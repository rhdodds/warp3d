c ***************************************************************
c *                                                             *
c * didriv.f                                                    *
c *                                                             *
c * domain_driver - drive the computation of j-integral values  *
c *                 using the domain integral technique.        *
c *                                                             *
c *                 or:                                         *
c *                                                             *
c *                 drive the computation of mixed-mode         *
c *                 stress intensity factors and t-stress       *
c *                 using the interaction integral.             *
c *                                                             *
c *                 last modified: 4/21/22 rhd                  *
c *                                                             *
c ***************************************************************
c
c
      subroutine didriv
c
      use global_data, only : out, nonode, noelem, nodof,
     a  dstmap, props, iprops, lprops, bits, nelblk,
     b  ltmstp, lsldnm, stname,
     c  scoords => c, sdispl => u, svelocities => v,
     d  saccel => a! old common.main
c
      use constants
c
      use main_data, only : incmap, incid, crdmap, output_packets,
     a                      packet_file_no, inverse_incidences,
     b                      fgm_node_values,
     c                      elems_to_blocks, temper_nodes,
     d                      temper_elems, initial_state_option,
     e                      fgm_node_values_used, initial_state_step,
     f                      initial_stresses_input
c
      use j_data, only : debug_driver, domain_type, crack_plane_normal,
     a  front_order, q_values, expanded_front_nodes,
     b  num_front_nodes, qvals_given, last_compr, front_element_list,
     c  comput_i, comput_j, front_q_area, front_length, compr_q_list,
     d  e_front, front_nodes, domain_min_i, front_list_length,
     e  domain_min_j, nu_front, fgm_e, fgm_nu, process_initial_state,
     f  block_seg_curves, seg_snode_e, seg_snode_nu, domain_origin,
     g  snode_alpha_ij, print_elem_values, domain_max_j, print_totals,
     h  domain_max_j, domain_avg_j, domain_max_i, domain_avg_i,
     i  static_min, static_max, static_avg, first_domain, 
     j  nowring, j_storage, i_storage, j_from_ks, ks_from_j, ring_list,
     k  max_domain_rings, q_element_maps, num_auto_rings, static_j,
     l  j_geonl_formulation, j_linear_formulation,
     m  temperatures_on_model, domain_node_sets, 
     n  compute_domain_extents

c
      implicit none
c
c                 locals. allocatables automatically released
c
      integer :: num_rings, block, iword, list_entry, snode, iout, i,
     &           q_map_len, ring, dof, orig_node, elemno, last_ring
      integer, allocatable :: q_node_map(:), q_new_map(:,:),
     &                        q_old_map(:,:)
      logical :: error, bad_domain, last_domain, user_def_ct, ldum,
     &           setup_node_props, bad
      double precision :: nx, ny, nz
      real :: rword
      real, parameter :: e_tol=0.001
      equivalence ( rword, iword ) ! still needed for
!                                    integers-reals mixed in same vector
c
      if( debug_driver ) write(out,9000) ! set by user input
c
c               0. process user-defined initial state ?
c
      process_initial_state = initial_state_option .and.
     &                        ltmstp >= initial_state_step
c
c               1. check validity of domain before starting computations.
c
      bad_domain = .false.
c
      if( domain_type .lt. 1 ) then
         bad_domain = .true.
         write(out,9110)
      end if
c
      if( front_order .lt. 1 .or. front_order .gt. 2 ) then
         bad_domain = .true.
         write(out,9115)
      end if
c
      nx = crack_plane_normal(1)
      ny = crack_plane_normal(2)
      nz = crack_plane_normal(3)
      if( abs(sqrt(nx*nx+ny*ny+nz*nz)-1.0d0) .gt. 0.001d0 ) then
         bad_domain = .true.
         write(out,9120)
      end if
c
c               2. perform exhaustive check on consistency of
c                  1) number of front nodes, 2) front interpolation
c                  order, and 3) domain type.
c
      call dickdm( bad_domain, out, scoords, crdmap, nonode, noelem )
      if( bad_domain ) then
        write(out,9140); call didrive_release( out ); return
      end if
c
c               3. allocate space for a q-value at each
c                  structure node. allocate space for expanded lists
c                  of coincident front nodes at each front location.
c
      allocate( q_values(nonode) )
      q_values(1:nonode) = rzero
      if( allocated( expanded_front_nodes ) ) 
     &   deallocate( expanded_front_nodes )
      allocate( expanded_front_nodes(num_front_nodes) )
      expanded_front_nodes(1:num_front_nodes)%node_count = 0
c
c               4. output header information for domain
c
      call diheadr
c
c               5. when q-values are given by the user,
c                  expand the compressed list of q-values
c                  into a list of length number of structure
c                  nodes.
c
      if( qvals_given ) then
        do list_entry = 1, last_compr
          snode           = compr_q_list(list_entry,1)
          iword           = compr_q_list(list_entry,2)
          q_values(snode) = rword
        end do
      end if
c
c               6. set up the domain. for automatic domains,
c                  set q-values at front nodes. then find coincident
c                  front nodes and set their q-values.
c
      call distup( scoords, crdmap, nonode, debug_driver, out )
c
c               7a. create vector with list of elements on the
c                   crack front for this domain
c
      call di_cf_elem( num_front_nodes, front_list_length )
c
c               7b. use crack front element to set if J computations
c                   will use small or large strain formulation.
c                   all elements in domain must match be
c                   the same formulation.
c
      elemno = front_element_list(1) ! element on front
      j_geonl_formulation  = lprops(18,elemno)
      j_linear_formulation = .not. j_geonl_formulation
c
c               8. at point on front where integral is being computed,
c                  build the global -> crack rotation matrix.
c                  gather coordinates and displacements of
c                  crack-front nodes, and rotate them to local
c                  crack-front system.
c
      call dimrot( scoords, crdmap, sdispl, dstmap, debug_driver, out )
c
c               8c. calculate strain e33 at node at domain origin.
c                   this is for T-stress calculations using the
c                   interaction integral. linear-elastic analysis only
c
      if( comput_i ) call di_calc_e33( out )
c
c               8c. calculate properties of a curve passing through
c                   the front nodes. these will be used to compute
c                   distance 'r' from integration points to a curved
c                   crack front. linear-elastic analysis for
c                   interaction integrals
c
      if( comput_i ) call di_calc_curvature( debug_driver, out )
c
c               9. compute area under the q-function over
c                  that part of crack front for this domain. the
c                  area must be > 0 else fatal error in domain
c                  (user forgot to set q-values on front nodes)
c
      call di_front_q_area( scoords, crdmap, out, debug_driver,
     &                      error, bad_domain )
      if( bad_domain ) then
         write(iout,9210)
         call didrive_release( out )
         return
      end if
      write(out,9200) front_q_area, front_length
c
c               10. formerly set logical flags to indicate if the nodal
c                   velocities and accelerations are all zero for
c                   this load step.
c                   * no longer used *
c
c
c               11. Need to do this when temperatures have been
c                   specified in model. Possibly leading to computation
c                   of J_6

c                   Not needed if the user has input initial stresses
c                   or set an initial state for  J-processing to accommodate
c                   prior thermo-mechanical processing. In such
c                   cases those effects, FGMs, thermal, .. will be
c                   included thru J_7, J_8
c
c                   Build the node average value of thermal expansion
c                   coefficient. for temperature-dependent material
c                   properties, also build the node average value of
c                   young's modulus and poisson's ratio. for
c                   temperature-independent material properties, values
c                   of e and nu are obtained within dicmj.f
c                   nodal properties are needed for domain
c                   integral computations to compute spatial derivatives
c                   within the domain.
c
c
      temperatures_on_model = any( temper_nodes .ne. zero ) .or.
     &                        any( temper_elems .ne. 0 )
      setup_node_props =  temperatures_on_model
      if( process_initial_state .or. initial_stresses_input ) 
     &      setup_node_props = .false.
c
      if( setup_node_props ) call di_node_props_setup( 1 ) !  MPI parallel
c
c              12. Needed when FGM nodal properties are actually used
c                  for elements in model, when the user
c                  has defined an initial state, or temperature
c                  dependent stress-strain curves are used (effectively
c                  makes the material an FGM), or input initial
c                  stresses.
c
c                  Needed for I integrals with temperatures in model.
c
c                  Even if not needed, the code builds small, dummy
c                  arrays to satisfy checks.
c
c                  Build the nodal averages of strain energy
c                  density (W) and nodal values of displacement
c                  gradient (3x3) at t_n relative to coordinates at t=0.
c                  These terms are used to calculate the
c                  derivative of W wrt crack local X and derivatives
c                  of the gradients wrt crack local X -- both at
c                  integration points for J(7) and J(8).
c
      call di_setup_J7_J8( 1 ) ! MPI parallel
c
c              12b. at point on front where integral is being computed,
c                   collect young's modulus and poisson's ratio. this
c                   assumes that all elements connected to this crack-front
c                   node have identical, homogeneous material properties,
c                   or that fgm material properties have been assigned
c                   to the model. for homogeneous material, "props" contains
c                   material data. for fgms, read data from "fgm_node_values."
c                   for temperature-dependent properties, segmental data
c                   arrays contain the properties.
c
      orig_node = front_nodes(domain_origin)
      elemno    = inverse_incidences(orig_node)%element_list(1)
      block     = elems_to_blocks(elemno,1)
      e_front   = props(7,elemno)
      nu_front  = props(8,elemno)
      if( fgm_e )  e_front  = fgm_node_values(orig_node,1)
      if( fgm_nu ) nu_front = fgm_node_values(orig_node,2)
      if( temperatures_on_model .and.
     &    allocated(block_seg_curves) ) then
           if( block_seg_curves(block) ) then
              e_front     = seg_snode_e(orig_node)
              nu_front    = seg_snode_nu(orig_node)
           end if
      end if
c
c               12c. at point on front where integral is being computed,
c                    props(7,elemno) and props(8,elemno) may be zero
c                    if the element has been killed. in this case,
c                    use the first nonzero values of e and nu.
c
      if( e_front .lt. e_tol ) then
         bad = .true.
         do i = 1, noelem
            e_front  = props(7,i)
            nu_front = props(8,i)
            if( e_front .gt. e_tol ) then
              bad = .false.
              write(out,9220) e_front, nu_front
              exit
            end if
         end do
         if( bad ) then
           write(out,9230)
           return
         end if
      end if
c
c               12d. output info for domain about temperature at crack
c                    front and alpha values. Also E, nu at front
c
      if( temperatures_on_model ) then
       write(out,9905) orig_node,
     &                 temper_nodes(orig_node) + temper_elems(elemno)
       if( allocated( snode_alpha_ij ) ) then
        write(out,9910) orig_node, (snode_alpha_ij(orig_node,i),i=1,6)
       end if
      end if
c
      write(out,9900) orig_node, e_front, orig_node, nu_front
c
      if( print_elem_values .or. print_totals ) then
         if( comput_i ) write(out,9950)
      end if
c
c              13a. init various tracking values for for J and I
c
      domain_min_j       =  1.0d30
      domain_max_j       = -1.0d30
      domain_avg_j       =  zero
      domain_min_i(1:10) =  1.0d30
      domain_max_i(1:10) = -1.0d30
      domain_avg_i(1:10) =  zero
      static_min         =  1.0d30
      static_max         = -1.0d30
      static_avg         =  zero
c
c              13b. separate drivers for comput values
c                   on 1 domain for user-defined q-values, and
c                   for defining-computing over automatic domains

      if( qvals_given ) then
         call didrive_user_q_values ! done with J, I calcs. time to output
      else
         call didrive_auto_domains  !   14. ....
      end if
c
c             15a. write j-integral and i-integral data
c                  to standard output
c
      call di_write_std_out( ltmstp, stname, lsldnm, out )
c
c             15b. write j-integral and i-integral data
c                   to packets
c
      if( output_packets ) call di_write_packets(ltmstp,stname,lsldnm)
c
c
      write(out,9160)
c
c              16. release arrays used for both user defined and
c                  automatic domains. call routines to get
c                  releases on onther than rank 0 for MPI
c
      call didrive_release( out )  ! only does rank 0 for MPI
c
      if( temperatures_on_model ) call di_node_props_setup( 2 )
      call di_setup_J7_J8( 2)
c
 9000 format(//,'>> Entered domain_driver...')
 9110 format(//,'>>> Domain type (a->d) not specified')
 9115 format(//,'>>> Front interpolation order not specified')
 9120 format(//,'>>> Crack plane normal not unit length')
 9140 format(//,'>>> Invalid domain definition. Skipped..')
 9060 format(/,1x,i5,2x,10(2x,e11.4),3x,'(',i3,')')
 9160 format(/,'>>> Completed domain integral calculations',/)
 9200 format(/,' area under q along front:       ',e13.6,
     & /,      ' length of crack front segment:  ',e13.6)
 9210 format(//,'>>> Area under q-function along the crack front',
     &          ' is not >0.',/,
     &          '    J-value cannot be computed. Check that q-values',
     &          ' are properly specified',
     &        /,'    for the front nodes. Domain integration procedure',
     &          ' terminated....',//)
 9220 format(/,'>>> Warning: a zero E exists at node on front where',
     &      /, '             J or I is being computed. The alternative',
     &      /, '             values used for E and nu are: ',2d14.6,//)
 9230 format(/,'>>> Warning: a zero E exists at node on front where',
     &      /, '             J or I is being computed. A non-zero E',
     &      /, '             cannot be found in element properties.',
     &      /, '             J-I computations skipped.',//)
 9900 format(' Young''s modulus at crack front node  ',i8,' : ',e13.6,
     &    /, ' Poisson''s ratio at crack front node  ',i8,' : ',e13.6 )
 9905 format(' temperature at crack front node      ',i8,' : ',e13.6 )
 9910 format(' alpha values at crack front node     ',i8,' : ',
     &       6(e13.6,2x) )
 9950 format(/,1x,'output sequence:',
     &       /,5x,'interaction integral value',5x,'auxiliary field',
     &       /,5x,'I_KI  ',25x,'plane stress',
     &       /,5x,'I_KI  ',25x,'plane strain',
     &       /,5x,'I_KII ',25x,'plane stress',
     &       /,5x,'I_KII ',25x,'plane strain',
     &       /,5x,'I_KIII',25x,'anti-plane shear',
     &       /,5x,'I_T11 ',25x,'plane stress',
     &       /,5x,'I_T11 ',25x,'plane strain',
     &       /,5x,'I_T13 ',25x,'anti-plane shear',/)
c
      return
c
      contains
c     ========
c
c              contains routines for didriv to simplify logic
c
      subroutine didrive_user_q_values
      implicit none
c
      first_domain = .true.
      nowring = 0
      allocate( j_storage(1,11) )
      allocate( i_storage(1,13,8) )
      allocate( j_from_ks(1,2) )
      allocate( ks_from_j(1,3) )
      j_storage(1,1:11)     = zero
      i_storage(1,1:13,1:8) = zero
      j_from_ks(1,1:2)      = zero
      ks_from_j(1,1:3)      = zero
      call dicmj   !   do the computation over elements - uses MPI
      last_domain = .true.
c
      return
c
      end subroutine didrive_user_q_values
c
      subroutine didrive_auto_domains
      implicit none
c
      double precision :: domain_extents(max_domain_rings,6)
c
c
c              14. user wants automatic construction of domains.
c
c              14a. get last ring at which output will be
c                   printed. domains are always generated starting
c                   at ring 1 but j-values may not be computed for
c                   every domain.
c
      last_ring = max_domain_rings
      do i = max_domain_rings, 1, -1
         if( ring_list(i) .eq. 1 ) then
            last_ring = i
            exit
         end if
      end do
c
c              14b. allocate arrays needed to support construction/
c                   definition of the domains. for type 4, we need
c                   only a nodal bit map. for types 1-3,
c                   we need two sets of nodal bit maps. each set
c                   has 3 maps of length to record all structure nodes.
c                   element list stored in the common vector.
c
      if( .not. allocated( q_element_maps ) ) then
          allocate( q_element_maps(noelem/30+1) )
          q_element_maps(1:noelem/30+1) = 0
      end if
      q_map_len = nonode/30 + 1 ! 30 bits/word used
      if( domain_type .ne. 4 ) then
         allocate( q_new_map(q_map_len,3) )
         allocate( q_old_map(q_map_len,3) )
         q_new_map(1:q_map_len,1:3) = 0
         q_old_map(1:q_map_len,1:3) = 0
      else
         allocate( q_node_map(q_map_len) )
      end if
c
c              14c. set up to accumulate statistics for
c                   computed domain values.
c
      allocate( j_storage(num_auto_rings,11) )
      allocate( i_storage(num_auto_rings,13,8) )
      allocate( j_from_ks(num_auto_rings,2) )
      allocate( ks_from_j(num_auto_rings,3) )
      j_storage(1:num_auto_rings,1:11)     = zero
      i_storage(1:num_auto_rings,1:13,1:8) = zero
      j_from_ks(1:num_auto_rings,1:2)      = zero
      ks_from_j(1:num_auto_rings,1:3)      = zero
      first_domain = .true.
      static_j     = .false.
c
c              14d. loop over all domains. construct definition of
c                   the domain (q-values, element list). call driver
c                   to actually calculate value for domain.
c
      do ring = 1, last_ring
         if( domain_type .eq. 4 ) then
            call diexp4( ring, q_node_map, nonode, noelem, incmap,
     &                   iprops, incid, out, q_map_len, bits,
     &                   compute_domain_extents, domain_extents  )
         else
            call diexp13( ring, nonode, noelem, incmap, iprops, incid,
     &                    out, bits, q_new_map, q_old_map, q_map_len,
     &                    compute_domain_extents, domain_extents )
        end if
         if( ring_list(ring) .eq. 1 ) then
            last_domain = ring .eq. last_ring
            nowring = ring 
            call dicmj  !  do calcs over elements, runs MPI parallel
            first_domain = .false.
         end if
      end do
c
      if( compute_domain_extents ) then 
        write(out,9000)
        write(out,9005)
        do ring = 1, last_ring
          write(out,9010) ring, domain_extents(ring,1:6)
        end do
      end if
c
      return
c
 9000 format(/," Physical Extents of Domains in Global Coordinates")
 9005 format("    Domain        X-min          Y-min         Z-min",
     &  "          X-max          Y-max          Z-max")
 9010 format(5x, i3, 6f15.5)     
c
      end subroutine didrive_auto_domains
c
      end subroutine didriv
c
c
c
c
c              independent support routine used by didriv.
c
c ***************************************************************
c *                                                             *
c *                      didrive_release                        *
c *                                                             *
c *                 last modified: 1/6/22 rhd                   *
c *                                                             *
c ***************************************************************
c
      subroutine didrive_release( iout )
      use j_data
      implicit none
c
      integer :: iout, i
      integer, parameter :: num_flags=19
      integer :: flags(num_flags)
c
      flags = 0
c
c              16a. integers, logicals
c
      if( allocated( compr_q_list ) )
     &  deallocate( compr_q_list, stat=flags(1) )
c
      if( allocated( domain_node_sets ) )
     &  deallocate( domain_node_sets, stat=flags(2) )
c
      if( allocated( expanded_front_nodes ) )
     &  deallocate( expanded_front_nodes, stat=flags(3) )
c
      if( allocated( q_element_maps ) )
     &  deallocate( q_element_maps, stat=flags(4) )
c
      if( allocated( count_alpha ) )
     &  deallocate( count_alpha, stat=flags(5) )
c
      if( allocated( extrap_counts ) )
     &  deallocate( extrap_counts, stat=flags(6) )
c
      if( allocated( front_element_list ) )
     &  deallocate( front_element_list, stat=flags(7) )
c
      if( allocated( block_seg_curves ) )
     &  deallocate( block_seg_curves, stat=flags(8))
c
c              16b. reals
c
      if( allocated( q_values ) )
     &  deallocate( q_values, stat=flags(9) )
c
      if( allocated( seg_snode_e ) )
     &  deallocate( seg_snode_e, stat=flags(10) )
c
      if( allocated( seg_snode_nu ) )
     &  deallocate( seg_snode_nu, stat=flags(11) )
c
      if( allocated( snode_alpha_ij ) )
     &  deallocate( snode_alpha_ij, stat=flags(12) )
c
c              16c. doubles
c
      if( allocated( swd_at_nodes ) )
     &  deallocate( swd_at_nodes, stat=flags(13) )
c
      if( allocated( strain_at_nodes ) )
     &  deallocate( strain_at_nodes, stat=flags(14) )
c
      if( allocated( j_storage ) )
     &  deallocate( j_storage, stat=flags(15) )
c
      if( allocated( j_from_ks ) )
     &  deallocate( j_from_ks, stat=flags(16) )
c
      if( allocated( ks_from_j ) )
     &  deallocate( ks_from_j, stat=flags(17))
c
      if( allocated( displ_grad_at_nodes ))
     &  deallocate( displ_grad_at_nodes, stat=flags(18) )
c
      if( allocated( i_storage ) )
     &  deallocate( i_storage, stat=flags(19) )
c
      if( any( flags .ne. 0 ) ) then
        write(iout,9000)
        do i = 1, num_flags
          if( flags(i) .ne. 0 ) write(iout,9010) i, flags(i)
        end do
        write(iout,9020)
        call die_abort
      end if
c
      return
c
 9000 format(/,"FATAL ERROR in didrive_release")
 9010 format(10x, "flag no, deallocate error code: ",2i6)
 9020 format(/,".... Job terminated" )
c
      end
c ***************************************************************
c
c ***************************************************************
c *                                                             *
c * domain_header - output info at start of domain computations *
c *                                                             *
c *                 last modified: 3/8/2018 rhd                 *
c *                                                             *
c ***************************************************************
c
      subroutine diheadr
      use global_data, only : stname, lsldnm, ltmstp, total_model_time,
     &                        out
      use j_data, only : domain_id
      implicit none
c
c                 locals
c
      real, external :: wcputime
c
      write(out,1116)
      write(out,1017) stname
      write(out,1118) lsldnm, ltmstp, total_model_time
      write(out,1119) domain_id(1:24), wcputime(1)
c
      return
c
 1116 format(//,25x,'domain-integral post-processing',/,25x,
     &   31('-') )
 1017 format(/,1x,'structure: ',a8)
 1118 format(1x,'loading: ',a8 ,5x, 'step: ', i5, 8x,
     & 'model analysis time: ',e14.6 )
 1119 format(1x,'domain id: ',a24,6x,
     &   'current elapsed wall time: ',f8.1)
c
      end
c ***************************************************************
c *                                                             *
c * domain_check  - exhaustive testing of domain defintiion for *
c *                 consistency                                 *
c *                                                             *
c *                last modified: 1/10/22 rhd                   *
c *                                                             *
c ***************************************************************
c
      subroutine dickdm( bad_domain, iout, scoord, coord_map,
     &                   nonode, noelem )
      use constants
      use j_data, only : front_nodes, num_front_nodes, rings_given,
     &                   qvals_given, front_order, domain_node_sets,
     &                   compr_q_list, q_element_maps, domain_type
      implicit none
c
c              parameters
c
      integer :: nonode, noelem, iout, coord_map(*)
      double precision :: scoord(*)
      logical :: bad_domain
c
c              locals
c
      integer :: numfrn, frtint, node_id, i, j, node_err, dtype, ftype,
     &           node_set_id 
      double precision :: x, y, z, x1, y1, z1, distance, cum_distance
      integer inttbl(3,7)
      data  ((inttbl(i,j),i=1,3),j=1,7) 
     & / 4,0,0, 1,0,0,  3,1,0,  4,0,1,  4,3,0,  4,0,0,  4,4,3 /
      character(len=1), parameter :: labs(4) = [ 'a','b','c','d' ]
      logical :: consis, linear, quad, cubic, header, ok, user_def_ct
c
c             check consistency of the domain definition
c
      header = .true.
      user_def_ct = front_nodes(1) .lt. 0
      if( num_front_nodes .eq. 0 ) then
         write(iout,9900)
         header = .false.
         write(iout,9910)
         bad_domain = .true.
      end if
c
      if( rings_given .and. qvals_given ) then
         if( header ) then
             write(iout,9900)
             header = .false.
         end if
         write(iout,9920)
         bad_domain = .true.
         return
      end if
c
c             check consistency of number of front nodes specified
c             with interpolation order.
c
      numfrn = num_front_nodes
      frtint = front_order
      consis = .false.
      if( frtint .eq. 1 ) consis = .true.
      if( frtint .eq. 2 )
     &   consis = ( numfrn .ge. 3 ) .and. ( mod(numfrn,2) .eq. 1 )
      if( frtint .eq. 3 )
     &   consis = ( numfrn .ge. 4 ) .and.
     &            ( numfrn .eq. ( (numfrn/3) +1 + 2*(numfrn/3) ) )
      if( .not. consis ) then
         if( header ) then
             write(iout,9900)
             header = .false.
         end if
         write(iout,9980)
         bad_domain = .true.
         return
      end if
c
c             verify that the order of front nodes given
c             defines a sequence of increasing distances
c             from the first node.
c
      if( numfrn .gt. 1 ) then
        if( user_def_ct ) then
          node_set_id = -front_nodes(1)
          node_id = domain_node_sets(node_set_id)%node_list(1)
          x1 = scoord(coord_map(node_id))
          y1 = scoord(coord_map(node_id)+1)
          z1 = scoord(coord_map(node_id)+2)
        else
          x1 = scoord(coord_map(front_nodes(1)))
          y1 = scoord(coord_map(front_nodes(1))+1)
          z1 = scoord(coord_map(front_nodes(1))+2)
        end if
        cum_distance = zero
        do i = 2, numfrn
          if( user_def_ct ) then
            node_set_id = -front_nodes(i)
            node_id = domain_node_sets(node_set_id)%node_list(1)
            x = scoord(coord_map(node_id))
            y = scoord(coord_map(node_id)+1)
            z = scoord(coord_map(node_id)+2)
          else
            x = scoord(coord_map(front_nodes(i)))
            y = scoord(coord_map(front_nodes(i))+1)
            z = scoord(coord_map(front_nodes(i))+2)
          end if
          distance = sqrt( (x-x1)**2 + (y-y1)**2 + (z-z1)**2 )
          if( distance .le. cum_distance ) then
            if( header ) then
              write(iout,9900)
              header = .false.
            end if
            if( user_def_ct ) then
               node_err = node_id
            else
               node_err = front_nodes(i)
            end if
            write(iout,9994) node_err
            bad_domain = .true.
            return
          end if
          cum_distance = distance
        end do
      end if
c
c             make checks for user specification of q-values
c
      if( qvals_given  ) then
          ok = .false.
          do i = 1, nonode
            if( compr_q_list(i,2) .gt. rzero ) ok = .true.
          end do
          if( .not. ok ) then
            if( header ) then
               write(iout,9900)
               header = .false.
            end if
            write(iout,9930)
            bad_domain = .true.
          end if
          ok = .false.
          do i = 1, noelem/30+1
            if( q_element_maps(i) .gt. 0 ) ok = .true.
          end do
          if( .not. ok ) then
            if( header ) then
               write(iout,9900)
               header = .false.
            end if
            write(iout,9940)
            bad_domain = .true.
          end if
          return
      end if
c
c             make checks for automatic generation of
c             q-values.
c
c             for automatic generation, check consistency of
c             front node interpolation with the function
c             type specified. if user does not specify value,
c             assume function type based on front order
c             interpolation and number of front nodes. if we get
c             this far in the checks, the number of front nodes
c             is consistent with the order of interpolation.
c
      numfrn = num_front_nodes
      frtint = front_order
      ftype  = domain_type
c
      if( ftype .eq. 0 ) then
         if( numfrn .gt. 7 ) then
            dtype = 4
          else
            dtype = inttbl(frtint,numfrn)
         end if
         domain_type = dtype
         write(iout,9990) labs(dtype)
         return
      end if
c
c             user specified a function type. check it for
c             consistency with the front nodes and interpolation.
c
      ftype  = domain_type
      linear = frtint .eq. 1
      quad   = frtint .eq. 2
      cubic  = frtint .eq. 3
      if( ftype .eq. 1 .or. ftype .eq. 3 ) then
         consis = (numfrn .eq. 2 .and. linear) .or.
     &            (numfrn .eq. 3 .and. quad )  .or.
     &            (numfrn .eq. 4 .and. cubic )
      end if
      if( ftype .eq. 2 ) then
         consis = (numfrn .eq. 3 .and. linear) .or.
     &            (numfrn .eq. 5 .and. quad )  .or.
     &            (numfrn .eq. 7 .and. cubic )
      end if
      if( ftype .eq. 4 ) consis = .true.
      if( .not. consis ) then
         write(iout,9992)
         bad_domain = .true.
         return
      end if
c
c
      return
c
 9900 format(/,' >>>>> domain definition error:')
 9910 format(/, ' >>> no front nodes specified',// )
 9920 format(/, ' >>> cannot specify q-values & automatic',//)
 9930 format(/, ' >>> no q-values specified',//)
 9940 format(/, ' >>> no element list specified',//)
 9970 format(/, ' >>> cannot use normal orientation for function',
     &                   ' type:d',//)
 9980 format(/, ' >>>the number of front nodes and the order of'
     &   /,     '    interpolation are inconsistent',//)
 9990 format(/, ' >>> default value adopted for function type: ',
     &              a1,//)
 9992 format(/, ' >>> the specified function type and the number of',
     &   /,     '       front nodes/interpolation order are '
     &   /,     '       inconsistent',//)
 9994 format(/, ' >>> the front nodes are not in order of monotonically',
     &   /,     '       increasing distance from first front node...'
     &   /,     '       first front node in error: ',i7,
     &   /,     "       Did you forget the keyword 'sets' in the ",
     &          "'front nodes ...'",
     &   /,     "       command for domains that use multiple nodes ",
     &          "at each front location?", /)
c
      end
c
c **********************************************************************
c *                                                                    *
c * di_cf_elem - create a list of elements incident on the crack front *
c *                                                                    *
c *              written by:    mcw                                    *
c *              last modified: 11/23/22 rhd                           *
c *                                                                    *
c **********************************************************************
c
      subroutine di_cf_elem( num_front_nodes, front_list_length )
c
      use global_data, only : out, noelem
      use main_data, only: inverse_incidences
      use j_data, only : front_element_list, expanded_front_nodes 
c
      implicit none
c
c             parameters
c
      integer :: num_front_nodes, front_list_length       
      logical, parameter :: local_debug = .false.
c
c             local arguments
c
      integer :: i, j, elem, fnode_index
      integer, allocatable, dimension(:) :: elem_flags(:)
c
      if( local_debug ) write(out,9005)
      if( allocated( front_element_list ) ) 
     &    deallocate( front_element_list )
      allocate( elem_flags(noelem) )
      elem_flags = 0
      if( local_debug ) then
         write(out,9010) 1
         write(out,*) "     ... num_front_nodes: ", num_front_nodes
      end if
c
c             build list of all elements attached a front node in the
c             list. run over front node is parallel. for domains
c             with thousands of front nodes this saves time.
c
c$OMP PARALLEL DO PRIVATE( fnode_index ) SCHEDULE( dynamic )
      do fnode_index = 1, num_front_nodes
         call di_add_elem_to_list( fnode_index, elem_flags )
      end do
c$OMP END PARALLEL DO
      if( local_debug ) write(out,9010) 2
c
c             save final list of front elements in list of
c             exact size required.
c
      front_list_length = count( elem_flags > 0, dim=1 )
      allocate( front_element_list(1:front_list_length) )
      if( local_debug ) write(out,9015) 3, front_list_length
c
      j = 0 
      do elem = 1, noelem
       if( elem_flags(elem) > 0 ) then
          j = j + 1
          front_element_list(j) = elem
          cycle
        end if
      end do 
      deallocate( elem_flags ) 
      if( local_debug ) write(out,9010) 4
c
      if( j .ne. front_list_length ) then  ! confirming
        write(out,9000) front_list_length, j
        call die_gracefully
      end if
c
      return
 9000 format(/,">>>>>> FATAL INTERNAL ERROR in di_cf_elem",
     &   10x, 'front_list_length, j: ',2i10)
 9005 format(/,".... entered di_cf_elem .....")
 9010 format(5x,"@ ",i3)
 9015 format(5x,"@, list length ",i3,i8)
c
      end
c
c **********************************************************************
c *                                                                    *
c *         di_add_elem_to_list - support for di_cf_elem               *
c *                                                                    *
c *              written by:    rhd                                    *
c *              last modified: 11/23/22 rhd                           *
c *                                                                    *
c **********************************************************************
c
c
      subroutine di_add_elem_to_list( fnode_index, elem_flags )
c
      use main_data, only: inverse_incidences
      use j_data, only :  expanded_front_nodes 
c
      implicit none
c
      integer :: elem_flags(*), fnode_index
      integer :: j, k, fnode, numelems_connected, elem, 
     &           numexpanded_nodes
c
      numexpanded_nodes = expanded_front_nodes(fnode_index)%node_count
      do j = 1, numexpanded_nodes
        fnode = expanded_front_nodes(fnode_index)%node_list(j)
        numelems_connected = inverse_incidences(fnode)%element_count
        do k = 1, numelems_connected
           elem = inverse_incidences(fnode)%element_list(k)
c$OMP ATOMIC UPDATE
           elem_flags(elem) = elem_flags(elem) + 1
        end do
      end do
      return 
c
      end

c
c **********************************************************************
c *                                                                    *
c * di_cf_elem - create a list of elements incident on the crack front *
c *                                                                    *
c *              written by:    mcw                                    *
c *              last modified: 1/7/22 rhd                             *
c *                                                                    *
c **********************************************************************
c
      subroutine old_di_cf_elem( num_front_nodes, front_list_length )
c
      use main_data, only: inverse_incidences
      use j_data, only : front_element_list, expanded_front_nodes 
c
      implicit none
c
c             parameters
c
      integer :: num_front_nodes, front_list_length       
c
c             local arguments
c
      integer :: i, j, k, numexpanded_nodes, fnode, numelems_connected,
     &           elem, now_size, size
      integer, parameter :: start_list_size = 50
      integer, allocatable, dimension(:) :: new_list
c
      if( allocated( front_element_list ) ) 
     &    deallocate( front_element_list )
      allocate( front_element_list(1:start_list_size) )
      size = start_list_size
      now_size = 0
c
c             find and add non-duplicate front elements to list
c             that resizes as needed
c
      do i = 1, num_front_nodes
         numexpanded_nodes = expanded_front_nodes(i)%node_count
         do j = 1, numexpanded_nodes
            fnode = expanded_front_nodes(i)%node_list(j)
            numelems_connected = inverse_incidences(fnode)%element_count
            do k = 1, numelems_connected
               elem = inverse_incidences(fnode)%element_list(k)
               call old_di_cf_elem_add
            end do
         end do
      end do
c
c             save final list of front elements in list of
c             exact size required.
c
      front_list_length = now_size
      if( now_size .ne. size ) then
        size = now_size
        allocate( new_list(size) )
        new_list(1:now_size) = front_element_list(1:now_size)
        call move_alloc( new_list, front_element_list )
      end if
c
      return
c
      contains
c     ========
c
      subroutine old_di_cf_elem_add
      implicit none
c
      integer :: i
c
      if( now_size .eq. 0 ) then
         now_size = 1
         front_element_list(1) = elem
         return
      end if

      do i = 1, now_size
        if( front_element_list(i) == elem ) return
      end do
c
c              add to list. resize if needed. move_alloc frees old
c              front_element_list
c
      if( now_size == size ) then
        size = size * 2
        allocate( new_list(size) )
        new_list(1:now_size) = front_element_list(1:now_size)
        call move_alloc( new_list, front_element_list )
      end if
      now_size = now_size + 1
      front_element_list(now_size) = elem
c
      return
      end subroutine old_di_cf_elem_add
      end subroutine old_di_cf_elem
c
c **********************************************************************
c *                                                                    *
c *                    di_chk_proportionality (for J_cutoff)           *
c *                                                                    *
c *      solutions for load step now_step completed. cumulative        *
c *      factors for loading patterns and constraints were printed     *
c *      by calling routine.                                           *
c *                                                                    *
c *      store loading pattern names and factors for step 1 in a       *
c *      persistent table. for steps 2, 3, ... confirm                 *
c *      continued proportional loading -- can be complex !            *
c *                                                                    *
c *              written by:    rhd                                    *
c *              last modified: 6/26/22 rhd                            *
c *                                                                    *
c **********************************************************************
c
      subroutine di_chk_proportionality( now_step, pattern_count,
     &            all_pattern_ids, all_pattern_factors,
     &            sum_constraints )
c
      use global_data, only: out, nodof
      use main_data, only: cnstrn_in
      use j_data, only : J_cutoff_now_frnt_position, J_cutoff_exceeded,
     &  J_cutoff_step_1_constraint_factor, 
     &  J_cutoff_step_1_num_patterns, patterns_step_1,
     &  J_load_ratio_this_step 
      use constants
c
      implicit none
c
c             parameters
c
      integer, intent(in) :: now_step, pattern_count
      character(len=*), intent(in) :: all_pattern_ids(*)
      double precision, intent(in) :: all_pattern_factors(*),
     &                                sum_constraints
c
c             locals
c
      integer :: i, j
      double precision :: cons_sum_value, step_1_factor, now_factor,
     &                    constraints_ratio, tol, comparison_ratio
      logical :: cons_all_zero, first,
     &           mpcs_found, bad_ratio, now_factor_not_zero,
     &           step_1_factor_is_zero, constraints_only,
     &           step_1_factor_not_zero, now_factor_is_zero 
      logical, parameter :: local_debug = .false.
c
      if( local_debug ) then
       write(out,*)
     &   ".... entered di_chk_proportionality ...."
       write(out,*) '      pattern_count: ', pattern_count
      end if
c
c             initialize J cutoff should domains be computed now that 
c             step solution is converged. (not adaptive sub-steps).
c
      J_cutoff_now_frnt_position = 0
      J_load_ratio_this_step = zero 
c
c             for step 1, store the patterns and factors in the
c             table for J cutoff use. all defined
c             loading patterns are in the list for step 1, even
c             if never mentioned in a step definition.
c
      if( now_step == 1 ) then
        J_cutoff_step_1_constraint_factor = sum_constraints
        J_cutoff_step_1_num_patterns = pattern_count
        do i = 1, pattern_count ! could = 0, no loop
           patterns_step_1(i)%id = all_pattern_ids(i) 
           patterns_step_1(i)%factor = all_pattern_factors(i) 
        end do
        if( local_debug ) then
         write(out,9000) J_cutoff_step_1_constraint_factor,
     &       J_cutoff_step_1_num_patterns
         do i = 1, pattern_count ! could = 0, no loop
            write(out,9005) patterns_step_1(i)%id,
     &          patterns_step_1(i)%factor
         end do
        end if
        return
      end if
c
c             for step 2, 3, .. is loading for just completed step
c             proportional to loading in step 1 ? can get complex !
c   
      if( pattern_count .ne. J_cutoff_step_1_num_patterns ) then
        write(out,9010) 
        write(out,9015) pattern_count, J_cutoff_step_1_num_patterns
        write(out,9020)
        call die_gracefully
      end if
c
c             current step has 1 or more patterns. all defined
c             loading patterns on the model before step1 
c             are in the list for step 1, even
c             if never mentioned in a step definition.
c
      J_load_ratio_this_step = zero
      tol = thousandth * thousandth 
      first = .true.
      bad_ratio = .false.
      comparison_ratio = zero
c 
      do i = 1, pattern_count ! could = 0
        step_1_factor = patterns_step_1(i)%factor
        step_1_factor_is_zero = abs(step_1_factor) < tol
        step_1_factor_not_zero = abs(step_1_factor) > tol
        now_factor = all_pattern_factors(i)
        now_factor_not_zero =  abs(now_factor) > tol
        now_factor_is_zero =  abs(now_factor) < tol
        if( local_debug ) write(out,9130) i, step_1_factor,
     &                    now_factor
        if( step_1_factor_is_zero .and. now_factor_not_zero ) then
            bad_ratio = .true.
            write(out,9010)
            write(out,9100) patterns_step_1(i)%id
            cycle
        end if
        if( step_1_factor_not_zero .and. now_factor_is_zero ) then
            bad_ratio = .true.
            write(out,9010)
            write(out,9140) patterns_step_1(i)%id
            cycle
        end if
        if( step_1_factor_is_zero .or. now_factor_is_zero ) cycle
        if( first ) then
             comparison_ratio = now_factor / step_1_factor
             if( local_debug ) write(out,*) '.. i, comp ratio: ',i,
     &              comparison_ratio
             first = .false.
             cycle
        end if
        if( abs( now_factor/step_1_factor - comparison_ratio )
     &          > tol ) then
             bad_ratio = .true.
             write(out,9010)
             write(out,9110) patterns_step_1(i)%id, 
     &                       now_factor
             cycle
        end if
      end do ! on i
c
      if( bad_ratio ) then
            write(out,9020)
            call die_gracefully
      end if
      J_load_ratio_this_step = comparison_ratio
c
c
c             pattern loads thru this step are proportional to
c             those in step 1. if any constraints are non-zero,
c             the constraint factor for this step relative to
c             step 1 must have the same ratio as
c             to step 1.
c
      if( local_debug ) then
        write(out,9045) J_load_ratio_this_step
      end if
c 
      constraints_only = abs( J_load_ratio_this_step ) < tol
      if( constraints_only ) then ! only constraints loading
         constraints_ratio = sum_constraints / 
     &                      J_cutoff_step_1_constraint_factor
         J_load_ratio_this_step = constraints_ratio
         if( local_debug ) then
           write(out,9050) J_load_ratio_this_step
         end if
         write(out,9200)  J_load_ratio_this_step
         return
      end if
c             
c             patterns & constraints must be proportional to
c             step 1. if all constraints = 0, they don't matter
c             & we already have proportional patterns
c
      tol = thousandth * thousandth 
      call cons_sum( cons_sum_value, mpcs_found )
      cons_all_zero = cons_sum_value < tol
      if( local_debug ) write(out,9055) cons_sum_value, cons_all_zero 
      if( mpcs_found ) write(out,9120)
      if( cons_all_zero ) then
        write(out,9200)  J_load_ratio_this_step
        return
      end if
c
      constraints_ratio = sum_constraints / 
     &                      J_cutoff_step_1_constraint_factor
      if( abs(constraints_ratio - J_load_ratio_this_step) > tol ) then
        write(out,9010) 
        write(out,9060) constraints_ratio, J_load_ratio_this_step
        write(out,9020)
        call die_gracefully
      end if
      write(out,9200)  J_load_ratio_this_step
      return
c
c
 9000 format(" ... routine di_chk_proportionality ...",
     &   /,  "       step = 1, constraint factor: ",f8.5,
     &   /,  "       num patterns: ",i3 )
 9005 format(15x,a8,2x,f10.5)
 9010 format(/1x,'>>>>> Fatal Error:  non-proportional loading ',
     &          ' encountered ....',
     &       /1x,'                    J cutoff feature requires ',
     &      'that each step loading be proportional to step 1')
 9015 format(/,5x,'number of load patterns: ',i5,' for this step',
     & /,5x,'does not match the number for step 1: ',i5, //)
 9020 format(5x,'Job terminated ...')
 9035 format(5x,'non-equal ratios these patterns to step 1 values',
     &  /,   5x,'step pattern: ',a8,' bad ratio: ',f10.5)
 9045 format(5x,'after pattern checking. load ratio: ',f10.5)
 9050 format(5x,'constraints only. loading ratio: ',f10.5)
 9055 format(5x,'cons_sum, flag: ',d14.6,l3)
 9060 format(/,5x,'this step constraints ratio to step 1:  ',f10.5,
     &     /,5x   'does not match pattern ratio to step 1: ',f10.5,//)  
 9100 format(/,5x,'pattern: ',a8,' has 0.0 factor in step 1 and',
     &     ' cannot be changed to maintain proportional loading',// ) 
 9110 format(5x,'pattern: ',a8,' factor: ',f10.5,2x,
     &     'creates non-proportional ',
     &     ' loading relative to step 1',// ) 
 9120 format(/1x,'>>>>> Warning:  multi-point constraints are present',
     &       /1x,'                WARP3D cannot check that loads ',
     &  'remain proportional to step 1',
     &       /1x,'                for J cutoff capability',//)
 9130 format('... i, s1 factor, now factor: ',i3,2f10.5)
 9140 format(/,5x,'pattern: ',a8,' has non-zero factor in step 1 and',
     &     ' and is now zero.',// ) 
 9200 format('>>>>> INFO: step load ratio to step 1 is: ',f10.5,
     &     ' for use in J cutoff processing'//)
c
      end  
        
         

      

