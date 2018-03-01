c ***************************************************************
c *                                                             *
c * didriv.f                                                    *
c *    calling tree:                                            *
c *                                                             *
c *    didriv.f                                                 *
c *       -dickdm                                               *
c *       -diheadr                                              *
c *       -distup.f                                             *
c *       -di_cf_elem                                           *
c *       -dimrot.f                                             *
c *            -difrtn.f                                        *
c *                 -di1dsf.f                                   *
c *                 -difrts                                     *
c *       -di_calc_e33                                          *
c *       -difrar.f                                             *
c *            -di1dsf.f                                        *
c *       -di_expan_coeff_setup                                 *
c *            -di_node_expan_setup                             *
c *       -di_fgm_setup                                         *
c *            -di_nod_vals                                     *
c *                 -di_extrap_to_nodes                         *
c *                      -ndpts1                                *
c *                      -oulgf                                 *
c *       -dicmj.f                                              *
c *            -dielem_a.f                                      *
c *            -dielem_b.f                                      *
c *            -dielem_c.f                                      *
c *       -diexp4.f                                             *
c *       -diexp13.f                                            *
c *                                                             *
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
c *                 last modified: 9/10/15  by rhd              *
c *                                                             *
c ***************************************************************
c
c
      subroutine didriv
      use global_data ! old common.main
c
      use main_data, only : incmap, incid, crdmap, output_packets,
     &                      packet_file_no, inverse_incidences,
     &                      fgm_node_values, temperatures_ref,
     &                      elems_to_blocks, temper_nodes, temper_elems
      use j_data
c
      implicit integer (a-z)
c
c                 local declarations
c
      logical error, bad_domain, last_domain, user_def_ct, ldum
c
      double precision
     & nx, ny, nz, zero
      real rword, rzero, e_tol
c
      allocatable q_node_map(:), q_new_map(:,:), q_old_map(:,:)
c
      integer num_rings, block
      equivalence ( rword, iword )
c
      data rzero, zero, e_tol / 0.0, 0.0d0, 0.001 /
c
      debug_driver = .false.   ! in module J-data
      if ( debug_driver ) write(out,9000)
c
c               1. check validity of domain before starting computations.
c
      bad_domain = .false.
c
      if ( domain_type .lt. 1 ) then
         bad_domain = .true.
         write(out,9110)
      end if
c
      if ( front_order .lt. 1 .or. front_order .gt. 2 ) then
         bad_domain = .true.
         write(out,9115)
      end if
c
      nx = crack_plane_normal(1)
      ny = crack_plane_normal(2)
      nz = crack_plane_normal(3)
      if ( abs(sqrt(nx*nx+ny*ny+nz*nz)-1.0) .gt. 0.001 ) then
         bad_domain = .true.
         write(out,9120)
      end if
c
c               2. perform exhaustive check on consistency of
c                  1) number of front nodes, 2) front interpolation
c                  order, and 3) domain type.
c
      call dickdm( bad_domain, out, c, crdmap, nonode, noelem )
      if ( bad_domain ) then
        write(out,9140)
        return
      end if
c
c               3. allocate space for a q-value at each
c                  structure node. allocate space for expanded lists
c                  of coincident front nodes at each front location.
c
      allocate( q_values(nonode) )
      q_values(1:nonode) = rzero
      max_exp_front      = 200
      allocate( expanded_front_nodes(0:max_exp_front,num_front_nodes) )
      expanded_front_nodes(0:max_exp_front,1:num_front_nodes) = 0
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
      if ( qvals_given ) then
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
      call distup( c, crdmap, nonode, debug_driver, out )
c
c               7. allocate a vector of logicals and assign .true. for
c                  each element connected to a crack front node.
c
      allocate( crack_front_elem(1:noelem) )
      crack_front_elem(1:noelem) = .false.
      call di_cf_elem( crack_front_elem, num_front_nodes,
     &                 max_exp_front, expanded_front_nodes )
c
c               8. at point on front where integral is being computed,
c                  build the global->crack rotation matrix.
c                  gather coordinates and displacements of
c                  crack-front nodes, and rotate them to local
c                  crack-front system.
c
      call dimrot( c, crdmap, u, dstmap, debug_driver, out )
c
c               8c. calculate strain e33 at node at domain origin.
c                   this is for T-stress calculations using the
c                   interaction integral
c
      if( comput_i ) call di_calc_e33( out )
c
c               8c. calculate properties of a curve passing through
c                   the front nodes. these will be used to compute
c                   distance 'r' from integration points to a curved
c                   crack front.
c
      if( comput_i ) call di_calc_curvature( debug_driver, out )
c
c               9. compute area under the q-function over
c                  that part of crack front for this domain. the
c                  area must be >0 else fatal error in domain
c                  (user forgot to set q-values on front nodes)
c
      call di_front_q_area( c, crdmap, out, debug_driver,
     &                      error, bad_domain )
      if ( bad_domain ) then
         write(iout,9210)
         if ( allocated( q_values ) ) deallocate( q_values )
         if ( allocated( expanded_front_nodes ) )
     &                   deallocate( expanded_front_nodes )
         if ( allocated( crack_front_elem ) )
     &                   deallocate( crack_front_elem )
         return
      end if
      write(out,9200) front_q_area, front_length
c
c               10. set logical flags to indicate if the nodal
c                   velocities and accelerations are all zero for
c                   this load step. if so, some later computations
c                   can be skipped.
c
      process_velocities = .false.
      process_accels     = .false.
      do dof = 1, nodof
         if ( v(dof) .ne. zero ) process_velocities = .true.
         if ( a(dof) .ne. zero ) process_accels     = .true.
      end do
c
c            di_expan_coeff_setup
c
c               11. Build the node average value of thermal expansion
c                   coefficient. for temperature-dependent material
c                   properties, also build the node average value of
c                   young's modulus and poisson's ratio. for
c                   temperature-independent material properties, values
c                   of e and nu are obtained within dicmj.f
c                   nodal properties are needed for domain
c                   integral computations to compute spatial derivatives
c                   within the domain.
c
      if( temperatures .or. temperatures_ref )
     &       call di_node_props_setup (1, nonode, nelblk )
c
c              12. Build the nodal averages of strain energy
c                  density (stress work density) and strains.
c                  These terms are used to calculate the
c                  derivative of the strain energy density, which
c                  appears in the domain integral when material
c                  properties vary spatially (e.g. fgms). The nodal
c                  values calculated in di_fgm_setup will be used
c                  to compute their spatial derivatives at integration
c                  points.
c
      call di_fgm_setup ( 1, nonode, out, temperatures )
c
c               12b. at point on front where integral is being computed,
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
      if( fgm_e ) then
         e_front  = fgm_node_values(orig_node,1)
      end if
      if( fgm_nu ) then
         nu_front = fgm_node_values(orig_node,2)
      end if
      if( temperatures .or. temperatures_ref ) then
         if( block_seg_curves(block) ) then
            e_front     = seg_snode_e(orig_node)
            nu_front    = seg_snode_nu(orig_node)
         end if
      end if
c
c          props(7,elemno) and props(8,elemno) may be zero
c          if the element has been killed. in this case,
c          use the first nonzero values of e and nu.
c
      if( e_front .lt. e_tol ) then
         do i = 1, noelem
            e_front  = props(7,i)
            nu_front = props(8,i)
            if( e_front .gt. e_tol ) exit
         end do
      end if
c
      if( temperatures .or. temperatures_ref ) then
         write(out,9905) orig_node,
     &                   temper_nodes(orig_node) + temper_elems(elemno)
         write(out,9910) orig_node, (snode_alpha_ij(orig_node,i),i=1,6)
      end if
c
      write(out,9900) orig_node, e_front, orig_node, nu_front
c
      if( print_elem_values .or. print_totals ) then
         if( comput_i ) write(out,9950)
      end if
c
c              13. if q-values given by user, we compute domain
c                  integral right now. otherwise, we set up a
c                  loop to generate q-values for automatic
c                  domains and their computation.
c
      domain_min_j       =  1.0e30
      domain_max_j       = -1.0e30
      domain_avg_j       =  zero
      domain_min_i(1:10) =  1.0e30
      domain_max_i(1:10) = -1.0e30
      domain_avg_i(1:10) =  zero
      static_min         =  1.0e30
      static_max         = -1.0e30
      static_avg         =  zero
c
      if ( qvals_given ) then
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
          call dicmj
          last_domain = .true.
          go to 8000
      end if
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
         if ( ring_list(i) .eq. 1 ) then
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
      if ( .not. allocated( q_element_maps ) ) then
          allocate( q_element_maps(noelem/30+1) )
          q_element_maps(1:noelem/30+1) = 0
      end if
      q_map_len = nonode/30 + 1 ! 30 bits/word used
      if ( domain_type .ne. 4 ) then
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
         if ( domain_type .eq. 4 ) then
            call diexp4( ring, q_node_map, nonode, noelem, incmap,
     &                   iprops, incid, out, q_map_len, bits )
         else
            call diexp13( ring, nonode, noelem, incmap, iprops, incid,
     &                    out, bits, q_new_map, q_old_map, q_map_len )
         end if
         if ( ring_list(ring) .eq. 1 ) then
            last_domain = ring .eq. last_ring
            nowring = ring
            call dicmj
            first_domain = .false.
         end if
      end do
c
c              14e. release allocatable arrays for automatic domains
c
      if ( domain_type .ne. 4 ) then
         deallocate( q_new_map )
         deallocate( q_old_map )
      else
         deallocate( q_node_map )
      end if
c
 8000 continue
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
c                  automatic domains
c
      deallocate( q_values )
      deallocate( expanded_front_nodes )
      deallocate( crack_front_elem )
      if ( allocated( j_storage ) ) deallocate( j_storage )
      if ( allocated( i_storage ) ) deallocate( i_storage )
      if ( allocated( j_from_ks ) ) deallocate( j_from_ks )
      if ( allocated( ks_from_j ) ) deallocate( ks_from_j )
      if ( allocated( compr_q_list ) )    deallocate( compr_q_list )
      if ( allocated( q_element_maps ) )  deallocate( q_element_maps )
      if ( allocated( node_set ) )        deallocate( node_set )
      if ( allocated( block_seg_curves ) ) deallocate( block_seg_curves)
      if( temperatures .or. temperatures_ref )
     &    call di_node_props_setup ( 2, nonode, nelblk )
      call di_fgm_setup( 2, nonode, out, ldum )
c
 9000 format(//,'>> Entered domain_driver...')
 9110 format(//,'>>> Domain type (a->d) not specified')
 9115 format(//,'>>> Front interpolation order not specified')
 9120 format(//,'>>> Crack plane normal not unit length')
 9140 format(//,'>>> Invalid domain definition. Skipped..')
 9060 format(/,1x,i5,2x,10(2x,e11.4),3x,'(',i3,')')
 9160 format(/,'>>> Completed domain integral calculations')
 9200 format(/,' area under q along front:       ',e13.6,
     & /,      ' length of crack front segment:  ',e13.6)
 9210 format(//,'>>> Area under q-function along the crack front',
     &          ' is not >0.',/,
     &          '    J-value cannot be computed. Check that q-values',
     &          ' are properly specified',
     &        /,'    for the front nodes. Domain integration procedure',
     &          ' terminated....',//)
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
      end
c
c
c           support routines used by didriv.
c
c
c ***************************************************************
c *                                                             *
c * domain_header - output info at start of domain computations *
c *                                                             *
c ***************************************************************
c
c
      subroutine diheadr
      use global_data ! old common.main
      use j_data
      implicit integer (a-z)
c
c                 local declarations
c
      real wcputime
      external wcputime
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
 1119 format(1x,'domain id: ',a24,6x,'current elapsed wall time: ',f8.1)
c
      end
c ***************************************************************
c *                                                             *
c * domain_check  - exhaustive testing of domain defintiion for *
c *                 consistency                                 *
c *                                                             *
c ***************************************************************
c
      subroutine dickdm( bad_domain, iout, scoord, coord_map,
     &                   nonode, noelem )
      use j_data
      implicit integer (a-z)
c
      double precision
     & scoord(*)
      dimension  coord_map(*)
c
      double precision
     &  x, y, z, x1, y1, z1, zero, distance, cum_distance
      real rzero
c
      dimension  inttbl(3,7)
      character(len=1) :: labs(4)
      logical    consis, linear, quad, cubic, header, bad_domain,
     &           ok, user_def_ct
      data       labs / 'a','b','c','d' /
      data       inttbl / 4,0,0,  1,0,0,  3,1,0, 4,0,1,
     &                    4,3,0,  4,0,0,  4,4,3 /
      data zero, rzero / 0.0, 0.0 /
c
c             check consistency of the domain definition
c
      header = .true.
      user_def_ct = front_nodes(1) .lt. 0
      if ( num_front_nodes .eq. 0 ) then
         write(iout,9900)
         header = .false.
         write(iout,9910)
         bad_domain = .true.
      end if
c
      if ( rings_given .and. qvals_given ) then
         if ( header ) then
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
      if ( frtint .eq. 1 ) consis = .true.
      if ( frtint .eq. 2 )
     &   consis = ( numfrn .ge. 3 ) .and. ( mod(numfrn,2) .eq. 1 )
      if ( frtint .eq. 3 )
     &   consis = ( numfrn .ge. 4 ) .and.
     &             ( numfrn .eq. ( (numfrn/3) +1 + 2*(numfrn/3) ) )
      if ( .not. consis ) then
         if ( header ) then
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
      if ( numfrn .gt. 1 ) then
        if ( user_def_ct ) then
          node_id = node_set(-front_nodes(1),1)
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
          if ( user_def_ct ) then
            node_id = node_set(-front_nodes(i),1)
            x = scoord(coord_map(node_id))
            y = scoord(coord_map(node_id)+1)
            z = scoord(coord_map(node_id)+2)
          else
            x = scoord(coord_map(front_nodes(i)))
            y = scoord(coord_map(front_nodes(i))+1)
            z = scoord(coord_map(front_nodes(i))+2)
          end if
          distance = sqrt( (x-x1)**2 + (y-y1)**2 + (z-z1)**2 )
          if ( distance .le. cum_distance ) then
            if ( header ) then
              write(iout,9900)
              header = .false.
            end if
            if ( user_def_ct ) then
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
      if ( qvals_given  ) then
          ok = .false.
          do i = 1, nonode
            if ( compr_q_list(i,2) .gt. rzero ) ok = .true.
          end do
          if ( .not. ok ) then
            if ( header ) then
               write(iout,9900)
               header = .false.
            end if
            write(iout,9930)
            bad_domain = .true.
          end if
          ok = .false.
          do i = 1, noelem/30+1
            if ( q_element_maps(i) .gt. 0 ) ok = .true.
          end do
          if ( .not. ok ) then
            if ( header ) then
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
      if ( ftype .eq. 0 ) then
         if ( numfrn .gt. 7 ) then
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
      if ( ftype .eq. 1 .or. ftype .eq. 3 ) then
         consis = (numfrn .eq. 2 .and. linear) .or.
     &            (numfrn .eq. 3 .and. quad )  .or.
     &            (numfrn .eq. 4 .and. cubic )
      end if
      if ( ftype .eq. 2 ) then
         consis = (numfrn .eq. 3 .and. linear) .or.
     &            (numfrn .eq. 5 .and. quad )  .or.
     &            (numfrn .eq. 7 .and. cubic )
      end if
      if ( ftype .eq. 4 ) consis = .true.
      if ( .not. consis ) then
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
c * di_cf_elem - create a logical vector whose entries are .true. for  *
c *              elements incident on the crack tip, and .false. for   *
c *              those that are not. dicmj will use this info to set a *
c *              flag for each element that is analyzed by dielem. if  *
c *              a user includes the domain integral command           *
c *              'omit crack front elements for fgms yes', the flag    *
c *              will cause terms7 and 8 to be set to zero.            *
c *                                                                    *
c *              written by:    mcw                                    *
c *              last modified: 09/06/01 by mcw                        *
c *                                                                    *
c **********************************************************************
c
      subroutine di_cf_elem( crack_front_elem, num_front_nodes,
     &                       max_exp_front, expanded_front_nodes )
      use main_data, only: inverse_incidences
      implicit none
c
c             dummy arguments
c
      integer num_front_nodes, max_exp_front,
     &        expanded_front_nodes(0:max_exp_front,1:num_front_nodes)
      logical crack_front_elem(*)
c
c             local arguments
c
      integer i, j, k, numexpanded_nodes, fnode, numelems_connected,
     &        elem
c
      do i = 1, num_front_nodes
         numexpanded_nodes = expanded_front_nodes(0,i)
         do j = 1, numexpanded_nodes
            fnode = expanded_front_nodes(j,i)
            numelems_connected = inverse_incidences(fnode)%element_count
            do k = 1, numelems_connected
               elem = inverse_incidences(fnode)%element_list(k)
               crack_front_elem(elem) = .true.
            end do
         end do
      end do
c
      return
      end
c
c **********************************************************************
c *                                                                    *
c * di_node_props_setup - obtain alpha values at nodes. for            *
c *                       temperature-dependent properties, also       *
c *                       compute e and nu values at nodes. this       *
c *                       routine replaces di_expan_coeff_setup,       *
c *                       and the routine it calls, di_node_props,     *
c *                       replaces di_node_expan_coeff.                *
c *                                                                    *
c *                       written by:    mcw                           *
c *                       last modified: 06/11 by mcm                  *
c *                                 fixes various MPI issues           *
c *                                                                    *
c **********************************************************************
c
      subroutine di_node_props_setup( do_itin, nonodeaa, nelblkin )
      use global_data ! old common.main
      use j_data, only: count_alpha, snode_alpha_ij, seg_snode_e,
     &                  seg_snode_nu, block_seg_curves,
     &                  process_temperatures, front_nodes, domain_origin
      implicit integer (a-z)
c
c             dummy arguments
c
      integer do_itin, nonodeaa, nelblkin
c
c             local arguments
c
      integer orig_node
      real zero
      data zero / 0.0 /
c
      integer nonodea,do_it
c
      orig_node = front_nodes(domain_origin)
c
c          if we are running in MPI:
c             alert slaves to enter this routine and allocate data
c             structures for handling temperature loads in the J
c             calculations.
c          if we are running in serial:
c             this is just a dummy routine which immediately returns
c
      call wmpi_alert_slaves( 30 )
c
c          broadcast the do_it variable
c
c           Bugfix, need to make a local copy of nelblk, nonodea, and do_it
      if (myid .eq. 0) then
            do_it = do_itin
            nelblk = nelblkin
            nonodea = nonodeaa
      end if
c
      call wmpi_bcast_int ( do_it )
      call wmpi_bcast_int ( nonodea )
c     It's this one below
      call wmpi_bcast_int ( nelblk )
      call wmpi_bcast_int ( orig_node )
c
c          do_it = 1 -- allocate and fill J temp data structures
c          do_it = 2 -- deallocate J temp data structures
c
      if ( do_it .eq. 1 ) then
         allocate( count_alpha(nonodea) )
         allocate( snode_alpha_ij(nonodea,6) )
         allocate( seg_snode_e(nonodea) )
         allocate( seg_snode_nu(nonodea) )
         allocate( block_seg_curves(nelblk) )
         count_alpha(1:nonodea)        = 0
         snode_alpha_ij(1:nonodea,1:6) = zero
         seg_snode_e(1:nonodea)        = zero
         seg_snode_nu(1:nonodea)       = zero
         block_seg_curves(1:nelblk)   = .false.
         call di_node_props ( count_alpha, snode_alpha_ij, seg_snode_e,
     &                        seg_snode_nu, block_seg_curves,
     &                        process_temperatures, orig_node )
c
      else
         if( allocated( count_alpha) )     deallocate( count_alpha )
         if( allocated( snode_alpha_ij ) ) deallocate( snode_alpha_ij )
         if( allocated( seg_snode_e ) )    deallocate( seg_snode_e )
         if( allocated( seg_snode_nu ) )   deallocate( seg_snode_nu )
         if( allocated( block_seg_curves)) deallocate( block_seg_curves)
      end if
c
      return
      end
c
c **********************************************************************
c *                                                                    *
c * di_node_props - build average nodal values of material properties  *
c *                 at each node in the model. all values are          *
c *                 single-precision reals because double-precision    *
c *                 accuracy is unnecessary.                           *
c *                                                                    *
c *                 written by: mcw                                    *
c *                 last modified: 06/11 by mcm                        *
c *                                 (fixes various MPI issues)         *
c *                                                                    *
c **********************************************************************
c
      subroutine di_node_props ( count_alpha, snode_alpha_ij,
     &                           seg_snode_e, seg_snode_nu,
     &                           block_seg_curves, process_temperatures,
     &                           orig_node )
      use global_data ! old common.main
      use main_data, only : incmap, incid, fgm_node_values,
     &                      elems_to_blocks, temper_nodes,
     &                      temper_nodes_ref, temper_elems,
     &                      inverse_incidences
c      use j_data, only : e_front, nu_front, alpha_front
      use segmental_curves, only: seg_curve_table, seg_curves_type
      implicit integer (a-z)
c
c             dummy arguments
c
      integer count_alpha(*), orig_node
      real snode_alpha_ij(nonode,*), seg_snode_e(*), seg_snode_nu(*)
      logical block_seg_curves(*), process_temperatures
c
c             local arguments
c
      integer block, first_elem, processor, block_props(nelblk),
     &        bit_flags, curve_set, first_curve, curve_set_type, snode,
     &        incptr, i, span, elemno, num_enodes, enode,
     &        num_curves_in_set, curve_no, block_num
      double precision
     &     elem_uniform_temp, enode_temper, alpha_temper, e_temper,
     &     nu_temper, linear_interpolate
      double precision,
     &  dimension (:), allocatable :: curve_temp, curve_es, curve_nus,
     &                 curve_alphas
      real alphax, alphay, alphaz, alphaxy, alphaxz, alphayz, neg_99,
     &     fgm_tol, sum, zero, one, rcount
      external linear_interpolate
      logical debug, seg_alphas, fgm_alphas, constant_alphas, flag
c          I seemingly have no choice but to do this to fix a bug
      real, dimension(:), allocatable :: alpha_copy
c
      data  neg_99, fgm_tol, zero, one / -99.0, 1.0, 0.0, 1.0 /
c
c
      debug                 = .false.
c
      constant_alphas       = .false.
      fgm_alphas            = .false.
      seg_alphas            = .false.
      process_temperatures  = .false.
      block_props(1:nelblk) = 0
c
c             loop through element blocks.
c             check type of property assignment in each block.
c             property assignment: element nodes (fgm)   = 1
c                                  constant in element   = 2
c                                  temperature dependent = 3
c
      do block = 1, nelblk
         span       = elblks(0,block)
         first_elem = elblks(1,block)
         processor  = elblks(2,block)
         block_num  = elems_to_blocks(first_elem,1)
c
c             skip block if not handled by current processor
c
         if( processor .ne. myid ) cycle
         if( debug ) write(out,*) "processor ", processor, " block ",
     &                             block, " first_elem ", first_elem
c
c             1. check for fgm alphas. the identifier 'fgm_mark' for
c             fgms assigned in inmat.f is -99.0. elprp.f then assigns
c             this value to props(9,*), props(13,*) and props(34,*).
c
         fgm_alphas = abs( props(9,first_elem) - neg_99 ) .le. fgm_tol
         if( fgm_alphas ) then
            process_temperatures = .true.
            block_props(block)   = 1
            call di_fgm_alphas( span, first_elem, count_alpha,
     &                          snode_alpha_ij )
            goto 100
         end if
c
c             2. check for alphas that are constant throughout element.
c                (this is check #2 because fgm assignment of alpha
c                places '-99' in the props array.)
c
         alphax  = props(9,first_elem)
         alphay  = props(13,first_elem)
         alphaz  = props(34,first_elem)
         alphaxy = props(35,first_elem)
         alphaxz = props(36,first_elem)
         alphayz = props(37,first_elem)
         sum     = abs(alphax)  + abs(alphay)  + abs(alphaz)
     &           + abs(alphaxy) + abs(alphaxz) + abs(alphayz)
         if ( sum .gt. zero ) then
            constant_alphas      = .true.
            process_temperatures = .true.
            block_props(block)   = 2
            call di_constant_alphas( span, first_elem, count_alpha,
     &                               snode_alpha_ij )
            goto 100
         end if
c
c             3. check for temperature-dependent material properties
c                assigned using segmental curves. curve_set_type
c                curve_set_type =0, temperature-independent properties
c                               =1, temperature-dependent properties
c                               =2, strain-rate dependent properties
c
         bit_flags  = iprops(24,first_elem)
         seg_alphas = .true.
         if ( iand(bit_flags,4) .eq. 0 ) seg_alphas = .false.
         if( seg_alphas ) then
            seg_alphas = .false.
            curve_set         = iprops(21,first_elem)
            num_curves_in_set = seg_curve_table(1,curve_set)
            first_curve       = seg_curve_table(2,curve_set)
            curve_set_type    = seg_curves_type(first_curve)
            if( curve_set_type .eq. 1 ) then
               seg_alphas              = .true.
               process_temperatures    = .true.
               block_seg_curves(block) = .true.
               block_props(block)      = 3
               call di_seg_alpha_e_nu( curve_set, num_curves_in_set,
     &                   first_curve, curve_set_type, span, first_elem,
     &                   count_alpha, snode_alpha_ij,
     &                   seg_snode_e, seg_snode_nu )
            end if
         end if
c
 100     continue
         if( debug ) then
            write(out,*)"block ",block," prop type ",block_props(block)
         end if
      end do
c
c             average sum of alpha values at each node in structure.
c
c             for parallel execution:
c             nodes lying outside domain managed by current
c             processor will have zero alpha values, but this won't
c             affect computations.
c
      do snode = 1,nonode
         if ( count_alpha(snode) .ne. 0 ) then
            rcount                  = one / real(count_alpha(snode))
            snode_alpha_ij(snode,1) = snode_alpha_ij(snode,1) * rcount
            snode_alpha_ij(snode,2) = snode_alpha_ij(snode,2) * rcount
            snode_alpha_ij(snode,3) = snode_alpha_ij(snode,3) * rcount
            snode_alpha_ij(snode,4) = snode_alpha_ij(snode,4) * rcount
            snode_alpha_ij(snode,5) = snode_alpha_ij(snode,5) * rcount
            snode_alpha_ij(snode,6) = snode_alpha_ij(snode,6) * rcount
         end if
      end do
c
c             send flags to root processor. if the flag is true on any
c             slave processor, set the flag on root to true.
c
      call wmpi_redlog( process_temperatures, numprocs )
      do block = 1, nelblk
         call wmpi_redlog( block_seg_curves(block), numprocs )
      end do
c
c             for temperature-dependent material properties, we need to
c             make sure that root has the array data from which to obtain
c             e_front, nu_front, and alpha_front.
c
c           Fix a wacky bug where snode_alpha_ij gets messed up
      allocate(alpha_copy(6))
      alpha_copy(1:6) = snode_alpha_ij(orig_node,1:6)
c      call wmpi_send_real( seg_snode_e(orig_node), 1 )
      call wmpi_send_real_new( seg_snode_e( orig_node) , 1)
c      call wmpi_send_real( seg_snode_nu(orig_node), 1 )
      call wmpi_send_real_new( seg_snode_nu(orig_node), 1 )
c      call wmpi_send_real( snode_alpha_ij(orig_node,1), 1 )
c      call wmpi_send_real( snode_alpha_ij(orig_node,2), 1 )
c      call wmpi_send_real( snode_alpha_ij(orig_node,3), 1 )
c      call wmpi_send_real( snode_alpha_ij(orig_node,4), 1 )
c      call wmpi_send_real( snode_alpha_ij(orig_node,5), 1 )
c      call wmpi_send_real( snode_alpha_ij(orig_node,6), 1 )
       call wmpi_send_real_new(alpha_copy,6)
      if (myid .eq. 0) then
            snode_alpha_ij(orig_node,1:6) = alpha_copy(1:6)
      end if
c
      deallocate(alpha_copy)
c
 200  if( debug ) then
         write(out, 500)
         if( fgm_alphas ) write(out,501)
         if( constant_alphas ) write(out,502)
         if( seg_alphas ) write(out,503)
         do snode = 1,nonode
            write(out,504) snode, (snode_alpha_ij(snode,i), i=1,6)
         end do
      end if
c
      return
c
 500  format(//,'***alpha values at structure nodes***')
 501  format('   --fgm alphas detected')
 502  format('   --constant alphas detected')
 503  format('   --temperature-dependent alphas detected')
 504  format('node',2x,i7,2x,6(e14.4,2x))
c
      end
c
c **********************************************************************
c *                                                                    *
c * di_fgm_alphas - retrieve nodal values of alpha assigned through    *
c *                 'fgm' input. this routine was separated from       *
c *                 di_node_props for clarity.                         *
c *                                                                    *
c *                 written by: mcw                                    *
c *                 last modified: 02/05 by mcw                        *
c *                                                                    *
c **********************************************************************
c
      subroutine di_fgm_alphas( span, first_elem, count_alpha,
     &                          snode_alpha_ij )
      use global_data ! old common.main
      use main_data, only : incmap, incid, fgm_node_values
      implicit integer (a-z)
c
c             dummy variables
c
      integer span, first_elem, count_alpha(*)
      real snode_alpha_ij(nonode,6)
c
c             local variables
c
      integer i, elemno, num_enodes, incptr, enode, snode
      real zero
      data zero / 0.0 /
c
c             assign fgm alpha values to nodes of elements in current
c             block from fgm_alpha_values
c
c             loop through elements in block
c
      do i = 1, span
         elemno     = first_elem + ( i - 1 )
         num_enodes = iprops(2,elemno)
         incptr     = incmap(elemno)
c
c             loop through nodes on current element and assign fgm alpha.
c             we sum element contributions to make logic easier for
c             averaging values when some blocks use fgm alphas and some
c             use constant alphas. fgm alphas must be isotropic.
c
         do enode = 1, num_enodes
            snode = incid( incptr + enode - 1 )
            count_alpha(snode)      = count_alpha(snode) + 1
            snode_alpha_ij(snode,1) = snode_alpha_ij(snode,1)
     &                              + fgm_node_values(snode,3)
            snode_alpha_ij(snode,2) = snode_alpha_ij(snode,2)
     &                              + fgm_node_values(snode,3)
            snode_alpha_ij(snode,3) = snode_alpha_ij(snode,3)
     &                              + fgm_node_values(snode,3)
            snode_alpha_ij(snode,4) = zero
            snode_alpha_ij(snode,5) = zero
            snode_alpha_ij(snode,6) = zero
         end do
      end do
      return
      end
c
c **********************************************************************
c *                                                                    *
c * di_constant_alphas - retrieve element values of alpha assigned     *
c *                      through element definitions, and add values   *
c *                      to nodal sums. the average nodal value is     *
c *                      then computed in the calling routine.         *
c *                      this routine was separated from               *
c *                      di_node_props for clarity.                    *
c *                                                                    *
c *                 written by: mcw                                    *
c *                 last modified: 02/05 by mcw                        *
c *                                                                    *
c **********************************************************************
c
      subroutine di_constant_alphas( span, first_elem, count_alpha,
     &                               snode_alpha_ij )
      use global_data ! old common.main
      use main_data, only : incmap, incid
      implicit integer (a-z)
c
c             dummy variables
c
      integer span, first_elem, count_alpha(*)
      real snode_alpha_ij(nonode,6)
c
c             local variables
c
      integer i, elemno, num_enodes, incptr, enode, snode
      real alphax, alphay, alphaz, alphaxy, alphaxz, alphayz
c
c             assign alpha values from element constant values defined
c             in the material property definition for the block
c
c             loop through elements in block
c
      do i = 1, span
         elemno     = first_elem + ( i - 1 )
         num_enodes = iprops(2,elemno)
         incptr     = incmap(elemno)
         alphax     = props(9,elemno)
         alphay     = props(13,elemno)
         alphaz     = props(34,elemno)
         alphaxy    = props(35,elemno)
         alphaxz    = props(36,elemno)
         alphayz    = props(37,elemno)
c
c             loop through nodes on current element and assign them the
c             alpha value of the element.
c
         do enode = 1, num_enodes
            snode                   = incid( incptr + enode - 1 )
            count_alpha(snode)      = count_alpha(snode) + 1
            snode_alpha_ij(snode,1) = snode_alpha_ij(snode,1)
     &                              + alphax
            snode_alpha_ij(snode,2) = snode_alpha_ij(snode,2)
     &                              + alphay
            snode_alpha_ij(snode,3) = snode_alpha_ij(snode,3)
     &                              + alphaz
            snode_alpha_ij(snode,4) = snode_alpha_ij(snode,4)
     &                              + alphaxy
            snode_alpha_ij(snode,5) = snode_alpha_ij(snode,5)
     &                              + alphaxz
            snode_alpha_ij(snode,6) = snode_alpha_ij(snode,6)
     &                              + alphayz
         end do
      end do
      return
      end
c
c **********************************************************************
c *                                                                    *
c * di_seg_alpha_e_nu - retrieve nodal values of alpha, e and nu       *
c *                     assigned through temperature-dependent         *
c *                     segmental curves. to simplify logic, nodal     *
c *                     alpha values are added to sum and averaged     *
c *                     in calling routine. e and nu values at nodes   *
c *                     are not averaged. this routine                 *
c *                     was separated from di_node_props for clarity.  *
c *                                                                    *
c *                 written by: mcw                                    *
c *                 last modified: 02/05 by mcw                        *
c *                                                                    *
c **********************************************************************
c
      subroutine di_seg_alpha_e_nu( curve_set, num_curves_in_set,
     &                  first_curve, curve_set_type, span, first_elem,
     &                  count_alpha, snode_alpha_ij, seg_snode_e,
     &                  seg_snode_nu )
      use global_data ! old common.main
      use main_data, only : incmap, incid, temper_elems, temper_nodes,
     &                      temper_nodes_ref
      use segmental_curves, only : seg_curve_table, seg_curves_value,
     &                             seg_curves_ym, seg_curves_nu,
     &                             seg_curves_alpha
      implicit integer (a-z)
c
c             dummy variables
c
      integer curve_set, curve_set_type, num_curves_in_set, first_curve,
     &        span, first_elem, count_alpha(*)
      real snode_alpha_ij(nonode,6), seg_snode_e(*), seg_snode_nu(*)
c
c             local variables
c
      integer i, curve_no, elemno, num_enodes, incptr, enode, snode
      double precision
     &     elem_uniform_temp, enode_temper, alpha_temper, e_temper,
     &     nu_temper, linear_interpolate
      double precision,
     &  dimension (:), allocatable :: curve_temp, curve_es, curve_nus,
     &                 curve_alphas
      real zero
      logical debug
      external linear_interpolate
      data zero / 0.0 /
c
      debug = .false.
c
      if( debug ) then
         write(out,1000) curve_set, curve_set_type, num_curves_in_set,
     &                   first_curve
      end if
c
c             allocate temporary arrays.
c             get alpha, e and nu values from segmental curves table.
c
      if( allocated(curve_temp) ) deallocate( curve_temp )
      if( allocated(curve_alphas) ) deallocate( curve_alphas )
      if( allocated(curve_es) ) deallocate( curve_es )
      if( allocated(curve_nus) ) deallocate( curve_nus )
      allocate( curve_temp(num_curves_in_set),
     &          curve_alphas(num_curves_in_set),
     &          curve_es(num_curves_in_set),
     &          curve_nus(num_curves_in_set) )
c
      do i = 1, num_curves_in_set
         curve_no        = seg_curve_table(i+1,curve_set)
         curve_temp(i)   = seg_curves_value(curve_no)
         curve_es(i)     = seg_curves_ym(curve_no)
         curve_nus(i)    = seg_curves_nu(curve_no)
         curve_alphas(i) = seg_curves_alpha(curve_no)
c
         if( debug ) then
            write(out,2000) i, curve_no, curve_temp(i), curve_es(i),
     &                   curve_nus(i), curve_alphas(i)
         end if
c
      end do
c
c             loop through elements in current block
c
      do i = 1, span
         elemno            = first_elem + ( i - 1 )
         num_enodes        = iprops(2,elemno)
         incptr            = incmap(elemno)
         elem_uniform_temp = temper_elems(elemno)
c
c             loop through nodes on each element. get temperature at
c             each node. do not remove reference temperatures of model
c             nodes (temper_nodes_ref) -- we need absolute temperatures
c             for temperature-dependent values.
c
         do enode = 1, num_enodes
            snode               = incid( incptr + enode - 1 )
            count_alpha(snode)  = count_alpha(snode) + 1
            enode_temper        = temper_nodes(snode)
     &                          + elem_uniform_temp
c
c             use nodal temperature to interpolate material properties
c             from segmental curves.
c
            alpha_temper = linear_interpolate( enode_temper,
     &                     num_curves_in_set, curve_temp, curve_alphas )
            e_temper     = linear_interpolate( enode_temper,
     &                     num_curves_in_set, curve_temp, curve_es )
            nu_temper    = linear_interpolate( enode_temper,
     &                     num_curves_in_set, curve_temp, curve_nus )
c
c             store temperature-dependent material values. we sum element
c             alpha contributions to make logic easier for averaging.
c             temperature-dependent e and nu values are not summed for
c             averaging.
c
            seg_snode_e(snode)      = e_temper
            seg_snode_nu(snode)     = nu_temper
            snode_alpha_ij(snode,1) = snode_alpha_ij(snode,1)
     &                              + alpha_temper
            snode_alpha_ij(snode,2) = snode_alpha_ij(snode,2)
     &                              + alpha_temper
            snode_alpha_ij(snode,3) = snode_alpha_ij(snode,3)
     &                              + alpha_temper
            snode_alpha_ij(snode,4) = zero
            snode_alpha_ij(snode,5) = zero
            snode_alpha_ij(snode,6) = zero
c
            if( debug ) then
               write(out,3000) elemno, enode, snode,temper_nodes(snode),
     &         elem_uniform_temp, temper_nodes_ref(snode), enode_temper,
     &         alpha_temper, e_temper, nu_temper
            end if
c
         end do
      end do
      if( allocated(curve_temp) ) deallocate( curve_temp )
      if( allocated(curve_alphas) ) deallocate( curve_alphas )
      if( allocated(curve_es) ) deallocate( curve_es )
      if( allocated(curve_nus) ) deallocate( curve_nus )
c
      return
 1000 format(/,'curve_set         = ', i2,/,'curve_set_type    = ', i2,
     &       /,'num_curves_in_set = ', i2,/,'first_curve       = ', i2 )
 2000 format(/,'i               = ', i2,   /,'curve_no        = ', i2,
     &       /,'curve_temp(i)   = ', e11.4,/,'curve_es(i)     = ',e11.4,
     &       /,'curve_nus(i)    = ', e11.4,/,'curve_alphas(i) = ',e11.4)
 3000 format(/,'elemno                  = ', i8,
     &       /,'enode                   = ', i8,
     &       /,'snode                   = ', i8,
     &       /,'temper_nodes(snode)     = ', e11.4,
     &       /,'elem_uniform_temp       = ', e11.4,
     &       /,'temper_nodes_ref(snode) = ', e11.4,
     &       /,'enode_temper(snode)     = ', e11.4,
     &       /,'alpha_temper(snode)     = ', e11.4,
     &       /,'e_temper(snode)         = ', e11.4,
     &       /,'nu_temper(snode)        = ', e11.4 )
      end
c
c **********************************************************************
c *                                                                    *
c * di_fgm_setup - allocate data structures for two terms used         *
c * in the calculation of the derivative of the stress work density.   *
c * these are: nodal values of stress work density and strain.         *
c *                                                                    *
c *                   written by: mcw                                  *
c *                last modified: 09/01 by mcw                         *
c *                                                                    *
c **********************************************************************
c
      subroutine di_fgm_setup( do_it, nonode, out, temperatures )
      use main_data, only: fgm_node_values, temperatures_ref
      use j_data, only: extrap_counts, swd_at_nodes, strain_at_nodes,
     &                  fgm_e, fgm_nu, comput_i, comput_j
      implicit none
c
c          dummy variables
c
      integer do_it, nonode, out
      logical temperatures
c
c          local variables
c
      double precision
     &     zero
      integer i, node, flags(3), ios
      logical debug, callit
c
      data zero /0.0d0/
c
      debug = .false.
c
c          if we are running in MPI:
c             alert slaves to enter this routine and allocate data
c             structures for handling temperature loads in the J
c             calculations.
c          if we are running in serial:
c             this is just a dummy routine which immediately returns
c
      call wmpi_alert_slaves( 15 )
c
c          broadcast variables
c
      call wmpi_bcast_int ( do_it )
      call wmpi_bcast_int ( nonode )
      call wmpi_bcast_int ( out )
      call wmpi_bcast_log ( temperatures )
      call wmpi_bcast_log ( temperatures_ref )
      call wmpi_bcast_log ( comput_j )
      call wmpi_bcast_log ( comput_i )
c
c          if do_it = 2 -- deallocate fgm node-value data structures
c
      if( do_it .eq. 2 ) then
         flags(1:3) = 0
         if( allocated(extrap_counts) )
     &            deallocate( extrap_counts, stat = flags(1) )
         if( allocated(swd_at_nodes) )
     &            deallocate( swd_at_nodes, stat = flags(2) )
         if( allocated(strain_at_nodes) )
     &            deallocate( strain_at_nodes, stat = flags(3) )
         if( any(flags(1:3) .ne. 0) ) then
            write(out, 100)
            call die_abort
         end if
         return
      end if
c
c          do_it = 1 -- allocate and fill fgm node-value data structures
c
      flags(1:3) = 0
      if( allocated( extrap_counts ) )
     &    deallocate( extrap_counts, stat = flags(1) )
      if( allocated( swd_at_nodes ) )
     &    deallocate( swd_at_nodes, stat = flags(2) )
      if( allocated( strain_at_nodes ) )
     &    deallocate( strain_at_nodes, stat = flags(3) )
      if( any(flags(1:3) .ne. 0) ) then
         write(out, 100)
         call die_abort
      end if
c
c           even if fgm properties have not been assigned,
c           extrap_counts, swd_at_nodes and strain_at_nodes
c           must be allocated for use in dielem_*.f.
c
      flags(1:3) = 0
      allocate( extrap_counts(1:nonode),       stat = flags(1) )
      allocate( swd_at_nodes(1:nonode),        stat = flags(2) )
      allocate( strain_at_nodes(1:6,1:nonode), stat = flags(3) )
      if( any(flags(1:3) .ne. 0) ) then
         write(out, 200)
         call die_abort
      end if
      extrap_counts(1:nonode)       = zero
      swd_at_nodes(1:nonode)        = zero
      strain_at_nodes(1:6,1:nonode) = zero
c
c           determine if 'fgm' values of e or nu have been assigned
c           at the nodes. if not, terms 7 and 8 of j will be set to zero,
c           so skip extrapolation process in di_nod_vals. dicmj.f also
c           uses these flags for calculating auxiliary fields.
c
      fgm_e  = .false.
      if( allocated(fgm_node_values) ) then
         do node=1,nonode
            if( fgm_node_values(node,1) .ne. zero ) then
               fgm_e  = .true.
               exit
            end if
         end do
      end if
c
      fgm_nu = .false.
      if( allocated(fgm_node_values) ) then
         do node=1,nonode
            if( fgm_node_values(node,2) .ne. zero ) then
               fgm_nu = .true.
               exit
            end if
         end do
      end if
c
c             it is necessary to call di_nod_vals:
c
c                1) for J-integral computations when e and nu have
c                   been assigned using fgm material properties.
c
c                2) for I-integral computations involving thermal
c                   loads.
c
      callit = .false.
      if( fgm_e .and. comput_j )            callit = .true.
      if( fgm_nu .and. comput_j )           callit = .true.
      if( comput_i .and. temperatures )     callit = .true.
      if( comput_i .and. temperatures_ref ) callit = .true.
c
c          extrapolate stress work density and (total) strains to nodes.
c
      if( callit ) call di_nod_vals( extrap_counts, swd_at_nodes,
     &                               strain_at_nodes )
c
      return
 100  format('>>>>>deallocation error in di_fgm_setup')
 200  format('>>>>>allocation error in di_fgm_setup')
      end
c
c **********************************************************************
c *                                                                    *
c * di_nod_vals - build average nodal values of the strain             *
c *          and stress work density for all nodes in the model.       *
c *                                                                    *
c *                   written by: mcw                                  *
c *                last modified: 09/01 by mcw                         *
c *                                                                    *
c **********************************************************************
c
      subroutine di_nod_vals( extrap_counts, swd_at_nodes,
     &                        strain_at_nodes )
      use global_data ! old common.main
c
      use main_data, only : incmap, incid, elems_to_blocks
      use elem_block_data, only : urcs_n_blocks, eps_n_blocks
c
      implicit integer (a-z)
c
c
c         dummy variables
c
      dimension extrap_counts(*)
c
      double precision
     &     swd_at_nodes(*), strain_at_nodes(6,*)
c
c         local variables
c
      double precision
     &     swd_at_gpts(mxgp), strain_at_gpts(6, mxgp), swd_sum,
     &     strain_sum(6), elem_ave_swd, elem_ave_strain(6),
     &     enode_swd(mxndel), enode_strains(6,mxndel), rnum_gpts,
     &     rcount
c
      double precision,
     &    dimension(:), pointer :: urcs_n, eps_n
c
      logical debug, chk_killed
      data one, zero /1.0d0, 0.0d0/
c
      debug = .false.
c
c                build average nodal values of the stress work density
c                (swd) and strains in the model. for each element,
c                shape functions will be used to extrapolate integration-
c                point values to each of the element's nodes where they
c                are then averaged.
c
c                loop over all structure elements. if element is involved
c                in J computations, call element dependent routine. handle
c                only elements whose data is owned by this processor.
c
      do elemno = 1, noelem
         blk = elems_to_blocks(elemno,1)
         if( elblks(2,blk) .ne. myid ) then
            if( debug ) then
               write(out,*) 'element ', elemno, 'skipped: not myid'
            end if
            cycle
         end if
         if( chk_killed( elemno ) ) then
            skipped_killed = skipped_killed + 1
            cycle
         end if
c
c                get element properties
c
         span             = 1
         etype            = iprops(1,elemno)
         num_enodes       = iprops(2,elemno)
         int_order        = iprops(5,elemno)
         num_gpts         = iprops(6,elemno)
      rnum_gpts       = dble( num_gpts )
         incptr           = incmap(elemno)
         rel_elem         = elems_to_blocks(elemno,2)
c
c                gather gauss point swd and strains for current element.
c
         swd_at_gpts(1:num_gpts)        = zero
         strain_at_gpts(1:6,1:num_gpts) = zero
c
         urcs_n => urcs_n_blocks(blk)%ptr
         eps_n => eps_n_blocks(blk)%ptr
         sig_offset = (rel_elem - 1) * nstrs * num_gpts
         eps_offset = (rel_elem - 1) * nstr  * num_gpts
c
         do gpn = 1, num_gpts
            swd_at_gpts(gpn)      = urcs_n(sig_offset + 7)
            strain_at_gpts(1,gpn) = eps_n(eps_offset + 1)
            strain_at_gpts(2,gpn) = eps_n(eps_offset + 2)
            strain_at_gpts(3,gpn) = eps_n(eps_offset + 3)
            strain_at_gpts(4,gpn) = eps_n(eps_offset + 4)
            strain_at_gpts(5,gpn) = eps_n(eps_offset + 5)
            strain_at_gpts(6,gpn) = eps_n(eps_offset + 6)
            sig_offset = sig_offset + nstrs
            eps_offset = eps_offset + nstr
         end do
c
         if( debug ) then
            write(out,100) elemno
            do gpt = 1,num_gpts
               write(out,200) gpt, swd_at_gpts(gpt),
     &              (strain_at_gpts(i,gpt), i=1,6)
            end do
         end if
c
c                extrapolate Gauss point values of swd and strain
c                for the current element to the element nodes.
c
         enode_swd(1:num_enodes)         = zero
         enode_strains(1:6,1:num_enodes) = zero
c
         call di_extrap_to_nodes( nonode, span, elemno, etype,
     &        num_enodes, int_order, num_gpts, extrap_counts,
     &        swd_at_gpts, enode_swd, strain_at_gpts, enode_strains,
     &        out, mxvl, mxoupr, mxndel, mxgp )
c
c                add element node values to structure node-value
c                arrays.
c
         do enode = 1, num_enodes
            snode                    = incid(incptr + enode - 1)
            extrap_counts(snode)     = extrap_counts(snode) + 1
            swd_at_nodes(snode)      = swd_at_nodes(snode)
     &                               + enode_swd(enode)
            strain_at_nodes(1,snode) = strain_at_nodes(1,snode)
     &                               + enode_strains(1,enode)
            strain_at_nodes(2,snode) = strain_at_nodes(2,snode)
     &                               + enode_strains(2,enode)
            strain_at_nodes(3,snode) = strain_at_nodes(3,snode)
     &                               + enode_strains(3,enode)
            strain_at_nodes(4,snode) = strain_at_nodes(4,snode)
     &                               + enode_strains(4,enode)
            strain_at_nodes(5,snode) = strain_at_nodes(5,snode)
     &                               + enode_strains(5,enode)
            strain_at_nodes(6,snode) = strain_at_nodes(6,snode)
     &                               + enode_strains(6,enode)
         end do
      end do
c
c                average the values at structure nodes.
c
      do snode = 1, nonode
         if ( extrap_counts(snode) .ne. 0 ) then
        rcount                  = one / dble( extrap_counts(snode) )
            swd_at_nodes(snode)      = swd_at_nodes(snode)      * rcount
            strain_at_nodes(1,snode) = strain_at_nodes(1,snode) * rcount
            strain_at_nodes(2,snode) = strain_at_nodes(2,snode) * rcount
            strain_at_nodes(3,snode) = strain_at_nodes(3,snode) * rcount
            strain_at_nodes(4,snode) = strain_at_nodes(4,snode) * rcount
            strain_at_nodes(5,snode) = strain_at_nodes(5,snode) * rcount
            strain_at_nodes(6,snode) = strain_at_nodes(6,snode) * rcount
         end if
      end do
c
      return
 100  format(/,'***before: gauss pt values for elem ',i4,/,'gpt',5x,
     &       'swd',15x,'strains')
 200  format(i3,1x,f10.5,1x,6(1x,f10.5))
      end
c
c ****************************************************************
c *                                                              *
c *     subroutine to extrapolate strain and stress work density *
c *     values from integration points to nodes for one element. *
c *                                                              *
c *                   written by: mcw                            *
c *                last modified: 02/01 by mcw                   *
c *                                                              *
c ****************************************************************
c
      subroutine di_extrap_to_nodes( nonode, span, elemno, etype,
     &     num_enodes, int_order, num_gpts, extrap_counts, swd_at_gpts,
     &     enode_swd, strain_at_gpts, enode_strains, out,
     &     mxvl, mxoupr, mxndel, mxgp )
c
      implicit none
c
c              dummy variables
c
      integer nonode, span, elemno, etype, num_enodes, int_order,
     &     num_gpts, extrap_counts(*), out, mxvl, mxoupr, mxndel, mxgp
c
      double precision
     &     swd_at_gpts(*), enode_swd(*), strain_at_gpts(6,*),
     &     enode_strains(6,*)
c
c              local variables
c
      integer val, gpn, enode, i, j, k, eps, node, idummy_vec(1),
     &        value
c
      double precision
     &     lg(mxgp), rdummy_vec(1), rnum_gpts, node_value,
     &     xi, eta, zeta, zero, one, one_by_ngp
c
      logical lagrangian_extrap, debug
      data one, zero /1.0d0, 0.0d0/
c
      debug = .false.
      enode_strains(1:6,1:num_enodes) = zero
      enode_swd(1:num_enodes)         = zero
c
c                       loop over element nodes. for 8 and 20-noded
c                       hex elements using 2x2x2 integration, we
c                       extrapolate to the nodes using the lagrangian
c                       polynomials. otherwise we just average the
c                       integration-point values and use that value
c                       at every node. warp3d uses the following
c                       arguments to describe elements and integration
c                       order:
c
c                       etype = 1: 20-noded brick
c                       etype = 2:  8-noded brick
c                       etype = 3: 12-noded brick
c                       etype = 4: 15-noded brick
c                       etype = 5:  9-noded brick
c
c                       etype = 1, int_order = 1: 27 point rule (not used)
c                       etype = 1, int_order = 8:  8 point rule
c                       etype = 1, int_order = 9: 14 point rule
c                       etype = 2, int_order = 1:  8 point rule
c                       etype = 2, int_order = 2:  6 point rule (not used)

c
      lagrangian_extrap = .false.
      if( etype.eq.2 ) lagrangian_extrap = .true.
      if( int_order.eq.1 .or. int_order.eq.8 )
     &         lagrangian_extrap = .true.
c
      if ( lagrangian_extrap ) go to 1000
c
c                       to get nodal values for the elements, extrapolation
c                       from integration points to node points is
c                       not possible. instead we use the average of all
c                       integration points at each node.
c
      do gpn = 1,num_gpts
         enode_strains(1,1) = enode_strains(1,1) + strain_at_gpts(1,gpn)
         enode_strains(2,1) = enode_strains(2,1) + strain_at_gpts(2,gpn)
         enode_strains(3,1) = enode_strains(3,1) + strain_at_gpts(3,gpn)
         enode_strains(4,1) = enode_strains(4,1) + strain_at_gpts(4,gpn)
         enode_strains(5,1) = enode_strains(5,1) + strain_at_gpts(5,gpn)
         enode_strains(6,1) = enode_strains(6,1) + strain_at_gpts(6,gpn)
         enode_swd(1)       = enode_swd(1)       + swd_at_gpts(gpn)
      end do
c
      if( debug ) then
         write(out,100)
         write(out,200) gpn, (enode_strains(i,1),i=1,6),
     &                  enode_swd(1)
      end if
c
c                       assign element integration point averages
c                       to the nodes.
c
      rnum_gpts  = num_gpts
      one_by_ngp = one / rnum_gpts
c
      do enode = 1, num_enodes
         enode_strains(1,enode) = enode_strains(1,1) * one_by_ngp
         enode_strains(2,enode) = enode_strains(2,1) * one_by_ngp
         enode_strains(3,enode) = enode_strains(3,1) * one_by_ngp
         enode_strains(4,enode) = enode_strains(4,1) * one_by_ngp
         enode_strains(5,enode) = enode_strains(5,1) * one_by_ngp
         enode_strains(6,enode) = enode_strains(6,1) * one_by_ngp
         enode_swd(enode)       = enode_swd(1)       * one_by_ngp
      end do
c
      if( debug) then
         write(out,300) elemno
         do enode = 1, num_enodes
            write(out,400) enode, enode_swd(enode),
     &           (enode_strains(i,enode), i=1,6)
         end do
      end if
      return
c
 1000 continue
c
c                       extrapolate integration point values to
c                       nodes using lagrangian shape functions.
c
      do enode = 1,num_enodes
c
c                       find the lagrangian shape functions for
c                       each gauss point at the current element node.
c
         call ndpts1( idummy_vec, 0, rdummy_vec, etype, enode,
     &                xi, eta, zeta )
         call oulgf( etype, xi, eta, zeta, lg, int_order )
c
c                       extrapolate strains and swd from integration
c                       points to current element node.
c
         do gpn = 1, num_gpts
            enode_strains(1,enode) = enode_strains(1,enode)
     &                             + strain_at_gpts(1,gpn) * lg(gpn)
            enode_strains(2,enode) = enode_strains(2,enode)
     &                             + strain_at_gpts(2,gpn) * lg(gpn)
            enode_strains(3,enode) = enode_strains(3,enode)
     &                             + strain_at_gpts(3,gpn) * lg(gpn)
            enode_strains(4,enode) = enode_strains(4,enode)
     &                             + strain_at_gpts(4,gpn) * lg(gpn)
            enode_strains(5,enode) = enode_strains(5,enode)
     &                             + strain_at_gpts(5,gpn) * lg(gpn)
            enode_strains(6,enode) = enode_strains(6,enode)
     &                             + strain_at_gpts(6,gpn) * lg(gpn)
            enode_swd(enode)       = enode_swd(enode)
     &                             + swd_at_gpts(gpn)      * lg(gpn)
         end do
c
      end do
c
      if( debug ) then
         write(out,300) elemno
         do enode = 1,num_enodes
            write(out,400) enode, enode_swd(enode),
     &           (enode_strains(i,enode),i=1,6)
         end do
      end if
c
      return
c
 100  format(/,'***summation of values at gpts')
 200  format(i5,7f10.5)
 300  format(/,'***node values from extrap routine',
     & ' elem ',i5,/,'enode ',3x,'swd',8x,'strains')
 400  format(i3,1x,f10.5,1x,6(1x,f10.5))
 500  format(/,'***extrapolated nodal values from extrap routine',/,
     &'enode ',15x,'swd',15x,'strains')
c
      end
c
c****************************************************************
c                                                               *
c      subroutine to write j and i-integral data to             *
c      standard output                                          *
c                    written by: mcw                            *
c                 last modified: 12/02 by mcw                   *
c                                                               *
c****************************************************************
c
      subroutine di_write_std_out( ltmstp, stname, lsldnm, out )
      use j_data
      implicit none
c
c             dummy variables
c
      integer ltmstp, out
      character(len=8) :: stname, lsldnm
c
c             local variables
c
      integer i, j, ring, skipped_killed
      double precision
     & rg_count
c
      rg_count = dble(ring_count)
c
      write( out,9000 )
c
c             write j-integral results to standard output
c
      if( comput_j ) then
         write(out,9005)
         write(out,9010)
         do i=1,ring_count
            ring           = int( j_storage(i,1) )
            skipped_killed = int( j_storage(i,11) )
            write(out,9015) ring, (j_storage(i,j),j=2,10),skipped_killed
         end do
c
c             print average, min, max of static and dynamic
c             domain values for j-integral computations.
c
         if( ring_count .gt. 1 ) then
            write(out,9020) domain_avg_j / rg_count, domain_min_j,
     &           domain_max_j
            if ( .not. static_j ) then
               write(out,9025)
               write(out,9020) static_avg / rg_count, static_min,
     &              static_max
            end if
         end if
         write(out,9030)
         if ( symmetric_domain ) write(out,9035)
         if ( rings_given ) write(out,9040)
c
c             print stress intensity factors calculated from J assuming
c             either pure mode-I, pure mode-II, or pure mode-III loading.
c
         if( j_to_k ) then
            write(out,9041)
            write(out,9042)
            write(out,9043)
            write(out,9044)
            do i=1,ring_count
               ring = int( j_storage(i,1) )
               write(out,9045) ring, ks_from_j(i,1), ks_from_j(i,2),
     &                         ks_from_j(i,1), ks_from_j(i,2),
     &                         ks_from_j(i,3)
            end do

         end if
c
      end if
c
c             write I-integral data to standard output
c
      if( comput_i ) then
         write(out,9050)
         if( symmetric_domain ) write(out,9051)
         if( process_temperatures ) write(out,9053)
c
c             KI plane stress, plane strain
c
         if( out_pstress ) then
            write(out,9055)
            do i=1,ring_count
               ring           = int( i_storage(i,1,1) )
               skipped_killed = int( i_storage(i,13,1) )
               write(out,9060) ring, (i_storage(i,j,1),j=2,11),
     &            skipped_killed
            end do
            if( ring_count .gt. 1 ) then
               write(out,9020) domain_avg_i(1) / rg_count,
     &              domain_min_i(1), domain_max_i(1)
            end if
         end if
c
         if( out_pstrain ) then
            write(out,9075)
            do i=1,ring_count
               ring           = int( i_storage(i,1,2) )
               skipped_killed = int( i_storage(i,13,2) )
               write(out,9060) ring, (i_storage(i,j,2),j=2,11),
     &            skipped_killed
            end do
            if( ring_count .gt. 1 ) then
               write(out,9020) domain_avg_i(2) / rg_count,
     &              domain_min_i(2), domain_max_i(2)
            end if
         end if
c
c             if symmetric option was used, assume mode-I loading
c             and skip output for KII, KIII.
c
         if( symmetric_domain ) goto 111
c
c             KII plane stress, plane strain
c
         if( out_pstress ) then
            write(out,9085)
            do i=1,ring_count
               ring           = int( i_storage(i,1,3) )
               skipped_killed = int( i_storage(i,13,3) )
               write(out,9060) ring, (i_storage(i,j,3),j=2,11),
     &            skipped_killed
            end do
            if( ring_count .gt. 1 ) then
               write(out,9020) domain_avg_i(3) / rg_count,
     &              domain_min_i(3), domain_max_i(3)
            end if
         end if
c
         if( out_pstrain ) then
            write(out,9095)
            do i=1,ring_count
               ring           = int( i_storage(i,1,4) )
               skipped_killed = int( i_storage(i,13,4) )
               write(out,9060) ring, (i_storage(i,j,4),j=2,11),
     &            skipped_killed
            end do
            if( ring_count .gt. 1 ) then
               write(out,9020) domain_avg_i(4) / rg_count,
     &              domain_min_i(4), domain_max_i(4)
            end if
         end if
c
c             KIII
c
         write(out,9105)
         do i=1,ring_count
            ring           = int( i_storage(i,1,5) )
            skipped_killed = int( i_storage(i,13,5) )
            write(out,9060) ring, (i_storage(i,j,5),j=2,11),
     &           skipped_killed
         end do
         if( ring_count .gt. 1 ) then
            write(out,9020) domain_avg_i(5) / rg_count,
     &           domain_min_i(5), domain_max_i(5)
         end if
c
 111     continue
c
c             print J-values calculated from the mixed-mode
c             stress intensity factors computed using the
c             interaction integral.
c
         write(out,9200)
         if( out_pstress .and. .not. out_pstrain ) then
            write(out,9202)
            write(out,9210)
            do i=1,ring_count
               ring = int( i_storage(i,1,1) )
               write(out,9230) ring, j_from_ks(i,1)
            end do
         end if
         if( out_pstrain .and. .not. out_pstress ) then
            write(out,9204)
            write(out,9215)
            do i=1,ring_count
               ring = int( i_storage(i,1,1) )
               write(out,9230) ring, j_from_ks(i,2)
            end do
         end if
         if( out_pstress .and. out_pstrain ) then
            write(out,9206)
            write(out,9220)
            do i=1,ring_count
               ring = int( i_storage(i,1,1) )
               write(out,9240) ring, j_from_ks(i,1), j_from_ks(i,2)
            end do
         end if
         write(out,9260)
c
c             T-stress T11 & T33, plane stress and plane strain
c
         if( face_loading ) write(out,9054)
c
         if( out_pstress ) then
            write(out,9107)
            do i=1,ring_count
               ring           = int( i_storage(i,1,6) )
               skipped_killed = int( i_storage(i,13,6) )
               write(out,9065) ring, (i_storage(i,j,6),j=2,12),
     &            skipped_killed
            end do
            if( ring_count .gt. 1 ) then
               write(out,9115) domain_avg_i(6) / rg_count,
     &              domain_min_i(6), domain_max_i(6),
     &              domain_avg_i(9) / rg_count, domain_min_i(9),
     &              domain_max_i(9)
            end if
         end if
c
         if( out_pstrain ) then
            write(out,9110)
            do i=1,ring_count
               ring           = int( i_storage(i,1,7) )
               skipped_killed = int( i_storage(i,13,7) )
               write(out,9065) ring, (i_storage(i,j,7),j=2,12),
     &            skipped_killed
            end do
            if( ring_count .gt. 1 ) then
               write(out,9115) domain_avg_i(7) / rg_count,
     &              domain_min_i(7), domain_max_i(7),
     &              domain_avg_i(10) / rg_count, domain_min_i(10),
     &              domain_max_i(10)
            end if
         end if
c
c             if symmetric option was used, assume mode-I loading
c             and skip output for T13.
c
         if( symmetric_domain ) goto 222
c
c             T-stress T13, anti-plane shear
c
         write(out,9112)
         do i=1,ring_count
            ring           = int( i_storage(i,1,8) )
            skipped_killed = int( i_storage(i,13,8) )
            write(out,9060) ring, (i_storage(i,j,8),j=2,11),
     &           skipped_killed
         end do
         if( ring_count .gt. 1 ) then
            write(out,9020) domain_avg_i(8) / rg_count,
     &           domain_min_i(8), domain_max_i(8)
         end if
c
 222     continue
c
         write(out,9125)
         if ( symmetric_domain ) write(out,9130)
         if ( rings_given ) write(out,9040)
      end if
c
      return
c
 9000 format(///,4x,'---- TABULARIZED RESULTS ----',/)
 9005 format(/,15x,'Total J-values')
 9010 format(/,1x,30x,'J-integral components',
     &  /,      1x,27x,'---------------------------',
     & //,1x,'domain',8x,'dm1',10x,'dm2',10x,'dm3',10x,'dm4',
     &    10x,'dm5',10x,'dm6',10x,'dm7',10x,'dm8',7x,'total J',3x,
     &    'killed ele' )
 9015 format(/,1x,i5,2x,9(2x,e11.4),3x,'(',i3,')')
 9020 format(/,3x,' average  ',3x,' minimum  ',3x,' maximum',
     &       /,1x, e11.4,2x,e11.4,2x,e11.4 )
 9025 format(/,10x,'Static J-values' )
 9030 format(/,1x,'* dm1: stress work density term',
     &       /,1x,'* dm2: traction-displacement derivative term',
     &       /,1x,'* dm3: kinetic energy term',
     &       /,1x,'* dm4: acceleration term',
     &       /,1x,'* dm5: crack face loading term',
     &       /,1x,'* dm6: thermal loading term (zero for an fgm)',
     &       /,1x,'* dm7: stress times derivative of strain',
     &         1x, '(zero if not an fgm)'
     &       /,1x,'* dm8: derivative of stress work density',
     &         1x, '(zero if not an fgm)'
     &       /,1x,'* ( ): number of killed elements skipped' )
 9035 format(/,1x,'* note: J-values doubled for symmetry *')
 9040 format(/,1x,'* note: automatic domain builder cannot detect',
     &       /,1x,'        domains that reach a boundary',/)
 9041 format(//,10x,'Stress intensity factors obtained from J-values',
     &       ' assuming pure mode I, mode II or mode III loading')
 9042 format(/,16x,'pure mode I loading',17x,'pure mode II loading',12x,
     &       'pure mode III loading')
 9043 format(10x,'---------------------------------',3x,
     &       '-----------------------------------',3x,
     &       '---------------------')
 9044 format(1x,'domain',3x,'KI plane stress',3x,'KI plane strain',
     &       3x,'KII plane stress',3x,'KII plane strain',3x,
     &       'KIII anti-plane shear' )
 9045 format(/,1x,i5,5x,e11.4,7x,e11.4,8x,e11.4,8x,e11.4,10x,e11.4)
 9050 format(///,15x,'Static I-values')
 9051 format(//,1x,'>>>>> warning: ''symmetric'' option detected.',
     &       /,16x,'Values for KII, KIII and T13 set to zero.',
     &       ' results valid for mode-I loading only.' )
 9053 format(//,1x,'>>>>> warning: thermal loading detected.',
     &       /,16x,'Interaction Integral results below for K and T'
     &       ' do not currently include thermal effects.')
 9054 format(//,1x,'>>>>> warning: crack-face loading detected.',
     &       /,16x,'Interaction Integral results below for T'
     &       ' do not include crack-face loading effects.')
 9055 format(/,1x,30x,'I-integral components',
     &       /,26x,'-- KI = 1, KII = 0, KIII = 0, plane stress --',
     &       //,1x,'domain',8x,'dm1',10x,'dm2',10x,'dm3',10x,'dm4',
     &       10x,'dm5',10x,'dm6',10x,'dm7',10x,'dm8',7x,'total I',5x,
     &       'KI (pstrs)',1x,'killed ele' )
 9060 format(/,1x,i5,2x,10(2x,e11.4),3x,'(',i3,')')
 9065 format(/,1x,i5,2x,11(2x,e11.4),3x,'(',i3,')')
 9075 format(/,1x,30x,'I-integral components',
     &       /,26x,'-- KI = 1, KII = 0, KIII = 0, plane strain --',
     &       //,1x,'domain',8x,'dm1',10x,'dm2',10x,'dm3',10x,'dm4',
     &       10x,'dm5',10x,'dm6',10x,'dm7',10x,'dm8',7x,'total I',5x,
     &       'KI (pstrn)',1x,'killed ele' )
 9085 format(/,1x,30x,'I-integral components',
     &       /,26x,'-- KI = 0, KII = 1, KIII = 0, plane stress --',
     &       //,1x,'domain',8x,'dm1',10x,'dm2',10x,'dm3',10x,'dm4',
     &       10x,'dm5',10x,'dm6',10x,'dm7',10x,'dm8',7x,'total I',5x,
     &       'KII (pstrs)',1x,'killed ele' )
 9095 format(/,1x,30x,'I-integral components',
     &       /,26x,'-- KI = 0, KII = 1, KIII = 0, plane strain --',
     &       //,1x,'domain',8x,'dm1',10x,'dm2',10x,'dm3',10x,'dm4',
     &       10x,'dm5',10x,'dm6',10x,'dm7',10x,'dm8',7x,'total I',5x,
     &       'KII (pstrn)',1x,'killed ele' )
 9105 format(/,1x,30x,'I-integral components',
     &       /,26x,'-- KI = 0, KII = 0, KIII = 1 --',
     &       //,1x,'domain',8x,'dm1',10x,'dm2',10x,'dm3',10x,'dm4',
     &       10x,'dm5',10x,'dm6',10x,'dm7',10x,'dm8',7x,'total I',8x,
     &       'KIII',4x,'killed ele' )
 9200 format(///,10x,'J-values computed from mixed-mode stress',
     &       ' intensity factors KI, KII and KIII.')
 9202 format(10x,'KI- and KII-values correspond to plane stress',
     &       ' conditions.')
 9204 format(10x,'KI- and KII-values correspond to plane strain',
     &       ' conditions.')
 9206 format(10x,'KI- and KII-values correspond to either',
     &       ' plane stress or plane strain conditions.')
 9210 format(/,1x,'domain',3x,'J plane stress')
 9215 format(/,1x,'domain',3x,'J plane strain')
 9220 format(/,1x,'domain',3x,'J plane stress',3x,'J plane strain')
 9230 format(/,1x,i5,5x,e11.4)
 9240 format(/,1x,i5,5x,e11.4,6x,e11.4)
 9260 format(/)
 9107 format(/,1x,30x,'I-integral components',
     &       /,26x,'--T-stress (plane stress auxiliary fields)--',
     &       //,1x,'domain',8x,'dm1',10x,'dm2',10x,'dm3',10x,'dm4',
     &       10x,'dm5',10x,'dm6',10x,'dm7',10x,'dm8',7x,'total I',4x,
     &       'T11 (pstrs)',6x,'T33',5x,'killed ele' )
 9110 format(/,1x,30x,'I-integral components',
     &       /,26x,'--T-stress ("plane strain" auxiliary fields)--',
     &       //,1x,'domain',8x,'dm1',10x,'dm2',10x,'dm3',10x,'dm4',
     &       10x,'dm5',10x,'dm6',10x,'dm7',10x,'dm8',7x,'total I',4x,
     &       'T11 (pstrn)',6x,'T33',5x,'killed ele' )
 9112 format(/,1x,30x,'I-integral components',
     &       /,26x,'--T-stress (anti-plane shear auxiliary fields)--',
     &       //,1x,'domain',8x,'dm1',10x,'dm2',10x,'dm3',10x,'dm4',
     &       10x,'dm5',10x,'dm6',10x,'dm7',10x,'dm8',7x,'total I',8x,
     &       'T13',5x,'killed ele' )
 9115 format(/,3x,' average  ',3x,' minimum  ',3x,' maximum',
     &       /,1x, 3(e11.4,2x),
     &       /,1x, 3(e11.4,2x))
 9125 format(/,1x,'* dm1: stress * deriv of aux. displ.',
     &       /,1x,'* dm2: aux stress * derivative of displacement',
     &       /,1x,'* dm3: mixed strain energy density',
     &       ' (aux stress * strain)',
     &       /,1x,'* dm4: stress * 2nd deriv of aux displ',
     &       /,1x,'* dm5: stress * deriv of aux strain)',
     &       /,1x,'* dm6: deriv of constitutive tensor * strain * aux',
     &       'strain',
     &       /,1x,'* dm7: thermal strains (not yet implemented)',
     &       /,1x,'* dm8: traction * deriv of aux displ',
     &       /,1x,'* ( ): number of killed elements skipped' )
 9130 format(/,1x,'* note: I-values doubled for symmetry *')
c
      end
c
c
c****************************************************************
c                                                               *
c      subroutine to write J and I-integral data to packets     *
c                                                               *
c                    written by: mcw                            *
c                 last modified: 12/17 by mcw                   *
c                                                               *
c****************************************************************
c
      subroutine di_write_packets( ltmstp, stname, lsldnm )
      use main_data, only : output_packets, packet_file_no
      use j_data
      implicit none
c
c             dummy variables
c
      integer ltmstp
      character(len=8) :: stname, lsldnm
c
c             local variables
c
      integer i, j, ring, skipped_killed, num_lines
      double precision
     & rg_count
c
      rg_count = dble(ring_count)
c
c             write J-integral data to packet 17
c
      if( output_packets .and. comput_j .and. output_packet_j ) then
         if( ring_count.eq.1 ) then
            num_lines = 3
            write( packet_file_no ) 17, num_lines, ltmstp, 0
            write( packet_file_no ) 1, domain_id, stname, lsldnm, ltmstp
         else
            num_lines = 1 + ring_count + 1 + ring_count
            write( packet_file_no ) 17, num_lines, ltmstp, 0
            write( packet_file_no ) ring_count, domain_id, stname,
     &           lsldnm, ltmstp
         end if
         do i=1,ring_count
            ring           = int( j_storage(i,1) )
            skipped_killed = int( j_storage(i,11) )
            write( packet_file_no ) ring, (j_storage(i,j),j=2,10),
     &                              skipped_killed
         end do
         if( ring_count.gt.1 ) then
            write( packet_file_no ) domain_avg_j/rg_count, domain_min_j,
     &           domain_max_j
         end if
c
c             write stress intensity factors computed from J
c
         do i=1,ring_count
            ring = int( j_storage(i,1) )
            write( packet_file_no ) ring, ks_from_j(i,1),
     &           ks_from_j(i,2), ks_from_j(i,1), ks_from_j(i,2),
     &           ks_from_j(i,3)
         end do
      end if
c
c             write I-integral stress intensity factor data to packet 27
c
      if( output_packets .and. comput_i .and. output_packet_i ) then
         if( ring_count.eq.1 ) then
            num_lines = 7
            write( packet_file_no ) 27, num_lines, ltmstp, 0
            write( packet_file_no ) 1, domain_id, stname, lsldnm, ltmstp
         else
            num_lines = 1 + 5 * ( ring_count + 1 ) + ring_count
            write( packet_file_no ) 27, num_lines, ltmstp,0
            write( packet_file_no ) ring_count, domain_id, stname,
     &           lsldnm, ltmstp
         end if
c
c             plane stress, plane strain KI
c
         do i=1,ring_count
            ring           = int( i_storage(i,1,1) )
            skipped_killed = int( i_storage(i,13,1) )
            write( packet_file_no ) ring, (i_storage(i,j,1),j=2,11),
     &                              skipped_killed
         end do
         if( ring_count .gt. 1 ) then
            write( packet_file_no ) domain_avg_i(1) / rg_count,
     &           domain_min_i(1), domain_max_i(1)
         end if
c
         do i=1,ring_count
            ring           = int( i_storage(i,1,2) )
            skipped_killed = int( i_storage(i,13,2) )
            write( packet_file_no ) ring, (i_storage(i,j,2),j=2,11),
     &                              skipped_killed
         end do
         if( ring_count .gt. 1 ) then
            write( packet_file_no ) domain_avg_i(2) / rg_count,
     &           domain_min_i(2), domain_max_i(2)
         end if
c
c             plane stress, plane strain KII
c
         do i=1,ring_count
            ring           = int( i_storage(i,1,3) )
            skipped_killed = int( i_storage(i,13,3) )
            write( packet_file_no ) ring, (i_storage(i,j,3),j=2,11),
     &                              skipped_killed
         end do
         if( ring_count .gt. 1 ) then
            write( packet_file_no ) domain_avg_i(3) / rg_count,
     &           domain_min_i(3), domain_max_i(3)
         end if
c
         do i=1,ring_count
            ring           = int( i_storage(i,1,4) )
            skipped_killed = int( i_storage(i,13,4) )
            write( packet_file_no ) ring, (i_storage(i,j,4),j=2,11),
     &                              skipped_killed
         end do
         if( ring_count .gt. 1 ) then
            write( packet_file_no ) domain_avg_i(4) / rg_count,
     &           domain_min_i(4), domain_max_i(4)
         end if
c
c             KIII
c
         do i=1,ring_count
            ring           = int( i_storage(i,1,5) )
            skipped_killed = int( i_storage(i,13,5) )
            write( packet_file_no ) ring, (i_storage(i,j,5),j=2,11),
     &                              skipped_killed
         end do
         if( ring_count .gt. 1 ) then
            write( packet_file_no ) domain_avg_i(5) / rg_count,
     &           domain_min_i(5), domain_max_i(5)
         end if
c
c             J-values computed from K-values
c
         do i=1,ring_count
            ring = int( i_storage(i,1,1) )
            write( packet_file_no ) ring, j_from_ks(i,1), j_from_ks(i,2)
         end do
c
      end if
c
c             write I-integral T-stress data to packet 28
c
      if( output_packets .and. comput_i .and. output_packet_i ) then
         if( ring_count.eq.1 ) then
            num_lines = 4
            write( packet_file_no ) 28, num_lines, ltmstp, 0
            write( packet_file_no ) 1, domain_id, stname, lsldnm, ltmstp
         else
            num_lines = 1 + 3 * ( ring_count + 1 )
            write( packet_file_no ) 28, num_lines, ltmstp,0
            write( packet_file_no ) ring_count, domain_id, stname,
     &           lsldnm, ltmstp
         end if
c
c             T-stresses T11 and T33 for plane stress and "plane-strain"
c             auxiliary fields
c
         do i=1,ring_count
            ring           = int( i_storage(i,1,6) )
            skipped_killed = int( i_storage(i,13,6) )
            write( packet_file_no ) ring, (i_storage(i,j,6),j=2,12),
     &                              skipped_killed
         end do
         if( ring_count .gt. 1 ) then
            write( packet_file_no ) domain_avg_i(6) / rg_count,
     &           domain_min_i(6), domain_max_i(6),
     &           domain_avg_i(9) / rg_count,
     &           domain_min_i(9), domain_max_i(9)
         end if
c
         do i=1,ring_count
            ring           = int( i_storage(i,1,7)  )
            skipped_killed = int( i_storage(i,13,7) )
            write( packet_file_no ) ring, (i_storage(i,j,7),j=2,12),
     &                              skipped_killed
         end do
         if( ring_count .gt. 1 ) then
            write( packet_file_no ) domain_avg_i(7) / rg_count,
     &           domain_min_i(7), domain_max_i(7),
     &           domain_avg_i(10) / rg_count,
     &           domain_min_i(10), domain_max_i(10)
         end if
      end if
c
c             T-stress T13 for anti-plane shear auxiliary fields
c
      if( output_packets .and. comput_i .and. output_packet_i ) then
         do i=1,ring_count
            ring           = int( i_storage(i,1,8) )
            skipped_killed = int( i_storage(i,13,8) )
            write( packet_file_no ) ring, (i_storage(i,j,8),j=2,11),
     &                              skipped_killed
         end do
         if( ring_count .gt. 1 ) then
            write( packet_file_no ) domain_avg_i(8) / rg_count,
     &           domain_min_i(8), domain_max_i(8)
         end if
      end if
c
      return
 9060 format(/,1x,i5,2x,10(2x,e11.4),3x,'(',i3,')')
      end
c
c *******************************************************************
c *                                                                 *
c *   calculate strain e33 at domain origin for T-stress calcs.     *
c *   calculate strain e33 as the difference between the            *
c *   deformed and undeformed crack-front lengths delta_L / L       *
c *                                                                 *
c *                                    written by: mcw              *
c *                                 last modified: 12/20/12/ rhd    *
c *                                                                 *
c *******************************************************************
c
c
      subroutine di_calc_e33( out)
      use j_data
      implicit none
c
c             parameter declarations
c
      integer out
c
c             local declarations
c
      integer i, j, nfn, ngp, nfelem, inc, elem, gp, node
      double precision
     &  sf(4), dsf(4), jacob, jacobi, xsi, displ_coords(3,30),
     &  length_undisp, length_disp, j1, j2, j3, detj, dieldp, w,
     &  dummy, zero, half, one, two, x1, z1, x2, z2, x3, z3
      logical debug
c
      data zero, half, one, two / 0.0d0, 0.5d0, 1.0d0, 2.0d0 /
c
      debug = .false.
c
      if( debug ) write(out,9000) domain_type, num_front_nodes
c
c             generate array with local coordinates of displaced front nodes
c
      do i = 1, num_front_nodes
         displ_coords(1,i) = front_coords(1,i) + front_node_displ(1,i)
         displ_coords(2,i) = front_coords(2,i) + front_node_displ(2,i)
         displ_coords(3,i) = front_coords(3,i) + front_node_displ(3,i)
         if( debug ) then
            write(out,9049) 'front node ', i
            write(out,9050) 'undispl: ',front_coords(1,i),
     &                      front_coords(2,i), front_coords(3,i)
            write(out,9050) 'delta:   ',front_node_displ(1,i),
     &                      front_node_displ(2,i), front_node_displ(3,i)
            write(out,9050) 'displ:   ',displ_coords(1,i),
     &                      displ_coords(2,i), displ_coords(3,i)
            write(out,9051) ' '
         end if
      end do
c
c             nfelem -- number of elements along crack-front in domain
c             nfn    -- number of nodes along an element edge
c             ngp    -- number of integration points
c                           = 1 for linear elements (exact)
c                           = 4 for quadratic elements (approximate)
c             inc    -- increment to find first crack-front node of
c                       next element along crack front
c
      if( front_order.eq.1 ) then
         nfelem = num_front_nodes - 1
         nfn    = 2
         ngp    = 1
         inc    = 1
      else if( front_order.eq.2 ) then
         nfelem = (num_front_nodes - 1)/2
         nfn    = 3
         ngp    = 4
         inc    = 2
      else
         write(out,9100) front_order
         write(out,9150) domain_type
         call die_abort
         stop
      end if
      if( debug ) write(out,9200) nfn, ngp, inc
c
      length_undisp = zero
      length_disp   = zero
      node          = 1
c
      do elem = 1, nfelem
         do gp = 1, ngp
            call getgpts( 17, ngp, gp, xsi, dummy, dummy, w )
            call di1dsf( xsi, dsf, sf, nfn )
            j1            = dieldp(dsf(1),front_coords(1,node),nfn,1,3)
            j2            = dieldp(dsf(1),front_coords(2,node),nfn,1,3)
            j3            = dieldp(dsf(1),front_coords(3,node),nfn,1,3)
            detj          = sqrt( j1**2 + j2**2 + j3**2 )
            length_undisp = length_undisp + detj * w
c
            j1            = dieldp(dsf(1),displ_coords(1,node),nfn,1,3)
            j2            = dieldp(dsf(1),displ_coords(2,node),nfn,1,3)
            j3            = dieldp(dsf(1),displ_coords(3,node),nfn,1,3)
            detj          = sqrt( j1**2 + j2**2 + j3**2 )
            length_disp   = length_disp + detj * w
            if( debug ) write(out,9300) elem, gp, xsi, w
         end do
         node = node + inc
      end do
c
c             compute strain for front segment as (change in length) / length
c
      e33_front = (length_disp - length_undisp) / length_undisp
c
      if( debug ) write(out,9800) domain_type, num_front_nodes,
     &                            length_undisp, length_disp,
     &                            length_disp-length_undisp, e33_front
c
      if( comput_i .and. out_pstrain ) write(out,9900) e33_front
c
      return
c
 9000 format('>>>>> entered di_calc_e33',/,5x,'domain_type:     ',i2,
     &       /,5x,'num_front_nodes: ',i2)
 9049 format(a,i7)
 9050 format(a9,5x,3(2x,e26.16))
 9051 format(a)
 9100 format('>>>>> SYSTEM ERROR: in di_calc_e33 domain_type has',
     &       ' illegal value: ',i10,/,'          job terminated.')
 9150 format(5x,'domain_type: ',i2)
 9200 format(5x,'nfn: ',i2,2x,'ngp: ',i2,2x,'inc: ',i2)
 9300 format(5x,'elem: ',i8,2x,'point: ',i2,2x,'xsi: ',e11.4,2x,'w: ',
     &       e11.4)
 9800 format(/,'      > leaving di_calc_e33',
     & /,20x,'domain_type:           ',i10,
     & /,20x,'num_front_nodes:       ',i10,
     & /,20x,'length_undisp:         ',e26.16,
     & /,20x,'length_disp:           ',e26.16,
     & /,20x,'change in length:      ',e26.16,
     & /,20x,'e33_front:             ',e13.6)
 9900 format(//,' crack front tangential strain',
     &       /, ' epsilon_33 for computation of',
     &       /, ' ''plane strain'' T11, T33:        ',e13.6)
c
      end
c
c *******************************************************************
c *                                                                 *
c *   calculate coefficients of curve described by crack front      *
c *   nodes.                                                        *
c *                                                                 *
c *                                    written by: mcw              *
c *                                 last modified: 04/26/04         *
c *                                                                 *
c *******************************************************************
c
c
      subroutine di_calc_curvature( debug, out )
      use j_data
      implicit none
c
c             dummy variables
c
      integer out
      logical debug
c
c             local variables
c
      integer i, j
      double precision
     & zero, half, one, two, toler, x1, z1, x2, z2, x3, z3, x4, z4,
     & x5, z5, a, b, c, d, e, x_center, z_center, circle_radius, q
c
      data zero, half, one, two, toler
     & / 0.0d0, 0.5d0, 1.0d0, 2.0d0, 1.0d-04 /
c
c      debug = .true.
c
c             return if elements are not quadratic
c
      if( debug ) then
         write(out,*)
         write(out,*) "debug di_calc_curvature"
         write(out,*) "front_order        = ", front_order
         write(out,*) "num_front_nodes    = ", num_front_nodes
      end if

      if( front_order .ne. 2 ) return
c
c             determine if crack front is curved.
c             find coefficients of a quadratic expression defined by
c             three crack front nodes: x = az^2 + bz + c
c             if 'a' equals zero +/- tolerance, consider the front
c             nodes to be colinear, and exit.
c
      if( debug ) write(out,*) "@1"

      if( num_front_nodes .eq. 3 ) then
        if( debug ) write(out,*) "@2"

         x1 = front_coords(1,1)
         z1 = front_coords(3,1)
         x2 = front_coords(1,2)
         z2 = front_coords(3,2)
         x3 = front_coords(1,3)
         z3 = front_coords(3,3)
         q =   ( x3 - x2 ) / ( ( z3 - z2 ) * ( z3 - z1 ) )
     &       - ( x1 - x2 ) / ( ( z1 - z2 ) * ( z3 - z1 ) )
         if( abs(q) .gt. toler ) crack_curvature(1) = one
      end if
c
c             if num_front_nodes is greater than 5, we have a domain
c             type 'd'. in this case, estimate curvature based on
c             the first two elements along the crack in the domain.
c
      if( debug ) write(out,*) "@3"

      if( num_front_nodes .ge. 5 ) then
        if( debug ) write(out,*) "@4"

         x1 = front_coords(1,1)
         z1 = front_coords(3,1)
         x2 = front_coords(1,2)
         z2 = front_coords(3,2)
         x3 = front_coords(1,3)
         z3 = front_coords(3,3)
         x4 = front_coords(1,4)
         z4 = front_coords(3,4)
         x5 = front_coords(1,5)
         z5 = front_coords(3,5)
         q =   ( x3 - x2 ) / ( ( z3 - z2 ) * ( z3 - z1 ) )
     &       - ( x1 - x2 ) / ( ( z1 - z2 ) * ( z3 - z1 ) )
         if( abs(q) .gt. toler ) crack_curvature(1) = one
         if( debug ) write(out,*) "q first segment  = ", q
         q =   ( x5 - x4 ) / ( ( z5 - z4 ) * ( z5 - z3 ) )
     &       - ( x3 - x4 ) / ( ( z3 - z4 ) * ( z5 - z3 ) )
         if( abs(q) .gt. toler ) crack_curvature(1) = one
         if( debug ) write(out,*) "q second segment = ", q
      end if
c
c     12/19/2007 GVT: use the "nint" function to check if the value
c                     is integer zero; may help avoid floating point
c                     truncation
c
      if( debug ) write(out,*) "@5"

      if( nint(crack_curvature(1)) .eq. 0 ) go to 1111
c     BUG FIX BY GREG T. if( int(crack_curvature(1)) .eq. 0 ) go to 1111
c
c             contents of array 'crack_curvature':
c
c             crack_curvature(1) = 0 for straight front
c                                = 1 for curved front
c             crack_curvature(2) = 0 for unknown curvature
c                                = 1 for for a circular curve
c             crack_curvature(3) = local x coordinate for circle
c                                = coefficient 'a' for polynomial
c             crack_curvature(4) = local z coordinate for circle
c                                = coefficient 'b' for polynomial
c             crack_curvature(5) = circle radius for circle
c                                = coefficient 'c' for polynomial
c             crack_curvature(6) = coefficient 'd' for polynomial
c             crack_curvature(7) = coefficient 'e' for polynomial
c
c             for a circular crack, find center and radius of a
c             circle defined by the first three crack-front nodes.
c
      if( debug ) write(out,*) "@6"

      if( nint(crack_curvature(2)) .eq. 1 ) then
c     BUG FIX BY GREG T. if( int(crack_curvature(2)) .eq. 1 ) then
        if( debug ) write(out,*) "@7"

         x_center = (   ( x1*x1 + z1*z1 ) * ( z3 - z2 )
     &                + ( x2*x2 + z2*z2 ) * ( z1 - z3 )
     &                + ( x3*x3 + z3*z3 ) * ( z2 - z1 ) )
     &            / ( two * (   x1 * ( z3 - z2 )
     &                        + x2 * ( z1 - z3 )
     &                        + x3 * ( z2 - z1 ) ) )
         z_center = (   ( x1*x1 + z1*z1 ) * ( x3 - x2 )
     &                + ( x2*x2 + z2*z2 ) * ( x1 - x3 )
     &                + ( x3*x3 + z3*z3 ) * ( x2 - x1 ) )
     &            / ( two * (   z1 * ( x3 - x2 )
     &                        + z2 * ( x1 - x3 )
     &                        + z3 * ( x2 - x1 ) ) )
         circle_radius = (   (x1 - x_center)**two
     &                     + (z1 - z_center)**two )**half
c
         crack_curvature(3) = x_center
         crack_curvature(4) = z_center
         crack_curvature(5) = circle_radius
         go to 1111
      end if
c
c             for crack fronts of unknown curvature:
c
      if( debug ) write(out,*) "@8"

      if( nint( crack_curvature(2) ) .eq. 0 ) then
c     BUG FIX BY GREG T. if( int( crack_curvature(2) ) .eq. 0 ) then
        if( debug ) write(out,*) "@9"
c
c             for three nodes on front, use quadratic fit:
c             x = az^2 + bz + c
c
         if( num_front_nodes .eq. 3 ) then
            if( debug ) write(out,*) "@10"

            a =   ( x3 - x2 ) / ( ( z3 - z2 ) * ( z3 - z1 ) )
     &          - ( x1 - x2 ) / ( ( z1 - z2 ) * ( z3 - z1 ) )
            b = ( x1 - x2 + a * ( z2*z2 - z1*z1 ) ) / ( z1 - z2 )
            c =   x1 - a * z1*z1 - b * z1
            crack_curvature(3) = a
            crack_curvature(4) = b
            crack_curvature(5) = c
            go to 1111
         end if
c
c             for five nodes on front, fit nodes to a fourth-order
c             polynomial:  x = az^4 + bz^3 + cz^2 + dz + e.
c             for more than five nodes on front, use curvature
c             determined from first five nodes.
c
         if( debug ) write(out,*) "@11"

         if( num_front_nodes .ge. 5 ) then
           if( debug ) write(out,*) "@12"
c
c             call least-squares regression subroutine to compute
c             coefficients of polynomial.
c
            call di_calc_coefficients( x1, z1, x2, z2, x3, z3, x4, z4,
     &                               x5, z5, a, b, c, d, e, out, debug )
            crack_curvature(3) = a
            crack_curvature(4) = b
            crack_curvature(5) = c
            crack_curvature(6) = d
            crack_curvature(7) = e
         end if
      end if

      if( debug ) write(out,*) "@13"
c
 1111 continue
      if( debug ) then
         write(out,*) "curvature 'q'      = ", q
         write(out,*) "tolerance          = ", toler
         write(out,*) "crack_curvature(1) = ", crack_curvature(1)
         write(out,*) "crack_curvature(2) = ", crack_curvature(2)
         write(out,*) "crack_curvature(3) = ", crack_curvature(3)
         write(out,*) "crack_curvature(4) = ", crack_curvature(4)
         write(out,*) "crack_curvature(5) = ", crack_curvature(5)
         write(out,*) "crack_curvature(6) = ", crack_curvature(6)
         write(out,*) "crack_curvature(7) = ", crack_curvature(7)
      end if
c
      return
c
      end
c
c
c ********************************************************************
c *                                                                  *
c * subroutine to compute coefficients of a fourth-order             *
c * polynomial by a linear least-squares regression of               *
c * the five data points of the five crack-front nodes.              *
c *                                                                  *
c *                                written by: mcw                   *
c *                             last modified: 4/24/04               *
c *                                                                  *
c *    subroutine:  di_calc_coefficients                             *
c *                                                                  *
c *    purpose: solve a linear set of simultaneous equations using   *
c *             gauss elimination with partial pivoting. one rhs     *
c *             is permitted. the coefficient matrix may be non-symm *
c *             etric. the coefficient matrix is replaced by it's    *
c *             triangulated form.                                   *
c *             solves: [a] {x} = {b}                                *
c *                                                                  *
c *    written by:  r. h. dodds                                      *
c *    last revision: rhd 10/14/96 -- original coding                *
c *    parameters: (double precision unless noted)                   *
c *        a       -- coefficient matrix of                          *
c *        x       -- vector to receive computed solution            *
c *        b       -- vector containing the right hand side          *
c *        nrow_a  -- (integer) dimensioned number of rows for       *
c *                   matrix a                                       *
c *        neqns   -- (integer) number of linear equations           *
c *        iwork   -- (integer) work vector of length neqns          *
c *        dwork   -- work vector of length neqns                    *
c *        error   -- (logical) returned .true. if the equations     *
c *                    could not be solved due to a zero diagonal    *
c *                                                                  *
c *        the work vectors must be allocated by the calling program *
c *        for use by this routine to perform the pivoting           *
c *        operations.                                               *
c *                                                                  *
c *    local variables: (double precision unless noted)              *
c *        workmax, rmax, r, sum xmult -- temporaries                *
c *        i, j, k, ll, lk  -- (integer) indexes for loops           *
c *                                                                  *
c *    global variables: none                                        *
c *    routines called:  none                                        *
c *    intrinsics called: max, abs                                   *
c *                                                                  *
c ********************************************************************
c
      subroutine di_calc_coefficients( x1, z1, x2, z2, x3, z3, x4, z4,
     &                               x5, z5, aa, bb, cc, dd, ee,
     &                               out, debug )
      implicit double precision (a-h, o-z)
c
c             dummy arguments
c
      integer out
      double precision
     & x1, z1, x2, z2, x3, z3, x4, z4, x5, z5, aa, bb, cc, dd, ee
      logical debug
c
c             local arguments
c
      logical error
      integer nrow_a, neqns
      dimension a(5,5), x(5), b(5), iwork(5), dwork(5)
      data zero, half, one, two, tolerance
     & / 0.0d0, 0.5d0, 1.0d0, 2.0d0, 1.0d-20 /
c
      nrow_a   = 5
      neqns    = 5
c
c             fill the vector of known data points
c
      b(1)     = x1
      b(2)     = x2
      b(3)     = x3
      b(4)     = x4
      b(5)     = x5
c
c             build the Vandermonde matrix
c
      a(1,1)   = one
      a(2,1)   = one
      a(3,1)   = one
      a(4,1)   = one
      a(5,1)   = one
      a(1,2)   = z1
      a(2,2)   = z2
      a(3,2)   = z3
      a(4,2)   = z4
      a(5,2)   = z5
      a(1,3)   = z1*z1
      a(2,3)   = z2*z2
      a(3,3)   = z3*z3
      a(4,3)   = z4*z4
      a(5,3)   = z5*z5
      a(1,4)   = z1*z1*z1
      a(2,4)   = z2*z2*z2
      a(3,4)   = z3*z3*z3
      a(4,4)   = z4*z4*z4
      a(5,4)   = z5*z5*z5
      a(1,5)   = z1*z1*z1*z1
      a(2,5)   = z2*z2*z2*z2
      a(3,5)   = z3*z3*z3*z3
      a(4,5)   = z4*z4*z4*z4
      a(5,5)   = z5*z5*z5*z5
c
c            set up the pivoting indexes and vector to store
c            maximum term on each row
c
      error = .false.
      do i = 1, neqns
        iwork(i) = i
        workmax = 0.0
        do j = 1, neqns
          workmax = max( workmax, abs(a(i,j)) )
        end do
        dwork(i) = workmax
      end do
c
c           perform triangulation of the coefficient matrix
c           with parital pivoting. terminate if a divide by zero
c           diagonal occurs even after the pivoting.
c
      do k = 1, neqns - 1
        rmax = 0.0
        do i = k, neqns
          ll = iwork(i)
          if ( dwork(ll) .le. tolerance ) then
            error = .true.
            return
          end if
          r  = abs( a(ll,k) ) / dwork(ll)
          if ( r .gt. rmax ) then
            rmax = r
            j    = i
          end if
        end do
c
        lk       = iwork(j)
        iwork(j) = iwork(k)
        iwork(k) = lk
c
        do i = k + 1, neqns
          ll = iwork(i)
          if ( abs( a(lk,k) ) .le. tolerance ) then
            error = .true.
            return
          end if
          xmult = a(ll,k) / a(lk,k)
          do j = k + 1, neqns
            a(ll,j) = a(ll,j) - xmult * a(lk,j)
          end do
          a(ll,k) = xmult
        end do
      end do
c
c           perform forward reduction of the rhs.
c
      do k = 1, neqns - 1
        do i =  k + 1, neqns
          ll = iwork(i)
          b(ll) = b(ll) - a(ll,k) * b(iwork(k))
        end do
      end do
c
c           perform back substitution to recover solution vector
c
      x(neqns) = b(iwork(neqns)) / a(iwork(neqns),neqns)
      do i = neqns - 1, 1, -1
        ll = iwork(i)
        sum = b(ll)
        do j = i + 1, neqns
          sum = sum - a(ll,j) * x(j)
        end do
        x(i) = sum / a(ll,i)
      end do
c
c             assign solution vector components to coefficients
c
      aa = x(5)
      bb = x(4)
      cc = x(3)
      dd = x(2)
      ee = x(1)
c
      return
c
      end
