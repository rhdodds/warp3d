
c
c
c      various supporting routines for domain integrals
c
c
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
c *                       called only if temperatures have been        *
c *                       at some point over the loading history       *
c *                       see didriv for more discussion               *
c *                                                                    *
c *                       *** supports MPI execution ***               *
c *                       written by:    mcw                           *
c *                       last modified: 3/25/2018 rhd                 *
c *                                                                    *
c **********************************************************************
c
      subroutine di_node_props_setup( do_it_local )
c
      use global_data, only : out, myid, nonode, nelblk
      use j_data, only: count_alpha, snode_alpha_ij, seg_snode_e,
     &                  seg_snode_nu, block_seg_curves,
     &                  process_temperatures, front_nodes, domain_origin
c
      implicit none
c
c             parameters
c
      integer :: do_it_local
c
c             locals
c
      integer :: orig_node, do_it
      real, parameter :: zero=0.0
c
c
c          if we are running in MPI:
c            alert workers to enter this routine. process as described
c          otherwise do serial
c
      call wmpi_alert_slaves( 30 )
c
      if( myid .eq. 0 ) then
            do_it   = do_it_local
            orig_node = front_nodes(domain_origin)
      end if
c
      call wmpi_bcast_int( do_it )
      call wmpi_bcast_int( orig_node )
c
c          do_it = 1 -- allocate and fill J thermal data structures
c          do_it = 2 -- deallocate data structures
c
      if( do_it .eq. 1 ) then
         allocate( count_alpha(nonode) )
         allocate( snode_alpha_ij(nonode,6) )
         allocate( seg_snode_e(nonode) )
         allocate( seg_snode_nu(nonode) )
         allocate( block_seg_curves(nelblk) )
         count_alpha    = 0
         snode_alpha_ij = zero
         seg_snode_e    = zero
         seg_snode_nu   = zero
         block_seg_curves = .false.
         call di_node_props( count_alpha, snode_alpha_ij, seg_snode_e,
     &                       seg_snode_nu, block_seg_curves,
     &                       process_temperatures, orig_node )
         deallocate( count_alpha )
         return
      end if
c
c           release module j_data arrays on any rank
c
      if( do_it .eq. 2 ) call didrive_release( out )
c
      return
      end

c
c **********************************************************************
c *                                                                    *
c * di_node_props - build average nodal values of material properties  *
c *                 at each node in the model. all values are          *
c *                 single-precision reals because double-precision    *
c *                 is overkill                                        *
c *                                                                    *
c *                 written by: mcw                                    *
c *                 last modified: 4/17/2018 rhd                       *
c *                                                                    *
c **********************************************************************
c
      subroutine di_node_props ( count_alpha, snode_alpha_ij,
     &                           seg_snode_e, seg_snode_nu,
     &                           block_seg_curves, process_temperatures,
     &                           orig_node )
c
      use global_data, only : nelblk, elblks, myid, out, iprops,
     a                        props, numprocs, nonode  ! old common.main
      use main_data, only : incmap, incid, fgm_node_values,
     a                      elems_to_blocks, temper_nodes,
     b                      temper_nodes_ref, temper_elems,
     c                      inverse_incidences
      use segmental_curves, only: seg_curve_table, seg_curves_type
c
      implicit none
c
c             parameters
c
      integer :: count_alpha(*), orig_node
      real :: snode_alpha_ij(nonode,6), seg_snode_e(*), seg_snode_nu(*)
      logical :: block_seg_curves(*), process_temperatures
c
c             locals
c
      integer :: block, first_elem, rank, bit_flags, curve_set,
     &           first_curve, curve_set_type, snode,
     &           incptr, i, span, elemno, num_enodes, enode,
     &           num_curves_in_set, curve_no, block_num
      integer, allocatable :: block_props(:)
      double precision ::  elem_uniform_temp, enode_temper,
     &                     alpha_temper, e_temper, nu_temper
      double precision, external :: linear_interpolate
      double precision, dimension (:), allocatable :: curve_temp,
     &          curve_es, curve_nus, curve_alphas
      real :: alphax, alphay, alphaz, alphaxy, alphaxz, alphayz,
     &        sum, rc
      real, parameter :: neg_99=-99.0, fgm_tol=1.0, zero=0.0, one=1.0
      logical :: debug, seg_alphas, fgm_alphas, constant_alphas, flag
      real, dimension(:) :: alpha_copy(6)
c
      debug                 = .false.
c
      constant_alphas       = .false.
      fgm_alphas            = .false.
      seg_alphas            = .false.
      process_temperatures  = .false.
      allocate( block_props(nelblk) )
      block_props           = 0
c
c             loop through element blocks.
c             check type of property assignment in each block.
c             property assignment: element nodes (fgm)   = 1
c                                  constant in element   = 2
c                                  temperature dependent = 3
c
      do block = 1, nelblk ! change to threaded at some point
         span       = elblks(0,block)
         first_elem = elblks(1,block)
         rank       = elblks(2,block)
         block_num  = elems_to_blocks(first_elem,1)
c
c             skip block if not handled by current processor
c
         if( rank .ne. myid ) cycle
         if( debug ) write(out,*) "rank ", rank, " block ",
     &                             block, " first_elem ", first_elem
c
c             1. check for fgm alphas. the identifier 'fgm_mark' for
c                fgms assigned in inmat.f is -99.0. elprp.f then assigns
c                this value to props(9,*), props(13,*) and props(34,*).
c
         fgm_alphas = abs( props(9,first_elem) - neg_99 ) .le. fgm_tol
         if( fgm_alphas ) then
            process_temperatures = .true.
            block_props(block)   = 1
            call di_fgm_alphas( span, first_elem, count_alpha,
     &                          snode_alpha_ij )
            cycle
         end if
c
c             2. check for alphas that are constant throughout element.
c                (this is check #2 because fgm assignment of alpha
c                places '-99' in the props array.)
c
c                ordering must match strains, stresses. see inalpha.f
c
         alphax  = props(9,first_elem)
         alphay  = props(13,first_elem)
         alphaz  = props(34,first_elem)
         alphaxy = props(35,first_elem)
         alphayz = props(36,first_elem)
         alphaxz = props(37,first_elem)
         sum     = abs(alphax)  + abs(alphay)  + abs(alphaz)
     &           + abs(alphaxy) + abs(alphaxz) + abs(alphayz)
         if( sum .gt. zero ) then
            constant_alphas      = .true.
            process_temperatures = .true.
            block_props(block)   = 2
            call di_constant_alphas( span, first_elem, count_alpha,
     &                               snode_alpha_ij )
            cycle
         end if
c
c             3. check for temperature-dependent material properties
c                assigned using segmental curves. curve_set_type
c                curve_set_type =0, temperature-independent properties
c                               =1, temperature-dependent properties
c                               =2, strain-rate dependent properties
c
         bit_flags  = iprops(24,first_elem)
         seg_alphas = iand( bit_flags,4 ) > 0
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
      end do !  over all blocks
c
c             average sum of alpha values at each node in structure.
c
c             for MPI:
c             reduce counts to root and re-broadcast all
c             reduce nodal alpha_ij values to root and re-broadcast all
c             for simplicity, each rank then computes nodal averages
c             for all nodes in model
c
      call wmpi_allreduce_int( count_alpha, nonode )
      call wmpi_allreduce_real( snode_alpha_ij, nonode*6 )
c
      do snode = 1,nonode
         if( count_alpha(snode) .ne. 0 ) then
            rc                      = one / real(count_alpha(snode))
            snode_alpha_ij(snode,:) = snode_alpha_ij(snode,:) * rc
         end if
      end do
c
c             send flags to root processor. if the flag is true on any
c             worker processor, set the flag on root to true. send
c             reduced block_seg_curves logical vector to all ranks for
c             later use during actual domain computations
c
      call wmpi_allreduce_vec_log( process_temperatures, 1 )
      call wmpi_allreduce_vec_log( block_seg_curves, nelblk )
c
c             for temperature-dependent material properties, we need to
c             make sure that root has the array data from which to obtain
c             e_front, nu_front, and alpha_front.
c
c             alpha_ij values are reduced to root, bcast and averaged on all
c             ranks above - so root has the full, correct alpha_ij array.
c
c             for temperature dependent stress-strain, each block
c             above puts its e, nu value into the global snode
c             locations - it does not add to what is there. We can't use
c             sum overall ranks, bcast and average since not all ranks
c             likely do not contribute to seg_snode_e, seg_snode_nu.
c
c             we use the MPI reduction operatore MAX which keeps on root
c             the maximum value of each vector location sent by the workers.
c             a neat solution here.
c
c             But we only need to do this for the crack front orig_node,
c             not the whole nonode vector.
c
      if( process_temperatures ) then
         call wmpi_reduce_real_max( seg_snode_e(orig_node), 1)
         call wmpi_reduce_real_max( seg_snode_nu(orig_node), 1 )
      end if
c
      if( debug ) then
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
      end subroutine di_node_props
c
c **********************************************************************
c *                                                                    *
c * di_fgm_alphas - retrieve nodal values of alpha assigned through    *
c *                 'fgm' input. this routine was separated from       *
c *                 di_node_props for clarity.                         *
c *                                                                    *
c *                 written by: mcw                                    *
c *                 last modified: 4/17/2018 rhd                       *
c *                                                                    *
c **********************************************************************
c
      subroutine di_fgm_alphas( span, first_elem, count_alpha,
     &                          snode_alpha_ij )
c
      use global_data, only : iprops, nonode ! old common.main
      use main_data, only : incmap, incid, fgm_node_values
c
      implicit none
c
c             parameters
c
      integer :: span, first_elem, count_alpha(*)
      real :: snode_alpha_ij(nonode,6)
c
c             local variables
c
      integer :: i, elemno, num_enodes, incpos, enode, snode
      real :: x
      real, parameter :: zero=0.0
c
c             assign fgm alpha values to nodes of elements in current
c             block from fgm_alpha_values.
c
c             'fgm' alphas are isotropic. col 3 of fgm values
c
c             loop through elements in block
c
      do i = 1, span
         elemno     = first_elem + ( i - 1 )
         num_enodes = iprops(2,elemno)
         incpos     = incmap(elemno) - 1
c
c             loop through nodes on current element and assign fgm alpha.
c             we sum element contributions to make logic easier for
c             averaging values when some blocks use fgm alphas and some
c             use constant alphas. fgm alphas must be isotropic.
c
         do enode = 1, num_enodes
            snode = incid( incpos + enode )
            count_alpha(snode) = count_alpha(snode) + 1
            x = fgm_node_values(snode,3)
            snode_alpha_ij(snode,1:3) = snode_alpha_ij(snode,1:3) + x
            snode_alpha_ij(snode,4:6) = zero
         end do
      end do
c
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
c *                 last modified: 4/17/2018 rhd                       *
c *                                                                    *
c **********************************************************************
c
      subroutine di_constant_alphas( span, first_elem, count_alpha,
     &                               snode_alpha_ij )
      use global_data, only : iprops, props, nonode ! old common.main
      use main_data, only : incmap, incid
      implicit integer (a-z)
c
c             parameters
c
      integer :: span, first_elem, count_alpha(*)
      real :: snode_alpha_ij(nonode,6)
c
c             local variables
c
      integer :: i, elemno, num_enodes, incpos, enode, snode
      real :: a(6)
c
c             assign alpha values from element constant values defined
c             in the material property definition for the block
c
c             loop through elements in block
c
      do i = 1, span
         elemno     = first_elem + ( i - 1 )
         num_enodes = iprops(2,elemno)
         incpos     = incmap(elemno) - 1
         a(1) = props(9,elemno)   ! alphax
         a(2) = props(13,elemno)  ! alphay
         a(3) = props(34,elemno)  ! alphaz
         a(4) = props(35,elemno)  ! alphaxy
         a(5) = props(36,elemno)  ! alphayz
         a(6) = props(37,elemno)  ! alphaxz
c
c             loop through nodes on current element and assign them the
c             alpha value of the element.
c
         do enode = 1, num_enodes
            snode                   = incid(incpos + enode )
            count_alpha(snode)      = count_alpha(snode) + 1
            snode_alpha_ij(snode,:) = snode_alpha_ij(snode,:) + a
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
c *                 last modified: 4/17/2018 rhd                       *
c *                                                                    *
c **********************************************************************
c
      subroutine di_seg_alpha_e_nu( curve_set, num_curves_in_set,
     &                  first_curve, curve_set_type, span, first_elem,
     &                  count_alpha, snode_alpha_ij, seg_snode_e,
     &                  seg_snode_nu )
c
      use global_data, only : out, iprops, nonode  ! old common.main
      use main_data, only : incmap, incid, temper_elems, temper_nodes,
     &                      temper_nodes_ref
      use segmental_curves, only : seg_curve_table, seg_curves_value,
     &                             seg_curves_ym, seg_curves_nu,
     &                             seg_curves_alpha
      implicit none
c
c             parameters
c
      integer :: curve_set, curve_set_type, num_curves_in_set,
     &           first_curve, span, first_elem, count_alpha(*)
      real :: snode_alpha_ij(nonode,6), seg_snode_e(*), seg_snode_nu(*)
c
c             local variables
c
      integer :: i, curve_no, elemno, num_enodes, incpos, enode, snode
      double precision :: elem_uniform_temp, enode_temper, alpha_temper,
     &                    e_temper, nu_temper, x
      double precision, external :: linear_interpolate
      double precision, dimension (:), allocatable :: curve_temp,
     &                    curve_es, curve_nus, curve_alphas
      real, parameter :: zero=0.0
      logical :: debug
c
      debug = .false.
c
      if( debug ) write(out,1000) curve_set, curve_set_type,
     &            num_curves_in_set, first_curve
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
         if( debug )
     &       write(out,2000) i, curve_no, curve_temp(i), curve_es(i),
     &                   curve_nus(i), curve_alphas(i)
c
      end do
c
c             loop through elements in current block
c
      do i = 1, span
         elemno            = first_elem + ( i - 1 )
         num_enodes        = iprops(2,elemno)
         incpos            = incmap(elemno) - 1
         elem_uniform_temp = temper_elems(elemno)
c
c             loop through nodes on each element. get temperature at
c             each node. do not remove reference temperatures of model
c             nodes (temper_nodes_ref) -- we need absolute temperatures
c             for temperature-dependent values.
c
         do enode = 1, num_enodes
            snode               = incid( incpos + enode )
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
            x                       = alpha_temper
            snode_alpha_ij(snode,1:3) = snode_alpha_ij(snode,1:3) + x
            snode_alpha_ij(snode,4:6) = zero
c
            if( debug )
     &         write(out,3000) elemno, enode, snode,temper_nodes(snode),
     &         elem_uniform_temp, temper_nodes_ref(snode), enode_temper,
     &         alpha_temper, e_temper, nu_temper
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
c *                             di_setup_J7_J8                         *
c *                                                                    *
c * allocate data structures and drive computations                    *
c * to enable subsequent computations of J(7) and J(8) terms.          *
c *                                                                    *
c * if computations not actually required, create small, dummy arrays  *
c * to simplify checking                                               *                                                                           *
c *                                                                    *
c * get a stress work value (W) at all model nodes by extrapoaltion    *
c * from element integration point values, then node averaging         *
c *                                                                    *
c * also get the displacement gradiaents partial(u_i/x_j) (3x3)        *
c * at model nodes by extrapolation from element integration points    *
c *                                                                    *
c *                *** supports MPI execution ***                      *
c *                   written by: mcw                                  *
c *                last modified: 9/15/2018 rhd                        *
c *                                                                    *
c **********************************************************************
c
      subroutine di_setup_J7_J8( do_it_local )
c
      use main_data, only: fgm_node_values, temperatures_ref,
     &                     fgm_node_values_defined,
     &                     fgm_node_values_used,
     &                     initial_state_option,
     &                     initial_stresses_input
      use global_data, only : nonode, out, temperatures, myid, nelblk
      use j_data, only: extrap_counts, swd_at_nodes, strain_at_nodes,
     &                  fgm_e, fgm_nu, comput_i, comput_j,
     &                  displ_grad_at_nodes, j_linear_formulation,
     &                  j_geonl_formulation, block_seg_curves,
     &                  process_initial_state
      implicit none
c
c          parameters
c
      integer :: do_it_local
c
c          local variables
c
      double precision, parameter :: zero=0.0d0
      integer, parameter :: num_flags = 4
      integer :: i, node, flags(num_flags), ios, do_it
      logical :: debug, build_J7_J8_data
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
      if( myid .eq. 0 ) do_it = do_it_local
      call wmpi_bcast_int( do_it )
      call wmpi_bcast_log( temperatures ) ! just in case
      call wmpi_bcast_log( temperatures_ref ) !    "
      call wmpi_bcast_log( initial_state_option ) ! "
      call wmpi_bcast_log( comput_j )
      call wmpi_bcast_log( comput_i )
      call wmpi_bcast_log( j_linear_formulation )
      call wmpi_bcast_log( process_initial_state )
      j_geonl_formulation = .not. j_linear_formulation
c
c          if do_it = 2 -- deallocate all module J-data.
c                          releases any rank
c
      if( do_it .eq. 2 ) then
         call didrive_release( out )
         return
      end if
c
c          do_it = 1 -- allocate and fill node-value data structures
c                       deallocate here is belts + suspenders
c
      flags = 0
      if( allocated( extrap_counts ) )
     &    deallocate( extrap_counts, stat = flags(1) )
      if( allocated( swd_at_nodes ) )
     &    deallocate( swd_at_nodes, stat = flags(2) )
      if( allocated( strain_at_nodes ) )
     &    deallocate( strain_at_nodes, stat = flags(3) )
      if( allocated(displ_grad_at_nodes) )
     &    deallocate( displ_grad_at_nodes )
      if( any(flags .ne. 0) ) then
         write(out, 100) flags
         call die_abort
      end if
c
c           even if node average terms not required,
c           extrap_counts, swd_at_nodes and strain_at_nodes
c           must be allocated for use in dielem_*.f.
c
c           allocate minimal size arrays. routines that need these
c           data will chhck sizes first.
c
      flags = 0
      allocate( extrap_counts(1), stat = flags(1) )
      allocate( swd_at_nodes(1), stat = flags(2) )
      allocate( strain_at_nodes(1,1), stat = flags(3) )
      allocate( displ_grad_at_nodes(1,1), stat = flags(4) )
      if( any(flags .ne. 0) ) then
         write(out, 200) flags
         call die_abort
      end if
      extrap_counts       = 0
      swd_at_nodes        = zero
      strain_at_nodes     = zero
      displ_grad_at_nodes = zero
c
c           determine if 'fgm' values of e or nu are assigned in any
c           element properties.
c
c           interaction integrals need this info
c
      fgm_e  = .false. ;  fgm_nu = .false.
      if( fgm_node_values_used ) then
          fgm_e  = any( fgm_node_values(1:nonode,1) .ne. zero )
          fgm_nu = any( fgm_node_values(1:nonode,2) .ne. zero )
      end if
c
c             determine when it is required to build nodal averaged
c             values of terms needed for J7, J8 or for
c             interaction integrals. block seg_curves -> temperature
c             dependent stress-strain curves appear in model.
c
      build_J7_J8_data = .false.
      if( fgm_node_values_used )            build_J7_J8_data = .true.
      if( comput_i .and. temperatures )     build_J7_J8_data = .true.
      if( comput_i .and. temperatures_ref ) build_J7_J8_data = .true.
      if( process_initial_state )           build_J7_J8_data = .true.
      if( initial_stresses_input )          build_J7_J8_data = .true.
            if( allocated( block_seg_curves ) ) then
        if( any( block_seg_curves ) )       build_J7_J8_data = .true.
      end if
c
      if( .not. build_J7_J8_data ) return
c
c             will compute J7 and J8 terms. allocate model node
c             arrays based on small or large displacement formulation
c             for teh domain. allocates dummies for unsed
c             data for passing as parameters and releasing.
c
      deallocate( extrap_counts, swd_at_nodes, strain_at_nodes,
     &            displ_grad_at_nodes )
      flags = 0
      allocate( extrap_counts(nonode), stat = flags(1) )
      allocate( swd_at_nodes(nonode), stat = flags(2) )
      if( j_linear_formulation ) then
         allocate( strain_at_nodes(6,nonode), stat = flags(3) )
         allocate( displ_grad_at_nodes(1,1), stat = flags(4) )
      end if
      if( j_geonl_formulation ) then
         allocate( strain_at_nodes(1,1), stat = flags(3) )
         allocate( displ_grad_at_nodes(9,nonode), stat = flags(4) )
      end if
      if( any(flags(1:4) .ne. 0) ) then
         write(out, 200) flags
         call die_abort
      end if
      extrap_counts       = 0
!DIR$ VECTOR ALIGNED
      swd_at_nodes        = zero
!DIR$ VECTOR ALIGNED
      strain_at_nodes     = zero
!DIR$ VECTOR ALIGNED
      displ_grad_at_nodes = zero
c
      call di_node_vals( extrap_counts, swd_at_nodes, strain_at_nodes,
     &                   displ_grad_at_nodes )
      deallocate( extrap_counts )
c
      if( debug .and. myid == 0 ) then
         write(out,*) ' ... extrapolation to model nodes ...'
         do i = 1, nonode
c           write(out,9000) i, extrap_counts(i), swd_at_nodes(i)
           write(out,9100) i, strain_at_nodes(1:3,i)
           write(out,9110) displ_grad_at_nodes(1,i),
     &                     displ_grad_at_nodes(5,i),
     &                     displ_grad_at_nodes(9,i)
         end do
      end if
c
      return
c
 100  format('>>>>> deallocation error in di_fgm_setup:',10i6)
 200  format('>>>>> allocation error in di_fgm_setup:', 10i6)
 9000 format(2x,i7,i7,d14.6)
 9100 format(2x,i7,6d14.6)
 9110 format(2x,7x,6d14.6)
c
      end
c
c **********************************************************************
c *                                                                    *
c * di_node_vals - build average nodal values of the stress work       *
c *          denisty (W) and (1) strains or (2)  displacement          *
c *          gradients at model  nodes by extrapolation from element   *
c *          values and averaging                                      *
c *                                                                    *
c *                   written by: rhd                                  *
c *                last modified: 4/17/2018 rhd                        *
c *                                                                    *
c **********************************************************************
c
      subroutine di_node_vals( extrap_counts, swd_at_nodes,
     &                         strain_at_nodes, displ_grad_at_nodes  )
      use global_data, only : noelem, nonode, myid, out, elblks
      use main_data, only : elems_to_blocks
      use j_data, only : j_linear_formulation, j_geonl_formulation
      implicit none
c
c         parameters
c
      integer :: extrap_counts(nonode)
      double precision :: swd_at_nodes(nonode),
     &                    strain_at_nodes(6,nonode),
     &                    displ_grad_at_nodes(9,nonode)
c
c         locals (mxgp, mxndel are parameters in global_data
c
      integer :: i, elemno, blk, enode, snode, skipped_killed
      double precision :: rc, time_start, time_end
      double precision, external :: omp_get_wtime
      double precision, parameter :: one=1.0d0
      logical, parameter :: debug = .false.
      logical, external :: chk_killed
c
c                build average nodal values of the stress work density
c                (swd) and displacement derivatives in the model.
c                for each element,
c                shape functions will be used to extrapolate integration-
c                point values to each of the element's nodes where they
c                are then averaged.
c
c                loop over all structure elements. handle
c                only elements whose data is owned by this processor.
c
      time_start = omp_get_wtime()
c$OMP PARALLEL DO PRIVATE( elemno, blk )
      do elemno = 1, noelem
         blk = elems_to_blocks(elemno,1)
         if( elblks(2,blk) .ne. myid ) cycle
         if( chk_killed( elemno ) ) then ! function call
c$OMP ATOMIC
            skipped_killed = skipped_killed + 1 ! init in dicmj
            cycle
         end if
         call di_nod_vals_one_elem( elemno, blk, extrap_counts,
     &                              swd_at_nodes, strain_at_nodes,
     &                              displ_grad_at_nodes  )
      end do  ! over structure elements
c$ OMP END PARALLEL DO
      time_end = omp_get_wtime()
c      write(out,*) '... eleme loop time: ', time_end - time_start
c
c                reduce vectors to root then broadcast to all.
c                compute nodal average values. each rank thus has a
c                complete copy for all model nodes.
c
      call wmpi_allreduce_int( extrap_counts, nonode )
      call wmpi_allreduce_dble( swd_at_nodes, nonode )
c
      do snode = 1, nonode
       if( extrap_counts(snode) .eq. 0 ) cycle
       rc = one / dble( extrap_counts(snode) )
       swd_at_nodes(snode) = swd_at_nodes(snode) * rc
      end do

c
      if( j_linear_formulation ) then
        call wmpi_allreduce_dble( strain_at_nodes, 6*nonode )
        do snode = 1, nonode
          if( extrap_counts(snode) .eq. 0 ) cycle
          rc = one / dble( extrap_counts(snode) )
!DIR$ VECTOR ALIGNED
          strain_at_nodes(1:6,snode) = strain_at_nodes(1:6,snode) * rc
        end do
      end if

c
      if( j_geonl_formulation ) then
       call wmpi_allreduce_dble( displ_grad_at_nodes, 9*nonode )
        do snode = 1, nonode
          if( extrap_counts(snode) .eq. 0 ) cycle
          rc = one / dble( extrap_counts(snode) )
!DIR$ VECTOR ALIGNED
          displ_grad_at_nodes(1:9,snode) =
     &        displ_grad_at_nodes(1:9,snode) * rc
        end do
      end if
c
      return
c
      end
c **********************************************************************
c *                                                                    *
c * di_nod_vals_one_elem - process 1 element to create element nodal   *
c *         values &  add to model model nodes to average computations *
c *                                                                    *
c *                   written by: rhd                                  *
c *                last modified: 6/15/2018 rhd                        *
c *                                                                    *
c **********************************************************************
c
      subroutine di_nod_vals_one_elem( elemno, blk, extrap_counts,
     &      swd_at_nodes, strain_at_nodes, displ_grad_at_nodes  )
c
      use global_data, only : out, nstrs, nstr, iprops, lprops, nonode,
     &                        mxndel, mxgp, mxndel, mxvl, mxoupr,
     &                        scoords => c, sdispl => u, dstmap
      use main_data, only : incmap, crdmap, incid, elems_to_blocks,
     &                      trn, trnmat
      use j_data, only : j_linear_formulation, j_geonl_formulation,
     &                   process_initial_state
      use elem_block_data, only : urcs_n_blocks, eps_n_blocks,
     &                            initial_state_data
c
      implicit none
c
c         parameters
c
      integer :: elemno, blk, extrap_counts(nonode)
      double precision :: swd_at_nodes(nonode),
     &                    strain_at_nodes(6,nonode),
     &                    displ_grad_at_nodes(9,nonode)
c
c         locals (mxgp, mxndel are parameters in global_data
c
      integer :: i, j, k, ptno, etype, num_enodes, int_order,
     &           num_gpts, incpos, rel_elem, gpn, sig_offset,
     &           eps_offset, enode, snode, gpt, ierr, cpos,
     &           e_snodes(mxndel), last
      integer, parameter :: nrow_dsf = mxndel
      integer, save :: message_count2=0
      double precision :: swd_at_gpts(mxgp), strain_at_gpts(6,mxgp),
     &     displ_grad_at_gpts(9,mxgp), e_coords(3,mxndel),
     &     e_displ(3,mxndel), dsf(nrow_dsf,3), jacob(3,3), jacobi(3,3),
     &     enode_swd(mxndel), enode_strains(6,mxndel), rnum_gpts,
     &     rcount, enode_displ_grad(9,mxndel), xsi, eta,
     &     zeta, weight, lg(mxgp), R(3,3), t(3), W0, eps_0(6),
     &     e_coords_n(3,mxndel), detvol_0, detvol_n, vol_0, vol_n,
     &     f_vec(9), f_tens(3,3), detf, factor, j_bar
      equivalence (f_vec, f_tens)
      double precision, parameter :: zero=0.0d0, one=1.0d0, two=2.d0
      double precision, pointer :: urcs_n(:), eps_n(:), z(:,:)
      logical :: geonl, linearform, ok
      logical, parameter :: debug = .false.
c
c                get element properties
c
      etype            = iprops(1,elemno)
      num_enodes       = iprops(2,elemno)
      int_order        = iprops(5,elemno)
      num_gpts         = iprops(6,elemno)
      geonl            = lprops(18,elemno)
      rnum_gpts        = dble( num_gpts )
      rel_elem         = elems_to_blocks(elemno,2)
      ierr             = 0
      linearform       = .not. geonl
c
      incpos = incmap(elemno); last = incpos+num_enodes-1
      e_snodes(1:num_enodes) = incid(incpos:last)
c
c                formulation must match that for domain
c
      ok = .true.
      if( j_linear_formulation ) then  ! element must match domain
         if( geonl ) ok = .false.
      end if
      if( j_geonl_formulation ) then
         if( linearform ) ok = .false.
      end if
      if( .not. ok ) then
         message_count2 = message_count2 + 1
         if( message_count2 <=10 ) then
           write(out,9230) elemno
           if( message_count2 == 10 ) write(out,9240)
         end if
         return
      end if
      urcs_n => urcs_n_blocks(blk)%ptr
      eps_n  => eps_n_blocks(blk)%ptr
c
      sig_offset = (rel_elem - 1) * nstrs * num_gpts
      eps_offset = (rel_elem - 1) * nstr  * num_gpts
c
c             1. gather integration point work densities.
c                adjust for plastic value at the user-defined
c                initial state if required
c
!DIR$ VECTOR ALIGNED
      do gpn = 1, num_gpts
         swd_at_gpts(gpn) = urcs_n(sig_offset + 7)
         sig_offset       = sig_offset + nstrs
      end do
c
      if( process_initial_state ) then
        associate( z => initial_state_data(blk)%W_plastic_nis_block )
!DIR$ VECTOR ALIGNED
        do gpn = 1, num_gpts
           W0 = z(rel_elem,gpn)
           swd_at_gpts(gpn) = swd_at_gpts(gpn)  - W0
        end do
        end associate
      end if
c
c             2. for small eps, gather strains at points.
c                strains are in model global coordinates.
c
      if( j_linear_formulation ) then
        do gpn = 1, num_gpts
         i = eps_offset + 1; j = i + 5
         strain_at_gpts(1:6,gpn) = eps_n(i:j)
         eps_offset = eps_offset + nstr
        end do
      end if
c
c             3. for large eps, compute displacement gradient tensor
c                at each integration point. grad (3x3) = F - I.
c                for 8-node hex, F is repalced by F-bar.
c
c                convert local node system {u} -> global as needed
c                global = trans(R) local. displacement grads
c                are in model global coordinates
c
      if( j_geonl_formulation ) then
       do enode = 1, num_enodes
         snode = e_snodes(enode)
         cpos  = crdmap(snode)
         e_coords(1,enode) = scoords(cpos+0)
         e_coords(2,enode) = scoords(cpos+1)
         e_coords(3,enode) = scoords(cpos+2)
         e_displ(1,enode)  = sdispl(dstmap(snode)+0)
         e_displ(2,enode)  = sdispl(dstmap(snode)+1)
         e_displ(3,enode)  = sdispl(dstmap(snode)+2)
         if( trn(snode) ) then
           R = trnmat(snode)%mat  ! 3x3
           t = e_displ(1:3,enode)
           e_displ(1,enode) = R(1,1)*t(1) + R(2,1)*t(2) + R(3,1)*t(3)
           e_displ(2,enode) = R(1,2)*t(1) + R(2,2)*t(2) + R(3,2)*t(3)
           e_displ(3,enode) = R(1,3)*t(1) + R(2,3)*t(2) + R(3,3)*t(3)
         end if
       end do ! over enode
c
       call di_displ_grad_one_elem( elemno, etype, num_enodes,
     &                              num_gpts, int_order, e_coords,
     &                              e_displ, displ_grad_at_gpts )
c
      end if ! on nlgeom
c
      if( debug  .and. (elemno == 77 .or. elemno == 18000) ) then
         write(out,100) elemno
         do gpt = 1, num_gpts
            write(out,210) gpt, strain_at_gpts(1:6,gpt)
            write(out,210) gpt, displ_grad_at_gpts(1:9,gpt)
         end do
      end if
c
c             4. extrapolate integration point values to element nodes
c
!DIR$ VECTOR ALIGNED
      enode_swd        = zero  !  all terms    always
!DIR$ VECTOR ALIGNED
      enode_strains    = zero  !     "         small eps only
!DIR$ VECTOR ALIGNED
      enode_displ_grad = zero  !     "         large eps only
c
      call di_extrap_to_nodes( elemno, etype, num_enodes, int_order,
     &     num_gpts, swd_at_gpts, enode_swd, 1, lg, out )
c
      if( j_linear_formulation )
     &  call di_extrap_to_nodes( elemno, etype, num_enodes, int_order,
     &     num_gpts, strain_at_gpts, enode_strains, 6, lg, out )
c
      if( j_geonl_formulation )
     &  call di_extrap_to_nodes( elemno, etype, num_enodes, int_order,
     &     num_gpts, displ_grad_at_gpts, enode_displ_grad, 9, lg, out )
c
      if( debug .and. (elemno == 77 .or. elemno == 18000) ) then
         write(out,110) elemno
         do enode = 1, num_enodes
            write(out,210) enode, enode_strains(1:6,enode)
            write(out,210) enode, enode_displ_grad(1:9,enode)
         end do
      end if
c
c             add element node values to global structure node-value
c             arrays.
c
      do enode = 1, num_enodes
         snode                = e_snodes(enode)
c$OMP ATOMIC
         extrap_counts(snode) = extrap_counts(snode) + 1
c$OMP ATOMIC
         swd_at_nodes(snode)  = swd_at_nodes(snode) + enode_swd(enode)
      end do
c
      if( j_linear_formulation ) then
        do enode = 1, num_enodes
          snode = e_snodes(enode)
!DIR$ VECTOR ALIGNED
          do k = 1, 6
c$OMP ATOMIC
           strain_at_nodes(k,snode) = strain_at_nodes(k,snode) +
     &                                enode_strains(k,enode)
          end do
        end do
      end if
c
      if( j_geonl_formulation ) then
      do enode = 1, num_enodes
         snode = e_snodes(enode)
!DIR$ VECTOR ALIGNED
         do k = 1, 9
c$OMP ATOMIC
          displ_grad_at_nodes(k,snode) =
     &       displ_grad_at_nodes(k,snode) + enode_displ_grad(k,enode)
         end do
        end do
      end if
c
      return
c
 100  format(/,'*** gauss pt values for elem ',i8,/,'gpt',5x,
     &       15x,'values')
 110  format(/,'*** element node values. elem ',i8,/,15x,'values')
 210  format(i3,1x,9d14.6)
 220  format(i3,1x,9d14.6)
 230  format(4x,9d14.6)
 9230 format(/,'>> Warning: all elements in domain must have the same',
     &  ' formulation (small strain, larger strain. element: ',i8,
     & /,    '           element skipped.')
 9240 format(/,'>> Warning: *** no more messages about formulation',
     &   ' mismatch ***')
c
      end subroutine di_nod_vals_one_elem
c
c ****************************************************************
c *                                                              *
c *     subroutine di_displ_grad_one_elem                        *
c *                                                              *
c *     compute displacement gradient tensor at each integration *
c *     point of one element. for 8-node element, these are      *
c *     taken from F-bar rather than F                           *
c *                                                              *
c *                   written by: rhd                            *
c *                last modified: 04/22/2018 rhd                 *
c *                                                              *
c ****************************************************************
c
      subroutine di_displ_grad_one_elem( elemno, etype, num_enodes,
     &        num_gpts, int_order, e_coords, e_displ, displ_grads )
c
      use global_data, only : out, mxndel, mxgp
c
      implicit none
c
c              parameters
c
      integer :: elemno, etype, num_enodes, num_gpts, int_order
      double precision :: e_coords(3,mxndel), e_displ(3,mxndel),
     &                    displ_grads(9,mxgp)
c
c              locals
c
      integer :: ptno, ierr, enode, i
c
      double precision :: jacob(3,3), jacobi(3,3), F_vec(9,mxgp),
     &                    F_tensor(3,3,mxgp), dsf_1(mxndel),
     &                    dsf_2(mxndel), dsf_3(mxndel)
      equivalence( F_vec, F_tensor )
c
      double precision :: xsi, eta, zeta, weight, detvol_0,
     &                    jbar, vol_0, dux, dvx, dwx,
     &                    duy, dvy, dwy, duz, dvz, dwz, nx, ny, nz,
     &                    J_bar, detf(mxgp), factor, sum_detF_dvel
      double precision, parameter :: third = 1.0d0/3.0d0, zero = 0.0d0,
     &                               one = 1.0d0
c
c              compute the displacement gradient tensor F-I at
c              each integration point.
c
c              for 8-node hex, we compute F-bar and then gradient
c              tensor = (F-bar) - I
c
c              1. get F at each integration point
c                  1a. get integration point parametric coords
c                  1b. parametric derivatives at point for each enode
c                  1c. coordinate Jacobian, inverse, det for undeformed
c                      element shape
c                  1d. add to element undeformed volume
c                  1e. compute each enode contribution to displacement
c                      gradents at this integration point. derivatives
c                      are wrt undeformed shape
c                  1f. add 1.0 to diagonals to make F
c
c                  note equivalencing of 3x3 tensors and 9x1 vectors
c                  for convenience
c
      sum_detF_dvel = zero
      vol_0         = zero
      associate( F => F_tensor )
c
      do ptno = 1, num_gpts
c
        call getgpts( etype, int_order, ptno, xsi, eta, zeta, weight )
        call derivs( etype, xsi, eta, zeta, dsf_1, dsf_2, dsf_3 )
        call dielcj( dsf_1, dsf_2, dsf_3, e_coords, num_enodes,
     &               jacob, jacobi, detvol_0, ierr, out, .false. )
        if( ierr .ne. 0 ) then
          write(out,9100) elemno, ptno; call die_abort
        end if
c
        vol_0 = vol_0 + detvol_0  ! weight = 1 for 2x2x2
c
        dux   = zero
        dvx   = zero
        dwx   = zero
        duy   = zero
        dvy   = zero
        dwy   = zero
        duz   = zero
        dvz   = zero
        dwz   = zero
!DIR$ VECTOR ALIGNED
        do enode = 1, num_enodes
          nx = dsf_1(enode) * jacobi(1,1) + dsf_2(enode) * jacobi(1,2) +
     &         dsf_3(enode) * jacobi(1,3)
          ny = dsf_1(enode) * jacobi(2,1) + dsf_2(enode) * jacobi(2,2) +
     &         dsf_3(enode) * jacobi(2,3)
          nz = dsf_1(enode) * jacobi(3,1) + dsf_2(enode) * jacobi(3,2) +
     &         dsf_3(enode) * jacobi(3,3)
          dux = dux + nx * e_displ(1,enode)
          dvx = dvx + nx * e_displ(2,enode)
          dwx = dwx + nx * e_displ(3,enode)
          duy = duy + ny * e_displ(1,enode)
          dvy = dvy + ny * e_displ(2,enode)
          dwy = dwy + ny * e_displ(3,enode)
          duz = duz + nz * e_displ(1,enode)
          dvz = dvz + nz * e_displ(2,enode)
          dwz = dwz + nz * e_displ(3,enode)
        end do ! over enodes
c
        F_tensor(1,1,ptno) = dux + one
        F_tensor(2,1,ptno) = dvx
        F_tensor(3,1,ptno) = dwx
        F_tensor(1,2,ptno) = duy
        F_tensor(2,2,ptno) = dvy + one
        F_tensor(3,2,ptno) = dwy
        F_tensor(1,3,ptno) = duz
        F_tensor(2,3,ptno) = dvz
        F_tensor(3,3,ptno) = dwz + one
        detf(ptno) =  F(1,1,ptno) * F(2,2,ptno) * F(3,3,ptno) +
     &                F(2,1,ptno) * F(3,2,ptno) * F(1,3,ptno)
     &              + F(3,1,ptno) * F(1,2,ptno) * F(2,3,ptno) -
     &                F(1,1,ptno) * F(3,2,ptno) * F(2,3,ptno)
     &              - F(2,1,ptno) * F(1,2,ptno) * F(3,3,ptno) -
     &                F(3,1,ptno) * F(2,2,ptno) * F(1,3,ptno)
        sum_detF_dvel = sum_detF_dvel + detF(ptno) * detvol_0
c
      end do ! over points
c
c              3. Unless element is 8-node hex, we'ere done.
c                 Set displacement gradient tensor as F - I at
c                 each integration point and leave.
c
      if( etype .ne. 2 ) then ! 8-node hex is type 2
!DIR$ VECTOR ALIGNED
      do ptno = 1, num_gpts
!DIR$ VECTOR ALIGNED
      displ_grads(1:9,ptno) = F_vec(1:9,ptno)
           displ_grads(1,ptno)   = F_vec(1,ptno) - one
           displ_grads(5,ptno)   = F_vec(5,ptno) - one
           displ_grads(9,ptno)   = F_vec(9,ptno) - one
        end do
        return
      end if
c
c              4. for 8-node hex, replace F with F-bar
c
      J_bar = sum_detF_dvel / vol_0
c
!DIR$ VECTOR ALIGNED
      do ptno = 1, num_gpts
        factor = (J_bar/ detf(ptno) ) ** third
        displ_grads(1,ptno) = (F_vec(1,ptno) * factor) - one
        displ_grads(2,ptno) =  F_vec(2,ptno) * factor
        displ_grads(3,ptno) =  F_vec(3,ptno) * factor
        displ_grads(4,ptno) =  F_vec(4,ptno) * factor
        displ_grads(5,ptno) = (F_vec(5,ptno) * factor) - one
        displ_grads(6,ptno) =  F_vec(6,ptno) * factor
        displ_grads(7,ptno) =  F_vec(7,ptno) * factor
        displ_grads(8,ptno) =  F_vec(8,ptno) * factor
        displ_grads(9,ptno) = (F_vec(9,ptno) * factor) - one
      end do
c
      end associate
c
      return
c
 9100 format('>>>>> Fatal error: in di_displ_grad_one_elem.',
     &  /,   '                   bad coordinate [J] ele, pt: ',2i8,
     &     /,'                   job terminated.' )
 9110 format('>>>>> Fatal error: in di_displ_grad_one_elem.',
     &  /,   '                   bad coordinate [J] ele,pt: ',2i8,
     &     /,'                   deformed shape. job terminated.' )
c
      end
c
c ****************************************************************
c *                                                              *
c *     extrapolate integration point values to nodes for 1 elem *
c *                                                              *
c *                   written by: mcw                            *
c *                last modified: 04/15/2018 rhd                 *
c *                                                              *
c ****************************************************************
c
      subroutine di_extrap_to_nodes( elemno, etype, num_enodes,
     &     int_order, num_gpts, gpt_vals, enode_vals, nvalues,
     &     lg, out )
c
      implicit none
c
c              parameters
c
      integer :: elemno, etype, num_enodes, int_order, num_gpts, out,
     &           nvalues
      double precision :: gpt_vals(nvalues,*), enode_vals(nvalues,*),
     &                    lg(*)
c
c              locals
c
      integer :: gpn, enode, idummy_vec(1)
      double precision :: rdummy_vec(1), xi, eta, zeta
      double precision, parameter ::  one=1.0d0, zero=0.d0
      logical :: lagrangian_extrap, debug
c
      debug = .false.
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
      if( etype .eq. 2 ) lagrangian_extrap = .true.
      if( int_order .eq. 1 .or. int_order .eq. 8 )
     &         lagrangian_extrap = .true.
c
      if ( lagrangian_extrap ) then
          call  di_extrap_to_nodes_lag
      else  ! node values are all = integration pt averages
          call di_extrap_to_nodes_not_lag
      end if
c
      return
c
      contains
c     ========
c
      subroutine di_extrap_to_nodes_not_lag
      implicit none
      double precision :: x
c
c              to get nodal values for the element, extrapolation
c              from integration points to node points is
c              not possible. instead we use the average of all
c              integration points at each node. node arrays zeroed
c              before call
c
      x = one / dble( num_gpts )
c
      do gpn = 1, num_gpts
!DIR$ VECTOR ALIGNED
      enode_vals(1:nvalues,1) = enode_vals(1:nvalues,1) +
     &                              gpt_vals(1:nvalues,gpn)
      end do
c
      enode_vals(1:nvalues,1) = enode_vals(1:nvalues,1) * x
c
      do enode = 2, num_enodes
!DIR$ VECTOR ALIGNED
      enode_vals(1:nvalues,enode) = enode_vals(1:nvalues,1)
      end do
c
      return
c
      end subroutine di_extrap_to_nodes_not_lag
c
      subroutine di_extrap_to_nodes_lag
      implicit none
c
c              extrapolate integration point values to nodes using
c              lagrangian shape functions.
c                a) find the lagrangian shape functions for
c                   each gauss point at the current element node.
c                b) extrapolate values from integration
c                   points to an element node.
c
      do enode = 1, num_enodes
c
       call ndpts1( idummy_vec, 0, rdummy_vec, etype, enode,
     &              xi, eta, zeta )
       call oulgf( etype, xi, eta, zeta, lg, int_order )
c
       do gpn = 1, num_gpts
!DIR$ VECTOR ALIGNED
       enode_vals(1:nvalues,enode) =  enode_vals(1:nvalues,enode) +
     &                        gpt_vals(1:nvalues,gpn) * lg(gpn)
         end do
c
      end do  ! over enodes
c
      return
c
      end subroutine di_extrap_to_nodes_lag
c
      end subroutine di_extrap_to_nodes
c
c****************************************************************
c                                                               *
c      subroutine to write j and i-integral data to             *
c      standard output                                          *
c                    written by: mcw                            *
c                 last modified: 3/6/2018 rhd                   *
c                                                               *
c****************************************************************
c
      subroutine di_write_std_out( ltmstp, stname, lsldnm, out )
c
      use j_data, only : ring_count, comput_j, comput_i,
     a    output_packet_j, domain_id, output_packet_i,
     b    ks_from_j, j_storage, domain_min_j, j_from_ks,
     c    domain_max_j, i_storage, domain_avg_i, domain_avg_j,
     d    domain_max_i, domain_min_i, static_j, static_avg, static_max,
     e    face_loading, symmetric_domain, out_pstrain, out_pstress,
     f    rings_given, j_to_k, static_min, process_temperatures
c
      implicit none
c
c             parameters
c
      integer :: ltmstp, out
      character(len=8) :: stname, lsldnm
c
c             local variables
c
      integer :: i, j, ring, skipped_killed
      double precision :: rg_count
c
      rg_count = dble(ring_count)
c
      write( out,9000 )
c
c             write j-integral results to standard output
c
      if( comput_j ) then
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
 9000 format(/)
 9005 format(/,3x,'*****   Total J-values   *****')
 9010 format(/,1x,45x,'J-integral components',
     &  /,      1x,42x,27('='),
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
     &       /,1x,'* dm4: acceleration & velocity gradient term',
     &       /,1x,'* dm5: crack face loading term',
     &       /,1x,'* dm6: thermal loading term (zero for an fgm)',
     &       /,1x,'* dm7: stress times displacement gradient term',
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
     &       //,1x,'domain',8x,'II1',10x,'II2',10x,'II3',10x,'II4',
     &       10x,'II5',10x,'II6',10x,'II7',10x,'II8',7x,'total I',5x,
     &       'KI (pstrs)',1x,'killed ele' )
 9060 format(/,1x,i5,2x,10(2x,e11.4),3x,'(',i3,')')
 9065 format(/,1x,i5,2x,11(2x,e11.4),3x,'(',i3,')')
 9075 format(/,1x,30x,'I-integral components',
     &       /,26x,'-- KI = 1, KII = 0, KIII = 0, plane strain --',
     &       //,1x,'domain',8x,'II1',10x,'II2',10x,'II3',10x,'II4',
     &       10x,'II5',10x,'II6',10x,'II7',10x,'II8',7x,'total I',5x,
     &       'KI (pstrn)',1x,'killed ele' )
 9085 format(/,1x,30x,'I-integral components',
     &       /,26x,'-- KI = 0, KII = 1, KIII = 0, plane stress --',
     &       //,1x,'domain',8x,'II1',10x,'II2',10x,'II3',10x,'II4',
     &       10x,'II5',10x,'II6',10x,'II7',10x,'II8',7x,'total I',5x,
     &       'KII (pstrs)',1x,'killed ele' )
 9095 format(/,1x,30x,'I-integral components',
     &       /,26x,'-- KI = 0, KII = 1, KIII = 0, plane strain --',
     &       //,1x,'domain',8x,'II1',10x,'II2',10x,'II3',10x,'II4',
     &       10x,'II5',10x,'II6',10x,'II7',10x,'II8',7x,'total I',5x,
     &       'KII (pstrn)',1x,'killed ele' )
 9105 format(/,1x,30x,'I-integral components',
     &       /,26x,'-- KI = 0, KII = 0, KIII = 1 --',
     &       //,1x,'domain',8x,'II1',10x,'II2',10x,'II3',10x,'II4',
     &       10x,'II5',10x,'II6',10x,'II7',10x,'II8',7x,'total I',8x,
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
     &       //,1x,'domain',8x,'II1',10x,'II2',10x,'II3',10x,'II4',
     &       10x,'II5',10x,'II6',10x,'II7',10x,'II8',7x,'total I',4x,
     &       'T11 (pstrs)',6x,'T33',5x,'killed ele' )
 9110 format(/,1x,30x,'I-integral components',
     &       /,26x,'--T-stress ("plane strain" auxiliary fields)--',
     &       //,1x,'domain',8x,'II1',10x,'II2',10x,'II3',10x,'II4',
     &       10x,'II5',10x,'II6',10x,'II7',10x,'II8',7x,'total I',4x,
     &       'T11 (pstrn)',6x,'T33',5x,'killed ele' )
 9112 format(/,1x,30x,'I-integral components',
     &       /,26x,'--T-stress (anti-plane shear auxiliary fields)--',
     &       //,1x,'domain',8x,'II1',10x,'II2',10x,'II3',10x,'II4',
     &       10x,'II5',10x,'II6',10x,'II7',10x,'II8',7x,'total I',8x,
     &       'T13',5x,'killed ele' )
 9115 format(/,3x,' average  ',3x,' minimum  ',3x,' maximum',
     &       /,1x, 3(e11.4,2x),
     &       /,1x, 3(e11.4,2x))
 9125 format(/,1x,'* II1: stress * deriv of aux. displ.',
     &       /,1x,'* II2: aux stress * derivative of displacement',
     &       /,1x,'* II3: mixed strain energy density',
     &       ' (aux stress * strain)',
     &       /,1x,'* II4: stress * 2nd deriv of aux displ',
     &       /,1x,'* II5: stress * deriv of aux strain)',
     &       /,1x,'* II6: deriv of constitutive tensor * strain * aux',
     &       'strain',
     &       /,1x,'* II7: thermal strains (not yet implemented)',
     &       /,1x,'* II8: traction * deriv of aux displ',
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
c                 last modified: 4/6/2018 rhd                   *
c                                                               *
c****************************************************************
c
      subroutine di_write_packets( ltmstp, stname, lsldnm )
c
      use main_data, only : output_packets, packet_file_no
      use j_data, only : ring_count, comput_j, comput_i,
     a    output_packet_j, domain_id, output_packet_i,
     b    ks_from_j, j_storage, domain_min_j, j_from_ks,
     c    domain_max_j, i_storage, domain_avg_i, domain_avg_j,
     d    domain_max_i, domain_min_i
c
      implicit none
c
c             parameters
c
      integer :: ltmstp
      character(len=8) :: stname, lsldnm
c
c             local variables
c
      integer :: i, j, ring, skipped_killed, num_lines
      double precision :: rg_count
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
c *                                 last modified: 4/6/2019 rhd     *
c *                                                                 *
c *******************************************************************
c
c
      subroutine di_calc_e33( out )
      use j_data, only : domain_type, num_front_nodes, front_coords,
     a                   front_order, e33_front, comput_i, out_pstrain,
     b                   front_node_displ
      implicit none
c
c             parameter declarations
c
      integer :: out
c
c             local declarations
c
      integer :: i, j, nfn, ngp, nfelem, inc, elem, gp, node
      double precision :: sf(4), dsf(4), jacob, jacobi, xsi,
     &                    displ_coords(3,30), length_undisp,
     &                    length_disp, j1, j2, j3, detj, dieldp, w,
     &                    dummy, x1, z1, x2, z2, x3, z3
      double precision, parameter ::  zero=0.d0, half=0.5d0, one=1.0d0,
     &                                two=2.0d0
      logical :: debug
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
            call di1dsfa( xsi, dsf, sf, nfn )
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
 9900 format(//,' crack front tangential strain (e_33) for ',
     & 'T11, T33: ',e13.6)
c
      end
c
c *******************************************************************
c *                                                                 *
c *   calculate coefficients of curve described by crack front      *
c *   nodes.                                                        *
c *                                                                 *
c *                                    written by: mcw              *
c *                                 last modified: 04/6/2018 rhd    *
c *                                                                 *
c *******************************************************************
c
c
      subroutine di_calc_curvature( debug, out )
      use j_data, only : front_order, num_front_nodes, front_coords,
     a                   crack_curvature
      implicit none
c
c             parameters
c
      integer :: out
      logical :: debug
c
c             local variables
c
      integer :: i, j
      double precision :: x1, z1, x2, z2, x3, z3, x4, z4,
     & x5, z5, a, b, c, d, e, x_center, z_center, circle_radius, q
      double precision, parameter :: zero=0.0d0, half=0.5d0, one=1.0d0,
     &                               two=2.0d0, toler=1.0d-04
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
      if( debug ) write(out,*) "@5"

      if( nint(crack_curvature(1)) .eq. 0 ) go to 1111
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
     &                                 x5, z5, aa, bb, cc, dd, ee,
     &                                 out, debug )
      implicit none
c
c             parameters
c
      integer :: out
      double precision :: x1, z1, x2, z2, x3, z3, x4, z4, x5, z5, aa,
     &                    bb, cc, dd, ee
      logical :: debug
c
c             locals
c
      logical :: error
      integer :: nrow_a, neqns, i, j, k, ll, lk
      integer :: iwork(5)
      double precision :: workmax, r, rmax, xmult, sum
      double precision :: a(5,5), x(5), b(5), dwork(5)
      double precision, parameter :: zero=0.0d0, half=0.5d0, one=1.0d0,
     &                               two=2.0d0, tolerance=1.0d-20
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
        workmax = zero
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
c *******************************************************************
c *                                                                 *
c *   subroutine di1dsfa ---- calculates shape function values and  *
c *                           derivatives for 1-dimension           *
c *                                                                 *
c *******************************************************************
c
      subroutine di1dsfa( xsi, dsf, sf, nlnode )
c
c              parameter declarations
c
      implicit none
c
      integer :: nlnode
      double precision  :: xsi, dsf(*), sf(*)
c
c              local declarations
c
      double precision :: xsisqr
      double precision, parameter :: half = 0.5d0, one = 1.0d0,
     &                               two = 2.0d0
c
      select case( nlnode )
c
        case(1, 2)  ! linear
        sf(1)  = half * ( one - xsi )
        sf(2)  = half * ( one + xsi )
        dsf(1) = -half
        dsf(2) =  half
c
        case( 3 )  ! quadratic
        xsisqr = xsi * xsi
        sf(1)  = half * ( xsisqr - xsi )
        sf(2)  = one - xsisqr
        sf(3)  = half * ( xsisqr + xsi )
        dsf(1) = xsi - half
        dsf(2) = -two*xsi
        dsf(3) = xsi + half
c
      end select
c
      return
      end
