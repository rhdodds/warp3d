c ***************************************************************
c *                                                             *
c * domain_compute - drive execution of element routines to     *
c *                  compute j and i-integrals for a single     *
c *                  domain                                     *
c *                                                             *
c *                  -> drive MPI to process all elements       *
c *                                                             *
c *                  last update:  5/31/2018 rhd                *
c *                                                             *
c ***************************************************************
c
c
      subroutine dicmj
c
      use global_data, only : myid, iout=>out, nelblk, elblks,
     a                        worker_processor, numprocs, num_threads
c
      use j_data, only : debug_driver, ring_count, symmetric_domain,
     a  print_elem_values, q_element_maps, q_values,
     b  front_element_list, comput_i, e_front, nu_front, comput_j,
     d  cf_traction_flags, expanded_front_nodes, front_q_area,
     e  face_loading, static_j,
     f  domain_min_j, domain_max_j,
     g  domain_min_i, domain_max_i, domain_avg_i, static_min,
     h  static_max, static_avg, ks_from_j, nowring, domain_avg_j,
     i  e33_front, domain_id, j_from_ks, print_totals, j_storage,
     j  i_storage, size_j_values, size_i_values,
     k  j_geonl_formulation, j_linear_formulation,
     l  temperatures_on_model
c
      implicit none
c
c          local declarations.
c
      integer :: now_thread
      integer, external ::  omp_get_thread_num
c
      double precision ::
     &   jtotal, symm_factor, di_value_j, static_di,
     &   static_total, Ki, T11, T33, T13, sign
c
      double precision :: diterms(size_i_values), itotal(size_i_values),
     &   iiterms(size_i_values,size_i_values),
     &   di_value_i(size_j_values)
c
      double precision, allocatable :: diterms_thrds(:,:),
     &                                 iiterms_thrds(:,:,:)
c
      integer :: i, j, skipped_killed, ierr, blk
      integer, save :: message_count=0, message_count2=0
      logical ::  chk_killed
      double precision, parameter :: zero=0.0d0, one=1.0d0,
     &                               two=2.0d0, toler_static_j=1.0d-5
c
c          if running parallel MPI:
c             alert all workers to enter this routine, then send all
c             data required for J-I calculations that is only on root.
c
      call wmpi_alert_slaves( 27 )
      call wmpi_send_jint
c
      if( debug_driver ) write(iout,*) ' >>> entered element driver'
c
c          keep track of ring_count for storing and retrieving output
c
      ring_count = ring_count + 1
c
c          initialize totals for domain and set symmetry factor.
c          print necessary column headers. use arrays to reduce
c          per thread contibutions for J and II terms. saves
c          lots of atomic updates inside block driver
c          routine.
c
      diterms          = zero ! all
      iiterms          = zero !  "
      skipped_killed   = 0
      message_count    = 0
      message_count2   = 0
      symm_factor      = one
      if( symmetric_domain ) symm_factor = two
      if( myid .eq. 0 ) then
         if( print_elem_values ) write(iout,9100)
         if( print_elem_values ) write(iout,9110)
      end if
c
      allocate( diterms_thrds(size_j_values,num_threads),
     &   iiterms_thrds(size_i_values,size_i_values,num_threads) )
      diterms_thrds = zero
      iiterms_thrds = zero
c
c          loop over all element blocks. inside block checks
c          if element is in current domain
c
c$OMP PARALLEL DO PRIVATE( now_thread, blk ) ! all else shared
      do blk = 1, nelblk
         if( elblks(2,blk) .ne. myid ) cycle
         now_thread = omp_get_thread_num() + 1
         call dicmj_do_a_blk( blk, now_thread,
     &      diterms_thrds(1,now_thread), iiterms_thrds(1,1,now_thread),
     &      message_count, message_count2  )
      end do
c$OMP END PARALLEL DO
c
c          reduce all thread contributions to J and II terms
c
      do i = 1, num_threads
        diterms = diterms + diterms_thrds(:,i)
        iiterms = iiterms + iiterms_thrds(:,:,i)
      end do
c
      call wmpi_reduce_vec( diterms, size_j_values )
      call wmpi_reduce_vec( iiterms(1,1), size_i_values**2 )
      call wmpi_redint( skipped_killed )
      call wmpi_redlog( face_loading )
c
      if( worker_processor ) then
        if( allocated( q_values ) ) deallocate( q_values )
        if( allocated( q_element_maps ) ) deallocate( q_element_maps )
        if( allocated( front_element_list ) )
     &      deallocate( front_element_list )
        if( allocated( expanded_front_nodes ) )
     &      deallocate( expanded_front_nodes )
        return
      end if
c
c             done with this domain computation. print values
c             as required based on user specified flags.
c
c             only rank 0 runs this code. see deallocate/return
c             above for workers
c
c             J-integral results:
c
      if( comput_j ) call dicmj_finalize_j_terms
c
c             interaction-integral results:
c             calculate K from relationship with interaction integral
c
c             the 8 terms for I-integral results are stored
c             in array e_iresults(8,7) as follows:
c
c                               value    auxiliary field
c
c                  iiterms(i,1): KI       plane stress
c                  iiterms(i,2): KI       plane stress
c                  iiterms(i,3): KII      plane stress
c                  iiterms(i,4): KII      plane stress
c                  iiterms(i,5): KIII     anti-plane shear
c                  iiterms(i,6): T11,T33  plane stress
c                  iiterms(i,7): T11,T33  plane strain
c                  iiterms(i,8): T13      anti-plane strain
c
      if( comput_i ) call dicmj_finalize_i_terms
c
c             compute J-values from K-values...
c
c             including KI plane stress, KII plane stress, and KIII
c
      j_from_ks(ring_count,1) = i_storage(ring_count,11,1)**2 / e_front
     &                        + i_storage(ring_count,11,3)**2 / e_front
     &                        + i_storage(ring_count,11,5)**2
     &                        * (one + nu_front ) / e_front
c
c             including KI plane strain, KII plane strain, and KIII
c
      j_from_ks(ring_count,2) = i_storage(ring_count,11,2)**2
     &                          * (one - nu_front**2 ) / e_front
     &                        + i_storage(ring_count,11,4)**2
     &                          * (one - nu_front**2 ) / e_front
     &                        + i_storage(ring_count,11,5)**2
     &                        * (one + nu_front ) / e_front
c
c             print sum of values from all elements in domain
c
      call dicmj_print_values
c
      return
c
 9100 format(/,1x,25x,'domain integral components',
     &  /,      1x,25x,'--------------------------')
 9110 format(7x,'element',5x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8' )
c
      contains
c     ========
c
c***************************************************************
c                                                              *
c      dicmj -> contains dicmj_finalize_j_terms                *
c                                                              *
c        written by: rhd                                       *
c        last modified: 3/30/2018 rhd                          *
c                                                              *
c***************************************************************
c
      subroutine dicmj_finalize_j_terms
      implicit none
c
      jtotal = diterms(1) + diterms(2) + diterms(3) + diterms(4) +
     &         diterms(5) + diterms(6) + diterms(7) + diterms(8)
      if( abs(diterms(3)+diterms(4)).lt.toler_static_j*abs(jtotal) )
     &     static_j = .true.
      static_total = diterms(1) + diterms(2) + diterms(5)
     &             + diterms(6) + diterms(7) + diterms(8)
      di_value_j   = jtotal / front_q_area
      static_di    = static_total / front_q_area
      domain_min_j = min( domain_min_j, di_value_j )
      domain_max_j = max( domain_max_j, di_value_j )
      domain_avg_j = domain_avg_j + di_value_j
      static_min   = min( static_min, static_di  )
      static_max   = max( static_max, static_di  )
      static_avg   = static_avg + static_di
c
c             calculate stress intensity factor K, from J, for:
c
c               1. plane stress where loading is pure mode I or mode II.
c               2. plane strain where loading is pure mode I or mode II.
c               3. anti-plane shear where loading is pure mode III.
c
c             if J is negative, calculate K from |J|, and then change sign.
c
      sign = di_value_j / abs(di_value_j)
c
      ks_from_j(ring_count,1) = sqrt(abs(di_value_j) * e_front)*sign
c
      ks_from_j(ring_count,2) =
     &     sqrt(abs(di_value_j) * e_front/(one - nu_front**2))*sign
c
      ks_from_j(ring_count,3) =
     &     sqrt(abs(di_value_j) * e_front/(one + nu_front)) * sign
c
c          store J-integral values for output in didriv.f
c
      j_storage(ring_count,1)   = dble( nowring )
      j_storage(ring_count,2:9) = diterms(1:8)
      j_storage(ring_count,10)  = di_value_j
      j_storage(ring_count,11)  = dble(skipped_killed)
c
      return
      end subroutine dicmj_finalize_j_terms
c
c***************************************************************
c                                                              *
c      dicmj -> contains dicmj_finalize_i_terms                *
c                                                              *
c        written by: rhd                                       *
c        last modified: 3/30/2018 rhd                          *
c                                                              *
c***************************************************************
c
      subroutine dicmj_finalize_i_terms
      implicit none
c
      do j = 1, 8
         itotal(j) = iiterms(1,j) + iiterms(2,j) + iiterms(3,j)
     &             + iiterms(4,j) + iiterms(5,j) + iiterms(6,j)
     &             + iiterms(7,j) + iiterms(8,j)
         di_value_i(j) = itotal(j) / front_q_area
c
c          calculate stress intensity factors and T-stress
c
         Ki  = zero
         T11 = zero
         T33 = zero
         T13 = zero
c
c          KI (plane stress auxiliary fields)
c
         if( j .eq. 1 ) Ki = di_value_i(j) * e_front / two
c
c          KI (plane strain auxiliary fields)
c
         if( j .eq. 2 ) Ki = di_value_i(j) * e_front
     &                     / ( two * (one - nu_front**2) )
c
c          KII (plane stress auxiliary fields)
c
         if( j .eq. 3 ) Ki = di_value_i(j) * e_front / two
c
c          KII (plane strain auxiliary fields)
c
         if( j .eq. 4 ) Ki = di_value_i(j) * e_front
     &                     / ( two * (one - nu_front**2) )
c
c          KIII (anti-plane shear auxiliary fields)
c
         if( j .eq. 5 ) Ki = di_value_i(j) * e_front
     &                     / (two * (one + nu_front))
c
c          T11, T33 (plane stress auxiliary fields)
c
         if( j .eq. 6 ) then
            T11 = e_front * di_value_i(j)
            Ki  = T11
            T33 = zero
            domain_min_i(9) = T33
            domain_max_i(9) = T33
            domain_avg_i(9) = T33
         end if
c
c          T11, T33 (plane strain auxiliary fields)
c
         if( j .eq. 7 ) then
           T11 = e_front / (one - nu_front**2 )
     &         * ( di_value_i(j) + nu_front * e33_front )
           T33 = e33_front * e_front + nu_front * T11
           Ki  = T11
           domain_min_i(10) = min( domain_min_i(10), T33 )
           domain_max_i(10) = max( domain_max_i(10), T33 )
           domain_avg_i(10) = domain_avg_i(10) + T33
         end if
c
c             T13 (anti-plane shear auxiliary fields)
c
         if( j .eq. 8 ) then
            T13 = di_value_i(j) * e_front / (two * (one + nu_front ))
            Ki  = T13
         end if
c
c             determine maximum and minimum values
c
         domain_min_i(j) = min( domain_min_i(j), Ki )
         domain_max_i(j) = max( domain_max_i(j), Ki )
         domain_avg_i(j) = domain_avg_i(j) + Ki
c
c             store i-integral values for output in didriv.f
c
         i_storage(ring_count,1,j)   = dble( nowring )
         i_storage(ring_count,2:9,j) = iiterms(1:8,j)
         i_storage(ring_count,10,j)  = di_value_i(j)
         i_storage(ring_count,11,j)  = Ki
         if( j.eq.7 ) i_storage(ring_count,12,j) = T33
         i_storage(ring_count,13,j)  = dble(skipped_killed)
      end do  !  on i
c
      return
      end subroutine dicmj_finalize_i_terms
c***************************************************************
c                                                              *
c      dicmj -> contains dicmj_print_values                    *
c                                                              *
c        written by: rhd                                       *
c        last modified: 3/30/2018 rhd                          *
c                                                              *
c***************************************************************
c
      subroutine dicmj_print_values
      implicit none


      if ( print_elem_values ) then
         write(iout,9160)
         if( comput_j ) write(iout,9170) "J     ", (diterms(i),i=1,8),
     &        jtotal
         if( comput_i ) then
            write(iout,9170) "I_KI  ", (iiterms(i,1),i=1,8), itotal(1)
            write(iout,9170) "I_KI  ", (iiterms(i,2),i=1,8), itotal(2)
            write(iout,9170) "I_KII ", (iiterms(i,3),i=1,8), itotal(3)
            write(iout,9170) "I_KII ", (iiterms(i,4),i=1,8), itotal(4)
            write(iout,9170) "I_KIII", (iiterms(i,5),i=1,8), itotal(5)
            write(iout,9170) "I_T11 ", (iiterms(i,6),i=1,8), itotal(6)
            write(iout,9170) "I_T11 ", (iiterms(i,7),i=1,8), itotal(7)
            write(iout,9170) "I_T13 ", (iiterms(i,8),i=1,8), itotal(8)
         end if
      end if
c
c             print sum of values from all elements in domain,
c             and J, I for domain
c
      if( print_totals ) then   ! should not be used now (3/9/2018, rhd)
         write(iout,9180)
         if( nowring .eq. 0 ) then
            if( comput_j .and. comput_i ) write(iout,9185)
            if( comput_j .and. .not. comput_i ) write(iout,9186)
            if( comput_i .and. .not. comput_j ) write(iout,9187)
            if( comput_j ) then
               write(iout,9220) "J     ", domain_id(1:8),
     &              (diterms(j),j=1,8), di_value_j, skipped_killed
            end if
            if( comput_i ) then
               write(iout,9220) "I_KI  ", domain_id(1:8),
     &              (iiterms(i,1),i=1,8), di_value_i(1), skipped_killed
               write(iout,9220) "I_KI  ", domain_id(1:8),
     &              (iiterms(i,2),i=1,8), di_value_i(2), skipped_killed
               write(iout,9220) "I_KII ", domain_id(1:8),
     &              (iiterms(i,3),i=1,8), di_value_i(3), skipped_killed
               write(iout,9220) "I_KII ", domain_id(1:8),
     &              (iiterms(i,4),i=1,8), di_value_i(4), skipped_killed
               write(iout,9220) "I_KIII", domain_id(1:8),
     &              (iiterms(i,5),i=1,8), di_value_i(5), skipped_killed
               write(iout,9220) "I_T11 ", domain_id(1:8),
     &              (iiterms(i,6),i=1,8), di_value_i(6), skipped_killed
               write(iout,9220) "I_T11 ", domain_id(1:8),
     &              (iiterms(i,7),i=1,8), di_value_i(7), skipped_killed
               write(iout,9220) "I_T13 ", domain_id(1:8),
     &              (iiterms(i,8),i=1,8), di_value_i(8), skipped_killed
            end if
         end if
c
         if( nowring.gt.0 ) then
            if( comput_j .and.       comput_i ) write(iout,9250)
            if( comput_j .and. .not. comput_i ) write(iout,9270)
            if( comput_i .and. .not. comput_j ) write(iout,9290)
            if( comput_j ) then
               write(iout,9300) "J     ", nowring, (diterms(i),i=1,8),
     &              di_value_j, skipped_killed
            end if
            if( comput_i ) then
               write(iout,9300) "I_KI  ", nowring,
     &              (iiterms(i,1),i=1,8), di_value_i(1), skipped_killed
               write(iout,9300) "I_KI  ", nowring,
     &              (iiterms(i,2),i=1,8), di_value_i(2), skipped_killed
               write(iout,9300) "I_KII ", nowring,
     &              (iiterms(i,3),i=1,8), di_value_i(3), skipped_killed
               write(iout,9300) "I_KII ", nowring,
     &              (iiterms(i,4),i=1,8), di_value_i(4), skipped_killed
               write(iout,9300) "I_KIII", nowring,
     &              (iiterms(i,5),i=1,8), di_value_i(5), skipped_killed
               write(iout,9300) "I_T11 ", nowring,
     &              (iiterms(i,6),i=1,8), di_value_i(6), skipped_killed
               write(iout,9300) "I_T11 ", nowring,
     &              (iiterms(i,7),i=1,8), di_value_i(7), skipped_killed
               write(iout,9300) "I_T13 ", nowring,
     &              (iiterms(i,8),i=1,8), di_value_i(8), skipped_killed
            end if
         end if
      end if
c
      return
c
 9100 format(/,1x,25x,'domain integral components',
     &  /,      1x,25x,'--------------------------')
 9110 format(7x,'element',5x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8' )
 9130 format(/,5x,' gauss point strains for element ',i7,', gpt ',i2,
     &       ':',/, 3(10x,3(e11.4,2x),/))
 9140 format(1x,a,1x,i7,8(1x,e11.4))
 9150 format('')
 9160 format(/,1x,'element totals:',3x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',8x,'total')
 9170 format(1x,a,7x,9(1x,e11.4))
 9180 format(/,1x,25x,'domain integral components',
     &  /,      1x,25x,'--------------------------')
 9185 format(8x,'domain',8x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',5x,'total J,I',
     &    2x,'killed ele' )
 9186 format(8x,'domain',8x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',7x,'total J',
     &    2x,'killed ele' )
 9187 format(8x,'domain',8x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',7x,'total I',
     &    2x,'killed ele' )
 9220 format(1x,a,1x,a8,9(1x,e11.4),2x,'(',i3,')')
 9250 format(8x,'domain',6x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',5x,'total J,I',
     &    2x,'killed ele' )
 9270 format(8x,'domain',6x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',7x,'total J',
     &    2x,'killed ele' )
 9290 format(8x,'domain',6x,'dm1',9x,'dm2',9x,'dm3',9x,'dm4',
     &    9x,'dm5',9x,'dm6',9x,'dm7',9x,'dm8',7x,'total I',
     &    2x,'killed ele' )
 9300 format(1x,a,1x,i7,9(1x,e11.4),2x,'(',i3,')')
 9400 format(/,1x,'area under q-function along crack front:  ',e11.4,
     &       /,1x,'length along crack front for this domain: ',e11.4,/)
c
      end subroutine dicmj_print_values
c
      end subroutine dicmj

c     ****************************************************************
c     *                                                              *
c     *                      dicmj_do_a_blk                          *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/23/2018 rhd              *
c     *                                                              *
c     *      Processing 1 domain and output values. Run MPI          *
c     *      parallel over (mpi) domains, threads over blocks        *
c     *      in each rank. In large models, many (MPI) domains       *
c     *      and element blocks have no elements within the J-II     *
c     *      domain and are skipped                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine dicmj_do_a_blk( blk, now_thread, diterms,
     &                           iiterms, message_count,
     &                           message_count2 )
c
      use global_data, only : myid, iout=>out, noelem, nonode, bits,
     a                        iprops, props, lprops, elblks, mxndel,
     b                        nstrs, mxgp, mxelpr, nstr,
     c                        worker_processor, numprocs,
     d                        scoord=>c, sdispl=>u, sveloc=>v,
     e                        saccel=>a
      use main_data, only   : incmap, elems_to_blocks, incid,
     a                        fgm_node_values, eq_node_force_indexes,
     b                        trn, trnmat,
     c                        temper_elems, temper_nodes,
     d                        temper_nodes_ref, eq_node_forces,
     e                        eq_node_force_len
      use elem_block_data, only : rot_n1_blocks,
     a                            urcs_n_blocks, cdest_blocks,
     b                            edest_blocks, eps_n_blocks,
     c                            initial_state_data
      use j_data, only : debug_driver, symmetric_domain,
     a  print_elem_values, q_element_maps, q_values, fgm_e, fgm_nu,
     b  block_seg_curves, seg_snode_nu, swd_at_nodes, strain_at_nodes,
     c  ignore_face_loads, omit_crack_front_elems, front_element_list,
     d  comput_i, cf_tractions, domain_rot, debug_elements,
     e  one_point_rule, e_front, nu_front, comput_j,
     f  cf_traction_flags, face_loading, process_temperatures,
     g  snode_alpha_ij, seg_snode_e, displ_grad_at_nodes,
     h  size_j_values, size_i_values, front_list_length,
     i  j_geonl_formulation, j_linear_formulation,
     j  process_initial_state, temperatures_on_model
c
      implicit none
c
c          parameters
c
      integer :: blk, now_thread, message_count, message_count2
      double precision :: diterms(size_j_values),
     &                    iiterms(size_i_values,size_i_values)
c
c          local declarations. includes pointers to simplify
c          addressing into blocks.
c
      double precision ::
     &   e_coord(3,mxndel), e_displ(3,mxndel), e_vel(3,mxndel),
     &   e_accel(3,mxndel), q_element(mxndel),
     &   e_stress(nstrs*mxgp), e_rots(9*mxgp), e_force(3,mxndel),
     &   symm_factor, e_node_temps(mxndel),
     &   elem_uniform_temp, e_alpha_ij(6,mxndel), eq_load_modifier,
     &   dummy, elem_nod_strains(6,mxndel), elem_nod_swd(mxndel),
     &   elem_nod_displ_derivs(9,mxndel), e_strain(9,mxgp), sign,
     &   e(mxndel), nu(mxndel), cf_load(3), e_W_is(mxgp)
c
      double precision :: e_jresults(size_j_values),
     &   e_iresults(size_i_values,size_i_values)
c
      integer :: eps_offset, gpn, node, node_id, etype,
     &           skipped_killed, ierr, elemno, num_enodes,
     &           num_gpts, rel_elem, i, j, i_elem,
     &           offset, element, numrow_sig, span, felem
      integer :: snodes(mxndel)
      integer, dimension (:,:), pointer :: edest, cdest
      real :: e_props(mxelpr)
      logical, external :: dibmck
      logical ::  geonl, chk_killed, omit_J7_J8_front_elems, ok,
     &            elem_temps, face_loads_possible, is_front_element,
     &            front_elem_flag, seg_curves_flag, brick, linearform,
     &            process_swd_derivs, process_strain_derivs,
     &            process_grad_derivs, process_blk
      double precision, parameter :: zero=0.0d0, one=1.0d0,
     &                               two=2.0d0, toler_static_j=1.0d-5
c
      span  = elblks(0,blk)
      felem = elblks(1,blk)
c
c          loop over all elements in this block. is there work to
c          do: (1) element is in domain, and (2) not killed
c
      process_blk = .false.
c
      do i_elem = 1, span
         elemno = felem + i_elem - 1
         if( .not. dibmck( elemno, q_element_maps, bits ) ) cycle
         if( chk_killed( elemno ) ) cycle
         process_blk = .true.
         exit
      end do
      if( .not. process_blk ) return
c
      symm_factor      = one
      if( symmetric_domain ) symm_factor = two
      process_swd_derivs    = size( swd_at_nodes ) == nonode
      process_strain_derivs = size( strain_at_nodes, 2 ) == nonode
      process_grad_derivs   = size( displ_grad_at_nodes, 2 ) == nonode
c
      do i_elem = 1, span
c
c          handle only elements in J-II domain
c
         elemno = felem + i_elem - 1
         if( .not. dibmck( elemno, q_element_maps, bits ) ) cycle
         if( chk_killed( elemno ) ) then
c$OMP ATOMIC UPDATE
            skipped_killed = skipped_killed + 1
            cycle
         end if
c
         etype      = iprops(1,elemno)
         brick      = etype .ge. 1  .and. etype .le. 5
         if( .not. brick ) then
c$OMP ATOMIC UPDATE
             message_count = message_count + 1
             if( message_count <=10 ) then
               write(iout,9200) elemno
               if( message_count == 10 ) write(iout,9210)
             end if
             cycle
         end if
c
         num_enodes = iprops(2,elemno)
         num_gpts   = iprops(6,elemno)
         geonl      = lprops(18,elemno)
         linearform = .not. geonl
         ok = .true.
         if( j_linear_formulation ) then  ! element must match domain
            if( geonl ) ok = .false.
         end if
         if( j_geonl_formulation ) then
            if( linearform ) ok = .false.
         end if
         if( .not. ok ) then
c$OMP ATOMIC UPDATE
             message_count2 = message_count2 + 1
             if( message_count2 <=10 ) then
               write(iout,9230) elemno
               if( message_count2 == 10 ) write(iout,9240)
             end if
             cycle
         end if
c
         rel_elem   = elemno - felem + 1
         cdest      => cdest_blocks(blk)%ptr
         edest      => edest_blocks(blk)%ptr
c
c          build copy of element coordinates, displacements,
c          velocities, accelerations, q-values from global data.
c          zero load terms that correspond to constrained dof.
c          if we have temperature loadings, get temps for element
c          nodes (node values + element uniform value). rotate
c          displacements, velocities and accelerations from
c          constraint-compatible to global coordinates. obtain
c          element nodal values of stress work density and displ
c          gradients structure node-value arrays.
c
         call dicmj_element_setup
c
c          crack-face tractions. the contribution to J arising
c          from crack-face tractions is calculated element-by-element
c          using the equivalent nodal loads on each element face
c          computed during setup of load step for applied pressures,
c          surface tractions and body forces. here we just get
c          the element equiv. nodal forces. they are passed to
c          element J routine which sorts out the loaded face.
c
c          a similar procedure is used to compute the contribution
c          to I.
c
         face_loads_possible = .not. ignore_face_loads .and.
     &                               eq_node_force_len .gt. 0
         if( face_loads_possible ) then
            if( eq_node_force_indexes(elemno) .ne. 0 ) then
               call vec_ops( e_force,  !    copy routine
     &              eq_node_forces(eq_node_force_indexes(elemno)),
     &              dummy, 3*num_enodes, 5 )
            end if
         end if
c
         is_front_element = .false.
         do i = 1, front_list_length
           if( front_element_list(i) == elemno ) then
             is_front_element = .true.
             exit
            end if
         end do
         omit_J7_J8_front_elems = omit_crack_front_elems .and.
     &                            is_front_element
c
c          gather element stresses
c
         offset = (rel_elem-1) * nstrs * num_gpts
         associate( urcs_n =>  urcs_n_blocks(blk)%ptr )
         do i = 1, nstrs * num_gpts
            e_stress(i) = urcs_n(offset+i)
         end do
         end associate

c          gather element properties
c
         e_props(1:mxelpr) = props(1:mxelpr,elemno)
c
c          gather 3x3 rotation matrices from polar decompositions
c          at each gauss point of element for geonl
c
         if( geonl ) then
           offset =  (rel_elem-1) * 9 * num_gpts
           associate( rots => rot_n1_blocks(blk)%ptr )
           do i = 1, 9 * num_gpts
             e_rots(i) = rots(offset+i)
           end do
           end associate
         end if
c
         e_strain = zero ! all terms
c
c          for I computations,
c          gather element strains at integration points.
c          and store in tensor form. convert engineering
c          (total) strains of off-diagonal terms to tensor strains.
c
c          for crack-face loading, make a copy of the tractions
c          input by the user in the domain definition. if none
c          was input, dielem_c.f automatically uses the equivalent
c          nodal loads used to solve the boundary value problem.
c
         if( comput_i ) call dicmj_more_i_setup
c
c          set flag for crack-front element.
c
         front_elem_flag = is_front_element
c
c          call the element routine to compute contribution
c          to the j-integral and or i-integral for domain
c
         element    = elemno  !  protect variables
         numrow_sig = nstrs   !       "
         call dielem( e_coord, q_element, e_displ, e_vel, e_accel,
     a                e_force, e_stress, e_props, e_props, e_props,
     b                e_rots, domain_rot, e_jresults, e_node_temps,
     c                elem_temps, e_alpha_ij, ierr, element, geonl,
     d                numrow_sig, snodes, e_W_is, elem_nod_swd,
     e                elem_nod_strains, elem_nod_displ_derivs,
     f                omit_J7_J8_front_elems, fgm_e,
     f                fgm_nu, e_strain, e, e_front, nu, nu_front,
     g                e_iresults, cf_traction_flags, cf_load,
     h                front_elem_flag, seg_curves_flag )
c
         if( comput_j ) then
            do i = 1, size_j_values
               e_jresults(i) = e_jresults(i) * symm_factor
               diterms(i)    = diterms(i) + e_jresults(i)
            end do
            if( print_elem_values .and. comput_j )
     &        write(iout,9140) "J     ", elemno,
     &                         e_jresults(1:size_j_values)
         end if
c
         if( comput_i ) then
            do i = 1, size_i_values
               do j = 1, size_i_values
                  e_iresults(i,j) = e_iresults(i,j) * symm_factor
                  iiterms(i,j)    = iiterms(i,j) + e_iresults(i,j)
               end do
            end do
         end if
c
         if( comput_i .and. symmetric_domain ) then
            iiterms(1:size_i_values,3:5) = zero
            iiterms(1:size_i_values,8)   = zero
         end if
c
         if( print_elem_values .and. comput_i ) then
            write(iout,9140) "I_KI  ", elemno, (e_iresults(i,1),i=1,8)
            write(iout,9140) "I_KI  ", elemno, (e_iresults(i,2),i=1,8)
            write(iout,9140) "I_KII ", elemno, (e_iresults(i,3),i=1,8)
            write(iout,9140) "I_KII ", elemno, (e_iresults(i,4),i=1,8)
            write(iout,9140) "I_KIII", elemno, (e_iresults(i,5),i=1,8)
            write(iout,9140) "I_T11 ", elemno, (e_iresults(i,6),i=1,8)
            write(iout,9140) "I_T11 ", elemno, (e_iresults(i,7),i=1,8)
            write(iout,9140) "I_T13 ", elemno, (e_iresults(i,8),i=1,8)
            write(iout,9150)
         end if
      end do  ! **** on all elements in block ****


      return
 9140 format(1x,a,1x,i7,8(1x,e11.4))
 9150 format('')
 9200 format(/,'>> Warning: only brick elements supported for domain',
     &  ' integrals. element: ',i8,
     & /,    '           element skipped.')

 9210 format(/,'>> Warning: *** no more messages about bricks ***')
 9230 format(/,'>> Warning: all elements in domain must have the same',
     &  ' formulation (small strain, larger strain. element: ',i8,
     & /,    '           element skipped.')
 9240 format(/,'>> Warning: *** no more messages about formulation',
     &   ' mismatch ***')


      contains
c     ========
c
c***************************************************************
c                                                              *
c      dicmj_do_a_blk -> contains dicmj_element_setup          *
c                                                              *
c        written by: rhd                                       *
c        last modified: 6/23/2018 rhd                          *
c                                                              *
c***************************************************************
c
      subroutine dicmj_element_setup
      implicit none
c
      integer :: incpos, k1, k2, i1, i2, i3, enode, snode
      logical, parameter :: local_debug = .true.
c
      seg_curves_flag = .false.
c
      incpos     = incmap(elemno)
      k1         = num_enodes
      k2         = 2*num_enodes
c
!DIR$ VECTOR ALIGNED
      do enode = 1, num_enodes
         e_coord(1,enode) = scoord(cdest(enode,rel_elem))
         e_coord(2,enode) = scoord(cdest(enode+k1,rel_elem))
         e_coord(3,enode) = scoord(cdest(enode+k2,rel_elem))
         snode            = incid(incpos+enode-1)
         snodes(enode)    = snode
         q_element(enode) = q_values(snode)
         i1               = edest(enode,   rel_elem)
         i2               = edest(k1+enode,rel_elem)
         i3               = edest(k2+enode,rel_elem)
         e_displ(1,enode) = sdispl(i1)
         e_displ(2,enode) = sdispl(i2)
         e_displ(3,enode) = sdispl(i3)
         e_vel(1,enode)   = sveloc(i1)
         e_vel(2,enode)   = sveloc(i2)
         e_vel(3,enode)   = sveloc(i3)
         e_accel(1,enode) = saccel(i1)
         e_accel(2,enode) = saccel(i2)
         e_accel(3,enode) = saccel(i3)
         e(enode)         = props(7,elemno)
         nu(enode)        = props(8,elemno)
         if( fgm_e )  e(enode)  = fgm_node_values(snode,1)
         if( fgm_nu ) nu(enode) = fgm_node_values(snode,2)
         if( temperatures_on_model .and. 
     &       allocated( block_seg_curves ) ) then
            if( block_seg_curves(blk) ) then
               seg_curves_flag = .true.
               e(enode)        = seg_snode_e(snode)
               nu(enode)       = seg_snode_nu(snode)
            end if
         end if
c
c         if necessary, transform nodal values from
c         constraint-compatible coordinates to global
c         coordinates.
c
         if( trn(snode) ) then
            call di_trans_nodvals( e_displ(1,enode),
     &                             trnmat(snode)%mat(1,1))
            call di_trans_nodvals( e_vel(1,enode),
     &                             trnmat(snode)%mat(1,1))
            call di_trans_nodvals( e_accel(1,enode),
     &                             trnmat(snode)%mat(1,1))
         end if
c
         e_force(1,enode)          = zero
         e_force(2,enode)          = zero
         e_force(3,enode)          = zero
         e_node_temps(enode)       = zero
         e_alpha_ij(1:6,enode)     = zero
         elem_nod_swd(enode)       = zero
         elem_nod_strains(1:6,enode) = zero
         elem_nod_displ_derivs(1:9,enode) = zero
c
      end do ! on enode
c

      if( process_swd_derivs ) then
!DIR$ VECTOR ALIGNED
        do enode = 1, num_enodes
          snode = incid(incpos+enode-1)
          elem_nod_swd(enode) = swd_at_nodes(snode)
        end do
      end if
c
      if( process_strain_derivs ) then
        do enode = 1, num_enodes
          snode = incid(incpos+enode-1)
!DIR$ VECTOR ALIGNED
          elem_nod_strains(1:6,enode) = strain_at_nodes(1:6,snode)
        end do
      end if
c
      if( process_grad_derivs ) then
        do enode = 1, num_enodes
          snode = incid(incpos+enode-1)
!DIR$ VECTOR ALIGNED
          elem_nod_displ_derivs(1:9,enode) =
     &               displ_grad_at_nodes(1:9,snode)
        end do
      end if
c
c        for problems with temperature loads, pull out
c        current temperatures of element nodes and set the thermal
c        expansion coefficients at element nodes. remove reference
c        temperatures of model nodes. the global temperatures flag
c        refers to increment temps for step just completed. cannot use
c        it since early application of temps may be followed by only
c        mechanical loads
c
      elem_temps = .false.
      if( temperatures_on_model ) then
        elem_uniform_temp = temper_elems(elemno)
!DIR$ VECTOR ALIGNED
        do enode = 1, num_enodes
         snode = incid(incpos + enode - 1)
         e_node_temps(enode) = temper_nodes(snode) + elem_uniform_temp
     &                        - temper_nodes_ref(snode)
         if( abs(e_node_temps(enode)) > zero ) elem_temps = .true.
        end do
        if( allocated( snode_alpha_ij ) ) then
         do enode = 1, num_enodes
          snode = incid( incpos + enode - 1 )
!DIR$ VECTOR ALIGNED
          e_alpha_ij(1:6,enode) = snode_alpha_ij(snode,1:6)
         end do
        end if
      end if
c
c        setup vector of energy density values at element
c        integration points at user-defined initial state
c        if required. otherwise set zero vector to eliminate
c        logic in element computational routines.
c
      if( process_initial_state ) then
        do gpn = 1, num_gpts
         e_W_is(gpn) =
     &     initial_state_data(blk)%W_plastic_nis_block(i_elem,gpn)
        end do

        if( local_debug .and.
     &      (elemno == 8955 .or. elemno == 9060 ) ) then
          write(iout,*) ' @ 1. dicmj_element_setup. elem: ',elemno
          do gpn = 1, num_gpts
           write(iout,9000)  gpn, e_W_is(gpn)
          end do
        end if
      else
          e_W_is = zero
      end if
 9000 format(10x,i3,f15.8)
c
      return
      end subroutine dicmj_element_setup
c
c***************************************************************
c                                                              *
c      dicmj -> contains dicmj_more_i_setup                    *
c                                                              *
c        written by: rhd                                       *
c        last modified: 3/30/2018 rhd                          *
c                                                              *
c***************************************************************
c
      subroutine dicmj_more_i_setup
      implicit none

c
c          gather element strains at integration points.
c          and store in tensor form. convert engineering
c          (total) strains of off-diagonal terms to tensor strains.
c
      eps_offset = (rel_elem - 1) * nstr * num_gpts
      associate( eps => eps_n_blocks(blk)%ptr )
c
!DIR$ VECTOR ALIGNED
      do gpn = 1, num_gpts
         e_strain(1,gpn) = eps(eps_offset + 1)
         e_strain(2,gpn) = eps(eps_offset + 4)/two
         e_strain(3,gpn) = eps(eps_offset + 6)/two
         e_strain(4,gpn) = eps(eps_offset + 4)/two
         e_strain(5,gpn) = eps(eps_offset + 2)
         e_strain(6,gpn) = eps(eps_offset + 5)/two
         e_strain(7,gpn) = eps(eps_offset + 6)/two
         e_strain(8,gpn) = eps(eps_offset + 5)/two
         e_strain(9,gpn) = eps(eps_offset + 3)
         eps_offset      = eps_offset + nstr
      end do
c
      end associate
c
      if( debug_driver ) then
         do gpn = 1, num_gpts
            write(iout,9130) elemno, gpn, (e_strain(i,gpn),i=1,9)
         end do
      end if
c
c          for crack-face loading, make a copy of the tractions
c          input by the user in the domain definition. if none
c          was input, dielem_c.f automatically uses the equivalent
c          nodal loads used to solve the boundary value problem.
c
      cf_load(1) = cf_tractions(1)
      cf_load(2) = cf_tractions(2)
      cf_load(3) = cf_tractions(3)
c
      return
 9130 format(/,5x,' gauss point strains for element ',i7,', gpt ',i2,
     &       ':',/, 3(10x,3(e11.4,2x),/))

      end subroutine  dicmj_more_i_setup
c
      end subroutine dicmj_do_a_blk


c***************************************************************
c                                                              *
c subroutine to transform a 3x1 vector of constraint-          *
c compatible-coordinate-system values to global-coordinate-    *
c system values. the routine premultiplies the incoming        *
c constraint-compatible values by the transpose of the stored  *
c global-to-constraint coordinate system rotation matrix.      *
c                                                              *
c                written by: mcw                               *
c                last modified: 3/26/2018 rhd                  *
c                                                              *
c***************************************************************
c
      subroutine di_trans_nodvals( vec, transmat )
      implicit none
      double precision :: vec(3), transmat(3,3)
      double precision :: tempvec(3)
c
      tempvec(1) = transmat(1,1) * vec(1) + transmat(2,1) * vec(2)
     &           + transmat(3,1) * vec(3)
      tempvec(2) = transmat(1,2) * vec(1) + transmat(2,2) * vec(2)
     &           + transmat(3,2) * vec(3)
      tempvec(3) = transmat(1,3) * vec(1) + transmat(2,3) * vec(2)
     &           + transmat(3,3) * vec(3)
c
      vec(1) = tempvec(1)
      vec(2) = tempvec(2)
      vec(3) = tempvec(3)
c
      return
      end
