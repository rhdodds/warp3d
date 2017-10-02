c
c **********************************************************************
c *                                                                    *
c *    mm10.f                                                          *
c *                                                                    *
c *         original code written by : mcm                             *
c *         extensively revised by Tim Truster and R. Dodds            *
c *                                                                    *
c *         Stress/strain update routines for                          *
c *         crystal plasticity material modeled via Beaudoin et al.    *
c *                                                                    *
c **********************************************************************
c

c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine mm10                         *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/8/2016 rhd               *
c     *                                                              *
c     *              crystal plasticity stress-strain update         *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10( gp, span, ncrystals, hist_sz, history_n,
     &                 history_np1, local_work, uddt, gp_temps,
     &                 gp_temp_inc, iout, display_matl_messages,
     &                 do_nonlocal, nonlocal_state, maxnonlocal,
     &                 iter_0_extrapolate_off )
c
      use segmental_curves, only: max_seg_points
      use mm10_defs
      use mm10_constants
c
      implicit none
      include 'include_sig_up'
c
c                 parameter definitions
c
      integer, intent(in) :: gp, span, hist_sz, iout, maxnonlocal,
     &                       ncrystals(mxvl)
      logical, intent(in) :: display_matl_messages, do_nonlocal,
     &                       iter_0_extrapolate_off
c
      double precision, intent(in) ::
     &      uddt(mxvl,nstr), gp_temps(mxvl), gp_temp_inc(mxvl)
      double precision, intent(inout) ::
     &      history_n(span,hist_sz), history_np1(span,hist_sz),
     &      nonlocal_state(mxvl,maxnonlocal)
c
c                 locals
c
      type(crystal_props) :: cc_props
      type(crystal_state) :: cc_n, cc_np1
c
      integer :: i, c, now_element, iloop, number_crystals, crys_no
      double precision :: sig_avg(6), se(6), p_strain_ten_c(6),
     &                    p_strain_ten(6),
     &                    tang_avg(6,6), tang_avg_vec(36), ! see equiv
     &                    slip_avg(length_comm_hist(5))
      double precision :: t_work_inc, p_work_inc,p_strain_inc,
     &                    n_avg, p_strain_avg
      logical :: debug, gpall, locdebug, mat_debug
      equivalence( tang_avg, tang_avg_vec )
c
      debug = .false.
c
      if( debug ) write (iout,*) ".... entering mm10"
c
c              we have data passed for integration point # gp for
c              all (span) elements in this block. the CP updating
c              routines process only a single int point. loop to
c              process that point for all elements in block.
c
c              for small strain analysis, rot_blk_n1 set to identity
c              now by warp3d.
c
      local_work%material_cut_step = .false.
c
      do i = 1, span
c
        iloop = i ! protect index i
        now_element = local_work%felem + i - 1
        if( debug ) write(iout,*) 'Updating element: ', now_element
c
c              for step = 1, init element history
c
        if( local_work%step .eq. 1 ) then
           if( debug ) write(iout,*) "Init GP history"
           call mm10_init_general_hist( span, history_n(iloop,1) )
           call mm10_init_uout_hist( span, history_n(iloop,1) )
           call mm10_init_slip_hist( span, history_n(iloop,1) )
        end if
c
        sig_avg      = zero  ! 6 x 1
        slip_avg     = zero  ! length_comm_hist(5) x 1
        t_work_inc   = zero
        p_work_inc   = zero
        p_strain_inc = zero
        tang_avg_vec = zero  ! 36 x 1
        p_strain_ten = zero
        n_avg        = zero
c
        number_crystals = ncrystals(i)
c
c              loop on all crystals at integration point
c
        do c = 1, number_crystals
          crys_no   = c  ! protect index c
          debug     = local_work%debug_flag(i)
          locdebug  = .false.
          mat_debug = locdebug
c
c              initialize G,H arrays for certain CP constitutive
c              models. Note: all crystals in the element block
c              will have the same hardening model and slip systems
c
      call mm10_set_cons( local_work, cc_props, 1, iloop, c )
          if( locdebug ) write(iout,*) "Setting up properties"
          call mm10_a_do_crystal
          if( local_work%material_cut_step ) then
            return
          end if
c              release allocated G & H arrays for
c              certain CP constitutive models
c
      call mm10_set_cons( local_work, cc_props, 2, iloop, c )
        end do ! over crystals
c
c              finalize averages over all crystals at point.
c              p_strain_avg -> average effective strain increment
c                              from the average tensor
c
        call mm10_a_crystal_avgs
c
c              nonlocal state values returned are creep rate
c              wrt real time and total creep strain;
c              only for first crystal
c               nonlocal(1) -> effective creep rate wrt real time
c               nonlocal(2) -> (total) effective creep strain
c               nonlocal(3) -> effective creep eponent (n)
c               nonlocal(4) -> effective Norton constant B
c
        if( do_nonlocal ) call mm10_a_compute_local
c
c              store results for integration point into
c              variables passed by warp3d
c
        call mm10_a_store_crystal
c
      end do ! over span
c
c
      return
c
      contains ! for mm10
c
c              ****************************************************
c              *  contains: mm10_a_crystal_avgs                   *
c              ****************************************************
c
      subroutine mm10_a_crystal_avgs
      implicit none
c
      double precision :: rncry, s_trace, t1, t2

      rncry        = dble( ncrystals(iloop) )
      sig_avg      = sig_avg / rncry         ! vector
      tang_avg_vec = tang_avg_vec / rncry    ! vector
      slip_avg     = slip_avg / rncry        ! vector
      t_work_inc   = t_work_inc / rncry      ! scalar
      p_work_inc   = p_work_inc / rncry      ! scalar
      p_strain_inc = p_strain_inc / rncry    ! scalar
      p_strain_ten = p_strain_ten / rncry    ! scalar
      t1           = p_strain_ten(1)**2 + p_strain_ten(2)**2 +
     &               p_strain_ten(3)**2
      t2           = p_strain_ten(4)**2 + p_strain_ten(5)**2 +
     &               p_strain_ten(6)**2
      s_trace = (two/three) * ( t1 + t2/two )
c
      if( (exponent(t1) .lt. -400) .or. (exponent(t1) .gt. 600)  .or.
     &    (s_trace .eq. zero) ) then
        p_strain_avg = zero
        n_avg = zero
      else
        p_strain_avg = sqrt( s_trace )
        n_avg        = n_avg / p_strain_avg / rncry
      end if
      if( locdebug ) write(iout,*) "stress", sig_avg
      if( locdebug ) write(iout,*) "tang", tang_avg
c
      return
      end subroutine mm10_a_crystal_avgs
c
c              ****************************************************
c              *  contains: mm10_a_do_crystal                     *
c              ****************************************************
c
      subroutine mm10_a_do_crystal
      implicit none
c
c              we use some work vectors to eliminate
c              overhead of compiler generated temporaries
c              and to-from copying
c
      integer :: len, sh2, sh3, eh2, eh3
      double precision :: work_vec1(9), work_vec2(6), work_gradfe(27),
     &                    work_R(9), user_initial_stresses(6)
c
      call mm10_init_cc_props( local_work%c_props(iloop,c),
     &              local_work%angle_type(iloop),
     &              local_work%angle_convention(iloop),
     &              local_work%debug_flag(iloop), cc_props )
      cc_props%out = iout
c
      if( local_work%step .eq. 1 ) then
         user_initial_stresses(1:6) =
     &            local_work%urcs_blk_n(iloop,1:6,gp)
         call mm10_init_cc_hist0( cc_props,
     &           local_work%c_props(iloop,c)%init_angles(1),
     &           history_n(iloop,1), user_initial_stresses,
     &           span, crys_no, hist_sz )
      end if
c
      if( locdebug ) write(iout,*) "Copying n to struct"
      sh2  = indexes_common(2,1)
      eh2  = indexes_common(2,2)
      sh3  = indexes_common(3,1)
      eh3  = indexes_common(3,2)
      work_gradfe(1:27) = history_n(iloop,sh2:eh2)
      work_R(1:9) = history_n(iloop,sh3:eh3)
      call mm10_copy_cc_hist( crys_no, span, history_n(iloop,1),
     &         work_gradfe,  work_R, ! both readonly,
     &         cc_props, cc_n )
c
      work_vec1(1:9) = local_work%rot_blk_n1(iloop,1:9,gp)
      work_vec2(1:6) = uddt(iloop,1:6)
      call mm10_setup_np1(
     &        work_vec1, work_vec2, ! read only in subroutine
     &        local_work%dt, gp_temps(iloop), local_work%step,
     &        iloop-1+local_work%felem, local_work%iter,
     &        local_work%gpn, cc_np1 )
c
      if( locdebug ) write(iout,*) "Updating crystal ", c
      call mm10_solve_crystal( cc_props, cc_np1, cc_n,
     &        local_work%material_cut_step, iout, .false., 0,
     &        p_strain_ten_c, iter_0_extrapolate_off )
      if( local_work%material_cut_step ) then
          call mm10_set_cons( local_work, cc_props, 2, i, c )
          return
      end if
c
c                  accumulate sums for subsequent averaging
c                    cp_strain_inc -> effective plastic increment
c                    p_strain_ten -> plastic strain increment tensor
c                    n_avg -> effective creep exponent
c
      sig_avg      = sig_avg + cc_np1%stress     ! 6x1 vector
      tang_avg     = tang_avg + cc_np1%tangent   ! 6x6 matrix
      len = length_comm_hist(5)
      slip_avg(1:len) = slip_avg(1:len) + cc_np1%slip_incs(1:len)
      t_work_inc   = t_work_inc + cc_np1%work_inc
      p_work_inc   = p_work_inc + cc_np1%p_work_inc
      p_strain_inc = p_strain_inc + cc_np1%p_strain_inc
      p_strain_ten = p_strain_ten + p_strain_ten_c
      n_avg        = n_avg + cc_np1%p_strain_inc*cc_np1%u(12)
c
c                  store the CP history for this crystal
c
      call mm10_store_cryhist( crys_no, span, cc_props, cc_np1,
     &                         cc_n, history_np1(iloop,1) )
c
      return
c
      end subroutine mm10_a_do_crystal
c
c              ****************************************************
c              *  contains: mm10_a_store_crystal                  *
c              ****************************************************
c
      subroutine mm10_a_store_crystal
      implicit none
c
      integer :: sh, eh, sh3, eh3, sh4, eh4, sh5,eh5, len
      sh   = indexes_common(1,1)
      eh   = indexes_common(1,2)
      sh3  = indexes_common(3,1)
      eh3  = indexes_common(3,2)
      sh4  = indexes_common(4,1)
      eh4  = indexes_common(4,2)
      sh5  = indexes_common(5,1)
      eh5  = indexes_common(5,2)
      len  = eh5 - sh5 + 1
c
      local_work%urcs_blk_n1(iloop,1:6,gp) = sig_avg(1:6)
      history_np1(iloop,sh3:eh3) = local_work%rot_blk_n1(iloop,1:9,gp)
      history_np1(iloop,sh:eh) = tang_avg_vec(1:36)
      history_np1(iloop,sh5:eh5) = history_n(iloop,sh5:eh5) +
     &                             slip_avg(1:len)
      local_work%urcs_blk_n1(iloop,7,gp) =
     &                           local_work%urcs_blk_n(iloop,7,gp)
     &                           + t_work_inc
      local_work%urcs_blk_n1(iloop,8,gp) =
     &                           local_work%urcs_blk_n(iloop,8,gp)
     &                           + p_work_inc
      local_work%urcs_blk_n1(iloop,9,gp) =
     &                           local_work%urcs_blk_n(iloop,9,gp)
     &                           + p_strain_inc
      history_np1(iloop,sh4+0) = history_n(iloop,sh4+0) + t_work_inc
      history_np1(iloop,sh4+1) = history_n(iloop,sh4+1) + p_work_inc
      history_np1(iloop,sh4+2) = history_n(iloop,sh4+2) + p_strain_inc
c
      return
      end subroutine mm10_a_store_crystal

c
c              ****************************************************
c              *  contains: mm10_a_compute_local                  *
c              ****************************************************
c
      subroutine mm10_a_compute_local
      implicit none
c
      double precision :: s_trace, se(6), t1, t2, B_eff
c
      nonlocal_state(iloop,1) = p_strain_avg / local_work%dt
      nonlocal_state(iloop,2) = local_work%urcs_blk_n(iloop,9,gp) +
     &                          p_strain_avg
      if( n_avg .lt. one ) then ! limit range of values
        nonlocal_state(iloop,3) = one
      elseif( n_avg .gt. ten ) then
        nonlocal_state(iloop,3) = ten
      else
        nonlocal_state(iloop,3) = n_avg
      endif
c
      s_trace = (sig_avg(1) + sig_avg(2) + sig_avg(3)) / three
      se(1:6) = sig_avg(1:6)
      se(1)   = se(1) - s_trace
      se(2)   = se(2) - s_trace
      se(3)   = se(3) - s_trace
      t1      = se(1)**2 + se(2)**2 + se(3)**2
      t2      = se(4)**2 + se(5)**2 + se(6)**2
      s_trace = sqrt( onept5 * ( t1 + two*t2 ) )
      B_eff   = p_strain_avg / local_work%dt/(s_trace**n_avg)
      nonlocal_state(iloop,4) = B_eff
c
      return
      end subroutine mm10_a_compute_local

      end subroutine mm10
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_set_cons                     *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 9/30/2016 rhd               *
c     *                                                              *
c     *       set interaction matrices G & H for slip system type    *
c     *       used for this element block                            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine mm10_set_cons( local_work, cc_props, isw, i, c )
c
      use mm10_defs ! to get definition of cc_props
      use mm10_constants
c
      implicit none
      include 'include_sig_up'
c
      type(crystal_props) :: cc_props
      integer :: isw, s_type1, n_hard, h_type, allocate_status
      integer :: i, c
      logical :: process_G_H
c
      h_type  = local_work%c_props(i, c)%h_type
      s_type1 = local_work%c_props(i, c)%s_type
      n_hard  = local_work%c_props(i, c)%num_hard
c
      process_G_H = ( h_type .eq. 4 ) .or.
     &              ( h_type .eq. 7 ) .or.
     &              ( h_type .eq. 8 ) .or.
     &              ( h_type .eq. 9 )
      if( .not. process_G_H ) return
c
      if( (h_type .eq. 4) .or. (h_type .eq. 7) ) then
       select case( isw )
         case( 1 ) ! allocate G,H interaction matrices
           allocate( cc_props%Gmat(n_hard,n_hard),
     &               cc_props%Hmat(n_hard,n_hard),
     &               stat=allocate_status)
           if( allocate_status .ne. 0 ) then
              write(*,*) ' error allocating G matrix'
              call die_gracefully
           end if
           call mm10_mrr_GH( s_type1, n_hard, cc_props%Gmat,
     &                     cc_props%Hmat,i,c)
         case( 2 ) ! deallocate G,H matrices
           deallocate( cc_props%Gmat, cc_props%Hmat )
         case default
            write(*,*) '>>>> FATAL ERROR. invalid isw, mm10_set_cons'
            write(*,*) '                  job terminated'
            call die_abort
         end select
         return
      end if
c       Armstrong-Frederick
      if( h_type .eq. 8 ) then
       select case( isw )
         case( 1 ) ! allocate G=q,H=many_params matrices
           allocate( cc_props%Gmat(n_hard,n_hard),
     &               cc_props%Hmat(7,n_hard), stat=allocate_status)
           if( allocate_status .ne. 0 ) then
              write(*,*) ' error allocating G matrix'
              call die_gracefully
           end if
           call mm10_DJGM_GH(local_work,s_type1,n_hard,
     &               cc_props%Gmat, cc_props%Hmat,i,c)
         case( 2 ) ! deallocate G,H matrices
           deallocate( cc_props%Gmat, cc_props%Hmat )
         case default
           write(*,*) '>>>> FATAL ERROR. invalid isw, mm10_set_cons'
           write(*,*) '                  job terminated'
           call die_abort
       end select
      end if
c
      if( h_type .eq. 9 ) then
       select case( isw )
         case( 1 ) ! allocate G=q,H=many_params matrices
           allocate( cc_props%Gmat(n_hard,n_hard),
     &               cc_props%Hmat(7,n_hard), stat=allocate_status)
           if( allocate_status .ne. 0 ) then
              write(*,*) ' error allocating G matrix'
              call die_gracefully
           end if
           call mm10_DJGM_GH(local_work,s_type1,n_hard,
     &               cc_props%Gmat, cc_props%Hmat,i,c)
         case( 2 ) ! deallocate G,H matrices
           deallocate( cc_props%Gmat, cc_props%Hmat )
         case default
           write(*,*) '>>>> FATAL ERROR. invalid isw, mm10_set_cons'
           write(*,*) '                  job terminated'
           call die_abort
       end select
      end if
c
      return
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_init_cc_hist0                *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 9/30/2017 rhd               *
c     *                                                              *
c     *    Initialize a crystal history variables                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_init_cc_hist0( props, angles, history,
     &                               user_initial_stresses, span,
     &                               crys_no, hist_size )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      integer :: span, crys_no, hist_size
      double precision :: history(span,*), angles(*),
     &                    user_initial_stresses(*)
c
      integer :: a, sh, eh, sh2, eh2, len1, len2
      double precision :: work_hist1(hist_size),
     &                    work_hist2(hist_size) ! automatics
c
c              Stress
c
      sh = index_crys_hist(crys_no,1,1)
      eh = index_crys_hist(crys_no,1,2)
      history(1,sh:eh) = zero
      history(1,sh+0) = user_initial_stresses(1)
      history(1,sh+1) = user_initial_stresses(2)
      history(1,sh+2) = user_initial_stresses(3)
      history(1,sh+3) = user_initial_stresses(4)
      history(1,sh+4) = user_initial_stresses(5)
      history(1,sh+5) = user_initial_stresses(6)
c
c              Angles
c
      sh = index_crys_hist(crys_no,2,1)
      history(1,sh+0) = angles(1)
      history(1,sh+1) = angles(2)
      history(1,sh+2) = angles(3)
c
c              Rotation
c
      sh = index_crys_hist(crys_no,3,1)
      history(1,sh+0) = one
      history(1,sh+1) = zero
      history(1,sh+2) = zero
      history(1,sh+3) = zero
      history(1,sh+4) = one
      history(1,sh+5) = zero
      history(1,sh+6) = zero
      history(1,sh+7) = zero
      history(1,sh+8) = one
c
c              D
c
      sh = index_crys_hist(crys_no,4,1)
      eh = index_crys_hist(crys_no,4,2)
      history(1,sh:eh) = zero
c
c              eps
c
      sh = index_crys_hist(crys_no,5,1)
      eh = index_crys_hist(crys_no,5,2)
      history(1,sh:eh) = zero
c
c              slip_incs
c
      sh = index_crys_hist(crys_no,6,1)
      eh = index_crys_hist(crys_no,6,2)
      history(1,sh:eh) = zero
c
c              Hardening
c
      sh   = index_crys_hist(crys_no,7,1)
      eh   = index_crys_hist(crys_no,7,2)
      sh2  = index_crys_hist(crys_no,8,1)
      eh2  = index_crys_hist(crys_no,8,2)
      len1 = eh - sh + 1
      len2 = eh2 - sh2 + 1
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type ) ! hardening model
        case( 1 ) ! simple Voche
          call mm10_init_voche( props, work_hist1, work_hist2 )
        case( 2 ) ! MTS
          call mm10_init_mts( props, work_hist1, work_hist2 )
        case( 3 ) ! User
         call mm10_init_user( props, work_hist1, work_hist2 )
        case( 4 ) ! ORNL
         call mm10_init_ornl( props, work_hist1, work_hist2 )
        case( 7 ) ! mrr
         call mm10_init_mrr( props, work_hist1, work_hist2 )
        case( 8 ) ! Armstrong-Frederick
         call mm10_init_arfr( props, work_hist1, work_hist2 )
        case( 9 ) ! DJGM
         call mm10_init_djgm( props, work_hist1, work_hist2 )
        case default
         call mm10_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
      history(1,sh:eh)   = work_hist1(1:len1)
      history(1,sh2:eh2) = work_hist2(1:len2)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_unknown_hard_error           *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 1/27/14                     *
c     *                                                              *
c     *     A common error message for the general hardening setup   *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_unknown_hard_error( props )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
c
      write(props%out,101) props%h_type
 101  format(
     &      10x,'>> Error: unknown hardening type ', 'i6', '.',
     &    /,10x, 'Aborting...')
      call die_gracefully
c
      end
cc
c **********************************************************************
c *                                                                    *
c *         warp3d called routines for sizes                           *
c *                                                                    *
c **********************************************************************
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_set_sizes_special            *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 06/16/2016 tjt              *
c     *                                                              *
c     *    called by warp3d for each material model to obtain        *
c     *    various sizes of data for the model                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_set_sizes_special( size_data, local_el  )
      use main_data, only : imatprp
      use mm10_defs, only : one_crystal_hist_size, common_hist_size
      use mm10_constants
      use global_data ! old common.main
      implicit none
c
      integer :: size_data(*), local_el, matnum, ncrystals
c
c        size_data(1)  :  no. of words of history data for each
c                         integration point
c        in this case sizeof(__)*number of crystals
c
      matnum    = iprops(38,local_el)
      ncrystals = imatprp(101,matnum)
c
c              total history size is going to be:
c
      size_data(1) = common_hist_size +
     &               ncrystals*(one_crystal_hist_size)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_tangent                      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/10/2016 rhd              *
c     *                                                              *
c     *     get consistent tangent after a converged stress update.  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_tangent( props, np1, n, gaspt, Jmat, nrJmat,
     &                         ncJmat )
      use mm10_defs
      use mm10_constants
      implicit none
c
c                 parameter definitions
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      integer :: gaspt, nrJmat, ncJmat  ! nr=nc = 6 + props%num_hard
      double precision, dimension(nrJmat,ncJmat) :: Jmat
c
c
c                 locals
c
      integer :: len, i, info, gpp
      logical :: debug, gpall, locdebug
      double precision, dimension(6,6) :: J11, JJ, JR, JA, JB, JK
      double precision, dimension(6) :: d_mod, d_barp, tw, work_vec
      double precision, dimension(6,props%num_hard) :: ed, J12
      double precision, dimension(props%num_hard,6) :: J21
      double precision, dimension(props%num_hard,12) :: beta

c
c                 automatic vectors-arrays
c

      integer, dimension(props%num_hard) :: ipiv
      double precision, dimension(6,props%nslip) :: symtqmat
      double precision, dimension(6,props%nslip) :: dgammadd
      double precision,
     &    dimension(props%num_hard,props%num_hard) :: J22, alpha
      double precision :: alpha1
c!DIR$ ASSUME_ALIGNED Jmat:64
c
      debug    = .false.     ! props%debug
      gpall    = props%gpall ! print iteration norms for all int points
      gpp      = props%gpp   ! set one particular G.P. to print
      locdebug = .false.
c
      call mm10_a_zero_vec( np1%tangent, 36 )
c
c              pull Jacobian from solver evaluation
c
      len              = props%num_hard
      J11(1:6,1:6)     = Jmat(1:6,1:6)
      J12(1:6,1:len)   = Jmat(1:6,7:6+len)
      J21(1:len,1:6)   = Jmat(7:6+len,1:6)
      J22(1:len,1:len) = Jmat(7:6+len,7:6+len)
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type  )
        case( 1 )
          call mm10_ed_voche( props, np1, n, np1%stress, np1%tau_tilde,
     &                        ed )
        case( 2 )
          call mm10_ed_mts( props, np1, n, np1%stress, np1%tau_tilde,
     &                        ed )
        case( 3 )
          call mm10_ed_user( props, np1, n, np1%stress, np1%tau_tilde,
     &                        ed )
        case( 4 )
          call mm10_ed_ornl( props, np1, n, np1%stress, np1%tau_tilde,
     &                        ed )
        case( 7 )
          call mm10_ed_mrr( props, np1, n, np1%stress, np1%tau_tilde,
     &                        ed )
        case( 8 )
          call mm10_ed_arfr( props, np1, n, np1%stress, np1%tau_tilde,
     &                        ed )
        case( 9 )
          call mm10_ed_djgm( props, np1, n, np1%stress, np1%tau_tilde,
     &                        ed )
        case default
          call mm10_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
      call mm10_symSWmat( np1%stress, np1%qc, props%nslip, symtqmat )
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type )
        case( 1 )
          call mm10_dgdd_voche( props, np1, n, np1%stress,
     &                          np1%tau_tilde, np1%D, dgammadd )
        case( 2 )
         call mm10_dgdd_mts( props, np1, n, np1%stress,
     &                          np1%tau_tilde, np1%D, dgammadd )
        case( 3 )
         call mm10_dgdd_user( props, np1, n, np1%stress,
     &                          np1%tau_tilde, np1%D, dgammadd )
        case( 4 )
         call mm10_dgdd_ornl( props, np1, n, np1%stress,
     &                          np1%tau_tilde, np1%D, dgammadd )
        case( 7 )
         call mm10_dgdd_mrr( props, np1, n, np1%stress,
     &                          np1%tau_tilde, np1%D, dgammadd )
        case( 8 )
         call mm10_dgdd_arfr( props, np1, n, np1%stress,
     &                          np1%tau_tilde, np1%D, dgammadd )
        case( 9 )
         call mm10_dgdd_djgm( props, np1, n, np1%stress,
     &                          np1%tau_tilde, np1%D, dgammadd )
        case default
          call mm10_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
      call mm10_a_zero_vec( JA, 36 )
      do i = 1, props%nslip
c              work_vec = stiffness * ms + 2 *  symtqmat
       call mm10_a_mult_type_4( work_vec, props%stiffness,
     &                          np1%ms(1,i), symtqmat, two  )
c              JA = JA +  work_vec * trans(dgammadd(:,i))
       call mm10_a_mult_type_5( JA, work_vec, dgammadd(1,i) )
      end do
c
c              compute tangent matrix T where
c                     T = JJ^-1*JR
c                    JR = (Cijkl - JA - JB)
c                    JJ = J11 - J12*J22^-1*J21
c                    JB = J12*J22^-1*ed
c
      call mm10_a_copy_vector( JJ, J11, 36 )  ! JJ = J11 (6 x 6)
      call mm10_a_copy_vector( alpha, J22, len*len )   ! len x len
      call mm10_a_copy_vector( beta, J21, len*6 )  ! beta = J21
      call mm10_a_copy_vector( beta(1,7), ed, len*6 )
      if( locdebug ) write (*,*) "beta", beta(1:6,1:12) ! beta = ed
c
c              compute J22^-1*J21 and J22^-1*ed
c
      call DGESV( len, 2*6, alpha, len, ipiv, beta, len, info )
c
c              JJ = JJ - J12*beta(*,1:6)   [ JJ 6x6 ]
c
      if( locdebug ) write (*,*) "beta", beta(1:6,1:12)
      call DGEMM ('N','N', 6, 6, len, -one, J12, 6, beta, len,
     &             one, JJ, 6 )
c
c
      call mm10_a_zero_vec( JB, 36 )
      call DGEMM ('N','N', 6 ,6, len, one, J12, 6, beta(1,7),
     &             len, one, JB, 6 )  ! make JB 6 x
      if( debug ) write (*,*) "JB", JB(1,2)
c
      JR = props%stiffness - JA - JB  ! 6 x 6
c
c              np1%tangent = matmul(JJ, JR)
c              avoid explicitly computing the inverse
c
      call DGESV( 6, 6, JJ, 6, ipiv, JR, 6, info )
      call mm10_a_copy_vector( np1%tangent, JR, 36 )
      if( locdebug ) write (*,*) "JR", JR(1:6,1:6)
c
      return
      end subroutine

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_setup                        *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 9/29/2016 rhd               *
c     *                                                              *
c     *     setup hardening for a particular state np1               *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_setup( props, np1, n )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
c              locals
c
      integer :: i, j, k, s, t
      double precision :: RE(6,6), rw(3,3), rwc(3,3), trans_nRp(3,3),
     &                    curv(3,3,3), cn(3), tm(3), const, t1, t2,
     &                    work(3,3), work_vec(3)
      double precision, parameter :: alpha = 1.0d0/3.0d0, ! geometric hardening
     &                               twothirds = 2.0d0/3.0d0
c
c              effective strain increment
c
      t1 = np1%D(1)*np1%D(1) + np1%D(2)*np1%D(2) + np1%D(3)*np1%D(3)
      t2 = np1%D(4)*np1%D(4) + np1%D(5)*np1%D(5) + np1%D(6)*np1%D(6)
      np1%dg = dsqrt(twothirds*(t1+half*t2))
c
c              compute current m and q tensors
c              yes, these are supposed to be transposes.  We actually
c              need the backwards rotation from the lattice state.
c
      trans_nRp(1,1) = n%Rp(1,1)
      trans_nRp(2,1) = n%Rp(1,2)
      trans_nRp(3,1) = n%Rp(1,3)
      trans_nRp(1,2) = n%Rp(2,1)
      trans_nRp(2,2) = n%Rp(2,2)
      trans_nRp(3,2) = n%Rp(2,3)
      trans_nRp(1,3) = n%Rp(3,1)
      trans_nRp(2,3) = n%Rp(3,2)
      trans_nRp(3,3) = n%Rp(3,3)
c
      call mm10_RT2RVE( trans_nRp, RE ) ! makes RE(6,6)
      call mm10_RT2RVW( trans_nRp, RW ) ! makes RW(3,3)
      call mm10_a_mult_type_1( work, np1%R(1,1), trans_nRp )
      call mm10_RT2RVW( work, RWC ) ! makes RWC(3,3)
c
      do i = 1, props%nslip
        call mm10_a_mult_type_2( np1%ms(1,i), RE, props%ms(1,i) )
        call mm10_a_mult_type_3( np1%qs(1,i), RW, props%qs(1,i) )
        call mm10_a_mult_type_3( np1%qc(1,i), RWC, props%qs(1,i) )
      end do
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type )
        case( 1 ) ! voche
          call mm10_setup_voche(props, np1, n)
        case( 2 ) ! MTS
          call mm10_setup_mts(props, np1, n)
        case( 3 ) ! user
          call mm10_setup_user(props, np1, n)
        case( 4 ) ! ORNL
          call mm10_setup_ornl(props, np1, n)
        case( 7 )  ! MRR
          call mm10_setup_mrr(props, np1, n)
        case( 8 )  ! Armstrong-Frederick
          call mm10_setup_arfr(props, np1, n)
        case( 9 )  ! DJDM
          call mm10_setup_DJGM(props, np1, n)
        case default ! error
          call mm10_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
c
c              compute quatities related to backstress, gradients, etc.
c
      select case( props%h_type )

        case( 1, 2, 3 ) ! voche, MTS, user
c
c                  tau lambdas for geometric hardening
c                  lattice curvature
c
!DIR$ IVDEP
          do i = 1, 3
            curv(i,1,1) = n%gradFeinv(i,1,1) - n%gradFeinv(i,1,1)
            curv(i,1,2) = n%gradFeinv(i,1,2) - n%gradFeinv(i,2,1)
            curv(i,1,3) = n%gradFeinv(i,1,3) - n%gradFeinv(i,3,1)
            curv(i,2,1) = n%gradFeinv(i,2,1) - n%gradFeinv(i,1,2)
            curv(i,2,2) = n%gradFeinv(i,2,2) - n%gradFeinv(i,2,2)
            curv(i,2,3) = n%gradFeinv(i,2,3) - n%gradFeinv(i,3,2)
            curv(i,3,1) = n%gradFeinv(i,3,1) - n%gradFeinv(i,1,3)
            curv(i,3,2) = n%gradFeinv(i,3,2) - n%gradFeinv(i,2,3)
            curv(i,3,3) = n%gradFeinv(i,3,3) - n%gradFeinv(i,3,3)
          end do
c
c                  Acharya's large strain definition of lambda
c
          const = props%k_0 * props%burgers * alpha * alpha *
     &        np1%mu_harden * np1%mu_harden / two / props%theta_0
c
          do t = 1, props%nslip
c
c                  normal into current coordinates (This could be
c                  the problem...)
c
           call mm10_a_mult_type_3( work_vec, trans_nRp,
     &                              props%ns(1,t) )
           call mm10_a_mult_type_3( cn, n%R(1,1), work_vec )
c
c                  large-strain lambda
c
           tm(1) = zero; tm(2) = zero; tm(3) = zero
!DIR$ IVDEP
           do i = 1, 3
               tm(i) = tm(i) + half * (
     &           - curv(i,1,2) * cn(3)
     &           + curv(i,1,3) * cn(2)
     &           + curv(i,2,1) * cn(3)
     &           - curv(i,2,3) * cn(1)
     &           - curv(i,3,1) * cn(2)
     &           + curv(i,3,2) * cn(1) )
           end do ! on i
           np1%tau_l(t) = const *
     &                 dsqrt( tm(1)*tm(1) + tm(2)*tm(2) + tm(3)*tm(3) )
          end do ! on t
c
        case( 4 ) ! ORNL
        case( 7 ) ! MRR
        case( 8 ) ! Armstrong-Frederick
        case( 9 ) ! DJGM
        case default
          call mm10_unknown_hard_error(props)
      end select
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_store_cryhist                *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/3/2016 rhd               *
c     *                                                              *
c     *          Copy the state np1 struct to the history            *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_store_cryhist( crys_no, span, props, np1,
     &                               n, history )
      use mm10_defs
      use mm10_constants
      implicit none
c
      integer :: span, crys_no
      double precision :: history(span,*)
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      integer :: sh, eh, len1, len2
c!DIR$ ASSUME_ALIGNED history:64

c
      sh = index_crys_hist(crys_no,1,1)
      eh = index_crys_hist(crys_no,1,2)
      history(1,sh:eh) = np1%stress(1:6)
c
c
      sh = index_crys_hist(crys_no,2,1)
      history(1,sh+0) = np1%euler_angles(1)
      history(1,sh+1) = np1%euler_angles(2)
      history(1,sh+2) = np1%euler_angles(3)
c
      sh = index_crys_hist(crys_no,3,1)
      history(1,sh+0) = np1%Rp(1,1)
      history(1,sh+1) = np1%Rp(2,1)
      history(1,sh+2) = np1%Rp(3,1)
      history(1,sh+3) = np1%Rp(1,2)
      history(1,sh+4) = np1%Rp(2,2)
      history(1,sh+5) = np1%Rp(3,2)
      history(1,sh+6) = np1%Rp(1,3)
      history(1,sh+7) = np1%Rp(2,3)
      history(1,sh+8) = np1%Rp(3,3)
c
      sh = index_crys_hist(crys_no,4,1)
      eh = index_crys_hist(crys_no,4,2)
      history(1,sh:eh) = np1%D(1:6)
c
      sh = index_crys_hist(crys_no,5,1)
      eh = index_crys_hist(crys_no,5,2)
      history(1,sh:eh) = np1%eps(1:6)
c
      sh = index_crys_hist(crys_no,6,1)
      eh = index_crys_hist(crys_no,6,2)
      len2 = length_crys_hist(6)
      if( len2 .ne. eh-sh+1 ) then ! sanity check
         write(props%out,9000) 1
         call die_abort
      end if
      history(1,sh:eh) = np1%slip_incs(1:len2)
c
      sh   = index_crys_hist(crys_no,7,1)
      len1 = props%num_hard
      history(1,sh:sh-1+len1) = np1%tau_tilde(1:len1)
c
      sh   = index_crys_hist(crys_no,8,1)
      eh   = index_crys_hist(crys_no,8,2)
      len2 = length_crys_hist(8)
      if( len2 .ne. eh-sh+1 ) then ! sanity check
         write(props%out,9000) 2
         call die_abort
      end if
      history(1,sh:eh) = np1%u(1:len2)
c
      sh = index_crys_hist(crys_no,9,1)
      history(1,sh:sh-1+len1) = np1%tt_rate(1:len1)
c
      sh = index_crys_hist(crys_no,10,1)
      history(1,sh:sh-1+6) = np1%ep(1:6)
c
      sh = index_crys_hist(crys_no,11,1)
      history(1,sh:sh-1+6) = np1%ed(1:6)
c
      return
c
 9000 format('>> FATAL ERROR: mm10_store_cryhist @ ',i2 )

      end


c
c --------------------------------------------------------------------
c
c     Operational subroutines (do not require modification for new
c     constitutive models)
c
c --------------------------------------------------------------------
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_solve_crystal                *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/5/2016 rhd               *
c     *                                                              *
c     *     Advance a crystal from n to np1, store tangent, and      *
c     *     store other output                                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_solve_crystal( props, np1, n, cut, iout, fat,
     &                               gp, p_strain_ten_c,
     &                               iter_0_extrapolate_off )
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      use main_data, only: asymmetric_assembly
      implicit none
c
c              parameters
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      logical :: cut, fat, iter_0_extrapolate_off
      integer :: iout, gp
      double precision :: p_strain_ten_c(6)
c
c              locals -- max_uhard is in param_def
c
      integer :: nrJmat, ncJmat
      double precision :: vec1(max_uhard), vec2(max_uhard)
      double precision, dimension(max_uhard,max_uhard) :: arr1,arr2
      complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
      double precision, allocatable :: Jmat(:,:)
      logical :: no_load

      if( cut ) then ! somehow we got back into this routine
                     ! after one crystal failed already
         write(iout,*) "mm10 returned for crystal solve after
     & another crystal failed!"
            call die_gracefully ! segmentation faults will follow if not stopped here
      end if

      nrJmat = 6 + props%num_hard
      ncJmat = 2 * nrJmat ! extra space needed for Broyden solver
      allocate( Jmat(nrJmat,ncJmat) )

      call mm10_solve_strup( props, np1, n, vec1, vec2, arr1, arr2,
     &   ivec1, ivec2, cut, gp, iter_0_extrapolate_off, no_load,
     &   Jmat, nrJmat, ncJmat )
c
      if( cut ) then  ! global step size reduction needed
        write(iout,*) "mm10 stress update failed"
        deallocate( Jmat )
        return
      end if
c
c              special update with extrapolate off: return stress
c              due to creep & linear elastic stiffness matrix
c              (from mm10_solve_strup)
c
      if( iter_0_extrapolate_off .or. no_load) then
          deallocate( Jmat )
          return
      end if
c
      call mm10_tangent( props, np1, n, gp, Jmat, nrJmat, ncJmat )
c
c              make tangent symmetric
c
      if( .not. asymmetric_assembly .and. .not. fat )
     &   call mm10_a_make_symm_1( np1%tangent )
c
      call mm10_formvecs( props, np1, n, np1%stress, np1%tau_tilde,
     &                    vec1, vec2 )
c
      call mm10_update_rotation( props, np1, n, vec1, vec2 )
c
c              compute and store lots of other state and values for
c              eventual output
c
      call mm10_output( props, np1, n, vec1, vec2, p_strain_ten_c )
c
      deallocate( Jmat )
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_update_euler_angles          *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 9/30/2016 rhd               *
c     *                                                              *
c     *               update euler angles to the new rotation        *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_update_euler_angles( props, np1, n )
      use mm10_defs
      use mm10_constants
c
      implicit none
      type (crystal_props) :: props
      type (crystal_state) :: n, np1
c
      double precision, dimension(3,3) :: full_rot, work1
      double precision :: psiK, phiK, thetaK, psi, phi, theta,
     &                    pps, pms
      double precision, external :: mm10_atan2
      double precision, parameter :: tol=1.0d-16
c
c           Note: This subroutine needs a major fix if
c           ever going to support anything other than degrees+
c           Kocks convention
c
      call mm10_a_mult_type_3t( work1, np1%Rp, np1%R )
      call mm10_a_mult_type_1( full_rot, props%g, work1 )
c
      psiK = mm10_atan2( full_rot(3,2), full_rot(3,1) )
      phiK = mm10_atan2( full_rot(2,3), full_rot(1,3) )
      if( full_rot(3,3) .gt. one ) full_rot(3,3) = one
      thetaK = dacos( full_rot(3,3) )

      if( props%angle_convention .eq. 1 ) then
        psi = psiK
        phi = phiK
        theta = thetaK
      else
        write (props%out,*) "Angle convention not implemented."
        call die_gracefully
      end if
c
      if( props%angle_type .eq. 1) then
        np1%euler_angles(1) = one_eighty/pi*psi
        np1%euler_angles(2) = one_eighty/pi*theta
        np1%euler_angles(3) = one_eighty/pi*phi
      elseif( props%angle_type .eq. 2 ) then
        np1%euler_angles(1) = psi
        np1%euler_angles(2) = theta
        np1%euler_angles(3) = phi
      else
        write (props%out,*) "Unrecognized angle convention."
        call die_gracefully
      end if
c
      return
      end subroutine
c
c              helper for the above, atan2 with range 0 to 2*pi
c
      double precision function mm10_atan2(a, b) ! inlined
      use mm10_constants
      implicit none
      double precision ::  a, b
c
      mm10_atan2 = datan2(a, b)
      if( mm10_atan2 .lt. zero ) mm10_atan2 = mm10_atan2 + two*pi
c
      return
      end function
c

c
c
c
c **********************************************************************
c *                                                                    *
c *                    mm10_invasym                                    *
c *                                                                    *
c *         written by : mcm                                           *
c *         last modified : 9/30/2016 rhd                              *
c *                                                                    *
c *        inversion of a non-symmetric matrix using LAPACK            *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_invasym( A, n )
      implicit none
c
      integer, intent(in) :: n
      double precision, intent(inout), dimension(n,n) :: A

      integer :: info, lwork
      integer, allocatable :: ipivt(:)
      double precision, allocatable :: work(:)
c!DIR$ ASSUME_ALIGNED a:64
c
c              allocate storage, factor, inverse
c
      lwork = n * n
      allocate( ipivt(n), work(lwork) )
      call DGETRF( n, n, A, n, ipivt, info )
      call DGETRI( n, A, n, ipivt, work, lwork, info )
c
c              done. free storage
c
      deallocate( ipivt, work )
c
      return
      end
c
c **********************************************************************
c *                                                                    *
c *            mm10_rotation_matrix                                    *
c *                                                                    *
c *         written by : mcm                                           *
c *         last modified : 9/30/2016 rhd                              *
c *                                                                    *
c *     given euler angles, an angle convention, and an angle type     *
c *     send back the correct rotation matrix.                         *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_rotation_matrix( angles, aconv, atype, r, out )
      use mm10_constants
      implicit none
c
      double precision, dimension(3), intent(in) :: angles
      character, intent(in) :: aconv*5
      character, intent(in) :: atype*7
      integer, intent(in) :: out
c
      double precision, dimension(3,3), intent(out) :: r
      double precision :: a, b, c, psi, theta, phi
c!DIR$ ASSUME_ALIGNED angles:64, r:64
c
      a = angles(1)
      b = angles(2)
      c = angles(3)
c
      if( atype .eq. 'degrees' ) then
         a = a*pi/one_eighty
         b = b*pi/one_eighty
         c = c*pi/one_eighty
      elseif (atype .eq. 'radians') then ! no conversion needed
      else
         write(out,9000)
      end if

      if( aconv .eq. 'kocks' ) then
         psi   = a
         theta = b
         phi   = c
      elseif( aconv .eq. 'bunge' ) then
         psi   = a - pi/two
         theta = b
         phi   = pi/two - c
      elseif( aconv .eq. 'roe' ) then
         psi   = a
         theta = b
         phi   = three*pi/two - c
      else
         write (out,9001)
      end if
c
      r(1,1) = -sin(psi)*sin(phi)-cos(psi)*cos(phi)*cos(theta)
      r(1,2) = cos(psi)*sin(phi)-sin(psi)*cos(phi)*cos(theta)
      r(1,3) = cos(phi)*sin(theta)
      r(2,1) = sin(psi)*cos(phi)-cos(psi)*sin(phi)*cos(theta)
      r(2,2) = -cos(psi)*cos(phi)-sin(psi)*sin(phi)*cos(theta)
      r(2,3) = sin(phi)*sin(theta)
      r(3,1) = cos(psi)*sin(theta)
      r(3,2) = sin(psi)*sin(theta)
      r(3,3) = cos(theta)
c
      return
c
 9000 format(/'danger: unknown angle type passed to rotation_matrix'/)
 9001 format(/'danger: unknown angle convention passed to',
     &        ' rotation_matrix'/)
c
      end
c
c **********************************************************************
c *                                                                    *
c *                    mm10_invsym                                     *
c *                                                                    *
c *         written by : mcm                                           *
c *         last modified : 9/30/2016 rhd                              *
c *                                                                    *
c *         inversion of a symmetric matrix using LAPACK               *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_invsym( A, n )
      implicit none
c
      integer, intent(in) :: n
      double precision, intent(inout), dimension(n,n) :: A
c
      integer :: i, j, info, lwork
      integer, allocatable :: ipiv(:)
      double precision, allocatable :: work(:)
c!DIR$ ASSUME_ALIGNED a:64
c
c              allocate storage, factor, inverse, convert
c              symmetric to full storage
c
      lwork = n * n
      allocate( ipiv(n), work(lwork) )
      call DSYTRF('U', n, A, n, ipiv, work, lwork, info)
      call DSYTRI('U', n, A, n, ipiv, work, info )
c
      do i = 1, n
        do j = 1, i-1
          A(i,j) = A(j,i)
        end do
      end do
c
      deallocate( ipiv, work )
c
      return
      end
c
c **********************************************************************
c *                                                                    *
c *                    mm10_rt2rve                                     *
c *                                                                    *
c *         written by : mcm                                           *
c *         last modified : 9/30/2016 rhd                              *
c *                                                                    *
c *         takes a 3x3 rotation tensor and returns it in a 6x6 form   *
c *         suitable for rotating voigt-type strain vectors            *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_rt2rve( rt, rv )
      use mm10_constants
      implicit none
c
      double precision, dimension(3,3), intent(in) :: rt
      double precision, dimension(6,6), intent(out) :: rv
c
c!DIR$ ASSUME_ALIGNED rt:64, rv:64
c
      rv(1,1) = rt(1,1)**2
      rv(1,2) = rt(1,2)**2
      rv(1,3) = rt(1,3)**2
      rv(1,4) = two*rt(1,1)*rt(1,2)
      rv(1,5) = two*rt(1,3)*rt(1,2)
      rv(1,6) = two*rt(1,1)*rt(1,3)
      rv(2,1) = rt(2,1)**2
      rv(2,2) = rt(2,2)**2
      rv(2,3) = rt(2,3)**2
      rv(2,4) = two*rt(2,1)*rt(2,2)
      rv(2,5) = two*rt(2,3)*rt(2,2)
      rv(2,6) = two*rt(2,1)*rt(2,3)
      rv(3,1) = rt(3,1)**2
      rv(3,2) = rt(3,2)**2
      rv(3,3) = rt(3,3)**2
      rv(3,4) = two*rt(3,1)*rt(3,2)
      rv(3,5) = two*rt(3,3)*rt(3,2)
      rv(3,6) = two*rt(3,1)*rt(3,3)
      rv(4,1) = rt(1,1)*rt(2,1)
      rv(4,2) = rt(1,2)*rt(2,2)
      rv(4,3) = rt(1,3)*rt(2,3)
      rv(4,4) = rt(1,1)*rt(2,2)  + rt(2,1)*rt(1,2)
      rv(4,5) = rt(1,2)*rt(2,3)  + rt(1,3)*rt(2,2)
      rv(4,6) = rt(1,1)*rt(2,3)  + rt(1,3)*rt(2,1)
      rv(5,1) = rt(2,1)*rt(3,1)
      rv(5,2) = rt(3,2)*rt(2,2)
      rv(5,3) = rt(2,3)*rt(3,3)
      rv(5,4) = rt(2,1)*rt(3,2)  + rt(2,2)*rt(3,1)
      rv(5,5) = rt(2,2)*rt(3,3)  + rt(3,2)*rt(2,3)
      rv(5,6) = rt(2,1)*rt(3,3)  + rt(2,3)*rt(3,1)
      rv(6,1) = rt(1,1)*rt(3,1)
      rv(6,2) = rt(1,2)*rt(3,2)
      rv(6,3) = rt(1,3)*rt(3,3)
      rv(6,4) = rt(1,1)*rt(3,2)  + rt(1,2)*rt(3,1)
      rv(6,5) = rt(1,2)*rt(3,3)  + rt(1,3)*rt(3,2)
      rv(6,6) = rt(1,1)*rt(3,3)  + rt(3,1)*rt(1,3)
c
      return
      end subroutine
c
c **********************************************************************
c *                                                                    *
c *                     mm10_rt2rvw                                    *
c *                                                                    *
c *         written by : mcm                                           *
c *         last modified : 10/1/2016 rhd                              *
c *                                                                    *
c *         takes a 3x3 rotation tensor and returns it in a 3x3 form   *
c *         suitable for rotating my 3x1 skew vectors                  *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_rt2rvw( rt, rv )
      implicit none
c
      double precision, dimension(3,3), intent(in) :: rt
      double precision, dimension(3,3), intent(out) :: rv
c!DIR$ ASSUME_ALIGNED rt:64, rv:64
c
      rv(1,1) = rt(2,2)*rt(3,3) - rt(2,3)*rt(3,2)
      rv(1,2) = rt(2,1)*rt(3,3) - rt(2,3)*rt(3,1)
      rv(1,3) = rt(2,1)*rt(3,2) - rt(2,2)*rt(3,1)
      rv(2,1) = rt(1,2)*rt(3,3) - rt(1,3)*rt(3,2)
      rv(2,2) = rt(1,1)*rt(3,3) - rt(1,3)*rt(3,1)
      rv(2,3) = rt(1,1)*rt(3,2) - rt(1,2)*rt(3,1)
      rv(3,1) = rt(1,2)*rt(2,3) - rt(1,3)*rt(2,2)
      rv(3,2) = rt(1,1)*rt(2,3) - rt(1,3)*rt(2,1)
      rv(3,3) = rt(1,1)*rt(2,2) - rt(1,2)*rt(2,1)
c
      return
      end subroutine
c
c **********************************************************************
c *                                                                    *
c *                mm10_et2ev                                          *
c *                                                                    *
c *         written by : mcm                                           *
c *         last modified : 10/1/2016 rhd                              *
c *                                                                    *
c *         strain tensor to strain vector -- note ordering            *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_et2ev( et, ev )
      use mm10_constants
      implicit none
c
      double precision, dimension(3,3), intent(in) :: et
      double precision, dimension(6), intent(out) :: ev
c
c!DIR$ ASSUME_ALIGNED et:64, ev:64
c
      ev(1) = et(1,1)
      ev(2) = et(2,2)
      ev(3) = et(3,3)
      ev(4) = two*et(1,2)
      ev(5) = two*et(2,3)
      ev(6) = two*et(1,3)
c
      return
      end subroutine
c
c **********************************************************************
c *                                                                    *
c *                      mm10_wt2wv                                    *
c *                                                                    *
c *         written by : mcm                                           *
c *         last modified : 10/1/2016 rhd                              *
c *                                                                    *
c *         skew  tensor to skew vector                                *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_wt2wv( wt, wv )
      implicit none
c
      double precision, dimension(3,3), intent(in) :: wt
      double precision, dimension(3), intent(out) :: wv
c!DIR$ ASSUME_ALIGNED wt:64, wv:64
c
      wv(1) = wt(2,3)
      wv(2) = wt(1,3)
      wv(3) = wt(1,2)
c
      return
      end subroutine

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_calc_grads                   *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/12/2016 rhd              *
c     *                                                              *
c     *    calculate the gradient of Re.T (Fe) through a linear fit  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_calc_grads( ngp, elem_type, order, geonl,
     &                            rot_blk, jac, Rps, gradFes, iout,
     &                            now_elem, span, hist_size )
      use mm10_constants
      implicit none
      include 'param_def'
c
c                 parameter definitions
c
      integer :: ngp, elem_type, order, iout, now_elem, span, hist_size
      logical :: geonl
      double precision :: rot_blk(mxvl,9,ngp), Rps(span,hist_size,ngp),
     &                    gradFes(span,hist_size,ngp)
      double precision, dimension(mxvl,3,3) :: jac
c
c                 locals
c
      integer :: i, igp, a, b, lwork, info
      logical :: ok
      double precision :: Rt3d(3,3,mxgp), jacinv(3,3), RHS(mxgp),
     &                    work(8), grads(3,3,3), intermat(mxgp,4),
     &                    onerot(9), local_jac(3,3), local_Rps(9),
     &                    grads_vec(27)
      equivalence( grads, grads_vec )
      double precision :: weight, fact
c
c                 constants
c
c!DIR$ ASSUME_ALIGNED rot_blk:64, Rps:64, gradFes:64, jac:64
c
c                get R components and stick in the right place
c
      if( geonl ) then
        local_jac(1:3,1:3) = jac(now_elem,1:3,1:3)
        call mm10_a_invert_33( jacinv, local_jac, ok )
        do igp = 1, ngp
          onerot = rot_blk(now_elem,1:9,igp)
          local_Rps = Rps(now_elem,1:9,igp)
          call mm10_a_mult_type_3t( Rt3d(1,1,igp), local_Rps, onerot )
        end do
      else
        call mm10_a_copy_vector( Rt3d, Rps, ngp*9 )
      end if
c
c              for each Rt component create an interpolation,
c              solve for the coefficients, and store the gradient
c
      do a = 1, 3
        do b = 1, 3
          intermat = zero
          RHS      = zero
          do i = 1, ngp
c              1-3 are the coordinates
            call getgpts( elem_type, order, i, intermat(i,1),
     &                    intermat(i,2), intermat(i,3), weight )
            intermat(i,4) = one
            RHS(i) = Rt3d(a,b,i)
          end do
c              solve with LAPACK
          lwork = 8
          call DGELS('N',  ngp, 4, 1, intermat, mxgp, RHS, ngp, work,
     &                lwork, info)
          if( info .ne. 0 ) then
            write(iout,9000) 2
            call die_abort
          end if
c
c              get the gradient
c
          grads(a,b,1) = RHS(1)
          grads(a,b,2) = RHS(2)
          grads(a,b,3) = RHS(3)
c
c              take to the current coordinates
c
          if( geonl ) then
c            grads(a,b,1:3) = matmul(jacinv,grads(a,b,1:3))
            grads(a,b,1) = jacinv(1,1)*RHS(1) + jacinv(1,2)*RHS(2) +
     &                     jacinv(1,3)*RHS(3)
            grads(a,b,2) = jacinv(2,1)*RHS(1) + jacinv(2,2)*RHS(2) +
     &                     jacinv(2,3)*RHS(3)
            grads(a,b,3) = jacinv(3,1)*RHS(1) + jacinv(3,2)*RHS(2) +
     &                     jacinv(3,3)*RHS(3)
          end if
        end do
      end do
c
c              flatten and store
c
      do igp = 1, ngp
       gradFes(now_elem,1:27,igp) =  grads_vec(1:27)
      end do
c
      return
9000  format('*** FATAL ERROR: mm10_calc_grads @ ',i2,
     &   /,   '                 execution terminated')
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_setup_np1                    *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 9/29/2016 rhd               *
c     *                                                              *
c     *    Initialize the state np1 structure with new strain/temp/  *
c     *    time increment                                            *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_setup_np1( Rur, dstrain, dt, T, step, elem,
     &                           iter, gp, np1 )
      use mm10_defs
      use mm10_constants
      implicit none
c
      integer :: step, elem, gp, iter
      type(crystal_state) :: np1
      double precision :: Rur(3,3), dstrain(6)
      double precision :: dt, T
c
      integer :: i, j
c
c!DIR$ ASSUME_ALIGNED Rur:64, dstrain:64
c
c              scalars
c
      np1%temp         = T
      np1%tinc         = dt
      np1%step         = step
      np1%elem         = elem
      np1%iter         = iter
      np1%gp           = gp
      np1%dg           = zero
      np1%tau_v        = zero
      np1%tau_y        = zero
      np1%mu_harden    = zero
      np1%work_inc     = zero
      np1%p_work_inc   = zero
      np1%p_strain_inc = zero
c
c              vectors
c
      do i = 1, 6
        np1%D(i)      = dstrain(i)
        np1%stress(i) = zero
        np1%eps(i)    = zero
      end do
c
      np1%euler_angles(1) = zero
      np1%euler_angles(2) = zero
      np1%euler_angles(3) = zero
c
      do i = 1, max_slip_sys
        np1%tau_l(i)     = zero
        np1%slip_incs(i) = zero
      end do
c
      do i = 1, max_uhard
        np1%tau_tilde(i) = zero
        np1%tt_rate(i)   = zero
        np1%u(i)         = zero
      end do
c
c              matrices
c
      do j = 1, 3
        np1%R(1,j)  = Rur(1,j)
        np1%R(2,j)  = Rur(2,j)
        np1%R(3,j)  = Rur(3,j)
        np1%Rp(1,j) = zero
        np1%Rp(2,j) = zero
        np1%Rp(3,j) = zero
      end do
c
      call mm10_a_zero_vec( np1%gradFeinv, 27 )
      call mm10_a_zero_vec( np1%tangent, 36 )
      call mm10_a_zero_vec( np1%ms, 6*max_slip_sys )
      call mm10_a_zero_vec( np1%qs, 3*max_slip_sys )
      call mm10_a_zero_vec( np1%qc, 3*max_slip_sys )
c
      return
c
      end

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_init_uout_hist               *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 9/30/2016                   *
c     *                                                              *
c     *    initialize the user output                                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_init_uout_hist( span, history )
      use mm10_defs, only : indexes_common
      use mm10_constants
c
      implicit none
      integer :: span
      double precision :: history(span,*)
c
      integer :: sh
c!DIR$ ASSUME_ALIGNED history:64
c
      sh  = indexes_common(4,1)
      history(1,sh+0) = zero
      history(1,sh+1) = zero
      history(1,sh+2) = zero
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_init_slip_hist               *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *    initialize the slip totals (output variable)              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_init_slip_hist( span, history )
      use mm10_defs, only : indexes_common
      use mm10_constants
c
      implicit none
      integer :: span
      double precision :: history(span,*)
c
      integer :: sh, eh
c!DIR$ ASSUME_ALIGNED history:64
c
      sh  = indexes_common(5,1)
      eh  = indexes_common(5,2)
      history(1,sh:eh) = zero
c
      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_init_cc_props                *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 9/30/2016 rhd               *
c     *                                                              *
c     *    Copy properties over from local_work into the update      *
c     *    structure                                                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_init_cc_props( inc_props, atype, aconv, debug,
     &                               cc_props )
      use mm10_defs
      use mm10_constants
      implicit none
      include 'include_sig_up'
      integer :: atype, aconv
      logical :: debug
      type(crystal_properties) :: inc_props
      type(crystal_props) :: cc_props
c
c              scalars
c
      cc_props%rate_n      = inc_props%rateN
      cc_props%tau_hat_y   = inc_props%tauHat_y
      cc_props%G_0_y       = inc_props%Go_y
      cc_props%burgers     = inc_props%burgers
      cc_props%p_v         = inc_props%p_v
      cc_props%q_v         = inc_props%q_v
      cc_props%boltzman    = inc_props%boltzman
      cc_props%theta_0     = inc_props%theta_o
      cc_props%eps_dot_0_v = inc_props%eps_dot_o_v
      cc_props%eps_dot_0_y = inc_props%eps_dot_o_y
      cc_props%p_y         = inc_props%p_y
      cc_props%q_y         = inc_props%q_y
      cc_props%tau_a       = inc_props%tau_a
      cc_props%tau_hat_v   = inc_props%tauHat_v
      cc_props%G_0_v       = inc_props%Go_v
      cc_props%k_0         = inc_props%k_o
      cc_props%mu_0        = inc_props%mu_o
      cc_props%D_0         = inc_props%D_o
      cc_props%T_0         = inc_props%t_o
      cc_props%tau_y       = inc_props%tau_y
      cc_props%tau_v       = inc_props%tau_v
      cc_props%voche_m     = inc_props%voche_m
      cc_props%iD_v        = inc_props%iD_v
      cc_props%u1          = inc_props%u1
      cc_props%u2          = inc_props%u2
      cc_props%u3          = inc_props%u3
      cc_props%u4          = inc_props%u4
      cc_props%u5          = inc_props%u5
      cc_props%u6          = inc_props%u6
      cc_props%u7          = inc_props%u7
      cc_props%u8          = inc_props%u8
      cc_props%u9          = inc_props%u9
      cc_props%u10         = inc_props%u10
      cc_props%solver      = inc_props%solver
      cc_props%strategy    = inc_props%strategy
      cc_props%gpall       = inc_props%gpall
      cc_props%gpp         = inc_props%gpp
      cc_props%method      = inc_props%method
      cc_props%miter       = inc_props%miter
      cc_props%atol        = inc_props%atol
      cc_props%atol1       = inc_props%atol1
      cc_props%rtol        = inc_props%rtol
      cc_props%rtol1       = inc_props%rtol1
      cc_props%xtol        = inc_props%xtol
      cc_props%xtol1       = inc_props%xtol1
      cc_props%alter_mode  = inc_props%alter_mode
c
      cc_props%cp_001 = inc_props%cp_001
      cc_props%cp_002 = inc_props%cp_002
      cc_props%cp_003 = inc_props%cp_003
      cc_props%cp_004 = inc_props%cp_004
      cc_props%cp_005 = inc_props%cp_005
      cc_props%cp_006 = inc_props%cp_006
      cc_props%cp_007 = inc_props%cp_007
      cc_props%cp_008 = inc_props%cp_008
      cc_props%cp_009 = inc_props%cp_009
      cc_props%cp_010 = inc_props%cp_010
      cc_props%cp_011 = inc_props%cp_011
      cc_props%cp_012 = inc_props%cp_012
      cc_props%cp_013 = inc_props%cp_013
      cc_props%cp_014 = inc_props%cp_014
      cc_props%cp_015 = inc_props%cp_015
      cc_props%cp_016 = inc_props%cp_016
      cc_props%cp_017 = inc_props%cp_017
      cc_props%cp_018 = inc_props%cp_018
      cc_props%cp_019 = inc_props%cp_019
      cc_props%cp_020 = inc_props%cp_020
      cc_props%cp_021 = inc_props%cp_021
      cc_props%cp_022 = inc_props%cp_022
      cc_props%cp_023 = inc_props%cp_023
      cc_props%cp_024 = inc_props%cp_024
      cc_props%cp_025 = inc_props%cp_025
      cc_props%cp_026 = inc_props%cp_026
      cc_props%cp_027 = inc_props%cp_027
      cc_props%cp_028 = inc_props%cp_028
      cc_props%cp_029 = inc_props%cp_029
      cc_props%cp_030 = inc_props%cp_030
      cc_props%cp_031 = inc_props%cp_031
      cc_props%cp_032 = inc_props%cp_032
      cc_props%cp_033 = inc_props%cp_033
      cc_props%cp_034 = inc_props%cp_034
      cc_props%cp_035 = inc_props%cp_035
      cc_props%cp_036 = inc_props%cp_036
      cc_props%cp_037 = inc_props%cp_037
      cc_props%cp_038 = inc_props%cp_038
      cc_props%cp_039 = inc_props%cp_039
      cc_props%cp_040 = inc_props%cp_040
      cc_props%cp_041 = inc_props%cp_041
      cc_props%cp_042 = inc_props%cp_042
      cc_props%cp_043 = inc_props%cp_043
      cc_props%cp_044 = inc_props%cp_044
      cc_props%cp_045 = inc_props%cp_045
      cc_props%cp_046 = inc_props%cp_046
      cc_props%cp_047 = inc_props%cp_047
      cc_props%cp_048 = inc_props%cp_048
      cc_props%cp_049 = inc_props%cp_049
      cc_props%cp_050 = inc_props%cp_050
      cc_props%cp_051 = inc_props%cp_051
      cc_props%cp_052 = inc_props%cp_052
      cc_props%cp_053 = inc_props%cp_053
      cc_props%cp_054 = inc_props%cp_054
      cc_props%cp_055 = inc_props%cp_055
      cc_props%cp_056 = inc_props%cp_056
      cc_props%cp_057 = inc_props%cp_057
      cc_props%cp_058 = inc_props%cp_058
      cc_props%cp_059 = inc_props%cp_059
      cc_props%cp_060 = inc_props%cp_060
      cc_props%cp_061 = inc_props%cp_061
      cc_props%cp_062 = inc_props%cp_062
      cc_props%cp_063 = inc_props%cp_063
      cc_props%cp_064 = inc_props%cp_064
      cc_props%cp_065 = inc_props%cp_065
      cc_props%cp_066 = inc_props%cp_066
      cc_props%cp_067 = inc_props%cp_067
      cc_props%cp_068 = inc_props%cp_068
      cc_props%cp_069 = inc_props%cp_069
      cc_props%cp_070 = inc_props%cp_070
      cc_props%cp_071 = inc_props%cp_071
      cc_props%cp_072 = inc_props%cp_072
      cc_props%cp_073 = inc_props%cp_073
      cc_props%cp_074 = inc_props%cp_074
      cc_props%cp_075 = inc_props%cp_075
      cc_props%cp_076 = inc_props%cp_076
      cc_props%cp_077 = inc_props%cp_077
      cc_props%cp_078 = inc_props%cp_078
      cc_props%cp_079 = inc_props%cp_079
      cc_props%cp_080 = inc_props%cp_080
      cc_props%cp_081 = inc_props%cp_081
      cc_props%cp_082 = inc_props%cp_082
      cc_props%cp_083 = inc_props%cp_083
      cc_props%cp_084 = inc_props%cp_084
      cc_props%cp_085 = inc_props%cp_085
      cc_props%cp_086 = inc_props%cp_086
      cc_props%cp_087 = inc_props%cp_087
      cc_props%cp_088 = inc_props%cp_088
      cc_props%cp_089 = inc_props%cp_089
      cc_props%cp_090 = inc_props%cp_090
      cc_props%cp_091 = inc_props%cp_091
      cc_props%cp_092 = inc_props%cp_092
      cc_props%cp_093 = inc_props%cp_093
      cc_props%cp_094 = inc_props%cp_094
      cc_props%cp_095 = inc_props%cp_095
      cc_props%cp_096 = inc_props%cp_096
      cc_props%cp_097 = inc_props%cp_097
      cc_props%cp_098 = inc_props%cp_098
      cc_props%cp_099 = inc_props%cp_099
      cc_props%cp_100 = inc_props%cp_100
c
      cc_props%angle_type       = atype
      cc_props%angle_convention = aconv
      cc_props%nslip            = inc_props%nslip
c
      cc_props%h_type           = inc_props%h_type
      cc_props%num_hard         = inc_props%num_hard
      cc_props%tang_calc        = inc_props%tang_calc
      cc_props%debug            = debug
      cc_props%s_type           = inc_props%s_type
      cc_props%cnum             = inc_props%cnum
c
c                    vectors and arrays
c
      cc_props%st_it(1) = inc_props%st_it(1)
      cc_props%st_it(2) = inc_props%st_it(2)
      cc_props%st_it(3) = inc_props%st_it(3)
      call mm10_a_copy_vector( cc_props%g, inc_props%rotation_g, 9 )
      call mm10_a_copy_vector( cc_props%ms, inc_props%ms,
     &                         6*max_slip_sys )
      call mm10_a_copy_vector( cc_props%qs, inc_props%qs,
     &                         3*max_slip_sys )
      call mm10_a_copy_vector( cc_props%ns, inc_props%ns,
     &                         3*max_slip_sys )
      call mm10_a_copy_vector( cc_props%stiffness,
     &                         inc_props%init_elast_stiff, 36 )
c
      return
      end subroutine

c
c *********************************************************************
c *                                                                   *
c *         Initialize USER hardening model                           *
c *                                                                   *
c *********************************************************************
c
      subroutine mm10_init_user( props, tau_tilde, uhist )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      double precision, dimension(props%num_hard) :: tau_tilde
      double precision, dimension(max_uhard) :: uhist
c
      write (props%out,*) "Not implemented"
      call die_gracefully
c
      return
      end
c
      subroutine mm10_setup_user( props, np1, n )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      write(props%out,*) "Not implemented"
      call die_gracefully
c
      return
      end

c *********************************************************************
c *                                                                   *
c *         Initialize simple Voche hardening model                   *
c *                                                                   *
c *********************************************************************
c
      subroutine mm10_init_voche( props, tau_tilde, uhist )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      double precision, dimension(max_uhard) :: tau_tilde
      double precision, dimension(max_uhard) :: uhist
      double precision, parameter :: init_hard = 1.0d-5
c
      tau_tilde(1) = props%tau_y+init_hard
c
      return
      end
c
      subroutine mm10_setup_voche( props, np1, n )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
c              define a mu_harden at state np1 for the
c              sake of the CP model
c
      np1%mu_harden = props%stiffness(6,6)
c
c
c              alternate mode for Voce: use fixed pre-factor
c
      if( props%alter_mode ) np1%dg =  props%eps_dot_0_y * np1%tinc
c
      return
 9000 format(/1x,'>>>>> Warning: routine mm10_init_voche has',
     &  /,1x,    '      coded patch to eliminate \dot gamma^0',
     &  /,1x,    '      change code if desired, comment this',
     &  /,1x,    '      message and run....',
     &  /,1x,    '      terminated.'/)
c
      end

c *********************************************************************
c *                                                                   *
c *         Initialize MTS hardening model                            *
c *                                                                   *
c *********************************************************************
c
      subroutine mm10_init_mts( props, tau_tilde, uhist )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      double precision, dimension(max_uhard) :: tau_tilde
      double precision, dimension(max_uhard) :: uhist
c
      double precision, parameter :: mone = -1.0d0
c
      tau_tilde(1) = mone     ! This only works because these
c                               are actually flags
      uhist(1) = mone
      uhist(2) = mone
c
      return
      end
c
      subroutine mm10_setup_mts( props, np1, n )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      double precision :: dgc
      double precision, parameter :: init_hard = 0.1d0
c
c              continuum effective rate
c
      dgc = np1%dg / np1%tinc
c
c              new shear modulus
c
      if( np1%temp .eq. zero ) then
        np1%mu_harden = props%mu_0
      else
        np1%mu_harden = props%mu_0 - props%D_0 /
     &                  ( exp(props%T_0/np1%temp)- one )
      end if
c
c              threshold contributions
c
      if( dgc .eq. zero) then ! happens during initial stiffness
c                               matrix setup only
         np1%tau_v = props%tau_hat_v
         np1%tau_y = props%tau_hat_y
      else
         np1%tau_v = props%tau_hat_v*(one-(props%boltzman*np1%temp
     &       /(np1%mu_harden*(props%burgers**three)*props%G_0_v)*
     &       dlog(props%eps_dot_0_v/dgc))**(one/props%q_v))
     &       **(one/props%p_v)
         np1%tau_y = props%tau_hat_y*(one-(props%boltzman*np1%temp
     &       /(np1%mu_harden*(props%burgers**three)*props%G_0_y)*
     &       dlog(props%eps_dot_0_y/dgc))**(one/props%q_y))
     &       **(one/props%p_y)
      end if
c
c              use existing labels as a convenience,
c              actually get/set the history
c
      np1%u(1) = np1%tau_y
      if( n%u(1) .lt. zero ) then
        n%tau_y = np1%tau_y
      else
        n%tau_y = n%u(1)
      end if
c
      np1%u(2) = np1%mu_harden
      if( n%u(2) .lt. zero ) then
        n%mu_harden = np1%mu_harden
      else
        n%mu_harden = n%u(2)
      end if
c
c              same here -- check for previous step flag
c
      if( n%tau_tilde(1) .lt. zero ) then
        n%tau_tilde(1) = props%tau_a +
     &     (np1%mu_harden/props%mu_0)*np1%tau_y + init_hard
      end if
c
      return
      end
c
c **********************************************************************
c *                                                                    *
c *         Initialize Ma-Roters-Raabe hardening model                 *
c *                                                                    *
c **********************************************************************
c
c Variable conversion table:
c      WARP3D        Matlab
c      rate_n
c      tau_hat_y     c7
c      G_0_y         v_attack
c      burgers,
c      p_v
c      q_v
c      boltzman
c      theta_0       rho_initial
c      eps_dot_0_v
c      eps_dot_0_y
c      p_y
c      q_y
c      tau_a         Qbulk
c      tau_hat_v     c8
c      G_0_v         Qslip
c      k_0
c      mu_0
c      D_0
c      T_0
c      voche_m
c      u1            c1
c      u2            c2
c      u3            c3
c      u4            c4
c      u5            c5
c      u6            c6
c
      subroutine mm10_init_mrr( props, tau_tilde, uhist )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      double precision :: tau_tilde(props%num_hard)
      double precision, dimension(max_uhard) :: uhist
c
c
      if( props%theta_0 .eq. zero ) then
         tau_tilde(1:props%num_hard) = 1.0d8 ! initial densities for
c                                              edges
      else
         tau_tilde(1:props%num_hard) = props%theta_0 ! initial densities
c                                                      for edges
      end if
c
      return
      end
c
      subroutine mm10_setup_mrr( props, np1, n )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision :: time
c
c              increment the total time
c
      time     = n%u(1) + np1%tinc
      np1%u(1) = time
c
      return
      end
c
c **********************************************************************
c *                                                                    *
c *         Initialize ORNL ferritic-martensitic steel hardening model *
c *                                                                    *
c **********************************************************************
c
c Variable conversion table:
c      WARP3D        Matlab
c      rate_n
c      tau_hat_y     c7
c      G_0_y         v_attack
c      burgers,
c      p_v
c      q_v
c      boltzman
c      theta_0       rho_initial
c      eps_dot_0_v
c      eps_dot_0_y
c      p_y
c      q_y
c      tau_a         Qbulk
c      tau_hat_v     c8
c      G_0_v         Qslip
c      k_0
c      mu_0
c      D_0
c      T_0
c      voche_m
c      u1            c1
c      u2            c2
c      u3            c3
c      u4            c4
c      u5            c5
c      u6            c6
c
c
      subroutine mm10_init_ornl( props, tau_tilde, uhist )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      double precision :: tau_tilde(props%num_hard)
      double precision, dimension(max_uhard) :: uhist
c
c              initial densities for edges
c
      if( props%theta_0 .eq. zero ) then
         tau_tilde(1:props%num_hard) = 1.0d8
      else
         tau_tilde(1:props%num_hard) = props%theta_0
      end if
c
      return
      end subroutine
c
      subroutine mm10_setup_ornl( props, np1, n )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision :: time
c
c              increment the total time
c
      time     = n%u(1) + np1%tinc
      np1%u(1) = time
c
      return
      end
c
c
c **********************************************************************
c *                                                                    *
c *         Initialize  Armstrong-Frederick hardening model            *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_init_arfr( props, tau_tilde, uhist )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      double precision :: tau_tilde(props%num_hard)
      double precision, dimension(max_uhard) :: uhist
c
      if( props%s_type .eq. 9 ) then
                tau_tilde(1:3) = props%G_0_y ! Initial g_0 (MPa)
                tau_tilde(4:6) = props%eps_dot_0_y
      elseif( props%s_type .eq. 10 ) then
                tau_tilde(1:3) = props%G_0_y ! Initial g_0 (MPa)
                tau_tilde(4:6) = props%eps_dot_0_y
                tau_tilde(7:18) = props%G_0_v
      else
          write(props%out,101) props%s_type
          call die_gracefully
      end if
c
      return
c
 101  format(
     &  10x,'>> Error: initial values not defined for Ti6242 for ',
     & 'i6', '.', /,10x, 'Aborting...')
c
      end
c
      subroutine mm10_setup_arfr( props, np1, n )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision :: time
c
c              increment the total time
c
      time     = n%u(1) + np1%tinc
      np1%u(1) = time
c
      return
      end
c
c
c **********************************************************************
c *                                                                    *
c *         Initialize  AFRL Ti-6242 high temperature hardening model  *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_init_DJGM( props, tau_tilde, uhist )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      double precision :: tau_tilde(props%num_hard)
      double precision, dimension(max_uhard) :: uhist
c
      if( props%s_type .eq. 7 ) then !bcc12
                tau_tilde(1:12) = props%G_0_y ! Initial g_0 (MPa)
      elseif( props%s_type .eq. 9 ) then
                tau_tilde(1:3) = props%G_0_y ! Initial g_0 (MPa)
                tau_tilde(4:6) = props%eps_dot_0_y
      elseif( props%s_type .eq. 10 ) then
                tau_tilde(1:3) = props%G_0_y ! Initial g_0 (MPa)
                tau_tilde(4:6) = props%eps_dot_0_y
                tau_tilde(7:18) = props%G_0_v
      else
          write(props%out,101) props%s_type
          call die_gracefully
      end if
c
      return
c
 101  format(
     &  10x,'>> Error: initial values not defined for Ti6242 for ',
     & 'i6', '.', /,10x, 'Aborting...')
c
      end
c
      subroutine mm10_setup_DJGM( props, np1, n )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision :: time
c
c              increment the total time
c
      time     = n%u(1) + np1%tinc
      np1%u(1) = time
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10c_unknown_hard_error          *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10c_unknown_hard_error( props )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
c
      write(props%out,101) props%h_type
 101  format(
     &      10x,'>> Error: unknown hardening type ', 'i6', '.',
     &    /,10x, 'Aborting...')
      call die_gracefully
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_copy_cc_hist                 *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *    Initialize the state n structure                          *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_copy_cc_hist( crys_no, span, history, gradfe,
     &                              R, props, n)
      use mm10_defs
      use mm10_constants
      implicit none
c
      integer :: crys_no, span, sh, eh
c
      double precision :: history(span,*), gradfe(3,3,3), R(3,3)
      type(crystal_props) :: props
      type(crystal_state) :: n
c
      integer :: len1, len2
c!DIR$ ASSUME_ALIGNED history:64, gradfe:64, R:64
c
      len1 = props%num_hard
c
      call mm10_a_copy_vector( n%R, R, 9 )
      call mm10_a_copy_vector( n%gradFeinv, gradfe, 27 )
c
c              scalars only used at n+1
c
      n%temp      = zero
      n%tinc      = zero
      n%dg        = zero
      n%tau_v     = zero
      n%tau_y     = zero
      n%mu_harden = zero
c
      sh = index_crys_hist(crys_no,1,1)
      eh = index_crys_hist(crys_no,1,2)
      n%stress(1:6) = history(1,sh:eh)
c
      sh = index_crys_hist(crys_no,2,1)
      eh = index_crys_hist(crys_no,2,2)
      n%euler_angles(1) = history(1,sh+0)
      n%euler_angles(2) = history(1,sh+1)
      n%euler_angles(3) = history(1,sh+2)

c
      sh = index_crys_hist(crys_no,3,1)
      eh = index_crys_hist(crys_no,3,2)
      n%Rp(1,1) = history(1,sh+0)
      n%Rp(2,1) = history(1,sh+1)
      n%Rp(3,1) = history(1,sh+2)
      n%Rp(1,2) = history(1,sh+3)
      n%Rp(2,2) = history(1,sh+4)
      n%Rp(3,2) = history(1,sh+5)
      n%Rp(1,3) = history(1,sh+6)
      n%Rp(2,3) = history(1,sh+7)
      n%Rp(3,3) = history(1,sh+8)
c
      sh = index_crys_hist(crys_no,4,1)
      eh = index_crys_hist(crys_no,4,2)
      n%D(1:6) = history(1,sh:eh)
c
      sh = index_crys_hist(crys_no,5,1)
      eh = index_crys_hist(crys_no,5,2)
      n%eps(1:6) = history(1,sh:eh)
c
      sh = index_crys_hist(crys_no,6,1)
      eh = index_crys_hist(crys_no,6,2)
      len2 = eh - sh + 1
      if( len2 .ne. length_crys_hist(6) ) then ! sanity check
         write(props%out,9000) 1
         call die_abort
      end if
      n%slip_incs(1:len2-1) = history(1,sh:eh)
c
      sh = index_crys_hist(crys_no,7,1)
      eh = index_crys_hist(crys_no,7,2)
      n%tau_tilde(1:len1) = history(1,sh:sh-1+len1)
c
      sh = index_crys_hist(crys_no,8,1)
      eh = index_crys_hist(crys_no,8,2)
      len2 = eh - sh + 1
      if( len2 .ne. length_crys_hist(8) ) then ! sanity check
         write(props%out,9000) 2
         call die_abort
      end if
      n%u(1:len2-1) = history(1,sh:eh)
c
      sh = index_crys_hist(crys_no,9,1)
      eh = index_crys_hist(crys_no,9,2)
      n%tt_rate(1:len1) = history(1,sh:sh-1+len1)
c
      sh = index_crys_hist(crys_no,10,1)
      eh = index_crys_hist(crys_no,10,2)
      n%ep(1:6) = history(1,sh:eh)
c
      sh = index_crys_hist(crys_no,11,1)
      eh = index_crys_hist(crys_no,11,2)
      n%ed(1:6) = history(1,sh:eh)
c
      return
c
 9000 format('>> FATAL ERROR: mm10_copy_cc_hist @ ',i2 )
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_init_general_hist            *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 9/30/2016 rhd               *
c     *                                                              *
c     *    initialize general GP history (grad Fe, tangent, R)       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_init_general_hist( span, history )
      use mm10_defs, only : indexes_common
      use mm10_constants
      implicit none
c
      integer :: span
      double precision :: history(span,*)
c
      integer :: sh, eh
c!DIR$ ASSUME_ALIGNED history:64
c
c              cep
c
      sh  = indexes_common(1,1)
      eh  = indexes_common(1,2)
      history(1,sh:eh) = zero
c
c              grad_fe
c
      sh  = indexes_common(2,1)
      eh  = indexes_common(2,2)
      history(1,sh:eh) = zero
c
c              R from F=R*U
c
      sh  = indexes_common(3,1)
      history(1,sh+0) = one
      history(1,sh+1) = zero
      history(1,sh+2) = zero
      history(1,sh+3) = zero
      history(1,sh+4) = one
      history(1,sh+5) = zero
      history(1,sh+6) = zero
      history(1,sh+7) = zero
      history(1,sh+8) = one
c
      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_solve_strup                  *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/5/2016 rhd               *
c     *                                                              *
c     *     Solve the stress update adaptively (if required)         *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_solve_strup( props, np1, n, vec1, vec2, arr1,
     &                             arr2, ivec1, ivec2, fail, gp,
     &                             iter_0_extrapolate_off,
     &                             no_load, Jmat, nrJmat, ncJmat )
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
c              parameters  -- mm10_defs has param_def -> maxuhard)
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      integer :: gp, nrJmat, ncJmat
      logical :: fail, iter_0_extrapolate_off, no_load
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
      complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
      double precision :: Jmat(nrJmat,ncJmat)
c
c              locals
c
      type(crystal_state) :: curr
      integer :: cuts, i, gpp,  faili(10), len1
      double precision :: stress(6), ostress(6), R1(6), frac, step,
     &                    failr(10), temp1, temp2
      logical :: debug, gpall, locdebug
c
c              automatic vectors
c
      double precision :: tt(props%num_hard), ott(props%num_hard)
c
      double precision, parameter :: mult = 0.5d0
      integer, parameter :: mcuts = 4
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, arr1:64, arr2:64, Jmat:64
c
      debug    = props%debug
      gpall    = props%gpall ! true to print iteration norms all GPs
      gpp      = props%gpp   ! set one particular G.P. to print
      locdebug = .false.
c
      frac = zero
      step = one
      cuts = 0
c
      len1         = props%num_hard
      stress(1:6)  = n%stress(1:6)
      ostress(1:6) = stress(1:6)
      call mm10_a_copy_vector( tt, n%tau_tilde, len1 )
      call mm10_a_copy_vector( ott, tt, len1 )
      if( locdebug ) write(props%out,*) " n%tt_rate =", n%tt_rate(1)
c
      call mm10_setup( props, np1, n )
      fail  = .false.
      temp1 = np1%D(1)**2 + np1%D(2)**2 + np1%D(3)**2 + np1%D(4)**2 +
     &        np1%D(5)**2 + np1%D(6)**2
c
      temp2 = stress(1)**2 + stress(2)**2 +
     &        stress(3)**2 + stress(4)**2 +
     &        stress(5)**2 + stress(6)**2
      no_load = (temp2 .eq. zero) .and. (temp1 .eq. zero)
c
c              update for stress: either elasticity, or plasticity
c
      if( iter_0_extrapolate_off .or. no_load ) then
c
c              assign elastic stiffness as the tangent matrix.
c              check whether to update the stress accounting for
c              evolving creep according to
c                 np1%stress = n%stress - R1(n%stress,uddt)
c              where
c                 R1 = stress - n%stress - matmul( props%stiffness,
c                      np1%D - dbarp ) + 2.0d0 * symTW
c
        call mm10_solve_strup_a
c
      else
c
c              update stress, state variables using usual material
c              N-R solvers
c
        call mm10_solve_strup_iterate
        if( fail ) return
c
      end if
c
      np1%stress(1:6) = stress(1:6)
      call mm10_a_copy_vector( np1%tau_tilde, tt, len1 )
      call mm10_a_copy_vector( np1%tt_rate, curr%tt_rate, max_uhard )
c
      if( locdebug ) then
           write(props%out,*) " stress(2)=", stress(2),
     &                        " np1%stress=", np1%stress(2)
           write(props%out,*) "  tt=", tt(1), " np1%tau_tilde =",
     &                        np1%tau_tilde(1)
           write(props%out,*) " np1%tt_rate =", np1%tt_rate(1)
      end if
c
      return

      contains  ! for mm10_solve_strup
c
c              ****************************************************
c              *  contains: mm10_solve_strup_a                    *
c              ****************************************************
c
      subroutine mm10_solve_strup_a
      implicit none
c
      call mm10_a_copy_vector( np1%tangent, props%stiffness, 36 )
      if( .not. no_load ) then
          call mm10_formvecs( props, np1, n, stress, tt, vec1, vec2 )
          call mm10_formR1( props, np1, n, vec1, vec2, stress, tt,
     &                      R1, gp)
          stress(1:6) = stress(1:6) - R1(1:6)
      end if
c
c              ensure that zero rates are set, though this value
c              should not be kept in warp3d history during
c              extrapolaion anyway
        curr%tt_rate = zero
c
      return
c
      end subroutine mm10_solve_strup_a
c
c              ****************************************************
c              *  contains: mm10_solve_strup_iterate              *
c              ****************************************************
c
      subroutine mm10_solve_strup_iterate
      implicit none
c
      logical :: solver_failed
      double precision :: D_work(6), tinc_work, temp_work
c
      do while( frac .lt. one )
c
        D_work(1:6) =  np1%D*(step+frac)
        tinc_work = np1%tinc*(step+frac)
        temp_work = (np1%temp-n%temp)*(step+frac)+n%temp
        call mm10_setup_np1( np1%R, D_work, tinc_work, temp_work,
     &                   np1%step, np1%elem, np1%iter, np1%gp, curr)
c
        call mm10_setup( props, curr, n )
c
c              some subroutines modify n%tau_tilde in setup
c
        call mm10_a_copy_vector(  tt, n%tau_tilde, len1 )
c
        if( props%solver ) then
          call mm10_solve( props, curr, n, vec1, vec2, arr1, arr2,
     &         ivec1, ivec2, stress, tt, fail, faili, failr, gp,
     &         np1%tinc*step, Jmat, nrJmat, ncJmat )
        else
          call mm10_solveB( props, curr, n, vec1, vec2, arr1, arr2,
     &             ivec1, ivec2, stress, tt, fail, faili, failr, gp,
     &             np1%tinc*step, Jmat )
        end if
        if( fail ) then
          if( locdebug ) write(props%out,*) "Adapting"
          stress(1:6) = ostress(1:6)
          call mm10_a_copy_vector( tt, ott, len1 )
          step = step * mult
          cuts = cuts + 1
          if( cuts .gt. mcuts ) exit
          fail = .false.
        else
          ostress(1:6) = stress(1:6)
          call mm10_a_copy_vector( ott, tt, len1 )
          frac = frac + step
        end if
c
      end do ! over iterations
c
c              solver succeeded ?
c
      solver_failed = fail .or. any(isnan(tt)) .or. any(isnan(stress))
      if( .not. solver_failed ) return ! good to go
c
c              solver failed. messages. new stresses = old stresses
c
      write(props%out,*)" >>> Warning: mm10 implicit solution failed."
      if( faili(4) .eq. 1 ) then
          write(props%out,*)" Stress prediction failed at iter=",
     &                        faili(1), " for miter=", faili(2)
      else
          write(props%out,*)" Material update failed at iter=",
     &                        faili(1), " for miter=", faili(2)
      end if
      if( faili(3) .eq.1 ) then
          write(props%out,*)" Reason: absolute residual norm"
      else if( faili(3) .eq. 2 ) then
          write(props%out,*)" Reason: relative residual norm"
      else
          write(props%out,*)" Reason: encountered NaN"
      end if
      write(props%out,*) vec1(1:6)
      write(props%out,*) tt(1:6)
      write(props%out,*) stress, np1%D
c      write(props%out,*) props%Hmat(1:7,1)
      write(props%out,*) np1%ms(1:6,1)
      write(props%out,*)" AbsNorm=", failr(1), " AbsTol=", failr(2)
      write(props%out,*)" RelNorm=", failr(3), " RelTol=", failr(4)
      write(props%out,*)" Error occured in: element=",
     &  np1%elem, " GaussPoint=", np1%gp
      write(props%out,*)" Fractional step frac=",
     &  frac, " step=", step, " cuts=", cuts
c
      fail = .true.
      np1%stress(1:6) = n%stress(1:6)
      call mm10_a_copy_vector( np1%tau_tilde, n%tau_tilde, len1 )
      return

      end subroutine mm10_solve_strup_iterate
c
      end subroutine mm10_solve_strup
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_solve                        *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/6/2016 rhd               *
c     *                                                              *
c     *     Solve a stress update to prescribed strain state         *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_solve( props, np1, n, vec1, vec2, arr1, arr2,
     &                       ivec1, ivec2, stress, tt, fail, faili,
     &                       failr, gp, dtinc, J, nrJ, ncJ )
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
c                 parameter definitions
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      integer :: gp, faili(10), nrJ, ncJ
      double precision :: dtinc
      double precision :: stress(6), tt(props%num_hard), failr(10)
      double precision, dimension(max_uhard) :: vec1,vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
      complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
      double precision :: J(nrJ,ncJ)
c
      logical :: fail
c
c                 locals
c
      integer :: iter, miter, info, ls, gpp, ttind, maxfev,
     &           nfev, njev, i
      double precision :: nR, inR, atol, rtol, uB, alpha, ls1, ls2,
     &                    nlsx, nRs, dt, cos_ang, xetol, xtol1,
     &                    zerotol, dxerr, nR1, inR1, atol1, rtol1,
     &                    factor, mm10_enorm, t1, t2, dtrace,
     &                    stepmx
c
      double precision, dimension(6,6) :: J11
c
      double precision, dimension(6) :: R1, x1, dx1, xnew1, d1, d2
      logical :: debug, gpall, solver, strategy, locdebug
c
c                 automatic vectors-arrays
c
      double precision, dimension(6+props%num_hard) :: R, x, dx,
     &                xnew, g, work_vec1
      integer, dimension(6+props%num_hard) :: ipiv
      double precision, dimension(props%num_hard) :: x2
      double precision :: trans_J(ncJ,nrJ), minus_J(nrJ,ncJ)
c
c                 constants
c
      double precision, parameter :: c = 1.0d-4, red = 0.5d0,
     &                               xtol = 0.001d0
      integer, parameter :: mls = 10, mmin = 1
c
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, arr1:64, arr2:64
c!DIR$ ASSUME_ALIGNED ivec1:64, ivec2:64, stress:64, tt:64, J:64
c!DIR$ ASSUME_ALIGNED trans_J:64, minus_J:64
c
      locdebug = .false.
      if( locdebug ) write(*,*) "Entering mm10_solve"
c
      atol    = props%atol
      atol1   = props%atol1
      rtol    = props%rtol
      rtol1   = props%rtol1
      xetol   = props%xtol
      xtol1   = props%xtol1
      miter   = props%miter
      ttind   = 1 ! index of tt to print while debugging
c
c              trust region parameters
c
      maxfev = 3*miter
c
c              debug flags for printing iteration behavior
c
      debug = props%debug
      gpall = props%gpall ! true to print iteration norms all GPs
      gpp   = props%gpp   ! set one particular G.P. to print
c
c              solver flags
c
      solver   = props%solver ! true for Mark's N.R. routine,
c                               false for trust-region
      strategy = props%strategy ! true for geometric l.s.,
c                                 false for cubic l.s.
c
      x(1:6) = stress(1:6)
      x(7:props%num_hard+6) = tt(1:props%num_hard)
c
      if( locdebug ) write(props%out,*)" Guess syy=", x(2)
c
c              prediction of yield stress to initialize the
c              integration algorithm; helps for Orowan flow rule
c              type models
c
      if( props%h_type .gt. 0 ) then
          call mm10_solve_predict
      else
          inR1 = one
      end if
c
c              material update algorithm
c
      call mm10_solve_update
c
      stress(1:6) = x(1:6)
      tt(1:props%num_hard) = x(7:6+props%num_hard)
c
      return
c
      contains
c
c              ****************************************************
c              *  contains: mm10_solve_predict                    *
c              ****************************************************

      subroutine mm10_solve_predict
      implicit none
c
      integer :: iout
      logical :: tan_mat_implemented
      double precision :: dotR1, work_vec(6), minus_J11(6,6)
c
      iout = props%out
      if( debug ) write(iout,*) "Entering stress prediction"
      iter = 0
      x1   = x(1:6)
c
c              extrapolate hardening variables.
c              predict the new hardening variable by extrapolation.
c              use cosine of the "angle" between the new and old
c                displacement increment to indicate continuity
c                of load direction.
c
      dt      = dtinc !np1%tinc
      d1(1:6) = np1%D(1:6)
      dtrace  = (d1(1) + d1(2) + d1(3))/three
      d1(1)   = d1(1) - dtrace
      d1(2)   = d1(2) - dtrace
      d1(3)   = d1(3) - dtrace
      t1      = d1(1)**2 + d1(2)**2 + d1(3)**2
      t2      = d1(4)**2 + d1(5)**2 + d1(6)**2
      if( t1+t2 .eq. zero ) then
          d1(1:6) = zero
        else
          d1(1:6) = d1(1:6)/dsqrt( t1+t2 )
      end if
c
      d2(1:6) = n%D(1:6)
      dtrace  = (d2(1) + d2(2) + d2(3))/three
      d2(1)   = d2(1) - dtrace
      d2(2)   = d2(2) - dtrace
      d2(3)   = d2(3) - dtrace
      t1      = d2(1)**2 + d2(2)**2 + d2(3)**2
      t2      = d2(4)**2 + d2(5)**2 + d2(6)**2
      if( t1+t2 .eq. zero ) then
          d2(1:6) = zero
        else
          d2(1:6) = d2(1:6)/sqrt( t1+t2 )
      end if
      t1      = d1(1)*d2(1) + d1(2)*d2(2) + d1(3)*d2(3)
      t2      = d1(4)*d2(4) + d1(5)*d2(5) + d1(6)*d2(6)
      cos_ang = dmax1( t1+t2, zero )
      x2 = x(7:6+props%num_hard) +
     &     cos_ang * n%tt_rate(1:props%num_hard) * dt
c
      if( locdebug ) then
          write(iout,*)" Stress prediction module, G.P.=", gp
          write(iout,*)" Extrapol tt6=",
     &                x2(6), " Previous tt6=", x(6+6)
      end if
c
      call mm10_formvecs( props, np1, n, x1, x2, vec1, vec2 )
      call mm10_formR1( props, np1, n, vec1, vec2, x1,x2, R1, gp )
      nR1 = dsqrt( R1(1)**2 + R1(2)**2 + R1(3)**2 + R1(4)**2 +
     &             R1(5)**2 + R1(6)**2 )
      inR1 = nR1
c
c              run the Newton-Raphson loop
c
      if( locdebug ) write(iout,'("Iter ",i3," norm ",E10.3)')
     &                iter, nR1
      tan_mat_implemented = (props%tang_calc .eq. 0) .or.
     &                      (props%tang_calc .eq. 3) .or.
     &                      (props%tang_calc .eq. 4)
c
      do while( (nR1 .gt. atol1) .and. (nR1/inR1 .gt. rtol1)  )
c
c              Jacobian
c
        if( tan_mat_implemented ) then
           call mm10_formarrs( props, np1, n, x1, x2, vec1, vec2,
     &                         arr1, arr2, 1 )
           call mm10_formJ11( props, np1, n, vec1, vec2, arr1, arr2,
     &                        x1, x2, J11 )
        else if( props%tang_calc.eq. 2 ) then ! complex difference
           call mm10_formJ11i( props, np1, n, ivec1, ivec2, x1,
     &                         x2, J11 )
        else
          write(iout,*) 'real variable finite difference',
     &                ' not available in mm10_solve'
          call die_gracefully
        end if
c
        dx1(1:6) = R1(1:6)
        minus_J11 = -J11
        call DGESV( 6, 1, minus_J11, 6, ipiv, dx1, 6, info )
        if( locdebug ) write(iout,*)" Iter=",
     &                 iter, " dx1=", dx1(2)
c
c              line search
c
        alpha = one
        dotR1 =( R1(1)**2 + R1(2)**2 + R1(3)**2 + R1(4)**2 +
     &           R1(5)**2 + R1(6)**2 )
        ls1   = half *dotR1
        call mm10_a_mult_type_2t( work_vec, J11, R1 ) ! trans(J11)xR1
        ls2 = c * ( dx1(1)*work_vec(1) + dx1(2)*work_vec(2) +
     &              dx1(3)*work_vec(3) + dx1(4)*work_vec(4) +
     &              dx1(5)*work_vec(5) + dx1(6)*work_vec(6) )
        ls    = 0
c
        do  ! LS loop
          nlsx  = ls1 + ls2*alpha
          xnew1 = x1 + alpha*dx1
c
c                   residual
c
          call mm10_formvecs( props, np1, n, xnew1, x2, vec1, vec2 )
          call mm10_formR1( props, np1, n, vec1, vec2,
     &                      xnew1, x2, R1, gp )
          if( locdebug ) write(iout,*) "R1 ls", R1(2)
          dotR1 =( R1(1)**2 + R1(2)**2 + R1(3)**2 + R1(4)**2 +
     &             R1(5)**2 + R1(6)**2 )
          nR1 = dsqrt( dotR1 )
          nRs = half*dotR1
c
c                   ls convergence test
c
          if( (nRs .le. nlsx) .or. (ls .gt. mls) ) then
            x1 = xnew1
            if( locdebug )  write(iout,*)" Iter=",
     &         iter, " syy=", x1(2), " AbsNorm=", nR1, " alpha=", alpha
            exit
          else
            alpha = red*alpha
            ls = ls + 1
          end if
        end do ! line search loop
c
c              next NR iteration
c
        iter = iter + 1
        if( locdebug ) write(iout,
     &           '("Iter ",i3," norm ",E10.3," ls",F10.3," dx",F10.3)')
     &            iter, nR1, alpha, dxerr
c
c              check for failure and record reason
c
        if( (iter .gt. miter) .or. any(isnan(x1)) ) then
          fail     = .true.
          faili(1) = iter
          faili(2) = miter
          if( nR1 .gt. atol1) then
            faili(3) = 1
          else if( nR1/inR1 .gt. rtol1 ) then
            faili(3) = 2
          else if( any(isnan(x1)) ) then
            faili(3) = 3
          end if
          faili(4) = 1
          failr(1) = nR1
          failr(2) = atol1
          failr(3) = nR1/inR1
          failr(4) = rtol1
          return
        end if
c
      end do ! NR loop
c
c              optional statistics from N-R algorithm
c
      if( locdebug ) then
          write(iout,*)" Stress pred conv iter=", iter
          write(iout,*)" AbsNorm=", nR1, " AbsTol=", atol1
          write(iout,*)" RelNorm=", nR1/inR1, " RelTol=", rtol1
          write(iout,*)" Guess syy=",x(2), " actual syy=", x1(2)
      end if
c
c              put predicted stress and hardening back to primary
c              variable x
c
      x(1:6)                = x1(1:6)
      x(7:6+props%num_hard) = x2(1:props%num_hard)
c
      return

      end subroutine mm10_solve_predict
c
c              ****************************************************
c              *  contains: mm10_solve_update                     *
c              ****************************************************
c
      subroutine mm10_solve_update
      implicit none
c
      integer :: iout, isize
      double precision :: dotR
      logical :: tan_mat_implemented
c
      iout  = props%out
      isize = 6 + props%num_hard
c
c              material update algorithm
c
      if( locdebug ) then
         write(iout,*)" Material update module, G.P.=", gp
         write(iout,*)" Guess syy=",x(2), " Guess tt=", x(6+ttind)
      end if
c
      iter = 0
      call mm10_formR( props, np1, n, vec1, vec2, x(1), x(7), R, gp )
      nR  = dsqrt( dot_product(R,R) )   ! R not 6x1
      inR = nR
      if( inR. eq. zero ) inR = inR1
c
c              NR loop
c
      if( locdebug ) write(iout,'("Iter ",i3," norm ",E10.3)') iter, nR
      tan_mat_implemented = (props%tang_calc .eq. 0) .or.
     &                      (props%tang_calc .eq. 3) .or.
     &                      (props%tang_calc .eq. 4)

c
      do while( ((nR > atol) .and. (nR/inR > rtol) ) .or.(iter < mmin) )
c
c              Jacobian -- not 6x6
c
        if( tan_mat_implemented ) then
            call mm10_formJ( props, np1, n, vec1, vec2, arr1, arr2,
     &                       x(1),  x(7), J )

        else if( props%tang_calc.eq. 2 ) then ! complex difference
            call mm10_formJi( props, np1, n, ivec1, ivec2, x(1),
     &                        x(7), J)
        else
          write(iout,*) 'mm10_solve_update real variable finite ',
     &             'difference not available in mm10_solve'
          call die_gracefully
        end if
c
        dx      = R
        minus_J = -J  ! eliminate repeated temporaries
        call DGESV( isize, 1, minus_J, isize, ipiv, dx, isize, info )
c
        if( locdebug) write(iout,*) " Iter=", iter, " dx1=", dx(2),
     &                           " dx2=", dx(6+ttind)
c
c              line search
c
       alpha     = one
       dotR      = dot_product( R, R )  ! R is not 6x1
       ls1       = half * dotR
       trans_J   = transpose( J )  !   J is not 6x6
       work_vec1 = matmul( trans_J, R ) ! eliminates repeated temps
       ls2       = c * dot_product( dx, work_vec1 )
c
       ls = 0
       do
          nlsx = ls1 + ls2*alpha
          xnew = x + alpha*dx
c
c              residual
c
          call mm10_formR( props, np1, n, vec1, vec2, xnew(1),
     &                     xnew(7), R, gp)
          dotR = dot_product( R, R ) ! not 6x 1
          nR   = dsqrt(dotR)
          nRs  = half * dotR
c
c              LS convergence test
c
          if( (nRs .le. nlsx) .or. (ls .gt. mls) ) then
            x = xnew
            if( locdebug ) then
               write(iout,*) " L.S. converged, Iter=", iter, " syy=",
     &              x(2)," tt(5)=", x(6+ttind), " nR=", nR
              write(iout,*) " alpha=", alpha
            end if
            exit ! done with LS iterations
          else
            alpha = red*alpha
            ls = ls + 1
          end if
        end do  ! on line search
c
c              increment and check for failure.
c              record data and reason for failure
c
        iter = iter + 1
        if( locdebug ) write(iout,
     &     '("Iter ",i3," norm ",E10.3," ls",F10.3)') iter, nR, alpha
        if( (iter .gt. miter) .or. any(isnan(x)) ) then
          fail     = .true.
          faili(1) = iter
          faili(2) = miter
          if( nR .gt. atol) then
            faili(3) = 1
          else if( nR/inR .gt. rtol ) then
            faili(3) = 2
          else if( any(isnan(x)) ) then
            faili(3) = 3
          endif
          faili(4) = 2
          failr(1) = nR
          failr(2) = atol
          failr(3) = nR/inR
          failr(4) = rtol
          return
        end if

      end do  ! N-R loop
c
c              optional statistics from N-R algorithm
c
      if( locdebug ) then
        write(iout,*)" Material upd conv, iter=", iter
        write(iout,*)" AbsNorm=", nR, " AbsTol=", atol
        write(iout,*)" RelNorm=", nR/inR, " RelTol=", rtol
        write(iout,*)" Guess syy=",stress(2), " actual syy=", x(2)
        write(iout,*)" Guess tt6=",x2(6), " actual tt6=",x(6+ttind)
      end if

      return
      end subroutine mm10_solve_update

      end subroutine mm10_solve

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_update_rotation              *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 1/14/14                     *
c     *                                                              *
c     *               Update the plastic rotation                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_update_rotation( props, np1, n, vec1, vec2 )
      use mm10_defs
      use mm10_constants
c
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      double precision :: wbarp(3), wbarp_full(3,3), expw(3,3),
     &                    vec1(max_uhard), vec2(max_uhard), w(3,3)
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64
c
      call mm10_form_wbarp( props, np1, n, vec1, vec2, np1%stress,
     &                      np1%tau_tilde, wbarp )
      call mm10_WV2WT( wbarp, wbarp_full )
      call mm10_expmw3x3( wbarp_full, expw )
c
      call mm10_a_mult_type_1( np1%Rp, expw, n%Rp ) ! A = B *C
c
      return
      end

c
c **********************************************************************
c *                                                                    *
c *                     mm10_wv2wt                                     *
c *                                                                    *
c *         written by : mcm                                           *
c *         last modified : 10/1/2016 rhd                              *
c *                                                                    *
c *         skew vector to skew tensor                                 *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_wv2wt( wv, wt )
      use mm10_constants
      implicit none
c
      double precision, dimension(3,3), intent(out) :: wt
      double precision, dimension(3), intent(in) :: wv
c
c
c!DIR$ ASSUME_ALIGNED wv:64, wt:64
c
      wt(1,1) = zero
      wt(1,2) = wv(3)
      wt(1,3) = wv(2)
      wt(2,1) = -wv(3)
      wt(2,2) = zero
      wt(2,3) = wv(1)
      wt(3,1) = -wv(2)
      wt(3,2) = -wv(1)
      wt(3,3) = zero
c
      return
      end
c
c
c **********************************************************************
c *                                                                    *
c *                                mm10_expmw3x3                       *
c *                                                                    *
c *         written by : mcm                                           *
c *         last modified : 3/22/12 mcm                                *
c *                                                                    *
c *         calculates exp(w) where w is a 3x3 skew matrix.            *
c *         returns full 3x3 because the result                        *
c *         is only orthogonal (not skew or symmetric)                 *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_expmw3x3( w, a )
      use mm10_constants
      implicit none
c
      double precision, dimension(3,3), intent(in) :: w
      double precision, dimension(3,3), intent(out) :: a
      double precision :: alpha
      integer :: i
c!DIR$ ASSUME_ALIGNED w:64, a:64
c
c              compute alpha
c
      alpha = dsqrt( w(2,3)**2 + w(1,3)**2 + w(1,2)**2 )
c
c              algorithm will fail with alpha = 0 (=> w=0)
c              correct result is expm(0) = i, so program that in
c
      if( alpha .lt. 1.0d-16 ) then
        a = zero
      else
        a = w
        call dgemm( 'N', 'N', 3, 3, 3, (one-dcos(alpha))/alpha**2,
     &              w, 3, w, 3, sin(alpha)/alpha, a, 3 )
      end if
c
c              add the identity
c
      a(1,1) = a(1,1) + one
      a(2,2) = a(2,2) + one
      a(3,3) = a(3,3) + one
c
      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_output                       *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/22/2016 rhd              *
c     *                                                              *
c     *     Calculate various other user output quantities           *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_output( props, np1, n, vec1, vec2,
     &                        p_strain_ten_c )
      use mm10_defs
      use mm10_constants
      implicit none
c
c              parameters
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: p_strain_ten_c
      double precision, dimension(max_uhard) :: vec1, vec2
c
c              locals
c
      integer :: i, info, sysID, numAct, nslip
      double precision, dimension(6) :: ee, dbarp, ewwe, ep,
     &                                  eeunrot, se, ed, work_vec1
      double precision, dimension(3) :: wp
      double precision, dimension(6,6) :: S, erot
      double precision :: maxslip, ec_dot, n_eff, rs, ec_slip, ep_dot,
     &                    mm10_rs,curslip, B_eff, s_trace, ed_dot
      double precision, dimension(max_slip_sys) :: dgammadtau, dif_slp
      double precision ::  t1, t2
c
c ***** START: Add new Constitutive Models into this block *****
      nslip = props%nslip
      select case( props%h_type )
        case( 1 ) ! Voche
          do i = 1, nslip
           call mm10_slipinc( props, np1, n, np1%stress,
     &                        np1%tau_tilde, i, np1%slip_incs(i) )
c                    addition for diffusion
c           rs = mm10_rs( props, np1, n, np1%stress, np1%tau_tilde, i )
c           np1%slip_incs(i) = np1%slip_incs(i) +
c     &                        rs*np1%tinc*props%iD_v
          end do
       case( 2 ) ! MTS
         do i = 1, nslip
           call mm10_slipinc( props, np1, n, np1%stress,
     &                        np1%tau_tilde, i, np1%slip_incs(i) )
         end do
       case(3)  ! User
         do i =1, nslip
           call mm10_slipinc_user( props, np1, n, np1%stress,
     &                        np1%tau_tilde, i, np1%slip_incs(i) )
         end do
       case( 4 )  ! ORNL
         np1%slip_incs(1:nslip)  = vec1(1:nslip)
       case( 7 ) ! MRR
         np1%slip_incs(1:nslip)  = vec1(1:nslip)
       case( 8 ) ! Armstrong-Frederick
         np1%slip_incs(1:nslip)  = vec1(1:nslip)
       case( 9 )  ! DJGM
         np1%slip_incs(1:nslip)  = vec1(1:nslip)
       case default
         call mm10_unknown_hard_error( props )
      end select
c
      call mm10_update_euler_angles( props, np1, n )
c
c     Compute diffusion rate
      dif_slp = zero
      select case( props%h_type )
        case( 1 ) ! Voche
          do i = 1, nslip
           rs = mm10_rs( props, np1, n, np1%stress, np1%tau_tilde, i )
           dif_slp(i) = rs*np1%tinc*props%iD_v
          end do
          case( 2 )! MTS
          case( 3 )  ! User
          case( 4 ) ! ORNL
          case( 7 ) ! MRR
          case( 8 ) ! Armstrong-Frederick
          case( 9 ) ! DJGM
          case default
            call mm10_unknown_hard_error( props )
      end select
c     Compute manually
          ed = zero
          do i = 1, nslip
            ed = ed + dif_slp(i)*np1%ms(1:6,i)
          end do
c     Isotropic diffusion
        work_vec1 = zero
	if ( abs(props%cp_031-one)<1.0e-5 ) then
	  call mm10_halite_formRpp( props, work_vec1, np1%stress, np1%tinc )
          ed = ed - work_vec1 ! Ran's value is negative
	end if
      np1%ed = ed / np1%tinc
c          store the plastic strain increment for nonlocal averages
      t1 = ed(1)**2 + ed(2)**2 + ed(3)**2
      t2 = ed(4)**2 + ed(5)**2 + ed(6)**2
      ed_dot = dsqrt( two/three*(t1+ half*t2) )
c
c           compute diffusion rate
c
      ed_dot = ed_dot / np1%tinc
      np1%u(15) = ed_dot
c
c              miscellaneous things
c
      np1%work_inc = np1%stress(1)*np1%D(1) + np1%stress(2)*np1%D(2) +
     &               np1%stress(3)*np1%D(3) + np1%stress(4)*np1%D(4) +
     &               np1%stress(5)*np1%D(5) + np1%stress(6)*np1%D(6)
c              plastic strain and work
      call mm10_a_copy_vector( S, props%stiffness, 36 )
      eeunrot(1:6) = np1%stress(1:6)
      call DPOSV( 'U', 6, 1, S, 6, eeunrot, 6, INFO )
      call mm10_RT2RVE( np1%R, erot )
      call mm10_a_mult_type_2( ee, erot, eeunrot )
c
c      call mm10_form_dbarp( props, np1, n,vec1, vec2, np1%stress,
c     &                      np1%tau_tilde, dbarp )
c     Compute manually
          dbarp = zero
          do i = 1, nslip
            dbarp = dbarp + np1%slip_incs(i)*np1%ms(1:6,i)
          end do
c
      call mm10_form_wp( props, np1, n,vec1, vec2, np1%stress,
     &                   np1%tau_tilde, wp )
      call mm10_symSW( ee, wp, ewwe )
c
      ep = dbarp + ewwe
      np1%ep = ep / np1%tinc
c           compute plastic creep rate
      t1 = ep(1)**2 + ep(2)**2 + ep(3)**2
      t2 = ep(4)**2 + ep(5)**2 + ep(6)**2
      ep_dot = dsqrt( two/three*(t1+ half*t2) )
      ep_dot = ep_dot / np1%tinc
      np1%u(11) = ep_dot
c          store the plastic strain increment for nonlocal averages
      ep = ep + ed
      p_strain_ten_c(1:6) = ep(1:6)
      t1 = ep(1)**2 + ep(2)**2 + ep(3)**2
      t2 = ep(4)**2 + ep(5)**2 + ep(6)**2
      np1%p_strain_inc = dsqrt( two/three*(t1+ half*t2) )
c      np1%p_work_inc = dot_product(np1%stress, ep)
      np1%p_work_inc = np1%stress(1)*ep(1) + np1%stress(2)*ep(2) +
     &               np1%stress(3)*ep(3) + np1%stress(4)*ep(4) +
     &               np1%stress(5)*ep(5) + np1%stress(6)*ep(6)
c           store the lattice strains
      np1%eps = ee
c
c           compute effective creep exponent
c
      ec_dot = np1%p_strain_inc / np1%tinc
      if( ec_dot > zero ) then  !   AAA
        ep(1:6) = ep(1:6) / np1%tinc
c
c           generalization of CP model implementation for other
c           slip rate eqns, requiring other forms of d_gamma/d_tau.
c           vector dgammadtau should be 1 x n_slip
c
c ***** START: Add new Constitutive Models into this block *****
        select case( props%h_type  )
          case( 1 ) ! Voche
            call mm10_dgdt_voche( props, np1, n, np1%stress,
     &                            np1%tau_tilde, dgammadtau )
          case( 2 )! MTS
            call mm10_dgdt_mts( props, np1, n, np1%stress,
     &                          np1%tau_tilde, dgammadtau )
          case( 3 )  ! User
            call mm10_dgdt_user( props, np1, n, np1%stress,
     &                           np1%tau_tilde, dgammadtau )
          case( 4 ) ! ORNL
            call mm10_dgdt_ornl( props, np1, n, np1%stress,
     &                           np1%tau_tilde, dgammadtau )
          case( 7 ) ! MRR
            call mm10_dgdt_mrr( props, np1, n, np1%stress,
     &                          np1%tau_tilde, dgammadtau )
          case( 8 ) ! Armstrong-Frederick
            call mm10_dgdt_arfr( props, np1, n, np1%stress,
     &                          np1%tau_tilde, dgammadtau )
          case( 9 ) ! DJGM
            call mm10_dgdt_djgm( props, np1, n, np1%stress,
     &                          np1%tau_tilde, dgammadtau)
          case default
            call mm10_unknown_hard_error( props )
        end select
c ****** END: Add new Constitutive Models into this block ******
        n_eff = zero
        do i = 1, nslip
           rs = mm10_rs( props, np1, n, np1%stress, np1%tau_tilde, i )
           t1 = np1%ms(1,i)*ep(1) + np1%ms(2,i)*ep(2) +
     &          np1%ms(3,i)*ep(3)
           t2 = np1%ms(4,i)*ep(4) + np1%ms(5,i)*ep(5) +
     &          np1%ms(6,i)*ep(6)
           ec_slip = t1 + half*t2
c
c           n_eff = n_eff + two/three/ec_dot/ec_dot*rs*
c     &             dgammadtau(i)/np1%tinc*ec_slip ! from D. Parks
c
c           re-order to suppress floating overflow
c
           t1 = rs * dgammadtau(i) / ec_dot
           t2 = ec_slip / ec_dot / np1%tinc
           n_eff = n_eff + (two/three)*t1*t2
         end do
      else
        n_eff = 1.0d10
      end if  !  AAA
c
      np1%u(12) = n_eff
      s_trace = (np1%stress(1)+np1%stress(2)+np1%stress(3))/three
      se(1:6) = np1%stress(1:6)
      se(1) = se(1) - s_trace
      se(2) = se(2) - s_trace
      se(3) = se(3) - s_trace
      t1    = se(1)**2 + se(2)**2 + se(3)**2
      t2    = se(4)**2 + se(5)**2 + se(6)**2
      s_trace = sqrt( onept5 * ( t1 + two*t2 ) )
      np1%u(13) = s_trace
      if( ec_dot < 1.d-100 ) then ! elastic
        B_eff = zero
      elseif( n_eff > 100.d0 ) then ! likely overflow
        B_eff = -one
      else
        B_eff = ec_dot / s_trace**n_eff
      end if
      np1%u(14) = B_eff
c
c     Combine creep and diffusion slip rates
      np1%slip_incs(1:nslip) = np1%slip_incs(1:nslip) + dif_slp(1:nslip)
c
c              compute user history output
c
c              6:8 is about the active slip systems: maximum slip rate,
c              identifier of that system, and how many systems have rates
c              on the same order of magnitude as that one.
c
      maxslip = zero
      sysID = 0
      do i = 1, nslip
        curslip = dabs(np1%slip_incs(i))
        if( curslip > maxslip ) then
          maxslip = curslip
          sysID = i
        end if
      end do
      np1%u(6) = maxslip/np1%tinc
      np1%u(7) = dble(sysID)
c
      numAct = 0
      maxslip = pt1*maxslip
      do i = 1, nslip
        curslip = dabs(np1%slip_incs(i))
        if( curslip >= maxslip ) numAct = numAct + 1
      end do
      np1%u(8) = dble(numAct)
c
      return
      end






c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_a_mult...                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 9/29/2016 rhd               *
c     *                                                              *
c     *     various matrix multiplication routines for mm10_a        *
c     *     these should be all inlined                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_a_mult_type_1( a, b, c )
      implicit none
      double precision :: a(3,3), b(3,3), c(3,3)
c!DIR$ ASSUME_ALIGNED a:64, b:64, c:64
c
c                     [a] = [b] * [c]
c
      a(1,1) = b(1,1)*c(1,1) + b(1,2)*c(2,1) + b(1,3)*c(3,1)
      a(2,1) = b(2,1)*c(1,1) + b(2,2)*c(2,1) + b(2,3)*c(3,1)
      a(3,1) = b(3,1)*c(1,1) + b(3,2)*c(2,1) + b(3,3)*c(3,1)
c
      a(1,2) = b(1,1)*c(1,2) + b(1,2)*c(2,2) + b(1,3)*c(3,2)
      a(2,2) = b(2,1)*c(1,2) + b(2,2)*c(2,2) + b(2,3)*c(3,2)
      a(3,2) = b(3,1)*c(1,2) + b(3,2)*c(2,2) + b(3,3)*c(3,2)
c
      a(1,3) = b(1,1)*c(1,3) + b(1,2)*c(2,3) + b(1,3)*c(3,3)
      a(2,3) = b(2,1)*c(1,3) + b(2,2)*c(2,3) + b(2,3)*c(3,3)
      a(3,3) = b(3,1)*c(1,3) + b(3,2)*c(2,3) + b(3,3)*c(3,3)
c
      return
      end
c
      subroutine mm10_a_mult_type_1a( a, b )
      implicit none
      double precision :: a(3,3), b(3,3)
c
      double precision :: w(3,3)
c!DIR$ ASSUME_ALIGNED a:64, b:64, w:64
c
c                     [b] = [a] * [b]
c
      w(1,1) = a(1,1)*b(1,1) + a(1,2)*b(2,1) + a(1,3)*b(3,1)
      w(2,1) = a(2,1)*b(1,1) + a(2,2)*b(2,1) + a(2,3)*b(3,1)
      w(3,1) = a(3,1)*b(1,1) + a(3,2)*b(2,1) + a(3,3)*b(3,1)
c
      w(1,2) = a(1,1)*b(1,2) + a(1,2)*b(2,2) + a(1,3)*b(3,2)
      w(2,2) = a(2,1)*b(1,2) + a(2,2)*b(2,2) + a(2,3)*b(3,2)
      w(3,2) = a(3,1)*b(1,2) + a(3,2)*b(2,2) + a(3,3)*b(3,2)
c
      w(1,3) = a(1,1)*b(1,3) + a(1,2)*b(2,3) + a(1,3)*b(3,3)
      w(2,3) = a(2,1)*b(1,3) + a(2,2)*b(2,3) + a(2,3)*b(3,3)
      w(3,3) = a(3,1)*b(1,3) + a(3,2)*b(2,3) + a(3,3)*b(3,3)
c
      b = w   !   a = b * a
c
      return
      end
c
      subroutine mm10_a_mult_type_2( a, b, c )
      implicit none
      double precision :: a(6), b(6,6), c(6)
      integer :: j
c!DIR$ ASSUME_ALIGNED a:64, b:64, c:64
c
c                     a = [b] * c
c
      a(1) = b(1,1)*c(1)
      a(2) = b(2,1)*c(1)
      a(3) = b(3,1)*c(1)
      a(4) = b(4,1)*c(1)
      a(5) = b(5,1)*c(1)
      a(6) = b(6,1)*c(1)

      do j = 2, 6
        a(1) = a(1) + b(1,j)*c(j)
        a(2) = a(2) + b(2,j)*c(j)
        a(3) = a(3) + b(3,j)*c(j)
        a(4) = a(4) + b(4,j)*c(j)
        a(5) = a(5) + b(5,j)*c(j)
        a(6) = a(6) + b(6,j)*c(j)
      end do
c
      return
      end

      subroutine mm10_a_mult_type_4( a, b, c, d, const )
      implicit none
      double precision :: a(6), b(6,6), c(6), d(6), const
      integer :: j
c!DIR$ ASSUME_ALIGNED a:64, b:64, c:64, d:64
c
c                     a = [b] * c + const * d
c
      a(1) = const * d(1)
      a(2) = const * d(2)
      a(3) = const * d(3)
      a(4) = const * d(4)
      a(5) = const * d(5)
      a(6) = const * d(6)

      do j = 1, 6
        a(1) = a(1) + b(1,j)*c(j)
        a(2) = a(2) + b(2,j)*c(j)
        a(3) = a(3) + b(3,j)*c(j)
        a(4) = a(4) + b(4,j)*c(j)
        a(5) = a(5) + b(5,j)*c(j)
        a(6) = a(6) + b(6,j)*c(j)
      end do


c
      return
      end


      subroutine mm10_a_mult_type_2t( a, b, c )
      implicit none
      double precision :: a(6), b(6,6), c(6)
      integer :: j
c!DIR$ ASSUME_ALIGNED a:64, b:64, c:64
c
c                     a = trans( [b] ) * c
c
      a(1) = b(1,1)*c(1) + b(2,1)*c(2) + b(3,1)*c(3) + b(4,1)*c(4) +
     &       b(5,1)*c(5) + b(6,1)*c(6)
c
      a(2) = b(1,2)*c(1) + b(2,2)*c(2) + b(3,2)*c(3) + b(4,2)*c(4) +
     &       b(5,2)*c(5) + b(6,2)*c(6)
c
      a(3) = b(1,3)*c(1) + b(2,3)*c(2) + b(3,3)*c(3) + b(4,3)*c(4) +
     &       b(5,3)*c(5) + b(6,3)*c(6)
c
      a(4) = b(1,4)*c(1) + b(2,4)*c(2) + b(3,4)*c(3) + b(4,4)*c(4) +
     &       b(5,4)*c(5) + b(6,4)*c(6)
c
      a(5) = b(1,5)*c(1) + b(2,5)*c(2) + b(3,5)*c(3) + b(4,5)*c(4) +
     &       b(5,5)*c(5) + b(6,5)*c(6)
c
      a(6) = b(1,6)*c(1) + b(2,6)*c(2) + b(3,6)*c(3) + b(4,6)*c(4) +
     &       b(5,6)*c(5) + b(6,6)*c(6)
c
      return
      end

      subroutine mm10_a_mult_type_3t( a, b, c )
      implicit none
      double precision :: a(3,3), b(3,3), c(3,3)
c!DIR$ ASSUME_ALIGNED a:64, b:64, c:64
c
c                     [a] = [b] * trans( [c] )
c
      a(1,1) = b(1,1)*c(1,1) +  b(1,2)*c(1,2) +  b(1,3)*c(1,3)
      a(1,2) = b(1,1)*c(2,1) +  b(1,2)*c(2,2) +  b(1,3)*c(2,3)
      a(1,3) = b(1,1)*c(3,1) +  b(1,2)*c(3,2) +  b(1,3)*c(3,3)

      a(2,1) = b(2,1)*c(1,1) +  b(2,2)*c(1,2) +  b(2,3)*c(1,3)
      a(2,2) = b(2,1)*c(2,1) +  b(2,2)*c(2,2) +  b(2,3)*c(2,3)
      a(2,3) = b(2,1)*c(3,1) +  b(2,2)*c(3,2) +  b(2,3)*c(3,3)

      a(3,1) = b(3,1)*c(1,1) +  b(3,2)*c(1,2) +  b(3,3)*c(1,3)
      a(3,2) = b(3,1)*c(2,1) +  b(3,2)*c(2,2) +  b(3,3)*c(2,3)
      a(3,3) = b(3,1)*c(3,1) +  b(3,2)*c(3,2) +  b(3,3)*c(3,3)

c
      return
      end

      subroutine mm10_a_mult_type_5( a, b, c )
      implicit none
      double precision :: a(6,6), b(6), c(6)
      integer :: i, j
c!DIR$ ASSUME_ALIGNED a:64, b:64, c:64
c
c                    [a] = [a] + b * trans( c )
c
      do i = 1 , 6
       do j = 1, 6
         a(i,j) = a(i,j) + b(i)*c(j)
       end do
      end do
c
      return
      end


      subroutine mm10_a_mult_type_3( a, b, c )
      implicit none
      double precision :: a(3), b(3,3), c(3)
c!DIR$ ASSUME_ALIGNED a:64, b:64, c:64
c
c                     a = [b] * c
c
      a(1) = b(1,1)*c(1) + b(1,2)*c(2) + b(1,3)*c(3)
      a(2) = b(2,1)*c(1) + b(2,2)*c(2) + b(2,3)*c(3)
      a(3) = b(3,1)*c(1) + b(3,2)*c(2) + b(3,3)*c(3)
c
      return
      end


c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_a_zero_vec                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 9/29/2016 rhd               *
c     *                                                              *
c     *    zero a vector - should be inlined                         *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_a_zero_vec( vec, nterms )
      use mm10_constants
      implicit none
      integer :: nterms
      double precision :: vec(nterms)
c!DIR$ ASSUME_ALIGNED vec:64
c
      vec = zero
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_a_make_symm_1                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 9/29/2016 rhd               *
c     *                                                              *
c     *                 make a symmetric 6x6 matrix                  *
c     *                                                              *
c     ****************************************************************
c

      subroutine mm10_a_make_symm_1( a )
      use mm10_constants
      implicit none
c
      integer :: i, j
      double precision :: a(6,6), work(6,6)
c!DIR$ ASSUME_ALIGNED a:64, work(3,3):64
c
c                     [a] = ( [a] + trans[a] ) / 2.0
c
      do j = 1, 6
        do i = 1, 6
          work(i,j) = ( a(i,j) + a(j,i) ) * half
        end do
      end do
c
      a = work
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_a_copy_vec                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 9/29/2016 rhd               *
c     *                                                              *
c     *     copy a vector - should be inlined                        *
c     *                                                              *
c     ****************************************************************
c

      subroutine mm10_a_copy_vector( a, b, nterms )
      implicit none
      integer :: nterms
      double precision :: a(nterms), b(nterms), c(nterms)
c!DIR$ ASSUME_ALIGNED a:64, b:64, c:64
c
      a = b
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_a_invert_33                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 10/11/2016 rhd              *
c     *                                                              *
c     *                     invert a general 3 x 3 matrix            *
c     *                                                              *
c     ****************************************************************
c


      subroutine mm10_a_invert_33( ainv, a,  ok_flag )
      implicit none
c
      double precision, dimension(3,3), intent(in)  :: a
      double precision, dimension(3,3), intent(out) :: ainv
      logical, intent(out) :: ok_flag
c
      double precision, parameter :: eps = 1.0d-10
      double precision :: det
      double precision, dimension(3,3) :: cofactor
c!DIR$ ASSUME_ALIGNED a:64, ainv:64
c
      det =    a(1,1)*a(2,2)*a(3,3)
     1       - a(1,1)*a(2,3)*a(3,2)
     2       - a(1,2)*a(2,1)*a(3,3)
     3       + a(1,2)*a(2,3)*a(3,1)
     4       + a(1,3)*a(2,1)*a(3,2)
     5       - a(1,3)*a(2,2)*a(3,1)

      if( abs(det) .le. eps ) then
         ainv = 0.0d0
         ok_flag = .false.
         return
      end if
c
      cofactor(1,1) = +(a(2,2)*a(3,3)-a(2,3)*a(3,2))
      cofactor(1,2) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
      cofactor(1,3) = +(a(2,1)*a(3,2)-a(2,2)*a(3,1))
      cofactor(2,1) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
      cofactor(2,2) = +(a(1,1)*a(3,3)-a(1,3)*a(3,1))
      cofactor(2,3) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
      cofactor(3,1) = +(a(1,2)*a(2,3)-a(1,3)*a(2,2))
      cofactor(3,2) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
      cofactor(3,3) = +(a(1,1)*a(2,2)-a(1,2)*a(2,1))
c
c                ainv = transpose(cofactor) / det
c
      ainv(1,1) = cofactor(1,1) / det
      ainv(2,1) = cofactor(1,2) / det
      ainv(3,1) = cofactor(1,3) / det
c
      ainv(1,2) = cofactor(2,1) / det
      ainv(2,2) = cofactor(2,2) / det
      ainv(3,2) = cofactor(2,3) / det
c
      ainv(1,3) = cofactor(3,1) / det
      ainv(2,3) = cofactor(3,2) / det
      ainv(3,3) = cofactor(3,3) / det
c
      ok_flag = .true.
c
      return
c
      end


