c
c **********************************************************************
c *                                                                    *
c *    mm10.f                                                          *
c *                                                                    *
c *         written by : mcm                                           *
c *         last modified : 5/30/2016 rhd                              *
c *                                                                    *
c *         Stress/strain update routines AND HELPERS for              *
c *         crystal plasticity material modeled via Beaudoin et al.    *
c *                                                                    *
c **********************************************************************
c

c
c     ****************************************************************
c     *                                                              *
c     *                       module mm10_defs                       *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 5/22/2016 rhd               *
c     *                                                              *
c     *       small module to hold crystal update data structs       *
c     *       also hold integer indexes into history vector for      *
c     *       an integration point                                   *
c     *                                                              *
c     ****************************************************************
c
      module mm10_defs
c
      implicit integer (a-z)
$add param_def
c              includes all key size limits for WARP3D (e.g. mxvl)
c
c              Massive list of properties
c
      type :: crystal_props
        double precision :: rate_n, tau_hat_y, G_0_y, burgers,
     &        p_v, q_v, boltzman, theta_0, eps_dot_0_v, eps_dot_0_y,
     &        p_y, q_y, tau_a, tau_hat_v, G_0_v,
     &        k_0, mu_0, D_0, T_0, tau_y, tau_v, voche_m,
     &        u1, u2, u3, u4, u5, u6
        double precision :: atol, atol1, rtol, rtol1, xtol, xtol1
        double precision, dimension(3,3) :: g
        double precision :: ms(6,max_slip_sys), qs(3,max_slip_sys),
     &                      ns(3,max_slip_sys)
        double precision, dimension(6,6) :: stiffness
        integer :: angle_type, angle_convention, nslip,
     &        h_type, miter, gpp, s_type, cnum, method,
     &        st_it(3), num_hard, out
        logical :: real_tang, solver, strategy, debug, gpall
        ! constants for use in material models
        double precision, dimension(:,:), allocatable :: Gmat,Hmat
      end type
c
      type :: crystal_state
        double precision, dimension(3,3) :: R, Rp
        double precision, dimension(6) :: stress, D, eps
        double precision, dimension(3) :: euler_angles
        double precision, dimension(max_slip_sys) :: tau_l, slip_incs
        double precision, dimension(3,3,3) :: gradFeinv
        double precision, dimension(6,6) :: tangent
        double precision, dimension(max_uhard) :: tau_tilde
        double precision, dimension(max_uhard) :: tt_rate
        double precision :: temp, tinc, dg, tau_v, tau_y,
     &                      mu_harden, work_inc, p_work_inc,
     &                      p_strain_inc
        double precision :: ms(6,max_slip_sys), qs(3,max_slip_sys),
     &                      qc(3,max_slip_sys)
        double precision, dimension(max_uhard) :: u
        integer :: step, elem, gp, iter
      end type
c
c              store integer indexes into history vector for a an
c              integration point.
c              WARP3D makes sure an mm10 routine is called to
c              create the arraays & set values. No need to save
c              across restarts. These arrays are red-only during
c              threaded processing of element blocks.
c
      integer, allocatable :: indexes_common(:,:),
     &                        index_crys_hist(:,:,:)
c
      end module
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine mm10                         *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 4/6/2016 tjt                *
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
c
      implicit integer (a-z)
$add include_sig_up
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
      integer :: i, c, co, cn, now_element
      double precision, dimension(6) :: sig_avg, se
      double precision, dimension(6) :: p_strain_ten_c, p_strain_ten
      double precision, dimension(6,6) :: tang_avg
      double precision, dimension(max_slip_sys) :: slip_avg
      double precision :: t_work_inc, p_work_inc,p_strain_inc, zero,
     &                    one, two, three, ten, t1, t2
      double precision :: n_avg,p_strain_avg, B_eff, s_trace
      type(crystal_props) :: cc_props
      type(crystal_state) :: cc_n, cc_np1
      logical :: debug, gpall, locdebug
c
      data zero, one, two, three, ten / 0.0d0, 1.0d00, 2.0d00, 3.0d00,
     &    10.0d00 /
c
      debug = .false.
c
      if( debug ) write (iout,*) "In mm10"
c
c                 initialize G,H arrays for certain CP constitutive
c                 models. Note: all crystals in the element block
c                 will have the same hardening model and slip systems
c
      call mm10_set_cons(local_work,cc_props,1)
c
c                 we have data passed for integration point # gp for
c                 all (span) elements in this block. the CP updating
c                 routines process only a single int point. loop to
c                 process that point for all elements in block.
c
c                 for small strain analysis, rot_blk_n1 set to identity
c                 now by warp3d.
c
      local_work%material_cut_step = .false.
c
      do i = 1, span
c
        now_element = local_work%felem + i - 1
c
c                 for step = 1, init element history
c
        if( local_work%step .eq. 1 ) then
              if( debug ) write(iout,*) "Init GP history"
              call mm10_init_general_hist( history_n(i,1:72) )
              call mm10_init_uout_hist( history_n(i,73:75) )
              call mm10_init_slip_hist( history_n(i,76:
     &              76+max_slip_sys-1) )
        end if
c
        if( debug ) write(iout,*) 'Updating element: ', now_element
        sig_avg      = zero
        slip_avg     = zero
        t_work_inc   = zero
        p_work_inc   = zero
        p_strain_inc = zero
        tang_avg     = zero
        p_strain_ten = zero
        n_avg        = zero
c
c                 loop on all crystals at integration point
c
        do c = 1, ncrystals(i)
          debug = local_work%debug_flag(i)
          locdebug = (debug
     &               .and.((local_work%c_props(i,c)%st_it(3)
     &               .eq.-2).or.
     &               (local_work%c_props(i,c)%st_it(3).eq.
     &                i-1+local_work%felem))
     &               .and.((local_work%c_props(i,c)%st_it(1)
     &               .eq.-2).or.
     &               (local_work%c_props(i,c)%st_it(1).eq.
     &               local_work%step))
     &               .and.((local_work%c_props(i,c)%st_it(2)
     &               .eq.-2).or.
     &               (local_work%c_props(i,c)%st_it(2)
     &               .eq.local_work%iter)))
          mat_debug = locdebug
          locdebug = .false.
          co = 76+max_slip_sys+(c-1)*(30+max_slip_sys
     &           +3*max_uhard)
          cn = 76+max_slip_sys+(c)*(30+max_slip_sys
     &           +3*max_uhard)-1
          if( locdebug ) write(iout,*) "Setting up properties"
          call mm10_init_cc_props( local_work%c_props(i,c),
     &              local_work%angle_type(i),
     &              local_work%angle_convention(i),
     &              local_work%debug_flag(i), cc_props )
          cc_props%out = iout
c
          if( local_work%step .eq. 1 ) then
            if( locdebug ) write(iout,*) "Init history 0"
            call mm10_init_cc_hist0( cc_props,
     &           local_work%c_props(i,c)%init_angles(1:3),
     &           history_n(i,co:cn) )
          end if
          if( locdebug ) write(iout,*) "Copying n to struct"
          call mm10_copy_cc_hist( history_n(i,co:cn),
     &         history_n(i,37:63),
     &         history_n(i,64:72),
     &         cc_props, cc_n )
          call mm10_setup_np1(
     &        local_work%rot_blk_n1(i,1:9,gp), uddt(i,1:6),
     &        local_work%dt, gp_temps(i), local_work%step,
     &        i-1+local_work%felem, local_work%iter,
     &        local_work%gpn,cc_np1 )
          if( locdebug) write(iout,*) "Updating crystal ", c
          call mm10_solve_crystal( cc_props, cc_np1, cc_n,
     &        local_work%material_cut_step, iout, .false., 0,
     &        p_strain_ten_c, iter_0_extrapolate_off )
          if( local_work%material_cut_step ) then
            call mm10_set_cons( local_work, cc_props, 2)
            return
          endif
c
c                 accumulate sums for subsequent averaging
c                    cp_strain_inc -> effective plastic increment
c                    p_strain_ten -> plastic strain increment tensor
c                    n_avg -> effective creep exponent
c
          sig_avg      = sig_avg + cc_np1%stress
          tang_avg     = tang_avg + cc_np1%tangent
          slip_avg     = slip_avg + cc_np1%slip_incs
          t_work_inc   = t_work_inc + cc_np1%work_inc
          p_work_inc   = p_work_inc + cc_np1%p_work_inc
          p_strain_inc = p_strain_inc + cc_np1%p_strain_inc
          p_strain_ten = p_strain_ten + p_strain_ten_c
          n_avg        = n_avg + cc_np1%p_strain_inc*cc_np1%u(12)
c
c                 store the CP history for this crystal
c
          call mm10_store_cryhist( cc_props, cc_np1, cc_n,
     &                             history_np1(i,co:cn) )
c
        end do ! over all crystals at int point
c
c                 finalize averages over all crystals at point.
c                    p_strain_avg -> average effective strain increment
c                                    from the average tensor
c
        rncry        = dble( ncrystals(i) )
        sig_avg      = sig_avg / rncry
        tang_avg     = tang_avg /  rncry
        slip_avg     = slip_avg / rncry
        t_work_inc   = t_work_inc / rncry
        p_work_inc   = p_work_inc / rncry
        p_strain_inc = p_strain_inc / rncry
        p_strain_ten = p_strain_ten / rncry
        t1           = p_strain_ten(1)**2 + p_strain_ten(2)**2 +
     &                 p_strain_ten(3)**2
        t2           = p_strain_ten(4)**2 + p_strain_ten(5)**2 +
     &                 p_strain_ten(6)**2
        p_strain_avg = sqrt( two/three * ( t1 + t2/two ) )
        n_avg        = n_avg / p_strain_avg / rncry
        if( locdebug ) write(iout,*) "stress", sig_avg
        if( locdebug ) write(iout,*) "tang", tang_avg
c
c                 nonlocal state values returned are creep rate
c                 wrt real time and total creep strain;
c                 ONLY for first crystal
c                     nonlocal(1) -> effective creep rate wrt real time
c                     nonlocal(2) -> (total) effective creep strain
c                     nonlocal(3) -> effective creep eponent (n)
c                     nonlocal(4) -> effective Norton constant B
c
        if( do_nonlocal ) then
          nonlocal_state(i,1) = p_strain_avg / local_work%dt
          nonlocal_state(i,2) = local_work%urcs_blk_n(i,9,gp) +
     &                          p_strain_avg
          if( n_avg .lt. one ) then ! limti range of values
            nonlocal_state(i,3) = one
          elseif( n_avg .gt. ten ) then
            nonlocal_state(i,3) = ten
          else
            nonlocal_state(i,3) = n_avg
          endif
          s_trace = (sig_avg(1) + sig_avg(2) + sig_avg(3)) / three
          se(1:6) = sig_avg(1:6)
          se(1)   = se(1) - s_trace
          se(2)   = se(2) - s_trace
          se(3)   = se(3) - s_trace
          t1      = se(1)**2 + se(2)**2 + se(3)**2
          t2      = se(4)**2 + se(5)**2 + se(6)**2
          s_trace = sqrt( three/two * ( t1 + two*t2 ) )
          B_eff   = p_strain_avg / local_work%dt/(s_trace**n_avg)
          nonlocal_state(i,4) = B_eff
        end if
c
c                 store results for integration point into
c                 variables passed by warp3d
c
        call mm10_store_gp (sig_avg,
     &     local_work%urcs_blk_n1(i,1:6,gp),
     &     tang_avg, history_np1(i,1:36),
     &     slip_avg, history_n(i,76:76+max_slip_sys-1),
     &     history_np1(i,76:76+max_slip_sys-1),
     &     t_work_inc, p_work_inc, p_strain_inc,
     &     local_work%urcs_blk_n(i,7:9,gp),
     &     local_work%urcs_blk_n1(i,7:9,gp),
     &     local_work%rot_blk_n1(i,1:9,gp),
     &     history_np1(i,64:72) )
c
      end do ! over span
c
c                  release allocated G & H arrays for
c                  certain CP constitutive models
c
      call mm10_set_cons( local_work, cc_props, 2 )
c
      return
c
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_set_cons                     *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 4/6/2016 rhd                *
c     *                                                              *
c     *       set interaction matrices G & H for slip system type    *
c     *       used for this element block                            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine mm10_set_cons( local_work, cc_props, isw )
c
      use segmental_curves, only: max_seg_points
      use mm10_defs ! to get definition of cc_props
c
      implicit integer (a-z)
$add include_sig_up
c
      type(crystal_props) :: cc_props
      integer :: isw, s_type1, n_hard
      logical :: process_G_H
c
      process_G_H = ( local_work%c_props(1,1)%h_type .eq. 4 ) .or.
     &              ( local_work%c_props(1,1)%h_type .eq. 7 )
      if( .not. process_G_H ) return
c
      select case( isw )

      case( 1 ) ! allocate G,H interaction matrices
c
         s_type1 = local_work%c_props(1,1)%s_type
         n_hard  = local_work%c_props(1,1)%num_hard
         allocate( cc_props%Gmat(n_hard,n_hard),
     &             cc_props%Hmat(n_hard,n_hard), stat=allocate_status)
         if( allocate_status .ne. 0 ) then
            write(*,*) ' error allocating G matrix'
            call die_gracefully
         endif
         call mm10_mrr_GH( s_type1, n_hard, cc_props%Gmat,
     &                     cc_props%Hmat)
c
      case( 2 ) ! deallocate G,H matrices
c
         deallocate( cc_props%Gmat, cc_props%Hmat )
c
      case default
c
          write(*,*) '>>>> FATAL ERROR. invalid isw, mm10_set_cons'
          write(*,*) '                  job terminated'
          call die_abort
c
      end select
c
      return
c
      end subroutine
c
c --------------------------------------------------------------------
c
c   Plugin Subroutines:
c   EVERY new constitutive model must be embedded into each of these
c   routines
c
c --------------------------------------------------------------------
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_init_cc_hist0                *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *    Initialize the crystal history variables                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_init_cc_hist0(props, angles, history)
      use mm10_defs
      implicit integer(a-z)
c
      type(crystal_props) :: props
      integer :: a
      double precision, dimension(3) :: angles
      double precision, dimension(3,3) :: I
      double precision,
     &  dimension(24+max_uhard+max_uhard) :: history
c
      I = 0.0d0
      do a=1,3
        I(a,a) = 1.0d0
      end do
c           Stress
      history(1:6) = 0.0d0
c           Angles
      history(7:9) = angles(1:3)
c           Rotation
      history(10:18) = reshape(I, (/9/))
c           D
      history(18+1:24) = 0.0d0
c           eps
      history(24+1:30) = 0.0d0
c           slip_incs
      history(30+1:30+max_slip_sys) = 0.0d0
c           Hardening
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! Simple voche
            call mm10_init_voche(props,
     &            history(30+max_slip_sys+1:30+max_slip_sys+max_uhard),
     & history(30+max_slip_sys+max_uhard+1:30+max_slip_sys+2*max_uhard))
      elseif (props%h_type .eq. 2) then ! MTS
            call mm10_init_mts(props,
     &            history(30+max_slip_sys+1:30+max_slip_sys+max_uhard),
     & history(30+max_slip_sys+max_uhard+1:30+max_slip_sys+2*max_uhard))
      elseif (props%h_type .eq. 3) then ! User
            call mm10_init_user(props,
     &            history(30+max_slip_sys+1:30+max_slip_sys+max_uhard),
     & history(30+max_slip_sys+max_uhard+1:30+max_slip_sys+2*max_uhard))
      elseif (props%h_type .eq. 4) then ! ORNL
            call mm10_init_ornl(props,
     &            history(30+max_slip_sys+1:30+max_slip_sys+max_uhard),
     & history(30+max_slip_sys+max_uhard+1:30+max_slip_sys+2*max_uhard))
      elseif (props%h_type .eq. 7) then ! MRR
            call mm10_init_mrr(props,
     &            history(30+max_slip_sys+1:30+max_slip_sys+max_uhard),
     & history(30+max_slip_sys+max_uhard+1:30+max_slip_sys+2*max_uhard))
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********

      return
      end subroutine

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_unknown_hard_error          *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 1/27/14                     *
c     *                                                              *
c     *     A common error message for the general hardening setup   *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_unknown_hard_error(props)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
c
      write(props%out,101) props%h_type
 101  format(
     &      10x,'>> Error: unknown hardening type ', 'i6', '.',
     &    /,10x, 'Aborting...')
      call die_gracefully

      end subroutine
c
c *****************************************************************************
c *                                                                           *
c *         WARP ROUTINES (set up or external calculations)                   *
c *                                                                           *
c *****************************************************************************
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_set_sizes_special            *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 04/01/2015 tjt              *
c     *                                                              *
c     *    called by warp3d for each material model to obtain        *
c     *    various sizes of data for the model                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_set_sizes_special( size_data, local_el  )
      use main_data, only: imatprp
      implicit integer (a-z)
$add common.main
      dimension size_data(*)
      integer :: local_el, matnum, ncrystals
c
c        size_data(1)  :  no. of words of history data for each
c                         integration point
c
c
c        in this case sizeof(__)*number of crystals
c
c
      matnum = iprops(38,local_el)
      ncrystals = imatprp(101,matnum)
c
c       So total history size is going to be:
c
c
      size_data(1) = 76+max_slip_sys+
     &               ncrystals*(30+max_slip_sys+max_uhard
     &               +max_uhard+max_uhard)
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_tangent                      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/27/13                    *
c     *                                                              *
c     *     Calculate the consistent tangent after a converged       *
c     *     stress update.                                           *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_tangent(props, np1, n, vec1, vec2,arr1, arr2,
     &     ivec1, ivec2, gaspt)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      double precision, dimension(6,6) :: J11, JJ, JR, JA, JB,JK
      double precision, dimension(6) :: d_mod, d_barp, tw
      double precision, dimension(6,props%num_hard) :: ed
      double precision, dimension(6,props%num_hard) :: J12
      double precision, dimension(props%num_hard,6) :: J21
      double precision, dimension(props%num_hard,12) :: beta
      double precision, dimension(3) :: w_p
      double precision, dimension(props%nslip,6) :: symtqmat
      double precision, dimension(6,props%nslip) :: dgammadd
      double precision,
     &    dimension(props%num_hard,props%num_hard) :: J22, alpha
      double precision :: alpha1
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
      double complex, dimension(max_uhard) :: ivec1, ivec2

      logical :: debug, gpall, locdebug
      integer, dimension(props%num_hard) :: ipiv
      integer :: i, info, gpp, gaspt
c
      debug = .false.!props%debug
      gpall = props%gpall ! true to print iteration norms for all Gauss points
      gpp = props%gpp ! set one particular G.P. to print
      locdebug = (debug .and.(gpall.or.(gaspt.eq.gpp))
     & .and.((props%st_it(3).eq.-2).or.
     &       (props%st_it(3).eq.np1%elem))
     & .and.((props%st_it(1).eq.-2).or.
     &       (props%st_it(1).eq.np1%step))
     & .and.((props%st_it(2).eq.-2).or.
     &       (props%st_it(2).eq.np1%iter)))
c
      np1%tangent = 0.0d0
c
      if(props%real_tang) then ! tangent matrix implemented
      call mm10_formvecs(props, np1, n,
     &     np1%stress, np1%tau_tilde, vec1, vec2)
      call mm10_formarrs(props, np1, n,
     &     np1%stress, np1%tau_tilde, vec1, vec2, arr1, arr2,2)
      call mm10_formJ11(props, np1, n, vec1, vec2, arr1, arr2,
     & np1%stress, np1%tau_tilde, J11)
      call mm10_formJ12(props, np1, n, vec1, vec2, arr1, arr2,
     & np1%stress, np1%tau_tilde, J12)
      call mm10_formJ21(props, np1, n, vec1, vec2, arr1, arr2,
     & np1%stress, np1%tau_tilde, J21)
      call mm10_formJ22(props, np1, n, vec1, vec2, arr1, arr2,
     & np1%stress, np1%tau_tilde, J22)
      else ! compute tangent my complex method
      call mm10_formJ11i(props, np1, n, ivec1, ivec2,
     & np1%stress, np1%tau_tilde, J11)
      call mm10_formJ12i(props, np1, n, ivec1, ivec2,
     & np1%stress, np1%tau_tilde, J12)
      call mm10_formJ21i(props, np1, n, ivec1, ivec2,
     & np1%stress, np1%tau_tilde, J21)
      call mm10_formJ22i(props, np1, n, ivec1, ivec2,
     & np1%stress, np1%tau_tilde, J22)
      endif

        if (locdebug) write (*,*) "J11", J11(1:6,1:6)
        if (locdebug) write (*,*) "J12", J12(1:6,1:props%num_hard)
        if (locdebug) write (*,*) "J21", J21(1:props%num_hard,1:6)
        if (locdebug) write (*,*) "J22", J22(1:props%num_hard
     % ,1:props%num_hard)
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
        call mm10_ed_voche(props, np1, n, np1%stress, np1%tau_tilde, ed)
      elseif (props%h_type .eq. 2) then ! MTS
        call mm10_ed_mts(props, np1, n, np1%stress, np1%tau_tilde, ed)
      elseif (props%h_type .eq. 3) then ! User
        call mm10_ed_user(props, np1, n, np1%stress, np1%tau_tilde, ed)
      elseif (props%h_type .eq. 4) then ! ORNL
        call mm10_ed_ornl(props, np1, n, np1%stress, np1%tau_tilde, ed)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_ed_mrr(props, np1, n, np1%stress, np1%tau_tilde, ed)
c        if (debug) write (*,*) "ed", ed(1:6,1:12)
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
      call mm10_symSWmat(np1%stress, np1%qc, props%nslip, symtqmat)
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
        call mm10_dgdd_voche(props, np1, n, np1%stress,
     &       np1%tau_tilde, np1%D, dgammadd)
      elseif (props%h_type .eq. 2) then ! MTS
        call mm10_dgdd_mts(props, np1, n, np1%stress,
     &       np1%tau_tilde, np1%D, dgammadd)
      elseif (props%h_type .eq. 3) then ! User
        call mm10_dgdd_user(props,np1, n, np1%stress,
     &       np1%tau_tilde, np1%D, dgammadd)
      elseif (props%h_type .eq. 4) then ! ORNL
        call mm10_dgdd_ornl(props,np1, n, np1%stress,
     &       np1%tau_tilde, np1%D, dgammadd)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_dgdd_mrr(props,np1, n, np1%stress,
     &       np1%tau_tilde, np1%D, dgammadd)
c        if (debug) write (*,*) "dgammadd", dgammadd(1:6,1:12)
      else
        call mm10_unknown_hard_error(props)
      endif
c ******* END: Add new Constitutive Models into this block *********
c
      JA = 0.0d0
      do i=1,props%nslip
        call DGER(6,6,1.d0,
     &      matmul(props%stiffness, np1%ms(1:6,i))
     &      + 2.0*symtqmat(1:6,i), 1, dgammadd(1:6,i),1,JA,6)
      end do
c        if (debug) write (*,*) "JA", JA(1,2)
c        if (debug) write (*,*) "JA", JA(1:6,1:6)
c
c Compute tangent matrix T where
c T = JJ^-1*JR
c JR = (Cijkl - JA - JB)
c JJ = J11 - J12*J22^-1*J21
c JB = J12*J22^-1*ed
c
      JJ = J11
      alpha = J22
c      call mm10_invasym(alpha, 1)
c      call DGEMM ('N','N',1,6,1,-1.d0,alpha,1,J21,1,0.d0,beta,1)
c     JJ = JJ + J12*beta
c      call DGEMM ('N','N',6,6,1,1.d0,J12,6,beta,1,1.d0,
c     &                 JJ,6)
c      call mm10_invasym(JJ, 6)
c        if (debug) write (*,*) "JJ", JJ(1,1)
c     Avoid explicitly computing the inverse
      call dcopy(props%num_hard*6,J21,1,beta,1)
      call dcopy(props%num_hard*6,ed,1,beta(1,7),1)
        if (locdebug) write (*,*) "beta", beta(1:6,1:12)
c Compute J22^-1*J21 and J22^-1*ed
      call DGESV(props%num_hard, 2*6, alpha, props%num_hard, ipiv,
     & beta, props%num_hard, info)
c     JJ = JJ - J12*beta(*,1:6)
        if (locdebug) write (*,*) "beta", beta(1:6,1:12)
      call DGEMM ('N','N',6,6,props%num_hard,-1.d0,J12,6,beta,
     &                 props%num_hard,1.d0,JJ,6)
        if (locdebug) write (*,*) "JJ", JJ(1:6,1:6)
c
c
      JB = 0.0
c      call DGEMM ('N','N',1,6,1,1.d0,alpha,1,ed,1,0.d0,beta,1)
c      call DGEMM ('N','N',6,6,1,1.d0,J12,6,beta,1,1.d0,JB,6)
c     JB = J12*beta(*,7:12)
      call DGEMM ('N','N',6,6,props%num_hard,1.d0,J12,6,beta(1,7),
     &                 props%num_hard,1.d0,JB,6)
c        if (debug) write (*,*) "JB", JB(1,2)
c
      JR = props%stiffness - JA - JB
c
c      np1%tangent = matmul(JJ, JR)
c     Avoid explicitly computing the inverse
       call DGESV(6, 6, JJ, 6, ipiv, JR,
     & 6, info)
      call dcopy(6*6,JR,1,np1%tangent,1)
        if (locdebug) write (*,*) "JR", JR(1:6,1:6)
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
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *     setup hardening for a particular state np1               *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_setup(props, np1, n)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      integer :: i, j, k, s, t
      integer, external :: mm10_l_c
      double precision, dimension(6,6) :: RE
      double precision, dimension(3,3) :: RW, RWC
      double precision, dimension(3,3,3) :: curv
      double precision, dimension(3) :: cn, tm
      double precision :: alpha
c           The alpha for geometric hardening
      parameter(alpha=1.0d0/3.0d0)
c
c           Calculate effective strain increment
      np1%dg = dsqrt(2.0d0/3.0d0*(dot_product(np1%D(1:3),np1%D(1:3))+
     &      0.5d0*dot_product(np1%D(4:6),np1%D(4:6))))
c
c           Calculate the current m and q tensors
c           Yes, these are supposed to be transposes.  We actually need the
c           backwards rotation from the lattice state.
      call mm10_RT2RVE(transpose(n%Rp), RE)
      call mm10_RT2RVW(transpose(n%Rp), RW)
      call mm10_RT2RVW(matmul(np1%R,transpose(n%Rp)), RWC)
      do i=1,props%nslip
        np1%ms(1:6,i) = matmul(RE, props%ms(1:6,i))
        np1%qs(1:3,i) = matmul(RW, props%qs(1:3,i))
        np1%qc(1:3,i) = matmul(RWC, props%qs(1:3,i))
      end do
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
        call mm10_setup_voche(props, np1, n)
      elseif (props%h_type .eq. 2) then ! MTS
        call mm10_setup_mts(props, np1, n)
      elseif (props%h_type .eq. 3) then ! User
        call mm10_setup_user(props, np1, n)
      elseif (props%h_type .eq. 4) then ! ORNL
        call mm10_setup_ornl(props, np1, n)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_setup_mrr(props, np1, n)
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
c
c    Compute quatities related to backstress, gradients, etc.
c
      if ((props%h_type .eq. 1).or.(props%h_type .eq. 2)
     &   .or.(props%h_type .eq. 3)) then ! voche, MTS, user
c           Calculate the tau lambdas for geometric hardening
c           Lattice curvature
c
      curv = 0.0
      do i=1,3
        do k=1,3
          do s=1,3
            curv(i,k,s) = n%gradFeinv(i,k,s) - n%gradFeinv(i,s,k)
          end do
        end do
      end do
c           Use Acharya's large strain definition of lambda
      do t=1,props%nslip
        ! Normal into current coordinates (This could be the problem...)
        cn = matmul(n%R, matmul(transpose(n%Rp), props%ns(1:3,t)))
c       Calculate the large-strain lambda
        tm = 0.0
        do i=1,3
          do j=1,3
            do k=1,3
              do s=1,3
                tm(i) = tm(i) + curv(i,j,k)*0.5d0 *
     &             dble(mm10_l_c(j,s,k))*cn(s)
              end do
            end do
          end do
          np1%tau_l(t) = props%k_0*props%burgers*alpha**2.0d0
     &                       *np1%mu_harden**2.0d0/
     &      (2.0d0*props%theta_0)*dsqrt(dot_product(tm,tm))
        end do
      end do
c
      elseif (props%h_type .eq. 4) then ! ORNL
c add back stress calculations here
c
      elseif (props%h_type .eq. 7) then ! MRR
      else
        call mm10_unknown_hard_error(props)
      endif
c
c
      return
      end subroutine
c
c     Helper for the above
c
      integer function mm10_l_c(i,j,k)
        implicit none
        integer :: i,j,k
        mm10_l_c = (i-j)*(j-k)*(k-i)/2
        return
      end function
cc
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_store_cryhist                *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *    Copy the state np1 struct to the history                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_store_cryhist(props, np1, n, history)
      use mm10_defs
      implicit integer(a-z)
c
      double precision,
     &  dimension(30+3*max_uhard+max_slip_sys) :: history
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      history(1:6) = np1%stress
      history(7:9) = np1%euler_angles
      history(10:18) = reshape(np1%Rp, (/9/))
      history(18+1:24) = np1%D
      history(24+1:30) = np1%eps
      history(30+1:30+max_slip_sys) = np1%slip_incs
     &       (1:max_slip_sys)
      history(30+max_slip_sys+1:30+max_slip_sys
     &        +props%num_hard) =
     &   np1%tau_tilde(1:props%num_hard)
      history(30+max_slip_sys+max_uhard+1:
     &   30+max_slip_sys+2*max_uhard) =
     &   np1%u(1:max_uhard)
      history(30+max_slip_sys+2*max_uhard+1:
     &   30+max_slip_sys+2*max_uhard+props%num_hard) =
     &   np1%tt_rate(1:props%num_hard)
c
      return
      end subroutine

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
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *     Advance a crystal from n to np1, store tangent, and      *
c     *     store other output                                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_solve_crystal(props, np1, n, cut, iout, fat,gp,
     &                  p_strain_ten_c,iter_0_extrapolate_off)
      use mm10_defs
      use main_data, only: asymmetric_assembly
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(max_uhard) :: vec1,vec2
      double precision, dimension(6) :: p_strain_ten_c
      double precision, dimension(max_uhard,max_uhard) :: arr1,arr2
      double complex, dimension(max_uhard) :: ivec1,ivec2
      logical :: cut, fat, iter_0_extrapolate_off, no_load
      integer :: iout,gp
c
      call mm10_solve_strup(props, np1, n, vec1, vec2, arr1, arr2,
     &   ivec1, ivec2, cut, gp, iter_0_extrapolate_off, no_load)
c
      if (cut) then
        write(iout,*) "mm10 stress update failed"
        return
      end if
c
c      write(*,*) 'iter_0_extrapolate_off in solv_cry',
c     &    iter_0_extrapolate_off
c     Special update with extrapolate off: return stress due to creep
c     and linear elastic stiffness matrix (from mm10_solve_strup)
      if (iter_0_extrapolate_off.or.no_load) then
          return
      end if
c
      call mm10_tangent(props, np1, n, vec1, vec2, arr1, arr2,
     &        ivec1, ivec2,gp)
c
c      call mm10_ur_tangent(props, np1, n)
c
c      call mm10_num_tangent(props, np1, n)
c
      if (.not. asymmetric_assembly .and. .not. fat) then
        np1%tangent = 0.5d0*(np1%tangent + transpose(np1%tangent))
      end if
c
      call mm10_formvecs(props, np1, n, np1%stress, np1%tau_tilde,
     &  vec1, vec2)
c
      call mm10_update_rotation(props, np1, n, vec1, vec2)
c
      call mm10_output(props, np1, n, vec1, vec2,p_strain_ten_c)
c
      return
      end subroutine
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_update_euler_angles          *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 3/28/12                     *
c     *                                                              *
c     *    update euler angles to the new rotation                   *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_update_euler_angles(props,np1,n)
            use mm10_defs
            implicit none
            type (crystal_props) :: props
            type (crystal_state) :: n, np1
c
            double precision, dimension(3,3) :: full_rot
            double precision :: psiK, phiK, thetaK, psi, phi, theta,
     &            pi, pps, pms, tol
            double precision, external :: mm10_atan2
            parameter(tol=1.0d-16)
c
            pi = 2.0d0*acos(0.0d0)
c
c                 Note: This subroutine needs a major fix if I'm
c                 ever going to support anything other than degrees+
c                 Kocks convention
            full_rot = matmul(props%g, matmul(np1%Rp, transpose(np1%R)))
c
            psiK = mm10_atan2(full_rot(3,2),full_rot(3,1))
            phiK = mm10_atan2(full_rot(2,3),full_rot(1,3))
            thetaK = dacos(full_rot(3,3))

            if (props%angle_convention .eq. 1) then
                  psi = psiK
                  phi = phiK
                  theta = thetaK
            else
                  write (*,*) "Angle convention not implemented."
                  call die_gracefully
            end if
c
            if (props%angle_type .eq. 1) then
                  np1%euler_angles(1) = 180.0d0/pi*psi
                  np1%euler_angles(2) = 180.0d0/pi*theta
                  np1%euler_angles(3) = 180.0d0/pi*phi
            elseif (props%angle_type .eq. 2) then
                  np1%euler_angles(1) = psi
                  np1%euler_angles(2) = theta
                  np1%euler_angles(3) = phi
            else
                  write (*,*) "Unrecognized angle convention."
                  call die_gracefully
            end if
c
            return
      end subroutine
c
c     Helper for the above, atan2 with range 0 to 2*pi
c
      function mm10_atan2(a, b)
            implicit none
            double precision :: mm10_atan2, a, b, pi
c
            pi = 4.d0*datan(1.d0)
            mm10_atan2 = datan2(a, b)
            if( mm10_atan2 .lt. 0.0d0 )
     &          mm10_atan2 = mm10_atan2 + 2.0d0*pi
c
            return
      end function
c
c
c ****************************************************************************
c *                                                                          *
c *    mm10_invasym                                                          *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 3/22/12 mcm                                      *
c *                                                                          *
c *         Calculate the inversion of a non-symmetric matrix using LAPACK   *
c *                                                                          *
c ****************************************************************************
c
      subroutine mm10_invasym(A,n)
            implicit none
            double precision, intent(inout), dimension(n,n) :: A
            integer, intent(in) :: n

            integer :: i,j,info,lwork
            integer, allocatable :: ipivt(:)
            double precision, allocatable :: work(:)
c
c           Allocate storage
            allocate(ipivt(n))
            lwork = n*n
            allocate(work(lwork))
c           Factor
            call DGETRF(n,n,A,n,ipivt,info)
c           Inverse
            call DGETRI(n,A,n,ipivt,work,lwork,info)
c           Free storage
            deallocate(ipivt)
            deallocate(work)

            return
       end subroutine
c
c
c ****************************************************************************
c *                                                                          *
c *    mm10_rotation_matrix                                                       *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 3/22/12 mcm                                      *
c *                                                                          *
c *         Given euler angles, an angle convention, and an angle type       *
c *         send back the correct rotation matrix.                           *
c *                                                                          *
c ****************************************************************************
c
      subroutine mm10_rotation_matrix(angles, aconv, atype, R, out)
            implicit none
            double precision, dimension(3), intent(in) :: angles
            character, intent(in) :: aconv*5
            character, intent(in) :: atype*7
            integer, intent(in) :: out
c
            double precision, dimension(3,3), intent(out) :: R
c
            double precision :: a, b, c, psi, theta, phi, pi
c
            pi = 2d0*dacos(0.0d0)
c
            a = angles(1)
            b = angles(2)
            c = angles(3)

            if (atype .eq. 'degrees') then
                  a = a*pi/180.d0
                  b = b*pi/180.d0
                  c = c*pi/180.d0
            elseif (atype .eq. 'radians') then
            else
                  write (out,9000)
            end if

            if (aconv .eq. 'kocks') then
                  psi = a
                  theta = b
                  phi = c
            elseif (aconv .eq. 'bunge') then
                  psi = a - pi/2.d0
                  theta = b
                  phi = pi/2.d0 - c
            elseif (aconv .eq. 'roe') then
                  psi = a
                  theta = b
                  phi = 3.d0*pi/2.d0-c
            else
                  write (out,9001)
            end if


            R(1,1) = -sin(psi)*sin(phi)-cos(psi)*cos(phi)*cos(theta)
            R(1,2) = cos(psi)*sin(phi)-sin(psi)*cos(phi)*cos(theta)
            R(1,3) = cos(phi)*sin(theta)
            R(2,1) = sin(psi)*cos(phi)-cos(psi)*sin(phi)*cos(theta)
            R(2,2) = -cos(psi)*cos(phi)-sin(psi)*sin(phi)*cos(theta)
            R(2,3) = sin(phi)*sin(theta)
            R(3,1) = cos(psi)*sin(theta)
            R(3,2) = sin(psi)*sin(theta)
            R(3,3) = cos(theta)

            return
 9000 format(/'Danger: Unknown angle type passed to rotation_matrix'/)
 9001 format(/'Danger: Unknown angle convention passed to',
     &        ' rotation_matrix'/)
      end subroutine
c
c
c ****************************************************************************
c *                                                                          *
c *    mm10_invsym                                                                *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 3/22/12 mcm                                      *
c *                                                                          *
c *         Calculate the inversion of a symmetric matrix using LAPACK       *
c *                                                                          *
c ****************************************************************************
c
      subroutine mm10_invsym(A,n)
            implicit none
            double precision, intent(inout), dimension(n,n) :: A
            integer, intent(in) :: n

            integer :: i,j,info,lwork
            integer, allocatable :: ipiv(:)
            double precision, allocatable :: work(:)
c
c           Allocate storage
            allocate(ipiv(n))
            lwork = n*n
            allocate(work(lwork))
c           Factor
            call DSYTRF('U',n,A,n,ipiv,work,lwork,info)
c           Inverse
            call DSYTRI('U',n,A,n,ipiv,work,info)
c           Sym -> Full
            do i=1,n
                  do j=1,i-1
                        A(i,j) = A(j,i)
                  end do
             end do
c           Free storage
            deallocate(ipiv)
            deallocate(work)

            return
       end subroutine
c
c
c ****************************************************************************
c *                                                                          *
c *    mm10_RT2RVE                                                                *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 3/22/12 mcm                                      *
c *                                                                          *
c *         Takes a 3x3 rotation tensor and returns it in a 6x6 form         *
c *         suitable for rotating Voigt-type strain vectors                  *
c *                                                                          *
c ****************************************************************************
c
      subroutine mm10_RT2RVE( RT, RV )
            implicit none
            double precision, dimension(3,3), intent(in) :: RT
            double precision, dimension(6,6), intent(out) :: RV

            RV(1,1)=RT(1,1)**2.d0
            RV(1,2)=RT(1,2)**2.d0
            RV(1,3)=RT(1,3)**2.d0
            RV(1,4)=2.d0*RT(1,1)*RT(1,2)
            RV(1,5)=2.d0*RT(1,3)*RT(1,2)
            RV(1,6)=2.d0*RT(1,1)*RT(1,3)
            RV(2,1)=RT(2,1)**2.d0
            RV(2,2)=RT(2,2)**2.d0
            RV(2,3)=RT(2,3)**2.d0
            RV(2,4)=2*RT(2,1)*RT(2,2)
            RV(2,5)=2*RT(2,3)*RT(2,2)
            RV(2,6)=2*RT(2,1)*RT(2,3)
            RV(3,1)=RT(3,1)**2.d0
            RV(3,2)=RT(3,2)**2.d0
            RV(3,3)=RT(3,3)**2.d0
            RV(3,4)=2.d0*RT(3,1)*RT(3,2)
            RV(3,5)=2.d0*RT(3,3)*RT(3,2)
            RV(3,6)=2.d0*RT(3,1)*RT(3,3)
            RV(4,1)=RT(1,1)*RT(2,1)
            RV(4,2)=RT(1,2)*RT(2,2)
            RV(4,3)=RT(1,3)*RT(2,3)
            RV(4,4)=RT(1,1)*RT(2,2)+RT(2,1)*RT(1,2)
            RV(4,5)=RT(1,2)*RT(2,3)+RT(1,3)*RT(2,2)
            RV(4,6)=RT(1,1)*RT(2,3)+RT(1,3)*RT(2,1)
            RV(5,1)=RT(2,1)*RT(3,1)
            RV(5,2)=RT(3,2)*RT(2,2)
            RV(5,3)=RT(2,3)*RT(3,3)
            RV(5,4)=RT(2,1)*RT(3,2)+RT(2,2)*RT(3,1)
            RV(5,5)=RT(2,2)*RT(3,3)+RT(3,2)*RT(2,3)
            RV(5,6)=RT(2,1)*RT(3,3)+RT(2,3)*RT(3,1)
            RV(6,1)=RT(1,1)*RT(3,1)
            RV(6,2)=RT(1,2)*RT(3,2)
            RV(6,3)=RT(1,3)*RT(3,3)
            RV(6,4)=RT(1,1)*RT(3,2)+RT(1,2)*RT(3,1)
            RV(6,5)=RT(1,2)*RT(3,3)+RT(1,3)*RT(3,2)
            RV(6,6)=RT(1,1)*RT(3,3)+RT(3,1)*RT(1,3)

            return
      end subroutine
c
c
c ****************************************************************************
c *                                                                          *
c *    mm10_RT2RVW                                                                *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 3/22/12 mcm                                      *
c *                                                                          *
c *         Takes a 3x3 rotation tensor and returns it in a 3x3 form         *
c *         suitable for rotating my 3x1 skew vectors                        *
c *                                                                          *
c ****************************************************************************
c
      subroutine mm10_RT2RVW( RT, RV )
            implicit none
            double precision, dimension(3,3), intent(in) :: RT
            double precision, dimension(3,3), intent(out) :: RV

            RV(1,1)=RT(2,2)*RT(3,3)-RT(2,3)*RT(3,2)
            RV(1,2)=RT(2,1)*RT(3,3)-RT(2,3)*RT(3,1)
            RV(1,3)=RT(2,1)*RT(3,2)-RT(2,2)*RT(3,1)
            RV(2,1)=RT(1,2)*RT(3,3)-RT(1,3)*RT(3,2)
            RV(2,2)=RT(1,1)*RT(3,3)-RT(1,3)*RT(3,1)
            RV(2,3)=RT(1,1)*RT(3,2)-RT(1,2)*RT(3,1)
            RV(3,1)=RT(1,2)*RT(2,3)-RT(1,3)*RT(2,2)
            RV(3,2)=RT(1,1)*RT(2,3)-RT(1,3)*RT(2,1)
            RV(3,3)=RT(1,1)*RT(2,2)-RT(1,2)*RT(2,1)

            return
      end subroutine
c
c ****************************************************************************
c *                                                                          *
c *    mm10_ET2EV                                                                 *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 3/22/12 mcm                                      *
c *                                                                          *
c *         Strain tensor to strain vector                                   *
c *                                                                          *
c ****************************************************************************
c
      subroutine mm10_ET2EV(ET, EV)
            implicit none
            double precision, dimension(3,3), intent(in) :: ET
            double precision, dimension(6), intent(out) :: EV

            EV(1) = ET(1,1)
            EV(2) = ET(2,2)
            EV(3) = ET(3,3)
            EV(4) = 2.d0*ET(1,2)
            EV(5) = 2.d0*ET(2,3)
            EV(6) = 2.d0*ET(1,3)

            return
      end subroutine
c
c
c ****************************************************************************
c *                                                                          *
c *    mm10_WT2WV                                                                 *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 3/22/12 mcm                                      *
c *                                                                          *
c *         Skew   tensor to skew   vector                                   *
c *                                                                          *
c ****************************************************************************
c
      subroutine mm10_WT2WV(WT, WV)
            implicit none
            double precision, dimension(3,3), intent(in) :: WT
            double precision, dimension(3), intent(out) :: WV

            WV(1) = WT(2,3)
            WV(2) = WT(1,3)
            WV(3) = WT(1,2)

            return
      end subroutine
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_store_gp                     *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *    Store all required gauss point data                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_store_gp(stress_in, stress_out, tang_in,
     &      tang_out, slip_inc, slip_n, slip_np1, t_work, p_work,
     &      p_strain, u_old, u_new, R_in, R_out)
      use mm10_defs
      implicit integer(a-z)
c
      double precision, dimension(6) :: stress_in, stress_out
      double precision, dimension(9) :: R_in, R_out
      double precision, dimension(6,6) :: tang_in
      double precision, dimension(36) :: tang_out
      double precision, dimension(max_slip_sys) :: slip_inc, slip_n,
     &      slip_np1
      double precision :: t_work, p_work, p_strain
      double precision, dimension(3) :: u_old, u_new
c
      stress_out(1:6) = stress_in(1:6)
      tang_out(1:36) = reshape(tang_in, (/ 36 /))
      slip_np1 = slip_n + slip_inc
c
      u_new(1) = u_old(1) + t_work
      u_new(2) = u_old(2) + p_work
      u_new(3) = u_old(3) + p_strain
c
      R_out(1:9) = R_in(1:9)
      return
      end subroutine
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_calc_grads                   *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 12/20/13                    *
c     *                                                              *
c     *    calculate the gradient of Re.T (Fe) through a linear      *
c     *    fit                                                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_calc_grads(ngp, elem_type, order, geonl, rot_blk,
     &            jac, Rps, gradFes)
      implicit none
      integer :: ngp, elem_type, order
      double precision, dimension(9,ngp) :: rot_blk, Rps
      double precision, dimension(27,ngp) :: gradFes
      double precision, dimension(3,3) :: jac
      logical :: geonl
c
      integer :: i, a, b, lwork, info
      double precision, dimension(ngp,3,3) :: Rt
      double precision, dimension(3,3) :: jacinv
      double precision, dimension(ngp,4) :: intermat
      double precision, dimension(ngp) :: RHS
      double precision :: weight, fact
      double precision, dimension(8) :: work
      double precision, dimension(3,3,3) :: grads
c
c
c           Get R components and stick in the right place
      if (geonl) then
        jacinv = jac
        call mm10_invasym(jacinv, 3)
        do i=1, ngp
          Rt(i,1:3,1:3) = matmul(reshape(Rps(1:9,i),(/3,3/)),
     &        transpose(reshape(rot_blk(1:9,i), (/3,3/))))
        end do
      else
        do i=1, ngp
          Rt(i,1:3,1:3) = reshape(Rps(1:9,i), (/3,3/))
        end do
      end if
c
c     For each Rt component create an interpolation, solve for the
c     coefficients, and store the gradient
      do a=1,3
        do b=1,3
          intermat = 0.0d0
          RHS = 0.0d0
          do i=1,ngp
c           1-3 are the coordinates
            call getgpts( elem_type, order, i, intermat(i,1),
     &            intermat(i,2), intermat(i,3), weight)
            intermat(i,4) = 1.0D0
            RHS(i) = Rt(i,a,b)
          end do
c           Solve with LAPACK
          lwork = 8
          call DGELS('N',  ngp, 4, 1, intermat, ngp, RHS, ngp, work,
     &            lwork, info)
c           Extract coefs
          if (info .ne. 0) then
            write (*,*) "Error finding least squares mm10 grad."
          end if
c           Get the gradient
          grads(a,b,1) = RHS(1)
          grads(a,b,2) = RHS(2)
          grads(a,b,3) = RHS(3)
c           Take to the current coordinates
          if (geonl) then
            grads(a,b,1:3) = matmul(jacinv,grads(a,b,1:3))
          end if
        end do
      end do
c
c     Flatten and store
      do i=1,ngp
        gradFes(1:27, i) = reshape(grads(1:3,1:3,1:3), (/27/))
      end do
c
      return
      end subroutine

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_setup_np1                   *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *    Initialize the state np1 structure with new strain/temp/  *
c     *     time increment                                           *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_setup_np1(Rur, dstrain, dt, T, step, elem,
     &      iter,gp,np1)
      use mm10_defs
      implicit none
c
      double precision, dimension(9) :: Rur
      double precision, dimension(6) :: dstrain
      double precision :: dt, T
      integer :: step, elem, gp, iter
      type(crystal_state) :: np1
c
      np1%R = reshape(Rur(1:9), (/3,3/))
      np1%D = dstrain(1:6)
      np1%temp = T
      np1%tinc = dt
      np1%step = step
      np1%elem = elem
      np1%iter = iter
      np1%gp = gp
c
      end subroutine

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_init_uout_hist               *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *    initialize the user output                                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_init_uout_hist(history)
      implicit none
      double precision :: history(3)
c
      history(1:3) = 0.0d0
c
      return
      end subroutine

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
      subroutine mm10_init_slip_hist(history)
      implicit integer (a-z)
$add param_def
      double precision :: history(max_slip_sys)
c
      history(1:max_slip_sys) = 0.0d0
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_init_cc_props                *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *    Copy properties over from local_work into the update      *
c     *    structure                                                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_init_cc_props(inc_props, atype, aconv, debug,
     &                              cc_props)
      use mm10_defs
      implicit integer (a-z)
$add include_sig_up
      integer :: atype, aconv
      logical :: debug
      type(crystal_properties) :: inc_props
      type(crystal_props) :: cc_props
c
c     Just a whole lot of copying
      cc_props%rate_n = inc_props%rateN
      cc_props%tau_hat_y = inc_props%tauHat_y
      cc_props%G_0_y = inc_props%Go_y
      cc_props%burgers = inc_props%burgers
      cc_props%p_v = inc_props%p_v
      cc_props%q_v = inc_props%q_v
      cc_props%boltzman = inc_props%boltzman
      cc_props%theta_0 = inc_props%theta_o
      cc_props%eps_dot_0_v = inc_props%eps_dot_o_v
      cc_props%eps_dot_0_y = inc_props%eps_dot_o_y
      cc_props%p_y = inc_props%p_y
      cc_props%q_y = inc_props%q_y
      cc_props%tau_a = inc_props%tau_a
      cc_props%tau_hat_v = inc_props%tauHat_v
      cc_props%G_0_v = inc_props%Go_v
      cc_props%k_0 = inc_props%k_o
      cc_props%mu_0 = inc_props%mu_o
      cc_props%D_0 = inc_props%D_o
      cc_props%T_0 = inc_props%t_o
      cc_props%tau_y = inc_props%tau_y
      cc_props%tau_v = inc_props%tau_v
      cc_props%voche_m = inc_props%voche_m
      cc_props%u1 = inc_props%u1
      cc_props%u2 = inc_props%u2
      cc_props%u3 = inc_props%u3
      cc_props%u4 = inc_props%u4
      cc_props%u5 = inc_props%u5
      cc_props%u6 = inc_props%u6
      cc_props%solver = inc_props%solver
      cc_props%strategy = inc_props%strategy
      cc_props%gpall = inc_props%gpall
      cc_props%gpp = inc_props%gpp
      cc_props%st_it(1:3) = inc_props%st_it(1:3)
      cc_props%method = inc_props%method
      cc_props%miter = inc_props%miter
      cc_props%atol = inc_props%atol
      cc_props%atol1 = inc_props%atol1
      cc_props%rtol = inc_props%rtol
      cc_props%rtol1 = inc_props%rtol1
      cc_props%xtol = inc_props%xtol
      cc_props%xtol1 = inc_props%xtol1
c
      cc_props%g(1:3,1:3) = inc_props%rotation_g(1:3,1:3)
c
      cc_props%ms(1:6,1:max_slip_sys) =
     &      inc_props%ms(1:6,1:max_slip_sys)
      cc_props%qs(1:3,1:max_slip_sys) =
     &      inc_props%qs(1:3,1:max_slip_sys)
      cc_props%ns(1:3,1:max_slip_sys) =
     &      inc_props%ns(1:3,1:max_slip_sys)
c
      cc_props%stiffness(1:6,1:6) =
     &      inc_props%init_elast_stiff(1:6,1:6)
c
      cc_props%angle_type = atype
      cc_props%angle_convention = aconv
      cc_props%nslip = inc_props%nslip

      cc_props%h_type = inc_props%h_type
      cc_props%num_hard = inc_props%num_hard
      cc_props%real_tang = inc_props%real_tang
      cc_props%debug = debug
      cc_props%s_type = inc_props%s_type
      cc_props%cnum = inc_props%cnum
c
      return
      end subroutine

c
c *****************************************************************************
c *                                                                           *
c *         User hardening routines                                           *
c *                                                                           *
c *****************************************************************************
c
c           Initialize history
      subroutine mm10_init_user(props, tau_tilde, uhist)
      use mm10_defs
      implicit integer(a-z)
c
      type(crystal_props) :: props
      double precision, dimension(props%num_hard) :: tau_tilde
      double precision, dimension(max_uhard) :: uhist
c
      write (*,*) "Not implemented"
      call die_gracefully

      return
      end subroutine
c
c           Setup user hardening
      subroutine mm10_setup_user(props, np1, n)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      write(*,*) "Not implemented"
      call die_gracefully
c
      return
      end subroutine
c -----------------------
c     Simple voche:
c
c -----------------------
c           Initialize history
      subroutine mm10_init_voche(props, tau_tilde, uhist)
      use mm10_defs
      implicit integer(a-z)
c
      type(crystal_props) :: props
      double precision, dimension(max_uhard) :: tau_tilde
      double precision, dimension(max_uhard) :: uhist
      double precision :: init_hard
      parameter(init_hard = 0.1d0)
c
      tau_tilde(1) = props%tau_y+init_hard
      ! Nothing with the user history
      return
      end subroutine
c
c           Setup voche law hardening
      subroutine mm10_setup_voche(props, np1, n)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      ! No setup actually required, but define a mu_harden at state np1
      ! for the sake of the CP model
      np1%mu_harden = props%stiffness(6,6)
c
      return
      end subroutine
c -------------
c     MTS:
c
c -------------
c           Initialize history
      subroutine mm10_init_mts(props, tau_tilde, uhist)
      use mm10_defs
      implicit integer(a-z)
c
      type(crystal_props) :: props
      double precision, dimension(max_uhard) :: tau_tilde
      double precision, dimension(max_uhard) :: uhist
c
      tau_tilde(1) = -1.0d0 ! This only works because these are actually flags
      uhist(1) = -1.0d0
      uhist(2) = -1.0d0
c
      return
      end subroutine
c
c           Setup MTS hardening
      subroutine mm10_setup_mts(props, np1, n)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      double precision :: dgc, init_hard
      parameter(init_hard=0.1d0)
c
c     Continuum effective rate
      dgc = np1%dg / np1%tinc
c
c     New shear modulus
      np1%mu_harden = props%mu_0 - props%D_0 / (exp(props%T_0/np1%temp)
     &      - 1.0d0)
c
c     Get the threshold contributions
      np1%tau_v = props%tau_hat_v*(1.d0-(props%boltzman*np1%temp
     &       /(np1%mu_harden*(props%burgers**3.d0)*props%G_0_v)*
     &       dlog(props%eps_dot_0_v/dgc))**(1.d0/props%q_v))
     &       **(1.d0/props%p_v)
      np1%tau_y = props%tau_hat_y*(1.d0-(props%boltzman*np1%temp
     &       /(np1%mu_harden*(props%burgers**3.d0)*props%G_0_y)*
     &       dlog(props%eps_dot_0_y/dgc))**(1.d0/props%q_y))
     &       **(1.d0/props%p_y)
c
c     Used existing labels as a convenience, actually get/set the history
      np1%u(1) = np1%tau_y
      if (n%u(1) .lt. 0.0d0) then
        n%tau_y = np1%tau_y
      else
        n%tau_y = n%u(1)
      end if

      np1%u(2) = np1%mu_harden
      if (n%u(2) .lt. 0.0d0) then
        n%mu_harden = np1%mu_harden
      else
        n%mu_harden = n%u(2)
      end if
c
c     Same here -- check for previous step flag
      if (n%tau_tilde(1) .lt. 0.0d0) then
        n%tau_tilde(1) = props%tau_a +
     &     (np1%mu_harden/props%mu_0)*np1%tau_y + init_hard
      end if
c
      return
      end subroutine
c
c *****************************************************************************
c *                                                                           *
c *         Ma-Roters-Raabe hardening routines                                *
c *                                                                           *
c *****************************************************************************
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
c           Initialize history
      subroutine mm10_init_mrr(props, tau_tilde, uhist)
      use mm10_defs
      implicit integer(a-z)
c
      type(crystal_props) :: props
      double precision :: tau_tilde(props%num_hard)
      double precision, dimension(max_uhard) :: uhist
c
      if(props%theta_0.eq.0.d0) then
      tau_tilde(1:props%num_hard) = 1.0d8 ! Initial densities for edges
      else
      tau_tilde(1:props%num_hard) = props%theta_0 ! Initial densities for edges
      endif
c
      return
      end subroutine
c
c           Setup mrr hardening
      subroutine mm10_setup_mrr(props, np1, n)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision :: time
c increment the total time
      time = n%u(1) + np1%tinc
      np1%u(1) = time
c
      return
      end subroutine
c
c *****************************************************************************
c *                                                                           *
c *         ORNL ferritic-martensitic steel hardening routines                *
c *                                                                           *
c *****************************************************************************
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
c           Initialize history
      subroutine mm10_init_ornl(props, tau_tilde, uhist)
      use mm10_defs
      implicit integer(a-z)
c
      type(crystal_props) :: props
      double precision :: tau_tilde(props%num_hard)
      double precision, dimension(max_uhard) :: uhist
c
      if(props%theta_0.eq.0.d0) then
      tau_tilde(1:props%num_hard) = 1.0d8 ! Initial densities for edges
      else
      tau_tilde(1:props%num_hard) = props%theta_0 ! Initial densities for edges
      endif
c
      return
      end subroutine
c
c           Setup ornl hardening
      subroutine mm10_setup_ornl(props, np1, n)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision :: time
c increment the total time
      time = n%u(1) + np1%tinc
      np1%u(1) = time
c
      return
      end subroutine
c
c ****************************************************************************
c *                                                                          *
c *    mm10.f                                                                *
c *                                                                          *
c *         written by : mcm                                                 *
c *                                                                          *
c *         Solver functions                                                 *
c *                                                                          *
c ****************************************************************************
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10c_unknown_hard_error          *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 1/27/14                     *
c     *                                                              *
c     *     A common error message for the general hardening setup   *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10c_unknown_hard_error(props)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
c
      write(props%out,101) props%h_type
 101  format(
     &      10x,'>> Error: unknown hardening type ', 'i6', '.',
     &    /,10x, 'Aborting...')
      call die_gracefully

      end subroutine
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
      subroutine mm10_copy_cc_hist(history, gradfe, R, props, n)
      use mm10_defs
      implicit integer(a-z)
c
      double precision,
     &  dimension(24+max_uhard+max_uhard) :: history
      double precision, dimension(27) :: gradfe
      double precision, dimension(9) :: R
      type(crystal_props) :: props
      type(crystal_state) :: n
c
c           Not provided, could be useful...
      n%R(1:3,1:3) = reshape(R, (/3,3/))
      n%gradFeinv(1:3,1:3,1:3) = reshape(gradfe, (/3,3,3/))
      n%temp = 0.0d0
c           Only used at n+1
      n%tinc = 0.0d0
      n%dg = 0.0d0
      n%tau_v = 0.0d0
      n%tau_y = 0.0d0
      n%mu_harden = 0.0d0
c
      n%stress = history(1:6)
      n%euler_angles(1:3) = history(7:9)
c
      n%Rp = reshape(history(10:18), (/3,3/))
      n%D(1:6) = history(18+1:24)
      n%eps(1:6) = history(24+1:30)
      n%slip_incs(1:max_slip_sys) =
     &  history(30+1:30+max_slip_sys)
      n%tau_tilde(1:props%num_hard) =
     &  history(30+max_slip_sys+1:30+max_slip_sys
     &          +props%num_hard)
c
      n%u(1:max_uhard) =
     &  history(30+max_slip_sys+max_uhard+1:30+max_slip_sys
     &          +2*max_uhard)
      n%tt_rate(1:props%num_hard) =
     &  history(30+max_slip_sys+2*max_uhard+1:30+max_slip_sys
     &          +2*max_uhard+props%num_hard)
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_init_general_hist            *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *    initialize general GP history (grad Fe, tangent, R)       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_init_general_hist(history)
      implicit none
      double precision :: history(72)
      double precision :: eye(3,3)
c
      history(1:63) = 0.0d0
c
      eye = 0.0d0
      eye(1,1) = 1.0d0
      eye(2,2) = 1.0d0
      eye(3,3) = 1.0d0
      history(64:72) = reshape(eye, (/9/))
      return
      end subroutine

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_solve_strup                  *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *     Solve the stress update adaptively (if required)         *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_solve_strup(props, np1, n, vec1, vec2, arr1,
     &   arr2, ivec1, ivec2, fail, gp, iter_0_extrapolate_off,
     &    no_load)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      logical :: fail
c
      type(crystal_state) :: curr
      double precision, dimension(6) :: stress, ostress, R1
      double precision :: tt(props%num_hard), ott(props%num_hard)
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
      double complex, dimension(max_uhard) :: ivec1, ivec2
      integer :: cuts, mcuts, i, gp, gpp
      double precision :: frac, step, mult
      integer, dimension(10) :: faili
      double precision, dimension(10) :: failr
      double precision :: temp1,temp2
      logical :: debug, gpall, locdebug, iter_0_extrapolate_off,
     &    no_load
c
      parameter(mcuts = 4)
      parameter(mult = 0.5d0)
      debug = props%debug
      gpall = props%gpall ! true to print iteration norms for all Gauss points
      gpp = props%gpp ! set one particular G.P. to print
c      Set flag for specific step, element, iteration, ...
      locdebug = (debug .and.(gpall.or.(gp.eq.gpp))
     & .and.((props%st_it(3).eq.-2).or.
     &       (props%st_it(3).eq.np1%elem))
     & .and.((props%st_it(1).eq.-2).or.
     &       (props%st_it(1).eq.np1%step))
     & .and.((props%st_it(2).eq.-2).or.
     &       (props%st_it(2).eq.np1%iter)))
c
      frac = 0.0d0
      step = 1.0d0
      cuts = 0
c
      stress = n%stress
      tt = n%tau_tilde
      ostress = stress
      ott = tt
      if(locdebug) then ! print statement for debugging
      write(props%out,*) " n%tt_rate =", n%tt_rate(1)
      endif
c
      call mm10_setup(props, np1, n)
      fail = .false.
c     check if zero increments or zero stress
      temp1 = np1%D(1)*np1%D(1)+np1%D(2)*np1%D(2)+
     &        np1%D(3)*np1%D(3)+np1%D(4)*np1%D(4)+
     &        np1%D(5)*np1%D(5)+np1%D(6)*np1%D(6)
      temp2 = stress(1)*stress(1)+stress(2)*stress(2)+
     &        stress(3)*stress(3)+stress(4)*stress(4)+
     &        stress(5)*stress(5)+stress(6)*stress(6)
      no_load = (temp2.eq.0.d0).and.(temp1.eq.0.d0)
c       write(*,*) 'noload', no_load

c     Branch on update for stress: either elasticity, or plasticity
      if(iter_0_extrapolate_off.or.no_load) then
c
        np1%tangent = props%stiffness ! assign elastic stiffness as the tangent matrix
c       check whether to update the stress accounting for evolving creep
c       according to np1%stress = n%stress - R1(n%stress,uddt)
c       where R1 = stress - n%stress - matmul(props%stiffness, np1%D - dbarp)
c     &      + 2.0d0 * symTW
        if(.not.no_load) then
          call mm10_formvecs(props, np1, n, stress, tt, vec1, vec2)
          call mm10_formR1(props, np1, n, vec1, vec2, stress, tt,
     &                     R1, gp)
          stress = stress - R1
        end if
c ensure that zero rates are set, though this value should not be kept in warp3d history during extrapolaion anyway
        curr%tt_rate = 0.d0
c
      else ! update stress, state variables using usual material N-R solvers
c
      do while (frac .lt. 1.0d0)
        call mm10_setup_np1(reshape(np1%R, (/9/)), np1%D*(step+frac),
     &      np1%tinc*(step+frac), (np1%temp-n%temp)*(step+frac)+n%temp,
     &      np1%step, np1%elem, np1%iter, np1%gp, curr)
        call mm10_setup(props, curr, n)
        tt = n%tau_tilde ! some subroutines modify n%tau_tilde in setup
        if(props%solver) then
        call mm10_solve(props, curr, n, vec1, vec2, arr1, arr2,
     &       ivec1, ivec2, stress, tt, fail, faili, failr, gp,
     &       np1%tinc*step)
        else
        call mm10_solveB(props, curr, n, vec1, vec2, arr1, arr2,
     &       ivec1, ivec2, stress, tt, fail, faili, failr, gp,
     &       np1%tinc*step)
        endif
        if (fail) then
          if (locdebug) write(*,*) "Adapting"
          stress = ostress
          tt = ott
          step = step * mult
          cuts = cuts + 1
          if (cuts .gt. mcuts) exit
          fail = .false.
        else
          ostress = stress
          ott = tt
          frac = frac + step
        end if
      end do
c
c
      if (fail .or. any(isnan(tt)) .or.
     &   any(isnan(stress))) then
      write(props%out,*)" >>> Warning: mm10 implicit solution failed."
        if(faili(4).eq.1) then
      write(props%out,*)" Stress prediction failed at iter=",
     &  faili(1), " for miter=", faili(2)
        else
      write(props%out,*)" Material update failed at iter=",
     &  faili(1), " for miter=", faili(2)
        endif
        if(faili(3).eq.1) then
      write(props%out,*)" Reason: absolute residual norm"
        elseif(faili(3).eq.2) then
      write(props%out,*)" Reason: relative residual norm"
        else
      write(props%out,*)" Reason: encountered NaN"
        endif
      write(props%out,*)" AbsNorm=",
     &  failr(1), " AbsTol=", failr(2)
      write(props%out,*)" RelNorm=",
     &  failr(3), " RelTol=", failr(4)
      write(props%out,*)" Error occured in: element=",
     &  np1%elem, " GaussPoint=", np1%gp
      write(props%out,*)" Fractional step frac=",
     &  frac, " step=", step, " cuts=", cuts
      fail = .true.
      np1%stress = n%stress
      np1%tau_tilde = n%tau_tilde
      return
      end if ! failure from solver
c
      end if
c
c
      np1%stress = stress
      np1%tau_tilde = tt
      np1%tt_rate = curr%tt_rate ! store rates of hardening too
      if(locdebug) then ! print statement for debugging
      write(props%out,*)" stress(2)=",
     &  stress(2), " np1%stress=", np1%stress(2)
      write(props%out,*)"  tt=",
     &  tt(1), " np1%tau_tilde =", np1%tau_tilde(1)
      write(props%out,*) " np1%tt_rate =", np1%tt_rate(1)
      endif
c
      return
      end subroutine
c
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_solve                        *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *     Solve a stress update to prescribed strain state         *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_solve(props, np1, n, vec1, vec2, arr1, arr2,
     &          ivec1, ivec2, stress, tt, fail,faili,failr, gp,dtinc)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      logical :: fail
c
      double precision, dimension(6+props%num_hard) :: R,x,dx, xnew,g
      double precision,
     & dimension(6+props%num_hard,6+props%num_hard) :: J
      double precision :: nR, inR, atol, rtol, uB, alpha, ls1, ls2,
     &      nlsx, nRs, c, red, dt, cos_ang, xetol, xtol1, zerotol,
     &      dxerr
      integer :: iter, miter, info, ls, mls, mmin, gp, gpp, ttind
      integer, dimension(6+props%num_hard) :: ipiv
      logical :: debug, gpall, solver, strategy, locdebug
      double precision, dimension(6) :: R1, x1, dx1, xnew1, d1,d2
      double precision, dimension(props%num_hard) :: x2
      double precision, dimension(6,6) :: J11
      double precision, dimension(max_uhard) :: vec1,vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1,arr2
      double complex, dimension(max_uhard) :: ivec1, ivec2
      double precision :: nR1, inR1, atol1, rtol1, dtinc
      integer, dimension(10) :: faili
      double precision, dimension(10) :: failr
c      Trust region solver
      integer maxfev,nfev,njev,lr, i
      double precision xtol,factor, mm10_enorm
      double precision rout(6+props%num_hard,6+props%num_hard)
      double precision, dimension(6+props%num_hard) :: diag,qtf,
     *      wa1,wa2,wa3,wa4
c      Cubic line search
      double precision stepmx
c
c      Convergence parameters: Newton with geometric line search
      parameter(c = 1.0d-4)
      parameter(red = 0.5d0)
      parameter(mls = 10)
      parameter(mmin = 1)
      atol = props%atol
      atol1 = props%atol1
      rtol = props%rtol
      rtol1 = props%rtol1
      xetol = props%xtol
      xtol1 = props%xtol1
      miter = props%miter
      zerotol = 1.d-12
      ttind = 1 ! index of tt to print while debugging
c       Trust region parameters
      xtol = 1.0d-2
      maxfev = 3*miter
      lr=((6+props%num_hard)*(6+props%num_hard+1))/2
c      Debug flags for printing iteration behavior
      debug = props%debug
      gpall = props%gpall ! true to print iteration norms for all Gauss points
      gpp = props%gpp ! set one particular G.P. to print
c      Solver flags
      solver = props%solver ! true for Mark's N.R. routine, false for trust-region
      strategy = props%strategy ! true for geometric l.s., false for cubic l.s.
c      Set flag for specific step, element, iteration, ...
      locdebug = (debug .and.(gpall.or.(gp.eq.gpp))
     & .and.((props%st_it(3).eq.-2).or.
     &       (props%st_it(3).eq.np1%elem))
     & .and.((props%st_it(1).eq.-2).or.
     &       (props%st_it(1).eq.np1%step))
     & .and.((props%st_it(2).eq.-2).or.
     &       (props%st_it(2).eq.np1%iter)))
c
      if (locdebug) write(*,*) "Entering solution routine"
c
      x(1:6) = stress
      x(7:props%num_hard+6) = tt
c
      if(locdebug) then ! print statement for debugging
      write(props%out,*)" Guess syy=",
     &  x(2)
      endif
c
c Prediction of yield stress to initialize the integration algorithm; helps
c for Orowan flow rule type models
      if (props%h_type .gt. 3) then
c
c      if (debug) write(*,*) "Entering stress prediction"
      iter = 0
      x1 = x(1:6)
c
      ! Module to extrapolate the hardening variables
      dt = dtinc !np1%tinc
      ! Predict the new hardening variable by extrapolation
      ! Use cosine of the "angle" between the new and old
      ! displacement increment to indicate continuity of load
      ! direction
      d1(1:6) = np1%D(1:6)
      alpha = d1(1) + d1(2) + d1(3)
      d1(1) = d1(1) - 1.d0/3.d0*alpha
      d1(2) = d1(2) - 1.d0/3.d0*alpha
      d1(3) = d1(3) - 1.d0/3.d0*alpha
      d1(1:6) = d1(1:6)/dsqrt(dot_product(d1,d1))
      d2(1:6) = n%D(1:6)
      alpha = d2(1) + d2(2) + d2(3)
      d2(1) = d2(1) - 1.d0/3.d0*alpha
      d2(2) = d2(2) - 1.d0/3.d0*alpha
      d2(3) = d2(3) - 1.d0/3.d0*alpha
      d2(1:6) = d2(1:6)/dsqrt(dot_product(d2,d2))
      cos_ang = dmax1(dot_product(d1,d2),0.d0)
      !write(*,*) cos_ang, n%tt_rate(2)
      x2 = x(7:6+props%num_hard) + cos_ang*
     &        n%tt_rate(1:props%num_hard)*dt

      if(locdebug)
     & write(props%out,*)" Stress prediction module, G.P.=", gp
      if(locdebug) then ! print statement for debugging
      write(props%out,*)" Extrapol tt6=",
     &  x2(6), " Previous tt6=", x(6+6)
      endif
      call mm10_formvecs(props, np1, n, x1, x2, vec1, vec2)
      call mm10_formR1(props, np1, n, vec1, vec2, x1,x2, R1, gp)
      nR1 = dsqrt(dot_product(R1,R1))
      inR1 = nR1
c       Newton-Raphson loop
      if (locdebug) write(*,'("Iter ",i3," norm ",E10.3)') iter, nR1
c      dxerr = 1.1*xtol1
c      if(locdebug) then ! print statement for debugging
c      write(props%out,*)" xtol1=",
c     &  xtol1, " dxerr=", dxerr, 'T/f', ((nR1 .gt. atol1)
c     & .and. (nR1/inR1 .gt. rtol1)
c     &   .and. (dxerr .gt. xtol1))
c      endif
c      do while ((nR1 .gt. atol1) .and. (nR1/inR1 .gt. rtol1)
c     &   .and. (dxerr .gt. xtol1))
      do while ((nR1 .gt. atol1) .and. (nR1/inR1 .gt. rtol1)
     &   )
c           Jacobian
        if(props%real_tang) then ! tangent matrix implemented
        call mm10_formarrs(props, np1, n, x1, x2, vec1, vec2,
     &       arr1, arr2,1)
        call mm10_formJ11(props, np1, n, vec1, vec2,
     &       arr1, arr2, x1, x2, J11)
        else
        call mm10_formJ11i(props, np1, n, ivec1, ivec2,
     &       x1, x2, J11)
        endif
c        if(locdebug) then
c          write(*,*) 'vec1', vec1(1:6)
c        endif
c        if (locdebug) write (*,*) "R1", R1(1:6)
c        if (locdebug) write (*,*) "J11", J11(1:6,1:6)
c           Increment
        dx1 = R1
        call DGESV(6, 1, -J11, 6, ipiv, dx1, 6, info)
      if(locdebug) then ! print statement for debugging
      write(props%out,*)" Iter=",
     &  iter, " dx1=", dx1(2)
      endif

c           Line search
        alpha = 1.0d0
        ls1 = 0.5d0*dot_product(R1,R1)
        ls2 = c*dot_product(dx1, matmul(transpose(J11),R1))
        ls = 0
        do
          nlsx = ls1 + ls2*alpha
          xnew1 = x1 + alpha*dx1
c             Residual
          call mm10_formvecs(props, np1, n, xnew1, x2, vec1, vec2)
          call mm10_formR1(props, np1, n, vec1, vec2,
     &         xnew1,x2, R1,gp)
        if (locdebug) write (*,*) "R1 ls", R1(2)
          nR1 = dsqrt(dot_product(R1,R1))
          nRs = 0.5d0*dot_product(R1,R1)
c             Line search convergence test
          if ((nRs .le. nlsx) .or. (ls .gt. mls)) then
c         compute incremental solution vector error
c        dxerr = 0.d0
c        do i = 1,6
c          if (dabs(x1(i)).gt.zerotol) then
c            dxerr = dxerr + dabs(alpha*dx1(i)/x1(i))
c          endif
c        end do
c        do i = 1,props%num_hard
c          if (dabs(x2(i)).gt.zerotol) then
c            dxerr = dxerr + dabs(alpha*dx2(i)/x2(i))
c          endif
c        end do
            x1 = xnew1
      if(locdebug) then ! print statement for debugging
      write(props%out,*)" Iter=",
     &  iter, " syy=", x1(2), " AbsNorm=", nR1, " alpha=", alpha
      endif
            exit
          else
            alpha = red*alpha
            ls = ls + 1
          end if
        end do
c           Increment and check for failure
        iter = iter + 1
        if (locdebug) write(*,'("Iter ",i3," norm ",E10.3," ls",F10.3,
     &   " dx",F10.3)')
     &            iter, nR1, alpha, dxerr
c
        if ((iter .gt. miter) .or. any(isnan(x1))) then
c         Record data and reason for failure
          fail = .true.
          faili(1) = iter
          faili(2) = miter
          if((nR1 .gt. atol1)) then
          faili(3) = 1
          elseif((nR1/inR1 .gt. rtol1)) then
          faili(3) = 2
          elseif(any(isnan(x1))) then
          faili(3) = 3
          endif
          faili(4) = 1
          failr(1) = nR1
          failr(2) = atol1
          failr(3) = nR1/inR1
          failr(4) = rtol1
          return
        end if
c
      end do
c       Output statistics from N-R algorithm
      if(locdebug) then ! print statement for debugging
      write(props%out,*)" Stress pred conv iter=",
     &  iter
      write(props%out,*)" AbsNorm=",
     &  nR1, " AbsTol=", atol1
      write(props%out,*)" RelNorm=",
     &  nR1/inR1, " RelTol=", rtol1
      write(props%out,*)" Guess syy=",
     &  x(2), " actual syy=", x1(2)
      endif
c        Copy predicted stress and hardening back to primary variable x
      x(1:6) = x1(1:6)
      x(7:6+props%num_hard) = x2(1:props%num_hard)
c
      else
c
        inR1 = 1.d0
c
      endif ! stress prediction

c
c Material update algorithm
c
      if(locdebug) ! print statement for debugging
     & write(props%out,*)" Material update module, G.P.=", gp
      if(locdebug) then ! print statement for debugging
      write(props%out,*)" Guess syy=",
     &  x(2), " Guess tt=", x(6+ttind)
      endif
c
c
      iter = 0
      call mm10_formR(props, np1, n, vec1, vec2, x(1:6),
     & x(7:6+props%num_hard), R, gp)
c      write (*,*) "R1", R(1:6)
c      if (debug) write (*,*) "R2", R(7:6+props%num_hard)
      nR = dsqrt(dot_product(R,R))
      inR = nR
      if(inR.eq.0.d0) then
        inR = inR1
      endif
c       Newton-Raphson loop
      if (locdebug) write(*,'("Iter ",i3," norm ",E10.3)') iter, nR
c      dxerr = 1.1*xtol
c      do while (((nR .gt. atol) .and. (nR/inR .gt. rtol)
c     &   .and. (dxerr .gt. xetol)) .or.
      do while (((nR .gt. atol) .and. (nR/inR .gt. rtol)
     &   ) .or.
     &          (iter.lt.mmin))
c           Jacobian
        if(props%real_tang) then ! tangent matrix implemented
        call mm10_formJ(props, np1, n, vec1, vec2, arr1, arr2, x(1:6),
     & x(7:6+props%num_hard), J)
        else
        call mm10_formJi(props, np1, n, ivec1, ivec2, x(1:6),
     & x(7:6+props%num_hard), J)
        endif
c        if (debug) write (*,*) "J11", J(1:6,1:6)
c        if (debug) write (*,*) "J12", J(1:6,7:6+props%num_hard)
c        if (debug) write (*,*) "J21", J(7:6+props%num_hard,1:6)
c        if (debug) write (*,*) "J22",
c     &    J(7:6+props%num_hard,7:6+props%num_hard)
c           Increment
        dx = R
        call DGESV(6+props%num_hard, 1, -J, 6+props%num_hard, ipiv, dx,
     & 6+props%num_hard, info)
c
      if(locdebug) then ! print statement for debugging
      write(props%out,*)" Iter=",
     &  iter, " dx1=", dx(2), " dx2=", dx(6+ttind)
      endif
c
c           Line search
        alpha = 1.0d0
        ls1 = 0.5d0*dot_product(R,R)
        ls2 = c*dot_product(dx, matmul(transpose(J),R))
        ls = 0
        do
          nlsx = ls1 + ls2*alpha
          xnew = x + alpha*dx
c             Residual
          call mm10_formR(props, np1, n, vec1, vec2, xnew(1:6),
     & xnew(7:6+props%num_hard), R, gp)
c      if (debug) write (*,*) "R1", R(1:6)
c      if (debug) write (*,*) "R2", R(7:6+props%num_hard)
          nR = dsqrt(dot_product(R,R))
          nRs = 0.5d0*dot_product(R,R)
c             Line search convergence test
          if ((nRs .le. nlsx) .or. (ls .gt. mls)) then
c         compute incremental solution vector error
c        dxerr = 0.d0
c        do i = 1,6+props%num_hard
c          if (dabs(x(i)).gt.zerotol) then
c            dxerr = dxerr + dabs(alpha*dx(i)/x(i))
c          endif
c        end do
            x = xnew
      if(locdebug) then ! print statement for debugging
      write(props%out,*)" L.S. converged, Iter=",
     &  iter, " syy=", x(2), " tt(5)=", x(6+ttind), " nR=", nR
      write(props%out,*)" alpha=",
     &  alpha
c      write(props%out,*)" dxerr=",
c     &  dxerr
      endif
            exit
          else
            alpha = red*alpha
            ls = ls + 1
          end if
        end do
c
c           Increment and check for failure
        iter = iter + 1
        if (locdebug) write(*,'("Iter ",i3," norm ",E10.3," ls",F10.3)')
     &            iter, nR, alpha
c         Record data and reason for failure
        if ((iter .gt. miter) .or. any(isnan(x))) then
          fail = .true.
          faili(1) = iter
          faili(2) = miter
          if((nR .gt. atol)) then
          faili(3) = 1
          elseif((nR/inR .gt. rtol)) then
          faili(3) = 2
          elseif(any(isnan(x))) then
          faili(3) = 3
          endif
          faili(4) = 2
          failr(1) = nR
          failr(2) = atol
          failr(3) = nR/inR
          failr(4) = rtol
          return
        end if

      end do
c
c       Output statistics from N-R algorithm
      if(locdebug) then ! print statement for debugging
      write(props%out,*)" Material upd conv, iter=",
     &  iter
      write(props%out,*)" AbsNorm=",
     &  nR, " AbsTol=", atol
      write(props%out,*)" RelNorm=",
     &  nR/inR, " RelTol=", rtol
      write(props%out,*)" Guess syy=",
     &  stress(2), " actual syy=", x(2)
      write(props%out,*)" Guess tt6=",
     &  x2(6), " actual tt6=", x(6+ttind)
      endif
c           Set for return
      stress = x(1:6)
      tt = x(7:6+props%num_hard)
c      if (debug) write (*,*) "stress", x(1:6)
c      write (*,*) "fail", fail
c      write (*,*) "iter", iter

      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_update_rotation              *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 1/14/14                     *
c     *                                                              *
c     *     Update the plastic rotation                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_update_rotation(props, np1, n, vec1, vec2)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      double precision, dimension(3) :: wbarp
      double precision, dimension(3,3) :: wbarp_full, expw
      double precision, dimension(max_uhard) :: vec1, vec2
c
      call mm10_form_wbarp(props, np1, n, vec1, vec2,
     &      np1%stress, np1%tau_tilde,
     &      wbarp)
      call mm10_WV2WT(wbarp, wbarp_full)
c
      call mm10_expmw3x3(wbarp_full, expw)
c
      np1%Rp = matmul(expw, n%Rp)
c
      end subroutine
c
c ****************************************************************************
c *                                                                          *
c *    mm10_WV2WT                                                                 *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 3/22/12 mcm                                      *
c *                                                                          *
c *         Skew vector to skew tensor                                       *
c *                                                                          *
c ****************************************************************************
c
      subroutine mm10_WV2WT(WV, WT)
            implicit none
            double precision, dimension(3,3), intent(out) :: WT
            double precision, dimension(3), intent(in) :: WV

            WT(1,1) = 0.0d0
            WT(1,2) = WV(3)
            WT(1,3) = WV(2)
            WT(2,1) = -WV(3)
            WT(2,2) = 0.d0
            WT(2,3) = WV(1)
            WT(3,1) = -WV(2)
            WT(3,2) = -WV(1)
            WT(3,3) = 0.d0

            return
      end subroutine
c
c ****************************************************************************
c *                                                                          *
c *    mm10_expmw3x3                                                              *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 3/22/12 mcm                                      *
c *                                                                          *
c *         Calculates exp(W) where W is a 3x3 skew matrix.                  *
c *         Returns full 3x3 because the result                              *
c *         is only orthogonal (not skew or symmetric)                       *
c *                                                                          *
c ****************************************************************************
c
      subroutine mm10_expmw3x3(W, A)
            implicit none
            double precision, dimension(3,3), intent(in) :: W
            double precision, dimension(3,3), intent(out) :: A
            double precision :: alpha
            integer :: i
c
c           Compute alpha
            alpha = DSQRT(W(2,3)**2.d0+W(1,3)**2.d0+W(1,2)**2.d0)
c
c           Algorithm will fail with alpha = 0 (=> W=0)
c           Correct result is expm(0) = I, so program that in
            if (alpha .lt. 1.0d-16) then
                  A = 0.0d0
            else
                  A=W
                  call DGEMM('n','n',3,3,3,(1.d0-dcos(alpha))/
     &                 (alpha**2.d0),W,3,W,3,sin(alpha)/alpha,A,3)
            end if

c           Add the identity
            do i=1,3
                  A(i,i) = A(i,i) + 1.0D0
            end do

            return
      end subroutine
c234567
