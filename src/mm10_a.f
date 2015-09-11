c
c ****************************************************************************
c *                                                                          *
c *    mm10.f                                                                *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 1/3/2015 rhd                                     *
c *                                                                          *
c *         Stress/strain update routines AND HELPERS for crystal plasticity *
c *         material modeled via Beaudoin et al.                             *
c *                                                                          *
c ****************************************************************************
c

c
c     ****************************************************************
c     *                                                              *
c     *                 module mm10_defs                             *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 12/11/13                    *
c     *                                                              *
c     *    small module to hold crystal update data structs          *
c     *                                                              *
c     ****************************************************************
c
      module mm10_defs
            implicit integer (a-z)
$add param_def
c           Massive list of properties
            type :: crystal_props
                  double precision :: rate_n, tau_hat_y, G_0_y, burgers,
     &                  p_v, q_v, boltzman, theta_0, eps_dot_0_v,
     &                  eps_dot_0_y,
     &                  p_y, q_y, 
     &                  tau_a, tau_hat_v, G_0_v,
     &                  k_0, mu_0, D_0, T_0, tau_y, tau_v, voche_m,
     &                  u1, u2, u3, u4, u5, u6
                  double precision :: atol, atol1, rtol, rtol1
                  double precision, dimension(3,3) :: g
                  double precision :: ms(6,max_slip_sys),
     &                  qs(3,max_slip_sys), ns(3,max_slip_sys)
                  double precision, dimension(6,6) :: stiffness
                  integer :: angle_type, angle_convention, nslip, 
     &                  h_type, miter, gpp, s_type, cnum
                  integer :: num_hard
                  logical :: real_tang, solver, strategy, debug, gpall
                  integer :: out
            end type
            type :: crystal_state
                  double precision, dimension(3,3) :: R, Rp
                  double precision, dimension(6) :: stress, D, eps
                  double precision, dimension(3) :: euler_angles
                  double precision, dimension(max_slip_sys) :: tau_l,
     &                  slip_incs
                  double precision, dimension(3,3,3) :: gradFeinv
                  double precision, dimension(6,6) :: tangent
                  double precision, dimension(max_uhard) :: tau_tilde
                  double precision, dimension(max_uhard) :: tt_rate
                  double precision :: temp, tinc,
     &                   dg, tau_v, tau_y,
     &                   mu_harden, work_inc, p_work_inc, p_strain_inc
                  double precision :: ms(6,max_slip_sys),
     &                  qs(3,max_slip_sys), qc(3,max_slip_sys)
                  double precision, dimension(max_uhard) :: u
                  integer :: step, elem, gp
            end type
      end module

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10                              *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 12/11/13                     *
c     *                                                              *
c     *    crystal plasticity stress-strain update                   *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10(gp, span, ncrystals, hist_sz,
     &            history_n, history_np1,
     &            local_work, uddt, gp_temps, 
     &            gp_temp_inc, iout)
            use segmental_curves, only: max_seg_points
            use mm10_defs
            implicit integer (a-z)
$add include_sig_up
            double precision, intent(in) :: uddt(mxvl,nstr)
            double precision, intent(in) :: gp_temps(mxvl),
     &                                      gp_temp_inc(mxvl)
            integer, intent(in) :: iout
            integer :: span, ncrystals(mxvl), hist_sz, gp
            double precision :: history_n(span,hist_sz)
            double precision :: history_np1(span,hist_sz)
c     
            logical :: debug
            integer :: i,c,co,cn
            double precision, dimension(6) :: sig_avg
            double precision, dimension(6,6) :: tang_avg
            double precision, dimension(max_slip_sys) :: slip_avg
            double precision :: t_work_inc, p_work_inc,p_strain_inc
            type(crystal_props) :: cc_props
            type(crystal_state) :: cc_n, cc_np1
c
            debug = .false.
c
            if (debug) write (*,*) "In mm10"
c
c                 Loop on gauss points
            local_work%material_cut_step = .false.
            do i=1,span
c                       Initialize the element history if it's step
                  if (local_work%step .eq. 1) then
                        if (debug) write(*,*) "Init GP history"
                        call mm10_init_general_hist(history_n(i,1:72))
                        call mm10_init_uout_hist(history_n(i,73:75))
                        call mm10_init_slip_hist(history_n(i,76:
     &                        76+max_slip_sys-1))
                  end if
c                       Fix a problem with the rotations
c                       I feel this really should be done elsewhere
                  if (.not. local_work%geo_non_flg) then
                    local_work%rot_blk_n1(i, 1:9, gp) = 
     &                  (/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 
     &                    0.d0, 0.d0, 1.d0/)
                  end if

                  if (debug) write (*,*) "Updating local ", i
                  sig_avg = 0.d0
                  slip_avg = 0.d0
                  t_work_inc = 0.d0
                  p_work_inc = 0.d0
                  p_strain_inc = 0.d0
                  tang_avg = 0.d0
c                       Loop on crystals
                  do c=1,ncrystals(i)
                      co = 76+max_slip_sys+(c-1)*(30+max_slip_sys
     &                     +3*max_uhard)
                      cn = 76+max_slip_sys+(c)*(30+max_slip_sys
     &                     +3*max_uhard)-1
                      if (debug) write(*,*) "Setting up properties"
                      call mm10_init_cc_props(local_work%c_props(i,c),
     &                        local_work%angle_type(i), 
     &                        local_work%angle_convention(i),
     &                        local_work%debug_flag,cc_props)
                      cc_props%out = iout
c
                      if (local_work%step .eq. 1) then
                        if (debug) write(*,*) "Init history 0"
                        call mm10_init_cc_hist0(cc_props, 
     &                     local_work%c_props(i,c)%init_angles(1:3),
     &                     history_n(i,co:cn))
                      end if
                      if (debug) write(*,*) "Copying n to struct"
                      call mm10_copy_cc_hist(history_n(i,co:cn),
     &                   history_n(i,37:63),
     &                   history_n(i,64:72),
     &                   cc_props, cc_n)
c
                      call mm10_setup_np1(
     &                  local_work%rot_blk_n1(i,1:9,gp), uddt(i,1:6), 
     &                  local_work%dt, gp_temps(i), local_work%step,
     &                  i+local_work%felem, gp, cc_np1)
c
                      if (debug) write (*,*) "Updating crystal ", c
                      call mm10_solve_crystal(cc_props, cc_np1, cc_n,
     &                  local_work%material_cut_step, iout, .false.,gp)
c
                      if (local_work%material_cut_step) return
c                       
c                     Add stuff into the average
                      sig_avg = sig_avg + cc_np1%stress
                      tang_avg = tang_avg + cc_np1%tangent
                      slip_avg = slip_avg + cc_np1%slip_incs
                      t_work_inc = t_work_inc + cc_np1%work_inc
                      p_work_inc = p_work_inc + cc_np1%p_work_inc
                      p_strain_inc = p_strain_inc + cc_np1%p_strain_inc
c
c                     Store the CP history for this crystal
                      call mm10_store_cryhist(cc_props, cc_np1, cc_n, 
     &                  history_np1(i,co:cn))
                  end do
c
c                 Do the division for the averages
                  sig_avg = sig_avg / DBLE(ncrystals(i))
                  tang_avg = tang_avg / DBLE(ncrystals(i))
                  slip_avg = slip_avg / DBLE(ncrystals(i))
                  t_work_inc = t_work_inc / DBLE(ncrystals(i))
                  p_work_inc = p_work_inc / DBLE(ncrystals(i))
                  p_strain_inc = p_strain_inc / DBLE(ncrystals(i))
c
c                 Actually store the GP data
c
                  call mm10_store_gp(sig_avg, 
     &               local_work%urcs_blk_n1(i,1:6,gp),
     &               tang_avg, history_np1(i,1:36),
     &               slip_avg, history_n(i,76:76+max_slip_sys-1),
     &               history_np1(i,76:76+max_slip_sys-1),
     &               t_work_inc, p_work_inc, p_strain_inc, 
     &               local_work%urcs_blk_n(i,7:9,gp),
     &               local_work%urcs_blk_n1(i,7:9,gp),
     &               local_work%rot_blk_n1(i,1:9,gp),
     &               history_np1(i,64:72))
            end do
c
            return
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
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_output                       *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/27/13                    *
c     *                                                              *
c     *     Calculate various other user output quantities           *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_output(props, np1, n,vec1,vec2)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      integer :: i, info, sysID, numAct
      double precision, dimension(6) :: ee, dbarp, ewwe, ep, eeunrot
      double precision, dimension(3) :: wp
      double precision, dimension(6,6) :: S, erot
      double precision, dimension(max_uhard) :: vec1,vec2
      double precision :: maxslip, ec_dot, n_eff, rs, ec_slip, 
     &                    mm10_rs
      double precision, dimension(max_slip_sys) :: dgammadtau
c
c     Store the slip increments
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
        do i=1,props%nslip
           call mm10_slipinc(props, np1, n, np1%stress,
     &                        np1%tau_tilde, i, np1%slip_incs(i))
        end do
      elseif (props%h_type .eq. 2) then ! MTS
        do i=1,props%nslip
           call mm10_slipinc(props, np1, n, np1%stress,
     &                        np1%tau_tilde, i, np1%slip_incs(i))
        end do
      elseif (props%h_type .eq. 3) then ! User
        do i=1,props%nslip
           call mm10_slipinc_user(props, np1, n, np1%stress,
     &                        np1%tau_tilde, i, np1%slip_incs(i))
        end do
      elseif (props%h_type .eq. 7) then ! MRR
c        do i=1,props%nslip
c           call mm10_slipinc_mrr(props, np1, n, 
c     &           np1%stress,np1%tau_tilde, i, np1%slip_incs(i))
c        end do
        np1%slip_incs(1:props%nslip)  = vec1(1:props%nslip)
      else
        call mm10_unknown_hard_error(props)
      end if
c ********* END: Add new Constitutive Models into this block *********
c
c     Call a function to store the Euler angles
c
      call mm10_update_euler_angles(props, np1, n)
c
c     Compute user history output
c
c     6:8 is about the active slip systems: maximum slip rate,
c     identifier of that system, and how many systems have rates
c     on the same order of magnitude as that one.
      maxslip = 0.d0
      sysID = 0
      do i=1,props%nslip
        if(dabs(np1%slip_incs(i)).gt.maxslip) then
          maxslip = dabs(np1%slip_incs(i))
          sysID = i
        endif
      end do
      np1%u(6) = maxslip/np1%tinc
      np1%u(7) = dble(sysID)
c
      numAct = 0
      do i=1,props%nslip
        if(dabs(np1%slip_incs(i)).ge.0.1d0*maxslip) then
          numAct = numAct + 1
        endif
      end do
      np1%u(8) = dble(numAct)
c
c     Take care of miscellaneous things
c
      np1%work_inc = dot_product(np1%stress, np1%D)
c     Plastic strain and work
      S = props%stiffness
c      call mm10_invsym(S, 6)
c     Avoid forming the inverse explicitly
      call dcopy(6, np1%stress, 1, eeunrot, 1)
      call DPOSV( 'U', 6, 6, S, 6, eeunrot, 6, INFO )
      call mm10_RT2RVE(np1%R, erot)
      ee = matmul(erot, eeunrot)
c      ee = matmul(erot, matmul(S, np1%stress))
      call mm10_form_dbarp(props, np1, n,vec1,vec2,
     &      np1%stress, np1%tau_tilde, 
     &      dbarp)
      call mm10_form_wp(props, np1, n,vec1,vec2,
     &      np1%stress, np1%tau_tilde,
     &      wp)
      call mm10_symSW(ee, wp, ewwe)
c
      ep = dbarp + ewwe
      np1%p_strain_inc = dsqrt(2.0d0/3.0d0*
     &      (dot_product(ep(1:3),ep(1:3))+
     &      0.5d0*dot_product(ep(4:6),ep(4:6))))
      np1%p_work_inc = dot_product(np1%stress, ep)
c
c     And finally, store the lattice strains
c
      np1%eps = ee
c
c     Compute plastic creep rate and effective creep exponent
c
      ec_dot = np1%p_strain_inc/np1%tinc
      np1%u(11) = ec_dot
      if(ec_dot.gt.0.d0) then
      ep(1:6) = ep(1:6)/np1%tinc
c
c Generalization of CP model implementation for other slip rate
c equations, requiring other forms of d_gamma/d_tau
c Vector dgammadtau should be 1 x n_slip
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
        call mm10_dgdt_voche(props,np1, n, np1%stress,
     &         np1%tau_tilde, dgammadtau)
      elseif (props%h_type .eq. 2) then ! MTS
        call mm10_dgdt_mts(props, np1,n, np1%stress,
     & np1%tau_tilde, dgammadtau)
      elseif (props%h_type .eq. 3) then ! User
        call mm10_dgdt_user(props,np1, n, np1%stress,
     & np1%tau_tilde, dgammadtau)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_dgdt_mrr(props,np1, n, np1%stress,
     & np1%tau_tilde, dgammadtau)
      else
        call mm10_unknown_hard_error(props)
      endif
c ******* END: Add new Constitutive Models into this block *********
      n_eff = 0.d0
        do i=1,props%nslip
           rs = mm10_rs(props, np1, n, np1%stress,
     & np1%tau_tilde, i)
           ec_slip = (dot_product(np1%ms(1:3,i),ep(1:3))+
     &      0.5d0*dot_product(np1%ms(4:6,i),ep(4:6)))
           n_eff = n_eff + 2.d0/3.d0/ec_dot/ec_dot*rs*
     &             dgammadtau(i)/np1%tinc*ec_slip ! formula according to David Parks
        end do
      else
      n_eff = 1.d0
      endif
      np1%u(12) = n_eff

      return
      end subroutine
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
     &     ivec1, ivec2)
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

      logical :: debug
      integer, dimension(props%num_hard) :: ipiv
      integer :: i, info
c
      debug = .false.
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

        if (debug) write (*,*) "J11", J11(1:6,1:6)
        if (debug) write (*,*) "J12", J12(1:6,1:props%num_hard)
        if (debug) write (*,*) "J21", J21(1:props%num_hard,1:6)
        if (debug) write (*,*) "J22", J22(1:props%num_hard
     % ,1:props%num_hard)
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
        call mm10_ed_voche(props, np1, n, np1%stress, np1%tau_tilde, ed)
      elseif (props%h_type .eq. 2) then ! MTS
        call mm10_ed_mts(props, np1, n, np1%stress, np1%tau_tilde, ed)
      elseif (props%h_type .eq. 3) then ! User
        call mm10_ed_user(props, np1, n, np1%stress, np1%tau_tilde, ed)
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
        if (debug) write (*,*) "beta", beta(1:6,1:12)
c Compute J22^-1*J21 and J22^-1*ed
      call DGESV(props%num_hard, 2*6, alpha, props%num_hard, ipiv, 
     & beta, props%num_hard, info)
c     JJ = JJ - J12*beta(*,1:6)
        if (debug) write (*,*) "beta", beta(1:6,1:12)
      call DGEMM ('N','N',6,6,props%num_hard,-1.d0,J12,6,beta,
     &                 props%num_hard,1.d0,JJ,6)
        if (debug) write (*,*) "JJ", JJ(1:6,1:6)
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
        if (debug) write (*,*) "JR", JR(1:6,1:6)
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
      subroutine mm10_solve_crystal(props, np1, n, cut, iout, fat,gp)
      use mm10_defs
      use main_data, only: asymmetric_assembly
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(max_uhard) :: vec1,vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1,arr2
      double complex, dimension(max_uhard) :: ivec1,ivec2
      logical :: cut, fat
      integer :: iout,gp
c
      call mm10_solve_strup(props, np1, n, vec1, vec2, arr1, arr2,
     &   ivec1, ivec2, cut,gp)
c      
      if (cut) then
        write(iout,*) "mm10 stress update failed"
        return
      end if
c
      call mm10_tangent(props, np1, n, vec1, vec2, arr1, arr2,
     &        ivec1, ivec2)
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
      call mm10_output(props, np1, n, vec1, vec2)
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
c *    mm10_invasym                                                               *
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
c
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
c *****************************************************************************
c *                                                                           *
c *         WARP ROUTINES (set up or external calculations)                   *
c *                                                                           *
c *****************************************************************************
c
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
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine cnst10                            *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 3/22/12                     *
c     *                                                              *
c     *    gets the consistent tangent for a block                   *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine cnst10(gp, span, hist_sz,
     &            history_n, history_np1,
     &            local_work, iout)
      use segmental_curves, only: max_seg_points
      implicit integer (a-z)
$add param_def
$add include_tan_ek
      integer, intent(in) :: iout
      integer :: span, ncrystals(mxvl), hist_sz, gp
      double precision :: history_n(span,hist_sz)
      double precision :: history_np1(span,hist_sz)
c
      integer :: i
c
      do i=1,span
        local_work%cep(i,:,:) = 
     &            reshape(history_np1(i,1:36), (/6,6/))
        local_work%cep(i,:,:) = local_work%cep(i,:,:)*
     &                  local_work%det_jac_block(i,gp)*
     &                  local_work%weights(gp) 
      end do
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine lnstff10                          *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 3/22/12                     *
c     *                                                              *
c     *   Get the linearized tangent for a block                     *
c     *   TODO: include the elastic rotation                         *
c     *                                                              *
c     ****************************************************************
c
      subroutine lnstff10(gp, span, local_work)
c
      implicit integer (a-z)
$add param_def
$add include_lin_ek
c
      integer :: gp, span
c
      integer :: i, ci, tc, a, b
      double precision, dimension(6,6) :: totalC, Cci, Srot, Ct
      double precision, dimension(3,3) :: g
c
      do i = 1, span
      tc = 0
      totalC = 0.0d0
      do ci = 1, local_work%ncrystals(i)
            g = local_work%cp_g_rot(i,1:3,1:3,ci)
            Ct = local_work%cp_stiff(i,1:6,1:6,ci)
            Srot = 0.0d0
            call mm10_RT2RVE(transpose(g), Srot)
            Cci = matmul(Ct, transpose(Srot))
            Cci = matmul(Srot, Cci)
            totalC = totalC + Cci
            tc = tc + 1
      end do
      totalC = totalC / DBLE(tc)
      do b = 1, 6
       do a = 1, 6
        local_work%cep(i,a,b) = totalC(a,b) * local_work%weights(gp) *
     &      local_work%det_jac_block(i, gp)
       end do
      end do
      end do


      return

      end subroutine
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_set_e_nu                     *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 7/11/12                     *
c     *                                                              *
c     *    An annoying extra helper to enable the linear stiffness   *
c     *    routines to play nicely with the CP model                 *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_set_e_nu(matnum,elnum,e,nu)
      use main_data, only : matprp, lmtprp, imatprp, dmatprp, smatprp
      use crystal_data, only : c_array, angle_input, crystal_input,
     &                              data_offset
c
      implicit integer (a-z)
$add param_def
c
      integer :: matnum, elnum
      real, intent(out) :: e, nu
c
c     Local
      integer :: c, cnum, ncry, osn

      e = 0.0d0
      nu = 0.0d0
      ncry = imatprp(101,matnum)
      do c=1,ncry
c                 Get the local crystal number
            if (imatprp(104,matnum) .eq. 1) then
                  cnum = imatprp(105,matnum)
            elseif (imatprp(104,matnum) .eq. 2) then
                  osn = data_offset(elnum) 
                  cnum = crystal_input(osn,c)
c                 Couldn't do this earlier, so check here
                  if ((cnum .gt. max_crystals) .or. 
     &                  (cnum .lt. 0)) then
                   write (*,'("Crystal ", i3, " not valid")')
     &                  cnum
                        call die_gracefully
                  end if
            else
                  write(*,*) "INSANITY IN mm10" 
                  call die_gracefully
            end if
                  e = e + 
     &             SNGL(c_array(cnum)%e)
                  nu = nu
     &             + SNGL(c_array(cnum)%nu)
      end do
      e = e / 
     &       SNGL(ncry)
      nu = nu /
     &       SNGL(ncry)
      return
      end subroutine
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
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_set_sizes                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 12/14/14 rhd                *
c     *                                                              *
c     *    called by warp3d for each material model to obtain        *
c     *    various sizes of data for the model                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_set_sizes( info_vector )
      implicit integer (a-z)
$add common.main
      integer, dimension(*) :: info_vector
c
c        set infor_data
c
c         1        number of history values per integration 
c                  point. Abaqus calles these "statev". Values
c                  double or single precsion based on hardware.
c    
c         2        number of values in the symmetric part of the 
c                  [D] for each integration point. for solid
c                  elements this is 21, for cohesive elements this 6.
c
c         3        = 0, the material model returns "unrotated"
c                       Cauchy stresses at n+1
c                  = 1, the material model returns the standard
c                       Cauchy stresses at n+1
c
c         4        number of state variables per point to be output
c                  when user requests this type of results
c
      info_vector(1) = -100     ! set by special version above
      info_vector(2) = -100    ! set by special version above
      info_vector(3) = -100    ! set by special version above 
      info_vector(4) = 41+max_slip_sys+max_uhard
c
      return
      end

      
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_labels                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 12/15/2014 (rhd)               *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, outi,
     &      comment_lines, max_comment_lines, num_comment_lines )
      implicit integer (a-z)
$add common.main
c
c                       parameters
c
      integer :: size_state, num_states, outi, max_comment_lines,
     &           num_comment_lines
      character(len=8)  :: state_labels(size_state)
      character(len=60) :: state_descriptors(size_state)
      character(len=80) :: comment_lines(max_comment_lines)
c
c                       locals
c
      integer :: i
      logical, save :: do_print = .false.
c
      num_states = 41 + max_slip_sys + max_uhard
      num_comment_lines = 0
c     
      state_labels(1) = "euler-1"
      state_labels(2) = "euler-2"
      state_labels(3) = "euler-3"

      state_descriptors(1) = "convention as input"
      state_descriptors(2) = "convention as input"
      state_descriptors(3) = "convention as input"

      state_labels(4) = "Rp-11"
      state_labels(5) = "Rp-21"
      state_labels(6) = "Rp-31"
      state_labels(7) = "Rp-12"
      state_labels(8) = "Rp-22"
      state_labels(9) = "Rp-32"
      state_labels(10) = "Rp-13"
      state_labels(11) = "Rp-23"
      state_labels(12) = "Rp-33"

      state_descriptors(4:12) = "plastic rotation component"

      state_labels(13) = "mcFei-11"
      state_labels(14) = "mcFei-21"
      state_labels(15) = "mcFei-31"
      state_labels(16) = "mcFei-12"
      state_labels(17) = "mcFei-22"
      state_labels(18) = "mcFei-32"
      state_labels(19) = "mcFei-13"
      state_labels(20) = "mcFei-23"
      state_labels(21) = "mcFei-33"

      state_descriptors(13:21) = "-curl(Fe^-1) pseudo Nye tensor"
      
      state_labels(22) = "eps-11"
      state_labels(23) = "eps-22"
      state_labels(24) = "eps-33"
      state_labels(25) = "eps-13"
      state_labels(26) = "eps-23"
      state_labels(27) = "eps-12"

      state_descriptors(22:27) = "lattice strain"

      state_labels(28) = "R-11"
      state_labels(29) = "R-21"
      state_labels(30) = "R-31"
      state_labels(31) = "R-12"
      state_labels(32) = "R-22"
      state_labels(33) = "R-32"
      state_labels(34) = "R-13"
      state_labels(35) = "R-23"
      state_labels(36) = "R-33"

      state_descriptors(28:36) = "RU rotation"

      state_labels(37) = "max-slip-rate"
      state_descriptors(37) = "max-slip-rate"

      state_labels(38) = "max-sys-ID"
      state_descriptors(38) = "max-sys-ID"

      state_labels(39) = "number-active"
      state_descriptors(39) = "number-active"

      state_labels(40) = ""
      state_descriptors(40) = "eff. creep rate"

      state_labels(41) = ""
      state_descriptors(41) = "n_eff"

      do i = 1, max_slip_sys
            write(state_labels(i+41), 9000) i
            state_descriptors(i+41) = "integrated slip"
      end do

      do i = 1, max_uhard
            write(state_labels(i+max_slip_sys+41), 9020)
     & i
            state_descriptors(i+max_slip_sys+41) =
     &  "inter. hardening var."
      end do
c
      if( do_print ) then
        do i = 1, num_states 
          write(outi,9010) i, state_labels(i), state_descriptors(i)
        end do
        do_print = .false.   
      end if
c          
      return
 9000 format("slip-",i2.2)    
 9010 format(2x,i3,2x,a8,2x,a) 
 9020 format("hard.-",i2.2)      
      end
      

c
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_values                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 12/5/2014 (rhd)                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_values( itype, elem_states_output,
     &                                 nrow_states, num_states  )
c
c                       access some global data structures
c
      use elem_block_data, only: history_blocks, history_blk_list
      use main_data, only: elems_to_blocks
      implicit integer (a-z)
$add common.main
c
c                       parameters
c
      integer :: nrow_states, itype, num_states
#dbl      double precision :: elem_states_output(nrow_states,*)
#sgl      real  :: elem_states_output(nrow_states,*)
c
c                       locals
c
#dbl      double precision, 
#sgl      real,
     & allocatable :: history_dump(:,:,:), one_elem_states(:)
      integer :: relem, elnum, hist_size, blockno
      logical :: do_a_block
      double precision :: zero
      data zero / 0.0d00 /
c      
c           build CP states values output.
c
c              itype > 0 => this is the block number. do all elements
c                           in the block
c
c              itype < 0 => this is an element number. put state
c                           values into column 1 of results.
c 
      do_block = .true.
      if( itype. gt. 0 ) then
         do_a_block = .true.
         blockno = itype
      else
         do_a_block = .false.
         elnum = -itype
         blockno = elems_to_blocks(elnum,1)
      end if          
c
      local_debug = .false.      
      felem       = elblks(1,blockno)
      elem_type   = iprops(1,felem)
      mat_type    = iprops(25,felem)
      int_points  = iprops(6,felem)
      span        = elblks(0,blockno)
      hist_size   = history_blk_list(blockno)
      if( local_debug ) write(out,9050) blockno, felem, elem_type,         
     &         mat_type, int_points, span, hist_size
c
c           temporary block of history so it can be re-organized
c
      allocate( one_elem_states(nrow_states) )
      allocate( history_dump(hist_size,int_points,span) )
      history_dump = reshape( history_blocks(blockno)%ptr,
     &           (/hist_size,int_points,span/) )
c      
      if( do_a_block ) then    
        do relem = 1, span
           elnum = felem + relem - 1  ! absolute element number
           one_elem_states(1:nrow_states) = zero 
           call mm10_states_values_a
           elem_states_output(1:nrow_states,relem) = 
     &                one_elem_states(1:nrow_states)
        end do
      else
           relem = elnum + 1 - felem
           one_elem_states(1:nrow_states) = zero
           call mm10_states_values_a
           elem_states_output(1:nrow_states,1) =
     &                one_elem_states(1:nrow_states)
      end if  
c        
      deallocate( history_dump, one_elem_states )
c        
      return
c      
 9050 format(10x,"block, felem, etype, mtype:  ",4i7,
     &  /,10x,   "int_pts, span, hist_size:    ",3i7 )
c
      contains
c     ========      
     
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm10_states_values_a                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 5/30/2015 (tjt)            *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_states_values_a
c      
      use main_data, only: matprp, imatprp
      use crystal_data, only : c_array, angle_input, crystal_input,
     &                         data_offset
c
      implicit integer (a-z)
$add common.main
c
c                       locals
c
      integer :: matnum, nslip, s, e, d, z, b, w,
     &           cnum, osn, kk, sc
#dbl      double precision, 
#sgl      real,
     & dimension(3,3,3) :: gradFe
#dbl      double precision, 
#sgl      real,
     & dimension(3,3) :: nye
c
c
c           Goal for state output:
c           1:3                     euler angles      Cry 7:9
c           4:12                    Rp                Cry 10:18
c           13:21                   curl Fe^-1        ****
c                                   grad Fe^-1 is     37:63
c           22:27                   lattice strain    
c           28:36                   R    
c           37:39                   active slip systems   
c           40                      effective creep rate
c           41                      n_eff for creep
c           41+1:41+max_slip_sys    slip history      76:76+max_slip_sys-1
c           41+max_slip_sys+1:end    hardening      76:76+max_slip_sys-1
c          
c           Unfortunately we don't store ncrystals in the history, so
c           we need to access it in the material properties
c           If there are multiple crystals per GP then just look 
c           at the first first one
c
      matnum = iprops(38,felem)
      if( imatprp(104,matnum) .eq. 1 ) then
          cnum = imatprp(105,matnum)
      elseif( imatprp(104,matnum) .eq. 2 ) then
          osn  = data_offset(elnum)
          cnum = crystal_input(osn,1)
      else
          write(out,9000) elnum
          call die_gracefully
      end if
c
      nslip = c_array(cnum)%nslip
c
c           Start of the first crystal
      sc = 76+max_slip_sys - 1
c
c
c           Results 1-3: euler angles of the first crystal, 
c           first GP
c
      s = 7 + sc
      e = 9 + sc
      one_elem_states(1:3) = 
     &          ( history_dump(s:e,1,relem) )
c
c           Results 4-12: plastic rotation of the first crystal, 
c           first GP
c
      s = 10 + sc
      e = 18 + sc
      one_elem_states(4:12) = 
     &          ( history_dump(s:e,1,relem) )
c
c           Results 13-21: Nye tensor, averaged over gauss points
c
      s = 37
      e = 63
      gradFe = reshape( sum(history_dump(s:e,1:int_points,relem),2 ) /
     &            dble(int_points), (/3,3,3/) )
      nye(1:3,1:3) = zero
      do d = 1,3
          do z = 1,3
            do b = 1,3
                do w = 1,3
                  kk = (z-b)*(b-w)*(w-z)/2
                  nye(d,z) = nye(d,z) - dble(kk)*gradFe(d,b,w)
                end do
             end do
          end do
      end do
      one_elem_states(13:21) = reshape(nye, (/9/))
c
c           Results 22:27: the lattice strain, averaged over Gauss points
c    
      s = sc+24+1
      e = sc+30
      one_elem_states(22:27) = sum( 
     &      history_dump(s:e,1:int_points,relem), 2 ) / dble(int_points)
c
c           Results 28:36: R, first Gauss point
c    
      s = 64
      e = 72
      one_elem_states(28:36) = ( 
     &      history_dump(s:e,1,relem) )
c
c           Results 37:39: active slip systems, first Gauss point
c    
      s = sc + 30 + max_slip_sys + max_uhard + 6
      e = sc + 30 + max_slip_sys + max_uhard + 8
      one_elem_states(37:39) = ( 
     &      history_dump(s:e,1,relem) )
c
c           Results 40:41: effective creep, first Gauss point
c    
      s = sc + 30 + max_slip_sys + max_uhard + 11
      e = sc + 30 + max_slip_sys + max_uhard + 12
      one_elem_states(40:41) = ( 
     &      history_dump(s:e,1,relem) )
c
c           Results 39+1:39+max_slip_sys : the slip totals, 
c           padded with zeros as 
c           required  and averaged over Gauss points
c
      s = 76
      e = 76 + nslip - 1
      one_elem_states(41+1:41+max_slip_sys) = 
     &         sum( history_dump(s:e,1:int_points,relem),2 ) / 
     &         dble(int_points)

c
c           Results 13: tau_tilde of first crystal, avged
c
      s = sc + 30 + max_slip_sys + 1
      e = sc + 30 + max_slip_sys + max_uhard
      one_elem_states(41+1+max_slip_sys:41+max_slip_sys+max_uhard) = 
     &         sum( history_dump(s:e,1:int_points,relem),2 ) /
     &         dble(int_points)
c
      return     
c
 9000 format(/1x,
     &'>>>>> Error: invalid crystal detected. element: ',i7,
     & /,14x,'routine oustates_values_mm10_a',
     & /,14x,'job terminated....'/)
c
      end subroutine mm10_states_values_a
      end subroutine mm10_states_values
c
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
      subroutine mm10_setup_np1(Rur, dstrain, dt, T, step, elem, gp,
     &      np1)
      use mm10_defs
      implicit none
c
      double precision, dimension(9) :: Rur
      double precision, dimension(6) :: dstrain
      double precision :: dt, T
      integer :: step, elem, gp
      type(crystal_state) :: np1
c
      np1%R = reshape(Rur(1:9), (/3,3/))
      np1%D = dstrain(1:6)
      np1%temp = T
      np1%tinc = dt
      np1%step = step
      np1%elem = elem
      np1%gp = gp
c
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
      cc_props%miter = inc_props%miter
      cc_props%atol = inc_props%atol
      cc_props%atol1 = inc_props%atol1
      cc_props%rtol = inc_props%rtol
      cc_props%rtol1 = inc_props%rtol1
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
