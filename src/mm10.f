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
                  double precision, dimension(3,3) :: g
                  double precision :: ms(6,max_slip_sys),
     &                  qs(3,max_slip_sys), ns(3,max_slip_sys)
                  double precision, dimension(6,6) :: stiffness
                  integer :: angle_type, angle_convention, nslip, 
     &                  h_type
                  integer :: num_hard
                  logical :: real_tang
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
     &                        local_work%angle_convention(i),cc_props)
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
      double precision :: maxslip
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
c     1:3 is about the active slip systems: maximum slip rate,
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
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formR2                       *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *     Form R2                                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formR2(props, np1, n,vec1,vec2, 
     &           stress, tt, R2, gp)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: R2, h
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1,vec2
      integer :: gp
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
        call mm10_h_voche(props, np1, n, stress, tt, h)
      elseif (props%h_type .eq. 2) then ! MTS
        call mm10_h_mts(props, np1, n, stress, tt, h)
      elseif (props%h_type .eq. 3) then ! User
        call mm10_h_user(props, np1, n, stress, tt, h)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_h_mrr(props, np1, n, vec1, vec2,stress,tt,h,gp)
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
      R2 = tt - h
c Compute time rate of change of hardening parameters, store for
c possible prediction during next time/load step
      np1%tt_rate(1:props%num_hard) = (h(1:props%num_hard)
     &     - n%tau_tilde(1:props%num_hard))/np1.tinc
c      if(gp.ne.0) then
c        write(*,*) "R2 =", R2(6), " h= ", h(6), " ttrate= ", 
c     &  np1%tt_rate(6)
c      endif
c      
      return
      end subroutine
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formR2i                      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *     Form R2                                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formR2i(props, np1, n,ivec1,ivec2, 
     &           stress, tt, R2)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double complex, dimension(6) :: stress
      double complex, dimension(props%num_hard) :: h
      double complex, dimension(props%num_hard) :: R2
      double complex, dimension(props%num_hard) :: tt
      double complex, dimension(max_uhard) :: ivec1,ivec2
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
      elseif (props%h_type .eq. 2) then ! MTS
      elseif (props%h_type .eq. 3) then ! User
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_hi_mrr(props, np1, n, ivec1, ivec2, stress, tt, h)
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
      R2 = tt - h
c      
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ11                      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 03/31/15                    *
c     *                                                              *
c     *     Form the stress varying with stress part                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ11(props, np1, n, vec1, vec2,
     &      arr1, arr2, stress, tt, J11)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(6,6) :: J11
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1,vec2
      double precision, dimension(max_uhard,max_uhard) :: 
     &  arr1,arr2
c
      double precision :: mm10_rs
c
      double precision, dimension(6,props%nslip) :: symtqmat
      double precision, dimension(props%nslip) :: dgammadtau
      double precision, dimension(3) :: wp
      double precision, dimension(6,6) :: Iw
      double precision :: rs
      integer :: i
      logical :: debug
c
c
            debug = .false.
c
            if (debug) write (*,*) "In mm10"
c
      J11 = 0.0
      call mm10_symSWmat(stress, np1%qc, props%nslip, symtqmat)
c
c Generalization of CP model implementation for other slip rate
c equations, requiring other forms of d_gamma/d_tau
c Vector dgammadtau should be 1 x n_slip
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
        call mm10_dgdt_voche(props,np1, n, stress, tt, dgammadtau)
      elseif (props%h_type .eq. 2) then ! MTS
        call mm10_dgdt_mts(props, np1,n, stress, tt, dgammadtau)
      elseif (props%h_type .eq. 3) then ! User
        call mm10_dgdt_user(props,np1, n, stress, tt, dgammadtau)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_dgdt_mrr(props,np1, n, stress, tt, dgammadtau)
      else
        call mm10_unknown_hard_error(props)
      endif
c ******* END: Add new Constitutive Models into this block *********
c
      do i=1,props%nslip
        call DGER(6,6,dgammadtau(i),
     &      matmul(props%stiffness, np1%ms(1:6,i))
     &      + 2.0d0*symtqmat(1:6,i), 1, np1%ms(1:6,i),1,J11,6)
        if (debug) write (*,*) "J11", J11(1,1)
      end do
c
      call mm10_form_wp(props, np1, n, vec1, vec2, stress, tt, wp)
      call mm10_IW(wp, Iw)
      J11 = J11 + Iw
c
      do i=1,6
        J11(i,i) = J11(i,i) + 1.0d0
      end do
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ12                      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 03/31/15                    *
c     *                                                              *
c     *     Form the stress varying with hardening part              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ12(props, np1, n, vec1, vec2, arr1,
     &            arr2,stress, tt, J12)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(6,props%num_hard) :: J12
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1,vec2
      double precision, dimension(max_uhard,max_uhard) ::
     &     arr1,arr2
c
      double precision, dimension(6) :: symTW
      double precision, dimension(6,props%nslip) :: symtqmat
      double precision, dimension(6) :: tempmat
      double precision, dimension(props%nslip,props%num_hard)
     & :: dgammadtt
      double precision :: dgam, mm10_slipinc
      integer :: i
      logical :: debug
c
c
            debug = .false.
c
      J12 = 0.0d0
      call mm10_symSWmat(stress, np1%qc, props%nslip, symtqmat)
c        if (debug) write (*,*) "symtqmat", symtqmat(1:6,1:props%nslip)
c
c Generalization of CP model implementation for other slip rate
c equations, requiring other forms of d_gamma/d_hardening
c Vector dgammadtt should be n_slip x n_hardening, which is the
c derivative of slip rate alpha with respect to hardening variable
c beta.
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
        call mm10_dgdh_voche(props,np1, n, stress, tt, dgammadtt)
      elseif (props%h_type .eq. 2) then ! MTS
        call mm10_dgdh_mts(props, np1,n, stress, tt, dgammadtt)
      elseif (props%h_type .eq. 3) then ! User
        call mm10_dgdh_user(props,np1, n, stress, tt, dgammadtt)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_dgdh_mrr(props,np1, n, stress, tt, dgammadtt)
      else
        call mm10_unknown_hard_error(props)
      endif
c      if (debug) write(*,*) "dgammadtt", dgammadtt(1:12,1:12)
c      if (debug) write(*,*) "np1%ms(1:6,i)", np1%ms(1:6,1:props%nslip)
c      if (debug) write(*,*) "np1%qc(1:3,i)", np1%qc(1:3,1:props%nslip)
c ******* END: Add new Constitutive Models into this block *********
c
      do i=1,props%nslip
        tempmat(1:6) = matmul(props%stiffness, np1%ms(1:6,i))
     &      + 2.0d0*symtqmat(1:6,i)
        if (debug) write (*,*) "tempmat", tempmat(1:6)
        call DGER(6,props%num_hard,1.d0,
     &      tempmat(1:6), 1, dgammadtt(i,1:props%num_hard),1,
     &      J12(1:6,1:props%num_hard),6)
      end do !props%nslip,
c
      return
      end subroutine

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ21                      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/27/13                    *
c     *                                                              *
c     *     Form the hardening varying with stress part              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ21(props, np1, n, vec1, vec2,
     &     arr1, arr2, stress, tt, J21)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard,6) :: J21
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1,vec2
      double precision, dimension(max_uhard,max_uhard) :: 
     &    arr1,arr2
c
      double precision, dimension(props%num_hard,6) :: estr
c
      J21 = 0.0d0
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
        call mm10_estress_voche(props, np1, n, stress, tt, estr)
      elseif (props%h_type .eq. 2) then ! MTS
        call mm10_estress_mts(props, np1, n, stress, tt, estr)
      elseif (props%h_type .eq. 3) then ! User
        call mm10_estress_user(props, np1, n, stress, tt, estr)
      elseif (props%h_type .eq. 7) then ! User
        call mm10_estress_mrr(props, np1, n, vec1, vec2, arr1, arr2,
     &            stress, tt, estr)
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
      J21 = -estr
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ22                      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/27/13                    *
c     *                                                              *
c     *     Form the hardening varying with hardening part           *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ22(props, np1, n, vec1, vec2,
     &     arr1, arr2, stress, tt, J22)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard
     &         ,props%num_hard) :: J22, etau
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1,vec2
      double precision, dimension(max_uhard,max_uhard) :: 
     &   arr1,arr2
c
c
      J22 = 0.0d0
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
        call mm10_ehard_voche(props, np1, n, stress, tt, etau)
      elseif (props%h_type .eq. 2) then ! MTS
        call mm10_ehard_mts(props, np1, n, stress, tt, etau)
      elseif (props%h_type .eq. 3) then ! User
        call mm10_ehard_user(props, np1, n, stress, tt, etau)
      elseif (props%h_type .eq. 7) then ! User
        call mm10_ehard_mrr(props, np1, n, vec1, vec2, arr1, arr2,
     &          stress, tt, etau)
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
      J22 = etau
c
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formvecs                     *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 05/30/15                    *
c     *                                                              *
c     *     Form intermediate vectors which are repeatedly used by   *
c     *     other constitutive routines (e.g. precompute slip rates) *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formvecs(props, np1, n, stress, tt, vec1, vec2)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1,vec2
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
      elseif (props%h_type .eq. 2) then ! MTS
      elseif (props%h_type .eq. 3) then ! User
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_v_mrr(props, np1, n, stress, tt, 
     &   vec1, vec2)
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
c
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formarrs                     *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 05/30/15                    *
c     *                                                              *
c     *     Form intermediate arrays which are repeatedly used by    *
c     *     other constitutive routines                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formarrs(props, np1, n, stress, tt, vec1, vec2,
     &          arr1, arr2, both)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1,vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1,arr2
      integer both
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
      elseif (props%h_type .eq. 2) then ! MTS
      elseif (props%h_type .eq. 3) then ! User
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_a_mrr(props, np1, n, stress, tt, vec1, vec2,
     &          arr1, arr2, both)
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
c
      end subroutine

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ11i                     *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 03/31/15                    *
c     *                                                              *
c     *     Form the stress varying with stress part                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ11i(props, np1, n, ivec1, ivec2,
     &      stress, tt, J11)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, zeroA
      double precision, dimension(6,6) :: J11
      double precision, dimension(props%num_hard) :: tt, zeroB
      double complex, dimension(max_uhard) :: ivec1,ivec2
c
      integer :: i, k
      logical :: debug
c
      double complex, dimension(6) :: A, Ri
      double complex, dimension(props%num_hard) :: B
c
      double precision :: h, hi
      double complex :: i1
c
c
            debug = .false.
c
            if (debug) write (*,*) "In mm10"
      h = 1.0d-12
      i1 = (0.d0, 1.d0)
c
      J11 = 0.0d0
      zeroA = 0.d0
      zeroB = 0.d0
            do k = 1,6
                A = dcmplx (stress, zeroA)
                B = dcmplx (tt, zeroB)
                if(stress(k).eq.0.d0) then
                hi = h
                else
                hi = h*dabs(stress(k))
                endif
                A(k) = A(k) + i1*hi
                call mm10_formvecsi(props,np1,n,A,B,ivec1,ivec2)
                call mm10_formR1i(props,np1,n,ivec1,ivec2,A,B,Ri)
                J11(1:6,k) = 1.d0/hi*aimag(Ri)
            enddo
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ12i                     *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 03/31/15                    *
c     *                                                              *
c     *     Form the stress varying with hardening part              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ12i(props, np1, n, ivec1, ivec2,
     &            stress, tt, J12)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, zeroA
      double precision, dimension(6,props%num_hard) :: J12
      double precision, dimension(props%num_hard) :: tt, zeroB
      double complex, dimension(max_uhard) :: ivec1,ivec2
c
      integer :: i, k
      logical :: debug
c
      double complex, dimension(6) :: Ri, A
      double complex, dimension(props%num_hard) :: B
c
      double precision :: h, hi
      double complex :: i1
c
c
            debug = .false.
c
            if (debug) write (*,*) "In mm10"
      h = 1.0d-12
      i1 = (0.d0, 1.d0)
c
      J12 = 0.0d0
      zeroA = 0.d0
      zeroB = 0.d0
            do k = 1,props%num_hard
                A = dcmplx (stress, zeroA)
                B = dcmplx (tt, zeroB)
                if(tt(k).eq.0.d0) then
                hi = h
                else
                hi = h*dabs(tt(k))
                endif
                B(k) = B(k) + i1*hi
                call mm10_formvecsi(props,np1,n,A,B,ivec1,ivec2)
                call mm10_formR1i(props,np1,n,ivec1,ivec2,A,B,Ri)
                J12(1:6,k) = 1.d0/hi*aimag(Ri)
            enddo
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ21i                     *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/27/13                    *
c     *                                                              *
c     *     Form the hardening varying with stress part              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ21i(props, np1, n, ivec1, ivec2,
     &     stress, tt, J21)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, zeroA
      double precision, dimension(props%num_hard,6) :: J21
      double precision, dimension(props%num_hard) :: tt, zeroB
      double complex, dimension(max_uhard) :: ivec1,ivec2
      integer :: i, k
      logical :: debug
c
      double complex, dimension(6) :: A
      double complex, dimension(props%num_hard) :: B, Ri
c
      double precision :: h, hi
      double complex :: i1
c
c
            debug = .false.
c
            if (debug) write (*,*) "In mm10"
      h = 1.0d-12
      i1 = (0.d0, 1.d0)
c
      J21 = 0.0d0
      zeroA = 0.d0
      zeroB = 0.d0
            do k = 1,6
                A = dcmplx (stress, zeroA)
                B = dcmplx (tt, zeroB)
                if(stress(k).eq.0.d0) then
                hi = h
                else
                hi = h*dabs(stress(k))
                endif
                A(k) = A(k) + i1*hi
                call mm10_formvecsi(props,np1,n,A,B,ivec1,ivec2)
                call mm10_formR2i(props,np1,n,ivec1,ivec2,A,B,Ri)
                J21(1:props%num_hard,k) = 1.d0/hi*aimag(Ri)
            enddo
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ22i                     *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/27/13                    *
c     *                                                              *
c     *     Form the hardening varying with hardening part           *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ22i(props, np1, n, ivec1, ivec2,
     &     stress, tt, J22)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, zeroA
      double precision, 
     &   dimension(props%num_hard,props%num_hard) :: J22
      double precision, dimension(props%num_hard) :: tt, zeroB
      double complex, dimension(max_uhard) :: ivec1,ivec2
c
      integer :: i, k
      logical :: debug
c
      double complex, dimension(6) :: A
      double complex, dimension(props%num_hard) :: Ri, B
c
      double precision :: h, hi
      double complex :: i1
c
c
            debug = .false.
c
            if (debug) write (*,*) "In mm10"
      h = 1.0d-12
      i1 = (0.d0, 1.d0)
c
      J22 = 0.0d0
      zeroA = 0.d0
      zeroB = 0.d0
            do k = 1,props%num_hard
                A = dcmplx (stress, zeroA)
                B = dcmplx (tt, zeroB)
                if(tt(k).eq.0.d0) then
                hi = h
                else
                hi = h*dabs(tt(k))
                endif
                B(k) = B(k) + i1*hi
                call mm10_formvecsi(props,np1,n,A,B,ivec1,ivec2)
                call mm10_formR2i(props,np1,n,ivec1,ivec2,A,B,Ri)
                J22(1:props%num_hard,k) = 1.d0/hi*aimag(Ri)
            enddo
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formvecsi                    *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 05/30/15                    *
c     *                                                              *
c     *     Form intermediate vectors which are repeatedly used by   *
c     *     other constitutive routines (e.g. precompute slip rates) *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formvecsi(props, np1, n, stress, tt, ivec1, ivec2)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double complex, dimension(6) :: stress
      double complex, dimension(props%num_hard) :: tt
      double complex, dimension(max_uhard) :: ivec1,ivec2
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
      elseif (props%h_type .eq. 2) then ! MTS
      elseif (props%h_type .eq. 3) then ! User
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_vi_mrr(props, np1, n, stress, tt, 
     &   ivec1, ivec2)
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
c
      end subroutine
c
c *****************************************************************************
c *                                                                           *
c *         Small constitutive functions                                      *
c *                                                                           *
c *****************************************************************************
c
c
c           Form d_bar_p
      subroutine mm10_form_dbarp(props, np1, n, vec1, vec2,
     &     stress, tt, dbar)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(6) :: dbar
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1,vec2
c
      integer :: i
c
      double precision :: slipinc
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
      dbar = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
        dbar = dbar + slipinc*np1%ms(1:6,i)
      end do
      elseif (props%h_type .eq. 2) then ! MTS
      dbar = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
        dbar = dbar + slipinc*np1%ms(1:6,i)
      end do
      elseif (props%h_type .eq. 3) then ! User
      dbar = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc_user(props, np1, n, stress, tt, i, slipinc)
        dbar = dbar + slipinc*np1%ms(1:6,i)
      end do
      elseif (props%h_type .eq. 7) then ! MRR
c      dbar = 0.0
c      do i=1,props%nslip
c c        call mm10_slipinc_mrr(props, np1, n, stress, tt, i, slipinc)
c        slipinc = vec1(i)
c        dbar = dbar + slipinc*np1%ms(1:6,i)
c      end do
        dbar(1:6) = matmul(np1%ms(1:6,1:props%nslip),
     &                     vec1(1:props%nslip))
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
      return
      end subroutine
c
c
c           Form w_bar_p
      subroutine mm10_form_wbarp(props, np1, n, vec1, vec2,
     &   stress, tt, wbar)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(3) :: wbar
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1,vec2
c
      integer :: i
c
      double precision :: slipinc
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
      wbar = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
        wbar = wbar + slipinc*np1%qs(1:3,i)
      end do
      elseif (props%h_type .eq. 2) then ! MTS
      wbar = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
        wbar = wbar + slipinc*np1%qs(1:3,i)
      end do
      elseif (props%h_type .eq. 3) then ! User
      wbar = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc_user(props, np1, n, stress, tt, i, slipinc)
        wbar = wbar + slipinc*np1%qs(1:3,i)
      end do
      elseif (props%h_type .eq. 7) then ! MRR
c      wbar = 0.0
c      do i=1,props%nslip
c c        call mm10_slipinc_mrr(props, np1, n, stress, tt, i, slipinc)
c        slipinc = vec1(i)
c        wbar = wbar + slipinc*np1%qs(1:3,i)
c      end do
        wbar(1:3) = matmul(np1%qs(1:3,1:props%nslip),
     &                     vec1(1:props%nslip))
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
      return
      end subroutine

c           Form w_p (in the current configuration)
      subroutine mm10_form_wp(props, np1, n, vec1, vec2, stress, tt, w)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(3) :: w
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1,vec2
c
      integer :: i
c
      double precision :: slipinc
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
      w = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
        w = w + slipinc*np1%qc(1:3,i)
      end do
      elseif (props%h_type .eq. 2) then ! MTS
      w = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
        w = w + slipinc*np1%qc(1:3,i)
      end do
      elseif (props%h_type .eq. 3) then ! User
      w = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc_user(props, np1, n, stress, tt, i, slipinc)
        w = w + slipinc*np1%qc(1:3,i)
      end do
      elseif (props%h_type .eq. 7) then ! MRR
      w = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc_mrr(props, np1, n, stress, tt, i, slipinc)
        w = w + slipinc*np1%qc(1:3,i)
      end do
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
      return
      end subroutine
c
c           Form d_bar_p
      subroutine mm10_form_dbarpi(props, np1, n, ivec1, ivec2,
     &     stress, tt, dbar)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double complex, dimension(6) :: stress
      double complex, dimension(6) :: dbar
      double complex, dimension(props%num_hard) :: tt
      double complex, dimension(max_uhard) :: ivec1,ivec2
c
      integer :: i
c
      double complex :: slipinc
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
      elseif (props%h_type .eq. 2) then ! MTS
      elseif (props%h_type .eq. 3) then ! User
      elseif (props%h_type .eq. 7) then ! MRR
      dbar = (0.d0,0.d0)
      do i=1,props%nslip
c        call mm10_slipinc_mrr(props, np1, n, stress, tt, i, slipinc)
        slipinc = ivec1(i)
        dbar = dbar + slipinc*np1%ms(1:6,i)
      end do
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
      return
      end subroutine
c
c           Form w_p (in the current configuration)
      subroutine mm10_form_wpi(props, np1, n, ivec1, ivec2, 
     &     stress, tt, w)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double complex, dimension(6) :: stress
      double complex, dimension(3) :: w
      double complex, dimension(props%num_hard) :: tt
      double complex, dimension(max_uhard) :: ivec1,ivec2
c
      integer :: i
c
      double complex :: slipinc
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
      elseif (props%h_type .eq. 2) then ! MTS
      elseif (props%h_type .eq. 3) then ! User
      elseif (props%h_type .eq. 7) then ! MRR
      w = (0.d0,0.d0)
      do i=1,props%nslip
        call mm10_slipinci_mrr(props, np1, n, stress, tt, i, slipinc)
        w = w + slipinc*np1%qc(1:3,i)
      end do
      else
        call mm10_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
      return
      end subroutine
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
      subroutine mm10_init_cc_props(inc_props, atype, aconv, cc_props)
      use mm10_defs
      implicit integer (a-z)
$add include_sig_up
      integer :: atype, aconv
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
     &   arr2, ivec1, ivec2, fail,gp)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      logical :: fail
c
      type(crystal_state) :: curr
      double precision, dimension(6) :: stress, ostress
      double precision :: tt(props%num_hard), ott(props%num_hard)
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
      double complex, dimension(max_uhard) :: ivec1, ivec2
      integer :: cuts, mcuts, i, gp
      double precision :: frac, step, mult
      integer, dimension(10) :: faili
      double precision, dimension(10) :: failr
      logical :: debug
      parameter(mcuts = 8)
      parameter(mult = 0.5d0)
      parameter(debug = .false.)
c
      frac = 0.0d0
      step = 1.0d0
      cuts = 0
c
      stress = n%stress
      tt = n%tau_tilde
      ostress = stress
      ott = tt
c
      call mm10_setup(props, np1, n)
      fail = .false.

      do while (frac .lt. 1.0d0)
        call mm10_setup_np1(reshape(np1%R, (/9/)), np1%D*(step+frac),
     &      np1%tinc*(step+frac), (np1%temp-n%temp)*(step+frac)+n%temp,
     &      np1%step, np1%elem, np1%gp, curr)
        call mm10_setup(props, curr, n)
        tt = n%tau_tilde ! some subroutines modify n%tau_tilde in setup
        call mm10_solve(props, curr, n, vec1, vec2, arr1, arr2,
     &       ivec1, ivec2, stress, tt, fail, faili, failr, gp,
     &       np1%tinc*step)
        if (fail) then
          if (debug) write(*,*) "Adapting"
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
      end if
c      
      np1%stress = stress
      np1%tau_tilde = tt
      np1%tt_rate = curr%tt_rate ! store rates of hardening too
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
     &      nlsx, nRs, c, red, dt, cos_ang
      integer :: iter, miter, info, ls, mls, mmin, gp, gpp
      integer, dimension(6+props%num_hard) :: ipiv
      logical :: debug, gpall, solver, strategy
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
      parameter(atol = 1.0d-5)
      parameter(atol1 = 1.0d-5)
      parameter(rtol = 5.0d-5)
      parameter(rtol1 = 1.0d-5)
      parameter(miter = 30)
      parameter(c = 1.0d-4)
      parameter(red = 0.5d0)
      parameter(mls = 10)
      parameter(mmin = 1)
c       Trust region parameters
      xtol = 1.0d-2
      maxfev = 3*miter
      lr=((6+props%num_hard)*(6+props%num_hard+1))/2
c      Debug flags for printing iteration behavior
      debug = .false.
      gpall = .false. ! true to print iteration norms for all Gauss points
      gpp = 0 ! set one particular G.P. to print
c      Solver flags
      solver = .true. ! true for Mark's N.R. routine, false for trust-region
      strategy = .true. ! true for geometric l.s., false for cubic l.s.
c
      if (debug) write(*,*) "Entering solution routine"
c
      x(1:6) = stress
      x(7:props%num_hard+6) = tt
c
      if(debug .and.(gpall.or.(gp.eq.gpp))) then ! print statement for debugging
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
      d1(1:6) = np1%D(1:6)/dsqrt(dot_product(np1%D,np1%D))
      d2(1:6) = n%D(1:6)/dsqrt(dot_product(n%D,n%D))
      cos_ang = dmax1(dot_product(d1,d2),0.d0)
      !write(*,*) cos_ang, n%tt_rate(2)
      x2 = x(7:6+props%num_hard) + cos_ang*
     &        n%tt_rate(1:props%num_hard)*dt

      if(debug)
     & write(props%out,*)" Stress prediction module, G.P.=", gp
      if(debug .and.(gpall.or.(gp.eq.gpp))) then ! print statement for debugging
      write(props%out,*)" Extrapol tt6=",
     &  x2(6), " Previous tt6=", x(6+6)
      endif
      call mm10_formvecs(props, np1, n, x1, x2, vec1, vec2)
      call mm10_formR1(props, np1, n, vec1, vec2, x1,x2, R1, gp)
      nR1 = dsqrt(dot_product(R1,R1))
      inR1 = nR1
c       Newton-Raphson loop
      if (debug) write(*,'("Iter ",i3," norm ",E10.3)') iter, nR1
      do while ((nR1 .gt. atol1) .and. (nR1/inR1 .gt. rtol1))
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
c        if (debug) write (*,*) "J11", J11(1:6,1:6)
c           Increment
        dx1 = R1
        call DGESV(6, 1, -J11, 6, ipiv, dx1, 6, info)
      if(debug .and.(gpall.or.(gp.eq.gpp))) then ! print statement for debugging
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
          nR1 = dsqrt(dot_product(R1,R1))
          nRs = 0.5d0*dot_product(R1,R1)
c             Line search convergence test
          if ((nRs .le. nlsx) .or. (ls .gt. mls)) then
            x1 = xnew1
      if(debug .and.(gpall.or.(gp.eq.gpp))) then ! print statement for debugging
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
        if (debug) write(*,'("Iter ",i3," norm ",E10.3," ls",F10.3)')
     &            iter, nR1, alpha
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
      if(debug .and.(gpall.or.(gp.eq.gpp))) then ! print statement for debugging
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
      if(debug) ! print statement for debugging
     & write(props%out,*)" Material update module, G.P.=", gp
      if(debug .and.(gpall.or.(gp.eq.gpp))) then ! print statement for debugging
      write(props%out,*)" Guess syy=",
     &  x(2), " Guess tt=", x(6+6)
      endif
c
      if(solver) then ! Newton Raphson with line search
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
      if (debug) write(*,'("Iter ",i3," norm ",E10.3)') iter, nR
      do while (((nR .gt. atol) .and. (nR/inR .gt. rtol)) .or. 
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
      if(debug .and.(gpall.or.(gp.eq.gpp))) then ! print statement for debugging
      write(props%out,*)" Iter=",
     &  iter, " dx1=", dx(2), " dx2=", dx(6+6)
      endif
c
        if(strategy) then ! geometric line search
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
            x = xnew
      if(debug .and.(gpall.or.(gp.eq.gpp))) then ! print statement for debugging
      write(props%out,*)" L.S. converged, Iter=",
     &  iter, " syy=", x(2), " tt(5)=", x(6+6), " nR=", nR
      write(props%out,*)" alpha=",
     &  alpha
      endif
            exit
          else
            alpha = red*alpha
            ls = ls + 1
          end if
        end do
c
        else ! cubic line search
        ls1 = 0.5d0*dot_product(R,R)
        g = matmul(transpose(J),R)
        stepmx = 1.d0
        xtol = 1.d-3
        wa1 = 1.d0 ! scalex
        i = 1 ! priter
          call mm10_nwclsh(6+props%num_hard,x,ls1,dx,g,stepmx,xtol,
     *                  wa1,xnew,R,nR,wa2,info,ls,i,iter,
     &                  props, np1, n, vec1, vec2,gp)
            if(info.eq.0) then
            x = xnew
            nR = dsqrt(2.0d0*nR)
            alpha = 1.d0
            endif
        endif
c           Increment and check for failure
        iter = iter + 1
        if (debug) write(*,'("Iter ",i3," norm ",E10.3," ls",F10.3)')
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
      else ! trust region dogleg solver
c         Get typical stress and tt
c           alpha = mm10_enorm(6,stress)
c           if(alpha.eq.0.d0) then
c             diag(1:6) = 1.d0
c           else
c             do i = 1,6
c             diag(i) = 1.d5/dabs(x(i))
c             enddo
c           endif
c           alpha = mm10_enorm(6,tt)
c           if(alpha.eq.0.d0) then
c             diag(7:6+props%num_hard) = 1.d0
c           else
c             diag(7:6+props%num_hard) = 10.d0
c     *           /dabs(tt(1:props%num_hard))
c           endif
         call mm10_hybrj(6+props%num_hard,x,R,J,6+props%num_hard,
     *             xtol,maxfev,diag,
     *             1,100.0d0,maxfev,info,nfev,njev,rout,lr,qtf,wa1,wa2,
     *                 wa3,wa4,atol,rtol,   nR,inR, gpall, debug,gpp,
     &      props, np1, n, vec1, vec2, arr1, arr2, ivec1, ivec2,gp)
c         Record data and reason for failure
        if (info.ne.1) then
          if(debug .and.(gpall.or.(gp.eq.gpp))) 
     &      write(*,*) "info=",info, "gp=",gp
          fail = .true.
          faili(1) = nfev
          faili(2) = maxfev
          faili(3) = info
          faili(4) = 2
          return
        end if
      endif
c
c       Output statistics from N-R algorithm
      if(debug .and.(gpall.or.(gp.eq.gpp))) then ! print statement for debugging
      write(props%out,*)" Material upd conv, iter=",
     &  iter
      write(props%out,*)" AbsNorm=",
     &  nR, " AbsTol=", atol
      write(props%out,*)" RelNorm=",
     &  nR/inR, " RelTol=", rtol
      write(props%out,*)" Guess syy=",
     &  stress(2), " actual syy=", x(2)
      write(props%out,*)" Guess tt6=",
     &  x2(6), " actual tt6=", x(6+6)
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
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_hybrj                        *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 7/03/15                     *
c     *                                                              *
c     *     Set of subroutines for trust region dogleg solver        *
c     *                                                              *
c     ****************************************************************
c
c 
      subroutine mm10_hybrj(n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,
     *             mode,factor,nprint,info,nfev,njev,r,lr,qtf,wa1,wa2,
     *                 wa3,wa4,atol,rtol,  fnorm,inR,gpall,debug,gpp,
     &      props, np1, np0, vec1, vec2, arr1, arr2, ivec1, ivec2,gp)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, np0
      double precision, dimension(max_uhard) :: vec1,vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1,arr2
      double complex, dimension(max_uhard) :: ivec1, ivec2
      integer n,ldfjac,maxfev,mode,nprint,info,nfev,njev,lr,gp,gpp
      logical debug,gpall
      double precision xtol,factor,atol,rtol
      double precision x(n),fvec(n),fjac(ldfjac,n),diag(n),r(lr),
     *                 qtf(n),wa1(n),wa2(n),wa3(n),wa4(n)
c     **********
c
c     subroutine hybrj
c
c     the purpose of hybrj is to find a zero of a system of
c     n nonlinear functions in n variables by a modification
c     of the powell hybrid method. the user must provide a
c     subroutine which calculates the functions and the jacobian.
c
c     the subroutine statement is
c
c       subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,
c                        mode,factor,nprint,info,nfev,njev,r,lr,qtf,
c                        wa1,wa2,wa3,wa4)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions and the jacobian. fcn must
c         be declared in an external statement in the user
c         calling program, and should be written as follows.
c
c         subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
c         integer n,ldfjac,iflag
c         double precision x(n),fvec(n),fjac(ldfjac,n)
c         ----------
c         if iflag = 1 calculate the functions at x and
c         return this vector in fvec. do not alter fjac.
c         if iflag = 2 calculate the jacobian at x and
c         return this matrix in fjac. do not alter fvec.
c         ---------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of hybrj.
c         in this case set iflag to a negative integer.
c
c       n is a positive integer input variable set to the number
c         of functions and variables.
c
c       x is an array of length n. on input x must contain
c         an initial estimate of the solution vector. on output x
c         contains the final estimate of the solution vector.
c
c       fvec is an output array of length n which contains
c         the functions evaluated at the output x.
c
c       fjac is an output n by n array which contains the
c         orthogonal matrix q produced by the qr factorization
c         of the final approximate jacobian.
c
c       ldfjac is a positive integer input variable not less than n
c         which specifies the leading dimension of the array fjac.
c
c       xtol is a nonnegative input variable. termination
c         occurs when the relative error between two consecutive
c         iterates is at most xtol.
c
c       maxfev is a positive integer input variable. termination
c         occurs when the number of calls to fcn with iflag = 1
c         has reached maxfev.
c
c       diag is an array of length n. if mode = 1 (see
c         below), diag is internally set. if mode = 2, diag
c         must contain positive entries that serve as
c         multiplicative scale factors for the variables.
c
c       mode is an integer input variable. if mode = 1, the
c         variables will be scaled internally. if mode = 2,
c         the scaling is specified by the input diag. other
c         values of mode are equivalent to mode = 1.
c
c       factor is a positive input variable used in determining the
c         initial step bound. this bound is set to the product of
c         factor and the euclidean norm of diag*x if nonzero, or else
c         to factor itself. in most cases factor should lie in the
c         interval (.1,100.). 100. is a generally recommended value.
c
c       nprint is an integer input variable that enables controlled
c         printing of iterates if it is positive. in this case,
c         fcn is called with iflag = 0 at the beginning of the first
c         iteration and every nprint iterations thereafter and
c         immediately prior to return, with x and fvec available
c         for printing. fvec and fjac should not be altered.
c         if nprint is not positive, no special calls of fcn
c         with iflag = 0 are made.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set as follows.
c
c         info = 0   improper input parameters.
c
c         info = 1   relative error between two consecutive iterates
c                    is at most xtol.
c
c         info = 2   number of calls to fcn with iflag = 1 has
c                    reached maxfev.
c
c         info = 3   xtol is too small. no further improvement in
c                    the approximate solution x is possible.
c
c         info = 4   iteration is not making good progress, as
c                    measured by the improvement from the last
c                    five jacobian evaluations.
c
c         info = 5   iteration is not making good progress, as
c                    measured by the improvement from the last
c                    ten iterations.
c
c       nfev is an integer output variable set to the number of
c         calls to fcn with iflag = 1.
c
c       njev is an integer output variable set to the number of
c         calls to fcn with iflag = 2.
c
c       r is an output array of length lr which contains the
c         upper triangular matrix produced by the qr factorization
c         of the final approximate jacobian, stored rowwise.
c
c       lr is a positive integer input variable not less than
c         (n*(n+1))/2.
c
c       qtf is an output array of length n which contains
c         the vector (q transpose)*fvec.
c
c       wa1, wa2, wa3, and wa4 are work arrays of length n.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... dogleg,dpmpar,enorm,
c                            qform,qrfac,r1mpyq,r1updt
c
c       fortran-supplied ... dabs,dmax1,dmin1,mod
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,iflag,iter,j,jm1,l,ncfail,ncsuc,nslow1,nslow2
      integer iwa(1)
      logical jeval,sing
      double precision actred,delta,epsmch,fnorm,fnorm1,one,pnorm,
     *                 prered,p1,p5,p001,p0001,ratio,sum,temp,xnorm,
     *                 zero,inR,temp2
      double precision mm10_dpmpar,mm10_enorm
      data one,p1,p5,p001,p0001,zero
     *     /1.0d0,1.0d-1,5.0d-1,1.0d-3,1.0d-4,0.0d0/
c
c     epsmch is the machine precision.
c
      epsmch = mm10_dpmpar(1)
c
      info = 0
      iflag = 0
      nfev = 0
      njev = 0
c
c     check the input parameters for errors.
c
      if (n .le. 0 .or. ldfjac .lt. n .or. xtol .lt. zero
     *    .or. maxfev .le. 0 .or. factor .le. zero
     *    .or. lr .lt. (n*(n + 1))/2) then
       write(*,*) n,ldfjac,xtol,maxfev,factor,lr,mode
       go to 300
      endif
      if (mode .ne. 2) go to 20
      do 10 j = 1, n
         if (diag(j) .le. zero) go to 300
   10    continue
   20 continue
c
c     evaluate the function at the starting point
c     and calculate its norm.
c
      iflag = 1
c      call fcn(n,x,fvec,fjac,ldfjac,iflag)
      call mm10_formR(props, np1, np0, vec1, vec2, x(1:6),
     & x(7:6+props%num_hard), fvec, gp)
      nfev = 1
      if (iflag .lt. 0) go to 300
      fnorm = mm10_enorm(n,fvec)
      inR = fnorm
            if (debug .and.(gpall.or.(gp.eq.gpp))) then
               temp = mm10_enorm(6,fvec(1:6))
               temp2 = mm10_enorm(props%num_hard,
     *                 fvec(7:6+props%num_hard))
               write(*,*) "iter = 0",
     *                    "AbsNorm(R1) = ", temp,
     *                    "AbsNorm(R2) = ", temp2
            endif
c
c     initialize iteration counter and monitors.
c
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
c
c     beginning of the outer loop.
c
   30 continue
         jeval = .true.
c
c        calculate the jacobian matrix.
c
         iflag = 2
c         call fcn(n,x,fvec,fjac,ldfjac,iflag)
        if(props%real_tang) then ! tangent matrix implemented
        call mm10_formJ(props, np1, np0, vec1, vec2, arr1, arr2, x(1:6),
     & x(7:6+props%num_hard), fjac)
        else
        call mm10_formJi(props, np1, np0, ivec1, ivec2, x(1:6),
     & x(7:6+props%num_hard), fjac)
        endif
         njev = njev + 1
         if (iflag .lt. 0) go to 300
c
c        compute the qr factorization of the jacobian.
c
         call mm10_qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1
     &    ,wa2,wa3)
c
c        on the first iteration and if mode is 1, scale according
c        to the norms of the columns of the initial jacobian.
c
         if (iter .ne. 1) go to 70
         if (mode .eq. 2) go to 50
         do 40 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) .eq. zero) diag(j) = one
   40       continue
   50    continue
c
c        on the first iteration, calculate the norm of the scaled x
c        and initialize the step bound delta.
c
         do 60 j = 1, n
            wa3(j) = diag(j)*x(j)
   60       continue
         xnorm = mm10_enorm(n,wa3)
         delta = factor*xnorm
         if (delta .eq. zero) delta = factor
   70    continue
c
c        form (q transpose)*fvec and store in qtf.
c
         do 80 i = 1, n
            qtf(i) = fvec(i)
   80       continue
         do 120 j = 1, n
            if (fjac(j,j) .eq. zero) go to 110
            sum = zero
            do 90 i = j, n
               sum = sum + fjac(i,j)*qtf(i)
   90          continue
            temp = -sum/fjac(j,j)
            do 100 i = j, n
               qtf(i) = qtf(i) + fjac(i,j)*temp
  100          continue
  110       continue
  120       continue
c
c        copy the triangular factor of the qr factorization into r.
c
         sing = .false.
         do 150 j = 1, n
            l = j
            jm1 = j - 1
            if (jm1 .lt. 1) go to 140
            do 130 i = 1, jm1
               r(l) = fjac(i,j)
               l = l + n - i
  130          continue
  140       continue
            r(l) = wa1(j)
            if (wa1(j) .eq. zero) sing = .true.
  150       continue
c
c        accumulate the orthogonal factor in fjac.
c
         call mm10_qform(n,n,fjac,ldfjac,wa1)
c
c        rescale if necessary.
c
         if (mode .eq. 2) go to 170
         do 160 j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
  160       continue
  170    continue
c
c        beginning of the inner loop.
c
  180    continue
c
c           if requested, call fcn to enable printing of iterates.
c
            if (nprint .le. 0) go to 190
            iflag = 0
            if (debug .and.(gpall.or.(gp.eq.gpp))) then
               temp = mm10_enorm(6,fvec(1:6))
               temp2 = mm10_enorm(props%num_hard,
     *                 fvec(7:6+props%num_hard))
               write(*,*) "iter = ", iter,
     *                    "AbsNorm(R1) = ", temp,
     *                    "AbsNorm(R2) = ", temp2
            endif
c     *         call fcn(n,x,fvec,fjac,ldfjac,iflag)
            if (iflag .lt. 0) go to 300
  190       continue
c
c           determine the direction p.
c
            call mm10_dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
c
c           store the direction p and x + p. calculate the norm of p.
c
            do 200 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  200          continue
            pnorm = mm10_enorm(n,wa3)
c
c           on the first iteration, adjust the initial step bound.
c
            if (iter .eq. 1) delta = dmin1(delta,pnorm)
c
c           evaluate the function at x + p and calculate its norm.
c
            iflag = 1
c            call fcn(n,wa2,wa4,fjac,ldfjac,iflag)
      call mm10_formR(props, np1, np0, vec1, vec2, wa2(1:6),
     & wa2(7:6+props%num_hard), wa4, gp)
            nfev = nfev + 1
            if (iflag .lt. 0) go to 300
            fnorm1 = mm10_enorm(n,wa4)
            if (debug .and.(gpall.or.(gp.eq.gpp))) then
               temp = mm10_enorm(6,wa4(1:6))
               temp2 = mm10_enorm(props%num_hard,
     *                 wa4(7:6+props%num_hard))
               write(*,*) "nslow1 = ", nslow1,
     *                    "AbsNorm(R1) = ", temp,
     *                    "AbsNorm(R2) = ", temp2
            endif
c
c           compute the scaled actual reduction.
c
            actred = -one
            if (fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
c
c           compute the scaled predicted reduction.
c
            l = 1
            do 220 i = 1, n
               sum = zero
               do 210 j = i, n
                  sum = sum + r(l)*wa1(j)
                  l = l + 1
  210             continue
               wa3(i) = qtf(i) + sum
  220          continue
            temp = mm10_enorm(n,wa3)
            prered = zero
            if (temp .lt. fnorm) prered = one - (temp/fnorm)**2
c
c           compute the ratio of the actual to the predicted
c           reduction.
c
            ratio = zero
            if (prered .gt. zero) ratio = actred/prered
c
c           update the step bound.
c
            if (ratio .ge. p1) go to 230
               ncsuc = 0
               ncfail = ncfail + 1
               delta = p5*delta
               go to 240
  230       continue
               ncfail = 0
               ncsuc = ncsuc + 1
               if (ratio .ge. p5 .or. ncsuc .gt. 1)
     *            delta = dmax1(delta,pnorm/p5)
               if (dabs(ratio-one) .le. p1) delta = pnorm/p5
  240       continue
c
c           test for successful iteration.
c
            if (ratio .lt. p0001) go to 260
c
c           successful iteration. update x, fvec, and their norms.
c
            do 250 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
               fvec(j) = wa4(j)
  250          continue
            xnorm = mm10_enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  260       continue
c
c           determine the progress of the iteration.
c
            nslow1 = nslow1 + 1
            if (actred .ge. p001) nslow1 = 0
            if (jeval) nslow2 = nslow2 + 1
            if (actred .ge. p1) nslow2 = 0
c
c           test for convergence.
c
            if (delta .le. xtol*xnorm .and. 
     *         ((fnorm .le. atol).or.(fnorm/inR .le. rtol))) info = 1
            if (debug .and.(gpall.or.(gp.eq.gpp))) then
               write(*,*) "delta = ", delta,
     *                    "xx = ", xtol*xnorm
            endif
            if (info .ne. 0) go to 300
c
c           tests for termination and stringent tolerances.
c
            if (nfev .ge. maxfev) info = 2
            if (p1*dmax1(p1*delta,pnorm) .le. epsmch*xnorm) info = 3
            if (nslow2 .eq. 5) info = 4
            if (nslow1 .eq. 10) info = 5
            if (info .ne. 0) go to 300
c
c           criterion for recalculating jacobian.
c
            if (ncfail .eq. 2) go to 290
c
c           calculate the rank one modification to the jacobian
c           and update qtf if necessary.
c
            do 280 j = 1, n
               sum = zero
               do 270 i = 1, n
                  sum = sum + fjac(i,j)*wa4(i)
  270             continue
               wa2(j) = (sum - wa3(j))/pnorm
               wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
               if (ratio .ge. p0001) qtf(j) = sum
  280          continue
c
c           compute the qr factorization of the updated jacobian.
c
            call mm10_r1updt(n,n,r,lr,wa1,wa2,wa3,sing)
            call mm10_r1mpyq(n,n,fjac,ldfjac,wa2,wa3)
            call mm10_r1mpyq(1,n,qtf,1,wa2,wa3)
c
c           end of the inner loop.
c
            jeval = .false.
            go to 180
  290    continue
c
c        end of the outer loop.
c
         go to 30
  300 continue
c
c     termination, either normal or user imposed.
c
      if (iflag .lt. 0) info = iflag
      iflag = 0
      if (debug .and.(gpall.or.(gp.eq.gpp))) then !call fcn(n,x,fvec,fjac,ldfjac,iflag)
           
               temp = mm10_enorm(6,fvec(1:6))
               temp2 = mm10_enorm(props%num_hard,
     *                 fvec(7:6+props%num_hard))
               write(*,*) "AbsNorm(R1) = ", temp,
     *                    "AbsNorm(R2) = ", temp2
      endif
      return
c
c     last card of subroutine hybrj.
c
      end subroutine
c
c
      double precision function mm10_dpmpar(i)
      integer i
c     **********
c
c     Function dpmpar
c
c     This function provides double precision machine parameters
c     when the appropriate set of data statements is activated (by
c     removing the c from column 1) and all other data statements are
c     rendered inactive. Most of the parameter values were obtained
c     from the corresponding Bell Laboratories Port Library function.
c
c     The function statement is
c
c       double precision function dpmpar(i)
c
c     where
c
c       i is an integer input variable set to 1, 2, or 3 which
c         selects the desired machine parameter. If the machine has
c         t base b digits and its smallest and largest exponents are
c         emin and emax, respectively, then these parameters are
c
c         dpmpar(1) = b**(1 - t), the machine precision,
c
c         dpmpar(2) = b**(emin - 1), the smallest magnitude,
c
c         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
c
c     Argonne National Laboratory. MINPACK Project. November 1996.
c     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'
c
c     **********
      integer mcheps(4)
      integer minmag(4)
      integer maxmag(4)
      double precision dmach(3)
      equivalence (dmach(1),mcheps(1))
      equivalence (dmach(2),minmag(1))
      equivalence (dmach(3),maxmag(1))
c
c     Machine constants for the IBM 360/370 series,
c     the Amdahl 470/V6, the ICL 2900, the Itel AS/6,
c     the Xerox Sigma 5/7/9 and the Sel systems 85/86.
c
c     data mcheps(1),mcheps(2) / z34100000, z00000000 /
c     data minmag(1),minmag(2) / z00100000, z00000000 /
c     data maxmag(1),maxmag(2) / z7fffffff, zffffffff /
c
c     Machine constants for the Honeywell 600/6000 series.
c
c     data mcheps(1),mcheps(2) / o606400000000, o000000000000 /
c     data minmag(1),minmag(2) / o402400000000, o000000000000 /
c     data maxmag(1),maxmag(2) / o376777777777, o777777777777 /
c
c     Machine constants for the CDC 6000/7000 series.
c
c     data mcheps(1) / 15614000000000000000b /
c     data mcheps(2) / 15010000000000000000b /
c
c     data minmag(1) / 00604000000000000000b /
c     data minmag(2) / 00000000000000000000b /
c
c     data maxmag(1) / 37767777777777777777b /
c     data maxmag(2) / 37167777777777777777b /
c
c     Machine constants for the PDP-10 (KA processor).
c
c     data mcheps(1),mcheps(2) / "114400000000, "000000000000 /
c     data minmag(1),minmag(2) / "033400000000, "000000000000 /
c     data maxmag(1),maxmag(2) / "377777777777, "344777777777 /
c
c     Machine constants for the PDP-10 (KI processor).
c
c     data mcheps(1),mcheps(2) / "104400000000, "000000000000 /
c     data minmag(1),minmag(2) / "000400000000, "000000000000 /
c     data maxmag(1),maxmag(2) / "377777777777, "377777777777 /
c
c     Machine constants for the PDP-11. 
c
c     data mcheps(1),mcheps(2) /   9472,      0 /
c     data mcheps(3),mcheps(4) /      0,      0 /
c
c     data minmag(1),minmag(2) /    128,      0 /
c     data minmag(3),minmag(4) /      0,      0 /
c
c     data maxmag(1),maxmag(2) /  32767,     -1 /
c     data maxmag(3),maxmag(4) /     -1,     -1 /
c
c     Machine constants for the Burroughs 6700/7700 systems.
c
c     data mcheps(1) / o1451000000000000 /
c     data mcheps(2) / o0000000000000000 /
c
c     data minmag(1) / o1771000000000000 /
c     data minmag(2) / o7770000000000000 /
c
c     data maxmag(1) / o0777777777777777 /
c     data maxmag(2) / o7777777777777777 /
c
c     Machine constants for the Burroughs 5700 system.
c
c     data mcheps(1) / o1451000000000000 /
c     data mcheps(2) / o0000000000000000 /
c
c     data minmag(1) / o1771000000000000 /
c     data minmag(2) / o0000000000000000 /
c
c     data maxmag(1) / o0777777777777777 /
c     data maxmag(2) / o0007777777777777 /
c
c     Machine constants for the Burroughs 1700 system.
c
c     data mcheps(1) / zcc6800000 /
c     data mcheps(2) / z000000000 /
c
c     data minmag(1) / zc00800000 /
c     data minmag(2) / z000000000 /
c
c     data maxmag(1) / zdffffffff /
c     data maxmag(2) / zfffffffff /
c
c     Machine constants for the Univac 1100 series.
c
c     data mcheps(1),mcheps(2) / o170640000000, o000000000000 /
c     data minmag(1),minmag(2) / o000040000000, o000000000000 /
c     data maxmag(1),maxmag(2) / o377777777777, o777777777777 /
c
c     Machine constants for the Data General Eclipse S/200.
c
c     Note - it may be appropriate to include the following card -
c     static dmach(3)
c
c     data minmag/20k,3*0/,maxmag/77777k,3*177777k/
c     data mcheps/32020k,3*0/
c
c     Machine constants for the Harris 220.
c
c     data mcheps(1),mcheps(2) / '20000000, '00000334 /
c     data minmag(1),minmag(2) / '20000000, '00000201 /
c     data maxmag(1),maxmag(2) / '37777777, '37777577 /
c
c     Machine constants for the Cray-1.
c
c     data mcheps(1) / 0376424000000000000000b /
c     data mcheps(2) / 0000000000000000000000b /
c
c     data minmag(1) / 0200034000000000000000b /
c     data minmag(2) / 0000000000000000000000b /
c
c     data maxmag(1) / 0577777777777777777777b /
c     data maxmag(2) / 0000007777777777777776b /
c
c     Machine constants for the Prime 400.
c
c     data mcheps(1),mcheps(2) / :10000000000, :00000000123 /
c     data minmag(1),minmag(2) / :10000000000, :00000100000 /
c     data maxmag(1),maxmag(2) / :17777777777, :37777677776 /
c
c     Machine constants for the VAX-11.
c
c     data mcheps(1),mcheps(2) /   9472,  0 /
c     data minmag(1),minmag(2) /    128,  0 /
c     data maxmag(1),maxmag(2) / -32769, -1 /
c
c     Machine constants for IEEE machines.
c
      data dmach(1) /2.22044604926d-16/
      data dmach(2) /2.22507385852d-308/
      data dmach(3) /1.79769313485d+308/
c
      mm10_dpmpar = dmach(i)
      return
c
c     Last card of function dpmpar.
c
      end
c
c
      subroutine mm10_dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
      integer n,lr
      double precision delta
      double precision r(lr),diag(n),qtb(n),x(n),wa1(n),wa2(n)
c     **********
c
c     subroutine dogleg
c
c     given an m by n matrix a, an n by n nonsingular diagonal
c     matrix d, an m-vector b, and a positive number delta, the
c     problem is to determine the convex combination x of the
c     gauss-newton and scaled gradient directions that minimizes
c     (a*x - b) in the least squares sense, subject to the
c     restriction that the euclidean norm of d*x be at most delta.
c
c     this subroutine completes the solution of the problem
c     if it is provided with the necessary information from the
c     qr factorization of a. that is, if a = q*r, where q has
c     orthogonal columns and r is an upper triangular matrix,
c     then dogleg expects the full upper triangle of r and
c     the first n components of (q transpose)*b.
c
c     the subroutine statement is
c
c       subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
c
c     where
c
c       n is a positive integer input variable set to the order of r.
c
c       r is an input array of length lr which must contain the upper
c         triangular matrix r stored by rows.
c
c       lr is a positive integer input variable not less than
c         (n*(n+1))/2.
c
c       diag is an input array of length n which must contain the
c         diagonal elements of the matrix d.
c
c       qtb is an input array of length n which must contain the first
c         n elements of the vector (q transpose)*b.
c
c       delta is a positive input variable which specifies an upper
c         bound on the euclidean norm of d*x.
c
c       x is an output array of length n which contains the desired
c         convex combination of the gauss-newton direction and the
c         scaled gradient direction.
c
c       wa1 and wa2 are work arrays of length n.
c
c     subprograms called
c
c       minpack-supplied ... dpmpar,enorm
c
c       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,j,jj,jp1,k,l
      double precision alpha,bnorm,epsmch,gnorm,one,qnorm,sgnorm,sum,
     *                 temp,zero
      double precision mm10_dpmpar,mm10_enorm
      data one,zero /1.0d0,0.0d0/
c
c     epsmch is the machine precision.
c
      epsmch = mm10_dpmpar(1)
c
c     first, calculate the gauss-newton direction.
c
      jj = (n*(n + 1))/2 + 1
      do 50 k = 1, n
         j = n - k + 1
         jp1 = j + 1
         jj = jj - k
         l = jj + 1
         sum = zero
         if (n .lt. jp1) go to 20
         do 10 i = jp1, n
            sum = sum + r(l)*x(i)
            l = l + 1
   10       continue
   20    continue
         temp = r(jj)
         if (temp .ne. zero) go to 40
         l = j
         do 30 i = 1, j
            temp = dmax1(temp,dabs(r(l)))
            l = l + n - i
   30       continue
         temp = epsmch*temp
         if (temp .eq. zero) temp = epsmch
   40    continue
         x(j) = (qtb(j) - sum)/temp
   50    continue
c
c     test whether the gauss-newton direction is acceptable.
c
      do 60 j = 1, n
         wa1(j) = zero
         wa2(j) = diag(j)*x(j)
   60    continue
      qnorm = mm10_enorm(n,wa2)
      if (qnorm .le. delta) go to 140
c
c     the gauss-newton direction is not acceptable.
c     next, calculate the scaled gradient direction.
c
      l = 1
      do 80 j = 1, n
         temp = qtb(j)
         do 70 i = j, n
            wa1(i) = wa1(i) + r(l)*temp
            l = l + 1
   70       continue
         wa1(j) = wa1(j)/diag(j)
   80    continue
c
c     calculate the norm of the scaled gradient and test for
c     the special case in which the scaled gradient is zero.
c
      gnorm = mm10_enorm(n,wa1)
      sgnorm = zero
      alpha = delta/qnorm
      if (gnorm .eq. zero) go to 120
c
c     calculate the point along the scaled gradient
c     at which the quadratic is minimized.
c
      do 90 j = 1, n
         wa1(j) = (wa1(j)/gnorm)/diag(j)
   90    continue
      l = 1
      do 110 j = 1, n
         sum = zero
         do 100 i = j, n
            sum = sum + r(l)*wa1(i)
            l = l + 1
  100       continue
         wa2(j) = sum
  110    continue
      temp = mm10_enorm(n,wa2)
      sgnorm = (gnorm/temp)/temp
c
c     test whether the scaled gradient direction is acceptable.
c
      alpha = zero
      if (sgnorm .ge. delta) go to 120
c
c     the scaled gradient direction is not acceptable.
c     finally, calculate the point along the dogleg
c     at which the quadratic is minimized.
c
      bnorm = mm10_enorm(n,qtb)
      temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta)
      temp = temp - (delta/qnorm)*(sgnorm/delta)**2
     *       + dsqrt((temp-(delta/qnorm))**2
     *               +(one-(delta/qnorm)**2)*(one-(sgnorm/delta)**2))
      alpha = ((delta/qnorm)*(one - (sgnorm/delta)**2))/temp
  120 continue
c
c     form appropriate convex combination of the gauss-newton
c     direction and the scaled gradient direction.
c
      temp = (one - alpha)*dmin1(sgnorm,delta)
      do 130 j = 1, n
         x(j) = temp*wa1(j) + alpha*x(j)
  130    continue
  140 continue
      return
c
c     last card of subroutine dogleg.
c
      end subroutine
c
c
      double precision function mm10_enorm(n,x)
      integer n
      double precision x(n)
c     **********
c
c     function enorm
c
c     given an n-vector x, this function calculates the
c     euclidean norm of x.
c
c     the euclidean norm is computed by accumulating the sum of
c     squares in three different sums. the sums of squares for the
c     small and large components are scaled so that no overflows
c     occur. non-destructive underflows are permitted. underflows
c     and overflows do not occur in the computation of the unscaled
c     sum of squares for the intermediate components.
c     the definitions of small, intermediate and large components
c     depend on two constants, rdwarf and rgiant. the main
c     restrictions on these constants are that rdwarf**2 not
c     underflow and rgiant**2 not overflow. the constants
c     given here are suitable for every known computer.
c
c     the function statement is
c
c       double precision function enorm(n,x)
c
c     where
c
c       n is a positive integer input variable.
c
c       x is an input array of length n.
c
c     subprograms called
c
c       fortran-supplied ... dabs,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i
      double precision agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs,
     *                 x1max,x3max,zero
      data one,zero,rdwarf,rgiant /1.0d0,0.0d0,3.834d-20,1.304d19/
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = n
      agiant = rgiant/floatn
      do 90 i = 1, n
         xabs = dabs(x(i))
         if (xabs .gt. rdwarf .and. xabs .lt. agiant) go to 70
            if (xabs .le. rdwarf) go to 30
c
c              sum for large components.
c
               if (xabs .le. x1max) go to 10
                  s1 = one + s1*(x1max/xabs)**2
                  x1max = xabs
                  go to 20
   10          continue
                  s1 = s1 + (xabs/x1max)**2
   20          continue
               go to 60
   30       continue
c
c              sum for small components.
c
               if (xabs .le. x3max) go to 40
                  s3 = one + s3*(x3max/xabs)**2
                  x3max = xabs
                  go to 50
   40          continue
                  if (xabs .ne. zero) s3 = s3 + (xabs/x3max)**2
   50          continue
   60       continue
            go to 80
   70    continue
c
c           sum for intermediate components.
c
            s2 = s2 + xabs**2
   80    continue
   90    continue
c
c     calculation of norm.
c
      if (s1 .eq. zero) go to 100
         mm10_enorm = x1max*dsqrt(s1+(s2/x1max)/x1max)
         go to 130
  100 continue
         if (s2 .eq. zero) go to 110
            if (s2 .ge. x3max)
     *         mm10_enorm = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))
            if (s2 .lt. x3max)
     *         mm10_enorm = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
            go to 120
  110    continue
            mm10_enorm = x3max*dsqrt(s3)
  120    continue
  130 continue
      return
c
c     last card of function enorm.
c
      end
c
c
      subroutine mm10_qform(m,n,q,ldq,wa)
      integer m,n,ldq
      double precision q(ldq,m),wa(m)
c     **********
c
c     subroutine qform
c
c     this subroutine proceeds from the computed qr factorization of
c     an m by n matrix a to accumulate the m by m orthogonal matrix
c     q from its factored form.
c
c     the subroutine statement is
c
c       subroutine qform(m,n,q,ldq,wa)
c
c     where
c
c       m is a positive integer input variable set to the number
c         of rows of a and the order of q.
c
c       n is a positive integer input variable set to the number
c         of columns of a.
c
c       q is an m by m array. on input the full lower trapezoid in
c         the first min(m,n) columns of q contains the factored form.
c         on output q has been accumulated into a square matrix.
c
c       ldq is a positive integer input variable not less than m
c         which specifies the leading dimension of the array q.
c
c       wa is a work array of length m.
c
c     subprograms called
c
c       fortran-supplied ... min0
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,j,jm1,k,l,minmn,np1
      double precision one,sum,temp,zero
      data one,zero /1.0d0,0.0d0/
c
c     zero out upper triangle of q in the first min(m,n) columns.
c
      minmn = min0(m,n)
      if (minmn .lt. 2) go to 30
      do 20 j = 2, minmn
         jm1 = j - 1
         do 10 i = 1, jm1
            q(i,j) = zero
   10       continue
   20    continue
   30 continue
c
c     initialize remaining columns to those of the identity matrix.
c
      np1 = n + 1
      if (m .lt. np1) go to 60
      do 50 j = np1, m
         do 40 i = 1, m
            q(i,j) = zero
   40       continue
         q(j,j) = one
   50    continue
   60 continue
c
c     accumulate q from its factored form.
c
      do 120 l = 1, minmn
         k = minmn - l + 1
         do 70 i = k, m
            wa(i) = q(i,k)
            q(i,k) = zero
   70       continue
         q(k,k) = one
         if (wa(k) .eq. zero) go to 110
         do 100 j = k, m
            sum = zero
            do 80 i = k, m
               sum = sum + q(i,j)*wa(i)
   80          continue
            temp = sum/wa(k)
            do 90 i = k, m
               q(i,j) = q(i,j) - temp*wa(i)
   90          continue
  100       continue
  110    continue
  120    continue
      return
c
c     last card of subroutine qform.
c
      end
c
c
      subroutine mm10_qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
      integer m,n,lda,lipvt
      integer ipvt(lipvt)
      logical pivot
      double precision a(lda,n),rdiag(n),acnorm(n),wa(n)
c     **********
c
c     subroutine qrfac
c
c     this subroutine uses householder transformations with column
c     pivoting (optional) to compute a qr factorization of the
c     m by n matrix a. that is, qrfac determines an orthogonal
c     matrix q, a permutation matrix p, and an upper trapezoidal
c     matrix r with diagonal elements of nonincreasing magnitude,
c     such that a*p = q*r. the householder transformation for
c     column k, k = 1,2,...,min(m,n), is of the form
c
c                           t
c           i - (1/u(k))*u*u
c
c     where u has zeros in the first k-1 positions. the form of
c     this transformation and the method of pivoting first
c     appeared in the corresponding linpack subroutine.
c
c     the subroutine statement is
c
c       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
c
c     where
c
c       m is a positive integer input variable set to the number
c         of rows of a.
c
c       n is a positive integer input variable set to the number
c         of columns of a.
c
c       a is an m by n array. on input a contains the matrix for
c         which the qr factorization is to be computed. on output
c         the strict upper trapezoidal part of a contains the strict
c         upper trapezoidal part of r, and the lower trapezoidal
c         part of a contains a factored form of q (the non-trivial
c         elements of the u vectors described above).
c
c       lda is a positive integer input variable not less than m
c         which specifies the leading dimension of the array a.
c
c       pivot is a logical input variable. if pivot is set true,
c         then column pivoting is enforced. if pivot is set false,
c         then no column pivoting is done.
c
c       ipvt is an integer output array of length lipvt. ipvt
c         defines the permutation matrix p such that a*p = q*r.
c         column j of p is column ipvt(j) of the identity matrix.
c         if pivot is false, ipvt is not referenced.
c
c       lipvt is a positive integer input variable. if pivot is false,
c         then lipvt may be as small as 1. if pivot is true, then
c         lipvt must be at least n.
c
c       rdiag is an output array of length n which contains the
c         diagonal elements of r.
c
c       acnorm is an output array of length n which contains the
c         norms of the corresponding columns of the input matrix a.
c         if this information is not needed, then acnorm can coincide
c         with rdiag.
c
c       wa is a work array of length n. if pivot is false, then wa
c         can coincide with rdiag.
c
c     subprograms called
c
c       minpack-supplied ... dpmpar,enorm
c
c       fortran-supplied ... dmax1,dsqrt,min0
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,j,jp1,k,kmax,minmn
      double precision ajnorm,epsmch,one,p05,sum,temp,zero
      double precision dpmpar,mm10_enorm
      data one,p05,zero /1.0d0,5.0d-2,0.0d0/
c
c     epsmch is the machine precision.
c
      epsmch = mm10_dpmpar(1)
c
c     compute the initial column norms and initialize several arrays.
c
      do 10 j = 1, n
         acnorm(j) = mm10_enorm(m,a(1,j))
         rdiag(j) = acnorm(j)
         wa(j) = rdiag(j)
         if (pivot) ipvt(j) = j
   10    continue
c
c     reduce a to r with householder transformations.
c
      minmn = min0(m,n)
      do 110 j = 1, minmn
         if (.not.pivot) go to 40
c
c        bring the column of largest norm into the pivot position.
c
         kmax = j
         do 20 k = j, n
            if (rdiag(k) .gt. rdiag(kmax)) kmax = k
   20       continue
         if (kmax .eq. j) go to 40
         do 30 i = 1, m
            temp = a(i,j)
            a(i,j) = a(i,kmax)
            a(i,kmax) = temp
   30       continue
         rdiag(kmax) = rdiag(j)
         wa(kmax) = wa(j)
         k = ipvt(j)
         ipvt(j) = ipvt(kmax)
         ipvt(kmax) = k
   40    continue
c
c        compute the householder transformation to reduce the
c        j-th column of a to a multiple of the j-th unit vector.
c
         ajnorm = mm10_enorm(m-j+1,a(j,j))
         if (ajnorm .eq. zero) go to 100
         if (a(j,j) .lt. zero) ajnorm = -ajnorm
         do 50 i = j, m
            a(i,j) = a(i,j)/ajnorm
   50       continue
         a(j,j) = a(j,j) + one
c
c        apply the transformation to the remaining columns
c        and update the norms.
c
         jp1 = j + 1
         if (n .lt. jp1) go to 100
         do 90 k = jp1, n
            sum = zero
            do 60 i = j, m
               sum = sum + a(i,j)*a(i,k)
   60          continue
            temp = sum/a(j,j)
            do 70 i = j, m
               a(i,k) = a(i,k) - temp*a(i,j)
   70          continue
            if (.not.pivot .or. rdiag(k) .eq. zero) go to 80
            temp = a(j,k)/rdiag(k)
            rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2))
            if (p05*(rdiag(k)/wa(k))**2 .gt. epsmch) go to 80
            rdiag(k) = mm10_enorm(m-j,a(jp1,k))
            wa(k) = rdiag(k)
   80       continue
   90       continue
  100    continue
         rdiag(j) = -ajnorm
  110    continue
      return
c
c     last card of subroutine qrfac.
c
      end
c
c
      subroutine mm10_r1mpyq(m,n,a,lda,v,w)
      integer m,n,lda
      double precision a(lda,n),v(n),w(n)
c     **********
c
c     subroutine r1mpyq
c
c     given an m by n matrix a, this subroutine computes a*q where
c     q is the product of 2*(n - 1) transformations
c
c           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
c
c     and gv(i), gw(i) are givens rotations in the (i,n) plane which
c     eliminate elements in the i-th and n-th planes, respectively.
c     q itself is not given, rather the information to recover the
c     gv, gw rotations is supplied.
c
c     the subroutine statement is
c
c       subroutine r1mpyq(m,n,a,lda,v,w)
c
c     where
c
c       m is a positive integer input variable set to the number
c         of rows of a.
c
c       n is a positive integer input variable set to the number
c         of columns of a.
c
c       a is an m by n array. on input a must contain the matrix
c         to be postmultiplied by the orthogonal matrix q
c         described above. on output a*q has replaced a.
c
c       lda is a positive integer input variable not less than m
c         which specifies the leading dimension of the array a.
c
c       v is an input array of length n. v(i) must contain the
c         information necessary to recover the givens rotation gv(i)
c         described above.
c
c       w is an input array of length n. w(i) must contain the
c         information necessary to recover the givens rotation gw(i)
c         described above.
c
c     subroutines called
c
c       fortran-supplied ... dabs,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,j,nmj,nm1
      double precision cos,one,sin,temp
      data one /1.0d0/
c
c     apply the first set of givens rotations to a.
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 50
      do 20 nmj = 1, nm1
         j = n - nmj
         if (dabs(v(j)) .gt. one) cos = one/v(j)
         if (dabs(v(j)) .gt. one) sin = dsqrt(one-cos**2)
         if (dabs(v(j)) .le. one) sin = v(j)
         if (dabs(v(j)) .le. one) cos = dsqrt(one-sin**2)
         do 10 i = 1, m
            temp = cos*a(i,j) - sin*a(i,n)
            a(i,n) = sin*a(i,j) + cos*a(i,n)
            a(i,j) = temp
   10       continue
   20    continue
c
c     apply the second set of givens rotations to a.
c
      do 40 j = 1, nm1
         if (dabs(w(j)) .gt. one) cos = one/w(j)
         if (dabs(w(j)) .gt. one) sin = dsqrt(one-cos**2)
         if (dabs(w(j)) .le. one) sin = w(j)
         if (dabs(w(j)) .le. one) cos = dsqrt(one-sin**2)
         do 30 i = 1, m
            temp = cos*a(i,j) + sin*a(i,n)
            a(i,n) = -sin*a(i,j) + cos*a(i,n)
            a(i,j) = temp
   30       continue
   40    continue
   50 continue
      return
c
c     last card of subroutine r1mpyq.
c
      end
c
c
      subroutine mm10_r1updt(m,n,s,ls,u,v,w,sing)
      integer m,n,ls
      logical sing
      double precision s(ls),u(m),v(n),w(m)
c     **********
c
c     subroutine r1updt
c
c     given an m by n lower trapezoidal matrix s, an m-vector u,
c     and an n-vector v, the problem is to determine an
c     orthogonal matrix q such that
c
c                   t
c           (s + u*v )*q
c
c     is again lower trapezoidal.
c
c     this subroutine determines q as the product of 2*(n - 1)
c     transformations
c
c           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
c
c     where gv(i), gw(i) are givens rotations in the (i,n) plane
c     which eliminate elements in the i-th and n-th planes,
c     respectively. q itself is not accumulated, rather the
c     information to recover the gv, gw rotations is returned.
c
c     the subroutine statement is
c
c       subroutine r1updt(m,n,s,ls,u,v,w,sing)
c
c     where
c
c       m is a positive integer input variable set to the number
c         of rows of s.
c
c       n is a positive integer input variable set to the number
c         of columns of s. n must not exceed m.
c
c       s is an array of length ls. on input s must contain the lower
c         trapezoidal matrix s stored by columns. on output s contains
c         the lower trapezoidal matrix produced as described above.
c
c       ls is a positive integer input variable not less than
c         (n*(2*m-n+1))/2.
c
c       u is an input array of length m which must contain the
c         vector u.
c
c       v is an array of length n. on input v must contain the vector
c         v. on output v(i) contains the information necessary to
c         recover the givens rotation gv(i) described above.
c
c       w is an output array of length m. w(i) contains information
c         necessary to recover the givens rotation gw(i) described
c         above.
c
c       sing is a logical output variable. sing is set true if any
c         of the diagonal elements of the output s are zero. otherwise
c         sing is set false.
c
c     subprograms called
c
c       minpack-supplied ... dpmpar
c
c       fortran-supplied ... dabs,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more,
c     john l. nazareth
c
c     **********
      integer i,j,jj,l,nmj,nm1
      double precision cos,cotan,giant,one,p5,p25,sin,tan,tau,temp,
     *                 zero
      double precision dpmpar
      data one,p5,p25,zero /1.0d0,5.0d-1,2.5d-1,0.0d0/
c
c     giant is the largest magnitude.
c
      giant = mm10_dpmpar(3)
c
c     initialize the diagonal element pointer.
c
      jj = (n*(2*m - n + 1))/2 - (m - n)
c
c     move the nontrivial part of the last column of s into w.
c
      l = jj
      do 10 i = n, m
         w(i) = s(l)
         l = l + 1
   10    continue
c
c     rotate the vector v into a multiple of the n-th unit vector
c     in such a way that a spike is introduced into w.
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 nmj = 1, nm1
         j = n - nmj
         jj = jj - (m - j + 1)
         w(j) = zero
         if (v(j) .eq. zero) go to 50
c
c        determine a givens rotation which eliminates the
c        j-th element of v.
c
         if (dabs(v(n)) .ge. dabs(v(j))) go to 20
            cotan = v(n)/v(j)
            sin = p5/dsqrt(p25+p25*cotan**2)
            cos = sin*cotan
            tau = one
            if (dabs(cos)*giant .gt. one) tau = one/cos
            go to 30
   20    continue
            tan = v(j)/v(n)
            cos = p5/dsqrt(p25+p25*tan**2)
            sin = cos*tan
            tau = sin
   30    continue
c
c        apply the transformation to v and store the information
c        necessary to recover the givens rotation.
c
         v(n) = sin*v(j) + cos*v(n)
         v(j) = tau
c
c        apply the transformation to s and extend the spike in w.
c
         l = jj
         do 40 i = j, m
            temp = cos*s(l) - sin*w(i)
            w(i) = sin*s(l) + cos*w(i)
            s(l) = temp
            l = l + 1
   40       continue
   50    continue
   60    continue
   70 continue
c
c     add the spike from the rank 1 update to w.
c
      do 80 i = 1, m
         w(i) = w(i) + v(n)*u(i)
   80    continue
c
c     eliminate the spike.
c
      sing = .false.
      if (nm1 .lt. 1) go to 140
      do 130 j = 1, nm1
         if (w(j) .eq. zero) go to 120
c
c        determine a givens rotation which eliminates the
c        j-th element of the spike.
c
         if (dabs(s(jj)) .ge. dabs(w(j))) go to 90
            cotan = s(jj)/w(j)
            sin = p5/dsqrt(p25+p25*cotan**2)
            cos = sin*cotan
            tau = one
            if (dabs(cos)*giant .gt. one) tau = one/cos
            go to 100
   90    continue
            tan = w(j)/s(jj)
            cos = p5/dsqrt(p25+p25*tan**2)
            sin = cos*tan
            tau = sin
  100    continue
c
c        apply the transformation to s and reduce the spike in w.
c
         l = jj
         do 110 i = j, m
            temp = cos*s(l) + sin*w(i)
            w(i) = -sin*s(l) + cos*w(i)
            s(l) = temp
            l = l + 1
  110       continue
c
c        store the information necessary to recover the
c        givens rotation.
c
         w(j) = tau
  120    continue
c
c        test for zero diagonal elements in the output s.
c
         if (s(jj) .eq. zero) sing = .true.
         jj = jj + (m - j + 1)
  130    continue
  140 continue
c
c     move w back into the last column of the output s.
c
      l = jj
      do 150 i = n, m
         s(l) = w(i)
         l = l + 1
  150    continue
      if (s(jj) .eq. zero) sing = .true.
      return
c
c     last card of subroutine r1updt.
c
      end
c
c
c    End of trust region dogleg solver
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_nwclsh                       *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 7/03/15                     *
c     *                                                              *
c     *     Cubic line search from NLEQSLV package                   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine mm10_nwclsh(n,xc,fcnorm,d,g,stepmx,xtol,scalex,
     *                  xp,fp,fpnorm,xw,retcd,gcnt,priter,iter,
     &                  props, np1, np0, vec1, vec2, gp)
c
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, np0
      double precision, dimension(max_uhard) :: vec1,vec2
      integer n,retcd,gcnt
      double precision  stepmx,xtol,fcnorm,fpnorm
      double precision  xc(*)
      double precision  d(*),g(*),xp(*),fp(*),xw(*)
      double precision  scalex(*)
c
      integer priter,iter
c
c-------------------------------------------------------------------------
c
c     Find a next acceptable iterate using a safeguarded cubic line search
c     along the newton direction
c
c     Arguments
c
c     In       n       Integer          dimension of problem
c     In       xc      Real(*)          current iterate
c     In       fcnorm  Real             0.5 * || f(xc) ||**2
c     In       d       Real(*)          newton direction
c     In       g       Real(*)          gradient at current iterate
c     In       stepmx  Real             maximum stepsize
c     In       xtol    Real             relative step size at which
c                                       successive iterates are considered
c                                       close enough to terminate algorithm
c     In       scalex  Real(*)          diagonal scaling matrix for x()
c     In       fvec    Name             name of routine to calculate f()
c     In       xp      Real(*)          new x()
c     In       fp      Real(*)          new f(x)
c     In       fpnorm  Real             .5*||fp||**2
c     Out      xw      Real(*)           workspace for unscaling x(*)
c
c     Out      retcd   Integer          return code
c                                         0 new satisfactory x() found
c                                         1 no  satisfactory x() found
c                                           sufficiently distinct from xc()
c
c     Out      gcnt    Integer          number of steps taken
c     In       priter  Integer           >0 unit if intermediate steps to be printed
c                                        -1 if no printing
c
c-------------------------------------------------------------------------

      integer i, gp
      double precision  alpha,slope,rsclen,oarg(4)
      double precision  lambda,lamhi,lamlo,t
      double precision  ddot,dnrm2, nudnrm, ftarg
      double precision  dlen
      double precision a, b, disc, fpt, fpt0, fpnorm0, lambda0
      integer idamax, miter
      logical firstback

      parameter (alpha = 1.0d-4)

      double precision Rhalf, Rone, Rtwo, Rthree, Rten, Rzero
      parameter(Rzero=0.0d0)
      parameter(Rhalf=0.5d0, Rone=1.0d0, Rtwo=2.0d0, Rten=10.0d0)
      parameter(Rthree=3.0d0)
      double precision t1,t2
      double precision dlamch

c     silence warnings issued by ftncheck

      lambda0 = Rzero
      fpnorm0 = Rzero

c     safeguard initial step size

      dlen = dnrm2(n,d,1)
c      if( dlen .gt. stepmx ) then
c          lamhi  = stepmx / dlen
c      else
          lamhi  = Rone
c      endif

c     compute slope  =  g-trans * d

      slope = ddot(n,g,1,d,1)

c     compute the smallest value allowable for the damping
c     parameter lambda ==> lamlo

      rsclen = nudnrm(n,d,xc)
      lamlo  = xtol / rsclen

c     initialization of retcd and lambda (linesearch length)

      retcd  = 2
      lambda = lamhi
      gcnt   = 0
      firstback = .true.
      iter = 0
      miter = 200

      do while( retcd .eq. 2 .and. iter .lt. miter)
         iter = iter + 1

c        compute next x
         do i=1,n
            xp(i) = xc(i) + lambda*d(i)
         enddo

c        evaluate functions and the objective function at xp

c         call nwfvec(xp,n,scalex,fvec,fp,fpnorm,xw)
c     NOTE: in nwfvec, they unscale x first; thus, needs special
c           treatment if scalex ~= 1.d0
      call mm10_formR(props, np1, np0, vec1, vec2, xp(1:6),
     & xp(7:6+props%num_hard), fp, gp)
      fpnorm = 0.5d0*dot_product(fp(1:6+props%num_hard)
     &     ,fp(1:6+props%num_hard))
         gcnt = gcnt + 1
         ftarg = fcnorm + alpha * lambda * slope

         if( priter .gt. 0) then
            oarg(1) = lambda
            oarg(2) = ftarg
            oarg(3) = fpnorm
            oarg(4) = abs(fp(idamax(n,fp,1)))
c        write(*,*) "iter=", iter, "gcnt=", gcnt, 
c     &      "ftarg=",ftarg,"fpnorm=",fpnorm,
c     &      "fcnorm=",fcnorm,"slope=",slope
         endif

c        first is quadratic
c        test whether the standard step produces enough decrease
c        of the objective function.
c        If not update lambda and compute a new next iterate

         if( fpnorm .le. ftarg ) then
             retcd = 0
         else
            if( fpnorm .gt. lamlo**2 * sqrt(dlamch('O')) ) then
c               safety against overflow in what follows (use lamlo**2 for safety)
                lambda = lambda/Rten
                firstback = .true.
            else
                if( firstback ) then
                   t = ((-lambda**2)*slope/Rtwo)/
     *                        (fpnorm-fcnorm-lambda*slope)
                   firstback = .false.
                else
                   fpt  = fpnorm -fcnorm - lambda*slope
                   fpt0 = fpnorm0-fcnorm - lambda0*slope
                   a = fpt/lambda**2 - fpt0/lambda0**2
                   b = -lambda0*fpt/lambda**2 + lambda*fpt0/lambda0**2
                   a = a /(lambda - lambda0)
                   b = b /(lambda - lambda0)
                   if( abs(a) .le. dlamch('E') ) then
c                      not quadratic but linear
                       t = -slope/(2*b)
                   else
c                      use Higham procedure to compute roots acccurately
c                      Higham: Accuracy and Stability of Numerical Algorithms, second edition,2002, page 10.    
c                      Actually solving 3*a*x^2+2*b*x+c=0 ==> (3/2)*a*x^2+b*x+(c/2)=0
                       disc = b**2 - Rthree * a * slope
                       t1 = -(b+sign(Rone,b)*sqrt(disc))/(Rthree*a)
                       t2 = slope/(Rthree*a)/t1
                       if(a .gt. Rzero ) then
c                          upward opening parabola ==> rightmost is solution
                           t = max(t1,t2)
                       else
c                          downward opening parabola ==> leftmost is solution
                           t = min(t1,t2)
                       endif
                   endif
                   t = min(t, Rhalf*lambda)
                endif
                lambda0 = lambda
                fpnorm0 = fpnorm
                lambda = max(t,lambda/Rten)
                if(lambda .lt. lamlo) then
                   retcd = 1
                endif
            endif
         endif
      enddo

      if(iter.eq.miter) then
        retcd = 1
      endif

      return
      end
c
c-----------------------------------------------------------------------

      function nudnrm(n, d, x)
      integer n
      double precision  d(*), x(*)
      double precision nudnrm

c-------------------------------------------------------------------------
c
c     calculate  max( abs(d[*]) / max(x[*],1) )
c
c     Arguments
c
c     In   n        Integer       number of elements in d() and x()
c     In   d        Real(*)       vector d
c     In   x        Real(*)       vector x
c
c-------------------------------------------------------------------------

      integer i
      double precision  t

      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)

      t = Rzero
      do i=1,n
         t = max(t, abs(d(i)) / max(abs(x(i)),Rone))
      enddo
      nudnrm = t

      return
      end
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
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_num_tangent                  *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/27/13                    *
c     *                                                              *
c     *     Calculate the consistent tangent after a converged       *
c     *     stress update with a numerical derivative (for check)    *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_num_tangent(props, np1, n)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      type(crystal_state) :: np1star
      integer :: i
      double precision, dimension(6) :: dD
      double precision :: eps
      logical :: cut
      parameter (eps=1d-8)
c
      np1%tangent = 0.0d0
      do i=1,6
            dD = 0.0d0
            dD(i) = eps
            call mm10_setup_np1(reshape(np1%R, (/9/)), np1%D+dD, 
     &            np1%tinc, 
     &            np1%temp, np1%step, np1%elem, np1%gp, np1star)
            call mm10_solve_strup(props, np1star, n, cut)
            if (cut) then
              write (*,*) "Numerical tangent failed"
              call die_gracefully
            end if
            np1%tangent(:,i) = (np1star%stress - np1%stress) / eps
      end do

      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_ur_tangent                   *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 1/16/14                     *
c     *                                                              *
c     *     Calculate the consistent tangent after a converged       *
c     *     stress update.                                           *
c     *                                                              *
c     *     This routine uses the old, obsolete unrolling method     *
c     *     for computing the tangent.  It's just for debug/         *
c     *     comparison purposes.                                     *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_ur_tangent(props, np1, n)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      double precision, dimension(6,6) :: A
      double precision, dimension(6) :: d_mod, d_barp, tw
      double precision, dimension(6) :: b
      double precision, dimension(3) :: w_p
      double precision :: alpha
      double precision, dimension(42,42) :: sys
      double precision, dimension(42) :: rhs
      double precision, dimension(7,6) :: AA
      double precision, dimension(7,7) :: Jac
      integer :: i, j, k, r_ind, c_ind, nes, nrhs, lda, ldb, info
      integer, dimension(42) :: ipiv
c
      call mm10_formJ(props, np1, n, np1%stress, np1%tau_tilde, Jac)
      call mm10_form_dbarp(props, np1, n, np1%stress, np1%tau_tilde,
     &      d_barp)
      call mm10_form_wp(props, np1, n, np1%stress, np1%tau_tilde,
     &      w_p)
      call mm10_symSW(np1%stress, w_p, tw)
c
      if (props%h_type .eq. 1) then ! voche
        call mm10_ed_voche(props, np1, n, np1%stress, np1%tau_tilde, b)
      elseif (props%h_type .eq. 2) then ! MTS
        call mm10_ed_mts(props, np1, n, np1%stress, np1%tau_tilde, b)
      elseif (props%h_type .eq. 3) then ! User
        call mm10_ed_user(props, np1, n, np1%stress, np1%tau_tilde, b)
      else
        call mm10_unknown_hard_error(props)
      end if
c
      d_mod = np1%D
      d_mod(4:6) = 0.5d0 * d_mod(4:6)
c
      b = -b
      A = 0.0d0
      A = -props%stiffness
      alpha = 2.0d0/(3.0d0*np1%dg**2.0d0)
      call DGER(6, 6, alpha, matmul(props%stiffness,d_barp) + 2.0d0*tw,
     &      1, d_mod, 1, A, 6)
c
      np1%tangent = 0.0d0
c
c     Unroll
c
      AA(1:6,1:6) = -A
      AA(7,1:6) = -b
      rhs = 0.0d0
      sys = 0.0d0
c
      do i=1,7
        do j=1,6
          r_ind = (i-1)*6+j
          rhs(r_ind) = AA(i,j)
          do k=1,7
            c_ind = (k-1)*6+j
            sys(r_ind, c_ind) = Jac(i,k)
          end do
        end do
      end do

c     Solve the equation
      nes = 42
      nrhs = 1
      lda = 42
      ldb = 42
      call DGESV(nes, nrhs, sys, lda, ipiv, rhs, ldb, info)
c
      np1%tangent = 0.0

      do i=1,6
        np1%tangent(i,1:6) = rhs(((i-1)*6+1):i*6)
      end do


      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ                        *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *     Form the jacobian from lower subroutines                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ(props, np1, n, vec1, vec2, arr1, arr2,
     &           stress, tt, J)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision,
     & dimension(6+props%num_hard,6+props%num_hard) :: J
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      call mm10_formarrs(props, np1, n, stress, tt, vec1, vec2,
     &       arr1, arr2,2)
c
      call mm10_formJ11(props, np1, n, vec1, vec2, arr1, arr2,
     &   stress, tt, J(1:6,1:6))
c      write (*,*) "J11", J(1,1)
      call mm10_formJ12(props, np1, n, vec1, vec2, arr1, arr2,
     & stress, tt,
     & J(1:6,7:6+props%num_hard))
c      write (*,*) "J12", J(1,7)
      call mm10_formJ21(props, np1, n, vec1, vec2, arr1, arr2,
     & stress, tt,
     & J(7:6+props%num_hard, 1:6))
c      write (*,*) "J21", J(7,1)
      call mm10_formJ22(props, np1, n, vec1, vec2, arr1, arr2,
     & stress, tt,
     & J(7:6+props%num_hard,7:6+props%num_hard))
c      write (*,*) "J22", J(7,7)
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJi                       *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *     Form the jacobian from lower subroutines                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJi(props, np1, n, ivec1, ivec2,
     &           stress, tt, J)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision,
     & dimension(6+props%num_hard,6+props%num_hard) :: J
      double precision, dimension(props%num_hard) :: tt
      double complex, dimension(max_uhard) :: ivec1, ivec2
c
      call mm10_formJ11i(props, np1, n, ivec1, ivec2, 
     &   stress, tt, J(1:6,1:6))
c      write (*,*) "J11", J(1,1)
      call mm10_formJ12i(props, np1, n, ivec1, ivec2, 
     & stress, tt,
     & J(1:6,7:6+props%num_hard))
c      write (*,*) "J12", J(1,7)
      call mm10_formJ21i(props, np1, n, ivec1, ivec2, 
     & stress, tt,
     & J(7:6+props%num_hard, 1:6))
c      write (*,*) "J21", J(7,1)
      call mm10_formJ22i(props, np1, n, ivec1, ivec2, 
     & stress, tt,
     & J(7:6+props%num_hard,7:6+props%num_hard))
c      write (*,*) "J22", J(7,7)
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_form_numJ                    *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *     Form the jacobian from numerical differentiation         *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_form_numJ(props, np1, n, stress, tt, J)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(7,7) :: J
      double precision :: tt
c
      double precision :: eps, pt
      double precision, dimension(7) :: R, pR
      double precision, dimension(6) :: pS
      integer :: i
      parameter(eps = 1d-6)
c
      J = 0.0d0
      call mm10_formR(props, np1, n, stress, tt, R)
      do i=1,7
        if (i .ne. 7) then
          pS = stress
          pS(i) = pS(i) + eps
          call mm10_formR(props, np1, n, pS, tt, pR)
        else
          pt = tt + eps
          call mm10_formR(props, np1, n, stress, pt, pR)
        end if
        J(:,i) = (pR-R)/eps
      end do
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formR                        *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *     Form the residual from lower subroutines                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formR(props, np1, n, vec1, vec2,
     &   stress, tt, R, gp)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(6+props%num_hard) :: R
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
      integer :: gp
c
      call mm10_formvecs(props, np1, n, stress, tt,
     &    vec1, vec2)
      call mm10_formR1(props, np1, n, vec1, vec2,stress, tt, R(1:6), gp)
      call mm10_formR2(props, np1, n, vec1, vec2,stress, tt,
     & R(7:6+props%num_hard),gp)
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formR1                       *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *     Form R1                                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formR1(props, np1, n, vec1, vec2,
     &     stress, tt, R1,gp)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(6) :: R1
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
c
      double precision, dimension(6) :: dbarp
      double precision, dimension(3) :: wp
      double precision, dimension(6) :: symTW,temp
      integer :: gp
c
      call mm10_form_dbarp(props, np1, n, vec1, vec2, 
     &   stress, tt, dbarp)
      call mm10_form_wp(props, np1, n, vec1, vec2,stress, tt, wp)
      call mm10_symSW(stress, wp, symTW)
c
      R1 = stress - n%stress - matmul(props%stiffness, np1%D - dbarp) 
     &      + 2.0d0 * symTW
c      if(gp.eq.0) then
c       write(*,*) "R1=", R1(2), " de=", np1%D(2) - dbarp(2)
c       write(*,*) "dp=", dbarp(2), " d=", np1%D(2)
c       write(*,*) "wp=", wp(2), "STW=", symTW(2)
c       temp = matmul(props%stiffness, np1%D - dbarp)
c       write(*,*) "C*de= ", temp(2)
c       temp = stress - n%stress
c       write(*,*) "sn-sn_1= ", temp(2)
c       write(*,*) "C= ", props%stiffness(1,1)
c       write(*,*) "slipinc= ", vec1(6)
c      endif
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formR1i                      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 11/26/13                    *
c     *                                                              *
c     *     Form R1                                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formR1i(props, np1, n, ivec1, ivec2,
     &     stress, tt, R1)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double complex, dimension(6) :: stress
      double complex, dimension(6) :: R1
      double complex, dimension(props%num_hard) :: tt
      double complex, dimension(max_uhard) :: ivec1, ivec2
c
      double complex, dimension(6) :: dbarp, temp
      double complex, dimension(6,6) :: stiff2
      double precision, dimension(6,6) :: zeroff
      double complex, dimension(3) :: wp
      double complex, dimension(6) :: symTW
c
      call mm10_form_dbarpi(props, np1, n, ivec1, ivec2, 
     &   stress, tt, dbarp)
      call mm10_form_wpi(props, np1, n, ivec1, ivec2,stress, tt, wp)
      call mm10_symSWi(stress, wp, symTW)
c
      zeroff = 0.d0
      stiff2 = dcmplx(props%stiffness, zeroff)
      temp = np1%D - dbarp
      R1 = stress - n%stress - matmul(stiff2, np1%D - dbarp) 
     &      + 2.0d0 * symTW
c
      return
      end subroutine
c
c *****************************************************************************
c *                                                                           *
c *         Internal utility routines -- initializing history, setting gauss  *
c *               point data, etc.                                            *
c *                                                                           *
c *****************************************************************************
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_dump_state                   *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 2/25/14                     *
c     *                                                              *
c     *     Dump a state variable to try to figure out what's going  *
c     *     on with these random errors.                             *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_dump_state(state)
      use mm10_defs
      implicit none
c
      type(crystal_state) :: state
c
      write(*,*)
      write(*,*) "step, element, gp"
      write(*,*) state%step, state%elem, state%gp
      write(*,*) "R"
      write(*,*) transpose(state%R)
      write(*,*) "Rp"
      write(*,*) transpose(state%Rp)
      write(*,*) "stress"
      write(*,*) state%stress
      write(*,*) "D"
      write(*,*) state%D
      write(*,*) "angles"
      write(*,*) state%euler_angles
      write(*,*) "tau_l"
      write(*,*) state%tau_l
      write(*,*) "slip_incs"
      write(*,*) state%slip_incs
      write(*,*) "gradFeinv"
      write(*,*) state%gradFeInv
      write(*,*) "tangent"
      write(*,*) state%tangent
      write(*,*) "tau_tilde, temp, tinc, dg, tau_v, tau_y"
      write(*,*) state%tau_tilde, state%temp, state%tinc, state%dg,
     &            state%tau_v, state%tau_y
      write(*,*) "mu_harden, work_inc, p_work_inc, p_strain_inc"
      write(*,*) state%mu_harden, state%work_inc, state%p_work_inc,
     &            state%p_strain_inc
      write(*,*) "ms"
      write(*,*) transpose(state%ms)
      write(*,*) "qs"
      write(*,*) transpose(state%qs)
      write(*,*) "qc"
      write(*,*) transpose(state%qc)
      write(*,*) "u"
      write(*,*) state%u
      write(*,*)
c
      end subroutine
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_dump_props                   *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 2/25/14                     *
c     *                                                              *
c     *     Dump a props variable to try to figure out what's going  *
c     *     on with these random errors.                             *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_dump_props(props)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
c
      write(*,*) "PROPS"
      write(*,*) props%rate_n, props%tau_hat_y, props%G_0_y, 
     &                  props%burgers,
     &                  props%p_v, props%q_v, props%boltzman, 
     &                  props%theta_0, props%eps_dot_0_v,
     &                  props%eps_dot_0_y,
     &                  props%p_y, props%q_y,
     &                  props%tau_a, props%tau_hat_v, props%G_0_v,
     &                  props%k_0, props%mu_0, props%D_0, props%T_0, 
     &                  props%tau_y, props%tau_v, props%voche_m,
     &                  props%u1, props%u2, props%u3, props%u4, 
     &                  props%u5, props%u6
      write(*,*) transpose(props%g)
      write(*,*) props%ms
      write(*,*) props%qs
      write(*,*) props%ns
      write(*,*) transpose(props%stiffness)
      write(*,*) props%angle_type, props%angle_convention, props%nslip,
     &      props%h_type
c
      end subroutine

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
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_setup_np1                    *
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
      info_vector(4) = 39+max_slip_sys+max_uhard
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
      num_states = 39 + max_slip_sys + max_uhard
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

      do i = 1, max_slip_sys
            write(state_labels(i+39), 9000) i
            state_descriptors(i+39) = "integrated slip"
      end do

      do i = 1, max_uhard
            write(state_labels(i+max_slip_sys+39), 9020)
     & i
            state_descriptors(i+max_slip_sys+39) =
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
c           39+1:39+max_slip_sys    slip history      76:76+max_slip_sys-1
c           39+max_slip_sys+1:end    hardening      76:76+max_slip_sys-1
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
c           Results 39+1:39+max_slip_sys : the slip totals, 
c           padded with zeros as 
c           required  and averaged over Gauss points
c
      s = 76
      e = 76 + nslip - 1
      one_elem_states(39+1:39+max_slip_sys) = 
     &         sum( history_dump(s:e,1:int_points,relem),2 ) / 
     &         dble(int_points)

c
c           Results 13: tau_tilde of first crystal, avged
c
      s = sc + 30 + max_slip_sys + 1
      e = sc + 30 + max_slip_sys + max_uhard
      one_elem_states(40+max_slip_sys:39+max_slip_sys+max_uhard) = 
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
c *****************************************************************************
c *                                                                           *
c *         START HELPER ROUTINES                                             *
c *                                                                           *
c *****************************************************************************
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
c *    mm10_symSW                                                            *    *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 12/11/12 mcm                                     *
c *                                                                          *
c *   Take the symmetric part of S*W for some stress and skew tensor         *
c *                                                                          *
c ****************************************************************************
c
      subroutine mm10_symSW(S, W, SW)
            implicit none
            double precision, dimension(6), intent(in) :: S
            double precision, dimension(3), intent(in) :: W
            double precision, dimension(6), intent(out) :: SW
c
            SW = 0.0d0
            SW(1) = S(4)*W(3) - S(6)*W(2)
            SW(2) = S(4)*W(3) - S(5)*W(1)
            SW(3) = S(6)*W(2) + S(5)*W(1)
            SW(4) = 0.5d0*(W(3)*(S(1)-S(2)) + W(1)*S(6) - W(2)*S(5))
            SW(5) = 0.5d0*(W(1)*(S(2)-S(3)) + W(2)*S(4) + W(3)*S(6))
            SW(6) = 0.5d0*(W(2)*(S(1)-S(3)) + W(1)*S(4) - W(3)*S(5))

            return

      end subroutine
c
c ****************************************************************************
c *                                                                          *
c *    mm10_symSWi                                                           *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 12/11/12 mcm                                     *
c *                                                                          *
c *   Take the symmetric part of S*W for some stress and skew tensor         *
c *                                                                          *
c ****************************************************************************
c
      subroutine mm10_symSWi(S, W, SW)
            implicit none
            double complex, dimension(6), intent(in) :: S
            double complex, dimension(3), intent(in) :: W
            double complex, dimension(6), intent(out) :: SW
c
            SW = 0.d0
            SW(1) = S(4)*W(3) - S(6)*W(2)
            SW(2) = S(4)*W(3) - S(5)*W(1)
            SW(3) = S(6)*W(2) + S(5)*W(1)
            SW(4) = 0.5d0*(W(3)*(S(1)-S(2)) + W(1)*S(6) - W(2)*S(5))
            SW(5) = 0.5d0*(W(1)*(S(2)-S(3)) + W(2)*S(4) + W(3)*S(6))
            SW(6) = 0.5d0*(W(2)*(S(1)-S(3)) + W(1)*S(4) - W(3)*S(5))

            return

      end subroutine
c
c ****************************************************************************
c *                                                                          *
c *    mm10_symSWmat                                                         *
c *                                                                          *
c *         written by : tjt                                                 *
c *         last modified : 03/31/15 tjt                                     *
c *                                                                          *
c *   Take the symmetric part of S*W for some stress and skew tensor         *
c *                                                                          *
c ****************************************************************************
c
      subroutine mm10_symSWmat(S, W, n, SW)
            implicit none
            integer n
            double precision, dimension(6), intent(in) :: S
            double precision, dimension(3,n), intent(in) :: W
            double precision, dimension(6,n), intent(out) :: SW
c
            SW = 0.0d0
            SW(1,1:n) = S(4)*W(3,1:n) - S(6)*W(2,1:n)
            SW(2,1:n) = S(4)*W(3,1:n) - S(5)*W(1,1:n)
            SW(3,1:n) = S(6)*W(2,1:n) + S(5)*W(1,1:n)
            SW(4,1:n) = 0.5d0*(W(3,1:n)*(S(1)-S(2))
     & + W(1,1:n)*S(6) - W(2,1:n)*S(5))
            SW(5,1:n) = 0.5d0*(W(1,1:n)*(S(2)-S(3))
     & + W(2,1:n)*S(4) + W(3,1:n)*S(6))
            SW(6,1:n) = 0.5d0*(W(2,1:n)*(S(1)-S(3))
     & + W(1,1:n)*S(4) - W(3,1:n)*S(5))

            return

      end subroutine
c
c
c ****************************************************************************
c *                                                                          *
c *    mm10_IW                                                               *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 12/11/12 mcm                                     *
c *                                                                          *
c *   Get Iik*Wlj - Wik*Ijl in our voigt notation                            *
c *                                                                          *
c ****************************************************************************
c
      subroutine mm10_IW(W, IW)
            implicit none
            double precision, dimension(3), intent(in) :: W
            double precision, dimension(6,6), intent(out) :: IW
c
            IW = 0.0d0
            IW(1,4) =  2.0d0*W(3)
            IW(1,6) = -2.0d0*W(2)
            IW(2,4) =  2.0d0*W(3)
            IW(2,5) = -2.0d0*W(1)
            IW(3,5) =  2.0d0*W(1)
            IW(3,6) =  2.0d0*W(2)
            IW(4,1) =  W(3)
            IW(4,2) = -W(3)
            IW(4,5) = -W(2)
            IW(4,6) =  W(1)
            IW(5,2) =  W(1)
            IW(5,3) = -W(1)
            IW(5,4) =  W(2)
            IW(5,6) =  W(3)
            IW(6,1) =  W(2)
            IW(6,3) = -W(2)
            IW(6,4) =  W(1)
            IW(6,5) = -W(3)

            return
      end subroutine
c
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
c
c -----------------------------------------------------------------------------
c
c    Slip rate and hardening functions
c
c -----------------------------------------------------------------------------
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
c
c           Actual user sliprate function
      subroutine mm10_slipinc_user(props, np1, n, stress, tt, i,
     &                             slipinc)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      double precision :: slipinc
      integer :: i
c
      write (*,*) "Not implemented"
      call die_gracefully
c
      return
      end subroutine
c
c           Actual user hardening function
      subroutine mm10_h_user(props, np1, n, stress, tt, h)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt, h
c
      write (*,*) "Not implemented"
      call die_gracefully
c
      return
      end subroutine
c
c           Derivative of hardening fn wrt stress
      subroutine mm10_estress_user(props, np1, n, stress, tt, et)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(props%num_hard,6) :: et
c
      write (*,*) "Not implemented"
      call die_gracefully
c
      return
      end subroutine
c
c           Derivative of hardening fn wrt hardening
      subroutine mm10_ehard_user(props, np1, n, stress, tt, etau)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      double precision, 
     &   dimension(props%num_hard,props%num_hard) :: etau
c
      write (*,*) "Not implemented"
      call die_gracefully
c
      return
      end subroutine
c
c
c           Derivative of hardening fn wrt strain
      subroutine mm10_ed_user(props, np1, n, stress, tt, ed)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(6,props%num_hard) :: ed
c
      write (*,*) "Not implemented"
      call die_gracefully
c
      return
      end subroutine
c
c           Derivative of sliprate wrt resolved shear stress
      subroutine mm10_dgdt_user(props, np1, n, stress, tt, dgammadtau)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: rs
      double precision, dimension(props%nslip) :: dgammadtau
      double precision, dimension(props%num_hard) :: tt
c
c
      write (*,*) "Not implemented"
      call die_gracefully
c
      return
      end subroutine
c
c           Derivative of sliprate wrt hardening variables
      subroutine mm10_dgdh_user(props, np1, n, stress, tt, dgammadtt)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(props%nslip,props%num_hard)
     &        :: dgammadtt
c
c
      write (*,*) "Not implemented"
      call die_gracefully
c
      return
      end subroutine
c
c           Derivative of sliprate wrt strain increment
      subroutine mm10_dgdd_user(props, np1, n, stress, tt, dgammadd)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(props%nslip,6) :: dgammadd
c
c
      write (*,*) "Not implemented"
      call die_gracefully
c
      return
      end subroutine
c
c *****************************************************************************
c *                                                                           *
c *         Built in hardening routines                                       *
c *                                                                           *
c *****************************************************************************
c
c
c           Calculate the slip increment along system i
      subroutine mm10_slipinc(props, np1, n, stress, tt, i,
     &                             slipinc)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt
      integer :: i
c
      double precision :: slipinc, mm10_rs
      double precision :: rs
c
      rs = mm10_rs(props, np1, n, stress, tt, i)
      slipinc = np1%dg/tt * dabs(rs/tt)**(props%rate_n-1.0)*rs
      return
      end subroutine
c      
c           Calculate the resolved shear along system i     
      function mm10_rs(props, np1, n, stress, tt, i)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt
      integer :: i
c
      double precision :: mm10_rs
c
      mm10_rs = dot_product(stress, np1%ms(1:6,i))
c
      return
      end function
c      
c           Calculate the resolved shear along system i     
      function mm10_rsi(props, np1, n, stress, tt, i)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double complex, dimension(6) :: stress
      double complex :: tt
      integer :: i
c
      double complex :: mm10_rsi
c
      mm10_rsi = stress(1)*np1%ms(1,i)
     &         + stress(2)*np1%ms(2,i)
     &         + stress(3)*np1%ms(3,i)
     &         + stress(4)*np1%ms(4,i)
     &         + stress(5)*np1%ms(5,i)
     &         + stress(6)*np1%ms(6,i)
c
      return
      end function
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
c
c           Actual voche law hardening function
      subroutine mm10_h_voche(props, np1, n, stress, tt, h)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(1) :: tt, h
      double precision :: dg
      integer :: i
c
      double precision :: slipinc
c
      h = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
        h = h + (1.0d0-(tt-props%tau_y)/props%tau_v+np1%tau_l(i)/(
     &      tt-props%tau_y))**(props%voche_m)*
     &      dabs(slipinc)
      end do
      h(1) = n%tau_tilde(1) + props%theta_0*h(1)
c
      return
      end subroutine
c
c           Derivative of hardening fn wrt stress
      subroutine mm10_estress_voche(props, np1, n, stress, tt, et)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(1) :: tt
      double precision, dimension(6) :: et
c
      double precision :: mm10_rs
      double precision :: rs
      integer :: i
c
      et = 0.0d0
c
      do i=1,props%nslip
        rs = mm10_rs(props, np1, n, stress, tt, i)
        et = et + (1.0d0-(tt(1)-props%tau_y)/props%tau_v+np1%tau_l(i)
     &      /(tt(1)-props%tau_y))**(props%voche_m)*
     &      dabs(rs)**(props%rate_n-2.d0)*rs*np1%ms(1:6,i)
      end do
c
      et = props%theta_0*np1%dg*props%rate_n/tt(1)**
     &      (props%rate_n)*et
c
      return
      end subroutine
c
c           Derivative of hardening fn wrt hardening
      subroutine mm10_ehard_voche(props, np1, n, stress, tt, etau)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt, etau
c
      double precision :: slipinc
      integer :: i
c
      etau = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
        etau = etau + (props%voche_m*(1.0d0/props%tau_v+np1%tau_l(i)/
     &      (tt-props%tau_y)**2.0d0)*(1.0d0-(tt-props%tau_y)/
     &       props%tau_v+
     &      np1%tau_l(i)/(tt-props%tau_y))**(-1.0d0)+props%rate_n/tt)*
     &      (1.0d0-(tt-props%tau_y)/props%tau_v+np1%tau_l(i)/
     &      (tt-props%tau_y))**(props%voche_m)*
     &      dabs(slipinc)
      end do

      etau = -props%theta_0*etau
      etau = 1.0d0 - etau
c
      return
      end subroutine
c
c           Derivative of hardening fn wrt strain
      subroutine mm10_ed_voche(props, np1, n, stress, tt, ed)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(1) :: tt
      double precision, dimension(6) :: ed
c
      double precision, dimension(1) :: h
      double precision, dimension(6) :: d_mod
c
      ed = 0.0d0
c
      call mm10_h_voche(props, np1, n, stress, tt, h)
      d_mod = np1%D
      d_mod(4:6) = 0.5d0 * d_mod(4:6)
c
      ed = 2.0d0*(h(1) - n%tau_tilde(1))/
     &     (3.0d0*np1%dg**2.0d0)*d_mod
c
      return
      end subroutine
c
c           Derivative of sliprate wrt resolved shear stress
      subroutine mm10_dgdt_voche(props, np1, n, stress, tt, dgammadtau)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt, rs
      double precision, dimension(props%nslip) :: dgammadtau
c
      double precision :: mm10_rs
c
      integer :: s
c
      do s=1,props%nslip
        rs = mm10_rs(props, np1, n, stress, tt, s)
        dgammadtau(s) = dabs(rs)**(props%rate_n-1.0d0)
        dgammadtau(s) = np1%dg*props%rate_n/tt**(props%rate_n)
     &     *dgammadtau(s)
      end do
c
      return
      end subroutine
c
c           Derivative of sliprate wrt hardening variables
      subroutine mm10_dgdh_voche(props, np1, n, stress, tt, dgammadtt)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt, dgam
      double precision, dimension(props%nslip,1) :: dgammadtt
c
      double precision :: slipinc
c
      integer :: s
c
      do s=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, s, dgam)
        dgammadtt(s,1) = -props%rate_n/tt * dgam
      end do
c
      return
      end subroutine
c
c           Derivative of sliprate wrt strain increment
      subroutine mm10_dgdd_voche(props, np1, n, stress, tt, D, 
     &                         dgammadd)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, D, d_mod
      double precision :: tt, alpha, dgam
      double precision, dimension(props%nslip,6) :: dgammadd
c
      double precision :: slipinc
c
      integer :: s
c
      d_mod = D
      d_mod(4:6) = 0.5d0 * d_mod(4:6)
      alpha = 2.0d0/(3.0d0*np1%dg**2.0d0)
      do s=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, s, dgam)
        dgammadd(s,1:6) = alpha * dgam * d_mod(1:6)
      end do
c
      return
      end subroutine
c
c
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
c           Actual MTS hardening function
      subroutine mm10_h_mts(props, np1, n, stress, tt, h)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(1) :: tt, h
      integer :: i
c
      double precision :: slipinc
      double precision :: ct, cta
c
      cta = (props%mu_0/
     &      np1%mu_harden)*tt(1) - (props%mu_0/np1%mu_harden)*
     &     props%tau_a - np1%tau_y
      ct = 1.0d0 - cta/np1%tau_v
      h = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
        h(1) = h(1) + (ct + np1%tau_l(i)/cta)**(props%voche_m)*
     &      dabs(slipinc)
      end do
c
      h(1) = props%tau_a*(1.0d0 - np1%mu_harden/n%mu_harden) + 
     &      (np1%mu_harden/props%mu_0)*(np1%tau_y - n%tau_y) + 
     &      (np1%mu_harden/n%mu_harden)*n%tau_tilde(1) + 
     &      props%theta_0 * (np1%mu_harden/props%mu_0)*h(1)

      return
      end subroutine
c
c           Derivative of hardening fn wrt stress
      subroutine mm10_estress_mts(props, np1, n, stress, tt, et)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt
      double precision, dimension(6) :: et
c
      double precision :: mm10_rs
c
      double precision :: rs, cta, ct
      integer :: i
c
      cta = (props%mu_0/
     &      np1%mu_harden)*tt - (props%mu_0/np1%mu_harden)*props%tau_a -
     &      np1%tau_y
      ct = 1.0d0 - cta/np1%tau_v
      et = 0.0d0
      do i = 1, props%nslip
        rs = mm10_rs(props, np1, n, stress, tt, i)
        et = et + (ct+np1%tau_l(i)/cta)**(props%voche_m)*
     &      dabs(rs)**(props%rate_n-2.0d0)*rs*np1%ms(1:6,i)
      end do

      et =  props%theta_0 * 
     &      (np1%mu_harden/props%mu_0)*et*props%rate_n*
     &      np1%dg/tt**props%rate_n

c
      return
      end subroutine
c
c           Derivative of hardening fn wrt hardening
      subroutine mm10_ehard_mts(props, np1, n, stress, tt, etau)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt, etau
c
      double precision :: slipinc
c
      double precision :: ur, s, ct, cta
      integer :: i
c
      cta = (props%mu_0/
     &      np1%mu_harden)*tt - (props%mu_0/np1%mu_harden)*props%tau_a -
     &      np1%tau_y
      ct = 1.0d0 - cta/np1%tau_v
c
      ur = np1%mu_harden / props%mu_0
c
      etau = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
        etau = etau+(props%voche_m*(1.d0/np1%tau_v+np1%tau_l(i)/
     &         cta**2.0d0)*
     &      (ct+np1%tau_l(i)/cta)**(-1.0d0) + ur*props%rate_n/tt)*
     &      (ct+np1%tau_l(i)/cta)**(props%voche_m)*
     &      dabs(slipinc)
      end do
c
      etau = -props%theta_0*etau
      etau = 1.0d0 - etau

      return
      end subroutine
c
c           Derivative of hardening fn wrt strain
      subroutine mm10_ed_mts(props, np1, n, stress, tt, ed)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt
      double precision, dimension(6) :: ed
c
      double precision :: slipinc
c
      double precision :: dslip, lnv, lny, dgc, ty, tv, mnp0, m0np, sc
      double precision, dimension(6) :: d_mod, dydd, dvdd
      integer :: s
c
c     Form a bunch of simple constitutive things
c
c
      d_mod = np1%D
      d_mod(4:6) = 0.5d0 * d_mod(4:6)
c
      dgc = np1%dg / np1%tinc
c
c     Form d_tauy/d_deltad by steps
c
      lny = dlog(props%eps_dot_0_y/dgc)
      ty = props%boltzman*np1%temp
     &       /(np1%mu_harden*(props%burgers**3)*props%G_0_y)*
     &       lny
      dydd = 2.0d0*props%tau_hat_y/(3.0d0*np1%dg**2.0d0*
     &        props%q_y*props%p_y*
     &      lny)*(1.0d0-ty**(1.0d0/props%q_y))**(1.0d0/props%p_y-1.0d0)*
     &      ty**(1.0d0/props%q_y)*d_mod

c
c     Form d_tauv/d_deltad by steps
c
      lnv = dlog(props%eps_dot_0_v/dgc)
      tv = props%boltzman*np1%temp
     &       /(np1%mu_harden*(props%burgers**3)*props%G_0_v)*
     &       lnv
      dvdd = 2.0d0*props%tau_hat_v/(3.0d0*np1%dg**2.0d0*
     &        props%q_v*props%p_v*
     &      lnv)*(1.0d0-tv**(1.0d0/props%q_v))**(1.0d0/props%p_v-1.0d0)*
     &      tv**(1.0d0/props%q_v)*d_mod

c
c     Form a couple more common components
c
      mnp0 = np1%mu_harden / props%mu_0
      sc = tt/mnp0 - props%tau_a/mnp0 - np1%tau_y
c
c     Glue everything together
c
      ed = 0.0d0
      do s=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, s, slipinc)
        ed = ed + (props%voche_m*(1.d0/np1%tau_v
     &      +np1%tau_l(s)/sc**2.0d0)*
     &      (1.0d0-sc/np1%tau_v+np1%tau_l(s)/sc)
     &         **(props%voche_m-1.0d0)*dydd
     &      + props%voche_m/(np1%tau_v)**2.0d0*sc*
     &      (1.0d0-sc/np1%tau_v+np1%tau_l(s)/sc)
     &         **(props%voche_m-1.0)*dvdd
     &      + 2.0d0/(3.0d0*np1%dg**2.0d0)*(1.0d0-sc/np1%tau_v
     &         +np1%tau_l(s)/sc)
     &      **(props%voche_m)*d_mod) * 
     &      dabs(slipinc)
      end do

      ed = props%theta_0 * mnp0 * ed + mnp0*dydd

      return
      end subroutine
c
c           Derivative of sliprate wrt resolved shear stress
      subroutine mm10_dgdt_mts(props, np1, n, stress, tt, dgammadtau)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt, rs
      double precision, dimension(props%nslip) :: dgammadtau
c
      double precision :: mm10_rs
c
      integer :: s
c
      do s=1,props%nslip
        rs = mm10_rs(props, np1, n, stress, tt, s)
        dgammadtau(s) = dabs(rs)**(props%rate_n-1.0d0)
        dgammadtau(s) = np1%dg*props%rate_n/tt**(props%rate_n)
     &     *dgammadtau(s)
      end do
c
      return
      end subroutine
c
c           Derivative of sliprate wrt hardening variables
      subroutine mm10_dgdh_mts(props, np1, n, stress, tt, dgammadtt)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt, dgam
      double precision, dimension(props%nslip,1) :: dgammadtt
c
      double precision :: slipinc
c
      integer :: s
c
      do s=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, s, dgam)
        dgammadtt(s,1) = -props%rate_n/tt * dgam
      end do
c
      return
      end subroutine
c
c           Derivative of sliprate wrt strain increment
      subroutine mm10_dgdd_mts(props, np1, n, stress, tt, D, 
     &                         dgammadd)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, D, d_mod
      double precision :: tt, alpha, dgam
      double precision, dimension(6,props%nslip) :: dgammadd
c
      double precision :: slipinc
c
      integer :: s
c
      d_mod = D
      d_mod(4:6) = 0.5d0 * d_mod(4:6)
      alpha = 2.0d0/(3.0d0*np1%dg**2.0d0)
      do s=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, s, dgam)
        dgammadd(1:6,s) = alpha * dgam * d_mod(1:6)
      end do
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
c           Form intermediate vectors for faster calculations
      subroutine mm10_v_mrr(props, np1, n, stress, tt, 
     &   vec1, vec2)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(12) :: tt, h
      double precision, dimension(max_uhard) :: vec1, vec2
      integer :: alpha
      double precision :: slipinc
c
      do alpha = 1,props%nslip
         call mm10_slipinc_mrr(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
         vec1(alpha) = slipinc
      enddo
c
      end subroutine
c
c           Form intermediate arrays for faster calculations
      subroutine mm10_a_mrr(props, np1, n, stress, tt, 
     &   vec1, vec2, arr1, arr2, both)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(12) :: tt, h
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
      integer both
c
         call mm10_dgdt_mrr(props, np1, n, stress, tt, 
     &                            arr1(1:props%num_hard,1))
      if(both.eq.2) then
         call mm10_dgdh_mrr(props, np1, n, stress, tt, 
     &               arr2(1:props%nslip,1:props%num_hard))
      endif
c
      end subroutine
c
c           Form intermediate vectors for faster calculations
      subroutine mm10_vi_mrr(props, np1, n, stress, tt, 
     &   ivec1, ivec2)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double complex, dimension(6) :: stress
      double complex, dimension(12) :: tt, h
      double complex, dimension(max_uhard) :: ivec1, ivec2
      integer :: alpha
      double complex :: slipinc
c
      do alpha = 1,props%nslip
         call mm10_slipinci_mrr(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
         ivec1(alpha) = slipinc
      enddo
c
      end subroutine
c
c           Actual mrr sliprate function
      subroutine mm10_slipinc_mrr(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(12) :: tt
      double precision :: h, slipinc, mm10_rs
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m
      double precision, dimension(12,12) :: Gmat, Hmat
c
c Load some material parameters
        dt = np1%tinc
        k = props%boltzman
        theta = np1%temp
        G = props%mu_0
        b = props%burgers
        c1 = props%u1
        c2 = props%u2
        c3 = props%u3
        p_e = props%p_v
        q_e = props%q_v
        Qslip = props%G_0_v
        v_attack = props%G_0_y
c        
c Compute the shear modulus using Roter's function
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      call mm10_mrr_GH(props,Gmat,Hmat)
c        
c      ms = np1.ms(1:6,alpha);
c      rs = stress*ms; % tau^a
      rs = mm10_rs(props, np1, n, stress, tt, alpha)
c        
c         [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          rhoF = dot_product(Gmat(alpha,1:12),tt(1:12))
          rhoP = dot_product(Hmat(alpha,1:12),tt(1:12))
c
c Compute some stresses and rates
        gamma_0 = v_attack*k*theta/(c1*c3*G*b*b)*dsqrt(rhoP) ! (15)
        tpass = c1*G*b*dsqrt(rhoP) ! (16)
        tcut = Qslip/(c2*c3*b*b)*dsqrt(rhoF) ! (17)
        fract = ((dabs(rs)-tpass)/tcut)
c
c Evaluate the slip rate equation
c        if(alpha.eq.6) then
c          write(*,*) "fract=", fract, "rs=", rs
c          write(*,*) "tcut=", tcut, "tpass=", tpass
c          write(*,*) "gamma_0=", gamma_0, "dt=", dt
c          write(*,*) "exp=", dexp (-(Qslip/k/theta)*(1.d0 - fract))
c        endif

        if(fract.gt.1.d0) then
            ! linear extrapolation past the too-high stress (rs) value
            b = gamma_0
            x = fract
            m = b * (-q_e*(Qslip/k/theta)*(1.d0 - 1.d0)
     &  **(q_e-1.d0))
     &         * (- p_e*(1.d0)**(p_e-1.d0)) * dsign(1.d0,rs)/tcut
            y = m*x + b
            slipinc = dt * y
        elseif(dabs(rs).gt.0.d0) then
        slipinc = dt * gamma_0 * dsign(1.d0,rs)
     & * dexp (-(Qslip/k/theta)*(1.d0 - fract**p_e)**q_e) !(14)
        else
            slipinc = 0.d0
        endif
c       if(np1.step.gt.21) then
c       write(*,*) alpha, "rs", rs, "fract", fract, "slipinc", slipinc
c       endif
c
      return
      end subroutine
c
c           Actual mrr hardening function
      subroutine mm10_h_mrr(props, np1, n, vec1, vec2, 
     & stress, tt, h,gp)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(12) :: tt, h
      double precision, dimension(max_uhard) :: vec1, vec2
      integer :: alpha, gp
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, rho,
     &  rho_n, pi, c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  tem1, tem2, tem3, tem4
      double precision, dimension(12,12) :: Gmat, Hmat
c Load some material parameters
        Qbulk = props%tau_a
        k = props%boltzman
        theta = np1%temp
        if(theta.le.0.d0) then
          write (*,*) "Roters model requires non-zero 
     &   nodal temperatures"
      call die_gracefully
        endif
        c1 = props%u1
        c2 = props%u2
        c3 = props%u3
        c4 = props%u4
        c5 = props%u5
        c6 = props%u6
        c7 = props%tau_hat_y
        c8 = props%tau_hat_v
        G = props%mu_0
        b = props%burgers
        v = 0.3d0 !props%nu
        dt = np1%tinc
        PI=4.D0*DATAN(1.D0)
c      write(*,*) "pi", pi
c        
c Compute the shear modulus using Roter's function
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      call mm10_mrr_GH(props,Gmat,Hmat)
c
c
c      write(*,*) "Gmat", Gmat(1,1)
c      write(*,*) "G", G
      do alpha = 1,12
c
c Get dislocation density
          rho = tt(alpha) ! rho^a_SSD
          rho_n = n%tau_tilde(alpha) ! rho^a_SSD
c
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
c        write(*,*) "rs", rs
c          
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          rhoF = dot_product(Gmat(alpha,1:12),tt(1:12))
          rhoP = dot_product(Hmat(alpha,1:12),tt(1:12))
c        write(*,*) "rhoF", rhoF
c        write(*,*) "rhoP", rhoP
c          
          tpass = c1*G*b*dsqrt(rhoP) ! (16)
          ddipole = dsqrt(3.d0)*G*b/(16.d0*pi*(1.d0-v))/
     &       (dabs(rs)) ! (42)
          rhoM = (2.d0*k/(c1*c2*c3*G*b**3.d0))*
     &       theta*dsqrt(rhoF*rhoP) ! (13)
c        write(*,*) "tpass", tpass
c        write(*,*) "ddipole", ddipole
c        write(*,*) "rhoM", rhoM
c          
c          call mm10_slipinc_mrr(props, np1, n, stress, tt, alpha, 
c     &     slipinc)
          slipinc = vec1(alpha)
c       if(np1.step.gt.21) then
c       write(*,*) alpha, "slipinc", slipinc
c       endif
          gammadot = dabs(slipinc/dt)
c      if((alpha.eq.6).and.(gp.eq.1)) then
c      write(*,*) "rhoF=", rhoF, "rhoP=", rhoP
c      write(*,*) "tpass=", tpass, "dd=", ddipole
c      write(*,*) "rhoM=", rhoM, "gdot=", gammadot
c      endif
c        write(*,*) "gammadot", gammadot
c              
c Evaluate the hardening equation
c      write(*,*) "c4", c4
c      write(*,*) "c5", c5
c      write(*,*) "c6", c6
c      write(*,*) "c7", c7
c      write(*,*) "c8", c8
c      write(*,*) "Qbulk", Qbulk
c      write(*,*) "rho", rho
c      write(*,*) "b", b
      tem1 = c4/b*dsqrt(rhoP)*gammadot
      tem2 = c6*ddipole/b*rhoM*gammadot
      tem3 = c5*rho*gammadot
      tem4 = c7*dexp(-Qbulk/k/theta)*dabs(rs)/(k*theta)
     &        *rho*rho*gammadot**c8
c      if((alpha.eq.6).and.(gp.eq.1)) then
c      write(*,*) "h1=", tem1, "h2=", tem2
c      write(*,*) "h3=", tem3, "h4=", tem4
c      endif
          h(alpha) = rho_n + dt*(tem1
     &        + tem2 - tem3
     &        - tem4) ! (18)
      enddo
c
      return
      end subroutine
c
c           Imaginary mrr sliprate function
      subroutine mm10_slipinci_mrr(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double complex, dimension(6) :: stress
      double complex, dimension(12) :: tt, temp
      double complex :: h, slipinc, mm10_rsi
      integer :: alpha, i
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44
      double complex :: rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m
      double precision, dimension(12,12) :: Gmat, Hmat
c
c Load some material parameters
        dt = np1%tinc
        k = props%boltzman
        theta = np1%temp
        G = props%mu_0
        b = props%burgers
        c1 = props%u1
        c2 = props%u2
        c3 = props%u3
        p_e = props%p_v
        q_e = props%q_v
        Qslip = props%G_0_v
        v_attack = props%G_0_y
c        
c Compute the shear modulus using Roter's function
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      call mm10_mrr_GH(props,Gmat,Hmat)
c        
c      ms = np1.ms(1:6,alpha);
c      rs = stress*ms; % tau^a
      rs = mm10_rsi(props, np1, n, stress, tt, alpha)
c        
c         [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          temp = (Gmat(alpha,1:12)*tt(1:12))
          rhoF = sum(temp)
          temp = (Hmat(alpha,1:12)*tt(1:12))
          rhoP = sum(temp)
c
c Compute some stresses and rates
        gamma_0 = v_attack*k*theta/(c1*c3*G*b*b)*cdsqrt(rhoP) ! (15)
        tpass = c1*G*b*cdsqrt(rhoP) ! (16)
        tcut = Qslip/(c2*c3*b*b)*cdsqrt(rhoF) ! (17)
          if(dreal(rs).lt.0.d0) then
        fract = (-rs-tpass)/tcut
          else
        fract = (rs-tpass)/tcut
          endif
c
c Evaluate the slip rate equation
        if(dreal(fract).gt.1.d0) then
            ! linear extrapolation past the too-high stress (rs) value
            b = gamma_0
            x = fract
            m = b * (-q_e*(Qslip/k/theta)*(1.d0 - 1.d0)
     &  **(q_e-1.d0))
     &         * (- p_e*(1.d0)**(p_e-1.d0)) * dsign(1.d0,dreal(rs))/tcut
            y = m*x + b
            slipinc = dt * y
        elseif(dabs(dreal(rs)).gt.0.d0) then
        slipinc = dt * gamma_0 * dsign(1.d0,dreal(rs))
     & * cdexp (-(Qslip/k/theta)*(1.d0 - fract**p_e)**q_e) !(14)
        else
            slipinc = 0.d0
        endif
c
      return
      end subroutine
c
c           Imaginary mrr hardening function
      subroutine mm10_hi_mrr(props, np1, n, ivec1, ivec2, 
     & stress, tt, h)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double complex, dimension(6) :: stress
      double complex, dimension(12) :: tt, h, temp
      double complex, dimension(max_uhard) :: ivec1, ivec2
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, rho_n, 
     &  pi, c4, c5, c6, c7, c8, v, Qbulk
      double complex :: rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, rho,
     &  mm10_rsi,
     &  ddipole, rhoM, slipinc, gammadot, 
     &  tem1, tem2, tem3
      double precision, dimension(12,12) :: Gmat, Hmat
c Load some material parameters
        Qbulk = props%tau_a
        k = props%boltzman
        theta = np1%temp
        if(theta.le.0.d0) then
          write (*,*) "Roters model requires non-zero 
     &   nodal temperatures"
      call die_gracefully
        endif
        c1 = props%u1
        c2 = props%u2
        c3 = props%u3
        c4 = props%u4
        c5 = props%u5
        c6 = props%u6
        c7 = props%tau_hat_y
        c8 = props%tau_hat_v
        G = props%mu_0
        b = props%burgers
        v = 0.3d0 !props%nu
        dt = np1%tinc
        PI=4.D0*DATAN(1.D0)
c      write(*,*) "pi", pi
c        
c Compute the shear modulus using Roter's function
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      call mm10_mrr_GH(props,Gmat,Hmat)
c
c
c      write(*,*) "Gmat", Gmat(1,1)
c      write(*,*) "G", G
      do alpha = 1,12
c
c Get dislocation density
          rho = tt(alpha) ! rho^a_SSD
          rho_n = n%tau_tilde(alpha) ! rho^a_SSD
c
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rsi(props, np1, n, stress, tt, alpha)
          rs = dsign(1.d0,dreal(rs))*rs
c        write(*,*) "rs", rs
c          
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          temp = (Gmat(alpha,1:12)*tt(1:12))
          rhoF = sum(temp)
          temp = (Hmat(alpha,1:12)*tt(1:12))
          rhoP = sum(temp)
c        write(*,*) "rhoF", rhoF
c        write(*,*) "rhoP", rhoP
c          
          tpass = c1*G*b*cdsqrt(rhoP) ! (16)
          ddipole = dsqrt(3.d0)*G*b/(16.d0*pi*(1.d0-v))/
     &       (rs) ! (42)
          rhoM = (2.d0*k/(c1*c2*c3*G*b**3.d0))*
     &       theta*cdsqrt(rhoF*rhoP) ! (13)
c        write(*,*) "tpass", tpass
c        write(*,*) "ddipole", ddipole
c        write(*,*) "rhoM", rhoM
c          
c          call mm10_slipinc_mrr(props, np1, n, stress, tt, alpha, 
c     &     slipinc)
          slipinc = ivec1(alpha)
          slipinc = dsign(1.d0,dreal(slipinc))*slipinc
          gammadot = slipinc/dt
c        write(*,*) "gammadot", gammadot
c              
c Evaluate the hardening equation
c      write(*,*) "c4", c4
c      write(*,*) "c5", c5
c      write(*,*) "c6", c6
c      write(*,*) "c7", c7
c      write(*,*) "c8", c8
c      write(*,*) "Qbulk", Qbulk
c      write(*,*) "rho", rho
c      write(*,*) "b", b
      tem1 = c4/b*cdsqrt(rhoF)*gammadot
      tem2 = c6*ddipole/b*rhoM*gammadot
      tem3 = c5*rho*gammadot
c      write(*,*) "c1", tem1
c      write(*,*) "c2", tem2
c      write(*,*) "c3", tem3
      tem1 = c7*dexp(-Qbulk/k/theta)*rs/(k*theta)
     &        *rho*rho*gammadot**c8
c      write(*,*) "c1", tem1
          h(alpha) = rho_n + dt*(c4/b*cdsqrt(rhoP)*gammadot
     &        + c6*ddipole/b*rhoM*gammadot - c5*rho*gammadot
     &        - c7*dexp(-Qbulk/k/theta)*rs/(k*theta)
     &        *rho*rho*gammadot**c8) ! (18)
      enddo
c
      return
      end subroutine
c
c           Wrapper version, mrr slipinc function
      subroutine mm10_slipinc_mrrW(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(12) :: tt, zerosV
      double precision :: h, slipinc, mm10_rs
      integer :: alpha
c
      double complex, dimension(6) :: stressi
      double complex, dimension(12) :: tti
      double complex :: slipinci
c
      zerosV = 0.d0
      stressi = dcmplx(stress,zerosV(1:6))
      tti = dcmplx(tt,zerosV)
c
      call mm10_slipinci_mrr(props, np1, n, stressi, tti, 
     &                            alpha, slipinci)
c
      slipinc = dreal(slipinci)
c
      return
      end subroutine
c
c           Derivative of hardening fn wrt stress
      subroutine mm10_estress_mrr(props, np1, n, vec1, vec2, 
     &        arr1, arr2, stress, tt, et)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, dtdstress
      double precision, dimension(12) :: tt
      double precision, dimension(12,6) :: et
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, rho,
     &  rho_n, pi, c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  dddipole, dslipinc, badterm
      double precision, dimension(12,12) :: Gmat, Hmat
      double precision, dimension(props%nslip) :: dslip
      integer :: alpha
c Load some material parameters
        Qbulk = props%tau_a
        k = props%boltzman
        theta = np1%temp
        c1 = props%u1
        c2 = props%u2
        c3 = props%u3
        c4 = props%u4
        c5 = props%u5
        c6 = props%u6
        c7 = props%tau_hat_y
        c8 = props%tau_hat_v
        G = props%mu_0
        b = props%burgers
        v = 0.3d0 !props%nu
        dt = np1%tinc
        PI=4.D0*DATAN(1.D0)
c
c compute derivatives of slip increments with respect to resolved
c shear stress
c        call mm10_dgdt_mrr(props, np1, n, stress, 
c     &         tt, dslip)
        dslip(1:props%num_hard) = arr1(1:props%num_hard,1)
c        
c Compute the shear modulus using Roter's function
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      call mm10_mrr_GH(props,Gmat,Hmat)
c
      do alpha = 1,12

          ! Get dislocation density
          rho = tt(alpha) ! rho^a_SSD
c
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
c
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          rhoF = dot_product(Gmat(alpha,1:12),tt(1:12))
          rhoP = dot_product(Hmat(alpha,1:12),tt(1:12))
c
          tpass = c1*G*b*dsqrt(rhoP) ! (16)
          ddipole = dsqrt(3.d0)*G*b/(16.d0*pi*(1.d0-v))/
     &       (dabs(rs)) ! (42)
          dddipole = -dsqrt(3.d0)*G*b/(16.d0*pi*(1.d0-v))/
     &          (dabs(rs))**2.d0*dsign(1.d0,rs)
          rhoM = (2.d0*k/(c1*c2*c3*G*b**3.d0))*
     &       theta*dsqrt(rhoF*rhoP) ! (13)
c          
c          call mm10_slipinc_mrr(props, np1, n, stress, tt, alpha, 
c     &     slipinc)
          slipinc = vec1(alpha)
          gammadot = dabs(slipinc/dt)
          dtdstress(1:6) = np1%ms(1:6,alpha)
          dslipinc = dslip(alpha)/dt ! always positive

          ! Evaluate the equation
          if(gammadot.eq.0.d0) then
          badterm = 0.d0
          else
          badterm = c8*dabs(rs)*gammadot**(c8-1.d0)*dslipinc
          endif
          et(alpha,1:6) = dt*(c4/b*dsqrt(rhoF)*dsign(1.d0,rs)
     &         *dslipinc + c6/b*rhoM*(ddipole*dsign(1.d0,rs)
     &         *dslipinc + dddipole*gammadot)
     &        - c5*rho*dsign(1.d0,rs)*dslipinc
     &    - c7*dexp(-Qbulk/k/theta)/(k*theta)*rho**2.d0*
     &      (badterm + dsign(1.d0,rs)*gammadot**c8))
     &    *dtdstress(1:6) ! d(18)
      
      enddo
c
      return
      end subroutine
c
c           Derivative of hardening fn wrt hardening
      subroutine mm10_ehard_mrr(props, np1, n, vec1, vec2,
     &       arr1, arr2, stress, tt, etau)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(12) :: tt
      double precision, dimension(12,12) :: etau
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, rho,
     &  rho_n, pi, c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  dddipole, dslipinc, badterm, deltaij,
     &  drhoF, drhoP, drhoM
      double precision, dimension(12,12) :: Gmat, Hmat
      double precision, dimension(12,12) :: dslip
      integer :: alpha, beta
c      
c compute derivatives of slip increments with respect to densities
c        call mm10_dgdh_mrr(props, np1, n, stress, tt, dslip)
        dslip(1:props%nslip,1:props%num_hard) = 
     &   arr2(1:props%nslip,1:props%num_hard)
c
c Load some material parameters
        Qbulk = props%tau_a
        k = props%boltzman
        theta = np1%temp
        c1 = props%u1
        c2 = props%u2
        c3 = props%u3
        c4 = props%u4
        c5 = props%u5
        c6 = props%u6
        c7 = props%tau_hat_y
        c8 = props%tau_hat_v
        G = props%mu_0
        b = props%burgers
        v = 0.3d0 !props%nu
        dt = np1%tinc
        PI=4.D0*DATAN(1.D0)
c       
c Compute the shear modulus using Roter's function
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      call mm10_mrr_GH(props,Gmat,Hmat)
c
c Compute drho_alpha/drho_beta
c loop over numerator hardening variable
      do alpha = 1,12

c Get dislocation density
        rho = tt(alpha) ! rho^a_SSD
c
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
          
c         [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          rhoF = dot_product(Gmat(alpha,1:12),tt(1:12))
          rhoP = dot_product(Hmat(alpha,1:12),tt(1:12))
c
c          call mm10_slipinc_mrr(props, np1, n, stress, tt, alpha, 
c     &     slipinc)
          slipinc = vec1(alpha)
          gammadot = dabs(slipinc/dt)
          
        tpass = c1*G*b*dsqrt(rhoP) ! (16)
        ddipole = dsqrt(3.d0)*G*b/(16.d0*pi*(1.d0-v))/
     &       (dabs(rs)) ! (42)
        rhoM = (2.d0*k/(c1*c2*c3*G*b**3.d0))*
     &       theta*dsqrt(rhoF*rhoP) ! (13)
c        write(*,*) "tpass", tpass
c        write(*,*) "ddipole", ddipole
c        write(*,*) "rhoM", rhoM
c          
c loop over denominator hardening variable
        do beta = 1,12
c          
c           [drhoF,drhoP] = mm10_drhoFP_mrr(props, np1, n, tt, alpha, beta);
          drhoF = Gmat(alpha,beta)
          drhoP = Hmat(alpha,beta)
          
          dddipole = 0.d0
          drhoM = 0.5d0*(2.d0*k/(c1*c2*c3*G*b**3.d0))*
     &       theta/dsqrt(rhoF*rhoP) * (drhoF*rhoP + rhoF*drhoP)
          
          dslipinc = dslip(alpha,beta)/dt
c        write(*,*) "dslipinc", dslipinc
          
          if(alpha.eq.beta) then
              deltaij = 1.d0
          else
              deltaij = 0.d0
          endif

c Evaluate the equation
          if(gammadot.eq.0.d0) then
          badterm = 0.d0
          else
          badterm = c8*gammadot**(c8-1.d0)*dslipinc
          endif
          etau(alpha,beta) = deltaij - dt*(c4/b*(0.5d0*drhoF/
     &        dsqrt(rhoF)*gammadot*dsign(1.d0,rs) 
     &     + dsqrt(rhoF)*dslipinc) + c6/b*(rhoM*ddipole*dslipinc
     &     + rhoM*dddipole*gammadot*dsign(1.d0,rs)
     &     + drhoM*ddipole*gammadot*dsign(1.d0,rs))
     &     - c5*(rho*dslipinc + deltaij*gammadot*dsign(1.d0,rs))
     &     - c7*dexp(-Qbulk/k/theta)/(k*theta)*dabs(rs)*
     &      (rho**2.d0*badterm + 2.d0*rho*dsign(1.d0,rs)
     &      *deltaij*gammadot**c8))*dsign(1.d0,rs)
      
c          write(*,*) "etau", etau(alpha,beta)
        enddo !beta
      
      enddo !alpha
c
      return
      end subroutine
c
c
c           Derivative of hardening fn wrt strain
      subroutine mm10_ed_mrr(props, np1, n, stress, tt, ed)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(12) :: tt
      double precision, dimension(6,12) :: ed
c
      ed = 0.d0
c
      return
      end subroutine
c
c           Derivative of sliprate wrt resolved shear stress
      subroutine mm10_dgdt_mrr(props, np1, n, stress, tt, dgammadtau)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: rs
      double precision, dimension(props%nslip) :: dgammadtau
      double precision, dimension(12) :: tt
      double precision :: h, slipinc, mm10_rs
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, dfract,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m,
     &  dslipinc, slipexp
      double precision, dimension(12,12) :: Gmat, Hmat
c
c Load some material parameters
        dt = np1%tinc
        k = props%boltzman
        theta = np1%temp
        G = props%mu_0
        b = props%burgers
        c1 = props%u1
        c2 = props%u2
        c3 = props%u3
        p_e = props%p_v
        q_e = props%q_v
        Qslip = props%G_0_v
        v_attack = props%G_0_y
c        
c Compute the shear modulus using Roter's function
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      call mm10_mrr_GH(props,Gmat,Hmat)
c        
      do alpha = 1,12
c        
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
c        
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          rhoF = dot_product(Gmat(alpha,1:12),tt(1:12))
          rhoP = dot_product(Hmat(alpha,1:12),tt(1:12))
c
c Compute one dependency
        gamma_0 = v_attack*k*theta/(c1*c3*G*b*b)*sqrt(rhoP) ! (15)
        tpass = c1*G*b*dsqrt(rhoP) ! (16)
        tcut = Qslip/(c2*c3*b*b)*dsqrt(rhoF) ! (17)
        fract = ((dabs(rs)-tpass)/tcut)

c Evaluate the equation
        dfract = dsign(1.d0,rs)/tcut
        if(fract.gt.1.d0) then
c linear extrapolation past the too-high stress (rs) value
            b = gamma_0
            m = b * (-q_e*(Qslip/k/theta)*(1.d0 - 1.d0)**(q_e-1.d0))
     &         * (- p_e*(1.d0)**(p_e-1.d0)) * dsign(1.d0,rs)/tcut
            dslipinc = dt * m
        elseif(dabs(rs).gt.0.d0) then
        slipexp = dexp (-(Qslip/k/theta)*(1.d0 - fract**p_e)**q_e)
     &          * dsign(1.d0,rs)
        dslipinc = dt * (gamma_0 * slipexp * -(Qslip/k/theta)*q_e*
     &      (1.d0 - fract**p_e)**(q_e-1.d0)
     &      * -p_e*fract**(p_e-1.d0) * dfract) !(14)
        else
            dslipinc = 0.d0
        endif
        
        dgammadtau(alpha) = dslipinc

      enddo
c
      return
      end subroutine
c
c           Derivative of sliprate wrt hardening variables
      subroutine mm10_dgdh_mrr(props, np1, n, stress, tt, dgammadtt)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(12) :: tt
      double precision, dimension(props%nslip,12) :: dgammadtt
      double precision :: mm10_rs, rs
      integer :: alpha, beta
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, dfract,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m,
     &  dslipinc, slipexp, drhoF, drhoP, dgamma_0,
     &  dtcut, dtpass
      double precision, dimension(12,12) :: Gmat, Hmat
c
c Load some material parameters
        dt = np1%tinc
        k = props%boltzman
        theta = np1%temp
        G = props%mu_0
        b = props%burgers
        c1 = props%u1
        c2 = props%u2
        c3 = props%u3
        p_e = props%p_v
        q_e = props%q_v
        Qslip = props%G_0_v
        v_attack = props%G_0_y
c        
c Compute the shear modulus using Roter's function
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      call mm10_mrr_GH(props,Gmat,Hmat)
        
c Compute derivative of slip rate alpha w.r.t. density beta
c loop over slip rate
      do alpha = 1,12
c        
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
c        
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          rhoF = dot_product(Gmat(alpha,1:12),tt(1:12))
          rhoP = dot_product(Hmat(alpha,1:12),tt(1:12))
c          
c Compute one dependency
        gamma_0 = v_attack*k*theta/(c1*c3*G*b*b)*sqrt(rhoP) ! (15)
        tpass = c1*G*b*dsqrt(rhoP) ! (16)
        tcut = Qslip/(c2*c3*b*b)*dsqrt(rhoF) ! (17)
        fract = ((dabs(rs)-tpass)/tcut)
c        
c loop over density
        do beta = 1,12
        
c           [drhoF,drhoP] = mm10_drhoFP_mrr(props, np1, n, tt, alpha, beta);
          drhoF = Gmat(alpha,beta)
          drhoP = Hmat(alpha,beta)
        
          dgamma_0 = 0.5d0*v_attack*k*theta/(c1*c3*G*b*b)
     &              /dsqrt(rhoP)*drhoP ! (15)
          dtpass = 0.5d0*c1*G*b/dsqrt(rhoP)*drhoP ! (16)
          dtcut = 0.5d0*Qslip/(c2*c3*b*b)/dsqrt(rhoF)*drhoF ! (17)
          dfract = ((-dtpass)*tcut - 
     &             (dabs(rs)-tpass)*dtcut)/tcut**2.d0

c Evaluate the equation
          if(fract.gt.1.d0) then
              ! linear extrapolation past the too-high stress (rs) value
              dslipinc = dt * (-q_e*(Qslip/k/theta)*(1.d0 - 1.d0)
     &           **(q_e-1.d0)) * (- p_e*(1.d0)**(p_e-1.d0))
     &            * dsign(1.d0,rs) * (dgamma_0/tcut 
     &             - gamma_0*dtcut/tcut**2.d0)
          elseif(dabs(rs).gt.0.d0) then
          slipexp = dexp (-(Qslip/k/theta)*
     &             (1.d0 - fract**p_e)**q_e) !(14)
          dslipinc = dt * (dgamma_0 * slipexp * dsign(1.d0,rs)
     &          + gamma_0 * dsign(1.d0,rs) * slipexp *
     &          -(Qslip/k/theta)*q_e*(1.d0 - fract**p_e)**(q_e-1.d0)
     &        * -p_e*fract**(p_e-1.d0) * dfract)
          else
              dslipinc = 0.d0
          endif
        
          dgammadtt(alpha,beta) = dslipinc

        enddo !beta
      
      enddo !alpha
c
      return
      end subroutine
c
c           Derivative of sliprate wrt strain increment
      subroutine mm10_dgdd_mrr(props, np1, n, stress, tt, D, 
     &              dgammadd)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, D
      double precision, dimension(6,props%nslip) :: dgammadd
      double precision, dimension(12) :: tt
c
      dgammadd = 0.d0
c
      return
      end subroutine
c
c     Interaction matrices for parallel and forest dislocations
      subroutine mm10_mrr_GH(props,G,H)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      double precision, dimension(12,12) :: G, H

      G = reshape((/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00
     &  /), shape(G))
      H = reshape((/ 
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00
     &    /), shape(H))
c
      return
      end subroutine

