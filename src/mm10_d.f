

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
      elseif (props%h_type .eq. 4) then ! ORNL
c        do i=1,props%nslip
c           call mm10_slipinc_ornl(props, np1, n, 
c     &           np1%stress,np1%tau_tilde, i, np1%slip_incs(i))
c        end do
        np1%slip_incs(1:props%nslip)  = vec1(1:props%nslip)
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
      elseif (props%h_type .eq. 4) then ! ORNL
        call mm10_dgdt_ornl(props,np1, n, np1%stress,
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