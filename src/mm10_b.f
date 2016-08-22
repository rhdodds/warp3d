c
c ****************************************************************************
c *                                                                          *
c *    mm10.f                                                                *
c *                                                                          *
c *         written by : mcm                                                 *
c *                                                                          *
c *         Slip rate and hardening functions                                *
c *                                                                          *
c ****************************************************************************
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10b_unknown_hard_error          *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 1/27/14                     *
c     *                                                              *
c     *     A common error message for the general hardening setup   *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10b_unknown_hard_error(props)
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
      elseif (props%h_type .eq. 4) then ! ORNL
        call mm10_h_ornl(props, np1, n, vec1, vec2,stress,tt,h,gp)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_h_mrr(props, np1, n, vec1, vec2,stress,tt,h,gp)
      elseif (props%h_type .eq. 9) then ! DJGM
        call mm10_h_DJGM(props, np1, n, vec1, vec2,stress,tt,h)
      else
        call mm10b_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
      R2(1:props%num_hard) = tt(1:props%num_hard) - h(1:props%num_hard)
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
      elseif (props%h_type .eq. 4) then ! ORNL
        call mm10_hi_ornl(props, np1, n, ivec1, ivec2, stress, tt, h)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_hi_mrr(props, np1, n, ivec1, ivec2, stress, tt, h)
      elseif (props%h_type .eq. 9) then ! DJGM
      else
        call mm10b_unknown_hard_error(props)
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
      elseif (props%h_type .eq. 4) then ! ORNL
        call mm10_dgdt_ornl(props,np1, n, stress, tt, dgammadtau)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_dgdt_mrr(props,np1, n, stress, tt, dgammadtau)
           if (debug) write (*,*) "dgammadtau", dgammadtau(1:3)
           if (debug) write (*,*) "stress", stress(1:3)
      elseif (props%h_type .eq. 9) then ! DJGM
        call mm10_dgdt_DJGM(props,np1, n, stress, tt, dgammadtau)
      else
        call mm10b_unknown_hard_error(props)
      endif
c ******* END: Add new Constitutive Models into this block *********
c
      do i=1,props%nslip
        call DGER(6,6,dgammadtau(i),
     &      matmul(props%stiffness, np1%ms(1:6,i))
     &      + 2.0d0*symtqmat(1:6,i), 1, np1%ms(1:6,i),1,J11,6)
c        if (debug) write (*,*) "J11", J11(1,1)
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
      elseif (props%h_type .eq. 4) then ! ORNL
        call mm10_dgdh_ornl(props,np1, n, stress, tt, dgammadtt)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_dgdh_mrr(props,np1, n, stress, tt, dgammadtt)
      elseif (props%h_type .eq. 9) then ! DJGM
        call mm10_dgdh_DJGM(props,np1, n, stress, tt, dgammadtt)
      else
        call mm10b_unknown_hard_error(props)
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
      elseif (props%h_type .eq. 4) then ! ORNL
        call mm10_estress_ornl(props, np1, n, vec1, vec2, arr1, arr2,
     &            stress, tt, estr)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_estress_mrr(props, np1, n, vec1, vec2, arr1, arr2,
     &            stress, tt, estr)
      elseif (props%h_type .eq. 9) then ! DJGM
        call mm10_estress_DJGM(props, np1, n, vec1, vec2, arr1, arr2,
     &            stress, tt, estr)
      else
        call mm10b_unknown_hard_error(props)
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
      elseif (props%h_type .eq. 4) then ! ORNL
        call mm10_ehard_ornl(props, np1, n, vec1, vec2, arr1, arr2,
     &          stress, tt, etau)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_ehard_mrr(props, np1, n, vec1, vec2, arr1, arr2,
     &          stress, tt, etau)
      elseif (props%h_type .eq. 9) then ! DJGM
        call mm10_ehard_DJGM(props, np1, n, vec1, vec2, arr1, arr2,
     &          stress, tt, etau)
      else
        call mm10b_unknown_hard_error(props)
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
      double precision :: zero
c
      zero = 0.d0
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
         vec1 = zero
         vec2 = zero
      elseif (props%h_type .eq. 2) then ! MTS
         vec1 = zero
         vec2 = zero
      elseif (props%h_type .eq. 3) then ! User
         vec1 = zero
         vec2 = zero
      elseif (props%h_type .eq. 4) then ! ORNL
        call mm10_v_ornl(props, np1, n, stress, tt, 
     &   vec1, vec2)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_v_mrr(props, np1, n, stress, tt, 
     &   vec1, vec2)
      elseif (props%h_type .eq. 9) then ! DJGM
        call mm10_v_DJGM(props, np1, n, stress, tt, 
     &   vec1, vec2)
      else
        call mm10b_unknown_hard_error(props)
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
      elseif (props%h_type .eq. 4) then ! ORNL
        call mm10_a_ornl(props, np1, n, stress, tt, vec1, vec2,
     &          arr1, arr2, both)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_a_mrr(props, np1, n, stress, tt, vec1, vec2,
     &          arr1, arr2, both)
      elseif (props%h_type .eq. 9) then ! DJGM
        call mm10_a_DJGM(props, np1, n, stress, tt, vec1, vec2,
     &          arr1, arr2, both)
      else
        call mm10b_unknown_hard_error(props)
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
      elseif (props%h_type .eq. 4) then ! ORNL
        call mm10_vi_ornl(props, np1, n, stress, tt, 
     &   ivec1, ivec2)
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_vi_mrr(props, np1, n, stress, tt, 
     &   ivec1, ivec2)
      else
        call mm10b_unknown_hard_error(props)
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
        call mm10_slipinc(props, np1, n, stress, tt(1), i, slipinc)
        dbar = dbar + slipinc*np1%ms(1:6,i)
      end do
      elseif (props%h_type .eq. 2) then ! MTS
      dbar = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt(1), i, slipinc)
        dbar = dbar + slipinc*np1%ms(1:6,i)
      end do
      elseif (props%h_type .eq. 3) then ! User
      dbar = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc_user(props, np1, n, stress, tt, i, slipinc)
        dbar = dbar + slipinc*np1%ms(1:6,i)
      end do
      elseif (props%h_type .eq. 4) then ! ORNL
c      dbar = 0.0
c      do i=1,props%nslip
c c        call mm10_slipinc_ornl(props, np1, n, stress, tt, i, slipinc)
c        slipinc = vec1(i)
c        dbar = dbar + slipinc*np1%ms(1:6,i)
c      end do
        dbar(1:6) = matmul(np1%ms(1:6,1:props%nslip),
     &                     vec1(1:props%nslip))
      elseif (props%h_type .eq. 7) then ! MRR
c      dbar = 0.0
c      do i=1,props%nslip
c c        call mm10_slipinc_mrr(props, np1, n, stress, tt, i, slipinc)
c        slipinc = vec1(i)
c        dbar = dbar + slipinc*np1%ms(1:6,i)
c      end do
        dbar(1:6) = matmul(np1%ms(1:6,1:props%nslip),
     &                     vec1(1:props%nslip))
      elseif (props%h_type .eq. 9) then ! DJGM
c      dbar = 0.0
c      do i=1,props%nslip
c c        call mm10_slipinc_DJGM(props, np1, n, stress, tt, i, slipinc)
c        slipinc = vec1(i)
c        dbar = dbar + slipinc*np1%ms(1:6,i)
c      end do
        dbar(1:6) = matmul(np1%ms(1:6,1:props%nslip),
     &                     vec1(1:props%nslip))
      else
        call mm10b_unknown_hard_error(props)
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
      elseif (props%h_type .eq. 4) then ! ORNL
c      wbar = 0.0
c      do i=1,props%nslip
c c        call mm10_slipinc_ornl(props, np1, n, stress, tt, i, slipinc)
c        slipinc = vec1(i)
c        wbar = wbar + slipinc*np1%qs(1:3,i)
c      end do
        wbar(1:3) = matmul(np1%qs(1:3,1:props%nslip),
     &                     vec1(1:props%nslip))
      elseif (props%h_type .eq. 7) then ! MRR
c      wbar = 0.0
c      do i=1,props%nslip
c c        call mm10_slipinc_mrr(props, np1, n, stress, tt, i, slipinc)
c        slipinc = vec1(i)
c        wbar = wbar + slipinc*np1%qs(1:3,i)
c      end do
        wbar(1:3) = matmul(np1%qs(1:3,1:props%nslip),
     &                     vec1(1:props%nslip))
      elseif (props%h_type .eq. 9) then ! DJGM
c      wbar = 0.0
c      do i=1,props%nslip
c c        call mm10_slipinc_DJGM(props, np1, n, stress, tt, i, slipinc)
c        slipinc = vec1(i)
c        wbar = wbar + slipinc*np1%qs(1:3,i)
c      end do
        wbar(1:3) = matmul(np1%qs(1:3,1:props%nslip),
     &                     vec1(1:props%nslip))
      else
        call mm10b_unknown_hard_error(props)
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
      elseif (props%h_type .eq. 4) then ! ornl
      w = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc_ornl(props, np1, n, stress, tt, i, slipinc)
        w = w + slipinc*np1%qc(1:3,i)
      end do
      elseif (props%h_type .eq. 7) then ! MRR
      w = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc_mrr(props, np1, n, stress, tt, i, slipinc)
        w = w + slipinc*np1%qc(1:3,i)
      end do
      elseif (props%h_type .eq. 9) then ! DJGM
      w = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc_DJGM(props, np1, n, stress, tt, i, slipinc)
        w = w + slipinc*np1%qc(1:3,i)
      end do
      else
        call mm10b_unknown_hard_error(props)
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
      elseif (props%h_type .eq. 4) then ! ORNL
      dbar = (0.d0,0.d0)
      do i=1,props%nslip
c        call mm10_slipinc_ornl(props, np1, n, stress, tt, i, slipinc)
        slipinc = ivec1(i)
        dbar = dbar + slipinc*np1%ms(1:6,i)
      end do
      elseif (props%h_type .eq. 7) then ! MRR
      dbar = (0.d0,0.d0)
      do i=1,props%nslip
c        call mm10_slipinc_mrr(props, np1, n, stress, tt, i, slipinc)
        slipinc = ivec1(i)
        dbar = dbar + slipinc*np1%ms(1:6,i)
      end do
      else
        call mm10b_unknown_hard_error(props)
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
      elseif (props%h_type .eq. 4) then ! ORNL
      w = (0.d0,0.d0)
      do i=1,props%nslip
        call mm10_slipinci_ornl(props, np1, n, stress, tt, i, slipinc)
        w = w + slipinc*np1%qc(1:3,i)
      end do
      elseif (props%h_type .eq. 7) then ! MRR
      w = (0.d0,0.d0)
      do i=1,props%nslip
        call mm10_slipinci_mrr(props, np1, n, stress, tt, i, slipinc)
        w = w + slipinc*np1%qc(1:3,i)
      end do
      else
        call mm10b_unknown_hard_error(props)
      end if
c ******* END: Add new Constitutive Models into this block *********
c
      return
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
c
c
c *****************************************************************************
c *                                                                           *
c *         User hardening routines                                           *
c *                                                                           *
c *****************************************************************************
c
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
c      write(*,*) stress
c      write(*,*) i
c      write(*,*) tt
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
c           Form intermediate vectors for faster calculations
      subroutine mm10_v_mrr(props, np1, n, stress, tt, 
     &   vec1, vec2)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt, h
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
      double precision, dimension(props%num_hard) :: tt, h
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
      double complex, dimension(props%num_hard) :: tt, h
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
      double precision, dimension(props%num_hard) :: tt
      double precision :: h, slipinc, mm10_rs
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
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
        if(G.lt.0.d0) then
        G = -G
        else
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
        endif
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
c        
c      ms = np1.ms(1:6,alpha);
c      rs = stress*ms; % tau^a
      rs = mm10_rs(props, np1, n, stress, tt, alpha)
c        
c         [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          rhoF = dot_product(props%Gmat(alpha,1:props%num_hard),
     &   tt(1:props%num_hard))
          rhoP = dot_product(props%Hmat(alpha,1:props%num_hard),
     &   tt(1:props%num_hard))
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

            call splunj
        if(fract.gt.1.d0) then
            ! linear extrapolation past the too-high stress (rs) value
c          write(*,*) 'rs too high', fract, rs, 'rhoP', rhoP 
c     &   , tpass, alpha, c1, G
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
      double precision, dimension(props%num_hard) :: tt, h,
     &       rhoFs,rhoPs
      double precision, dimension(max_uhard) :: vec1, vec2
      integer :: alpha, gp
      logical :: mat_debug
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, rho,
     &  rho_n, pi, c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  tem1, tem2, tem3, tem4
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat

       mat_debug = .false.

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
        if(G.lt.0.d0) then
        G = -G
        else
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
        endif
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Gmat,props%num_hard,tt,1,0.d0,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Hmat,props%num_hard,tt,1,0.d0,rhoPs,1)
c
c
c      write(*,*) "Gmat", Gmat(1,1)
c      write(*,*) "G", G
      do alpha = 1,props%num_hard
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
c          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
c          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
          rhoF = rhoFs(alpha)
          rhoP = rhoPs(alpha)
        if(mat_debug) then
        write(*,*) "rhoF", rhoF
        write(*,*) "rhoP", rhoP
        endif
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
      double complex, dimension(props%num_hard) :: tt, temp
      double complex :: h, slipinc, mm10_rsi
      integer :: alpha, i
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44
      double complex :: rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
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
        if(G.lt.0.d0) then
        G = -G
        else
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
        endif
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
c        
c      ms = np1.ms(1:6,alpha);
c      rs = stress*ms; % tau^a
      rs = mm10_rsi(props, np1, n, stress, tt, alpha)
c        
c         [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          temp = (props%Gmat(alpha,1:props%num_hard)*
     &     tt(1:props%num_hard))
          rhoF = sum(temp)
          temp = (props%Hmat(alpha,1:props%num_hard)*
     &     tt(1:props%num_hard))
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
      double complex, dimension(props%num_hard) :: tt, h, temp
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
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
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
        if(G.lt.0.d0) then
        G = -G
        else
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
        endif
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
c
c
c      write(*,*) "Gmat", Gmat(1,1)
c      write(*,*) "G", G
      do alpha = 1,props%num_hard
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
          temp = (props%Gmat(alpha,1:props%num_hard)
     &      *tt(1:props%num_hard))
          rhoF = sum(temp)
          temp = (props%Hmat(alpha,1:props%num_hard)
     &      *tt(1:props%num_hard))
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
      double precision, dimension(props%num_hard) :: tt, zerosV
      double precision :: h, slipinc, mm10_rs
      integer :: alpha
c
      double complex, dimension(6) :: stressi
      double complex, dimension(props%num_hard) :: tti
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
      double precision, dimension(props%num_hard) :: tt,rhoPs,rhoFs
      double precision, dimension(props%num_hard,6) :: et
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, rho,
     &  rho_n, pi, c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  dddipole, dslipinc, badterm
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
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
        if(G.lt.0.d0) then
        G = -G
        else
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
        endif
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Gmat,props%num_hard,tt,1,0.d0,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Hmat,props%num_hard,tt,1,0.d0,rhoPs,1)
c
      do alpha = 1,props%num_hard

          ! Get dislocation density
          rho = tt(alpha) ! rho^a_SSD
c
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
c
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
c          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
c          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
          rhoF = rhoFs(alpha)
          rhoP = rhoPs(alpha)
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
      double precision, dimension(props%num_hard) :: tt,rhoPs,rhoFs
      double precision, dimension(props%num_hard,props%num_hard) :: etau
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
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
      double precision, dimension(props%num_hard,props%num_hard)
     &    :: dslip
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
        if(G.lt.0.d0) then
        G = -G
        else
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
        endif
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Gmat,props%num_hard,tt,1,0.d0,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Hmat,props%num_hard,tt,1,0.d0,rhoPs,1)
c
c Compute drho_alpha/drho_beta
c loop over numerator hardening variable
      do alpha = 1,props%num_hard

c Get dislocation density
        rho = tt(alpha) ! rho^a_SSD
c
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
          
c         [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
c          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
c          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
          rhoF = rhoFs(alpha)
          rhoP = rhoPs(alpha)
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
        do beta = 1,props%num_hard
c          
c           [drhoF,drhoP] = mm10_drhoFP_mrr(props, np1, n, tt, alpha, beta);
          drhoF = props%Gmat(alpha,beta)
          drhoP = props%Hmat(alpha,beta)
          
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
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(6,props%num_hard) :: ed
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
      double precision, dimension(props%num_hard) :: tt,rhoFs,rhoPs
      double precision :: h, slipinc, mm10_rs
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, dfract,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m,
     &  dslipinc, slipexp
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
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
        if(G.lt.0.d0) then
        G = -G
        else
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
        endif
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Gmat,props%num_hard,tt,1,0.d0,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Hmat,props%num_hard,tt,1,0.d0,rhoPs,1)
c        
      do alpha = 1,props%num_hard
c        
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
c        
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
c          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
c          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
          rhoF = rhoFs(alpha)
          rhoP = rhoPs(alpha)
c
c Compute one dependency
        gamma_0 = v_attack*k*theta/(c1*c3*G*b*b)*sqrt(rhoP) ! (15)
        tpass = c1*G*b*dsqrt(rhoP) ! (16)
        tcut = Qslip/(c2*c3*b*b)*dsqrt(rhoF) ! (17)
        fract = ((dabs(rs)-tpass)/tcut)

c Evaluate the equation
        dfract = dsign(1.d0,rs)/tcut
        if(fract.gt.1.d0) then
c            write(*,*) G, b, theta, c1, rhoP
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
      double precision, dimension(props%num_hard) :: tt,rhoFs,rhoPs
      double precision, dimension(props%nslip,props%num_hard)
     &    :: dgammadtt
      double precision :: mm10_rs, rs
      integer :: alpha, beta
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, dfract,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m,
     &  dslipinc, slipexp, drhoF, drhoP, dgamma_0,
     &  dtcut, dtpass
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
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
        if(G.lt.0.d0) then
        G = -G
        else
        K11=123.323d0+6.7008d-8*theta**3.d0
     &     -1.1342d-4*theta**2.d0-7.8788d-3*theta
        K12=70.6512d0+4.4105d-8*theta**3.d0
     &     -7.5498d-5*theta**2.d0-3.9992d-3*theta
        K44=31.2071d0+7.0477d-9*theta**3.d0
     &     -1.2136d-5*theta**2.d0-8.3274d-3*theta
        G = 1.d0/3.d0*(K11-K12+K44)*1d9
        endif
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Gmat,props%num_hard,tt,1,0.d0,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Hmat,props%num_hard,tt,1,0.d0,rhoPs,1)
        
c Compute derivative of slip rate alpha w.r.t. density beta
c loop over slip rate
      do alpha = 1,props%num_hard
c        
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
c        
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
c          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
c          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
          rhoF = rhoFs(alpha)
          rhoP = rhoPs(alpha)
c          
c Compute one dependency
        gamma_0 = v_attack*k*theta/(c1*c3*G*b*b)*sqrt(rhoP) ! (15)
        tpass = c1*G*b*dsqrt(rhoP) ! (16)
        tcut = Qslip/(c2*c3*b*b)*dsqrt(rhoF) ! (17)
        fract = ((dabs(rs)-tpass)/tcut)
c        
c loop over density
        do beta = 1,props%num_hard
        
c           [drhoF,drhoP] = mm10_drhoFP_mrr(props, np1, n, tt, alpha, beta);
          drhoF = props%Gmat(alpha,beta)
          drhoP = props%Hmat(alpha,beta)
        
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
      double precision, dimension(props%num_hard) :: tt
c
      dgammadd = 0.d0
c
      return
      end subroutine
c
c *****************************************************************************
c *                                                                           *
c *         ORNL ferric-martensitic steel hardening routines                  *
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
c           Form intermediate vectors for faster calculations
      subroutine mm10_v_ornl(props, np1, n, stress, tt, 
     &   vec1, vec2)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt, h
      double precision, dimension(max_uhard) :: vec1, vec2
      integer :: alpha
      double precision :: slipinc
c
      do alpha = 1,props%nslip
         call mm10_slipinc_ornl(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
         vec1(alpha) = slipinc
      enddo
c
      end subroutine
c
c           Form intermediate arrays for faster calculations
      subroutine mm10_a_ornl(props, np1, n, stress, tt, 
     &   vec1, vec2, arr1, arr2, both)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt, h
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
      integer both
c
         call mm10_dgdt_ornl(props, np1, n, stress, tt, 
     &                            arr1(1:props%num_hard,1))
      if(both.eq.2) then
         call mm10_dgdh_ornl(props, np1, n, stress, tt, 
     &               arr2(1:props%nslip,1:props%num_hard))
      endif
c
      end subroutine
c
c           Form intermediate vectors for faster calculations
      subroutine mm10_vi_ornl(props, np1, n, stress, tt, 
     &   ivec1, ivec2)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double complex, dimension(6) :: stress
      double complex, dimension(props%num_hard) :: tt, h
      double complex, dimension(max_uhard) :: ivec1, ivec2
      integer :: alpha
      double complex :: slipinc
c
      do alpha = 1,props%nslip
         call mm10_slipinci_ornl(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
         ivec1(alpha) = slipinc
      enddo
c
      end subroutine
c
c           Actual ornl sliprate function
      subroutine mm10_slipinc_ornl(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      double precision :: h, slipinc, mm10_rs
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, tau0, 
     &  p_e, q_e, Qslip, v_s, K11, K12, K44, rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m,
     &  fM, lamda, G0, rhoM
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
c
c Load some material parameters
        dt = np1%tinc
        k = props%boltzman
        theta = np1%temp
        G0 = props%mu_0
        b = props%burgers
        c1 = props%u1
        c2 = props%u2
        tau0 = props%u3
        p_e = props%p_v
        q_e = props%q_v
        Qslip = props%G_0_v
        v_s = props%G_0_y
        if( v_s.lt.0.d0 ) then
            v_s = exp(-v_s)
        endif
        fM = 0.1d0
        lamda = c2*b
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - 1.d0)
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
c        
c      ms = np1.ms(1:6,alpha);
c      rs = stress*ms; % tau^a
      rs = mm10_rs(props, np1, n, stress, tt, alpha)
c        
c         [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
c          rhoF = dot_product(props%Gmat(alpha,1:props%num_hard),
c     &   tt(1:props%num_hard))
       rhoP = dot_product(props%Hmat(alpha,1:props%num_hard),
     &   tt(1:props%num_hard))
c
c Compute some stresses and rates
        rhoM = fM*tt(alpha)
        gamma_0 = rhoM*b*v_s/b*lamda
        tpass = c1*G*b*dsqrt(rhoP) ! (16)
        tcut = tau0*G/G0
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
            m = b * (-q_e*(Qslip/k/theta)*(1.d0))
     &         * (- p_e*(1.d0)**(p_e-1.d0)) * dsign(1.d0,rs)/tcut
            y = m*x + b
            slipinc = dt * y
        elseif(dabs(rs).eq.0.d0) then
            slipinc = 0.d0
        else
          if(fract.lt.0.d0) then
            p_e = 1.d0 ! deal with low stresses by allowing the fraction to be small
          endif
          slipinc = dt * gamma_0 * dsign(1.d0,rs)
     & * dexp (-(Qslip/k/theta)*(1.d0 - fract**p_e)**q_e) !(14)x
        endif
c       if(np1.step.gt.21) then
c       write(*,*) alpha, "rs", rs, "fract", fract, "slipinc", slipinc
c       endif
c
      return
      end subroutine
c
c           Actual ornl hardening function
      subroutine mm10_h_ornl(props, np1, n, vec1, vec2, 
     & stress, tt, h,gp)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt, h,
     &   rhoFs,rhoPs
      double precision, dimension(max_uhard) :: vec1, vec2
      integer :: alpha, gp
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, rho,
     &  rho_n, pi, c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  tem1, tem2, tem3, tem4, G0
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
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
        G0 = props%mu_0
        b = props%burgers
        v = 0.3d0 !props%nu
        dt = np1%tinc
        PI=4.D0*DATAN(1.D0)
c      write(*,*) "pi", pi
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - 1.d0)
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Gmat,props%num_hard,tt,1,0.d0,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Hmat,props%num_hard,tt,1,0.d0,rhoPs,1)
c
c
c      write(*,*) "Gmat", Gmat(1,1)
c      write(*,*) "G", G
      do alpha = 1,props%num_hard
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
c          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
c          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
          rhoF = rhoFs(alpha)
          rhoP = rhoPs(alpha)
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
c          call mm10_slipinc_ornl(props, np1, n, stress, tt, alpha, 
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
c           Imaginary ornl sliprate function
      subroutine mm10_slipinci_ornl(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double complex, dimension(6) :: stress
      double complex, dimension(props%num_hard) :: tt, temp
      double complex :: h, slipinc, mm10_rsi
      integer :: alpha, i
c
      double precision :: dt, k, theta, G, b, c1, c2, tau0, 
     &  p_e, q_e, Qslip, v_s, K11, K12, K44
      double complex :: rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m,
     &  fM, lamda, G0, rhoM
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
c
c Load some material parameters
        dt = np1%tinc
        k = props%boltzman
        theta = np1%temp
        G0 = props%mu_0
        b = props%burgers
        c1 = props%u1
        c2 = props%u2
        tau0 = props%u3
        p_e = props%p_v
        q_e = props%q_v
        Qslip = props%G_0_v
        v_s = props%G_0_y
        if( v_s.lt.0.d0 ) then
            v_s = exp(-v_s)
        endif
        fM = 0.1d0
        lamda = c2*b
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - 1.d0)
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
c        
c      ms = np1.ms(1:6,alpha);
c      rs = stress*ms; % tau^a
      rs = mm10_rsi(props, np1, n, stress, tt, alpha)
c        
c         [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          temp = (props%Hmat(alpha,1:props%num_hard)*
     &     tt(1:props%num_hard))
          rhoP = sum(temp)
c
c Compute some stresses and rates
        rhoM = fM*tt(alpha)
        gamma_0 = rhoM*b*v_s/b*lamda
        tpass = c1*G*b*cdsqrt(rhoP) ! (16)
        tcut = tau0*G/G0
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
            m = b * (-q_e*(Qslip/k/theta)*(1.d0))
     &         * (- p_e*(1.d0)**(p_e-1.d0)) * dsign(1.d0,dreal(rs))/tcut
            y = m*x + b
            slipinc = dt * y
        elseif(dabs(dreal(rs)).eq.0.d0) then
            slipinc = 0.d0
        else
          if(dreal(fract).lt.0.d0) then
            p_e = 1.d0 ! deal with low stresses by allowing the fraction to be small
          endif
          slipinc = dt * gamma_0 * dsign(1.d0,dreal(rs))
     & * cdexp (-(Qslip/k/theta)*(1.d0 - fract**p_e)**q_e) !(14)
        endif
c
      return
      end subroutine
c
c           Imaginary ornl hardening function
      subroutine mm10_hi_ornl(props, np1, n, ivec1, ivec2, 
     & stress, tt, h)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double complex, dimension(6) :: stress
      double complex, dimension(props%num_hard) :: tt, h, temp
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
     &  tem1, tem2, tem3, G0
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
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
        G0 = props%mu_0
        b = props%burgers
        v = 0.3d0 !props%nu
        dt = np1%tinc
        PI=4.D0*DATAN(1.D0)
c      write(*,*) "pi", pi
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - 1.d0)
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
c
c
c      write(*,*) "Gmat", Gmat(1,1)
c      write(*,*) "G", G
      do alpha = 1,props%num_hard
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
          temp = (props%Gmat(alpha,1:props%num_hard)
     &      *tt(1:props%num_hard))
          rhoF = sum(temp)
          temp = (props%Hmat(alpha,1:props%num_hard)
     &      *tt(1:props%num_hard))
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
c          call mm10_slipinc_ornl(props, np1, n, stress, tt, alpha, 
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
c           Wrapper version, ornl slipinc function
      subroutine mm10_slipinc_ornlW(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt, zerosV
      double precision :: h, slipinc, mm10_rs
      integer :: alpha
c
      double complex, dimension(6) :: stressi
      double complex, dimension(props%num_hard) :: tti
      double complex :: slipinci
c
      zerosV = 0.d0
      stressi = dcmplx(stress,zerosV(1:6))
      tti = dcmplx(tt,zerosV)
c
      call mm10_slipinci_ornl(props, np1, n, stressi, tti, 
     &                            alpha, slipinci)
c
      slipinc = dreal(slipinci)
c
      return
      end subroutine
c
c           Derivative of hardening fn wrt stress
      subroutine mm10_estress_ornl(props, np1, n, vec1, vec2, 
     &        arr1, arr2, stress, tt, et)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, dtdstress
      double precision, dimension(props%num_hard) :: tt,rhoFs,rhoPs
      double precision, dimension(props%num_hard,6) :: et
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 

     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, rho,
     &  rho_n, pi, c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  dddipole, dslipinc, badterm, G0
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
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
        G0 = props%mu_0
        b = props%burgers
        v = 0.3d0 !props%nu
        dt = np1%tinc
        PI=4.D0*DATAN(1.D0)
c
c compute derivatives of slip increments with respect to resolved
c shear stress
c        call mm10_dgdt_ornl(props, np1, n, stress, 
c     &         tt, dslip)
        dslip(1:props%num_hard) = arr1(1:props%num_hard,1)
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - 1.d0)
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Gmat,props%num_hard,tt,1,0.d0,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Hmat,props%num_hard,tt,1,0.d0,rhoPs,1)
c
      do alpha = 1,props%num_hard

          ! Get dislocation density
          rho = tt(alpha) ! rho^a_SSD
c
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
c
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
c          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
c          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
          rhoF = rhoFs(alpha)
          rhoP = rhoPs(alpha)
c
          tpass = c1*G*b*dsqrt(rhoP) ! (16)
          ddipole = dsqrt(3.d0)*G*b/(16.d0*pi*(1.d0-v))/
     &       (dabs(rs)) ! (42)
          dddipole = -dsqrt(3.d0)*G*b/(16.d0*pi*(1.d0-v))/
     &          (dabs(rs))**2.d0*dsign(1.d0,rs)
          rhoM = (2.d0*k/(c1*c2*c3*G*b**3.d0))*
     &       theta*dsqrt(rhoF*rhoP) ! (13)
c          
c          call mm10_slipinc_ornl(props, np1, n, stress, tt, alpha, 
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
      subroutine mm10_ehard_ornl(props, np1, n, vec1, vec2,
     &       arr1, arr2, stress, tt, etau)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt,rhoPs,rhoFs
      double precision, dimension(props%num_hard,props%num_hard) :: etau
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, rho,
     &  rho_n, pi, c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  dddipole, dslipinc, badterm, deltaij,
     &  drhoF, drhoP, drhoM,G0
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
      double precision, dimension(props%num_hard,props%num_hard)
     &    :: dslip
      integer :: alpha, beta
c      
c compute derivatives of slip increments with respect to densities
c        call mm10_dgdh_ornl(props, np1, n, stress, tt, dslip)
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
        G0 = props%mu_0
        b = props%burgers
        v = 0.3d0 !props%nu
        dt = np1%tinc
        PI=4.D0*DATAN(1.D0)
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - 1.d0)
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Gmat,props%num_hard,tt,1,0.d0,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Hmat,props%num_hard,tt,1,0.d0,rhoPs,1)

c
c Compute drho_alpha/drho_beta
c loop over numerator hardening variable
      do alpha = 1,props%num_hard

c Get dislocation density
        rho = tt(alpha) ! rho^a_SSD
c
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
          
c         [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
c          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
c          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
          rhoF = rhoFs(alpha)
          rhoP = rhoPs(alpha)
c
c          call mm10_slipinc_ornl(props, np1, n, stress, tt, alpha, 
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
        do beta = 1,props%num_hard
c          
c           [drhoF,drhoP] = mm10_drhoFP_ornl(props, np1, n, tt, alpha, beta);
          drhoF = props%Gmat(alpha,beta)
          drhoP = props%Hmat(alpha,beta)
          
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
      subroutine mm10_ed_ornl(props, np1, n, stress, tt, ed)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(6,props%num_hard) :: ed
c
      ed = 0.d0
c
      return
      end subroutine
c
c           Derivative of sliprate wrt resolved shear stress
      subroutine mm10_dgdt_ornl(props, np1, n, stress, tt, dgammadtau)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: rs
      double precision, dimension(props%nslip) :: dgammadtau
      double precision, dimension(props%num_hard) :: tt,rhoFs,rhoPs
      double precision :: h, slipinc, mm10_rs
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, tau0, 
     &  p_e, q_e, Qslip, v_s, K11, K12, K44, dfract,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m,
     &  dslipinc, slipexp, G0, fM, lamda, rhoM, p2
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
c
c Load some material parameters
        dt = np1%tinc
        k = props%boltzman
        theta = np1%temp
        G0 = props%mu_0
        b = props%burgers
        c1 = props%u1
        c2 = props%u2
        tau0 = props%u3
        p_e = props%p_v
        q_e = props%q_v
        Qslip = props%G_0_v
        v_s = props%G_0_y
        if( v_s.lt.0.d0 ) then
            v_s = exp(-v_s)
        endif
        fM = 0.1d0
        lamda = c2*b
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - 1.d0)
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Gmat,props%num_hard,tt,1,0.d0,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Hmat,props%num_hard,tt,1,0.d0,rhoPs,1)
c        
      do alpha = 1,props%num_hard
c        
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
c        
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
c          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
c          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
          rhoF = rhoFs(alpha)
          rhoP = rhoPs(alpha)
c
c Compute one dependency
        rhoM = fM*tt(alpha)
        gamma_0 = rhoM*b*v_s/b*lamda ! (15)
        tpass = c1*G*b*dsqrt(rhoP) ! (16)
        tcut = tau0*G/G0 ! (17)
        fract = ((dabs(rs)-tpass)/tcut)

c Evaluate the equation
        dfract = dsign(1.d0,rs)/tcut
        if(fract.gt.1.d0) then
c linear extrapolation past the too-high stress (rs) value
            b = gamma_0
            m = b * (-q_e*(Qslip/k/theta)*(1.d0))
     &         * (- p_e*(1.d0)**(p_e-1.d0)) * dsign(1.d0,rs)/tcut
            dslipinc = dt * m
        elseif(dabs(rs).eq.0.d0) then
            dslipinc = 0.d0
        else
          if(fract.lt.0.d0) then
            p2 = 1.d0 ! deal with low stresses by allowing the fraction to be small
          else
            p2 = p_e
          endif
        slipexp = dexp (-(Qslip/k/theta)*(1.d0 - fract**p2)**q_e)
     &          * dsign(1.d0,rs)
        dslipinc = dt * (gamma_0 * slipexp * -(Qslip/k/theta)*q_e*
     &      (1.d0 - fract**p2)**(q_e-1.d0)
     &      * -p2*fract**(p2-1.d0) * dfract) !(14)
        endif
        
        dgammadtau(alpha) = dslipinc

      enddo
c
      return
      end subroutine
c
c           Derivative of sliprate wrt hardening variables
      subroutine mm10_dgdh_ornl(props, np1, n, stress, tt, dgammadtt)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt,rhoPs,rhoFs
      double precision, dimension(props%nslip,props%num_hard)
     &    :: dgammadtt
      double precision :: mm10_rs, rs
      integer :: alpha, beta
c
      double precision :: dt, k, theta, G, b, c1, c2, tau0, 
     &  p_e, q_e, Qslip, v_s, K11, K12, K44, dfract,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m,
     &  dslipinc, slipexp, drhoF, drhoP, dgamma_0,
     &  dtcut, dtpass, G0, fM, rhoM, lamda, drhoM, p2
       !double precision, dimension(props%num_hard,props%num_hard)
       !&   :: Gmat, Hmat
c
c Load some material parameters
        dt = np1%tinc
        k = props%boltzman
        theta = np1%temp
        G0 = props%mu_0
        b = props%burgers
        c1 = props%u1
        c2 = props%u2
        tau0 = props%u3
        p_e = props%p_v
        q_e = props%q_v
        Qslip = props%G_0_v
        v_s = props%G_0_y
        if( v_s.lt.0.d0 ) then
            v_s = exp(-v_s)
        endif
        fM = 0.1d0
        lamda = c2*b
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - 1.d0)
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Gmat,props%num_hard,tt,1,0.d0,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,1.d0,
     &           props%Hmat,props%num_hard,tt,1,0.d0,rhoPs,1)
        
c Compute derivative of slip rate alpha w.r.t. density beta
c loop over slip rate
      do alpha = 1,props%num_hard
c        
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
c        
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
c          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
c          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
c     &    tt(1:props%num_hard))
          rhoF = rhoFs(alpha)
          rhoP = rhoPs(alpha)
c          
c Compute one dependency
        rhoM = fM*tt(alpha)
        gamma_0 = rhoM*b*v_s/b*lamda ! (15)
        tpass = c1*G*b*dsqrt(rhoP) ! (16)
        tcut = tau0*G/G0 ! (17)

        fract = ((dabs(rs)-tpass)/tcut)
c        
c loop over density
        do beta = 1,props%num_hard
        
c           [drhoF,drhoP] = mm10_drhoFP_mrr(props, np1, n, tt, alpha, beta);
          drhoF = props%Gmat(alpha,beta)
          drhoP = props%Hmat(alpha,beta)
          if(alpha.eq.beta) then
             drhoM = fM
          else
             drhoM = 0.d0
          endif
        
          dgamma_0 = drhoM*b*v_s/b*lamda ! (15)
          dtpass = 0.5d0*c1*G*b/dsqrt(rhoP)*drhoP ! (16)
          dtcut = 0.d0 ! (17)
          dfract = ((-dtpass)*tcut - 
     &             (dabs(rs)-tpass)*dtcut)/tcut**2

c Evaluate the equation
          if(fract.gt.1.d0) then
              ! linear extrapolation past the too-high stress (rs) value
              dslipinc = dt * (-q_e*(Qslip/k/theta)*(1.d0)
     &            * (- p_e*(1.d0)**(p_e-1.d0)))
     &            * dsign(1.d0,rs) * (dgamma_0/tcut 
     &             - gamma_0*dtcut/tcut**2)
          elseif(dabs(rs).eq.0.d0) then
              dslipinc = 0.d0
          else
            if(fract.lt.0.d0) then
              p2 = 1.d0
            else
              p2 = p_e
            endif
          slipexp = dexp (-(Qslip/k/theta)*
     &             (1.d0 - fract**p2)**q_e) !(14)
          dslipinc = dt * (dgamma_0 * slipexp * dsign(1.d0,rs)
     &          + gamma_0 * dsign(1.d0,rs) * slipexp *
     &          -(Qslip/k/theta)*q_e*(1.d0 - fract**p2)**(q_e-1.d0)
     &        * -p2*fract**(p2-1.d0) * dfract)
          endif
        
          dgammadtt(alpha,beta) = dslipinc

        enddo !beta
      
      enddo !alpha
c
      return
      end subroutine
c
c           Derivative of sliprate wrt strain increment
      subroutine mm10_dgdd_ornl(props, np1, n, stress, tt, D, 
     &              dgammadd)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, D
      double precision, dimension(6,props%nslip) :: dgammadd
      double precision, dimension(props%num_hard) :: tt
c
      dgammadd = 0.d0
c
      return
      end subroutine
c
c
c *****************************************************************************
c *                                                                           *
c *         AFRL Ti-6242 high temperature hardening routines                  *
c *                                                                           *
c *****************************************************************************
c
c
c
c
c     DJGM:
c
c           Form intermediate vectors for faster calculations
      subroutine mm10_v_DJGM(props, np1,n, stress, tt,
     &   vec1, vec2)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt, h
      double precision, dimension(max_uhard) :: vec1, vec2
      integer :: alpha
      double precision :: slipinc


      do alpha = 1,props%nslip
         call mm10_slipinc_DJGM(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
         vec1(alpha) = slipinc
      enddo

      return
      end
c
c
c     DJGM:
c
c           Form intermediate arrays for faster calculations
       subroutine mm10_a_DJGM(props, np1, n, stress, tt, 
     &   vec1, vec2, arr1, arr2, both)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt, h
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
      integer both
c
         call mm10_dgdt_DJGM(props, np1, n, stress, tt, 
     &                            arr1(1:props%num_hard,1))
      if(both.eq.2) then
         call mm10_dgdh_DJGM(props, np1, n, stress, tt, 
     &               arr2(1:props%nslip,1:props%num_hard))
      endif
      
      return
      end
c
c           Calculate the slip increment along system i  
c           mrr model
      subroutine mm10_slipinc_DJGM(props, np1, n, stress, tt, 
     &                             alpha, slipinc)
c
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      double precision :: slipinc, mm10_rs
      integer :: alpha
c
      double precision :: dt, tau, g_alpha
      double precision :: h_0_alpha, gamma_dot_tilde,
     &     g_tilde, r_alpha, n_alpha, m_alpha, g_0_alpha
        
        ! Load hard coded material parameters
        ! q=Gmat matrix is loaded at top of mm10
        !h_0_alpha = props%Hmat(1,alpha)
        gamma_dot_tilde = props%Hmat(2,alpha)
        !g_tilde = props%Hmat(3,alpha)
        !r_alpha = props%Hmat(4,alpha)
        !n_alpha = props%Hmat(5,alpha)
        m_alpha = props%Hmat(6,alpha)
        !g_0_alpha = props%Hmat(7,alpha)
        
        tau = mm10_rs(props, np1, n, stress, tt, alpha)
        g_alpha = tt(alpha)
        dt = np1%tinc
        
        ! Equation [4]
        slipinc = dt * gamma_dot_tilde * 
     &  abs(tau/g_alpha)**(1.d0/m_alpha)*dsign(1.d0,tau)
      
      return
      end
c
c
c     Hardening function for DJGM model
c
      subroutine mm10_h_DJGM(props, np1,
     &             n,vec1,vec2, stress, tt, g)
c
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt, g, h
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision :: slipinc, mm10_rs
      integer :: slip_a, slip_b
c
      double precision :: dt, tau, g_alpha
      double precision :: h_0_alpha, gamma_dot_tilde,
     &     g_tilde, r_alpha, n_alpha, m_alpha, g_0_alpha,
     &     g_s, gamma_dot, temp, g_n, g_dot
c Deka, Dhyanjyoti, et al. "Crystal plasticity modeling of deformation and 
c creep in polycrystalline Ti-6242." Metallurgical and materials 
c transactions A 37.5 (2006): 1371-1388.

c evolution of hardening parameter 'g'

        ! Load material parameters
        ! q=Gmat matrix is loaded at top of mm10
        
        dt = np1%tinc
        
        ! Equation [6], g_s equation has to be modified
        do slip_b = 1,props%nslip
            h_0_alpha = props%Hmat(1,slip_b)
            gamma_dot_tilde = props%Hmat(2,slip_b)
            g_tilde = props%Hmat(3,slip_b)
            r_alpha = props%Hmat(4,slip_b)
            n_alpha = props%Hmat(5,slip_b)
            m_alpha = props%Hmat(6,slip_b)
            g_0_alpha = props%Hmat(7,slip_b)
            slipinc = vec1(slip_b)
            gamma_dot = abs(slipinc)/dt
            temp = g_tilde * (gamma_dot/gamma_dot_tilde)**n_alpha
            g_s = max(temp, 2.d0/3.d0*g_0_alpha) ! threshold to ensure no slip during elastic response
            h(slip_b) = h_0_alpha * abs(1.d0-tt(slip_b)/g_s)
     &      **r_alpha * sign(1.d0,1.d0-tt(slip_b)/g_s)
        enddo
        
        ! Equation [5]
        do slip_a = 1,props%nslip
            g_dot = 0.d0
            do slip_b = 1,props%nslip
                slipinc = vec1(slip_b)
                gamma_dot = abs(slipinc)/dt
                g_dot = g_dot + props%Gmat(slip_a, slip_b) * 
     &                         h(slip_b) * abs(gamma_dot)
            enddo
            g_n = n%tau_tilde(slip_a)
            g(slip_a) = g_n + g_dot*dt
        enddo
        
      return
      end
c
c
c     DJGM:
c
c           Derivative of hardening fn wrt stress
      subroutine mm10_estress_DJGM(props,
     & np1, n,vec1,vec2,arr1,arr2, stress, tt, et)
 
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, dtdstress, ms
      double precision, dimension(props%num_hard) :: tt,h,dh
      double precision, dimension(props%num_hard,6) :: et
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, slipinc
      double precision :: h_0_alpha, gamma_dot_tilde,
     &     g_tilde, r_alpha, n_alpha, m_alpha, g_0_alpha,
     &     g_s, gamma_dot, temp, g_n, g_dot, dg_s, dslip
      double precision, dimension(props%nslip) :: dslipinc
      integer :: slip_a, slip_b
      
      et = 0.d0

      ! compute derivatives of slip increments with respect to resolved
      ! shear stress
        dslipinc(1:props%num_hard) = arr1(1:props%num_hard,1)

        ! Load material parameters
        dt = np1%tinc       
        
        ! Equation [6], g_s equation has to be modified
        do slip_b = 1,props%nslip
            h_0_alpha = props%Hmat(1,slip_b)
            gamma_dot_tilde = props%Hmat(2,slip_b)
            g_tilde = props%Hmat(3,slip_b)
            r_alpha = props%Hmat(4,slip_b)
            n_alpha = props%Hmat(5,slip_b)
            m_alpha = props%Hmat(6,slip_b)
            g_0_alpha = props%Hmat(7,slip_b)
            slipinc = vec1(slip_b)
            gamma_dot = abs(slipinc)/dt
            dslip = dslipinc(slip_b)/dt
            temp = g_tilde * (gamma_dot/gamma_dot_tilde)**n_alpha
c            write(*,*) 'slip_b', slip_b, 'dslip',dslip
c            write(*,*) 'gamma_dot', gamma_dot, 'temp',temp
            if( temp .gt. 2.d0/3.d0*g_0_alpha ) then ! threshold to ensure no slip during elastic response
                g_s = temp
                h(slip_b) = h_0_alpha * abs(1.d0-tt(slip_b)/g_s)
     &          **r_alpha * sign(1.d0,1.d0-tt(slip_b)/g_s)
                dg_s = g_tilde * n_alpha / gamma_dot * 
     &              (gamma_dot/gamma_dot_tilde)**n_alpha * dslip
                dh(slip_b) = h_0_alpha * r_alpha * 
     &              abs(1.d0-tt(slip_b)/g_s)**(r_alpha-1.d0)
     &              *(dg_s*tt(slip_b)/g_s**2)
            else
                g_s = 2.d0/3.d0*g_0_alpha
                h(slip_b) = h_0_alpha * abs(1.d0-tt(slip_b)/g_s)
     &          **r_alpha * sign(1.d0,1.d0-tt(slip_b)/g_s)
                dh(slip_b) = 0.d0
            endif
c            write(*,*) 'h(slip_b)', h(slip_b), 'g_s',g_s
c            write(*,*) 'dh(slip_b)', dh(slip_b)
        enddo
        
        ! Equation [5]
        do slip_a = 1,props%nslip
            do slip_b = 1,props%nslip
                dtdstress(1:6) = np1%ms(1:6,slip_b)
                slipinc = vec1(slip_b)
                gamma_dot = abs(slipinc)/dt
                dslip = dslipinc(slip_b)/dt
                et(slip_a,1:6) = et(slip_a,1:6) + dt*
     &          props%Gmat(slip_a, slip_b)*
     &              (dh(slip_b)*abs(gamma_dot) +
     &          h(slip_b)*sign(1.d0,gamma_dot)*dslip)*dtdstress(1:6)
            enddo
        enddo
c
       return
        end
c
c
c     DJGM:
c
c           Derivative of hardening fn wrt hardening
        subroutine mm10_ehard_DJGM(props,
     &      np1, n,vec1,vec2,arr1,arr2, stress, tt, etau)
 
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, dtdstress, ms
      double precision, dimension(props%num_hard) :: tt,h,dh
      double precision, dimension(props%num_hard,props%num_hard) :: etau
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, slipinc, deltaij
      double precision :: h_0_alpha, gamma_dot_tilde,
     &     g_tilde, r_alpha, n_alpha, m_alpha, g_0_alpha,
     &     g_s, gamma_dot, temp, g_n, g_dot, dg_s, dslip
      double precision, dimension(props%nslip,props%num_hard)
     &        :: dslipinc
      integer :: slip_a, slip_b
      
       etau = 0.d0

       dslipinc = arr2(1:props%nslip,1:props%num_hard)

        ! Load material parameters
        dt = np1%tinc
        
        ! Equation [6], g_s equation has to be modified
        do slip_b = 1,props%nslip
        
            h_0_alpha = props%Hmat(1,slip_b)
            gamma_dot_tilde = props%Hmat(2,slip_b)
            g_tilde = props%Hmat(3,slip_b)
            r_alpha = props%Hmat(4,slip_b)
            n_alpha = props%Hmat(5,slip_b)
            m_alpha = props%Hmat(6,slip_b)
            g_0_alpha = props%Hmat(7,slip_b)
            slipinc = vec1(slip_b)
            gamma_dot = abs(slipinc)/dt
            dslip = dslipinc(slip_b,slip_b)/dt
            temp = g_tilde * (gamma_dot/gamma_dot_tilde)**n_alpha
            if( temp .gt. 2.d0/3.d0*g_0_alpha ) then ! threshold to ensure no slip during elastic response
                g_s = temp
                h(slip_b) = h_0_alpha * abs(1.d0-tt(slip_b)/g_s)
     &          **r_alpha * sign(1.d0,1.d0-tt(slip_b)/g_s)
                dg_s = g_tilde * n_alpha / gamma_dot * 
     &              (gamma_dot/gamma_dot_tilde)**n_alpha * dslip
                dh(slip_b) = h_0_alpha * r_alpha * 
     &              abs(1.d0-tt(slip_b)/g_s)**(r_alpha-1.d0)
     &              *(-1.d0/g_s + dg_s*tt(slip_b)/g_s**2)
            else
                g_s = 2.d0/3.d0*g_0_alpha
                h(slip_b) = h_0_alpha * abs(1.d0-tt(slip_b)/g_s)
     &          **r_alpha * sign(1.d0,1.d0-tt(slip_b)/g_s)
                dh(slip_b) = 0.d0
            endif
            
        enddo
        
        ! Equation [5]
        do slip_a = 1,props%nslip
            do slip_b = 1,props%nslip
                  if( slip_a.eq.slip_b ) then
                      deltaij = 1.d0
                  else
                      deltaij = 0.d0
                  endif
                slipinc = vec1(slip_b)
                gamma_dot = abs(slipinc)/dt
                dslip = dslipinc(slip_b,slip_b)/dt
                etau(slip_a,slip_b) = deltaij - dt*
     &          props%Gmat(slip_a, slip_b)*
     &              (dh(slip_b)*abs(gamma_dot) + 
     &               h(slip_b)*sign(1.d0,gamma_dot)*dslip)
            enddo
        enddo

       return
       end
c
c
c           Derivative of hardening fn wrt strain
      subroutine mm10_ed_DJGM(props, np1, n, stress, tt, ed)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(6,props%num_hard) :: ed
c
      ed = 0.d0
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_dgdt_DJGM                    *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 06/10/16                    *
c     *                                                              *
c     *     Calculate partial gamma(s) w.r.t. tau(s) for each slip   *
c     *     system, for use in J11 in material integration. DJGM     *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_dgdt_DJGM(props, np1,
     &      n, stress, tt, dgdt)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt, dgdt
      double precision :: slipinc, mm10_rs
      integer :: slip_a
c
      double precision :: dt, tau, g_alpha
      double precision :: h_0_alpha, gamma_dot_tilde,
     &     g_tilde_alpha, r_alpha, n_alpha, m_alpha, g_0_alpha
        
c Deka, Dhyanjyoti, et al. "Crystal plasticity modeling of deformation and 
c creep in polycrystalline Ti-6242." Metallurgical and materials 
c transactions A 37.5 (2006): 1371-1388.

        ! Load hard coded material parameters
        ! q=Gmat matrix is loaded at top of mm10
        
        dgdt = 0.d0

        dt = np1%tinc
        
        ! Equation [4], slip rate vector
        do slip_a = 1,props%nslip
            gamma_dot_tilde = props%Hmat(2,slip_a)
            m_alpha = props%Hmat(6,slip_a)
            tau = mm10_rs(props, np1, n, stress, tt, slip_a)
            g_alpha = tt(slip_a)
            dgdt(slip_a) = dt * (1.d0/m_alpha) * gamma_dot_tilde
     &      * abs(tau/g_alpha)**((1.d0/m_alpha) - 1.d0) * 
     &        (1.d0/g_alpha)
        enddo
        
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_dgdh_DJGM                    *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 06/10/16                    *
c     *                                                              *
c     *     Calculate partial gamma(s) w.r.t. tt for each slip       *
c     *     system, for use in J12 in material integration. DJGM     *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_dgdh_DJGM(props, np1,
     &      n, stress, tt, dgammadtt)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      double precision, dimension(props%num_hard,props%num_hard)
     &                 :: dgammadtt
      double precision :: slipinc, mm10_rs, dslipinc
      integer :: slip_a, slip_b
c
      double precision :: dt, tau, g_alpha
      double precision :: h_0_alpha, gamma_dot_tilde,
     &     g_tilde_alpha, r_alpha, n_alpha, m_alpha, g_0_alpha
        
        dgammadtt = 0.d0
        
        dt = np1%tinc
        
        ! Equation [4]
        slipinc = dt * gamma_dot_tilde * 
     &  abs(tau/g_alpha)**(1.d0/m_alpha) * dsign(1.d0,tau)

        
        ! Equation [4], slip rate vector
      do slip_a = 1,props%num_hard
        slip_b = slip_a
            gamma_dot_tilde = props%Hmat(2,slip_a)
            m_alpha = props%Hmat(6,slip_a)
            tau = mm10_rs(props, np1, n, stress, tt, slip_a)
            g_alpha = tt(slip_a)
        dslipinc = dt * (1.d0/m_alpha) * gamma_dot_tilde
     &  * abs(tau/g_alpha)**(1.d0/m_alpha)
     &  * (-1.d0/g_alpha)*sign(1.d0,tau)
          dgammadtt(slip_a,slip_b) = dslipinc
      enddo
        
      return
      end
c
c
c           Derivative of sliprate wrt strain increment
      subroutine mm10_dgdd_DJGM(props, np1, n, stress, tt, D, 
     &              dgammadd)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, D
      double precision, dimension(6,props%nslip) :: dgammadd
      double precision, dimension(props%num_hard) :: tt
c
      dgammadd = 0.d0
c
      return
      end subroutine
