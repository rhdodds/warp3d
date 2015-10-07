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
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_h_mrr(props, np1, n, vec1, vec2,stress,tt,h,gp)
      else
        call mm10b_unknown_hard_error(props)
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
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_dgdt_mrr(props,np1, n, stress, tt, dgammadtau)
      else
        call mm10b_unknown_hard_error(props)
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
      elseif (props%h_type .eq. 7) then ! User
        call mm10_estress_mrr(props, np1, n, vec1, vec2, arr1, arr2,
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
      elseif (props%h_type .eq. 7) then ! User
        call mm10_ehard_mrr(props, np1, n, vec1, vec2, arr1, arr2,
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
c
c ******* START: Add new Constitutive Models into this block *********
      if (props%h_type .eq. 1) then ! voche
      elseif (props%h_type .eq. 2) then ! MTS
      elseif (props%h_type .eq. 3) then ! User
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_v_mrr(props, np1, n, stress, tt, 
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
      elseif (props%h_type .eq. 7) then ! MRR
        call mm10_a_mrr(props, np1, n, stress, tt, vec1, vec2,
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
      elseif (props%h_type .eq. 7) then ! MRR
      w = 0.0d0
      do i=1,props%nslip
        call mm10_slipinc_mrr(props, np1, n, stress, tt, i, slipinc)
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
      double precision, dimension(props%num_hard,props%num_hard)
     &   :: Gmat, Hmat
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
      call mm10_mrr_GH(props,Gmat,Hmat)
c        
c      ms = np1.ms(1:6,alpha);
c      rs = stress*ms; % tau^a
      rs = mm10_rs(props, np1, n, stress, tt, alpha)
c        
c         [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
     &   tt(1:props%num_hard))
          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
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
      double precision, dimension(props%num_hard) :: tt, h
      double precision, dimension(max_uhard) :: vec1, vec2
      integer :: alpha, gp
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, rho,
     &  rho_n, pi, c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  tem1, tem2, tem3, tem4
      double precision, dimension(props%num_hard,props%num_hard)
     &    :: Gmat, Hmat
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
      call mm10_mrr_GH(props,Gmat,Hmat)
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
          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
     &    tt(1:props%num_hard))
          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
     &    tt(1:props%num_hard))
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
      double precision, dimension(props%num_hard,props%num_hard)
     &      :: Gmat, Hmat
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
      call mm10_mrr_GH(props,Gmat,Hmat)
c        
c      ms = np1.ms(1:6,alpha);
c      rs = stress*ms; % tau^a
      rs = mm10_rsi(props, np1, n, stress, tt, alpha)
c        
c         [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          temp = (Gmat(alpha,1:props%num_hard)*
     &     tt(1:props%num_hard))
          rhoF = sum(temp)
          temp = (Hmat(alpha,1:props%num_hard)*
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
      double precision, dimension(props%num_hard,props%num_hard)
     &     :: Gmat, Hmat
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
      call mm10_mrr_GH(props,Gmat,Hmat)
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
          temp = (Gmat(alpha,1:props%num_hard)
     &      *tt(1:props%num_hard))
          rhoF = sum(temp)
          temp = (Hmat(alpha,1:props%num_hard)
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
      double precision, dimension(props%num_hard) :: tt
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
      double precision, dimension(props%num_hard,props%num_hard)
     &     :: Gmat, Hmat
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
      call mm10_mrr_GH(props,Gmat,Hmat)
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
          rhoF = dot_product(Gmat(alpha,1:props%num_hard)
     &     ,tt(1:props%num_hard))
          rhoP = dot_product(Hmat(alpha,1:props%num_hard)
     &      ,tt(1:props%num_hard))
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
      double precision, dimension(props%num_hard) :: tt
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
      double precision, dimension(props%num_hard,props%num_hard)
     &    :: Gmat, Hmat
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
      call mm10_mrr_GH(props,Gmat,Hmat)
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
          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
     &     tt(1:props%num_hard))
          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
     &     tt(1:props%num_hard))
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
      double precision, dimension(props%num_hard) :: tt
      double precision :: h, slipinc, mm10_rs
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, K11, K12, K44, dfract,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m,
     &  dslipinc, slipexp
      double precision, dimension(props%num_hard,props%num_hard)
     &      :: Gmat, Hmat
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
      call mm10_mrr_GH(props,Gmat,Hmat)
c        
      do alpha = 1,props%num_hard
c        
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
c        
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
     &   tt(1:props%num_hard))
          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
     &   tt(1:props%num_hard))
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
      double precision, dimension(props%num_hard) :: tt
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
      double precision, dimension(props%num_hard,props%num_hard)
     &      :: Gmat, Hmat
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
      call mm10_mrr_GH(props,Gmat,Hmat)
        
c Compute derivative of slip rate alpha w.r.t. density beta
c loop over slip rate
      do alpha = 1,props%num_hard
c        
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rs(props, np1, n, stress, tt, alpha)
c        
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          rhoF = dot_product(Gmat(alpha,1:props%num_hard),
     &      tt(1:props%num_hard))
          rhoP = dot_product(Hmat(alpha,1:props%num_hard),
     &      tt(1:props%num_hard))
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
      double precision, dimension(props%num_hard) :: tt
c
      dgammadd = 0.d0
c
      return
      end subroutine
c
c     Interaction matrices for parallel and forest dislocations, general
      subroutine mm10_mrr_GH(props,G,H)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      double precision, dimension(props%num_hard,props%num_hard)
     &    :: G, H
c
      if(props%s_type.eq.1) then ! FCC
         call mm10_mrr_GHfcc(G,H)
      elseif(props%s_type.eq.2) then ! BCC
         call mm10_mrr_GHbcc(G,H)
      elseif(props%s_type.eq.6) then ! FCC12
         call mm10_mrr_GHfcc12(G,H)
      elseif(props%s_type.eq.7) then ! BCC12
         call mm10_mrr_GHbcc12(G,H)
      elseif(props%s_type.eq.8) then ! BCC48
         call mm10_mrr_GHbcc48(G,H)
      else ! calculate manually
         call mm10_mrr_GHman(props,G,H)
      endif
c
      return
      end subroutine
c
c     Manual calculation of interactions
      subroutine mm10_mrr_GHman(props,G,H)
c
      use mm10_defs
      use crystal_data, only : c_array
      integer :: cnum
c
      type(crystal_props) :: props
      double precision, dimension(props%num_hard,props%num_hard)
     &    :: G, H
      double precision :: temp
      double precision, dimension(3) :: bs,ns
      double precision, dimension(3,props%num_hard) :: tvecs
      integer :: i,s,zeta,xi
c
c     Initialize crystal
c      call initialize_new_crystal(1, 0)
c      c_array(1)%slip_type = props%s_type
c      call finalize_new_crystal(1, 0)
c  ccc The crystal does not have to be re-initialized, it is stored in memory
c     BUT: we need to know which one to access
      cnum = props%cnum ! get crystal ID so that we can find the slip system
c     Now calculate with the slip vectors
c Get tangent vectors
      do i = 1,props%num_hard
         bs = c_array(cnum)%bi(i,:)
         ns = c_array(cnum)%ni(i,:)
         tvecs(1,i) = ns(2)*bs(3) - ns(3)*bs(2)
         tvecs(2,i) = ns(3)*bs(1) - ns(1)*bs(3)
         tvecs(3,i) = ns(1)*bs(2) - ns(2)*bs(1)
      enddo
c Fill in template with vector norms w/o the scaling coefficients
      do zeta = 1,props%num_hard
        do xi = 1,props%num_hard
          temp = c_array(cnum)%ni(xi,1)*tvecs(1,zeta) +
     &           c_array(cnum)%ni(xi,2)*tvecs(2,zeta) +
     &           c_array(cnum)%ni(xi,3)*tvecs(3,zeta)
          G(xi,zeta) = abs(temp)
          temp = 1.d0 - temp*temp
          if(temp.lt.0.d0) temp = 0.d0
          H(xi,zeta) = sqrt(temp)
        enddo
      enddo
c
      return
      end subroutine
c
c     Interaction matrices for FCC system, unitary interactions
      subroutine mm10_mrr_GHfcc(G,H)
c
      double precision, dimension(12,12)
     &    :: G, H
c
      G(1:12,1) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01
     &    /)
      G(1:12,2) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01
     &    /)
      G(1:12,3) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01
     &    /)
      G(1:12,4) = (/
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01
     &    /)
      G(1:12,5) = (/
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01
     &    /)
      G(1:12,6) = (/
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01
     &    /)
      G(1:12,7) = (/
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01
     &    /)
      G(1:12,8) = (/
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01
     &    /)
      G(1:12,9) = (/
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.4280904158206347D-01, +9.4280904158206347D-01
     &    /)
      G(1:12,10) = (/
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00
     &    /)
      G(1:12,11) = (/
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00
     &    /)
      G(1:12,12) = (/
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00
     &    /)
      H(1:12,1) = (/
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01
     &    /)
      H(1:12,2) = (/
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01
     &    /)
      H(1:12,3) = (/
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01
     &    /)
      H(1:12,4) = (/
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01
     &    /)
      H(1:12,5) = (/
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01
     &    /)
      H(1:12,6) = (/
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01
     &    /)
      H(1:12,7) = (/
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01
     &    /)
      H(1:12,8) = (/
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01
     &    /)
      H(1:12,9) = (/
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +3.3333333333333309D-01, +3.3333333333333309D-01
     &    /)
      H(1:12,10) = (/
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00
     &    /)
      H(1:12,11) = (/
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00
     &    /)
      H(1:12,12) = (/
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00
     &    /)
c
      return
      end subroutine
c
c     Interaction matrices for BCC system, unitary interactions
      subroutine mm10_mrr_GHbcc(G,H)
c
      double precision, dimension(12,12)
     &    :: G, H
c
      G(1:12,1) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +2.8867513459481287D-01
     &    /)
      G(1:12,2) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +0.0000000000000000D+00, +2.8867513459481287D-01
     &    /)
      G(1:12,3) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +8.6602540378443860D-01, +5.7735026918962573D-01
     &    /)
      G(1:12,4) = (/
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01
     &    /)
      G(1:12,5) = (/
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +5.7735026918962573D-01, +2.8867513459481287D-01
     &    /)
      G(1:12,6) = (/
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +2.8867513459481287D-01, +5.7735026918962573D-01
     &    /)
      G(1:12,7) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +8.6602540378443860D-01
     &    /)
      G(1:12,8) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +5.7735026918962573D-01, +8.6602540378443860D-01
     &    /)
      G(1:12,9) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +2.8867513459481287D-01, +0.0000000000000000D+00
     &    /)
      G(1:12,10) = (/
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01
     &    /)
      G(1:12,11) = (/
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +0.0000000000000000D+00, +8.6602540378443860D-01
     &    /)
      G(1:12,12) = (/
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443860D-01, +0.0000000000000000D+00
     &    /)
      H(1:12,1) = (/
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +9.5742710775633810D-01
     &    /)
      H(1:12,2) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +1.0000000000000000D+00, +9.5742710775633810D-01
     &    /)
      H(1:12,3) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +5.0000000000000011D-01, +8.1649658092772603D-01
     &    /)
      H(1:12,4) = (/
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01
     &    /)
      H(1:12,5) = (/
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +8.1649658092772603D-01, +9.5742710775633810D-01
     &    /)
      H(1:12,6) = (/
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +9.5742710775633810D-01, +8.1649658092772603D-01
     &    /)
      H(1:12,7) = (/
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +5.0000000000000011D-01
     &    /)
      H(1:12,8) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +8.1649658092772603D-01, +5.0000000000000011D-01
     &    /)
      H(1:12,9) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.5742710775633810D-01, +1.0000000000000000D+00
     &    /)
      H(1:12,10) = (/
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01
     &    /)
      H(1:12,11) = (/
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +1.0000000000000000D+00, +5.0000000000000011D-01
     &    /)
      H(1:12,12) = (/
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +5.0000000000000011D-01, +1.0000000000000000D+00
     &    /)
c
      return
      end subroutine
c
c     Interaction matrices for FCC system, unitary interactions
      subroutine mm10_mrr_GHfcc12(G,H)
c
      double precision, dimension(12,12)
     &    :: G, H
c
      G(1:12,1) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01
     &    /)
      G(1:12,2) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01
     &    /)
      G(1:12,3) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01
     &    /)
      G(1:12,4) = (/
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01
     &    /)
      G(1:12,5) = (/
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01
     &    /)
      G(1:12,6) = (/
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01
     &    /)
      G(1:12,7) = (/
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103173D-01, +4.7140452079103173D-01
     &    /)
      G(1:12,8) = (/
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.7140452079103168D-01, +4.7140452079103168D-01
     &    /)
      G(1:12,9) = (/
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.4280904158206347D-01, +9.4280904158206347D-01
     &    /)
      G(1:12,10) = (/
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00
     &    /)
      G(1:12,11) = (/
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00
     &    /)
      G(1:12,12) = (/
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103173D-01, +4.7140452079103173D-01,
     & +4.7140452079103168D-01, +4.7140452079103168D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +9.4280904158206347D-01, +9.4280904158206347D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00
     &    /)
      H(1:12,1) = (/
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01
     &    /)
      H(1:12,2) = (/
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01
     &    /)
      H(1:12,3) = (/
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01
     &    /)
      H(1:12,4) = (/
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01
     &    /)
      H(1:12,5) = (/
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01
     &    /)
      H(1:12,6) = (/
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01
     &    /)
      H(1:12,7) = (/
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819676D-01, +8.8191710368819676D-01
     &    /)
      H(1:12,8) = (/
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.8191710368819687D-01, +8.8191710368819687D-01
     &    /)
      H(1:12,9) = (/
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +3.3333333333333309D-01, +3.3333333333333309D-01
     &    /)
      H(1:12,10) = (/
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00
     &    /)
      H(1:12,11) = (/
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00
     &    /)
      H(1:12,12) = (/
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819676D-01, +8.8191710368819676D-01,
     & +8.8191710368819687D-01, +8.8191710368819687D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +3.3333333333333309D-01, +3.3333333333333309D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00
     &    /)
c
      return
      end subroutine
c
c     Interaction matrices for BCC system, unitary interactions
      subroutine mm10_mrr_GHbcc12(G,H)
c
      double precision, dimension(12,12)
     &    :: G, H
c
      G(1:12,1) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01
     &    /)
      G(1:12,2) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01
     &    /)
      G(1:12,3) = (/
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01
     &    /)
      G(1:12,4) = (/
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01
     &    /)
      G(1:12,5) = (/
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01
     &    /)
      G(1:12,6) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01
     &    /)
      G(1:12,7) = (/
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01
     &    /)
      G(1:12,8) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01
     &    /)
      G(1:12,9) = (/
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962573D-01, +5.7735026918962573D-01
     &    /)
      G(1:12,10) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962573D-01, +5.7735026918962573D-01
     &    /)
      G(1:12,11) = (/
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00
     &    /)
      G(1:12,12) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00
     &    /)
      H(1:12,1) = (/
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01
     &    /)
      H(1:12,2) = (/
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01
     &    /)
      H(1:12,3) = (/
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01
     &    /)
      H(1:12,4) = (/
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01
     &    /)
      H(1:12,5) = (/
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01
     &    /)
      H(1:12,6) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01
     &    /)
      H(1:12,7) = (/
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01
     &    /)
      H(1:12,8) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01
     &    /)
      H(1:12,9) = (/
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772603D-01, +8.1649658092772603D-01
     &    /)
      H(1:12,10) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772603D-01, +8.1649658092772603D-01
     &    /)
      H(1:12,11) = (/
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00
     &    /)
      H(1:12,12) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00
     &    /)
c
      return
      end subroutine
c
c     Interaction matrices for BCC system, unitary interactions
      subroutine mm10_mrr_GHbcc48(G,H)
c
      double precision, dimension(48,48)
     &    :: G, H
c
      G(1:48,1) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +6.6666666666666674D-01, +6.6666666666666674D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +8.3333333333333337D-01, +8.3333333333333337D-01,
     & +5.0000000000000000D-01, +5.0000000000000000D-01,
     & +8.3333333333333337D-01, +8.3333333333333337D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +2.1821789023599236D-01, +2.1821789023599236D-01,
     & +8.7287156094396956D-01, +8.7287156094396956D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +8.7287156094396956D-01, +8.7287156094396956D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +7.6376261582597338D-01, +1.0910894511799618D-01
     &    /)
      G(1:48,2) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +6.6666666666666674D-01, +6.6666666666666674D-01,
     & +3.3333333333333337D-01, +3.3333333333333337D-01,
     & +5.0000000000000000D-01, +5.0000000000000000D-01,
     & +1.6666666666666663D-01, +1.6666666666666663D-01,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +1.6666666666666663D-01, +1.6666666666666663D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +2.1821789023599236D-01, +2.1821789023599236D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +3.2732683535398865D-01, +3.2732683535398860D-01
     &    /)
      G(1:48,3) = (/
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +3.3333333333333337D-01, +3.3333333333333337D-01,
     & +6.6666666666666674D-01, +6.6666666666666674D-01,
     & +1.6666666666666663D-01, +1.6666666666666663D-01,
     & +5.0000000000000000D-01, +5.0000000000000000D-01,
     & +1.6666666666666663D-01, +1.6666666666666663D-01,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +2.1821789023599236D-01, +2.1821789023599236D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +1.0910894511799618D-01, +7.6376261582597338D-01
     &    /)
      G(1:48,4) = (/
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +6.6666666666666674D-01, +6.6666666666666674D-01,
     & +8.3333333333333337D-01, +8.3333333333333337D-01,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +8.3333333333333337D-01, +8.3333333333333337D-01,
     & +5.0000000000000000D-01, +5.0000000000000000D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +8.7287156094396956D-01, +8.7287156094396956D-01,
     & +2.1821789023599236D-01, +2.1821789023599236D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +8.7287156094396956D-01, +8.7287156094396956D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +3.2732683535398860D-01, +3.2732683535398865D-01
     &    /)
      G(1:48,5) = (/
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +5.0000000000000000D-01, +5.0000000000000000D-01,
     & +1.6666666666666663D-01, +1.6666666666666663D-01,
     & +6.6666666666666674D-01, +6.6666666666666674D-01,
     & +6.6666666666666674D-01, +6.6666666666666674D-01,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +5.0000000000000000D-01, +5.0000000000000000D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +2.1821789023599242D-01, +2.1821789023599242D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +2.1821789023599236D-01, +2.1821789023599236D-01,
     & +6.5465367070797720D-01, +1.3877787807814457D-17
     &    /)
      G(1:48,6) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +1.6666666666666663D-01, +1.6666666666666663D-01,
     & +5.0000000000000000D-01, +5.0000000000000000D-01,
     & +6.6666666666666674D-01, +6.6666666666666674D-01,
     & +6.6666666666666674D-01, +6.6666666666666674D-01,
     & +5.0000000000000000D-01, +5.0000000000000000D-01,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +2.1821789023599242D-01, +2.1821789023599242D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +2.1821789023599236D-01, +8.7287156094396956D-01
     &    /)
      G(1:48,7) = (/
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +8.3333333333333337D-01, +8.3333333333333337D-01,
     & +3.3333333333333331D-01, +3.3333333333333331D-01,
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +8.3333333333333348D-01, +8.3333333333333348D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.7287156094396967D-01, +8.7287156094396967D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +1.3877787807814457D-17, +1.3877787807814457D-17,
     & +8.7287156094396956D-01, +2.1821789023599236D-01
     &    /)
      G(1:48,8) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.3333333333333337D-01, +8.3333333333333337D-01,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +3.3333333333333331D-01, +3.3333333333333331D-01,
     & +8.3333333333333348D-01, +8.3333333333333348D-01,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +8.7287156094396967D-01, +8.7287156094396967D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +8.7287156094396956D-01, +8.7287156094396956D-01,
     & +1.3877787807814457D-17, +6.5465367070797720D-01
     &    /)
      G(1:48,9) = (/
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +5.0000000000000000D-01, +5.0000000000000000D-01,
     & +8.3333333333333337D-01, +8.3333333333333337D-01,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +8.3333333333333348D-01, +8.3333333333333348D-01,
     & +6.6666666666666674D-01, +6.6666666666666674D-01,
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +8.7287156094396967D-01, +8.7287156094396967D-01,
     & +2.1821789023599236D-01, +2.1821789023599236D-01,
     & +8.7287156094396956D-01, +8.7287156094396956D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +9.8198050606196574D-01, +3.2732683535398865D-01
     &    /)
      G(1:48,10) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +1.6666666666666663D-01, +1.6666666666666663D-01,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +5.0000000000000000D-01, +5.0000000000000000D-01,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +6.6666666666666674D-01, +6.6666666666666674D-01,
     & +3.3333333333333331D-01, +3.3333333333333331D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +2.1821789023599242D-01, +2.1821789023599242D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +1.3877787807814457D-17, +1.3877787807814457D-17,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +3.2732683535398865D-01, +9.8198050606196574D-01
     &    /)
      G(1:48,11) = (/
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +1.6666666666666663D-01, +1.6666666666666663D-01,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +5.0000000000000000D-01, +5.0000000000000000D-01,
     & +3.3333333333333331D-01, +3.3333333333333331D-01,
     & +6.6666666666666674D-01, +6.6666666666666674D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +2.1821789023599242D-01, +2.1821789023599242D-01,
     & +1.3877787807814457D-17, +1.3877787807814457D-17,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +7.6376261582597338D-01, +5.4554472558998102D-01
     &    /)
      G(1:48,12) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +2.8867513459481287D-01, +2.8867513459481287D-01,
     & +5.7735026918962573D-01, +5.7735026918962573D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.3333333333333337D-01, +8.3333333333333337D-01,
     & +5.0000000000000000D-01, +5.0000000000000000D-01,
     & +8.3333333333333348D-01, +8.3333333333333348D-01,
     & +1.6666666666666669D-01, +1.6666666666666669D-01,
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +6.6666666666666674D-01, +6.6666666666666674D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +3.2732683535398860D-01, +3.2732683535398860D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +8.7287156094396967D-01, +8.7287156094396967D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +8.7287156094396956D-01, +8.7287156094396956D-01,
     & +2.1821789023599236D-01, +2.1821789023599236D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +7.6376261582597338D-01, +7.6376261582597338D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +5.4554472558998102D-01, +7.6376261582597338D-01
     &    /)
      G(1:48,13) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962595D-01, +5.7735026918962595D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +9.4491118252306827D-01, +1.8898223650461365D-01
     &    /)
      G(1:48,14) = (/
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.7735026918962595D-01, +5.7735026918962595D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +1.8898223650461365D-01, +9.4491118252306827D-01
     &    /)
      G(1:48,15) = (/
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.7735026918962595D-01, +5.7735026918962595D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +1.8898223650461365D-01, +9.4491118252306827D-01
     &    /)
      G(1:48,16) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962595D-01, +5.7735026918962595D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +9.4491118252306827D-01, +1.8898223650461365D-01
     &    /)
      G(1:48,17) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +3.7796447300922731D-01, +7.5592894601845462D-01
     &    /)
      G(1:48,18) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +5.7735026918962595D-01, +5.7735026918962595D-01,
     & +5.7735026918962595D-01, +5.7735026918962595D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +7.5592894601845462D-01, +3.7796447300922731D-01
     &    /)
      G(1:48,19) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +3.7796447300922731D-01, +7.5592894601845462D-01
     &    /)
      G(1:48,20) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +5.7735026918962595D-01, +5.7735026918962595D-01,
     & +5.7735026918962595D-01, +5.7735026918962595D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +7.5592894601845462D-01, +3.7796447300922731D-01
     &    /)
      G(1:48,21) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962595D-01, +5.7735026918962595D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01
     &    /)
      G(1:48,22) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +5.7735026918962595D-01, +5.7735026918962595D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01
     &    /)
      G(1:48,23) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +5.7735026918962595D-01, +5.7735026918962595D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01
     &    /)
      G(1:48,24) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +1.0000000000000002D+00, +1.0000000000000002D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +2.8867513459481298D-01, +2.8867513459481298D-01,
     & +5.7735026918962595D-01, +5.7735026918962595D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +3.7796447300922731D-01, +3.7796447300922731D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01,
     & +5.6694670951384096D-01, +5.6694670951384096D-01,
     & +1.8898223650461365D-01, +1.8898223650461365D-01
     &    /)
      G(1:48,25) = (/
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +1.8898223650461368D-01, +1.8898223650461368D-01,
     & +6.9293486718358344D-01, +6.9293486718358344D-01,
     & +1.2598815766974242D-01, +1.2598815766974242D-01,
     & +8.8191710368819709D-01, +8.8191710368819709D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +9.4491118252306838D-01, +9.4491118252306838D-01,
     & +1.3877787807814457D-17, +1.3877787807814457D-17,
     & +6.5982887907385812D-01, +6.5982887907385812D-01,
     & +2.0619652471058070D-01, +2.0619652471058070D-01,
     & +7.8354679390020654D-01, +7.8354679390020654D-01,
     & +3.7115374447904520D-01, +3.7115374447904520D-01,
     & +7.0106818401597426D-01, +7.0106818401597426D-01,
     & +4.1239304942116146D-02, +4.1239304942116146D-02,
     & +9.4850401366867110D-01, +9.4850401366867110D-01,
     & +5.3611096424750981D-01, +5.3611096424750981D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +3.2991443953692917D-01, +3.2991443953692917D-01,
     & +9.8974331861078724D-01, +2.4743582965269689D-01
     &    /)
      G(1:48,26) = (/
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +4.4095855184409849D-01, +4.4095855184409849D-01,
     & +6.2994078834871209D-02, +6.2994078834871209D-02,
     & +7.5592894601845473D-01, +7.5592894601845473D-01,
     & +2.5197631533948489D-01, +2.5197631533948489D-01,
     & +8.1892302485332602D-01, +8.1892302485332602D-01,
     & +3.1497039417435613D-01, +3.1497039417435613D-01,
     & +4.1239304942116128D-01, +4.1239304942116128D-01,
     & +2.4743582965269678D-01, +2.4743582965269678D-01,
     & +6.1858957413174198D-01, +6.1858957413174198D-01,
     & +3.7115374447904526D-01, +3.7115374447904526D-01,
     & +4.5363235436327742D-01, +4.5363235436327742D-01,
     & +1.2371791482634840D-01, +1.2371791482634840D-01,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +1.2371791482634842D-01, +1.2371791482634842D-01,
     & +7.0106818401597437D-01, +7.0106818401597437D-01,
     & +3.7115374447904526D-01, +3.7115374447904526D-01,
     & +9.0726470872655496D-01, +9.0726470872655496D-01,
     & +2.4743582965269689D-01, +9.8974331861078724D-01
     &    /)
      G(1:48,27) = (/
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +6.9293486718358344D-01, +6.9293486718358344D-01,
     & +1.8898223650461368D-01, +1.8898223650461368D-01,
     & +8.8191710368819709D-01, +8.8191710368819709D-01,
     & +1.2598815766974242D-01, +1.2598815766974242D-01,
     & +9.4491118252306838D-01, +9.4491118252306838D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +6.5982887907385812D-01, +6.5982887907385812D-01,
     & +1.3877787807814457D-17, +1.3877787807814457D-17,
     & +7.8354679390020654D-01, +7.8354679390020654D-01,
     & +2.0619652471058070D-01, +2.0619652471058070D-01,
     & +7.0106818401597426D-01, +7.0106818401597426D-01,
     & +3.7115374447904520D-01, +3.7115374447904520D-01,
     & +9.4850401366867110D-01, +9.4850401366867110D-01,
     & +4.1239304942116146D-02, +4.1239304942116146D-02,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +5.3611096424750981D-01, +5.3611096424750981D-01,
     & +9.8974331861078724D-01, +9.8974331861078724D-01,
     & +3.2991443953692917D-01, +9.0726470872655496D-01
     &    /)
      G(1:48,28) = (/
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +6.2994078834871209D-02, +6.2994078834871209D-02,
     & +4.4095855184409849D-01, +4.4095855184409849D-01,
     & +2.5197631533948489D-01, +2.5197631533948489D-01,
     & +7.5592894601845473D-01, +7.5592894601845473D-01,
     & +3.1497039417435613D-01, +3.1497039417435613D-01,
     & +8.1892302485332602D-01, +8.1892302485332602D-01,
     & +2.4743582965269678D-01, +2.4743582965269678D-01,
     & +4.1239304942116128D-01, +4.1239304942116128D-01,
     & +3.7115374447904526D-01, +3.7115374447904526D-01,
     & +6.1858957413174198D-01, +6.1858957413174198D-01,
     & +1.2371791482634840D-01, +1.2371791482634840D-01,
     & +4.5363235436327742D-01, +4.5363235436327742D-01,
     & +1.2371791482634842D-01, +1.2371791482634842D-01,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +3.7115374447904526D-01, +3.7115374447904526D-01,
     & +7.0106818401597437D-01, +7.0106818401597437D-01,
     & +2.4743582965269689D-01, +2.4743582965269689D-01,
     & +9.0726470872655496D-01, +3.2991443953692917D-01
     &    /)
      G(1:48,29) = (/
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +1.2598815766974242D-01, +1.2598815766974242D-01,
     & +2.5197631533948489D-01, +2.5197631533948489D-01,
     & +1.8898223650461371D-01, +1.8898223650461371D-01,
     & +6.2994078834871209D-02, +6.2994078834871209D-02,
     & +4.4095855184409855D-01, +4.4095855184409855D-01,
     & +3.1497039417435618D-01, +3.1497039417435618D-01,
     & +2.0619652471058070D-01, +2.0619652471058070D-01,
     & +3.7115374447904520D-01, +3.7115374447904520D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +2.4743582965269678D-01, +2.4743582965269678D-01,
     & +4.1239304942116140D-02, +4.1239304942116140D-02,
     & +1.2371791482634842D-01, +1.2371791482634842D-01,
     & +3.7115374447904514D-01, +3.7115374447904514D-01,
     & +1.2371791482634839D-01, +1.2371791482634839D-01,
     & +3.2991443953692917D-01, +3.2991443953692917D-01,
     & +2.4743582965269689D-01, +2.4743582965269689D-01,
     & +5.3611096424750981D-01, +5.3611096424750981D-01,
     & +3.7115374447904526D-01, +8.6602540378443893D-01
     &    /)
      G(1:48,30) = (/
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +8.8191710368819720D-01, +8.8191710368819720D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +6.9293486718358355D-01, +6.9293486718358355D-01,
     & +8.1892302485332591D-01, +8.1892302485332591D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +6.1858957413174198D-01, +6.1858957413174198D-01,
     & +7.8354679390020654D-01, +7.8354679390020654D-01,
     & +4.1239304942116134D-01, +4.1239304942116134D-01,
     & +6.5982887907385812D-01, +6.5982887907385812D-01,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +9.4850401366867110D-01, +9.4850401366867110D-01,
     & +4.5363235436327742D-01, +4.5363235436327742D-01,
     & +7.0106818401597426D-01, +7.0106818401597426D-01,
     & +9.0726470872655507D-01, +9.0726470872655507D-01,
     & +9.8974331861078735D-01, +9.8974331861078735D-01,
     & +7.0106818401597437D-01, +7.0106818401597437D-01,
     & +8.6602540378443893D-01, +3.7115374447904526D-01
     &    /)
      G(1:48,31) = (/
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +2.5197631533948489D-01, +2.5197631533948489D-01,
     & +1.2598815766974242D-01, +1.2598815766974242D-01,
     & +6.2994078834871209D-02, +6.2994078834871209D-02,
     & +1.8898223650461371D-01, +1.8898223650461371D-01,
     & +3.1497039417435618D-01, +3.1497039417435618D-01,
     & +4.4095855184409855D-01, +4.4095855184409855D-01,
     & +3.7115374447904520D-01, +3.7115374447904520D-01,
     & +2.0619652471058070D-01, +2.0619652471058070D-01,
     & +2.4743582965269678D-01, +2.4743582965269678D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.2371791482634842D-01, +1.2371791482634842D-01,
     & +4.1239304942116140D-02, +4.1239304942116140D-02,
     & +1.2371791482634839D-01, +1.2371791482634839D-01,
     & +3.7115374447904514D-01, +3.7115374447904514D-01,
     & +2.4743582965269689D-01, +2.4743582965269689D-01,
     & +3.2991443953692917D-01, +3.2991443953692917D-01,
     & +3.7115374447904526D-01, +3.7115374447904526D-01,
     & +5.3611096424750981D-01, +7.0106818401597437D-01
     &    /)
      G(1:48,32) = (/
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +8.8191710368819720D-01, +8.8191710368819720D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +6.9293486718358355D-01, +6.9293486718358355D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +8.1892302485332591D-01, +8.1892302485332591D-01,
     & +7.8354679390020654D-01, +7.8354679390020654D-01,
     & +6.1858957413174198D-01, +6.1858957413174198D-01,
     & +6.5982887907385812D-01, +6.5982887907385812D-01,
     & +4.1239304942116134D-01, +4.1239304942116134D-01,
     & +9.4850401366867110D-01, +9.4850401366867110D-01,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +7.0106818401597426D-01, +7.0106818401597426D-01,
     & +4.5363235436327742D-01, +4.5363235436327742D-01,
     & +9.8974331861078735D-01, +9.8974331861078735D-01,
     & +9.0726470872655507D-01, +9.0726470872655507D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +7.0106818401597437D-01, +5.3611096424750981D-01
     &    /)
      G(1:48,33) = (/
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +1.8898223650461368D-01, +1.8898223650461368D-01,
     & +4.4095855184409849D-01, +4.4095855184409849D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +8.1892302485332602D-01, +8.1892302485332602D-01,
     & +1.2598815766974242D-01, +1.2598815766974242D-01,
     & +7.5592894601845473D-01, +7.5592894601845473D-01,
     & +3.7115374447904520D-01, +3.7115374447904520D-01,
     & +4.5363235436327742D-01, +4.5363235436327742D-01,
     & +5.3611096424750981D-01, +5.3611096424750981D-01,
     & +7.0106818401597437D-01, +7.0106818401597437D-01,
     & +1.3877787807814457D-17, +1.3877787807814457D-17,
     & +4.1239304942116128D-01, +4.1239304942116128D-01,
     & +3.2991443953692917D-01, +3.2991443953692917D-01,
     & +9.0726470872655496D-01, +9.0726470872655496D-01,
     & +2.0619652471058070D-01, +2.0619652471058070D-01,
     & +6.1858957413174198D-01, +6.1858957413174198D-01,
     & +4.1239304942116146D-02, +4.1239304942116146D-02,
     & +8.6602540378443882D-01, +1.2371791482634842D-01
     &    /)
      G(1:48,34) = (/
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +6.9293486718358344D-01, +6.9293486718358344D-01,
     & +6.2994078834871209D-02, +6.2994078834871209D-02,
     & +9.4491118252306838D-01, +9.4491118252306838D-01,
     & +3.1497039417435613D-01, +3.1497039417435613D-01,
     & +8.8191710368819709D-01, +8.8191710368819709D-01,
     & +2.5197631533948489D-01, +2.5197631533948489D-01,
     & +7.0106818401597426D-01, +7.0106818401597426D-01,
     & +1.2371791482634840D-01, +1.2371791482634840D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +3.7115374447904526D-01, +3.7115374447904526D-01,
     & +6.5982887907385812D-01, +6.5982887907385812D-01,
     & +2.4743582965269678D-01, +2.4743582965269678D-01,
     & +9.8974331861078724D-01, +9.8974331861078724D-01,
     & +2.4743582965269689D-01, +2.4743582965269689D-01,
     & +7.8354679390020654D-01, +7.8354679390020654D-01,
     & +3.7115374447904526D-01, +3.7115374447904526D-01,
     & +9.4850401366867110D-01, +9.4850401366867110D-01,
     & +1.2371791482634842D-01, +8.6602540378443882D-01
     &    /)
      G(1:48,35) = (/
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +4.4095855184409849D-01, +4.4095855184409849D-01,
     & +1.8898223650461368D-01, +1.8898223650461368D-01,
     & +8.1892302485332602D-01, +8.1892302485332602D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +7.5592894601845473D-01, +7.5592894601845473D-01,
     & +1.2598815766974242D-01, +1.2598815766974242D-01,
     & +4.5363235436327742D-01, +4.5363235436327742D-01,
     & +3.7115374447904520D-01, +3.7115374447904520D-01,
     & +7.0106818401597437D-01, +7.0106818401597437D-01,
     & +5.3611096424750981D-01, +5.3611096424750981D-01,
     & +4.1239304942116128D-01, +4.1239304942116128D-01,
     & +1.3877787807814457D-17, +1.3877787807814457D-17,
     & +9.0726470872655496D-01, +9.0726470872655496D-01,
     & +3.2991443953692917D-01, +3.2991443953692917D-01,
     & +6.1858957413174198D-01, +6.1858957413174198D-01,
     & +2.0619652471058070D-01, +2.0619652471058070D-01,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +4.1239304942116146D-02, +9.4850401366867110D-01
     &    /)
      G(1:48,36) = (/
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +6.2994078834871209D-02, +6.2994078834871209D-02,
     & +6.9293486718358344D-01, +6.9293486718358344D-01,
     & +3.1497039417435613D-01, +3.1497039417435613D-01,
     & +9.4491118252306838D-01, +9.4491118252306838D-01,
     & +2.5197631533948489D-01, +2.5197631533948489D-01,
     & +8.8191710368819709D-01, +8.8191710368819709D-01,
     & +1.2371791482634840D-01, +1.2371791482634840D-01,
     & +7.0106818401597426D-01, +7.0106818401597426D-01,
     & +3.7115374447904526D-01, +3.7115374447904526D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +2.4743582965269678D-01, +2.4743582965269678D-01,
     & +6.5982887907385812D-01, +6.5982887907385812D-01,
     & +2.4743582965269689D-01, +2.4743582965269689D-01,
     & +9.8974331861078724D-01, +9.8974331861078724D-01,
     & +3.7115374447904526D-01, +3.7115374447904526D-01,
     & +7.8354679390020654D-01, +7.8354679390020654D-01,
     & +1.2371791482634842D-01, +1.2371791482634842D-01,
     & +9.4850401366867110D-01, +4.1239304942116146D-02
     &    /)
      G(1:48,37) = (/
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +3.1497039417435613D-01, +3.1497039417435613D-01,
     & +1.8898223650461368D-01, +1.8898223650461368D-01,
     & +6.2994078834871181D-02, +6.2994078834871181D-02,
     & +1.2598815766974242D-01, +1.2598815766974242D-01,
     & +2.5197631533948484D-01, +2.5197631533948484D-01,
     & +5.3611096424750981D-01, +5.3611096424750981D-01,
     & +3.7115374447904531D-01, +3.7115374447904531D-01,
     & +3.7115374447904520D-01, +3.7115374447904520D-01,
     & +1.2371791482634842D-01, +1.2371791482634842D-01,
     & +3.2991443953692917D-01, +3.2991443953692917D-01,
     & +2.4743582965269689D-01, +2.4743582965269689D-01,
     & +2.7755575615628914D-17, +2.7755575615628914D-17,
     & +2.4743582965269681D-01, +2.4743582965269681D-01,
     & +4.1239304942116140D-02, +4.1239304942116140D-02,
     & +1.2371791482634836D-01, +1.2371791482634836D-01,
     & +2.0619652471058067D-01, +2.0619652471058067D-01,
     & +3.7115374447904514D-01, +6.1858957413174198D-01
     &    /)
      G(1:48,38) = (/
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +9.4491118252306838D-01, +9.4491118252306838D-01,
     & +8.1892302485332591D-01, +8.1892302485332591D-01,
     & +6.9293486718358344D-01, +6.9293486718358344D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +8.8191710368819709D-01, +8.8191710368819709D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +7.0106818401597437D-01, +7.0106818401597437D-01,
     & +7.0106818401597426D-01, +7.0106818401597426D-01,
     & +4.5363235436327748D-01, +4.5363235436327748D-01,
     & +9.8974331861078735D-01, +9.8974331861078735D-01,
     & +9.0726470872655507D-01, +9.0726470872655507D-01,
     & +6.5982887907385812D-01, +6.5982887907385812D-01,
     & +4.1239304942116128D-01, +4.1239304942116128D-01,
     & +9.4850401366867110D-01, +9.4850401366867110D-01,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +7.8354679390020654D-01, +7.8354679390020654D-01,
     & +6.1858957413174198D-01, +3.7115374447904514D-01
     &    /)
      G(1:48,39) = (/
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +3.1497039417435613D-01, +3.1497039417435613D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +6.2994078834871181D-02, +6.2994078834871181D-02,
     & +1.8898223650461368D-01, +1.8898223650461368D-01,
     & +2.5197631533948484D-01, +2.5197631533948484D-01,
     & +1.2598815766974242D-01, +1.2598815766974242D-01,
     & +3.7115374447904531D-01, +3.7115374447904531D-01,
     & +5.3611096424750981D-01, +5.3611096424750981D-01,
     & +1.2371791482634842D-01, +1.2371791482634842D-01,
     & +3.7115374447904520D-01, +3.7115374447904520D-01,
     & +2.4743582965269689D-01, +2.4743582965269689D-01,
     & +3.2991443953692917D-01, +3.2991443953692917D-01,
     & +2.4743582965269681D-01, +2.4743582965269681D-01,
     & +2.7755575615628914D-17, +2.7755575615628914D-17,
     & +1.2371791482634836D-01, +1.2371791482634836D-01,
     & +4.1239304942116140D-02, +4.1239304942116140D-02,
     & +3.7115374447904514D-01, +3.7115374447904514D-01,
     & +2.0619652471058067D-01, +7.8354679390020654D-01
     &    /)
      G(1:48,40) = (/
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +8.1892302485332591D-01, +8.1892302485332591D-01,
     & +9.4491118252306838D-01, +9.4491118252306838D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +6.9293486718358344D-01, +6.9293486718358344D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +8.8191710368819709D-01, +8.8191710368819709D-01,
     & +7.0106818401597437D-01, +7.0106818401597437D-01,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +4.5363235436327748D-01, +4.5363235436327748D-01,
     & +7.0106818401597426D-01, +7.0106818401597426D-01,
     & +9.0726470872655507D-01, +9.0726470872655507D-01,
     & +9.8974331861078735D-01, +9.8974331861078735D-01,
     & +4.1239304942116128D-01, +4.1239304942116128D-01,
     & +6.5982887907385812D-01, +6.5982887907385812D-01,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +9.4850401366867110D-01, +9.4850401366867110D-01,
     & +6.1858957413174198D-01, +6.1858957413174198D-01,
     & +7.8354679390020654D-01, +2.0619652471058067D-01
     &    /)
      G(1:48,41) = (/
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +1.2598815766974242D-01, +1.2598815766974242D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +4.4095855184409855D-01, +4.4095855184409855D-01,
     & +8.1892302485332591D-01, +8.1892302485332591D-01,
     & +1.8898223650461371D-01, +1.8898223650461371D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +4.1239304942116140D-02, +4.1239304942116140D-02,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +3.2991443953692917D-01, +3.2991443953692917D-01,
     & +9.0726470872655507D-01, +9.0726470872655507D-01,
     & +2.0619652471058070D-01, +2.0619652471058070D-01,
     & +6.1858957413174198D-01, +6.1858957413174198D-01,
     & +5.3611096424750981D-01, +5.3611096424750981D-01,
     & +7.0106818401597437D-01, +7.0106818401597437D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.1239304942116134D-01, +4.1239304942116134D-01,
     & +3.7115374447904514D-01, +3.7115374447904514D-01,
     & +4.5363235436327742D-01, +7.0106818401597426D-01
     &    /)
      G(1:48,42) = (/
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +2.5197631533948489D-01, +2.5197631533948489D-01,
     & +8.8191710368819720D-01, +8.8191710368819720D-01,
     & +3.1497039417435618D-01, +3.1497039417435618D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +6.2994078834871209D-02, +6.2994078834871209D-02,
     & +6.9293486718358355D-01, +6.9293486718358355D-01,
     & +1.2371791482634842D-01, +1.2371791482634842D-01,
     & +9.4850401366867110D-01, +9.4850401366867110D-01,
     & +2.4743582965269689D-01, +2.4743582965269689D-01,
     & +9.8974331861078735D-01, +9.8974331861078735D-01,
     & +3.7115374447904520D-01, +3.7115374447904520D-01,
     & +7.8354679390020654D-01, +7.8354679390020654D-01,
     & +3.7115374447904526D-01, +3.7115374447904526D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +2.4743582965269678D-01, +2.4743582965269678D-01,
     & +6.5982887907385812D-01, +6.5982887907385812D-01,
     & +1.2371791482634839D-01, +1.2371791482634839D-01,
     & +7.0106818401597426D-01, +4.5363235436327742D-01
     &    /)
      G(1:48,43) = (/
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +1.2598815766974242D-01, +1.2598815766974242D-01,
     & +8.1892302485332591D-01, +8.1892302485332591D-01,
     & +4.4095855184409855D-01, +4.4095855184409855D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +1.8898223650461371D-01, +1.8898223650461371D-01,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +4.1239304942116140D-02, +4.1239304942116140D-02,
     & +9.0726470872655507D-01, +9.0726470872655507D-01,
     & +3.2991443953692917D-01, +3.2991443953692917D-01,
     & +6.1858957413174198D-01, +6.1858957413174198D-01,
     & +2.0619652471058070D-01, +2.0619652471058070D-01,
     & +7.0106818401597437D-01, +7.0106818401597437D-01,
     & +5.3611096424750981D-01, +5.3611096424750981D-01,
     & +4.1239304942116134D-01, +4.1239304942116134D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +4.5363235436327742D-01, +4.5363235436327742D-01,
     & +3.7115374447904514D-01, +1.2371791482634839D-01
     &    /)
      G(1:48,44) = (/
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +8.8191710368819720D-01, +8.8191710368819720D-01,
     & +2.5197631533948489D-01, +2.5197631533948489D-01,
     & +9.4491118252306827D-01, +9.4491118252306827D-01,
     & +3.1497039417435618D-01, +3.1497039417435618D-01,
     & +6.9293486718358355D-01, +6.9293486718358355D-01,
     & +6.2994078834871209D-02, +6.2994078834871209D-02,
     & +9.4850401366867110D-01, +9.4850401366867110D-01,
     & +1.2371791482634842D-01, +1.2371791482634842D-01,
     & +9.8974331861078735D-01, +9.8974331861078735D-01,
     & +2.4743582965269689D-01, +2.4743582965269689D-01,
     & +7.8354679390020654D-01, +7.8354679390020654D-01,
     & +3.7115374447904520D-01, +3.7115374447904520D-01,
     & +8.6602540378443893D-01, +8.6602540378443893D-01,
     & +3.7115374447904526D-01, +3.7115374447904526D-01,
     & +6.5982887907385812D-01, +6.5982887907385812D-01,
     & +2.4743582965269678D-01, +2.4743582965269678D-01,
     & +7.0106818401597426D-01, +7.0106818401597426D-01,
     & +1.2371791482634839D-01, +3.7115374447904514D-01
     &    /)
      G(1:48,45) = (/
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +9.4491118252306838D-01, +9.4491118252306838D-01,
     & +1.2598815766974242D-01, +1.2598815766974242D-01,
     & +8.8191710368819709D-01, +8.8191710368819709D-01,
     & +1.8898223650461368D-01, +1.8898223650461368D-01,
     & +6.9293486718358344D-01, +6.9293486718358344D-01,
     & +3.2991443953692917D-01, +3.2991443953692917D-01,
     & +9.8974331861078735D-01, +9.8974331861078735D-01,
     & +4.1239304942116140D-02, +4.1239304942116140D-02,
     & +9.4850401366867110D-01, +9.4850401366867110D-01,
     & +5.3611096424750981D-01, +5.3611096424750981D-01,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +2.0619652471058067D-01, +2.0619652471058067D-01,
     & +7.8354679390020654D-01, +7.8354679390020654D-01,
     & +3.7115374447904520D-01, +3.7115374447904520D-01,
     & +7.0106818401597426D-01, +7.0106818401597426D-01,
     & +2.7755575615628914D-17, +2.7755575615628914D-17,
     & +6.5982887907385812D-01, +4.1239304942116128D-01
     &    /)
      G(1:48,46) = (/
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +3.1497039417435613D-01, +3.1497039417435613D-01,
     & +8.1892302485332591D-01, +8.1892302485332591D-01,
     & +2.5197631533948484D-01, +2.5197631533948484D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +6.2994078834871181D-02, +6.2994078834871181D-02,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +2.4743582965269689D-01, +2.4743582965269689D-01,
     & +9.0726470872655507D-01, +9.0726470872655507D-01,
     & +1.2371791482634836D-01, +1.2371791482634836D-01,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +3.7115374447904531D-01, +3.7115374447904531D-01,
     & +7.0106818401597437D-01, +7.0106818401597437D-01,
     & +3.7115374447904514D-01, +3.7115374447904514D-01,
     & +6.1858957413174198D-01, +6.1858957413174198D-01,
     & +1.2371791482634842D-01, +1.2371791482634842D-01,
     & +4.5363235436327748D-01, +4.5363235436327748D-01,
     & +2.4743582965269681D-01, +2.4743582965269681D-01,
     & +4.1239304942116128D-01, +6.5982887907385812D-01
     &    /)
      G(1:48,47) = (/
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +9.4491118252306838D-01, +9.4491118252306838D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +8.8191710368819709D-01, +8.8191710368819709D-01,
     & +1.2598815766974242D-01, +1.2598815766974242D-01,
     & +6.9293486718358344D-01, +6.9293486718358344D-01,
     & +1.8898223650461368D-01, +1.8898223650461368D-01,
     & +9.8974331861078735D-01, +9.8974331861078735D-01,
     & +3.2991443953692917D-01, +3.2991443953692917D-01,
     & +9.4850401366867110D-01, +9.4850401366867110D-01,
     & +4.1239304942116140D-02, +4.1239304942116140D-02,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +5.3611096424750981D-01, +5.3611096424750981D-01,
     & +7.8354679390020654D-01, +7.8354679390020654D-01,
     & +2.0619652471058067D-01, +2.0619652471058067D-01,
     & +7.0106818401597426D-01, +7.0106818401597426D-01,
     & +3.7115374447904520D-01, +3.7115374447904520D-01,
     & +6.5982887907385812D-01, +6.5982887907385812D-01,
     & +2.7755575615628914D-17, +2.4743582965269681D-01
     &    /)
      G(1:48,48) = (/
     & +3.2732683535398865D-01, +3.2732683535398865D-01,
     & +5.4554472558998102D-01, +5.4554472558998102D-01,
     & +4.3643578047198484D-01, +4.3643578047198484D-01,
     & +6.5465367070797720D-01, +6.5465367070797720D-01,
     & +9.8198050606196585D-01, +9.8198050606196585D-01,
     & +1.0910894511799618D-01, +1.0910894511799618D-01,
     & +8.1892302485332591D-01, +8.1892302485332591D-01,
     & +3.1497039417435613D-01, +3.1497039417435613D-01,
     & +7.5592894601845462D-01, +7.5592894601845462D-01,
     & +2.5197631533948484D-01, +2.5197631533948484D-01,
     & +4.4095855184409860D-01, +4.4095855184409860D-01,
     & +6.2994078834871181D-02, +6.2994078834871181D-02,
     & +9.0726470872655507D-01, +9.0726470872655507D-01,
     & +2.4743582965269689D-01, +2.4743582965269689D-01,
     & +8.6602540378443882D-01, +8.6602540378443882D-01,
     & +1.2371791482634836D-01, +1.2371791482634836D-01,
     & +7.0106818401597437D-01, +7.0106818401597437D-01,
     & +3.7115374447904531D-01, +3.7115374447904531D-01,
     & +6.1858957413174198D-01, +6.1858957413174198D-01,
     & +3.7115374447904514D-01, +3.7115374447904514D-01,
     & +4.5363235436327748D-01, +4.5363235436327748D-01,
     & +1.2371791482634842D-01, +1.2371791482634842D-01,
     & +4.1239304942116128D-01, +4.1239304942116128D-01,
     & +2.4743582965269681D-01, +2.7755575615628914D-17
     &    /)
      H(1:48,1) = (/
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +7.4535599249992979D-01, +7.4535599249992979D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +5.5277079839256660D-01, +5.5277079839256660D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +5.5277079839256660D-01, +5.5277079839256660D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +1.8898223650461357D-01, +1.8898223650461357D-01,
     & +9.7590007294853320D-01, +9.7590007294853320D-01,
     & +4.8795003647426655D-01, +4.8795003647426655D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +1.8898223650461357D-01, +1.8898223650461357D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +4.8795003647426655D-01, +4.8795003647426655D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +6.4549722436790280D-01, +9.9402979738800490D-01
     &    /)
      H(1:48,2) = (/
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +7.4535599249992979D-01, +7.4535599249992979D-01,
     & +9.4280904158206336D-01, +9.4280904158206336D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.7590007294853320D-01, +9.7590007294853320D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01
     &    /)
      H(1:48,3) = (/
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.4280904158206336D-01, +9.4280904158206336D-01,
     & +7.4535599249992979D-01, +7.4535599249992979D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.7590007294853320D-01, +9.7590007294853320D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.9402979738800490D-01, +6.4549722436790280D-01
     &    /)
      H(1:48,4) = (/
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +7.4535599249992979D-01, +7.4535599249992979D-01,
     & +5.5277079839256660D-01, +5.5277079839256660D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +5.5277079839256660D-01, +5.5277079839256660D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +1.8898223650461357D-01, +1.8898223650461357D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +4.8795003647426655D-01, +4.8795003647426655D-01,
     & +9.7590007294853320D-01, +9.7590007294853320D-01,
     & +1.8898223650461357D-01, +1.8898223650461357D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +4.8795003647426655D-01, +4.8795003647426655D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01
     &    /)
      H(1:48,5) = (/
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +7.4535599249992979D-01, +7.4535599249992979D-01,
     & +7.4535599249992979D-01, +7.4535599249992979D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +9.7590007294853320D-01, +9.7590007294853320D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.7590007294853320D-01, +9.7590007294853320D-01,
     & +7.5592894601845440D-01, +1.0000000000000000D+00
     &    /)
      H(1:48,6) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +7.4535599249992979D-01, +7.4535599249992979D-01,
     & +7.4535599249992979D-01, +7.4535599249992979D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +9.7590007294853320D-01, +9.7590007294853320D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +9.7590007294853320D-01, +4.8795003647426655D-01
     &    /)
      H(1:48,7) = (/
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +5.5277079839256660D-01, +5.5277079839256660D-01,
     & +9.4280904158206336D-01, +9.4280904158206336D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +5.5277079839256649D-01, +5.5277079839256649D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +4.8795003647426627D-01, +4.8795003647426627D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +1.8898223650461357D-01, +1.8898223650461357D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +1.8898223650461357D-01, +1.8898223650461357D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +4.8795003647426655D-01, +9.7590007294853320D-01
     &    /)
      H(1:48,8) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.5277079839256660D-01, +5.5277079839256660D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.4280904158206336D-01, +9.4280904158206336D-01,
     & +5.5277079839256649D-01, +5.5277079839256649D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +4.8795003647426627D-01, +4.8795003647426627D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.8898223650461357D-01, +1.8898223650461357D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +1.8898223650461357D-01, +1.8898223650461357D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +4.8795003647426655D-01, +4.8795003647426655D-01,
     & +1.0000000000000000D+00, +7.5592894601845440D-01
     &    /)
      H(1:48,9) = (/
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +5.5277079839256660D-01, +5.5277079839256660D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +5.5277079839256649D-01, +5.5277079839256649D-01,
     & +7.4535599249992979D-01, +7.4535599249992979D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +4.8795003647426627D-01, +4.8795003647426627D-01,
     & +9.7590007294853320D-01, +9.7590007294853320D-01,
     & +4.8795003647426655D-01, +4.8795003647426655D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +1.8898223650461357D-01, +1.8898223650461357D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +1.8898223650461357D-01, +9.4491118252306805D-01
     &    /)
      H(1:48,10) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +7.4535599249992979D-01, +7.4535599249992979D-01,
     & +9.4280904158206336D-01, +9.4280904158206336D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.7590007294853320D-01, +9.7590007294853320D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +9.4491118252306805D-01, +1.8898223650461357D-01
     &    /)
      H(1:48,11) = (/
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +9.4280904158206336D-01, +9.4280904158206336D-01,
     & +7.4535599249992979D-01, +7.4535599249992979D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.7590007294853320D-01, +9.7590007294853320D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +6.4549722436790280D-01, +8.3808170984752572D-01
     &    /)
      H(1:48,12) = (/
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +5.0000000000000011D-01, +5.0000000000000011D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +8.1649658092772603D-01, +8.1649658092772603D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +5.5277079839256660D-01, +5.5277079839256660D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +5.5277079839256649D-01, +5.5277079839256649D-01,
     & +9.8601329718326935D-01, +9.8601329718326935D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +7.4535599249992979D-01, +7.4535599249992979D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +4.8795003647426627D-01, +4.8795003647426627D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +4.8795003647426655D-01, +4.8795003647426655D-01,
     & +9.7590007294853320D-01, +9.7590007294853320D-01,
     & +1.8898223650461357D-01, +1.8898223650461357D-01,
     & +6.4549722436790280D-01, +6.4549722436790280D-01,
     & +1.8898223650461357D-01, +1.8898223650461357D-01,
     & +8.3808170984752572D-01, +6.4549722436790280D-01
     &    /)
      H(1:48,13) = (/
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772592D-01, +8.1649658092772592D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398799D-01, +9.8198050606196574D-01
     &    /)
      H(1:48,14) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.1649658092772592D-01, +8.1649658092772592D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.8198050606196574D-01, +3.2732683535398799D-01
     &    /)
      H(1:48,15) = (/
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.1649658092772592D-01, +8.1649658092772592D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.8198050606196574D-01, +3.2732683535398799D-01
     &    /)
      H(1:48,16) = (/
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772592D-01, +8.1649658092772592D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398799D-01, +9.8198050606196574D-01
     &    /)
      H(1:48,17) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +9.2582009977255142D-01, +6.5465367070797686D-01
     &    /)
      H(1:48,18) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +8.1649658092772592D-01, +8.1649658092772592D-01,
     & +8.1649658092772592D-01, +8.1649658092772592D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +6.5465367070797686D-01, +9.2582009977255142D-01
     &    /)
      H(1:48,19) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +9.2582009977255142D-01, +6.5465367070797686D-01
     &    /)
      H(1:48,20) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +8.1649658092772592D-01, +8.1649658092772592D-01,
     & +8.1649658092772592D-01, +8.1649658092772592D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +6.5465367070797686D-01, +9.2582009977255142D-01
     &    /)
      H(1:48,21) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772592D-01, +8.1649658092772592D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01
     &    /)
      H(1:48,22) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.1649658092772592D-01, +8.1649658092772592D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01
     &    /)
      H(1:48,23) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +8.1649658092772592D-01, +8.1649658092772592D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01
     &    /)
      H(1:48,24) = (/
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +8.6602540378443860D-01, +8.6602540378443860D-01,
     & +0.0000000000000000D+00, +0.0000000000000000D+00,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.5742710775633810D-01, +9.5742710775633810D-01,
     & +8.1649658092772592D-01, +8.1649658092772592D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.2582009977255142D-01, +9.2582009977255142D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +8.2375447104791388D-01, +8.2375447104791388D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01
     &    /)
      H(1:48,25) = (/
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +7.2100018712984359D-01, +7.2100018712984359D-01,
     & +9.9203174552379325D-01, +9.9203174552379325D-01,
     & +4.7140452079103129D-01, +4.7140452079103129D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +3.2732683535398766D-01, +3.2732683535398766D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +7.5141589705045231D-01, +7.5141589705045231D-01,
     & +9.7851059943021512D-01, +9.7851059943021512D-01,
     & +6.2133277860475644D-01, +6.2133277860475644D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +7.1309424437485402D-01, +7.1309424437485402D-01,
     & +9.9914929801701369D-01, +9.9914929801701369D-01,
     & +3.1676511180119155D-01, +3.1676511180119155D-01,
     & +8.4414751910646835D-01, +8.4414751910646835D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.4401083817138132D-01, +9.4401083817138132D-01,
     & +1.4285714285714138D-01, +9.6890428330360967D-01
     &    /)
      H(1:48,26) = (/
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +8.9752746785575066D-01, +8.9752746785575066D-01,
     & +9.9801390072069940D-01, +9.9801390072069940D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.6773340156674170D-01, +9.6773340156674170D-01,
     & +5.7390337110447320D-01, +5.7390337110447320D-01,
     & +9.4910149657117848D-01, +9.4910149657117848D-01,
     & +9.1100602236709483D-01, +9.1100602236709483D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +7.8571428571428570D-01, +7.8571428571428570D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +8.9118891772442388D-01, +8.9118891772442388D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +7.1309424437485391D-01, +7.1309424437485391D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +4.2056004125370655D-01, +4.2056004125370655D-01,
     & +9.6890428330360967D-01, +1.4285714285714138D-01
     &    /)
      H(1:48,27) = (/
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +7.2100018712984359D-01, +7.2100018712984359D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +4.7140452079103129D-01, +4.7140452079103129D-01,
     & +9.9203174552379325D-01, +9.9203174552379325D-01,
     & +3.2732683535398766D-01, +3.2732683535398766D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +7.5141589705045231D-01, +7.5141589705045231D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +6.2133277860475644D-01, +6.2133277860475644D-01,
     & +9.7851059943021512D-01, +9.7851059943021512D-01,
     & +7.1309424437485402D-01, +7.1309424437485402D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +3.1676511180119155D-01, +3.1676511180119155D-01,
     & +9.9914929801701369D-01, +9.9914929801701369D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +8.4414751910646835D-01, +8.4414751910646835D-01,
     & +1.4285714285714138D-01, +1.4285714285714138D-01,
     & +9.4401083817138132D-01, +4.2056004125370655D-01
     &    /)
      H(1:48,28) = (/
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.9801390072069940D-01, +9.9801390072069940D-01,
     & +8.9752746785575066D-01, +8.9752746785575066D-01,
     & +9.6773340156674170D-01, +9.6773340156674170D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.4910149657117848D-01, +9.4910149657117848D-01,
     & +5.7390337110447320D-01, +5.7390337110447320D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +9.1100602236709483D-01, +9.1100602236709483D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +7.8571428571428570D-01, +7.8571428571428570D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +8.9118891772442388D-01, +8.9118891772442388D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +7.1309424437485391D-01, +7.1309424437485391D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +4.2056004125370655D-01, +9.4401083817138132D-01
     &    /)
      H(1:48,29) = (/
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.9203174552379325D-01, +9.9203174552379325D-01,
     & +9.6773340156674170D-01, +9.6773340156674170D-01,
     & +9.8198050606196563D-01, +9.8198050606196563D-01,
     & +9.9801390072069940D-01, +9.9801390072069940D-01,
     & +8.9752746785575066D-01, +8.9752746785575066D-01,
     & +9.4910149657117848D-01, +9.4910149657117848D-01,
     & +9.7851059943021512D-01, +9.7851059943021512D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +9.9914929801701369D-01, +9.9914929801701369D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +9.4401083817138132D-01, +9.4401083817138132D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +8.4414751910646835D-01, +8.4414751910646835D-01,
     & +9.2857142857142849D-01, +4.9999999999999956D-01
     &    /)
      H(1:48,30) = (/
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +4.7140452079103107D-01, +4.7140452079103107D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +7.2100018712984359D-01, +7.2100018712984359D-01,
     & +5.7390337110447343D-01, +5.7390337110447343D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +7.8571428571428570D-01, +7.8571428571428570D-01,
     & +6.2133277860475644D-01, +6.2133277860475644D-01,
     & +9.1100602236709471D-01, +9.1100602236709471D-01,
     & +7.5141589705045231D-01, +7.5141589705045231D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +3.1676511180119155D-01, +3.1676511180119155D-01,
     & +8.9118891772442388D-01, +8.9118891772442388D-01,
     & +7.1309424437485402D-01, +7.1309424437485402D-01,
     & +4.2056004125370633D-01, +4.2056004125370633D-01,
     & +1.4285714285714060D-01, +1.4285714285714060D-01,
     & +7.1309424437485391D-01, +7.1309424437485391D-01,
     & +4.9999999999999956D-01, +9.2857142857142849D-01
     &    /)
      H(1:48,31) = (/
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.6773340156674170D-01, +9.6773340156674170D-01,
     & +9.9203174552379325D-01, +9.9203174552379325D-01,
     & +9.9801390072069940D-01, +9.9801390072069940D-01,
     & +9.8198050606196563D-01, +9.8198050606196563D-01,
     & +9.4910149657117848D-01, +9.4910149657117848D-01,
     & +8.9752746785575066D-01, +8.9752746785575066D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +9.7851059943021512D-01, +9.7851059943021512D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +9.9914929801701369D-01, +9.9914929801701369D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +9.4401083817138132D-01, +9.4401083817138132D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +8.4414751910646835D-01, +7.1309424437485391D-01
     &    /)
      H(1:48,32) = (/
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +4.7140452079103107D-01, +4.7140452079103107D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +7.2100018712984359D-01, +7.2100018712984359D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +5.7390337110447343D-01, +5.7390337110447343D-01,
     & +6.2133277860475644D-01, +6.2133277860475644D-01,
     & +7.8571428571428570D-01, +7.8571428571428570D-01,
     & +7.5141589705045231D-01, +7.5141589705045231D-01,
     & +9.1100602236709471D-01, +9.1100602236709471D-01,
     & +3.1676511180119155D-01, +3.1676511180119155D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +7.1309424437485402D-01, +7.1309424437485402D-01,
     & +8.9118891772442388D-01, +8.9118891772442388D-01,
     & +1.4285714285714060D-01, +1.4285714285714060D-01,
     & +4.2056004125370633D-01, +4.2056004125370633D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +7.1309424437485391D-01, +8.4414751910646835D-01
     &    /)
      H(1:48,33) = (/
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +8.9752746785575066D-01, +8.9752746785575066D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +5.7390337110447320D-01, +5.7390337110447320D-01,
     & +9.9203174552379325D-01, +9.9203174552379325D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +8.9118891772442388D-01, +8.9118891772442388D-01,
     & +8.4414751910646835D-01, +8.4414751910646835D-01,
     & +7.1309424437485391D-01, +7.1309424437485391D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.1100602236709483D-01, +9.1100602236709483D-01,
     & +9.4401083817138132D-01, +9.4401083817138132D-01,
     & +4.2056004125370655D-01, +4.2056004125370655D-01,
     & +9.7851059943021512D-01, +9.7851059943021512D-01,
     & +7.8571428571428570D-01, +7.8571428571428570D-01,
     & +9.9914929801701369D-01, +9.9914929801701369D-01,
     & +4.9999999999999967D-01, +9.9231742781784316D-01
     &    /)
      H(1:48,34) = (/
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +7.2100018712984359D-01, +7.2100018712984359D-01,
     & +9.9801390072069940D-01, +9.9801390072069940D-01,
     & +3.2732683535398766D-01, +3.2732683535398766D-01,
     & +9.4910149657117848D-01, +9.4910149657117848D-01,
     & +4.7140452079103129D-01, +4.7140452079103129D-01,
     & +9.6773340156674170D-01, +9.6773340156674170D-01,
     & +7.1309424437485402D-01, +7.1309424437485402D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +7.5141589705045231D-01, +7.5141589705045231D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +1.4285714285714138D-01, +1.4285714285714138D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +6.2133277860475644D-01, +6.2133277860475644D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +3.1676511180119155D-01, +3.1676511180119155D-01,
     & +9.9231742781784316D-01, +4.9999999999999967D-01
     &    /)
      H(1:48,35) = (/
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.9752746785575066D-01, +8.9752746785575066D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +5.7390337110447320D-01, +5.7390337110447320D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.9203174552379325D-01, +9.9203174552379325D-01,
     & +8.9118891772442388D-01, +8.9118891772442388D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +7.1309424437485391D-01, +7.1309424437485391D-01,
     & +8.4414751910646835D-01, +8.4414751910646835D-01,
     & +9.1100602236709483D-01, +9.1100602236709483D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +4.2056004125370655D-01, +4.2056004125370655D-01,
     & +9.4401083817138132D-01, +9.4401083817138132D-01,
     & +7.8571428571428570D-01, +7.8571428571428570D-01,
     & +9.7851059943021512D-01, +9.7851059943021512D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +9.9914929801701369D-01, +3.1676511180119155D-01
     &    /)
      H(1:48,36) = (/
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +9.9801390072069940D-01, +9.9801390072069940D-01,
     & +7.2100018712984359D-01, +7.2100018712984359D-01,
     & +9.4910149657117848D-01, +9.4910149657117848D-01,
     & +3.2732683535398766D-01, +3.2732683535398766D-01,
     & +9.6773340156674170D-01, +9.6773340156674170D-01,
     & +4.7140452079103129D-01, +4.7140452079103129D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +7.1309424437485402D-01, +7.1309424437485402D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +7.5141589705045231D-01, +7.5141589705045231D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +1.4285714285714138D-01, +1.4285714285714138D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +6.2133277860475644D-01, +6.2133277860475644D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +3.1676511180119155D-01, +9.9914929801701369D-01
     &    /)
      H(1:48,37) = (/
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +9.4910149657117848D-01, +9.4910149657117848D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.9801390072069940D-01, +9.9801390072069940D-01,
     & +9.9203174552379325D-01, +9.9203174552379325D-01,
     & +9.6773340156674170D-01, +9.6773340156674170D-01,
     & +8.4414751910646835D-01, +8.4414751910646835D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +9.4401083817138132D-01, +9.4401083817138132D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +9.9914929801701369D-01, +9.9914929801701369D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +9.7851059943021512D-01, +9.7851059943021512D-01,
     & +9.2857142857142860D-01, +7.8571428571428570D-01
     &    /)
      H(1:48,38) = (/
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +3.2732683535398766D-01, +3.2732683535398766D-01,
     & +5.7390337110447343D-01, +5.7390337110447343D-01,
     & +7.2100018712984359D-01, +7.2100018712984359D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +4.7140452079103129D-01, +4.7140452079103129D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +7.1309424437485391D-01, +7.1309424437485391D-01,
     & +7.1309424437485402D-01, +7.1309424437485402D-01,
     & +8.9118891772442388D-01, +8.9118891772442388D-01,
     & +1.4285714285714060D-01, +1.4285714285714060D-01,
     & +4.2056004125370633D-01, +4.2056004125370633D-01,
     & +7.5141589705045231D-01, +7.5141589705045231D-01,
     & +9.1100602236709483D-01, +9.1100602236709483D-01,
     & +3.1676511180119155D-01, +3.1676511180119155D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +6.2133277860475644D-01, +6.2133277860475644D-01,
     & +7.8571428571428570D-01, +9.2857142857142860D-01
     &    /)
      H(1:48,39) = (/
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +9.4910149657117848D-01, +9.4910149657117848D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +9.9801390072069940D-01, +9.9801390072069940D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +9.6773340156674170D-01, +9.6773340156674170D-01,
     & +9.9203174552379325D-01, +9.9203174552379325D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +8.4414751910646835D-01, +8.4414751910646835D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +9.4401083817138132D-01, +9.4401083817138132D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +9.9914929801701369D-01, +9.9914929801701369D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +9.7851059943021512D-01, +6.2133277860475644D-01
     &    /)
      H(1:48,40) = (/
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +5.7390337110447343D-01, +5.7390337110447343D-01,
     & +3.2732683535398766D-01, +3.2732683535398766D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +7.2100018712984359D-01, +7.2100018712984359D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +4.7140452079103129D-01, +4.7140452079103129D-01,
     & +7.1309424437485391D-01, +7.1309424437485391D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +8.9118891772442388D-01, +8.9118891772442388D-01,
     & +7.1309424437485402D-01, +7.1309424437485402D-01,
     & +4.2056004125370633D-01, +4.2056004125370633D-01,
     & +1.4285714285714060D-01, +1.4285714285714060D-01,
     & +9.1100602236709483D-01, +9.1100602236709483D-01,
     & +7.5141589705045231D-01, +7.5141589705045231D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +3.1676511180119155D-01, +3.1676511180119155D-01,
     & +7.8571428571428570D-01, +7.8571428571428570D-01,
     & +6.2133277860475644D-01, +9.7851059943021512D-01
     &    /)
      H(1:48,41) = (/
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.9203174552379325D-01, +9.9203174552379325D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +8.9752746785575066D-01, +8.9752746785575066D-01,
     & +5.7390337110447343D-01, +5.7390337110447343D-01,
     & +9.8198050606196563D-01, +9.8198050606196563D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +9.9914929801701369D-01, +9.9914929801701369D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +9.4401083817138132D-01, +9.4401083817138132D-01,
     & +4.2056004125370633D-01, +4.2056004125370633D-01,
     & +9.7851059943021512D-01, +9.7851059943021512D-01,
     & +7.8571428571428570D-01, +7.8571428571428570D-01,
     & +8.4414751910646835D-01, +8.4414751910646835D-01,
     & +7.1309424437485391D-01, +7.1309424437485391D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +9.1100602236709471D-01, +9.1100602236709471D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +8.9118891772442388D-01, +7.1309424437485402D-01
     &    /)
      H(1:48,42) = (/
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.6773340156674170D-01, +9.6773340156674170D-01,
     & +4.7140452079103107D-01, +4.7140452079103107D-01,
     & +9.4910149657117848D-01, +9.4910149657117848D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.9801390072069940D-01, +9.9801390072069940D-01,
     & +7.2100018712984359D-01, +7.2100018712984359D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +3.1676511180119155D-01, +3.1676511180119155D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +1.4285714285714060D-01, +1.4285714285714060D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +6.2133277860475644D-01, +6.2133277860475644D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +7.5141589705045231D-01, +7.5141589705045231D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +7.1309424437485402D-01, +8.9118891772442388D-01
     &    /)
      H(1:48,43) = (/
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.9203174552379325D-01, +9.9203174552379325D-01,
     & +5.7390337110447343D-01, +5.7390337110447343D-01,
     & +8.9752746785575066D-01, +8.9752746785575066D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +9.8198050606196563D-01, +9.8198050606196563D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +9.9914929801701369D-01, +9.9914929801701369D-01,
     & +4.2056004125370633D-01, +4.2056004125370633D-01,
     & +9.4401083817138132D-01, +9.4401083817138132D-01,
     & +7.8571428571428570D-01, +7.8571428571428570D-01,
     & +9.7851059943021512D-01, +9.7851059943021512D-01,
     & +7.1309424437485391D-01, +7.1309424437485391D-01,
     & +8.4414751910646835D-01, +8.4414751910646835D-01,
     & +9.1100602236709471D-01, +9.1100602236709471D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +8.9118891772442388D-01, +8.9118891772442388D-01,
     & +9.2857142857142860D-01, +9.9231742781784316D-01
     &    /)
      H(1:48,44) = (/
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +4.7140452079103107D-01, +4.7140452079103107D-01,
     & +9.6773340156674170D-01, +9.6773340156674170D-01,
     & +3.2732683535398799D-01, +3.2732683535398799D-01,
     & +9.4910149657117848D-01, +9.4910149657117848D-01,
     & +7.2100018712984359D-01, +7.2100018712984359D-01,
     & +9.9801390072069940D-01, +9.9801390072069940D-01,
     & +3.1676511180119155D-01, +3.1676511180119155D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +1.4285714285714060D-01, +1.4285714285714060D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +6.2133277860475644D-01, +6.2133277860475644D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +4.9999999999999956D-01, +4.9999999999999956D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +7.5141589705045231D-01, +7.5141589705045231D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +7.1309424437485402D-01, +7.1309424437485402D-01,
     & +9.9231742781784316D-01, +9.2857142857142860D-01
     &    /)
      H(1:48,45) = (/
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +3.2732683535398766D-01, +3.2732683535398766D-01,
     & +9.9203174552379325D-01, +9.9203174552379325D-01,
     & +4.7140452079103129D-01, +4.7140452079103129D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +7.2100018712984359D-01, +7.2100018712984359D-01,
     & +9.4401083817138132D-01, +9.4401083817138132D-01,
     & +1.4285714285714060D-01, +1.4285714285714060D-01,
     & +9.9914929801701369D-01, +9.9914929801701369D-01,
     & +3.1676511180119155D-01, +3.1676511180119155D-01,
     & +8.4414751910646835D-01, +8.4414751910646835D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +9.7851059943021512D-01, +9.7851059943021512D-01,
     & +6.2133277860475644D-01, +6.2133277860475644D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +7.1309424437485402D-01, +7.1309424437485402D-01,
     & +1.0000000000000000D+00, +1.0000000000000000D+00,
     & +7.5141589705045231D-01, +9.1100602236709483D-01
     &    /)
      H(1:48,46) = (/
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.4910149657117848D-01, +9.4910149657117848D-01,
     & +5.7390337110447343D-01, +5.7390337110447343D-01,
     & +9.6773340156674170D-01, +9.6773340156674170D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.9801390072069940D-01, +9.9801390072069940D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +4.2056004125370633D-01, +4.2056004125370633D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +7.1309424437485391D-01, +7.1309424437485391D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +7.8571428571428570D-01, +7.8571428571428570D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +8.9118891772442388D-01, +8.9118891772442388D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +9.1100602236709483D-01, +7.5141589705045231D-01
     &    /)
      H(1:48,47) = (/
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +3.2732683535398766D-01, +3.2732683535398766D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +4.7140452079103129D-01, +4.7140452079103129D-01,
     & +9.9203174552379325D-01, +9.9203174552379325D-01,
     & +7.2100018712984359D-01, +7.2100018712984359D-01,
     & +9.8198050606196574D-01, +9.8198050606196574D-01,
     & +1.4285714285714060D-01, +1.4285714285714060D-01,
     & +9.4401083817138132D-01, +9.4401083817138132D-01,
     & +3.1676511180119155D-01, +3.1676511180119155D-01,
     & +9.9914929801701369D-01, +9.9914929801701369D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +8.4414751910646835D-01, +8.4414751910646835D-01,
     & +6.2133277860475644D-01, +6.2133277860475644D-01,
     & +9.7851059943021512D-01, +9.7851059943021512D-01,
     & +7.1309424437485402D-01, +7.1309424437485402D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +7.5141589705045231D-01, +7.5141589705045231D-01,
     & +1.0000000000000000D+00, +9.6890428330360967D-01
     &    /)
      H(1:48,48) = (/
     & +9.4491118252306805D-01, +9.4491118252306805D-01,
     & +8.3808170984752572D-01, +8.3808170984752572D-01,
     & +8.9973541084243724D-01, +8.9973541084243724D-01,
     & +7.5592894601845440D-01, +7.5592894601845440D-01,
     & +1.8898223650461299D-01, +1.8898223650461299D-01,
     & +9.9402979738800490D-01, +9.9402979738800490D-01,
     & +5.7390337110447343D-01, +5.7390337110447343D-01,
     & +9.4910149657117848D-01, +9.4910149657117848D-01,
     & +6.5465367070797686D-01, +6.5465367070797686D-01,
     & +9.6773340156674170D-01, +9.6773340156674170D-01,
     & +8.9752746785575055D-01, +8.9752746785575055D-01,
     & +9.9801390072069940D-01, +9.9801390072069940D-01,
     & +4.2056004125370633D-01, +4.2056004125370633D-01,
     & +9.6890428330360967D-01, +9.6890428330360967D-01,
     & +4.9999999999999967D-01, +4.9999999999999967D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +7.1309424437485391D-01, +7.1309424437485391D-01,
     & +9.2857142857142849D-01, +9.2857142857142849D-01,
     & +7.8571428571428570D-01, +7.8571428571428570D-01,
     & +9.2857142857142860D-01, +9.2857142857142860D-01,
     & +8.9118891772442388D-01, +8.9118891772442388D-01,
     & +9.9231742781784316D-01, +9.9231742781784316D-01,
     & +9.1100602236709483D-01, +9.1100602236709483D-01,
     & +9.6890428330360967D-01, +1.0000000000000000D+00
     &    /)
c
      return
      end subroutine

