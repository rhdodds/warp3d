c
c **********************************************************************
c *                                                                    *
c *    mm10.f                                                          *
c *                                                                    *
c *         written by : mcm                                           *
c *                                                                    *
c *         Slip rate and hardening functions                          *
c *         Updated extensively by TJT and RHD                         *
c *                                                                    *
c **********************************************************************
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
      subroutine mm10b_unknown_hard_error( props )
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

      end
c
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
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     *                         Form R2                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formR2( props, np1, n, vec1, vec2, 
     &                        stress, tt, R2, gp )
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      integer :: gp
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: R2, h
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1,vec2
c
      integer :: len
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, R2:64
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type )
        case( 1 ) ! Voche
           call mm10_h_voche( props, np1, n, stress, tt, h ) 
        case( 2 ) ! MTS
           call mm10_h_mts( props, np1, n, stress, tt, h )
        case( 3 ) ! User
           call mm10_h_user( props, np1, n, stress, tt, h )
        case( 4 ) ! ORNL
           call mm10_h_ornl( props, np1, n, vec1, vec2, stress,
     &                       tt, h, gp )
        case( 7 ) !MRR
           call mm10_h_mrr( props, np1, n, vec1, vec2, stress,
     &                      tt, h, gp )
        case( 8 ) !Armstrong-Frederick
           call mm10_h_arfr( props, np1, n, vec1, vec2, stress,
     &                      tt, h )
        case( 9 ) ! DJGM
           call mm10_h_DJGM( props, np1, n, vec1, vec2, stress,
     &                      tt, h )
        case default
           call mm10b_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
c              compute time rate of change of hardening parameters, 
c              store for possible prediction during next time/load step
c
      len = props%num_hard
      R2(1:len) = tt(1:len) - h(1:len)
      np1%tt_rate(1:len) = ( h(1:len) - n%tau_tilde(1:len) ) / np1%tinc
c      
      return
      end 
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formR2i                      *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/18/2016 rhd              *
c     *                                                              *
c     *                       Form R2                                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formR2i( props, np1, n, ivec1, ivec2, 
     &                         stress, tt, R2 )
      use iso_Fortran_env
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      complex(kind=real64), dimension(6) :: stress
      complex(kind=real64), dimension(size_num_hard) :: R2
      complex(kind=real64), dimension(size_num_hard) :: tt
      complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
c
c              locals - automatics
c
      complex(kind=real64), dimension(size_num_hard) :: h
c!DIR$ ASSUME_ALIGNED ivec1:64, ivec2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, R2:64
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type )
        case( 1, 2, 3, 8, 9 ) ! Voche, MTS, USER, Armstrong-Frederick, DJGM
          continue
        case( 4 )  ! ORNL
          call mm10_hi_ornl( props, np1, n, ivec1, ivec2, stress, tt, h)
        case( 7 ) ! MRR
          call mm10_hi_mrr( props, np1, n, ivec1, ivec2, stress, tt, h)
        case default
          call mm10b_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
      R2 = tt - h
c      
      return
      end 
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ11                      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 1/15/2016 rhd              *
c     *                                                              *
c     *     Form the stress varying with stress part                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ11( props, np1, n, vec1, vec2,
     &                         arr1, arr2, stress, tt, J11 )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(6,6) :: J11
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: 
     &                             arr1,arr2
c
      integer :: i
      integer, save :: passno = 0
      logical :: debug
      double precision, allocatable :: symtqmat(:,:), dgammadtau(:)
      double precision :: wp(3), work_vec(6), Iw(6,6)
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, arr1:64, arr2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, J11:64
c
      passno = passno + 1
      debug = .false.
      if( debug ) then
          write(props%out,*) "In mm10_formJ11"
          write(props%out,*) "props%nslip, size_nslip, passno: ",
     &                        props%nslip, size_nslip, passno
      end if
c
      allocate( symtqmat(6,size_nslip) )
      allocate( dgammadtau(size_nslip) )
c      
      J11 = zero
      call mm10_symSWmat( stress, np1%qc, props%nslip, symtqmat )
c
c              generalization of CP model implementation for other 
c              slip rate equations, requiring other forms of
c              d_gamma/d_tau.
c              vector dgammadtau should be 1 x n_slip
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type )
       case( 1 ) ! voche
        call mm10_dgdt_voche( props,np1, n, stress, tt, dgammadtau )
       case( 2 ) ! MTS
        call mm10_dgdt_mts( props, np1,n, stress, tt, dgammadtau )
       case( 3 ) ! User
        call mm10_dgdt_user( props,np1, n, stress, tt, dgammadtau )
       case( 4 ) ! ORNL
        call mm10_dgdt_ornl( props,np1, n, stress, tt, dgammadtau )
       case( 7 ) ! MRR
        call mm10_dgdt_mrr( props,np1, n, stress, tt, dgammadtau )
       case( 8 ) ! Armstrong-Frederick
        call mm10_dgdt_arfr( props,np1, n, stress, tt, dgammadtau )
       case( 9 ) ! DJGM
        call mm10_dgdt_DJGM( props,np1, n, stress, tt, dgammadtau )
       case default
        call mm10b_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
      do i = 1, props%nslip
       call mm10_b_mult_type_4( work_vec, props%stiffness, np1%ms(1,i),
     &                          symtqmat(1,i), two )   
       call DGER( 6, 6, dgammadtau(i), work_vec, 1, np1%ms(1,i),
     &            1, J11, 6 )
     
c        call DGER(6,6,dgammadtau(i),
c     &      matmul(props%stiffness, np1%ms(1:6,i))
c     &      + 2.0d0*symtqmat(1:6,i), 1, np1%ms(1:6,i),1,J11,6)
      end do
      deallocate( symtqmat )
      deallocate( dgammadtau )
c
      call mm10_form_wp( props, np1, n, vec1, vec2, stress, tt, wp )
      call mm10_IW (wp, Iw )
      J11 = J11 + Iw
c
c      add J of creep strain due to pressure precipitation (Ran)
      if ( abs(props%cp_031-one)<1.0e-5 ) then
            call mm10_halite_formJpp( props, J11, np1%tinc )
      end if
c
      do i = 1, 6
        J11(i,i) = J11(i,i) + one
      end do
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ12                      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/14/2016 rhd              *
c     *                                                              *
c     *     Form the stress varying with hardening part              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ12( props, np1, n, vec1, vec2, arr1,
     &                         arr2, stress, tt, J12 )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      double precision, dimension(6) :: stress
      double precision, dimension(6,size_num_hard) :: J12
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) ::
     &     arr1,arr2
c
      integer :: i, len
      logical :: debug
      double precision, dimension(6) :: tempmat
      double precision :: work_vec(max_uhard)
c
c              automatics
c      
      double precision, dimension(6,size_nslip) :: symtqmat
      double precision, dimension(size_nslip,size_num_hard)
     &                  :: dgammadtt
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, arr1:64, arr2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, J12:64     
c
      debug = .false.
      J12 = zero
      call mm10_symSWmat( stress, np1%qc, props%nslip, symtqmat )
c
c
c              generalization of CP model implementation for other 
c              slip rate equations, requiring other forms of 
c              d_gamma/d_hardening. Vector dgammadtt should be 
c              n_slip x n_hardening, which is the
c              derivative of slip rate alpha with respect to 
c              hardening variable beta.
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type )
       case( 1 )  ! voche
        call mm10_dgdh_voche( props,np1, n, stress, tt, dgammadtt )
       case( 2 )  ! MTS
        call mm10_dgdh_mts( props, np1,n, stress, tt, dgammadtt )
       case( 3 ) ! User
        call mm10_dgdh_user( props,np1, n, stress, tt, dgammadtt )
       case( 4 ) ! ORNL
        call mm10_dgdh_ornl( props,np1, n, stress, tt, dgammadtt )
       case( 7 ) ! MRR
        call mm10_dgdh_mrr( props,np1, n, stress, tt, dgammadtt )
       case( 8 ) ! Armstrong-Frederick
        call mm10_dgdh_arfr( props,np1, n, stress, tt, dgammadtt )
       case( 9)  ! DJGM
        call mm10_dgdh_DJGM( props,np1, n, stress, tt, dgammadtt )
       case default 
        call mm10b_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
      len = props%num_hard
c      
      do i = 1, props%nslip
        call mm10_b_mult_type_4( tempmat, props%stiffness, np1%ms(1,i),
     &                           symtqmat(1,i), two )             
        work_vec(1:len) = dgammadtt(i,1:len) 
        call DGER( 6, len, one, tempmat, 1, work_vec, 1, J12, 6 )
      end do  
c
c              original code
c
c        tempmat(1:6) = matmul(props%stiffness, np1%ms(1:6,i))
c     &      + two*symtqmat(1:6,i) 
c        call DGER( 6, props%num_hard, one,
c     &      tempmat(1:6), 1, dgammadtt(i,1:props%num_hard),1,
c     &      J12(1:6,1:props%num_hard),6)
c
      return
      end 


c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ21                      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     *     Form the hardening varying with stress part              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ21( props, np1, n, vec1, vec2,
     &                         arr1, arr2, stress, tt, J21 )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard,6) :: J21
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
c              automatics
c
      double precision, dimension(size_num_hard,6) :: estr
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, arr1:64, arr2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, J21:64     
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type ) 
        case( 1 ) ! Voche
          call mm10_estress_voche( props, np1, n, stress, tt, estr )
        case( 2 )  ! MTS
          call mm10_estress_mts( props, np1, n, stress, tt, estr )
        case( 3 ) ! User
          call mm10_estress_user( props, np1, n, stress, tt, estr )
        case( 4 ) ! ORNL
          call mm10_estress_ornl( props, np1, n, vec1, vec2, arr1, 
     &                            arr2, stress, tt, estr )
        case( 7 ) ! MRR
          call mm10_estress_mrr( props, np1, n, vec1, vec2, arr1,
     &                           arr2, stress, tt, estr )
        case( 8 ) ! Armstrong-Frederick
          call mm10_estress_arfr( props, np1, n, vec1, vec2, arr1,
     &                           arr2, stress, tt, estr )
        case( 9 ) ! DJGM
          call mm10_estress_DJGM( props, np1, n, vec1, vec2, arr1,
     &                            arr2, stress, tt, estr )
        case default
          call mm10b_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
      J21 = -estr
c
      return
      end 
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ22                      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/16/2016 rhd              *
c     *                                                              *
c     *     Form the hardening varying with hardening part           *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ22( props, np1, n, vec1, vec2,
     &                         arr1, arr2, stress, tt, J22 )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, 
     &       dimension(size_num_hard,size_num_hard) :: J22
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
c              automatics
c
      double precision, 
     &      dimension(size_num_hard,size_num_hard) :: etau
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, arr1:64, arr2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, J22:64     
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type ) 
        case( 1 ) ! Voche
          call mm10_ehard_voche( props, np1, n, stress, tt, etau )
        case( 2 ) ! MTS
          call mm10_ehard_mts( props, np1, n, stress, tt, etau )
        case( 3 ) ! User
          call mm10_ehard_user( props, np1, n, stress, tt, etau )
        case( 4 ) ! ORNL
          call mm10_ehard_ornl( props, np1, n, vec1, vec2, arr1, arr2,
     &                          stress, tt, etau )
        case( 7 ) ! MRR
          call mm10_ehard_mrr( props, np1, n, vec1, vec2, arr1, arr2,
     &                         stress, tt, etau )
        case( 8 ) ! Armstrong-Frederick
          call mm10_ehard_arfr( props, np1, n, vec1, vec2, arr1, arr2,
     &                         stress, tt, etau )
        case( 9 ) ! DJGM
          call mm10_ehard_DJGM( props, np1, n, vec1, vec2, arr1, arr2,
     &                          stress, tt, etau )
        case default
          call mm10b_unknown_hard_error(props)
      end select
c ****** END: Add new Constitutive Models into this block ******
c
      J22 = etau
      return
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formvecs                     *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     *     Form intermediate vectors which are repeatedly used by   *
c     *     other constitutive routines (e.g. precompute slip rates) *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formvecs( props, np1, n, stress, tt, vec1, vec2 )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1,vec2
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, stress:64, tt:64
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type )
        case( 1, 2, 3 ) ! Voche, MTS, User
         vec1 = zero
         vec2 = zero
        case( 4 ) ! ORNL
          call mm10_v_ornl( props, np1, n, stress, tt, vec1, vec2 )
        case( 7 ) ! MRR
          call mm10_v_mrr( props, np1, n, stress, tt, vec1, vec2 )
        case( 8 ) ! Armstrong-Frederick
          call mm10_v_arfr( props, np1, n, stress, tt, vec1, vec2 )
        case( 9 ) ! DJGM
          call mm10_v_DJGM( props, np1, n, stress, tt, vec1, vec2 )
        case default
          call mm10b_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
      return
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formarrs                     *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     *     Form intermediate arrays which are repeatedly used by    *
c     *     other constitutive routines                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formarrs( props, np1, n, stress, tt, vec1, vec2,
     &                          arr1, arr2, both )
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      integer :: both
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1,vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1,arr2
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type )
        case( 1, 2, 3 ) ! Voche, MTS, User
          continue
        case( 4 ) ! ORNL
          call mm10_a_ornl( props, np1, n, stress, tt, vec1, vec2,
     &                      arr1, arr2, both )
        case( 7 ) ! MRR
          call mm10_a_mrr( props, np1, n, stress, tt, vec1, vec2,
     &                     arr1, arr2, both )
        case( 8 ) ! Armstrong-Frederick
          call mm10_a_arfr( props, np1, n, stress, tt, vec1, vec2,
     &                     arr1, arr2, both )
        case( 9 ) ! DJGM
         call mm10_a_DJGM( props, np1, n, stress, tt, vec1, vec2,
     &                     arr1, arr2, both)
        case default
         call mm10b_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
      return      
c
      end

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
      subroutine mm10_formJ11i( props, np1, n, ivec1, ivec2,
     &                          stress, tt, J11)
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(6,6) :: J11
      double precision, dimension(size_num_hard) :: tt
      complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
c
      integer :: k
      logical :: debug
c
      complex(kind=real64), dimension(6) :: A, Ri
      double precision :: h, hi, zeroA(6)
      complex(kind=real64) :: i1
c
c              automatics
c
      complex(kind=real64), dimension(size_num_hard) :: B
      double precision, dimension(size_num_hard) :: zeroB
c!DIR$ ASSUME_ALIGNED ivec1:64, ivec2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, J11:64     
c
      debug = .false.
      if( debug ) write (props%out,*) "In mm10_formJ11i"
      h = 1.0d-12
      i1 = (zero, one)
c
      J11   = zero
      zeroA = zero
      zeroB = zero
c      
      do k = 1, 6
       A = dcmplx( stress, zeroA )
       B = dcmplx( tt, zeroB )
       if( stress(k) == zero) then
         hi = h
       else
         hi = h*dabs(stress(k))
       end if
       A(k) = A(k) + i1*hi
       call mm10_formvecsi( props, np1, n, A, B, ivec1, ivec2 )
       call mm10_formR1i( props, np1, n, ivec1, ivec2, A, B, Ri )
       J11(1:6,k) = one / hi*aimag(Ri)
      end do
c
      return
      end 
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ12i                     *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/18/2016 rhd              *
c     *                                                              *
c     *     Form the stress varying with hardening part              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ12i( props, np1, n, ivec1, ivec2,
     &                          stress, tt, J12 )
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c      
      double precision, dimension(6) :: stress
      double precision, dimension(6,size_num_hard) :: J12
      double precision, dimension(size_num_hard) :: tt
      complex(kind=real64), dimension(max_uhard) :: ivec1,ivec2
c
      integer :: k
      logical :: debug
c
      complex(kind=real64), dimension(6) :: Ri, A
      double precision, dimension(6) :: zeroA
      double precision :: h, hi
      complex(kind=real64) :: i1
c
c              automatics
c
      complex(kind=real64), dimension(size_num_hard) :: B
      double precision, dimension(size_num_hard) :: zeroB
c!DIR$ ASSUME_ALIGNED ivec1:64, ivec2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, J12:64     
c
      debug = .false.
      if( debug ) write (props%out,*) "In mm10_formJ12i"
      h = 1.0d-12
      i1 = (zero, one)
c
      J12   = zero
      zeroA = zero
      zeroB = zero
c      
      do k = 1,props%num_hard
        A = dcmplx (stress, zeroA)
        B = dcmplx (tt, zeroB)
        if( tt(k) == zero ) then
          hi = h
        else
          hi = h*dabs(tt(k))
        end if
        B(k) = B(k) + i1*hi
        call mm10_formvecsi( props, np1, n, A, B, ivec1, ivec2 )
        call mm10_formR1i( props, np1, n, ivec1, ivec2, A, B, Ri )
        J12(1:6,k) = one/ hi*aimag(Ri)
      end do
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ21i                     *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/18/2016 rhd              *
c     *                                                              *
c     *     Form the hardening varying with stress part              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ21i( props, np1, n, ivec1, ivec2,
     &                          stress, tt, J21 )
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard,6) :: J21
      double precision, dimension(size_num_hard) :: tt
      complex(kind=real64), dimension(max_uhard) :: ivec1,ivec2
      integer :: k
      logical :: debug
c
      complex(kind=real64), dimension(6) :: A
      double precision :: h, hi, zeroA(6)
      complex(kind=real64) :: i1
c
c              automatics
c
      complex(kind=real64), dimension(size_num_hard) :: B, Ri
      double precision, dimension(size_num_hard) :: zeroB
c!DIR$ ASSUME_ALIGNED ivec1:64, ivec2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, J21:64     
c
      debug = .false.
      if( debug ) write (props%out,*) "In mm10_formJ121i"
      h = 1.0d-12
      i1 = (zero, one)
c
      J21   = zero
      zeroA = zero
      zeroB = zero
      do k = 1, 6
        A = dcmplx (stress, zeroA)
        B = dcmplx (tt, zeroB)
        if( stress(k) == zero ) then
          hi = h
        else
          hi = h*dabs(stress(k))
        endif
        A(k) = A(k) + i1*hi
        call mm10_formvecsi( props, np1, n, A, B, ivec1, ivec2 )
        call mm10_formR2i( props, np1, n, ivec1, ivec2, A, B, Ri )
        J21(1:props%num_hard,k) = one / hi*aimag(Ri)
      end do
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJ22i                     *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/18/2016 rhd              *
c     *                                                              *
c     *     Form the hardening varying with hardening part           *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ22i( props, np1, n, ivec1, ivec2,
     &                          stress, tt, J22 )
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, 
     &   dimension(size_num_hard,size_num_hard) :: J22
      double precision, dimension(size_num_hard) :: tt
      complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
c
      integer :: k
      logical :: debug
c
      complex(kind=real64), dimension(6) :: A
      double precision :: h, hi, zeroA(6)
      complex(kind=real64) :: i1
c
c              automatics
c
      double precision, dimension(size_num_hard) :: zeroB
      complex(kind=real64), dimension(size_num_hard) :: Ri, B
c!DIR$ ASSUME_ALIGNED ivec1:64, ivec2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, J22:64     
c
      debug = .false.
      if( debug ) write(props%out,*) "In mm10_formJ22i"
      h = 1.0d-12
      i1 = (zero, one)
c
      J22   = zero
      zeroA = zero
      zeroB = zero
c      
      do k = 1, props%num_hard
        A = dcmplx (stress, zeroA)
        B = dcmplx (tt, zeroB)
        if( tt(k) == zero ) then
          hi = h
        else
          hi = h*dabs(tt(k))
        end if
        B(k) = B(k) + i1*hi
        call mm10_formvecsi( props, np1, n, A, B, ivec1, ivec2 )
        call mm10_formR2i( props, np1, n, ivec1, ivec2, A, B, Ri )
        J22(1:props%num_hard,k) = one / hi*aimag(Ri)
      end do
c
      return
      end 
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formvecsi                    *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     *     Form intermediate vectors which are repeatedly used by   *
c     *     other constitutive routines (e.g. precompute slip rates) *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formvecsi( props, np1, n, stress, tt, 
     &                           ivec1, ivec2 )
      use iso_Fortran_env
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      complex(kind=real64), dimension(6) :: stress
      complex(kind=real64), dimension(size_num_hard) :: tt
      complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type )
       case( 1, 2, 3, 8, 9 )  ! voche, MTS, User
         continue
       case( 4 ) ! ORNL
        call mm10_vi_ornl( props, np1, n, stress, tt, ivec1, ivec2 )
       case( 7 ) ! MRR
        call mm10_vi_mrr( props, np1, n, stress, tt, ivec1, ivec2 )
       case default
        call mm10b_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
      return
c
      end 
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJnew                     *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/15/2016                  *
c     *                                                              *
c     *     Form the jacobian from lower subroutines                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJ( props, np1, n, vec1, vec2, arr1, arr2,
     &                       stress, tt, J )
      use mm10_defs
      implicit none
c
c              parameters
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision,
     & dimension(six_plus_num_hard,six_plus_num_hard) :: J
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
c              locals
c
      integer :: len
      double precision :: local_J11(6,6), local_J12(6,max_uhard)
      double precision, allocatable, dimension(:,:) :: local_J21,
     &                                                 local_J22
c
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, arr1:64, arr2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, J:64
c
c
c              use local arrays with size of J submatrices to 
c              support calls.  no copyin req'd. only copyout
c
      len = props%num_hard
      allocate( local_J21(len,6), local_J22(len,len) )
c
      call mm10_formarrs( props, np1, n, stress, tt, vec1, vec2,
     &                    arr1, arr2, 2 )
c
      call mm10_formJ11( props, np1, n, vec1, vec2, arr1, arr2,
     &                   stress, tt, local_J11 ) ! J(1:6,1:6))
      call mm10_formJ12( props, np1, n, vec1, vec2, arr1, arr2,
     &                   stress, tt, local_J12 ) !
      call mm10_formJ21( props, np1, n, vec1, vec2, arr1, arr2,
     &                   stress, tt, local_J21 ) !
      call mm10_formJ22( props, np1, n, vec1, vec2, arr1, arr2,
     &                   stress, tt, local_J22 ) 
c
      J(1:6,1:6)         = local_J11
      J(1:6,7:6+len)     = local_J12(1:6,1:len)
      J(7:6+len,1:6)     = local_J21(1:len,1:6)
      J(7:6+len,7:6+len) = local_J22(1:len,1:len)
c
      deallocate( local_J21, local_J22 )
c
      return
      end    

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formJi                       *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     *     Form the jacobian from lower subroutines                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formJi( props, np1, n, ivec1, 
     &                        ivec2, stress, tt, J )
      use iso_Fortran_env
      use mm10_defs
      implicit none
c
c              parameters
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision,
     & dimension(six_plus_num_hard,six_plus_num_hard) :: J
      double precision, dimension(size_num_hard) :: tt
      complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
c
c              locals
c
      integer :: len
      double precision :: local_J11(6,6), local_J12(6,max_uhard)
      double precision, allocatable, dimension(:,:) :: local_J21,
     &                                                 local_J22
c
c!DIR$ ASSUME_ALIGNED ivec1:64, ivec2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, J:64
c
      len = props%num_hard
      allocate( local_J21(len,6), local_J22(len,len) )
c
      call mm10_formJ11i( props, np1, n, ivec1, ivec2, 
     &                    stress, tt, local_J11 )
      call mm10_formJ12i( props, np1, n, ivec1, ivec2, 
     &                    stress, tt, local_J12 )
      call mm10_formJ21i( props, np1, n, ivec1, ivec2, 
     &                    stress, tt, local_J21 )
      call mm10_formJ22i( props, np1, n, ivec1, ivec2, 
     &                    stress, tt, local_J22 )
c
      J(1:6,1:6)         = local_J11
      J(1:6,7:6+len)     = local_J12(1:6,1:len)
      J(7:6+len,1:6)     = local_J21(1:len,1:6)
      J(7:6+len,7:6+len) = local_J22(1:len,1:len)
c
      deallocate( local_J21, local_J22 )
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
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     *     Form the residual from lower subroutines                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formR( props, np1, n, vec1, vec2, stress, tt, 
     &                       R, gp)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      integer :: gp      
      double precision, dimension(6) :: stress
      double precision, dimension(six_plus_num_hard) :: R
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, stress:64, tt:64, R:64
c
      call mm10_formvecs( props, np1, n, stress, tt, vec1, vec2 )
      call mm10_formR1( props, np1, n, vec1, vec2, stress, 
     &                  tt, R(1), gp )
      call mm10_formR2( props, np1, n, vec1, vec2, stress, tt,
     &                  R(7), gp )
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formR1                       *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     *     Form R1                                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formR1( props, np1, n, vec1, vec2,
     &                        stress, tt, R1, gp )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(6) :: R1
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
c
      integer :: gp
      double precision :: dbarp(6), wp(3), symTW(6), work_vec1(6),
     &                    work_vec2(6) 
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, stress:64, tt:64, R1:64
c
      call mm10_form_dbarp( props, np1, n, vec1, vec2, 
     &                      stress, tt, dbarp )
      call mm10_form_wp( props, np1, n, vec1, vec2, stress, tt, wp )
      call mm10_symSW( stress, wp, symTW) 
c
      work_vec1 =  np1%D - dbarp
      if ( abs(props%cp_031-one)<1.0e-5 ) then
        call mm10_halite_formRpp( props, work_vec1, stress, np1%tinc )
      end if
      call mm10_b_mult_type_2( work_vec2, props%stiffness, work_vec1 )
      R1 = stress - n%stress - work_vec2 + two * symTW
c
c            original code
c      
c      R1 = stress - n%stress - matmul(props%stiffness, work_vec1) 
c     &      + two * symTW
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_formR1i                      *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     *     Form R1                                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_formR1i( props, np1, n, ivec1, ivec2,
     &                         stress, tt, R1 )
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n  
      complex(kind=real64), dimension(6) :: stress
      complex(kind=real64), dimension(6) :: R1
      complex(kind=real64), dimension(size_num_hard) :: tt
      complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
c
      complex(kind=real64), dimension(6) :: dbarp, temp
      complex(kind=real64), dimension(6,6) :: stiff2
      double precision, dimension(6,6) :: zeroff
      complex(kind=real64), dimension(3) :: wp
      complex(kind=real64), dimension(6) :: symTW
c!DIR$ ASSUME_ALIGNED ivec1:64, ivec2:64, stress:64, tt:64, R1:64
c
      call mm10_form_dbarpi( props, np1, n, ivec1, ivec2, 
     &                       stress, tt, dbarp )
      call mm10_form_wpi( props, np1, n, ivec1, ivec2 ,stress, tt, wp )
      call mm10_symSWi( stress, wp, symTW )
c
      zeroff = zero
      stiff2 = dcmplx( props%stiffness, zeroff )
      temp = np1%D - dbarp
      R1 = stress - n%stress - matmul(stiff2, np1%D - dbarp) 
     &      + two * symTW
c
      return
      end 

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_form_dbarp                   *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_form_dbarp( props, np1, n, vec1, vec2,
     &                            stress, tt, dbar )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(6) :: dbar
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
c
      integer :: i, nslip
      double precision :: slipinc, rs
      double precision, external :: mm10_rs
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, dbar:64
c
      nslip = props%nslip
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type )
        case( 1 ) ! Voche
          dbar = zero
          do i = 1, nslip
            call mm10_slipinc( props, np1, n, stress, tt, i, slipinc )
            rs = mm10_rs( props, np1, n, stress, tt, i )
            dbar = dbar + (rs*np1%tinc*props%iD_v + slipinc) *
     &             np1%ms(1:6,i)
          end do
        case( 2 ) ! MTS     
          dbar = zero
          do i = 1, nslip
            call mm10_slipinc( props, np1, n, stress, tt, i, slipinc )
            dbar = dbar + slipinc*np1%ms(1:6,i)
          end do
        case( 3 ) ! User
          dbar = zero
          do i = 1, nslip
            call mm10_slipinc_user( props, np1, n, stress, tt, i,
     &                              slipinc )
            dbar = dbar + slipinc*np1%ms(1:6,i)
          end do
        case( 4 ) ! ORNL
          call mm10_b_mult_type_2a( dbar, np1%ms(1,1), vec1, nslip  )
        case( 7 ) ! MRR
          call mm10_b_mult_type_2a( dbar, np1%ms(1,1), vec1, nslip  )
        case( 8 ) ! Armstrong-Frederick
          call mm10_b_mult_type_2a( dbar, np1%ms(1,1), vec1, nslip  )
        case( 9 ) !DJGM
          call mm10_b_mult_type_2a( dbar, np1%ms(1,1), vec1, nslip  )
        case default
          call mm10b_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_form_wbarp                   *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_form_wbarp( props, np1, n, vec1, vec2,
     &                            stress, tt, wbar )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(3) :: wbar
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
c
      integer :: i, nslip
      double precision :: slipinc, rs
      double precision, external :: mm10_rs
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, wbar:64

      nslip = props%nslip
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type )
        case( 1 )  ! voche
          wbar = zero
          do i = 1, nslip
            call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
            wbar = wbar + slipinc*np1%qs(1:3,i)
c                       addition for diffusion
            rs = mm10_rs( props, np1, n, stress, tt, i )
            wbar = wbar + rs*np1%tinc*props%iD_v*np1%qs(1:3,i)
          end do
        case( 2 ) ! MTS
          wbar = zero
          do i = 1, nslip
            call mm10_slipinc( props, np1, n, stress, tt, i, slipinc )
            wbar = wbar + slipinc*np1%qs(1:3,i)
          end do
        case( 3 ) ! User
          wbar = zero
          do i = 1, nslip
            call mm10_slipinc_user( props, np1, n, stress, tt, i, 
     &                              slipinc )
            wbar = wbar + slipinc*np1%qs(1:3,i)
          end do
        case( 4 ) ! ORNL
          call mm10_b_mult_type_2b( wbar, np1%qs(1,1), vec1, nslip  )          
        case( 7 )  ! MRR
          call mm10_b_mult_type_2b( wbar, np1%qs(1,1), vec1, nslip  )           
        case( 8 )  ! Armstrong-Frederick
          call mm10_b_mult_type_2b( wbar, np1%qs(1,1), vec1, nslip  )         
        case( 9 )  ! DJGM
          call mm10_b_mult_type_2b( wbar, np1%qs(1,1), vec1, nslip  )          
         case default
          call mm10b_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
      return
      end subroutine

c           Form w_p (in the current configuration)
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_form_wp                      *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_form_wp( props, np1, n, vec1, vec2, stress, 
     &                         tt, w )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(3) :: w
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
c
      integer :: i, nslip
      double precision :: slipinc, rs
      double precision, external :: mm10_rs
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, w:64
c
      nslip = props%nslip
c      
      select case ( props%h_type )
        case( 1 ) ! Voche
          w = zero
          do i = 1, nslip
            call mm10_slipinc( props, np1, n, stress, tt, i, slipinc )
            w = w + slipinc*np1%qc(1:3,i)
c             addition for diffusion
            rs = mm10_rs ( props, np1, n, stress, tt, i )
            w = w + rs*np1%tinc*props%iD_v*np1%qc(1:3,i)
          end do
        case( 2 ) ! MTS
          w = zero
          do i = 1, nslip
            call mm10_slipinc( props, np1, n, stress, tt, i, slipinc )
            w = w + slipinc*np1%qc(1:3,i)
          end do
        case( 3 ) ! User
          w = zero
          do i = 1, nslip
            call mm10_slipinc_user( props, np1, n, stress, tt, i, 
     &                              slipinc )
            w = w + slipinc*np1%qc(1:3,i)
          end do
        case( 4 ) ! ornl
          w = zero
          do i = 1, nslip
            call mm10_slipinc_ornl( props, np1, n, stress, tt, i,
     &                              slipinc)
            w = w + slipinc*np1%qc(1:3,i)
          end do
        case( 7 ) ! MRR
          w = zero
          do i = 1, nslip
            call mm10_slipinc_mrr( props, np1, n, stress, tt, i,
     &                             slipinc )
            w = w + slipinc*np1%qc(1:3,i)
          end do
        case( 8 ) ! Armstrong-Frederick
          w = zero
          do i = 1, nslip
            call mm10_slipinc_arfr( props, np1, n, stress, tt, i,
     &                             slipinc )
            w = w + slipinc*np1%qc(1:3,i)
          end do
        case( 9 ) ! DJGM
          w = zero
          do i = 1, nslip
            call mm10_slipinc_DJGM( props, np1, n, stress, tt, i, 
     &                              slipinc )
            w = w + slipinc*np1%qc(1:3,i)
          end do
        case default
          call mm10b_unknown_hard_error( props )
      end select
c
      return
      end 
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_form_dbarpi                  *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/18/2016 rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_form_dbarpi( props, np1, n, ivec1, ivec2,
     &                             stress, tt, dbar)
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      complex(kind=real64), dimension(6) :: stress
      complex(kind=real64), dimension(6) :: dbar
      complex(kind=real64), dimension(size_num_hard) :: tt
      complex(kind=real64), dimension(max_uhard) :: ivec1,ivec2
c
      integer :: i
      complex(kind=real64) :: slipinc
c!DIR$ ASSUME_ALIGNED ivec1:64, ivec2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, dbar:64
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type )
        case( 1, 2, 3, 8, 9 ) ! Voche, MTS, User
          continue
        case(4 ) ! ORNL
          dbar = (zero,zero)
          do i = 1, props%nslip
            slipinc = ivec1(i)
            dbar = dbar + slipinc*np1%ms(1:6,i)
          end do
        case( 7 ) ! MRR
          dbar = (zero,zero)
          do i = 1, props%nslip
            slipinc = ivec1(i)
            dbar = dbar + slipinc*np1%ms(1:6,i)
          end do
        case default
          call mm10b_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c      
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_form_wbarpi                  *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     *           Form w_p (in the current configuration)            *
c     *                                                              *
c     ****************************************************************
c
       subroutine mm10_form_wpi( props, np1, n, ivec1, ivec2, 
     &                           stress, tt, w)
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      complex(kind=real64), dimension(6) :: stress
      complex(kind=real64), dimension(3) :: w
      complex(kind=real64), dimension(size_num_hard) :: tt
      complex(kind=real64), dimension(max_uhard) :: ivec1,ivec2
c
      integer :: i
      complex(kind=real64) :: slipinc
c!DIR$ ASSUME_ALIGNED ivec1:64, ivec2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64, w:64
c
c ***** START: Add new Constitutive Models into this block *****
      select case( props%h_type )
         case( 1, 2, 3, 8, 9 )  !Voche, MTS, User
           continue
         case( 4 ) ! ORNL
           w = (zero,zero)
           do i = 1, props%nslip
            call mm10_slipinci_ornl( props, np1, n, stress, tt, i, 
     &                               slipinc)
            w = w + slipinc*np1%qc(1:3,i)
           end do
        case( 7 )  ! MRR
          w = (zero,zero)
          do i = 1, props%nslip
           call mm10_slipinci_mrr( props, np1, n, stress, tt, i, 
     &                             slipinc )
           w = w + slipinc*np1%qc(1:3,i)
          end do
        case default
          call mm10b_unknown_hard_error( props )
      end select
c ****** END: Add new Constitutive Models into this block ******
c
      return
      end subroutine

c
c **********************************************************************
c *                                                                    *
c *           mm10_symSW                                               * 
c *                                                                    *
c *         written by : mcm                                           *
c *         last modified : 10/15/2016 rhd                             *
c *                                                                    *
c *   Take the symmetric part of S*W for some stress and skew tensor   *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_symsw(s, w, sw)
      use mm10_constants
      implicit none
c      
      double precision, dimension(6), intent(in) :: s
      double precision, dimension(3), intent(in) :: w
      double precision, dimension(6), intent(out) :: sw
!dir$ assume_aligned s:64, w:64, sw:64     
c
      sw = zero
      sw(1) = s(4)*w(3) - s(6)*w(2)
      sw(2) = s(4)*w(3) - s(5)*w(1)
      sw(3) = s(6)*w(2) + s(5)*w(1)
      sw(4) = half*(w(3)*(s(1)-s(2)) + w(1)*s(6) - w(2)*s(5))
      sw(5) = half*(w(1)*(s(2)-s(3)) + w(2)*s(4) + w(3)*s(6))
      sw(6) = half*(w(2)*(s(1)-s(3)) + w(1)*s(4) - w(3)*s(5))
c
      return
c
      end
c
c **********************************************************************
c *                                                                    *
c *         written by : tjt                                           *
c *         last modified : 10/15/2016 rhd                             *
c *                                                                    *
c *   Take the symmetric part of S*W for some stress and skew tensor   *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_symswi(s, w, sw)
      use iso_Fortran_env
      use mm10_constants
      implicit none
c      
      complex(kind=real64), dimension(6), intent(in) :: s
      complex(kind=real64), dimension(3), intent(in) :: w
      complex(kind=real64), dimension(6), intent(out) :: sw
c!DIR$ ASSUME_ALIGNED S:64, W:64, SW:64     
c
      sw = zero
      sw(1) = s(4)*w(3) - s(6)*w(2)
      sw(2) = s(4)*w(3) - s(5)*w(1)
      sw(3) = s(6)*w(2) + s(5)*w(1)
      sw(4) = half*(w(3)*(s(1)-s(2)) + w(1)*s(6) - w(2)*s(5))
      sw(5) = half*(w(1)*(s(2)-s(3)) + w(2)*s(4) + w(3)*s(6))
      sw(6) = half*(w(2)*(s(1)-s(3)) + w(1)*s(4) - w(3)*s(5))
c
      return
c
      end 
c
c **********************************************************************
c *                                                                    *
c *                      mm10_symSWmat                                 *
c *                                                                    *
c *         written by : tjt                                           *
c *         last modified : 10/15/2016 rhd                             *
c *                                                                    *
c *   Take the symmetric part of S*W for some stress and skew tensor   *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_symswmat(s, w, n, sw)
      use mm10_constants
      implicit none
c      
      integer :: n
      double precision, dimension(6), intent(in) :: s
      double precision, dimension(3,n), intent(in) :: w
      double precision, dimension(6,n), intent(out) :: sw
c!DIR$ ASSUME_ALIGNED S:64, W:64, SW:64     
c
      sw = zero
      sw(1,1:n) = s(4)*w(3,1:n) - s(6)*w(2,1:n)
      sw(2,1:n) = s(4)*w(3,1:n) - s(5)*w(1,1:n)
      sw(3,1:n) = s(6)*w(2,1:n) + s(5)*w(1,1:n)
      sw(4,1:n) = half*(w(3,1:n)*(s(1)-s(2))
     &             + w(1,1:n)*s(6) - w(2,1:n)*s(5))
      sw(5,1:n) = half*(w(1,1:n)*(s(2)-s(3))
     &           + w(2,1:n)*s(4) + w(3,1:n)*s(6))
                  sw(6,1:n) = half*(w(2,1:n)*(s(1)-s(3))
     &             + w(1,1:n)*s(4) - w(3,1:n)*s(5))
c
      return
c
      end 
c
c
c **********************************************************************
c *                                                                    *
c *                             mm10_IW                                *
c *                                                                    *
c *         written by : mcm                                           *
c *         last modified : 10/15/2016 rhd                             *
c *                                                                    *
c *        Get Iik*Wlj - Wik*Ijl in our voigt notation                 *
c *                                                                    *
c **********************************************************************
c
      subroutine mm10_iw(w, iw)
      use mm10_constants
      implicit none
c            
      double precision, dimension(3), intent(in) :: w
      double precision, dimension(6,6), intent(out) :: iw
c!DIR$ ASSUME_ALIGNED w:64, iw:64    
c
      iw = zero
      iw(1,4) =  two*w(3)
      iw(1,6) = -two*w(2)
      iw(2,4) =  two*w(3)
      iw(2,5) = -two*w(1)
      iw(3,5) =  two*w(1)
      iw(3,6) =  two*w(2)
      iw(4,1) =  w(3)
      iw(4,2) = -w(3)
      iw(4,5) = -w(2)
      iw(4,6) =  w(1)
      iw(5,2) =  w(1)
      iw(5,3) = -w(1)
      iw(5,4) =  w(2)
      iw(5,6) =  w(3)
      iw(6,1) =  w(2)
      iw(6,3) = -w(2)
      iw(6,4) =  w(1)
      iw(6,5) = -w(3)
c
      return
      end
c
c **********************************************************************
c *                                                                    *
c *         User hardening routines                                    *
c *                                                                    *
c **********************************************************************
c
c
      subroutine mm10_slipinc_user( props, np1, n, stress, tt, i,
     &                              slipinc )
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision :: slipinc
      integer :: i
c
      write(props%out,*) "Not implemented: mm10_slipinc_user"
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
      double precision, dimension(size_num_hard) :: tt, h
c
      write(props%out,*) "Not implemented: mm10_h_user"
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
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(size_num_hard,6) :: et
c
      write(props%out,*) "Not implemented: mm10_estress_user"
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
      double precision, dimension(size_num_hard) :: tt
      double precision, 
     &   dimension(size_num_hard,size_num_hard) :: etau
c
      write(props%out,*) "Not implemented: mm10_ehard_user"
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
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(6,size_num_hard) :: ed
c
      write(props%out,*) "Not implemented: mm10_ed_user"
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
      double precision, dimension(size_nslip) :: dgammadtau
      double precision, dimension(size_num_hard) :: tt
c
c
      write(props%out,*) "Not implemented: mm10_dgdt_user"
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
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(size_nslip,size_num_hard)
     &        :: dgammadtt
c
c
      write(props%out,*) "Not implemented: mm10_dgdh_user"
      call die_gracefully
c
      return
      end subroutine
c
c           Derivative of sliprate wrt strain increment
      subroutine mm10_dgdd_user( props, np1, n, stress, tt, D,
     &                           dgammadd)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, D
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(size_nslip,6) :: dgammadd
c
c
      write(props%out,*) "Not implemented: mm10_dgdd_user"
      call die_gracefully
c
      return
      end subroutine
c
c **********************************************************************
c *                                                                    *
c *         Built in hardening routines                                *
c *                                                                    *
c **********************************************************************
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_slipinc                      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     *       Calculate the slip increment along system i            *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_slipinc( props, np1, n, stress, tt, i, slipinc )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt(*), slipinc
      integer :: i
c
      double precision :: rs
      double precision, external :: mm10_rs
c!DIR$ ASSUME_ALIGNED stress:64
c
      rs = mm10_rs( props, np1, n, stress, tt, i )
      slipinc = np1%dg/tt(1) * dabs(rs/tt(1))**(props%rate_n-one)*rs
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 function mm10_rs                             *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     &          Calculate the resolved shear along system i         *  
c     *                                                              *
c     ****************************************************************
c      
      function mm10_rs( props, np1, n, stress, tt, i )
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
c
      integer :: i      
      double precision :: stress(6), tt(*)
c
      double precision :: mm10_rs
c!DIR$ ASSUME_ALIGNED stress:64
c
      mm10_rs =  stress(1)*np1%ms(1,i)
     &         + stress(2)*np1%ms(2,i)
     &         + stress(3)*np1%ms(3,i)
     &         + stress(4)*np1%ms(4,i)
     &         + stress(5)*np1%ms(5,i)
     &         + stress(6)*np1%ms(6,i)
c
      return
      end
c      
      function mm10_rsi( props, np1, n, stress, tt, i )
      use iso_Fortran_env
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      integer :: i
      complex(kind=real64) :: stress(6), tt(*)
c
      complex(kind=real64) :: mm10_rsi
c!DIR$ ASSUME_ALIGNED stress:64
c
      mm10_rsi = stress(1)*np1%ms(1,i)
     &         + stress(2)*np1%ms(2,i)
     &         + stress(3)*np1%ms(3,i)
     &         + stress(4)*np1%ms(4,i)
     &         + stress(5)*np1%ms(5,i)
     &         + stress(6)*np1%ms(6,i)
c
      return
      end 
c -----------------------
c     Simple voche:
c
c -----------------------
c
c           Actual voche law hardening function
      subroutine mm10_h_voche(props, np1, n, stress, tt, h)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(1) :: tt, h
      double precision :: h_term
      integer :: i
c
      double precision :: slipinc
c
      h = zero
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
        h_term = one - (tt(1)-props%tau_y)/props%tau_v
     &           + np1%tau_l(i)/(tt(1)-props%tau_y)
        h(1) = h(1) + abs(h_term)**(props%voche_m) * 
     &      sign(one,h_term)
     &      * abs(slipinc)
      end do
      h(1) = n%tau_tilde(1) + props%theta_0*h(1)
c
      return
      end subroutine
c
c           Derivative of hardening fn wrt stress
      subroutine mm10_estress_voche(props, np1, n, stress, tt, et)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(1) :: tt
      double precision, dimension(6) :: et
c
      double precision :: mm10_rs
      double precision :: rs, h_term
      integer :: i
c
      et = zero
c
      do i=1,props%nslip
        rs = mm10_rs(props, np1, n, stress, tt, i)
        h_term = one - (tt(1)-props%tau_y)/props%tau_v
     &           + np1%tau_l(i)/(tt(1)-props%tau_y)
        et(1:6) = et(1:6) + abs(h_term)**(props%voche_m) * 
     &      sign(one,h_term) * 
     &      dabs(rs)**(props%rate_n-two)*rs*np1%ms(1:6,i)
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
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(*) :: tt
      double precision :: h_term
      double precision, 
     &      dimension(size_num_hard,size_num_hard) :: etau
c
      double precision :: slipinc
      integer :: i
c
      etau(1,1) = zero
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
        h_term = one - (tt(1)-props%tau_y)/props%tau_v
     &           + np1%tau_l(i)/(tt(1)-props%tau_y)
        etau(1,1) = etau(1,1) + ( props%voche_m * 
     &      (-one/props%tau_v - np1%tau_l(i)/
     &      (tt(1)-props%tau_y)**2) * abs(slipinc)/abs(h_term)
     &      - abs(slipinc) * props%rate_n/tt(1) * 
     &      sign(one,h_term)
     &      * sign(one,slipinc) ) * (abs(h_term)**props%voche_m)
      end do

      etau(1,1) = props%theta_0 * etau(1,1)
      etau(1,1) = one - etau(1,1)
c
      return
      end subroutine
c
c           Derivative of hardening fn wrt strain
      subroutine mm10_ed_voche(props, np1, n, stress, tt, ed)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(*) :: tt
      double precision, dimension(6) :: ed
c
      ed = zero
c
      return
      end subroutine
c
c           Derivative of sliprate wrt resolved shear stress
      subroutine mm10_dgdt_voche(props, np1, n, stress, tt, dgammadtau)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt(*), rs
      double precision, dimension(size_nslip) :: dgammadtau
c
      double precision, external :: mm10_rs
c
      integer :: s
c
      do s=1,props%nslip
        rs = mm10_rs(props, np1, n, stress, tt, s)
        dgammadtau(s) = dabs(rs)**(props%rate_n-one)
        dgammadtau(s) = np1%dg*props%rate_n/tt(1)**(props%rate_n)
     &     *dgammadtau(s)
c   additional term for diffusion
        dgammadtau(s) = dgammadtau(s) + np1%tinc*props%iD_v
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
      double precision :: tt(*), dgam
      double precision, dimension(size_nslip,1) :: dgammadtt
      integer :: s
c
      do s=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, s, dgam)
        dgammadtt(s,1) = -props%rate_n/tt(1) * dgam
      end do
c
      return
      end subroutine
c
c           Derivative of sliprate wrt strain increment
      subroutine mm10_dgdd_voche(props, np1, n, stress, tt, D, 
     &                         dgammadd)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, D
      double precision :: tt(*)
      double precision, dimension(size_nslip,6) :: dgammadd
c
      dgammadd = zero
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
      use mm10_constants
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
      ct = one - cta/np1%tau_v
      h = zero
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
        h(1) = h(1) + (ct + np1%tau_l(i)/cta)**(props%voche_m)*
     &      dabs(slipinc)
      end do
c
      h(1) = props%tau_a*(one - np1%mu_harden/n%mu_harden) + 
     &      (np1%mu_harden/props%mu_0)*(np1%tau_y - n%tau_y) + 
     &      (np1%mu_harden/n%mu_harden)*n%tau_tilde(1) + 
     &      props%theta_0 * (np1%mu_harden/props%mu_0)*h(1)

      return
      end subroutine
c
c           Derivative of hardening fn wrt stress
      subroutine mm10_estress_mts(props, np1, n, stress, tt, et)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt(*)
      double precision, dimension(6) :: et
c
      double precision, external :: mm10_rs
c
      double precision :: rs, cta, ct
      integer :: i
c
      cta = (props%mu_0/
     &      np1%mu_harden)*tt(1) -
     &      (props%mu_0/np1%mu_harden)*props%tau_a -
     &      np1%tau_y
      ct = one - cta/np1%tau_v
      et = zero
      do i = 1, props%nslip
        rs = mm10_rs(props, np1, n, stress, tt, i)
        et = et + (ct+np1%tau_l(i)/cta)**(props%voche_m)*
     &      dabs(rs)**(props%rate_n-two)*rs*np1%ms(1:6,i)
      end do

      et =  props%theta_0 * 
     &      (np1%mu_harden/props%mu_0)*et*props%rate_n*
     &      np1%dg/tt(1)**props%rate_n

c
      return
      end subroutine
c
c           Derivative of hardening fn wrt hardening
      subroutine mm10_ehard_mts(props, np1, n, stress, tt, etau)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt(*)
      double precision, 
     &      dimension(size_num_hard,size_num_hard) :: etau
c
      integer :: i
      double precision :: slipinc, ur, ct, cta
c
      cta = (props%mu_0/
     &      np1%mu_harden)*tt(1) -
     &      (props%mu_0/np1%mu_harden)*props%tau_a -
     &      np1%tau_y
      ct = one - cta/np1%tau_v
c
      ur = np1%mu_harden / props%mu_0
c
      etau(1,1) = zero
      do i=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, i, slipinc)
        etau(1,1) = etau(1,1) + (props%voche_m*
     &    (one/np1%tau_v+np1%tau_l(i)/ cta**2)*
     &      (ct+np1%tau_l(i)/cta)**(-one) + ur*props%rate_n/tt(1))*
     &      (ct+np1%tau_l(i)/cta)**(props%voche_m)*
     &      dabs(slipinc)
      end do
c
      etau(1,1) = -props%theta_0 * etau(1,1)
      etau(1,1) = one - etau(1,1)

      return
      end subroutine
c
c           Derivative of hardening fn wrt strain
      subroutine mm10_ed_mts(props, np1, n, stress, tt, ed)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt(*)
      double precision, dimension(6) :: ed
c
      double precision :: slipinc
c
      double precision :: lnv, lny, dgc, ty, tv, mnp0, sc
      double precision, dimension(6) :: d_mod, dydd, dvdd
      integer :: s
c
c     Form a bunch of simple constitutive things
c
c
      d_mod = np1%D
      d_mod(4:6) = half * d_mod(4:6)
c
      dgc = np1%dg / np1%tinc
c
c     Form d_tauy/d_deltad by steps
c
      lny = dlog(props%eps_dot_0_y/dgc)
      ty = props%boltzman*np1%temp
     &       /(np1%mu_harden*(props%burgers**3)*props%G_0_y)*
     &       lny
      dydd = two*props%tau_hat_y/(three*np1%dg**2*
     &        props%q_y*props%p_y*
     &      lny)*(one-ty**(one/props%q_y))**(one/props%p_y-one)*
     &      ty**(one/props%q_y)*d_mod

c
c     Form d_tauv/d_deltad by steps
c
      lnv = dlog(props%eps_dot_0_v/dgc)
      tv = props%boltzman*np1%temp
     &       /(np1%mu_harden*(props%burgers**3)*props%G_0_v)*
     &       lnv
      dvdd = two*props%tau_hat_v/(three*np1%dg**2*
     &        props%q_v*props%p_v*
     &      lnv)*(one-tv**(one/props%q_v))**(one/props%p_v-one)*
     &      tv**(one/props%q_v)*d_mod

c
c     Form a couple more common components
c
      mnp0 = np1%mu_harden / props%mu_0
      sc = tt(1)/mnp0 - props%tau_a/mnp0 - np1%tau_y
c
c     Glue everything together
c
      ed = zero
      do s=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, s, slipinc)
        ed = ed + (props%voche_m*(one/np1%tau_v
     &      +np1%tau_l(s)/sc**2)*
     &      (one-sc/np1%tau_v+np1%tau_l(s)/sc)
     &         **(props%voche_m-one)*dydd
     &      + props%voche_m/(np1%tau_v)**2*sc*
     &      (one-sc/np1%tau_v+np1%tau_l(s)/sc)
     &         **(props%voche_m-one)*dvdd
     &      + two/(three*np1%dg**2)*(one-sc/np1%tau_v
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
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: tt(*), rs
      double precision, dimension(size_nslip) :: dgammadtau
c
      double precision, external :: mm10_rs
c
      integer :: s
c
      do s=1,props%nslip
        rs = mm10_rs(props, np1, n, stress, tt, s)
        dgammadtau(s) = dabs(rs)**(props%rate_n-one)
        dgammadtau(s) = np1%dg*props%rate_n/tt(1)**(props%rate_n)
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
      double precision :: tt(*), dgam
      double precision, dimension(size_nslip,1) :: dgammadtt
      integer :: s
c
      do s=1,props%nslip
        call mm10_slipinc(props, np1, n, stress, tt, s, dgam)
        dgammadtt(s,1) = -props%rate_n/tt(1) * dgam
      end do
c
      return
      end subroutine
c
c           Derivative of sliprate wrt strain increment
      subroutine mm10_dgdd_mts(props, np1, n, stress, tt, D, 
     &                         dgammadd)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, D, d_mod
      double precision :: tt(*), alpha, dgam
      double precision, dimension(6,size_nslip) :: dgammadd
      integer :: s
c
      d_mod = D
      d_mod(4:6) = half * d_mod(4:6)
      alpha = two/(three*np1%dg**2)
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
      double precision, dimension(size_num_hard) :: tt
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

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_a_mrr                        *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     *     Form the stress varying with hardening part              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_a_mrr( props, np1, n, stress, tt, vec1, vec2,
     &                       arr1, arr2, both )
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      integer both      
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      integer :: nslip, num_hard  
      double precision, allocatable :: temp_arr2(:,:)    
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, arr1:64, arr2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64     
c
c              arr1 treated as vector inside
c
      nslip    = props%nslip
      num_hard = props%num_hard
c      
      call mm10_dgdt_mrr( props, np1, n, stress, tt, arr1(1,1) )
c      
c              arr2 must have lead dimension = nslip for compatibility 
c              with _dgdh_mrr. values defined inside the routine - no
c              need to copyin. only copyout needed
c
      if( both == 2 ) then
         allocate( temp_arr2(nslip,num_hard) )
         call mm10_dgdh_mrr( props, np1, n, stress, tt, temp_arr2 )
         arr2(1:nslip,1:num_hard) = temp_arr2
         deallocate( temp_arr2 ) 
      end if
c
      return
      end 
c
c           Form intermediate vectors for faster calculations
      subroutine mm10_vi_mrr(props, np1, n, stress, tt, 
     &   ivec1, ivec2)
      use iso_Fortran_env
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      complex(kind=real64), dimension(6) :: stress
      complex(kind=real64), dimension(size_num_hard) :: tt
      complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
      integer :: alpha
      complex(kind=real64) :: slipinc
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
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision :: slipinc
      double precision, external :: mm10_rs
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, rs,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m
       !double precision, dimension(size_num_hard,size_num_hard)
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
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
c        
c      ms = np1.ms(1:6,alpha);
c      rs = stress*ms; % tau^a
      rs = mm10_rs( props, np1, n, stress, tt, alpha )
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
        if(fract.gt.one) then
            ! linear extrapolation past the too-high stress (rs) value
            b = gamma_0
            x = fract
            m = b * (-q_e*(Qslip/k/theta))
     &         * (- p_e) * dsign(one,rs)/tcut
            y = m*x + b
            slipinc = dt * y
        elseif(dabs(rs).gt.zero) then
        slipinc = dt * gamma_0 * dsign(one,rs)
     & * dexp (-(Qslip/k/theta)*(one - fract**p_e)**q_e) !(14)
        else
            slipinc = zero
        endif
c
      return
      end 
c
c           Actual mrr hardening function
      subroutine mm10_h_mrr(props, np1, n, vec1, vec2, 
     & stress, tt, h,gp)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt, h,
     &       rhoFs,rhoPs
      double precision, dimension(max_uhard) :: vec1, vec2
      integer :: alpha, gp
      logical :: mat_debug
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  rs,
     &  rhoF, rhoP, tpass, rho,
     &  rho_n, c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  tem1, tem2, tem3, tem4
       !double precision, dimension(size_num_hard,size_num_hard)
       !&   :: Gmat, Hmat

       mat_debug = .false.

c Load some material parameters
        Qbulk = props%tau_a
        k = props%boltzman
        theta = np1%temp
        if(theta.le.zero) then
          write (props%out,*) "Roters model requires non-zero 
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
c        
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Gmat,props%num_hard,tt,1,zero,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Hmat,props%num_hard,tt,1,zero,rhoPs,1)
c
c
      do alpha = 1,props%num_hard
c
c Get dislocation density
          rho = tt(alpha) ! rho^a_SSD
          rho_n = n%tau_tilde(alpha) ! rho^a_SSD
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
        if(mat_debug) then
        write(*,*) "rhoF", rhoF
        write(*,*) "rhoP", rhoP
        endif
c          
          tpass = c1*G*b*dsqrt(rhoP) ! (16)
          ddipole = dsqrt(three)*G*b/(16.d0*pi*(one-v))/
     &       (dabs(rs)) ! (42)
          rhoM = (two*k/(c1*c2*c3*G*b**3))*
     &       theta*dsqrt(rhoF*rhoP) ! (13)
          slipinc = vec1(alpha)
          gammadot = dabs(slipinc/dt)
      tem1 = c4/b*dsqrt(rhoP)*gammadot
      tem2 = c6*ddipole/b*rhoM*gammadot
      tem3 = c5*rho*gammadot
      tem4 = c7*dexp(-Qbulk/k/theta)*dabs(rs)/(k*theta)
     &        *rho*rho*gammadot**c8
          h(alpha) = rho_n + dt*(tem1
     &        + tem2 - tem3
     &        - tem4) ! (18)
      enddo
c
      return
      end 
c
c           Imaginary mrr sliprate function
      subroutine mm10_slipinci_mrr(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      complex(kind=real64), dimension(6) :: stress
      complex(kind=real64), dimension(size_num_hard) :: tt, temp
      complex(kind=real64) :: slipinc
      complex(kind=real64), external :: mm10_rsi
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack
      complex(kind=real64) :: rs, bb,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, x, y, m
       !double precision, dimension(size_num_hard,size_num_hard)
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
c                Compute some stresses and rates
        gamma_0 = v_attack*k*theta/(c1*c3*G*b*b)*cdsqrt(rhoP) ! (15)
        tpass = c1*G*b*cdsqrt(rhoP) ! (16)
        tcut = Qslip/(c2*c3*b*b)*cdsqrt(rhoF) ! (17)
          if(dreal(rs).lt.zero) then
        fract = (-rs-tpass)/tcut
          else
        fract = (rs-tpass)/tcut
          endif
c
c            Evaluate the slip rate equation
        if(dreal(fract).gt.one) then
            ! linear extrapolation past the too-high stress (rs) value
            bb = gamma_0
            x = fract
            m = bb * (-q_e*(Qslip/k/theta))
     &         * (- p_e) * dsign(one,dreal(rs))/tcut
            y = m*x + bb
            slipinc = dt * y
        elseif(dabs(dreal(rs)).gt.zero) then
        slipinc = dt * gamma_0 * dsign(one,dreal(rs))
     & * cdexp (-(Qslip/k/theta)*(one - fract**p_e)**q_e) !(14)
        else
            slipinc = zero
        endif
c
      return
      end 
c
! c           Imaginary mrr hardening function
c
      subroutine mm10_hi_mrr( props, np1, n, ivec1, ivec2, 
     &                        stress, tt, h)
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      complex(kind=real64), dimension(6) :: stress
      complex(kind=real64), dimension(size_num_hard) :: tt, h, temp
      complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
      integer :: alpha
c
      double precision :: 
     &  dt, k_boltz, theta, G, b, c1, c2, c3, rho_n, c4, c5, c6, 
     &  c7, c8, v, Qbulk, tem0, tem4, chk_size, tem5
      double precision, parameter :: zero_tol = 1.0d-30
c     
      complex(kind=real64) ::
     &  rs, rhoF, rhoP, tpass, rho, mm10_rsi, ddipole, rhoM, slipinc, 
     &  gammadot, tem1, tem2, tem3
c
      logical :: zero_flag
c     
       !double precision, dimension(size_num_hard,size_num_hard)
       !&   :: Gmat, Hmat
c        
c                        Load some material parameters
c
      Qbulk   = props%tau_a
      k_boltz = props%boltzman
      theta   = np1%temp
c      
      if( theta .le. zero ) then
          write (props%out,*) "Roters model requires non-zero ",
     &                        "nodal temperatures"
          call die_gracefully
      end if
c
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
c        
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
c
c
      do alpha = 1,props%num_hard
c
c Get dislocation density
          rho = tt(alpha) ! rho^a_SSD
          rho_n = n%tau_tilde(alpha) ! rho^a_SSD
c
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rsi(props, np1, n, stress, tt, alpha)
          rs = dsign(one,dreal(rs))*rs
c          
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          temp = (props%Gmat(alpha,1:props%num_hard)
     &      *tt(1:props%num_hard))
          rhoF = sum(temp)
          temp = (props%Hmat(alpha,1:props%num_hard)
     &      *tt(1:props%num_hard))
          rhoP = sum(temp)
c          
          tpass = c1*G*b*cdsqrt(rhoP) ! (16)
          ddipole = root3*G*b/(16.d0*pi*(one-v))/
     &       (rs) ! (42)
          rhoM = (two*k_boltz/(c1*c2*c3*G*b**3))*
     &       theta*cdsqrt(rhoF*rhoP) ! (13)
c          
c          call mm10_slipinc_mrr(props, np1, n, stress, tt, alpha, 
c     &     slipinc)
          slipinc = ivec1(alpha)
          slipinc = dsign(one,dreal(slipinc)) * slipinc
          gammadot = slipinc/dt
          zero_flag = .false.
          chk_size = dsqrt( dble(gammadot)**2 + dimag(gammadot)**2 )
          if( chk_size .lt. zero_tol ) then
            gammadot = dcmplx( zero, zero )
            zero_flag = .true.
          end if  
c              
c             Evaluate the hardening equation
c
          tem0 = c4/b*cdsqrt(rhoF)*gammadot
          tem2 = c6*ddipole/b*rhoM*gammadot
          tem3 = c5*rho*gammadot
          tem4 = dexp(-Qbulk/k_boltz/theta)
          tem5 = dcmplx( zero, zero )
          if( .not. zero_flag ) tem5 = gammadot**c8
          tem1 = c7*tem4*rs/(k_boltz*theta)*rho*rho*tem5
          h(alpha) = rho_n + dt*( tem0 + tem2 - tem3 - tem1 ) 
      end do
c
      return
      end
c
      
c
c           Wrapper version, mrr slipinc function
      subroutine mm10_slipinc_mrrW(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt, zerosV
      double precision :: slipinc
      integer :: alpha
c
      complex(kind=real64), dimension(6) :: stressi
      complex(kind=real64), dimension(size_num_hard) :: tti
      complex(kind=real64) :: slipinci
c
      zerosV = zero
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
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, dtdstress
      double precision, dimension(size_num_hard) :: tt,rhoPs,rhoFs
      double precision, dimension(size_num_hard,6) :: et
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  rs,
     &  rhoF, rhoP, tpass, rho,
     &  c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  dddipole, dslipinc, badterm
       !double precision, dimension(size_num_hard,size_num_hard)
       !&   :: Gmat, Hmat
      double precision, dimension(size_nslip) :: dslip
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
c
c compute derivatives of slip increments with respect to resolved
c shear stress
c        call mm10_dgdt_mrr(props, np1, n, stress, 
c     &         tt, dslip)
        dslip(1:props%num_hard) = arr1(1:props%num_hard,1)
c        
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Gmat,props%num_hard,tt,1,zero,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Hmat,props%num_hard,tt,1,zero,rhoPs,1)
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
          ddipole = root3*G*b/(16.d0*pi*(one-v))/
     &       (dabs(rs)) ! (42)
          dddipole = -root3*G*b/(16.d0*pi*(one-v))/
     &          (dabs(rs))**2*dsign(one,rs)
          rhoM = (two*k/(c1*c2*c3*G*b**3))*
     &       theta*dsqrt(rhoF*rhoP) ! (13)
c          
c          call mm10_slipinc_mrr(props, np1, n, stress, tt, alpha, 
c     &     slipinc)
          slipinc = vec1(alpha)
          gammadot = dabs(slipinc/dt)
          dtdstress(1:6) = np1%ms(1:6,alpha)
          dslipinc = dslip(alpha)/dt ! always positive

          ! Evaluate the equation
          if(gammadot.eq.zero) then
          badterm = zero
          else
          badterm = c8*dabs(rs)*gammadot**(c8-one)*dslipinc
          endif
          et(alpha,1:6) = dt*(c4/b*dsqrt(rhoF)*dsign(one,rs)
     &         *dslipinc + c6/b*rhoM*(ddipole*dsign(one,rs)
     &         *dslipinc + dddipole*gammadot)
     &        - c5*rho*dsign(one,rs)*dslipinc
     &    - c7*dexp(-Qbulk/k/theta)/(k*theta)*rho**2*
     &      (badterm + dsign(one,rs)*gammadot**c8))
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
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt,rhoPs,rhoFs
      double precision, dimension(size_num_hard,size_num_hard) :: etau
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  rs,
     &  rhoF, rhoP, tpass, rho,
     &  c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  dddipole, dslipinc, badterm, deltaij,
     &  drhoF, drhoP, drhoM
       !double precision, dimension(size_num_hard,size_num_hard)
       !&   :: Gmat, Hmat
      double precision, dimension(size_num_hard,size_num_hard)
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
c       
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Gmat,props%num_hard,tt,1,zero,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Hmat,props%num_hard,tt,1,zero,rhoPs,1)
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
        ddipole = root3*G*b/(16.d0*pi*(one-v))/
     &       (dabs(rs)) ! (42)
        rhoM = (two*k/(c1*c2*c3*G*b**3))*
     &       theta*dsqrt(rhoF*rhoP) ! (13)
c          
c loop over denominator hardening variable
        do beta = 1,props%num_hard
c          
c           [drhoF,drhoP] = mm10_drhoFP_mrr(props, np1, n, tt, alpha, beta);
          drhoF = props%Gmat(alpha,beta)
          drhoP = props%Hmat(alpha,beta)
          dddipole = zero
          drhoM = half*(two*k/(c1*c2*c3*G*b**3))*
     &       theta/dsqrt(rhoF*rhoP) * (drhoF*rhoP + rhoF*drhoP)
          dslipinc = dslip(alpha,beta)/dt
          if(alpha.eq.beta) then
              deltaij = one
          else
              deltaij = zero
          endif

c Evaluate the equation
          if(gammadot.eq.zero) then
          badterm = zero
          else
          badterm = c8*gammadot**(c8-one)*dslipinc
          endif
          etau(alpha,beta) = deltaij - dt*(c4/b*(half*drhoF/
     &        dsqrt(rhoF)*gammadot*dsign(one,rs) 
     &     + dsqrt(rhoF)*dslipinc) + c6/b*(rhoM*ddipole*dslipinc
     &     + rhoM*dddipole*gammadot*dsign(one,rs)
     &     + drhoM*ddipole*gammadot*dsign(one,rs))
     &     - c5*(rho*dslipinc + deltaij*gammadot*dsign(one,rs))
     &     - c7*dexp(-Qbulk/k/theta)/(k*theta)*dabs(rs)*
     &      (rho**2*badterm + two*rho*dsign(one,rs)
     &      *deltaij*gammadot**c8))*dsign(one,rs)
      
        end do !beta
      
      end do !alpha
c
      return
      end 
c
c
c           Derivative of hardening fn wrt strain
      subroutine mm10_ed_mrr(props, np1, n, stress, tt, ed)   
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(6,size_num_hard) :: ed
c
      ed = zero
c
      return
      end subroutine
c
c           Derivative of sliprate wrt resolved shear stress
      subroutine mm10_dgdt_mrr(props, np1, n, stress, tt, dgammadtau)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: rs
      double precision, dimension(size_nslip) :: dgammadtau
      double precision, dimension(size_num_hard) :: tt,rhoFs,rhoPs
      double precision, external ::  mm10_rs
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, dfract,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, m,
     &  dslipinc, slipexp
       !double precision, dimension(size_num_hard,size_num_hard)
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
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Gmat,props%num_hard,tt,1,zero,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Hmat,props%num_hard,tt,1,zero,rhoPs,1)
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
            m = b * (-q_e*(Qslip/k/theta))
     &         * (- p_e) * dsign(one,rs)/tcut
            dslipinc = dt * m
        elseif(dabs(rs).gt.zero) then
        slipexp = dexp (-(Qslip/k/theta)*(one - fract**p_e)**q_e)
     &          * dsign(one,rs)
        dslipinc = dt * (-gamma_0 * slipexp *(Qslip/k/theta)*q_e*
     &      (one - fract**p_e)**(q_e-one) * 
     &      (-one)*p_e*fract**(p_e-one) * dfract) !(14)
        else
            dslipinc = zero
        endif
        
        dgammadtau(alpha) = dslipinc

      end do
c
      return
      end 
c
c           Derivative of sliprate wrt hardening variables
      subroutine mm10_dgdh_mrr(props, np1, n, stress, tt, dgammadtt)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt,rhoFs,rhoPs
      double precision, dimension(size_nslip,size_num_hard)
     &    :: dgammadtt
      double precision, external :: mm10_rs
      double precision :: rs
      integer :: alpha, beta
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  p_e, q_e, Qslip, v_attack, dfract,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract,
     &  dslipinc, slipexp, drhoF, drhoP, dgamma_0,
     &  dtcut, dtpass
       !double precision, dimension(size_num_hard,size_num_hard)
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
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Gmat,props%num_hard,tt,1,zero,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Hmat,props%num_hard,tt,1,zero,rhoPs,1)
        
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
        
c       [drhoF,drhoP] = mm10_drhoFP_mrr(props, np1, n, tt, alpha, beta);
          drhoF = props%Gmat(alpha,beta)
          drhoP = props%Hmat(alpha,beta)
        
          dgamma_0 = half*v_attack*k*theta/(c1*c3*G*b*b)
     &              /dsqrt(rhoP)*drhoP ! (15)
          dtpass = half*c1*G*b/dsqrt(rhoP)*drhoP ! (16)
          dtcut = half*Qslip/(c2*c3*b*b)/dsqrt(rhoF)*drhoF ! (17)
          dfract = ((-dtpass)*tcut - 
     &             (dabs(rs)-tpass)*dtcut)/tcut**2

c Evaluate the equation
          if(fract.gt.one) then
              ! linear extrapolation past the too-high stress (rs) value
              dslipinc = dt * (-q_e*(Qslip/k/theta)) * (- p_e)
     &            * dsign(one,rs) * (dgamma_0/tcut 
     &             - gamma_0*dtcut/tcut**2)
          elseif(dabs(rs).gt.zero) then
          slipexp = dexp (-(Qslip/k/theta)*
     &             (one - fract**p_e)**q_e) !(14)
          dslipinc = dt * (dgamma_0 * slipexp * dsign(one,rs)
     &          + gamma_0 * dsign(one,rs) * slipexp *(-one)*
     &          (Qslip/k/theta)*q_e*(one - fract**p_e)**(q_e-one)
     &        * (-one)*p_e*fract**(p_e-one) * dfract)
          else
              dslipinc = zero
          endif
        
          dgammadtt(alpha,beta) = dslipinc

        enddo !beta
      
      enddo !alpha
c
      return
      end
c
c           Derivative of sliprate wrt strain increment
      subroutine mm10_dgdd_mrr(props, np1, n, stress, tt, D, 
     &              dgammadd)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, D
      double precision, dimension(6,size_nslip) :: dgammadd
      double precision, dimension(size_num_hard) :: tt
c
      dgammadd = zero
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
      double precision, dimension(size_num_hard) :: tt
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
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_a_ornl                       *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_a_ornl( props, np1, n, stress, tt, vec1, vec2,
     &                        arr1, arr2, both )
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      integer :: both
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      integer :: nslip, num_hard  
      double precision, allocatable :: temp_arr2(:,:)    
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, arr1:64, arr2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64     
c
c              arr1 treated as vector inside
c
      nslip    = props%nslip
      num_hard = props%num_hard
c      
      call mm10_dgdt_ornl( props, np1, n, stress, tt, arr1(1,1) )
c      
c              arr2 must have lead dimension = nslip for compatibility 
c              with _dgdh_ornl. values defined inside the routine - no
c              need to copyin. only copyout needed
c
      if( both == 2 ) then
         allocate( temp_arr2(nslip,num_hard) )
         call mm10_dgdh_ornl( props, np1, n, stress, tt, temp_arr2 )
         arr2(1:nslip,1:num_hard) = temp_arr2
         deallocate( temp_arr2 ) 
      end if
c
      return
      end 
c
c           Form intermediate vectors for faster calculations
      subroutine mm10_vi_ornl(props, np1, n, stress, tt, 
     &   ivec1, ivec2)
      use iso_Fortran_env
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      complex(kind=real64), dimension(6) :: stress
      complex(kind=real64), dimension(size_num_hard) :: tt
      complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
      integer :: alpha
      complex(kind=real64) :: slipinc
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
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision :: slipinc, mm10_rs
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, tau0, 
     &  p_e, q_e, Qslip, v_s, rs,
     &  rhoP, gamma_0, tpass, tcut, fract, x, y, m,
     &  fM, lamda, G0, rhoM
       !double precision, dimension(size_num_hard,size_num_hard)
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
        if( v_s.lt.zero ) then
            v_s = exp(-v_s)
        endif
        fM = 0.1d0
        lamda = c2*b
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - one)
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

        if(fract.gt.one) then
            ! linear extrapolation past the too-high stress (rs) value
            b = gamma_0
            x = fract
            m = b * (-q_e*(Qslip/k/theta))
     &         * (- p_e) * dsign(one,rs)/tcut
            y = m*x + b
            slipinc = dt * y
        elseif(dabs(rs).eq.zero) then
            slipinc = zero
        else
          if(fract.lt.zero) then
            p_e = one ! deal with low stresses by allowing the fraction to be small
          endif
          slipinc = dt * gamma_0 * dsign(one,rs)
     & * dexp (-(Qslip/k/theta)*(one - fract**p_e)**q_e) !(14)x
        endif
c
      return
      end 
c
c           Actual ornl hardening function
      subroutine mm10_h_ornl(props, np1, n, vec1, vec2, 
     & stress, tt, h,gp)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt, h,
     &   rhoFs,rhoPs
      double precision, dimension(max_uhard) :: vec1, vec2
      integer :: alpha, gp
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  rs,
     &  rhoF, rhoP, tpass, rho,
     &  rho_n, c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  tem1, tem2, tem3, tem4, G0
       !double precision, dimension(size_num_hard,size_num_hard)
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
c        
c              New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - one)
c       Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Gmat,props%num_hard,tt,1,zero,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Hmat,props%num_hard,tt,1,zero,rhoPs,1)
c
c
      do alpha = 1,props%num_hard
c
c Get dislocation density
          rho = tt(alpha) ! rho^a_SSD
          rho_n = n%tau_tilde(alpha) ! rho^a_SSD
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
          ddipole = root3*G*b/(16.d0*pi*(one-v))/
     &       (dabs(rs)) ! (42)
          rhoM = (two*k/(c1*c2*c3*G*b**3))*
     &       theta*dsqrt(rhoF*rhoP) ! (13)
c          
c          call mm10_slipinc_ornl(props, np1, n, stress, tt, alpha, 
c     &     slipinc)
          slipinc = vec1(alpha)
          gammadot = dabs(slipinc/dt)
      tem1 = c4/b*dsqrt(rhoP)*gammadot
      tem2 = c6*ddipole/b*rhoM*gammadot
      tem3 = c5*rho*gammadot
      tem4 = c7*dexp(-Qbulk/k/theta)*dabs(rs)/(k*theta)
     &        *rho*rho*gammadot**c8
          h(alpha) = rho_n + dt*(tem1
     &        + tem2 - tem3
     &        - tem4) ! (18)
      enddo
c
      return
      end 
c
c           Imaginary ornl sliprate function
      subroutine mm10_slipinci_ornl(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      complex(kind=real64), dimension(6) :: stress
      complex(kind=real64), dimension(size_num_hard) :: tt, temp
      complex(kind=real64) :: slipinc
      complex(kind=real64), external :: mm10_rsi
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, tau0, 
     &  p_e, q_e, Qslip, v_s, G0
      complex(kind=real64) :: rs,
     &  rhoP, gamma_0, tpass, tcut, fract, x, y, m, bb,
     &  fM, lamda, rhoM
       !double precision, dimension(size_num_hard,size_num_hard)
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
        if( v_s.lt.zero ) then
            v_s = exp(-v_s)
        endif
        fM = 0.1d0
        lamda = c2*b
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - one)
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
          if(dreal(rs).lt.zero) then
        fract = (-rs-tpass)/tcut
          else
        fract = (rs-tpass)/tcut
          endif
c
c Evaluate the slip rate equation
        if(dreal(fract).gt.one) then
            ! linear extrapolation past the too-high stress (rs) value
            bb = gamma_0
            x = fract
            m = bb * (-q_e*(Qslip/k/theta))* (- p_e) *
     &           dsign(one,dreal(rs))/tcut
            y = m*x + bb
            slipinc = dt * y
        elseif(dabs(dreal(rs)).eq.zero) then
            slipinc = zero
        else
          if(dreal(fract).lt.zero) then
            p_e = one ! deal with low stresses by allowing the fraction to be small
          endif
          slipinc = dt * gamma_0 * dsign(one,dreal(rs))
     & * cdexp (-(Qslip/k/theta)*(one - fract**p_e)**q_e) !(14)
        endif
c
      return
      end 
c
c           Imaginary ornl hardening function
      subroutine mm10_hi_ornl(props, np1, n, ivec1, ivec2, 
     & stress, tt, h)
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      complex(kind=real64), dimension(6) :: stress
      complex(kind=real64), dimension(size_num_hard) :: tt, h, temp
      complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  rho_n, G0,
     &  c4, c5, c6, c7, c8, v, Qbulk
      complex(kind=real64) :: rs,
     &  rhoF, rhoP, tpass, rho,
     &  mm10_rsi,
     &  ddipole, rhoM, slipinc, gammadot, 
     &  tem1, tem2, tem3
       !double precision, dimension(size_num_hard,size_num_hard)
       !&   :: Gmat, Hmat
c Load some material parameters
        Qbulk = props%tau_a
        k = props%boltzman
        theta = np1%temp
        if(theta.le.zero) then
            write (props%out,9000)
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
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - one)
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
c
c
      do alpha = 1,props%num_hard
c
c                Get dislocation density
          rho = tt(alpha) ! rho^a_SSD
          rho_n = n%tau_tilde(alpha) ! rho^a_SSD
c
c          ms = np1.ms(1:6,alpha);
c          rs = stress*ms; % tau^a
          rs = mm10_rsi(props, np1, n, stress, tt, alpha)
          rs = dsign(1.d0,dreal(rs))*rs
c          
c           [rhoF,rhoP] = mm10_rhoFP_mrr(props, np1, n, tt, alpha);
          temp = (props%Gmat(alpha,1:props%num_hard)
     &      *tt(1:props%num_hard))
          rhoF = sum(temp)
          temp = (props%Hmat(alpha,1:props%num_hard)
     &      *tt(1:props%num_hard))
          rhoP = sum(temp)
c          
          tpass = c1*G*b*cdsqrt(rhoP) ! (16)
          ddipole = root3*G*b/(16.d0*pi*(one-v))/
     &       (rs) ! (42)
          rhoM = (two*k/(c1*c2*c3*G*b**3))*
     &       theta*cdsqrt(rhoF*rhoP) ! (13)
c          
c          call mm10_slipinc_ornl(props, np1, n, stress, tt, alpha, 
c     &     slipinc)
          slipinc = ivec1(alpha)
          slipinc = dsign(one,dreal(slipinc))*slipinc
          gammadot = slipinc/dt
      tem1 = c4/b*cdsqrt(rhoF)*gammadot
      tem2 = c6*ddipole/b*rhoM*gammadot
      tem3 = c5*rho*gammadot
      tem1 = c7*dexp(-Qbulk/k/theta)*rs/(k*theta)
     &        *rho*rho*gammadot**c8
          h(alpha) = rho_n + dt*(c4/b*cdsqrt(rhoP)*gammadot
     &        + c6*ddipole/b*rhoM*gammadot - c5*rho*gammadot
     &        - c7*dexp(-Qbulk/k/theta)*rs/(k*theta)
     &        *rho*rho*gammadot**c8) ! (18)
      enddo
c
      return
c
9000  format('>> FATAL ERROR: Roters model requires non-zero nodal',
     & /,     '                temperatures. terminated' )
      end 

c           Wrapper version, ornl slipinc function
      subroutine mm10_slipinc_ornlW(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt, zerosV
      double precision :: slipinc
      integer :: alpha
c
      complex(kind=real64), dimension(6) :: stressi
      complex(kind=real64), dimension(size_num_hard) :: tti
      complex(kind=real64) :: slipinci
c
      zerosV = zero
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
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, dtdstress
      double precision, dimension(size_num_hard) :: tt,rhoFs,rhoPs
      double precision, dimension(size_num_hard,6) :: et
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  rs,
     &  rhoF, rhoP, tpass, rho,
     &  c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  dddipole, dslipinc, badterm, G0
       !double precision, dimension(size_num_hard,size_num_hard)
       !&   :: Gmat, Hmat
      double precision, dimension(size_nslip) :: dslip
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
c
c compute derivatives of slip increments with respect to resolved
c shear stress
c        call mm10_dgdt_ornl(props, np1, n, stress, 
c     &         tt, dslip)
        dslip(1:props%num_hard) = arr1(1:props%num_hard,1)
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - one)
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Gmat,props%num_hard,tt,1,zero,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Hmat,props%num_hard,tt,1,zero,rhoPs,1)
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
          ddipole = root3*G*b/(16.d0*pi*(one-v))/
     &       (dabs(rs)) ! (42)
          dddipole = -root3*G*b/(16.d0*pi*(one-v))/
     &          (dabs(rs))**2*dsign(one,rs)
          rhoM = (two*k/(c1*c2*c3*G*b**3))*
     &       theta*dsqrt(rhoF*rhoP) ! (13)
c          
c          call mm10_slipinc_ornl(props, np1, n, stress, tt, alpha, 
c     &     slipinc)
          slipinc = vec1(alpha)
          gammadot = dabs(slipinc/dt)
          dtdstress(1:6) = np1%ms(1:6,alpha)
          dslipinc = dslip(alpha)/dt ! always positive

          ! Evaluate the equation
          if(gammadot.eq.zero) then
          badterm = zero
          else
          badterm = c8*dabs(rs)*gammadot**(c8-one)*dslipinc
          endif
          et(alpha,1:6) = dt*(c4/b*dsqrt(rhoF)*dsign(one,rs)
     &         *dslipinc + c6/b*rhoM*(ddipole*dsign(one,rs)
     &         *dslipinc + dddipole*gammadot)
     &        - c5*rho*dsign(one,rs)*dslipinc
     &    - c7*dexp(-Qbulk/k/theta)/(k*theta)*rho**2*
     &      (badterm + dsign(one,rs)*gammadot**c8))
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
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt,rhoPs,rhoFs
      double precision, dimension(size_num_hard,size_num_hard) :: etau
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, k, theta, G, b, c1, c2, c3, 
     &  rs,
     &  rhoF, rhoP, tpass, rho,
     &  c4, c5, c6, c7, c8, v, mm10_rs,
     &  ddipole, rhoM, slipinc, gammadot, Qbulk,
     &  dddipole, dslipinc, badterm, deltaij,
     &  drhoF, drhoP, drhoM,G0
       !double precision, dimension(size_num_hard,size_num_hard)
       !&   :: Gmat, Hmat
      double precision, dimension(size_num_hard,size_num_hard)
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
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - one)
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Gmat,props%num_hard,tt,1,zero,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Hmat,props%num_hard,tt,1,zero,rhoPs,1)

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
        ddipole = root3*G*b/(16.d0*pi*(one-v))/
     &       (dabs(rs)) ! (42)
        rhoM = (two*k/(c1*c2*c3*G*b**3))*
     &       theta*dsqrt(rhoF*rhoP) ! (13)
c          
c               loop over denominator hardening variable
        do beta = 1,props%num_hard
c          
c           [drhoF,drhoP] = mm10_drhoFP_ornl(props, np1, n, tt, alpha, beta);
          drhoF = props%Gmat(alpha,beta)
          drhoP = props%Hmat(alpha,beta)
          
          dddipole = zero
          drhoM = half*(two*k/(c1*c2*c3*G*b**3))*
     &       theta/dsqrt(rhoF*rhoP) * (drhoF*rhoP + rhoF*drhoP)
          
          dslipinc = dslip(alpha,beta)/dt
          
          if(alpha.eq.beta) then
              deltaij = one
          else
              deltaij = zero
          endif

c Evaluate the equation
          if(gammadot.eq.zero) then
          badterm = zero
          else
          badterm = c8*gammadot**(c8-one)*dslipinc
          endif
          etau(alpha,beta) = deltaij - dt*(c4/b*(half*drhoF/
     &        dsqrt(rhoF)*gammadot*dsign(one,rs) 
     &     + dsqrt(rhoF)*dslipinc) + c6/b*(rhoM*ddipole*dslipinc
     &     + rhoM*dddipole*gammadot*dsign(one,rs)
     &     + drhoM*ddipole*gammadot*dsign(one,rs))
     &     - c5*(rho*dslipinc + deltaij*gammadot*dsign(one,rs))
     &     - c7*dexp(-Qbulk/k/theta)/(k*theta)*dabs(rs)*
     &      (rho**2*badterm + two*rho*dsign(one,rs)
     &      *deltaij*gammadot**c8))*dsign(one,rs)
      
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
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(6,size_num_hard) :: ed
c
      ed = zero
c
      return
      end subroutine
c
c           Derivative of sliprate wrt resolved shear stress
      subroutine mm10_dgdt_ornl(props, np1, n, stress, tt, dgammadtau)
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision :: rs
      double precision, dimension(size_nslip) :: dgammadtau
      double precision, dimension(size_num_hard) :: tt,rhoFs,rhoPs
      double precision, external ::  mm10_rs
      integer :: alpha
c
      double precision :: dt, k, theta, G, b, c1, c2, tau0, 
     &  p_e, q_e, Qslip, v_s, dfract,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, m,
     &  dslipinc, slipexp, G0, fM, lamda, rhoM, p2
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
        if( v_s.lt.zero ) then
            v_s = exp(-v_s)
        endif
        fM = 0.1d0
        lamda = c2*b
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - one)
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Gmat,props%num_hard,tt,1,zero,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Hmat,props%num_hard,tt,1,zero,rhoPs,1)
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
        dfract = dsign(one,rs)/tcut
        if(fract.gt.one) then
c linear extrapolation past the too-high stress (rs) value
            b = gamma_0
            m = b * (-q_e*(Qslip/k/theta))
     &         * (- p_e) * dsign(one,rs)/tcut
            dslipinc = dt * m
        elseif(dabs(rs).eq.zero) then
            dslipinc = zero
        else
          if(fract.lt.zero) then
            p2 = one ! deal with low stresses by allowing the fraction to be small
          else
            p2 = p_e
          endif
        slipexp = dexp (-(Qslip/k/theta)*(one - fract**p2)**q_e)
     &          * dsign(one,rs)
        dslipinc = dt * (gamma_0 * slipexp * 
     &      (-one)*(Qslip/k/theta)*q_e*
     &      (one - fract**p2)**(q_e-one)
     &      * (-one)*p2*fract**(p2-one) * dfract) !(14)
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
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt,rhoPs,rhoFs
      double precision, dimension(size_nslip,size_num_hard)
     &    :: dgammadtt
      double precision, external:: mm10_rs
      double precision :: rs
      integer :: alpha, beta
c
      double precision :: dt, k, theta, G, b, c1, c2, tau0, 
     &  p_e, q_e, Qslip, v_s, dfract,
     &  rhoF, rhoP, gamma_0, tpass, tcut, fract, 
     &  dslipinc, slipexp, drhoF, drhoP, dgamma_0,
     &  dtcut, dtpass, G0, fM, rhoM, lamda, drhoM, p2
       !double precision, dimension(size_num_hard,size_num_hard)
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
        if( v_s.lt.zero ) then
            v_s = exp(-v_s)
        endif
        fM = 0.1d0
        lamda = c2*b
c        
c New shear modulus
        G = G0 - props%D_0 / (exp(props%T_0/theta) - one)
c Load the interaction matrices for parallel and forest dislocs
c        [Gmat,Hmat] = mm10_mrr_GH(props);
      ! call mm10_mrr_GH(props,Gmat,Hmat)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Gmat,props%num_hard,tt,1,zero,rhoFs,1)
      call dgemv('N',props%num_hard,props%num_hard,one,
     &           props%Hmat,props%num_hard,tt,1,zero,rhoPs,1)
        
c        Compute derivative of slip rate alpha w.r.t. density beta
c        loop over slip rate
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
             drhoM = zero
          endif
        
          dgamma_0 = drhoM*b*v_s/b*lamda ! (15)
          dtpass = half*c1*G*b/dsqrt(rhoP)*drhoP ! (16)
          dtcut = zero ! (17)
          dfract = ((-dtpass)*tcut - 
     &             (dabs(rs)-tpass)*dtcut)/tcut**2

c Evaluate the equation
          if(fract.gt.one) then
              ! linear extrapolation past the too-high stress (rs) value
              dslipinc = dt * (-q_e*(Qslip/k/theta)
     &            * (- p_e))
     &            * dsign(one,rs) * (dgamma_0/tcut 
     &             - gamma_0*dtcut/tcut**2)
          elseif(dabs(rs).eq.zero) then
              dslipinc = zero
          else
            if(fract.lt.zero) then
              p2 = one
            else
              p2 = p_e
            endif
          slipexp = dexp (-(Qslip/k/theta)*
     &             (one - fract**p2)**q_e) !(14)
          dslipinc = dt * (dgamma_0 * slipexp * dsign(one,rs)
     &          + gamma_0 * dsign(one,rs) * slipexp *
     &          (-one)*(Qslip/k/theta)*q_e*(one - fract**p2)**(q_e-one)
     &        * (-one)*p2*fract**(p2-one) * dfract)
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
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, D
      double precision, dimension(6,size_nslip) :: dgammadd
      double precision, dimension(size_num_hard) :: tt
c
      dgammadd = zero
c
      return
      end subroutine
c
c
c *****************************************************************************
c *                                                                           *
c *         Armstrong-Frederick hardening model                               *
c *                                                                           *
c *****************************************************************************
c
c
c
c
c     ARFR:
c
c           Form intermediate vectors for faster calculations
      subroutine mm10_v_arfr(props, np1,n, stress, tt,
     &   vec1, vec2)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
      integer :: alpha
      double precision :: slipinc


      do alpha = 1,props%nslip
         call mm10_slipinc_arfr(props, np1, n, stress, tt, 
     &                            alpha, slipinc)
         vec1(alpha) = slipinc
      enddo

      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_a_arfr                       *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 05/18/2017 tjt              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_a_arfr( props, np1, n, stress, tt, vec1, vec2,
     &                        arr1, arr2, both )
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      integer :: both
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      integer :: nslip, num_hard  
      double precision, allocatable :: temp_arr2(:,:)    
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, arr1:64, arr2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64     
c
c              arr1 treated as vector inside
c
      nslip    = props%nslip
      num_hard = props%num_hard
c      
      call mm10_dgdt_arfr( props, np1, n, stress, tt, arr1(1,1) )
c
c              arr2 must have lead dimension = nslip for compatibility 
c              with _dgdh_arfr. values defined inside the routine - no
c              need to copyin. only copyout needed
c
      if( both == 2 ) then
         allocate( temp_arr2(nslip,num_hard) )
         call mm10_dgdh_arfr( props, np1, n, stress, tt, temp_arr2 )
         arr2(1:nslip,1:num_hard) = temp_arr2
         deallocate( temp_arr2 ) 
      end if
c      
      return
      end
c
c           Calculate the slip increment along system i  
c           arfr model
      subroutine mm10_slipinc_arfr(props, np1, n, stress, tt, 
     &                             alpha, slipinc)
c
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision, external :: mm10_rs
      double precision :: slipinc
      integer :: alpha
c
      double precision :: dt, tau, g_alpha
      double precision ::  gamma_dot_tilde, m_alpha
        
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
     &  abs(tau/g_alpha)**(one/m_alpha)*dsign(one,tau)
      
      return
      end
c
c
c     Hardening function for arfr model
c
      subroutine mm10_h_arfr( props, np1,
     &             n, vec1, vec2, stress, tt, g )
c
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt, g, h
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision :: slipinc
      integer :: slip_a, slip_b
c
      double precision :: dt
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
            g_s = max(temp, one/four*g_0_alpha) ! threshold to ensure no slip during elastic response
            h(slip_b) = h_0_alpha * abs(one-tt(slip_b)/g_s)
     &      **r_alpha * sign(one,one-tt(slip_b)/g_s)
        enddo
        
        ! Equation [5]
        do slip_a = 1,props%nslip
            g_dot = zero
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
c     ARFR:
c
c           Derivative of hardening fn wrt stress
      subroutine mm10_estress_arfr(props,
     & np1, n,vec1,vec2,arr1,arr2, stress, tt, et)
 
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, dtdstress
      double precision, dimension(size_num_hard) :: tt,h,dh
      double precision, dimension(size_num_hard,6) :: et
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, slipinc
      double precision :: h_0_alpha, gamma_dot_tilde,
     &     g_tilde, r_alpha, n_alpha, m_alpha, g_0_alpha,
     &     g_s, gamma_dot, temp, dg_s, dslip
      double precision, dimension(size_nslip) :: dslipinc
      integer :: slip_a, slip_b
      
      et = zero

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
            if( temp .gt. one/four*g_0_alpha ) then ! threshold to ensure no slip during elastic response
                g_s = temp
                h(slip_b) = h_0_alpha * abs(one-tt(slip_b)/g_s)
     &          **r_alpha * sign(one,one-tt(slip_b)/g_s)
                dg_s = g_tilde * n_alpha / gamma_dot * 
     &              (gamma_dot/gamma_dot_tilde)**n_alpha * dslip
                dh(slip_b) = h_0_alpha * r_alpha * 
     &              abs(one-tt(slip_b)/g_s)**(r_alpha-one)
     &              *(dg_s*tt(slip_b)/g_s**2)
            else
                g_s = one/four*g_0_alpha
                h(slip_b) = h_0_alpha * abs(one-tt(slip_b)/g_s)
     &          **r_alpha * sign(one,one-tt(slip_b)/g_s)
                dh(slip_b) = zero
            endif
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
c     ARFR:
c
c           Derivative of hardening fn wrt hardening
        subroutine mm10_ehard_arfr(props,
     &      np1, n,vec1,vec2,arr1,arr2, stress, tt, etau)
 
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt,h,dh
      double precision, dimension(size_num_hard,size_num_hard) :: etau
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, slipinc, deltaij
      double precision :: h_0_alpha, gamma_dot_tilde,
     &     g_tilde, r_alpha, n_alpha, m_alpha, g_0_alpha,
     &     g_s, gamma_dot, temp, dg_s, dslip
      double precision, dimension(size_nslip,size_num_hard)
     &        :: dslipinc
      integer :: slip_a, slip_b
      
       etau = zero

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
            if( temp .gt. one/four*g_0_alpha ) then ! threshold to ensure no slip during elastic response
                g_s = temp
                h(slip_b) = h_0_alpha * abs(one-tt(slip_b)/g_s)
     &          **r_alpha * sign(one,one-tt(slip_b)/g_s)
                dg_s = g_tilde * n_alpha / gamma_dot * 
     &              (gamma_dot/gamma_dot_tilde)**n_alpha * dslip
                dh(slip_b) = h_0_alpha * r_alpha * 
     &              abs(one-tt(slip_b)/g_s)**(r_alpha-one)
     &              *(-one/g_s + dg_s*tt(slip_b)/g_s**2)
            else
                g_s = one/four*g_0_alpha
                h(slip_b) = h_0_alpha * abs(one-tt(slip_b)/g_s)
     &          **r_alpha * sign(one,one-tt(slip_b)/g_s)
                dh(slip_b) = zero
            endif
            
        enddo
        
        ! Equation [5]
        do slip_a = 1,props%nslip
            do slip_b = 1,props%nslip
                  if( slip_a.eq.slip_b ) then
                      deltaij = one
                  else
                      deltaij = zero
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
      subroutine mm10_ed_arfr(props, np1, n, stress, tt, ed)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(6,size_num_hard) :: ed
c
      ed = 0.d0
c
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_dgdt_arfr                    *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 05/18/17                    *
c     *                                                              *
c     *     Calculate partial gamma(s) w.r.t. tau(s) for each slip   *
c     *     system, for use in J11 in material integration. ARFR     *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_dgdt_arfr(props, np1,  n, stress, tt, dgdt)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt, dgdt
      double precision, external:: mm10_rs
      integer :: slip_a
c
      double precision :: dt, tau, g_alpha
      double precision :: gamma_dot_tilde, m_alpha
        
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
c     *                 subroutine mm10_dgdh_arfr                    *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 05/18/17                    *
c     *                                                              *
c     *     Calculate partial gamma(s) w.r.t. tt for each slip       *
c     *     system, for use in J12 in material integration. ARFR     *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_dgdh_arfr(props, np1,
     &      n, stress, tt, dgammadtt)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(size_num_hard,size_num_hard)
     &                 :: dgammadtt
      double precision, external :: mm10_rs
      double precision :: dslipinc
      integer :: slip_a, slip_b
c
      double precision :: dt, tau, g_alpha
      double precision ::  gamma_dot_tilde, m_alpha
        
        dgammadtt = 0.d0
        
        dt = np1%tinc
        
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
      subroutine mm10_dgdd_arfr(props, np1, n, stress, tt, D, 
     &              dgammadd)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, D
      double precision, dimension(6,size_nslip) :: dgammadd
      double precision, dimension(size_num_hard) :: tt
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
      double precision, dimension(size_num_hard) :: tt
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
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_a_DJGM                       *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/15/2016 rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_a_DJGM( props, np1, n, stress, tt, vec1, vec2,
     &                        arr1, arr2, both )
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      integer :: both
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      integer :: nslip, num_hard  
      double precision, allocatable :: temp_arr2(:,:)    
c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, arr1:64, arr2:64, stress:64
c!DIR$ ASSUME_ALIGNED tt:64     
c
c              arr1 treated as vector inside
c
      nslip    = props%nslip
      num_hard = props%num_hard
c      
      call mm10_dgdt_DJGM( props, np1, n, stress, tt, arr1(1,1) )
c
c              arr2 must have lead dimension = nslip for compatibility 
c              with _dgdh_DJGM. values defined inside the routine - no
c              need to copyin. only copyout needed
c
      if( both == 2 ) then
         allocate( temp_arr2(nslip,num_hard) )
         call mm10_dgdh_DJGM( props, np1, n, stress, tt, temp_arr2 )
         arr2(1:nslip,1:num_hard) = temp_arr2
         deallocate( temp_arr2 ) 
      end if
c      
      return
      end
c
c           Calculate the slip increment along system i  
c           mrr model
      subroutine mm10_slipinc_DJGM(props, np1, n, stress, tt, 
     &                             alpha, slipinc)
c
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt
      double precision, external :: mm10_rs
      double precision :: slipinc
      integer :: alpha
c
      double precision :: dt, tau, g_alpha
      double precision :: gamma_dot_tilde, m_alpha
        
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
     &  abs(tau/g_alpha)**(one/m_alpha)*dsign(one,tau)
      
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
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt, g, h
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision :: slipinc
      integer :: slip_a, slip_b
c
      double precision :: dt
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
            g_s = max(temp, one/four*g_0_alpha) ! threshold to ensure no slip during elastic response
            h(slip_b) = h_0_alpha * abs(one-tt(slip_b)/g_s)
     &      **r_alpha * sign(one,one-tt(slip_b)/g_s)
        enddo
        
        ! Equation [5]
        do slip_a = 1,props%nslip
            g_dot = zero
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
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress, dtdstress
      double precision, dimension(size_num_hard) :: tt,h,dh
      double precision, dimension(size_num_hard,6) :: et
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, slipinc
      double precision :: h_0_alpha, gamma_dot_tilde,
     &     g_tilde, r_alpha, n_alpha, m_alpha, g_0_alpha,
     &     g_s, gamma_dot, temp, dg_s, dslip
      double precision, dimension(size_nslip) :: dslipinc
      integer :: slip_a, slip_b
      
      et = zero

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
            if( temp .gt. one/four*g_0_alpha ) then ! threshold to ensure no slip during elastic response
                g_s = temp
                h(slip_b) = h_0_alpha * abs(one-tt(slip_b)/g_s)
     &          **r_alpha * sign(one,one-tt(slip_b)/g_s)
                dg_s = g_tilde * n_alpha / gamma_dot * 
     &              (gamma_dot/gamma_dot_tilde)**n_alpha * dslip
                dh(slip_b) = h_0_alpha * r_alpha * 
     &              abs(one-tt(slip_b)/g_s)**(r_alpha-one)
     &              *(dg_s*tt(slip_b)/g_s**2)
            else
                g_s = one/four*g_0_alpha
                h(slip_b) = h_0_alpha * abs(one-tt(slip_b)/g_s)
     &          **r_alpha * sign(one,one-tt(slip_b)/g_s)
                dh(slip_b) = zero
            endif
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
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt,h,dh
      double precision, dimension(size_num_hard,size_num_hard) :: etau
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
c
      double precision :: dt, slipinc, deltaij
      double precision :: h_0_alpha, gamma_dot_tilde,
     &     g_tilde, r_alpha, n_alpha, m_alpha, g_0_alpha,
     &     g_s, gamma_dot, temp, dslip, dg_s
      double precision, dimension(size_nslip,size_num_hard)
     &        :: dslipinc
      integer :: slip_a, slip_b
      
       etau = zero

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
            if( temp .gt. one/four*g_0_alpha ) then ! threshold to ensure no slip during elastic response
                g_s = temp
                h(slip_b) = h_0_alpha * abs(one-tt(slip_b)/g_s)
     &          **r_alpha * sign(one,one-tt(slip_b)/g_s)
                dg_s = g_tilde * n_alpha / gamma_dot * 
     &              (gamma_dot/gamma_dot_tilde)**n_alpha * dslip
                dh(slip_b) = h_0_alpha * r_alpha * 
     &              abs(one-tt(slip_b)/g_s)**(r_alpha-one)
     &              *(-one/g_s + dg_s*tt(slip_b)/g_s**2)
            else
                g_s = one/four*g_0_alpha
                h(slip_b) = h_0_alpha * abs(one-tt(slip_b)/g_s)
     &          **r_alpha * sign(one,one-tt(slip_b)/g_s)
                dh(slip_b) = zero
            endif
            
        enddo
        
        ! Equation [5]
        do slip_a = 1,props%nslip
            do slip_b = 1,props%nslip
                  if( slip_a.eq.slip_b ) then
                      deltaij = one
                  else
                      deltaij = zero
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
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(6,size_num_hard) :: ed
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
      subroutine mm10_dgdt_DJGM(props, np1, n, stress, tt, dgdt)
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, n
      double precision, dimension(6) :: stress
      double precision, dimension(size_num_hard) :: tt, dgdt
      double precision, external :: mm10_rs
      integer :: slip_a
c
      double precision :: dt, tau, g_alpha
      double precision ::  gamma_dot_tilde, m_alpha
        
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
      double precision, dimension(size_num_hard) :: tt
      double precision, dimension(size_num_hard,size_num_hard)
     &                 :: dgammadtt
      double precision, external:: mm10_rs
      double precision :: dslipinc
      integer :: slip_a, slip_b
c
      double precision :: dt, tau, g_alpha
      double precision :: gamma_dot_tilde, m_alpha
        
        dgammadtt = 0.d0
        
        dt = np1%tinc
        
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
      double precision, dimension(6,size_nslip) :: dgammadd
      double precision, dimension(size_num_hard) :: tt
c
      dgammadd = 0.d0
c
      return
      end subroutine
c
c
c     Halite/diffusion functions
c
c            formR for pressure-precipitation model
      subroutine mm10_halite_formRpp( props, work_vec1, stress, tinc )
      use mm10_defs
      use mm10_constants
      implicit none
c
c      global variables
      type(crystal_props) :: props
      double precision, dimension(6) :: stress, work_vec1
      double precision :: tinc
c      local variables
      double precision :: I_dev(6, 6), eps_norton(6), temp(6)
c
      I_dev(1:3, 4:6) = zero
      I_dev(4:6, 1:3) = zero
      I_dev(1, 1) = two/three
      I_dev(1, 2) = -one/three
      I_dev(1, 3) = -one/three
      I_dev(2, 1) = -one/three
      I_dev(2, 2) = two/three
      I_dev(2, 3) = -one/three
      I_dev(3, 1) = -one/three
      I_dev(3, 2) = -one/three
      I_dev(3, 3) = two/three
      I_dev(4, 4) = one
      I_dev(4, 5) = zero
      I_dev(4, 6) = zero
      I_dev(5, 4) = zero
      I_dev(5, 5) = one
      I_dev(5, 6) = zero
      I_dev(6, 4) = zero
      I_dev(6, 5) = zero
      I_dev(6, 6) = one
      call mm10_b_mult_type_2( temp, I_dev, stress )
      eps_norton = three/two * props%cp_033 * tinc * temp
      eps_norton(4:6) = eps_norton(4:6) * two
      work_vec1 = work_vec1 - eps_norton
      end subroutine
c
c            formJ for pressure-precipitation model
      subroutine mm10_halite_formJpp( props, J11, tinc )
      use mm10_defs
      use mm10_constants
      implicit none
c
      type(crystal_props) :: props
      double precision :: J11(6, 6)
      double precision :: tinc
c      local variables
      double precision :: I_dev(6, 6), temp(6, 6)
      integer :: i
c
      I_dev(1:3, 4:6) = zero
      I_dev(4:6, 1:3) = zero
      I_dev(1, 1) = two/three
      I_dev(1, 2) = -one/three
      I_dev(1, 3) = -one/three
      I_dev(2, 1) = -one/three
      I_dev(2, 2) = two/three
      I_dev(2, 3) = -one/three
      I_dev(3, 1) = -one/three
      I_dev(3, 2) = -one/three
      I_dev(3, 3) = two/three
      I_dev(4, 4) = one
      I_dev(4, 5) = zero
      I_dev(4, 6) = zero
      I_dev(5, 4) = zero
      I_dev(5, 5) = one
      I_dev(5, 6) = zero
      I_dev(6, 4) = zero
      I_dev(6, 5) = zero
      I_dev(6, 6) = one
      do i = 1, 6
        call mm10_b_mult_type_2(temp(1:6, i), 
     &            props%stiffness, I_dev(1:6, i))
      end do
      temp(4:6, 1:6) = temp(4:6, 1:6) * two
      J11 = J11 + three/two*props%cp_033*temp*tinc
c
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_b_mult...                    *
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
      subroutine mm10_b_mult_type_1( a, b, c )
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
      subroutine mm10_b_mult_type_1a( a, b )
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
      subroutine mm10_b_mult_type_2( a, b, c )
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

c
      subroutine mm10_b_mult_type_2a( a, b, c, n )
      implicit none
      integer :: j, n
      double precision :: a(6), b(6,n), c(n)
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
      if ( n == 1 ) return

      do j = 2, n
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

      subroutine mm10_b_mult_type_2b( a, b, c, n )
      implicit none
      integer :: j, n
      double precision :: a(3), b(3,n), c(n)
c!DIR$ ASSUME_ALIGNED a:64, b:64, c:64      
c
c                     a = [b] * c
c
      a(1) = b(1,1)*c(1)
      a(2) = b(2,1)*c(1)
      a(3) = b(3,1)*c(1)
      if ( n == 1 ) return

      do j = 2, n
        a(1) = a(1) + b(1,j)*c(j)
        a(2) = a(2) + b(2,j)*c(j)
        a(3) = a(3) + b(3,j)*c(j)
      end do
c
      return
      end


      subroutine mm10_b_mult_type_4( a, b, c, d, const )
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


      subroutine mm10_b_mult_type_2t( a, b, c )
      implicit none
      double precision :: a(6), b(6,6), c(6)
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

      subroutine mm10_b_mult_type_3t( a, b, c )
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

      subroutine mm10_b_mult_type_5( a, b, c )
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


      subroutine mm10_b_mult_type_3( a, b, c )
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
c     *                 subroutine mm10_b_zero_vec                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 9/29/2016 rhd               *
c     *                                                              *
c     *    zero a vector - should be inlined                         *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_b_zero_vec( vec, nterms )
      implicit none
      integer :: nterms
      double precision :: vec(nterms)
      double precision, parameter :: zero = 0.0d0
c!DIR$ ASSUME_ALIGNED vec:64       
c
      vec = zero
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_b_make_symm_1                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 9/29/2016 rhd               *
c     *                                                              *
c     *                 make a symmetric 6x6 matrix                  *
c     *                                                              *
c     ****************************************************************
c

      subroutine mm10_b_make_symm_1( a )
      implicit none
c
      integer :: i, j      
      double precision :: a(6,6), work(6,6)
      double precision, parameter :: half = 0.5d0
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
c     *                 subroutine mm10_b_copy_vec                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 9/29/2016 rhd               *
c     *                                                              *
c     *     copy a vector - should be inlined                        *
c     *                                                              *
c     ****************************************************************
c
 
      subroutine mm10_b_copy_vector( a, b, nterms )
      implicit none
      integer :: nterms
      double precision :: a(nterms), b(nterms)
c!DIR$ ASSUME_ALIGNED a:64, b:64
c
      a = b
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm10_b_invert_33                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 10/11/2016 rhd              *
c     *                                                              *
c     *                     invert a general 3 x 3 matrix            *
c     *                                                              *
c     ****************************************************************
c
 
 
      subroutine mm10_b_invert_33( ainv, a,  ok_flag )
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

      
