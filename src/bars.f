c     ****************************************************************
c     *                                                              *
c     *                      subroutine bar_mass                     *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 02/07/17                   *
c     *                                                              *
c     *     compute bar mass matrices                                *
c     *                                                              *
c     ****************************************************************
c
      subroutine bar_mass(span, props, emass, mel, volume, totdof)
      implicit none
      
      include 'param_def'

      integer :: span, totdof
      real :: props(mxelpr,*)
      double precision :: volume(*), mel(totdof,*), emass(*)

      integer :: i,j

      do i=1,span
         emass(i) = props(44, i)
         do j=1,totdof
            mel(j,i) = props(44,i) / 2
         end do
      end do


      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine bar_stiffness_symm           *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 02/07/17                   *
c     *                                                              *
c     *     compute symmetric bar stiffness                          *
c     *                                                              *
c     ****************************************************************
c
      subroutine bar_stiffness_symm(local_work, ek_symm, span, nrow_ek)
      implicit none
      include 'param_def'
      include 'include_tan_ek'
      
      integer :: span, nrow_ek
      double precision :: ek_symm(span, nrow_ek)
      double precision, allocatable, dimension(:,:,:) :: ek_full
      
      integer :: i,j,k,l

      ! Fix this later, but for now be very lazy
      allocate(ek_full(span, 6, 6))

      call bar_stiffness_asymm(local_work, ek_full, span, 36)

      do i=1,span
            ! Arg, nowhere do we write down how this is supposed
            ! to get unrolled
            l = 1
            do j=1,6    ! col
                  do k=1,j    ! row
                        ek_symm(i,l) = ek_full(i,k,j)
                        l = l + 1
                  end do
            end do
      end do

      deallocate(ek_full)

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine bar_stiffness_asymm          *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 02/07/17                   *
c     *                                                              *
c     *     compute asymmetric bar stiffness                         *
c     *                                                              *
c     ****************************************************************
c
      subroutine bar_stiffness_asymm(local_work, ek_asymm, span,
     &            nrow_ek)
      implicit none
      include 'param_def'
      include 'include_tan_ek'

      integer :: span, nrow_ek
      double precision :: ek_asymm(span, nrow_ek)
      
      integer :: i

      do i=1,span
            call bar_stiffness_single(local_work%ce(i,1:6), 
     &            local_work%bar_stiffness(i), ek_asymm(i,:))
      end do
      
      end subroutine


c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine bar_stiffness_single         *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 02/07/17                   *
c     *                                                              *
c     *     compute a single  6x6 bar stiffness matrix               *
c     *                                                              *
c     ****************************************************************
c
      subroutine bar_stiffness_single(coords, k, ek)
      implicit none
      
      double precision :: coords(2,3), ek(6,6)
      real :: k
c
c           "Real" bars
c

      double precision :: n(3)
      double precision :: T(2,6)
      double precision :: lK(2,2)
      integer :: i

c      n = coords(2,:) - coords(1,:)
c      n = n / sqrt(dot_product(n,n))
c      
c      T = 0.0
c      T(1,1) = n(1)
c      T(2,2) = n(1)
c      T(1,3) = n(2)
c      T(2,4) = n(2)
c      T(1,5) = n(3)
c      T(2,6) = n(3)
c
c      lK(1,1) = k
c      lK(1,2) = -k
c      lK(2,1) = -k
c      lK(2,2) = k
c      
c      ek = matmul(transpose(T),matmul(lK,T))

c           "Penalty" bars
      ek = 0.0
      ek(1,1) =  k
      ek(1,2) = -k
      ek(2,1) = -k
      ek(2,2) =  k

      ek(3,3) =  k
      ek(3,4) = -k
      ek(4,3) = -k
      ek(4,4) =  k

      ek(5,5) =  k
      ek(5,6) = -k
      ek(6,5) = -k
      ek(6,6) =  k


      end subroutine

c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine bar_forces                   *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 02/07/17                   *
c     *                                                              *
c     *     computes the bar element force vectors                   *
c     *                                                              *
c     ****************************************************************
c
      subroutine bar_forces(local_work, forces, volumes, span)
      implicit none
      include 'param_def'
      include 'include_sig_up'      

      integer :: span
      double precision :: forces(span, 6)
      double precision :: volumes(span)

      integer :: i
      double precision, dimension(6,6) :: k_elem
      double precision, dimension(6):: u_elem

      ! Be faster later
      do i=1,span
            ! Volume still zero
            volumes(i) = 0.0

            ! Get element stiffness
            call bar_stiffness_single(local_work%ce_0(i,1:6),
     &            local_work%bar_stiffness(i), k_elem)
            u_elem = local_work%ue(i,1:6) + local_work%due(i,1:6)
            forces(i,1:6) = matmul(k_elem, u_elem)
      end do

      end subroutine


      
