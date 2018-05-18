      implicit none
c
      integer :: i, node
      double precision :: gamma, beta, alpha, pi, toler, gamma_r,
     &                    beta_r,  alpha_r, c, s, mag
      double precision, parameter :: zero = 0.0d0, one=1.0d0, two=2.0d0
      double precision :: rotx(3,3), roty(3,3), rotz(3,3), q(3,3),
     &                    qchk(3,3),  x_old(3), x_new(3), unit_y(3),
     &                    p1(3), p2(3), p3(3), v1(3), v2(3), v3(3),
     &                    q_local_global(3,3), q_global_local(3,3)
c
c                   given rotated coordinates, compute the local->global
c                   rotation, then transpose is global->local for
c                   constraint transformation matrix
c
c                   coordinates of 3 nodes in rotated system
c
      p1(1) =-0.118811416E+01
      p1(2) = 0.672382423E+00
      p1(3) = 0.326016658E+00


      p2(1) = -0.115275882E+01
      p2(2) = 0.707737768E+00
      p2(3) =  0.412619212E+00

      p3(1) =-0.127550476E+01
      p3(2) = 0.717884314E+00
      p3(3) = 0.343117650E+00

      v1 = p2 - p1
      v2 = p3 - p1
      call normalize( v1, mag )
      call normalize( v2, mag )
      call cross_prod( v1, v2, v3 )
      call normalize( v3, mag )
      q_local_global(1:3,1) = v1
      q_local_global(1:3,2) = v2
      q_local_global(1:3,3) = v3
      q_global_local = transpose( q_local_global )
       write(*,*) '... q_global_local ...'
      do i = 1, 3
         write(*,9100) q_global_local(i,1:3)
      end do
      unit_y(1) = zero
      unit_y(2) = one
      unit_y(3) = zero
c
      unit_y = matmul( q_global_local, unit_y )
      write(*,*) '... unit_y rotated ...'
      write(*,*) unit_y

c
c              compute rotation to original coordinates to
c              give general position in space
c

      gamma = 20.0d0
      beta  = -60.0d0
      alpha = 45.0d0
      pi = atan(one)*4.0d0
c
      toler = 0.00001d0
c
      rotx = zero
      roty = zero
      rotz = zero
      q    = zero
      qchk = zero

      gamma_r = gamma * pi/180.d0
      beta_r  = beta *  pi/180.d0
      alpha_r = alpha * pi/180.d0
c
      c = cos(gamma_r)
      s = sin(gamma_r)
      rotx(1,1) = one
      rotx(2,2) = c
      rotx(3,3) = c
      rotx(3,2) = s
      rotx(2,3) = -s

      c = cos(beta_r)
      s = sin(beta_r)
      roty(1,1) = c
      roty(1,3) = s
      roty(2,2) = one
      roty(3,1) = -s
      roty(3,3) = c

      c = cos(alpha_r)
      s = sin(alpha_r)
      rotz(1,1) = c
      rotz(1,2) = -s
      rotz(2,1) = s
      rotz(2,2) = c
      rotz(3,3) = one

      q = matmul( roty, rotx )
      q = matmul( rotz, q )
      write(*,*) '... Q ...'
      do i = 1, 3
         write(*,9100) q(i,1:3)
      end do



      write(*,*) '... trans(q) * q ...'
      qchk = matmul( transpose( q ), q )
      do i = 1, 3
         write(*,9000) qchk(i,1:3)
      end do

      do i = 1, 3020
       read(*,*,end=100) node, x_old
       x_new = matmul( q, x_old )
       if( abs(x_new(1) ) <= toler ) x_new(1) = zero
       if( abs(x_new(2) ) <= toler ) x_new(2) = zero
       if( abs(x_new(3) ) <= toler ) x_new(3) = zero
       write(*,1000) node, x_new
      end do
c
 100  stop
 1000 format(i8,3e17.9)
 9000 format(10x,3f10.3)
 9100 format(10x,3e16.8)
      end

      subroutine cross_prod (vec1, vec2, vec_out)
c
c
      double precision
     &     vec1(3), vec2(3), vec_out(3)
      data zero /0.0/
c
      vec_out(1) = vec1(2)*vec2(3) -
     &     vec2(2)*vec1(3)
      vec_out(2) = vec2(1)*vec1(3) -
     &     vec1(1)*vec2(3)
      vec_out(3) = vec1(1)*vec2(2) -
     &     vec2(1)*vec1(2)
c
      return
      end
c
c
      subroutine normalize (vec, mag)
c
c
      double precision
     &     vec(3), zero, mag
      integer i, j
      data zero /0.0/
c
      mag = sqrt ( vec(1)**2 + vec(2)**2 + vec(3)**2 )
c
      do i=1, 3
         vec(i) = vec(i) / mag
      enddo
c
      return
      end
