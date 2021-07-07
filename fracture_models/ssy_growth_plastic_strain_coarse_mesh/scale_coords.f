      implicit double precision (a-h,o-z)
c
      dimension x(3000), y(3000), z(3000)
      integer :: bound_nodes(500)
      double precision :: nu, KI, kappa

      factor = 25.4d0 * 0.9842517962599695d0

      radius_mm = 2000.d0
      radius_m = 2.0d0
      pi = 3.14159265359d0
c
      read(*,*) nnode

      do i = 1, nnode
       read(*,*,end=100) j, x(j), y(j), z(j)
       x(i) = x(i) * factor
       y(i) = y(i) * factor
       z(i) = z(i) * factor
      end do

      read(*,*) nboundary
      do i = 1, nboundary
        read(*,*) bound_nodes(i)
      end do

      do i = 1, nboundary
        node = bound_nodes(i)
        xb = x(node)
        yb = y(node)
        angle = atan2(yb,xb)
        x(node) = radius_mm * cos(angle)
        y(node) = radius_mm * sin(angle)
      end do

      do i = 1, nnode
       write(*,1000)  i, x(i), y(i), z(i)
      end do

      stop  !   <<<<<<<<<<<<<<<<<<<<<<<<

      ymod = 88000.d0
      nu = 0.34d0
      gmod = ymod / 2.0d0 / (1.0d0 + nu )
      KI = 1.0d0
      applied_k = 1.0d0
      write(*,*) ' '
      write(*,*) ' '

c
c            make KI displacements on boundary for 1.0 MPa-sqrt(m)
c            Zr4 properties.
c            compute in meters, final in mm

c
      do i = 1, nboundary
        node = bound_nodes(i)
        xb = x(node)
        yb = y(node)
        theta = atan2( yb, xb ) 
        ct2 = cos(theta*0.5d0)
        st2 = sin(theta*0.5d0)
        ct  = cos(theta)
        st  = sin(theta)
        term1 = (3.0-4.0*nu-ct) * ct2
        term2 = (3.0-4.0*nu-ct) * st2
        term3 = sqrt( radius_m * 0.5d0 / pi )
        term4 = (1.0d0 + nu ) * applied_k / ymod
        u_k = 1000.0d0 * term4 * term3 * term1
        v_k = 1000.0d0 * term4 * term3 * term2
        write(*,2000) node, u_k, v_k
      end do
      
c
 100  stop
 1000 format(i6,3e17.9)
 2000 format(i9,1x,'u ',e17.9,3x,'v ',e17.9)
      end
