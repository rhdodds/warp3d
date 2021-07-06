      implicit double precision (a-h,o-z)
c
      integer, parameter :: max_nodes = 100000
      dimension x(max_nodes), y(max_nodes), z(max_nodes)
      integer :: bound_nodes(max_nodes)
      double precision :: nu, KI, kappa, rad, zero
      logical, parameter :: debug = .false.

      factor = 1.0d0  ! already in mm

      radius_mm = 2000.d0
      radius_m = 2.0d0
      pi = 3.14159265359d0
      zero = 0.0d0
c
      read(*,*) nnode
c
      do i = 1, nnode
       read(*,*,end=100) j, x(j), y(j), z(j)
       x(i) = x(i) * factor
       y(i) = y(i) * factor
       z(i) = z(i) * factor
      end do
      write(*,*) '.... coordinates read '
c
c             find boundary nodes (R = 2000.)
c
      nboundary = 0
      do i = 1, nnode
       rad = sqrt( x(i)*x(i) + y(i)*y(i) )
       if( rad < 0.9d0*radius_mm ) cycle
       nboundary = nboundary + 1
       bound_nodes(nboundary) = i
      end do
      write(*,*) '.... # boundary nodes: ', nboundary
c
c             fix radius on boundary nodes
c
      do i = 1, nboundary
        node = bound_nodes(i)
        xb = x(node)
        yb = y(node)
        angle = atan2(yb,xb)
        x(node) = radius_mm * cos(angle)
        y(node) = radius_mm * sin(angle)
      end do
c
      open(unit=11,file="coords_mm.inp")
      write(11,1002) "*echo off"
      write(11,1002) " coordinates"
      do i = 1, nnode
       write(11,1000) i, x(i), y(i), z(i)
      end do
      write(11,1002) "*echo on"
      close(unit=11)
c
c             du, dv displacements for K_I = 1 MPa sqrt(m) on bound
c
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
        write(*,2000) node, u_k, v_k ! mm
      end do
      
c
 100  stop
 1000 format(i6,3e17.9)
 1002 format(a)
 2000 format(i9,1x,'u ',e17.9,3x,'v ',e17.9)
      end
