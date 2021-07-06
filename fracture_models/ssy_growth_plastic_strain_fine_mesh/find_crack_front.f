      implicit none


      double precision, allocatable :: coords(:,:), xbar(:), ybar(:),
     &       zbar(:), stresses(:,:), states(:,:), max_x(:), max_y(:),
     &       max_z(:)
      double precision :: zero, mises, x_0, crack_extension,
     &                    mesh_limit_xbar, max_front, minimum_y_bar,
     &                    mises_min, large
      integer, allocatable :: incid(:,:)
      integer :: nnode, nelem, i, j, out, step, lnode, snode, 
     &           num_sig_values, num_states_values
      character(len=30) :: fname1, fname2
      logical :: test, debug

      out = 6
      zero = 0.0d0
      x_0 = 0.0d0  ! global X-coord of straight crack front
      num_sig_values = 26
      num_states_values = 9
      debug = .false.
      mesh_limit_xbar = 11.0d0 !   of uniform mesh
      minimum_y_bar = 0.04d0 ! ignore element w/ Ybar >
      mises_min = 0.1d0
      large = 1.0d20
c
      open(unit=10,file='../raw_coords.dat')
      read(10,*) nnode
      allocate( coords(3,nnode) )
      do i = 1, nnode
        read(10,*) j, coords(1,i), coords(2,j), coords(3,i)
      end do
      close(unit=10)
c
      open(unit=10,file='../raw_incid.dat')
      read(10,*) nelem
      allocate( incid(8,nelem) ) 
      do i = 1, nelem
        read(10,*) j, incid(1:8,j)
      end do
      close(unit=10)
      write(out,*) '>> coords & incids read ....'
c
c              find x-bar, y-bar, z-bar of every element
c
      allocate( xbar(nelem), ybar(nelem), zbar(nelem),
     &          max_x(nelem), max_y(nelem), max_z(nelem) )
c
      do i = 1, nelem
        xbar(i) = zero
        ybar(i) = zero
        zbar(i) = zero
        max_x(i) = -large
        max_y(i) = -large
        max_z(i) = -large
        do lnode = 1, 8
          snode = incid(lnode,i)
          xbar(i) = xbar(i) + coords(1,snode) * 0.125d0
          ybar(i) = ybar(i) + coords(2,snode) * 0.125d0
          zbar(i) = zbar(i) + coords(3,snode) * 0.125d0
          max_x(i) = max( coords(1,snode), max_x(i) )
          max_y(i) = max( coords(2,snode), max_y(i) )
          max_z(i) = max( coords(3,snode), max_z(i) )
        end do
      end do
      write(out,*) '>> element center coords computed ....'
c
      allocate( stresses(num_sig_values,nelem), 
     &          states(num_states_values,nelem) )
c
c              loop over all possible steps. process steps
c              with results files. find max crack extension for step

      do step = 1, 10000
c
       write(fname1,9000) step
!       write(fname2,9002) step
       inquire(file=fname1,exist=test)
       if( .not. test ) cycle
!       inquire(file=fname2,exist=test) ! both results files must exist
!       if( .not. test ) cycle
c
       open( unit=10, file=fname1, status='old', access="stream",
     &         form="unformatted" )
!       open( unit=11, file=fname2, status='old', access="stream",
!     &         form="unformatted" )
       if( debug ) write(out,*) '      ... processing step: ', step
c
       do i = 1, nelem
        read(10) stresses(1:num_sig_values,i)
       end do
!       do i = 1, nelem
!        read(11) states(1:num_states_values,i)
!       end do
c
       close(unit=10)
!       close(unit=11)
       if( debug ) write(out,*) '         ... stresses read'
c
       max_front = zero
       do i = 1, nelem
        if( ybar(i) > minimum_y_bar ) cycle
        mises = stresses(8,i)
        if( mises > mises_min ) cycle ! killed elements have zero mises
        if( max_x(i) > max_front ) max_front = max_x(i)
       end do ! nelem
c
       crack_extension = max_front - x_0
       if( crack_extension < zero ) crack_extension = zero
       if( max_front >= mesh_limit_xbar ) then
         write(out,9912) step, mesh_limit_xbar
       else
         write(out,9910) step, crack_extension
       end if
c        
      end do ! step
      
  


      stop
 9000 format('wes',i7.7,'_stream')
 9002 format('wem',i7.7,'_stream_smcs')
 9910 format('>> step, max da: ',i5,f10.3)
 9912 format('>> step, extension beyond mesh limit of: ',i5,f10.3)
      end
