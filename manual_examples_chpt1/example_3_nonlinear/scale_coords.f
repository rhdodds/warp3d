      implicit double precision (a-h,o-z)
c
      scale = 1000.0d0

      read(*,*) nnode, noelem
      do i = 1, nonode
       read(*,*,end=100) x, y, z
       ynew = y
       if( y > 1.0d0 )then 
          ynew = 1.0d0
        end if
       xnew = x ! - 2220.d0
       write(*,1000)  xnew, ynew, z
      end do
c
 100  stop
 1000 format(3e17.9)
      end
