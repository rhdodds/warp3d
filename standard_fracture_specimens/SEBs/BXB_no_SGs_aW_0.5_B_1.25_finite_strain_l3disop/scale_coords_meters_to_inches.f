      implicit double precision (a-h,o-z)
c
c         convert meters to inches
c
c
      scale = 39.37007874d00
      do
       read(*,*,end=100) node, x, y, z
       write(*,1000) node, x*scale, y*scale, z*scale
      end do
c
 100  stop
 1000 format(i8,3f17.10)
      end
