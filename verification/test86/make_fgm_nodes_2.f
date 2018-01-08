      implicit none

      integer :: i, inode, nnodes

      double precision :: coords(3,100000), v1(3), v2(3), v3(3)
      double precision :: e, yld_pt, n_power
      character(len=72) :: line



      nnodes = 45603

      open(unit=10,file="coordinates.inp")
      read(10,fmt="(a)") line
      read(10,fmt="(a)") line
      read(10,fmt="(a)") line
      do i = 1, nnodes
         read(10,*) inode, coords(1:3,inode)
      end do
      close(10)
c
      v2(1) = -0.5d0
      v2(2) = 0.5d0
      v2(3) = 0.0d0
c
      do i = 1, nnodes
        v1(1) = coords(1,i) -0.5d0
        v1(2) = coords(2,i)
        v1(3) = 0.0d0
        call cross_prod( v1, v2, v3 )
        if( v3(3) < 0.0d0 ) then
                e = 30000.0d0
                yld_pt = 60.0d0
                n_power = 10.0d0
        else
                e = 10000.0d0
                yld_pt = 30.0d0
                n_power = 5.0d0
        end if
        write(*,9000) i,e, yld_pt, n_power
      end do

 9000 format(2x,i6,2x,'e',f10.1,' yld_pt ',f5.1,' n_power ',f5.1)
      stop
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine cross_prod                   *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 07/02/98                   *
c     *                                                              *
c     *        This routine computes the cross product of the two    *
c     *        input vectors.                                        *
c     *                                                              *
c     ****************************************************************
c
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
