c     ****************************************************************
c     *                                                              *
c     *                      subroutine setelb                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 08/11/2017 rhd             *
c     *                                                              *
c     *     this subroutine sets the element library to include all  *
c     *     of the elements to be implemented by the program.        *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine setelb( nlibel, elelib, outmap, mxlbel, out )
c
      implicit none
c
      integer :: nlibel, outmap(mxlbel,*), mxlbel, out
      character(len=8) :: elelib(mxlbel)
c
      integer :: i
c
      nlibel = 19
      if( nlibel > mxlbel ) then
         write(out,9000)
         call die_abort
      end if
c
      elelib(1:mxlbel) = ' '
c
c                       element type one: q3disop
c
      elelib(1) = 'q3disop '
      do i = 1, 30
        outmap(1,i) = i
      end do
c
c                       element type two: l3disop
c
      elelib(2) = 'l3disop '
      do i = 1, 30
        outmap(2,i) = i
      end do
c
c                       element type 3: ts12isop
c
      elelib(3) = 'ts12isop'
      do i = 1, 30
        outmap(3,i) = i
      end do
c
c                       element type 4: ts15isop
c
      elelib(4) = 'ts15isop'
      do i = 1, 30
        outmap(4,i) = i
      end do
c
c                       element type 5: ts9isop
c
      elelib(5) = 'ts9isop'
      do i = 1, 30
        outmap(5,i) = i
      end do
c
c                       element type 6: tet10
c
      elelib(6) = 'tet10'
      do i = 1, 30
        outmap(6,i) = i
      end do
c
c                       element type 7: wedge15
c
      elelib(7) = 'wedge15'
      do i = 1, 30
        outmap(7,i) = i
      end do
c
c                       element type 8: tri6
c
      elelib(8) = 'tri6'
      do i = 1, 30
        outmap(8,i) = i
      end do
c
c                       element type 9: quad8
c
      elelib(9) = 'quad8'
      do i = 1, 30
        outmap(9,i) = i
      end do
c
c                       element type 10: axiquad8
c
      elelib(10) = 'axiquad8'
      do i = 1, 30
        outmap(10,i) = i
      end do
c
c                       element type 11: axitri6
c
      elelib(11) = 'axitri6'
      do i = 1, 30
        outmap(11,i) = i
      end do
c
c                       element type 12: 'inter_8'
c
      elelib(12) = 'inter_8'
      do i = 1,30
        outmap(12,i) = i
      end do
c
c                       element type 13: tet4
c
      elelib(13) = 'tet4'
      do i = 1, 30
        outmap(13,i) = i
      end do
c
c                       element type 14: 'trint6'
c
      elelib(14) = 'trint6'
      do i = 1,30
        outmap(14,i) = i
      end do
c
c                       element type 15: 'trint12'
c
      elelib(15) = 'trint12'
      do i = 1,30
        outmap(15,i) = i
      end do
c
c                       element type 16: not a real element.
c                                        used to provide integration
c                                        points for 4-node quad
c
      elelib(16) = 'not_real'
      do i = 1,30
        outmap(16,i) = i
      end do
c
c                       element type 17: not a real element.
c                                        used to provide integration
c                                        points for 3 node line.
c                                        used in domain integral work
c
      elelib(17) = 'line3'
      do i = 1,30
        outmap(17,i) = i
      end do
c
c
c                       element type 18: 'bar2'
c
      elelib(18) = 'bar2'
      do i = 1,30
        outmap(18,i) = i
      end do
c
c
c                       element type 19: 'link2'
c
      elelib(19) = 'link2'
      do i = 1,30
        outmap(19,i) = i
      end do
c      
      return
c
 9000 format(/,'>>>> FATAL ERROR: setelb. incorrect (nlibel,mxlebl)',
     & ' values',
     &       /,'                  job terminated',//)
      end
