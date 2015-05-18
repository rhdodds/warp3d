c     ****************************************************************
c     *                                                              *
c     *                      subroutine oueff1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 04/18/90                   *
c     *                                                              *
c     *     this subroutine computes the effective strain measure    *
c     *     for a block of similar, non-conflicting q3disop or       *
c     *     l3disop elements.                                        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine oueff1( span, str, efeps, mxvl )
      implicit integer (a-z)
#dbl      double precision
#sgl      real
     &     str(mxvl,*), efeps(*), root23, onep5
      data root23, onep5 / 0.471404, 1.5 /
c
c                       compute the effective strain measure.
c
      do i = 1, span
       efeps(i) = root23 *
     &      sqrt( (str(i,1)-str(i,2))**2+(str(i,2)-str(i,3))**2+
     &            (str(i,1)-str(i,3))**2+
     &            onep5*(str(i,4)**2+str(i,5)**2+str(i,6)**2) )
      end do
c
      return
      end 
