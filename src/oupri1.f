c     ****************************************************************
c     *                                                              *
c     *                      subroutine oupri1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 07/04/91                   *
c     *                                                              *
c     *     this subroutine computes the principal cauchy stresses   *
c     *     or almansi strains and the direction cosines of their    *
c     *     corresponding normals at a strain point, either a gauss  *
c     *     point or a node point, for a block of similar, non-      *
c     *     conflicting q3disop or l3disop elements.                 *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine oupri1( span, str, prstr, angles, stress )
      implicit integer (a-z)
$add param_def
#dbl      double precision
#sgl      real
     &     str(mxvl,*), prstr(mxvl,*), angles(mxvl,ndim,*)
      logical stress
c
c                    locally allocated
c
#dbl      double precision
#sgl      real
     &     tstr(nstr,mxvl), wk(ndim,mxvl), ev(nstr,mxvl),
     &     evec(ndim,ndim,mxvl), factor, one, half
      data one, half / 1.0, 0.5 /
c
      factor = one
      if ( .not. stress ) factor = half
      do i = 1, span
         tstr(1,i) = str(i,1)
         tstr(2,i) = str(i,4) * factor
         tstr(3,i) = str(i,2)
         tstr(4,i) = str(i,6) * factor
         tstr(5,i) = str(i,5) * factor
         tstr(6,i) = str(i,3) 
      end do
c             
      do i = 1, span
         call ou3dpr( tstr(1,i), ndim, 1, ev(1,i), evec(1,1,i),
     &                ndim, wk(1,i), ier )
      end do
c                               
c
      do i = 1, span
         prstr(i,1)    = ev(1,i)
         prstr(i,2)    = ev(2,i)
         prstr(i,3)    = ev(3,i)
         angles(i,1,1) = evec(1,1,i)
         angles(i,1,2) = evec(1,2,i)
         angles(i,1,3) = evec(1,3,i)
         angles(i,2,1) = evec(2,1,i)
         angles(i,2,2) = evec(2,2,i)
         angles(i,2,3) = evec(2,3,i)
         angles(i,3,1) = evec(3,1,i)
         angles(i,3,2) = evec(3,2,i)
         angles(i,3,3) = evec(3,3,i)
      end do
c
c
      return
      end


