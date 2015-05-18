c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouinv1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 07/16/89                   *
c     *                                                              *
c     *     this subroutine computes the invariants of the cauchy    *
c     *     stress or the almansi strain at a strain point, either a *
c     *     gauss point or a node point, for a block of similar,     *
c     *     non-conflicting q3disop or l3disop elements.             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouinv1( span, str, inv, stress )
      implicit integer (a-z)
$add param_def
#dbl      double precision
#sgl      real
     &     str(mxvl,*), inv(mxvl,*)
      logical stress
c
c                    locally allocated
c
#dbl      double precision
#sgl      real
     &     tstr(mxvl,nstr), factor, one, half
      data one, half / 1.0, 0.5 /
c
      factor = one
      if ( .not. stress ) factor = half
      do i = 1, span
         tstr(i,1) = str(i,1)
         tstr(i,2) = str(i,2)
         tstr(i,3) = str(i,3)
         tstr(i,4) = str(i,4) * factor
         tstr(i,5) = str(i,5) * factor
         tstr(i,6) = str(i,6) * factor
      enddo
c
c                       compute the principal invariants.
c
      do i = 1, span                        
       inv(i,1) = tstr(i,1)+tstr(i,2)+tstr(i,3)
       inv(i,2) = tstr(i,4)*tstr(i,4)+tstr(i,5)*tstr(i,5) +
     &            tstr(i,6)*tstr(i,6)-tstr(i,1)*tstr(i,2) -
     &            tstr(i,2)*tstr(i,3)-tstr(i,1)*tstr(i,3)
       inv(i,3)=
     &   tstr(i,1)*( tstr(i,2)*tstr(i,3)-tstr(i,5)*tstr(i,5) ) -
     &   tstr(i,4)*( tstr(i,4)*tstr(i,3)-tstr(i,5)*tstr(i,6) ) +
     &   tstr(i,6)*( tstr(i,4)*tstr(i,5)-tstr(i,2)*tstr(i,6) )
      end do
c
c
      return
      end
