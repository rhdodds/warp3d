c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouyld1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 06/19/02 rhd               *
c     *                                                              *
c     *     this subroutine computes the mises equiv. stress at a    *
c     *     given strain point for a block of elements               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouyld1( span, gpstr, yf, mxvl )
      implicit integer (a-z)
      double precision
     &     gpstr(mxvl,*), yf(*), iroot2, six
      data iroot2, six / 0.70711, 6.0 /
c
c                       compute the von-mises stress.
c
      do i = 1, span
         yf(i) = sqrt( (gpstr(i,1)-gpstr(i,2))**2+
     &                 (gpstr(i,2)-gpstr(i,3))**2+
     &                 (gpstr(i,1)-gpstr(i,3))**2+
     &             six*(gpstr(i,4)**2+gpstr(i,5)**2+
     &                  gpstr(i,6)**2) )*iroot2
      end do
c
      return
      end




