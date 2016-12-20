c     ****************************************************************
c     *                                                              *
c     *                      subroutine yield_function               *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 10/04/94                   *
c     *                                   03/05/95 kck               *
c     *                                                              *
c     *     this subroutine computes mises yield function value      *
c     *     at the nodes/elements after the primary values have      *
c     *     been averaged                                            * 
c     *                                                              *
c     ****************************************************************
c
      subroutine yield_function ( results, maxnum, num )
      implicit integer (a-z)
      double precision
     &     results(maxnum,*), iroot2, six
      data iroot2, six / 0.70711, 6.0 /
c
      do i = 1, num
        results(i,8) = sqrt( 
     &     ( results(i,1) - results(i,2) ) ** 2 +
     &     ( results(i,2) - results(i,3) ) ** 2 +
     &     ( results(i,1) - results(i,3) ) ** 2 +
     &   six*( results(i,4)**2 + results(i,5)**2+
     &         results(i,6)**2 ) )*iroot2
c
      end do
c
      return 
      end
