c     ****************************************************************
c     *                                                              *
c     *                      subroutine princ_inv_strain             *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 10/04/94                   *
c     *                                   03/05/95 kck               *        
c     *                                   06/11/97 rhd               *        
c     *                                                              *
c     *     this subroutine computes principal invariants            *
c     *     at the nodes/elem after the primary values have          *
c     *     been averaged                                            * 
c     *                                                              *
c     ****************************************************************
c
      subroutine princ_inv_strain( results, nrowd, num )
      implicit integer (a-z)
#dbl      double precision
#sgl      real
     &     results(nrowd,*), half, t1, t2, t3
      data  half  / 0.5  / 
c   
       do i = 1, num
         t1 =  half * results(i,4)
         t2 =  half * results(i,5)
         t3 =  half * results(i,6)
         results(i,8) = results(i,1) +  results(i,2) +
     &                        results(i,3)
         results(i,9) = - t1 * t1 - t2 * t2 - t3 * t3 +
     &                        results(i,1) * results(i,2) +
     &                        results(i,2) * results(i,3) +
     &                        results(i,1) * results(i,3)
         results(i,10) =
     &     results(i,1) * ( results(i,2) * results(i,3) -
     &              t2 * t2 ) -  t1 * ( t1 * results(i,3) -
     &              t2 * t3 ) + t3 * ( t1 * t2 -
     &              results(i,2) * t3 ) 
       end do
c
       return
       end
