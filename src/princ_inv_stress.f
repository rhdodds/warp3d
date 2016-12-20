c     ****************************************************************
c     *                                                              *
c     *                      subroutine princ_inv_stress             *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 10/04/94                   *
c     *                                   03/05/95 kck               *        
c     *                                   06/15/97 rhd               *        
c     *                                                              *
c     *     this subroutine computes principal invariants            *
c     *     at the nodes/elements after the primary values have      * 
c     *     been averaged                                            * 
c     *                                                              *
c     ****************************************************************
c
      subroutine princ_inv_stress( results, mxvl, num )
      implicit integer (a-z)
      double precision
     &     results(mxvl,*)
        
       do i = 1, num
         results(i,12) = results(i,1) +  results(i,2) +
     &                        results(i,3)
         results(i,13) = results(i,1) * results(i,2) + 
     &                        results(i,2) * results(i,3) +
     &                        results(i,1) * results(i,3) -
     &                        results(i,4) * results(i,4) -
     &                        results(i,5) * results(i,5) -
     &                        results(i,6) * results(i,6)
         results(i,14) = results(i,1) * 
     &            ( results(i,2) * results(i,3) -
     &              results(i,5) * results(i,5) ) -
     &                        results(i,4) *
     &            ( results(i,4) * results(i,3) -
     &              results(i,5) * results(i,6) ) +
     &                        results(i,6) *
     &            ( results(i,4) * results(i,5) -
     &              results(i,2) * results(i,6) ) 
       end do
c
       return
       end








