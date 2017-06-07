c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine equiv_strain                 *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 10/04/94                   *          
c     *                                   06/02/95 kck               *          
c     *                                                              *          
c     *     this subroutine computes mises yield function value      *          
c     *     at the nodes/elem after the primary values have          *          
c     *     been averaged                                            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine equiv_strain ( results, maxnum, num )                          
      implicit integer (a-z)                                                    
      double precision                                                          
     &     results(maxnum,*),  root23, onep5                                    
      data root23, onep5 / 0.471404, 1.5 /                                      
c                                                                               
      do i = 1, num                                                             
        results(i,7) = root23 * sqrt(                                           
     &     ( results(i,1) - results(i,2) ) ** 2 +                               
     &     ( results(i,2) - results(i,3) ) ** 2 +                               
     &     ( results(i,1) - results(i,3) ) ** 2 +                               
     &   onep5 * ( results(i,4)**2 + results(i,5)**2+                           
     &         results(i,6)**2 ) )                                              
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
