c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine equiv_strain                 *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 12/14/2018 rhd             *          
c     *                                                              *          
c     *     this subroutine computes mises yield function value      *          
c     *     at the nodes/elem after the primary values have          *          
c     *     been averaged                                            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine equiv_strain( results, maxnum, num )                          
      implicit none  
c
      integer :: maxnum, num                                                  
      double precision :: results(maxnum,*)
c
      integer :: i
      double precision, parameter :: root23=0.471404d0, onep5=1.5d0                                      
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
