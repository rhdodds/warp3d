c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine thyme                        *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 06/15/90                   *          
c     *                                                              *          
c     *     this subroutine handles the timing of various routines.  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine thyme( calc, flag )    
c                                              
      use global_data, only : times, strtm
c      
      implicit none
c      
      integer :: calc, flag      
      real ::  t1
      real, external :: wcputime                                                        
c
      if ( flag .eq. 1 ) then  ! initial call                                                                               
         strtm = wcputime ( 1 )                                                 
      else  ! accumulate wall time and counts.
         t1            = wcputime ( 1 )                                         
         times(calc,1) = times(calc,1) + t1 - strtm                             
         times(calc,2) = times(calc,2) + 1.0                                    
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
