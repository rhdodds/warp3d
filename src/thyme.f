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
      use global_data ! old common.main
      implicit integer (a-z)                                                    
      real  t1, wcputime                                                        
      external wcputime                                                         
c                                                                               
c                                                                               
      if ( slave_processor ) return                                             
c                                                                               
      if ( flag .eq. 1 ) then                                                   
c                                                                               
c                       initial call.                                           
c                                                                               
         strtm = wcputime ( 1 )                                                 
c                                                                               
      else                                                                      
c                                                                               
c                       final call                                              
c                                                                               
         t1            = wcputime ( 1 )                                         
         times(calc,1) = times(calc,1) + t1 - strtm                             
         times(calc,2) = times(calc,2) + 1.0                                    
c                                                                               
      end if                                                                    
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
