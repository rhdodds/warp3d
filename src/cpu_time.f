c     ****************************************************************          
c     *                                                              *          
c     *               functions wcputime, wwalltime                  *          
c     *                      (both are same )                        *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 6/10/2013 rhd              *          
c     *                                                              *          
c     *        returns wall time for process in units of seconds     *          
c     *        previously returned CPU time before all jobs ran      *          
c     *        parallel. now wall time makes more sense.             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      real function wcputime( arg )                                             
      use performance_data, only : start_wall_time                              
      implicit none                                                             
      integer :: arg                                                            
      double precision :: pt, omp_get_wtime, et                                 
c                                                                               
      pt = omp_get_wtime()                                                      
      et = pt - start_wall_time                                                 
      wcputime = real(et)                                                       
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
      real function wwalltime( arg )                                            
      use performance_data, only : start_wall_time                              
      implicit none                                                             
      integer :: arg                                                            
      double precision :: pt, omp_get_wtime, et                                 
c                                                                               
      pt = omp_get_wtime()                                                      
      et = pt - start_wall_time                                                 
      wwalltime = real(et)                                                      
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      function setstarttime                   *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified : 09/28/2011                 *          
c     *                                                              *          
c     *        initialize timing system with the starting wall time  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine setstarttime                                                   
      use performance_data, only : start_wall_time                              
      double precision, external:: omp_get_wtime                                
c                                                                               
      start_wall_time =   omp_get_wtime()                                       
c                                                                               
      return                                                                    
      end                                                                       
