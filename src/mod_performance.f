c     ****************************************************************          
c     *                                                              *          
c     *                      module performance_data                 *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified : 04/27/2017 rhd             *          
c     *                                                              *          
c     *                      stores various profiling data           *          
c     *                                                              *          
c     ****************************************************************          
                                                                                
      module performance_data                                                   
      implicit none                                                             
c                                                                               
      double precision, save :: start_wall_time                                 
      logical :: time_assembly                                                  
      double precision, save :: start_assembly_step, assembly_total             
      integer, save :: ntimes_assembly                                          
c                                                                               
      real, save, private :: time_pardiso, time_warp,                           
     &                       start_run_pardiso                                  
c                                                                               
      contains                                                                  
c                                                                               
        subroutine t_start_assembly( tstart )                                   
        implicit none                                                           
        double precision :: omp_get_wtime, tstart                               
c                                                                               
        tstart = omp_get_wtime()                                                
        return                                                                  
c                                                                               
        end subroutine                                                          
c                                                                               
        subroutine t_end_assembly( tstore, tstart )                             
        implicit none                                                           
        double precision :: omp_get_wtime, tstore, tstart                       
c                                                                               
        tstore = tstore + (omp_get_wtime()-tstart)                              
c                                                                               
        return                                                                  
        end subroutine                                                          
                                                                                
        subroutine t_init_performance                                           
        implicit none                                                           
c                                                                               
        time_pardiso = 0.0                                                      
        time_warp = 0.0                                                         
c                                                                               
        return                                                                  
        end subroutine                                                          
c                                                                               
        subroutine t_performance_eoj( t1 )                                      
        implicit none                                                           
        real :: t1, t2                                                          
c                                                                               
        call cpu_time( t2 )                                                     
        t1 = t2 - time_warp                                                     
c                                                                               
        return                                                                  
        end subroutine                                                          
                                                                                
        subroutine t_performance_start_pardiso                                  
        implicit none                                                           
c                                                                               
        call cpu_time( start_run_pardiso )                                      
c                                                                               
        return                                                                  
        end subroutine                                                          
                                                                                
        subroutine t_performance_end_pardiso                                    
        implicit none                                                           
        real :: t1                                                              
c                                                                               
        call cpu_time( t1 )                                                     
        time_pardiso = time_pardiso + (t1 - start_run_pardiso)                  
c                                                                               
        return                                                                  
        end subroutine                                                          
                                                                                
        subroutine t_performance_eoj_pardiso( t1 )                              
        implicit none                                                           
        real :: t1                                                              
c                                                                               
        t1 = time_pardiso                                                       
c                                                                               
        return                                                                  
        end subroutine                                                          
                                                                                
                                                                                
                                                                                
      end module performance_data                                               
