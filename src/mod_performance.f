c     ****************************************************************          
c     *                                                              *          
c     *                      module performance_data                 *          
c     *                                                              *          
c     *                       written by : mcm                       *          
c     *                                                              *          
c     *                   last modified : 3/8/2026 rhd               *          
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
      double precision, external :: omp_get_wtime            
      integer, save :: ntimes_assembly                                          
c                                                                               
      real, save, private :: time_pardiso, time_warp,  time_mumps,                        
     &                       start_run_pardiso,  start_run_mumps                                
c                                                                               
      contains                                                                  
c                                                                               
        subroutine t_start_assembly( tstart )                                   
        implicit none                                                           
        double precision ::  tstart                               
c                                                                               
        tstart = omp_get_wtime()                                                
        return                                                                  
c                                                                               
        end subroutine                                                          
c                                                                               
        subroutine t_end_assembly( tstore, tstart )                             
        implicit none                                                           
        double precision :: tstore, tstart                       
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
        time_mumps = 0.0                                                       
c                                                                               
        return                                                                  
        end subroutine                                                          
c                                                                               
        subroutine t_performance_eoj( t1 )                                      
        implicit none                                                           
        real :: t1                                                     
c                                                                               
        t1 = omp_get_wtime() - time_warp                                                     
c                                                                               
        return                                                                  
        end subroutine                                                          
                                                                                
        subroutine t_performance_start_pardiso                                  
        implicit none                                                           
c                                                                               
        start_run_pardiso =  omp_get_wtime()                                   
c                                                                               
        return                                                                  
        end subroutine                                                          
                                                                                
        subroutine t_performance_end_pardiso                                    
        implicit none                                                           
c
        time_pardiso = time_pardiso + ( omp_get_wtime() - 
     &                 start_run_pardiso )                  
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
c                                                                                
      end module performance_data                                               
