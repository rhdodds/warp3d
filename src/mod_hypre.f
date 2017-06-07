c **********************************************************************        
c     Module hypre_paramters                                                    
c                                                                               
c     Just stores parameters for the solver.                                    
c                                                                               
c     created: mcm 5/11                                                         
c                                                                               
c **********************************************************************        
      Module hypre_parameters                                                   
            implicit none                                                       
                                                                                
c           1 = parasails - default                                             
c           2 = boomeramg                                                       
            integer :: precond_type                                             
c           1 = pcg - default                                                   
            integer :: hsolver_type                                             
c           per hypre - default = 0                                             
            integer :: precond_printlevel                                       
c           per hypre - default = 0                                             
            integer :: solver_printlevel                                        
                                                                                
c           Solver properties                                                   
            double precision :: hypre_tol                                       
            integer :: hypre_max                                                
                                                                                
c           Parasails properties                                                
c           default = 1                                                         
            integer :: levels                                                   
c           default = 0.1                                                       
            double precision :: threshold                                       
c           default = 0.1                                                       
            double precision :: filter                                          
c           default = 1                                                         
            integer :: symme                                                    
c           default = 0                                                         
            double precision :: loadbal                                         
c                                                                               
c           BoomerAMG properties                                                
c                                                                               
c           maximum number of mg levels                                         
            integer :: max_levels                                               
c           multi grid drop threshold                                           
            double precision :: mg_threshold                                    
c           coarsening type                                                     
            integer :: coarsening                                               
c           number of levels of aggressive coarsening                           
            integer :: agg_levels                                               
c           interpolation type                                                  
            integer :: interpolation                                            
c           truncation tolerance                                                
            double precision :: truncation                                      
c           Relaxation type                                                     
            integer :: relaxation                                               
c           Relaxation weight and outer weight                                  
            double precision :: relax_wt, relax_outer_wt                        
c           Number of sweeps (both ways                                         
            integer :: sweeps                                                   
c           CF relaxation or not                                                
            integer :: cf                                                       
c           Cycle type                                                          
            integer :: cycle_type                                               
                                                                                
                                                                                
                                                                                
c           I do not remember what this is for                                  
            integer :: error_count                                              
c           Count the number of preconditioning errors.  If 2 in a a row, give u
            integer :: precond_fail_count                                       
c           This logical determines if the manager!!! needs to call an adaptive 
            logical :: hyp_trigger_step                                         
                                                                                
c           Important - need for allocation/waking/sleeping procs               
            logical :: hyp_first_solve                                          
                                                                                
      end module hypre_parameters                                               
