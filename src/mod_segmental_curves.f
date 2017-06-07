c     ****************************************************************          
c     *                                                              *          
c     *              f-90 module segmental_curves                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                    last modified : 5/20/2017 rhd             *          
c     *                                                              *          
c     *     define the variables and data structures to support      *          
c     *     segmental stress-strain curves                           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      module segmental_curves                                                   
c                                                                               
      parameter (max_seg_points=20, max_seg_curves=20,                          
     &           max_seg_curve_sets=10)                                         
c                                                                               
c          This module has two sets of variables. The first                     
c          set down to the next comments are used to store the                  
c          data from user input to define various segmental                     
c          stress-strain curves. These variables are saved and                  
c          restored during setup for an execution of restart.                   
c          None of these variables are changed during                           
c          problem solution once the input data has                             
c          been verified. they are shared read-only across threads              
c                                                                               
c                     double precision/reals                                    
c                                                                               
!dir$ attributes align: 64 :: seg_curves, seg_curves_min_stress,                
     & seg_curves_value, seg_curves_ym, seg_curves_nu,                          
     & seg_curves_alpha, seg_curves_gp_sigma_0,                                 
     & seg_curves_gp_h_u,seg_curves_gp_beta_u,                                  
     & seg_curves_gp_delta_u                                                    
c                                                                               
       double precision ::                                                      
     & seg_curves(max_seg_points,2,max_seg_curves),                             
     & seg_curves_min_stress(max_seg_curves),                                   
     & seg_curves_value(max_seg_curves),                                        
     & seg_curves_ym(max_seg_curves),                                           
     & seg_curves_nu(max_seg_curves),                                           
     & seg_curves_alpha(max_seg_curves),                                        
     & seg_curves_gp_sigma_0(max_seg_curves),                                   
     & seg_curves_gp_h_u(max_seg_curves),                                       
     & seg_curves_gp_beta_u(max_seg_curves),                                    
     & seg_curves_gp_delta_u(max_seg_curves)                                    
c                                                                               
c                     integers                                                  
c                                                                               
!dir$ attributes align: 64 :: num_seg_points, seg_curves_type,                  
     &   seg_curve_table                                                        
c                                                                               
       integer ::                                                               
     &   num_seg_points(max_seg_curves),                                        
     &   seg_curves_type(max_seg_curves),                                       
     &   max_current_pts,                                                       
     &   max_current_curves, num_points, num_curve, num_seg_curve_sets,         
     &   seg_curve_table(max_seg_curves+1,max_seg_curve_sets)                   
c                                                                               
c                     logicals                                                  
c                                                                               
      logical ::                                                                
     &  seg_curve_def(max_seg_curves)                                           
c                                                                               
c          These variable below are used during problem                         
c          solution to support block-by-block computations.                     
c          The module is used to reduced the large number                       
c          of dummy arguments that would be passed many levels                  
c          down thru block computation rouitnes.                                
c                                                                               
c          For OMP threaded parallel, these variables must be                   
c          private for each thread since different threads could                
c          write concurrently on a single, global instance of this              
c          module. We use the threadprivate declaration so                      
c          each thread gets their own copy. This works even for                 
c          allocated arrays                                                     
c                                                                               
       double precision, dimension (:), allocatable ::                          
     &  sigma_curve_min_values,                                                 
     &  curve_temps, curve_e_values, curve_nu_values,                           
     &  curve_alpha_values, curve_rates,                                        
     &  curve_gp_sig_0_values, curve_gp_h_u_values,                             
     &  curve_gp_beta_u_values, curve_gp_delta_u_values                         
                                                                                
!dir$ attributes align: 64 :: curve_plseps_values                               
        double precision, dimension(max_seg_points) ::                          
     &   curve_plseps_values                                                    
c                                                                               
       double precision, dimension (:,:), allocatable ::                        
     &    sigma_curves, sigma_inter_table                                       
c                                                                               
      integer :: active_curve_set=0, now_blk_relem=0                            
      data  curve_plseps_values(1:max_seg_points)/max_seg_points*0.0/           
c                                                                               
c$OMP THREADPRIVATE(                                                            
c$OMP&  sigma_curve_min_values, curve_temps, curve_e_values,                    
c$OMP&  curve_nu_values, curve_alpha_values, curve_rates,                       
c$OMP&  sigma_curves, sigma_inter_table, active_curve_set,                      
c$OMP&  now_blk_relem, curve_plseps_values, curve_gp_sig_0_values,              
c$OMP&  curve_gp_h_u_values, curve_gp_beta_u_values,                            
c$OMP&  curve_gp_delta_u_values )                                               
c                                                                               
c                                                                               
      end module                                                                
                                                                                
                                                                                
