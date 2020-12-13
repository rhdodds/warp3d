c     ****************************************************************          
c     *                                                              *          
c     *                   f-90 module damage_data                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *              last modified : 12/12/20 rhd                    *          
c     *                                                              *          
c     *     define the variables and data structures to support      *          
c     *     crack growth using damage parameters (e.g., the Gurson   *          
c     *     model, ctoa, etc.)                                       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      module damage_data                                                        
c                                                                               
      parameter (mxstp_store = 10)                                              
c                                                                               
c                     arrays                                                    
c                                                                               
      integer, save, allocatable :: dam_ptr(:),
     &          smcs_states_intlst(:), deleted_elist_to_stop(:)   
c                            
      double precision ::                                                         
     &   del_poros(mxstp_store), del_deff(mxstp_store)   
c
      integer :: user_kill_list_now(100)                      
c                                                                               
c                     scalar double precision/reals                             
c                                                                               
       double precision ::                                                         
     &   porosity_limit, gurson_cell_size,                                      
     &   crack_plane_coord, release_fraction,                                   
     &   critical_angle, release_height,                                        
     &   crack_plane_sign, char_length,                                         
     &   init_crit_ang, 
     &   smcs_gamma, smcs_beta_1, smcs_beta_2, smcs_a_plus,
     &   smcs_a_minus, smcs_kappa, smcs_alpha_1, smcs_alpha_2, 
     &   smcs_cutoff_triaxiality, max_eps_critical, smcs_type_4_A,
     &   smcs_type_4_n,
     &   smcs_type_4_c1, smcs_type_4_c2, smcs_type_4_c3,
     &   smcs_alpha, smcs_beta,                               
     &   control_load_fact, old_load_fact,                                      
     &   min_load_fact, overshoot_limit, CTOA_range,                            
     &   perm_load_fact, max_porosity_change,                                   
     &   max_plast_strain_change,                                               
     &   init_ctoa_dist, ctoa_dist,                                             
     &   crkpln_srch_tol, max_deff_change,                                      
     &   critical_cohes_deff_fract,                                             
     &   ppr_kill_displ_fraction
                                               
c                                                                               
c                     scalar integers                                           
c                                                                               
      integer :: crack_growth_type,                                                
     &   num_kill_elem, max_dam_state, csttail, num_print_list,                 
     &   num_kill_order_list, release_type, smcs_type,                                     
     &   crk_pln_normal_idx, num_crack_plane_nodes, crack_front_start,          
     &   crack_front_end, crkfrnt_garbage_start, crkfrnt_garbage_end,           
     &   min_steps_for_release, num_nodes_thick, num_crack_fronts,              
     &   num_nodes_back, num_nodes_grwinc, num_steps_min,                       
     &   num_elements_killed, stop_killed_elist_length,                                                  
     &   num_elements_in_force_release, num_ctoa_released_nodes,
     &   num_user_kill_elems, killed_element_limit, num_top_list               
c                                                                               
c                     scalar logicals                                           
c                                                                               
      logical :: no_killed_elems, print_status, kill_order,                        
     &  kill_order_now, no_released_nodes, list_crkpln_nodes,                   
     &  list_crkfrnt_nodes, growth_by_kill, growth_by_release,                  
     &  growth_by_cohesive, enforce_node_release,                               
     &  overshoot_control_crk_grth, overshoot_allocated,                        
     &  load_size_control_crk_grth, g_stp_cntrl_allocated,                      
     &  const_front, master_lines_set, load_reduced, all_elems_killed,
     &  print_top_list, smcs_states, smcs_stream, smcs_text,
     &  smcs_list_file_flag           
c                                                                               
c                     character                             
c      
      character(len=80) :: smcs_list_file_name                                                                         
c                                                                               
      end module                                                                
