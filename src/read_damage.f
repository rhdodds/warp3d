c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine read_damage                  *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 4/8/26 rhd                 *          
c     *                                                              *          
c     *              reads damage data from restart file             * 
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                             
      subroutine read_damage( action, fileno ) 
c                      
      use global_data, only : out, mxedof, nonode, mxconn
c                                                                               
      use elem_extinct_data                                                     
      use node_release_data                                                     
      use damage_data                                                           
c                                                                               
      implicit none
c
      integer, intent(in) :: action, fileno
c
      integer :: isize, count, i, j, nrow_ek, np
      logical :: std_kill                                                  
c                                                                               
      std_kill = .not. use_mesh_regularization
c
      select case( action )
c
c              read dam_state, dam_ifv
c                                                                               
      case( 1 )
      call rdbk_i4( fileno, dam_state, num_kill_elem )                             
      call rd2d_r8( fileno, dam_ifv, mxedof, mxedof, num_kill_elem )                                                
c                                                                               
c              read dam_print list                             
c                                                                               
      case( 2 )
      call rdbk_i4( fileno, dam_print_list, num_print_list )                       
c                                                                               
c              read kill_order_list                            
c                                                                               
      case( 3 )
      call rdbk_i4( fileno, kill_order_list, num_kill_order_list )                 
c                                                                               
c              read dam_node_elecnt                            
c                                                                               
      case( 4 )
      call rdbk_i4( fileno, dam_node_elecnt, nonode )                              
c                                                                               
c              read dam_face_nodes, dam_dbar_elems                                  
c                                                                               
      case( 5 )
      call rdbk_i4( fileno, dam_face_nodes, 4*num_kill_elem )                      
      call rdbk_r8( fileno, dam_dbar_elems, 2*num_kill_elem )          
c                                                                               
c              read all of the node release crack growth variables                                
c                                                                               
      case( 6 )
      call rdbk_i4( fileno, crack_plane_nodes, num_crack_plane_nodes )             
      call rdbk_i4( fileno, inv_crkpln_nodes, nonode )                             
      call rdbk_i4( fileno, num_neighbors, num_crack_plane_nodes )                 
      call rd2d_i4( fileno, neighbor_nodes, mxconn, mxconn,                        
     &           num_crack_plane_nodes )                                        
      call rdbk_i4( fileno, crkpln_nodes_state, num_crack_plane_nodes )            
      call rdbk_r8( fileno, crkpln_nodes_react,                                    
     &              num_crack_plane_nodes )                            
      call rd2d_i4( fileno, crack_front_nodes, num_crack_plane_nodes,              
     &              num_crack_plane_nodes, 2 )                                      
c                                                                               
c              read the info for traction separation           
c                                                                               
      case( 7 )
      call rdbk_r8( fileno, node_release_frac,                                     
     &           num_crack_plane_nodes )                            
c                                                                               
c              read the old angles for overshoot control       
c                                                                               
      case( 8 ) 
      call rd2d_r8( fileno, old_angles_at_front, num_crack_plane_nodes,             
     &              num_crack_plane_nodes, mxconn )           
c                                                                               
c              read the old damage values needed for           
c              load size control during crack growth           
c                                                                               
      case( 9 )
      if( crack_growth_type .eq. 1 ) then                                       
         call rdbk_r8( fileno, gt_old_porosity, num_kill_elem )           
         call rdbk_r8( fileno, del_poros, mxstp_store )                
      else if( crack_growth_type .eq. 3 ) then                                  
         call rdbk_r8( fileno, old_plast_strain, num_kill_elem )                                                 
      else if( crack_growth_type .eq. 4 ) then                                  
         call rdbk_r8( fileno, cohes_old_deff, num_kill_elem )               
         call rdbk_r8( fileno, del_deff, mxstp_store )                 
      end if                                                                    
c                                                                               
c              read the data for constant front growth         
c                                                                               
      case( 10 )
      call rdbk_i4( fileno, master_nodes, num_crack_fronts )                       
      call rd2d_i4( fileno, crack_front_list,                                      
     &              num_crack_fronts*num_nodes_grwinc,                             
     &              num_crack_fronts*num_nodes_grwinc, 
     &              num_nodes_thick )           
      call rd2d_i4( fileno, master_lines, num_crack_fronts,                        
     &              num_crack_fronts, num_nodes_back + 1 )                         
c                                                                               
c              read the data smcs         
c                                                                               
      case( 11 )
      call rdbk_r8( fileno, smcs_old_epsplas, num_kill_elem ) 
      call rdbk_r8( fileno, smcs_tear_param, num_kill_elem )           
      read(fileno) isize
      if( isize > 0 ) then
          allocate(smcs_states_intlst(isize)) 
          read(fileno)  smcs_states_intlst(1:isize)   
      end if 
      if( use_weighted ) then
        isize = num_kill_elem
        call rdbk_r8( fileno, smcs_weighted_T, isize )
        call rdbk_r8( fileno, smcs_weighted_zeta, isize )
        call rdbk_r8( fileno, smcs_weighted_bar_theta, isize )
      end if  
c                                                                               
c              read the mesh regularization data
c                                                                               
      case( 12 )
      call rdbk_r8( fileno, smcs_d_values, num_kill_elem )           
      call rdbk_r8( fileno, smcs_eps_plas_at_death, num_kill_elem )    
      call rdbk_r8( fileno, smcs_stress_at_death, num_kill_elem )
      call rdbk_i4( fileno, smcs_start_kill_step, num_kill_elem )   
c
c              read the Oddy distortion metrics
c                                                                               
      case( 13 )
      if( .not. use_distortion_metric ) return
      if( .not. allocated( Oddy_metrics ) ) 
     &        allocate( Oddy_metrics(num_kill_elem,2) )
      read(fileno) Oddy_metrics
c
      case default
         write(out,9000)
         call die_abort
c
      end select
c                                                                               
      return    
c
 9000 format('>> FATAL ERROR: routine read_damage'                                     
     &  /,   '                job terminated' )      
 9100 format(10x,8f5.2)                          
c
      end subroutine read_damage                                                                       
                                                                                
                                                                                
