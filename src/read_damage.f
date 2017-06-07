c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine read_damage                  *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 10/16/94                   *          
c     *                   last modified : 07/28/95 (rhd)             *          
c     *                                                              *          
c     *     this subroutine reads in the damage routines from        *          
c     *     a file.                                                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine read_damage( status, fileno, prec_fact )                       
      use global_data ! old common.main
c                                                                               
      use elem_extinct_data                                                     
      use node_release_data                                                     
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      go to (100,200,300,400,500,600,700,800,900,1000), status                  
c                                                                               
c                               read dam_state, dam_ifv, dam_blk_killed         
c                                                                               
 100  continue                                                                  
      call rd2d( fileno, dam_ifv, prec_fact*mxedof, prec_fact*mxedof,           
     &           num_kill_elem )                                                
      call rdbk( fileno, dam_state, num_kill_elem )                             
      call rdbk( fileno, dam_blk_killed, nelblk )                               
      go to 9999                                                                
c                                                                               
c                               read dam_print list                             
c                                                                               
 200  continue                                                                  
      call rdbk( fileno, dam_print_list, num_print_list )                       
      call rdbk( fileno, old_mises, num_print_list*prec_fact )                  
      call rdbk( fileno, old_mean, num_print_list*prec_fact )                   
      go to 9999                                                                
c                                                                               
c                               read kill_order_list                            
c                                                                               
 300  continue                                                                  
      call rdbk( fileno, kill_order_list, num_kill_order_list )                 
      go to 9999                                                                
c                                                                               
c                               read dam_node_elecnt                            
c                                                                               
 400  continue                                                                  
      call rdbk( fileno, dam_node_elecnt, nonode )                              
      go to 9999                                                                
c                                                                               
c                               read dam_face_nodes,                            
c                               dam_dbar_elems                                  
c                                                                               
 500  continue                                                                  
      call rdbk( fileno, dam_face_nodes, 4*num_kill_elem )                      
      call rdbk( fileno, dam_dbar_elems, (2*prec_fact)*num_kill_elem )          
      go to 9999                                                                
c                                                                               
c                               read all of the node release crack              
c                               growth variables                                
c                                                                               
 600  continue                                                                  
      call rdbk( fileno, crack_plane_nodes, num_crack_plane_nodes )             
      call rdbk( fileno, inv_crkpln_nodes, nonode )                             
      call rdbk( fileno, num_neighbors, num_crack_plane_nodes )                 
      call rd2d( fileno, neighbor_nodes, mxconn, mxconn,                        
     &           num_crack_plane_nodes )                                        
      call rdbk( fileno, crkpln_nodes_state, num_crack_plane_nodes )            
      call rdbk( fileno, crkpln_nodes_react,                                    
     &           num_crack_plane_nodes * prec_fact )                            
      call rd2d( fileno, crack_front_nodes, num_crack_plane_nodes,              
     &           num_crack_plane_nodes,2 )                                      
      go to 9999                                                                
c                                                                               
c                               read the info for traction separation           
c                                                                               
 700  continue                                                                  
      call rdbk( fileno, node_release_frac,                                     
     &           num_crack_plane_nodes * prec_fact )                            
      go to 9999                                                                
c                                                                               
c                               read the old angles for overshoot control       
c                                                                               
 800  continue                                                                  
      call rd2d( fileno, old_angles_at_front, num_crack_plane_nodes             
     &       * prec_fact, num_crack_plane_nodes * prec_fact, mxconn )           
      go to 9999                                                                
c                                                                               
c                               read the old damage values needed for           
c                               load size control during crack growth           
c                                                                               
 900  continue                                                                  
      if( crack_growth_type .eq. 1 ) then                                       
         call rdbk( fileno, old_porosity, num_kill_elem * prec_fact )           
         call rdbk( fileno, del_poros, mxstp_store * prec_fact )                
      else if( crack_growth_type .eq. 3 ) then                                  
         call rdbk( fileno, old_plast_strain, num_kill_elem *                   
     &              prec_fact )                                                 
      else if( crack_growth_type .eq. 4 ) then                                  
         call rdbk( fileno, old_deff, num_kill_elem * prec_fact )               
         call rdbk( fileno, del_deff, mxstp_store * prec_fact )                 
      end if                                                                    
      go to 9999                                                                
c                                                                               
c                               read the data for constant front growth         
c                                                                               
 1000 continue                                                                  
      call rdbk( fileno, master_nodes, num_crack_fronts )                       
      call rd2d( fileno, crack_front_list,                                      
     &           num_crack_fronts*num_nodes_grwinc,                             
     &           num_crack_fronts*num_nodes_grwinc, num_nodes_thick )           
      call rd2d( fileno, master_lines, num_crack_fronts,                        
     &           num_crack_fronts, num_nodes_back + 1 )                         
      go to 9999                                                                
c                                                                               
c                                                                               
 9999 continue                                                                  
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
