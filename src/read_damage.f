c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine read_damage                  *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 9/5/23 rhd                 *          
c     *                                                              *          
c     *              reads damage data from restart file             * 
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine read_damage( action, fileno, prec_fact ) 
c                      
      use global_data, only : out, mxedof, nonode, mxconn
c                                                                               
      use elem_extinct_data                                                     
      use node_release_data                                                     
      use damage_data                                                           
c                                                                               
      implicit none
c
      integer, intent(in) :: action, fileno, prec_fact
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
      call rdbk( fileno, dam_state, num_kill_elem )                             
      call rd2d( fileno, dam_ifv, prec_fact*mxedof,
     &        prec_fact*mxedof, num_kill_elem )                                                
c                                                                               
c              read dam_print list                             
c                                                                               
      case( 2 )
      call rdbk( fileno, dam_print_list, num_print_list )                       
c                                                                               
c              read kill_order_list                            
c                                                                               
      case( 3 )
      call rdbk( fileno, kill_order_list, num_kill_order_list )                 
c                                                                               
c              read dam_node_elecnt                            
c                                                                               
      case( 4 )
      call rdbk( fileno, dam_node_elecnt, nonode )                              
c                                                                               
c              read dam_face_nodes, dam_dbar_elems                                  
c                                                                               
      case( 5 )
      call rdbk( fileno, dam_face_nodes, 4*num_kill_elem )                      
      call rdbk( fileno, dam_dbar_elems, (2*prec_fact)*num_kill_elem )          
c                                                                               
c              read all of the node release crack growth variables                                
c                                                                               
      case( 6 )
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
c                                                                               
c              read the info for traction separation           
c                                                                               
      case( 7 )
      call rdbk( fileno, node_release_frac,                                     
     &           num_crack_plane_nodes * prec_fact )                            
c                                                                               
c              read the old angles for overshoot control       
c                                                                               
      case( 8 ) 
      call rd2d( fileno, old_angles_at_front, num_crack_plane_nodes             
     &       * prec_fact, num_crack_plane_nodes * prec_fact, mxconn )           
c                                                                               
c              read the old damage values needed for           
c              load size control during crack growth           
c                                                                               
      case( 9 )
      if( crack_growth_type .eq. 1 ) then                                       
         call rdbk( fileno, gt_old_porosity, num_kill_elem * prec_fact )           
         call rdbk( fileno, del_poros, mxstp_store * prec_fact )                
      else if( crack_growth_type .eq. 3 ) then                                  
         call rdbk( fileno, old_plast_strain, num_kill_elem *                   
     &              prec_fact )                                                 
      else if( crack_growth_type .eq. 4 ) then                                  
         call rdbk( fileno, cohes_old_deff, num_kill_elem * prec_fact )               
         call rdbk( fileno, del_deff, mxstp_store * prec_fact )                 
      end if                                                                    
c                                                                               
c              read the data for constant front growth         
c                                                                               
      case( 10 )
      call rdbk( fileno, master_nodes, num_crack_fronts )                       
      call rd2d( fileno, crack_front_list,                                      
     &           num_crack_fronts*num_nodes_grwinc,                             
     &           num_crack_fronts*num_nodes_grwinc, num_nodes_thick )           
      call rd2d( fileno, master_lines, num_crack_fronts,                        
     &           num_crack_fronts, num_nodes_back + 1 )                         
c                                                                               
c              read the data smcs         
c                                                                               
      case( 11 )
      call rdbk( fileno, smcs_old_epsplas, num_kill_elem * prec_fact ) 
      call rdbk( fileno, smcs_tear_param, num_kill_elem * prec_fact )           
      read(fileno) isize
      if( isize > 0 ) then
          allocate(smcs_states_intlst(isize)) 
          read(fileno)  smcs_states_intlst(1:isize)   
      end if 
      if( use_weighted ) then
        isize = num_kill_elem * prec_fact 
        call rdbk( fileno, smcs_weighted_T, isize )
        call rdbk( fileno, smcs_weighted_zeta, isize )
        call rdbk( fileno, smcs_weighted_bar_theta, isize )
      end if  
c                                                                               
c              read the mesh regularization data
c                                                                               
      case( 12 )
      call rdbk( fileno, smcs_d_values, num_kill_elem * prec_fact )           
      call rdbk( fileno, smcs_eps_plas_at_death, 
     &           num_kill_elem * prec_fact )    
      call rdbk( fileno, smcs_stress_at_death, 
     &           num_kill_elem * prec_fact )
      call rdbk( fileno, smcs_start_kill_step, num_kill_elem )   
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
      contains
c     ========
c
      subroutine read_damage_dvector( vec, n )
      implicit none
c
      integer, intent(in) :: n
      double precision, intent(out) :: vec(n)
c
      read(fileno) vec(1:n)
c
      return
c
      end subroutine read_damage_dvector                                                                
      end subroutine read_damage                                                                       
                                                                                
                                                                                
