c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine allocate_damage              *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 10/16/94                   *          
c     *                   last modified :  804/97 ag                 *          
c     *                                                              *          
c     *     this subroutine allocates all the information for        *          
c     *     the damage routines                                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine allocate_damage ( status )                                     
      use global_data ! old common.main
      use elem_extinct_data                                                     
      use node_release_data                                                     
      use damage_data                                                           
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
c                                                                               
c                                                                               
      double precision                                                          
     &     zero, dumd1, dumd2, dumd3, dumd4, dumd5, dumd6,                      
     &     porosity, plast_strain,                                              
     &     values(20)                                                           
      logical debug, duml                                                       
      data zero /0.0/                                                           
      real dumr                                                                 
c                                                                               
      debug = .false.                                                           
c                                                                               
c                                                                               
      go to (100,200,300,400,500,600,700,800,900,1000,1100,                     
     &       1200 ), status                                                     
c                                                                               
c                               allocate dam_state, dam_ifv,                    
c                               dam_blk_killed                                  
c                                                                               
 100  continue                                                                  
c                                                                               
      if( allocated(dam_state) ) then                                           
         call errmsg(215,dum,'dam_state',dumr,dumd1)                            
         go to 9999                                                             
      end if                                                                    
      allocate(dam_state(num_kill_elem))                                        
      dam_state(1:num_kill_elem) = 0                                            
c                                                                               
      if( allocated(dam_ifv) ) then                                             
         call errmsg(215,dum,'dam_ifv',dumr,dumd1)                              
         go to 9999                                                             
      end if                                                                    
      allocate( dam_ifv(mxedof,num_kill_elem))                                  
      dam_ifv(1:mxedof,1:num_kill_elem) = zero                                  
c                                                                               
      if( allocated(dam_blk_killed) ) then                                      
         call errmsg(215,dum,'dam_blk_killed',dumr,dumd1)                       
         go to 9999                                                             
      end if                                                                    
      allocate( dam_blk_killed (nelblk))                                        
      dam_blk_killed(1:nelblk) = .false.                                        
c                                                                               
      go to 9999                                                                
c                                                                               
c                               allocate dam_print_list                         
c                               also allocate old_mises and                     
c                               old_mean to calculate change in                 
c                               stress over step.                               
c                                                                               
 200  continue                                                                  
      if( allocated(dam_print_list) ) deallocate( dam_print_list )              
      allocate( dam_print_list(num_print_list) )                                
c                                                                               
      if( allocated(old_mises) ) deallocate( old_mises )                        
      allocate( old_mises( num_print_list ) )                                   
      old_mises(1:num_print_list) = zero                                        
c                                                                               
      if( allocated( old_mean ) ) deallocate( old_mean )                        
      allocate( old_mean( num_print_list ) )                                    
      old_mean(1:num_print_list) = zero                                         
c                                                                               
      go to 9999                                                                
c                                                                               
c                                                                               
c                                allocate kill_order_list                       
c                                                                               
 300  continue                                                                  
      if( allocated(kill_order_list) ) deallocate(kill_order_list)              
      allocate(kill_order_list(num_kill_order_list))                            
      go to 9999                                                                
c                                                                               
c                               allocate dam_node_elecnt                        
c                                                                               
 400  continue                                                                  
      if( allocated(dam_node_elecnt) ) then                                     
         call errmsg(215,dum,'dam_node_elecnt',dumr,dumd1)                      
         go to 9999                                                             
      end if                                                                    
      allocate(dam_node_elecnt(nonode))                                         
      dam_node_elecnt(1:nonode) = 0                                             
      go to 9999                                                                
c                                                                               
c                                                                               
c                               allocate dam_face_nodes and                     
c                               dam_dbar_elems. these store the                 
c                               nodes for each element to compute               
c                               average extension normal to crack plane,        
c                               the average extension when the                  
c                               element reaches critical porosity,              
c                               and the average extension on the                
c                               previous load step.                             
c                                                                               
 500  continue                                                                  
      if( allocated(dam_face_nodes) ) then                                      
         call errmsg(215,dum,'dam_face_nodes_elecnt',dumr,dumd1)                
         go to 9999                                                             
      end if                                                                    
      if( allocated(dam_dbar_elems) ) then                                      
         call errmsg(215,dum,'dam_dbar_elems',dumr,dumd1)                       
         go to 9999                                                             
      end if                                                                    
      allocate(dam_face_nodes(4,num_kill_elem))                                 
      allocate(dam_dbar_elems(2,num_kill_elem))                                 
      go to 9999                                                                
c                                                                               
c                                allocate: (for node_release crack growth)      
c                                  crack_plane_nodes -- nodes on crack plane    
c                                  inv_crkpln_nodes -- inverse of               
c                                      crack_plane_nodes                        
c                                  num_neighbors -- # of neighboring nodes      
c                                      for each crack plane node                
c                                  neighbor_nodes -- the neighboring nodes      
c                                  crack_front_nodes -- linked list of          
c                                      the nodes on the crack front             
c                                  crkpln_nodes_react -- an array to hold       
c                                      the reaction force at the released       
c                                      node at release time                     
c                                  crkpln_nodes_state -- the current release    
c                                      state for the killed nodes               
c                                                                               
 600  continue                                                                  
      if( allocated(crack_plane_nodes) ) deallocate( crack_plane_nodes )        
      allocate( crack_plane_nodes(num_crack_plane_nodes) )                      
      if( allocated(inv_crkpln_nodes) ) deallocate( inv_crkpln_nodes )          
      allocate( inv_crkpln_nodes(nonode) )                                      
c                                                                               
      if( allocated(num_neighbors) ) deallocate( num_neighbors )                
      allocate( num_neighbors(num_crack_plane_nodes) )                          
c                                                                               
      if( allocated(neighbor_nodes) ) deallocate(neighbor_nodes)                
      allocate( neighbor_nodes(mxconn,num_crack_plane_nodes) )                  
c                                                                               
      if( allocated(crack_front_nodes) ) deallocate( crack_front_nodes )        
      allocate( crack_front_nodes(num_crack_plane_nodes,2) )                    
c                                                                               
      if( allocated(crkpln_nodes_react) )                                       
     &         deallocate( crkpln_nodes_react )                                 
      allocate( crkpln_nodes_react(num_crack_plane_nodes) )                     
      if( allocated(crkpln_nodes_state) )                                       
     &         deallocate( crkpln_nodes_state )                                 
      allocate( crkpln_nodes_state(num_crack_plane_nodes) )                     
      go to 9999                                                                
c                                                                               
c                            allocate: (for node_release crack growth           
c                                       and traction separation law)            
c                                node_release_frac(crack_plane_nodes) --        
c                                     stores the last fraction of release       
c                                     for a released node under traction        
c                                     separation.  We initialize it to          
c                                     zero in this routine                      
c                                                                               
 700  continue                                                                  
      if( allocated(node_release_frac) ) deallocate( node_release_frac )        
      allocate( node_release_frac(num_crack_plane_nodes) )                      
      node_release_frac(1:num_crack_plane_nodes) = zero                         
      go to 9999                                                                
c                                                                               
c                            allocate: (for overshoot control algorithm)        
c                                  old_angles_at_front -- Stores the old        
c                                       angle at each crack front node for      
c                                       each of its neighbors.                  
c                                  overshoot_allocated -- logical variable      
c                                       indicating if overshoot data            
c                                       structures have been allocated          
c                                                                               
 800  continue                                                                  
      if( allocated(old_angles_at_front) )                                      
     &         deallocate( old_angles_at_front )                                
      allocate( old_angles_at_front(num_crack_plane_nodes,                      
     &         mxconn) )                                                        
      old_angles_at_front(1:num_crack_plane_nodes,1:mxconn) = zero              
      overshoot_allocated = .true.                                              
      go to 9999                                                                
c                                                                               
c                                                                               
c                            allocate: (for step size control)                  
c                                g_stp_cntrl_allocated: indicates that          
c                                     the structures for step size control      
c                                     for crack growth or have been             
c                                     allocated.                                
c                                - for gurson crack growth:                     
c                                     old_porosity(num_kill_elem) --            
c                                       stores the porosity at the previous     
c                                       step for comparison to current          
c                                       porosity                                
c                                - for smcs crack growth:                       
c                                     old_plast_strain(num_kill_elem) --        
c                                       stores the plastic strain at the        
c                                       previous step for comparison to         
c                                       current plast strain                    
c                                - for cohesive crack growth:                   
c                                     old_deff(num_kill_elem) --                
c                                       stores the effective (relative)         
c                                       displacement at the previous step       
c                                       for comparison to the current value     
c                                                                               
 900  continue                                                                  
      select case ( crack_growth_type )                                         
      case( 1 )                                                                 
c                                                                               
c            gurson growth, initialize storage of old porosities                
c                                                                               
      if ( allocated( old_porosity ) ) deallocate( old_porosity )               
      allocate( old_porosity( num_kill_elem ) )                                 
c                                                                               
c             Now initialize the old_porosity with the current                  
c             porosity values...  If we are reallocating due to a               
c             restart, then don't bother to fill the array;                     
c             we will be refilling from the stored values. If                   
c             we are here before the computation of the first                   
c             load step, then fill array with the initial porosity.             
c                                                                               
      if ( g_stp_cntrl_allocated ) go to 9999                                   
      if ( debug ) write ( out, * ) '>>> initializing porosity'                 
      if ( .not. incflg ) then                                                  
        do elem  = 1, noelem                                                    
           elem_ptr = dam_ptr(elem)                                             
           if ( elem_ptr .eq. 0 ) cycle                                         
               old_porosity(elem_ptr) = props(26,elem)                          
        end do                                                                  
      else                                                                      
        do elem  = 1, noelem                                                    
           elem_ptr = dam_ptr(elem)                                             
           if ( elem_ptr .eq. 0 ) cycle                                         
           call dam_param( elem, duml, debug, porosity, dumd1,                  
     &                     dumd2, dumd3, dumd4, duml, dumd5,                    
     &                     dumd6 )                                              
           old_porosity( elem_ptr ) = porosity                                  
           if ( debug ) write (*,'("   Poros. for:",i5,"=",e13.6)')             
     &              elem, porosity                                              
        end do                                                                  
      end if                                                                    
c                                                                               
      case( 3 )                                                                 
c                                                                               
c         crack growth is smcs, initialize storage of old plastic strains       
c                                                                               
      if ( allocated( old_plast_strain ) )                                      
     &           deallocate( old_plast_strain )                                 
      allocate( old_plast_strain( num_kill_elem ) )                             
c                                                                               
c             Now initialize old_plast_strain with the current                  
c             plastic strain values...  if we are reallocating                  
c             due to a restart, then don't bother to fill the array;            
c             we will be refilling from the stored values. If                   
c             we are here before the computation of the first                   
c             load step, then fill array with zeroes.                           
c                                                                               
c                                                                               
      if ( g_stp_cntrl_allocated ) go to 9999                                   
      if ( .not. incflg ) then                                                  
        do elem_ptr = 1, num_kill_elem                                          
           old_plast_strain(elem_ptr) = zero                                    
        end do                                                                  
      else                                                                      
        if ( debug ) write(out,*) '>>> init plast_strain'                       
        do elem = 1, noelem                                                     
           elem_ptr = dam_ptr( elem )                                           
           if ( elem_ptr .eq. 0 ) cycle                                         
           call dam_param( elem, duml, debug, dumd5,                            
     &                     plast_strain, dumd2, dumd3, dumd4, duml,             
     &                     dumd5, dumd6 )                                       
           old_plast_strain( elem_ptr ) = plast_strain                          
           if ( debug ) write (*,'("   Pl. e for:",i5,"is",e13.6)')             
     &              elem, plast_strain                                          
        end do                                                                  
      end if                                                                    
c                                                                               
      case( 4 )                                                                 
c                                                                               
c             crack growth using cohesive elements                              
c                                                                               
      if ( allocated(old_deff) ) deallocate(old_deff)                           
      allocate( old_deff(num_kill_elem) )                                       
c                                                                               
c             Now initialize the old_deff with the current                      
c             deff values...  If we are reallocating due to a                   
c             restart, then don't bother to fill the array;                     
c             we will be refilling from the stored values. If                   
c             we are here before the computation of the first                   
c             load step, then fill array with zero                              
c                                                                               
      if ( g_stp_cntrl_allocated ) go to 9999                                   
      if ( debug ) write(out,*) '>>> initializing deff'                         
      if ( .not. incflg ) then                                                  
        do elem  = 1, noelem                                                    
           elem_ptr = dam_ptr(elem)                                             
           if ( elem_ptr .eq. 0 ) cycle                                         
               old_deff(elem_ptr) = zero                                        
        end do                                                                  
      else                                                                      
        do elem  = 1, noelem                                                    
           elem_ptr = dam_ptr(elem)                                             
           if ( elem_ptr .eq. 0 ) cycle                                         
           call dam_param_cohes( elem, duml, debug, values, 1 )                 
           old_deff(elem_ptr) = values(6)                                       
           if ( debug ) write (*,'("   deff. for:",i5,"=",e13.6)')              
     &              elem, values(6)                                             
        end do                                                                  
      end if                                                                    
c                                                                               
      case default                                                              
        write(*,*) '>>> invalid case in allocate_damage, type 9'                
        call die_abort                                                          
        stop                                                                    
      end select                                                                
c                                                                               
      g_stp_cntrl_allocated = .true.                                            
      go to 9999                                                                
c                                                                               
c                                                                               
c                            allocate: (constant front growth)                  
c                                  master_nodes -- holds the list of the        
c                                       master nodes, one per crack front,      
c                                       which will control the growth of        
c                                       its front                               
c                                                                               
 1000  continue                                                                 
c                                                                               
      if ( allocated(master_nodes) )                                            
     &         deallocate( master_nodes )                                       
      allocate( master_nodes(num_crack_fronts) )                                
      master_nodes(1:num_crack_fronts) = zero                                   
c                                                                               
      goto 9999                                                                 
c                                                                               
c                                                                               
c                            allocate: (constant front growth)                  
c                                  crack_front_list -- holds the list of the    
c                                       nodes on each of the crack fronts.      
c                                  master_lines -- holds a list of              
c                                       num_nodes_back nodes of the nodes       
c                                       on the master line behind each          
c                                       master node.                            
c                                                                               
 1100 continue                                                                  
c                                                                               
c           calculate number of nodes in between growth increments.             
c                                                                               
c   ==========> better comments needed <=================                       
c                                                                               
      num_nodes_grwinc = max(init_ctoa_dist,ctoa_dist)/char_length + 1          
      num_nodes_grwinc = max(num_nodes_grwinc,100)                              
                                                                                
c                                                                               
      if ( allocated(crack_front_list) )                                        
     &     deallocate( crack_front_list )                                       
      allocate( crack_front_list(num_crack_fronts * num_nodes_grwinc,           
     &     num_nodes_thick) )                                                   
      crack_front_list(1:num_crack_fronts*num_nodes_grwinc,                     
     &                      1:num_nodes_thick) = 0                              
c                                                                               
      if ( allocated(master_lines) )                                            
     &     deallocate( master_lines )                                           
      allocate( master_lines(num_crack_fronts, num_nodes_back + 1))             
      master_lines(1:num_crack_fronts, 1:num_nodes_back + 1) = 0                
c                                                                               
      goto 9999                                                                 
c                                                                               
c                            allocate: ( dam_ptr )                              
c                                                                               
 1200 continue                                                                  
      if( allocated(dam_ptr) )                                                  
     &     deallocate( dam_ptr )                                                
      allocate( dam_ptr(noelem) )                                               
      dam_ptr(1:noelem) = 0                                                     
      go to 9999                                                                
c                                                                               
 9999 continue                                                                  
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
