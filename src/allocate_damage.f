c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine allocate_damage              *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 4/21/23 rhd                *          
c     *                                                              *          
c     *     allocates information for the damage routines as needed  *   
c     *                                                              *          
c     ****************************************************************          
c  
c                                                                               
      subroutine allocate_damage ( dowhat )       
c                              
      use global_data, only: out, mxedof, nelblk, nonode, noelem,
     &                       mxconn, props, incflg
      use elem_extinct_data                                                     
      use node_release_data                                                     
      use damage_data   
      use constants
c                                                        
      implicit none 
c
      integer :: dowhat, ido                                                   
c                                                                               
      integer :: dum, elem, elem_ptr, i
      double precision ::                                                          
     &      dumd1, dumd2, dumd3, dumd4, dumd5, dumd6, dumd7,                      
     &     dumd8, dumd9, dumd10, plast_strain,                                              
     &     values(20)         
      logical :: debug, duml, standard_kill_method, get_princ                                                      
      real :: dumr
c
      include 'include_damage_values'   
      type(values_T) :: elem_values
c                                                                               
      debug = .false.
      standard_kill_method = .true.
      if( use_mesh_regularization ) standard_kill_method = .false. 
c                                                                               
c                               allocate dam_state, dam_ifv
c  
      select case( dowhat )       
c                                                                      
      case( 1 ) ! dowhat
c                                                                               
         if( allocated(dam_state) ) then                                           
            call errmsg(215,dum,'dam_state',dumr,dumd1)                            
            return
         end if        
         allocate( dam_state(num_kill_elem) )                                        
         dam_state(1:num_kill_elem) = 0                                            
c                                                                                  
         if( allocated(dam_ifv) ) then                                             
            call errmsg(215,dum,'dam_ifv',dumr,dumd1)                              
            return
         end if                                                                    
         allocate( dam_ifv(mxedof,num_kill_elem) )                                  
         dam_ifv = zero    
c                                                                                  
c                               allocate dam_print_list                         
c                               also allocate old_mises and                     
c                               old_mean to calculate change in                 
c                               stress over step.                               
c                                                                               
      case( 2 ) ! dowhat
         if( allocated(dam_print_list) ) deallocate( dam_print_list )              
         allocate( dam_print_list(num_print_list) )                                
c                                                                                 
         if( allocated(gt_old_mises) ) deallocate( gt_old_mises )                        
         allocate( gt_old_mises( num_print_list ) )                                   
         gt_old_mises(1:num_print_list) = zero                                        
c                                                                                 
         if( allocated( gt_old_mean ) ) deallocate( gt_old_mean )                        
         allocate( gt_old_mean( num_print_list ) )                                    
         gt_old_mean(1:num_print_list) = zero                                         
c                                                                               
c                                                                               
c                                allocate kill_order_list                       
c                                                                               
      case( 3 ) ! dowhat
         if( allocated(kill_order_list) ) deallocate(kill_order_list)              
         allocate(kill_order_list(num_kill_order_list))                            
c                                                                               
c                               allocate dam_node_elecnt                        
c                                                                               
      case( 4 ) ! dowhat
         if( allocated(dam_node_elecnt) ) then                                     
            call errmsg(215,dum,'dam_node_elecnt',dumr,dumd1)                      
            return
         end if                                                                    
         allocate(dam_node_elecnt(nonode))                                         
         dam_node_elecnt = 0                                                                               
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
      case( 5 ) ! dowhat                                                                
         if( allocated(dam_face_nodes) ) then                                      
            call errmsg(215,dum,'dam_face_nodes_elecnt',dumr,dumd1)                
            return
         end if                                                                    
         if( allocated(dam_dbar_elems) ) then                                      
            call errmsg(215,dum,'dam_dbar_elems',dumr,dumd1)                       
            return
         end if                                                                    
         allocate( dam_face_nodes(4,num_kill_elem) )                                 
         allocate( dam_dbar_elems(2,num_kill_elem) )                                 
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
      case( 6 ) ! dowhat                                                
         if( allocated(crack_plane_nodes) )
     &       deallocate( crack_plane_nodes )        
         allocate( crack_plane_nodes(num_crack_plane_nodes) )                      
         if( allocated(inv_crkpln_nodes) ) 
     &       deallocate( inv_crkpln_nodes )          
         allocate( inv_crkpln_nodes(nonode) )                                      
c                                                                               
         if( allocated(num_neighbors) ) deallocate( num_neighbors )                
         allocate( num_neighbors(num_crack_plane_nodes) )                          
c                                                                               
         if( allocated(neighbor_nodes) ) deallocate(neighbor_nodes)                
         allocate( neighbor_nodes(mxconn,num_crack_plane_nodes) )                  
c                                                                               
         if( allocated(crack_front_nodes) )
     &       deallocate( crack_front_nodes )        
         allocate( crack_front_nodes(num_crack_plane_nodes,2) )                    
c                                                                               
         if( allocated(crkpln_nodes_react) )                                       
     &         deallocate( crkpln_nodes_react )                                 
         allocate( crkpln_nodes_react(num_crack_plane_nodes) )                     
         if( allocated(crkpln_nodes_state) )                                       
     &         deallocate( crkpln_nodes_state )                                 
         allocate( crkpln_nodes_state(num_crack_plane_nodes) )                     
c                                                                               
c                            allocate: (for node_release crack growth           
c                                       and traction separation law)            
c                                node_release_frac(crack_plane_nodes) --        
c                                     stores the last fraction of release       
c                                     for a released node under traction        
c                                     separation.  We initialize it to          
c                                     zero in this routine                      
c                                                                               
      case( 7 ) ! dowhat
         if( allocated(node_release_frac) )
     &       deallocate( node_release_frac )        
         allocate( node_release_frac(num_crack_plane_nodes) )                      
         node_release_frac(1:num_crack_plane_nodes) = zero                         
c                                                                               
c                            allocate: (for overshoot control algorithm)        
c                                  old_angles_at_front -- Stores the old        
c                                       angle at each crack front node for      
c                                       each of its neighbors.                  
c                                  overshoot_allocated -- logical variable      
c                                       indicating if overshoot data            
c                                       structures have been allocated          
c                                                                               
      case( 8 ) ! dowhat
         if( allocated(old_angles_at_front) )                                      
     &       deallocate( old_angles_at_front )                                
         allocate( old_angles_at_front(num_crack_plane_nodes,                      
     &             mxconn) )                                                        
         old_angles_at_front(1:num_crack_plane_nodes,1:mxconn) = zero              
         overshoot_allocated = .true.                                              
c                                                                               
c   
c            element deletion. type 1 is Gurson. type 3 forms
c            of smcs. type 4 cohesive. type 5 user requested.'
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
      case( 9 ) ! dowhat
         select case( crack_growth_type )                                         
           case( 1 )   !   crack_growth_type  gurson                                                           
c                                                                               
c            gurson growth, initialize storage of old porosities                
c                                                                               
             if ( allocated( gt_old_porosity ) ) 
     &                deallocate( gt_old_porosity )               
             allocate( gt_old_porosity( num_kill_elem ) )                                 
c                                                                               
c             Now initialize the old_porosity with the current                  
c             porosity values...  If we are reallocating due to a               
c             restart, then don't bother to fill the array;                     
c             we will be refilling from the stored values. If                   
c             we are here before the computation of the first                   
c             load step, then fill array with the initial porosity.             
c                                                                               
             if ( g_stp_cntrl_allocated ) return
             if ( debug ) write(out,*) '>>> initializing porosity'                 
             if ( .not. incflg ) then                                                  
                do elem  = 1, noelem                                                    
                   elem_ptr = dam_ptr(elem)                                             
                   if( elem_ptr .eq. 0 ) cycle                                         
                   gt_old_porosity(elem_ptr) = props(26,elem)                          
                end do                                                                  
             else                                                                      
                do elem  = 1, noelem                                                    
                   elem_ptr = dam_ptr(elem)                                             
                   if ( elem_ptr .eq. 0 ) cycle                                         
                   call mm_return_values("avg_porosity", values )                             
                   gt_old_porosity(elem_ptr) = values(1)
                   if( debug )
     &                  write(out,'("   Poros. for:",i5,"=",e13.6)')             
     &                  elem, values(1)
                end do   ! elem                                                                
             end if ! incflg                                                                    
c                                                                               
           case( 3, 5 )  !   crack_growth_type  smcs, user                                                          
c                                                                               
c            crack growth is smcs, initialize storage of 
c            old plastic strains       
c                                                                               
              if( allocated( old_plast_strain ) )                                      
     &            deallocate( old_plast_strain )                                 
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
              if( g_stp_cntrl_allocated ) return
              if( .not. incflg ) then                                                  
                do elem_ptr = 1, num_kill_elem                                          
                   old_plast_strain(elem_ptr) = zero                                    
                end do                                                                  
              else                                                                      
                if( debug ) write(out,*) '>>> init plast_strain'                       
                do elem = 1, noelem                                                     
                   elem_ptr = dam_ptr( elem )                                           
                   if( elem_ptr .eq. 0 ) cycle                                         
                   ido = 1
                   get_princ = .false. 
                   call dam_param_smcs_get_values( elem, ido,  
     &                   elem_values, get_princ )    
                   old_plast_strain(elem_ptr) = elem_values%eps_plas
                   if( debug )
     &                  write(out,'("   Pl. e for:",i5,"is",e13.6)')             
     &                       elem, plast_strain                                          
                end do                                                                  
              end if  ! incflg                                                                    
c                                                                               
           case( 4 ) ! crack_growth_type  cohesive                                                             
c          
c             crack growth using cohesive elements                              
c                                                                               
              if( allocated(cohes_old_deff) ) deallocate(cohes_old_deff)                           
              allocate( cohes_old_deff(num_kill_elem) )                                       
c                                                                               
c             Now initialize the old_deff with the current                      
c             deff values...  If we are reallocating due to a                   
c             restart, then don't bother to fill the array;                     
c             we will be refilling from the stored values. If                   
c             we are here before the computation of the first                   
c             load step, then fill array with zero                              
c                                                                               
              if( g_stp_cntrl_allocated ) return
              if( debug ) write(out,*) '>>> initializing deff'                         
              if( .not. incflg ) then                                                  
                 do elem  = 1, noelem                                                    
                   elem_ptr = dam_ptr(elem)                                             
                   if( elem_ptr .eq. 0 ) cycle                                         
                   cohes_old_deff(elem_ptr) = zero                                        
                 end do                                                                  
              else                                                                      
                 do elem  = 1, noelem                                                    
                    elem_ptr = dam_ptr(elem)                                             
                    if( elem_ptr .eq. 0 ) cycle                                         
                    call dam_param_cohes( elem, duml, debug, values, 1 )                 
                    cohes_old_deff(elem_ptr) = values(6)                                       
                    if( debug )
     &                 write(out,'("   deff. for:",i5,"=",e13.6)')              
     &                       elem, values(6)                                             
                 end do                                                                  
              end if    ! incflg                                                                
c                                                                               
           case default  !  crack_growth_type                                                            
             write(out,*) '>>> invalid case in allocate_damage, type 9'                
             call die_abort         
c                                                 
         end select   !   crack_growth_type    
c                                                                               
      g_stp_cntrl_allocated = .true.                                            
c                                                                               
c                                                                               
c                            allocate: (constant front growth)                  
c                                  master_nodes -- holds the list of the        
c                                       master nodes, one per crack front,      
c                                       which will control the growth of        
c                                       its front                               
c                                                                               
      case( 10 ) !  dowhat                                                                              
         if( allocated(master_nodes) ) deallocate( master_nodes )                                       
         allocate( master_nodes(num_crack_fronts) )                                
         master_nodes(1:num_crack_fronts) = zero                                   
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
      case( 11 )  ! dowhat                                                              
c                                                                               
c           calculate number of nodes in between growth increments.             
c              
         num_nodes_grwinc = max( init_ctoa_dist,ctoa_dist ) /
     &                      char_length + 1   
         num_nodes_grwinc = max( num_nodes_grwinc,200 )                              
c                                                                               
         if( allocated(crack_front_list) )                                        
     &       deallocate( crack_front_list )                                       
         allocate( crack_front_list(num_crack_fronts * num_nodes_grwinc,           
     &             num_nodes_thick) )   
         crack_front_list(1:num_crack_fronts*num_nodes_grwinc,                     
     &                    1:num_nodes_thick) = 0                              
c                                                                               
         if( allocated(master_lines) ) deallocate( master_lines )                                           
         allocate( master_lines(num_crack_fronts, num_nodes_back + 1) )             
         master_lines(1:num_crack_fronts, 1:num_nodes_back + 1) = 0                
c                                                                               
c                            allocate: ( dam_ptr )                              
c                                                                               
       case( 12 )  ! dowhat
          if( allocated(dam_ptr) )  deallocate( dam_ptr )                                                
          allocate( dam_ptr(noelem) )                                               
          dam_ptr(1:noelem) = 0                                                     
c
c                            SMCS data                             
c                                                                               
       case( 13 )  ! dowhat
c
c                            data structures not needed for SMCS that
c                            may have been allocated.
c
          if( allocated( dam_dbar_elems ) )
     &        deallocate( dam_dbar_elems )               
          if( allocated( dam_face_nodes ) )
     &        deallocate( dam_face_nodes )     
          if( allocated( gt_old_porosity ) )
     &        deallocate( gt_old_porosity )    
c 
          if( allocated( smcs_old_epsplas ) ) 
     &       deallocate( smcs_old_epsplas )               
          if( allocated( smcs_tear_param ) ) 
     &        deallocate( smcs_tear_param  )  
c
          allocate( smcs_old_epsplas(num_kill_elem) )
          allocate( smcs_tear_param(num_kill_elem) )
c
          do i = 1, num_kill_elem
            smcs_old_epsplas(i) = zero  
            smcs_tear_param(i)  = zero  
          end do
c
c                            mesh regularization                             
c                                                                               
       case( 14 )  ! dowhat
c
c                            data structures not needed for SMCS that
c                            may have been allocated.
c
          if( allocated( dam_ifv ) )
     &        deallocate( dam_ifv )               
          if( allocated( dam_dbar_elems ) )
     &        deallocate( dam_dbar_elems )               
          if( allocated( dam_face_nodes ) )
     &        deallocate( dam_face_nodes )     
          if( allocated( gt_old_porosity ) )
     &        deallocate( gt_old_porosity )    
c 
          if( allocated( smcs_d_values ) ) 
     &        deallocate( smcs_d_values  )  
          if( allocated( smcs_eps_plas_at_death ) )
     &        deallocate( smcs_eps_plas_at_death  )  
          if( allocated( smcs_stress_at_death ) )
     &        deallocate( smcs_stress_at_death  )  
          if( allocated( killed_estiffs ) )
     &        deallocate( killed_estiffs )               
          if( allocated( smcs_start_kill_step ) )
     &        deallocate( smcs_start_kill_step )               
c
          allocate( smcs_d_values(num_kill_elem) )
          allocate( smcs_eps_plas_at_death(num_kill_elem) )
          allocate( smcs_stress_at_death(num_kill_elem) )
          allocate( killed_estiffs(num_kill_elem) )
          allocate( smcs_start_kill_step(num_kill_elem) )
c
          do i = 1, num_kill_elem
             smcs_d_values(i)            = zero
             smcs_eps_plas_at_death(i)   = zero
             smcs_stress_at_death(i)      = zero
             killed_estiffs(i)%num_terms = 0
             smcs_start_kill_step(i) = 0
          end do
c
c                            allocate Oddy distortion metrics
c                                                                               
       case( 15 )  ! dowhat
c
         if( allocated( Oddy_metrics ) ) deallocate( Oddy_metrics )
         if( use_distortion_metric ) then
             allocate( Oddy_metrics(num_kill_elem,2) )
             Oddy_metrics = 100000.0  ! single precision, large number
         end if  
c
c                            de-allocate Oddy distortion metrics
c                                                                               
       case( 16 )  ! dowhat
c
         if( allocated( Oddy_metrics ) ) deallocate( Oddy_metrics )
         if( allocated( Oddy_metrics_initial ) ) 
     &           deallocate( Oddy_metrics_initial )
c
c                            allocate Oddy distortion metrics at
c                            time = 0 for output to a file
c                                                                               
       case( 17 )  ! dowhat
c  
           if( Oddy_print_initial ) then
                if( allocated( Oddy_metrics_initial ) ) 
     &               deallocate(  Oddy_metrics_initial ) 
                allocate( Oddy_metrics_initial(num_kill_elem,2) )
                Oddy_metrics_initial(1:num_kill_elem,1) = 100000.0
                Oddy_metrics_initial(1:num_kill_elem,2) = -100000.0
           end if
c
      end select  !  dowhat   
c                                                                               
      return                                                                    
      end                                                                       
