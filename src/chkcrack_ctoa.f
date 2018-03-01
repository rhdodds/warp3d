c ********************************************************************          
c *                                                                  *          
c *  routines to support crack growth by critical ctoa criterion     *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_print_node               *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 04/22/02                   *          
c     *                                                              *          
c     *     This routine prints out the status of the crack front    *          
c     *     nodes for a crack growth by node release model at the    *          
c     *     beginning of a load step.                                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dam_print_node( step, iter )                                   
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : crack_front_nodes, inv_crkpln_nodes,        
     &     num_neighbors, neighbor_nodes, crkpln_nodes_state                    
      use main_data, only : cnstrn, cnstrn_in, output_packets,                  
     &                      packet_file_no                                      
c                                                                               
      use damage_data                                                           
      implicit integer (a-z)                                                    
      parameter (max_local_list=300)                                            
      double precision                                                          
     &     d32460, angle, zero, max_angle, two,                                 
     &     local_angles_list(max_local_list), denom_angle                       
      dimension local_node_list(max_local_list,3)                               
      character(len=1) :: star                                                  
      data d32460, zero, two, star / 32460.0, 0.0, 2.0, '*'/                    
      logical first_entry                                                       
c                                                                               
c           print the header                                                    
c                                                                               
      write(out,*) ' '                                                          
      write(out,*) ' ********************************************* '            
      write(out,*) ' ***         Crack Front Status            *** '            
      write(out,*) ' ********************************************* '            
      write(out,*)                                                              
c                                                                               
      write (out,*) 'crack front node  neighbor  CTOA(degrees)',                
     &     '  CTOA/critical CTOA'                                               
      write (out,*) '----------------  --------  -------------',                
     &     '  ------------------'                                               
c                                                                               
c          set up for parsing the crack front node list                         
c                                                                               
      pointer = crack_front_start                                               
      if( pointer .eq. -1 ) then                                                
         write (out,*) '>>>>> No crack front nodes'                             
         go to 9999                                                             
      end if                                                                    
      max_angle      = zero                                                     
      max_node       = 0                                                        
      max_neighbor   = 0                                                        
      num_local_list = 0                                                        
c                                                                               
c          check next crack front node                                          
c                                                                               
 10   continue                                                                  
      first_entry     = .true.                                                  
      node            = crack_front_nodes(pointer,1)                            
      node_data_entry = inv_crkpln_nodes(node)                                  
c                                                                               
      do neighbor = 1, num_neighbors(node_data_entry)                           
c                                                                               
c                find out if current neighbor is constrained -- if so,          
c                skip to next node.  if not, must compute the angle             
c                                                                               
         neighbor_node = neighbor_nodes(neighbor,node_data_entry)               
         if( cnstrn(dstmap(neighbor_node)+crk_pln_normal_idx-1)                 
     &        .ne. d32460 ) cycle                                               
c                                                                               
c               neighbor node is unconstrained.  call routine to get            
c               the angle between it and the original node. Print info.         
c                                                                               
c               note that the angle returned is between the line and            
c               the crack plane -- to make it relate to the CTOA, we            
c               must multiply it by two.                                        
c                                                                               
         call get_slope( neighbor_node, node, crk_pln_normal_idx,               
     &                   angle )                                                
         angle = angle * two                                                    
         if( angle .gt. max_angle ) then                                        
            max_angle    = angle                                                
            max_node     = node                                                 
            max_neighbor = neighbor_node                                        
         end if                                                                 
c                                                                               
         num_local_list = num_local_list + 1                                    
         if ( num_local_list .gt. max_local_list ) then                         
           write(out,9030)                                                      
           call die_gracefully                                                  
         end if                                                                 
c                                                                               
         local_node_list(num_local_list,1) = node                               
         local_node_list(num_local_list,2) = neighbor_node                      
         local_node_list(num_local_list,3) = 0                                  
         local_angles_list(num_local_list) = angle                              
                                                                                
                                                                                
         if( crkpln_nodes_state(inv_crkpln_nodes(neighbor_node) )               
     &        .eq. 0 ) then                                                     
            if( first_entry ) then                                              
               write(out,"(6x,i5,9x,i5,7x,f7.2,a1,11x,f4.1,3x,'<=')")           
     &              node, neighbor_node,                                        
     &              angle, star, angle/(init_crit_ang*two)                      
               first_entry = .false.                                            
            else                                                                
               write(out,"(20x,i5,7x,f7.2,a1,11x,f4.1,3x,'<=')")                
     &              neighbor_node,                                              
     &              angle, star, angle/(init_crit_ang*two)                      
            end if                                                              
         else                                                                   
            if( first_entry ) then                                              
               write(out,"(6x,i5,9x,i5,7x,f7.2,12x,f4.1,3x,'<=')")              
     &              node, neighbor_node,                                        
     &              angle,angle/(critical_angle*two)                            
               first_entry = .false.                                            
               local_node_list(num_local_list,3) = 1                            
            else                                                                
               write(out,"(20x,i5,7x,f7.2,12x,f4.1,3x,'<=')")                   
     &              neighbor_node,                                              
     &              angle, angle/(critical_angle*two)                           
               local_node_list(num_local_list,3) = 1                            
            end if                                                              
         end if                                                                 
c                                                                               
      end do                                                                    
c                                                                               
c         now go to next entry in list                                          
c                                                                               
      pointer = crack_front_nodes(pointer,2)                                    
      if( pointer .ne. -1 ) go to 10                                            
c                                                                               
      write(out,9000) max_angle, max_node, max_neighbor                         
      write(out,9050) num_ctoa_released_nodes                                   
                                                                                
c                                                                               
      if( output_packets .and. num_local_list .gt. 0  ) then                    
        write(packet_file_no) 9, num_local_list, step, iter                     
        do i = 1, num_local_list                                                
          node        = local_node_list(i,1)                                    
          neighbor    = local_node_list(i,2)                                    
          use_init    = local_node_list(i,3)                                    
          angle       = local_angles_list(i)                                    
          denom_angle = critical_angle*two                                      
          if ( use_init .eq. 0 ) denom_angle = init_crit_ang*two                
          write(packet_file_no) node, neighbor, angle, use_init,                
     &                          angle/(denom_angle)                             
        end do                                                                  
      end if                                                                    
c                                                                               
 9999 continue                                                                  
      write (out,*)                                                             
      return                                                                    
c                                                                               
 9000 format (/1x,'** Max CTOA: ',f6.2,' degrees',                              
     &        ' @ crack front node: ',i7,                                       
     &        ' connected to node: ',i7,                                        
     &      /,1x,'   NOTE: * indicates critical CTOA is for',                   
     &        ' initiation',                                                    
     &        / )                                                               
 9030 format(/,'FATAL ERROR: list length exceeded in',                          
     & ' dam_print_node. job aborted.' )                                        
 9050 format(8x,'** Total nodes released in model: ',i5,' **')                  
c                                                                               
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_print_front              *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 12/18/2016 rhd             *          
c     *                                                              *          
c     *     This routine prints out the status of the crack front    *          
c     *     nodes for a crack growth by node release model at the    *          
c     *     beginning of a load step if it uses a constant front.    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dam_print_front( step, iter )                                  
      use global_data ! old common.main
      use node_release_data, only : master_lines                                
      use main_data, only : output_packets, packet_file_no                      
      use damage_data                                                           
c                                                                               
      implicit none                                                             
                                                                                
      integer :: step, iter                                                     
                                                                                
      integer, parameter ::  max_local_list = 200                               
      double precision ::                                                       
     &     d32460, angle, zero, max_angle, two,                                 
     &     local_angles_list(max_local_list)                                    
      integer :: i, local_master_list(max_local_list), max_node,                
     &           max_neighbor, num_local_list, num_line, node, idummy           
      character(len=3) :: special_char                                          
      data d32460, zero, two / 32460.0d0, 0.0d0, 2.0d0 /                        
      logical :: use_init                                                       
c                                                                               
c           print the header                                                    
c                                                                               
      write(out,*) ' '                                                          
      write(out,*) ' ********************************************* '            
      write(out,*) ' ***         Crack Front Status            *** '            
      write(out,*) ' ********************************************* '            
      write(out,*)                                                              
c                                                                               
      write(out,*) 'master node       CTOA(degrees)',                           
     &     '  CTOA/critical CTOA'                                               
      write(out,*) '-----------       -------------',                           
     &     '  ------------------'                                               
c                                                                               
c          loop over the entries in master_lines                                
c                                                                               
      max_angle      = zero                                                     
      max_node       = 0                                                        
      max_neighbor   = 0                                                        
      num_local_list = 0                                                        
c                                                                               
      do num_line = 1, num_crack_fronts                                         
c                                                                               
         call get_slope_master_line( num_line, use_init, angle,                 
     &                               idummy )                                   
         if ( angle .eq. zero ) cycle                                           
         node  = master_lines( num_line, 1)                                     
         angle = angle * two                                                    
         if( angle .gt. max_angle ) then                                        
            max_angle = angle                                                   
            max_node = node                                                     
         end if                                                                 
         special_char(1:3) = ' '                                                
         if( use_init ) special_char(1:3) = '(*)'                               
         if( use_init ) then                                                    
          write(out,"(6x,i5,7x,f7.2,a3,11x,f4.1,3x,'<=')")                      
     &        node, angle, special_char, angle/(init_crit_ang*two)              
         else                                                                   
          write(out,"(6x,i5,7x,f7.2,a3,11x,f4.1,3x,'<=')")                      
     &        node, angle, special_char, angle/(critical_angle*two)             
         end if                                                                 
c                                                                               
c          save into vector for packet output after loop                        
c                                                                               
         num_local_list = num_local_list + 1                                    
         if ( num_local_list .gt. max_local_list ) then                         
           write(out,9030)                                                      
           call die_gracefully                                                  
         end if                                                                 
         local_master_list(num_local_list) = node                               
         local_angles_list(num_local_list) = angle                              
c                                                                               
      end do                                                                    
c                                                                               
      write(out,9000) max_angle, max_node                                       
      write(out,9050) num_ctoa_released_nodes                                   
      write(out,*)                                                              
c                                                                               
      if( output_packets .and. num_local_list .gt. 0  ) then                    
        write(packet_file_no) 8, num_local_list, step, iter                     
        do i = 1, num_local_list                                                
          node  = local_master_list(i)                                          
          angle = local_angles_list(i)                                          
          if( use_init ) then                                                   
           write(packet_file_no) node, angle, 0,                                
     &                           angle/(init_crit_ang*two)                      
          else                                                                  
           write(packet_file_no) node, angle, 1,                                
     &                           angle/(critical_angle*two)                     
          end if                                                                
        end do                                                                  
      end if                                                                    
                                                                                
      return                                                                    
 9000 format (/1x,'** Max CTOA: ',f6.2,' degrees',                              
     &        ' @ crack front node: ',i7,                                       
     &      /,1x,'   NOTE: * indicates critical CTOA is for',                   
     &        ' initiation',                                                    
     &        / )                                                               
 9030 format(/,'FATAL ERROR: list length exceeded in',                          
     & ' dam_print_front. job aborted.' )                                       
 9050 format(8x,'** Total nodes released in model: ',i5,' **')                  
c                                                                               
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine chk_node_release             *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 06/03/94                   *          
c     *                   last modified : 07/28/95                   *          
c     *                                                              *          
c     *        This routine updates previously released nodes and    *          
c     *        checks the current crack front for crack growth       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine chk_node_release( debug, step, iter )                          
      use global_data ! old common.main
      use node_release_data, only: crkpln_nodes_state                           
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
      logical debug, killed_this_time                                           
c                                                                               
c                                                                               
c     This routine:                                                             
c        1) updates information for previously released nodes by                
c           - incrementing the state variable for each released node            
c           - eliminiating any new constraints placed on previously             
c              released degrees of freedom                                      
c                                                                               
c        2) finds new nodes for releaseing and executes release through         
c            two passes of the linked list for the crack front:                 
c               Pass 1: compute CTOA, mark node for release if                  
c                         CTOA > critical angle and alter data structures       
c                  for the released node                                        
c               Pass 2: update crack front list to reflect newly released       
c                       nodes                                                   
c                                                                               
c                                                                               
      if( debug ) write(out,*) '>>>> entering chk_node_release'                 
c                                                                               
c                                                                               
c        1) Update previously released nodes                                    
c                                                                               
c             loop through all crack plane nodes -- if node                     
c             has been released, increment state by one. This                   
c             keeps a record of the load steps since release                    
c             of the node, used by the load size reduction                      
c             mechanism.                                                        
c                                                                               
c             Even though the exact state value is not used for                 
c             traction separation, we still increment the state                 
c             value since it shows if a node has been released                  
c             or not.                                                           
c                                                                               
      if( debug ) write(out,*) '>>>>>> Checking released nodes '                
c                                                                               
      do i = 1, num_crack_plane_nodes                                           
         if ( crkpln_nodes_state(i).gt.0 )                                      
     &        crkpln_nodes_state(i) = crkpln_nodes_state(i) + 1                 
      end do                                                                    
c                                                                               
c                                                                               
c        2) Check nodes on crack front for release.  If any are released,       
c           then update the crack front list to reflect the changes.            
c           Note that if we are using constant front growth                     
c           we call a different check front routine.                            
c                                                                               
      if( debug ) write(out,*) '>>>>>> Checking crack front nodes'              
      if( const_front ) then                                                    
         call chk_crack_front_const( killed_this_time, debug,                   
     &                               step, iter )                               
      else                                                                      
         call chk_crack_front( killed_this_time, debug, step, iter )            
      end if                                                                    
c                                                                               
      if ( killed_this_time ) growth_k_flag = .true.                            
c                                                                               
      enforce_node_release = .false.                                            
c                                                                               
      if( debug ) write(out,*) '<<<< leaving chk_node_release'                  
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine chk_reconstraint             *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 10/10/95                   *          
c     *                                                              *          
c     *        This routine checks the current constraints to see    *          
c     *        if a previously released node has been re-constrained *          
c     *        in the direction normal to the crack plane            *          
c     *        by user-input constraints.  If a released node has    *          
c     *        been re-constrained, unconstrain it.                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine chk_reconstraint                                               
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : inv_crkpln_nodes,                           
     &     crkpln_nodes_state, crkpln_nodes_react                               
      use main_data, only : invdst, cnstrn_in, cnstrn                           
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      character(len=1) :: dums                                                  
      real dumr                                                                 
      double precision                                                          
     &     d32460, dumd                                                         
      data d32460 / 32460.0 /                                                   
c                                                                               
c             traverse the singly-linked list of constraints in the model. If   
c             the user has input new constraints, then some of the previously   
c             released nodes may be re-constrained.  If a released node         
c             has been re-constrained, remove the constraint.                   
c                                                                               
c                set up pointers for constraint linked list                     
c                                                                               
      cst_ptr   = csthed                                                        
      above_ptr = -1                                                            
c                                                                               
c                enter top of constraint linked list loop.                      
c                                                                               
 10   continue                                                                  
      node       = invdst(cst_ptr)                                              
      crkpln_ptr = inv_crkpln_nodes(node)                                       
c                                                                               
c                    if node is not on crack plane, or the dof is not in        
c                    the crack plane direction, or the node is not yet          
c                    killed, then go to next entry                              
c                                                                               
c      if ( (crkpln_ptr .eq. 0) .or.                                            
c     &     (cst_ptr .ne. dstmap(node) + crk_pln_normal_idx -1) .or.            
c     &     (crkpln_nodes_state(crkpln_ptr) .eq. 0) ) go to 20                  
      if (crkpln_ptr .eq. 0) go to 20                                           
      if (cst_ptr .ne. dstmap(node) + crk_pln_normal_idx -1) go to 20           
      if (crkpln_nodes_state(crkpln_ptr) .eq. 0) go to 20                       
c                                                                               
c                    the node has been re-constrained and must be               
c                    re-released.  write a warning message,                     
c                    set the constraint value to null (32460),                  
c                    and remove it from the linked list.  Removing a            
c                    constraint from linked list automatically                  
c                    updates the traversal pointer to point to the next         
c                    entry                                                      
c                                                                               
      call errmsg( 248, node, dums, dumr, dumd )                                
      cnstrn(cst_ptr)    = d32460                                               
      cnstrn_in(cst_ptr) = d32460                                               
c                                                                               
      if ( above_ptr .eq.-1 ) then                                              
c                                                                               
c                           at top of list.  move head pointer.                 
c                                                                               
         csthed          = cstmap(cst_ptr)                                      
         cstmap(cst_ptr) = 0                                                    
         cst_ptr         = csthed                                               
      else                                                                      
c                                                                               
c                           in middle of list. correct pointers                 
c                                                                               
         cstmap(above_ptr) = cstmap(cst_ptr)                                    
         cstmap(cst_ptr)   = 0                                                  
         if ( csttail .eq. cst_ptr ) csttail = above_ptr                        
         cst_ptr = cstmap(above_ptr)                                            
      end if                                                                    
      go to 30                                                                  
c                                                                               
c                    the constraint was not corrected.  we need to              
c                    move to the next entry in the list.                        
c                                                                               
 20   continue                                                                  
      above_ptr = cst_ptr                                                       
      cst_ptr   = cstmap(cst_ptr)                                               
c                                                                               
c                    check to see if we are at end of list. if not, go to       
c                    next value.                                                
c                                                                               
 30   continue                                                                  
      if ( cst_ptr .ne. -1 ) go to 10                                           
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine chk_crack_front              *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 5/29/02, rhd               *          
c     *                                                              *          
c     *        This routine traverses the nodes on the current       *          
c     *        crack front, releasing any nodes that have (1) CTOA   *          
c     *        larger than the given critical angle, or (2) the      *          
c     *        user has specifically requested a node release in     *          
c     *        this step.                                            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine chk_crack_front( killed_this_time, debug, step, iter )         
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : inv_crkpln_nodes, num_neighbors,            
     &     neighbor_nodes, crack_front_nodes, crkpln_nodes_state                
      use main_data, only : cnstrn, output_packets, packet_file_no              
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
      parameter (max_local_list=200)                                            
c                                                                               
      double precision                                                          
     &     d32460, angle, two, one, crit_angle, hundred,                        
     &     local_killed_angles(max_local_list)                                  
      dimension local_killed_list(max_local_list)                               
      data d32460, two, one, hundred / 32460.0, 2.0, 1.0, 100.0/                
      logical debug, kill, killed_this_time                                     
c                                                                               
c          traverse the linked list of crack front nodes and                    
c          check if nodes need to be released.                                  
c                                                                               
c             set pointer to top of list, and the previous node                 
c             pointer to null.  Also, set the flag that shows if                
c             any nodes have been released this step to false.                  
c                                                                               
      killed_this_time = .false.                                                
      num_killed_now   = 0                                                      
      prev_node_ptr    = -1                                                     
      node_ptr         = crack_front_start                                      
      if( node_ptr .eq. -1 ) return                                             
c                                                                               
c             for current node, find out if it is to be                         
c             released. loop over the neighbors, calculating                    
c             the opening angle (from the crack plane to the connected          
c             element edges).                                                   
c                                                                               
 10   continue                                                                  
c                                                                               
      node            = crack_front_nodes(node_ptr,1)                           
      node_data_entry = inv_crkpln_nodes(node)                                  
      kill            = .false.                                                 
      if( debug ) then                                                          
        write(out,*) ' >> Checking node ', node                                 
        write(out,*) '    > enforced release flag: ',                           
     &                   enforce_node_release                                   
      end if                                                                    
c                                                                               
      do neighbor = 1, num_neighbors(node_data_entry)                           
c                                                                               
c                find out if current neighbor is constrained -- if so,          
c                skip to next node.  if not, compute the angle.  If             
c                angle is larger than the critical angle (within a              
c                range of CTOA_range of the critical angle), signal             
c                that the node should be killed.  Note that we have two         
c                critical angles -- one for initiation and one for              
c                crack growth.  We use the initiation angle on any              
c                neighboring nodes that are unconstrained but have a state      
c                of 0, indicating that they were initially unconstrained.       
c                                                                               
         neighbor_node = neighbor_nodes(neighbor,node_data_entry)               
         dof           = dstmap(neighbor_node)+crk_pln_normal_idx-1             
         if( cnstrn(dof) .ne. d32460 ) cycle                                    
c                                                                               
         call get_slope( neighbor_node, node, crk_pln_normal_idx,               
     &        angle )                                                           
         if( debug ) write (out,'("     Angle is: ",e13.6)') angle              
c                                                                               
         if( crkpln_nodes_state(inv_crkpln_nodes(neighbor_node))                
     &        .eq. 0) then                                                      
            if(debug) write (out,*) ' using initiation angle'                   
            crit_angle = init_crit_ang * ( one - CTOA_range )                   
         else                                                                   
            if(debug) write (out,*) ' using crack growth angle'                 
            crit_angle = critical_angle * ( one - CTOA_range )                  
         end if                                                                 
c                                                                               
         if ( angle .ge. crit_angle .or. enforce_node_release ) then            
            kill = .true.                                                       
            exit                                                                
         end if                                                                 
c                                                                               
      end do                                                                    
c                                                                               
c           if node should be released, then print message and release it       
c                                                                               
c           note that the angle as computed in WARP3D is between the            
c           crack plane and the appropriate element edge.  When we output       
c           we need to multiply by two so we get the total opening angle        
c           for the specimen.                                                   
c                                                                               
      if( kill ) then                                                           
         if( debug ) write(out,*) ' -> Crit.CTOA reached. Kill node.'           
         if( .not. killed_this_time ) write (out,9000)                          
         write(out,9010) node, angle * two                                      
         killed_this_time = .true.                                              
         call release_node( node, node_data_entry, debug )                      
         num_killed_now = num_killed_now + 1                                    
         if ( num_killed_now .gt. max_local_list ) then                         
           write(out,9030)                                                      
           call die_gracefully                                                  
         end if                                                                 
         local_killed_list(num_killed_now) = node                               
         local_killed_angles(num_killed_now) = angle * two                      
      end if                                                                    
c                                                                               
c            move to the next entry in the list.  If we are not at the end,     
c            process the next node.                                             
c                                                                               
      prev_node_ptr = node_ptr                                                  
      node_ptr      = crack_front_nodes(node_ptr,2)                             
      if( node_ptr .ne. -1 ) goto 10                                            
c                                                                               
c            print critical CTOA if any elements were killed.  Note that        
c            again we multiply the angle by two before we output so that        
c            we are outputting the total critical CTOA.  Also update the        
c            crack front to refect killed elements.                             
c                                                                               
c                                                                               
      if( killed_this_time ) then                                               
         if( enforce_node_release ) then                                        
           write(out,9025)                                                      
         else                                                                   
           write(out,9020) critical_angle * two,                                
     &          init_crit_ang * two, CTOA_range * hundred                       
         end if                                                                 
         call update_crack_front( debug )                                       
         write(out,9050)  num_ctoa_released_nodes                               
      end if                                                                    
c                                                                               
      if( output_packets .and. killed_this_time ) then                          
        write(packet_file_no) 18, num_killed_now+1, step, iter                  
        write(packet_file_no) critical_angle*two,                               
     &        init_crit_ang*two, CTOA_range*hundred                             
        do i = 1, num_killed_now                                                
          write(packet_file_no) local_killed_list(i),                           
     &         local_killed_angles(i)                                           
        end do                                                                  
      end if                                                                    
                                                                                
      if( debug ) write(out,*) '>> leaving chk_crack_front'                     
c                                                                               
c                                                                               
      return                                                                    
 9000 format(/,' >> node release option invoked for the following:',/)          
 9010 format(  '        node: ',i7,'  CTOA: ',f6.2,' (degrees)')                
 9020 format(/,'        * critical CTOA:',                                      
     &       /,'            for growth    : ',f6.2,' (degrees)',                
     &       /,'            for initiation: ',f6.2,' (degrees)',                
     &       /,'            nodes released within ',f4.1,' % of ',              
     &                        'critical angle.',                                
     &       /)                                                                 
 9025 format(/,'        *** user requested node release ***', / )               
 9030 format(/,'FATAL ERROR: list length exceeded in chk_crack_front',          
     & /,      '             job aborted.' )                                    
 9050 format(8x,'** Total nodes released in model: ',i5,' **')                  
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                subroutine chk_crack_front_const              *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 5/29/02  rhd               *          
c     *                                                              *          
c     *        This routine checks the crack front for the constant  *          
c     *        front algorithm, where the CTOA is measured a given   *          
c     *        distance behind the crack front.  The algorithm just  *          
c     *        cycles through the master_line lists, evaluating the  *          
c     *        CTOA for each master node.                            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine chk_crack_front_const( killed_this_time, debug, step,          
     &                                  iter )                                  
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : master_lines                                
      use main_data, only :  output_packets, packet_file_no                     
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
      parameter (max_local_list=200)                                            
c                                                                               
      double precision                                                          
     &     d32460, angle, two, one, crit_angle, hundred, dist,                  
     &     targ_dist, local_killed_angles(max_local_list)                       
      dimension local_killed_list(max_local_list)                               
      data d32460, two, one, hundred                                            
     &     / 32460.0, 2.0, 1.0, 100.0/                                          
      logical debug, kill, killed_this_time, use_init                           
c                                                                               
      if( debug ) write(*,*) '>>> in chk_front_master_line'                     
c                                                                               
      killed_this_time = .false.                                                
      num_killed_now   = 0                                                      
c                                                                               
c          loop over the entries in master_lines                                
c                                                                               
      do num_line = 1, num_crack_fronts                                         
c                                                                               
c             get CTOA for each master node                                     
c                                                                               
         if( debug ) write (out,'("    check front:",i7)') num_line             
         call get_slope_master_line( num_line, use_init, angle,                 
     &        dum)                                                              
         if( debug ) write (out,'("     Angle is: ",e13.6)') angle              
c                                                                               
         if( use_init ) then                                                    
            if( debug ) write(out,*) ' using initiation angle'                  
            crit_angle = init_crit_ang * ( one - CTOA_range )                   
         else                                                                   
            if( debug ) write(out,*) ' using crack growth angle'                
            crit_angle = critical_angle * ( one - CTOA_range )                  
         end if                                                                 
c                                                                               
c           if node should be released, then print message and release it       
c                                                                               
c           note that the angle as computed in WARP3D is between the            
c           crack plane and the appropriate location on the crack face.         
c           When we output we need to multiply by two so we get the total       
c           opening angle for the specimen.                                     
c                                                                               
c           To eliminate mesh dependence of the CTOA results, the CTOA is       
c           measured at a given distance behind the crack front.  This          
c           allows an arbitrary number of elements between the crack tip and    
c           the point at which the CTOA is measured.  However, the large        
c           amount of deformation present at the crack tip may cause a number   
c           of nodes to release consecutively after the first.  To avoid this   
c           "spurt growth", we automatically release all the fronts ahead of    
c           the first release until the crack has grown the distance used to    
c           evaluate the CTOA.  This significantly improves results.            
c                                                                               
c                                                                               
         if( angle .ge. crit_angle .or. enforce_node_release ) then             
c                                                                               
c              loop over nodes within critical distance                         
c                                                                               
            if( debug ) write(out,*) ' ->Crit.CTOA reached.Kill node.'          
            loop = 0                                                            
            orig_node = master_lines(num_line,1)                                
c                                                                               
            do while (.true.)                                                   
c                                                                               
c                  set the killed_this_time flag to indicate nodes have been    
c                  released.                                                    
c                                                                               
               loop = loop + 1                                                  
               if( .not. killed_this_time ) then                                
                  write(out,9000)                                               
                  killed_this_time = .true.                                     
               end if                                                           
c                                                                               
c                  print the releasing master nodes                             
c                                                                               
               num_killed_now = num_killed_now + 1                              
               if( num_killed_now .gt. max_local_list ) then                    
                  write(out,9030)                                               
                  call die_gracefully                                           
               end if                                                           
               local_killed_angles(num_killed_now) = angle * two                
c                                                                               
               if( loop .eq. 1 ) then                                           
                  write(out,9010) master_lines(num_line,1), angle * two         
                  local_killed_list(num_killed_now) =                           
     &                 master_lines(num_line,1)                                 
               else                                                             
                  write(out,9015) master_lines(num_line,1), angle * two         
                  local_killed_list(num_killed_now) =                           
     &               -master_lines(num_line,1)                                  
               end if                                                           
                                                                                
c                                                                               
c                  release the front related to the master node, then           
c                  update the crack front.                                      
c                                                                               
               call release_front( master_lines(num_line,1), debug )            
c                                                                               
               call update_crack_front( debug )                                 
c                                                                               
c                  if the master_lines entry is now zero, then the              
c                  crack has coalesed or grown to a free surface:               
c                exit the release loop.                                         
c                                                                               
           if( master_lines(num_line,1) .eq. 0 ) exit                           
c                                                                               
c                  now get the distance between current master node and         
c                  the old master node which initiated the growth increment.    
c                  If releasing the next increment would get us closer to       
c                  the target distance, then continue releasing.  Otherwise     
c                  stop releasing.                                              
c                                                                               
               call get_dist_to_node ( num_line, orig_node, dist )              
c                                                                               
               if( use_init ) then                                              
                  targ_dist = init_ctoa_dist                                    
               else                                                             
                  targ_dist = ctoa_dist                                         
               end if                                                           
               if( abs(dist - targ_dist) .lt.                                   
     &             abs(dist + char_length - targ_dist)) exit                    
c                                                                               
            end do                                                              
c                                                                               
         end if                                                                 
c                                                                               
      end do                                                                    
c                                                                               
c            print critical CTOA if any nodes were released.  Note that         
c            again we multiply the angle by two before we output so that        
c            we are outputting the total critical CTOA.  Also update the        
c            crack front to refect killed elements.                             
c                                                                               
      if( killed_this_time ) then                                               
         if( enforce_node_release ) then                                        
           write(out,9025)                                                      
         else                                                                   
           write(out,9020) critical_angle * two,                                
     &        init_crit_ang * two, CTOA_range * hundred                         
         end if                                                                 
         write(out,9050)  num_ctoa_released_nodes                               
      end if                                                                    
c                                                                               
      if( output_packets .and. killed_this_time ) then                          
        write(packet_file_no) 19, num_killed_now+1, step, iter                  
        write(packet_file_no) critical_angle*two,                               
     &        init_crit_ang*two, CTOA_range*hundred                             
        do i = 1, num_killed_now                                                
          write(packet_file_no) local_killed_list(i),                           
     &         local_killed_angles(i)                                           
        end do                                                                  
      end if                                                                    
c                                                                               
      if ( debug ) write(out,*) '>> leaving chk_front_master_line'              
c                                                                               
c                                                                               
      return                                                                    
 9000 format(/,' >> node release option invoked for the following:',/)          
 9010 format(  '        master node:  ',i7,'  CTOA: ',f6.2,' (degrees)')        
 9015 format(  '        interim node: ',i7,'  CTOA: ',f6.2,' (degrees)')        
 9020 format(/,'        * critical CTOA:',                                      
     &       /,'            for growth    : ',f6.2,' (degrees)',                
     &       /,'            for initiation: ',f6.2,' (degrees)',                
     &       /,'            nodes released within ',f4.1,' % of ',              
     &                        'critical angle.',                                
     &       /)                                                                 
 9025 format(/,'        *** user requested node release ***', / )               
 9030 format(/,'FATAL ERROR: list length exceeded in',                          
     &       /,'             chk_crack_front_const. job aborted.' )             
 9050 format(8x,'** Total nodes released in model: ',i5,' **')                  
c                                                                               
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine release_node                 *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 10/10/95                   *          
c     *                                                              *          
c     *    This routine releases a node -- removes the constraint    *          
c     *    and sets up reaction force release information.           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine release_node( node, node_data_entry, debug )                   
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : crkpln_nodes_state,                         
     &     crkpln_nodes_react                                                   
      use main_data, only : cnstrn_in, cnstrn                                   
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      double precision                                                          
     &     d32460, point_two, point_eight                                       
      data d32460, point_two, point_eight / 32460.0, 0.2, 0.8/                  
      logical debug                                                             
c                                                                               
c                                                                               
      if (debug) write (*,*) '>> releasing node..'                              
c                                                                               
c           To release the node:                                                
c                   a) set global flag that says "no elements                   
c                      have been released" to false.                            
c                                                                               
      no_released_nodes = .false.                                               
c                                                                               
c                                                                               
c                   b) remove constraint from node.  Set the                    
c                      constraint value to null (32460) and remove it           
c                      from the constraint linked list (cstmap).                
c                                                                               
      dof            = dstmap(node)+crk_pln_normal_idx-1                        
      cnstrn(dof)    = d32460                                                   
      cnstrn_in(dof) = d32460                                                   
      cst_ptr        = csthed                                                   
c                                                                               
      if ( cst_ptr .eq. dof ) then                                              
c                                                                               
c                           we are killing the first entry in the list.         
c                           move the list head pointer.                         
c                                                                               
         csthed      = cstmap(dof)                                              
         cstmap(dof) = 0                                                        
      else                                                                      
c                                                                               
c                           cycle through the list to find the entry.  When     
c                           found, make the entry that points to the            
c                           killed entry point to the next valid entry.         
c                                                                               
         do while (.true.)                                                      
            if ( cstmap(cst_ptr) .eq. dof ) then                                
               cstmap(cst_ptr) = cstmap(dof)                                    
               if ( cstmap(cst_ptr) .eq. -1 ) csttail = cst_ptr                 
               cstmap(dof) = 0                                                  
               exit                                                             
            end if                                                              
            cst_ptr = cstmap(cst_ptr)                                           
            if ( cst_ptr .eq. -1 ) then                                         
               write (out,*) '>>>> Fatal Error: failed to find',                
     &              'constraint in constraint linked list.'                     
               call die_gracefully                                              
               stop                                                             
            end if                                                              
         end do                                                                 
      end if                                                                    
c                                                                               
c                  c) Add newly released node to the node release algortihm.    
c                     Store the internal force at the node, and initialize      
c                     the node's state to 1.                                    
c                                                                               
      crkpln_nodes_state(node_data_entry) = 1                                   
      if ( release_type .eq. 1 ) then                                           
        crkpln_nodes_react(node_data_entry) = ifv(dof)                          
      end if                                                                    
      if ( release_type .eq. 2 ) then                                           
        load(dof) = load(dof) - point_two * ifv(dof)                            
        crkpln_nodes_react ( node_data_entry ) = point_eight * ifv(dof)         
      end if                                                                    
c                                                                               
      num_ctoa_released_nodes = num_ctoa_released_nodes + 1                     
c                                                                               
 9999 continue                                                                  
      if (debug) write (*,*) '<< finished releasing node..'                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine release_front                *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 7/31/97                    *          
c     *                                                              *          
c     *    Used in the const_front algorithm. This routine, given a  *          
c     *    master node, finds the attatched crack front, stores it   *          
c     *    for use in the force relaxation algorithm, and releases   *          
c     *    the crack front constraints on each crack front node.     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine release_front( master_node, debug )                            
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : crkpln_nodes_state,                         
     &     crkpln_nodes_react, crack_front_list, inv_crkpln_nodes,              
     &     num_neighbors, neighbor_nodes                                        
      use main_data, only : cnstrn                                              
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
      character(len=1) :: dums                                                  
      real dumr                                                                 
c                                                                               
      double precision                                                          
     &     d32460, dumd                                                         
      data d32460 / 32460.0/                                                    
      logical debug, crack_node, same_front                                     
c                                                                               
c                                                                               
      if (debug) write (*,*) '>> releasing front...'                            
c                                                                               
c           check if master node has already been released. This can happen     
c           if two cracks are just about to coalese, such that there is         
c           a single line of constrained nodes surrounded by free nodes.        
c           In this case, return.                                               
c                                                                               
      dof = dstmap(master_node)+crk_pln_normal_idx-1                            
      if (cnstrn(dof) .eq. d32460) goto 9999                                    
c                                                                               
c           initialize const_front_list for holding the crack front as          
c           we find it.  This will be used by the reaction force                
c           relaxation algorithm.                                               
c                                                                               
      entry = 0                                                                 
      do while (.true.)                                                         
         entry = entry + 1                                                      
c                                                                               
c               if no free space in the const_front_list                        
c               is found, then stop execution                                   
c                                                                               
         if (entry .gt. num_crack_fronts * num_nodes_grwinc) then               
            call errmsg (285, dum, dums, dumr, dumd)                            
            call die_gracefully                                                 
            stop                                                                
         endif                                                                  
c                                                                               
c               free space found; start the list                                
c                                                                               
         if (crack_front_list(entry,1) .eq. 0) then                             
            list_entry = entry                                                  
            crack_front_list(list_entry,1) = master_node                        
            exit                                                                
         endif                                                                  
      enddo                                                                     
c                                                                               
c           loop over the crack front node list.  For each node, release        
c           it, then check each of its neighbors to find out if they are        
c           nodes on the crack front.  If a neighbor is a crackfront node,      
c           and it is not on the crack_front_node list, then add it to the      
c           list.  In this way, the crack_front_node list will be both          
c           processed such that each node on it is released, and the list       
c           will be expanded as new crack front nodes are found.                
c                                                                               
      if (debug) write (out,*) '> check crack_front_list'                       
      ptr = 0                                                                   
      ptr_last = 1                                                              
      do while (.true.)                                                         
         ptr = ptr + 1                                                          
         if (debug) write (out,*) '  - ptr = ', ptr                             
c                                                                               
c              check if we have processed all the nodes on the list; exit       
c              if we have.                                                      
c                                                                               
         if (ptr .gt. num_nodes_thick) exit                                     
         node = crack_front_list ( entry, ptr )                                 
         if (ptr .gt. ptr_last) exit                                            
         if (node .eq. 0) exit                                                  
         node_data_entry = inv_crkpln_nodes ( node )                            
c                                                                               
c              node has not yet been processed.  release it.                    
c                                                                               
         call release_node ( node, node_data_entry, debug )                     
         if (node .ne. master_node) write (out,9000) node,master_node           
c                                                                               
c              now search the node's neighbors, looking for crack front         
c              nodes.  A neighbor node is a crack front node if it is           
c              a constrained node and one of it's neighbors is an               
c              unconstrained node.  Once neighbor nodes have been               
c              located, add them to crack_front_list if they are not            
c              already in the list.                                             
c                                                                               
         do i = 1, num_neighbors(node_data_entry)                               
c                                                                               
            neighbor_node = neighbor_nodes (i,node_data_entry)                  
            dof = dstmap(neighbor_node)+crk_pln_normal_idx-1                    
            if (cnstrn(dof) .eq. d32460) cycle                                  
c                                                                               
c                       check if node is already assigned in                    
c                       const_front_list; skip it if it is                      
c                                                                               
            do j = 1, ptr_last                                                  
               if ( crack_front_list(list_entry, j)                             
     &             .eq. neighbor_node ) go to 100                               
            end do                                                              
c                                                                               
c                       check the neighbors of the neighbors                    
c                                                                               
            crack_node = .false.                                                
            neighbor_data_entry = inv_crkpln_nodes(neighbor_node)               
c                                                                               
            do j = 1, num_neighbors(neighbor_data_entry)                        
               neigh2 = neighbor_nodes (j,neighbor_data_entry)                  
               dof2   = dstmap(neigh2)+crk_pln_normal_idx-1                     
               if ( cnstrn(dof2) .ne. d32460 ) cycle                            
c                                                                               
c                          check if node was released just this step            
c                                                                               
               do k = 1, ptr_last                                               
                  if ( crack_front_list(list_entry, k)                          
     &                .eq. neighbor_node ) cycle                                
               end do                                                           
c                                                                               
c                          skip node if it is the node we are searching         
c                          from, or if node found is not on same crack          
c                          front.                                               
c                                                                               
               if (neigh2 .eq. node) cycle                                      
           if (.not. same_front(neighbor_node, neighbor_data_entry,             
     &               node) ) cycle                                              
c                                                                               
               crack_node = .true.                                              
               exit                                                             
            end do                                                              
c                                                                               
c                       if node is a crack node, then put it in                 
c                       crack_front_list                                        
c                                                                               
            if (crack_node) then                                                
               ptr_last = ptr_last + 1                                          
               if (debug) write (out,*) '     cfn neighbor:',                   
     &                neighbor_node,' into ',ptr_last                           
               if (ptr_last .gt. num_nodes_thick) then                          
                  call errmsg (285, dum, dums, dumr, dumd)                      
                  call die_gracefully                                           
                  stop                                                          
               endif                                                            
               crack_front_list(list_entry,ptr_last) = neighbor_node            
            endif                                                               
c                                                                               
 100        continue                                                            
         end do                                                                 
c                                                                               
c              loop over to next entry in crack_front_list                      
c                                                                               
      end do                                                                    
c                                                                               
c              if we are using release by steps, then we don't a                
c              need the nodal list anymore.  Clear it out.                      
c                                                                               
      if (release_type .eq. 1) then                                             
         do i=1, num_nodes_thick                                                
            crack_front_list(list_entry,i) = 0                                  
         end do                                                                 
      endif                                                                     
c                                                                               
c                                                                               
 9999 continue                                                                  
      if (debug) write (*,*) '<< finished releasing front..'                    
c                                                                               
      return                                                                    
 9000 format(  '             node: ',i7,'  slave of master node: ',i6)          
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine update_crack_front           *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 7/29/97                    *          
c     *                                                              *          
c     *      This routine updates the crack front list, removing     *          
c     *      entries of newly released nodes and adding to the       *          
c     *      to the list their neighbors that are still constrained  *          
c     *      in the direction normal to the crack plane.             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine update_crack_front( debug )                                    
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : inv_crkpln_nodes, num_neighbors,            
     &     neighbor_nodes, crack_front_nodes, crkpln_nodes_state,               
     &     master_nodes                                                         
      use main_data, only : cnstrn                                              
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      character(len=1) :: dums                                                  
      real dumr                                                                 
      double precision                                                          
     &     d32460, dumd, zero                                                   
      data d32460, zero  / 32460.0, 0.0 /                                       
      logical debug, next_node_found                                            
c                                                                               
      integer master                                                            
c                                                                               
c                                                                               
c            Traverse linked list of the crack front nodes.                     
c            If we find a node that has been newly killed:                      
c                   a) delete the node from the crack front list                
c                   b) add any neighbors of this node that are still            
c                      constrained to the crack front list                      
c                                                                               
      if (debug) write (out,*) '>>>>>> Updating linked list.'                   
c                                                                               
c            initialize the pointers for traversing the list.                   
c                                                                               
      pointer = crack_front_start                                               
      above_node_ptr = -1                                                       
c                                                                               
c            top of list traversal loop                                         
c                                                                               
 10   continue                                                                  
c                                                                               
c            a node is newly killed if its state is 1.  get the                 
c            state and skip node if the state is not 1.                         
c                                                                               
      node = crack_front_nodes(pointer,1)                                       
      node_data_entry = inv_crkpln_nodes(node)                                  
      if (crkpln_nodes_state(node_data_entry).ne.1) goto 20                     
c                                                                               
c                  node is newly killed.                                        
c                                                                               
c                    a) remove it from the list                                 
c                                                                               
c                       Note: by doing this, we won't need to move              
c                             the pointers for the next list entry              
c                                                                               
      temp_ptr = crack_front_nodes(pointer,2)                                   
      call rm_from_list (pointer, above_node_ptr)                               
      pointer = temp_ptr                                                        
c                                                                               
c                    b) add any still constrained neighbors to                  
c                       crack_front_nodes                                       
c                                                                               
      next_node_found = .false.                                                 
      do neighbor = 1, num_neighbors(node_data_entry)                           
c                                                                               
         neighbor_node = neighbor_nodes(neighbor,node_data_entry)               
         dof= dstmap(neighbor_node)+crk_pln_normal_idx-1                        
         if (cnstrn(dof).ne.d32460) then                                        
c                                                                               
c                         neighbor is constrained. check if it is               
c                         already on the list.  If it isn't, add to             
c                         the list.                                             
c                                                                               
c                         if we are doing const_growth, then replace the        
c                         old master node in the list of master nodes with      
c                         the new node.                                         
c                                                                               
            new_node = neighbor_node                                            
            next_node_found = .true.                                            
           if (find_in_list(neighbor_node,dumi).eq.-1) then                     
               call add_to_list (neighbor_node)                                 
               master_entry = master(node)                                      
               if (master_entry .gt. 0) then                                    
                  master_nodes(master_entry) = neighbor_node                    
               endif                                                            
            endif                                                               
         end if                                                                 
      end do                                                                    
c                                                                               
c           If we are using the constant front algorithm:                       
c            Since we have removed a node, we need to update the master         
c            line list.  If we removed a node, but didn't add a new             
c            one, then a crack has coalesed, so we zero out the                 
c            corresponding entry.                                               
c                                                                               
      if ( const_front ) call update_master_line (node, new_node,               
     &     next_node_found, debug )                                             
c                                                                               
      goto 30                                                                   
c                                                                               
c            node was not released in this step, so we need to                  
c            move the pointers to the next entry in the list.                   
c                                                                               
 20   continue                                                                  
      above_node_ptr = pointer                                                  
      pointer = crack_front_nodes(pointer,2)                                    
c                                                                               
c            if we are not at end of list, check the next node.                 
c                                                                               
 30   continue                                                                  
      if (pointer .ne. -1) goto 10                                              
c                                                                               
      if (debug) write (out,*) '>>>>>> Linked list is now updated.'             
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine update_master_line           *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 10/16/97                   *          
c     *                                                              *          
c     *        This routine updates the master line list after       *          
c     *        all the crack front releases have been calculated.    *          
c     *        If the node in front of the old master node is        *          
c     *        unconstrained, then a crack has coalesed. In this     *          
c     *        case we zero out the entry in the master_line list    *          
c     *        to indicate the crack front no longer exists.         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine update_master_line( old_node, new_node,                        
     &     next_node_found, debug )                                             
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : master_lines, old_angles_at_front           
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
      double precision                                                          
     &     zero                                                                 
      data zero  / 0.0 /                                                        
c                                                                               
      logical next_node_found, debug                                            
c                                                                               
c          if no replacement node was found, then the crack is coalesing        
c          or ending.  Zero out the corresponding entry in the master_lines     
c          structure.                                                           
c                                                                               
      if ( .not. next_node_found ) then                                         
c                                                                               
         do i = 1, num_crack_fronts                                             
            if ( master_lines(i,1) .eq. old_node ) then                         
               if (debug) write (*,*) '-> dead node, zero master_line:',        
     &              old_node                                                    
               do j = 1, num_nodes_back + 1                                     
                  master_lines(i,j) = 0                                         
               enddo                                                            
            endif                                                               
         enddo                                                                  
c                                                                               
c           if the crack is growing normally, then new_node holds the           
c           new crack front node.  move all of the entries in the               
c           master_line list over one to make room for the new crack            
c           tip node.                                                           
c                                                                               
      else                                                                      
c                                                                               
         do num_line = 1, num_crack_fronts                                      
            if ( master_lines(num_line,1) .eq. old_node ) then                  
               do i = num_nodes_back,1,-1                                       
                  master_lines( num_line, i+1 ) =                               
     &                 master_lines( num_line, i )                              
               enddo                                                            
               exit                                                             
            endif                                                               
         enddo                                                                  
         master_lines( num_line, 1 ) = new_node                                 
c                                                                               
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine get_slope                    *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 08/28/95                   *          
c     *                                                              *          
c     *        This routine gets the angle between a crack front     *          
c     *        node and an unconstrained neighboring crack face      *          
c     *        node.  This is used in the node_release algorithm.    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine get_slope( node1, node2, normal, angle )                       
      use global_data ! old common.main
c                                                                               
      use main_data, only : crdmap                                              
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      double precision                                                          
     &      angle, dist_x, dist_y, dist_z, plane_dist, pi                       
      data pi /3.14159/                                                         
c                                                                               
c               calculate the distances in the three coordinate                 
c               directions between the two nodes.                               
c                                                                               
      dist_x = c(crdmap(node1)) + u(dstmap(node1))                              
     &      - ( c(crdmap(node2)) + u(dstmap(node2)))                            
      dist_y = c(crdmap(node1)+1) + u(dstmap(node1)+1)                          
     &      - ( c(crdmap(node2)+1) + u(dstmap(node2)+1))                        
      dist_z = c(crdmap(node1)+2) + u(dstmap(node1)+2)                          
     &      - ( c(crdmap(node2)+2) + u(dstmap(node2)+2))                        
c                                                                               
      if ( normal .eq. 1 ) then                                                 
c                                                                               
c               x is direction of crack plane normal                            
c                                                                               
         plane_dist = sqrt (dist_y ** 2 + dist_z ** 2)                          
         angle = atan( dist_x / plane_dist ) * 180.0 / pi                       
c                                                                               
      else if (normal .eq. 2) then                                              
c                                                                               
c               y is direction of crack plane normal                            
c                                                                               
         plane_dist = sqrt (dist_x ** 2 + dist_z ** 2)                          
         angle = atan( dist_y / plane_dist ) * 180.0 / pi                       
c                                                                               
      else if ( normal .eq. 3 ) then                                            
c                                                                               
c               z is direction of crack plane normal                            
c                                                                               
         plane_dist = sqrt (dist_x ** 2 + dist_y ** 2)                          
         angle = atan( dist_z / plane_dist ) * 180.0 / pi                       
c                                                                               
      end if                                                                    
c                                                                               
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine get_slope_master_line        *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 10/16/97                   *          
c     *                                                              *          
c     *        This routine gets the CTOA at a master node at a      *          
c     *        defined distance behind the crack front, as defined   *          
c     *        by the master_line for the node.                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine get_slope_master_line( num_line, use_init, angle,              
     &     num_elems_back)                                                      
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : master_lines, crkpln_nodes_state,           
     &     inv_crkpln_nodes                                                     
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      double precision                                                          
     &     angle, dist, d_dist, height, pi, zero, factor, d_height,             
     &     one_eighty, one, four, ctoa_dist_back, dumd                          
      data zero, one_eighty, one, four / 0.0, 180.0, 1.0, 4.0/                  
      real dumr                                                                 
      character(len=1) :: dums                                                  
      logical debug, use_init                                                   
c                                                                               
      debug = .false.                                                           
      pi = four * atan(one)                                                     
      num_elems_back = 0                                                        
c                                                                               
c           if master_lines entry is zeroed, then set angle to zero             
c           and leave.                                                          
c                                                                               
      if (master_lines(num_line,1) .eq. 0) then                                 
         angle = zero                                                           
         return                                                                 
      endif                                                                     
c             determine if we should use the initial angle or continuation      
c             angle. If the second node in the master_line list has never       
c             been released, then the crack has not yet moved, so we can        
c             use the initiation angle and distance.                            
c                                                                               
      if ( crkpln_nodes_state( inv_crkpln_nodes(                                
     &     master_lines( num_line, 2 ) ) ) .eq.0 ) then                         
         use_init = .true.                                                      
         ctoa_dist_back = init_ctoa_dist                                        
         if (debug) write (*,*) 'Using initiation angle'                        
      else                                                                      
         use_init = .false.                                                     
         ctoa_dist_back = ctoa_dist                                             
         if (debug) write (*,*) 'Using release angle'                           
      endif                                                                     
c                                                                               
c           calculate arc distance along the master line until we               
c           find the nodes which bound our distance, and interpolate.           
c                                                                               
      dist = zero                                                               
      height = zero                                                             
      base_node = master_lines(num_line,1)                                      
      next_node = master_lines(num_line,2)                                      
      node_idx = 2                                                              
      old_height = zero                                                         
c                                                                               
c               travel along the master line                                    
c                                                                               
      do while (.true.)                                                         
c                                                                               
         num_elems_back = num_elems_back + 1                                    
c                                                                               
c                  get distance                                                 
c                                                                               
         call get_dist ( next_node, base_node, crk_pln_normal_idx,              
     &        d_dist, d_height)                                                 
c                                                                               
c                  check distance -- if further than we want to go,             
c                  calculate the angle                                          
c                                                                               
         if (dist + d_dist .ge. ctoa_dist_back) then                            
c                                                                               
            factor = ( ctoa_dist_back - dist ) / d_dist                         
            dist = dist + d_dist * factor                                       
            height = height + d_height * factor                                 
            angle = atan( height / dist ) * one_eighty / pi                     
            exit                                                                
c                                                                               
         else                                                                   
c                                                                               
c                  we haven't passed our distance yet.  Continue travelling     
c                  along the master line.                                       
c                                                                               
            dist = dist + d_dist                                                
            height = height + d_height                                          
            node_idx = node_idx + 1                                             
            if (node_idx .gt. num_nodes_back + 1) then                          
               call errmsg( 290, num_nodes_back, dums, dumr, dumd)              
               call die_gracefully                                              
               stop                                                             
            endif                                                               
            base_node = next_node                                               
            next_node = master_lines(num_line,node_idx)                         
            if (next_node .eq. 0) then                                          
               call errmsg( 291, master_lines(num_line,1), dums, dumr,          
     &              dumd)                                                       
               call die_gracefully                                              
               stop                                                             
            endif                                                               
c                                                                               
         endif                                                                  
c                                                                               
      end do                                                                    
      if (debug) then                                                           
         write (*,*) ' nodes;',base_node, next_node                             
         write (*,*) ' number of elements back;',num_elems_back                 
         write (*,*) ' dist:',dist, '  height:',height                          
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
                                                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine get_dist_to_node             *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 5/5/98                     *          
c     *                                                              *          
c     *        This routine gets the distance from a given old       *          
c     *        master node to the current master node.               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine get_dist_to_node( num_line, old_node, dist)                    
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : master_lines                                
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      double precision                                                          
     &     dist, d_dist, height, zero, factor, d_height,                        
     &     one_eighty, one, four, ctoa_dist_back, dumd                          
      data zero, one_eighty, one, four / 0.0, 180.0, 1.0, 4.0/                  
      real dumr                                                                 
      character(len=1) :: dums                                                  
      logical debug, use_init                                                   
c                                                                               
      debug = .false.                                                           
c                                                                               
c           if master_lines entry is zeroed, then set distance to zero          
c           and leave.                                                          
c                                                                               
      if (master_lines(num_line,1) .eq. 0) then                                 
         dist = zero                                                            
         return                                                                 
      endif                                                                     
c                                                                               
c           calculate arc distance along the master line until we               
c           find the nodes which bound our distance, and interpolate.           
c                                                                               
      dist = zero                                                               
      base_node = master_lines(num_line,1)                                      
      next_node = master_lines(num_line,2)                                      
      node_idx = 2                                                              
c                                                                               
c               travel along the master line                                    
c                                                                               
      do while (.true.)                                                         
c                                                                               
c                  get distance                                                 
c                                                                               
         call get_dist ( next_node, base_node, crk_pln_normal_idx,              
     &        d_dist, d_height)                                                 
c                                                                               
c                  check distance -- if further than we want to go,             
c                  calculate the angle                                          
c                                                                               
         dist = dist + d_dist                                                   
         node_idx = node_idx + 1                                                
         if (node_idx .gt. num_nodes_back + 1) then                             
            call errmsg( 290, num_nodes_back, dums, dumr, dumd)                 
            call die_gracefully                                                 
            stop                                                                
         endif                                                                  
         base_node = next_node                                                  
         next_node = master_lines(num_line,node_idx)                            
         if ( base_node .eq. old_node ) then                                    
            exit                                                                
         else if (next_node .eq. 0) then                                        
            call errmsg( 291, master_lines(num_line,1), dums, dumr,             
     &           dumd)                                                          
            call die_gracefully                                                 
            stop                                                                
         endif                                                                  
c                                                                               
      end do                                                                    
c                                                                               
      if (debug) then                                                           
         write (*,*) '>>>>> in get_dist_to_node: dist to old node:',dist        
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine get_dist                     *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 10/16/97                   *          
c     *                                                              *          
c     *        This routine calculates the distance between two      *          
c     *        nodes originally on the crack plane.  Distance is     *          
c     *        given in terms of planar distance and height.         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine get_dist( node1, node2, normal, plane_dist, height)            
      use global_data ! old common.main
c                                                                               
      use main_data, only : crdmap                                              
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      double precision                                                          
     &      height, dist_x, dist_y, dist_z, plane_dist                          
c                                                                               
c               calculate the distances in the three coordinate                 
c               directions between the two nodes.                               
c                                                                               
      dist_x = c(crdmap(node1)) + u(dstmap(node1))                              
     &      - ( c(crdmap(node2)) + u(dstmap(node2)))                            
      dist_y = c(crdmap(node1)+1) + u(dstmap(node1)+1)                          
     &      - ( c(crdmap(node2)+1) + u(dstmap(node2)+1))                        
      dist_z = c(crdmap(node1)+2) + u(dstmap(node1)+2)                          
     &      - ( c(crdmap(node2)+2) + u(dstmap(node2)+2))                        
c                                                                               
      if ( normal .eq. 1 ) then                                                 
c                                                                               
c               x is direction of crack plane normal                            
c                                                                               
         plane_dist = sqrt (dist_y ** 2 + dist_z ** 2)                          
         height = dist_x                                                        
c                                                                               
      else if (normal .eq. 2) then                                              
c                                                                               
c               y is direction of crack plane normal                            
c                                                                               
         plane_dist = sqrt (dist_x ** 2 + dist_z ** 2)                          
         height = dist_y                                                        
c                                                                               
      else if ( normal .eq. 3 ) then                                            
c                                                                               
c               z is direction of crack plane normal                            
c                                                                               
         plane_dist = sqrt (dist_x ** 2 + dist_y ** 2)                          
         height = dist_z                                                        
c                                                                               
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      function same_front                     *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 7/29/97                    *          
c     *                                                              *          
c     *   This function, given two nodes on crack fronts, returns    *          
c     *   as true if the nodes are neighbors on the same crack front.*          
c     *   This function is used in the const_front algorithm.        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      logical function same_front (new_node, new_data_entry,                    
     &    check_node )                                                          
c       
      use global_data                                                                         
      use node_release_data, only : inv_crkpln_nodes, num_neighbors,            
     &     neighbor_nodes                                                       
      use main_data, only : cnstrn, cnstrn_in                                   
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      double precision                                                          
     &     d32460                                                               
      data d32460  / 32460.0 /                                                  
c                                                                               
c           This function checks a newly found crack front node to see if       
c           it is on the same front as the node from which the new node         
c           was found.  This becomes an issue when two crack fronts are         
c           back to back; in this case, the crack front location algorithm      
c           can break down, finding both crack fronts.                          
c                                                                               
c           The basic idea behind this routine is to loop back over the         
c           unconstrained nodes adjoining the crack front to make sure          
c           that there is a path from the new_node to the check_node            
c           consisting of precisely two unconstrained nodes.                    
c                                                                               
c           This routine operates as follows:                                   
c                                                                               
c           - find an unconstrained neighbor to the new node.                   
c           - loop over this neighbor's neighbors.  for each that is            
c             unconstrained:                                                    
c                - check the constrained neighbors of this new node.  If        
c                  one is the check_node, then new_node and check_node          
c                  are on the same crack front.                                 
c                                                                               
      same_front = .false.                                                      
      do i = 1, num_neighbors(new_data_entry)                                   
         node1 = neighbor_nodes(i,new_data_entry)                               
         if( node1 .eq. check_node ) cycle                                      
         node1_data_entry = inv_crkpln_nodes(node1)                             
         dof = dstmap(node1)+crk_pln_normal_idx-1                               
         if( cnstrn(dof) .ne. d32460 ) cycle                                    
         do j = 1, num_neighbors(node1_data_entry)                              
            node2 = neighbor_nodes(j,node1_data_entry)                          
            node2_data_entry = inv_crkpln_nodes(node2)                          
            dof = dstmap(node2)+crk_pln_normal_idx-1                            
            if( cnstrn(dof) .ne. d32460 ) cycle                                 
            do k = 1, num_neighbors(node2_data_entry)                           
               node3 = neighbor_nodes(k,node2_data_entry)                       
               if( node3 .eq. check_node ) then                                 
                  same_front = .true.                                           
                  go to 9999                                                    
               end if                                                           
            end do                                                              
         end do                                                                 
      end do                                                                    
c                                                                               
 9999 continue                                                                  
      return                                                                    
      end                                                                       
