c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine overshoot_CTOA                  *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 4/8/96                     *          
c     *                                                              *          
c     *        This routine traverses the nodes on the current       *          
c     *        crack front, extrapolating the CTOA for each angle    *          
c     *        between the crack front and it's neighbors. If any    *          
c     *        of the extrapolated angles are larger than            *          
c     *        overshoot_limit of of the critical release CTOA,      *          
c     *        reduce the size of the load step to come closer       *          
c     *        to the release CTOA.                                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine overshoot_CTOA ( new_load_fact, mf, mf_nm1, debug )            
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : inv_crkpln_nodes, num_neighbors,            
     &     neighbor_nodes, crack_front_nodes, crkpln_nodes_state                
      use main_data, only : cnstrn                                              
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      double precision                                                          
     &     d32460, angle, two, new_load_fact, crit_angle, one, mf,              
     &     mf_nm1, zero, dtemp                                                  
      character(len=1) :: dums                                                  
      data d32460, two, one, zero / 32460.0, 2.0, 1.0, 0.0 /                    
      logical debug                                                             
c                                                                               
c          If this step starts a new loading condition, then we do              
c          not have enough information to change the load                       
c          step size; return to calling routine.                                
c          If there are no more crack plane                                     
c         nodes, then return to calling routine.                                
c                                                                               
      if ( debug ) write (out,*) '>>>>>> in overshoot_control'                  
      if (mf_nm1 .eq. zero) then                                                
         if ( debug ) write (out,*) ' new loading condition; skip'              
         goto 9999                                                              
      endif                                                                     
c                                                                               
c          traverse the linked list of crack front nodes and                    
c          check angles.                                                        
c                                                                               
c             set pointer to top of list, and the previous node                 
c             pointer to null.                                                  
c                                                                               
      prev_node_ptr = -1                                                        
      node_ptr = crack_front_start                                              
      if (node_ptr .eq. -1) goto 9999                                           
c                                                                               
c             traverse the list of crack front nodes.  Extrapolate the angles   
c             between the crack front node and each of its neighbors for the    
c             forthcoming step.  If prediction indicates the angle              
c             will be more than overshoot_limit of the critical angle,          
c             set a factor which reduces the load step size to an               
c             acceptible level.                                                 
c                                                                               
 10   continue                                                                  
c                                                                               
      node = crack_front_nodes( node_ptr, 1 )                                   
      node_data_entry = inv_crkpln_nodes( node )                                
c                                                                               
      if ( debug ) write (out,*) ' >> Checking node ', node                     
c                                                                               
c                loop over neighbor nodes                                       
c                                                                               
      do neighbor = 1, num_neighbors( node_data_entry )                         
c                                                                               
c                if neighbor is constrained in the direction normal             
c                to the crack plane, skip it                                    
c                                                                               
         neighbor_node = neighbor_nodes( neighbor, node_data_entry )            
         dof= dstmap( neighbor_node ) + crk_pln_normal_idx - 1                  
         if ( cnstrn(dof).ne.d32460 ) cycle                                     
         if ( debug ) write (out,'("     neighbor is: ",i6)')                   
     &        neighbor_node                                                     
c                                                                               
c                find angle and appropriate critical angle ( initial            
c                growth or continued growth )                                   
c                                                                               
         call get_slope( neighbor_node, node, crk_pln_normal_idx,               
     &        angle )                                                           
         if ( crkpln_nodes_state( inv_crkpln_nodes( neighbor_node ) )           
     &        .eq.0 ) then                                                      
            crit_angle = init_crit_ang                                          
         else                                                                   
            crit_angle = critical_angle                                         
         endif                                                                  
         if ( debug ) write (out,'("     Angle is: ",f6.3)') angle              
     &        * two                                                             
c                                                                               
c                if current angle is less than critical angle, check            
c                for potential overshoot in upcoming step                       
c                                                                               
         if ( angle .lt. crit_angle * ( one - CTOA_range ) )                    
     &        call check_CTOA_overshoot (                                       
     &        node_ptr, neighbor, angle, crit_angle, new_load_fact,             
     &        mf, mf_nm1, debug )                                               
c                                                                               
      end do                                                                    
c                                                                               
c            move to the next entry in the list, unless we are at the end.      
c                                                                               
      prev_node_ptr = node_ptr                                                  
      node_ptr = crack_front_nodes( node_ptr, 2 )                               
      if ( node_ptr .ne. -1 ) goto 10                                           
c                                                                               
c            if load factor has been changed, alert the user                    
c                                                                               
      if ( new_load_fact .ne. one ) then                                        
        dtemp =  max(new_load_fact, min_load_fact)                              
        call errmsg ( 270, dum, dums, real(overshoot_limit), dtemp )            
      end if                                                                    
c                                                                               
c                                                                               
 9999 continue                                                                  
      if ( debug ) write(out,*) '<<<<< leaving overshoot_reduction'             
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine over_CTOA_const                 *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 4/8/96                     *          
c     *                                                              *          
c     *        This routine traverses the nodes on the current       *          
c     *        crack front, extrapolating the CTOA for each angle    *          
c     *        between the crack front and it's neighbors. If any    *          
c     *        of the extrapolated angles are larger than            *          
c     *        overshoot_limit of of the critical release CTOA,      *          
c     *        reduce the size of the load step to come closer       *          
c     *        to the release CTOA.                                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine over_CTOA_const ( new_load_fact, mf, mf_nm1, debug )           
      use global_data ! old common.main
      use damage_data                                                           
      implicit integer (a-z)                                                    
c                                                                               
      double precision                                                          
     &     d32460, angle, two, new_load_fact, crit_angle, one, mf,              
     &     mf_nm1, zero, dtemp, dtol                                            
      character(len=1) :: dums                                                  
      data d32460, two, one, zero, dtol / 32460.0, 2.0, 1.0, 0.0, 0.001/        
      logical debug, use_init                                                   
c                                                                               
      if ( debug ) write(out,*) '>>>>> entering over_CTOA_master'               
c                                                                               
c         loop over master_lines                                                
c                                                                               
      do num_line = 1, num_crack_fronts                                         
c                                                                               
         if (debug) write (out,'("    check front:",i7)') num_line              
         call get_slope_master_line( num_line, use_init, angle, idummy )        
         if (angle .le. dtol) cycle                                             
c                                                                               
         if ( use_init ) then                                                   
            crit_angle = init_crit_ang                                          
         else                                                                   
            crit_angle = critical_angle                                         
         endif                                                                  
         if ( debug ) write (out,'("     Angle is: ",f6.3)') angle              
     &        * two                                                             
c                                                                               
c                if current angle is less than critical angle, check            
c                for potential overshoot in upcoming step                       
c                                                                               
         if ( angle .lt. crit_angle * ( one - CTOA_range ) )                    
     &        call check_CTOA_overshoot (                                       
     &        num_line, 1, angle, crit_angle, new_load_fact,                    
     &        mf, mf_nm1, debug )                                               
c                                                                               
      end do                                                                    
c                                                                               
c            if load factor has been changed, alert the user                    
c                                                                               
      if ( new_load_fact .ne. one ) then                                        
        dtemp =  max(new_load_fact, min_load_fact)                              
        call errmsg ( 270, dum, dums, real(overshoot_limit), dtemp )            
      end if                                                                    
c                                                                               
c                                                                               
 9999 continue                                                                  
      if ( debug ) write(out,*) '<<<<< leaving over_CTOA_master'                
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                subroutine check_CTOA_overshoot               *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 04/05/96                   *          
c     *                                                              *          
c     *        This routine finds the extrapolated opening angle     *          
c     *        for given node and neighbor for next step.  If this   *          
c     *        is more than allowed overshoot_limit of the           *          
c     *        critical release angle, then it calculates the        *          
c     *        reduction in load step size required to release the   *          
c     *        node closer to the critical angle.                    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine check_CTOA_overshoot ( node_ptr, neighbor, angle,              
     &     crit_angle, new_load_fact, mf, mf_nm1, debug )                       
      use global_data ! old common.main
      use node_release_data, only : old_angles_at_front                         
      use damage_data                                                           
      implicit integer (a-z)                                                    
c                                                                               
      double precision                                                          
     &     angle, crit_angle, new_load_fact, mf, mf_nm1,                        
     &     new_angle, load_fact, two, one                                       
      data two, one /2.0, 1.0/                                                  
      logical debug                                                             
c                                                                               
c             find projected angle. If last load step was reduced by            
c             the step reduction algorithm, then divide change in angle by      
c             the load factor to approxiate angle as if the load step           
c             had not been reduced.  If multiplication factors                  
c             have changed over the step, include change in estimate.           
c             If load factor has been reduced permanently, reduce the           
c             estimate by the same factor.                                      
c                                                                               
      new_angle = angle + ( angle - old_angles_at_front( node_ptr,              
     &     neighbor ) ) / control_load_fact  * mf / mf_nm1                      
     &     * perm_load_fact                                                     
      if ( debug ) write (out,'("     prediction:",f6.3)')                      
     &     new_angle * two                                                      
c                                                                               
c             if projected angle is larger than the max accepted angle          
c             (crit_angle * (1 + % allowed error)), calculate the new           
c             load factor term -- if it is smaller than the one calculated      
c             so far, store the new one.  In calculating the load factor,       
c             include effects of permanent load size reductions and             
c             changes in the multiplication factor.                             
c                                                                               
      if ( new_angle .gt. ( ( one + overshoot_limit ) *                         
     &     crit_angle ) ) then                                                  
         load_fact = ( crit_angle - angle ) *  control_load_fact /              
     &        ( angle - old_angles_at_front( node_ptr, neighbor ) )             
     &        * ( mf_nm1 / mf ) / perm_load_fact                                
         if ( debug ) write (out,'("     new load fact:",e13.6)')               
     &        load_fact                                                         
         new_load_fact = min( new_load_fact, load_fact )                        
      endif                                                                     
c                                                                               
c             store current angle for use in next step load reduction           
c                                                                               
      old_angles_at_front( node_ptr, neighbor ) = angle                         
c                                                                               
c                                                                               
 9999 continue                                                                  
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine CTOA_cut_step                   *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 4/8/96                     *          
c     *                                                              *          
c     *         This routine checks if the load step size is too     *          
c     *         large. We find all neighbors of crack front nodes    *          
c     *         which produce CTOAs larger than the critical angle.  *          
c     *         If the number of load steps since release of one     *          
c     *         of these neighbors is smaller than the minimum       *          
c     *         number ( min_steps_for_release ), then cut the load  *          
c     *         step size by half.                                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine CTOA_cut_step ( debug )                                        
      use global_data ! old common.main
c                                                                               
      use node_release_data, only : inv_crkpln_nodes, num_neighbors,            
     &     neighbor_nodes, crack_front_nodes, crkpln_nodes_state                
      use damage_data                                                           
      use main_data, only : cnstrn                                              
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      double precision                                                          
     &     d32460, angle, two, crit_angle, one                                  
      character(len=1) :: dums                                                  
      real dumr                                                                 
      data d32460, two, one / 32460.0, 2.0, 1.0 /                               
      logical debug                                                             
c                                                                               
c                                                                               
      if ( debug ) write (out,*) '>>>>>> in CTOA_cut_step'                      
c                                                                               
c          traverse the linked list of crack front nodes and                    
c          check angles.                                                        
c                                                                               
c             set pointer to top of list, and the previous node                 
c             pointer to null.  If there are no more crack plane                
c            nodes, then skip rest of routine.                                  
c                                                                               
      prev_node_ptr = -1                                                        
      node_ptr = crack_front_start                                              
      if (node_ptr .eq. -1) goto 9999                                           
c                                                                               
c             traverse the list of crack front nodes. If node is about          
c             to be released, check each neighbor who creates angle             
c             larger than the critical angle, and find how many steps           
c             since the neighbor was released.  If less than                    
c             min_steps_for_release, cut step size in half.                     
c                                                                               
 10   continue                                                                  
c                                                                               
      node = crack_front_nodes( node_ptr, 1 )                                   
      node_data_entry = inv_crkpln_nodes( node )                                
c                                                                               
c                loop over neighbor nodes                                       
c                                                                               
      do neighbor = 1, num_neighbors( node_data_entry )                         
c                                                                               
c                if neighbor is constrained in the direction normal             
c                to the crack plane, skip it                                    
c                                                                               
         neighbor_node = neighbor_nodes( neighbor, node_data_entry )            
         dof = dstmap( neighbor_node ) + crk_pln_normal_idx - 1                 
         if ( cnstrn(dof).ne.d32460 ) cycle                                     
c                                                                               
c                find angle and appropriate critical angle ( initial            
c                growth or continued growth ).  See if node will be             
c                released.                                                      
c                                                                               
         call get_slope( neighbor_node, node, crk_pln_normal_idx,               
     &        angle )                                                           
         if ( crkpln_nodes_state( inv_crkpln_nodes( neighbor_node ) )           
     &        .eq.0 ) then                                                      
            crit_angle = init_crit_ang * (one - CTOA_range)                     
         else                                                                   
            crit_angle = critical_angle * (one - CTOA_range)                    
         end if                                                                 
c                                                                               
c                if current angle is larger than critical angle, then           
c                node will be released this step. Check                         
c                if number of steps since release of neighbor is less           
c                the min_steps_for_release.                                     
c                                                                               
         if ( angle .gt. crit_angle )then                                       
            neighbor_state = crkpln_nodes_state( inv_crkpln_nodes               
     &                       (neighbor_node) )                                  
            if ( neighbor_state.gt.0 .and.                                      
     &           neighbor_state.lt.min_steps_for_release) then                  
c                                                                               
c                    too few steps between release. cut rate in                 
c                    half.  Note that we change the perm_load_fact              
c                    variable -- unlike the temp_load_fact, which               
c                    changes the load step size for one step, this              
c                    variable changes the load step size permenantly.           
c                                                                               
               perm_load_fact = perm_load_fact / two                            
               call errmsg ( 273, dum, dums, dumr, perm_load_fact )             
               goto 9999                                                        
            end if                                                              
         end if                                                                 
c                                                                               
      end do                                                                    
c                                                                               
c            move to the next entry in the list, unless we are at the end.      
c                                                                               
      prev_node_ptr = node_ptr                                                  
      node_ptr = crack_front_nodes( node_ptr, 2 )                               
      if ( node_ptr .ne. -1 ) go to 10                                          
c                                                                               
 9999 continue                                                                  
c                                                                               
      if ( debug ) then                                                         
         write (out,'("  ===== perm load fact:",e13.6)') perm_load_fact         
         write(out,*) '<<<<< leaving CTOA_cut_step'                             
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine CTOA_cut_step_const             *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 4/8/96                     *          
c     *                                                              *          
c     *         This routine checks if the load step size is too     *          
c     *         large. We find all master nodes with CTOAs           *          
c     *         larger than the critical angle.                      *          
c     *         If the number of load steps since release of         *          
c     *         of previous master node is smaller than the minimum  *          
c     *         number ( min_steps_for_release ), then cut the load  *          
c     *         step size by half.                                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine CTOA_cut_step_const ( debug )                                  
      use global_data ! old common.main
      use node_release_data, only : inv_crkpln_nodes,                           
     &     crkpln_nodes_state, master_lines                                     
      use damage_data                                                           
      implicit integer (a-z)                                                    
c                                                                               
      double precision                                                          
     &     d32460, angle, two, crit_angle, one, zero                            
      character(len=1) :: dums                                                  
      real dumr                                                                 
      data d32460, two, one, zero / 32460.0, 2.0, 1.0, 0.0 /                    
      logical debug, use_init                                                   
c                                                                               
c                                                                               
      if ( debug ) write (out,*) '>>>>>> in CTOA_cut_step'                      
c                                                                               
c          loop over the master lines for the master_nodes. If                  
c          node is about to be released, find how many steps                    
c          since nearest neighbor was released.  If less than                   
c          min_steps_for_release, cut step size in half.                        
c                                                                               
      do num_line = 1, num_crack_fronts                                         
c                                                                               
         if (debug) write (out,'("    check front:",i7)') num_line              
         call get_slope_master_line( num_line, use_init, angle, idummy )        
         if (angle .eq. zero) cycle                                             
         if ( debug ) write (out,'("     Angle is: ",e13.6)') angle             
c                                                                               
         if ( use_init ) then                                                   
            crit_angle = init_crit_ang * (one - CTOA_range)                     
         else                                                                   
            crit_angle = critical_angle * (one - CTOA_range)                    
         end if                                                                 
c                                                                               
c                if current angle is larger than critical angle, then           
c                node will be released this step. Check                         
c                if number of steps since release of neighbor is less           
c                the min_steps_for_release.                                     
c                                                                               
         if ( angle .gt. crit_angle )then                                       
            neighbor_state = crkpln_nodes_state( inv_crkpln_nodes               
     &                       (master_lines(num_line,2)) )                       
            if ( neighbor_state.gt.0 .and.                                      
     &           neighbor_state.lt.min_steps_for_release) then                  
c                                                                               
c                    too few steps between release. cut rate in                 
c                    half.  Note that we change the perm_load_fact              
c                    variable -- unlike the temp_load_fact, which               
c                    changes the load step size for one step, this              
c                    variable changes the load step size permenantly.           
c                                                                               
               perm_load_fact = perm_load_fact / two                            
               call errmsg ( 273, dum, dums, dumr, perm_load_fact )             
               goto 9999                                                        
            end if                                                              
         end if                                                                 
c                                                                               
      end do                                                                    
c                                                                               
 9999 continue                                                                  
c                                                                               
      if ( debug ) then                                                         
         write (out,'("  ===== perm load fact:",e13.6)') perm_load_fact         
         write(out,*) '<<<<< leaving CTOA_cut_step'                             
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
