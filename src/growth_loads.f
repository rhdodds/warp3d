c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine growth_loads                 *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 01/03/96                   *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine growth_loads                                                   
      use global_data ! old common.main
c                                                                               
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
      double precision                                                          
     &   zero                                                                   
      logical debug                                                             
      data zero, debug / 0.0, .false./                                          
c                                                                               
c                                                                               
c            If crack growth is present in this analysis, slowly release        
c            any loads imposed on the structure due to killed elements          
c            or released nodes.  This release ensures smooth convergence        
c            of the solution while the geometry of the model is changed         
c            due to crack growth.                                               
c                                                                               
c                                                                               
      if ( debug ) write(*,*) '>>>> in growth_loads'                            
c                                                                               
      if ( growth_by_kill ) then                                                
         call killed_elem_loads( debug )                                        
      else if ( growth_by_release ) then                                        
         call released_node_loads( debug )                                      
      end if                                                                    
c                                                                               
c      if ( debug ) then                                                        
c         write (*,*) ' These are non-zero terms in new load_vector:'           
c         do i = 1, nodof                                                       
c            if ( load(i) .ne. zero ) write(*,'(" load(",i7,")=",e13.6)')       
c     &           i,load(i)                                                     
c         end do                                                                
c      end if                                                                   
c                                                                               
      if ( debug ) write(*,*) '<<<< leaving growth_loads'                       
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine killed_elem_loads            *          
c     *                                                              *          
c     *                       written by : asg                       *          
c     *                                                              *          
c     *                   last modified : 10/11/94                   *          
c     *                                                              *          
c     *     this subroutine releases the nodal loads that are        *          
c     *     imposed on the structure when an element is killed       *          
c     *     through element extinction (gurson) crack growth.        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine killed_elem_loads( debug )                                     
      use global_data ! old common.main
c                                                                               
      use elem_extinct_data, only : dam_ifv, dam_state, dam_dbar_elems          
      use main_data,         only : dload                                       
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      double precision                                                          
     &  dbar_now, refer_deform, now_fraction, last_fraction,                    
     &  fraction, dbar_zero, ldincr, zero, one, dumd                            
      real dumr                                                                 
      logical debug, no_print                                                   
      integer edest(mxedof)                                                     
      character(len=1) :: dums                                                  
      data zero, one / 0.0, 1.0 /                                               
c                                                                               
c          - For Gurson Crack Growth:                                           
c                                                                               
c            add a fraction of the internal load vector saved from              
c            killed elements to the dload vector to simulate the slow           
c            release of the killed element forces. there are two                
c            release options: (1) fixed no. of load steps, (2) a                
c            linear traction separation law.                                    
c                                                                               
c                                                                               
c            loop over all elements. skip non-killable elements                 
c            immediately. if element is killable, determine                     
c            if the internal forces present at extinction have                  
c            already been 100% released.                                        
c                                                                               
      if ( debug ) write(*,*) '>> dealing w/ killed element loads'              
      write (out,9200)                                                          
      no_print = .true.                                                         
      rcount   = 0                                                              
c                                                                               
      do elem = 1, noelem                                                       
        elem_ptr = dam_ptr(elem)                                                
        if ( elem_ptr .eq. 0 ) cycle                                            
        if ( dam_state(elem_ptr) .eq. 0 ) cycle                                 
        if ( release_type .eq. 1 ) then                                         
          if ( dam_state(elem_ptr) .gt. max_dam_state ) cycle                   
        end if                                                                  
        if ( release_type .eq. 2 ) then                                         
          if ( dam_dbar_elems(2,elem_ptr) .lt. zero ) cycle                     
        end if                                                                  
c                                                                               
c            element not yet 100% released. compute increment of                
c            internal forces at extinction to subtract off this                 
c            load step (notice we actually add due to force                     
c            sign convention).                                                  
c                                                                               
c            for fixed number of steps, the fraction is trivial.                
c                                                                               
        if ( release_type .eq. 1 ) then                                         
         fraction = one / dble(max_dam_state)                                   
             now_fraction =  dble(dam_state(elem_ptr)) /                        
     &                       dble(max_dam_state)                                
        end if                                                                  
c                                                                               
c            for traction-separation law we have to get current                 
c            elongation normal to crack plane and compute                       
c            relative change to start of release and the                        
c            relative to total fraction release up to now.                      
c            row 2 of dam_dbar_elems stores the released                        
c            fraction up to now.                                                
c                                                                               
        if ( release_type .eq. 2 ) then                                         
          call growth_set_dbar( elem, elem_ptr, debug, -2, dbar_now )           
          dbar_zero      = dam_dbar_elems(1,elem_ptr)                           
          refer_deform   = release_fraction * gurson_cell_size                  
          now_fraction   = (dbar_now-dbar_zero)/refer_deform                    
          last_fraction  = dam_dbar_elems(2,elem_ptr)                           
c                                                                               
          if ( now_fraction .gt. one ) then                                     
            now_fraction = one                                                  
            dam_dbar_elems(2,elem_ptr) = -one                                   
          else                                                                  
            dam_dbar_elems(2,elem_ptr) = now_fraction                           
          end if                                                                
          fraction = now_fraction - last_fraction                               
c                                                                               
        end if                                                                  
c                                                                               
c            output message about release of element forces                     
c                                                                               
      write (out,9300) elem, now_fraction * 100.0                               
      no_print = .false.                                                        
        rcount   = rcount + 1                                                   
c                                                                               
        num_dof = iprops(2,elem) * iprops(4,elem)                               
c                                                                               
c            modify incrmental load vector for structure for this               
c            step by the incremental force release for the                      
c            this element. also update the total load on structure.             
c                                                                               
        call get_single_edest_terms( edest, elem )                              
        do dof = 1, num_dof                                                     
          sdof        = edest(dof)                                              
          ldincr      = fraction * dam_ifv(dof,elem_ptr)                        
          load(sdof)  = load(sdof) + ldincr                                     
        end do                                                                  
c                                                                               
      end do                                                                    
c                                                                               
c         if no elements are undergoing release, output a message.              
c                                                                               
      if ( no_print ) write (out,9400)                                          
      if ( rcount .gt. 0 ) write(out,9500) rcount                               
      num_elements_in_force_release = rcount                                    
c                                                                               
      return                                                                    
 9000 format('  > current gurson release type 2 parameters: ',                  
     & /,10x,'element:     ',i7,                                                
     & /,10x,'elem_ptr:    ',i7,                                                
     & /,10x,'dbar_now,  dbar_zero: ',2f10.6,                                   
     & /,10x,'refer_deform, now_fraction, last_fraction, fraction',             
     & /,10x,5x,4f10.6 )                                                        
 9100 format('  > updated load for element:')                                   
 9150 format(10x,i3,i8,e14.6)                                                   
 9200 format(/1x,'  >> force release information for killed elements:')         
 9300 format(1x,'       element: ',i7,'. forces released (%): ',f5.1)           
 9400 format(1x,'        *no forces are currently being released*')             
 9500 format(1x,'       total elements in active release: ',i6)                 
                                                                                
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine released_node_loads          *          
c     *                                                              *          
c     *                       written by : asg                       *          
c     *                                                              *          
c     *                   last modified : 10/11/94                   *          
c     *                                                              *          
c     *     this subroutine releases the nodal loads that are        *          
c     *     imposed on the structure when a node is released         *          
c     *     through nodal release (discrete) crack growth.           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine released_node_loads( debug )                                   
      use global_data ! old common.main
      use node_release_data, only : crack_plane_nodes,                          
     &     crkpln_nodes_state, crkpln_nodes_react, node_release_frac,           
     &     crack_front_list, inv_crkpln_nodes                                   
      use damage_data                                                           
      implicit integer (a-z)                                                    
c                                                                               
      double precision                                                          
     &  now_fraction, fraction, new_height, ldincr, zero, one                   
      logical debug, no_print                                                   
      real dumr                                                                 
      character(len=1) :: dums                                                  
      data zero, one / 0.0, 1.0 /                                               
c                                                                               
c          - For Node Release Crack Growth:                                     
c                                                                               
c            add a fraction of the reaction forces saved from                   
c            released nodes to the dload vector to simulate the slow            
c            release of the reaction forces. there are two                      
c            release options: (1) fixed no. of load steps, (2) a                
c            linear traction separation law.                                    
c                                                                               
c            loop over all of the crack plane nodes.   If a node                
c            has not yet been released, then skip it.  Also, skip the           
c            node if it has been 100% released.                                 
c                                                                               
      if ( debug ) write (*,*) '>> dealing with released nodes'                 
      write (out, 9200)                                                         
      no_print = .true.                                                         
c                                                                               
c         1) if we are using the const_front algorithm and we are using the     
c            traction separation law, then we skip the normal procedure and     
c            release the other crack front nodes based on the release angle     
c            of the master node.                                                
c                                                                               
c           In this case, all of the crack front nodes that are                 
c           undergoing release are listed in the crack_front_list               
c           structure.  The first node in each front list in this               
c           structure is the master node for the front.  First                  
c           determine its release angle.  Then loop through the                 
c           other crack front nodes, and release them based on the              
c           master node's release fraction.                                     
c                                                                               
c                                                                               
      if ( const_front .and. release_type .eq. 2) then                          
         if (debug) write(*,*) 'Const_front and trac sep.'                      
c                                                                               
c               loop over each set of fronts in crack_front_list                
c                                                                               
         do i = 1, num_crack_fronts * num_nodes_grwinc                          
c                                                                               
c                   if first node is zero, skip this front                      
c                                                                               
            node = crack_front_list(i,1)                                        
            if (node .eq. 0) cycle                                              
          node_data_entry = inv_crkpln_nodes(node)                              
c                                                                               
c                   get release fraction for master node (see                   
c                below for explaination of traction-seperation                  
c                  algorithm)                                                   
c                                                                               
            new_height = u(dstmap(node)+crk_pln_normal_idx-1)                   
            now_fraction = crack_plane_sign*new_height/release_height           
            now_fraction = aint(now_fraction * 100.0) / 100.0                   
            if ( now_fraction.lt.zero ) now_fraction = - now_fraction           
            if ( now_fraction.gt.one ) now_fraction = one                       
c                                                                               
            fraction = now_fraction -                                           
     &           node_release_frac(node_data_entry)                             
c                                                                               
            if ( debug ) write(*,'("new_height,now_frac,frac:",                 
     &           3e16.9)') new_height, now_fraction, fraction                   
c                                                                               
            if ( fraction.lt.zero ) then                                        
               call errmsg( 253, node, dums, dumr, new_height )                 
               fraction = zero                                                  
            else                                                                
               node_release_frac(node_data_entry) = now_fraction                
            endif                                                               
c                                                                               
c                  write out information about the releasing force              
c                                                                               
          write (out,9300) node, now_fraction * 100.0                           
          no_print = .false.                                                    
c                                                                               
c                 for all nodes on the crack front, calculate the               
c                 amount of force to release based on the master                
c                  node release angle, and subtract that from the               
c                 total load vector.                                            
c                                                                               
            do j = 1, num_nodes_thick                                           
               node = crack_front_list(i,j)                                     
               node_data_entry = inv_crkpln_nodes(node)                         
             if (node .eq. 0) exit                                              
               dof = dstmap(node) + crk_pln_normal_idx - 1                      
               if ( debug )                                                     
     &            write(*,'(" change dof",i7," from load:",e13.6)')             
     &            dof, load(dof)                                                
               ldincr = fraction*crkpln_nodes_react(node_data_entry)            
               load(dof)  = load(dof) - ldincr                                  
               if ( debug )                                                     
     &            write (*,'(20x,"new load:",e13.6)') load(dof)                 
c                                                                               
c                 if the force has been fully released on the master node,      
c                 then zero out the terms in the crack front list so we can     
c                 reuse the space.                                              
c                                                                               
               if (now_fraction .eq. one) crack_front_list(i,j) = 0             
c                                                                               
            enddo                                                               
c                                                                               
         enddo                                                                  
         goto 9999                                                              
      endif                                                                     
c                                                                               
c           2)  handle cases of general growth or const. front growth           
c               with release steps                                              
c                                                                               
      do node_loop = 1, num_crack_plane_nodes                                   
         node = crack_plane_nodes(node_loop)                                    
         if ( crkpln_nodes_state(node_loop).eq.0 ) cycle                        
         if ( release_type .eq. 1 .and.                                         
     &        crkpln_nodes_state(node_loop).gt.max_dam_state ) cycle            
         if ( release_type .eq. 2 ) then                                        
           if ( node_release_frac(node_loop) .ge. one ) cycle                   
         end if                                                                 
c                                                                               
c            node not yet 100% released. compute increment of                   
c            reaction force at release to subtract off this                     
c            load step (notice we actually add due to force                     
c            sign convention).                                                  
c                                                                               
c            for fixed number of steps, the fraction is trivial. compute        
c            the total % released thus far for printing.                        
c                                                                               
        if ( release_type .eq. 1 ) then                                         
         fraction = one / dble(max_dam_state)                                   
         now_fraction = dble(crkpln_nodes_state(node_loop)) /                   
     &                      dble(max_dam_state)                                 
c                                                                               
c            if traction-separation law, get distance from crack plane.         
c            the height at which the force is completely released was           
c            calculated in dam_init_release2 as follows:                        
c                consider a line segment with length characteristic_length      
c                (which is specified in the user input) at an angle to the      
c                crack plane of release_fraction * critical_angle. One end      
c                point of the line rests on the crack plane.  The release       
c                height is the distance between the other endpoint and          
c                the crack plane.                                               
c            use the fraction distance_from_crack_plane/release_height          
c            to calculate the fraction of the reaction force to release.        
c                                                                               
c            There may be a problem with the node travelling in the             
c            reverse direction, moving through the crack plane. This            
c            causes the fraction to be negative. We simply change the           
c            sign of the fraction in this case, so we still release             
c            part of the force. round fraction to 2 significant                 
c            figures to remove effect of any small round-off in                 
c            displacements.                                                     
c                                                                               
c                                                                               
        else if ( release_type .eq. 2 ) then                                    
c                                                                               
           new_height = u(dstmap(node)+crk_pln_normal_idx-1)                    
           now_fraction = crack_plane_sign*new_height/release_height            
           now_fraction = aint(now_fraction * 100.0) / 100.0                    
           if ( now_fraction.lt.zero ) now_fraction = - now_fraction            
           if ( now_fraction.gt.one ) now_fraction = one                        
           fraction = now_fraction - node_release_frac(node_loop)               
c                                                                               
           if ( debug ) write(*,'("new_height,now_frac,frac:",3e16.9)')         
     &          new_height, now_fraction, fraction                              
c                                                                               
           if ( fraction.lt.zero ) then                                         
              call errmsg( 253, node, dums, dumr, new_height )                  
              fraction = zero                                                   
           else                                                                 
              node_release_frac(node_loop) = now_fraction                       
           endif                                                                
c                                                                               
        end if                                                                  
c                                                                               
c            write out information about the releasing force                    
c                                                                               
      write (out,9300) node, now_fraction * 100.0                               
      no_print = .false.                                                        
c                                                                               
c            modify incrmental load vector for structure for this               
c            step by the incremental force release for the                      
c            this node in the crack plane normal direction.                     
c                                                                               
        dof = dstmap(node) + crk_pln_normal_idx - 1                             
        if ( debug ) write(*,'(" change dof",i7," from load:",e13.6)')          
     &       dof, load(dof)                                                     
        ldincr = fraction*crkpln_nodes_react(node_loop)                         
        load(dof)  = load(dof) - ldincr                                         
        if ( debug ) write (*,'(20x,"new load:",e13.6)') load(dof)              
                                                                                
c                                                                               
      end do                                                                    
c                                                                               
 9999 continue                                                                  
c                                                                               
c          if no forces released, then output a message.                        
c                                                                               
      if ( no_print) write (out,9400)                                           
c                                                                               
c                                                                               
      return                                                                    
 9200 format(/1x,'  >> Force release information for released nodes:')          
 9300 format(1x,'     % of forces released from node ',i7,                      
     &        ' : ',f6.1)                                                       
 9400 format(1x,'     no forces are currently being released.')                 
      if ( debug ) write (*,*) '<< finished with released nodes'                
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
