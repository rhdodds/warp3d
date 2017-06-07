c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine find_neighbors               *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified : 08/22/95                   *          
c     *                                                              *          
c     *     this subroutine finds all of the neighboring crack       *          
c     *     plane nodes to a given crack plane node. This routine    *          
c     *     assumes that all of the elements are l3disop w/ 8 nodes  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine find_neighbors (node, num_neighbors, node_list,                
     &     crack_node, inv_crkpln_nodes, crack_plane_normal )                   
      use global_data ! old common.main
c                                                                               
      use main_data, only : incmap, incid, cnstrn_in,                           
     &                      inverse_incidences                                  
c                                                                               
c              passed parameters:                                               
c                                                                               
c                   node <int> -- input                                         
c                       node for which we are finding the neighbors             
c                   inv_crkpln_nodes (array) <int> -- input                     
c                       points from global node number to entry                 
c                       in crack_plane data structures                          
c                   crack_plane_normal <int> -- input                           
c                       direction of the normal to the crack plane              
c                                                                               
c                   num_neighbors <int> -- output                               
c                       number of neighboring nodes to 'node'                   
c                   node_list(mxconn)  <int> -- output                          
c                       list of neighboring nodes to 'node' on the              
c                       crack plane                                             
c                   crack_node <logical> -- output                              
c                       true if node is one the crack front                     
c                       (if it is constrained and has one                       
c                       or more neighboring nodes that are                      
c                       unconstrained)                                          
c                                                                               
c                   NOTE: we assume 8 nodes per element here!                   
c                                                                               
c                                                                               
      implicit integer (a-z)                                                    
      real dumr                                                                 
      character(len=1) :: dums                                                  
      double precision                                                          
     &     d32460, dumd                                                         
      data d32460 / 32460.0/                                                    
      dimension node_list(mxconn), adj_node_inc(3,8),                           
     &     inv_crkpln_nodes(*)                                                  
      logical crack_node, debug, already_in_list                                
      data adj_node_inc /2,4,5,1,3,6,2,4,7,1,3,8,1,6,8,2,5,7,3,6,8,             
     &     4,5,7/                                                               
      data debug /.false./                                                      
c                                                                               
      if ( debug ) write (*,*) '>>>> entering find_neighbors'                   
c                                                                               
      do i = 1, mxconn                                                          
         node_list(i) = 0                                                       
      enddo                                                                     
      num_neighbors = 0                                                         
      crack_node = .false.                                                      
c                                                                               
c               get the number of connected elements to the given node          
c               from the inverse incidences.  For each element,                 
c               use the incidence list relations to find the 3 nodes            
c               on the element that share a side with the given element.        
c               If the current element does not have 8 nodes,                   
c               skip the element.                                               
c                                                                               
c                                                                               
c                     NOTE: we assume 8 nodes per element here!                 
c                                                                               
      num_elems = inverse_incidences(node)%element_count                        
      do elem_list = 1, num_elems                                               
         elem = inverse_incidences(node)%element_list(elem_list)                
         if ( iprops(2,elem).ne. 8 ) cycle                                      
         incptr = incmap(elem) - 1                                              
c                                                                               
c                     The neighboring nodes on the element can be found         
c                     from the incidence list.                                  
c                                                                               
c                     find which element node the given node is                 
c                                                                               
         do i = 1, 8                                                            
            if (node .eq. incid(incptr+i)) then                                 
               node_inc = i                                                     
               exit                                                             
            endif                                                               
         enddo                                                                  
c                                                                               
c                     using the incidence relations, loop over the              
c                     3 nodes adjacent to the node tfor wich we are             
c                     locating the neighbors                                    
c                                                                               
         do adj_loop = 1, 3                                                     
            adj_node = incid(incptr+adj_node_inc(adj_loop,node_inc))            
c                                                                               
c                                                                               
c                         check if adj_node is on crack plane.                  
c                         if it is not on the crack plane, go to next           
c                         adjacent node.  If it is on the crack plane,          
c                         and it is the first neighbor found, start the list.   
c                                                                               
            if ( inv_crkpln_nodes(adj_node).eq.0 ) cycle                        
            if ( num_neighbors .eq. 0 ) then                                    
               num_neighbors = 1                                                
               node_list(1) = adj_node                                          
            else                                                                
c                                                                               
c                         the node is on the crack plane, but                   
c                         is not the first node found. check                    
c                         to see if it is already in the list.                  
c                                                                               
               already_in_list = .false.                                        
               do i = 1, num_neighbors                                          
                  if ( node_list(i) .eq. adj_node ) then                        
                     already_in_list = .true.                                   
                     exit                                                       
                  endif                                                         
               enddo                                                            
               if ( already_in_list ) cycle                                     
c                                                                               
c                                                                               
c                         node is not on list -- add to list.                   
c                                                                               
               num_neighbors = num_neighbors + 1                                
               if ( num_neighbors .gt.mxconn ) then                             
                  call errmsg(252,node,dums,dumr,dumd)                          
                  call die_abort                                                
                  stop                                                          
               endif                                                            
               node_list(num_neighbors) = adj_node                              
            endif                                                               
c                                                                               
c                         We have a new neighbor node.  Check if this           
c                         new node makes the check node a crack front           
c                         node. It does if:                                     
c                            a) the check node is constrained in the            
c                                  crack plane normal direction, and            
c                            b) new neighbor node is unconstrained in the       
c                                  crack plane normal direction,                
c                         then 'node' is on the crack front;                    
c                         set crack_node as true.                               
c                                                                               
            axis = crack_plane_normal-1                                         
            if ( (cnstrn_in(dstmap(node)+axis).ne.d32460) .and.                 
     &           cnstrn_in(dstmap(adj_node)+axis) .eq.d32460 )                  
     &           crack_node = .true.                                            
c                                                                               
         enddo                                                                  
      enddo                                                                     
c                                                                               
c                                                                               
      if ( debug ) then                                                         
         write (*,*) '>>>> For node ',node                                      
         write (*,*) '  # of neighbors:',num_neighbors                          
         write (*,*) '  list of neighbors:'                                     
         do i=1, num_neighbors                                                  
            write (*,'(6x,i5)') node_list(i)                                    
         enddo                                                                  
         write (*,'("   crack_node is ",l5)') crack_node                        
         write (*,*)                                                            
         write (*,*) '<<<<< leaving find_neighbors'                             
      endif                                                                     
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
