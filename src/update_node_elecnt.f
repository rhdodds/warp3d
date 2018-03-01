c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine update_node_elecnt           *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 06/30/94                   *          
c     *                                                              *          
c     *        Updates the vector that holds the number of           *          
c     *        elements connected to each node after a node is       *          
c     *        killed.                                               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine update_node_elecnt( ielem, debug )                             
      use global_data ! old common.main
c                                                                               
      use elem_extinct_data, only : dam_node_elecnt                             
      use main_data, only : incid, incmap, inverse_incidences                   
c                                                                               
      implicit integer (a-z)                                                    
      logical debug                                                             
c                                                                               
      if ( debug ) write (*,*) '>>>> in update_node_elecnt '                    
c                                                                               
c                loop over nodes connected to killed element.  Get node         
c                number.  Then decrease the count of the elements connected     
c                to that node by one.                                           
c                                                                               
      nnode = iprops(2,ielem)                                                   
      do i = 1, nnode                                                           
         node = incid(incmap(ielem)+(i-1))                                      
         elem = inverse_incidences(node)%element_list(1)                        
         ndof = iprops(4,elem)                                                  
         dam_node_elecnt(node)= dam_node_elecnt(node) - 1                       
      end do                                                                    
      if (debug) write (*,*) '>>>> leaving update_node_elecnt'                  
c                                                                               
      return                                                                    
 9000 format ('>> node ',i2,' of elem ',i7,' is ',i6)                           
 9010 format ('>>   elems connected to node is ',i2)                            
 9020 format('old constraints for node ',i7,' :',3(1x,e14.6))                   
      end                                                                       
