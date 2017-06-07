c     ****************************************************************          
c     *                                                              *          
c     *           extract indexes into the global                    *          
c     *           displacmenent, velocity, acceleration vectors      *          
c     *           for equations for a list of elements. this         *          
c     *           essentially de-blocks the indexes for a set of     *          
c     *           elements                                           *          
c     *                                                              *          
c     *                       written by  : rhd                      *          
c     *                   last modified : 02/18/98                   *          
c     *                                                              *          
c     ****************************************************************          
                                                                                
                                                                                
      subroutine get_edest_terms( table, elem_list, list_length )               
      use global_data ! old common.main
      use elem_block_data, only :  edest_blocks                                 
      use main_data,       only :  elems_to_blocks                              
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      dimension table(mxedof,*), elem_list(*)                                   
      integer, dimension (:,:), pointer :: edest                                
c                                                                               
      do i = 1, list_length                                                     
       elem = elem_list(i)                                                      
       if ( elem .le. 0 ) cycle                                                 
       totdof   = iprops(2,elem) * iprops(4,elem)                               
       blk      = elems_to_blocks(elem,1)                                       
       rel_elem = elems_to_blocks(elem,2)                                       
       edest    => edest_blocks(blk)%ptr                                        
       do dof = 1, totdof                                                       
         table(dof,i) = edest(dof,rel_elem)                                     
       end do                                                                   
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *           extract indexes into the global                    *          
c     *           displacmenent, velocity, acceleration vectors      *          
c     *           for equations for a single element.                *          
c     *                                                              *          
c     *                       written by  : rhd                      *          
c     *                   last modified : 02/18/98                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine get_single_edest_terms( edest_vec, elem )                      
      use global_data ! old common.main
      use elem_block_data, only :  edest_blocks                                 
      use main_data,       only :  elems_to_blocks                              
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      dimension edest_vec(*)                                                    
      integer, dimension (:,:), pointer :: edest                                
c                                                                               
      if ( elem .le. 0 ) return                                                 
      totdof   = iprops(2,elem) * iprops(4,elem)                                
      blk      = elems_to_blocks(elem,1)                                        
      rel_elem = elems_to_blocks(elem,2)                                        
      edest    => edest_blocks(blk)%ptr                                         
      do dof = 1, totdof                                                        
         edest_vec(dof) = edest(dof,rel_elem)                                   
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
