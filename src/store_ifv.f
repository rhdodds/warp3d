c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine store_ifv                    *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 06/06/94                   *          
c     *                   last modified : 07/29/95 (rhd)             *          
c     *                                                              *          
c     *        continue to reduce forces a killed element imposes    *          
c     *        on the structure                                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine store_ifv( elem, elem_ptr, debug )                             
      use global_data ! old common.main
c                                                                               
      use elem_extinct_data, only : dam_ifv, dam_state                          
      use elem_block_data,   only : edest_blocks                                
      use main_data,         only : elems_to_blocks                             
      use damage_data, only : release_type                                      
c                                                                               
      implicit integer (a-z)                                                    
      double precision                                                          
     &     zero, one                                                            
      logical debug                                                             
      integer, dimension(:,:), pointer :: edest                                 
      data zero, one /0.0,1.0/                                                  
c                                                                               
      if ( debug ) write (*,*) '>>>>> Entering store_ifv'                       
c                                                                               
c            initialize state of release for element                            
c                                                                               
      dam_state(elem_ptr) = 1                                                   
c                                                                               
c            subtract the element internal forces from the                      
c            structure load vector. for all killable elements                   
c            which are not in the process of releasing, we                      
c            grab their internal forces at start of each load                   
c            step in another routine so they are available                      
c            here.                                                              
c                                                                               
      num_edof = iprops(2,elem) * iprops(4,elem)                                
c                                                                               
c            for the traction separation release, provide                       
c            an option to immediately                                           
c            drop the internal forces by a fraction, then release the           
c            remaining forces linearly with opening displacement.               
c            current reduction is zero.                                         
c                                                                               
      if ( release_type .eq. 2 ) then                                           
        do i = 1, num_edof                                                      
            dam_ifv(i,elem_ptr) = one * dam_ifv(i,elem_ptr)                     
        end do                                                                  
      end if                                                                    
c                                                                               
      blk         = elems_to_blocks(elem,1)                                     
      rel_elem    = elems_to_blocks(elem,2)                                     
      edest       => edest_blocks(blk)%ptr                                      
      do i = 1, num_edof                                                        
        load(edest(i,rel_elem)) = load(edest(i,rel_elem)) -                     
     &                            dam_ifv(i,elem_ptr)                           
      end do                                                                    
c                                                                               
      if ( debug ) then                                                         
         write (*,*) '>>>>> load in storeifv for elem:',elem                    
         do i = 1, num_edof                                                     
            write (*,'("new load:",i5,1x,e14.6)') edest(i,rel_elem),            
     &        load(edest(i,rel_elem))                                           
         end do                                                                 
         write (*,*) '<<<<< Leaving store_ifv'                                  
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
