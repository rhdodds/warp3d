c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oupads                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 3/3/2016 rhd               *          
c     *                                                              *          
c     *     this subroutine adds the nodal results from a block of   *          
c     *     non-conflicting, similar elements into the total nodal   *          
c     *     results data structure.                                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oupads( span, belinc, nnode, num_vals, outmap,                 
     &                   num_struct_nodes, nodal_values )                       
      use elblk_data, only : elestr                                             
      implicit integer (a-z)                                                    
c                                                                               
      include 'param_def'                                                       
      dimension belinc(nnode,*), outmap(*)                                      
c                                                                               
      type :: node_entry                                                        
         integer :: count                                                       
         double precision, dimension(:), pointer :: node_values                 
      end type node_entry                                                       
c                                                                               
      type (node_entry), dimension (num_struct_nodes) :: nodal_values           
         double precision, dimension(:), pointer :: snode_values                
c                                                                               
c                        local declarations                                     
c                                                                               
      double precision                                                          
     &   zero                                                                   
      data zero / 0.0d0 /                                                       
c                                                                               
c                       look at nodes present on each element of                
c                       this block. if the vector of node average               
c                       values does not exist for a node, create and            
c                       zero it. update the counter for each structural node    
c                       accessed in the processing of the element               
c                       block.                                                  
c                                                                               
      do j = 1, nnode                                                           
       do i = 1, span                                                           
         snode = belinc(j,i)                                                    
         if ( .not. associated(nodal_values(snode)%node_values) ) then          
             allocate( nodal_values(snode)%node_values(mxstmp) )                
             nodal_values(snode)%node_values(1:mxstmp) = zero                   
         end if                                                                 
         nodal_values(snode)%count = nodal_values(snode)%count + 1              
       end do                                                                   
      end do                                                                    
c                                                                               
c                       add the results from each node of each                  
c                       element in this block into the global nodal results     
c                       data structure.                                         
c                                                                               
      do k = 1, num_vals                                                        
         map = outmap(k)                                                        
         do j = 1, nnode                                                        
          do i = 1, span                                                        
            snode = belinc(j,i)                                                 
            snode_values => nodal_values(snode)%node_values                     
            snode_values(map) = snode_values(map) + elestr(i,k,j)               
          end do                                                                
         end do                                                                 
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
