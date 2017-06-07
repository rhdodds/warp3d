c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dupmas                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 03/13/04 rhd               *          
c     *                                                              *          
c     *     this subroutine creates a separate copy of all element   *          
c     *     data necessary for the mass computations of each         *          
c     *     element in a block of similar, non-conflicting elements. *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dupmas( span, nnode, bcdst, c, totdof, ce_block,               
     &                   mxvl )                                                 
      implicit integer (a-z)                                                    
      double precision                                                          
     &  c(*), ce_block(mxvl,*)                                                  
      integer bcdst(totdof,*)                                                   
c                                                                               
c                       gather element coordinates.                             
c                                                                               
      do j = 1, 3*nnode                                                         
!DIR$ VECTOR ALIGNED                                                            
         do i = 1, span                                                         
            ce_block(i,j) = c(bcdst(j,i))                                       
         end do                                                                 
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dupmas_fgm                   *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 3/13/04 rhd                *          
c     *                                                              *          
c     *     this subroutine creates a separate copy of all element   *          
c     *     data necessary for the mass computations of each         *          
c     *     element in a block of similar, non-conflicting elements. *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dupmas_fgm( span, nnode, belinc, fgm_node_values,              
     &                       num_struct_nodes, rho_block, mxndel )              
      implicit integer (a-z)                                                    
c                                                                               
c                        parameter declarations                                 
c                                                                               
      real  fgm_node_values(num_struct_nodes,*)                                 
      double precision                                                          
     &   rho_block(mxndel,*)                                                    
      integer belinc(nnode,*)                                                   
c                                                                               
c                 if the model has fgm properties at the model nodes,           
c                 build a table of values for nodes of elements in the          
c                 block                                                         
c                                                                               
       do i = 1, nnode                                                          
!DIR$ VECTOR ALIGNED                                                            
         do k = 1, span                                                         
            rho_block(i,k) = fgm_node_values(belinc(i,k),5)                     
         end do                                                                 
       end do                                                                   
c                                                                               
      return                                                                    
      end                                                                       
