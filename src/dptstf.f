c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dptstf                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 07/29/2017 rhd             *          
c     *                                                              *          
c     *     this subroutine creates a separate copy of all element   *          
c     *     data necessary for the tangent stiffness computation of  *          
c     *     each element in a block of similar, non-conflicting      *          
c     *     elements.                                                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dptstf( span, bedst, bcdst, belinc, felem, ngp,                
     &                   nnode, ndof, geonl, totdof, mat_type,                  
     &                   trn_e_flags,                                           
     &                   trn_e_block,                                           
     &                   ce,  ce_0,                                                  
     &                   trne,                                                  
     &                   trnmte,                                                
     &                   ue,                                                    
     &                   due,                                                   
     &                   local_cp,                                              
     &                   local_icp, trn, elem_type, surface,                    
     &                   cohesive_elem )                                        
      use global_data ! old common.main
      use elem_extinct_data, only : dam_state                                   
      use damage_data, only : dam_ptr, growth_by_kill                           
                                                                                
      implicit integer (a-z)                                                    
      dimension bedst(totdof,*), bcdst(totdof,*), belinc(nnode,*),              
     &          local_cp(*),                                                    
     &          local_icp(mxutsz,*)                                             
      logical   trn_e_flags(*), trn_e_block, trne(mxvl,*),                      
     &          geonl, trn(*), cohesive_elem,                                   
     &          middle_surface                                                  
      double precision                                                          
     & ce(mxvl,*), trnmte(mxvl,mxedof,*), ue(mxvl,*), due(mxvl,*),              
     & ce_orig(mxvl,mxecor), djcoh(mxvl), zero, ce_0(mxvl,mxecor)                                  
      data zero / 0.0 /                                                         
c                                                                               
c           transformation matrices are used to define a "local"                
c           coordinate systen for imposition of constraints, e.g.,              
c           skewed. if transformations are present, they are used               
c           to rotate the computed element stiffness terms for those            
c           nodes with constraint coordinate systems. also, for                 
c           geometric nonlinear elements, the element coordinates must          
c           be updated and so the local nodal displacements                     
c           must be rotated back to global -- all element computations          
c           are performed in global coordinates.                                
c                                                                               
c           build table of element nodal coordinates and                        
c           transformation flags for block.                                     
c           this is done for each node of each element in the block.            
c           we then set a unique flag for each element                          
c           and for the whole block (which can be skipped most times).          
c                                                                               
      do i = 1, span                                                            
        trn_e_flags(i) = .false.                                                
      end do                                                                    
c                                                                               
      trn_e_block         = .false.                                             
      k = 1                                                                     
      do j = 1, nnode                                                           
!DIR$ VECTOR ALIGNED                                                            
         do i = 1, span                                                         
            ce(i,k)        = c(bcdst(k,i))                                      
            ce(i,k+1)      = c(bcdst(k+1,i))                                    
            ce(i,k+2)      = c(bcdst(k+2,i))                                    
            ce_0(i,k)      = c(bcdst(k,i))                                      
            ce_0(i,k+1)    = c(bcdst(k+1,i))                                    
            ce_0(i,k+2)    = c(bcdst(k+2,i))                                    
        end do                                                                 
         do i = 1, span                                                         
            trne(i,j)      = trn(belinc(j,i))                                   
            trn_e_flags(i) = trn_e_flags(i) .or. trne(i,j)                      
            trn_e_block    = trn_e_block .or. trne(i,j)                         
         end do                                                                 
         k = k + 3                                                              
      end do                                                                    
c                                                                               
c           gather element transformation matrices -                            
c           all dof of all elements in block. skip dof that have                
c           no transformation  - their row in transformation matrix             
c           table is 0.                                                         
c                                                                               
      if ( trn_e_block ) call duptrans(  span, felem, trnmte )                  
c                                                                               
c           gather data for geometrically nonlinear formulation.                
c                                                                               
c             a) get element nodal displacements at start of                    
c                step (ue) and increment over step (due)                        
c             b) rotate from constraint (local nodal) to global                 
c                for elements with nodes that have the local                    
c                system. can possibly skip whole block.                         
c             c) update global coordinates of all element nodes                 
c                for this element block to current estimate for                 
c                end of step configuration.                                     
c             d) for cohesive elements, compute the reference surface           
c                for computations                                               
c             e) for killed elements in the block, set coordinates              
c                back to initial (t=0) values.                                  
c             f) gather unrotated cauchy stresses at n+1                        
c             g) gather [R,n+1] to transform unrotated cauchy                   
c                stresses to cauchy stresses at n+1                             
c                                                                               
c                                                                               
      if ( .not. geonl ) go to 300                                              
c                                                                               
      do  j = 1, totdof                                                         
!DIR$ VECTOR ALIGNED                                                            
       do i = 1, span                                                           
          ue(i,j)  = u(bedst(j,i))                                              
          due(i,j) = du(bedst(j,i))                                             
       end do                                                                   
      end do                                                                    
c                                                                               
      if ( trn_e_block ) then                                                   
       do i = 1, span                                                           
         if ( trn_e_flags(i) ) then                                             
           call trnvec( ue, trnmte, trne, ndof, nnode, i, 2 )                   
           call trnvec( due, trnmte, trne, ndof, nnode, i, 2 )                  
         end if                                                                 
        end do                                                                  
      end if                                                                    
c                                                                               
      do j = 1, totdof                                                          
!DIR$ VECTOR ALIGNED                                                            
       do i = 1, span                                                           
         ce_orig(i,j) = ce(i,j)                                                 
         ce(i,j)      = ce(i,j) + ue(i,j) + due(i,j)                            
       end do                                                                   
      end do                                                                    
c                                                                               
c           set up for the cohesive elements.                                   
c                                                                               
      middle_surface = surface .eq. 2                                           
!DIR$ VECTOR ALIGNED                                                            
      djcoh(1:span)  = zero                                                     
      if( cohesive_elem ) then                                                  
        if( middle_surface ) call chk_cohes_penetrate( span, mxvl,              
     &               felem, mxndel, nnode, elem_type, ce, djcoh )               
       call cohes_ref_surface( span, mxvl, mxecor, surface, nnode,              
     &                          totdof, ce_orig, ce, djcoh )                    
      end if                                                                    
c                                                                               
c           all elements in a block must be killable or                         
c           not killable. skip non-killable blocks here.                        
c           use staged checks since lower data                                  
c           exists only if we use crack growth by killing elements              
c           and elements have been killed. if dam_ptr for                       
c           first element in block = 0, the block has no                        
c           killable elements (all elements in a block must be                  
c           killable or non-killable).                                          
c                                                                               
      if( .not. growth_by_kill ) go to 300                                      
      if(  dam_ptr(felem) .eq. 0 ) go to 300                                    
      do i = 1, span                                                            
        element = felem + i - 1                                                 
        if ( dam_ptr(element) .eq. 0 ) cycle                                    
        if ( dam_state(dam_ptr(element)) .ne. 0 ) then                          
          do j = 1, 3*nnode                                                     
             ce(i,j) = c(bcdst(j,i))                                            
          end do                                                                
        end if                                                                  
      end do                                                                    
c                                                                               
c                 make copies of vectors of subscripts                          
c                 used in matrix multiplies for element                         
c                 stiffnesses. this reduces access to shared                    
c                 global variables during parallel processing                   
c                 of element blocks.                                            
c                                                                               
 300  continue                                                                  
      do i = 1, mxedof                                                          
       local_cp(i) = cp(i)                                                      
      end do                                                                    
c                                                                               
      do i = 1, mxutsz                                                          
       local_icp(i,1) = icp(i,1)                                                
       local_icp(i,2) = icp(i,2)                                                
      end do                                                                    
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
                                                                                
                                                                                
