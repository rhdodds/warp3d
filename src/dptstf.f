c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dptstf                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 03/17/21 rhd               *          
c     *                                                              *          
c     *     creates a separate (local) copy of all element data      *
c     *     necessary for the tangent stiffness computation of       *          
c     *     each element in the block                                *          
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
c        
      use global_data, only : out, mxvl, mxecor, mxedof, mxutsz,
     &                        mxndel, scoords => c, icp, sdispl=>u,
     &                        sdu => du, cp   
      use elem_extinct_data, only : dam_state, smcs_d_values                                   
      use damage_data, only : dam_ptr, growth_by_kill, 
     &                        use_mesh_regularization, tol_regular =>
     &                        tolerance_mesh_regularization  
      use constants  
c
      implicit none
c
c                      parameter declarations
c
      integer, intent(in) :: span, felem, ngp, nnode, ndof, totdof, 
     &                       mat_type, elem_type, surface
      logical :: trn_e_flags(*), trn_e_block, trne(mxvl,*),                     
     &           geonl, trn(*), cohesive_elem,                                   
     &           middle_surface      
      integer, intent(in) :: bedst(totdof,*), bcdst(totdof,*),
     &         belinc(nnode,*)                
      integer, intent(out) :: local_cp(*), local_icp(mxutsz,*)                                             
      double precision :: ce(mxvl,*), trnmte(mxvl,mxedof,*), 
     &                    ue(mxvl,*), due(mxvl,*), 
     &                    ce_orig(mxvl,mxecor), djcoh(mxvl),
     &                    ce_0(mxvl,mxecor)     
c
c                       locally defined variables
c
      integer :: i, j, k, element, elem_ptr
      logical :: standard_kill
c                                                                               
c                 make copies of vectors of subscripts                          
c                 used in matrix multiplies for element                         
c                 stiffnesses. this reduces access to shared                    
c                 global variables during parallel processing                   
c                 of element blocks.                                            
c                                                                               
       local_cp(1:mxedof) = cp(1:mxedof)                                                      
       local_icp(1:mxutsz,1:2) = icp(1:mxutsz,1:2)                                                
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
      trn_e_flags(1:span) = .false.                                                
c                                                                               
      trn_e_block         = .false.                                             
      k = 1                                                                     
      do j = 1, nnode                                                           
         do i = 1, span                                                         
            ce(i,k)        = scoords(bcdst(k,i))                                      
            ce(i,k+1)      = scoords(bcdst(k+1,i))                                    
            ce(i,k+2)      = scoords(bcdst(k+2,i))                                    
            ce_0(i,k)      = scoords(bcdst(k,i))                                      
            ce_0(i,k+1)    = scoords(bcdst(k+1,i))                                    
            ce_0(i,k+2)    = scoords(bcdst(k+2,i))                                    
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
      if( trn_e_block ) call duptrans( span, felem, trnmte )                  
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
c             e) for (totally) killed elements in the block, set coordinates              
c                back to initial (t=0) values.                                  
c             f) gather unrotated cauchy stresses at n+1                        
c             g) gather [R,n+1] to transform unrotated cauchy                   
c                stresses to cauchy stresses at n+1                             
c                                                                               
c                                                                               
      if( .not. geonl ) return
c                                                                               
      do  j = 1, totdof                                                         
       do i = 1, span                                                           
          ue(i,j)  = sdispl(bedst(j,i))                                              
          due(i,j) = sdu(bedst(j,i))                                             
       end do                                                                   
      end do                                                                    
c                                                                               
      if( trn_e_block ) then                                                   
       do i = 1, span                                                           
         if( trn_e_flags(i) ) then                                             
           call trnvec( ue, trnmte, trne, ndof, nnode, i, 2 )                   
           call trnvec( due, trnmte, trne, ndof, nnode, i, 2 )                  
         end if                                                                 
       end do                                                                  
      end if                                                                    
c                                                                               
      do j = 1, totdof                                                          
       do i = 1, span                                                           
         ce_orig(i,j) = ce(i,j)                                                 
         ce(i,j)      = ce(i,j) + ue(i,j) + due(i,j)                            
       end do                                                                   
      end do                                                                    
c                                                                               
c           set up for the cohesive elements.                                   
c                                                                               
      middle_surface = surface .eq. 2                                           
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
      if( .not. growth_by_kill ) return
      if( dam_ptr(felem) .eq. 0 ) return ! no killable elements this block
      standard_kill = .not. use_mesh_regularization
c
      if( standard_kill ) then
         do i = 1, span                                                            
            element = felem + i - 1   
            elem_ptr = dam_ptr(element)                                              
            if( dam_state(elem_ptr) == 0 ) cycle ! not yet killed
            do j = 1, 3*nnode   ! killed                                                
               ce(i,j) = scoords(bcdst(j,i))                                            
            end do                                                                
         end do 
         return
      end if
c
      if( use_mesh_regularization ) then
         do i = 1, span                                                            
            element = felem + i - 1   
            elem_ptr = dam_ptr(element)                                              
            if( smcs_d_values(elem_ptr) > tol_regular ) then ! fully killed
               do j = 1, 3*nnode                                                     
                 ce(i,j) = scoords(bcdst(j,i))                                            
               end do                                                                
            end if                                                                  
         end do 
         return
      end if
c
      write(out,9000)
      call die_abort                                                                   
      return 
c
 9000 format(">>> FATAL ERROR: inconsistent state in dptstf"
     & /,".... job terminated ....",//)
c                                                                               
      end                                                                       
                                                                                
                                                                                
