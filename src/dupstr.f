c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dupstr                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 08/25/23 rhd               *          
c     *                                                              *          
c     *     creates a local(block) copy of element                   *          
c     *     data necessary for global stress vector recovery for     *          
c     *     each element in a block of similar elements.             *          
c     *     process only element data stored globally                *          
c     *     in non-blocked structures.                               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dupstr( span, bedst, belinc, felem,                            
     &                   nnode, ndof, totdof,                                   
     &                   trn_e_flags,                                           
     &                   trn_e_block,                                           
     &                   trne,                                                  
     &                   trnmte,                                                
     &                   ue,                                                    
     &                   due, trn,
     &                   killed_status  )                                             
c
      use global_data, only : out, mxvl, mxedof, sdispl=>u,
     &                        sdu => du   
      use elem_extinct_data, only : dam_state                                 
      use damage_data, only : dam_ptr, growth_by_kill, 
     &                        use_mesh_regularization
      use constants  
c                                                                               
      implicit none                                                    
c                                                                               
c           parameter declarations                                              
c                  
      integer, intent(in) :: span, felem, nnode, ndof, totdof                                                              
      integer, intent(in) :: bedst(totdof,*), belinc(nnode,*)                                
      logical :: trn_e_flags(*), trn_e_block, trne(mxvl,*), trn(*),
     &           killed_status(*)                                  
      double precision :: trnmte(mxvl,mxedof,*)
      double precision, intent(out) :: ue(mxvl,*), due(mxvl,*)                           
c                                                                               
c           local declarations                                                  
c                     
      integer :: k, i, j, element, elem_ptr                                                          
      logical :: standard_kill                                                      
      logical, parameter :: local_debug = .false.
c                                                                               
      if( local_debug ) write(out,9100)                                        
c                                                                               
c           gather element transformation flags.                                
c           then element transformation matrices.                               
c           Note: extensive comments on this are provided                       
c                 in dptstf.f                                                   
c                                                                               
      if( local_debug ) write(out,9200)                                       
c                                                                               
      trn_e_flags(1:span) = .false.                                                
c                                                                               
      trn_e_block         = .false.                                             
      k = 1                                                                     
      do j = 1, nnode                                                           
         do i = 1, span                                                         
            trne(i,j)      = trn(belinc(j,i))                                   
            trn_e_flags(i) = trn_e_flags(i) .or. trne(i,j)                      
            trn_e_block    = trn_e_block .or. trne(i,j)                         
         end do                                                                 
         k = k + 3                                                              
      end do                                                                    
c                                                                               
c           gather element transformation matrices                              
c                                                                               
      if( trn_e_block )  call duptrans( span, felem, trnmte )                  
c                                                                               
c            gather element displacements at state n and                        
c            the total displacement increment from n ->                         
c            n+1.                                      
c                                                                               
      do  j = 1, totdof                                                         
!DIR$ IVDEP
         do i = 1, span                                                         
            ue(i,j)  = sdispl(bedst(j,i))                                            
            due(i,j) = sdu(bedst(j,i))                                           
         end do                                                                 
      end do                                                                    
c          
c            set up for killable elements.
c                                                                     
c            zero total and incremental displacements                           
c            for (totally) killed elements. use staged checks since lower        
c            data exists only if we use crack growth by killing elements             
c            and elements have been killed. if dam_ptr for                      
c            first element in block = 0, the block has no                       
c            killable elements (all elements in a block must be                 
c            killable or non-killable).                                         
c                                                                                                                                                              
      if( .not. growth_by_kill )  go to 200                                     
      if( dam_ptr(felem) .eq. 0 ) go to 200  ! no killable elements in blk
      standard_kill = .not. use_mesh_regularization
c
c            setup so computed strains, stress, internal forces = 0.
c            the internal forces were grabbed and saved prior to this
c            at the moment the code killed the element.
c
      if( standard_kill ) then
        do i = 1, span                                                            
          element = felem + i - 1
          elem_ptr = dam_ptr(element)                                                        
          if( dam_state(elem_ptr) == 0 ) cycle ! not killable or not yet killed
!DIR$ IVDEP
          do j = 1, totdof                                                     
            ue(i,j)  = zero                                                    
            due(i,j) = zero                                                    
          end do                                                               
        end do  
        go to 200
      end if
c             
c            compute strains, stresses, internal forces assuming
c            no damage.  internal forces will be reduced later for 
c            damage. if damage exceeds limit (element fully removed),
c            setup so zero strains, stresses will be computed.
c            computed internal forces will be zeroed for fully 
c            damaged element
c                                                    
      if( use_mesh_regularization ) then
         do i = 1, span                                                            
            element = felem + i - 1   
            elem_ptr = dam_ptr(element) 
            if( dam_state(elem_ptr) == 0 ) cycle ! not killable                                           
            if( killed_status(i) ) then ! been fully removed
!DIR$ IVDEP
              do j = 1, totdof                                                     
                ue(i,j)  = zero                                                    
                due(i,j) = zero  
              end do                                                  
            end if   
         end do
         go to 200
      end if    
      write(out,9000) ! sanity check
      call die_gracefully
      return 
c                                                                               
c           if required, rotate total and incremental                           
c           displacements at specific nodes from constraint                     
c           compatible coordinates to global coordinates.                       
c                                                                               
 200  continue                                                                  
      if( trn_e_block ) then                                                   
        if( local_debug )  write(out,9500)                                     
        do i = 1, span                                                          
          if( trn_e_flags(i) ) then                                            
           call trnvec( ue, trnmte, trne, ndof, nnode, i, 2 )                   
           call trnvec( due, trnmte, trne, ndof, nnode, i, 2 )                  
          end if                                                                
        end do                                                                  
      end if                                                                    
c                                                                               
      if( local_debug ) write(out,9150)                                        
c                                                                               
      return                                                                    
c  
c
 9000 format(">>> FATAL ERROR: inconsistent state in dupstr"
     & /,".... job terminated ....",//)
 9100 format(8x,'>> entered dupstr...' )                                        
 9150 format(8x,'>> leaving dupstr...' )                                        
 9200 format(12x,'>> gather element transformation flags, matrices...' )        
 9400 format(12x,'>> gather element nodal displacement increments...' )         
 9500 format(12x,'>> transform element displ. to global coord. sys...' ) 
 9600 format(1x,"@0.1 routine dupstr. fully removed elem: ",i6)       
c                                                                               
      end                                                                       
