c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dupstr                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 01/20/12 rhd               *          
c     *                                                              *          
c     *     this subroutine creates a separate copy of  element      *          
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
     &                   due, trn )                                             
      use global_data ! old common.main
c                                                                               
      use elem_extinct_data, only : dam_state                                   
      use damage_data, only : dam_ptr, growth_by_kill                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c           parameter declarations                                              
c                                                                               
      dimension bedst(totdof,*), belinc(nnode,*)                                
      logical trn_e_flags(*),                                                   
     &        trn_e_block, trne(mxvl,*), trn(*)                                 
      double precision                                                          
     & trnmte(mxvl,mxedof,*), ue(mxvl,*), due(mxvl,*)                           
c                                                                               
c           local declarations                                                  
c                                                                               
      logical local_debug                                                       
      double precision                                                          
     &   half, zero                                                             
      data local_debug, half, zero / .false., 0.5, 0.0 /                        
c                                                                               
c                                                                               
      if ( local_debug ) write(out,9100)                                        
c                                                                               
c           gather element transformation flags.                                
c           then element transformation matrices.                               
c           Note: extensive comments on this are provided                       
c                 in dptstf.f                                                   
c                                                                               
      if ( local_debug )  write(out,9200)                                       
c                                                                               
      do i = 1, span                                                            
        trn_e_flags(i) = .false.                                                
      end do                                                                    
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
      if ( trn_e_block )  call duptrans( span, felem, trnmte )                  
c                                                                               
c            gather element displacements at state n and                        
c            the total displacement increment from n ->                         
c            n+1. zero displacements for killed elements. if iter=0             
c            and displacement increments are all zero, just leave.              
c            there is nothing to do.                                            
c                                                                               
      do  j = 1, totdof                                                         
!DIR$ VECTOR ALIGNED                                                            
         do i = 1, span                                                         
            ue(i,j)  = u(bedst(j,i))                                            
            due(i,j) = du(bedst(j,i))                                           
         end do                                                                 
      end do                                                                    
c                                                                               
c            zero total and incremental displacements                           
c            for killed elements. use staged checks since lower data            
c            exists only if we use crack growth by killing elements             
c            and elements have been killed. if dam_ptr for                      
c            first element in block = 0, the block has no                       
c            killable elements (all elements in a block must be                 
c            killable or non-killable).                                         
c                                                                               
c                                                                               
      if( .not. growth_by_kill )  go to 200                                     
      if( dam_ptr(felem) .eq. 0 ) go to 200                                     
      do i = 1, span                                                            
        element = felem + i - 1                                                 
        if ( dam_ptr(element) .eq. 0 ) cycle                                    
        if ( dam_state(dam_ptr(element)) .ne. 0 ) then                          
           do j = 1, totdof                                                     
             ue(i,j)  = zero                                                    
             due(i,j) = zero                                                    
           end do                                                               
        end if                                                                  
      end do                                                                    
c                                                                               
c           if required, rotate total and incremental                           
c           displacements at specific nodes from constraint                     
c           compatible coordinates to global coordinates.                       
c                                                                               
 200  continue                                                                  
      if ( trn_e_block ) then                                                   
        if ( local_debug )  write(out,9500)                                     
        do i = 1, span                                                          
          if ( trn_e_flags(i) ) then                                            
           call trnvec( ue, trnmte, trne, ndof, nnode, i, 2 )                   
           call trnvec( due, trnmte, trne, ndof, nnode, i, 2 )                  
          end if                                                                
        end do                                                                  
      end if                                                                    
c                                                                               
      if ( local_debug ) write(out,9150)                                        
c                                                                               
      return                                                                    
c                                                                               
 9100 format(8x,'>> entered dupstr...' )                                        
 9150 format(8x,'>> leaving dupstr...' )                                        
 9200 format(12x,'>> gather element transformation flags, matrices...' )        
 9400 format(12x,'>> gather element nodal displacement increments...' )         
 9500 format(12x,'>> transform element displ. to global coord. sys...' )        
c                                                                               
      end                                                                       
