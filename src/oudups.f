c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oudups                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 4/26/2017 rhd              *          
c     *                                                              *          
c     *     gathers data for a block of elements to support          *          
c     *     generation of a patran, packet or hard copy output.      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oudups( span, felem, ngp, geonl, stress, is_cohesive )         
      use global_data ! old common.main
      use main_data, only : elems_to_blocks                                     
      use elem_block_data, only : history_blocks, rot_n1_blocks,                
     &                            eps_n_blocks, urcs_n_blocks,                  
     &                            history_blk_list                              
      use elblk_data, only : elem_hist, blk_size_hist, urcs_blk_n,              
     &                       rot_blk_n1, ddtse, blk_size_gp                     
c                                                                               
      implicit none                                                             
c                                                                               
      integer :: span, felem, ngp                                               
      logical :: geonl, stress, is_cohesive                                     
c                                                                               
c             local declarations                                                
c                                                                               
      integer :: blk, rel_elem, hist_size, hist_offset, rot_offset,             
     &           eps_offset, sig_offset, mxhist, mxngp                          
      logical :: is_solid                                                       
      logical, parameter ::  local_debug = .false.                              
c                                                                               
      blk         = elems_to_blocks(felem,1)                                    
      rel_elem    = elems_to_blocks(felem,2)                                    
      hist_size   = history_blk_list(blk)                                       
      hist_offset = (rel_elem-1)*hist_size*ngp + 1                              
      rot_offset  = (rel_elem-1)*9*ngp + 1                                      
      eps_offset  = (rel_elem-1)*nstr*ngp + 1                                   
      sig_offset  = (rel_elem-1)*nstrs*ngp + 1                                  
      is_solid    = .not. is_cohesive                                           
c                                                                               
      if( local_debug ) then                                                    
         write(out,*) '..... oudups....'                                        
         write(out,*) '....    span, felem, ngp, geonl, stress'                 
         write(out,*) span, felem, ngp, geonl, stress                           
         write(out,*) '.... blk, rel_elm, eps_offset: '                         
         write(out,*) blk, rel_elem, eps_offset                                 
         write(out,*) '.... is_cohesive: ', is_cohesive                         
      end if                                                                    
c                                                                               
c             gather history data. careful: the local                           
c             block size may be larger than stored block size.                  
c             uses non-standard gastr routine !                                 
c                                                                               
c             history data:                                                     
c              o The global blocks are sized(hist_size,ngp,span)                
c              o The local block is sized (mxvl,mxhist,mxngp).                  
c                -> mxhist: for all element blocks, the maximum                 
c                           no. of words of history data per                    
c                           gauss point                                         
c                -> mxngp:  for all elements blocks, the maximum                
c                           no. of integration points for an element            
c                                                                               
c              This makes it possible to pass a 2-D array slice for             
c              all elements of the block for a single gauss point.              
c                                                                               
      mxhist = blk_size_hist                                                    
      mxngp  = blk_size_gp                                                      
      call ou_gastr( elem_hist(1,1,1),                                          
     &               history_blocks(blk)%ptr(hist_offset),                      
     &               ngp, mxhist, mxngp, hist_size, span, mxvl )                
c                                                                               
c             gather stresses. for geonl, gather [Rot], the current             
c             rotation for transforming unrotated cauchy stresses               
c             to cauchy stresses (skip interface-cohesive elements).            
c             if not stresses, gather strain data.                              
c                                                                               
      if( geonl .and. is_solid )                                                
     &   call tanstf_gastr( rot_blk_n1,                                         
     &      rot_n1_blocks(blk)%ptr(rot_offset), ngp, 9, span )                  
c                                                                               
      if( stress ) then                                                         
        call tanstf_gastr( urcs_blk_n,                                          
     &         urcs_n_blocks(blk)%ptr(sig_offset), ngp, nstrs, span )           
      else                                                                      
        call tanstf_gastr( ddtse, eps_n_blocks(blk)%ptr(eps_offset),            
     &         ngp, nstr, span )                                                
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ou_gastr                     *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 03/9/13 rhd                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine ou_gastr( history_local, history_global, ngp,                  
     &                     mxhist, mxngp, hist_size, span, mxvl )               
      implicit none                                                             
c                                                                               
c               parameter declarations                                          
c                                                                               
      integer ngp, mxhist, mxngp, hist_size, mxvl, span                         
      double precision                                                          
     & history_local(mxvl,mxhist,mxngp),                                        
     & history_global(hist_size,ngp,span)                                       
c                                                                               
c               local declarations                                              
c                                                                               
      integer i, j, k                                                           
c                                                                               
      if( ngp .ne. 8 ) then                                                     
        do k = 1, ngp                                                           
         do  j = 1, hist_size                                                   
            do  i = 1, span                                                     
               history_local(i,j,k) = history_global(j,k,i)                     
            end do                                                              
         end do                                                                 
        end do                                                                  
        return                                                                  
      end if                                                                    
c                                                                               
c                number of gauss points = 8, unroll.                            
c                                                                               
      do  j = 1, hist_size                                                      
        do  i = 1, span                                                         
            history_local(i,j,1) = history_global(j,1,i)                        
            history_local(i,j,2) = history_global(j,2,i)                        
            history_local(i,j,3) = history_global(j,3,i)                        
            history_local(i,j,4) = history_global(j,4,i)                        
            history_local(i,j,5) = history_global(j,5,i)                        
            history_local(i,j,6) = history_global(j,6,i)                        
            history_local(i,j,7) = history_global(j,7,i)                        
            history_local(i,j,8) = history_global(j,8,i)                        
        end do                                                                  
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
