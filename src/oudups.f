c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oudups                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 8/4/25 rhd                 *          
c     *                                                              *          
c     *     gathers data for a block of elements to support          *          
c     *     generation of a patran, packet or hard copy output.      *  
c     *                                                              *
c     *     zero results for element being or already killed         *        
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oudups( span, felem, ngp, geonl, do_stress, 
     &                   is_cohesive )  
c       
      use global_data, only : mxvl, nstr, nstrs, out
      use main_data, only : elems_to_blocks                                     
      use elem_block_data, only : history_blocks, rot_n1_blocks,                
     &                            eps_n_blocks, urcs_n_blocks,                  
     &                            history_blk_list                              
      use elblk_data, only : elem_hist, blk_size_hist, urcs_blk_n,              
     &                       rot_blk_n1, ddtse, blk_size_gp  
      use damage_data, only : dam_ptr, growth_by_kill,
     &                        use_mesh_regularization
      use elem_extinct_data, only : dam_state
      use constants
c                                                                               
      implicit none                                                             
c                                                                               
      integer, intent(in) :: span, felem, ngp                                               
      logical, intent(in) :: geonl, do_stress, is_cohesive                                     
c                                                                               
c             local declarations                                                
c                                                                               
      integer :: blk, rel_elem, hist_size, hist_offset, rot_offset,             
     &           eps_offset, sig_offset, mxhist, mxngp, relem,
     &           element, elem_ptr                          
      logical :: is_solid, process_blk, do_strains                                                       
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
      do_strains  = .not. do_stress                                       
c                                                                               
      if( local_debug ) then                                                    
         write(out,*) '..... oudups....'                                        
         write(out,*) '....    span, felem, ngp, geonl, stress'                 
         write(out,*) span, felem, ngp, geonl, do_stress                           
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
      if( do_stress ) call tanstf_gastr( urcs_blk_n,                                          
     &         urcs_n_blocks(blk)%ptr(sig_offset), ngp, nstrs, span )           
      if( do_strains ) call tanstf_gastr( ddtse, 
     &            eps_n_blocks(blk)%ptr(eps_offset), ngp, nstr, span )                                                
c  
c              for killed elements in block, zero local copy of results.
c
      process_blk = is_solid .and. growth_by_kill .and.
     &              dam_ptr(felem) > 0 .and. use_mesh_regularization
      if( .not. process_blk ) return
c
      do relem = 1, span                                                            
         element  = felem + relem - 1       
         elem_ptr = dam_ptr(element)    
         if( dam_state(elem_ptr) == 0 ) cycle ! element not yet killing
         call oudups_1( relem, elem_hist(1,1,1), ngp, 
     &                  mxhist, mxngp, hist_size, span, mxvl    )
         if( do_stress) call oudups_2( relem, urcs_blk_n, ngp, 
     &                              nstrs, span, mxvl )      
         if( do_strains ) call oudups_2( relem, ddtse, ngp, 
     &                                     nstr, span, mxvl )
      end do
c
      return
c
      contains              
c     ========
c     ****************************************************************          
c     *                                                              *          
c     *              internal subroutine oudups_1                    *          
c     *                                                              *          
c     ****************************************************************          
c
      subroutine oudups_1( relem, history_local, ngp,                  
     &                     mxhist, mxngp, hist_size, span, mxvl )               
      implicit none                                                             
c                                                                               
c               parameter declarations                                          
c                                                                               
      integer, intent(in) :: relem, ngp, mxhist, mxngp, hist_size, 
     &                       mxvl, span                         
      double precision, intent(out) :: history_local(mxvl,mxhist,mxngp)
c                                                                               
c               local declarations                                              
c                                                                               
      integer :: j, k                                                           
c                                                                               
      do k = 1, ngp                                                           
         do  j = 1, hist_size                                                   
            history_local(relem,j,k) = zero
         end do                                                                 
      end do    
c                                                                    
      return
      end subroutine oudups_1
c        
c     ****************************************************************          
c     *                                                              *          
c     *              internal subroutine oudups_2                    *          
c     *                                                              *          
c     ****************************************************************          
c
      subroutine oudups_2( relem, ml, ngp, nprm, span, mxvl )                        
      implicit none                                                             
c                                                                               
c               parameter declarations                                          
c                                                                               
      integer :: relem, ngp, nprm, span, mxvl                                                
      double precision :: ml(mxvl,nprm,ngp)
c                                                                               
      integer :: j, k                                                        
c                                                                               
      do k = 1, ngp                                                           
        do  j = 1, nprm                                                        
          ml(relem,j,k) = zero
        end do                                                              
      end do                                                                 
c
      return       
      end subroutine oudups_2
c
      end subroutine oudups                                                                 

                                                                  
