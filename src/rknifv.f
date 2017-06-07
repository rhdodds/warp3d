c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine rknifv                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 7/2/12 rhd                 *          
c     *                                                              *          
c     *     drives comptuation of internal force vectors for a       *          
c           block of similar elements. integral trans B * sigma      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine rknifv( eleifv, updated_element_volumes, nrow_ifv,             
     &                   local_work )                                           
     &                                                                          
      use segmental_curves, only : max_seg_points, max_seg_curves               
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
      include 'include_sig_up'                                                  
c                                                                               
c                       parameter declarations                                  
c                                                                               
      double precision                                                          
     &  eleifv(nrow_ifv,*), updated_element_volumes(*)                          
c                                                                               
c                       local variables                                         
c                                                                               
      double precision                                                          
     &  zero, xi, eta, zeta, element_volumes_for_blk(mxvl)                      
      logical geonl, bbar, local_debug                                          
      data local_debug, zero / .false., 0.0 /                                   
c                                                                               
      span            = local_work%span                                         
      felem           = local_work%felem                                        
      elem_type       = local_work%elem_type                                    
      order           = local_work%int_order                                    
      ngp             = local_work%num_int_points                               
      nnode           = local_work%num_enodes                                   
      ndof            = local_work%num_enode_dof                                
      geonl           = local_work%geo_non_flg                                  
      step            = local_work%step                                         
      iter            = local_work%iter                                         
      bbar            = local_work%bbar_flg                                     
      mat_type        = local_work%mat_type                                     
      totdof          = local_work%totdof                                       
      surf            = local_work%surface                                      
c                                                                               
c                       initiialize ifv's for the block.                        
c                                                                               
!DIR$ VECTOR ALIGNED                                                            
      eleifv(1:span,1:totdof)         = zero                                    
!DIR$ VECTOR ALIGNED                                                            
      element_volumes_for_blk(1:span) = zero                                    
c                                                                               
c                       compute all the shape function derivates, and           
c                       inverse jacobians.  also calculate volume               
c                       terms if using bbar. use element shape                  
c                       at n+1 for geonl.                                       
c                                                                               
      if( bbar .and. elem_type .eq. 2 )                                         
     &  call zero_vol( local_work%vol_block, local_work%volume_block,           
     &                 span, mxvl )                                             
c                                                                               
c             this subroutine is called only for a block of                     
c             cohesive elements. here the element coordinates and               
c             the global displacements are rotated to a coordinate              
c             system in which the normal axis (Z rotated) is                    
c             perpendicular to the surface ot the cohesive element              
c                                                                               
      if( local_work%is_cohes_elem ) then                                       
           call cohes_rot_mat( span, felem, nnode, elem_type,                   
     &                         local_work%ce_n1,                                
     &                         local_work%cohes_rot_block )                     
        if ( geonl )                                                            
     &      call cohes_mirror_refsurf( span, mxvl, totdof, nnode,               
     &                                 local_work%ce_n1 )                       
      end if                                                                    
c                                                                               
      do gpn = 1, ngp                                                           
        call getgpts( elem_type, order, gpn, xi, eta, zeta,                     
     &                local_work%weights(gpn) )                                 
        call derivs( elem_type, xi, eta, zeta, local_work%nxi(1,gpn),           
     &               local_work%neta(1,gpn),local_work%nzeta(1,gpn) )           
        call jacob1( elem_type, span, felem, gpn, local_work%jac,               
     &    local_work%det_j(1,gpn), local_work%gama(1,1,1,gpn),                  
     &    local_work%cohes_rot_block,                                           
     &    local_work%nxi(1,gpn), local_work%neta(1,gpn),                        
     &    local_work%nzeta(1,gpn), local_work%ce_n1, nnode )                    
c                                                                               
        if ( local_work%is_cohes_elem )                                         
     &     call shapef( elem_type, xi, eta, zeta,                               
     &                  local_work%shape(1,gpn) )                               
c                                                                               
        if ( bbar .and. elem_type .eq. 2 ) then                                 
          call vol_terms( local_work%gama(1,1,1,gpn),                           
     &                    local_work%det_j(1,gpn),                              
     &                    local_work%vol_block,                                 
     &                    local_work%nxi(1,gpn), local_work%neta(1,gpn),        
     &                    local_work%nzeta(1,gpn),                              
     &                    local_work%volume_block, span, mxvl )                 
        end if                                                                  
      end do                                                                    
c                                                                               
c                                                                               
      if ( bbar .and. elem_type .eq. 2 )                                        
     &  call vol_avg( local_work%vol_block, local_work%volume_block,            
     &                span, mxvl )                                              
c                                                                               
c                       compute ifv's, one gauss point at                       
c                       a time for all elements in this block.                  
c                       volumetric type elements have their current             
c                       volume computed and saved.                              
c                                                                               
      do gpn = 1, ngp                                                           
         call gpifv1( eleifv, span, gpn, local_work%weights(gpn),               
     &                local_work%det_j(1,gpn), local_work%b,                    
     &                local_work%urcs_blk_n1(1,1,gpn), local_work,              
     &                element_volumes_for_blk )                                 
      end do                                                                    
                                                                                
                                                                                
c      if(local_work%step .eq. 44 .and. felem .eq. 1) then                      
c       write(*,9200) felem                                                     
c       do i = 1, span                                                          
c         write(*,*) ' '                                                        
c         write(*,*) 'element: ', felem+i-1                                     
c         write(*,9300) eleifv(i,1:totdof)                                      
c       end do                                                                  
c      end if                                                                   
c 9200 format(/,2x,'... rknifv for block with first element: ',i6)              
c 9300 format(5x,8e14.6)                                                        
                                                                                
c                                                                               
c                                                                               
c                       transform element internal force                        
c                       vectors to constraint compatible coordinates.           
c                                                                               
      if ( local_work%trn_e_block ) then                                        
         do i = 1, span                                                         
           if (  local_work%trn_e_flags(i) ) then                               
             call trnvecs( eleifv, local_work%trnmte, local_work%trne,          
     &                     ndof, nnode, i, 1, span )                            
           end if                                                               
         end do                                                                 
      end if                                                                    
c                                                                               
c                       complete processing updated element volumes.            
c                       if element has been killed, don't change its            
c                       volume. we freeze the volume as the value               
c                       when element is killed.                                 
c                                                                               
      if ( .not. local_work%is_cohes_elem ) then                                
      if ( local_debug )                                                        
     &    write(*,*) '.. calling update_element_volumes ..'                     
       call update_element_volumes( felem, span,                                
     &                               updated_element_volumes,                   
     &                               element_volumes_for_blk )                  
      end if                                                                    
c                                                                               
      go to 9999                                                                
c                                                                               
c                                                                               
 9999 continue                                                                  
      if ( local_debug ) write(local_work%iout,*)                               
     &     "    >>> leaving rknifv ...."                                        
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine update_element_volumes               *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 01/08/2010                 *          
c     *                     removed increment of dam_state [bug]     *          
c     *                                                              *          
c     *     this subroutine saves the volume for elements not        *          
c     *     killed                                                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
                                                                                
      subroutine update_element_volumes(                                        
     &  felem, span, updated_element_volumes, element_volumes_for_blk )         
      use global_data ! old common.main
c                                                                               
      use elem_extinct_data, only : dam_state                                   
      use damage_data, only : dam_ptr, max_dam_state, growth_by_kill            
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                       parameter declarations                                  
c                                                                               
      double precision                                                          
     & updated_element_volumes(*), element_volumes_for_blk(*)                   
c                                                                               
c                       local variables                                         
c                                                                               
      logical update_volume(mxvl), local_debug                                  
      data local_debug / .false. /                                              
c                                                                               
c                       if crack growth by element killing is                   
c                       not being used, all element volumes are                 
c                       updated.                                                
c                                                                               
      if ( local_debug )                                                        
     &   write(*,*)                                                             
     &       '... in update_ele...growth_by_killm felem, span: ',               
     &       growth_by_kill, felem, span                                        
      if ( .not. growth_by_kill ) then                                          
!DIR$ VECTOR ALIGNED                                                            
        do i = 1, span                                                          
           updated_element_volumes(i) = element_volumes_for_blk(i)              
        end do                                                                  
        return                                                                  
      end if                                                                    
c                                                                               
c                       don't update volume for killed elements.                
c                       build list of ones to do.                               
c                                                                               
      do i = 1, span                                                            
        update_volume(i) = .false.                                              
      end do                                                                    
c                                                                               
      if( local_debug )  write(*,*)                                             
     & '      .. running over block. felem, span: ', felem,span                 
      do i = 1, span                                                            
       abselem     = felem + i - 1                                              
       elem_ptr    = dam_ptr(abselem)                                           
c                                                                               
c                       zero elem_ptr means element is not killable             
c                                                                               
       if ( elem_ptr .eq. 0 ) then                                              
         update_volume(i) = .true.                                              
         cycle                                                                  
       end if                                                                   
c                                                                               
c                       if already killed, skip                                 
c                                                                               
       if ( dam_state(elem_ptr) .gt. 0 ) cycle                                  
c                                                                               
c                       element is killable but not yet killed.                 
c                                                                               
       update_volume(i) = .true.                                                
c                                                                               
      end do                                                                    
c                                                                               
c                       update volume of non-killed elements.                   
!DIR$ VECTOR ALIGNED                                                            
      do i = 1, span                                                            
       if( update_volume(i) )  updated_element_volumes(i) =                     
     &                         element_volumes_for_blk(i)                       
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
