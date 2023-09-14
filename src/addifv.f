c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine addifv                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 8/25/23 rhd                *          
c     *                                                              *          
c     *     assembles the internal force vectors for a block of      *          
c     *     similar, elements into the global internal force vector. *          
c     *     this is a classic scatter operation                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine addifv( span, bedst, totdof, ifv,                      
     &                   felem, sum_ifv, num_term_ifv,                          
     &                   eleifv )
      use global_data, only : out, mxelpr, mxedof, iprops                      
      use damage_data, only : dam_ptr, growth_by_kill,
     &                        use_mesh_regularization,
     &                        tol_regular =>
     &                        tolerance_mesh_regularization  
      use elem_extinct_data, only : dam_ifv, dam_state,
     &                        smcs_d_values
      use constants
c
      implicit none                                                    
c                     
      integer :: span, totdof, felem, num_term_ifv
      integer :: bedst(totdof,*)
      double precision :: ifv(*), sum_ifv, eleifv(span,*)                     
      logical, parameter :: debug = .false.                                                          
c                                                                               
c             locals                                                            
c   
      logical :: elems_in_blk_killable, standard_process
c                                                                                 
c      if( debug ) write (out,*) '>>>>  inside addifv'                           
c
c              with crack growth by element death, more logic
c              is required. if 1st element in block is killable
c              they are all killable.
c
      standard_process = .true.
c
      If( growth_by_kill ) then  ! is blk of killable elems?
        elems_in_blk_killable = dam_ptr(felem) > 0
        if( elems_in_blk_killable ) standard_process = .false.
      end if
c
      if( standard_process ) then 
        call addifv_standard ! no killable elements in block
        return
      end if 
c
c              solution uses crack growth by element deletion &
c              this block has killable elements.
 
      if( use_mesh_regularization ) then 
          call addifv_regularization
          return
      end if
      if( .not. use_mesh_regularization ) then
          call addifv_no_regularization
          return
      end if 
c
      write(out,9000) 
      call die_abort
c
      return
9000  format('>> FATAL ERROR: addifv. inconsistent condition.',
     & ' job terminated.'//)
      
      contains
c     ========                                                    
c     ****************************************************************          
c     *                                                              *          
c     *               internal subroutine addifv_standard            *          
c     *                                                              *          
c     *                   last modified : 3/6/2021 rhd               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine addifv_standard
      implicit none
c
      integer :: i, j
c
c              no killable elements in block. simple scatter update
c
       do j = 1, totdof ! for each element in blk                                                         
         do i = 1, span                                                         
c$OMP ATOMIC UPDATE                                                             
            ifv(bedst(j,i)) = ifv(bedst(j,i)) + eleifv(i,j)                     
         end do                                                                 
      end do 
c           
c             thread private value passed for sum_ifv, num_term_ifv
c                                                 
      do j = 1, totdof                                                          
!DIR$ IVDEP                                                                     
         do i = 1, span                                                         
            sum_ifv = sum_ifv + abs(eleifv(i,j))                                
         end do                                                                 
      end do                                                                    
      num_term_ifv = num_term_ifv + totdof * span                               
c                                                                               
      if( debug ) then                                                          
       write(out,9200) felem                                                    
       do i = 1, span                                                           
         write(out,*) ' '                                                       
         write(out,*) 'element: ', felem+i-1                                    
         write(out,9300) eleifv(i,1:totdof)                                     
       end do                                                                   
      end if   
c
      return
c
 9200 format(/,2x,'... addifv for block with first element: ',i7)               
 9300 format(5x,8e14.6)     
c
      end subroutine  addifv_standard      
c     ****************************************************************          
c     *                                                              *          
c     *       internal subroutine addifv_no_regularization           *          
c     *                                                              *          
c     *                   last modified : 3/6/2021 rhd               *          
c     *                                                              *          
c     ****************************************************************          
c
      subroutine addifv_no_regularization
c
      implicit none
c     
      integer :: i, j, relem, element, elem_ptr
c
c              elements in block that are being killed have eleifv = 0.
c              for killable but not yet killed, we save eleifv 
c              for subsequent reduction. need to save since we don't know
c              yet when element will be killed.
c
      do j = 1, totdof ! for each element in blk                                                         
         do i = 1, span  ! simple scatter update                                                       
c$OMP ATOMIC UPDATE                                                             
            ifv(bedst(j,i)) = ifv(bedst(j,i)) + eleifv(i,j)                     
         end do                                                                 
      end do 
c           
c             thread private value passed for sum_ifv, num_term_ifv
c                                                 
      do j = 1, totdof                                                          
!DIR$ IVDEP                                                                     
         do i = 1, span                                                         
            sum_ifv = sum_ifv + abs(eleifv(i,j))                                
         end do                                                                 
      end do                                                                    
      num_term_ifv = num_term_ifv + totdof * span                               
c                                                                               
      if( debug ) then                                                          
       write(out,9200) felem                                                    
       do i = 1, span                                                           
         write(out,*) ' '                                                       
         write(out,*) 'element: ', felem+i-1                                    
         write(out,9300) eleifv(i,1:totdof)                                     
       end do                                                                   
      end if                                                                    
c                                                                               
c             if we have NOT started releasing the internal forces,             
c             store the element contribution to the ifv in the                  
c             crack growth data structures.                                     
c             this is not real efficient but we don't know here                 
c             which elements have just been killed.                             
c                                                                               
      do relem = 1, span                                                            
         element  = felem + relem - 1       
         elem_ptr =  dam_ptr(element)    
         if( dam_state(elem_ptr) > 0 ) cycle  ! already killed       
!DIR$ IVDEP   
         do i = 1, totdof
            dam_ifv(i,elem_ptr) = eleifv(relem,i)   
         end do                        
      end do                                                                    
c                                                                               
      if( debug ) write(out,*) '<<<<  leaving  addifv_no_regularization'                         
      return                                                                    
c                                                                               
 9000 format(1x,'dof:',i2,' element ifv:',e14.6,' total ifv:',e14.6)            
 9100 format(1x,'FATAL ERROR: in addifv_no_regularization. contact ',
     &   'WARP3D developers. Job terminated'//)             
 9200 format(/,2x,'... addifv for block with first element: ',i7)               
 9300 format(5x,8e14.6)     
                                                        
      end subroutine addifv_no_regularization       

c     ========                                                    
c     ****************************************************************          
c     *                                                              *          
c     *        internal subroutine addifv_regularization             *          
c     *                                                              *          
c     *                   last modified : 8/25/2023 rhd              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine addifv_regularization
c
      implicit none
c
      integer :: i, j, relem, element, elem_ptr
      logical, parameter :: local_debug = .false.
c
c              element strains, stresses, internal forces were 
c              computed assuming *no* damage
c
c              all elements in block are killable. those not yet in
c              death have dam_state = 0. zero ifv for elements in 
c              killing or killed state. then do a regular scatter 
c              update
c
      if( local_debug ) write(out,*) 
     &          "@10.5 entered addifv_regularization "           
c
c             if we have NOT started releasing the internal forces,             
c             store the element contribution to the ifv in the                  
c             crack growth data structures.                                     
c             this is not real efficient but we don't know here                 
c             which elements have just been killed.    
c             for already in death, zero anyway for safety in case a bug                         
c                 
      do relem = 1, span                                                            
         element  = felem + relem - 1       
         elem_ptr =  dam_ptr(element)    
         if( dam_state(elem_ptr) .eq. 0  ) then! not yet killed 
!DIR$ IVDEP   
           do i = 1, totdof
              dam_ifv(i,elem_ptr) = eleifv(relem,i)   
           end do   
         end if                                
         if( dam_state(elem_ptr) > 0 ) then
!DIR$ IVDEP   
           do i = 1, totdof
              eleifv(relem,i) = zero
           end do   
         end if        
      end do  ! relem                                                                   
c               
       do j = 1, totdof ! for each element in blk                                                         
         do i = 1, span                                                         
c$OMP ATOMIC UPDATE                                                             
            ifv(bedst(j,i)) = ifv(bedst(j,i)) + eleifv(i,j)                     
         end do                                                                 
      end do 
c           
c             thread private value passed for sum_ifv, num_term_ifv
c                                                 
      do j = 1, totdof                                                          
!DIR$ IVDEP                                                                     
         do i = 1, span                                                         
            sum_ifv = sum_ifv + abs(eleifv(i,j))                                
         end do                                                                 
      end do                                                                    
      num_term_ifv = num_term_ifv + totdof * span                               
c                                                                               
      if( debug ) then                                                          
       write(out,9200) felem                                                    
       do i = 1, span                                                           
         write(out,*) ' '                                                       
         write(out,*) 'element: ', felem+i-1                                    
         write(out,9300) eleifv(i,1:totdof)                                     
       end do                                                                   
      end if   
c
      return
c
 9100 format(1x,'FATAL ERROR: in addifv_regularization. contact ',
     &   'WARP3D developers. Job terminated'//)             
 9200 format(/,2x,'... addifv for block with first element: ',i7)               
 9300 format(5x,8e14.6)
c
      end subroutine addifv_regularization   

      end subroutine addifv                                                 
