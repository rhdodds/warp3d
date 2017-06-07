c                                                                               
c                                                                               
c        mm08 is the placeholder for the UMAT material                          
c        interface in WARP3D.                                                   
c                                                                               
c        internally, WARP3D processes the UMAT as material                      
c        model #8 but the the special name UMAT rather than                     
c        mm08                                                                   
c                                                                               
c        the UMAT processes 1 integration point per call rather                 
c        than a block of elements for the integration point as                  
c        the other WARP3D models. the code in stiffness and stress              
c        update handles the looping over elements in the block to hide          
c        the element blocking from the UMAT                                     
c                                                                               
c        actual code here is called by WARP3D to support UMAT                   
c        operations, e.g., drive to get states values for output                
c                                                                               
c                                                                               
      subroutine mm08  !   never called by WARP3D                               
      return                                                                    
      end                                                                       
                                                                                
c                                                                               
c       routines umat_output_states and umat_states_labels                      
c       must be defined in the user_routine.f file                              
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *     subroutine mm08_states_labels   (warp3d model 8 )        *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *               last modified : 1/11/2015 (rhd)                *          
c     *                                                              *          
c     ****************************************************************          
                                                                                
                                                                                
      subroutine mm08_states_labels( size_state,                                
     &      num_states, state_labels, state_descriptors, out,                   
     &      comment_lines, max_comment_lines, num_comment_lines )               
      implicit none                                                             
c                                                                               
c                       parameters                                              
c                                                                               
      integer :: size_state, num_states, out, max_comment_lines,                
     &           num_comment_lines                                              
      character(len=8)  :: state_labels(size_state)                             
      character(len=60) :: state_descriptors(size_state)                        
      character(len=80) :: comment_lines(max_comment_lines)                     
c                                                                               
c                       locals                                                  
c                                                                               
c                                                                               
c                       call the UMAT routine to get labels                     
c                                                                               
      call umat_states_labels( size_state,                                      
     &      num_states, state_labels, state_descriptors, out,                   
     &      comment_lines, max_comment_lines, num_comment_lines )               
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *     subroutine mm08_states_values     (warp3d model 8 )      *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *               last modified : 12/14/2014 (rhd)               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm08_states_values( itype, elem_states_output,                 
     &                                 nrow_states, num_states  )               
      use global_data ! old common.main
c                                                                               
c                       access some global data structures                      
c                                                                               
      use elem_block_data, only: history_blocks, history_blk_list               
      use main_data, only: elems_to_blocks                                      
      implicit integer (a-z)                                                    
c                                                                               
c                       parameters                                              
c                                                                               
      integer :: nrow_states, itype, num_states                                 
      double precision :: elem_states_output(nrow_states,*)                     
c                                                                               
c                       locals                                                  
c                                                                               
      double precision,                                                         
     & allocatable :: history_dump(:,:,:), one_elem_states(:),                  
     &                avgs(:), ip_state_values(:)                               
      integer :: relem, elnum, hist_size, blockno                               
      logical :: do_a_block, local_debug                                        
      double precision :: zero                                                  
      data zero / 0.0d00 /                                                      
c                                                                               
c           build umat states values output.                                    
c                                                                               
c              itype > 0 => this is the block number. do all elements           
c                           in the block                                        
c                                                                               
c              itype < 0 => this is an element number. put state                
c                           values into column 1 of results.                    
c                                                                               
      do_a_block = .true.                                                       
      if( itype. gt. 0 ) then                                                   
         do_a_block = .true.                                                    
         blockno = itype                                                        
      else                                                                      
         do_a_block = .false.                                                   
         elnum = -itype                                                         
         blockno = elems_to_blocks(elnum,1)                                     
      end if                                                                    
c                                                                               
      local_debug = .false.                                                     
      felem       = elblks(1,blockno)                                           
      elem_type   = iprops(1,felem)                                             
      mat_type    = iprops(25,felem)                                            
      int_points  = iprops(6,felem)                                             
      span        = elblks(0,blockno)                                           
      hist_size   = history_blk_list(blockno)                                   
      if( local_debug ) write(out,9050) blockno, felem, elem_type,              
     &         mat_type, int_points, span, hist_size                            
c                                                                               
c           temporary block of history so it can be re-organized.               
c           working vectors.                                                    
c                                                                               
      allocate( one_elem_states(nrow_states) )                                  
      allocate( avgs(nrow_states), ip_state_values(nrow_states) )               
      allocate( history_dump(hist_size,int_points,span) )                       
c                                                                               
      history_dump = reshape( history_blocks(blockno)%ptr,                      
     &           (/hist_size,int_points,span/) )                                
c                                                                               
      if( do_a_block ) then                                                     
        do relem = 1, span                                                      
           elnum = felem + relem - 1  ! absolute element number                 
           one_elem_states(1:num_states) = zero                                 
           call mm08_states_values_a                                            
           elem_states_output(1:num_states,relem) =                             
     &                one_elem_states(1:num_states)                             
        end do                                                                  
      else                                                                      
           relem = elnum + 1 - felem                                            
           one_elem_states(1:num_states) = zero                                 
           call mm08_states_values_a                                            
           elem_states_output(1:num_states,1) =                                 
     &                one_elem_states(1:num_states)                             
      end if                                                                    
c                                                                               
      deallocate( history_dump, one_elem_states, avgs,                          
     &            ip_state_values  )                                            
c                                                                               
      return                                                                    
c                                                                               
 9050 format(10x,"block, felem, etype, mtype:  ",4i7,                           
     &  /,10x,   "int_pts, span, hist_size:    ",3i7 )                          
c                                                                               
      contains                                                                  
c     ========                                                                  
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *             subroutine mm08_states_values_a                  *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 12/15/2014 (rhd)           *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm08_states_values_a                                           
c                                                                               
      implicit none                                                             
c                                                                               
c                       locals                                                  
c                                                                               
      integer :: ipt, num_states_returned                                       
      double precision ::                                                       
     & a_bar, b_bar                                                             
c                                                                               
      avgs(1:num_states) = zero                                                 
c                                                                               
      do ipt = 1, int_points                                                    
        call umat_output_states( history_dump(1,ipt,relem),                     
     &                   ip_state_values, num_states_returned )                 
        if( num_states .ne. num_states_returned ) then                          
            write(out,9000) elnum                                               
            call die_abort                                                      
        end if                                                                  
        avgs(1:num_states) = avgs(1:num_states) +                               
     &                       ip_state_values(1:num_states)                      
      end do                                                                    
c                                                                               
      one_elem_states(1:num_states) = avgs(1:num_states) /                      
     &                                dble(int_points)                          
c                                                                               
      return                                                                    
c                                                                               
 9000 format(/1x,                                                               
     &'>>>>> Error: umat routines did not return correct number of '            
     & /,14x,'state variables for element: ',i7,                                
     & /,14x,'detected by routine oustates_values_umat_a',                      
     & /,14x,'job terminated....'/)                                             
                                                                                
      end subroutine mm08_states_values_a                                       
      end subroutine mm08_states_values                                         
                                                                                
                                                                                
                                                                                
