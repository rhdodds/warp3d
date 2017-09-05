c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oupstr_elem                  *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 8/20/2017 rhd              *          
c     *                                                              *          
c     *  drive output of stress or strain element results to         *          
c     *  (1) a Patran file in either binary or formatted forms       *          
c     *  or (2) a flat file with text or stream format               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oupstr_elem( stress, oubin, ouasc, flat_file,                  
     &                      stream_file, text_file, compressed )                
      use global_data ! old common.main
      use main_data, only : cohesive_ele_types, bar_types, link_types                                  
      use elblk_data, only : urcs_blk_n  ! readonly here for debugging          
                                                                                
      implicit none                                                             
c                                                                               
      logical :: stress, oubin, ouasc, flat_file,                              
     &           stream_file, text_file, compressed                             
c                                                                               
c                       locally allocated. the big one is                       
c                       elem_results.                                           
c                                                                               
      integer :: i, blk, span, felem, elem_type, int_order, mat_type,           
     &           num_enodes, num_enode_dof, totdof, num_int_points,             
     &           output_loc, num_short_strain,  num_short_stress,               
     &           num_vals, iblk, ifelem                                         
      integer :: elem_out_map(mxelmp)                                           
                                                                                
      logical ::  bbar_flg, geo_non_flg, long_out_flg, nodpts_flg,              
     &            center_output, first_block, last_block, cohesive_elem,        
     &            local_debug, is_bar_elem, is_link_elem                                                   
                                                                                
      integer :: elem, j                                                        
c                                                                               
      double precision, parameter :: zero = 0.d0, small_tol = 1.d-80            
      double precision :: elem_results(mxvl,mxstmp)                             
c                                                                               
      first_block = .true.                                                      
      last_block  = .false.                                                     
      local_debug = .false.                                                     
c                                                                               
c                       compute the stress/strain element output                
c                       vectors and assemble them into the global               
c                       element results for each block of similar,              
c                       non-conflicting elements.                               
c                                                                               
      do blk = 1, nelblk                                                        
         if( elblks(2,blk) .ne. myid ) cycle                                    
         elem_results      = zero    ! all entries                              
         span              = elblks(0,blk)                                      
         felem             = elblks(1,blk)                                      
         elem_type         = iprops(1,felem)                                    
         int_order         = iprops(5,felem)                                    
         mat_type          = iprops(25,felem)                                   
         num_enodes        = iprops(2,felem)                                    
         num_enode_dof     = iprops(4,felem)                                    
         totdof            = num_enodes * num_enode_dof                         
         geo_non_flg       = lprops(18,felem)                                   
         bbar_flg          = lprops(19,felem)                                   
         num_int_points    = iprops(6,felem)                                    
         long_out_flg      = lprops(16,felem)                                   
         output_loc        = iprops(12,felem)                                   
         nodpts_flg        = .true.                                             
         center_output     = .true.                                             
         num_short_stress  = 11                                                 
         num_short_strain  = 7                                                  
         num_vals          = num_short_strain + 15                              
         if( stress ) num_vals = num_short_stress + 15                          
         cohesive_elem     = cohesive_ele_types(elem_type)  
         is_bar_elem       = bar_types(elem_type)
         is_link_elem      = link_types(elem_type)
         if( local_debug ) write(*,*) '.. block, span, felem: ',                
     &                      blk, span, felem                                    
c                                                                               
c                       for cohesive, bar & link elements, output 
c                       zeroes for element result values.                                  
c                                                                               
         if( cohesive_elem .or. is_bar_elem .or. is_link_elem ) then                                               
           call oust_elem( stress, oubin, ouasc, num_vals,                      
     &                     elem_results, mxvl, span, first_block,               
     &                     felem, last_block, flat_file,                        
     &                     stream_file, text_file, compressed )                 
           first_block = .false.                                                
           cycle                                                                
         end if                                                                 
c                                                                               
c                       process all solid element types                         
c                                                                               
         elem_out_map(1:mxelmp) = outmap(elem_type,1:mxelmp)                    
c                                                                               
c                       duplicate necessary element block data.                 
c                                                                               
         call oudups( span, felem, num_int_points, geo_non_flg, stress,         
     &                  .false. )                                               
                                                                                
         if( local_debug ) then                                                 
            do i =1, span                                                       
                elem = felem + i -1                                             
                write(*,*) '... GP results for element: ', elem                 
                do j = 1, 8  ! gp number                                        
                    write(*,9000)  j, urcs_blk_n(i,1:6,j)                       
                end do                                                          
            end do                                                              
         end if                                                                 
                                                                                
c                       compute the element block center stress/strain          
c                       data.                                                   
c                                                                               
         iblk = blk                                                             
         ifelem = felem                                                         
         call ouprks( span, iblk, ifelem, elem_type, int_order,                 
     &                num_int_points, num_enodes, geo_non_flg,                  
     &                stress, mat_type, center_output,                          
     &                num_short_stress, num_short_strain, .true. )              
c                                                                               
c                       copy element block element stress/strain data           
c                       into the global element results array.                  
c                                                                               
         call oupele( span, num_short_strain, num_short_stress,                 
     &                stress, elem_results(1,1), mxvl )                         
c                                                                               
c                       calculate the extended values of stress and             
c                       strain for each element                                 
c                                                                               
         call ouext2( elem_results, mxvl, span, stress )                        
c                                                                               
c                       output the element results                              
c                       to a file compatable with patran for post               
c                       processing.                                             
c                                                                               
c                       zero small values to prevent 3-digit exponents          
c                       in formated output                                      
c                                                                               
         where( abs(elem_results) .lt. small_tol ) elem_results = zero          
c                                                                               
         call oust_elem( stress, oubin, ouasc, num_vals, elem_results,          
     &                   mxvl, span, first_block, felem, last_block,            
     &                   flat_file,stream_file, text_file,                      
     &                   compressed )                                           
         first_block = .false.                                                  
c                                                                               
      end do                                                                    
c                                                                               
c                                                                               
      where( abs(elem_results) .lt. small_tol ) elem_results = zero             
c                                                                               
      last_block = .true.                                                       
      call oust_elem( stress, oubin, ouasc, num_vals, elem_results,             
     &                mxvl, span, first_block, felem, last_block,               
     &                flat_file, stream_file, text_file, compressed )           
c                                                                               
      return                                                                    
c                                                                               
 9000 format(5x,i2,6e14.6)                                                      
c                                                                               
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
