c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oupstr_elem                  *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 9/18/2025 rhd              *          
c     *                                                              *          
c     *  drive output of stress or strain element results to         *          
c     *  (1) a Patran file in either binary or formatted forms       *          
c     *  or (2) a flat file with text or stream format               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oupstr_elem( do_stress, oubin, ouasc, flat_file,                  
     &                        stream_file, text_file, compressed )
c     
      use global_data, only : nelblk, elblks, myid, iprops, lprops,
     &                        outmap, mxelmp, mxstmp, mxvl
      use main_data, only : cohesive_ele_types, bar_types, link_types                                  
      use elblk_data, only : urcs_blk_n  ! readonly here for debugging      
      use constants    
      use output_value_indexes, only : num_short_strain, 
     &            num_short_stress, num_long_strain, num_long_stress
                                                                                
      implicit none                                                             
c                                                                               
      logical :: do_stress, oubin, ouasc, flat_file,                              
     &           stream_file, text_file, compressed                             
c                                                                               
c                       locally allocated                                           
c                                                                               
      integer :: i, blk, span, felem, elem_type, int_order, mat_type,           
     &           num_enodes, num_enode_dof, totdof, num_int_points,             
     &           output_loc, num_vals, iblk, ifelem                                         
      integer :: elem_out_map(mxelmp)                                           
      logical :: bbar_flg, geo_non_flg, long_out_flg, nodpts_flg,              
     &           center_output, first_block, last_block, cohesive_elem,        
     &           local_debug, is_bar_elem, is_link_elem, do_strains,
     &           do_elem_output                                                  
                                                                                
      integer :: elem, j                                                        
c                                                                               
      double precision, parameter :: small_tol = 1.d-80            
      double precision :: elem_results(mxvl,mxstmp)                             
c         
      do_strains  = .not. do_stress    
      do_elem_output = .true.                                                                  
      first_block = .true.                                                      
      last_block  = .false.                                                     
      local_debug = .false.                                                     
c                                                                               
c                       process all elements by blocks. for threads only,
c                       the order of elements is sequential over blocks.
c                       can output to file block by block.
c
c                       for MPI, process blocks owned by this domain.
c                       an external program must be run to combine block
c                       results into element sequential order before
c                       writing the final results file.
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
         if( do_strains ) num_vals = num_long_strain                              
         if( do_stress )  num_vals = num_long_stress                          
         cohesive_elem     = cohesive_ele_types(elem_type)  
         is_bar_elem       = bar_types(elem_type)
         is_link_elem      = link_types(elem_type)
         if( local_debug ) write(*,*) '.. block, span, felem: ',                
     &                      blk, span, felem   
c                                                                               
c                                                                               
c                       for cohesive, bar & link elements, output 
c                       zeroes for element result values.                                  
c                                                                               
         if( cohesive_elem .or. is_bar_elem .or. is_link_elem ) then                                               
           call oust_elem( do_stress, oubin, ouasc, num_vals,                      
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
c                       duplicate necessary element block data. these
c                       are put into elestr array in module elblk_data.  
c                       a working array shared by lower level routines
c                       thru the module. elestr is mxvl x max output
c                       values supported by code                             
c                                                                               
         call oudups( span, felem, num_int_points, geo_non_flg, 
     &                do_stress, .false. )   
c                                            
c                       compute element center stress/strain values         
c                       for elements in block. values put into elestr.
c                                                                               
         iblk = blk                                                             
         ifelem = felem 
         call ouprks( span, iblk, ifelem, elem_type, int_order,                 
     &                num_int_points, num_enodes, geo_non_flg,                  
     &                do_stress, mat_type, center_output,
     &                do_elem_output )              
c                                                                               
c                       copy elestr values into local elem_results
c                       array for the block. Use of local elem_results
c                       and "hidden" elestr could have been designed
c                       better. 
c                                                                               
         call oupele( span, num_short_strain, num_short_stress,                 
     &                do_stress, elem_results(1,1), mxvl )                         
c                                                                               
c                       calculate the extended (long) values of
c                       stress/strain for each element                                 
c                                                                               
         call ouext2( elem_results, mxvl, span, do_stress )                        
c                                                                               
c                       output the element results for block                            
c                       to a file compatable with patran or our flat
c                       file structure
c                                                                               
c                       zero small values to prevent 3-digit exponents          
c                       in formated output                                      
c                                                                               
         where( abs(elem_results) .lt. small_tol ) elem_results = zero          
c                                                                               
         call oust_elem( do_stress, oubin, ouasc, num_vals,         
     &                   elem_results,  mxvl, span, first_block, 
     &                   felem, last_block,            
     &                   flat_file, stream_file, text_file,                      
     &                   compressed )                                           
         first_block = .false.                                                  
c                                                                               
      end do                                                                    
c                                                                                                                                                              
c                       close the patran or flat results file. 
c                                                                               
      where( abs(elem_results) .lt. small_tol ) elem_results = zero             
c                                                                               
      last_block = .true.                                                       
      call oust_elem( do_stress, oubin, ouasc, num_vals, elem_results,             
     &                mxvl, span, first_block, felem, last_block,               
     &                flat_file, stream_file, text_file, compressed )           
c                                                                               
      return                                                                    
c                                                                               
 9000 format(5x,i2,6e14.6)                                                      
c                                                                               
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
