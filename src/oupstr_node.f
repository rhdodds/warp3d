c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oupstr_node                  *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 1/29/2017 rhd              *          
c     *                                                              *          
c     *     drives output of stress or strain nodal results to       *          
c     *     (1) patran files in either binary or formatted forms or  *          
c     *     (2) flat text or stream files                            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oupstr_node( stress, oubin, ouasc, flat_file,                  
     &                      stream_file, text_file, compressed  )               
      use global_data ! old common.main
      use main_data, only : incmap, incid, cohesive_ele_types                   
      implicit integer (a-z)                                                    
      logical stress, oubin, ouasc, flat_file,                                  
     &        stream_file, text_file, compressed                                
c                                                                               
c                local declarations                                             
c                                                                               
      logical  bbar_flg, geo_non_flg, long_out_flg, nodpts_flg,                 
     &         center_output, cohesive_elem, do_average, serial                 
      integer elem_out_map(mxelmp)                                              
c                                                                               
c                data structure for averaged nodal results. vector              
c                of derived types. global vector allocated for all              
c                structure nodes. derived type allocated only for               
c                nodes touched by elements in this domain.                      
c                                                                               
c                patran nodal result file contains values only                  
c                for nodes touched by elements in this domain.                  
c                they are averaged by the number of elements that               
c                touch then in this domain.                                     
c                                                                               
      type :: node_entry                                                        
         integer :: count                                                       
         double precision, dimension(:), pointer :: node_values                 
      end type node_entry                                                       
c                                                                               
      type (node_entry), dimension (:), pointer :: nodal_values                 
         double precision, dimension(:), pointer :: snode_values                
c                                                                               
c                                                                               
c                       create structure size vector of derived types.          
c                       set count for each node to zero. nullify                
c                       pointers to all nodal vectors.                          
c                                                                               
      allocate( nodal_values(nonode) )                                          
      do node = 1, nonode                                                       
       nodal_values(node)%count = 0                                             
       nullify( nodal_values(node)%node_values )                                
      end do                                                                    
c                                                                               
c                       compute the stress/strain nodal output                  
c                       vectors and assemble them into the global               
c                       nodal results for each element block                    
c                                                                               
      do blk = 1, nelblk                                                        
c                                                                               
         if ( elblks(2,blk) .ne. myid ) cycle                                   
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
         center_output     = output_loc .eq. 3                                  
         num_short_stress  = 11                                                 
         num_short_strain  = 7                                                  
         cohesive_elem     = cohesive_ele_types(elem_type)                      
c                                                                               
                                                                                
         num_vals = 22                                                          
         if( stress ) num_vals = 26                                             
c                                                                               
c                       skip cohesive elements for now. they do not             
c                       contribute to nodal results at structure level.         
c                                                                               
         if ( cohesive_elem ) cycle                                             
c                                                                               
         do i = 1, mxelmp                                                       
            elem_out_map(i) = outmap(elem_type,i)                               
         end do                                                                 
c                                                                               
c                       duplicate necessary element block data.                 
c                                                                               
         call oudups( span, felem, num_int_points, geo_non_flg, stress,         
     &                  .false. )                                               
c                                                                               
c                       compute the element block nodal stress/strain           
c                       data.                                                   
c                                                                               
         iblk   = blk                                                           
         ifelem = felem                                                         
         call ouprks( span, iblk, ifelem, elem_type, int_order,                 
     &                num_int_points, num_enodes, geo_non_flg, stress,          
     &                mat_type, center_output, num_short_stress,                
     &                num_short_strain, .false. )                               
c                                                                               
c                       add the element block nodal stress/strain data          
c                       into the global nodal results data structure.           
c                                                                               
         call oupads( span, incid(incmap(felem)), num_enodes, num_vals,         
     &                elem_out_map, nonode, nodal_values(1) )                   
c                                                                               
      end do                                                                    
c                       for scalar processing and/or mpi with just              
c                       one process, average the total results at               
c                       each node, using node count for the model.              
c                       then calculate the extended values of stress            
c                       and strain at each node in the model.                   
c                                                                               
c                       we can have nodes with only interface-cohesive          
c                       elements attached. make sure they have a vector         
c                       of zero values to write in the file.                    
c                                                                               
c                       for multi-process mpi jobs, we just write               
c                       the summed results + the node count to the              
c                       simplified result file (for the domain).                
c                                                                               
c                       the external program to combine mpi result              
c                       files handles nodes with missing values.                
c                                                                               
c                                                                               
      serial = .not. use_mpi                                                    
      do_average =  serial .or. (use_mpi .and. numprocs .eq. 1)                 
c                                                                               
      if( do_average ) then                                                     
       do snode = 1, nonode                                                     
         count = nodal_values(snode)%count                                      
         if( count .eq. 0 ) then  ! probably node w/ only cohesive              
            allocate( nodal_values(snode)%node_values(mxstmp) )                 
            nodal_values(snode)%node_values(1:mxstmp) = zero                    
            nodal_values(snode)%count = 1                                       
            cycle                                                               
         end if                                                                 
         rn_count = dble(count)                                                 
         snode_values => nodal_values(snode)%node_values                        
         do j = 1, num_vals                                                     
           map = elem_out_map(j)                                                
           snode_values(map) = snode_values(map) / rn_count                     
         end do                                                                 
         call ouext2( snode_values(1), 1, 1, stress )                           
       end do ! on snode                                                        
      end if                                                                    
c                                                                               
c                       output the averaged total nodal results                 
c                       to a (1) file compatable with patran for post           
c                       processing in binary or formatted file type             
c                       or (2) flat file in text or stream types.               
c                                                                               
c                       only results for nodes that appear in this              
c                       domain are written. see notes                           
c                       above for values actually written.                      
c                                                                               
      call oustpa( stress, oubin, ouasc, num_vals, nodal_values(1),             
     &             nonode, flat_file, stream_file, text_file,                   
     &             compressed )                                                 
c                                                                               
c                       deallocate all instances of the derived type            
c                       then the structure size vector                          
c                                                                               
      do node = 1, nonode                                                       
       if ( nodal_values(node)%count .ne. 0 ) then                              
        snode_values => nodal_values(node)%node_values                          
         deallocate(nodal_values(node)%node_values)                             
       end if                                                                   
      end do                                                                    
      deallocate( nodal_values )                                                
c                                                                               
      return                                                                    
      end                                                                       
