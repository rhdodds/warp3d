c     ****************************************************************
c     *                                                              *
c     *                      find_patran_work                        *
c     *                                                              *
c     ****************************************************************
c
      subroutine find_patran_work( combine_binary, combine_formatted, 
     &                             local_count, matl_states_name )  
c
      use global_data, only : termout, nodal_stresses_step_list,
     &   nodal_strains_step_list, element_stresses_step_list,
     &   element_strains_step_list, element_states_step_list,
     &   found_node_stresses, found_node_strains, found_ele_stresses,
     &   found_ele_strains, found_ele_states, last_step_no
      implicit none
c
c                parameter declarations
c                ----------------------
c
      integer :: local_count
      logical :: combine_binary, combine_formatted   
      character(len=*) :: matl_states_name
c
c                local variables
c                ---------------    
c
      integer :: step
      logical :: l1, l2, l3, l4, l5  
      logical, parameter :: local_debug = .false.    
      character(len=80) :: node_stress_fname, node_strain_fname,
     &                     element_stress_fname, element_strain_fname,
     &                     element_states_fname
      character(len=80) :: dummy_file_name
c
c                examine all files in directory to build lists
c                of step numbers to be combined for each type
c                of result file
c
      do step = 1, last_step_no
        nodal_stresses_step_list(step)   = 0
        nodal_strains_step_list(step)    = 0
        element_stresses_step_list(step) = 0
        element_strains_step_list(step)  = 0
        element_states_step_list(step)   = 0
      end do
c
      local_count = 0
      found_node_stresses = .false.
      found_node_strains  = .false.
      found_ele_stresses  = .false.
      found_ele_strains   = .false.
      found_ele_states    = .false.
c
      do step = 1, last_step_no
c
        call build_patran_file_names( step, 0, node_stress_fname,
     &                         node_strain_fname, element_stress_fname,
     &                         element_strain_fname, matl_states_name,
     &                         element_states_fname, termout,
     &                         combine_binary, combine_formatted,
     &                         dummy_file_name, 1, 1 )
        inquire( file=node_stress_fname,exist=l1 )
        inquire( file=node_strain_fname,exist=l2 )
        inquire( file=element_stress_fname,exist=l3 )
        inquire( file=element_strain_fname,exist=l4 )
        inquire( file=element_states_fname,exist=l5 )
c
        if( local_debug .and. step <=  2 ) then
            write(termout,9000) step, node_stress_fname,
     &                         node_strain_fname, element_stress_fname,
     &                         element_strain_fname, matl_states_name,
     &                         element_states_fname 
            write(termout,9010) l1, l2, l3, l4, l5
        end if     
c
        if( l1 ) then
            nodal_stresses_step_list(step) = 1
            found_node_stresses = .true.
        end if
        if( l2 ) then
            nodal_strains_step_list(step) = 1
            found_node_strains = .true.
        end if
        if( l3 ) then
            element_stresses_step_list(step) = 1
            found_ele_stresses = .true.
        end if
        if( l4 ) then
            element_strains_step_list(step)  = 1
            found_ele_strains = .true.
        end if
        if( l5 ) then
            element_states_step_list(step)  = 1
            found_ele_states = .true.
        end if
        if( l1 .or. l2 .or. l3 .or. l4 .or. l5 ) 
     &            local_count = local_count + 1
c
      end do
c
      return  
c
 9000 format("... Patran file names for step: ",i7,/,10(10x,a,/))  
 9010 format("    Exit flags for step: ",10l2)        
c
      end subroutine find_patran_work
c **********************************************************************
c *                                                                    *
c *     support routine:      build_patran_file_names                  *
c *                                                                    *
c **********************************************************************
c
c
      subroutine build_patran_file_names( step, myid, node_stress_fname,
     &                         node_strain_fname, element_stress_fname,
     &                         element_strain_fname, matl_states_name, 
     &                         element_states_fname, termout,
     &                         combine_binary, combine_formatted,
     &                         file_name, result_type, data_type  )
      implicit none
c
c                parameter declarations
c                ----------------------
c
      integer :: termout, step, myid, result_type, data_type
      logical ::  combine_binary, combine_formatted
      character(len=*) :: node_stress_fname, node_strain_fname,
     &                    element_stress_fname, element_strain_fname,
     &                    element_states_fname, file_name,
     &                    matl_states_name
c
c                local variables
c                ---------------    
c
      integer :: nchars
      logical :: local_debug
      character :: step_id*5, proc_id*4, format*1
c
      local_debug              = .false.
      node_stress_fname(1:)    = ' '  
      node_strain_fname(1:)    = ' '
      element_stress_fname(1:) = ' '
      element_strain_fname(1:) = ' '
      element_states_fname(1:) = ' '
c
      write(step_id,9000) step
      write(proc_id,9100) myid
c
      if( combine_binary ) then
         format = 'b'
      elseif( combine_formatted ) then
         format = 'f'
      else
         write(termout,*) 
     &       '>> FATAL ERROR: routine build_patran_file_names'
         write(termout,*) '>>   invalid format flag...'
         write(termout,*) '>>   job terminated...'
         stop
      end if
c
      node_stress_fname(1:)    = 'wn' // format // 's' // 
     &                           step_id // '.' // proc_id
      node_strain_fname(1:)    = 'wn' // format // 'e' //
     &                           step_id // '.' // proc_id
      element_stress_fname(1:) = 'we' // format // 's' //
     &                           step_id // '.' // proc_id
      element_strain_fname(1:) = 'we' // format // 'e' // 
     &                           step_id // '.' // proc_id
      nchars = len_trim( matl_states_name )
      element_states_fname(1:) = 'we' // format // 'm' // 
     &                           step_id // '_' //
     &                           matl_states_name(1:nchars) // '.' //
     &                           proc_id
c
c          data_type = 0   (element results)
c                    = 1   (nodal results)
c
c          result_type  = 0 (stress results)
c                       = 1 (strain results)  
c                       = 2 (states results)
c
c
      if( data_type .eq. 0 ) then
        if( result_type .eq. 1 ) then
           file_name(1:) = element_strain_fname(1:)
        elseif( result_type .eq. 0 ) then
           file_name(1:) = element_stress_fname(1:)
        elseif( result_type .eq. 2 ) then
           file_name(1:) = element_states_fname(1:)
        else
           write(termout,9300)
           stop
        end if
      elseif( data_type .eq. 1 ) then 
        if( result_type .eq. 1 ) then
           file_name(1:) = node_strain_fname(1:)
        elseif( result_type .eq. 0 ) then
           file_name(1:) = node_stress_fname(1:)
        else
           write(termout,9300)
           stop
        end if
      else
        write(termout,9400)
        stop
      end if
c       
      if( local_debug ) then
        write(termout,*) ' >> local debug in build_patran_file_names...'
        write(termout,9200) node_stress_fname, node_strain_fname,
     &                      element_stress_fname, element_strain_fname,
     &                      element_states_fname
      end if
c 
      return      
c
 9000 format( i5.5 )
 9100 format( i4.4 )
 9200 format( 5x,' > file names: ', 5(/10x,a) )
 9300 format( '>> FATAL ERROR: routine build_patran_file_names...',
     &  /,    '                invalid results_type',
     &  /,    '                job terminated...' )
 9400 format( '>> FATAL ERROR: routine build_patran_file_names...',
     &  /,    '                invalid data_type',
     &  /,    '                job terminated...' )
c
      end
c **********************************************************************
c *                                                                    *
c *          do_patran_nodal_step                                      *
c *                                                                    *
c **********************************************************************
c
c
      subroutine do_patran_nodal_step( step, combine_binary, 
     &          combine_formatted, termout, sig_eps_type, file_no,
     &          missing_flg, did_patran_binary, did_patran_formatted,
     &          consistent_num_values )
c
      use global_data, only : last_proc_id, num_model_nodes
      implicit none
c
c
c                parameter declarations
c                ----------------------
c
      integer :: step, termout, sig_eps_type, file_no,
     &           consistent_num_values
      logical :: combine_binary, combine_formatted, missing_flg,
     &           did_patran_binary, did_patran_formatted 
c
c                local variables
c                ---------------    
c
      integer:: myid, open_status, num_values, allocate_status,
     &          now_num_values, num_nodes, now_num_nodes, node_id, 
     &          miss_count
      integer, allocatable, dimension(:)  :: node_list
      logical :: local_debug, print_header, print, bad1, bad2
      double precision, allocatable,  dimension(:,:) :: node_values
      double precision, parameter :: zero = 0.0d00
      character(len=80) :: node_stress_fname, node_strain_fname,
     &                     element_stress_fname, element_strain_fname,
     &                     dummy_fname, dummy_name
      character(len=14) :: file_name
      character(len=4) :: title1_binary(80), title2_binary(80)
      character(len=1) :: title1_formatted(80), title2_formatted(80)
c
c
      local_debug = .false.
c
c                 open the root processor file for this load step.
c                 get information about step then close down this
c                 file.
c
      dummy_name(1:) = "..dummy.."
      call build_patran_file_names( step, 0, node_stress_fname,
     &   node_strain_fname, element_stress_fname,
     &   element_strain_fname, dummy_name, dummy_fname, termout,
     &   combine_binary, combine_formatted, file_name, sig_eps_type, 1 )
c
      if( local_debug ) write(termout,9700) step, file_name(1:14),
     &                   file_no 
c
c                open the node result file for processor
c                zero. read the header records to get the 
c                number of data values per node and the number
c                of nodes. each result file for other processors
c                must have the same values else bad.
c
      call open_old_patran_result_file( file_no, file_name, 
     &   combine_binary, combine_formatted, termout, 1, open_status, 0 )
c
      call read_patran_node_file_header( 
     &   file_no, title1_binary, num_values, title2_binary, num_nodes, 
     &   title1_formatted, title2_formatted, termout, combine_binary,
     &   combine_formatted )
      close( unit=file_no )
c  
c                check consistency
c    
      if(  local_debug ) write(termout,9710) num_values, num_nodes
      if( consistent_num_values == -1 )
     &          consistent_num_values = num_values
      bad1 = num_nodes .ne. num_model_nodes
      bad2 = num_values .ne. consistent_num_values
      if( bad1 .or. bad2 ) then
        write(termout,9730) num_nodes, num_model_nodes, num_values,
     &                      consistent_num_values,
     &                      file_name(1:len_trim(file_name))  
        stop
      end if
c
c                 allocate space to hold the data values for all
c                 nodes in the model. keep track of which nodes
c                 have been read from the files for the step.
c
      allocate( node_values(num_nodes,num_values),
     &          stat=allocate_status)
      call check_allocate_status( termout, 1, allocate_status )
c
      allocate( node_list(num_nodes), stat=allocate_status )
      call check_allocate_status( termout, 2, allocate_status )
c
      node_values(1:num_nodes,1:num_values) = zero
      node_list(1:num_nodes)                = 0
      if(  local_debug ) write(termout,9720) num_nodes
c
c                loop over all the separate rank files of 
c                node results for this load step.
c                open the file and read the
c                node results into the single array here. keep
c                track of how many times a node appears in result files 
c                to enable averaging.
c
      do myid = 0, last_proc_id ! MPI ranks
c
        call build_patran_file_names( step, myid, node_stress_fname,
     &     node_strain_fname, element_stress_fname, 
     &     element_strain_fname, dummy_name, dummy_fname, termout,
     &     combine_binary, combine_formatted, file_name, sig_eps_type, 
     &     1 )
c
c                open the node result file for the processor.
c                read all data for a binary or a formmated file. if
c                the open fails, keep looking for processor files
c                until last_proc_id in case user blew away a
c                processor file.
c
        call open_old_patran_result_file(
     &     file_no, file_name, combine_binary, combine_formatted,
     &     termout, 2, open_status, 1 )
        if( open_status .ne. 0 ) cycle
        if( local_debug ) write(termout,9705) myid, file_name(1:14),
     &                                        file_no
c
c               read node results for this processor. check
c               header info for consistency with other processor
c               files for this step.
c
        call read_patran_node_file_header_dummy( file_no, 
     &     now_num_values, now_num_nodes, termout, combine_binary,
     &     combine_formatted )
        bad1 = now_num_nodes .ne. num_model_nodes
        bad2 = now_num_values .ne. consistent_num_values
        if( bad1 .or. bad2 ) then
          write(termout,9730) now_num_nodes, num_model_nodes, 
     &                        now_num_values, consistent_num_values,
     &                        file_name(1:len_trim(file_name))  
          stop
        end if
c
c                read data values until end of file. stuff into the
c                single node results table. update count of
c                node appearances.
c
       call read_patran_node_values( file_no, node_values, num_values,
     &    num_nodes, node_list, termout, combine_binary, 
     &    combine_formatted ) 
     &   
       close( unit=file_no,status='keep' )
c
c              end loop over processor files for step
c
      end do
c
c                issue warning messages for nodes that have no results.
c                disable printing of messages in curent version based
c                on a global message in calling routine.
c
      print_header = .true.
      miss_count   = 0
      print        = .false.
c
      do node_id = 1, num_nodes
         if( node_list(node_id) .eq. 0 ) missing_flg = .true.
      end do
c
      if( print ) then
       do node_id = 1, num_nodes
         if( node_list(node_id) .eq. 0 ) then
           if( print_header ) write(termout,9300)
           print_header = .false.
           write(termout,9400) node_id
           miss_count = miss_count + 1
           if( miss_count .eq. 5 ) then
              write(termout,9301)
              exit
           end if
         end if ! on node_list
       end do
      end if ! on print
c
c               node results data looks ok. create a new,
c               single patran results file for this step.
c               write in all the data.
c
      call open_new_patran_result_file( file_no, file_name, 
     &   combine_binary, combine_formatted, termout, 2 )
c
      if( combine_binary ) did_patran_binary = .true.
      if( combine_formatted ) did_patran_formatted = .true.
      call write_patran_node_results( file_no, title1_binary, 
     &   num_values, title2_binary, num_nodes, node_values, node_list,
     &   title1_formatted, title2_formatted, termout, combine_binary, 
     &   combine_formatted, sig_eps_type )
c    
      close( unit=file_no )
c
      deallocate( node_values, node_list, stat=allocate_status )
      call check_deallocate_status( termout, 20, allocate_status )
c
      return    
c
 9120 format( '>> FATAL ERROR: routine do_patran_nodal_step...',
     &  /,    '                internal error @ 4',
     &  /,    '                job terminated...' )
 9300 format(10x,'> Note: The following nodes have no results',
     &     /,10x '        in the Patran result files...')
 9301 format(10x,'> Note: Further error messages of this nature',
     &     /,10x '        are suppressed...')
 9400 format(15x,i7)
c
 9700 format(5x,'> local debug in do_node_step: ',
     & /,10x,'step: ',i5,/,10x,'file_name: ',a14,
     & /,10x,'file_no: ',i5 )
 9705 format(5x,'> local debug in do_node_step: ',
     & /,10x,'myid: ',i5,/,10x,'file_name: ',a14,
     & /,10x,'file_no: ',i5 )
 9710 format(10x,'number of data values, nodes: ',i5,i9)
 9720 format(5x,'> data arrays allocated. no. nodes: ',i7)
 9730 format( '>> FATAL ERROR: routine do_patran_nodal_step...',
     &  /,    '                inconsistent values.',
     &  /,    '                num_nodes, num_model_nodes:',2i10,
     &  /,    '                num_values, num_consistent: ',2i10,
     &  /,    '                file: ',a,
     &  /,    '                job terminated...' )
c
      end

c **********************************************************************
c *                                                                    *
c *          do_patran_element_step                                    *
c *                                                                    *
c **********************************************************************
c
c
      subroutine do_patran_element_step( step, combine_binary,
     &  combine_formatted, termout, result_type, file_no, missing_flg,
     &  matl_state_name, did_patran_binary, did_patran_formatted,
     &  consistent_num_values )
c
      use global_data, only : num_model_elems, last_proc_id 

      implicit none
c
c                parameter declarations
c                ----------------------
c
      integer :: step, termout, result_type, file_no, 
     &           consistent_num_values
      logical :: combine_binary, combine_formatted, missing_flg,
     &           did_patran_binary, did_patran_formatted  
      character(len=*) :: matl_state_name
c
c                local variables
c                ---------------    
c
      integer :: myid, open_status, num_values, allocate_status,
     &           now_num_values, num_elements, element_id, 
     &           now_num_elements, miss_count
      character(len=80) :: node_stress_fname, node_strain_fname,
     &               element_stress_fname, element_strain_fname,
     &               element_states_fname
      character :: file_name*80, title1_binary(80)*4, 
     &             title2_binary(80)*4, title1_formatted(80)*1,
     &             title2_formatted(80)*1
      logical :: local_debug, fatal, print_header, print, bad1, bad2
c
      real, allocatable,   dimension(:,:) :: element_values
      integer, allocatable, dimension(:)  :: element_list
c
      local_debug = .false.
c
c                 open the root processor file for this load step.
c                 get information about step then close down this
c                 file.
c
      call build_patran_file_names( step, 0, node_stress_fname,
     &   node_strain_fname, element_stress_fname,
     &   element_strain_fname, matl_state_name, element_states_fname, 
     &   termout, combine_binary, combine_formatted, file_name,
     &   result_type, 0 )
      if( local_debug ) write(termout,9700) step,
     &                  file_name(1:len_trim(file_name)), file_no 
c
c                open the element result file for processor
c                zero. read the header records to get the 
c                number of data values per element and the number
c                of elements. each result file for other processors
c                must have the same values else bad.
c
      call open_old_patran_result_file( file_no, file_name,
     &   combine_binary, combine_formatted, termout, 1, open_status, 0 )
c
      call read_patran_element_file_header( file_no, title1_binary,
     &   num_values, title2_binary, num_elements, title1_formatted,
     &   title2_formatted, termout, combine_binary, combine_formatted )
      close( unit=file_no )
c      
      if( local_debug ) write(termout,9710) num_values, num_elements
      if( consistent_num_values == -1 )
     &          consistent_num_values = num_values
      bad1 = num_elements .ne. num_model_elems
      bad2 = num_values .ne. consistent_num_values
      if( bad1 .or. bad2 ) then
        write(termout,9730) num_elements, num_model_elems, num_values,
     &                      consistent_num_values,
     &                      file_name(1:len_trim(file_name))  
 
        stop
      end if
c
c                 allocate space to hold the data values for all
c                 elements in the mode. keep track of which element
c                 values have been read from the files for the step.
c
      allocate( element_values(num_values,num_elements),
     &          stat=allocate_status)
      call check_allocate_status( termout, 1, allocate_status )
      allocate( element_list(num_elements), stat=allocate_status )
      call check_allocate_status( termout, 2, allocate_status )
      element_list(1:num_elements) = 0
      element_values(1:num_values,1:num_elements) = 0.0
      if( local_debug ) write(termout,9720) num_elements
c
c                loop over all the separate processor files of 
c                element results for this load step.
c                open the file and read the
c                element results into the single array here. keep
c                track if an element appears other than once... bad.
c
      do myid = 0, last_proc_id
c
        call build_patran_file_names( step, myid, node_stress_fname,
     &     node_strain_fname, element_stress_fname,
     &     element_strain_fname, matl_state_name, 
     &     element_states_fname, termout, combine_binary, 
     &     combine_formatted, file_name, result_type, 0 )
        if( local_debug ) write(termout,9705) myid,
     &                     file_name(1:len_trim(file_name)), file_no 
c
c                open the element result file for the processor.
c                read all data for a binary or a formmated file. if
c                the open fails, keep looking for processor files
c                until max_procs in case user blew away a
c                processor file.
c
        call open_old_patran_result_file( file_no, file_name, 
     &       combine_binary,combine_formatted, termout, 2, 
     &       open_status, 1 )
        if( open_status .ne. 0 ) cycle
c
c               read element results for this processor. check
c               header info for consistency with other processor
c               files for this step.
c
        call read_patran_element_file_header_dummy( file_no, 
     &     now_num_values, termout, combine_binary, combine_formatted,
     &     now_num_elements )
        bad1 = now_num_elements .ne. num_model_elems
        bad2 = now_num_values .ne. consistent_num_values
        if( bad1 .or. bad2 ) then
          write(termout,9730) num_elements, num_model_elems, num_values,
     &                        consistent_num_values, 
     &                        file_name(1:len_trim(file_name))  
          stop
        end if

c
c                read data values until end of file. stuff into the
c                single element results table. update count of
c                element appearances.
c
       call read_patran_element_values( file_no, element_values,
     &     num_values, num_elements, element_list, termout,
     &     combine_binary, combine_formatted ) 
     &   
       close( unit=file_no )
c
c              end loop over processor files for step
      end do
c
c                check that each element has exactly 1 set of 
c                values read from all the processor file
c
      fatal        = .false.
      print_header = .true.
      miss_count   = 0
      print        = .false.
c
      do element_id = 1, num_elements
         if( element_list(element_id) .eq. 0 ) missing_flg = .true.
      end do

      if( print ) then
        do element_id = 1, num_elements
         if(  element_list(element_id) .ne. 1 ) then
           if(  print_header ) write(termout,9300)
           print_header = .false.
           write(termout,9400) element_id
           miss_count = miss_count + 1
           if(  miss_count .eq. 5 ) then
              write(termout,9301)
              exit
           end if

         end if
        end do
      end if
c
c               element results data look ok. create a new,
c               single patran results file for this step.
c               write in all the data.
c
      if( combine_binary )    did_patran_binary = .true.
      if( combine_formatted ) did_patran_formatted = .true.
      call open_new_patran_result_file( file_no, file_name, 
     &             combine_binary, combine_formatted, termout, 1 )
c
      call write_patran_element_results( 
     &         file_no, title1_binary, num_values, title2_binary,
     &         num_elements, element_values, title1_formatted,
     &         title2_formatted, termout, combine_binary, 
     &         combine_formatted, element_list )
c    
      close( unit=file_no )
c
      deallocate( element_values, element_list, stat=allocate_status)
      call check_deallocate_status( termout, 1, allocate_status )
c
      return    
c
 9120 format( '>> FATAL ERROR: routine do_patran_element_step...',
     &  /,    '                internal error @ 4',
     &  /,    '                job terminated...' )
 9300 format(15x,'> Note: The following elements have an appearance',
     &     /,15x '        count not equal to 1 in Patran result',
     &     /,15x '        files. This could result from any',
     &     /,15x,'        elements which do not emit results.')
 9301 format(15x,'> Note: Further messages of this nature',
     &     /,15x '        are suppressed...')
 9400 format(15x,i7)
c
 9700 format(5x,'> local debug in do_patran_element_step: ',
     & /,10x,'step: ',i5,/,10x,'file_name: ',a,
     & /,10x,'file_no: ',i5 )
 9705 format(5x,'> local debug in do_patran_element_step: ',
     & /,10x,'myid: ',i5,/,10x,'file_name: ',a,
     & /,10x,'file_no: ',i5 )
 9710 format(10x,'number of data values, elements: ',i5,i10)
 9720 format(5x,'> data arrays allocated. no. elements: ',i7)
 9730 format( '>> FATAL ERROR: routine do_patran_element_step...',
     &  /,    '                inconsistent values.',
     &  /,    '                num_elements, num_model_elements:',2i10,
     &  /,    '                num_values, num_consistent: ',2i10,
     &  /,    '                file: ',a,
     &  /,    '                job terminated...' )
c
      end
c **********************************************************************
c *                                                                    *
c *          open_old_patran_result_file                               *
c *                                                                    *
c **********************************************************************
c
c
      subroutine open_old_patran_result_file( file_no, file_name, 
     &             combine_binary, combine_formatted, termout, 
     &             location, open_status, dowhat )
      implicit none
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, termout, location, open_status, dowhat
      character(len=*) :: file_name
      logical :: combine_binary, combine_formatted
c
c                  local declarations
c                  ------------------
c
c
      if( combine_binary ) then
        open( unit=file_no, file=file_name, status='old', 
     &        access='sequential', form='unformatted',
     &        iostat= open_status )
      elseif( combine_formatted ) then 
        open( unit=file_no, file=file_name, status='old', 
     &        access='sequential', form='formatted',
     &        iostat= open_status )
      else 
        write(termout,9200) location
        stop
      end if
c
c                 dowhat = 0 : check open status and stop job if
c                              open failed
c                        = 1 : don't check open status
c
      if( dowhat .eq. 1 ) return
c
      if( open_status .ne. 0 ) then
        write(termout,9100) location
        stop
      end if        
c
      return
c
 9100 format( '>> FATAL ERROR: failed file open @ ',i2,
     &  /,    '                job terminated...' )
 9200 format( '>> FATAL ERROR: invalid system state @ ',i2,
     &  /,    '                job terminated...' )
c
      end
c **********************************************************************
c *                                                                    *
c *          read_patran_element_file_header                           *
c *                                                                    *
c **********************************************************************
c
c
      subroutine read_patran_element_file_header( file_no, 
     &                   title1_binary,
     &                   num_values, title2_binary, num_elements, 
     &                   title1_formatted, title2_formatted, termout,
     &                   combine_binary, combine_formatted )
      implicit none
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, termout, num_values, num_elements
      character(len=4) :: title1_binary(80), title2_binary(80)
      character(len=1) :: title1_formatted(80), title2_formatted(80)
      logical :: combine_binary, combine_formatted
c
c                  local declarations
c                  ------------------
c
      integer :: read_status
c              
      if( combine_binary ) then 
         read(unit=file_no,iostat=read_status) title1_binary, num_values
         call check_read_status( file_no, read_status, termout, 1 )
         read(unit=file_no,iostat=read_status) title2_binary
         call check_read_status( file_no, read_status, termout, 2 )
         read(unit=file_no,iostat=read_status) title2_binary
         call check_read_status( file_no, read_status, termout, 3 )
         read(unit=file_no,iostat=read_status) num_elements
         call check_read_status( file_no, read_status, termout, 3 )
      elseif( combine_formatted ) then  
         read(unit=file_no,iostat=read_status,
     &        fmt=9300) title1_formatted
         call check_read_status( file_no, read_status, termout, 4 )
         read(unit=file_no,iostat=read_status,fmt=*)
     &        num_values
         call check_read_status( file_no, read_status, termout, 5 )
         read(unit=file_no,iostat=read_status,
     &        fmt=9300) title2_formatted
         call check_read_status( file_no, read_status, termout, 6 )
         read(unit=file_no,iostat=read_status,
     &        fmt=9300) title2_formatted
         call check_read_status( file_no, read_status, termout, 7 )
         read(unit=file_no,iostat=read_status,fmt=*) num_elements
         call check_read_status( file_no, read_status, termout, 8 )
      else
        write(termout,9200)
        stop
      end if
c
      return
c
 9200 format( '>> FATAL ERROR: invalid system state in ',
     &  /,    '                routine: read_element_file_header',
     &  /,    '                job terminated...' )
 9300 format(80a1)
 9320 format(i10)
c
      end
c **********************************************************************
c *                                                                    *
c *          read_patran_element_file_header_dummy                     *
c *                                                                    *
c **********************************************************************
c
c
      subroutine read_patran_element_file_header_dummy(
     &          file_no, now_num_values, termout,
     &          combine_binary, combine_formatted, now_num_elements )
      implicit none
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, termout, now_num_values, now_num_elements
      logical :: combine_binary, combine_formatted
c
c                  local declarations
c                  ------------------
c
      integer :: read_status
      character(len=4) :: title1_binary(80), title2_binary(80)
      character(len=1) :: title1_formatted(80), title2_formatted(80)
c
      if(  combine_binary ) then 
         read(unit=file_no,iostat=read_status) title1_binary, 
     &        now_num_values
         call check_read_status( file_no, read_status, termout, 1 )
         read(unit=file_no,iostat=read_status) title2_binary
         call check_read_status( file_no, read_status, termout, 2 )
         read(unit=file_no,iostat=read_status) title2_binary
         call check_read_status( file_no, read_status, termout, 3 )
         read(unit=file_no,iostat=read_status) now_num_elements
         call check_read_status( file_no, read_status, termout, 3 )
      elseif( combine_formatted ) then  
         read(unit=file_no,iostat=read_status,
     &        fmt=9300) title1_formatted
         call check_read_status( file_no, read_status, termout, 4 )
         read(unit=file_no,iostat=read_status,fmt=*)
     &         now_num_values
         call check_read_status( file_no, read_status, termout, 5 )
         read(unit=file_no,iostat=read_status,
     &        fmt=9300) title2_formatted
         call check_read_status( file_no, read_status, termout, 6 )
         read(unit=file_no,iostat=read_status,
     &        fmt=9300) title2_formatted
         call check_read_status( file_no, read_status, termout, 7 )
         read(unit=file_no,iostat=read_status,fmt=*) now_num_elements
         call check_read_status( file_no, read_status, termout, 8 )
      else
        write(termout,9200)
        stop
      end if
c
      return
c
 9200 format( '>> FATAL ERROR: invalid system state in ',
     &  /,    '                routine: read_element_file_header_dummy',
     &  /,    '                job terminated...' )
 9300 format(80a1)
 9320 format(i10)
c
      end


c **********************************************************************
c *                                                                    *
c *          read_patran_element_values                                *
c *                                                                    *
c **********************************************************************
c
c
      subroutine read_patran_element_values(
     &           file_no, element_values, num_values, num_elements,
     &           element_list, termout, combine_binary, 
     &           combine_formatted  ) 
      implicit none
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, num_values, num_elements, element_list(*),
     &           termout
      logical :: combine_binary, combine_formatted
      real :: element_values(num_values,num_elements)
c
c                  local declarations
c                  ------------------
c
      integer :: read_status, element_id, elem_family, i
      real :: element_data(100)
c
c
      do
        if( combine_binary ) then 
          read(unit=file_no,iostat=read_status) element_id,
     &         elem_family, 
     &        ( element_data(i), i = 1, num_values )
        elseif( combine_formatted ) then
          read(unit=file_no,iostat=read_status,fmt=9330)
     &         element_id, elem_family,
     &        ( element_data(i), i = 1, num_values )
        else
         write(termout,9200)
         stop
        end if
c
        if( read_status .ne. 0 ) return
        if( element_id .le. num_elements ) then
          element_values(1:num_values,element_id) =
     &                    element_data(1:num_values)
          element_list(element_id) = element_list(element_id) + 1
        else
          write(termout,9130) num_values, num_elements, element_id
          stop
        end if
c
      end do
c
      return
c
 9130 format( '>> FATAL ERROR: routine read_patran_element_data...',
     &  /,    '                internal error @ 5',
     &  /,    '                num_values, num_elements, element_id: ',
     & 3i10,
     &  /,    '                job terminated...' )
 9330 format(2i8,/,(6e13.6))
 9200 format( '>> FATAL ERROR: invalid system state in ',
     &  /,    '                routine: read_element_data',
     &  /,    '                job terminated...' )
c
      end
c **********************************************************************
c *                                                                    *
c *          open_new_patran_result_file                               *
c *                                                                    *
c **********************************************************************
c
      subroutine open_new_patran_result_file( file_no, file_name, 
     &             combine_binary, combine_formatted, termout, 
     &             location )
      implicit none
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, termout, location
      character(len=*) :: file_name
      logical :: combine_binary, combine_formatted
c
c                  local declarations
c                  ------------------
c
      integer :: open_status, dot_pos
      character(len=80) :: local_fname
c
c                  remove the .proc_id from file name
c
      dot_pos = index( file_name, ".")
      if( dot_pos .eq. 0 ) then
        write(termout,9300) file_name(1:len_trim(file_name))
        stop
      end if
      local_fname(1:) = file_name(1:dot_pos-1)
c      
      if( combine_binary ) then
        open( unit=file_no, file=local_fname, status='unknown', 
     &        access='sequential', form='unformatted',
     &        iostat= open_status )
      elseif( combine_formatted ) then 
        open( unit=file_no, file=local_fname, status='unknown', 
     &        access='sequential', form='formatted',
     &        iostat= open_status )
      else 
        write(termout,9200) location
        stop
      end if
c
      if( open_status .ne. 0 ) then
        write(termout,9100) location
        stop
      end if        
c
      return
c
 9100 format( '>> FATAL ERROR: failed  open for new file @ ',i2,
     &  /,    '                open_new_patran_result_file',
     &  /,    '                job terminated...' )
 9200 format( '>> FATAL ERROR: invalid system state @ ',i2,
     &  /,    '                open_new_patran_result_file',
     &  /,    '                job terminated...' )
 9300 format( '>> FATAL ERROR: invalid system state in routine:',
     &  /,    '                open_new_patran_result_file',
     &  /,    '                bad file name: ',a,
     &  /,    '                job terminated...' )
c
      end
c **********************************************************************
c *                                                                    *
c *          write_patran_element_results                              *
c *                                                                    *
c **********************************************************************
c
c
      subroutine write_patran_element_results( 
     &         file_no, title1_binary, num_values, title2_binary,
     &         num_elements, element_values, title1_formatted,
     &         title2_formatted, termout, combine_binary, 
     &         combine_formatted, element_list ) 
      implicit none 
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, termout, num_values, num_elements
      character(len=4) :: title1_binary(80), title2_binary(80)
      character(len=1) :: title1_formatted(80), title2_formatted(80)
      logical :: combine_binary, combine_formatted
      real :: element_values(num_values,num_elements)
      integer :: element_list(num_elements)
c
c                  local declarations
c                  ------------------
c
      integer ::  element_id, write_status, i


c       
      if( combine_binary ) then 
c
         write(unit=file_no,iostat=write_status) title1_binary,
     &                                           num_values
         call check_write_status( file_no, write_status, termout, 1 )
         write(unit=file_no,iostat=write_status) title2_binary
         call check_write_status( file_no, write_status, termout, 2 )
         write(unit=file_no,iostat=write_status) title2_binary
         call check_write_status( file_no, write_status, termout, 3 )
         do element_id = 1, num_elements
          if( element_list(element_id) .eq. 1 )
     &        write(unit=file_no,iostat=write_status)
     &        element_id, 8,
     &        (element_values(i,element_id), i = 1, num_values)
          call check_write_status( file_no, write_status, termout, 10 )
         end do
c
      elseif( combine_formatted ) then
c
         write(unit=file_no,iostat=write_status,
     &        fmt=9300) title1_formatted
         call check_write_status( file_no, write_status, termout, 4 )
          write(unit=file_no,iostat=write_status,
     &        fmt=9315) num_values
         call check_write_status( file_no, write_status, termout, 5 )
         write(unit=file_no,iostat=write_status,
     &        fmt=9300) title2_formatted
         call check_write_status( file_no, write_status, termout, 6 )
         write(unit=file_no,iostat=write_status,
     &        fmt=9300) title2_formatted
         call check_write_status( file_no, write_status, termout, 7 )
         do element_id = 1, num_elements
          if(  element_list(element_id) .eq. 1 )
     &       write(unit=file_no,iostat=write_status,fmt=9330)
     &       element_id, 8,
     &       ( element_values(i,element_id), i = 1, num_values )
          call check_write_status( file_no, write_status, termout, 8 )
         end do
      else
        write(termout,9200)
        stop
c           
      end if
c
      return
c
 9200 format( '>> FATAL ERROR: invalid system state in routine:',
     &  /,    '                write_element_results',
     &  /,    '                job terminated...' )
 9300 format(80a1)
 9315 format(i5)
 9330 format(2i8,/,(6e13.6))
c
      end

c **********************************************************************
c *                                                                    *
c *          read_patran_node_file_header                              *
c *                                                                    *
c **********************************************************************
c
c
      subroutine read_patran_node_file_header( file_no, title1_binary,
     &                   num_values, title2_binary, num_nodes, 
     &                   title1_formatted, title2_formatted, termout,
     &                   combine_binary, combine_formatted )
      implicit none
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, termout, num_values, num_nodes
      character(len=4) :: title1_binary(80), title2_binary(80)
      character(len=1) :: title1_formatted(80), title2_formatted(80)
      logical :: combine_binary, combine_formatted
c
c                  local declarations
c                  ------------------
c
      integer :: read_status,  dum1, dum3
      real :: dum2
c
c              
      if( combine_binary ) then 
         read(unit=file_no,iostat=read_status) title1_binary,
     &                 num_nodes, dum1, dum2, dum3, num_values 
         call check_read_status( file_no, read_status, termout, 1 )
         read(unit=file_no,iostat=read_status) title2_binary
         call check_read_status( file_no, read_status, termout, 2 )
         read(unit=file_no,iostat=read_status) title2_binary
         call check_read_status( file_no, read_status, termout, 3 )
      elseif( combine_formatted ) then  
         read(unit=file_no,iostat=read_status,
     &        fmt=9300) title1_formatted
         call check_read_status( file_no, read_status, termout, 4 )
         read(unit=file_no,iostat=read_status,fmt=*)
     &         num_nodes, dum1, dum2, dum3, num_values
         call check_read_status( file_no, read_status, termout, 5 )
         read(unit=file_no,iostat=read_status,
     &        fmt=9300) title2_formatted
         call check_read_status( file_no, read_status, termout, 6 )
         read(unit=file_no,iostat=read_status,
     &        fmt=9300) title2_formatted
         call check_read_status( file_no, read_status, termout, 7 )
      else
        write(termout,9200)
        stop
      end if
c
      return
c
 9200 format( '>> FATAL ERROR: invalid system state in ',
     &  /,    '                routine: read_node_file_header',
     &  /,    '                job terminated...' )
 9300 format(80a1)
c
      end
c **********************************************************************
c *                                                                    *
c *          read_node_file_header_dummy                               *
c *                                                                    *
c **********************************************************************
c
c
      subroutine read_patran_node_file_header_dummy( file_no, 
     &  num_values, num_nodes, termout, combine_binary, 
     &  combine_formatted )
      implicit none
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, termout, num_values, num_nodes
      character(len=4) :: title1_binary(80), title2_binary(80)
      character(len=1) :: title1_formatted(80), title2_formatted(80)
      logical :: combine_binary, combine_formatted
c
c                  local declarations
c                  ------------------
c
      integer :: read_status, dum1, dum3
      real :: dum2
c
c              
      if( combine_binary ) then 
         read(unit=file_no,iostat=read_status) title1_binary,
     &                 num_nodes, dum1, dum2, dum3, num_values 
         call check_read_status( file_no, read_status, termout, 1 )
         read(unit=file_no,iostat=read_status) title2_binary
         call check_read_status( file_no, read_status, termout, 2 )
         read(unit=file_no,iostat=read_status) title2_binary
         call check_read_status( file_no, read_status, termout, 3 )
      elseif( combine_formatted ) then  
         read(unit=file_no,iostat=read_status,
     &        fmt=9300) title1_formatted
         call check_read_status( file_no, read_status, termout, 4 )
         read(unit=file_no,iostat=read_status,fmt=*)
     &         num_nodes, dum1, dum2, dum3, num_values
         call check_read_status( file_no, read_status, termout, 5 )
         read(unit=file_no,iostat=read_status,
     &        fmt=9300) title2_formatted
         call check_read_status( file_no, read_status, termout, 6 )
         read(unit=file_no,iostat=read_status,
     &        fmt=9300) title2_formatted
         call check_read_status( file_no, read_status, termout, 7 )
      else
        write(termout,9200)
        stop
      end if
c
      return
c
 9200 format( '>> FATAL ERROR: invalid system state in ',
     &  /,    '                routine: read_node_file_header_dummy',
     &  /,    '                job terminated...' )
 9300 format(80a1)
c
      end
c **********************************************************************
c *                                                                    *
c *          read_patran_node_values                                   *
c *                                                                    *
c **********************************************************************
c
c
      subroutine read_patran_node_values(
     &           file_no, node_values, num_values, num_nodes,
     &           node_list, termout, combine_binary, 
     &           combine_formatted  ) 
      implicit none
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, num_values, num_nodes, node_list(*),
     &           termout
      logical :: combine_binary, combine_formatted
      double precision :: node_values(num_nodes,num_values)
c
c                  local declarations
c                  ------------------
c
      integer :: read_status, node_id, lcount
      double precision :: node_data(100)
      integer :: bcount
c
c
      bcount = 0
      do
        if( combine_binary ) then 
          read(unit=file_no,iostat=read_status) node_id, lcount,
     &        node_data(1:num_values)
        elseif( combine_formatted ) then
          read(unit=file_no,iostat=read_status,fmt=9330) node_id,
     &                             lcount
          read(unit=file_no,iostat=read_status,fmt=9332)
     &          node_data(1:num_values)
          bcount = bcount + num_values + 1
        else
         write(termout,9200)
         stop
        end if
c
        if( read_status .ne. 0 ) return
        if( node_id .le. num_nodes ) then
          node_values(node_id,1:num_values) =
     &                 node_values(node_id,1:num_values) +
     &                 node_data(1:num_values)
          node_list(node_id) = node_list(node_id) + lcount
        else
          write(termout,9130)
          stop
        end if
c
      end do
c
      return
c
 9130 format( '>> FATAL ERROR: routine read_node_data...',
     &  /,    '                internal error @ 5',
     &  /,    '                job terminated...' )
 9330 format(2i8)
 9332 format(e23.15)
 9200 format( '>> FATAL ERROR: invalid system state in ',
     &  /,    '                routine: read_node_data',
     &  /,    '                job terminated...' )
c
      end

c **********************************************************************
c *                                                                    *
c *          write_patran_node_results                                 *
c *                                                                    *
c **********************************************************************
c
c
      subroutine write_patran_node_results( 
     &         file_no, title1_binary, num_values, title2_binary,
     &         num_nodes, node_values, node_list, title1_formatted,
     &         title2_formatted, termout, combine_binary, 
     &         combine_formatted, sig_eps_type ) 
      implicit none 
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, termout, num_values, num_nodes, node_list(*),
     &           sig_eps_type
      character(len=4) :: title1_binary(80), title2_binary(80)
      character(len=1) :: title1_formatted(80), title2_formatted(80)
      logical :: combine_binary, combine_formatted
      double precision :: node_values(num_nodes,num_values)
c
c                  local declarations
c                  ------------------
c
      integer :: node_id, write_status, dum1, dum3
      real :: dum2
      double precision :: scale
      real :: single_values(100)
c
      dum1 = 0
      dum2 = 0.0
      dum3 = 0
c
      do node_id = 1, num_nodes
        scale = 0.0d0
        if(  node_list(node_id) .gt. 0 ) 
     &       scale = 1.0d0 /  real( node_list(node_id) )
        node_values(node_id,1:num_values) =
     &            node_values(node_id,1:num_values) * scale
      end do      
c
      call compute_extra_sig_eps_values( node_values, num_nodes,
     &                                   sig_eps_type )
c       
      if( combine_binary ) then 
c
         write(unit=file_no,iostat=write_status) title1_binary,
     &                 num_nodes, dum1, dum2, dum3, num_values 
         call check_write_status( file_no, write_status, termout, 1 )
         write(unit=file_no,iostat=write_status) title2_binary
         call check_write_status( file_no, write_status, termout, 2 )
         write(unit=file_no,iostat=write_status) title2_binary
         call check_write_status( file_no, write_status, termout, 3 )
         do node_id = 1, num_nodes
          single_values(1:num_values) = node_values(node_id,
     &                                  1:num_values)
          if( node_list(node_id) .gt. 0 )
     &       write(unit=file_no,iostat=write_status)
     &       node_id, single_values(1:num_values)
          call check_write_status( file_no, write_status, termout, 12 )
         end do
c
      elseif( combine_formatted ) then  
c
         write(unit=file_no,iostat=write_status,
     &        fmt=9300) title1_formatted
         call check_write_status( file_no, write_status, termout, 4 )
         write(unit=file_no,iostat=write_status,
     &        fmt=9315) num_nodes, dum1, dum2, dum3, num_values
         call check_write_status( file_no, write_status, termout, 5 )
         write(unit=file_no,iostat=write_status,
     &        fmt=9300) title2_formatted
         call check_write_status( file_no, write_status, termout, 6 )
         write(unit=file_no,iostat=write_status,
     &        fmt=9300) title2_formatted
         call check_write_status( file_no, write_status, termout, 7 )
         do node_id = 1, num_nodes
          if( node_list(node_id) .gt. 0 )
     &        write(unit=file_no,iostat=write_status,fmt=9330)
     &        node_id, node_values(node_id,1:num_values)
          call check_write_status( file_no, write_status, termout, 8 )
         end do
c
      else
c
        write(termout,9200)
        stop
c
      end if
c
      return
c
 9200 format( '>> FATAL ERROR: invalid system state in routine:',
     &  /,    '                write_node_results',
     &  /,    '                job terminated...' )
 9300 format(80a1)
 9315 format(2i9,e15.6,2i9)
 9330 format(i8,(5e13.6))
c
      end
