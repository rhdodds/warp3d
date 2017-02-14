c     ****************************************************************
c     *                                                              *
c     *                      find_flat_work                          *
c     *                                                              *
c     ****************************************************************
c
      subroutine find_flat_work( combine_stream, combine_text, 
     &                           local_count, matl_states_name )
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
      logical :: combine_stream, combine_text 
      character(len=*) :: matl_states_name
c
c                local variables
c                ---------------    
c
      integer :: step
      logical :: l1, l2, l3, l4, l5      
      logical, parameter :: local_debug = .false.,
     &                      local_debug2 = .false. 
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
        call build_flat_file_names( step, 0, node_stress_fname,
     &                         node_strain_fname, element_stress_fname,
     &                         element_strain_fname, matl_states_name,
     &                         element_states_fname, termout,
     &                         combine_stream, combine_text,
     &                         dummy_file_name, 1, 1 )
        inquire( file=node_stress_fname,exist=l1 )
        inquire( file=node_strain_fname,exist=l2 )
        inquire( file=element_stress_fname,exist=l3 )
        inquire( file=element_strain_fname,exist=l4 )
        inquire( file=element_states_fname,exist=l5 )
c
        if( local_debug2 .and. step <=  30 ) then
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
        if( l1 .or. l2 .or. l3 .or. l4 .or. l5 ) then
                 local_count = local_count + 1
        end if
c
      end do
c
      if( .not. local_debug ) return
c
      write(termout,*) "... list of flat files to process ..." 
      write(termout,*) '      total number of files: ', local_count
      do step = 1, last_step_no
        if( nodal_stresses_step_list(step) == 1 )
     &      write(termout,*) '        step: ',step,' nodal stresses'
        if( nodal_strains_step_list(step) == 1 )
     &      write(termout,*) '        step: ',' nodal strains'
        if( element_stresses_step_list(step) == 1 )
     &      write(termout,*) '        step: ',step,' element stresses'
       if( element_strains_step_list(step) == 1 )
     &      write(termout,*) '        step: ',step,' element strains'
       if( element_states_step_list(step) == 1 )
     &      write(termout,*) '        step: ',step,' element states'
      end do 
c
      return      
c
 9000 format("... Flat file names for step: ",i7,/,10(10x,a,/))  
 9010 format("    Exit flags for step: ",10l2)        
c
      end subroutine find_flat_work
c **********************************************************************
c *                                                                    *
c *     support routine:      build_flat_file_names                    *
c *                                                                    *
c **********************************************************************
c
c
      subroutine build_flat_file_names( step, myid, node_stress_fname,
     &                         node_strain_fname, element_stress_fname,
     &                         element_strain_fname, matl_states_name, 
     &                         element_states_fname, termout,
     &                         combine_stream, combine_text,
     &                         file_name, result_type, data_type  )
      implicit none
c
c
c                parameter declarations
c                ----------------------
c
      integer :: termout, step, myid, result_type, data_type
      logical :: combine_stream, combine_text
      character(len=*) :: node_stress_fname, node_strain_fname,
     &                    element_stress_fname, element_strain_fname,
     &                    element_states_fname, file_name,
     &                    matl_states_name
c
c                local variables
c                ---------------    
c
      integer :: fcount, nchars
      character(len=5) :: step_id
      character(len=4) :: proc_id
      character(len=6) :: format
      logical :: local_debug
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
      if( combine_stream ) then
         format = 'stream'
         fcount = 6
      elseif( combine_text ) then
         format = 'text'
         fcount = 4
      else
         write(termout,*)
     &       '>> FATAL ERROR: routine build_flat_file_names'
         write(termout,*) '>>   invalid format flag...'
         write(termout,*) '>>   job terminated...'
         stop
      end if
c
      node_stress_fname(1:)    = 'wn' // 's' // 
     &                           step_id // '_' //
     &                           format(1:fcount) // "." // proc_id     
      node_strain_fname(1:)    = 'wn' // 'e' // 
     &                           step_id // '_' //
     &                           format(1:fcount) // "." // proc_id     
c
      element_stress_fname(1:) = 'we' // 's' //
     &                           step_id // '_' // 
     &                           format(1:fcount) // "." // proc_id
      element_strain_fname(1:) = 'we' // 'e' //
     &                           step_id // '_' // 
     &                           format(1:fcount) // "." // proc_id
      nchars = len_trim( matl_states_name )
      element_states_fname(1:) = 'we' // 'm' // 
     &                           step_id // '_' // format(1:fcount) //
     &                           '_' // matl_states_name(1:nchars) //
     &                           '.' // proc_id
c
c          data_type = 0   (element results)
c                    = 1   (nodal results)
c
c          result_type  = 0 (stress results)
c                       = 1 (strain results)  
c                       = 2 (states results)
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
        write(termout,*) ' >> local debug in build_flat_file_names...'
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
 9300 format( '>> FATAL ERROR: routine build_flat_file_names...',
     &  /,    '                invalid results_type',
     &  /,    '                job terminated...' )
 9400 format( '>> FATAL ERROR: routine build_file_names...',
     &  /,    '                invalid data_type',
     &  /,    '                job terminated...' )
c
      end
c **********************************************************************
c *                                                                    *
c *          do_flat_element_step                                      *
c *                                                                    *
c **********************************************************************
c
c
      subroutine do_flat_element_step( step, combine_stream,
     &   combine_text, termout,  results_type, file_no, missing_flg,
     &   matl_state_name, did_flat_stream, did_flat_text,
     &   consistent_num_values )
c
      use global_data, only :  last_proc_id, num_model_elems
      implicit none
c
c
c                parameter declarations
c                ----------------------
c
      integer :: step, termout, results_type, file_no,
     &           consistent_num_values
      logical :: combine_stream, combine_text, missing_flg,  
     &           did_flat_stream, did_flat_text
      character(len=*) :: matl_state_name
c
c                local variables
c                ---------------    
c
      integer :: proc_id, open_status, num_values,
     &           allocate_status, num_elements,
     &           element_id,miss_count, num_header_lines
      integer, parameter :: max_header_lines = 50 
      character(len=80) :: node_stress_fname, node_strain_fname,
     &                     element_stress_fname, element_strain_fname,
     &                     element_states_fname,
     &                     header_lines(max_header_lines)
      character(len=80) :: file_name
      logical :: local_debug, fatal, print_header, print
c
      double precision, allocatable,   dimension(:,:) :: element_values
      double precision, parameter :: zero = 0.0d00
      integer, allocatable, dimension(:)  :: element_list
c
      local_debug = .false.
c
c                 open the root processor file for this load step.
c                 get information about step then close down this
c                 file.
c
      file_name(1:) = " "
      call build_flat_file_names( step, 0, node_stress_fname,
     &                       node_strain_fname, element_stress_fname,
     &                       element_strain_fname, matl_state_name,
     &                       element_states_fname, termout,
     &                       combine_stream, combine_text,
     &                       file_name, results_type, 0 )
      if( local_debug ) write(termout,9700) step,
     &                  file_name(1:len_trim(file_name)), file_no 
c
c                open the element result file for processor
c                zero. read the header records to get the 
c                number of data values per element and the number
c                of elements. each result file for other processors
c                must have the same values else bad.
c
      call open_old_flat_result_file( file_no, file_name,
     &      combine_stream, combine_text, termout, 1, open_status, 0 )
c
      num_header_lines = 0
      if( combine_text ) call read_flat_header( file_no, header_lines,
     &               num_header_lines, max_header_lines, termout )
c
      if( local_debug ) 
     &     write(termout,*) '.... num_header_lines: ', num_header_lines
c
c                 allocate space to hold the data values for all
c                 elements in the mode. keep track of which element
c                 values have been read from the files for the step.
c
c                 results_type = 0,1,2 (stresses, strains, states)
c
      if( combine_stream ) read(file_no)   num_values
      if( combine_text )   read(file_no,*) num_values
      close( unit=file_no )
      if( consistent_num_values == -1 ) 
     &    consistent_num_values = num_values
      if( num_values .ne. consistent_num_values ) then
        write(termout,9730) num_values, consistent_num_values,
     &                      file_name(1:len_trim(file_name))  
        stop
      end if
      num_elements = num_model_elems
      if( local_debug ) write(termout,9710) num_values 
c      
      allocate( element_values(num_values,num_elements),
     &          stat=allocate_status)
      call check_allocate_status( termout, 1, allocate_status )
      allocate( element_list(num_elements), stat=allocate_status )
      call check_allocate_status( termout, 2, allocate_status )
      element_list(1:num_elements) = 0
      element_values(1:num_values,1:num_elements) = zero
      if( local_debug ) write(termout,9720) num_elements
c
c                loop over all the separate processor files of 
c                element results for this load step.
c                open the file and read the
c                element results into the single array here. keep
c                track if an element appears other than once... bad.
c
      do proc_id = 0, last_proc_id
c
        call build_flat_file_names( step, proc_id, node_stress_fname,
     &                       node_strain_fname, element_stress_fname,
     &                       element_strain_fname, matl_state_name,
     &                       element_states_fname, termout,
     &                       combine_stream, combine_text,
     &                       file_name, results_type, 0 )
        if( local_debug ) write(termout,9705) proc_id,
     &                    file_name(1:len_trim(file_name)), file_no  
c
c                open the element result file for the processor.
c                read all data for a binary or a formmated file. if
c                the open fails, keep looking for processor files
c                until max_procs in case user blew away a
c                processor file.
c
        call open_old_flat_result_file( file_no, file_name, 
     &       combine_stream, combine_text, termout, 2, open_status, 1 )
        if( open_status .ne. 0 ) cycle
c
c               skip element results header for this processor.
c
        if( combine_text ) call skip_comment_lines( file_no )
c
c                read data values until end of file. stuff into the
c                single element results table. update count of
c                element appearances.
c
       call read_flat_element_values( file_no, element_values, 
     &      num_values, num_elements, element_list, termout, 
     &      combine_stream, combine_text ) 
     &   
       close( unit=file_no )
c
c              end loop over processor files for step
      end do
c
c                check that each element has exactly 1 set of 
c                values read from all the processor file
c
      fatal = .false.
      print_header = .true.
      miss_count = 0
      print = .true.
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
c               element results data looks ok. create a new,
c               single flat results file for this step.
c               write in all the data.
c
c
      if( combine_text )   did_flat_text   = .true.
      if( combine_stream ) did_flat_stream = .true.
c
      call open_new_flat_result_file( file_no, file_name, 
     &     combine_stream, combine_text, termout, 1 )
c
      call write_flat_element_results( file_no, num_header_lines, 
     &    header_lines, num_values, num_elements, element_values,
     &    termout, combine_stream, combine_text, element_list )
c    
      close(unit=file_no,status='keep')
c
      deallocate( element_values, element_list, stat=allocate_status)
      call check_deallocate_status( termout, 1, allocate_status )
c
      return    
c
 9300 format(15x,'> Note: The following elements have an appearance',
     &     /,15x '        count not equal to 1 in Patran result',
     &     /,15x '        files. This could result from any',
     &     /,15x,'        elements which do not emit results.')
 9301 format(15x,'> Note: Further messages of this nature',
     &     /,15x '        are suppressed...')
 9400 format(15x,i7)
c
 9700 format(5x,'> local debug in do_flat_element_step: ',
     & /,10x,'step: ',i5,/,10x,'file_name: ',a,
     & /,10x,'file_no: ',i5 )
 9705 format(5x,'> local debug in do_flat_element_step: ',
     & /,10x,'proc_id: ',i5,/,10x,'file_name: ',a,
     & /,10x,'file_no: ',i5 )
 9710 format(10x,'number of data values: ',i5)
 9720 format(5x,'> data arrays allocated. no. elements: ',i7)
 9730 format( '>> FATAL ERROR: routine do_flat_element_step...',
     &  /,    '                inconsistent values.',
     &  /,    '                num_values, num_consistent: ',2i10,
     &  /,    '                file: ',a,
     &  /,    '                job terminated...' )
c
      end

c **********************************************************************
c *                                                                    *
c *          open_old_flat_result_file                                 *
c *                                                                    *
c **********************************************************************
c
c
      subroutine open_old_flat_result_file( file_no, file_name, 
     &             combine_stream, combine_text, termout, 
     &             location, open_status, dowhat )
      implicit none
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, termout, location, open_status, dowhat
      character(len=*) :: file_name
      logical :: combine_stream, combine_text
c
c                  local declarations
c                  ------------------
c
c
      if( combine_stream ) then
        open( unit=file_no, file=file_name, status='old', 
     &        access='stream', form='unformatted',
     &        iostat= open_status )
      elseif( combine_text ) then 
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
 9100 format( '>> FATAL ERROR: failed flat file open @ ',i2,
     &  /,    '                job terminated...' )
 9200 format( '>> FATAL ERROR: invalid system state @ ',i2,
     &  /,    '                job terminated...' )
c
      end
c **********************************************************************
c *                                                                    *
c *          read_flat_header                                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine read_flat_header( file_no, header_lines,
     &                             num_header_lines, max_header_lines,
     &                             termout )
      implicit none
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, termout, num_header_lines, max_header_lines
      character(len=80) :: header_lines(max_header_lines)
c
c                  local declarations
c                  ------------------
c
      character :: first_char*1, line*80
c
      header_lines(1:max_header_lines) = " "
      num_header_lines = 0
c      
      do 
        read(file_no,fmt="(a80)") line
        first_char = line(1:1)
        if( first_char .eq. "#"  .or.
     &      first_char .eq. "!"  .or.
     &      first_char .eq. "c"  .or.
     &      first_char .eq. "C"  ) then
          num_header_lines = num_header_lines + 1
          header_lines(num_header_lines) = line
          cycle
        end if     
        backspace(file_no) ! next line has number of values for states
c                            processing      
        return
       end do
c
      return
      end

c **********************************************************************
c *                                                                    *
c *          read_flat_element_values                                  *
c *                                                                    *
c **********************************************************************
c
c
      subroutine read_flat_element_values(
     &           file_no, element_values, num_values, num_elements,
     &           element_list, termout, combine_stream, 
     &           combine_text ) 
      implicit none
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, num_values, num_elements, element_list(*),
     &           termout
      logical :: combine_stream, combine_text
      double precision :: element_values(num_values,num_elements)
c
c                  local declarations
c                  ------------------
c
      integer :: read_status, element_id, now_num_values
      double precision, allocatable :: element_data(:)
c
c              next line (value) on both file types should be
c              the number of state variables. This value must
c              be same as value passed in or we have an
c              error somewhere.
c
      if( combine_stream ) read(file_no,iostat=read_status) 
     &                       now_num_values
      if( combine_text )   read(file_no,*,iostat=read_status)
     &                       now_num_values
      if( read_status .ne. 0 ) then
          write(termout,9140)
          stop
      end if     
      if( now_num_values .ne. num_values ) then
          write(termout,9150) num_values, now_num_values
        stop
      end if
c      
      allocate( element_data( num_values) )
c      
      do
        if( combine_stream ) then 
          read(unit=file_no,iostat=read_status) element_id,
     &        element_data
        elseif( combine_text ) then
          read(unit=file_no,iostat=read_status,fmt=9330)
     &         element_id, element_data
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
          write(termout,9130) element_id
          stop
        end if
c
      end do
c
      deallocate( element_data )
      return
c
 9130 format( '>> FATAL ERROR: routine read_flat_element_data...',
     &  /,    '                internal error. bad element #: ',i10,
     &  /,    '                job terminated...' )
 9140 format( '>> FATAL ERROR: routine read_flat_element_data...',
     &  /,    '                internal error. ',
     &  /,    '                job terminated...' )
 9150 format( '>> FATAL ERROR: routine read_flat_element_data...',
     &  /,    '                internal error. number of values',
     &  /,    '                inconsistent: ',2i10,
     &  /,    '                job terminated...' )
 9330 format(i8,30e15.6)
 9200 format( '>> FATAL ERROR: invalid system state in ',
     &  /,    '                routine: read_flat_element_data',
     &  /,    '                job terminated...' )
c
      end

c **********************************************************************
c *                                                                    *
c *          open_new_flat_result_file                                 *
c *                                                                    *
c **********************************************************************
c
c
      subroutine open_new_flat_result_file( file_no, file_name, 
     &         combine_stream, combine_text, termout, location )
      implicit none
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, termout, location
      character(len=*) :: file_name
      logical :: combine_stream, combine_text
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
      if( combine_stream ) then
        open( unit=file_no, file=local_fname, status='unknown', 
     &        access='stream', form='unformatted', iostat= open_status )
      elseif( combine_text ) then 
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
 9100 format( '>> FATAL ERROR: failed open_new_flat_result_file @ ',i2,
     &  /,    '                job terminated...' )
 9200 format( '>> FATAL ERROR: invalid system state @ ',i2,
     &  /,    '                file: open_new_flat_result_file',
     &  /,    '                job terminated...' )
 9300 format( '>> FATAL ERROR: invalid system state in routine:',
     &  /,    '                open_new_flat_result_file',
     &  /,    '                bad file name: ',a,
     &  /,    '                job terminated...' )
c
      end
c **********************************************************************
c *                                                                    *
c *          write_flat_element_results                                *
c *                                                                    *
c **********************************************************************
c
c
      subroutine write_flat_element_results( 
     &         file_no, num_header_lines, header_lines, num_values, 
     &         num_elements, element_values, termout, combine_stream, 
     &         combine_text, element_list ) 
      implicit none 
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, termout, num_values, num_elements,
     &           num_header_lines 
      logical :: combine_stream, combine_text
      double precision :: element_values(num_values,num_elements)
      integer :: element_list(num_elements)
      character(len=80) :: header_lines(*)
c
c                  local declarations
c                  ------------------
c
      integer ::  element_id, write_status, i
c       
      if( combine_stream ) then 
c
         do element_id = 1, num_elements
          if( element_list(element_id) .eq. 1 )
     &        write(unit=file_no,iostat=write_status)
     &        (element_values(i,element_id), i = 1, num_values)
          call check_write_status( file_no, write_status, termout, 10 )
         end do
c
      elseif( combine_text ) then
c
         do i = 1, num_header_lines
          write(file_no,9000) trim(header_lines(i))
         end do
         do element_id = 1, num_elements
          if( element_list(element_id) .eq. 1 )
     &       write(unit=file_no,iostat=write_status,fmt=9200)
     &         ( element_values(i,element_id), i = 1, num_values )
          call check_write_status( file_no, write_status, termout, 8 )
         end do
      else
        write(termout,9300)
        stop
c           
      end if
c
      return
c
 9000 format(a)
 9300 format( '>> FATAL ERROR: invalid system state in routine:',
     &  /,    '                write_element_results',
     &  /,    '                job terminated...' )
 9200 format(30e15.6)
c
      end

c **********************************************************************
c *                                                                    *
c *          do_flat_nodal_step                                        *
c *                                                                    *
c **********************************************************************
c
c
      subroutine do_flat_nodal_step( step, combine_stream, combine_text,
     &      termout, sig_eps_type, file_no, missing_flg,
     &      did_flat_stream, did_flat_text, consistent_num_values )
c
      use global_data, only :  num_model_nodes, num_model_elems,
     &                         last_proc_id
      implicit none
c
c
c                parameter declarations
c                ----------------------
c
      integer :: step, termout, sig_eps_type, file_no, 
     &           consistent_num_values
      logical :: combine_stream, combine_text, missing_flg,
     &           did_flat_stream, did_flat_text 
c
c                local variables
c                ---------------    
c
      integer:: proc_id, open_status, num_values, allocate_status, 
     &          num_nodes, node_id, miss_count, num_header_lines
      integer, allocatable, dimension(:)  :: count_list
      integer, parameter :: max_header_lines = 50 
c
      logical ::  local_debug, print_header, print
c
      double precision, allocatable,  dimension(:,:) :: node_values
      double precision, parameter :: zero = 0.0d00
c
      character(len=80) :: node_stress_fname, node_strain_fname,
     &               element_stress_fname, element_strain_fname,
     &               element_states_fname, matl_states_fname,
     &               header_lines(max_header_lines), file_name
c
      local_debug = .false.
c
c                 get name of root processor file for this load step.
c
      matl_states_fname(1:) = "..dummy.."
      call build_flat_file_names( step, 0, node_stress_fname,
     &    node_strain_fname, element_stress_fname,
     &    element_strain_fname, matl_states_fname, 
     &    element_states_fname, termout,  combine_stream, combine_text,
     &    file_name, sig_eps_type, 1 )
c
      if( local_debug ) write(termout,9700) step, file_name(1:20),
     &                                      file_no 
c
c                open the node result file for processor
c                zero. read/store header lines for text result
c
      call open_old_flat_result_file( file_no, file_name, 
     &    combine_stream, combine_text, termout, 1, open_status, 0 )
      num_header_lines = 0
      if( combine_text ) call read_flat_header( file_no, header_lines,
     &                       num_header_lines, max_header_lines,
     &                       termout )
c
      if( local_debug ) 
     &     write(termout,*) '.... num_header_lines: ', num_header_lines
c
c                 allocate space to hold the data values for all
c                 nodes in the model. keep track of which nodes
c                 have been read from the files for the step.
c
      if( combine_stream ) read(file_no)   num_values
      if( combine_text )   read(file_no,*) num_values
      close( unit=file_no )
      if( consistent_num_values == -1 ) 
     &    consistent_num_values = num_values
      if( num_values .ne. consistent_num_values ) then
        write(termout,9730) num_values, consistent_num_values,
     &                      file_name(1:len_trim(file_name))  
        stop
      end if
 
      num_nodes = num_model_nodes
      allocate( node_values(num_nodes,num_values), stat=allocate_status)
      call check_allocate_status( termout, 10001, allocate_status )
c
      allocate( count_list(num_nodes), stat=allocate_status )
      call check_allocate_status( termout, 1002, allocate_status )
c
      node_values(1:num_nodes,1:num_values) = zero
      count_list(1:num_nodes)               = 0
      if( local_debug ) write(termout,9720) num_nodes
c
c                loop over all the separate rank files of 
c                node results for this load step. open the file and read 
c                the node results into the single array here. keep
c                track of how many times a node appears in result files 
c                to enable averaging.
c
      do proc_id = 0, last_proc_id ! MPI ranks
c
        call build_flat_file_names( step, proc_id, node_stress_fname,
     &     node_strain_fname, element_stress_fname,
     &     element_strain_fname,  matl_states_fname, 
     &     element_states_fname, termout, combine_stream, combine_text,
     &     file_name, sig_eps_type, 1 )
c
c                open the node result file for the processor.
c                read all data for a stream or text file. if
c                the open fails, keep looking for processor files
c                until last_proc_id in case user blew away a
c                processor file.
c
        call open_old_flat_result_file( file_no, file_name, 
     &     combine_stream, combine_text, termout, 2, open_status, 1 )
        if( open_status .ne. 0 ) cycle
        if( local_debug ) write(termout,9705) proc_id, file_name(1:20),
     &                                        file_no
c
c               skip header lines for text results files
c
       if( combine_text ) call skip_comment_lines( file_no )
c
c                read node results for this processor. check
c                read data values until end of file. stuff into the
c                single node results table. update count of
c                node appearances.
c
       call read_flat_node_values( file_no, node_values, num_values,
     &    num_nodes, count_list, termout, combine_stream, combine_text ) 
     &   
       close( unit=file_no )
c
      end do  !  over processor files
c
c                issue warning messages for nodes that have no results.
c                disable printing of messages in curent version based
c                on a global message in calling routine.
c
      print_header = .true.
      miss_count   = 0
      print        = .true.
c
      do node_id = 1, num_nodes
         if( count_list(node_id) .eq. 0 ) missing_flg = .true.
      end do
c
      if( print ) then
       do node_id = 1, num_nodes
         if( count_list(node_id) .eq. 0 ) then
           if( print_header ) write(termout,9300)
           print_header = .false.
           write(termout,9400) node_id
           miss_count = miss_count + 1
           if( miss_count .eq. 50 ) then
              write(termout,9301)
              exit
           end if
         end if ! on node_list
       end do
      end if ! on print
c
c               node results data looks ok. create a new,
c               single flat results file for this step.
c               write in all the data.
c
      if( combine_stream ) did_flat_stream = .true.
      if( combine_text )   did_flat_text   = .true.
c      
      call open_new_flat_result_file( file_no, file_name, 
     &     combine_stream, combine_text, termout, 2 )
c
      call write_flat_node_results( file_no, num_header_lines,
     &      header_lines, num_values, num_nodes, node_values,  
     &      count_list, termout, combine_stream, combine_text,
     &      sig_eps_type )
c    
      close( unit=file_no )
c
      deallocate( node_values, count_list, stat=allocate_status )
      call check_deallocate_status( termout, 1003, allocate_status )
c
      return    
c
 9300 format(10x,'> Note: The following nodes have no results',
     &     /,10x '        in the flat result files...')
 9301 format(10x,'> Note: Further error messages of this nature',
     &     /,10x '        are suppressed...')
 9400 format(15x,i7)
c
 9700 format(5x,'> local debug in do_flat_node_step: ',
     & /,10x,'step: ',i5,/,10x,'file_name: ',a20,
     & /,10x,'file_no: ',i5 )
 9705 format(5x,'> local debug in do_flat_node_step: ',
     & /,10x,'proc_id: ',i5,/,10x,'file_name: ',a20,
     & /,10x,'file_no: ',i5 )
 9720 format(5x,'> data arrays allocated. no. nodes: ',i7)
 9730 format( '>> FATAL ERROR: routine do_flat_nodal_step...',
     &  /,    '                inconsistent values.',
     &  /,    '                num_values, num_consistent: ',2i10,
     &  /,    '                file: ',a,
     &  /,    '                job terminated...' )
c
      end
c **********************************************************************
c *                                                                    *
c *          read_flat_node_values                                     *
c *                                                                    *
c **********************************************************************
c
c
      subroutine read_flat_node_values( file_no, node_values, 
     &             num_values, num_nodes, node_counts, termout, 
     &             combine_stream, combine_text ) 
      implicit none
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, num_values, num_nodes, node_counts(*),
     &           termout
      logical :: combine_stream, combine_text
      double precision :: node_values(num_nodes,num_values)
c
c                  local declarations
c                  ------------------
c
      integer :: read_status, node_id, lcount, now_num_values
      double precision :: node_data(100)
c
c              next line (value) on both file types should be
c              the number of state variables. This value must
c              be same as value passed in or we have an
c              error somewhere.
c
      if( combine_stream ) read(file_no,iostat=read_status) 
     &                       now_num_values
      if( combine_text )   read(file_no,*,iostat=read_status)
     &                       now_num_values
      if( read_status .ne. 0 ) then
          write(termout,9140)
          stop
      end if     
      if( now_num_values .ne. num_values ) then
          write(termout,9150) num_values, now_num_values
        stop
      end if
c
      do
        if( combine_stream ) then 
          read(unit=file_no,iostat=read_status) node_id, lcount,
     &          node_data(1:num_values)
        elseif( combine_text ) then
          read(unit=file_no,iostat=read_status,fmt=9330) node_id,
     &          lcount, node_data(1:num_values)
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
          node_counts(node_id) = node_counts(node_id) + lcount
        else
          write(termout,9130)
          stop
        end if
c
      end do
c
      return
c
 9130 format( '>> FATAL ERROR: routine read_flat_node_values...',
     &  /,    '                internal error @ 5',
     &  /,    '                job terminated...' )
 9140 format( '>> FATAL ERROR: routine read_flat_nodal_data...',
     &  /,    '                internal error. ',
     &  /,    '                job terminated...' )
 9150 format( '>> FATAL ERROR: routine read_flat_nodal_data...',
     &  /,    '                internal error. number of values',
     &  /,    '                inconsistent: ',2i10,
     &  /,    '                job terminated...' )
 9330 format(2i9,30d15.6)
 9200 format( '>> FATAL ERROR: invalid system state in ',
     &  /,    '                routine: read_flat_node_values',
     &  /,    '                job terminated...' )
c
      end

c **********************************************************************
c *                                                                    *
c *          write_flat_node_results                                   *
c *                                                                    *
c **********************************************************************
c
c
      subroutine write_flat_node_results( file_no, num_header_lines,
     &         header_lines, num_values, num_nodes, node_values, 
     &         node_counts,termout, combine_stream, combine_text, 
     &         sig_eps_type ) 
      implicit none 
c
c                  parameter declarations
c                  ----------------------
c
      integer :: file_no, termout, num_values, num_nodes,
     &           sig_eps_type, num_header_lines
      integer :: node_counts(*)
      logical :: combine_stream, combine_text
      double precision :: node_values(num_nodes,num_values)
      character(len=80) :: header_lines(*)
c
c                  local declarations
c                  ------------------
c
      integer :: node_id, write_status, i
      double precision :: scale
      double precision, parameter :: zero = 0.0d0, one = 1.0d0
c
      do node_id = 1, num_nodes
        scale = zero
        if( node_counts(node_id) .gt. 0 ) 
     &       scale = one /  dble( node_counts(node_id) )
        node_values(node_id,1:num_values) =
     &            node_values(node_id,1:num_values) * scale
      end do      
c
      call compute_extra_sig_eps_values( node_values, num_nodes,
     &                                   sig_eps_type )
c       
      if( combine_stream ) then 
c
         do node_id = 1, num_nodes
           write(unit=file_no,iostat=write_status)
     &          node_values(node_id,1:num_values)
           call check_write_status( file_no, write_status, termout,
     &                              200 )  
         end do
c
      elseif( combine_text ) then  
c
         do i = 1, num_header_lines
          write(file_no,9000) trim( header_lines(i) )
         end do
c         
         do node_id = 1, num_nodes
            write(unit=file_no,iostat=write_status,fmt=9330)
     &            node_values(node_id,1:num_values)
            call check_write_status( file_no, write_status, termout,
     &                               201 )
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
 9000 format(a)
 9200 format( '>> FATAL ERROR: invalid system state in routine:',
     &  /,    '                write_flat_node_results',
     &  /,    '                job terminated...' )
 9330 format(30e15.6)
c
      end
