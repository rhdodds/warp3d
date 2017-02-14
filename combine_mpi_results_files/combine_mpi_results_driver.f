c     ****************************************************************
c     *                                                              *
c     *                    f-90 module main_data                     *
c     *                                                              *
c     ****************************************************************
c
c
c
      module size_parameter_data
c
      integer, parameter :: max_states_matls = 100 
c
      end module size_parameter_data

      module global_data
      use size_parameter_data, only : max_states_matls
c

      integer :: num_model_nodes, num_model_elems, last_step_no,
     &           last_proc_id
      logical :: combine_nodal_results, combine_element_results, 
     &           combine_patran_binary, 
     &           combine_patran_formatted,
     &           found_ele_stresses, found_ele_strains,
     &           found_node_stresses, found_node_strains,
     &           found_ele_states, did_patran_binary,
     &           did_patran_formatted, did_flat_stream, did_flat_text,     
     &           nodes_with_no_patran_values, 
     &           elements_with_no_patran_values,
     &           nodes_with_no_flat_values, 
     &           elements_with_no_flat_values, combine_flat_stream,
     &           combine_flat_text

      integer :: termin, termout, count, states_count
      integer, allocatable, dimension(:) ::
     &           nodal_stresses_step_list,  
     &           nodal_strains_step_list,
     &           element_stresses_step_list,
     &           element_strains_step_list,
     &           element_states_step_list
      character(len=80) :: states_materials(max_states_matls) 
c
      end module global_data

c     ****************************************************************
c     *                                                              *
c     *    main program for the program to combine flat & patran     *
c     *    result files for each step produced by mpi processes into *
c     *    single,                                                   *
c     *                                                              *
c     *               can run in parallel using OpenMP               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 2/11/2017 rhd              *
c     *                                                              *
c     ****************************************************************
c
      program combine_results 
      use size_parameter_data
      use global_data
      implicit none
c
      integer, external :: omp_get_max_threads
c
      termin  = 5
      termout = 6
c
      write(termout,9000)
c
c                command line must have number of model nodes, elements 
c                as parameters
c
c                <command> -nodes ### -elements ###   
c                <command> -elements ### -nodes ###
c      
      call get_model_sizes( termout, num_model_nodes, num_model_elems )
c
c                find states_header_... files in curent directory.
c                save names of material models in a list with count.
c                needed to make sure we process all states output
c                for different material models
c
c                also find the (1) highest numbered step with results
c                files, (2) highest numbered MPI rank
c
      call extract_directory_info( termout, states_count,
     &   states_materials, max_states_matls, last_step_no,
     &   last_proc_id )
      if( last_step_no == 0 ) then ! no result files
         write(termout,*) " "
         write(termout,9150)
         stop
      end if 
      if( last_proc_id == 0 ) then ! no multi rank MPI results
         write(termout,*) " "
         write(termout,9160)
         stop
      end if 

c      
      write(termout,9010) omp_get_max_threads() 
      write(termout,9020) num_model_nodes, num_model_elems,
     &                    last_step_no, last_proc_id
      write(termout,*) ' '
c
      combine_nodal_results          = .true.
      combine_element_results        = .true.

      allocate( nodal_stresses_step_list(last_step_no),  
     &          nodal_strains_step_list(last_step_no), 
     &          element_stresses_step_list(last_step_no), 
     &          element_strains_step_list(last_step_no),  
     &          element_states_step_list(last_step_no) )  
c      
c                process Patran files
c
      nodes_with_no_patran_values    = .false.
      elements_with_no_patran_values = .false. 
      did_patran_binary              = .false.
      did_patran_formatted           = .false.
      call process_patran_binary
      call process_patran_formatted
c
c                process flat files
c
      nodes_with_no_flat_values      = .false.
      elements_with_no_flat_values   = .false.
      did_flat_stream                = .false.
      did_flat_text                  = .false.
      call process_flat_stream
      call process_flat_text
c      
c                all done
c
      write(termout,9100)
      if( nodes_with_no_patran_values ) write(termout,9110)
      if( elements_with_no_patran_values ) write(termout,9120)
c      
      if( nodes_with_no_flat_values ) write(termout,9130)
      if( elements_with_no_flat_values ) write(termout,9140)
c
      write(termout,*) ' '
      write(termout,*) '>> Normal termination...'
      write(termout,*) ' '
      stop
c
 9000 format(/,
     &       ' ****************************************************',
     &  /,   ' *                                                  *',
     &  /,   ' *     Combine Nodal and Element Result Files       *',
     &  /,   ' *     From File Parts Created During Parallel      *',
     &  /,   ' *     MPI Execution of WARP3D                      *',
     &  /,   ' *                                                  *',
     &  /,   ' *             build date:  2-10-2017               *',
     &  /,   ' *                                                  *',
     &  /,   ' *                                                  *',
     &  /,   ' ****************************************************'
     &  // )
 9010 format('... Using # of threads:                ',i8)
 9020 format('... Number of model nodes, elements, steps, MPI ranks: ',
     &   3i8,i5)
 9100 format(//,'.... All results processed ....' )
 9110 format(/,'>>>>> Warning: 1 or more nodes in the Patran nodal',
     & ' files',
     & /,15x,'had no results. This could be, for example, ',
     &   'from cohesive',
     & /,15x,'elements with nodes on a symmetry plane -- since ',
     &    'cohesive',
     & /,15x,'elements make no contributions to nodal strain/stress',
     & /,15x,'results. Patran result values for these nodes are zero.')
 9120 format(/,'>>>>> Warning: 1 or more elements in the Patran',
     & '  element files',
     & /,15x,'had no results. Patran result values for these',
     & /,15x,'elements are zero.')
 9130 format(/,'>>>>> Warning: 1 or more nodes in the flat nodal',
     & ' files',
     & /,15x,'had no results. This could be, for example, ',
     &   'from cohesive',
     & /,15x,'elements with nodes on a symmetry plane -- since ',
     &    'cohesive',
     & /,15x,'elements make no contributions to nodal strain/stress',
     & /,15x,'results. flat result values for these nodes are zero.')
 9140 format(/,'>>>>> Warning: 1 or more elements in the flat element',
     & ' files',
     & /,15x,'had no results. flat result values for these',
     & /,15x,'elements are zero.')
 9150 format(/,'>>>>> Error: no files in current directory have a step',
     & /,13x,'number in position 5-9 of the name cocnistent with',
     & /,13x,'WARP3D naming of Patran and flat results files.',
     & /,13x,'This program must be executed in directory having',
     & /,13x,'the results.',/,13x,'Program terminated...'//)
 9160 format(/,'>>>>> Error: no files in current directory have an MPI',
     & /,13x,'rank number > 0. If only 1 MPI rank, no combination',
     & /,13x,'oerpations are possible...',
     & /,13x,'Program terminated...'//)
c
      contains
c     ========
c
c **********************************************************************
c *      contains:       process_flat_text                             *
c **********************************************************************
c
      subroutine process_flat_text
      implicit none
c      
      character(len=80) :: dummy_name 
c
      write(termout,*) ' '
      write(termout,*) '>> Processing flat text result files...'
c      
      combine_flat_text   = .true.
      combine_flat_stream = .not. combine_flat_text
c      
      dummy_name(1:) = "..dummy.."
      call find_flat_work( combine_flat_stream,  combine_flat_text,
     &                     count, dummy_name )
c      
      if( count > 0 ) then
        if( combine_nodal_results ) 
     &      call combine_flat_nodal(
     &           combine_flat_stream, combine_flat_text )
        if( combine_element_results ) 
     &      call combine_flat_element_sigeps(
     &           combine_flat_stream, combine_flat_text )
      end if
c
      if( states_count > 0 .and. combine_element_results )
     &    call combine_flat_states( combine_flat_stream,
     &                              combine_flat_text )
c
      if( .not. did_flat_text ) 
     &  write(termout,*) '      > no text files to process'
c
      return     
c
      end subroutine process_flat_text
c
c **********************************************************************
c *      contains:       process_flat_stream                           *
c **********************************************************************
c
      subroutine process_flat_stream
      implicit none
c
      character(len=80) ::  dummy_name 

      write(termout,*) ' '
      write(termout,*) '>> Processing flat stream result files...'
c      
      combine_flat_text   = .false.
      combine_flat_stream = .not. combine_flat_text
c
      dummy_name(1:) = "..dummy.."
      call find_flat_work( combine_flat_stream,  combine_flat_text,
     &                     count, dummy_name )
      if( count > 0 ) then
        if( combine_nodal_results ) 
     &      call combine_flat_nodal(
     &           combine_flat_stream, combine_flat_text )
        if( combine_element_results ) 
     &      call combine_flat_element_sigeps(
     &           combine_flat_stream, combine_flat_text )
      end if
c
      if( states_count > 0 .and. combine_element_results )
     &    call combine_flat_states( combine_flat_stream,
     &                              combine_flat_text )
c
      if( .not. did_flat_stream ) 
     &  write(termout,*) '      > no text files to process'
c
      return
c
      end subroutine process_flat_stream
c
c **********************************************************************
c *      contains:       process_patran_binary                         *
c **********************************************************************
c
      subroutine process_patran_binary
      implicit none
c      
      character(len=80) :: dummy_name 
c
c                process Patran binary node and element results files.
c
      write(termout,*) '>> Processing Patran binary result files...'
      combine_patran_binary    = .true.
      combine_patran_formatted = .false.
c      
      dummy_name(1:) = "..dummy.."
      call find_patran_work( combine_patran_binary, 
     &                  combine_patran_formatted, count, dummy_name )      
c      
      if( count > 0 ) then
        if( combine_nodal_results ) 
     &      call combine_patran_nodal(
     &           combine_patran_binary, combine_patran_formatted )
c     
        if( combine_element_results ) 
     &      call combine_patran_element_sig_eps(
     &           combine_patran_binary, combine_patran_formatted )
      end if
c     
      if( states_count > 0 .and. combine_element_results ) 
     &      call combine_patran_states( combine_patran_binary,
     &                                  combine_patran_formatted )
c   
      if( .not. did_patran_binary ) 
     &  write(termout,*) '      > no binary files to process'
c     
      return
      end subroutine process_patran_binary

c
c **********************************************************************
c *      contains:       process_patran_formatted                      *
c **********************************************************************
c
      subroutine process_patran_formatted
c      
      character(len=80) :: dummy_name 
c 
c                process Patran formatted node and element results
c                files.
c
      write(termout,*) ' '
      write(termout,*) '>> Processing Patran formatted result files...'
      combine_patran_binary    = .false.
      combine_patran_formatted = .true.
c 
      dummy_name(1:) = "..dummy.."
c     
      call find_patran_work( combine_patran_binary, 
     &                 combine_patran_formatted, count, dummy_name )
c
      if( count > 0 ) then
        if( combine_nodal_results ) 
     &      call combine_patran_nodal(
     &           combine_patran_binary, combine_patran_formatted )
c     
        if( combine_element_results ) 
     &      call combine_patran_element_sig_eps( 
     &           combine_patran_binary, combine_patran_formatted )
      end if
c
      if( states_count > 0 .and. combine_element_results ) 
     &         call combine_patran_states( combine_patran_binary,
     &                                     combine_patran_formatted )
c     
      if( .not. did_patran_formatted ) 
     &  write(termout,*) '      > no formatted files to process'
c     
      return
      end subroutine process_patran_formatted

c
c **********************************************************************
c *      contains:       combine_patran_nodal                          *
c **********************************************************************
c
      subroutine combine_patran_nodal( combine_binary,
     &                                 combine_formatted )
      implicit none
c
c                parameter declarations
c                ----------------------
c
      logical :: combine_binary, combine_formatted
c
c                local variables
c                ---------------    
c
      integer, external ::  OMP_GET_THREAD_NUM  
      integer :: step, file_no, my_thread, consistent_num_values
      logical :: found_nodes_no_values
c
      if( found_node_stresses ) then
c
      write(termout,*) '   > Processing nodal stresses...'
      consistent_num_values = -1      
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_nodes_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, last_step_no
         if( nodal_stresses_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_nodes_no_values = .false.
         call do_patran_nodal_step( step, combine_binary,
     &      combine_formatted, termout, 0, file_no, 
     &      found_nodes_no_values, did_patran_binary,
     &      did_patran_formatted, consistent_num_values )
c$OMP ATOMIC UPDATE
         nodes_with_no_patran_values = nodes_with_no_patran_values .or.
     &                                 found_nodes_no_values
         write(termout,*) '      > completed step:   ',step
        end do
c$OMP END DO
c$OMP END PARALLEL

      end if ! processed nodal stresses
c
c
      if( found_node_strains ) then
c
        write(termout, *) '   > Processing nodal strains...'
        consistent_num_values = -1
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_nodes_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, last_step_no
         if( nodal_strains_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_nodes_no_values = .false.
         call do_patran_nodal_step( step, combine_binary, 
     &      combine_formatted, termout, 1, file_no,
     &      found_nodes_no_values, did_patran_binary,
     &      did_patran_formatted, consistent_num_values )
c$OMP ATOMIC UPDATE
         nodes_with_no_patran_values = nodes_with_no_patran_values .or.
     &                         found_nodes_no_values
         write(termout,*) '      > completed step: ',step
        end do
c$OMP END DO
c$OMP END PARALLEL
      end if   
c
      return

      end subroutine combine_patran_nodal
c **********************************************************************
c *      contains:       combine_patran_element_sig_eps                *
c **********************************************************************
c
      subroutine combine_patran_element_sig_eps( combine_binary,
     &                                           combine_formatted )
      implicit none
c
c                parameter declarations!
c                ----------------------
c
      logical :: combine_binary, combine_formatted
c
c                local variables
c                ---------------    
c
      integer, external ::  OMP_GET_THREAD_NUM  
      integer :: step, file_no, my_thread, matls_count, local_count,
     &           nchars, consistent_num_values
      logical :: found_elements_no_values
      character(len=80) :: matl_states_name
c
      matl_states_name(1:) = "..dummy.."
c      
      if( found_ele_stresses ) then
        write(termout,*) '   > Processing element stresses...'
        consistent_num_values = -1
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_elements_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, last_step_no
         if( element_stresses_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_elements_no_values = .false.
         call do_patran_element_step( step, combine_binary,
     &       combine_formatted, termout, 0, file_no,
     &       found_elements_no_values,  matl_states_name,
     &       did_patran_binary, did_patran_formatted,
     &       consistent_num_values )
c$OMP ATOMIC UPDATE
         elements_with_no_patran_values = 
     &            elements_with_no_patran_values .or.
     &            found_elements_no_values
         write(termout,*) '      > completed step: ',step
        end do
c$OMP END DO
c$OMP END PARALLEL
      end if  
c
      if( found_ele_strains ) then
        write(termout,*) '   > Processing element strains...'
        consistent_num_values = -1
cc
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_elements_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, last_step_no
         if( element_strains_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_elements_no_values = .false.
         call do_patran_element_step( step, combine_binary,
     &       combine_formatted, termout, 1, file_no,
     &       found_elements_no_values, matl_states_name, 
     &       did_patran_binary, did_patran_formatted,
     &       consistent_num_values )
c$OMP ATOMIC UPDATE
         elements_with_no_patran_values =
     &              elements_with_no_patran_values .or.
     &              found_elements_no_values
         write(termout,*) '      > completed step: ',step
        end do
c$OMP END DO
c$OMP END PARALLEL
      end if
c
      return
c
      end subroutine combine_patran_element_sig_eps
c **********************************************************************
c *      contains:       combine_patran_states                         *
c **********************************************************************
c
      subroutine combine_patran_states( combine_binary,
     &                                  combine_formatted )
      implicit none
c
c                parameter declarations
c                ----------------------
c
      logical :: combine_binary, combine_formatted
c
c                local variables
c                ---------------    
c
      integer, external ::  OMP_GET_THREAD_NUM  
      integer :: step, file_no, my_thread, matls_count, local_count,
     &           nchars, consistent_num_values
      logical :: found_elements_no_values, found_binary, process,
     &           found_formatted, check_binary, check_formatted,
     &           write_hdr
      logical, parameter :: ldebug = .false.
      character(len=80) :: matl_states_name
c
      if( ldebug ) write(termout,9010)
      write_hdr = .true.
c
c              loop over the list of material models appearing in
c              states header files. run the entire combine process 
c              for each states material.
c
c              Example: FE model has mises_gurson and  cohesive
c                 materials and there are states_headers file for each.
c                 combine all mises_gurson results from MPI ranks
c                 for each load step. Do same for cohesive.
c
      do matls_count = 1, states_count
        matl_states_name(1:) = states_materials(matls_count)(1:)
        if( ldebug ) write(termout,*)
     &    '.. loop. matls_count, matl name: ',
     &            matls_count, matl_states_name   
        nchars = len_trim( matl_states_name ) 
c            
c              see if real work to do. if processing binary but no
c              binary states results exist for this material, skip.
c              same process if we are doinf formatted results but
c              none exists for this material.
c
        found_binary = .false.
        found_formatted = .false.
        check_binary = .true.
        check_formatted = .false.
        call find_patran_work( check_binary, check_formatted, 
     &                         local_count, matl_states_name )
        if( found_ele_states ) found_binary = .true.
        check_binary = .false.
        check_formatted = .true.
        call find_patran_work( check_binary, check_formatted, 
     &                         local_count, matl_states_name )
        if( found_ele_states ) found_formatted = .true.
        process = ( combine_binary .and. found_binary ) .or.
     &            ( combine_formatted .and. found_formatted ) 
        if( .not. process ) cycle
c        
c              binary or formmated states result are to be
c              combined for this material type. get the list of
c              load steps that have required results.
c
        call find_patran_work( combine_binary, combine_formatted, 
     &                         local_count, matl_states_name )
        local_count = 0
        do step = 1, last_step_no
          if( element_states_step_list(step) .ne. 0 )
     &          local_count = local_count + 1
        end do
        if( ldebug ) write(termout,9020) local_count
        if( write_hdr ) then
         write(termout,*) '   > Processing element states...'
         write_hdr = .false.
        end if
c        
c              for this material states, eg., mises_gurson, combine
c              MPI processor results for each step into a single file
c              containing results for all elements
c
        consistent_num_values = -1
c
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_elements_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, last_step_no
         if( element_states_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_elements_no_values = .false.
         call do_patran_element_step( step, combine_binary,
     &       combine_formatted, termout, 2, file_no,
     &       found_elements_no_values, matl_states_name,
     &       did_patran_binary, did_patran_formatted, 
     &       consistent_num_values )
c$OMP ATOMIC UPDATE
         elements_with_no_patran_values =
     &              elements_with_no_patran_values .or.
     &              found_elements_no_values
         write(termout,*) '      > completed step: ',step
        end do
c$OMP END DO
c$OMP END PARALLEL
c    
      end do ! matls_count

      if( ldebug ) write(termout,9030)
      return
c
 9005 format(1x,"FATAL ERROR. inconsistency in routine",
     & /,1x,    "             combine_patran_states",
     & /,1x,    "             states_count = 0",
     & /,1x,    "             Job Aborted.") 
 9010 format(1x,"... enter combine_patran_states....") 
 9020 format(10x,'local_count: ',i10)     
 9030 format(1x,"... leave combine_patran_states.... ")
 9050 format(1x,"WARNING: inconsistent Patran elements states..",
     &  /,    1x,"        no results found for material: ",a,
     &  /,    1x,"        but states_header file exists..",
     &  /,    1x,"        skipping this material states",/)          

      end subroutine combine_patran_states
c
c **********************************************************************
c *      contains:       combine_flat_element_sigeps                   *
c **********************************************************************
c
      subroutine combine_flat_element_sigeps( combine_stream,
     &                                        combine_text )
      implicit none
c      
      logical :: combine_stream, combine_text
c
      integer, external ::  OMP_GET_THREAD_NUM  
      integer :: step, file_no, my_thread,  consistent_num_values
      logical :: found_elements_no_values
      character(len=80) :: matl_states_name
c
      matl_states_name(1:) = "..dummy.."
c
      if( found_ele_stresses ) then
        consistent_num_values = -1
        write(termout,*) '   > Processing element stresses...'
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_elements_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, last_step_no
         if( element_stresses_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_elements_no_values = .false.
         call do_flat_element_step( step, combine_stream,
     &      combine_text, termout, 0, file_no, found_elements_no_values,
     &      matl_states_name, did_flat_stream, did_flat_text,
     &      consistent_num_values )
c$OMP ATOMIC UPDATE
         elements_with_no_flat_values = 
     &            elements_with_no_flat_values .or.
     &            found_elements_no_values
         write(termout,*) '      > completed step: ',step
        end do
c$OMP END DO
c$OMP END PARALLEL
      end if  
c
      if( found_ele_strains ) then
        consistent_num_values = -1
        write(termout,*) '   > Processing element strains...'
c
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_elements_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, last_step_no
         if( element_strains_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_elements_no_values = .false.
         call do_flat_element_step( step, combine_stream,
     &      combine_text, termout, 1, file_no, found_elements_no_values,
     &      matl_states_name, did_flat_stream, did_flat_text,
     &      consistent_num_values )
c$OMP ATOMIC UPDATE
         elements_with_no_flat_values =
     &              elements_with_no_flat_values .or.
     &              found_elements_no_values
         write(termout,*) '      > completed step: ',step
        end do
c$OMP END DO
c$OMP END PARALLEL
      end if
c    
      return
      end subroutine combine_flat_element_sigeps

c **********************************************************************
c *      contains:      combine_flat_nodal                             *
c **********************************************************************
c
      subroutine combine_flat_nodal( combine_stream, combine_text )
      implicit none
c
      logical :: combine_stream, combine_text
c      
      integer, external ::  OMP_GET_THREAD_NUM  
      integer :: step, file_no, my_thread, consistent_num_values
      logical :: found_nodes_no_values
c
      if( found_node_stresses ) then
        consistent_num_values = -1
c
        write(termout,*) '   > Processing nodal stresses...'
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_nodes_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, last_step_no
         if( nodal_stresses_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_nodes_no_values = .false.
         call do_flat_nodal_step( step, combine_stream, combine_text,
     &        termout, 0, file_no, found_nodes_no_values,
     &        did_flat_stream, did_flat_text, consistent_num_values )
c$OMP ATOMIC UPDATE
         nodes_with_no_flat_values = nodes_with_no_flat_values .or.
     &                                 found_nodes_no_values
         write(termout,*) '      > completed step:   ',step
        end do
c$OMP END DO
c$OMP END PARALLEL

      end if ! processed nodal stresses
c
c
      if( found_node_strains ) then
        consistent_num_values = -1
c
        write(termout, *) '   > Processing nodal strains...'
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_nodes_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, last_step_no
         if( nodal_strains_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_nodes_no_values = .false.
         call do_flat_nodal_step( step, combine_stream, 
     &      combine_text, termout, 1, file_no, found_nodes_no_values,
     &      did_flat_stream, did_flat_text, consistent_num_values )
c$OMP ATOMIC UPDATE
         nodes_with_no_flat_values = nodes_with_no_flat_values .or.
     &                         found_nodes_no_values
         write(termout,*) '      > completed step: ',step
        end do
c$OMP END DO
c$OMP END PARALLEL
      end if   
c
      return

      end subroutine combine_flat_nodal
c **********************************************************************
c *      contains:       combine_flat_states                           *
c **********************************************************************
c
      subroutine combine_flat_states( combine_stream, combine_text )
      implicit none
c
c                parameter declarations
c                ----------------------
c
      logical :: combine_stream, combine_text
c
c                local variables
c                ---------------    
c
      integer, external ::  OMP_GET_THREAD_NUM  
      integer :: step, file_no, my_thread, matls_count, local_count,
     &           nchars, consistent_num_values
      logical :: found_elements_no_values, found_stream, process,
     &           found_text, check_stream, check_text, write_hdr
      logical, parameter :: ldebug = .false.
      character(len=80) :: matl_states_name
c
      if( ldebug ) write(termout,9010)
      write_hdr = .true.
c
c              loop over the list of material models appearing in
c              states header files. run the entire combine process 
c              for each states material.
c
c              Example: FE model has mises_gurson and  cohesive
c                 materials and there are states_headers file for each.
c                 combine all mises_gurson results from MPI ranks
c                 for each load step. Do same for cohesive.
c
      do matls_count = 1, states_count
c        
        matl_states_name(1:) = states_materials(matls_count)(1:)
        if( ldebug ) write(termout,*)
     &    '.. loop. matls_count, matl name: ', matls_count, 
     &      matl_states_name(1:len_trim(matl_states_name))   
        nchars = len_trim( matl_states_name ) 
c            
c              see if real work to do. if processing binary but no
c              binary states results exist for this material, skip.
c              same process if we are doinf formatted results but
c              none exists for this material.
c
        found_stream  = .false.
        found_text    = .false.
        check_stream  = .true.
        check_text    = .false.
        call find_flat_work( check_stream, check_text, local_count,
     &                       matl_states_name )
        if( found_ele_states ) found_stream = .true.
        check_stream = .false.
        check_text   = .true.
        call find_flat_work( check_stream, check_text, local_count,
     &                       matl_states_name )
        if( found_ele_states ) found_text = .true.
        process = ( combine_stream .and. found_stream ) .or.
     &            ( combine_text .and. found_text )

        write(*,*) '.. check_stream: ', check_stream
        write(*,*) '.. check_text: ',   check_text
        write(*,*) '.. found_stream: ', found_stream
        write(*,*) '.. found_text: ',   found_text
        write(*,*) '.. process: ', process

        if( .not. process ) cycle
c        
c              stream or text states result are to be
c              combined for this material type. get the list of
c              load steps that have required results.
c
        call find_flat_work( combine_stream, combine_text, 
     &                       local_count, matl_states_name )
        local_count = 0
        do step = 1, last_step_no
          if( element_states_step_list(step) .ne. 0 )
     &          local_count = local_count + 1
        end do
        if( ldebug ) write(termout,9020) local_count
        if( write_hdr ) then
         write(termout,*) '   > Processing element states...'
         write_hdr = .false.
        end if
        consistent_num_values  = -1      
c        
c              for this material states, eg., mises_gurson, combine
c              MPI processor results for each step into a single file
c              containing results for all elements
c
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_elements_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, last_step_no
         if( element_states_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_elements_no_values = .false.
         call do_flat_element_step( step, combine_stream,
     &       combine_text, termout, 2, file_no,
     &       found_elements_no_values, matl_states_name,
     &       did_flat_stream, did_flat_text, consistent_num_values )
c$OMP ATOMIC UPDATE
         elements_with_no_flat_values =
     &              elements_with_no_flat_values .or.
     &              found_elements_no_values
         write(termout,*) '      > completed step: ',step
        end do
c$OMP END DO
c$OMP END PARALLEL
c    
      end do ! matls_count
c
      if( ldebug ) write(termout,9030)
      return
c
 9005 format(1x,"FATAL ERROR. inconsistency in routine",
     & /,1x,    "             combine_flat_states",
     & /,1x,    "             states_count = 0",
     & /,1x,    "             Job Aborted.") 
 9010 format(//,1x,"... enter combine_flat_states....") 
 9020 format(10x,'local_count: ',i10)     
 9030 format(1x,"... leave combine_flat_states.... ")
 9050 format(1x,"WARNING: inconsistent flat elements states..",
     &  /,    1x,"        no results found for material: ",a,
     &  /,    1x,"        but states_header file exists..",
     &  /,    1x,"        skipping this material states",/)          

      end subroutine combine_flat_states
c
      end program combine_results
