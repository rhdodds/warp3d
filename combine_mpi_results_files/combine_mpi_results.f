c     ****************************************************************
c     *                                                              *
c     *                    f-90 module main_data                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     ****************************************************************
c
c
c
      module size_parameter_data
c
      integer, parameter :: max_steps    = 100000
      integer, parameter :: max_elements = 8000000
      integer, parameter :: max_nodes    = 8000000
      integer, parameter :: max_procs    = 128  ! MPI ranks

      integer :: num_model_nodes, num_model_elems,
     &           num_flat_elem_strains, num_flat_elem_stresses
c
      end module

c     ****************************************************************
c     *                                                              *
c     *    main program for the program to combine flat & patran     *
c     *    result files for each step produced by mpi processes into *
c     *    single, standard result files for patran processing       *
c     *                                                              *
c     *               can run in parallel using OpenMP               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/30/2017 rhd              *
c     *                                                              *
c     ****************************************************************
c
      program combine_results
      implicit none
      call driver
      stop
      end


      subroutine driver 
      use size_parameter_data
      implicit none
c
      logical :: combine_nodal_results, combine_element_results, 
     &           combine_patran_binary, 
     &           combine_patran_formatted,
     &           found_ele_stresses, found_ele_strains,
     &           found_node_stresses, found_node_strains,
     &           nodes_with_no_patran_values, 
     &           elements_with_no_patran_values,
     &           nodes_with_no_flat_values, 
     &           elements_with_no_flat_values, combine_flat_stream,
     &           combine_flat_text

      integer :: termin, termout, count,
     &           nodal_stresses_step_list(max_steps),  
     &           nodal_strains_step_list(max_steps),
     &           element_stresses_step_list(max_steps),
     &           element_strains_step_list(max_steps)
      integer, external :: omp_get_max_threads
c
      termin  = 5
      termout = 6
      num_flat_elem_strains  = 22
      num_flat_elem_stresses = 26
c
      write(termout,9000)
c
c                command line must have number of model nodes, elements 
c                as parameters
c
c                <command> -nodes ### -elements ###   
c                <command> -elements ### -nodes ###
c      
      call get_model_sizes
      write(termout,9010) omp_get_max_threads() 
      write(termout,9020) num_model_nodes, num_model_elems
      write(termout,*) ' '
c
c
c                process Patran binary node and element results files.
c
c
      combine_nodal_results          = .true.
      combine_element_results        = .true.
      nodes_with_no_patran_values    = .false.
      elements_with_no_patran_values = .false. 
c      
      write(termout,*) '>> Processing Patran binary result files...'
      combine_patran_binary    = .true.
      combine_patran_formatted = .false.
      call driver_find_patran_work( combine_patran_binary, 
     &                              combine_patran_formatted, count )
      if( count .eq. 0 ) then ! no Patran binary files to combine 
        write(termout,*) '      > no Patran binary files found...'
        write(termout,*) ' '
      else
        if( combine_nodal_results ) 
     &      call driver_combine_patran_nodal_results(
     &           combine_patran_binary, combine_patran_formatted )
        if( combine_element_results ) 
     &      call driver_combine_patran_element_results(
     &           combine_patran_binary, combine_patran_formatted )
       end if
c
c
c                process Patran formatted node and element results
c                files.
c
c
      write(termout,*) ' '
      write(termout,*) '>> Processing Patran formatted result files...'
      combine_patran_binary    = .false.
      combine_patran_formatted = .true.
      call driver_find_patran_work( combine_patran_binary, 
     &                              combine_patran_formatted, count )
      if( count .eq. 0 ) then ! no Patran formated files to combine 
         write(termout,*) '      > no Patran formatted files found...'
         write(termout,*) ' '
      else
        if( combine_nodal_results ) 
     &      call driver_combine_patran_nodal_results(
     &           combine_patran_binary, combine_patran_formatted )
        if( combine_element_results ) 
     &      call driver_combine_patran_element_results( 
     &           combine_patran_binary, combine_patran_formatted )
      end if
c
c
c                process flat (formatted) result files of node/element 
c                results
c
      combine_nodal_results          = .true.
      combine_element_results        = .true.
      nodes_with_no_flat_values      = .false.
      elements_with_no_flat_values   = .false. 
c      
      write(termout,*) '>> Processing flat text result files...'
      combine_flat_text   = .true.
      combine_flat_stream = .not. combine_flat_text
      call driver_find_flat_work( combine_flat_stream, 
     &                            combine_flat_text, count )
      if( count .eq. 0 ) then ! no flat text files to combine 
        write(termout,*) '      > no flat text files found...'
        write(termout,*) ' '
      else
        if( combine_nodal_results ) 
     &      call driver_combine_flat_nodal_results(
     &           combine_flat_stream, combine_flat_text )
        if( combine_element_results ) 
     &      call driver_combine_flat_element_results(
     &           combine_flat_stream, combine_flat_text )
       end if
c
c
c                process flat (stream) result files of node/element 
c                results
c
      combine_nodal_results          = .true.
      combine_element_results        = .true.
c      
      write(termout,*) ' '
      write(termout,*) '>> Processing flat stream result files...'
      combine_flat_text   = .false.
      combine_flat_stream = .not. combine_flat_text
      call driver_find_flat_work( combine_flat_stream, 
     &                            combine_flat_text, count )
      if( count .eq. 0 ) then ! no flat stream files to combine 
        write(termout,*) '      > no flat stream files found...'
        write(termout,*) ' '
      else
        if( combine_nodal_results ) 
     &      call driver_combine_flat_nodal_results(
     &           combine_flat_stream, combine_flat_text )
        if( combine_element_results ) 
     &      call driver_combine_flat_element_results(
     &           combine_flat_stream, combine_flat_text )
       end if
c
c
c                all done
c
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
     &  /,   ' *             build date: 1-29-2017                *',
     &  /,   ' *                                                  *',
     &  /,   ' *                                                  *',
     &  /,   ' ****************************************************'
     &  // )
 9010 format('... using # of threads:             ',i8)
 9020 format('... number of model nodes, elements:',2i8)
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
c
      contains
c     ========

c **********************************************************************
c *      contains:       get_model_sizes                               *
c **********************************************************************
c
      subroutine get_model_sizes
      implicit none
c      
      integer :: num_command_args, nchars, status
      character(len=20) :: chars
c
      num_command_args =  command_argument_count ()
      if( num_command_args .ne. 4 ) then  
         write(termout,9000)
         stop
      end if
c      
      num_model_nodes = 0
      num_model_elems = 0
c      
      call get_command_argument( 1, chars, nchars, status )
      if( chars(1:5) .eq. '-node' ) then
        call get_command_argument( 2, chars, nchars, status ) 
        read(chars,*) num_model_nodes
      elseif( chars(1:5) .eq. '-elem' ) then
        call get_command_argument( 2, chars, nchars, status ) 
        read(chars,*) num_model_elems
      else
        write(termout,9000)
        stop         
      end if
c      
      call get_command_argument( 3, chars, nchars, status )
      if( chars(1:5) .eq. '-node' ) then
        call get_command_argument( 4, chars, nchars, status ) 
        read(chars,*) num_model_nodes
      elseif( chars(1:5) .eq. '-elem' ) then
        call get_command_argument( 4, chars, nchars, status ) 
        read(chars,*) num_model_elems
      else
        write(termout,9000)
        stop         
      end if
c
      if( num_model_nodes < 1  .or. num_model_elems < 1 ) then  
         write(termout,9000)
         stop
      end if
c
      return
c
 9000 format('>>> invalid command line. required form:',
     & /,10x,' <command>  -nodes ### -elements ###',
     & /,10x,' <command>  -elements ### -nodes ###',
     & //,10x,'... job terminated ...',//)     
c 
      end subroutine  get_model_sizes     


c **********************************************************************
c *      contains:       driver_find_flat_work                         *
c **********************************************************************
c
      subroutine driver_find_flat_work( combine_stream, combine_text, 
     &                                  local_count )
      implicit none
c
      integer :: local_count
      logical :: combine_stream, combine_text   
c
      integer :: step
      logical :: l1, l2, l3, l4      
      logical, parameter :: local_debug = .false.
      character(len=80) :: node_stress_fname, node_strain_fname,
     &                     element_stress_fname, element_strain_fname
      character(len=14) :: dummy_file_name
c
c                examine all files in directory to build lists
c                of step numbers to be combined for each type
c                of result file
c
      do step = 1, max_steps
        nodal_stresses_step_list(step)   = 0
        nodal_strains_step_list(step)    = 0
        element_stresses_step_list(step) = 0
        element_strains_step_list(step)  = 0
      end do
c
      local_count = 0
      found_node_stresses = .false.
      found_node_strains  = .false.
      found_ele_stresses  = .false.
      found_ele_strains   = .false.
c
      do step = 1, max_steps
c
        call build_flat_file_names( step, 0, node_stress_fname,
     &                         node_strain_fname, element_stress_fname,
     &                         element_strain_fname, termout,
     &                         combine_stream, combine_text,
     &                         dummy_file_name, 1, 1 )
        inquire( file=node_stress_fname,exist=l1 )
        inquire( file=node_strain_fname,exist=l2 )
        inquire( file=element_stress_fname,exist=l3 )
        inquire( file=element_strain_fname,exist=l4 )
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
        if( l1 .or. l2 .or. l3 .or. l4 ) local_count = local_count + 1
c
      end do
c
      if( .not. local_debug ) return
c
      write(termout,*) "... list of flat files to process ..." 
      write(termout,*) '      total number of files: ', local_count
      do step = 1, max_steps
        if( nodal_stresses_step_list(step) == 1 )
     &      write(termout,*) '        step: ',step,' nodal stresses'
        if( nodal_strains_step_list(step) == 1 )
     &      write(termout,*) '        step: ',' nodal strains'
        if( element_stresses_step_list(step) == 1 )
     &      write(termout,*) '        step: ',step,' element stresses'
       if( element_strains_step_list(step) == 1 )
     &      write(termout,*) '        step: ',step,' element strains'
      end do 
c
      return      
c
      end subroutine driver_find_flat_work

c **********************************************************************
c *      contains:       driver_find_patran_work                       *
c **********************************************************************
c
      subroutine driver_find_patran_work( combine_binary,
     &                                    combine_formatted, 
     &                                    local_count )
      implicit none
c
      integer :: local_count
      logical :: combine_binary, combine_formatted   
c
      integer :: step
      logical :: l1, l2, l3, l4      
      character(len=80) :: node_stress_fname, node_strain_fname,
     &                     element_stress_fname, element_strain_fname
      character(len=14) :: dummy_file_name
cc
c                examine all files in directory to build lists
c                of step numbers to be combined for each type
c                of result file
c
      do step = 1, max_steps
        nodal_stresses_step_list(step)   = 0
        nodal_strains_step_list(step)    = 0
        element_stresses_step_list(step) = 0
        element_strains_step_list(step)  = 0
      end do
c
      local_count = 0
      found_node_stresses = .false.
      found_node_strains  = .false.
      found_ele_stresses  = .false.
      found_ele_strains   = .false.
c
      do step = 1, max_steps
c
        call build_patran_file_names( step, 0, node_stress_fname,
     &                         node_strain_fname, element_stress_fname,
     &                         element_strain_fname, termout,
     &                         combine_binary, combine_formatted,
     &                         dummy_file_name, 1, 1 )
        inquire( file=node_stress_fname,exist=l1 )
        inquire( file=node_strain_fname,exist=l2 )
        inquire( file=element_stress_fname,exist=l3 )
        inquire( file=element_strain_fname,exist=l4 )
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
        if( l1 .or. l2 .or. l3 .or. l4 ) local_count = local_count + 1
c
      end do
c
      return      
c
      end subroutine driver_find_patran_work

c **********************************************************************
c *      contains:       driver_combine_patran_nodal_results           *
c **********************************************************************
c
      subroutine driver_combine_patran_nodal_results( combine_binary,
     &                   combine_formatted )
      implicit none
c
      logical :: combine_binary, combine_formatted
c      
      integer, external ::  OMP_GET_THREAD_NUM  
      integer :: step, file_no, my_thread
      logical :: found_nodes_no_values
c
      if( found_node_stresses ) then
c
      write(termout,*) '   > Processing nodal stresses...'
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_nodes_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, max_steps
         if( nodal_stresses_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_nodes_no_values = .false.
         call do_patran_nodal_step( step, combine_binary,
     &                       combine_formatted, termout, 0, file_no,
     &                       found_nodes_no_values )
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
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_nodes_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, max_steps
         if( nodal_strains_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_nodes_no_values = .false.
         call do_patran_nodal_step( step, combine_binary, 
     &                       combine_formatted, termout, 1, file_no,
     &                       found_nodes_no_values )
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

      end subroutine driver_combine_patran_nodal_results
c **********************************************************************
c *      contains:       driver_combine_patran_element_results         *
c **********************************************************************
c
      subroutine driver_combine_patran_element_results( combine_binary,
     &                   combine_formatted )
      implicit none
c      
      logical :: combine_binary, combine_formatted
c
      integer, external ::  OMP_GET_THREAD_NUM  
      integer :: step, file_no, my_thread
      logical :: found_elements_no_values
c
c
      if( found_ele_stresses ) then
        write(termout,*) '   > Processing element stresses...'
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_elements_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, max_steps
         if( element_stresses_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_elements_no_values = .false.
         call do_patran_element_step( step, combine_binary,
     &                         combine_formatted, termout, 0, file_no,
     &                         found_elements_no_values )
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
c
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_elements_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, max_steps
         if( element_strains_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_elements_no_values = .false.
         call do_patran_element_step( step, combine_binary,
     &                         combine_formatted, termout, 1, file_no,
     &                         found_elements_no_values )
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
      end subroutine driver_combine_patran_element_results
c
c **********************************************************************
c *      contains:       driver_combine_flat_element_results           *
c **********************************************************************
c
      subroutine driver_combine_flat_element_results( combine_stream,
     &                                                combine_text )
      implicit none
c      
      logical :: combine_stream, combine_text
c
      integer, external ::  OMP_GET_THREAD_NUM  
      integer :: step, file_no, my_thread
      logical :: found_elements_no_values
c
c
      if( found_ele_stresses ) then
        write(termout,*) '   > Processing element stresses...'
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_elements_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, max_steps
         if( element_stresses_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_elements_no_values = .false.
         call do_flat_element_step( step, combine_stream,
     &                         combine_text, termout, 0, file_no,
     &                         found_elements_no_values )
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
        write(termout,*) '   > Processing element strains...'
c
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_elements_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, max_steps
         if( element_strains_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_elements_no_values = .false.
         call do_flat_element_step( step, combine_stream,
     &                         combine_text, termout, 1, file_no,
     &                         found_elements_no_values )
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
      end subroutine driver_combine_flat_element_results

c **********************************************************************
c *      contains:       driver_combine_flat_nodal_results             *
c **********************************************************************
c
      subroutine driver_combine_flat_nodal_results( combine_stream,
     &                                              combine_text )
      implicit none
c
      logical :: combine_stream, combine_text
c      
      integer, external ::  OMP_GET_THREAD_NUM  
      integer :: step, file_no, my_thread
      logical :: found_nodes_no_values
c
      if( found_node_stresses ) then
c
      write(termout,*) '   > Processing nodal stresses...'
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_nodes_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, max_steps
         if( nodal_stresses_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_nodes_no_values = .false.
         call do_flat_nodal_step( step, combine_stream,
     &                       combine_text, termout, 0, file_no,
     &                       found_nodes_no_values )
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
c
        write(termout, *) '   > Processing nodal strains...'
c$OMP PARALLEL
c$OMP& PRIVATE(step,my_thread,file_no,found_nodes_no_values)
c
c$OMP DO SCHEDULE(DYNAMIC)
        do step = 1, max_steps
         if( nodal_strains_step_list(step) .eq. 0 ) cycle
         my_thread = OMP_GET_THREAD_NUM()
         file_no = 10 + my_thread
         found_nodes_no_values = .false.
         call do_flat_nodal_step( step, combine_stream, 
     &                       combine_text, termout, 1, file_no,
     &                       found_nodes_no_values )
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

      end subroutine driver_combine_flat_nodal_results

c
      end subroutine driver

c **********************************************************************
c *                                                                    *
c *     support routine:      build_flat_file_names                    *
c *                                                                    *
c **********************************************************************
c
c
      subroutine build_flat_file_names( step, myid, node_stress_fname,
     &                         node_strain_fname, element_stress_fname,
     &                         element_strain_fname, termout,
     &                         combine_stream, combine_text,
     &                         file_name, sig_eps_type, data_type  )
      implicit none
c
      integer :: termout, step, myid, sig_eps_type, data_type
      logical :: combine_stream, combine_text
      character(len=*) :: node_stress_fname, node_strain_fname,
     &                    element_stress_fname, element_strain_fname
c
      integer :: fcount
      character(len=5) :: step_id
      character(len=4) :: proc_id
      character(len=6) :: format
      character(len=*) :: file_name
      logical :: local_debug
c
      local_debug              = .false.
      node_stress_fname(1:)    = ' '  
      node_strain_fname(1:)    = ' '
      element_stress_fname(1:) = ' '
      element_strain_fname(1:) = ' '
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
     &       '>> FATAL ERROR: routine build_text_file_names'
         write(termout,*) '>>   invalid format flag...'
         write(termout,*) '>>   job terminated...'
         stop
      end if
c
      node_stress_fname(1:)    = 'wn' // 's' // 
     &                           step_id // '_' //
     &                           format(1:fcount) // "." //proc_id     
      node_strain_fname(1:)    = 'wn' // 'e' // 
     &                           step_id // '_' //
     &                           format(1:fcount) // "." //proc_id     
c
      element_stress_fname(1:) = 'we' // 's' //
     &                           step_id // '_' // 
     &                           format(1:fcount) // "." //proc_id
      element_strain_fname(1:) = 'we' // 'e' //
     &                           step_id // '_' // 
     &                           format(1:fcount) // "." //proc_id
c
c          data_type = 0   (element results)
c                    = 1   (nodal results)
c
c          sig_eps_type = 0 (stress results)
c                       = 1 (strain results)  
c
      if( data_type .eq. 0 ) then
        if( sig_eps_type .eq. 1 ) then
           file_name(1:) = element_strain_fname(1:)
        elseif( sig_eps_type .eq. 0 ) then
           file_name(1:) = element_stress_fname(1:)
        else
           write(termout,9300)
           stop
        end if
      elseif( data_type .eq. 1 ) then 
        if( sig_eps_type .eq. 1 ) then
           file_name(1:) = node_strain_fname(1:)
        elseif( sig_eps_type .eq. 0 ) then
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
     &                      element_stress_fname, element_strain_fname
      end if
c 
      return      
c
 9000 format( i5.5 )
 9100 format( i4.4 )
 9200 format( 5x,' > file names: ', 4(/10x,a14) )
 9300 format( '>> FATAL ERROR: routine build_file_names...',
     &  /,    '                invalid sig_eps_type',
     &  /,    '                job terminated...' )
 9400 format( '>> FATAL ERROR: routine build_file_names...',
     &  /,    '                invalid data_type',
     &  /,    '                job terminated...' )

c
      end
c **********************************************************************
c *                                                                    *
c *          do_flat_element_step                                    *
c *                                                                    *
c **********************************************************************
c
c
      subroutine do_flat_element_step( step, combine_stream,
     &                                 combine_text, termout, 
     &                                 sig_eps_type, file_no,
     &                                 missing_flg )
c
      use size_parameter_data
      implicit none
c
c
c                parameter declarations
c                ----------------------
c
      integer :: step, termout, sig_eps_type, file_no
      logical :: combine_stream, combine_text, missing_flg
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
     &                     header_lines(max_header_lines)
      character(len=20) :: file_name
      logical :: local_debug, fatal, print_header, print
c
      double precision, allocatable,   dimension(:,:) :: element_values
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
     &                       element_strain_fname, termout,
     &                       combine_stream, combine_text,
     &                       file_name, sig_eps_type, 0 )
      if( local_debug ) write(termout,9700) step, file_name(1:20),
     &                                      file_no 
c
c                open the element result file for processor
c                zero. read the header records to get the 
c                number of data values per element and the number
c                of elements. each result file for other processors
c                must have the same values else bad.
c
      call open_old_flat_result_file( file_no, file_name,
     &                           combine_stream,
     &                           combine_text, termout, 1,
     &                           open_status, 0 )
c
      num_header_lines = 0
      if( combine_text ) call read_flat_header( file_no, header_lines,
     &                       num_header_lines, max_header_lines,
     &                       termout )
c
      if( local_debug ) 
     &     write(*,*) '.... num_header_lines: ', num_header_lines
c
c
c                 allocate space to hold the data values for all
c                 elements in the mode. keep track of which element
c                 values have been read from the files for the step.
c
      if( sig_eps_type .eq. 0 ) num_values = num_flat_elem_stresses
      if( sig_eps_type .eq. 1 ) num_values = num_flat_elem_strains
      num_elements = num_model_elems
      if( local_debug ) write(termout,9710) num_values 
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
      do proc_id = 0, max_procs
c
        call build_flat_file_names( step, proc_id, node_stress_fname,
     &                       node_strain_fname, element_stress_fname,
     &                       element_strain_fname, termout,
     &                       combine_stream, combine_text,
     &                       file_name, sig_eps_type, 0 )
        if( local_debug ) write(termout,9705) proc_id, file_name(1:20),
     &                                        file_no 
c
c                open the element result file for the processor.
c                read all data for a binary or a formmated file. if
c                the open fails, keep looking for processor files
c                until max_procs in case user blew away a
c                processor file.
c
        call open_old_flat_result_file(
     &           file_no, file_name, combine_stream,
     &           combine_text, termout, 2, open_status, 1 )
        if( open_status .ne. 0 ) cycle
c
c               skip element results header for this processor.
c
       if( combine_text ) call skip_flat_header( file_no )
c
c                read data values until end of file. stuff into the
c                single element results table. update count of
c                element appearances.
c
       call read_flat_element_values(
     &           file_no, element_values, num_values, num_elements,
     &           element_list, termout, combine_stream, 
     &           combine_text ) 
     &   
       close( unit=file_no,status='keep' )
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
        if(  element_list(element_id) .eq. 0 ) missing_flg = .true.
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
      if( combine_text ) file_name(14:) = " "  ! trim off the proc id
      if( combine_stream )  file_name(16:) = " " 
      call open_new_flat_result_file( file_no, file_name, 
     &             combine_stream, combine_text, termout, 1 )
c
      call write_flat_element_results( file_no, num_header_lines, 
     &      header_lines, num_values, num_elements, element_values,
     &      termout, combine_stream, combine_text, element_list )
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
 9700 format(5x,'> local debug in do_element_step: ',
     & /,10x,'step: ',i5,/,10x,'file_name: ',a20,
     & /,10x,'file_no: ',i5 )
 9705 format(5x,'> local debug in do_element_step: ',
     & /,10x,'proc_id: ',i5,/,10x,'file_name: ',a20,
     & /,10x,'file_no: ',i5 )
 9710 format(10x,'number of data values: ',i5)
 9720 format(5x,'> data arrays allocated. no. elements: ',i7)
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
     &        recl=3500, iostat= open_status )
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
c *          skip_flat_header                                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine skip_flat_header( file_no )
      implicit none
c
      integer :: file_no
c
      character :: first_char*1
c
      do 
        read(file_no,fmt="(a1)") first_char
        if( first_char .eq. "#"  .or.
     &      first_char .eq. "!"  .or.
     &      first_char .eq. "c"  .or.
     &      first_char .eq. "C"  ) cycle
        backspace(file_no)
        return
       end do
c
      return
      end
c
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
      integer :: file_no, termout, num_header_lines, max_header_lines
      character(len=80) :: header_lines(max_header_lines)
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
        backspace(file_no)
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
      integer :: read_status, element_id
      double precision, allocatable :: element_data(:)
c
      allocate( element_data(num_values) )
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
      return
c
 9130 format( '>> FATAL ERROR: routine read_flat_element_data...',
     &  /,    '                internal error. bad element #: ',i10,
     &  /,    '                job terminated...' )
 9330 format(i8,30e15.6)
 9200 format( '>> FATAL ERROR: invalid system state in ',
     &  /,    '                routine: read_element_data',
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
     &             combine_stream, combine_text, termout, 
     &             location )
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
      integer :: open_status
c
      if( combine_stream ) then
        open( unit=file_no, file=file_name, status='unknown', 
     &        access='stream', form='unformatted', iostat= open_status )
      elseif( combine_text ) then 
        open( unit=file_no, file=file_name, status='unknown', 
     &        access='sequential', form='formatted',
     &        recl=3500, iostat= open_status )
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
     &  /,    '                job terminated...' )
 9200 format( '>> FATAL ERROR: invalid system state @ ',i2,
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
     &                               termout, sig_eps_type, file_no, 
     &                               missing_flg )
c
      use size_parameter_data
      implicit none
c
c
c                parameter declarations
c                ----------------------
c
      integer :: step, termout, sig_eps_type, file_no
      logical :: combine_stream, combine_text, missing_flg
c
c                local variables
c                ---------------    
c
      integer:: proc_id, open_status, num_values,
     &          allocate_status, num_nodes,
     &          node_id, miss_count, num_header_lines
      integer, allocatable, dimension(:)  :: count_list
      integer, parameter :: max_header_lines = 50 
c
      logical ::  local_debug, print_header, print
c
      character(len=80) :: node_stress_fname, node_strain_fname,
     &               element_stress_fname, element_strain_fname,
     &               header_lines(max_header_lines)
      character(len=20) :: file_name
c
      double precision, allocatable,  dimension(:,:) :: node_values
      double precision, parameter :: zero = 0.0d00
c
c
      local_debug = .true.
c
c                 get name of root processor file for this load step.
c
      call build_flat_file_names( step, 0, node_stress_fname,
     &                       node_strain_fname, element_stress_fname,
     &                       element_strain_fname, termout,
     &                       combine_stream, combine_text,
     &                       file_name, sig_eps_type, 1 )
c
      if( local_debug ) write(termout,9700) step, file_name(1:20),
     &                                      file_no 
c
c                open the node result file for processor
c                zero. read/store header lines for text result
c
      call open_old_flat_result_file( file_no, file_name, 
     &                           combine_stream, combine_text, termout,
     &                           1, open_status, 0 )
      num_header_lines = 0
      if( combine_text ) call read_flat_header( file_no, header_lines,
     &                       num_header_lines, max_header_lines,
     &                       termout )
c
      if( local_debug ) 
     &     write(*,*) '.... num_header_lines: ', num_header_lines
c
c                 allocate space to hold the data values for all
c                 nodes in the model. keep track of which nodes
c                 have been read from the files for the step.
c
      if( sig_eps_type .eq. 0 ) num_values = num_flat_elem_stresses
      if( sig_eps_type .eq. 1 ) num_values = num_flat_elem_strains
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
      do proc_id = 0, max_procs ! MPI ranks
c
        call build_flat_file_names( step, proc_id, node_stress_fname,
     &                       node_strain_fname, element_stress_fname,
     &                       element_strain_fname, termout,
     &                       combine_stream, combine_text,
     &                       file_name, sig_eps_type, 1 )
c
c                open the node result file for the processor.
c                read all data for a stream or text file. if
c                the open fails, keep looking for processor files
c                until max_procs in case user blew away a
c                processor file.
c
        call open_old_flat_result_file( file_no, file_name, 
     &                combine_stream, combine_text, termout, 2, 
     &                open_status, 1 )
        if( open_status .ne. 0 ) cycle
        if( local_debug ) write(termout,9705) proc_id, file_name(1:20),
     &                                        file_no
c
c               skip header lines for text results files
c
       if( combine_text ) call skip_flat_header( file_no )
c
c                read node results for this processor. check
c                read data values until end of file. stuff into the
c                single node results table. update count of
c                node appearances.
c
       call read_flat_node_values( file_no, node_values, num_values,
     &                             num_nodes, count_list, termout, 
     &                             combine_stream, combine_text ) 
     &   
       close( unit=file_no,status='keep' )
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
      if( combine_text ) file_name(14:) = " "  ! trim off the proc id
      if( combine_stream )  file_name(16:) = " " 
      call open_new_flat_result_file( file_no, file_name, 
     &                     combine_stream, combine_text, termout, 2 )
c
      call write_flat_node_results( file_no, num_header_lines,
     &      header_lines, num_values, num_nodes, node_values,  
     &      count_list, termout, combine_stream, combine_text,
     &      sig_eps_type )
c    
      close(unit=file_no,status='keep')
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
      integer :: read_status, node_id, lcount
      double precision :: node_data(100)
c
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
c **********************************************************************
c *                                                                    *
c *     support routine:      build_patran_file_names                  *
c *                                                                    *
c **********************************************************************
c
c
      subroutine build_patran_file_names( step, myid, node_stress_fname,
     &                         node_strain_fname, element_stress_fname,
     &                         element_strain_fname, termout,
     &                         combine_binary, combine_formatted,
     &                         file_name, sig_eps_type, data_type  )
      implicit none
c
      integer :: termout, step, myid, sig_eps_type, data_type
      character(len=*) :: node_stress_fname, node_strain_fname,
     &                    element_stress_fname, element_strain_fname
      character(len=5) :: step_id
      character(len=4) :: proc_id
      character(len=1) :: format
      character(len=14) :: file_name
      logical ::  combine_binary, combine_formatted, local_debug
c
      local_debug              = .false.
      node_stress_fname(1:)    = ' '  
      node_strain_fname(1:)    = ' '
      element_stress_fname(1:) = ' '
      element_strain_fname(1:) = ' '
c
      write(step_id,9000) step
      write(proc_id,9100) myid
c
      if(  combine_binary ) then
         format = 'b'
      elseif ( combine_formatted ) then
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
c
c          data_type = 0   (element results)
c                    = 1   (nodal results)
c
c          sig_eps_type = 0 (stress results)
c                       = 1 (strain results)  
c
      if(  data_type .eq. 0 ) then
        if(  sig_eps_type .eq. 1 ) then
           file_name(1:) = element_strain_fname(1:)
        elseif( sig_eps_type .eq. 0 ) then
           file_name(1:) = element_stress_fname(1:)
        else
           write(termout,9300)
           stop
        end if
      elseif ( data_type .eq. 1 ) then 
        if(  sig_eps_type .eq. 1 ) then
           file_name(1:) = node_strain_fname(1:)
        elseif ( sig_eps_type .eq. 0 ) then
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
     &                      element_stress_fname, element_strain_fname
      end if
c 
      return      
c
 9000 format( i5.5 )
 9100 format( i4.4 )
 9200 format( 5x,' > file names: ', 4(/10x,a14) )
 9300 format( '>> FATAL ERROR: routine build_file_names...',
     &  /,    '                invalid sig_eps_type',
     &  /,    '                job terminated...' )
 9400 format( '>> FATAL ERROR: routine build_file_names...',
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
     &                          combine_formatted,
     &                          termout, sig_eps_type, file_no,
     &                          missing_flg )
c
      use size_parameter_data
      implicit none
c
c
c                parameter declarations
c                ----------------------
c
      integer :: step, termout, sig_eps_type, file_no
      logical :: combine_binary, combine_formatted, missing_flg
c
c                local variables
c                ---------------    
c
      integer:: myid, open_status, num_values,
     &          allocate_status, now_num_values, num_nodes,
     &          node_id, miss_count
      character(len=80) :: node_stress_fname, node_strain_fname,
     &               element_stress_fname, element_strain_fname
      character(len=14) :: file_name
      character(len=4) :: title1_binary(80), title2_binary(80)
      character(len=1) :: title1_formatted(80), title2_formatted(80)
      logical :: local_debug, print_header, print
c
      double precision, allocatable,  dimension(:,:) :: node_values
      integer, allocatable, dimension(:)  :: node_list
c   
      double precision, parameter :: zero = 0.0d00
c
c
      local_debug = .false.
c
c                 open the root processor file for this load step.
c                 get information about step then close down this
c                 file.
c
      call build_patran_file_names( step, 0, node_stress_fname,
     &                       node_strain_fname, element_stress_fname,
     &                       element_strain_fname, termout,
     &                       combine_binary, combine_formatted,
     &                       file_name, sig_eps_type, 1 )
c
      if(  local_debug ) write(termout,9700) step, file_name(1:14),
     &                   file_no 
c
c                open the node result file for processor
c                zero. read the header records to get the 
c                number of data values per node and the number
c                of nodes. each result file for other processors
c                must have the same values else bad.
c
      call open_old_patran_result_file( file_no, file_name, 
     &                           combine_binary,
     &                           combine_formatted, termout, 1,
     &                           open_status, 0 )
c
      call read_patran_node_file_header( 
     &        file_no, title1_binary,
     &        num_values, title2_binary, num_nodes, 
     &        title1_formatted, title2_formatted, termout,
     &        combine_binary, combine_formatted )
      close( unit=file_no,status='keep' )
      if(  local_debug ) write(termout,9710) num_values, num_nodes
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
      do myid = 0, max_procs ! MPI ranks
c
        call build_patran_file_names( step, myid, node_stress_fname,
     &                       node_strain_fname, element_stress_fname,
     &                       element_strain_fname, termout,
     &                       combine_binary, combine_formatted,
     &                       file_name, sig_eps_type, 1 )
c
c                open the node result file for the processor.
c                read all data for a binary or a formmated file. if
c                the open fails, keep looking for processor files
c                until max_procs in case user blew away a
c                processor file.
c
        call open_old_patran_result_file(
     &           file_no, file_name, combine_binary,
     &           combine_formatted, termout, 2, open_status, 1 )
        if( open_status .ne. 0 ) cycle
        if( local_debug ) write(termout,9705) myid, file_name(1:14),
     &                                         file_no
c
c               read node results for this processor. check
c               header info for consistency with other processor
c               files for this step.
c
        call read_patran_node_file_header_dummy( 
     &        file_no, now_num_values, termout,
     &        combine_binary, combine_formatted )
        if( num_values .ne. now_num_values ) then
         write(termout,9120)
         stop
        end if
c
c                read data values until end of file. stuff into the
c                single node results table. update count of
c                node appearances.
c
       call read_patran_node_values(
     &           file_no, node_values, num_values, num_nodes,
     &           node_list, termout, combine_binary, 
     &           combine_formatted ) 
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
c
      call open_new_patran_result_file( 
     &             file_no, file_name, 
     &             combine_binary, combine_formatted, termout, 2 )
c
      call write_patran_node_results( 
     &         file_no, title1_binary, num_values, title2_binary,
     &         num_nodes, node_values, node_list,  title1_formatted,
     &         title2_formatted, termout, combine_binary, 
     &         combine_formatted, sig_eps_type )
c    
      close(unit=file_no,status='keep')
c
      deallocate( node_values, node_list, stat=allocate_status )
      call check_deallocate_status( termout, 20, allocate_status )
c
      return    
c
 9120 format( '>> FATAL ERROR: routine do_nodal_step...',
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
     &                            combine_formatted,
     &                            termout, sig_eps_type, file_no,
     &                            missing_flg )
c
      use size_parameter_data
      implicit none
c
c
c                parameter declarations
c                ----------------------
c
      integer :: step, termout, sig_eps_type, file_no
      logical :: combine_binary, combine_formatted, missing_flg
c
c                local variables
c                ---------------    
c
      integer :: myid, open_status, num_values,
     &           allocate_status, now_num_values, num_elements,
     &           element_id, now_num_elements, miss_count
      character(len=80) :: node_stress_fname, node_strain_fname,
     &               element_stress_fname, element_strain_fname
      character(len=14) :: file_name
      character(len=4) :: title1_binary(80), title2_binary(80)
      character(len=1) :: title1_formatted(80), title2_formatted(80)
      logical :: local_debug, fatal, print_header, print
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
     &                       node_strain_fname, element_stress_fname,
     &                       element_strain_fname, termout,
     &                       combine_binary, combine_formatted,
     &                       file_name, sig_eps_type, 0 )
      if(  local_debug ) write(termout,9700) step, file_name(1:14),
     &                   file_no 
c
c                open the element result file for processor
c                zero. read the header records to get the 
c                number of data values per element and the number
c                of elements. each result file for other processors
c                must have the same values else bad.
c
      call open_old_patran_result_file( file_no, file_name,
     &                           combine_binary,
     &                           combine_formatted, termout, 1,
     &                           open_status, 0 )
c
      call read_patran_element_file_header( 
     &        file_no, title1_binary,
     &        num_values, title2_binary, num_elements, 
     &        title1_formatted, title2_formatted, termout,
     &        combine_binary, combine_formatted )
      close( unit=file_no,status='keep' )
      if(  local_debug ) write(termout,9710) num_values 
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
      if(  local_debug ) write(termout,9720) num_elements
c
c                loop over all the separate processor files of 
c                element results for this load step.
c                open the file and read the
c                element results into the single array here. keep
c                track if an element appears other than once... bad.
c
      do myid = 0, max_procs
c
        call build_patran_file_names( step, myid, node_stress_fname,
     &                       node_strain_fname, element_stress_fname,
     &                       element_strain_fname, termout,
     &                       combine_binary, combine_formatted,
     &                       file_name, sig_eps_type, 0 )
        if(  local_debug ) write(termout,9705) myid, file_name(1:14),
     &                                         file_no 
c
c                open the element result file for the processor.
c                read all data for a binary or a formmated file. if
c                the open fails, keep looking for processor files
c                until max_procs in case user blew away a
c                processor file.
c
        call open_old_patran_result_file(
     &           file_no, file_name, combine_binary,
     &           combine_formatted, termout, 2, open_status, 1 )
        if(  open_status .ne. 0 ) cycle
c
c               read element results for this processor. check
c               header info for consistency with other processor
c               files for this step.
c
        call read_patran_element_file_header_dummy( 
     &        file_no, now_num_values, termout,
     &        combine_binary, combine_formatted, now_num_elements )
        if(  num_values .ne. now_num_values .or.
     &       num_elements .ne. now_num_elements ) then
         write(termout,9120)
         stop
       end if   
c
c                read data values until end of file. stuff into the
c                single element results table. update count of
c                element appearances.
c
       call read_patran_element_values(
     &           file_no, element_values, num_values, num_elements,
     &           element_list, termout, combine_binary, 
     &           combine_formatted ) 
     &   
       close( unit=file_no,status='keep' )
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
      print = .false.
c
      do element_id = 1, num_elements
         if(  element_list(element_id) .eq. 0 ) missing_flg = .true.
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
c               single patran results file for this step.
c               write in all the data.
c
c
      call open_new_patran_result_file( 
     &             file_no, file_name, 
     &             combine_binary, combine_formatted, termout, 1 )
c
      call write_patran_element_results( 
     &         file_no, title1_binary, num_values, title2_binary,
     &         num_elements, element_values, title1_formatted,
     &         title2_formatted, termout, combine_binary, 
     &         combine_formatted, element_list )
c    
      close(unit=file_no,status='keep')
c
      deallocate( element_values, element_list, stat=allocate_status)
      call check_deallocate_status( termout, 1, allocate_status )
c
      return    
c
 9120 format( '>> FATAL ERROR: routine do_element_step...',
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
 9700 format(5x,'> local debug in do_element_step: ',
     & /,10x,'step: ',i5,/,10x,'file_name: ',a14,
     & /,10x,'file_no: ',i5 )
 9705 format(5x,'> local debug in do_element_step: ',
     & /,10x,'myid: ',i5,/,10x,'file_name: ',a14,
     & /,10x,'file_no: ',i5 )
 9710 format(10x,'number of data values: ',i5)
 9720 format(5x,'> data arrays allocated. no. elements: ',i7)
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
     &        recl=350, iostat= open_status )
      elseif( combine_formatted ) then 
        open( unit=file_no, file=file_name, status='old', 
     &        access='sequential', form='formatted',
     &        recl=350, iostat= open_status )
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
         read(unit=file_no,iostat=read_status,
     &        fmt=9320) num_elements
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
         read(unit=file_no,iostat=read_status,
     &        fmt=9320) now_num_elements
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
          write(termout,9130)
          stop
        end if
c
      end do
c
      return
c
 9130 format( '>> FATAL ERROR: routine read_element_data...',
     &  /,    '                internal error @ 5',
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
      integer :: open_status
c
      if( combine_binary ) then
        open( unit=file_no, file=file_name(1:9), status='unknown', 
     &        access='sequential', form='unformatted',
     &        recl=350, iostat= open_status )
      elseif( combine_formatted ) then 
        open( unit=file_no, file=file_name(1:9), status='unknown', 
     &        access='sequential', form='formatted',
     &        recl=350, iostat= open_status )
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
     &  /,    '                job terminated...' )
 9200 format( '>> FATAL ERROR: invalid system state @ ',i2,
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
     &  num_values, termout, combine_binary, combine_formatted )
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
c **********************************************************************
c *                                                                    *
c *          check_read_status                                         *
c *                                                                    *
c **********************************************************************
c
c
      subroutine check_read_status( file_no, read_status, termout, 
     &                              location )
      implicit none
c
      integer :: file_no, read_status, termout, location
c
      if( read_status .eq. 0 ) return
      write(termout,*) '>> FATAL ERROR: result file read'
      write(termout,*) '                failed at location: ',
     &                  location
      write(termout,*) '>> Job terminated...'
      stop
      end
c **********************************************************************
c *                                                                    *
c *          check_write_status                                        *
c *                                                                    *
c **********************************************************************
c
c
      subroutine check_write_status( file_no, write_status, termout,
     &                              location )
      implicit none
c
      integer :: file_no, write_status, termout, location
c
      if( write_status .eq. 0 ) return
      write(termout,*) '>> FATAL ERROR: result file write'
      write(termout,*) '                failed at location: ',
     &                  location
      write(termout,*) '                status code: ', write_status
      write(termout,*) '>> Job terminated...'
      stop
      end
c **********************************************************************
c *                                                                    *
c *          check_allocate_status                                     *
c *                                                                    *
c **********************************************************************
c
c
      subroutine check_allocate_status( termout, location, 
     &                                  allocate_status )
      implicit none
c
      integer ::  termout, location, allocate_status
c
      if( allocate_status .eq. 0 ) return
      write(termout,*) '>> FATAL ERROR: allocation failure'
      write(termout,*) '                at location: ',
     &                  location
      write(termout,*) '>> Job terminated...'
      stop
      end
c **********************************************************************
c *                                                                    *
c *          check_deallocate_status                                   *
c *                                                                    *
c **********************************************************************
c
c
      subroutine check_deallocate_status( termout, location, 
     &                                  allocate_status )
      implicit none
c
      integer ::  termout, location, allocate_status
c
      if( allocate_status .eq. 0 ) return
      write(termout,*) '>> FATAL ERROR: deallocation failure'
      write(termout,*) '                at location: ',
     &                  location
      write(termout,*) '>> Job terminated...'
      stop
      end
c     ****************************************************************
c     *                                                              *
c     *           subroutine compute_extra_sig_eps_values            *
c     *                                                              *
c     ****************************************************************
c
      subroutine compute_extra_sig_eps_values( node_values, num_nodes,
     &                                         sig_eps_type )
      implicit none
c
      integer :: num_nodes, sig_eps_type
      double precision :: node_values(num_nodes,*)
      logical :: stress 
c
      stress = sig_eps_type .eq. 0
c
      if( stress ) then
         call princ_inv_stress( node_values, num_nodes )
         call princ_stress( node_values, num_nodes )
      else
         call princ_inv_strain( node_values, num_nodes )
         call princ_strain( node_values, num_nodes )
      end if
c
      return
      end
    
c     ****************************************************************
c     *                                                              *
c     *                  subroutine princ_inv_stress                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/18/2017 rhd              *
c     *                                                              *
c     *     this subroutine computes the mises stress, principal     *
c     *     invariants at the nodes/elements                         * 
c     *                                                              *
c     ****************************************************************
c
      subroutine princ_inv_stress( results, num_nodes )
      implicit none
c
      integer :: num_nodes
      double precision :: results(num_nodes,*)
c
      integer :: i
      double precision, parameter :: six=6.0d0, iroot2=0.70710678118d0
c        
      do i = 1, num_nodes
         results(i,8) = sqrt( (results(i,1)-results(i,2))**2+
     &                 (results(i,2)-results(i,3))**2+
     &                 (results(i,1)-results(i,3))**2+
     &             six*(results(i,4)**2+results(i,5)**2+
     &                  results(i,6)**2) )*iroot2
         results(i,12) = results(i,1) +  results(i,2) +
     &                        results(i,3)
         results(i,13) = results(i,1) * results(i,2) + 
     &                        results(i,2) * results(i,3) +
     &                        results(i,1) * results(i,3) -
     &                        results(i,4) * results(i,4) -
     &                        results(i,5) * results(i,5) -
     &                        results(i,6) * results(i,6)
         results(i,14) = results(i,1) * 
     &            ( results(i,2) * results(i,3) -
     &              results(i,5) * results(i,5) ) -
     &                        results(i,4) *
     &            ( results(i,4) * results(i,3) -
     &              results(i,5) * results(i,6) ) +
     &                        results(i,6) *
     &            ( results(i,4) * results(i,5) -
     &              results(i,2) * results(i,6) ) 
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine princ_inv_strain                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/29/2017 rhd              *
c     *                                                              *
c     *     this subroutine computes equivalent strain and           *
c     *     principal invariants  at the nodes/elem                  * 
c     *                                                              *
c     ****************************************************************
c
      subroutine princ_inv_strain( results, num_nodes )
      implicit none
c
      integer :: num_nodes
      double precision :: results(num_nodes,*)
c
      integer :: i
      double precision :: t1, t2, t3
      double precision, parameter ::  half = 0.5d0,
     &                                root23 = 0.47140452079d0,
     &                                onep5 = 1.5d0
c   
      do i = 1, num_nodes
         t1 =  half * results(i,4)
         t2 =  half * results(i,5)
         t3 =  half * results(i,6)
         results(i,7) = root23 * sqrt( 
     &     ( results(i,1) - results(i,2) ) ** 2 +
     &     ( results(i,2) - results(i,3) ) ** 2 +
     &     ( results(i,1) - results(i,3) ) ** 2 +
     &         onep5 * ( results(i,4)**2 + results(i,5)**2+
     &         results(i,6)**2 ) )
         results(i,8) = results(i,1) +  results(i,2) +
     &                        results(i,3)
         results(i,9) = - t1 * t1 - t2 * t2 - t3 * t3 +
     &                        results(i,1) * results(i,2) +
     &                        results(i,2) * results(i,3) +
     &                        results(i,1) * results(i,3)
         results(i,10) =
     &     results(i,1) * ( results(i,2) * results(i,3) -
     &              t2 * t2 ) -  t1 * ( t1 * results(i,3) -
     &              t2 * t3 ) + t3 * ( t1 * t2 -
     &              results(i,2) * t3 ) 
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine princ_strain                 *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 1/18/2016 rhd              *
c     *                                                              *
c     *     this subroutine computes principal strain values         *
c     *     at the nodes/elem after the primary values have          *
c     *     been averaged                                            * 
c     *                                                              *
c     ****************************************************************
c
      subroutine princ_strain( results, num_nodes )
      implicit none
c
      integer :: num_nodes   
      double precision ::  results(num_nodes,*)
c
c                    locals
c
      integer :: i, ier
      integer, parameter :: nstr=6,ndim=3
      double precision ::
     &  temp_strain(nstr), wk(ndim), ev(nstr), evec(ndim,ndim)
      double precision, parameter :: half = 0.5d00
c
c
c        calculate the principal strains and there direction cosines by
c        passing off strains into dummy array, then calling an eigenvalue
c        eigenvector routine, and then passing the values from this routine
c        back to the appropriate places within results
c      
       do i = 1, num_nodes
          temp_strain(1) = results(i,1)
          temp_strain(2) = results(i,4) * half
          temp_strain(3) = results(i,2)
          temp_strain(4) = results(i,6) * half
          temp_strain(5) = results(i,5) * half
          temp_strain(6) = results(i,3)
          call ou3dpr( temp_strain, ndim, 1, ev, evec, ndim, wk, ier )   
          results(i,11) = ev(1)     
          results(i,12) = ev(2)
          results(i,13) = ev(3)
          results(i,14) = evec(1,1)
          results(i,15) = evec(2,1)
          results(i,16) = evec(3,1)
          results(i,17) = evec(1,2)
          results(i,18) = evec(2,2)
          results(i,19) = evec(3,2)
          results(i,20) = evec(1,3)
          results(i,21) = evec(2,3)
          results(i,22) = evec(3,3)
      end do
      return
      end  
c     ****************************************************************
c     *                                                              *
c     *                      subroutine princ_stress                 *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 10/04/94                   *
c     *                                   03/05/95 kck               *        
c     *                                   06/10/97 rhd               *        
c     *                                                              *
c     *     this subroutine computes principal stress values         *
c     *     at the nodes/elem after the primary values have          *
c     *     been averaged                                            * 
c     *                                                              *
c     ****************************************************************
c
      subroutine princ_stress( results, num_nodes )
      implicit none
c
      integer :: num_nodes   
      double precision ::  results(num_nodes,*)
c
c                    locals
c
      integer :: i, ier
      integer, parameter :: nstr=6, ndim=3
      double precision :: temp_stress(nstr), wk(ndim), ev(nstr),
     &                    evec(ndim,ndim)
c
c
c        calculate the principal stresses and there direction cosines by
c        passing off stresses into dummy array, then calling an eigenvalue
c        eigenvector routine, and then passing the values from this routine
c        back to the appropriate places within results
c      
       do i = 1, num_nodes
          temp_stress(1) = results(i,1)
          temp_stress(2) = results(i,4)
          temp_stress(3) = results(i,2)
          temp_stress(4) = results(i,6)
          temp_stress(5) = results(i,5)
          temp_stress(6) = results(i,3)
          call ou3dpr( temp_stress, ndim, 1, ev, evec, ndim, wk, ier )   
          results(i,15) = ev(1)     
          results(i,16) = ev(2)
          results(i,17) = ev(3)
          results(i,18) = evec(1,1)
          results(i,19) = evec(2,1)
          results(i,20) = evec(3,1)
          results(i,21) = evec(1,2)
          results(i,22) = evec(2,2)
          results(i,23) = evec(3,2)
          results(i,24) = evec(1,3)
          results(i,25) = evec(2,3)
          results(i,26) = evec(3,3)
      end do
      return
      end  
c *******************************************************************
c *                                                                 *
c *      element output service routine -- ou3dpr                   *
c *                                                                 *
c *******************************************************************
c
c
      subroutine ou3dpr( a,n,ijob,d,z,iz,wk,ier )
      implicit double precision (a-h,o-z)
c
c   purpose             - eigenvalues and (optionally) eigenvectors of
c                           a real symmetric matrix in symmetric
c                           storage mode
c
c   arguments    a      - input real symmetric matrix of order n,
c                           stored in symmetric storage mode,
c                           whose eigenvalues and eigenvectors
c                           are to be computed. input a is
c                           destroyed if ijob is equal to 0 or 1.
c                n      - input order of the matrix a.
c                ijob   - input option parameter, when
c                           ijob = 0, compute eigenvalues only
c                           ijob = 1, compute eigenvalues and eigen-
c                             vectors.
c                           ijob = 2, compute eigenvalues, eigenvectors
c                             and performance index.
c                           ijob = 3, compute performance index only.
c                           if the performance index is computed, it is
c                           returned in wk(1). the routines have
c                           performed (well, satisfactorily, poorly) if
c                           wk(1) is (less than 1, between 1 and 100,
c                           greater than 100).
c                d      - output vector of length n,
c                           containing the eigenvalues of a.
c                z      - output n by n matrix containing
c                           the eigenvectors of a.
c                           the eigenvector in column j of z corres-
c                           ponds to the eigenvalue d(j).
c                           if ijob = 0, z is not used.
c                iz     - input row dimension of matrix z exactly as
c                           specified in the dimension statement in the
c                           calling program.
c                wk     - work area, the length of wk depends
c                           on the value of ijob, when
c                           ijob = 0, the length of wk is at least n.
c                           ijob = 1, the length of wk is at least n.
c                           ijob = 2, the length of wk is at least
c                             n(n+1)/2+n.
c                           ijob = 3, the length of wk is at least 1.
c                ier    - error parameter (output)
c                         terminal error
c                           ier = 128+j, indicates that ourt2s failed
c                             to converge on eigenvalue j. eigenvalues
c                             and eigenvectors 1,...,j-1 have been
c                             computed correctly, but the eigenvalues
c                             are unordered. the performance index
c                             is set to 1000.0
c                         warning error (with fix)
c                           ier = 66, indicates ijob is less than 0 or
c                             ijob is greater than 3. ijob set to 1.
c                           ier = 67, indicates ijob is not equal to
c                             zero, and iz is less than the order of
c                             matrix a. ijob is set to zero.
c
c                                  specifications for arguments
      integer            n,ijob,iz,ier
      double precision
     &                   a(*),d(*),wk(*),z(iz,*)
c                                  specifications for local variables
      integer            jer,na,nd,iiz,ibeg,il,kk,lk,i,j,k,l
      double precision
     &                   anorm,asum,pi,sumz,sumr,an,s,ten,rdelp,zero,
     &                   one,thous
      data               rdelp /0.745058d-08/
      data               zero,one/0.d0,1.d0/,ten/10.d0/,thous/1000.d0/
c
c                                  initialize error parameters
c                                  first executable statement
      ier = 0
      jer = 0
      if( ijob .ge. 0 .and. ijob .le. 3) go to 5
c
c                                  warning error - ijob is not in the
c                                    range
      ier = 66
      ijob = 1
      go to 10
    5 if( ijob .eq. 0) go to 20
   10 if( iz .ge. n) go to 15
c                                  warning error - iz is less than n
c                                    eigenvectors can not be computed,
c                                    ijob set to zero
      ier = 67
      ijob = 0
   15 if( ijob .eq. 3) go to 65
   20 na = (n*(n+1))/2
      if( ijob .ne. 2) go to 35
      do 30 i=1,na
         wk(i) = a(i)
   30 continue
c                                  save input a if ijob = 2
   35 nd = 1
      if( ijob .eq. 2) nd = na+1
c                                 wnfe00001.0000 reduce a to symmetric tridiagonal
c                                    form
      call ouhous (a,n,d,wk(nd),wk(nd))
      iiz = 1
      if( ijob .eq. 0) go to 50
      iiz = iz
c                                  set z to the identity matrix
      do 45 i=1,n
         do 40 j=1,n
            z(i,j) = zero
   40    continue
         z(i,i) = one
   45 continue
c                                  compute eigenvalues and eigenvectors
   50 call ourt2s (d,wk(nd),n,z,iiz,jer)
      if( ijob .eq. 0) go to 9000
      if( jer .gt. 128) go to 55
c                                  back transform eigenvectors
      call ouobks (a,n,1,n,z,iz)
   55 if( ijob .le. 1) go to 9000
c                                  move input matrix back to a
      do 60 i=1,na
         a(i) = wk(i)
   60 continue
      wk(1) = thous
      if( jer .ne. 0) go to 9000
c                                  compute 1 - norm of a
   65 anorm = zero
      ibeg = 1
      do 75 i=1,n
         asum = zero
         il = ibeg
         kk = 1
         do 70 l=1,n
            asum =asum+abs(a(il))
            if( l .ge. i) kk = l
            il = il+kk
   70    continue
         anorm = dmax1(anorm,asum)
         ibeg = ibeg+i
   75 continue
      if( anorm .eq. zero) anorm = one
c                                  compute performance index
      pi = zero
      do 90 i=1,n
         ibeg = 1
         s = zero
         sumz = zero
         do 85 l=1,n
            lk = ibeg
            kk = 1
            sumz = sumz+abs(z(l,i))
            sumr = -d(i)*z(l,i)
            do 80 k=1,n
               sumr = sumr+a(lk)*z(k,i)
               if( k .ge. l) kk = k
               lk = lk+kk
   80       continue
            s = s+abs(sumr)
            ibeg = ibeg+l
   85    continue
         if( sumz .eq. zero) go to 90
         pi = dmax1(pi,s/sumz)
   90 continue
      an = n
      pi = pi/(anorm*ten*an*rdelp)
      wk(1) = pi
 9000 continue
      if(  ier .ne. 0 ) return
      if( jer .eq. 0) go to 9005
      ier = jer
      return
 9005 return
      end
c
c
c *******************************************************************
c *                                                                 *
c *      element output service routine -- ouobks                   *
c *                                                                 *
c *******************************************************************
c
c
      subroutine ouobks( a,n,m1,m2,z,iz )
      implicit double precision (a-h,o-z)
c
c   purpose             - back transformation to form the eigenvectors
c                           of the original symmetric matrix from the
c                           eigenvectors of the tridiagonal matrix
c
c   arguments    a      - the array contains the details of the house-
c                           holder reduction of the original matrix a
c                           as generated by routine ouhous.(input)
c                n      - order of the real symmetric matrix.(input)
c                m1     - m1 and m2 are two input scalars such that
c                           eigenvectors m1 to m2 of the tridiagonal
c                           matrix a have been found and normalized
c                           according to the euclidean norm.
c                m2     - see above - m1
c                z      - a two dimensional array of size n x (m2-m1+1)
c                           which contains eigenvectors m1 to m2 of
c                           tridiagonal matrix t, normalized according
c                           to euclidean norm. input z can be produced
c                           by routine ourt2s, the resultant
c                           matrix overwrites the input z.(input/output)
c                iz     - row dimension of matrix z exactly as
c                           specified in the dimension statement in the
c                           calling program. (input)
c
c
      double precision
     &                   zero
      dimension          a(*),z(iz,*)
      data zero /0.0/

c                                  first executable statement
      if( n.eq.1) go to 30
      do 25 i=2,n
         l = i-1
         ia = (i*l)/2
         h = a(ia+i)
         if( h.eq.zero) go to 25
c                                  derives eigenvectors m1 to m2 of
c                                  the original matrix from eigenvectors
c                                  m1 to m2 of the symmetric
c                                  tridiagonal matrix
         do 20 j = m1,m2
            s = zero
            do 10 k = 1,l
               s = s+a(ia+k)*z(k,j)
   10       continue
            s = s/h
            do 15 k=1,l
               z(k,j) = z(k,j)-s*a(ia+k)
   15       continue
   20    continue
   25 continue
   30 return
      end
c *******************************************************************
c *                                                                 *
c *      element output service routine -- ouhous                   *
c *                                                                 *
c *******************************************************************
c
c
      subroutine ouhous( a,n,d,e,e2 )
      implicit double precision (a-h,o-z)
c
c   purpose             - reduction of a symmetric matrix to symmetric
c                           tridiagonal form using a householder
c                           reduction
c
c   arguments    a      - the given n x n, real symmetric matrix a,
c                           where a is stored in symmetric storage mode.
c                           the input a is replaced by the details of
c                           the householder reduction of a.
c                n      - input order of a and the length of d, e, and
c                           e2.
c                d      - the output array of length n, giving the
c                           diagonal elements of the tridiagonal matrix.
c                e      - the output array of length n, giving the sub-
c                           diagonal in the last (n-1) elements, e(1) is
c                           set to zero.
c                e2     - output array of length n.  e2(i) = e(i)**2.
c
c
      dimension          a(*),d(n),e(n),e2(n)
      double precision
     &                   a,d,e,e2,zero,h,scale,one,scale1,f,g,hh
      data               zero/0.d0/,one/1.d0/
c                                  first executable statement
      np1 = n+1
      nn = (n*np1)/2-1
      nbeg = nn+1-n
      do 70 ii = 1,n
         i = np1-ii
         l = i-1
         h = zero
         scale = zero
         if( l .lt. 1) go to 10
c                                  scale row (algol tol then not needed)
         nk = nn
         do 5 k = 1,l
            scale = scale+abs(a(nk))
            nk = nk-1
    5    continue
         if( scale .ne. zero) go to 15
   10    e(i) = zero
         e2(i) = zero
         go to 65
   15    nk = nn
         scale1 = one/scale
         do 20 k = 1,l
            a(nk) = a(nk)*scale1
            h = h+a(nk)*a(nk)
            nk = nk-1
   20    continue
         e2(i) = scale*scale*h
         f = a(nn)
         g = -sign(sqrt(h),f)
         e(i) = scale*g
         h = h-f*g
         a(nn) = f-g
         if( l .eq. 1) go to 55
         f = zero
         jk1 = 1
         do 40 j = 1,l
            g = zero
            ik = nbeg+1
            jk = jk1
c                                  form element of a*u
            do 25 k = 1,j
               g = g+a(jk)*a(ik)
               jk = jk+1
               ik = ik+1
   25       continue
            jp1 = j+1
            if( l .lt. jp1) go to 35
            jk = jk+j-1
            do 30 k = jp1,l
               g = g+a(jk)*a(ik)
               jk = jk+k
               ik = ik+1
   30       continue
c                                  form element of p
   35       e(j) = g/h
            f = f+e(j)*a(nbeg+j)
            jk1 = jk1+j
   40    continue
         hh = f/(h+h)
c                                  form reduced a
         jk = 1
         do 50 j = 1,l
            f = a(nbeg+j)
            g = e(j)-hh*f
            e(j) = g
            do k = 1,j
               a(jk) = a(jk)-f*e(k)-g*a(nbeg+k)
               jk = jk+1
            enddo
   50    continue
   55    do 60 k = 1,l
            a(nbeg+k) = scale*a(nbeg+k)
   60    continue
   65    d(i) = a(nbeg+i)
         a(nbeg+i) = h*scale*scale
         nbeg = nbeg-i+1
         nn = nn-i
   70 continue
      return
      end
c *******************************************************************
c *                                                                 *
c *      element output service routine -- ourt2s
c *                                                                 *
c *******************************************************************
c
c
      subroutine ourt2s(d,e,n,z,iz,ier)
      implicit double precision (a-h,o-z)
c   purpose             - eigenvalues and (optionally) eigenvectors of
c                           a symmetric tridiagonal matrix using the
c                           ql method.
c
c   arguments    d      - on input, the vector d of length n contains
c                           the diagonal elements of the symmetric
c                           tridiagonal matrix t.
c                           on output, d contains the eigenvalues of
c                           t in ascending order.
c                e      - on input, the vector e of length n contains
c                           the sub-diagonal elements of t in position
c                           2,...,n. on output, e is destroyed.
c                n      - order of tridiagonal matrix t.(input)
c                z      - on input, z contains the identity matrix of
c                           order n.
c                           on output, z contains the eigenvectors
c                           of t. the eigenvector in column j of z
c                           corresponds to the eigenvalue d(j).
c                iz     - input row dimension of matrix z exactly as
c                           specified in the dimension statement in the
c                           calling program. if iz is less than n, the
c                           eigenvectors are not computed. in this case
c                           z is not used.
c                ier    - error parameter
c                         terminal error
c                           ier = 128+j, indicates that ourt2s failed
c                             to converge on eigenvalue j. eigenvalues
c                             and eigenvectors 1,...,j-1 have been
c                             computed correctly, but the eigenvalues
c                             are unordered.
c
c
      dimension          d(*),e(*),z(iz,*)
      data               rdelp /0.745058e-08/
      data               zero,one/0.0,1.0/
c                                  move the last n-1 elements
c                                  of e into the first n-1 locations
c                                  first executable statement
      ier  = 0
      if( n .eq. 1) go to 9005
      do 5  i=2,n
         e(i-1) = e(i)
    5 continue
      e(n) = zero
      b = zero
      f = zero
      do  60  l=1,n
         j = 0
         h = rdelp*(abs(d(l))+abs(e(l)))
         if( b.lt.h) b = h
c                                  look for small sub-diagonal element
         do 10  m=l,n
            k=m
            if( abs(e(k)) .le. b) go to 15
   10    continue
   15    m = k
         if( m.eq.l) go to 55
   20    if( j .eq. 30) go to 85
         j = j+1
         l1 = l+1
         g = d(l)
         p = (d(l1)-g)/(e(l)+e(l))
         r = sqrt(p*p+one)
         d(l) = e(l)/(p+sign(r,p))
         h = g-d(l)
         do 25 i = l1,n
            d(i) = d(i)-h
   25    continue
         f = f+h
c                                  ql transformation
         p = d(m)
         c = one
         s = zero
         mm1 = m-1
         mm1pl = mm1+l
         if( l.gt.mm1) go to 50
         do 45 ii=l,mm1
            i = mm1pl-ii
            G = C*E(I)
            H = C*P
            if( abs(p).lt.abs(e(i))) go to 30
            c = e(i)/p
            r = sqrt(c*c+one)
            e(i+1) = s*p*r
            s = c/r
            c = one/r
            go to 35
   30       c = p/e(i)
            r = sqrt(c*c+one)
            e(i+1) = s*e(i)*r
            s = one/r
            c = c*s
   35       p = c*d(i)-s*g
            d(i+1) = h+s*(c*g+s*d(i))
            if( iz .lt. n) go to 45
c                                  form vector
            do 40 k=1,n
               h = z(k,i+1)
               z(k,i+1) = s*z(k,i)+c*h
               z(k,i) = c*z(k,i)-s*h
   40       continue
   45    continue
   50    e(l) = s*p
         d(l) = c*p
         if(  abs(e(l)) .gt.b) go to 20
   55    d(l) = d(l) + f
   60 continue
c                                  order eigenvalues and eigenvectors
      do  80  i=1,n
         k = i
         p = d(i)
         ip1 = i+1
         if( ip1.gt.n) go to 70
         do 65  j=ip1,n
            if( d(j) .ge. p) go to 65
            k = j
            p = d(j)
   65    continue
   70    if( k.eq.i) go to 80
         d(k) = d(i)
         d(i) = p
         if( iz .lt. n) go to 80
         do 75 j = 1,n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
   75    continue
   80 continue
      go to 9005
   85 ier = 128+l
 9000 continue
 9005 return
      end

