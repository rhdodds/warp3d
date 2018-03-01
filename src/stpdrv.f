c     ****************************************************************
c     *                                                              *
c     *                      subroutine stpdrv                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 2/8/2018 rhd               *
c     *                                                              *
c     *     drive the solution process list of steps specified       *
c     *     on the compute command. compute intermediate             *
c     *     steps as required to satisfy request.                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine stpdrv( stplst, nsteps, ldnum )
      use global_data ! old common.main
c
      use main_data, only : cnstrn, cnstrn_in, temp_nodmap, temp_nodlod,
     &                      output_packets, run_user_solution_routine,
     &                      output_command_file, max_step_limit,
     &                      output_step_bitmap_list, stpchk,
     &                      user_cnstrn_stp_factors
c
      use damage_data, only : growth_by_kill, growth_by_release
c
      implicit none
c
      integer :: stplst(*), nsteps, ldnum
c
c              declare all locals here - visible in contains
c
      integer :: icn, iplist, diff, istep, step
      real :: dumr, con_stp_factor
      character :: dums
      logical :: mf_ratio_change, stpdrv_error
      double precision :: mf, mf_nm1, dumd, load_reduce_fact
      logical, parameter :: local_debug = .true.,msg_flag = .false.
      double precision, parameter :: d32460 = 32460.0d0
c
c          at least one convergence test for global Newton iterations
c          must be defined to proceed
c
      call errchk( 10, 0, .false. )
c
c          clean up any temporary data structures left
c          over from parsing the input.
c
      if ( allocated( temp_nodmap ) ) deallocate( temp_nodmap )
      if ( allocated( temp_nodlod ) ) deallocate( temp_nodlod )
c
c          before analysis of load (time) step 1, perform
c          system-wide array initializations.
c
      if( (ltmstp .eq. 0) .and. (.not. incflg) ) call incomp
c
c          if we are running MPI:
c          check if either new constraints or new dynamic
c          analysis parameters were input during last step.
c          If so, then send the required data to the slave
c          processors. If this is the very first time step,
c          then skip these routines because they were called
c          inside of incomp.  Also, if we are using crack
c          growth by element extinction, then we need to check
c          to see if the data structures have been allocated.
c
      if( ltmstp .ge. 0) then
         if ( new_constraints ) call wmpi_send_const
         if ( new_analysis_param ) call wmpi_send_analysis
      end if
      if( growth_by_kill ) call wmpi_growth_init
c
c          before this compute displ... command, the user may have
c          changed the constraints.  if analysis uses crack
c          growth by node release, we must also check for
c          any previously released nodes are now reconstrained.
c
      if( growth_by_release .and. new_constraints )
     &    call chk_reconstraint
c
c          find the range of time steps for which this
c          loading is defined.
c
      lowstp = stprng(ldnum,1)
      histep = stprng(ldnum,2)
c
c          -----   loop over steps in user defined list -----
c          extract each step from user step list and process.
c          intermediate steps maybe required, e.g., compute
c          for step 20 when 10-19 have not been computed yet.
c
      icn    = 0; iplist = 1
      do while ( iplist .ne. 0 )
        call trxlst( stplst, nsteps, iplist, icn, step )
c
c          is time step > max steps defined ?
c          is time step requested < 0 ?
c          is time step < previous step solved ?
c
        if( step .gt. max_step_limit ) then
          call errmsg( 125, step, dums, dumr, dumd )
          cycle
        end if
        if( step .lt. 0 ) then
          call errmsg( 69, step, dums, dumr, dumd )
          cycle
        end if
        if( step .le. ltmstp ) then
          call errmsg( 123, step, dums, dumr, dumd )
          cycle
        end if
c
        diff = step - ltmstp
        if( diff .gt. 1 ) then ! fill in intermediate steps
         do while ( diff .ge. 1 )
           istep = ltmstp + 1
           diff  = diff - 1
           call stpdrv_one_step( istep, stpdrv_error )
           if( stpdrv_error ) return
           call stpdrv_output( istep ) ! output commands file
         end do
         cycle
        else  ! no fill in needed
          call stpdrv_one_step( step, stpdrv_error )
          if( stpdrv_error ) return
          call stpdrv_output( step )  ! output commands file
        end if   !  diff .gt. 1
c
      end do
      return
c

c
      contains     ! ***** note contains here *****


c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine stpdrv_output                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/10/2018 rhd             *
c     *                                                              *
c     *     drive the optional, user-define output commands stored   *
c     *     in a file completion of a load (time) step               *
c     *                                                              *
c     ****************************************************************
c
      subroutine stpdrv_output( now_step )
      implicit none
c
      integer :: now_step
c
c
c          locals
c
      logical :: matchs, sflag_1, sflag_2, here_debug, is_iostat_end,
     &           output_this_step
      logical, external :: ouchk_map_entry
      character(len=80) :: line
      integer :: read_status, word, bit
c
c          check to see if the user has chosen to use the
c
c             output commands use file .. after steps ....
c
c          command.  if so, drive the output system to execute
c          commands in the user defined file.
c
c          we call oudrive directly from here to process each (logical)
c          line from the file once we pre-read the line to determine
c          it is an output command.
c
c          only output commands and comment lines are permitted in
c          the specified file of commands.
c
c          This logic is a bit tricky since we are mingling
c          peeks at the next input lines before scan, then
c          backspacing and invoking scan it all looks ok.
c
c          if the user has defined the file and we just completed
c          a load (time) step indicated in the command, we open anc
c          run commands in the output file. Then close the
c          output commands file.
c

      here_debug = .false.
      if( here_debug ) write(out,*) " entered stpdrv_output"
      if( .not. ouchk_map_entry(now_step) ) return
c
c          open the file of output commands and make scan aware
c          of it. the global variable "in" has the i/o device number
c
      call infile_stpdrv_open( output_command_file )
      if( here_debug ) write (out,*) "opened output commands"
c
      do  ! over lines in output commands file
c
c            1. read a line from output commands file.
c            2. if eof, close and return to stpdrv
c            3. skip comment lines
c            4. backspace and have scan read the line
c            5. if output command, call oudrive to complete
c               scanning of line and perform output
c            6. if not an "output" command, kill job. we're
c               deep in the solution logic, not overall command
c               processing logic in main program
c
c            this effort is required since readsc will not let us
c            check for an eof. it sees an eof, pops the file stack
c            & starts reading again from prior file. We don't want
c            that to happen here in the code driving the solution
c            over load steps. the comments is check is needed here
c            as well. support he last lines of the output commands
c            file has comment lines. readsc seems them, skips over them,
c            hits and eof and pops the stack - not what we want here.
c            the eof check here and the comment checks make it work.
c
        read(in,fmt="(a80)",iostat=read_status) line
        if( is_iostat_end( read_status ) ) then
          if( here_debug ) write(out,*) " @ 5 eof condition"
          call infile_stpdrv_close( output_command_file )
          return
        end if
        if( line(1:2) .eq. "! " ) cycle
        if( line(1:2) .eq. "c " ) cycle
        if( line(1:2) .eq. "C " ) cycle
        if( line(1:) .eq. " " ) cycle
        backspace (unit=in)
        call readsc
        if( here_debug ) write(out,*) " @ 1 just readsc"
        if( matchs("output",4) ) then
           call oudrive( sflag_1, sflag_2, stname, ltmstp )
        else  ! only output commands & comments allowed
           write(out,9100)
           call die_gracefully
        end if
      end do
c
 9100 format(
     & /1x, '>>>>> error: while processing file of output commands.',
     & /14x,'only output commands and comment lines allowed.',
     & /14x,'job terminated...',//)
 9200 format(
     & /1x, '>>>>> error: no list of steps exists for output commands',
     & /14x,'file: ',a80,
     & /14x,'job terminated...',//)
c
      end subroutine stpdrv_output
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine stpdrv_one_step                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 8/21/2016 rhd              *
c     *                                                              *
c     *            oversee setting up solution for one step          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine stpdrv_one_step( now_step, stpdrv_error  )
      implicit none
c
      integer :: now_step
      logical :: stpdrv_error
c
      integer :: i
c
      stpdrv_error = .false.
c
      if( output_packets ) then
        call close_packets_file( msg_flag )
        call open_packets_file( msg_flag )
      end if
c
c          check to make sure the load-time step to be com-
c          puted is defined in the loading given. if
c          not, then cease processing of the loading
c          list given.
c
      if( (now_step .lt. lowstp) .or. (now_step .gt. histep) ) then
        write(out,9121) now_step, lodnam(ldnum)
        num_error = num_error + 1
        stpdrv_error = .true.
        call die_abort
      else if ( .not. stpchk(now_step) ) then
        write(out,9121) now_step, lodnam(ldnum)
        num_error = num_error + 1
        stpdrv_error = .true.
        call die_abort
      end if
c
c          the step is a valid one. will load step
c          require more than wall time allowed.
c          Yes - make restart & quit.
c
c          if option is in effect, call the user routine to
c          modify nonlinear solution parameters and loadings
c          for next step.
c
      call steptime( now_step, 2 )
c
c          if option is in effect, call the user routine to
c          modify nonlinear solution parameters and loadings
c          for next step.
c
      if( run_user_solution_routine )
     &   call stpdrv_user_solution_parms( now_step )
c
c          compute: (1) increment in nodal and element
c          temps specified by use for step, (2) the
c          total nodal forces specified by user
c          and equivalent element nodal forces
c          for end of step conditions. put into
c          user specified nodal coordinate systems.
c
      if( now_step .eq. 1 ) then
        mf = 1.0
        mf_nm1 = 1.0
        mf_ratio_change = .false.
      end if
c
      call eqivld( mf, mf_nm1, mf_ratio_change, now_step, ldnum )
c
c          check if any algorithms require a smaller load step
c          size.  If so, reduce the constraints, incremental
c          loads, incremental temperatures and accumulated
c          pattern multipliers accordingly.
c          set the constraint vector for
c          the step from the user defined constraint values.
c          remember that users can specifiy constraints
c          in the step loading with a multiplier stored
c          in cnstrn_stp_factors - default value is 1.0.
c          make sure nodes w/o constraints are not changed.
c
c          calls for check.. and modify.. needed even for step 1
c          to set up other data structures/global variables
c
      con_stp_factor = user_cnstrn_stp_factors(now_step)
      do i = 1, nodof
         cnstrn(i) = cnstrn_in(i) * con_stp_factor
         if ( cnstrn_in(i) .eq. d32460 ) cnstrn(i) = d32460
      end do
      if( now_step == 1 ) then
             load_reduce_fact = 1.0d0
      else
        call check_for_step_reduction( load_reduce_fact, mf,
     &                                mf_nm1 )
      end if
      call modify_load( load_reduce_fact, mf, mf_nm1, now_step )
c
c          flags for new constraints and new dynamic
c          analysis parameters
c
      new_constraints    = .false.
      new_analysis_param = .false.
c
c          execute the iteration driver to compute
c          displacements for this load-time step.
c
      call mnralg( mf, mf_nm1, mf_ratio_change, now_step, ldnum )
c
c          if global load step size was reduced via  reduction
c          checks, then change it back to the original size.
c
      call original_step_size( mf, mf_nm1, now_step )
c
c          possibly update/output wall time info
c
      call steptime( now_step, 3 )

c
      return
c
 9121 format(/1x,'>>>>> FATAL ERROR: the load step to be solved: ',i7,
     &       /1x,'                   is not defined for loading: ',a8,
     &       /1x,'                   job terminated....')

 9160 format(7x,
     & '>> computing first element stiffness matrices (@ t=0)')
c
      end subroutine stpdrv_one_step
c
      end subroutine stpdrv

c
c
c     ****************************************************************
c     *                                                              *
c     *              subroutine stpdrv_user_solution_parms           *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/6/2015 rhd              *
c     *                                                              *
c     *     drive execution of user_solution_parms routine           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine stpdrv_user_solution_parms( now_step )
      use global_data ! old common.main
c
      use damage_data, only :  perm_load_fact
      use main_data, only : step_load_data, convergence_history,
     &                      extrapolate, extrap_off_next_step,
     &                      divergence_check, diverge_check_strict,
     &                      line_search, user_cnstrn_stp_factors
      implicit none
c
      integer :: now_step
c
      integer :: num_patt, i, ld_cond_number
      include 'include_usr_parm'
      type(solution_parameters) :: usolution_parms
      type(step_definition) :: next_step_loading
c

c
c          put current values of solution parameters into local data
c          structure for passing to user solution routines. Reset
c          system values to those returned by the user routine
c
      usolution_parms%simulation_time = total_model_time
      usolution_parms%time_step = dt
      usolution_parms%newmark_beta = nbeta
c
      usolution_parms%maximum_iterations = mxiter
      usolution_parms%minimum_iterations = mniter
      usolution_parms%nonconvergent_solutions_flag = halt
      usolution_parms%adaptive_solution = adaptive_flag
      usolution_parms%extrapolate_solution = extrapolate
      usolution_parms%extrap_off_next_step = extrap_off_next_step
      usolution_parms%divergence_check = divergence_check
      usolution_parms%diverge_check_strict = diverge_check_strict
      usolution_parms%line_search = line_search

c
      usolution_parms%convergence_test_1 = convrg(1)
      usolution_parms%convergence_test_1_tolerance = tol(1)
      usolution_parms%convergence_test_2 = convrg(2)
      usolution_parms%convergence_test_2_tolerance = tol(2)
      usolution_parms%convergence_test_3 = convrg(3)
      usolution_parms%convergence_test_3_tolerance = tol(3)
      usolution_parms%convergence_test_4 = convrg(4)
      usolution_parms%convergence_test_4_tolerance = tol(4)
c
      usolution_parms%batch_messages = batch_messages
      usolution_parms%wall_time_limit_seconds = time_limit
      usolution_parms%material_messages =  signal_flag
      usolution_parms%trace_solution = trace(1)
c
      usolution_parms%bbar_stabiliation_factor = eps_bbar
      usolution_parms%consistent_q_matrix = qbar_flag
c
      usolution_parms%reset_load_reduction_factor = perm_load_fact
c
c          put definition of current load (time) step into local data
c          structure for passing to user solution routines. Reset
c          system values to those returned by the user routine
c
      num_patt = step_load_data(now_step)%num_load_patterns
      next_step_loading%number_load_patt = num_patt
      if( num_patt .gt. 0 ) then
        do i = 1, num_patt
          ld_cond_number = step_load_data(now_step)%load_patt_num(i)
          next_step_loading%load_patt_nums(i) = ld_cond_number
          next_step_loading%load_patt_ids(i) = lodnam(ld_cond_number)
          next_step_loading%load_patt_multipliers(i) =
     &       step_load_data(now_step)%load_patt_factor(i)
        end do
        deallocate( step_load_data(now_step)%load_patt_num )
        deallocate( step_load_data(now_step)%load_patt_factor )
      end if
c
      next_step_loading%constraints_multiplier =
     &   user_cnstrn_stp_factors(now_step)
c
c          call the user supplied routine to potentially modify
c          parameters, loading definition for next step
c
      call user_solution_parameters(
     &     now_step, out, usolution_parms, next_step_loading,
     &     convergence_history )
c
c          set WARP3D values to those possibly modified by the
c          user routine
c
      dt    =  usolution_parms%time_step
      nbeta =  usolution_parms%newmark_beta
c
      mxiter = usolution_parms%maximum_iterations
      mniter = usolution_parms%minimum_iterations
      halt   = usolution_parms%nonconvergent_solutions_flag
c
      adaptive_flag = usolution_parms%adaptive_solution
      extrapolate   = usolution_parms%extrapolate_solution
      extrap_off_next_step = usolution_parms%extrap_off_next_step
      divergence_check = usolution_parms%divergence_check
      diverge_check_strict = usolution_parms%diverge_check_strict
      line_search = usolution_parms%line_search
c
      convrg(1) = usolution_parms%convergence_test_1
      tol(1)    = usolution_parms%convergence_test_1_tolerance
      convrg(2) = usolution_parms%convergence_test_2
      tol(2)    = usolution_parms%convergence_test_2_tolerance
      convrg(3) = usolution_parms%convergence_test_3
      tol(3)    = usolution_parms%convergence_test_3_tolerance
      convrg(4) = usolution_parms%convergence_test_4
      tol(4)    = usolution_parms%convergence_test_4_tolerance
c
      batch_messages = usolution_parms%batch_messages
      time_limit = usolution_parms%wall_time_limit_seconds
      signal_flag = usolution_parms%material_messages
      trace(1) = usolution_parms%trace_solution
c
      eps_bbar = usolution_parms%bbar_stabiliation_factor
      qbar_flag = usolution_parms%consistent_q_matrix
c
      perm_load_fact = usolution_parms%reset_load_reduction_factor
c
c          update loading definition for next load (time) step. delete
c          current vector, re-allocate to possibly new size and load
c          updated values from user routine.
c
      user_cnstrn_stp_factors(now_step) =
     &    next_step_loading%constraints_multiplier
c
      num_patt = next_step_loading%number_load_patt
      step_load_data(now_step)%num_load_patterns = num_patt
c
      if( num_patt .gt. 0 ) then
        allocate( step_load_data(now_step)%load_patt_num(num_patt) )
        allocate( step_load_data(now_step)%load_patt_factor(num_patt) )
        step_load_data(now_step)%load_patt_num(1:num_patt) =
     &     next_step_loading%load_patt_nums(1:num_patt)
        step_load_data(now_step)%load_patt_factor(1:num_patt) =
     &       next_step_loading%load_patt_multipliers(1:num_patt)
       end if
c
      return
      end

