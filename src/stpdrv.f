c     ****************************************************************
c     *                                                              *
c     *                      subroutine stpdrv                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 7/24/22 rhd                *
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
      use main_data, only : cnstrn, cnstrn_in, temp_nodmap,
     &  temp_nodlod, output_packets, run_user_solution_routine,
     &  output_command_file, max_step_limit, stpchk, 
     &  extrap_off_next_step, user_cnstrn_stp_factors,
     &  user_cnstrn_stp_factors, last_step_adapted,
     &  last_step_num_iters
c
      use damage_data, only : growth_by_kill, growth_by_release
      use constants
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
      logical :: matchs, sflag_1, sflag_2, here_debug, is_iostat_end
      logical, external :: ouchk_map_entry
      character(len=80) :: line
      integer :: read_status
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
c     *                   last modified : 11/14/22  rhd              *
c     *                                                              *
c     *            oversee setting up solution for one step          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine stpdrv_one_step( now_step, stpdrv_error  )
c
      use j_data, only : J_cutoff_active, J_compute_step_2_automatic,
     &                   J_max_step_1, J_cutoff_e, J_cutoff_nu,
     &                   J_max_now_step, J_auto_step_2_delta_K
c
      implicit none
c
      integer :: now_step
      logical :: stpdrv_error
c
      integer :: i
c
c          optional processing for J adaptive/limiting definitions of 
c          this and future load steps.
c
c          features:
c
c             - J_cutoff_active: monitor the value of J/J_elastic at 
c                   the local of maximum J in each step. if value exceeds
c                   user limit, end job before starting this step.
c                   optionally write a restart file.
c
c             - J_ratio_adaptive_steps: implies Kr adaptive step loads.
c                   Early in loading compute the change in Kr over
c                   the steps. Increases/decrease step sizes
c                   to best maintain the user-specified change.
c                   at more plastic deformation, switch to use the change in
c                   in J/J_elastic ratio over the just completed step
c                   for adaptive.
c                   if the change exceeds a user specified value,
c                   e.g. 0.5, reduce this and subsequent step sizes.
c                   if the change is too small, increase this
c                   and subsequent load step sizes.
c                   the result is a dynamically varying size of each
c                   step aiming to maintain the user specified change
c                   (e.g. 0.5)
c                   see detail description of decisions to set the
c                   adaptive step increment size in Section 2.18 of manual
c 
c             - J_compute_step_2_automatic: the user specifies a loading
c                   for step 1 sufficiently small to insure a linear-
c                   elastic response. With large parametric studies,
c                   the definition of sizes for load steps 2, 3, 4, ..
c                   can become challenging for manual processing.
c                   strategy here: compute max K_I found along front
c                   in step 1. Set pattern factors for steps 2, 3, ...
c                   to increase K_I by a specified amount assuming the
c                   solution remains linear. E.g \Delta K_I in step 2
c                   set to be 25 based on step 1 solution.
c                   This feature combined with J_ratio_adaptive_steps
c                   that kicks in after step 3
c                   eliminates the need for careful definition of
c                   load steps sizes
c
c             The Kr-J adaptive step loading works only for load
c             control simulations. All displacement constraints
c             must = 0.
c
c             Solutions should *not* use the nlgeom formulation. After some
c             amount of displacements, the equivalent force vector entries
c             will no longer be proportional to those in load step 1. This
c             invalidates scaling of KI in step 1 to step n.
c  
c          now_step is the step number we are about to compute
c          displacements
c
      if( J_compute_step_2_automatic .and. now_step == 2 )
     &    call stpdrv_J_auto_size_step_2 ! just finished step 1
c
      if( J_cutoff_active ) call stpdrv_J_cutoff( now_step ) ! may just return
c
c          end of J cutoff, adaptive J-Jr based processing
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
        mf = one
        mf_nm1 = one
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
      current_load_time_step = now_step
      last_step_adapted = .false.
      last_step_num_iters = 0
      call mnralg( mf, mf_nm1, mf_ratio_change, now_step, ldnum )
c
c          if global load step size was reduced via reduction
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
9121  format(/1x,'>>>>> FATAL ERROR: the load step to be solved: ',i7,
     &       /1x,'                   is not defined for loading: ',a8,
     &       /1x,'                   job terminated....')
c
      end subroutine stpdrv_one_step
c
c     ****************************************************************
c     *                                                              *
c     *            subroutine stpdrv_J_auto_size_step_2              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 7/14/2022 rhd              *
c     *                                                              *
c     *            set size of step 2 based on K_I max along front   *
c     *            found in step 1.                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine stpdrv_J_auto_size_step_2
c
      use j_data, only : J_max_step_1, J_cutoff_e, J_cutoff_nu,
     &                   J_auto_step_2_delta_K
      use main_data, only : step_load_data
c
      implicit none
c
      integer :: i, j, num_defined_steps, npatt
      double precision :: step_1_patt_factor, K_max_step_1,
     &                    step_2_factor 
      logical, parameter :: here_debug = .false.
c
      K_max_step_1 = sqrt( J_cutoff_e * J_max_step_1 /
     &                     (one - J_cutoff_nu**2) )
      step_2_factor = J_auto_step_2_delta_K / K_max_step_1
c
c              modify pattern factors for steps 2->last defined
c              to make delta K_I = wanted value assuming linear
c              solution in step 2. most likely to
c              be combined with J_ratio_adaptive_steps feature
c              to subsequently adjust sizes for step3 3, 4, 5 ,..
c              for increased nonlinear behavior.
c
      num_defined_steps = size( step_load_data )
      do i = 2, num_defined_steps
        npatt = step_load_data(i)%num_load_patterns
        do j = 1, npatt
          step_load_data(i)%load_patt_factor(j) = 
     &        step_load_data(1)%load_patt_factor(j) * step_2_factor
        end do
        user_cnstrn_stp_factors(i) = user_cnstrn_stp_factors(1) *
     &                                 step_2_factor 
      end do
c
      write(out,9100)
      write(out,8900) J_max_step_1
      write(out,9000) K_max_step_1
      write(out,9010) J_auto_step_2_delta_K
      write(out,9020) step_2_factor
c
      return
c   
 8900 format(   '      maximum J on front for step 1:   ',e14.6)
 9000 format(   '      maximum K_J on front for step 1: ',e14.6)
 9010 format(   '      requested delta K in step 2:     ',e14.6)
 9020 format(   '      multiplier on step 1 loading =   ',e14.6)  
 9100 format(//,'>>>>> J-adaptive loading setup for step 2...')
c
      end subroutine stpdrv_J_auto_size_step_2
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine stpdrv_J_cutoff                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/2/23 rhd                 *
c     *                                                              *
c     *       has J_/J-elastic ratio reached user limit ?            *
c     *       and/or adapt this and future steps sizes to maintain   *
c     *       a user-specified target for the ratio change in        *
c     *       each step                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine stpdrv_J_cutoff( now_step )
c
      use j_data, only :
     &     J_cutoff_active, J_cutoff_exceeded,
     &     J_cutoff_restart_file, J_count_exceeded,
     &     J_cutoff_num_frnt_positions, J_max_now_step,
     &     J_cutoff_max_value, J_cutoff_ratio, 
     &     J_cutoff_frnt_pos_max_ratio, J_ratio_last_step,
     &     J_target_diff, J_ratio_adaptive_steps,
     &     J_limit_ratio_increase, J_limit_ratio_decrease,
     &     J_diff_at_2_set, J_diff_at_2, Kr_min_limit,
     &     Kr_target_diff
      implicit none
c
      integer :: now_step
c
      integer :: step_just_completed  
      logical :: ldummy1, ldummy2, ratio_set,
     &           convergence_caused_reduction
      logical, parameter :: here_debug = .false. 
      double precision :: J_ratio_diff, J_load_factor,
     &                    diff_ratio, slope, now_target_diff,
     &                    Kr_now, Kr_last_step, now_target_diff1
c
c                now_step is the step number about to be computed.
c
c                we first examine the step change in J/J_elastic
c                values once step 3 is completed, i.e., the first
c                difference is ratio @ 3 - ratio @ 2.
c
c                J_cutoff_max_value: J/Je value at location on front
c                  where J max occurs. Th
c                  front location at which this value occurs
c                  often moves during early loading as plastic
c                  deformation evolves
c                J_cutoff_exceeded (logical) set by di_process_J_cutoff
c                  if J_cutoff_max_value exceed user-specified limit
c                J_count_exceeded set by di_process_J_cutoff. count of
c                number frnt positions at which J/Je exceeds user limit
c
      step_just_completed = now_step - 1
      if( step_just_completed >= 3 ) then ! write summary
        J_ratio_diff = J_cutoff_max_value - J_ratio_last_step
        Kr_now = one / sqrt( J_cutoff_max_value ) ! save
        Kr_last_step = one / sqrt( J_ratio_last_step ) ! save
        write(out,9200) step_just_completed, J_cutoff_ratio, 
     &         J_count_exceeded, J_cutoff_num_frnt_positions,
     &         J_max_now_step, J_cutoff_frnt_pos_max_ratio,
     &         J_cutoff_max_value,
     &         J_ratio_diff, one/sqrt(J_cutoff_max_value) ! FAD vlaue
      end if
      if( step_just_completed >= 2 ) then ! update ratio for last step
        J_ratio_last_step = J_cutoff_max_value
        if( .not. J_ratio_adaptive_steps ) write(out,9235)
     &          step_just_completed, J_ratio_last_step
      end if
c
      if( J_cutoff_exceeded ) then ! stop analysis
         write(out,9205) 
         if( J_cutoff_restart_file ) then
            write(out,9210)
            call store( ' ','J_ratio_limit_exceeded.db', 
     &                 ldummy1, ldummy2 )
         end if
         write(out,9220)
         call warp3d_normal_stop
      end if 
c
c              continue to next load step with analysis.
c
      if( .not. J_ratio_adaptive_steps ) return 
c
c              Start using change in Kr over the step to set the
c              step increment multiplier for next step.
c
c              Once Kr decreases to about 0.6-0.8 (user-specified),
c              we switch to use changes in J/J_e to set step 
c              increment multiplier for next step.
c
c              see figures in Section 2.18 for details of
c              decision process.
c                
c              These conditions seem to prevent reducing or 
c              increasing load steps sizes to rapidly   
c
c              The diff_ratio = J_ratio_diff - J_target_diff can be
c              slightly negative in early steps when the position of
c              max J/J_e changes along the front.
c
c              Kr = 1.0 / sqrt( J/Je ) continually decreases from 1
c                                      with loading. some increases
c                                      above 1.0 to say 1.1-1.2 are
c                                      observed for cerftain crack
c                                      geometries/loading
c              J/Je = 1.0 / Kr**2
c         
      if( step_just_completed < 3 ) return ! no adaptive loading yet
c
c               use change in Kr over each step for adaptive loading 
c               until Kr decreases to the specified level, e.g. 0.6
c               this scheme is very simple
c
      if( Kr_now > Kr_min_limit ) then ! still under Kr size control
         call stpdrv_J_cutoff_Kr( now_step, Kr_now, Kr_last_step,
     &                            Kr_target_diff )
         return
      end if
c       
c               under J/Je control at higher deformations
c
c               save J/Je ratio when we switched to J/Je control.
c
      if( .not. J_diff_at_2_set ) then ! time to set J_diff_at_2
          J_diff_at_2 = J_ratio_diff
          J_diff_at_2_set = .true.
          write(out,9240)  ! switching to J/Je control
      end if
c
c               increases linearly the target delta(J/Je) for next
c               step is we are in a simple transition region of 
c               delta(JKr) control to delta(J.Je) control.
c               tries to limit a big jump in the target change for
c               delta(J/Je).
c
       now_target_diff = J_target_diff ! user-specified input value 
       if( J_cutoff_max_value < eight ) then ! interpolate to get target
        slope = ( J_target_diff - J_diff_at_2 ) / (eight - two )
        now_target_diff1 = J_diff_at_2 + (J_cutoff_max_value-two)*slope
        now_target_diff = min( J_target_diff, now_target_diff1 )
        if( here_debug ) then
           write(out,*) "..  J_target_diff ,J_diff_at_2: ",
     &       J_target_diff ,J_diff_at_2
           write(out,*) ".. slope, now_tar_diff1: ",slope,
     &           now_target_diff1 
           write(out,*) ".. now_tar_diff: ", now_target_diff 
        end if
      end if
c
c               basic scheme now applied that now_target_diff
c               for delta(J/Je) set.
c
      ratio_set = .false.  ! simplifies logic below
      diff_ratio = J_ratio_diff - now_target_diff
      if( diff_ratio < zero ) then ! max position on frnt likely moved
         J_load_factor = J_limit_ratio_increase
         ratio_set = .true.
      else
         J_load_factor = now_target_diff / J_ratio_diff
      end if
c
c               limit algorithm predicted adjustments in load factor
c               to dampen/accelerate loading rate (based on experience)
c
      if( .not. ratio_set ) then
        if( J_load_factor < one ) then  ! decrease load step size
          ratio_set = .true.
          if( J_load_factor < ptone ) J_load_factor = ptone ! limit decrease
          if( J_load_factor > pt75 ) J_load_factor = pt75
        end if
      end if
c
      if( .not. ratio_set ) then
        if( J_load_factor > J_limit_ratio_increase ) ! increase step sizes
     &       J_load_factor = J_limit_ratio_increase ! default = 1.1
      end if  
c
c               if last load step caused issues with global Newton
c               convergence, reduce next step size by 75%.
c
      convergence_caused_reduction = .false.
      if( J_load_factor > quarter ) then
        if( last_step_adapted .or. last_step_num_iters > 8 ) then 
          J_load_factor = quarter
          convergence_caused_reduction = .true.
        end if
      end if
c
      write(out,9230) now_target_diff, J_load_factor
      if( convergence_caused_reduction ) write(out,9232)
c
      if( here_debug ) then
        write(out,*) ' @1   step_just_completed: ',step_just_completed
        write(out,*) '      J_cutoff_max_value: ', J_cutoff_max_value
        write(out,*) '      J_ratio_diff: ', J_ratio_diff
      end if  
      call stpdrv_J_adapt_scale_loads( now_step, J_load_factor  )
c
      return
c
9200  format(//,'>>>>> Summary for J-cutoff after step: ',i6,
     & /,'               user limit: ',f5.1,
     &   ' exceeded at: ',i3, ' of: ',i3, ' crack front positions',
     & /,'               max J on front ',e14.6,
     & ' at position: ',i5 ,' J/Je ratio: ', f12.4,
     & /,'               change in J-ratio over previous step: ',f7.3,
     & /,'               Kr = K_I / K_J (FAD): ',f6.3)
9205  format(//,'>>>>> User-specified limit on J/J_elastic ',
     &    ' exceeded ...' )
9210  format(/, '      Writing restart file: ',
     &        ' J_ratio_limit_exceeded.db' )
9220  format(//, '>>>> Job terminated normally...',//)
9230  format( '               current target J-ratio change'
     &   ' for next step: ',f7.3,
     &        /,'               multiplier applied to next & future ',
     & '  steps: ',f6.2)
9232  format(   '               convergence behavior of global Newton ',
     & 'iterations for prior step governs multiplier')
9235  format(//,'>>>>> INFO: Max J-ratio for step: ',i6,2x,f7.2)
9240  format('               switching to J/Je load control from Kr')
c
      end subroutine stpdrv_J_cutoff

c     ****************************************************************
c     *                                                              *
c     *                 subroutine stpdrv_J_cutoff_Kr                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 4/11/24 rhd                *
c     *                                                              *
c     *       Adaptive load factors based on increments of Kr        *
c     *       over load steps to capture early part of FADs          *
c     *                                                              *
c     ****************************************************************

      subroutine stpdrv_J_cutoff_Kr( now_step, Kr_now, Kr_last_step,
     &                               Kr_target_diff )
c
      implicit none
c
      integer :: now_step
      double precision :: Kr_now, Kr_last_step, Kr_target_diff
c
      logical, parameter :: here_debug = .true. 
      double precision ::  Kr_diff, Kr_load_factor
c
      Kr_diff = - ( Kr_now - Kr_last_step ) ! so we can use + values
      write(out,9100) Kr_last_step, -Kr_diff,
     &                -Kr_target_diff
c
c              if last Kr > 1.0 + tol, decrease the step size. This is
c              likely a localized unloading situation on the crack front.
c              when Kr > 1, J_el > J_total, where J_el is the scaled value
c              from linear-elastic solution in load step 1. J_el > J_total
c              (from the domain integral) is not realistic. It means
c              J_total (domain) is decreasing due to near tip (elastic)
c              unloading.
c              Current scheme reduces next load steps size a small amount.
c
      if( Kr_now >= 1.005d0 ) then
        Kr_load_factor = 0.9d0 ! small decrease in load step size
        write(out,9120) Kr_load_factor
        call stpdrv_J_adapt_scale_loads( now_step, Kr_load_factor  )
        return
      end if         
c
c              may want to increase load step sizes a bit. max increase
c              factor 1.1 for now
c
      if( Kr_diff <= Kr_target_diff ) then
        if( Kr_diff > pt75 * Kr_target_diff ) then
           write(out,9120) one
           return ! no change
        end if
        Kr_load_factor = oneptone ! increase load steps
        write(out,9120) Kr_load_factor
        call stpdrv_J_adapt_scale_loads( now_step, Kr_load_factor  )
        return
      end if
c
c              decrease load step sizes. minimum decrease factor 
c              is 0.8 even if 0.9 indicated
c
      Kr_load_factor = Kr_target_diff / Kr_diff   ! reduce load steps
      Kr_load_factor = min( Kr_load_factor, point_eight ) 
      write(out,9120) Kr_load_factor
      call stpdrv_J_adapt_scale_loads( now_step, Kr_load_factor  )
c
      return    
c
 9100 format('               Kr last step, Kr_diff, target Kr diff: ',
     & 3f8.4)
 9130 format('               no change in load factor next step')
 9120 format('               load factor next step: ',f5.2)   
c
      end subroutine stpdrv_J_cutoff_Kr
c
c     ****************************************************************
c     *                                                              *
c     *            subroutine stpdrv_J_adapt_scale_loads             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 8/8/2022 rhd               *
c     *                                                              *
c     *       scale existing pattern factors for                     *
c     *       about-to-be-computed and all subsequent steps by       * 
c     *       the passed factor                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine stpdrv_J_adapt_scale_loads( now_step, factor )
c
      use main_data, only : step_load_data
      implicit none
c
      integer :: now_step, i, j, num_defined_steps, npatt
      double precision :: factor
      logical, parameter :: here_debug = .false.
c
      num_defined_steps = size( step_load_data )
      if( here_debug ) then
       write(*,*) '.... entered stpdrv_J_adapt_scale_loads'
       write(*,*) '        now_step, factor: ', now_step, factor
       write(*,*) '        num_defined_steps: ',num_defined_steps 
       write(*,*) '        before updates:'
       do i = 1, 6
        npatt = step_load_data(i)%num_load_patterns
        write(*,*) "      > i, npatt: ",i, npatt 
        do j = 1, npatt
          write(*,*) '       j, stepdata: ', j,
     &     step_load_data(i)%load_patt_factor(j) 
        end do
      end do
      end if
c
      do i = now_step, num_defined_steps
        npatt = step_load_data(i)%num_load_patterns
        do j = 1, npatt
          step_load_data(i)%load_patt_factor(j) = 
     &        step_load_data(now_step-1)%load_patt_factor(j) * factor
        end do
        user_cnstrn_stp_factors(i) = 
     &      user_cnstrn_stp_factors(now_step-1) * factor 
      end do
c
      if( here_debug ) then
        write(*,*) '        after updates:'
        do i = 1, 6
         npatt = step_load_data(i)%num_load_patterns
         write(*,*) "      > i, npatt: ",i, npatt 
          do j = 1, npatt
            write(*,*) '       j, stepdata: ', j,
     &         step_load_data(i)%load_patt_factor(j) 
          end do
        end do
      end if
c
      if( now_step == 4 .and. here_debug ) then
         write(*,*) '  @ 3 now_step = 4'
          do i = now_step, num_defined_steps
            npatt = step_load_data(i)%num_load_patterns
              write(*,*) "          npatt: ", npatt
              do j = 1, npatt
                write(*,*) "      ",
     &                     step_load_data(i)%load_patt_factor(j)
              end do
           end do
      end if
         
      return
      end subroutine stpdrv_J_adapt_scale_loads
c
      end subroutine stpdrv
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
     &    real( next_step_loading%constraints_multiplier )
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

