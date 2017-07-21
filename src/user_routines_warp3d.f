c
c           user_routines_warp3d.f   Distribution version
c
c           Updated:  2/6/2015
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine user_solution_parmeters           *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/7/2015 rhd              *
c     *                                                              *
c     *     user routine to potentially alter solution parameters    *
c     *     and loading before the next load (time) step             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine user_solution_parameters( now_step, iout,
     &    sol_parms, next_step_loading, convergence_history )
      implicit none
c
      integer now_step, iout ! others are user-defined types
c
c           declaration of data structure used to communicate between
c           WARP3D and the user_solution_params routine
c
c           See Section 2.10 of User Manual for more details. The
c           variables below have 1-to-1 correspondence with
c           user input commands.
c
      include 'include_usr_parm'
      type :: step_convergence_data
        logical :: step_converged
        logical :: adaptive_used
        integer :: iterations_for_convergence
        integer :: adapt_substeps
      end type
c
      type(solution_parameters) :: sol_parms
      type(step_definition) :: next_step_loading
      type(step_convergence_data), dimension(5) ::
     &                         convergence_history
c
c            locals
c
      integer i, num_patts, step_to_write
      logical debug_local
      data debug_local  / .false. /
c
      if( debug_local ) then ! values provided by WARP3D
        write(iout,*)
     &   '... called user_solution_parameters before step: ',now_step
        write(iout,9000) sol_parms%time_step, sol_parms%newmark_beta,
     &   sol_parms%maximum_iterations, sol_parms%minimum_iterations,
     &   sol_parms%nonconvergent_solutions_flag,
     &   sol_parms%adaptive_solution,
     &   sol_parms%extrapolate_solution
        write(iout,9010)
     &   sol_parms%convergence_test_1,
     &   sol_parms%convergence_test_1_tolerance,
     &   sol_parms%convergence_test_2,
     &   sol_parms%convergence_test_2_tolerance,
     &   sol_parms%convergence_test_3,
     &   sol_parms%convergence_test_3_tolerance,
     &   sol_parms%convergence_test_4,
     &   sol_parms%convergence_test_4_tolerance
        write(iout,9020)
     &   sol_parms%batch_messages,
     &   sol_parms%wall_time_limit_seconds,
     &   sol_parms%material_messages,
     &   sol_parms%trace_solution,
     &   sol_parms%bbar_stabiliation_factor,
     &   sol_parms%consistent_q_matrix,
     &   sol_parms%reset_load_reduction_factor
c
        num_patts =  next_step_loading%number_load_patt
        write(iout,9150) next_step_loading%constraints_multiplier,
     &             next_step_loading%number_load_patt
        if( num_patts .ne. 0 ) then
          do i = 1, num_patts
            write(iout,9160) next_step_loading%load_patt_nums(i),
     &      next_step_loading%load_patt_ids(i),
     &      next_step_loading%load_patt_multipliers(i)
          end do
        end if
c
c            convergence history over last 5 steps. step n is in 5th
c            row, step n-1 in row 4, etc
c
        if( now_step .gt. 1 ) write(iout,9300) now_step
        step_to_write = now_step - 1
        do i = 5, 1, -1
           if( step_to_write .lt. 1 ) exit
           write(iout,9310) step_to_write,
     &        convergence_history(i)%step_converged,
     &        convergence_history(i)%adaptive_used,
     &        convergence_history(i)%iterations_for_convergence
           if( convergence_history(i)%adaptive_used ) then
             write(iout,9320) convergence_history(i)%adapt_substeps
           end if
          step_to_write = step_to_write - 1
          end do
      end if
c
c            modify nonlinear solution parameters as desired
c            and loading/constraints for next step.
c
      if( debug_local ) then
        write(iout,*) ' '
        write(iout,*) '... leaving user_solution_parameters ...'
      end if

      return
 9000 format(/,6x,'Nonlinear solution parameters: ',
     & /,8x,'time step:                   ',e14.6,
     & /,8x,'newmark beta:                ',f8.3,
     & /,8x,'max, min global iterations:  ',2i3,
     & /,8x,'stop on nonconvergent iters: ',l2,
     & /,8x,'adaptive solution:           ',l2,
     & /,8x,'extrapolate displacements:   ',l2 )
 9010  format(
     &   8x,'convergence test 1, tol:     ',l2,1x,f10.4,
     & /,8x,'convergence test 2, tol:     ',l2,1x,f10.4,
     & /,8x,'convergence test 3, tol:     ',l2,1x,f10.4,
     & /,8x,'convergence test 4, tol:     ',l2,1x,f10.4 )
 9020 format(
     & /,8x,'batch messages:              ',l2,
     & /,8x,'wall time limits (secs):     ',f10.0,
     & /,8x,'material messages:           ',l2,
     & /,8x,'trace solution:              ',l2,
     & /,8x,'bbar stabilize factor:       ',f5.2,
     & /,8x,'consistent [Q] matrix:       ',l2,
     & /,8x,'reset load reduction factor: ',f5.2)
 9150 format(/,6x,'Next step definition:',
     & /,8x,'constraint multiplier:       ',f10.5
     & /,8x,'number of loading patterns:  ',i2 )
 9160 format(15x,'pattern #: ',i5,'   id: ',a8,3x,'multiplier: ',
     & f10.5 )
 9200 format(/ )
 9300 format(/,6x,'Convergence history. Next step is: ',i8)
 9310 format(8x,'Step:', i8,'  Converged: ',l2,
     &  '  Adaptive used: ',l2,
     &  ' global iterations: ', i4 )
 9320 format(23x,'Adaptive substeps used: ',i4)
      end
c
c
c *******************************************************************
c *                                                                 *
c *  User nodal loads routine (optional)                            *
c *                                                                 *
c *******************************************************************
c
c
      subroutine user_nodal_loads(
     &  load_name, user_file_name, pattern_values, initial_temps,
     &  nnode, step_np1, time_n, dtime,
     &  node_coords, forces_set, temps_set, use_linear_k, kout )
      implicit none
c
      character :: load_name*8, user_file_name*80
      integer  nnode, step_np1, kout
      double precision pattern_values(nnode,4), initial_temps(*),
     &                 node_coords(nnode,3),
     &                 time_n, dtime
      logical forces_set, temps_set, use_linear_k
c
c       Routine for direct specification of incremental nodal
c       forces and temperatures for a loading condition (pattern).
c
c       In the input file...
c
c           loading right_end   <--  loading name is right_end
c              nodal loads
c               user_routine (<file name>)   <-- call this routine
c                                                for loads. file name
c                                                is optional
c           loading left
c              element loads
c                  .
c                  .
c
c  Arguments:  floating points are double precision
c  ----------
c    load_name       --  character * 8 (contains string right_end in
c                        above example)
c    user_file_name  --  the character string/label (if any) provided
c                        by user in input after keyword user_routine
c    pattern_values  --  number of model nodes x 4 array. Put
c                        values of force x, y, z and
c                        temperature in cols 1-4.
c                        array was zeroed before call to this routine
c    initial_temps   --  user specified initial temperatures
c    nnode           --  number of nodes in structural model
c    step_np1        --  value of n+1 where solution is being advanced
c                        from step n to n+1 (value = 1 for first
c                        step of analysis)
c    time_n          --  model simulation time at start of step (=0 for
c                        step 1)
c    dtime           --  increment of model simulation time from
c                        n -> n+1
c    node_coords     --  number of model nodes x 3 array. contains
c                        coordinates of all structure nodes (col 1
c                        is x, 2 is y, 3 is z).
c                        These are undeformed coordinates.
c    forces_set      --  (logical) set to .true. if any nodal forces
c                        are defined. initialized to .false.
c                        before call
c    temps_set       --  (logical) set to .true. if any nodal
c                        temperatures are defined.  initialized to .false.
c                        before call
c    use_linear_k    --  return as false or true. return as true
c                        if new step loads are sufficiently non-
c                        proportional that (1) linear stiffness should
c                        be used for iteration 1 of step and (2) no
c                        extrapolation of displacements at start of
c                        step.
c    kout            --  Fortran device number for error messages,
c                        debugging output
c
c               local variables
c               ---------------
c
      double precision x
      integer node
c
c                code below is just a simple example for illustration
c
      forces_set = .true.
      temps_set = .false.
      use_linear_k = .false.
c
      do node = 1, nnode
       x = node_coords(node,1)
       if( abs(x-5.0d0) .lt. 0.001d0 ) then
         pattern_values(node,4) = 0.0d00
       end if
      end do
      return
c
      end

