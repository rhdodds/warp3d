c
c     ****************************************************************
c     *                                                              *
c     *              subroutine check_for_step_reduction             *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 2/12/2018 rhd              *
c     *                                                              *
c     *        This routine branches on the type of load step        *
c     *        reduction to check before computing the step.         *
c     *        This allows algorithms which require finer            *
c     *        granularity of step size to request it.               *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine check_for_step_reduction( temp_load_fact, mf, mf_nm1 )
      use global_data ! old common.main
      use damage_data
      implicit none
c
      double precision ::  temp_load_fact, mf, mf_nm1
c
      logical, parameter ::  debug = .false.
      double precision, parameter ::  one = 1.0d0, zero = 1.0d0
c
      if ( debug ) then
         write(*,*) ' > entered check_for_step_reduction'
         write(*,*) '   load_size_...: ',load_size_control_crk_grth
         write(*,*) '   crack_growth_type: ',crack_growth_type
      end if
c
      temp_load_fact = one
c
c        1) check for permanent step size reductions:
c
c           - if we are using crack growth by node release, and
c             automatic load step sizing is enabled, then check the
c             number of load steps between crack growth increments.
c             If this is less than the user-defined value of
c             min_steps_for_release, then cut the loading in half.
c
c           - if we are using crack growth by element extinction
c             with the gurson material model, and automatic load step
c             sizing is enabled, then check the change in porosity
c             for the last load step. subsequent load steps
c             can be increased or decreased in size.
c
c           - for cohesive elements, the load step sizes can
c             be adjusted up and down using same algorithms as
c             the gurson elements
c
      if ( load_size_control_crk_grth ) then
         if ( crack_growth_type .eq. 1 ) call gurson_cut_step ( debug )
         if ( crack_growth_type .eq. 2 ) then
            if (const_front) then
               call ctoa_cut_step_const ( debug )
            else
               call ctoa_cut_step ( debug )
            endif
         endif
         if ( crack_growth_type .eq. 3 ) call smcs_cut_step ( debug )
         if ( crack_growth_type .eq. 4 ) call cohes_cut_step ( debug )
      endif
c
c
c        2) check for temporary step size reductions
c
c           - overshoot control for crack growth routines; reduces
c             step size to prevent large overshoots of criterion
c             for crack growth. Branch on type of crack growth if
c             overshoot control is enabled.
c
      if ( overshoot_control_crk_grth ) then
         if ( crack_growth_type .eq. 2 ) then
c
c                 if we are measuring a defined distance for the CTOA,
c                 call a different routine
c
            if (const_front) then
               call over_CTOA_const ( temp_load_fact, mf, mf_nm1,
     &              debug )
            else
               call overshoot_CTOA ( temp_load_fact, mf, mf_nm1,
     &              debug )
            endif
         endif
      endif
c
c
 9999 continue
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine modify_load                  *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 2/12/2018 rhd              *
c     *                                                              *
c     *        Use the load factor calculated by the load step       *
c     *        reduction criterion (from damage, crack growth, etc.) *
c     *        to reduce the incremental loads and constraints.      *
c     *        This redefines the user's specified incremental loads *
c     *        for the step - permanently or transiently. We change  *
c     *        the applied nodal loads, equiv. element loads, temps  *
c     *        for the step but not the tables that hold the user's  *
c     *        step definitions.                                     *
c     *                                                              *
c     ****************************************************************
c
      subroutine modify_load ( temp_load_fact, mf, mf_nm1, nowstep )
      use global_data ! old common.main
c
      use main_data,  only : cnstrn, dload, rload, dtemp_nodes,
     &                       load_pattern_factors, rload_nm1,
     &                       dtemp_elems, actual_cnstrn_stp_factors,
     &                       user_cnstrn_stp_factors
      use damage_data
c
      implicit none
c
      integer :: nowstep
      double precision :: temp_load_fact, mf, mf_nm1
c
      integer :: i
      double precision :: step_factor, total_factor
      double precision, parameter :: one = 1.0d0
      logical, parameter :: debug = .false.
c
      if ( debug ) write (out,9000)
c
c          find the new load factor for the step, and update the old
c          load factor.  There are two factors to include in the new
c          load factor: a permanent part (perm_load_fact) for permanent
c          load reductions, and a one-step-only reduction factor
c          (temp_load_fact).
c
      old_load_fact     = control_load_fact
      control_load_fact = temp_load_fact * perm_load_fact
c
c          if load factor is too small, set it to mimimum value.
c
      if ( control_load_fact .lt. min_load_fact )
     &     control_load_fact = min_load_fact
      if ( debug ) write(out,9010) control_load_fact, old_load_fact
c
c          scale load factors for extrapolation algorithm
c
      if ( debug ) write(out,9020) mf, mf_nm1
      mf     = mf * control_load_fact
      mf_nm1 = mf_nm1 * old_load_fact
      if ( debug ) write(out,9030) mf, mf_nm1
c
c          skip scaling the loads-temps if control_load_fact is 1. just set
c          dload for step, update the global load vector (which has
c          inertia and crack growth loads added later). update the
c          actual, accumulated multipliers for each loading pattern.
c          initialize rload n-1 for next step. The rload and rload_nm1
c          vectors contain only the user specified nodal forces and
c          equivalent nodal forces from specified element pressures,
c          tractions.
c
      if ( control_load_fact .eq. one ) then
        do i = 1, nodof
          dload(i)     = rload(i) - rload_nm1(i)
          load(i)      = load(i) + dload(i)
          rload_nm1(i) = rload(i)
        end do
        actual_cnstrn_stp_factors(nowstep) =
     &            user_cnstrn_stp_factors(nowstep)
        do i = 1, numlod
         step_factor  = load_pattern_factors(i,2)
         total_factor = load_pattern_factors(i,1)
         load_pattern_factors(i,1) = total_factor + step_factor
        end do
        if ( debug ) then
          write(out,9040)
          do i = 1, numlod
            step_factor  = load_pattern_factors(i,2) * control_load_fact
            total_factor = load_pattern_factors(i,1)
            write(out,9050) lodnam(i),  step_factor, total_factor
          end do
        end if
        go to 9999
      end if
     &
c
c          scale all of the non-zero constraints and the incremental
c          loading by the scaling factor from load step size control.
c          see comments above for vectors. We must also scale the
c          user specified increment of nodal and element temperatures
c          for the step.
c
      do i = 1, nodof
        if ( cstmap(i) .ne. 0 ) cnstrn(i) = cnstrn(i) *
     &                                      control_load_fact
        dload(i)     = ( rload(i) - rload_nm1(i) ) * control_load_fact
        load(i)      = load(i) + dload(i)
        rload(i)     = rload_nm1(i) + dload(i)
        rload_nm1(i) = rload(i)
      end do
      actual_cnstrn_stp_factors(nowstep) =
     &       user_cnstrn_stp_factors(nowstep) * control_load_fact
      dtemp_nodes(1:nonode) = dtemp_nodes(1:nonode) * control_load_fact
      dtemp_elems(1:noelem) = dtemp_elems(1:noelem) * control_load_fact
      do i = 1, numlod
        step_factor  = load_pattern_factors(i,2) *  control_load_fact
        total_factor = load_pattern_factors(i,1)
        load_pattern_factors(i,1) = total_factor + step_factor
      end do
      if ( debug ) then
        write(out,9040) control_load_fact
        do i = 1, numlod
          step_factor  = load_pattern_factors(i,2) * control_load_fact
          total_factor = load_pattern_factors(i,1)
          write(out,9050) lodnam(i),  step_factor, total_factor
        end do
      end if
c
 9999 continue
      if ( debug ) write (out,*) '<<<<<<<< leaving modify load'
      return
c
 9000 format('>>> Entering modify_loads...')
 9010 format('>>> Load facts(new, old):',2e13.6)
 9020 format('    original mf, mf_nm1:',2e13.6)
 9030 format('    new mf, mf_nm1:',2e13.6)
 9040 format('... Current pattern loading factors. control factor: ',
     &      f8.3 )
 9050 format(10x, a8, 2f10.3 )
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine original_step_size                *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 2/12/2018                  *
c     *                                                              *
c     *        This routine restores the original loading if         *
c     *        the reduction algorithms reduced the step size        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine original_step_size ( mf, mf_nm1, step )
      use global_data ! old common.main
      use damage_data
      implicit none
c
      integer :: step
      double precision :: mf, mf_nm1
c
c          change back load multipliers. On step 1 the old
c          load factor is zero, so don't change mf_nm1.
c
      mf = mf / control_load_fact
      if( step .gt. 1 ) mf_nm1 = mf_nm1 / old_load_fact
c
      return
      end
