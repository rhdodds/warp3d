c     ****************************************************************
c     *                                                              *
c     *                      subroutine eqivld                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 11/17/2016 rhd             *
c     *                                                              *
c     *     compute the applied total load vector at end of current  *
c     *     step based on user specified definition of the load step.*
c     *     This includes specified nodal forces, specified          *
c     *     element tractions, pressures. The step increment may     *
c     *     also include specified nodal and element (constant)      *
c     *     temperatrures. This code enables approximate processing  *
c     *     of deformation dependent, specified element tractions,   *
c     *     pressures by using the current geometry.                 *
c     *     Convert vector of total applied nodal loads acting on    *
c     *     model at end of step to constraint compatable            *
c     *     global coordinates if necessary.                         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine eqivld( mf_tot, mf_nm1_tot, mf_ratio_change, step,
     &                   lodnum )
      use global_data ! old common.main
c
      use main_data, only : dtemp_nodes, dtemp_elems, temper_elems,
     &                      rload, dload, inverse_incidences,
     &                      load_pattern_factors, elem_eq_loads,
     &                      eq_node_force_indexes, eq_node_forces,
     &                      node_load_defs, crdmap, cnstrn_in,
     &                      temper_nodes_ref, total_user_nodal_forces,
     &                      user_cnstrn_stp_factors
c
      implicit integer (a-z)
c
c
c      step       - load step number for which loads are required
c      lodnum     - loading number for the nonlinear loading
c                   condition (loading patterns and nonlinear
c                   loading definitions are stored together)\
c      numlod     - total number of loading patterns and nonlinear
c                   loading conditions defined by user
c      mf_tot
c      mf_nm1_tot         multipliers and flags used for displaement
c      mf_ratio_change    extrapolation and adaptive sub-stepping
c                         if req'd
c      rload      - total of all user specified nodal forces and
c                   equivalent nodal forces from specified tractions-
c                   pressures. values here are based on accumulated
c                   multipliers of step loading patterns (as adjusted
c                   by global load reductions during damage-crack
c                   growth) + the user specified patterns-multipliers
c                   for this step.
c      rload_nm1  - rload at start of this load step.
c      load_pattern_factors - col. 2 stores user specified pattern
c                   multipliers for this step. col. 1 stores accumulated
c                   pattern multipliers prior to this step (as adjusted
c                   by any global loading reductions)
c
c      total_user_nodal_forces - exists only if there are loading
c                   patterns that reference the user_routine. contains
c                   the accumulated increments of nodal forces supplied
c                   by the user_routine thru step n. Updated here to
c                   n+1. More than 1 load pattern can employ the
c                   user_routine. The total_user_nodal_forces are
c                   saved across restarts. user_routines may
c                   provide any spatial or temporal varition of
c                   incremental (nodal) forces and nodal temperatures.
c
      double precision ::
     &  mf, mf_nm1, mf_tot, mf_nm1_tot, mf_ratio,
     &  step_factor, total_factor, dummy
      logical :: mf_ratio_change, user_mf_ratio_change
      character(len=80) :: user_file_name
      double precision, parameter :: zero = 0.0d00
      logical, parameter :: debug = .false.
c
c        1) initialization:
c
      dtemp_nodes(1:nonode) = zero
      dtemp_elems(1:noelem) = zero
      temperatures          = .false.
      load_pattern_factors(1:mxlc,2) = zero
c
      if( allocated( total_user_nodal_forces ) ) then
       rload(1:nodof) = total_user_nodal_forces(1:nodof)
      else
       rload(1:nodof)= zero
      end if
c
      call mem_allocate( 21 )
c
c        2) resolve the various multiplication factors and logical flag
c           used to support possible displacement extrapolation and
c           adaptive sub-incrementation during the load step. If current
c           step is not proportional to last step, we likely will
c           turn off displacement extrapolation for the next step
c           and use the linear stiffness on iteration 1 for next
c           step.  This can markedly improve convergence rates.
c
c
      call eqivld_resolve_extrapolate(
     &   mf_tot, mf_nm1_tot, mf_ratio_change, step, lodnum, numlod,
     &   user_cnstrn_stp_factors, debug, out, nodof, cstmap,
     &   cnstrn_in )
c
c        3) loop over all loading conditions (patterns). if that pattern
c           contributes to the incremental load for this step,
c           set the step increment factor for the pattern.
c
      load_pattern_factors(1:mxlc,2) = zero
c
      do ldcond = 1, numlod
         call get_step_factor( ldcond, step, mf )
         if ( mf .ne. zero ) mf_tot = mf
         if ( mf .eq. zero ) cycle
         load_pattern_factors(ldcond,2) = load_pattern_factors(ldcond,2)
     &                                    + mf
      end do
c
      if ( debug ) then
        write(out,*) '.. updated pattern factors. incr and total...'
        do ldcond = 1, numlod
          write(out,*) '  ',lodnam(ldcond),'  ',
     &                    load_pattern_factors(ldcond,2),
     &                    load_pattern_factors(ldcond,1)
        end do
      end if
c
c        4) loop over all loading conditions (patterns). compute the
c           total multiplier for that pattern thru the end of this step.
c           get applied nodal forces and equiv. nodal forces from
c           specified element loads for the total pattern value.
c           get the incremental nodal temps and element temps for
c           user specified step factors. at end of loop, rload
c           contains the total applied nodal loads + equiv. nodal
c           loads. note that contributions to rload from specified
c           element tractions-pressures are based on the current
c           deformed geometry for large displacement analysis.
c           once rload_nm1 is subtracted from rload, the incremental
c           nodal loads (dload) for the step reflect an approximate
c           treatement of deformation dependent imposed element
c           tractions-pressures.
c
c           for user_routine we keep an accumulation of incremental
c           nodal forces from the user routine to include in rload.
c           by keeping the acummulation vector (and across restarts),
c           the user_routine can truly have varying incremental forces
c           over time, temperature, etc.
c
c
      do ldcond = 1, numlod
         how_defined = node_load_defs(ldcond)%how_defined
         step_factor  = load_pattern_factors(ldcond,2)
         total_factor = load_pattern_factors(ldcond,1) + step_factor
         select case ( how_defined )
         case( 0 )  ! user defined list in input file
            if ( total_factor .ne. zero ) then
               call nodal_loads( total_factor, ldcond, dstmap,
     &                           rload, debug )
               call driv_eload( ldcond, total_factor, rload )
            end if
            if ( step_factor .ne. zero )
     &         call nodal_incr_temps( step_factor, ldcond,
     &                                temperatures, debug )
         case ( 1 ) ! user subroutine
            if( .not. allocated( total_user_nodal_forces ) ) then
              allocate( total_user_nodal_forces(nodof) )
              total_user_nodal_forces(1:nodof) = zero
            end if
            user_mf_ratio_change = .false.
            user_file_name = node_load_defs(ldcond)%user_file_name
            call eqivld_drive_user_nodal( step,
     &            step_factor, ldcond, temper_nodes_ref,
     &            total_user_nodal_forces, rload, dtemp_nodes,
     &            crdmap, user_mf_ratio_change, user_file_name,
     &            debug  )
            if( user_mf_ratio_change ) mf_ratio_change = .true.
         case default
            write(iout,*) '>>> invalid loading condition type'
            write(iout,*) '    in routine eqivld'
            call die_abort
         end select
      end do
c
c        5) rotate end of step, total applied nodal forces + equiv.
c           nodal forces into constraint compatible coordinates (alternate
c           corodinate systems at nodes specified to impose constraints).
c
      call rotate_loads( rload, debug )
c
c
c        6) convert pointer-vector data structure for element
c           equivloads to a packed, single vector of values with
c           an integer indexing vector. makes mpi communication
c           much simpler. throw away point-vector structure as
c           we copy & compact.
c
      call mem_allocate( 22 )
      count = 1
      do i = 1, noelem
       ncols = elem_eq_loads(i)%ncols
       if( ncols .gt. 0 ) then
         eq_node_force_indexes(i) = count
         call vec_ops( eq_node_forces(count),
     &                 elem_eq_loads(i)%forces(1,1), dummy, 3*ncols,
     &                 5 )
         count = count + 3*ncols
         deallocate( elem_eq_loads(i)%forces )
       end if
      end do
      deallocate( elem_eq_loads )

      if ( debug ) write(*,*) '<<<< leaving eqivld.f'
c
      return
c
      end
c     ****************************************************************
c     *                                                              *
c     *                subroutine eqivld_drive_user_nodal            *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *                 last modified : 11/19/2016 rhd               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine eqivld_drive_user_nodal(
     &  step, step_factor, ldcond, initial_temps,
     &  tunf, rload, warp_incr_temperatures, crdmap,
     &  user_mf_ratio_change, user_file_name, debug )
      use global_data ! old common.main
      implicit none
c
c         parameters
c
      integer :: step, ldcond, crdmap(*)
      logical :: debug, user_mf_ratio_change
      double precision ::
     &  step_factor, tunf(nodof), rload(nodof),
     &  warp_incr_temperatures(*), initial_temps(*)
      character(len=80) :: user_file_name
c
c        local variables
c
      integer :: num_model_nodes, step_np1, kout, node, dloc, i
      character(len=8) :: load_name
      logical :: forces_set, temps_set, process
      double precision ::
     & time_n, dtime, f1, f2, f3
      double precision, parameter :: zero = 0.0d0
c
      double precision,
     &       dimension(:,:), allocatable :: incr_values,
     &                                      node_coords
c
c        set local copies of key variables just in case the user rouitne
c        decides to change them
c
      if( debug ) write(out,*) '... inside eqivld_drive_user_nodal ..'
      num_model_nodes = nonode          ! common.main
      step_np1 = step
      time_n =  total_model_time        ! common.main
      dtime = dt
      load_name = lodnam(ldcond)        ! common.main
      kout = out
      forces_set = .false.
      temps_set  = .false.
      process = step_factor .ne. zero
      if( .not. process ) return
c
c        local work tables for user routine to store load values
c        in flat data structures. build flat structure for nodal
c        coordinates (undeformed). Zero incremental loads. Columns
c        are delta Fx, Fy, Fz, temp
c
      allocate( incr_values(nonode,4) )
      allocate( node_coords(nonode,3) )
c
      do node = 1, nonode
        node_coords(node,1) = c(crdmap(node)+0)
        node_coords(node,2) = c(crdmap(node)+1)
        node_coords(node,3) = c(crdmap(node)+2)
        incr_values(node,1:4) = zero
      end do
c
      if( debug ) then
         write(out,9010) nonode, step, time_n, dt, step_factor
         write(out,*) '... load_name: ', load_name
      end if
c
c        invoke user routine. it sets two logical flags
c        about what it did. returns nodal forces and temperatures
c        interpreted here as patttern values to be scaled by user factor
c        in nonlinear step definition
c
      call user_nodal_loads(
     &  load_name, user_file_name, incr_values, initial_temps,
     &  num_model_nodes, step_np1, time_n, dtime,
     &  node_coords, forces_set, temps_set, user_mf_ratio_change,
     &  kout )
c
      if( debug ) write(out,*) '... forces_set, temps_set',
     &                        forces_set, temps_set
c
c        tunf -> total_user_nodal_forces (shorthand). rload was
c        preloaded with tunf at step n before call. include
c        incremental forces from user routine in both.
c
c        this scheme allows user routine to provide non-proportional
c        loads across load(time) steps.
c
      if( forces_set ) then
       do node = 1, nonode
          dloc = dstmap(node)
          f1 = incr_values(node,1) * step_factor
          f2 = incr_values(node,2) * step_factor
          f3 = incr_values(node,3) * step_factor
          tunf(dloc+0)  = tunf(dloc+0) + f1
          tunf(dloc+1)  = tunf(dloc+1) + f2
          tunf(dloc+2)  = tunf(dloc+2) + f3
          rload(dloc+0) = rload(dloc+0) + f1
          rload(dloc+1) = rload(dloc+1) + f2
          rload(dloc+2) = rload(dloc+2) + f3
       end do
      end if
c
      if( temps_set ) then
        temperatures = .true.  ! common.main
        do node = 1, nonode
          warp_incr_temperatures(node) = warp_incr_temperatures(node) +
     &                              incr_values(node,4) * step_factor
        end do
      end if
c
c        all done. release work arrays
c
      deallocate( incr_values,  node_coords )
c
      return
c
 9010 format(' ... nonode, step, time_n, dt, step_factor: ',i10,i7,
     & 3e14.6)
c
      end

c     ****************************************************************
c     *                                                              *
c     *                subroutine get_step_factor                    *
c     *                                                              *
c     *                    written by : rd                           *
c     *                                                              *
c     *                 last modified : 05/1/02                      *
c     *                                                              *
c     *     for a specified load step number and a loading pattern   *
c     *     number, return the user specified multiplier (could be   *
c     *     = 0)                                                     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine get_step_factor( ldcond, step, mf )
      use main_data, only : step_load_data
      implicit none
c
c         ldcond   -  one of possibly many pattern load condition nos.
c         step     -  load step being processed
c         mf       -  the current nonlinear load step may reference
c                     the loading pattern ldcond. if so, return
c                     the user specified multipler for the step.
c                     otherwise, mf = 0.0
c
      integer ldcond, step
      double precision
     &  mf
c
      integer i, num_patt
      double precision
     &  zero
       data zero / 0.0d00 /
c
c           step_load_data lists the user specified patterns and multipliers
c           for all load steps.
c
      mf       = zero
      num_patt = step_load_data(step)%num_load_patterns
      if ( num_patt .eq. 0 ) return
c
      do i = 1, num_patt
        if ( step_load_data(step)%load_patt_num(i) .eq. ldcond ) then
           mf = step_load_data(step)%load_patt_factor(i)
        end if
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *              subroutine eqivld_resolve_extrapolate           *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 09/10/12                   *
c     *                                                              *
c     *     resolve values for the various multipliers, flags used   *
c     *     to support displacement extrapolation at the start of    *
c     *     this load step                                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine eqivld_resolve_extrapolate(
     &     mf_tot, mf_nm1_tot, mf_ratio_change, now_step,
     &     nonlinear_loading,
     &     num_loading_conds, user_cnstrn_stp_factors, debug, iout,
     &     nodof, cstmap, cnstrn_in )
      implicit none
c
      integer now_step, nonlinear_loading, num_loading_conds,
     &        iout, nodof, cstmap(*)
      double precision
     &  mf, mf_nm1, mf_tot, mf_nm1_tot, mf_ratio, cnstrn_in(*)
      real user_cnstrn_stp_factors(*)
      logical debug, mf_ratio_change

      integer ldcond, i, exit_point
      double precision
     &   zero, one, cfnm1, cfn, con_ratio, toler
      logical  non_zero_mf, non_zero_mf_nm1, all_zero_cons,
     &         loading_patterns_exist
      data zero, one, toler / 0.0d00, 1.0d00, 1.0d-05 /
c
c           mf_tot and mf_nm1_tot are the variables used
c           by extrapolation and the adaptive step sizing algorithm.
c
c           Simple case: step n has 1 loading pattern; step n-1 has
c           1 (same) loading pattern. Multiplers may be different.
c           Set the multipler ratio between n-1 and n for subsequent
c           extrapolation.
c
c           More complex: n-1 and/or n have more than one pattern
c           in the definiton. If the patterns are same and mutipliers
c           from n-1 to n all have same ratio, the incremental loading
c           for step n is proportional. Extrpaolation will be fine
c           using just a single scaling factor.
c
c           Extrapolation will likely not work well if the
c           patterns at n-1 and n are different or they are the
c           same but multiplier ratios are not all the same, i.e.,
c           incremental loading from n-1 to n is non-proportional.
c           We set flag for eventual warning message.
c
      if( debug )
     &   write(iout,9100) now_step, num_loading_conds, nonlinear_loading
c
      mf_tot          = zero
      mf_nm1_tot      = zero
      mf_ratio        = zero
      mf_ratio_change = .false.
      non_zero_mf     = .false.
      non_zero_mf_nm1 = .false.
c
      do ldcond = 1, num_loading_conds ! defined for model
         if( ldcond .eq. nonlinear_loading ) cycle ! non. lin. loading
         call get_step_factor( ldcond, now_step, mf )
         if( debug ) write(iout,9110) ldcond, mf
         if( mf .ne. zero ) then
            mf_tot = mf
            non_zero_mf =.true.
         end if
         if( now_step .eq. 1 ) then
           mf_ratio = one
           cycle
         end if
         call get_step_factor( ldcond, now_step-1, mf_nm1 )
         if ( debug ) write(iout,9120) mf_nm1
         if( abs(mf_nm1) .gt. zero ) then
              mf_nm1_tot = mf_nm1
              non_zero_mf_nm1 = .true.
         end if
c
c                   if this loading pattern was not present at
c                   n-1, or has been dropped at n,
c                   or is in n-1 and n but the multiplier is different,
c                   set logical flag that will eventually
c                   cause a warning message to be issued about
c                   non-proportional incremental loads with
c                   extrapolation.
c
          if( mf_nm1 .ne. zero ) then ! pattern present at n-1
c
c                     1) pattern present or not at n?
c
               if( mf .eq. zero ) then
c
c                         a) pattern not present at n
c
                  mf_ratio_change = .true.
                  if( debug ) write(iout,9130)
               else
c
c                         b) pattern present at n-1 & n, check ratio
c                            of multipliers.
c
                  if( mf_ratio .eq. zero ) then
                     mf_ratio = mf / mf_nm1
                  else
                     if( mf / mf_nm1 .ne. mf_ratio ) then
c
c                               multiplier changed from n-1 to n
c
                        mf_ratio_change = .true.
                        if( debug ) write(iout,9150)
                     end if
                  end if
               end if
            else
c
c                     2) pattern present at n but not n-1
c
               if( mf .ne. zero ) then
c
c                            pattern appears at n but not n-1
c
                  mf_ratio_change = .true.
                  if( debug ) write(iout,9140)
               end if
c
            end if
c

      end do
c
c                   check the user specified constraints
c                   for this step (n) vs prior step (n-1) against the
c                   changes in other loading patterns specifed for the
c                   two stpes. we're looking to see of incremental load -
c                   including - constraints is proportional from n-1
c                   to n.
c
c                   Eliminate simple cases to start:
c
c                    - if step = 1 nothing to check further, return
c
c                    - if mf_ratio_change = .true. from above, return.
c                      the loading patterns (w/ multipliers) for step n
c                      are non-proportional w/o consideration of
c                      constraints. no further checks needed. return
c
c                    - if imposed constraints are all = 0.0, then
c                      they have no effect on proportinality
c                      of the incremental loading. just return.
c
c                    - we can still have no loading patterns,
c                      just constraints only loading
c
      if( now_step .eq. 1 ) then
         exit_point = 1
         go to 8000
      end if
      if( mf_ratio_change ) then
          mf_tot = one
          mf_nm1_tot = one
          exit_point = 1
         go to 8000
      end if
c
      all_zero_cons = .true.
      do i = 1, nodof
        if( cstmap(i) .eq. 0 ) cycle
        if( abs( cnstrn_in(i) ) .gt. zero ) all_zero_cons = .false.
      end do
      if( all_zero_cons ) then
            exit_point = 2
            go to 8000
      end if
c
c                   we have non-zero constraints at current step.
c                   we have to assume the previous constraints
c                   are same values.
c
c                   do we also have loading patterns included at
c                   n-1 and n? If not, just simple constraint
c                   loading of model and proportionality
c                   not an issue.
c
c                   if yes, they have already been
c                   verified proportional above. Does the specified
c                   constraint multiplier at n-1 and n match
c                   proportionailty of the loading patterns so that
c                   the incremental loading "vector" (including
c                   constraints) is in same "direction"?
c
      loading_patterns_exist =  non_zero_mf .and. non_zero_mf_nm1
      cfnm1 = user_cnstrn_stp_factors(now_step-1)
      cfn   = user_cnstrn_stp_factors(now_step)
      con_ratio = zero
      if( abs(cfnm1) .gt. toler ) con_ratio = cfn / cfnm1
c
      if( loading_patterns_exist ) then
         if( abs(con_ratio - mf_ratio) .le. toler ) then
            exit_point = 3
            go to 8000
         end if
         mf_ratio_change = .true.
         exit_point = 4
         go to 8000
      end if
c
c                constraints is the only patterm listed
c                for steps n-1 and n. just set multipliers
c                for use by mnralg.
c
       mf_tot     = cfn
       mf_nm1_tot = cfnm1
       exit_point = 5
       go to 8000
c
 8000  continue
       if( debug ) then
         write(iout,9200) exit_point, now_step,
     &                    mf_tot, mf_nm1_tot,  mf_ratio_change
         if( exit_point .ge. 2 ) write(iout,9210)
     &              all_zero_cons
         if( exit_point .gt. 2 )
     &           write(iout,9220) cfnm1, cfn, con_ratio, non_zero_mf,
     &               non_zero_mf_nm1
         if( exit_point .eq. 5 )
     &           write(iout,9230)  mf_tot, mf_nm1_tot
         write(iout,*) ' '
      end if
c
      return
 9200 format(/,1x,'... leave eqivld_resolve_extrapolation ...',
     & /,10x,'exit_point, now_step: ',2i10,
     & /,10x,'mf_tot, mf_nm1_tot,  mf_ratio_change: ',2f10.4,l5)
 9210 format(10x,'all_zero_cons: ',l5)
 9220 format(/,10x,'cfnm1, cfn, con_ratio: ',3f10.4,
     & /,10x,'non_zero_mf, non_zero_mf_nm1: ', 3l2 )
 9230 format(/,10x,'mf_tot, mf_nm1_tot (just cons): ',2f10.4 )
 9100 format( /,1x,'... enter eqivld_resolve_extrapolation ...',
     & /,5x,'load step, no. model loadings, nonlinear loading no: ',
     & i7,2i5)
 9110 format(10x,'.. load pattern, returned mf: ',i3,f10.6)
 9120 format(20x,'multiplier for prior step: ',f10.6)
 9130 format(20x,'pattern present in prior step but not this one...')
 9140 format(20x,'pattern not present in prior step...')
 9150 format(20x,'pattern multiplers not-proportional to prior step...')

      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine nodal_incr_temps                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 05/1/02                    *
c     *                                                              *
c     *     extract user-specified, incremental nodal temperatures   *
c     *     for a pattern, scale by step factor                      *
c     *     and assemble into structure temperture vector            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine nodal_incr_temps( mf, ldcond, temperatures, debug )
      use main_data, only : dtemp_nodes, node_load_defs
      implicit integer (a-z)
c
      double precision
     &  mf, zero
      logical debug, temperatures
c
      real node_vals(10)
      integer, dimension(:,:), pointer :: node_lod_data
      data zero / 0.0 /
c
      if ( debug ) write(*,*) '>> in nodal_incr_temps'
c
c               check for nodal loads and process if they are defined.
c
      node_count = node_load_defs(ldcond)%node_count
      if ( node_count .eq. 0 ) return
      node_lod_data => node_load_defs(ldcond)%nodal_loads
c
      do col = 1, node_count
c
c               set the node at which loads occur and the index
c               into the table containing the magnitudes of the loads
c               at each dof of the node
c
         node       = node_lod_data(col,1)
         loddat_col = node_lod_data(col,2)
c
c               extract pattern load on node, multiply by
c               step factor. assemble into in the structure size vector.
c
         call loddat_ops( 3, node_vals, loddat_col )
         dtemp_nodes(node) = dtemp_nodes(node) + node_vals(4)*mf
         if ( abs(dtemp_nodes(node)) .ne. zero ) temperatures = .true.
c
      end do
c
      if ( debug ) write(*,*) '<< leave nodal_incr_temps'
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine nodal_loads                  *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 05/1/02 rhd                *
c     *                                                              *
c     *     extract nodal loads for a pattern, scale by step factor  *
c     *     and assemble into structure load vector                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine nodal_loads( mf, ldcond, dstmap, rload, debug )
      use main_data, only : node_load_defs
      implicit integer (a-z)
c
      double precision
     &  mf, rload(*)
      integer dstmap(*)
      logical debug
c
      real node_vals(10)
      integer, dimension(:,:), pointer :: node_lod_data
c
      if ( debug ) write(*,*) '>> in nodal_loads'
c
c               check for nodal loads and process if they are defined.
c               this code works only for 3 dof per node.
c
      node_count = node_load_defs(ldcond)%node_count
      if ( node_count .eq. 0 ) return
      node_lod_data => node_load_defs(ldcond)%nodal_loads
c
      do col = 1, node_count
c
c               set the node at which loads occur and the index
c               into the table containing the magnitudes of the loads
c               at each dof of the node
c
         node       = node_lod_data(col,1)
         loddat_col = node_lod_data(col,2)
c
c               extract pattern load on node, multiply by
c               step factor. assemble into in the structure size vector.
c
         dloc = dstmap(node)
         call loddat_ops( 3, node_vals, loddat_col )
         rload(dloc+0) = rload(dloc+0) + node_vals(1)*mf
         rload(dloc+1) = rload(dloc+1) + node_vals(2)*mf
         rload(dloc+2) = rload(dloc+2) + node_vals(3)*mf
c
      end do
c
      if ( debug ) write(*,*) '<< leave nodal_loads'
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine rotate_loads                     *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 8/21/2016 rhd              *
c     *                                                              *
c     *     this subroutine rotates the global nodal load vector     *
c     *     into constraint compatable global coordinates.           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rotate_loads( rload, debug )
      use global_data ! old common.main
c
      use main_data, only : trn, trnmat, inverse_incidences
c
      implicit none
c
      double precision :: rload(*)
      logical :: debug
c
      integer :: node, dptr, ndof, i
      double precision ::
     &  ndlod(mxvl,mxndof), zero, trnmte(mxvl,mxedof,mxndof)
      logical :: loads_found, trne(mxvl,mxndel)
      data zero / 0.0d0 /
c
      if( debug ) write(*,*) '>> in rotate_loads -- rotating loads'
c
c           loop over all the nodes in the structure. If the node
c           has a local coordinate system for definition of
c           constraints, then check for existence of applied forces
c           at the node. rotate the forces into constraint compatible
c           coordinates.
c
c           Note that we copy the dload entry for a given node into
c           ndlod, which has a size of (mxvl,nxndof), even though we
c           only ever use a one in the first index.  This is because
c           trnvec usually works on a block structure, and thus requires
c           the vector to be vectorizable. same for trnmte and
c           trne (local) arrays.
c
c           the value of ndof = 3 is hard coded here ....
c
      do node = 1, nonode
       if( .not. trn(node) ) cycle
       dptr        = dstmap(node)
       ndof        = iprops(4,inverse_incidences(node)%element_list(1))
       ndlod(1,1)  = rload(dptr+0)
       ndlod(1,2)  = rload(dptr+1)
       ndlod(1,3)  = rload(dptr+2)
       loads_found = abs(ndlod(1,1)) + abs(ndlod(1,2)) +
     &               abs(ndlod(1,3)) .ne. zero
       if( debug ) then
         write(*,*) '>> node, dptr, ndof, loads_found: ',node, dptr,
     &    ndof, loads_found
         if( loads_found ) write(*,9000) ndlod(1,1:3)
       end if
       if( .not. loads_found ) cycle
c
c           extract trans. matrix to constraint compatable
c           global coordinates for this node and transform the
c           nodal incremental load vector to constraint compatable
c           form.
c
       trne(1,1)     = .true.
       trnmte(1,1,1) = trnmat(node)%mat(1,1)
       trnmte(1,1,2) = trnmat(node)%mat(1,2)
       trnmte(1,1,3) = trnmat(node)%mat(1,3)
       trnmte(1,2,1) = trnmat(node)%mat(2,1)
       trnmte(1,2,2) = trnmat(node)%mat(2,2)
       trnmte(1,2,3) = trnmat(node)%mat(2,3)
       trnmte(1,3,1) = trnmat(node)%mat(3,1)
       trnmte(1,3,2) = trnmat(node)%mat(3,2)
       trnmte(1,3,3) = trnmat(node)%mat(3,3)
       call trnvec( ndlod, trnmte, trne, ndof, 1, 1, 1 )
       if ( debug ) write(*,9000) ndlod(1,1:3)
       rload(dptr+0) = ndlod(1,1)
       rload(dptr+1) = ndlod(1,2)
       rload(dptr+2) = ndlod(1,3)
c
      end do
c
      return
c
 9000 format('>> ndlod: ',3f15.6)
c
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine eqiv_out_patters                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 5/21/02 rhd                *
c     *                                                              *
c     *     at completion of a load step, output a summary of        *
c     *     loading pattern multipliers acting on model. this aids   *
c     *     analysis interpretation when gloabl load step reductions *
c     *     occur during crack growth and other types of analyses    *
c     *                                                              *
c     ****************************************************************
c
      subroutine eqiv_out_patterns
      use global_data ! old common.main
c
      use main_data,  only : load_pattern_factors, output_packets,
     &                       packet_file_no, actual_cnstrn_stp_factors
c
      implicit integer (a-z)
      double precision
     &     step_factor, total_factor, zero, sum
      character(len=12) pattern_names(3)
      real tfacts(3)
      data zero / 0.0 /
c
      write(out,9000) ltmstp+1
c
c           sum up the constraint factors through the current load
c           step
c
      sum = zero
      do i = 1, ltmstp+1
        sum = sum + actual_cnstrn_stp_factors(i)
      end do
c
c           write out the loading patterns and accumulated multipliers
c           3 to the line. insert constraints multiplier at end.
c
      line_count = 0
      pattern_count = 0
      do i = 1, numlod
        if ( line_count .eq. 3 ) then
          write(out,9010) pattern_names(1), tfacts(1),
     &                    pattern_names(2), tfacts(2),
     &                    pattern_names(3), tfacts(3)
          line_count = 0
        end if
        step_factor  = load_pattern_factors(i,2)
        total_factor = load_pattern_factors(i,1)
        if ( total_factor .ne. zero ) then
          line_count                = line_count + 1
          pattern_names(line_count) = lodnam(i)
          tfacts(line_count)        = total_factor
          pattern_count = pattern_count + 1
        end if
      end do
c
      line_count = line_count + 1
      pattern_count = pattern_count + 1
      pattern_names(line_count) = 'constraints'
      tfacts(line_count) = sum
      write(out,9010) (pattern_names(i), tfacts(i), i=1,line_count)
      write(out,*) ' '
c
c
c           write binary packets for the same information if packet
c           output is in effect.
c
      if( .not. output_packets ) return
      write(packet_file_no) 24, pattern_count, ltmstp+1, 0
      do i = 1, numlod
        total_factor = load_pattern_factors(i,1)
        if ( total_factor .ne. zero ) then
           pattern_names(1) = lodnam(i)
           write(packet_file_no) pattern_names(1), total_factor
        end if
      end do
      pattern_names(1) = 'constraints'
      write(packet_file_no) pattern_names(1), sum


      return
 9000 format(/1x,'>> total applied load pattern factors through step: ',
     & i7 )
 9010 format(6x,'> ',3(a12,' ',f10.3,2x))
      end
