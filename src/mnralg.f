c     ****************************************************************
c     *                                                              *
c     *                      subroutine mnralg                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified: 5/16/2018 rhd               *
c     *                                                              *
c     *     supervises advancing the solution from                   *
c     *     step n to n+1 using a newton iteration process.          *
c     *     a solver (direct or iterative)                           *
c     *     resolves the linear equations at each newton iteration.  *
c     *     extrapolation of displacements to accelerate convergence *
c     *     and adaptive solution options are available for static   *
c     *     and dynamic analyses.                                    *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine mnralg( mf, mf_nm1, mf_ratio_change, step, ldnum )
      use global_data ! old common.main
c
      use main_data, only : mdiag, pbar, cnstrn, dload,
     &     max_loddat_blks, convergence_history, umat_used,
     &     divergence_check, asymmetric_assembly, creep_model_used,
     &     non_zero_imposed_du, extrapolated_du, extrapolate,
     &     extrap_off_next_step, line_search, ls_details,
     &     ls_min_step_length, ls_max_step_length, ls_rho,
     &     ls_slack_tol
      use adaptive_steps, only : adapt_result, adapt_disp_fact,
     &                           adapt_load_fact
      use hypre_parameters, only : hyp_trigger_step
      use distributed_stiffness_data, only: parallel_assembly_used,
     &                               distributed_stiffness_used,
     &                               parallel_assembly_allowed
      use performance_data
      use mod_mpc, only: mpcs_exist, tied_con_mpcs_constructed
      use stiffness_data, only :  total_lagrange_forces,
     &                            d_lagrange_forces,
     &                            i_lagrange_forces
c
      implicit integer (a-z)
      double precision :: mf, mf_nm1, dt_original
      logical :: mf_ratio_change
      double precision :: sumrr, sumdd
c
c           local variables
c
      double precision ::
     &  mgres1, magdu1, zero, one, mgload, ext_tol,
     &  scaling_load, mag
      double precision, allocatable, dimension(:) ::u_n_local
c
      logical :: dynamic, material_cut_step, emit_extrap_msg,
     & emit_forced_linear_k_for_step, ls_request_adaptive,
     & diverging_flag, user_extrapolation_on, user_extrapolation_off
c
      character (len=1) :: dums
      logical :: local_debug, first_solve, first_subinc,
     &        hypre_solver, cnverg, adaptive, local_mf_ratio_change,
     &        check_crk_growth
c
      type :: info_mnralg
        logical :: adaptive_used
        logical :: step_converged
        integer :: adapt_substeps
        integer :: iters_at_convergence
        integer :: adapt_iters_at_convergence
      end type
      type( info_mnralg ) info
c
      data zero, one, ext_tol, local_debug / 0.0d00, 1.0d00,
     &      1.0d-20, .false. /
c
c          predct:          .true. if user requested extrapolate on
c          adaptive_flag    .true. if user requested adaptive solution
c          first_subinc     .true. if we are executing the the first
c                           subincrement after adaptive has reduced
c                           or restired the step size
c          mf_ratio_change  .true. if the equivloads process says this
c                           step load increment is NOT proportional to
c                           step n laod increment
c
      if( local_debug ) then
         write(out,*) " "
         write(out,*) "..... entered mnralg ....."
      end if
c
      now_step = step ! just for protection
      iout = out
c
c          store model displacement vector at start of step (n).
c          "u" will be updated to n+1 in this routine.
c          at bottom here we set du = u(n+1) - u(n) to get
c          ready for extrapolation in next load(time) step.
c          Needed since adaptive will set du to be just a fraction
c          of total du over n -> n+1.
c
      allocate( u_n_local(nodof) )
      u_n_local(1:nodof) = u(1:nodof)
c
c          used to track how the step converged for updating
c          global onvergence_history info
c
      info%adaptive_used = .false.
      info%step_converged = .false.
      info%adapt_substeps = 0
      info%iters_at_convergence = 1
      info%adapt_iters_at_convergence = 0
c
c          forcing of linear stiffness at start of load step now
c          means to turn off extrapolation (if on) for just the new
c          step. the strain-stress updating routines for iter = 0
c          will return the linear [D] for integration points
c          if extrapolation = .false.
c
      user_extrapolation_on = extrapolate ! better name for local use
      user_extrapolation_off = .not. extrapolate ! convenience
      emit_extrap_msg = .true.
      emit_forced_linear_k_for_step = .true.
      hypre_solver = solver_flag .eq. 9
      local_mf_ratio_change = mf_ratio_change
c
c          initialize adaptive solution stack for the load step. we try
c          to solve the entire step in one "increment" first.
c
      if( local_debug ) write(iout,9000) step, mf, mf_nm1,
     &                  local_mf_ratio_change
c
      dynamic      = total_mass .gt. zero .and. dt .gt. zero
      adaptive     = adaptive_flag .and. step .gt. 1
      first_subinc = .false.  ! 1st substep for adaptive
      call adapt_check( scaling_adapt, 3, step, out )
      dt_original  = dt
c
c           calculate the scaling value for extrapolated displacements.
c           this value is used for real steps only (not subincrements).
c           the scaling value is modified if the step is subincremented.
c           mf_ratio_change = .true. means the loading routines believe
c           the next step is not proportional to the last step
c           (then turn off extrapolation for the step which forces [D]
c           linear to be used on iteration 1)
c
      if( abs(mf_nm1) .le. ext_tol ) mf_nm1 = one
      scaling_load = mf / mf_nm1
      if( local_mf_ratio_change ) scaling_load = zero
      if( local_debug ) write(out,9002) mf_nm1, scaling_load
c
      if( local_debug ) write (iout,9500) dynamic,adaptive,dt_original
c
c           adjust global load vector (load) for the gradual release of
c           any user-specified absolute constraints. The reduction
c           is independent of scaling, adaptivity.
c
c           These releases are generated by the user thru the
c           release constraints command and are not part of the
c           built in crack growth by node release facilities.
c
      call mnralg_release_con_forces
c
c           initialize data structures which handle the temperature
c           loadings in the case of adaptive step reduction.
c
c           for MPI::
c               send workers temperature information
c               and element equiv. nodal loads for this step.
c
      call wmpi_send_temp_eqloads
c
c           save user-defined temperature increments for nodes and
c           (uniform) element increments in local save vectors
c           in case an adaptive solution is executed for step.
c
      if( temperatures ) call mnralg_scale_temps( 1, iout )
c
c
c          initialize solution for this load step  (iter = 0 part)
c          -------------------------------------------------------
c
c          this is target for starting an adaptive subincrement
c          within the step. save critical values for solution state
c          needed to support adaptive restart of the step if req'd.
c
 1000 continue
c
c          incorporate changes in load vector due to node/element
c          releases for crack growth. skip for step 1. skip when
c          this is the first subincrement in adaptive for step
c          since checking/loads adjustment will have already been done
c          for step. no crack growth of any type allowed during
c          adaptive solution since we're trying to resolve only
c          the Newton.
c
      check_crk_growth = step .gt. 1  .and.  (.not. first_subinc)
      if( check_crk_growth ) then
         call chkcrack( step, 0 )
         call growth_loads
      end if
      if( show_details ) write(out,9150) step
c
c          set up the displacement increment to start the step or
c          adaptive subincrement. vector du contains converged
c          displacement increment from previous step
c
c          options for extrapolating (predct) converged displacement
c          increments from the last step. none, (du=0) or a
c          computed factor based on the load step factors s
c          specified by the user.
c
c          the scaling factor for the step is a function of the
c          scaling factor due to changes in the load mulitplier and
c          the scaling factor due to subincrementation.
c
c          extrapolation on is now the default with new solution
c          algorithms in Oct. 2015.
c
c          the linear stiffness on iteration 1 (extrapolation=off)
c          is enforced if
c
c            - user requested extrapolation is off
c            - the user requested it just for this step with the
c              nonlinear parameters command (but not
c               substeps if adative kicks in)
c            - the loading routines determined that the current step
c              is not proportional to last one (mf_ratio_change
c              = true.). Not for substeps if adaptive kicks in
c            - we just reduced the step size for adaptive. use linear
c              [D] for 1st substep after reduction
c
c          extrapolated_du keeps track of what
c          to do about du extrapolation for just this step and adaptive
c          substeps as needed.
c
c          no extrapolation is possible for step 1. further, step 1
c          cannot have adaptive solution in this version.
c
c          user specified extrapolation value (prdmlt) no longer
c          supported
c
c          note local variables user_extrapolation_on,
c          user_extrapolation_off for simpler code
c
      scaling_factor = scaling_load * scaling_adapt
c
      if( local_debug ) write (iout,9600)
     &      scaling_load, scaling_adapt, scaling_factor
c
      if( first_subinc ) then ! easiest case
         extrapolated_du = .false.
         du(1:nodof) = zero
         if( show_details ) write(out,9154) step
      elseif( extrap_off_next_step ) then
         extrapolated_du = .false.
         du(1:nodof) = zero
         if( show_details ) write(out,9154) step
         extrap_off_next_step = .false.
      else ! regular load step or continuation of adaptive substep
         extrapolated_du = .false.
         if( user_extrapolation_on .and. step .gt. 1  ) then
           extrapolated_du = .true.
           if( local_mf_ratio_change  ) then ! non-proportional loading
             extrapolated_du = .false.
             call errmsg ( 255, dum, dums, dumr, dumd)
             emit_extrap_msg = .false. ! no need for repeated messages
           else
             if( abs( scaling_factor ) > 5.0d0 .and. show_details ) then
               if( show_details ) write(out,9158) step, scaling_factor
             else
               if( show_details ) write(out,9152) step, scaling_factor
             end if
             if( scaling_factor < zero ) then
               scaling_factor = zero
               if( show_details ) write(out,9153) scaling_factor
             end if
             if( scaling_factor > 5.0d0 ) then
                scaling_factor = one
                if( show_details ) write(out,9153) scaling_factor
             end if
             extrapolated_du = .true.
             du(1:nodof) =  scaling_factor * du(1:nodof)
             if( local_debug ) write (iout,9530) scaling_factor
           end if
        else ! no extrapolation to start this step.
           du(1:nodof) = zero
           if( show_details ) write(out,9154) step
        end if
      end if
c
c          apply constraints to the starting displacement change for
c          step. non-zero constraints are scaled by current adaptive
c          factor (<= 1). the extrapolated displacements are thus
c          scaled back if the load step has been subdivided.
c          set global variable if user-imposed constrains have non-zero
c          values -- used during iter=0 as part of stress update
c          to get loading for step
c
      non_zero_imposed_du = .false.
      do i = 1, nodof
         if( cstmap(i) .ne. 0 ) then
            du(i) = cnstrn(i) * adapt_disp_fact
            if( du(i) .ne. zero ) non_zero_imposed_du = .true.
         end if
      end do
c
c          if adaptive for a (real) dynamic solution, then cut time
c          step and save original for restoration at end of step.
c
      if ( adaptive ) then
           dt = dt_original* adapt_disp_fact
           if ( local_debug ) write (iout,9510) dt_original,
     &                        adapt_disp_fact, dt
      end if
c
c          MPI solution::
c            send worker processes information needed before the first
c            Newton iteration. Workers inform root the stress-strain
c            values it has for elements it does not own that
c            are no longer valid.
c
c          run uexternaldb for Abaqus support
c
      call wmpi_send_step
      call wmpi_send_itern
      call wmpi_get_str ( 1 )
      douextdb = 3  ! in common.main. tells uexternal what to do
      call wmpi_do_uexternaldb
c
c          re-define the step temperature increment for nodes
c          and elements in case we are in adaptive solution
c
      if( temperatures ) call mnralg_scale_temps( 2, iout )
c
c
c          If contact has been enabled for this analysis, then search
c          for contact against all the nodes.  If contact is found,
c          compute the contact force.
c
      call contact_find
c
c          perform the strain-stress recovery and compute the internal
c          forces to reflect non-zero displacements applied before
c          the step (either extrapolate or non-zero constraints). this is
c          a simple way to get the "applied" forces equivalent to the
c          imposed, incremental displacements for step.
c
      if( show_details ) write(out,9155) step
      material_cut_step = .false.
      call drive_eps_sig_internal_forces( step, 0,
     &          material_cut_step )
c
c          element stiffness matrices. (stifup gets only new element
c          stiffnesss - not a structure stiff).
c          We're processing stiffness at start of load step (iter=0)
c
      call stifup( step, 1, out, newstf, show_details )
c
c          setup nodal forces that enforce for MPCs. we call them
c          Lagrange forces. See comments in module stiffness_data.
c
      if( mpcs_exist .or. tied_con_mpcs_constructed ) then
        if( .not. allocated( total_lagrange_forces ) ) then
           allocate( total_lagrange_forces(nodof) )
           total_lagrange_forces = zero
        end if
        if( .not. allocated( d_lagrange_forces ) )
     &            allocate( d_lagrange_forces(nodof) )
        d_lagrange_forces = zero
        if( .not. allocated( i_lagrange_forces ) )
     &            allocate( i_lagrange_forces(nodof) )
        i_lagrange_forces = zero
      end if
c
c          compute the effective load increment to start load step
c          including real applied (nodal) load increments, estimated
c          inertia and estimated load due to imposed displacement/temp
c          increments. result is stuffed into the "res" vector. that
c          vector determines incremental displacements for iteration
c          1 of step.
c
      call uppbar( pbar, mdiag, a, v, load, dload, nbeta, dt,
     &             adapt_load_fact, nodof, out )
      call upres_iter_0( mgres1, nbeta, dt, nodof, out, cstmap,
     &                   pbar, ifv, mdiag, du, res, dstmap )
c
c          if mgres1 (the square root of the sum of the squares of the
c          res) is zero then no load has been applied for this step,
c          including directly applied loads, inertial forces, contact
c          forces, imposed displacements, temperature increments.
c          This means that the effective load increment for the step
c          = {0}.  This  implies n solution is needed for step.
c
c          This needs modification to support creep loading at
c          fixed external loads. In this case, mgres1 = 0
c          and imposed displacement increments = 0. We should
c          let solution continue here so that "fake" v. small loads
c          are not needed in such creep solutions.
c
      if( mgres1.eq.zero .and. zrocon ) then
         call errmsg(230,step,dums,dumr,dumd)
         iter = 0
         go to 100
      endif
c
c          begin netwon iterations loop to reduce residual forces
c          ------------------------------------------------------
c
      iter              = 1
      material_cut_step = .false.
      first_solve       = .true.
c
c          target for top of iteration loop. update element [K]s.
c          we add nodal mass (with newmark beta and dt) to
c          diagonal terms of element stiffnesses.
c
 25   continue
      msg_count_1 = 0 ! used to prevent excessive messages in low level routines
      msg_count_2 = 0
      if( iter .gt. 1 ) ! start of step done above
     &         call stifup( step, iter, out, newstf, show_details )
c
c          if we're using hypre solver, we need to distribute the
c          stiffness to MPI ranks (that's the second flag).
c          Gets set below either way, but all
c          parallel assembly must be distributed
c          some features used in the model may disable
c          distributed assembly in MPI.
c
c          MPCs/tied contact  -> no hypre, symmetric only
c          hypre -> only with MPI. ok without parallel
c                   assembly. equations assembled on root
c                   then distributed to ranks.
c
      parallel_assembly_used = .false.
      if( .not. (mpcs_exist .or. tied_con_mpcs_constructed)
     &     .and. hypre_solver .and. parallel_assembly_allowed ) then
            parallel_assembly_used = .true.
            distributed_stiffness_used = .true.
      end if
c
c          we can use hypre (MPI) w/o distributed assembly.
c
      distributed_stiffness_used = .false.
      if( hypre_solver ) distributed_stiffness_used = .true.
c
      if( dynamic ) call inclmass
c
c           Do some checks for things we haven't implemented yet
c
      if( asymmetric_assembly .and. parallel_assembly_used ) then
            write(out,9720); write(out,9718)
            call die_gracefully
      end if
c
      if( asymmetric_assembly .and. (mpcs_exist .or.
     &      tied_con_mpcs_constructed) ) then
            write(out,9715); write(out,9718)
            call die_gracefully
      end if
c
c          MPI:
c            default is to gather all of the element stiffnesses back
c            to the root processor so that we can conduct assembly of
c            model stiffness matrix.
c
c            with hypre solver and under right conditions (see above)
c            we can use the faster, less memory or root distributed
c            assembly.
c
      if( .not. parallel_assembly_used ) then
            call t_start_assembly( start_assembly_step )
            call wmpi_combine_stf
            call t_end_assembly( assembly_total, start_assembly_step )
      end if
c
c          solve the linearized equations to compute the displacement
c          increment idu and update the total displacement increment
c          for the load step (du).
c
c          Enforce MPCs when:
c              - symmetric MKL solvers (direct/iterative)
c              - no MPI
c
      call eqn_solve( iter, step, first_solve,
     &                nodof, solver_flag, use_mpi, show_details,
     &                du, idu, iout )
c
c          This applies only to hypre:  If we return from the solver
c          and hypre has indicated that it needs an adaptive step
c          (recall issue with solution of iterations immediately
c          preceding an adaptive call) then get the next scaling
c          factor, call for a reset and return to the start
c          of the step.
c
      if( hyp_trigger_step ) then
               call adapt_check( scaling_adapt, 1, step, out )
               call adaptive_reset
               first_subinc = .true.
               go to 1000
      end if
c
c          run strain-stress-internal force update. include line
c          search in this process if option in on.
c          an "instrumented" line search is available for research
c          work on algortihms
c
c      call mnralg_ls_instrumented
      call mnralg_ls
c
c          did a faiure occur in strain computation for finite
c          strains (e.g. inverted element) or a material stress
c          update process faile ? then use adaptive immediately.
c
      if(  material_cut_step ) then ! from a material model across all ranks
       if( adaptive  ) then
         if( show_details ) write(out,9420)
         go to 50
       else
         call abort_job
       end if
      end if
c
      if( ls_request_adaptive ) then ! ls search failed
       if ( adaptive  ) then
         if ( show_details ) write(out,9425)
         go to 50
       else
         call abort_job
       end if
      end if

c
      if( prnres ) call oures( ldnum, iter )
c
c          perform convergence tests for newton iterations. even
c          if apparently converged, the user may want a set minimum
c          number of iterations.
c
c          also check for possible divergence of the iterations
c          (if global flag is on). can trigger immediate adaptive or
c          termination.
c
      call cvtest( cnverg, step, iter, magdu1, mgload, adapt_load_fact,
     &             diverging_flag  )
      if( cnverg .and. iter .ge. mniter ) go to 100
      iter = iter + 1
      info%iters_at_convergence = info%iters_at_convergence + 1
      if( iter .le. mxiter ) then
        if( .not. divergence_check ) go to 25 ! next iteration
        if( .not. diverging_flag ) go to 25   ! next iteration
      end if
c
c          iteration number beyond prescribed limit. computations maybe
c          halted or the step subdivided using the adaptive strategy.
c          we automatically save for restart if the solution does not
c          converge.
c
      if( show_details ) write(out,9410) step
 50   continue
      if( adaptive ) then
          info%adaptive_used = .true.
          call adapt_check( scaling_adapt, 1, step, out )
          if( adapt_result .eq. 1 ) go to 100
          if( adapt_result .eq. 2 ) then
               call adaptive_reset
               first_subinc = .true.
               go to 1000  ! do over with smaller load inrement
          end if
          if( adapt_result .eq. 3 ) call abort_job
      end if
      iter = iter -1
      info%iters_at_convergence = info%iters_at_convergence -1
      if( halt ) call abort_job
c
c          iterations for the current step produced a converged
c          solution or the user wants to keep going if non-converged.
c          update global solution, n -> n+1. save critical values
c          for solution state needed to support adaptive
c          restart of the step.
c
 100  continue
      call update
      call adaptive_save
c
c          reset time step to original size
c
      total_model_time = total_model_time + dt
c
c          uexternaldb for Abaqus compatible support
c
      douextdb = 4  ! common.main. tells uexternaldb what to do
      call wmpi_do_uexternaldb
c
      if( adaptive ) dt = dt_original
c
c          MPI:
c            if we are using crack growth, then the slave processors
c            must send to the root processor the pertinent element
c            information for killable elements so that root can
c            assess the current growth criterion.
c
      call wmpi_get_grow
c
c          output convergence status.
c
      call mnralg_ouconv( 1, iter, step, cnverg, out,
     &                    total_model_time, adapt_load_fact  )
      adapt_result = 1
      scaling_load = one
      if( adaptive ) then
          call adapt_check( scaling_adapt, 2, step, out )
          info%adapt_iters_at_convergence =
     &        info%adapt_iters_at_convergence + iter
          info%adapt_substeps =
     &        info%adapt_substeps + 1
      end if
      if( cnverg ) call energy ( step, adapt_result, mdiag )
      first_subinc = .false.
      local_mf_ratio_change = .false.
      if( adapt_result .eq. 2 ) go to 1000
c
      if( cnverg ) info%step_converged = .true.
      if( local_debug ) write(iout,9800) info%adaptive_used,
     &     info%adapt_substeps, info%iters_at_convergence,
     &     info%adapt_iters_at_convergence, info%step_converged
c
c          update globally stored convergence history with info
c          for this step
c
      call update_convergence_history( step, iout, local_debug,
     &                                 info )
c
c          output the summary of loading pattern multipliers
c          imposed on model thru the step that just converged.
c
      call eqiv_out_patterns
c
c          update the last time step computed and the last nonlinear
c          load name successfully used.
c
      ltmstp = ltmstp + 1
      lsldnm = lodnam(ldnum)
c
c          release space used by additional temperature vectors needed
c          to support adaptive processing
c
      if( temperatures ) call mnralg_scale_temps( 3, iout )
c
c          print element killing information if requested.
c
      call dam_print( step, iter )
c
c          get total displacement change over n -> n+1. du
c          will not be this if adaptive was used over n -> n+1
c
      du(1:nodof) = u(1:nodof) - u_n_local(1:dof)
      deallocate( u_n_local )
c
      return  ! back to stpdrv
c
 9000 format(/,1x,'... mnralg: step, mf, mf_nm1,'
     & ' mf_ratio_change: ',i7,2f10.6, l2)
 9002 format(/,1x,'... mnralg:  mf_nm1, scaling_load',
     & 2f10.6)
 9100 format(7x,
     & '>> updating strains/stress/forces step, iteration: ',i7,i3)
 9110 format(7x,
     & '>> updating internal forces step, iteration:       ',i7,i3)
 9150 format(/,7x,
     & '>> starting global newton solution for step:       ',i7)
 9152 format(7x,
     & '>> extrapolating displacements to start step:      ',i7,
     & ' (* ',f7.3,')')
 9153 format(7x,
     & '>> extrapolation scale factor reset to :           ',f7.3)
 9154 format(7x,
     & '>> no extrapolating displacements to start step:   ',i7)
 9155 format(7x,
     & '>> delta forces for loads/displ/temps/creep. step: ',i7)
 9158 format(7x,
     & '>> extrapolating displacements to start step:      ',i7,
     & ' (* ',e10.3,')')
 9410 format(/1x,'>> iteration limit exceeded for current step:',i7,
     & /,1x       '   (or adaptive sub-increment) or the solution ',
     &   'appears to be diverging....')
 9420 format(/,2x,'>>> material model stress update computations',
     &     /,2x,'    requested load step reduction...')
 9425 format(/,2x,
     &'>>> line search has requested load step reduction...')
 9500 format (/,'... mnralg:',
     & /8x,'Dynamic parameter = ',l2,
     & /8x,'Adaptive solution = ',l2,
     & /8x,'Original time Step (before subincrementing) = ',
     &        e16.6 )
 9510 format(/,'... mnralg:',
     & /,8x,'Original time step = ',e16.6,
     & /,8x,'Adaptive factor = ', e16.6,
     & /,8x,'Current time step = ',e16.6)
 9530 format(/,'... mnralg:',
     & /,8x,'>>>> Scaling previous displacements by   ',e16.6)
 9600 format(/,'... mnralg:',
     & /,8x,'scaling_load, scaling_adapt, scaling_factor ',
     &             3e16.6,// )
 9610 format(3x,i8,2x,e14.6)
 9700 format(10x,i4, 2x, e14.6, 2x, e14.6)
 9715 format(/,'>> ERROR: tied contact and/or multi-point constraints',
     &    /, '          cannot be used with the asymmetric assembly-',
     &    /, '          solver process' )
 9718 format(/,'>> Terminating analysis at this point',//)
 9720 format(/,'>> ERROR: asymmetric stiffness not implemented with'
     &    /, '          distributed (parallel) assembly' )
 9800 format(2x,'... leaving mnralg info data... ',
     & /,10x,'adaptive_used:              ',l2,
     & /,10x,'adapt_substeps:             ',i2,
     & /,10x,'iters_at_convergence:       ',i2,
     & /,10x,'adapt_iters_at_convergence: ',i2,
     & /,10x,'step_converged:             ',l2, //)
c
      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *                     subroutine mnralg_ls                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/3/2015 rhd              *
c     *                                                              *
c     *    run line search for this global Newton iteration if       *
c     *    user requested . instrumented version to gather behavior  *
c     *    for all LS iterations to support algorithm tuning         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine mnralg_ls
      implicit none
c
c              local variables
c
      double precision ::
     &      s_value, s0_value, r_value, alpha, rho, slack_toler,
     &      alpha_min, s_values(0:30), r_values(0:30), zero, one,
     &      ls_reduce_fraction, alpha_values(30)
      double precision, allocatable :: du0(:)

      logical :: ls_debug, no_line_search, line_search_details
      integer :: ls_ell, i
      character(len=1) :: ls_flag
      data zero, one / 0.0d00, 1.0d00 /
c
      if( show_details ) write(iout,9210) step, iter
c
      ls_debug = .false.
      line_search_details = ls_details
      if( ls_debug ) write(iout,8900) step, iter
      rho = ls_rho
      no_line_search = .not. line_search
      slack_toler = ls_slack_tol
      alpha_min = ls_min_step_length
      ls_request_adaptive = .false.
      ls_reduce_fraction = 0.9d00 ! could make user input
      ls_flag = "*"
c
c              1. we always compute s0 even if the user does not
c                 request line search. if s0 > 0 we issue a
c                 warning message that the global NR is in
c                 in an uphill situation.
c
c                 if s0 >=0 start next Newton iteration
c                 even if user requested line search.
c
c                 may want to consider immediate adaptive.
c
c                 skip MPC dof if present in s computations.
c                 use separate loops for efficiency.
c
      call mnralg_ls_get_s( s0_value )
      if( s0_value > zero .and. show_details )
     &       write(iout,9010) s0_value, step, iter
      if( s0_value >= zero .or. no_line_search ) then
        call mnralg_ls_update_res_0
        return
      end if
      if( ls_debug ) write(iout,9000) s0_value
      if( line_search_details ) write(iout,8930)
c
c              2. loop to find acceptable search length (alpha) if
c                 one exists
c
      alpha = ls_max_step_length  ! defaults to 1.0
      ls_ell = 1
      allocate( du0(1:nodof) )
      du0(1:nodof) = du(1:nodof)  ! since residual/contact/... use du
      s_values = zero
      r_values = zero
      s_values(0) = s0_value
c
      do ! step length search
c
c              2a. new alpha (step length). get s-value.
c                   material routines may want immediate adaptive.
c
         if( ls_debug ) write(iout,8910) ls_ell, alpha
         call mnralg_ls_update_res( alpha, du0 )
         if( material_cut_step ) then
           deallocate( du0 )
           return
         end if
         call mnralg_ls_get_s( s_value )
         s_values(ls_ell) = s_value
c
c              2b. compute r-ratio. check for convergence of step
c                  length
c
         r_value = abs( s_values(ls_ell) / s_values(0) )
         r_values(ls_ell) = r_value
         alpha_values(ls_ell) = alpha
         if( ls_debug ) write(iout,9030) s_value, r_value
c
         if( r_value < slack_toler ) then ! converged
             if( ls_debug ) write(iout,8920) slack_toler
             if( line_search_details ) write(iout,9205)  ls_ell,
     &                                    alpha, r_value
             if( show_details ) write(iout,9200) alpha, ls_flag,
     &                                           step, iter
             deallocate( du0 )
             return
         end if
         if( line_search_details ) write(iout,9205)  ls_ell,
     &                 alpha_values(ls_ell), r_values(ls_ell)
c
c              2c. if r increases or if r is not decreasing
c                  fast enough, then stop ls. use step length
c                  = 1.0
c
         if( ls_ell .ge. 2 ) then ! increasing r ?
           if( r_values(ls_ell) .gt. ls_reduce_fraction *
     &         r_values(ls_ell-1) ) then ! r-value not decreasing
               alpha = one
               ls_request_adaptive = .false.
               call mnralg_ls_update_res( alpha, du0 )
               deallocate( du0 )
               ls_flag = "#"
               if( show_details ) write(iout,9200) alpha, ls_flag,
     &                            step, iter
               return
           end if
         end if
c
c              2d. get new step length. simple back track for now.
c                  check if smaller than min step length.
c                  up to top of search loop.
c
         alpha = rho * alpha  ! simple back track.
         if( alpha .lt. alpha_min ) then
             ls_flag = "%"
             if( show_details ) write(iout,9100) ls_ell,
     &            step, iter, r_value, alpha, ls_flag
             ls_request_adaptive = .false.
             alpha = one
             call mnralg_ls_update_res( alpha, du0 )
             deallocate( du0 )
             return
         end if
         ls_ell = ls_ell + 1
      end do
c
      return
c
 8900 format(6x,'.... entered line search routine. step, iter: ',
     &  i7,i3)
 8910 format(10x,'.... start new ls iteration, alpha: ',i2,f10.4)
 8920 format(20x,'satisfied search tolerance of: ',f10.3)
 8930 format(7x,
     & '>> begin line search ' )
 9000 format(3x,'       ... s0 value: ', e14.6)
 9010 format(7x,
     &      '>> line search s0 value:',f8.3,' > 0.',14x,i7,i3,
     & /,7x,'   using step length 1.0 this iteration')
 9030 format(20x,'new s, r:  ',2f10.4)
 9100 format(7x,
     &    '>> line search did not converge. iterations: ',i3,3x,i7,i3,
     &/7x,'      r-ratio: ',f8.2,' step length:', f6.3,a2)
 9200 format(7x,
     & '>> line search completed. step length: ',f6.3,a2,4x,i7,i3)
 9205 format(7x,
     & '     line search iteration, alpha, r-value: ',i2,f6.3,f7.3)
 9210 format(7x,
     & '>> updating strains/stress/forces step, iteration: ',i7,i3)
 9220 format(7x,
     & '      ... finished line search. iterations, step length: ',
     &   i2,1x,f7.4)
c
      end subroutine mnralg_ls
c     ****************************************************************
c     *                                                              *
c     *              subroutine mnralg_ls_instrumented               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/3/2015 rhd              *
c     *                                                              *
c     *    run line search for this global Newton iteration if       *
c     *    user requested . instrumented version to gather behavior  *
c     *    for all LS iterations to support algorithm tuning         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine mnralg_ls_instrumented
      implicit none
c
c              local variables
c
      double precision ::
     &      s_value, s0_value, r_value, alpha, rho, slack_toler,
     &      alpha_min, s_values(0:30), r_values(0:30), zero, one,
     &      ls_reduce_fraction, alpha_values(30)
      double precision, allocatable :: du0(:)

      logical :: ls_debug, no_line_search, line_search_details
      integer :: ls_ell, i, curve_type, bug
      character(len=1) :: ls_flag
      data zero, one / 0.0d00, 1.0d00 /
c
      if( show_details ) write(iout,9210) step, iter
c
      ls_debug = .false.
      line_search_details = .true.
      if( ls_debug ) write(iout,8900) step, iter
      rho = ls_rho
      no_line_search = .not. line_search
      slack_toler = ls_slack_tol
      alpha_min   = 0.01 ! ls_min_step_length
      ls_request_adaptive = .false.
      ls_reduce_fraction = 0.9d00
      ls_flag = "*"
c
c              we always compute s0 even if the user does not
c              request line search. if s0 > 0 we issue a
c              warning message that the global NR is in in an uphill
c              situation.
c
c              if s0 >=0 start next Newton iteration
c              even if user requested line search.
c
c              may want to consider immediate adaptive.
c
c              skip MPC dof if present in s computations.
c              use separate loops for efficiency.
c
      call mnralg_ls_get_s( s0_value )
      if( s0_value > zero .and. show_details )
     &       write(iout,9010) s0_value, step, iter
      if( s0_value >= zero .or. no_line_search ) then
        call mnralg_ls_update_res_0
        return
      end if
      if( line_search_details ) write(iout,8930)
c
c              loop over all step lengths. compute r-value
c              and print for studies of how r-values change
c              with step lengths
c
      alpha = 1.0d00
      ls_ell = 1
      allocate( du0(1:nodof) )
      du0(1:nodof) = du(1:nodof)  ! since residual/contact/... use du
      s_values = zero
      r_values = zero
      s_values(0) = s0_value
c
      do ! step length search
c
         call mnralg_ls_update_res( alpha, du0 )
         if( material_cut_step ) then
           deallocate( du0 )
           return
         end if
         call mnralg_ls_get_s( s_value )
         s_values(ls_ell) = s_value
c
c              2. compute r-ratio. check for convergence of step
c                 length
c
         r_value = abs( s_values(ls_ell) / s_values(0) )
         r_values(ls_ell) = r_value
         alpha_values(ls_ell) = alpha
         if( line_search_details ) write(iout,9205)  ls_ell,
     &                 alpha_values(ls_ell), r_values(ls_ell)
         alpha = rho * alpha  ! simple back track.
         if( alpha .lt. alpha_min ) then
             ls_flag = "%"
             ls_request_adaptive = .false.
             bug = 0
             if ( step .eq. 5 ) bug = 1
             call mnralg_ls_instrumented_a( ls_ell,
     &          r_values(1), alpha_values(1), curve_type, iout, bug )
              write(iout,9300) curve_type
             deallocate( du0 )
             return
         end if
         ls_ell = ls_ell + 1
      end do
      return

 9300 format(7x,
     & '      ... LS r-values curve type: ',i3)
 8900 format(6x,'.... entered line search routine. step, iter: ',
     &  i7,i3)
 8910 format(10x,'.... start new ls iteration, alpha: ',i2,f10.4)
 8920 format(20x,'satisfied search tolerance of: ',f10.3)
 8930 format(//,7x,
     & '>> begin line search ' )
 9000 format(3x,'       ... s0 value: ', e14.6)
 9010 format(7x,
     &      '>> line search s0 value:',f8.3,' > 0.',14x,i7,i3,
     & /,7x,'   using step length 1.0 this iteration')
 9030 format(20x,'new s, r:  ',2f10.4)
 9100 format(7x,
     &    '>> line search did not converge. iterations: ',i3,3x,i7,i3,
     &/7x,'      r-ratio: ',f8.2,' step length:', f6.3,a2)
 9200 format(7x,
     & '>> line search completed. step length: ',f6.3,a2,4x,i7,i3)
 9205 format(7x,
     & '     line search iteration, alpha, r-value: ',i2,f6.3,f7.3)
 9210 format(7x,
     & '>> updating strains/stress/forces step, iteration: ',i7,i3)
 9220 format(7x,
     & '      ... finished line search. iterations, step length: ',
     &   i2,1x,f7.4)
c
      end subroutine mnralg_ls_instrumented
c
      subroutine mnralg_ls_instrumented_a( npts, yvalues, xvalues,
     &                      curve_type, iout, bug )
      implicit none
      integer npts,curve_type, i, first_min, iout, bug
      double precision yvalues(*), xvalues(*), alpha,
     &                 ymin
      if( bug .eq. 1 ) then
      write(iout,*) '<<< num points: ', npts
      do i = 1, npts
         write(iout,*) i, xvalues(i), yvalues(i)
      end do
      end if

      curve_type = 1
      do i = 2, npts
        if( yvalues(i) < yvalues(i-1) ) cycle
          go to 200
      end do
      return

 200  curve_type = 2
      do i = 2, npts
        if( yvalues(i) > yvalues(i-1) ) cycle
          go to 300
      end do
      return

 300  curve_type = 3
      return
      end subroutine mnralg_ls_instrumented_a

c     ****************************************************************
c     *                                                              *
c     *  subroutines: mnralg_ls_update_res_0, mnralg_ls_update_res   *
c     *                mnralg_ls_get_s                               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/9/2016 rhd               *
c     *                                                              *
c     *    support routines to make ls code clean. these will be     *
c     *    inlined by compiler                                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mnralg_ls_update_res_0
      implicit none
c
      du(1:nodof) = du(1:nodof) + idu(1:nodof)
c
c          MPI:
c            send all the worker processors information determined from
c            the equation solve to help them continue with the newton
c            iterations. just du at present. drive ... supervises
c            ranks to compute strains, stresses, new internal forces.
c
      call wmpi_send_itern
c
      call drive_eps_sig_internal_forces( step, iter,
     &                                    material_cut_step )
      if( material_cut_step ) return ! reduced from all ranks
      call contact_find
      call upres( iter, out, nodof, dt, nbeta, num_term_loads,
     &            sum_loads, mgload, mdiag, pbar, du, v, a,
     &            ifv, res, cstmap, dstmap, load )
      return
      end subroutine mnralg_ls_update_res_0

      subroutine mnralg_ls_update_res( alpha, du0 )
      implicit none
      double precision :: alpha, du0(*)
c
      du(1:nodof) = du0(1:nodof) + alpha*idu(1:nodof)
c
c          MPI:
c            send all the worker processors information determined from
c            the equation solve to help them continue with the newton
c            iterations. just du at present. drive ... supervises
c            ranks to compute strains, stresses, new internal forces.
c            iterations.
c
      call wmpi_send_itern
c
      call drive_eps_sig_internal_forces( step, iter,
     &                                    material_cut_step )
      if( material_cut_step ) return ! reduced from all ranks
      call contact_find
      call upres( iter, out, nodof, dt, nbeta, num_term_loads,
     &            sum_loads, mgload, mdiag, pbar, du, v, a,
     &            ifv, res, cstmap, dstmap, load )
      return
      end subroutine mnralg_ls_update_res

      subroutine mnralg_ls_get_s( svalue )
      implicit none
      double precision :: svalue
      integer :: i
c
      svalue = zero
      if( allocated( d_lagrange_forces ) ) then
        do i = 1, nodof
         if( cstmap(i) .ne. 0 ) cycle ! abs constraint
         if( d_lagrange_forces(1) .ne. zero ) cycle
         svalue = svalue  - res(i)*idu(i)
        end do
      else  ! no mpcs active
        do i = 1, nodof   !   fix for MPCs
         if( cstmap(i) .ne. 0 ) cycle ! abs constraint
         svalue = svalue  - res(i)*idu(i)
        end do
      end if
      return
      end subroutine mnralg_ls_get_s
c
      end subroutine mnralg





c
c     ****************************************************************
c     *                                                              *
c     *                subroutine update_convergence_history         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 08/15/2013                 *
c     *                                                              *
c     * load step completed. update global convergence history       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine update_convergence_history( step, iout, local_debug,
     &                                       info )
c
      use main_data, only : convergence_history
      implicit integer (a-z)
c
      logical local_debug
c
      type :: info_mnralg
        logical :: adaptive_used
        logical :: step_converged
        integer :: adapt_substeps
        integer :: iters_at_convergence
        integer :: adapt_iters_at_convergence
      end type
      type( info_mnralg ) info
c
c          we save convergence history for 5 most recent steps.
c          the step just completed is in position 5. bump prior
c          steps up in vector to make room.
c
      do i = 2, 5
       j = i - 1
       convergence_history(j)%step_converged =
     &       convergence_history(i)%step_converged
       convergence_history(j)%adaptive_used =
     &       convergence_history(i)%adaptive_used
       convergence_history(j)%iterations_for_convergence =
     &       convergence_history(i)%iterations_for_convergence
       convergence_history(j)%adapt_substeps =
     &       convergence_history(i)%adapt_substeps
      end do
c
      convergence_history(5)%step_converged = info%step_converged
      convergence_history(5)%adaptive_used = info%adaptive_used
      if( info%adaptive_used ) then
          convergence_history(5)%iterations_for_convergence =
     &          info%adapt_iters_at_convergence
          convergence_history(5)%adapt_substeps = info%adapt_substeps
      else
          convergence_history(5)%iterations_for_convergence =
     &          info%iters_at_convergence
          convergence_history(5)%adapt_substeps = 1
      end if
c
      if( local_debug ) then
        write(iout,9000)
        do i = 1, 5
          write(iout,9010) i, convergence_history(1)%step_converged,
     &       convergence_history(i)%adaptive_used,
     &       convergence_history(i)%iterations_for_convergence,
     &       convergence_history(i)%adapt_substeps
        end do
      end if
c
      return
 9000 format(/,1x,'... update_convergence_history ...')
 9010 format(/,5x,'(i): ',i2,
     & /,10x,'step_converged:             ', l2,
     & /,10x,'adaptive_used:              ', l2,
     & /,10x,'iterations_for_convergence: ', i2,
     & /,10x,'adapt_substeps:             ', i2,// )
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine adaptive_save                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 05/27/96                   *
c     *                                                              *
c     * save key solution variables in case the analysis requires an *
c     * adaptive restart with a smaller load step size.              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine adaptive_save
      use global_data ! old common.main
c
      use main_data, only       : du_nm1
      use elem_block_data, only : rot_n_blocks, rot_n1_blocks,
     &                            rot_blk_list
c
      implicit integer (a-z)
      double precision
     &    dummy(1)
c
c                      MPI:
c                        workers run this subroutine for the
c                        data they own.
c
       call wmpi_alert_slaves ( 18 )
c
c
c                      we also save the current 3x3 rotation matrices
c                      from the most recent polar decompositions as
c                      the values for start of new step.
c
c                      Note that for each block, we handle the data
c                      only if the list entry for the block is 1, for
c                      instance rot_blk_list(blk)=1.  This provides a simple
c                      mechanism to skip unallocated blocks due either to
c                      processor assignment or material model.
c
c
      if ( .not. allocated(du_nm1) ) then
        allocate( du_nm1(nodof), stat=alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 1
           call die_abort
           stop
        end if
      end if
      call vec_ops( du_nm1, du, dummy, nodof, 5 )
c
      if ( .not. allocated(rot_n1_blocks) ) go to 200
      do blk = 1, nelblk
        if ( rot_blk_list(blk) .eq. 1 ) then
          felem       = elblks(1,blk)
          span       = elblks(0,blk)
          ngp        = iprops(6,felem)
          block_size = span * ngp * 9
          call vec_ops( rot_n_blocks(blk)%ptr(1),
     &                  rot_n1_blocks(blk)%ptr(1), dummy,
     &                  block_size, 5 )
        end if
      end do
c
 200  continue
      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 adaptive_saves: @',i2,/,
     &       '>>> Job terminated....')
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine adaptive_reset               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified :  3/4/13 rhd                *
c     *                                                              *
c     * the adaptive processor has decided to re-start solution for  *
c     * the load step. reset key data arrays to those at the         *
c     * beginning of the step. this is based on analysis type        *
c     * and the type of material model                               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine adaptive_reset
      use global_data ! old common.main
c
      use main_data,       only : du_nm1, nonlocal_analysis
      use elem_block_data, only : history_blocks, history1_blocks,
     &                            rot_n_blocks, rot_n1_blocks,
     &                            eps_n_blocks, eps_n1_blocks,
     &                            urcs_n_blocks, urcs_n1_blocks,
     &                            history_blk_list, rot_blk_list,
     &                            eps_blk_list,
     &                            urcs_blk_list, nonlocal_flags,
     &                            nonlocal_data_n, nonlocal_data_n1
      implicit integer (a-z)
c
c                      lcoal values
c
      double precision
     &    dummy(1)
      logical chk
c
c                      MPI:
c                        tell the slaves to run this subroutine for the
c                        data they own.
c
      call wmpi_alert_slaves ( 19 )
c
c                      the displacement increment for the previous
c                      step or converged sub-increment due to
c                      adaptive solution.
c
      call vec_ops( du, du_nm1, dummy, nodof, 5 )
c
c                      reset:
c                       1) the structural history data,
c                       2) 3x3 rotations for material points (geo.
c                          nonlin elements) and
c                       3) strains
c                       4) unrotated cauchy stresses
c
c                      Note that for each block, we do the copy of data
c                      only if the list entry for the block is .ne. 0, for
c                      instance history_blk_list(blk)>0. This provides a
c                      simple mechanism to skip unallocated blocks due
c                      either to processor assignment or material model.
c
      do blk = 1, nelblk
        felem      = elblks(1,blk)
        span       = elblks(0,blk)
        ngp        = iprops(6,felem)
        hist_size  = span * ngp * history_blk_list(blk)
        rot_size   = span * ngp * 9
        eps_size   = span * ngp * nstr
        sig_size   = span * ngp * nstrs
        if ( hist_size .gt. 0 )
     &      call vec_ops( history1_blocks(blk)%ptr(1),
     &                    history_blocks(blk)%ptr(1), dummy,
     &                    hist_size, 5 )
        if ( rot_blk_list(blk) .eq. 1 )
     &      call vec_ops( rot_n1_blocks(blk)%ptr(1),
     &                    rot_n_blocks(blk)%ptr(1), dummy,
     &                    rot_size, 5 )
        if ( eps_blk_list(blk) .eq. 1 )
     &      call vec_ops( eps_n1_blocks(blk)%ptr(1),
     &                    eps_n_blocks(blk)%ptr(1), dummy,
     &                    eps_size, 5 )
        if ( urcs_blk_list(blk) .eq. 1 )
     &      call vec_ops( urcs_n1_blocks(blk)%ptr(1),
     &                    urcs_n_blocks(blk)%ptr(1), dummy,
     &                    sig_size, 5 )
      end do
c
c                      nonlocal shared state values if they exist.
c                      reset values at n+1 to those at n.
c
      if( .not. nonlocal_analysis ) return
      n = nonlocal_shared_state_size
      do i = 1, noelem
         if( .not. nonlocal_flags(i) ) cycle
         chk = allocated( nonlocal_data_n(i)%state_values ) .and.
     &         allocated( nonlocal_data_n1(i)%state_values )
        if( chk ) then
           nonlocal_data_n1(i)%state_values(1:n) =
     &              nonlocal_data_n(i)%state_values(1:n)
        else
           write(out,9100) i
           call die_abort
        end if
      end do
c
      return
c
 9100 format(">>>>> FATAL ERROR. adaptive_reset. nonlocal. elem: ",i8,
     &      /,"      Job terminated." )
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine abort_job                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 04/12/95                   *
c     *                                                              *
c     *     write a non-converged restart file and terminate         *
c     *     execution                                                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine abort_job
      use global_data ! old common.main
      implicit integer (a-z)
c
      call iodevn( idummy, iout, dummy, 1 )
      write(out,9000)
      call die_abort
      stop
c
 9000 format(//,' >>>>> analysis terminated.')
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine eqn_solve                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : rhd 11/21/2015 rhd         *
c     *                                                              *
c     *     solution of the linearized equilibrium equations         *
c     *     using an available equation solver.                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine eqn_solve( iter, step, first_solve,
     &                      nodof, solver_flag, use_mpi, show_details,
     &                      du, idu, iout )
c
      use main_data, only : extrapolated_du
      use hypre_parameters, only: hyp_trigger_step
      use stiffness_data, only: d_lagrange_forces,
     &                          i_lagrange_forces
c
      implicit none
c
c          parameters
c
      integer :: iter, step, solver_flag, iout, nodof
      logical :: first_solve, use_mpi, show_details
      double precision :: du(nodof), idu(nodof)
c
c          locals
c
      logical :: hypre, not_hypre, no_extrapolation, new_pre_cond,
     &           cpardiso, pardiso
c
c          solve for the change in the displacements for this iteration.
c          use either the Pardiso threaded solver (direct or iterative)
c          or the MPI_threads hypre solver from LLNL.
c
c          the EBE solver is deprecated  and removed.
c
c          the direct solver re-triangulates each time a solution
c          is obtained. it also maintains data structures for assembly
c          of [k-eff] structure dependent on constraints. for
c          predict on (iter=1) we skip the solution and go
c          straight to stress updating with predicted displacment
c          vector. set idu equal to the extrapolated displacements.
c
c          the first_solve flag tells the direct solver if it has done
c          a solution already for this step independent of what iter
c          equals. that way, it can change the solver data structure for
c          modified constraints even if we skipped the solution on iter=1
c          due to predict flag. first_solve = .true. if we have not actually
c          performed a triangulation yet for the load step. first_solve
c          is just local to this routine.
c
      hypre      = solver_flag .eq. 9
      pardiso    = solver_flag .eq. 7  .or.  solver_flag .eq. 8
      cpardiso   = solver_flag .eq. 10  .or.  solver_flag .eq. 11
      not_hypre = .not. hypre
c
c          MPI:
c            1) hypre -
c                  just solve using all mpi ranks + threads
c                  assembly on ranks or root
c            2) pardiso (threads only) -
c                   we still allow MPI + usual
c                   threaded pardiso. Workers
c                   perform element level computations. We bring
c                   [Ke]s back to root for assembly & solve.
c                   MPI workers ruinning on same node with root will steal
c                   spin-wait cycles from the solver use
c            3) cpardiso -
c                  just solve using all mpi ranks + threads.
c                  element [ks]s on ranks. assembly on root.
c                  cpardiso runs on root and workers and handles
c                  distribution of equations to ranks.
c
c          drive the equation solver. Code to allow changing
c          solver during a run has been removed.
c
c          give some help for an iterative solver about
c          when is a good point to build a new preconditioner
c
c          if using linear stiffness at beginning of a step
c          (no extrapolation of du),
c          the old conditioner is likely from a tangent stiffness
c          and not very good. iteration 2 of such a step should
c          then also have a new pre-cond for the tangent K
c          rather than linear K used in iteration one.
c
c          for now iter = 1 with extrapolation gets hint of no
c          new pre-conditioner.
c          may want to make this a user option. when the user knows
c          that large stiffness changes can occur from step n -> n+1
c          forcing a pre-conditioner update may be beneficial.
c
c          the solvers may choose to ignore the hint of new_pre_cond
c
      if ( not_hypre .and. show_details ) write(iout,9142) step, iter
      no_extrapolation = .not. extrapolated_du
      new_pre_cond = .false.
      if( iter .gt. 2 ) then
         new_pre_cond = .false.
      elseif( iter .eq. 1 )then
         new_pre_cond = .true.  ! linear K being used
      elseif( iter .eq. 2 .and. no_extrapolation ) then
         new_pre_cond = .true.  ! first tangent K after linear K
      end if
c
      call drive_assemble_solve( first_solve, iter, new_pre_cond )
      first_solve = .false.
c
c          check if hypre wants an adaptive step.
c          We can only deal with this in the main mnralg, so
c          simply keep going up the stack.
c
      if( hyp_trigger_step ) return
c
c          update after this iteration:
c
c           - du is estimated change in displacement from n -> n+1
c             now updated in line search routine here
c           - estimated change in Lagrange forces from n -> n+1
c                for MPCs (always self-equilibrating)
c
c
      if( allocated( d_lagrange_forces ) ) then
         if( .not. allocated( i_lagrange_forces ) ) then
            write(iout,9200); call die_gracefully
         end if
         d_lagrange_forces(1:nodof) = d_lagrange_forces(1:nodof) +
     &                                i_lagrange_forces(1:nodof)
      end if
c
      return
c
 9142 format(7x,
     & '>> running sparse solver step, iteration:          ',i7,i3)
 9200 format(1x,'>> FATAL ERROR: Job Aborted.',
     & /,5x,'i_lagrange_forces not allocated. eqn_solve')
c
      end

c     ****************************************************************
c     *                                                              *
c     *               subroutine mnralg_scale_temps                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/4/13  rhd                *
c     *                                                              *
c     *     save-scale-restore nodal and element temperature values  *
c     *     to support adaptive solutions                            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine mnralg_scale_temps ( process_type, iout )
      use global_data ! old common.main
      use main_data, only : dtemp_nodes, dtemp_elems
      use adaptive_steps, only : adapt_temper_fact ! set in adapt_check
      implicit integer (a-z)
c
      double precision
     &   mag, zero, f
      double precision, dimension (:), save, allocatable ::
     & dtemp_nodes_step, dtemp_elems_step
       data zero /0.0d00/
c
c             MPI:
c               alert MPI slave processors that we are scaling
c               temperatures. send the process_type.
c
      call wmpi_alert_slaves ( 21 )
      call wmpi_bcast_int ( process_type )
      call wmpi_bcast_int ( iout )
c
      select case ( process_type )
      case( 1 )
c
c          save the user-defined increment of nodal and element
c          temperatures for the load step. the dynamically allocated
c          vectors are local (save) to this routine. we use allocation
c          checks just to make sure we've not missed a logic problem
c          in the allocate/deallocate cycle.
c
      if ( allocated( dtemp_nodes_step ) ) then
        write(iout,9100) 1
        call die_abort
        stop
      end if
      allocate( dtemp_nodes_step(nonode) )
      if ( allocated( dtemp_elems_step ) ) then
        write(iout,9100) 2
        call die_abort
        stop
      end if
      allocate( dtemp_elems_step(noelem), stat = alloc_stat )
      if ( alloc_stat .ne. 0 ) then
        write(iout,9100) 5
        call die_abort
        stop
      end if
       dtemp_nodes_step(1:nonode) = dtemp_nodes(1:nonode)
       dtemp_elems_step(1:noelem) = dtemp_elems(1:noelem)
c
      case( 2 )
c
c          we could be doing adaptive solution for step. the vectors
c          dtemp_nodes, dtemp_elems must have the actual temperature
c          increments for this full step or a fraction for the
c          adaptive substep. scale the saved total step increment
c          (above) that was determined before mnralg by eqivlds.
c
        f = adapt_temper_fact
        dtemp_nodes(1:nonode) = f * dtemp_nodes_step(1:nonode)
        dtemp_elems(1:noelem) = f * dtemp_elems_step(1:noelem)
c
      case( 3 )
c
c          release allocated space that saved user-defined temperature
c          increment for step.
c
      if ( allocated( dtemp_nodes_step ) ) then
          deallocate( dtemp_nodes_step )
      else
        write(iout,9100) 3
        call die_abort
        stop
      end if
      if ( allocated( dtemp_elems_step ) ) then
          deallocate( dtemp_elems_step )
      else
        write(iout,9100) 4
        call die_abort
        stop
      end if
c
      case default
c
        write(iout,9100) 10
        call die_abort
c
      end select
      return
c
 9100 format('>> FATAL ERROR: scale_temps reports error: ',i2,
     &     /,'                job aborted' )
      end
c     ****************************************************************
c     *                                                              *
c     *               subroutine mnralg_release_con_forces           *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 2/24/2014 rhd              *
c     *                                                              *
c     *     adjust global load vector if there are user-specified    *
c     *     constrained dof being released. see routine              *
c     *     release_constraints and data structure in mod_main       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine mnralg_release_con_forces
      use global_data ! old common.main
      use main_data, only : release_cons_table, dload
c
      implicit integer (a-z)
c
c           local variables
c
      double precision ::
     &  one, factor, reaction, dload_before, zero
      integer :: inode, idof, sdof, n1, n2
c
      logical local_debug, found
      data zero, one, local_debug / 0.0d00, 1.0d00, .false. /
c
      if( local_debug ) write(out,9000)
c
      if( .not. allocated( release_cons_table ) ) go to 9999
c
      found = .false.
c
c          the release_cons_table exists so there must still be 1
c          or more structural sdof being gradually released from an
c          absolute constraint.
c
c          loop thru table to find entries. compute the force by
c          which the former reaction force should be reduced this step.
c          the fraction < 1.0 and depends on the number of user
c          requested release steps and how many release steps have
c          already occurred.
c
c          during pass thru table, see if there are any more releases
c          to be done. if not, delete the release_cons_table.
c
      if( local_debug ) write(out,9100)
      do inode = 1, nonode
        do idof = 1, 3  ! over u, v, w components
         n1 = release_cons_table(idof,inode)%num_release_steps
         n2 =
     &     release_cons_table(idof,inode)%remaining_steps_for_release
         if( n2 .ne. 0 ) then
           factor = one / real(n1)
           reaction = release_cons_table(idof,inode)%reaction_force
           sdof = 3 * (inode-1) + idof
          load(sdof) = load(sdof) - factor * reaction
             dload_before =  dload(sdof)
             dload(sdof) = dload(sdof) - factor * reaction
           release_cons_table(idof,inode)%remaining_steps_for_release =
     &         n2 - 1
           if( n2 - 1 .gt. 0 ) found = .true.
           if( local_debug )
     &       write(out,9110) inode, idof, n1, n2,
     &                      sdof, reaction, found, factor,
     &                      dload_before, dload(sdof)
          end if
         end do
       end do
c
       if( .not. found ) then
         if( local_debug ) write(out,9120)
         deallocate(release_cons_table)
      end if
c
 9999 continue
      if( local_debug ) write(out,9010)
      return

 9000 format(1x,'.... entered  mnralg_release_con_forces ....' )
 9010 format(1x,'.... leaving  mnralg_release_con_forces ....' )
 9100 format(1x,'inode',3x'idof',2x,'n1',3x,
     & 'n2',2x,'sdof',4x,'reaction',4x,'found',4x,'factor',4x,
     & 'dload before', 4x,'dload after')
 9110 format(1x, i7, i5, i5, i5, i7, e14.6, l5, f10.3, f15.3, f15.3 )
 9120 format(10x,'.... delete release_cons_table on return .....')
c
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine mnralg_ouconv                *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 5/13/2014 rhd              *
c     *                                                              *
c     *     outputs convergence status of global Newton iterations   *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine mnralg_ouconv( itpr, lstitr, step, cnvflg, out,
     &                   total_model_time, adapt_load_fact )
      implicit none
c
      integer :: itpr, lstitr, step, out
      logical :: cnvflg
      double precision ::
     & total_model_time, adapt_load_fact
c
      character (len=10) :: fraction
      double precision,
     &  parameter :: one = 1.0d00
c
c                       newton nonlinear solution algorithm for
c                       a load/time step
c
      if( itpr .ne. 1 ) return
c
      write(fraction,9000) adapt_load_fact

      if( cnvflg ) then
        if( adapt_load_fact .eq. one ) then
          write(out,9100) step, lstitr, total_model_time
        else
          write(out,9102) step, fraction(6:10), lstitr,
     &                    total_model_time
        end if

      else
            write(out,9110) step, lstitr, total_model_time
      end if
c
      return
c
 9000 format(f10.4)
 9100 format(/1x,'>> solution for step: ',i7,
     &           ' converged:',i3,' iters.  ',
     &           ' model time:',e14.6)
 9102 format(/1x,'>> solution for step: ',i7,a4,
     &           ' converged:',i3,' iters.  ',
     &           ' model time:',e14.6)
 9110 format(/1x,'>> solution for step: ',i7,
     &           ' failed to converge:',i2,' iters.  ',
     &           ' model time:',e14.6)
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine stifup                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 10/27/2015 rhd             *
c     *                                                              *
c     *     invoke the process to compute element [K]s               *
c     *                                                              *
c     ****************************************************************
c
      subroutine stifup( step, iter, out, stiff_update_flg,
     &                   show_details )
      implicit none
      integer :: step, iter, out
      logical :: stiff_update_flg, show_details
c
      integer :: local_step, local_iter
c
      stiff_update_flg = .true.
      local_step = step   ! just protect values
      local_iter = iter
      if ( show_details ) write(out,9110) step, iter
      call tanstf( .false., local_step, local_iter )
      return
c
 9110 format(7x,
     & '>> computing tangent stiffness step, iteration:    ',i7,i3)
c
      end
