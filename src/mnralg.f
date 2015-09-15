c     ****************************************************************
c     *                                                              *
c     *                      subroutine mnralg                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified: 9/14/2015 rhd               *
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
c
      use main_data, only : mdiag, pbar, cnstrn, dload,
     &     max_loddat_blks, convergence_history, umat_used,
     &     divergence_check, asymmetric_assembly
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
$add common.main
#dbl      double precision
#sgl      real
     &  mf, mf_nm1, dt_original
      logical mf_ratio_change
      double precision sumrr, sumdd
c
c           local variables
c
#dbl      double precision
#sgl      real
     &  mgres1, magdu1, zero, one, mgload, ext_tol,
     &  scaling_load, ls_s0, ls_s1, mag
c     
      logical dynamic, material_cut_step, emit_extrap_msg,
     & emit_forced_linear_k_for_step,
     & forced_linear_k_for_iter_1, diverging_flag
c
      character * 1 dums
      logical local_debug, first_solve, first_subinc, check_crk_growth,
     &        hypre_solver, cnverg, adaptive, local_mf_ratio_change
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
     &      0.00001d00, .false. /
c
c          predct:          .true. if use requested extrapolate on
c          adaptive_flag    .true. if user requested adaptive solution
c          first_subinc     .true. if we are executing the the first
c                           subincrement after adaptive has reduced
c                           the step size
c          mf_ratio_change  .true. if the equivloads process says this
c                           step load increment is NOT proportional to
c                           step n laod increment
c
      if( local_debug ) then
         write(out,*) " "
         write(out,*) "..... entered mnralg ....."
      end if   
c
      now_step = step
      iout = out
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
      call adapt_check( scaling_adapt, 3, step )
      dt_original  = dt
c
c           calculate the scaling value for extrapolated displacements.
c           this value is used for real steps only (not subincrements).
c           the scaling value is modified if the step is subincremented.
c           mf_ratio_change = .true. means the loading routines believe
c           the next step is not proportional to the last step
c           (then turn off extrapolation for the step and use the
c           linear stiffness on iteration 1)
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
c          adaptive subincrement.
c
c          options for extrapolating (predct) converged displacement
c          increments from the last step. none, (du=0), user
c          defined scaling factor (prdmlt >0), or a computed scaling
c          factor based on the load step factors specified by the user.
c
c          the scaling factor for the step is a function of the
c          scaling factor due to changes in the load mulitplier and
c          the scaling factor due to subincrementation.
c

c          the linear stiffness on iteration 1 is enforced if
c
c            - the user requested it just for this step (but not
c               substeps if adative kicks in)
c            - the loading routines determined that the current step
c              is not proportional to last one (mf_ratio_change
c              = true.). Not for substeps if adaptive kicks in
c            - we just reduced the step size for adaptive. use liner
c              [K] for 1st substep after reduction
c
c          the stiffness update code takes care of the case of
c          when the user requested linear stiffness for iteration
c          1 until further changes by user.
c
      forced_linear_k_for_iter_1 = .false.
      if( linstf_nxt_step ) forced_linear_k_for_iter_1 = .true.
      if( first_subinc ) forced_linear_k_for_iter_1 = .true.
c
      scaling_factor = scaling_load * scaling_adapt
c
      if( local_debug ) write (iout,9600)
     &      scaling_load, scaling_adapt, scaling_factor
c
      if( predct .and. step .gt. 1  ) then
         if( local_mf_ratio_change .or. linstf_nxt_step ) then
           forced_linear_k_for_iter_1 = .true.
           if( emit_extrap_msg ) then
             call errmsg ( 255, dum, dums, dumr, dumd)
             emit_extrap_msg = .false.
         end if
         end if
         if( show_details ) write(out,9152) step, scaling_factor
         if( prdmlt .ne. zero ) then
           du(1:nodof) = prdmlt * du(1:nodof)
         else
           du(1:nodof) =  scaling_factor  * du(1:nodof)
           if(local_debug) write (iout,9530) scaling_factor
         end if
      else
         du(1:nodof) = zero
         if(  local_mf_ratio_change .or. linstf_nxt_step ) then
            forced_linear_k_for_iter_1 = .true.
            if( emit_forced_linear_k_for_step ) then
               emit_forced_linear_k_for_step = .false.
               if( step .ne. 1 ) call errmsg ( 321, dum, 
     &                           dums, dumr, dumd)
            end if
         end if
      end if
c
c          apply constraints to the starting displacement change for
c          step. non-zero constraints are scaled by current adaptive
c          factor (<= 1). the extrapolated displacements are thus
c          scaled back if the load step has been subdivided.
c
      do i = 1, nodof
         if( cstmap(i) .ne. 0 ) du(i) = cnstrn(i) * adapt_disp_fact
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
c          element stiffness matrices. (stifup gets only new element
c          stiffnesss - not a structure stiff).  a
c          linear stiffness on iteration 1 is supported.
c          We're processing stiffness at start of load step (iter=0)
c
      call stifup( step, 0, forced_linear_k_for_iter_1 )
c
c          perform the strain-stress recovery and compute the internal
c          forces to reflect non-zero displacements applied before
c          the step (either extrapolate or non-zero constraints). this is
c          a simple way to get the "applied" forces equivalent to the
c          imposed, incremental displacements for step.
c
      if( predct .or. (.not.zrocon) .or. first_subinc .or.
     &    temperatures .or. umat_used ) then
           if( show_details ) write(out,9155) step
           material_cut_step = .false.
           call drive_eps_sig_internal_forces( step, 0,
     &          material_cut_step )
      end if
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
c          target for top of iteration loop. update stiffness.
c          (stifup gets only new element
c          stiffnesss - not a structure stiff).  a
c          linear stiffness on iteration 1 is supported. right after
c          stifup, we add nodal mass (with newmark beta and dt) to
c          diagonal terms of element stiffnesses.
c
 25   continue
      if( iter .gt. 1 ) ! start of step done above
     &   call stifup( step, iter, .false. )
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
      call eqn_solve( iter, step, first_subinc, first_solve,
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
               call adapt_check( scaling_adapt, 1, step )
               call adaptive_reset
               force_k_update = 2
               first_subinc = .true.
               go to 1000
      end if
c
c          MPI:
c            send all the worker processors information determined from
c            the equation solve to help them continue with the newton
c            iterations.
c
      call wmpi_send_itern
c
c          compute the line search factor s0 for the corrective
c          displacements and the residual loads used to compute
c          them. Line search not implemented. These factors
c          for FYI at this point.
c
      sumrr = zero ; sumdd = zero; ls_s0 = zero
      do i = 1, nodof
         if( cstmap(i) .eq. 0 ) then ! no abs constraint
          ls_s0 = ls_s0 - res(i)*du(i)
          sumrr = sumrr + abs(res(i)); sumdd = sumdd + abs(du(i))
         end if
      end do
      if( local_debug ) write(iout,9710) sumrr, sumdd
c
c          update strains/stresses to define a new estimate of
c          the solution at n+1. compute/assemble a new internal
c          force vector (ifv) using the updates estimate of the
c          stresses and nodal displacements (large displ).
c   
c          a material model may request aborting the current
c          iteration and an immediate adaptive step reduction.
c
      if( show_details ) write(out,9100) step, iter
      call drive_eps_sig_internal_forces( step, iter,
     &                                    material_cut_step )
      if(  material_cut_step ) then ! from a material model
        if ( adaptive  ) then
          if ( show_details ) write(out,9420)
          go to 50
        else
          call abort_job
        end if
      end if
c
c          compute the residual force vector at the end of current
c          ith iteration to advance solution from n -> n+1. this
c          will be driver to compute corrective displacements of
c          iteration i+1.
c          get norm of current total load (applied forces +
c          reactions + inertia - used for convergence check). residuals
c          at constrained dof are set zero. print residuals if
c          requested. compute the line search factor for the residuals
c          after iteration 1.
c
c          compute new contact forces if contact included.
c
      call contact_find
      call upres( iter, out, nodof, dt, nbeta, num_term_loads,
     &            sum_loads, mgload, mdiag, pbar, du, v, a, 
     &            ifv, res, cstmap, dstmap, load )
c
      sumrr = zero; sumdd = zero
      ls_s1 = zero
      do i = 1, nodof
        if( cstmap(i) .eq. 0 ) then ! no abs constraint
          sumrr = sumrr + abs(res(i)); sumdd = sumdd + abs(du(i))
          ls_s1 = ls_s1 - res(i)*idu(i)
        end if
      end do
      if( local_debug ) write(iout,9710) sumrr, sumdd
c
      if( abs(ls_s0) .gt. 1.0d-06 .and. show_details )
     &               write(out,9120) ls_s1/ls_s0
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
      if ( adaptive ) then
          info%adaptive_used = .true.
          call adapt_check( scaling_adapt, 1, step )
          if( adapt_result .eq. 1 ) go to 100
          if( adapt_result .eq. 2 ) then
               call adaptive_reset
               force_k_update = 2
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
          call adapt_check( scaling_adapt, 2, step )
          info%adapt_iters_at_convergence = 
     &        info%adapt_iters_at_convergence + iter
          info%adapt_substeps = 
     &        info%adapt_substeps + 1
      end if
      if( cnverg ) call energy ( step, adapt_result, mdiag )
      first_subinc = .false.
      linstf_nxt_step = .false.
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
      return  ! back to stpdrv
c
 9000 format(/,1x,'... mnralg: step, mf, mf_nm1,'
     & ' mf_ratio_change: ',i6,2f10.6, l2)
 9002 format(/,1x,'... mnralg:  mf_nm1, scaling_load',
     & 2f10.6)
 9100 format(7x,
     & '>> updating strains/stress/forces step, iteration: ',i5,i3)
 9110 format(7x,
     & '>> updating internal forces step, iteration:       ',i5,i3)
 9120 format(7x,
     & '>> line search beta ratio:                         ',f7.3)
 9150 format(/,7x,
     & '>> starting global newton solution for step:       ',i5)
 9152 format(7x,
     & '>> extrapolating displacements to start step:      ',i5,
     & ' (* ',f6.2,')')
 9155 format(7x,
     & '>> stresses/forces for imposed displacements. step:',i5)
 9410 format(/1x,'>> iteration limit exceeded for current step:',i6,
     & /,1x       '   (or adaptive sub-increment) or the solution ',
     &   'appears to be diverging....')
 9420 format(/,2x,'>>> material model stress update computations',
     &     /,2x,'    requested load step reduction...')
 9500 format (/,'... mnralg:',
     & /8x,'Dynamic parameter = ',l2,
     & /8x,'Adaptive solution = ',l2,
     & /8x,'Original time Step (before subincrementing) = ',
     &        e16.6 )
 9510 format(/,'... mnralg:',
     & /,8x,'Original time step = ',e16.6,
     & /,8x,'Adaptive factor = ', e16.6,
     & /,8x,'Current time step = ',e16.6)
 9520 format(7x,
     & '>> automatically forcing a second iteration:       ',i5)
 9530 format(/,'... mnralg:',
     & /,8x,'>>>> Scaling previous displacements by   ',e16.6)
 9600 format(/,'... mnralg:',
     & /,8x,'scaling_load, scaling_adapt, scaling_factor ',
     &             3e16.6,// )
 9610 format(3x,i8,2x,e14.6)
 9700 format(10x,i4, 2x, e14.6, 2x, e14.6)
 9710 format(10x,'sumrr: ',e14.6,'  sumdd: ',e14.6)
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
      end
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
c
      use main_data, only       : du_nm1
      use elem_block_data, only : rot_n_blocks, rot_n1_blocks,
     &                            rts_blocks, rts_nm1_blocks,
     &                            rot_blk_list, rts_blk_list
c
      implicit integer (a-z)
$add common.main
#dbl      double precision
#sgl      real
     &    dummy(1)
c
c                      MPI:
c                        workers run this subroutine for the
c                        data they own.
c
       call wmpi_alert_slaves ( 18 )
c
c                      reset the rts vectors as necessary.
c                      these are deviators of trial elastic stresses
c                      computed at the gauss point during the most
c                      recent stress update. they are needed to
c                      form the consistent tangent matrix. here we
c                      save the just computed values for the converged
c                      load step. note that we only allocate blocks
c                      actually used.
c
c                      we also save the current 3x3 rotation matrices
c                      from the most recent polar decompositions as
c                      the values for start of new step.
c
c                      Note that for each block, we handle the data
c                      only if the list entry for the block is 1, for
c                      instance rts_blk_list(blk)=1.  This provides a simple
c                      mechanism to skip unallocated blocks due either to
c                      processor assignment or material model.
c
      if ( .not. allocated(rts_blocks) ) go to 100
      if ( .not. allocated(rts_nm1_blocks) ) then
        call rts_init( 0, 3 )
      end if
      do blk = 1, nelblk
        if ( rts_blk_list(blk) .eq. 1 ) then
          felem      = elblks(1,blk)
          span       = elblks(0,blk)
          ngp        = iprops(6,felem)
          block_size = span * ngp * nstr
          call vec_ops( rts_nm1_blocks(blk)%ptr(1),
     &                  rts_blocks(blk)%ptr(1), dummy, block_size, 5 )
        end if
      end do
c
 100  continue
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
c
      use main_data,       only : du_nm1, nonlocal_analysis
      use elem_block_data, only : history_blocks, history1_blocks,
     &                            rot_n_blocks, rot_n1_blocks,
     &                            rts_blocks, rts_nm1_blocks,
     &                            eps_n_blocks, eps_n1_blocks,
     &                            urcs_n_blocks, urcs_n1_blocks,
     &                            history_blk_list, rot_blk_list,
     &                            rts_blk_list, eps_blk_list,
     &                            urcs_blk_list, nonlocal_flags,
     &                            nonlocal_data_n, nonlocal_data_n1
      implicit integer (a-z)
$add common.main
c
c                      lcoal values
c
#dbl      double precision
#sgl      real
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
c                       3) rts vector for mm01, mm03. these are
c                          deviators of trial elastic stresses
c                          computed at the gauss point during the most
c                          recent stress update. they are needed to
c                          form the consistent tangent matrix. here we
c                          reload the values present at the start
c                          of this load step
c                       4) strains
c                       5) unrotated cauchy stresses
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
        rts_size   = span * ngp * nstr
        if ( hist_size .gt. 0 )
     &      call vec_ops( history1_blocks(blk)%ptr(1),
     &                    history_blocks(blk)%ptr(1), dummy,
     &                    hist_size, 5 )
        if ( rot_blk_list(blk) .eq. 1 )
     &      call vec_ops( rot_n1_blocks(blk)%ptr(1),
     &                    rot_n_blocks(blk)%ptr(1), dummy,
     &                    rot_size, 5 )
        if ( rts_blk_list(blk) .eq. 1 )
     &      call vec_ops( rts_blocks(blk)%ptr(1),
     &                    rts_nm1_blocks(blk)%ptr(1), dummy,
     &                    rts_size, 5 )
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
      implicit integer (a-z)
$add common.main
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
c     *                   last modified : rhd 4/19/2015 rhd          *
c     *                                                              *
c     *     solution of the linearized equilibrium equations         *
c     *     using an available equation solver.                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine eqn_solve( iter, step, first_subinc, first_solve,
     &                      nodof, solver_flag, use_mpi, show_details,
     &                      du, idu, iout ) 
c
      use hypre_parameters, only: hyp_trigger_step
      use stiffness_data, only: d_lagrange_forces, 
     &                          i_lagrange_forces
c
      implicit none
c
c          parameters
c
      integer :: iter, step, solver_flag, iout, nodof 
      logical :: first_subinc, first_solve, use_mpi, show_details
#dbl      double precision ::
#sgl      real ::
     &  du(nodof), idu(nodof)            
c
c          locals
c
      logical :: hypre, not_hypre
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
      not_hypre = .not. hypre
c
c          MPI:
c            1) hypre - just solve
c            2) Pardiso - we still allow MPI + Pardiso. Workers
c               perform element level computations. We bring 
c               [Ke]s back to root for assembly & solve.
c               Put MPI procs to sleep while Pardiso runs.
c
      if( use_mpi .and. not_hypre ) call wmpi_suspend( 1 )
c
c          drive the equation solver. Code to allow changing
c          solver during a run has been removed.
c
      if ( not_hypre .and. show_details ) write(iout,9142) step, iter
      call direct_routine_sparse( first_solve, iter)  ! run solve
      first_solve = .false.
c
c          MPI: wake up worker processes
c
      if( use_mpi .and. not_hypre ) call wmpi_suspend( 2 )
c
c          check if hypre wants an adaptive step.
c          We can only deal with this in the main mnralg, so
c          simply keep going up the stack.
c
      if( hyp_trigger_step ) return
c
c          update after this iteration:
c
c           - estimated change in displacement from n -> n+1
c           - estimated change in Lagrange forces from n -> n+1
c                for MPCs (always self-equilibrating)
c           
      du(1:nodof) = du(1:nodof) + idu(1:nodof)
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
     & '>> running sparse solver step, iteration:          ',i5,i3)
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
      use main_data, only : dtemp_nodes, dtemp_elems
      use adaptive_steps, only : adapt_temper_fact ! set in adapt_check
      implicit integer (a-z)
$add common.main
c
#dbl      double precision
#sgl      real
     &   mag, zero, f
#dbl      double precision, dimension (:), save, allocatable ::
#sgl      real, dimension (:), save, allocatable ::
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
      use main_data, only : release_cons_table, dload
c
      implicit integer (a-z)
$add common.main
c
c           local variables
c
#dbl      double precision ::
#sgl      real ::
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
 9110 format(1x, i6, i5, i5, i5, i6, e14.6, l5, f10.3, f15.3, f15.3 )
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
#dbl      double precision ::
#sgl      real ::
     & total_model_time, adapt_load_fact
c
      character (len=6) :: fraction
#dbl      double precision,
#sgl      real,
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
          write(out,9102) step, fraction(3:6), lstitr,
     &                    total_model_time 
        end if
        
      else
            write(out,9110) step, lstitr, total_model_time
      end if
c
      return
c
 9000 format(f6.3)
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
