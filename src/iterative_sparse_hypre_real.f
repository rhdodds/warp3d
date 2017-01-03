c     ***********************************************************************
c     *                                                                     *
c     *     Actual routines for the hypre solver library                    *
c     *     Includes local_stiffness_mod, used to distribute rows (equally) *
c     *     to processors from assembled stiffness.                         *
c     *     Note: this means this file uses MPI calls.                      *
c     *                                                                     *
c     ***********************************************************************
c
c
c     **********************************************************************
c     *                                                                    *
c     *     iterative_sparse_hypre                                         *
c     *                                                                    *
c     *     All ranks enter this, calls conversion                         *
c     *     and distribution functions.  Then it calls the                 *
c     *     hypre functions to setup the hypre data structs.  Finally it   *
c     *     alls a separate function to actually solve the system.         *
c     *                                                                    *
c     *     On return it gathers the solution vector to root and           *
c     *     frees the memory that can/should be freed                      *
c     *                                                                    *
c     *     written by: mcm 2/11                                           *
c     *     modified: mcm 5/11                                             *
c     *                                                                    *
c     **********************************************************************
c
       subroutine iterative_sparse_hypre(sol_vec,out,error_code )
c
      use distributed_stiffness_data
      use local_stiffness_mod
      use hypre_parameters
      use performance_data
      implicit none    
      include 'HYPREf.h'
c
c           input
      integer :: out, error_code
      double precision :: sol_vec
      dimension :: sol_vec(*)
c
c           Local variables
      logical :: debug
      integer :: ierr, i,j
      real :: wcputime
      external :: wcputime
c           Hypre things - I make a note of this because they are integer*8
c           used as pointers.
      integer*8 :: hypre_mat
      integer*8 :: hypre_par_mat
      integer*8 :: hypre_rhs
      integer*8 :: hypre_par_rhs
      integer*8 :: hypre_soln
      integer*8 :: hypre_par_soln
      integer :: hypre_err
c
c
      debug = .false.
c
      if (debug) write (*,'("Rank ",i3," entering hypre soln fn")') 
     &             local_k%my_rank
c
      if (hyp_first_solve) then
            error_count = 0
      end if
c
c           Insert to HYPRE structs.  Refer to the hypre manual if you
c           want to know the details of each call.
c           Matrix
      call HYPRE_IJMatrixCreate(local_k%comm,local_k%local_start,
     &      local_k%local_start+local_k%local_n-1, local_k%local_start,
     &      local_k%local_start+local_k%local_n-1, hypre_mat
     &      , hypre_err)
      call HYPRE_IJMatrixSetObjectType(hypre_mat,HYPRE_PARCSR,
     &      hypre_err)
      call HYPRE_IJMatrixSetRowSizes(hypre_mat,local_k%cols_per_row,
     &      hypre_err)
      call HYPRE_IJMatrixInitialize(hypre_mat,hypre_err)
      call HYPRE_IJMatrixSetValues(hypre_mat,local_k%local_n,
     &      local_k%cols_per_row,local_k%vec_indexes, 
     &      local_k%col_indexes,local_k%coefs, hypre_err)
      call HYPRE_IJMatrixAssemble(hypre_mat,hypre_err)
      call HYPRE_IJMatrixGetObject(hypre_mat,hypre_par_mat,hypre_err)
c
c           RHS
      call HYPRE_IJVectorCreate(local_k%comm, local_k%local_start,
     &      local_k%local_start+local_k%local_n-1, hypre_rhs, hypre_err)
      call HYPRE_IJVectorSetObjectType(hypre_rhs,HYPRE_PARCSR,hypre_err)
      call HYPRE_IJVectorInitialize(hypre_rhs,hypre_err)
      call HYPRE_IJVectorSetValues(hypre_rhs,local_k%local_n,
     &      local_k%vec_indexes, local_k%rhs, hypre_err)
      call HYPRE_IJVectorAssemble(hypre_rhs,hypre_err)
      call HYPRE_IJVectorGetObject(hypre_rhs,hypre_par_rhs,hypre_err)
c
c           Solution vector
      call HYPRE_IJVectorCreate(local_k%comm,local_k%local_start,
     &      local_k%local_start+local_k%local_n-1, hypre_soln,
     &       hypre_err)
      call HYPRE_IJVectorSetObjectType(hypre_soln, HYPRE_PARCSR,
     &      hypre_err)
      call HYPRE_IJVectorInitialize(hypre_soln,hypre_err)
      call HYPRE_IJVectorSetValues(hypre_soln,local_k%local_n,
     &      local_k%vec_indexes,local_k%soln, hypre_err)
      call HYPRE_IJVectorAssemble(hypre_soln, hypre_err)
      call HYPRE_IJVectorGetObject(hypre_soln,hypre_par_soln,hypre_err)
c
      if (local_k%my_rank .eq. 0) then
            write(out,
     &      '(15x,"HYPRE structs created         @ ", f10.2)')
     &      wcputime(1)
      end if
c
c           Enter our solution routine
      call solve_hypre(hypre_par_mat,hypre_par_soln,hypre_par_rhs,
     &                              error_code, local_k%comm, out)
c           Output HYPRE to fortran vector
c           First take it back to fortran vectors, then gather
c           it to root
      call HYPRE_IJVectorGetValues(hypre_soln,local_k%local_n,
     &      local_k%vec_indexes,local_k%soln,hypre_err)
      call merge_soln(local_k, sol_vec)
      call MPI_Barrier(local_k%comm,ierr)
c
c           Free HYPRE structs
      call HYPRE_IJMatrixDestroy(hypre_mat,hypre_err)
      call HYPRE_IJVectorDestroy(hypre_rhs, hypre_err)
      call HYPRE_IJVectorDestroy(hypre_soln,hypre_err)
c
c           We may as well free our allocated data as we cannot
c           keep it between solves (see above comment on hypre-driver)
c      call clear_dist(local_k)
c
      if (local_k%my_rank .eq. 0) then
            write(out,
     &      '(15x,"Solver finished               @ ",f10.2)')
     &      wcputime(1)
      end if
c
c      
c      
      return
      end 
c
c
c     **********************************************************************
c     *                                                                    *
c     *     solve_hypre                                                    *
c     *                                                                    *
c     *     Solves Ax=b using various hypre solvers.                       *
c     *                                                                    *
c     *     Must provide handles to hypre data structs.  (HANDLES, not the *
c     *     actual struct.  See hypre documentation)                       *
c     *                                                                    *
c     *     written by: mcm 2/11                                           *
c     *     modified: mcm 5/11                                             *
c     *                                                                    *
c     **********************************************************************
c
      subroutine solve_hypre(A,x,b, error_code, usecomm, out)
      use hypre_parameters
      implicit none    
      include 'HYPREf.h'
      include 'HYPRE_error_f.h'
      include 'mpif.h'
c
c           Input
      integer*8 :: A, x, b
      integer :: error_code, usecomm, out
c
c           Local variables
      integer*8 :: precond, solver
      integer :: hypre_err, ierr
      integer :: rank, procs
      Integer :: total_iters, zero_err
      real :: wcputime
      external :: wcputime

c
      hypre_err = 0
c
c           Get MPI info and sync all the top level (i.e. which solver/precond
c           we want to use) parameters (as root only will have the actual
c           values).
      call MPI_Comm_rank(usecomm, rank, ierr)
      call MPI_Comm_size(usecomm, procs, ierr)
      call MPI_Bcast(precond_type,1,MPI_INTEGER,0,usecomm,
     &      ierr)
      call MPI_Bcast(hsolver_type,1,MPI_INTEGER,0,usecomm,
     &      ierr)
      call MPI_Bcast(precond_printlevel,1,MPI_INTEGER,0,usecomm,
     &       ierr)
      call MPI_Bcast(solver_printlevel,1,MPI_INTEGER,0,usecomm,
     &       ierr)
c
c
c           Create and setup the preconditioner.  Branch on the preconditioner
c           type flag:
c
c           1 - Parasails
c           2 - BoomerAMG (not implemented)
c
      if (precond_type .eq. 1) then
            call HYPRE_ParaSailsCreate(usecomm, precond,
     &            hypre_err)
c
c                 sync parameters
            call MPI_Bcast(threshold,1,MPI_DOUBLE_PRECISION,0,
     &            usecomm, ierr)
            call MPI_Bcast(levels,1,MPI_INTEGER,0,
     &            usecomm,ierr)
            call MPI_Bcast(symme,1,MPI_INTEGER,0,
     &            usecomm,ierr)
            call MPI_Bcast(loadbal,1,MPI_DOUBLE_PRECISION,0,
     &            usecomm,ierr)
            call MPI_Bcast(filter,1,MPI_DOUBLE_PRECISION,0,
     &            usecomm,ierr)
c
            if (levels .eq. 0) then
                  levels = 1
            end if
c                 Set the params
            call HYPRE_ParaSailsSetParams(precond,threshold,levels,
     &            hypre_err)
            call HYPRE_ParaSailsSetSym(precond,1,hypre_err)
            call HYPRE_ParaSailsSetLoadbal(precond,loadbal,hypre_err)
            call HYPRE_ParaSailsSetFilter(precond, filter,hypre_err)
            call HYPRE_ParaSailsSetLogging(precond,precond_printlevel,
     &                  hypre_err)
      else if (precond_type .eq. 2) then
c           Sync parameters
            call MPI_Bcast(mg_threshold,1,MPI_DOUBLE_PRECISION,0,
     &            usecomm, ierr)
            call MPI_Bcast(max_levels,1,MPI_INTEGER,0,
     &            usecomm,ierr)
            call MPI_Bcast(interpolation,1,MPI_INTEGER,0,
     &            usecomm,ierr)
            call MPI_Bcast(coarsening,1,MPI_INTEGER,0,
     &            usecomm,ierr)
            call MPI_Bcast(agg_levels,1,MPI_INTEGER,0,
     &            usecomm,ierr)
            call MPI_Bcast(truncation,1,MPI_DOUBLE_PRECISION,0,
     &            usecomm, ierr)
            call MPI_Bcast(relaxation,1,MPI_INTEGER,0,
     &            usecomm,ierr)
            call MPI_Bcast(relax_wt,1,MPI_DOUBLE_PRECISION,0,
     &            usecomm, ierr)
            call MPI_Bcast(relax_outer_wt,1,MPI_DOUBLE_PRECISION,0,
     &            usecomm, ierr)
            call MPI_Bcast(sweeps,1,MPI_INTEGER,0,
     &            usecomm,ierr)
            call MPI_Bcast(cf,1,MPI_INTEGER,0,
     &            usecomm,ierr)
            call MPI_Bcast(cycle_type,1,MPI_INTEGER,0,
     &            usecomm,ierr)
c           Set up and set params
            call HYPRE_BoomerAMGCreate(precond,
     &            hypre_err)
            call HYPRE_BoomerAMGSetMaxLevels(precond, max_levels, 
     &            hypre_err)
            call HYPRE_BoomerAMGSetStrongThrshld(precond, mg_threshold,
     &            hypre_err)
            call HYPRE_BoomerAMGSetCoarsenType(precond, coarsening,
     &            hypre_err)
            call HYPRE_BoomerAMGSetCycleType(precond,cycle_type,
     &            hypre_err)
            call HYPRE_BoomerAMGSetNumSweeps(precond, sweeps, hypre_err)
            call HYPRE_BoomerAMGSetRelaxType(precond, relaxation,
     &            hypre_err)
            call HYPRE_BoomerAMGSetRelaxOrder(precond, cf, hypre_err)
            call HYPRE_BoomerAMGSetRelaxWt(precond, relax_wt, hypre_err)
            call HYPRE_BoomerAMGSetOuterWt(precond, relax_outer_wt,
     &            hypre_err)
            call HYPRE_BoomerAMGSetTruncFactor(precond, truncation,
     &            hypre_err)
            call HYPRE_BoomerAMGSetInterpType(precond, interpolation, 
     &            hypre_err)
            call HYPRE_BoomerAMGSetAggNumLevels(precond, agg_levels,
     &            hypre_err)
            call HYPRE_BoomerAMGSetPrintLevel(precond, 
     &            precond_printlevel, hypre_err)


      else
c           Input reader should not let you get here, but just in case...
            if (rank .eq. 0) then
                  write (out,*) "Error: Unknown preconditioner."
            end if
            call MPI_Abort(usecomm,2,ierr)
      end if
c
c           Now branch on solver type, setup, and solve.  Solver types:
c     
c           1 - parallel PCG
c     
      if (hsolver_type .eq. 1) then
            call HYPRE_ParCSRPCGCreate(usecomm,solver,hypre_err)
c           Sync and set the parameters
            call MPI_Bcast(hypre_tol,1,MPI_DOUBLE_PRECISION,0,
     &            usecomm,ierr)
            call MPI_Bcast(hypre_max, 1, MPI_INTEGER,0,
     &            usecomm,ierr)
            call HYPRE_ParCSRPCGSetMaxIter(solver,hypre_max,hypre_err)
            call HYPRE_ParCSRPCGSetTol(solver, hypre_tol, hypre_err)
            call HYPRE_ParCSRPCGSetPrintLevel(solver,solver_printlevel,
     &                  hypre_err)
c
c                 Set preconditioner.  Note hypre uses codes here, they are:
c                 1 - DS? (can't find)
c                 2 - BoomerAMG
c                 3 - Pilut
c                 4 - ParaSails
            if (precond_type .eq. 1) then
c                 parasails
                  call HYPRE_ParCSRPCGSetPrecond(solver,4,precond,
     &                  hypre_err)
            else if (precond_type .eq. 2) then
c                 boomerAMG
                  call HYPRE_ParCSRPCGSetPrecond(solver,2,precond,
     &                  hypre_err)
            else
c                 There is no way to be here, we call errors above
            end if
c
c                 Setup and get our error code
            call HYPRE_ParCSRPCGSetup(solver,A,b,x,hypre_err)
c
c                 A non-zero error code implies that the solver failed to setup
c                 the preconditioner.  This commonly means that we need an
c                 adaptive step.  Note we count the number of
c                 adaptive calls we make in a row as we are unlikely
c                 to recover from 2 in a row.
c           
c                 If we do need a step, set the flag accordingly, deallocate
c                 and return.
c
c                 BoomerAMG seems to think it fails even though it actually
c                 doesn't.  Take this into account (error code 256)
            if (((precond_type .eq. 1) .and. (hypre_err .gt. 0)) .or.
     &                 ((precond_type .eq. 2) .and. (hypre_err .gt. 0)
     &                  .and. (hypre_err .ne. 256))) then
                  precond_fail_count=precond_fail_count+1
                  if (rank .eq. 0) then
                        error_count = error_count+1
                        write (*,*) 
     &                  ">>> Preconditioning failed."
                        write (*,*) hypre_err
                        error_code = 1
                  end if
                  if (precond_type .eq. 1) then
                        call HYPRE_ParaSailsDestroy(precond,hypre_err)
                  else if (precond_type .eq. 2) then
                        call HYPRE_BoomerAMGDestroy(precond,hypre_err)
                  end if
                  call HYPRE_ParCSRPCGDestroy(solver,hypre_err)
                  return
            end if
c          
            if (rank .eq. 0) then
                  write(out,
     &            '(15x,"Preconditioner created        @ ",f10.2)')
     &            wcputime(1)
            end if
c                 Solve and display stats.
            call HYPRE_ParCSRPCGSolve(solver,A,b,x,hypre_err)
            if (rank .eq. 0) then
                  write(out,
     &            '(15x,"System solved                 @ ",f10.2)')
     &            wcputime(1)
            end if
            call HYPRE_ParCSRPCGGetNumIterations(solver,
     &            total_iters, hypre_err)
            if (rank .eq. 0) then
                  write(out,
     &            '(15x,"Iterations                        ",i8)')
     &            total_iters
            end if
      else
c                 If we had other solver types, here is where they
c                 would be.
            if (rank .eq. 0) then
                  write (*,*) "Unknown solver type"
            end if
            call MPI_Abort(usecomm,2,ierr)
      end if
c           Clean up
      if (precond_type .eq. 1) then
            call HYPRE_ParaSailsDestroy(precond,hypre_err)
      else if (precond_type .eq. 2) then
            call HYPRE_BoomerAMGDestroy(precond,hypre_err)
      end if
      if (hsolver_type .eq. 1) then
            call HYPRE_ParCSRPCGDestroy(solver,hypre_err)
      end if
c           Reset our error code and "number of failures in a row count"
      if (rank .eq. 0) then
            error_code = 0
      end if
      precond_fail_count = 0
c
      return
      end subroutine solve_hypre
