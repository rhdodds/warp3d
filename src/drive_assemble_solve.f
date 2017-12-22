c
c     ****************************************************************
c     *                                                              *
c     *  assemble & solve linear equations for a Newton iteration    *
c     *                                                              *
c     *                       written by  : rhd                      *
c     *                   last modified : 1/5/2017 rhd               *
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_assemble_solve( first_solve, now_iteration,
     &                                 suggested_new_precond )
      use global_data ! old common.main
c
      use elem_block_data, only : edest_blocks
      use main_data, only : repeat_incid, modified_mpcs,
     &                      asymmetric_assembly
      use stiffness_data, only : ncoeff, big_ncoeff, k_coeffs,
     &                           k_indexes,
     &                           ncoeff_from_assembled_profile
      use mod_mpc, only : tied_con_mpcs_constructed, mpcs_exist
      use hypre_parameters, only: precond_fail_count, hyp_trigger_step
      use performance_data
      use distributed_stiffness_data, only: parallel_assembly_used,
     &                                   distributed_stiffness_used
c
      implicit none
c
c                    parameter declarations
c
      integer :: now_iteration
      logical :: first_solve, suggested_new_precond
c
c                    locals
c
      integer :: i, k, error_code, nnz, mkl_threads, num_enode_dof,
     &           num_struct_dof, iresult, itype, num_terms, i_offset,
     &    eqn_num
      integer, save :: neqns, old_neqns, old_ncoeff
      integer, external :: curr_neqns
      integer :: code_vec(mxedof,mxvl), edest(mxedof,mxconn)
      integer, allocatable :: k_ptrs(:)
      integer, allocatable, save :: dof_eqn_map(:), eqn_node_map(:),
     &                              save_k_indexes(:), save_k_ptrs(:)

      logical :: new_size
      logical, save :: cpu_stats, save_solver
      logical, parameter :: local_debug = .false.,
     &                      local_debug2 = .false.
c
      real, external :: wcputime
      double precision, parameter :: zero = 0.0d00
      double precision, allocatable :: p_vec(:), k_diag(:), u_vec(:)
c
      data old_neqns, old_ncoeff, cpu_stats, save_solver
     &     / 0, 0, .true., .false. /
c
      if( local_debug ) write(*,*) '... drive_assem_solve ... @ 1'
      if( .not. show_details ) cpu_stats = .false.
c
      num_enode_dof  = iprops(4,1) ! always = 3 in warp3d
      num_struct_dof = nonode * num_enode_dof
      new_size       = .false.
c
c              stiffness data structure development. can skip
c              if not 1st solution for this load(time) step.
c              can also skip if # eqns and sparsity unchanged
c              from prior load(time) step.

      if( first_solve ) then ! 1st for this load(time) step
c
c           1. build sparsity data structures
c              ------------------------------
c
c              based on decisions/checks made in mnralg,
c              setup for threaded or MPI parallel assembly by
c              generating the sparsity data structures
c
        call t_start_assembly( start_assembly_step ) ! timing
        if( .not. parallel_assembly_used ) then
          call ds_setup_sparsity_local( iresult, use_mpi, numprocs )
          if( iresult .eq. 1 ) then ! no equations to solve
             idu(1:nodof) = zero
             return
          end if
        end if
        if( local_debug ) write(*,*) '... drive_assem_solve  @ 2'
        if( parallel_assembly_used ) then
          call ds_setup_sparsity_distributed( iresult )
          if( iresult .eq. 0 ) then  ! no equations to solve
             idu(1:nodof) = zero
             return
          end if
        end if
        call t_end_assembly( assembly_total, start_assembly_step )
        if( local_debug ) write(*,*) '... drive_assem_solve  @ 3'

c
      end if   ! on first_solve
c
c          2.  assemble the equations and rh side here
c              ---------------------------------------
c
c             => threads only assembly
c
c                - symmetric assembly with sparse (MKL)
c                  direct/iterative solvers. Can have MPC
c                  from user or tied-contact.
c
c                - asymmetric assembly with sparse (MKL)
c                  direct/iterative solvers. The sparsity is symmetric
c                  but not coefficients, e.g., from crystal plasticity
c                  material model
c
c                - a symmetric/asymmetric assembly on root followed by
c                  distribution to MPI ranks for solution by hypre
c
c                or
c
c                - a symmetric/asymmetric assembly on root follwed
c                  by use of MPI cpardiso direct solver on ranks
c
c             => distributed assembly across MPI ranks followed by
c                hypre or CPardiso solve.
c
      call t_start_assembly( start_assembly_step )
      ntimes_assembly = ntimes_assembly + 1
      if( .not. parallel_assembly_used ) then
        call ds_root_assembly
        itype = 2
        if( new_size ) itype = 1
       end if
      if( parallel_assembly_used ) call ds_distributed_assembly
c
c              distributed assembly across MPI ranks followed by
c              do we need a distributed stiffness -- but did not use
c              parallel assembly?  If so, take care of that now.
c              The "true" flag indicates full (not LT) assembly and the
c             'blockrow' parameter means agglomerate by rows and map
c              by blocks of equal nnz size.
c
      if( .not. parallel_assembly_used
     &    .and. distributed_stiffness_used ) then
            call wmpi_alert_slaves(39)
            call distribute_from_assembled(neqns, ncoeff,
     &                  k_diag, p_vec, u_vec, k_coeffs, k_ptrs,
     &                  k_indexes, .true., "blockrow", itype, out )
      end if
c
      if( local_debug ) write(*,*) '... drive_assem_solve @ 4'
      call t_end_assembly( assembly_total, start_assembly_step )
c
c
c          3.  solve the equations
c              -------------------
c
c              computed u_vec has the corrections to du.
c              du is current estimate of total displacement
c              change from n -> n+1
c
c              p_vec was set above and is the set of residual
c              nodal forces at n+1. they include MPC effects as
c              needed for symmetric MKL solvers.
c
c               hypre_solver        => solver_flag .eq. 9
c               pardiso_asymmetric  => solver_flag .eq. 8
c               pardiso_symmetric   => solver_flag .eq. 7
c               cpardiso_symmetric  => solver_flag .eq. 10
c               cpardiso_asymmetric => solver_flag .eq. 11
c
      select case( solver_flag )

      case( 7 ) ! symmetric pardiso
c
        if( asymmetric_assembly ) then
          write(out,9120); call die_gracefully
        end if
        if( local_debug ) write(*,*) '... drive_assem_solve  @ 5a'
        call pardiso_symmetric( neqns, ncoeff, k_diag, p_vec,
     &                          u_vec, k_coeffs, k_ptrs, k_indexes,
     &                          cpu_stats, itype, out,
     &                          solver_out_of_core, solver_memory,
     &                          solver_scr_dir, solver_mkl_iterative,
     &                          suggested_new_precond )
        if( local_debug ) write(*,*) '... drive_assem_solve @ 5b'

      case( 8 ) ! asymmetric pardiso
c
        if( .not. asymmetric_assembly ) then
          write(out,9110); call die_gracefully
        end if
        if( local_debug ) write(*,*) '... drive_assem_solve  @ 5c'
        call pardiso_unsymmetric( neqns, nnz, k_ptrs, k_indexes,
     &            k_coeffs,  p_vec, u_vec, cpu_stats, itype, out,
     &            solver_mkl_iterative )
        if( local_debug ) write(*,*) '... drive_assem_solve @ 5d'
c

      case( 9 ) ! hypre_solver
c
        call ds_drive_hypre
        if( hyp_trigger_step ) return ! failed and wants Newton
c
      case( 10 ) ! symmetric cpardiso
c
        if( asymmetric_assembly ) then
          write(out,9120); call die_gracefully
        end if
        if( local_debug ) write(*,*) '... drive_assem_solve @ 5e'
        call wmpi_alert_slaves( 3 )
        call cpardiso_symmetric( neqns, ncoeff, k_diag, p_vec,
     &                          u_vec, k_coeffs, k_ptrs, k_indexes,
     &                          cpu_stats, itype, out, myid )
        if( local_debug ) write(*,*) '... drive_assem_solve @ 5f'
c
      case( 11 ) ! asymmetric cpardiso
c
        if( .not. asymmetric_assembly ) then
             write(out,9125); call die_gracefully
        end if
        if( local_debug ) write(*,*) '... drive_assem_solve @ 5g'
        call wmpi_alert_slaves( 31 )
        call cpardiso_unsymmetric( neqns, nnz, k_ptrs, k_indexes,
     &            k_coeffs,  p_vec, u_vec, cpu_stats, itype, out,
     &            myid )
        if( local_debug ) write(*,*) '... drive_assem_solve @ 5h'
c
      case default ! bad solver type
c
            write(out,9301) solver_flag
            call die_gracefully
c
      end select
c
      if( local_debug ) write(out,*) ' @ 6'
      if( .not. parallel_assembly_used )
     &  deallocate( k_diag, k_coeffs, k_ptrs, k_indexes, p_vec )
c
c          4.  for distributed assembly/solve, reorder solution vector
c              -------------------------------------------------------
c
c              to eliminate the effects of the remapping
c              to improve load balancing
c
      call t_start_assembly( start_assembly_step )
      call reorder_soln_vec( u_vec )
      call t_end_assembly( assembly_total, start_assembly_step )
c
c          5.  apply mpc equations
c              -------------------
c
c              if they exist to compute the corrective displacements
c              to du at dependent dofs.
c
c              compute corrections estimate of increment of
c              lagrange multipliers over the step
c
c              these are nodal forces at all MPC dof to enforce the
c              multi-points constraints
c
c              only for symmetric systems. put into their
c              positions in u_vec.
c
      if( tied_con_mpcs_constructed .or. mpcs_exist ) then
         if( local_debug ) write(out,*) ' @ 7'
         call mpcs_apply( u_vec, neqns, nodof, cstmap )
      end if
c
c
c          6.  expand compressed u_vec to structure nodof
c              ------------------------------------------
c
c              u_vec from solver is the reduced number of equations
c              from absolute constraints and MPCs. expand to
c              nodof size for structure.
c
c              Again - these are the corrective displacements to
c              be added to du.
c
      do i = 1, nodof
       if( dof_eqn_map(i) .eq. 0 ) then
          idu(i) = zero
       else
          idu(i) = u_vec(dof_eqn_map(i))
       end if
      end do
c
      if( local_debug2 ) write(out,9900) (i, u_vec(i), i=1, neqns)
c
      deallocate( u_vec ) ! scratch vec used by solvers
c
c          7.  re-define number of non-zero coefficients in the
c              equations as solved.
c              -------------------------------------------------
c
c              needed only if we inserted mpcs. The next time thru
c              the sover this now "correct" count of the number of
c              terms. saves a lot of work during re-insertion of the
c              mpc coefficients.
c
      if( tied_con_mpcs_constructed .or. mpcs_exist )
     &    ncoeff = big_ncoeff
c
      cpu_stats = .false.   ! only print them 1st time thru solver
c
      return
c
 9100 format('>> constrained load vector: ',
     & 1000(/3x,8f12.6) )
 9110 format(1x,'>> FATAL ERROR: Job Aborted.',
     & /,5x,'asymmetric solver requested for symmetric equations')
 9115 format(1x,'>> FATAL ERROR: Job Aborted.',
     & /,5x,'asymmetric solver requested for symmetric equations')
 9120 format(1x,'>> FATAL ERROR: Job Aborted.',
     & /,5x,'symmetric solver requested for asymmetric equations')
 9125 format(1x,'>> FATAL ERROR: Job Aborted.',
     & /,5x,'asymmetric solver requested w/o asymmetric assembly')
 9301 format(1x,
     &'>> FATAL ERROR: Job Aborted.',
     &5x,' Requested equation solver cannot be used',
     &/5x,'Solver_flag is set to:',i12 )
 9602 format('>> Terminating analysis at this point',//)
 9900 format('>> incremental displacements: ',
     & 1000(/4x,i8, f12.8) )
c
c
       contains
c      ========
c
c
c     ****************************************************************
c     *                                                              *
c     *   set up equation sparsity for local assembly: symmetric     *
c     *                                                              *
c     *                       written by  : rhd                      *
c     *                   last modified : 6/6/2017 rhd               *
c     *                                                              *
c     ****************************************************************
c
      subroutine ds_setup_sparsity_local( ireturn, using_mpi,
     &                                    num_ranks )
      implicit none
c
      integer :: ireturn, num_ranks  ! local to this routine
      logical :: using_mpi
c
c                    locals
c
      integer :: i, j, k, l, srow, start_srow
      integer, allocatable :: start_kindex_locs(:), edest(:,:,:),
     &                        scol_flags(:,:), scol_lists(:,:)
      character(len=200) mkl_string
      integer :: mkl_num_thrds, next_space, count_previous, count_now,
     &           nrow_lists, safety_factor
      integer, external :: mkl_get_max_threads
      logical, parameter :: local_debug_2 = .false.
c
c              1. generate equation numbers for all unconstrained
c                 nodal dof. dof_eqn_map(i) gives equation
c                 number for structure dof number. eqn_node_map
c                 gives the node number for an equation number.
c                 set the actual number of equations we will solve
c                 (neqns). if this turns out = 0, all nodes in
c                 model are fully constrained and we just leave.
c                 the solution is known from constraints.
c
c
      if( local_debug ) write(out,*) '.. @ (1)  '
      call thyme( 21, 1 )
      if( cpu_stats .and. show_details ) then
          call mkl_get_version_string( mkl_string )
          mkl_num_thrds = mkl_get_max_threads()
          if( using_mpi ) then
             write(out,9300) num_threads, num_ranks, mkl_string(38:45),
     &                       mkl_string(61:68), wcputime(1)
           else
             write(out,9400) num_threads, mkl_string(38:45),
     &                       mkl_string(61:68), wcputime(1)
          end if
      end if
      if( .not. allocated ( dof_eqn_map ) ) then
          allocate( dof_eqn_map(num_struct_dof) )
          allocate( eqn_node_map(num_struct_dof) )
      end if
      call dof_map( dof_eqn_map, cstmap, nonode, num_enode_dof,
     &              eqn_node_map, neqns )
      if( cpu_stats .and. show_details ) write(out,9409) wcputime(1)
c
      ireturn = 1
      if( neqns .le. 0 ) return
c
c              2. decide if this set of equations has a different
c                 sparsity from the previous time we came thru here.
c                 right now, this is based on the actual number
c                 of equations. in the future, this test may need
c                 to actual verify the same sparsity structure, e.g.,
c                 the same numberof equations and the same
c                 number of off diagonal terms.
c
      if( neqns .ne. old_neqns )  new_size = .true.
      if( local_debug ) write(out,*) '....  neqns, old_neqns: ',
     &                             neqns, old_neqns
      old_neqns  = neqns
      ireturn = 2
      if( .not. new_size ) return
c
c              3. compute the number of non-zero terms in the
c                 upper-triangle of the assembled structure stiffness.
c                 this is done by simulating a row-by-row
c                 sparse format assembly (rather than an element-by-
c                 element assembly). this way we never need to allocate
c                 memory for the profile of the upper-triangle.
c                 we save the count of the number terms in the
c                 "assembled" profile - the count can become
c                 changed later in this routine if, for example, there
c                 are mpc equations which modify the starting
c                 "assembled" stiffness.
c
c                 we also get the number of non-zero terms on each row
c                 to right of diagonal ( start_kindex_locs )
c
      safety_factor = 10
      allocate( start_kindex_locs(neqns) )
      allocate( scol_flags(neqns,num_threads) )
      allocate( edest(mxedof,mxconn,num_threads) )
      edest = 0
      nrow_lists = mxconn * mxedof * safety_factor
      allocate( scol_lists(nrow_lists,num_threads) )
c
      call count_profile_symmetric( neqns, ncoeff, num_threads,
     &                              eqn_node_map, iprops, dof_eqn_map,
     &                              start_kindex_locs, edest,
     &                              scol_flags, scol_lists,
     &                              nrow_lists  )
c
      num_terms  = neqns + ncoeff
      old_ncoeff = ncoeff
      ncoeff_from_assembled_profile = ncoeff
      call thyme( 21, 2 )
      call thyme( 22, 1 )
c
      if( cpu_stats .and. show_details ) then
          write(out,9410) wcputime(1)
          write(out,9411) neqns
          write(out,9412) num_terms
      end if
c
c              4. build sparsity of equilibrium equations
c                 as defined by the k_indexes, k_ptrs
c                 vectors. see diagram at end for example.
c                 The work is done is sparse format so that
c                 we only need a workspace vector of neqns
c                 in length. Refer to comments inside called
c                 routine.
c
c                 start_kindex_locs from count.. routine has k_ptrs.
c                 convert to the starting locating in k_indexes for
c                 storage of column indexes for each equation.
c
      start_kindex_locs(neqns) = 0
c
      if( local_debug_2 )  then
          write(out,*) '...  start_kindex_locs ...'
          do i = 1, neqns
             write(out,*) '        ',i,start_kindex_locs(i)
          end do
      end if
c
      select case( neqns )
        case( 1 )
        case( 2 )
          start_kindex_locs(1) = 1
          if( ncoeff_from_assembled_profile == 0 ) ! 2x2 no (1,2) term
     &           start_kindex_locs(1) = 0
        case( 3: )
          start_srow = 0
          do srow = 1, neqns-1
            if( start_kindex_locs(srow) > 0 ) then
               start_srow = srow
               exit
            end if
          end do
          if( start_srow == 0 ) then
               write(out,9210)
               call die_abort
          end if
          count_previous = start_kindex_locs(start_srow)
          start_kindex_locs(start_srow) = 1
          do i = start_srow+1, neqns-1
           count_now = start_kindex_locs(i)
           start_kindex_locs(i) = start_kindex_locs(i-1) +
     &                            count_previous
           count_previous = count_now
          end do
          if( local_debug_2 ) then
              write(out,*) '... start_k_indexes_locs after update'
              do i = 1, neqns
                 write(out,*) '        ',i,start_kindex_locs(i)
              end do
          end if
      end select
c
      if( cpu_stats .and. show_details ) write(out,9419) wcputime(1)
      if( allocated( save_k_indexes ) ) deallocate( save_k_indexes )
      if( allocated( save_k_ptrs ) )    deallocate( save_k_ptrs )
      allocate( save_k_indexes(ncoeff_from_assembled_profile) )
      allocate( save_k_ptrs(neqns) )
      if( local_debug ) write(out,*)
     &                 ' @ 2.2 ... saving K-indexes, ptrs'
c
      call build_col_sparse_symm(
     &      neqns, num_threads, eqn_node_map, dof_eqn_map,
     &      save_k_indexes, save_k_ ptrs,
     &      start_kindex_locs, edest, scol_lists, nrow_lists, ncoeff )
c
      deallocate( start_kindex_locs, scol_flags, edest, scol_lists )

      if( cpu_stats .and. show_details ) write(out,9420) wcputime(1)
      call thyme( 22, 2 )
c
      ireturn = 3
      return
c
 9210 format('>> FATAL ERROR: computations of dense storage data',
     & /,    '                failed in equation solver.',
     & /,    '                inconsistencies in data structure',
     & /,    '                job terminated' )
 9300  format (
     &  10x, '>> solver wall time statistics (secs):'
     & /,15x,'number of MKL threads used      ',i10,
     & /,15x,'number of MPI ranks:            ',i10,
     & /,15x,'MKL version, build: ',5x,a8,1x,a8,
     & /,15x,'starting work                 @ ',f10.2 )
 9400  format (
     &  10x, '>> solver wall time statistics (secs):'
     & /,15x,'number of MKL threads used      ',i10,
     & /,15x,'MKL version, build: ',5x,a8,1x,a8,
     & /,15x,'starting work                 @ ',f10.2 )
 9409  format(
     &  15x, 'finished dof setup            @ ',f10.2 )
 9410  format(
     &  15x, 'non-zero term counts finished @ ',f10.2 )
 9411  format(
     &  15x, 'number of equations             ',i10)
 9412  format(
     &  15x, 'non-zero terms in profile       ',i10)
 9419  format(
     &  15x, 'building sparse ptrs, indexes @ ',f10.2 )
 9420  format(
     &  15x, 'k_ptrs, k_indexes done        @ ',f10.2 )
c
      end subroutine ds_setup_sparsity_local

c     ****************************************************************
c     *                                                              *
c     *     set up equation sparsity for distributed assembly        *
c     *                                                              *
c     *                       written by  : rhd                      *
c     *                   last modified : 04/13/2015                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine ds_setup_sparsity_distributed( ireturn )
      implicit none
c
      integer :: ireturn  ! local to this routine
c
c           If we have the same sparsity pattern we can skip
c           all of this and simply assemble the coefficients
c
      if( .not. allocated ( dof_eqn_map ) ) then
          allocate( dof_eqn_map(num_struct_dof) )
          allocate( eqn_node_map(num_struct_dof) )
      end if
      call dof_map( dof_eqn_map, cstmap, nonode, num_enode_dof,
     &              eqn_node_map, neqns )
c
      ireturn = 1
      if( neqns .le. 0 ) return
c
      if( neqns .ne. old_neqns ) new_size = .true.
      old_neqns = neqns
      deallocate(dof_eqn_map)
      deallocate(eqn_node_map)
c
      ireturn = 2
      if( .not. new_size ) return
c
c           Determine the sparsity structure of your blocks of elements
c
      call wmpi_alert_slaves( 41 )
      call determine_local_sparsity
c
c           Determine the initial map for assembling the sparse
c           structure this should just be a simple heuristic to keep
c           a good amount of locality
c
      call wmpi_alert_slaves( 42 )
      call determine_initial_map
c
c           Assemble the sparse structure for the above map, now we know
c           global nnz and full sparsity.  We can make a parmetis graph
c           to balance.
c
      call wmpi_alert_slaves( 43 )
      call assemble_sparsity
c
c           Determine a better mapping for load balance, pass rows as
c           necessary and make a note of the new equation numbers
c
      call wmpi_alert_slaves( 44 )
      call determine_ordering
c
c           Move the nnz structure around to match the final map
c
      call wmpi_alert_slaves( 45 )
      call move_sparsity
c
c           all done
c
      ireturn = 3
c
      return
c
      end subroutine ds_setup_sparsity_distributed

c     ****************************************************************
c     *                                                              *
c     *                        ds_root_assembly                      *
c     *                                                              *
c     *                       written by  : rhd                      *
c     *                   last modified : 09/14/2015                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine ds_root_assembly
      implicit none
c
      logical :: hypre_solver, rebuild_mpcs
      double precision :: px, py, pz
c
c              1. sparsity of equilibrium equations already defined on
c                 this or previous solution. neqns = 0 means all dof in
c                 model are constrained. just leave.
c
      if( local_debug ) write(out,*) ' @ 3'
      if( neqns .le. 0 ) return
c
c              2. make a working copy of the k_ptrs and k_indexes from
c                 saved values for use by the sparse solvers
c                 (some sprase solvers may alter them and/or
c                 the values/sizes are changed during insertion of
c                 mpcs). allocate space for
c                 the diagonal terms, the off-diagonal, non-zero
c                 coefficients and the effective load vector. initialize
c                 terms as needed.
c
c                 The dimension of the arrays are inflated in order to
c                 allow using different sparse solvers. The format of
c                 these vectors in the assembly process follow
c                 arrays in the now old old NASA-VSS sparse solver.
c
c                 Other solvers
c                 will have to perform mapping/transformation of these
c                 arrays to their native format, e.g if the diagonal
c                 terms are snot stored in a separate vector.
c                 Therefore, the sizes
c                 of these arrays are larger than the needed sizes
c                 in the NASA-VSS solver.
c
c                 NASA-VSS solver uses:
c
c                       k_ptrs(neqns)
c                       k_indexes(ncoeff)
c                       k_diag(neqns)
c                       k_coeffs(ncoeff)
c                       p_vec(neqns)
c                  u_vec(neqns)
c
c                 The hypre solver needs a non-symmetric CSR
c                 format matrix.  The hypre functions take care of
c                 conversion, but we need to overallocate the
c                 arrays in order to fit all the
c                 data.  Note for non-symmetric CSR we need:
c
c                       k_ptrs(neqns+1)
c                       k_indexes(2*ncoeff+neqns)
c                       k_coeffs(2*ncoeff+neqns)
c                       p_vec(neqns)
c                       u_vec(neqns)
c
c                 k_diag is not used in the CSR format, but we'll
c                 need it for assembly.
c
      hypre_solver =  solver_flag .eq. 9
      call thyme( 18, 1 )
      if( hypre_solver .or. asymmetric_assembly  ) then
            allocate( k_ptrs(neqns + 1 ),
     &                k_indexes(2*ncoeff + neqns),
     &                k_coeffs(2*ncoeff + neqns) )
      else
            allocate( k_ptrs(neqns + 1),
     &                k_indexes(ncoeff + neqns),
     &                k_coeffs(ncoeff + neqns) )
      end if
      allocate( k_diag(neqns), p_vec(neqns), u_vec(neqns) )
      if( local_debug ) write(out,*) '  @ 3.1'
      do i = 1, neqns
        k_ptrs(i) = save_k_ptrs(i)
        k_diag(i) = zero
        p_vec(i)  = zero
        u_vec(i)  = zero
      end do
      if( local_debug ) then
         write(out,*)  '  @ 3.2, ncoeff: ',ncoeff
         write(out,*)  '         ncoeff_from_assembled_profile: ',
     &                        ncoeff_from_assembled_profile
         write(out,*)  '         size k_indexes: ', size(k_indexes)
      end if
      k_indexes(1:ncoeff_from_assembled_profile) =
     &            save_k_indexes(1:ncoeff_from_assembled_profile)
      if( local_debug ) write(out,*) ' @ 4'
c
c              3. Here we branch for full (asymmetric) or
c                 upper-triangle symmetric assembly.
c
c                 the code in assem_by_row runs parallel with threads
c                 where a thread builds an entire equation row
c
      if( asymmetric_assembly ) then
c
c             3a. The only difference is we convert the UT
c                 VSS ptrs and indexes to a full, structurally
c                 symmetric CSR format and use a modified
c                 assem_by_row which picks up the LT terms.
c                 rows are assembled in parallel
c
        nnz = 0  ! total number of non-zero terms in assembled eqns.
        call convert_vss_csr( neqns, ncoeff, save_k_ptrs(1),
     &      save_k_indexes(1), nnz, k_ptrs(1), k_indexes(1) )
        if( cpu_stats .and. show_details ) write(out,9999) wcputime(1)
        k_coeffs = zero
        call assem_by_row_asymmetric( neqns, num_threads, eqn_node_map,
     &                                dof_eqn_map, k_ptrs, k_indexes,
     &                                iprops, k_coeffs )

      else
c             3b. assemble symmetric  equilibrium equations
c                 in sparse format using a row-by-row algorithm
c                 rather than a conventional element-by-element.
c                 rows are assembled in parallel.
c
        k_coeffs = zero
        call assem_by_row( neqns, num_threads, eqn_node_map,
     &                     dof_eqn_map, k_diag, k_coeffs,
     &                     k_indexes, k_ptrs, iprops,
     &                     dcp, noelem )
c
      end if ! for asymmetric/symmetric assembly
c
      if( cpu_stats .and. show_details ) write(out,9470) wcputime(1)
      if( local_debug ) write(out,*) ' @ 5'
c
c              4. set the rhs of the equations. p_vec has the
c                 neqns terms for unconstrained dof from res.
c
      if( local_debug2 ) then
         write(out,9500)
         i_offset = 0
         do i = 1, nonode
           px = res(i_offset+1)
           py = res(i_offset+2)
           pz = res(i_offset+3)
           i_offset = i_offset + 3
           write(out,9510) i, px, py, pz
         end do
      end if
c
      do i = 1, nodof
        eqn_num = dof_eqn_map(i) ! dof # -> eqn # map
        if( eqn_num .ne. 0 ) p_vec(eqn_num) = res(i)
      end do
c
      call thyme( 18, 2 )
      if( local_debug ) write(out,*) ' @ 6'
      if( cpu_stats .and. show_details ) write(out,9480) wcputime(1)
c
c              5. if any multi-point constraint equations exist,
c                 modify the stiffness matrix and rhs accordingly.
c                 if the # eqns has not changed, a much faster set
c                 of routines runs, taking advantage of info learned
c                 during the previous solve. note that the mpc
c                 insertion process for a new set of equations
c                 (new_size is .true.) resets the value of ncoeff !
c                 Thats why we save the "assembled" profile count
c                 in step 3 of the process above to construct the
c                 equation sparsity.
c
c                 We also need to recreate MPCs if they've
c                 been modified! Added a flag for this case
c
c                 ==>> MPCs/tied cons are supported only for symmetric
c                      assembly
c
c                 The re-insert code has a bug. minimize effects
c                 by rebuilding MPC data structures on 1st iteration
c                 of each step.
c
      if( local_debug2 ) then
          write(out,*) '... before MPC enforcement ...'
          write(out,*) '       neqns, ncoeff: ', neqns, ncoeff
          write(out,*) '       .... k_ptrs ....'
          write(out,9010) (i, k_ptrs(i), i = 1, neqns+1)
          k_ptrs(neqns+1) = 0
          write(out,*) '       .... k_indexes ....'
          write(out,*) '      sizeof k_indexes: ', sizeof(k_indexes)
          write(out,9010) (i, k_indexes(i), i = 1, ncoeff+neqns)
      end if


      if( tied_con_mpcs_constructed .or. mpcs_exist ) then
          rebuild_mpcs = now_iteration .eq. 1
          if( new_size .or. modified_mpcs .or. rebuild_mpcs ) then
            if( local_debug ) write(out,*) '.... @ 6.1'
            call mpc_insert_terms(neqns, k_ptrs, k_diag, dstmap,
     &                            dof_eqn_map)
            if( local_debug ) write(out,*) '.... @ 6.2'
            call mpc_modify_stiffness(neqns, k_diag, p_vec)
            if( local_debug ) write(out,*) '.... @ 6.3'
            call mpc_remove_dep_eqns(neqns, k_ptrs)
            if( local_debug ) write(out,*) '.... @ 6.4'
            modified_mpcs = .false.
         else
            if( local_debug ) write(out,*) '.... @ 6.5'
            call mpc_reinsert_terms(neqns, k_ptrs, dstmap, dof_eqn_map)
            if( local_debug ) write(out,*) '.... @ 6.6'
            call mpc_modify_stiffness(neqns, k_diag, p_vec)
            if( local_debug ) write(out,*) '.... @ 6.7'
            call mpc_remove_dep_eqns(neqns, k_ptrs)
            if( local_debug ) write(out,*) '.... @ 6.8'
         end if
         if( cpu_stats .and. show_details ) then
            write(out,9490) wcputime(1)
            write(out,9412) ncoeff+neqns
         end if
         if( local_debug ) then
            write(out,*) "... pvec after MPC mods ..."
            do k = 1, neqns
               write(out,9300) k, p_vec(k)
            end do
         end if

      end if ! on tied_con_mpcs_constructed .or. mpcs_exist
c
      if( local_debug2 ) then
          write(out,*) '... after MPC enforcement ...'
          write(out,*) '       neqns, ncoeff: ', neqns, ncoeff
          write(out,*) '       .... k_ptrs ....'
          write(out,9010) (i, k_ptrs(i), i = 1, neqns+1)
          write(out,*) '       .... k_indexes ....'
          write(out,9010) (i, k_indexes(i), i = 1, ncoeff+neqns)
      end if
c
      if( local_debug ) write(out,*) ' @ 7'

c
c              6. write equations ready to solve in various formats
c                 as requested by user. makes our equations available
c                 to drive solvers w/o dragging warp3d along
c
      if( sparse_stiff_output ) then
          call save_solver_sparse( neqns, ncoeff,
     &                           k_diag, p_vec, k_coeffs, k_ptrs,
     &                           k_indexes, dstmap, dof_eqn_map,
     &                           eqn_node_map, num_struct_dof,
     &                           nonode, sparse_stiff_output,
     &                           sparse_stiff_binary,
     &                           sparse_stiff_file_name, out )
          if( sparse_research )
     &         call store_solver( sparse_stiff_binary,
     &                            sparse_stiff_file_name,
     &                            neqns, num_struct_dof,
     &                            dof_eqn_map, eqn_node_map )
      end if
c
c              7. deallocate all element stiffness matrices.
c                 6/23/2017 - we now run w/o deleteing all those
c                 blocks on each global Newton iteration.
c
c      call estiff_allocate( 5 )
      if( local_debug ) write(*,*) ' @ 8'
c
      return
c
 9010 format(6x,12i8)
 9300 format(10x,i5,3(1x,f10.6))
 9412  format(
     &  15x, 'non-zero terms in profile       ',i10)
 9470  format(
     &  15x, 'sparse [k] assembly done      @ ',f10.2 )
 9480  format(
     &  15x, 'load vector assembly done     @ ',f10.2 )
 9490  format(
     &  15x, 'multi-point constraints added @ ',f10.2 )
 9500 format( 2x,"... Force vector for equation solving: ",
     &   /,   2x,"   node       /----   px, py, pz    ----/")
 9510 format( 2x,i7, 3e14.6 )
 9999 format(
     &  15x, 'upper triangle to full done   @ ',f10.2)
c
      end subroutine ds_root_assembly

c     ****************************************************************
c     *                                                              *
c     *                     ds_distributed_assembly                  *
c     *                                                              *
c     *                       written by  : rhd                      *
c     *                   last modified : 04/14/2015                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine ds_distributed_assembly
      implicit none
c
c
c           The idea here is to move the sparsity to the correct processors
c           and then assemble the actual data.
c           The first part (the sparsity) would be done above and would save
c           memory, but because we might spawn processors must be done here.
c           Move the nnz structure around to match the final map
c
      call wmpi_alert_slaves( 46 )
      call assemble_coefs( new_size )
c
c           and get and distribute the load vector
c
      call wmpi_alert_slaves( 47 )
      call assem_load_vec
      call wmpi_alert_slaves( 48 )
      call dist_final_setup
c
c           We need to allocate the solution vector on root
c
      neqns = curr_neqns()   !   a function in mpi code
      allocate( u_vec(neqns) )
c
c           It turns out that we need this
c
      if( .not. allocated ( dof_eqn_map ) ) then
          allocate( dof_eqn_map(num_struct_dof) )
          allocate( eqn_node_map(num_struct_dof) )
      end if
      call dof_map( dof_eqn_map, cstmap, nonode, num_enode_dof,
     &              eqn_node_map, neqns )
c
      return
c
c           Checksum matches.  We have the right distributed equations
c
      end subroutine ds_distributed_assembly


c     ****************************************************************
c     *                                                              *
c     *                     ds_drive_hypre                           *
c     *                                                              *
c     *                       written by  : rhd                      *
c     *                   last modified : 04/14/2015                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine ds_drive_hypre
      implicit none

       if( asymmetric_assembly ) then
          write(out,*) "Note to Mark.  This should work, but check"
        end if
c
        call wmpi_alert_slaves( 40 )
c
c           run the hyper solver on MPI ranks
c
        call iterative_sparse_hypre( u_vec, out, error_code )
c
c           The returned error_code means one of two things:
c
c           1. we've run into the condition "hypre can't solve the
c              increment before an adaptive reset is necessary", or
c
c           2. "these equations are too ill-conditioned for hypre
c              to solve."
c
c           If the adaptive Newton solution flag is not set
c           assume (2) and just exit.
c
c           Else hypre will keep count of how many adaptive
c           calls in a row it has triggered.  If this is the
c           first hypre-triggered adaptive step, then deallocate
c           the structs and return to mnalgr, requesting a usual Newton
c           adaptive sub-step.  If this is the second (-in a row-),
c           then quit out, as we're likely not going to recover.
c
      if( error_code .eq. 1 ) then
         if( adaptive_flag ) then
           if( precond_fail_count .ge. 2 ) then
              write(out,9600)
              write(out,9602)
              call die_gracefully
           else
              if( local_debug ) write (out,*) '... @ 8.1'
              write(out,9605)
              hyp_trigger_step = .true.
              if( .not. (parallel_assembly_used) ) then
                deallocate( k_diag )
                deallocate( k_coeffs )
                deallocate( k_ptrs )
                deallocate( k_indexes )
                deallocate( p_vec )
              end if
              deallocate( u_vec )
              return
           end if
         else ! user does not allow adaptive Newton.
               write(out,9607)
               write(out,9602)
               call die_gracefully
         end if
      end if
c
c           this means we solved the equations and do not need an
c           adaptive call, so set the flag accordingly.
c
      hyp_trigger_step = .false.
c
      return
c
 9600 format('>> ERROR: hypre has requested 2 adaptive steps',
     &          ' in a row due to non-convergence.')
 9602 format('>> Terminating analysis at this point',//)
 9605 format('>> Note: hypre has requested adaptive step reduction'/)
 9607 format('>> ERROR: hypre failed to converge. Try adjusting load',
     &  /    '          step size on hypre parameters.'//)
c
      end subroutine ds_drive_hypre

      end subroutine drive_assemble_solve




c ------------------------------------------------------------------------
c                        NASA - VSS
c
c               sparse matrix storage format
c
c   example:
c
c                1    2    3    4     5     6
c
c          1  | 100   1    2                5  |  | x1 |     | 201 |
c          2  |     200    6    7           9  |  | x2 |     | 202 |
c          3  |          300   10    11    12  |  | x3 |     | 203 |
c      a = 4  |                400   13    14  |  | x4 |  =  | 204 |
c          5  |                     500    15  |  | x5 |     | 205 |
c          6  |                           600  |  | x6 |     | 206 |
c
c     number of equations    = 6
c
c     number of coefficients = 12
c
c
c     k_ptrs    = { 3, 3, 3, 2, 1, 0}
c
c     k_indesxs = { 2, 3, 6,  3, 4, 6,  4, 5, 6,  5, 6,  6}
c
c     k_coefs   = { 1, 2, 5,  6, 7, 9, 10,11,12, 13,14, 15}
c
c     k_diag    = { 100, 200, 300, 400, 500, 600}
c
c     k_rhs     = { 201, 202, 203, 204, 205, 206}
c
c ------------------------------------------------------------------------
c     ****************************************************************
c     *                                                              *
c     *           subroutines save_solver_sparse                     *
c     *           write the equations in sparse format onto a        *
c     *           stream, unformatted file or a formatted file       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/15/2014                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine save_solver_sparse( neqns, ncoeff, k_diag, p_vec,
     &    k_coeffs, k_ptrs, k_indexes, dstmap, dof_eqn_map,
     &    eqn_node_map, num_struct_dof, nonode, sparse_stiff_output,
     &    sparse_stiff_binary, sparse_stiff_file_name, iout )
c
      implicit none
      integer  neqns, ncoeff, k_ptrs(*), k_indexes(*), dstmap(*),
     &         dof_eqn_map(*), eqn_node_map(*), num_struct_dof, nonode
      double precision
     &  k_diag(*), p_vec(*), k_coeffs(*)
      logical sparse_stiff_output, sparse_stiff_binary, connected
      character(len=*) :: sparse_stiff_file_name
      integer fileno, open_result, iout, i
c
      sparse_stiff_output = .false.
c
c                  find an available unit number to use
c
      do fileno = 11, 99
        inquire(unit=fileno, opened=  connected )
        if ( .not. connected ) go to 100
      end do
      write(iout,9000)
      return
c
 100  continue
      if ( sparse_stiff_binary ) then
        open(unit=fileno, file=sparse_stiff_file_name,
     &       status='replace', access='stream',form='unformatted',
     &       iostat=open_result  )
      else
        open(unit=fileno,file=sparse_stiff_file_name,
     &       status='replace', access='sequential', form='formatted',
     &       iostat=open_result  )
      end if
      if ( open_result .ne. 0 ) then
         write(iout,9010)
         return
      end if
c
      if ( sparse_stiff_binary ) then
          write(iout,9020)
          write(fileno) neqns, ncoeff
          write(iout,9022)
          write(fileno) (k_diag(i), i=1,neqns)
          write(iout,9024)
          write(fileno) (p_vec(i), i=1,neqns)
          write(iout,9026)
          write(fileno) (k_coeffs(i), i=1,ncoeff)
          write(iout,9028)
          write(fileno) (k_ptrs(i), i=1,neqns)
          write(iout,9030)
          write(fileno) (k_indexes(i), i=1,ncoeff)
          write(iout,9032)
          write(fileno) num_struct_dof, nonode
          write(fileno) dstmap(1:nonode)
          write(fileno) dof_eqn_map(1:num_struct_dof)
          write(fileno) eqn_node_map(1:neqns)
          close( fileno,status='keep' )
          write(iout,9034)
          return
      end if
c
      write(iout,9020)
      write(fileno,9500) neqns, ncoeff
      write(iout,9022)
      write(fileno,9510) (k_diag(i), i=1,neqns)
      write(iout,9024)
      write(fileno,9510) (p_vec(i), i=1,neqns)
      write(iout,9026)
      write(fileno,9510) (k_coeffs(i), i=1,ncoeff)
      write(iout,9028)
      write(fileno,9520) (k_ptrs(i), i=1,neqns)
      write(iout,9030)
      write(fileno,9520) (k_indexes(i), i=1,ncoeff)
      write(fileno,9500) num_struct_dof, nonode
      write(fileno,9520) dstmap(1:nonode)
      write(fileno,9520) dof_eqn_map(1:num_struct_dof)
      write(fileno,9520) eqn_node_map(1:num_struct_dof)
      write(iout,9032)
      close( fileno,status='keep' )
      write(iout,9034)
      return
c
 9000 format('>> WARNING: could not find a unit number to write',
     & /     '            sparse matrix. action skipped...',/)
 9010 format('>> WARNING: could not open file to save sparse [K]',
     & /     '            action skipped...',/)
 9020 format(15x,'writing sparse stiffness in binary format...')
 9022 format(17x,'> matrix sizes written...')
 9024 format(17x,'> diagonal terms written...')
 9026 format(17x,'> rhs written...')
 9028 format(17x,'> non-zero upper triangle terms written...')
 9030 format(17x,'> row pointers written...')
 9032 format(17x,'> non-zero column numbers written...')
 9034 format(17x,'> file closed...')
 9500 format(2i10)
 9510 format(4e20.12)
 9520 format(7i12)
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine store_solver                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/15/02 rhd                *
c     *                                                              *
c     *     write a binary file of key structure data after assembly *
c     *     and before solve by sparse solver. this file is used to  *
c     *     support external programs for solver research.           *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine store_solver( do_binary, file_name, neqns,
     &                         num_struct_dof, dof_eqn_map,
     &                         eqn_node_map )
      use global_data ! old common.main
      use main_data
      implicit integer (a-z)
      intrinsic size
      logical do_binary
      character(len=*) :: file_name
      dimension dof_eqn_map(*), eqn_node_map(*)
c
c               parameter and local variable declarations
c
      allocatable ele_info(:)
      real dumr, rzero
      double precision
     &     dumd
      character(len=1) :: dums
      data check_data_key, rzero / 2147483647, 0.0 /
c
      if ( .not. do_binary ) then
         write(out,*) '>>> WARNING: solver research data cannot'
         write(out,*) '             be written to an ascii file'
         write(out,*) '             request ignored....'
         return
      end if
c
      prec_fact = 2
      fileno = 11
        open(unit=fileno, file=file_name,
     &       status='old', access='sequential', form='unformatted',
     &       iostat=open_result, position='append'  )
      if ( open_result .ne. 0 ) then
         write(out,9070)
         return
      end if
      write(out,9050)
c
c                       write out sizes and integer scalars
c
      write(fileno) noelem, nonode, nodof, csthed, inctop, nelblk
      write(fileno) check_data_key
      write(out,9000)
c
c                       write out integer vectors
c
      call wrtbk( fileno, invdst, nodof )
      call wrtbk( fileno, incmap, noelem )
      call wrtbk( fileno, dstmap, nonode )
      call wrtbk( fileno, cstmap, nodof )
      call wrtbk( fileno, incid,  inctop )
      allocate( ele_info(noelem) )
      do i = 1, noelem
        ele_info(i) = iprops(2,i)
      end do
      call wrtbk( fileno, ele_info, noelem )
      deallocate( ele_info )
      write (fileno) check_data_key
c
c                       write more integer vectors, arrays
c
      call wrt2d( fileno, elblks(0,1) , 4, 4, nelblk  )
      call wrtbk( fileno, cp, mxedof )
      call wrtbk( fileno, dcp, mxedof )
      call wrt2d( fileno, icp, mxutsz, mxutsz, 2  )
      write (fileno) check_data_key
      write(out,9010)
c
c                       constraint data (double)
c
      call wrtbk( fileno, cnstrn, prec_fact*nodof )
      write (fileno) check_data_key
      write(out,9040)
c
c                       inverse incidences
c
      call store_solver2( fileno )
      write (fileno) check_data_key
c
c                       element stiffness blocks and element
c                       destination vectors
c
      call store_solver3( fileno )
      write(out,9080)
c
c                       number of actual equations to solve, the
c                       structure dof -> equation no. map vector, and
c                       equation no. -> structure node map vector
c
      write(fileno) neqns, num_struct_dof
      write(fileno) dof_eqn_map(1:num_struct_dof)
      write(fileno) eqn_node_map(1:num_struct_dof)
      write(fileno) check_data_key
c
c                       close the file
c
      close( fileno, status='keep' )
      write(out,9060)
c
c
c
 1000 format ( 3x, e16.6 )
 9000 format(17x,'> scalars written...')
 9010 format(17x,'> integer arrays written...')
 9040 format(17x,'> double precision arrays written...')
 9050 format(15x,'appending solver research data to binary file')
 9060 format(15x,'solver research file written and closed...')
 9070 format('>> WARNING: could not open file to save sparse [K]',
     & /     '            action skipped...',/)
 9080 format(17x,'> element stiffness and destinations written...')

      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine store_solver_2                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/15/02                   *
c     *                                                              *
c     *     write data into restart file for the requested data      *
c     *     structures which additional complexity                   *
c     *                                                              *
c     ****************************************************************
c
      subroutine store_solver2( fileno )
      use global_data ! old common.main
c
      use main_data, only : inverse_incidences, inverse_dof_map
c
      implicit integer (a-z)
      data check_data_key / 2147483647 /
c
      write(fileno) ( inverse_incidences(stnd)%element_count, stnd = 1,
     &                nonode )
      do stnd = 1, nonode
        ecount = inverse_incidences(stnd)%element_count
        write(fileno) ( inverse_incidences(stnd)%element_list(i),
     &                  i = 1, ecount )
      end do
      write(out,9000)
c
      do stnd = 1, nonode
        ecount = inverse_incidences(stnd)%element_count
        write(fileno) (( inverse_dof_map(stnd)%edof_table(i,j),
     &                  i = 1, ecount ), j = 1, 3 )
      end do
      write(out,9100)
c
      return
 9000 format(17x,'> inverse_incidences written...')
 9100 format(17x,'> inverse_dof_maps written...')
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine store_solver_3                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/15/02                   *
c     *                                                              *
c     *     write data into restart file for the requested data      *
c     *     structures which additional complexity                   *
c     *                                                              *
c     ****************************************************************
c
      subroutine store_solver3( fileno )
      use global_data ! old common.main
c
      use elem_block_data, only : edest_blocks, estiff_blocks
c
      implicit integer (a-z)
      data check_data_key / 2147483647 /
c
      do blk = 1, nelblk
        felem         = elblks(1,blk)
        span          = elblks(0,blk)
        nnode         = iprops(2,felem)
        num_enode_dof = iprops(4,felem)
        totdof        = nnode * num_enode_dof
        write(fileno) nnode, totdof, span
        write(fileno)  edest_blocks(blk)%ptr
        write(fileno)  estiff_blocks(blk)%ptr
        write(fileno) check_data_key
      end do
      return
      end


c     ****************************************************************
c     *                                                              *
c     *                 subroutine store_csr_formatted               *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 12/18/2016                 *
c     *                                                              *
c     *     write assembled equations IN CSR in the following format:*
c     *     n, nnz  (2I8)                                            *
c     *     value(i), col_index(i) (E25.16, 2X, I8)                  *
c     *     row_ptrs (I8)                                            *
c     *     RHS (E25.16)                                             *
c     *                                                              *
c     ****************************************************************
c
      subroutine store_csr_formatted( filename, n, nnz, row_ptrs,
     &                        col_indexes, values, rhs, diagonal )
      implicit none

      integer :: n, nnz, row_ptrs(*), col_indexes(*)
      double precision ::  diagonal(*), rhs(*), values(*)
      character(len=*) :: filename

      integer :: i, fileno
      logical :: connected
c
c          Find a file number
c
      do fileno = 11, 99
        inquire(unit=fileno, opened=connected )
        if ( .not. connected ) go to 100
      end do
      write(*,*) "Couldn't find a file number to write to, skipping"
      return
c
 100  continue
c
c           Basically the plan is to convert to CSR (with existing fn)
c           and write out in
c           the format specified above.
c
      call pardiso_symmetric_map( n, nnz, diagonal, rhs,values,
     &                   row_ptrs, col_indexes )
      nnz = nnz + n
      open (unit = fileno, file=filename)
      write(*,*) "Writing"
      write(fileno,5) n,nnz
5     format(2I10)
      do i=1, nnz
            write(fileno,10) values(i), col_indexes(i)
      end do
10    format(E25.16,2X,I10)
      do i = 1, n+1
            write(fileno,15) row_ptrs(i)
      end do
15    format(I10)
      do i = 1, n
            write(fileno,20) rhs(i)
      end do
20    format(E25.16)
      close(fileno)
      write (*,*) "Done writing"
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine convert_vss_csr                   *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 12/10/13                   *
c     *                                                              *
c     *     Convert the upper triangle vss sparsity pattern to       *
c     *     a full CSR format.                                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine convert_vss_csr( neqns, ncoeff, vss_ptrs,
     &                            vss_indexes, nnz, ptrs, indexes )
      implicit none

      integer, intent(in) :: neqns, ncoeff, vss_ptrs(neqns+1),
     &                        vss_indexes(neqns+ncoeff)
      integer, intent(out) :: nnz, ptrs(neqns+1),
     &                        indexes(2*ncoeff+neqns)
      integer :: i, j, s, o, k, c, r, work(neqns)
c
      nnz = neqns + 2*ncoeff
c
c           Copy the number of terms per row
c
      ptrs = 0
      ptrs(1:neqns) = vss_ptrs(1:neqns)
c
c           Add the diagonal terms
c
      ptrs(1:neqns) = ptrs(1:neqns) + 1
c
c           Add the number of terms per column, store for later use
c
      work = 0
      do i = 1, ncoeff
        ptrs(vss_indexes(i)) = ptrs(vss_indexes(i)) + 1
        work(vss_indexes(i)) = work(vss_indexes(i)) + 1
      end do
c
c           Convert to offsets
c
      do i = 2, neqns
        ptrs(i) = ptrs(i-1) + ptrs(i)
      end do
      do i = neqns+1, 2, -1
        ptrs(i) = ptrs(i-1) + 1
      end do
      ptrs(1) = 1
c
c           Insert the UT and diagonal indexes
c
      indexes = 0
      s = 1
      do i = 1, neqns
        o = ptrs(i) + work(i)
             ! o is pointing at the diagonal
        indexes(o) = i
        o = o + 1
        do j = s, s+vss_ptrs(i)-1
             ! j indexes the original structure
          indexes(o) = vss_indexes(j)
          o = o + 1
        end do
        s = j
      end do
c
c           Insert the LT indexes
c
      work = 0
      s = 1
      do c = 1, neqns
        do j = s, s+vss_ptrs(c) - 1
          r = vss_indexes(j)
          ! row r column c, columns will be strictly increasing
          o = ptrs(r) + work(r)
          indexes(o) = c
          work(r) = work(r) + 1
        end do
        s = j
      end do
c
c     All done algorithm is O(nnz), don't think we can do better
c
      return
      end
