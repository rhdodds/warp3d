c     ****************************************************************
c     *                                                              *
c     *       drive MUMPS solver for symmetric equations             *
c     *                                                              *
c     *      last modified : 13/7/2026 rhd                           *
c     *                                                              *
c     ****************************************************************
c
c
       module dmumps_struc_def      
c      =======================
c
c  this file is part of mumps 5.8.2, released
c  on mon jan 12 15:17:08 utc 2026
c
c
c  copyright 1991-2026 cerfacs, cnrs, ens lyon, inp toulouse, inria,
c  mumps technologies, university of bordeaux.
c
c  this version of mumps is provided to you free of charge. it is
c  released under the cecill-c license
c  (see doc/cecill-c_v1-en.txt, doc/cecill-c_v1-fr.txt, and
c  https://cecill.info/licences/licence_cecill-c_v1-en.html)
c

!       include 'dmumps_struc.h'
!
!  this file is part of mumps 5.8.2, released
!  on mon jan 12 15:17:08 utc 2026
!
!
!  copyright 1991-2026 cerfacs, cnrs, ens lyon, inp toulouse, inria,
!  mumps technologies, university of bordeaux.
!
!  this version of mumps is provided to you free of charge. it is
!  released under the cecill-c license 
!  (see doc/cecill-c_v1-en.txt, doc/cecill-c_v1-fr.txt, and
!  https://cecill.info/licences/licence_cecill-c_v1-en.html)
!
      type dmumps_struc
        sequence
!
! this structure contains all parameters 
! for the interface to the user, plus internal
! information from the solver
!
! *****************
! input parameters
! *****************
!    -----------------
!    mpi communicator
!    -----------------
        integer :: comm
!    ------------------
!    problem definition
!    ------------------
!    solver (sym=0 unsymmetric,sym=1 symmetric positive definite, 
!        sym=2 general symmetric)
!    type of parallelism (par=1 host working, par=0 host not working)
        integer ::  sym, par
        integer ::  job 
!    --------------------
!    order of input matrix 
!    --------------------
        integer ::  n
!
!    ----------------------------------------
!    assembled input matrix : user interface
!    ----------------------------------------
        integer    :: nz  ! standard integer input + bwd. compat.
        integer(8) :: nnz ! 64-bit integer input
        double precision, dimension(:), pointer :: a
        integer, dimension(:), pointer :: irn, jcn
!    --------------
!    scaling arrays
!    --------------
        double precision, dimension(:), pointer :: colsca, rowsca
        double precision, dimension(:), pointer :: colsca_loc
        double precision, dimension(:), pointer :: rowsca_loc
        integer, dimension(:), pointer :: rowind, colind
        double precision, dimension(:), pointer :: pivots
!
!       ------------------------------------
!       case of distributed assembled matrix
!       matrix on entry:
!       ------------------------------------
        integer    :: nz_loc  ! standard integer input + bwd. compat.
        integer    :: pad1
        integer(8) :: nnz_loc ! 64-bit integer input
        integer, dimension(:), pointer :: irn_loc, jcn_loc
        double precision, dimension(:), pointer :: a_loc, pad2
!
!    ----------------------------------------
!    unassembled input matrix: user interface
!    ----------------------------------------
        integer :: nelt, pad3
        integer, dimension(:), pointer :: eltptr
        integer, dimension(:), pointer :: eltvar
        double precision, dimension(:), pointer :: a_elt, pad4
!
!    ---------------------------------------------
!    symmetric permutation : 
!               perm_in if given by user (optional)
!    ---------------------------------------------
        integer, dimension(:), pointer :: perm_in
!
!    ----------------
!    format by blocks
!    ----------------
        integer :: nblk, pad5
        integer, dimension(:), pointer :: blkptr
        integer, dimension(:), pointer :: blkvar
!
! ******************
! input/output data 
! ******************
!    --------------------------------------------------------
!    rhs / sol_loc
!    -------------
!       right-hand side and solution
!    -------------------------------------------------------
        double precision, dimension(:), pointer :: rhs, redrhs
        double precision, dimension(:), pointer :: rhs_sparse
        double precision, dimension(:), pointer :: sol_loc
        double precision, dimension(:), pointer :: rhs_loc
        integer, dimension(:), pointer :: irhs_sparse
        integer, dimension(:), pointer :: irhs_ptr
        integer, dimension(:), pointer :: isol_loc
        integer, dimension(:), pointer :: irhs_loc
        integer :: lrhs, nrhs, nz_rhs, nloc_rhs, lrhs_loc, lredrhs
        integer :: lsol_loc, nsol_loc
        integer :: ld_rhsintr, pad6
!    ----------------------------
!    control parameters,
!    statistics and output data
!    ---------------------------
        integer ::  icntl(60)
        integer ::  info(80) 
        integer :: infog(80)
        double precision ::  cost_subtrees
        double precision ::  cntl(15)
        double precision ::  rinfo(40)
        double precision ::  rinfog(40)
! the options array for metis/parmetis
        integer ::  metis_options(40)
!    ---------------------------------------------------------
!    permutations computed during analysis:
!       sym_perm: symmetric permutation 
!       uns_perm: column permutation (optional)
!    ---------------------------------------------------------
        integer, dimension(:), pointer :: sym_perm, uns_perm
! 
!    -----
!    schur
!    -----
        integer ::  nprow, npcol, mblock, nblock
        integer ::  schur_mloc, schur_nloc, schur_lld
        integer ::  size_schur
        double precision, dimension(:), pointer :: schur
        double precision, dimension(:), pointer :: schur_cinterface
        integer, dimension(:), pointer :: listvar_schur
!    -------------------------------------
!    case of distributed matrix on entry:
!    dmumps potentially provides mapping
!    -------------------------------------
        integer, dimension(:), pointer :: mapping
!    --------------
!    version number
!    --------------
        character(len=30) ::  version_number
!    -----------
!    out-of-core
!    -----------
        character(len=1023) :: ooc_tmpdir
        character(len=255) :: ooc_prefix
!    ------------------------------------------
!    name of file to dump a matrix/rhs to disk
!    ------------------------------------------
        character(len=1023) ::  write_problem
!    -----------
!    save/restore
!    -----------
        character(len=1023) :: save_dir
        character(len=255)  :: save_prefix
        character(len=7)   ::  pad7  
!
!
! **********************
! internal working data
! *********************
        integer(8) :: keep8(150), max_surf_master
        integer ::  inst_number
!       for mpi
        integer ::  comm_nodes, myid_nodes, comm_load
        integer ::  myid, nprocs, nslaves
        integer ::  ass_irecv
!       is is used for the factors + workspace for contrib. blocks
        integer, dimension(:), pointer :: is
        integer ::  keep(500)
!       the following data/arrays are computed during the analysis
!       phase and used during the factorization and solve phases.
        integer ::  lna
        integer ::  nbsa
        integer,pointer,dimension(:) :: step, ne_steps, nd_steps
        integer,pointer,dimension(:) :: frere_steps, dad_steps
        integer,pointer,dimension(:) :: fils, frtptr, frtelt
        integer(8),pointer,dimension(:) :: ptrar, ptr8arr
        integer,pointer,dimension(:) :: nincolarr,ninrowarr,ptrdebarr
        integer,pointer,dimension(:) :: na, procnode_steps
!       info for pruning tree 
        integer,pointer,dimension(:) :: step2node
!       ptlust_s and ptrfac are two pointer arrays computed during
!       factorization and used by the solve
        integer, dimension(:), pointer :: ptlust_s
        integer(8), dimension(:), pointer :: ptrfac
!       main real working arrays for factorization/solve phases
        double precision, dimension(:), pointer :: s
        real(kind(0.e0)), dimension(:), pointer :: lps
!       information on mapping
        integer, dimension(:), pointer :: procnode
!       input matrix ready for numerical assembly 
!           -arrowhead format in case of assembled matrix
!           -element format otherwise
!       element entry: internal data
        integer :: nelt_loc, leltvar
        integer, dimension(:), pointer :: eltproc
!       candidates and node partitionning
        integer, dimension(:,:), pointer :: candidates
        integer, dimension(:),   pointer :: istep_to_iniv2
        integer, dimension(:),   pointer :: future_niv2
        integer, dimension(:,:), pointer :: tab_pos_in_pere 
        logical, dimension(:),   pointer :: i_am_cand
!       for heterogeneous architecture
        integer, dimension(:), pointer :: mem_dist
!       compressed rhs
        integer, dimension(:),   pointer :: glob2loc_rhs
        logical  :: glob2loc_sol_alloc, pad11
        integer, dimension(:),   pointer :: glob2loc_sol
        double precision, dimension(:),   pointer :: rhsintr
!       info on the subtrees to be used during factorization
        double precision, dimension(:), pointer :: mem_subtree
        double precision, dimension(:), pointer :: cost_trav
        integer, dimension(:),   pointer :: my_root_sbtr
        integer, dimension(:),   pointer :: my_first_leaf
        integer, dimension(:),   pointer :: my_nb_leaf
        integer, dimension(:),   pointer :: depth_first
        integer, dimension(:),   pointer :: depth_first_seq
        integer, dimension(:),   pointer :: sbtr_id
        integer, dimension(:),   pointer :: sched_dep
        integer, dimension(:),   pointer :: sched_grp
        integer, dimension(:),   pointer :: sched_sbtr
        integer, dimension(:),   pointer :: croix_manu
        double precision, dimension(:),   pointer :: wk_user
        integer :: nbsa_local
        integer :: lwk_user
!    internal control array
        double precision ::  dkeep(230)
!    for simulating parallel out-of-core stack.
        double precision, dimension(:),pointer :: cb_son_size
!    instance number used/managed by the c/f77 interface
        integer ::  instance_number
!    ooc management data that must persist from factorization to solve.
        integer ::  ooc_max_nb_nodes_for_zone
        integer, dimension(:,:),   pointer :: ooc_inode_sequence
        integer(8),dimension(:,:), pointer :: ooc_size_of_block
        integer(8), dimension(:,:),   pointer :: ooc_vaddr
        integer,dimension(:), pointer :: ooc_total_nb_nodes
        integer,dimension(:), pointer :: ooc_nb_files
        integer :: ooc_nb_file_type,pad12
        integer,dimension(:), pointer :: ooc_file_name_length
        character,dimension(:,:), pointer :: ooc_file_names  
!    indices of nul pivots
        integer,dimension(:), pointer :: pivnul_list
!    array needed to manage additionnal candidate processor 
        integer, dimension(:,:), pointer :: sup_proc, pad14
!    lists of nodes where processors work. built/used in solve phase.
        integer, dimension(:), pointer :: iptr_working, working
!    internal data structures accessor
        character, dimension(:), pointer :: intr_encoding
!    low-rank
        integer, pointer, dimension(:) :: lrgroups
        integer :: nbgrp,pad13
!    pointer encoding for fdm_f data
        character, dimension(:), pointer :: fdm_f_encoding
!    pointer array encoding blr factors pointers
        character, dimension(:), pointer :: blrarray_encoding
!    multicore
        integer :: lpool_a_l0_omp, lpool_b_l0_omp
        integer :: l_phys_l0_omp
        integer :: l_virt_l0_omp
        integer :: ll0_omp_mapping, ll0_omp_factors
        integer(8) :: thread_la
! estimates before l0_omp
        integer, dimension(:,:), pointer    :: i4_l0_omp
        integer(8), dimension(:,:), pointer :: i8_l0_omp
! pool before l0_omp
        integer, dimension(:), pointer :: ipool_b_l0_omp
! pool after l0_omp
        integer, dimension(:), pointer :: ipool_a_l0_omp
! subtrees
        integer, dimension(:), pointer :: phys_l0_omp
! amalgamated subtrees
        integer, dimension(:), pointer :: virt_l0_omp
! mapping of amalgamated subtrees
        integer, dimension(:), pointer :: virt_l0_omp_mapping
! from heaviest to lowest subtree
        integer, dimension(:), pointer :: perm_l0_omp
! to get leafs in global pool
        integer, dimension(:), pointer :: ptr_leafs_l0_omp
! mapping of the subtree nodes
        integer, dimension(:), pointer :: l0_omp_mapping
! mpi to omp - mumps agile
        integer, dimension(:), pointer :: mtko_procs_map
! for rr on root
        double precision, dimension(:), pointer :: singular_values
        integer ::  nb_singular_values
        integer ::  deficiency, pad16
! to know if ooc files are associated to a saved and so if they should be removed.
        logical :: associated_ooc_files
      end type dmumps_struc
!
      end module dmumps_struc_def
c     ===========================
c
c     ****************************************************************
c     *                                                              *
c     *         drive MUMPS solver for symmetric equations           *
c     *            uses MKL LAPACK/BLAS for ifort/ifx                *
c     *            uses OpenBLAS for gfortran                        *
c     *                                                              *
c     *          last modified : 3/7/2026 rhd                        *
c     *                                                              *
c     ****************************************************************
c
       subroutine drive_mumps( neq, ncoeff, diagonals, rhs,  
     &                        solution_vec, off_diagonals, row_counts,
     &                        col_indexes, print_time_stats,
     &                        itype, out )
c
      use performance_data, only : t_performance_start_mumps,
     &                             t_performance_end_mumps
      use dmumps_struc_def
      use global_data, only : ltmstp
      use constants, only : zero
c
      implicit none
c
c              parameter declarations
c
      double precision ::  diagonals(neq), rhs(neq), solution_vec(neq),
     &                     off_diagonals(ncoeff)
      integer :: row_counts(neq), col_indexes(ncoeff)
      logical :: print_time_stats
      integer :: itype, out, neq, ncoeff
c
c              locally defined.
c
      integer :: nnz, count, k, col_count, eqn, p, num_off_diag, 
     &           off_count, A_count, i9, negpiv, local_now_step
      integer ::  n_total_threads, num_procs, omp_solver, mkl_solver, 
     &           save_omp_threads, save_mkl_threads
      integer, save :: num_calls = 0
      integer, external :: omp_get_num_procs, mkl_get_max_threads,
     &                     omp_get_max_threads
      logical :: prn 
      double precision :: fact_entries, fact_gb, gflop_elim
c
      integer, target, save, allocatable :: w_irn(:), w_jrn(:)
      double precision, target, save, allocatable :: aval(:), lrhs(:)
      real, external :: wwalltime
c
      type(DMUMPS_STRUC), save :: id
c
c              solution types (itype):
c                
c                 1 - first time solution for a matrix:
c                               setup ordering method and perform
c                               pre-processing steps.
c                 2 - Solution of above same matrix equations
c                     with a new set of coefficients but same
c                     sparsity
c                 3 - no solution. just release data.
c
      prn =  print_time_stats
      call t_performance_start_mumps
      nnz = ncoeff + neq
c      
      if( itype == 3 ) then
        deallocate( w_jrn, w_irn, aval, lrhs )
        id%JOB = -2
        call DMUMPS( id) 
        return
      end if  
c      
c              Save the number of omp and mkl threads 
c              to restore after mumps.
c              Use our allocator routine to set the number of omp and
c              mkl threads to use during MUMPS.
c
c              Turn off mkl dynamic so it cannot change # threads
c
      mkl_solver = 0
      save_omp_threads = omp_get_max_threads()
#ifdef MKL      
      save_mkl_threads = mkl_get_max_threads()
#endif    
c
      call get_mumps_thread_split( save_omp_threads, omp_solver,
     &                             mkl_solver )
c
#ifdef MKL      
      call mkl_set_num_threads( mkl_solver )
      call mkl_set_dynamic( 0 )
#endif      
      call omp_set_num_threads( omp_solver )
      if( prn ) write(out,9404) omp_solver, mkl_solver
c
c              Build triplet storage format from WARP3D 
c              NASA-VSS upper triangle      
c
      if( itype == 1 ) then
        if( allocated( w_jrn) ) deallocate( w_jrn )
        if( allocated( w_irn) ) deallocate( w_irn )
        if( allocated( aval) )  deallocate( aval )
        if( allocated( lrhs) )  deallocate( lrhs )
        allocate( w_jrn(nnz), w_irn(nnz), aval(nnz), lrhs(neq) )
        p         = 1
        col_count = 1
c      
        do eqn = 1, neq
         w_irn(p) = eqn
         w_jrn(p) = eqn
         p = p + 1
         num_off_diag = row_counts(eqn)
         do k = 1, num_off_diag
            w_irn(p) = eqn
            w_jrn(p) = col_indexes(col_count)
            p = p + 1
            col_count = col_count + 1
           end do
         end do
      end if  
c              
c              Load coefficients into aval vector
c           
      if( itype == 1 .or. itype == 2 ) then
          p = 1
          off_count = 1
c      
          do eqn = 1, neq
             aval(p)  = diagonals(eqn)
             p = p + 1
             num_off_diag = row_counts(eqn)
             do k = 1, num_off_diag
               aval(p)  = off_diagonals(off_count)
               p = p + 1
               off_count = off_count + 1
             end do
          end do
      end if         
c
      if( prn ) write(out,9200) wwalltime( 1 )
      call mumps_check_diagonals( 0 )
c
c              initialize MUMPS for first solve of these equations
c
      if( itype == 1 ) then
        id%COMM = -987654   ! centralized, no MPI
        id%SYM  = 2         ! SPD matrix =1, indefinite = 2
        id%PAR  = 1         ! host participates
        id%JOB  = -1
        call DMUMPS(id)
        id%N     = neq
        id%NZ    = nnz
        id%IRN   => w_jrn  !irn
        id%JCN   => w_irn  !jcn
        id%A     => aval
        id%NRHS  = 1
        id%RHS   => lrhs
c
        id%ICNTL(1) = -1  ! output controls
        id%ICNTL(2) = -1  !
        id%ICNTL(3) = -1  !
        id%ICNTL(4) = 0   !
        id%ICNTL(7) = 0   ! force AMD
        if( neq > 1000 ) id%ICNTL(7) = 4    ! force PORD
        id%ICNTL(14) = 0   ! increase working memory n%
        id%ICNTL(24) = 1   ! null pivot row detection
        if( prn .and. id%ICNTL(7)==0 ) write(out,9202)  
        if( prn .and. id%ICNTL(7)==4 ) write(out,9204)  
c
c              Phase 1: Analysis (ordering, symbolic factorization)
c
        call thyme( 23, 1 )
        id%JOB = 1
        call DMUMPS( id )
        if( id%INFO(1) /= 0 ) then
         write(out,*) 'MUMPS error in analysis: ', id%INFO(1)
         call die_abort
        end if
        write(out,9210) wwalltime( 1 )
        call thyme( 23, 2 )
c    
      end if  ! itype == 1
c
c              Numeric factorization & solve
c
      num_calls = num_calls + 1
      call thyme( 26, 1 )
      id%JOB = 2
      call DMUMPS( id )
      if( id%INFO(1) /= 0 ) then
         write(out,*) 'MUMPS error in factorization: ', id%INFO(1)
         call die_abort
      end if
      if( prn ) write(out,9220) wwalltime( 1 )
      i9 = id%INFOG(9)! [L] factor storage 
      if( i9 .lt. 0 ) then
         fact_entries = dble(-i9) * 1.0d6
      else
         fact_entries = dble(i9)
      end if
      fact_gb = fact_entries * 8.0d0 / (1024.0d0**3)
      gflop_elim = id%RINFOG(3) / 1.0d9   ! after factorization: elimination flops
      negpiv  = id%INFOG(12)   ! for SYM=1 or 2: total number of negative pivots
c
      if( prn ) write(out,9240) fact_gb
      if( prn ) write(out,9250) gflop_elim
      if( prn ) write(out,9270) negpiv
      if( .not. prn ) write(out,9280) negpiv
c
      id%JOB = 3
      lrhs(1:neq) = rhs(1:neq) ! make a target version
      call DMUMPS( id )
      if( id%INFO(1) /= 0 ) then
         write(out,*) 'MUMPS error in solve: ', id%INFO(1)
         call die_abort
      end if
      if( prn ) write(out,9230) wwalltime( 1 )
c
      solution_vec(1:neq) = lrhs(1:neq)
c
      if( prn ) write(out,9260) dble(id%INFOG(17))/1024.d0
c
      call thyme( 26, 2 )
      call t_performance_end_mumps
c
c              reset the number of omp, mkl threasds to values
c              before MUMPS solve
c      
#ifdef MKL
      call mkl_set_num_threads( save_mkl_threads )
#endif      
      call omp_set_num_threads( save_omp_threads )
c      
      return
      
 9200 format(15x,'VSS map -> MUMPS completed    @ ',f10.2 )
 9202 format(15x,'reordering by : AMD')
 9204 format(15x,'reordering by : PORD')
 9210 format(15x,'reorder + symbolic fact done  @ ',f10.2 )
 9220 format(15x,'factorization done            @ ',f10.2 )
 9230 format(15x,'forward/backward done         @ ',f10.2 )
 9240 format(15x,'no. terms in [L] factor:',8x,f10.2,' x 10**9')
 9250 format(15x,'factorization GFLOP: :',10x,f10.2 )
 9260 format(15x,'peak memory usage (GB) :',8x,f10.2 )
 9270 format(15x,'negative pivots found :',9x,i10 )
 9280 format(15x,'negative pivots found :',i8 )
 9404  format (15x,'number of OMP threads used      ',i10,
     & /,15x,'number of MKL threads used      ',i10 )      
c
      return
c
      contains
c     ========
c     ******************************************************************
c     *       contains:   mumps_check_diagonals                        *
c     ******************************************************************
c
      subroutine mumps_check_diagonals( message_type )
c
      implicit none
c
      integer :: i, message_type
      double precision :: min_diagonal, max_diagonal
c
      min_diagonal =  1.0d50
      max_diagonal = -1.0d50
      do i = 1, neq
       if( diagonals(i) .eq. zero ) cycle
       min_diagonal = min( min_diagonal, diagonals(i) )   
       max_diagonal = max( max_diagonal, diagonals(i) )   
      end do
c
      write(out,9020) min_diagonal, max_diagonal
      if( min_diagonal < zero ) write(out,9010)
c
      return
c                
 9010 format(15x,'*** negative diagonal(s) found ***')
 9020 format(15x,'min diagonal term: ',d10.3,
     &     /,15x,'max diagonal term: ',d10.3 )
c
      end subroutine mumps_check_diagonals      

      end subroutine drive_mumps
      
      subroutine get_mumps_thread_split( n_total,
     &                                   n_omp_solver,
     &                                   n_mkl_solver )

c     ------------------------------------------------------------
c     Determine thread split for MUMPS solver phase
c
c     Input:
c        n_total       total threads requested by user
c
c     Output:
c        n_omp_solver  OMP threads during solver
c        n_mkl_solver  MKL threads during solver
c
c     Conservative default rule based on benchmarking:
c
c        for n_total <= 5:
c           n_mkl = 1
c
c        for n_total  > 5:
c           n_mkl = nint( 0.30 * n_total )
c           n_mkl = max( n_mkl, 1 )
c           n_mkl = min( n_mkl, 4 )
c
c        n_omp = n_total - n_mkl
c
c     This favors keeping most threads in outer MUMPS/OpenMP
c     parallelism while giving MKL a modest number of threads
c     for dense BLAS/LAPACK work.
c     ------------------------------------------------------------

      implicit none

      integer          n_total
      integer          n_omp_solver
      integer          n_mkl_solver

      integer          mkl_cap
      integer          small_cut
      double precision frac

      parameter      ( mkl_cap   = 4 )
      parameter      ( small_cut = 5 )
      parameter      ( frac      = 0.30d0 )

c     ---- degenerate/serial case

      if( n_total .le. 1 ) then
         n_omp_solver = 1
         n_mkl_solver = 1
         return
      endif

c     ---- small thread counts: keep MKL serial

      if( n_total .le. small_cut ) then
         n_mkl_solver = 1
         n_omp_solver = n_total - n_mkl_solver
         if( n_omp_solver .lt. 1 ) n_omp_solver = 1
         return
      endif

c     ---- general conservative rule

      n_mkl_solver = nint( frac * dble(n_total) )

      if( n_mkl_solver .lt. 1 ) n_mkl_solver = 1
      if( n_mkl_solver .gt. mkl_cap ) n_mkl_solver = mkl_cap

      n_omp_solver = n_total - n_mkl_solver

c     ---- guards

      if( n_omp_solver .lt. 1 ) n_omp_solver = 1

      return
      end      