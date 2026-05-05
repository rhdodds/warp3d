c     ****************************************************************
c     *                                                              *
c     *       drive MUMPS solver for symmetric equations             *
c     *                                                              *
c     *      last modified : 5/4/2026 rhd                            *
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
c     *            uses MKL LAPACK/BLAS for ifort/ifx (x86-64)       *
c     *            uses MKL LAPACK/BLAS for gfortran (x86-64)        *
c     *            uses Apple Accelerate w/ gfortran (arm64)         *
c     *                                                              *
c     *          last modified : 5/4/26  rhd                         *
c     *                                                              *
c     ****************************************************************
c
       subroutine drive_mumps( neq, ncoeff, diagonals, rhs,  
     &                        solution_vec, off_diagonals, row_counts,
     &                        col_indexes, print_time_stats,
     &                        itype, out )
c
      use dmumps_struc_def, only : DMUMPS_STRUC
      use global_data, only : mumps_solver_type
      use constants, only : zero
c
      implicit none
c
c              parameter declarations
c
      double precision :: diagonals(neq), rhs(neq), solution_vec(neq),
     &                    off_diagonals(ncoeff)
      integer :: row_counts(neq), col_indexes(ncoeff)
      logical :: print_time_stats
      integer :: itype, out, neq, ncoeff
c
c              locally defined.
c
      integer :: nnz, count, A_count, i9, negpiv, local_now_step,
     &           omp_threads, mkl_solver
      integer, save :: num_calls = 0
      integer, external :: omp_get_num_procs, mkl_get_max_threads,
     &                     omp_get_max_threads
      integer, save :: solver_method = 1
      logical :: prn, failed_factorization
      logical, parameter :: ldb = .false.
      double precision :: fact_entries, fact_gb, gflop_elim
c
      integer, target, save, allocatable :: w_irn(:), w_jrn(:)
      double precision, target, save, allocatable :: aval(:), lrhs(:)
      real :: sfactor_start, sfactor_end
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
c              if we enter with solver type = 1 (symmetric positive
c              definite) and factorization fails, switch to symmetric
c              indefinite and try again. Set global flag to
c              use type 2 from this point on. Also saved in restart.   
c              if solver type = 2 and factorization fails, stop.
c              if we have multi-point constraints, the assembler 
c              sets solver type = 2 (system init type = 1)
c
      prn =  print_time_stats
      nnz = ncoeff + neq
c      
      if( itype == 3 ) then
        deallocate( w_jrn, w_irn, aval, lrhs )
        id%JOB = -2
        call DMUMPS( id) 
        return
      end if  
c  
 100  continue     !  in case solver fails and change type
      if( ldb ) write(out,*) "....solution type: ", itype
c
c              Turn off mkl dynamic so it cannot change # threads
c
      omp_threads = omp_get_max_threads()
c
#ifdef MKL     
      mkl_solver = omp_threads
      call mkl_set_num_threads( mkl_solver )
      call mkl_set_dynamic( 0 )
      if( prn ) write(out,9404) omp_threads, mkl_solver
#else
      if( prn ) write(out,9405) omp_threads
#endif      
c
c              Check diagonals.
c              Build triplet storage format from WARP3D 
c              NASA-VSS upper triangle. Initialize MUMPS for first 
c              solve of these equations. reordering and
c              symbolic factorization.
c 
      call warp3d_mumps_check_diagonals( 0 )
c
      if( itype == 1 ) then
         call warp3d_mumps_vss_map()
         call warp3d_mumps_setup
      end if         
c
c              Load coefficients. factorize     
c       
      call warp3d_mumps_load_coeffs()
      id%A => aval
c      
      num_calls = num_calls + 1
      call thyme( 26, 1 )
      if( prn ) then
        if( mumps_solver_type == 1 ) write(out,9500)
        if( mumps_solver_type == 2 ) write(out,9505)
       end if 
      id%JOB = 2
      if( ldb ) sfactor_start =  wwalltime( 1 )  
      call DMUMPS( id )
      if( ldb ) then
         sfactor_end =  wwalltime( 1 )  
         write(out,*) "..... factor time: ", sfactor_end-sfactor_start
      end if  
c  
      failed_factorization = id%INFO(1) /= 0
      if( failed_factorization ) then
         if( mumps_solver_type == 2 ) then
           write(out,9510) 
           call die_abort
         end if 
         id%JOB = -2        !   reset MUMPS for indefinite solver
         call DMUMPS( id) 
         itype = 1          !   restart solver
         mumps_solver_type = 2
         write(out,9515)
         go to 100
      end if          
c      
      if( prn ) then
         write(out,9220) wwalltime( 1 )
         i9 = id%INFOG(9)! [L] factor storage 
         if( i9 .lt. 0 ) then
            fact_entries = dble(-i9) * 1.0d6
         else
            fact_entries = dble(i9)
         end if
         fact_gb = fact_entries * 8.0d0 / (1024.0d0**3)
         gflop_elim = id%RINFOG(3) / 1.0d9   ! after factorization: elimination flops
       end if  
       negpiv  = id%INFOG(12)   ! for SYM=1 or 2: total number of negative pivots
c
      if( prn ) then
         write(out,9240) fact_gb
         write(out,9250) gflop_elim
         write(out,9270) negpiv
      else
         write(out,9280) negpiv
      end if         
c
c              Forward/backward solve
c                  
      id%NRHS  = 1
      id%RHS   => lrhs ! loaded earlier
      id%JOB = 3
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
c
      return
      
 9200 format(15x,'map -> MUMPS completed        @ ',f10.2 )
 9220 format(15x,'factorization done            @ ',f10.2 )
 9230 format(15x,'forward/backward solve        @ ',f10.2 )
 9240 format(15x,'no. terms in [L] factor:',8x,f10.2,' x 10**9')
 9250 format(15x,'factorization GFLOP: :',10x,f10.2 )
 9260 format(15x,'peak memory usage (GB) :',8x,f10.2 )
 9270 format(15x,'negative pivots found :',9x,i10 )
 9280 format(15x,'negative pivots found :',i8 )
 9404 format(15x,'number of OMP threads used      ',i10,
     & /,15x,'number of MKL threads used      ',i10 )      
 9405 format(15x,'OMP threads for MUMPS           ',i10 ) 
 9500 format(15x,'use symm. pos definite factorization')
 9505 format(15x,'use symm. indefinite factorization')
 9510 format(/1x,'>>>>> Fatal error: the indefinite ',
     &           'factorization failed.',
     &  /20x,'solution terminated...',//)
 9515 format(/1x,'>>>>> warning: the symmetric positive definite ',
     &       'factorization failed.',
     &  /16x,'resetting the solver to try indefinite factorization',
     &  /16x,'for all further solves this analysis',// )
c
      return
c
      contains
c     ========
c     ******************************************************************
c     *       contains:   warp3d_mumps_check_diagonals                 *
c     ******************************************************************
c
      subroutine warp3d_mumps_check_diagonals( message_type )
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
      if( min_diagonal <= zero ) then
          write(out,9010)
          write(out,9012)
          mumps_solver_type = 2
      end if          
c
      return
c                
 9010 format(15x,'*** negative diagonal(s) found ***')
 9012 format(15x,'*** switching to indefinite solver ***')
 9020 format(15x,'min diagonal term: ',d10.3,
     &     /,15x,'max diagonal term: ',d10.3 )
c
      end subroutine warp3d_mumps_check_diagonals
c
c     ******************************************************************
c     *       contains:   warp3d_mumps_vss_map                         *
c     ******************************************************************
c
      subroutine warp3d_mumps_vss_map()
c
      implicit none
c
      integer :: p, eqn, num_off_diag, k, col_count      
c
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
c      
      return
      end subroutine warp3d_mumps_vss_map
c
c
c     ******************************************************************
c     *       contains:   mumps_load_coeffs                           *
c     ******************************************************************
c
      subroutine warp3d_mumps_load_coeffs()
c
      implicit none
      
      integer :: p, off_count, eqn, k, num_off_diag     
c
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
        lrhs(eqn) = rhs(eqn) ! make a target version
      end do
c
      return
      end subroutine warp3d_mumps_load_coeffs
c
c
c     ******************************************************************
c     *       contains:   mumps_setup                                  *
c     ******************************************************************
c
      subroutine warp3d_mumps_setup()
c
      implicit none
c 
c             initial call to MUMPS for new set of equations.
c            
      id%COMM = -987654       ! centralized, no MPI
      id%SYM  = mumps_solver_type !solver_method ! SPD matrix =1, indefinite = 2
      id%PAR  =  1            ! host participates
      id%JOB  = -1
      call DMUMPS(id)
      id%N     = neq
      id%NZ    = nnz
      id%IRN   => w_jrn  !irn
      id%JCN   => w_irn  !jcn
c
      id%cntl(1) = 0.01 
      id%ICNTL(1) = -1  ! output controls
      id%ICNTL(2) = -1  !
      id%ICNTL(3) = -1  !
      id%ICNTL(4) = 0   !
      id%ICNTL(7) = 0   ! force AMD
      if( neq > 1000 ) id%ICNTL(7) = 5   ! 4 = force PORD
c                                           5 = metis        
      id%ICNTL(14) = 0   ! increase working memory n%
      id%ICNTL(24) = 1   ! null pivot row detection
      if( prn .and. id%ICNTL(7)==0 ) write(out,9202)  
      if( prn .and. id%ICNTL(7)==4 ) write(out,9204)  
      if( prn .and. id%ICNTL(7)==5 ) write(out,9206)  
c
c            Phase 1: Analysis (ordering, symbolic factorization)
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
      return
c       
 9202 format(15x,'reordering by : AMD')
 9204 format(15x,'reordering by : PORD')
 9206 format(15x,'reordering by : METIS')
 9210 format(15x,'reorder + symbolic fact done  @ ',f10.2 )
c         
      end subroutine warp3d_mumps_setup      
c         
      end subroutine drive_mumps
