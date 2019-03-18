c
c     ****************************************************************
c     *                                                              *
c     *  Driver to call Pardiso (threaded) asymmetric solver with    *
c     *  full CSR matrix. [K] is structurally symmetric but values   *
c     *  are not symmetric                                           *
c     *                                                              *
c     *  written by: mcm                                             *
c     *  last modified : 5/13/2017 rhd                               *
c     *                                                              *
c     ****************************************************************
c
      subroutine pardiso_unsymmetric( n, nnz, k_ptrs, k_indexes,
     &            k_coeffs, p_vec, u_vec, cpu_stats, itype, out,
     &            use_iterative )
      use performance_data, only : t_performance_start_pardiso,
     &                             t_performance_end_pardiso
c
      implicit none
c
c                parameter declarations
c
      integer :: n, nnz, k_ptrs(n+1), k_indexes(nnz), itype, out
      double precision :: k_coeffs(nnz), p_vec(n), u_vec(n)
      logical :: cpu_stats, use_iterative
c
c                local for pardiso
c
      integer :: maxfct, mnum, mtype, perm, nrhs, msglvl, error, phase,
     &           mkl_ooc_flag 
      integer, save :: pt(64), iparm(64), num_calls
      logical :: direct_solve
      logical, save :: pardiso_mat_defined
      data num_calls, pardiso_mat_defined / 0, .false. /
c
c                set driver and Pardiso parameters
c
      maxfct = 1
      mnum   = 1
      mtype  = 1 ! structurally symmetric, values not symmetric
      perm   = 0 ! will be dummy
      nrhs   = 1
      msglvl = 0 ! no printed stats from inside Pardiso
      mkl_ooc_flag = 0 ! not supported for unsymmetric
      direct_solve = .not. use_iterative
c
c                solution types (itype):
c                 1 - first time solution for a matrix:
c                               setup ordering method and perform
c                               pre-processing steps.
c                 2 - Solution of above same matrix equations
c                     with a new set of coefficients but same
c                     sparsity
c                 3 - no solution. just release data.
c
      call t_performance_start_pardiso
c
      select case( itype )
       case( 1 )
        call pardiso_unsymmetric_setup
        if( direct_solve )  call pardiso_unsymmetric_direct
        if( use_iterative ) call pardiso_unsymmetric_iterative
       case( 2 )
        if( direct_solve )  call pardiso_unsymmetric_direct
        if( use_iterative ) call pardiso_unsymmetric_iterative
       case(3 )
        call pardiso_unsymmetric_release
       case default
        call warp3d_pardiso_mess( 11, out, error, mkl_ooc_flag,
     &                            cpu_stats, iparm )
      end select
c
      call t_performance_end_pardiso
c
      return
c
      contains
c     ========
c
c     ******************************************************************
c     *       contains:   pardiso_unsymmetric_setup                    *
c     ******************************************************************
c
      subroutine pardiso_unsymmetric_setup
      implicit none
c
c
c              solve first time during execution or solve equations
c              with different sparisty
c
      call thyme(23, 1)
c
      iparm(1:64) = 0
      iparm(1) = 1 ! No solver default...
      iparm(2) = 3 ! Parallel fill-in reducing
      iparm(3) = 0 ! numbers of processors. MKL_NUM_THREADS overrides
      iparm(4) = 0 ! no iterative-direct algorithm
      if( use_iterative ) iparm(4) = 81 ! Unless we want it less strict
      iparm(5) = 0 ! No user ordering
      iparm(6) = 0 ! Write separate x
      iparm(7) = 0 ! not in use
      iparm(8) = 30 ! Do a lot of refinement
      iparm(10) = 13 ! Permute with 10E-13
      iparm(11) = 1 ! scaling.
      iparm(12) = 0 ! olve standard  equation: Ax = b
      iparm(13) = 1 ! Add some extra permutation for
c                     non-symmetric matrices
      iparm(14) = 0 ! Output: number of perturbed pivots
      iparm(15) = 0 ! not in use
      iparm(16) = 0 ! not in use
      iparm(17) = 0 ! not in use
      iparm(18) = -1 ! return: number of nonzeros in the factor LU
      iparm(19) = -1 ! return: Mflops for LU factorization
      iparm(20) = -1 ! return: Numbers of CG Iterations
      iparm(21) = 0 ! Different pivoting not available
      iparm(24) = 1 ! use 2-level parallelism for triangulation
      iparm(25) = 0 ! Parallel backsolve
      iparm(27) = 0 ! Don't check matrices
      iparm(28) = 0 ! double precision
      iparm(31) = 0 ! Full solve
      iparm(35) = 0 ! Fortran
      iparm(60) = 0 ! Do everything in core for now...
c
      error  = 0 ! make sure it gets set
      pt(1:64) = 0 ! init internal solver memory pointer only
c                    necessary for first call to pardiso.
c
c             if previously sovled a set of equations,release memory.
c
      if( pardiso_mat_defined ) then
        phase = -1
        call pardiso( pt, maxfct, mnum, mtype, phase, n, k_coeffs,
     &       k_ptrs, k_indexes, perm, nrhs, iparm, msglvl, p_vec,
     &       u_vec, error )
        pardiso_mat_defined = .false.
        call warp3d_pardiso_mess( 6, out, error, mkl_ooc_flag,
     &                            cpu_stats, iparm )
      end if
      call warp3d_pardiso_mess( 1, out, error, mkl_ooc_flag,
     &                         cpu_stats, iparm )
c
      phase = 11 ! reordering and symbolic factorization
      call pardiso( pt, maxfct, mnum, mtype, phase, n, k_coeffs,
     &     k_ptrs, k_indexes, perm, nrhs, iparm, msglvl, p_vec,
     &     u_vec, error )
      pardiso_mat_defined = .true.
      call warp3d_pardiso_mess( 2, out, error, mkl_ooc_flag,
     &                          cpu_stats, iparm )
c
      call thyme( 23, 2 )
      return
c
      end subroutine pardiso_unsymmetric_setup


c     ******************************************************************
c     *       contains:   pardiso_unsymmetric_release                  *
c     ******************************************************************
c
      subroutine pardiso_unsymmetric_release
      implicit none
c
      if( num_calls == 0 ) return
      if( .not. pardiso_mat_defined ) ! job aborted
     &      call warp3d_pardiso_mess( 7, out, error, mkl_ooc_flag,
     &                               cpu_stats, iparm)
      phase = -1
      call pardiso(pt, maxfct, mnum, mtype, phase, n, k_coeffs,
     &       k_ptrs, k_indexes, perm, nrhs, iparm, msglvl, p_vec,
     &       u_vec, error)
      pardiso_mat_defined = .false.
      num_calls = 0
      return
c
      end subroutine pardiso_unsymmetric_release

c     ******************************************************************
c     *       contains:   pardiso_unsymmetric_direct                   *
c     ******************************************************************
c
      subroutine pardiso_unsymmetric_direct
      implicit none
c
c              direct solve: factorization, forward/backward pass
c
      num_calls = num_calls + 1
      call thyme( 25, 1)
      phase = 22
      call warp3d_pardiso_mess( 4, out,  error, mkl_ooc_flag,
     &                          cpu_stats, iparm )
      call pardiso( pt, maxfct, mnum, mtype, phase, n, k_coeffs,
     &       k_ptrs, k_indexes, perm, nrhs, iparm, msglvl, p_vec,
     &       u_vec, error )
       call warp3d_pardiso_mess( 3, out, error, mkl_ooc_flag,
     &                           cpu_stats, iparm )
      call thyme( 25, 2 )
c
c             solve: forward and backward substitution
c
      phase = 33
      call thyme( 26, 1 )
      call pardiso(pt, maxfct, mnum, mtype, phase, n, k_coeffs,
     &     k_ptrs, k_indexes, perm, nrhs, iparm, msglvl, p_vec,
     &     u_vec, error)
      call warp3d_pardiso_mess( 5, out, error, mkl_ooc_flag,
     &                          cpu_stats, iparm )
      call thyme( 26, 2 )
c
      return
c
      end subroutine pardiso_unsymmetric_direct
c
c     ******************************************************************
c     *       contains:   pardiso_unsymmetric_iterative                *
c     ******************************************************************
c
      subroutine pardiso_unsymmetric_iterative
      implicit none
c
c              direct solve: factorization, forward/backward pass
c
c
      num_calls = num_calls + 1
      call thyme( 25, 1 )
      phase = 23 ! solve with iterative refinement
      call warp3d_pardiso_mess( 9, out, error, mkl_ooc_flag,
     &                          cpu_stats, iparm )
      call pardiso(pt, maxfct, mnum, mtype, phase, n, k_coeffs,
     &     k_ptrs, k_indexes, perm, nrhs, iparm, msglvl, p_vec,
     &     u_vec, error)
      call warp3d_pardiso_mess( 8, out, error, mkl_ooc_flag,
     &                          cpu_stats, iparm )
      call thyme( 25, 2 )
      return

      end subroutine pardiso_unsymmetric_iterative

      end subroutine pardiso_unsymmetric




