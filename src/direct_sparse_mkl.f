c
c     ****************************************************************
c     *                                                              *
c     *  Driver to call the MKL structurally symmetric solver with a *
c     *  full CSR matrix.  Very similar to direct_sparse_mkl but a   *
c     *  bit simpler because the matrix is already CSR.              *
c     *                                                              *
c     *  written by: mcm                                             *
c     *  last modified : 11/5/13                                     *
c     *                                                              *
c     ****************************************************************
c
      subroutine mkl_unsymmetric( n, nnz, k_ptrs, k_indexes, 
     &            k_coeffs, p_vec, u_vec, cpu_stats, itype, out, 
     &            use_iterative)
      use main_data, only: pardiso_first
      implicit none
c
      integer :: n, nnz, k_ptrs(n+1), k_indexes(nnz), itype, out
      double precision :: k_coeffs(nnz), p_vec(n), u_vec(n)
      logical :: cpu_stats, use_iterative
c
c     Local stuff for pardiso
c
      integer :: maxfct, mnum, mtype, perm, nrhs, iparm(64),
     &      msglvl, error, phase
      integer, save :: pt(64)
      logical :: mkl_ooc_flag
      integer, save :: num_calls
      logical, save :: pardiso_mat_defined
      data num_calls / 0 /
      data pardiso_mat_defined / .false. /
      data mkl_ooc_flag / .false. /
c
c     Set pardiso parameters
c
      maxfct = 1
      mnum = 1
      mtype = 1 ! Structurally symmetric, but not actually symmetric
      perm = 0 ! Will be dummy
      nrhs = 1
      if (pardiso_first) then
        pt = 0
        pardiso_first = .false.
      end if
c
      iparm = 0
      iparm(1) = 1 ! Non default...
      iparm(2) = 3 ! Parallel fill-in reducing
      iparm(4) = 0 ! no iterative-direct algorithm
      if( use_iterative ) iparm(4) = 81 ! Unless we want it
      iparm(5) = 0 ! No user ordering
      iparm(6) = 0 ! Write separate x
      iparm(8) = 30 ! Do a lot of refinement
      iparm(10) = 13 ! Permute with 10E-13
      iparm(11) = 1 ! Do matching
      iparm(12) = 0 ! Ax = b
      iparm(13) = 1 ! Add some extra permutation for non-symmetric matrices
      iparm(18) = -1 ! output nnz in factorization
      iparm(19) = -1 ! output MFLOPS
      iparm(21) = 0 ! Different pivoting not available
      iparm(24) = 1 ! More parallelism
      iparm(25) = 0 ! Parallel backsolve
      iparm(27) = 0 ! Don't check matrices
      iparm(28) = 0 ! double precision
      iparm(31) = 0 ! Full solve
      iparm(35) = 0 ! Fortran
      iparm(60) = 0 ! Do everything in core for now...
c
      msglvl = 0 ! Don't print stuff
c
      error = 0 ! Make sure it gets set
c
c     3 itype possibilities:
c           1 -- solve a system with a new sparsity
c           2 -- solve a system with a stored sparsity
c           3 -- release memory
c
      if (itype .eq. 3) then
        if (.not. pardiso_mat_defined) call warp3d_pardiso_mess(
     &            7, out, error, mkl_ooc_flag, cpu_stats, iparm)
        phase = -1
        call pardiso(pt, maxfct, mnum, mtype, phase, n, k_coeffs,
     &       k_ptrs, k_indexes, perm, nrhs, iparm, msglvl, p_vec,
     &       u_vec, error)
        pardiso_mat_defined = .false.
        num_calls = 0
        return
      end if

c     If we are in phase 2 check to make sure we've defined the matrix,
c     if we did skip ahead to numeric factorization and backsolve
      if ((itype .eq. 2) .and. pardiso_mat_defined) then
        goto 1000
      elseif (itype .eq. 2) then
        write (*,*) "Need the appropriate error message for no matrix"
        call die_gracefully
      end if
c
c     So factorize
c
      call thyme(23, 1)
c     Need to free memory, as we have a new matrix structure
      if (pardiso_mat_defined) then
        phase = -1
        call pardiso(pt, maxfct, mnum, mtype, phase, n, k_coeffs,
     &       k_ptrs, k_indexes, perm, nrhs, iparm, msglvl, p_vec,
     &       u_vec, error) 
        pardiso_mat_defined = .false.
        call warp3d_pardiso_mess( 6, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
      end if
c     No idea what this message says
      call warp3d_pardiso_mess( 1, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
c     Symbolic factor!
      phase = 11
      call pardiso(pt, maxfct, mnum, mtype, phase, n, k_coeffs,
     &     k_ptrs, k_indexes, perm, nrhs, iparm, msglvl, p_vec,
     &     u_vec, error)
      pardiso_mat_defined = .true.
      call warp3d_pardiso_mess( 2, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
      call thyme( 23, 2 )

c     At this point we have a factorized matrix
 1000 continue
c     Now lets factorize and solve.  If we're using the iterative version
c     this takes 1 call, if not it takes 2 calls
      num_calls = num_calls + 1
      call thyme( 25, 1)
      if (use_iterative) then
        phase = 23 ! Factorize and solve with iterative refinement
        call warp3d_pardiso_mess( 9, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
        call pardiso(pt, maxfct, mnum, mtype, phase, n, k_coeffs,
     &     k_ptrs, k_indexes, perm, nrhs, iparm, msglvl, p_vec,
     &     u_vec, error)
        call warp3d_pardiso_mess( 8, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
        call thyme( 25, 2 )
      else
        phase = 22 ! Factorize
        call warp3d_pardiso_mess( 4, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
        call pardiso(pt, maxfct, mnum, mtype, phase, n, k_coeffs,
     &       k_ptrs, k_indexes, perm, nrhs, iparm, msglvl, p_vec,
     &       u_vec, error)
        call warp3d_pardiso_mess( 3, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
        call thyme( 25, 2 )
        phase = 33 ! Actual forward/backsolve
        call thyme( 26 , 1)
        call pardiso(pt, maxfct, mnum, mtype, phase, n, k_coeffs,
     &     k_ptrs, k_indexes, perm, nrhs, iparm, msglvl, p_vec,
     &     u_vec, error)
        call warp3d_pardiso_mess( 5, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
        call thyme( 26, 2 )
      end if
c
c     Temp error printing
c
      if (error .ne. 0) then
            write (*,*) "ERROR"
            write (*,*) error
            call die_gracefully
      end if

      return

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *  Intel MKL Paradiso Sparse Solver Driver: A driver to call   *
c     *  the serial or parallel sparse (direct) symmetric solver     *
c     *  routines. The routines are part of the Intel MKL product    *
c     *                                                              *
c     *  written by: rhd   modified by: rhd                          *
c     *  last modified : 11/21/2015 rhd. hint about new_precond      *
c     *                                                              *
c     ****************************************************************
c
       subroutine direct_sparse_mkl( neq, ncoeff, k_diag, rhs,
     &                        sol_vec, eqn_coeffs, k_pointers,
     &                        k_indices, cpu_stats, itype, out, 
     &                        solver_out_of_core, solver_memory,
     &                        solver_scr_dir, solver_mkl_iterative,
     &                        suggested_new_precond )
c
      implicit none      
c
c                parameter declarations
c
      double precision ::  k_diag(*), rhs(*), sol_vec(*), 
     &                     eqn_coeffs(*)
      integer :: k_pointers(*), k_indices(*)
      logical :: cpu_stats, solver_out_of_core, solver_mkl_iterative,
     &           suggested_new_precond  
      integer :: itype, out, neq, ncoeff, solver_memory
      character (len=*) ::  solver_scr_dir
c
c                locally defined.
c
      real, external :: wcputime
      character(len=150) solver_directory 
      logical ::  pardiso_mat_defined, use_iterative
      integer ::  mkl_ooc_flag
      save num_calls, pardiso_mat_defined
c
c                local for pardiso
c
      external :: pardiso
      integer(kind=8), save :: pt(64)
      integer, save :: handle, iparm(64), msglvl, mtype
      integer :: maxfct, mnum, phase, nrhs, error,
     &           idum, num_calls 
      double precision :: ddum
c
      data  nrhs /1/, maxfct /1/, mnum /1/, num_calls / 0 /
      data  pardiso_mat_defined / .false. /
c
c                solution types (itype):
c                 1 - first time solution for a matrix:
c             		      setup ordering method and perform
c		                   pre-processing steps.
c                 2 - Solution of above same matrix equations
c                     with a new set of coefficients but same
c                     sparsity
c                 3 - no solution. just release data.
c
      use_iterative = solver_mkl_iterative
      if ( itype .eq. 3 .and. num_calls .ne. 0 ) then
        if( .not. pardiso_mat_defined ) 
     &   call warp3d_pardiso_mess( 7, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
        phase = -1 ! release internal memory
        call pardiso( pt, maxfct, mnum, mtype, phase, neq, ddum, 
     &                idum, idum, idum, nrhs, iparm, msglvl, 
     &                ddum, ddum, error)
        pardiso_mat_defined = .false.
        num_calls = 0
        return
      end if
c
c		            1.-  Map the input arrays from NASA-VSS
c                   format to the MKL format.
c
c                   if we're solving with the out-of-core option,
c                   the initialization, reordering steps must be
c                   repeated. for in-core solutions jump to
c                   factorization for itype = 2
c
      call map_vss_mkl( neq, ncoeff, k_diag, rhs, eqn_coeffs,
     &                  k_pointers, k_indices )
      if ( itype .eq. 2 .and. .not. solver_out_of_core ) go to 1000
c
c	            	2.-  Initialization. 
c
      call thyme( 23, 1 )
      if ( pardiso_mat_defined ) then
        phase = -1 ! release internal memory
        call pardiso( pt, maxfct, mnum, mtype, phase, neq, ddum, 
     &                idum, idum, idum, nrhs, iparm, msglvl, 
     &                ddum, ddum, error )
        pardiso_mat_defined = .false.
        call warp3d_pardiso_mess( 6, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
      end if
c
c              2a.-  create and write the out-of-core configuration
c                    file if necessary.
c
      mkl_ooc_flag = 0
      if( .not. use_iterative ) then
         if( solver_out_of_core ) then
            mkl_ooc_flag = 1
            call warp3d_dss_ooc( out, solver_memory, solver_scr_dir )
         end if 
      end if 
c  
      call warp3d_pardiso_mess( 1, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
c
c	             3.-  initialize, input  sparsity structure, reorder,
c                   symbolic factorization.
c
      iparm(1:64) = 0
      iparm(1) = 1 ! no solver default
      iparm(2) = 3 ! parallel reordering
      iparm(3) = 1 ! numbers of processors. MKL_NUM_THREADS overrides
      iparm(4) = 0 ! no iterative-direct algorithm
      if( use_iterative ) iparm(4) = 52
      iparm(5) = 0 ! no user fill-in reducing permutation
      iparm(6) = 0 ! =0 solution on the first n compoments of x
      iparm(7) = 0 ! not in use
      iparm(8) = 0 ! numbers of iterative refinement steps
      iparm(9) = 0 ! not in use
      iparm(10) = 13 ! perturb the pivot elements with 1E-13
      iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
      iparm(12) = 0 ! not in use
      iparm(13) = 0 ! maximum weighted matching algorithm is 
c                     switched-off (default for symmetric). 
c                     Try iparm(13) = 1 in case of inappropriate accuracy
      iparm(14) = 0 ! Output: number of perturbed pivots
      iparm(15) = 0 ! not in use
      iparm(16) = 0 ! not in use
      iparm(17) = 0 ! not in use
      iparm(18) = -1 ! Output: number of nonzeros in the factor LU
      iparm(19) = -1 ! Output: Mflops for LU factorization
      iparm(20) = -1 ! Output: Numbers of CG Iterations
      iparm(24) = 1 ! two-level factorization for better performance
      iparm(25) = 0 ! parallel forward-backward solve
      iparm(60) = mkl_ooc_flag 
      error = 0 ! initialize error flag
      msglvl = 0 ! print statistical information
      mtype = -2 ! symmetric, indefinite
      pt(1:64) = 0 ! init internal solver memory pointer only
c                    necessary for first call to pardiso.
c
      phase = 11 ! reordering and symbolic factorization
      call pardiso( pt, maxfct, mnum, mtype, phase, neq, eqn_coeffs,
     &              k_pointers, k_indices, idum, nrhs, iparm, msglvl,
     &              ddum, ddum, error )
      pardiso_mat_defined = .true.
      call warp3d_pardiso_mess( 2, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
      call thyme( 23, 2 )
c
c            		4.-  Numeric Factorization:
c                      -> input numerical values of coeffs
c                      -> Factor a matrix into L and L transpose
c                   For iterative, we do possible factoriztion and
c                   solve together in one call.
c
 1000 continue
      if( use_iterative ) then
        num_calls = num_calls + 1; call thyme( 25, 1 )
        phase = 23 ! iterative solve for displ
        iparm(4) = 52
        if( suggested_new_precond ) iparm(4) = 0 ! new factorization
        call warp3d_pardiso_mess( 9, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
        call pardiso( pt, maxfct, mnum, mtype, phase, neq,
     &     eqn_coeffs, k_pointers, k_indices, idum, nrhs, iparm, 
     &     msglvl, rhs, sol_vec, error )
        call warp3d_pardiso_mess( 8, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
        call thyme( 25, 2 )
        return
      end if
c
      num_calls = num_calls + 1
      call thyme( 25, 1 )
      phase = 22 ! only factorization
      call warp3d_pardiso_mess( 4, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
      call pardiso( pt, maxfct, mnum, mtype, phase, neq,
     &              eqn_coeffs, k_pointers, k_indices,
     &              idum, nrhs, iparm, msglvl, ddum, ddum, error )
      call warp3d_pardiso_mess( 3, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
      call thyme( 25, 2 )
c
c		            5.-  Solve:
c		                Forward and backward substitution
c
      call thyme( 26, 1 )
      iparm(8) = 0 ! max numbers of iterative refinement steps
      phase = 33   ! only factorization
      CALL pardiso( pt, maxfct, mnum, mtype, phase, neq,
     &              eqn_coeffs, k_pointers, k_indices, idum, nrhs, 
     &              iparm, msglvl, rhs, sol_vec, error )
      call warp3d_pardiso_mess( 5, out, 
     &        error, mkl_ooc_flag, cpu_stats, iparm )
      call thyme( 26, 2 )
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *  Messages during equation solving with the Intel MKL pkg.    *
c     *                                                              *
c     *  written by: rh   modified by: rhd  last modified: 7/27/10   *
c     *                                                              *
c     ****************************************************************
c
      subroutine warp3d_pardiso_mess( mess_no, iout, 
     &               ier, ooc_flag, cpu_stats, iparm )
      implicit none
c
      logical::  cpu_stats
      real, external :: wcputime
      integer :: iparm(*), mess_no, iout, ier, ooc_flag
c
      real :: start_factor_cpu_time
c
      select case ( mess_no )
c
        case( 1 )
         if( cpu_stats ) then 
            write(iout,9480) wcputime(1)
            if ( ooc_flag .ne. 0 ) write(iout,9580)
         end if
c
        case( 2 )
         if( ier .ne. 0 ) then
           write(iout,9470)
           write(iout,*) '       >> @ phase 11, ier: ',ier
           call die_gracefully
         end if
         if( cpu_stats )  then
           write(iout,9482) wcputime(1)
           write(iout,2010) dble(iparm(18))/1000.d0/1000.d0/1000.d0
           write(iout,2020) dble(iparm(19))/1000.d0
           write(iout,2022) dble(iparm(15))/1000.d0/1000.d0
         end if
c
        case( 3 )
         if( ier .ne. 0 ) then
            write(iout,9470)
            write(iout,*) '       >> @ phase 22, ier: ',ier
            call die_gracefully
         end if
         if( cpu_stats ) then
           write(iout,9490) wcputime(1)
           write(iout,2030) dble(iparm(17))/1000.d0/1000.d0
         end if

        case( 4 )
         if( cpu_stats ) then
           if( ooc_flag .eq. 0 ) write(iout,9184)
           if( ooc_flag .gt. 0 ) write(iout,9284)
           if( iparm(4) .gt. 0 ) write(iout,9286)
           start_factor_cpu_time = wcputime(1)
         end if
c
        case( 5 )
         if( ier .ne. 0 ) then
           write(iout,9470)
           write(iout,*) '       >> @ phase 33, ier: ',ier
           call die_gracefully
         end if
         if( cpu_stats ) write(iout,9492) wcputime(1)
c
        case( 6 )
         if( ier .ne. 0 ) then
           write(iout,9470)
           write(iout,*) '       >>  phase -1, ier: ',ier
           write(iout,*) '       >>  releasing solver data'
           call die_gracefully
         end if
c
        case( 7 )
           write(iout,9470)
           write(iout,*) '  >>  inconsistency @ 1'
           write(iout,*) '  >>  releasing solver data'
           call die_gracefully
c
        case( 8 )
         if( ier .ne. 0 ) then
           write(iout,9472)
           write(iout,*) '       >>  phase 52, ier: ',ier
           write(iout,*) '       >>  iterative solve'
           call die_gracefully
         end if
         if( cpu_stats ) then
           write(iout,9600) iparm(20)
         else 
           write(iout,9610) iparm(20)
         end if
c
        case( 9 )
          if( cpu_stats ) write(iout,9286)
c
      end select
      return
c
 2010 format(
     &  15x,'terms in factored matrix (B):    ', f9.2)
 2020 format(
     &  15x,'factorization op count (GFlop):  ', f9.2)
 2022 format(
     &  15x,'reordering memory (GB):          ', f9.2)
 2030 format(
     &  15x,'factorization memory (GB):       ', f9.2)
 9470  format(
     &  15x, 'FATAL ERRROR: mkl sparse solver' )
 9472  format(
     &  15x, 'FATAL ERRROR: mkl sparse -iterative- solver' )
 9480  format(
     &  15x, 'mkl solver initializing       @ ',f10.2 )
 9580  format(
     &  15x, ' -> out-of-memory operation' )
 9482  format(
     &  15x, 'reorder-symbolic factor. done @ ',f10.2 )
 9490  format(
     &  15x, 'numeric factorization done    @ ',f10.2 )
 9492  format(
     &  15x, 'numeric loadpass done         @ ',f10.2 )
 9184  format(
     &  15x, 'start in-memory factorization')
 9284  format(
     &  15x, 'start out-of-memory factorization')
 9286  format(
     &    15x, 'using conjugate gradient-krylov method',
     & /, 15x, 'with direct solver when needed and as the',
     & /, 15x, 'pre-conditioner: -- requires in-memory operation --' )
 9600  format(
     &  15x,'conjugate gradient-Krylov iters  ', i9)
 9610  format(7x,
     & '>> conjugate gradient-Krylov iterations:           ',i5)
c
      end

c
c     ****************************************************************
c     *                                                              *
c     *  map the default vss sparse solver format that warp3d        *
c     *  assembles into to the MKL sparse solver format              *
c     *                                                              *
c     *  written by: rh   modified by: rhd  last modified: 2/5/09    *
c     *                                                              *
c     ****************************************************************
c
      subroutine map_vss_mkl( neq, ncoeff, diag, rhs, amat, kpt, kind )
c
      implicit none
c
      double precision :: diag(*), rhs(*), amat(*)
      integer :: neq, ncoeff, kpt(*), kind(*)         
c
c	This sub maps arrays needed for the NASA sparse solver to
c	arrays in the format that the MKL 7 solver needs ( Harwell-Boeing
c	format)
c		See format example at the end of this subroutine
c
c	Definitions:
c		ncoeff: # of non-zero terms excluding the diagonal
c		amat  : has the dimension of (neq+ncoeff), however,
c			the input vector contains the vss format (
c			i.e. only the ncoeff number of terms); the output
c			will have all terms filled. Same with indices
c			array.
c
c   Note:    It is assumed that the last term in kpt(neq) == 0
c            is always zero.  Therefore, the second loop will not
c            be performed when i=neq ( See bellow )
c
c
c	Input(vss Format)   -------->     Output (MKL format)
c	     Terms Used                        Terms Used
c
c	       amat(ncoeff)                  k_coeff( ncoeff +neq )
c	       kind(ncoeff)                  k_indices ( ncoeff + neq )
c              kpt(neq)                      k_pointers ( neq + 1)
c
     	double precision ::  zero
	     data zero / 0.0 /
	     integer, allocatable :: kpt_l(:)
	     integer :: i_old, i_new, i, j, nonzt_sum
c
	     allocate( kpt_l( neq + 1) ); kpt_l(1:neq) = kpt(1:neq)
c
c                 Map amat and kind arrays
c      
	     i_old = ncoeff + 1; i_new = ncoeff + neq + 1
c
     	do i = neq, 1, -1
	      do j = 1, kpt(i)
	        i_new = i_new - 1; i_old = i_old - 1
	        amat(i_new) = amat(i_old); kind(i_new) = kind(i_old)
	      end do
	      if( diag(i) .ne. zero ) then
	        i_new = i_new - 1; amat(i_new) = diag(i); kind(i_new) = i
	      end if
    	end do
c
c                Map kpt ---> k_pointers
c
	     nonzt_sum = 1
     	kpt(1) = nonzt_sum
	     do i =1 ,neq
	      if( diag(i). ne . zero ) then
	        nonzt_sum = nonzt_sum + kpt_l(i) + 1; kpt(i+1) = nonzt_sum
	      end if
	     end do
c
c               Note: always pointers(neq+1) = # of non-zero terms  + 1
c
     	deallocate( kpt_l  )
     	return
      end

c     ****************************************************************
c     *                                                              *
c     *  Set-up configuration file for out-of-core solution          *
c     *                                                              *
c     *  written by: rh   modified by: rhd  last modified: 2/19/09   *
c     *                                                              *
c     ****************************************************************
c
      subroutine warp3d_dss_ooc( iout, solver_memory, solver_scr_dir )
      implicit none
c
      integer :: iout, solver_memory, ifileno, ii
      character(len=*) :: solver_scr_dir
      integer, external :: warp3d_get_device_number
c
c                 find an unused file number to use
c
      ifileno = warp3d_get_device_number()
      if( ifileno .eq. -1 ) then
         write(iout,*) '>>> Fatal error in warp3d_dss_ooc...'
         call die_gracefully
      end if
c
      open(unit=ifileno,dispose='keep',file='pardiso_ooc.cfg',
     &         status='unknown')
      ii = len_trim( solver_scr_dir )
      write(ifileno,9000) solver_scr_dir(1:ii)
      write(ifileno,*) 'MKL_PARDISO_OOC_MAX_CORE_SIZE = ',
     &             solver_memory
      write(ifileno,*) 'MKL_PARDISO_OOC_KEEP_FILE = 1'
      close(unit=ifileno)
c
 9000 format(1x, 'MKL_PARDISO_OOC_PATH = ',a) 
      return
      end
