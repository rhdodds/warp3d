c     ****************************************************************
c     *                                                              *
c     *  drive Pardiso solver for symmetric equations: direct or     *
c     *  iterative. Pardiso is threads only. CPardiso is MPI +       *
c     *  threads                                                     *
c     *                                                              *
c     *      last modified : 11/18/2019 rhd                          *
c     *                                                              *
c     ****************************************************************
c
       subroutine pardiso_symmetric( neq, ncoeff, k_diag, rhs,
     &                        sol_vec, eqn_coeffs, k_pointers,
     &                        k_indices, print_cpu_stats, itype, out,
     &                        solver_out_of_core, solver_memory,
     &                        solver_scr_dir, solver_mkl_iterative )
c
      use performance_data, only : t_performance_start_pardiso,
     &                             t_performance_end_pardiso
      use  global_data, only : ltmstp, solver_threads, num_threads
c
      implicit none
c
c                parameter declarations
c
      double precision ::  k_diag(*), rhs(*), sol_vec(*),
     &                     eqn_coeffs(*)
      integer :: k_pointers(*), k_indices(*)
      logical :: print_cpu_stats, solver_out_of_core,
     &           solver_mkl_iterative
      integer :: itype, out, neq, ncoeff, solver_memory
      character (len=*) ::  solver_scr_dir
c
c                locally defined.
c
      real, external :: wcputime
      logical ::  pardiso_mat_defined, use_iterative, direct_solve
      logical :: use_non_pardiso
      integer ::  mkl_ooc_flag, local_now_step
      double precision, parameter :: zero = 0.d0
      double precision :: min_diagonal, max_diagonal
      save num_calls, pardiso_mat_defined
c
c                local for pardiso
c
      integer(kind=8), save :: pt(64)
      integer, save :: iparm(64), msglvl, mtype
      integer :: maxfct, mnum, phase, nrhs, error, idum, num_calls
      double precision :: ddum
c
      data  nrhs /1/, maxfct /1/, mnum /1/, num_calls / 0 /,
     &      pardiso_mat_defined / .false. /
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
      use_iterative = solver_mkl_iterative
      direct_solve = .not. use_iterative
c
      local_now_step = ltmstp
      use_non_pardiso = local_now_step > 100000000 ! 181  ! step number
      if( use_non_pardiso ) then
        call pardiso_dsptrf
        return
      end if
c
c                  set number of threads Pardiso must use. Turn off
c                  dynamic adjust of thread count.
c
      call mkl_set_num_threads( solver_threads )
      call omp_set_num_threads( solver_threads )
      call mkl_set_dynamic( 1 )
c
      select case( itype )
      case( 1 )
        call pardiso_symmetric_setup
        call pardiso_check_diagonals( 0 ) 
        if( use_iterative ) call pardiso_symmetric_iterative
        if( direct_solve )  call pardiso_symmetric_direct
      case( 2 )
        call pardiso_check_diagonals( 1 )
        call pardiso_symmetric_map( neq, ncoeff, k_diag,
     &                              eqn_coeffs, k_pointers, k_indices )
        if( use_iterative ) call pardiso_symmetric_iterative
        if( direct_solve )  call pardiso_symmetric_direct
      case( 3 )
        call pardiso_symmetric_release
      case default ! then die
        call warp3d_pardiso_mess( 11, out, error, mkl_ooc_flag,
     &                            print_cpu_stats, iparm )
      end select
c      
      call omp_set_num_threads( num_threads )
c
      call t_performance_end_pardiso
      return
c
      contains
c     ========
c     ******************************************************************
c     *       contains:   pardiso_check_diagonals                      *
c     ******************************************************************
c
      subroutine pardiso_check_diagonals( message_type )
      implicit none
c
      integer :: i, message_type
c
      min_diagonal =  1.0d50
      max_diagonal = -1.0d50
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
      do i = 1, neq
       if( k_diag(i) .eq. zero ) cycle
       min_diagonal = min( min_diagonal, k_diag(i) )   
       max_diagonal = max( max_diagonal, k_diag(i) )   
      end do
c
      if( message_type == 1 ) write(out,9000) min_diagonal, max_diagonal
      if( message_type == 0 ) write(out,9020) min_diagonal, max_diagonal
      if( min_diagonal < zero ) write(out,9010)
c
      return
c                
 9000 format(12x,'min, max diagonal terms: ',2d10.3)
 9010 format(12x,'*** negative diagonal(s) found ***')
 9020 format(17x,'min diagonal term: ',d10.3,
     &     /,17x,'max diagonal term: ',d10.3 )
c
      end subroutine pardiso_check_diagonals      
      
c     ******************************************************************
c     *       contains:   pardiso_dsptrf                               *
c     ******************************************************************
c
      subroutine pardiso_dsptrf
c
c          convert to full matrix form and use the dsptrf 
c          solver for symmetric indefinite systems.
c          runs parallel but takes much longer than
c          pardiso
c          
c
      implicit none
c
      integer :: i, j, row, col, ap_position, next, nsinfo,
     &           count, nsymm, num_indexes, nterms_prior_row,
     &           neqns, nterms_row
c
      integer, allocatable, dimension(:) ::
     &   ipiv(:),  row_pos(:)
      double precision, allocatable, dimension(:) ::
     &              ap(:), solution_vec(:)
      double precision, parameter :: zero = 0.0d0
      real :: time_result, solver_time_start, solver_time
      real, external :: wwalltime
      logical, parameter :: ldebug = .false.

      call pardiso_symmetric_release
c
      neqns = neq
      nsymm = neqns*(neqns+1)/2
      num_indexes = nsymm - neqns 
      write(out,*) '... neqns, num_indexes: ', neqns, num_indexes
      allocate( ipiv(10*neqns), ap(nsymm), row_pos(neqns) )
      ap = zero
c
c          insert diagonals. build locations of diagonal terms in ap
c 
      row_pos(1) = 1
      nterms_prior_row = neqns
      do i = 2, neqns
        row_pos(i) = row_pos(i-1) + nterms_prior_row
        nterms_prior_row = nterms_prior_row - 1
      end do
      if( row_pos(neqns) .ne. nsymm ) stop ' ..bad last row_pos'
c
      do row = 1, neqns
       if( ldebug ) then
         write(out,*) ' row, row_pos, k_diag: ', row, row_pos(row),
     &            k_diag(row)
       end if
       ap(row_pos(row)) = k_diag(row)
      end do
c
      count = 0
      do i = 1, nsymm
        if( ap(i) /= zero ) count = count + 1
      end do
      write(out,*) '... ap count non-zero terms: ',count
c
c          insert off diagonals
c  
      next = 1
      do row = 1, neqns - 1
       nterms_row = k_pointers(row)
       if( ldebug ) then 
        write(out,*) '>> start row, nrterms_row, row_position: ', 
     &         row, nterms_row, row_pos(row)
       end if
       if( nterms_row == 0 ) then
          write(out,*) '... row no off diagonal terms: ', row
          cycle 
       end if
       do j = 1, nterms_row
        col = k_indices(next)
        ap_position = row_pos(row) + ( col - row )
        ap(ap_position) = eqn_coeffs(next)
        if( ldebug ) then
         write(out,9010) j, col, ap_position, eqn_coeffs(next)
        end if  
        next = next + 1
       end do
      end do
c
      count = 0
      do i = 1, nsymm
        if( ap(i) /= zero ) count = count + 1
      end do
      write(out,*) '... ap count non-zero terms: ',count

      nsinfo = 0

      write(out,*) '... call dsptrf ...'
      solver_time_start = wwalltime(1)
      call dsptrf( 'L', neqns, ap, ipiv, nsinfo )
      write(out,*) '.. nsinfo after factor: ', nsinfo
      sol_vec(1:neqns) = rhs(1:neqns)
      call dsptrs( 'L', neqns, 1, ap, ipiv, sol_vec, neqns, nsinfo )
      write(out,*) '.. nsinfo after solve: ', nsinfo
      solver_time = wwalltime( 1 ) - solver_time_start
      write(out,*) '... wall time for solution only: ',solver_time
c
      deallocate( row_pos )
      deallocate( ap )
      deallocate( ipiv )

      return
 9010 format(5x,'.. j, col, ap_position, k_coefs: ',3i5,f14.4)


      end subroutine pardiso_dsptrf







c     ******************************************************************
c     *       contains:   pardiso_symmetric_release                    *
c     ******************************************************************
c
      subroutine pardiso_symmetric_release
      implicit none
c
      if ( num_calls .eq. 0 ) return
      if( .not. pardiso_mat_defined ) ! job will be aborted
     &  call warp3d_pardiso_mess( 7, out,
     &        error, mkl_ooc_flag, print_cpu_stats, iparm )
        phase = -1 ! release internal memory
        call pardiso( pt, maxfct, mnum, mtype, phase, neq, ddum,
     &                idum, idum, idum, nrhs, iparm, msglvl,
     &                ddum, ddum, error)
        pardiso_mat_defined = .false.
        num_calls = 0
      return
c
      end subroutine pardiso_symmetric_release
c
c     ******************************************************************
c     *       contains:   pardiso_symmetric_setup                      *
c     ******************************************************************
c
      subroutine pardiso_symmetric_setup
      implicit none
c
c              map vss to csr equation storage
c
      call pardiso_symmetric_map( neq, ncoeff, k_diag, eqn_coeffs,
     &                            k_pointers, k_indices )
c
c             if previously sovled a set of equations,release memory.
c
      call thyme( 23, 1 )
      if ( pardiso_mat_defined ) then
        phase = -1 ! release internal memory
        call pardiso( pt, maxfct, mnum, mtype, phase, neq, ddum, idum,
     &                idum, idum, nrhs, iparm, msglvl, ddum, ddum,
     &                error )
        pardiso_mat_defined = .false.
        call warp3d_pardiso_mess( 6, out, error, mkl_ooc_flag,
     &                           print_cpu_stats, iparm )
      end if
c
c              create and write the out-of-core configuration
c              file if necessary.
c
      mkl_ooc_flag = 0
      if( direct_solve .and. solver_out_of_core ) then
        mkl_ooc_flag = 1
        call warp3d_dss_ooc( out, solver_memory, solver_scr_dir )
      end if
c
      call warp3d_pardiso_mess( 1, out, error, mkl_ooc_flag,
     &                          print_cpu_stats, iparm )
c
c              initialize, input  sparsity structure, reorder,
c              symbolic factorization.
c
      iparm(1:64) = 0
      iparm(1) = 1 ! no solver default
      iparm(2) = 3 ! 3 ! parallel reordering
      iparm(3) = 0 ! numbers of processors. MKL_NUM_THREADS overrides
      iparm(4) = 0 ! no iterative-direct algorithm
      if( use_iterative ) iparm(4) = 52
      iparm(5) = 0 ! no user fill-in reducing permutation
      iparm(6) = 0 ! =0 solution on the first n compoments of x
      iparm(7) = 0 ! not in use
      iparm(8) = 0 ! numbers of iterative refinement steps
      iparm(9) = 0 ! not in use
      iparm(10) = 8 ! perturb the pivot elements with 1E-13
      iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
      iparm(12) = 0 ! not in use
      iparm(13) = 1 ! maximum weighted matching algorithm is
c                     switched-off (default for symmetric).
c                     Try iparm(13) = 1 in case of inappropriate accuracy
      iparm(14) = 0 ! Output: number of perturbed pivots
      iparm(15) = 0 ! not in use
      iparm(16) = 0 ! not in use
      iparm(17) = 0 ! not in use
      iparm(18) = -1 ! return: number of nonzeros in the factor LU
      iparm(19) = -1 ! return: Mflops for LU factorization
      iparm(20) = -1 ! return: Numbers of CG Iterations
      iparm(21) = 0  ! Bunch-Kaufman pivoting symmetric indefinite
      iparm(23) = 0  ! will have # negative eigenvalues
      iparm(24) = 1 ! use 2 level factorization
      iparm(25) = 2 ! parallel forward-backward solve
      iparm(27) = 0 !  check input matrix for errors (=1)
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
      call warp3d_pardiso_mess( 2, out, error, mkl_ooc_flag,
     &                          print_cpu_stats, iparm )
      call thyme( 23, 2 )
c
      return
c
      end subroutine pardiso_symmetric_setup
c
c     ******************************************************************
c     *       contains:   pardiso_symmetric_iterative                  *
c     ******************************************************************
c
      subroutine pardiso_symmetric_iterative
      implicit none
c
c             iterative solution- we do possible factoriztion and
c             solve together in one call.
c
c             for now, ignore the provided hint about when
c             are good times to refactor
c
      if( solver_out_of_core ) then ! invalid condition
        call warp3d_pardiso_mess( 10, out,  error, mkl_ooc_flag,
     &                            print_cpu_stats, iparm )
      end if
c
      num_calls = num_calls + 1
      call thyme( 25, 1 )
      phase     = 23 ! iterative solve for displ
      iparm(4)  = 52
      call warp3d_pardiso_mess( 9, out,  error, mkl_ooc_flag,
     &                         print_cpu_stats, iparm )
      call pardiso( pt, maxfct, mnum, mtype, phase, neq,
     &   eqn_coeffs, k_pointers, k_indices, idum, nrhs, iparm,
     &   msglvl, rhs, sol_vec, error )
      call warp3d_pardiso_mess( 8, out, error, mkl_ooc_flag,
     &                          print_cpu_stats, iparm )
      call thyme( 25, 2 )
      return
c
      end subroutine pardiso_symmetric_iterative
c
c     ******************************************************************
c     *       contains:   pardiso_symmetric_diect                      *
c     ******************************************************************
c
      subroutine pardiso_symmetric_direct
      implicit none
c
c              direct solve: factorization, forward/backward pass
c              combine into one operation.
c
      num_calls = num_calls + 1
      call thyme( 26, 1 )
      iparm(8) = 0 ! max numbers of iterative refinement steps
      phase = 23   ! only forward/backward solve
      CALL pardiso( pt, maxfct, mnum, mtype, phase, neq,
     &              eqn_coeffs, k_pointers, k_indices, idum, nrhs,
     &              iparm, msglvl, rhs, sol_vec, error )
      call warp3d_pardiso_mess( 5, out, error, mkl_ooc_flag,
     &                          print_cpu_stats, iparm )
      call thyme( 26, 2 )
c
      return
c
      end subroutine pardiso_symmetric_direct
c
      end subroutine pardiso_symmetric
c

c
c     ****************************************************************
c     *                                                              *
c     *  Messages during equation solving with Pardiso               *
c     *                                                              *
c     *  last modified: 1/12/2016 rhd                                *
c     *                                                              *
c     ****************************************************************
c
      subroutine warp3d_pardiso_mess( mess_no, iout, ier, ooc_flag,
     &                                cpu_stats, iparm )
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
           write(iout,9900)
           call die_abort
         end if
         if( cpu_stats )  then
           write(iout,9482) wcputime(1)
           if( iparm(18) < 0 ) then
             write(iout,2012)
           else
             write(iout,2010) dble(iparm(18))/1000.d0/1000.d0
           end if
           write(iout,2020) dble(iparm(19))/1000.d0
           write(iout,2022) dble(iparm(15))/1000.d0/1000.d0
         end if
c
        case( 3 )
         if( ier .ne. 0 ) then
            write(iout,9470)
            if( ier .eq. -2 ) write(iout,9700) 
            if( ier .eq. -4 ) write(iout,9710)
            write(iout,9900)
            call die_abort
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
           if( ier .eq. -2 ) write(iout,9700)
           if( ier .eq. -4 ) write(iout,9710)
           write(iout,9900)
           call die_abort
         end if
         if( cpu_stats ) write(iout,9492) wcputime(1)
         if( iparm(23) > 0 ) write(iout,9493) iparm(23)
c
        case( 6 )
         if( ier .ne. 0 ) then
           write(iout,9470)
           write(iout,*) '       >>  phase -1, ier: ',ier
           write(iout,*) '       >>  releasing solver data'
           write(iout,9900)
           call die_abort
         end if
c
        case( 7 )
           write(iout,9470)
           write(iout,*) '  >>  inconsistency @ 1'
           write(iout,*) '  >>  releasing solver data'
           write(iout,9900)
           call die_abort
c
        case( 8 )
         if( ier .ne. 0 ) then
           write(iout,9472)
           write(iout,*) '       >>  phase 52, ier: ',ier
           write(iout,*) '       >>  iterative solve'
           write(iout,9900)
           call die_abort
         end if
         if( cpu_stats ) then
           write(iout,9600) iparm(20)
           write(iout,9602) wcputime(1)
         else
           write(iout,9610) iparm(20)
         end if
c
        case( 9 )
          if( cpu_stats ) write(iout,9286)
c
        case( 10 )
           write(iout,9720)
           write(iout,9900)
           call die_abort
c
        case( 11 )
           write(iout,9740)
           write(iout,9900)
           call die_abort
c
        case default
           write(iout,9730)
           write(iout,9900)
           call die_abort

      end select
      return
c
 2010 format(
     &  15x,'terms in factored matrix (M):    ', f9.2)
 2012 format(
     &  15x,'terms in factored matrix:          > 2.15B')
 2020 format(
     &  15x,'factorization op count (GFlop):  ', f9.2)
 2022 format(
     &  15x,'reordering memory (GB):          ', f9.2)
 2030 format(
     &  15x,'factorization memory (GB):       ', f9.2)
 9470  format(/,
     &  10x, '>>>>>>> FATAL ERRROR: mkl sparse solver' )
 9472  format(/,
     &  10x, '>>>>>>> FATAL ERRROR: mkl sparse -iterative- solver' )
 9480  format(
     &  15x, 'Pardiso solver initializing   @ ',f10.2 )
 9580  format(
     &  15x, ' -> out-of-memory operation' )
 9482  format(
     &  15x, 'reorder-symbolic factor. done @ ',f10.2 )
 9490  format(
     &  15x, 'numeric factorization done    @ ',f10.2 )
 9492  format(
     &  15x, 'numeric factor/solve done     @ ',f10.2 )
 9493  format(/1x,'>>>>> Warning: negative eigenvalues found',
     & /16x,'during equation solve. count: ',i5,/)
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
 9602  format(
     &  15x,'Krylov iterations done        @ ',f10.2 )
 9610  format(7x,
     & '>> conjugate gradient-Krylov iterations:           ',i5)
 9700  format(
     & /,7x,'Pardiso factorization/solve step reported error code -2',
     & /,7x,'This is an out-of-memory error. Usually this means',
     & /,7x,'the factorization step has attempted to allocate more',
     & /,7x,'virtual memory than available on your machine.' )
 9710  format(/,
     & /,7x,'Pardiso factorization/solve step reported error code: -4',
     & /,7x,'There is a zero pivot that cannot be locally re-ordered.')
 9720  format(/,
     &   10x, '>>>>>>> FATAL ERRROR: iterative solver cannot be',
     &     ' used with',
     & /,10x, '              out-of-core option',
     & /,10x  '              Job aborted.')
 9730  format(/,
     &   10x, '>>>>>>> FATAL ERRROR: invalid message condition in',
     & /,10x, '              warp3d_pardiso_mess',
     & /,10x  '              Job aborted.')
 9740  format(/,
     &   10x, '>>>>>>> FATAL ERRROR: invalid itype in',
     &   ' pardiso_symmetric',
     & /,10x  '              or pardiso_unsymmetric. Job aborted.')
 9900 format(/,
     &   10x, '>>>>>>> Job terminated.....',//)
c
      end

c
c     ****************************************************************
c     *                                                              *
c     *  map the default vss sparse solver format that warp3d        *
c     *  assembles into to the MKL sparse solver format              *
c     *                                                              *
c     *  last modified: 11/18/2019                                   *
c     *                                                              *
c     ****************************************************************
c
      subroutine pardiso_symmetric_map( neq, ncoeff, diag, amat,
     &                                  kpt, kind )
c
      implicit none
c
      double precision :: diag(*), amat(*), min_diagonal, max_diagonal
      integer :: neq, ncoeff, kpt(*), kind(*)
c
c      This sub maps arrays needed for the NASA sparse solver to
c      arrays in the format that the MKL solver needs ( Harwell-Boeing
c      format)
c            See format example at the end of this subroutine
c
c      Definitions:
c            ncoeff: # of non-zero terms excluding the diagonal
c            amat  : has the dimension of (neq+ncoeff), however,
c                  the input vector contains the vss format (
c                  i.e. only the ncoeff number of terms); the output
c                  will have all terms filled. Same with indices
c                  array.
c
c       Note: It is assumed that the last term in kpt(neq) == 0
c             is always zero.  Therefore, the second loop will not
c             be performed when i=neq ( See bellow )
c
c
c      Input(vss Format)   -------->     Output (MKL/CSR format)
c           Terms Used                        Terms Used
c
c             amat(ncoeff)                  k_coeff( ncoeff +neq )
c             kind(ncoeff)                  k_indices ( ncoeff + neq )
c              kpt(neq)                      k_pointers ( neq + 1)
c
      double precision ::  zero
      integer, allocatable :: kpt_l(:)
      integer :: i_old, i_new, i, j, nonzt_sum
      data zero / 0.0d0 /
c
      allocate( kpt_l(neq+1) )
      kpt_l(1:neq) = kpt(1:neq)
c
c                 Map amat and kind arrays
c
      i_old = ncoeff + 1; i_new = ncoeff + neq + 1
c
      do i = neq, 1, -1
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
        do j = 1, kpt(i)
          i_new = i_new - 1
          i_old = i_old - 1
          amat(i_new) = amat(i_old)
          kind(i_new) = kind(i_old)
        end do
        if( diag(i) .ne. zero ) then
          i_new = i_new - 1
          amat(i_new) = diag(i)
          kind(i_new) = i
        end if
      end do
c
c                Map kpt ---> k_pointers
c
      nonzt_sum = 1
      kpt(1) = nonzt_sum
      do i = 1, neq
       if( diag(i). ne . zero ) then
         nonzt_sum = nonzt_sum + kpt_l(i) + 1
         kpt(i+1) = nonzt_sum
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
      open(unit=ifileno,file='pardiso_ooc.cfg',
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
