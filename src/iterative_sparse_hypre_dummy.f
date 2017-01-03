c     *****************************************************************
c     *                                                               *
c     *     iterative_sparse_hypre (dummy)                            *
c     *                                                               *
c     *     created by: mcm 2/11                                      *
c     *                                                               *
c     *****************************************************************
      subroutine iterative_sparse_hypre( neq, ncoeff, k_diag, rhs,
     &                        sol_vec, eqn_coeffs, k_pointers,
     &                        k_indices )
      implicit none
c           Input
      integer neq, ncoeff, k_pointers, k_indices
      double precision  k_diag, rhs, sol_vec, eqn_coeffs
      dimension k_diag(*), rhs(*), eqn_coeffs(*), k_pointers(*),
     &          k_indices(*), sol_vec(*)
c           Actually, this error message is not that useful 
c           (not in control flow).  But just in case...
      write(*,*) ' '
      write(*,*) '>> The hypre solvers are not available on this system'
      write(*,*) '>> in this version of WARP3D. please use one of the'
      write(*,*) '>> following keywords with the command: '
      write(*,*) '>> solution technique...'
      write(*,*) '>>   direct iterative, direct sparse  ',
     &  '  = Intel Paradiso solver'
      write(*,*) ' '
      write(*,*) '>> please see the README file in the'   
      write(*,*) '>> ../hypre_solver_dir/ directory for more'
      write(*,*) '>> information on this sparse solver...'
      write(*,*) ' '
      write(*,*) '>> job terminated...'
      write(*,*) ' '
      call die_gracefully
      stop
      end

