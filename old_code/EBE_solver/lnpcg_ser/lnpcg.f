c     ****************************************************************
c     *                                                              *
c     *                      subroutine lnpcg                        *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 06/23/91                   *
c     *                                   11/17/97 AG                *
c     *                                                              *
c     *     this subroutine supervises the computation of the        *
c     *     incremental displacements using a linear pcg             * 
c     *     algorithm for a given iteration of a given time step     * 
c     *     until convergence or the iterative limit has been        * 
c     *     reached. it is also used to calculate the initial        * 
c     *     accelerations for unconstrained dof.                     * 
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine lnpcg( svec, step, stpitr, pcg_precond_time )
c
      use pcg_data, only: pdiag, sdir, lres, lpsres
c
      implicit integer (a-z)              
$add common.main
#dbl      double precision
#sgl      real
     &  svec(*), alpha, beta, z1tr1, magres, magztr
      logical cnverg, trcpcg, cpu_stats, debug
      real t1, pcg_start_time, pcg_precond_time, pcg_end_time,
     &     wcputime
c
#dbl      double precision
#sgl      real
     &     zero
      external wcputime
c
      data cpu_stats, zero, debug / .true., 0.0, .false. /
c
      pcg_start_time = wcputime( 2 )
      trcpcg         = trace(5)
c
      if ( debug ) write (*,*) 'Lets start lnpcg.'
c
c                       initialize the total iteration counter. 
c                             
      iter = 1
c
c                       set the iteration counter for the first iter-
c                       ation after the last restart.
c
 10   rsiter = 1
c
c                       begin iteration loop
c                                     
c                       check for restart. if the iteration number
c                       is greater than the restart iteration, restart
c                       the algorithm.
c                       
 15   if( rsiter .gt. restrt ) go to 10
c
      if (debug) write (*,*) '-> iter ', iter
c
c                       solve for the pseudo residual from the last
c                       iteration, used in computing the new search
c                       direction.
c    
      call thyme( 12, 1 )
c
c                       multiply the new search direction by the 
c                       preconditioner
c
      call ebeslv ( lres, lpsres )
      call thyme( 12, 2 )
c
c                       compute the dot product of the linear residual and
c                       the pseudo residual, and the mutliplier of the the
c                       old search direction, if necessary.
c
      call thyme(13,1)
      call lfindb( step, stpitr, rsiter, iter, trcpcg, z1tr1,
     &             beta, lres, lpsres )  
      call thyme( 13, 2 )
c 
c                       compute the new search direction.
c    
      call thyme( 15, 1 )
      if( rsiter .eq. 1 ) then
        sdir(1:nodof) = lpsres(1:nodof)
      else
        sdir(1:nodof) = lpsres(1:nodof) + beta*sdir(1:nodof)
      end if
      call thyme( 15, 2 )
c
c                       compute the step length for this iteration.
c           
      call thyme( 9, 1 )
      call lfinda( step, stpitr, iter, trcpcg, z1tr1, alpha, sdir )
      call thyme( 9, 2 )
c 
c                       compute the new estimate of the solution 
c                       vector.
c    
      call thyme( 16, 1 )
      svec(1:nodof) = svec(1:nodof) + alpha*sdir(1:nodof)
      call thyme(16,2)
c
c                       update the residual vector for the linear
c                       ebe pcg algorithm.                
c                       
      call thyme( 14, 1 )
      call lupres( step, stpitr, iter, trcpcg, magres, alpha, lres )
      call thyme( 14, 2 )
c
c                       perform any and all convergence tests.
c                       
      call lcvtst( cnverg, step, stpitr, iter, magres, magztr,
     &             z1tr1, lres ) 
c
c                       branch on whether or not convergence has been
c                       met.
c
      if( cnverg ) go to 9999
c
c                       convergence has not been met. update the
c                       iteration numbers and test to see if the to-
c                       tal exceeds the maximum number of iterations
c                       allowed. 
c                       
      iter   = iter + 1      
      rsiter = rsiter + 1      
      if( iter .le. mxlitr ) then
c
c                       total iteration number within prescribed 
c                       limit. return to perform computations for 
c                       the next iteration.
c
         go to 15       
c
      else
c
         call ouconv( 2, iter, stpitr, cnverg, out ) 
c
c                       output timings
c
         t1 = wcputime (1)
         call outime
         write(*,*)
         write(*,*)
         write(*,*) 'total elapsed time (cpu-secs): ',t1
         write(*,*)
         write(*,*) 
         call die_gracefully
         stop
c
      end if
c
c                       output convergence status.
c
 9999 continue
c
      call ouconv( 2, iter, stpitr, cnverg, out )
      pcg_end_time = wcputime( 2 ) - pcg_start_time
      if ( cpu_stats ) then 
         cpu_stats = .false.
         if (pcg_end_time .eq. zero) then
            write(out,9100) pcg_precond_time, pcg_end_time
         else
            cpu_stats = .false.
            write(out,9200) pcg_precond_time, pcg_end_time,
     &           real(iter) / pcg_end_time, 
     &           int(real(iter*noelem) / pcg_end_time) 
         endif
      end if
c
 9100  format (/,
     &  10x, '>> pcg solver cpu statistics:'
     & /15x, 'preconditioner (secs):  ',f10.3,
     & /15x, 'cg iterations (secs):   ',f10.3,/)
 9200  format (/,
     &  10x, '>> pcg solver cpu statistics:'
     & /15x, 'preconditioner (secs):  ',f10.3,
     & /15x, 'cg iterations (secs):   ',f10.3,
     & /15x, 'cg iterations/sec:      ',f10.2,
     & /15x, 'elements/sec:           ',i10 /  )
c
c
      return
      end
