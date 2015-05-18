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
c     *     This is the parallel version, where each processor       *
c     *     does the operations for section of the total mesh        *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine lnpcg( svec, step, stpitr, pcg_precond_time )
      use pcg_data, only: pdiag, sdir, lres, lpsres, diag
      use mpi_lnpcg
      implicit integer (a-z)              
$add common.main
#dbl      double precision
#sgl      real
     &   svec(*), alpha, beta, z1tr1, magres, magztr, svec_local(mxdof)
      logical cnverg, trcpcg, cpu_stats, debug, ebe_stats
      real t1, pcg_start_time, pcg_precond_time, pcg_end_time,
     &     wcputime, time1, time2, ebetime, alltot
      external wcputime
c
#dbl      double precision
#sgl      real
c     &     zero, diag_local(mxdof), lres_local(mxdof)
     &     zero, lres_local(mxdof) 
c
      data cpu_stats, zero, debug, ebe_stats
     &     / .true., 0.0, .false., .true. /
c
c
      pcg_start_time = wcputime( 2 )
      trcpcg = trace(5)
      ebetime = zero
c
      if (debug) write (*,*) 'Lets start lnpcg.'
c
c                       tell the slave processors to start lnpcg.  Then
c                       call dd_pcg_init to intialize the processor
c                       local variables for the lnpcg solver.
c
      call wmpi_alert_slaves ( 5 ) 
      call dd_pcg_init( lres_local, svec_local, svec,
     &     mydof, refdof)
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
c     
      if (debug) write (*,*) myid,': -> iter ', iter
c
c                       solve for the pseudo residual from the last
c                       iteration, used in computing the new search
c                       direction.
c    
      call thyme( 12, 1 )
c
c                       multiply the new search direction by the 
c                         preconditioner
c
      if (debug) write (*,*) myid,':-> ebeslv '
c
      if (ebe_stats) time1 = wcputime (1)
c
      call ebeslv (lres_local,lpsres,diag,refdof,mydof,iter)
c
      if (ebe_stats) then
         time2 = wcputime (1)
         ebetime = ebetime + time2 - time1
      endif
c
      call thyme( 12, 2 )
c
c                       compute the dot product of the linear residual and
c                       the pseudo residual, and the mutliplier of the the
c                       old search direction, if necessary.
c
      call thyme(13,1)
      if (debug) write (*,*) myid,':-> lfindb '
      call lfindb( step, stpitr, rsiter, iter, trcpcg, z1tr1,
     &             beta, lres_local, lpsres, mydof)  
      if (debug) write (*,*) myid,':==> beta:', beta
      call thyme( 13, 2 )
c 
c                       compute the new search direction.
c    
      call thyme( 15, 1 )
c     
      if( rsiter .eq. 1 ) then
         sdir(1:refdof) = lpsres(1:refdof)
      else
         sdir(1:refdof) = lpsres(1:refdof) + beta*sdir(1:refdof)
      end if
c    
      call thyme( 15, 2 )
c
c                       compute the step length for this iteration.
c           
      call thyme( 9, 1 )
c                                            
      if (debug) write (*,*) myid,':-> lfinda '
c
      call lfinda( step, stpitr, iter, trcpcg, z1tr1, alpha, sdir, 
     &     mydof, refdof)
c
      if (debug) write (*,*) myid,':    alpha is ', alpha
c    
c
      call thyme( 9, 2 )
c 
c                       compute the new estimate of the solution 
c                       vector.
c    
      call thyme( 16, 1 )
c                       
      if ( debug) write (*,*) myid,':>>>> doing svec_local'
c
      svec_local(1:mydof) = svec_local(1:mydof) + alpha*sdir(1:mydof)
c    
      call thyme(16,2)
c
c                       update the residual vector for the linear
c                       ebe pcg algorithm.                
c                       
      call thyme( 14, 1 )
c     
      if ( debug) write (*,*) myid,':>>>> do lupres'
c
      call lupres( step, stpitr, iter, trcpcg, magres, alpha,
     &             lres_local, mydof, refdof, local_nodes%local2global )
c     
      call thyme( 14, 2 )
c
c                       perform any and all convergence tests.
c                       
      if ( debug) write (*,*) myid,':>>>> do lcvtst'
c
      call lcvtst( cnverg, step, stpitr, iter, magres, magztr,
     &     z1tr1, lres_local, mydof ) 
c
c                       branch on whether or not convergence has been
c                       met.
c
      if ( cnverg ) go to 9999
c
c                       convergence has not been met. update the
c                       iteration numbers and test to see if the to-
c                       tal exceeds the maximum number of iterations
c                       allowed. 
c                       
      iter   = iter + 1      
      rsiter = rsiter + 1      
c
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
         if (myid .eq. 0) then
         
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
         else
            return
         endif
c
      end if
c
c
 9999 continue
c
c                       we have converged. collect svec terms all back
c                       to processor zero.
c
      if ( debug)write (*,*) '>>>>>>>>>>>>>>> Hey, we made it!'
      call dd_combine (svec_local, svec)
c
c                       output convergence status.
c
      if ( slave_processor ) return
c
      call ouconv( 2, iter, stpitr, cnverg, out )
      pcg_end_time = wcputime( 2 ) - pcg_start_time
      if ( cpu_stats ) then 
c         cpu_stats = .false.
         if (pcg_end_time .eq. zero) then
            write(out,9100) pcg_precond_time, pcg_end_time
         else
c            cpu_stats = .false.
            write(out,9200) pcg_precond_time, pcg_end_time,
     &           real(iter) / pcg_end_time, 
     &           int(real(iter*noelem) / pcg_end_time) ,
     &           ebetime/pcg_end_time
         endif
      end if
c     
 9100 format (/,
     &  10x, '>> pcg solver cpu statistics:'
     & /15x, 'preconditioner (secs):  ',f10.3,
     & /15x, 'cg iterations (secs):   ',f10.3,/)
 9200 format (/,
     &  10x, '>> pcg solver cpu statistics:'
     & /15x, 'preconditioner (secs):  ',f10.3,
     & /15x, 'cg iterations (secs):   ',f10.3,
     & /15x, 'cg iterations/sec:      ',f10.2,
     & /15x, 'elements/sec:           ',i10,
     & /15x, '% time for h-w prec:    ',f10.5,/)
c
c
      return
      end






