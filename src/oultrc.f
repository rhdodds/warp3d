c     ****************************************************************
c     *                                                              *
c     *                      subroutine oultrc                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 10/07/90                   *
c     *                                 : 02/20/94                   *
c     *                                                              *
c     *     this subroutine outputs the status of the convergence    *
c     *     test in order to trace the soluton using the linear      * 
c     *     pcg algorithm.                                           * 
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine oultrc( step, stpitr, iter, out, lcnvrg, magres,
     &                  magztr, testvl, ratio, cnverg )
      implicit integer (a-z)
#dbl      double precision
#sgl      real
     & magres, magztr, testvl(*), ratio(*)
      logical lcnvrg(*), cnverg
      real time, wcputime 
      external wcputime
c
c                      we only write a lpcg trace message
c                      every frequency iterations
c
      time = wcputime( 1 )
      frequency = 50
      if ( mod(iter,frequency) .ne. 0 .and. .not. cnverg) return
c
      write(out,9000) time, step, stpitr, iter
c
c                       test comparing the ratio of the norm of the    
c                       linear residual vector to the norm of initial
c                       linear residual.
c
      if( lcnvrg(1) ) then
            write(out,9010) testvl(1),magres,ratio(1)
      end if
c                       
c                       test comparing the ratio of the norm of the    
c                       dot product of the pseudo linear residual
c                       and the linear residual to the norm of the initial 
c                       dot product.
c  
      if( lcnvrg(2) ) then
            write(out,9020) testvl(2),magztr,ratio(2)
      end if
c
 9000 format(/10x,'lpcg convergence tests @ wall time:',f8.1,
     &      /,10x,28('-'),/,13x,
     &       'step: ',i8,1x,'nonlinear iter: ',i3,
     &       1x,'lpcg iter: ',i4 )
c
 9010 format(13x,'norm of linear residual vector:          ',e12.5,
     &      /13x,'norm of original linear residual vector: ',e12.5,
     &      /13x,'ratio*100: ',f15.5)
c
 9020 format(13x,'norm of linear dot product:          ',e12.5,
     &      /13x,'norm of original linear dot product: ',e12.5,
     &      /13x,'ratio*100: ',f15.5)
c
      return
      end
