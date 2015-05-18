c     ****************************************************************
c     *                                                              *
c     *                      subroutine lcvtst                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 06/23/91                   *
c     *                                                              *
c     *     this subroutine performs the convergence test for the    *
c     *     linear ebe pcg algorithm. if the convergence criteria is * 
c     *     met, the convergence flag is set.                        * 
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine lcvtst( cnverg, step, stpitr, iter, magres, magztr,
     &                   z1tr1, lres ) 
c
      implicit integer (a-z)
$add common.main
#dbl      double precision
#sgl      real
     &  magres, magztr, z1tr1, testvl(2), ratio(2), lres(*)
      logical cnverg,trcsol
c
#dbl      double precision
#sgl      real
     &  zero, hundred
      data zero, hundred / 0.0, 100.0 /
c
c                       initialize and set flags.
c
      trcsol = trace(4)
      cnverg = .false.
c
c                       check for displacement control in the 
c                       convergence tests of the nl solution
c                       algorithm. 
c
      if( convrg(5) ) then
         cnverg = .true.
         go to 9999     
      end if
c
c                       loop over all possible tests.
c
      do 10 i = 1, 2
         if( lcnvrg(i) ) then
            go to (100,200) i
c
c                       test comparing the ratio of the norm of the    
c                       linear residual to the norm of the initial 
c                       linear residual.
c
 100        continue
c
c                       compute the magnitude of linear residual 
c                       vector and compare it to magres.
c
            testvl(1) = zero
            do j = 1, nodof
               testvl(1) = testvl(1) + lres(j)*lres(j)
            end do
            testvl(1) = sqrt(testvl(1))
c
c                       perform convergence test.
c                       
            ratio(1) = (testvl(1)/magres)*hundred
            if (ratio(1) .lt. ltol(1) ) then
               cnverg = .true.
            else
               cnverg = .false.
            end if
c
            go to 10
c
c                       test comparing the ratio of the norm of the    
c                       dot product of the pseudo linear residual
c                       and the linear residual to the norm of the initial 
c                       dot product.
c
 200        continue
c
c                       set the initial dot product norm.
c
            if( iter .eq. 1 ) magztr = sqrt(z1tr1)
c
c                       compute the magnitude of the current dot 
c                       product and compare it to magztr.
c
            testvl(2) = sqrt(z1tr1)
c
c                       perform convergence test.
c                       
            ratio(2) = (testvl(2)/magztr)*hundred
            if( ratio(2). lt. ltol(2) ) cnverg = .true.
c
            go to 10
c
         end if
c
 10   continue
c
c                       output the results of the convergence tests
c                       if desired by the user.
c
      if( trcsol ) then
         call oultrc( step, stpitr, iter, out, lcnvrg, magres,
     &                magztr, testvl, ratio, cnverg )
      end if
c
 9999 return
      end




























