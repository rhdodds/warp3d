c     ****************************************************************
c     *                                                              *
c     *                      subroutine lupres                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 07/08/91                   *
c     *                                                              *
c     *     this subroutine updates the residual for use in the      *
c     *     linear pcg solution algorithm.                           * 
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine lupres( step, stpitr, iter, trcpcg, magres, alpha,
     &                   lres ) 
      implicit integer (a-z)
$add common.main
#dbl      double precision
#sgl      real
     &  alpha, magres, lres(*)
      logical trcpcg
c
#dbl      double precision
#sgl      real
     & zero
      data zero / 0.0 /
c
c
      if( iter .eq. 1 ) then
         magres = zero
         do i = 1,nodof     
            magres = magres+ lres(i)*lres(i)
         end do
         magres = sqrt(magres)
      end if
c
      do i = 1, nodof
         if( cstmap(i) .eq. 0 ) then       
            lres(i) = lres(i) - alpha*ifv(i)
         else
            lres(i) = zero
         end if
      end do
c
c                       output the linear residual if the trace 
c                       flag is on, and if it is desired that  
c                       linear residuals be output.
c                           
      if( trcpcg ) then
         if( prlres ) call oulres( step, stpitr, iter )
      end if
c
      return
      end
