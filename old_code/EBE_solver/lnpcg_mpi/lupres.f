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
c     *     this is the parallel version                             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine lupres( step, stpitr, iter, trcpcg, magres, alpha,
     &     lres_local, mydof, refdof, local2global)
      implicit integer (a-z)
$add common.main
#dbl      double precision
#sgl      real
     &  alpha, magres, lres_local(*), magres_tmp
      logical trcpcg
      dimension local2global(*)
c
#dbl      double precision
#sgl      real
     & zero
      data zero / 0.0 /
c
c
      if( iter .eq. 1 ) then
c
         magres = zero
         do i = 1,mydof
            magres = magres + lres_local(i)*lres_local(i)
         end do
c
         call wmpi_dotprod ( magres )
         magres = sqrt(magres)
c
      end if
c
      do i = 1, refdof
         if( cstmap(local2global(i)) .eq. 0 ) then       
            lres_local(i) = lres_local(i) - alpha*ifv(i)
         else
            lres_local(i) = zero
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
 1000 continue
      return
      end

