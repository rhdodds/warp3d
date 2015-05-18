c     ****************************************************************
c     *                                                              *
c     *                      subroutine lfindb                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 10/07/90                   *
c     *                                                              *
c     *     this subroutine computes multiplier of the previous      *
c     *     search direction for use in computing the new search     *
c     *     direction in the linear pcg algorithm.                   *
c     *                                                              *
c     ****************************************************************
c
c                       
c
      subroutine lfindb( step, stpitr, rsiter, iter, trcpcg, z1tr1,
     &                   beta, lres, lpsres ) 
      implicit integer (a-z)
#dbl      double precision
#sgl      real
     &    beta, z1tr1, z2tr2, lres(*), lpsres(*)     
      logical trcpcg
$add common.main
c
#dbl      double precision
#sgl      real
     &  zero
      data zero / 0.0 /
c
c                       if it is iteration one, perform only the
c                       computation of z1tr1. otherwise, perform all
c                       computations necessary to find beta.
c                       
      if( rsiter .gt. 1 ) then 
c
c                       assign dot product from the previous iteration.
c
         z2tr2 = z1tr1
c
c                       compute the dot product for the current iteration.
c                            
         z1tr1 = zero
         do i = 1, nodof
            z1tr1 = z1tr1 + lres(i)*lpsres(i)
         end do
c
c                       compute beta.
c                       
         beta = z1tr1 / z2tr2
c
c                       output the value of beta if the trace flag
c                       is on.
c                       
         if( trcpcg ) then
            call oulbta( step, stpitr, iter, beta )
         end if
c
      else
c
c                       compute the dot product for the current iteration.
c                            
         z1tr1 = zero
         beta = zero
         do i = 1, nodof
            z1tr1 = z1tr1 + lres(i)*lpsres(i)
         end do
c
      end if
c
c
      return
      end

