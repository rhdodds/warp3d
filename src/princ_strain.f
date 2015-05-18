c     ****************************************************************
c     *                                                              *
c     *                      subroutine princ_strain                 *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 10/04/94                   *
c     *                                   03/05/95 kck               *        
c     *                                   06/10/97 rhd               *        
c     *                                                              *
c     *     this subroutine computes principal strain values         *
c     *     at the nodes/elem after the primary values have          *
c     *     been averaged                                            * 
c     *                                                              *
c     ****************************************************************
c
      subroutine princ_strain( results, nrowd, num )
      implicit integer (a-z)
$add param_def
#dbl      double precision
#sgl      real
     &     results(nrowd,*)
c
c                    locally allocated
c
#dbl      double precision
#sgl      real
     &  temp_strain(nstr), wk(ndim), ev(nstr), evec(ndim,ndim), half
c
      data  half / 0.5 / 
c
c
c        calculate the principal strains and there direction cosines by
c        passing off strains into dummy array, then calling an eigenvalue
c        eigenvector routine, and then passing the values from this routine
c        back to the appropriate places within results
c      
       do i = 1, num
          temp_strain(1) = results(i,1)
          temp_strain(2) = results(i,4) * half
          temp_strain(3) = results(i,2)
          temp_strain(4) = results(i,6) * half
          temp_strain(5) = results(i,5) * half
          temp_strain(6) = results(i,3)
          call ou3dpr( temp_strain, ndim, 1, ev, evec, ndim, wk, ier )   
          results(i,11) = ev(1)     
          results(i,12) = ev(2)
          results(i,13) = ev(3)
          results(i,14) = evec(1,1)
          results(i,15) = evec(2,1)
          results(i,16) = evec(3,1)
          results(i,17) = evec(1,2)
          results(i,18) = evec(2,2)
          results(i,19) = evec(3,2)
          results(i,20) = evec(1,3)
          results(i,21) = evec(2,3)
          results(i,22) = evec(3,3)
      end do
      return
      end  









