c     ****************************************************************
c     *                                                              *
c     *                      subroutine princ_stress                 *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 10/04/94                   *
c     *                                   03/05/95 kck               *        
c     *                                   06/10/97 rhd               *        
c     *                                                              *
c     *     this subroutine computes principal stress values         *
c     *     at the nodes/elem after the primary values have          *
c     *     been averaged                                            * 
c     *                                                              *
c     ****************************************************************
c
      subroutine princ_stress( results, nrowd, num )
      implicit integer (a-z)
      include 'param_def'
      double precision
     &     results(nrowd,*) 
c
c                    locally allocated
c
      double precision
     &  temp_stress(nstr), wk(ndim), ev(nstr), evec(ndim,ndim)
c
c
c        calculate the principal stresses and there direction cosines by
c        passing off stresses into dummy array, then calling an eigenvalue
c        eigenvector routine, and then passing the values from this routine
c        back to the appropriate places within results
c      
       do i = 1, num
          temp_stress(1) = results(i,1)
          temp_stress(2) = results(i,4)
          temp_stress(3) = results(i,2)
          temp_stress(4) = results(i,6)
          temp_stress(5) = results(i,5)
          temp_stress(6) = results(i,3)
          call ou3dpr( temp_stress, ndim, 1, ev, evec, ndim, wk, ier )   
          results(i,15) = ev(1)     
          results(i,16) = ev(2)
          results(i,17) = ev(3)
          results(i,18) = evec(1,1)
          results(i,19) = evec(2,1)
          results(i,20) = evec(3,1)
          results(i,21) = evec(1,2)
          results(i,22) = evec(2,2)
          results(i,23) = evec(3,2)
          results(i,24) = evec(1,3)
          results(i,25) = evec(2,3)
          results(i,26) = evec(3,3)
      end do
      return
      end  






