c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine rscmp1                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 06/30/91                   *          
c     *                                                              *          
c     *     this subroutine computes the polar decompostion of the   *          
c     *     deformation gradient into the right stretch tensor.      *          
c     *     the computations are for a gauss point of a q3disop      *          
c     *     or l3disop element.                                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine rscmp1( span, f, u )                                           
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
c                                                                               
c                      parameter declarations                                   
c                                                                               
      double precision                                                          
     & f(mxvl,ndim,*), u(mxvl,*)                                                
c                                                                               
c                      locally allocated arrays & constants                     
c                                                                               
      double precision                                                          
     &  c(mxvl,nstr), cc(mxvl,nstr), iu(mxvl), iiu(mxvl),                       
     &  iiiu(mxvl), a1(mxvl), b1(mxvl), c1(mxvl), one                           
      data one / 1.0 /                                                          
c                                                                               
c                       u is in symmetric upper triangular form.                
c                                                                               
c                                                                               
c                       compute the invariants of the right                     
c                       stretch tensor, the metric tensor, and                  
c                       its square.                                             
c                                                                               
      call ivcmp1( span, f, c, cc, iu, iiu, iiiu )                              
c                                                                               
c                       compute multipliers.                                    
c                                                                               
      do i = 1, span                                                            
         a1(i)= one/(iu(i)*iiu(i)-iiiu(i))                                      
         b1(i)= iu(i)*iiiu(i)                                                   
         c1(i)= iu(i)*iu(i)-iiu(i)                                              
      end do                                                                    
c                                                                               
c                       compute the right stretch tensor.                       
c                                                                               
      do i = 1, span                                                            
         u(i,1)= a1(i) * ( b1(i) + c1(i)*c(i,1) - cc(i,1) )                     
         u(i,2)= a1(i) * (         c1(i)*c(i,2) - cc(i,2) )                     
         u(i,3)= a1(i) * ( b1(i) + c1(i)*c(i,3) - cc(i,3) )                     
         u(i,4)= a1(i) * (         c1(i)*c(i,4) - cc(i,4) )                     
         u(i,5)= a1(i) * (         c1(i)*c(i,5) - cc(i,5) )                     
         u(i,6)= a1(i) * ( b1(i) + c1(i)*c(i,6) - cc(i,6) )                     
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
