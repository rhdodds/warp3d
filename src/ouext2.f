c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine ouext2                       *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 1/19/2017 rhd              *          
c     *                                                              *          
c     *     this subroutine computes derived stress/strain values    *          
c     *     after the primary values have been averaged              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine ouext2 ( results, nrowd, noval, do_stresses )                  
      implicit none                                                             
c                                                                               
      integer :: nrowd, noval                                                   
      double precision :: results(nrowd,*)                                      
      logical :: do_stresses                                                    
c                                                                               
      if( do_stresses ) then                                                    
         call princ_inv_stress ( results, nrowd, noval )                        
         call princ_stress ( results, nrowd, noval )                            
         call yield_function( results, nrowd, noval )                           
      else                                                                      
         call princ_inv_strain ( results, nrowd, noval )                        
         call princ_strain ( results, nrowd, noval )                            
         call equiv_strain( results, nrowd, noval )                             
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine princ_inv_strain             *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 10/04/94                   *          
c     *                                   03/05/95 kck               *          
c     *                                   06/11/97 rhd               *          
c     *                                                              *          
c     *     this subroutine computes principal invariants            *          
c     *     at the nodes/elem after the primary values have          *          
c     *     been averaged                                            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine princ_inv_strain( results, nrowd, num )                        
      implicit integer (a-z)                                                    
      double precision                                                          
     &     results(nrowd,*), half, t1, t2, t3                                   
      data  half  / 0.5  /                                                      
c                                                                               
       do i = 1, num                                                            
         t1 =  half * results(i,4)                                              
         t2 =  half * results(i,5)                                              
         t3 =  half * results(i,6)                                              
         results(i,8) = results(i,1) +  results(i,2) +                          
     &                        results(i,3)                                      
         results(i,9) = - t1 * t1 - t2 * t2 - t3 * t3 +                         
     &                        results(i,1) * results(i,2) +                     
     &                        results(i,2) * results(i,3) +                     
     &                        results(i,1) * results(i,3)                       
         results(i,10) =                                                        
     &     results(i,1) * ( results(i,2) * results(i,3) -                       
     &              t2 * t2 ) -  t1 * ( t1 * results(i,3) -                     
     &              t2 * t3 ) + t3 * ( t1 * t2 -                                
     &              results(i,2) * t3 )                                         
       end do                                                                   
c                                                                               
       return                                                                   
       end                                                                      
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine princ_inv_stress             *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 10/04/94                   *          
c     *                                   03/05/95 kck               *          
c     *                                   06/15/97 rhd               *          
c     *                                                              *          
c     *     this subroutine computes principal invariants            *          
c     *     at the nodes/elements after the primary values have      *          
c     *     been averaged                                            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine princ_inv_stress( results, mxvl, num )                         
      implicit integer (a-z)                                                    
      double precision                                                          
     &     results(mxvl,*)                                                      
                                                                                
       do i = 1, num                                                            
         results(i,12) = results(i,1) +  results(i,2) +                         
     &                        results(i,3)                                      
         results(i,13) = results(i,1) * results(i,2) +                          
     &                        results(i,2) * results(i,3) +                     
     &                        results(i,1) * results(i,3) -                     
     &                        results(i,4) * results(i,4) -                     
     &                        results(i,5) * results(i,5) -                     
     &                        results(i,6) * results(i,6)                       
         results(i,14) = results(i,1) *                                         
     &            ( results(i,2) * results(i,3) -                               
     &              results(i,5) * results(i,5) ) -                             
     &                        results(i,4) *                                    
     &            ( results(i,4) * results(i,3) -                               
     &              results(i,5) * results(i,6) ) +                             
     &                        results(i,6) *                                    
     &            ( results(i,4) * results(i,5) -                               
     &              results(i,2) * results(i,6) )                               
       end do                                                                   
c                                                                               
       return                                                                   
       end                                                                      
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
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
      include 'param_def'                                                       
      double precision                                                          
     &     results(nrowd,*)                                                     
c                                                                               
c                    locally allocated                                          
c                                                                               
      double precision                                                          
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
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine yield_function               *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 10/04/94                   *          
c     *                                   03/05/95 kck               *          
c     *                                                              *          
c     *     this subroutine computes mises yield function value      *          
c     *     at the nodes/elements after the primary values have      *          
c     *     been averaged                                            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine yield_function ( results, maxnum, num )                        
      implicit integer (a-z)                                                    
      double precision                                                          
     &     results(maxnum,*), iroot2, six                                       
      data iroot2, six / 0.70711, 6.0 /                                         
c                                                                               
      do i = 1, num                                                             
        results(i,8) = sqrt(                                                    
     &     ( results(i,1) - results(i,2) ) ** 2 +                               
     &     ( results(i,2) - results(i,3) ) ** 2 +                               
     &     ( results(i,1) - results(i,3) ) ** 2 +                               
     &   six*( results(i,4)**2 + results(i,5)**2+                               
     &         results(i,6)**2 ) )*iroot2                                       
c                                                                               
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine equiv_strain                 *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 10/04/94                   *          
c     *                                   06/02/95 kck               *          
c     *                                                              *          
c     *     this subroutine computes mises yield function value      *          
c     *     at the nodes/elem after the primary values have          *          
c     *     been averaged                                            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine equiv_strain ( results, maxnum, num )                          
      implicit integer (a-z)                                                    
      double precision                                                          
     &     results(maxnum,*),  root23, onep5                                    
      data root23, onep5 / 0.471404, 1.5 /                                      
c                                                                               
      do i = 1, num                                                             
        results(i,7) = root23 * sqrt(                                           
     &     ( results(i,1) - results(i,2) ) ** 2 +                               
     &     ( results(i,2) - results(i,3) ) ** 2 +                               
     &     ( results(i,1) - results(i,3) ) ** 2 +                               
     &   onep5 * ( results(i,4)**2 + results(i,5)**2+                           
     &         results(i,6)**2 ) )                                              
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c *      element output service routine -- ou3dpr                   *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine ou3dpr( a,n,ijob,d,z,iz,wk,ier )                               
      implicit double precision (a-h,o-z)                                       
c   purpose             - eigenvalues and (optionally) eigenvectors of          
c                           a real symmetric matrix in symmetric                
c                           storage mode                                        
c                                                                               
c   arguments    a      - input real symmetric matrix of order n,               
c                           stored in symmetric storage mode,                   
c                           whose eigenvalues and eigenvectors                  
c                           are to be computed. input a is                      
c                           destroyed if ijob is equal to 0 or 1.               
c                n      - input order of the matrix a.                          
c                ijob   - input option parameter, when                          
c                           ijob = 0, compute eigenvalues only                  
c                           ijob = 1, compute eigenvalues and eigen-            
c                             vectors.                                          
c                           ijob = 2, compute eigenvalues, eigenvectors         
c                             and performance index.                            
c                           ijob = 3, compute performance index only.           
c                           if the performance index is computed, it is         
c                           returned in wk(1). the routines have                
c                           performed (well, satisfactorily, poorly) if         
c                           wk(1) is (less than 1, between 1 and 100,           
c                           greater than 100).                                  
c                d      - output vector of length n,                            
c                           containing the eigenvalues of a.                    
c                z      - output n by n matrix containing                       
c                           the eigenvectors of a.                              
c                           the eigenvector in column j of z corres-            
c                           ponds to the eigenvalue d(j).                       
c                           if ijob = 0, z is not used.                         
c                iz     - input row dimension of matrix z exactly as            
c                           specified in the dimension statement in the         
c                           calling program.                                    
c                wk     - work area, the length of wk depends                   
c                           on the value of ijob, when                          
c                           ijob = 0, the length of wk is at least n.           
c                           ijob = 1, the length of wk is at least n.           
c                           ijob = 2, the length of wk is at least              
c                             n(n+1)/2+n.                                       
c                           ijob = 3, the length of wk is at least 1.           
c                ier    - error parameter (output)                              
c                         terminal error                                        
c                           ier = 128+j, indicates that ourt2s failed           
c                             to converge on eigenvalue j. eigenvalues          
c                             and eigenvectors 1,...,j-1 have been              
c                             computed correctly, but the eigenvalues           
c                             are unordered. the performance index              
c                             is set to 1000.0                                  
c                         warning error (with fix)                              
c                           ier = 66, indicates ijob is less than 0 or          
c                             ijob is greater than 3. ijob set to 1.            
c                           ier = 67, indicates ijob is not equal to            
c                             zero, and iz is less than the order of            
c                             matrix a. ijob is set to zero.                    
c                                                                               
c                                  specifications for arguments                 
      integer            n,ijob,iz,ier                                          
      double precision                                                          
     &                   a(*),d(*),wk(*),z(iz,*)                                
c                                  specifications for local variables           
      integer            jer,na,nd,iiz,ibeg,il,kk,lk,i,j,k,l                    
      double precision                                                          
     &                   anorm,asum,pi,sumz,sumr,an,s,ten,rdelp,zero,           
     &                   one,thous                                              
      data               rdelp /0.745058d-08/                                   
      data               zero,one/0.d0,1.d0/,ten/10.d0/,thous/1000.d0/          
c                                  initialize error parameters                  
c                                  first executable statement                   
      ier = 0                                                                   
      jer = 0                                                                   
      if (ijob .ge. 0 .and. ijob .le. 3) go to 5                                
c                                  warning error - ijob is not in the           
c                                    range                                      
      ier = 66                                                                  
      ijob = 1                                                                  
      go to 10                                                                  
    5 if (ijob .eq. 0) go to 20                                                 
   10 if (iz .ge. n) go to 15                                                   
c                                  warning error - iz is less than n            
c                                    eigenvectors can not be computed,          
c                                    ijob set to zero                           
      ier = 67                                                                  
      ijob = 0                                                                  
   15 if (ijob .eq. 3) go to 65                                                 
   20 na = (n*(n+1))/2                                                          
      if (ijob .ne. 2) go to 35                                                 
      do 30 i=1,na                                                              
         wk(i) = a(i)                                                           
   30 continue                                                                  
c                                  save input a if ijob = 2                     
   35 nd = 1                                                                    
      if (ijob .eq. 2) nd = na+1                                                
c                                  reduce a to symmetric tridiagonal            
c                                    form                                       
      call ouhous (a,n,d,wk(nd),wk(nd))                                         
      iiz = 1                                                                   
      if (ijob .eq. 0) go to 50                                                 
      iiz = iz                                                                  
c                                  set z to the identity matrix                 
      do 45 i=1,n                                                               
         do 40 j=1,n                                                            
            z(i,j) = zero                                                       
   40    continue                                                               
         z(i,i) = one                                                           
   45 continue                                                                  
c                                  compute eigenvalues and eigenvectors         
   50 call ourt2s (d,wk(nd),n,z,iiz,jer)                                        
      if (ijob .eq. 0) go to 9000                                               
      if (jer .gt. 128) go to 55                                                
c                                  back transform eigenvectors                  
      call ouobks (a,n,1,n,z,iz)                                                
   55 if (ijob .le. 1) go to 9000                                               
c                                  move input matrix back to a                  
      do 60 i=1,na                                                              
         a(i) = wk(i)                                                           
   60 continue                                                                  
      wk(1) = thous                                                             
      if (jer .ne. 0) go to 9000                                                
c                                  compute 1 - norm of a                        
   65 anorm = zero                                                              
      ibeg = 1                                                                  
      do 75 i=1,n                                                               
         asum = zero                                                            
         il = ibeg                                                              
         kk = 1                                                                 
         do 70 l=1,n                                                            
            asum =asum+abs(a(il))                                               
            if (l .ge. i) kk = l                                                
            il = il+kk                                                          
   70    continue                                                               
         anorm = dmax1(anorm,asum)                                              
         ibeg = ibeg+i                                                          
   75 continue                                                                  
      if (anorm .eq. zero) anorm = one                                          
c                                  compute performance index                    
      pi = zero                                                                 
      do 90 i=1,n                                                               
         ibeg = 1                                                               
         s = zero                                                               
         sumz = zero                                                            
         do 85 l=1,n                                                            
            lk = ibeg                                                           
            kk = 1                                                              
            sumz = sumz+abs(z(l,i))                                             
            sumr = -d(i)*z(l,i)                                                 
            do 80 k=1,n                                                         
               sumr = sumr+a(lk)*z(k,i)                                         
               if (k .ge. l) kk = k                                             
               lk = lk+kk                                                       
   80       continue                                                            
            s = s+abs(sumr)                                                     
            ibeg = ibeg+l                                                       
   85    continue                                                               
         if (sumz .eq. zero) go to 90                                           
         pi = dmax1(pi,s/sumz)                                                  
   90 continue                                                                  
      an = n                                                                    
      pi = pi/(anorm*ten*an*rdelp)                                              
      wk(1) = pi                                                                
 9000 continue                                                                  
      if ( ier .ne. 0 ) return                                                  
      if (jer .eq. 0) go to 9005                                                
      ier = jer                                                                 
      return                                                                    
 9005 return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *      element output service routine -- ouobks                   *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine ouobks( a,n,m1,m2,z,iz )                                       
      implicit double precision (a-h,o-z)                                       
c                                                                               
c   purpose             - back transformation to form the eigenvectors          
c                           of the original symmetric matrix from the           
c                           eigenvectors of the tridiagonal matrix              
c                                                                               
c   arguments    a      - the array contains the details of the house-          
c                           holder reduction of the original matrix a           
c                           as generated by routine ouhous.(input)              
c                n      - order of the real symmetric matrix.(input)            
c                m1     - m1 and m2 are two input scalars such that             
c                           eigenvectors m1 to m2 of the tridiagonal            
c                           matrix a have been found and normalized             
c                           according to the euclidean norm.                    
c                m2     - see above - m1                                        
c                z      - a two dimensional array of size n x (m2-m1+1)         
c                           which contains eigenvectors m1 to m2 of             
c                           tridiagonal matrix t, normalized according          
c                           to euclidean norm. input z can be produced          
c                           by routine ourt2s, the resultant                    
c                           matrix overwrites the input z.(input/output)        
c                iz     - row dimension of matrix z exactly as                  
c                           specified in the dimension statement in the         
c                           calling program. (input)                            
c                                                                               
c                                                                               
      double precision                                                          
     &                   zero                                                   
      dimension          a(*),z(iz,*)                                           
      data zero /0.0/                                                           
                                                                                
c                                  first executable statement                   
      if (n.eq.1) go to 30                                                      
      do 25 i=2,n                                                               
         l = i-1                                                                
         ia = (i*l)/2                                                           
         h = a(ia+i)                                                            
         if (h.eq.zero) go to 25                                                
c                                  derives eigenvectors m1 to m2 of             
c                                  the original matrix from eigenvectors        
c                                  m1 to m2 of the symmetric                    
c                                  tridiagonal matrix                           
         do 20 j = m1,m2                                                        
            s = zero                                                            
            do 10 k = 1,l                                                       
               s = s+a(ia+k)*z(k,j)                                             
   10       continue                                                            
            s = s/h                                                             
            do 15 k=1,l                                                         
               z(k,j) = z(k,j)-s*a(ia+k)                                        
   15       continue                                                            
   20    continue                                                               
   25 continue                                                                  
   30 return                                                                    
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c *      element output service routine -- ouhous                   *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine ouhous( a,n,d,e,e2 )                                           
      implicit double precision (a-h,o-z)                                       
c                                                                               
c   purpose             - reduction of a symmetric matrix to symmetric          
c                           tridiagonal form using a householder                
c                           reduction                                           
c                                                                               
c   arguments    a      - the given n x n, real symmetric matrix a,             
c                           where a is stored in symmetric storage mode.        
c                           the input a is replaced by the details of           
c                           the householder reduction of a.                     
c                n      - input order of a and the length of d, e, and          
c                           e2.                                                 
c                d      - the output array of length n, giving the              
c                           diagonal elements of the tridiagonal matrix.        
c                e      - the output array of length n, giving the sub-         
c                           diagonal in the last (n-1) elements, e(1) is        
c                           set to zero.                                        
c                e2     - output array of length n.  e2(i) = e(i)**2.           
c                                                                               
c                                                                               
      dimension          a(*),d(n),e(n),e2(n)                                   
      double precision                                                          
     &                   a,d,e,e2,zero,h,scale,one,scale1,f,g,hh                
      data               zero/0.d0/,one/1.d0/                                   
c                                  first executable statement                   
      np1 = n+1                                                                 
      nn = (n*np1)/2-1                                                          
      nbeg = nn+1-n                                                             
      do 70 ii = 1,n                                                            
         i = np1-ii                                                             
         l = i-1                                                                
         h = zero                                                               
         scale = zero                                                           
         if (l .lt. 1) go to 10                                                 
c                                  scale row (algol tol then not needed)        
         nk = nn                                                                
         do 5 k = 1,l                                                           
            scale = scale+abs(a(nk))                                            
            nk = nk-1                                                           
    5    continue                                                               
         if (scale .ne. zero) go to 15                                          
   10    e(i) = zero                                                            
         e2(i) = zero                                                           
         go to 65                                                               
   15    nk = nn                                                                
         scale1 = one/scale                                                     
         do 20 k = 1,l                                                          
            a(nk) = a(nk)*scale1                                                
            h = h+a(nk)*a(nk)                                                   
            nk = nk-1                                                           
   20    continue                                                               
         e2(i) = scale*scale*h                                                  
         f = a(nn)                                                              
         g = -sign(sqrt(h),f)                                                   
         e(i) = scale*g                                                         
         h = h-f*g                                                              
         a(nn) = f-g                                                            
         if (l .eq. 1) go to 55                                                 
         f = zero                                                               
         jk1 = 1                                                                
         do 40 j = 1,l                                                          
            g = zero                                                            
            ik = nbeg+1                                                         
            jk = jk1                                                            
c                                  form element of a*u                          
            do 25 k = 1,j                                                       
               g = g+a(jk)*a(ik)                                                
               jk = jk+1                                                        
               ik = ik+1                                                        
   25       continue                                                            
            jp1 = j+1                                                           
            if (l .lt. jp1) go to 35                                            
            jk = jk+j-1                                                         
            do 30 k = jp1,l                                                     
               g = g+a(jk)*a(ik)                                                
               jk = jk+k                                                        
               ik = ik+1                                                        
   30       continue                                                            
c                                  form element of p                            
   35       e(j) = g/h                                                          
            f = f+e(j)*a(nbeg+j)                                                
            jk1 = jk1+j                                                         
   40    continue                                                               
         hh = f/(h+h)                                                           
c                                  form reduced a                               
         jk = 1                                                                 
         do 50 j = 1,l                                                          
            f = a(nbeg+j)                                                       
            g = e(j)-hh*f                                                       
            e(j) = g                                                            
            do k = 1,j                                                          
               a(jk) = a(jk)-f*e(k)-g*a(nbeg+k)                                 
               jk = jk+1                                                        
            enddo                                                               
   50    continue                                                               
   55    do 60 k = 1,l                                                          
            a(nbeg+k) = scale*a(nbeg+k)                                         
   60    continue                                                               
   65    d(i) = a(nbeg+i)                                                       
         a(nbeg+i) = h*scale*scale                                              
         nbeg = nbeg-i+1                                                        
         nn = nn-i                                                              
   70 continue                                                                  
      return                                                                    
      end                                                                       
c *******************************************************************           
c *                                                                 *           
c *      element output service routine -- ourt2s                               
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine ourt2s(d,e,n,z,iz,ier)                                         
      implicit double precision (a-h,o-z)                                       
c   purpose             - eigenvalues and (optionally) eigenvectors of          
c                           a symmetric tridiagonal matrix using the            
c                           ql method.                                          
c                                                                               
c   arguments    d      - on input, the vector d of length n contains           
c                           the diagonal elements of the symmetric              
c                           tridiagonal matrix t.                               
c                           on output, d contains the eigenvalues of            
c                           t in ascending order.                               
c                e      - on input, the vector e of length n contains           
c                           the sub-diagonal elements of t in position          
c                           2,...,n. on output, e is destroyed.                 
c                n      - order of tridiagonal matrix t.(input)                 
c                z      - on input, z contains the identity matrix of           
c                           order n.                                            
c                           on output, z contains the eigenvectors              
c                           of t. the eigenvector in column j of z              
c                           corresponds to the eigenvalue d(j).                 
c                iz     - input row dimension of matrix z exactly as            
c                           specified in the dimension statement in the         
c                           calling program. if iz is less than n, the          
c                           eigenvectors are not computed. in this case         
c                           z is not used.                                      
c                ier    - error parameter                                       
c                         terminal error                                        
c                           ier = 128+j, indicates that ourt2s failed           
c                             to converge on eigenvalue j. eigenvalues          
c                             and eigenvectors 1,...,j-1 have been              
c                             computed correctly, but the eigenvalues           
c                             are unordered.                                    
c                                                                               
c                                                                               
      dimension          d(*),e(*),z(iz,*)                                      
      data               rdelp /0.745058e-08/                                   
      data               zero,one/0.0,1.0/                                      
c                                  move the last n-1 elements                   
c                                  of e into the first n-1 locations            
c                                  first executable statement                   
      ier  = 0                                                                  
      if (n .eq. 1) go to 9005                                                  
      do 5  i=2,n                                                               
         e(i-1) = e(i)                                                          
    5 continue                                                                  
      e(n) = zero                                                               
      b = zero                                                                  
      f = zero                                                                  
      do  60  l=1,n                                                             
         j = 0                                                                  
         h = rdelp*(abs(d(l))+abs(e(l)))                                        
         if (b.lt.h) b = h                                                      
c                                  look for small sub-diagonal element          
         do 10  m=l,n                                                           
            k=m                                                                 
            if (abs(e(k)) .le. b) go to 15                                      
   10    continue                                                               
   15    m = k                                                                  
         if (m.eq.l) go to 55                                                   
   20    if (j .eq. 30) go to 85                                                
         j = j+1                                                                
         l1 = l+1                                                               
         g = d(l)                                                               
         p = (d(l1)-g)/(e(l)+e(l))                                              
         r = sqrt(p*p+one)                                                      
         d(l) = e(l)/(p+sign(r,p))                                              
         h = g-d(l)                                                             
         do 25 i = l1,n                                                         
            d(i) = d(i)-h                                                       
   25    continue                                                               
         f = f+h                                                                
c                                  ql transformation                            
         p = d(m)                                                               
         c = one                                                                
         s = zero                                                               
         mm1 = m-1                                                              
         mm1pl = mm1+l                                                          
         if (l.gt.mm1) go to 50                                                 
         do 45 ii=l,mm1                                                         
            i = mm1pl-ii                                                        
            G = C*E(I)                                                          
            H = C*P                                                             
            if (abs(p).lt.abs(e(i))) go to 30                                   
            c = e(i)/p                                                          
            r = sqrt(c*c+one)                                                   
            e(i+1) = s*p*r                                                      
            s = c/r                                                             
            c = one/r                                                           
            go to 35                                                            
   30       c = p/e(i)                                                          
            r = sqrt(c*c+one)                                                   
            e(i+1) = s*e(i)*r                                                   
            s = one/r                                                           
            c = c*s                                                             
   35       p = c*d(i)-s*g                                                      
            d(i+1) = h+s*(c*g+s*d(i))                                           
            if (iz .lt. n) go to 45                                             
c                                  form vector                                  
            do 40 k=1,n                                                         
               h = z(k,i+1)                                                     
               z(k,i+1) = s*z(k,i)+c*h                                          
               z(k,i) = c*z(k,i)-s*h                                            
   40       continue                                                            
   45    continue                                                               
   50    e(l) = s*p                                                             
         d(l) = c*p                                                             
         if ( abs(e(l)) .gt.b) go to 20                                         
   55    d(l) = d(l) + f                                                        
   60 continue                                                                  
c                                  order eigenvalues and eigenvectors           
      do  80  i=1,n                                                             
         k = i                                                                  
         p = d(i)                                                               
         ip1 = i+1                                                              
         if (ip1.gt.n) go to 70                                                 
         do 65  j=ip1,n                                                         
            if (d(j) .ge. p) go to 65                                           
            k = j                                                               
            p = d(j)                                                            
   65    continue                                                               
   70    if (k.eq.i) go to 80                                                   
         d(k) = d(i)                                                            
         d(i) = p                                                               
         if (iz .lt. n) go to 80                                                
         do 75 j = 1,n                                                          
            p = z(j,i)                                                          
            z(j,i) = z(j,k)                                                     
            z(j,k) = p                                                          
   75    continue                                                               
   80 continue                                                                  
      go to 9005                                                                
   85 ier = 128+l                                                               
 9000 continue                                                                  
 9005 return                                                                    
      end                                                                       
                                                                                
