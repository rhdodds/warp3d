c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine blcmp1                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 04/26/2017 rhd             *          
c     *                                                              *          
c     *     computes the linear strain-displacement matrices at a    *          
c     *     given gauss point for element block                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine blcmp1( span, b, gama, nxi, neta, nzeta, nnode )               
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c                       parameter declarations                                  
c                                                                               
      integer :: span, nnode                                                    
      double precision ::                                                       
     &  b(mxvl,mxedof,*), gama(mxvl,ndim,*), nxi(*), neta(*), nzeta(*)          
c                                                                               
c                       locals                                                  
c                                                                               
      integer :: j, i, bpos1, bpos2                                             
      double precision, allocatable :: btemp(:,:,:)                             
c                                                                               
      allocate( btemp(mxvl,mxndel,ndim) )                                       
c                                                                               
c                       compute building blocks of b.                           
c                       btemp - j,1 = NX for node j at gpn                      
c                       btemp - j,2 = NY for node j at gpn                      
c                       btemp - j,3 = NZ for node j at gpn                      
c                                                                               
      do j = 1, nnode                                                           
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
         do i = 1, span                                                         
            btemp(i,j,1)= gama(i,1,1)*nxi(j)+gama(i,1,2)*neta(j)+               
     &                    gama(i,1,3)*nzeta(j)                                  
            btemp(i,j,2)= gama(i,2,1)*nxi(j)+gama(i,2,2)*neta(j)+               
     &                    gama(i,2,3)*nzeta(j)                                  
            btemp(i,j,3)= gama(i,3,1)*nxi(j)+gama(i,3,2)*neta(j)+               
     &                    gama(i,3,3)*nzeta(j)                                  
         end do                                                                 
      end do                                                                    
c                                                                               
c                       set position parameters                                 
c                                                                               
      bpos1 = nnode                                                             
      bpos2 = 2 * nnode                                                         
c                                                                               
c                       compute the linear strain-                              
c                       displacement matrices, using btemp.                     
c                                                                               
      do j = 1, nnode                                                           
!DIR$ VECTOR ALIGNED                                                            
            b(1:span,j,1)=       btemp(1:span,j,1)                              
!DIR$ VECTOR ALIGNED                                                            
            b(1:span,j,4)=       btemp(1:span,j,2)                              
!DIR$ VECTOR ALIGNED                                                            
            b(1:span,j,6)=       btemp(1:span,j,3)                              
c                                                                               
!DIR$ VECTOR ALIGNED                                                            
            b(1:span,bpos1+j,2)= btemp(1:span,j,2)                              
!DIR$ VECTOR ALIGNED                                                            
            b(1:span,bpos1+j,4)= btemp(1:span,j,1)                              
!DIR$ VECTOR ALIGNED                                                            
            b(1:span,bpos1+j,5)= btemp(1:span,j,3)                              
c                                                                               
!DIR$ VECTOR ALIGNED                                                            
            b(1:span,bpos2+j,3)= btemp(1:span,j,3)                              
!DIR$ VECTOR ALIGNED                                                            
            b(1:span,bpos2+j,5)= btemp(1:span,j,2)                              
!DIR$ VECTOR ALIGNED                                                            
            b(1:span,bpos2+j,6)= btemp(1:span,j,1)                              
      end do                                                                    
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine blcmp_cohes                  *          
c     *                                                              *          
c     *                       written by : aroy                      *          
c     *                                                              *          
c     *                   last modified : 04/21/2016 rhd             *          
c     *                                                              *          
c     *     this subroutine computes the B (=RLN) matrices           *          
c     *     at a given gauss point for a block of similar,           *          
c     *     non-conflicting interface elements.                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine blcmp_cohes( span, b, rot, shape, etype, nnode )               
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c                       parameter declarations                                  
c                                                                               
      integer :: span, etype, nnode                                             
      double precision ::                                                       
     &  b(mxvl,mxedof,*), rot(mxvl,ndim,*), shape(*), sh(mxndel)                
c                                                                               
      integer :: i, j, bpos1, bpos2                                             
      double precision :: zero                                                  
      data zero  / 0.0d00 /                                                     
c                                                                               
c            compute sh = L*N (N -> shape fn. array)                            
c                                                                               
c      do i = 1, nnode/2                                                        
c       sh(i) = -shape(i)                                                       
c      end do                                                                   
c      do i = nnode/2+1, nnode                                                  
c       sh(i) = shape(i)                                                        
c      end do                                                                   
c                                                                               
c           compute B = R*sh                                                    
c                                                                               
c                                                                               
      if( nnode .eq. 12 ) then                                                  
!DIR$ IVDEP                                                                     
       sh(1:6) = -shape(1:6)                                                    
!DIR$ IVDEP                                                                     
       sh(7:12) = shape(7:12)                                                   
       do j = 1, 12                                                             
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
           do i = 1, span                                                       
              b(i,j,1) = sh(j)*rot(i,1,1)                                       
              b(i,j,2) = sh(j)*rot(i,2,1)                                       
              b(i,j,3) = sh(j)*rot(i,3,1)                                       
c                                                                               
              b(i,12+j,1) = sh(j)*rot(i,1,2)                                    
              b(i,12+j,2) = sh(j)*rot(i,2,2)                                    
              b(i,12+j,3) = sh(j)*rot(i,3,2)                                    
c                                                                               
              b(i,24+j,1) = sh(j)*rot(i,1,3)                                    
              b(i,24+j,2) = sh(j)*rot(i,2,3)                                    
              b(i,24+j,3) = sh(j)*rot(i,3,3)                                    
c                                                                               
           end do                                                               
         end do                                                                 
        return                                                                  
      end if                                                                    
c                                                                               
      if( nnode .eq. 6 ) then                                                   
        sh(1:3) = -shape(1:3)                                                   
        sh(4:6) = shape(4:6)                                                    
        do j = 1, 6                                                             
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
           do i = 1, span                                                       
              b(i,j,1) = sh(j)*rot(i,1,1)                                       
              b(i,j,2) = sh(j)*rot(i,2,1)                                       
              b(i,j,3) = sh(j)*rot(i,3,1)                                       
c                                                                               
              b(i,6+j,1) = sh(j)*rot(i,1,2)                                     
              b(i,6+j,2) = sh(j)*rot(i,2,2)                                     
              b(i,6+j,3) = sh(j)*rot(i,3,2)                                     
c                                                                               
              b(i,12+j,1) = sh(j)*rot(i,1,3)                                    
              b(i,12+j,2) = sh(j)*rot(i,2,3)                                    
              b(i,12+j,3) = sh(j)*rot(i,3,3)                                    
c                                                                               
           end do                                                               
         end do                                                                 
        return                                                                  
      end if                                                                    
c                                                                               
      if( nnode .eq. 8 ) then                                                   
        sh(1:4) = -shape(1:4)                                                   
        sh(5:8) = shape(5:8)                                                    
        do j = 1, 8                                                             
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
           do i = 1, span                                                       
              b(i,j,1) = sh(j)*rot(i,1,1)                                       
              b(i,j,2) = sh(j)*rot(i,2,1)                                       
              b(i,j,3) = sh(j)*rot(i,3,1)                                       
c                                                                               
              b(i,8+j,1) = sh(j)*rot(i,1,2)                                     
              b(i,8+j,2) = sh(j)*rot(i,2,2)                                     
              b(i,8+j,3) = sh(j)*rot(i,3,2)                                     
c                                                                               
              b(i,16+j,1) = sh(j)*rot(i,1,3)                                    
              b(i,16+j,2) = sh(j)*rot(i,2,3)                                    
              b(i,16+j,3) = sh(j)*rot(i,3,3)                                    
c                                                                               
           end do                                                               
         end do                                                                 
        return                                                                  
      end if                                                                    
c                                                                               
c              general number of nodes                                          
c                                                                               
      bpos1 = nnode                                                             
      bpos2 = 2 * nnode                                                         
c                                                                               
      do i = 1, nnode/2                                                         
       sh(i) = -shape(i)                                                        
      end do                                                                    
      do i = nnode/2+1, nnode                                                   
       sh(i) = shape(i)                                                         
      end do                                                                    
      do j=1,nnode                                                              
!DIR$ IVDEP                                                                     
!DIR$ VECTOR ALIGNED                                                            
           do i = 1, span                                                       
              b(i,j,1) = sh(j)*rot(i,1,1)                                       
              b(i,j,2) = sh(j)*rot(i,2,1)                                       
              b(i,j,3) = sh(j)*rot(i,3,1)                                       
c                                                                               
              b(i,bpos1+j,1) = sh(j)*rot(i,1,2)                                 
              b(i,bpos1+j,2) = sh(j)*rot(i,2,2)                                 
              b(i,bpos1+j,3) = sh(j)*rot(i,3,2)                                 
c                                                                               
              b(i,bpos2+j,1) = sh(j)*rot(i,1,3)                                 
              b(i,bpos2+j,2) = sh(j)*rot(i,2,3)                                 
              b(i,bpos2+j,3) = sh(j)*rot(i,3,3)                                 
c                                                                               
           end do                                                               
       end do                                                                   
c                                                                               
       return                                                                   
       end                                                                      
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine blcmp1_axisymm               *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 07/21/2016 rhd             *          
c     *                                                              *          
c     *     computes the linear strain-displacement matrices at a    *          
c     *     given gauss point for a block of axisymmetric elements   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine blcmp1_axisymm( span, b, gama, nxi, neta, nzeta,               
     &                       shape, ce, radius, etype, nnode )                  
      implicit none                                                             
      include 'param_def'                                                       
c                                                                               
c                       parameter declarations                                  
c                                                                               
      integer :: span, etype, nnode                                             
      double precision ::                                                       
     &  b(mxvl,mxedof,*), gama(mxvl,ndim,*), nxi(*), neta(*),                   
     &  nzeta(*), shape(*), ce(mxvl,*), radius(*)                               
c                                                                               
c                       locals                                                  
c                                                                               
      integer :: i, j, row, col, bpos1, bpos2                                   
      double precision ::                                                       
     &  btemp(mxvl,mxndel,ndim), zero ! on stack                                
      logical ::  local_debug, axisym                                           
      data zero, local_debug / 0.0d0, .false. /                                 
c                                                                               
c                  compute building blocks of b for axisymmetric                
c                  elements.                                                    
c                       btemp - j,1 = NX for node j at gpn                      
c                       btemp - j,2 = NY for node j at gpn                      
c                       btemp - j,3 = nj/ri for node j at gpn                   
c                  radius(i) = radius to the current Gauss point for each       
c                              element                                          
c                  shape(j)  = shape function evaluated at the current          
c                              gauss point                                      
c                  Note that btemp(i,j,3) is for the axisymmetric hoop          
c                  strain = nj/ri                                               
c                                                                               
c                  When computing the b building blocks for                     
c                  axisymmetric elements only the upper left 2x2 of             
c                  the Jacobian matrix is used.                                 
c                                                                               
        do j = 1, nnode                                                         
!DIR$ IVDEP                                                                     
          do i = 1, span                                                        
              btemp(i,j,1) = gama(i,1,1)*nxi(j)+gama(i,1,2)*neta(j)             
              btemp(i,j,2) = gama(i,2,1)*nxi(j)+gama(i,2,2)*neta(j)             
              btemp(i,j,3) = shape(j) / radius(i)                               
           end do                                                               
        end do                                                                  
c                                                                               
        bpos1 = nnode                                                           
        bpos2 = 2 * nnode                                                       
c                                                                               
c                       compute the linear strain-                              
c                       displacement [b] matrices, using btemp                  
c                       for axisymmetric elements. The third group is all       
c                       zeros since there is no z-dof; btemp(i,j,3) is used     
c                       only for the third row in the first group to            
c                       give the hoop strain for the axisymmetric               
c                       element; the bottom two rows are zeros for no           
c                       yz or zx shear strain.                                  
c                                                                               
        do  j = 1, nnode                                                        
!DIR$ IVDEP                                                                     
           do i = 1, span                                                       
c                                                                               
              b(i,j,1)=       btemp(i,j,1)                                      
              b(i,j,3)=       btemp(i,j,3)                                      
              b(i,j,4)=       btemp(i,j,2)                                      
c                                                                               
              b(i,bpos1+j,2)= btemp(i,j,2)                                      
              b(i,bpos1+j,4)= btemp(i,j,1)                                      
c                                                                               
           end do                                                               
        end do                                                                  
c                                                                               
      return                                                                    
c                                                                               
      end                                                                       
