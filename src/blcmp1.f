c     ****************************************************************
c     *                                                              *
c     *                      subroutine blcmp1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 04/18/2016 rhd             *
c     *                                                              *
c     *     this subroutine computes the linear strain-              *
c     *     displacement matrices at a given gauss point for a       *
c     *     block of similar, non-conflicting elements.              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine blcmp1_srt( span, b, gama, nxi, neta, nzeta, shape,
     &                       ce, radius, etype, nnode )
      implicit integer (a-z)
$add param_def
c
c                       parameter declarations
c       
#dbl      double precision
#sgl      real
     &  b(mxvl,mxedof,*), gama(mxvl,ndim,*), nxi(*), neta(*),
     &  nzeta(*), shape(*), ce(mxvl,*), radius(*)
c
c                       locally allocated arrays
c       
#dbl      double precision ::
#sgl      real ::
     &  btemp(mxvl,mxndel,ndim), zero ! on stack
      logical ::  local_debug, axisym
      data zero, local_debug / 0.0, .false. /
@!DIR$ ASSUME_ALIGNED b:64, gama:64, nxi:64, neta:64, nzeta:64  
@!DIR$ ASSUME_ALIGNED btemp:64  
c
c
      if( local_debug ) then
        write(*,500)
        do i = 1, span
          write(*,*)'  Inverse Jacobian, elem #',i
          do row = 1, ndim
            write(*,510) (gama(i,row,col),col=1,ndim)
          end do
          write(*,505) radius(i)
        end do
        write(*,*) '[n, nxi, neta, nzeta]'
        do i = 1, nnode
          write(*,520) shape(i), nxi(i), neta(i), nzeta(i)
        end do
      end if
c
c                  compute building blocks of b, then form the
c                  full b matrix.  Choose between axisymmetric elements
c                  or volume elements.
c
      axisym = etype .eq. 10  .or.  etype .eq. 11
      if( axisym ) then
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
@!DIR$ LOOP COUNT MAX=###  
@!DIR$ IVDEP
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
@!DIR$ LOOP COUNT MAX=###  
@!DIR$ IVDEP
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
      end if
c
      if ( .not. axisym ) then
c
c                  Compute building blocks of b for volume elements.
c                       btemp - j,1 = NX for node j at gpn        
c                       btemp - j,2 = NY for node j at gpn        
c                       btemp - j,3 = NZ for node j at gpn        
c                     
        do  j = 1, nnode
@!DIR$ LOOP COUNT MAX=###  
@!DIR$ IVDEP
           do i = 1, span
              btemp(i,j,1)= gama(i,1,1)*nxi(j)+gama(i,1,2)*neta(j)+
     &                      gama(i,1,3)*nzeta(j)   
              btemp(i,j,2)= gama(i,2,1)*nxi(j)+gama(i,2,2)*neta(j)+
     &                      gama(i,2,3)*nzeta(j)   
              btemp(i,j,3)= gama(i,3,1)*nxi(j)+gama(i,3,2)*neta(j)+
     &                      gama(i,3,3)*nzeta(j)   
           end do
        end do
c
c                       set position parameters necessary for
c                       the computation of the tangent nonlinear 
c                       strain-displacement matrices.
c
        bpos1 = nnode
        bpos2 = 2 * nnode
c
c                       compute the linear strain-
c                       displacement matrices, using btemp.
c                       [b] was zeroed at creation in tanstf or
c                       drive_eps_sig_internal_forces. The zero
c                       locations in [b] never change
c       
        do  j = 1, nnode
@!DIR$ LOOP COUNT MAX=###  
@!DIR$ IVDEP
           do i = 1, span
c
              b(i,j,1)=       btemp(i,j,1)
              b(i,j,4)=       btemp(i,j,2)
              b(i,j,6)=       btemp(i,j,3)
c
              b(i,bpos1+j,2)= btemp(i,j,2)
              b(i,bpos1+j,4)= btemp(i,j,1)
              b(i,bpos1+j,5)= btemp(i,j,3)
c
              b(i,bpos2+j,3)= btemp(i,j,3)
              b(i,bpos2+j,5)= btemp(i,j,2)
              b(i,bpos2+j,6)= btemp(i,j,1) 
           end do
        end do

      end if
c
      if( local_debug ) then
        do i=1,span
          write(*,*)
          write(*,*)'   element #',i
          write(*,*)'   3 btemp arrays:'
          do row=1,3
            write(*,600)(btemp(i,col,row),col=1,nnode)
          end do
          write(*,*)'  [b] matrix:'
          do row=1,6
            write(*,600)(b(i,col,row),col=1,3*nnode)
          end do
        end do
      end if
c
      return
c
c
500   format(1x,//,'>>>>> blcmp1_srt',/,
     &               '     form the b matrix',/)
505   format(1x,'  Radius to Gauss point = r = ',f12.6)
510   format(1x,3f12.6)
520   format(1x,4f12.6)
600   format(1x,60f12.6)
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine blcmp1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 04/21/2016                 *
c     *                                                              *
c     *     this subroutine computes the linear strain-              *
c     *     displacement matrices at a given gauss point for a       *
c     *     block of similar, non-conflicting elements.              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine blcmp1( span, b, gama, nxi, neta, nzeta, nnode )
      implicit none
$add param_def
c
c                       parameter declarations
c 
      integer :: span, nnode      
#dbl      double precision ::
#sgl      real ::
     &  b(mxvl,mxedof,*), gama(mxvl,ndim,*), nxi(*), neta(*), nzeta(*)
c
c                       locals
c      
      integer :: j, i, bpos1, bpos2 
#dbl      double precision ::
#sgl      real ::
     &  btemp(mxvl,mxndel,ndim), zero  ! on stack
      data zero / 0.0d0 /
@!DIR$ ASSUME_ALIGNED b:64, gama:64, nxi:64, neta:64, nzeta:64  
@!DIR$ ASSUME_ALIGNED btemp:64  
c
c                       compute building blocks of b. 
c                       btemp - j,1 = NX for node j at gpn        
c                       btemp - j,2 = NY for node j at gpn        
c                       btemp - j,3 = NZ for node j at gpn        
c                     
      do j = 1, nnode
@!DIR$ LOOP COUNT MAX=###  
@!DIR$ IVDEP
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
c                       set position parameters necessary for
c                       the computation of the tangent nonlinear 
c                       strain-displacement matrices.
c
      bpos1 = nnode
      bpos2 = 2 * nnode
c
c                       compute the linear strain-
c                       displacement matrices, using btemp.
c       
      do  j = 1, nnode
@!DIR$ LOOP COUNT MAX=###  
@!DIR$ IVDEP
         do i = 1, span
c
            b(i,j,1)=       btemp(i,j,1)
            b(i,j,4)=       btemp(i,j,2)
            b(i,j,6)=       btemp(i,j,3)
c
            b(i,bpos1+j,2)= btemp(i,j,2)
            b(i,bpos1+j,4)= btemp(i,j,1)
            b(i,bpos1+j,5)= btemp(i,j,3)
c
            b(i,bpos2+j,3)= btemp(i,j,3)
            b(i,bpos2+j,5)= btemp(i,j,2)
            b(i,bpos2+j,6)= btemp(i,j,1) 
         end do
      end do
c
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
$add param_def
c
c                       parameter declarations
c  
      integer :: span, etype, nnode     
#dbl      double precision ::
#sgl      real ::
     &  b(mxvl,mxedof,*), rot(mxvl,ndim,*), shape(*), sh(mxndel)
c
      integer :: i, j, bpos1, bpos2     
      double precision :: zero
      data zero  / 0.0d00 /
@!DIR$ ASSUME_ALIGNED b:64, sh:64, rot:64, shape:64
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
       sh(1:6) = -shape(1:6)
       sh(7:12) = shape(7:12)
        do j = 1, 12
@!DIR$ LOOP COUNT MAX=###  
@!DIR$ IVDEP
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
@!DIR$ LOOP COUNT MAX=###  
@!DIR$ IVDEP
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
@!DIR$ LOOP COUNT MAX=###  
@!DIR$ IVDEP
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
@!DIR$ LOOP COUNT MAX=###  
@!DIR$ IVDEP
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
