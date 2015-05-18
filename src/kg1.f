c     ****************************************************************
c     *                                                              *
c     *                      subroutine kg1                          *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 07/02/91                   *
c     *                                                              *
c     *     this subroutine computes the geometric tangent stiff-    *
c     *     nesses at a gauss point in uniform global coordinates    *
c     *     at state (n+1) for a block of similar, non-conflicting   *
c     *     q3disop or l3disop elements.                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine kg1( span, cp, icp, gama, nxi, neta, nzeta, nnode,
     &                sig, dj, w, ek, vol, bbar, nsz, totdof )
      use main_data, only: asymmetric_assembly
      implicit integer (a-z)
$add param_def
c
c                       parameter declarations
c
#dbl      double precision
#sgl      real
     &    gama(mxvl,ndim,*), nxi(*), neta(*), nzeta(*), sig(mxvl,*),
     &    dj(*), ek(nsz,*), w, vol(mxvl,8,*)
      integer cp(*), icp(mxutsz,*)
      logical bbar
c
c                       locally allocated arrays - on stack
c
#dbl      double precision
#sgl      real
     &    gtg(mxvl,mxnusz), gxi(mxvl,mxndel),
     &    geta(mxvl,mxndel), gzeta(mxvl,mxndel)
c
c                       calculate the arrays that are the building 
c                       blocks of the geometric stiffness.
c    
      if ( bbar ) then                   
        do j = 1, nnode
          do i = 1, span
           gxi(i,j)   = vol(i,j,1)
           geta(i,j)  = vol(i,j,2)
           gzeta(i,j) = vol(i,j,3)
          end do
        end do
      else
        do j = 1, nnode
          do i = 1, span
            gxi(i,j)  =  gama(i,1,1)*nxi(j)+gama(i,1,2)*neta(j)+
     &                   gama(i,1,3)*nzeta(j)
            geta(i,j) =  gama(i,2,1)*nxi(j)+gama(i,2,2)*neta(j)+
     &                   gama(i,2,3)*nzeta(j)
            gzeta(i,j) = gama(i,3,1)*nxi(j)+gama(i,3,2)*neta(j)+
     &                   gama(i,3,3)*nzeta(j)
          end do
        end do
      end if
c
c                       compute the geometric tangent stiffnesses for
c                       one (x,y,z) direction.
c
      do j = 1, cp(nnode)+nnode
         do i = 1, span
            gtg(i,j)= (gxi(i,icp(j,1))*gxi(i,icp(j,2))*sig(i,1)+
     &                 geta(i,icp(j,1))*geta(i,icp(j,2))*sig(i,2)+
     &                 gzeta(i,icp(j,1))*gzeta(i,icp(j,2))*sig(i,3)+
     &                (gxi(i,icp(j,1))*geta(i,icp(j,2))+
     &                 gxi(i,icp(j,2))*geta(i,icp(j,1)))*sig(i,4)+
     &                (geta(i,icp(j,1))*gzeta(i,icp(j,2))+
     &                 geta(i,icp(j,2))*gzeta(i,icp(j,1)))*sig(i,5)+
     &                (gxi(i,icp(j,1))*gzeta(i,icp(j,2))+
     &                 gxi(i,icp(j,2))*gzeta(i,icp(j,1)))*sig(i,6))*
     &                 dj(i)*w
         end do
      end do
c
c                       incorporate the 1d geometric tangent stiffnesses
c                       into the element stiffnesses.
c
      if (.not. asymmetric_assembly) then ! symmetric
      do l = 1, nnode
c
         cp1 = cp(l)
         cp2 = cp(nnode+l)+nnode
         cp3 = cp(2*nnode+l)+2*nnode
c
         do j = 1, l
            do  i = 1, span
               ek(cp1+j,i) = ek(cp1+j,i) + gtg(i,cp1+j)
               ek(cp2+j,i) = ek(cp2+j,i) + gtg(i,cp1+j)
               ek(cp3+j,i) = ek(cp3+j,i) + gtg(i,cp1+j)
            end do
         end do
      end do

      else ! asymmetric, need to do both sides of matrix
      do l = 1, nnode
c
         cp1 = cp(l)
         cp2 = cp(nnode+l)+nnode
         cp3 = cp(2*nnode+l)+2*nnode
c
         do j = 1, l
            do  i = 1, span
               ! col l, row j
               r = j
               c = l
               k = (c-1)*totdof + r
               ek(k,i) = ek(k,i) + gtg(i,cp1+j)
               if (r .ne. c) then
                 k = (r-1)*totdof + c
                 ek(k,i) = ek(k,i) + gtg(i,cp1+j)
               end if

               ! col nnode+l, row nnode+j
               r = nnode + j
               c = nnode + l
               k = (c-1)*totdof + r
               ek(k,i) = ek(k,i) + gtg(i,cp1+j)
               if (r .ne. c) then
                 k = (r-1)*totdof + c
                 ek(k,i) = ek(k,i) + gtg(i,cp1+j)
               end if

               ! col 2*nnode+l, row 2*nnode + j
               r = 2*nnode + j
               c = 2*nnode + l
               k = (c-1)*totdof + r
               ek(k,i) = ek(k,i) + gtg(i,cp1+j)
               if (r .ne. c) then
                 k = (r-1)*totdof + c
                 ek(k,i) = ek(k,i) + gtg(i,cp1+j)
               end if

            end do
         end do
      end do

      end if

      return
      end
