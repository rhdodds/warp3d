c     ****************************************************************
c     *                                                              *
c     *                      subroutine gtlsn1                       *
c     *                                                              *
c     *             -- incremental strain for solid elements --      *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 4/30/2016 rhd              *
c     *                                                              *
c     *      processes a material (gauss) point for a block of       *
c     *      identical type solid elements. compute the [B] matrix   *
c     *      at the point for each element in the block. multiply    *
c     *      [B] into incremental displacement vector for element    *
c     *      nodes: (n+1) - n to define a strain increment for       *
c     *      n->n+1. used for both small and finite strain           *
c     *      theory. for finite strains, the [B] is evaluated using  *
c     *      mid-step (n+1/2) element configuration. the strain      *
c     *      increment is then most often noted D (or the rate D     *
c     *      * dt).                                                  *
c     *                                                              *
c     ****************************************************************
c
c           
      subroutine gtlsn1( span, nnode,
     &                   due, deps, gama, nxi, neta,
     &                   nzeta, vol_block, bbar, eps_bbar, b )
      implicit integer (a-z)
      include 'param_def'
c
c                      parameter declarations
c
      double precision ::
     & due(mxvl,*), deps(mxvl,nstr), gama(mxvl,3,3),
     & nxi(*), neta(*), nzeta(*), vol_block(mxvl,8,*), eps_bbar,
     & b(mxvl,mxedof,*)
      logical :: bbar
c
c                      locals
c
c!DIR$ ASSUME_ALIGNED due:64, deps:64, gama:64, nxi:64, neta:64
c!DIR$ ASSUME_ALIGNED nzeta:64, vol_block:64, b:64  
c      
c                       compute linear strain-displacement
c                       [B] matrix for this material (gauss) point at
c                       all elements in the block. modify
c                       for B-bar as needed.

      call blcmp1( span, b, gama, nxi, neta, nzeta, nnode )  
      if ( bbar ) call bmod( b, vol_block, span, mxvl, eps_bbar,
     &                       mxedof )
c
c                       multiply [B] x displacement increment. this is
c                       done for all elements in the block for this
c                       material (gauss) point nuber. take advantage of
c                       sparsity in [B] during multiply.
c   
      deps = 0.0d00
c      
      bpos1 = nnode
      bpos2 = 2*nnode 
c     
      do j = 1, nnode
!DIR$ LOOP COUNT MAX=128  
!DIR$ IVDEP
         do i = 1,span
           deps(i,1) = deps(i,1) +  b(i,j,1) * due(i,j) +
     &                              b(i,bpos1+j,1) * due(i,bpos1+j) +
     &                              b(i,bpos2+j,1) * due(i,bpos2+j)
           deps(i,2) = deps(i,2) +  b(i,j,2) * due (i,j) +
     &                              b(i,bpos1+j,2) * due(i,bpos1+j)+
     &                              b(i,bpos2+j,2) * due(i,bpos2+j)
           deps(i,3) = deps(i,3) +  b(i,j,3) * due (i,j) +
     &                              b(i,bpos1+j,3) * due(i,bpos1+j)+
     &                              b(i,bpos2+j,3) * due(i,bpos2+j)
           deps(i,4) = deps(i,4) +  b(i,j,4)*due(i,j)+      
     &                              b(i,bpos1+j,4)*due(i,bpos1+j)
           deps(i,5) = deps(i,5) +  b(i,bpos1+j,5)*due(i,bpos1+j)+      
     &                              b(i,bpos2+j,5)*due(i,bpos2+j)
           deps(i,6) = deps(i,6) +  b(i,j,6)*due(i,j)+      
     &                              b(i,bpos2+j,6)*due(i,bpos2+j)
         end do
      end do
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine gtlsn2                       *
c     *                                                              *
c     *         -- incremental "strain" for cohesive elements --     *
c     *                                                              *
c     *                       written by : aroy                      *
c     *                                                              *
c     *                   last modified : 06/12/99                   *
c     *                                                              *
c     *    this subroutine computes the incremental linear           *
c     *    "strain" of a given gauss point for use in the stress     *
c     *    recovery routine for an element in a block of similar,    *
c     *    non-conflicting cohesive elements excluding geometric     *
c     *    nonlinearity. the accumulated linear strains are updated  *
c     *    "strains" are just displacement jumps across the          *
c     *    interface for cohesive elements                           *
c     *                                                              *
c     ****************************************************************
c
c           
      subroutine gtlsn2( span, nnode, due, dgstrn, dgstrs, rot, shape, 
     &                   etype, gpn, felem, iout )
      implicit integer (a-z)
      include 'param_def'
c
c                      parameter declarations
c
      double precision ::
     & due(mxvl,*), dgstrn(mxvl,nstr), dgstrs(mxvl,*), rot(mxvl,3,3),
     & shape(*), b(mxvl,mxedof,nstr)
c
c                      locals
c
      double precision, parameter :: zero = 0.0d00
      logical :: local_debug
      data local_debug / .false. /
c!DIR$ ASSUME_ALIGNED due:64, dgstrn:64, dgstrs:64, rot:64
c!DIR$ ASSUME_ALIGNED shape:64, b:64  
c      
c         compute linear relative displacement jump
c         [b] matrix for this gauss point at all elements in 
c         the block. 
c
      call blcmp_cohes( span, b, rot, shape, etype,  nnode )  
c   
      dgstrn = zero
c      
      bpos1 = nnode
      bpos2 = 2*nnode
c     
      do j = 1, nnode
!DIR$ LOOP COUNT MAX=128  
!DIR$ IVDEP
         do i = 1, span
           dgstrn(i,1)= dgstrn(i,1)+b(i,j,1) * due(i,j) +
     &                              b(i,bpos1+j,1) * due(i,bpos1+j) +
     &                              b(i,bpos2+j,1) * due(i,bpos2+j)
           dgstrn(i,2)= dgstrn(i,2) + b(i,j,2) * due (i,j) +
     &                                b(i,bpos1+j,2) * due(i,bpos1+j)+
     &                                b(i,bpos2+j,2) * due(i,bpos2+j)
           dgstrn(i,3)= dgstrn(i,3)+ b(i,j,3) * due (i,j) +
     &                               b(i,bpos1+j,3) * due(i,bpos1+j)+
     &                               b(i,bpos2+j,3) * due(i,bpos2+j)
         end do
      end do
c
c                       update the accumulated displacement jumps.
c                       dgstrs is total "strain" at end of step.
c                       dgstrn is total "strain" increment over step.
c       
!DIR$ LOOP COUNT MAX=128  
!DIR$ IVDEP
      do i = 1, span
         dgstrs(i,1) = dgstrs(i,1) + dgstrn(i,1)
         dgstrs(i,2) = dgstrs(i,2) + dgstrn(i,2)
         dgstrs(i,3) = dgstrs(i,3) + dgstrn(i,3)
         dgstrs(i,4) = zero
         dgstrs(i,5) = zero
         dgstrs(i,6) = zero
      end do
c
      if ( local_debug .and. gpn .eq. 1 ) then
        write(iout,9000) 
        do i = 1, span
          write(iout,9100) felem + i - 1
          do j = 1, nnode
            write(iout,9200) j, due(i,j), due(i,bpos1+j), due(i,bpos2+j)
          end do
        end do
      end if
c
      return
c
 9000 format('>>>> Debug in gtlsn2 (cohesive element)...')
 9100 format('      Element: ',i10)
 9200 format(15x,i3,3f20.10)
      end

