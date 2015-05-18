c     ****************************************************************
c     *                                                              *
c     *  driver for the sparse matrix re-ordering. this code         *
c     *  had f90 allocatable arrays. it calls lower-level, f77       *
c     *  to perform the actual work.                                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/28/95                   *
c     *                                                              *
c     ****************************************************************
c
      subroutine sparse_reorder( neq, ncoeff, diag, rhs, coefs, kpt,
     &                           kind, perm, xadj, adjncy, itype  )
#dbl      implicit double precision (a-h,o-z)
#sgl      implicit real (a-h,o-z)
c
      integer xadj, adjncy, perm, deg, qsize, qlink
#sgl      real xdjncy
#dbl      double precision xdjncy
c
c		these allocatables are released before returning.
c               this set of routines retains no space between calls.
c
      allocatable jtemp(:), maxan(:), itemp(:),
     &            xtemp(:), xdjncy(:), invp(:), deg(:),
     &            marker(:), qsize(:), qlink(:), btemp(:)
c
      dimension coefs(*), diag(*), rhs(*), kind(*), kpt(*), perm(*)
      data nb / 2 /
c
c
      allocate( jtemp(neq+1) )
      allocate( maxan(neq+1) )
c  
      call sparse_reorder2( neq, kpt, kind, jtemp, maxan,
     &                      ncof, maxn )
c
      incof        = max(ncof,neq)
#sgl      isize_xdjncy = ((ncoeff * nb) + neq) * 2
#dbl      isize_xdjncy = ((ncoeff * nb) + neq)
      allocate( itemp(incof) )
      allocate( xtemp(ncof) )
      allocate( xdjncy(isize_xdjncy) )
c
c               perform minimum degree reordering if required. skip
c               if re-solving a matrix with identical sparisty.
c
      if ( itype .eq. 1 ) then 
        call sparse_reorder3( neq, kpt, kind, itemp, jtemp,
     &                        maxn, nt, xadj, adjncy, xtemp, coefs,
     &                        xdjncy, jmax, maxan )
        nterms = nt
        maxcon = maxn
c
c		minimum degree reordering procedure
c
        allocate( invp(neq+1)   )
        allocate( deg(neq+1)    )
        allocate( marker(neq+1) )
        allocate( qsize(neq+1)  )
        allocate( qlink(neq+1)  )
        call sparse_reorder5( neq, xadj, adjncy, perm, invp, deg,
     &                        marker, qsize, qlink )
        deallocate( invp   )
        deallocate( deg    )
        deallocate( marker )
        deallocate( qsize  )
        deallocate( qlink  )
      end if
c
c                recompute adjncy, re-order non-zero coefficients.
c
      call sparse_reorder3( neq, kpt, kind, itemp, jtemp,
     &                      maxn, nt, xadj, adjncy, xtemp, coefs,
     &                      xdjncy, jmax, maxan )
      deallocate( jtemp )
      deallocate( maxan )
c
      allocate( btemp(neq) )
      call sparse_reorder4( neq, xadj, adjncy, perm,
     &                      diag, rhs, coefs, kpt, kind,
     &                      xdjncy, itemp, btemp )
      deallocate( btemp )
      deallocate( itemp )
      deallocate( xtemp )
      deallocate( xdjncy )
c
      return
      end
