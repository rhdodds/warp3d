c     ****************************************************************
c     *                                                              *
c     *  serice routines for the sparse matrix re-ordering.          *
c     *  this code is straightforward f77                            *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/28/95                   *
c     *                                                              *
c     ****************************************************************
c
      subroutine sparse_reorder4( neqns, xadj, adjncy, perm, maxcon,
     &                            nterms, ncoeff, d, b, a, irow, icoln,
     &                            xdjncy, itemp, btemp )
      return
      end

      subroutine sparse_reorder3( neq, irow, icoln, ncoeff, itemp,
     &                            jtemp, maxn, nt, xadj, adjncy,
     &                            xtemp, a, xdjncy, jmax, maxan )
      return
      end


      subroutine sparse_reorder2( neq, irow, icoln, ncoeff, jtemp,
     &                            maxan, ncof, maxn )
      return
      end

      subroutine  sparse_reorder5(
     &  neqns, xadj, adjncy, perm, invp, dhead, qsize, llist,
     &  marker, nofsub, maxcon, nterms )
      return
      end

      subroutine  sparse_reorder6(
     &  ehead, neqns, xadj, adjncy, delta, mdeg, dhead, dforw,
     &  dbakw, qsize, llist, marker, maxint, tag )
      return
      end

      subroutine sparse_reorder7(
     &  mdnode, xadj, adjncy, dhead, dforw, dbakw, qsize, llist,
     &  marker, maxint, tag )
      return
      end

      subroutine sparse_reorder8(
     &  neqns, xadj, adjncy, dhead, dforw, dbakw, qsize, llist,
     &  marker )
      return
      end

      subroutine sparse_reorder9(
     &  neqns, perm, invp, qsize )
      return
      end
