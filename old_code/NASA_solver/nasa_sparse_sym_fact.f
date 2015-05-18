c                 sparse_sym_fact
c
c     purpose - this routine performs symbolic factorization
c        on a permuted linear system and it also sets up the
c        compressed data structure for the system.
c
c     input parameters -
c        neqns  - number of equations.
c        (xadj,adjncy) - the adjacency structure.
c        (perm,invp) - the permutation vector and its inverse.
c
c     updated parameter -
c        maxsub - size of the subscript array nzsub.  on return,
c                 it contains the number of subscripts used
c
c     output parameters -
c        xlnz   - index into the nonzero storage vector lnz.
c        (xnzsub,nzsub) - the compressed subscript vectors.
c        maxlnz - the number of nonzeros found.
c        flag   - error flag.  positive value indicates that
c                 nzsub array is too small.
c
c     working parameters -
c        mrglnk - a vector of size neqns.  at the kth step,
c                 mrglnk(k), mrglnk(mrglnk(k)) , ..., is a list
c                 containing all those columns l(*,j) with j
c                 less than k, such that its first off-diagonal
c                 nonzero is l(k,j).  thus, the nonzero structure
c                 of column l(*,k) can be found by merging that
c                 of such columns l(*,j) with the structure of
c                 a(*,k).
c        rchlnk - a vector of size neqns.  it is used to
c                 accumulate the structure of each column l(*,k).
c                 at the end of the kth step, rchlnk(k),
c                 rchlnk(rchlnk(k)), ..., is the list of
c                 positions of nonzeros in column k of the
c                 factor l.
c        marker - an integer vector of length neqns.  it is used
c                 to test if mass symbolic elimination can be
c                 performed.  that is, it is used to check
c                 whether the structure of the current column k
c                 being processed is completely determined by
c                 the single column mrglnk(k).
c
c****************************************************************
c
c
      subroutine sparse_sym_fact( neqns, xadj, adjncy, xlnz, maxlnz,
     &                            xnzsub, nzsub, maxsub, flag )
c
c****************************************************************
c
      implicit integer (a-z)
      dimension xadj(*), adjncy(*), xlnz(*), xnzsub(*), nzsub(*)
      allocatable perm(:), invp(:), rchlnk(:), mrglnk(:),
     &            marker(:)
c
c****************************************************************
c
c
c        ------------------
c        initialization ...
c        ------------------
c
      allocate( perm(neqns+1) )
      allocate( invp(neqns+1) )
      allocate( rchlnk(neqns+1) )
      allocate( mrglnk(neqns+1) )
      allocate( marker(neqns+1) )
c
      nzbeg = 1
      nzend = 0
      do k = 1, neqns
        perm(k)   = k
        invp(k)   = k
        mrglnk(k) = 0
        marker(k) = 0
        xlnz(k)   = 0
      end do
      xlnz(1) = 1
c
c        ----------------------------------------------
c        for each column ... knz counts the number
c        of nonzeros in column k accumulated in rchlnk.
c        ----------------------------------------------
c
         np1 = neqns + 1
         do  1700  k = 1, neqns
             knz = 0
             mrgk = mrglnk(k)
             mrkflg = 0
             marker(k) = k
             if  ( mrgk .ne. 0 )  marker(k) = marker(mrgk)
             xnzsub(k) = nzend
             node = perm(k)
             jstrt = xadj(node)
             jstop = xadj(node+1) - 1
             if  ( jstrt .gt. jstop )  go to 1600
c                -----------------------------------
c                use rchlnk to link through the
c                structure of a(*,k) below diagonal.
c                -----------------------------------
                 rchlnk(k) = np1
                 do  300  j = jstrt, jstop
                     nabor = adjncy(j)
                     nabor = invp(nabor)
                     if  ( nabor .le. k )  go to 300
                         rchm = k
  200                    continue
                             m = rchm
                             rchm = rchlnk(m)
                             if  ( rchm .le. nabor )  go to 200
                         knz = knz+1
                         rchlnk(m) = nabor
                         rchlnk(nabor) = rchm
                         if  ( marker(nabor) .ne. marker(k) )
     1                       mrkflg = 1
  300            continue
c                --------------------------------------
c                test for mass symbolic elimination ...
c                --------------------------------------
                 lmax = 0
                 if  ( mrkflg .ne. 0  .or.  mrgk .eq. 0 )  go to 400
                 if  ( mrglnk(mrgk) .ne. 0 )  go to 400
                     xnzsub(k) = xnzsub(mrgk) + 1
                     knz = xlnz(mrgk+1) - (xlnz(mrgk) + 1)
                     go to 1500
  400            continue
c                -----------------------------------------------
c                link through each column i that affects l(*,k).
c                -----------------------------------------------
                 i = k
  500            continue
                     i = mrglnk(i)
                     if  ( i .eq. 0 )  go to 900
                         inz = xlnz(i+1) - (xlnz(i) + 1)
                         jstrt = xnzsub(i) +  1
                         jstop = xnzsub(i) + inz
                         if  ( inz .le. lmax )  go to 600
                             lmax = inz
                             xnzsub(k) = jstrt
  600                    continue
c                        ----------------------------------
c                        merge structure of l(*,i) in nzsub
c                        into rchlnk.
c                        ----------------------------------
                         rchm = k
                         do  800  j = jstrt, jstop
                             nabor = nzsub(j)
  700                        continue
                                 m = rchm
                                 rchm = rchlnk(m)
                                 if  ( rchm .lt. nabor )  go to 700
                             if  ( rchm .eq. nabor )  go to 800
                                 knz = knz+1
                                 rchlnk(m) = nabor
                                 rchlnk(nabor) = rchm
                                 rchm = nabor
  800                    continue
                         go to 500
  900            continue
c                -----------------------------------
c                check if subscripts duplicate those
c                of another column.
c                -----------------------------------
                 if  ( knz .eq. lmax )  go to 1500
c                    -----------------------------------------------
c                    or if tail of k-1st column matches head of kth.
c                    -----------------------------------------------
                     if  ( nzbeg .gt. nzend )  go to 1300
                         i = rchlnk(k)
                         do  1000  jstrt = nzbeg, nzend
                             if  ( nzsub(jstrt) - i )  1000, 1100, 1300
 1000                    continue
                         go to 1300
 1100                    continue
                         xnzsub(k) = jstrt
                         do  1200  j = jstrt, nzend
                             if  ( nzsub(j) .ne. i )  go to 1300
                                 i = rchlnk(i)
                                 if  ( i .gt. neqns )  go to 1500
 1200                    continue
                         nzend = jstrt - 1
 1300                continue
c                    ----------------------------------------
c                    copy the structure of l(*,k) from rchlnk
c                    to the data structure (xnzsub, nzsub).
c                    ----------------------------------------
                     nzbeg = nzend +  1
                     nzend = nzend + knz
                     if  ( nzend .gt. maxsub )  go to 1800
                         i = k
                         do  1400  j = nzbeg, nzend
                             i = rchlnk(i)
                             nzsub(j) = i
                             marker(i) = k
 1400                    continue
                         xnzsub(k) = nzbeg
                         marker(k) = k
 1500            continue
c                ---------------------------------------------
c                update the vector mrglnk.  note column l(*,k)
c                just found is required to determine column
c                l(*,j), where l(j,k) is the first nonzero in
c                l(*,k) below diagonal.
c                ---------------------------------------------
                 if  ( knz .le. 1 )  go to 1600
                     kxsub = xnzsub(k)
                     i = nzsub(kxsub)
                     mrglnk(k) = mrglnk(i)
                     mrglnk(i) = k
 1600        continue
             xlnz(k+1) = xlnz(k) + knz
 1700    continue
c
c              flag = 0 means sufficient memory available for fill.
c              maxlnz = number of non-zero terms in factored matrix.
c
         maxlnz = xlnz(neqns) - 1
         maxsub = xnzsub(neqns)
         xnzsub(neqns+1) = xnzsub(neqns)
         flag = 0
         deallocate( perm )
         deallocate( invp )
         deallocate( rchlnk )
         deallocate( mrglnk )
         deallocate( marker )
         return
c
 1800    continue
c        ----------------------------------------------------
c        error - insufficient storage for nonzero subscripts.
c                maxlnz not defined on exit
c        ----------------------------------------------------
         flag = 1
         deallocate( perm )
         deallocate( invp )
         deallocate( rchlnk )
         deallocate( mrglnk )
         deallocate( marker )
         return
c
      end
