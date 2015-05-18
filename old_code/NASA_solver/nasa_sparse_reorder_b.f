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
      subroutine sparse_reorder4( neqns, xadj, adjncy, perm,
     &                            d, b, a, irow, icoln,
     &                            xdjncy, itemp, btemp )
c
#dbl      implicit double precision (a-h,o-z)
#sgl      implicit real (a-h,o-z)
c
      dimension adjncy(*), perm(*), itemp(*), xadj(*), icoln(*),
     &          irow(*), a(*), b(*), d(*), xdjncy(*),  btemp(*)
      integer adjncy, perm, xadj
      data ibig / 1000000 /
c
      neq = neqns
      do ii = 1, neq
        itemp(perm(ii)) = ii
      end do
c
      lcont = 1
      do i = 1, neq
        jq = perm(i)
        is = xadj(perm(i))
        iend = xadj(jq+1)
        ncont = 0
           do j = is, iend-1
                if( itemp(adjncy(j)) .gt. i ) then
                  icoln(lcont) = itemp(adjncy(j))
                  a(lcont) = xdjncy(j)
                  lcont = lcont + 1
                  ncont = ncont + 1
                end if
           end do
         irow(i) = ncont
      end do
c
c           reorder the matrix with small to large
c
      mycont = 1
      jcont = 1
      do i = 1, neq
        jmin = ibig
        do k = 1, irow(i)
         icont = jcont
         jmin = ibig
         do j = 1, irow(i)
            if( icoln(icont) .lt. jmin ) then     
              jmin = icoln(icont) 
              jj = icont
              ak11 = a(icont)
            end if
            icont = icont + 1
         end do
         xdjncy(mycont) = ak11
         adjncy(mycont) = jmin
         mycont = mycont + 1
         icoln(jj) = ibig
        end do
       jcont = jcont + irow(i)
      end do
c    
      do i = 1 , lcont - 1
         a(i) = xdjncy(i)
         icoln(i) = adjncy(i)
      end do
c
      do ii = 1, neq
        btemp(ii) = b(perm(ii))
      end do
c
      do ii = 1, neq
        b(ii) = btemp(ii)
      end do
c
      do ii = 1, neq
        btemp(ii) = d(perm(ii))
      end do
c
      do ii = 1, neq
        d(ii) = btemp(ii)
      end do
c
      return
      end

      subroutine sparse_reorder3( neq, irow, icoln, itemp,
     &                            jtemp, maxn, nt, xadj, adjncy,
     &                            xtemp, a, xdjncy, jmax, maxan )
#dbl      implicit double precision (a-h,o-z)
#sgl      implicit real (a-h,o-z)
      dimension xtemp(*), a(*), xdjncy(*), irow(*), itemp(*),
     &          jtemp(*), maxan(*), icoln(*)
      integer adjncy(*), xadj(*)
c 
      do ii = 1, neq
        jtemp(ii) = 0
      end do
      icont = 1
      ncont = 1
      lcont = 1
      jmax  = 0
      xadj(1) = lcont
c
      do i = 1, neq
       kcont = 0
        do j = 1, jtemp(i)
          adjncy(ncont) = itemp(maxan(i)+j-1)
          xdjncy(ncont) = xtemp(maxan(i)+j-1)
          ncont = ncont + 1
          kcont = kcont + 1
        end do
        do j = 1, irow(i)
          iix = icoln(icont)
          xdjncy(ncont) = a(icont)
          adjncy(ncont) = icoln(icont)
          ncont = ncont + 1
          itemp(maxan(i)+jtemp(i)) = iix
          itemp(maxan(iix)+jtemp(iix)) = i
          xtemp(maxan(iix)+jtemp(iix)) = a(icont)
          if( jtemp(iix) .ge. maxn ) then
            write(*,*)' increase the dimension of maxn to ....',maxn
            call die_gracefully
            stop
          end if
          jtemp(iix) = jtemp(iix) + 1
          jtemp(i) = jtemp(i) + 1
          icont = icont + 1
          kcont = kcont + 1
        end do
        lcont = lcont + kcont
        xadj(i+1) = lcont
        jmax = max(jmax,jtemp(i))
      end do
c
      nt = lcont - 1
      xadj(neq+1) = lcont
c
      return
      end


      subroutine sparse_reorder2( neq, irow, icoln, jtemp,
     &                            maxan, ncof, maxn )
#dbl      implicit double precision (a-h,o-z)
#sgl      implicit real (a-h,o-z)
      integer neq, irow(*), jtemp(*), maxan(*), icoln(*)
c 
c              nt number of terms in adjncy
c
      do ii = 1, neq
        jtemp(ii) = 0
      end do
c
      icont = 1
      do i = 1 , neq
        do j = 1, irow(i)
           iix        = icoln(icont)
           jtemp(iix) = jtemp(iix) + 1
           jtemp(i)   = jtemp(i) + 1
           icont      = icont + 1
        end do
      end do
c
      jmax = 0
      maxan(1) = 1
      do i = 1, neq
       jmax       = max(jmax,jtemp(i))
       maxan(i+1) = maxan(i) + jtemp(i)
      end do
      maxn = jmax
      ncof = maxan(neq+1)
c
      return
      end

c     sparse_reorder5
c
c     purpose - this routine implements the minimum degree
c        algorithm.  it makes use of the implicit representation
c        of elimination graphs by quotient graphs, and the
c        notion of indistinguishable nodes.  it also implements
c        the modifications by multiple elimination and minimum
c        external degree.
c        ---------------------------------------------
c        caution - the adjacency vector adjncy will be
c        destroyed.
c        ---------------------------------------------
c
c     input parameters -
c        neqns  - number of equations.
c        (xadj,adjncy) - the adjacency structure.
c        delta  - tolerance value for multiple elimination.
c        maxint - maximum machine representable (short) integer
c                 (any smaller estimate will do) for marking
c                 nodes.
c
c     output parameters -
c        perm   - the minimum degree ordering.
c        invp   - the inverse of perm.
c        nofsub - an upper bound on the number of nonzero
c                 subscripts for the compressed storage scheme.
c
c     working parameters -
c        dhead  - vector for head of degree lists.
c        invp   - used temporarily for degree forward link.
c        perm   - used temporarily for degree backward link.
c        qsize  - vector for size of supernodes.
c        llist  - vector for temporary linked lists.
c        marker - a temporary marker vector.
c
c***************************************************************
c
c
      subroutine  sparse_reorder5(
     &  neqns, xadj, adjncy, perm, invp, dhead, marker,
     &  qsize, llist )
c***************************************************************
c
      integer    adjncy(1), dhead(1) , invp(1)  , llist(1) ,
     1           marker(1), perm(1)  , qsize(1)
      integer    xadj(1)
      integer    delta , ehead , i     , maxint, mdeg  ,
     1           mdlmt , mdnode, neqns , nextmd, num, tag
c
c***************************************************************
c
         maxint = 10000000
         delta = 0
         if  ( neqns .le. 0 )  return
c
c        ------------------------------------------------
c        initialization for the minimum degree algorithm.
c        ------------------------------------------------
c         nofsub = 0
         call sparse_reorder8( neqns, xadj, dhead, invp,
     &                         perm, qsize, llist, marker )
c
c        ----------------------------------------------
c        num counts the number of ordered nodes plus 1.
c        ----------------------------------------------
         num = 1
c
c        -----------------------------
c        eliminate all isolated nodes.
c        -----------------------------
         nextmd = dhead(1)
  100    continue
             if  ( nextmd .le. 0 )  go to 200
                 mdnode = nextmd
                 nextmd = invp(mdnode)
                 marker(mdnode) = maxint
                 invp(mdnode) = - num
                 num = num + 1
                 go to 100
c
  200    continue
c        ----------------------------------------
c        search for node of the minimum degree.
c        mdeg is the current minimum degree;
c        tag is used to facilitate marking nodes.
c        ----------------------------------------
         if  ( num .gt. neqns )  go to 1000
         tag = 1
         dhead(1) = 0
         mdeg = 2
  300    continue
             if  ( dhead(mdeg) .gt. 0 )  go to 400
                 mdeg = mdeg + 1
                 go to 300
  400        continue
c            -------------------------------------------------
c            use value of delta to set up mdlmt, which governs
c            when a degree update is to be performed.
c            -------------------------------------------------
             mdlmt = mdeg + delta
             ehead = 0
c
  500        continue
                 mdnode = dhead(mdeg)
                 if  ( mdnode .gt. 0 )  go to 600
                     mdeg = mdeg + 1
                     if  ( mdeg .gt. mdlmt )  go to 900
                         go to 500
  600            continue
c                ----------------------------------------
c                remove mdnode from the degree structure.
c                ----------------------------------------
                 nextmd = invp(mdnode)
                 dhead(mdeg) = nextmd
                 if  ( nextmd .gt. 0 )  perm(nextmd) = - mdeg
                 invp(mdnode) = - num
c                 nofsub = nofsub + mdeg + qsize(mdnode) - 2
                 if  ( num+qsize(mdnode) .gt. neqns )  go to 1000
c                ----------------------------------------------
c                eliminate mdnode and perform quotient graph
c                transformation.  reset tag value if necessary.
c                ----------------------------------------------
                 tag = tag + 1
                 if  ( tag .lt. maxint )  go to 800
                     tag = 1
                     do  700  i = 1, neqns
                         if  ( marker(i) .lt. maxint )  marker(i) = 0
  700                continue
  800            continue
                 call  sparse_reorder7(
     &               mdnode, xadj, adjncy, dhead, invp, perm, qsize,
     &               llist, marker, maxint, tag )
                 num = num + qsize(mdnode)
                 llist(mdnode) = ehead
                 ehead = mdnode
                 if  ( delta .ge. 0 )  go to 500
  900        continue
c            -------------------------------------------
c            update degrees of the nodes involved in the
c            minimum degree nodes elimination.
c            -------------------------------------------
             if  ( num .gt. neqns )  go to 1000
             call  sparse_reorder6(
     &         ehead, neqns, xadj, adjncy, delta, mdeg,
     &         dhead, invp, perm, qsize, llist, marker,
     &         maxint, tag )
             go to 300
c
 1000    continue
         call sparse_reorder9( neqns, perm, invp, qsize )
         return
c
      end
c ****************************************************************
c     sparse_reorder6
c
c     purpose - this routine updates the degrees of nodes
c        after a multiple elimination step.
c
c     input parameters -
c        ehead  - the beginning of the list of eliminated
c                 nodes (i.e., newly formed elements).
c        neqns  - number of equations.
c        (xadj,adjncy) - adjacency structure.
c        delta  - tolerance value for multiple elimination.
c        maxint - maximum machine representable (short)
c                 integer.
c
c     updated parameters -
c        mdeg   - new minimum degree after degree update.
c        (dhead,dforw,dbakw) - degree doubly linked structure.
c        qsize  - size of supernode.
c        llist  - working linked list.
c        marker - marker vector for degree update.
c        tag    - tag value.
c
c***************************************************************
c
c
      subroutine  sparse_reorder6(
     &  ehead, neqns, xadj, adjncy, delta, mdeg, dhead, dforw,
     &  dbakw, qsize, llist, marker, maxint, tag )
c
c***************************************************************
c
         integer    adjncy(1), dbakw(1) , dforw(1) , dhead(1) ,
     1              llist(1) , marker(1), qsize(1)
         integer    xadj(1)
         integer    deg   , deg0  , delta , ehead , elmnt ,
     1              enode , fnode , i     , iq2   , istop ,
     1              istrt , j     , jstop , jstrt , link  ,
     1              maxint, mdeg  , mdeg0 , mtag  , nabor ,
     1              neqns , node  , q2head, qxhead, tag
c
c***************************************************************
c
         mdeg0 = mdeg + delta
         elmnt = ehead
  100    continue
c            -------------------------------------------------------
c            for each of the newly formed element, do the following.
c            (reset tag value if necessary.)
c            -------------------------------------------------------
             if  ( elmnt .le. 0 )  return
             mtag = tag + mdeg0
             if  ( mtag .lt. maxint )  go to 300
                 tag = 1
                 do  200  i = 1, neqns
                     if  ( marker(i) .lt. maxint )  marker(i) = 0
  200            continue
                 mtag = tag + mdeg0
  300        continue
c            ---------------------------------------------
c            create two linked lists from nodes associated
c            with elmnt: one with two nabors (q2head) in
c            adjacency structure, and the other with more
c            than two nabors (qxhead).  also compute deg0,
c            number of nodes in this element.
c            ---------------------------------------------
             q2head = 0
             qxhead = 0
             deg0 = 0
             link = elmnt
  400        continue
                 istrt = xadj(link)
                 istop = xadj(link+1) - 1
                 do  700  i = istrt, istop
                     enode = adjncy(i)
                     link = - enode
                     if  ( enode )  400, 800, 500
c
  500                continue
                     if  ( qsize(enode) .eq. 0 )  go to 700
                         deg0 = deg0 + qsize(enode)
                         marker(enode) = mtag
c                        ----------------------------------
c                        if enode requires a degree update,
c                        then do the following.
c                        ----------------------------------
                         if  ( dbakw(enode) .ne. 0 )  go to 700
c                            ---------------------------------------
c                            place either in qxhead or q2head lists.
c                            ---------------------------------------
                             if  ( dforw(enode) .eq. 2 )  go to 600
                                 llist(enode) = qxhead
                                 qxhead = enode
                                 go to 700
  600                        continue
                             llist(enode) = q2head
                             q2head = enode
  700            continue
  800        continue
c            --------------------------------------------
c            for each enode in q2 list, do the following.
c            --------------------------------------------
             enode = q2head
             iq2 = 1
  900        continue
                 if  ( enode .le. 0 )  go to 1500
                 if  ( dbakw(enode) .ne. 0 )  go to 2200
                     tag = tag + 1
                     deg = deg0
c                    ------------------------------------------
c                    identify the other adjacent element nabor.
c                    ------------------------------------------
                     istrt = xadj(enode)
                     nabor = adjncy(istrt)
                     if  ( nabor .eq. elmnt )  nabor = adjncy(istrt+1)
c                    ------------------------------------------------
c                    if nabor is uneliminated, increase degree count.
c                    ------------------------------------------------
                     link = nabor
                     if  ( dforw(nabor) .lt. 0 )  go to 1000
                         deg = deg + qsize(nabor)
                         go to 2100
 1000                continue
c                        --------------------------------------------
c                        otherwise, for each node in the 2nd element,
c                        do the following.
c                        --------------------------------------------
                         istrt = xadj(link)
                         istop = xadj(link+1) - 1
                         do  1400  i = istrt, istop
                             node = adjncy(i)
                             link = - node
                             if  ( node .eq. enode )  go to 1400
                             if  ( node )  1000, 2100, 1100
c
 1100                        continue
                             if  ( qsize(node) .eq. 0 )  go to 1400
                             if  ( marker(node) .ge. tag )  go to 1200
c                                -------------------------------------
c                                case when node is not yet considered.
c                                -------------------------------------
                                 marker(node) = tag
                                 deg = deg + qsize(node)
                                 go to 1400
 1200                        continue
c                            ----------------------------------------
c                            case when node is indistinguishable from
c                            enode.  merge them into a new supernode.
c                            ----------------------------------------
                             if  ( dbakw(node) .ne. 0 )  go to 1400
                             if  ( dforw(node) .ne. 2 )  go to 1300
                                 qsize(enode) = qsize(enode) +
     1                                          qsize(node)
                                 qsize(node) = 0
                                 marker(node) = maxint
                                 dforw(node) = - enode
                                 dbakw(node) = - maxint
                                 go to 1400
 1300                        continue
c                            --------------------------------------
c                            case when node is outmatched by enode.
c                            --------------------------------------
                             if  ( dbakw(node) .eq.0 )
     1                             dbakw(node) = - maxint
 1400                    continue
                         go to 2100
 1500            continue
c                ------------------------------------------------
c                for each enode in the qx list, do the following.
c                ------------------------------------------------
                 enode = qxhead
                 iq2 = 0
 1600            continue
                     if  ( enode .le. 0 )  go to 2300
                     if  ( dbakw(enode) .ne. 0 )  go to 2200
                         tag = tag + 1
                         deg = deg0
c                        ---------------------------------
c                        for each unmarked nabor of enode,
c                        do the following.
c                        ---------------------------------
                         istrt = xadj(enode)
                         istop = xadj(enode+1) - 1
                         do  2000  i = istrt, istop
                             nabor = adjncy(i)
                             if  ( nabor .eq. 0 )  go to 2100
                             if  ( marker(nabor) .ge. tag )  go to 2000
                                 marker(nabor) = tag
                                 link = nabor
c                                ------------------------------
c                                if uneliminated, include it in
c                                deg count.
c                                ------------------------------
                                 if  ( dforw(nabor) .lt. 0 )  go to 1700
                                     deg = deg + qsize(nabor)
                                     go to 2000
 1700                            continue
c                                    -------------------------------
c                                    if eliminated, include unmarked
c                                    nodes in this element into the
c                                    degree count.
c                                    -------------------------------
                                     jstrt = xadj(link)
                                     jstop = xadj(link+1) - 1
                                     do  1900  j = jstrt, jstop
                                         node = adjncy(j)
                                         link = - node
                                         if  ( node )  1700, 2000, 1800
c
 1800                                    continue
                                         if  ( marker(node) .ge. tag )
     1                                         go to 1900
                                             marker(node) = tag
                                             deg = deg + qsize(node)
 1900                                continue
 2000                    continue
 2100                continue
c                    -------------------------------------------
c                    update external degree of enode in degree
c                    structure, and mdeg (min deg) if necessary.
c                    -------------------------------------------
                     deg = deg - qsize(enode) + 1
                     fnode = dhead(deg)
                     dforw(enode) = fnode
                     dbakw(enode) = - deg
                     if  ( fnode .gt. 0 )  dbakw(fnode) = enode
                     dhead(deg) = enode
                     if  ( deg .lt. mdeg )  mdeg = deg
 2200                continue
c                    ----------------------------------
c                    get next enode in current element.
c                    ----------------------------------
                     enode = llist(enode)
                     if  ( iq2 .eq. 1 )  go to 900
                         go to 1600
 2300        continue
c            -----------------------------
c            get next element in the list.
c            -----------------------------
             tag = mtag
             elmnt = llist(elmnt)
             go to 100
c
      end
c  ***************************************************
c   sparse_reorder7   
c
c     purpose - this routine eliminates the node mdnode of
c        minimum degree from the adjacency structure, which
c        is stored in the quotient graph format.  it also
c        transforms the quotient graph representation of the
c        elimination graph.
c
c     input parameters -
c        mdnode - node of minimum degree.
c        maxint - estimate of maximum representable (short)
c                 integer.
c        tag    - tag value.
c
c     updated parameters -
c        (xadj,adjncy) - updated adjacency structure.
c        (dhead,dforw,dbakw) - degree doubly linked structure.
c        qsize  - size of supernode.
c        marker - marker vector.
c        llist  - temporary linked list of eliminated nabors.
c
c***************************************************************
c
c
      subroutine sparse_reorder7(
     &  mdnode, xadj, adjncy, dhead, dforw, dbakw, qsize, llist,
     &  marker, maxint, tag )
c
c***************************************************************
c
         integer    adjncy(1), dbakw(1) , dforw(1) , dhead(1) ,
     1              llist(1) , marker(1), qsize(1)
         integer    xadj(1)
         integer    elmnt , i     , istop , istrt , j     ,
     1              jstop , jstrt , link  , maxint, mdnode,
     1              nabor , node  , npv   , nqnbrs, nxnode,
     1              pvnode, rlmt  , rloc  , rnode , tag   ,
     1              xqnbr
c
c***************************************************************
c
c        -----------------------------------------------
c        find reachable set and place in data structure.
c        -----------------------------------------------
         marker(mdnode) = tag
         istrt = xadj(mdnode)
         istop = xadj(mdnode+1) - 1
c        -------------------------------------------------------
c        elmnt points to the beginning of the list of eliminated
c        nabors of mdnode, and rloc gives the storage location
c        for the next reachable node.
c        -------------------------------------------------------
         elmnt = 0
         rloc = istrt
         rlmt = istop
         do  200  i = istrt, istop
             nabor = adjncy(i)
             if  ( nabor .eq. 0 )  go to 300
                 if  ( marker(nabor) .ge. tag )  go to 200
                     marker(nabor) = tag
                     if  ( dforw(nabor) .lt. 0 )  go to 100
                         adjncy(rloc) = nabor
                         rloc = rloc + 1
                         go to 200
  100                continue
                     llist(nabor) = elmnt
                     elmnt = nabor
  200    continue
  300    continue
c            -----------------------------------------------------
c            merge with reachable nodes from generalized elements.
c            -----------------------------------------------------
             if  ( elmnt .le. 0 )  go to 1000
                 adjncy(rlmt) = - elmnt
                 link = elmnt
  400            continue
                     jstrt = xadj(link)
                     jstop = xadj(link+1) - 1
                     do  800  j = jstrt, jstop
                         node = adjncy(j)
                         link = - node
                         if  ( node )  400, 900, 500
  500                    continue
                         if  ( marker(node) .ge. tag  .or.
     1                         dforw(node) .lt. 0 )  go to 800
                             marker(node) = tag
c                            ---------------------------------
c                            use storage from eliminated nodes
c                            if necessary.
c                            ---------------------------------
  600                        continue
                                 if  ( rloc .lt. rlmt )  go to 700
                                     link = - adjncy(rlmt)
                                     rloc = xadj(link)
                                     rlmt = xadj(link+1) - 1
                                     go to 600
  700                        continue
                             adjncy(rloc) = node
                             rloc = rloc + 1
  800                continue
  900            continue
                 elmnt = llist(elmnt)
                 go to 300
 1000    continue
         if  ( rloc .le. rlmt )  adjncy(rloc) = 0
c        --------------------------------------------------------
c        for each node in the reachable set, do the following ...
c        --------------------------------------------------------
         link = mdnode
 1100    continue
             istrt = xadj(link)
             istop = xadj(link+1) - 1
             do  1700  i = istrt, istop
                 rnode = adjncy(i)
                 link = - rnode
                 if  ( rnode )  1100, 1800, 1200
 1200            continue
c                --------------------------------------------
c                if rnode is in the degree list structure ...
c                --------------------------------------------
                 pvnode = dbakw(rnode)
                 if  ( pvnode .eq. 0  .or.
     1                 pvnode .eq. (-maxint) )  go to 1300
c                    -------------------------------------
c                    then remove rnode from the structure.
c                    -------------------------------------
                     nxnode = dforw(rnode)
                     if  ( nxnode .gt. 0 )  dbakw(nxnode) = pvnode
                     if  ( pvnode .gt. 0 )  dforw(pvnode) = nxnode
                     npv = - pvnode
                     if  ( pvnode .lt. 0 )  dhead(npv) = nxnode
 1300            continue
c                ----------------------------------------
c                purge inactive quotient nabors of rnode.
c                ----------------------------------------
                 jstrt = xadj(rnode)
                 jstop = xadj(rnode+1) - 1
                 xqnbr = jstrt
                 do  1400  j = jstrt, jstop
                     nabor = adjncy(j)
                     if  ( nabor .eq. 0 )  go to 1500
                         if  ( marker(nabor) .ge. tag )  go to 1400
                             adjncy(xqnbr) = nabor
                             xqnbr = xqnbr + 1
 1400            continue
 1500            continue
c                ----------------------------------------
c                if no active nabor after the purging ...
c                ----------------------------------------
                 nqnbrs = xqnbr - jstrt
                 if  ( nqnbrs .gt. 0 )  go to 1600
c                    -----------------------------
c                    then merge rnode with mdnode.
c                    -----------------------------
                     qsize(mdnode) = qsize(mdnode) + qsize(rnode)
                     qsize(rnode) = 0
                     marker(rnode) = maxint
                     dforw(rnode) = - mdnode
                     dbakw(rnode) = - maxint
                     go to 1700
 1600            continue
c                --------------------------------------
c                else flag rnode for degree update, and
c                add mdnode as a nabor of rnode.
c                --------------------------------------
                 dforw(rnode) = nqnbrs + 1
                 dbakw(rnode) = 0
                 adjncy(xqnbr) = mdnode
                 xqnbr = xqnbr + 1
                 if  ( xqnbr .le. jstop )  adjncy(xqnbr) = 0
c
 1700        continue
 1800    continue
         return
c
      end
c  **********************************************************
c  sparse_reorder8
c
c     purpose - this routine performs initialization for the
c        multiple elimination version of the minimum degree
c        algorithm.
c
c     input parameters -
c        neqns  - number of equations.
c        (xadj,adjncy) - adjacency structure.
c
c     output parameters -
c        (dhead,dforw,dbakw) - degree doubly linked structure.
c        qsize  - size of supernode (initialized to one).
c        llist  - linked list.
c        marker - marker vector.
c
c***************************************************************
c
c
      subroutine sparse_reorder8(
     &  neqns, xadj, dhead, dforw, dbakw, qsize, llist,
     &  marker )
c
c***************************************************************
c
         integer    dbakw(1) , dforw(1) , dhead(1) ,
     1              llist(1) , marker(1), qsize(1)
         integer    xadj(1)
         integer    fnode , ndeg  , neqns , node
c
c***************************************************************
c
         do  node = 1, neqns
             dhead(node) = 0
             qsize(node) = 1
             marker(node) = 0
             llist(node) = 0
         end do
c        ------------------------------------------
c        initialize the degree doubly linked lists.
c        ------------------------------------------
         do node = 1, neqns
             ndeg = xadj(node+1) - xadj(node) + 1
             fnode = dhead(ndeg)
             dforw(node) = fnode
             dhead(ndeg) = node
             if  ( fnode .gt. 0 )  dbakw(fnode) = node
             dbakw(node) = - ndeg
         end do
         return
c
      end
c **********************************************************
c sparse_reorder9
c
c     purpose - this routine performs the final step in
c        producing the permutation and inverse permutation
c        vectors in the multiple elimination version of the
c        minimum degree ordering algorithm.
c
c     input parameters -
c        neqns  - number of equations.
c        qsize  - size of supernodes at elimination.
c
c     updated parameters -
c        invp   - inverse permutation vector.  on input,
c                 if qsize(node)=0, then node has been merged
c                 into the node -invp(node); otherwise,
c                 -invp(node) is its inverse labelling.
c
c     output parameters -
c        perm   - the permutation vector.
c
c***************************************************************
c
c
      subroutine sparse_reorder9(
     &  neqns, perm, invp, qsize )
c
c***************************************************************
c
         integer    invp(1)  , perm(1)  , qsize(1)
         integer    father, neqns , nextf , node  , nqsize,
     1              num   , root
c
c***************************************************************
c
         do node = 1, neqns
             nqsize = qsize(node)
             if  ( nqsize .le. 0 )  perm(node) = invp(node)
             if  ( nqsize .gt. 0 )  perm(node) = - invp(node)
         end do
c        ------------------------------------------------------
c        for each node which has been merged, do the following.
c        ------------------------------------------------------
         do  500  node = 1, neqns
             if  ( perm(node) .gt. 0 )  go to 500
c                -----------------------------------------
c                trace the merged tree until one which has
c                not been merged, call it root.
c                -----------------------------------------
                 father = node
  200            continue
                     if  ( perm(father) .gt. 0 )  go to 300
                         father = - perm(father)
                         go to 200
  300            continue
c                -----------------------
c                number node after root.
c                -----------------------
                 root = father
                 num = perm(root) + 1
                 invp(node) = - num
                 perm(root) = num
c                ------------------------
c                shorten the merged tree.
c                ------------------------
                 father = node
  400            continue
                     nextf = - perm(father)
                     if  ( nextf .le. 0 )  go to 500
                         perm(father) = - root
                         father = nextf
                         go to 400
  500    continue
c        ----------------------
c        ready to compute perm.
c        ----------------------
         do  node = 1, neqns
             num = - invp(node)
             invp(node) = num
             perm(num) = node
         end do
         return
c
      end
