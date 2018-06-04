c ***************************************************************
c *                                                             *
c * domain expand 4 - expand type 4 automatic domain            *
c *                                                             *
c ***************************************************************
c
c
      subroutine diexp4( ring, node_map, nonode, noelem, incmap,
     &                    eprops, incid, out, q_map_len, bits )
c
      use j_data, only : debug_driver, num_front_nodes, q_values,
     &                   expanded_front_nodes, q_element_maps
      use main_data, only : inverse_incidences
      implicit none
      include 'param_def'
c
c                     parameter declarations
c
      integer :: ring, nonode, noelem, out, q_map_len
      integer :: node_map(*), incmap(*), eprops(mxelpr,*), incid(*),
     &           bits(*)
c
c                     local declarations
c
      integer :: front_pos, i, ii, snode, elem, incptr, nnode, elnd,
     &           count, index
      integer :: local_list(10)
      real, parameter :: rzero=0.0, rone=1.0
      logical, external :: dibmck
c
c             build list of all nodes for elements appearing in
c             the previous domain. if ring = 1 we use the list of
c             crack front nodes to initialize the process.
c
      if( debug_driver ) write(out,9000) ring
      node_map(1:q_map_len) = 0
c
      if( ring .eq. 1 ) then
         do front_pos = 1, num_front_nodes
           do ii = 1, expanded_front_nodes(0,front_pos)
             snode = expanded_front_nodes(ii,front_pos)
             call diadit( node_map, snode, bits )
           end do
         end do
      else
         do elem = 1, noelem
          if( dibmck( elem, q_element_maps, bits ) ) then
             incptr = incmap(elem)-1
             nnode  = eprops(2,elem)
             do elnd = 1, nnode
               snode = incid(incptr+elnd)
               call diadit( node_map, snode, bits )
             end do
           end if
         end do
       end if
c
c              debug for list of nodes just created
c
       if( debug_driver ) then
           write(out,9010)
           count = 0
           index = 0
           do snode = 1, nonode
             if( dibmck( snode, node_map, bits ) ) then
               index = index + 1
               local_list(index) = snode
               count = count + 1
               if( mod(index,10) .eq. 0 ) then
                   write(out,9040) local_list
                   index = 0
               end if
             end if
           end do
           if( index .gt. 0 ) write(out,9040) (local_list(ii),
     &                                          ii=1,index)
           write(out,9015) count
       end if
c
c             build list of all elements incident on nodes in
c             current list. use nodal incidences to greatly
c             speed up this process.
c
      do snode = 1, nonode
         if( dibmck( snode, node_map, bits ) ) then
           do i = 1, inverse_incidences(snode)%element_count
             elem = inverse_incidences(snode)%element_list(i)
             call diadit( q_element_maps, elem, bits )
           end do
         end if
      end do
c
c              debug for list of elements to be in domain.
c
       if( debug_driver ) then
           write(out,9030)
           index = 0
           count = 0
           do elem = 1, noelem
             if( dibmck( elem, q_element_maps, bits ) ) then
               index = index + 1
               local_list(index) = elem
               count = count + 1
               if( mod(index,10) .eq. 0 ) then
                  write(out,9040) local_list
                  index = 0
               end if
             end if
           end do
           if( index .gt. 0 ) write(out,9040) (local_list(ii),
     &                                          ii=1,index)
           write(out,9035) count
       end if
c
c             for all nodes in the node map, set the q-value = 1
c
      do snode = 1, nonode
         if( dibmck( snode, node_map, bits ) ) then
           q_values(snode) = rone
         else
           q_values(snode) = rzero
         end if
      end do
c
      return
c
 9000 format(//,1x,'>>> entered expand_4, for ring: ',i3,' ...')
 9010 format(2x,'>>> nodes in current map:')
 9015 format(2x,'>>> number of nodes in current map:',i5)
 9020 format(1x,i10)
 9030 format(2x,'>>> elements to be in domain:')
 9035 format(2x,'>>> number of elements to be in domain:',i5)
 9040 format(10x,10i7)
c
      end
c ***************************************************************
c *                                                             *
c * domain expand 13 - expand type 1-3 automatic domain         *
c *                                                             *
c ***************************************************************
c
      subroutine diexp13( ring, nonode, noelem, incmap, eprops, incid,
     &                    out, bits, newmap, nodmap, q_map_len )
      use j_data, only : debug_driver, q_values, domain_type,
     a  q_element_maps, num_front_nodes, front_order,
     b  expanded_front_nodes
      use main_data, only : inverse_incidences
      implicit none
      include 'param_def'
c
c                     parameter declarations
c
      integer :: ring, nonode, noelem, out, q_map_len
      integer :: incmap(*), eprops(mxelpr,*), incid(*), bits(*),
     &           nodmap(q_map_len,*), newmap(q_map_len,*)
c
c                     local declarations
c
      integer :: ii, elem, cornno, frnpos, count, fnode, snode, elenum,
     &           incptr, jj, edge, nenode, ecorn1, ecorn2, scorn1,
     &           scorn2, i, index, nnode, numedge, mapno, nownod
      integer :: local_list(10), corner(3), enodes(30),
     &           edge_table(0:3,20)
      logical, external :: dibmck
      logical :: flags1(3), flags2(3), fl1, fl2
      real, parameter :: rzero=0.0, rone=1.0
c
      if( debug_driver ) write(out,9000) ring
      if( ring == 1 ) call diexp13_ring_1
      if( ring > 1 )  call diexp13_not_ring_1
c
c             create definition for next ring. for ring 1, the
c             front node map just created starts the process.
c
c             set q-values and element list for the domain.
c
      do snode = 1, nonode
       q_values(snode) = rzero
       if( .not. dibmck( snode, newmap(1,domain_type), bits ) ) cycle
       q_values(snode) = rone
       do i = 1, inverse_incidences(snode)%element_count
         elem = inverse_incidences(snode)%element_list(i)
         call diadit( q_element_maps, elem, bits )
       end do
      end do
c
      if( debug_driver ) call diexp13_debug
c
c             copy the 3 node maps for current ring into old maps
c             to set up generation of next ring.
c
      do ii = 1, q_map_len
        nodmap(ii,1) = newmap(ii,1)
        nodmap(ii,2) = newmap(ii,2)
        nodmap(ii,3) = newmap(ii,3)
      end do
c
c             all done building domains.
c
      return
c
 9000 format(//,1x,'>>> entered expand_13, for ring: ',i3,' ...')
c
      contains
c     ========
c
      subroutine diexp13_ring_1
      implicit none
c
c             convert list of coincident nodes at each corner
c             node on the front to a bit map. there are 2 or three
c             corner nodes. for 2 corner nodes we just leave
c             one of the maps all zero which makes code logic
c             simpler in later steps. set which positions in the list
c             of front nodes are the corner nodes.
c
      if( domain_type .eq. 1 ) then
         corner(1) = 1
         corner(2) = num_front_nodes
         corner(3) = 0
        else if ( domain_type .eq. 2 ) then
         corner(1) = 1
         corner(2) = front_order + 1
         corner(3) = num_front_nodes
        else
         corner(1) = 0
         corner(2) = 1
         corner(3) = num_front_nodes
      end if
c
      do ii = 1, q_map_len
        newmap(ii,1) = 0
        newmap(ii,2) = 0
        newmap(ii,3) = 0
      end do
      if( debug_driver )  write(out,9400) corner
c
      do cornno = 1, 3
        frnpos = corner(cornno)
        if( frnpos .eq. 0 ) cycle
        count = expanded_front_nodes(0,frnpos)
        do fnode = 1, count
          snode = expanded_front_nodes(fnode,frnpos)
          call diadit( newmap(1,cornno), snode, bits )
        end do
      end do
      if( debug_driver )  write(out,9500)
      return
c
 9400 format(/,'       corner(1-3): ',4i3 )
 9500 format(/,'       initial 3 node maps created')
c
      end subroutine diexp13_ring_1
c
      subroutine diexp13_not_ring_1
      implicit none
c
      do elenum = 1, noelem
        if( .not. dibmck( elenum, q_element_maps, bits ) ) cycle
c
c             get element edge topology data and incidences.
c
        call digete( eprops(1,elenum), nnode, numedge, edge_table )
        incptr = incmap(elenum)
        jj = incptr - 1
        do ii = 1, nnode
         enodes(ii) = incid(jj+ii)
        end do
        if( debug_driver ) write(out,9700) elenum, numedge
c
c             loop over all edges of this element. get the two
c             structure nodes corresponding to the two corner nodes
c             for the edge. classify the edge based on it's
c             connectivity to nodes listed in the three node maps.
c             build 3 logical flags for each corner node, marking
c             it's appearance in the 3 node maps. based on connectivity
c             of edge nodes on the 3 nodes lists, add a corner
c             node to the lists.  logic tree:
c
c               0 = both corner nodes appear in the 3 node lists
c               0 = neither corner node appears in the 3 node lists
c               if corner node 1 only appears in the 3 lists,
c                  add corner node 2 to the node list in which corner
c                  node 1 appears.
c               if corner node 2 only appears in the 3 lists,
c                  add corner node 1 to the node list in which corner
c                  node 2 appears.
c
        if( debug_driver ) write(out,9800)
        do edge = 1, numedge
           nenode = edge_table(0,edge)
           ecorn1 = edge_table(1,edge)
           ecorn2 = edge_table(1+nenode-1,edge)
           scorn1 = enodes(ecorn1)
           scorn2 = enodes(ecorn2)
           flags1(1) = dibmck( scorn1, nodmap(1,1), bits )
           flags1(2) = dibmck( scorn1, nodmap(1,2), bits )
           flags1(3) = dibmck( scorn1, nodmap(1,3), bits )
           flags2(1) = dibmck( scorn2, nodmap(1,1), bits )
           flags2(2) = dibmck( scorn2, nodmap(1,2), bits )
           flags2(3) = dibmck( scorn2, nodmap(1,3), bits )
           fl1 = flags1(1) .or. flags1(2) .or. flags1(3)
           fl2 = flags2(1) .or. flags2(2) .or. flags2(3)
           if ( debug_driver ) write(out,9800) edge, nenode, ecorn1,
     &                         ecorn2, scorn1, scorn2, flags1, flags2,
     &                         fl1, fl2
           if( fl1 .and. fl2 ) cycle
           if( (.not. fl1) .and. (.not. fl2) ) cycle
           if( fl1 ) then
               if ( flags1(1) ) mapno = 1
               if ( flags1(2) ) mapno = 2
               if ( flags1(3) ) mapno = 3
               nownod = scorn2
             else
               if ( flags2(1) ) mapno = 1
               if ( flags2(2) ) mapno = 2
               if ( flags2(3) ) mapno = 3
               nownod = scorn1
           end if
           call diadit( newmap(1,mapno), nownod, bits )
           if( debug_driver ) write(out,9900) nownod, mapno
        end do
      end do
c
      return
c
 9700 format(/,'      element set up data: ',
     &  /,     '       ielem, elenum, edaddr, neterm:',4i8,
     &  /,     '       numedg, incvec, incpos, nenode: ',4i8 )
 9800 format(/,'      processing element edge: ',i4,
     &  /,     '       nenode,ecorn1,ecorn2,scorn1,scorn2: ',5i5,
     &  /,     '       flags1, flags2, fl1, fl2: ',8l2 )
 9900 format(/,'       saving node. mapno, node: ',2i4)
      end subroutine diexp13_not_ring_1
c
      subroutine diexp13_debug
      implicit none
c
      write(out,9010)
      count = 0
      index = 0
      do snode = 1, nonode
        if( dibmck( snode, newmap(1,domain_type), bits ) ) then
          index = index + 1
          local_list(index) = snode
          count = count + 1
          if( mod(index,10) .eq. 0 ) then
              write(out,9040) local_list
              index = 0
          end if
        end if
      end do
      if( index .gt. 0 ) write(out,9040) (local_list(ii),
     &                                     ii=1,index)
      write(out,9015) count
c
c         debug for list of elements to be in domain.
c
      write(out,9030)
      index = 0
      count = 0
      do elem = 1, noelem
        if( dibmck( elem, q_element_maps, bits ) ) then
          index = index + 1
          local_list(index) = elem
          count = count + 1
          if( mod(index,10) .eq. 0 ) then
             write(out,9040) local_list
             index = 0
          end if
        end if
      end do
      if( index .gt. 0 ) write(out,9040) (local_list(ii),
     &                                     ii=1,index)
      write(out,9035) count
c
      return
c
 9010 format(2x,'>>> nodes in current map:')
 9015 format(2x,'>>> number of nodes in current map:',i5)
 9020 format(i10)
 9030 format(2x,'>>> elements to be in domain:')
 9035 format(2x,'>>> number of elements to be in domain:',i5)
 9040 format(10x,10i7)
 9500 format(/,'       initial 3 node maps created')
 9600 format(/,' >>> building domain for ring: ',i4)
c
      end subroutine diexp13_debug
      end subroutine diexp13
c *******************************************************************
c *                                                                 *
c *   support action -- add to bit map from specified item          *
c *                                                                 *
c *******************************************************************
c
c
      subroutine diadit( mapvec, item, bits )
      implicit none
      integer :: mapvec(*), bits(*), item
      integer :: word, bit
c
c             mapvec contains a bit map. item is an entity
c             representing a node, element, etc.
c             turn on bit in bit map for integer.
c
      word = ( item - 1 ) / 30  + 1
      bit  = item - ( word-1 ) * 30
      mapvec(word) = ior( mapvec(word),bits(bit) )
c
      return
      end
c *******************************************************************
c *                                                                 *
c *   support action -- test for entry presence in a bit map        *
c *                                                                 *
c *******************************************************************
c
c
      function dibmck( entry, mapvec, bits ) result( ans )
      implicit none
      logical :: ans
      integer, intent(in) :: mapvec(*), bits(*), entry
      integer :: word, bit
c
c             mapvec contains a bit map representing a set of
c             nodes, elements, etc. entry is a node number, element
c             number, etc. determine if the entry is present in the
c             bit map.
c
      word = ( entry-1 ) / 30  + 1
      bit  = entry - ( word-1 ) * 30
      ans  = iand( bits(bit),mapvec(word) ) .ne. 0
c
      return
      end
c *******************************************************************
c *                                                                 *
c *   support action -- get element edge data                       *
c *                                                                 *
c *******************************************************************
c
      subroutine digete( eprops, nnode, numedg, edge_table )
      implicit none
c
c             paramter declarations
c
      integer :: nnode, numedg,  eprops(*), edge_table(0:3,*)
c
c             local declarations
c
      integer :: etype, j
      integer :: edge_8_table(0:2,12),edge_20_table(0:3,12),
     &           edge_12_table(0:3,12),edge_15_table(0:3,12),
     &           edge_9_table(0:3,12)
      data edge_8_table /
     & 2,1,5,  2,5,8,  2,8,4,  2,4,1,  2,2,6,  2,6,7,  2,7,3,
     & 2,3,2,  2,7,8,  2,3,4,  2,6,5,  2,2,1 /
      data edge_20_table /
     & 3,1,17,5,  3,5,16,8,  3,8,20,4,  3,4,12,1,  3,2,18,6,
     & 3,6,14,7,  3,7,19,3,  3,3,10,2,  3,7,15,8,  3,3,11,4,
     & 3,6,13,5,  3,2,9,1 /
      data edge_12_table /
     & 2,1,5,0,   2,5,8,0,   2,8,4,0,   3,4,12,1,  2,2,6,0,
     & 2,6,7,0,   2,7,3,0,   3,3,10,2,  2,7,8,0,   3,3,11,4,
     & 2,6,5,0,   3,2,9,1 /
      data edge_9_table /
     & 2,1,5,0,   2,5,8,0,   2,8,4,0,   2,4,1,0,   2,2,6,0,
     & 2,6,7,0,   2,7,3,0,   2,3,2,0,   2,7,8,0,   2,3,4,0,
     & 2,6,5,0,   3,2,9,1 /
      data edge_15_table /
     & 3,1,14,5,  2,5,8,0,   2,8,4,0,   3,4,12,1,  3,2,15,6,
     & 2,6,7,0,   2,7,3,0,   3,3,10,2,  2,7,8,0,   3,3,11,4,
     & 3,6,13,5,  3,2,9,1 /
c
      etype = eprops(1)
      nnode = eprops(2)
c
      select case( etype )

        case( 1 ) ! 20-node brick
        numedg = 12
        do j = 1, 12
              edge_table(0,j) = edge_20_table(0,j)
              edge_table(1,j) = edge_20_table(1,j)
              edge_table(2,j) = edge_20_table(2,j)
              edge_table(3,j) = edge_20_table(3,j)
        end do

        case( 2)  ! 8-node brick
        numedg = 12
        do j = 1, 12
           edge_table(0,j) = edge_8_table(0,j)
           edge_table(1,j) = edge_8_table(1,j)
           edge_table(2,j) = edge_8_table(2,j)
        end do

        case( 3 ) ! 12-node brick
        numedg = 12
        do j = 1, 12
           edge_table(0,j) = edge_12_table(0,j)
           edge_table(1,j) = edge_12_table(1,j)
           edge_table(2,j) = edge_12_table(2,j)
           edge_table(3,j) = edge_12_table(3,j)
        end do
c

        case( 4 ) ! 15-node brick
        numedg = 12
        do j = 1, 12
           edge_table(0,j) = edge_15_table(0,j)
           edge_table(1,j) = edge_15_table(1,j)
           edge_table(2,j) = edge_15_table(2,j)
           edge_table(3,j) = edge_15_table(3,j)
        end do
        return
c
        case( 5 ) !  9-node brick
        numedg = 12
        do j = 1, 12
           edge_table(0,j) = edge_9_table(0,j)
           edge_table(1,j) = edge_9_table(1,j)
           edge_table(2,j) = edge_9_table(2,j)
           edge_table(3,j) = edge_9_table(3,j)
        end do
c
        case default
          write(*,9000)
          call die_abort
c
      end select
c
      return
 9000 format('** FATAL ERROR: routine digete. bad etype',
     & /,    '                job terminated.' )
      end



