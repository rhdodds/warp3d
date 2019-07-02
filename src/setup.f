c     ****************************************************************
c     *                                                              *
c     *                      subroutine setup                        *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 7/1/2019 rhd               *
c     *                                                              *
c     *     this subroutine sets necessary mapping vectors and       *
c     *     arrays and various other necessary data in preparation   *
c     *     for the solution of the current structure.               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine setup( fatal )
      use global_data ! old common.main
c
      use main_data,       only : invdst, incid, incmap,
     &                            elems_to_blocks, rload,
     &                            repeat_incid, inverse_incidences,
     &                            inverse_dof_map
c
      implicit none
c
      logical :: fatal
c
      integer :: new_inctop, nnode, old_inc_map, enode, dof, ndof,
     &           stnd, elem, nodf, incptr, elnd, chknod, pos, 
     &           ecount, what_node, ldofn, i, j, type,
     &           nstrou, stnod
      integer, external :: clmloc
      logical :: message_flag
      real :: dumr
      character :: dums
      double precision, parameter :: zero = 0.0d0
      double precision :: dumd
      integer, allocatable, dimension(:) :: temp_invdst, new_incid,
     &                                      new_incmap
c
c              incidences have been checked for holes, duplicates.
c              the user may have given data with elements not in
c              sequence. rebuild the incmap and incid vectors 
c              so that incidences for element n+1 follow immediately
c              those for element n.
c
c              there are places in the code where incidences for
c              a block of elements are *assumed* to be contiguous in
c              incid. 
c
c              there are places like this:
c
c              call xxx( incid(incmap(felem)), num_nodes_elem, ..
c
c              subroutine xxx( belinc, num_enodes, ...
c              integer :: belinc(num_enodes,*) .....)
c
c              where felem is the first element in a block.
c
      allocate( new_incmap(noelem) )
      allocate( new_incid(inctop) )
      new_inctop = 0
c
      do elem = 1, noelem
       nnode = iprops(2,elem)
       old_inc_map = incmap(elem)
       new_incmap(elem) = new_inctop + 1
       do enode = 1, nnode
        new_inctop = new_inctop + 1
        new_incid(new_inctop) = incid(old_inc_map+enode-1)
       end do
      end do
      if( inctop .ne. new_inctop ) then
            write(*,*) '.... bad new_inctop:', new_inctop
            call die_gracefully
      end if
      call move_alloc( new_incmap, incmap )
      call move_alloc( new_incid, incid )
c
c              set the inverse mappings for the structure.
c              run only by manager process.
c
c              allocate data structure for nodal incidences,
c              i.e. the list of elements connected to
c              each model node. zero element counts for nodes
c
      call init_maps( 0, 4 )
      do stnd = 1, nonode
         inverse_incidences(stnd)%element_count = 0
      end do
c
c              allocate the logical array indicating whether
c              given element has repeated nodes in its
c              incidences (collapsed elements).
c
      call init_maps( 0, 3 )
c
c              find the maximum element to node connectivity.
c              the maximum degree of connectivity allowed
c              is fixed by a parameter variable (mxconn) used
c              size quite a few, smaller arrays throughout
c              the code
c
      do elem = 1, noelem
         repeat_incid(elem) = .false.
         nnode  = iprops(2,elem)
         ndof   = iprops(4,elem)
         incptr = incmap(elem)-1
         do elnd = 1, nnode
            stnd = incid(incptr+elnd)
c
c              check to see if the current node occurred earlier 
c              in the incidence list-- if so, skip this node
c
            do chknod = 1, elnd - 1
               if( incid(incptr+chknod) .eq. stnd ) then
                  repeat_incid(elem) = .true.
                  go to 10
               end if
            end do
            ecount = inverse_incidences(stnd)%element_count
            ecount = ecount + 1
            inverse_incidences(stnd)%element_count = ecount
            if( ecount .gt. mxconn ) then
               write(out,9100) stnd, mxconn
               call die_gracefully
               stop
            end if
 10         continue
         end do
      end do
c
c              find node with largest number of connected
c              elements. issue information message to user
c
      lgnmcn = 0
      do stnd = 1, nonode
         ecount = inverse_incidences(stnd)%element_count
         if( ecount .gt. lgnmcn ) then
            lgnmcn = ecount
            what_node = stnd
         end if
      end do
      write(out,9000) lgnmcn, what_node
c
c              allocate the inverse dof destination mapping, allocate
c              the list of elements connected to each model node.
c              zero list counts again. they are filled to final 
c              values again as we fill lists of elements on the node
c
      call init_maps( 0, 1 )
      call init_maps( 0, 5 )
      call init_maps( 0, 6 )
      call init_maps( 0, 7 )
      do stnd = 1, nonode
         inverse_incidences(stnd)%element_count = 0
      end do
c
c              build the lists of elements connected to each node.
c              build the inverse dof maps. a set of maps is built for
c              each node. for a node, the ecount x 3 array is built.
c              ecount = number of elements connected to the model 
c              node. let i be the ith element connected to the node.
c              then the inverse dof map (i,j) sets the element
c              local dof number for that element for dof j (j=1,2,3)
c
      do elem = 1,noelem
         nnode  = iprops(2,elem)
         ndof   = iprops(4,elem)
         incptr = incmap(elem)-1
         do elnd = 1, nnode
            stnd = incid(incptr+elnd)
c
c
c              check to see if the current node occurred earlier in the
c              incidence list-- if so, skip this node
c
            do chknod = 1, elnd - 1
               if( incid(incptr+chknod) .eq. stnd ) go to 20
            end do
            ecount = inverse_incidences(stnd)%element_count
            ecount = ecount + 1
            inverse_incidences(stnd)%element_count = ecount
            inverse_incidences(stnd)%element_list(ecount) = elem
            do dof = 1, 3
               inverse_dof_map(stnd)%edof_table(ecount,dof) =
     &                     nnode*(dof-1)+elnd
            end do
 20         continue
         end do
      end do
c
c
c              make sure all nodes have at least one element
c              attached. if not set fatal error flag to prevent 
c              solution.  write list of nodes
c
      message_flag = .true.
      do stnd = 1, nonode
       if( inverse_incidences(stnd)%element_count .eq. 0 ) then
         if( message_flag ) write(out,9200)
         message_flag = .false.
         write(out,*) '   ',stnd
       end if
      end do
c
      if( .not. message_flag ) then
         num_fatal = num_fatal + 1
         write(out,*) ' '
      end if
c
c
c              set the destination arrays and the number of total 
c              degrees of freedom. build list in temporary vector of
c              maximum possible length. allocate permanant vector 
c              of correct length and copy.
c
      allocate( temp_invdst(nonode*mxndof) )
      ldofn = 0
      do stnod = 1, nonode
         elem = inverse_incidences(stnod)%element_list(1)
         ndof         = iprops(4,elem)
         dstmap(stnod) = ldofn + 1
         do i = 1, ndof
            temp_invdst(ldofn+i) = stnod
         end do
         ldofn = ldofn + ndof
      end do
      nodof = ldofn
c
      call init_maps( 0, 2 )
      do i = 1, nodof
        invdst(i) = temp_invdst(i)
      end do
      deallocate( temp_invdst )
c
c
c              building of cdest, edest now done by incomp.cinitialize
c              all necessary vectors which depend on nodof
c
      do i = 1, nodof
         u(i)     = zero
         du(i)    = zero
         v(i)     = zero
         a(i)     = zero
         ifv(i)   = zero
         load(i)  = zero
         rload(i) = zero
      end do
c
c
c              set the maximum output map entry to be accessed.
c
c
      lgoump = 0
      do j = 1, noelem
         type   = iprops(1,j)
         nstrou = iprops(3,j)/two16
         do i = 1, nstrou
            lgoump = max(outmap(type,i),lgoump)
         end do
      end do
c
c
c              create the column position data structures for use 
c              with any upper triangular operations.
c
      do j = 1, mxedof
         cp(j)  = clmloc(j)
         dcp(j) = cp(j)+j
         do i = 1, j
            pos        = cp(j)+i
            icp(pos,1) = i
            icp(pos,2) = j
         end do
      end do
c
 9999 return
 9000 format(/,'>> Maximum element-to-node connectivity: ',i3,
     &       ' @ node: ',i7,//)
 9100 format(/,'>> Fatal Error. Internal table overflow.',
     &       /,'   routine setup, element-to-node connectivity',
     &       /,'   at structure node: ',i7,' exceeds fixed',
     &       /,'   limit (mxconn) of: ',i5,
     &       /,'   Job terminated....',//)
 9200 format(/,'>> Fatal Error. the following nodes have no elements',
     &             ' attached: ' )
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine setup_slave                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/15/01 rhd               *
c     *                                                              *
c     *     this subroutine sets up the inverse mapping and          *
c     *     inverse dof mapping data structures. for a slave         *
c     *     processor                                                *
c     *                                                              *
c     ****************************************************************
c
      subroutine setup_slave( element_node_counts )
      use global_data ! old common.main
c
      use main_data,       only : incid, incmap, inverse_incidences,
     &                            inverse_dof_map
c
      implicit integer (a-z)
      dimension element_node_counts(*)
c
c                       set the inverse mappings for the structure.
c                       allocate data structure for nodal incidences,
c                       i.e. the list of elements connected to
c                       each model node. zero element counts for nodes
c
      call init_maps( 0, 4 )
      do stnd = 1, nonode
         inverse_incidences(stnd)%element_count = 0
      end do
c
c                       build count of number of elements connected to
c                       each model node. if element has repeated node
c                       numbers the element counts just once.
c
      do elem = 1, noelem
         nnode  = element_node_counts(elem)
         ndof   = 3
         incptr = incmap(elem)-1
         do elnd = 1, nnode
            stnd = incid(incptr+elnd)
c
c                           check to see if the current node
c                           occurred earlier in the incidence
c                           list-- if so, skip this node
c
            do chknod = 1, elnd - 1
               if ( incid(incptr+chknod) .eq. stnd ) go to 10
            end do
            ecount = inverse_incidences(stnd)%element_count
            ecount = ecount + 1
            inverse_incidences(stnd)%element_count = ecount
 10         continue
         end do
      end do
c
c                         allocate the inverse dof destination
c                         mapping, allocate the list of elements
c                         connected to each model node.
c                         zero list counts again. they are filled
c                         to final values again as we fill lists of
c                         elements on the node
c
      call init_maps( 0, 5 )
      call init_maps( 0, 6 )
      call init_maps( 0, 7 )
      do stnd = 1, nonode
         inverse_incidences(stnd)%element_count = 0
      end do
c
c                         build the lists of elements connected to
c                         each node. build the inverse dof maps.
c                         a set of maps is built for each node.
c                         for a node, the ecount x 3 array is built.
c                         ecount = number of elements connected to the
c                         model node. let i be the ith element
c                         connected to the node. then the inverse
c                         dof map (i,j) sets the element local dof number
c                         for that element for dof j (j=1,2,3)
c
      do elem = 1, noelem
         nnode  = element_node_counts(elem)
         ndof   = 3
         incptr = incmap(elem)-1
         do elnd = 1, nnode
            stnd = incid(incptr+elnd)
c
c                           check to see if the current node
c                           occurred earlier in the incidence
c                           list-- if so, skip this node
c
            do chknod = 1, elnd - 1
               if ( incid(incptr+chknod) .eq. stnd) go to 20
            end do
            ecount = inverse_incidences(stnd)%element_count
            ecount = ecount + 1
            inverse_incidences(stnd)%element_count = ecount
            inverse_incidences(stnd)%element_list(ecount) = elem
            do dof = 1, 3
               inverse_dof_map(stnd)%edof_table(ecount,dof) =
     &                     nnode*(dof-1)+elnd
            end do
 20         continue
         end do
      end do
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                 function clmloc                              *
c     *                                                              *
c     ****************************************************************
c

      integer function clmloc( col )
      implicit none
c
      integer :: col, stpos, i
      stpos = 0
      do i = 1, col-1
         stpos = stpos + col - i
      end do
c
      clmloc= stpos
c
      return
      end

