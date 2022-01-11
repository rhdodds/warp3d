c ***************************************************************
c *                                                             *
c * set_up_domain  -- set up the domain for processing          *
c *                   find coincident front nodes and set their *
c *                   q-values                                  *
c *                                                             *
c *            last updated:  1/7/22 rhd                        *
c *                                                             *
c ***************************************************************
c
c
      subroutine distup( scoord, coord_map, num_struct_nodes, debug,
     &                   iout )
c
      use j_data, only : front_nodes, domain_node_sets,
     & expanded_front_nodes, rings_given, domain_type, q_values,
     & verify_front, num_front_nodes, front_order, box_tol_relat,
     & q_vals_linear
c
      implicit none
c
c                     parameter declarations
c
      integer :: coord_map(*), num_struct_nodes, iout
      double precision ::  scoord(*)
      logical :: debug
c
c                     local declarations
c
      integer, allocatable :: coinc_list(:)
      integer :: count,fnode,
     &           front_pos, param, i, ii, pos, snode1, snode2, snode3,
     &           node1, node2, coord_loc, snode, jj
      double precision ::  x, y, z, dist, x1, y1, z1, dx, dy, dz,
     &                    x2, y2, z2, xp, xm, yp, ym, zp, zm, dumd,
     &                    box_tol
      real :: dumr
      real, parameter :: rone=1.0, rhalf=0.5
      logical :: inside_box, user_def_ct
      character(len=1) :: dums
      double precision :: zero=0.0d0, toler=0.001d0
c
      if( debug ) then
         write(iout,*) ' >>> entered set_up_domain'
         write(iout,*) '     > no. front nodes, interpolation: ',
     &                 num_front_nodes, front_order
      end if
c
      allocate( coinc_list(1000) )
c
c             determine if the user has specificed the crack front nodes
c             using node sets. this is done when the model is not
c             collapsed at the crack tip or the user wants to define a 
c             non-std starting domain shape as domain 1. 
c             For this user defined crack tip, load
c             expanded_front_nodes with node_set specified by the
c             user. check to be sure that each front position has the
c             same number of crack tip nodes.
c
      user_def_ct = front_nodes(1) .lt. 0
      if( user_def_ct ) call distup_user_ct
c
c             for automatic domain definitions, set the q-values at
c             user specified front nodes based on the domain type
c             and the user specified front interpolation order.
c
      if( rings_given ) then
        if( domain_type .eq. 4 ) then
          do ii = 1, num_front_nodes
            q_values(front_nodes(ii)) = rone
          end do
        end if
        if( domain_type .eq. 1 ) q_values(front_nodes(1)) = rone
        if( domain_type .eq. 3 )
     &               q_values(front_nodes(num_front_nodes)) = rone
        if( domain_type .eq. 2 ) then
          if( front_order .eq. 1 ) pos = 2
          if( front_order .eq. 2 ) pos = 3
          q_values(front_nodes(pos)) = rone
        end if
        q_vals_linear = .true.
      end if
c
c             if linear interpolation of q-values over the domain
c             is on, we also need to interpolate q-values along
c             front for quadratic elements.
c
      if( q_vals_linear .and. front_order .eq. 2 ) then
        do i = 1, num_front_nodes-2, 2
          snode1 = front_nodes(i)
          snode2 = front_nodes(i+1)
          snode3 = front_nodes(i+2)
          q_values(snode2) = rhalf * ( q_values(snode1) +
     &                                 q_values(snode3) )
       end do
      end if
c
c             for the special user-defined front nodes, we're done.
c             no need to search for coincident nodes at each front
c             location
c
      if( user_def_ct ) then
         deallocate( coinc_list )
         return
      end if
c
c             look at first 1 or 2 nodes specified to be on the
c             crack front. compute a distance to be used for locating
c             all other nodes conincident with each specified
c             crack front node. set the "box" size for proximity
c             checking as a fraction of the distance measure.
c
      node1 = front_nodes(1)
      node2 = 0
      if( num_front_nodes .eq. 1 ) then
         coord_loc = coord_map(node1) - 1
         x         = scoord(coord_loc+1)
         y         = scoord(coord_loc+2)
         z         = scoord(coord_loc+3)
         dist      = sqrt( x*x + y*y + z*z )
         if( dist .lt. toler ) dist = toler
       else
         node1      = front_nodes(1)
         node2      = front_nodes(2)
         coord_loc  = coord_map(node1) - 1
         x1         = scoord(coord_loc+1)
         y1         = scoord(coord_loc+2)
         z1         = scoord(coord_loc+3)
         coord_loc  = coord_map(node2) - 1
         x2         = scoord(coord_loc+1)
         y2         = scoord(coord_loc+2)
         z2         = scoord(coord_loc+3)
         dx         = abs(x1-x2)
         dy         = abs(y1-y2)
         dz         = abs(z1-z2)
         dist       = sqrt( dx*dx + dy*dy + dz*dz )
      end if
      box_tol = box_tol_relat * dist
      if( debug ) then
         write(iout,*) '     > first 2 nodes on front:     ',node1,node2
         write(iout,*) '     > distance between nodes:     ',dist
         write(iout,*) '     > front-node box tolerance:   ',box_tol
      end if
c
c             loop over each specified front node. locate all other
c             nodes in model that are conincident with the front node.
c             (use the box test). set the q-value of each coincident
c             node to be same as user specified value for the specified
c             front node. dump verification list if user requested it.
c
      if( verify_front ) then
         write(iout,8000) box_tol_relat
         write(iout,8500) box_tol
         write(iout,9000)
      end if
c
      do front_pos = 1, num_front_nodes
        fnode = front_nodes(front_pos)
        count = 1
        coinc_list(count) = fnode
        coord_loc  = coord_map(fnode) - 1
        x          = scoord(coord_loc+1)
        y          = scoord(coord_loc+2)
        z          = scoord(coord_loc+3)
        xp         = x + box_tol
        xm         = x - box_tol
        yp         = y + box_tol
        ym         = y - box_tol
        zp         = z + box_tol
        zm         = z - box_tol
        do snode = 1, num_struct_nodes
          if( snode .eq. fnode ) cycle
          coord_loc  = coord_map(snode) - 1
          x          = scoord(coord_loc+1)
          y          = scoord(coord_loc+2)
          z          = scoord(coord_loc+3)
          inside_box = .true.
          if( x .lt. xm .or. x .gt. xp ) inside_box = .false.
          if( y .lt. ym .or. y .gt. yp ) inside_box = .false.
          if( z .lt. zm .or. z .gt. zp ) inside_box = .false.
          if( inside_box ) then
             count = count + 1
             if( count .gt. 1000 ) then
              write(iout,9300) count, 1000
              call die_abort
              stop
             end if
             coinc_list(count) = snode
             q_values(snode)   = q_values(fnode)
          end if
        end do ! snode
        if( verify_front ) then
         if( count .gt. 1 ) write(iout,9100) fnode,
     &                       (coinc_list(ii),ii=2,count)
         if( count .eq. 1 ) write(iout,9200) fnode
        end if
        expanded_front_nodes(front_pos)%node_count = count
        allocate( expanded_front_nodes(front_pos)%node_list(count) )
        expanded_front_nodes(front_pos)%node_list(1:count) = 
     &           coinc_list(1:count)
      end do ! front_pos
c
      deallocate( coinc_list )
      return
c
 8000 format(1x,'box tolerance ratio: ',e11.4)
 8500 format(1x,'box tolerance:       ',e11.4)
 9000 format(/,' coincident nodes on crack front:',
     &  //,'   listed node           /--- coincident nodes ---/',/)
 9100 format(i10,5x,20(7i8,/,15x))
 9200 format(i16,/,14x,'none')
 9300 format('** FATAL ERROR: routine distup. number of coincident',
     & /,    '                front nodes: ',i4,' is > dimensioned',
     & /,    '                size of = ',i4,
     & /,    '                job terminated.' )
c
      contains
c     --------
c
c ***************************************************************
c *                                                             *
c *                  internal routine distup_user_ct            *
c *                                                             *
c *                 last updated:  1/7/22 rhd                   *
c *                                                             *
c ***************************************************************
c
c
      subroutine distup_user_ct
      implicit none
c
      integer :: num_front_sets, node_set_id, nnl, num_node_chk, fs
      logical :: message
c
      message = .true.
c
      num_front_sets = num_front_nodes
c
      do fs = 1, num_front_sets
        node_set_id = -front_nodes(fs)
        associate( lst => domain_node_sets(node_set_id)%node_list ) 
          nnl = domain_node_sets(node_set_id)%node_count
          allocate( expanded_front_nodes(fs)%node_list(nnl) )
          expanded_front_nodes(fs)%node_count = nnl
          expanded_front_nodes(fs)%node_list(1:nnl) = lst(1:nnl)
          front_nodes(fs) = lst(1)
        end associate
      end do
c
c             each front set must have same # nodes
c
      num_node_chk = expanded_front_nodes(1)%node_count
      do fs = 2, num_front_sets
        if( expanded_front_nodes(fs)%node_count .ne. num_node_chk
     &         .and. message ) then
         call errmsg( 264, param, dums, dumr, dumd )
         message = .false.
        end if
      end do
c
      return
c
      end subroutine distup_user_ct
      end subroutine distup
