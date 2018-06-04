c     ****************************************************************
c     *                                                              *
c     *                  subroutine dimrot                           *
c     *                                                              *
c     *                  written by : rhd                            *
c     *                                                              *
c     *             last modified : 2/17/2018 rhd                    *
c     *                                                              *
c     *    compute the 3x3 global -> crack front local rotation      *
c     *                                                              *
c     ****************************************************************
c
      subroutine dimrot( scoord, coord_map, u, dstmap, debug, iout )
c
      use main_data, only : trn, trnmat
      use j_data, only : num_front_nodes, front_order, front_nodes,
     &                   front_coords, front_node_displ,
     &                   domain_origin, domain_type, domain_rot,
     &                   crack_plane_normal, crack_front_tangent,
     &                   tangent_vector_defined
      implicit none
c
c                     parameter declarations
c
      integer :: iout, coord_map(*), dstmap(*)
      double precision :: scoord(*), u(*)
      logical :: debug
c
c                     local declarations
c
      integer :: node, node_id, caseno, frnint, status, i, ftype,
     &           nfnode
      double precision :: nx, ny, nz, xcrack(3), tvec(3),
     &                    length, x, y, z, unode, vnode, wnode
      double precision, parameter :: zero = 0.0d0, one = 1.0d0
c
      if( debug ) then
         write(iout,*) ' '
         write(iout,*) ' >>> entered make domain rotation'
         write(iout,*) '     > no. front nodes, interpolation: ',
     &                 num_front_nodes, front_order
      end if
c
c             gather coordinates and displacements of crack-front nodes
c
      do node = 1, num_front_nodes
         node_id                  = front_nodes(node)
         front_coords(1,node)     = scoord(coord_map(node_id))
         front_coords(2,node)     = scoord(coord_map(node_id)+1)
         front_coords(3,node)     = scoord(coord_map(node_id)+2)
         front_node_displ(1,node) = u(dstmap(node_id)+0)
         front_node_displ(2,node) = u(dstmap(node_id)+1)
         front_node_displ(3,node) = u(dstmap(node_id)+2)
      end do
c
      domain_rot = zero
      nx = crack_plane_normal(1)
      ny = crack_plane_normal(2)
      nz = crack_plane_normal(3)
c
c
c             if only one node on the crack front. we have the
c             simple 2-D case. rotation matrix determined only by
c             direction cosines. global-z and crack-z are assumed
c             to coincide.
c
      if( num_front_nodes .eq. 1 ) then
        nx = crack_plane_normal(1)
        ny = crack_plane_normal(2)
        domain_rot(1,1) = ny
        domain_rot(2,1) = nx
        domain_rot(1,2) = -one * nx
        domain_rot(2,2) = ny
        domain_rot(3,3) = one
        call dimrot_coord_displ
        if( .not. debug ) return
        write(iout,9100) nx, ny, domain_rot(1,1), domain_rot(1,2),
     &                   domain_rot(1,3), domain_rot(2,1),
     &                   domain_rot(2,2), domain_rot(2,3),
     &                   domain_rot(3,1), domain_rot(3,2),
     &                   domain_rot(3,3)
        return
      end if
c
c             more than one node on front specified. we have a
c             3-D situation.
c             get the position of the front-node to be used as
c             the domain origin. for domain types 'a', 'b' and
c             'c', (types 1, 2 and 3, respectively), take the
c             domain origin as the location where the crack-front
c             advance is maximum (where the q-value equals 1.0).
c             for domain type 'd' (uniform crack advance) take the
c             domain origin as the first front node.
c
      domain_origin = 0
      caseno        = 0
c
      select case( domain_type )
       case( 1 )
         domain_origin = 1
         if( num_front_nodes.eq.2 ) caseno = 1
         if( num_front_nodes.eq.3 ) caseno = 3
         if( num_front_nodes.eq.4 ) caseno = 4
       case( 2 )
         if( num_front_nodes .eq. 3 ) then
            domain_origin = 2
            caseno        = 2
         end if
         if( num_front_nodes .eq. 5 ) then
            domain_origin = 3
            caseno        = 5
         end if
         if( num_front_nodes .eq. 6 ) then
            domain_origin = 3
            caseno        = 5
         end if
       case( 3 )
         caseno = 1
         if( num_front_nodes.eq.2 ) domain_origin = 2
         if( num_front_nodes.eq.3 ) domain_origin = 3
       case( 4 )
         domain_origin = 1
      end select
c
c             based on the number of front nodes,
c             interpolation order, and q-function type construct a
c             unit tangent vector to the front that lies in the
c             crack plane (crack-z direction). crack-y is given
c             by dir cosines of normal to crack plane. construct
c             crack-x and final 3x3 structure-to-crack rotation
c             by cross-product.
c
c             for q-function type 'd' (=4), we just use the first
c             two nodes on the front to find the tangent vector.
c             this case is for a uniform crack advance along the
c             full front.
c
c             alternatively, the user may have defined components of
c             the crack front tangent vector rather than using the
c             computed value. this feature helps, for example, when the
c             mesh has the crack front intersecting a symmetry plane
c             at other than 90-degrees. A correct user-defined tangent
c             vector restores domain independent J-values. difrtn
c             takes care of this case as well.
c
c
      frnint = front_order
      ftype  = domain_type
      nfnode = num_front_nodes
c
      if( domain_type .eq. 4 ) then
c
c             get crack front (unit) tangent vector (crack-z) using
c             first 2 nodes listed on front.
c
        call dimrott( 2, front_nodes, 1, scoord, coord_map, tvec, 1,
     &                crack_front_tangent, tangent_vector_defined,
     &                iout, debug, status )
      else
        if ( nfnode .gt. 5 ) then
           write(iout,*) '>> too many nodes on front'
           call die_abort
           stop
        end if
        if( caseno .eq. 0 ) then
           write(iout,*) '>> invalid front node list/element type'
           write(iout,*) '>> combination'
           call die_abort
           stop
        end if
        if( debug ) then
           write(iout,*) '     > rotation by dir. cosines '
           write(iout,*) '        front interpolation: ', frnint
           write(iout,*) '        domain type:         ', ftype
        end if
c
c             get crack front (unit) tangent vector (crack-z).
c
        call dimrott( num_front_nodes, front_nodes, caseno, scoord,
     &               coord_map, tvec, domain_type,
     &               crack_front_tangent, tangent_vector_defined,
     &               iout, debug, status )
      end if
      if( status .eq. 1 ) then
         write(iout,9300)
         domain_rot(1,1) = one
         domain_rot(2,2) = one
         domain_rot(3,3) = one
         return
      end if
c
c             compute crack-x unit vector from z and y. fill in
c             global-to-crack rotation matrix.
c
      xcrack(1) = ny * tvec(3)  -  nz * tvec(2)
      xcrack(2) = nz * tvec(1)  -  nx * tvec(3)
      xcrack(3) = nx * tvec(2)  -  ny * tvec(1)
c
      length          = sqrt( xcrack(1)**2 + xcrack(2)**2 +
     &                        xcrack(3)**2 )
      xcrack(1)       = xcrack(1) / length
      xcrack(2)       = xcrack(2) / length
      xcrack(3)       = xcrack(3) / length
      domain_rot(1,1) = xcrack(1)
      domain_rot(1,2) = xcrack(2)
      domain_rot(1,3) = xcrack(3)
      domain_rot(2,1) = nx
      domain_rot(2,2) = ny
      domain_rot(2,3) = nz
      domain_rot(3,1) = tvec(1)
      domain_rot(3,2) = tvec(2)
      domain_rot(3,3) = tvec(3)
c
      call dimrot_coord_displ
c
      if ( .not. debug ) return
      write(iout,9200) xcrack, nx, ny, nz, tvec,
     &                 domain_rot(1,1), domain_rot(1,2),
     &                 domain_rot(1,3), domain_rot(2,1),
     &                 domain_rot(2,2), domain_rot(2,3),
     &                 domain_rot(3,1), domain_rot(3,2),
     &                 domain_rot(3,3)
      write(iout,9006) (node,(front_coords(i,node),
     &                     i=1,3), node=1,num_front_nodes )
      write(iout,9010) (node,(front_node_displ(i,node),
     &                     i=1,3), node=1,num_front_nodes )
c
      return
c
 9006 format (' >>>  front-node coordinates in crack reference frame',
     & //,5x,'node',15x,'x',19x,'y',19x,'z',15(/,i5,3f20.6))
 9010 format (' >>>  front-node displacements in crack reference frame',
     & //,5x,'node',15x,'u',19x,'v',19x,'w',15(/,i5,3f20.6))
 9100 format('     > debug for crack rotation by 2-d dir cosines.',
     & /,    '       nx -- ', f8.4,
     & /,    '       ny -- ', f8.4,
     & /,    '       rotation matrix -- ',
     &  3(/,10x,3f20.6)// )
 9200 format('      > debug for crack rotation by 3-d dir cosines.',
     & /,    '        xcrack -- ', 3f8.4,
     & /,    '        ycrack -- ', 3f8.4,
     & /,    '        zcrack -- ', 3f8.4
     & /,    '        rotation matrix -- ',
     &  3(/,10x,3f20.6)// )
 9300 format(  ' >>>>> ERROR: domain computation of global-to-'
     &       /,'       crack system rotation failed. identity',
     &       ' matrix assumed.',// )
c
      contains
c     ========
c
      subroutine dimrot_coord_displ
      implicit none
c
c              rotate coordinates of crack-front
c              nodes to local system.
c
      do node = 1, num_front_nodes
         x = domain_rot(1,1) * front_coords(1,node)
     &     + domain_rot(1,2) * front_coords(2,node)
     &     + domain_rot(1,3) * front_coords(3,node)
         y = domain_rot(2,1) * front_coords(1,node)
     &     + domain_rot(2,2) * front_coords(2,node)
     &     + domain_rot(2,3) * front_coords(3,node)
         z = domain_rot(3,1) * front_coords(1,node)
     &     + domain_rot(3,2) * front_coords(2,node)
     &     + domain_rot(3,3) * front_coords(3,node)
         front_coords(1,node) = x
         front_coords(2,node) = y
         front_coords(3,node) = z
c
c            if necessary, transform nodal displacement values
c            from constraint-compatible coordinates to global
c            coordinates. then rotate displacements from global
c            to local coordinates
c
         node_id = front_nodes(node)
         if( trn(node_id) ) then
            call di_trans_nodvals( front_node_displ(1,node),
     &           trnmat(node_id)%mat(1,1) )
         end if
c
         unode = domain_rot(1,1) * front_node_displ(1,node)
     &         + domain_rot(1,2) * front_node_displ(2,node)
     &         + domain_rot(1,3) * front_node_displ(3,node)
         vnode = domain_rot(2,1) * front_node_displ(1,node)
     &         + domain_rot(2,2) * front_node_displ(2,node)
     &         + domain_rot(2,3) * front_node_displ(3,node)
         wnode = domain_rot(3,1) * front_node_displ(1,node)
     &         + domain_rot(3,2) * front_node_displ(2,node)
     &         + domain_rot(3,3) * front_node_displ(3,node)
         front_node_displ(1,node) = unode
         front_node_displ(2,node) = vnode
         front_node_displ(3,node) = wnode
      end do
c
      return
      end subroutine dimrot_coord_displ
c
      end subroutine dimrot
c
c     ****************************************************************
c     *                                                              *
c     *                   support routine dimrott                    *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 2/17/2018 rhd              *
c     *                                                              *
c     *        unit tangent vector to front directed along crack     *
c     *        local-Z direction                                     *
c     *                                                              *
c     ****************************************************************
c
      subroutine dimrott( nfnode, fnodes, caseno, scoord,
     &                    coord_map, tvec, fnctyp, crack_front_tangent,
     &                    tangent_vector_defined, iout, debug, status )
      implicit none
c
c               parameter declarations
c
      integer :: nfnode, fnodes(*), caseno, fnctyp, iout, status,
     &           coord_map(*)
      double precision :: scoord(*), tvec(*), crack_front_tangent(*)
      logical :: tangent_vector_defined, debug
c
c               local declarations
c
      integer :: ii, jj, coord_loc, nfn
      double precision :: dsf(4), sf(4), tvec1(3), tvec2(3),
     &                    dx, dy, dz, xsi
      double precision, parameter :: zero = 0.0d0, half = 0.5d0,
     &                               one = 1.0d0
c
c             compute or load unit tangent vector to the crack
c             front according to the situation indicated by the
c             case number or the user-defined flag.
c             for cases 2,5,6 we average two tangents
c             computed for adjacent line segments. for cases 1,3,4 we
c             compute a single tangent at start or end of line segment.
c             fnodes is list of structure nodes defining the segment(s).
c
c             tangent_vector_defined = .T. --> the user specified the
c                                              tangent vector to use.
c                                              components are in
c                                              crack_front_tangent
c                                              skip computations.
c
c
c             case no.           situation
c             --------           ---------
c                1               linear; 2-nodes on front
c                2               linear; 3-nodes on front
c                3               quadratic; 3-nodes on front
c                4               cubic; 4-nodes on front
c                5               quadratic; 5-nodes on front
c                6               cubic; 7-nodes on front
c
c             fnctyp: 1-4, used for cases 3,4 to determine if
c                     tangent is needed for first or last node
c                     on front. = 1 (first node), =3 (last node).
c
c             status: returned value = 0 (ok) = 1 (bad data)
c
      if( tangent_vector_defined ) then
        dx = crack_front_tangent(1)
        dy = crack_front_tangent(2)
        dz = crack_front_tangent(3)
        call difrts( dx, dy, dz, tvec, status )
        if ( .not. debug ) return
        write(iout,9200)  tangent_vector_defined,  caseno, fnctyp,
     &                  status, tvec(1), tvec(2), tvec(3)
      end if
c
c             case 1,3,4: 2, 3 or four nodes on front but only a single
c                         element (linear, quadratic or cubic). compute
c                         tangent at first or last node (2,3 or 4). use
c                         utility routine to get derivatives of the 1-D
c                         isoparametric shape functions.
c
c             case 2,5,6: two linear, quadratic or cubic segments
c                         connecting 3, 5 or 6 front nodes.
c                         compute average tangent at center (common)
c                         node. compute two unit tangents, average
c                         components, then restore unit length.
      select case( caseno )
c
      case( 1, 3, 4 )
      if( fnctyp .eq. 1 ) then
         xsi = -one
       elseif( fnctyp .eq. 3 ) then
         xsi = one
       else
         write(iout,9100) fnctyp
         call die_gracefully
         stop
      end if
      call di1dsfc( xsi, dsf, sf, nfnode )
      dx = zero
      dy = zero
      dz = zero
      do ii = 1, nfnode
        coord_loc  = coord_map(fnodes(ii)) - 1
        dx = dx + dsf(ii)*scoord(coord_loc+1)
        dy = dy + dsf(ii)*scoord(coord_loc+2)
        dz = dz + dsf(ii)*scoord(coord_loc+3)
      end do
      call difrts( dx, dy, dz, tvec, status )
c
      case( 2, 5, 6 )
      nfn = 2
      if( caseno .eq. 5) nfn = 3
      if( caseno .eq. 6) nfn = 4
      xsi = one
      call di1dsfc( xsi, dsf, sf, nfn )
      dx = zero
      dy = zero
      dz = zero
      do ii = 1, nfn
        coord_loc = coord_map(fnodes(ii)) - 1
        dx = dx + dsf(ii)*scoord(coord_loc+1)
        dy = dy + dsf(ii)*scoord(coord_loc+2)
        dz = dz + dsf(ii)*scoord(coord_loc+3)
      end do
      call difrts( dx, dy, dz, tvec1, status )
      if( status .eq. 1 ) return
      xsi = -one
      call di1dsfc( xsi, dsf, sf, nfn )
      dx = zero
      dy = zero
      dz = zero
      jj = nfn - 1
      do ii = 1, nfn
        coord_loc = coord_map(fnodes(jj+ii)) - 1
        dx = dx + dsf(ii)*scoord(coord_loc+1)
        dy = dy + dsf(ii)*scoord(coord_loc+2)
        dz = dz + dsf(ii)*scoord(coord_loc+3)
      end do
      call difrts( dx, dy, dz, tvec2, status )
      if( status .eq. 1 ) return
      dx = (tvec1(1) + tvec2(1)) * half
      dy = (tvec1(2) + tvec2(2)) * half
      dz = (tvec1(3) + tvec2(3)) * half
      call difrts( dx, dy, dz, tvec, status )
c
      end select
c
c             all done.
c
      if ( .not. debug ) return
      write(iout,9200)  tangent_vector_defined,  caseno, fnctyp,
     &                  status, tvec(1), tvec(2), tvec(3)
      return
c
 9100 format('>>>>> SYSTEM ERROR:in difrtn. fnctyp has illegal',
     &  /,   '                   value: ',i10,
     &  /,   '                   job terminated.')
 9200 format('      > leaving dimrott...',
     & /,20x,'tangent_vector_defined: ',l1,
     & /,20x,'caseno:       ',i10,
     & /,20x,'domain_type:  ',i10,
     & /,20x,'status:       ',i10,
     & /,20x,'tangent vec:  ',3f10.6)
c
      end
c *******************************************************************
c *                                                                 *
c *   support routine: called by difrts                             *
c *                                                                 *
c *******************************************************************
c
c
      subroutine difrts( dx, dy, dz, tvec, status )
      implicit none
c
      integer :: status
      double precision :: dx, dy, dz, tvec(3)
c
      double precision :: length
c
      length = sqrt( dx*dx + dy*dy + dz*dz )
      status = 1
      if( length .lt. 1.0d-07 ) return
      tvec(1) = dx / length
      tvec(2) = dy / length
      tvec(3) = dz / length
      status = 0
c
      return
      end
c *******************************************************************
c *                                                                 *
c *   subroutine di1dsfc ---- calculates shape function values and  *
c *                           derivatives for 1-dimension           *
c *                                                                 *
c *******************************************************************
c
      subroutine di1dsfc( xsi, dsf, sf, nlnode )
c
c              parameter declarations
c
      implicit none
c
      integer :: nlnode
      double precision  :: xsi, dsf(*), sf(*)
c
c              local declarations
c
      double precision :: xsisqr
      double precision, parameter :: half = 0.5d0, one = 1.0d0,
     &                               two = 2.0d0
c
      select case( nlnode )
c
        case(1, 2)  ! linear
        sf(1)  = half * ( one - xsi )
        sf(2)  = half * ( one + xsi )
        dsf(1) = -half
        dsf(2) =  half
c
        case( 3 )  ! quadratic
        xsisqr = xsi * xsi
        sf(1)  = half * ( xsisqr - xsi )
        sf(2)  = one - xsisqr
        sf(3)  = half * ( xsisqr + xsi )
        dsf(1) = xsi - half
        dsf(2) = -two*xsi
        dsf(3) = xsi + half
c
      end select
c
      return
      end

