c *******************************************************************
c *                                                                 *
c *   domain_rotation                                               *
c *                                                                 *
c *******************************************************************
c
c
      subroutine dimrot( scoord, coord_map, debug, iout )
      use main_data
      use j_data
      implicit integer (a-z)
      include 'common.main'
c
c                     parameter declarations
c
      integer iout
      dimension  scoord(*), coord_map(*)
      double precision
     & scoord
      logical debug
c
c                     local declarations
c
      integer node, node_id, caseno
      double precision
     & nx, ny, nz, one, zero, xcrack(3), tvec(3), length, x, y, z,
     & unode, vnode, wnode
      data zero, one / 0.0d00, 1.0d00 /
c
      if ( debug ) then
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
      if( domain_type.eq.1 ) then
         domain_origin = 1
         if( num_front_nodes.eq.2 ) caseno = 1
         if( num_front_nodes.eq.3 ) caseno = 3
         if( num_front_nodes.eq.4 ) caseno = 4
      end if
      if( domain_type.eq.2 ) then
         if( num_front_nodes.eq.3 ) then
            domain_origin = 2
            caseno        = 2
         end if
         if( num_front_nodes.eq.5 ) then
            domain_origin = 3
            caseno        = 5
         end if
         if( num_front_nodes.eq.6 ) then
            domain_origin = 3
            caseno        = 5
         end if
      end if
      if( domain_type.eq.3 ) then
         caseno = 1
         if( num_front_nodes.eq.2 ) domain_origin = 2
         if( num_front_nodes.eq.3 ) domain_origin = 3
      end if
      if( domain_type.eq.4 ) domain_origin = 1
c
c             if only one node on the crack front. we have the
c             simple 2-D case. rotation matrix determined only by
c             direction cosines. global-z and crack-z are assumed
c             to coincide.
c
      domain_rot(1,1) = zero
      domain_rot(2,1) = zero
      domain_rot(3,1) = zero
      domain_rot(1,2) = zero
      domain_rot(2,2) = zero
      domain_rot(3,2) = zero
      domain_rot(1,3) = zero
      domain_rot(2,3) = zero
      domain_rot(3,3) = zero
c
      nx = crack_plane_normal(1)
      ny = crack_plane_normal(2)
      nz = crack_plane_normal(3)
c
      if ( num_front_nodes .eq. 1 ) then
        nx = crack_plane_normal(1)
        ny = crack_plane_normal(2)
        domain_rot(1,1) = ny
        domain_rot(2,1) = nx
        domain_rot(1,2) = -one * nx
        domain_rot(2,2) = ny
        domain_rot(3,3) = one
        if ( .not. debug ) return
        write(iout,9100) nx, ny, domain_rot(1,1), domain_rot(1,2),
     &                   domain_rot(1,3), domain_rot(2,1),
     &                   domain_rot(2,2), domain_rot(2,3),
     &                   domain_rot(3,1), domain_rot(3,2),
     &                   domain_rot(3,3)
        return
      end if
c
c             more than one node on front specified. we have a
c             3-D situation. based on the number of front nodes,
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
c             the new feature enables the user to specify components of
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
      if ( domain_type .eq. 4 ) then
c
c             get crack front (unit) tangent vector (crack-z) using
c             first 2 nodes listed on front.
c
        call difrtn( 2, front_nodes, 1, scoord, coord_map, tvec, 1,
     &               crack_front_tangent, tangent_vector_defined,
     &               iout, debug, status )
      else
        if ( nfnode .gt. 5 ) then
           write(iout,*) '>> too many nodes on front'
           call die_abort
           stop
        end if
        if ( caseno .eq. 0 ) then
           write(iout,*) '>> invalid front node list/element type'
           write(iout,*) '>> combination'
           call die_abort
           stop
        end if
        if ( debug ) then
           write(iout,*) '     > rotation by dir. cosines '
           write(iout,*) '        front interpolation: ', frnint
           write(iout,*) '        domain type:         ', ftype
        end if
c
c             get crack front (unit) tangent vector (crack-z).
c
        call difrtn( num_front_nodes, front_nodes, caseno, scoord,
     &               coord_map, tvec, domain_type,
     &               crack_front_tangent, tangent_vector_defined,
     &               iout, debug, status )
      end if
      if ( status .eq. 1 ) then
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
c
c              rotate coordinates of crack-front
c              nodes to local system.
c
      do node = 1,num_front_nodes
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
      end
