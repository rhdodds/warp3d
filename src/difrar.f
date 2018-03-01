c     ****************************************************************
c     *                                                              *
c     *                subroutine di_front_q_area                    *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 2/17/2018 rhd              *
c     *                                                              *
c     * compute area under q-function along front for this domain    *
c     *                                                              *
c     ****************************************************************
c
      subroutine di_front_q_area( scoord, coord_map, iout,
     &                            debug, error, bad_domain )
      use j_data, only : num_front_nodes, q_values, front_q_area,
     &                   front_order, front_length, front_nodes
      implicit none
c
c             parameter declarations
c
      integer :: iout, coord_map(*)
      double precision :: scoord(*)
      logical :: debug, error, bad_domain
c
c             local declarations
c
      integer :: i, inode, jnode, coord_loc, ngpnts, nlnode, nsegmn,
     &           isgmnt, iloc, igpnt, lnode, snode
      double precision :: sf(3), dsf(3), xsi(3), weight(3),
     &                    x1, y1, z1, x2, y2, z2, dl, dx, dy, dz,
     &                    qvalue, x, y, z
      double precision, parameter :: zero = 0.0d0, half = 0.5d0,
     &  gp_loc = sqrt(0.6d0), weight13 = 5.d0/9.d0, weight2 = 8.d0/9.d0
c
c             if no nodes on the crack front return with error.
c             check for the special case of only one node on the
c             crack front (i.e. 2-dimensions )
c
      if( debug ) write(iout,*) ' >>> entered compute front area'
      error = .true.
c
      if( num_front_nodes .eq. 1 )  then
        front_q_area = q_values(front_nodes(1))
        front_length = zero
        error        = .false.
        if( front_q_area .le. zero ) then
          bad_domain = .true.
          return
        end if
      end if

      select case( front_order ) ! linear, quadratic for now
c
      case( 1 )  ! linear elements at front this domain
        front_q_area  = zero
        front_length  = zero
        do i = 1, num_front_nodes - 1
         inode      = front_nodes(i)
         jnode      = front_nodes(i+1)
         coord_loc  = coord_map(inode) - 1
         x1         = scoord(coord_loc+1)
         y1         = scoord(coord_loc+2)
         z1         = scoord(coord_loc+3)
         coord_loc  = coord_map(jnode) - 1
         x2         = scoord(coord_loc+1)
         y2         = scoord(coord_loc+2)
         z2         = scoord(coord_loc+3)
         dl         = sqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )
         front_q_area  = front_q_area + dl*half*(q_values(inode) +
     &                   q_values(jnode))
         front_length = front_length + dl
        end do
c
c             front order = 2: quadratic elements on front
c             check to see if the number of crack front nodes
c             are valid for the interpolation order.
c
      case( 2 )  ! quadratic front elements this domain
        ngpnts = 3
        nlnode = 3
        nsegmn = ( num_front_nodes - 1 )/( nlnode - 1 )
        if( mod( num_front_nodes - 1, nlnode - 1 ) .ne. 0 ) then
          write(iout,9000)
          call die_gracefully
        end if
c
c             loop over the number of gauss points and calculate
c             the gaussian weighted integrand, accumulating sum.
c
        front_q_area = zero
        front_length = zero
        xsi(1)       = -gp_loc
        xsi(2)       = zero
        xsi(3)       = gp_loc
        weight(1)    = weight13
        weight(2)    = weight2
        weight(3)    = weight13
        do isgmnt = 1, nsegmn
         iloc  = ( isgmnt - 1 )*( nlnode - 1 )
         do igpnt = 1, ngpnts
           call di1dsf( xsi(igpnt), dsf, sf, nlnode )
           dx     = zero
           dy     = zero
           dz     = zero
           qvalue = zero
           do lnode = 1, nlnode
            snode      = front_nodes(iloc+lnode)
            coord_loc  = coord_map(snode) - 1
            x          = scoord(coord_loc+1)
            y          = scoord(coord_loc+2)
            z          = scoord(coord_loc+3)
            dx         = dx + dsf(lnode)*x
            dy         = dy + dsf(lnode)*y
            dz         = dz + dsf(lnode)*z
            qvalue     = qvalue  +  sf(lnode)*q_values(snode)
          end do
          dl           = sqrt( dx**2 + dy**2 + dz**2 )
          front_q_area = front_q_area + qvalue*dl*weight(igpnt)
          front_length = front_length + dl*weight(igpnt)
         end do
        end do
c
      case default
        write(iout,9010)
        call die_gracefully
c
      end select
c
      error = .false.
      bad_domain = front_q_area <= zero
      return
c
 9000 format('>> FATAL INTERNAL ERROR: invalid number of front nodes',
     &   /,  '                         for domain integral. stop.' )
 9010 format('>> FATAL INTERNAL ERROR: invalid front interpolation',
     &   /,  '                         for domain integral. stop.' )
c
      end

