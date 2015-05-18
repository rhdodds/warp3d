c *******************************************************************
c *                                                                 *
c *    subroutine domain_frnt_area --- calculates the area under    *
c *    the q-function along the crack front                         *
c *                                                                 *
c *******************************************************************
c
      subroutine difrar( scoord, coord_map, iout,
     &                   debug, error, bad_domain )
      use j_data
      implicit integer(a-z)
c
c             parameter declarations
c
#dbl      double precision
#sgl      real
     &  scoord(*)
      integer coord_map(*)
      logical  debug, error, bad_domain
c
c             local declarations
c
#dbl      double precision
#sgl      real
     & sf(3), dsf(3), xsi(3), weight(3),
     & zero, x1, y1, z1, x2, y2, z2, half, dl, gp_loc,
     & dx, dy, dz, qvalue, weight13, weight2, x, y, z
      data zero, half, gp_loc, weight13, weight2
     &   / 0.0, 0.5, 0.774596669241483, 0.555555555555556,
     &     0.888888888888889 /
c
c             if no nodes on the crack front return with error.
c             check for the special case of only one node on the
c             crack front (i.e. 2-dimensions )
c
      if ( debug ) then
         write(iout,*) ' >>> entered compute front area'
      end if
c
      error = .true.
      if ( num_front_nodes .eq. 1 )  then
        front_q_area = q_values(front_nodes(1))
        front_length = zero
        error        = .false.
        if ( debug ) write(iout,9000) front_q_area, front_length 
        if ( front_q_area .le. zero ) then
          bad_domain = .true.
          write(iout,9020)
          return
        endif
      endif
c
c             front order = 1: linear q variation along front
c             integrate simple linear q-function variation
c             by trapezoidal rule. loop over each segment and
c             sum contribution.
c
      if ( front_order .eq. 1 ) then
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
        if ( debug ) write(iout,9000) front_q_area, front_length 
        error = .false.
        if ( front_q_area .le. zero ) then
          bad_domain = .true.
          write(iout,9020)
          return
        endif
        write(iout,9000) front_q_area, front_length 
        return
      end if
c
c             quadratic elements on front
c             check to see if the number of crack front nodes
c             are valid for the interpolation order. 
c
      ngpnts = 3
      nlnode = 3
      if ( debug )  write( iout, 9010 ) ngpnts, nlnode
c
      nsegmn = ( num_front_nodes - 1 )/( nlnode - 1 )
      if ( mod( num_front_nodes - 1, nlnode - 1 ) .ne. 0 ) return
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
      if ( front_q_area .le. zero ) then
          bad_domain = .true.
          write(iout,9020)
          return
      endif
c
      write(iout,9000) front_q_area, front_length
c
      error = .false.
      return
c
c
 9000 format(/,' area under q along front:       ',e13.6,
     & /,      ' length of crack front segment:  ',e13.6)
 9010 format( '      ngpnts: ', i6, '     nlnode: ', i6 ) 
 9020 format(//,'>>> Area under q-function along the crack front',
     &          ' is not >0.',/,
     &          '    J-value cannot be computed. Check that q-values',
     &          ' are properly specified',
     &        /,'    for the front nodes. Domain integration procedure',
     &          ' terminated....',//)
c
      end

