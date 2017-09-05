c     ****************************************************************
c     *                                                              *
c     *                      subroutine ounds1                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 8/12/2017 rhd              *
c     *                                                              *
c     *     elestr on entry has element strains or stresses at       *
c     *     all integration points for all elements in block.        *
c     *     replace elestr with values at element nodes              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ounds1( span, elem_type, int_order, num_enode, ngp,
     &                   do_stresses, num_short_stress,
     &                   num_short_strain )
      use elblk_data, only : elestr ! mxvl,mxoupr,mxoupt
      implicit none
      include 'param_def'
      integer :: span, elem_type, int_order, num_enode, ngp,
     &           num_short_stress, num_short_strain
      logical :: do_stresses
c
c                       local declarations
c
      integer :: idumvec(1), num_vals, shift, elnod, tetpt, j, k, i
      logical :: lagrangian_extrap, tet10, tet4, threed_solid_elem,
     &           hex_elem, wedge_elem, tet_elem, twod_elem, quad_elem,
     &           triangle_elem, axisymm_elem, cohesive_elem, bar_elem,
     &           link_elem
      double precision :: temstr(mxvl,mxoupr,mxndel), lg(mxgp),
     &                    dum_vec(1), rngpts, xi, eta, zeta
      double precision, parameter :: zero = 0.0d0
c
c                       set number of values to define at element nodes.
c                       determine if hex/tet element and that the order
c                       of integration allows lagrangian extrapolation
c                       to element nodes.
c
      num_vals = num_short_strain
      if( do_stresses ) num_vals = num_short_stress
c
      call set_element_type( elem_type, threed_solid_elem,
     &                       hex_elem, wedge_elem, tet_elem,
     &                       twod_elem, quad_elem,
     &                       triangle_elem, axisymm_elem,
     &                       cohesive_elem, bar_elem, link_elem )

      tet4  = tet_elem .and. num_enode .eq. 4
      tet10 = tet_elem .and. num_enode .eq. 10
      lagrangian_extrap = .false.
      if( hex_elem ) then
         if( elem_type .eq. 2 ) then
            lagrangian_extrap = .true.
         else if( int_order .eq. 1 .or. int_order .eq. 8 ) then
            lagrangian_extrap = .true.
         end if
      end if
c
      if( lagrangian_extrap ) then
          call ounds1_lagrag
          return
      end if
c
      if( tet10 ) then
          call ounds1_tet10
          return
      end if
c
c              this code handles elements for which we cannot
c              devise an acceptable extrapolation procedure to
c              the element nodes. for example, 14 pt integration
c              rule for hexes, tet4 elements, etc. we use the
c              average of all integration points at each node.
c
      temstr = zero  ! all terms
c
      do j = 1, ngp
       do k = 1, num_vals
         do i = 1, span
            temstr(i,k,1)= temstr(i,k,1) + elestr(i,k,j)
         end do
       end do
      end do
c
      rngpts = ngp
      elestr = zero ! all terms
c
      do j = 1, num_enode
       elestr(1:span,1:num_vals,j) =
     &                 temstr(1:span,1:num_vals,1) / rngpts
      end do
c
      return

      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *                  subroutine ounds1_tet10                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/3/2016 rhd               *
c     *                                                              *
c     *     conjure up node values from integration point values     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ounds1_tet10
      implicit none
c
c              for 10 node tet, we set corner node values equal
c              to the nearest integration point. then average
c              corner node values to define mid-side node values.
c
c              int order   element node, int point pairs
c               4 pts      (1,1), (2,2), (3,3), (4,4)
c               5 pts      (1,2), (2,3), (3,4), (4,5)
c
      temstr = zero
      shift = 0
      if( ngp .eq. 5 ) shift = 1
      do elnod = 1, 4
        tetpt = shift + elnod
        temstr(1:span,1:num_vals,elnod) =
     &         elestr(1:span,1:num_vals,tetpt)
      end do
c
      call ounds1_avg(  5, 1, 2 )
      call ounds1_avg(  6, 2, 3 )
      call ounds1_avg(  7, 3, 1 )
      call ounds1_avg(  8, 1, 4 )
      call ounds1_avg(  9, 4, 2 )
      call ounds1_avg( 10, 4, 3 )
c
      elestr = zero
      elestr(1:span,1:num_vals,1:10) = temstr(1:span,1:num_vals,1:10)
c
      return
c
      end subroutine ounds1_tet10
c

c     ****************************************************************
c     *                                                              *
c     *                  subroutine ounds1_lagrag                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/3/2016 rhd               *
c     *                                                              *
c     *     lagrangian extrapolation from Gauss point values to      *
c     *     element nodes for hex elements                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ounds1_lagrag
      implicit none
c
      temstr = zero
c
      do elnod = 1, num_enode
c
c              find the lagrangian shape functions for
c              each gauss point at the current element node.
c
         call ndpts1( idumvec, 0, dum_vec, elem_type, elnod,
     &                xi, eta, zeta )
         call oulgf( elem_type, xi, eta, zeta, lg, int_order )
c
c              extrapolate stresses/strains to
c              nodal values using the the lagrangian
c              shape functions.
c
         do j = 1, ngp
          do k = 1, num_vals
            do i = 1, span
             temstr(i,k,elnod) = temstr(i,k,elnod) +
     &                           elestr(i,k,j)*lg(j)
            end do
          end do
         end do
c
      end do ! over element nodes
c
c              replace the element stresses/strains at
c              the gauss points with those at the nodes.
c
      elestr = zero
      elestr(1:span,1:num_vals,1:num_enode) =
     &                    temstr(1:span,1:num_vals,1:num_enode)
c
      return
      end subroutine ounds1_lagrag

c     ****************************************************************
c     *                                                              *
c     *                      subroutine outetavg                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/3/2016 rhd               *
c     *                                                              *
c     *     average corner node strain-stress values to define       *
c     *     values at a mid-side node                                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ounds1_avg( nm, n1, n2)
      implicit none
      double precision :: half
      integer :: nm, n1, n2
      data half / 0.5d00 /
c
      do k = 1, num_vals
       temstr(1:span,k,nm) = half * ( temstr(1:span,k,n1) +
     &                                temstr(1:span,k,n2 ) )
      end do
c
      return
      end subroutine ounds1_avg
c
      end subroutine ounds1
