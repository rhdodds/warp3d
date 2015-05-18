c     ****************************************************************
c     *                                                              *
c     *                      subroutine ounds1                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/11/04 (rhd)              *
c     *                                                              *
c     *     this subroutine extrapolates the element stresses or     *
c     *     strains from the integration points used for stiffness   *
c     *     strains-stresses to the nodes of the element for a       *
c     *     block of similar solid elements.                         *
c     *                                                              *
c     ****************************************************************
c
c           
      subroutine ounds1( span, type, order, nnode, ngp, stress,
     &                   num_short_stress, num_short_strain )
      use elblk_data, only : elestr
      implicit integer (a-z)
$add param_def
      logical stress
c
c                       local declarations
c
#dbl      double precision
#sgl      real
     &     temstr(mxvl,mxoupr,mxndel), lg(mxgp), zero,
     &     dum_vec(1), rngpts, xi, eta, zeta
      dimension idumvec(1)
      logical lagrangian_extrap, hex_elem, tet10, tet4
      data zero /0.0/
c
c                       set number of values to define at element nodes.
c                       determine if hex element and that the order
c                       of integration allows lagrangian extrapolation
c                       to elemetn nodes. 
c                       
      num_vals = num_short_strain
      if ( stress ) num_vals = num_short_stress
c
      lagrangian_extrap = .false.
      hex_elem = .false.
      hex_elem = type .ge. 1 .and. type .le. 5
      tet4 = type .eq. 13
      tet10 = type .eq. 6
      if( hex_elem ) then
         if( type .eq. 2 ) then
            lagrangian_extrap = .true.
         else if( order .eq. 1 .or. order .eq. 8 ) then
            lagrangian_extrap = .true.
         end if
      end if
c
      if ( lagrangian_extrap ) go to 1000
c
      if ( tet10 ) then
c
c                       for 10 node tet, we set corner node values equal
c                       to the nearest integration point. then average
c                       corner node values to define mid-side node values.
c               
        shift = 0
        if ( ngp .eq. 5 ) shift = 1
        do elnod = 1, 4
          tetpt = shift + elnod
          temstr(1:span,1:num_vals,elnod) =
     &                 elestr(1:span,1:num_vals,tetpt)
        end do
        call outetavg( span, 5, 1, 2, temstr, mxvl, mxoupr, num_vals )
        call outetavg( span, 6, 2, 3, temstr, mxvl, mxoupr, num_vals )
        call outetavg( span, 7, 3, 1, temstr, mxvl, mxoupr, num_vals )
        call outetavg( span, 8, 1, 4, temstr, mxvl, mxoupr, num_vals )
        call outetavg( span, 9, 4, 2, temstr, mxvl, mxoupr, num_vals )
        call outetavg( span, 10, 4, 3, temstr, mxvl, mxoupr, num_vals )
        elestr(1:span,1:num_vals,1:10) = temstr(1:span,1:num_vals,1:10)
        return
      end if
c
c                       this code handles elements for which we cannot
c                       devise an acceptable extrapolation procedure to
c                       the element nodes. for example, 14 pt integration
c                       rule for hexes, tet4 elements, etc. we use the
c                       average of all integration points at each node.
c
      temstr(1:span,1:num_vals,1) = zero
c
      do j = 1, ngp
       do k = 1, num_vals
         do i = 1, span 
            temstr(i,k,1)= temstr(i,k,1) + elestr(i,k,j)
         end do
       end do
      end do
c        
c                       replace the element stresses/strains at 
c                       the gauss points with those at the nodes.
c
      rngpts = ngp
      do j = 1, nnode
       elestr(1:span,1:num_vals,j) =
     &                 temstr(1:span,1:num_vals,1) / rngpts
      end do
      return
c
c                       lagrangian extrapolation from
c                       Gauss point values to element nodes for
c                       hex elements
c
 1000 continue
c                        
      temstr(1:span,1:num_vals,1:nnode) = zero
      do elnod = 1, nnode
c
c                       find the lagrangian shape functions for
c                       each gauss point at the current element node.
c       
         call ndpts1( idumvec, 0, dum_vec, type, elnod,
     &                xi, eta, zeta ) 
         call oulgf( type, xi, eta, zeta, lg, order ) 
c
c                       extrapolate stresses/strains to current
c                       nodal values using the the lagrangian
c                       shape functions.        
c          
         do j = 1, ngp
          do k = 1, num_vals
            do i = 1, span 
             temstr(i,k,elnod)= temstr(i,k,elnod) + elestr(i,k,j)*lg(j)
            end do
          end do
         end do
c
      end do
c        
c                       replace the element stresses/strains at 
c                       the gauss points with those at the nodes.
c
      elestr(1:span,1:num_vals,1:nnode) =
     &                    temstr(1:span,1:num_vals,1:nnode)
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine outetavg                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/26/01                   *
c     *                                                              *
c     *     average corner node strain-stress values to define       *
c     *     values at a mid-side node                                *
c     *                                                              *
c     ****************************************************************
c
c           
      subroutine outetavg( span, nm, n1, n2, temstr, mxvl, mxoupr,
     &                     num_vals )
      implicit integer (a-z)
#dbl      double precision
#sgl      real
     &     temstr(mxvl,mxoupr,*), half
#dbl      data half / 0.5d00 /
#sgl      data half / 0.5 /
c
      do k = 1, num_vals
       temstr(1:span,k,nm) = half * ( temstr(1:span,k,n1) +
     &                                temstr(1:span,k,n2 ) )
      end do
c
      return
      end
