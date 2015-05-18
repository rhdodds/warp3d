c     ****************************************************************
c     *                                                              *
c     *                      subroutine gpmas1                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/20/2015 rhd             *
c     *                                                              *
c     *     this subroutine computes the lumped mass matrix terms    *
c     *     for elements in the block. uses undeformed coordinates.  *
c     *     compute contribution just for the specified gauss point. *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine gpmas1( span, felem, type, order, gpn, nnode,
     &                   volume, mel, totdof, fgm_props, rho, emass,
     &                   iout, rho_fgm_flags, axisymm_elem,
     &                   implemented_elem, ce_block, rho_block )
      implicit none
$add param_def
c
c                       parameter declarations
c
      integer :: span, felem, type, order, gpn, nnode,
     &                 totdof, iout
#dbl      double precision ::
#sgl      real ::
     &  volume(*), mel(totdof,*), rho(*), emass(*), rho_fgm_flags(*),
     &  ce_block(mxvl,*), rho_block(mxndel,*)
      logical :: fgm_props, axisymm_elem, implemented_elem
c
c                       locals
c
      integer :: i, j
#dbl      double precision ::
#sgl      real ::
     &   shape(mxndel), radius(mxvl), dummy(mxvl,3,3),
     &   jac(mxvl,ndim,ndim),  gama(mxvl,ndim,ndim),
     &   xi, eta ,zeta, wfactor, dj(mxvl),
     &   nxi(mxndel), neta(mxndel), nzeta(mxndel),
     &   one
c
      data one / 1.0d00 /
c                       
c                       message for unsupported elements                       
c
      if( .not. implemented_elem ) then
        write(iout,9500) type,span,felem,order,gpn,nnode,totdof
        do j = 1, nnode
          do i = 1, span
            mel(j,i) = one
          end do
        end do
        do i = 1, span
          volume(i) = one
        end do
        return
      end if
c
c                       compute the shape functions and the 
c                       shape function derivatives.
c
      call getgpts( type, order, gpn, xi, eta, zeta, wfactor )
      call derivs(  type, xi, eta, zeta, nxi, neta, nzeta )
      call shapef(  type, xi, eta, zeta, shape )
c
c                       for axisymmetric elements, compute the radius
c                       to the current gauss point for each element
c                       in the block.
c
      if( axisymm_elem ) then
       call get_radius( radius, nnode, span, shape, ce_block, mxvl )
      end if
c
c                       compute the determinate of the coordinate
c                       jacobian (use initial coordinates).
c
      call jacob1( type, span, felem, gpn, jac, dj, gama, dummy,
     &             nxi, neta, nzeta, ce_block, nnode )
c
c                       some elements may have the mass density
c                       defined using the fgm capabilities, i.e,
c                       thru nodal values. interpolate density
c                       at this gauss point for all elements in
c                       the block.
c
      if ( fgm_props )
     &    call set_fgm_solid_props_for_block(
     &                span, felem, type, gpn, nnode, rho, shape, 
     &                rho_block, 1, rho_fgm_flags )
c                    
c                       compute the gauss point contribution to the
c                       the nodal mass term.
c
      call mass_sum( type, totdof, nnode, mxvl, span, radius,
     &               shape, dj, wfactor, mel, volume, rho, emass )
c
c
      return
c
9500    format(1x,//,'>>>>>  in gpmas1, Mass matrix not yet defined',/,
     &               '       for element type = ',i3,/,
     &               '   Skipping mass matrix computation, set',/,
     &               '   unit values for now.',//,
     &               '     span =',i6,/,
     &               '    felem =',i6,/,
     &               '    order =',i6,/,
     &               '      gpn =',i6,/,
     &               '    nnode =',i6,/,
     &               '   totdof =',i6)
      end
c     ****************************************************************
c     *                                                              *
c     *                     subroutine mass_sum                      *
c     *                                                              *
c     *                       written by: rau                        *
c     *                    last modified: 07/11/00                   *
c     *                                                              *
c     *     this routine computes the current gauss point            *
c     *     contribution to the terms of the consistent              *
c     *     diagonal mass matrix for one (any) dof. the              *
c     *     matrix will be scaled later to include mass              *
c     *     density.                                                 *
c     *                                                              *
c     ****************************************************************
c
c         for triangular elements, the correct area is given by
c         0.5*|J|. for tetrahedral elements, the correct volume is
c         given by (1/6)*|J|.
c
c         For axisymmetric elements, include the 2*pi*radius scalar term
c         in the element stiffness summation.
c
c
c         Variables:
c
c         elem_type = integer flag for the element type
c         totdof    = total number of degrees of freedom, 
c                     dimension for mel()
c         mxvl      = maximum number of elements in a block
c         span      = number of elements in the block
c         radius()  = radius to the current Gauss point for each element
c                     in the current block.
c         shape()   = shape function array evaluated at current Gauss point
c         dj()      = Jacobian determinant for each element in the block
c         w         = Gauss integration weight
c         mel()     = mass contribution at a node for each element
c                     in the block at the current Gauss integration point
c         volume()  = volume contribution at a node for each element in the
c                     block at the current Gauss integration point
c
c
      subroutine mass_sum( elem_type, totdof, nnode, mxvl, span, radius,
     &                     shape, dj, w, mel, volume, rho, emass )
      implicit none
c
c               parameter declarations
c
      integer :: elem_type, totdof, nnode, mxvl, span
#dbl      double precision ::
#sgl      real ::
     &         radius(*), 
     &         shape(*), 
     &         dj(*), 
     &         w, 
     &         mel(totdof,*), 
     &         volume(*),
     &         rho(*),
     &         emass(*)
c
c               local variables 
c
      integer :: i,j
#dbl      double precision ::
#sgl      real ::
     &         two, 
     &         half, 
     &         pi, 
     &         scalar,
     &         weight
c
      data two, half, pi 
     & / 2.d0, 0.5d0, 3.14159265358979323846d0 /
c
c               the element type determines the calculations 
c               performed. compute the current mass and volume
c               contributions for the different element types.
c
      go to ( 100, 100, 100, 100, 100, 600, 100, 800, 100,
     &        1000, 1100, 100, 600 ), elem_type
c
c               element numbers 1-5, 7, 9, 12: hex elements
c
 100  continue
c
      do j = 1, nnode                                                
         do i = 1, span                                               
           mel(j,i) = mel(j,i) + shape(j)*shape(j)*dj(i)*w*rho(i)            
         end do                                                       
      end do                                                         
c                                                                      
      do i = 1, span                                                 
         volume(i) = volume(i) + dj(i)*w
         emass(i) = emass(i) + dj(i)*w*rho(i)
      end do      
c
      return
c
c               element number 6, 13: tet10, tet4
c                                                   
 600   continue
c
       scalar = 1.0D0/6.0D0
       weight = scalar*w
       do j = 1, nnode
          do i = 1, span
            mel(j,i) = mel(j,i) + shape(j)*shape(j)*dj(i)*weight*rho(i)
          end do
       end do
c
       do i = 1, span
          volume(i) = volume(i) + dj(i)*weight
          emass(i) = emass(i) + dj(i)*weight*rho(i)
       end do
c     
       return
c
c               element number 8: tri6
c
 800   continue
c
       scalar = half
       weight = scalar*w
       do j = 1, nnode
          do i = 1, span
            mel(j,i) = mel(j,i) + shape(j)*shape(j)*dj(i)*weight*rho(i)
          end do
       end do
c
       do i = 1, span
          volume(i) = volume(i) + dj(i)*weight
          emass(i) = emass(i) + dj(i)*weight*rho(i)
       end do
c     
       return
c
c               element numbers 10: axisymmetric quad 
c
 1000  continue
c
c               axisymmetric elements, multiply by 2*pi*radius
c
      scalar = two * pi
       do j = 1, nnode
          do i = 1, span
            mel(j,i) = mel(j,i) + 
     &                 scalar*radius(i)*shape(j)*shape(j)*dj(i)*w*rho(i)
          end do
       end do
       do i = 1, span
          volume(i) = volume(i) + scalar*radius(i)*dj(i)*w
          emass(i) = emass(i) + dj(i)*w*rho(i)
       end do
c
       return
c
c               element number 11: axisymmetric tri 
c               note that for the axisymmetric triangle 2*pi*0.5 = pi.
c                                                   
 1100  continue
c
       scalar = pi
       do j = 1, nnode
          do i = 1, span
            mel(j,i) = mel(j,i) + 
     &                 scalar*radius(i)*shape(j)*shape(j)*dj(i)*w*rho(i)
          end do
       end do
       do i = 1, span
          volume(i) = volume(i) + scalar*radius(i)*dj(i)*w
          emass(i)  = emass(i) + dj(i)*w*rho(i)
       end do
c
       return
       end








