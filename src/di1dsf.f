c *******************************************************************
c *                                                                 *
c *    subroutine di1dsf ---- calculates shape function values and  *
c *                           derivatives for 1-dimension           *
c *                                                                 *
c *******************************************************************
c
      subroutine di1dsf( xsi, dsf, sf, nlnode )
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
