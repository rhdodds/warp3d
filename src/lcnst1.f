c     ****************************************************************
c     *                                                              *
c     *  subroutines to compute the linear-elastic constitutive      *
c     *  matrices for all material models in WARP3D.                 *
c     *                                                              *
c     ****************************************************************
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine lcnst1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 10/25/89                   *
c     *                                 : 05/27/99 aroy              *
c     *                                 : 11/21/2010, RHD            *
c     *                                                              *
c     *      this subroutine determines the linear constituitive     *
c     *      matrices for 3D solid elements                          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine lcnst1( span, cep, nu, e, dj, w )
      implicit integer (a-z)
$add param_def
#dbl      double precision
#sgl      real
     &     cep(mxvl,nstr,*), e(*), nu(*), dj(*), c1(mxvl), c2(mxvl),
     &     c3(mxvl), c4(mxvl) ,w, zero, one, two
#dbl      data zero,one,two / 0.0d00, 1.0d00, 2.0d00 /
#sgl      data zero,one,two /0.0, 1.0, 2.0/
c
c
      do i= 1,span
         cep(i,1,4)= zero
         cep(i,1,5)= zero
         cep(i,1,6)= zero
         cep(i,2,4)= zero
         cep(i,2,5)= zero
         cep(i,2,6)= zero
         cep(i,3,4)= zero
         cep(i,3,5)= zero
         cep(i,3,6)= zero
         cep(i,4,1)= zero
         cep(i,4,2)= zero
         cep(i,4,3)= zero
         cep(i,4,5)= zero
         cep(i,4,6)= zero
         cep(i,5,1)= zero
         cep(i,5,2)= zero
         cep(i,5,3)= zero
         cep(i,5,4)= zero
         cep(i,5,6)= zero
         cep(i,6,1)= zero
         cep(i,6,2)= zero
         cep(i,6,3)= zero
         cep(i,6,4)= zero
         cep(i,6,5)= zero
c
         c1(i)= (e(i)/((one+nu(i))*(one-two*nu(i))))*dj(i)*w
         c2(i)= (one-nu(i))*c1(i)
         c3(i)= ((one-two*nu(i))/two)*c1(i)
         c4(i)= nu(i)*c1(i)
c
         cep(i,1,1)= c2(i)
         cep(i,2,2)= c2(i)
         cep(i,3,3)= c2(i)
         cep(i,4,4)= c3(i)
         cep(i,5,5)= c3(i)
         cep(i,6,6)= c3(i)
         cep(i,1,2)= c4(i)
         cep(i,1,3)= c4(i)
         cep(i,2,1)= c4(i)
         cep(i,3,1)= c4(i)
         cep(i,2,3)= c4(i)
         cep(i,3,2)= c4(i)
c
      end do
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                     Subroutine adjust_cnst                   *
c     *                                                              *
c     *                       written by: gvt                        *
c     *                   last modified : 08/25/98                   *
c     *                                                              *
c     ****************************************************************
c
c
c         Include other scalars with the linear material matrix required for
c         axisymmetric or planar elements.
c
c         For triangular elements, the correct area is given by
c         0.5*|J|.
c
c         For axisymmetric elements, include the 2*pi*radius scalar term
c         in the element stiffness summation.
c
c         Determine the correct scalar value and then multiply all the
c         terms in the material matrix [D] = cep_block.
c
c         Variables:
c
c         elem_type = integer flag for the element type, should get one of
c                     8=tri6, 10=axiquad8, 11=axitri6
c         nstr  = number of strains, should be 6
c         mxvl  = maximum number of elements in a block
c         span  = number of elements in the block
c         rad() = radius to the current Gauss point for each element in the
c                 current block.
c         cep() = material matrix, should already have been updated, extra
c                 scalar values in the Gauss integration
c                 loop are included here into the material matrix.
c
      subroutine adjust_cnst( elem_type, nstr, mxvl, span, rad, cep )
      implicit none
c
c                  parameters
c
      integer elem_type, nstr, mxvl, span
#dbl      double precision
#sgl      real
     &         rad(*),
     &         cep(mxvl,nstr,*)
c
c                   local variables
c
      integer i, row
c
#dbl      double precision
#sgl      real
     &         scalar
c
c                    the element type determines the scalar
c                    multiple to get the correct adjustment to |J|.
c
      call adjust_scalar_weights( elem_type, scalar )
      do row = 1, nstr
          do i = 1, span
              cep(i,row,1) = scalar*cep(i,row,1)
              cep(i,row,2) = scalar*cep(i,row,2)
              cep(i,row,3) = scalar*cep(i,row,3)
              cep(i,row,4) = scalar*cep(i,row,4)
              cep(i,row,5) = scalar*cep(i,row,5)
              cep(i,row,6) = scalar*cep(i,row,6)
          end do
      end do
c
      return
      end

