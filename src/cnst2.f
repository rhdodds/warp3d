c *********************************************************************
c *                                                                   *
c *   cnst2 -- tangent modulus matrix for nonlinear elastic model     *
c *            vectorized version                                     *
c *                                                                   *
c *********************************************************************
c
      subroutine cnst2( felem, gpn, e, nu, sigyld,
     &                  exp, strain_n1, history, cep, span, dj, w )
      implicit integer (a-z)
$add param_def
c
c                   parameter declarations
c   
#dbl      double precision
#sgl      real
     & strain_n1(mxvl,*), history(span,*), cep(mxvl,6,6), e(*),
     & nu(*), sigyld(*), exp(*), dj(*), w
c
c                   locally defined 
c
#dbl      double precision
#sgl      real
     & k1, k2(mxvl), epsyld(mxvl), epslim(mxvl), sigeff(mxvl),
     & epseff(mxvl), l2, c1, c2, c3, c4, theta, pi2,
     & root2, zero, one, two, three,
     & epsnc(mxvl), signc(mxvl), rnc(mxvl)
      logical debug, nonlinear_points, nonlinear_flags(mxvl)       
c
      data pi2, root2
     & / 1.570795, 1.4142135623730  /
      data zero, one, two, k1, three
     & / 0.0, 1.0, 2.0, 0.95, 3.0/  
c
c              pull out model properties
c
      nonlinear_points = .false.
      nonlin_point     = 0
      do i = 1, span
        epsyld(i) = sigyld(i) / e(i)
        k2(i)     = root2 * ( one-k1 ) / exp(i) + one
        epslim(i) = k1 * epsyld(i) * (two/three) * (one+nu(i))
        sigeff(i) = history(i,2)
        epseff(i) = history(i,3)
      end do
      do i = 1, span
        nonlinear_flags(i) = .false.
        if ( epseff(i) .gt. epslim(i) ) then
           nonlinear_points = .true.
           if ( nonlin_point .eq. 0 ) nonlin_point = i
           nonlinear_flags(i) = .true.
        end if
      end do
c
      debug = .false.
      if ( debug ) then
        call iodevn( innum, iout, dummy, 1 )
        write(iout,*) ' '
        i = nonlin_point
        write(iout,9000) felem+i-1, gpn, e(i), nu(i), sigyld(i), exp(i),
     &                   (strain_n1(i,j),j=1,6), (history(i,j),j=1,3)
      end if
c
      do  i = 1, span
       if ( .not. nonlinear_flags(i) )  then
c
c              strain point is currently linear. return linear
c              elastic, isotropic [d].
c
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
          c1 = (e(i)/((one+nu(i))*(one-two*nu(i))))*dj(i)*w
          c2 = (one-nu(i))*c1   
          c3 = ((one-two*nu(i))/two)*c1
          c4 = nu(i)*c1
          cep(i,1,1)= c2
          cep(i,2,2)= c2
          cep(i,3,3)= c2
          cep(i,4,4)= c3
          cep(i,5,5)= c3
          cep(i,6,6)= c3
          cep(i,1,2)= c4
          cep(i,1,3)= c4
          cep(i,2,1)= c4
          cep(i,3,1)= c4
          cep(i,2,3)= c4
          cep(i,3,2)= c4
       end if
      end do
c      
      if ( .not. nonlinear_points ) return
c
c             some points are nonlinear
c
      do i = 1, span
       if ( nonlinear_flags(i) ) then
        theta    = atan( k2(i)**(1.-exp(i))/exp(i) )
        l2       = -tan(pi2-theta)
        c1       = one / (l2 + one)
        c2       = k2(i)**exp(i)
        epsnc(i) = ( two*k1 + l2*c2 - k2(i) ) * c1
        signc(i) = ( two*l2*k1 - l2*c2 + k2(i) ) * c1
        rnc(i)   = root2 * ( k1 - l2*k1 + l2*c2 - k2(i) ) * c1
       end if
      end do
c
      call cnst2a( cep, e, nu, sigeff, epseff, strain_n1,
     &             sigyld, exp, signc, rnc, k2,
     &             span, dj, w, nonlinear_flags )
c
      if ( .not. debug ) return
        write(iout,9034)
        i = nonlin_point
        do j = 1, 6
          write(iout,9008) cep(i,j,1), cep(i,j,2), cep(i,j,3), 
     &                     cep(i,j,4), cep(i,j,5), cep(i,j,6)
        end do
      write(iout,*) ' '
      return
c
 9000 format('>> debug from cnst2. elem, gpn : ',i8,i2,
     & /,    '    e, nu, sigyld, exp : ',4f10.2,
     & /,    '    strains @ n+1 :',
     & /,10x,3e15.6,/,10x,3e15.6,
     & /,    '    history : ',f4.1,f10.3,f10.6 )
 9008 format(3x,6e14.4)
 9020 format(//,5x,'>> update [d] element, point: ',2i6)
 9034 format('    elasto-plastic d matrix: ' )
 9038 format('    point is elastic.  epseff, sigeff: ',e15.6,
     & f10.3 )
c
      end
c ********************************************************************
c *                                                                  *
c *   cnst2a -- tangent modulus matrix                               *
c *                                                                  *
c ********************************************************************
c
c
      subroutine cnst2a( cep, e, nu, sigeff, epseff, strain, 
     &                   sigyld, exp, signc, rnc, k2,
     &                   span, dj, w, nonlinear_flags )
      implicit integer (a-z)
$add param_def
c
#dbl      double precision
#sgl      real            
     & cep(mxvl,6,6), e(*), nu(*), sigeff(*), epseff(*),
     & strain(mxvl,*), sigyld(*), exp(*), signc(*),
     & rnc(*), k2(*), dj(*), w
      logical  nonlinear_flags(*)
c
c
#dbl      double precision
#sgl      real            
     & ev(mxvl), e1(mxvl), e2(mxvl), e3(mxvl),
     & e4(mxvl), e5(mxvl), e6(mxvl), sign(mxvl), radical,
     & stiff, const, g(mxvl), expm1, c1,
     & c2(mxvl), c3(mxvl), temp, denom, wf, 
     & twothd, fnine, third, half, one, two 
      data  twothd, fnine, third, half, one, two
     & / 0.6666666667, 0.444444444444444, 0.3333333333333, 0.5,
     &   1.0, 2.0 /
c
c             compute the tangent modulus matrix
c
      do i = 1, span
       if ( nonlinear_flags(i) ) then
         ev(i) = third * ( strain(i,1) + strain(i,2) + strain(i,3) )
         e1(i) = strain(i,1) - ev(i)
         e2(i) = strain(i,2) - ev(i)
         e3(i) = strain(i,3) - ev(i)
         e4(i) = strain(i,4) * half
         e5(i) = strain(i,5) * half
         e6(i) = strain(i,6) * half
         sign(i) = sigeff(i) / sigyld(i)
       end if
      end do
c
      do i = 1, span
      if ( nonlinear_flags(i) ) then
       if ( sign(i) .le. k2(i) ) then
        radical = sqrt( rnc(i)*rnc(i) - (sign(i)-signc(i))**2 )
        stiff   = (sign(i) - signc(i))/radical -
     &                twothd*(half-nu(i))
        const   = one/epseff(i) - e(i)/sigeff(i)/stiff
        g(i)    = -fnine * sigeff(i) * const / epseff(i)
     &                /epseff(i)
       else
        expm1 = exp(i) - one
        c1    = sign(i)**expm1
        denom = twothd * ( -half + nu(i) ) + exp(i) * c1
        temp  = -fnine * sigeff(i) * expm1 * c1 / epseff(i)
     &              / epseff(i) /epseff(i)
        g(i)  = temp / denom
       end if
       c2(i) = twothd * sigeff(i) / epseff(i)
       c3(i) = third * ( -c2(i) + e(i)/(one-two*nu(i)) )
      end if
      end do
c
c             3-D (row ordering x, y, z, xy, yz, xz)
c
      do i = 1, span
       if ( nonlinear_flags(i) ) then
          wf = dj(i) * w
          cep(i,1,1) = ( c3(i) + g(i)*e1(i)*e1(i) + c2(i) ) * wf
          cep(i,2,1) = ( c3(i) + g(i)*e2(i)*e1(i) ) * wf
          cep(i,3,1) = ( c3(i) + g(i)*e3(i)*e1(i) ) * wf
          cep(i,2,2) = ( c3(i) + g(i)*e2(i)*e2(i) + c2(i) ) * wf
          cep(i,3,2) = ( c3(i) + g(i)*e3(i)*e2(i) ) * wf
          cep(i,3,3) = ( c3(i) + g(i)*e3(i)*e3(i) + c2(i) ) * wf
          cep(i,1,2) = cep(i,2,1)
          cep(i,1,3) = cep(i,3,1)
          cep(i,2,3) = cep(i,3,2)
          cep(i,1,4) = ( g(i)*e1(i)*e4(i) ) * wf
          cep(i,1,5) = ( g(i)*e1(i)*e5(i) ) * wf
          cep(i,1,6) = ( g(i)*e1(i)*e6(i) ) * wf
          cep(i,2,4) = ( g(i)*e2(i)*e4(i) ) * wf
          cep(i,2,5) = ( g(i)*e2(i)*e5(i) ) * wf
          cep(i,2,6) = ( g(i)*e2(i)*e6(i) ) * wf
          cep(i,3,4) = ( g(i)*e3(i)*e4(i) ) * wf
          cep(i,3,5) = ( g(i)*e3(i)*e5(i) ) * wf
          cep(i,3,6) = ( g(i)*e3(i)*e6(i) ) * wf
          cep(i,4,1) = cep(i,1,4)
          cep(i,5,1) = cep(i,1,5)
          cep(i,6,1) = cep(i,1,6)
          cep(i,4,2) = cep(i,2,4)
          cep(i,5,2) = cep(i,2,5)
          cep(i,6,2) = cep(i,2,6)
          cep(i,4,3) = cep(i,3,4)
          cep(i,5,3) = cep(i,3,5)
          cep(i,6,3) = cep(i,3,6)
          cep(i,4,4) = ( g(i)*e4(i)*e4(i) + c2(i)*half ) * wf
          cep(i,4,5) = ( g(i)*e4(i)*e5(i) ) * wf
          cep(i,4,6) = ( g(i)*e4(i)*e6(i) ) * wf
          cep(i,5,5) = ( g(i)*e5(i)*e5(i) + c2(i)*half ) * wf
          cep(i,5,6) = ( g(i)*e5(i)*e6(i) ) * wf
          cep(i,6,6) = ( g(i)*e6(i)*e6(i) + c2(i)*half ) * wf
          cep(i,5,4) = cep(i,4,5)
          cep(i,6,4) = cep(i,4,6)
          cep(i,6,5) = cep(i,5,6)
       end if
      end do
c
      return
      end

