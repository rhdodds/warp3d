c *********************************************************************
c *                                                                   *
c *   cnst3 -- tangent modulus matrix for gurson model                *
c *            vectorized version                                     *
c *                                                                   *
c *********************************************************************
c
c
c
c      parameter definitions:
c      ======================
c
c                  mode/type
c
c   element          1,1         -- element number being processed
c
c   gpn              1,1         -- gauss point (or strain point) being
c                                   processed
c
c   first            1,3         -- a logical flag indicating if this
c                                   is the very first call to the routine
c                                   the linear elastic [d] is returned
c                                   whenever first = .true.
c
c   iter             1,1         -- equilibrium iteration for the load
c                                   step. iter 1 is for application
c                                   of the "real" load increment.
c
c   e                1,2         -- young's modulus
c 
c   nu               1,2         -- poisson's ratio
c
c   q1, q2,          1,2         -- q1, q2, q3 constants used in gurson 
c   q3                              model
c
c   nucleation       1,3         -- .true. if void nucleation is to be
c                                   included in the stress update
c
c   nuc_s_n          1,2         -- nucleation constants for gurson
c   nuc_e_n,                        model (s_n, e_n, f_n)
c   nuc_f_n  
c
c   stress_trial     1,2         -- 6x1 vector of the trial elastic
c                                   stress state for n+1. not needed
c                                   if first = .true. computed by
c                                   mm03 for use here. for first 
c                                   iteration of a step (application
c                                   of real load), this should be
c                                   final trial elastic stress state 
c                                   computed in resolution of previous
c                                   load step. see mm03 for ordering. 
c
c   history          1,2         -- 9x1 vector of history data for the
c                                   gauss point. passed in as history
c                                   state at 'n'. not needed if
c                                   first=.true.
c
c   history1         1,2         -- 9x1 vector of history data for the
c                                   gauss point. passed in as history
c                                   state for current estimate
c                                   of the solution at 'n+1'. not needed
c                                   if first=.true.
c
c   cep              2,2         -- 6x6 cepent [d]
c
c   iout             1,1         -- output device number for error messages.
c
c   span             1,1         -- number of elements in this block
c
c   dj               1,2         -- determinant of coordinate jacobian
c
c   w                1,2         -- gauss integration weight value at
c                                   point
c
c
c   routine names:
c   =============
c
c         all gurson routines are named with the prefix mm03
c
c   mode keys:
c   =========
c         1 (input)
c         2 (output)
c         3 (input by calling routine, updated by this routine and
c            returned)
c
c
c   type keys:
c   =========
c         1 (integer)
c         2 (floating point; real or double depending on computer)
c         3 (logical)
c
c
c   double-single precision:
c   =======================
c
c         we use a #dbl or #sgl at beginning of lines which
c         depend on single/double precision implementations.
c
c
c
c
c   large-strain issues:
c   ====================
c
c         this routine is unaware of finite-strain vs. small-strain
c         issues. we assume that "rotation neutral" stresses
c         are passed as arguments. this can be accomplished
c         in a variety of ways including:
c           a) a jaumann rate approach with the hughes-winget technique
c              to handle the rotation
c           b) using the "unrotated" reference configuration of
c              dienes, halquist, taylor & flanagan, dodds & healy...
c
      subroutine cnst3( element,  gpn, first, iter, e, nu, q1,
     &                  q2, q3,
     &                  nucleation, nuc_s_n, nuc_e_n, nuc_f_n,
     &                  stress_trial, history, history1, cep,
     &                  span, dj, w )
      implicit integer (a-z)
$add param_def
c
c                   parameter declarations
c 
#dbl      double precision
#sgl      real
     & stress_trial(mxvl,*), history(span,*),
     & history1(span,*), cep(mxvl,6,6), e(*), nu(*), q1(*), 
     & q2(*), q3(*), nuc_s_n(*), nuc_e_n(*), nuc_f_n(*), dj(*), w
      logical nucleation(*), first  
c
c                   locally defined array, variables
c
      logical debug, ldummy
#dbl      double precision
#sgl      real
     & dword(mxvl), mpi(mxvl), mpn(mxvl), c1(mxvl), c2(mxvl),
     & c3(mxvl), c4(mxvl), ebarp(mxvl), sbar(mxvl),
     & pn1(mxvl), hprime(mxvl), f(mxvl), dep(mxvl), deq(mxvl),
     & mqi(mxvl), mqn(mxvl), qn1(mxvl), ebarpn(mxvl),
     & fn(mxvl), shear_mod(mxvl), bulk_mod(mxvl),
     & twog(mxvl), threeg(mxvl),
     & debarp(mxvl), sbar2(mxvl), sbar3(mxvl), p(mxvl),
     & q(mxvl), sm(mxvl), sx(mxvl), sy(mxvl), sz(mxvl),
     & qe(mxvl), pe(mxvl), sl(mxvl), n(mxvl,6), a1(mxvl),
     & a2(mxvl), a3(mxvl), a4(mxvl), a5(mxvl), anuc(mxvl),
     & anuc_prime(mxvl) 
#dbl      double precision
#sgl      real
     & twothd, third, half, root32, zero, one, two, three, root_2_pi,
     & term1, term2, temp, b1, b2, d1, d2, d3, d4, h1, h2, h3, h4,
     & h5, h6, h7, h8, beta, ch, sh, pgp, pgq, pgsbar, pgf,
     & cap_a11, cap_a12, cap_b1, cap_c1, c5, c6, e10, e11, e12,
     & e13, m1, m2, c7, c8, c9, d10, d11, d12, d13,
     & cap_a21, cap_a22, cap_b2, cap_c2, denom, h10, h11, h12, h13,
     & con_1, con_2, con_3, wf
      logical   nonlinear_points
      integer   iword(mxvl*2), state(mxvl)
      equivalence (iword, dword )
      data      twothd, third / 0.6666666667, 0.333333333 /
      data      half / 0.5 /, root32 / 1.224744871 /
      data      zero, one, two, three / 0.0, 1.0, 2.0, 3.0 /
      data      root_2_pi / 2.50663 /
c
c
c        cep (output)       -- 6x6 update elastic-plastic [d].
c        pn1    (input)     -- (-) mean macrostress at n+1
c        qn1    (input)     -- equivalent macrostress at n+1
c        hprime (input)     -- current (matrix) plastic modulus (n+1)
c        sbar   (input)     -- current (matrix) equivalent stress (n+1)
c        ebarp  (input)     -- current plastic strain in matrix (n+1)
c        ebarpn (input)     -- plastic strain in matrix at start of
c                              step (n)
c        f      (input)     -- current void volume fraction (n+1)
c        fn     (input)     -- void volume fraction at start
c                              of step (n)
c        dep    (input)     -- increment of (macro) plastic volume
c                              strain over step
c        deq    (input)     -- increment of (macro) plastic deviatoric
c                              strain over step
c        shear_mod (input)  -- elastic shear modulus
c        bulk_mod (input)   -- bulk modulus
c        q1, q2, q3 (input) -- Gurson model constants
c        nucleation (input) -- logical true if nucleation is to
c                              be modeled
c        nuc_s_n
c        nuc_e_n (input)    -- constants for nucleation in
c        nuc_f_n               Gurson's model
c        sigx, y ...        -- trial elastic stresses at n+1
c
c        iter   (input)     -- current iteration number
c
c
c              set the debugging output level for routine. dump
c              key parameters if debugging.
c
      call iodevn( idummy, iout, ldummy, 1 )
      if( first ) then
        write(iout,9600)
        call die_abort
      end if
c
      debug = .false.
      if ( debug ) then
        write(iout,9000) element, gpn, iter, e(1), nu(1), dj(1), w,
     &                   q1(1), q2(1), q3(1), nucleation(1),
     &                   nuc_s_n(1), nuc_e_n(1), nuc_f_n(1), first,
     &                   (stress_trial(1,i),i=1,6) 
      end if
c
c              get material state variable for gauss point from
c              history vector. if first = .t., this is the first
c              time ever in here (step=iter=1).
c              two states are possible: = -1 implies strain point
c              is currently in a linear state. = 1 implies 
c              strain point is currently undergoing plastic loading.
c
      nonlinear_points = .false.
      nonlin_point    = 0
c
      do i = 1, span
       dword(i) = history1(i,6)
      end do
c
#dbl      inc_factor = 2
#sgl      inc_factor = 1
      j = 1
      do i = 1, span
       state(i) = iword(j)
       if ( first ) state(i) = -1
       if ( state(i) .eq. 1 ) then
          nonlinear_points = .true.
          if ( nonlin_point .eq. 0 ) nonlin_point = i
       end if
       j = j + inc_factor
      end do
c
      do  i = 1, span
       if ( state(i) .eq. -1 ) then
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
          c1(i)= (e(i)/((one+nu(i))*(one-two*nu(i))))*dj(i)*w
          c2(i)= (one-nu(i))*c1(i)   
          c3(i)= ((one-two*nu(i))/two)*c1(i)
          c4(i)= nu(i)*c1(i)
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
       end if
      end do
c      
      if ( .not. nonlinear_points ) return
c
c              points are elastic-plastic.  pull data from history.
c              get the cepent [dep] for the strain point.
c              for iteration 1 of a step we have a special situation.
c              the two history vectors are identical. the trial
c              elastic stress vector was computed
c              as part of the last stress update of the previous step.
c              it is ->not<- the current contact stress.
c              we set want a "continuum" tangent in this case. we
c              set  dep, deq = 0 (pn1, qn1 scale the trial stress to
c              be the contact stress).
c              a continnum tangent at the contact stress point on the
c              yield surface is computed.
c
      do i = 1, span
        if ( state(i) .eq. 1 ) then
          ebarp(i)     = history1(i,1)
          sbar(i)      = history1(i,2)
          pn1(i)       = history1(i,3)
          hprime(i)    = history1(i,4)
          f(i)         = history1(i,5)
          dep(i)       = history1(i,7)
          deq(i)       = history1(i,8)
          qn1(i)       = history1(i,9)
          ebarpn(i)    = history(i,1)
          fn(i)        = history(i,5)
          shear_mod(i) = half * e(i) /(one+nu(i))
          bulk_mod(i)  = e(i) * third / ( one - two*nu(i) )
          if ( iter .le. 1 ) then
            dep(i) = zero
            deq(i) = zero
          end if
        end if
      end do
c
c            1) compute some frequently used constants
c
c            2) pressure and equivalent stress for macrostresses
c               at n+1. passed down from history state at (n).
c
c
      if ( debug ) then
        i = nonlin_point
        write(iout,9109) i
        write(iout,9110) shear_mod(i), hprime(i), sbar(i),
     &    ebarp(i), ebarpn(i), f(i), fn(i), dep(i), deq(i),
     &    q1(i), q2(i), q3(i),  nuc_s_n(i), nuc_e_n(i),
     &    nuc_f_n(i), nucleation(i), stress_trial(i,1),
     &    stress_trial(i,2), stress_trial(i,3),
     &    stress_trial(i,4), stress_trial(i,5), stress_trial(i,6)
      end if
c 
      do i = 1, span
        if ( state(i) .eq. 1 ) then
          twog(i)   = two * shear_mod(i)
          threeg(i) = three * shear_mod(i)
          debarp(i) = ebarp(i) - ebarpn(i)
          sbar2(i)  = sbar(i) * sbar(i)
          sbar3(i)  = sbar(i)**3
          p(i)      = pn1(i)
          q(i)      = qn1(i)
        end if
      end do
c
c            3) process terms of trial elastic stress state at n+1.
c               get pressure, equivalent stress and yield surface
c               normal.
c
      do i = 1, span
        if ( state(i) .eq. 1 ) then
         sm(i)  = ( stress_trial(i,1) + stress_trial(i,2) +
     &              stress_trial(i,3) ) * third
         sx(i)  = stress_trial(i,1) - sm(i)
         sy(i)  = stress_trial(i,2) - sm(i)
         sz(i)  = stress_trial(i,3) - sm(i)
         qe(i)  = root32 * sqrt( sx(i)*sx(i) + sy(i)*sy(i) +
     &            sz(i)*sz(i) + two *
     &            ( stress_trial(i,4)*stress_trial(i,4) +
     &            stress_trial(i,5)*stress_trial(i,5) +
     &            stress_trial(i,6) * stress_trial(i,6) ) )
         pe(i)  = -sm(i)
         sl(i)  = three / two / qe(i)
         n(i,1) = sx(i) * sl(i)
         n(i,2) = sy(i) * sl(i)
         n(i,3) = sz(i) * sl(i)
         n(i,4) = stress_trial(i,4) * sl(i)
         n(i,5) = stress_trial(i,5) * sl(i)
         n(i,6) = stress_trial(i,6) * sl(i)
        end if
      end do
c
      if ( debug ) then
         i = nonlin_point
         write(iout,*) ' '
         write(iout,9120) q(i), qe(i), p(i), pe(i), n(i,1), n(i,2),
     &                    n(i,3), n(i,4), n(i,5), n(i,6)
      end if
c
c            4)  evaluate various constants defined on pages CT-dg-3
c                CT-dg-7 of the notes. the names here are those used
c                in the notes. equation numbers in the notes are
c                referred to throughout.
c
c
c                  4a)  a1 and a1 from Eq. (4)
c
      do i = 1, span
        if ( state(i) .eq. 1 ) then
          a1(i) = (deq(i) * q(i) - dep(i) * p(i)) /
     &            ( one - f(i) )**2 / sbar(i)
          a2(i) = (deq(i) * q(i) - dep(i) * p(i)) * hprime(i)  /
     &            ( f(i) - one ) / sbar2(i)
c
c                  4b)  a3, a4, a5 from Eq. (9)
c
          a3(i) = one / sbar(i) / ( f(i) - one )
          a4(i) = p(i)  +  dep(i) * bulk_mod(i)
          a5(i) = threeg(i) * deq(i)  -  q(i)
        end if
      end do
c
c                  4c)  b1 and b2 from Eq. (15). need to evaluate
c                       A(ebarp) for nucleation and its derivative
c                       wrt to ebarp.
c
      do i = 1, span
        if ( state(i) .eq. 1 ) then
         if ( nucleation(i) ) then 
           term1   = nuc_f_n(i) / nuc_s_n(i) / root_2_pi
           term2   = ( ( ebarp(i) - nuc_e_n(i) ) / nuc_s_n(i) )**2
           anuc(i) = term1 * exp(-half*term2)
           temp    = -nuc_f_n(i) * (ebarp(i)-nuc_e_n(i)) /
     &                nuc_s_n(i)**3 / root_2_pi
           anuc_prime(i) = temp * exp(-half*term2) 
         else
           anuc(i)       = zero
           anuc_prime(i) = zero
         end if
        end if
      end do
c
      do i = 1, span
        if ( state(i) .eq. 1 ) then
        b1 = ( one  -  f(i) ) / ( one  +  dep(i) )
        b2 = ( anuc(i)  +  debarp(i) * anuc_prime(i) ) /
     &       ( one + dep(i) )
c
c                  4d)  d1 -> d4 from Eqs. (20) -> (24)
c
        d1 = ( a2(i) - one ) / ( a1(i) * b2  +  a2(i) -  one )
        d2 = -b2 / ( a1(i) * b2  +  a2(i)  -  one )
        d3 = one / ( one  -  a2(i)  -  a1(i) * b2 )
        d4 = a1(i)  / ( one  -  a2(i)  -  a1(i) * b2 ) 
c
c                  4e)  h1 -> h8 from Eqs. (27) and (29)
c
        h1 = d1 * b1  +  d2 * a3(i) * a4(i)
        h2 = d2 * a3(i)  * a5(i)
        h3 = a3(i) * dep(i) * d2
        h4 = a3(i) * deq(i) * d2
        h5 = a3(i) * a4(i)  * d3  +  b1 * d4
        h6 = a3(i) * a5(i)  * d3
        h7 = a3(i) * dep(i) * d3
        h8 = a3(i) * deq(i) * d3
c
c                  4f)  A11, A12, B1, C1 in Eq. (32). first we
c                       need numerical values for derivatives of 
c                       Gurson yield function at end of step. use
c                       cap_ pre-fix to denote uc symbols in text.
c
        beta = -(one + half) * q2(i) 
        ch   = cosh( beta*p(i)/sbar(i) )
        sh   = sinh( beta*p(i)/sbar(i) )
c             
c                       4f-1)  partial g / partial p
c
        pgp = two * beta * f(i) * q1(i) * sh / sbar(i)
c
c                       4f-2)  partial g / partial q
c
        pgq = two * q(i) / sbar2(i)
c
c                       4f-3)  partial g / partial sbar
c
        pgsbar = ( two / sbar3(i) ) *
     &           ( -q(i)*q(i) - beta * f(i) * p(i) * q1(i) *
     &            sbar(i) * sh )
c
c                       4f-4)  partial g / partial f
c
        pgf = two * ( -f(i) * q3(i) + q1(i) * ch )

c
        cap_a11 = h1 * pgf  +  bulk_mod(i) * pgp  + 
     &            hprime(i) * h5 * pgsbar
        cap_a12 = h2* pgf -  threeg(i) * pgq  +
     &            hprime(i) * h6 * pgsbar
        cap_b1  = ( h3 * pgf +  pgp  +  hprime(i) *
     &            h7 * pgsbar ) *  bulk_mod(i)
        cap_c1  = ( h4 * pgf - pgq + hprime(i) * 
     &            h8 * pgsbar) * twog(i) 
c      if ( debug ) then
c        write(iout,*) ' '
c        write(iout,9200) beta, ch, sh, pgp, pgq, pgsbar, pgf,
c     &                   cap_a11, cap_a12, cap_b1, cap_c1,
c     &                   h1,h2,h3,h4,h5,h6,h7,h8,d1,d2,d3,d4,
c     &                   anuc(1), anuc_prime(1),
c     &                   b1,b2, a1(1),a2(1),a3(1),a4(1),a5(1)
c      end if
c
c
c            5)  evaluate various constants defined on pages CT-R1-3
c                CT-R1-5 the notes. the names here are those used
c                in the notes. equation numbers in the notes are
c                referred to throughout.
c
c
c                  5a)  c5, c6 from Eq. (17)
c
        c5 = two / sbar2(i)
        c6 = two * two * q(i) * hprime(i) / sbar3(i)
c
c                  5b)  e10->e13 from Eq. (18)
c
        e10 = c6 * h7
        e11 = c5  +  c6 * h8
        e12 = c6 * h5
        e13 = threeg(i) * c5  +  c6 * h6
c
c                  5c)  m1, m2 and c7 -> c11 from Eqs.(6) and
c                       (8) -> (10). use ch and abs(sh) from above.
c
        m1 = two * q1(i)
        m2 = - three * q2(i) * half
        sh = sinh( m2 * p(i) / sbar(i) )
        ch = cosh( m2 * p(i) / sbar(i) )
        c7 = f(i) * m1 * m2 * m2 * ch / sbar2(i)
        c8 = m1 * m2 * sh / sbar(i)
        c9 = -f(i) * hprime(i) * ( m1 * m2 * m2 * p(i) * ch / sbar3(i) +
     &                      m1 * m2 * sh / sbar2(i) )
c
c                  5d)  d10 -> d13 from Eq. (20)
c      
        d10 = -c7  -  c8 * h3  -  c9 * h7
        d11 = -c9 * h8  -  c8 * h4
        d12 = c8 * h1  +  c9 * h5  +  c7 * bulk_mod(i)
        d13 = c8 * h2  +  c9 * h6
c
c                  5e)  A21, A22, B2, C2 from Eq. (24), (25)
c      
        cap_a21 = deq(i) * d12  -  dep(i) * e12  +  pgq
        cap_a22 = deq(i) * d13  -  dep(i) * e13  +  pgp
        cap_b2  = -( deq(i) * d10  +  dep(i) * e10 ) * bulk_mod(i)
        cap_c2  = -( deq(i) * d11  + dep(i) * e11 ) * twog(i)
c      if ( debug ) then
c        write(iout,9300) cap_a21, cap_a22, cap_b2, cap_c2,
c     &                   d10,d11,d12,d13, m1,m2, c7,c8,c9,
c     &                   e10,e11,e12,e13, c5,c6
c      end if
c
c
c            6)  solve the pair of linear equations to compute
c                numerical values for d(Dep) and d(deq). compute
c                constants in the solution h10 -> h13 using
c                Eqs. (26-27). The compute values for the 
c                terms mpi, mpn, mqi, mqn defined in Eq. (28 -> 32)
c
        denom = cap_a11 * cap_a22  -  cap_a12 * cap_a21
        h10   = cap_a22 / denom
        h11   = cap_a12 / denom
        h12   = cap_a11 / denom
        h13   = cap_a21 / denom
c
        mpi(i) = h10 * cap_b1  -  h11 * cap_b2
        mpn(i) = h10 * cap_c1  -  h11 * cap_c2
        mqi(i) = h12 * cap_b2  -  h13 * cap_b1
        mqn(i) = h12 * cap_c2  -  h13 * cap_c1
       end if
c      if ( debug ) then
c        write(iout,9400) denom, h10,h11,h12,h13,
c     &                   mpi(1), mpn(1), mqi(1), mqn(1)
c      end if
       end do
c
c
c            7)  the four most difficult scalars needed to form [D]
c                are now available (mpi, mpn, mqi, mqn).
c                the [D] ceps of 5, 6x6 matrices added together
c                with an appropriate constant in front of each 6x6.
c                add in each 6x6 in turn, following the notes on
c                pages CT-EF-##. note that we have a choice of
c                the full non-symmetric [D] or a symmetry version.
c                cep was zeroed in nldfep driver.
c
c
c                 7a) [D] (1) Eq. (17) and pg. CT-EF-7
c                 7b) [D] (2) Eq. (17) and pg. CT-EF-7
c                 7c) [D] (4) Eq. (17) and pg. CT-EF-9. exercise the
c                     symmetry option by default
c
      do i = 1, span
       if ( state(i) .eq. 1 ) then
        con_1 = twog(i) * ( one - threeg(i) * deq(i) / qe(i) )
        con_2 = bulk_mod(i) * ( one - mpi(i) )
        con_3 = half * ( twog(i) * mqi(i) + bulk_mod(i) * mpn(i) )
        cep(i,1,1) =  con_1 * twothd + con_2 - two * con_3 * n(i,1)
        cep(i,2,1) = -con_1 * third  + con_2 - con_3 * (n(i,2)+n(i,1))
        cep(i,3,1) = -con_1 * third  + con_2 - con_3 * (n(i,3)+n(i,1))
        cep(i,1,2) = -con_1 * third  + con_2 - con_3 * (n(i,1)+n(i,2))
        cep(i,2,2) =  con_1 * twothd + con_2 - two * con_3 * n(i,2)
        cep(i,3,2) = -con_1 * third  + con_2 - con_3 * (n(i,2)+n(i,3))
        cep(i,1,3) = -con_1 * third  + con_2 - con_3 * (n(i,1)+n(i,3))
        cep(i,2,3) = -con_1 * third  + con_2 - con_3 * (n(i,2)+n(i,3))
        cep(i,3,3) =  con_1 * twothd + con_2 - two * con_3 * n(i,3)
        cep(i,4,4) =  con_1 * half
        cep(i,5,5) =  con_1 * half
        cep(i,6,6) =  con_1 * half
        cep(i,4,1) = -con_3 * n(i,4)
        cep(i,5,1) = -con_3 * n(i,5)
        cep(i,6,1) = -con_3 * n(i,6)
        cep(i,4,2) = -con_3 * n(i,4)
        cep(i,5,2) = -con_3 * n(i,5)
        cep(i,6,2) = -con_3 * n(i,6)
        cep(i,4,3) = -con_3 * n(i,4)
        cep(i,5,3) = -con_3 * n(i,5)
        cep(i,6,3) = -con_3 * n(i,6)
        cep(i,1,4) = -con_3 * n(i,4)
        cep(i,1,5) = -con_3 * n(i,5)
        cep(i,1,6) = -con_3 * n(i,6)
        cep(i,2,4) = -con_3 * n(i,4)
        cep(i,2,5) = -con_3 * n(i,5)
        cep(i,2,6) = -con_3 * n(i,6)
        cep(i,3,4) = -con_3 * n(i,4)
        cep(i,3,5) = -con_3 * n(i,5)
        cep(i,3,6) = -con_3 * n(i,6)
        cep(i,4,5) =  zero
        cep(i,4,6) =  zero
        cep(i,5,4) =  zero
        cep(i,5,6) =  zero
        cep(i,6,4) =  zero
        cep(i,6,5) =  zero
       end if
      end do
c
c                 7d) [D] (3) Eq. (17) and pg. CT-EF-8. multiply in
c                     the weight factor and determinant of coord.
c                     jacobian.
c
      do i = 1, span
       if ( state(i) .eq. 1 ) then 
        wf = dj(i) * w
        con_1 = twog(i) * ( twog(i) * deq(i) / qe(i)  -  mqn(i) )
        cep(i,1,1) = (cep(i,1,1) + con_1 * n(i,1) * n(i,1)) * wf
        cep(i,2,1) = (cep(i,2,1) + con_1 * n(i,2) * n(i,1)) * wf
        cep(i,3,1) = (cep(i,3,1) + con_1 * n(i,3) * n(i,1)) * wf
        cep(i,4,1) = (cep(i,4,1) + con_1 * n(i,4) * n(i,1)) * wf
        cep(i,5,1) = (cep(i,5,1) + con_1 * n(i,5) * n(i,1)) * wf
        cep(i,6,1) = (cep(i,6,1) + con_1 * n(i,6) * n(i,1)) * wf
        cep(i,1,2) = (cep(i,1,2) + con_1 * n(i,1) * n(i,2)) * wf
        cep(i,2,2) = (cep(i,2,2) + con_1 * n(i,2) * n(i,2)) * wf
        cep(i,3,2) = (cep(i,3,2) + con_1 * n(i,3) * n(i,2)) * wf
        cep(i,4,2) = (cep(i,4,2) + con_1 * n(i,4) * n(i,2)) * wf
        cep(i,5,2) = (cep(i,5,2) + con_1 * n(i,5) * n(i,2)) * wf
        cep(i,6,2) = (cep(i,6,2) + con_1 * n(i,6) * n(i,2)) * wf
        cep(i,1,3) = (cep(i,1,3) + con_1 * n(i,1) * n(i,3)) * wf
        cep(i,2,3) = (cep(i,2,3) + con_1 * n(i,2) * n(i,3)) * wf
        cep(i,3,3) = (cep(i,3,3) + con_1 * n(i,3) * n(i,3)) * wf
        cep(i,4,3) = (cep(i,4,3) + con_1 * n(i,4) * n(i,3)) * wf
        cep(i,5,3) = (cep(i,5,3) + con_1 * n(i,5) * n(i,3)) * wf
        cep(i,6,3) = (cep(i,6,3) + con_1 * n(i,6) * n(i,3)) * wf
        cep(i,1,4) = (cep(i,1,4) + con_1 * n(i,1) * n(i,4)) * wf
        cep(i,2,4) = (cep(i,2,4) + con_1 * n(i,2) * n(i,4)) * wf
        cep(i,3,4) = (cep(i,3,4) + con_1 * n(i,3) * n(i,4)) * wf
        cep(i,4,4) = (cep(i,4,4) + con_1 * n(i,4) * n(i,4)) * wf
        cep(i,5,4) = (cep(i,5,4) + con_1 * n(i,5) * n(i,4)) * wf
        cep(i,6,4) = (cep(i,6,4) + con_1 * n(i,6) * n(i,4)) * wf
        cep(i,1,5) = (cep(i,1,5) + con_1 * n(i,1) * n(i,5)) * wf
        cep(i,2,5) = (cep(i,2,5) + con_1 * n(i,2) * n(i,5)) * wf
        cep(i,3,5) = (cep(i,3,5) + con_1 * n(i,3) * n(i,5)) * wf
        cep(i,4,5) = (cep(i,4,5) + con_1 * n(i,4) * n(i,5)) * wf
        cep(i,5,5) = (cep(i,5,5) + con_1 * n(i,5) * n(i,5)) * wf
        cep(i,6,5) = (cep(i,6,5) + con_1 * n(i,6) * n(i,5)) * wf
        cep(i,1,6) = (cep(i,1,6) + con_1 * n(i,1) * n(i,6)) * wf
        cep(i,2,6) = (cep(i,2,6) + con_1 * n(i,2) * n(i,6)) * wf
        cep(i,3,6) = (cep(i,3,6) + con_1 * n(i,3) * n(i,6)) * wf
        cep(i,4,6) = (cep(i,4,6) + con_1 * n(i,4) * n(i,6)) * wf
        cep(i,5,6) = (cep(i,5,6) + con_1 * n(i,5) * n(i,6)) * wf
        cep(i,6,6) = (cep(i,6,6) + con_1 * n(i,6) * n(i,6)) * wf
       end if
      end do
c
c               all done building 6x6 cep tangent.
c
      if ( debug ) then
        write(iout,9500)
        i = nonlin_point
        do j = 1, 6
          write(iout,9510) cep(i,j,1),cep(i,j,2),cep(i,j,3),
     &                     cep(i,j,4),cep(i,j,5),cep(i,j,6)
        end do
        write(iout,9520)
      end if
c
      return
c
 9000 format(/,'>> debug from cnst3 elem, gpn, iter : ',i8,2i2,
     & /,    '    e, nu, dj(1), w : ',4f15.7,
     & /,    '    q1, q2, q3, nucleation: ',3f10.3, l10,
     & /,    '    nuc_s_n, nuc_e_n, nuc_f_n, first  : ',3f10.3,l5,
     & /,    '    trial elastic stresses @ n+1 :',
     & /,10x,3e15.6,/,10x,3e15.6 )
 9030 format(//,5x,'>> point is elastic, [cep] = [delas]')
 9500 format('  >>>> [cep]: ', /)
 9510 format(2x,6e14.6)
 9520 format(//)
 9100 format('  >> update [cep].  current trial stress data: ',
     & /,    '      sigx    :',e14.6,
     & /,    '      sigy    :',e14.6,
     & /,    '      sigz    :',e14.6,
     & /,    '      sigxy   :',e14.6,
     & /,    '      sigyz   :',e14.6,
     & /,    '      sigxz   :',e14.6,
     & /,    '      sm      :',e14.6,
     & /,    '      sx      :',e14.6,
     & /,    '      sy      :',e14.6,
     & /,    '      sz      :',e14.6,
     & /,    '      sl      :',e14.6,
     & /,    '      nx      :',e14.6,
     & /,    '      ny      :',e14.6,
     & /,    '      nz      :',e14.6,
     & /,    '      nxy     :',e14.6,
     & /,    '      nyz     :',e14.6,
     & /,    '      nxz     :',e14.6,
     & /,    '      q       :',e14.6,
     & /,    '      p       :',e14.6, // )
 9109 format('  >> debugging for element in block: ',i4)
 9110 format('  >> arguments   passed: ',
     & /,    '      shear_mod    :',e14.6,
     & /,    '      hprime       :',e14.6,
     & /,    '      sbar         :',e14.6,
     & /,    '      ebarp        :',e14.6,
     & /,    '      ebarpn       :',e14.6,
     & /,    '      f            :',e14.6,
     & /,    '      fn           :',e14.6,
     & /,    '      dep          :',e14.6,
     & /,    '      deq          :',e14.6,
     & /,    '      q1           :',e14.6,
     & /,    '      q2           :',e14.6,
     & /,    '      q3           :',e14.6,
     & /,    '      nuc_s_n      :',e14.6,
     & /,    '      nuc_e_n      :',e14.6,
     & /,    '      nuc_f_n      :',e14.6,
     & /,    '      nucleation   :',l2,
     & /,    '      stress_el(1) :',e14.6,
     & /,    '      stress_el(2) :',e14.6,
     & /,    '      stress_el(3) :',e14.6,
     & /,    '      stress_el(4) :',e14.6,
     & /,    '      stress_el(5) :',e14.6,
     & /,    '      stress_el(6) :',e14.6,// )
 9120 format('  >> q, qe, p, pe :',4e14.6,
     & /,    '      nx      :',e14.6,
     & /,    '      ny      :',e14.6,
     & /,    '      nz      :',e14.6,
     & /,    '      nxy     :',e14.6,
     & /,    '      nyz     :',e14.6,
     & /,    '      nxz     :',e14.6, // )
 9200 format('  >> terms in section (4): ',
     & /,    '      beta         :',e14.6,
     & /,    '      ch           :',e14.6,
     & /,    '      sh           :',e14.6,
     & /,    '      pgp          :',e14.6,
     & /,    '      pgq          :',e14.6,
     & /,    '      pgsbar       :',e14.6,
     & /,    '      pgf          :',e14.6,
     & /,    '      cap_a11      :',e14.6,
     & /,    '      cap_a12      :',e14.6,
     & /,    '      cap_b1       :',e14.6,
     & /,    '      cap_c1       :',e14.6,
     & /,    '      h1           :',e14.6,
     & /,    '      h2           :',e14.6,
     & /,    '      h3           :',e14.6,
     & /,    '      h4           :',e14.6,
     & /,    '      h5           :',e14.6,
     & /,    '      h6           :',e14.6,
     & /,    '      h7           :',e14.6,
     & /,    '      h8           :',e14.6,
     & /,    '      d1           :',e14.6,
     & /,    '      d2           :',e14.6,
     & /,    '      d3           :',e14.6,
     & /,    '      d4           :',e14.6,
     & /,    '      anuc         :',e14.6,
     & /,    '      anuc_prime   :',e14.6,
     & /,    '      b1           :',e14.6,
     & /,    '      b2           :',e14.6,
     & /,    '      a1           :',e14.6,
     & /,    '      a2           :',e14.6,
     & /,    '      a3           :',e14.6,
     & /,    '      a4           :',e14.6,
     & /,    '      a5           :',e14.6, // )
 9300 format('  >> terms in section (5): ',
     & /,    '      cap_a21      :',e14.6,
     & /,    '      cap_a22      :',e14.6,
     & /,    '      cap_b2       :',e14.6,
     & /,    '      cap_c2       :',e14.6,
     & /,    '      d10          :',e14.6,
     & /,    '      d11          :',e14.6,
     & /,    '      d12          :',e14.6,
     & /,    '      d13          :',e14.6,
     & /,    '      m1           :',e14.6,
     & /,    '      m2           :',e14.6,
     & /,    '      c7           :',e14.6,
     & /,    '      c8           :',e14.6,
     & /,    '      c9           :',e14.6,
     & /,    '      e10          :',e14.6,
     & /,    '      e11          :',e14.6,
     & /,    '      e12          :',e14.6,
     & /,    '      e13          :',e14.6,
     & /,    '      c5           :',e14.6,
     & /,    '      c6           :',e14.6, //)
 9400 format('  >> terms in section (6): ',
     & /,    '      denom        :',e14.6,
     & /,    '      h10          :',e14.6,
     & /,    '      h11          :',e14.6,
     & /,    '      h12          :',e14.6,
     & /,    '      h13          :',e14.6,
     & /,    '      mpi          :',e14.6,
     & /,    '      mpn          :',e14.6,
     & /,    '      mqi          :',e14.6,
     & /,    '      mqn          :',e14.6, //)
 9600 format('>>>> FATAL ERROR: cnst3 has inconsistent state. ',
     & /,    '                  please contact WARP3D developers.',
     & /,    '                  Job terminated....' )
c
        end


