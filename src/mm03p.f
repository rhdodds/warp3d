c ************************************************************************
c *                                                                      *
c *    routine  mm03p -- vectorized pre-processor for block of mises or  *
c *                      gurson elements                                 *
c *                                                                      *
c ************************************************************************
c
c
      subroutine mm03p(
     &  step, iter, span, gpn, deps, e, nu, sigma_o, f0, q1, q2, q3,
     &  n_power, h_fixed, null_point, nucleation, stress_n, stress_n1,
     &  stress_n1_elas, p_trial, q_trial, history, history1,
     &  yld_func, nonlinear_flag, process_block, iout, segmental,
     &  curve_type, felem, ym_n, nu_at_n, hist_size )
      use segmental_curves, only : now_blk_relem

      implicit integer (a-z)
$add param_def
c
c               parameter declarations
c
#dbl      double precision
#sgl      real
     &  deps(mxvl,6), e(*), nu(*), sigma_o(*), f0(*), q1(*), q2(*),
     &  q3(*), stress_n(nstrs,*), stress_n1(nstrs,*),
     &  stress_n1_elas(mxvl,6,*), ym_n(*), nu_at_n(*),
     &  p_trial(*), q_trial(*), history(span,hist_size,*),
     &  history1(span,hist_size,*), yld_func(*), n_power(*), h_fixed(*)
c
      logical nucleation(*), nonlinear_flag(*), null_point(*), 
     &        process_block, segmental
c
c               locally allocated arrays for vectorization
c
#dbl      double precision 
#sgl      real
     &  a(mxvl), b(mxvl), c(mxvl), d(mxvl), shear_mod(mxvl),
     &  term1(mxvl), term2(mxvl), term3(mxvl),
     &  root22, third, zero, one, two, six, onehalf, toler,
     &  dword,  half, chk_yld_eps, h, mm03sc, plastic_strain,
     &  flow_stress(mxvl), hnow, ddumy, e_n, nu_n, g_n,
     &  een1, een2, een3, een4, een5, een6, e1, e2, e3, e4, e5,
     &  e6
      integer iword(2)
      logical    debug, gurson(mxvl), previously_linear, now_nonlinear
      equivalence ( dword, iword )
c
      data     root22 / 0.70710678 /, third / 0.3333333333333 /
      data     zero, one, two, six, onehalf / 0.0, 1.0, 2.0, 6.0, 1.5 /
      data     toler / 1.0e-10 /, half / 0.5 /
      data     chk_yld_eps / 0.002 /
c
c               initialize history on step 1
c
      debug = .false.
c
      if ( step .eq. 1 ) then
        iword(1)    = -1
        iword(2)    = 0
        if( .not. segmental ) then
           do i = 1, span
              if ( n_power(i) .gt. zero ) then
                 h = e(i) * (0.9999*e(i)) / ( e(i)-0.9999*e(i))
              else
                 h  = h_fixed(i)
              end if
              history(i,4,gpn)  = h
           end do
        endif
        do i = 1, span 
          history(i,1,gpn)  = zero
          history(i,2,gpn)  = sigma_o(i)
          history(i,3,gpn)  = zero
          history(i,5,gpn)  = f0(i)
          history(i,6,gpn)  = dword
          history(i,7,gpn)  = zero
          history(i,8,gpn)  = zero
          history(i,9,gpn)  = zero
          history(i,10,gpn) = zero
        end do
      end if
c
c               check for null gauss points
c
      if( debug ) write (iout,*) ' null point calcs'
      do i = 1, span
        null_point(i) = abs( deps(i,1) ) + abs( deps(i,2) ) +
     &                  abs( deps(i,3) ) + abs( deps(i,4) ) + 
     &                  abs( deps(i,5) ) + abs( deps(i,6) )
     &                  .le. toler * chk_yld_eps
      end do
c
      do i = 1, span
       if ( .not. null_point(i) ) go to 100
      end do
c
c               all points are null (no strain increment)
c
      state    = -1
      iword(1) = state
      do i = 1, span
        stress_n1(1,i)          = stress_n(1,i)
        stress_n1(2,i)          = stress_n(2,i)
        stress_n1(3,i)          = stress_n(3,i)
        stress_n1(4,i)          = stress_n(4,i)
        stress_n1(5,i)          = stress_n(5,i)
        stress_n1(6,i)          = stress_n(6,i)
        stress_n1(7,i)          = stress_n(7,i)
        stress_n1(8,i)          = stress_n(8,i)
        stress_n1(9,i)          = stress_n(9,i)
        stress_n1_elas(i,1,gpn) = zero   
        stress_n1_elas(i,2,gpn) = zero   
        stress_n1_elas(i,3,gpn) = zero   
        stress_n1_elas(i,4,gpn) = zero   
        stress_n1_elas(i,5,gpn) = zero   
        stress_n1_elas(i,6,gpn) = zero
c
        history1(i,1,gpn)  = history(i,1,gpn)
        history1(i,2,gpn)  = history(i,2,gpn)
        history1(i,3,gpn)  = history(i,3,gpn)
        history1(i,4,gpn)  = history(i,4,gpn)
        history1(i,5,gpn)  = history(i,5,gpn)
        history1(i,6,gpn)  = dword
        history1(i,7,gpn)  = zero
        history1(i,8,gpn)  = zero
        history1(i,9,gpn)  = zero
        history1(i,10,gpn) = zero
      end do
      process_block              = .false.  
      if ( debug ) write(iout,9000)
      return
c
c               we must process this block of elements.
c
 100  continue
      if ( debug ) write(iout,9010) step, iter, span, gpn
c
c
c               1.  trial elastic stress at n+1. set actual
c                   stress at n+1 as trial stress. 
c
      do i = 1, span
c
c                      elastic strain at n from stresses at n
c                      and elastic constants at n. a zero modulus
c                      at n means element must have been killed
c
        e_n  = ym_n(i)
        if ( e_n .gt. zero ) then
           nu_n = nu_at_n(i)
           g_n  = e_n/two/(one+nu_n)
           een1 = (stress_n(1,i)-nu_n*(stress_n(2,i)+stress_n(3,i)))/e_n
           een2 = (stress_n(2,i)-nu_n*(stress_n(1,i)+stress_n(3,i)))/e_n
           een3 = (stress_n(3,i)-nu_n*(stress_n(1,i)+stress_n(2,i)))/e_n
           een4 = stress_n(4,i) / g_n 
           een5 = stress_n(5,i) / g_n 
           een6 = stress_n(6,i) / g_n 
        else
           een1 = zero
           een2 = zero
           een3 = zero
           een4 = zero
           een5 = zero
           een6 = zero
        end if  
c
c                      strain for trial stress computation at n+1
c                      is elastic strain at n + increment of strain
c                      over step.
c        
        e1 = een1 + deps(i,1)
        e2 = een2 + deps(i,2)
        e3 = een3 + deps(i,3)
        e4 = een4 + deps(i,4)
        e5 = een5 + deps(i,5) 
        e6 = een6 + deps(i,6)
c
c                      trial elastic stress using strain at n+1 and
c                      elastic constants at n+1
c
        c(i) = e(i) / ( ( one + nu(i) ) * ( one - two * nu(i) ) )
        a(i) = c(i) * ( one - nu(i) )
        b(i) = c(i) * nu(i)
        shear_mod(i) = e(i) / ( two*(one+nu(i)))
c
        stress_n1(1,i) = a(i)*e1 + b(i)*(e2+e3)
        stress_n1(2,i) = a(i)*e2 + b(i)*(e1+e3)
        stress_n1(3,i) = a(i)*e3 + b(i)*(e1+e2)
        stress_n1(4,i) = shear_mod(i) * e4
        stress_n1(5,i) = shear_mod(i) * e5
        stress_n1(6,i) = shear_mod(i) * e6
c
        stress_n1_elas(i,1,gpn) = stress_n1(1,i)
        stress_n1_elas(i,2,gpn) = stress_n1(2,i)
        stress_n1_elas(i,3,gpn) = stress_n1(3,i)
        stress_n1_elas(i,4,gpn) = stress_n1(4,i)
        stress_n1_elas(i,5,gpn) = stress_n1(5,i)
        stress_n1_elas(i,6,gpn) = stress_n1(6,i)
c
      end do
c
      if ( debug ) then
        write(iout,9020)
        do i = 1, span
          write(iout,9030) i, (deps(i,j), 
     &                     stress_n(j,i), stress_n1(j,i), j = 1, 6 )
        end do
      end if
c
c
c              2.  hydrostatic and mises components of trial
c                  elastic stress at n+1
c
      do i = 1, span
        p_trial(i) = - ( stress_n1(1,i) + stress_n1(2,i) +
     &                   stress_n1(3,i) ) * third
        a(i)       = stress_n1(1,i) - stress_n1(3,i)
        b(i)       = stress_n1(1,i) - stress_n1(2,i)
        c(i)       = stress_n1(2,i) - stress_n1(3,i)
        d(i)       = stress_n1(4,i)*stress_n1(4,i) +
     &               stress_n1(5,i)*stress_n1(5,i) +
     &               stress_n1(6,i)*stress_n1(6,i)
        q_trial(i) = root22 * sqrt( a(i)*a(i) + b(i)*b(i) + 
     &                              c(i)*c(i) + six*d(i) )
       end do
c
      if ( debug ) then
        write(iout,9040)
        do i = 1, span
         write(iout,9042) i, p_trial(i), q_trial(i) 
        end do
      end if
c
c              3.  set the current flow stress for material to check
c                  loading/unloading status. this is needed only for
c                  temperature dependent segmental curves for which
c                  the flow properties may have changed from
c                  previous load step. using the current flow
c                  curve at this temperature and the total plastic
c                  strain in the history, get the new uniaxial
c                  stress. save in history as well so update
c                  routine knows about change in curve.
c
c
      do i = 1, span
        flow_stress(i) =  history(i,2,gpn)
      end do
c
      if ( segmental .and. curve_type .eq. 1 ) then
          do i = 1, span
            plastic_strain = history(i,1,gpn)
            now_blk_relem = i
            flow_stress(i) = mm03sc( plastic_strain, hnow, ddumy,
     &                               ddumy, 1 )
            history(i,2,gpn)  = flow_stress(i)
            history1(i,2,gpn) = flow_stress(i)
          end do
      end if 
c
c
c              4.  evaluate either mises yield function or the
c                  gurson- tvergaard yield function for the trial
c                  elastic stress state. current matrix stress is
c                  history(2), current void fraction is history(5).
c
      do i = 1, span
        gurson(i) = f0(i) .gt. zero .or. nucleation(i)
        if ( gurson(i) ) then
          term1(i)       = q_trial(i) / flow_stress(i)
          term2(i)       = two * q1(i) * history(i,5,gpn) *
     &                     cosh( onehalf*q2(i)*p_trial(i)/
     &                     flow_stress(i) )
          term3(i)       = one + q3(i) * history(i,5,gpn) *
     &                      history(i,5,gpn)
          yld_func(i) = term1(i)*term1(i) + term2(i) - term3(i)
        else
          yld_func(i) = q_trial(i) - flow_stress(i)
       end if
      end do
c
      if ( debug ) then
        write(iout,9050)
        do i = 1, span
         write(iout,9052) i, yld_func(i)
        end do
      end if
c
c              5.  set a flag to indicate if gauss point for
c                  each element in block is still linear. for iter=0,
c                  we are doing a stress estimate for imposed/
c                  extrapolated displacements. linear points must
c                  remain linear during this update.
c
      process_block = .false.
      do i = 1, span
        previously_linear = history(i,1,gpn) .eq. zero
        now_nonlinear     = (history(i,1,gpn) .gt. zero) .or.
     &                      (yld_func(i) .gt. zero)
        nonlinear_flag(i) = now_nonlinear
        process_block = process_block .or. nonlinear_flag(i)
      end do
c
      if ( debug ) then
         write(iout,9060)
         do i = 1, span
          write(iout,9062) i, nonlinear_flag(i)
         end do
         write(iout,9064) process_block
      end if
c
c              6.  for elements that are still linear, compute 
c                  updated strain energy density.
c
      do i= 1, span
        if ( .not. nonlinear_flag(i) ) then
          stress_n1(7,i) = stress_n(7,i) +
     &      half * ( 
     &      deps(i,1) * (stress_n1(1,i) + stress_n(1,i))
     &    + deps(i,2) * (stress_n1(2,i) + stress_n(2,i))
     &    + deps(i,3) * (stress_n1(3,i) + stress_n(3,i))
     &    + deps(i,4) * (stress_n1(4,i) + stress_n(4,i))
     &    + deps(i,5) * (stress_n1(5,i) + stress_n(5,i))
     &    + deps(i,6) * (stress_n1(6,i) + stress_n(6,i)) )
          stress_n1(8,i) = stress_n(8,i)
          stress_n1(9,i) = stress_n(9,i)
        end if
      end do    
c
      if ( debug ) then
         write(iout,9070)
         do i = 1, span
          if (.not.  nonlinear_flag(i) ) write(iout,9072) i,
     &            stress_n1(7,i)
         end do
      end if
c
c              7.  set element histories at n+1 to those at n
c                  those for nonlinear elements will be modified.
c
      do i= 1, span
        history1(i,1,gpn)  = history(i,1,gpn)
        history1(i,2,gpn)  = history(i,2,gpn)
        history1(i,3,gpn)  = history(i,3,gpn)
        history1(i,4,gpn)  = history(i,4,gpn)
        history1(i,5,gpn)  = history(i,5,gpn)
        history1(i,6,gpn)  = history(i,6,gpn)
        history1(i,7,gpn)  = history(i,7,gpn)
        history1(i,8,gpn)  = history(i,8,gpn)
        history1(i,9,gpn)  = history(i,9,gpn)
        history1(i,10,gpn) = history(i,10,gpn)
      end do 

c
      return
c
 9000 format(/,'>> all points in element block are null')
 9010 format(/,'>> model 3 pre-processor for step, iter, span, gpn: ',
     &  4i5)
 9020 format(/,'>> old and trial elastic stresses...',
     &       /,'      deps     stress_n  stress_n1')
 9030 format(/,i3, f11.8, 2f10.3,
     &         5(/,3x, f11.8, 2f10.3) )
 9040 format(/,'>>    p_trial   q_trial')
 9042 format(i3,2f10.3)
 9050 format(/,'>>    yield_function')
 9052 format(i3,f10.4)
 9060 format(/,'>>    nonlinear_flag')
 9062 format(i3,l4)
 9064 format(/,'>> process block flag: ',l1)
 9070 format(/,'>> updated energy densities for linear elements:')
 9072 format(i4,f10.6)
 9100 format('>> trial elastic state computation: ',
     &  /,   '      e, nu: ',f10.2,f10.3 )
 9110 format('      deps     dsigel       sigel ',
     & 6(/,1x,3f11.4) )
c
      end
c *******************************************************************
c *                                                                 *
c *      material model # 03 routine -- mm03f                       *
c *                                                                 *
c *******************************************************************
c
c
      function mm03f( stress, sbar, f, q1, q2, q3 )
c
c
c              compute the value of the gurson yield
c              function for the state of stress tau.
c
#dbl      implicit double precision (a-h,o-z)
#dbl      double precision
#sgl      real
     &   mm03f, mm03g
      dimension stress(*)
      data     root22 / 0.70710678 /, third / 0.3333333333333 /
      data     six / 6.0 /
c
      p  = - ( stress(1) + stress(2) + stress(3) ) * third
      a  = stress(1) - stress(3)
      b  = stress(1) - stress(2)
      c  = stress(2) - stress(3)
      d  = stress(4)*stress(4) + stress(5)*stress(5) +
     &     stress(6)*stress(6)
      q  = root22 * sqrt( a*a + b*b + c*c + six*d )
c
      mm03f = mm03g( q, sbar, f, p, q1, q2, q3 )
c
      return
      end
c *******************************************************************
c *                                                                 *
c *      material model # 03 routine -- mm03g                       *
c *                                                                 *
c *******************************************************************
c
c
      function mm03g( q, sbar, f, p, q1, q2, q3 )
c
c
c              compute the value of the Gurson yield
c              given the scaler parameters
c
#dbl      implicit double precision (a-h,o-z)
#dbl      double precision
#sgl      real
     &  mm03g
      data one, two, onehalf / 1.0, 2.0, 1.5 /
c
      term1  = (q/sbar)*(q/sbar)
      term2  = two * q1 * f * cosh(onehalf*q2*p/sbar)
      term3  = one + q3 * f * f
      mm03g  = term1 + term2 - term3
c
      return
      end
c ************************************************************************
c *                                                                      *
c *    routine  mm03us -- update state variables f, ebarp, sbar given    *
c *                       estimates for dep, deq. this is the fast       *
c *                       version called first. if this version fails    *
c *                       to converge, the slow one is called that uses  *
c *                       a newton-bisection method after bracketing the *
c *                       root first.                                    *
c *                                                                      *
c ************************************************************************
c
c
      subroutine mm03us( p_new, dep, q_new, deq, ebarp, f, sbar,
     &                   consth, sigma_o, eps_o, h_constant, 
     &                   e, n_power, m_power, dt, eps_ref, rate_depen,
     &                   mpoweri, nucleation, nuc_s_n, nuc_e_n,
     &                   nuc_f_n, ebarp_new, sbar_new, f_new, h_new,
     &                   iout, debug, segmental,
     &                   power_law, converge, first_time )
#dbl      implicit double precision (a-h,o-z)
#dbl      double precision
#sgl      real
     &  n_power, m_power, mpoweri, nuc_s_n, nuc_e_n, nuc_f_n
      logical debug, consth, nucleation, rate_depen,
     &         segmental, power_law, converge, first_time
      data    max_iterations / 10 /
      data    zero, one, tol / 0.0, 1.0, 0.000001 /
c
c           build iteration loop to find ebarp_new. in the process
c           we find f_new and sbar_new. if nucleation is neglected the
c           expressions are simpler. use a newton's method with iterations
c           on ebarp_new. an alternative estimate for ebarp_new is
c           just  ebarp_new = ebarp 
c
      t_numer        = -p_new * dep + q_new * deq
      ebarp_new      = ebarp + t_numer / ( one - f ) / sbar
      strain_toler   = tol * max( eps_o, ebarp_new )
      converge       = .true.
c
c
      do loop = 1, max_iterations
      if ( debug ) write(iout,9200) loop
c
c           get the update matrix equivalent stress and current
c           plastic modulus (sbar_new, h_new). the matrix response
c           can be very complex. see details of solution inside 
c           lower routine.
c
      call mm03sb( sbar_new, h_new, ebarp_new, ebarp, consth,
     &             sigma_o, h_constant, e, n_power,m_power,
     &             dt, eps_ref, rate_depen, mpoweri, iout, 
     &             debug, segmental, power_law, first_time )
c
c           for nucleation, get the A function and the derivative
c           of A wrt ebarp at the current estimate for ebarp_new.
c
      a_nuc   = zero
      a_prime = zero
      if ( nucleation ) then
        call mm03ap( a_nuc, a_prime, ebarp_new, nuc_s_n, nuc_e_n,
     &               nuc_f_n )
      end if
c
c           evaluate the residual of netwon loop and check
c           for convergence.
c
      d_ebarp = ebarp_new - ebarp
      term1   = dep + f + a_nuc*d_ebarp
      term2   = sbar_new * ( one - term1/(one + dep) )
      r       = d_ebarp - t_numer / term2
c
c          compute correction for matrix plastic strain and check for
c          convergnce of matrix plastic strain.
c
      term1    = (one-dep) * t_numer * ( a_nuc + d_ebarp * a_prime )
      term2    = (one - f - ebarp * a_nuc + ebarp_new * a_nuc) 
      term3    = (one+dep) * t_numer * h_new
      dr_ebarp = one - term1 / term2 / term2 / sbar_new +
     &           term3 / term2 / sbar_new / sbar_new
      eps_correction = r / dr_ebarp 
      if ( abs(eps_correction) .le. strain_toler ) then 
c
c           converged values for ebarp_new and sbar_new found. get
c           f_new and return. make sure f_new is positive or
c           the load increments are probably too big.
c
         f_new = (f + dep + a_nuc * d_ebarp) / ( one + dep )
         if ( f_new .lt. zero ) then 
            converge = .false.
            return
         end if
         if ( ebarp_new .lt. ebarp ) then
            converge = .false.
            return
         end if
         if ( debug ) then 
            write(iout,9100) ebarp_new, sbar_new, f_new, h_new,
     &                       a_nuc, a_prime,loop
         end if
         return
c
      end if
c
c          residual of netwon loop still too large. try again.
c
      ebarp_new = ebarp_new - eps_correction
c
      end do
c
c           iterations to find ebarp_new did not converge
c
      converge = .false.
      return
c
 9100 format(/,8x,' >> returning from mm03us:',
     &   /,12x,'ebarp_new, sbar_new: ',f10.7,f10.3,
     &   /,12x,'f_new, h_new:        ',f10.5,e10.3,
     &   /,12x,'A_nuc, a_prime:      ',f10.5,f10.5,
     &   /,12x,'iterations required: ',i5 )
 9200 format(/,6x,'>> Starting iteration: ',i3,
     &        ' to find ebarp, f, sbar...' )
c
      end
c ************************************************************************
c *                                                                      *
c *    routine  mm03ss -- update state variables f, ebarp, sbar given    *
c *                       estimates for dep, deq. a slower but better    *
c *                       version called if the simple-fast one fails    *
c *                                                                      *
c ************************************************************************
c
c
      subroutine mm03ss( p_new, dep, q_new, deq, ebarp, f, sbar,
     &                   consth, sigma_o, eps_o, h_constant, 
     &                   e, n_power, m_power, dt, eps_ref, rate_depen,
     &                   mpoweri, nucleation, nuc_s_n, nuc_e_n,
     &                   nuc_f_n, ebarp_new, sbar_new, f_new, h_new,
     &                   elenum, ptno, iout, debug, segmental,
     &                   power_law, converge, first_time )
#dbl      implicit double precision (a-h,o-z)
#dbl      double precision
#sgl      real
     &  n_power, m_power, mpoweri, nuc_s_n, nuc_e_n, nuc_f_n,
     &  sfactors(5)
      integer  elenum, ptno
      logical  local_debug, debug, consth, nucleation, rate_depen,
     &         segmental, power_law, converge, debug_mm03sb,
     &         first_time
      data    local_debug, debug_mm03sb
     &        / .false., .false. /
      data    zero, half, one, two, tol / 0.0, 0.5, 1.0, 2.0, 0.00001 /
      data    sfactors / 1.5, 2.0, 2.0, 4.0, 4.0 /
      data    max_iterations / 12 /
c
c           find ebarp_new and in the process we find f_new and sbar_new.
c           if nucleation is neglected the expressions are simpler.
c           using the development described on pages 16-17 of the
c           Gurson notes, the problem resolves to solving a
c           single nonlinear algebraic equation for ebarp_new.
c           the solution has proven fatal for a pure newton
c           method as there are "flat" spots that cause newton
c           failure.
c
c           here we use a combination of bisection-newton as outlined
c           in "numerical recipes" p. 359 (fortran, 2nd. ed.). bisection
c           simply keeps the newton method within the bounding limits
c           on ebarp_new and when convergence of newton seems slow,
c           a bisection is performed and netown resumed. so far
c           the algortihm has eliminated previous difficulties
c           encountered, especially when nucleation is invoked.
c
c
c           for names of variables and equations, just refer to the
c           gurson notes. this routine looks long because we make
c           a sequence of subroutine calls each time the "residual"
c           function in newton's method is required.
c    
      t1             = -p_new * dep + q_new * deq
      t3             = one + dep
      t2             = ( f + dep ) / t3
      ebarp_new      = ebarp + t1 / ( one - f ) / sbar
      strain_tol     = tol * max( eps_o, ebarp_new )
c
c           find the upper-bound on new ebarp - a new ebarp
c           large enough to change sign of the residual function.
c           assume current ebarp (start of step) is lower bound.
c           see comments below about calling mm03sb, mm03ap. 
c           increase estimate of upper-bound ebarp using ratios
c           stored in sfactors for convenience.
c
      ebarp_big =  ebarp
      do loop_bound = 1, 5
        call mm03sb( sbar_new, h_new, ebarp_big, ebarp, consth,
     &             sigma_o, h_constant, e, n_power, m_power,
     &             dt, eps_ref, rate_depen, mpoweri, iout,
     &             debug_mm03sb, segmental, power_law, first_time )
        a_nuc   = zero
        a_prime = zero
        if ( nucleation ) then
           call mm03ap( a_nuc, a_prime, ebarp_big, nuc_s_n, nuc_e_n,
     &                  nuc_f_n )
        end if
        d_ebarp = ebarp_big - ebarp
        term1   = sbar_new * ( one - t2 - a_nuc*d_ebarp/t3 )
        r       = d_ebarp - t1 / term1
        if ( local_debug ) then
          write(iout,9460) loop_bound, ebarp, ebarp_big, r
        end if
        if ( loop_bound .eq. 1 ) then
          r_low  = r
          if ( ebarp .eq. zero ) ebarp_big = eps_o
        else if ( r * r_low .le. zero ) then
          r_high = r
          go to 100
        end if
        ebarp_big = ebarp_big * sfactors(loop_bound)
      end do
c
c           we failed to bound ebarp_new - this is really bad!
c
c      write(iout,9400) elenum, ptno, ebarp, p_new, dep, q_new, deq,
c     &                 ebarp + t1 / ( one - f ) / sbar
      converge = .false.
      return
c
c           ebarp_new is bracketed between ebarp and ebarp_big.
c           use a combination of newton's method and bisection
c           to locate the value of ebarp_new. 
c
 100  continue
      loop = 0
      if ( r_low .eq. zero ) then
        ebarp_new = ebarp
        go to 1000
      else if ( r_high .eq. zero ) then
        ebarp_new = ebarp_big
        go to 1000
      else if ( r_low .lt. zero ) then
        ebarp_low = ebarp
        ebarp_high = ebarp_big
      else
        ebarp_high = ebarp
        ebarp_low = ebarp_big
      end if
      ebarp_new = ( ebarp_low + ebarp_high ) * half
      ebarp_old = ebarp_new
      deps_old  = abs( ebarp_low - ebarp_high )
      deps      = deps_old
      if ( local_debug ) then
        write(iout,9420) r_low, r_high, ebarp_big, ebarp_high,
     &                   ebarp_low, ebarp_new, deps, ebarp,
     &                   loop_bound-1
      end if
c
c           get the updated matrix equivalent stress and current
c           plastic modulus (sbar_new, h_new). the matrix response
c           can be very complex. see details of solution inside 
c           lower routine.
c
      call mm03sb( sbar_new, h_new, ebarp_new, ebarp, consth,
     &             sigma_o, h_constant, e, n_power,m_power,
     &             dt, eps_ref, rate_depen, mpoweri, iout,
     &             debug_mm03sb, segmental, power_law, first_time )
c
c           for nucleation, get the A function and the derivative
c           of A wrt ebarp at the current estimate for ebarp_new.
c
      a_nuc   = zero
      a_prime = zero
      if ( nucleation ) then
        call mm03ap( a_nuc, a_prime, ebarp_new, nuc_s_n, nuc_e_n,
     &               nuc_f_n )
      end if
c
c           evaluate the "residual" function and it's derivative
c           at the current estimate of ebarp_new
c
      d_ebarp  = ebarp_new - ebarp
      term1    = sbar_new * ( one - t2 - a_nuc*d_ebarp/t3 )
      r        = d_ebarp - t1 / term1
      term1    = t1 * ( -a_nuc/t3 - a_prime*d_ebarp/t3 )
      term2    = sbar_new * (one - t2 - a_nuc*d_ebarp/t3)**2
      term3    = t1 * h_new
      term4    = sbar_new * sbar_new * (one - t2 - a_nuc*d_ebarp/t3)
      dr_ebarp = one + term1/term2 + term3/term4
      if ( local_debug ) then
         write(iout,9430) ebarp_new, sbar_new, h_new, a_nuc, a_prime,
     &                    r, dr_ebarp
      end if
c
c           we must perform the bisection-newtion loop to improve
c           current estimate of ebarp_new
c
      do loop = 1, max_iterations
      if ( local_debug ) write(iout,9200) loop
c
      term1 = ( ebarp_new - ebarp_high ) * dr_ebarp - r  
      term2 = ( ebarp_new - ebarp_low  ) * dr_ebarp - r  
      if ( term1*term2 .ge. zero .or.
     &     abs(two*r) .gt. abs(deps_old*dr_ebarp) ) then
        deps_old = deps
        deps = half * ( ebarp_high - ebarp_low ) 
        ebarp_new = ebarp_low + deps
        if ( local_debug ) then
          write(iout,9440) loop, term1, term2, deps, ebarp_low,
     &                     ebarp_new, r
        end if
        if ( ebarp_low .eq. ebarp_new ) go to 1000
      else
        deps_old = deps
        deps = r / dr_ebarp
        temp = ebarp_new
        ebarp_new = ebarp_new - deps
        if ( local_debug ) then
          write(iout,9450) loop, term1, term2, deps, r, ebarp_new
        end if
        if ( temp .eq. ebarp_new )  go to 1000
      end if
c
c           convergence check on ebarp_new
c
      if ( abs(ebarp_new-ebarp_old) .le. strain_tol ) go to 1000
      ebarp_old = ebarp_new
c
c           get the update matrix equivalent stress and current
c           plastic modulus (sbar_new, h_new). the matrix response
c           can be very complex. see details of solution inside 
c           lower routine.
c
      call mm03sb( sbar_new, h_new, ebarp_new, ebarp, consth,
     &             sigma_o, h_constant, e, n_power,m_power,
     &             dt, eps_ref, rate_depen, mpoweri, iout,
     &             debug_mm03sb, segmental, power_law, first_time )
c
c           for nucleation, get the A function and the derivative
c           of A wrt ebarp at the current estimate for ebarp_new.
c
      a_nuc   = zero
      a_prime = zero
      if ( nucleation ) then
        call mm03ap( a_nuc, a_prime, ebarp_new, nuc_s_n, nuc_e_n,
     &               nuc_f_n )
      end if
c
c           evaluate the residual of netwon loop and the derivative.
c           set new lower and upper bounds.
c
      d_ebarp  = ebarp_new - ebarp
      term1    = sbar_new * ( one - t2 - a_nuc*d_ebarp/t3 )
      r        = d_ebarp - t1 / term1
      term1    = t1 * ( -a_nuc/t3 - a_prime*d_ebarp/t3 )
      term2    = sbar_new * (one - t2 - a_nuc*d_ebarp/t3)**2
      term3    = t1 * h_new
      term4    = sbar_new * sbar_new * (one - t2 - a_nuc*d_ebarp/t3)
      dr_ebarp = one + term1/term2 + term3/term4
      if ( r .lt. zero ) then
         ebarp_low  = ebarp_new
      else
         ebarp_high = ebarp_new
      end if
c
      end do
c
c           iterations to find ebarp_new did not converge
c
      converge = .false.
      write(iout,9000) elenum, ptno
      return
c
 1000 continue 
c
c           converged values for ebarp_new and sbar_new found. get
c           f_new and return. make sure f_new is positive or
c           the load increments are probably too big. if loop = 0
c           we skipped the newton-iteration after the bracketing.
c           need to update state for ebarp_new.
c
      converge = .true.
      if ( loop .eq. 0 ) then 
       call mm03sb( sbar_new, h_new, ebarp_new, ebarp, consth,
     &             sigma_o, h_constant, e, n_power,m_power,
     &             dt, eps_ref, rate_depen, mpoweri, iout, 
     &             debug_mm03sb, segmental, power_law, first_time )
       a_nuc   = zero
       a_prime = zero
       if ( nucleation ) then
         call mm03ap( a_nuc, a_prime, ebarp_new, nuc_s_n, nuc_e_n,
     &                nuc_f_n )
       end if
      end if
      d_ebarp  = ebarp_new - ebarp
      f_new = t2 + a_nuc * d_ebarp / t3
      if ( f_new .lt. zero ) then 
          write(iout,9300) elenum, ptno
          converge = .false.
      end if
      if ( debug ) then 
         write(iout,9100) ebarp_new, sbar_new, f_new, h_new,
     &                       a_nuc, a_prime, loop
      end if
      return
c
 9000 format(10x,i7,i3,' newton loop in mm03ss failed to converge') 
 9100 format(/,8x,' >> returning from mm03ss:',
     &   /,12x,'ebarp_new, sbar_new: ',f10.7,f10.3,
     &   /,12x,'f_new, h_new:        ',f10.5,e10.3,
     &   /,12x,'A_nuc, a_prime:      ',f10.5,f10.5,
     &   /,12x,'iterations required: ',i5 )
 9200 format(/,6x,'>> Starting iteration: ',i3,
     &        ' to find ebarp, f, sbar...' )
 9300 format(10x,i4,i3,' f_new < 0 during update') 
 9400 format(10x,i7,i3,' failed to bracket ebarp_new.',
     &   /,12x,'ebarp, p_new, dep, q_new, deq: ',f12.8, f10.3,
     &         f12.8, f10.3, f12.8,
     &   /,12x,'starting estimate for ebarp_new: ',f12.8 )
 9420 format(/,8x,' >> mm03ss. ebarp_new bracketed., setup of',
     &       /,8x,'    newton-bisection iterations...',
     &       /,8x,'    r_low, r_high: ',2f10.6,
     &       /,8x,'    ebarp_big, ebarp_high, ebarp_low: ',3f12.9,
     &       /,8x,'    ebarp_new, deps, ebarp: ',3f12.9,
     &       /,8x,'    iterations to find bound: ',i4 )
 9430 format(/,8x,' >> ready to start newton-bisection loop:',
     &       /,8x,'    ebarp_new, sbar_new, h_new: ',f12.9,2f10.3,
     &       /,8x,'    a_nuc, a_prime, r, dr_ebarp: ',4f12.6)
 9440 format(
     & /,8x,'     > @ 1: loop, term1, term2, deps: ',i2,2f10.3,f13.9,
     & /,8x,'            ebarp_low, ebarp_new, r: ',2f12.6,f10.6 )
 9450 format(
     & /,8x,'     > @ 2: loop, term1, term2, deps: ',i2,2f10.3,f13.9,
     & /,8x,'            r, ebarp_new: ',f10.6,f12.6 )
 9460 format(/,8x,'  > loop_bound, ebarp, ebarp_big, r: ',
     & i2, 2f13.9, f10.5 )
c
      end
c ************************************************************************
c *                                                                      *
c *    routine  mm03sb -- update the matrix equivalent stress - sbar -   *
c *                       given a value of the total plastic strain      *
c *                       in the matrix. the matrix can be inviscid or   * 
c *                       visco-plastic                                  * 
c *                                                                      *
c ************************************************************************
c
      subroutine mm03sb( sbar_new, h_new, ebarp_new, ebarp, consth,
     &                   sigma_o, h_constant, e, n_power, 
     &                   m_power, dt, eps_ref, rate_depen, mpoweri, 
     &                   iout, debug, segmental, power_law,
     &                   first_time )
#dbl      implicit double precision (a-h,o-z)
#dbl      double precision
#sgl      real
     &   n_power, mpoweri, mm03is, m_power, mm03sc
      logical consth, debug, rate_depen, segmental, power_law,
     &        first_time
      data one / 1.0 / 
c
c                 1. get the updated inviscid equivalent stress
c                    and plastic modulus given the total plastic
c                    strain in the matrix. matrix can have
c                    constant hardening, a linear + power-law or
c                    a piecewise linear (segmental).
c                    for linear + power law or segmental,
c                    call a lower-level routine which iterates
c                    to find new stress. the segmental response can
c                    be rate dependent thru a set of user input curves
c                    of stress vs. plastic strain at various plastic
c                    strain rate. mm03sc returns the rate dependent
c                    plastic modulus only when passed a case=2 since this
c                    requires some significant work.
c
      if ( consth ) then
           sbar_inviscid = sigma_o + h_constant * ebarp_new
           h_inviscid    = h_constant
      else if (power_law) then
           sbar_inviscid = mm03is( ebarp_new, sigma_o, e, n_power,
     &                             h_inviscid, iout )
        else if ( segmental ) then
           d_ebarp = ebarp_new - ebarp
           sbar_inviscid = mm03sc( ebarp_new, h_inviscid, d_ebarp,
     &                             dt, 2 )   
      end if
c
c                 2. if matrix is power-law rate sensitive, scale inviscid
c                    equivalent stress using power-law viscosity. 
c                    convert inviscid plastic modulus to the
c                    rate sensitive value.
c
      if ( rate_depen ) then
          d_ebarp = ebarp_new - ebarp
          rate_multiplier = ( eps_ref * d_ebarp/dt + one )**mpoweri
          sbar_new = sbar_inviscid * rate_multiplier
          term0 = ( sbar_new / sbar_inviscid ) ** (one - m_power )
          term1 = term0 * eps_ref * sbar_inviscid / m_power / dt
          term2 = h_inviscid * sbar_new / sbar_inviscid
          h_new = term1 + term2
      else
          sbar_new = sbar_inviscid
          h_new    = h_inviscid
      end if
c
      if ( .not. debug ) return
c
      write(iout,*) ' '
      write(iout,*) '        > mm03sb: update sbar, h'
      write(iout,9000)  ebarp_new, ebarp, sbar_new, h_new
      write(iout,9100)  sigma_o, e, rate_depen
      if ( consth ) then
         write(iout,9110) h_constant
      else
         write(iout,9120) n_power
      end if
      if ( rate_depen ) then
         write(iout,9200) m_power, dt, eps_ref, mpoweri,
     &                    rate_multiplier, sbar_inviscid, h_inviscid
      end if
      return
c
 9000 format(
     & /,15x,'ebarp_new, ebarp                :',2f10.6,
     & /,15x,'sbar_new, h_new:                :',f10.3,e10.3 )
 9100 format(
     &   15x,'sigma_o, e, rate_depen          :',2f10.3,l3 )
 9110 format(
     &   15x,'h_constant                      :',e10.3 )
 9120 format(
     &   15x,'n_power                         :',f10.3 )
 9200 format(
     & /,15x,'mpower, dt, eps_ref, mpoweri    :',4f10.3,
     & /,15x,'rate_mult., sbar_inviscid       :',2f10.3,
     & /,15x,'h_inviscid                      :',e10.3 )
c
      end
c ************************************************************************
c *                                                                      *
c *    routine  mm03is -- linear + power law model for uniaxial          *
c *                       response (cubic transition)                    * 
c *                                                                      *
c ************************************************************************
c
c
      function mm03is( eps_pls, sigma_o, e, power, hprime, iout )
c
#dbl      implicit double precision (a-h,o-z)
#dbl      double precision
#sgl      real
     & mm03is, n1, n2, n3, n4
      logical local_debug, converge
      data local_debug / .false. /
      data zero, one, two, three / 0.0, 1.0, 2.0, 3.0 /
      data beta, alpha, toler, max_itr
     &     / 0.95, 1.1, 0.000001, 20 /
c
c              given the total plastic strain, return the stress,
c              and new plastic modulus. an iterative solution is
c              required to find the stress. a newton iteration
c              scheme converges very fast.
c
      if ( eps_pls .le. zero ) then
       mm03is = sigma_o
       hprime = e * (0.9999*e) / ( e-0.9999*e)
       return
      end if 
c
c              set factors that define start and end values of strain
c              (relative to eps_o) for the cubic transition between the
c              linear and power-law regions. the sigma_o value
c              passed is the unreduced value of yield-stress used in
c              yield function evaluations, etc.
c              the small transition region is defined by factors
c              beta and alpha. this is needed as we define the transition
c              part of curve relative to the user specified yield point.
c              'a' and 'b' refer to the start and end points of the cubic
c              transition. changing beta and alpha can destroy the validity
c              of the cubic transition.  see notes of derivation.
c
      sig_o   = sigma_o / beta
      eps_o   = sig_o / e
      poweri  = one / power
      eps_a   = beta  * eps_o
      eps_b   = alpha * eps_o
      h       = eps_b - eps_a
      h2      = h * h
      h3      = h2 * h
      sig_a   = beta * sig_o
      sig_b   = sig_o * alpha**poweri
      epower  = e / power
      pm1     = poweri-one
      et_b    = epower * alpha**pm1
c
      iterno  = 1
      strain  = eps_pls + eps_o
      stress  = zero
c
c             use Newton iterations to correct starting value.
c
400   continue 
      if ( strain .ge. eps_b ) then
         stress_new = sig_o * ( strain/eps_o )**poweri  
         fprime     = epower * ( strain/eps_o )**pm1  
      else
         ta = strain - eps_a
         tb = strain - eps_b
         n1 = -tb*tb * (-h - two*ta) / h3
         n2 =  ta * tb * tb / h2
         n3 =  ta * ta * (h - two*tb) / h3
         n4 =  ta * ta * tb / h2
         stress_new = n1*sig_a + n3*sig_b + n2*e + n4*et_b
         pn1 = two*(tb)*(three*strain-two*eps_a-eps_b+h)/h3
         pn3 = two*(ta)*(-three*strain+eps_a+two*eps_b+h)/h3
         pn2 = -tb*(-three*strain+two*eps_a+eps_b)/h2
         pn4 = -ta*(-three*strain+eps_a+two*eps_b)/h2
         fprime = pn1*sig_a + pn3*sig_b + pn2*e + pn4*et_b
      end if
      resid      = strain - stress_new/e - eps_pls
      deps       = -resid / ( one - fprime/e )
      strain_new = strain + deps
c      if ( local_debug ) then
c        write (iout,9020) iterno, strain, stress_new, fprime,
c     &                    resid, deps, strain_new
c      end if  
      converge =  abs( strain_new - strain ) .lt. toler * strain_new
     &                              .and. 
     &            abs( stress_new - stress ) .lt. toler * stress_new          
      if ( .not. converge ) then
          iterno = iterno + 1
          if ( iterno .gt. max_itr ) then
            mm03is = sig_o
            hprime = e * e
            return
          end if
          strain = strain_new
          stress = stress_new
          go to 400
      end if                                      
c
c             newton iteration converged. set final stress
c             and compute tangent modulus, plastic modulus
c             from that value.
c
      et     = fprime
      hprime = ( e * et ) / ( e - et ) 
      mm03is =  stress_new
c
      return
c
 9020 format (1x,'>>  Iteration number : ',i2,
     &   /,10x, 'strain, stress_new, fprime: ',f10.6,f10.3,f10.3,
     &   /,10x, 'resid, deps, strain_new:    ',f10.6,f10.6,f10.6 )
c
      end
cc *******************************************************************
c *                                                                 *
c *      material model # 03 routine -- mm03ap                      *
c *                                                                 *
c *******************************************************************
c
c
      subroutine mm03ap( a_nuc, a_prime, ebarp, nuc_s_n, nuc_e_n,
     &                   nuc_f_n )
#dbl      implicit double precision (a-h,o-z)
#dbl      double precision
#sgl      real
     &   nuc_s_n, nuc_e_n,  nuc_f_n
      data root_2_pi / 2.50663 /
      data half / 0.5 /
c
      term1   = nuc_f_n / nuc_s_n / root_2_pi
      term2   = ( ( ebarp - nuc_e_n ) / nuc_s_n )**2
      a_nuc   = term1 * exp(-half*term2)
      a_prime = -nuc_f_n * (ebarp-nuc_e_n) / nuc_s_n**3 / root_2_pi
      a_prime = a_prime * exp(-half*term2) 
c
      return
      end
c *******************************************************************
c *                                                                 *
c *      material model # 03 routine -- mm03dr                      *
c *                                                                 *
c *******************************************************************
c
c
      subroutine mm03dr( phi, gp, gq, gsbar, gf, q, sbar,
     &                   p, f, ebarp, f_old, ebarp_old, nucleation, 
     &                   shear_mod, bulk_mod, h, q1,
     &                   q2, q3, nuc_s_n, nuc_e_n, nuc_f_n, dep, deq, 
     &                   d11, d12, d21, d22, debug, iout )
c
c              compute the value of the Gurson yield function (phi)
c              given the scaler parameters. compute derivatives
c              of the yield function wrt to the state variables.
c              compute partial derivatives of the 2 residual
c              functions (r1, e2) wrt dep and deq.
c
#dbl      implicit double precision (a-h,o-z)
#dbl      double precision
#sgl      real
     &  mm03g, nuc_s_n, nuc_e_n, nuc_f_n
      logical nucleation, debug
      data onehalf / 1.5 /
      data one, two, three, four / 1.0, 2.0, 3.0, 4.0 /
      data zero / 0.0 /
c
c              1. the yield function is call phi or g.
c
      phi = mm03g( q, sbar, f, p, q1, q2, q3 )
c
c             1a. the hyperbolic sin and cos of the argument pressure
c                 term that appears repeatedly.
c
      beta = -onehalf * q2 
      ch   = cosh( beta*p/sbar )
      sh   = sinh( beta*p/sbar )
c
c              2. partial g / partial p
c
      gp = two * beta * f * q1 * sh / sbar
c
c              3. partial g / partial q
c
      gq = two * q / sbar**2
c
c              4. partial g / partial sbar
c
      gsbar = (two / sbar**3 ) *
     &        ( -q*q - beta * f * p * q1 * sbar * sh )
c
c              5. partial g / partial f
c
      gf = two * ( -f * q3 + q1 * ch)
c
c              6. for nucleation, get the A function and the
c                 derivative of 'a' wrt ebarp at the current
c                 estimate for ebarp_new.
c
      a_nuc   = zero
      a_prime = zero
      if ( nucleation ) then
        call mm03ap( a_nuc, a_prime, ebarp, nuc_s_n, nuc_e_n,
     &               nuc_f_n )
      end if
c
c              7. increment of matrix plastic strain.
c
      d_ebarp = ebarp - ebarp_old
c
c              8. partial f / partial dep, 
c                 partial d_ebarp / partial dep
c
      a1 = one / ( one + dep ) - 
     &     ( dep + f_old + d_ebarp * a_nuc ) / ( one + dep )**2
      a1 = -a1
      a2 = d_ebarp * a_prime  / ( one + dep )
      b1 = ( -p * dep + q * deq ) / sbar / ( one - f )**2
      b2 = ( -p - bulk_mod * dep ) / sbar / ( one - f )
      b3 = ( -p * dep + q * deq ) * h / ( one - f ) / sbar**2
c
      df_dep     = ( a1 - a2 * b2 + a1 * b3 ) / ( a2 * b1 - b3 - one )
      debarp_dep = ( b2 - a1 * b1 ) / ( one - a2 * b1 + b3 )
c
c              9. partial f / partial deq, 
c                 partial d_ebarp / partial deq
c
      b4 = d_ebarp * a_prime / ( one + dep )
      b5 = ( q - three * shear_mod * deq ) / ( one - f ) / sbar
c
      df_deq     = b4 * b5 / ( one + b3 - b1 * b4 )
      debarp_deq = b5 / ( one + b3 - b1 * b4 )
c
c             10. d21 = partial r2 / partial dep
c                 d22 = partial r2 / partial deq
c
      d21 = gp * bulk_mod   +    gf * df_dep   +   
     &      h * gsbar * debarp_dep
      d22 = -three * shear_mod * gq + h * gsbar * debarp_deq +
     &      gf *  df_deq 
c
c             11. d11 = partial r1 / partial dep
c
      beta  = -onehalf * q2
      alpha = two * beta * q1
      term1 = two * q / sbar**2
      term2 = alpha * deq * sh * df_dep / sbar
      term3 = ( ( four * dep * q / sbar**3  ) + 
     &        alpha * deq * f * sh / sbar**2 ) * h * debarp_dep
      term4 = alpha * beta * f * deq * ch *
     &        ( bulk_mod * sbar - p * h * debarp_dep ) / sbar**3
      d11   = term1 + term2 - term3 + term4
c
c             12. d12 = partial r1 / partial deq
c
      term1 = three * three * shear_mod * dep / sbar**2
      term2 = alpha * f * sh / sbar
      term3 = alpha * deq * sh * df_deq / sbar
      term4 = ( four * dep * q / sbar**3 + 
     &        alpha * deq * f * sh / sbar**2 ) *
     &        h * debarp_deq
      term5 = alpha * beta * deq * f * p * ch * h *
     &        debarp_deq / sbar**3
      d12   = -term1 + term2 + term3 - term4 - term5
c
      if ( debug ) then
        write(iout,9100) phi, gp, gq, gsbar, gf,
     &                   d_ebarp, a_nuc, a_prime,
     &                   df_dep, debarp_dep, df_deq, debarp_deq,
     &                   d11, d12, d21, d22 
        write(iout,*) ' '
      end if
c
      return
c
 9100 format(/,8x,' >> returning from mm03dr:',
     &   /,12x,'phi                          : ',f10.4,
     &   /,12x,'gp, gq, gsbar, gf            : ',4f10.4,
     &   /,12x,'d_ebarp, a_nuc, a_prime      : ',f10.6,2f10.4,
     &   /,12x,'df_dep, debarp_dep           : ',2f10.4,
     &   /,12x,'df_deq, debarp_deq           : ',2f10.4,
     &   /,12x,'d11, d12                     : ',2e12.4,
     &   /,12x,'d21, d22                     : ',2e12.4 )
c
      end
c ****************************************************************
c *                                                              *
c *      function to lookup stress for given plastic strain on   *
c *      segmentally defined stress-strain curve                 *
c *                                                              *
c ****************************************************************
c
c
      function mm03sc( ebarp, hprime, deplas, dtime, caseh )
c
c                      parameter declarations
c
      use segmental_curves, only :  now_blk_relem, curve_plseps_values,
     &                              seg_curve_table, active_curve_set,
     &                              sigma_curves, num_seg_points,
     &                              max_seg_points, seg_curves_type,
     &                              curve_rates, sigma_inter_table
#dbl      implicit double precision (a-h,o-z)
#dbl      double precision mm03sc
#sgl      real mm03sc
      integer caseh     
c
c
c                      local declarations
c
      integer first_curve, curve_set_type, pt_high, pt_low
#dbl      double precision 
#sgl      real 
     & stress_values(max_seg_points), mm03lint,
     & curve_high(max_seg_points), curve_low(max_seg_points)
c
      data zero, one, deriv_factor, two / 0.0, 1.0, 0.02, 2.0 /
c
      first_curve       = seg_curve_table(2,active_curve_set)
      curve_set_type    = seg_curves_type(first_curve)
      num_curves_in_set = seg_curve_table(1,active_curve_set)
      numpts            = num_seg_points(first_curve)
c
c               check for input point before first point on curve
c               or beyond last point on curve.
c
      if ( curve_set_type .eq. 2 ) then
        eps_pls_dot = zero
        if ( dtime .gt. zero ) eps_pls_dot =  deplas / dtime
        do i = 1, numpts
         stress_values(i) = mm03lint( eps_pls_dot, num_curves_in_set,
     &                                curve_rates,
     &                                sigma_inter_table(1,i) ) 
        end do               
      else
        do i = 1, numpts
          stress_values(i) = sigma_curves(i,now_blk_relem)
        end do
      end if
c      
      if ( ebarp .lt. curve_plseps_values(1) ) then
        mm03sc = stress_values(1)
        hprime = 1.0e10
        return
      end if
c
c                 check for input point beyond right end of curve
c
      if ( ebarp .ge. curve_plseps_values(numpts) ) then
        mm03sc = stress_values(numpts)
        hprime = zero
        return
      end if
c
c            point is actually on the curve between points i and i-1.
c            interpolate linearly to get stress. set plastic modulus. 
c
      do i = 2, numpts
        if ( ebarp .lt. curve_plseps_values(i) ) then
          sig_high = stress_values(i)
          sig_low  = stress_values(i-1)
          eps_high = curve_plseps_values(i)
          eps_low  = curve_plseps_values(i-1) 
          hprime   = (sig_high - sig_low) / ( eps_high - eps_low )
          mm03sc   = sig_low + hprime * ( ebarp - eps_low )
          pt_high  = i
          pt_low   = i-1
          go to 1000
        end if
      end do
c
c            if we get here things are very, very bad. stop job
c
      write(*,*) '>> FATAL ERROR: routine mm03sc reached an impossible'
      write(*,*) '>>              condition. job terminated.'
      write(*,*) '>> ebarp: ', ebarp
      write(*,*) '>> first_curve, numpts: ', first_curve, numpts
      call abort_job
c
c            for rate dependent, segmental material, we must modify 
c            the hprime computed above to include the effects
c            of the change in equivalent stress at this plastic strain
c            due to a change in plastic strain rate. do only if calling
c            routine really needs value since this requires some work.
c
 1000 continue
      if ( curve_set_type .ne. 2 ) return 
      if ( caseh .eq. 1 ) return
      if ( deplas .le. zero ) return
c      
c                      for rate dependent response, we now know
c                      hprime for the current plastic strain rate.
c                      the hprime value for use in computation of the
c                      consistent [D] matrix needs another term.
c                      this term provides  the rate of change of stress
c                      w.r.t. plastic strain rate at this plastic
c                      strain level. 
c
c     d sigbar = ( partial sigbar / partial ebarp ) * d ebarp +
c                ( partial sigbar / partial eps_plas_dot ) *  d eps_plas_dot
c
c     but d eps_plas_dot = d eps_plas / d time
c
c     so
c  
c     d sigbar = ( partial sigbar / partial eps_plas  +
c                  (1/ d time ) partial sigbar / partial eps_plas_dot )
c                 *  d eps_plas
c
c                      this defines the modfied, rate dependent hprime
c                      for use in consistent [D] computation
c
c     d sigbar = ( modified hprime ) * d eps_plas
c
c                      since the rate dependent, uniaxial stress vs.
c                      plastic strain response is available only in
c                      segmental form, we must compute an approximate
c                      derivative,  partial sigbar / partial eps_plas_dot,
c                      using a central differnce method about the
c                      current plastic strain rate (passed in).
c
c                      the code below constructs interpolated segmental
c                      stress vs. plastic strain curves at a slightly
c                      larger and slower plastic strain rate.
c                      The stress value on those two curves at the
c                      current plastic strain (include mixed hardening
c                      beta effect) gives the data to approximate
c                      partial sigbar / partial eps_plas_dot. 
c
c
      rate_high = (one + deriv_factor ) * eps_pls_dot
      rate_low  = (one - deriv_factor ) * eps_pls_dot
c
c                      construct the two interpolated stress vs.
c                      plastic strain curves at slightly larger and
c                      lower plastic strain rates. 
c
      do i = 1, numpts
        curve_high(i) = mm03lint( rate_high, num_curves_in_set,
     &                            curve_rates, sigma_inter_table(1,i) ) 
        curve_low(i)  = mm03lint( rate_low, num_curves_in_set,
     &                            curve_rates, sigma_inter_table(1,i) ) 
      end do
c
c                      compute hprime for each of the two slightly
c                      different plastic rate curves
c
      hprime_high = ( curve_high(pt_high) - curve_high(pt_low) ) /
     &              ( eps_high - eps_low )
      hprime_low  = ( curve_low(pt_high) - curve_low(pt_low) ) /
     &              ( eps_high - eps_low ) 
      dep = ebarp - eps_low 
      stress_high = curve_high(pt_low) +  hprime_high * dep
      stress_low  = curve_low(pt_low)  +  hprime_low  * dep
c
c                      get approximate:
c                           partial sigbar / partial eps_plas_dot
c                      and modify the passed in hprime value for
c                      use in consistent [D] computation
c
      derivative  = (stress_high - stress_low ) /
     &              (two * deriv_factor * eps_pls_dot)
      hprime      = hprime + derivative / dtime
      return
c
      end
c     ****************************************************************
c     *                                                              *
c     *                        function mm03lint                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 05/2/00                     *
c     *                                                              *
c     *    execute linear interpolation on a tabular function        *
c     *    of a single variable where the x values are sorted in     *
c     *    increasing value but are not req'd to be uniformly spaced *
c     *                                                              *
c     ****************************************************************
c

      function mm03lint( xvalue, n, x, y )
      implicit integer (a-z)
#dbl      double precision
#sgl      real
     & xvalue, x(n), y(n), x1, x2, y1, y2, mm03lint
c
      if ( xvalue .le. x(1) ) then
        mm03lint = y(1)
        return
      end if
c
      if ( xvalue .ge. x(n) ) then
        mm03lint = y(n)
        return
      end if
c
      do point = 2, n
        if ( xvalue .gt. x(point) ) cycle
        x1 = x(point-1)
        x2 = x(point)
        y1 = y(point-1)
        y2 = y(point)
        mm03lint = y1 + (xvalue-x1)*(y2-y1)/(x2-x1)
        return
      end do
c
      write(*,*) '>> FATAL ERROR: mm03lint'
      write(*,*) '                job aborted'
      write(*,*) '>> xvalue: ',xvalue
      do i = 1, n
        write(*,*) '  i,x,y: ',i,x(i),y(i)
      end do
      stop
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm03_set_sizes                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 1/11/2015  rhd              *
c     *                                                              *
c     *    called by warp3d for each material model to obtain        *
c     *    various sizes of data for the model                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm03_set_sizes( info_vector )
      dimension info_vector(*)
c
c        set infor_data
c
c         1        number of history values per integration 
c                  point. Abaqus calles these "statev". Values
c                  double or single precsion based on hardware.
c    
c         2        number of values in the symmetric part of the 
c                  [D] for each integration point. for solid
c                  elements this is 21, for cohesive elements this 6.
c
c         3        = 0, the material model returns "unrotated"
c                       Cauchy stresses at n+1
c                  = 1, the material model returns the standard
c                       Cauchy stresses at n+1
c
c         4        number of state variables per point to be output
c                  when user requests this type of results
c
      info_vector(1) = 11
      info_vector(2) = 21
      info_vector(3) = 0
      info_vector(4) = 6
c
      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm03_states_values                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 1/11/15 (rhd)                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm03_states_values( itype, elem_states_output,
     &                                nrow_states, num_states  )
c
c                       access some global data structures
c
      use elem_block_data, only: history_blocks, history_blk_list
      use main_data, only: elems_to_blocks
c      
      implicit integer (a-z)
$add common.main
c
c                       parameters
c
      integer :: nrow_states, itype, num_states
#dbl      double precision :: elem_states_output(nrow_states,*)
#sgl      real  :: elem_states_output(nrow_states,*)
c
c                       locals
c
#dbl      double precision, 
#sgl      real,
     & allocatable :: history_dump(:,:,:), one_elem_states(:)
      integer :: relem, elnum, hist_size, blockno
      logical :: do_a_block
      double precision :: zero
      data zero / 0.0d00 /
c      
c           build deformation plasticity states values output.
c
c              itype > 0 => this is the block number. do all elements
c                           in the block
c
c              itype < 0 => this is an element number. put state
c                           values into column 1 of results.
c 
      do_block = .true.
      if( itype. gt. 0 ) then
         do_a_block = .true.
         blockno = itype
      else
         do_a_block = .false.
         elnum = -itype
         blockno = elems_to_blocks(elnum,1)
      end if          
c
      local_debug = .false.
      felem       = elblks(1,blockno)
      elem_type   = iprops(1,felem)
      mat_type    = iprops(25,felem)
      int_points  = iprops(6,felem)
      span        = elblks(0,blockno)
      hist_size   = history_blk_list(blockno)
      
      if( local_debug ) write(out,9050) blockno, felem, elem_type,         
     &         mat_type, int_points, span, hist_size
c
c           temporary block of history so it can be re-organized
c
      allocate( one_elem_states(nrow_states) )
      allocate( history_dump(hist_size,int_points,span) )
      history_dump = reshape( history_blocks(blockno)%ptr,
     &           (/hist_size,int_points,span/) )
c      
      if( do_a_block ) then    
        do relem = 1, span
           elnum = felem + relem - 1  ! absolute element number
           one_elem_states(1:nrow_states) = zero 
           call mm03_states_values_a
           elem_states_output(1:nrow_states,relem) = 
     &                one_elem_states(1:nrow_states)
        end do
      else
        relem = elnum + 1 - felem
        one_elem_states(1:nrow_states) = zero
        call mm03_states_values_a
        elem_states_output(1:nrow_states,1) =
     &                one_elem_states(1:nrow_states)
      end if  
c        
      deallocate( history_dump, one_elem_states )
c        
      return
c      
 9050 format(10x,"block, felem, etype, mtype:  ",4i7,
     &  /,10x,   "int_pts, span, hist_size:    ",3i7 )
c 
      contains
c     ========      
     
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm03_states_values_a              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/11/15 (rhd)              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm03_states_values_a
c   
      implicit none      
c
c                       locals
c
      integer :: ipt   
#dbl      double precision :: 
#sgl      real ::
     & epspls, state, dword, sigeff, hprime, porosity_f, rate_yield 
       integer :: iword(2)
       equivalence ( dword, iword )
c
      epspls     = zero
      sigeff     = zero
      hprime     = zero
      porosity_f = zero
      state      = zero
      rate_yield = zero
c       
      do ipt = 1, int_points
        epspls = epspls + history_dump(1,ipt,relem)
        sigeff = sigeff + history_dump(2,ipt,relem)
        hprime = hprime + history_dump(4,ipt,relem)
        porosity_f = porosity_f + history_dump(5,ipt,relem)
        dword  = history_dump(6,ipt,relem)
        state  = state + dble(iword(1))
        rate_yield = rate_yield + history_dump(9,ipt,relem)
      end do
c
      one_elem_states(1) = epspls / dble(int_points)   
      one_elem_states(2) = sigeff / dble(int_points)
      one_elem_states(3) = hprime / dble(int_points)
      one_elem_states(4) = porosity_f / dble(int_points)
      one_elem_states(5) = state / dble(int_points)
      one_elem_states(6) = rate_yield / dble(int_points)
c
      return     
c
      end subroutine mm03_states_values_a
      end subroutine mm03_states_values
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm03_states_labels                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/11/15 (rhd)              *
c     *                                                              *
c     ****************************************************************
c
      
      subroutine mm03_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
      implicit none
c
c                       parameters
c
      integer :: size_state, num_states, out, max_comment_lines,
     &           num_comment_lines
      character(len=8)  :: state_labels(size_state)
      character(len=60) :: state_descriptors(size_state)
      character(len=80) :: comment_lines(max_comment_lines)
c
c                       locals
c
      integer :: i
c
      num_states = 6    
      state_labels(1) = "epspls"
      state_labels(2) = "sig_flow"
      state_labels(3) = "H'"
      state_labels(4) = "porosity"
      state_labels(5) = "state"
      state_labels(6) = "q"
            
      state_descriptors(1) = "Eq. plastic strain"
      state_descriptors(2) = "See states_header file"
      state_descriptors(3) = "Plastic modulus"
      state_descriptors(4) = "Porosity (f)"
      state_descriptors(5) = "See states_header file"
      state_descriptors(6) = "See states header file"
c   

c      
      comment_lines(1) = "Notes on Mises material model:"
      comment_lines(2) = "  sig_flow - current equivalent flow stress"
      comment_lines(3) = "  porosity - 0.0"
      comment_lines(4) = "  state: -1 actively linear elastic"
      comment_lines(5) = "          1 actively plastic"
      comment_lines(6) = "          value is average over" //
     &                          " integration points"
      comment_lines(7) = "  q - not defined for mises"
c  
      i = 8
      comment_lines(i+0) = " "
      comment_lines(i+1) = "Notes on Gurson material model:"
      comment_lines(i+2) = "  sig_flow - dense matrix equivalent stress"
      comment_lines(i+3) = "  state - same as above for mises"
      comment_lines(i+4) = "  q - mises equiv stress for continuum"
c      
      num_comment_lines = 12
c   
      return
      end




