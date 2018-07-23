c ************************************************************************
c *                                                                      *
c *    routine  mm03p -- vectorized pre-processor for block of mises or  *
c *                      gurson elements                                 *
c *                      updated: 6/12/2018 rhd                          *
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
c
      use segmental_curves, only : now_blk_relem

      implicit none
      include 'param_def'
c
c               parameter declarations
c
      integer :: step, iter, span, gpn, curve_type, felem, iout,
     &           hist_size
      double precision ::
     &  deps(mxvl,6), e(*), nu(*), sigma_o(*), f0(*), q1(*), q2(*),
     &  q3(*), stress_n(nstrs,*), stress_n1(nstrs,*),
     &  stress_n1_elas(mxvl,6,*), ym_n(*), nu_at_n(*),
     &  p_trial(*), q_trial(*), history(span,hist_size,*),
     &  history1(span,hist_size,*), yld_func(*), n_power(*), h_fixed(*)
c
      logical :: nucleation(*), nonlinear_flag(*), null_point(*),
     &           process_block, segmental
c
c               locals
c
      integer :: i, j, iword(2), state
      double precision ::
     &  a, b, c, d, shear_mod, term1, term2, term3,
     &  dword, h, plastic_strain,
     &  flow_stress(mxvl), hnow, ddumy, e_n, nu_n, g_n,
     &  een1, een2, een3, een4, een5, een6, trial_eps(6)
      double precision, external :: mm03sc
      equivalence ( dword, iword )
      logical ::  debug, gurson, previously_linear, now_nonlinear,
     &            active_point(mxvl)
c
      double precision, parameter ::
     &    root22 = dsqrt(2.d0)/2.d0, third = 1.d0/3.d0, zero = 0.d0,
     &    one = 1.d0, two = 2.d0, six = 6.d0, onehalf = 1.5d0,
     &    toler = 1.0d-10, half = 0.5d0, chk_yld_eps = 0.002d0
c
c               initialize history on step 1, time = 0
c
      debug = .false.
      if( step .eq. 1 ) call mm03p_a
c
c               check for null gauss points. set last n1 values zero
c               to prevent susequent uninitialzed variable computations
c               should stress update fail (eg, Gurson does not converge)
c
      if( debug ) write (iout,*) ' null point calcs'
!DIR$ VECTOR ALIGNED
      do i = 1, span
        null_point(i) = abs( deps(i,1) ) + abs( deps(i,2) ) +
     &                  abs( deps(i,3) ) + abs( deps(i,4) ) +
     &                  abs( deps(i,5) ) + abs( deps(i,6) )
     &                  .le. toler * chk_yld_eps
       active_point(i) = .not. null_point(i)
!DIR$ VECTOR ALIGNED
       stress_n1(7:nstrs,i) = zero
      end do
c
      process_block = any( active_point(1:span) )
      if( .not. process_block ) then !  all points are null
        state    = -1
        iword(1) = state
!DIR$ VECTOR ALIGNED
        do i = 1, span
         stress_n1(:,i)          = stress_n(:,i)
         stress_n1_elas(i,:,gpn) = zero
         history1(i,1:5,gpn)     = history(i,1:5,gpn)
         history1(i,6,gpn)       = dword
         history1(i,7:9,gpn)     = zero
!DIR$ VECTOR ALIGNED
         history1(i,10:15,gpn)   = history(i,10:15,gpn) ! no eps elastic change
       end do
       if ( debug ) write(iout,9000)
       return
      end if
c
c               we must process this block of elements.
c
      if( debug ) write(iout,9010) step, iter, span, gpn
c
c
c               1.  trial elastic stress at n+1. set actual
c                   stress at n+1 as trial stress.
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
c
c                      strain for trial stress computation at n+1
c                      is elastic strain at n + increment of
c                      mechaical strain over step over step.
c
!DIR$ VECTOR ALIGNED
        trial_eps(1:6) = history(i,10:15,gpn) + deps(i,1:6)
c
c                      trial elastic stress using strain at n+1 and
c                      elastic constants at n+1
c
        c = e(i) / ( ( one + nu(i) ) * ( one - two * nu(i) ) )
        a = c * ( one - nu(i) )
        b = c * nu(i)
        shear_mod = e(i) / ( two*(one+nu(i)))
c
        stress_n1(1,i) = a*trial_eps(1) + b*(trial_eps(2)+trial_eps(3))
        stress_n1(2,i) = a*trial_eps(2) + b*(trial_eps(1)+trial_eps(3))
        stress_n1(3,i) = a*trial_eps(3) + b*(trial_eps(1)+trial_eps(2))
        stress_n1(4,i) = shear_mod * trial_eps(4)
        stress_n1(5,i) = shear_mod * trial_eps(5)
        stress_n1(6,i) = shear_mod * trial_eps(6)
c
!DIR$ VECTOR ALIGNED
        stress_n1_elas(i,1:6,gpn) = stress_n1(1:6,i)
c
      end do ! over span

c
      if( debug ) then
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
!DIR$ VECTOR ALIGNED
      do i = 1, span
        p_trial(i) = - ( stress_n1(1,i) + stress_n1(2,i) +
     &                   stress_n1(3,i) ) * third
        a       = stress_n1(1,i) - stress_n1(3,i)
        b       = stress_n1(1,i) - stress_n1(2,i)
        c       = stress_n1(2,i) - stress_n1(3,i)
        d       = stress_n1(4,i)*stress_n1(4,i) +
     &               stress_n1(5,i)*stress_n1(5,i) +
     &               stress_n1(6,i)*stress_n1(6,i)
        q_trial(i) = root22 * sqrt( a*a + b*b +
     &                              c*c + six*d )
       end do
c
      if( debug ) then
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
!DIR$ VECTOR ALIGNED
      flow_stress(1:span) =  history(1:span,2,gpn)
c
      if ( segmental .and. curve_type .eq. 1 ) then
!DIR$ VECTOR ALIGNED
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
!DIR$ VECTOR ALIGNED
      do i = 1, span
        gurson = f0(i) .gt. zero .or. nucleation(i)
        if( gurson ) then
          term1 = q_trial(i) / flow_stress(i)
          term2 = two * q1(i) * history(i,5,gpn) *
     &               cosh( onehalf*q2(i)*p_trial(i)/
     &               flow_stress(i) )
          term3 = one + q3(i) * history(i,5,gpn) * history(i,5,gpn)
          yld_func(i) = term1*term1 + term2 - term3
        else
          yld_func(i) = q_trial(i) - flow_stress(i)
       end if
      end do
c
      if( debug ) then
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
!DIR$ VECTOR ALIGNED
      do i = 1, span
        previously_linear = history(i,1,gpn) .eq. zero
        now_nonlinear     = (history(i,1,gpn) .gt. zero) .or.
     &                      (yld_func(i) .gt. zero)
        nonlinear_flag(i) = now_nonlinear
        process_block = process_block .or. nonlinear_flag(i)
      end do
c
      if( debug ) then
         write(iout,9060)
         do i = 1, span
          write(iout,9062) i, nonlinear_flag(i)
         end do
         write(iout,9064) process_block
      end if
c
c              6.  for elements that are still linear, compute
c                  updated strain energy density. update accumulative
c                  elastic strains
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
        if( .not. nonlinear_flag(i) ) then
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
!DIR$ VECTOR ALIGNED
          history1(i,10:15,gpn) = history(i,10:15,gpn) + deps(i,1:6)
        end if
      end do
c
      if( debug ) then
         write(iout,9070)
         do i = 1, span
          if( .not.  nonlinear_flag(i) ) write(iout,9072) i,
     &            stress_n1(7,i)
         end do
      end if
c
c              7.  set element histories at n+1 to those at n
c                  those for nonlinear elements will be modified.
c                  ** except elastic strains **
c
      do j = 1, 9
!DIR$ VECTOR ALIGNED
        history1(1:span,j,gpn)  = history(1:span,j,gpn)
      end do

c
      return
c
 9000 format(/,'>> all points in element block are null')
 9010 format(/,'>> model 3 pre-processor for step, iter, span, gpn: ',
     &  4i6)
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
      contains
c     ========
      subroutine mm03p_a  ! setup step 1
      implicit none
c
      if( .not. segmental ) then
!DIR$ VECTOR ALIGNED
           do i = 1, span
              if( n_power(i) .gt. zero ) then
                 h = e(i) * (0.9999d0*e(i)) / ( e(i)-0.9999d0*e(i))
              else
                 h  = h_fixed(i)
              end if
              history(i,4,gpn)  = h
           end do
      endif
c
      iword(1) = -1
      iword(2) = 0
!DIR$ VECTOR ALIGNED
      do i = 1, span
          history(i,1,gpn)  = zero        ! equiv plastic strain
          history(i,2,gpn)  = sigma_o(i)  ! current static yld stress
          history(i,3,gpn)  = zero        ! current p_new  GT model
c                   4                     ! current H-prime
          history(i,5,gpn)  = f0(i)       ! current porosity GT model
          history(i,6,gpn)  = dword       ! current rate dependent yld stress
          history(i,7,gpn)  = zero        ! dep GT model
          history(i,8,gpn)  = zero        ! deq GT model
          history(i,9,gpn)  = zero        ! q   GT model
!DIR$ VECTOR ALIGNED
          history(i,10:15,gpn) = zero     ! accumulated elastic strains
      end do
c
c               compute initial elastic strain values using user-defined
c               initial stresses pre-loaded into stress_n.
c               these are zero unless model has user-initial stresses
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
        e_n  = ym_n(i)  ! at time 0 of simulation
        if( e_n .le. zero ) cycle ! sanity for bad user-input
        nu_n = nu_at_n(i)
        g_n  = e_n/two/(one+nu_n)
        een1 = (stress_n(1,i)-nu_n*(stress_n(2,i)+stress_n(3,i)))/e_n
        een2 = (stress_n(2,i)-nu_n*(stress_n(1,i)+stress_n(3,i)))/e_n
        een3 = (stress_n(3,i)-nu_n*(stress_n(1,i)+stress_n(2,i)))/e_n
        een4 = stress_n(4,i) / g_n
        een5 = stress_n(5,i) / g_n
        een6 = stress_n(6,i) / g_n
        history(i,10,gpn) = een1
        history(i,11,gpn) = een2
        history(i,12,gpn) = een3
        history(i,13,gpn) = een4
        history(i,14,gpn) = een5
        history(i,15,gpn) = een6
      end do
c
      return
c
      end subroutine mm03p_a
      end subroutine mm03p
c *******************************************************************
c *                                                                 *
c *      material model # 03 routine -- mm03f                       *
c *                                                                 *
c *******************************************************************
c
      function mm03f( stress, sbar, f, q1, q2, q3 ) result( value )
c
c              compute the value of the gurson yield
c              function for the state of stress tau.
c
      implicit none
c
      double precision :: stress(*), sbar, f, q1, q2, q3, value
c
      double precision :: p, a, b, c, d, q
      double precision, external :: mm03g
      double precision, parameter :: root22 = dsqrt(2.d0)/2.d0,
     &                               third = 1.d0/3.d0, six = 6.d0
c
      p  = - ( stress(1) + stress(2) + stress(3) ) * third
      a  = stress(1) - stress(3)
      b  = stress(1) - stress(2)
      c  = stress(2) - stress(3)
      d  = stress(4)*stress(4) + stress(5)*stress(5) +
     &     stress(6)*stress(6)
      q  = root22 * sqrt( a*a + b*b + c*c + six*d )
c
      value = mm03g( q, sbar, f, p, q1, q2, q3 )
c
      return
      end function mm03f
c *******************************************************************
c *                                                                 *
c *      material model # 03 routine -- mm03g                       *
c *                                                                 *
c *******************************************************************
c
      function mm03g( q, sbar, f, p, q1, q2, q3 ) result( value )
c
c
c              compute the value of the Gurson yield
c              given the scaler parameters
c
      implicit none
c
      double precision :: q, sbar, f, p, q1, q2, q3, value
c
      double precision :: term1, term2, term3
      double precision, parameter :: one = 1.d0, two = 2.d0,
     &                               onehalf = 1.5d0
c
      term1  = (q/sbar)*(q/sbar)
      term2  = two * q1 * f * cosh(onehalf*q2*p/sbar)
      term3  = one + q3 * f * f
      value  = term1 + term2 - term3
c
      return
      end function mm03g
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
      implicit double precision (a-h,o-z)
      double precision
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
      implicit double precision (a-h,o-z)
      double precision
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
      implicit none
c
c                 parameters
c
      integer :: iout
      double precision :: sbar_new, h_new, ebarp_new, ebarp,
     &                    sigma_o, h_constant, e, n_power,
     &                    m_power, dt, eps_ref, mpoweri
      logical :: consth, debug, rate_depen, segmental, power_law,
     &           first_time
c
c                locals
c
      double precision :: term0, term1, term2, sbar_inviscid,
     &                    h_inviscid, d_ebarp, rate_multiplier
      double precision, external :: mm03is, mm03sc
      double precision, parameter :: one = 1.d0
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
     &              result( value )
c
      implicit none
c
c              parameters
c
      integer :: iout
      double precision :: eps_pls, sigma_o, e, power, hprime, value
c
c              locals
c
      integer :: iterno
      integer, parameter :: max_itr = 20
      logical :: converge
      logical, parameter :: local_debug = .false.
      double precision :: sig_o, eps_o, poweri, eps_a, eps_b, h, h2,
     &                    h3, sig_a, sig_b, epower, pm1, et_b, strain,
     &                    stress, stress_new, fprime, ta, tb, n1, n2,
     &                    n3, n4, pn1, pn2, pn3, pn4, resid, deps,
     &                    et, strain_new
      double precision, parameter ::
     &   zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.d0,
     &   beta = 0.95d0, alpha = 1.1d0, toler = 0.000001d0
c
c              given the total plastic strain, return the stress,
c              and new plastic modulus. an iterative solution is
c              required to find the stress. a newton iteration
c              scheme converges very fast.
c
      if( eps_pls <= zero ) then
       value = sigma_o
       hprime = e * (0.9999d0*e) / ( e-0.9999d0*e)
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
      if( strain >= eps_b ) then
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
      converge =  abs( strain_new - strain ) .lt. toler * strain_new
     &                              .and.
     &            abs( stress_new - stress ) .lt. toler * stress_new
      if( .not. converge ) then
          iterno = iterno + 1
          if( iterno > max_itr ) then
            value = sig_o
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
      value =  stress_new
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
      implicit none
c
      double precision :: a_nuc, a_prime, ebarp, nuc_s_n, nuc_e_n,
     &                    nuc_f_n
      double precision :: term1, term2
      double precision, parameter :: root_2_pi = 2.50663d0,
     &                               half = 0.5d0
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
      implicit double precision (a-h,o-z)
      double precision
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
     &          result( value )
c
c                      parameter declarations
c
      use segmental_curves, only :  now_blk_relem, curve_plseps_values,
     &                              seg_curve_table, active_curve_set,
     &                              sigma_curves, num_seg_points,
     &                              max_seg_points, seg_curves_type,
     &                              curve_rates, sigma_inter_table
      implicit none
c
c                      parameter declarations
c
      integer :: caseh
      double precision :: ebarp, hprime, deplas, dtime, value
c
c                      local declarations
c
      integer :: i, first_curve, curve_set_type, pt_high, pt_low,
     &           num_curves_in_set, numpts
      double precision :: stress_values(max_seg_points),
     &   curve_high(max_seg_points), curve_low(max_seg_points),
     &   eps_pls_dot, sig_high, sig_low, eps_high, eps_low,
     &   rate_high, rate_low, hprime_high, hprime_low, dep,
     &   stress_low, stress_high, derivative
      double precision, external :: mm03lint
      double precision, parameter :: zero = 0d0, one = 1.d0,
     &           two = 2.d0, deriv_factor = 0.02d0
c
      first_curve       = seg_curve_table(2,active_curve_set)
      curve_set_type    = seg_curves_type(first_curve)
      num_curves_in_set = seg_curve_table(1,active_curve_set)
      numpts            = num_seg_points(first_curve)
c
c               check for input point before first point on curve
c               or beyond last point on curve.
c
      if( curve_set_type .eq. 2 ) then
        eps_pls_dot = zero
        if( dtime .gt. zero ) eps_pls_dot =  deplas / dtime
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
      if( ebarp .lt. curve_plseps_values(1) ) then
        value = stress_values(1)
        hprime = 1.0e10
        return
      end if
c
c                 check for input point beyond right end of curve
c
      if( ebarp .ge. curve_plseps_values(numpts) ) then
        value = stress_values(numpts)
        hprime = zero
        return
      end if
c
c            point is actually on the curve between points i and i-1.
c            interpolate linearly to get stress. set plastic modulus.
c
      do i = 2, numpts
        if( ebarp .lt. curve_plseps_values(i) ) then
          sig_high = stress_values(i)
          sig_low  = stress_values(i-1)
          eps_high = curve_plseps_values(i)
          eps_low  = curve_plseps_values(i-1)
          hprime   = (sig_high - sig_low) / ( eps_high - eps_low )
          value    = sig_low + hprime * ( ebarp - eps_low )
          pt_high  = i
          pt_low   = i-1
          go to 1000
        end if
      end do
c
c            if we get here things are very bad. stop job
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
      if( curve_set_type .ne. 2 ) return
      if( caseh .eq. 1 ) return
      if( deplas .le. zero ) return
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
c     *                   last modified: 6/13/2018 rhd               *
c     *                                                              *
c     *    execute linear interpolation on a tabular function        *
c     *    of a single variable where the x values are sorted in     *
c     *    increasing value but are not req'd to be uniformly spaced *
c     *                                                              *
c     ****************************************************************
c
      function mm03lint( xvalue, n, x, y ) result( value )
      implicit none
c
      integer :: n
      double precision :: xvalue, x(n), y(n), value
c
      integer :: i, point
      double precision :: x1, x2, y1, y2
c
      if( xvalue .le. x(1) ) then
        value = y(1)
        return
      end if
c
      if( xvalue .ge. x(n) ) then
        value = y(n)
        return
      end if
c
      do point = 2, n
        if( xvalue .gt. x(point) ) cycle
        x1 = x(point-1)
        x2 = x(point)
        y1 = y(point-1)
        y2 = y(point)
        value = y1 + (xvalue-x1)*(y2-y1)/(x2-x1)
        return
      end do
c
      write(*,*) '>> FATAL ERROR: mm03lint'
      write(*,*) '                job aborted'
      write(*,*) '>> xvalue: ',xvalue
      do i = 1, n
        write(*,*) '  i,x,y: ',i,x(i),y(i)
      end do
      call die_abort
c
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm03_set_sizes                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 6/13/2018 rhd               *
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
      info_vector(1) = 15
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
      use global_data ! old common.main
c
c                       access some global data structures
c
      use elem_block_data, only: history_blocks, history_blk_list
      use main_data, only: elems_to_blocks
c
      implicit integer (a-z)
c
c                       parameters
c
      integer :: nrow_states, itype, num_states
      double precision :: elem_states_output(nrow_states,*)
c
c                       locals
c
      double precision,
     & allocatable :: history_dump(:,:,:), one_elem_states(:)
      integer :: relem, elnum, hist_size, blockno
      logical :: do_a_block, local_debug
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
      do_a_block = .true.
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
      double precision ::
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




c *******************************************************************
c *                                                                 *
c *        material model # 3 --  stress update for gurson model    *
c *                               and mises model using backward    *
c *                               euler integration                 *
c *                                                                 *
c *******************************************************************
c
c
c      parameter definitions:
c      ======================
c
c                  mode/type
c
c   ae               1,2         -- young's modulus
c
c   anu              1,2         -- poisson's ratio
c
c   af0              1,2         -- initial porosity for gurson material
c
c   aeps_ref         1,2         -- reference strain rate for power-law
c                                   visco-plastic model
c
c   asigyld          1,2         -- uniaxial yield stress (inviscid)
c
c   am_power         1,2         -- exponent for power-law visco-plastic
c
c   an_power         1,2         -- exponent for the power law part of
c                                   the inviscid, uniaxial stress-strain
c                                   curve
c
c   ah_fixed         1,2         -- inviscid plastic modulus for a material
c                                   with constant strain hardening,
c                                   initial plastic modulus for segemental
c                                   uniaxial response
c
c   aq1, aq2,        1,2         -- q1, q2, q3 constants used in gurson
c   aq3                             model
c
c   anucleation      1,3         -- .true. if void nucleation is to be
c                                   included in the stress update
c
c   anuc_s_n         1,2         -- nucleation constants for gurson
c   anuc_e_n,                       model (s_n, e_n, f_n)
c   anuc_f_n
c
c   stress_n         1,2         -- 7x1 vector of stresses at start
c                                   of load step (i.e. at 'n')
c                                   ordering: x, y, z, xy, yz, xz,
c                                   work density
c
c   stress_n1        1,2         -- updated stress state at n+1 (7x1)
c                                   loaded with trial elastic stress
c                                   by model pre-processor
c
c   deps             1,2         -- mxvl x 6 array of total strain increment
c                                   for the load step (deps advanes the
c                                   solution from n->n+1). ordering is above.
c
c   history          1,2         -- span x hist_size x #gpns array of
c                                   history data for state 'n' (start of
c                                   step). the history vector
c                                   must be made available to the
c                                   consistent tangent modulus routine.
c
c   history1         3,2         -- span x hist_size x #gpns array of
c                                   history data for state 'n+1'
c                                   (end of step). prior to calling,
c                                   contains history at state 'n'.
c                                   update to n+1 here.
c
c   f1               1,2         -- yield function evaluated for trial
c                                   elastic stress state at n+1
c
c   aq_trial         1,2         -- mises equivalent stress for the trial
c                                   elastic stress state at n+1
c
c   ap_trial         1,2         -- negative of mean stress for the trial
c                                   elastic stress state at n+1
c
c   mxvl             1,1         -- dimensioned number of rows for (most)
c                                   element block data structures
c
c   hist_size        1,1         -- number of history words per gauss point
c                                   required by this material model
c
c
c        ** these parameters are passed through the args derived **
c        ** type to reduce overhead in calling routine           **
c
c
c                  mode/type
c   step             1,1         -- load step number in the analysis
c
c   aiter            1,1         -- equilibrium iteration for the load
c                                   step.
c
c   element          1,1         -- element number being processed
c
c   relem            1,1         -- relative element number of the
c                                   block being processed
c
c   gpn              1,1         -- gauss point (or strain point) being
c                                   processed
c
c   atime_incr       1,2         -- time increment for load step. used
c                                   for viscoplastic computations
c
c   jout             1,1         -- output device number for error messages.
c
c   allow_cut        1,3         -- specifies if material model will be allowed
c                                   to request an adaptive step size cut due to
c                                   reversed plasticity
c
c   cut_step         2,3         -- set true by model if a reduced global
c                                   step size is needed in material model
c                                   computations
c   signal_flag      1,3         -- .true. if material model is allowed
c                                   to write messages about state changes
c                                   in material at the gauss point
c
c   asegmental       1,3         -- .true. if stress-strain curve is
c                                   represented with segmental model
c
c   apower_law       1,3         -- .true. if stress-strain curve is
c                                   represented with linear+power-law
c                                   model
c
c
c   routine names:
c   =============
c
c         all support routines are named with the prefix mm03
c
c   mode keys:
c   =========
c         1 (input)
c         2 (output)
c         3 (input by calling routine, updated by this routine and
c            returned)
c
c   type keys:
c   =========
c         1 (integer)
c         2 (floating point; real or double depending on computer)
c         3 (logical)
c
c   stress-strain curve models:
c   ===========================
c
c         the inviscid stress-strain curve can be modeled as linear+
c         constant hardening, linear+power-law hardening or
c         piecewise linear (segmental).
c         for linear+constant hardening:
c             define e, nu, sigyld, h_fixed. set n_power=0.0
c         for linear+power-law hardening:
c             define e, nu, sigyld, n_power. set h_fixed=0.0
c         for segmental curve:
c             define e, nu, sigyld, points on curve. set h_fixed=0.0
c
c         to include power-law viscoplasticity for matrix material:
c            define inviscid stress-strain curve above. define
c            eps_ref, m_power and time_incr. set these = 0.0 for
c            an inviscid response.
c
c         to provide general viscoplastic response, use a set of stress
c         vs. plastic strain curves specified at various plastic strain
c         rate.
c
c         temperature dependence of the uniaxial (matrix) response is
c         modeled through a user specified set of stress vs. plastic
c         strain curves at various temperatures. the calling routine
c         has interpolated the curves to build teh one for this gauss
c         point.e
c
c   double-single precision:
c   =======================
c
c         depend on single/double precision implementations.
c
c
c   large-strain issues:
c   ====================
c
c         this routine is unaware of finite-strain vs. small-strain
c         issues. we assume that "rotation neutral" stresses and strain
c         increments are passed as arguments. this can be accomplished
c         in a variety of ways including:
c           a) a jaumann rate approach with the hughes-winget technique
c              to handle the rotation
c           b) using the "unrotated" reference configuration of
c              dienes, halquist, taylor & flanagan, dodds & healy...
c
      subroutine mm03(
     &  args, ae, anu, af0, aeps_ref, asigyld, am_power, an_power,
     &  ah_fixed, aq1, aq2, aq3, anucleation, anuc_s_n, anuc_e_n,
     &  anuc_f_n, stress_n, stress_n1, stress_n1_elas, deps, history,
     &  history1, f1, ap_trial, aq_trial, mxvl, hist_size,
     &  asig_cur_min_val, span )
c
      implicit double precision (a-h,o-z)
c
c                   parameter declarations
c
      integer hist_size, span
      dimension
     &    deps(mxvl,*), stress_n(*), stress_n1(*),
     &    history(span,hist_size,*), history1(span,hist_size,*),
     &    stress_n1_elas(mxvl,6,*)
      logical anucleation
c
      type :: arguments
        integer :: iter, abs_element, relem, ipoint, iout
        logical :: allow_cut, segmental, power_law,
     &             rate_depend_segmental, signal_flag, cut_step
        double precision :: dtime, step_scale_fact
      end type
      type (arguments) ::args
c
c              the following variables are local to this
c              routine but are accessible by all the subsequent
c              routines within the "contains" section. This is
c              a fortran 90 feature to eliminate a
c              common area. none of these are initialized when
c              this routine is called. these variables go on the
c              stack and thus do not create problems with
c              execution of this code within a threaded
c              parallel construct.
c
c        ---- start of the shared locals with contains -------
c
      double precision
     &  e, nu, f0, eps_ref, sigma_o, m_power, n_power, h_fixed,
     &  q1, q2, q3, nuc_s_n, nuc_e_n, nuc_f_n,
     &  stime, ftime, sbar, dmeps(6), ebarp, ebarp_new,
     &  dtauel(6), f2, smel, std(6), deps_sub(6), stress_start(6),
     &  p_new, q_new, shear_mod, dep_last_sub, deq_last_sub,
     &  et, h, g, f, p_trial, q_trial, deq, dep,
     &  h_new, f_new, time_incr, newyld,
     &  senerg, threeg, mpoweri, dt_eref, h_constant,
     &  dsenrg, oldsig(6), dep_old, deq_old,
     &  bulk_mod, sbar_new, dep_n, deq_n, equiv_eps,
     &  eps_o, h_inviscid, inviscid_stress, dword, zero, half,
     &  one, two, three, third, six, hprime_init,
     &  dsig(6), deps_plas(6), factor1, factor2, deps_plas_bar,
     &  sig_cur_min_val

c
      logical  debug, plastic_loading, nucleation, gurson,
     &         unload, signal, uloadp, rate_depen,
     &         consth, lword(2), null_point,
     &         power_law, segmental, first_time
      logical  allow_cut
c
      integer  element, gpn, relem
      integer  state, dbele, dbptno, iword(2), elenum, ptno,
     &         gt_converge, vm_converge, iout, iterno
c
      real  rword(2)
      equivalence  ( rword,iword,lword,dword )
      data zero, half, one, two, three, third, six
     &  / 0.0, 0.5, 1.0, 2.0, 3.0, 0.3333333333, 6.0 /
c
c        ---- end of the shared locals with contains -------
c
      double precision
     &  root23, ptone, root22, eps_tol, root2
      data root23, ptone, root22, root2 / 0.4714, 0.1, 0.70710678,
     &                                    1.41421356   /
      data  eps_tol / 100.0 /
c
c                   put parameters into shared local variables
c                   for use by the contains routines
c
      relem      = args%relem
      allow_cut  = args%allow_cut
      iout       = args%iout
      e          = ae
      nu         = anu
      f0         = af0
      eps_ref    = aeps_ref
      sigma_o    = asigyld
      m_power    = am_power
      n_power    = an_power
      h_fixed    = ah_fixed
      q1         = aq1
      q2         = aq2
      q3         = aq3
      nucleation = anucleation
      nuc_s_n    = anuc_s_n
      nuc_e_n    = anuc_e_n
      nuc_f_n    = anuc_f_n
      time_incr  = args%dtime
      elenum     = args%abs_element
      element    = args%abs_element
      ptno       = args%ipoint
      gpn        = args%ipoint
      gurson     = f0 .gt. zero .or. nucleation
      p_trial    = ap_trial
      q_trial    = aq_trial
      signal     = args%signal_flag
      iterno     = args%iter
      segmental  = args%segmental
      power_law  = args%power_law
      step_scale = args%step_scale_fact
      first_time = .true.
      if ( segmental       ) sig_cur_min_val = asig_cur_min_val
      if ( .not. segmental ) sig_cur_min_val = sigma_o
c
      debug = .false.
      if ( debug ) then
        write(iout,9000) element, gpn, e, nu, f0, eps_ref,
     &                   sigma_o, m_power, n_power, h_fixed,
     &                   q1, q2, q3, nucleation, nuc_s_n,
     &                   nuc_e_n, nuc_f_n, eps_ref,
     &                   time_incr,
     &                   (deps(relem,k),k=1,6),
     &                   (stress_n(k),k=1,7), gurson, f1, p_trial,
     &                   q_trial
      end if
c
c              we have a constant hardening modulus if the
c              hardening exponent is zero. we have a power-law rate
c              dependent response if m_power > 0 and time step > 0.
c              the segmental response can be rate dependent through
c              a set of curves or through a single curve + the power-law
c              rate dependence
c
      if ( segmental ) then
         consth       = .false.
         eps_o        = sigma_o / e
         shear_mod    = half * e /(one+nu)
         bulk_mod     = e * third / ( one - two*nu )
         threeg       = three * shear_mod
         hprime_init  = h_fixed
         rate_depen   = time_incr .gt. zero .and. m_power .gt. zero
         if ( rate_depen ) mpoweri = one / m_power
      else
         consth       = n_power .eq. zero
         eps_o        = sigma_o / e
         shear_mod    = half * e /(one+nu)
         bulk_mod     = e * third / ( one - two*nu )
         threeg       = three * shear_mod
         hprime_init  = zero
         rate_depen   = time_incr .gt. zero .and. m_power .gt. zero
         if ( n_power .eq. zero ) h_constant = h_fixed
         if ( rate_depen ) then
              mpoweri = one / m_power
c              dt_eref = time_incr / eps_ref
         end if
      end if
c
c              set state variable values to those at start of the
c              load step. starting values for dep, deq newton
c              loop are later set to a fraction of those converged
c              at end of last step.
c
      ebarp      = history(relem,1,gpn)
      sbar       = history(relem,2,gpn)
      h          = history(relem,4,gpn)
      f          = history(relem,5,gpn)
      dword      = history(relem,6,gpn)
      state      = iword(1)
      dep_n      = history(relem,7,gpn)
      deq_n      = history(relem,8,gpn)
      senerg     = stress_n(7)
      dtauel(1)  = stress_n1_elas(relem,1,gpn) - stress_n(1)
      dtauel(2)  = stress_n1_elas(relem,2,gpn) - stress_n(2)
      dtauel(3)  = stress_n1_elas(relem,3,gpn) - stress_n(3)
      dtauel(4)  = stress_n1_elas(relem,4,gpn) - stress_n(4)
      dtauel(5)  = stress_n1_elas(relem,5,gpn) - stress_n(5)
      dtauel(6)  = stress_n1_elas(relem,6,gpn) - stress_n(6)
c
c
c              trial stress state was computed by model
c              pre-processor.
c
      if ( debug ) write(iout,9012) sbar, ebarp, state, e, nu,
     &                              f, h, dep_n, deq_n
c
c               yield function value (f1) computed by model
c               pre-processor. determine if point is still linear,
c               plastically loading or elastically unloading.
c               use value of yield function at trial elastic state.
c               state =  1 (plastic at start of step)
c                     = -1 (elastic at start of step)
c
      if ( state .eq. -1 ) then
         plastic_loading = f1 .gt. zero
      else
         plastic_loading = f1 .ge. zero
      end if
c
      if ( plastic_loading ) go to 500
c
c              material point is still linear or is
c              elastically unloading.
c
      dsenrg = deps(relem,1) * ( stress_n(1)+stress_n1(1) )  +
     &         deps(relem,2) * ( stress_n(2)+stress_n1(2) )  +
     &         deps(relem,3) * ( stress_n(3)+stress_n1(3) )  +
     &         deps(relem,4) * ( stress_n(4)+stress_n1(4) )  +
     &         deps(relem,5) * ( stress_n(5)+stress_n1(5) )  +
     &         deps(relem,6) * ( stress_n(6)+stress_n1(6) )
      senerg       = senerg + dsenrg * half
      stress_n1(7) = senerg
      stress_n1(8) = stress_n(8)
      stress_n1(9) = stress_n(9)
      if ( debug ) write(iout,9016) f1, senerg
c
c              update state and history variables. leave
c              porosity unchanged.
c
      if ( state .eq. 1 ) then
        if ( debug ) write(iout,9036) sbar, f1, senerg
        if ( signal .and.
     &       ( abs(f1) .gt. 0.01 * sigma_o)  ) then
          write(iout,9030) element, gpn, f1
        end if
      end if
      state    = -1
      iword(1) = state
      history1(relem,6,gpn)  = dword
      history1(relem,7,gpn)  = zero
      history1(relem,8,gpn)  = zero
      history1(relem,9,gpn)  = zero
      history1(relem,10:15,gpn) = history(relem,10:15,gpn) +
     &                            deps(relem,1:6) ! elas eps
      return
c
c              trial elastic stresses are outside the yield surface.
c              if just yielded from an elastic state, notify user.
c              for the gurson material, if strain increment is
c              considered really large, notify user and set
c              step cut back flag. numerical testing
c              indicates that the integration procedure deteriorates
c              for strain increments > 20 * yield strain. the updated
c              stresses contain enough errors to possibly cause spurious
c              unloadings in follow-on load steps.
c
c
  500 continue
      if ( gurson ) then
        equiv_eps = (deps(relem,1)-deps(relem,2))**2 +
     &              (deps(relem,2)-deps(relem,3))**2 +
     &              (deps(relem,1)-deps(relem,3))**2 +
     &              half * ( deps(relem,4)**2 +
     &              deps(relem,5)**2 + deps(relem,6)**2 )
        equiv_eps = root23 * sqrt ( equiv_eps )
        if ( equiv_eps .gt. eps_tol * eps_o ) then
         if ( signal  ) then
           write(iout,9021) element, gpn, equiv_eps/eps_o
         end if
         if ( allow_cut ) then
               args%cut_step  = .true.
               write(iout,9040) element, gpn, equiv_eps/eps_o
               return
            else
               write(iout,9042)
               call abort_job
         end if
        end if
      end if
c
      if ( state .ne. 1 .and. signal  )
     &         write(iout,9017) element, gpn, f1
c
c              given the material state at start of step, update
c              stresses using implicit, return mapping. stress_n1
c              contains updated, unrotated cauchy stress at end of step
c              on return. we branch to separate routines
c              to handle (1) a real gurson material with voids or
c              (2) a mises material. the mises update procedure
c              is much simpler. for mises we use a single-step
c              procedure. for gurson, we use subincrements
c              if the strain increment is large. for a viscoplastic
c              response we use 1 subincrement since I can't get
c              the algorithm to work with subincrements.
c
      if ( gurson ) then
c
        if ( ebarp .gt. zero ) then
          nsubin    = max( idint(one+ptone*equiv_eps/eps_o), 1 )
        else
          nsubin = 1
        end if
        if ( rate_depen ) nsubin = 1
        scale     = one / dble(nsubin)
        total_dep = zero
        total_deq = zero
        dep_n     = dep_n * scale * step_scale
        deq_n     = deq_n * scale * step_scale
        do i = 1, 6
          dtauel(i)       = dtauel(i) * scale
          deps_sub(i)     = deps(relem,i) * scale
          stress_start(i) = stress_n(i)
        end do
c
        do isub = 1, nsubin
           stress_n1(1) = stress_start(1) + dtauel(1)
           stress_n1(2) = stress_start(2) + dtauel(2)
           stress_n1(3) = stress_start(3) + dtauel(3)
           stress_n1(4) = stress_start(4) + dtauel(4)
           stress_n1(5) = stress_start(5) + dtauel(5)
           stress_n1(6) = stress_start(6) + dtauel(6)
           p_trial = - ( stress_n1(1) + stress_n1(2) +
     &                   stress_n1(3) ) * third
           a       = stress_n1(1) - stress_n1(3)
           b       = stress_n1(1) - stress_n1(2)
           c       = stress_n1(2) - stress_n1(3)
           d       = stress_n1(4)*stress_n1(4) +
     &               stress_n1(5)*stress_n1(5) +
     &               stress_n1(6)*stress_n1(6)
           q_trial = root22 * sqrt( a*a + b*b + c*c + six*d )
           call mm03s( stress_n1, isub )
           if ( gt_converge .gt. 0 ) then
            if ( allow_cut ) then
               args%cut_step = .true.
               write(iout,9050) element, gpn, gt_converge
               return
            else
               write(iout,9060) element, gpn, gt_converge
               call abort_job
            end if
           end if
           dsenrg = deps_sub(1) * ( stress_start(1)+stress_n1(1) )  +
     &              deps_sub(2) * ( stress_start(2)+stress_n1(2) )  +
     &              deps_sub(3) * ( stress_start(3)+stress_n1(3) )  +
     &              deps_sub(4) * ( stress_start(4)+stress_n1(4) )  +
     &              deps_sub(5) * ( stress_start(5)+stress_n1(5) )  +
     &              deps_sub(6) * ( stress_start(6)+stress_n1(6) )
           senerg  = senerg + dsenrg * half
           ebarp     = ebarp_new
           sbar      = sbar_new
           h         = h_new
           f         = f_new
           total_dep = total_dep + dep
           total_deq = total_deq + deq
           stress_start(1) = stress_n1(1)
           stress_start(2) = stress_n1(2)
           stress_start(3) = stress_n1(3)
           stress_start(4) = stress_n1(4)
           stress_start(5) = stress_n1(5)
           stress_start(6) = stress_n1(6)
        end do
        deq = total_deq
        dep = total_dep
      end if
c
      if ( .not. gurson ) then
           call mm03q( stress_n1 )
           if ( vm_converge .eq. 1 .or. vm_converge .eq. 2 ) then
              if ( allow_cut ) then
                 args%cut_step = .true.
                 write(iout,9052) element, gpn, vm_converge
                 return
              else
                 write(iout,9062) element, gpn
                 call abort_job
              end if
           end if
           if ( vm_converge .eq. 3 ) write(iout,9064) element, gpn
           dsenrg = deps(relem,1) * ( stress_n(1)+stress_n1(1) )  +
     &         deps(relem,2) * ( stress_n(2)+stress_n1(2) )  +
     &         deps(relem,3) * ( stress_n(3)+stress_n1(3) )  +
     &         deps(relem,4) * ( stress_n(4)+stress_n1(4) )  +
     &         deps(relem,5) * ( stress_n(5)+stress_n1(5) )  +
     &         deps(relem,6) * ( stress_n(6)+stress_n1(6) )
           senerg = senerg + dsenrg * half
      end if
c
c              update history for strain point to reflect that
c              material is plastic. for the matrix we save
c              the plastic strain, equivalent stress, plastic
c              modulus, porosity, the state flag. we save the
c              porosity and plastic strain at start of step for
c              use by the consistent [cep] routine. the converged
c              values of dep and deq are saved for use in starting
c              the next iteration of this step. we save the pressure
c              and equivalent stress for updated stress state (needed
c              by consistent tangent routine). the plastic modulus
c              includes rate dependent modifications.
c
      stress_n1(7) = senerg
      state       = 1
      history1(relem,1,gpn)  = ebarp_new
      history1(relem,2,gpn)  = sbar_new
      history1(relem,3,gpn)  = p_new
      history1(relem,4,gpn)  = h_new
      history1(relem,5,gpn)  = f_new
      iword(1)               = state
      history1(relem,6,gpn)  = dword
      history1(relem,7,gpn)  = dep
      history1(relem,8,gpn)  = deq
      history1(relem,9,gpn)  = q_new
c
      dsig(1) = stress_n1(1) - stress_n(1)
      dsig(2) = stress_n1(2) - stress_n(2)
      dsig(3) = stress_n1(3) - stress_n(3)
      dsig(4) = stress_n1(4) - stress_n(4)
      dsig(5) = stress_n1(5) - stress_n(5)
      dsig(6) = stress_n1(6) - stress_n(6)
c
      deps_plas(1) = deps(relem,1) - (dsig(1) - nu*(dsig(2)+dsig(3)))/e
      deps_plas(2) = deps(relem,2) - (dsig(2) - nu*(dsig(1)+dsig(3)))/e
      deps_plas(3) = deps(relem,3) - (dsig(3) - nu*(dsig(1)+dsig(2)))/e
      deps_plas(4) = deps(relem,4) - dsig(4) / shear_mod
      deps_plas(5) = deps(relem,5) - dsig(5) / shear_mod
      deps_plas(6) = deps(relem,6) - dsig(6) / shear_mod
c
      history1(relem,10:15,gpn) = history(relem,10:15,gpn) +
     &                            deps(relem,1:6) - deps_plas(1:6)! elas eps
c
      stress_n1(8) = stress_n(8)  +  half * (
     &       deps_plas(1) * (stress_n1(1) + stress_n(1))
     &     + deps_plas(2) * (stress_n1(2) + stress_n(2))
     &     + deps_plas(3) * (stress_n1(3) + stress_n(3))
     &     + deps_plas(4) * (stress_n1(4) + stress_n(4))
     &     + deps_plas(5) * (stress_n1(5) + stress_n(5))
     &     + deps_plas(6) * (stress_n1(6) + stress_n(6)) )
c
      factor1 = ( deps_plas(1) - deps_plas(2) )**2  +
     &          ( deps_plas(2) - deps_plas(3) )**2  +
     &          ( deps_plas(1) - deps_plas(3) )**2
      factor2 = deps_plas(4)**2 +  deps_plas(5)**2 +
     &          deps_plas(6)**2
      deps_plas_bar =  (root2/three) * sqrt( factor1 +
     &                 (three/two)*factor2 )
      stress_n1(9) = stress_n(9) + deps_plas_bar
      if ( debug ) then
         write(iout,9070) deps_plas(1:6)
         write(iout,9072) deps(relem,1:6)
         write(iout,9074) dsig(1:6)
      end if
c
      return
c
c
 9000 format(/,'>> debug from mm03. elem, gpn : ',i8,i2,
     & /,    '    e, nu, f0, eps_ref : ',4f10.3,
     & /,    '    sigma_o, m_power, n_power, h_fixed : ',4f10.3,
     & /,    '    q1, q2, q3, nucleation: ',3f10.3, l10,
     & /,    '    nuc_s_n, nuc_e_n, nuc_f_n : ',3f10.3,
     & /,    '    eps_ref, time_incre : ',f10.3,e10.3,
     & /,    '    deps :',
     & /,10x,3e15.6,/,10x,3e15.6,
     & /,    '    stresses @ n :',
     & /,10x,3e15.6,/,10x,4e15.6
     & /,    '    gurson, f1, p_trial : ',l1,2f10.3,
     & /,    '    q_trial: ',f10.3 )
 9001 format('Updated plastic strain rates. ele, gpn, old, new:',
     & i6,i3,2f15.5)
 9012 format(//, 5x, 'history vector data:',
     &        /, 8x, 'current yield stress', f10.3,
     &        /, 8x, 'total plastic strain', f10.7,
     &        /, 8x, 'state variable      ', i5,
     &        /, 8x, 'youngs modulus      ', f10.1,
     &        /, 8x, 'poissons ratio      ', f10.3,
     &        /, 8x, 'void fraction       ', f10.3,
     &        /, 8x, 'plastic modulus     ', e10.2
     &        /, 8x, 'dep @ start of step ', f10.6,
     &        /, 8x, 'deq @ start of step ', f10.6 )
 9016 format( //, 5x, 40hpoint remains elastic.  yield function = ,
     &        f12.3 ,
     &  /,5x, 24hstrain energy density = ,f10.4/)
 9017 format(10x,i5,i3,' point yields.  f    = ',f12.3 )
 9020 format(10x,i5,i3,' reversed plastic yielding. dot, f =',2f12.3 )
 9021 format(10x,i5,i3,' large strain incr. / eps_o =',f12.3 )
 9022 format(16x,' excessive reversed plasticity in step' )
 9030 format(10x,i5,i3,' elastic unloading f = ',f12.3 )
 9036 format( //, 5x, 'point is unloading'
     &         /, 8x, 'current yield stress = ', f10.2,
     &         /, 8x, 'yield function       = ', f10.3 ,
     &         /, 8x, 'strain energy        = ',f12.9///)
 9040 format(
     &/,3x,'>> Warning: strain increment too large.',
     &/,3x,'            material model requesting step size reduction.',
     &/,3x,'            element, gauss point, deps / eps_o: ',
     &             i8,i3, f12.3 )
 9042 format(
     &/,3x,'>> Warning: strain increment too large.',
     &/,3x,'            load step reduction not enabled.,'
     &/,3x,'            analysis terminated.' )
 9044 format(
     &/,3x,'>> Warning: strain increment too large.',
     &/,3x,'            load step reduction not enabled for',
     &/,3x,'            loads computation due to imposed',
     &/,3x,'            displacements at start of step.',
     &/,3x,'            this may cause problems later...')
 9050 format(
     &/,3x,'>> Warning: iterations for gurson model failed.',
     &/,3x,'            material model requesting step size reduction.',
     &/,3x,'            element, gauss point, gt-flag: ',i8,i3,i3 )
 9052 format(
     &/,3x,'>> Warning: iterations for mises model failed.',
     &/,3x,'            material model requesting step size reduction.',
     &/,3x,'            element, gauss point, gt-flag: ',i8,i3,i3 )
 9060 format(
     &/,3x,'>> FATAL ERROR: iterations for gurson model failed.',
     &/,3x,'                did not converge in 30 iterations.',
     &/,3x,'                load step reduction not enabled.,'
     &/,3x,'                analysis terminated.',
     &/,3x,'                element, gauss point, gt-flag: ',i8,i3,i3 )
 9062 format(
     &/,3x,'>> FATAL ERROR: iterations for mises model failed.',
     &/,3x,'                did not converge.',
     &/,3x,'                load step reduction not enabled.,'
     &/,3x,'                analysis terminated.',
     &/,3x,'                element, gauss point, gt-flag: ',i8,i3,i3 )
 9064 format(
     &/,3x,'>> WARNING: during mises update, the increment of',
     &/,3x,'                plastic strain is negative.',
     &/,3x,'                element, gauss point ',i8,i3 )
c
 9070 format(
     & /,    '    deps plastic :',
     & /,10x,3e15.6,/,10x,3e15.6 )
c
 9072 format(
     & /,    '    deps :',
     & /,10x,3e15.6,/,10x,3e15.6 )
c
 9074 format(
     & /,    '    dsig :',
     & /,10x,3e15.6,/,10x,3e15.6 )

      contains

c ************************************************************************
c *                                                                      *
c *    routine  mm03s --    3-d radial return procedure Gurson model     *
c *                         (for a single subincrement of strain)        *
c *                                                                      *
c ************************************************************************
c
c
      subroutine mm03s( newsig, isubincr )
c
      implicit none
      double precision
     &  newsig(*), mm03f
      integer isubincr
c
c              This routine accesses the variables defined in the
c              calling routine thru the f-90 host by association
c              feature.
c
c                   local variables in this routine
c
      logical   converge
c
      double precision
     & toler, phi, gp, gq, gsbar, gf, d11, d12, d21, d22, r1,
     & r2, del_dep, del_deq, eps_compare, chk1, tol_yf
c
      integer max_iterations, loop
c
      data  toler, tol_yf, max_iterations /  0.00001, 0.001, 30 /
c
c            compute scalars required to correct
c            trial elastic stress vector to the yield surface:
c
c              dep       -- increment of macroscopic, plastic volume
c                           change
c              deq       -- increment of macroscopic, deviatoric
c                           plastic strain
c              f_new     -- updated void fraction
c              ebarp_new -- updated plastic strain in matrix
c              sbar_new  -- updated equivslent stress in matrix
c
c
c            A. set initial estimates for dep and deq. currently,
c               we use 0.0 for the first subincrement of strain. use
c               the previously converged values from previous
c               subincrement for a multi-subincrement update process.
c               these seem to be the most stable values, found after
c               much experimentation.
c
      gt_converge = 0
      if ( isubincr .gt. 1 ) then
        dep = dep_last_sub
        deq = deq_last_sub
      else
        dep = zero
        deq = zero
      end if
      dep_old = dep
      deq_old = deq
c
c            B. start the main iteration loop to update state variables.
c               we run at least two iterations to catch any small
c               problems with false convergence upon start up early in
c               the loading.
c
      do loop = 1, max_iterations
c
      p_new = p_trial + bulk_mod * dep
      q_new = q_trial - three * shear_mod * deq
      if ( debug ) write(iout,9100) loop, dep, deq, p_new, q_new
c
c            C. execute a local newton iteration to update
c               ebarp, sbar, and f -> ebarp_new, sbar_new, f_new,
c               h_new. call a fast routine first that sometimes
c               does not converge. if it fails we call a slower,
c               but very reliable routine.
c
      call mm03us( p_new, dep, q_new, deq, ebarp, f, sbar,
     &             consth, sigma_o, eps_o, h_constant, e,
     &             n_power, m_power, time_incr, eps_ref, rate_depen,
     &             mpoweri, nucleation, nuc_s_n, nuc_e_n,
     &             nuc_f_n, ebarp_new, sbar_new, f_new, h_new,
     &             iout, debug, segmental,
     &             power_law, converge, first_time )
      if ( .not. converge ) then
         call mm03ss( p_new, dep, q_new, deq, ebarp, f, sbar,
     &                consth, sigma_o, eps_o, h_constant, e,
     &                n_power, m_power, time_incr, eps_ref, rate_depen,
     &                mpoweri, nucleation, nuc_s_n, nuc_e_n,
     &                nuc_f_n, ebarp_new, sbar_new, f_new, h_new,
     &                elenum, ptno, iout, debug, segmental,
     &                power_law, converge, first_time )
      end if
      if ( .not. converge ) then
       gt_converge = 1
       return
      end if
c
c            D. evaluate the yield function and derivatives of
c               the yield function wrt the state variables. evaluate
c               the four (4) partial derivatives of R1 and
c               R2 that define th Jacobian for the 2 nonlinear
c               scalar equations.
c
      call mm03dr( phi, gp, gq, gsbar, gf, q_new, sbar_new,
     &             p_new, f_new, ebarp_new, f, ebarp, nucleation,
     &             shear_mod, bulk_mod, h_new, q1,
     &             q2, q3, nuc_s_n, nuc_e_n, nuc_f_n, dep, deq,
     &             d11, d12, d21, d22, debug, iout )
c
c            E. evaluate the residual functions R1 and R2 using
c               current state variables.
c
      r1 = dep * gq + deq * gp
      r2 = phi
c
c            G. solve the 2x2 set of linear equations to define
c               correction terms for dep and deq.
c
      del_dep = (d12*r2 - d22*r1 ) / (d11*d22 - d12*d21)
      del_deq = (d21*r1 - d11*r2) / (d11*d22 - d12*d21)
      dep     = dep + del_dep
      deq     = deq + del_deq
c
c            F. check for convergence. we used to force 2
c               iterations no matter what on the if ( converge )
c               stm. but that can cause problems if the load step
c               is incredibly small. convergence tests in mm03us,
c               mm03sb say it is non-coverged just because
c               changes are so small.
c
c
      converge    = .true.
      eps_compare = toler * max( abs(dep), abs(deq), eps_o )
      if ( debug ) then
         write(iout,9200) loop, r1, r2, del_dep, del_deq,
     &                    dep, deq, dep_old, deq_old, eps_compare
      end if
      if ( abs( dep - dep_old ) .gt. eps_compare ) converge = .false.
      if ( abs( deq - deq_old ) .gt. eps_compare ) converge = .false.
      if ( abs(r2) .gt. tol_yf*sigma_o ) converge = .false.
c      if ( converge .and. loop .gt. 1 ) go to 1000
      if ( converge ) go to 1000
c
      dep_old = dep
      deq_old = deq
c
      end do
c
c            the state update for dep, deq failed to converge in
c            the number of iterations allowed.
c
      gt_converge = 2
      return
c
c            compute final stress state on yield surface using
c            radial return along deviator direction of trial
c            elastic stress. converged values of p_new and
c            q_new above are returned thru the model common
c            for later use by consistent tangent routine.
c
 1000 continue
      dep_last_sub = dep
      deq_last_sub = deq
      smel   = -one * p_trial
      std(1) = newsig(1) - smel
      std(2) = newsig(2) - smel
      std(3) = newsig(3) - smel
      std(4) = newsig(4)
      std(5) = newsig(5)
      std(6) = newsig(6)
      if ( deq .gt. zero ) then
        scale = three * shear_mod * deq / q_trial
        newsig(1) = newsig(1) - scale * std(1)
        newsig(2) = newsig(2) - scale * std(2)
        newsig(3) = newsig(3) - scale * std(3)
        newsig(4) = newsig(4) - scale * std(4)
        newsig(5) = newsig(5) - scale * std(5)
        newsig(6) = newsig(6) - scale * std(6)
      end if
      newsig(1) = newsig(1) - bulk_mod * dep
      newsig(2) = newsig(2) - bulk_mod * dep
      newsig(3) = newsig(3) - bulk_mod * dep
c
c            new stresses on yield surface are in newsig. check that
c            the yield function is satisfied for the new stress state.
c            (the final check for bonehead mistakes in all the
c             iteration processes above).
c
      chk1 = mm03f( newsig, sbar_new, f_new, q1, q2, q3 )
      if ( abs(chk1) .gt. tol_yf * sigma_o ) then
         write(iout,9300) elenum, ptno, chk1
      end if
      if ( .not. debug ) return
c
      write(iout,9130) ( newsig(i), i = 1, 6 )
      write(iout,9140) chk1
      return
c
 9100 format(/,3x,'>> Starting iteration: ',i3,' to find dep, deq',
     & //,7x,'dep, deq, p_new, q_new: ',2f10.6,2f10.3 )
 9130 format(/,' >> Updated stresses:',6(/,3x,f10.3) )
 9140 format(/,' >> Gurson yield function: ',f12.3 )
 9200 format(/,3x,'>> Finished iteration: ',i3,' to find dep, deq',
     & //,7x,'r1, r2                   : ',2e14.6,
     &  /,7x,'del_dep, del_deq         : ',2e14.6,
     &  /,7x,'dep, deq                 : ',2e14.6,
     &  /,7x,'dep_old, deq_old         : ',2e14.6,
     &  /,7x,'strain tolerance         : ',e14.6 )
 9300 format(/,'>>>> Warning: routine mm03s. invalid update of',
     &       /,'              Gurson model. yield surface not',
     &       /,'              satisfied. element, point: ',2i8,
     &       /,'              yield function: ',e14.6)
c
      end subroutine mm03s
c ********************************************************************
c *                                                                  *
c *    routine  mm03q --  mises update  driver                       *
c *                                                                  *
c ********************************************************************
c
c
      subroutine mm03q( newsig )
c
      implicit none
c
c                   parameter declarations
c
      double precision
     &    newsig(*)
c
c              This routine accesses the variables defined in the
c              calling routine thru the f-90 host by association
c              feature.
c
c                   local variables in this routine
c
      double precision
     &  newyld_high, newyld_mid, toler, sig_tol, unused, eps_compare,
     &  deplas_low, deplas_high, resid_high,
     &  fl, fh, deplas, deplas_riddr, deplas_mid, resid_mid,
     &  sridder, deplas_new, resid_new, resid, term0, term1, term2,
     &  sm, resid_low
c
      integer loop, iconverge
c
      data toler, sig_tol, unused / 0.000001, 0.000001, -1.11e30 /
c
c            mises material with or without viscoplasticity
c
c            find the uniaxial plastic strain increment consistent
c            with the trial stress (yt) and the uniaxial stress (y)
c            vs. plastic strain (epspls) strain curve.
c
c            we have three models available for the inviscid uniaxial
c            response: (1) constant tangent modulus, (2) linear+power
c            law and (3) piece-wise linear. the uniaxial rsponse can
c            also be rate sensitive if the time increment for the step
c            is > 0 and the power-law rate parameters are set. The response
c            can be also rate dependent if the user has specified a set
c            of segmental curves (stress vs. plastic strain at various
c            plastic strain rates).
c
c            for temperature dependence, the stress-strain curve for
c            the specific gauss point temperature has been computed
c            before entering the material model.
c
c            an iterative scheme is required if the unixial hardening
c            response is a function of the accumulated plastic
c            strain. this occurs for any rate dependence or the
c            power-law or piece-wise (segmental) inviscid model.
c
c            we must solve for the plastic strain increment
c            using the scalar eqn:
c
c              resid function = q_trial - 3*g*deplas - yield(eplas) = 0
c
c            where yield is the uniaxial stress at the specified
c            level of temperature, plastic strain and plastic strain rate.
c            for constant hardening (rate independent), we have
c            yield = old-yield + deplas * h_constant.
c
c            we first bracket the value of deplas between 0 and
c            the value based on the smallest stress value on the
c            stress-strain curve and whether or not the
c            initial plastic modulus is >0 or <0. The segmental
c            curve can decrease then increase again to model
c            upper-lower yield behavior.
c
c            with these bounds in hand, we use a variant
c            of regula-falsi called Ridder's method (see Numerical
c            Recipes). The procedure keeps the updated estimates
c            for deplas bracketed while achieving very nearly
c            the quadratic convernce of a Newton method.
c
c            this algorithm guarantees that deplas can be found
c            for points very near breaks on the segemental
c            curve, segmental curves that 'stiffen' after a
c            softening response, for viscoplastic response with
c            a sudden rate change. a pure newton fails miserably in
c            these cases.
c
c            note: points on the segmental curve are plastic
c                  strain vs. stress. they are not req'd to be
c                  monotonically increasing
c
c            local debug segments are commented out to eliminate
c            repeated activation checks.
c
c            we call a subroutine to handle the different stress-strain
c            curve models. it returns the args listed and sets
c            a plastic modulus (h_inviscid) in the shared
c            contains variables. for rate dependent segmental, this
c            is the rate dependent value for the consistent tangent.
c
c            1. set the lower and upper bounds on deplas. verify the
c               root (deplas) is bracketed. determine if the upper/
c               lower bounds are in fact the solution. note the special
c               estimate if the initial plastic modulus is negative,
c               i.e., the segmental curve stress decreases immediately
c               from the yield stress to model an upper-lower yield
c               point response.
c
      vm_converge = 0
      loop        = 0
      eps_compare = toler * max( ebarp, eps_o )
      deplas_low  = zero
      resid_low   = q_trial - sbar
      if ( hprime_init .ge. zero ) then
        deplas_high = ( q_trial - sig_cur_min_val ) / threeg
      else
        deplas_high = ( q_trial - sig_cur_min_val ) /
     &                ( threeg + hprime_init )
      end if
      call mm03qq( deplas_high, newyld_high, resid_high, first_time, 1 )
c
      fl = resid_low
      fh = resid_high
      if ( (fl.gt.zero .and. fh.lt.zero) .or.
     &     (fl.lt.zero .and. fh.gt.zero) ) go to 100
      if ( abs(fl) .le. sig_tol*sbar  ) then
        newyld     = sbar
        deplas     = zero
        iconverge  = 1
      elseif ( abs(fh) .le. sig_tol*sbar ) then
        deplas     = deplas_high
        newyld     = newyld_high
        iconverge  = 2
      else
        write(iout,9300)
        vm_converge = 1
        return
      end if
      go to 1000
c
c            2. execute a maximum of 20 iterations of regula-falsi
c               with Ridder's improvements. debugs are commented
c               to eliminate repeated checks.
c
 100  continue
      deplas_riddr = unused
c
      do loop = 1, 20
c
c                     2a. evaluate the residual function at the
c                         current mid-point of bracketed range.
c                         compute Ridder's "magic" factor:)
c
      deplas_mid = half * ( deplas_low + deplas_high )
      call mm03qq( deplas_mid, newyld_mid, resid_mid, first_time, 1 )
      sridder = sqrt( resid_mid**2 - resid_low*resid_high )
      if ( sridder .eq. zero ) then
        write(iout,9310)
        vm_converge = 1
        return
      end if
c                     2b. compute the new estimate of deplas
c                         using Ridder's special update. check
c                         for convergence of deplas
c
      deplas_new = deplas_mid + (deplas_mid-deplas_low) *
     &             (sign(one,resid_low-resid_high)*resid_mid/sridder)
      if ( abs(deplas_new-deplas_riddr) .le. eps_compare ) then
         iconverge = 3
         deplas    = deplas_new
         go to 1000
      end if
c                     2c. evaluate residual function using the
c                         Ridder estimate for deplas. check
c                         convergence of the actual residual
c                         function.
c
      deplas_riddr =  deplas_new
      call mm03qq( deplas_riddr, newyld, resid_new, first_time, 1 )
      if ( abs(resid_new) .le. sig_tol*newyld ) then
         iconverge = 4
         deplas = deplas_new
         go to 1000
      end if
c                     2d. must do another iteration. bookeeping
c                         of high/low values to keep the root
c                         bracketed.
c
      if ( sign(resid_mid,resid_new) .ne. resid_mid ) then
        deplas_low = deplas_mid
        resid_low = resid_mid
        deplas_high = deplas_riddr
        resid_high = resid_new
      elseif ( sign(resid_low,resid_new) .ne. resid_low ) then
        deplas_high = deplas_riddr
        resid_high = resid_new
      elseif ( sign(resid_high,resid_new) .ne. resid_high ) then
        deplas_low = deplas_riddr
        resid_low  = resid_new
      else
        write(iout,9420)
        vm_converge = 1
        return
      end if
c
c                     2e. check convergence of new brackets on
c                         deplas - they may have collapsed to within
c                         the tolerance.
c
      if ( abs(deplas_high-deplas_low) .le. eps_compare ) then
        iconverge = 5
        deplas    = ( deplas_high + deplas_low ) * half
        go to 1000
      end if
c
c                     2f. all done with this iteration. fall thru
c                         at limit of iterations, set flag and
c                         return.
c
      end do
c
      vm_converge = 2
      return
c
c
c            3. iterations converged on deplas. except for convergence
c               checks at 3 and 5 we have the most current newyld
c               and plastic modulus. do a final using deplas
c               for them. then get rate dependent plastic modulus if
c               needed for power-law rate dependency
c               (mm03qq puts rate dependent stress in newyld,
c               static stress in inviscid_stress and static plastic
c               modulus into h_inviscid, all in common). for segmental
c               curve, get the final plastic modulus which could be
c               rate dependent. a final check in case deplas came
c               out negative - not a good outcome...
c
c
 1000 continue
      if ( iconverge .eq. 3 .or. iconverge .eq. 5 .or. segmental )
     &  call mm03qq( deplas, newyld, resid, first_time, 2 )
      if ( deplas .le. zero ) then
        if ( abs(deplas) .gt. 0.001*eps_o ) vm_converge = 3
      end if
      h      = h_inviscid
      if ( rate_depen ) then
         term0  = ( newyld / inviscid_stress ) ** (one - m_power )
         term1  = term0 * eps_ref*inviscid_stress / m_power / time_incr
         term2  = h_inviscid * newyld / inviscid_stress
         h      = term1 + term2
      end if
c
c               plastic strain increment available.  compute
c               the scale factor used to find the final
c               stresses. update the equivalent plastic strain
c               by the converged increment. set updated state
c               variables for calling routine. the matrix state
c               and the macroscopic state are same for a mises
c               (no void) material.
c
      h_new     = h
      ebarp_new = ebarp + deplas
      sbar_new  = newyld
      f_new     = zero
      q_new     = sbar_new
      p_new     = - ( newsig(1) + newsig(2) + newsig(3) ) * third
      dep       = zero
      deq       = deplas
      scale     = threeg * deplas / q_trial
      sm        = - p_new
      newsig(1) = newsig(1) - scale * ( newsig(1) - sm )
      newsig(2) = newsig(2) - scale * ( newsig(2) - sm )
      newsig(3) = newsig(3) - scale * ( newsig(3) - sm )
      newsig(4) = newsig(4) - scale * newsig(4)
      newsig(5) = newsig(5) - scale * newsig(5)
      newsig(6) = newsig(6) - scale * newsig(6)
c
      if ( .not. debug ) return
      write(iout,9120) scale, deplas, newyld, ebarp_new, h_new,
     &                 h_inviscid
      write(iout,9130) ( newsig(i), i = 1, 6 )
c
 9120 format(/,' >> final factors.  scale, deplas: ',f10.6,f12.9,
     & /,      '    newyld, ebarp_new: ',f10.3,f12.9,
     & /,      '    h_new, h_inviscid: ',2e10.3 )
 9130 format(/,' >> Updated stresses:',6(/,3x,f10.3) )
 9200 format(/,6x,'>> Starting iteration: ',i3,
     &        ' to find deplas...' )
 9300 format(/,3x,'>> Error: deplas not properly bracketed in',
     &  /,7x,'mises stress update.')
 9310 format(/,3x,'>> Error: divide by zero in mises update')
 9420 format(/,3x,'>> Error: failed to make new bracket',
     &       /,3x,'          for mises update of deplas')
 9500 format(2x,'   >> resid, deplas_new, newyld: ',f10.6,2f13.8)
 9700 format('** loops, iconverge, resid: ',i3,i3,f20.7)
c
      end subroutine mm03q
c ********************************************************************
c *                                                                  *
c *    routine  mm03qq --  mises update  driver                      *
c *                                                                  *
c ********************************************************************
c
c
      subroutine mm03qq( deplas, sig_bar, resid_now,
     &                   first_iter, caseh )
      implicit none
c
c                    parameter declarations (and functions called)

      double precision
     &   mm03is, mm03sc, deplas, sig_bar, resid_now
      integer caseh
      logical first_iter
c
c
c              This routine accesses the variables defined in the
c              calling routine thru the f-90 host by association
c              feature.
c
c                   local variables in this routine
c
      double precision
     & factor
c
c              caseh tells the segmental routine to compute (=2) or
c              not compute (=1) the plastic modulus. This takes
c              quite a bit of work for rate dependent response wo
c              we only do it of really necessary.
c
      if ( consth ) then
        inviscid_stress = sigma_o + h_constant * (ebarp+deplas)
        h_inviscid      = h_constant
      elseif ( power_law ) then
        inviscid_stress = mm03is( ebarp+deplas, sigma_o, e,
     &                            n_power, h_inviscid, iout )
      else
        inviscid_stress = mm03sc( ebarp+deplas, h_inviscid,
     &                            deplas, time_incr, caseh )
      end if
c
      if ( rate_depen ) then
        factor      = eps_ref * deplas/time_incr + one
        sig_bar     = inviscid_stress * factor**mpoweri
        resid_now   = q_trial - threeg * deplas - sig_bar
      else
        sig_bar    = inviscid_stress
        resid_now  = q_trial - threeg * deplas - sig_bar
      end if
c
      return
      end subroutine mm03qq
      end subroutine mm03

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
c                                   stress state for n+1.
c                                   of real load), this should be
c                                   final trial elastic stress state
c                                   computed in resolution of previous
c                                   load step. see mm03 for ordering.
c
c   history          1,2         -- 9x1 vector of history data for the
c                                   gauss point. passed in as history
c                                   state at 'n'.
c
c   history1         1,2         -- 9x1 vector of history data for the
c                                   gauss point. passed in as history
c                                   state for current estimate
c                                   of the solution at 'n+1'.
c
c   cep              2,2         -- 6x6 cepent [d]
c
c   iout             1,1         -- output device number for error messages.
c
c   span             1,1         -- number of elements in this block
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
c         depend on single/double precision implementations.
c
c
c   large-strain issues:
c   ====================
c
c         this routine is unaware of finite-strain vs. small-strain
c         issues. we assume that "rotation neutral" stresses
c         are passed as arguments.
c
      subroutine cnst3( element, gpn, iter, e, nu, q1,
     &                  q2, q3, nucleation, nuc_s_n, nuc_e_n, nuc_f_n,
     &                  stress_trial, history, history1, cep,
     &                  span, iout )
      implicit none
      include 'param_def'
c
c                   parameter declarations
c
      integer :: element, gpn, iter, iout, span
      double precision ::
     & stress_trial(mxvl,*), history(span,*),
     & history1(span,*), cep(mxvl,6,6), e(*), nu(*), q1(*),
     & q2(*), q3(*), nuc_s_n(*), nuc_e_n(*), nuc_f_n(*)
      logical :: nucleation(*)
c
c                   locally defined array, variables
c
      double precision ::
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
      double precision::
     & twothd, third, half, root32, zero, one, two, three, root_2_pi,
     & term1, term2, temp, b1, b2, d1, d2, d3, d4, h1, h2, h3, h4,
     & h5, h6, h7, h8, beta, ch, sh, pgp, pgq, pgsbar, pgf,
     & cap_a11, cap_a12, cap_b1, cap_c1, c5, c6, e10, e11, e12,
     & e13, m1, m2, c7, c8, c9, d10, d11, d12, d13,
     & cap_a21, cap_a22, cap_b2, cap_c2, denom, h10, h11, h12, h13,
     & con_1, con_2, con_3
      logical ::  nonlinear_points, debug, ldummy
      integer ::  iword(mxvl*2), state(mxvl), i, j, nonlin_point,
     &            inc_factor
      equivalence (iword, dword )
      data      twothd, third / 0.6666666667d00, 0.333333333d00 /
      data      half / 0.5d00 /, root32 / 1.224744871d00 /
      data      zero, one, two, three
     &           / 0.0d00, 1.0d00, 2.0d00, 3.0d00 /
      data      root_2_pi / 2.50663d00 /
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
      debug = .false.
      if( debug ) then
        write(iout,9000) element, gpn, iter, e(1), nu(1),
     &                   q1(1), q2(1), q3(1), nucleation(1),
     &                   nuc_s_n(1), nuc_e_n(1), nuc_f_n(1),
     &                   (stress_trial(1,i),i=1,6)
      end if
c
c              get material state variable for gauss point from
c              history vector.
c              two states are possible: = -1 implies strain point
c              is currently in a linear state. = 1 implies
c              strain point is currently undergoing plastic loading.
c
      nonlinear_points = .false.
      nonlin_point    = 0
!DIR$ VECTOR ALIGNED
      dword(1:span) = history1(1:span,6)
c
      inc_factor = 2
      j = 1
      do i = 1, span
       state(i) = iword(j)
       if( state(i) .eq. 1 ) then
          nonlinear_points = .true.
          if( nonlin_point .eq. 0 ) nonlin_point = i
       end if
       j = j + inc_factor
      end do
c
c              process linear strain points
c
      cep = zero
c
!DIR$ VECTOR ALIGNED
      do  i = 1, span
       if( state(i) .ne. -1 ) cycle
       c1(i)= (e(i)/((one+nu(i))*(one-two*nu(i))))
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
      end do
c
      if( .not. nonlinear_points ) return
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
!DIR$ VECTOR ALIGNED
      do i = 1, span
        if( state(i) .ne. 1 ) cycle
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
      end do
c
      if( iter .le. 1 ) then
!DIR$ VECTOR ALIGNED
        do i = 1, span
          if( state(i) .ne. 1 ) cycle
            dep(i) = zero
            deq(i) = zero
        end do
      end if
c
c
c
c            1) compute some frequently used constants
c
c            2) pressure and equivalent stress for macrostresses
c               at n+1. passed down from history state at (n).
c
c
      if( debug ) then
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
!DIR$ VECTOR ALIGNED
      do i = 1, span
        if( state(i) .ne. 1 ) cycle
        twog(i)   = two * shear_mod(i)
        threeg(i) = three * shear_mod(i)
        debarp(i) = ebarp(i) - ebarpn(i)
        sbar2(i)  = sbar(i) * sbar(i)
        sbar3(i)  = sbar(i)**3
        p(i)      = pn1(i)
        q(i)      = qn1(i)
      end do
c
c            3) process terms of trial elastic stress state at n+1.
c               get pressure, equivalent stress and yield surface
c               normal.
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
        if( state(i) .ne. 1 ) cycle
        sm(i)  = ( stress_trial(i,1) + stress_trial(i,2) +
     &             stress_trial(i,3) ) * third
        sx(i)  = stress_trial(i,1) - sm(i)
        sy(i)  = stress_trial(i,2) - sm(i)
        sz(i)  = stress_trial(i,3) - sm(i)
        qe(i)  = root32 * sqrt( sx(i)*sx(i) + sy(i)*sy(i) +
     &           sz(i)*sz(i) + two *
     &           ( stress_trial(i,4)*stress_trial(i,4) +
     &           stress_trial(i,5)*stress_trial(i,5) +
     &           stress_trial(i,6) * stress_trial(i,6) ) )
        pe(i)  = -sm(i)
        sl(i)  = three / two / qe(i)
        n(i,1) = sx(i) * sl(i)
        n(i,2) = sy(i) * sl(i)
        n(i,3) = sz(i) * sl(i)
        n(i,4) = stress_trial(i,4) * sl(i)
        n(i,5) = stress_trial(i,5) * sl(i)
        n(i,6) = stress_trial(i,6) * sl(i)
      end do
c
      if( debug ) then
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
!DIR$ VECTOR ALIGNED
      do i = 1, span
        if( state(i) .ne. 1 ) cycle
        a1(i) = (deq(i) * q(i) - dep(i) * p(i)) /
     &          ( one - f(i) )**2 / sbar(i)
        a2(i) = (deq(i) * q(i) - dep(i) * p(i)) * hprime(i)  /
     &          ( f(i) - one ) / sbar2(i)
c
c                4b)  a3, a4, a5 from Eq. (9)
c
        a3(i) = one / sbar(i) / ( f(i) - one )
        a4(i) = p(i)  +  dep(i) * bulk_mod(i)
        a5(i) = threeg(i) * deq(i)  -  q(i)
      end do
c
c                  4c)  b1 and b2 from Eq. (15). need to evaluate
c                       A(ebarp) for nucleation and its derivative
c                       wrt to ebarp.
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
        if( state(i) .ne. 1 ) cycle
        if( nucleation(i) ) then
          term1   = nuc_f_n(i) / nuc_s_n(i) / root_2_pi
          term2   = ( ( ebarp(i) - nuc_e_n(i) ) / nuc_s_n(i) )**2
          anuc(i) = term1 * exp(-half*term2)
          temp    = -nuc_f_n(i) * (ebarp(i)-nuc_e_n(i)) /
     &               nuc_s_n(i)**3 / root_2_pi
          anuc_prime(i) = temp * exp(-half*term2)
        else
          anuc(i)       = zero
          anuc_prime(i) = zero
        end if
      end do
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
        if( state(i) .ne. 1 ) cycle
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
!DIR$ VECTOR ALIGNED
      do i = 1, span
       if( state(i) .ne. 1 ) cycle
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
      end do
c
c                 7d) [D] (3) Eq. (17) and pg. CT-EF-8. multiply in
c                     the weight factor and determinant of coord.
c                     jacobian.
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
       if( state(i) .ne. 1 ) cycle
       con_1 = twog(i) * ( twog(i) * deq(i) / qe(i)  -  mqn(i) )
       cep(i,1,1) = (cep(i,1,1) + con_1 * n(i,1) * n(i,1))
       cep(i,2,1) = (cep(i,2,1) + con_1 * n(i,2) * n(i,1))
       cep(i,3,1) = (cep(i,3,1) + con_1 * n(i,3) * n(i,1))
       cep(i,4,1) = (cep(i,4,1) + con_1 * n(i,4) * n(i,1))
       cep(i,5,1) = (cep(i,5,1) + con_1 * n(i,5) * n(i,1))
       cep(i,6,1) = (cep(i,6,1) + con_1 * n(i,6) * n(i,1))
       cep(i,1,2) = (cep(i,1,2) + con_1 * n(i,1) * n(i,2))
       cep(i,2,2) = (cep(i,2,2) + con_1 * n(i,2) * n(i,2))
       cep(i,3,2) = (cep(i,3,2) + con_1 * n(i,3) * n(i,2))
       cep(i,4,2) = (cep(i,4,2) + con_1 * n(i,4) * n(i,2))
       cep(i,5,2) = (cep(i,5,2) + con_1 * n(i,5) * n(i,2))
       cep(i,6,2) = (cep(i,6,2) + con_1 * n(i,6) * n(i,2))
       cep(i,1,3) = (cep(i,1,3) + con_1 * n(i,1) * n(i,3))
       cep(i,2,3) = (cep(i,2,3) + con_1 * n(i,2) * n(i,3))
       cep(i,3,3) = (cep(i,3,3) + con_1 * n(i,3) * n(i,3))
       cep(i,4,3) = (cep(i,4,3) + con_1 * n(i,4) * n(i,3))
       cep(i,5,3) = (cep(i,5,3) + con_1 * n(i,5) * n(i,3))
       cep(i,6,3) = (cep(i,6,3) + con_1 * n(i,6) * n(i,3))
       cep(i,1,4) = (cep(i,1,4) + con_1 * n(i,1) * n(i,4))
       cep(i,2,4) = (cep(i,2,4) + con_1 * n(i,2) * n(i,4))
       cep(i,3,4) = (cep(i,3,4) + con_1 * n(i,3) * n(i,4))
       cep(i,4,4) = (cep(i,4,4) + con_1 * n(i,4) * n(i,4))
       cep(i,5,4) = (cep(i,5,4) + con_1 * n(i,5) * n(i,4))
       cep(i,6,4) = (cep(i,6,4) + con_1 * n(i,6) * n(i,4))
       cep(i,1,5) = (cep(i,1,5) + con_1 * n(i,1) * n(i,5))
       cep(i,2,5) = (cep(i,2,5) + con_1 * n(i,2) * n(i,5))
       cep(i,3,5) = (cep(i,3,5) + con_1 * n(i,3) * n(i,5))
       cep(i,4,5) = (cep(i,4,5) + con_1 * n(i,4) * n(i,5))
       cep(i,5,5) = (cep(i,5,5) + con_1 * n(i,5) * n(i,5))
       cep(i,6,5) = (cep(i,6,5) + con_1 * n(i,6) * n(i,5))
       cep(i,1,6) = (cep(i,1,6) + con_1 * n(i,1) * n(i,6))
       cep(i,2,6) = (cep(i,2,6) + con_1 * n(i,2) * n(i,6))
       cep(i,3,6) = (cep(i,3,6) + con_1 * n(i,3) * n(i,6))
       cep(i,4,6) = (cep(i,4,6) + con_1 * n(i,4) * n(i,6))
       cep(i,5,6) = (cep(i,5,6) + con_1 * n(i,5) * n(i,6))
       cep(i,6,6) = (cep(i,6,6) + con_1 * n(i,6) * n(i,6))
      end do
c
c               all done building 6x6 cep tangent.
c
      if( debug ) then
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
     & /,    '    e, nu : ',2f15.7,
     & /,    '    q1, q2, q3, nucleation: ',3f10.3, l10,
     & /,    '    nuc_s_n, nuc_e_n, nuc_f_n : ',3f10.3,
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
c
        end


