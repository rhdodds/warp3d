c *******************************************************************
c *                                                                 *
c *        material model # 5 -- cyclic plasticity model            *
c *                                                                 *
c *        mm05.f  -- last modified 9/29/2017 rhd                   *
c *                                                                 *
c *******************************************************************
c
      subroutine mm05(
     &  step, iter, felem, gpn, mxvl, hist_size, nstrs, nstrn, span,
     &  iout, signal_flag, adaptive_possible, cut_step_size_now,
     &  nonlin_hard, generalized_pl, sig_tol, mm_props, e_vec_np1,
     &  e_vec_n,
     &  nu_vec_np1, nu_vec_n, sigyld_gp_vec_np1, sigyld_gp_vec_n,
     &  h_gp_np1, h_gp_n, beta_gp_np1, beta_gp_n,
     &  delta_gp_np1, delta_gp_n, tau,
     &  sigyld_nl_vec_np1, sigyld_nl_vec_n, qu_nl_np1, qu_nl_n,
     &  bu_nl_np1, bu_nl_n, hu_nl_np1, hu_nl_n, gu_nl_np1, gu_nl_n,
     &  trial_elas_stress_np1, stress_n, stress_np1,
     &  deps, history_n, history_np1, do_nonlocal, nonlocal_state,
     &  maxnonlocal )
      implicit none
c
c                   parameter declarations
c                   ----------------------
      integer
     &  step, iter, felem, gpn, mxvl, hist_size, span,
     &  iout, nstrs, nstrn, maxnonlocal
c
      logical
     &   signal_flag, adaptive_possible, cut_step_size_now,
     &   nonlin_hard, generalized_pl, do_nonlocal
c
      double precision
     & sig_tol,  mm_props(mxvl,6), e_vec_np1(mxvl), e_vec_n(mxvl),
     & nu_vec_np1(mxvl), nu_vec_n(mxvl),
     & sigyld_gp_vec_np1(mxvl), sigyld_gp_vec_n(mxvl),
     & h_gp_np1(mxvl), h_gp_n(mxvl),
     & beta_gp_np1(mxvl), beta_gp_n(mxvl),
     & delta_gp_np1(mxvl), delta_gp_n(mxvl), tau(mxvl),
     & sigyld_nl_vec_np1(mxvl), sigyld_nl_vec_n(mxvl),
     & qu_nl_np1(mxvl), qu_nl_n(mxvl), bu_nl_np1(mxvl), bu_nl_n(mxvl),
     & hu_nl_np1(mxvl), hu_nl_n(mxvl), gu_nl_np1(mxvl), gu_nl_n(mxvl),
     & trial_elas_stress_np1(mxvl,nstrn),
     & stress_n(mxvl,nstrs), stress_np1(mxvl,nstrs), deps(mxvl,nstrn),
     & history_n(span,hist_size), history_np1(span,hist_size),
     & nonlocal_state(mxvl,maxnonlocal)
c
c               description of parameters
c               -------------------------
c
c     step              : current load step number
c     iter              : current newton iteration number
c     felem             : first element of the current block
c     gpn               : gauss point number being processed for block
c     mxvl              : maximum no. elements per block
c     hist_size         : number of history words per gauss point
c     nstrs             : number of stress terms (6 + extras)
c     nstrn             : number of incremental strain terms (6)
c     span              : number of elements in current block
c     iout              : write messates to this device number
c     signal_flag       : user wants notification messages for key
c                         events in material response
c     adaptive_possible : .true. if the material model may request
c                         immediate reduction of global load step size.
c                         no stress-histroy update required
c (*) cut_step_size_now : set .true. if material model wants immediate
c                         reduction of global load step size.
c                         no stress-histroy update required
c     nonlin_hard       : .true. if nonlinear_hardening option selected
c     generalized_pl    : .true. if generalized_plasticity option selected
c     mm_props          : material parameter values input by user for
c                         each element in block
c     e_vec_np1         : Young's modulus for each element in block at end
c                         of load step (n+1)
c     e_vec_n           : Young's modulus for each element in block at start
c                         of load step (n)
c     nu_vec_np1        : Poisson's ratio for each element in block at end
c                         of load step (n+1)
c     nu_vec_n          : Poisson's ratio for each element in block at start
c                         of load step (n)
c (*) trial_elas_stress_np1 : trial elastic stress vector to be used later by
c                         consistent tangent routine for model
c (**)stress_n          : stresses at start of load step (n) for all
c                         elements in block for this gauss point
c (*) stress_np1        : stresses at end of load step (n+1) for all
c                         elements in block for this gauss point
c     deps              : current estimate of strain increment over the
c                         load step (minus increment of thermal strain)
c (**)history_n         : history values at start of load step (n) for all
c                         elements in block for this gauss point
c (*) history_np1       : history values at end of load step (n+1) for all
c                         elements in block for this gauss point
c
c     Arrays used by the generalized plasticity model:
c     -----------------------------------------------
c
c     sigyld_gp_vec_np1 : yield stress for each element in block at end
c                         of load step (n+1)
c     sigyld_gp_vec_n   : yield stress for each element in block at start
c                         of load step (n)
c     h_gp_np1          : terminal hardening modulus of the gp model for
c                         each element in block at end of load step (n+1)
c     h_gp_n            : terminal hardening modulus of the gp model for
c                         each element in block at start of load step (n)
c     beta_gp_np1       : beta parameter of the gp model for each eleemnt in
c                         block at end of load step (n+1)
c     beta_gp_n         : beta parameter of the gp model for each eleemnt in
c                         block at start of load step (n)
c     delta_gp_np1      : delta parameter of the gp model for each element
c                         in block at end of load step (n+1)
c     delta_gp_n        : delta parameter of the gp model for each element
c                         in block at start of load step (n)
c     tau               : measure of kinematic/isotropic hardening in the
c                         gp model for each element in block
c
c     Arrays used by the nonlinear hardening model:
c     --------------------------------------------
c
c     sigyld_nl_vec_np1 : yield stress for each element in the block at end
c                         of load step (n+1)
c     sigyld_nl_vec_n   : yield stress for each element in the block at start
c                         of load step (n)
c     qu_nl_np1         : uniaxial q-value for each element in the block at
c                         end of load step (n+1)
c     qu_nl_n           : uniaxial q-value for each element in the block at
c                         start of load step (n)
c     bu_nl_np1         : uniaxial b-value for each element in the block at
c                         end of load step (n+1)
c     bu_nl_n           : uniaxial b-value for each element in the block at
c                         start of load step (n)
c     hu_nl_np1         : uniaxial h-value for each element in the block at
c                         end of load step (n+1)
c     hu_nl_n           : uniaxial h-value for each element in the block at
c                         start of load step (n)
c     gu_nl_np1         : uniaxial gamma value for each element in the block
c                         at end of the load step (n+1)
c     gu_nl_n           : uniaxial gamma value for each element in the block
c                       : at start of the load step (n)
c
c    (*)  values to be updated by this material model
c    (**) needs to be initialized on step 1
c
c   strain ordering:
c     deps-xx, deps-yy, deps-zz, gamma-xy, gamma-yz, gamma-xz
c
c   stress ordering (at n and n+1):
c     (1) sig-xx
c     (2) sig-yy
c     (3) sig-zz
c     (4) tau-xy
c     (5) tau-yz
c     (6) tau-xz
c     (7) total work density
c     (8) total plastic work density
c     (9) total plastic strain
c
c    mm_props ordering:
c      FA - Frederick-Armstrong, now referred to as the
c           nonlinear hardening option
c      GP - the Generalized Plasticity option
c
c      In August 2009 we changed the manual write up extensively
c      for the model to include the Generalized Plasticity option
c      and also modified the model input
c      parameters (names and values) to all conform to their
c      definitions in the 1-D formulations -- the reason is that the
c      user will invariably calibrate model parameters from a
c      1-D test result.
c
c      The computational code here has not been changed. It uses
c      a mix of 1-D and 3-D definitions for the material parameters.
c      The table below shows current status:
c
c      Nonlinear_hardening option
c
c         Code uses  Definition  Old input      New (1-D) input
c           q_bar      1-D        q_bar            q_u
c           b          3-D        b                b_u
c           h_bar      3-D        h_bar            h_u
c           gamma      3-D        gamma            gamma_u
c
c      The input translators (inmat.f) convert h_u, b_u, and
c      gamma_u to their 3-D values before storing in the material
c      properties table. The values pulled out here by the computational
c      code are the values cooresponding to columns 1 and 2 in the above
c      table on which the code was written.
c
c      Conversion formulas:  b = sqrt(2/3) * b_u
c                            h_bar = (2/3) * h_u
c                            gamma = sqrt(2/3) * gamma_u
c
c
c     (1) FA -> Q_bar,
c     (2) FA -> b,
c     (3) FA -> H_bar,
c     (4) FA and GP -> type (if type < 0 use GP model, otherwise use FA)
c     (5) FA -> gamma,
c     (6) FA and GP -> sigtol
c
c   history ordering (at n and n+1)
c     (1) lambda (consistency parameter)
c     (2) k(e_bar_p)  (size of yield surface). mises equiv stress
c                     is sqrt(3/2) * k
c     (3) e_bar_p     (accumulated plastic strain). tensorial
c                     definition. mises equivalent
c                     uniaxial strain = sqrt(2/3)*e_bar_p
c     (4) state  (integer indicating elastic or plastic state)
c     FA: (5) H      (kinematic hardening modulus)
c     GP: (5) (|| sigma-alpha || at end of step, used by tangent subroutine)
c     (6)-(11)   (backstress for kinematic hardening)
c     (12) s     (number of subincrements)
c
c     GP: (13) k_n   (size of yield surface at previous subincrement)
c     GP: (14) dHi   (change of isotropic hardening modulus btwn subincrements)
c
c                       declare local variables
c                       -----------------------
c
      logical :: debug, yield(mxvl), prior_linear(mxvl)
      integer :: i, j, iostat(mxvl), instat(mxvl)
      double precision :: shear_mod, c, a, b, delastic,
     &     shear_mod_vec(mxvl), alpha_n(mxvl,nstrn),
     &     trace_eps_np1(mxvl), yf_vec(mxvl), mises_equiv_stress,
     &     mises_eqiv_eps_pls
      double precision, parameter :: zero = 0.0d0, one = 1.0d0,
     &     two = 2.0d0, half = 0.5d0, root32 = dsqrt(3.0d0/2.0d0),
     &     root23 = dsqrt(2.0d0/3.0d0)
c
      double precision :: type ! see set from mm_props
c
      if( generalized_pl ) then
          call mm05_gp(
     &      step, iter, felem, gpn, mxvl, hist_size, nstrs, nstrn, span,
     &      iout, signal_flag, adaptive_possible, cut_step_size_now,
     &      mm_props, e_vec_np1, e_vec_n, nu_vec_np1, nu_vec_n,
     &      sigyld_gp_vec_np1, sigyld_gp_vec_n, h_gp_np1, h_gp_n,
     &      beta_gp_np1, beta_gp_n, delta_gp_np1, delta_gp_n, tau,
     &      trial_elas_stress_np1, stress_n, stress_np1,
     &      deps, history_n, history_np1 )
      end if
c
      if ( nonlin_hard ) then
          call mm05_fa(
     &      step, iter, felem, gpn, mxvl, hist_size, nstrs, nstrn, span,
     &      iout, signal_flag, adaptive_possible, cut_step_size_now,
     &      mm_props, e_vec_np1, nu_vec_np1, sigyld_nl_vec_np1,
     &      trial_elas_stress_np1, stress_n, stress_np1,
     &      deps, history_n, history_np1 )
      end if
c
c                       if needed, set nonlocal state variables for
c                       use by cohesive material model
c    (2) k(e_bar_p)  (size of yield surface). mises equiv stress
c                     is sqrt(3/2) * k
c     (3) e_bar_p     (accumulated plastic strain). tensorial
c                     definition. mises equivalent
c                     uniaxial strain = sqrt(2/3)*e_bar_p

c
       if( .not. do_nonlocal ) return
       do i = 1, span
        mises_equiv_stress = root32 * history_np1(i,2)
        mises_eqiv_eps_pls = root23 * history_np1(i,3)
        nonlocal_state(i,1) = mises_equiv_stress
        nonlocal_state(i,2) = mises_eqiv_eps_pls
       end do
c
       return
c
       end

c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_fa.f  -- last modified 2/05/04 by kbc
c *                                                                 *
c *******************************************************************
c
c
      subroutine mm05_fa(
     &  step, iter, felem, gpn, mxvl, hist_size, nstrs, nstrn, span,
     &  iout, signal_flag, adaptive_possible, cut_step_size_now,
     &  mm_props, e_vec, nu_vec, sigyld_vec,
     &  trial_elas_stress_np1, stress_n, stress_np1,
     &  deps, history_n, history_np1 )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &  step, iter, felem, gpn, mxvl, hist_size, span,
     &  iout, nstrs, nstrn
c
      logical
     &   signal_flag, adaptive_possible, cut_step_size_now
c
      double precision
     & mm_props(mxvl,5), e_vec(mxvl), nu_vec(mxvl),
     & sigyld_vec(mxvl), stress_n(mxvl,nstrs),
     & stress_np1(mxvl,nstrs), deps(mxvl,nstrn),
     & trial_elas_stress_np1(mxvl,nstrn), history_n(span,hist_size),
     & history_np1(span,hist_size)
c
c   mm_props ordering:
c     (1) Q_bar
c     (2) b
c     (3) H_bar
c     (4) type (switch for FA vs GP, not used at this level)
c     (5) gamma
c     (6) sig_tol
c
c   history ordering (at n and n+1)
c     (1) lambda (consistency parameter)
c     (2) k(ep)  (size of yield surface)
c     (3) ep     (accumulated plastic strain)
c     (4) state  (integer indicating elastic or plastic state)
c     (5) H      (kinematic hardening modulus)
c     (6)-(11)   (backstress for kinematic hardening)
c     (12) s     (number of subincrements)
c
c                       declare local variables
c                       -----------------------
c
       logical debug, yield(mxvl), prior_linear(mxvl)
       integer i, j, iostat(mxvl), instat(mxvl)
      double precision
     &     shear_mod_vec(mxvl), alpha_n(mxvl, nstrn),
     &     trace_eps_np1(mxvl),  yf_vec(mxvl), zero
      data zero / 0.0d00/
c
c       write(*,*) 'entering fa code'
c       write(*,*) (deps(1,j), j=1,6)
c
c               initialize stress_n and history_n on step 1
c                -------------------------------------------
c
      if ( step .eq. 1 ) then
        call mm05_fa_step1_set( stress_n, history_n, mm_props,
     &                      sigyld_vec, span, mxvl )
      end if
c
c    compute deviatoric stress, trial elastic state, evaluate mat. state
c    --------------------------------------------------------------------
c
       call mm05_fa_init(span, mxvl, iout, deps,iostat, felem,
     &   prior_linear, history_n, stress_n, e_vec, nu_vec, signal_flag,
     &   shear_mod_vec, alpha_n,
     &   trial_elas_stress_np1, debug, gpn,sigyld_vec, instat, yield,
     &   iter, trace_eps_np1, yf_vec, step )
c
c
c              compute stresses and internal variables at n+1
c              ----------------------------------------------
c
c       write(*,*) 'just before fa compute'
c       write(*,*) (deps(1,j), j=1,6)
c       write(*,*) 'sigy: ', sigyld_vec(1)
       call mm05_fa_compute(
     &    span, mxvl, debug, iout, history_n, history_np1,
     &    mm_props, sigyld_vec, yield, e_vec, iter,
     &    stress_n, stress_np1, shear_mod_vec,
     &    adaptive_possible, cut_step_size_now, signal_flag,
     &    gpn, felem, step, deps, prior_linear, trial_elas_stress_np1 )
c
       if ( cut_step_size_now ) then
          return
       end if
c
c                   update elastic elements
c                   -----------------------
c
       do i = 1, span
          if (yield(i) ) cycle
          history_np1(i,1)    = zero
          history_np1(i,2)    = history_n(i,2)
          history_np1(i,3)    = history_n(i,3)
          history_np1(i,5:11) = history_n(i,5:11)
          history_np1(i,12)   = zero
          stress_np1(i,1:6)   = trial_elas_stress_np1(i,1:6)
       end do
c
c
c    compute total update stresses, calculate the energy density
c    -----------------------------------------------------------
c
       call mm05_final(
     &                  span, mxvl, stress_n, stress_np1, e_vec,
     &                  nu_vec, shear_mod_vec, deps, instat,
     &                  history_np1, trace_eps_np1, debug, iout )
c
c
c                        plastic work update
c                        -------------------
c
      if ( iter .gt. 0 )
     &     call mm05_plastic_work( iout, span, mxvl, stress_n,
     &     stress_np1, yield, deps, nu_vec, e_vec, shear_mod_vec )
c
c       write(*,*) 'exiting fa code'
       return
c
       end

c *******************************************************************
c *                                                                 *
c *        material model # 5 -- cyclic plasticity model            *
c *                                                                 *
c *        mm05_fa_step1_set  -- last modified 1/23/04  by kbc      *
c *                                                                 *
c *        for step 1, initialize the stress and history vectors    *
c *                                                                 *
c *******************************************************************
      subroutine mm05_fa_step1_set( stress_n, history_n, mm_props,
     &                      sigyld_vec, span, mxvl )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &   mxvl, span
c
      double precision
     & mm_props(mxvl,5), sigyld_vec(mxvl),  stress_n(mxvl,*),
     &  history_n(span,*)
c
c                  local parameters
c                  ----------------
      integer i, iword(2)
c
      double precision
     &  root2third, zero, dword
c
      equivalence( dword, iword )
c
      data root2third, zero / 0.81649658, 0.0 /
c
      iword(1) = 3
      iword(2) = 0
c
c                 no longer zero _n for step 1. Could be pre-loaded
c                 w/ user-defined initial stresses
c
        do i = 1, span
            stress_n(i,7) = zero
            stress_n(i,8) = zero
            stress_n(i,9) = zero
c
            history_n(i,1) = zero
            history_n(i,2) = sigyld_vec(i)*root2third
            history_n(i,3) = zero
            history_n(i,4) = dword
            history_n(i,5) = mm_props(i,3)
            history_n(i,6) = zero
            history_n(i,7) = zero
            history_n(i,8) = zero
            history_n(i,9) = zero
            history_n(i,10) = zero
            history_n(i,11) = zero
            history_n(i,12) = zero
         end do
c
      return
c
      end


c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_fa_init  -- last modified 1/15/04 by kbc               *
c *                                                                 *
c *        compute deviatoric stress, trial elastic state,          *
c *        and evaluate material state at the beginning of each     *
c *        step                                                     *
c *                                                                 *
c *******************************************************************
      subroutine mm05_fa_init( span, mxvl, iout, deps, iostat, felem,
     & prior_linear, history_n, stress_n, e_vec, nu_vec, signal_flag,
     & shear_mod_vec, alpha_n, trial_elas_stress_np1, debug, gpn,
     & sigyld_vec, instat, yield, iter, trace_eps_np1, yf_vec, step )
      implicit none
c
c               parameter declarations
c               ----------------------
      integer span, mxvl, iout, iostat(*), instat(*), iter, step,
     &              gpn, felem
      logical debug, prior_linear(*), yield(*), signal_flag
      double precision
     & deps(mxvl, *), history_n(span,*), stress_n(mxvl,*),  yf_vec(*),
     & e_vec(*), shear_mod_vec(*), alpha_n( mxvl, *), nu_vec(*),
     & trial_elas_stress_np1(mxvl, *), sigyld_vec(*), trace_eps_np1(*)
c
c               local parameters
c               -----------------
      double precision
     & dword, deps_mean, de1, de2, de3, de4, de5, de6, trace_deps,
     & e, nu, g, een1, een2, een3, een4, een5, een6, e1, e2, e3,
     & e4, e5, e6, zero, one, two, three, yld_tol, rtse(6),
     & mrts( mxvl ), eps_mean, k_n(mxvl)
c
       integer iword(2), i, j
       equivalence ( dword, iword )
       data zero, one, two, three, yld_tol
     & / 0.0, 1.0, 2.0, 3.0, 0.00001 /
c
       do i = 1, span
          dword  = history_n(i, 4)
          iostat(i) = iword(1)
          prior_linear(i) = iostat(i) .eq. 3
c
c      deviatoric strain components
c
          trace_deps = deps(i,1) + deps(i,2) + deps(i,3)
          deps_mean  = trace_deps / three
          de1        = deps(i,1) - deps_mean
          de2        = deps(i,2) - deps_mean
          de3        = deps(i,3) - deps_mean
          de4        = deps(i,4)
          de5        = deps(i,5)
          de6        = deps(i,6)
c
c     elastic components of total strain at start of step
c
          e    = e_vec(i)
          nu   = nu_vec(i)
          g    = e/two/(one+nu)
          een1 = (stress_n(i,1)-nu*(stress_n(i,2)+ stress_n(i,3)))/e
          een2 = (stress_n(i,2)-nu*(stress_n(i,1)+ stress_n(i,3)))/e
          een3 = (stress_n(i,3)-nu*(stress_n(i,1)+ stress_n(i,2)))/e
          een4 = stress_n(i,4) / g
          een5 = stress_n(i,5) / g
          een6 = stress_n(i,6) / g
c
c     deviatoric elastic strain components at start of step plus
c     deviatoric strain increment over step
c
          trace_eps_np1(i) = een1 + een2 + een3 + trace_deps
          eps_mean         = (een1 + een2 + een3 ) / three
          e1               = (een1 - eps_mean) + de1
          e2               = (een2 - eps_mean) + de2
          e3               = (een3 - eps_mean) + de3
          e4               = een4 + de4
          e5               = een5 + de5
          e6               = een6 + de6
c
c    deviatoric components of trial elastic stress state
c
          shear_mod_vec(i)           = g
          trial_elas_stress_np1(i,1) = two*shear_mod_vec(i)*e1
          trial_elas_stress_np1(i,2) = two*shear_mod_vec(i)*e2
          trial_elas_stress_np1(i,3) = two*shear_mod_vec(i)*e3
          trial_elas_stress_np1(i,4) = shear_mod_vec(i)*e4
          trial_elas_stress_np1(i,5) = shear_mod_vec(i)*e5
          trial_elas_stress_np1(i,6) = shear_mod_vec(i)*e6
c
c     pull out backstress and other parameters at start of step
c
          alpha_n(i,1) = history_n(i,6)
          alpha_n(i,2) = history_n(i,7)
          alpha_n(i,3) = history_n(i,8)
          alpha_n(i,4) = history_n(i,9)
          alpha_n(i,5) = history_n(i,10)
          alpha_n(i,6) = history_n(i,11)
          k_n(i)       = history_n(i,2)
c
c     calculate relative trial stress ( deviatoric )
c
          rtse(1) = trial_elas_stress_np1(i,1) - alpha_n(i,1)
          rtse(2) = trial_elas_stress_np1(i,2) - alpha_n(i,2)
          rtse(3) = trial_elas_stress_np1(i,3) - alpha_n(i,3)
          rtse(4) = trial_elas_stress_np1(i,4) - alpha_n(i,4)
          rtse(5) = trial_elas_stress_np1(i,5) - alpha_n(i,5)
          rtse(6) = trial_elas_stress_np1(i,6) - alpha_n(i,6)
c
          mrts(i) = sqrt(rtse(1)**2+rtse(2)**2+rtse(3)**2 +
     &              two*(rtse(4)**2+rtse(5)**2+rtse(6)**2 ))

          yf_vec(i) = mrts(i) - k_n(i)
c
c     set various flags to indicate yielding
c         state variable = history_n(i,4)
c                        = 1, if gauss point is yield
c                        = 3, if not yielding
c
c     N.B. iteration zero is a pseudo iteration used to compute
c          stresses due to imposed displacement and temperature
c          changes.  yielding from a previously linear state is
c          not allowed for this iteration.
c
          instat(i) = 3
          yield(i) = .false.
c
          if ( yf_vec(i) .ge. yld_tol*k_n(i) ) then
              yield(i) = .true.
              instat(i) = 1
          end if
c
          if ( iter .gt. 0 ) then
             if(signal_flag) then
                if( history_n(i,3) .eq. zero .and. yield(i)) then
                   write(iout, 1010) felem+i-1, gpn
                end if
             end if
             cycle
          end if
c
          if ( prior_linear(i) )then
              yield(i) = .false.
              instat(i) = 3
              cycle
          end if

       end do
c

       return
c
 1000  format( 3e15.6, /, 3e15.6 )
 1010  format( 'element ', i6, ' at gauss point ', i3,
     &      '   begins yielding' )
       end
c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_fa_compute  -- last modified 12/18/2016 rhd         *
c *                                                                 *
c *        subincrementation based on truncation error with         *
c *        extrapolation                                            *
c *                                                                 *
c *******************************************************************
       subroutine mm05_fa_compute(
     &    span, mxvl, debug, iout, history_n, history_np1,
     &    mm_props, sigyld_vec, yield, e_vec, iter,
     &    stress_n, stress_np1, shear_mod_vec,
     &    adaptive_possible, cut_step_size_now, signal_flag,
     &    gpn, felem, step, deps, prior_linear, trial_stress )
       implicit none
c
c                    parameter declarations
c                    ----------------------
       integer span, mxvl, iout, gpn, felem, step, iter
       logical yield(*), debug, adaptive_possible, cut_step_size_now,
     &         signal_flag, prior_linear(*)
      double precision
     & history_n(span,*), history_np1(span,*), mm_props(mxvl,*),
     & sigyld_vec(*), stress_n(mxvl,*),deps(mxvl, *),
     & stress_np1(mxvl,*), shear_mod_vec(*),
     & trial_stress(mxvl, *), e_vec(*)
c
c                    local parameters
c                    ----------------
c
       integer i, j, s, nsubinc, l
       logical local_debug
      double precision
     & gamma, Q_bar, b, H_bar,
     & g,sig_tol, k_n, rel_tol, eta, error,a, bb,
     & mean_stress,  one, two, three, root2third,
     & lambda_1step, lambda_2step, lambda_nstep, mean_deps,
     & aps_n, aps_np1, aps_1step, aps_2step, aps_nstep,
     & alpha_n(6), alpha_np1(6), alpha_1step(6),  stress_1step(6),
     & alpha_2step(6), stress_2step(6), alpha_nstep(6), stress_nstep(6),
     & dev_stress_n(6), dev_stress_np1(6),diff(6), devdeps(6), sig_t(6),
     & norm_equiv_sm1, alpha_sm1(6)
c
       data  two, one, three,  root2third
     &   /  2.0, 1.0, 3.0, 0.81649658  /
c
       do i = 1, span
         local_debug = .false.
c
         if( .not. yield(i) ) cycle
c
         g       = shear_mod_vec(i)
         aps_n   = history_n(i, 3)
         k_n = history_n(i,2)
         Q_bar   = mm_props(i,1)
         b       = mm_props(i,2)
         H_bar   = mm_props(i,3)
         gamma   = mm_props(i,5)
         sig_tol = mm_props(i,6)
c
        rel_tol = sig_tol * sigyld_vec(i)
c
        mean_stress = ( stress_n(i, 1) + stress_n(i,2) +
     &                   stress_n(i, 3) )/three
        dev_stress_n(1) = stress_n(i,1) - mean_stress
        dev_stress_n(2) = stress_n(i,2) - mean_stress
        dev_stress_n(3) = stress_n(i,3) - mean_stress
        dev_stress_n(4) = stress_n(i,4)
        dev_stress_n(5) = stress_n(i,5)
        dev_stress_n(6) = stress_n(i,6)
c
        alpha_n(1:6) = history_n(i,6:11)
c
        mean_deps = ( deps(i, 1) + deps(i,2) +
     &                   deps(i, 3) )/three
c
        devdeps(1) = deps(i,1) - mean_deps
        devdeps(2) = deps(i,2) - mean_deps
        devdeps(3) = deps(i,3) - mean_deps
        devdeps(4) = deps(i,4)
        devdeps(5) = deps(i,5)
        devdeps(6) = deps(i,6)
c
c       compute the fraction of the step that is purely elastic and
c       update the stress_s vector to include the additional elastic
c       components.  also update devdeps.
c
c        write(*,*) (deps(i,j), j=1,6)
c
        call mm05_fa_elastic_fraction(  g, dev_stress_n, alpha_n,
     &      devdeps, eta, k_n, iout, felem, gpn, i, prior_linear(i))
c
c     use only one substep and return (for debugging)
        if(local_debug) then
           write(*,*) 'elastic fraction: ', eta
           write(iout,9005) (devdeps(j), j=1,6)
           write(iout,9004) (trial_stress(i,j), j=1,6)
        call mm05_fa_nsteps( 1, devdeps, dev_stress_n, alpha_n, aps_n,
     &      g, Q_bar, b, H_bar, gamma,  sigyld_vec(i),
     &      stress_1step, alpha_1step, aps_1step, lambda_1step,
     &      sig_t, alpha_sm1, norm_equiv_sm1, iout, gpn, felem, i, step,
     &      iter, signal_flag, adaptive_possible, cut_step_size_now )
c
         stress_np1(i,1:6) = stress_1step(1:6)
         history_np1(i,1) = lambda_1step
         history_np1(i,2) = root2third*(sigyld_vec(i) +
     &                  Q_bar*(one-exp( -b*aps_1step )))
         history_np1(i,3) = aps_1step
         history_np1(i,5) = H_bar
         history_np1(i,6:11) = alpha_1step(1:6)
         trial_stress(i,1:6) = sig_t(1:6)
         write(iout,9002) (stress_np1(i,j), j=1,6)
         write(iout,*) 'lambda ', lambda_1step
         return
        end if
c
        call mm05_fa_nsteps(1, devdeps, dev_stress_n, alpha_n, aps_n,
     &      g, Q_bar, b, H_bar, gamma,  sigyld_vec(i),
     &      stress_1step, alpha_1step, aps_1step, lambda_1step,
     &      sig_t, alpha_sm1, norm_equiv_sm1, iout,gpn, felem, i, step,
     &      iter, signal_flag, adaptive_possible, cut_step_size_now)
c
c        write(*,*) 'return from 1step calcs'
c
        call mm05_fa_nsteps(2, devdeps, dev_stress_n, alpha_n, aps_n,
     &      g, Q_bar, b, H_bar, gamma,  sigyld_vec(i),
     &      stress_2step, alpha_2step, aps_2step, lambda_2step,
     &      sig_t, alpha_sm1, norm_equiv_sm1, iout,gpn,felem, i,step,
     &      iter, signal_flag, adaptive_possible, cut_step_size_now)
c
c        write(*,*) 'return from 2step calcs'
c
         diff(1) = abs(stress_2step(1) - stress_1step(1))
         diff(2) = abs(stress_2step(2) - stress_1step(2))
         diff(3) = abs(stress_2step(3) - stress_1step(3))
         diff(4) = abs(stress_2step(4) - stress_1step(4))
         diff(5) = abs(stress_2step(5) - stress_1step(5))
         diff(6) = abs(stress_2step(6) - stress_1step(6))
         error = maxval(diff)
         nsubinc = ceiling(1.1*error/rel_tol)
         if (nsubinc .lt. two) nsubinc = two
         if(local_debug) then
           write(iout,9000) felem+i-1, gpn, step
           write(iout,9001) error, nsubinc
           write(iout,9002) (stress_1step(l), l=1,6)
           write(iout,9003) (stress_2step(l), l=1,6)
         end if
c
         history_np1(i,12) = nsubinc*one
         if(signal_flag) then
            write(iout,9006) felem+i-1, gpn, step, nsubinc
         end if
c
c        if error between 1 step and and 2 step solutions is small enough
c        that the necessary number of substeps is less than 3 then solution
c        is accurate enough; extrapolate and finish.  Otherwise
c        recompute solution using number of substeps estimated from
c        truncation error (nsubinc) then extrapolate 2 step and nsubinc step
c        solutions and finish
c
         if(nsubinc .lt. three ) then
c
          dev_stress_np1(1) = 2*stress_2step(1) - stress_1step(1)
          dev_stress_np1(2) = 2*stress_2step(2) - stress_1step(2)
          dev_stress_np1(3) = 2*stress_2step(3) - stress_1step(3)
          dev_stress_np1(4) = 2*stress_2step(4) - stress_1step(4)
          dev_stress_np1(5) = 2*stress_2step(5) - stress_1step(5)
          dev_stress_np1(6) = 2*stress_2step(6) - stress_1step(6)
          alpha_np1(1) = 2*alpha_2step(1) - alpha_1step(1)
          alpha_np1(2) = 2*alpha_2step(2) - alpha_1step(2)
          alpha_np1(3) = 2*alpha_2step(3) - alpha_1step(3)
          alpha_np1(4) = 2*alpha_2step(4) - alpha_1step(4)
          alpha_np1(5) = 2*alpha_2step(5) - alpha_1step(5)
          alpha_np1(6) = 2*alpha_2step(6) - alpha_1step(6)
          aps_np1 = 2*aps_2step - aps_1step
          lambda_nstep = lambda_2step
c
         else
c
           call mm05_fa_nsteps(nsubinc, devdeps, dev_stress_n, alpha_n,
     &         aps_n, g, Q_bar, b, H_bar, gamma,  sigyld_vec(i),
     &         stress_nstep, alpha_nstep, aps_nstep, lambda_nstep,
     &         sig_t, alpha_sm1, norm_equiv_sm1, iout, gpn, felem,
     &         i, step, iter, signal_flag, adaptive_possible,
     &         cut_step_size_now )
c
            a = (one*nsubinc)/(nsubinc-2)
            bb = (one*2)/(nsubinc-2)
            dev_stress_np1(1) = a*stress_nstep(1) - bb*stress_2step(1)
            dev_stress_np1(2) = a*stress_nstep(2) - bb*stress_2step(2)
            dev_stress_np1(3) = a*stress_nstep(3) - bb*stress_2step(3)
            dev_stress_np1(4) = a*stress_nstep(4) - bb*stress_2step(4)
            dev_stress_np1(5) = a*stress_nstep(5) - bb*stress_2step(5)
            dev_stress_np1(6) = a*stress_nstep(6) - bb*stress_2step(6)
            alpha_np1(1) = a*alpha_nstep(1) - bb*alpha_2step(1)
            alpha_np1(2) = a*alpha_nstep(2) - bb*alpha_2step(2)
            alpha_np1(3) = a*alpha_nstep(3) - bb*alpha_2step(3)
            alpha_np1(4) = a*alpha_nstep(4) - bb*alpha_2step(4)
            alpha_np1(5) = a*alpha_nstep(5) - bb*alpha_2step(5)
            alpha_np1(6) = a*alpha_nstep(6) - bb*alpha_2step(6)
            aps_np1  = a*aps_nstep - bb*aps_2step
c
         end if
c
         stress_np1(i,1:6) = dev_stress_np1(1:6)
c
c       update history
c
         history_np1(i,1) = lambda_nstep
         history_np1(i,2) = root2third*(sigyld_vec(i) +
     &                  Q_bar*(one-exp( -b*aps_np1 )))
         history_np1(i,3) = aps_np1
         history_np1(i,5) = H_bar
         history_np1(i,6:11) = alpha_np1(1:6)
         trial_stress(i,1:6) = sig_t(1:6)
c
       end do
c
       return
c
 9000  format('Debug mm05_fa_compute:  element ', i7, ' at gp ', i3,
     &   ' step number ', i7 )
 9001  format('computed error: ', e15.6, '  nsubinc: ', i6)
 9002  format('stress_1step: ', / 3e15.6 / 3e15.6)
 9003  format('stress_2step: ', / 3e15.6 / 3e15.6)
 9004  format('trial stress: ', / 3e15.6 / 3e15.6)
 9005  format('devdeps: ', / 3e15.6 / 3e15.6)
 9006  format('element ', i7, ' at gp ', i3, ' step ', i7,
     &        ' used ', i6, ' substeps')
c
       end
c
c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_fa_nsteps  -- last modified 3/22/05 by kbc          *
c *                                                                 *
c *                                                                 *
c *******************************************************************
c
        subroutine mm05_fa_nsteps(nsubinc, devdeps, dev_stress_n,
     &      alpha_n, aps_n, g,  Q_bar, b, H_bar, gamma, sigyld,
     &      stress_j, alpha_j, aps_j, lambda,
     &      sig_t, alpha_sm1, norm_equiv_j, iout,gpn,felem, i, step,
     &      iter, signal_flag, adaptive_possible, cut_step_size_now)
      implicit none
c
c    parameter declarations
c    ----------------------
c
      integer iout, gpn, felem, step, iter, i, nsubinc
      logical signal_flag, adaptive_possible, cut_step_size_now
      double precision
     &  devdeps(*), dev_stress_n(*), alpha_n(*), aps_n, g, H_bar,
     &  Q_bar, sigyld, stress_j(*), alpha_j(*), aps_j, lambda, gamma,
     &  sig_t(*), alpha_sm1(*), b, norm_equiv_j
c
c    local parameters
c    ----------------
c
      integer s, l,j
      logical local_debug
      double precision
     &  one, equiv_t(6), norm_equiv_t, k_j, two, d, root2third,
     &  equiv(6), k_np1, H_np1, lg, kappa, beta, ro
c
       data  two,  one, root2third
     &   /   2.0, 1.0, 0.81649658  /
c
           stress_j(1:6) = dev_stress_n(1:6)
           alpha_j(1:6)  = alpha_n(1:6)
           aps_j = aps_n
           local_debug = .false.
c
           if(local_debug) then
              write(iout,*) 'devdeps: ', (devdeps(j), j=1,6)
              write(iout,*) 'dev_stress_n: ', (dev_stress_n(j), j=1,6)
              write(iout,*) 'alpha_n: ', (alpha_n(j), j=1,6)
              write(iout,*) 'aps_n: ', aps_n
              write(iout,*) 'g, Q_bar, b, H_bar, gamma, sigyld: ',
     &             g, Q_bar, b, H_bar, gamma, sigyld
           end if
c
           d = one/nsubinc
c
           do s = 1, nsubinc
c
              sig_t(1) = stress_j(1) + two*g*devdeps(1)*d
              sig_t(2) = stress_j(2) + two*g*devdeps(2)*d
              sig_t(3) = stress_j(3) + two*g*devdeps(3)*d
              sig_t(4) = stress_j(4) +     g*devdeps(4)*d
              sig_t(5) = stress_j(5) +     g*devdeps(5)*d
              sig_t(6) = stress_j(6) +     g*devdeps(6)*d
c
c              write(*,*) 'about to call lambda solve'
c              write(*,*) 'sigy: ', sigyld
c             write(*,*) 'sigt_t before: ', (sig_t(j), j=1,6)
c
              call  mm05_fa_lambda_solve( iout, sigyld, alpha_j,
     &          sig_t, equiv, aps_j,gamma, lambda, k_np1, Q_bar, b,
     &          H_bar, lg, g, adaptive_possible, cut_step_size_now,
     &          signal_flag, gpn, felem, step, i, H_np1, iter, s )
c
               kappa = two*g*lambda/k_np1
               beta  = lambda*H_np1/k_np1
               ro    = lg*(one+kappa)+beta
c
               stress_j(1) = sig_t(1) - kappa/ro * equiv(1)
               stress_j(2) = sig_t(2) - kappa/ro * equiv(2)
               stress_j(3) = sig_t(3) - kappa/ro * equiv(3)
               stress_j(4) = sig_t(4) - kappa/ro * equiv(4)
               stress_j(5) = sig_t(5) - kappa/ro * equiv(5)
               stress_j(6) = sig_t(6) - kappa/ro * equiv(6)
c
               alpha_j(1) = (alpha_j(1) + beta/ro *equiv(1))/lg
               alpha_j(2) = (alpha_j(2) + beta/ro *equiv(2))/lg
               alpha_j(3) = (alpha_j(3) + beta/ro *equiv(3))/lg
               alpha_j(4) = (alpha_j(4) + beta/ro *equiv(4))/lg
               alpha_j(5) = (alpha_j(5) + beta/ro *equiv(5))/lg
               alpha_j(6) = (alpha_j(6) + beta/ro *equiv(6))/lg
c
            end do
c
9000  format('Debug mm05_fa_compute:  element ', i7, ' at gp ', i3,
     &   ' step number ', i6 )

       end


c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_fa_lambda_solve  -- last modified 2/06/04 by kbc    *
c *                                                                 *
c *                                                                 *
c *******************************************************************
c
       subroutine mm05_fa_lambda_solve(
     &          iout, sig_y, alpha_s, sig_t, equiv, aps_n, gamma,
     &          lambda, k_np1, Q_bar, b, H_bar, lg, g,
     &          adaptive_possible, cut_step_size_now, signal_flag,
     &          gpn, felem, step, i, H_np1, iter, sub )
       implicit none
c
c    parameter declarations
c    ----------------------
c
      integer iout, gpn, felem, step, i, iter, sub
      logical signal_flag, adaptive_possible, cut_step_size_now
      double precision
     &    sig_y, alpha_s(*), sig_t(*), equiv(*), aps_n,
     &    lambda, k_np1,  Q_bar, b, H_bar, lg, gamma, g, H_np1
c
c  the following parameters should be current upon exiting this
c  subroutine, even though they are functions of the current lambda
c  only, they are not recalculated before the stress and history
c  variables are updated:  equiv, aps_n, lambda, k_np1, lg, H_np1
c
c
c    local  parameter declarations
c    -----------------------------
c
      logical debug
      integer j, l, max_iter

      double precision
     &    norm_equiv, res, dres, eq_dot_sigt, tol, aps_np1,
     &    dk_np1, dH_np1, zero, one, two, root2third
c
      data max_iter, zero, one, two, root2third
     &   / 15, 0.0, 1.0, 2.0, 0.8164956809 /
      debug = .false.
c
      if( iter .gt. 1 ) then
        tol = 0.000001
      else
        tol = 0.001
      end if
c
      if( debug ) then
          write(iout, *) 'inside lambda solve'
c          write(iout, 9000)  (equiv(l), l=1,6)
          write(iout, 9000)  (sig_t(l), l=1,6)
c          write(iout, 9000)  (alpha_s(l), l=1,6)
c          write(*,*) 'sigy ', sig_y
      end if
c
      lambda = zero
      do j = 1, max_iter
c
c       compute current values of isotropic and kinematic
c       hardening functions
c
        aps_np1 = aps_n + lambda
        k_np1  = root2third*(sig_y + Q_bar*(one-exp( -b*aps_np1 )))
        dk_np1 = root2third*b*Q_bar*(exp( -b*aps_np1 ))
        H_np1  = H_bar
        dH_np1 = zero
        lg     = one + lambda*gamma
c
        equiv(1) = lg*sig_t(1)-alpha_s(1)
        equiv(2) = lg*sig_t(2)-alpha_s(2)
        equiv(3) = lg*sig_t(3)-alpha_s(3)
        equiv(4) = lg*sig_t(4)-alpha_s(4)
        equiv(5) = lg*sig_t(5)-alpha_s(5)
        equiv(6) = lg*sig_t(6)-alpha_s(6)
c
        norm_equiv = sqrt(equiv(1)**2+equiv(2)**2+equiv(3)**2 +
     &               two*(equiv(4)**2+equiv(5)**2+equiv(6)**2 ))
c
        res = lg*( k_np1 + two*g*lambda) + lambda*H_np1 - norm_equiv
c        write(*,*) 'residual: ', res, abs(res), tol*sig_y
c
c       test for convergence, updating history variables
c       and exiting loop if calculations have converged.
c
        if ( abs(res) .lt. tol*sig_y ) then
           if ( lambda .lt. zero ) then
             call mm05_find_pos_lambda(iout, sig_y, sig_t, equiv,
     &              aps_n, gamma, lambda, k_np1, Q_bar, b, H_bar, lg,
     &              g, alpha_s, adaptive_possible, cut_step_size_now,
     &              H_np1, tol, aps_np1 )
           end if
           aps_n = aps_np1
           return
        end if
c
        eq_dot_sigt = equiv(1)*sig_t(1) +
     &                equiv(2)*sig_t(2) +
     &                equiv(3)*sig_t(3) +
     &            two*equiv(4)*sig_t(4) +
     &            two*equiv(5)*sig_t(5) +
     &            two*equiv(6)*sig_t(6)
c
        dres = gamma*(k_np1+two*g*lambda)+lg*(dk_np1 + two*g) +
     &         H_np1 + lambda*dH_np1 - gamma*eq_dot_sigt/norm_equiv
        lambda = lambda - res/dres
c
      end do
c
c     deal with convergence failure (exceeding max. iterations)
c     call for a step size reduction if allowed, otherwise kill
c     the analysis
c
      if ( j .gt. max_iter ) then
         if ( adaptive_possible ) then
             cut_step_size_now  = .true.
             write( iout, 9020 ) felem+i-1, gpn
             return
         else
             write(iout, 9200 )
             call die_abort
             end if
         end if
c
      return
c
 9000  format( 3e15.6, /, 3e15.6 )
 9020  format( 'element ', i7, ' at gauss point ', i3,
     &  ' requests step size reduction:  mm05 failed to converge' // )
 9200  format( '>> Fatal Error: routine mm05_fa_lambda_solve.' //
     &         ' newton iterations failed to converge.',
     &         ' Job terminated.' // )
c
      end
c
c
c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_find_pos_lambda  -- last modified 2/05/04 by kbc    *
c *                                                                 *
c *                                                                 *
c *******************************************************************
c
       subroutine mm05_find_pos_lambda(iout, sigy, sig_t, equiv,
     &            aps_n, gamma, lambda, k_np1, Q, b, H_bar, lg, g,
     &            alpha_n, adaptive_possible, cut_step_size_now,
     &            H_np1, tol, aps_np1 )

       implicit none
c
c    parameter declarations
c    ----------------------
c
      integer iout
      logical adaptive_possible, cut_step_size_now
      double precision
     &    sigy, alpha_n(*), sig_t(*), equiv(*), aps_n, tol,
     &    lambda, k_np1,  Q, b, H_bar, lg, gamma, g, H_np1,
     &    aps_np1
c
c    local  parameter declarations
c    -----------------------------
c
      logical local_debug
      integer sc_count, sc_max, bs_count, bs_max, nr_count, nr_max

      double precision
     &    dk_np1, dH_np1,  res_i, res_f, lambda_i, lambda_f, res,
     &    zero, one, two, root2third
c
      data zero, one, two, root2third, sc_max, bs_max, nr_max
     &   / 0.0, 1.0, 2.0, 0.8164956809, 4, 10, 10 /
c
      local_debug = .false.
c
c   set starting interval: lambda_i = 0 and find lambda_f such
c   the residual at the interval endpoints have opposite signs
c
      lambda_i = zero
      sc_count = 0
      bs_count = 0
      nr_count = 0
c
      call mm05_compute_residual( lambda_i, aps_n, sig_t,
     &     alpha_n, gamma, sigy, H_bar, Q, b, g, k_np1,
     &     H_np1, lg, res_i, equiv, aps_np1, dk_np1, dH_np1, iout)
c
      lambda_f = -1*lambda
c
      call mm05_compute_residual( lambda_f, aps_n, sig_t,
     &     alpha_n, gamma, sigy, H_bar, Q, b, g, k_np1,
     &     H_np1, lg, res_f, equiv, aps_np1, dk_np1, dH_np1, iout)
c
       do while ( sign( one, res_i)*sign(one, res_f) .gt. zero )
         sc_count = sc_count+1
         if( sc_count > sc_max ) then
           write(iout, 9000)
           if( adaptive_possible ) then
             cut_step_size_now = .true.
             write(iout, 9010)
             return
           else
             write(iout, 9020 )
             call die_abort
           end if
         end if
         lambda_f = 2*lambda_f
      call mm05_compute_residual( lambda_f, aps_n, sig_t,
     &     alpha_n, gamma, sigy, H_bar, Q, b, g, k_np1,
     &     H_np1, lg, res_f, equiv, aps_np1, dk_np1, dH_np1, iout)
      end do
c
c   find lambda > 0 using Newton-Raphson iterations with interval
c   bisection if the lambda at any iteration falls outside of the
c   known interval
c
      do while ( bs_count < bs_max )
        bs_count = bs_count+1
        lambda = lambda_i + 0.5*(lambda_f - lambda_i)
        do while (nr_count < nr_max )
c
           call mm05_compute_residual( lambda, aps_n, sig_t,
     &          alpha_n, gamma, sigy, H_bar, Q, b, g, k_np1,
     &          H_np1, lg, res, equiv, aps_np1, dk_np1, dH_np1,
     &          iout)
c
           if ( abs(res) < tol*sigy ) then
              return
           end if
c
           nr_count = nr_count +1
           call mm05_update_lambda(lambda, res, equiv, gamma,
     &          k_np1, g, H_np1, dk_np1, dH_np1, sig_t )
           if( lambda < lambda_i .or. lambda > lambda_f ) then
             exit
           end if
        end do
        lambda = lambda_i + 0.5*(lambda_f - lambda_i)
c
         call mm05_compute_residual( lambda, aps_n, sig_t,
     &        alpha_n, gamma, sigy, H_bar, Q, b, g, k_np1,
     &        H_np1, lg, res, equiv, aps_np1, dk_np1, dH_np1, iout)
c
        if( (sign( one, res_i)*sign(one, res)) .gt. zero ) then
           lambda_i = lambda
c
           call mm05_compute_residual( lambda_i, aps_n, sig_t,
     &        alpha_n, gamma, sigy, H_bar, Q, b, g, k_np1,
     &        H_np1, lg, res_i, equiv, aps_np1, dk_np1, dH_np1, iout)
c
        else
           lambda_f = lambda
c
           call mm05_compute_residual( lambda_f, aps_n, sig_t,
     &        alpha_n, gamma, sigy, H_bar, Q, b, g, k_np1,
     &        H_np1, lg, res_f, equiv, aps_np1, dk_np1, dH_np1, iout)
c
        end if
c
        nr_count = 0
      end do
c
      if( adaptive_possible ) then
        cut_step_size_now = .true.
        write(iout, 9040)
        return
      else
         write(iout, 9030 )
         call die_abort
      end if

      return
c
 9000 format ('Exceeded maximum attempts to find suitable bracket
     &         for interval bisection' )
 9010 format ('requesting step size reduction' )
 9020 format ('Fatal Error: can not find lambda > 0' )
 9030 format ('Fatal Error: exceeded maximum interval bisection
     &         iterations without finding lambda > 0')
 9040 format ('Exceeded maximum interval bisection attempts,
     &         requesting step size reduction. ' // )
      end

c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *      mm05_compute_residual --last modified 2/05/04 by kbc       *
c *                                                                 *
c *                                                                 *
c *******************************************************************
c
      subroutine mm05_compute_residual( lambda, aps_n, sig_t,
     &     alpha_n, gamma, sigy, H_bar, Q_bar, b, g, k_np1, H_np1,
     &     lg, res, equiv, aps_np1, dk_np1, dH_np1, iout )
      implicit none
c
c    parameter declarations
c    ----------------------
      integer iout
      double precision
     & lambda, aps_n, sig_t(*), alpha_n(*),
     & gamma, sigy, H_bar, b, Q_bar, g, lg,
     & k_np1, dk_np1, H_np1, dH_np1, res, aps_np1, equiv(*)
c
c    local parameters
c    ----------------
c
      double precision
     &    zero, one, two, root2third, norm_equiv
c
      data zero, one, two, root2third
     &   / 0.0, 1.0, 2.0, 0.8164956809 /
c
       aps_np1 = aps_n + lambda
       k_np1  = root2third*(sigy + Q_bar*(one-exp( -b*aps_np1 )))
       dk_np1 = root2third*b*Q_bar*(exp( -b*aps_np1 ))
       H_np1  = H_bar
       dH_np1 = zero
       lg     = one + lambda*gamma
c
       equiv(1) = lg*sig_t(1)-alpha_n(1)
       equiv(2) = lg*sig_t(2)-alpha_n(2)
       equiv(3) = lg*sig_t(3)-alpha_n(3)
       equiv(4) = lg*sig_t(4)-alpha_n(4)
       equiv(5) = lg*sig_t(5)-alpha_n(5)
       equiv(6) = lg*sig_t(6)-alpha_n(6)
c
       norm_equiv = sqrt(equiv(1)**2+equiv(2)**2+equiv(3)**2 +
     &               two*(equiv(4)**2+equiv(5)**2+equiv(6)**2 ))
c
       res = lg*( k_np1 + two*g*lambda) + lambda*H_np1 - norm_equiv
c
      return
c
      end

c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *      mm05_update_lambda --last modified 2/05/04 by kbc          *
c *                                                                 *
c *                                                                 *
c *******************************************************************
c
      subroutine mm05_update_lambda(lambda, res, equiv, gamma,
     &          k_np1, g, H_np1, dk_np1, dH_np1, sig_t )
c
c    parameter declarations
c    ----------------------
      double precision
     & lambda, res, equiv(*), gamma, k_np1, g, H_np1, dk_np1, dH_np1,
     & norm_equiv, eq_dot_sigt, one, two, sig_t(*)
c
      data one, two  / 1.0, 2.0 /
c
       norm_equiv = sqrt(equiv(1)**2+equiv(2)**2+equiv(3)**2 +
     &               two*(equiv(4)**2+equiv(5)**2+equiv(6)**2 ))

       eq_dot_sigt = equiv(1)*sig_t(1)  +
     &                equiv(2)*sig_t(2) +
     &                equiv(3)*sig_t(3) +
     &            two*equiv(4)*sig_t(4) +
     &            two*equiv(5)*sig_t(5) +
     &            two*equiv(6)*sig_t(6)

       dres = gamma*(k_np1+two*g*lambda)+
     &        (one+lambda*gamma)*(dk_np1 + two*g) +  H_np1 +
     &        lambda*dH_np1 - gamma*eq_dot_sigt/norm_equiv
c
      lambda = lambda - res/dres
c
       return
c
       end

c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *     mm05_fa_elastic_fraction  -- last modified 1/23/04 by kbc   *
c *                                                                 *
c *    this subroutine computes the fraction of the the strain      *
c *    increment that produces purely elastic response and updates  *
c *    the stress_n vector to include the stress caused by this     *
c *    purely elastic strain fraction                               *
c *                                                                 *
c *******************************************************************
c
       subroutine mm05_fa_elastic_fraction( g, stress_n, alpha_n,
     &                                 devdeps, eta, k_n, iout,
     &                                 felem, gpn, i, prior_linear )
       implicit none
c
c    parameter declarations
c    ----------------------
c
      logical prior_linear
      integer iout, felem, gpn, i
      double precision
     &    g, stress_n(*), alpha_n(*), devdeps(*), eta, k_n
c
c    local  parameter declarations
c    -----------------------------
c
      logical debug
      double precision
     &    a, b, c, zero, two, temp, tol
c
      data zero, two, tol / 0.0, 2.0, 0.00001/
      debug = .false.
c
        eta = zero
c
c        write(iout, 9000) (devdeps(i), i=1,6)
c
           a = two*g*((stress_n(1)-alpha_n(1))*devdeps(1) +
     &                (stress_n(2)-alpha_n(2))*devdeps(2) +
     &                (stress_n(3)-alpha_n(3))*devdeps(3) +
     &                (stress_n(4)-alpha_n(4))*devdeps(4) +
     &                (stress_n(5)-alpha_n(5))*devdeps(5) +
     &                (stress_n(6)-alpha_n(6))*devdeps(6) )
c
           b = g*g*(4*(devdeps(1)**2 + devdeps(2)**2 + devdeps(3)**2 )+
     &              2*(devdeps(4)**2 + devdeps(5)**2 + devdeps(6)**2 ))
c
c   the quantity 'c' is -(||sigma_n - alpha_n||^2 - k_n^2) = -yf,
c   i.e, the negative of the yield function.  The following uses c = 0 if the
c   previous step was plastic, or if it was elastic but 0 < yf < tol*k_n
c   where tol is the yield tolerance used to determine yielding in subroutine
c   mm05_init.
c
c
           c = 0
c
           if( prior_linear ) then
              c =  ((stress_n(1)-alpha_n(1))**2 +
     &                      (stress_n(2)-alpha_n(2))**2 +
     &                      (stress_n(3)-alpha_n(3))**2 +
     &                 two*((stress_n(4)-alpha_n(4))**2 +
     &                      (stress_n(5)-alpha_n(5))**2 +
     &                      (stress_n(6)-alpha_n(6))**2 ) )
            temp = k_n - sqrt(c)
            c = k_n*k_n - c
            if( c .lt. zero ) then
               if( temp .lt. tol*k_n) then
                  c  = zero
               else
                 write(iout, 9200 )  felem+i-1, gpn
                 call die_abort
               end if
             end if
          end if

          if(debug) then
              write(iout, *) a, b, c, k_n
          end if

c
           if( c .eq. zero ) then
             if( a .gt. zero ) then
                eta = zero
             else
                eta = -two*a/b
             end if
           else
              eta = -a/b + sqrt( (a/b)**2 + c/b )
              if( eta .lt. zero ) then
                 write(iout, 9300 ) felem+i-1, gpn
                 call die_abort
              end if
           end if

c
c          update deviatoric stress to be sig_n + eta*2G*deps
c
           stress_n(1) = stress_n(1) + eta*2*g*devdeps(1)
           stress_n(2) = stress_n(2) + eta*2*g*devdeps(2)
           stress_n(3) = stress_n(3) + eta*2*g*devdeps(3)
           stress_n(4) = stress_n(4) + eta*g*devdeps(4)
           stress_n(5) = stress_n(5) + eta*g*devdeps(5)
           stress_n(6) = stress_n(6) + eta*g*devdeps(6)
c
c          update deviatoric strain increment to exclude elastic part
c
           devdeps(1:6) = (1-eta)*devdeps(1:6)
c
      return
c
 9000  format( 3e15.6, /, 3e15.6 )
 9200  format( '>> Fatal Error: routine mm05_elastic_fraction:',
     &        '  element ', i6, ' gauss point ', i3, ' c < 0 ',
     &         ' Job terminated.' // )
 9300  format( '>> Fatal Error: routine mm05_elastic_fraction:',
     &        '  element ', i6, ' gauss point ', i3, ' eta < 0',
     &         ' Job terminated.' // )
c
      end

c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *         subroutine mm05_debug -- prints out variables           *
c *                                                                 *
c *******************************************************************
c
c
      subroutine mm05_debug(
     &  step, iter, felem, gpn, mxvl, hist_size, nstrs, nstrn, span,
     &  iout, signal_flag, adaptive_possible, cut_step_size_now,
     &  mm_props, e_vec, nu_vec, sigyld_vec,
     &  trial_elas_stress_np1, stress_n, stress_np1,
     &  deps, history_n, history_np1 )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &  step, iter, felem, gpn, mxvl, hist_size,
     &  span, iout, nstrs, nstrn
c
      logical
     &   signal_flag, adaptive_possible, cut_step_size_now
c
      double precision
     & mm_props(mxvl,5), e_vec(mxvl), nu_vec(mxvl),
     & sigyld_vec(mxvl), stress_n(mxvl,nstrs),
     & stress_np1(mxvl,nstrs), deps(mxvl,nstrn),
     & trial_elas_stress_np1(mxvl,nstrn), history_n(span,hist_size),
     & history_np1(span,hist_size)
c
c                          local parameters
c
      integer i, j
c
c                         print out all input variables
c
          write( iout, * ) 'debug mode - echoing input variables'
          write( iout, 9000 ) step, iter, felem, gpn, mxvl, hist_size,
     &          span, nstrs, nstrn,
     &          signal_flag, adaptive_possible, cut_step_size_now
c          do i =1, span
           i = 1
           write( iout, 9010 )  e_vec(i), nu_vec(i),
     &        ( mm_props(i,j), j=1,5 ),
     &        ( stress_n(i,j), j=1, nstrs ),
     &        ( stress_np1(i,j), j=1, nstrs ),
     &        ( deps(i,j), j=1, nstrn ),
     &        ( trial_elas_stress_np1(i,j), j=1, nstrn ),
     &        ( history_n(i,j), j=1, hist_size ),
     &        ( history_np1(i,j), j=1, hist_size )
c          end do
c
       return
c
c
 9000 format( '  step = ', i3, ',  iteration = ', i4,',  felem = ', i3,
     & /, '  gpn = ', i3,  ',  mxvl = ', i3, ',  hist_size = ', i3,
     & /, '  span = ', i3, ',  nstrs = ', i3, ',  nstrn = ', i3,
     & /, '  signal_flag = ', l3, '  adapt. = ', l3, '  cut_step = ',
     &    l3 )
 9010 format( '  e_vec, nu_vec = ', 2e12.3,
     & /,     '  mm_props = ', 5e12.3,
     & /,     '  stress @ n   = ',    3e15.6, / , 17x, 3e15.6,
     & /,                                         17x, 3e15.6,
     & /,     '  stress @ n+1 = ',  3e15.6, / , 17x, 3e15.6,
     & /,                                         17x, 3e15.6,
     & /,     '  deps         = ',  3e15.6, / , 17x, 3e15.6,
     & /,     '  trial stress = ',  3e15.6, / , 17x, 3e15.6,
     & /,     '  history @ n  = ', 4e15.6, / , 17x, 4e15.6,
     & /,                                      17x, 3e15.6,
     & /,     '  history @ n+1= ', 4e15.6, / , 17x, 4e15.6,
     & /,                                      17x, 4e15.6,
     & /,                                      17x,  e15.6 )
c
       end


c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        subroutine cnst5 -- computes consistent tangent          *
c *                                                                 *
c *******************************************************************
c
c
      subroutine cnst5(
     &  span, felem, gpn, iter, iout, mxvl, nstrn,
     &  e_vec, nu_vec, mm_props, sig_trial, history_n,
     &  history_np1, stress_np1, dmat,
     &  h_gp_np1, beta_gp_np1, delta_gp_np1, tau )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &  span, felem, gpn, iter, iout, mxvl, nstrn
c
      double precision
     & mm_props(mxvl,6), e_vec(mxvl), nu_vec(mxvl),
     & sig_trial(mxvl,nstrn), history_n(span,*),
     & history_np1(span,*), dmat(mxvl,nstrn,nstrn),
     & stress_np1(mxvl,nstrn), h_gp_np1(mxvl), beta_gp_np1(mxvl),
     & delta_gp_np1(mxvl), tau(mxvl)
c
c               description of parameters
c               -------------------------
c
c     step              : current load step number
c     iter              : current newton iteration number. iter 1 is
c                         for application of the "real" load increment.
c     felem             : first element of the current block
c     gpn               : gauss point number being processed for block
c     mxvl              : maximum no. elements per block
c     nstrn             : number of strain-stress components (=6)
c     span              : number of elements in current block
c     iout              : write messages to this device number
c     mm_props          : material parameter values input by user for
c                         each element in block
c     e_vec             : Young's modulus for each element in block
c     nu_vec            : Poisson's ratio for each element in block
c (#) trial_elas_stress_np1 : trial elastic stress vector defined by stress
c                             update routine
c                         consistent tangent routine for model
c     history_n         : history values at start of load step (n) for all
c                         elements in block for this gauss point
c     history_np1       : history values at end of load step (n+1) for all
c                         elements in block for this gauss point
c (!) stress_np1        : current estimate of 6 stress components for
c                         end of step (see ordering below)
c (*) dmat              : 6x6 (symmetric) tangent (consistent) for
c                         this gauss point for each element of block
c                         (see stress ordering below)
c
c    (*)  values to be updated by this material model
c    (!)  for finite strain computations, these are unrotated
c         Cauchy stresses
c    (#)  used by constitutive update procedures based on some
c         form of elastic predictor - return mapping algorithm.
c         the contents of this array are set by the corresponding
c         stress update routine for the model. for finite
c         strain computations, the contents will be unrotated
c         Cauchy stress terms of some form as set by the stress
c         update routine.
c
c   Finite strain issues:
c     this constitutive model only sees stresses that have already
c     been rotation "neutralized" by WARP3D. this routine can
c     thus operate simply as "small strain" theory. WARP3D will
c     "rotate" the compute [D] matrices on return for finite strain-
c     rotation effects.
c
c   strain ordering:
c     deps-xx, deps-yy, deps-zz, gamma-xy, gamma-yz, gamma-xz
c
c   stress ordering (at n and n+1):
c     sig-xx, sig-yy, sig-zz, tau-xy, tau-yz, tau-xz
c
c                      local variables
c                      ----------------
c
      integer ::i, iword(2), state, l, m , j, t
      logical :: yield(mxvl), debug
      double precision ::
     &     c1, c2, c3, c4, fact, zero, one, two, dword,
     &     n(6), k, dk, H, dH, kappa, dkappa, beta, dbeta,
     &     ro, zeta, a, b, c, equiv(6), norm_equiv, g, bulk,
     &     twothird, onethird, mean_n, half, Q_bar, H_bar,
     &     gamma, lambda, aps, lg, temp, rt2third,
     &     three, four, mean_sig, type

       data zero,one,two, twothird, onethird, half,
     &         rt2third, three, four  /0.0d0, 1.0d0, 2.0d0,
     &     0.66666666d0, 0.33333333d0, 0.5d0, 0.81649658d0, 3.0d0,
     &     4.0d0 /
c
       equivalence(dword, iword)
c
       type = mm_props(1,4)
c
c                if GP model is indicated call subroutine for GP
c                tangent and return. otherwise - compute tangent
c                for FA model
c
       if( type .lt. zero ) then
        call  cnst5_gp(
     &  span, felem, gpn, iter, iout, mxvl, nstrn,
     &  e_vec, nu_vec, h_gp_np1, beta_gp_np1,
     &  delta_gp_np1, tau, sig_trial, history_n, history_np1,
     &  stress_np1, dmat )
        return
       end if
c
       do i = 1, span
         dword = history_np1(i,4)
         state = iword(1)
         yield(i) = .false.
         if( state .eq. 1 ) yield(i) = .true.
      end do
c
      do i = 1, span
         if( yield(i) ) cycle
         dmat(i,1,4) = zero
         dmat(i,1,5) = zero
         dmat(i,1,6) = zero
         dmat(i,2,4) = zero
         dmat(i,2,5) = zero
         dmat(i,2,6) = zero
         dmat(i,3,4) = zero
         dmat(i,3,5) = zero
         dmat(i,3,6) = zero
         dmat(i,4,1) = zero
         dmat(i,4,2) = zero
         dmat(i,4,3) = zero
         dmat(i,4,5) = zero
         dmat(i,4,6) = zero
         dmat(i,5,1) = zero
         dmat(i,5,2) = zero
         dmat(i,5,3) = zero
         dmat(i,5,4) = zero
         dmat(i,5,6) = zero
         dmat(i,6,1) = zero
         dmat(i,6,2) = zero
         dmat(i,6,3) = zero
         dmat(i,6,4) = zero
         dmat(i,6,5) = zero
c
         c1 = e_vec(i)/((one+nu_vec(i))*(one-two*nu_vec(i)))
         c2 = (one-nu_vec(i))*c1
         c3 = ((one-two*nu_vec(i))/two)*c1
         c4 = nu_vec(i)*c1
c
         dmat(i,1,1)= c2
         dmat(i,2,2)= c2
         dmat(i,3,3)= c2
         dmat(i,4,4)= c3
         dmat(i,5,5)= c3
         dmat(i,6,6)= c3
         dmat(i,1,2)= c4
         dmat(i,1,3)= c4
         dmat(i,2,1)= c4
         dmat(i,3,1)= c4
         dmat(i,2,3)= c4
         dmat(i,3,2)= c4
      end do
c
      do i = 1, span
       if( .not. yield(i) ) cycle
       fact = 1.0d00
       g     = e_vec(i)/(two*(1+nu_vec(i)))
       bulk  = e_vec(i)/(three*(one-two*nu_vec(i)))
       Q_bar = mm_props(i,1)
       b     = mm_props(i,2)
       H_bar = mm_props(i,3)
       gamma = mm_props(i,5)
c
       lambda = history_np1(i,1)
       k      = history_np1(i, 2)
       aps    = history_np1(i,3)
       H      = history_np1(i, 5)
       dk     = rt2third*b*Q_bar*(exp( -b*aps ))
       dH     = zero
       lg     = one + lambda*gamma
c
       kappa = two*g*lambda/k
       beta  = lambda*H/k
       ro    = lg*(one+kappa)+beta
c
       mean_sig = ( stress_np1(i,1)+stress_np1(i,2) +
     &                    stress_np1(i,3) )/three
c
       equiv(1) = stress_np1(i,1)- mean_sig - history_np1(i,6)
       equiv(2) = stress_np1(i,2)- mean_sig - history_np1(i,7)
       equiv(3) = stress_np1(i,3)- mean_sig - history_np1(i,8)
       equiv(4) = stress_np1(i,4)-history_np1(i,9)
       equiv(5) = stress_np1(i,5)-history_np1(i,10)
       equiv(6) = stress_np1(i,6)-history_np1(i,11)
c
       norm_equiv = sqrt(equiv(1)**2+equiv(2)**2+equiv(3)**2 +
     &              two*(equiv(4)**2+equiv(5)**2+equiv(6)**2 ))
c
       n(1) = equiv(1)/norm_equiv
       n(2) = equiv(2)/norm_equiv
       n(3) = equiv(3)/norm_equiv
       n(4) = equiv(4)/norm_equiv
       n(5) = equiv(5)/norm_equiv
       n(6) = equiv(6)/norm_equiv
c
       dkappa = two*g*(k-lambda*dk)/k**2
       dbeta  = H/k + (lambda/k**2)*(dH*k - H*dk)
c
       temp = n(1)*sig_trial(i,1) +
     &        n(2)*sig_trial(i,2) +
     &        n(3)*sig_trial(i,3) +
     &        two*n(4)*sig_trial(i,4) +
     &        two*n(5)*sig_trial(i,5) +
     &        two*n(6)*sig_trial(i,6)
c
       zeta   = (two*g*lg) / (dk + two*g + H + lambda*dH +
     &           gamma*(k + four*g*lambda + dk*lambda - temp))
c
       a = two*g/ro*(lg+beta)
       b = k*zeta/ro*(gamma*(kappa*(one+kappa) - lambda*dkappa )
     &     + dbeta*kappa - beta*dkappa - dkappa )
       c = gamma*zeta*kappa/ro
c
c                      diagonal terms
c
       mean_n  = onethird*(n(1) + n(2) + n(3) )
       dmat(i,1,1) = (bulk + twothird*a + (b*n(1)-c*sig_trial(i,1))*
     &                                     (n(1)-mean_n) )
       dmat(i,2,2) = (bulk + twothird*a + (b*n(2)-c*sig_trial(i,2))*
     &                                     (n(2)-mean_n) )
       dmat(i,3,3) = (bulk + twothird*a + (b*n(3)-c*sig_trial(i,3))*
     &                                     (n(3)-mean_n) )
       dmat(i,4,4) = (half*a + (b*n(4) - c*sig_trial(i,4))*n(4))
       dmat(i,5,5) = (half*a + (b*n(5) - c*sig_trial(i,5))*n(5))
       dmat(i,6,6) = (half*a + (b*n(6) - c*sig_trial(i,6))*n(6))
c
c                     off diagonal terms, symmetrized
c
       dmat(i,1,2) = (bulk - onethird*a +
     &                half*((b*n(1)-c*sig_trial(i,1))*(n(2)-mean_n)+
     &                      (b*n(2)-c*sig_trial(i,2))*(n(1)-mean_n)))
       dmat(i,1,3) = (bulk - onethird*a +
     &                half*((b*n(1)-c*sig_trial(i,1))*(n(3)-mean_n)+
     &                      (b*n(3)-c*sig_trial(i,3))*(n(1)-mean_n)))
       dmat(i,1,4) = half*( (b*n(1) - c*sig_trial(i,1))*n(4) +
     &                  (b*n(4) - c*sig_trial(i,4))*(n(1)-mean_n))
       dmat(i,1,5) = half*( (b*n(1) - c*sig_trial(i,1))*n(5) +
     &                  (b*n(5) - c*sig_trial(i,5))*(n(1)-mean_n))
       dmat(i,1,6) = half*( (b*n(1) - c*sig_trial(i,1))*n(6) +
     &                  (b*n(6) - c*sig_trial(i,6))*(n(1)-mean_n))
c
       dmat(i,2,3) = (bulk - onethird*a +
     &              half*((b*n(2)-c*sig_trial(i,2))*(n(3)-mean_n)+
     &                    (b*n(3)-c*sig_trial(i,3))*(n(2)-mean_n)))
       dmat(i,2,4) = half*((b*n(2) - c*sig_trial(i,2))*n(4) +
     &                 (b*n(4) - c*sig_trial(i,4))*(n(2)-mean_n))
       dmat(i,2,5) = half*((b*n(2) - c*sig_trial(i,2))*n(5) +
     &                 (b*n(5) - c*sig_trial(i,5))*(n(2)-mean_n))
       dmat(i,2,6) = half*((b*n(2) - c*sig_trial(i,2))*n(6) +
     &                 (b*n(6) - c*sig_trial(i,6))*(n(2)-mean_n))
c
       dmat(i,3,4) = half*((b*n(3) - c*sig_trial(i,3))*n(4) +
     &                 (b*n(4) - c*sig_trial(i,4))*(n(3)-mean_n))

       dmat(i,3,5) = half*((b*n(3) - c*sig_trial(i,3))*n(5) +
     &                 (b*n(5) - c*sig_trial(i,5))*(n(3)-mean_n))

       dmat(i,3,6) = half*((b*n(3) - c*sig_trial(i,3))*n(6) +
     &                 (b*n(6) - c*sig_trial(i,6))*(n(3)-mean_n))
c
       dmat(i,4,5) = half*((b*n(4) - c*sig_trial(i,4))*n(5) +
     &                        (b*n(5) - c*sig_trial(i,5))*n(4))
       dmat(i,4,6) = half*((b*n(4) - c*sig_trial(i,4))*n(6) +
     &                        (b*n(6) - c*sig_trial(i,6))*n(4))
c
       dmat(i,5,6) = half*((b*n(5) - c*sig_trial(i,5))*n(6) +
     &                        (b*n(6) - c*sig_trial(i,6))*n(5))
c
c                  symmetry
c
       dmat(i,2,1) = dmat(i,1,2)
       dmat(i,3,1) = dmat(i,1,3)
       dmat(i,4,1) = dmat(i,1,4)
       dmat(i,5,1) = dmat(i,1,5)
       dmat(i,6,1) = dmat(i,1,6)
c
       dmat(i,3,2) = dmat(i,2,3)
       dmat(i,4,2) = dmat(i,2,4)
       dmat(i,5,2) = dmat(i,2,5)
       dmat(i,6,2) = dmat(i,2,6)
c
       dmat(i,4,3) = dmat(i,3,4)
       dmat(i,5,3) = dmat(i,3,5)
       dmat(i,6,3) = dmat(i,3,6)
c
       dmat(i,5,4) = dmat(i,4,5)
       dmat(i,6,4) = dmat(i,4,6)
c
       dmat(i,6,5) = dmat(i,5,6)
c
      end do
c
      return
c
 1000 format ( 6e12.4)
      end
c
c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_gp.f  -- last modified 2/05/04 by jcs
c *                                                                 *
c *******************************************************************
c
c
      subroutine mm05_gp(
     &  step, iter, felem, gpn, mxvl, hist_size, nstrs, nstrn, span,
     &  iout, signal_flag, adaptive_possible, cut_step_size_now,
     &  mm_props, e_vec_np1, e_vec_n, nu_vec_np1,
     &  nu_vec_n, sigyld_vec_np1, sigyld_vec_n, h_gp_np1, h_gp_n,
     &  beta_gp_np1, beta_gp_n, delta_gp_np1, delta_gp_n, tau,
     &  trial_elas_stress_np1, stress_n, stress_np1,
     &  deps, history_n, history_np1 )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &  step, iter, felem, gpn, mxvl, hist_size, span,
     &  iout, nstrs, nstrn
c
      logical
     &   signal_flag, adaptive_possible, cut_step_size_now
c
      double precision
     & mm_props(mxvl,5), e_vec_np1(mxvl), e_vec_n(mxvl),
     & nu_vec_np1(mxvl), nu_vec_n(mxvl),
     & sigyld_vec_np1(mxvl), sigyld_vec_n(mxvl),
     & h_gp_np1(mxvl), h_gp_n(mxvl), beta_gp_np1(mxvl), beta_gp_n(mxvl),
     & delta_gp_np1(mxvl), delta_gp_n(mxvl), tau(mxvl),
     & stress_n(mxvl,nstrs), stress_np1(mxvl,nstrs), deps(mxvl,nstrn),
     & trial_elas_stress_np1(mxvl,nstrn), history_n(span,hist_size),
     & history_np1(span,hist_size)
c
c               description of parameters
c               -------------------------
c   mm_props ordering:
c     (6) sig_tol
c
c   history ordering (at n and n+1)
c     (1) lambda (consistency parameter)
c     (2) k(ep)  (size of yield surface)
c     (3) ep     (accumulated plastic strain)
c     (4) state  (integer indicating elastic or plastic state)
c     (5) (|| sigma-alpha || at start of step, used by tangent subroutine)
c     (6)-(11)   (backstress for kinematic hardening)
c     (12) number of subincrements
c     (13) k_n   (size of yield surface at previous subincrement)
c     (14) dHi   (change of isotropic hardening modulus btwn subincrements)
c
c                       declare local variables
c                       -----------------------
c
       logical debug, yield(mxvl), prior_linear(mxvl)
       integer i, j, iostat(mxvl), instat(mxvl)
      double precision
     &     g_vec_n(mxvl), g_vec_np1(mxvl), alpha_n(mxvl, nstrn),
     &     trace_eps_np1(mxvl), zero, one, rse(6),
     &     hi_n, hi_n1, hk_n, hk_n1, lk, eps_n, htol,
     &     two, three, twothird, root2third
       data zero, one, htol, two, three, twothird, root2third
     &      / 0.0, 1.0, 0.000001, 2.0, 3.0, 0.66666666, 0.81649658 /
c
c               initialize stress_n and history_n on step 1
c                -------------------------------------------
c
      if ( step .eq. 1 ) then
        call mm05_gp_step1_set( stress_n, history_n, mm_props,
     &    sigyld_vec_n, span, mxvl )
      end if
c
c    compute deviatoric stress, trial elastic state, evaluate mat. state
c    --------------------------------------------------------------------
c
       call mm05_gp_init(span, mxvl, iout, deps, iostat, felem,
     &   prior_linear, history_n, stress_n, signal_flag,
     &   g_vec_np1, g_vec_n, alpha_n, trial_elas_stress_np1, debug,
     &   gpn, sigyld_vec_np1, sigyld_vec_n, instat, yield, iter,
     &   trace_eps_np1, step, h_gp_np1, h_gp_n, tau,
     &   e_vec_np1, e_vec_n, nu_vec_np1, nu_vec_n )
c
c
c              compute stresses and internal variables at n+1
c              ----------------------------------------------
c
       call mm05_gp_compute(
     &   span, mxvl, debug, iout, history_np1, history_n,
     &   mm_props, sigyld_vec_np1, sigyld_vec_n, yield,
     &   iter, stress_np1, stress_n, adaptive_possible,
     &   cut_step_size_now, signal_flag, g_vec_np1, g_vec_n,
     &   h_gp_np1, h_gp_n, delta_gp_np1, delta_gp_n,
     &   beta_gp_np1, beta_gp_n, tau, gpn, felem, step, deps,
     &   prior_linear, trial_elas_stress_np1 )
c
       if ( cut_step_size_now ) then
          return
       end if
c
c                   update elastic elements
c                   -----------------------
c
       do i = 1, span
          if ( yield(i) ) cycle
c
          eps_n   = history_n(i,3)
          hi_n    = twothird*tau(i)*h_gp_n(i)
          hi_n1   = twothird*tau(i)*h_gp_np1(i)
          hk_n    = twothird*(one-tau(i))*h_gp_n(i)
          hk_n1   = twothird*(one-tau(i))*h_gp_np1(i)
          lk    = one
          if ( abs( hk_n ) .gt. htol ) lk = hk_n1/hk_n
c
          history_np1(i,1)    = zero
          history_np1(i,2)    = root2third*sigyld_vec_np1(i)+hi_n1*eps_n
          history_np1(i,3)    = history_n(i,3)
          history_np1(i,6:11) = lk*history_n(i,6:11)
          history_np1(i,12)   = zero
          history_np1(i,13)   = root2third*sigyld_vec_n(i)+hi_n*eps_n
          history_np1(i,14)   = hi_n1 - hi_n
          stress_np1(i,1:6)   = trial_elas_stress_np1(i,1:6)
c
          rse(1) = stress_np1(i,1) - history_np1(i,6)
          rse(2) = stress_np1(i,2) - history_np1(i,7)
          rse(3) = stress_np1(i,3) - history_np1(i,8)
          rse(4) = stress_np1(i,4) - history_np1(i,9)
          rse(5) = stress_np1(i,5) - history_np1(i,10)
          rse(6) = stress_np1(i,6) - history_np1(i,11)
          history_np1(i,5) = sqrt( rse(1)**2 + rse(2)**2 + rse(3)**2
     &         + two*(rse(4)**2 + rse(5)**2 + rse(6)**2 ))
c
       end do
c
c      call mm05_gp_debug(
c     &  step, iter, felem, gpn, mxvl, hist_size, nstrs, nstrn, span,
c     &  iout, signal_flag, adaptive_possible, cut_step_size_now,
c     &  mm_props, e_vec_np1, nu_vec_np1, sigyld_vec_np1,
c     &  e_vec_n, nu_vec_n, sigyld_vec_n, h_gp_np1, h_gp_n,
c     &  beta_gp_np1, beta_gp_n, delta_gp_np1, delta_gp_n, tau,
c     &  trial_elas_stress_np1, stress_n, stress_np1,
c     &  deps, history_n, history_np1 )
c
c
c
c
c
c    compute total update stresses, calculate the energy density
c    -----------------------------------------------------------
c
       call mm05_final(
     &                  span, mxvl, stress_n, stress_np1, e_vec_np1,
     &                  nu_vec_np1, g_vec_np1, deps, instat,
     &                  history_np1, trace_eps_np1, debug, iout )
c

c
c
c                        plastic work update
c                        -------------------
c
      if ( iter .gt. 0 )
     &     call mm05_plastic_work( iout, span, mxvl, stress_n,
     &      stress_np1, yield, deps, nu_vec_np1, e_vec_np1, g_vec_np1 )
c
c
c
       return
c
       end

c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_gp_step1_set  -- last modified 1/23/04  by kbc      *
c *                                                                 *
c *        for step 1, initialize the stress and history vectors    *
c *                                                                 *
c *******************************************************************
c
      subroutine mm05_gp_step1_set( stress_n, history_n, mm_props,
     &                      sigyld_vec_n, span, mxvl )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &   mxvl, span
c
      double precision
     & mm_props(mxvl,5), sigyld_vec_n(mxvl),  stress_n(mxvl,*),
     &  history_n(span,*)
c
c                  local parameters
c                  ----------------
      integer i, iword(2)
c
      double precision
     &  root2third, zero, dword
c
      equivalence( dword, iword )
c
      data root2third, zero / 0.81649658, 0.0 /
c
      iword(1) = 3
      iword(2) = 0
c
c                 no longer zero _n for step 1. Could be pre-loaded
c                 w/ user-defined initial stresses
c
        do i = 1, span
            stress_n(i,7) = zero
            stress_n(i,8) = zero
            stress_n(i,9) = zero
c
            history_n(i,1) = zero
            history_n(i,2) = sigyld_vec_n(i)*root2third
            history_n(i,3) = zero
            history_n(i,4) = dword
            history_n(i,5) = zero
            history_n(i,6) = zero
            history_n(i,7) = zero
            history_n(i,8) = zero
            history_n(i,9) = zero
            history_n(i,10) = zero
            history_n(i,11) = zero
            history_n(i,12) = zero
            history_n(i,13) = zero
            history_n(i,14) = zero
         end do
c
      return
c
      end


c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_gp_init  -- last modified  7/11/11 by jcs           *
c *                                                                 *
c *        compute deviatoric stress, trial elastic state,          *
c *        and evaluate material state at the beginning of each     *
c *        step                                                     *
c *                                                                 *
c *******************************************************************
c
       subroutine mm05_gp_init(span, mxvl, iout, deps, iostat, felem,
     &   prior_linear, history_n, stress_n, signal_flag,
     &   g_vec_np1, g_vec_n, alpha_n, trial_elas_stress_np1, debug,
     &   gpn, sigyld_vec_np1, sigyld_vec_n, instat, yield, iter,
     &   trace_eps_np1, step, h_gp_np1, h_gp_n, tau,
     &   e_vec_np1, e_vec_n, nu_vec_np1, nu_vec_n )
      implicit none
c
c               parameter declarations
c               ----------------------
      integer span, mxvl, iout, iostat(*), instat(*), iter, step,
     &              gpn, felem
      logical debug, prior_linear(*), yield(*), signal_flag
      double precision
     & deps(mxvl, *), history_n(span,*), stress_n(mxvl,*),
     & g_vec_n(*), g_vec_np1(*), alpha_n( mxvl, *),
     & trial_elas_stress_np1(mxvl, *), sigyld_vec_np1(*),
     & sigyld_vec_n(*), trace_eps_np1(*),
     & e_vec_n(*), e_vec_np1(*), nu_vec_n(*), nu_vec_np1(*),
     & h_gp_n(*), h_gp_np1(*), tau(*)
c
c               local parameters
c               -----------------
      integer iword(2), i, j, l
      logical local_debug
      double precision
     & dword, deps_mean, de1, de2, de3, de4, de5, de6, trace_deps,
     & een1, een2, een3, een4, een5, een6, e1, e2, e3, e4, e5, e6,
     & e_n, e_n1, nu_n, nu_n1, g_n, g_n1, hk_n, hk_n1, hi_n1, lk,
     & f_n, eps_n, k_n, kb_n1, f_n1, F_trial, rtse(6),
     & zero, one, two, three, root2third, twothird, yld_tol,
     & bar_alpha_n(6), rsen(6), norm_rtse, norm_rsen,
     & eps_mean, norm_deps, stress_n_mean
c
       equivalence ( dword, iword )
       data zero, one, two, three, twothird, yld_tol, root2third
     & / 0.0, 1.0, 2.0, 3.0, 0.66666666, 0.0000001, 0.8164965809 /
       local_debug = .false.
c
       do i = 1, span
          dword  = history_n(i, 4)
          iostat(i) = iword(1)
          prior_linear(i) = iostat(i) .eq. 3
c
c      deviatoric strain components
c
          trace_deps = deps(i,1) + deps(i,2) + deps(i,3)
          deps_mean  = trace_deps / three
          de1        = deps(i,1) - deps_mean
          de2        = deps(i,2) - deps_mean
          de3        = deps(i,3) - deps_mean
          de4        = deps(i,4)
          de5        = deps(i,5)
          de6        = deps(i,6)
c
c     elastic components of total strain at start of step
c
          e_n  = e_vec_n(i)
          nu_n = nu_vec_n(i)
          g_n  = e_n/(two*(one+nu_n))
          een1 = (stress_n(i,1)-nu_n*(stress_n(i,2)
     &                              + stress_n(i,3)))/e_n
          een2 = (stress_n(i,2)-nu_n*(stress_n(i,1)
     &                              + stress_n(i,3)))/e_n
          een3 = (stress_n(i,3)-nu_n*(stress_n(i,1)
     &                              + stress_n(i,2)))/e_n
          een4 = stress_n(i,4) / g_n
          een5 = stress_n(i,5) / g_n
          een6 = stress_n(i,6) / g_n
c
          g_vec_n(i) = g_n
c
c     deviatoric elastic strain components at start of step plus
c     deviatoric strain increment over step
c
          trace_eps_np1(i) = een1 + een2 + een3 + trace_deps
          eps_mean         = (een1 + een2 + een3 ) / three
          e1               = (een1 - eps_mean) + de1
          e2               = (een2 - eps_mean) + de2
          e3               = (een3 - eps_mean) + de3
          e4               = een4 + de4
          e5               = een5 + de5
          e6               = een6 + de6
c
c    deviatoric components of trial elastic stress state
c
          e_n1  = e_vec_np1(i)
          nu_n1 = nu_vec_np1(i)
          g_n1  = e_n1/(two*(one+nu_n1))
          trial_elas_stress_np1(i,1) = two*g_n1*e1
          trial_elas_stress_np1(i,2) = two*g_n1*e2
          trial_elas_stress_np1(i,3) = two*g_n1*e3
          trial_elas_stress_np1(i,4) = g_n1*e4
          trial_elas_stress_np1(i,5) = g_n1*e5
          trial_elas_stress_np1(i,6) = g_n1*e6
c
          g_vec_np1(i) = g_n1
c
c     pull out backstress and other parameters at start of step
c
          alpha_n(i,1) = history_n(i,6)
          alpha_n(i,2) = history_n(i,7)
          alpha_n(i,3) = history_n(i,8)
          alpha_n(i,4) = history_n(i,9)
          alpha_n(i,5) = history_n(i,10)
          alpha_n(i,6) = history_n(i,11)
c
c     determine trial backstress value for an elastic step
c
          hk_n  = (one-tau(i))*twothird*h_gp_n(i)
          hk_n1 = (one-tau(i))*twothird*h_gp_np1(i)
          lk    = one
          if ( abs( hk_n ) .gt. yld_tol ) lk = hk_n1/hk_n
c
          bar_alpha_n(1) = lk*alpha_n(i,1)
          bar_alpha_n(2) = lk*alpha_n(i,2)
          bar_alpha_n(3) = lk*alpha_n(i,3)
          bar_alpha_n(4) = lk*alpha_n(i,4)
          bar_alpha_n(5) = lk*alpha_n(i,5)
          bar_alpha_n(6) = lk*alpha_n(i,6)
c
c     calculate relative trial stress ( deviatoric )
c
          rtse(1) = trial_elas_stress_np1(i,1) - bar_alpha_n(1)
          rtse(2) = trial_elas_stress_np1(i,2) - bar_alpha_n(2)
          rtse(3) = trial_elas_stress_np1(i,3) - bar_alpha_n(3)
          rtse(4) = trial_elas_stress_np1(i,4) - bar_alpha_n(4)
          rtse(5) = trial_elas_stress_np1(i,5) - bar_alpha_n(5)
          rtse(6) = trial_elas_stress_np1(i,6) - bar_alpha_n(6)
c
          norm_rtse = sqrt( rtse(1)**2 + rtse(2)**2 + rtse(3)**2 +
     &         two*( rtse(4)**2 + rtse(5)**2 + rtse(6)**2 ) )

c
c     calculate relative deviatoric stress at step n
c     from total stress at n and backstress
c
          stress_n_mean = ( stress_n(i,1) + stress_n(i,2) +
     &                     stress_n(i,3) ) / three
          rsen(1) = stress_n(i,1) - stress_n_mean - alpha_n(i,1)
          rsen(2) = stress_n(i,2) - stress_n_mean - alpha_n(i,2)
          rsen(3) = stress_n(i,3) - stress_n_mean - alpha_n(i,3)
          rsen(4) = stress_n(i,4) - alpha_n(i,4)
          rsen(5) = stress_n(i,5) - alpha_n(i,5)
          rsen(6) = stress_n(i,6) - alpha_n(i,6)
c
          norm_rsen = sqrt( rsen(1)**2 + rsen(2)**2 +
     &         rsen(3)**2 + two*( rsen(4)**2 + rsen(5)**2 +
     &         rsen(6)**2 ))
c
c     estimate the values for the yield function of the current step
c     and the real value for the yield function of the previous step
c
          k_n   = history_n(i,2)
          hi_n1 = tau(i)*( twothird*h_gp_np1(i) )
          eps_n = history_n(i,3)
          kb_n1 = root2third*sigyld_vec_np1(i) + hi_n1*eps_n
c
          f_n   = norm_rsen - k_n
          f_n1  = norm_rtse - kb_n1
          F_trial = f_n1 - f_n
c
          if ( local_debug ) then
             write(iout,9000) felem+i-1, gpn, step
             write(iout,9010) f_n, norm_rsen, k_n
             write(iout,9020) f_n1, norm_rtse, kb_n1
             write(iout,9030) f_n1, F_trial, yld_tol*k_n
          end if
c
c
c     set various flags to indicate yielding
c         state variable = history_n(i,4)
c                        = 1, if gauss point is yield
c                        = 3, if not yielding
c
c     N.B. iteration zero is a pseudo iteration used to compute
c          stresses due to imposed displacement and temperature
c          changes.  yielding from a previously linear state is
c          not allowed for this iteration.
c
          instat(i) = 1
          yield(i) = .true.
c
          if ( ( f_n1 .lt. yld_tol*k_n ) .or.
     &         ( F_trial .lt. yld_tol*k_n ) ) then
             yield(i) = .false.
             instat(i) = 3
          end if
c
          if ( local_debug ) then
             if( instat(i) .eq. 3 ) then
                write(iout, *) 'state is elastic'
c                write(iout, *) 'F_trial w/o h ', F_trial
c                write(iout, 9040)  (rtse(l), l=1,6)
c                write(iout, 9040)  (rsen(l), l=1,6)
             else
                write(iout, *) 'state is plastic'
             end if
          end if
c
          if ( iter .gt. 0 ) then
             if(signal_flag) then
                if( history_n(i,3) .eq. zero .and. yield(i)) then
                   write(iout, 9050) felem+i-1, gpn
                end if
             end if
             cycle
          end if
c
          if ( prior_linear(i) )then
             yield(i) = .false.
             instat(i) = 3
             cycle
          end if

       end do
c

       return
c
 9000  format('Debug mm05_gp_init:  element ', i7, ' at gp ', i3,
     &   ' step number ', i6 )
 9010  format( 'f_n,  norm_rsen, kn    =', 3e15.6 )
 9020  format( 'f_n1, norm_rtse, kb_n1 =', 3e15.6 )
 9030  format( 'f_n1, F_trial, tol*k_n =', 3e15.6 )
 9040  format( 3e15.6, /, 3e15.6 )
 9050  format( 'element ', i7, ' at gauss point ', i3,
     &      '   begins yielding' )
       end

c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_gp_compute  -- last modified 7/05/11 by jcs         *
c *                                                                 *
c *        currently includes:                                      *
c *           !*automatic subincrementation revision                *
c *           !*subincrementation of strain increment               *
c *                                                                 *
c *******************************************************************
c
       subroutine mm05_gp_compute(
     &   span, mxvl, debug, iout, history_np1, history_n,
     &   mm_props, sigyld_vec_np1, sigyld_vec_n, yield,
     &   iter, stress_np1, stress_n, adaptive_possible,
     &   cut_step_size_now, signal_flag, g_vec_np1, g_vec_n,
     &   h_gp_np1, h_gp_n, delta_gp_np1, delta_gp_n,
     &   beta_gp_np1, beta_gp_n, tau, gpn, felem, step, deps,
     &   prior_linear, trial_stress )
       implicit none
c
c                    parameter declarations
c                    ----------------------
       integer span, mxvl, iout, gpn, felem, step, iter
       logical yield(*), debug, adaptive_possible, cut_step_size_now,
     &         signal_flag, prior_linear(*)
      double precision
     & history_n(span,*), history_np1(span,*), mm_props(mxvl,*),
     & sigyld_vec_np1(*), sigyld_vec_n(*), stress_np1(mxvl,*),
     & stress_n(mxvl,*), deps(mxvl,*), g_vec_np1(*), g_vec_n(*),
     & trial_stress(mxvl,*), h_gp_n(mxvl), h_gp_np1(mxvl),
     & delta_gp_np1(mxvl), delta_gp_n(mxvl),
     & beta_gp_n(mxvl), beta_gp_np1(mxvl), tau(mxvl)
c
c                    local parameters
c                    ----------------
c
       integer i, j, s, nsubinc, mxsubinc, l
       logical local_debug, onestep
      double precision
     & delta_n, Hi_n, Hk_n, beta_n, g_n, k0_n, k_n,
     & delta_n1, Hi_n1, Hk_n1, beta_n1, g_n1, k0_n1,
     & g_s, Hi_s, Hk_s, beta_s, delta_s, k0_s,
     & sig_tol, rel_tol, eta, error, a, b, twothird,
     & mean_stress, norm_equiv_sm1, one, two, three, root2third,
     & lambda_1step, lambda_2step, lambda_nstep, mean_deps,
     & eps_n, eps_np1, eps_1step, eps_2step, eps_nstep, dHis,
     & alpha_n(6), alpha_np1(6), alpha_1step(6), alpha_2step(6),
     & alpha_nstep(6), dev_stress_n(6), dev_stress_np1(6),
     & stress_1step(6), stress_2step(6), stress_nstep(6),
     & diff(6), devdeps(6), sig_t(6), equiv(6),
     & norm_equiv, lk, equiv_nr(6), norm_equiv_nr
c
       data  one, two, three, twothird, root2third
     &   /  1.0, 2.0, 3.0, 0.66666666, 0.81649658  /
       data mxsubinc / 10 /
       local_debug = .false.
c       onestep = .false.
c
       do i = 1, span
c
          if( .not. yield(i) ) cycle
c
c     set material vectors to scalar values and scale them from uniaxial
c     values to 3d values
c
          g_n      = g_vec_n(i)
          g_n1     = g_vec_np1(i)
          Hi_n     = tau(i)*twothird*h_gp_n(i)
          Hi_n1    = tau(i)*twothird*h_gp_np1(i)
          Hk_n     = (one-tau(i))*twothird*h_gp_n(i)
          Hk_n1    = (one-tau(i))*twothird*h_gp_np1(i)
          beta_n   = root2third*beta_gp_n(i)
          beta_n1  = root2third*beta_gp_np1(i)
          delta_n  = twothird*delta_gp_n(i)
          delta_n1 = twothird*delta_gp_np1(i)
          k0_n     = root2third*sigyld_vec_n(i)
          k0_n1    = root2third*sigyld_vec_np1(i)
c
          sig_tol  = one*mm_props(i,6)
          k_n      = history_n(i,2)
          eps_n    = history_n(i,3)
c
          rel_tol = sig_tol * k_n
c
          mean_stress = ( stress_n(i, 1) + stress_n(i,2) +
     &         stress_n(i, 3) )/three
          dev_stress_n(1) = stress_n(i,1) - mean_stress
          dev_stress_n(2) = stress_n(i,2) - mean_stress
          dev_stress_n(3) = stress_n(i,3) - mean_stress
          dev_stress_n(4) = stress_n(i,4)
          dev_stress_n(5) = stress_n(i,5)
          dev_stress_n(6) = stress_n(i,6)
c
          alpha_n(1:6) = history_n(i,6:11)
c
          mean_deps = ( deps(i, 1) + deps(i,2) +
     &         deps(i, 3) )/three
c
          devdeps(1) = deps(i,1) - mean_deps
          devdeps(2) = deps(i,2) - mean_deps
          devdeps(3) = deps(i,3) - mean_deps
          devdeps(4) = deps(i,4)
          devdeps(5) = deps(i,5)
          devdeps(6) = deps(i,6)
c
c     compute the fraction of the step that is purely elastic and
c     update the stress_s vector to include the additional elastic
c     components.  also update devdeps.
c
          call mm05_gp_elastic_fraction( dev_stress_n, alpha_n,
     &         devdeps, eta, k_n, iout, felem, gpn, i, step,
     &         g_n, g_n1, g_s, Hi_n, Hi_n1, Hi_s, Hk_n, Hk_n1, Hk_s,
     &         beta_n, beta_n1, beta_s, delta_n, delta_n1, delta_s,
     &         k0_n, k0_n1, k0_s, eps_n, adaptive_possible,
     &         cut_step_size_now, onestep )
          if ( cut_step_size_now ) return
c
c     calculate the number of subincrements
c     1) determine stresses for 1 subincrement
c
          call mm05_gp_nsteps(1, devdeps, dev_stress_n, alpha_n,
     &         eps_n, g_s, Hi_s, Hk_s, beta_s, delta_s, k0_s,
     &         g_n1, Hi_n1, Hk_n1, beta_n1, delta_n1, k0_n1,
     &         stress_1step, alpha_1step, eps_1step, lambda_1step, dHis,
     &         sig_t, norm_equiv_sm1, iout, gpn, felem, i, step,
     &         iter, signal_flag, adaptive_possible, cut_step_size_now)
          if ( cut_step_size_now ) return
c          if ( onestep ) goto 1000
c
c     2) determine stresses for 2 subincrements
c
          call mm05_gp_nsteps(2, devdeps, dev_stress_n, alpha_n, eps_n,
     &         g_s, Hi_s, Hk_s, beta_s, delta_s, k0_s,
     &         g_n1, Hi_n1, Hk_n1, beta_n1, delta_n1, k0_n1,
     &         stress_2step, alpha_2step, eps_2step, lambda_2step, dHis,
     &         sig_t, norm_equiv_sm1, iout, gpn, felem, i, step,
     &         iter, signal_flag, adaptive_possible, cut_step_size_now)
          if ( cut_step_size_now ) return
c
c     3) compute error estimate
c
          diff(1) = abs(stress_2step(1) - stress_1step(1))
          diff(2) = abs(stress_2step(2) - stress_1step(2))
          diff(3) = abs(stress_2step(3) - stress_1step(3))
          diff(4) = abs(stress_2step(4) - stress_1step(4))
          diff(5) = abs(stress_2step(5) - stress_1step(5))
          diff(6) = abs(stress_2step(6) - stress_1step(6))
          error = maxval(diff)
c
c     4) determine number of subincrements needed to drive error
c        below some tolerance
c
          nsubinc = ceiling(1.1*error/rel_tol)
          if (nsubinc .lt. two) nsubinc = two
          if(local_debug) then
             write(iout,9000) felem+i-1, gpn, step
             write(iout,9001) error, nsubinc
             write(iout,9002) (stress_1step(l), l=1,6)
             write(iout,9003) (stress_2step(l), l=1,6)
          end if
c
c     5) limit number of subincrements to some reasonable number
c
          if (nsubinc .gt. mxsubinc) then
             if (signal_flag) write(iout,9007) nsubinc, mxsubinc
             nsubinc = mxsubinc
          end if
c
          history_np1(i,12) = nsubinc*one
          if(signal_flag) then
             write(iout,9006) felem+i-1, gpn, step, nsubinc
          end if
c
c     if the number of subincrements is less then 3, extrapolate
c
          if( nsubinc .lt. three ) then
c
             dev_stress_np1(1) = two*stress_2step(1) - stress_1step(1)
             dev_stress_np1(2) = two*stress_2step(2) - stress_1step(2)
             dev_stress_np1(3) = two*stress_2step(3) - stress_1step(3)
             dev_stress_np1(4) = two*stress_2step(4) - stress_1step(4)
             dev_stress_np1(5) = two*stress_2step(5) - stress_1step(5)
             dev_stress_np1(6) = two*stress_2step(6) - stress_1step(6)
             alpha_np1(1) = two*alpha_2step(1) - alpha_1step(1)
             alpha_np1(2) = two*alpha_2step(2) - alpha_1step(2)
             alpha_np1(3) = two*alpha_2step(3) - alpha_1step(3)
             alpha_np1(4) = two*alpha_2step(4) - alpha_1step(4)
             alpha_np1(5) = two*alpha_2step(5) - alpha_1step(5)
             alpha_np1(6) = two*alpha_2step(6) - alpha_1step(6)
             eps_np1 = two*eps_2step - eps_1step
             lambda_nstep = lambda_2step
c
c     o.w., compute stresses using nsubinc subincrements and then extrapolate
c
          else
c
             call mm05_gp_nsteps(
     &            nsubinc, devdeps, dev_stress_n, alpha_n, eps_n,
     &            g_s, Hi_s, Hk_s, beta_s, delta_s, k0_s,
     &            g_n1, Hi_n1, Hk_n1, beta_n1, delta_n1, k0_n1,
     &            stress_nstep, alpha_nstep, eps_nstep, lambda_nstep,
     &            dHis, sig_t, norm_equiv_sm1,
     &            iout, gpn, felem, i, step, iter, signal_flag,
     &            adaptive_possible, cut_step_size_now)
c
             if ( cut_step_size_now ) return
c
             a = (one*nsubinc)/(nsubinc-two)
             b = (one*two)/(nsubinc-two)
             dev_stress_np1(1) = a*stress_nstep(1) - b*stress_2step(1)
             dev_stress_np1(2) = a*stress_nstep(2) - b*stress_2step(2)
             dev_stress_np1(3) = a*stress_nstep(3) - b*stress_2step(3)
             dev_stress_np1(4) = a*stress_nstep(4) - b*stress_2step(4)
             dev_stress_np1(5) = a*stress_nstep(5) - b*stress_2step(5)
             dev_stress_np1(6) = a*stress_nstep(6) - b*stress_2step(6)
             alpha_np1(1) = a*alpha_nstep(1) - b*alpha_2step(1)
             alpha_np1(2) = a*alpha_nstep(2) - b*alpha_2step(2)
             alpha_np1(3) = a*alpha_nstep(3) - b*alpha_2step(3)
             alpha_np1(4) = a*alpha_nstep(4) - b*alpha_2step(4)
             alpha_np1(5) = a*alpha_nstep(5) - b*alpha_2step(5)
             alpha_np1(6) = a*alpha_nstep(6) - b*alpha_2step(6)
             eps_np1  = a*eps_nstep - b*eps_2step
c
             if(local_debug) then
                write(iout,9004) (stress_nstep(l), l=1,6)
             end if
c
          end if
c
          if(local_debug) then
             write(iout,9005) (dev_stress_np1(l), l=1,6)
          end if
c
c     compute the deviatoric part of the updated stress
c     and store in the stress vector - add the mean stress
c     later
c
c
c 1000     if ( onestep ) then
c             dev_stress_np1(1:6) = stress_1step(1:6)
c             alpha_np1(1:6) = alpha_1step(1:6)
c             eps_np1 = eps_1step
c             lambda_nstep = lambda_1step
c          end if
c
          stress_np1(i,1:6) = dev_stress_np1(1:6)
c
c     update history
c
          history_np1(i,1) = lambda_nstep
          history_np1(i,2) = k0_n1+Hi_n1*eps_np1
          history_np1(i,3) = eps_np1
          history_np1(i,5) = norm_equiv_sm1
          history_np1(i,6:11) = alpha_np1(1:6)
          history_np1(i,13) = k_n
          history_np1(i,14) = dHis
c
          trial_stress(i,1:6) = sig_t(1:6)
c
       end do
c
       return
c
 9000  format('Debug mm05_gp_compute:  element ', i7, ' at gp ', i3,
     &   ' step number ', i7 )
 9001  format('computed error: ', e15.6, '  nsubinc: ', i6)
 9002  format('stress_1step: ', / 3e15.6 / 3e15.6)
 9003  format('stress_2step: ', / 3e15.6 / 3e15.6)
 9004  format('stress_nstep: ', / 3e15.6 / 3e15.6)
 9005  format('extrapolated stresses: ', / 3e15.6 / 3e15.6)
 9006  format('element ', i7, ' at gp ', i3, ' step ', i6,
     &        ' used ', i3, ' substeps')
 9007  format('>>>>> warning: error tolerance for adaptive integration',
     &        'leads to ', i5, ' subincrements. exceeds limit.',/15x,
     &        'setting number of subincrements to ', i3, /)
c
       end
c
c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_gp_nsteps  -- last modified 7/05/11 by jcs          *
c *                                                                 *
c *                                                                 *
c *******************************************************************
c
      subroutine mm05_gp_nsteps(
     &     nsubinc, devdeps, dev_stress_s, alpha_s, eps_s,
     &     g_s, Hi_s, Hk_s, beta_s, delta_s, k0_s,
     &     g_n1, Hi_n1, Hk_n1, beta_n1, delta_n1, k0_n1,
     &     stress_j, alpha_j, eps_j, lambda, dHi,
     &     sig_t, norm_equiv_j, iout, gpn, felem,
     &     i, step, iter, signal_flag, adaptive_possible,
     &     cut_step_size_now)

      implicit none
c
c     parameter declarations
c     ----------------------
c
      integer iout, gpn, felem, step, iter, i, nsubinc
      logical signal_flag, adaptive_possible, cut_step_size_now
      double precision
     &     devdeps(*), dev_stress_s(*), alpha_s(*), eps_s,
     &     stress_j(*), alpha_j(*), eps_j, lambda, dHi,
     &     sig_t(*),  norm_equiv_j, g_s, Hi_s, Hk_s, beta_s, delta_s,
     &     k0_s, g_n1, Hi_n1, Hk_n1, beta_n1, delta_n1, k0_n1
c
c    local parameters
c    ----------------
c
      integer j, l
      logical local_debug
      double precision
     &     one, equiv_t(6), norm_equiv_t, k_j, two, d, root2third,
     &     dHk, dk0, dg, dbeta, ddel, lk, mk, zero, tol, dnsubinc,
     &     alpha_j1(6), Hk_j, Hk_j1, Hi_j, Hi_j1, g_j, g_j1,
     &     beta_j, beta_j1, delta_j, delta_j1, k0_j, k0_j1, alp_t(6),
     &     eps_j1, stress_j1(6)
c
      data  zero, one, two, root2third, tol
     &     /   0.0, 1.0, 2.0, 0.81649658, 0.00000001  /
      local_debug = .false.
c
c     determine size of subincrement and the associated change in the
c     material parameters per subincrement
c
      dnsubinc = dble( nsubinc )
      d = one/dnsubinc
      if (local_debug) write(iout,*) nsubinc, dnsubinc, d
c
      dHk   = d*(Hk_n1 - Hk_s)
      dHi   = d*(Hi_n1 - Hi_s)
      dg    = d*(g_n1 - g_s)
      dk0   = d*(k0_n1 - k0_s)
      dbeta = d*(beta_n1 - beta_s)
      ddel  = d*(delta_n1 - delta_s)
c
c     initialize the material properties for the 0'th subincrement
c
      Hk_j    = Hk_s
      Hi_j    = Hi_s
      g_j     = g_s
      k0_j    = k0_s
      beta_j  = beta_s
      delta_j = delta_s
c
c     initialize subincrementation vectors and scalars
c
      stress_j(1:6) = dev_stress_s(1:6)
      alpha_j(1:6)  = alpha_s(1:6)
      eps_j         = eps_s
c
c     perform subincrementation loop
c
      do j = 1, nsubinc
c
c     increment material properties from substep j --> j+1
c
         Hk_j1    = Hk_j + dHk
         Hi_j1    = Hi_j + dHi
         g_j1     = g_j + dg
         k0_j1    = k0_j + dk0
         beta_j1  = beta_j + dbeta
         delta_j1 = delta_j + ddel
c
c     compute temperature dependent material ratios to adjust stress
c     and backstress tensors at "j" to their updated values at "j+1"
c     also determine the yield surface size at the previous subincrement
c
         mk = g_j1/g_j
         lk = one
         if ( abs( Hk_j ) .gt. tol ) lk = Hk_j1/Hk_j
         k_j = k0_j + Hi_j*eps_j
c
c     compute a trial stress state (deviatoric) and trial backstress
c
         sig_t(1) = mk*stress_j(1) + two*g_j1*devdeps(1)*d
         sig_t(2) = mk*stress_j(2) + two*g_j1*devdeps(2)*d
         sig_t(3) = mk*stress_j(3) + two*g_j1*devdeps(3)*d
         sig_t(4) = mk*stress_j(4) +     g_j1*devdeps(4)*d
         sig_t(5) = mk*stress_j(5) +     g_j1*devdeps(5)*d
         sig_t(6) = mk*stress_j(6) +     g_j1*devdeps(6)*d
c
         alp_t(1:6) = lk*alpha_j(1:6)
c
c     determine relative trial stress w/ thermal effects on alpha
c     and its norm
c
         equiv_t(1:6) = sig_t(1:6) - alp_t(1:6)
         norm_equiv_t = sqrt(
     &        equiv_t(1)**2 + equiv_t(2)**2 + equiv_t(3)**2 +
     &   two*(equiv_t(4)**2 + equiv_t(5)**2 + equiv_t(6)**2))
c
c     determine the yield surface radius at the substep j
c
         norm_equiv_j = sqrt(
     &            (stress_j(1)-alpha_j(1))**2 +
     &            (stress_j(2)-alpha_j(2))**2 +
     &            (stress_j(3)-alpha_j(3))**2 +
     &        two*(stress_j(4)-alpha_j(4))**2 +
     &        two*(stress_j(5)-alpha_j(5))**2 +
     &        two*(stress_j(6)-alpha_j(6))**2 )
c
         if ( local_debug ) then
            write(iout,9000) felem+i-1, gpn, step, j, nsubinc
            write(iout,9001) (stress_j(l), l=1,6)
            write(iout,9002) (alpha_j(l), l=1,6)
            write(iout,9007) g_j,Hk_j,Hi_j,k0_j,beta_j,delta_j
            write(iout,9008) g_j1,Hk_j1,Hi_j1,k0_j1,beta_j1,delta_j1
            write(iout,9011) (sig_t(l), l=1,6)
            write(iout,9012) (alp_t(l), l=1,6)
            write(iout,9013) norm_equiv_t, norm_equiv_j
         end if
c
c     determine consistency parameter lambda
c

         call mm05_gp_lambda_solve(
     &        iout, norm_equiv_t, norm_equiv_j, k_j, g_j1,
     &        delta_j1, beta_j1, Hk_j1, Hi_j1, dHi, k0_j1,
     &        adaptive_possible, cut_step_size_now,
     &        signal_flag, gpn, felem, step, i, iter, eps_j, lambda )
         if ( cut_step_size_now ) return
c
c     update subincrementation values of stress, backstress, & eps
c
         eps_j1 = eps_j+lambda
c
         stress_j1(1) = sig_t(1)-2*g_j1*lambda*equiv_t(1)/norm_equiv_t
         stress_j1(2) = sig_t(2)-2*g_j1*lambda*equiv_t(2)/norm_equiv_t
         stress_j1(3) = sig_t(3)-2*g_j1*lambda*equiv_t(3)/norm_equiv_t
         stress_j1(4) = sig_t(4)-2*g_j1*lambda*equiv_t(4)/norm_equiv_t
         stress_j1(5) = sig_t(5)-2*g_j1*lambda*equiv_t(5)/norm_equiv_t
         stress_j1(6) = sig_t(6)-2*g_j1*lambda*equiv_t(6)/norm_equiv_t
c
         alpha_j1(1) = alp_t(1) +
     &        Hk_j1*lambda*equiv_t(1)/norm_equiv_t
         alpha_j1(2) = alp_t(2) +
     &        Hk_j1*lambda*equiv_t(2)/norm_equiv_t
         alpha_j1(3) = alp_t(3) +
     &        Hk_j1*lambda*equiv_t(3)/norm_equiv_t
         alpha_j1(4) = alp_t(4) +
     &        Hk_j1*lambda*equiv_t(4)/norm_equiv_t
         alpha_j1(5) = alp_t(5) +
     &        Hk_j1*lambda*equiv_t(5)/norm_equiv_t
         alpha_j1(6) = alp_t(6) +
     &        Hk_j1*lambda*equiv_t(6)/norm_equiv_t
c
         if ( local_debug ) then
            write(iout,9021) (stress_j1(l), l=1,6)
            write(iout,9022) (alpha_j1(l), l=1,6)
         end if
c
c     increment j+1 -> j to initialize for next subincrement
c
         alpha_j(1:6) = alpha_j1(1:6)
         stress_j(1:6) = stress_j1(1:6)
         eps_j = eps_j1
c
         Hk_j    = Hk_j1
         Hi_j    = Hi_j1
         g_j     = g_j1
         k0_j    = k0_j1
         beta_j  = beta_j1
         delta_j = delta_j1
c
      end do
c
 9000 format('Debug mm05_gp_nsteps:  element ', i7, ' at gp ', i3,
     &   ' step number ', i6, ' subincrement', i3, ' of', i3 )
 9001 format('stress_j : ', / 3e15.6 / 3e15.6)
 9002 format('alpha_j  : ', / 3e15.6 / 3e15.6)
 9007 format('g_j, Hk_j, Hi_j:         ', 3f15.6 /
     &       'k0_j, beta_j, delta_j:   ', 3f15.6 )
 9008 format('g_j1, Hk_j1, Hi_j1:      ', 3f15.6 /
     &       'k0_j1, beta_j1, delta_j1:', 3f15.6 )
 9011 format('sig_t    : ', / 3e15.6 / 3e15.6)
 9012 format('alp_t    " ', / 3e15.6 / 3e15.6)
 9013 format('norm_equiv_t, norm_equiv_j:', 2e15.6)
 9021 format('stress_j1: ', / 3e15.6 / 3e15.6)
 9022 format('alpha_j1 : ', / 3e15.6 / 3e15.6)
       end
c
c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_gp_lambda_solve                                     *
c *                 -- last modified 7/05/11 by jcs                 *
c *                                                                 *
c *                                                                 *
c *******************************************************************
c
      subroutine mm05_gp_lambda_solve(
     &     iout, norm_equiv_t, norm_equiv_j, k_j, g_n1,
     &     delta_n1, beta_n1, Hk_n1, Hi_n1, dHi, k0_n1,
     &     adaptive_possible, cut_step_size_now,
     &     signal_flag, gpn, felem, step, i, iter, eps_j, lambda )

      implicit none
c
c     parameter declarations
c     ----------------------
c
      integer iout, gpn, felem, step, i, iter, sub
      logical signal_flag, adaptive_possible, cut_step_size_now
      double precision
     &     norm_equiv_t, norm_equiv_j, k_j, g_n1, delta_n1,
     &     beta_n1, Hk_n1, Hi_n1, dHi, k0_n1, eps_j, lambda
c
c
c     local parameters
c     ----------------
      double precision
     &     a1, a2, a3, a4, a5, g1, H_n1, f_j, kb_j1, a, b, c,
     &     zero, half, one, two, temp, l1, l2
c
      data zero, half, one, two / 0.0, 0.5, 1.0, 2.0 /
c
c     initialize variables and constants prior to starting loop
c
c
      H_n1   = Hi_n1 + Hk_n1
      f_j    = norm_equiv_j - k_j
      kb_j1  = k0_n1 + Hi_n1*eps_j
c
c     determine quantities needed to compute quadratic coefficients
c
      g1 = half*( two*g_n1 + H_n1 )
      a1 = norm_equiv_t - kb_j1
      a2 = a1 - f_j
      a3 = delta_n1 - ( two*g_n1 + dHi )
      a4 = beta_n1*( delta_n1 + H_n1 )
c
c     calculate quadratic coefficients: a*lam^2 + b*lam + c = 0
c
      a = two*g1*a3
      b = two*g1*a2 + a4 - a1*a3
      c = - a1*a2
c
c     compute discriminant, check that it is positive
c
      temp = b**2 -4*a*c
      if( temp .lt. zero) then
         write(iout,7000) felem+i-1, step, gpn
         if ( adaptive_possible ) then
            cut_step_size_now = .true.
            return
         end if
         write(iout,7010)
         write(iout,7001) a, b, c
         call die_abort
      end if
c
c     compute roots of the quadratic equation, check that at least
c     one of the roots is positive
c
      l1 = half/a*(-b-sqrt(temp))
      l2 = half/a*(-b+sqrt(temp))
c
      if(l1 .lt. zero .and. l2 .lt. zero) then
         write(iout,7002) felem+i-1, step, gpn
         write(iout,7011) H_n1, f_j, k_j, kb_j1, a1, a2, a3, a4,
     &     g1, a, b, c, l1, l2, norm_equiv_t, norm_equiv_j
         if ( adaptive_possible ) then
            cut_step_size_now = .true.
            return
         end if
         write(iout,7010)
         call die_abort
      end if
c
c     set lambda as the smallest positive root
c
      if(l1 .gt. zero) then
         lambda = l1
         if(lambda .gt. l2 .and. l2 .gt. zero) then
            lambda = l2
         end if
      else
         lambda = l2
      end if
c
c     final check to ensure that lambda solves the quadratic equation
c     within some tolerance
c
      temp  = a*lambda**2 + b*lambda + c
      temp  = abs(temp)
      if( temp > 0.000001) then
         write(iout,7003) felem+i-1, step, gpn, lambda
         write(iout,7011) H_n1, f_j, k_j, kb_j1, a1, a2, a3, a4,
     &     g1, a, b, c, l1, l2, norm_equiv_t, norm_equiv_j
         write(iout,7012) temp
         if ( adaptive_possible ) then
            cut_step_size_now = .true.
            return
         end if
         write(iout,7010)
         call die_abort
      end if
c
 7000 format('>> ERROR: no real value for lambda in generalized',
     &     ' plasticity algorithm',/,' for element ', i7, ' at step '
     &     , i6,' gauss point ', i2)
 7001 format('(a,b,c) = ', 3(1x,e13.6), / )
 7002 format('>> ERROR: no positive value for lambda in',
     &     ' generalized plasticity algorithm',/,' for element ', i7,
     &     ' at step ', i6, ' gauss point ', i2)
 7003 format('>> ERROR: unknown error with quadratic solution in',
     &     ' generalized plasticity algorithm',/,' for element ', i7,
     &     ' at step ', i6, ' gauss point ', i2, ' lambda ', e13.6)
 7010 format( '>> Fatal Error: cannot requst step size reduction',
     &        ' Job terminated.' // )
 7011 format('H_n1,   f_j,   k_j, kb_j1', 4(1x,e13.6), /,
     &       '  a1,    a2,    a3,    a4', 4(1x,e13.6), /,
     &       '  g1,     a,     b,     c', 4(1x,e13.6), /,
     &       '  l1,    l2,   net,   nej', 4(1x,e13.6) )
 7012 format('temp: ', e13.6)
c
      end

c
c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_gp_elastic_fraction                                 *
c *                 -- last modified 8/18/11 by jcs                 *
c *                                                                 *
c *    this subroutine computes the fraction of the the strain      *
c *    increment that produces purely elastic response and updates  *
c *    the stress_n vector to include the stress caused by this     *
c *    purely elastic strain fraction                               *
c *                                                                 *
c *******************************************************************
c
       subroutine mm05_gp_elastic_fraction( stress_n, alpha_n,
     &    devdeps, eta, k_n, iout, felem, gpn, i, step,
     &    g_n, g_n1, g_s, Hi_n, Hi_n1, Hi_s, Hk_n, Hk_n1, Hk_s,
     &    beta_n, beta_n1, beta_s, delta_n, delta_n1, delta_s,
     &    k0_n, k0_n1, k0_s, eps_n, adaptive_possible,
     &    cut_step_size_now, onestep )
       implicit none
c
c    parameter declarations
c    ----------------------
c
      logical adaptive_possible, cut_step_size_now, onestep
      integer iout, felem, gpn, i, step
      double precision
     &    stress_n(*), alpha_n(*), devdeps(*), eta, k_n,
     &    g_n, g_n1, g_s, Hi_n, Hi_n1, Hi_s, Hk_n, Hk_n1, Hk_s,
     &    beta_n, beta_n1, beta_s, delta_n, delta_n1, delta_s,
     &    k0_n, k0_n1, k0_s, eps_n
c
c    local  parameter declarations
c    -----------------------------
c
      logical local_debug
      integer l, loop
      double precision
     &    zero, half, one, two, tol, stol, eta_low, eta_high,
     &    eta_mid, eta_new, eta_ridder, f_low, f_high, f_mid, f_new,
     &    sridder, eta_comp, unused, t1, t2, t3, t4, t5, t6, lk,
     &    eta_pt, f_pt, ptol, df, a, b, c, x1, x2, fx1, fx2
c
      data zero, half, one, two, tol, stol, eta_comp, unused, ptol, c
     &     / 0.0, 0.5, 1.0, 2.0, 1.0e-10, 1.0e-10, 1.0e-10, -1.11e30,
     &       0.0001, 0.6180339887 /
      local_debug = .false.
c
c     determine the six constant scalar variables from contracting
c     two tensors
c
      t1 =    stress_n(1)**2 + stress_n(2)**2 + stress_n(3)**2 +
     &  two*( stress_n(4)**2 + stress_n(5)**2 + stress_n(6)**2 )
c
      t2 =    stress_n(1)*alpha_n(1) +
     &        stress_n(2)*alpha_n(2) +
     &        stress_n(3)*alpha_n(3) +
     &  two*( stress_n(4)*alpha_n(4) +
     &        stress_n(5)*alpha_n(5) +
     &        stress_n(6)*alpha_n(6) )
c
      t3 =    alpha_n(1)**2 + alpha_n(2)**2 + alpha_n(3)**2 +
     &  two*( alpha_n(4)**2 + alpha_n(5)**2 + alpha_n(6)**2 )
c
      t4 =    devdeps(1)*stress_n(1) +
     &        devdeps(2)*stress_n(2) +
     &        devdeps(3)*stress_n(3) +
     &        devdeps(4)*stress_n(4) +
     &        devdeps(5)*stress_n(5) +
     &        devdeps(6)*stress_n(6)
c
      t5 =    devdeps(1)*alpha_n(1) +
     &        devdeps(2)*alpha_n(2) +
     &        devdeps(3)*alpha_n(3) +
     &        devdeps(4)*alpha_n(4) +
     &        devdeps(5)*alpha_n(5) +
     &        devdeps(6)*alpha_n(6)
c
      t6 =    devdeps(1)**2 + devdeps(2)**2 + devdeps(3)**2 +
     & half*( devdeps(4)**2 + devdeps(5)**2 + devdeps(6)**2 )
c
      if ( local_debug ) then
         write(iout,9000) felem+i-1, gpn, step
         write(iout,9001) t1, t2, t3, t4, t5, t6
      end if
c
c      if ( onestep ) then
c         eta = zero
c         goto 1000
c      end if
c
c     the yield function for an entirely elastic strain
c     increment is nonlinear for temperature-dependent
c     material properties. the value of eta remains bounded
c     between 0 (fully elastic-plastic) and 1 (fully elastic).
c     the yield function is negative for eta = 0 and positive
c     for eta = 1. the (relatively) small material property changes
c     with temperature suggest that only one root exists between
c     0 < eta < 1.
c
c     consequently, we adopt a robust variant of relgula-falsi called
c     Ridder's method (see Numerical Recipes). previous experience
c     with this nonlinear solution technique suggests that the value
c     of eta should converge at a nearly quadratic rate. Ridder's method
c     also does not require a derivative of the yield function.
c
c     Tte subroutine mm05_gp_felas_eta determines the yield function
c     value for a given function of eta.
c
c     1. set the lower and upper bounds on eta. verify the
c        root (deplas) is bracketed. determine if the upper/
c        lower bounds are in fact the solution.
c
      eta_low  = zero
      eta_high = one
c      f_low = func( eta_low, ... )
      call mm05_gp_felas_eta(
     &  eta_low, f_low, iout, felem, gpn, i,
     &  g_n, g_n1, Hi_n, Hi_n1, Hk_n, Hk_n1, k0_n, k0_n1,
     &  t1, t2, t3, t4, t5, t6, eps_n )
c      f_high = func( eta_high, ... )
      call mm05_gp_felas_eta(
     &  eta_high, f_high, iout, felem, gpn, i,
     &  g_n, g_n1, Hi_n, Hi_n1, Hk_n, Hk_n1, k0_n, k0_n1,
     &  t1, t2, t3, t4, t5, t6, eps_n )
c
      if ( local_debug ) then
         write(iout,9002) eta_low, f_low
         write(iout,9004) eta_high, f_high
      end if
c
c     if the root is properly bracketed, proceed to Ridder's method
c
      if ( (f_low .lt. zero .and. f_high .gt. zero) ) go to 100
c
c     check if stress state is outside of the yield surface at the
c     start of the load step. this can occur in the GP model.
c     need to check the secant of the yield function near eta_low
c     if it's negative, then we perform a golden section search
c     to determine the minimum value of the yield function.
c     otherwise, eta is zero
c
      if ( f_low .gt. zero .and. f_high .gt. zero ) then
         eta_pt = eta_low + ptol
c      f_pt = func( etapt, ... )
         call mm05_gp_felas_eta(
     &   eta_pt, f_pt, iout, felem, gpn, i,
     &   g_n, g_n1, Hi_n, Hi_n1, Hk_n, Hk_n1, k0_n, k0_n1,
     &   t1, t2, t3, t4, t5, t6, eps_n )
         df = (f_pt-f_low)/ptol
c
         if ( df .gt. zero ) then
            eta = eta_low
            go to 1000
         else
            go to 200
         end if
      end if
c
c     check remaining cases where the root is improperly bracketed
c
      if ( abs(f_low) .le. stol*k_n  ) then
         eta = eta_low
      elseif ( abs(f_high) .le. stol*k_n ) then
         write(iout, 9200) felem+i-1, gpn
         call die_abort
      else
c
c     general catch all case
c
         write(iout, 9201) felem+i-1, gpn
         write(iout, 9002) eta_low,  f_low
         write(iout, 9004) eta_high, f_high
         call die_abort
      end if
      go to 1000
c
c     2. execute a maximum of 50 iterations of regula-falsi
c        with Ridder's improvements.
c
 100  continue
      eta_ridder = unused
c
      do loop = 1, 50
c
c     2a. evaluate the residual function at the
c         current mid-point of bracketed range.
c         compute Ridder's "magic" factor:)
c
         eta_mid = half * ( eta_low + eta_high )
c         f_mid = func( eta_mid )
c
         call mm05_gp_felas_eta(
     &        eta_mid, f_mid, iout, felem, gpn, i,
     &        g_n, g_n1, Hi_n, Hi_n1, Hk_n, Hk_n1, k0_n, k0_n1,
     &        t1, t2, t3, t4, t5, t6, eps_n )
c
         sridder = sqrt( f_mid**2 - f_low*f_high )
         if ( sridder .eq. zero ) then
            write(iout,9202)
            if ( adaptive_possible ) then
               cut_step_size_now = .true.
               return
            else
               write(iout,9205)
               call die_abort
            end if
         end if
c
c     2b. compute the new estimate of eta
c         using Ridder's special update.
c         check for convergence of eta
c
         eta_new = eta_mid + (eta_mid-eta_low) *
     &             (sign(one,f_low-f_high)*f_mid/sridder)
         if ( abs(eta_new-eta_ridder) .le. eta_comp ) then
            eta    = eta_new
            go to 1000
         end if
c
c     2c. evaluate residual function using the
c         Ridder estimate for eta. check
c         convergence of the actual residual function.
c
         eta_ridder =  eta_new
c         f_new = func( eta_ridder, ... )
c
         call mm05_gp_felas_eta(
     &        eta_ridder, f_new, iout, felem, gpn, i,
     &        g_n, g_n1, Hi_n, Hi_n1, Hk_n, Hk_n1, k0_n, k0_n1,
     &        t1, t2, t3, t4, t5, t6, eps_n )
c
         if ( abs(f_new) .le. stol*k_n ) then
            eta = eta_new
            go to 1000
         end if
c
c     2d. must do another iteration. bookeeping
c         of high/low values to keep the root bracketed.
c
         if ( sign(f_mid,f_new) .ne. f_mid ) then
            eta_low = eta_mid
            f_low = f_mid
            eta_high = eta_ridder
            f_high = f_new
         elseif ( sign(f_low,f_new) .ne. f_low ) then
            eta_high = eta_ridder
            f_high = f_new
         elseif ( sign(f_high,f_new) .ne. f_high ) then
            eta_low = eta_ridder
            f_low  = f_new
         else
            write(iout,9203)
            cut_step_size_now = .true.
            return
         end if
c
c     2e. check convergence of new brackets on
c         eta - they may have collapsed to within
c         the tolerance.
c
         if ( abs(eta_high-eta_low) .le. eta_comp ) then
            eta    = ( eta_high + eta_low ) * half
            go to 1000
         end if
c
c     2f. all done with this iteration. fall thru
c         at limit of iterations, set flag, and
c         return.
c
      end do
c
      write(iout,9204)
      if ( adaptive_possible ) then
         cut_step_size_now = .true.
         return
      else
         write(iout,9205)
         call die_abort
      end if
c
c     3. perform a golden section search to detemine the minimum
c        value of the yield function between eta_low and eta_high
c
 200  continue
c
c     3a. set the initial variables prior to starting loop
c
      a = eta_low
      b = eta_high
      x1 = c*a + (one-c)*b
      call mm05_gp_felas_eta(
     &     x1, fx1, iout, felem, gpn, i,
     &     g_n, g_n1, Hi_n, Hi_n1, Hk_n, Hk_n1, k0_n, k0_n1,
     &     t1, t2, t3, t4, t5, t6, eps_n )
      x2 = (one-c)*a + c*b
      call mm05_gp_felas_eta(
     &     x2, fx2, iout, felem, gpn, i,
     &     g_n, g_n1, Hi_n, Hi_n1, Hk_n, Hk_n1, k0_n, k0_n1,
     &     t1, t2, t3, t4, t5, t6, eps_n )
c
c     3b. begin the loop that will eventually determine the
c         value of eta where f is a minimum value.
c
      do loop = 1,10000
c
c     3c. based on the value of fx1 vs fx2, update the variables
c
         if ( fx1 .lt. fx2 ) then
            b   = x2
            x2  = x1
            fx2 = fx1
            x1  = c*a + (one-c)*b
            call mm05_gp_felas_eta(
     &           x1, fx1, iout, felem, gpn, i,
     &           g_n, g_n1, Hi_n, Hi_n1, Hk_n, Hk_n1, k0_n, k0_n1,
     &           t1, t2, t3, t4, t5, t6, eps_n )
c
         else
            a   = x1
            x1  = x2
            fx1 = fx2
            x2  = (one-c)*a + c*b
            call mm05_gp_felas_eta(
     &           x2, fx2, iout, felem, gpn, i,
     &           g_n, g_n1, Hi_n, Hi_n1, Hk_n, Hk_n1, k0_n, k0_n1,
     &           t1, t2, t3, t4, t5, t6, eps_n )
         end if
c
c     3d. check if the difference between a and b meets the tolerance
c         for eta. if so, set eta to the mean value of a and b and exit
c         the loop
c
         if ( abs(b-a) .lt. eta_comp ) then
            eta = half*(a+b)
            go to 1000
         end if
c
c     3e. all done with this iteration. fall thru
c         at limit of iterations, set flag, and
c         return.
c
      end do
c
      write(iout,9206)
      if ( adaptive_possible ) then
         cut_step_size_now = .true.
         return
      else
         write(iout,9205)
         call die_abort
      end if

c
c     4. iterations converged on eta. compute the new values
c        of the stress, backstress, and material properties
c        based on the current value of eta.
c
 1000 continue
c
      if ( local_debug ) then
         call mm05_gp_felas_eta(
     &        eta, f_new, iout, felem, gpn, i,
     &        g_n, g_n1, Hi_n, Hi_n1, Hk_n, Hk_n1, k0_n, k0_n1,
     &        t1, t2, t3, t4, t5, t6, eps_n )
         write(iout,9005) eta, f_new
      end if
c
      g_s     = eta*g_n1 + (one-eta)*g_n
      Hi_s    = eta*Hi_n1 + (one-eta)*Hi_n
      Hk_s    = eta*Hk_n1 + (one-eta)*Hk_n
      k0_s    = eta*k0_n1 + (one-eta)*k0_n
      beta_s  = eta*beta_n1 + (one-eta)*beta_n
      delta_s = eta*delta_n1 + (one-eta)*delta_n
c
      stress_n(1) = g_s/g_n*stress_n(1) + eta*2*g_s*devdeps(1)
      stress_n(2) = g_s/g_n*stress_n(2) + eta*2*g_s*devdeps(2)
      stress_n(3) = g_s/g_n*stress_n(3) + eta*2*g_s*devdeps(3)
      stress_n(4) = g_s/g_n*stress_n(4) + eta*g_s*devdeps(4)
      stress_n(5) = g_s/g_n*stress_n(5) + eta*g_s*devdeps(5)
      stress_n(6) = g_s/g_n*stress_n(6) + eta*g_s*devdeps(6)
c
      lk = one
      if ( abs( Hk_n ) .gt. tol ) lk = Hk_s/Hk_n
      alpha_n(1:6) = lk*alpha_n(1:6)
c
      devdeps(1:6) = (1-eta)*devdeps(1:6)
c
      if ( local_debug ) then
         write(iout,9010) (stress_n(l), l=1,6)
         write(iout,9011) (alpha_n(l), l=1,6)
         write(iout,9012) (devdeps(l), l=1,6)
      end if
c
      return
c
 9000 format('Debug mm05_gp_elastic_fraction:  element ', i7, ' at gp ',
     &   i3, ' step number ', i6  )
 9001 format( 't1 t2 t3:', 3e15.6, /, 't4 t5 t6:', 3e15.6 )
 9002 format( 'eta_low  f_low: ', 2e15.6 )
 9003 format( 'eta_mid  f_mid: ', 2e15.6 )
 9004 format( 'eta_high f_high:', 2e15.6 )
 9005 format( 'eta f           ', 2e15.6 )
 9010 format( 'stress_n = ', 3e15.6, /, '           ', 3e15.6 )
 9011 format( 'alpha_n =  ', 3e15.6, /, '           ', 3e15.6 )
 9012 format( 'devdeps =  ', 3e15.6, /, '           ', 3e15.6 )
 9200 format( '>> Fatal Error: routine mm05_gp_elastic_fraction:',
     &        '  element ', i7, ' gauss point ', i3, ' eta = 1 ',
     &         ' Job terminated.' // )
 9201 format( '>> Fatal Error: routine mm05_gp_elastic_fraction:',
     &         '   element ', i7, ' gauss point ', i3, /,
     &         '   eta is not properly bracketed.', /,
     &         ' Job terminated.' // )
 9202 format(/,3x,'>> Error: divide by zero in ',
     &         'mm05_gp_elastic_fraction')
 9203 format(/,3x,'>> Error: failed to make new bracket',
     &       /,3x,'          for Ridder update of eta')
 9204 format(/,3x,'>> Error: Ridder method failed to converge to eta')
 9205 format( '>> Fatal Error: cannot requst step size reduction',
     &        ' Job terminated.' // )
 9206 format(/,3x,'>> Error: Golden section earch fails to converge',
     &       ' to eta')
c
      end


c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        subroutine mm05_gp_felas_eta                             *
c *                 -- last modified 7/05/11 by jcs                 *
c *                                                                 *
c *        Computes the temperature-dependent yield function for    *
c *        the GP model, assuming an entirely elastic load step     *
c *        with fraction eta                                        *
c *                                                                 *
c *******************************************************************
c
c
      subroutine mm05_gp_felas_eta(
     &  eta, felas, iout, felem, gpn, i,
     &  g_n, g_n1, Hi_n, Hi_n1, Hk_n, Hk_n1, k0_n, k0_n1,
     &  t1, t2, t3, t4, t5, t6, eps_n )
      implicit none
c
c    parameter declarations
c    ----------------------
c
      integer iout, felem, gpn, i
      double precision
     &    eta, felas, g_n, g_n1, Hi_n, Hi_n1, Hk_n, Hk_n1,
     &    k0_n, k0_n1, t1, t2, t3, t4, t5, t6, eps_n

c
c    local  parameter declarations
c    -----------------------------
c
      double precision
     &    g_b, Hi_b, Hk_b, k0_b, c1, c2, a1, a2, a3, a4,
     &    one, four, tol
c
      data one, four, tol / 1.0, 4.0, 0.00000001 /
c
c     initialize material properties for the temperature at eta
c
      g_b     = eta*g_n1  + (one-eta)*g_n
      Hi_b    = eta*Hi_n1 + (one-eta)*Hi_n
      Hk_b    = eta*Hk_n1 + (one-eta)*Hk_n
      k0_b    = eta*k0_n1 + (one-eta)*k0_n
c
      c1 = g_b/g_n
      c2 = one
      if ( abs( Hk_n ) .gt. tol ) c2 = Hk_b/Hk_n
c
c     determine the quantities a1, a2, a3, and a4 that define
c     the yield function at temperature eta
c
      a1 = (c1**2)*t1 - 2*c1*c2*t2 + (c2**2)*t3
      a2 = four*g_b*eta*( c1*t4 - c2*t5 )
      a3 = four*(g_b**2)*(eta**2)*t6
      a4 = k0_b + Hi_b*eps_n
c
c     evaluate the yield function
c
      felas = sqrt( a1 + a2 + a3 ) - a4

      return

      end

c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *         subroutine mm05_gp_debug -- prints out variables        *
c *                                                                 *
c *******************************************************************
c
c
      subroutine mm05_gp_debug(
     &  step, iter, felem, gpn, mxvl, hist_size, nstrs, nstrn, span,
     &  iout, signal_flag, adaptive_possible, cut_step_size_now,
     &  mm_props, e_vec_np1, nu_vec_np1, sigyld_vec_np1,
     &  e_vec_n, nu_vec_n, sigyld_vec_n, h_gp_np1, h_gp_n,
     &  beta_gp_np1, beta_gp_n, delta_gp_np1, delta_gp_n, tau,
     &  trial_elas_stress_np1, stress_n, stress_np1,
     &  deps, history_n, history_np1 )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &  step, iter, felem, gpn, mxvl, hist_size,
     &  span, iout, nstrs, nstrn
c
      logical
     &   signal_flag, adaptive_possible, cut_step_size_now
c
      double precision
     & mm_props(mxvl,5), e_vec_np1(mxvl), nu_vec_np1(mxvl),
     & sigyld_vec_np1(mxvl), stress_n(mxvl,nstrs),
     & e_vec_n(mxvl), nu_vec_n(mxvl), sigyld_vec_n(mxvl),
     & h_gp_np1(mxvl), h_gp_n(mxvl), beta_gp_np1(mxvl), beta_gp_n(mxvl),
     & delta_gp_np1(mxvl), delta_gp_n(mxvl), tau(mxvl),
     & stress_np1(mxvl,nstrs), deps(mxvl,nstrn),
     & trial_elas_stress_np1(mxvl,nstrn), history_n(span,hist_size),
     & history_np1(span,hist_size)
c
c                          local parameters
c
      integer i, j
c
          if ( gpn .ne. 1 ) return
c
c                         print out all input variables
c
          write( iout, * ) 'debug mode - echoing input variables'
          write( iout, 9000 ) step, iter, felem, gpn, mxvl, hist_size,
     &          span, nstrs, nstrn,
     &          signal_flag, adaptive_possible, cut_step_size_now
c          do i =1, span
           i = 1
           write( iout, 9010 )  e_vec_n(i), nu_vec_n(i),
     &        sigyld_vec_n(i), h_gp_n(i), beta_gp_n(i),
     &        delta_gp_n(i), e_vec_np1(i), nu_vec_np1(i),
     &        sigyld_vec_np1(i), h_gp_np1(i), beta_gp_np1(i),
     &        delta_gp_np1(i), tau(i),
     &        ( stress_n(i,j), j=1, nstrs ),
     &        ( stress_np1(i,j), j=1, nstrs ),
     &        ( deps(i,j), j=1, nstrn ),
     &        ( trial_elas_stress_np1(i,j), j=1, nstrn ),
     &        ( history_n(i,j), j=1, hist_size ),
     &        ( history_np1(i,j), j=1, hist_size )
c          end do
c
       return
c
c
 9000 format( '  step = ', i3, ',  iteration = ', i4,',  felem = ', i3,
     & /, '  gpn = ', i3,  ',  mxvl = ', i3, ',  hist_size = ', i3,
     & /, '  span = ', i3, ',  nstrs = ', i3, ',  nstrn = ', i3,
     & /, '  signal_flag = ', l3, '  adapt. = ', l3, '  cut_step = ',
     &    l3 )
 9010 format( '  @n  : e,   nu,     yld_pt  = ', 3e15.6,
     & /,     '        h_u, beta_u, delta_u = ', 3e15.6,
     & /,     '  @n+1: e,   nu,     yld_pt  = ', 3e15.6,
     & /,     '        h_u, beta_u, delta_u = ', 3e15.6,
     & /,     '  tau:', e12.4,
     & /,     '  stress @ n   = ',    3e15.6, / , 17x, 3e15.6,
     & /,                                         17x, 3e15.6,
     & /,     '  stress @ n+1 = ',  3e15.6, / , 17x, 3e15.6,
     & /,                                         17x, 3e15.6,
     & /,     '  deps         = ',  3e15.6, / , 17x, 3e15.6,
     & /,     '  trial stress = ',  3e15.6, / , 17x, 3e15.6,
     & /,     '  history @ n  = ', 4e15.6, / , 17x, 4e15.6,
     & /,                                      17x, 4e15.6,
     & /,                                      17x, 2e15.6,
     & /,     '  history @ n+1= ', 4e15.6, / , 17x, 4e15.6,
     & /,                                      17x, 4e15.6,
     & /,                                      17x, 2e15.6 )
c
       end


c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        subroutine cnst5_gp -- computes consistent tangent       *
c *             -- last modified 1/1/2016 rhd                       *
c *                                                                 *
c *******************************************************************
c
c
      subroutine cnst5_gp(
     &  span, felem, gpn, iter, iout, mxvl, nstrn,
     &  e_vec_np1, nu_vec_np1, h_gp_np1, beta_gp_np1,
     &  delta_gp_np1, tau, sig_trial, history_n, history_np1,
     &  stress_np1, dmat )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer :: span, felem, gpn, iter, iout, mxvl, nstrn
c
      double precision ::
     & e_vec_np1(mxvl), nu_vec_np1(mxvl), h_gp_np1(mxvl),
     & beta_gp_np1(mxvl), delta_gp_np1(mxvl), tau(mxvl),
     & sig_trial(mxvl,nstrn), history_n(span,*),
     & history_np1(span,*), dmat(mxvl,nstrn,nstrn),
     & stress_np1(mxvl,nstrn)
c
c               description of parameters
c               -------------------------
c
c     step              : current load step number
c     iter              : current newton iteration number. iter 1 is
c                         for application of the "real" load increment.
c     felem             : first element of the current block
c     gpn               : gauss point number being processed for block
c     mxvl              : maximum no. elements per block
c     nstrn             : number of strain-stress components (=6)
c     span              : number of elements in current block
c     iout              : write messages to this device number
c     e_vec_np1         : Young's modulus for each element in block at end
c                         of load step (n+1)
c     nu_vec_np1        : Poisson's ratio for each element in block at end
c                         of load step (n+1)
c     h_gp_np1          : terminal hardening modulus of the gp model for
c                         each element in block at end of load step (n+1)
c     beta_gp_np1       : beta parameter of the gp model for each eleemnt in
c                         block at end of load step (n+1)
c     delta_gp_np1      : delta parameter of the gp model for each element
c                         in block at end of load step (n+1)
c     tau               : measure of kinematic/isotropic hardening in the
c                         gp model for each element in block
c (#) trial_elas_stress_np1 : trial elastic stress vector defined by stress
c                             update routine
c                         consistent tangent routine for model
c     history_n         : history values at start of load step (n) for all
c                         elements in block for this gauss point
c     history_np1       : history values at end of load step (n+1) for all
c                         elements in block for this gauss point
c     det_jac_block     : |J| at this gauss point for each element in
c                         block
c (!) stress_np1        : current estimate of 6 stress components for
c                         end of step (see ordering below)
c     weight            : integration point weight factor
c (*) dmat              : 6x6 (symmetric) tangent (consistent) for
c                         this gauss point for each element of block
c                         (see stress ordering below)
c
c    (*)  values to be updated by this material model
c    (!)  for finite strain computations, these are unrotated
c         Cauchy stresses
c    (#)  used by constitutive update procedures based on some
c         for of elastic predictor - return mapping algorithm.
c         the contents of this array are set by the corresponding
c         stress update routine for the model. for finite
c         strain computations, the contents will be unrotated
c         Cauchy stress terms of some form as set by the stress
c         update routine.
c
c   Finite strain issues:
c     this constitutive model only sees stresses that have already
c     been rotation "neutralized" by WARP3D. this routine can
c     thus operate simply as "small strain" theory. WARP3D will
c     "rotate" the compute [D] matrices on return for finite strain-
c     rotation effects.
c
c   strain ordering:
c     deps-xx, deps-yy, deps-zz, gamma-xy, gamma-yz, gamma-xz
c
c   stress ordering (at n and n+1):
c     sig-xx, sig-yy, sig-zz, tau-xy, tau-yz, tau-xz
c
c                      local variables
c                      ----------------
c
      integer :: i, iword(2), state, l, m , j, t
      logical :: yield(mxvl), debug
      double precision ::
     &     c1, c2, c3, c4, fact, dword,
     &     n(6), equiv(6), norm_equiv_np1, g, bulk,
     &     A_GP, C_GP, mean_n, norm_equiv_n, lambda,
     &     mean_sig, G1, BB1, BB2, BB3, BB4,
     &     k_n, k_n1, Hk_n1, H_n1, dHi, f_n1, f_n, del_n1, beta_n1,
     &     zero, one, two, twothird, onethird, half, rt2third,
     &     three, four

          data zero, one, two, twothird, onethird, half,
     &         rt2third, three, four  /0.0d0, 1.0d0, 2.0d0,
     &     0.66666666d0, 0.33333333d0, 0.5d0, 0.81649658d0, 3.0d0,
     &     4.0d0 /
c
       equivalence(dword, iword)
c
       do i = 1, span
         dword = history_np1(i,4)
         state = iword(1)
         yield(i) = .false.
         if( state .eq. 1 ) yield(i) = .true.
       end do
c
      do i = 1, span
       if( yield(i) ) cycle
       dmat(i,1,4) = zero
       dmat(i,1,5) = zero
       dmat(i,1,6) = zero
       dmat(i,2,4) = zero
       dmat(i,2,5) = zero
       dmat(i,2,6) = zero
       dmat(i,3,4) = zero
       dmat(i,3,5) = zero
       dmat(i,3,6) = zero
       dmat(i,4,1) = zero
       dmat(i,4,2) = zero
       dmat(i,4,3) = zero
       dmat(i,4,5) = zero
       dmat(i,4,6) = zero
       dmat(i,5,1) = zero
       dmat(i,5,2) = zero
       dmat(i,5,3) = zero
       dmat(i,5,4) = zero
       dmat(i,5,6) = zero
       dmat(i,6,1) = zero
       dmat(i,6,2) = zero
       dmat(i,6,3) = zero
       dmat(i,6,4) = zero
       dmat(i,6,5) = zero
c
       c1 = e_vec_np1(i)/((one+nu_vec_np1(i)) *
     &         (one-two*nu_vec_np1(i)))
       c2 = (one-nu_vec_np1(i))*c1
       c3 = ((one-two*nu_vec_np1(i))/two)*c1
       c4 = nu_vec_np1(i)*c1
c
       dmat(i,1,1) = c2
       dmat(i,2,2) = c2
       dmat(i,3,3) = c2
       dmat(i,4,4) = c3
       dmat(i,5,5) = c3
       dmat(i,6,6) = c3
       dmat(i,1,2) = c4
       dmat(i,1,3) = c4
       dmat(i,2,1) = c4
       dmat(i,3,1) = c4
       dmat(i,2,3) = c4
       dmat(i,3,2) = c4
      end do
c
      do i = 1, span
       if( .not. yield(i) ) cycle
       fact     = 1.0d00
       g        = e_vec_np1(i)/(two*(one+nu_vec_np1(i)))
       bulk     = e_vec_np1(i)/(three*(one-two*nu_vec_np1(i)))
       del_n1   = twothird*delta_gp_np1(i)
       beta_n1  = rt2third*beta_gp_np1(i)
       H_n1     = twothird*h_gp_np1(i)
       Hk_n1    = ( one - tau(i) )*H_n1
       G1       = g+half*H_n1
c
       lambda   = history_np1(i,1)
       k_n1     = history_np1(i,2)
       dHi      = history_np1(i,14)
c
       norm_equiv_n = history_np1(i,5)
       k_n          = history_np1(i,13)
c
       mean_sig = ( stress_np1(i,1)+stress_np1(i,2) +
     &         stress_np1(i,3) )/three
c
       equiv(1) = stress_np1(i,1)-history_np1(i, 6)-mean_sig
       equiv(2) = stress_np1(i,2)-history_np1(i, 7)-mean_sig
       equiv(3) = stress_np1(i,3)-history_np1(i, 8)-mean_sig
       equiv(4) = stress_np1(i,4)-history_np1(i, 9)
       equiv(5) = stress_np1(i,5)-history_np1(i,10)
       equiv(6) = stress_np1(i,6)-history_np1(i,11)
c
       norm_equiv_np1 = sqrt(equiv(1)**2+equiv(2)**2+equiv(3)**2 +
     &         two*(equiv(4)**2+equiv(5)**2+equiv(6)**2 ))
c
       f_n  = norm_equiv_n   - k_n
       f_n1 = norm_equiv_np1 - k_n1
c
       BB1 = H_n1 + del_n1 - dHi
       BB2 = ( del_n1 + H_n1 )*beta_n1
       BB3 = two*f_n1 - f_n + BB1*lambda
       BB4 = BB2 - BB1*f_n1
c
       A_GP = two*g*BB3/( two*G1*BB3 + BB4 )
       C_GP = two*g*lambda/( norm_equiv_np1 + (two*g+Hk_n1)*lambda )
c
c               elasitc-plastic weighting terms
c
       c1 = two*g*(1-C_GP)
       c2 = two*g*(C_GP-A_GP)
c
       n(1) = equiv(1)/norm_equiv_np1
       n(2) = equiv(2)/norm_equiv_np1
       n(3) = equiv(3)/norm_equiv_np1
       n(4) = equiv(4)/norm_equiv_np1
       n(5) = equiv(5)/norm_equiv_np1
       n(6) = equiv(6)/norm_equiv_np1
c
c               diagonal terms
c
       mean_n  = onethird*(n(1) + n(2) + n(3) )
       dmat(i,1,1) = bulk + twothird*c1 + c2*n(1)*n(1)
       dmat(i,2,2) = bulk + twothird*c1 + c2*n(2)*n(2)
       dmat(i,3,3) = bulk + twothird*c1 + c2*n(3)*n(3)
       dmat(i,4,4) = half*c1 + c2*n(4)*n(4)
       dmat(i,5,5) = half*c1 + c2*n(5)*n(5)
       dmat(i,6,6) = half*c1 + c2*n(6)*n(6)
c
c               off diagonal terms
c
       dmat(i,1,2) = (bulk - onethird*c1 + c2*n(1)*n(2))
       dmat(i,1,3) = (bulk - onethird*c1 + c2*n(1)*n(3))
       dmat(i,1,4) = (c2*n(1)*n(4))
       dmat(i,1,5) = (c2*n(1)*n(5))
       dmat(i,1,6) = (c2*n(1)*n(6))
c
       dmat(i,2,3) = bulk - onethird*c1 + c2*n(2)*n(3)
       dmat(i,2,4) = c2*n(2)*n(4)
       dmat(i,2,5) = c2*n(2)*n(5)
       dmat(i,2,6) = c2*n(2)*n(6)
c
       dmat(i,3,4) = c2*n(3)*n(4)
       dmat(i,3,5) = c2*n(3)*n(5)
       dmat(i,3,6) = c2*n(3)*n(6)
c
       dmat(i,4,5) = c2*n(4)*n(5)
       dmat(i,4,6) = c2*n(4)*n(6)
c
       dmat(i,5,6) = c2*n(5)*n(6)
c
c                symmetry
c
       dmat(i,2,1) = dmat(i,1,2)
       dmat(i,3,1) = dmat(i,1,3)
       dmat(i,4,1) = dmat(i,1,4)
       dmat(i,5,1) = dmat(i,1,5)
       dmat(i,6,1) = dmat(i,1,6)
c
       dmat(i,3,2) = dmat(i,2,3)
       dmat(i,4,2) = dmat(i,2,4)
       dmat(i,5,2) = dmat(i,2,5)
       dmat(i,6,2) = dmat(i,2,6)
c
       dmat(i,4,3) = dmat(i,3,4)
       dmat(i,5,3) = dmat(i,3,5)
       dmat(i,6,3) = dmat(i,3,6)
c
       dmat(i,5,4) = dmat(i,4,5)
       dmat(i,6,4) = dmat(i,4,6)
c
       dmat(i,6,5) = dmat(i,5,6)
      end do
c
      return
c
 9000 format( '  iteration = ', i4,',  felem = ', i3,
     & /, '  gpn = ', i3,  ', ','  span = ', i3  )
 9010 format( '  @n+1: e,   nu,             = ', 2e15.6,
     & /,     '        h_u, beta_u, delta_u = ', 3e15.6,
     & /,     '  tau:', e12.4,
     & /,     '  stress @ n+1 = ',  3e15.6, / , 17x, 3e15.6,
     & /,     '  history @ n+1= ', 4e15.6, / , 17x, 4e15.6,
     & /,                                      17x, 4e15.6,
     & /,                                      17x, 2e15.6 )
 9020 format('Consistent tangent')
 9030 format(6(1x,e15.8))
c
      end

c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_final  -- last modified 11/14/02 by kbc             *
c *                                                                 *
c *        calculate total update stress (deviatoric component      *
c *        only was computed in mm05_compute) and calculate the     *
c *        energy density from a trapezoidal numerical integration  *
c *        of increments of strain and average stresses             *
c *                                                                 *
c *******************************************************************
       subroutine mm05_final(
     &                  span, mxvl, stress_n, stress_np1, e_vec,
     &                  nu_vec, shear_mod_vec, deps, instat,
     &                  history_np1, trace_eps_np1, debug, iout )

       implicit none
c
c                    parameter declarations
c                    ----------------------
       integer span, mxvl, instat(*), iout
       logical debug
      double precision
     & stress_n(mxvl, *), stress_np1(mxvl, *), e_vec(*),
     & nu_vec(*), shear_mod_vec(*), deps(mxvl, *),
     & history_np1(span, *), trace_eps_np1(*)
c
c                   local parameter declarations
c                   ----------------------------
       integer i, iword(2)
      double precision
     & sig_mean_np1, one, two, three, half, dword
       equivalence ( dword, iword )
       data one, two, three, half / 1.0, 2.0, 3.0, 0.5 /
c
c
c
       do i = 1, span

c          sig_mean_np1 = trace_eps_np1(i)*e_vec(i)/
c     &                   ( three*(one-2*nu_vec(i)))

          sig_mean_np1 = trace_eps_np1(i)
     &        *(three*e_vec(i)*nu_vec(i)/((one+nu_vec(i))*
     &       (one-two*nu_vec(i)))+two*shear_mod_vec(i))/three
c
          stress_np1(i,1) = stress_np1(i,1) + sig_mean_np1
          stress_np1(i,2) = stress_np1(i,2) + sig_mean_np1
          stress_np1(i,3) = stress_np1(i,3) + sig_mean_np1
c
          stress_np1(i,7) = stress_n(i,7) + half*(
     &               deps(i,1)*(stress_n(i,1)+stress_np1(i,1))
     &             + deps(i,2)*(stress_n(i,2)+stress_np1(i,2))
     &             + deps(i,3)*(stress_n(i,3)+stress_np1(i,3))
     &             + deps(i,4)*(stress_n(i,4)+stress_np1(i,4))
     &             + deps(i,5)*(stress_n(i,5)+stress_np1(i,5))
     &             + deps(i,6)*(stress_n(i,6)+stress_np1(i,6)))
c
          iword(1) = instat(i)
          history_np1(i,4) = dword
c
       end do

       return
c
       end

c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity model       *
c *                                                                 *
c *        mm05_plastic_work  -- last modified 11/15/02 by kbc      *
c *                                                                 *
c *        calculate the total plastic work density and the         *
c *        total plastic strain and store in stress vector          *
c *                                                                 *
c *******************************************************************

      subroutine  mm05_plastic_work( iout, span, mxvl, stress_n,
     &     stress_np1, yield, deps, nu_vec, e_vec, shear_mod_vec )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer  span, mxvl, iout
c
      logical yield(*)
c
      double precision
     &   stress_n(mxvl,*), stress_np1(mxvl,*), deps(mxvl,*),
     &   nu_vec(*), e_vec(*), shear_mod_vec(*)
c
c                       local parameters
c
       integer i
c
      double precision
     & zero, deps_plas_bar, dsig(6), deps_plas(6), half, root2,
     & three, two, factor1, factor2
c
       data zero, half, root2, three, two / 0.0, 0.5, 1.414213562373095,
     & 3.0, 2.0 /
c
c
c
       do i = 1, span
        stress_np1(i,8)   =    stress_n(i,8)
        stress_np1(i,9)   =    stress_n(i,9)
        if ( yield(i)) then
          deps_plas_bar = zero
          dsig(1) = stress_np1(i,1) - stress_n(i,1)
          dsig(2) = stress_np1(i,2) - stress_n(i,2)
          dsig(3) = stress_np1(i,3) - stress_n(i,3)
          dsig(4) = stress_np1(i,4) - stress_n(i,4)
          dsig(5) = stress_np1(i,5) - stress_n(i,5)
          dsig(6) = stress_np1(i,6) - stress_n(i,6)
          deps_plas(1) = deps(i,1) -
     &                   (dsig(1)-nu_vec(i)*(dsig(2)+dsig(3)))/e_vec(i)
          deps_plas(2) = deps(i,2) -
     &                   (dsig(2)-nu_vec(i)*(dsig(1)+dsig(3)))/e_vec(i)
          deps_plas(3) = deps(i,3) -
     &                   (dsig(3)-nu_vec(i)*(dsig(2)+dsig(1)))/e_vec(i)

          deps_plas(4) = deps(i,4) - dsig(4)/shear_mod_vec(i)
          deps_plas(5) = deps(i,5) - dsig(5)/shear_mod_vec(i)
          deps_plas(6) = deps(i,6) - dsig(6)/shear_mod_vec(i)
          stress_np1(i,8) = stress_n(i,8) + half * (
     &        deps_plas(1) * (stress_np1(i,1) + stress_n(i,1))
     &   +    deps_plas(2) * (stress_np1(i,2) + stress_n(i,2))
     &   +    deps_plas(3) * (stress_np1(i,3) + stress_n(i,3))
     &   +    deps_plas(4) * (stress_np1(i,4) + stress_n(i,4))
     &   +    deps_plas(5) * (stress_np1(i,5) + stress_n(i,5))
     &   +    deps_plas(6) * (stress_np1(i,6) + stress_n(i,6)))
          factor1 = ( deps_plas(1) - deps_plas(2) )**2 +
     &              ( deps_plas(1) - deps_plas(3) )**2 +
     &              ( deps_plas(3) - deps_plas(2) )**2
          factor2 = deps_plas(4)**2 + deps_plas(5)**2 +
     &              deps_plas(6)**2
          deps_plas_bar = (root2/three)*sqrt( factor1 +
     &                    (three/two)*factor2 )
          stress_np1(i,9) = stress_n(i,9)+deps_plas_bar
        end if
       end do
c
       return
c
       end
c *******************************************************************
c *                                                                 *
c *        material model # 5 -- adv. cyclic plasticity             *
c *                                                                 *
c *           set 3 material model dependent output values          *
c *                                                                 *
c *******************************************************************
c
c
      subroutine oumm05( gpn, mxvl, span, iout, elestr,
     &                   stress, history )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &  gpn, mxvl, span, iout
c
c
      double precision
     & stress(mxvl,*), elestr(mxvl,*), history(mxvl,*)
c
c               description of parameters
c               -------------------------
c
c     gpn               : gauss point number being processed for block
c     mxvl              : maximum no. elements per block
c     span              : number of elements in current block
c     iout              : write messages to this device number
c     stress            : current stresses for all
c                         elements in block for this gauss point
c (*) elestr            : stresses to be output for elements
c     history           : current history values for all
c                         elements in block for this gauss point
c
c    (*)  values to be updated by this material model
c
c
c   stress ordering                elestr ordering
c     (1) sig-xx                   (1) sig-xx
c     (2) sig-yy                   (2) sig-yy
c     (3) sig-zz                   (3) sig-zz
c     (4) tau-xy                   (4) tau-xy
c     (5) tau-yz                   (5) tau-yz
c     (6) tau-xz                   (6) tau-xz
c     (7) total work density       (7) total work density
c                                  (8) mises equiv. stress
c                             (*)  (9) mat_val1
c                             (*) (10) mat_val2
c                             (*) (11) mat_val3
c
c  NOTE:  do not modify "stress" array or columns 1-8 of the "elestr"
c         array. only modify columns 9-11 of "elestr". These are the
c         3 "material model" dependent values output in the stress
c         values for a gauss point. See Section 2.12 of manual and
c         description of each material model. The output labels for
c         columns 9-11 are "c1", "c2", "c3".
c
       integer i, iword(2)
c
      double precision
     &      dword, one, three
c
       equivalence(dword, iword)
c
       data one, three / 1.0, 3.0 /
c
c
c   mat_va11 = c1 = accumulated plastic strain as used in material model
c                   (differs from other definitions of total plastic
c                    strain by a factor of sqrt(3/2))
c   mat_va11 = c2 = current radius of yield surface
c   mat_va11 = c3 = state flag
c                 = 1 for active plastic loading
c                 = 3 otherwise
c
c
       do i= 1, span
         elestr(i,9)  = history(i,3)
         elestr(i,10) = history(i,2)
         elestr(i,11) = three
         dword = history(i,4)
         if (iword(i) .eq. 1 ) then
           elestr(i,11) = one
         end if
       end do
c
c
       return
       end


c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm05_set_sizes                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 12/14/14 rhd                *
c     *                                                              *
c     *    called by warp3d for each material model to obtain        *
c     *    various sizes of data for the model                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm05_set_sizes( info_vector )
      dimension info_vector(*)
c
c        set infor_data
c
c         1        number of history values per integration
c                  point. Abaqus calls these "statev". Values
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
      info_vector(1) = 14
      info_vector(2) = 21
      info_vector(3) = 0
      info_vector(4) = 0
c
      return
      end
c            dummy routines for model not yet supporting
c            states output
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm05_states_values                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 1/3/2015 (rhd))                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm05_states_values( itype, elem_states_output,
     &                                 nrow_states, num_states  )
      use global_data ! old common.main
c
c                       access some global data structures
c
      use elem_block_data, only: history_blocks, history_blk_list
      use main_data, only: elems_to_blocks, cohesive_ele_types
c
      implicit integer (a-z)
c
c                       parameters
c
      integer :: nrow_states, itype, num_states
      double precision :: elem_states_output(nrow_states,*)
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm05_states_labels                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 1/11/2015 (rhd)                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm05_states_labels( size_state,
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
      num_states = 0
      num_comment_lines = 0
      state_labels(1) = "..."
      state_descriptors(1) = "...."
c
      return
      end

