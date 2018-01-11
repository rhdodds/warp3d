c     ****************************************************************
c     *                                                              *
c     *                      subroutine mm01                         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/21/12                    *
c     *                                                              *
c     *     this subroutine performs the recovery of the             *
c     *     unrotated cauchy stress for the rate-independent         *
c     *     mises material model for a block of similar, non-        *
c     *     conflicting elements. mixed isotropic-kinematic          *
c     *     hardening is supported. the material uniaxial stress-    *
c     *     strain curve is bilinear. material response can be       *
c     *     temperature dependent (e, nu, yield, hardening modulus)  *
c     *     temperature dependent thermal expansion coefficients do  *
c     *     not come into play here for stress updates - thermal     *
c     *     strain increments are computed and subtracted out        *
c     *     before we get here                                       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine mm01( span, felem, gpn, step, iter, ym_n1, nu_n1,
     &                 beta, hprime_n1, lnelas, yld_n1, cgn, cgn1,
     &                 deps, history, history1, rtse, dtemps,
     &                 ym_n, nu_n )
      implicit integer (a-z)
      include 'param_def'
c
c          update the stresses and internal state variables at gauss
c          point "gpn" for "span" elements in the current block.
c          values marked with (*) are updated by this routine
c
c          span     -- number of elements in this block
c          felem    -- actual element number for first element in blk
c          gpn      -- gauss point number being processed on
c                      this call for elements in the block
c          step     -- current load step number
c          iter     -- current (global) netwon iteration number (>1)
c          ym_n1    -- young's modulus for elements
c                      in the block at temperature for step n+1
c          nu_n1    -- poisson's ratio for elements
c                      in the block at temperature for step n+1
c          beta     -- isotropic/kinematic fractional factor for elements
c                      in the block
c         hprime_n1 -- plastic hardening modulus for elements
c                      in the block at temperature for step n+1
c          lnelas   -- logical flag indicating element response
c                      is always linear elastic and temperature
c                      independent
c          yld_n1   -- uniaxial yield stress for elements
c                      in the block at temperature for step n+1
c          cgn      -- cartesian stresses at step n for elements
c                      in the block
c    (*)   cgn1     -- cartesian stresses at end of this iteration for
c                      load step n+1 (output) for elements in the block
c          deps     -- cartesian increments of "mechanical" strain
c                      between n and n+1 for elements in the block
c          history  -- history data at step n for elements in the block
c    (*)   history1 -- history data at step n+1 for elements in the block
c    (*)   rtse     -- "relative" trial elstic stress state for elements
c                      in the block (output). computed and used here and
c                      later by consistent tangent generator (routine cnst1).
c          dtemps   -- temperature change over step at this gp for
c                      each element in block
c          ym_n     -- young's modulus for elements
c                      in the block at temperature for step n
c          nu_n     -- poisson's ratio for elements
c                      in the block at temperature for step n
c
c      Notes:
c      ------
c
c        o   the variable mxvl is defined as a parameter in the
c            include file param_def. it sets the maximum possible
c            number of elements in a block allowed in analyses.
c            most arrays are sized based on this variable
c
c        o   most parameters are vectors of length mxvl
c
c        o   the stresses for elements in the block (cgn, cgn1)
c            are arrays sized mxvl by *. The ordering of terms is
c            xx, yy, zz, xy, yz, xz, energy density, plastic
c            work density, acummulated (incremental) plastic
c            strain
c
c        o   the strain increment (deps) for elements in the block
c            is an array sized mxvl by *. The ordering of terms is
c            xx, yy, zz, gam-xy, gam-yz, gam-xz. The thermal
c            strain contribution has been subtracted before this
c            routine is called.
c
c        o   the rtse array (mxvl by 6) stores the trial elastic stresses
c            at n+1. They are computed here and returned for use
c            later by the consistent tangent routine.
c
c        o   the ordering of terms in the history vectors is
c            described in mm01_set_history
c
c
c                     parameters
c
      double precision
     &     nu_n1(*), beta(*), hprime_n1(*), cgn(mxvl,*),
     &     cgn1(mxvl,*), deps(mxvl,*), ym_n1(*), history(span,*),
     &     yld_n1(*), history1(span,*), rtse(mxvl,*), dtemps(*),
     &     ym_n(*), nu_n(*)
      logical lnelas(*)
c
c                     locals on the stack: these are not dynamically
c                     allocated because both mxvl and nstr
c                     are defined as parameters in the include file
c                     param_def
c
      double precision
     &     yfn(mxvl), mrts(mxvl), alpha_n(mxvl,6), devstr_n1(mxvl,nstr),
     &     shear_mod_n1(mxvl), lk(mxvl), kbar(mxvl), eps_vol_n1(mxvl),
     &     zero
      integer iostat(mxvl), instat(mxvl)
      logical yield(mxvl), trcmat, local_debug, prior_linear(mxvl),
     &        isothermal
c
      data zero / 0.0d00 /
c
c                       get the current output device for any meesages.
c                       set up history if this is an iteration during
c                       load step 1.
c
      call iodevn( local_in, local_out, trcmat, 2 )
      local_debug = .false.
      if ( local_debug ) then
        write(local_out,9000) span, felem, gpn, step, iter
      end if
c
      if ( step .eq. 1 )
     &   call mm01_set_history( history, history1, cgn, yld_n1,
     &                          hprime_n1, span, mxvl, ym_n1,
     &                          nu_n1 )
c
c                       do the basic setup of the trial elastic stress
c                       state, deviators at n, pull backstress at n from
c                       history, etc. to compute elastic trial
c                       state we use e, nu at n+1. evaluate the
c                       material state as elastic or currently plastic.
c
      call mm01_init( span, mxvl, local_out, deps, iostat,
     &                prior_linear, history, cgn,
     &                ym_n1, nu_n1, shear_mod_n1, alpha_n,
     &                rtse, isothermal, dtemps, local_debug,
     &                lk, beta, hprime_n1, yld_n1, kbar,
     &                mrts, yfn, instat, yield, iter, lnelas,
     &                eps_vol_n1, ym_n, nu_n )
c
c                       compute updated yield surface size, updated
c                       backstresses and save in history. compute
c                       deviators of updated stress stateat n+1. if all
c                       elements at this gauss point have no temperature
c                       change over the step, use the faster
c                       isothermal update procedure. these update
c                       routines use their own loops over span and skip
c                       linear elastic elements
c
c      if ( isothermal ) then
c           call mm01_simple1( span, mxvl, history, history1,
c     &                        kbar, mrts, shear_mod_n1,
c     &                        hprime_n1, beta, rtse, devstr_n1,
c     &                        yield, local_debug )
c      else
          call mm01_general( span, mxvl, history, history1,
     &                       kbar, mrts, shear_mod_n1,
     &                       hprime_n1, beta, rtse, devstr_n1,
     &                       yield, lk, local_debug, dtemps,
     &                       local_out )
c      end if
c
c                       update elements that are linear elastic at this
c                       point. note: we save the updated yield
c                       surface size and updated back stresses at n+1
c                       to reflect any changes due to temperature.
c                       add the backstresses at n to rtse to define
c                       updated deviators of the trial elastic stress
c                       at n+1
c
      do i = 1, span
         if( yield(i) ) cycle
         history1(i,1)    = zero
         history1(i,2)    = kbar(i)
         if( lnelas(i) ) history1(i,2) = zero
         history1(i,3)    = history(i,3)
         history1(i,5)    = hprime_n1(i)
         history1(i,6:11) = alpha_n(i,1:6) * lk(i)
         devstr_n1(i,1:6) = rtse(i,1:6) + alpha_n(i,1:6)
      end do
c
c                       compute the total updated stresses from their
c                       deviator values at state (n+1) and the
c                       mean stress (linear elastic) contribution.
c                       save the state variable, the elastic
c                       modulus and poisson's ratio at n+1 in
c                       updated history vector.
c                       calculate the energy density from a
c                       trapezoidal numerical integration of
c                       increments of strain and average stresses
c
      call mm01_sig_final( span, mxvl, cgn, cgn1,
     &                     ym_n1, nu_n1, shear_mod_n1,
     &                     devstr_n1, deps, instat, history1,
     &                     eps_vol_n1, local_debug )
c
c                       update the plastic
c                       work density for elements in the block.
c
      call mm01_plastic_work( span, mxvl, cgn, cgn1, yield, deps,
     &                              nu_n1, ym_n1, shear_mod_n1 )
c
      return
 9000 format( '>> debugging in mm01. ',
     &  /, 10x,'span, felem, gpn, step, iter: ',5i6 )
c
      end
c     ****************************************************************
c     *                                                              *
c     *              subroutine mm01_set_history                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 10/1/00                     *
c     *                                                              *
c     *    for step 1, initialize the history data for this gauss pt *
c     *    for all elements in the block                             *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm01_set_history( history, history1, cgn, sigma_o,
     &                             hprime, span, mxvl, ym_n1, nu_n1 )
      implicit integer (a-z)
c
c                     parameters
c
      double precision
     &     hprime(*), cgn(mxvl,*), history(span,*), sigma_o(*),
     &     history1(span,*), ym_n1(*), nu_n1(*)
c
c                     locals
c
      double precision
     &   kn, dword
      integer iword(2)
      equivalence ( dword, iword )
c
c                     numerical constants
c
      double precision
     &   root3, zero
c
      data zero, root3 / 0.0d00, 1.73205080756888d00 /
c
c           set up the history vector for the material.
c           (1) -- updated estimate for lamda * deltat
c                  over step (used by consisten tangent
c                  generator routine)
c           (2) -- equivalent (shear) stress, k, that sets radius
c                  of the yield surface. any amount ofisotropic
c                  hardening after yield changes value
c                  as can temperature dependent flow
c                  properties.
c           (3) -- accumulated (uniaxial) plastic strain
c           (4) -- state (integer)
c           (5) -- hprime at n on entry, at n+1 after
c                  updating (can vary with temperature).
c           (6)-(11) -- back stress for kinematic hardening
c
      iword(1) = 3
      iword(2) = 0
      do i = 1, span
         kn           = sigma_o(i) / root3
         history(i,1) = zero
         history(i,2) = kn
         history(i,3) = zero
         history(i,4) = dword
         history(i,5) = hprime(i)
         history(i,6) = zero
         history(i,7) = zero
         history(i,8) = zero
         history(i,9) = zero
         history(i,10) = zero
         history(i,11) = zero
c
         history1(i,1) = zero
         history1(i,2) = kn
         history1(i,3) = zero
         history1(i,4) = dword
         history1(i,5) = hprime(i)
         history1(i,6) = zero
         history1(i,7) = zero
         history1(i,8) = zero
         history1(i,9) = zero
         history1(i,10) = zero
         history1(i,11) = zero
c
         cgn(i,8)       = zero
         cgn(i,9)       = zero
c
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                    subroutine mm01_init                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 10/1/00                     *
c     *                                                              *
c     *    basic set up of trial elastic state at n+1, pull terms    *
c     *    history at n, stress deviators at n, etc. for all         *
c     *    elements in the block                                     *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm01_init( span, mxvl, iout, deps, iostat,
     &                      prior_linear, history, cgn,
     &                      ym_n1, nu_n1, shear_mod, alpha_n,
     &                      rtse, isothermal, dtemps, debug,
     &                      lk, beta, hprime_n1, yld_n1,
     &                      kbar, mrts, yf, instat, yield,
     &                      iter, lnelas, eps_vol_n1, ym_n, nu_nn )
      implicit none
c
c                     parameters
c
      integer span, mxvl, iostat(*), instat(*), iter, iout
      logical   prior_linear(*), isothermal, debug, yield(*),
     &          lnelas(*)
      double precision
     & deps(mxvl,*), history(span,*), cgn(mxvl,*), ym_n1(*),
     & shear_mod(*), alpha_n(mxvl,*), nu_n1(*), rtse(mxvl,*), dtemps(*),
     & lk(*), beta(*), hprime_n1(*), yld_n1(*), kbar(*), mrts(*), yf(*),
     & eps_vol_n1(*), ym_n(*), nu_nn(*)
c
c                     locals
c
      double precision
     & dword, three, eps_mean, de1, de2, de3, de4, de5, de6,
     & one, two, temper_tol, htol, deps_vol,
     & hbari_np1, hbark_np1, hbark_n, root3, root2, yld_tol,
     & e_n, nu_n, g_n, een1, een2, een3, een4, een5, een6, e1, e2,
     & e3, e4, e5, e6, e, nu, devstr_n1_elas(6), eps_mean_n, zero
c
      integer iword(2), i
      equivalence ( dword, iword )
      data one, two, three, temper_tol, htol, root3, root2, yld_tol,
     &     zero
     & / 1.0d00, 2.0d00, 3.0d00, 0.00001d00, 0.000001d00,
     &   1.7320508075688d00,
     &   1.414213562373095d00, 0.0000001d00, 0.0d00 /
c
c
c                       do the basic setup of the trial elastic stress
c                       state, check for yielding, set flags, etc.
c                       see mm01_init_history for map of history
c                       vector for a gauss point.
c
c                       we also track if the element block is isothermal
c                       over step (enables simpler update process
c                       later)
c
      if ( debug ) write(iout,*) ' ... debugging in mm01_init...'
      isothermal = .true.
c
      do i = 1, span
c
         dword           = history(i,4)
         iostat(i)       = iword(1)
         prior_linear(i) = iostat(i) .eq. 3
c
c                       deviatoric components of strain increment over
c                       step.
c
         deps_vol    = deps(i,1) + deps(i,2) + deps(i,3)
         eps_mean    = deps_vol / three
         de1         = deps(i,1) - eps_mean
         de2         = deps(i,2) - eps_mean
         de3         = deps(i,3) - eps_mean
         de4         = deps(i,4)
         de5         = deps(i,5)
         de6         = deps(i,6)
c
c                      elastic components of total strain at start of
c                      step.
c
         e_n         = ym_n(i)
         nu_n        = nu_nn(i)
         g_n         = e_n/two/(one+nu_n)
         een1        = (cgn(i,1)-nu_n*(cgn(i,2)+cgn(i,3)))/e_n
         een2        = (cgn(i,2)-nu_n*(cgn(i,1)+cgn(i,3)))/e_n
         een3        = (cgn(i,3)-nu_n*(cgn(i,1)+cgn(i,2)))/e_n
         een4        = cgn(i,4) / g_n
         een5        = cgn(i,5) / g_n
         een6        = cgn(i,6) / g_n
c
c                      deviatoric components of elastic strain
c                      at start of step plus deviatoric strain
c                      increment over the step. (volumetric term
c                      at n+1 saved for final update operation)
c
         eps_vol_n1(i) = een1 + een2 + een3 + deps_vol
         eps_mean_n   = (een1 + een2 + een3 ) / three
         e1           = (een1 -  eps_mean_n) + de1
         e2           = (een2 -  eps_mean_n) + de2
         e3           = (een3 -  eps_mean_n) + de3
         e4           = een4 + de4
         e5           = een5 + de5
         e6           = een6 + de6
c
c                       compute deviators for trial elastic
c                       stress state. Uses deviatoric elastic strain at n
c                       + the deviatoric strain increment over the step
c                       and temperature dependent moduli at
c                       n+1. this is way to get temperature effects
c                       on modulus and poisson's ratio properly
c                       included.
c
         shear_mod(i)  = ym_n1(i) / (two*(one+nu_n1(i)))
c
         devstr_n1_elas(1) = two * shear_mod(i) * e1
         devstr_n1_elas(2) = two * shear_mod(i) * e2
         devstr_n1_elas(3) = two * shear_mod(i) * e3
         devstr_n1_elas(4) = shear_mod(i) * e4
         devstr_n1_elas(5) = shear_mod(i) * e5
         devstr_n1_elas(6) = shear_mod(i) * e6
c
c                       pull out backstress at start of step
c                       at start of step) for ease of access
c
         alpha_n(i,1) =  history(i,6)
         alpha_n(i,2) =  history(i,7)
         alpha_n(i,3) =  history(i,8)
         alpha_n(i,4) =  history(i,9)
         alpha_n(i,5) =  history(i,10)
         alpha_n(i,6) =  history(i,11)
c
c                       set isotropic and kinematic plastic hardening
c                       moduli at start and end of step (can
c                       be different due to temperature).
c                       set parameter kbar - shear yield stress at
c                       n+1 based on yield stress at n+1, current
c                       plastic strain and isotropic (plastic) hardening
c                       modulus at n+1. this sets the size of yield
c                       surface to check for yielding with trial stress
c                       state. Also used later in stress update
c                       procedures.
c
c                       set the lk factor for temperature dependent
c                       kinematic hardening.
c
         hbari_np1  = beta(i) * hprime_n1(i)
         hbark_np1  = (one - beta(i)) * hprime_n1(i)
         hbark_n    = (one - beta(i)) * history(i,5)
         kbar(i)    = (yld_n1(i) + hbari_np1*history(i,3))/root3
         lk(i)      = one
         if ( abs( hbark_n ) .gt. htol ) lk(i) = hbark_np1 / hbark_n
c
c                       compute deviators of relative trial elastic
c                       stresses at n+1.
c                       note use of the backstress at n updated
c                       to the temperature at n+1. compute norm
c                       of relative trial stress and evaluate yield
c                       criterion.
c
         rtse(i,1) = devstr_n1_elas(1) - alpha_n(i,1)*lk(i)
         rtse(i,2) = devstr_n1_elas(2) - alpha_n(i,2)*lk(i)
         rtse(i,3) = devstr_n1_elas(3) - alpha_n(i,3)*lk(i)
         rtse(i,4) = devstr_n1_elas(4) - alpha_n(i,4)*lk(i)
         rtse(i,5) = devstr_n1_elas(5) - alpha_n(i,5)*lk(i)
         rtse(i,6) = devstr_n1_elas(6) - alpha_n(i,6)*lk(i)
c
         mrts(i) = sqrt( rtse(i,1)**2+rtse(i,2)**2+rtse(i,3)**2 +
     &              two*(rtse(i,4)**2+rtse(i,5)**2+rtse(i,6)**2) )
         yf(i) = mrts(i) - root2 * kbar(i)
c
c                      the relative trial (deviator) stress used
c                      in the actual update procedures does not
c                      reflect the back stress at n updated for
c                      the temperature at n+1.
c
         rtse(i,1) = devstr_n1_elas(1) - alpha_n(i,1)
         rtse(i,2) = devstr_n1_elas(2) - alpha_n(i,2)
         rtse(i,3) = devstr_n1_elas(3) - alpha_n(i,3)
         rtse(i,4) = devstr_n1_elas(4) - alpha_n(i,4)
         rtse(i,5) = devstr_n1_elas(5) - alpha_n(i,5)
         rtse(i,6) = devstr_n1_elas(6) - alpha_n(i,6)
         mrts(i) = sqrt( rtse(i,1)**2+rtse(i,2)**2+rtse(i,3)**2 +
     &              two*(rtse(i,4)**2+rtse(i,5)**2+rtse(i,6)**2) )
c
c                      update the isothermal flag to reflect the
c                      status of this element
c
         if ( abs( dtemps(i) ) .gt. temper_tol )  isothermal = .false.
c
c                      set various logical flags based on element
c                      status determined by the trial elastic state.
c
c                      if the lnelas flag is true, the element response
c                      must always be linear elastic no matter what.
c
c                      if this gauss point for element is yielding,
c                      set flags.
c
c                      state variable:
c                         = 1,  point is actively yielding
c                         = 3,  point is not actively yielding
c
c
         instat(i) = 3
         yield(i)  = .false.
         if( lnelas(i) ) cycle
c
         if( yf(i) .ge. yld_tol*root2*kbar(i) ) then
            yield(i) = .true.
            instat(i) = 1
         end if
c
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                    subroutine mm01_plastic_work              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 10/1/00                     *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm01_plastic_work( span, mxvl, cgn, cgn1, yield, deps,
     &                              nu, ym, shear_mod )
     &
      implicit none
c
c                     parameters
c
      integer span, mxvl
      double precision
     & cgn(mxvl,*), cgn1(mxvl,*), deps(mxvl,*), nu(*), ym(*),
     & shear_mod(*)
      logical yield(*)
c
c                     locals
c
      integer i
      double precision
     & zero, deps_plas_bar, dsig(6), deps_plas(6), half, root2, three,
     & two, factor1, factor2
c
      data zero, half, root2, three, two
     &  / 0.0d00, 0.5d00, 1.414213562373095d00, 3.0d00, 2.0d00 /

c
c                       for plastic points, compute the updated
c                       plastic work density and accumulated
c                       (uniaxial) plastic strain
c
c
      do i = 1, span
       cgn1(i,8)        = cgn(i,8)
       cgn1(i,9)        = cgn(i,9)
       if ( yield(i) ) then
         deps_plas_bar    = zero
         dsig(1) = cgn1(i,1) - cgn(i,1)
         dsig(2) = cgn1(i,2) - cgn(i,2)
         dsig(3) = cgn1(i,3) - cgn(i,3)
         dsig(4) = cgn1(i,4) - cgn(i,4)
         dsig(5) = cgn1(i,5) - cgn(i,5)
         dsig(6) = cgn1(i,6) - cgn(i,6)
         deps_plas(1) = deps(i,1) -
     &                  (dsig(1) - nu(i)*(dsig(2)+dsig(3)))/ym(i)
         deps_plas(2) = deps(i,2) -
     &                  (dsig(2) - nu(i)*(dsig(1)+dsig(3)))/ym(i)
         deps_plas(3) = deps(i,3) -
     &                  (dsig(3) - nu(i)*(dsig(1)+dsig(2)))/ym(i)
         deps_plas(4) = deps(i,4) - dsig(4) / shear_mod(i)
         deps_plas(5) = deps(i,5) - dsig(5) / shear_mod(i)
         deps_plas(6) = deps(i,6) - dsig(6) / shear_mod(i)
         cgn1(i,8) = cgn(i,8) +  half * (
     &       deps_plas(1) * (cgn1(i,1) + cgn(i,1))
     &     + deps_plas(2) * (cgn1(i,2) + cgn(i,2))
     &     + deps_plas(3) * (cgn1(i,3) + cgn(i,3))
     &     + deps_plas(4) * (cgn1(i,4) + cgn(i,4))
     &     + deps_plas(5) * (cgn1(i,5) + cgn(i,5))
     &     + deps_plas(6) * (cgn1(i,6) + cgn(i,6)) )
         factor1 = ( deps_plas(1) - deps_plas(2) )**2  +
     &             ( deps_plas(2) - deps_plas(3) )**2  +
     &             ( deps_plas(1) - deps_plas(3) )**2
         factor2 = deps_plas(4)**2 +  deps_plas(5)**2 +
     &             deps_plas(6)**2
         deps_plas_bar =  (root2/three) * sqrt( factor1 +
     &                    (three/two)*factor2 )
         cgn1(i,9) = cgn(i,9) + deps_plas_bar
        end if
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                    subroutine mm01_sig_final                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 10/1/00                     *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm01_sig_final( span, mxvl, cgn, cgn1,
     &                           ym, nu, shear_mod,
     &                           devstr_n1, deps, instat, history1,
     &                           eps_vol_n1, debug )
      implicit none
c
c                     parameters
c
      integer span, mxvl, instat(*)
      logical debug
      double precision
     & cgn(mxvl,*), cgn1(mxvl,*),
     & ym(*), nu(*), shear_mod(*), devstr_n1(mxvl,*), deps(mxvl,*),
     & history1(span,*), eps_vol_n1(*)
c
c                     locals
c
      integer i, iword(2), k
      double precision
     & sig_mean_np1, one, two, three, half, dword
      equivalence ( dword, iword )
      data one, two, three, half
     &  / 1.0d00, 2.0d00, 3.0d00, 0.5d00 /
c
c                       compute the updated stresses from their
c                       deviator values at state (n+1) and the
c                       (linear elastic) mean stress contribution.
c                       save the state variable, elastic modulus
c                       and poisson's ratio at n+1 in history.
c                       calculate the energy density from a
c                       trapezoidal numerical integration of
c                       increments of strain and average stresses
c
      do i = 1, span
         sig_mean_np1 = eps_vol_n1(i)
     &                  *(three*ym(i)*nu(i)/((one+nu(i))*
     &                  (one-two*nu(i))) + two*shear_mod(i))/three
         cgn1(i,1) = devstr_n1(i,1) + sig_mean_np1
         cgn1(i,2) = devstr_n1(i,2) + sig_mean_np1
         cgn1(i,3) = devstr_n1(i,3) + sig_mean_np1
         cgn1(i,4) = devstr_n1(i,4)
         cgn1(i,5) = devstr_n1(i,5)
         cgn1(i,6) = devstr_n1(i,6)
         cgn1(i,7) = cgn(i,7) + half * (
     &       deps(i,1) * (cgn1(i,1) + cgn(i,1))
     &     + deps(i,2) * (cgn1(i,2) + cgn(i,2))
     &     + deps(i,3) * (cgn1(i,3) + cgn(i,3))
     &     + deps(i,4) * (cgn1(i,4) + cgn(i,4))
     &     + deps(i,5) * (cgn1(i,5) + cgn(i,5))
     &     + deps(i,6) * (cgn1(i,6) + cgn(i,6)) )
c
         iword(1)       = instat(i)
         history1(i,4)  = dword
      end do
c
      if( debug ) then
          write(*,9100) 2, deps(2,1:6),cgn(2,1:6),cgn1(2,1:6)
      end if
 9100 format(i3,2x,6e14.6,/,5x,6e14.6,/,5x,6e14.6 )


c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm01_simple1                      *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *                   last modified: 10/1/00                     *
c     *                                                              *
c     *      stress update procedure when all elements in block      *
c     *      are isothermal over the load step                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm01_simple1( span, mxvl, history, history1,
     &                         kbar, mrts, shear_mod_n1,
     &                         hprime_n1, beta, rtse, devstr_n1,
     &                         yield, debug )
      implicit none
c
c                     parameters
c
      integer  mxvl, span
      logical debug, yield(*)
      double precision
     & history(span,*), history1(span,*), kbar(*),
     & mrts(*), shear_mod_n1(*), hprime_n1(*), beta(*),
     & rtse(mxvl,*), devstr_n1(mxvl,*)
c
c                     locals
c
      integer j, i
      double precision
     & lambda_deltat, k_np1, hbark_np1, const1, hbari_np1, const2
c
c                     numerical constants
c
      double precision
     & root2, twthrd, one, two, three, root23, root3
      data root2, twthrd, one, two, three, root23
     & / 1.414213562373095d00, 0.666666666666667d00, 1.0d00,
     &   2.0d00, 3.0d00, 0.816496580927d00 /
c
c                       see mm01_init_history for map of history
c                       vector for a gauss point.
c
      do i = 1, span
      if ( .not. yield(i) ) cycle
c
c                       compute isotropic and kinematic plastic
c                       moduli. get the plastic strain multiplier.
c                       update size of yield surface caused by
c                       isotropic hardening (if any).
c
      hbari_np1 = beta(i) * hprime_n1(i)
      hbark_np1 = (one - beta(i)) * hprime_n1(i)
c
      lambda_deltat = ( mrts(i) - root2*kbar(i) ) /
     &                ( twthrd*(three*shear_mod_n1(i) +
     &                  hprime_n1(i)) )
c
      k_np1 = kbar(i) + (root2/three) * hbari_np1 * lambda_deltat
c
c                       update scalars in history for element
c
      history1(i,1)   = lambda_deltat
      history1(i,2)   = k_np1
      history1(i,3)   = history(i,3) + lambda_deltat * root23
      history1(i,5)   = hprime_n1(i)
c
c                 updated backstresses and compute deviators for the
c                 updated stress state. updated backstresses
c                 saved into updated history
c
      const1 = twthrd * hbark_np1 * lambda_deltat / mrts(i)
      const2 = root2 * k_np1 / mrts(i)
c
      history1(i,6 ) = history(i,6 )  + const1 * rtse(i,1)
      history1(i,7 ) = history(i,7 )  + const1 * rtse(i,2)
      history1(i,8 ) = history(i,8 )  + const1 * rtse(i,3)
      history1(i,9 ) = history(i,9 )  + const1 * rtse(i,4)
      history1(i,10) = history(i,10)  + const1 * rtse(i,5)
      history1(i,11) = history(i,11)  + const1 * rtse(i,6)
c
      devstr_n1(i,1) = history1(i,6 ) + const2 * rtse(i,1)
      devstr_n1(i,2) = history1(i,7 ) + const2 * rtse(i,2)
      devstr_n1(i,3) = history1(i,8 ) + const2 * rtse(i,3)
      devstr_n1(i,4) = history1(i,9 ) + const2 * rtse(i,4)
      devstr_n1(i,5) = history1(i,10) + const2 * rtse(i,5)
      devstr_n1(i,6) = history1(i,11) + const2 * rtse(i,6)
c
      end do
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm01_general                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 10/1/00                     *
c     *                                                              *
c     *      stress update procedure when one or more elements       *
c     *      in the block undergo a temperature change over step.    *
c     *      the update procdedure handles additional terms in       *
c     *      isotropic and kinematic hardening due to temperature    *
c     *      changes over the step                                   *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm01_general( span, mxvl, history, history1,
     &                         kbar, mrts, shear_mod_n1,
     &                         hprime_n1, beta, rtse, devstr_n1,
     &                         yield, lk, debug, dtemps, iout )
      implicit none
c
c                     parameters
c
      integer  mxvl, span, iout
      logical debug, yield(*)
      double precision
     & history(span,*), history1(span,*), kbar(*),
     & mrts(*), shear_mod_n1(*), hprime_n1(*), beta(*),
     & rtse(mxvl,*), devstr_n1(mxvl,*), lk(*), dtemps(*)
c
c                     locals
c
      integer j, i
      double precision
     & lambda_deltat, hbark_np1, const1, term1, term2,
     & hbari_np1, const2, a, b, d, vbar, wbar, vwbar,
     & t1, t2, t3, hbark_n, discr, qroot1, qroot2, hbari_n
c
c                     constants
c
      double precision
     & root2, twthrd, one, two, three, root23, four, zero, root2o3
      data root2, twthrd, one, two, three, root23, four, zero,
     &     root2o3
     & / 1.414213562373095d00, 0.666666666666667d00, 1.0d00,
     &   2.0d00, 3.0d00,  0.816496580927d00, 4.0d00,
     &   0.0d00, 0.4714045207910d00 /
c
c                       see mm01_init_history for map of history
c                       vector for a gauss point.
c
      if ( debug ) write(iout,*) ' .. inside mm01_general'
      do i = 1, span
      if ( .not. yield(i) ) cycle
c
c                       compute isotropic and kinematic plastic
c                       moduli at start and end of step.
c
      hbari_np1 = beta(i) * hprime_n1(i)
      hbark_np1 = (one - beta(i)) * hprime_n1(i)
      hbari_n   = beta(i) * history(i,5)
      hbark_n   = (one - beta(i)) * history(i,5)
c
c                       set up terms of the quadratic equation
c                       to solve for the plastic multiplier
c                       lambda * deltat
c
c                       a) coefficients for the 3 tensors
c                          denoted u, v, w in writeup
c                       b) tensor products: v dot v, v dot w,
c                          w dot w. these are vbar, vwbar and wbar
c                       c) t1, t2, t3 are the resulting coefficients
c                          of terms for the quadratic equation.
c                          t1 multiplies lambda bar ** 2,
c                          t1 multiplies lambda bar
c                          t3 is the constant
c
      a = ( two*shear_mod_n1(i) + twthrd*hbark_np1 ) / mrts(i)
      b = one - lk(i)
      d = root2o3 * hbari_np1
c
      vbar  = mrts(i) * mrts(i)
      wbar  = history(i,6)**2 + history(i,7)**2 + history(i,8)**2 +
     &        two*( history(i,9)**2 + history(i,10)**2 +
     &        history(1,11)**2 )
      vwbar = rtse(i,1)*history(i,6) + rtse(i,2)*history(i,7)
     &        + rtse(i,3)*history(i,8) +
     &        two*( rtse(i,4)*history(i,9) + rtse(i,5)*history(i,10) +
     &        rtse(i,6)*history(i,11) )
c
      t1 = a*a*vbar - two*d*d
      t2 = -four*d*kbar(i) - two*a*vbar - two*a*b*vwbar
      t3 = b*b*wbar + two*b*vwbar + vbar -two*kbar(i)*kbar(i)
c
c                       compute discriminant of the quadratic.
c                       if it is negative, we hit a wierd case
c                       where the material is linear elastic. issue
c                       a warning, set platic multiplier to zero.
c                       compute the two roots, take the smaller
c                       root to define lambda bar = lambda * deltat
c
c
      discr = t2*t2 - four*t1*t3
      if ( discr .lt. zero ) then
        write(iout,9000)
        if ( debug )
     &     write(iout,9200) i, hbari_np1, hbark_np1, hbari_n, hbark_n,
     &                   a, b, d, vbar, wbar, vwbar, t1, t2, t3,
     &                   discr
        lambda_deltat = zero
      else
        qroot1 = (-t2 + sqrt( discr )) / (two*t1)
        qroot2 = (-t2 - sqrt( discr )) / (two*t1)
        lambda_deltat = min( qroot1, qroot2 )
      end if
c
      lambda_deltat = min( qroot1, qroot2 )
      if ( lambda_deltat .lt. zero ) then
        write(iout,9100)
        call die_abort
      end if
c
c                       update scalars in the history. lambda *deltat
c                       is used by the routine (cnst1) to compute the
c                       consistent tangent modulus. save the updated
c                       equivalent (shear) stress to set new size
c                       of the yield cylinder for any amount
c                       of isotropic hardening.
c                       updated the accumulated equivalent
c                       plastic strain (ebarp)
c
      history1(i,1) = lambda_deltat
      history1(i,2) = kbar(i) + (root2/three)*hbari_np1*lambda_deltat
      history1(i,3) = history(i,3) + lambda_deltat * root23
      history1(i,5) = hprime_n1(i)
c
c                 updated backstresses and compute deviators for the
c                 updated stress state.
c
      const1 = twthrd * hbark_np1 * lambda_deltat / mrts(i)
      const2 = one - lambda_deltat * ( twthrd *  hbark_np1 +
     &         two * shear_mod_n1(i) ) / mrts(i)
      b = one - lk(i)
c
      history1(i,6 ) = lk(i)*history(i,6 )  + const1 * rtse(i,1)
      history1(i,7 ) = lk(i)*history(i,7 )  + const1 * rtse(i,2)
      history1(i,8 ) = lk(i)*history(i,8 )  + const1 * rtse(i,3)
      history1(i,9 ) = lk(i)*history(i,9 )  + const1 * rtse(i,4)
      history1(i,10) = lk(i)*history(i,10)  + const1 * rtse(i,5)
      history1(i,11) = lk(i)*history(i,11)  + const1 * rtse(i,6)

c
      devstr_n1(i,1)  = history1(i,6) + b*history(i,6) +
     &                  const2 * rtse(i,1)
      devstr_n1(i,2)  = history1(i,7) + b*history(i,7) +
     &                  const2 * rtse(i,2)
      devstr_n1(i,3)  = history1(i,8) + b*history(i,8) +
     &                  const2 * rtse(i,3)
      devstr_n1(i,4)  = history1(i,9) + b*history(i,9) +
     &                  const2 * rtse(i,4)
      devstr_n1(i,5)  = history1(i,10) + b*history(i,10) +
     &                  const2 * rtse(i,5)
      devstr_n1(i,6)  = history1(i,11) + b*history(i,11) +
     &                  const2 * rtse(i,6)
c
      end do
c
      return
 9000 format(/,'>> Warning: solution for lambda-deltat',
     & /,      '            (plastic multiplier) < 0.',
     & /,      '            value set to zero',/)
 9100 format( '>> Fatal Error: routine mm01_general. lambda_deltat is',
     & /,     '                negative. Job terminated.',//)
 9200 format(' i, hbari_np1, hbark_np1, hbari_n, hbark_n: ',
     &     i3,4f10.3,
     &/,'  a, b, d, vbar, wbar, vwbar, t1, t2, t3: ',9e14.6,
     &/,'  discr: ',e14.6)

c
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm01_set_sizes                    *
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
      subroutine mm01_set_sizes( info_vector )
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
      info_vector(4) = 9
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm01_states_values                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 1/10/15 (rhd)                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm01_states_values( itype, elem_states_output,
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
           call mm01_states_values_a
           elem_states_output(1:nrow_states,relem) =
     &                one_elem_states(1:nrow_states)
        end do
      else
        relem = elnum + 1 - felem
        one_elem_states(1:nrow_states) = zero
        call mm01_states_values_a
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
c     *                 subroutine mm01_states_values_a              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/10/15 (rhd)              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm01_states_values_a
c
      implicit none
c
c                       locals
c
      integer :: ipt
      double precision ::
     & epspls, kbar, state, back_stress(6), dword
       integer :: iword(2)
       equivalence ( dword, iword )
c
      epspls =  zero
      kbar = zero
      state  = zero
      back_stress(1:6) = zero
c
      do ipt = 1, int_points
        kbar   = kbar + history_dump(2,ipt,relem)
        epspls = epspls + history_dump(3,ipt,relem)
        dword  = history_dump(4,ipt,relem)
        state  = state + dble(iword(1))
        back_stress(1:6) = back_stress(1:6) +
     &                     history_dump(6:11,ipt,relem)
      end do
c
      one_elem_states(1) = kbar / dble(int_points)
      one_elem_states(2) = epspls / dble(int_points)
      one_elem_states(3) = state  / dble(int_points)
      one_elem_states(4:9) = back_stress(1:6) / dble(int_points)
c
      return
c
      end subroutine mm01_states_values_a
      end subroutine mm01_states_values
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm01_states_labels                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 1/11/2015 (rhd)                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm01_states_labels( size_state,
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
      logical, save :: do_print = .false.
c
      num_states = 9
      state_labels(1) = "kbar"
      state_labels(2) = "epspls"
      state_labels(3) = "status"
      state_labels(4) = "alphaxx"
      state_labels(5) = "alphayy"
      state_labels(6) = "alphazz"
      state_labels(7) = "alphaxy"
      state_labels(8) = "alphayz"
      state_labels(9) = "alphaxz"
c
      state_descriptors(1) = "Radius of yield cylinder"
      state_descriptors(2) = "Plastic strain"
      state_descriptors(3) = "=1 active yield, 3=not active"
      state_descriptors(4:9) = "Backstress"
c
      num_comment_lines = 0
c
      if( do_print ) then
        do i = 1, 3
          write(out,9010) i, state_labels(i), state_descriptors(i)
        end do
        do_print = .false.
      end if
c
      return
 9010 format(2x,i3,2x,a8,2x,a)
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine cnst1                        *
c     *                                                              *
c     *                       written by : bh, rhd                   *
c     *                                                              *
c     *                   last modified : 12/21/2015 rhd             *
c     *                                                              *
c     *      computes the consistent tangent operator matrix for     *
c     *      bilinear material model                                 *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine cnst1( span, cep, rtsg, nu, e, kn1, hprime,
     &                  beta, ldt, dstates, felem, iout )
      implicit none
      include 'param_def'
c
c                       parameter declarations
c
      integer :: span, felem, iout
      double precision ::
     &   cep(mxvl,nstr,*), rtsg(mxvl,*), e(*), nu(*), kn1(*),
     &   hprime(*), beta(*), ldt(*), dstates(*)
c
c                      locals
c
      integer :: i, j, iestate, iword(2)
      double precision ::
     &   bb(mxvl), albar(mxvl),
     &   thbar(mxvl), c1(mxvl), c2(mxvl), c3(mxvl), c4(mxvl), g(mxvl),
     &   l(mxvl), k(mxvl), mrtsq(mxvl), gambar(mxvl), gamma(mxvl),
     &   gbar(mxvl), root2, zero, one, two, three, dword
c
      logical :: yield(mxvl), local_debug
      equivalence ( dword, iword )
c
      data zero, one, two, three, root2 / 0.0d0, 1.0d0, 2.0d0, 3.0d0,
     &      1.414213562d0/
c
      local_debug = .false.
c
      do i = 1, span
       dword   = dstates(i)
       iestate = iword(1)
       if( iestate .eq. 1 ) then
          yield(i) = .true.
       else
          yield(i) = .false.
       end if
      end do
c
      if( local_debug ) then
        write(iout,*) '>>> yield flags...'
        write(iout,*) i, (yield(j),j=1,span)
        do i = 1, span
         write(iout,fmt='(2x,i4,6f10.3)') i+felem-1,(rtsg(i,j),j=1,6)
         write(iout,fmt='(2x,i4,3f10.3)') i+felem-1, e(i),
     &                                 kn1(i), hprime(i)
        end do
      end if
c
      do i = 1, span
       if( yield(i) ) cycle
       cep(i,1,4) = zero
       cep(i,1,5) = zero
       cep(i,1,6) = zero
       cep(i,2,4) = zero
       cep(i,2,5) = zero
       cep(i,2,6) = zero
       cep(i,3,4) = zero
       cep(i,3,5) = zero
       cep(i,3,6) = zero
       cep(i,4,1) = zero
       cep(i,4,2) = zero
       cep(i,4,3) = zero
       cep(i,4,5) = zero
       cep(i,4,6) = zero
       cep(i,5,1) = zero
       cep(i,5,2) = zero
       cep(i,5,3) = zero
       cep(i,5,4) = zero
       cep(i,5,6) = zero
       cep(i,6,1) = zero
       cep(i,6,2) = zero
       cep(i,6,3) = zero
       cep(i,6,4) = zero
       cep(i,6,5) = zero
       c1(i) = (e(i)/((one+nu(i))*(one-two*nu(i))))
       c2(i) = (one-nu(i))*c1(i)
       c3(i) = ((one-two*nu(i))/two)*c1(i)
       c4(i) = nu(i)*c1(i)
       cep(i,1,1) = c2(i)
       cep(i,2,2) = c2(i)
       cep(i,3,3) = c2(i)
       cep(i,4,4) = c3(i)
       cep(i,5,5) = c3(i)
       cep(i,6,6) = c3(i)
       cep(i,1,2) = c4(i)
       cep(i,1,3) = c4(i)
       cep(i,2,1) = c4(i)
       cep(i,3,1) = c4(i)
       cep(i,2,3) = c4(i)
       cep(i,3,2) = c4(i)
      end do
c
      do i = 1, span
       if( .not. yield(i) ) cycle
       g(i) = e(i)/(two*(one+nu(i)))
       l(i) = (e(i)*nu(i))/((one+nu(i))*(one-two*nu(i)))
       k(i) = (three*l(i)+two*g(i))/three
       mrtsq(i) = rtsg(i,1)**2+rtsg(i,2)**2+rtsg(i,3)**2+two*
     &                (rtsg(i,4)**2+rtsg(i,5)**2+rtsg(i,6)**2)
       bb(i) = (root2*kn1(i)+(two/three)*(one-beta(i))*hprime(i)
     &              *ldt(i))/sqrt(mrtsq(i))
       gamma(i) = one/(one+hprime(i)/(three*g(i)))
       gambar(i) =  gamma(i)-one+bb(i)
       gbar(i) = g(i)*bb(i)
       albar(i) = k(i)- two*gbar(i)/three
       thbar(i) = two*g(i)*gambar(i)
       cep(i,1,1) =  (albar(i)+two*gbar(i)-thbar(i)*(rtsg(i,1)**2)/
     &                    mrtsq(i))
       cep(i,2,2) =  (albar(i)+two*gbar(i)-thbar(i)*(rtsg(i,2)**2)/
     &                    mrtsq(i))
       cep(i,3,3) =  (albar(i)+two*gbar(i)-thbar(i)*(rtsg(i,3)**2)/
     &                    mrtsq(i))
       cep(i,4,4) =  (gbar(i)-thbar(i)*(rtsg(i,4)**2)/mrtsq(i))
       cep(i,5,5) =  (gbar(i)-thbar(i)*(rtsg(i,5)**2)/mrtsq(i))
       cep(i,6,6) =  (gbar(i)-thbar(i)*(rtsg(i,6)**2)/mrtsq(i))
       cep(i,2,1) =  (albar(i)-thbar(i)*rtsg(i,1)*rtsg(i,2)/
     &                    mrtsq(i))
       cep(i,3,1) =  (albar(i)-thbar(i)*rtsg(i,1)*rtsg(i,3)/
     &                    mrtsq(i))
       cep(i,4,1) = -(thbar(i)*rtsg(i,1)*rtsg(i,4)/mrtsq(i))
       cep(i,5,1) = -(thbar(i)*rtsg(i,1)*rtsg(i,5)/mrtsq(i))
       cep(i,6,1) = -(thbar(i)*rtsg(i,1)*rtsg(i,6)/mrtsq(i))
       cep(i,3,2) =  (albar(i)-thbar(i)*rtsg(i,3)*rtsg(i,2)/
     &                    mrtsq(i))
       cep(i,4,2) = -(thbar(i)*rtsg(i,2)*rtsg(i,4)/mrtsq(i))
       cep(i,5,2) = -(thbar(i)*rtsg(i,2)*rtsg(i,5)/mrtsq(i))
       cep(i,6,2) = -(thbar(i)*rtsg(i,2)*rtsg(i,6)/mrtsq(i))
       cep(i,4,3) = -(thbar(i)*rtsg(i,3)*rtsg(i,4)/mrtsq(i))
       cep(i,5,3) = -(thbar(i)*rtsg(i,3)*rtsg(i,5)/mrtsq(i))
       cep(i,6,3) = -(thbar(i)*rtsg(i,3)*rtsg(i,6)/mrtsq(i))
       cep(i,5,4) = -(thbar(i)*rtsg(i,4)*rtsg(i,5)/mrtsq(i))
       cep(i,6,4) = -(thbar(i)*rtsg(i,4)*rtsg(i,6)/mrtsq(i))
       cep(i,6,5) = -(thbar(i)*rtsg(i,5)*rtsg(i,6)/mrtsq(i))
       cep(i,1,2) = cep(i,2,1)
       cep(i,1,3) = cep(i,3,1)
       cep(i,1,4) = cep(i,4,1)
       cep(i,1,5) = cep(i,5,1)
       cep(i,1,6) = cep(i,6,1)
       cep(i,2,3) = cep(i,3,2)
       cep(i,2,4) = cep(i,4,2)
       cep(i,2,5) = cep(i,5,2)
       cep(i,2,6) = cep(i,6,2)
       cep(i,3,4) = cep(i,4,3)
       cep(i,3,5) = cep(i,5,3)
       cep(i,3,6) = cep(i,6,3)
       cep(i,4,5) = cep(i,5,4)
       cep(i,4,6) = cep(i,6,4)
       cep(i,5,6) = cep(i,6,5)
      end do
c
      return
      end
