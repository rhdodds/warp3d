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
c         we use a #dbl or #sgl at beginning of lines which
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
#dbl      implicit double precision (a-h,o-z)
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
#r60        sequence
        integer :: iter, abs_element, relem, ipoint, iout
        logical :: allow_cut, segmental, power_law,
     &             rate_depend_segmental, signal_flag, cut_step
#sgl        real :: dtime, step_scale_fact
#dbl        double precision :: dtime, step_scale_fact
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
#dbl      double precision
#sgl      real
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
#dbl      double precision
#sgl      real
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
      history1(relem,10,gpn) = zero
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
#sgl          nsubin    = max( int(one+ptone*equiv_eps/eps_o), 1 )
#dbl          nsubin    = max( idint(one+ptone*equiv_eps/eps_o), 1 )
        else
          nsubin = 1
        end if
        if ( rate_depen ) nsubin = 1
#sgl        scale     = one / real(nsubin)
#dbl        scale     = one / dble(nsubin)
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
#dbl      double precision
#sgl      real
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
#dbl      double precision
#sgl      real
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
#dbl      double precision
#sgl      real
     &    newsig(*)
c
c              This routine accesses the variables defined in the
c              calling routine thru the f-90 host by association
c              feature.
c
c                   local variables in this routine
c
#dbl      double precision
#sgl      real
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

#dbl      double precision
#sgl      real
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
#dbl      double precision
#sgl      real
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
c         we use a #dbl or #sgl at beginning of lines which
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
$add param_def
c
c                   parameter declarations
c 
      integer :: element, gpn, iter, iout, span 
#dbl      double precision ::
#sgl      real ::
     & stress_trial(mxvl,*), history(span,*),
     & history1(span,*), cep(mxvl,6,6), e(*), nu(*), q1(*), 
     & q2(*), q3(*), nuc_s_n(*), nuc_e_n(*), nuc_f_n(*)
      logical :: nucleation(*)
c
c                   locally defined array, variables
c
#dbl      double precision ::
#sgl      real ::
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
#dbl      double precision::
#sgl      real ::
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
      dword(1:span) = history1(1:span,6)
c
#dbl      inc_factor = 2
#sgl      inc_factor = 1
      j = 1
@!DIR$ LOOP COUNT MAX=###  
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
@!DIR$ LOOP COUNT MAX=###  
      do  i = 1, span
       if( state(i) .ne. -1 ) cycle
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
@!DIR$ LOOP COUNT MAX=###  
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
        if( iter .le. 1 ) then
          dep(i) = zero
          deq(i) = zero
        end if
      end do
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
@!DIR$ LOOP COUNT MAX=###  
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
@!DIR$ LOOP COUNT MAX=###  
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
@!DIR$ LOOP COUNT MAX=###  
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
@!DIR$ LOOP COUNT MAX=###  
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
@!DIR$ LOOP COUNT MAX=###  
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
@!DIR$ LOOP COUNT MAX=###  
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
@!DIR$ LOOP COUNT MAX=###  
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


