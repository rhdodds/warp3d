c *******************************************************************
c *                                                                 *
c *        material model # 7 -- warp3d hydrogen model              *
c *                                                                 *
c *               written by: yueming liang                         *
c *                                                                 *
c *                last modified: 9/29/2016 rhd                     *
c *                                                                 *
c *                                                                 *
c *     this material model incorporates the equilibrium hydrogen   *
c * softeninging model by Sofronis et al. (sofronis, liang, and     *
c * aravas. hydrogen induced shear localization of the plastic      *
c *     flow in metals and alloys. eur. j. mech. a/solids. 20(2001) *
c *     857-872). the algorithm for updating stresses, strains, and *
c * the current hydrogen concentrations is the modified aravas      *
c * type of algorithm (n. aravas, on the numerical integration      *
c * of a class of pressure-dependent plasticity models. int. j.     *
c * num. method. engng. vol. 24, 1987, 1395-1416).                  *
c *                                                                 *
c *******************************************************************
      subroutine mm07(
     &  step, iter, felem, gpn, mxvl, hist_size, nstrs, nstrn, span,
     &  iout, signal_flag, adaptive_possible, cut_step_size_now,
     &  mm_props, e_vec, tan_e_vec, nu_vec, sigyld_vec,
     &  n_power_vec, trial_elas_stress_np1, stress_n, stress_np1,
     &  deps, history_n, history_np1 )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer ::
     &  step, iter, felem, gpn, mxvl, hist_size, span, iout,
     &  nstrs, nstrn
c
      logical ::
     &   signal_flag, adaptive_possible, cut_step_size_now
c
      double precision ::
     & mm_props(mxvl,10), e_vec(mxvl), tan_e_vec(mxvl), nu_vec(mxvl),
     & sigyld_vec(mxvl), n_power_vec(mxvl), stress_n(mxvl,nstrs),
     & stress_np1(mxvl,nstrs), deps(mxvl,nstrn),
     & trial_elas_stress_np1(mxvl,nstrn), history_n(span,hist_size),
     & history_np1(span,hist_size)
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
c     mm_props          : material parameter values input by user for
c                         each element in block
c     e_vec             : Young's modulus for each element in block
c     tan_e_vec         : constant tangent modulus after yield (may not
c                         be defined) for each element in block
c     nu_vec            : Poisson's ratio for each element in block
c     sigyld_vec        : yield stress for each element in block
c     n_power_vec       : power-law hardening exponent for each element
c                         in block (may not be defined)
c (*) trial_elas_stress_np1 : trial elastic stress vector to be used
c                          later by consistent tangent routine for model
c (**)stress_n          : stresses at start of load step (n) for all
c                         elements in block for this gauss point
c (*) stress_np1        : stresses at end of load step (n+1) for all
c                         elements in block for this gauss point
c     deps              : current estimate of strain increment over the
c                         load step (minus increment of thermal strain)
c (**)history_n         : history values at start of load step (n) for
c                         all elements in block for this gauss point
c (*) history_np1       : history values at end of load step (n+1) for
c                         all elements in block for this gauss point
c
c    (*)  values to be updated by this material model
c    (**) needs to be initialized on step 1
c
c
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
c   mm_props ordering
c (1) alpha, # of H sites per trap, C_T=alpha*theta_t*N_T
c (2) rho0, initial dislocation density, line length/m^3
c (3) gamma, rho=rho0+gamma*epsilon(p), line length/m^3
c (4) a, lattice parameter, m
c (5) W_B, trap binding energy, j/mole
c
c (6) beta, number of NILS per solvent atom, C_L=beta*theta_l*N_L
c (7) V_M, molar volume of host metal (m^3/mole)
c (8) V_H, partial molar volume of H in solution (m^3/mole)
c (9) c_0, initial hydrogen concentration (H/M)
c (10)rksi, softening parameter, sigY(c)=sig0*((rksi-1)*c+1)
c
c   history parameters ordering
c (1) plastic strain
c (2) c, current hydrogen concentration
c (3) etran, transformation strain due to hydrogen
c (4) yflag, indicator of yielding
c                        ( -1 = currently linear elasticity)
c                        (  1 = active yielding            )
c       (5) increment of plastic strain
c
c                   local variables
c                   ---------------
c
      integer :: i,j
      logical :: local_debug
      double precision ::
     & shear_mod, a,b,c, zero, one, two,half,three,onep5,
c
     & str_epd(6), str_epd_1(6), strdir(6),
     & deps_plas(6),stress_np1_1(6),
     & alpha, rho0, gamma, alat, wb, beta, vm, vh, c0, rksi,
     & rgas, rna, temper, rnl, alamda,
     & eps_n, hcon_n, etran_n, eps_0, sigma_0, bk, sigmaY,
     & pe, qe, fac1, tol, yflag, dekk, hcon_np1, eps_np1,
     & p_np1, etran_np1, delastic,aux1, dplastic,n_power,deleps,
     & e11,e22,e33,s11,s22,s33,s44,s55,s66,e44,e55,e66,et,
     & dsdep,dsdc,cl,ct,eps,cnow, thl, tht, dcdskk, dcdep,
     & dntdeps,skk
c
      data zero, one, two/0.0d00, 1.0d00, 2.0d00/
      data half, three, onep5 /0.5d00, 3.0d00, 1.5d00/
      data rgas,rna,temper,tol/8.3144d00,6.0232d23,300.d00,1.0d-6/
      data local_debug / .false. /
c
      if( step .eq. 1 ) then
          do i = 1, span
          if( gpn.eq.1 .and. i.eq.1. and. local_debug) then
            write(iout,*)
            write(iout,*) 'Mechanical properties are:'
            write(iout,*)
            write(iout,*) '  Young-s modulus: ', e_vec(i)
            write(iout,*) '  Poisson-s ratio: ', nu_vec(i)
            write(iout,*) '  Yield point    : ', sigyld_vec(i)
            write(iout,*) '  Hardening ratio: ', n_power_vec(i)
c
            write(iout,*)
            write(iout,*) 'Hydrogen related properties are:'
            write(iout,*)
            write(iout,*) '  alpha = ', mm_props(i,1)
            write(iout,*) '  rho0  = ', mm_props(i,2)
            write(iout,*) '  gamma = ', mm_props(i,3)
            write(iout,*) '  alat  = ', mm_props(i,4)
            write(iout,*) '  wb    = ', mm_props(i,5)
            write(iout,*) '  beta  = ', mm_props(i,6)
            write(iout,*) '  Vm    = ', mm_props(i,7)
            write(iout,*) '  Vh    = ', mm_props(i,8)
            write(iout,*) '  c0    = ', mm_props(i,9)
            write(iout,*) '  rksi  = ', mm_props(i,10)
          end if
            stress_n(i,7:9) = zero ! support user initial stresses now
            history_n(i,1) = zero
            history_n(i,2) = mm_props(i,9)
            history_n(i,3) = zero
            history_n(i,4) = -one
            history_n(i,5:7) = zero
          end do
      end if
c
      do i = 1, span
c
        alpha = mm_props(i,1)
        rho0  = mm_props(i,2)
        gamma = mm_props(i,3)
        alat  = mm_props(i,4)
        wb    = mm_props(i,5)
        beta  = mm_props(i,6)
        vm    = mm_props(i,7)
        vh    = mm_props(i,8)
        c0    = mm_props(i,9)
        rksi  = mm_props(i,10)
c
        alamda= vh/vm
        rnl   = rna/vm
c
        if(step.eq.1) then
           eps = zero
           skk = zero
           call mm07_hydrogen( iout, alpha, rho0, gamma, alat,
     &       wb, beta, vm, vh, c0, rksi, alamda, rnl, rna, rgas,
     &       temper, eps, skk, cnow, thl, tht, dcdskk, dcdep,
     &       dntdeps, cl, ct)
           history_n(i,2) = cnow
        end if
c
        eps_n   = history_n(i,1)
        hcon_n  = history_n(i,2)
        etran_n = history_n(i,3)
c

        sigma_0 = sigyld_vec(i)
        eps_0   = sigma_0 / e_vec(i)
        n_power = n_power_vec(i)
c
        bk = e_vec(i) / ( three * (one - two * nu_vec(i)))
c
        shear_mod = e_vec(i) /
     &         (two * ( one + nu_vec(i) ))
c
        a = two * shear_mod * ( one - nu_vec(i) )
     &         /( one - two * nu_vec(i))
c
        b = two * shear_mod * nu_vec(i)
     &         /( one - two * nu_vec(i))
c
        str_epd(1) = stress_n(i,1) + a * deps(i,1)
     &                       + b * deps(i,2) + b * deps(i,3)
        str_epd(2) = stress_n(i,2) + a * deps(i,2)
     &                       + b * deps(i,1) + b * deps(i,3)
        str_epd(3) = stress_n(i,3) + a * deps(i,3)
     &                       + b * deps(i,1) + b * deps(i,2)
        str_epd(4) = stress_n(i,4) + shear_mod * deps(i,4)
        str_epd(5) = stress_n(i,5) + shear_mod * deps(i,5)
        str_epd(6) = stress_n(i,6) + shear_mod * deps(i,6)
c
        trial_elas_stress_np1(i,1) = str_epd(1)
        trial_elas_stress_np1(i,2) = str_epd(2)
        trial_elas_stress_np1(i,3) = str_epd(3)
        trial_elas_stress_np1(i,4) = str_epd(4)
        trial_elas_stress_np1(i,5) = str_epd(5)
        trial_elas_stress_np1(i,6) = str_epd(6)
c
        pe = (str_epd(1) + str_epd(2) + str_epd(3))/three
c
        str_epd_1(1) = str_epd(1) - pe
        str_epd_1(2) = str_epd(2) - pe
        str_epd_1(3) = str_epd(3) - pe
        str_epd_1(4) = str_epd(4)
        str_epd_1(5) = str_epd(5)
        str_epd_1(6) = str_epd(6)
c
        fac1 = str_epd_1(1) * str_epd_1(1)
     &          + str_epd_1(2) * str_epd_1(2)
     &          + str_epd_1(3) * str_epd_1(3)
     &          + ( str_epd_1(4) * str_epd_1(4)
     &            + str_epd_1(5) * str_epd_1(5)
     &            + str_epd_1(6) * str_epd_1(6) ) * two
c
        qe = sqrt( onep5 * fac1 )
c
        strdir(1) = onep5 * str_epd_1(1) / qe
        strdir(2) = onep5 * str_epd_1(2) / qe
        strdir(3) = onep5 * str_epd_1(3) / qe
        strdir(4) = onep5 * str_epd_1(4) / qe
        strdir(5) = onep5 * str_epd_1(5) / qe
        strdir(6) = onep5 * str_epd_1(6) / qe
c
        call mm07_yield( hcon_n, rksi, eps_n, sigma_0, eps_0,
     &                   n_power, sigmaY, dsdep, dsdc)
c
        if( qe.le.sigmaY ) then
c
          call mm07_elastic( alpha, rho0, gamma, alat, wb, beta,
     &      vm, vh, c0, rksi, alamda, rnl, rna, rgas, temper, bk,
     &      shear_mod, pe, p_np1, hcon_n, hcon_np1, etran_n, eps_n,
     &      etran_np1, adaptive_possible, cut_step_size_now, iout,
     &      cl, ct )
c
          stress_np1(i,1) = str_epd_1(1) + p_np1
          stress_np1(i,2) = str_epd_1(2) + p_np1
          stress_np1(i,3) = str_epd_1(3) + p_np1
          stress_np1(i,4) = str_epd_1(4)
          stress_np1(i,5) = str_epd_1(5)
          stress_np1(i,6) = str_epd_1(6)
c
          history_np1(i,1) = eps_n
          history_np1(i,2) = hcon_np1
          history_np1(i,3) = etran_np1
          history_np1(i,4) = - one
          history_np1(i,5) = zero
          history_np1(i,6) = cl
          history_np1(i,7) = ct

          delastic =
     &        half * ( stress_np1(i,1) + stress_n(i,1) )*deps(i,1) +
     &        half * ( stress_np1(i,2) + stress_n(i,2) )*deps(i,2) +
     &        half * ( stress_np1(i,3) + stress_n(i,3) )*deps(i,3) +
     &        half * ( stress_np1(i,4) + stress_n(i,4) )*deps(i,4) +
     &        half * ( stress_np1(i,5) + stress_n(i,5) )*deps(i,5) +
     &        half * ( stress_np1(i,6) + stress_n(i,6) )*deps(i,6)
          stress_np1(i,7) = stress_n(i,7) +  delastic
          stress_np1(i,8) = zero
          stress_np1(i,9) = zero
        end if
c
        if( qe.gt.sigmaY ) then
c
          call mm07_plastic( alpha, rho0, gamma, alat, wb, beta,
     &       vm, vh, c0, rksi, alamda, rnl, rna, rgas, temper, bk,
     &       shear_mod, pe, p_np1, hcon_n, hcon_np1, etran_n,
     &       etran_np1, eps_n, eps_np1, adaptive_possible,
     &       cut_step_size_now, iout,
     &       sigma_0, eps_0, n_power, qe, cl,ct )
c
          deleps = eps_np1 - eps_n
          aux1 = two * shear_mod * deleps
c
          stress_np1_1(1) = str_epd_1(1) - aux1 * strdir(1)
          stress_np1_1(2) = str_epd_1(2) - aux1 * strdir(2)
          stress_np1_1(3) = str_epd_1(3) - aux1 * strdir(3)
          stress_np1_1(4) = str_epd_1(4) - aux1 * strdir(4)
          stress_np1_1(5) = str_epd_1(5) - aux1 * strdir(5)
          stress_np1_1(6) = str_epd_1(6) - aux1 * strdir(6)
c
          stress_np1(i,1) = stress_np1_1(1) + p_np1
          stress_np1(i,2) = stress_np1_1(2) + p_np1
          stress_np1(i,3) = stress_np1_1(3) + p_np1
          stress_np1(i,4) = stress_np1_1(4)
          stress_np1(i,5) = stress_np1_1(5)
          stress_np1(i,6) = stress_np1_1(6)
c
          history_np1(i,1) = eps_np1
          history_np1(i,2) = hcon_np1
          history_np1(i,3) = etran_np1
          history_np1(i,4) = one
          history_np1(i,5) = deleps
          history_np1(i,6) = cl
          history_np1(i,7) = ct
c
          delastic =
     &       half * ( stress_np1(i,1) + stress_n(i,1) ) *  deps(i,1) +
     &       half * ( stress_np1(i,2) + stress_n(i,2) ) *  deps(i,2) +
     &       half * ( stress_np1(i,3) + stress_n(i,3) ) *  deps(i,3) +
     &       half * ( stress_np1(i,4) + stress_n(i,4) ) *  deps(i,4) +
     &       half * ( stress_np1(i,5) + stress_n(i,5) ) *  deps(i,5) +
     &       half * ( stress_np1(i,6) + stress_n(i,6) ) *  deps(i,6)
          stress_np1(i,7) = stress_n(i,7) +  delastic
          stress_np1(i,8) = zero
          stress_np1(i,9) = zero
c
          deps_plas(1) = deleps * strdir(1)
          deps_plas(2) = deleps * strdir(2)
          deps_plas(3) = deleps * strdir(3)
          deps_plas(4) = deleps * strdir(4)
          deps_plas(5) = deleps * strdir(5)
          deps_plas(6) = deleps * strdir(6)
c
          dplastic =
     &       half * ( stress_np1(i,1) + stress_n(i,1) )*deps_plas(1) +
     &       half * ( stress_np1(i,2) + stress_n(i,2) )*deps_plas(2) +
     &       half * ( stress_np1(i,3) + stress_n(i,3) )*deps_plas(3) +
     &       one * ( stress_np1(i,4) + stress_n(i,4) )*deps_plas(4) +
     &       one * ( stress_np1(i,5) + stress_n(i,5) )*deps_plas(5) +
     &       one * ( stress_np1(i,6) + stress_n(i,6) )*deps_plas(6)
          stress_np1(i,8) = stress_n(i,8) + dplastic
          stress_np1(i,9) = stress_n(i,9) + deleps
c
        end if
      end do  ! on span
c
      return
      end
c *******************************************************************
       subroutine mm07_elastic( alpha, rho0, gamma, alat, wb, beta,
     &    vm, vh, c0, rksi, alamda, rnl, rna, rgas, temper, bk,
     &    shear_mod, pe, p_np1, hcon_n, hcon_np1, etran_n, eps_n,
     &    etran_np1, adaptive_possible, cut_step_size_now, iout,
     &    cl, ct )
       implicit none
c
       integer :: iout
       logical :: adaptive_possible, cut_step_size_now

       double precision ::
     &   alpha, rho0, gamma, alat, wb, beta, vm, vh,
     &   c0, rksi, alamda,
     &   rnl, rna, rgas, temper, bk, shear_mod, pe, p_np1, hcon_n,
     &   hcon_np1, etran_n, eps_n,etran_np1
c
c                  local variables
c
       integer :: j
      double precision ::
     &    zero, one, two, three, tol, cguess,
     &    delc, ddc, detkk, skk, cnow, thl, tht, dcdskk, dcdep,
     &    dntdeps, res, dskkdc, slope, eps, cl, ct,tol1

      data zero, one, two, three /0.0d00, 1.0d00, 2.0d00, 3.0d00/
      data tol,tol1, cguess/1.d-6,1.d-8, 1.d-3/
c
      delc = zero
      ddc  = zero
      do j = 1, 20
         delc = delc + ddc
         hcon_np1 = hcon_n + delc
         detkk = three * delc * alamda /
     &           ( three + etran_n )
         p_np1 = pe - bk * detkk
         skk = three * p_np1
         eps = eps_n
c
         call mm07_hydrogen( iout, alpha, rho0, gamma, alat,
     &      wb, beta, vm, vh, c0, rksi, alamda, rnl, rna, rgas,
     &      temper, eps, skk, cnow, thl, tht, dcdskk, dcdep,
     &      dntdeps, cl, ct )
c
         res = hcon_np1 - cnow
c
         if( abs(res/c0).lt.tol ) goto 1117
c
         dskkdc = - three * bk * alamda
     &          / ( one + etran_n / three )
c
         slope = one - dskkdc * dcdskk
         ddc = - res / slope
      end do
      write(iout,*) '*** no convergence in mm07 for elasticity ****'
      write(iout,*) '....need to cut the step size immediately.......'
      if( adaptive_possible ) cut_step_size_now=.true.

 1117  continue
       etran_np1 = etran_n + delc * alamda
       return
       end
c *******************************************************************
       subroutine  mm07_hydrogen( iout, alpha, rho0, gamma, alat,
     &       wb, beta, vm, vh, c0, rksi, alamda, rnl, rna, rgas,
     &       temper, eps, skk, cnow, thl, tht, dcdskk, dcdep,
     &       dntdeps, cl, ct)
       implicit none
c
       integer :: iout
c
      double precision ::
     & alpha, rho0, gamma, alat,
     & wb, beta, vm, vh, c0, rksi, alamda, rnl, rna, rgas,
     & temper, eps, skk, cnow, thl, tht, dcdskk, dcdep,
     & dntdeps, dthldskk, dthtdskk
c
c                local variables
c
      double precision ::
     & zero, one,two,three,half,fac,
     & aux, rkl, rkt, thl0, cl, ct, rho, epst2, tden,
     & n1,n2,n3,n4,n5,n6,n7,ten, tdenfac,k1,k2,k3
c
      data one, two, three /1.0d00, 2.0d00, 3.0d00/
      data zero , half,fac /0.0d00, 0.5d00,6932799.d0/
      data n1,n2,n3,n4/20.911d0,10.334d0,18.635d0,17.073d0/
      data n5,n6,n7,ten/8.302d0,2.034d0,0.197d0,10.d0/
      data k1,k2,k3/23.26d0,-2.33d0,-5.5d0/
c
      aux = vh /( three * rgas * temper )
      rkl = exp( skk * aux )
      rkt = exp( wb / (rgas * temper))
      thl0 = c0 / beta
      thl = thl0 * rkl / ( one - thl0 + thl0 * rkl )
      tht = thl * rkt / ( one - thl + thl * rkt )
      cl = beta * thl
c
      if( rho0.eq.zero .and .gamma.eq.zero )  then
c
c                 formula by Krom et al, 1999
c
            tden = k1 + k2 * exp(k3*eps)
            tden = exp ( tden * log(ten) )
            tdenfac = k2 * k3 * exp(k3 * eps)
            dntdeps = tden * tdenfac * log(ten)
c
c               formula by Taha and Sofronis
c
c           tden = n1 + n2 * eps
c     +          - n3 * eps * eps
c     +          + n4 * eps * eps * eps
c     +          - n5 * eps * eps * eps * eps
c     +          + n6 * eps * eps * eps * eps * eps
c     +          - n7 * eps * eps * eps * eps * eps * eps
c           tden = exp ( tden * log(ten) )
c           tdenfac =  n2
c     +          - n3 * eps
c     +          + n4 * eps * eps
c     +          - n5 * eps * eps * eps
c     +          + n6 * eps * eps * eps * eps
c     +          - n7 * eps * eps * eps * eps * eps
c           dntdeps = tdenfac * tden * log(ten)
c
      else
           epst2 = eps
           if( eps.gt.half ) then
              epst2 = half
              rho = rho0 + gamma * epst2
              tden = sqrt( two ) * rho / alat
              dntdeps = zero
           else
              rho = rho0 + gamma * epst2
              tden = sqrt( two ) * rho / alat
              dntdeps = sqrt ( two ) * gamma / alat
           endif
      endif
c
      ct = alpha * tht * tden / rnl
      cnow = cl + ct
      if( cnow.gt.one ) then
           write(iout,*) '>>>>>>>FATAL ERROR:'
           write(iout,*) 'c > 1 in mm07_hydrogen'
           stop
      end if
c
      dthldskk = thl*(one-thl0)*aux/(one-thl0+thl0*rkl)
      dthtdskk = dthldskk*rkt/( one -thl+thl*rkt)**two
      dcdskk   = beta*dthldskk+alpha*tden/rnl*dthtdskk
      dcdep    = alpha*tht/rnl*dntdeps
c
      return
      end

c *******************************************************************
       subroutine mm07_plastic( alpha, rho0, gamma, alat, wb, beta,
     &    vm, vh, c0, rksi, alamda, rnl, rna, rgas, temper, bk,
     &    shear_mod, pe, p_np1, hcon_n, hcon_np1, etran_n, etran_np1,
     &    eps_n, eps_np1, adaptive_possible, cut_step_size_now, iout,
     &    sigma_0, eps_0, n_power, qe, cl, ct )
       implicit none
c
       integer :: iout, j
       logical :: adaptive_possible, cut_step_size_now
c
      double precision ::
     & alpha, rho0, gamma, alat, wb, beta,
     & vm, vh, c0, rksi, alamda, rnl, rna, rgas, temper, bk,
     & shear_mod, pe, p_np1, hcon_n, hcon_np1, etran_n, etran_np1,
     & eps_n, eps_np1, sigma_0, eps_0, n_power, qe, cl, ct
c
c                   local variables
c
      double precision ::
     & tol, zero, one, two, three,
     & epsguess, ddeps, sigmaY, dsdep, ff1, slope,ff2,deleps,
     & delc, ddelc, detkk, dsdc, cnow, thl, tht, dcdskk,tol1,
     & dcdep, dntdeps, a11, a12, a21, a22, det1, skk, dskkdc
c
      data zero, one, two, three /0.0d00, 1.0d00, 2.0d00, 3.0d00/
      data tol,tol1, epsguess/1.d-6,1.d-8, 1.d-4/
c
      deleps = epsguess
      delc = zero
      ddeps = zero
      ddelc = zero
      do j = 1, 20
            delc = delc + ddelc
            deleps = deleps + ddeps
            hcon_np1 = hcon_n + delc
            eps_np1 = eps_n + deleps
            detkk = three * delc * alamda /
     &            ( three + etran_n )
            p_np1 = pe - bk * detkk
            skk = three * p_np1
c
            call mm07_yield( hcon_np1, rksi, eps_np1, sigma_0, eps_0,
     &                       n_power, sigmaY, dsdep, dsdc )
c
            ff1 = qe - three * shear_mod * deleps - sigmaY
c
            call mm07_hydrogen( iout, alpha, rho0, gamma, alat,
     &       wb, beta, vm, vh, c0, rksi, alamda, rnl, rna, rgas,
     &       temper, eps_np1, skk, cnow, thl, tht, dcdskk, dcdep,
     &       dntdeps, cl, ct )
c
            ff2 = hcon_np1 - cnow
c
            if(abs(ff1/sigma_0).lt.tol.and.abs(ff2/c0).lt.tol)
     &        goto 1121
c
            a11 = - dsdc
            a12 = - three * shear_mod - dsdep
            dskkdc = - three * bk * alamda
     &             / ( one + etran_n / three )
            a21 = one - dcdskk * dskkdc
            a22 = - alpha / rnl * dntdeps * tht
            det1 = a11 * a22 - a12 * a21
            ddelc = - ( a22 * ff1 - a12 * ff2 ) / det1
            ddeps = - ( a11 * ff2 - a21 * ff1 ) / det1
      end do
      write(iout,*) '*** no convergence in mm07 plasticiy!!!***'
      write(iout,*) '....need to cut the step size immediately...'
      if(adaptive_possible) cut_step_size_now=.true.
c
 1121 continue
c
      etran_np1 = etran_n + delc * alamda
c
      return
      end
c *******************************************************************

       subroutine mm07_yield( hcon, rksi, eps, sigma_0, eps_0,
     &    n_power, sigmaY, dsdep, dsdc)
       implicit none
c
      double precision ::
     & hcon, rksi, eps, sigma_0, eps_0, n_power,
     & sigmaY, dsdep, dsdc
c
c                   local variables
c
      double precision ::
     &    zero, one, aux1, aux2
      data zero, one / 0.d00, 1.d00 /
c
      aux2 = ( rksi - one )* hcon + one
      if( eps.le.zero ) then
          sigmaY = sigma_0 * aux2
          dsdep = sigma_0 / eps_0
      else
          aux1   = one + eps/eps_0
          sigmaY = sigma_0 * aux2 * aux1 ** ( one / n_power )
          dsdep  = sigmaY / n_power / ( eps + eps_0)
          dsdc   = sigma_0*(rksi - one )* aux1 **(one/n_power)
      end if
c
      return
      end
c *******************************************************************
c *                                                                 *
c *        material model # 7 -- adv. mises + hydrogen              *
c *                                                                 *
c *******************************************************************
c
c
      subroutine cnst7(
     &  span, felem, gpn, iter, iout, mxvl, nstrn,
     &  e_vec, nu_vec, sigyld_vec, n_power_vec, mm_props,
     &  trial_elas_stress, history_n,
     &  history_np1, stress_np1, dmat )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &  span, felem, gpn, iter, iout, mxvl, nstrn
c
      double precision
     & mm_props(mxvl,10), e_vec(mxvl), nu_vec(mxvl),
     & trial_elas_stress(mxvl,nstrn), history_n(span,*),
     & history_np1(span,*), dmat(mxvl,nstrn,nstrn),
     & stress_np1(mxvl,nstrn), sigyld_vec(mxvl), n_power_vec(mxvl)
c
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
c     sigyld_vec        : yield stress for each element in block
c     n_power_vec       : power-law hardening exponent
c     trial_elas_stress_np1 : trial elastic stress vector defined by
c                             stress update routine
c                         consistent tangent routine for model
c     history_n         : history values at start of load step (n) for
c                         all elements in block for this gauss point
c     history_np1       : history values at end of load step (n+1) for
c                         all elements in block for this gauss point
c     det_jac_block     : |J| at this gauss point for each element in
c                         block
c     stress_np1        : current estimate of 6 stress components for
c                         end of step (see ordering below)
c     weight            : integration point weight factor
c (*) dmat              : 6x6 (symmetric) tangent (consistent) for
c                         this gauss point for each element of block
c                         (see stress ordering below)
c
c    (*)  values to be updated by this material model
c
c   strain ordering:
c     deps-xx, deps-yy, deps-zz, gamma-xy, gamma-yz, gamma-xz
c
c   stress ordering (at n and n+1):
c     sig-xx, sig-yy, sig-zz, tau-xy, tau-yz, tau-xz
c
      integer :: i,jj,kk
      logical :: local_debug
      double precision ::
     &  zero, one, two,three, quarter, half, onep5,
     &  one_third, two_third,
c
     &  rna,rgas,temper,fact,c1,c2,c3,c4,bk,shear_mod,
     &  alpha,rho0,gamma,alat,wb,beta,vm,vh,c0,rksi,alamda,
     & rnl,hcon,delc,rlam,drlamdc,eps,skk,thl,tht,dthldskk,
     & dthtdskk,dcdskk,dcdep,fac1,fac2,sigmaY,dsdep,dsdc,
     & pe,qe,fac3,fac4,tden,sigma_0,hcon_np1,eps_np1,eps_0,
     & deleps,dntdeps,rbeta,rgamma,n_power,ctol,cnow,fac,
     &  cl, ct, etran,
c
     &  dcde(6),depde(6),delta(6),an(6),str_epd_1(6),ak(6,6)
c
      data zero, one, two, three /0.0d00,1.0d00,2.0d00,3.0d00/
      data quarter, half, onep5 /0.25d00, 0.5d00, 1.5d00/
      data one_third, two_third/0.333333333d00, 0.666666667d00/
      data rgas, rna, temper, ctol/8.3144d0,6.0232d23,3.d02,1.d-8/
c
      dcdep  = zero
      dsdc   = zero
      dcdskk = zero
      ak(:,:) = zero
c
      ak(1,1) = two_third
      ak(2,2) = two_third
      ak(3,3) = two_third
      ak(1,2) = - one_third
      ak(1,3) = - one_third
      ak(2,1) = - one_third
      ak(2,3) = - one_third
      ak(3,1) = - one_third
      ak(3,2) = - one_third
      ak(4,4) = half
      ak(5,5) = half
      ak(6,6) = half
c
c              dmat[] zeroed by warp3d
c
      do i = 1, span
c
       alpha = mm_props(i,1)
       rho0  = mm_props(i,2)
       gamma = mm_props(i,3)
       alat  = mm_props(i,4)
       wb    = mm_props(i,5)
       beta  = mm_props(i,6)
       vm    = mm_props(i,7)
       vh    = mm_props(i,8)
       c0    = mm_props(i,9)
       rksi  = mm_props(i,10)
c
       alamda= vh/vm
       rnl   = rna/vm
c
       fact = one
c
       sigma_0   = sigyld_vec(i)
       eps_0     = sigma_0 / e_vec(i)
       n_power   = n_power_vec(i)
c
       bk = e_vec(i) /
     &       ( three * ( one - two * nu_vec(i)))
       shear_mod = e_vec(i) /
     &       ( two * ( one + nu_vec(i)))
c
c               c1 is actually lamda
c
       c1 = nu_vec(i) * two * shear_mod
     &       /( one - two * nu_vec(i))
c
       c2 = c1 + two * shear_mod
       c3 = shear_mod
c
       dmat(i,1,1) = c2
       dmat(i,2,2) = c2
       dmat(i,3,3) = c2
       dmat(i,1,2) = c1
       dmat(i,1,3) = c1
       dmat(i,2,1) = c1
       dmat(i,2,3) = c1
       dmat(i,3,1) = c1
       dmat(i,3,2) = c1
       dmat(i,4,4) = c3
       dmat(i,5,5) = c3
       dmat(i,6,6) = c3
c
       hcon = history_np1(i,2)
       delc = history_np1(i,2) - history_n(i,2)
       etran= history_np1(i,3)
       rlam = three * alamda / ( three+ etran )
       drlamdc = -rlam * rlam / three
c
       eps = history_np1(i,1)
       skk = ( stress_np1(i,1)
     +          + stress_np1(i,2)
     +          + stress_np1(i,3) ) / three
       call mm07_hydrogen( iout, alpha, rho0, gamma, alat,
     &       wb, beta, vm, vh, c0, rksi, alamda, rnl, rna, rgas,
     &       temper, eps, skk, cnow, thl, tht, dcdskk, dcdep,
     &       dntdeps, cl, ct)
c
       if( history_np1(i,4) .lt. zero ) then
         fac1 = bk*(drlamdc*delc+rlam)*three *bk*dcdskk
     +       /(one + three*bk*(delc*drlamdc+rlam)*dcdskk)
         fac2 = fac1*fact
         dmat(i,1,1) = dmat(i,1,1)-fac2
         dmat(i,1,2) = dmat(i,1,2)-fac2
         dmat(i,1,3) = dmat(i,1,3)-fac2
         dmat(i,2,1) = dmat(i,2,1)-fac2
         dmat(i,2,2) = dmat(i,2,2)-fac2
         dmat(i,2,3) = dmat(i,2,3)-fac2
         dmat(i,3,1) = dmat(i,3,1)-fac2
         dmat(i,3,2) = dmat(i,3,2)-fac2
         dmat(i,3,3) = dmat(i,3,3)-fac2
       end if
c
       if( history_np1(i,4) .gt. zero ) then
          eps_np1 = history_np1(i,1)
          hcon_np1 = history_np1(i,2)
          deleps  = history_np1(i,5)
c
          pe = ( trial_elas_stress(i,1)
     &         + trial_elas_stress(i,2)
     &         + trial_elas_stress(i,3) ) / three
c
          str_epd_1(1) = trial_elas_stress(i,1) - pe
          str_epd_1(2) = trial_elas_stress(i,2) - pe
          str_epd_1(3) = trial_elas_stress(i,3) - pe
          str_epd_1(4) = trial_elas_stress(i,4)
          str_epd_1(5) = trial_elas_stress(i,5)
          str_epd_1(6) = trial_elas_stress(i,6)
c
          fac1 = str_epd_1(1) * str_epd_1(1)
     &         + str_epd_1(2) * str_epd_1(2)
     &         + str_epd_1(3) * str_epd_1(3)
     &        + (str_epd_1(4) * str_epd_1(4)
     &         + str_epd_1(5) * str_epd_1(5)
     &         + str_epd_1(6) * str_epd_1(6))* two
          qe = sqrt ( onep5 * fac1 )
c
          an(1) = onep5 * str_epd_1(1) / qe
          an(2) = onep5 * str_epd_1(2) / qe
          an(3) = onep5 * str_epd_1(3) / qe
          an(4) = onep5 * str_epd_1(4) / qe
          an(5) = onep5 * str_epd_1(5) / qe
          an(6) = onep5 * str_epd_1(6) / qe
c
          call mm07_yield( hcon_np1, rksi, eps_np1, sigma_0, eps_0,
     &                     n_power, sigmaY, dsdep, dsdc )
c
          fac1 = dsdep+three*shear_mod
          fac2 =  one + three * bk*(delc*drlamdc+rlam)*dcdskk
          fac3 = fac1*fac2+dsdc*dcdep
c
          delta(1) = one
          delta(2) = one
          delta(3) = one
          delta(4) = zero
          delta(5) = zero
          delta(6) = zero
c
          dcde(:)  =  zero
          depde(:) =  zero
c
          dcde(1) = (three *bk*fac1*dcdskk*delta(1)
     +          + two * shear_mod*dcdep*an(1))/fac3
c
          dcde(2) = (three *bk*fac1*dcdskk*delta(2)
     +          +two *shear_mod*dcdep*an(2))/fac3
c
          dcde(3) = (three *bk*fac1*dcdskk*delta(3)
     +          +two *shear_mod*dcdep*an(3))/fac3
c
          dcde(4) = (three *bk*fac1*dcdskk*delta(4)
     +          +two *shear_mod*dcdep*an(4))/fac3
c
          dcde(5) = (three *bk*fac1*dcdskk*delta(5)
     +          +two *shear_mod*dcdep*an(5))/fac3
c
          dcde(6) = (three *bk*fac1*dcdskk*delta(6)
     +          +two *shear_mod*dcdep*an(6))/fac3
c
          depde(1) = (two *shear_mod*fac2*an(1)
     +             -three *bk*dsdc*dcdskk*delta(1))/fac3
c
          depde(2) = (two *shear_mod*fac2*an(2)
     +             -three *bk*dsdc*dcdskk*delta(2))/fac3
c
          depde(3) = (two *shear_mod*fac2*an(3)
     +             -three *bk*dsdc*dcdskk*delta(3))/fac3
c
          depde(4) = (two *shear_mod*fac2*an(4)
     +             -three *bk*dsdc*dcdskk*delta(4))/fac3
c
          depde(5) = (two *shear_mod*fac2*an(5)
     +             -three *bk*dsdc*dcdskk*delta(5))/fac3
c
          depde(6) = (two *shear_mod*fac2*an(6)
     +             -three *bk*dsdc*dcdskk*delta(6))/fac3
c
          fac4 = two * shear_mod * shear_mod * deleps / qe
          do jj = 1, 6
            do kk = 1, 6
                 dmat(i,jj,kk) = dmat(i,jj,kk)
     +           +fac4*(two *an(jj)*an(kk) - three *ak(jj,kk))*fact
     +           -shear_mod*(an(jj)*depde(kk)+an(kk)*depde(jj))*fact
     +           -bk*(drlamdc*delc+rlam)*fact
     +           * half *(delta(jj)*dcde(kk)+delta(kk)*dcde(jj))
            end do
          end do
c
         end if  !  history_np1(i,4) .gt. zero
c
      end do
c
      return
      end
c *******************************************************************
c *                                                                 *
c *        material model # 7 -- mises + hydrogen model             *
c *                                                                 *
c *           set 3 material model dependent output values          *
c *                                                                 *
c *******************************************************************
c
c
      subroutine oumm07( gpn, mxvl, span, iout, elestr,
     &                   stress, history )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &  gpn, mxvl, span, iout,i
c
c
      double precision
     & stress(mxvl,*), elestr(mxvl,*), history(mxvl,*),
     & p, two, three, onep5, fac1
c
      data onep5, two, three/1.5d00, 2.d00, 3.d00/
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
c
c
c
       do i = 1, span
c          elestr(i,1) = stress(i,1)
c          elestr(i,2) = stress(i,2)
c          elestr(i,3) = stress(i,3)
c          elestr(i,4) = stress(i,4)
c          elestr(i,5) = stress(i,5)
c          elestr(i,6) = stress(i,6)
c          elestr(i,7) = stress(i,7)
c          p = ( stress(i,1) + stress(i,2) + stress(i,3))
c     &      / three
c          fac1 = ( stress(i,1) - p ) * ( stress(i,1) - p )
c     &         + ( stress(i,2) - p ) * ( stress(i,2) - p )
c     &         + ( stress(i,3) - p ) * ( stress(i,3) - p )
c     &         + ( stress(i,4) * stress(i,4)
c     &         +   stress(i,5) * stress(i,5)
c     &         +   stress(i,6) * stress(i,6) )* two
c
c          elestr(i,8) = sqrt( onep5 * fac1 )
c
          elestr(i,9) = history(i,1)
          elestr(i,10)= history(i,6)
          elestr(i,11)= history(i,7)
       enddo
c
       return
       end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm07_set_sizes                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 12/14/14  rhd               *
c     *                                                              *
c     *    called by warp3d for each material model to obtain        *
c     *    various sizes of data for the model                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm07_set_sizes( info_vector )
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
      info_vector(1) = 7
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
c     *             subroutine mm07_states_values                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 1/3/2015 (rhd))                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm07_states_values( itype, elem_states_output,
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
c     *                 subroutine mm07_states_labels                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 1/11/2015 (rhd)                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm07_states_labels( size_state,
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

