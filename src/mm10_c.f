c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine mm10_solveB                  *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/23/2016 rhd              *
c     *                                                              *
c     *              solve for stress/hardening using interface      *
c     *              to "nwnleq" solver package                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_solveB( props, np1, np0, vec1, vec2, arr1, arr2,
     &                        ivec1, ivec2, stress, tt, fail, faili,
     &                        failr, gaspt, dtinc, rjac)
c
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c
c              parameters
c
      type(crystal_props) :: props
      type(crystal_state) :: np1, np0
      double precision, dimension(max_uhard) :: vec1, vec2
      double precision, dimension(max_uhard,max_uhard) :: arr1, arr2
      complex(kind=real64), dimension(max_uhard) :: ivec1, ivec2
      double precision, dimension(6) :: stress
      double precision, dimension(props%num_hard) :: tt
      logical :: fail
      integer, dimension(10) :: faili
      double precision, dimension(10) :: failr
      integer :: gaspt
      double precision :: dtinc
c
c               locals
c
      include 'include_mm10'
      type(mm10_working_data) :: solve_work
c
      double precision :: nR, inR, atol, rtol, uB, alpha, ls1, ls2,
     &      nlsx, nRs, dt, cos_ang, xetol, xtol1, zerotol,
     &      dxerr, nR1, inR1, atol1, rtol1, dp1
      double precision, dimension(6) :: R1, x1, dx1, xnew1, d1, d2
      double precision, dimension(6,6) :: J11
      integer :: miter, info, ls, mmin, gpp, ttind
      logical :: debug, gpall, locdebug, jaccheck, numerjac, predict_ok,
     &           update_ok 
c
c               automatics
c
      double precision, dimension(6+props%num_hard) :: R, x, dx, 
     &                                                 xnew, g
      double precision, 
     & dimension(6+props%num_hard,6+props%num_hard) :: J
      double precision, dimension(props%num_hard) :: x2
      integer, dimension(6+props%num_hard) :: ipiv
c      
c               local for solver package
c
      integer ::  n,jacflg(5),maxit,njcnt,nfcnt,iter,termcd,method,
     &            global,xscalm,ldr,lrwork,qrwsiz,outopt(5)
      double precision :: xtol,ftol,btol,cndtol,stepmx,delta,sigma,
     &                    t1,t2,s_trace,dtrace
      double precision :: rjac1(6,2*6)
c      
c               automatics for solver package
c
      double precision ::  xp(6+props%num_hard),fp(6+props%num_hard),
     &                     gp(6+props%num_hard),
     &                     rjac(6+props%num_hard,2*(6+props%num_hard)),
     &                     rwork(9*(6+props%num_hard)),
     &                     rcdwrk(3*(6+props%num_hard)),
     &                     qrwork(props%num_hard*(6+props%num_hard)),
     &                     scalex(6+props%num_hard)
      integer  :: icdwrk(6+props%num_hard)
c
c               convergence parameters: Newton w/ geometric line search
c
      double precision, parameter :: c = 1.0d-4, red = 0.5d0
      integer, parameter :: mls = 10, min = 1

c!DIR$ ASSUME_ALIGNED vec1:64, vec2:64, arr1:64, arr2:64, ivec1:64
c!DIR$ ASSUME_ALIGNED ivec2:64, stress:64, tt:64



      atol = props%atol
      atol1 = props%atol1
      rtol = props%rtol
      rtol1 = props%rtol1
      xetol = props%xtol
      xtol1 = props%xtol1
      miter = props%miter
      zerotol = 1.0d-12
      ttind = 1 ! index of tt to print while debugging
c234567890123456
c              tang_calc: variable to denote Jacobian calculation
c              set this variable INSIDE mod_crystals.f
c                   0 for user supplied, no checking
c                   1 for real finite difference (mm10_solveB only)
c                   2 for complex finite difference
c                   3 for checking with real finite difference
c                   4 for checking with complex finite difference
c
c              flag for checking user supplied Jacobian numerically
c              flag for computing Jacobian numerically
c              debugging flags
c
      if( (props%tang_calc .eq. 3) .or.
     &    (props%tang_calc .eq. 4) ) then
        jaccheck = .true. 
      else
        jaccheck = .false. 
      end if
c     
      if( (props%tang_calc .eq. 1) ) then
        numerjac = .true.
      else
        numerjac = .false. ! complex version handled inside of J function
      end if
c
      debug = props%debug
      gpall = props%gpall ! true to print iteration norms all GPs
      gpp  = props%gpp ! set one particular G.P. to print
c
c      locdebug = (debug .and.(gpall.or.(gaspt.eq.gpp))
c     & .and.((props%st_it(3).eq.-2).or.
c     &       (props%st_it(3).eq.np1%elem))
c     & .and.((props%st_it(1).eq.-2).or.
c     &       (props%st_it(1).eq.np1%step))
c     & .and.((props%st_it(2).eq.-2).or.
c     &       (props%st_it(2).eq.np1%iter)))
      locdebug = .false.
      solve_work%locdebug = locdebug
c
      if( locdebug ) write(props%out,*) "Entering solution routine"
c
c              initial values for solution from last time step
c
      x(1:6) = stress(1:6)
      x(7:props%num_hard+6) = tt(1:props%num_hard)
c
c              prediction of yield stress to initialize
c              the integration algorithm; helps
c              for Orowan flow rule type models
c
c              extrapolation --
c              use cosine of the "angle" between the new and old
c              displacement increment to indicate continuity of load
c              direction
c
      if( props%h_type .gt. 3 ) then
          call mm10_solveB_stress_predict
          if( .not. predict_ok ) return
      else 
        inR1 = one
      end if 
c
c              material update algorithm
c
      if( locdebug ) then
         write(props%out,*)" Material update module, G.P.=", gaspt
         write(props%out,*)" Guess syy=",x(2), " Guess tt=", x(6+ttind)
      end if
c
      call mm10_solveB_stress_update
c
      return

      contains
c     ========

      subroutine  mm10_solveB_stress_predict
      implicit none
c
c              prediction of yield stress to initialize
c              the integration algorithm; helps
c              for Orowan flow rule type models
c
c              extrapolation --
c              use cosine of the "angle" between the new and old
c              displacement increment to indicate continuity of load
c              direction
c
      iter = 0
      x1(1:6) = x(1:6)
      dt = dtinc !np1%tinc
      d1(1:6) = np1%D(1:6)
      dtrace = (d1(1) + d1(2) + d1(3))/three
      d1(1) = d1(1) - dtrace
      d1(2) = d1(2) - dtrace
      d1(3) = d1(3) - dtrace
      t1 = d1(1)**2 + d1(2)**2 + d1(3)**2
      t2 = d1(4)**2 + d1(5)**2 + d1(6)**2
      if( t1+t2 .eq. zero ) then
        d1(1:6) = zero
      else
        d1(1:6) = d1(1:6)/sqrt( t1+t2 )
      end if
c
      d2(1:6) = np0%D(1:6)
      dtrace = (d2(1) + d2(2) + d2(3))/three
      d2(1) = d2(1) - dtrace
      d2(2) = d2(2) - dtrace
      d2(3) = d2(3) - dtrace
      t1 = d2(1)**2 + d2(2)**2 + d2(3)**2
      t2 = d2(4)**2 + d2(5)**2 + d2(6)**2
      if( t1+t2 .eq. zero ) then
        d2(1:6) = zero
      else
        d2(1:6) = d2(1:6)/sqrt( t1+t2 )
      end if
      t1 = d1(1)*d2(1) + d1(2)*d2(2) + d1(3)*d2(3)
      t2 = d1(4)*d2(4) + d1(5)*d2(5) + d1(6)*d2(6)
      cos_ang = dmax1( t1+t2, zero )
      x2(1:props%num_hard) = x(7:6+props%num_hard) + 
     &     cos_ang*np0%tt_rate(1:props%num_hard)*dt
c
      if( locdebug ) then
         write(props%out,*)" Stress prediction module, G.P.=", gaspt
         write(props%out,*)" Extrapol tt6=",
     &      x2(ttind), " Previous tt6=", x(6+ttind)
      end if 
c
c              copy stuff into solve_work. set solver arguments
c
      solve_work%props = props
      solve_work%np1 = np1
      solve_work%np0 = np0
      solve_work%vec1(1:max_uhard) = vec1(1:max_uhard)
      solve_work%vec2(1:max_uhard) = vec2(1:max_uhard)
      solve_work%arr1(1:max_uhard,1:max_uhard) = 
     &           arr1(1:max_uhard,1:max_uhard)
      solve_work%arr2(1:max_uhard,1:max_uhard) = 
     &           arr2(1:max_uhard,1:max_uhard)
      solve_work%ivec1(1:max_uhard) = ivec1(1:max_uhard)
      solve_work%ivec2(1:max_uhard) = ivec2(1:max_uhard)
      solve_work%gaspt = gaspt
      solve_work%solvfnc = 1
      solve_work%x2(1:props%num_hard) = x2(1:props%num_hard)
c
      n = 6
      t1      = stress(1)**2 + stress(2)**2 + stress(3)**2
      t2      = stress(4)**2 + stress(5)**2 + stress(6)**2
      s_trace = sqrt( three/two * ( t1 + two*t2 ) )
      if( s_trace .eq. zero ) then
         t1      = np1%D(1)**2 + np1%D(2)**2 + np1%D(3)**2
         t2      = np1%D(4)**2 + np1%D(5)**2 + np1%D(6)**2
         s_trace = dsqrt( onept5 * ( t1 + two*t2 ) )
         stepmx = props%mu_0*s_trace
         scalex(1:6) = one/stepmx
      else
         scalex(1:6) = one/s_trace
      end if
      maxit = miter
      if( numerjac ) then
        jacflg(1) = 0 ! finite difference 1 ! analytical 
      else
        jacflg(1) = 1 ! analytical 0 ! finite difference 
      endif
      jacflg(2) = -1
      jacflg(3) = -1
      jacflg(4) = 1
      xtol = xtol1 !x tolerance
      ftol = rtol1 !f tolerance
      btol = xtol1 !tolerance for backtracking
      cndtol = 1.0d-14 !tolerance of test for ill conditioning
      if( props%method .le. 6 ) then
        method = 0 ! Newton
        global = props%method ! linesearch/trust-region
      else
        method = 1 ! Broyden
        global = props%method-7 !linesearch/trust-region
      end if
      xscalm = 0 !for manual scaling
      stepmx = zero !maximum stepsize
      delta = -two !  ==> use min(Newton length, stepmx)
      sigma = red !reduction factor geometric linesearch
      ldr = 6
      lrwork = 9*(6+props%num_hard) 
c
c              setup size for QR decomposition
c
      call mm10_liqsiz( n, qrwsiz )
      if( locdebug ) then
        outopt(1) = 1
      else
        outopt(1) = 0
      end if
      if( jaccheck ) then
        outopt(2) = 1
      else
        outopt(2) = 0
      end if
      outopt(3) = 1 ! yes, output Jacobian
c
c              perform the actual solve for stress prediction
c
      call mm10_nwnleq(x1,n,scalex,maxit,
     *                  jacflg,xtol,ftol,btol,cndtol,method,global,
     *                  xscalm,stepmx,delta,sigma,rjac1,ldr,
     *                  rwork,lrwork,
     *                  rcdwrk,icdwrk,qrwork,qrwsiz,solve_work,outopt,
     *                  xp,fp,gp,njcnt,nfcnt,iter,termcd)
c
c              what happened; did we converge or not. set
c              causes of failure
c
      dp1 = fp(1)**2 + fp(2)**2 + fp(3)**2 + fp(4)**2 + fp(5)**2 +
     &      fp(6)**2 
      nR1 = dsqrt( dp1 )
      if( (iter > miter) .or. any(isnan(xp(1:6))) .or.
     &     (termcd ==5) .or. (termcd == 4) ) then  ! failed :(
          fail = .true.
          faili(1) = iter
          faili(2) = miter
          if( nR1 .gt. atol1 ) then
            faili(3) = 1
          elseif( any( isnan(x1) ) ) then
            faili(3) = 3
          end if
          faili(4) = 1
          failr(1) = nR1
          failr(2) = atol1
          failr(4) = rtol1
          if( locdebug ) then
            write(props%out,*) 'termcd', termcd, 'miter', miter, 
     &                          'iter', iter
          end if
          predict_ok = .false.
          return  ! Return out for failure to converge
      end if
c
c              copy stuff out of solver output.
c              optionally output statistics from N-R algorithm
c
      x1(1:6) = xp(1:6)
      if( locdebug ) then 
        write(props%out,*) 'termcd', termcd, 'miter', miter, 'iter',
     &                    iter
        write(props%out,*)" Stress pred conv iter=",  iter
        write(props%out,*)" AbsNorm=", nR1, " AbsTol=", atol1
        write(props%out,*)" Guess syy=",x(2), " actual syy=", x1(2)
      end if
c
c             copy predicted stress and hardening back
c             to primary variable x
c
      x(1:6) = x1(1:6)
      x(7:6+props%num_hard) = x2(1:props%num_hard)
c
      predict_ok = .true.
c      
      return
c
      end subroutine mm10_solveB_stress_predict

      subroutine mm10_solveB_stress_update
      implicit none
c      
      solve_work%props = props
      solve_work%np1 = np1
      solve_work%np0 = np0
      solve_work%vec1(1:max_uhard) = vec1(1:max_uhard)
      solve_work%vec2(1:max_uhard) = vec2(1:max_uhard)
      solve_work%arr1(1:max_uhard,1:max_uhard) = 
     &           arr1(1:max_uhard,1:max_uhard)
      solve_work%arr2(1:max_uhard,1:max_uhard) = 
     &           arr2(1:max_uhard,1:max_uhard)
      solve_work%ivec1(1:max_uhard) = ivec1(1:max_uhard)
      solve_work%ivec2(1:max_uhard) = ivec2(1:max_uhard)
      solve_work%gaspt = gaspt
      solve_work%solvfnc = 2
c
c         Set up solver arguments
c
      n = 6+props%num_hard
      t1      = x(1)**2 + x(2)**2 + x(3)**2
      t2      = x(4)**2 + x(5)**2 + x(6)**2
      s_trace = dsqrt( onept5 * ( t1 + two*t2 ) )
      if( s_trace .eq. zero ) then
        t1      = np1%D(1)**2 + np1%D(2)**2 + np1%D(3)**2
        t2      = np1%D(4)**2 + np1%D(5)**2 + np1%D(6)**2
        s_trace = dsqrt( onept5 * ( t1 + two*t2 ) )
        stepmx = props%mu_0*s_trace
        scalex(1:6) = one/stepmx
      else
        scalex(1:6) = one/s_trace
      end if
c
      scalex(7:6+props%num_hard) = one / x(7:6+props%num_hard)
      maxit = miter
      if( numerjac ) then
        jacflg(1) = 0 ! finite difference 1 ! analytical 
      else
        jacflg(1) = 1 ! analytical 0 ! finite difference 
      endif
      jacflg(2) = -1
      jacflg(3) = -1
      jacflg(4) = 1
      xtol = xetol !x tolerance
      ftol = rtol !f tolerance
      btol = xetol !tolerance for backtracking
      cndtol = 1.0d-14 !tolerance of test for ill conditioning
      if( props%method <= 6 ) then
        method = 0 ! Newton
        global = props%method ! linesearch/trust-region
      else
        method = 1 ! Broyden
        global = props%method-7 ! linesearch/trust-region
      end if
      xscalm = 0 !manual scaling
      stepmx = zero !maximum stepsize
      delta = -two !  ==> use min(Newton length, stepmx)
      sigma = red !reduction factor geometric linesearch
      ldr = 6+props%num_hard
      lrwork = 9*(6+props%num_hard)
c
c              setup size for QR decomposition
c
      call mm10_liqsiz(n,qrwsiz)
      if( locdebug ) then
        outopt(1) = 1
      else
        outopt(1) = 0
      end if
      if( jaccheck ) then
        outopt(2) = 1
      else
        outopt(2) = 0
      end if
      outopt(3) = 1 ! yes, output Jacobian
c
c              perform the actual solve for stress prediction
c
      call mm10_nwnleq(x,n,scalex,maxit,
     *                  jacflg,xtol,ftol,btol,cndtol,method,global,
     *                  xscalm,stepmx,delta,sigma,rjac,ldr,
     *                  rwork,lrwork,
     *                  rcdwrk,icdwrk,qrwork,qrwsiz,solve_work,outopt,
     *                  xp,fp,gp,njcnt,nfcnt,iter,termcd)
c
c              what happened; did we converge or not. set
c              causes of failure
c
      dp1 = fp(1)**2 + fp(2)**2 + fp(3)**2 + fp(4)**2 + fp(5)**2 +
     &      fp(6)**2 
      nR = dsqrt( dp1 )
      if( (iter > miter) .or. any(isnan(xp)) .or.
     &     (termcd == 5) .or. (termcd ==4) ) then
          fail = .true.
          faili(1) = iter
          faili(2) = miter
          if( nR > atol ) then
            faili(3) = 1
          elseif( any( isnan(x) ) ) then
            faili(3) = 3
          end if
          faili(4) = 2
          failr(1) = nR
          failr(2) = atol
          failr(4) = rtol
          if( locdebug ) then 
            write(props%out,*) 'termcd', termcd, 'miter', miter, 
     &               'iter', iter
          end if ! Return out for failure to converge
          update_ok = .false.
          return
      end if 
c
c              copy stuff out of solver output.
c              optionally output statistics from N-R algorithm
c
      if( locdebug ) then
        write(props%out,*) 'termcd', termcd, 'miter', miter, 
     &            'iter', iter
      end if
c      
      x(1:6+props%num_hard) = xp(1:6+props%num_hard)
      np1%tt_rate(1:props%num_hard) = 
     &     solve_work%np1%tt_rate(1:props%num_hard)
c
      if( locdebug ) then 
        write(props%out,*)" Material upd conv, iter=", iter
        write(props%out,*)" AbsNorm=", nR, " AbsTol=", atol
        write(props%out,*)" Guess syy=", stress(2), 
     &             " actual syy=", x(2)
        write(props%out,*)" Guess tt6=", x2(ttind), 
     &            " actual tt6=", x(6+ttind)
      end if
c
c              good stress update. all done
c
      stress(1:6) = x(1:6)
      tt(1:props%num_hard) = x(7:6+props%num_hard)
      update_ok = .true.
c
      return

      end subroutine mm10_solveB_stress_update
c      
      end subroutine mm10_solveB

c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine mm10_nwnleq                  *
c     *                                                              *
c     *              entry point for nwnleq solver package           *
c     *                                                              *
c     ****************************************************************
c
c ----------------------------------------
c
      subroutine mm10_nwnleq(x0,n,scalex,maxit,
     *                  jacflg,xtol,ftol,btol,cndtol,method,
     *                  global,xscalm,stepmx,delta,sigma,rjac,
     *                  ldr,rwork,lrwork,
     *                  rcdwrk,icdwrk,qrwork,qrwsiz,solve_work,
     *                  outopt,xp,fp,gp,njcnt,nfcnt,iter,termcd)
      use iso_Fortran_env
      use mm10_defs
      integer n,jacflg(5),maxit,njcnt,nfcnt,iter,termcd,method
      integer global,xscalm,ldr,lrwork,qrwsiz
      integer outopt(*)
      double precision  xtol,ftol,btol,cndtol,stepmx,delta,sigma
      double precision  xp(n),fp(n),gp(*),x0(*)
      double precision  rjac(ldr,*),rwork(*),rcdwrk(*),qrwork(*)
      double precision  scalex(*)
      integer           icdwrk(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work


c-------------------------------------------------------------------------
c
c     Solves systems of nonlinear equations using the Newton / Broyden
c     method with a global strategy either linesearch or double dogleg
c
c     In       x0      Real(*)         starting vector for x
c     In       n       Integer         dimension of problem
c     Inout    scalex  Real(*)         scaling factors x()
c     Inout    maxit   Integer         maximum number iterations
c     Inout    jacflg  Integer(*)      jacobian flag array
c                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
c                                                  3: user supplied banded
c                                      jacflg[2]: number of sub diagonals or -1 if not banded
c                                      jacflg[3]: number of super diagonals or -1 if not banded
c                                      jacflg[4]: 1 if adjusting step allowed when
c                                                   singular or illconditioned
c     Inout    xtol    Real            x tolerance
c     Inout    ftol    Real            f tolerance
c     Inout    btol    Real            x tolerance for backtracking
c     Inout    cndtol  Real            tolerance of test for ill conditioning
c     Inout    method  Integer         method to use
c                                        0 Newton
c                                        1 Broyden
c     In       global  Integer         global strategy to use
c                                        1 cubic linesearch
c                                        2 quadratic linesearch
c                                        3 geometric linesearch
c                                        4 double dogleg
c                                        5 powell dogleg
c                                        6 hookstep (More-Hebden Levenberg-Marquardt)
c     In       xscalm  Integer         scaling method
c                                        0 scale fixed and supplied by user
c                                        1 for scale from jac. columns a la Minpack
c     Inout    stepmx  Real            maximum stepsize
c     Inout    delta   Real            trust region radius
c                                        > 0.0 or special value for initial value
c                                        -1.0  ==> use min(Cauchy length, stepmx)
c                                        -2.0  ==> use min(Newton length, stepmx)
c     Inout    sigma   Real            reduction factor geometric linesearch
c     Inout    rjac    Real(ldr,*)     workspace jacobian
c                                         2*n*n for Broyden and n*n for Newton
c     In       ldr     Integer         leading dimension rjac
c     Out      rwork   Real(*)         real workspace (9n)
c     In       lrwork  Integer         size real workspace
c     In       rcdwrk  Real(*)         workspace for Dtrcon (3n)
c     In       icdwrk  Integer(*)      workspace for Dtrcon (n)
c     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
c     In       qrwsiz  Integer         size of qrwork
c     In       fjac    Name            optional name of routine to calculate
c                                      user supplied jacobian
c     In       fvec    Name            name of routine to calculate f(x)
c     In       outopt  Integer(*)      output options
c                                       outopt(1)
c                                         0 no output
c                                         1 output an iteration report
c                                       outopt(2)
c                                         0 do not check any user supplied jacobian
c                                         1 check user supplied jacobian if supplied
c     Out      xp      Real(*)         final values for x()
c     Out      fp      Real(*)         final values for f(x)
c     Out      gp      Real(*)         gradient of f() at xp()
c     Out      njcnt   Integer         number of jacobian evaluations
c     Out      nfcnt   Integer         number of function evaluations
c     Out      iter    Integer         number of (outer) iterations
c     Out      termcd  Integer         termination code
c                                       > 0 process terminated
c                                             1  function criterion near zero
c                                             2  no better point found
c                                             3  x-values within tolerance
c                                             4  iteration limit exceeded
c                                             5  singular/ill-conditioned jacobian
c                                             6  totally singular jacobian
c                                                (when allowSingular=TRUE)
c
c                                       < 0 invalid input parameters
c                                            -1  n not positive
c                                            -2  insufficient workspace rwork
c                                            -3  cannot check user supplied jacobian (not supplied)
c
c    The subroutine fvec must be declared as
c
c!        subroutine fvec(x,f,n,flag)
c         double precision x(*), f(*)
c         integer  n, flag
c
c         x() are the x values for which to calculate the function values f(*)
c         The dimension of these vectors is n
c         The flag argument is set to
c            0  for calculation of function values
c           >0  indicating that jacobian column <flag> is being computed
c               so that fvec can abort.
c
c    The subroutine fjac must be declared as
c
c!        subroutine mkjac(rjac,ldr,x,n)
c         integer ldr
c         double precision rjac(ldr,*), x(*)
c         integer  n
c
c         The routine calculates the jacobian in point x(*) of the
c         function. If any illegal values are encountered during
c         calculation of the jacobian it is the responsibility of
c         the routine to quit.

c-------------------------------------------------------------------------

      double precision epsm

c     check input parameters

      call mm10_nwpchk(n,lrwork,xtol,ftol,btol,cndtol,maxit,
     *            jacflg,method,global,stepmx,delta,sigma,
     *            epsm,outopt,scalex,xscalm,termcd)
      if(termcd .lt. 0) then
         return
      endif

c     first argument of nwsolv/brsolv is leading dimension 
c      of rjac in those routines
c     should be at least n

      if( method .eq. 0 ) then

         call mm10_nwsolv(ldr,x0,n,scalex,maxit,jacflg,
     *               xtol,ftol,btol,cndtol,global,xscalm,
     *               stepmx,delta,sigma,
     *               rjac,
     *               rwork(1    ),rwork(1+  n),
     *               rwork(1+2*n),rwork(1+3*n),
     *               rwork(1+4*n),rwork(1+5*n),
     *               rwork(1+6*n),rwork(1+7*n),
     *               rwork(1+8*n),rcdwrk,icdwrk,qrwork,qrwsiz,epsm,
     *               solve_work,outopt,xp,fp,gp,njcnt,nfcnt,iter,termcd)

      elseif( method .eq. 1 ) then

         call mm10_brsolv(ldr,x0,n,scalex,maxit,jacflg,
     *               xtol,ftol,btol,cndtol,global,xscalm,
     *               stepmx,delta,sigma,
     *               rjac,rjac(1,n+1),
     *               rwork(1    ),rwork(1+  n),
     *               rwork(1+2*n),rwork(1+3*n),
     *               rwork(1+4*n),rwork(1+5*n),
     *               rwork(1+6*n),rwork(1+7*n),
     *               rwork(1+8*n),rcdwrk,icdwrk,
     *               qrwork,qrwsiz,epsm,
     *               solve_work,outopt,xp,fp,gp,njcnt,nfcnt,iter,termcd)

      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwpchk(n,lrwk,
     *                  xtol,ftol,btol,cndtol,maxit,jacflg,method,
     *                  global,stepmx,delta,sigma,epsm,outopt,
     *                  scalex,xscalm,termcd)

      integer n,lrwk,jacflg(5)
      integer method,global,maxit,xscalm,termcd
      integer outopt(*)
      double precision  xtol,ftol,btol,cndtol,stepmx,delta,sigma
      double precision  scalex(*),epsm

c-------------------------------------------------------------------------
c
c     Check input arguments for consistency and modify if needed/harmless
c
c     Arguments
c
c     In       n       Integer         dimension of problem
c     In       lrwk    Integer         size real workspace
c     Inout    xtol    Real            x tolerance
c     Inout    ftol    Real            f tolerance
c     Inout    btol    Real            x tolerance for backtracking
c     Inout    cndtol  Real            tolerance of test for ill conditioning
c     Inout    maxit   Integer         maximum number iterations
c     Inout    jacflg  Integer(*)      jacobian flag
c     Inout    method  Integer         method to use (Newton/Broyden)
c     Inout    global  Integer         global strategy to use
c     Inout    stepmx  Real            maximum stepsize
c     Inout    delta     Real            trust region radius
c     Inout    sigma   Real            reduction factor geometric linesearch
c     Out      epsm                    machine precision
c     Inout    scalex  Real(*)         scaling factors x()
c     Inout    xscalm  integer         0 for fixed scaling, 1 for automatic scaling
c     Out      termcd  Integer         termination code (<0 on errors)
c
c-------------------------------------------------------------------------

      integer i,len
      double precision mm10_epsmch, mm10_dblhuge, Rhuge

      double precision Rzero, Rone, Rtwo, Rthree
      parameter(Rzero=0.0d0, Rone=1.0d0, Rtwo=2.0d0, Rthree=3.0d0)

      double precision Rhalf
      parameter(Rhalf = 0.5d0)

c     check that parameters only take on acceptable values
c     if not, set them to default values

c     initialize termcd to all ok

      termcd = 0

c     compute machine precision

      epsm = mm10_epsmch()

c     get largest double precision number
      Rhuge = mm10_dblhuge()

c     check dimensions of the problem

      if(n .le. 0) then
         termcd = -1
         return
      endif

c     check dimensions of workspace arrays

      len = 9*n
c      +2*n*n
      if(lrwk .lt. len) then
         termcd = -2
         return
      endif

c     check jacflg, method, and global

      if(jacflg(1) .gt. 3 .or. jacflg(1) .lt. 0) jacflg(1) = 0

      if(method .lt. 0 .or. method .gt. 1) method = 0

      if(global .lt. 0 .or. global .gt. 6) global = 4

c     set outopt to correct values

      if(outopt(1) .ne. 0 ) then
         outopt(1) = 1
      endif

      if(outopt(2) .ne. 0 ) then
         outopt(2) = 1
      endif

c     check scaling scale matrices

      if(xscalm .eq. 0) then
         do i = 1,n
            if(scalex(i) .lt. Rzero) scalex(i) = -scalex(i)
            if(scalex(i) .eq. Rzero) scalex(i) = Rone
         enddo
      else
         xscalm = 1
         do i = 1,n
            scalex(i) = Rone
         enddo
      endif
c     check step and function tolerances

      if(xtol .lt. Rzero) then
         xtol = epsm**(Rtwo/Rthree)
      endif

      if(ftol .lt. Rzero) then
         ftol = epsm**(Rtwo/Rthree)
      endif

      if( btol .lt. xtol ) btol = xtol

      cndtol = max(cndtol, epsm)

c     check reduction in geometric linesearch

      if( sigma .le. Rzero .or. sigma .ge. Rone ) then
         sigma = Rhalf
      endif

c     check iteration limit

      if(maxit .le. 0) then
         maxit = 150
      endif

c     set stepmx

      if(stepmx .le. Rzero) stepmx = Rhuge

c     check delta
      if(delta .le. Rzero) then
         if( delta .ne. -Rtwo ) then
            delta = -Rone
         endif
      elseif(delta .gt. stepmx) then
         delta = stepmx
      endif

      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine mm10_fjac                    *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/23/2016 rhd              *
c     *                                                              *
c     ****************************************************************

      subroutine mm10_fjac( solve_work, rjac, ldr, x, n )
c      
      use iso_Fortran_env
      use mm10_defs
      implicit none
c      
      include 'include_mm10'
      type(mm10_working_data) :: solve_work
      integer :: n, ldr
      double precision, dimension(n) :: x
      double precision :: rjac(ldr,n)
c!DIR$ ASSUME_ALIGNED x:64, rjac:64
c
      if( solve_work%solvfnc == 1 ) then
        if( solve_work%props%tang_calc == 2 ) then ! complex difference
          call mm10_formJ11i( solve_work%props, solve_work%np1, 
     &       solve_work%np0, solve_work%ivec1, solve_work%ivec2,
     &       x(1), solve_work%x2, rjac )
        else
          call mm10_formvecs( solve_work%props, solve_work%np1,
     &          solve_work%np0, x(1), solve_work%x2, 
     &          solve_work%vec1, solve_work%vec2 )
          call mm10_formarrs( solve_work%props, solve_work%np1, 
     &       solve_work%np0, x(1), solve_work%x2, 
     &       solve_work%vec1, solve_work%vec2,
     &       solve_work%arr1, solve_work%arr2,1 )
          call mm10_formJ11( solve_work%props, solve_work%np1,
     &        solve_work%np0, solve_work%vec1, solve_work%vec2,
     &        solve_work%arr1, solve_work%arr2, x(1), solve_work%x2,
     &        rjac)
        end if
      else
        if( solve_work%props%tang_calc == 2) then ! complex difference
          call mm10_formJi( solve_work%props, solve_work%np1, 
     &         solve_work%np0, solve_work%ivec1, solve_work%ivec2,
     &         x(1),x(7), rjac)
        else
          call mm10_formvecs( solve_work%props, solve_work%np1,
     &                     solve_work%np0, x(1), x(7), 
     &                     solve_work%vec1, solve_work%vec2)
          call mm10_formJ( solve_work%props, solve_work%np1, 
     &           solve_work%np0, solve_work%vec1, solve_work%vec2,
     &           solve_work%arr1, solve_work%arr2, x(1), x(7), rjac)
        end if
      end if
c
      return
      end
      
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine mm10_fvec                    *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 10/23/2016 rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_fvec( solve_work, x, fz, n, j )
      use iso_Fortran_env
      use mm10_defs
      implicit none
c      
      include 'include_mm10'
      type(mm10_working_data) :: solve_work
c      
      integer :: n, j
      double precision, dimension(n) :: x, fz
c!DIR$ ASSUME_ALIGNED x:64, fz:64
c
      if( solve_work%solvfnc == 1 ) then
        call mm10_formvecs( solve_work%props, solve_work%np1,
     &                   solve_work%np0, x(1:6), solve_work%x2, 
     &                     solve_work%vec1, solve_work%vec2 )
        call mm10_formR1( solve_work%props, solve_work%np1, 
     &                   solve_work%np0, solve_work%vec1, 
     &                   solve_work%vec2, x(1),solve_work%x2, fz, 
     &                   solve_work%gaspt )
      else
        call mm10_formR( solve_work%props, solve_work%np1,
     &                   solve_work%np0, solve_work%vec1,
     &                   solve_work%vec2, x(1),
     &                   x(7), fz, solve_work%gaspt )
        if( solve_work%locdebug ) then
          write(*,*) 'R', fz(9:12)
          write(*,*) 'str', x(1:6)
          write(*,*) 'tt', x(9:13)
        end if
      end if

      return
      end
      

c     ****************************************************************
c     *                                                              *
c     *                      subroutine mm10_fvec                    *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified: 11/3/2016 rhd               *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm10_fveci( solve_work, x, fz, n, j )
c
      use iso_Fortran_env
      use mm10_defs
      use mm10_constants
      implicit none
c      
      integer :: n, j
      complex(kind=real64), dimension(n) :: x, fz

      include 'include_mm10'
      type(mm10_working_data) :: solve_work
c
      integer :: length
      complex(kind=real64), dimension (solve_work%props%num_hard) :: x2i !  was dimension(n-6) :: x2i
      double precision, dimension(solve_work%props%num_hard) :: zero2 ! was dimension(n-6) :: zero2
c!DIR$ ASSUME_ALIGNED x:64, fz:64
      if( solve_work%solvfnc == 1 ) then
        zero2 = zero
        length = solve_work%props%num_hard
        x2i(1:length) = cmplx(solve_work%x2(1:length),zero2)
        call mm10_formvecsi( solve_work%props, solve_work%np1,
     &                    solve_work%np0, x(1), x2i, 
     &                    solve_work%ivec1, solve_work%ivec2 )
        call mm10_formR1i( solve_work%props, solve_work%np1, 
     &                     solve_work%np0, solve_work%ivec1, 
     &                     solve_work%ivec2, x(1), x2i, fz, 
     &                     solve_work%gaspt)
      else
        call mm10_formvecsi( solve_work%props, solve_work%np1,
     &                       solve_work%np0, x(1), x(7), 
     &                       solve_work%ivec1, solve_work%ivec2)
        call mm10_formR1i( solve_work%props, solve_work%np1, 
     &                     solve_work%np0, solve_work%ivec1, 
     &                     solve_work%ivec2, x(1),x(7), fz(1), 
     &                     solve_work%gaspt)
        call mm10_formR2i( solve_work%props, solve_work%np1, 
     &                     solve_work%np0, solve_work%ivec1, 
     &                     solve_work%ivec2, x(1), x(7), fz(7) )
        if( solve_work%locdebug ) then
          write(*,*) 'R', fz(9:12)
          write(*,*) 'str', x(1:6)
          write(*,*) 'tt', x(9:13)
        end if
      end if

      return
      end
c-----------------------------------------------------------------------------

      subroutine mm10_liqrfa(a, lda, n, tau, work, wsiz, info)
      integer  lda, n, wsiz, info
      double precision  a(lda,*), tau(*), work(*)

c-------------------------------------------------------------
c
c     QR decomposition of A(n,n)
c
c     Arguments
c
c      Inout A        Real(Lda,n)    Matrix to transform.
c      In    lda      Integer        Leading dimension of A
c      In    n        Integer        number of rows/cols A
c      Out   tau      Real(n)        Information for recovering
c      Out   work     Real(*)        workspace
c      In    wsiz     Integer        size of work()

c     Lapack blocked QR (much faster for larger n)
c
c-------------------------------------------------------------

      call dgeqrf(n,n,A,lda,tau,work,wsiz,info)

      return
      end

c=============================================================

      subroutine mm10_liqsiz(n,wrksiz)
      integer n, wrksiz

c-------------------------------------------------------------------------
c     Query the size of the double precision work array required
c     for optimal operation of the Lapack QR routines
c-------------------------------------------------------------------------

      double precision A(1), work(1)
      integer lwork, info

      lwork = -1
      call dgeqrf(n,n,A,n,work,work,lwork,info)
      if( info .ne. 0 ) then
          wrksiz = -1
      else
          wrksiz = int(work(1))
      endif

      return
      end

c=============================================================

      subroutine mm10_liqrqt(a, lda, n, tau, qty, work, wsiz, info)
      integer lda, n, wsiz, info
      double precision a(lda,*), tau(*), qty(*), work(*)

c-------------------------------------------------------------
c      Arguments
c
c      In    A     Real(Lda, n)    QR decomposition
c      In    Lda   Integer         Leading dimension A
c      In    n     Integer         Number of rows/columns in A
c      In    tau   Integer         Householder constants from QR
c      Inout qty   Real(n)         On input, vector y
c                                  On output, trans(Q)*y
c      Out   work  Real(*)         workspace
c      In    wsiz  Integer         size of work()
c
c     Liqrqt calculates trans(Q)*y from the QR decomposition
c
c     Lapack blocked
c-------------------------------------------------------------

      call dormqr('L','T',n,1,n,A,lda,tau,qty,n,work,wsiz,info)
      return
      end

c=============================================================

      subroutine mm10_liqrqq(q,ldq,tau,n,work,wsiz,info)
      integer n, ldq, wsiz, info
      double precision  q(ldq,*),tau(*),work(*)

c     Arguments
c
c     Inout  Q     Real(ldq,*)     On input, QR decomposition
c                                    On output, the full Q
c     In     ldq   Integer         leading dimension of Q
c     In     tau   Real(n)         Householder constants of
c                                     the QR decomposition
c     In     n     Integer         number of rows/columns in Q
c     Out    work  Real(n)         workspace of length n
c     In     wsiz  Integer         size of work()
c
c     Generate Q from QR decomposition Liqrfa (dgeqr2)
c
c     Lapack blocked
c-------------------------------------------------------------

      call dorgqr(n,n,n,q,ldq,tau,work,wsiz,info)

      return
      end

c-----------------------------------------------------------------------------

      subroutine mm10_nuzero( n, x )
      integer :: n
      double precision  :: x(*)

c     Parameters:
c
c     In    n        Integer           Number of elements.
c     In    x        Real(*)           Vector of reals.
c
c     Description:
c
c     Nuzero sets all elements of x to 0.
c     Does nothing when n <= 0

      double precision, parameter :: Rzero = 0.0d0
      integer i
c
      x(1:n) = Rzero
c
      return
      end

      subroutine mm10_limhpar(R, ldr, n, sdiag, qtf, dn, dnlen,
     *                   glen, delta, mu, d, work)
      integer ldr, n
      double precision R(ldr,*), sdiag(*)
      double precision qtf(*),dn(*), dnlen, glen, d(*),work(*)
      double precision delta, mu

c-----------------------------------------------------------------------------
c
c     Arguments
c
c       Inout  R      Real(ldr,*)     upper triangular matrix R from QR (unaltered)
c                                     strict lower triangle contains
c                                        transposed strict upper triangle of the upper
c                                        triangular matrix S.
c
c       In     n      Integer         order of R.
c
c       In     ldr    Integer         leading dimension of the array R.
c
c       Out    sdiag  Real(*)         vector of size n, containing the
c                                     diagonal elements of the upper
c                                     triangular matrix S.
c
c       In     qtr    Real(*)         trans(Q)*f vector of size n
c       In     dn     Real(*)         Newton step
c       In     dnlen  Real            length Newton step
c       In     glen   Real            length gradient vector
c
c       Inout  mu     Real            Levenberg-Marquardt parameter
c       In     delta  Real            size of trust region (euclidian norm)
c
c       Out    d      Real(*)         vector with step with norm very close to delta
c       Out    work   Real(*)         workspace of size n.
c
c     Description
c
c     determine Levenberg-Marquardt parameter mu such that
c     norm[(R**T R + mu*I)**(-1) * qtf] - delta approximately 0
c     See description in liqrev.f for further details
c
c     Algorithm comes from More: The Levenberg-Marquardt algorithm, Implementation and Theory
c     Lecture Notes in Mathematics, 1978, no. 630.
c     uses liqrev (in file liqrev.f) which is based on Nocedal's method (see comments in file)
c-----------------------------------------------------------------------------

      double precision phi, pnorm, qnorm, mulo, muhi,dmu, sqmu
      integer iter
      logical done
      double precision dnrm2

      double precision Rone
      parameter(Rone=1.0D0)

      phi = dnlen - delta
      muhi = glen/delta

      call dcopy(n,dn,1,d,1)
      call dscal(n, Rone/dnlen, d, 1)

c     solve R**T * x = dn
      call dtrsv("U","T","N",n,R,ldr,d,1)
      qnorm = dnrm2(n,d,1)
      mulo = (phi/dnlen)/qnorm**2
      mu = mulo

      iter = 0
      done = .false.
      do while( .not. done )
          iter = iter + 1
          sqmu = sqrt(mu)
          call mm10_liqrev(n, R, ldr, sqmu, qtf, d, sdiag, work)
          pnorm = dnrm2(n,d,1)
          call dcopy(n,d,1,work,1)
          call mm10_dtrstt(R, ldr, n, sdiag, work)
          done = abs(pnorm-delta) .le. .1d0*delta .or. iter .gt. 5
          if( .not. done ) then
              qnorm = dnrm2(n,work,1)
              if( pnorm .gt. delta ) then
                  mulo = max(mulo,mu)
              else if( pnorm .lt. delta ) then
                  muhi = min(muhi,mu)
              endif
              dmu = (pnorm-delta)/delta * (pnorm/qnorm)**2
              mu = max(mulo, mu + dmu)
          endif
      enddo
      return
      end

c-----------------------------------------------------------------------------

      subroutine mm10_liqrev(n,r,ldr,diag,b,x,sdiag,wrk)
      integer n,ldr
      double precision  r(ldr,*),b(*),x(*),sdiag(*),wrk(*)
      double precision diag

c-----------------------------------------------------------------------------
c
c     Arguments
c
c       In     n      Integer         order of R.
c       Inout  R      Real(ldr,*)     upper triangular matrix R from QR
c                                     unaltered
c                                     strict lower triangle contains
c                                        transposed strict upper triangle of the upper
c                                        triangular matrix S.
c
c       In     diag   Real            scalar for matrix D
c
c       In     ldr    Integer         leading dimension of the array R.
c       In     b      Real(*)         vector of size n
c
c       Out    x      Real(*)         vector of size n
c                                     on output contains least squares solution
c                                     of the system R*x = b, D*x = 0.
c
c       Out    sdiag  Real(*)         vector of size n, containing the
c                                     diagonal elements of the upper
c                                     triangular matrix S.
c
c       Out    wrk    Real(*)         workspace of size n.
c
c     Description
c
c     Given an n by n upper triangular matrix R, a diagonal matrix D with positive entries
c     and an n-vector b, determine an x which solves the system
c
c         |R*x| = |b|
c         |D*x| = |0|
c
c     in the least squares sense where D=diag*I.
c     This routine can be used for two different purposes.
c     The first is to provide a method of slightly modifying a singular or ill-conditioned matrix.
c     The second is for calculating a least squares solution to the above problem within
c     the context of e.g. a Levenberg-Marquardt algorithm combined with a More-Hebden algorithm
c     to determine a value of D (diagonal mu) such that x has a predetermined 2-norm.
c
c     The routine could also be used when the matrix R from the QR decomposition of a Jacobian
c     is ill-conditioned (or singular). Then it is difficult to calculate a Newton step
c     accurately (Dennis and Schnabel). D&S advise perturbing trans(J)*J with a positive
c     diagonal matrix.
c
c     The idea is  to solve (J^T * J + mu*I)x=b where mu is a small positive number.
c     Calculation of mu must be done in the calling routine.
c     Using a QR decomposition of J solving this system
c     is equivalent solving (R^T*R + mu*I)x=b, where R comes from the QR decomposition.
c     Solving this system is equivalent to solving the above least squares problem with the
c     elements of the matrix D set to sqrt(mu) which should be done in the calling routine.
c
c     On output the routine also provides an upper triangular matrix S such that
c     (see description of arguments above for the details)
c
c         (trans(R)*R + D*D) = trans(S)*S .
c
c     Method used here is described in
c     Nocedal and Wright, 2006, Numerical Optimization, Springer, ISBN 978-0-387-30303-1
c     page 258--261 (second edition)
c-----------------------------------------------------------------------------

      integer j,k
      double precision  bj,c,s,sum,temp
      double precision  ddot
      double precision Rzero
      parameter(Rzero=0.0d0)

c     copy R and b to preserve input and initialise S.
c     Save the diagonal elements of R in wrk.
c     Beware: the algorithm operates on an upper triangular matrix,
c     which is stored in lower triangle of R.
c
      do j=1,n
         call dcopy(n-j+1,r(j,j),ldr,r(j,j),1)
         wrk(j) = r(j,j)
      enddo
      call dcopy(n,b,1,x,1)

c     eliminate the diagonal matrix D using givens rotations.
c     Nocedal method: start at the bottom right
c     at end of loop R contains diagonal of S
c     save in sdiag and restore original diagonal of R

      do j=n,1,-1

c        initialise the row of D to be eliminated

         call mm10_nuzero(n-j+1,sdiag(j))
         sdiag(j) = diag

c        the transformations to eliminate the row of D

         bj = Rzero
         do k=j,n

c           determine a givens rotation which eliminates the
c           appropriate element in the current row of D.
c           accumulate the transformation in the row of S.

c           eliminate the diagonal element in row j of D
c           this generates fill-in in columns [j+1 .. n] of row j of D
c           successively eliminate the fill-in with givens rotations
c           for R[j+1,j+1] and D[j,j+1].
c           rows of R have been copied into the columns of R initially (see above)
c           perform all operations on those columns to preserve the original R

            if (sdiag(k) .ne. Rzero) then

               call mm10_nuvgiv(r(k,k),sdiag(k),c,s)
               if( k .lt. n ) then
                   call drot(n-k,r(k+1,k),1,sdiag(k+1),1,c,s)
               endif

c              compute the modified element of (b,0).

               temp =  c*x(k) + s*bj
               bj   = -s*x(k) + c*bj
               x(k) = temp

            endif

         enddo

      enddo

c     retrieve diagonal of S from diagonal of R
c     restore original diagonal of R

      do k=1,n
         sdiag(k) = r(k,k)
         r(k,k) = wrk(k)
      enddo

c     x now contains modified b
c     solve trans(S)*x = x
c     still to be done: guard against division by 0 to be absolutely safe
c     call dblepr('liqrev sdiag', 12, sdiag, n)
      x(n) = x(n) / sdiag(n)
      do j=n-1,1,-1
         sum  = ddot(n-j,r(j+1,j),1,x(j+1),1)
         x(j) = (x(j) - sum)/sdiag(j)
      enddo

      return
      end

c ----------------------------------------------------------------------

      subroutine mm10_dtrstt(S,ldr,n,sdiag,x)
      integer ldr, n
      double precision S(ldr,*), sdiag(*), x(*)
      integer j
      double precision sum, ddot

c     solve S*x = x where x is the result from subroutine liqrev
c     S is a lower triangular matrix with diagonal entries in sdiag()
c     and here it is in the lower triangular part of R as returned by liqrev

      x(1) = x(1) / sdiag(1)
      do j=2,n
         sum  = ddot(j-1,S(j,1),n,x,1)
         x(j) = (x(j) - sum)/sdiag(j)
      enddo

      return
      end

c ----------------------------------------------------------------------

      subroutine mm10_nuvgiv(x,y,c,s)
      double precision x,y,c,s

c     Parameters
c
c     Inout   x     Real       x input / c*x+s*y on output
c     Inout   y     Real       y input / 0       on output
c     Out     c     Real       c of tranformation (cosine)
c     Out     s     Real       s of tranformation (  sine)
c
c     Description
c
c     Nuvgiv calculates the givens rotator
c
c             |  c   s |
c         G = |        |
c             | -s   c |
c
c     with  c*c+s*s=1
c
c     for which G * | x | = | z |
c                   | y |   | 0 |
c
c     resulting in
c
c            c * x + s * y = z
c           -s * x + c * y = 0   ==>  s/c = y/x or c/s = x/y
c
c     Use Lapack dlartg routine
c     return c and s and the modified x and y

      double precision t

      double precision Rzero
      parameter(Rzero=0.0d0)

      call dlartg(x,y,c,s,t)
      x = t
      y = Rzero
      return
      end
      subroutine mm10_liqrup(Q,ldq,n,R,ldr,u,v,wk)
      integer ldq,n,ldr
      double precision Q(ldq,*),R(ldr,*),u(*),v(*),wk(*)

c-----------------------------------------------------------------------------
c
c     Arguments
c
c     Inout  Q       Real(ldq,n)      orthogonal matrix from QR
c     In     ldq     Integer          leading dimension of Q
c     In     n       Integer          order of Q and R
c     Inout  R       Real(ldr,n)      upper triangular matrix R from QR
c     In     ldr     Integer          leading dimension of R
c     In     u       Real(*)          vector u of size n
c     In     v       Real(*)          vector v of size n
c     Out    wk      Real(*)          workspace of size n
c
c     on return
c
c        Q       Q is the matrix with orthonormal columns in a QR
c                decomposition of the matrix B = A + u*v'
c
c        R       R is the upper triangular matrix in a QR
c                decomposition of the matrix B = A + u*v'
c
c     Description
c
c     The matrices Q and R are a QR decomposition of a square matrix
c     A = Q*R.
c     Given Q and R, qrupdt computes a QR decomposition of the rank one
c     modification B = A + u*trans(v) of A. Here u and v are vectors.
c
c     Source : procedure outlined in Dennis & Schnabel (Appendix A)
c              Algorithm 3.1.4 and 3.4.1a
c              modified (to use Lapack routines and more)
c
c-----------------------------------------------------------------------------

c     Local variables and functions

      integer k,i
      double precision  ddot

c     calculate wk = trans(Q)*u

      do i=1,n
         wk(i) = ddot(n,Q(1,i),1,u,1)
      enddo

c     zero components wk(n),wk(n-1)...wk(2)
c     and apply rotators to R and Q.

      do k=n-1,1,-1
         call mm10_jacrot(wk(k),wk(k+1),k,n,Q,ldq,R,ldr,k)
      enddo

c     r(1,1:n) += wk(1)*v(1:n)
      call daxpy(n,wk(1),v,1,R(1,1),ldr)

c     R is of upper hessenberg form. Triangularize R.
c      kr argument == k+1 to start applying rotation at column k+1
c      otherwise R(k,k) will be rotated twice and this way it also
c      avoids tiny roundoff errors.

      do k=1,n-1
         call mm10_jacrot(R(k,k),R(k+1,k),k,n,Q,ldq,R,ldr,k+1)
      enddo

      return
      end

c-----------------------------------------------------------------------------

      subroutine mm10_jacrot(a,b,k,n,Q,ldq,R,ldr,kr)

      double precision a, b
      integer k,n,ldr,ldq,kr
      double precision Q(ldq,*), R(ldr,*)

c-----------------------------------------------------------------------------
c
c     Arguments
c
c     Inout  a       Real             rotate argument
c     Inout  b       Real             rotate argument to rotate to zero
c     In     k       Integer          row/column number for rotation
c     In     n       Integer          order of Q and R
c     Inout  Q       Real(ldq,n)      orthogonal matrix from QR
c     In     ldq     Integer          leading dimension of Q
c     Inout  R       Real(ldr,n)      upper triangular matrix R from QR
c     In     ldr     Integer          leading dimension of R
c     In     u       Real(*)          vector u of size n
c     In     v       Real(*)          vector v of size n
c     In     kr      Integer          start R rotation in column kr
c                                     (should be k or k+1)
c
c-----------------------------------------------------------------------------

      double precision t
      double precision c,s
      double precision Rzero
      parameter(Rzero=0.0d0)

      call dlartg(a,b,c,s,t)
      a = t
      b = Rzero
      call drot(n-kr+1,R(k,kr),ldr,R(k+1,kr),ldr,c,s)
      call drot(n     ,Q(1,k) ,1  ,Q(1,k+1) ,1  ,c,s)

      return
      end

      subroutine mm10_nwbjac(rjac,r,ldr,n,xc,fc,fq,solve_work,epsm,jacflg,
     *                  wrk1,wrk2,wrk3,priter,
     *                  xscalm,scalex,gp,cndtol,rcdwrk,icdwrk,dn,
     *                  qtf,rcond,qrwork,qrwsiz,njcnt,iter,fstjac,ierr)

c-----------------------------------------------------------------------
c
c     Compute Jacobian matrix in xc, fc
c     scale it, compute gradient in xc and generate QR decomposition
c     calculate Newton step
c
c     Arguments
c
c     Out      rjac    Real(ldr,*)     jacobian (n columns) and used for storing full Q from Q
c     Out      r       Real(ldr,*)     used for storing R from QR factorization
c     In       ldr     Integer         leading dimension of rjac
c     In       n       Integer         dimensions of problem
c     In       xc      Real(*)         initial estimate of solution
c     Inout    fc      Real(*)         function values f(xc)
c     Wk       fq      Real(*)         workspace
c     In       fjac    Name            name of routine to calculate jacobian
c                                      (optional)
c     In       fvec    Name            name of routine to calculate f()
c     In       epsm    Real            machine precision
c     In       jacflg  Integer(*)      jacobian flag array
c                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
c                                                  3: user supplied banded
c                                      jacflg[2]: number of sub diagonals or -1 if not banded
c                                      jacflg[3]: number of super diagonals or -1 if not banded
c                                      jacflg[4]: 1 if adjusting step allowed when
c                                                   singular or illconditioned
c     Wk       wrk1    Real(*)         workspace
c     Wk       wrk2    Real(*)         workspace
c     Wk       wrk3    Real(*)         workspace
c     In       xscalm  Integer         x scaling method
c                                        1 from column norms of first jacobian
c                                          increased if needed after first iteration
c                                        0 scaling user supplied
c     Inout    scalex  Real(*)         scaling factors x(*)
c     Out      gp      Real(*)         gradient at xp()
c     In       cndtol  Real            tolerance of test for ill conditioning
c     Wk       rcdwrk  Real(*)         workspace
c     Wk       icdwrk  Integer(*)      workspace
c     Out      dn      Real(*)         Newton step
c     Out      qtf     Real(*)         workspace for nwnstp
c     Out      rcond   Real            estimated inverse condition of R from QR
c     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
c     In       qrwsiz  Integer         size of qrwork
c     Out      njcnt   Integer         number of jacobian evaluations
c     In       iter    Integer         iteration counter (used in scaling)
c     Inout    fstjac  logical         .true. if initial jacobian is available
c                                      on exit set to .false.
c     Out      ierr    Integer         error code
c                                        0 no error
c                                       >0 error in nwnstp (singular ...)
c
c-----------------------------------------------------------------------
      use iso_Fortran_env
      use mm10_defs
      integer ldr,n,iter, njcnt, ierr,priter
      integer jacflg(5),xscalm,qrwsiz
      logical fstjac
      double precision  epsm, cndtol, rcond
      double precision  rjac(ldr,*),r(ldr,*)
      double precision  xc(*),fc(*),dn(*)
      double precision  wrk1(n),wrk2(n),wrk3(*)
      double precision  qtf(*),gp(*),fq(*)
      double precision  scalex(*)
      double precision  rcdwrk(*),qrwork(*)
      integer           icdwrk(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work

      logical stepadj
      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)

c     evaluate the jacobian at the current iterate xc

      if( .not. fstjac ) then
         call mm10_nwfjac(xc,scalex,fc,fq,n,epsm,jacflg,solve_work,rjac,
     *               ldr,wrk1,wrk2,wrk3,priter)
         njcnt = njcnt + 1
      else
         fstjac = .false.
      endif

c     if requested calculate x scale from jacobian column norms a la Minpack

      if( xscalm .eq. 1 ) then
         call mm10_vunsc(n,xc,scalex)
         call mm10_nwcpsx(n,rjac,ldr,scalex,epsm,iter)
         call mm10_vscal(n,xc,scalex)
      endif

      call mm10_nwscjac(n,rjac,ldr,scalex)

c     evaluate the gradient at the current iterate xc
c     gp = trans(Rjac) * fc
      call dgemv('T',n,n,Rone,rjac,ldr,fc,1,Rzero,gp,1)

c     get broyden (newton) step
      stepadj = jacflg(4) .eq. 1
      call dcopy(n,fc,1,fq,1)
      call mm10_brdstp(rjac,r,ldr,fq,n,cndtol, stepadj,
     *            wrk1,dn,qtf,ierr,rcond,
     *            rcdwrk,icdwrk,qrwork,qrwsiz)

c     save some data about jacobian for later output
c      call mm10_nwsnot(0,ierr,rcond)

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_brdstp(rjac,r,ldr,fn,n,cndtol, stepadj,
     *                  qraux,dn,qtf,ierr,rcond,
     *                  rcdwrk,icdwrk,qrwork,qrwsiz)

      integer ldr,n,ierr,qrwsiz
      double precision  cndtol,rjac(ldr,*),r(ldr,*),qraux(*),fn(*)
      double precision  dn(*),qtf(*)
      double precision  rcdwrk(*),qrwork(*)
      integer           icdwrk(*)
      double precision  rcond
      logical           stepadj

c-----------------------------------------------------------------------
c
c     Calculate the newton step
c
c     Arguments
c
c     Inout    rjac    Real(ldr,*)     jacobian matrix at current iterate; on return full Q
c     Inout    r       Real(ldr,*)     jacobian matrix at current iterate; on return R fom QR
c     In       ldr     Integer         leading dimension of rjac
c     In       fn      Real(*)         function values at current iterate
c     In       n       Integer         dimension of problem
c     In       cndtol  Real            tolerance of test for ill conditioning
c     In       stepadj Logical         allow adjusting step for singular/illconditioned jacobian
c     Inout    qraux   Real(*)         QR info from liqrfa (calling Lapack dgeqrf)
c     Out      dn      Real(*)         Newton direction
c     Out      qtf     Real(*)         trans(Q)*f()
c     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
c                                      1 indicating Jacobian ill-conditioned
c                                      2 indicating Jacobian completely singular
c                                      3 indicating almost zero LM correction
c     Out      rcond   Real            inverse condition of upper triangular R of QR
c     Wk       rcdwrk  Real(*)         workspace
c     Wk       icdwrk  Integer(*)      workspace
c     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
c     In       qrwsiz  Integer         size of qrwork
c
c-----------------------------------------------------------------------

      integer info,k

      double precision Rone
      parameter(Rone=1.0d0)
      double precision mu

c     perform a QR factorization of rjac (simple Lapack routine)
c     check for singularity or ill conditioning
c     form qtf = trans(Q) * fn

      call mm10_liqrfa(rjac,ldr,n,qraux,qrwork,qrwsiz,ierr)
c     check for singularity or ill conditioning

      call mm10_cndjac(n,rjac,ldr,cndtol,rcond,rcdwrk,icdwrk,ierr)

c     compute qtf = trans(Q)*fn

      call dcopy(n,fn,1,qtf,1)
      call mm10_liqrqt(rjac, ldr, n, qraux, qtf, qrwork, qrwsiz, info)

c     copy the upper triangular part of the QR decomposition
c     contained in Rjac into R[*, 1..n].
c     form Q from the QR decomposition (taur/qraux in wrk1)

      call dlacpy('U',n,n,rjac,ldr,r,ldr)
      call mm10_liqrqq(rjac,ldr,qraux,n,qrwork,qrwsiz,info)

c     now Rjac[* ,1..n] holds expanded Q
c     now R[* ,1..n] holds full upper triangle R

      if( ierr .eq. 0 ) then
c         Normal Newton step
c         solve Jacobian*dn  =  -fn
c         ==> R*dn = - qtf

          call dcopy(n,qtf,1,dn,1)
          call dtrsv('U','N','N',n,r,ldr,dn,1)
          call dscal(n, -Rone, dn, 1)

      elseif( stepadj ) then
c         Adjusted Newton step
c         approximately from pseudoinverse(Jac+)
c         use mu to solve (trans(R)*R + mu*I*mu*I) * x = - trans(R) * fn
c         directly from the QR decomposition of R stacked with mu*I
c         a la Levenberg-Marquardt
          call mm10_compmu(r,ldr,n,mu,rcdwrk,ierr)
          if( ierr .eq. 0 ) then
             call mm10_liqrev(n,r,ldr,mu,qtf,dn,
     *                   rcdwrk(1+n),rcdwrk(2*n+1))
             call dscal(n, -Rone, dn, 1)

c            copy lower triangular R to upper triangular
             do k=1,n
                call dcopy (n-k+1,r(k,k),1,r(k,k),ldr)
                r(k,k) = rcdwrk(1+n+k-1)
             enddo
          endif
      endif

      return
      end


      subroutine mm10_brsolv(ldr,xc,n,scalex,maxit,
     *                  jacflg,xtol,ftol,btol,cndtol,global,xscalm,
     *                  stepmx,delta,sigma,
     *                  rjac,r,wrk1,wrk2,wrk3,wrk4,fc,fq,dn,d,qtf,
     *                  rcdwrk,icdwrk,qrwork,qrwsiz,epsm,
     *                  solve_work,outopt,xp,fp,gp,njcnt,nfcnt,iter,
     *                  termcd)
      use iso_Fortran_env
      use mm10_defs
      integer ldr,n,termcd,njcnt,nfcnt,iter
      integer maxit,jacflg(5),global,xscalm,qrwsiz
      integer outopt(*)
      double precision  xtol,ftol,btol,cndtol
      double precision  stepmx,delta,sigma,fpnorm,epsm
      double precision  rjac(ldr,*),r(ldr,*)
      double precision  xc(*),fc(*),xp(n),fp(n),dn(*),d(*)
      double precision  wrk1(n),wrk2(n),wrk3(*),wrk4(*)
      double precision  qtf(*),gp(*),fq(*)
      double precision  scalex(*)
      double precision  rcdwrk(*),qrwork(*)
      integer           icdwrk(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work

c-----------------------------------------------------------------------
c
c     Solve system of nonlinear equations with Broyden and global strategy
c
c
c     Arguments
c
c     In       ldr     Integer         leading dimension of rjac
c     In       xc      Real(*)         initial estimate of solution
c     In       n       Integer         dimensions of problem
c     Inout    scalex  Real(*)         scaling factors x(*)
c     In       maxit   Integer         maximum number of allowable iterations
c     In       jacflg  Integer(*)      jacobian flag array
c                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
c                                                  3: user supplied banded
c                                      jacflg[2]: number of sub diagonals or -1 if not banded
c                                      jacflg[3]: number of super diagonals or -1 if not banded
c                                      jacflg[4]: 1 if adjusting step allowed when
c                                                   singular or illconditioned
c     In       xtol    Real            tolerance at which successive iterates x()
c                                      are considered close enough to
c                                      terminate algorithm
c     In       ftol    Real            tolerance at which function values f()
c                                      are considered close enough to zero
c     Inout    btol    Real            x tolerance for backtracking
c     Inout    cndtol  Real            tolerance of test for ill conditioning
c     In       global  Integer         global strategy to use
c                                        1 cubic linesearch
c                                        2 quadratic linesearch
c                                        3 geometric linesearch
c                                        4 double dogleg
c                                        5 powell dogleg
c                                        6 hookstep (More-Hebden Levenberg-Marquardt)
c     In       xscalm  Integer         x scaling method
c                                        1 from column norms of first jacobian
c                                          increased if needed after first iteration
c                                        0 scaling user supplied
c     In       stepmx  Real            maximum allowable step size
c     In       delta   Real            trust region radius
c     In       sigma   Real            reduction factor geometric linesearch
c     Inout    rjac    Real(ldr,*)     jacobian (n columns)(compact QR decomposition/Q matrix)
c     Inout    r       Real(ldr,*)     stored R from QR decomposition
c     Wk       wrk1    Real(*)         workspace
c     Wk       wrk2    Real(*)         workspace
c     Wk       wrk3    Real(*)         workspace
c     Wk       wrk4    Real(*)         workspace
c     Inout    fc      Real(*)         function values f(xc)
c     Wk       fq      Real(*)         workspace
c     Wk       dn      Real(*)         workspace
c     Wk       d       Real(*)         workspace
c     Wk       qtf     Real(*)         workspace
c     Wk       rcdwrk  Real(*)         workspace
c     Wk       icdwrk  Integer(*)      workspace
c     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
c     In       qrwsiz  Integer         size of qrwork
c     In       epsm    Real            machine precision
c     In       fjac    Name            name of routine to calculate jacobian
c                                      (optional)
c     In       fvec    Name            name of routine to calculate f()
c     In       outopt  Integer(*)      output options
c     Out      xp      Real(*)         final x()
c     Out      fp      Real(*)         final f(xp)
c     Out      gp      Real(*)         gradient at xp()
c     Out      njcnt   Integer         number of jacobian evaluations
c     Out      nfcnt   Integer         number of function evaluations
c     Out      iter    Integer         number of (outer) iterations
c     Out      termcd  Integer         termination code
c
c-----------------------------------------------------------------------

      integer gcnt,retcd,ierr
      double precision  dum(2),dlt0,fcnorm,rcond
      logical fstjac
      logical jacevl,jacupd
      logical stepadj
      integer priter

      integer idamax

      double precision Rone
      parameter(Rone=1.0d0)

c     initialization

      retcd = 0
      iter  = 0
      njcnt = 0
      nfcnt = 0
      ierr  = 0

      dum(1) = 0
      dlt0 = delta

      if( outopt(1) .eq. 1 ) then
         priter = 1
      else
         priter = -1
      endif

c     evaluate function

      call mm10_vscal(n,xc,scalex)
      call mm10_nwfvec(xc,n,scalex,solve_work,fc,fcnorm,wrk1)

c     evaluate user supplied or finite difference jacobian and check user supplied
c     jacobian, if requested

      fstjac = .false.
      if(mod(jacflg(1),2) .eq. 1) then

        if( outopt(2) .eq. 1 ) then
           fstjac = .true.
           njcnt = njcnt + 1
           call mm10_nwfjac(xc,scalex,fc,fq,n,epsm,jacflg,solve_work,rjac,
     *                 ldr,wrk1,wrk2,wrk3,priter)
c      write(*,*) 'solve_work%solvfnc brsolv', solve_work%solvfnc
           call mm10_chkjac(rjac,ldr,xc,fc,n,epsm,jacflg,scalex,
     *                 fq,wrk1,wrk2,solve_work,termcd)
           if(termcd .lt. 0) then
c              copy initial values
               call dcopy(n,xc,1,xp,1)
               call dcopy(n,fc,1,fp,1)
               call mm10_vunsc(n,xp,scalex)
               fpnorm = fcnorm
               return
           endif
        endif

      endif

c     check stopping criteria for input xc

      call mm10_nwtcvg(xc,fc,xc,xtol,retcd,ftol,iter,maxit,n,ierr,
     &     termcd)

      if(termcd .gt. 0) then
          call dcopy(n,xc,1,xp,1)
          call dcopy(n,fc,1,fp,1)
          fpnorm = fcnorm
          if( outopt(3) .eq. 1 .and. .not. fstjac ) then
             njcnt = njcnt + 1
             call mm10_nwfjac(xp,scalex,fp,fq,n,epsm,jacflg,solve_work,
     *                   rjac,ldr,wrk1,wrk2,wrk3,priter)
          endif
          call mm10_vunsc(n,xp,scalex)
          return
      endif

      if( priter .gt. 0 ) then

         dum(1) = fcnorm
         dum(2) = abs(fc(idamax(n,fc,1)))

         if( global .eq. 0 ) then
            call mm10_nwprot(iter, -1, dum)
         elseif( global .le. 3 ) then
            call mm10_nwlsot(iter,-1,dum)
         elseif( global .eq. 4 ) then
            call mm10_nwdgot(iter,-1,0,dum)
         elseif( global .eq. 5 ) then
            call mm10_nwpwot(iter,-1,0,dum)
         elseif( global .eq. 6 ) then
            call mm10_nwmhot(iter,-1,0,dum)
         endif

      endif

      jacevl  = .true.
      stepadj = jacflg(4) .eq. 1
      if((maxit.gt.100).or.(maxit.lt.0)) then
        maxit = 100
      endif
      !write(*,*) 'maxit',maxit
      do while( termcd .eq. 0 )
         iter = iter+1

         if( jacevl ) then
            call mm10_nwbjac(rjac,r,ldr,n,xc,fc,fq,solve_work,epsm,jacflg,
     *                  wrk1,wrk2,wrk3,priter,
     *                  xscalm,scalex,gp,cndtol,rcdwrk,icdwrk,dn,
     *                  qtf,rcond,qrwork,qrwsiz,njcnt,iter,fstjac,ierr)

         else

c          - get broyden step
c          - calculate approximate gradient

            call dcopy(n,fc,1,fq,1)
            call mm10_brodir(rjac,ldr,r,fq,n,cndtol, stepadj,
     *                  dn,qtf,ierr,rcond,rcdwrk,icdwrk)

            if( ierr .eq. 0 ) then
               call dcopy(n,qtf,1,gp,1)
               call dtrmv('U','T','N',n,r,ldr,gp,1)
            endif
         endif
c      - choose the next iterate xp by a global strategy

         if( ierr .gt. 0 ) then
c           jacobian singular or too ill-conditioned
            call mm10_nweset(n,xc,fc,fcnorm,xp,fp,fpnorm,gcnt,priter,
     &         iter)
         elseif(global .eq. 0) then
            call mm10_nwpure(n,xc,dn,stepmx,scalex,
     *                  solve_work,xp,fp,fpnorm,wrk1,retcd,gcnt,
     *                  priter,iter)
         elseif(global .eq. 1) then
            call mm10_nwclsh0(n,xc,fcnorm,dn,gp,stepmx,btol,scalex,
     *                  solve_work,xp,fp,fpnorm,wrk1,retcd,gcnt,
     *                  priter,iter)
         elseif(global .eq. 2) then
            call mm10_nwqlsh(n,xc,fcnorm,dn,gp,stepmx,btol,scalex,
     *                  solve_work,xp,fp,fpnorm,wrk1,retcd,gcnt,
     *                  priter,iter)
         elseif(global .eq. 3) then
            call mm10_nwglsh(n,xc,fcnorm,dn,gp,sigma,stepmx,btol,
     *                  scalex,solve_work,xp,fp,fpnorm,wrk1,retcd,
     *                  gcnt,priter,iter)
         elseif(global .eq. 4) then
            call mm10_nwddlg(n,r,ldr,dn,gp,xc,fcnorm,stepmx,
     *                  btol,delta,qtf,scalex,
     *                  solve_work,d,fq,wrk1,wrk2,wrk3,wrk4,
     *                  xp,fp,fpnorm,retcd,gcnt,priter,iter)
         elseif(global .eq. 5) then
            call mm10_nwpdlg(n,r,ldr,dn,gp,xc,fcnorm,stepmx,
     *                  btol,delta,qtf,scalex,
     *                  solve_work,d,fq,wrk1,wrk2,wrk3,wrk4,
     *                  xp,fp,fpnorm,retcd,gcnt,priter,iter)
         elseif(global .eq. 6) then
            call mm10_nwmhlm(n,r,ldr,dn,gp,xc,fcnorm,stepmx,
     *                  btol,delta,qtf,scalex,
     *                  solve_work,d,fq,wrk1,wrk2,wrk3,wrk4,
     *                  xp,fp,fpnorm,retcd,gcnt,priter,iter)
         endif

         nfcnt = nfcnt + gcnt

c      - check stopping criteria for the new iterate xp

         call mm10_nwtcvg(xp,fp,xc,xtol,retcd,ftol,iter,maxit,n,ierr,
     &        termcd)

         if( termcd .eq. 3 .and. .not. jacevl ) then
c           global strategy failed but jacobian is out of date
c           try again with proper jacobian
c           reset trust region radius

            jacevl = .true.
            jacupd = .false.
            delta = dlt0
            termcd = 0

         elseif(termcd .gt. 0) then
            jacupd = .false.
         else
            jacupd = .true.
            jacevl = .false.
         endif

         if( jacupd ) then
c           perform Broyden update of current jacobian
c           update xc, fc, and fcnorm
            call mm10_brupdt(n,rjac,r,ldr,xc,xp,fc,fp,epsm,
     *                  wrk1,wrk2,wrk3)
            call dcopy(n,xp,1,xc,1)
            call dcopy(n,fp,1,fc,1)
            fcnorm = fpnorm
         endif

      enddo

      if( outopt(3) .eq. 1 ) then
c        final update of jacobian
         call mm10_brupdt(n,rjac,r,ldr,xc,xp,fc,fp,epsm,
     *               wrk1,wrk2,wrk3)
c        reconstruct Broyden matrix
c        calculate Q * R where Q is overwritten by result
c        Q is in rjac and R is in r
         call dtrmm('R','U','N','N',n,n,Rone,r,n,rjac,n)
c        unscale
         call mm10_nwunscjac(n,rjac,ldr,scalex)
      endif

      call mm10_vunsc(n,xp,scalex)

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_brupdt(n,q,r,ldr,xc,xp,fc,fp,epsm,dx,df,wa)
      integer n,ldr
      double precision  q(ldr,*),r(ldr,*)
      double precision  xc(*),xp(n),fc(*),fp(n),dx(*),df(*),wa(*)
      double precision  epsm

c-----------------------------------------------------------------------
c
c     Calculate new Q and R from rank-1 update with xp-xc and fp-fc
c     using Broyden method
c
c     Arguments
c
c     In       n       Integer         size of xc() etc.
c     Inout    Q       Real(ldr,n)     orthogonal matrix Q from QR
c                                       On output updated Q
c     Inout    R       Real(ldr,n)     upper triangular R from QR
c                                       On output updated R
c     In       ldr     Integer         leading dimension of Q and R
c     In       xc      Real(*)         current x() values
c     In       xp      Real(*)         new     x() values
c     In       fc      Real(*)         current f(xc)
c     In       fp      Real(*)         new     f(xp)
c     In       epsm    Real            machine precision
c     Wk       dx      Real(*)         workspace
c     Wk       df      Real(*)         workspace
c     Wk       wa      Real(*)         workspace
c
c-----------------------------------------------------------------------

      integer i
      double precision  eta,sts
      double precision  dnrm2
      logical doupdt

      double precision Rzero, Rone, Rtwo, Rhund
      parameter(Rzero=0.0d0, Rone=1.0d0, Rtwo=2.0d0, Rhund=100d0)

      eta    = Rhund * Rtwo * epsm
      doupdt = .false.

      do i=1,n
         dx(i) = xp(i) - xc(i)
         df(i) = fp(i) - fc(i)
      enddo

c     clear lower triangle

      do i=1,n-1
         call mm10_nuzero(n-i,r(i+1,i))
      enddo

c     calculate df - B*dx = df - Q*R*dx
c     wa = R*dx
c     df = df - Q*(R*dx) (!not really needed if qrupdt were to be changed)
c     do not update with noise

      call dcopy(n,dx,1,wa,1)
      call dtrmv('U','N','N',n,r,ldr,wa,1)
      call dgemv('N',n,n,-Rone,q,ldr,wa,1,Rone,df,1)

      do i=1,n
         if( abs(df(i)) .gt. eta*( abs(fp(i)) + abs(fc(i)) ) ) then
            doupdt = .true.
         else
            df(i)  = Rzero
         endif
      enddo

      if( doupdt ) then
c        equation 8.3.1 from Dennis and Schnabel (page 187)(Siam edition)
         sts = dnrm2(n,dx,1)
         call dscal(n,Rone/sts,dx,1)
         call dscal(n,Rone/sts,df,1)
         call mm10_liqrup(q,ldr,n,r,ldr,df,dx,wa)
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_brodir(q,ldr,r,fn,n,cndtol,stepadj,dn,qtf,
     *                  ierr,rcond,rcdwrk,icdwrk)

      integer ldr,n,ierr
      double precision  cndtol,q(ldr,*),r(ldr,*),fn(*)
      double precision  dn(*),qtf(*)
      double precision  rcdwrk(*)
      integer           icdwrk(*)
      double precision  rcond
      logical           stepadj

c-----------------------------------------------------------------------
c
c     Calculate the approximate newton direction
c
c     Arguments
c
c     Inout    Q       Real(ldr,*)     Q part from QR at current iterate
c     In       ldr     Integer         leading dimension of Q and R
c     In       R       Real(ldr,*)     upper triangular R from QR decomposition
c     In       fn      Real(*)         function values at current iterate
c     In       n       Integer         dimension of problem
c     In       cndtol  Real            tolerance of test for ill conditioning
c     In       stepadj Logical         allow adjusting step for singular/illconditioned jacobian
c     Out      dn      Real(*)         Newton direction
c     Out      qtf     Real(*)         trans(Q)*f()
c     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
c                                      1 indicating Jacobian ill-conditioned
c                                      2 indicating Jacobian completely singular
c                                      3 indicating almost zero LM correction
c     Out      rcond   Real            inverse condition of matrix
c     Wk       rcdwrk  Real(*)         workspace
c     Wk       icdwrk  Integer(*)      workspace
c
c     QR decomposition with no pivoting.
c
c-----------------------------------------------------------------------

      integer k
      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)
      double precision mu

c     check for singularity or ill conditioning

      call mm10_cndjac(n,r,ldr,cndtol,rcond,rcdwrk,icdwrk,ierr)

      if( ierr .eq. 0 ) then
c         form qtf = trans(Q) * fn

          call dgemv('T',n,n,Rone,q,ldr,fn,1,Rzero,qtf,1)

c         solve rjac*dn  =  -fn
c         ==> R*dn = - qtf

          call dcopy(n,qtf,1,dn,1)
          call dtrsv('U','N','N',n,r,ldr,dn,1)
          call dscal(n, -Rone, dn, 1)

      elseif( stepadj ) then
c         call intpr('ierr brodir', 12,ierr,1)
c         Adjusted Newton step
c         approximately from pseudoinverse(Jac+)
c         compute qtf = trans(Q)*fn

c         form qtf = trans(Q) * fn

          call dgemv('T',n,n,Rone,q,ldr,fn,1,Rzero,qtf,1)

c         use mu to solve (trans(R)*R + mu*I*mu*I) * x = - trans(R) * fn
c         directly from the QR decomposition of R stacked with mu*I
c         a la Levenberg-Marquardt
          call mm10_compmu(r,ldr,n,mu,rcdwrk,ierr)
          if( ierr .eq. 0 ) then
             call mm10_liqrev(n,r,ldr,mu,qtf,dn,
     *                   rcdwrk(1+n),rcdwrk(2*n+1))
             call dscal(n, -Rone, dn, 1)

c            copy lower triangular Rjac to upper triangular
             do k=1,n
                call dcopy (n-k+1,r(k,k),1,r(k,k),ldr)
                r(k,k) = rcdwrk(1+n+k-1)
             enddo
          endif
      endif
c      call mm10_nwsnot(1,ierr,rcond)

      return
      end

      subroutine mm10_nwclsh0(n,xc,fcnorm,d,g,stepmx,xtol,scalex,
     *             solve_work,xp,fp,fpnorm,xw,retcd,gcnt,priter,iter)
      use iso_Fortran_env
      use mm10_defs
      integer n,retcd,gcnt
      double precision  stepmx,xtol,fcnorm,fpnorm
      double precision  xc(*)
      double precision  d(*),g(*),xp(n),fp(n),xw(*)
      double precision  scalex(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work

      integer priter,iter

c-------------------------------------------------------------------------
c
c     Find a next acceptable iterate using a safeguarded cubic line search
c     along the newton direction
c
c     Arguments
c
c     In       n       Integer          dimension of problem
c     In       xc      Real(*)          current iterate
c     In       fcnorm  Real             0.5 * || f(xc) ||**2
c     In       d       Real(*)          newton direction
c     In       g       Real(*)          gradient at current iterate
c     In       stepmx  Real             maximum stepsize
c     In       xtol    Real             relative step size at which
c                                       successive iterates are considered
c                                       close enough to terminate algorithm
c     In       scalex  Real(*)          diagonal scaling matrix for x()
c     In       fvec    Name             name of routine to calculate f()
c     In       xp      Real(*)          new x()
c     In       fp      Real(*)          new f(x)
c     In       fpnorm  Real             .5*||fp||**2
c     Out      xw      Real(*)           workspace for unscaling x(*)
c
c     Out      retcd   Integer          return code
c                                         0 new satisfactory x() found
c                                         1 no  satisfactory x() found
c                                           sufficiently distinct from xc()
c
c     Out      gcnt    Integer          number of steps taken
c     In       priter  Integer           >0 unit if intermediate steps to be printed
c                                        -1 if no printing
c
c-------------------------------------------------------------------------

      integer i, gcntmax
      double precision  alpha,slope,rsclen,oarg(4)
      double precision  lambda,lamhi,lamlo,t
      double precision  ddot,dnrm2, mm10_nudnrm, ftarg
      double precision  dlen
      double precision a, b, disc, fpt, fpt0, fpnorm0, lambda0
      integer idamax
      logical firstback

      parameter (alpha = 1.0d-4)

      double precision Rhalf, Rone, Rtwo, Rthree, Rten, Rzero
      parameter(Rzero=0.0d0)
      parameter(Rhalf=0.5d0, Rone=1.0d0, Rtwo=2.0d0, Rten=10.0d0)
      parameter(Rthree=3.0d0)
      double precision t1,t2
      double precision dlamch

c     silence warnings issued by ftncheck

      lambda0 = Rzero
      fpnorm0 = Rzero

c     safeguard initial step size

      dlen = dnrm2(n,d,1)
      if( dlen .gt. stepmx ) then
          lamhi  = stepmx / dlen
      else
          lamhi  = Rone
      endif

c     compute slope  =  g-trans * d

      slope = ddot(n,g,1,d,1)

c     compute the smallest value allowable for the damping
c     parameter lambda ==> lamlo

      rsclen = mm10_nudnrm(n,d,xc)
      lamlo  = xtol / rsclen

c     initialization of retcd and lambda (linesearch length)

      retcd  = 2
      lambda = lamhi
      gcnt   = 0
      firstback = .true.
      gcntmax = 1000

      do while( (retcd .eq. 2) .and. (gcnt .le. gcntmax))

c        compute next x

         do i=1,n
            xp(i) = xc(i) + lambda*d(i)
         enddo

c        evaluate functions and the objective function at xp

         call mm10_nwfvec(xp,n,scalex,solve_work,fp,fpnorm,xw)
         gcnt = gcnt + 1
         ftarg = fcnorm + alpha * lambda * slope

         if( priter .gt. 0) then
            oarg(1) = lambda
            oarg(2) = ftarg
            oarg(3) = fpnorm
            oarg(4) = abs(fp(idamax(n,fp,1)))
            call mm10_nwlsot(iter,1,oarg)
         endif

c        first is quadratic
c        test whether the standard step produces enough decrease
c        of the objective function.
c        If not update lambda and compute a new next iterate

         if( fpnorm .le. ftarg ) then
             retcd = 0
         else
            if( fpnorm .gt. lamlo**2 * sqrt(dlamch('O')) ) then
c               safety against overflow in what follows (use lamlo**2 for safety)
                lambda = lambda/Rten
                firstback = .true.
            else
                if( firstback ) then
                   t = ((-lambda**2)*slope/Rtwo)/
     *                        (fpnorm-fcnorm-lambda*slope)
                   firstback = .false.
                else
                   fpt  = fpnorm -fcnorm - lambda*slope
                   fpt0 = fpnorm0-fcnorm - lambda0*slope
                   a = fpt/lambda**2 - fpt0/lambda0**2
                   b = -lambda0*fpt/lambda**2 + lambda*fpt0/lambda0**2
                   a = a /(lambda - lambda0)
                   b = b /(lambda - lambda0)
                   if( abs(a) .le. dlamch('E') ) then
c                      not quadratic but linear
                       t = -slope/(2*b)
                   else
c                      use Higham procedure to compute roots acccurately
c                      Higham: Accuracy and Stability of Numerical Algorithms, second edition,2002, page 10.    
c                      Actually solving 3*a*x^2+2*b*x+c=0 ==> (3/2)*a*x^2+b*x+(c/2)=0
                       disc = b**2 - Rthree * a * slope
                       t1 = -(b+sign(Rone,b)*sqrt(disc))/(Rthree*a)
                       t2 = slope/(Rthree*a)/t1
                       if(a .gt. Rzero ) then
c                          upward opening parabola ==> rightmost is solution
                           t = max(t1,t2)
                       else
c                          downward opening parabola ==> leftmost is solution
                           t = min(t1,t2)
                       endif
                   endif
                   t = min(t, Rhalf*lambda)
                endif
                lambda0 = lambda
                fpnorm0 = fpnorm
                lambda = max(t,lambda/Rten)
                if(lambda .lt. lamlo) then
                   retcd = 1
                endif
            endif
         endif
      enddo

      if(gcnt.ge.gcntmax) then
        retcd = 1
      endif

      return
      end

      subroutine mm10_nwddlg(n,rjac,ldr,dn,g,xc,fcnorm,stepmx,xtol,
     *                  delta,qtf,scalex,solve_work,d,xprev,
     *                  ssd,v,wa,fprev,xp,fp,fpnorm,retcd,gcnt,
     *                  priter,iter)
      use iso_Fortran_env
      use mm10_defs
      integer ldr, n, retcd, gcnt, priter, iter
      double precision  fcnorm, stepmx, xtol, fpnorm, delta
      double precision  rjac(ldr,*), dn(*), g(*), xc(*), qtf(*)
      double precision  scalex(*), d(*)
      double precision  xprev(*), xp(n), fp(n)
      double precision  ssd(*), v(*), wa(*), fprev(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work

c-------------------------------------------------------------------------
c
c     Find a next iterate xp by the double dogleg method
c
c     Arguments
c
c     In       n       Integer         size of problem: dimension x, f
c     In       Rjac    Real(ldr,*)     R of QR-factored jacobian
c     In       ldr     Integer         leading dimension of Rjac
c     Inout    dn      Real(*)         newton direction
c     Inout    g       Real(*)         gradient at current point
c                                      trans(jac)*f()
c     In       xc      Real(*)         current iterate
c     In       fcnorm  Real            .5*||f(xc)||**2
c     In       stepmx  Real            maximum stepsize
c     In       xtol    Real            x-tolerance (stepsize)
c     Inout    delta   Real            on input: initial trust region radius
c                                                if -1 then set to something
c                                                reasonable
c                                      on output: final value
c                                      ! Do not modify between calls while
c                                        still iterating
c     In       qtf     Real(*)         trans(Q)*f(xc)
c     In       scalex  Real(*)         scaling factors for x()
c     In       fvec    Name            name of subroutine to evaluate f(x)
c                                      ! must be declared external in caller
c     Wk       d       Real(*)         work vector
c     Wk       xprev   Real(*)         work vector
c     Wk       ssd     Real(*)         work vector
c     Wk       v       Real(*)         work vector
c     Wk       wa      Real(*)         work vector
c     Wk       fprev   Real(*)         work vector
c     Out      xp      Real(*)         new x()
c     Out      fp      Real(*)         new f(xp)
c     Out      fpnorm  Real            new .5*||f(xp)||**2
c     Out      retcd   Integer         return code
c                                       0  new satisfactory x() found
c                                       1  no  satisfactory x() found
c     Out      gcnt    Integer         number of steps taken
c     In       priter  Integer         print flag
c                                       -1 no intermediate printing
c                                       >0 yes for print of intermediate results
c     In       iter    Integer         current iteration (only used for above)
c
c     All vectors at least size n
c
c-------------------------------------------------------------------------

      integer i, gcntmax
      double precision  dnlen,ssdlen,alpha,beta,lambda,fpred
      double precision  sqalpha,eta,gamma,fpnsav,oarg(7)
      double precision  dnrm2,ddot, xpr(n)
      logical nwtstep
      integer dtype

      integer idamax

      double precision Rone, Rtwo, Rten, Rhalf, Rp2, Rp8
      parameter(Rhalf=0.5d0)
      parameter(Rone=1.0d0, Rtwo=2.0d0, Rten=10.0d0)
      parameter(Rp2 = Rtwo/Rten, Rp8 = Rone - Rp2)
      double precision Rzero
      parameter(Rzero=0.0d0)

c     length newton direction

      dnlen = dnrm2(n, dn, 1)

c     steepest descent direction and length

      sqalpha = dnrm2(n,g,1)
      alpha   = sqalpha**2

      call dcopy(n, g, 1, d, 1)
      call dtrmv('U','N','N',n,rjac,ldr,d,1)
      beta = dnrm2(n,d,1)**2

      call dcopy(n, g, 1, ssd, 1)
      call dscal(n, -(alpha/beta), ssd, 1)

      ssdlen = alpha*sqalpha/beta

c     set trust radius to ssdlen or dnlen if required

      if( delta .eq. -Rone ) then
         delta = min(ssdlen, stepmx)
      elseif( delta .eq. -Rtwo ) then
         delta = min(dnlen, stepmx)
      endif

c     calculate double dogleg parameter

      gamma = alpha*alpha/(-beta*ddot(n,g,1,dn,1))
c      call dgdbg(gamma, alpha*alpha, -beta*ddot(n,g,1,dn,1))
c     precautionary (just in case)
      eta = max(Rzero, min(Rone,Rp2 + Rp8*gamma))

      retcd = 4
      gcnt  = 0
      gcntmax = 1000

      do while( retcd .gt. 1 .and. (gcnt .le. gcntmax))
c        find new step by double dogleg algorithm

         call mm10_ddlgstp(n,dn,dnlen,delta,v,
     *               ssd,ssdlen,eta,d,dtype,lambda)
         nwtstep = dtype .eq. 4
c        compute the model prediction 0.5*||F + J*d||**2 (L2-norm)

         call dcopy(n,d,1,wa,1)
         call dtrmv('U','N','N',n,rjac,ldr,wa,1)
         call daxpy(n, Rone, qtf,1,wa,1)
         fpred = Rhalf * dnrm2(n,wa,1)**2

c        evaluate function at xp = xc + d

         do i=1,n
            xp(i) = xc(i) + d(i)
         enddo

         call mm10_nwfvec(xp,n,scalex,solve_work,fp,fpnorm,wa)
         gcnt = gcnt + 1

c        check whether the global step is acceptable

         oarg(2) = delta
         call mm10_nwtrup(n,fcnorm,g,d,nwtstep,stepmx,xtol,delta,
     *               fpred,retcd,xprev,fpnsav,fprev,xp,fp,fpnorm)

         if( priter .gt. 0 ) then
c            xpr(1:n) = xp(1:n)
c               call mm10_vunsc(n,xpr,scalex)
c            write(*,*) 'xp', xpr(1:n)
c            write(*,*) 'fp', fp(1:n)
            oarg(1) = lambda
            oarg(3) = delta
            oarg(4) = eta
            oarg(5) = fpnorm
            oarg(6) = abs(fp(idamax(n,fp,1)))
            call mm10_nwdgot(iter,dtype,retcd,oarg)
         endif

      enddo

      if(gcnt.ge.gcntmax) then
        retcd = 1
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_ddlgstp(n,dn,dnlen,delta,v,
     *                  ssd,ssdlen,eta,d,dtype,lambda)
      integer n
      double precision  dn(*), ssd(*), v(*), d(*)
      double precision  dnlen, delta, ssdlen, eta, lambda
      integer dtype

c-------------------------------------------------------------------------
c
c     Find a new step by the double dogleg algorithm
c     Internal routine for nwddlg
c
c     Arguments
c
c     In       n       Integer         size of problem
c     In       dn      Real(*)         current newton step
c     Out      dnlen   Real            length dn()
c     In       delta   Real            current trust region radius
c     Out      v       Real(*)         (internal) eta * dn() - ssd()
c     In       ssd     Real(*)         (internal) steepest descent direction
c     In       ssdlen  Real            (internal) length ssd
c     In       eta     Real            (internal) double dogleg parameter
c     Out      d       Real(*)         new step for x()
c     Out      dtype   Integer         steptype
c                                       1 steepest descent
c                                       2 combination of dn and ssd
c                                       3 partial newton step
c                                       4 full newton direction
c     Out      lambda  Real            weight of eta*dn() in d()
c                                      closer to 1 ==> more of eta*dn()
c
c-----------------------------------------------------------------------

      integer i
      double precision vssd, vlen
      double precision dnrm2, ddot

      if(dnlen .le. delta) then

c        Newton step smaller than trust radius ==> take it

         call dcopy(n, dn, 1, d, 1)
         delta = dnlen
         dtype = 4

      elseif(eta*dnlen .le. delta) then

c        take partial step in newton direction

         call dcopy(n, dn, 1, d, 1)
         call dscal(n, delta / dnlen, d, 1)
         dtype = 3

      elseif(ssdlen .ge. delta) then

c        take step in steepest descent direction

         call dcopy(n, ssd, 1, d, 1)
         call dscal(n, delta / ssdlen, d, 1)
         dtype = 1

      else

c        calculate convex combination of ssd and eta*dn with length delta

         do i=1,n
            v(i) = eta*dn(i) - ssd(i)
         enddo

         vssd = ddot(n,v,1,ssd,1)
         vlen = dnrm2(n,v,1)**2

         lambda =(-vssd+sqrt(vssd**2-vlen*(ssdlen**2-delta**2)))/vlen
         call dcopy(n, ssd, 1, d, 1)
         call daxpy(n, lambda, v, 1, d, 1)
         dtype = 2

      endif

      return
      end

      subroutine mm10_nwglsh(n,xc,fcnorm,d,g,sigma,stepmx,xtol,scalex,
     *                  solve_work,xp,fp,fpnorm,xw,retcd,gcnt,priter,
     *                      iter)
      use iso_Fortran_env
      use mm10_defs
      integer n,retcd,gcnt
      double precision  sigma,stepmx,xtol,fcnorm,fpnorm
      double precision  xc(*)
      double precision  d(*),g(*),xp(n),fp(n),xw(*)
      double precision  scalex(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work

      integer priter,iter

c-------------------------------------------------------------------------
c
c     Find a next acceptable iterate using geometric line search
c     along the newton direction
c
c     Arguments
c
c     In       n       Integer          dimension of problem
c     In       xc      Real(*)          current iterate
c     In       fcnorm  Real             0.5 * || f(xc) ||**2
c     In       d       Real(*)          newton direction
c     In       g       Real(*)          gradient at current iterate
c     In       sigma   Real             reduction factor for lambda
c     In       stepmx  Real             maximum stepsize
c     In       xtol    Real             relative step size at which
c                                       successive iterates are considered
c                                       close enough to terminate algorithm
c     In       scalex  Real(*)          diagonal scaling matrix for x()
c     In       fvec    Name             name of routine to calculate f()
c     In       xp      Real(*)          new x()
c     In       fp      Real(*)          new f(x)
c     In       fpnorm  Real             .5*||fp||**2
c     Out      xw      Real(*)          workspace for unscaling x
c
c     Out      retcd   Integer          return code
c                                         0 new satisfactory x() found
c                                         1 no  satisfactory x() found
c                                           sufficiently distinct from xc()
c
c     Out      gcnt    Integer          number of steps taken
c     In       priter  Integer           >0 if intermediate steps to be printed
c                                        -1 if no printing
c
c-------------------------------------------------------------------------

      integer i, gcntmax
      double precision  alpha,slope,rsclen,oarg(4)
      double precision  lambda,lamhi,lamlo
      double precision  ddot,dnrm2, mm10_nudnrm, ftarg
      double precision  dlen,xpr(n)

      integer idamax

      parameter (alpha = 1.0d-4)

      double precision Rone
      parameter(Rone=1.0d0)

c     safeguard initial step size

      dlen = dnrm2(n,d,1)
      if( dlen .gt. stepmx ) then
          lamhi  = stepmx / dlen
      else
          lamhi  = Rone
      endif

c     compute slope  =  g-trans * d

      slope = ddot(n,g,1,d,1)

c     compute the smallest value allowable for the damping
c     parameter lambda ==> lamlo

      rsclen = mm10_nudnrm(n,d,xc)
      lamlo  = xtol / rsclen

c     initialization of retcd and lambda (linesearch length)

      retcd  = 2
      lambda = lamhi
      gcnt   = 0
      gcntmax = 1000

         if( priter .gt. 0) then
c            xpr(1:n) = d(1:n)
c               call mm10_vunsc(n,xpr,scalex)
c            write(*,*) 'lamlo', lamlo
c            write(*,*) 'sigma', sigma
         endif

      do while( (retcd .eq. 2) .and. (gcnt .le. gcntmax))

c        compute next x

         do i=1,n
            xp(i) = xc(i) + lambda*d(i)
         enddo

c        evaluate functions and the objective function at xp

         call mm10_nwfvec(xp,n,scalex,solve_work,fp,fpnorm,xw)
         gcnt = gcnt + 1
         ftarg = fcnorm + alpha * lambda * slope

         if( priter .gt. 0) then
c            xpr(1:n) = xp(1:n)
c               call mm10_vunsc(n,xpr,scalex)
c            write(*,*) 'xp', xpr(1:n)
c            xpr(1:n) = fp(1:n)
c            write(*,*) 'fp', xpr(1:n)
            oarg(1) = lambda
            oarg(2) = ftarg
            oarg(3) = fpnorm
            oarg(4) = abs(fp(idamax(n,fp,1)))
            call mm10_nwlsot(iter,1,oarg)
         endif

c        test if the standard step produces enough decrease
c        of the objective function.
c        If not update lambda and compute a new next iterate


         if( priter .gt. 0) then
c            xpr(1:n) = d(1:n)
c               call mm10_vunsc(n,xpr,scalex)
c            write(*,*) 'lambda', lambda
         endif
         if( fpnorm .le. ftarg ) then
            retcd = 0
         else
            lambda  = sigma * lambda
            if(lambda .lt. lamlo) then
               retcd = 1
            endif
         endif

      enddo

      if(gcnt.ge.gcntmax) then
        retcd = 1
      endif

      return
      end

      subroutine mm10_nwmhlm(n,rjac,ldr,dn,g,xc,fcnorm,stepmx,xtol,
     *                  delta,qtf,scalex,solve_work,d,xprev,
     *                  ssd,v,wa,fprev,xp,fp,fpnorm,retcd,gcnt,
     *                  priter,iter)
      use iso_Fortran_env
      use mm10_defs
      integer ldr, n, retcd, gcnt, priter, iter
      double precision  fcnorm, stepmx, xtol, fpnorm, delta
      double precision  rjac(ldr,*), dn(*), g(*), xc(*), qtf(*)
      double precision  scalex(*), d(*)
      double precision  xprev(*), xp(n), fp(n)
      double precision  ssd(*), v(*), wa(*), fprev(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work

c-------------------------------------------------------------------------
c
c     Find a next iterate xp by the More-Hebden-Levenberg-Marquardt method
c
c     Arguments
c
c     In       n       Integer         size of problem: dimension x, f
c     In       Rjac    Real(ldr,*)     R of QR-factored jacobian
c     In       ldr     Integer         leading dimension of Rjac
c     Inout    dn      Real(*)         newton direction
c     Inout    g       Real(*)         gradient at current point
c                                      trans(jac)*f()
c     In       xc      Real(*)         current iterate
c     In       fcnorm  Real            .5*||f(xc)||**2
c     In       stepmx  Real            maximum stepsize
c     In       xtol    Real            x-tolerance (stepsize)
c     Inout    delta   Real            on input: initial trust region radius
c                                                if -1 then set to something
c                                                reasonable
c                                      on output: final value
c                                      ! Do not modify between calls while
c                                        still iterating
c     In       qtf     Real(*)         trans(Q)*f(xc)
c     In       scalex  Real(*)         scaling factors for x()
c     In       fvec    Name            name of subroutine to evaluate f(x)
c                                      ! must be declared external in caller
c     Wk       d       Real(*)         work vector
c     Wk       xprev   Real(*)         work vector
c     Wk       ssd     Real(*)         work vector
c     Wk       v       Real(*)         work vector
c     Wk       wa      Real(*)         work vector
c     Wk       fprev   Real(*)         work vector
c     Out      xp      Real(*)         new x()
c     Out      fp      Real(*)         new f(xp)
c     Out      fpnorm  Real            new .5*||f(xp)||**2
c     Out      retcd   Integer         return code
c                                       0  new satisfactory x() found
c                                       1  no  satisfactory x() found
c     Out      gcnt    Integer         number of steps taken
c     In       priter  Integer         print flag
c                                       -1 no intermediate printing
c                                       >0 yes for print of intermediate results
c     In       iter    Integer         current iteration (only used for above)
c
c     All vectors at least size n
c
c-------------------------------------------------------------------------

      integer i, gcntmax
      double precision  dnlen,glen,ssdlen,alpha,beta,mu,fpred
      double precision  fpnsav,oarg(6)
      double precision  dnrm2
      logical nwtstep
      integer dtype

      integer idamax

      double precision Rone, Rtwo, Rhalf
      parameter(Rhalf=0.5d0)
      parameter(Rone=1.0d0, Rtwo=2.0d0)
      gcntmax = 100

c     length newton direction

      dnlen = dnrm2(n, dn, 1)

c     gradient length and steepest descent direction and length

      glen  = dnrm2(n,g,1)
      alpha = glen**2

      call dcopy(n, g, 1, d, 1)
      call dtrmv('U','N','N',n,rjac,ldr,d,1)
      beta = dnrm2(n,d,1)**2

      call dcopy(n, g, 1, ssd, 1)
      call dscal(n, -(alpha/beta), ssd, 1)

      ssdlen = alpha*glen/beta

c     set trust radius to ssdlen or dnlen if required

      if( delta .eq. -Rone ) then
         delta = min(ssdlen, stepmx)
      elseif( delta .eq. -Rtwo ) then
         delta = min(dnlen, stepmx)
      endif

      retcd = 4
      gcnt  = 0

      do while( (retcd .gt. 1) .and. (gcnt.le.gcntmax) )
c        find new step by More Hebden LM  algorithm
c        reuse ssd as sdiag

         call mm10_nwmhstep(Rjac,ldr,n,ssd,qtf,dn,dnlen,glen,delta,mu,
     *                 d, v, dtype)
         nwtstep = dtype .eq. 2
c        compute the model prediction 0.5*||F + J*d||**2 (L2-norm)

         call dcopy(n,d,1,wa,1)
         call dtrmv('U','N','N',n,rjac,ldr,wa,1)
         call daxpy(n, Rone, qtf,1,wa,1)
         fpred = Rhalf * dnrm2(n,wa,1)**2

c        evaluate function at xp = xc + d

         do i=1,n
            xp(i) = xc(i) + d(i)
         enddo

         call mm10_nwfvec(xp,n,scalex,solve_work,fp,fpnorm,wa)
         gcnt = gcnt + 1
c         write(*,*) gcnt
c        check whether the global step is acceptable

         oarg(2) = delta
         call mm10_nwtrup(n,fcnorm,g,d,nwtstep,stepmx,xtol,delta,
     *               fpred,retcd,xprev,fpnsav,fprev,xp,fp,fpnorm)

         if( priter .gt. 0 ) then
            oarg(1) = mu
            oarg(3) = delta
            oarg(4) = dnrm2(n, d, 1)
            oarg(5) = fpnorm
            oarg(6) = abs(fp(idamax(n,fp,1)))
            call mm10_nwmhot(iter,dtype,retcd,oarg)
         endif

      enddo
      if(gcnt.ge.gcntmax) then
        retcd = 1
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwmhstep(R,ldr,n,sdiag,qtf,dn,dnlen,glen,
     *                  delta,mu,d, work, dtype)
      integer ldr, n
      double precision R(ldr,*)
      double precision sdiag(*), qtf(*), dn(*), d(*), work(*)
      double precision  dnlen, glen, delta, mu
      integer dtype

c-------------------------------------------------------------------------
c
c     Find a new step by the More Hebden Levemberg Marquardt algorithm
c     Internal routine for nwmhlm
c
c     Arguments
c
c     In       R       Real(ldr,*)     R of QR-factored jacobian
c     In       ldr     Integer         leading dimension of R
c     In       n       Integer         size of problem
c     Out      sdiag   Real(*)         diagonal of LM lower triangular modified R
c     In       qtf     Real(*)         trans(Q)*f(xc)
c     In       dn      Real(*)         current newton step
c     Out      dnlen   Real            length dn()
c     In       glen    Real            length gradient
c     In       delta   Real            current trust region radius
c     Inout    mu      Real            Levenberg-Marquardt parameter
c     Out      d       Real(*)         new step for x()
c     Work     work    Real(*)         work vector for limhpar
c     Out      dtype   Integer         steptype
c                                       1 LM step
c                                       2 full newton direction
c
c-----------------------------------------------------------------------

      double precision Rone
      parameter(Rone=1.0D0)

      if(dnlen .le. delta) then

c        Newton step smaller than trust radius ==> take it

         call dcopy(n, dn, 1, d, 1)
         delta = dnlen
         dtype = 2

      else

c        calculate LM step
         call mm10_limhpar(R, ldr, n, sdiag, qtf, dn, dnlen, glen, 
     *                delta,mu, d, work)
c        change sign of step d (limhpar solves for trans(R)*R+mu*I)=qtf instead of -qtf)
         call dscal(n,-Rone,d,1)
         dtype = 1
      endif

      return
      end
      subroutine mm10_nwnjac(rjac,ldr,n,xc,fc,fq,solve_work,epsm,jacflg,
     *                  wrk1,wrk2,wrk3,priter,
     *                  xscalm,scalex,gp,cndtol,rcdwrk,icdwrk,dn,
     *                  qtf,rcond,qrwork,qrwsiz,njcnt,iter,fstjac,ierr)

c-----------------------------------------------------------------------
c
c     Compute Jacobian matrix in xc, fc
c     scale it, compute gradient in xc and generate QR decomposition
c     calculate Newton step
c
c     Arguments
c
c     Out      rjac    Real(ldr,*)     jacobian (n columns)
c     In       ldr     Integer         leading dimension of rjac
c     In       n       Integer         dimensions of problem
c     In       xc      Real(*)         initial estimate of solution
c     Inout    fc      Real(*)         function values f(xc)
c     Wk       fq      Real(*)         workspace
c     In       fjac    Name            name of routine to calculate jacobian
c                                      (optional)
c     In       fvec    Name            name of routine to calculate f()
c     In       epsm    Real            machine precision
c     In       jacflg  Integer(*)      jacobian flag array
c                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
c                                                  3: user supplied banded
c                                      jacflg[2]: number of sub diagonals or -1 if not banded
c                                      jacflg[3]: number of super diagonals or -1 if not banded
c                                      jacflg[4]: 1 if adjusting step allowed when
c                                                   singular or illconditioned
c     Wk       wrk1    Real(*)         workspace
c     Wk       wrk2    Real(*)         workspace
c     Wk       wrk3    Real(*)         workspace
c     In       xscalm  Integer         x scaling method
c                                        1 from column norms of first jacobian
c                                          increased if needed after first iteration
c                                        0 scaling user supplied
c     Inout    scalex  Real(*)         scaling factors x(*)
c     Out      gp      Real(*)         gradient at xp()
c     In       cndtol  Real            tolerance of test for ill conditioning
c     Wk       rcdwrk  Real(*)         workspace
c     Wk       icdwrk  Integer(*)      workspace
c     Out      dn      Real(*)         Newton step
c     Out      qtf     Real(*)         workspace for nwnstp
c     Out      rcond   Real            estimated inverse condition of R from QR
c     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
c     In       qrwsiz  Integer         size of qrwork
c     Out      njcnt   Integer         number of jacobian evaluations
c     In       iter    Integer         iteration counter (used in scaling)
c     Inout    fstjac  logical         .true. if initial jacobian is available
c                                      on exit set to .false.
c     Out      ierr    Integer         error code
c                                        0 no error
c                                       >0 error in nwnstp (singular ...)
c
c-----------------------------------------------------------------------
      use iso_Fortran_env
      use mm10_defs
      integer ldr,n,iter, njcnt, ierr
      integer jacflg(5),xscalm,qrwsiz,priter
      logical fstjac
      double precision  epsm, cndtol, rcond
      double precision  rjac(ldr,*)
      double precision  xc(*),fc(*),dn(*),xpr(n)
      double precision  wrk1(n),wrk2(n),wrk3(*)
      double precision  qtf(*),gp(*),fq(*)
      double precision  scalex(*)
      double precision  rcdwrk(*),qrwork(*)
      integer           icdwrk(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work

      logical stepadj
      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)

c     evaluate the jacobian at the current iterate xc

c      if( priter .gt. 0 ) then

c            xpr(1:n) = xc(1:n)
c               call mm10_vunsc(n,xpr,scalex)
c            write(*,*) 'xp in jac', xpr(1:n)

c      endif

      if( .not. fstjac ) then
         call mm10_nwfjac(xc,scalex,fc,fq,n,epsm,jacflg,solve_work,
     *               rjac,ldr,wrk1,wrk2,wrk3,priter)
         njcnt = njcnt + 1
      else
         fstjac = .false.
      endif
c      if( priter .gt. 0 ) then
c        write (*,*) "ldr", ldr
c        write (*,*) "J11", rjac(1:n,1:n)

c      endif

c     if requested calculate x scale from jacobian column norms a la Minpack

      if( xscalm .eq. 1 ) then
         call mm10_vunsc(n,xc,scalex)
         call mm10_nwcpsx(n,rjac,ldr,scalex,epsm,iter)
         call mm10_vscal(n,xc,scalex)
      endif

      call mm10_nwscjac(n,rjac,ldr,scalex)

c     evaluate the gradient at the current iterate xc
c     gp = trans(Rjac) * fc
      call dgemv('T',n,n,Rone,rjac,ldr,fc,1,Rzero,gp,1)

c     get newton step
      stepadj = jacflg(4) .eq. 1
      call dcopy(n,fc,1,fq,1)
      call mm10_nwnstp(rjac,ldr,fq,n,cndtol, stepadj,
     *            wrk1,dn,qtf,ierr,rcond,
     *            rcdwrk,icdwrk,qrwork,qrwsiz)

c     save some data about jacobian for later output
c      call mm10_nwsnot(0,ierr,rcond)

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwnstp(rjac,ldr,fn,n,cndtol, stepadj,
     *                  qraux,dn,qtf,ierr,rcond,
     *                  rcdwrk,icdwrk,qrwork,qrwsiz)

      integer ldr,n,ierr,qrwsiz
      double precision  cndtol,rjac(ldr,*),qraux(*),fn(*)
      double precision  dn(*),qtf(*)
      double precision  rcdwrk(*),qrwork(*)
      integer           icdwrk(*)
      double precision  rcond
      logical           stepadj

c-----------------------------------------------------------------------
c
c     Calculate the newton step
c
c     Arguments
c
c     Inout    rjac    Real(ldr,*)     jacobian matrix at current iterate
c                                      overwritten with QR decomposition
c     In       ldr     Integer         leading dimension of rjac
c     In       fn      Real(*)         function values at current iterate
c     In       n       Integer         dimension of problem
c     In       cndtol  Real            tolerance of test for ill conditioning
c     In       stepadj Logical         allow adjusting step for singular/illconditioned jacobian
c     Inout    qraux   Real(*)         QR info from liqrfa (calling Lapack dgeqrf)
c     Out      dn      Real(*)         Newton direction
c     Out      qtf     Real(*)         trans(Q)*f()
c     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
c                                      1 indicating Jacobian ill-conditioned
c                                      2 indicating Jacobian completely singular
c                                      3 indicating almost zero LM correction
c     Out      rcond   Real            inverse condition of upper triangular R of QR
c     Wk       rcdwrk  Real(*)         workspace
c     Wk       icdwrk  Integer(*)      workspace
c     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
c     In       qrwsiz  Integer         size of qrwork
c
c-----------------------------------------------------------------------

      integer info,k

      double precision Rone
      parameter(Rone=1.0d0)
      double precision mu

c     perform a QR factorization of rjac (simple Lapack routine)
c     check for singularity or ill conditioning
c     form qtf = trans(Q) * fn

      call mm10_liqrfa(rjac,ldr,n,qraux,qrwork,qrwsiz,ierr)

c     check for singularity or ill conditioning

      call mm10_cndjac(n,rjac,ldr,cndtol,rcond,rcdwrk,icdwrk,ierr)

c     compute qtf = trans(Q)*fn

      call dcopy(n,fn,1,qtf,1)
      call mm10_liqrqt(rjac, ldr, n, qraux, qtf, qrwork, qrwsiz, info)

      if( ierr .eq. 0 ) then
c         Normal Newton step
c         solve rjac*dn  =  -fn
c         ==> R*dn = - qtf

          call dcopy(n,qtf,1,dn,1)
          call dtrsv('U','N','N',n,rjac,ldr,dn,1)
          call dscal(n, -Rone, dn, 1)

      elseif( stepadj ) then
c         Adjusted Newton step
c         approximately from pseudoinverse(Jac+)
c         use mu to solve (trans(R)*R + mu*I*mu*I) * x = - trans(R) * fn
c         directly from the QR decomposition of R stacked with mu*I
c         a la Levenberg-Marquardt
          call mm10_compmu(rjac,ldr,n,mu,rcdwrk,ierr)
          if( ierr .eq. 0 ) then
             call mm10_liqrev(n,rjac,ldr,mu,qtf,dn,
     *                   rcdwrk(1+n),rcdwrk(2*n+1))
             call dscal(n, -Rone, dn, 1)

c            copy lower triangular Rjac to upper triangular
             do k=1,n
                call dcopy (n-k+1,rjac(k,k),1,rjac(k,k),ldr)
                rjac(k,k) = rcdwrk(1+n+k-1)
             enddo
          endif
      endif

      return
      end

      subroutine mm10_nwsolv(ldr,xc,n,scalex,maxit,
     *                  jacflg,xtol,ftol,btol,cndtol,global,xscalm,
     *                  stepmx,delta,sigma,
     *                  rjac,wrk1,wrk2,wrk3,wrk4,fc,fq,dn,d,qtf,
     *                  rcdwrk,icdwrk,qrwork,qrwsiz,epsm,
     *                  solve_work,outopt,xp,fp,gp,njcnt,nfcnt,iter,
     *                  termcd)
      use iso_Fortran_env
      use mm10_defs
      integer ldr,n,termcd,njcnt,nfcnt,iter
      integer maxit,jacflg(5),global,xscalm,qrwsiz
      integer outopt(*)
      double precision  xtol,ftol,btol,cndtol
      double precision  stepmx,delta,sigma,fpnorm,epsm
      double precision  rjac(ldr,*)
      double precision  xc(*),fc(*),xp(n),fp(n),dn(*),d(*)
      double precision  wrk1(n),wrk2(n),wrk3(*),wrk4(*)
      double precision  qtf(*),gp(*),fq(*)
      double precision  scalex(*)
      double precision  rcdwrk(*),qrwork(*)
      integer           icdwrk(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work

c-----------------------------------------------------------------------
c
c     Solve system of nonlinear equations with Newton and global strategy
c
c
c     Arguments
c
c     In       ldr     Integer         leading dimension of rjac
c     In       xc      Real(*)         initial estimate of solution
c     In       n       Integer         dimensions of problem
c     Inout    scalex  Real(*)         scaling factors x(*)
c     In       maxit   Integer         maximum number of allowable iterations
c     In       jacflg  Integer(*)      jacobian flag array
c                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
c                                                  3: user supplied banded
c                                      jacflg[2]: number of sub diagonals or -1 if not banded
c                                      jacflg[3]: number of super diagonals or -1 if not banded
c                                      jacflg[4]: 1 if adjusting step allowed when
c                                                   singular or illconditioned
c     In       xtol    Real            tolerance at which successive iterates x()
c                                      are considered close enough to
c                                      terminate algorithm
c     In       ftol    Real            tolerance at which function values f()
c                                      are considered close enough to zero
c     Inout    btol    Real            x tolerance for backtracking
c     Inout    cndtol  Real            tolerance of test for ill conditioning
c     In       global  Integer         global strategy to use
c                                        1 cubic linesearch
c                                        2 quadratic linesearch
c                                        3 geometric linesearch
c                                        4 double dogleg
c                                        5 powell dogleg
c                                        6 hookstep (More-Hebden Levenberg-Marquardt)
c     In       xscalm  Integer         x scaling method
c                                        1 from column norms of first jacobian
c                                          increased if needed after first iteration
c                                        0 scaling user supplied
c     In       stepmx  Real            maximum allowable step size
c     In       delta   Real            trust region radius
c     In       sigma   Real            reduction factor geometric linesearch
c     Inout    rjac    Real(ldr,*)     jacobian (n columns)
c     Wk       wrk1    Real(*)         workspace
c     Wk       wrk2    Real(*)         workspace
c     Wk       wrk3    Real(*)         workspace
c     Wk       wrk4    Real(*)         workspace
c     Inout    fc      Real(*)         function values f(xc)
c     Wk       fq      Real(*)         workspace
c     Wk       dn      Real(*)         workspace
c     Wk       d       Real(*)         workspace
c     Wk       qtf     Real(*)         workspace
c     Wk       rcdwrk  Real(*)         workspace
c     Wk       icdwrk  Integer(*)      workspace
c     In       qrwork  Real(*)         workspace for Lapack QR routines (call liqsiz)
c     In       qrwsiz  Integer         size of qrwork
c     In       epsm    Real            machine precision
c     In       fjac    Name            name of routine to calculate jacobian
c                                      (optional)
c     In       fvec    Name            name of routine to calculate f()
c     In       outopt  Integer(*)      output options
c     Out      xp      Real(*)         final x()
c     Out      fp      Real(*)         final f(xp)
c     Out      gp      Real(*)         gradient at xp()
c     Out      njcnt   Integer         number of jacobian evaluations
c     Out      nfcnt   Integer         number of function evaluations
c     Out      iter    Integer         number of (outer) iterations
c     Out      termcd  Integer         termination code
c
c-----------------------------------------------------------------------

      integer gcnt,retcd,ierr
      double precision  dum(2),fcnorm,rcond, xpr(n)
      logical fstjac
      integer priter

      integer idamax

c     initialization

      retcd = 0
      iter  = 0
      njcnt = 0
      nfcnt = 0
      ierr  = 0

      dum(1) = 0

      if( outopt(1) .eq. 1 ) then
         priter = 1
      else
         priter = -1
      endif

c     evaluate function

      call mm10_vscal(n,xc,scalex)
      call mm10_nwfvec(xc,n,scalex,solve_work,fc,fcnorm,wrk1)
c      if( priter .gt. 0 ) then
c            xpr(1:n) = fc(1:n)
c            write(*,*) 'fp', xpr(1:n)
c      endif

c     evaluate user supplied or finite difference jacobian and check user supplied
c     jacobian, if requested

      fstjac = .false.
      if(mod(jacflg(1),2) .eq. 1) then

        if( outopt(2) .eq. 1 ) then
           fstjac = .true.
           njcnt = njcnt + 1
           call mm10_nwfjac(xc,scalex,fc,fq,n,epsm,jacflg,solve_work,
     *                 rjac,ldr,wrk1,wrk2,wrk3,priter)
c      write(*,*) 'solve_work%solvfnc fnwsolv', solve_work%solvfnc
           call mm10_chkjac(rjac,ldr,xc,fc,n,epsm,jacflg,scalex,
     *                 fq,wrk1,wrk2,solve_work,termcd)
           if(termcd .lt. 0) then
c              copy initial values
               call dcopy(n,xc,1,xp,1)
               call dcopy(n,fc,1,fp,1)
               call mm10_vunsc(n,xp,scalex)
               fpnorm = fcnorm
               return
           endif
        endif

      endif

c     check stopping criteria for input xc

      call mm10_nwtcvg(xc,fc,xc,xtol,retcd,ftol,iter,maxit,n,ierr,
     &       termcd)

      if(termcd .gt. 0) then
          call dcopy(n,xc,1,xp,1)
          call dcopy(n,fc,1,fp,1)
          fpnorm = fcnorm
          if( outopt(3) .eq. 1 .and. .not. fstjac ) then
             njcnt = njcnt + 1
             call mm10_nwfjac(xp,scalex,fp,fq,n,epsm,jacflg,solve_work,
     *                   rjac,ldr,wrk1,wrk2,wrk3,priter)
          endif
          call mm10_vunsc(n,xp,scalex)
          return
      endif

      if( priter .gt. 0 ) then

         dum(1) = fcnorm
         dum(2) = abs(fc(idamax(n,fc,1)))

         if( global .eq. 0 ) then
            call mm10_nwprot(iter, -1, dum)
         elseif( global .le. 3 ) then
            call mm10_nwlsot(iter,-1,dum)
         elseif( global .eq. 4 ) then
            call mm10_nwdgot(iter,-1,0,dum)
         elseif( global .eq. 5 ) then
            call mm10_nwpwot(iter,-1,0,dum)
         elseif( global .eq. 6 ) then
            call mm10_nwmhot(iter,-1,0,dum)
         endif

      endif
      if((maxit.gt.100).or.(maxit.lt.0)) then
        maxit = 100
      endif
      !write(*,*) 'maxit',maxit

      do while( termcd .eq. 0 )
         iter = iter + 1

c      if( priter .gt. 0 ) then

c            xpr(1:n) = xc(1:n)
c               call mm10_vunsc(n,xpr,scalex)
c            write(*,*) 'xp', xpr(1:n)

c      endif

         call mm10_nwnjac(rjac,ldr,n,xc,fc,fq,solve_work,epsm,jacflg,
     *               wrk1,wrk2,wrk3,priter,
     *               xscalm,scalex,gp,cndtol,rcdwrk,icdwrk,dn,
     *               qtf,rcond,qrwork,qrwsiz,njcnt,iter,fstjac,ierr)
c        - choose the next iterate xp by a global strategy

      if( priter .gt. 0 ) then

      endif
         if( ierr .gt. 0 ) then
c           jacobian singular or too ill-conditioned
            call mm10_nweset(n,xc,fc,fcnorm,xp,fp,fpnorm,gcnt,priter,
     &    iter)
         elseif(global .eq. 0) then
            call mm10_nwpure(n,xc,dn,stepmx,scalex,
     *                  solve_work,xp,fp,fpnorm,wrk1,retcd,gcnt,
     *                  priter,iter)
         elseif(global .eq. 1) then
            call mm10_nwclsh0(n,xc,fcnorm,dn,gp,stepmx,btol,scalex,
     *                  solve_work,xp,fp,fpnorm,wrk1,retcd,gcnt,
     *                  priter,iter)
         elseif(global .eq. 2) then
            call mm10_nwqlsh(n,xc,fcnorm,dn,gp,stepmx,btol,scalex,
     *                  solve_work,xp,fp,fpnorm,wrk1,retcd,gcnt,
     *                  priter,iter)
         elseif(global .eq. 3) then
            call mm10_nwglsh(n,xc,fcnorm,dn,gp,sigma,stepmx,btol,
     *                  scalex,solve_work,xp,fp,fpnorm,wrk1,retcd,
     *                  gcnt,priter,iter)
         elseif(global .eq. 4) then
            call mm10_nwddlg(n,rjac,ldr,dn,gp,xc,fcnorm,stepmx,
     *                  btol,delta,qtf,scalex,
     *                  solve_work,d,fq,wrk1,wrk2,wrk3,wrk4,
     *                  xp,fp,fpnorm,retcd,gcnt,priter,iter)
         elseif(global .eq. 5) then
            call mm10_nwpdlg(n,rjac,ldr,dn,gp,xc,fcnorm,stepmx,
     *                  btol,delta,qtf,scalex,
     *                  solve_work,d,fq,wrk1,wrk2,wrk3,wrk4,
     *                  xp,fp,fpnorm,retcd,gcnt,priter,iter)
         elseif(global .eq. 6) then
            call mm10_nwmhlm(n,rjac,ldr,dn,gp,xc,fcnorm,stepmx,
     *                  btol,delta,qtf,scalex,
     *                  solve_work,d,fq,wrk1,wrk2,wrk3,wrk4,
     *                  xp,fp,fpnorm,retcd,gcnt,priter,iter)
         endif

         nfcnt = nfcnt + gcnt

c        - check stopping criteria for the new iterate xp

         call mm10_nwtcvg(xp,fp,xc,xtol,retcd,ftol,iter,maxit,n,ierr,
     &    termcd)

         if(termcd .eq. 0) then
c           update xc, fc, and fcnorm
            call dcopy(n,xp,1,xc,1)
            call dcopy(n,fp,1,fc,1)
            fcnorm = fpnorm
         endif

      enddo

      if( outopt(3) .eq. 1 ) then
         call mm10_nwfjac(xp,scalex,fp,fq,n,epsm,jacflg,solve_work,rjac,
     *               ldr,wrk1,wrk2,wrk3,priter)
      endif

      call mm10_vunsc(n,xp,scalex)

      return
      end

      subroutine mm10_nwpdlg(n,rjac,ldr,dn,g,xc,fcnorm,stepmx,xtol,
     *                  delta,qtf,scalex,solve_work,d,xprev,
     *                  ssd,v,wa,fprev,xp,fp,fpnorm,retcd,gcnt,
     *                  priter,iter)
      use iso_Fortran_env
      use mm10_defs
      integer ldr, n, retcd, gcnt, priter, iter
      double precision  fcnorm, stepmx, xtol, fpnorm, delta
      double precision  rjac(ldr,*), dn(*), g(*), xc(*), qtf(*)
      double precision  scalex(*), d(*)
      double precision  xprev(*), xp(n), fp(n)
      double precision  ssd(*), v(*), wa(*), fprev(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work

c-------------------------------------------------------------------------
c
c     Find a next iterate xp by the Powell dogleg method
c
c     Arguments
c
c     In       n       Integer         size of problem: dimension x, f
c     In       Rjac    Real(ldr,*)     R of QR-factored jacobian
c     In       ldr     Integer         leading dimension of Rjac
c     Inout    dn      Real(*)         newton direction
c     Inout    g       Real(*)         gradient at current point
c                                      trans(jac)*f()
c     In       xc      Real(*)         current iterate
c     In       fcnorm  Real            .5*||f(xc)||**2
c     In       stepmx  Real            maximum stepsize
c     In       xtol    Real            x-tolerance (stepsize)
c     Inout    delta     Real            on input: initial trust region radius
c                                                if -1 then set to something
c                                                reasonable
c                                      on output: final value
c                                      ! Do not modify between calls while
c                                        still iterating
c     In       qtf     Real(*)         trans(Q)*f(xc)
c     In       scalex  Real(*)         scaling factors for x()
c     In       fvec    Name            name of subroutine to evaluate f(x)
c                                      ! must be declared external in caller
c     Wk       d       Real(*)         work vector
c     Wk       xprev   Real(*)         work vector
c     Wk       ssd     Real(*)         work vector
c     Wk       v       Real(*)         work vector
c     Wk       wa      Real(*)         work vector
c     Wk       fprev   Real(*)         work vector
c     Out      xp      Real(*)         new x()
c     Out      fp      Real(*)         new f(xp)
c     Out      fpnorm  Real            new .5*||f(xp)||**2
c     Out      retcd   Integer         return code
c                                       0  new satisfactory x() found
c                                       1  no  satisfactory x() found
c     Out      gcnt    Integer         number of steps taken
c     In       priter  Integer         print flag
c                                       -1 no intermediate printing
c                                       >0 yes for print of intermediate results
c     In       iter    Integer         current iteration (only used for above)
c
c     All vectors at least size n
c
c-------------------------------------------------------------------------

      integer i, gcntmax
      double precision  dnlen,ssdlen,alpha,beta,lambda,fpred
      double precision  sqalpha,fpnsav,oarg(5)
      double precision  dnrm2
      logical nwtstep
      integer dtype

      integer idamax

      double precision Rone, Rtwo, Rhalf
      parameter(Rhalf=0.5d0)
      parameter(Rone=1.0d0, Rtwo=2.0d0)

c     length newton direction

      dnlen = dnrm2(n, dn, 1)

c     steepest descent direction and length

      sqalpha = dnrm2(n,g,1)
      alpha   = sqalpha**2

      call dcopy(n, g, 1, d, 1)
      call dtrmv('U','N','N',n,rjac,ldr,d,1)
      beta = dnrm2(n,d,1)**2

      call dcopy(n, g, 1, ssd, 1)
      call dscal(n, -(alpha/beta), ssd, 1)

      ssdlen = alpha*sqalpha/beta

c     set trust radius to ssdlen or dnlen if required

      if( delta .eq. -Rone ) then
         delta = min(ssdlen, stepmx)
      elseif( delta .eq. -Rtwo ) then
         delta = min(dnlen, stepmx)
      endif

      retcd = 4
      gcnt  = 0
      gcntmax = 1000

      do while( (retcd .gt. 1) .and. (gcnt .le. gcntmax))

c        find new step by single dogleg algorithm

         call mm10_pwlstp(n,dn,dnlen,delta,v,
     *               ssd,ssdlen,d,dtype,lambda)
         nwtstep = dtype .eq. 3
c        compute the model prediction 0.5*||F + J*d||**2 (L2-norm)

         call dcopy(n,d,1,wa,1)
         call dtrmv('U','N','N',n,rjac,ldr,wa,1)
         call daxpy(n, Rone, qtf,1,wa,1)
         fpred = Rhalf * dnrm2(n,wa,1)**2

c        evaluate function at xp = xc + d

         do i=1,n
            xp(i) = xc(i) + d(i)
         enddo

         call mm10_nwfvec(xp,n,scalex,solve_work,fp,fpnorm,wa)
         gcnt = gcnt + 1

c        check whether the global step is acceptable

         oarg(2) = delta
         call mm10_nwtrup(n,fcnorm,g,d,nwtstep,stepmx,xtol,delta,
     *               fpred,retcd,xprev,fpnsav,fprev,xp,fp,fpnorm)

         if( priter .gt. 0 ) then
            oarg(1) = lambda
            oarg(3) = delta
            oarg(4) = fpnorm
            oarg(5) = abs(fp(idamax(n,fp,1)))
            call mm10_nwpwot(iter,dtype,retcd,oarg)
         endif

      enddo
      if(gcnt.ge.gcntmax) then
        retcd = 1
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_pwlstp(n,dn,dnlen,delta,v,
     *                  ssd,ssdlen,d,dtype,lambda)
      integer n
      double precision  dn(*), ssd(*), v(*), d(*)
      double precision  dnlen, delta, ssdlen, lambda
      integer dtype

c-------------------------------------------------------------------------
c
c     Find a new step by the Powell dogleg algorithm
c     Internal routine for nwpdlg
c
c     Arguments
c
c     In       n       Integer         size of problem
c     In       dn      Real(*)         current newton step
c     Out      dnlen   Real            length dn()
c     In       delta   Real            current trust region radius
c     Out      v       Real(*)         (internal) dn() - ssd()
c     In       ssd     Real(*)         (internal) steepest descent direction
c     In       ssdlen  Real            (internal) length ssd
c     Out      d       Real(*)         new step for x()
c     Out      dtype   Integer         steptype
c                                       1 steepest descent
c                                       2 combination of dn and ssd
c                                       3 full newton direction
c     Out      lambda  Real            weight of dn() in d()
c                                      closer to 1 ==> more of dn()
c
c-----------------------------------------------------------------------

      integer i
      double precision vssd, vlen
      double precision dnrm2, ddot

      if(dnlen .le. delta) then

c        Newton step smaller than trust radius ==> take it

         call dcopy(n, dn, 1, d, 1)
         delta = dnlen
         dtype = 3

      elseif(ssdlen .ge. delta) then

c        take step in steepest descent direction

         call dcopy(n, ssd, 1, d, 1)
         call dscal(n, delta / ssdlen, d, 1)
         dtype = 1

      else

c        calculate convex combination of ssd and dn with length delta

         do i=1,n
            v(i) = dn(i) - ssd(i)
         enddo

         vssd = ddot(n,v,1,ssd,1)
         vlen = dnrm2(n,v,1)**2

         lambda =(-vssd+sqrt(vssd**2-vlen*(ssdlen**2-delta**2)))/vlen
         call dcopy(n, ssd, 1, d, 1)
         call daxpy(n, lambda, v, 1, d, 1)
         dtype = 2

      endif

      return
      end

      subroutine mm10_nwpure(n,xc,d,stepmx,scalex,solve_work,
     *                  xp,fp,fpnorm,xw,retcd,gcnt,priter,iter)
      use iso_Fortran_env
      use mm10_defs
      integer n,retcd,gcnt
      double precision  stepmx,fpnorm
      double precision  xc(*)
      double precision  d(*),xp(n),fp(n),xw(*)
      double precision  scalex(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work

      integer priter,iter

c-------------------------------------------------------------------------
c
c     Find a next iterate using geometric line search
c     along the newton direction
c
c     Arguments
c
c     In       n       Integer          dimension of problem
c     In       xc      Real(*)          current iterate
c     In       d       Real(*)          newton direction
c     In       stepmx  Real             maximum stepsize
c     In       scalex  Real(*)          diagonal scaling matrix for x()
c     In       fvec    Name             name of routine to calculate f()
c     In       xp      Real(*)          new x()
c     In       fp      Real(*)          new f(x)
c     In       fpnorm  Real             .5*||fp||**2
c     Out      xw      Real(*)          workspace for unscaling x
c
c     Out      retcd   Integer          return code
c                                         0 new satisfactory x() found (!always)
c
c     Out      gcnt    Integer          number of steps taken
c     In       priter  Integer           >0 if intermediate steps to be printed
c                                        -1 if no printing
c
c-------------------------------------------------------------------------

      integer i
      double precision  oarg(3)
      double precision  lambda
      double precision  dnrm2
      double precision  dlen

      integer idamax

      double precision Rone
      parameter(Rone=1.0d0)

c     safeguard initial step size

      dlen = dnrm2(n,d,1)
      if( dlen .gt. stepmx ) then
          lambda = stepmx / dlen
      else
          lambda = Rone
      endif

      retcd  = 0
      gcnt   = 1

c     compute the next iterate xp

      do i=1,n
         xp(i) = xc(i) + lambda*d(i)
      enddo

c     evaluate functions and the objective function at xp

      call mm10_nwfvec(xp,n,scalex,solve_work,fp,fpnorm,xw)

      if( priter .gt. 0) then
         oarg(1) = lambda
         oarg(2) = fpnorm
         oarg(3) = abs(fp(idamax(n,fp,1)))
         call mm10_nwprot(iter,1,oarg)
      endif

      return
      end

      subroutine mm10_nwqlsh(n,xc,fcnorm,d,g,stepmx,xtol,scalex,
     *           solve_work,xp,fp,fpnorm,xw,retcd,gcnt,priter,iter)
      use iso_Fortran_env
      use mm10_defs
      integer n,retcd,gcnt
      double precision  stepmx,xtol,fcnorm,fpnorm
      double precision  xc(*)
      double precision  d(*),g(*),xp(n),fp(n),xw(*)
      double precision  scalex(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work

      integer priter,iter

c-------------------------------------------------------------------------
c
c     Find a next acceptable iterate using a safeguarded quadratic line search
c     along the newton direction
c
c     Arguments
c
c     In       n       Integer          dimension of problem
c     In       xc      Real(*)          current iterate
c     In       fcnorm  Real             0.5 * || f(xc) ||**2
c     In       d       Real(*)          newton direction
c     In       g       Real(*)          gradient at current iterate
c     In       stepmx  Real             maximum stepsize
c     In       xtol    Real             relative step size at which
c                                       successive iterates are considered
c                                       close enough to terminate algorithm
c     In       scalex  Real(*)          diagonal scaling matrix for x()
c     In       fvec    Name             name of routine to calculate f()
c     In       xp      Real(*)          new x()
c     In       fp      Real(*)          new f(x)
c     In       fpnorm  Real             .5*||fp||**2
c     Out      xw      Real(*)           workspace for unscaling x(*)
c
c     Out      retcd   Integer          return code
c                                         0 new satisfactory x() found
c                                         1 no  satisfactory x() found
c                                           sufficiently distinct from xc()
c
c     Out      gcnt    Integer          number of steps taken
c     In       priter  Integer           >0 unit if intermediate steps to be printed
c                                        -1 if no printing
c
c-------------------------------------------------------------------------

      integer i, gcntmax
      double precision  alpha,slope,rsclen,oarg(4)
      double precision  lambda,lamhi,lamlo,t
      double precision  ddot,dnrm2, mm10_nudnrm, ftarg
      double precision  dlen, xpr(n)

      integer idamax

      parameter (alpha = 1.0d-4)

      double precision Rone, Rtwo, Rten
      parameter(Rone=1.0d0, Rtwo=2.0d0, Rten=10.0d0)
      gcntmax = 100

c     safeguard initial step size

      dlen = dnrm2(n,d,1)
      if( dlen .gt. stepmx ) then
          lamhi  = stepmx / dlen
      else
          lamhi  = Rone
      endif

c     compute slope  =  g-trans * d

      slope = ddot(n,g,1,d,1)

c     compute the smallest value allowable for the damping
c     parameter lambda ==> lamlo

      rsclen = mm10_nudnrm(n,d,xc)
      lamlo  = xtol / rsclen

c     initialization of retcd and lambda (linesearch length)

      retcd  = 2
      lambda = lamhi
      gcnt   = 0

      do while( (retcd .eq. 2) .and. (gcnt.le.gcntmax) )

c        compute next x

         do i=1,n
            xp(i) = xc(i) + lambda*d(i)
         enddo

c        evaluate functions and the objective function at xp

         call mm10_nwfvec(xp,n,scalex,solve_work,fp,fpnorm,xw)
         gcnt = gcnt + 1
         ftarg = fcnorm + alpha * lambda * slope

         if( priter .gt. 0) then
            oarg(1) = lambda
            oarg(2) = ftarg
            oarg(3) = fpnorm
            oarg(4) = abs(fp(idamax(n,fp,1)))
            call mm10_nwlsot(iter,1,oarg)
         endif

c        test whether the standard step produces enough decrease
c        of the objective function.
c        If not update lambda and compute a new next iterate

         if( fpnorm .le. ftarg ) then
            retcd = 0
         else
            t = ((-lambda**2)*slope/Rtwo)/(fpnorm-fcnorm-lambda*slope)
            lambda  = max(lambda / Rten , t)
            if(lambda .lt. lamlo) then
               retcd = 1
            endif
         endif

      enddo
      if(gcnt.ge.gcntmax) then
        retcd = 1
      endif

      return
      end

      subroutine mm10_nwtrup(n,fcnorm,g,sc,nwtstep,stepmx,xtol,delta,
     *                  fpred,retcd,xprev,fpnsav,fprev,xp,fp,
     *                  fpnorm)

      integer n,retcd
      double precision  fcnorm,stepmx,xtol,delta,fpred,fpnsav,fpnorm
      double precision  xp(n),g(*)
      double precision  sc(*),xprev(*),fprev(*),fp(n)
      logical nwtstep

c-------------------------------------------------------------------------
c
c     Decide whether to accept xp=xc+sc as the next iterate
c     and updates the trust region delta
c
c     Arguments
c
c     In       n       Integer         size of xc()
c     In       fcnorm  Real            .5*||f(xc)||**2
c     In       g       Real(*)         gradient at xc()
c     In       sc      Real(*)         current step
c     In       nwtstep Logical         true if sc is newton direction
c     In       stepmx  Real            maximum step size
c     In       xtol    Real            minimum step tolerance
c     Inout    delta   Real            trust region radius
c     In       fpred   Real            predicted value of .5*||f()||**2
c
c     Inout    retcd   Integer         return code
c                                       0 xp accepted as next iterate;
c                                         delta trust region for next iteration.
c
c                                       1 xp unsatisfactory but
c                                         accepted as next iterate because
c                                         xp-xc .lt. smallest allowable
c                                         step length.
c
c                                       2 f(xp) too large.
c                                         continue current iteration with
c                                         new reduced delta.
c
c                                       3 f(xp) sufficiently small, but
c                                         quadratic model predicts f(xp)
c                                         sufficiently well to continue current
c                                         iteration with new doubled delta.
c
c                                      On first entry, retcd must be 4
c
c     Wk       xprev   Real(*)         (internal) work
c     Wk       fpnsav  Real            (internal) work
c     Wk       fprev   Real(*)         (internal) work
c     Inout    xp      Real(*)         new iterate x()
c     Inout    fp      Real(*)         new f(xp)
c     Inout    fpnorm  Real            new .5*||f(xp)||**2
c
c-------------------------------------------------------------------------

      double precision  ared,pred,slope,sclen,rln,dltmp
      double precision  dnrm2,ddot,mm10_nudnrm
      logical ret3ok

      double precision Rone, Rtwo, Rthree, Rfour, Rten
      double precision Rhalf, Rpten
      parameter(Rpten = 0.1d0)
      parameter(Rhalf=0.5d0)
      parameter(Rone=1.0d0, Rtwo=2.0d0, Rthree=3.0d0, Rfour=4.0d0)
      parameter(Rten=10.0d0)

      double precision Rp99,Rp4, Rp75
      parameter(Rp99=Rone-Rten**(-2), Rp4=Rten**(-4), Rp75=Rthree/Rfour)

      double precision alpha
      parameter(alpha = Rp4)

c     pred measures the predicted reduction in the function value

      ared  = fpnorm - fcnorm
      pred  = fpred  - fcnorm
      slope = ddot(n,g,1,sc,1)

      if(retcd .ne. 3) then
         ret3ok = .false.
      else
         ret3ok = fpnorm .ge. fpnsav .or. ared .gt. alpha * slope
      endif

      if(retcd .eq. 3 .and. ret3ok) then

c        reset xp,fp,fpnorm to saved values and terminate global step

         retcd = 0
         call dcopy(n,xprev,1,xp,1)
         call dcopy(n,fprev,1,fp,1)
         fpnorm = fpnsav
c        reset delta to initial value
c        but beware
c           if the trial step was a Newton step then delta is reset to
c           .5 * length(Newton step) which will be smaller
c           because delta is set to length(Newton step) elsewhere
c           see nwddlg.f and nwpdlg.f
         delta  = Rhalf*delta

      elseif(ared .gt. alpha * slope) then

c        fpnorm too large (decrease not sufficient)

         rln = mm10_nudnrm(n,sc,xp)
         if(rln .lt. xtol) then

c           cannot find satisfactory xp sufficiently distinct from xc

            retcd = 1

         else

c           reduce trust region and continue current global step

            retcd = 2
            sclen = dnrm2(n,sc,1)
            dltmp = -slope*sclen/(Rtwo*(ared-slope))

            if(dltmp .lt. Rpten*delta) then
               delta = Rpten*delta
            else
               delta = min(Rhalf*delta, dltmp)
            endif

         endif

      elseif(retcd .ne. 2 .and. (abs(pred-ared) .le. Rpten*abs(ared))
     *      .and. (.not. nwtstep) .and. (delta .le. Rp99*stepmx)) then

c        pred predicts ared very well
c        attempt a doubling of the trust region and continue global step
c        when not taking a newton step and trust region not at maximum

         call dcopy(n,xp,1,xprev,1)
         call dcopy(n,fp,1,fprev,1)
         fpnsav = fpnorm
         delta  = min(Rtwo*delta,stepmx)
         retcd  = 3

      else

c        fpnorm sufficiently small to accept xp as next iterate.
c        Choose new trust region.

         retcd = 0
         if(ared .ge. Rpten*pred) then

c           Not good enough. Decrease trust region for next iteration

            delta = Rhalf*delta
         elseif( ared .le. Rp75*pred ) then

c           Wonderful. Increase trust region for next iteration

            delta = min(Rtwo*delta,stepmx)
         endif

      endif

      return
      end

      subroutine mm10_nwtcvg(xplus,fplus,xc,xtol,retcd,ftol,iter,
     *                  maxit,n,ierr,termcd)

      integer n,iter,maxit,ierr,termcd,retcd
      double precision  xtol,ftol
      double precision  xplus(*),fplus(*),xc(*)

c-------------------------------------------------------------------------
c
c     Decide whether to terminate the nonlinear algorithm
c
c     Arguments
c
c     In       xplus   Real(*)         new x values
c     In       fplus   Real(*)         new f values
c     In       xc      Real(*)         current x values
c     In       xtol    Real            stepsize tolerance
c     In       retcd   Integer         return code from global search routines
c     In       ftol    Real            function tolerance
c     In       iter    Integer         iteration number
c     In       maxit   Integer         maximum number of iterations allowed
c     In       n       Integer         size of x and f
c     In       ierr    Integer         return code of cndjac (condition estimation)
c
c     Out      termcd  Integer         termination code
c                                        0 no termination criterion satisfied
c                                          ==> continue iterating
c                                        1 norm of scaled function value is
c                                          less than ftol
c                                        2 scaled distance between last
c                                          two steps less than xtol
c                                        3 unsuccessful global strategy
c                                          ==> cannot find a better point
c                                        4 iteration limit exceeded
c                                        5 Jacobian too ill-conditioned
c                                        6 Jacobian singular
c
c-------------------------------------------------------------------------

      double precision  fmax,rsx, mm10_nuxnrm
      integer idamax

c     check whether function values are within tolerance

      termcd = 0

      if( ierr .ne. 0 ) then
         termcd = 4 + ierr
         return
      endif

      fmax = abs(fplus(idamax(n,fplus,1)))
      if( fmax .le. ftol) then
         termcd = 1
         return
      endif

c     initial check at start so there is no xplus
c     so only a check of function values is useful
      if(iter .eq. 0) return

      if(retcd .eq. 1) then
         termcd = 3
         return
      endif

c     check whether relative step length is within tolerance
c     Dennis Schnabel Algorithm A7.2.3

      rsx = mm10_nuxnrm(n, xplus, xc)
c      write(*,*) 'rsx', rsx, xplus(1), xc(1), 'xtol', xtol,n
      if(rsx .le. xtol) then
        termcd = 2
        return
      endif

c     check iteration limit
c       write(*,*) 'got here'
      if(iter .ge. maxit) then
         termcd = 4
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nweset(n,xc,fc,fcnorm,xp,fp,fpnorm,gcnt,priter,
     &      iter)
      integer n, gcnt, priter, iter
      double precision xc(*),fc(*),fcnorm,xp(n),fp(n),fpnorm

c-------------------------------------------------------------------------
c
c     calling routine got an error in decomposition/update of Jacobian/Broyden
c     jacobian an singular or too ill-conditioned
c     prepare return arguments
c
c     Arguments
c
c     In       n       Integer         size of x
c     In       xc      Real(*)         current (starting) x values
c     In       fc      Real(*)         function values f(xc)
c     In       fcnorm  Real            norm fc
c     Out      xp      Real(*)         final x values
c     Out      fp      Real(*)         function values f(xp)
c     Out      fpnorm  Real            final norm fp
c     Out      gcnt    Integer         # of backtracking steps (here set to 0)
c     In       priter  Integer         flag for type of output
c     In       iter    Integer         iteration counter
c
c-------------------------------------------------------------------------

      call dcopy(n,xc,1,xp,1)
      call dcopy(n,fc,1,fp,1)
      fpnorm = fcnorm
      gcnt   = 0
      if( priter .gt. 0 ) then
         call mm10_nwjerr(iter)
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_chkjac1(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,
     & solve_work,termcd)
      use iso_Fortran_env
      use mm10_defs
      integer lda,n,termcd
      double precision  A(lda,*),xc(*),fc(*)
      double precision  epsm,scalex(*)
      double precision  fz(*),wa(*),xw(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work, solve_work1

c-------------------------------------------------------------------------
c
c     Check the user supplied jacobian against its finite difference approximation
c
c     Arguments
c
c     In       A       Real(lda,*)     user supplied jacobian
c     In       lda     Integer         leading dimension of ajanal
c     In       xc      Real(*)         vector of x values
c     In       fc      Real(*)         function values f(xc)
c     In       n       Integer         size of x
c     In       epsm    Real            machine precision
c     In       scalex  Real(*)         scaling vector for x()
c     Wk       fz      Real(*)         workspace
c     Wk       wa      Real(*)         workspace
c     Wk       xw      Real(*)         workspace
c     In       fvec    Name            name of routine to evaluate f(x)
c     Out      termcd  Integer         return code
c                                        0  user supplied jacobian ok
c                                      -10  user supplied jacobian NOT ok
c
c-------------------------------------------------------------------------

      integer i,j,errcnt
      double precision  ndigit,p,h,xcj,dinf
      double precision  tol
      double precision  mm10_rnudif
      integer idamax

      integer MAXERR
      parameter(MAXERR=10)

      double precision Rquart, Rten, t1, t2
      parameter(Rquart=0.25d0, Rten=10.0d0)

c      write(*,*) 'got to start'
      termcd = 0

c     compute the finite difference jacobian and check it against
c     the analytic one

      ndigit = -log10(epsm)
      t1 = sqrt( epsm )
      t2 = sqrt( 1.0d0 / Rten**ndigit )
      p = max( t1, t2 )
      tol    = epsm**Rquart

      errcnt = 0
      call dcopy(n,xc,1,xw,1)
      call mm10_vunsc(n,xw,scalex)

c
c Copy the work so that it is not over-written
c      write(*,*) 'got to copy'
      call mm10_copy_work(solve_work,solve_work1)
c      write(*,*) 'got past copy'
c      write(*,*) 'solve_work%solvfnc chk1', solve_work%solvfnc
c      write(*,*) 'solve_work1%solvfnc', solve_work1%solvfnc
      do j=1,n
         h = p + p * abs(xw(j))
         xcj   = xw(j)
         xw(j) = xcj + h

c        avoid (small) rounding errors
c        h = xc(j) - xcj but not here to avoid clever optimizers

         h = mm10_rnudif(xw(j), xcj)

c      write(*,*) 'before call', j
         call mm10_fvec(solve_work1,xw,fz,n,j)
c      write(*,*) 'after call', j
         xw(j) = xcj

         do i=1,n
            wa(i) = (fz(i)-fc(i))/h
         enddo

         dinf = abs(wa(idamax(n,wa,1)))

         do i=1,n
            if(abs(A(i,j)-wa(i)).gt.tol*dinf) then
               errcnt = errcnt + 1
               if( errcnt .gt. MAXERR ) then
                  termcd = -10
                  return
               endif
               call mm10_nwckot(i,j,A(i,j),wa(i))
            endif
         enddo
      enddo
      call mm10_delete_state(solve_work1%props)

c      call vscal(n,xc,scalex)

      if( errcnt .gt. 0 ) then
         termcd = -10
      endif
      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_chkjac2(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,
     *                   solve_work,termcd,dsub,dsuper)
      use iso_Fortran_env
      use mm10_defs
      integer lda,n,termcd,dsub,dsuper
      double precision  A(lda,*),xc(*),fc(*)
      double precision  epsm,scalex(*)
      double precision  fz(*),wa(*),xw(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work, solve_work1

c-------------------------------------------------------------------------
c
c     Check the user supplied jacobian against its finite difference approximation
c
c     Arguments
c
c     In       A       Real(lda,*)     user supplied jacobian
c     In       lda     Integer         leading dimension of ajanal
c     In       xc      Real(*)         vector of x values
c     In       fc      Real(*)         function values f(xc)
c     In       n       Integer         size of x
c     In       epsm    Real            machine precision
c     In       scalex  Real(*)         scaling vector for x()
c     Wk       fz      Real(*)         workspace
c     Wk       wa      Real(*)         workspace
c     Wk       xw      Real(*)         workspace
c     In       fvec    Name            name of routine to evaluate f(x)
c     Out      termcd  Integer         return code
c                                        0  user supplied jacobian ok
c                                      -10  user supplied jacobian NOT ok
c
c-------------------------------------------------------------------------

      integer i,j,k,dsum,errcnt
      double precision  ndigit,p,h,dinf
      double precision  tol
      double precision w(n),xstep(n)

      integer MAXERR
      parameter(MAXERR=10)

      double precision Rquart, Rten, Rzero, t1, t2
      parameter(Rquart=0.25d0, Rten=10.0d0, Rzero=0.0d0)

      dsum = dsub + dsuper + 1

      termcd = 0

c     compute the finite difference jacobian and check it against
c     the user supplied one

      ndigit = -log10(epsm)
      t1 = sqrt( epsm )
      t2 = sqrt( 1.0d0 / Rten**ndigit )
      p = max( t1, t2 )
      tol    = epsm**Rquart

      errcnt = 0
      call dcopy(n,xc,1,xw,1)
      call mm10_vunsc(n,xw,scalex)

      do j=1,n
          xstep(j) = p + p * abs(xw(j))
          w(j) = xw(j)
      enddo

c
c Copy the work so that it is not over-written
      call mm10_copy_work(solve_work,solve_work1)
      do k=1,dsum
         do j=k,n,dsum
            xw(j) = xw(j) + xstep(j)
         enddo

c        for non finite values error message will be wrong
         call mm10_fvec(solve_work1,xw,fz,n,n+k)

         do j=k,n,dsum
             h = xstep(j)
             xw(j) = w(j)
             dinf = Rzero
             do i=max(j-dsuper,1),min(j+dsub,n)
                wa(i) = (fz(i)-fc(i)) / h
                if(abs(wa(i)).gt.dinf) dinf = abs(wa(i))
             enddo

             do i=max(j-dsuper,1),min(j+dsub,n)
                if(abs(A(i,j)-wa(i)).gt.tol*dinf) then
                   errcnt = errcnt + 1
                   if( errcnt .gt. MAXERR ) then
                      termcd = -10
                      return
                   endif
                   call mm10_nwckot(i,j,A(i,j),wa(i))
                endif
             enddo
         enddo
      enddo
      call mm10_delete_state(solve_work1%props)

c      call vscal(n,xc,scalex)

      if( errcnt .gt. 0 ) then
         termcd = -10
      endif
      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_chkjaci(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,
     & solve_work,termcd)
      use iso_Fortran_env
      use mm10_defs
      integer lda,n,termcd
      double precision  A(lda,*),xc(*),fc(*)
      double precision  epsm,scalex(*)
      double precision  fz(*),wa(*),xw(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work, solve_work1

c-------------------------------------------------------------------------
c
c     Check the user supplied jacobian against its complex difference approximation
c
c     Arguments
c
c     In       A       Real(lda,*)     user supplied jacobian
c     In       lda     Integer         leading dimension of ajanal
c     In       xc      Real(*)         vector of x values
c     In       fc      Real(*)         function values f(xc)
c     In       n       Integer         size of x
c     In       epsm    Real            machine precision
c     In       scalex  Real(*)         scaling vector for x()
c     Wk       fz      Real(*)         workspace
c     Wk       wa      Real(*)         workspace
c     Wk       xw      Real(*)         workspace
c     In       fvec    Name            name of routine to evaluate f(x)
c     Out      termcd  Integer         return code
c                                        0  user supplied jacobian ok
c                                      -10  user supplied jacobian NOT ok
c
c-------------------------------------------------------------------------

      integer i,j,errcnt
      double precision  ndigit,p,h,xcj,dinf
      double precision, dimension(n) :: zeroN
      complex(kind=real64), dimension(n) :: xwc, fzc
      double precision  tol
      complex(kind=real64)  mm10_rnudifi
      integer idamax

      integer MAXERR
      parameter(MAXERR=10)

      double precision Rquart, Rten
      parameter(Rquart=0.25d0, Rten=10.0d0)

      complex(kind=real64) :: i1
      i1 = (0.d0, 1.d0)

      zeroN = 0.d0

c      write(*,*) 'got to start'
      termcd = 0

c     compute the complex difference jacobian and check it against
c     the analytic one


      ndigit = -log10(epsm)
      t1 = sqrt( epsm )
      t2 = sqrt( 1.0d0 / Rten**ndigit )
      p = max( t1, t2 )
c      p = (Rten**4)*(max(Rten**(-ndigit),epsm)) ! this one did not work, gave small errors
      tol    = epsm**Rquart
      
      errcnt = 0
c      call dcopy(n,xc,1,xw,1)
      xwc(1:n) = dcmplx (xc(1:n), zeroN)
      call mm10_vunsci(n,xwc,scalex)

c
c Copy the work so that it is not over-written
c      write(*,*) 'got to copy'
      call mm10_copy_work(solve_work,solve_work1)
c      write(*,*) 'got past copy'
c      write(*,*) 'solve_work%solvfnc chk1', solve_work%solvfnc
c      write(*,*) 'solve_work1%solvfnc', solve_work1%solvfnc
      do j=1,n
         h = p + p * abs(real(xwc(j)))
         xcj   = real(xwc(j))
         xwc(j) = cmplx(xcj, h)

c        avoid (small) rounding errors
c        h = xc(j) - xcj but not here to avoid clever optimizers

c         h = mm10_rnudifi(xwc(j), xcj)

c      write(*,*) 'before call', j
         call mm10_fveci(solve_work1,xwc,fzc,n,j)
c      write(*,*) 'after call', j
         xwc(j) = cmplx(xcj, 0.d0)

         do i=1,n
            wa(i) = aimag(fzc(i))/h
         enddo

         dinf = abs(wa(idamax(n,wa,1)))

         do i=1,n
            if(abs(A(i,j)-wa(i)).gt.tol*dinf) then
               errcnt = errcnt + 1
               if( errcnt .gt. MAXERR ) then
                  termcd = -10
                  return
               endif
               call mm10_nwckot(i,j,A(i,j),wa(i))
            endif
         enddo
      enddo
      call mm10_delete_state(solve_work1%props)

c      call vscal(n,xc,scalex)

      if( errcnt .gt. 0 ) then
         termcd = -10
      endif
      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_chkjac(A,lda,xc,fc,n,epsm,jacflg,scalex,fz,wa,xw,
     *                  solve_work,termcd)
      use iso_Fortran_env
      use mm10_defs
      integer lda,n,termcd,jacflg(5)
      double precision  A(lda,*),xc(*),fc(*)
      double precision  epsm,scalex(*)
      double precision  fz(*),wa(*),xw(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work

c-------------------------------------------------------------------------
c
c     Check the user supplied jacobian against its finite difference approximation
c
c     Arguments
c
c     In       A       Real(lda,*)     user supplied jacobian
c     In       lda     Integer         leading dimension of ajanal
c     In       xc      Real(*)         vector of x values
c     In       fc      Real(*)         function values f(xc)
c     In       n       Integer         size of x
c     In       epsm    Real            machine precision
c     In       jacflg  Integer(*)      indicates how to compute jacobian
c                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
c                                                  3: user supplied banded
c                                      jacflg[2]: number of sub diagonals or -1 if not banded
c                                      jacflg[3]: number of super diagonals or -1 if not banded
c                                      jacflg[4]: 1 if adjusting jacobian allowed when
c                                                   singular or illconditioned
c     In       scalex  Real(*)         scaling vector for x()
c     Wk       fz      Real(*)         workspace
c     Wk       wa      Real(*)         workspace
c     Wk       xw      Real(*)         workspace
c     In       fvec    Name            name of routine to evaluate f(x)
c     Out      termcd  Integer         return code
c                                        0  user supplied jacobian ok
c                                      -10  user supplied jacobian NOT ok
c
c-------------------------------------------------------------------------

c      write(*,*) 'solve_work%solvfnc chkjac', solve_work%solvfnc
      if(jacflg(1) .eq. 3) then
c        user supplied and banded
         call mm10_chkjac2(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,
     *                     solve_work,termcd,jacflg(2),jacflg(3))
      elseif( solve_work%props%tang_calc .eq. 4 ) then
c        check using complex derivative
         call mm10_chkjaci(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,
     &                     solve_work,termcd)
      else
         call mm10_chkjac1(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,
     &                     solve_work,termcd)
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_fdjac0(xc,fc,n,epsm,solve_work,fz,rjac,ldr)

      use iso_Fortran_env
      use mm10_defs
      integer ldr,n
      double precision  epsm
      double precision  rjac(ldr,*),fz(*),xc(*),fc(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work, solve_work1

c-------------------------------------------------------------------------
c
c     Compute the finite difference jacobian at the current point xc
c
c     Arguments
c
c     In       xc      Real(*)         current point
c     In       fc      Real(*)         function values at current point
c     In       n       Integer         size of x and f
c     In       epsm    Real            machine precision
c     In       fvec    Name            name of routine to evaluate f(x)
c     Wk       fz      Real(*)         workspace
c     Out      rjac    Real(ldr,*)     jacobian matrix at x
c                                        entry [i,j] is derivative of
c                                        f(i) wrt to x(j)
c     In       ldr     Integer         leading dimension of rjac
c
c-------------------------------------------------------------------------

      integer i,j
      double precision  ndigit,p,h,xcj
      double precision  mm10_rnudif

      double precision Rten, t1, t2
      parameter(Rten=10.d0)

      ndigit = -log10(epsm)
      t1 = sqrt( epsm )
      t2 = sqrt( 1.0d0 / Rten**ndigit )
      p = max( t1, t2 )

c
c Copy the work so that it is not over-written
      call mm10_copy_work(solve_work,solve_work1)
      do j=1,n
         h = p + p * abs(xc(j))

c        or as alternative h  = p * max(Rone, abs(xc(j)))

         xcj   = xc(j)
         xc(j) = xcj + h

c        avoid (small) rounding errors
c        h = xc(j) - xcj  but not here to avoid clever optimizers

         h = mm10_rnudif(xc(j), xcj)
         call mm10_fvec(solve_work1,xc,fz,n,j)
         xc(j) = xcj
         do i=1,n
            rjac(i,j) = (fz(i)-fc(i)) / h
         enddo
      enddo
      call mm10_delete_state(solve_work1%props)

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_fdjac2(xc,fc,n,epsm,solve_work,fz,rjac,ldr,dsub,
     *                  dsuper,w,xstep)
      use iso_Fortran_env
      use mm10_defs
      integer ldr,n,dsub,dsuper
      double precision  epsm
      double precision  rjac(ldr,*),fz(*),xc(*),fc(*)
      double precision  w(*), xstep(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work, solve_work1

c-------------------------------------------------------------------------
c
c     Compute a banded finite difference jacobian at the current point xc
c
c     Arguments
c
c     In       xc      Real(*)         current point
c     In       fc      Real(*)         function values at current point
c     In       n       Integer         size of x and f
c     In       epsm    Real            machine precision
c     In       fvec    Name            name of routine to evaluate f(x)
c     Wk       fz      Real(*)         workspace
c     Out      rjac    Real(ldr,*)     jacobian matrix at x
c                                        entry [i,j] is derivative of
c                                        f(i) wrt to x(j)
c     In       ldr     Integer         leading dimension of rjac
c     In       dsub    Integer         number of subdiagonals
c     In       dsuper  Integer         number of superdiagonals
c     Internal w       Real(*)         for temporary saving of xc
c     Internal xstep   Real(*)         stepsizes
c
c-------------------------------------------------------------------------

      integer i,j,k
      double precision  ndigit,p,h
      double precision  mm10_rnudif

      double precision Rten, t1, t2
      parameter(Rten=10.d0)

      integer dsum

      dsum = dsub + dsuper + 1

      ndigit = -log10(epsm)
      t1 = sqrt( epsm )
      t2 = sqrt( 1.0d0 / Rten**ndigit )
      p = max( t1, t2 )

      do k=1,n
         xstep(k) = p + p * abs(xc(k))
      enddo

c
c Copy the work so that it is not over-written
      call mm10_copy_work(solve_work,solve_work1)
c
      do k=1,dsum
         do j=k,n,dsum
            w(j) = xc(j)
            xc(j) = xc(j) + xstep(j)
         enddo

         call mm10_fvec(solve_work1,xc,fz,n,n+k)
         do j=k,n,dsum
             call mm10_nuzero(n,rjac(1,j))
c            fdjac0 for why
c            doing this ensures that results for fdjac2 and fdjac0 will be identical
             h = mm10_rnudif(xc(j),w(j))
             xc(j) = w(j)
             do i=max(j-dsuper,1),min(j+dsub,n)
                rjac(i,j) = (fz(i)-fc(i)) / h
             enddo
         enddo
      enddo
      call mm10_delete_state(solve_work1%props)

      return
      end

c-----------------------------------------------------------------------

      function mm10_nudnrm(n, d, x)
      integer n
      double precision  d(*), x(*)
      double precision mm10_nudnrm

c-------------------------------------------------------------------------
c
c     calculate  max( abs(d[*]) / max(x[*],1) )
c
c     Arguments
c
c     In   n        Integer       number of elements in d() and x()
c     In   d        Real(*)       vector d
c     In   x        Real(*)       vector x
c
c-------------------------------------------------------------------------

      integer i
      double precision  t

      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)

      t = Rzero
      do i=1,n
         t = max(t, abs(d(i)) / max(abs(x(i)),Rone))
      enddo
      mm10_nudnrm = t

      return
      end

c-----------------------------------------------------------------------

      function mm10_nuxnrm(n, xn, xc)
      integer n
      double precision  xn(*), xc(*)
      double precision mm10_nuxnrm

c-------------------------------------------------------------------------
c
c     calculate  max( abs(xn[*]-xc[*]) / max(xn[*],1) )
c
c     Arguments
c
c     In   n        Integer       number of elements in xn() and xc()
c     In   xn       Real(*)       vector xn
c     In   xc       Real(*)       vector xc
c
c-------------------------------------------------------------------------

      integer i
      double precision  t

      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)

      t = Rzero
      do i=1,n
         t = max(t, abs(xn(i)-xc(i)) / max(abs(xn(i)),Rone))
      enddo
      mm10_nuxnrm = t

      return
      end

c-----------------------------------------------------------------------

      function mm10_rnudif(x, y)
      double precision x, y
      double precision mm10_rnudif

c-------------------------------------------------------------------------
c
c     Return difference of x and y (x - y)
c
c     Arguments
c
c     In   x  Real      argument 1
c     In   y  Real      argument 2
c
c-------------------------------------------------------------------------

      mm10_rnudif = x - y
      return
      end

c-----------------------------------------------------------------------

      function mm10_rnudifi(x, y)
      use iso_Fortran_env
      complex(kind=real64) x, y
      complex(kind=real64) mm10_rnudifi

c-------------------------------------------------------------------------
c
c     Return difference of x and y (x - y)
c
c     Arguments
c
c     In   x  Comp      argument 1
c     In   y  Comp      argument 2
c
c-------------------------------------------------------------------------

      mm10_rnudifi = x - y
      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_compmu(r,ldr,n,mu,y,ierr)

      integer ldr,n,ierr
      double precision r(ldr,*),mu,y(*)

c-------------------------------------------------------------------------
c
c     Compute a small perturbation mu for the (almost) singular matrix R.
c     mu is used in the computation of the Levenberg-Marquardt step.
c
c     Arguments
c
c     In       R       Real(ldr,*)     upper triangular matrix from QR
c     In       ldr     Integer         leading dimension of R
c     In       n       Integer         column dimension of R
c     Out      mu      Real            sqrt(l1 norm of R * infinity norm of R
c                                      * n * epsm * 100) designed to make
c                                        trans(R)*R + mu * I not singular
c     Wk       y       Real(*)         workspace for dlange
c     Out      ierr    Integer         0 indicating mu ok
c                                      3 indicating mu much too small
c
c-------------------------------------------------------------------------

      double precision aifnrm,al1nrm,epsm
      double precision dlantr
      double precision mm10_epsmch

      double precision Rhund
      parameter(Rhund=100d0)

c     get the infinity norm of R
c     get the l1 norm of R
      ierr = 0
      aifnrm = dlantr('I','U','N',n,n,r,ldr,y)
      al1nrm = dlantr('1','U','N',n,n,r,ldr,y)
      epsm = mm10_epsmch()
      mu = sqrt(n*epsm*Rhund)*aifnrm*al1nrm
c     matrix consists of zero's or near zero's
c     LM correction in liqrev will not work
      if( mu .le. Rhund*epsm ) then
         ierr = 3
      endif
      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_cndjac(n,r,ldr,cndtol,rcond,rcdwrk,icdwrk,ierr)
      integer n,ldr,icdwrk(*),ierr
      double precision cndtol,rcond,r(ldr,*),rcdwrk(*)

c---------------------------------------------------------------------
c
c     Check r for singularity and/or ill conditioning
c
c     Arguments
c
c     In       n       Integer         dimension of problem
c     In       r       Real(ldr,*)     upper triangular R from QR decomposition
c     In       ldr     Integer         leading dimension of rjac
c     In       cndtol  Real            tolerance of test for ill conditioning
c                                       when rcond <= cndtol then ierr is set to 1
c                                       cndtol should be >= machine precision
c     Out      rcond   Real            inverse condition  of r
c     Wk       rcdwrk  Real(*)         workspace (for dtrcon)
c     Wk       icdwrk  Integer(*)      workspace (fordtrcon)
c     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
c                                      1 indicating Jacobian too ill-conditioned
c                                      2 indicating Jacobian completely singular
c
c---------------------------------------------------------------------

      integer i,info
      logical rsing
      double precision Rzero
      parameter(Rzero=0.0d0)

      ierr = 0

      rsing = .false.
      do i=1,n
         if( r(i,i) .eq. Rzero ) then
             rsing = .true.
         endif
      enddo

      if( rsing ) then
         ierr = 2
         rcond = Rzero
      else
         call dtrcon('1','U','N',n,r,ldr,rcond,rcdwrk,icdwrk,info)
         if( rcond .eq. Rzero ) then
             ierr = 2
         elseif( rcond .le. cndtol ) then
             ierr = 1
         endif
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwfjac(x,scalex,f,fq,n,epsm,jacflg,solve_work,
     *                  rjac,ldr,xw,w,xstep,priter)
      use iso_Fortran_env
      use mm10_defs
      integer ldr,n,jacflg(5), priter
      double precision  epsm
      double precision  x(*),f(*),scalex(*),xw(*),w(*),xstep(*)
      double precision  rjac(ldr,*),fq(*)
      include 'include_mm10'
      type(mm10_working_data) :: solve_work

c-------------------------------------------------------------------------
c
c     Calculate the jacobian  matrix
c
c     Arguments
c
c     In       x       Real(*)         (scaled) current x values
c     In       scalex  Real(*)         scaling factors x
c     In       f       Real(*)         function values f(x)
c     Wk       fq      Real(*)         (internal) workspace
c     In       n       Integer         size of x and f
c     In       epsm    Real            machine precision
c     In       jacflg  Integer(*)      indicates how to compute jacobian
c                                      jacflg[1]:  0 numeric; 1 user supplied; 2 numerical banded
c                                                  3: user supplied banded
c                                      jacflg[2]: number of sub diagonals or -1 if not banded
c                                      jacflg[3]: number of super diagonals or -1 if not banded
c                                      jacflg[4]: 1 if adjusting jacobian allowed when
c                                                   singular or illconditioned
c     In       fvec    Name            name of routine to evaluate f()
c     In       mkjac   Name            name of routine to evaluate jacobian
c     Out      rjac    Real(ldr,*)     jacobian matrix (unscaled)
c     In       ldr     Integer         leading dimension of rjac
c     Internal xw      Real(*)         used for storing unscaled x
c     Internal w       Real(*)         workspace for banded jacobian
c     Internal xstep   Real(*)         workspace for banded jacobian
c
c-------------------------------------------------------------------------

c     compute the finite difference or analytic jacobian at x

      call dcopy(n,x,1,xw,1)
      ! write(*,*) 'in J'
      call mm10_vunsc(n,xw,scalex)

c      if( priter .gt. 0 ) then

c            write(*,*) 'xw', xw(1:n)

c      endif
      if(jacflg(1) .eq. 0) then
         call mm10_fdjac0(xw,f,n,epsm,solve_work,fq,rjac,ldr)
      elseif(jacflg(1) .eq. 2) then
         call mm10_fdjac2(xw,f,n,epsm,solve_work,fq,rjac,ldr,jacflg(2),
     *               jacflg(3),w,xstep)
      else
         call mm10_fjac(solve_work,rjac,ldr,xw,n)
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwscjac(n,rjac,ldr,scalex)
      integer n, ldr
      double precision rjac(ldr,*), scalex(*)

c-------------------------------------------------------------------------
c
c     Scale jacobian
c
c     Arguments
c
c     In       n       Integer         size of x and f
c     Inout    rjac    Real(ldr,*)     jacobian matrix
c     In       ldr     Integer         leading dimension of rjac
c     In       scalex  Real(*)         scaling factors for x
c
c-------------------------------------------------------------------------

      integer j
      double precision t, Rone
      parameter(Rone=1.0d0)

      do j = 1,n
         t = Rone/scalex(j)
         call dscal(n,t,rjac(1,j),1)
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwunscjac(n,rjac,ldr,scalex)
      integer n, ldr
      double precision rjac(ldr,*), scalex(*)

c-------------------------------------------------------------------------
c
c     Unscale jacobian
c
c     Arguments
c
c     In       n       Integer         size of x and f
c     Inout    rjac    Real(ldr,*)     jacobian matrix
c     In       ldr     Integer         leading dimension of rjac
c     In       scalex  Real(*)         scaling factors for x
c
c-------------------------------------------------------------------------

      integer j
      double precision t

      do j = 1,n
         t = scalex(j)
         call dscal(n,t,rjac(1,j),1)
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwcpsx(n,rjac,ldr,scalex,epsm, mode)

      integer ldr,n,mode
      double precision  epsm
      double precision  scalex(*)
      double precision  rjac(ldr,*)

c-------------------------------------------------------------------------
c
c     Calculate scaling factors from the jacobian  matrix
c
c     Arguments
c
c     In       n       Integer         size of x and f
c     In       rjac    Real(ldr,*)     jacobian matrix
c     In       ldr     Integer         leading dimension of rjac
c     Out      scalex  Real(*)         scaling factors for x
c     In       epsm    Real            machine precision
c     In       mode    Integer         1: initialise, >1: adjust
c-------------------------------------------------------------------------

      integer k
      double precision  dnrm2

      if( mode .eq. 1 ) then
         do k=1,n
            scalex(k) = dnrm2(n,rjac(1,k),1)
            if( scalex(k) .le. epsm ) scalex(k) = 1
         enddo
      else if( mode .gt. 1 ) then
         do k=1,n
            scalex(k) = max(scalex(k),dnrm2(n,rjac(1,k),1))
         enddo
      endif
      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwcpmt(n, x, scalex, factor, wrk, stepsiz)
      double precision x(*), scalex(*), wrk(*)
      double precision factor, stepsiz
      integer n

c-------------------------------------------------------------------------
c
c     Calculate maximum stepsize
c
c     Arguments
c
c     In       n       Integer     size of x
c     In       x       Real(*)     x-values
c     In       scalex  Real(*)     scaling factors for x
c     In       factor  Real        multiplier
c     Inout    wrk     Real(*)     workspace
c     Out      stepsiz Real        stepsize
c
c     Currently not used
c     Minpack uses this to calculate initial trust region size
c     Not (yet) used in this code because it doesn't seem to help
c     Manually setting an initial trust region size works better
c
c-------------------------------------------------------------------------

      double precision Rzero
      parameter(Rzero=0.0d0)

      double precision  dnrm2

      call dcopy(n,x,1,wrk,1)
      call mm10_vscal(n,wrk,scalex)
      stepsiz = factor * dnrm2(n,wrk,1)
      if( stepsiz .eq. Rzero ) stepsiz = factor
      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_vscal(n,x,sx)

      integer n
      double precision  x(*),sx(*)

c-------------------------------------------------------------------------
c
c     Scale a vector x
c
c     Arguments
c
c     In       n       Integer         size of x
c     Inout    x       Real(*)         vector to scale
c     In       sx      Real(*)         scaling vector
c
c-------------------------------------------------------------------------

      integer i

      do i = 1,n
         x(i) = sx(i) * x(i)
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_vunsc(n,x,sx)

      integer n
      double precision  x(*),sx(*)

c-------------------------------------------------------------------------
c
c     Unscale a vector x
c
c     Arguments
c
c     In       n       Integer         size of x
c     Inout    x       Real(*)         vector to unscale
c     In       sx      Real(*)         scaling vector
c
c-------------------------------------------------------------------------

      integer i

      do i = 1,n
         x(i) = x(i) / sx(i)
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_vunsci(n,x,sx)
      use iso_Fortran_env
      integer n
      complex(kind=real64)  x(*)
      double precision  sx(*)

c-------------------------------------------------------------------------
c
c     Unscale a vector x
c
c     Arguments
c
c     In       n       Integer         size of x
c     Inout    x       Comp(*)         vector to unscale
c     In       sx      Real(*)         scaling vector
c
c-------------------------------------------------------------------------

      integer i

      do i = 1,n
         x(i) = x(i) / sx(i)
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwfvec(x,n,scalex,solve_work,f,fnorm,xw)
      use iso_Fortran_env
      use mm10_defs
      integer n
      double precision  x(*),xw(*),scalex(*),f(*),fnorm
      include 'include_mm10'
      type(mm10_working_data) :: solve_work

c-------------------------------------------------------------------------
c
c     Evaluate the function at current iterate x and scale its value
c
c     Arguments
c
c     In       x       Real(*)         x
c     In       n       Integer         size of x
c     In       scalex  Real(*)         scaling vector for x
c     In       fvec    Name            name of routine to calculate f(x)
c     Out      f       Real(*)         f(x)
c     Out      fnorm   Real            .5*||f(x)||**2
c     Internal xw      Real(*)         used for storing unscaled xc
c
c-------------------------------------------------------------------------

      double precision dnrm2

      double precision Rhalf
      parameter(Rhalf=0.5d0)

      call dcopy(n,x,1,xw,1)
      call mm10_vunsc(n,xw,scalex)
      call mm10_fvec(solve_work,xw,f,n,0)

      fnorm = Rhalf * dnrm2(n,f,1)**2

      return
      end

c-----------------------------------------------------------------------

      function mm10_epsmch()

c     Return machine precision
c     Use Lapack routine

      double precision mm10_epsmch
      double precision dlamch
      external dlamch

c     dlamch('e') returns negeps (1-eps)
c     dlamch('p') returns 1+eps

      mm10_epsmch = dlamch('p')

      return
      end

c-----------------------------------------------------------------------

      function mm10_dblhuge()

c     Return largest double precision number
c     Use Lapack routine

      double precision mm10_dblhuge
      double precision dlamch
      external dlamch

c     dlamch('o') returns max double precision

      mm10_dblhuge = dlamch('o')

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwckot(i,j,aij,wi)
      integer i,j
      double precision aij,wi

      write(*,*) 'Chkjac  possible error in jacobian[i,j] = ',
     &   i, j,  aij
      write(*,*) 'Estimated[i,j] = ', i, j, wi

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwlsot(iter,lstep,oarg)
      integer iter,lstep
      double precision oarg(4)
c{
c    /*
c     * Linesearch output
c     */

      if(lstep.le.0) then
        if(lstep.eq.-1) then
        write(*,*) '  Iter Jac Lambda Ftarg Fnorm Largest |f|'
        endif
        write(*,*) iter, '                ', oarg(1), oarg(2)
      else 
        write(*,'(" ",i3,"  ",E10.3,"  ",E10.3,
     &   "  ",E10.3,"  ",E10.3)')
     &    iter, oarg(1), oarg(2), oarg(3), oarg(4)
      end if

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwdgot(iter,lstep,retcd,oarg)
      integer iter,lstep,retcd
      double precision oarg(6)
c    /*
c     * Double dogleg output
c     */

c    char step;

c    /*
c     *  C gradient (cauchy) step
c     *  N newton step
c     *  P partial newton step
c     *  W convex combination of P and C
c     */

      if(lstep.le.0) then
        if(lstep.eq.-1) then
       write(*,*) '  Iter Jac Lambda Eta Dlt0 Dltn Fnorm Largest |f|'
        endif
        write(*,*) iter, '                ', oarg(1), oarg(2)
      else 
        if( lstep.eq.2 ) then
        write(*,'(" ",i3," CWPN",i1,"  ",d11.3)')
     &    iter, lstep, oarg(1)
        else
        write(*,'(" ",i3," CWPN",i1,"  ",d11.3,"  ",d11.3
     &,"  ",d11.3,"  ",d11.3,"  ",d11.3)')
     &    iter, lstep, oarg(4), oarg(2), oarg(3), oarg(5), oarg(6)
        endif
      end if

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwprot(iter,lstep,oarg)
      integer iter,lstep
      double precision oarg(3)
c    /*
c     * None global method output
c     */

      if(lstep.le.0) then
        if(lstep.eq.-1) then
        write(*,*) '  Iter Jac Lambda Fnorm Largest |f|'
        endif
        write(*,*) iter, '                ', oarg(1), oarg(2)
      else 
        write(*,'(" ",i3,"  ",E10.3,"  ",E10.3,
     &   "  ",E10.3)')
     &    iter, oarg(1), oarg(2), oarg(3)
      end if

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwpwot(iter,lstep,retcd,oarg)
      integer iter,lstep,retcd
      double precision oarg(5)
c    /*
c     * Single dogleg output
c     */

c    char step;

c    /*
c     *  C gradient (cauchy) step
c     *  N newton step
c     *  W convex combination of P and C
c     */

      if(lstep.le.0) then
        if(lstep.eq.-1) then
        write(*,*) '  Iter Jac Lambda Dlt0 Dltn Fnorm Largest |f|'
        endif
        write(*,*) iter, '                ', oarg(1), oarg(2)
      else 
        if( lstep.eq.2 ) then
        write(*,'(" ",i3," CWN",i1,"  ",d11.3)')
     &    iter, lstep, oarg(1)
        else
        write(*,'(" ",i3," CWN",i1,"  ",d11.3,"  ",d11.3
     &,"  ",d11.3,"  ",d11.3)')
     &    iter, lstep, oarg(2), oarg(3), oarg(4), oarg(5)
        endif
      end if

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwmhot(iter,lstep,retcd,oarg)
      integer iter,lstep,retcd
      double precision oarg(6)
c    /*
c     * More-Hebden-Levenberg-Marquardt output
c     */

c    char step;

c    /*
c     *  H MHLM (hook)step
c     *  N newton step
c     */

      if(lstep.le.0) then
        if(lstep.eq.-1) then
       write(*,*) '  Iter Jac mu dnorm Dlt0 Dltn Fnorm Largest |f|'
        endif
        write(*,*) iter, '                ', oarg(1), oarg(2)
      else 
        if( lstep.eq.1 ) then
        write(*,'(" ",i3," HN",i1,"  ",d11.3)')
     &    iter, lstep, oarg(1)
        else
        write(*,'(" ",i3," HN",i1,"  ",d11.3,"  ",d11.3
     &,"  ",d11.3,"  ",d11.3,"  ",d11.3)')
     &    iter, lstep, oarg(4), oarg(2), oarg(3), oarg(5), oarg(6)
        endif
      end if

      return
      end


c-----------------------------------------------------------------------

      subroutine mm10_nwjerr(iter)
      integer iter

      write(*,*) 'nwjerr'

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_nwrowhdr(iter)
      integer iter

      write(*,*) 'nwrowhdr'

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_copy_work(solve_work,solve_work1)
      use iso_Fortran_env
      use mm10_defs
      implicit none
      include 'include_mm10'
      type(mm10_working_data) :: solve_work, solve_work1
      
c      write(*,*) 'I found out we did not need this'
c      solve_work1 = solve_work
c      return

c      write(*,*) 'I found out we did not need this'
      call mm10_copy_props(solve_work1%props,solve_work%props)
      call mm10_copy_state(solve_work1%np1,solve_work%np1)
      call mm10_copy_state(solve_work1%np0,solve_work%np0)
      solve_work1%vec1 = solve_work%vec1
      solve_work1%vec2 = solve_work%vec2
      solve_work1%arr1 = solve_work%arr1
      solve_work1%arr2 = solve_work%arr2
      solve_work1%ivec1 = solve_work%ivec1
      solve_work1%ivec2 = solve_work%ivec2
      solve_work1%gaspt = solve_work%gaspt
      solve_work1%solvfnc = solve_work%solvfnc
      solve_work1%x2 = solve_work%x2
      solve_work1%locdebug = solve_work%locdebug


      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_copy_props( props1, props )

      use mm10_defs
      implicit none

      integer :: sizeArr, allocate_status, dim1, dim2
      type(crystal_props) :: props, props1

      props1%rate_n = props%rate_n
      props1%tau_hat_y = props%tau_hat_y
      props1%G_0_y = props%G_0_y
      props1%burgers = props%burgers
      props1%p_v = props%p_v
      props1%q_v = props%q_v
      props1%boltzman = props%boltzman
      props1%theta_0 = props%theta_0
      props1%eps_dot_0_v = props%eps_dot_0_v
      props1%eps_dot_0_y = props%eps_dot_0_y
      props1%p_y = props%p_y
      props1%q_y = props%q_y
      props1%tau_a = props%tau_a
      props1%tau_hat_v = props%tau_hat_v
      props1%G_0_v = props%G_0_v
      props1%k_0 = props%k_0
      props1%mu_0 = props%mu_0
      props1%D_0 = props%D_0
      props1%T_0 = props%T_0
      props1%tau_y = props%tau_y
      props1%tau_v = props%tau_v
      props1%voche_m = props%voche_m
      props1%u1 = props%u1
      props1%u2 = props%u2
      props1%u3 = props%u3
      props1%u4 = props%u4
      props1%u5 = props%u5
      props1%u6 = props%u6
      props1%atol = props%atol
      props1%atol1 = props%atol1
      props1%rtol = props%rtol
      props1%rtol1 = props%rtol1
      props1%xtol = props%xtol
      props1%xtol1 = props%xtol1
      props1%g = props%g
      props1%ms = props%ms
      props1%qs = props%qs
      props1%ns = props%ns
      props1%stiffness = props%stiffness
      props1%angle_type = props%angle_type
      props1%angle_convention = props%angle_convention
      props1%nslip = props%nslip
      props1%h_type = props%h_type
      props1%miter = props%miter
      props1%gpp = props%gpp
      props1%s_type = props%s_type
      props1%cnum = props%cnum
      props1%method = props%method
      props1%st_it = props%st_it
      props1%num_hard = props%num_hard
      props1%out = props%out
      props1%tang_calc = props%tang_calc
      props1%solver = props%solver
      props1%strategy = props%strategy
      props1%debug = props%debug
      props1%gpall = props%gpall
c
      props1%cp_001 = props%cp_001
      props1%cp_002 = props%cp_002
      props1%cp_003 = props%cp_003
      props1%cp_004 = props%cp_004
      props1%cp_005 = props%cp_005
      props1%cp_006 = props%cp_006
      props1%cp_007 = props%cp_007
      props1%cp_008 = props%cp_008
      props1%cp_009 = props%cp_009
      props1%cp_010 = props%cp_010
      props1%cp_011 = props%cp_011
      props1%cp_012 = props%cp_012
      props1%cp_013 = props%cp_013
      props1%cp_014 = props%cp_014
      props1%cp_015 = props%cp_015
      props1%cp_016 = props%cp_016
      props1%cp_017 = props%cp_017
      props1%cp_018 = props%cp_018
      props1%cp_019 = props%cp_019
      props1%cp_020 = props%cp_020
      props1%cp_021 = props%cp_021
      props1%cp_022 = props%cp_022
      props1%cp_023 = props%cp_023
      props1%cp_024 = props%cp_024
      props1%cp_025 = props%cp_025
      props1%cp_026 = props%cp_026
      props1%cp_027 = props%cp_027
      props1%cp_028 = props%cp_028
      props1%cp_029 = props%cp_029
      props1%cp_030 = props%cp_030
      props1%cp_031 = props%cp_031
      props1%cp_032 = props%cp_032
      props1%cp_033 = props%cp_033
      props1%cp_034 = props%cp_034
      props1%cp_035 = props%cp_035
      props1%cp_036 = props%cp_036
      props1%cp_037 = props%cp_037
      props1%cp_038 = props%cp_038
      props1%cp_039 = props%cp_039
      props1%cp_040 = props%cp_040
      props1%cp_041 = props%cp_041
      props1%cp_042 = props%cp_042
      props1%cp_043 = props%cp_043
      props1%cp_044 = props%cp_044
      props1%cp_045 = props%cp_045
      props1%cp_046 = props%cp_046
      props1%cp_047 = props%cp_047
      props1%cp_048 = props%cp_048
      props1%cp_049 = props%cp_049
      props1%cp_050 = props%cp_050
      props1%cp_051 = props%cp_051
      props1%cp_052 = props%cp_052
      props1%cp_053 = props%cp_053
      props1%cp_054 = props%cp_054
      props1%cp_055 = props%cp_055
      props1%cp_056 = props%cp_056
      props1%cp_057 = props%cp_057
      props1%cp_058 = props%cp_058
      props1%cp_059 = props%cp_059
      props1%cp_060 = props%cp_060
      props1%cp_061 = props%cp_061
      props1%cp_062 = props%cp_062
      props1%cp_063 = props%cp_063
      props1%cp_064 = props%cp_064
      props1%cp_065 = props%cp_065
      props1%cp_066 = props%cp_066
      props1%cp_067 = props%cp_067
      props1%cp_068 = props%cp_068
      props1%cp_069 = props%cp_069
      props1%cp_070 = props%cp_070
      props1%cp_071 = props%cp_071
      props1%cp_072 = props%cp_072
      props1%cp_073 = props%cp_073
      props1%cp_074 = props%cp_074
      props1%cp_075 = props%cp_075
      props1%cp_076 = props%cp_076
      props1%cp_077 = props%cp_077
      props1%cp_078 = props%cp_078
      props1%cp_079 = props%cp_079
      props1%cp_080 = props%cp_080
      props1%cp_081 = props%cp_081
      props1%cp_082 = props%cp_082
      props1%cp_083 = props%cp_083
      props1%cp_084 = props%cp_084
      props1%cp_085 = props%cp_085
      props1%cp_086 = props%cp_086
      props1%cp_087 = props%cp_087
      props1%cp_088 = props%cp_088
      props1%cp_089 = props%cp_089
      props1%cp_090 = props%cp_090
      props1%cp_091 = props%cp_091
      props1%cp_092 = props%cp_092
      props1%cp_093 = props%cp_093
      props1%cp_094 = props%cp_094
      props1%cp_095 = props%cp_095
      props1%cp_096 = props%cp_096
      props1%cp_097 = props%cp_097
      props1%cp_098 = props%cp_098
      props1%cp_099 = props%cp_099
      props1%cp_100 = props%cp_100

c         constants for use in material models
      if( allocated(props%Gmat) ) then
c     Gmat: First full determine size of allocated array
      dim1 = size(props%Gmat,1)
      dim2 = size(props%Gmat,2)
      sizeArr = dim1*dim2
c     Allocate G
      allocate( props1%Gmat(dim1,dim2), stat=allocate_status)
      if( allocate_status .ne. 0 ) then
            write(*,*) ' error allocating G matrix inside solverB'
            call die_gracefully
      endif
c     Copy the contents
      call mm10_c_copy_darray(props%Gmat,sizeArr,props1%Gmat)
      end if

      if( allocated(props%Hmat) ) then
c     Hmat: First full determine size of allocated array
      dim1 = size(props%Hmat,1)
      dim2 = size(props%Hmat,2)
      sizeArr = dim1*dim2
c     Allocate G
      allocate( props1%Hmat(dim1,dim2), stat=allocate_status)
      if( allocate_status .ne. 0 ) then
            write(*,*) ' error allocating G matrix inside solverB'
            call die_gracefully
      endif
c     Copy the contents
      call mm10_c_copy_darray(props%Hmat,sizeArr,props1%Hmat)
      end if

      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_delete_state( props1 )
c
      use mm10_defs
      implicit none
c
      type(crystal_props) :: props1
c
      if( allocated(props1%Gmat) ) then
      deallocate( props1%Gmat  )
      end if
      if( allocated(props1%Hmat) ) then
      deallocate( props1%Hmat )
      end if
c
      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_copy_state( n1, n )

      use mm10_defs
      implicit none
c
      type(crystal_state) :: n,n1

      n1%R = n%R
      n1%gradFeinv = n%gradFeinv
      n1%temp = n%temp
      n1%tinc = n%tinc
      n1%dg = n%dg
      n1%tau_v = n%tau_v
      n1%tau_y = n%tau_y
      n1%tau_l = n%tau_l
      n1%mu_harden = n%mu_harden
      n1%stress = n%stress
      n1%euler_angles = n%euler_angles
      n1%Rp = n%Rp
      n1%D = n%D
      n1%eps = n%eps
      n1%slip_incs = n%slip_incs
      n1%tau_tilde = n%tau_tilde
      n1%u = n%u
      n1%tt_rate = n%tt_rate
      n1%step = n%step
      n1%elem = n%elem
      n1%gp = n%gp
      n1%iter = n%iter
      n1%ms = n%ms
      n1%qs = n%qs
      n1%qc = n%qc
      n1%ed = n%ed
      n1%ep = n%ep
      return
      end

c-----------------------------------------------------------------------

      subroutine mm10_c_copy_darray(A,n,B)
c
      integer :: n, i
      double precision, dimension(n) :: A,B
c
      B(1:n) = A(1:n)
c
      return
      end subroutine

