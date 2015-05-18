c *********************************************************************
c *                                                                   *
c *                           subroutine cnst4                        *
c *                                                                   *
c *             written by : A. roy, S. roychowdhury & others         *
c *                                                                   *
c *                last modified by : 7/22/2014 rhd                   *
c *                                                                   *
c *     tangent modulus matrix for cohesive zone model (linear,       *
c *     exp1_intf, ppr, creep cavity options)                         *
c *                                                                   *
c *********************************************************************
c
      subroutine cnst4(
     &   step, iter, felem, gpn, iout, span, mxvl, time_n, dtime,
     &   nonlocal, numthreads, now_thread, w, c_type, intfprps,
     &   reladis, history, history1,
     &   cep, dj, temp_ref, dtemp, temp_n,
     &   top_surf_elems, bott_surf_elems, top_surf_stresses,
     &   bott_surf_stresses, top_surf_eps, bott_surf_eps,
     &   top_nonlocal_vars, bott_nonlocal_vars, top_solid_matl,
     &   bott_solid_matl )
c
      implicit none
c
c                   parameter declarations
c
      integer step, iter, felem, gpn, iout, span, mxvl, numthreads,
     a        now_thread, c_type, top_surf_elems(span),
     b        bott_surf_elems(span), top_solid_matl(span),
     c        bott_solid_matl(span)
c
#dbl      double precision
#sgl      real
     a reladis(mxvl,*), history(span,*),
     b history1(span,*), cep(mxvl,6,6),
     c intfprps(mxvl,*), dj(*), w, time_n, dtime,
     d temp_ref(*), dtemp(*), temp_n(*),
     e top_surf_stresses(mxvl,6), bott_surf_stresses(mxvl,6),
     f top_surf_eps(mxvl,6), bott_surf_eps(mxvl,6),
     g top_nonlocal_vars(mxvl,*), bott_nonlocal_vars(mxvl,*)
c
      logical nonlocal
c
c
c                [D] updating for cohesive material models.
c                routine is called with data for a full block of
c                elements for a single integration point.
c                this framework of processing is very similar
c                to that used by Abaqus Explicit.
c
c (input)     step:  global solution is advancing analysis from
c                    n -> n+1 (WARP3D load step, Abaqus standard
c                    increment, time step in dynamic analysis).
c                    In WARP3D, the first step in an analysis
c                    is 0 - > 1. Value of step is actuall n+1. At start
c                    of analysis, step = 1.
c
c (input)     iter:  equilibrium iteration number at current step.
c
c (input)    felem:  number of first element (global numbering) in
c                    this block. elements in model are numbered
c                    sequentially.
c
c (input)      gpn:  integration point number for this block of
c                    cohesive elements
c
c (input)     iout:  Fortran device number for output of any messages
c                    from the material model
c
c (inout)      span: number of coehsive elements in this block
c
c (input)      mxvl: maximum allowable number of elements per block
c                    (dimensioned number of rows for many matrices)
c
c (input)    time_n: simulation time at start of load increment
c                    (sum of all delta times)
c
c (input)    dtime:  time increment for this step
c
c (input)  nonlocal: .true. if the nonlocal formulation is in
c                    effect (means info about solid element connected
c                    top & bottom surfaces is passed)
c
c (input) numthreads: number of threads being used to process element
c                     blocks
c (input) now_thread: thread number executing this bock of cohesive
c                     elements
c
c (input)      w:    weight value for this integration point.
c                    Must be multiplied into [cep]
c
c (input)    c_type: the model supports various options for the cohesive
c                    formulation (= type )
c                      1 - linear elastic
c                      2 - bilinear -- not impelemented
c                      3 - ramp     -- not impelemented
c                      4 - exponential_1
c                      5 - exponential_2 -- not impelemented
c                      6 - PPR (Park, Paulino, Roesler)
c                      7 - Cavitation-based cohesive element, referred
c                          to as Cavit.
c                          Reference: Onck P, Van Der Giessen E.
c                          Growth of an initially sharp crack by
c                          grain boundary cavitation.
c                          J. Mech. Phys. Solids 1998;47:99–139
c
c (input) intfprps:  material properties and computational options
c                    at this integration point for all elements in the
c                    block. see mm04_init for detailed description of
c                    each column. most often the properties are
c                    identical for all elements in the block but
c                    this is not required.
c
c (input)   reladis: current estimate of displacement jumps between
c                    n=0 -> n+1 (ds1, ds2, dn). used for linear
c                    elastic and secant formulations
c                     ( ,1):  ds1 (shear sliding 1)
c                     ( ,2):  ds2 (shear sliding 2)
c                     ( ,3):  dn  (normal)
c                    ds1 and ds2 are orthogonal and in the tangent plane
c
c (input)   history: state variables at step n. 7 values per integration
c                    point per element. passed in to routine.
c
c (input)   history1: state variables at step n+1
c
c                      history content for exponential model (type = 4 )
c                      data for state 'n' and 'n+1'
c                          1 -- effective displacement
c                          2 -- maximum effective displacement
c                          3 -- effective traction
c                          4 -- not used
c
c                      history content for PPR (type = 6 )
c                      data for state 'n' and 'n+1'
c                          1 -- dn_limit: normal separation when
c                               normal traction falls to zero
c                          2 -- dt_limit: tangential separation when
c                               shear traction falls to zero
c                          3 -- deln_max: maximum normal separation
c                               during the loading history
c                          4 -- delt_max: maximum tangential separation
c                               during the loading history
c
c
c                      history content for Cavit (c_type = 7 ) data
c                      for state 'n+1'  -- see mm04.f
c
c                   ** updates to history, history1 are discarded by
c                      WARP3D **
c
c (output)   cep:   insert [D_T] for each element in the block.
c                   note dimensioning: mxvl x nstr x nstr
c                   WARP3D sets nstr = 6 even though the [D] is 3x3
c                   for cohesive -- just to keep sizes uniform
c                   with solid elements.
c
c (input)      dj:   determinant of Jacobian at this integration point
c                    for each element into block. Must be multiplied
c                    into [cep]
c
c (input)  temp_ref: user specified initial (reference) temperature for
c                    each element in  block
c
c (input)  dtemp:    temperature increment for step n -> n+1 for
c                    each element in  block
c
c (input)  temp_n:   temperature at start of step for each element in
c                    block. = initial temperature + all prior changes.
c
c Following data items are passed as actual values only if
c nonlocal = .true.. strains and stresses in solid elements are
c for the current estimate of the solution at n+1
c
c (input)  top_surf_elems: number of solid element connected to top
c                          surface for each element in  block
c
c (input) bott_surf_elems: number of solid element connected to bottom
c                          surface for each element in  block
c
c (input) top_surf_stresses: stresses for solid element attached to
c                          top surface of each element in block. 6
c                          components in model (global) coordinates.
c                          xx, yy, zz, zy, yz, xz. These are averaged
c                          centroidal values for solids.
c
c (input) bott_surf_stresses: stresses for solid element attached to
c                          bottom surface of each element in block. 6
c                          components in model (global) coordinates.
c                          xx, yy, zz, zy, yz, xz. These are averaged
c                          centroidal values for solids.
c
c (input) top_surf_eps:    strains for solid element attached to
c                          top surface of each element in block. 6
c                          components in model (global) coordinates.
c                          xx, yy, zz, zy, yz, xz. These are averaged
c                          centroidal values for solids. shears are
c                          gamma not epsilon
c
c (input) bott_surf_eps:   strains for solid element attached to
c                          bottom surface of each element in block. 6
c                          components in model (global) coordinates.
c                          xx, yy, zz, zy, yz, xz. These are averaged
c                          centroidal values for solids. shears are
c                          gamma not epsilon
c
c (input) top_nonlocal_vars: additional values from solid element
c                          material models passed to here for possible
c                          use by cohesive update. values and order of
c                          values set by the solid material model.
c                          These are average of integration point
c                          values for top solid elements.
c
c (input) bott_nonlocal_vars: additional values from solid element
c                          material models passed to here for possible
c                          use by cohesive update. values and order of
c                          values set by the solid material model.
c                          These are average of integration point
c                          values for bottom solid elements.
c
c (input) top_matl_model:  WARP3D material model type for top
c                          solid elements
c
c (input) bott_matl_model: WARP3D material model type for bottom
c                          solid elements
c
c                          Material model types:
c                          1 = bilinear, 2 = deform. plasticity
c                          3 = mises/gurson, 4 = this model
c                          5 = general cyclic, 6 = <none>
c                          7 = mises + H2, 8 = umat
c
c  Note: when a surface of a cohesive element is attached to a symmetry
c        plane that solid element = 0 above. then stresses have been
c        made equal for top and bottom solids. same for strains.
c
c      for a large displacement analysis, the area change of the
c      cohesive surface has already been computed by the corresponding
c      interface-cohesive element. this material model doesn't
c      see the difference between a small and large displacement
c      analysis - exactly as other material models in WARP3D
c
c ------------------------------------------------------------------
c
c                   locally defined
c
      integer i, info_vector(20), history_length
#dbl      double precision
#sgl      real
     & zero, e, tol, three, initial,
     & intfmat(mxvl,3,3), ppr_support(mxvl,20), one,
     & top_surf_mean_stress, top_surf_mean_eps,
     & bott_surf_mean_stress, bott_surf_mean_eps,
     & dn, comp_multiplier
c
      logical elem_killed(mxvl), local_debug, debug_ppr,
     &        debug_exp1
c
      data zero, e, tol, one, three, initial
     &  / 0.0d0, 2.71828182845904523536d0, 0.1d0, 1.0d0, 3.0d0,
     &  1.d-14 /
c
c
      local_debug = intfprps(1,29) .gt. zero
c      
      call mm04_set_sizes( info_vector ) 
      history_length = info_vector(1)  !  note this so we avoid bugs
c                                         if mm04 changes history size
c
c             debug output on entry to cnst4
c
      if( local_debug ) then
         write(iout,9500) step, iter, c_type, span, mxvl, gpn,
     &                    felem, time_n, dtime, w,  nonlocal,
     &                    numthreads, now_thread
         write(iout,9510)
         do i = 1, span
          write(iout,9512) i, felem+i-1, temp_ref(i), dtemp(i),
     1                     temp_n(i)
         end do
         write(iout,9520)
         do i = 1, span
           write(iout,9522) i, felem+i-1, dj(i)
         end do
      end if
c
      if( nonlocal .and. local_debug ) then
          write(iout,9620)
          do i = 1, span
            top_surf_mean_stress = (top_surf_stresses(i,1) +
     &                  top_surf_stresses(i,2) +
     &                  top_surf_stresses(i,3)) / three
            bott_surf_mean_stress = (bott_surf_stresses(i,1) +
     &                  bott_surf_stresses(i,2) +
     &                  bott_surf_stresses(i,3)) / three
            top_surf_mean_eps    = (top_surf_eps(i,1) +
     &                  top_surf_eps(i,2) +
     &                  top_surf_eps(i,3)) / three
            bott_surf_mean_eps   = (bott_surf_eps(i,1) +
     &                  bott_surf_eps(i,2) +
     &                  bott_surf_eps(i,3)) / three
            write(iout,9622) i, felem+i-1, top_surf_elems(i),
     &                  bott_surf_elems(i), top_surf_mean_stress,
     &                  top_surf_mean_eps,
     &                  bott_surf_mean_stress, bott_surf_mean_eps
          end do
          write(iout,9625)
          do i = 1, span
            write(iout,9627) i, felem+i-1, top_solid_matl(i),
     &                       bott_solid_matl(i)
            write(iout,9629) top_nonlocal_vars(i,1:5),
     &                       bott_nonlocal_vars(i,1:5)
          end do
          write(iout,*) " "
      end if
c
      do i = 1, span
        elem_killed(i) = intfprps(i,13) .gt. zero
        cep(i,1:3,1:3) = zero
      end do
c
c           handle each type of cohesive traction-separation type
c           separately
c
      select case ( c_type )
c     ----------------------
c
      case( 1 )
c     =========
c
c             Option 1:  linear-elastic option
c
c                    intfprps( ,2) - interface stiffness in the
c                                    longitudinal direction
c                    intfprps( ,3) - interface stiffness in the
c                                    transverse direction
c                    intfprps( ,4) - interface stiffness in the
c                                    normal direction
c
c                    (1,1) = (2,2) implies isotropic shear-sliding
c                                  stiffness.
c
c             For interface-cohesive elements connected to a symmetry
c             plane, the user should specify *doubled* values of the
c             stiffness.
c
      do i = 1,span
        cep(i,1,1) = intfprps(i,2) * dj(i) * w
        cep(i,2,2) = intfprps(i,3) * dj(i) * w
        cep(i,3,3) = intfprps(i,4) * dj(i) * w
        dn = reladis(i,3)        
        if( dn .lt. zero ) then
           comp_multiplier = intfprps(i,22)
           cep(i,3,3) = cep(i,3,3) * comp_multiplier
        end if
      end do
c
      case( 2, 3  )
c     =============
c
c             Option 2 & 3:  bilinear, ramp options -- not implemented
c
      write(iout,9200)
      call die_abort
c
c
      case( 4 )
c     =========
c
c             Option 4:  standard, 2 parameter exponential option
c             an exponential traction - separation law is used in this
c             cohesive zone model.
c
      debug_exp1 = .false.
      call cnst4_exp1_tan( span, mxvl, felem, elem_killed,
     &                        intfprps, reladis,
     &                        history, history1, intfmat,
     &                        debug_exp1, iout )
c
      do i = 1, span
        cep(i,1:3,1:3) = intfmat(i,1:3,1:3) * dj(i) * w
      end do
c
      case( 5 )
c     =========
c
c             Option 5:  alternate exponential option --
c             not implemented
c
      write(iout,9110)
      call die_abort
c
      case( 6 )
c     =========
c
c             PPR cohesive model.
c             Ref. JMPS 57 (6), 891-908 (2009)
c               1) Initialize the model parameters. put in ppr_support
c               2) Compute tangent [D] matrix
c
      debug_ppr = .false.
      call mm04_init_ppr( span, intfprps, ppr_support,
     &                    elem_killed, iout, mxvl )
      call cnst4_ppr( span, ppr_support, reladis, history,
     &                history1, intfmat, elem_killed,
     &                debug_ppr, felem, iout, mxvl )
      do i = 1, span
        cep(i,1:3,1:3) = intfmat(i,1:3,1:3) * dj(i) * w
      end do
c
      case ( 7 )
c     =========
c
c             Cavit cohesive model.
c             JMPS. 1998;47:99–139. starting point for refs
c               1) Initialize the model parameters.
c               2) Compute tangent [D] matrix
c
      call cnst4_cavit( step, iter, span, felem, iout, mxvl, 
     &                  dtime, nonlocal, history, history1, intfprps,
     &                  intfmat, reladis, elem_killed )
c
      do i = 1, span
        cep(i,1:3,1:3) = intfmat(i,1:3,1:3) * dj(i) * w
      end do
c      
      case default
c     ============
c
      write(iout,9120)
      call die_abort
c
      end select
c
      if( local_debug ) then
        write(iout,9530)
        do i = 1, span
          write(iout,9532) i, felem+i-1, cep(i,1,1:3),
     &                     cep(i,2,1:3), cep(i,3,1:3)
        end do
       end if
c
       return
c
 9200 format('>>>> FATAL ERROR: bilinear option not implemented',
     & /,    '                  use ppr or exponential option.',
     & /,    '                  job terminated by cnst4' )
 9110 format('>>>> FATAL ERROR: exponential 2 model not implemented',
     & /,    '                  use ppr or exponential option.',
     & /,    '                  job terminated by cnst4' )
 9120 format('>>>> FATAL ERROR: unknown cohesive formulation',
     & /,    '                  job terminated by cnst4' )
 9500 format(/," ...... entered cnst4 .....",
     1 /,10x,"step, iter, c_type, span, mxvl: ",i6,2i3,i4,i4,
     2 /,10x,"gpn, felem, time_n, dtime: ",i3,i8, e14.6,1x,e14.6,
     3 /,10x,"weight(w): ",f9.6,
     4 /,10x,"nonlocal, numthreads, now_thread: ",l1,2i4)
 9510 format(/,10x,"temperature data:",
     1 /,9x,"i  ",2x,"element",5x,"ref temp",8x,"dtemp",8x,"temp_n" )
 9512 format(7x,i3,i8,3f15.4)
 9520 format(/,10x,"Jacobian determinant:",
     1 /,9x,"i  ",2x,"element",10x,"dj" )
 9522 format(7x,i3,i8,f15.4)
 9530 format(/,10x,"[D] matrices with w and dj:",
     1 /,9x,"i  ",2x,"element",15x,"----- [cep] -----")
 9532 format(7x,i3,i8,2x,3f15.4,/20x,3f15.4,/,20x,3f15.4)
 9620 format(/,10x,"nonlocal data (@ estimate for end of step):",
     1 /,9x,"i  ",2x,"element",5x,"top solid",8x,"bottom solid",
     2  3x,"top sig mean",5x,"top eps mean",
     3  8x,"bott sig mean",5x,"bott eps mean" )
 9622 format(7x,i3,2x,i8,5x,i10,8x,i10,2x,2f15.6,7x,2f15.6)
 9625 format(/,10x,"more nonlocal data (material type, solid state",
     a " variables passed thru @ estimate for end of step):",
     b /,9x,"i  ",2x,"element", 2x,"top matl type",5x,
     c   "bottom matl type" )
 9627 format(7x,i3,2x,i8,2x,i8,5x,i10)
 9629 format(12x,"top solid vars @ n+1:    ",5f10.5,
     a /,12x,    "bottom solid vars @ n+1: ",5f10.5 )
c
      end
c
c    ****************************************************************
c    *                                                              *
c    *               subroutine cnst4_exp1_tan                      *
c    *                                                              *
c    *                       written by : A. Roy                    *
c    *                                    S. Roychowdhury           *
c    *                                                              *
c    *                    last modified : 7/3/2013 rhd              *
c    *                                                              *
c    *     computes the tangent [D] for a block of cohesive         *
c    *     elements for the exponential traction-seperation law     *
c    *     (homogeneous material)                                   *
c    *     Ref. IJNME 44,1267-1282 (1999)                           *
c    *                                                              *
c    ****************************************************************
c
      subroutine cnst4_exp1_tan( span, mxvl, felem, elem_killed,
     &                           intfprps, dis, history,
     &                           history1, tanmat, local_debug, iout )
      implicit integer (a-z)
c
c                   parameter declarations
c
#dbl      double precision
#sgl      real
     & intfprps(mxvl,*), tanmat(mxvl,3,3), dis(mxvl,*),
     & history(span,*), history1(span,*)
       logical local_debug, elem_killed(*)
c
c                   locally defined
c
#dbl      double precision
#sgl      real
     & zero, tol, dtol, one, t11, t12, t13, t21, t22, t23, t31, t32,
     & t33, effdis(mxvl), efftrac, e, maxeffdis(mxvl),
     & dtd1, dtd2, dtd3, b2, b4, eff3,
     & dn, dt1, dt2, beta, prior_max_d_eff, d_eff_at_n,
     & d_eff_at_peak, peak_intf_stress, comp_multiplier
       logical compression(mxvl), small_effdis(mxvl), loading(mxvl)
c
       data zero, one, e, dtol
     & / 0.0d0, 1.0d0, 2.71828182845904523536d0, 0.000001d0  /
c
      do i = 1, span
        tanmat(i,1:3,1:3) = zero
      end do
c
c         explanation of logical variables
c
c         compression = .true.   =>
c                       jump normal displacement (dn) from
c                       n (time) = 0 to current n (time) < 0
c
c         small_effdis = .true.  =>
c                       current effective displacement (deff) is very
c                       small compared to displacement at peak
c                       stress on traction-separation curve. used
c                       to prevent divide by zeros
c
c         loading = .true. if current deff > prior_max_d_eff
c                   .and. deff is > value at start of step.
c                   material point is being pushed farther out
c                   the traction separation curve
c
c         when small_effdis = .true.:
c           normal_stiffness = slope of normal t-s curve at origin
c           tangential_stiffness = slope of tangential t-s curve
c                                  at origin
c
c         when compression = .true.:
c             normal_stiffness = compression multiplier *
c                                slope of normal t-s curve at origin
c
c         referenced history variables at start of step
c                   1 - effective jump displacement
c                   2 - prior maximum effective displacement
c
      if( local_debug ) then
        write(iout,9000); write(iout,9005)
      end if
      do i = 1, span
        if( elem_killed(i) ) cycle
        d_eff_at_peak   = intfprps(i,11)
        beta            = intfprps(i,12)
        d_eff_at_n      = history(i,1)
        prior_max_d_eff = history(i,2)
        tol             = dtol * d_eff_at_peak
        dt1 = dis(i,1); dt2 = dis(i,2); dn = dis(i,3)
c
        compression(i) = dn .lt. zero
        if( compression(i) ) dn = zero
        effdis(i) = sqrt( beta*beta*(dt1*dt1+dt2*dt2) + dn*dn )
        maxeffdis(i) = prior_max_d_eff
        if( effdis(i) .gt. prior_max_d_eff ) maxeffdis(i) = effdis(i)
        small_effdis(i) = effdis(i) .le. tol
        loading(i) = effdis(i) .ge. d_eff_at_n   .and.
     &               effdis(i) .eq. maxeffdis(i)
        if( local_debug )
     &   write(iout,9010) felem+i-1, effdis(i),  compression(i),
     &          small_effdis(i), loading(i), d_eff_at_peak,
     &          d_eff_at_n, prior_max_d_eff, tol, dt1, dt2, dn
      end do
c
c
      do k = 1, span
c
        if( elem_killed(k) ) cycle
c
        peak_intf_stress = intfprps(k,5)
        d_eff_at_peak    = intfprps(k,11)
        beta             = intfprps(k,12)
        comp_multiplier  = intfprps(k,22)
c
        if( small_effdis(k) ) then
           tanmat(k,1,1) = e*beta**2*peak_intf_stress
     &                     /d_eff_at_peak
           tanmat(k,2,2) = tanmat(k,1,1)
           tanmat(k,3,3) = e*peak_intf_stress/d_eff_at_peak
           if( compression(k) )
     &       tanmat(k,3,3) = comp_multiplier * tanmat(k,3,3)
           cycle
        end if
c
c                unloading case -- unload towards origin.
c                can still have compression
c
        if( .not. loading(k) ) then
           tanmat(k,1,1) = e*beta**2*peak_intf_stress/d_eff_at_peak
     &                     *exp(-one*maxeffdis(k)/d_eff_at_peak)
           tanmat(k,2,2) = tanmat(k,1,1)
           tanmat(k,3,3) = tanmat(k,1,1)
           if( compression(k) ) tanmat(k,3,3) =
     &        comp_multiplier * e*peak_intf_stress/d_eff_at_peak
           cycle
        end if
c
c                loading case -- can still have compression (rare)
c
        efftrac = e*peak_intf_stress*effdis(k)/d_eff_at_peak*
     &            e**(-one*effdis(k)/d_eff_at_peak)
        dt1  = dis(k,1)  ! jump sliding 1
        dt2  = dis(k,2)  ! jump sliding 2
        dn   = dis(k,3)  ! jump normal
        b2   = beta * beta
        b4   = b2 * b2
        eff3 = effdis(k)**3
c
        dtd1 = e*b2*peak_intf_stress*dt1*
     &            exp(-one*effdis(k)/d_eff_at_peak)*
     &            (one-effdis(k)/d_eff_at_peak)/d_eff_at_peak
     &            /effdis(k)
        dtd2 = e*b2*peak_intf_stress*dt2*
     &            exp(-one*effdis(k)/d_eff_at_peak)*
     &            (one-effdis(k)/d_eff_at_peak)/d_eff_at_peak
     &            /effdis(k)
        dtd3 = e*peak_intf_stress*dn*
     &            exp(-one*effdis(k)/d_eff_at_peak)*
     &            (one-effdis(k)/d_eff_at_peak)/d_eff_at_peak
     &            /effdis(k)
        t11 = b2*efftrac/effdis(k)+ b2*dt1*dtd1/effdis(k)-
     &           b4*efftrac*dt1*dt1/eff3
        t12 = b2*dt1*dtd2/effdis(k)- b4*efftrac*dt1*dt2/eff3
        t21 = b2*dt2*dtd1/effdis(k)- b4*efftrac*dt2*dt1/eff3
        t22 = b2*efftrac/effdis(k)+ b2*dt2*dtd2/effdis(k)-
     &           b4*efftrac*dt2*dt2/eff3
        t13 = b2*dt1*dtd3/effdis(k)- b2*efftrac*dt1*dn/eff3
        t23 = b2*dt2*dtd3/effdis(k)- b2*efftrac*dt2*dn/eff3
        t31 = dn*dtd1/effdis(k)- b2*efftrac*dn*dt1/eff3
        t32 = dn*dtd2/effdis(k)- b2*efftrac*dn*dt2/eff3
        t33 = efftrac/effdis(k) + dn*dtd3/effdis(k)- efftrac*dn*dn/eff3
        if( compression(k) ) then
           t13 = zero; t23 = zero; t31 = zero; t32 = zero
           t33 = comp_multiplier * e*peak_intf_stress/d_eff_at_peak
        end if
        tanmat(k,1,1) = t11
        tanmat(k,1,2) = t12
        tanmat(k,2,1) = t21
        tanmat(k,2,2) = t22
        tanmat(k,1,3) = t13
        tanmat(k,2,3) = t23
        tanmat(k,3,1) = t31
        tanmat(k,3,2) = t32
        tanmat(k,3,3) = t33
c
      end do
c
      return
c
 9000 format(/,10x,'... inside cnst4_exp1_tan. after setup ...')
 9005 format(12x,'elem      effdis',6x,'comp flg',2x,'intf_comp flg',
     & 2x,'load flg',2x,'d_eff_at_peak',2x,'d_eff_at_n',2x,
     & 'prior_max_d_eff',4x,'tol', 12x, 'dt1', 10x,'dt2',11x,'dn' )
 9010 format(7x,i9, 2x,f12.8, l8, 4x,l8, 4x,l8,5x,f12.8,
     & 3x, f12.8, 2x, f12.8, 2x, f12.8, 2x, f12.8, 2x,f12.8,
     & 2x,f12.8)
c
      end
c    ****************************************************************
c    *                                                              *
c    *               subroutine cnst4_ppr                           *
c    *                                                              *
c    *          written by : Kyoungsoo Park 7/25/2010               *
c    *             updated : R. Dodds 8/31/2010; 11/21/2010         *
c    *                                                              *
c    *     computes the tangent stiffness [D_T] for a  block        *
c    *     cohesive elements for the unified potential-based        *
c    *     model (PPR)                                              *
c    *     Ref. JMPS 57 (6), 891-908 (2009)                         *
c    *                                                              *
c    ****************************************************************
c
      subroutine cnst4_ppr( span, ppr_support, del, history,
     &                      history1, tanmat, elem_killed,
     &                      local_debug, felem, iout, mxvl )
      implicit none
      integer :: i, span, mxvl, felem, iout
#dbl      double precision
#sgl      real
     &        ppr_support(mxvl,*), del(mxvl,*), tanmat(mxvl,3,3),
     &        history(span,*), history1(span,*)
      logical elem_killed(*), local_debug
c
#dbl      double precision
#sgl      real
     &        Gam_n, Gam_t, Tn_m, Tt_m, dn, dt, alph, beta, m, n,
     &        dGtn, dGnt, dnb, dtb, comp_multiplier,
     &        deln, delt, deln_max, delt_max, Dnn, Dnt, Dtn, Dtt, Tt,
     &        zero, tol, t(3,3), half, one, two, dtol
      logical :: compression
      data zero, tol, one, two, half
     & / 0.0d00, 0.00001d00, 1.0d00, 2.0d00, 0.5d00 /
c
      do i = 1, span
        if( .not. elem_killed(i) ) cycle
        tanmat(i,1:3,1:3) = zero
      end do
c
c            main loop to get [D]
c            --------------------
c
      do i = 1, span
       if ( elem_killed(i) ) cycle
c
c            extract material props, constants
c
        Tn_m  = ppr_support(i,4)
        Tt_m  = ppr_support(i,5)
        alph  = ppr_support(i,6)
        beta  = ppr_support(i,7)
        m     = ppr_support(i,10)
        n     = ppr_support(i,11)
        dn    = ppr_support(i,12)
        dt    = ppr_support(i,13)
        Gam_n = ppr_support(i,14)
        Gam_t = ppr_support(i,15)
        dGnt  = ppr_support(i,16)
        dGtn  = ppr_support(i,17)
        dnb   = ppr_support(i,18)
        dtb   = ppr_support(i,19)
        delt  = sqrt( del(i,1)**2 + del(i,2)**2 )
        deln  = del(i,3)
        dtol  = tol * dt
        deln_max = history(i,3)
        delt_max = history(i,4)
        compression = deln .lt. zero
c
        if( local_debug )
     &     write(iout,9000) felem + i - 1, del(i,1:3), dtol, deln_max,
     &                      delt_max, compression, dn, dt
c
c         Compute normal components (Dnn, Dnt): 4 cases
c         Case 1: Compression region
c
        if( compression ) then
          Dnn = -Gam_n/dn**2*(m/alph)**(m-one)*(alph+m)*
     & (Gam_t*(one-delt/dt)**beta*(n/beta+delt/dt)**n + dGtn)
          comp_multiplier  = ppr_support(i,20)
          Dnn = Dnn * comp_multiplier
          Dnt = zero
c
c         KP: note: - may need to consider different way to define Dnt
c
          deln = zero    !   ?????
c
c         Case 2: Complete failure
c
        elseif ((deln .ge. dn) .or. (delt .ge. dtb)
     &           .or. (deln_max .ge. dn)) then
          Dnn = zero
          Dnt = zero
c
c         Case 3: Softening region
c
        elseif( deln .ge. deln_max ) then
          Dnn =
     & (Gam_t*(one-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn**2 *
     & ((one-deln/dn)**(alph-two)*(alph-one)*alph*(deln/dn+m/alph)**m -
     & two*(one-deln/dn)**(alph-one)*alph*(deln/dn+m/alph)**(m-one)*m +
     & (one-deln/dn)**alph*(deln/dn+m/alph)**(m-two)*(m-one)*m)
          Dnt =
     & Gam_t/dt*(-(one-delt/dt)**(beta-one)*beta*(delt/dt+n/beta)**n +
     & (one-delt/dt)**beta*(delt/dt+n/beta)**(n-one)*n) *
     & Gam_n/dn*(-(one-deln/dn)**(alph-one)*alph*(deln/dn+m/alph)**m +
     & (one-deln/dn)**alph*(deln/dn+m/alph)**(m-one)*m)
c
c        Case 4: Unloading/reloading condition
        else
          Dnn =
     & (Gam_t*(one-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*((one-deln_max/dn)**alph*(deln_max/dn+m/alph)**(m-one)*m
     & -(one-deln_max/dn)**(alph-one)*alph*(deln_max/dn+m/alph)**m)
     & / deln_max
          Dnt =
     & Gam_t/dt*(-(one-delt/dt)**(beta-one)*beta*(delt/dt+n/beta)**n +
     & (one-delt/dt)**beta*(delt/dt+n/beta)**(n-one)*n) *
     & Gam_n/dn*(m*(one-deln_max/dn)**alph*(m/alph+deln_max/dn)**(m-one)
     &   -alph*(one-deln_max/dn)**(alph-one)*(m/alph+deln_max/dn)**m) *
     & deln/deln_max
        end if
c
c        Compute tangential components (Dtn, Dtt): 3 cases
c        Case 1: Complete failure
c
        if ((delt .ge. dt) .or. (deln .ge. dnb)
     &       .or. (delt_max .ge. dt))  then
          Tt = zero
          Dtt = zero
          Dtn = zero
c
c        Case 2: Softening region
c
        elseif (delt .ge. delt_max) then
          Tt = (Gam_n*(one-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(one-delt/dt)**beta*(delt/dt+n/alph)**(n-one)
     &       -beta*(one-delt/dt)**(beta-one)*(delt/dt+n/beta)**n)
          Dtt =
     & (Gam_n*(one-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt**2 *
     & ((one-delt/dt)**(beta-2)*(beta-one)*beta*(delt/dt+n/beta)**n -
     & two*(one-delt/dt)**(beta-one)*beta*(delt/dt+n/beta)**(n-one)*n +
     & (one-delt/dt)**beta*(delt/dt+n/beta)**(n-two)*(n-one)*n)
          Dtn =
     & Gam_t/dt*(-(one-delt/dt)**(beta-one)*beta*(delt/dt+n/beta)**n +
     & (one-delt/dt)**beta*(delt/dt+n/beta)**(n-one)*n) *
     & Gam_n/dn*(-(one-deln/dn)**(alph-one)*alph*(deln/dn+m/alph)**m +
     & (one-deln/dn)**alph*(deln/dn+m/alph)**(m-one)*m)
c
c       Case 3: Unloading/reloading condition
c
        else
          Tt = (Gam_n*(one-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(one-delt_max/dt)**beta*(delt_max/dt+n/alph)**(n-one)
     &   -beta*(one-delt_max/dt)**(beta-one)*(delt_max/dt+n/beta)**n) *
     & delt/delt_max
          Dtt =
     & (Gam_n*(one-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(one-delt_max/dt)**beta*(delt_max/dt+n/alph)**(n-one)
     &      -beta*(one-delt_max/dt)**(beta-one)*(delt_max/dt+n/beta)**n)
     & / delt_max
          Dtn =
     & Gam_n/dn*(-(one-deln/dn)**(alph-one)*alph*(deln/dn+m/alph)**m +
     & (one-deln/dn)**alph*(deln/dn+m/alph)**(m-one)*m) *
     & Gam_t/dt*(n*(one-delt_max/dt)**beta*(delt_max/dt+n/alph)**(n-one)
     &   -beta*(one-delt_max/dt)**(beta-one)*(delt_max/dt+n/beta)**n) *
     & delt/delt_max
        end if
c
c          Store tangent matrix terms into 3x3
c
        if( delt .lt. dtol ) then  ! becomes diagonal. no coupling
          tanmat(i,1,1) = Dtt
          tanmat(i,1,2) = zero
          tanmat(i,1,3) = zero
          tanmat(i,2,1) = zero
          tanmat(i,2,2) = Dtt
          tanmat(i,2,3) = zero
          tanmat(i,3,1) = zero
          tanmat(i,3,2) = zero
          tanmat(i,3,3) = Dnn
        else
          tanmat(i,1,1) = Dtt*(del(i,1)/delt)**2 +
     &                    Tt*del(i,2)**2/delt**3
          tanmat(i,1,2) = Dtt*del(i,1)*del(i,2)/delt**2 -
     &                    Tt*del(i,1)*del(i,2)/delt**3
          tanmat(i,1,3) = Dtn*del(i,1)/delt
          tanmat(i,2,1) = Dtt*del(i,1)*del(i,2)/delt**2 -
     &                    Tt*del(i,1)*del(i,2)/delt**3
          tanmat(i,2,2) = Dtt*(del(i,2)/delt)**2 +
     &                    Tt*del(i,1)**2/delt**3
          tanmat(i,2,3) = Dtn*del(i,2)/delt
          tanmat(i,3,1) = Dnt*del(i,1)/delt
          tanmat(i,3,2) = Dnt*del(i,2)/delt
          tanmat(i,3,3) = Dnn
        endif
      end do
c
c            end of main loop to get [D]
c            ---------------------------
c
c          The [D_T] is non-symmetric when the normal and shear
c          energies differ. Symmetrize the resulting (3x3).
c          Take care of Mode I only case and shear mode only case.
c
      if( abs( Gam_n - Gam_t) .gt. tol*max(Gam_n,Gam_t) ) then
        do i = 1, span
         t(1,1) = tanmat(i,1,1)
         t(1,2) = half * ( tanmat(i,1,2) + tanmat(i,2,1) )
         t(1,3) = half * ( tanmat(i,1,3) + tanmat(i,3,1) )
         t(2,2) = tanmat(i,2,2)
         t(2,1) = tanmat(i,1,2)
         t(2,3) = half * ( tanmat(i,3,2) + tanmat(i,2,3) )
         t(3,3) = tanmat(i,3,3)
         t(3,1) = tanmat(i,1,3)
         t(3,2) = tanmat(i,2,3)
         tanmat(i,1:3,1:3) = t(1:3,1:3)
        end do
      end if
c
      if( local_debug ) then
       do i = 1, span
         write(iout,9100) felem + i - 1, tanmat(i,1,1:3),
     &           tanmat(i,2,1:3),  tanmat(i,3,1:3)
       end do
      end if
c
      return
c
 9000 format('.... inside ppr cnst. element: ',i10,
     &   /,  '     del: ',3f15.10,
     &   /,  '     dtol: ',f15.10,
     &   /,  '     deln_max, delt_max: ',2f15.10,
     &   /,  '     compression: ',l1,
     &   /,  '     dn, dt: ', 2f15.9 )
 9100 format('.... tanmat for element: ',i10, 3(/,3f15.6) )
c
      end
c
c    ****************************************************************
c    *                                                              *
c    *               subroutine cnst4_cavit                         *
c    *                                                              *
c    *               written by rhd 7/22/2014                       *
c    *               updated 10/23/2104 rhd                         *
c    *                                                              *
c    *     computes the tangent stiffness [D_T] for a  block        *
c    *     cohesive elements for the cavitation-based cohesive      *
c    *     model (cavit)                                            *
c    *     Ref. JMPS. 1998, vol 47, pp 99-139                       *
c    *                                                              *
c    ****************************************************************
c
      subroutine cnst4_cavit(
     1 step, iter, span, felem, iout, mxvl, dtime, nonlocal,
     2 history, history1, intfprps, intfmat, reladis, elem_killed )
      implicit none
c
c             parameters
c
      integer :: step, iter, span, iout, mxvl, felem
#dbl      double precision ::
#sgl      real ::
     & history(span,*), history1(span,*), dtime, intfmat(mxvl,3,3),
     & reladis(mxvl,*), intfprps(mxvl,*)
      logical :: nonlocal, elem_killed(*)
c
c             locals 
c
      integer :: i
      logical :: here_debug
#dbl      double precision ::
#sgl      real ::
     & zero 
      data  zero / 0.0d0 /
c
      here_debug = .false.         
c
      do i = 1, span
        intfmat(i,1:3,1:3) = zero
        if( elem_killed(i) ) cycle
        intfmat(i,1,1) =  history1(i,6)  ! shear 
        intfmat(i,2,2) =  history1(i,6)  ! shear 
        intfmat(i,3,3) =  history1(i,7)  ! normal 
      end do   
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *        material model # 4 --  cohesive material model        *
c     *                                                              *
c     *                       written by : aroy                      *
c     *                                                              *
c     *                   last modified : 6/12/2014  rhd             *
c     *                                                              *4
c     *      linear constituitive  matrices [D_L] block of           *
c     *      interface-cohesive elements                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine lcnst4(
     &   step, cohes_type, span, mxvl, nstr, gpn, felem, iout, time_n,
     &   dtime, cep, intfprps, dj, w, temp_np1, temp_ref )
      implicit none
c
      integer step, cohes_type, span, mxvl, nstr, gpn, felem, iout
c
#dbl      double precision
#sgl      real
     & cep(mxvl,nstr,*), intfprps(mxvl,*), dj(*), w,
     & temp_np1(mxvl), temp_ref(mxvl), time_n, dtime
c
c                   locally defined. some are automatic arrays
c
      integer i
#dbl      double precision
#sgl      real
     & zero, e, beta, one, ppr_support(mxvl,20)
      logical elem_killed(mxvl), local_debug
c

      data zero, e, one / 0.0d0, 2.7182818284d0, 1.0d0 /
c
c
c (input)     step:  global solution is advancing analysis from
c                    n -> n+1 (WARP3D load step, Abaqus standard
c                    increment, time step in dynamic analysis).
c                    In WARP3D, the first step in an analysis
c                    is 0 - > 1. Value of step is actuall n+1. At start
c                    of analysis, step = 1.
c
c (input) cohes_type:  the model supports various options for
c                      the cohesive formulation (= type )
c                      1 - linear elastic
c                      2 - bilinear -- not impelemented
c                      3 - ramp     -- not impelemented
c                      4 - exponential_1
c                      5 - exponential_2 -- not impelemented
c                      6 - PPR (Park, Paulino, Roesler)
c                      7 - cavit (cavitation-based cohesive element)
c
c (inout)      span: number of coehsive elements in this block
c
c (input)      mxvl: maximum allowable number of elements per block
c                    (dimensioned number of rows for many matrices)
c
c (input)      nstr: dimension number of cols (WARP3D sets=6)
c
c (input)      gpn:  integration point number for this block of
c                    cohesive elements
c
c (input)     felem: number of first element (global numbering) in
c                    this block. elements in model are numbered
c                    sequentially.
c
c (input)      iout: Fortran device number for output of any messages
c                    from the material model
c
c (input)    time_n: simulation time at start of load increment
c                    (sum of all delta times)
c
c (input)    dtime:  time increment for this step
c
c (output)   cep:   insert [D_L] for each element in the block.
c                   note dimensioning: mxvl x nstr x nstr
c                   WARP3D sets nstr = 6 even though the [D] is 3x3
c                   for cohesive -- just to keep sizes uniform
c                   with solid elements.
c
c (input) intfprps:  material properties and computational options
c                    at this integration point for all elements in the
c                    block. see mm04_init for detailed description of
c                    each column. most often the properties are
c                    identical for all elements in the block but
c                    this is not required.
c
c (input)      dj:   determinant of Jacobian at this integration point
c                    for each element into block. Must be multiplied
c                    into [cep]
c
c (input)      w:    weight value for this integration point.
c                    Must be multiplied into [cep]
c
c (input) temp_ref:  user specified initial (reference) temperature for
c                    each element in  block
c
c (input) temp_np1:  temperature at end of current step for each
c                    element in block. = initial temperature + all
c                    prior changes.
c
c           Compute linear interface [D] (3x3) matrix for each element
c           in the block at this integration point. Linear [D] is
c           diagonal.
c
      local_debug = intfprps(1,29) .gt. zero
c
      if( local_debug ) then
         write(iout,9500) step, cohes_type, span, mxvl, nstr, gpn,
     &                    felem, time_n, dtime, w
         write(iout,9510)
         do i = 1, span
            write(iout,9512) i, felem+i-1, temp_ref(i), temp_np1(i)
         end do
         write(iout,9520)
         do i = 1, span
           write(iout,9522) i, felem+i-1, dj(i)
         end do
      end if
c
      do i = 1, span
        elem_killed(i) = intfprps(i,13) .gt. zero
        cep(i,1:3,1:3) = zero
      end do
c
c           handle each type of cohesive traction-separation type
c           separately
c
      select case ( cohes_type )
      case( 1 )
c     =========
c
c             Option 1:  linear-elastic option
c
c                    intfprps( ,2) - interface stiffness in the
c                                    longitudinal direction
c                    intfprps( ,3) - interface stiffness in the
c                                    transverse direction
c                    intfprps( ,4) - interface stiffness in the
c                                    normal direction
c
c                    (1,1) = (2,2) implies isotropic shear-sliding
c                                  stiffness.
c            For interface-cohesive elements connected ot a symmetry
c            plane, the user should specify *doubled* values of the
c            stiffness.
c
      do i = 1, span
          cep(i,1,1) = intfprps(i,2) * dj(i) * w
          cep(i,2,2) = intfprps(i,3) * dj(i) * w
          cep(i,3,3) = intfprps(i,4) * dj(i) * w
      end do
c
      case( 2, 3  )
c     =============
c
c             Option 2 & 3:  bilinear, ramp options -- not implemented
c
      write(iout,9200)
      call die_abort
c
c
      case( 4 )
c     =========
c
c             exponential function of equivalent displacement to
c             determine effective traction. IJNME 44,1267-1282 (1999)
c             uses properties:
c                   intfprps( ,5)  - effective peak stress on the
c                                    traction separation curve
c                   intfprps( ,11) - effective displacement jump
c                                    at peak stress
c                   intfprps( ,12) - (beta) a ratio to determine the
c                                    effective displacement jump under
c                                    mixed mode loading
c
c               (1,1) = (2,2) sets isotropic shear-sliding stiff. * cep
c
      do i = 1, span
        cep(i,1,1) = e * intfprps(i,5)/intfprps(i,11) * dj(i) * w *
     &                   intfprps(i,12)**2
        cep(i,2,2) = e * intfprps(i,5)/intfprps(i,11) * dj(i) * w *
     &                   intfprps(i,12)**2
        cep(i,3,3) = e * intfprps(i,5)/intfprps(i,11) * dj(i) * w
      end do
c
      case( 5 )
c     =========
c
c             Option 5:  alternate exponential option -- not implemented
c
      write(iout,9110)
      call die_abort
c
      case( 6 )
c     =========
c
c             Option 6: PPR cohesive model.
c             Ref. JMPS 57 (6), 891-908 (2009)
c
      call lcnst4_ppr( span, cep, intfprps, ppr_support,
     &                 elem_killed, dj, w, iout, mxvl, nstr )
c
c
      case( 7 )
c     =========
c
c             Option 7: cavitation-based cohesive model.
c             Ref. JMPS. 1998;47:99–139
c
      call lcnst4_cavit( span, cep, dtime, intfprps,
     &                   dj, w, iout, mxvl, nstr, gpn )
c
c
      case default
c     ============
c
        write(iout,9530)
        call die_abort
c
      end select
c
      if( local_debug ) then
        write(iout,9530)
        do i = 1, span
          write(iout,9532) i, felem+i-1, cep(i,1,1:3),
     &                     cep(i,2,1:3), cep(i,3,1:3)
        end do
       end if
c
      return
c
 9100 format('>>>> FATAL ERROR: ramp option not implemented',
     & /,    '                  use ppr or exponential option.',
     & /,    '                  job terminated by lcnst4' )
 9200 format('>>>> FATAL ERROR: bilinear option not implemented',
     & /,    '                  use ppr or exponential option.',
     & /,    '                  job terminated by lcnst4' )
 9110 format('>>>> FATAL ERROR: exponential 2 model not implemented',
     & /,    '                  use ppr or exponential option.',
     & /,    '                  job terminated by lcnst4' )
 9300 format(/," ...... calling lcnst4_cavit")
 9310 format(/," ...... returned from lcnst4_cavit")
 9500 format(/," ...... entered lcnst4 .....",
     1 /,10x,"step, cohes_type, span, mxvl, nstr: ",i6,i3,i4,i4,i3,
     2 /,10x,"gpn, felem, time_n, dtime: ",i3,i8, e14.6,1x,e14.6,
     3 /,10x,"weight(w): ",f9.6)
 9510 format(/,10x,"temperature data:",
     1 /,9x,"i  ",2x,"element",5x,"ref temp",8x,"temp end of step" )
 9512 format(7x,i3,i8,3f15.4)
 9520 format(/,10x,"Jacobian determinant:",
     1 /,9x,"i  ",2x,"element",10x,"dj" )
 9522 format(7x,i3,i8,f15.4)
 9530 format(/,10x,"[D] matrices:",
     1 /,9x,"i  ",2x,"element",15x,"----- [cep] -----")
 9532 format(7x,i3,i8,2x,3f15.4,/20x,3f15.4,/,20x,3f15.4)
c
      end
c
c    ****************************************************************
c    *                                                              *
c    *               subroutine lcnst4_ppr                          *
c    *                                                              *
c    *            written by : Kyoungsoo Park                       *
c    *                                                              *
c    *     computes the linear [D]                                  *
c    *     for a  block of similar non conflicting cohesive         *
c    *     elements for the unified potential-based model (PPR)     *
c    *     for the first iteration at the first step                *
c    *     (deln = 0, delt = 0)                                     *
c    *     Ref. JMPS 57 (6), 891-908 (2009)                         *
c    *                                                              *
c    ****************************************************************
c
      subroutine lcnst4_ppr( span, cep, intfprps, ppr_support,
     &                       elem_killed, dj, w, iout, mxvl, nstr )
      implicit none
c
      integer :: i, span, mxvl, nstr, iout
#dbl      double precision
#sgl      real
     &         cep(mxvl,nstr,*), intfprps(mxvl,*), dj(*), w,
     &         ppr_support(mxvl,20)
      logical elem_killed(*)
c
#dbl      double precision
#sgl      real
     &  Dnn, Dtt, Gam_n, Gam_t, dn, dt, m, n, alph, beta, dGtn, dGnt,
     &  zero, one, two
      data zero,one,two / 0.0d00, 1.0d00, 2.0d00 /
c
      call mm04_init_ppr( span, intfprps, ppr_support,
     &                    elem_killed, iout, mxvl )
c
c            the shear-sliding response is set to be isotropic
c            (1,1) = (2,2). cep zeroed by caller. note that the linear
c            [D] is diagonal.
c
      do i = 1, span
        alph  = ppr_support(i,6)
        beta  = ppr_support(i,7)
        m     = ppr_support(i,10)
        n     = ppr_support(i,11)
        dn    = ppr_support(i,12)
        dt    = ppr_support(i,13)
        Gam_n = ppr_support(i,14)
        Gam_t = ppr_support(i,15)
        dGnt  = ppr_support(i,16)
        dGtn  = ppr_support(i,17)
        Dnn = -Gam_n/dn**2*(m/alph)**(m-one)*(m+alph)*
     &         (Gam_t*(n/beta)**n+dGtn)
        Dtt = -Gam_t/dt**2*(n/beta)**(n-one)*(n+beta)*
     &        (Gam_n*(m/alph)**m+dGnt)
        cep(i,1,1) = Dtt * dj(i) * w
        cep(i,2,2) = Dtt * dj(i) * w
        cep(i,3,3) = Dnn * dj(i) * w
      end do
c
      return
      end
c
c    ****************************************************************
c    *                                                              *
c    *               subroutine lcnst4_cavit                        *
c    *                                                              *
c    *                    written by cj, xz                         *
c    *                    last modified: 10/23/2014 rhd             *
c    *                                                              *
c    *     computes the linear [D] (diagonal)                       *
c    *     for a  block of interface-cohesive                       *
c    *     elements for the cavitation-based model (cavit)          *
c    *     for the first iteration at the first step                *
c    *     Ref. JMPS. 1998, Vol. 47, pp 99-139                      *
c    *                                                              *
c    ****************************************************************
c

      subroutine lcnst4_cavit( span, cep, dtime, intfprps, 
     &                         dj, w, iout, mxvl, nstr, gpn )
      implicit none
c
c             parameters
c
      integer :: span, mxvl, nstr, iout, gpn
#dbl      double precision ::
#sgl      real ::
     & cep(mxvl,nstr,*), dj(*), w, dtime, intfprps(mxvl,*)
c     
c             locals - includes the defined type for cavity properties
c             that simplifies managing a large number of scalar values.
c
c             * change all other occurrences in 
c             * warp3d of the type def (mm04.f)
c
      integer :: i
      logical :: here_debug
#dbl      double precision ::
#sgl      real ::
     & three, one, half, four, stiff_shear, stiff_normal,
     & f_0, q  
c
      type :: props_for_cavit  ! multiple instances here, cnst4.f
        double precision
     &    eta_b, Sigma_0, D, a_0, b_0, V_0, N_I, F_N,     
     &    n_pow, N_max, psi_angle_radians, psi_angle_degrees, hpsi,
     &    compression_mult, Sthr, beta, v2dot_term2
      end type
      type( props_for_cavit) :: props
c
      data one, three, half, four
     & / 1.0d0, 3.0d0, 0.5d0, 4.0d0 /   
c
c            get properties for the cavity cohesive material.
c            all interface elements in block have same properties.
c            use derivatives evaluated for zero stress and zero
c            rates but with \Delta_t at t=0 known.
c
      here_debug = .false.
      call mm04_cavit_set_props( intfprps, mxvl, props, iout,
     &                           here_debug  )
      if( here_debug ) then
          write(iout,*) "... props%a_0: ", props%a_0 
          write(iout,*) "... props%b_0: ", props%b_0 
          write(iout,*) "... props%D: ", props%D
          write(iout,*) "... dtime: ", dtime 
      end if    
c      
      f_0 = (props%a_0 / props%b_0)**2
      q   = log(one/f_0) - half * (three-f_0) * (one-f_0)
      stiff_normal = props%b_0**2 * q / ( four * props%D * dtime )
      stiff_shear  = props%eta_b / dtime
c
      if( here_debug ) then
          write(iout,*) "... f_0: ", f_0
          write(iout,*) "... q: ", q
          write(iout,*) "... stiff_normal: ", stiff_normal
          write(iout,*) "... stiff_shear: ",stiff_shear 
      end if    
c
      do i = 1, span
         cep(i,1,1) = stiff_shear * dj(i) * w
         cep(i,2,2) = stiff_shear * dj(i) * w
         cep(i,3,3) = stiff_normal * dj(i) * w
      end do  
C      
      return
      end
