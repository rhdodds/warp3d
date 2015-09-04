c
c *******************************************************************
c *                                                                 *
c *        material model # 4 --  cohesive material model           *
c *                                                                 *
c *            written by  : aroy                                   *
c *       last modified by : cj,  cavitation-based cohesive element *
c *                               (referred to as Cavit) 8/1/2013   *
c *                                                                 *
c *       10/16/2014 rhd (convert cavit to use physical units)      *
c *        3/1/2015 rhd updated for new formulation that computes   *
c *                      NI and props set relative to NI not NR     *
c *                      remove R_I as input                        *
c *        3/5/2015 rhd add 2 more history. Tshear and max ever Ts  *
c *                                                                 *
c *******************************************************************
c
      subroutine mm04(
     1 step, iter, span, felem, gpn, iout, mxvl, time_n, dtime,
     2 nonlocal, numthreads, now_thread, intfprps, c_type,
     3 trac_n, trac_n1, reladis, delrlds, history,
     4 history1, temp_ref, dtemp, temp_n,
     5 top_surf_elems,
     6 bott_surf_elems, top_surf_stresses_n, bott_surf_stresses_n,
     7 top_surf_eps_n, bott_surf_eps_n, top_nonlocal_vars,
     8 bott_nonlocal_vars, top_solid_matl,
     9 bott_solid_matl )
c
      implicit none
c
c                   parameter declarations. arrays have sizes
c                   dimensioned to the maximum allowable elements
c                   per block during execution (mxvl) typically
c                   set to 128 or 256.
c
      integer step, iter, span, felem, gpn, iout, mxvl,
     1        numthreads, now_thread, c_type
c
#dbl      double precision
#sgl      real
     1 intfprps(mxvl,*), trac_n(mxvl,*), delrlds(mxvl,*),
     2 trac_n1(mxvl,*), reladis(mxvl,*),
     3 history(span,13), history1(span,13),
     4 time_n, dtime,
     5 top_surf_stresses_n(mxvl,6), bott_surf_stresses_n(mxvl,6),
     6 top_surf_eps_n(mxvl,6), bott_surf_eps_n(mxvl,6),
     7 temp_ref(span), dtemp(span), temp_n(span),
     8 top_nonlocal_vars(mxvl,*), bott_nonlocal_vars(mxvl,*)
      integer top_surf_elems(span), bott_surf_elems(span),
     1        top_solid_matl(span), bott_solid_matl(span)
      logical nonlocal
c
c                stress updating for cohesive material models.
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
c                    iter = 0 used to find stresses for imposed
c                    displacements and temperatures. History
c                    need not be updated for iter = 0.
c                    iter = 1, 2, 3, .... standard equilibrium
c                    iterations at specified, fixed loading.
c                    history1 should be updated.
c
c (inout)      span: number of coehsive elements in this block
c
c (input)     felem: number of first element (global numbering) in
c                    this block. elements in model are numbered
c                    sequentially.
c
c (input)      gpn:  integration point number for this block of
c                    cohesive elements
c (input)      iout: Fortran device number for output of any messages
c                    from the material model
c
c (input)      mxvl: maximum allowable number of elements per block
c                    (dimensioned number of rows for many matrices)
c
c (input)    time_n: simulation time at start of load increment
c                    (sum of all delta times)
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
c (input) intfprps:  material properties and computational options
c                    at this integration point for all elements in the
c                    block. see mm04_init for detailed description of
c                    each column. most often the properties are
c                    identical for all elements in the block but
c                    this is not required.
c
c (input)   c_type:  the model supports various options for the cohesive
c                    formulation
c                      1 - linear elastic
c                      2 - bilinear -- not impelemented
c                      3 - ramp     -- not impelemented
c                      4 - exponential_1
c                      5 - exponential_2 -- not impelemented
c                      6 - PPR (Park, Paulino, Roesler)
c                      7 - Cavitation-based cohesive element, referred
c                          to as Cavit
c                          Reference: Onck P, Van DerGiessen E.
c                          Growth of an initially sharp crack by grain
c                          boundary
c                          J. Mech. Phys. Solids 1998; 47, 99-139
c
c
c (input)    trac_n: cohesive tractions at time step n (input)
c (update)  trac_n1: cohesive tractions at time step n+1 (output)
c                     ( ,1):  cohesive traction t1 (shear 1)
c                     ( ,2):  cohesive traction t2 (shear 2)
c                     ( ,3):  cohesive traction tn (normal)
c                     ( ,4):  - not used by cohesive: set = 0
c                     ( ,5):  - not used by cohesive: set = 0
c                     ( ,6):  - not used by cohesive: set = 0
c                     ( ,7):  total cohesive energy at step n+1
c                     ( ,8):  unrecoverable cohesive energy at step n+1
c                    - t1 and t2 are orthogonal and in the tangent plane
c                    - positions 4, 5, 6 not used by cohesive models
c
c (input)   reladis: current estimate of displacement jumps between
c                    n=0 -> n+1 (ds1, ds2, dn). used for linear
c                    elastic and secant formulations
c                     ( ,1):  ds1 (shear sliding 1)
c                     ( ,2):  ds2 (shear sliding 2)
c                     ( ,3):  dn  (normal)
c                    ds1 and ds2 are orthogonal and in the tangent plane
c
c (input)   delrlds: current estimate of displacement jumps
c                    between n -> n+1
c                     ( ,1):  delta-ds1
c                     ( ,2):  delta-ds2
c                     ( ,3):  delta-dn
c                    delta-ds1 and delta-ds2 are orthogonal and in the
c                    tangent plane
c
c (input)   history: state variables at step n. The number of values per
c                    integration point per element was set in subroutine
c                    mm04_set_sizes.
c
c (update) history1: state variables at step n+1 (update these)
c
c               history content for exponential
c               model (c_type = 4 ) data for state 'n' and 'n+1'
c                   1 -- effective displacement
c                   2 -- maximum effective displacement
c                   3 -- effective traction
c                   4 to last -- not used
c
c               history content for PPR (c_type = 6 )
c               data for state 'n' and 'n+1'
c                   1 -- dn_limit: normal separation when
c                        normal traction falls to zero
c                   2 -- dt_limit: tangential separation when
c                        shear traction falls to zero
c                   3 -- deln_max: maximum normal separation
c                        during the loading history
c                   4 -- delt_max: maximum tangential separation
c                        during the loading history
c                   5 to last -- not used
c
c               history content for Cavit (c_type = 7 ) data
c               see comments in mm04.f. some values are normalized as
c               indicate. 
c               updated Nov 20, 2014 rhd
c                   1 -- N / N_I cavity density 
c                   2 -- a / a_0 cavity radius
c                   3 -- b / b_0 (2b) is distance between centers of 
c                        adjacent cavities) cavities.
c                   4 -- V / V_0 cavity volume 
c                   5 -- Tn traction normal to GB
c                   6 -- traction shear stiffness
c                   7 -- tractions normal stiffness
c                   8 -- critical state value. 0.0 or 1.0. not yet used.
c                   9 -- maximum normal traction (Tn) experienced over 
c                        loading history (not normalized)
c                  10 -- opening displacement value (delta_c)
c                        corresponding to above traction (not 
c                        normalized) 
c                  11 -- a_0 / b_0. used in to computing
c                        actual a / b ratio at output time 
c                        store in history just for convenience at 
c                        output time
c                  12 -- a/ Lnr. current cavity radius normalized
c                        by Needleman-Rice characteristic length
c                  13 -- material state for driving state table to 
c                        update stresses (= 0, 1, 2, 3). pull integer
c                        from left part of double precision word
c                  14 -- Ts shear traction on GB
c                  15 -- maximum shear traction (Ts) experienced over 
c                        loading history (not normalized)
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
c nonlocal = .true.
c
c (input)  top_surf_elems: number of solid element connected to top
c                          surface for each element in  block
c
c (input) bott_surf_elems: number of solid element connected to bottom
c                          surface for each element in  block
c
c (input) top_surf_stresses_n: stresses for solid element attached to
c                          top surface of each element in block. 6
c                          components in model (global) coordinates.
c                          xx, yy, zz, zy, yz, xz. These are averaged
c                          centroidal values for solids.
c
c (input) bott_surf_stresses_n: stresses for solid element attached to
c                          bottom surface of each element in block. 6
c                          components in model (global) coordinates.
c                          xx, yy, zz, zy, yz, xz. These are averaged
c                          centroidal values for solids.
c
c (input) top_surf_eps_n:  strains for solid element attached to
c                          top surface of each element in block. 6
c                          components in model (global) coordinates.
c                          xx, yy, zz, zy, yz, xz. These are averaged
c                          centroidal values for solids. shears are
c                          gamma not epsilon
c
c (input) bott_surf_eps_n: strains for solid element attached to
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
c (input) top_solid_matl: material model type in WARP3D system for
c                          surface for each solid element in a
c                          attached to top surface
c
c (input) bott_solid_matl: material model type in WARP3D system for
c                          surface for each solid element in a
c                          attached to bottom surface
c  Note: when a surface of a cohesive element is attached to a symmetry
c        plane that solid element = 0 above. then stresses have been
c        made equal for top and bottom solids. same for strains.
c
c
c      for a large displacement analysis, the area change of the
c      cohesive surface has already been computed by the corresponding
c      interface-cohesive element. this material model doesn't
c      see the difference between a small and large displacement
c      analysis - exactly as other material models in WARP3D
c
c
c ------------------------------------------------------------------
c
c                   locally defined
c
      integer i, history_length, info_vector(10)
#dbl      double precision
#sgl      real
     & half, zero, one, three, e, tol, forty, initial,
     & beta, ds1, ds2, dn, tn_n1, dn_n1, tshear,
     & intfmat(mxvl,3), ppr_support(mxvl,30),
     & top_surf_mean_stress, top_surf_mean_eps,
     & bott_surf_mean_stress, bott_surf_mean_eps, comp_multiplier
       logical elem_killed(mxvl), local_debug
       data zero, e, tol, one, three, forty, initial
     &  / 0.0d0, 2.71828182845904523536d0, 0.1d0, 1.0d0, 3.0d0,
     &    40.0d0, 1.d-14 /
      data half / 0.5d00 /
c
c ------------------------------------------------------------------
c
      local_debug = intfprps(1,29) .gt. zero
      local_debug = .false.
c 
      call mm04_set_sizes( info_vector ) 
      history_length = info_vector(1)  !  note this so we avoid bugs
c                                         if history changes size
c
c             debug output on entry to mm04
c
      if( local_debug ) then
        write(iout,9500) step, iter, span, felem, gpn, mxvl,
     1                   time_n, dtime, nonlocal, numthreads,
     2                   now_thread, c_type
        write(iout,9510)
        do i = 1, span
          write(iout,9512) i, felem+i-1, temp_ref(i), dtemp(i),
     1                     temp_n(i)
        end do
      end if
c
      if( nonlocal .and. local_debug ) then
          write(iout,9520)
          do i = 1, span
            top_surf_mean_stress = (top_surf_stresses_n(i,1) +
     &                  top_surf_stresses_n(i,2) +
     &                  top_surf_stresses_n(i,3)) / three
            bott_surf_mean_stress = (bott_surf_stresses_n(i,1) +
     &                  bott_surf_stresses_n(i,2) +
     &                  bott_surf_stresses_n(i,3)) / three
            top_surf_mean_eps    = (top_surf_eps_n(i,1) +
     &                  top_surf_eps_n(i,2) +
     &                  top_surf_eps_n(i,3)) / three
            bott_surf_mean_eps   = (bott_surf_eps_n(i,1) +
     &                  bott_surf_eps_n(i,2) +
     &                  bott_surf_eps_n(i,3)) / three
         write(iout,9522) i, felem+i-1, top_surf_elems(i),
     &      bott_surf_elems(i),  top_surf_mean_stress,
     &      top_surf_mean_eps,
     &       bott_surf_mean_stress, bott_surf_mean_eps
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
c             for step 1, init the history data. for each iteration of
c             step 1, we're getting the solution from n = 0 to the
c             current iteration of step 1. history1 zeroed before call.
c             some cohesive sub-types subsequently do their own
c             initialization.
c
      if( step .eq. 1 ) history(1:span,1:history_length) = zero
c
c             set flags for already killed elements (exhausted
c             the cohesive traction). zero the unused stress
c             locations (for cohesive materials)
c
      do i = 1, span
        trac_n1(i,4:6) = zero
        elem_killed(i) = intfprps(i,13) .gt. zero
      end do
c
c             branch on cohesive formulation
c
      select case ( c_type )
c     ----------------------
c
      case ( 1 )
c     ==========
c
c             linear traction - separation law is used in this
c             cohesive zone model. normal and shear deformations are
c             uncoupled. user sets the two shear stiffness values equal
c             to have isotropic shear response. for cohesive elements
c             on symmetry planes. user must set stiffness values at twice
c             values used in a full (w/o symmetry). see WARP3D
c             manual for the cohesive material model.
c
      do i = 1, span
        ds1 = reladis(i,1)
        ds2 = reladis(i,2)
        dn  = reladis(i,3)
        trac_n1(i,1) = intfprps(i,2) * ds1
        trac_n1(i,2) = intfprps(i,3) * ds2
        trac_n1(i,3) = intfprps(i,4) * dn
        if( dn .lt. zero ) then
           comp_multiplier = intfprps(i,22)
           trac_n1(i,3) = trac_n1(i,3) * comp_multiplier 
        end if
        trac_n1(i,7) = half * ( trac_n1(i,1)*ds1 +
     &                 trac_n1(i,2)*ds2 + trac_n1(i,3)*dn )
        trac_n1(i,8) = zero
        trac_n1(i,9) = zero
      end do
c
      if( local_debug ) then
        write(iout,9300) step, iter
        do i = 1, span
           write(iout,9310) i, reladis(i,1:3), trac_n1(i,1:3)
        end do
      end if
c
      case ( 2 )   !   bilinear not implemented
c     ==========
c
      write(iout,9200)
      call die_abort
c
      case ( 3 )   !   ramp not implemented
c
      write(iout,9100)
      call die_abort
c
      case( 4 )    ! standard, 2 parameter exponential option
c     =========
c
c             * the current estimate of total displacement jumps
c             at n+1. these secant terms are found by evaluating
c             traction components at n+1 using tractions separation curve
c             then dividing by current total displacement jump. this
c             could be re-designed to simplify w/o introducing the
c             secant stiffness.
c
       call mm04_exp1_secant( span, felem, intfprps, reladis,
     &                        history,  history1, intfmat,
     &                        elem_killed, mxvl )
c
c             diagonal secant stiffness avaiable.  use it to 
c             get new total tractions.
c             trac_n1( ,1:3): cohesive tractions at step n+1
c             trac_n1( ,7):   total cohesive energy at step n+1
c             trac_n1( ,8):   unrecoverable cohesive energy at step n+1
c
      do i = 1, span
          trac_n1(i,1) = intfmat(i,1)*reladis(i,1)
          trac_n1(i,2) = intfmat(i,2)*reladis(i,2)
          trac_n1(i,3) = intfmat(i,3)*reladis(i,3)
          trac_n1(i,7) = trac_n(i,7) +
     &       half * ( (trac_n1(i,1) + trac_n(i,1))*delrlds(i,1) +
     &              (trac_n1(i,2) + trac_n(i,2))*delrlds(i,2) +
     &              (trac_n1(i,3) + trac_n(i,3))*delrlds(i,3) )
          trac_n1(i,8) = trac_n1(i,7) - half *
     &                   ( trac_n1(i,1)*reladis(i,1) +
     &                     trac_n1(i,2)*reladis(i,2) +
     &                     trac_n1(i,3)*reladis(i,3) )
          trac_n1(i,9) = history1(i,2)
      end do
c
      if( local_debug ) then
        write(iout,9400) step, iter
        do i = 1, span
           write(iout,9310) i, reladis(i,1:3), trac_n1(i,1:3)
        end do
      end if
c
      case ( 5 ) ! alternate exponential - not implemented
      write(iout,9110)
      call die_abort
c
c             Paulino-Park-Roesler (PPR) option
c             1) user input material props have already been loaded into
c             the intfprps by WARP3D drivers
c             2) compute a local block of parameters from material
c             properties for PPR computations.
c             Note: each element in block must be PPR, but each one may
c                 have different user defined material properties.
c             3) compute cohesive traction
c
      case ( 6 )  ! ppr
c     =========
      call mm04_init_ppr( span, intfprps, ppr_support,
     &                    elem_killed, iout, mxvl )
      call mm04_traction_ppr( span, ppr_support, reladis, trac_n1,
     &             history, history1, elem_killed,
     &             local_debug, iout, mxvl )
      do i = 1, span
        trac_n1(i,7) = trac_n(i,7) +
     &       half * ( (trac_n1(i,1) + trac_n(i,1))*delrlds(i,1) +
     &             (trac_n1(i,2) + trac_n(i,2))*delrlds(i,2) +
     &             (trac_n1(i,3) + trac_n(i,3))*delrlds(i,3) )
        trac_n1(i,8) = trac_n1(i,7) - half*
     &                   ( trac_n1(i,1)*reladis(i,1) +
     &                     trac_n1(i,2)*reladis(i,2) +
     &                     trac_n1(i,3)*reladis(i,3) )
        trac_n1(i,9) = zero
      end do
c
c
c             Cavitation-based cohesive zone model option (Cavit)
c             user input material props have already been loaded into
c             the intfprps by WARP3D drivers
c
      case ( 7 )  ! cavit
c     =========
c
      if( .not. nonlocal ) then  ! only viable for nonlocal formulation
          write(iout,9600)
          call die_abort   
      end if
      call mm04_cavit_driver   !  within contains for simplicity
c
      case default
c     ============
c
      write(iout,9120)
      call die_abort
c
      end select
c
      return
 1111 format(6(F16.8))
 9100 format('>>>> FATAL ERROR: ramp model not implemented',
     & /,    '                  use exponential model.',
     & /,    '                  job terminated by mm04' )
 9120 format('>>>> FATAL ERROR: unknown cohesive formulation',
     & /,    '                  job terminated by mm04' )
 9200 format('>>>> FATAL ERROR: bilinear model not implemented',
     & /,    '                  use exponential model.',
     & /,    '                  job terminated by mm04' )
 9110 format('>>>> FATAL ERROR: exponential 2 model not implemented',
     & /,    '                  use exponential model.',
     & /,    '                  job terminated by mm04' )
 9300 format('... compute linear-elastic tractions:',
     &    /, '      step, iter: ',2i6)
 9310 format('... updated tractions, element in block: ', i4, /,
     &    10x, ' del(1-3):  ',3f15.10,/
     &    10x, ' tract(1-3):',3f15.4 )
 9400 format('... compute exp1_intf tractions:',
     &    /, '      step, iter: ',2i6)
 9500 format(/," ...... entered mm04 .....",
     1 /,10x,"step, iter, span, felem: ",i6,i3,i4,i7,
     2 /,10x,"gpn, mxvl, time_n, dtime: ",i2,i4,2e14.6,
     3 /,10x,"nonlocal, numthreads, nowthread, c_type: ",l2,i4,i4,i4 )
 9510 format(/,10x,"temperature data:",
     1 /,9x,"i  ",2x,"element",5x,"ref temp",2x,"dtemp over step",5x,
     2 "temp start of step" )
 9512 format(7x,i3,i8,2f15.4,5x,f15.4)
 9520 format(/,10x,"nonlocal data ( @ start of step ):",
     1 /,9x,"i  ",2x,"element",5x,"top solid",8x,"bottom solid",
     2  3x,"top sig mean",5x,"top eps mean",
     3  8x,"bott sig mean",5x,"bott eps mean" )
 9522 format(7x,i3,2x,i8,5x,i10,8x,i10,2x,2f15.6,7x,2f15.6)
 9600 format('>>>> FATAL ERROR: cavitation cohesive model',
     & /,    '                  requires the nonlocal option.',
     & /,    '                  job terminated by mm04' )
 9625 format(/,10x,"more nonlocal data (material type, solid state",
     a " variables passed thru @ start of step):",
     b /,9x,"i  ",2x,"element", 2x,"top matl type",5x,
     c   "bottom matl type" )
 9627 format(7x,i3,2x,i8,2x,i8,5x,i10)
 9629 format(12x,"top solid vars @ n:    ",5f10.5,
     a /,12x,    "bottom solid vars @ n: ",5f10.5 )
c
      
      contains   !     ....  note this  .... 
c     ========
c
c    ****************************************************************
c    *                                                              *
c    *          subroutine mm04_cavit_driver                        *
c    *                                                              *
c    *            written by : rhd   6/9/2014                       *
c    *            updated: 3/4/2014  rhd                            *	
c    *                                                              *
c    *     drive stress update of cavity cohesive model             *  
c    *                                                              *
c    ****************************************************************
c
      subroutine mm04_cavit_driver
      implicit none
c
c             variables declared here are local to this routine.
c      
      logical :: here_debug
c
#dbl      double precision ::
#sgl      real ::
     & l_nr, a_now, third, one, dword
      integer :: iword(2)
      equivalence (dword, iword)
c 
      type :: props_for_cavit  ! multiple instances here, cnst4.f
        double precision
     &    eta_b, Sigma_0, D, a_0, b_0, V_0, N_I, F_N, 
     &    n_pow, N_max, psi_angle_radians, psi_angle_degrees, hpsi,
     &    compression_mult, Sthr, beta, v2dot_term2
      end type
      type( props_for_cavit) :: props
c
      data one, third / 1.0d00, 0.33333333d00 /      
c  
c             see mm04_init for intfprps definition
c
      here_debug = felem .eq. 55 .and. gpn .eq. 1
      here_debug = .false.
c      
      if( here_debug ) write(iout,9000)
      call mm04_cavit_set_props( intfprps, mxvl, props, iout, 
     &                           .false. )
      if( here_debug ) write(iout,9005)
c      
c             see top of mm04 for history definition.
c             normalized state variables stored for convenience
c             of including those values in output
c
      if(  step .eq. 1 ) then
         do i = 1, span
           history(i,1)     = one  !  N / props%N_I
           history(i,2)     = one  !  a / props%a_0
           history(i,3)     = one  !  b / props%b_0
           history(i,4)     = one  !  V / props%V_0
           history(i,5:10)  = zero
           history(i,11)    = props%a_0 / props%b_0
           history(i,12)    = zero
           iword(1) = 0; iword(2) = 0
           history(i,13)    = dword
           history(i,14:15)    = zero
           history1(i,1:15) = history(i,1:15)
           trac_n(i,1:3)    = zero
         end do  
         if( here_debug ) write(iout,9010)
      end if    
c
c             material properties are identical for all
c             elements in block.
c
c             use the simple typedef above rather than making
c             a set of values for each element in block.
c          
      if( here_debug ) write(iout,9015)
      call mm04_cavit_sig_update( 
     &     step, iter, span, felem, gpn, iout, mxvl, dtime, nonlocal,
     &     props, intfprps,
     &     trac_n, trac_n1, reladis, delrlds, history, history1,
     &     top_surf_stresses_n, bott_surf_stresses_n,
     &     top_nonlocal_vars, bott_nonlocal_vars, 
     &     local_debug, history_length )
      if( here_debug ) write(iout,9020)
c      
c             updating completed. compute energy quantities and
c             other values in history used later in output.
c   
c             trac(*,7) = total work of tractions integrated over
c                         loading history (trapezoidal rule)
c
c             trac(*,8) = total work of tractions integrated over
c                         loading history (trapezoidal rule) -
c                         triangular area under a secant to origin
c                         unloading response
c
c             retain max norm & shear tractions ever reached 
c             and the opening displacement at max normal
c             traction for output and possible element death operations.
c                          
      do i = 1, span
        trac_n1(i,7) = trac_n(i,7) +
     &       half * ( (trac_n1(i,1) + trac_n(i,1))*delrlds(i,1) +
     &             (trac_n1(i,2) + trac_n(i,2))*delrlds(i,2) +
     &             (trac_n1(i,3) + trac_n(i,3))*delrlds(i,3) )
        trac_n1(i,8) = trac_n1(i,7) - half*
     &                   ( trac_n1(i,1)*reladis(i,1) +
     &                     trac_n1(i,2)*reladis(i,2) +
     &                     trac_n1(i,3)*reladis(i,3) )
        trac_n1(i,9)   = zero 
c
        history1(i,8) = zero ! reserved for a damage parameter
c        
        tn_n1          = trac_n1(i,3)
        dn_n1          = reladis(i,3)
        history1(i,9)  = history(i,9)
        history1(i,10) = history(i,10)
        if( tn_n1 .gt. history(i,9) ) then
           history1(i,9)  = tn_n1  ! see note 1
           history1(i,10) = dn_n1 ! length units in model
        end if
        tshear = history1(i,14)
        if( tshear .gt. history1(i,15) ) history1(i,15)  = tshear
        if( here_debug ) then
           write(iout,*) "...tn_n1, dn_n1: ",tn_n1, dn_n1
        write(iout,*) "...hist1(9,10): ",history1(i,9),history1(i,10)
        end if
c
      end do   !  over i for energy update, extra output values
c 
      if( here_debug ) write(iout,9050)
      return
      
 9000 format("... entered mm04_cavit_driver")
 9005 format("... props created by mm04_cavit_set_props")
 9010 format("... history initialized")
 9015 format("... calling mm04_cavit_sig_update")
 9020 format("... back from mm04_cavit_sig_update")
 9050 format("... leaving mm04_cavit_driver") 
c 
      end subroutine mm04_cavit_driver
      end subroutine mm04
c
c    ****************************************************************
c    *                                                              *
c    *          subroutine mm04_cavit_set_props                     *
c    *                                                              *
c    *            written by : rhd   6/12/2014                      *
c    *            updated: 3/4/2015 rhd                             *	
c    *                                                              *
c    *     set up properties for cavity cohesive option             *
c    *                                                              *
c    ****************************************************************
c
      subroutine mm04_cavit_set_props( intfprps, mxvl, props, iout,
     &                                 here_debug  )
      implicit none
c
c             parameters
c
      integer :: mxvl, iout
      logical :: here_debug
#dbl      double precision ::
#sgl      real ::
     & intfprps(mxvl,*)
c
c             locals - includes the defined type for cavity properties
c             that simplifies managing a large number of scalar values.
c             if changed, find all other instances in warp3d of this
c             typedef and change them too.
c
#dbl      double precision ::
#sgl      real ::
     & pi, one_eighty, zero, one, three, four, psi, x, half,
     & onept5
c 
      type :: props_for_cavit  ! multiple instances here, cnst4.f
        double precision
     &    eta_b, Sigma_0, D, a_0, b_0, V_0, N_I, F_N,     
     &    n_pow, N_max, psi_angle_radians, psi_angle_degrees, hpsi,
     &    compression_mult, Sthr, beta, v2dot_term2
      end type
      type( props_for_cavit) :: props
c      
      data  one,  pi, one_eighty, half, three, four, onept5
     &     / 1.0d00, 3.141592653589793d00, 180.0d00, 0.5d00, 
     &       3.0d00, 4.0d00, 1.5d00 /
c
c             all interface elements in the block have the same
c             cohesive properties.
c  
c             see mm04_init for intfprps definition
c
c             note changes 3/1/2015:
c
c               - R_I no longer allowed input
c               - a_0, b_0 must be the actual, physical values
c               - N_I now computed from b_0. 
c               - N_I no longer allowed as input property
c               - F_N and N_max multipliers are on N_I, not N_R
c               - N_R no longer used
c              
c
      props%eta_b   = intfprps(1,39)
      props%Sigma_0 = intfprps(1,41)
      props%D       = intfprps(1,43)
      props%a_0     = intfprps(1,44) 
      props%b_0     = intfprps(1,45) 
      props%psi_angle_degrees = intfprps(1,46)
      props%psi_angle_radians = intfprps(1,46) * pi / one_eighty
      x = props%psi_angle_radians ! makes next stm simpler to read
      props%hpsi    = ( one/(one+cos(x)) - half*cos(x) ) / 
     &                       sin(x)
      props%V_0     = (four/three) * pi * props%a_0**3 * props%hpsi 
      props%N_I     = one/ pi / props%b_0**2
      props%F_N     = intfprps(1,48) * props%N_I
      props%Sthr    = props%N_I / props%F_N
      props%n_pow   = intfprps(1,49)
      props%N_max   = intfprps(1,50) * props%N_I
      props%compression_mult = intfprps(1,22)
      props%beta    = (props%n_pow - one) * (props%n_pow + 0.4319d00) / 
     &                 props%n_pow**2
      props%v2dot_term2 = onept5 / props%n_pow
c
      if( here_debug ) then 
         write(iout,9000) props%eta_b, 
     &                    props%Sigma_0, props%D, props%a_0, 
     &                    props%b_0, props%psi_angle_degrees, 
     &                    props%psi_angle_radians, props%hpsi, 
     &                    props%V_0, props%N_I, props%F_N,   
     &                    props%n_pow, props%N_max, 
     &                    props%compression_mult
      end if
c    
      return
c
 9000 format("...  mm04_cavit_set_props ...",
     &  /,5x,"eta_b:          ",e14.5,
     &  /,5x,"Sigma_0:        ",f10.1
     &  /,5x,"D:              ",e14.5,
     &  /,5x,"a_0:            ",e14.5,
     &  /,5x,"b_0:            ",e14.5,
     &  /,5x,"psi (degrees):  ",f8.1,
     &  /,5x,"psi (radians):  ",f8.2,
     &  /,5x,"h(psi):         ",f8.3,
     &  /,5x,"V_0:            ",e14.5,
     &  /,5x,"N_I:            ",e14.5,
     &  /,5x,"F_N:            ",e14.5,
     &  /,5x,"n_pow:          ",f5.1,
     &  /,5x,"N_max:          ",e14.5,
     &  /,5x,"comp. factor:   ",f5.1 )
c    
      end       
c
c    ****************************************************************
c    *                                                              *
c    *          subroutine mm04_cavit_sig_update                    *
c    *                                                              *
c    *    updated: 8/28/2015 rhd added switch_to_linear feature     8
c    *                                                              *
c    *     computes the cohesive traction vector for the            *
c    *     cavitation-based cohesive zone model (cavit)             *
c    *     this code requires physical material property            *
c    *     values and returns physical tractions.                   *
c    *     Ref. JMPS 1998; 47, 99-139 (as starting point)           *
c    *                                                              *
c    ****************************************************************
c
      subroutine mm04_cavit_sig_update(
     & step, iter, span, felem, gpn, iout, mxvl, dtime, nonlocal,
     & props, intfprps,
     & trac_n, trac_n1, reladis, delrlds, history, history1,
     & top_surf_stresses_n, bott_surf_stresses_n,
     & top_nonlocal_vars, bott_nonlocal_vars, local_debug,
     & history_len )
c
      implicit none
c
c             parameter declarations
c
      integer :: step, iter, span, iout, mxvl, felem, gpn, history_len
#dbl      double precision ::
#sgl      real ::
     1 trac_n(mxvl,*), delrlds(mxvl,*), intfprps(mxvl,*),
     2 trac_n1(mxvl,*), reladis(mxvl,*),
     3 history(span,history_len), history1(span,history_len), dtime,
     4 top_surf_stresses_n(mxvl,6), bott_surf_stresses_n(mxvl,6),
     5 top_nonlocal_vars(mxvl,*), bott_nonlocal_vars(mxvl,*)
      logical nonlocal
c      
c             find/modify all instances of this typedef
c
      type :: props_for_cavit  ! multiple instances here, cnst4.f
        double precision
     &    eta_b, Sigma_0, D, a_0, b_0, V_0, N_I, F_N,      
     &    n_pow, N_max, psi_angle_radians, psi_angle_degrees, hpsi,
     &    compression_mult, Sthr, beta, v2dot_term2
      end type
      type( props_for_cavit) :: props
c
c             local variables
c
      integer :: i, abs_elem, iter_count, iword(2), old_state,
     &           new_state, current_state
#dbl      double precision ::
#sgl      real :
     & dword, zero, three, one, toler, six, 
     & half, four, two, third, onept5,
     & top_surf_mean_stress, top_surf_mean_eps,
     & bott_surf_mean_stress, bott_surf_mean_eps,
     & sigma_e, sigma_m, term1,term3, fraction,
     & v2_dot, pi, c_0, c_1, c_2, f_0, q,
     & d_strain, c_strain,
     & delta_c_dot, oldreldis, Tn_solid,
     & comp_multiplier, top_sigma_e, bott_sigma_e,
     & top_sigma_m, bott_sigma_m, sig_top(6), sig_bott(6), 
     & N_n, a_n, b_n, V_n, T_n, ratio_1, ratio_2, f_factor,
     & q_factor, Lnr, S, pm, delta_T, T_new, v1_dot, V_new, a_new,
     & N_new, b_new, N_dot, stiff_normal, stiff_shear,
     & max_ab_ratio, stiff_linear, mark, trial_T_new, LNR_toler,
     & tshear_1, tshear_2, f_sd, ab_ratio, ab_sd_ratio
c
      equivalence ( iword, dword )     
c
#dbl      double precision, 
#sgl      real,
     & external :: mm04_cavit_mises
c
      logical :: elem_killed(mxvl), local_debug, converged, 
     & nucleation_active, debug_span_loop, compute_solid_local,
     & iterative_solve, opening, closing, switch_to_linear,
     & debug_newton, open, neutral, penetrated
      logical :: degrade_shear, VVNT, modify_q, include_nucleation,
     &           include_cavity_growth  
c
      data zero, one, three, third, two, half, four,
     &  six, onept5, toler, pi, max_ab_ratio, mark, LNR_toler, 
     &  ab_sd_ratio
     & / 0.0d0, 1.0d0, 3.0d0, 0.3333333333333333d0, 2.0d0, 0.5d0,
     &   4.0d0, 6.0d0, 1.5d00, 1.0d-10, 3.14159265d0, 0.9999d00,
     &   1.0d40, 1.0d-06, 0.5d0 /
c
c             set options in formulation to include, debug options
c
      degrade_shear         = .false.
      VVNT                  = .true.
      modify_q              = .true.
      include_nucleation    = .false.
      debug_newton          = .false.
      include_cavity_growth = .true.
      compute_solid_local   = .false.
c      
      local_debug = felem .eq. 55 .and. gpn .eq. 1 .and. 
     &                  ( step .eq. 599 )
      local_debug = .false.
      if( local_debug ) write(iout,9100) 
c
      do i = 1, span  ! main loop over all elems in block
c
      abs_elem        = felem + i - 1
      debug_span_loop = .false. ! abs_elem .eq. 3182
c      
      if( debug_span_loop ) write(iout,9000) felem+i-1, gpn      
c
c             step 1: pull the cavity state variables from history
c                     at time t_n. back to un-normalized values
c
      N_n = history(i,1) * props%N_I !cavity density
      a_n = history(i,2) * props%a_0 !cavity radius
      b_n = history(i,3) * props%b_0 !half-spacing between cavities
      V_n = history(i,4) * props%V_0 !cavity volume
      T_n = history(i,5) !prior normal traction stress
c
c             step 2: obtain nonlocal variables passed from creep
c                     material model (umat at this time)
c                     
      top_sigma_m  = ( top_surf_stresses_n(i,1) +
     &                 top_surf_stresses_n(i,2) +
     &                 top_surf_stresses_n(i,3) ) / three
      bott_sigma_m = ( bott_surf_stresses_n(i,1) +
     &                 bott_surf_stresses_n(i,2) +
     &                 bott_surf_stresses_n(i,3) ) / three
      sigma_m = half * ( top_sigma_m + bott_sigma_m )
c
      sig_top(1:6)    = top_surf_stresses_n(i,1:6)
      top_sigma_e     = mm04_cavit_mises( sig_top )
      sig_bott(1:6)    = bott_surf_stresses_n(i,1:6)
      bott_sigma_e = mm04_cavit_mises( sig_bott )
      sigma_e      = half * ( top_sigma_e + bott_sigma_e )
      if( compute_solid_local ) call mm04_compute_solid_local( i )
c
      d_strain = half * ( top_nonlocal_vars(i,1) +  ! creep strain rate
     &                    bott_nonlocal_vars(i,1) )
      c_strain = half * ( top_nonlocal_vars(i,2) +  ! creep strain
     &                    bott_nonlocal_vars(i,2) )      
      if( debug_span_loop ) then
           write(iout,9200) N_n, a_n, b_n, V_n, T_n,
     &                   sigma_m, sigma_e, d_strain, c_strain,
     &                   reladis(i,3), delrlds(i,3)
      end if  
c          
c             step 3: compute normal traction at n+1. update steps
c                     based on current state of interface. 
c
      opening     = delrlds(i,3) .ge. zero
      closing     = .not. opening
      open        = reladis(i,3) .gt. zero
      neutral     = reladis(i,3) .eq. zero
      penetrated  = reladis(i,3) .lt. zero
c
c             set current state before entering update.
c
c               state
c                 1         load(time) step 1
c                 2         neutral and closing
c                 3         neutral and opening
c                 4         penetrated and further closing
c                 5         penetrated but re-opening
c                 6         open. either opening/closing.
c                           has a corner subcase of penetrated at n
c                           now open at n+1
c                 7         no cavity growth allowed
c                -1         logic failed to set a valid current state
c
      current_state = -1
      if( step .eq. 1 ) then
          current_state = 1
      elseif( neutral )   then
          if( closing ) current_state = 2
          if( opening ) current_state = 3            
      elseif( penetrated  ) then
          if( closing ) current_state = 4
          if( opening ) current_state = 5            
      elseif( open  ) then
          current_state = 6
          if( opening .and. open ) then
             oldreldis = reladis(i,3) - delrlds(i,3)
             if( oldreldis .lt. zero ) current_state = 8
          end if   
      end if    
      if( .not. include_cavity_growth ) current_state = 7
c
      stiff_normal = mark ! for error checking
      T_new        = mark ! for error checking
      new_state    = -1   ! saved in history for possible future use
c      
      select case( current_state )
      case( 1 ) 
          if( closing ) then
            new_state = 1
            call mm04_cavit_update_linear
          else
            new_state = 2
            call mm04_cavit_std_update
          end if   
      case( 2 ) 
          new_state = 3
          call mm04_cavit_update_linear
      case( 3 ) 
          new_state = 4
          call mm04_cavit_std_update
      case( 4 ) 
          new_state = 5
          call mm04_cavit_update_linear
      case( 5 )
          new_state = 6
          call  mm04_cavit_update_linear
      case( 6 ) 
          new_state = 7
          call mm04_cavit_std_update 
          if( switch_to_linear ) then
            call mm04_cavit_update_linear
            new_state = 10
          end if
      case( 7 ) 
          new_state = 8
          call  mm04_cavit_update_linear
      case( 8 ) 
          new_state = 9
          call  mm04_cavit_update_linear
      case( -1 ) ! no current state set
          write(iout,9300) felem+i-1, gpn
          write(iout,9310) reladis(i,3), delrlds(i,3),
     &                     opening, closing, neutral, penetrated      
          write(iout,9320) N_n, a_n, b_n, V_n, T_n,
     &                  sigma_m, sigma_e, d_strain, c_strain
          call die_abort
      case default ! really weird if gets here
         write(iout,*) '*** fatal error mm04 cavit @ 1'
         call die_abort
      end select
c
c             step 4: verify completion of update with a result.
c                     update stress table and history
c
      if( new_state .eq. -1 .or. abs(T_new)+one .gt. mark .or.
     &  stiff_normal .le. zero .or. stiff_normal+one .gt. mark ) then
          write(iout,9300) felem+i-1, gpn
          write(iout,*) '..Fatal error mm04 cavit @ 2'
          call die_abort
      end if 
c
c             step 5: compute tangential traction and current
c                     tangential stiffness. just viscous flow model
c                     but shear viscosity degrades with a/b 
c                     (if a/b > ab_sd_ratio)
c
      f_sd = one
      if( degrade_shear ) then
        ab_ratio = a_new / b_new
        if( ab_ratio .gt. ab_sd_ratio ) 
     &     f_sd = one-(ab_ratio-ab_sd_ratio)/(one-ab_sd_ratio)
      end if
c
      trac_n1(i,1)  = f_sd * props%eta_b * delrlds(i,1) / dtime
      trac_n1(i,2)  = f_sd * props%eta_b * delrlds(i,2) / dtime
      stiff_shear   = f_sd * props%eta_b / dtime 
c       
      trac_n1(i,3)  = T_new
c      
      history1(i,1) = N_new / props%N_I
      history1(i,2) = a_new / props%a_0
      history1(i,3) = b_new / props%b_0
      history1(i,4) = V_new / props%V_0
      history1(i,5) = T_new  ! for states output
      if( compute_solid_local ) history1(i,5) = Tn_solid 
C      if( abs(T_new) .gt. 1000. .and. step .gt. 1) then
C         write(*,9500) abs_elem, T_new, top_sigma_m , bott_sigma_m,
C     &        top_sigma_e, bott_sigma_e 
C         write(*,9510) top_surf_stresses_n(i,1:6),
C     &       bott_surf_stresses_n(i,1:6)
C         if( step .eq. 3 ) stop    
C9500  format(".. big Tn, abs_elem, Tn: ", i6,f10.2,
C     &  /,10x," top_sigma_m , bott_sigma_m: ",2f10.2,
C     &  /,10x," top_sigma_e, bott_sigma_e : ",2f10.2 )
C 9510 format(10x," top stresses: ",6f10.2,
C     & /, 10x," bottom stresses: ",6f10.2 )     
C      end if
       
c      write(iout,*) "   ... Tnew, Tn_solid: ", T_new, Tn_solid           
c      if(  abs(Tn_solid)  .gt. 500.d00 ) stop  
           
      history1(i,6) = stiff_shear
      history1(i,7) = stiff_normal
      iword(1) = new_state; iword(2) = 0
      history1(i,13) = dword
      history1(i,11) = props%a_0 / props%b_0 ! for output 
      tshear_1 = trac_n1(i,1)
      tshear_2 = trac_n1(i,2)
      history1(i,14) = sqrt( tshear_1**2 + tshear_2**2 )
      if( debug_span_loop ) 
     &     write(iout,*) "..... @ 1 new_state: ", new_state
c
c             step 6: a / L_nr (needleman & rice parameter)
c                     sigma_e can be = 0 or -> 0
c
      history1(i,12) = 1.0d04  ! for sigma_e -> 0
      if( sigma_e .gt. LNR_toler ) history1(i,12) =  a_new * 
     &     ( d_strain / props%D / sigma_e )**third
c      
      end do    ! i over span
c      
      return
c
9000  format(10x,'.... mm04 cavit update. elem, gpn: ',i6,i2)
9100  format('...... entered mm04_cavit_sig_update ....') 
9110  format(15x,'i, hist1(i,6), hist1(i,7): ', i4,2f20.9)   
9200  format(15x,"N_n:                               ",e14.5,
     &     /,15x, "a_n, b_n, V_n:                    ",3e14.5,
     &     /,15x, "T_n, sigma_m, sigma_e:            ",3e14.5,
     &     /,15x, "d_strain, c_strain, reladis(i,3): ",3e14.5,
     &     /,15x, "delrlds(i,3):                     ",e14.5 )   
9300  format('>>>> FATAL ERRROR: element, gpn: ', 2i6 )
9310  format(15x,'reladis(,3), delrdis(,3): ',2e14.6, 
     &    /,15x,'opening, closing, neutral, penetrated: ',4l2)
9320  format(15x,"N_n:                               ",e14.5,
     &     /,15x, "a_n, b_n, V_n:                    ",3e14.5,
     &     /,15x, "T_n, sigma_m, sigma_e:            ",3e14.5,
     &     /,15x, "d_strain, c_strain:               ",3e14.5 )
c
c
      contains   !   ....  note this  ....
c     ========
c
c    ****************************************************************
c    *                                                              *
c    *          subroutines: mm04_compute_solid_local               *
c    *                                                              *
c    *               last modified:  8/13/2015 rhd                  *
c    *                                                              *
c    *                                                              *
c    ****************************************************************
c
      subroutine mm04_compute_solid_local( ie )
      implicit none
      integer :: ie
c
c             locals
c
#dbl      double precision ::
#sgl      real :
     & work(6), rotate(3,3)
c
      work(1:6) = half * ( sig_top(1:6) + sig_bott(1:6) )
c      
      rotate(1:3,1)  = intfprps(ie,51:53)    ! global to local tangent
      rotate(1:3,2)  = intfprps(ie,54:56)    ! normal. dir 3 is local
      rotate(1:3,3)  = intfprps(ie,57:59)    ! normal
c
      call mm04_tn_solid( rotate, work, Tn_solid, iout, .false. )
c
      return
      end subroutine mm04_compute_solid_local
c
c    ****************************************************************
c    *                                                              *
c    *          subroutines: mm04_cavit_update_linear               *
c    *                                                              *
c    *             last modified:  7/21/2015 rhd                    *
c    *                                                              *
c    *                                                              *
c    ****************************************************************
c
      subroutine mm04_cavit_update_linear
      implicit none
c 
      call mm04_cavit_linear_stiff  ! stiff at t = 0, rate = 0 
      stiff_normal = stiff_linear 
      if( reladis(i,3) .lt. zero ) 
     &     stiff_normal = stiff_linear * props%compression_mult
      T_new = stiff_normal * reladis(i,3)
      N_new = props%N_I
      a_new = props%a_0
      b_new = props%b_0
      V_new = props%V_0  
c
      return
      end   subroutine mm04_cavit_update_linear         
c
c    ****************************************************************
c    *                                                              *
c    *          subroutines: mm04_cavit_std_update                  *
c    *                                                              *
c    *             last modified:  11/20/2014 rhd                   *
c    *                              3/10/2015 kbc (VVNT, qmod)      *
c    *                              7/21/2-15 rhd (error check for  *
c    *                                a < a_0 )                     *
c    *                                                              *
c    ****************************************************************
c
      subroutine mm04_cavit_std_update
      implicit none
#dbl      double precision, parameter ::
#sgl      real, parameter ::
     & local_tol_a_0 = 0.90d0
c
c             step 1: compute f, then q(f), then c1.
c                     c_0, c_1, c_2. these routines will be inlined
c                     use contains subs for simplicity in running 
c                     algorithm
c
      if( VVNT ) call mm04_cavit_dV_L_or_H
c
      call mm04_cavit_c1
      call mm04_cavit_c0
      call mm04_cavit_c2
c      
      delta_c_dot = delrlds(i,3) / dtime 
c      
      if( debug_span_loop ) then
           write(iout,9210) delta_c_dot, ratio_1
           if( sigma_e .gt. toler ) 
     &        write(iout,9215) Lnr, ratio_2, term1, props%v2dot_term2,
     &                         term3, pm,
     &                         fraction, v2_dot
           write(iout,9220) S, nucleation_active, c_0, c_1, c_2
      end if            
c
c             step 2: update normal traction T. simple linear solve
c                     (c_2 = 0) or with Newton iteration (c_2 > 0)            
c                     store updated normal traction and normal
c                     direction stiffness.
c
      iterative_solve = c_2 .gt. zero
      if( iterative_solve ) then
         call mm04_cavit_newton( iout, c_2, c_1, c_0, trac_n(i,3),
     &                           delta_c_dot, delta_T, dtime, 
     &                           props%Sigma_0, converged,
     &                           debug_newton, iter_count )
         if( .not. converged ) then
            write(iout,9800) abs_elem, gpn, c_0, c_1, c_2,
     &                       delta_c_dot, dtime, props%Sigma_0   
            call die_abort
         end if
         T_new = trac_n(i,3) + delta_T
         if( debug_span_loop ) write(iout,9230) delta_T, T_new, 
     &                                          iter_count
      end if
      if( .not. iterative_solve ) then
         T_new = ( delta_c_dot - c_0 ) / c_1      
         if( debug_span_loop ) write(iout,9232) T_new
      end if   
c
c             step 3: update remaining internal state variables and
c                     update the history vector (see top of mm04)
c                     V2_dot does not depend on T_new  
c 
c                     if we find that a_n < a_0 (within small tol),
c                     the stress updating/solution process is in
c                     an inconsistent state.
c 
      V1_dot = four * pi * props%D * T_new / q_factor
      V_new  = V_n + dtime * (v1_dot + v2_dot)
      a_new  = a_n + dtime*(v1_dot + v2_dot) / 
     &      ( four * pi * a_n**2 * props%hpsi )
c     
      switch_to_linear = .false.
      if( a_new .lt. props%a_0 ) then
        switch_to_linear = .true.
        if( a_new .gt. local_tol_a_0 * props%a_0 ) return
        write(iout,9290) 
        write(iout,9300) felem+i-1, gpn
        write(iout,9310) reladis(i,3), delrlds(i,3),
     &                   opening, closing, neutral, penetrated    
        write(iout,9315) c_0, c_1, c_2, delta_c_dot, T_new, a_n  
        write(iout,9320) v1_dot, v2_dot, V_new, a_new
        write(iout,9330)
        return
      end if     
c     
c                     a_new cannot be > (tol) * b_new. set limit
c                     on a_new. the stress/state variable updating
c                     with finite load/time increments would otherwise
c                     allow a_new to exceed b_new.
c
      N_new = N_n
      b_new = b_n
c
      N_dot = zero
      nucleation_active = nucleation_active .and. T_new .gt. zero
      if( nucleation_active .and. include_nucleation ) then
        N_dot = props%F_N * d_strain * (T_new/props%Sigma_0)**2
        N_new = N_n  +  N_dot * dtime
        b_new = b_n - half * b_n * N_dot / N_new 
      end if 
c
      if( a_new / b_new .gt. max_ab_ratio ) then
        a_new = max_ab_ratio * b_new 
      end if  
c      
      stiff_normal = one / (dtime*(c_1 + two * c_2 * T_new ))
      new_state    = 1
c      
      if( debug_span_loop ) write(iout,9240) v1_dot, V_new, a_new,     
     &                        b_new, N_new, N_dot, nucleation_active,
     &                        stiff_normal, new_state
c
      return
      
9210  format(15x, "delta_c_dot, ratio_1:             ",2e14.5 )
9215  format(15x, "Lnr, ratio_2, term1:              ",3e14.5,
     &    /,15x,  "term2, term3, pm, fraction, v2_dot: ",
     &    2e14.5,f4.1,2e14.5 )
9220  format(15x, "S, nuc active:                    ",e14.5, 5x, l1,
     &    /,15x,  "c_0, c_1, c_2:                    " 3e14.5 )
9230  format(15x, "delta_T, T_(n+1), Newton iters:   ",2e14.5,i5 )
9232  format(15x, "linear solve, T_(n+1):            ",e14.5)
9240  format(15x, "v1_dot, V_new:                    ",2e14.5,
     &    /,15x,  "a_new, b_new:                     ",2e14.5,
     &    /,15x,  "N_new, N_dot, nuc active:         ",2e14.5,5x,l1,
     &    /,15x,  "stiff_normal, new_state:          ",e14.5, i5 )
9290  format(">>>>> Warning:  mm04_cavit_std_update @ 3" )
9300  format(15x,'element, gpn: ', i7, i2)
9310  format(15x,'reladis(,3) @ n+1, delrdis(,3) n->n+1: ',2e14.5, 
     &    /,15x,'opening, closing, neutral, penetrated: ',4l2)
9315  format(15x,'c_0, c_1, c_2: ',3e14.5, 
     &    /,15x,'delta_c_dot, T_new, a_n: ',3e14.5)
9320  format(15x, "v1_dot, v2_dot, V_new, a_new:     ",4e14.5)
9330  format(15x, "switch to linear stiffness path")

9800  format('>>>> FATAL ERROR. routine mm04_cavit_sig_update.',
     & /,    '                  iterations in mm04_nr_itr did not',
     & /,    '                  converge. ',
     & /,15x,"element, gpn:                     " i6,i5,
     & /,15x,"c_0, c_1, c_2:                    " 3e14.5,
     & /,15x,"delta_c_dot, dtime:               ",2e14.5,
     & /,15x,"Sigma_0:                          ",e14.5,
     & /,    "                  job terminated.", // )
      
      end subroutine mm04_cavit_std_update
c
c    ****************************************************************
c    *                                                              *
c    *          subroutines: mm04_cavit_dV_L_or_H                   *
c    *                                                              *
c    *             last modified:  2/25/2015 kbc                    *
c    *                                                              *
c    ****************************************************************
c
      subroutine mm04_cavit_dV_L_or_H
      implicit none
c
c             determine if V_dot_L or V_dot_H eqns govern volume 
c             change per VVNT equations (van der Giessen, van der 
c             Burg, Needleman, and Tvergaard). 
c             See ORNL/LTR-2014/492 pg 14-15.
c
#dbl      double precision ::
#sgl      real :
     &      v1_dot_L, v1_dot_H,
     &      v2_dot_L, v2_dot_H,
     &      v_h_term, beta_signed
          integer :: L_or_H ! = 0 if L, 1 if H
c
c             compute v2_dot_L and v2_dot_H
c
      v2_dot_L = zero
      v2_dot_H = zero
c
      if( sigma_e .gt. toler ) then ! ok to compute v2_dot
c      
        term3   = two * pi * d_strain * a_n**3 * props%hpsi
        pm      = one
        beta_signed = props%beta
        if( sigma_m .lt. zero ) then 
          pm = -one
          beta_signed = (props%n_pow - one) * 
     &                  (props%n_pow + 0.4031d00)/ props%n_pow**2
        endif
c        
        fraction = sigma_m / sigma_e 
        v_h_term = one/(one-(0.87d0*a_n/b_n)**(three/props%n_pow))
c        
        if( abs(fraction) .gt. one ) then
          v2_dot_L = pm * term3 * ( props%v2dot_term2*abs(fraction)
     &             + beta_signed )**props%n_pow
          v2_dot_H = pm * term3 * 
     &               (v_h_term*(props%v2dot_term2*abs(fraction)
     &                 + pm/props%n_pow))**props%n_pow
        else
          v2_dot_L = fraction * term3 * 
     &             (props%v2dot_term2+ beta_signed )**props%n_pow
          v2_dot_H = fraction * term3 * 
     &               (v_h_term*(props%v2dot_term2 + 
     &                pm/props%n_pow))**props%n_pow
        end if
c   
      end if     ! on sigma_e test
c
c             compute v1_dot_L and v1_dot_H
c
      call mm04_cavit_qfactor( 0 )
      v1_dot_L = four*pi*props%D*T_n/q_factor
c      
      call mm04_cavit_qfactor( 1 )
      v1_dot_H = four*pi*props%D*T_n/q_factor
c        
      if( debug_span_loop )
     &   write(iout,9010) v1_dot_L, v1_dot_H,
     &                   v2_dot_L, v2_dot_H,a_n/b_n,
     &                   sigma_m, sigma_e, a_n, b_n,
     &                   beta_signed, v_h_term
c
      if( abs(v1_dot_L + v2_dot_L) .ge. abs(v1_dot_H + v2_dot_H) ) then
        call mm04_cavit_qfactor( 0 )
        v2_dot = v2_dot_L
      else
        call mm04_cavit_qfactor( 1 )
        v2_dot = v2_dot_H      
      end if
c      
      return    
c
9010  format(15x, "v1_dot_L, v1_dot_H,v2_dot_L, v2_dot_H,a/b:: ",5e14.5,
     &    /,15x,  "sigma_m, sigma_e,                           ",2e14.5,
     &    /,15x,  "a, b                                        ",2e14.5,
     &    /,15x,  "beta_signed, v_h_term                       ",2E14.5)
c
      end subroutine mm04_cavit_dV_L_or_H
c
c    ****************************************************************
c    *                                                              *
c    *          subroutines: mm04_cavit_qfactor                     *
c    *                                                              *
c    *             last modified:  2/25/2015 kbc                    *
c    *                                                              *
c    ****************************************************************
c
      subroutine mm04_cavit_qfactor( L_or_H )
      implicit none
c
c             computes q (used in c1).
c
c             if q modification is on then the base equation 
c             is modified to transition to a minimum 
c             value (rather than allow q to 
c             get abitrarily small) over a range of f values.
c
c     input params
      integer :: L_or_H ! = 0 if L (i.e. using v_dot_L), 
c                         = 1 if H (v_dot_H)
c
c     local parameters
#dbl      double precision ::
#sgl      real :
     & f_bar_i, f_bar_j, q_i,dq_i, q_j, dq_j,
     & del, N1, M1, N2, M2,
     & q_min
c     
      q_min   = 0.001d0
      f_bar_i = 0.25d0
      f_bar_j = 0.49d0
c
      q_i = log(one/f_bar_i) - half*(three-f_bar_i)*(one-f_bar_i)
      dq_i = (two - f_bar_i) - one/f_bar_i
c
      q_j  = q_min
      dq_j = zero
c
      if( a_n .ge. b_n ) then
         write(iout,9000) a_n, b_n
         call die_abort
      end if
c
      ratio_1 = ( a_n / b_n )**2   
      if( L_or_H .eq. 0 .and. sigma_e .gt. toler ) then 
         Lnr = ( props%D * sigma_e / d_strain ) ** third
         ratio_2 = ( a_n / ( a_n + 1.5d0 * Lnr ) )**2
         f_factor = max( ratio_1, ratio_2 )
      else
         f_factor = ratio_1
      end if 
c       
      q_factor = log(one/f_factor) -
     &           half*(three-f_factor)*(one-f_factor)
c
      if( modify_q ) then 
        if( f_factor .gt. f_bar_i) then
          if( f_factor .gt. f_bar_j ) then 
            q_factor = q_min
          else
            del = two*(f_factor-f_bar_i)/(f_bar_j - f_bar_i) - one
            N1 = (one-del)*(one-del)*(two+del)*0.25d0
            M1 = 0.125d0*(f_bar_j-f_bar_i)*(one-del)*(one-del)*(one+del)
            N2 = (one+del)*(one+del)*(two-del)*0.25d0
            M2 = 0.125d0*(f_bar_j-f_bar_i)*(one+del)*(one+del)*(del-one)
            q_factor = N1*q_i + M1*dq_i + N2*q_j + M2*dq_j
          endif
        end if
      end if
c    
      return
c
9000  format('>>>> FATAL ERROR. routine mm04_cavit_qfactor.',
     & /,    '                  a_n .ge. b_n',
     & /,    '                  job terminated.', // )
9001  format('elm, gpn, f, q_factor:', i4, 1x, i2, 1x, e14.5, 
     &   1x, e14.5)
c
      end subroutine mm04_cavit_qfactor
c
c    ****************************************************************
c    *                                                              *
c    *          subroutines: mm04_cavit_c0, c1, c2                  *
c    *                                                              *
c    *             last modified:  3/20/2015 kbc                    *
c    *                                                              *
c    ****************************************************************
c
c
      subroutine mm04_cavit_c0
      implicit none
c
c                     handle corner case of v. small sigma_e,
c                     different equations for tension/compression
c                     mean stress.
c
c                     if VVNT eqns are used, then v2_dot is already
c                     calculated during the subroutine that determines
c                     if vdot_L or vdot_H is used. If the van der 
c                     Geissen equations for vdot are used, then 
c                     v2_dot must be computed here before computing c_0
c
      if( VVNT ) then
          c_0 = v2_dot * pi * b_n**2 
          return
      endif     
c
      c_0    = zero
      v2_dot = zero
c      
      if( sigma_e .gt. toler ) then ! ok to compute
        term3   = two * pi * d_strain * a_n**3 * props%hpsi
        pm      = one
        if( sigma_m .lt. zero ) pm = -one
        fraction = sigma_m / sigma_e 
        if( abs(fraction) .gt. one ) then
          v2_dot = pm * term3 * ( props%v2dot_term2*abs(fraction)
     &           + props%beta )**props%n_pow
        else
          v2_dot = term3 * fraction * 
     &          ( props%v2dot_term2 + props%beta )**props%n_pow
        end if
      end if
c      
      c_0 = v2_dot * pi * b_n**2  
c
      return
      end subroutine mm04_cavit_c0  
c 
      subroutine mm04_cavit_c1
      implicit none
c
      if( .not. VVNT ) call mm04_cavit_qfactor( 0 )  
      c_1 = four * props%D / ( b_n * b_n * q_factor )
c      
      return
      end subroutine mm04_cavit_c1  
c
      subroutine mm04_cavit_c2
      implicit none
c
c                     non-zero only if active nucleation
c
      c_2 = zero
      S   = zero
      if( T_n .gt. zero ) S = c_strain * ( T_n / props%Sigma_0 )**2
      nucleation_active = S .gt. props%Sthr .and. 
     &                    N_n .lt. props%N_max    
      if( nucleation_active .and. include_nucleation ) 
     &   c_2 = V_n * props%F_N * d_strain /
     &       ( pi * b_n**2 * N_n * props%Sigma_0**2 )
c     
      return 
      end subroutine mm04_cavit_c2  
c      
c    ****************************************************************
c    *                                                              *
c    *          subroutines: mm04_cavit_linear-stiff                *
c    *                                                              *
c    *             last modified:  3/20/2015 kbc                    *
c    *                                                              *
c    ****************************************************************
c
c
      subroutine mm04_cavit_linear_stiff
      implicit none
c          
      f_0 = (props%a_0 / props%b_0)**2
      q   = log(one/f_0) - half * (three-f_0) * (one-f_0)
      stiff_linear =  props%b_0**2 * q / ( four * props%D * dtime )
      if( debug_span_loop ) then
        write(iout,*) "..... a_0, b_0, f_0, q: ",props%a_0,
     &        props%b_0, f_0, q
         write(iout,*) "..... K-linear: " , stiff_linear
      end if
c      
      return
      end  subroutine mm04_cavit_linear_stiff      
      
      end subroutine mm04_cavit_sig_update
c
c    ****************************************************************
c    *                                                              *
c    *               function mm04_cavit_mises                      *
c    *                                                              *
c    *               written by rhd 3/15/2013                       *
c    *                                                              *
c    ****************************************************************
c
c
#dbl      double precision
#sgl      real
     & function mm04_cavit_mises( sig )
       implicit none
c
c             parameters
c
#dbl      double precision ::
#sgl      real ::
     &  sig(6)
c
c             locals
c
#dbl      double precision ::
#sgl      real ::
     & t1, t2, t3, t4, t5, roothalf, six
      data roothalf, six / 0.70710678119d0, 6.0d00 /
c
      t1 = sig(1) - sig(2)
      t2 = sig(3) - sig(2)
      t3 = sig(3) - sig(1)
      t4 = t1*t1 + t2*t2 + t3*t3
      t5 = sig(4)*sig(4) + sig(5)*sig(5) + sig(6)*sig(6)
      mm04_cavit_mises = roothalf * sqrt( t4 + six*t5 )
c
      return
      end
c
c    ****************************************************************
c    *                                                              *
c    *               subroutine mm04_cavit_newton                   *
c    *                                                              *
c    *               written by cj 7/14/2013                        *
c    *               updated by rhd 11/13/2014                      *
c    *                                                              *
c    *     local newton iteration support for cavit option          *
c    *                                                              *
c    ****************************************************************
c
c
      subroutine mm04_cavit_newton( iout, c_2, c_1, c_0, T_n,
     &                        delta_c_dot, delta_T, dtime, Sigma_0, 
     &                        converged, local_debug, iter_count )
      implicit none
c
c             parameters
c
#dbl      double precision ::
#sgl      real ::
     & c_2, c_1, c_0, T_n, delta_c_dot, delta_T, dtime, zero, Sigma_0
      integer :: iout, iter_count
      logical :: converged, local_debug
c
c             locals
c
#dbl      double precision ::
#sgl      real ::
     & toler, delta_T_old, delta_T_new, f, df, two, f_small
      logical :: test1, test2
      integer :: max_newton, iter
      data max_newton, toler, zero, two / 20, 1.d-5, 0.0d0, 2.0d0 /
c
c      if( local_debug ) write(iout,9000)
      converged   = .true.
      delta_T_old = zero
      f_small     = toler * abs( c_0 - delta_c_dot + 
     &              c_1*(T_n + delta_T_old) + 
     &              c_2*(T_n + delta_T_old)**2 )   
c      
      do iter = 1, max_newton
        f = c_0 - delta_c_dot + c_1*(T_n + delta_T_old) +
     &      c_2*(T_n + delta_T_old)**2    
        df = c_1 + two * c_2 * (T_n + delta_T_old) 
        delta_T_new = delta_T_old - f / df
        test1 = abs( f ) .le. f_small
        test2 = abs( delta_T_new - delta_T_old ) .lt. 
     &           toler*Sigma_0
c        if( local_debug ) then
c           write(iout,9010) iter, delta_T_old, f, df, delta_T_new, 
c     &                      test1, test2, f_small, toler*Sigma_0
c        end if
        if( test1 .and. test2 )then
            delta_T = delta_T_new
            iter_count = iter
            return
        else
            delta_T_old = delta_T_new
        end if    
      end do
c
      converged = .false.
c
      return
c    
9000  format(10x,"... inside  mm04_cavit_newton ...") 
9010  format(15x, "iter, delta_T_old:                ",i3,e14.5,
     &    /,15x,  "f, df, delta_T_new:               ",3e14.5,
     &    /,15x,  "test1, test2, f_small:            ",2l3,e14.5
     &    /,15x,  "toler*Sigma_0:                    ",e14.5 )
c
      end
      

c
c    ****************************************************************
c    *                                                              *
c    *               subroutine mm04_init                           *
c    *                                                              *
c    *                last modified:  8/13/2015 rhd                 *
c    *                                                              *
c    *     WARP3D internal routine to help setup mm04               *
c    *     this subroutine extracts the necessary material          *
c    *     properties for cohesive zone models. called directly by  *
c    *     drivers in WARP3D to set up intfprps for passing to      *
c    *     mm04                                                     *
c    *                                                              *
c    ****************************************************************
c
c
c      type of cohesive zone models
c             =  1 - linear elastic
c             =  2 - bilinear
c             =  3 - ramp
c             =  4 - exponential_1
c             =  5 - exponential_2
c             =  6 - PPR
c             =  7 - cavit
c
c                    -- for the original cohesive options
c
c   intfprps( ,2)  - linear interface stiffness in the longitudinal direction
c   intfprps( ,3)  - linear interface stiffness in the transverse direction
c   intfprps( ,4)  - linear interface stiffness in the normal direction
c   intfprps( ,5)  - peak normal stress of the interface
c   intfprps( ,6)  - peak shear stress of the interface
c   intfprps( ,7)  - shape parameter for bilinear and ramp
c   intfprps( ,8)  - second ( additional) shape parameter for ramp
c   intfprps( ,9)  - separation distance at peak stress in sliding
c   intfprps( ,10) - separation distance at peak stress in opening
c   intfprps( ,11) - equivalent critical separation distance
c   intfprps( ,12) - for type 4:
c                       (beta) a ratio to determine the equivalent
c                              separation under mixed mode loading
c                    for type 6:
c                       = 0.0 for real mixed-mode response
c                       = 2.0 for Mode I (normal mode) only
c                             response
c                       = -2.0 for shear mode only response
c
c   intfprps( ,22) - compression (normal) stiffness multiplier
c
c
c   intfprps( ,13) - = 0.0 element is active, =1.0 element is killed
c
c                    -- below additional for PPR option
c
c   intfprps( ,23) - G_normal
c   intfprps( ,24) - G_shear
c   intfprps( ,25) - ratio_normal
c   intfprps( ,26) - ratio_shear
c   intfprps( ,27) - shape_normal
c   intfprps( ,28) - shape_shear
c   intfprps( ,29) - debug cohesive material computations (>0.0)
c
c
c                    -- below additional for cavit option
c
c   intfprps( ,29) - debug cohesive material computations (>0.0)
c   intfprps( ,30) - lambda_1  -- not used as of 11/15/2014 --
c   intfprps( ,31) - lambda_2     ""
c   intfprps( ,32) - lambda_3     ""
c   intfprps( ,33) - lambda_4     "" 
c   intfprps( ,34) - lambda_5     ""
c   intfprps( ,35) - lambda_6 (creep exponent n)  ""
c   intfprps( ,36) - delta_star (for embrittlement)  ""
c   intfprps( ,37) - sigma_star (for embrittlement)  ""
c   intfprps( ,38) - lambda_7 ( a_0 / R_i )   ""
c   intfprps( ,39) - eta_b
c   intfprps( ,40) - R_I -- not allowed on input as of 3/1/2015
c   intfprps( ,41) - Sigma_0
c   intfprps( ,42) - t_c
c   intfprps( ,43) - D
c   intfprps( ,44) - a_0 -- must be actual values 3/1/2015
c   intfprps( ,45) - b_0 -- must be actual values 3/1/2015
c   intfprps( ,46) - psi_angle
c   intfprps( ,47) - N_I  -- not allowed on input as of 2/21/2015
c                            N_I computed from b_0
c   intfprps( ,48) - F_n
c   intfprps( ,49) - n
c   intfprps( ,50) - N_max
c
c
      subroutine mm04_init( iout, span, first_elem_in_blk,
     &                      props, lprops, iprops,
     &                      cohes_type, intfprps, matprp,
     &                      global_to_element_rot )
      implicit none
c
$add param_def
c
c          need mxvl, mxelpr from the param_def file.
c          mxvl - max allowable number of elements/block in system
c          mxelpr - max allowable properties for elements
c          mxvl typically 128 or 256
c          mxelpr typically 200
c
c
      integer iout, span, first_elem_in_blk, cohes_type
      real props(mxelpr,*)      ! note it is single precision !!!
      logical lprops(mxelpr,*)  ! same space as props
      integer iprops(mxelpr,*)  ! same space as props
      real matprp(*)            ! note it is single precision !!!
#dbl      double precision
#sgl      real
     &    intfprps(mxvl,*), global_to_element_rot(mxvl,3,3)
c
c                     locals
c
      integer i, j, iand, dummy, innum, cavit_loc, start_col, nvalues
      logical is_ppr, local_debug, is_cavit, bad
#dbl      double precision
#sgl      real
     &    zero, one, value
      data zero, one / 0.0d00, 1.0d00 /
c
c
c             original cohesive material options had their
c             props stuffed in with element properties. just
c             load for all options.
c
c             row 1 of intfprps no longer used. just caused confusion
c
      do i = 1, span
           intfprps(i,2)  =  props(7,i)
           intfprps(i,3)  =  props(8,i)
           intfprps(i,4)  =  props(9,i)
           intfprps(i,5)  =  props(13,i)
           intfprps(i,6)  =  props(14,i)
           intfprps(i,7)  =  props(15,i)
           intfprps(i,8)  =  props(28,i)
           intfprps(i,9)  =  props(20,i)
           intfprps(i,10) =  props(21,i)
           intfprps(i,11) =  props(22,i)
           intfprps(i,12) =  props(23,i)
           intfprps(i,13) =  iprops(32,i)
           intfprps(i,14) =  props(34,i)
           intfprps(i,15) =  props(35,i)
           intfprps(i,16) =  props(36,i)
           intfprps(i,17) =  props(37,i)
           intfprps(i,18) =  props(39,i)
           intfprps(i,19) =  props(40,i)
           intfprps(i,20) =  props(41,i)
           intfprps(i,21) =  props(33,i)
           intfprps(i,22) =  props(29,i)
           intfprps(i,29) = zero
           if( iand( iprops(24,i), 1 ) .ne. 0 ) intfprps(i,29) = one
           intfprps(i,20) = 10.0
      end do
c
c             verify that all elements in block are the same
c             type of coehsive material option.
c             should not get her unless same but check.
c
      bad = .false.
      do i = 1, span
       if( iprops(27,i) .ne. cohes_type ) then
          write(iout,9200) first_elem_in_blk + i - 1
          bad = .true.
       end if
      end do
      if( bad ) then
       write(iout,9210)
       call die_abort
      end if
c
      is_ppr   = cohes_type .eq. 6
      is_cavit = cohes_type .eq. 7
c
c             set up for ppr formulation.
c
      if( is_ppr ) then
         do i = 1, span
           intfprps(i,23) =  props(35,i)
           intfprps(i,24) =  props(36,i)
           intfprps(i,25) =  props(37,i)
           intfprps(i,26) =  props(39,i)
           intfprps(i,27) =  props(40,i)
           intfprps(i,28) =  props(41,i)
         end do
      end if
c
c              set up for cavit formulation.
c
      if( is_cavit ) then
         cavit_loc = 200 ! must match inmat.f
         start_col = 29  ! see table listing above
         nvalues   = 21  ! see table above
         do j = 1, nvalues
           value = matprp(cavit_loc+j)
           intfprps(1:span,start_col+j) = value
          end do
      end if
c
c              get global to interface element local rotation
c              matrix and put in intfprops for mm04 to use
c              if needed
c
      do i = 1, span
           intfprps(i,51) = global_to_element_rot(i,1,1)
           intfprps(i,52) = global_to_element_rot(i,2,1)
           intfprps(i,53) = global_to_element_rot(i,3,1)
           intfprps(i,54) = global_to_element_rot(i,1,2)
           intfprps(i,55) = global_to_element_rot(i,2,2)
           intfprps(i,56) = global_to_element_rot(i,3,2)
           intfprps(i,57) = global_to_element_rot(i,1,3)
           intfprps(i,58) = global_to_element_rot(i,2,3)
           intfprps(i,59) = global_to_element_rot(i,3,3)
      end do
c      
      local_debug = intfprps(1,29) .ne. zero
      if ( .not. local_debug ) return
c
      write(iout,*) "... inside mm04_init ..."
      do i = 1, span
        write(iout,*) "   > element: ", first_elem_in_blk+i-1
        do j = 2, 22
          write(iout,9100) j, intfprps(i,j)
        end do
        if( is_ppr ) then
           write(iout,*) "     .... option is ppr ...."
           do j = 23, 28
             write(iout,9100) j, intfprps(i,j)
           end do
        end if
        if( is_cavit ) then
          write(iout,*) "     .... option is cavit ...."
          do j = 1, nvalues
            write(iout,9100) j, intfprps(i,start_col+j)
          end do
        end if
        write(iout,*) " "
        write(iout,*) "     .... global -> element local ]R] ...."
        write(iout,9220) intfprps(i,51), intfprps(i,54), intfprps(i,57)
        write(iout,9220) intfprps(i,52), intfprps(i,55), intfprps(i,58)
        write(iout,9220) intfprps(i,53), intfprps(i,56), intfprps(i,59)
      end do
c
      return
 9100 format( "     intfprps(",i2.2,")",5x, f12.5)
 9200 format( "... internal error mm04_init. for element: ", i8,
     &    /,  "    different coehsive option that element 1 of blk" )
 9210 format( ">>>>> FATAL ERROR: job terminated",//)
 9220 format(15x,3f10.5)
      end
c
c    ****************************************************************
c    *                                                              *
c    *                   subroutine mm04_exp1_secant                *
c    *                                                              *
c    *                       written by : A Roy,                    *
c    *                                    S. Roychowdhury           *
c    *                    last modified : 7/3/2013 rhd              *
c    *                                                              *
c    *     computes the secant stiffness for a block of cohesive    *
c    *     elements for the exponential traction-separation model   *
c    *     Ref. IJNME 44,1267-1282 (1999)                           *
c    *                                                              *
c    ****************************************************************
c
      subroutine mm04_exp1_secant( span, felem, intfprps, dis,
     &                             history, history1, cohmat,
     &                             elem_killed, mxvl )
      implicit integer (a-z)
c
c                   parameter declarations
c
#dbl      double precision
#sgl      real
     & intfprps(mxvl,*), cohmat(mxvl,3), dis(mxvl,*),
     & history(span,*), history1(span,*)
       logical elem_killed(*)
c
c                   locally defined
c
#dbl      double precision
#sgl      real
     & zero, tol, dtol, d_eff_at_peak, beta, dt1, dt2, dn, b2,
     & effdis(mxvl), efftrac(mxvl), e, maxeffdis(mxvl),
     & prior_max_d_eff, d_eff_at_n, peak_intf_stress, comp_multiplier
       logical compression(mxvl), small_effdis(mxvl), loading(mxvl),
     &         local_debug
c
       data zero, e, one, dtol
     & / 0.0d0, 2.71828182845904523536d0, 1.0d0, 0.000001d0 /
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
c         see also file cnst4.f
c
      local_debug = .false.
c
      do i = 1, span
        if( elem_killed(i) ) cycle
        d_eff_at_peak   = intfprps(i,11)
        beta            = intfprps(i,12)
        d_eff_at_n      = history(i,1)
        prior_max_d_eff = history(i,2)
        tol             = dtol * d_eff_at_peak
        dt1 = dis(i,1); dt2 = dis(i,2); dn = dis(i,3)
        compression(i) = dn .lt. zero
        if( compression(i) ) dn = zero
        effdis(i) = sqrt( beta*beta*(dt1*dt1+dt2*dt2) + dn*dn )
        maxeffdis(i) = prior_max_d_eff
        if( effdis(i) .gt. prior_max_d_eff ) maxeffdis(i) = effdis(i)
        small_effdis(i)  = effdis(i) .le. tol
        loading(i) = effdis(i) .ge. d_eff_at_n   .and.
     &               effdis(i) .eq. maxeffdis(i)
      end do
c
      do k = 1, span
c
        cohmat(k,1:3) = zero
        if( elem_killed(k) ) cycle
c
        peak_intf_stress = intfprps(k,5)
        d_eff_at_peak    = intfprps(k,11)
        beta             = intfprps(k,12)
        comp_multiplier  = intfprps(k,22)
c
        if( small_effdis(k) ) then
              cohmat(k,1) = e*beta**2*peak_intf_stress
     &                      /d_eff_at_peak
              cohmat(k,2) = cohmat(k,1)
              cohmat(k,3) = e*peak_intf_stress/d_eff_at_peak
              if( compression(k) )
     &          cohmat(k,3) = cohmat(k,3) * comp_multiplier
              efftrac(k) = zero
              cycle
        end if
c
        efftrac(k) = e*peak_intf_stress*effdis(k)/d_eff_at_peak*
     &               e**(-one*effdis(k)/d_eff_at_peak)
        b2   = beta * beta
        if( loading(k) ) then
           if ( compression(k) ) then
               cohmat(k,1) = efftrac(k)/effdis(k)*b2
               cohmat(k,2) = cohmat(k,1)
               cohmat(k,3) = e*peak_intf_stress/d_eff_at_peak
               cohmat(k,3) = cohmat(k,3) * comp_multiplier
           else
               cohmat(k,1) = efftrac(k)/effdis(k)*b2
               cohmat(k,2) = cohmat(k,1)
               cohmat(k,3) = efftrac(k)/effdis(k)
           end if
        else  ! unloading. can still be compression
           cohmat(k,1) = e*b2*peak_intf_stress
     &                   /d_eff_at_peak*exp(-one*maxeffdis(k)
     &                   /d_eff_at_peak)
           cohmat(k,2) = cohmat(k,1)
           if ( compression(k) ) then
             cohmat(k,3) = comp_multiplier *
     &                     e*peak_intf_stress/d_eff_at_peak
           else
             cohmat(k,3) = e*peak_intf_stress/d_eff_at_peak
     &                   *exp(-one*maxeffdis(k)/d_eff_at_peak)
           end if
        end if
c
      end do
c
c                update history
c
       do k = 1, span
         if( elem_killed(k) ) cycle
         history1(k,1) = effdis(k)
         history1(k,2) = maxeffdis(k)
         history1(k,3) = efftrac(k)
       end do
c
       return
       end
c
c    ****************************************************************
c    *                                                              *
c    *               subroutine mm04_traction_ppr                   *
c    *                                                              *
c    *         written by : Kyoungsoo Park                          *
c    *                                                              *
c    *     this subroutine computes the cohesive traction           *
c    *     based on the unified potential-based model (PPR)         *
c    *     Ref. JMPS 57 (6), 891-908 (2009)                         *
c    *                                                              *
c    ****************************************************************
c
      subroutine mm04_traction_ppr( span, ppr_support, del, trac_n1,
     &             history, history1, elem_killed,
     &             local_debug, iout, mxvl )
      implicit none
      integer :: i, span, mxvl, iout
#dbl      double precision
#sgl      real
     &        ppr_support(mxvl,*), del(mxvl,*), trac_n1(mxvl,*),
     &        history(span,*), history1(span,*)
      logical :: elem_killed(*), local_debug
c
c               local definitions
c
#dbl      double precision
#sgl      real
     &        Gam_n, Gam_t, Tn_m, Tt_m, dn, dt, alph, beta, m, n,
     &        dGtn, dGnt, dnb, dtb,
     &        deln, delt, deln_max, delt_max, Tn, Tt, dTn, zero, tol,
     &        dtol, one, two, Dtt, Dnn, comp_multiplier
      logical :: compression
      data zero, one, two, tol
     & / 0.0d00, 1.0d00, 2.0d00, 0.00001d00 /
c
c              regular (nonlinear) stress updating
c              -----------------------------------
c
 100  continue
      do i = 1, span
        trac_n1(i,1) = zero
        trac_n1(i,2) = zero
        trac_n1(i,3) = zero
        if( elem_killed(i) ) cycle
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
        dtol  = tol * ( abs( deln ) + delt )
        deln_max = history(i,3)
        delt_max = history(i,4)
        compression = deln .lt. zero
c
c        Compute normal traction (Tn): 4 cases
c        Case 1: Compression region
c
        if (compression) then
          dTn = -Gam_n/dn**2*(m/alph)**(m-1)*(alph+m)*
     & (Gam_t*(1-delt/dt)**beta*(n/beta+delt/dt)**n + dGtn)
          comp_multiplier = ppr_support(i,20)
          Tn = dTn*deln*comp_multiplier
          deln = zero
c
c        Case 2: Complete failure
c
        elseif ((deln .ge. dn) .or. (delt .ge. dtb)
     &           .or. (deln_max .ge. dn)) then
          Tn = zero
c
c        Case 3: Softening region
c
        elseif (deln .ge. deln_max) then
          Tn = (Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*(m*(1-deln/dn)**alph*(m/alph+deln/dn)**(m-1)
     &       -alph*(1-deln/dn)**(alph-1)*(m/alph+deln/dn)**m)
          deln_max = deln
c
c         Case 4: Unloading/reloading condition
c
        else
          Tn = (Gam_t*(1-delt/dt)**beta*(delt/dt+n/beta)**n+dGtn) *
     & Gam_n/dn*(m*(1-deln_max/dn)**alph*(m/alph+deln_max/dn)**(m-1)
     &   -alph*(1-deln_max/dn)**(alph-1)*(m/alph+deln_max/dn)**m) *
     & deln/deln_max
        end if
c
c        Compute tangential traction (Tt): 3 cases
c        Case 1: Complete failure
c
        if ((delt .ge. dt) .or. (deln .ge. dnb)
     &       .or. (delt_max .ge. dt))  then
          Tt = zero
c
c        Case 2: Softening region
c
        elseif (delt .ge. delt_max) then
          Tt = (Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt/dt)**beta*(delt/dt+n/alph)**(n-1)
     &       -beta*(1-delt/dt)**(beta-1)*(delt/dt+n/beta)**n)
          delt_max = delt
c
c         Case 3: Unloading/reloading condition
c
        else
          Tt = (Gam_n*(1-deln/dn)**alph*(deln/dn+m/alph)**m+dGnt) *
     & Gam_t/dt*(n*(1-delt_max/dt)**beta*(delt_max/dt+n/alph)**(n-1)
     &   -beta*(1-delt_max/dt)**(beta-1)*(delt_max/dt+n/beta)**n) *
     & delt/delt_max
        end if
c
c        Store cohesive tractions
c
        if( delt .lt. dtol ) then
          trac_n1(i,1) = zero
          trac_n1(i,2) = zero
        else
          trac_n1(i,1) = Tt*del(i,1)/delt
          trac_n1(i,2) = Tt*del(i,2)/delt
        end if
        trac_n1(i,3) = Tn
c
c          Store new history variables for n+1
c
        history1(i,1) = dn
        history1(i,2) = dt
        history1(i,3) = deln_max
        history1(i,4) = delt_max
      end do
c
      if( local_debug ) then
        do i = 1, span
           write(iout,9200) i, del(i,1:3), trac_n1(i,1:3)
        end do
      end if
c
      return
 9100 format('... compute PPR tractions, step 1, iter 0')
 9200 format('... updated tractions, element in block: ', i4, /,
     &    10x, ' del(1-3):  ',3f15.10,/
     &    10x, ' tract(1-3):',3f15.4 )
      end
c    ****************************************************************
c    *                                                              *
c    *               subroutine mm04_init_ppr                       *
c    *                                                              *
c    *            written by : Kyoungsoo Park                       *
c    *                                                              *
c    *     this subroutine extracts user defined properties for     *
c    *     PPR and creates a table of suppporting factors           *
c    *                                                              *
c    ****************************************************************
c
c   ppr_support( ,1)  - = 0.0 (active element), = 1.0 (already
c                            killed element)
c   ppr_support( ,2)  - Normal fracture energy (Gn)
c   ppr_support( ,3)  - Tangential fracture energy (Gt)
c   ppr_support( ,4)  - Normal cohesive strength (sig_max)
c   ppr_support( ,5)  - Tangential cohesive strength (tau_max)
c   ppr_support( ,6)  - Normal shape parameter (shape_n)
c   ppr_support( ,7)  - Tangential shape parameter (shape_t)
c   ppr_support( ,8)  - Normal ratio of the critical separation to the
c                       final opening width (ratio_n)
c   ppr_support( ,9)  - Tangential ratio of the critical separation
c                       to the final opening width (ratio_t)
c   ppr_support( ,10) - Exponent m in the PPR model
c   ppr_support( ,11) - Exponent n in the PPR model
c   ppr_support( ,12) - Normal final crack opening width
c   ppr_support( ,13) - Tangential final crack opening width
c   ppr_support( ,14) - Energy constant associated with the
c                       opening mode
c   ppr_support( ,15) - Energy constant associated with the
c                       shearing mode
c   ppr_support( ,16) - <Gn - Gt>
c   ppr_support( ,17) - <Gt - Gn>
c                       KP- Note that <.> is the Macaulay bracket
c   ppr_support( ,18) - Conjugate normal final crack opening width
c   ppr_support( ,19) - Conjugate tangential final crack opening width
c   ppr_support( ,20) - compression multiplier
c
      subroutine mm04_init_ppr( span, intfprps, ppr_support,
     &                          elem_killed, iout, mxvl )
      implicit none
#dbl      double precision
#sgl      real
     &      intfprps(mxvl,*), ppr_support(mxvl,*)
      logical elem_killed(*)
      integer :: i, span, mxvl, iout
c
#dbl      double precision
#sgl      real
     &      Gn, Gt, Tn_m, Tt_m, alph, beta, ln, lt, m, n, dn, dt,
     &      dnb, dtb, Gam_n, Gam_t, dGnt, dGtn, zero, one,
     &      comp_multiplier
      data zero, one / 0.0d00, 1.0d00 /
c
c             The user defined material properties are:
c                  Gn             = intfprps(,23)
c                  Gt             = intfprps(,24)
c                  Tn_m = sig_max = intfprps(,4)
c                  Tt_m = tau_max = intfprps(,5)
c                  alph = shape_n = intfprps(,27)   [ shape_normal ]
c                  beta = shape_t = intfprps(,28)   [ shape_shear ]
c                  ln = ratio_n   = intfprps(,25)   [ ratio_normal ]
c                  lt = ratio_t   = intfprps(,26)   [ ratio_shear ]
c
      do i = 1, span
c
      if ( elem_killed(i) ) cycle
      Gn   = intfprps(i,23)
      Gt   = intfprps(i,24)
      Tn_m = intfprps(i,5)
      Tt_m = intfprps(i,6)
      alph = intfprps(i,27)
      beta = intfprps(i,28)
      ln   = intfprps(i,25)
      lt   = intfprps(i,26)
      comp_multiplier = intfprps(i,22)
c
      m = (alph-one)*alph*ln**2/(one-alph*ln**2)
      n = (beta-one)*beta*lt**2/(one-beta*lt**2)
      dn = alph*Gn/(m*Tn_m)*(one-ln)**(alph-one)
     &     * (alph/m*ln+one)**(m-one)*(alph+m)*ln
      dt = beta*Gt/(n*Tt_m)*(one-lt)**(beta-one)
     &     * (beta/n*lt+one)**(n-one)*(beta+n)*lt
c
      if (Gt .GT. Gn) then
       dGnt = zero
       dGtn = Gt - Gn
      elseif (Gt .LT. Gn) then
       dGnt = Gn - Gt
       dGtn = zero
      else
       dGnt = zero
       dGtn = zero
      endif
c
      if (Gn .EQ. Gt) then
       Gam_n = -Gn*(alph/m)**(m)
       Gam_t = (beta/n)**(n)
      else
       Gam_n = (-Gn)**(dGnt/(Gn-Gt))*(alph/m)**(m)
       Gam_t = (-Gt)**(dGtn/(Gt-Gn))*(beta/n)**(n)
      endif
      if (Gt .GT. Gn) then
       call mm04_get_dn_dt_bar (Gam_t, dt, n, beta, dGtn, dtb, iout )
       dnb = dn
      elseif (Gt .LT. Gn) then
       call mm04_get_dn_dt_bar (Gam_n, dn, m, alph, dGnt, dnb, iout)
       dtb = dt
      else
       dnb = dn
       dtb = dt
      endif
c
        ppr_support(i,2) = Gn
        ppr_support(i,3) = Gt
        ppr_support(i,4) = Tn_m
        ppr_support(i,5) = Tt_m
        ppr_support(i,6) = alph
        ppr_support(i,7) = beta
        ppr_support(i,8) = ln
        ppr_support(i,9) = lt
        ppr_support(i,10) = m
        ppr_support(i,11) = n
        ppr_support(i,12) = dn
        ppr_support(i,13) = dt
        ppr_support(i,14) = Gam_n
        ppr_support(i,15) = Gam_t
        ppr_support(i,16) = dGnt
        ppr_support(i,17) = dGtn
        ppr_support(i,18) = dnb
        ppr_support(i,19) = dtb
        ppr_support(i,20) = comp_multiplier
      end do
c
      return
      end
c
c    ****************************************************************
c    *                                                              *
c    *               subroutine mm04_get_dn_dt_bar                  *
c    *                                                              *
c    *           written by : Kyoungsoo Park                        *
c    *                                                              *
c    *     this subroutine estimates the conjugate final crack      *
c    *     opening width in the PPR cohesive zone models            *
c    *                                                              *
c    ****************************************************************
c
      subroutine mm04_get_dn_dt_bar( Gam_n, dn, m, alph, dGnt, dnb,
     &                               iout )
      implicit none
      integer :: itr, iout
#dbl      double precision
#sgl      real
     &        Gam_n, dn, m, alph, dGnt, dnb, fx, dfx
#dbl      double precision
#sgl      real
     &        half, tol, one
      data half, tol, one / 0.5d00, 1.0d-08, 1.0d00 /
c
      itr = 0
      dnb = half*dn
      fx = Gam_n*(1-dnb/dn)**alph*(m/alph+dnb/dn)**m + dGnt
c
c      Newton Raphson method is used to solve a nonlinear equation fx
c      The solution of the nonlinear equation fx is unique,
c      and exists between 0 and dn. The initial guess is set to be
c      0.5*dn. Please, see Ref. JMPS 57 (6), 891-908 (2009) for more detail.
c
      do while ((abs(fx) .gt. tol) .and. (itr .lt. 200))
        fx = Gam_n*(one-dnb/dn)**alph*(m/alph+dnb/dn)**m + dGnt
        dfx = (-alph*(one-dnb/dn)**(alph-one)*(m/alph+dnb/dn)**m +
     &        m*(one-dnb/dn)**alph*(m/alph+dnb/dn)**(m-one))*Gam_n/dn
        dnb = dnb - fx/dfx
        itr = itr + 1
      end do
      if( abs(fx) .gt. tol ) then
         write(iout,9000)
         call die_abort
      end if
c
      return
 9000 format('>>>>> FATAL ERROR: ppr option of cohesive model.',
     & /,'                   failed to converge in mm04_get_dn_dt_bar',
     & /,'                   analysis terminated.' )
      end
      
      
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm04_set_sizes                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 12/14/114 rhd               *
c     *                                                              *
c     *    called by warp3d for each material model to obtain        *
c     *    various sizes of data for the model                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm04_set_sizes( info_vector )
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
      info_vector(2) = 6
      info_vector(3) = 0
      info_vector(4) = 10
c
      return
      end
      
      
      
      
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm04_states_values                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 12/5/2014 (rhd)                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm04_states_values( itype, elem_states_output,
     &                                 nrow_states, num_states  )
c
c                       access some global data structures
c
      use elem_block_data, only: history_blocks, history_blk_list
      use main_data, only: elems_to_blocks, cohesive_ele_types
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
      integer :: relem, elnum, hist_size, blockno, cohesive_elem
      logical :: do_a_block, cohesive_type
      double precision :: zero
      data zero / 0.0d00 /
c      
c           build cohesive states values output.
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
      cohesive_type  = iprops(27,felem)
      cohesive_elem  = cohesive_ele_types(elem_type)
      
      if( local_debug ) write(out,9050) blockno, felem, elem_type,         
     &         mat_type, int_points, span, hist_size, cohesive_type,
     &         cohesive_elem
c
      if( cohesive_type .ne. 7 ) then
         write(out,9000) felem
         call die_abort
      end if  
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
           call mm04_states_values_a
           elem_states_output(1:nrow_states,relem) = 
     &                one_elem_states(1:nrow_states)
        end do
      else
        relem = elnum + 1 - felem
        one_elem_states(1:nrow_states) = zero
        call mm04_states_values_a
        elem_states_output(1:nrow_states,1) =
     &                one_elem_states(1:nrow_states)
      end if  
c        
      deallocate( history_dump, one_elem_states )
c        
      return
c      
 9050 format(10x,"block, felem, etype, mtype:  ",4i7,
     &  /,10x,   "int_pts, span, hist_size:    ",3i7,
     &  /,10x,   "cohesive_type, cohesive_elem:    ",i7,l7 )
 9000 format(/1x,
     &'>>>>> Error: only cavit option supported for states output ',
     & /,14x,'routine oustates_values_mm04, element: ',i7,
     & /,14x,'job terminated....'/)
c 
      contains
c     ========      
     
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm04_states_values_a                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/5/2015   (rhd)           *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm04_states_values_a
c   
      implicit none      
c
c                       locals
c
      integer :: ipt   
#dbl      double precision :: 
#sgl      real :
     & N_ratio, a_ratio, b_ratio, V_ratio, lambda_4, a_b_ratio,
     & a_bar, b_bar, a_L_ratio, T_n, T_nmax, T_shear, T_shr_max
c
      N_ratio   = zero
      a_ratio   = zero
      b_ratio   = zero
      V_ratio   = zero
      a_b_ratio = zero
      a_L_ratio = zero
      T_n       = zero
      T_nmax    = zero
      T_shear   = zero
      T_shr_max = zero
c       
      do ipt = 1, int_points
        a_bar     = history_dump(2,ipt,relem)
        b_bar     = history_dump(3,ipt,relem)  
        N_ratio   = N_ratio + history_dump(1,ipt,relem)
        a_ratio   = a_ratio + a_bar
        b_ratio   = b_ratio + b_bar
        V_ratio   = V_ratio + history_dump(4,ipt,relem)
        lambda_4  = history_dump(11,ipt,relem)
        a_b_ratio = a_b_ratio + ( a_bar / b_bar ) * lambda_4
        a_L_ratio = a_L_ratio + history_dump(12,ipt,relem)
        T_n       = T_n +  history_dump(5,ipt,relem)
        T_nmax    = T_nmax +  history_dump(9,ipt,relem)
        T_shear   = T_shear +  history_dump(14,ipt,relem)
        T_shr_max = T_shr_max +  history_dump(15,ipt,relem)
      end do
c
      one_elem_states(1)  = N_ratio / dble(int_points)            
      one_elem_states(2)  = a_ratio / dble(int_points)            
      one_elem_states(3)  = b_ratio / dble(int_points)            
      one_elem_states(4)  = V_ratio / dble(int_points)            
      one_elem_states(5)  = a_b_ratio / dble(int_points)            
      one_elem_states(6)  = a_L_ratio / dble(int_points)            
      one_elem_states(7)  = T_n / dble(int_points)            
      one_elem_states(8)  = T_nmax / dble(int_points)            
      one_elem_states(9)  = T_shear/ dble(int_points)            
      one_elem_states(10) = T_shr_max / dble(int_points)            
c
      return     
c
      end subroutine mm04_states_values_a
      end subroutine mm04_states_values
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm01_states_labels                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 3/4/2015  (rhd)                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm04_states_labels( size_state,
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
      num_states = 10  
      num_comment_lines  = 0   
      state_labels(1)  = "N / N_I"
      state_labels(2)  = "a / a_0"
      state_labels(3)  = "b / b_0"
      state_labels(4)  = "V / V_0"
      state_labels(5)  = "a / b"
      state_labels(6)  = "a / Lnr"
      state_labels(7)  = "Tnormal"
      state_labels(8)  = "Tn (max)"
      state_labels(9)  = "Tshear"
      state_labels(10) = "Ts (max)"
c      
      state_descriptors(1)  = "Normalized cavity density"
      state_descriptors(2)  = "Normalized cavity radius"
      state_descriptors(3)  = "Normalized cavity center-center spacing"
      state_descriptors(4)  = "Normalized cavity volume"
      state_descriptors(5)  = "Cavity radius over spacing"
      state_descriptors(6)  = "Needleman-Rice diffusivity-creep ratio"
      state_descriptors(7)  = "Normal traction"
      state_descriptors(8)  = "Max normal traction over history"
      state_descriptors(9)  = "Shear traction"
      state_descriptors(10) = "Max shear traction over history"
c
      if( do_print ) then
        do i = 1, 10 
          write(out,9010) i, state_labels(i), state_descriptors(i)
        end do
        do_print = .false.   
      end if
c          
      return
 9010 format(2x,i3,2x,a8,2x,a)      
      end
      
      
c *******************************************************************
c *                                                                 *
c *                 subroutine mm04_tn_solid                        *
c *                                                                 *
c *                       written by : rhd                          *
c *                                                                 *
c *               last modified : 8/13/2015  (rhd)                  *
c *                                                                 *
c *  Compute traction normal to interface from average of solid     *
c *  element stresses. for nonlocal implementation of cohesive      *
c *  material model. can be useful for checking                     *
c *                                                                 *
c *******************************************************************
c
c
      subroutine mm04_tn_solid( rotate, gsig, tn_solid, iout, debug )
      implicit none
c
      integer :: iout      
#dbl      double precision :: 
#sgl      real :
     & rotate(3,3), gsig(6), tn_solid
      logical ::   debug
c
c                       locals
c
#dbl      double precision :: 
#sgl      real :
     & nx, ny, nz, tx, ty, tz, tn   
c
c      global(1,1) = gsig(1)   !  comments for reference to ordering
c      global(2,1) = gsig(4) ! xy
c      global(3,1) = gsig(6) ! xz
c      global(1,2) = gsig(4)
c      global(2,2) = gsig(2) 
c      global(3,2) = gsig(5) ! yz
c      global(1,3) = gsig(6)
c      global(2,3) = gsig(5)
c      global(3,3) = gsig(3)
c          
c                       unit normal to interface
c
      nx = rotate(3,1); ny = rotate(3,2); nz = rotate(3,3)
c
c                       treat solid stresses as acting at point on the
c                       boundary having (unit) outward
c                       normal nx i + ny j + nz k
c                       apply cauchy relation to get traction vector
c                       on the interface, then dot with normal
c
      Tx = gsig(1)*nx + gsig(4)*ny + gsig(6)*nz
      Ty = gsig(4)*nx + gsig(2)*ny + gsig(5)*nz
      Tz = gsig(6)*nx + gsig(5)*ny + gsig(3)*nz
c      
      Tn_solid = Tx*nx + Ty*ny + Tz*nz
c
      return
      end
