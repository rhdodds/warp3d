c     ****************************************************************
c     *                                                              *
c     *                      subroutine gptns1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *               last modified : 06/18/12 rhd                   *
c     *                                                              *
c     *     this subroutine computes the contributon to the tangent  *
c     *     stiffnes matrices for a block of similar elements in     *
c     *     uniform global coordinates for a single gauss point.     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine gptns1( cp, icp, gpn, props, iprops, ek, local_work, 
     &      sz)
      use main_data, only: asymmetric_assembly
      use elem_block_data, only : global_cep_blocks => cep_blocks
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                     parameter declarations
c
      real props(mxelpr,*)
      integer cp(*), icp(mxutsz,*), iprops(mxelpr,*)
#dbl      double precision
#sgl      real
     & ek(sz,*)
c
c                     locals
c
#dbl      double precision
#sgl      real
     & eps_bbar, weight, rad(mxvl), dummy, factors(mxvl), one
      logical include_qbar, geonl, bbar, first, qbar_flag,
     &        temps_to_process, iscp
      data  one
#sgl     &  / 1.0 /
#dbl     &  / 1.0d00 /
c
c                       set local versions of the data structure
c                       scalars. set logical to include/not include the
c                       so-called qbar modification of the
c                       constitutive matrix. see addtional
c                       comments in routine ctran1
c
      etype            = local_work%elem_type
      span             = local_work%span
      felem            = local_work%felem
      utsz             = local_work%utsz
      eps_bbar         = local_work%eps_bbar
      nnode            = local_work%num_enodes
      totdof           = local_work%totdof
      geonl            = local_work%geo_non_flg
      utsz             = local_work%utsz
      bbar             = local_work%bbar_flg
      mat_type         = local_work%mat_type
      weight           = local_work%weights(gpn)
      first            = local_work%first
      iter             = local_work%iter
      qbar_flag        = local_work%qbar_flag
      int_order        = local_work%int_order
      temps_to_process = local_work%temps_node_to_process
      include_qbar     = qbar_flag
c
c                      for axisymmetric elements, compute the radius to
c                      the current gauss point for each element in the block.
c
      if( local_work%is_axisymm_elem ) then
        call get_radius( rad, nnode, span, local_work%shape(1,gpn),
     &                   local_work%ce, mxvl )
      end if
c
c                     compute [B] matrices at this gauss
c                     point for all elements in block. for geonl,
c                      a) [B] is the linear form evaluated using
c                         the updated coordinates at n+1
c                      b) the unrotated cauchy stresses at n+1
c                         are rotated to define the cauchy stresses
c                         at n+1 for initial stress stiffness
c                         and stress modification of [D].
c                         the [R,n+1] is retrieved from block storage
c                         by getrm1 and used to form the 6x6
c                         transformation: {T} = [qn1] * {uc}.
c                         [qn1] is also used to rotated [D] matrices
c                         from unrotated to spatial configuration.
c                         [D*] = trans[ qn1 ] * [D] * [qn1 ]
c                     routine bmod applies the b-bar modifications
c                     for both small and geonl formulations.
c                     note: [D] and [cep] mean the same here.
c                     the ctran1 routine performs the [D] rotation.
c
      if( local_work%is_cohes_elem ) then
           call blcmp_cohes( span, local_work%b_block,
     &                       local_work%cohes_rot_block,
     &                       local_work%shape(1,gpn), etype, nnode )
           go to 1000
      end if
c
      if( geonl ) then
          call getrm1( span, local_work%qn1,
     &                 local_work%rot_blk_n1(1,1,gpn), 2 )
          call qmply1( span, mxvl, nstr, local_work%qn1,
     &                 local_work%urcs_blk_n1(1,1,gpn),
     &                 local_work%cs_blk_n1 )
          call blcmp1_srt
     &         ( span, local_work%b_block,
     &                 local_work%gama_block(1,1,1,gpn),
     &                 local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &                 local_work%nzeta(1,gpn), local_work%shape(1,gpn),
     &                 local_work%ce, rad, etype, nnode )
      else
          call blcmp1_srt
     &         ( span, local_work%b_block,
     &                 local_work%gama_block(1,1,1,gpn),
     &                 local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &                 local_work%nzeta(1,gpn), local_work%shape(1,gpn),
     &                 local_work%ce, rad, etype, nnode )
      end if
c
      if( bbar ) then
          call bmod ( local_work%b_block, local_work%vol_block,
     &                span, mxvl, eps_bbar, mxedof )
      end if
c
c		     branch on material type:
c                      1 = simple mises model- linear hardening
c                          (isotropic or mixed kinematic), has geonl option
c                          rate effects on flow properties. temperature
c                          dependent flow properties, modulus, nu,
c                          alpha.
c                      2 = nonlinear elastic model (deformation plasticity
c                          with linear + power-law stress-strain
c                          relation. no geonl option.
c                      3 = general mises and gurson model -
c                          rate dependent/independent, linear (iso)
c                          hardening or power-law hardening. gurson
c                          model with w/o nucleation (strain controlled).
c                          temperature dependent flow props, modulus,
c                          nu, alpha
c                      4 = interface constitutive models
c                          supports several cohesive zone models.
c                          no geometric stiffness matrix. set geonl
c                          false to bypass those additional computations
c                      5 = adv. cyclic plasticity model
c                      6 = adv. gurson model
c                      7 = mises + hydrogen
c                      8 = Abaqus UMAT
c                     10 = CP model
c
c                    [Dt] computation for all models is vectorized.  [Dt]
c                    computed in unrotated configuration.
c
 1000 continue
      call iodevn( idummy, iout, dummy, 1 )
      select case ( mat_type )
      case ( 1 )
        call drive_01_cnst( gpn, iout, local_work )
      case ( 2 )
        call drive_02_cnst( gpn, iout, local_work )
      case ( 3 )
        call drive_03_cnst( gpn, iout, local_work )
      case ( 4 )
        call drive_04_cnst( gpn, iout, local_work )
        geonl = .false.
      case ( 5 )
        call drive_05_cnst( gpn, iout, local_work )
      case ( 6 )
        call drive_06_cnst( gpn, iout, local_work )
      case ( 7 )
        call drive_07_cnst( gpn, iout, local_work )
      case ( 8 )
        call drive_umat_cnst( gpn, iout, local_work )
      case (10 )
        call drive_10_cnst( gpn, iout, local_work )
      case default
          write(*,9500)
          call die_abort
      end select
c
c                       Check to see if we're actually a damaged interface
c                       model.  If we are, apply damage to the tangent.
      if (local_work%is_inter_dmg) then
        call drive_11_cnst(gpn, iout, local_work)
      end if

c
c                       store cep (i.e. [Dt]) for all elements in
c                       block at this integration point. scale out
c                       det[J] * weight included above. put in
c                       global blocked data structure. these
c                       will be used to compute stresses for iter=0
c                       if we have imposed displacement and/or
c                       temperature change over this step.
c
c                       No need to do this for UMAT since it only
c                       has [Dt] in global cep_blocks
c
      if( .not. local_work%is_umat ) then
         do i = 1, span
           factors(i) = one / ( local_work%det_jac_block(i,gpn)
     &                    * weight )
         end do
         now_blk = local_work%blk
         call tanstf_store_cep( span, mxvl, gpn,
     &                          local_work%cep_sym_size, local_work%cep,
     &                          global_cep_blocks(now_blk)%vector,
     &                          factors )
      end if
c
c                       convert [Dt] from unrotated cauchy to cauchy
c                       at current deformed configuration for geometric
c                       nonlinear analysis. no computations
c                       for cohesive or deformation plasticity. for UMAT with
c                       hyperelastic formulations which use [F] to get strains, the
c                       [Dt] stored in WARP3D is really for Cauchy stress - not
c                       unrotated Cauchy stress. The code below skips the
c                       rotation but may include the [Q] modification as
c                       requested in user input.
c
      if( .not. geonl ) go to 2000 ! used to have an iter=0 check after
      if( local_work%is_cohes_elem ) go to 2000
c
      if( local_work%is_deform_plas ) then
            write(*,9000)
            call die_abort
            stop
      end if
      call ctran1( span, local_work%cep, local_work%qn1,
     &                   local_work%cs_blk_n1,
     &                   include_qbar, local_work%det_jac_block(1,gpn),
     &                   weight, local_work%is_umat,
     &                   local_work%umat_stress_type,
     &                   local_work%is_crys_pls )

c
c                      include other required scalars with the material
c                      matrix for axisymmetric (triangle or
c                      quadrillateral) or planar triangle elements.
c                      A scalar to adjust the jacobian determinant
c                      is also needed for tetrahedral elements.
c
 2000 continue
      if( local_work%adjust_const_elem ) then
        call adjust_cnst( etype, nstr, mxvl, span, rad,
     &                    local_work%cep )
      end if
c
c                     compute each part of the element tangent
c                     stiffness matrix and add it in to the
c                     total. for geonl, we include initial stress
c                     stiffness. separate code (unrolled) is used for
c                     the 8-node brick for speed.
c
c
      if (asymmetric_assembly) then
        call bdbt_asym( span, icp, local_work%b_block, 
     &                  local_work%bt_block, local_work%bd_block,
     &                  local_work%cep, ek, mxvl, mxedof, utsz,
     &                  nstr, totdof, mxutsz)
      else

      if( totdof .ne. 24 )
     &   call bdbtgen( span, icp, local_work%b_block,
     &                 local_work%bt_block, local_work%bd_block,
     &                 local_work%cep, ek, mxvl, mxedof, utsz, nstr,
     &                 totdof, mxutsz )
c
      if( totdof .eq. 24 )
     &  call bdbt( span, icp, local_work%b_block, local_work%bt_block,
     &             local_work%bd_block, local_work%cep, ek,
     &             mxvl, mxedof, utsz, nstr, mxutsz )
      end if

      if( geonl ) then
        if (.not. asymmetric_assembly) then
          call kg1( span, cp, icp, local_work%gama_block(1,1,1,gpn),
     &            local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &            local_work%nzeta(1,gpn), nnode,
     &            local_work%cs_blk_n1,
     &            local_work%det_jac_block(1,gpn), weight, ek,
     &            local_work%vol_block, bbar, utsz, totdof )
         else
          call kg1( span, cp, icp, local_work%gama_block(1,1,1,gpn),
     &            local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &            local_work%nzeta(1,gpn), nnode,
     &            local_work%cs_blk_n1,
     &            local_work%det_jac_block(1,gpn), weight, ek,
     &            local_work%vol_block, bbar, totdof*totdof, totdof )
         end if
      end if
c
      return
c
 9000 format('FATAL ERROR: the nonlinear elastic material model',
     &   /,  '             is not valid for use with elements that',
     &   /,  '             have geometric nonlinearity.',
     &   /,  '             job terminated.' )
 9500 format(1x,'>> Fatal Error: gptns1. invalid material type..',
     &    /, 1x,'                job terminated' )

      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine drive_01_cnst                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/28/02                   *
c     *                                                              *
c     *     drive [D] consistent update for bilinear mises model     *
c     *     (has temperature and loading rate dependence, cyclic)    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_01_cnst( gpn, iout, local_work )
c
      use main_data, only : matprp, lmtprp
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                     parameter declarations
c
c                     local variables
c
#dbl      double precision
#sgl      real
     & weight
c
      etype            = local_work%elem_type
      nnode            = local_work%num_enodes
      span             = local_work%span
      felem            = local_work%felem
      weight           = local_work%weights(gpn)
c
c
c              set material states and get [cep]. transform
c              unrotated (material) -> spatial using qn1.
c              get temperature dependent modulus and nu if needed.
c
      call set_e_nu_for_block(
     &    span, local_work%nu_v, local_work%e_v, local_work%segmental,
     &    felem, etype, local_work%int_order, gpn, nnode,
     &    local_work%temps_node_to_process, local_work%temps_node_blk )
c
      if ( local_work%fgm_enode_props ) then
             call set_fgm_solid_props_for_block(
     &                span, felem, etype, gpn, nnode,
     &                local_work%e_v, local_work%shape(1,gpn),
     &                local_work%enode_mat_props, 1,
     &                local_work%fgm_flags(1,1) )
             call set_fgm_solid_props_for_block(
     &                span, felem, etype, gpn, nnode,
     &                local_work%nu_v, local_work%shape(1,gpn),
     &                local_work%enode_mat_props, 2,
     &                local_work%fgm_flags(1,2) )
      end if
c
      call cnst1( span, local_work%cep, local_work%rtse(1,1,gpn),
     &            local_work%nu_v,
     &            local_work%e_v, local_work%elem_hist1(1,2,gpn),
     &            local_work%elem_hist1(1,5,gpn), local_work%beta_v,
     &            local_work%elem_hist1(1,1,gpn),
     &            local_work%det_jac_block(1,gpn), weight,
     &            local_work%elem_hist1(1,4,gpn), felem )
c
      return
      end


c     ****************************************************************
c     *                                                              *
c     *                      subroutine drive_02_cnst                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/28/02                   *
c     *                                                              *
c     *     drive [D] consistent update for deformation plasticity   *
c     *     model                                                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_02_cnst( gpn, iout, local_work )
c
      use main_data, only : matprp, lmtprp
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                     parameter declarations
c
c                     local variables
c
#dbl      double precision
#sgl      real
     & weight
c
c                      use linear+power law nonlinear elastic
c                      material model. not valid for geometric
c                      nonlinear formulation. [cep] computation
c                      is fully vectorized.
c
      etype            = local_work%elem_type
      nnode            = local_work%num_enodes
      span             = local_work%span
      felem            = local_work%felem
      weight           = local_work%weights(gpn)
c
      if ( local_work%fgm_enode_props ) then
             call set_fgm_solid_props_for_block(
     &                span, felem, etype, gpn, nnode,
     &                local_work%e_v, local_work%shape(1,gpn),
     &                local_work%enode_mat_props, 1,
     &                local_work%fgm_flags(1,1) )
             call set_fgm_solid_props_for_block(
     &                span, felem, etype, gpn, nnode,
     &                local_work%nu_v, local_work%shape(1,gpn),
     &                local_work%enode_mat_props, 2,
     &                local_work%fgm_flags(1,2) )
             call set_fgm_solid_props_for_block(
     &                span, felem, etype, gpn, nnode,
     &                local_work%sigyld_v, local_work%shape(1,gpn),
     &                local_work%enode_mat_props, 7,
     &                local_work%fgm_flags(1,7) )
             call set_fgm_solid_props_for_block(
     &                span, felem, etype, gpn, nnode,
     &                local_work%n_power_v, local_work%shape(1,gpn),
     &                local_work%enode_mat_props, 8,
     &                local_work%fgm_flags(1,8) )
      end if
c
      call cnst2( felem, gpn, local_work%e_v, local_work%nu_v,
     &            local_work%sigyld_v,
     &            local_work%n_power_v,
     &            local_work%ddtse(1,1,gpn),
     &            local_work%elem_hist1(1,1,gpn), local_work%cep,
     &            span, local_work%det_jac_block(1,gpn), weight )
c
      return
      end


c     ****************************************************************
c     *                                                              *
c     *                   subroutine drive_03_cnst                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/28/02                   *
c     *                                                              *
c     *     drive [D] consistent update for mises, gt plasticity     *
c     *     model                                                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_03_cnst( gpn, iout, local_work )
c
      use main_data, only : matprp, lmtprp
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                     parameter declarations
c
c                     local variables
c
#dbl      double precision
#sgl      real
     & weight
      logical first
c
      etype            = local_work%elem_type
      nnode            = local_work%num_enodes
      span             = local_work%span
      felem            = local_work%felem
      weight           = local_work%weights(gpn)
      first            = local_work%first
      iter             = local_work%iter
c
c               general mises & Gurson rate dependent model
c
c         general mises/gurson model material model. the [cep] routine
c         needs all history values at n+1 but only 2 values fron
c         history at start of step (n). computation of [cep] is now
c         vectorized.
c
      call set_e_nu_for_block(
     &   span, local_work%nu_v, local_work%e_v, local_work%segmental,
     &   felem, etype, local_work%int_order, gpn, nnode,
     &   local_work%temps_node_to_process, local_work%temps_node_blk )
c
      if( local_work%fgm_enode_props ) then
             call set_fgm_solid_props_for_block(
     &                span, felem, etype, gpn, nnode,
     &                local_work%e_v, local_work%shape(1,gpn),
     &                local_work%enode_mat_props, 1,
     &                local_work%fgm_flags(1,1) )
             call set_fgm_solid_props_for_block(
     &                span, felem, etype, gpn, nnode,
     &                local_work%nu_v, local_work%shape(1,gpn),
     &                local_work%enode_mat_props, 2,
     &                local_work%fgm_flags(1,2) )
      end if
c
      call cnst3( felem, gpn, first, iter, local_work%e_v,
     &            local_work%nu_v,
     &            local_work%q1_v, local_work%q2_v,
     &            local_work%q3_v, local_work%nuc_v,
     &            local_work%nuc_s_n_v, local_work%nuc_e_n_v,
     &            local_work%nuc_f_n_v, local_work%rtse(1,1,gpn),
     &            local_work%elem_hist(1,1,gpn),
     &            local_work%elem_hist1(1,1,gpn), local_work%cep,
     &            span, local_work%det_jac_block(1,gpn), weight )
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                   subroutine drive_04_cnst                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/27/13 rhd               *
c     *                                                              *
c     *     drive [D] consistent update for cohesive model           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_04_cnst( gpn, iout, local_work )
c
      use main_data, only : matprp, lmtprp
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                     parameter declarations
c
c                     local variables
c
#dbl      double precision
#sgl      real
     & weight, time_n, dtime, ddummy(mxvl)
      integer idummy(mxvl)  ! and ddummy could be just (1)
      logical  nonlocal, temperatures, temperatures_ref
c
      etype            = local_work%elem_type
      nnode            = local_work%num_enodes
      span             = local_work%span
      felem            = local_work%felem
      weight           = local_work%weights(gpn)
      step             = local_work%step
      iter             = local_work%iter
      blk              = local_work%blk
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      knumthreads       = local_work%num_threads
      kthread           = omp_get_thread_num() + 1
c
      time_n = local_work%time_n
      dtime  = local_work%dt
c
      nonlocal = local_work%is_cohes_nonlocal
      imxvl = mxvl
      igpn = gpn
      kout = iout
c
      if( nonlocal ) then
        call cnst4(
     1    step, iter, felem, igpn, kout, span, imxvl, time_n, dtime,
     2    nonlocal, knumthreads, kthread, weight,
     3    local_work%cohes_type,local_work%intf_prp_block,
     4    local_work%ddtse(1,1,gpn),
     5    local_work%elem_hist(1,1,gpn),
     6    local_work%elem_hist1(1,1,gpn), local_work%cep,
     7    local_work%det_jac_block(1,gpn),
     8    local_work%cohes_temp_ref(1),
     9    local_work%cohes_dtemp(1),
     a    local_work%cohes_temp_n(1),  ! nonlocal after here
     a    local_work%top_surf_solid_elements(1),
     b    local_work%bott_surf_solid_elements(1),
     c    local_work%top_surf_solid_stresses_n1(1,1),
     d    local_work%bott_surf_solid_stresses_n1(1,1),
     e    local_work%top_surf_solid_eps_n1(1,1),
     f    local_work%bott_surf_solid_eps_n1(1,1),
     g    local_work%nonlocal_stvals_top_n1(1,1),
     h    local_work%nonlocal_stvals_bott_n1(1,1),
     i    local_work%top_solid_matl(1),
     j    local_work%bott_solid_matl(1) )
      else
        call cnst4(
     1    step, iter, felem, igpn, kout, span, imxvl, time_n, dtime,
     2    nonlocal, knumthreads, kthread, weight,
     3    local_work%cohes_type,local_work%intf_prp_block,
     4    local_work%ddtse(1,1,gpn),
     5    local_work%elem_hist(1,1,gpn),
     6    local_work%elem_hist1(1,1,gpn), local_work%cep,
     7    local_work%det_jac_block(1,gpn),
     8    local_work%cohes_temp_ref(1),
     9    local_work%cohes_dtemp(1),
     a    local_work%cohes_temp_n(1),
     a    idummy(1),
     b    idummy(1),
     c    ddummy(1),
     d    ddummy(1),
     e    ddummy(1),
     f    ddummy(1),
     g    ddummy(1),
     h    ddummy(1),
     i    idummy(1),
     j    idummy(1) )
      end if
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine drive_05_cnst                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 08/5/2011                  *
c     *                                                              *
c     *    drive [D] consistent update for cyclic plasticity model   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_05_cnst( gpn, iout, local_work )
c
      use main_data, only : matprp, lmtprp
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                     parameter declarations
c
c                     local variables
c
#dbl      double precision
#sgl      real
     & weight, gp_tau_vec(mxvl)
      logical first, nonlin_hard, generalized_pl,
     &        local_debug
      data local_debug / .false. /

      parameter(zero=0)
c
      span             = local_work%span
      felem            = local_work%felem
      weight           = local_work%weights(gpn)
      first            = local_work%first
      iter             = local_work%iter
      etype            = local_work%elem_type
      nnode            = local_work%num_enodes
c
      matnum           = local_work%matnum
      nonlin_hard      = matprp(58,matnum) .gt. zero
      generalized_pl   = matprp(58,matnum) .lt. zero
c
c
c                  storage layout for mm05_props:
c
c             local_work  matl storage  FA option      GP option
c                 1        55             q_u             gp_h_u
c                 2        56             b_u             gp_tau
c                 3        57             h_u             gp_beta_u
c                 4        58             1.0              -1.0
c                 5        59            gamma_u         gp_delta_u
c                 6        60            sig_tol         sig_tol
C                      61-64 <available>

      call set_mm05_props_for_block(
     &   span, local_work%nu_v, local_work%e_v, local_work%mm05_props,
     &   nonlin_hard, generalized_pl, local_work%segmental,
     &   felem, etype, local_work%int_order, gpn, nnode,
     &   local_work%temps_node_to_process, local_work%temps_node_blk,
     &   iout )
c
c            build remaining local data vectors generalized_plasticity.
c
      do i = 1, span
        gp_tau_vec(i) = matprp(56,matnum)
      end do

c
      if( local_debug ) then
       write(iout,*) '>.... drive_05_cnst:'
       write(iout,9000) span, felem, iter, first, matnum, gpn
       write(iout,9010) nonlin_hard, generalized_pl
       write(iout,9020)
       write(iout,9030) (felem+i-1, local_work%e_v(i),
     &                   local_work%nu_v(i),
     &                   local_work%mm05_props(i,1:5), i=1,span)
      end if
c
      call cnst5( span, felem, gpn, first, iter, iout, mxvl, nstr,
     &            weight, local_work%e_v, local_work%nu_v,
     &            local_work%mm05_props,
     &            local_work%rtse(1,1,gpn),
     &            local_work%elem_hist(1,1,gpn),
     &            local_work%elem_hist1(1,1,gpn),
     &            local_work%urcs_blk_n1(1,1,gpn),
     &            local_work%cep,
     &            local_work%det_jac_block(1,gpn),
     &            local_work%mm05_props(1,1),
     &            local_work%mm05_props(1,3),
     &            locaL_work%mm05_props(1,5), gp_tau_vec)
c
 9000 format(3x,'span, felem, iter, first, matnum, gpn: ',3i4,l4,2i4)
 9010 format(3x,'nonlin_hard, generalized_pl: ',2l4)
 9020 format(3x,"Property values:",/,9x,'elem',5x,'E         ','nu',
     &   15x,'mm05_props(1-5)')
 9030 format(5x,i8,f10.1,2x,f8.3,5f10.3)
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine drive_06_cnst                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/2/03                    *
c     *                                                              *
c     *     drive [D] consistent update for adv. gurson model        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_06_cnst( gpn, iout, local_work )
c
      use main_data, only : matprp, lmtprp
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                     parameter declarations
c
c                     local variables
c
#dbl      double precision
#sgl      real
     & weight
      logical first
c
      span             = local_work%span
      felem            = local_work%felem
      weight           = local_work%weights(gpn)
      first            = local_work%first
      iter             = local_work%iter
c
      call cnst6( span, felem, gpn, first, iter, iout, mxvl, nstr,
     &            weight, local_work%e_v, local_work%nu_v,
     &            local_work%sigyld_v, local_work%n_power_v,
     &            local_work%f0_v,
     &            local_work%q1_v,
     &            local_work%q2_v, local_work%q3_v,
     &            local_work%nuc_v, local_work%nuc_s_n_v,
     &            local_work%nuc_e_n_v, local_work%nuc_f_n_v,
     &            local_work%rtse(1,1,gpn),
     &            local_work%mm06_props,
     &            local_work%elem_hist(1,1,gpn),
     &            local_work%elem_hist1(1,1,gpn),
     &            local_work%urcs_blk_n1(1,1,gpn),
     &            local_work%cep,
     &            local_work%det_jac_block(1,gpn),
     &            local_work%killed_status_vec )
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine drive_07_cnst                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/2/03                    *
c     *                                                              *
c     *     drive [D] consistent update for adv. mises model +       *
c     *     hydrogen effects                                         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_07_cnst( gpn, iout, local_work )
c
      use main_data, only : matprp, lmtprp
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                     parameter declarations
c
c                     local variables
c
#dbl      double precision
#sgl      real
     & weight
      logical first
c
      span             = local_work%span
      felem            = local_work%felem
      weight           = local_work%weights(gpn)
      first            = local_work%first
      iter             = local_work%iter
c
      call cnst7( span, felem, gpn, first, iter, iout, mxvl, nstr,
     &            weight, local_work%e_v, local_work%nu_v,
     &            local_work%sigyld_v, local_work%n_power_v,
     &            local_work%mm07_props,
     &            local_work%rtse(1,1,gpn),
     &            local_work%elem_hist(1,1,gpn),
     &            local_work%elem_hist1(1,1,gpn),
     &            local_work%urcs_blk_n1(1,1,gpn),
     &            local_work%cep,
     &            local_work%det_jac_block(1,gpn) )
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                   subroutine drive_umat_cnst                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 05/14/12                   *
c     *                                                              *
c     *          drive [D] consistent update for warp3d umat         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_umat_cnst( gpn, iout, local_work )
c
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
      implicit integer (a-z)
$add param_def
$add include_tan_ek  ! has local_work definition
c
c                     local variables
c
#dbl      double precision
#sgl      real
     & weight, symm_part_cep(21), factor
      logical local_debug, debug_now
c
c           1. pull a few values from work space for block
c
      span              = local_work%span
      felem             = local_work%felem
      weight            = local_work%weights(gpn)
      now_blk           = local_work%blk
      local_debug       =  .false.
c
c           2. the tangent [D] matrices are stored in the
c              global_cep_blocks. only lower symmetric terms
c              stored in a vector. stored there in rstgp1.
c
c              for this integration point (gpn), process all
c              elements in block.
c
      do ielem = 1, span
c
        noel = felem + ielem - 1
        debug_now = local_debug
        if( debug_now ) write(iout,9100) ielem, noel, gpn
c
c                2 (a) pull 21 terms of lower-triangle for this element
c                      global cep block is 21 x span x num integration
c                      points
c
        start_loc = ( 21 * span * (gpn-1) ) + 21 * (ielem-1)
        do k = 1, 21
          symm_part_cep(k) = gbl_cep_blocks(now_blk)%vector(start_loc+k)
        end do
c
c                2 (b) expand to 6 x 6 symmetric [D] and scale by
c                      integration weight factor
c
        factor = weight * local_work%det_jac_block(ielem,gpn)
        k = 1
        do i = 1, 6
         do j = 1, i
           local_work%cep(ielem,i,j) = symm_part_cep(k) * factor
           local_work%cep(ielem,j,i) = symm_part_cep(k) * factor
           k = k + 1
         end do
        end do
c
        if( debug_now ) write(iout,9110) symm_part_cep(1:4)
      end do
c
      if( debug_now ) write(iout,9900)
      return
c
 9000 format('>> Enter UMAT cnst driver...')
 9001 format(5x,'span, felem, gpn, iter: ',i4,i10,i3,i4 )
 9002 format(5x,'integration weight: ',f10.5)
 9006 format(5x,'num history terms: ',i4 )
 9100 format(5x,'... processing i, elem, gpn: ',i4,i10, i3)
 9110 format(5x,'symd 1-4: ',4e14.6)
 9900 format('>> Leave UMAT cnst driver...')
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine drive_10_cnst                *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 03/28/12                   *
c     *                                                              *
c     *     drive [D] consistent update for warp3d umat              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_10_cnst( gpn, iout, local_work )
c
      use main_data, only : matprp, lmtprp
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                     parameter declarations
c
c                     local variables
c
#dbl      double precision
#sgl      real
     & weight
      logical first
c
      span             = local_work%span
      felem            = local_work%felem
      weight           = local_work%weights(gpn)
      first            = local_work%first
      iter             = local_work%iter
c           This really is not a good thing, but it works b/c
c           the element history arrays are allocatable
      hist_size_for_blk = size(local_work%elem_hist,2)
c
      call cnst10( gpn, span,
     &           hist_size_for_blk,
     &           local_work%elem_hist(1:span,1:hist_size_for_blk,gpn),
     &           local_work%elem_hist1(1:span,1:hist_size_for_blk,gpn),
     &           local_work, iout)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine drive_11_cnst                *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 03/11/14                   *
c     *                                                              *
c     *     drive damage update to the tangent matrix                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_11_cnst( gpn, iout, local_work )
c
      use main_data, only : matprp, lmtprp
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                     parameter declarations
c
c                     local variables
c
#dbl      double precision
#sgl      real
     & weight
      logical first
c
      span             = local_work%span
      felem            = local_work%felem
      weight           = local_work%weights(gpn)
      first            = local_work%first
      iter             = local_work%iter
c
      dmg_loc = 1+local_work%macro_sz
c
      call cnst11( gpn, span, local_work%cep(1:span,1:6,1:6), 
     &            local_work%elem_hist(1:span,dmg_loc,gpn),
     &            local_work%sv, local_work%lv, local_work%tv, iout)
c
      return
      end



c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine bdbt                         *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 07/07/91                   *
c     *                                 : 01/31/96                   *
c     *                                                              *
c     *     this subroutine performs the multiplication of the       *
c     *     transpose of the strain-displacement matrix by the       *
c     *     constituitive matrix by the strain-displacement matrix   *
c     *     that is a building block of a stiffness matrix.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine bdbt( span, icp, b, bt, bd, d, ek,
     &                 mxvl, mxedof, utsz, nstr, mxutsz )
      implicit integer (a-z)
c
c                       parameter declarations
c
#dbl      double precision
#sgl      real
     &   b(mxvl,mxedof,*), ek(utsz,*), d(mxvl,nstr,*),
     &   bd(mxvl,mxedof,*), bt(mxvl,nstr,*)
      integer icp(mxutsz,*)

c                       set b transpose.
c
      do j = 1, 24
        do i = 1, span
         do k = 1, 6
          bt(i,k,j)= b(i,j,k)
         end do
       end do
      end do
c
c                       perform multiplication of b*d. do the
c                       full matrix. call a subroutine since
c                       D can be treated as 2-d inside.
c
c                       [D] for solids is 6x6. Also for cohesive,
c                       i: 4->6 and j: 4->6 should be zeroed by
c                       cnst.. routine.
c
      do j = 1, 24
       do i = 1, span
         do k = 1, 6
           bd(i,j,k) = d(i,1,k) * b(i,j,1)
     &               + d(i,2,k) * b(i,j,2)
     &               + d(i,3,k) * b(i,j,3)
     &               + d(i,4,k) * b(i,j,4)
     &               + d(i,5,k) * b(i,j,5)
     &               + d(i,6,k) * b(i,j,6)
         end do
       end do
      end do
c
c                       perform multiplication of bd*b . do
c                       only for upper triangular entries.
c                       do it in-line to reduce an enormous
c                       number of subroutine calls
c
      do j = 1, 300
        row = icp(j,1)
        col = icp(j,2)
        do i = 1, span
         ek(j,i) = ek(j,i)
     &         +   bt(i,1,col) * bd(i,row,1)
     &         +   bt(i,2,col) * bd(i,row,2)
     &         +   bt(i,3,col) * bd(i,row,3)
     &         +   bt(i,4,col) * bd(i,row,4)
     &         +   bt(i,5,col) * bd(i,row,5)
     &         +   bt(i,6,col) * bd(i,row,6)
        end do
      end do
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine bdbt_asym                    *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 12/9/13                    *
c     *                                                              *
c     *     B * D * B.T for the asymmetric case.  Handles all        *
c     *     elements, at least for now.                              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine bdbt_asym( span, icp, b, bt, bd, d, ek,
     &                    mxvl, mxedof, utsz, nstr, totdof, mxutsz )
      implicit integer (a-z)
c
c                       parameter declarations
c
#dbl      double precision
#sgl      real
     &   b(mxvl,mxedof,*), ek(totdof*totdof,*), d(mxvl,nstr,*),
     &   bd(mxvl,mxedof,*), bt(mxvl,nstr,*)
      integer icp(mxutsz,*)

      do i=1,span
        ek(:,i) = ek(:,i) + reshape(matmul(b(i,1:totdof,1:6),
     &      matmul(d(i,1:6,1:6), transpose(b(i,1:totdof,1:6)))),
     &      (/totdof*totdof/))
      end do
      
      return
      end subroutine


c     ****************************************************************
c     *                                                              *
c     *                      subroutine bdbtgen                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 07/10/97                   *
c     *                                                              *
c     *     this subroutine performs the multiplication of the       *
c     *     transpose of the strain-displacement matrix by the       *
c     *     constituitive matrix by the strain-displacement matrix   *
c     *     that is a building block of a stiffness matrix. this     *
c     *     routine handles any type of element.                     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine bdbtgen( span, icp, b, bt, bd, d, ek,
     &                    mxvl, mxedof, utsz, nstr, totdof, mxutsz )
      implicit integer (a-z)
c
c                       parameter declarations
c
#dbl      double precision
#sgl      real
     &   b(mxvl,mxedof,*), ek(utsz,*), d(mxvl,nstr,*),
     &   bd(mxvl,mxedof,*), bt(mxvl,nstr,*)
      integer icp(mxutsz,*)

c                       set b transpose.
c
      do j = 1, totdof
        do i = 1, span
         do k = 1, 6
          bt(i,k,j)= b(i,j,k)
         end do
       end do
      end do
c
c                       perform multiplication of b*d. do the
c                       full matrix. call a subroutine since
c                       D can be treated as 2-d inside.
c
c                       [D] for solids is 6x6. Also for cohesive,
c                       i: 4->6 and j: 4->6 should be zeroed by
c                       cnst.. routine.

c
      do j = 1, totdof
       do i = 1, span
         do k = 1, 6
           bd(i,j,k) = d(i,1,k) * b(i,j,1)
     &               + d(i,2,k) * b(i,j,2)
     &               + d(i,3,k) * b(i,j,3)
     &               + d(i,4,k) * b(i,j,4)
     &               + d(i,5,k) * b(i,j,5)
     &               + d(i,6,k) * b(i,j,6)
         end do
       end do
      end do
c
c                       perform multiplication of bd*b . do
c                       only for upper triangular entries.
c                       do it in-line to reduce an enormous
c                       number of subroutine calls
c
      do j = 1, utsz
        row = icp(j,1)
        col = icp(j,2)
        do i = 1, span
         ek(j,i) = ek(j,i)
     &         +   bt(i,1,col) * bd(i,row,1)
     &         +   bt(i,2,col) * bd(i,row,2)
     &         +   bt(i,3,col) * bd(i,row,3)
     &         +   bt(i,4,col) * bd(i,row,4)
     &         +   bt(i,5,col) * bd(i,row,5)
     &         +   bt(i,6,col) * bd(i,row,6)
        end do
      end do
c
      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine tanstf_store_cep                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 71//12 rhd                 *
c     *                                                              *
c     *        support for storing plain [D] w/o weights             *
c     *                                                              *
c     ****************************************************************
c
      subroutine tanstf_store_cep( span, mxvl, gpn,
     &   cep_sym_size, local_cep_block, global_cep_block, factors )
      implicit none
c
      integer span, mxvl, gpn, nrow_cep_for_matl,
     &        cep_sym_size, nrow, ielem, k, i, j
#dbl      double precision
#sgl      real
     &    local_cep_block(mxvl,6,6),
     &    global_cep_block(cep_sym_size,span,*), factors(*)
c
c         local_cep_block contains the 6x6 [D] for this
c         integration point (gpn) for all elements in the block.
c         it is symmetric and already
c         has the integration weight factor already included. cancel that
c         with factors to just leave the simple [D]. store the lower
c         triangle of [D] in global_cep_block.
c
c         Note: the block local [D] is always 6x6 allocated even for
c         cohesive elements.
c
c         Use two sets of loops to expose actual sizes to compiler
c
      if( cep_sym_size .eq. 6 ) go to 1000
c
      k = 1
      do i = 1, 6
        do j = 1, i
          do ielem = 1, span
            global_cep_block(k,ielem,gpn) = local_cep_block(ielem,i,j)
     &                                       * factors(ielem)
          end do
          k = k + 1
        end do
      end do
      return
c
 1000 continue       ! cohesive elements
      k = 1
      do i = 1, 3
        do j = 1, i
          do ielem = 1, span
            global_cep_block(k,ielem,gpn) = local_cep_block(ielem,i,j)
     &                                       * factors(ielem)
          end do
          k = k + 1
        end do
       end do
      return
c
      end


c     ****************************************************************
c     *                                                              *
c     *                      subroutine ctran1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 06/18/12 rhd               *
c     *                                                              *
c     *     transform [Dt] from a form relating the unrotated stress *
c     *     rate and the unrotated rate of deformation tensor to     *
c     *     one relating the cauchy stress rate and the rate of      *
c     *     (spatial) deformation strain for a block of elements.    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ctran1( span, cep, qn1, cs, qbar, dj, w, is_umat,
     &                   umat_stress_type, is_crys_pls  )
      implicit integer (a-z)
$add param_def
#dbl      double precision
#sgl      real
     &     cep(mxvl,nstr,*), qn1(mxvl,nstr,*), tc(mxvl,nstr,nstr),
     &     cs(mxvl,*), half, two, dj(*), w, wf, halfw
      double precision z(6,6), zt(6,6)
      logical qbar, is_umat, is_crys_pls
      data half, two / 0.5, 2.0 /
c
c
c             [cep] (mxvl x 6 x 6) relates increments
c             of unrotated cauchy stress to increments
c             of the unrotated deformation. transform [cep]
c             so it relates increments of cauchy stress to
c             increments of the deformation, both on the
c             spatial coordinates.
c
c             [cep*] = [qn1] * [cep] * trans([qn1])
c
c             [qn1] is a rotation matrix constructed from the
c             [R] obtained by polar decomposition of the deformation
c             gradient, [F] =[R][U].
c
c             for UMATs the UMAT computed [cep] may already refer to
c             the Cauchy stress. No rotation to be done.
c
c             For crystal plasticity model,
c
      do_transform = .true.
      if( is_crys_pls ) do_transform = .true.
      if( is_umat ) then
        if( umat_stress_type .eq. 1 ) do_transform = .false.
      end if
      if( .not. do_transform ) go to 1000
c
c
c             perform multiplication of [tc] = [qn1] * [cep]
c
      do j = 1, nstr
         do i = 1, span
c
            tc(i,j,1)= (qn1(i,j,1)*cep(i,1,1)+
     &                  qn1(i,j,2)*cep(i,2,1)+
     &                  qn1(i,j,3)*cep(i,3,1)+
     &                  qn1(i,j,4)*cep(i,4,1)+
     &                  qn1(i,j,5)*cep(i,5,1)+
     &                  qn1(i,j,6)*cep(i,6,1))
c
            tc(i,j,2)= (qn1(i,j,1)*cep(i,1,2)+
     &                  qn1(i,j,2)*cep(i,2,2)+
     &                  qn1(i,j,3)*cep(i,3,2)+
     &                  qn1(i,j,4)*cep(i,4,2)+
     &                  qn1(i,j,5)*cep(i,5,2)+
     &                  qn1(i,j,6)*cep(i,6,2))
c
            tc(i,j,3)= (qn1(i,j,1)*cep(i,1,3)+
     &                  qn1(i,j,2)*cep(i,2,3)+
     &                  qn1(i,j,3)*cep(i,3,3)+
     &                  qn1(i,j,4)*cep(i,4,3)+
     &                  qn1(i,j,5)*cep(i,5,3)+
     &                  qn1(i,j,6)*cep(i,6,3))
c
            tc(i,j,4)= (qn1(i,j,1)*cep(i,1,4)+
     &                  qn1(i,j,2)*cep(i,2,4)+
     &                  qn1(i,j,3)*cep(i,3,4)+
     &                  qn1(i,j,4)*cep(i,4,4)+
     &                  qn1(i,j,5)*cep(i,5,4)+
     &                  qn1(i,j,6)*cep(i,6,4))
c
            tc(i,j,5)= (qn1(i,j,1)*cep(i,1,5)+
     &                  qn1(i,j,2)*cep(i,2,5)+
     &                  qn1(i,j,3)*cep(i,3,5)+
     &                  qn1(i,j,4)*cep(i,4,5)+
     &                  qn1(i,j,5)*cep(i,5,5)+
     &                  qn1(i,j,6)*cep(i,6,5))
c
            tc(i,j,6)= (qn1(i,j,1)*cep(i,1,6)+
     &                  qn1(i,j,2)*cep(i,2,6)+
     &                  qn1(i,j,3)*cep(i,3,6)+
     &                  qn1(i,j,4)*cep(i,4,6)+
     &                  qn1(i,j,5)*cep(i,5,6)+
     &                  qn1(i,j,6)*cep(i,6,6))
c
         end do
      end do
c
c
c                       perform multiplication of
c                       [cep*] =  [tc] * transpose([qn1])
c
      do j = 1, nstr
         do i = 1, span
c
            cep(i,j,1)= tc(i,j,1)*qn1(i,1,1)+
     &                  tc(i,j,2)*qn1(i,1,2)+
     &                  tc(i,j,3)*qn1(i,1,3)+
     &                  tc(i,j,4)*qn1(i,1,4)+
     &                  tc(i,j,5)*qn1(i,1,5)+
     &                  tc(i,j,6)*qn1(i,1,6)
c
            cep(i,j,2)= tc(i,j,1)*qn1(i,2,1)+
     &                  tc(i,j,2)*qn1(i,2,2)+
     &                  tc(i,j,3)*qn1(i,2,3)+
     &                  tc(i,j,4)*qn1(i,2,4)+
     &                  tc(i,j,5)*qn1(i,2,5)+
     &                  tc(i,j,6)*qn1(i,2,6)
c
            cep(i,j,3)= tc(i,j,1)*qn1(i,3,1)+
     &                  tc(i,j,2)*qn1(i,3,2)+
     &                  tc(i,j,3)*qn1(i,3,3)+
     &                  tc(i,j,4)*qn1(i,3,4)+
     &                  tc(i,j,5)*qn1(i,3,5)+
     &                  tc(i,j,6)*qn1(i,3,6)
c
            cep(i,j,4)= tc(i,j,1)*qn1(i,4,1)+
     &                  tc(i,j,2)*qn1(i,4,2)+
     &                  tc(i,j,3)*qn1(i,4,3)+
     &                  tc(i,j,4)*qn1(i,4,4)+
     &                  tc(i,j,5)*qn1(i,4,5)+
     &                  tc(i,j,6)*qn1(i,4,6)
c
            cep(i,j,5)= tc(i,j,1)*qn1(i,5,1)+
     &                  tc(i,j,2)*qn1(i,5,2)+
     &                  tc(i,j,3)*qn1(i,5,3)+
     &                  tc(i,j,4)*qn1(i,5,4)+
     &                  tc(i,j,5)*qn1(i,5,5)+
     &                  tc(i,j,6)*qn1(i,5,6)
c
            cep(i,j,6)= tc(i,j,1)*qn1(i,6,1)+
     &                  tc(i,j,2)*qn1(i,6,2)+
     &                  tc(i,j,3)*qn1(i,6,3)+
     &                  tc(i,j,4)*qn1(i,6,4)+
     &                  tc(i,j,5)*qn1(i,6,5)+
     &                  tc(i,j,6)*qn1(i,6,6)
c
         end do
      end do
c
c            subtract the [Q-bar] matrix from the transformed
c            [cep]. this is the "initial stress" at the material
c            point level. this remains an option indicated by qbar.
c            note: we must multiply in the gauss weight factor and
c            gauss point det[J] for the subtracted terms. the [cep]
c            passed in had these factors included by the cnst...
c            routines. the [Q-bar] 6x6 comes from the tensor
c            expression -2 (de.De):s, where, s is the stress tensor,
c            de is the rate of deformation tensor and De is the virtual
c            rate of deformation tensor. this expression in matrix form
c            is: - trans([B]) * [Q-bar] * [B]. this modification of [cep]
c            is essential for convergence of nearly homogeneous
c            deformation problems.
c
 1000 continue
      if ( qbar ) then
        do i = 1, span
         wf    = dj(i) * w
         halfw = half * wf
         cep(i,1,1) = cep(i,1,1) - two * cs(i,1) * wf
         cep(i,2,2) = cep(i,2,2) - two * cs(i,2) * wf
         cep(i,3,3) = cep(i,3,3) - two * cs(i,3) * wf
         cep(i,4,1) = cep(i,4,1) - cs(i,4) * wf
         cep(i,6,1) = cep(i,6,1) - cs(i,6) * wf
         cep(i,4,2) = cep(i,4,2) - cs(i,4) * wf
         cep(i,5,2) = cep(i,5,2) - cs(i,5) * wf
         cep(i,5,3) = cep(i,5,3) - cs(i,5) * wf
         cep(i,6,3) = cep(i,6,3) - cs(i,6) * wf
         cep(i,4,4) = cep(i,4,4) - halfw * ( cs(i,1)+cs(i,2) )
         cep(i,5,5) = cep(i,5,5) - halfw * ( cs(i,2)+cs(i,3) )
         cep(i,6,6) = cep(i,6,6) - halfw * ( cs(i,1)+cs(i,3) )
         cep(i,5,4) = cep(i,5,4) - halfw * cs(i,6)
         cep(i,6,4) = cep(i,6,4) - halfw * cs(i,5)
         cep(i,6,5) = cep(i,6,5) - halfw * cs(i,4)
         cep(i,1,4) = cep(i,4,1)
         cep(i,1,6) = cep(i,6,1)
         cep(i,2,4) = cep(i,4,2)
         cep(i,2,5) = cep(i,5,2)
         cep(i,3,5) = cep(i,5,3)
         cep(i,3,6) = cep(i,6,3)
         cep(i,4,5) = cep(i,5,4)
         cep(i,4,6) = cep(i,6,4)
         cep(i,5,6) = cep(i,6,5)
c
c                      experiment with symmetrized version of the nonsymmetric term
c                      see Crisfield vol. 2, pg55.
c
c      z(1:6,1:6) = 0.0d0      
C      z(1,1:3) = cs(i,1)      
C      z(2,1:3) = cs(i,2)      
C      z(3,1:3) = cs(i,3)      
C      z(4,1:3) = cs(i,4)      
C      z(5,1:3) = cs(i,5)      
C      z(6,1:3) = cs(i,6)      
c      zt = transpose( z)
c      do k = 1, 6
c      do l = 1, 6
c        z(k,l) = 0.5d0*(z(k,l)+zt(l,k))
c        cep(i,k,l) = cep(i,k,l) + z(k,l)*wf
c      end do
c      end do
c         
        end do
      end if
c
      return
      end





