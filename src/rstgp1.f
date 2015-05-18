c     ****************************************************************
c     *                                                              *
c     *                      subroutine rstgp1                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/2/2012  rhd              *
c     *                                                              *
c     *     supervise the computation of strains, stresses and       *
c     *     accompaning stress data at an integration point          *
c     *     for a block of similar elements that use the same        *
c     *     material model code                                      *
c     *                                                              *
c     *              ** geometric nonlinear version **               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rstgp1( props, lprops, iprops, local_work )
      use segmental_curves, only : max_seg_points
c
      implicit integer (a-z)
$add param_def
c
c                      parameter declarations
c
      real    props(mxelpr,*)
      logical lprops(mxelpr,*)
      integer iprops(mxelpr,*)
$add include_sig_up
c
c                       locally defined variables
c
#dbl      double precision
#sgl      real
     &  internal_energy, beta_fact, eps_bbar, plastic_work
#dbl      double precision,
#sgl      real,
     & allocatable :: ddt(:,:), uddt(:,:), qnhalf(:,:,:),
     &                qn1(:,:,:)
c
      logical cut_step,
     &        adaptive, geonl, bbar, material_cut_step,
     &        local_debug, adaptive_flag
c
      data local_debug / .false. /
c
      internal_energy   = local_work%block_energy
      plastic_work      = local_work%block_plastic_work
      beta_fact         = local_work%beta_fact
      span              = local_work%span
      felem             = local_work%felem
      type              = local_work%elem_type
      order             = local_work%int_order
      gpn               = local_work%gpn
      ngp               = local_work%num_int_points
      nnode             = local_work%num_enodes
      ndof              = local_work%num_enode_dof
      geonl             = local_work%geo_non_flg
      step              = local_work%step
      iter              = local_work%iter
      bbar              = local_work%bbar_flg
      mat_type          = local_work%mat_type
      material_cut_step = local_work%material_cut_step
      adaptive_flag     = local_work%adaptive_flag
      eps_bbar          = local_work%eps_bbar
      adaptive          = adaptive_flag .and. step .gt. 1
      iout              = local_work%iout
      if( local_debug ) write(iout,*) '... in rstgp1'
c
      allocate( ddt(mxvl,nstr), uddt(mxvl,nstr),
     &          qnhalf(mxvl,nstr,nstr), qn1(mxvl,nstr,nstr) )
c
c        process cohesive elements separately
c
      if( local_work%is_cohes_elem ) then
        if ( local_debug ) write(*,*) '>> calling gtlsn2...'
        call gtlsn2( span, nnode,
     &               local_work%due, uddt,
     &               local_work%ddtse(1,1,gpn),
     &               local_work%cohes_rot_block,
     &               local_work%shape(1,gpn),
     &               local_work%elem_type,gpn, felem, iout )
        go to 7000
      end if
c
c        compute deformation gradients (F), perform polar decompositions,
c        F=RU, etc. to get strain increment over the step on the
c        unrotated configuration.
c
c        calculate the element displacements at n+1/2 and n+1
c
      call rstgp1_a( ndof, nnode, span, local_work%ue,
     &               local_work%due, local_work%uenh,
     &               local_work%uen1, mxvl )
c
c        find the deformation gradients F=RU and stress/strain
c        transformation matrices at n+1/2, n+1. [R,n+1] is computed
c        and stored in block structure rot_blk_n1. [qnhalf], [qn1]
c        and {dfn1} are returned. dfn1 is det[F,n+1] for energy
c        integration. if we get a bad deformation jacobian (det <= 0,
c        terminate strain computations and request
c        immediate step size reduction if possible.
c
c        we store F at n and n+1 in local_work for use by WARP3D UMAT
c        and crystal plasticity if they need them
c        (for large displacement analysis).
c
c        qn1 is computed at present only for UMAT and crystal plasticity
c
      error = 0
      call gtmat1( qnhalf, qn1, error, local_work )
      if ( error .eq. 1 ) then
         if ( adaptive ) then
            material_cut_step = .true.
            local_work%material_cut_step = material_cut_step
            go to 9999
         else
            write(iout,9820)
            call abort_job
         end if
      end if
c
c         compute the deformation tensor increment ddt (often called D).
c         ddt returned in vector (6x1) form. we use the linear [B] matrix
c         evaluated at n+1/2 configuration * the step displacement increment.
c         (displacement increment also equals the velocity * dt)
c
      call gtlsn1( span, nnode,
     &             local_work%due, ddt,
     &             local_work%gama_mid(1,1,1,gpn),
     &             local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &             local_work%nzeta(1,gpn),
     &             local_work%vol_block, bbar, eps_bbar,
     &             local_work%b )
c
c         compute the unrotated increment (uddt) of the deformation tensor
c         (also called delta-d) (the "d" rate * dt ). uddt = [qnhalf] * ddt.
c         the stress update is driven by uddt.
c         add delta-d to acumulated (integral) over all steps. The unrotated
c         increments and total all refer to the fixed (global) coordinate
c         axes.
c
      call qmply1( span, mxvl, nstr, qnhalf, ddt, uddt )
      call rstgp1_update_strains( span, mxvl, nstr, uddt,
     &                            local_work%ddtse(1,1,gpn) )
c
c------------------------------------------------------------------------
c
 7000 continue
      select case ( mat_type )
      case ( 1 )
c
c                vectorized mises plasticity model
c
       call drive_01_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
c
      case( 2 )
c
c                linear+power law deformation plasticity model.
c                not supported for finite strain solutions,
c
       write(iout,9000)
       call die_abort
       stop
c
      case( 3 )
c
c                general mises/gurson flow theory model.
c
       call drive_03_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
c
      case( 4 )
c
c                linear and non-linear cohesive model
c
       call drive_04_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
c
      case ( 5 )
c
c                cyclic plasticity model
c
       call drive_05_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
      case ( 6 )
c
c                adv. gurson-tvergaard model
c
       call drive_06_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
      case ( 7 )
c
c                mises + hydrogen effects
c
       call drive_07_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
c
      case ( 8 )
c
c                general UMAT
c
       call drive_umat_update( gpn, local_work, uddt, qn1, iout )
c
      case ( 9 )
c
c                ALCOA anisotropic plasticty
c
       call drive_09_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
c
      case ( 10 )
c
c               CP model
c
      call drive_10_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout)
c
      case default
        write(iout,*) '>>> invalid material model number'
        write(iout,*) '    in rstgp1'
        call die_abort
c
      end select
c
c     If we're actually an interface damage model, call the interface
c     damage calculations
c
      if (local_work%is_inter_dmg) then
            call drive_11_update(gpn, props, lprops, iprops,
     &            local_work, uddt, iout)
      end if

c
c --------------------------------------------------------------------
c
c          calculate the internal energy and plastic work
c          y integrating the densities over the deformed volume of the element.
c          urcs_blk_n1(..,7) is really the current (total) energy
c          density per unit deformed volume - look above...
c          increment of plastic work density stored in plastic_work_incr
c
      if ( iter .ne. 0 ) then
        call rstgp1_b( span, internal_energy, plastic_work,
     &                 local_work%urcs_blk_n1(1,7,gpn),
     &                 local_work%urcs_blk_n1(1,8,gpn),
     &                 local_work%det_j(1,gpn), local_work%dfn1, 1 )
        internal_energy = internal_energy * beta_fact *
     &                    local_work%weights(gpn)
        plastic_work    = plastic_work * beta_fact *
     &                    local_work%weights(gpn)
        local_work%block_energy       = internal_energy
        local_work%block_plastic_work = plastic_work
      end if
c
      if ( local_debug ) then
        write(iout,*) '>> rstgp1 .. gauss point: ', gpn
        write (iout,9500) internal_energy,
     &                 plastic_work
        write(iout,9110)
        do i = 1, span
         write(iout,9100) i, (local_work%urcs_blk_n(i,k,gpn),k=1,7),
     &                 (local_work%urcs_blk_n1(i,k,gpn),k=1,7),
     &                 (ddt(i,k),k=1,6),
     &                 (uddt(i,k),k=1,6)
        end do
      end if
c
 9999 continue
      deallocate( ddt, uddt, qnhalf, qn1 )
c
      return
c
 9000 format('>>> Fatal Error: the nonlinear elastic material',
     &     /,'                 is not compatible with large',
     &     /,'                 displacement elements.',
     &     /,'                 job terminated....' )
 9100 format(i5,7f15.6,/,5x,7f15.6,/,5x,6f15.6,/,5x,6f15.6)
 9110 format(1x,'Elem    /',20('-'),
     &        ' unrot. Cauchy @ n, unrot. Cauchy @ n+1,',
     &        ' ddt, uddt', 20('-'),'/')
 9500 format('  Internal energy inside of (rstgp1)    = ',e16.6,
     &     /,'  Plasstic work inside of (rstgp1)      = ',e16.6)
 9820 format(///,
     &       '>> FATAL ERROR: strain computation routines requested',
     &     /,'                immediate step size reduction due to',
     &     /,'                invalid determinant of a deformation',
     &     /,'                jacobian. the user has not allowed step',
     &     /,'                size reductions. analysis terminated.',
     &     /// )
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine rstgp2                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/20/2011 rhd             *
c     *                                                              *
c     *     supervise the computation of strains, stresses and       *
c     *     accompaning stress data at an integration point          *
c     *     for a block of similar elements that use the same        *
c     *     material model code                                      *
c     *                                                              *
c     *              ** geometric linear version **                  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rstgp2( props, lprops, iprops, local_work )
      use segmental_curves, only : max_seg_points, max_seg_curves
c
      implicit integer (a-z)
$add param_def
c
c                      parameter declarations
c
      real    props(mxelpr,*)
      logical lprops(mxelpr,*)
      integer iprops(mxelpr,*)
$add include_sig_up
c
c                       locally defined variables
c
#dbl      double precision
#sgl      real
     &  internal_energy, beta_fact, eps_bbar,
     &  uddt(mxvl,nstr), plastic_work, dummy_q(1)
c
      logical bbar, local_debug, signal_flag
c
      data local_debug / .false. /

      internal_energy   = local_work%block_energy
      plastic_work      = local_work%block_plastic_work
      beta_fact         = local_work%beta_fact
      span              = local_work%span
      felem             = local_work%felem
      type              = local_work%elem_type
      order             = local_work%int_order
      gpn               = local_work%gpn
      ngp               = local_work%num_int_points
      nnode             = local_work%num_enodes
      ndof              = local_work%num_enode_dof
      step              = local_work%step
      iter              = local_work%iter
      bbar              = local_work%bbar_flg
      mat_type          = local_work%mat_type
      signal_flag       = local_work%signal_flag
      eps_bbar          = local_work%eps_bbar
      number_points     = local_work%number_points
      curve_set         = local_work%curve_set_number
      iout              = local_work%iout
c
c          compute the strain increment in vector form. update the
c          accumulated strains. branch to call material model.
c              uddt        -> strain increment over step
c              ddtse       -> strain at n, updated to strain at
c                             n+1 here.
c              strain_n    -> strain at n for models that need it
c                             (e.g. Abaqus UMAT)
c              urcs_blk_n  -> stresses at n
c              urcs_blk_n1 -> updated stresses at n+1
c
      if ( local_work%is_cohes_elem ) then
        call gtlsn2( span, nnode,
     &               local_work%due, uddt,
     &               local_work%ddtse(1,1,gpn),
     &               local_work%cohes_rot_block,
     &               local_work%shape(1,gpn),
     &               local_work%elem_type, gpn, felem, iout )
      else
        call gtlsn1( span, nnode,
     &               local_work%due, uddt,
     &               local_work%gama(1,1,1,gpn), local_work%nxi(1,gpn),
     &               local_work%neta(1,gpn),
     &               local_work%nzeta(1,gpn), local_work%vol_block,
     &               bbar, eps_bbar, local_work%b )
        call rstgp1_update_strains( span, mxvl, nstr, uddt,
     &                              local_work%ddtse(1,1,gpn) )
      end if
c
c------------------------------------------------------------------------
c
      select case (  mat_type )
      case ( 1 )
c
c                vectorized mises plasticity model
c
       call drive_01_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
c
      case( 2 )
c
c                linear+power law deformation plasticity model,
c                which is mostly vectorized.
c
       call drive_02_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
c
      case( 3 )
c
c                general mises/gurson flow theory model.
c
       call drive_03_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
c
      case ( 4 )
c
c                linear and non-linear cohesive model
c
       call drive_04_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
      case ( 5 )
c
c                cyclic plasticity model
c
       call drive_05_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
      case ( 6 )
c
c                adv. gurson-tvergaard model
c
       call drive_06_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
      case ( 7 )
c
c                mises + hydrogen effects
c
       call drive_07_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
      case ( 8 )
c
c                general UMAT -- Abaqus compatible
c
       call drive_umat_update( gpn, local_work, uddt, dummy_q, iout )
c
      case ( 9 )
c
c                ALCOA anisotropic plasticity
c
       call drive_09_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
c
      case ( 10 )
c
c                 Crystal plasticity.  Note this call (linear geometry)
c                 makes very little sense, and I reserve the right to
c                 throw an error in the material model
       call drive_10_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
c      
      case default
        write(iout,*) '>>> invalid material model number'
        call die_abort
        stop
      end select
c
c
c     If we're actually an interface damage model, call the interface
c     damage calculations
c
      if (local_work%is_inter_dmg) then
            call drive_11_update(gpn, props, lprops, iprops,
     &            local_work, uddt, iout)
      end if
c
c --------------------------------------------------------------------
c
c            calculate the current (total) internal energy
c            by integrating the energy density over the volume
c            of the element. do only if iter > 0
c
      if ( iter .ne. 0 ) then
        call rstgp1_b( span, internal_energy, plastic_work,
     &                 local_work%urcs_blk_n1(1,7,gpn),
     &                 local_work%urcs_blk_n1(1,8,gpn),
     &                 local_work%det_j(1,gpn), local_work%dfn1, 2 )
        internal_energy = internal_energy * beta_fact *
     &                    local_work%weights(gpn)
        plastic_work    = plastic_work * beta_fact *
     &                    local_work%weights(gpn)
        local_work%block_energy       = internal_energy
        local_work%block_plastic_work = plastic_work
        if ( local_debug ) write (iout,9500) internal_energy
      end if
c
      return
c
 9500 format('  Internal energy inside of (rstgp1) = ',e16.6)
c
      end
c     ****************************************************************
c     *                                                              *
c     *                subroutine drive_01_update                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *              last modified : 06/24/12 RHD                    *
c     *                                                              *
c     *     this subroutine drives material model #1 to              *
c     *     update stresses and history for all elements in the      *
c     *     block for gauss point gpn                                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_01_update( gpn, props, lprops, iprops,
     &                            local_work, uddt, iout )
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
c
      implicit integer (a-z)
$add param_def
c
c                      parameter declarations
c

      real    props(mxelpr,*)
      logical lprops(mxelpr,*)
      integer iprops(mxelpr,*)
#dbl      double precision
#sgl      real
     &  uddt(mxvl,nstr)
$add include_sig_up
c
c                       locally defined variables
c
#dbl      double precision
#sgl      real
     &  dtime, gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero,  ddummy(1), gp_alpha, ymfgm, et
c
      logical geonl, local_debug, temperatures, segmental,
     &        temperatures_ref, fgm_enode_props
c
      data local_debug, zero / .false., 0.0 /
c
c           vectorized mises plasticity model with constant hardening
c           modulus. the model supports temperature dependence of
c           the elastic modulus, nu, hprime, and thermal
c           expansion alpha can vary. temperature dependent
c           properties enter through segmental curves.
c
      dtime             = local_work%dt
      span              = local_work%span
      felem             = local_work%felem
      type              = local_work%elem_type
      order             = local_work%int_order
      ngp               = local_work%num_int_points
      nnode             = local_work%num_enodes
      ndof              = local_work%num_enode_dof
      geonl             = local_work%geo_non_flg
      step              = local_work%step
      iter              = local_work%iter
      now_blk           = local_work%blk
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      segmental         = local_work%segmental
      number_points     = local_work%number_points
      curve_set         = local_work%curve_set_number
      fgm_enode_props   = local_work%fgm_enode_props
c
c          determine if the material elastic and properties are
c          described by segmental curve(s) and if they are temperature
c          or strain rate dependent [=0 no dependence, =1 temperature
c          dependent, =2 strain-rate dependent (not used here)]
c
      curve_type = -1
      if ( segmental ) call set_segmental_type( curve_set, curve_type,
     &                                          local_work%eps_curve )
c
c          determine if some material properties are specified using
c          values at model nodes to define fgms. interpolate values
c          at the current gauss point (overwrite the constant values).
c          the material cannot be both segmental & fgm (already
c          checked). interpolated values are the same at n and n+1.
c          at present only e, nu, tan_e, sig_yld and
c          isotropic alpha can be specified as fgm properties
c          interpolated from nodal values. note we have to
c          recompute h-prime after the correct e and tan_e are found.
c
      if ( fgm_enode_props ) then
c
          call  set_fgm_solid_props_for_block(
     &            span, felem, type, gpn, nnode,
     &            local_work%e_vec_n, local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 1,
     &            local_work%fgm_flags(1,1) )
          call  set_fgm_solid_props_for_block(
     &            span, felem, elem_type, gpn, nnode,
     &            local_work%nu_vec_n, local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 2,
     &            local_work%fgm_flags(1,2) )
          call  set_fgm_solid_props_for_block(
     &            span, felem, elem_type, gpn, nnode,
     &            local_work%alpha_vec_n(1,1), local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 3,
     &            local_work%fgm_flags(1,3) )
          call  set_fgm_solid_props_for_block(
     &            span, felem, elem_type, gpn, nnode,
     &            local_work%tan_e_vec(1), local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 6,
     &            local_work%fgm_flags(1,6) )
          call  set_fgm_solid_props_for_block(
     &            span, felem, elem_type, gpn, nnode,
     &            local_work%sigyld_vec(1), local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 7,
     &            local_work%fgm_flags(1,7) )
c
          do i = 1, span
            local_work%e_vec(i)  = local_work%e_vec_n(i)
            local_work%nu_vec(i) = local_work%nu_vec_n(i)
            gp_alpha  = local_work%alpha_vec_n(i,1)
            local_work%alpha_vec(i,1)   = gp_alpha
            local_work%alpha_vec(i,2)   = gp_alpha
            local_work%alpha_vec(i,3)   = gp_alpha
            local_work%alpha_vec(i,4)   = zero
            local_work%alpha_vec(i,5)   = zero
            local_work%alpha_vec(i,6)   = zero
            local_work%alpha_vec_n(i,2) = gp_alpha
            local_work%alpha_vec_n(i,3) = gp_alpha
            local_work%alpha_vec_n(i,4) = zero
            local_work%alpha_vec_n(i,5) = zero
            local_work%alpha_vec_n(i,6) = zero
            ymfgm = local_work%e_vec_n(i)
            et    = local_work%tan_e_vec(i)
            local_work%h_vec(i) = (ymfgm*et)/(ymfgm-et)
          end do
c
      end if
c
c           get increment of temperature at gauss point for elements
c           in the block, the temperature at end of step and the
c           reference temperature.
c
      call gauss_pt_temps(
     &        local_work%dtemps_node_blk, gpn, type, span, order,
     &        nnode, gp_dtemps, local_work%temps_node_blk,
     &        gp_temps, temperatures, local_work%temps_node_to_process,
     &        temperatures_ref, local_work%temps_ref_node_blk,
     &        gp_rtemps )
c
c          get temperature dependent young's modulus, poisson's
c          ratio, uniaxial yield stress and constant plastic
c          modulus for the temperature at n and n+1. Some
c          properties are not for this model but need to be passed
c          to satisfy syntax of call.
c
        if ( curve_type .eq. 0 .or. curve_type .eq. 1 ) then
          call set_up_segmental( span, gp_temps, local_work%e_vec,
     &        local_work%nu_vec, local_work%alpha_vec,
     &        local_work%e_vec_n, local_work%nu_vec_n,
     &        local_work%alpha_vec_n,
     &        local_work%gp_sig_0_vec,
     &        local_work%gp_h_u_vec,
     &        local_work%gp_beta_u_vec,
     &        local_work%gp_delta_u_vec,
     &        local_work%gp_sig_0_vec_n,
     &        local_work%gp_h_u_vec_n,
     &        local_work%gp_beta_u_vec_n,
     &        local_work%gp_delta_u_vec_n,
     &        gp_dtemps,
     &        ddummy, gpn, mxvl )
c
         call set_up_h_prime( span, local_work%h_vec,
     &                        local_work%sigyld_vec, felem )
        end if
c
c          subtract out the thermal strain increment from uddt (the
c          strain increment for step)
c
      if ( temperatures )
     &               call gp_temp_eps( span, uddt,
     &                    local_work%alpha_vec, gp_dtemps ,
     &                    gp_temps, gp_rtemps,
     &                    local_work%alpha_vec_n )
c
c          for iter = 0 we just use the stored [Dt] to compute
c          stress increment. These stresses are used to dompute
c          internal forces for applied nodal displacements and
c          temperatures - material history is unaffected by
c          iter = 0
c
      if( iter .eq. 0 ) then
        igpn = gpn
        call recstr_cep_uddt_for_block( mxvl, span,
     &    gbl_cep_blocks(now_blk)%vector,
     &    uddt, local_work%urcs_blk_n(1,1,gpn),
     &    local_work%urcs_blk_n1(1,1,gpn), 1, 21, igpn )
      else
       call mm01( span, felem, gpn, step, iter, local_work%e_vec,
     &           local_work%nu_vec, local_work%beta_vec,
     &           local_work%h_vec, local_work%lnelas_vec,
     &           local_work%sigyld_vec,
     &           local_work%urcs_blk_n(1,1,gpn),
     &           local_work%urcs_blk_n1(1,1,gpn),
     &           uddt, local_work%elem_hist(1,1,gpn),
     &           local_work%elem_hist1(1,1,gpn),
     &           local_work%rtse(1,1,gpn), gp_dtemps,
     &           local_work%e_vec_n, local_work%nu_vec_n  )
      endif
c
      return
c
 9600 format(' >> Fatal Error: loop to iterate on plastic-strain'
     & /,    '                  rates did not converge. Job aborted')
 9610 format(' >> rate iterations to converge: ',i3 )
c
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_02_update                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/24/12 rhd               *
c     *                                                              *
c     *     this subroutine drives material model 02 to              *
c     *     update stresses and history for all elements in the      *
c     *     block for an integration point                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_02_update( gpn, props, lprops, iprops,
     &                            local_work, uddt, iout )
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
c
      implicit integer (a-z)
$add param_def
c
c                      parameter declarations
c
      real    props(mxelpr,*)
      logical lprops(mxelpr,*)
      integer iprops(mxelpr,*)
#dbl      double precision
#sgl      real
     &  uddt(mxvl,nstr)
$add include_sig_up
c
c
c                       locally defined variables
c
#dbl      double precision
#sgl      real
     &  gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  ddummy(1), zero, gp_alpha, ddtse(mxvl,6), nowtemp
      logical signal_flag, fgm_enode_props, local_debug,
     &        temperatures, temperatures_ref
      data local_debug, zero / .false., 0.0 /
c
c
c          deformation plasticity model. properties are invariant of
c          temperature and loading rate. properties may vary spatially
c          for fgms - values here are interploated at the current
c          integration point.
c
c
      span              = local_work%span
      felem             = local_work%felem
      step              = local_work%step
      iter              = local_work%iter
      now_blk           = local_work%blk
      type              = local_work%elem_type
      order             = local_work%int_order
      nnode             = local_work%num_enodes
      signal_flag       = local_work%signal_flag
      fgm_enode_props   = local_work%fgm_enode_props
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
c
c          for fgms, interpolate values of material properties
c          at the current gauss point (overwrite the constant values).
c          at present only e, nu, sig_yld, n_power and
c          isotropic alpha can be specified
c          as fgm properties interpolated from nodal values.
c
      if ( fgm_enode_props ) then
c
          call  set_fgm_solid_props_for_block(
     &            span, felem, type, gpn, nnode,
     &            local_work%e_vec, local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 1,
     &            local_work%fgm_flags(1,1) )
          call  set_fgm_solid_props_for_block(
     &            span, felem, elem_type, gpn, nnode,
     &            local_work%nu_vec, local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 2,
     &            local_work%fgm_flags(1,2) )
          call  set_fgm_solid_props_for_block(
     &            span, felem, elem_type, gpn, nnode,
     &            local_work%alpha_vec(1,1), local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 3,
     &            local_work%fgm_flags(1,3) )
          call  set_fgm_solid_props_for_block(
     &            span, felem, elem_type, gpn, nnode,
     &            local_work%sigyld_vec(1), local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 7,
     &            local_work%fgm_flags(1,7) )
          call  set_fgm_solid_props_for_block(
     &            span, felem, elem_type, gpn, nnode,
     &            local_work%n_power_vec(1), local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 8,
     &            local_work%fgm_flags(1,8) )
          do i = 1, span
            gp_alpha  = local_work%alpha_vec(i,1)
            local_work%alpha_vec(i,1)   = gp_alpha
            local_work%alpha_vec(i,2)   = gp_alpha
            local_work%alpha_vec(i,3)   = gp_alpha
            local_work%alpha_vec(i,4)   = zero
            local_work%alpha_vec(i,5)   = zero
            local_work%alpha_vec(i,6)   = zero
          end do
      end if
c
c           get increment of temperature at gauss point for elements
c           in the block, the temperature at end of step and the
c           reference temperature.
c
      call gauss_pt_temps(
     &        local_work%dtemps_node_blk, gpn, type, span, order,
     &        nnode, gp_dtemps, local_work%temps_node_blk,
     &        gp_temps, temperatures, local_work%temps_node_to_process,
     &        temperatures_ref, local_work%temps_ref_node_blk,
     &        gp_rtemps )
c
c          for iter = 0, use the stored [Dt] from last stiffness update
c          to compute stress @ n+1 = stress @ n + [Dt]*uddt. Remove
c          temperature increment from uddt (deps) if needed.
c          iter = 0 is the process to compute stresses for equivalent loads
c          generation from imposed displacement and temperature loading.
c
      if( iter .gt. 0 ) go to 1000
c
      if( temperatures )
     &    call gp_temp_eps( span, uddt, local_work%alpha_vec, gp_dtemps,
     &                    gp_temps, gp_rtemps, local_work%alpha_vec )
      igpn = gpn
      call recstr_cep_uddt_for_block( mxvl, span,
     &      gbl_cep_blocks(now_blk)%vector,
     &      uddt, local_work%urcs_blk_n(1,1,gpn),
     &      local_work%urcs_blk_n1(1,1,gpn), 1, 21, igpn )
      return
c
c          process iter > 0. mm02 uses total strains adjusted for total
c          temperature strains to compute stresses @ n+1.
c
 1000 continue
      if( .not. temperatures ) then
         call mm02( step, iter, felem, gpn, local_work%e_vec,
     &           local_work%nu_vec, local_work%sigyld_vec,
     &           local_work%n_power_vec,local_work%urcs_blk_n(1,1,gpn),
     &           local_work%urcs_blk_n1(1,1,gpn),
     &           local_work%ddtse(1,1,gpn),  ! total strain @ n+1
     &           local_work%elem_hist(1,1,gpn),
     &           local_work%elem_hist1(1,1,gpn), span, iout,
     &           signal_flag )
         return
      end if
c
      do i = 1, span
           nowtemp = gp_temps(i) - gp_rtemps(i)
           ddtse(i,1) = local_work%ddtse(i,1,gpn) -
     &                  local_work%alpha_vec(i,1)*nowtemp
           ddtse(i,2) = local_work%ddtse(i,2,gpn) -
     &                  local_work%alpha_vec(i,2)*nowtemp
           ddtse(i,3) = local_work%ddtse(i,3,gpn) -
     &                  local_work%alpha_vec(i,3)*nowtemp
           ddtse(i,4) = local_work%ddtse(i,4,gpn) -
     &                  local_work%alpha_vec(i,4)*nowtemp
           ddtse(i,5) = local_work%ddtse(i,5,gpn) -
     &                  local_work%alpha_vec(i,5)*nowtemp
           ddtse(i,6) = local_work%ddtse(i,6,gpn) -
     &                  local_work%alpha_vec(i,6)*nowtemp
       end do
       call mm02( step, iter, felem, gpn, local_work%e_vec,
     &           local_work%nu_vec, local_work%sigyld_vec,
     &           local_work%n_power_vec,
     &           local_work%urcs_blk_n(1,1,gpn),
     &           local_work%urcs_blk_n1(1,1,gpn),
     &           ddtse, local_work%elem_hist(1,1,gpn),
     &           local_work%elem_hist1(1,1,gpn), span, iout,
     &           signal_flag )
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine drive_03_update                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *             last modified : 03/18/04 rhd                     *
c     *                                                              *
c     *     this subroutine drives material model 03 to              *
c     *     update stresses and history for all elements in the      *
c     *     for gauss point gpn                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_03_update( gpn, props, lprops, iprops,
     &                            local_work, uddt, iout )
      use segmental_curves, only : max_seg_points, max_seg_curves,
     &                             now_blk_relem, sigma_curve_min_values
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
c
      implicit integer (a-z)
$add param_def
c
c                      parameter declarations
c
      real    props(mxelpr,*)
      logical lprops(mxelpr,*)
      integer iprops(mxelpr,*)
#dbl      double precision
#sgl      real
     &  uddt(mxvl,nstr)
$add include_sig_up
c
c                       locally defined variables
c
#dbl      double precision
#sgl      real
     &  dtime, internal_energy, beta_fact, eps_bbar, ddt(mxvl,nstr),
     &  q(mxvl,nstr,nstr), stress_n(nstrs,mxvl),
     &  stress_n1(nstrs,mxvl), p_trial(mxvl), q_trial(mxvl), dfn1(mxvl),
     &  yld_func(mxvl), step_scale_fact, plastic_work, gp_temps(mxvl),
     &  gp_dtemps(mxvl), plastic_eps_rates(mxvl), gp_rtemps(mxvl),
     &  zero, gp_alpha, et, ymfgm, copy_sigyld_vec(mxvl),
     &  trans_factor, curve_min_value
c
      logical null_point(mxvl), cut_step, process_block,
     &        adaptive, geonl, bbar, material_cut_step,
     &        local_debug, signal_flag, adaptive_flag,
     &        power_law, temperatures, allow_cut, segmental,
     &        model_update, temperatures_ref, fgm_enode_props
c
      type :: arguments
#r60        sequence
        integer :: iter, abs_element, relem, ipoint, iout
        logical :: allow_cut, segmental, power_law,
     &             rate_depend_segmental, signal_flag, cut_step
#sgl        real :: dtime, step_scale_fact
#dbl        double precision :: dtime, step_scale_fact
      end type
c
      type (arguments) ::args

      data zero, local_debug / 0.0, .false. /
      data trans_factor / 0.95 /
c
      dtime             = local_work%dt
      internal_energy   = local_work%block_energy
      plastic_work      = local_work%block_plastic_work
      beta_fact         = local_work%beta_fact
      span              = local_work%span
      felem             = local_work%felem
      type              = local_work%elem_type
      order             = local_work%int_order
      ngp               = local_work%num_int_points
      nnode             = local_work%num_enodes
      ndof              = local_work%num_enode_dof
      geonl             = local_work%geo_non_flg
      step              = local_work%step
      iter              = local_work%iter
      now_blk           = local_work%blk
      bbar              = local_work%bbar_flg
      mat_type          = local_work%mat_type
      material_cut_step = local_work%material_cut_step
      signal_flag       = local_work%signal_flag
      adaptive_flag     = local_work%adaptive_flag
      eps_bbar          = local_work%eps_bbar
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      allow_cut         = local_work%allow_cut
      segmental         = local_work%segmental
      number_points     = local_work%number_points
      curve_set         = local_work%curve_set_number
      power_law         = local_work%power_law
      step_scale_fact   = local_work%step_scale_fact
      adaptive          = adaptive_flag .and. step .gt. 1
      fgm_enode_props   = local_work%fgm_enode_props
      hist_size_for_blk = local_work%hist_size_for_blk

c
      curve_type = -1
      if ( segmental ) call set_segmental_type( curve_set, curve_type,
     &                                          local_work%eps_curve )
c
c          for fgms, interpolate values of material properties
c          at the current gauss point (overwrite the constant values).
c          the material cannot be both segmental & fgm (already
c          checked). interpolated values are the same at n and n+1.
c          at present only e, nu, tan_e, sig_yld, n_power and
c          isotropic alpha can be specified
c          as fgm properties interpoalted from nodal values. note
c          we have to recompute h-prime after the correct e is found.
c
      if ( fgm_enode_props ) then
c
          call  set_fgm_solid_props_for_block(
     &            span, felem, type, gpn, nnode,
     &            local_work%e_vec_n, local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 1,
     &            local_work%fgm_flags(1,1) )
          call  set_fgm_solid_props_for_block(
     &            span, felem, elem_type, gpn, nnode,
     &            local_work%nu_vec_n, local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 2,
     &            local_work%fgm_flags(1,2) )
          call  set_fgm_solid_props_for_block(
     &            span, felem, elem_type, gpn, nnode,
     &            local_work%alpha_vec_n(1,1), local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 3,
     &            local_work%fgm_flags(1,3) )
          call  set_fgm_solid_props_for_block(
     &            span, felem, elem_type, gpn, nnode,
     &            local_work%tan_e_vec(1), local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 6,
     &            local_work%fgm_flags(1,6) )
          call  set_fgm_solid_props_for_block(
     &            span, felem, elem_type, gpn, nnode,
     &            local_work%sigyld_vec(1), local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 7,
     &            local_work%fgm_flags(1,7) )
          call  set_fgm_solid_props_for_block(
     &            span, felem, elem_type, gpn, nnode,
     &            local_work%n_power_vec(1), local_work%shape(1,gpn),
     &            local_work%enode_mat_props, 8,
     &            local_work%fgm_flags(1,8) )
c
          do i = 1, span
            local_work%e_vec(i)         = local_work%e_vec_n(i)
            local_work%nu_vec(i)        = local_work%nu_vec_n(i)
            gp_alpha                    = local_work%alpha_vec_n(i,1)
            local_work%alpha_vec(i,1)   = gp_alpha
            local_work%alpha_vec(i,2)   = gp_alpha
            local_work%alpha_vec(i,3)   = gp_alpha
            local_work%alpha_vec(i,4)   = zero
            local_work%alpha_vec(i,5)   = zero
            local_work%alpha_vec(i,6)   = zero
            local_work%alpha_vec_n(i,2) = gp_alpha
            local_work%alpha_vec_n(i,3) = gp_alpha
            local_work%alpha_vec_n(i,4) = zero
            local_work%alpha_vec_n(i,5) = zero
            local_work%alpha_vec_n(i,6) = zero
            ymfgm                       = local_work%e_vec_n(i)
            et                          = local_work%tan_e_vec(i)
            local_work%h_vec(i)         = (ymfgm*et)/(ymfgm-et)
          end do
c
      end if
c
c          get temperature increment at the gauss point for all elements
c          of block. compute the temperature at the end of the step
c          at the gauss point for all elements of block.
c
      call gauss_pt_temps(
     &        local_work%dtemps_node_blk, gpn, type, span, order,
     &        nnode, gp_dtemps, local_work%temps_node_blk,
     &        gp_temps, temperatures,
     &        local_work%temps_node_to_process,
     &        temperatures_ref, local_work%temps_ref_node_blk,
     &        gp_rtemps )
c
c          get temperature dependent elastic and flow properites.
c          build local copy of stress vs. plastic strain curve
c          for use at this gauss point for elements in block.
c          set yield stress at zero plastic strain.  Some
c          properties are not for this model but need to be passed
c          to satisfy syntax of call.
c
      if ( curve_type .eq. 0 .or. curve_type .eq. 1 ) then
         call set_up_segmental( span, gp_temps, local_work%e_vec,
     &        local_work%nu_vec, local_work%alpha_vec,
     &        local_work%e_vec_n, local_work%nu_vec_n,
     &        local_work%alpha_vec_n,
     &        local_work%gp_sig_0_vec,
     &        local_work%gp_h_u_vec,
     &        local_work%gp_beta_u_vec,
     &        local_work%gp_delta_u_vec,
     &        local_work%gp_sig_0_vec_n,
     &        local_work%gp_h_u_vec_n,
     &        local_work%gp_beta_u_vec_n,
     &        local_work%gp_delta_u_vec_n,
     &        gp_dtemps, plastic_eps_rates,
     &        gpn, mxvl )
c
         call set_up_h_prime( span, local_work%h_vec,
     &                        local_work%sigyld_vec, felem )
      end if
c
c          subtract out the thermal strain increment from uddt.
c
      if ( temperatures )
     &   call gp_temp_eps( span, uddt, local_work%alpha_vec,
     &                     gp_dtemps,
     &                     gp_temps, gp_rtemps,
     &                     local_work%alpha_vec_n )

      if( iter .eq. 0 ) then
        igpn = gpn
        call recstr_cep_uddt_for_block( mxvl, span,
     &      gbl_cep_blocks(now_blk)%vector(1),
     &      uddt, local_work%urcs_blk_n(1,1,gpn),
     &      local_work%urcs_blk_n1(1,1,gpn), 1, 21, igpn )
        return
       end if
c
c
c          for rate dependent response defined by a set of segmental
c          curves, set up the internal tables for interpolation in
c          the seg_curve module. we also load the initial stress,
c          initial plastic modulus and curve min stress value to
c          be the values for lowest defined strain rate.  Some
c          properties are not for this model but need to be passed
c          to satisfy syntax of call.
c
      if ( curve_type .eq. 2 ) then
         plastic_eps_rates(1:span) = zero
         call set_up_segmental( span, gp_temps, local_work%e_vec,
     &        local_work%nu_vec, local_work%alpha_vec,
     &        local_work%e_vec_n, local_work%nu_vec_n,
     &        local_work%alpha_vec_n,
     &        local_work%gp_sig_0_vec,
     &        local_work%gp_h_u_vec,
     &        local_work%gp_beta_u_vec,
     &        local_work%gp_delta_u_vec,
     &        local_work%gp_sig_0_vec_n,
     &        local_work%gp_h_u_vec_n,
     &        local_work%gp_beta_u_vec_n,
     &        local_work%gp_delta_u_vec_n,
     &        gp_dtemps,
     &        plastic_eps_rates, gpn, mxvl )
c
         call set_up_h_prime( span, local_work%h_vec,
     &                        local_work%sigyld_vec, felem )
      end if
c

c --------------------------------------------------------------------
c
c
c         general mises/gurson flow theory model.
c         we process element block in two phases.
c         in phase (1) we compute the trial elastic
c         stress state and evaluate the yield function
c         using vectorized code. elements in block
c         with this gp nonlinear are marked with a
c         logical. for elements remaining linear, we
c         are done and if all elements in block are
c         linear, we have finished block (the energy
c         densities for linear elements us updated)
c
c         in phase (2) we call the nonlinear stress update
c         routine at this gauss point only for elements
c         which are nonlinear. trial elastic state
c         info is passed down to minimize redundant
c         computations.
c
c         block data structures below are generated to
c         enable passing contiguous data in vectors to
c         model. the loops are fully unrolled to force
c         vectorization over span and not 1:6 or 1:7.
c
c         all data about segmental stress-strain curves are
c         in the segmental curves module. this saves
c         passing a large number of arguments down through
c         many levels of subroutines to where curves are needed
c         in stress update code.
c
c --------------------------------------------------------------------
c
c                set up local stress array at state 'n'
c
      numrows_stress = nstrs
      call rstgp1_d( span, numrows_stress, stress_n,
     &               local_work%urcs_blk_n(1,1,gpn), mxvl )
c
c                set up processing for all elements in block at
c                this gauss point. elements can be linear,
c                nonlinear or null (all strains = 0). set logical
c                if we have to process any elements in block.
c
c                build copy of yield stress at this gp for each element
c                in block. for power law hardening, this value must
c                be reduced to accomodate the cubic transition.
c                the mm03 routines all assume this reduction is
c                done here. we can't change local_work since the
c                next gp could be using the same vector.
c
      copy_sigyld_vec(1:span) = local_work%sigyld_vec(1:span)
      if ( .not. segmental ) then
        do i = 1, span
         if ( local_work%n_power_vec(i) .gt. zero )
     &    copy_sigyld_vec(i) = copy_sigyld_vec(i) * trans_factor
        end do
      end if
c
      call mm03p(
     &       step, iter, span, gpn, uddt, local_work%e_vec,
     &       local_work%nu_vec, copy_sigyld_vec,
     &       local_work%f0_vec, local_work%q1_vec, local_work%q2_vec,
     &       local_work%q3_vec, local_work%n_power_vec,
     &       local_work%h_vec, null_point, local_work%nuc_vec,
     &       stress_n, stress_n1, local_work%rtse,
     &       p_trial, q_trial, local_work%elem_hist(1,1,1),
     &       local_work%elem_hist1(1,1,1), yld_func,
     &       local_work%nonlinear_flag,
     &       process_block, iout, segmental, curve_type, felem,
     &       local_work%e_vec_n, local_work%nu_vec_n,
     &       hist_size_for_blk  )
c
      if ( .not. process_block ) go to 1000
c
c                 we use the args derived type here to reduce the number of
c                 parameters passed to mm03 since it is called many times.
c                 the variable now_blk_relem must be set as well (in
c                 segmental curve module)
c
      args%iter                  = iter
      args%ipoint                = gpn
      args%dtime                 = dtime
      args%iout                  = iout
      args%allow_cut             = allow_cut
      args%cut_step              = .false.
      args%signal_flag           = signal_flag
      args%segmental             = segmental
      args%power_law             = power_law
      args%step_scale_fact       = step_scale_fact
      args%rate_depend_segmental = curve_type .eq. 2
c
      do i = 1, span
        if ( null_point(i) ) cycle
        if ( .not. local_work%nonlinear_flag(i) ) cycle
        now_blk_relem = i
        args%abs_element = felem+i-1
        args%relem = i
        curve_min_value = zero
        If( segmental ) curve_min_value = sigma_curve_min_values(i)
        call mm03(
     &         args, local_work%e_vec(i),
     &         local_work%nu_vec(i), local_work%f0_vec(i),
     &         local_work%eps_ref_vec(i), copy_sigyld_vec(i),
     &         local_work%m_power_vec(i),
     &         local_work%n_power_vec(i), local_work%h_vec(i),
     &         local_work%q1_vec(i), local_work%q2_vec(i),
     &         local_work%q3_vec(i), local_work%nuc_vec(i),
     &         local_work%nuc_s_n_vec(i),
     &         local_work%nuc_e_n_vec(i),local_work%nuc_f_n_vec(i),
     &         stress_n(1,i), stress_n1(1,i), local_work%rtse, uddt,
     &         local_work%elem_hist(1,1,1),
     &         local_work%elem_hist1(1,1,1),
     &         yld_func(i), p_trial(i), q_trial(i), mxvl,
     &         hist_size_for_blk, curve_min_value, span )
        material_cut_step = material_cut_step .or. args%cut_step
      end do
c
      local_work%material_cut_step = material_cut_step
c
c                save stresses at n+1 into block data structure.
c                step cut flag for block if material model
c                failed to converge and it wants a load step
c                size reduction
c
 1000 continue
      call rstgp1_c( span, numrows_stress, stress_n1,
     &               local_work%urcs_blk_n1(1,1,gpn), mxvl )
c
c
      return
c
 9000 format('>>> Fatal Error: the nonlinear elastic material',
     &     /,'                 is not compatible with large',
     &     /,'                 displacement elements.',
     &     /,'                 job terminated....' )
 9100 format(i5,7f15.6,/,5x,7f15.6,/,5x,6f15.6,/,5x,6f15.6)
 9110 format(1x,'Elem    /',20('-'),
     &        ' unrot. Cauchy @ n, unrot. Cauchy @ n+1,',
     &        ' ddt, uddt', 20('-'),'/')
 9500 format('  Internal energy inside of (rstgp1)    = ',e16.6,
     &     /,'  Plasstic work inside of (rstgp1)      = ',e16.6)
 9600 format(' >> Fatal Error: loop to iterate on plastic-strain'
     & /,    '                  rates did not converge. Job aborted')
 9610 format(' >> rate iterations to converge: ',i3 )
 9820 format(///,
     &       '>> FATAL ERROR: strain computation routines requested',
     &     /,'                immediate step size reduction due to',
     &     /,'                invalid determinant of a deformation',
     &     /,'                jacobian. the user has not allowed step',
     &     /,'                size reductions. analysis terminated.',
     &     /// )
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_04_update                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 03/13/14 rhd               *
c     *                                                              *
c     *     this subroutine drives material model 04 to              *
c     *     update stresses and history for all elements in the      *
c     *     for gauss point gpn                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_04_update( gpn, props, lprops, iprops,
     &                            local_work, uddt, iout )
      use segmental_curves, only : max_seg_points
      use elem_block_data, only  : gbl_cep_blocks => cep_blocks
c
      implicit integer (a-z)
$add param_def
c
c                      parameter declarations
c
      real    props(mxelpr,*)
      logical lprops(mxelpr,*)
      integer iprops(mxelpr,*)
#dbl      double precision
#sgl      real
     &  uddt(mxvl,nstr)
$add include_sig_up
c
c                       locally defined variables
c
#dbl      double precision
#sgl      real
     &  time_n, dtime, ddummy(mxvl)
      integer idummy(mxvl)

      logical fgm_enode_props, nonlocal, temperatures,
     &        temperatures_ref
c
      span              = local_work%span
      felem             = local_work%felem
      step              = local_work%step
      iter              = local_work%iter
      now_blk           = local_work%blk
      nnode             = local_work%num_enodes
      fgm_enode_props   = local_work%fgm_enode_props
c
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
c
c          for iter = 0 we just use the stored [Dt] to compute
c          stress increment. These stresses are used to dompute
c          internal forces for applied nodal displacements and
c          temperatures - material history is unaffected by
c          iter = 0

      if( iter .eq. 0 ) then
        igpn = gpn
        call recstr_cep_uddt_for_block( mxvl, span,
     &      gbl_cep_blocks(now_blk)%vector,
     &      uddt, local_work%urcs_blk_n(1,1,gpn),
     &      local_work%urcs_blk_n1(1,1,gpn), 2, 6, igpn )
        return
       end if

       if( gpn .eq. 1 ) local_work%elem_hist1 = 0.0
c
c                linear and non-linear cohesive model
c
      if( nonlocal ) then
        call mm04(
     g    step, iter, span, felem, gpn, iout, imxvl, time_n, dtime,
     1    nonlocal, knumthreads, kthread,
     2    local_work%intf_prp_block,
     3    local_work%cohes_type,
     4    local_work%urcs_blk_n(1,1,gpn),
     5    local_work%urcs_blk_n1(1,1,gpn),
     6    local_work%ddtse(1,1,gpn), uddt(1,1),
     7    local_work%elem_hist(1,1,gpn),
     8    local_work%elem_hist1(1,1,gpn),
     9    local_work%cohes_temp_ref(1),
     h    local_work%cohes_dtemp(1),
     i    local_work%cohes_temp_n(1),
     a    local_work%top_surf_solid_elements(1),
     b    local_work%bott_surf_solid_elements(1),
     c    local_work%top_surf_solid_stresses_n(1,1),
     d    local_work%bott_surf_solid_stresses_n(1,1),
     e    local_work%top_surf_solid_eps_n(1,1),
     f    local_work%bott_surf_solid_eps_n(1,1),
     g    local_work%nonlocal_stvals_top_n(1,1),
     h    local_work%nonlocal_stvals_bott_n(1,1),
     i    local_work%top_solid_matl(1),
     j    local_work%bott_solid_matl(1) )
      else
        call mm04(
     g    step, iter, span, felem, gpn, iout, imxvl, time_n, dtime,
     1    nonlocal, knumthreads, kthread,
     2    local_work%intf_prp_block,
     3    local_work%cohes_type,
     4    local_work%urcs_blk_n(1,1,gpn),
     5    local_work%urcs_blk_n1(1,1,gpn),
     6    local_work%ddtse(1,1,gpn), uddt(1,1),
     7    local_work%elem_hist(1,1,gpn),
     8    local_work%elem_hist1(1,1,gpn),
     9    local_work%cohes_temp_ref(1),
     h    local_work%cohes_dtemp(1),
     i    local_work%cohes_temp_n(1),
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
c     *                 subroutine drive_05_update                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 07/31/2011                 *
c     *                                                              *
c     *     this subroutine drives material model 05 to              *
c     *     update stresses and history for all elements in the      *
c     *     for gauss point gpn                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_05_update( gpn, props, lprops, iprops,
     &                            local_work, uddt, iout )
      use main_data, only : matprp, lmtprp
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
c
      implicit integer (a-z)
$add param_def
c
c                      parameter declarations
c
      real    props(mxelpr,*)
      logical lprops(mxelpr,*)
      integer iprops(mxelpr,*)
#dbl      double precision
#sgl      real
     &  uddt(mxvl,nstr)
$add include_sig_up
c
c
c                       locally defined variables
c
#dbl      double precision
#sgl      real
     &  gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero, gp_alpha, dtime, sig_tol,
     &  nh_sigma_0_vec(mxvl), nh_q_u_vec(mxvl), nh_b_u_vec(mxvl),
     &  nh_h_u_vec(mxvl), nh_gamma_u_vec(mxvl), gp_tau_vec(mxvl)

c
      logical signal_flag, local_debug, temperatures,
     &        temperatures_ref, adaptive_possible, cut_step_size_now,
     &        segmental, nonlin_hard, generalized_pl
      data local_debug, zero / .false., 0.0 /
c
c                  NOTE:  at present, all elements in the block must be
c                         same cyclic material defined by the user.
c                         this restriction can be removed by tracking
c                         thru and allowing matnum to vary by element in
c                         the block
c
      matnum            = local_work%matnum
      dtime             = local_work%dt
      span              = local_work%span
      felem             = local_work%felem
      step              = local_work%step
      iter              = local_work%iter
      now_blk           = local_work%blk
      type              = local_work%elem_type
      order             = local_work%int_order
      nnode             = local_work%num_enodes
      signal_flag       = local_work%signal_flag
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      segmental         = local_work%segmental
      number_points     = local_work%number_points
      curve_set         = local_work%curve_set_number
      adaptive_possible = local_work%adaptive_flag .and.
     &                    step .gt. 1
      cut_step_size_now = .false.
      hist_size_for_blk = local_work%hist_size_for_blk
c
c          determine if the material elastic and properties are
c          temperature dependent. The temperature dependent props
c          are e, nu, alpha, gp_sigma_0, gp_h_u, gp_beta_u,
c          gp_delta_u. Points on a curve are defined in to
c          data to satisfy all error checks but are not used here.
c
      curve_type = -1
      if ( segmental ) call set_segmental_type( curve_set, curve_type,
     &                                          local_work%eps_curve )
c
c           get increment of temperature at gauss point for elements
c           in the block, the temperature at end of step and the
c           reference temperature.
c
      call gauss_pt_temps(
     &        local_work%dtemps_node_blk, gpn, type, span, order,
     &        nnode, gp_dtemps, local_work%temps_node_blk,
     &        gp_temps, temperatures, local_work%temps_node_to_process,
     &        temperatures_ref, local_work%temps_ref_node_blk,
     &        gp_rtemps )
c
c          get temperature dependent young's modulus, poisson's
c          ratio, alpha, gp_sigma_0, gp_h_u, gp_beta_u, gp_delta_u
c          for the temperature at n and n+1
c
      if ( curve_type .eq. 1 ) then
         call set_up_segmental( span, gp_temps, local_work%e_vec,
     &        local_work%nu_vec, local_work%alpha_vec,
     &        local_work%e_vec_n, local_work%nu_vec_n,
     &        local_work%alpha_vec_n,
     &        local_work%gp_sig_0_vec,
     &        local_work%gp_h_u_vec,
     &        local_work%gp_beta_u_vec,
     &        local_work%gp_delta_u_vec,
     &        local_work%gp_sig_0_vec_n,
     &        local_work%gp_h_u_vec_n,
     &        local_work%gp_beta_u_vec_n,
     &        local_work%gp_delta_u_vec_n,
     &        gp_dtemps,
     &        ddummy, gpn, mxvl )
      end if
c
c            subtract out the thermal strain increment from uddt (the
c            strain increment for step)
c
      if ( temperatures ) then
        call gp_temp_eps( span, uddt, local_work%alpha_vec, gp_dtemps,
     &                    gp_temps, gp_rtemps, local_work%alpha_vec_n )
      end if
c
c          for iter = 0 we just use the stored [Dt] to compute
c          stress increment. These stresses are used to dompute
c          internal forces for applied nodal displacements and
c          temperatures - material history is unaffected by
c          iter = 0

      if( iter .eq. 0 ) then
        igpn = gpn
        call recstr_cep_uddt_for_block( mxvl, span,
     &      gbl_cep_blocks(now_blk)%vector,
     &      uddt, local_work%urcs_blk_n(1,1,gpn),
     &      local_work%urcs_blk_n1(1,1,gpn), 1, 21, igpn )
        return
       end if
c
c
c            build local data vectors for nonlinear_hardening option
c            to maintain consistency with generalized_plasticity option.
c            the gp option has possibly temperature dependent properties.
c            the nh option could have e, nu and alpha defined as temperature
c            dependent by the user thru curves (but we really don't show that
c            option in the manual.
c
      do i = 1, span
        nh_sigma_0_vec(i) = local_work%sigyld_vec(i)
        nh_q_u_vec(i)     = matprp(55,matnum)
        nh_b_u_vec(i)     = matprp(56,matnum)
        nh_h_u_vec(i)     = matprp(57,matnum)
        nh_gamma_u_vec(i) = matprp(59,matnum)
      end do
c
      nonlin_hard      = matprp(58,matnum) .gt. zero
      generalized_pl   = matprp(58,matnum) .lt. zero
      sig_tol          = matprp(60,matnum)
c
c            build remaining local data vectors generalized_plasticity.
c
      do i = 1, span
        gp_tau_vec(i) = matprp(56,matnum)
      end do
c
c
c            available data for passing to cyclic model
c            ------------------------------------------
c
c   (*) updatable by material model
c
c     step              : current load step number
c     iter              : current newton iteration number
c     felem             : first element of the current block
c     gpn               : gauss point number being processed for block
c     mxvl              : maximum no. elements per block
c     hist_size_for_blk : number of history words per gauss point
c     nstrs             : number of stress terms (6 + extras)
c     nstr              : number of incremental strain terms (6)
c     span              : number of elements in current block
c     iout              : write messates to this device number
c     segmental         : .true. if user specified temperature dependent
c                         properties via curves for material. At present this
c                         is fully implemented only for generalized_plasticity.
c     nonlin_hard       : .true. if nonlinear_hardening option selected
c     generalized_pl    : .true. if generalized_plasticity option selected
c     signal_flag       : user wants notification messages for key
c                         events in material response
c     adaptive_possible : .true. if the material model may request
c                         immediate reduction of global load step size.
c                         no stress-histroy update required
c (*) cut_step_size_now : set .true. by material model if it wants immediate
c                         reduction of global load step size.
c                         no stress-history update required
c
c     uddt              : current estimate of strain increment over the
c                         load step (minus increment of thermal strain)
c
c     In local_work: same for both nonlinear_hardening (FA) and generalized_plasticity
c     -------------
c (*) rtse :            : trial elastic stress vector at n+1 to be used later by
c                         consistent tangent routine for model
c     urcs_blk_n        : stresses at start of load step (n) for all
c                         elements in block for this gauss point
c (*) urcs_blk_n1       : stresses at end of load step (n+1) for all
c                         elements in block for this gauss point
c     elem_hist         : history values at start of load step (n) for all
c                         elements in block for this gauss point
c (*) elem_hist1        : history values at end of load step (n+1) for all
c                         elements in block for this gauss point
c     e_vec             : Young's modulus for each element in block at n+1
c     nu_vec            : Poisson's ratio for each element in block at n+1
c
c                            props for each element in block at this gauss pt
c                            for gp_model. use for temp. dependent and
c                            temp. invariant cases
c
c     gp_sig_0_vec      : gp model. yield stress at n+1
c     gp_h_u_vec        : gp model. h_u at n+1
c     gp_beta_u_vec     : gp model. beta_u at n+1
c     gp_delta_u_vec    : gp model. delta_u at n+1
c     gp_sig_0_vec_n    : gp model. yield stress at n
c     gp_h_u_vec_n      : gp model. h_u at n
c     gp_beta_u_vec_n   : gp model. beta_u at n
c     gp_delta_u_vec_n  : gp model. delta_u at n
c
c     mm05_props(mxvl,10)    : values loaded from matprp by rknstr.f
c                              can we delete after James new code ???
c
c
c     In local space for this routine
c     -------------------------------
c     nh_sigma_0_vec     : nh model. yield stress at n+1
c     nh_q_u_vec         : nh model. q_u at n+1
c     nh_b_u_vec         : nh model. b_u at n+1
c     nh_h_u_vec         : nh model. h_u stress at n+1
c     nh_gamma_u_vec     : nh model. gamma_u at n+1
c
c     gp_tau_vec         : gp model tau values
c
c     sig_tol            : user defined tolerance for adaptive
c                          stress computations. Same for all
c                          elements in block.
c
c
c     NOTE:  the above data structures make it possible to have code
c            for both nonlinear_hardening and generalized_plasticity
c            written as though both have temperature dependent property
c            values. Only the gp model has temperature dependent values
c            at present.
c
c            All data values/vectors above exist and can be passed to
c            model routine for either nonlinear_hardening or
c            generalized plasticty. Some data values will be meaningless
c            (e.g. all the gp_.. vectors if the nonlinear_hardening
c            option is being used).
c
c
      call mm05( step, iter, felem, gpn, mxvl, hist_size_for_blk,
     &           nstrs, nstr, span, iout,
     &           signal_flag, adaptive_possible, cut_step_size_now,
     &           nonlin_hard, generalized_pl, sig_tol,
     &           local_work%mm05_props,
     &           local_work%e_vec,    !  at n+1
     &           local_work%e_vec_n,
     &           local_work%nu_vec,   !  at n+1
     &           local_work%nu_vec_n,
     &           local_work%gp_sig_0_vec,  ! at n+1
     &           local_work%gp_sig_0_vec_n,
     &           local_work%gp_h_u_vec,   ! at n+1
     &           local_work%gp_h_u_vec_n,
     &           local_work%gp_beta_u_vec,  ! at n+1
     &           local_work%gp_beta_u_vec_n,
     &           local_work%gp_delta_u_vec,   ! at n+1
     &           local_work%gp_delta_u_vec_n,
     &           gp_tau_vec,  ! constant over step
     &           nh_sigma_0_vec,     ! nh model. yield stress at n+1
     &           nh_sigma_0_vec,     ! nh model. yield stress at n
     &           nh_q_u_vec,         ! nh model. q_u at n+1
     &           nh_q_u_vec,         ! nh model. q_u at n
     &           nh_b_u_vec,         ! nh model. b_u at n+1
     &           nh_b_u_vec,         ! nh model. b_u at n
     &           nh_h_u_vec,         ! nh model. h_u stress at n+1
     &           nh_h_u_vec,         ! nh model. h_u stress at n
     &           nh_gamma_u_vec,     ! nh model. gama_u at n+1
     &           nh_gamma_u_vec,     ! nh model. gama_u at n
     &           local_work%rtse(1,1,gpn),
     &           local_work%urcs_blk_n(1,1,gpn),
     &           local_work%urcs_blk_n1(1,1,gpn),
     &           uddt,
     &           local_work%elem_hist(1,1,gpn),
     &           local_work%elem_hist1(1,1,gpn) )
c
      if ( adaptive_possible .and. cut_step_size_now ) then
          local_work%material_cut_step = .true.
      end if
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_06_update                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/18/02                   *
c     *                                                              *
c     *     this subroutine drives material model 06 to              *
c     *     update stresses and history for all elements in the      *
c     *     for gauss point gpn                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_06_update( gpn, props, lprops, iprops,
     &                            local_work, uddt, iout )
      use main_data, only : matprp, lmtprp
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
c
      implicit integer (a-z)
$add param_def
c
c                      parameter declarations
c
      real    props(mxelpr,*)
      logical lprops(mxelpr,*)
      integer iprops(mxelpr,*)
#dbl      double precision
#sgl      real
     &  uddt(mxvl,nstr)
$add include_sig_up
c
c
c                       locally defined variables
c
#dbl      double precision
#sgl      real
     &  gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero, gp_alpha, dtime
c
      logical signal_flag, local_debug, temperatures,
     &        temperatures_ref
      data local_debug, zero / .false., 0.0 /
c
      dtime             = local_work%dt
      span              = local_work%span
      felem             = local_work%felem
      step              = local_work%step
      iter              = local_work%iter
      now_blk           = local_work%blk
      type              = local_work%elem_type
      order             = local_work%int_order
      nnode             = local_work%num_enodes
      signal_flag       = local_work%signal_flag
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      hist_size_for_blk = local_work%hist_size_for_blk
c
c           get increment of temperature at gauss point for elements
c           in the block, the temperature at end of step and the
c           reference temperature.
c
      call gauss_pt_temps(
     &        local_work%dtemps_node_blk, gpn, type, span, order,
     &        nnode, gp_dtemps, local_work%temps_node_blk,
     &        gp_temps, temperatures, local_work%temps_node_to_process,
     &        temperatures_ref, local_work%temps_ref_node_blk,
     &        gp_rtemps )
c
c            subtract out the thermal strain increment from uddt (the
c            strain increment for step)
c
      if ( temperatures ) then
        call gp_temp_eps( span, uddt, local_work%alpha_vec, gp_dtemps,
     &                    gp_temps, gp_rtemps, local_work%alpha_vec )
      end if

c
c          for iter = 0 we just use the stored [Dt] to compute
c          stress increment. These stresses are used to dompute
c          internal forces for applied nodal displacements and
c          temperatures - material history is unaffected by
c          iter = 0
c
      if( iter .eq. 0 ) then
        igpn = gpn
        call recstr_cep_uddt_for_block( mxvl, span,
     &    gbl_cep_blocks(now_blk)%vector,
     &    uddt, local_work%urcs_blk_n(1,1,gpn),
     &    local_work%urcs_blk_n1(1,1,gpn), 1, 21, igpn )
        return
      end if

c
      call mm06( step, iter, felem, gpn, mxvl, hist_size_for_blk,
     &           nstrs, nstr, span, iout,
     &           signal_flag, adaptive_possible, cut_step_size_now,
     &           local_work%mm06_props,
     &           local_work%f0_vec, local_work%q1_vec,
     &           local_work%q2_vec, local_work%q3_vec,
     &           local_work%nuc_vec, local_work%nuc_s_n_vec,
     &           local_work%nuc_e_n_vec, local_work%nuc_f_n_vec,
     &           local_work%e_vec, local_work%tan_e_vec,
     &           local_work%nu_vec, local_work%sigyld_vec,
     &           local_work%n_power_vec, local_work%rtse(1,1,gpn),
     &           local_work%urcs_blk_n(1,1,gpn),
     &           local_work%urcs_blk_n1(1,1,gpn),
     &           uddt, local_work%elem_hist(1,1,gpn),
     &           local_work%elem_hist1(1,1,gpn),
     &           local_work%killed_status_vec )
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_07_update                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/18/02                   *
c     *                                                              *
c     *     this subroutine drives material model 07 to              *
c     *     update stresses and history for all elements in the      *
c     *     for gauss point gpn                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_07_update( gpn, props, lprops, iprops,
     &                            local_work, uddt, iout )
      use main_data, only : matprp, lmtprp
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
c
      implicit integer (a-z)
$add param_def
c
c                      parameter declarations
c
      real    props(mxelpr,*)
      logical lprops(mxelpr,*)
      integer iprops(mxelpr,*)
#dbl      double precision
#sgl      real
     &  uddt(mxvl,nstr)
$add include_sig_up
c
c
c                       locally defined variables
c
#dbl      double precision
#sgl      real
     &  gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero, gp_alpha, dtime
c
      logical signal_flag, local_debug, temperatures,
     &        temperatures_ref
      data local_debug, zero / .false., 0.0 /
c
      dtime             = local_work%dt
      span              = local_work%span
      felem             = local_work%felem
      step              = local_work%step
      iter              = local_work%iter
      now_blk           = local_work%blk
      type              = local_work%elem_type
      order             = local_work%int_order
      nnode             = local_work%num_enodes
      signal_flag       = local_work%signal_flag
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      hist_size_for_blk = local_work%hist_size_for_blk
c
c           get increment of temperature at gauss point for elements
c           in the block, the temperature at end of step and the
c           reference temperature.
c
      call gauss_pt_temps(
     &        local_work%dtemps_node_blk, gpn, type, span, order,
     &        nnode, gp_dtemps, local_work%temps_node_blk,
     &        gp_temps, temperatures, local_work%temps_node_to_process,
     &        temperatures_ref, local_work%temps_ref_node_blk,
     &        gp_rtemps )
c
c            subtract out the thermal strain increment from uddt (the
c            strain increment for step)
c
      if ( temperatures ) then
        call gp_temp_eps( span, uddt, local_work%alpha_vec, gp_dtemps,
     &                    gp_temps, gp_rtemps, local_work%alpha_vec )
      end if
c
c          for iter = 0 we just use the stored [Dt] to compute
c          stress increment. These stresses are used to dompute
c          internal forces for applied nodal displacements and
c          temperatures - material history is unaffected by
c          iter = 0

      if( iter .eq. 0 ) then
        igpn = gpn
        call recstr_cep_uddt_for_block( mxvl, span,
     &      gbl_cep_blocks(now_blk)%vector,
     &      uddt, local_work%urcs_blk_n(1,1,gpn),
     &      local_work%urcs_blk_n1(1,1,gpn), 1, 21, igpn )
        return
       end if
c
      call mm07( step, iter, felem, gpn, mxvl,  hist_size_for_blk,
     &           nstrs, nstr, span, iout,
     &           signal_flag, adaptive_possible, cut_step_size_now,
     &           local_work%mm07_props,
     &           local_work%e_vec, local_work%tan_e_vec,
     &           local_work%nu_vec, local_work%sigyld_vec,
     &           local_work%n_power_vec, local_work%rtse(1,1,gpn),
     &           local_work%urcs_blk_n(1,1,gpn),
     &           local_work%urcs_blk_n1(1,1,gpn),
     &           uddt, local_work%elem_hist(1,1,gpn),
     &           local_work%elem_hist1(1,1,gpn) )
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *            subroutine drive_umat_update  (umat)              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 03/4/2013 rhd              *
c     *                                                              *
c     *     drives material model 08 (Abaqus umat) to update         *
c     *     stresses and history for all elements in block at        *
c     *     integration point gpn                                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_umat_update( gpn, local_work, uddt, qn1, iout )
      use main_data, only : matprp, lmtprp, nonlocal_analysis
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : nonlocal_flags, nonlocal_data_n1,
     &                            gbl_cep_blocks => cep_blocks
c
      implicit integer (a-z)
$add param_def
c
c                      parameter declarations
c
#dbl      double precision
#sgl      real
     &  uddt(mxvl,nstr), qn1(mxvl,nstr,nstr)
$add include_sig_up
c
c                       locally defined variables
c
#dbl      double precision
#sgl      real
     &  gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero, one, gp_alpha, ddsddt(6), drplde(6), drpldt,
     &  big, pnewdt, predef(1), dpred(1), time(2), dtime,
     &  stress(6), stress_copy(6), stran(6), dstran(6),
     &  abq_props(50), temp_n, dtemp,
     &  statev(500), sse, spd, scd, coords(3), celent,
     &  dfgrd0(9), dfgrd1(9), dfgrd0_array(3,3), dfgrd1_array(3,3),
     &  drot(9), ddsdde(36),
     &  gp_coords(mxvl,3), symm_part_ddsdde(21), total_work_np1,
     &  plastic_work_np1, identity(9), check_key, t5, t6,
     &  temp_0, temp_n_0, s1, s2, ps(3), an(3,3),
     &  unrotated_cauchy(mxvl,6), real_npts,
     &  nonloc_ele_values(nonlocal_shared_state_size),
     &  sys_vals(nonlocal_shared_state_size)
c
      equivalence (dfgrd0, dfgrd0_array),  (dfgrd1, dfgrd1_array)
c
      logical signal_flag, local_debug, debug_now, temperatures,
     &        temperatures_ref, init_sig_eps, init_history,
     &        chk_umat_support, chk, chk2, do_nonlocal
      integer map(6)
      character * 8  cmname
      data zero, one, big, check_key / 0.0d00, 1.0d00, 1.0d06,
     &      -999999.9d00 /
      data identity /1.0d00, 0.0d00, 0.0d00, 0.0d00, 1.0d00,
     &   0.0d00, 0.0d00, 0.0d00, 1.0d00 /
      data map / 1, 2, 3, 4, 6, 5 /
c
c           1a. pull a few values from work space for block.
c               initialize selected data to make consistent
c               with usual assumptions of Abaqus umats.
c
      span              = local_work%span
      felem             = local_work%felem
      type              = local_work%elem_type
      order             = local_work%int_order
      ngp               = local_work%num_int_points
      nnode             = local_work%num_enodes
      step              = local_work%step
      iter              = local_work%iter
      now_blk           = local_work%blk
      signal_flag       = local_work%signal_flag
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      hist_size_for_blk = local_work%hist_size_for_blk
c
      knumthreads       = local_work%num_threads
      kthread           = omp_get_thread_num() + 1
c
      max_nstatv  = 500      ! dimensioned size above
      local_debug =  .false.
      chk_umat_support = .false.
c
      init_sig_eps = step .eq. 1 .and. gpn .eq. 1
      init_history = step .eq. 1 .and. gpn .eq. 1
c
      if( local_debug ) then
         write(iout,*) ' '
         write(iout,9000)
         write(iout,9001) span, felem, gpn
         write(iout,9002) step, iter
         write(iout,9004) temperatures, temperatures_ref
         write(iout,9006) hist_size_for_blk
      end if
c
c           1b. zero WARP3D data structure at start of steps so
c               values can be updated. history at n will have [D]
c               from linear stiffness if that's current computation
c               history was set in lnstff.
c
      if( init_sig_eps ) then
       call rstgp1_zero( local_work%urcs_blk_n, mxvl*nstrs*mxgp )
       call rstgp1_zero( local_work%strain_n, mxvl*nstr*mxgp )
      end if
c
      if( init_history ) then
        call rstgp1_zero( local_work%elem_hist,
     &                    span*hist_size_for_blk*ngp )
      end if
c
c           1c. values that remain constant over each element in
c               block. WARP3D uses last 21 entries
c               in history at each material point to store the
c               symmetric terms of the updated 6x6 [D]
c
      ndi    = 3
      nshr   = 3
      ntens  = 6
      npt    = gpn
      layer  = 1
      kspt   = 1
      kstep  = 1          !  WARP3D has no concept of Abaqus "step"
      kinc   = step       !  Abaqus increments are really WARP3D steps
      kiter  = iter
      kout   = iout
      nstatv = hist_size_for_blk
      if( nstatv .gt. max_nstatv ) then
         write(iout,9200) felem, max_nstatv, nstatv
         call die_abort
      end if
c
      dtime   = local_work%dt
      time(1) = local_work%time_n
      time(2) = time(1) + dtime
c
c           2. interpolate temperatures at material point for elements
c              in block
c
c           gp_rtemps: used defined reference temperature via initial
c                      conditions
c           gp_temps: total temperature at n+1
c           gp_dtemps: temperature at (n+1 - n)
c
      call gauss_pt_temps(
     &        local_work%dtemps_node_blk, gpn, type, span, order,
     &        nnode, gp_dtemps, local_work%temps_node_blk,
     &        gp_temps, temperatures, local_work%temps_node_to_process,
     &        temperatures_ref, local_work%temps_ref_node_blk,
     &        gp_rtemps )
c
c           3. subtract out the thermal strain increment from uddt (the
c              strain increment for step)
c
c              WARP3D takes care of removing thermal strains from total
c              and incremental strains when the user specifies coefficient
c              of thermal expansion (CTE) "alpha" as
c              part of the material properties outside the umat.
c
c              WARP3D has isotropic or anisotropic alpha values but
c              must be temeprature independent. If user does not specify
c              alpha values in material input, WARP3D assumes they are
c              zero. The umat must then take care of thermal strains
c              if temperature loadings are imposed.
c
c              We call the code to compute/subtract thermal strains
c              here even if alpha happens to be zero.
c
      if ( temperatures ) ! global flag set by loads processor
c                           to indicate user-specified temp changes over
c                           load step
     &    call gp_temp_eps( span, uddt, local_work%alpha_vec, gp_dtemps,
     &                    gp_temps, gp_rtemps, local_work%alpha_vec )
c
c           4. (x,y,z) coordinates at this integration point for all
c              elements in the block. umat specification requires the
c              current (n+1) coordinates for geonl. setup routines for
c              strain-stress updating of block set ce_n1 as ce @ n=0 or
c              ce @ n+1.
c
      call gauss_pt_coords( gpn, type, span, order, nnode, gp_coords,
     &                      local_work%ce_n1, iout )
c
c           5. set up identity deformation gradients F at n and n+1.
c              values for GEONL solution loaded inside element loop
c              below. Since WARP3D handles rotations for hypoelasticity
c              (like Abaqus Explicit), pass identity "drot" to umat.
c

      do k = 1, 9
       dfgrd0(k) = identity(k)
       dfgrd1(k) = identity(k)
       drot(k)   = identity(k)
      end do

c
      rpl         = zero
      ddsddt(1:6) = zero
      drplde(1:6) = zero
      drpldt      = zero

      predef(1) = zero; dpred(1) = zero

      cmname(1:) = "UMAT-WRP"

      if( local_debug ) write(iout,9008) time(1), time(2)
c
c               5. drive update over all elements in the block at this
c                  material point. the Abaqus umat processes 1
c                  point per call.
c
 1000 continue
c
      do ielem = 1, span

      noel = felem + ielem - 1
      debug_now = local_debug .and. ielem .eq. 1 .and. gpn .eq. 1
      if( debug_now )  write(iout,9100) ielem, noel, gpn
c
c               5.1 temperature at n (refernce + sum of
c                   all incremental changes)
c                   temperature increment over step
c                   reference temperature
c                   temperature at n - reference temperature
c                   pnewdt to big. if umat makes < 1.0, it wants load
c                     step size reduction
c
      dtemp    = gp_dtemps(ielem)
      temp_n   = gp_temps(ielem) - dtemp
      temp_0   = gp_rtemps(ielem)
      temp_n_0 = temp_n - temp_0
      pnewdt   = big
c
c               5.2 set properties and state(history) vectors.
c                   initialize statev at start up. load our key
c                   varaible to cehck for umat stomping over end of
c                   state variables.
c
      nprops = 50
      abq_props(1:nprops) = local_work%umat_props(ielem,1:nprops)
      statev(1:nstatv)    = local_work%elem_hist(ielem,1:nstatv,gpn)
      statev(nstatv+1)    = check_key
c
c               5.3 set stresses at n, strains at n and strain
c                   increment. adjust for themral components and
c                   swap order of xz, yz shear terms to match Abaqus.
c                   thermal part of incremental strain handled above.
c                   WARP3D processing here for strain at n requires
c                   temperature invariant CTEs (subtract reference
c                   temperature). User should not define "alpha" in
c                   WARP3D input if umat handles temperature effects.
c                   alpha below will then be zero.
c
      do j = 1, 6
       stress(map(j)) = local_work%urcs_blk_n(ielem,j,gpn)
       stress_copy(map(j)) = local_work%urcs_blk_n(ielem,j,gpn)
       stran(map(j))  = local_work%strain_n(ielem,j,gpn) -
     &                  local_work%alpha_vec(ielem,j) * temp_n_0
       dstran(map(j)) = uddt(ielem,j)
      end do
c
c               5.4 global coordinates of integration point and
c                   characteristic element length per Abaqus spec.
c                   both in deformed shape for large displacement
c                   analysis.
c
      coords(1:3) = gp_coords(ielem,1:3)
      celent = local_work%characteristic_length(ielem)
c
c               5.5 umats expect material stiffness to
c                   be initialized zero
c
      ddsdde(1:36) = zero
c
c               5.6 zero starting values of specific energy,
c                   dissipation. afterwards update WARP3D data
c                   to include increments from umat
c
      sse = zero; spd = zero; scd = zero
c
      if( debug_now ) then
       write(iout,9125) coords(1:3), celent
       write(iout,9115) dtemp, temp_n
       write(iout,9110) abq_props(1:10)
       write(iout,9105) stress(1:6), stran(1:6), dstran(1:6)
      end if
c
c               5.7 for NLGEOM, load deformation gradient F at n
c                   and n+1. otherwise pass identity.
c
      if( local_work%geo_non_flg ) then
       dfgrd0_array(1:3,1:3) = local_work%fn(ielem,1:3,1:3)
       dfgrd1_array(1:3,1:3) = local_work%fn1(ielem,1:3,1:3)
      end if
c
c               5.75 zero vector for umat to put nonlocal material
c                    state values if it wants to share with
c                    connected cohesive elements.
c
      nonloc_ele_values(1:nonlocal_shared_state_size) = zero
c
c               5.8 call the umat for material point. notice we
c                   add the current iteration number as a
c                   new last parameter. if kiter = 0, the umat
c                   is to subtract the thermal strain increment for
c                   step from dtran and immediately return (if the umat
c                   is handling temperature effects)
c
      call umat(   !   up means changeable by umat
     1   stress,   ! up
     2   statev, ! up
     3   ddsdde, ! up
     4   sse, ! up
     5   spd, ! up
     6   scd, ! up
     7   rpl, ddsddt, drplde, drpldt,
     8   stran, dstran, ! up
     9   time, dtime,
     a   temp_n, dtemp,
     b   predef, dpred,! up
     c   cmname,
     d   ndi, nshr, ntens, nstatv,
     e   abq_props, nprops,
     f   coords,
     g   drot,
     h   pnewdt, ! up
     i   celent,
     j   dfgrd0, dfgrd1,
     k   noel, npt,
     l   kslay, kspt, kstep, kinc, kiter, kout, kthread, knumthreads,
     m   nonloc_ele_values, nonlocal_shared_state_size )
c
c               5.9 make sure umat did not overwite the declared
c                   number of state variables.
c
      if( statev(nstatv+1) .ne. check_key ) then
       write(iout,9210) noel, gpn
       call die_abort
      end if
c
c              5.10 kiter = 0 processing. WARP3D uses iter = 0
c                   to compute stresses and then internal forces for
c                   imposed displacements and temperature
c                   increments. We use [Dt] from end of converged
c                   solution for previous step n * change in mechanical
c                   strain increment for imposed displ/temp increment.
c                   the umat removes thermal strain increments for
c                   imposed temp change when the umat handles thermal
c                   effects. otherwise WARP3D handled temp effects
c                   at top of this routine.
c
c                   pull [Dt] from global blocks of [Dt]s. expand to 6x6
c                   form from the 21 term vector symmetric form. multiply
c                   into dstran, add to starting stresses and cycle to
c                   next element in block. [cep] has WARP3D row order.
c                   stress returned has WARP3D row order.
c
c                   use stresses_copy since umat may have modified the
c                   stress vector passed to it (it should not if user reads
c                   our kiter - 0 info).

      if( kiter .eq. 0 ) then
c
        start_loc = ( 21 * span * (gpn-1) ) + 21 * (ielem-1)
        do k = 1, 21
          symm_part_ddsdde(k) =
     &             gbl_cep_blocks(now_blk)%vector(start_loc+k)
        end do
c
        call rstgp1_umat_cep_dstran( stress_copy, dstran,
     &                               symm_part_ddsdde )
        local_work%urcs_blk_n1(ielem,1:6,gpn) = stress_copy(1:6)
c
        if( debug_now ) then
          write(iout,9106) stress_copy(1:6)
        end if
c
        cycle ! to process next element in block
c
      end if
c
c               5.11 if wanted, check UMAT support routines
c
      if( debug_now .and. chk_umat_support ) then
         call sinv( stress, s1, s2, 3, 3 )
         call sprinc( stress, ps, 1, 3, 3 )
         call sprind( stress, ps, an, 1, 3, 3 )
         write(iout,9300) s1, s2
         write(iout,9310) ps
         call sprind( stress, ps, an, 1, 3, 3 )
         write(iout,9310) ps
         write(iout,9320) an
         write(iout,9330) pnewdt
      end if
c
c               5.11.1 check for umat wanting step reduction
c
      if ( local_work%allow_cut .and.  (pnewdt .lt. one) )
     &     local_work%material_cut_step = .true.
c
c               5.12 update WARP3D stress data structure at n+1
c                    copy updated state variables into WARP3D
c                    history @ n+1. sawp 5,6 from umat to match WARP3D
c                    make sure umat did not overwite the declared
c                    number of state variables.
c
      local_work%urcs_blk_n1(ielem,1:4,gpn)     = stress(1:4)
      local_work%urcs_blk_n1(ielem,5,gpn)       = stress(6)
      local_work%urcs_blk_n1(ielem,6,gpn)       = stress(5)
      local_work%elem_hist1(ielem,1:nstatv,gpn) = statev(1:nstatv)
c
c               5.13 update WARP3D material point work/dissipation
c                    values for increments from umat. total work
c                    in position 7, total dissipation on 8.
c
      total_work_np1   = local_work%urcs_blk_n(ielem,7,gpn) +
     &                   sse + spd + scd
      plastic_work_np1 = local_work%urcs_blk_n(ielem,8,gpn) + spd + scd
      local_work%urcs_blk_n1(ielem,7,gpn) =  total_work_np1
      local_work%urcs_blk_n1(ielem,8,gpn) =  plastic_work_np1
c
c               5.14 put updated algorithmic tangent in the global
c                    blocks of [Dt]s for the material point.
c                    we make a symmetric version
c                    and store only the 21 triangular terms. row
c                    order of stored terms is WARP3D. global cep block
c                    is 21 x span x num integration points
c
c                    Note: if UMAT returns Cauchy stress then the [D] it
c                    returns is also for Cauchy stress and deformation
c                    rate D. Will need to recognize this in stiffness
c                    update.
c
      call rstgp1_make_symmetric_store( ddsdde, symm_part_ddsdde )
      start_loc = ( 21 * span * (gpn-1) ) + 21 * (ielem-1)
      do k = 1, 21
       gbl_cep_blocks(now_blk)%vector(start_loc+k) = 
     &           symm_part_ddsdde(k)
      end do
c
c               5.15 umat had opportunity to provide an updated
c                    vector of nonlocal material values. save these
c                    to global data structure if nonlocal is
c                    active for this analysis. these are
c                    state variables to be shared with cohesive
c                    elements connected to solids.
c
c                    for first integration point of block, zero
c                    values in nonlocal data.
c
c                    for last integration point of block, average the
c                    values in nonlocal state data.
c
      do_nonlocal = .false.
      if( nonlocal_analysis ) then
        chk = nonlocal_flags(noel)
        if( chk ) then
          chk2 = allocated( nonlocal_data_n1(noel)%state_values )
          if( .not. chk2 ) then
             write(iout,9410) noel
             call die_abort
          end if
          do_nonlocal = .true.
        end if
      end if
c
      if( do_nonlocal ) then
        n = nonlocal_shared_state_size
        sys_vals(1:n) =
     &       nonlocal_data_n1(noel)%state_values(1:n)
        if( gpn .eq. 1 ) sys_vals(1:n) = zero
        sys_vals(1:n) = sys_vals(1:n) + nonloc_ele_values(1:n)
        if( gpn .eq. ngp ) then
          real_npts = ngp
          sys_vals(1:n) = sys_vals(1:n) / real_npts
        end if
         nonlocal_data_n1(noel)%state_values(1:n) =
     &            sys_vals(1:n)
      end if
c
      if( debug_now ) then
       write(iout,9120) stress(1:6)
       write(iout,9130) sse, spd, scd
       write(iout,9140) symm_part_ddsdde(1:4)
       write(iout,*) ' '
      end if
c
c               5.16 end of loop over all elements in block at this
c                    material point.
c
      end do

c                    UMAT may have computed Cauchy stresses @ n+1
c                    rather than unrotated Cauchy stresses. This
c                    would be the case for hypereleasticity type models.
c                    If so, convert t ounrotated Cauchy stresses for
c                    storage in WARP3D data structures. UMAT info
c                    routine indicates what definition the computed
c                    stressses follow.
c
      If( local_work%umat_stress_type .eq. 1 ) then  ! umat returns Cauchy
       call qmply1( span, mxvl, nstr, qn1,
     &          local_work%urcs_blk_n1(1,1,gpn), unrotated_cauchy(1,1) )
       local_work%urcs_blk_n1(1:span,1:6,gpn) =
     &               unrotated_cauchy(1:span,1:6)
      end if

c
      if( debug_now ) write(iout,9900)
c
 9000 format('>> Enter UMAT driver...')
 9001 format(5x,'span, felem, gpn: ',i4,i10,i3 )
 9002 format(5x,'step, iter: ',i8, i4 )
 9004 format(5x,'temps defined, ref temps defined: ', 2l3 )
 9006 format(5x,'num history terms: ',i4 )
 9008 format(5x,'time_n, delta_t: ', 2e20.9 )
 9100 format(5x,'... processing i, elem, gpn: ',i4,i10, i3)
 9105 format(5x,'Abauqs ordering...',
     & /,    5x,'stress_n:   ',6e14.6,/,5x,'stran_n:    ',6e14.6,
     & /,    5x,'dstran:     ',6e14.6 )
 9106 format(5x,'WARP3D ordering...',
     & /,    5x,'stress_n:   ',6e14.6 )
 9110 format(5x,'props: ',5e14.6,/,12x,5e14.6)
 9115 format(5x,'delta-temperature: ',f10.3,
     &    /, 5x,'temperature @ n:   ',f10.3 )
 9120 format(5x,'stress_n+1: ',6e14.6)
 9125 format(5x,'coords: ', 3e14.6,/,
     &       5x,'characteristic length: ',e14.6 )
 9130 format(5x,'sse, spd, scd: ',3e14.6)
 9140 format(5x,'symd 1-4: ',4e14.6)
 9200 format(5x,'... FATAL ERROR: ',
     & /10x,'umat state variable vector is too large.',
     & /,10x'current max is: ', i4,
     & /,10x'first element in block: ',i8,
     & /,10x,'aborting exectuion....' )
 9210 format(5x,'... FATAL ERROR: ',
     & /10x,'umat overwrote state variable vector',
     & /,10x'element, gpn: ',i8,i4,
     & /,10x,'aborting exectuion....' )
 9300 format("....  s1, s2: ", 2e14.6)
 9310 format("....  ps: ", 3e14.6)
 9320 format("....  an: ", 3(/,4x,3e14.6) )
 9330 format("....  pnewdt: ", e14.6 )
 9400 format(">>> FATAL ERROR: umat. @nonlocal. blk: ",i6,
     & /,".... job terminated ....",//)
 9410 format(">>> FATAL ERROR: umat. @nonlocal. noel: ",i6,
     & /,".... job terminated ....",//)
 9900 format('>> Leave UMAT driver...')
      return
      end

c     ****************************************************************
c     *                                                              *
c     *     subroutine drive_09_update  (available)                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/29/12                   *
c     *                                                              *
c     *     this subroutine drives material model 09 to              *
c     *     update stresses and history for all elements in the      *
c     *     for gauss point gpn                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_09_update( gpn, props, lprops, iprops,
     &                            local_work, uddt, iout )
      use main_data, only : matprp, lmtprp
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
c
      implicit integer (a-z)
$add param_def
c
c                      parameter declarations
c
      real    props(mxelpr,*)
      logical lprops(mxelpr,*)
      integer iprops(mxelpr,*)
#dbl      double precision
#sgl      real
     &  uddt(mxvl,nstr)
$add include_sig_up
c
c
c                       locally defined variables
c
#dbl      double precision
#sgl      real
     &  gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero, gp_alpha, dtime
c
      logical signal_flag, local_debug, temperatures,
     &        temperatures_ref
      data local_debug, zero / .false., 0.0 /
c
      dtime             = local_work%dt
      span              = local_work%span
      felem             = local_work%felem
      step              = local_work%step
      iter              = local_work%iter
      now_blk           = local_work%blk
      type              = local_work%elem_type
      order             = local_work%int_order
      nnode             = local_work%num_enodes
      signal_flag       = local_work%signal_flag
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      hist_size_for_blk = local_work%hist_size_for_blk
c
c           get increment of temperature at gauss point for elements
c           in the block, the temperature at end of step and the
c           reference temperature.
c
      call gauss_pt_temps(
     &        local_work%dtemps_node_blk, gpn, type, span, order,
     &        nnode, gp_dtemps, local_work%temps_node_blk,
     &        gp_temps, temperatures, local_work%temps_node_to_process,
     &        temperatures_ref, local_work%temps_ref_node_blk,
     &        gp_rtemps )
c
c            subtract out the thermal strain increment from uddt (the
c            strain increment for step)
c
      if ( temperatures ) then
        call gp_temp_eps( span, uddt, local_work%alpha_vec, gp_dtemps,
     &                    gp_temps, gp_rtemps, local_work%alpha_vec )
      end if
c
c          for iter = 0 we just use the stored [Dt] to compute
c          stress increment. These stresses are used to dompute
c          internal forces for applied nodal displacements and
c          temperatures - material history is unaffected by
c          iter = 0

      if( iter .eq. 0 ) then
        igpn = gpn
        call recstr_cep_uddt_for_block( mxvl, span,
     &      gbl_cep_blocks(now_blk)%vector,
     &      uddt, local_work%urcs_blk_n(1,1,gpn),
     &      local_work%urcs_blk_n1(1,1,gpn), 1, 21, igpn )
        return
       end if

c
      call mm09( step, iter, felem, gpn, mxvl, hist_size_for_blk,
     &           nstrs, nstr, span, iout,
     &           signal_flag, adaptive_possible, cut_step_size_now,
     &           local_work%umat_props,
     &           local_work%e_vec, local_work%tan_e_vec,
     &           local_work%nu_vec, local_work%sigyld_vec,
     &           local_work%n_power_vec, local_work%rtse(1,1,gpn),
     &           local_work%urcs_blk_n(1,1,gpn),
     &           local_work%urcs_blk_n1(1,1,gpn),
     &           uddt, local_work%elem_hist(1,1,gpn),
     &           local_work%elem_hist1(1,1,gpn),
     &           local_work%killed_status_vec )
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *     subroutine drive_10_update  (crystal     plasticity)     *
c     *                                                              *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 03/23/12                   *
c     *                                                              *
c     *     this subroutine drives material model 10 to              *
c     *     update stresses and history for all elements in the      *
c     *     for gauss point gpn                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_10_update( gpn, props, lprops, iprops,
     &                            local_work, uddt, iout )
      use main_data, only : matprp, lmtprp
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
c
      implicit integer (a-z)
$add param_def
c
c                      parameter declarations
c
      real    props(mxelpr,*)
      logical lprops(mxelpr,*)
      integer iprops(mxelpr,*)
#dbl      double precision
#sgl      real
     &  uddt(mxvl,nstr)
$add include_sig_up
c
c
c                       locally defined variables
c
#dbl      double precision
#sgl      real
     &  gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero, gp_alpha, dtime
c
      logical signal_flag, local_debug, temperatures,
     &        temperatures_ref
      data local_debug, zero / .false., 0.0 /
      integer :: ncrystals
c
      dtime             = local_work%dt
      span              = local_work%span
      felem             = local_work%felem
      step              = local_work%step
      iter              = local_work%iter
      type              = local_work%elem_type
      order             = local_work%int_order
      nnode             = local_work%num_enodes  
      signal_flag       = local_work%signal_flag
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      hist_size_for_blk = local_work%hist_size_for_blk
      now_blk           = local_work%blk
c
c           get increment of temperature at gauss point for elements
c           in the block, the temperature at end of step and the
c           reference temperature.
c
      call gauss_pt_temps(
     &        local_work%dtemps_node_blk, gpn, type, span, order,
     &        nnode, gp_dtemps, local_work%temps_node_blk,
     &        gp_temps, temperatures, local_work%temps_node_to_process,
     &        temperatures_ref, local_work%temps_ref_node_blk,
     &        gp_rtemps )
c
c            subtract out the thermal strain increment from uddt (the
c            strain increment for step) (We're actually going to do this
c            internally)
c
      if ( temperatures ) then
        call gp_temp_eps( span, uddt, local_work%alpha_vec, gp_dtemps,
     &                    gp_temps, gp_rtemps, local_work%alpha_vec )
      end if
c
c          for iter = 0 we just use the stored [Dt] to compute 
c          stress increment. These stresses are used to dompute 
c          internal forces for applied nodal displacements and 
c          temperatures - material history is unaffected by
c          iter = 0

      if( iter .eq. 0 ) then
        igpn = gpn
        call recstr_cep_uddt_for_block( mxvl, span,
     &      gbl_cep_blocks(now_blk)%vector,
     &      uddt, local_work%urcs_blk_n(1,1,gpn),
     &      local_work%urcs_blk_n1(1,1,gpn), 1, 21, igpn )
        return
       end if
c
      call mm10(  gpn, local_work%span, local_work%ncrystals,
     &            hist_size_for_blk,
     &            local_work%elem_hist(1,1,gpn),
     &            local_work%elem_hist1(1,1,gpn),
     &            local_work, uddt, gp_temps,
     &            gp_dtemps, iout)
c
      return
      end
c
c
c
c     ****************************************************************
c     *                                                              *
c     *     subroutine drive_11_update  (interface daamge)           *
c     *                                                              *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 03/10/14                   *
c     *                                                              *
c     *     this subroutine drives material model 11 to update       *
c     *     damage and the macroscale stress                         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_11_update( gpn, props, lprops, iprops,
     &                            local_work, uddt, iout )
      use main_data, only : matprp, lmtprp
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
c
      implicit integer (a-z)
$add param_def
c
c                      parameter declarations
c
      real    props(mxelpr,*)
      logical lprops(mxelpr,*)
      integer iprops(mxelpr,*)
#dbl      double precision
#sgl      real
     &  uddt(mxvl,nstr)
$add include_sig_up
c
c
c                       locally defined variables
c
#dbl      double precision
#sgl      real
     &  gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero, gp_alpha, dtime
c
      logical signal_flag, local_debug, temperatures,
     &        temperatures_ref
      data local_debug, zero / .false., 0.0 /
c
      dtime             = local_work%dt
      span              = local_work%span
      felem             = local_work%felem
      step              = local_work%step
      iter              = local_work%iter
      type              = local_work%elem_type
      order             = local_work%int_order
      nnode             = local_work%num_enodes  
      signal_flag       = local_work%signal_flag
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      hist_size_for_blk = local_work%hist_size_for_blk
      now_blk           = local_work%blk
c
c           get increment of temperature at gauss point for elements
c           in the block, the temperature at end of step and the
c           reference temperature.
c
      call gauss_pt_temps(
     &        local_work%dtemps_node_blk, gpn, type, span, order,
     &        nnode, gp_dtemps, local_work%temps_node_blk,
     &        gp_temps, temperatures, local_work%temps_node_to_process,
     &        temperatures_ref, local_work%temps_ref_node_blk,
     &        gp_rtemps )
c
c            subtract out the thermal strain increment from uddt (the
c            strain increment for step) (We're actually going to do this
c            internally)
c
      if ( temperatures ) then
        call gp_temp_eps( span, uddt, local_work%alpha_vec, gp_dtemps,
     &                    gp_temps, gp_rtemps, local_work%alpha_vec )
      end if
c
c          for iter = 0 we just use the stored [Dt] to compute 
c          stress increment. These stresses are used to dompute 
c          internal forces for applied nodal displacements and 
c          temperatures - material history is unaffected by
c          iter = 0

      if( iter .eq. 0 ) then
        igpn = gpn
        call recstr_cep_uddt_for_block( mxvl, span,
     &      gbl_cep_blocks(now_blk)%vector,
     &      uddt, local_work%urcs_blk_n(1,1,gpn),
     &      local_work%urcs_blk_n1(1,1,gpn), 1, 21, igpn )
        return
       end if
c     
      call mm11(  gpn, local_work%span,
     &            hist_size_for_blk,
     &            local_work%elem_hist(1,1,gpn),
     &            local_work%elem_hist1(1,1,gpn),
     &            local_work, uddt, gp_temps,
     &            gp_dtemps, iout)

      end subroutine

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine material_model_info               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 04/28/12                   *
c     *                                                              *
c     *     call the material model specific set up routine          *
c     *     to get a vector of various data sizes, parameters        *
c     *     return specific requested value to the calling routine   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine material_model_info( element_no, block_no, info_type,
     &                                 value )
      implicit integer (a-z)
$add common.main
c
c                      local data
c
      integer info_vector(10)
      integer :: inter_mat
      logical :: is_inter_dmg
c
      is_inter_dmg = .false.
c
c                      get the material model type associated with
c                      the element number or the elements in the
c                      block number. calling routine cannot set both
c                      element and block numbers.
c
c     info_vector:
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
      info_vector(1) = 0
      info_vector(2) = 21
      info_vector(3) = 0
      info_vector(4) = 0
c
      if( element_no .gt. 0 .and. block_no .gt. 0 ) then
         write(*,9000) 1
         call die_gracefully
      end if
c
      local_element_no = element_no
      if( element_no .le. 0 ) then
        if( block_no .gt. nelblk .or. block_no .le. 0 ) then
          write(*,9000) 2
          call die_gracefully
        else
          local_element_no = elblks(1,block_no)
        end if
      end if
c
      if( local_element_no .le. 0 .or.
     &    local_element_no .gt. noelem ) then
             write(*,9000) 3
             call die_gracefully
      end if
c
      mat_type = iprops(25,local_element_no)
c
c     See if we're actually a interface-damaged model
c
      if (iprops(42, local_element_no) .ne. -1) then
        inter_mat = iprops(42,local_element_no)
        is_inter_dmg = .true.
      end if
c
      select case( mat_type )
      case( 1 )
        call mm01_set_sizes( info_vector )
      case( 2 )
        call mm02_set_sizes( info_vector )
      case( 3 )
        call mm03_set_sizes( info_vector )
      case( 4 )
        call mm04_set_sizes( info_vector )
      case( 5 )
        call mm05_set_sizes( info_vector )
      case( 6 )
        call mm06_set_sizes( info_vector )
      case( 7 )
        call mm07_set_sizes( info_vector )
      case( 8 )
        call umat_set_features( info_vector )
      case( 9 )
        call mm09_set_sizes( info_vector )
      case(10 )
        call mm10_set_sizes_special( info_vector, local_element_no )
      case default
        write(*,9000) 4
        call die_gracefully
      end select
c
c     Change history length if we are actually an interface damaged material
c
      if (is_inter_dmg) then
        call mm11_set_sizes_special(inter_mat,info_vector, 
     &                              local_element_no)
      end if
c
      if( info_type .gt. 0 .and. info_type .le. 4 ) then
         value = info_vector(info_type)
      else
         write(*,9000) 5
         call die_gracefully
      end if
c
      return
c
 9000 format(/," SYSTEM Error: material_model_sizes.", /,
     &         "               error type: ", i5, /,
     &         "               Job aborted" )
      end


c     ****************************************************************
c     *                                                              *
c     *             subroutine rstgp1_update_strains                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/20/2011                 *
c     *                                                              *
c     *      computed total strains at n+1 by including the          *
c     *      increment over the current step: n+1 = n + deps         *
c     *      for geometric nonlinear theory, these are the unrotated *
c     *      strains. for linear theory, just the usual small-       *
c     *      strains. strains are in vector (not tensor) form at     *
c     *      at this point. compiler should inline this routine      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rstgp1_update_strains( span, mxvl, nstr, deps,
     &                                  strain_np1 )
      implicit integer (a-z)
c
c                      parameter declarations
c
#dbl      double precision
#sgl      real
     & deps(mxvl,*), strain_np1(mxvl,*)
c
      do i = 1, span
         strain_np1(i,1) = strain_np1(i,1) + deps(i,1)
         strain_np1(i,2) = strain_np1(i,2) + deps(i,2)
         strain_np1(i,3) = strain_np1(i,3) + deps(i,3)
         strain_np1(i,4) = strain_np1(i,4) + deps(i,4)
         strain_np1(i,5) = strain_np1(i,5) + deps(i,5)
         strain_np1(i,6) = strain_np1(i,6) + deps(i,6)
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine gauss_pt_coords                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 03/25/12                   *
c     *                                                              *
c     *     compute (x,y,z) coordinates for all elements in block    *
c     *     at this integration point                                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine gauss_pt_coords(
     &       gpn, etype, span, int_order,
     &       nnodel, gp_coords, node_coords, iout  )
      implicit integer (a-z)
$add param_def
c
c                      parameter declarations
c
#dbl      double precision
#sgl      real
     &  gp_coords(mxvl,3), node_coords(mxvl,*)
c
c                     locally defined arrays-variables
c
#dbl      double precision
#sgl      real
     &  sf(mxndel), xi, eta, zeta, weight, zero
      logical local_debug
#sgl      data zero, local_debug / 0.0, .false. /
#dbl      data zero, local_debug / 0.0d00, .false. /
c
      if( local_debug ) write(iout,*) '... in gauss_pt_coords'
c
c                     get the parametric coordinates for
c                     this integration point. then get the nodal
c                     shape functions evaluated at the point.
c
      call getgpts( etype, int_order, gpn, xi, eta, zeta, weight )
      call shapef( etype, xi, eta, zeta, sf(1) )
c
      do i = 1, span
         gp_coords(i,1) = zero
         gp_coords(i,2) = zero
         gp_coords(i,3) = zero
      end do
c
c                     interpolate (x,y,z) global coords at this gpn for
c                     each element in block. rows of node_coords have
c                     coords for element nodes -- all x-coord, then
c                     all y-coord, then all z-coord
c
      ky = nnodel
      kz = ky + nnodel
      do enode = 1, nnodel
        do i = 1, span
          gp_coords(i,1) = gp_coords(i,1)  +
     &                      sf(enode) * node_coords(i,enode)
          gp_coords(i,2) = gp_coords(i,2)  +
     &                      sf(enode) *  node_coords(i,ky+enode)
          gp_coords(i,3) = gp_coords(i,3)  +
     &                      sf(enode) *  node_coords(i,kz+enode)
        end do
      end do
c
      if ( .not. local_debug ) return
         write(iout,*) '>> in  gauss_pt_coords'
         write(iout,*) 'xi, eta, zeta:'
         write(iout,9000) xi, eta, zeta
c         write(iout,*) 'coords for element 1 of blk'
c         write(iout,9000) node_coords(1,1:3*nnodel)
      return
c
 9000 format(1x,f15.6 )
      end
c     ****************************************************************
c     *                                                              *
c     *             subroutines rstgp1_umat_cep_dstran               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 4/15/12                    *
c     *                                                              *
c     *    support for drive_umat. includes here for easy inlining   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rstgp1_umat_cep_dstran( stress, dstran, symm_vector )
      implicit none
c
c               parameters
c
#dbl      double precision
#sgl      real
     & stress(6), dstran(6), symm_vector(21)
c
c               locals
c
#dbl      double precision
#sgl      real
     & dnow(6,6), t5, t6
c
      integer i, j, k
c
c               make dstran and stress have WARP3D row ordering
c               for xz and yz shear terms
c
      t5 = dstran(5); t6 = dstran(6); dstran(5) = t6; dstran(6) = t5
      t5 = stress(5); t6 = stress(6); stress(5) = t6; stress(6) = t5
c
c               make 6x6 form of [Dt] from packed vector form of
c               symmetric storage. update stress by incremental
c               change. dstran is mechanical strain increment over
c               step.
c
      k = 1
      do i = 1, 6
        do j = 1, i
          dnow(i,j) = symm_vector(k)
          dnow(j,i) = dnow(i,j)
          k = k + 1
        end do
      end do
c
      do i = 1, 6
       do k = 1, 6
          stress(i) = stress(i) + dnow(i,k) * dstran(k)
       end do
      end do
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine rstgp1_make_symmetric_store           *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 03/26/12                   *
c     *                                                              *
c     *                   support for umat processing                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rstgp1_make_symmetric_store( matrix, symm_vector )
      implicit none
#dbl      double precision
#sgl      real
     & matrix(6,6), symm_vector(21)
c
#dbl      double precision
#sgl      real
     & transpose(6,6), symm_version(6,6), half
      integer i, j, k, map(6)
#sgl      data half / 0.5 /
#dbl      data half / 0.5d00 /
      data map / 1,2,3,4,6,5 /
c
c         1. compute transpose of 6 x 6 matrix
c         2. compute symmetrized version
c         3. swap rows, cols 5 & 6 to make shear ordering
c            compatible with WARP3D
c         4. store 21 terms in lower triangle by row
c
      do i = 1, 6
        do j = 1, 6
          transpose(i,j) = matrix(j,i)
        end do
      end do
c
      do j = 1, 6
        do i = 1, 6
          symm_version(map(i),map(j)) = half * ( matrix(i,j) +
     &                                  transpose(i,j) )
        end do
      end do
c
      k = 1
      do i = 1, 6
        do j = 1, i
          symm_vector(k) = symm_version(i,j)
          k = k + 1
        end do
      end do
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                     subroutines rstgp1s                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 4/1/12 (no folling)        *
c     *                                                              *
c     *  support routines for rstgp1. include here so they can be    *
c     *  inlined.                                                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rstgp1_a( ndof, nnode, span, ue, due, uenh, uen1,
     &                     mxvl )
      integer span
#dbl      double precision
#sgl      real
     &  ue(mxvl,*), due(mxvl,*), uenh(mxvl,*), uen1(mxvl,*),
     &  half
      data half / 0.5d00 /
c
         do j = 1, ndof*nnode
            do i = 1, span
               uenh(i,j) = ue(i,j) + half*due(i,j)
               uen1(i,j) = ue(i,j) + due(i,j)
            end do
         end do
c
      return
      end

      subroutine rstgp1_b( span, internal_energy, plastic_work,
     &                     gp_energies, gp_plast_work, det_j,
     &                     dfn1, itype )
      integer span
#dbl      double precision
#sgl      real
     &  internal_energy, plastic_work, gp_energies(*),
     &  det_j(*), dfn1(*), gp_plast_work(*)
c
      if ( itype .ne. 1 ) go to 100
      do i = 1, span
         internal_energy = internal_energy + gp_energies(i) *
     &                     dfn1(i) * det_j(i)
         plastic_work    = plastic_work +
     &                     gp_plast_work(i) * dfn1(i) * det_j(i)
      end do
      return
c
 100  continue
      do i = 1, span
         internal_energy = internal_energy + gp_energies(i) *
     &                       det_j(i)
         plastic_work    = plastic_work + gp_plast_work(i) * det_j(i)
      end do
      return
c
      end


      subroutine rstgp1_c( span, nstrs, stress_n1, urcs_blk_n1, mxvl )
      integer span
#dbl      double precision
#sgl      real
     &  stress_n1(nstrs,*), urcs_blk_n1(mxvl,*)
c
      do k = 1, nstrs
         do i = 1, span
           urcs_blk_n1(i,k) = stress_n1(k,i)
         end do
      end do
c
      return
      end

      subroutine rstgp1_d( span, nstrs, stress_n, urcs_blk_n, mxvl )
      integer span
#dbl      double precision
#sgl      real
     &  stress_n(nstrs,*), urcs_blk_n(mxvl,*)
c
      do k = 1, nstrs
         do i = 1, span
           stress_n(k,i) = urcs_blk_n(i,k)
         end do
      end do
c
      return
      end

      subroutine rstgp1_zero( vec, n )
      integer n
#dbl      double precision
#sgl      real
     &  vec(*)
#dbl      double precision
#sgl      real
     &  zero
      data zero / 0.0d00 /
c
      vec(1:n) = zero
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine recstr_cep_uddt_for_block         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/24/12                   *
c     *                                                              *
c     *     compute updated stresses for elements in this block at   *
c     *     this integration point. Uses the stored [Dt] * deps      *
c     *     already adjusted for incremental temperature             *
c     *                                                              *
c     ****************************************************************

      subroutine recstr_cep_uddt_for_block( mxvl, span, ceps_blk,
     &     deps_blk, stress_n, stress_np1, type, nrow_ceps_blk,
     &     gpn )
      implicit none
c
c                      parameter declarations
c
      integer mxvl, span, type, nrow_ceps_blk, gpn
#dbl      double precision
#sgl      real
     &  ceps_blk(nrow_ceps_blk,span,*), deps_blk(mxvl,*), 
     &  stress_n(mxvl,6), stress_np1(mxvl,6)
c
c                      locals

      integer ielem, i, j, k
#dbl      double precision
#sgl      real
     & full_cep(mxvl,6,6), zero
      data zero  / 0.0d00 /
c
c                      handle solid elements (type = 1) and
c                      cohesive elements (type = 2 ) to let
c                      compiler optimize loops.
c
      if( type .eq. 2 ) go to 1000
c
c                      expand compressed (symmetric) [Dts] to full 6x6
c                      for simplicity in coding next loop.
c
      k = 1
      do i = 1, 6
       do j = 1, i
         do ielem = 1, span
          full_cep(ielem,i,j) = ceps_blk(k,ielem,gpn)
          full_cep(ielem,j,i) = full_cep(ielem,i,j)
         end do
        k = k + 1
       end do
      end do
c
c                      compute stress @ n+1 = stress @ n + [Dt]* deps
c                      for each element in block at this
c                      integration point
c
      stress_np1(1:mxvl,1:6) = stress_n(1:mxvl,1:6)
c
      do i = 1, 6
       do k = 1, 6
         do ielem = 1, span
           stress_np1(ielem,i) = stress_np1(ielem,i) +
     &         full_cep(ielem,i,k) * deps_blk(ielem,k)
         end do
       end do
      end do
c
      return
c
c
c                      cohesive elements have 3x3 [Dt]
c
 1000 continue
      k = 1
      do i = 1, 3
       do j = 1, i
         do ielem = 1, span
          full_cep(ielem,i,j) = ceps_blk(k,ielem,gpn)
          full_cep(ielem,j,i) = full_cep(ielem,i,j)
         end do
        k = k + 1
       end do
      end do
c
c                      compute stress @ n+1 = stress @ n + [Dt]* deps
c                      for each element in block at this
c                      integration point
c
      stress_np1(1:mxvl,1:3) = stress_n(1:mxvl,1:3)
c      stress_np1(1:mxvl,1:3) = 0.0d00
c
      do i = 1, 3
       do k = 1, 3
         do ielem = 1, span
           stress_np1(ielem,i) =  stress_np1(ielem,i) +
     &         full_cep(ielem,i,k) * deps_blk(ielem,k)
         end do
       end do
      end do
c      
c      do i = 1, span
c         write(*,9000) i, stress_np1(i,1:3), deps_blk(i,1:3)
c      end do
c 9000  format('.. i, tracs: ',i4, 3f15.1,3x,3f15.6)      
c
      return
      end
