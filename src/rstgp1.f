c     ****************************************************************
c     *                                                              *
c     *                      subroutine rstgp1                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/23/2017 rhd              *
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
      use main_data, only : initial_stresses,
     &                      initial_stresses_user_routine,
     &                      initial_stresses_file
c
      implicit none
      include 'param_def'
c
c                      parameter declarations
c
      real    :: props(mxelpr,*)   ! all 3 are same by read-only
      logical :: lprops(mxelpr,*)
      integer :: iprops(mxelpr,*)
      include 'include_sig_up'
c
c                       locally defined variables
c
      integer :: i, k, span, felem, type, order, gpn, ngp, nnode, ndof,
     &           step, iter, mat_type, iout, error, ielem
      double precision :: internal_energy, beta_fact, eps_bbar,
     &                    plastic_work, bar_volumes(mxvl),
     &                    bar_areas_0(mxvl), bar_areas_mid(mxvl)
      double precision, parameter :: zero = 0.0d0
      double precision, allocatable :: ddt(:,:), uddt(:,:),
     &                                 qnhalf(:,:,:), qn1(:,:,:)
      logical :: cut_step, adaptive, geonl, bbar, material_cut_step,
     &           adaptive_flag, isnan
      logical, parameter :: local_debug = .false.
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
c        allocate and zero. only span rows are used but
c        array operators (e.g uddt = ..) will operate on
c        full content and access uninitialized values.
c
      allocate( ddt(mxvl,nstr), uddt(mxvl,nstr),
     &          qnhalf(mxvl,nstr,nstr), qn1(mxvl,nstr,nstr) )
!DIR$ VECTOR ALIGNED
      ddt    = zero
!DIR$ VECTOR ALIGNED
      uddt   = zero
!DIR$ VECTOR ALIGNED
      qnhalf = zero
!DIR$ VECTOR ALIGNED
      qn1    = zero
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
c        process bar elements separately
c
      if( local_work%is_bar_elem ) then
        if ( local_debug ) write(*,*) '>> calling gtlsn3...'
        bar_areas_0(1:span) = props(43,1:span)
        call gtlsn3_vols( span, mxvl, felem, iout, bar_areas_0(1),
     &                    bar_areas_mid(1),
     &                    local_work%ce_0, local_work%ce_mid,
     &                    bar_volumes(1) )
        call gtlsn3( span, local_work%due, uddt,
     &               local_work%ddtse(1,1,gpn), local_work%ce_mid,
     &               mxvl, nstr, local_work%elem_type, felem, iout )
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
      if( error .eq. 1 ) then
         if( adaptive ) then
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
c                creep model
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
c                 include initial stresses if present
c
      if( local_work%process_initial_stresses )  then
!DIR$ VECTOR ALIGNED
       do i = 1, span
         ielem = felem + i - 1
         local_work%urcs_blk_n1(i,1:6,gpn) =
     &     local_work%urcs_blk_n1(i,1:6,gpn) +
     &         initial_stresses(1:6,ielem)
       end do
      end if
c
c
c              If we're actually an interface damage model,
c              call the interface damage calculations
c
      if( local_work%is_inter_dmg )
     &       call drive_11_update(gpn, props, lprops, iprops,
     &            local_work, uddt, iout)
c
c --------------------------------------------------------------------
c
c          calculate the internal energy and plastic work
c          y integrating the densities over the deformed volume of the.
c          elment. urcs_blk_n1(..,7) is really the current (total) energy
c          density per unit deformed volume - look above...
c          increment of plastic work density stored in plastic_work_incr
c
      if( iter .ne. 0 ) then
        call rstgp1_b( span, internal_energy, plastic_work,
     &                 local_work%urcs_blk_n1(1,7,gpn),
     &                 local_work%urcs_blk_n1(1,8,gpn),
     &                 local_work%det_j(1,gpn), local_work%dfn1, 1,
     &                 local_work%is_bar_elem, local_work%is_link_elem,
     &                 bar_volumes(1) )
        internal_energy = internal_energy * beta_fact *
     &                    local_work%weights(gpn)
        plastic_work    = plastic_work * beta_fact *
     &                    local_work%weights(gpn)
        local_work%block_energy       = internal_energy
        local_work%block_plastic_work = plastic_work
      end if
c
      if( local_debug ) then
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
c     *                   last modified : 9/25/2017 rhd              *
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
      use main_data, only : initial_stresses,
     &                      initial_stresses_user_routine,
     &                      initial_stresses_file
      implicit none
      include 'param_def'
c
c                      parameter declarations
c
      real    :: props(mxelpr,*)  ! all 3 are same but read only
      logical :: lprops(mxelpr,*) ! 1st col is 1st elem of blk
      integer :: iprops(mxelpr,*)
      include 'include_sig_up'
c
c                       locally defined variables
c
      integer :: span, felem, type, order, gpn, ngp, nnode, ndof, step,
     &           iter, mat_type, number_points, curve_step, iout,
     &           curve_set, i, k
      double precision :: internal_energy, beta_fact, eps_bbar,
     &  uddt(mxvl,nstr), plastic_work, dummy_q(1), dummy_dfn1(1),
     &  bar_volumes(mxvl), bar_areas_0(mxvl), bar_areas_nx(mxvl)
      logical :: bbar, signal_flag, drive_material_model
      logical, save :: local_debug = .false.
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
      if( local_work%is_cohes_elem ) then
        call gtlsn2( span, nnode,
     &               local_work%due, uddt,
     &               local_work%ddtse(1,1,gpn),
     &               local_work%cohes_rot_block,
     &               local_work%shape(1,gpn),
     &               local_work%elem_type, gpn, felem, iout )
       else if( local_work%is_link_elem ) then
        call gtlsn4( span, mxelpr, props, local_work%due, uddt,
     &               local_work%ddtse(1,1,gpn),
     &               local_work%urcs_blk_n(1,1,gpn),
     &               local_work%urcs_blk_n1(1,1,gpn), mxvl, nstr,
     &               nstrs, local_work%elem_type, felem, iout )
      else if( local_work%is_bar_elem ) then
        bar_areas_0(1:span) = props(43,1:span)
        call gtlsn3_vols( span, mxvl, felem, iout, bar_areas_0(1),
     &                    bar_areas_nx(1),
     &                    local_work%ce_0, local_work%ce_0,
     &                    bar_volumes(1) )
        call gtlsn3( span, local_work%due, uddt,
     &               local_work%ddtse(1,1,gpn), local_work%ce_0,
     &               mxvl, nstr, local_work%elem_type, felem, iout )
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
      drive_material_model = .true.
      if( local_work%is_link_elem ) drive_material_model = .false.
c
c                 drive material code, include initial stresses if present
c
      if( drive_material_model ) then
         if( local_work%process_initial_stresses  )  then
!DIR$ VECTOR ALIGNED
            do i = 1, span
              local_work%urcs_blk_n(i,1:6,gpn) =
     &         local_work%initial_stresses(1:6,i)
c              if( felem+i-1 .eq. 2) then
c                write(*,*) '... loading initial stresses elem 2'
c                write(*,*) local_work%urcs_blk_n(i,1:6,gpn)
c              end if
            end do
         end if
         call rstgp2_drive_matls
c         write(*,*) '... updated stresses elem 2'
c         write(*,*) local_work%urcs_blk_n1(2,1:6,gpn)
      end if
c

c
c                 If we're actually an interface damage model,
c                 call the interface damage calculations
c
      if( local_work%is_inter_dmg )
     &       call drive_11_update(gpn, props, lprops, iprops,
     &            local_work, uddt, iout)
c
c --------------------------------------------------------------------
c
c            calculate the current (total) internal energy
c            by integrating the energy density over the volume
c            of the element. do only if iter > 0
c
      if( iter .ne. 0 ) then
        call rstgp1_b( span, internal_energy, plastic_work,
     &                 local_work%urcs_blk_n1(1,7,gpn),
     &                 local_work%urcs_blk_n1(1,8,gpn),
     &                 local_work%det_j(1,gpn), dummy_dfn1, 2,
     &                 local_work%is_bar_elem, local_work%is_link_elem,
     &                 bar_volumes(1)   )
        internal_energy = internal_energy * beta_fact *
     &                    local_work%weights(gpn)
        plastic_work    = plastic_work * beta_fact *
     &                    local_work%weights(gpn)
        local_work%block_energy       = internal_energy
        local_work%block_plastic_work = plastic_work
      end if
c
      if( local_debug ) then
        write(iout,*) '>> rstgp2 .. gauss point: ', gpn
        write (iout,9500) internal_energy,
     &                 plastic_work
        write(iout,9110)
        do i = 1, span
         write(iout,9100) i, (local_work%urcs_blk_n(i,k,gpn),k=1,7),
     &                 (local_work%urcs_blk_n1(i,k,gpn),k=1,7),
     &                 (uddt(i,k),k=1,6)
        end do
      end if

c
      return
c
 9100 format(i5,7f15.6,/,5x,7f15.6,/,5x,6f15.6,/,5x,6f15.6)
 9110 format(1x,'Elem    /',20('-'),
     &        ' unrot. Cauchy @ n, unrot. Cauchy @ n+1,',
     &        ' uddt', 20('-'),'/')
 9500 format('  Internal energy inside of (rstgp2)    = ',e16.6,
     &     /,'  Plasstic work inside of (rstgp2)      = ',e16.6)
 9510 format(2x,i3,6f10.2)
c
      contains
c     ========
c
      subroutine rstgp2_drive_matls
      implicit none
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
c                creep model
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
c                available
c
       call drive_09_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
c
      case ( 10 )
c
c                crystal plasticity.  make not make much sense to use
c                small displacement formulation
c
       call drive_10_update( gpn, props, lprops, iprops,
     &                       local_work, uddt, iout )
c
      case default
        write(iout,*) '>>> invalid material model number'
        call die_abort
        stop
      end select
c
      return
      end subroutine rstgp2_drive_matls
      end subroutine rstgp2
c     ****************************************************************
c     *                                                              *
c     *                subroutine drive_01_update                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *              last modified : 12/21/2015 rhd                  *
c     *                                                              *
c     *     drives material model #1 (bilinear) to                   *
c     *     update stresses and history for all elements in the      *
c     *     block for gauss point gpn                                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_01_update( gpn, props, lprops, iprops,
     &                            local_work, uddt_displ, iout )
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
      use main_data, only : extrapolated_du, non_zero_imposed_du
c
      implicit none
      include 'param_def'
c
c                      parameter declarations
c
      real ::    props(mxelpr,*)   ! all same but read only
      logical :: lprops(mxelpr,*)
      integer :: iprops(mxelpr,*)
      integer :: gpn, iout
      double precision ::  uddt_displ(mxvl,nstr)
      include 'include_sig_up'
c
c                       locally defined variables
c
      integer :: span, felem, type, order, ngp, nnode, ndof, step,
     &           iter, now_blk, mat_type, number_points, curve_set,
     &           hist_size_for_blk, curve_type, elem_type, i
c
      double precision ::
     &  dtime, gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero,  ddummy(1), gp_alpha, ymfgm, et, uddt_temps(mxvl,nstr),
     &  uddt(mxvl,nstr), cep(mxvl,6,6)
c
      logical :: geonl, local_debug, temperatures, segmental,
     &           temperatures_ref, fgm_enode_props
c
      data zero / 0.0d0 /
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
      mat_type          = local_work%mat_type
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      segmental         = local_work%segmental
      number_points     = local_work%number_points
      curve_set         = local_work%curve_set_number
      fgm_enode_props   = local_work%fgm_enode_props
      hist_size_for_blk = local_work%hist_size_for_blk
c
      local_debug       = .false. ! felem .eq. 1 .and. gpn .eq. 3
      if( local_debug ) then
        write(iout,9000) felem, gpn, span
        write(iout,9010) dtime, type, order, nnode, ndof, geonl, step,
     &                   iter, now_blk, mat_type,
     &                   temperatures, temperatures_ref, segmental,
     &                   number_points, curve_set,
     &                   fgm_enode_props, hist_size_for_blk
      end if
c
c          determine if the material elastic and properties are
c          described by segmental curve(s) and if they are temperature
c          or strain rate dependent [=0 no dependence, =1 temperature
c          dependent, =2 strain-rate dependent (not used here)]
c
      curve_type = -1
      if( segmental ) call set_segmental_type( curve_set, curve_type,
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
      if( fgm_enode_props ) then
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
          if( local_debug ) write(iout,9020)
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
      if( local_debug ) write(iout,9030)
c
c          get temperature dependent young's modulus, poisson's
c          ratio, uniaxial yield stress and constant plastic
c          modulus for the temperature at n and n+1. Some
c          properties are not for this model but need to be passed
c          to satisfy syntax of call.
c
      if( curve_type .eq. 0 .or. curve_type .eq. 1 ) then
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
         if( local_debug ) write(iout,9040)
      end if
c
c          get the thermal strain increment (actually negative of increment)
c
!DIR$ VECTOR ALIGNED
      uddt_temps = zero
      if ( temperatures ) then
        call gp_temp_eps( span, uddt_temps,
     &                    local_work%alpha_vec, gp_dtemps ,
     &                    gp_temps, gp_rtemps,
     &                    local_work%alpha_vec_n, type )
        if( local_debug ) write(iout,9050)
      end if
c
c            uddt_displ - strain increment due to displacement increment
c            uddt_temps - (negative) of strain increment just due
c                         to temperature change
c
!DIR$ VECTOR ALIGNED
      uddt = uddt_displ + uddt_temps
!DIR$ VECTOR ALIGNED
      cep  = zero
c
      do i = 1, span
       if( local_work%killed_status_vec(i) ) uddt(i,1:nstr) = zero
      end do
c
c            now standard update process. use nonlinear update and [D]
c            computation for iter = 0 and extrapolation or iter > 1
c            for iter = 0 and no extrapolation, use linear-elastic [D]
c            with props at n+1.

      if( iter >= 1 .or. extrapolated_du ) then !nonlinear update
       if( local_debug ) write(iout,9060)
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
       call cnst1( span, cep, local_work%rtse(1,1,gpn),
     &            local_work%nu_vec,
     &            local_work%e_vec, local_work%elem_hist1(1,2,gpn),
     &            local_work%elem_hist1(1,5,gpn), local_work%beta_vec,
     &            local_work%elem_hist1(1,1,gpn),
     &            local_work%elem_hist1(1,4,gpn), felem, iout )
      else  ! linear-elastic update
        if( local_debug ) write(iout,9070)
        call drive_01_update_a
      end if
c
c          save the [D] matrices (lower-triangle)
c
      call rstgp1_store_cep( span, mxvl, gpn,
     &         gbl_cep_blocks(now_blk)%vector, 21, cep )
      if( local_debug ) write(iout,9080)
c
      return
c
 9000 format(1x,'.... debug mm01. felem, gpn, span: ',i7,i3,i3)
 9010 format(10x,'...dtime, type, order, nnode, ndof:',e14.6,4i5,
     &     /,10x,'...geonl, step, iter, now_blk, mat_type: ',l2,4i5,
     &     /,10x,'...temperatures, temperatures_ref: ',
     &               2l2,
     &     /,10x,'...segmental, number_points, curve_set: ',l2,i3,i3,
     &     /,10x,'...fgm_enode_props, hist_size_for_blk: ',
     &    l3,i4 )
 9610 format(' >> rate iterations to converge: ',i3 )
 9020 format(10x,'...fgm properties determined...')
 9030 format(10x,'...temperatures computed at integration point...')
 9040 format(10x,'...temperatures dependent properties computed...')
 9050 format(10x,'...thermal strains computed...')
 9060 format(10x,'...update stresses nonlinear procedure...' )
 9070 format(10x,'...update stresses use linear [D]...' )
 9080 format(10x,'...[D]s saved to global structure...')
c
      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_01_update_a                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/21/2015 rhd             *
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_01_update_a
      implicit none
c
      integer :: k, m, i
      double precision :: one, two, e, nu, c1, c2, c3, c4
      data one, two / 1.0d00, 2.0d00 /
c
c              get linear-elastic [D] with potentially temperature
c              dependent properties
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
         if( local_work%killed_status_vec(i) ) cycle
         e  = local_work%e_vec(i)
         nu = local_work%nu_vec(i)
         c1 = (e/((one+nu)*(one-two*nu)))
         c2 = (one-nu)*c1
         c3 = ((one-two*nu)/two)*c1
         c4 = nu*c1
         cep(i,1,1)= c2
         cep(i,2,2)= c2
         cep(i,3,3)= c2
         cep(i,4,4)= c3
         cep(i,5,5)= c3
         cep(i,6,6)= c3
         cep(i,1,2)= c4
         cep(i,1,3)= c4
         cep(i,2,1)= c4
         cep(i,3,1)= c4
         cep(i,2,3)= c4
         cep(i,3,2)= c4
      end do
c
c              stresses at n+1 using linear-elastic [D]
c
       call drive_01_update_b( span, mxvl, uddt, cep,
     &                        local_work%urcs_blk_n(1,1,gpn),
     &                        local_work%urcs_blk_n1(1,1,gpn),
     &                        local_work%killed_status_vec )
c
      return
      end subroutine drive_01_update_a
c
      end subroutine drive_01_update
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_01_update_b                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/21/2015 rhd             *
c     *                                                              *
c     *     support routine for mm01 material driver.                *
c     *     should be inlined                                        *
c     *                                                              *
c     ****************************************************************

      subroutine drive_01_update_b( span, mxvl, uddt,
     &                              local_cep, stress_n, stress_np1,
     &                              killed_status )
      implicit none
c
      integer :: span, mxvl
      logical :: killed_status(*)
      double precision ::
     &  local_cep(mxvl,6,6), stress_n(mxvl,6), stress_np1(mxvl,6),
     &  uddt(mxvl,6), zero
      data zero / 0.0d00 /
c
      integer i, k, m
c
c              for each element in block, update stresses by
c              [D-elastic] * uddt. uddt contains thermal increment +
c              increment from imposed nodal displacements
c
!DIR$ VECTOR ALIGNED
      stress_np1 = stress_n
c
      do k = 1, 6
       do m = 1, 6
!DIR$ VECTOR ALIGNED
         do i = 1, span
           stress_np1(i,k) = stress_np1(i,k) +
     &                       local_cep(i,m,k) * uddt(i,m)
         end do
       end do
      end do
c
      do i = 1, span
        if( killed_status(i) ) stress_np1(i,1:6) = zero
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_02_update                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/27/2015 rhd             *
c     *                                                              *
c     *     this subroutine drives material model 02 to              *
c     *     update stresses and history for all elements in the      *
c     *     block for an integration point                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_02_update( gpn, props, lprops, iprops,
     &                            local_work, uddt_displ, iout )
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
      use main_data, only : extrapolated_du, non_zero_imposed_du
c
      implicit none
      include 'param_def'
c
c                      parameter declarations
c
      real ::    props(mxelpr,*)   ! all same but read only
      logical :: lprops(mxelpr,*)
      integer :: iprops(mxelpr,*)
      integer :: gpn, iout
      double precision ::  uddt_displ(mxvl,nstr)
      include 'include_sig_up'

c
c                       locally defined variables
c
      integer :: span, felem, type, order, ngp, nnode, ndof, step,
     &           iter, now_blk, elem_type, mat_type,
     &           hist_size_for_blk, i
c
      double precision ::
     &  dtime, gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero,  gp_alpha, cep(mxvl,6,6), ddtse(mxvl,6), nowtemp
c
      logical :: geonl, local_debug, temperatures,
     &           temperatures_ref, fgm_enode_props, signal_flag
c
      data zero / 0.0d0 /
c
c          deformation plasticity model. properties are invariant of
c          temperature and loading rate. properties may vary spatially
c          for fgms - values here are interploated at the current
c          integration point.
c
      span              = local_work%span
      felem             = local_work%felem
      step              = local_work%step
      iter              = local_work%iter
      now_blk           = local_work%blk
      type              = local_work%elem_type
      order             = local_work%int_order
      ndof              = local_work%num_enode_dof
      nnode             = local_work%num_enodes
      mat_type          = local_work%mat_type
      signal_flag       = local_work%signal_flag
      fgm_enode_props   = local_work%fgm_enode_props
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      hist_size_for_blk = local_work%hist_size_for_blk
c
      local_debug       = .false. ! felem .eq. 1 .and. gpn .eq. 3
      if( local_debug ) then
        write(iout,9000) felem, gpn, span
        write(iout,9010) dtime, type, order, nnode, ndof, geonl, step,
     &                   iter, now_blk, mat_type,
     &                   temperatures, temperatures_ref,
     &                   fgm_enode_props, hist_size_for_blk
      end if
c
c          for fgms, interpolate values of material properties
c          at the current gauss point (overwrite the constant values).
c          at present only e, nu, sig_yld, n_power and
c          isotropic alpha can be specified
c          as fgm properties interpolated from nodal values.
c
      if( fgm_enode_props ) then
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
!DIR$ VECTOR ALIGNED
          do i = 1, span
            gp_alpha  = local_work%alpha_vec(i,1)
            local_work%alpha_vec(i,1)   = gp_alpha
            local_work%alpha_vec(i,2)   = gp_alpha
            local_work%alpha_vec(i,3)   = gp_alpha
            local_work%alpha_vec(i,4)   = zero
            local_work%alpha_vec(i,5)   = zero
            local_work%alpha_vec(i,6)   = zero
          end do
          if( local_debug ) write(iout,9020)
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
      if( local_debug ) write(iout,9030)
c
c          process iter > 0 or iter=0 and extrapolated du. mm02
c          uses total strains adjusted for total
c          temperature strains to compute stresses @ n+1.
c
      if( temperatures ) then
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
         if( local_debug ) write(iout,9090)
      else
         ddtse = local_work%ddtse(1:span,1:6,gpn)
      end if
c
      if( iter >= 1 .or. extrapolated_du ) then !nonlinear update
         if( local_debug ) write(iout,9060)
         call mm02( step, iter, felem, gpn, local_work%e_vec,
     &           local_work%nu_vec, local_work%sigyld_vec,
     &           local_work%n_power_vec,local_work%urcs_blk_n(1,1,gpn),
     &           local_work%urcs_blk_n1(1,1,gpn),
     &           ddtse(1,1),  ! total strain @ n+1
     &           local_work%elem_hist(1,1,gpn),
     &           local_work%elem_hist1(1,1,gpn), span, iout,
     &           signal_flag )
         call cnst2( felem, gpn, local_work%e_vec, local_work%nu_vec,
     &               local_work%sigyld_vec, local_work%n_power_vec,
     &               ddtse, local_work%elem_hist1(1,1,gpn),
     &               cep, span, iout )
      else
        if( local_debug ) write(iout,9070)
!DIR$ VECTOR ALIGNED
        cep  = zero
        call drive_02_update_a
       if( felem .eq. 1 ) then
           write(*,*) '... stresses _n+1 element 1 from 02_update_a '
           write(*,*) local_work%urcs_blk_n1(1,1:6,gpn)
         end if
      end if
c
c          save the [D] matrices (lower-triangle)
c
      call rstgp1_store_cep( span, mxvl, gpn,
     &         gbl_cep_blocks(now_blk)%vector, 21, cep )
      if( local_debug ) write(iout,9080)
c
      return
c
 9000 format(1x,'.... debug mm02. felem, gpn, span: ',i7,i3,i3)
 9010 format(10x,'...dtime, type, order, nnode, ndof:',e14.6,4i5,
     &     /,10x,'...geonl, step, iter, now_blk, mat_type: ',l2,4i5,
     &     /,10x,'...temperatures, temperatures_ref: ',
     &               2l2,
     &     /,10x,'...fgm_enode_props, hist_size_for_blk: ',
     &    l3,i4 )
 9610 format(' >> rate iterations to converge: ',i3 )
 9020 format(10x,'... fgm properties determined...')
 9030 format(10x,'... temperatures computed at integration point ...')
 9040 format(10x,'... temperatures dependent properties computed ...')
 9050 format(10x,'... thermal strains computed ...')
 9060 format(10x,'... update stresses nonlinear procedure ...' )
 9070 format(10x,'... update stresses use linear [D] ...' )
 9080 format(10x,'...[ D]s saved to global structure ...')
 9090 format(10x,'... temperature effects subtracted ...')
c
      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_02_update_a                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/26/2015 rhd             *
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_02_update_a
      implicit none
c
      integer :: k, m, i
      double precision :: one, two, e, nu, c1, c2, c3, c4
      data one, two / 1.0d00, 2.0d00 /
c
c              get linear-elastic [D] with potentially temperature
c              dependent properties
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
         if( local_work%killed_status_vec(i) ) cycle
         e  = local_work%e_vec(i)
         nu = local_work%nu_vec(i)
         c1 = (e/((one+nu)*(one-two*nu)))
         c2 = (one-nu)*c1
         c3 = ((one-two*nu)/two)*c1
         c4 = nu*c1
         cep(i,1,1)= c2
         cep(i,2,2)= c2
         cep(i,3,3)= c2
         cep(i,4,4)= c3
         cep(i,5,5)= c3
         cep(i,6,6)= c3
         cep(i,1,2)= c4
         cep(i,1,3)= c4
         cep(i,2,1)= c4
         cep(i,3,1)= c4
         cep(i,2,3)= c4
         cep(i,3,2)= c4
      end do
c
c              stresses at n+1 using linear-elastic [D]
c
      call drive_02_update_b( span, mxvl, ddtse, cep,
     &                        local_work%initial_stresses,
     &                        local_work%urcs_blk_n1(1,1,gpn),
     &                        local_work%killed_status_vec )

c
      return
      end subroutine drive_02_update_a
c
      end subroutine drive_02_update
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_02_update_b                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/26/2015 rhd             *
c     *                                                              *
c     *     support routine for mm02 material driver.                *
c     *     should be inlined                                        *
c     *                                                              *
c     ****************************************************************

      subroutine drive_02_update_b( span, mxvl, ddtse, local_cep,
     &                              initial_stresses, stress_np1,
     &                              killed_status )
      implicit none
c
      integer :: span, mxvl
      logical :: killed_status(*)
      double precision :: initial_stresses(6,span),
     &                    local_cep(mxvl,6,6), stress_np1(mxvl,6),
     &                    ddtse(mxvl,6), zero
      data zero / 0.0d00 /
c
      integer i, k, m
c
c              for each element in block, update stresses by
c              [D-elastic] * uddt. uddt contains thermal increment +
c              increment from imposed nodal displacements
c
!DIR$ VECTOR ALIGNED
      stress_np1 = zero
c
      do k = 1, 6
       do m = 1, 6
!DIR$ VECTOR ALIGNED
         do i = 1, span
           stress_np1(i,k) = stress_np1(i,k) +
     &                       local_cep(i,m,k) * ddtse(i,m)
         end do
       end do
      end do
c
      do k = 1, 6
!DIR$ VECTOR ALIGNED
       do i = 1, span
         stress_np1(i,k) = stress_np1(i,k) + initial_stresses(k,i)
       end do
      end do
c
      do i = 1, span
        if( killed_status(i) ) stress_np1(i,1:6) = zero
      end do
c
      return
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
      subroutine drive_02_update_old( gpn, props, lprops, iprops,
     &                            local_work, uddt, iout )
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
c
      implicit integer (a-z)
      include 'param_def'
c
c                      parameter declarations
c
      real    props(mxelpr,*)
      logical lprops(mxelpr,*)
      integer iprops(mxelpr,*)
      double precision
     &  uddt(mxvl,nstr)
      include 'include_sig_up'
c
c
c                       locally defined variables
c
      double precision
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
     &                    gp_temps, gp_rtemps, local_work%alpha_vec,
     &                    type )
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
!DIR$ VECTOR ALIGNED
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
c     *             last modified : 11/9/2015 rhd                    *
c     *                                                              *
c     *     drives material model 03 to update stresses and history  *
c     *     for all elements in the block at 1 integration point     *
c     *     (gpn)                                                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_03_update( gpn, props, lprops, iprops,
     &                            local_work, uddt_displ, iout )
      use segmental_curves, only : max_seg_points, max_seg_curves,
     &                             now_blk_relem, sigma_curve_min_values
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
      use main_data, only : extrapolated_du, non_zero_imposed_du
c
      implicit none
      include 'param_def'
c
c                      parameter declarations
c
      real ::    props(mxelpr,*)   ! all same but read only
      logical :: lprops(mxelpr,*)
      integer :: iprops(mxelpr,*)
      integer :: gpn, iout
      double precision ::  uddt_displ(mxvl,nstr)
      include 'include_sig_up'
c
c                       locally defined variables
c
      double precision ::
     &  dtime, internal_energy, beta_fact, eps_bbar, ddt(mxvl,nstr),
     &  q(mxvl,nstr,nstr), stress_n(nstrs,mxvl),
     &  stress_n1(nstrs,mxvl), p_trial(mxvl), q_trial(mxvl), dfn1(mxvl),
     &  yld_func(mxvl), step_scale_fact, plastic_work, gp_temps(mxvl),
     &  gp_dtemps(mxvl), plastic_eps_rates(mxvl), gp_rtemps(mxvl),
     &  zero, gp_alpha, et, ymfgm, copy_sigyld_vec(mxvl),
     &  trans_factor, curve_min_value, uddt_temps(mxvl,nstr),
     &  uddt(mxvl,nstr), cep(mxvl,6,6)
c
      logical :: null_point(mxvl), cut_step, process_block,
     &           adaptive, geonl, bbar, material_cut_step,
     &           local_debug, signal_flag, adaptive_flag,
     &           power_law, temperatures, allow_cut, segmental,
     &           model_update, temperatures_ref, fgm_enode_props
c
      integer :: span, felem, type, order, ngp, nnode, ndof, step,
     &           iter, now_blk, mat_type, number_points, curve_set,
     &           hist_size_for_blk, curve_type, elem_type, i,
     &           numrows_stress
c
      data zero, trans_factor / 0.0d00, 0.95d00 /
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
      local_debug       = .false. ! felem .eq. 1 .and. gpn .eq. 1
      if( local_debug ) then
        write(iout,9000) felem, gpn, span
        write(iout,9010) dtime, type, order, nnode, ndof, geonl, step,
     &                   iter, now_blk, mat_type, adaptive_flag,
     &                   temperatures, temperatures_ref, segmental,
     &                   number_points, curve_set, power_law,
     &                   fgm_enode_props, hist_size_for_blk
      end if
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
      if( fgm_enode_props ) then
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
          if( local_debug ) write(iout,9020)
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
      if( local_debug ) write(iout,9030)
c
c          get temperature dependent elastic and flow properites.
c          build local copy of stress vs. plastic strain curve
c          for use at this gauss point for elements in block.
c          set yield stress at zero plastic strain.  Some
c          properties are not for this model but need to be passed
c          to satisfy syntax of call.
c
      if( curve_type .eq. 0 .or. curve_type .eq. 1 ) then
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
         if( local_debug ) write(iout,9040)
      end if
c
c          compute negative of incremental thermal strains
c
!DIR$ VECTOR ALIGNED
      uddt_temps = zero
      if( temperatures ) then
        call gp_temp_eps( span, uddt_temps, local_work%alpha_vec,
     &                     gp_dtemps, gp_temps, gp_rtemps,
     &                     local_work%alpha_vec_n, type )
        if( local_debug ) write(iout,9050)
      end if
c
c
c            uddt_displ - strain increment due to displacement increment
c            uddt_temps - (negative) of strain increment just due
c                         to temperature change
c
!DIR$ VECTOR ALIGNED
      uddt = uddt_displ + uddt_temps
!DIR$ VECTOR ALIGNED
      cep  = zero
      do i = 1, span
       if( local_work%killed_status_vec(i) ) uddt(i,1:nstr) = zero
      end do
c
c                set up local stress array at state 'n'. used by
c                nonlinear update.
c
      numrows_stress = nstrs
      call drive_03_update_d( span, numrows_stress, stress_n,
     &               local_work%urcs_blk_n(1,1,gpn), mxvl )
c
      if( iter >= 1 .or. extrapolated_du ) then !nonlinear update
        if( local_debug ) write(iout,9060)
        call drive_03_update_a
        call drive_03_update_c( span, numrows_stress, stress_n1,
     &               local_work%urcs_blk_n1(1,1,gpn), mxvl )
      else  ! linear-elastic update
        if( local_debug ) write(iout,9070)
        call drive_03_update_b
      end if
c
c          save the [D] matrices (lower-triangle)
c
      call rstgp1_store_cep( span, mxvl, gpn,
     &         gbl_cep_blocks(now_blk)%vector, 21, cep )
      if( local_debug ) write(iout,9080)
c
 9000 format(1x,'.... debug mm03. felem, gpn, span: ',i7,i3,i3)
 9010 format(10x,'...dtime, type, order, nnode, ndof:',e14.6,4i5,
     &     /,10x,'...geonl, step, iter, now_blk, mat_type: ',l2,4i5,
     &     /,10x,'...adaptive_flag, temperatures, temperatures_ref: ',
     &               3l2,
     &     /,10x,'...segmental, number_points, curve_set: ',l2,i3,i3,
     &     /,10x,'...power_law, fgm_enode_props, hist_size_for_blk: ',
     &    2l3,i4 )
 9020 format(10x,'...fgm properties determined...')
 9030 format(10x,'...temperatures computed at integration point...')
 9040 format(10x,'...temperatures dependent properties computed...')
 9050 format(10x,'...thermal strains computed...')
 9060 format(10x,'...update stresses nonlinear procedure...' )
 9070 format(10x,'...update stresses use linear [D]...' )
 9080 format(10x,'...[D]s saved to global structure...')
c
      return
c
      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_03_update_b                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/9/2015 rhd
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_03_update_b
      implicit none
c
      integer :: k, m, i
      double precision :: one, two, e, nu, c1, c2, c3, c4
      data one, two / 1.0d00, 2.0d00 /
c
c              get linear-elastic [D] with potentially temperature
c              dependent properties
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
         if( local_work%killed_status_vec(i) ) cycle
         e  = local_work%e_vec(i)
         nu = local_work%nu_vec(i)
         c1 = (e/((one+nu)*(one-two*nu)))
         c2 = (one-nu)*c1
         c3 = ((one-two*nu)/two)*c1
         c4 = nu*c1
         cep(i,1,1)= c2
         cep(i,2,2)= c2
         cep(i,3,3)= c2
         cep(i,4,4)= c3
         cep(i,5,5)= c3
         cep(i,6,6)= c3
         cep(i,1,2)= c4
         cep(i,1,3)= c4
         cep(i,2,1)= c4
         cep(i,3,1)= c4
         cep(i,2,3)= c4
         cep(i,3,2)= c4
      end do
c
c              stresses at n+1 using linear-elastic [D]
c
      call drive_03_update_e( span, mxvl, uddt, cep,
     &                        local_work%urcs_blk_n(1,1,gpn),
     &                        local_work%urcs_blk_n1(1,1,gpn),
     &                        local_work%killed_status_vec )
c
      return
      end subroutine drive_03_update_b
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_03_update_a                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/9/2015 rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_03_update_a
      implicit none
c
      type :: arguments
        integer :: iter, abs_element, relem, ipoint, iout
        logical :: allow_cut, segmental, power_law,
     &             rate_depend_segmental, signal_flag, cut_step
        double precision :: dtime, step_scale_fact
      end type
c
      type (arguments) ::args
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
!DIR$ VECTOR ALIGNED
      copy_sigyld_vec(1:span) = local_work%sigyld_vec(1:span)
      if( .not. segmental ) then
        do i = 1, span
         if( local_work%n_power_vec(i) .gt. zero )
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
      if( process_block ) then
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
        if( null_point(i) ) cycle
        if( .not. local_work%nonlinear_flag(i) ) cycle
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
      end do ! on span
c
      local_work%material_cut_step = material_cut_step
      if( material_cut_step ) return

      end if   ! on process_block
c
      call cnst3(
     &  felem, gpn, iter, local_work%e_vec,
     &  local_work%nu_vec, local_work%q1_vec, local_work%q2_vec,
     &  local_work%q3_vec, local_work%nuc_vec,
     &  local_work%nuc_s_n_vec, local_work%nuc_e_n_vec,
     &  local_work%nuc_f_n_vec,
     &  local_work%rtse(1,1,gpn), local_work%elem_hist(1,1,gpn),
     &  local_work%elem_hist1(1,1,gpn), cep, span, iout )
c
      do i = 1, span
        if( local_work%killed_status_vec(i) ) then
            cep(i,1:6,1:6) = zero
            stress_n1(1:nstrs,i) = zero
        end if
      end do
c
      return
c
      end subroutine drive_03_update_a
      end subroutine drive_03_update
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_03_update_e                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/9/2015 rhd              *
c     *                                                              *
c     *     support routine for mm03 material driver.                *
c     *     should be inlined                                        *
c     *                                                              *
c     ****************************************************************

      subroutine drive_03_update_e( span, mxvl, uddt,
     &                              local_cep, stress_n, stress_np1,
     &                              killed_status )
      implicit none
c
      integer :: span, mxvl
      logical :: killed_status(*)
      double precision ::
     &  local_cep(mxvl,6,6), stress_n(mxvl,6), stress_np1(mxvl,6),
     &  uddt(mxvl,6), zero
      data zero / 0.0d00 /
c
      integer i, k, m
c
c              for each element in block, update stresses by
c              [D-elastic] * uddt. uddt contains thermal increment +
c              increment from imposed nodal displacements
c
!DIR$ VECTOR ALIGNED
      stress_np1 = stress_n
c
      do k = 1, 6
       do m = 1, 6
!DIR$ VECTOR ALIGNED
         do i = 1, span
           stress_np1(i,k) = stress_np1(i,k) +
     &                       local_cep(i,m,k) * uddt(i,m)
         end do
       end do
      end do
c
      do i = 1, span
        if( killed_status(i) ) stress_np1(i,1:6) = zero
      end do
c
      return
      end


c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_03_update_c                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/9/2015 rhd
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_03_update_c( span, nstrs, stress_n1,
     &                              urcs_blk_n1, mxvl )
      implicit none
c
      integer :: span, nstrs, mxvl
      double precision ::
     &  stress_n1(nstrs,*), urcs_blk_n1(mxvl,*)
c
      integer :: k, i
c
      do k = 1, nstrs  !  not necessarily = 6
!DIR$ VECTOR ALIGNED
         do i = 1, span
           urcs_blk_n1(i,k) = stress_n1(k,i)
         end do
      end do
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_03_update_d                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/9/2015 rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_03_update_d( span, nstrs, stress_n,
     &                              urcs_blk_n, mxvl )
      implicit none
c
      integer ::  span, nstrs, mxvl
      double precision ::
     &  stress_n(nstrs,*), urcs_blk_n(mxvl,*)
      integer :: k, i
c
      do k = 1, nstrs    !  not necessarily = 6
!DIR$ VECTOR ALIGNED
         do i = 1, span
           stress_n(k,i) = urcs_blk_n(i,k)
         end do
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_04_update                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 10/26/2015 rhd             *
c     *                                                              *
c     *     drive model 04 interface-cohesive to                     *
c     *     update stresses, history, [D]s for all elements in the   *
c     *     block for gauss point gpn                                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_04_update( gpn, props, lprops, iprops,
     &                            local_work, uddt, iout )
      use main_data, only : matprp, lmtprp, extrapolated_du,
     &                      non_zero_imposed_du
      use segmental_curves, only : max_seg_points
      use elem_block_data, only  : gbl_cep_blocks => cep_blocks
c
      implicit none
      include 'param_def'
c
c              parameter declarations
c
      integer ::  gpn, iout
      real    ::  props(mxelpr,*)   ! all same but readonly
      logical ::  lprops(mxelpr,*)
      integer ::  iprops(mxelpr,*)
      double precision :: uddt(mxvl,nstr)
      include 'include_sig_up'
c
c              locally defined variables
c
      double precision ::
     &  time_n, dtime,  cep(mxvl,6,6), ddummy(mxvl), ds1,
     &  ds2, dn, tns1, tns2, tnn, zero
c
      integer :: i, idummy(mxvl), span, felem, step, now_blk, iter,
     &           nnode, knumthreads, kthread, imxvl

      logical :: fgm_enode_props, nonlocal, temperatures,
     &           temperatures_ref, local_debug
c
      data zero / 0.0d00 /
c
      local_debug       = .false.
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
      kthread           = local_work%now_thread
c
      time_n = local_work%time_n
      dtime  = local_work%dt
c
      nonlocal = local_work%is_cohes_nonlocal
      imxvl    = mxvl   ! protect mxvl
      if( gpn .eq. 1 ) then
!DIR$ VECTOR ALIGNED
        local_work%elem_hist1 = zero
      end if
c
c              get updated stresses and new [D]
c              - usual process for iter >=1 or extrapolated
c                incremental displacements
c              - iter = 0. get linear [D]s. update
c                stresses for user imposed incremental
c                displacements
c
c              model does not have temperature dependence as yet
c
      if( iter .ge. 1 .or. extrapolated_du ) then
         call drive_04_update_std
      else
         call drive_04_update_iter_0_no_extrapolate
         do i = 1, span
           ds1  = uddt(i,1)
           ds2  = uddt(i,2)
           dn   = uddt(i,3)
           tns1 = local_work%urcs_blk_n(i,1,gpn)
           tns2 = local_work%urcs_blk_n(i,2,gpn)
           tnn  = local_work%urcs_blk_n(i,3,gpn)
           local_work%urcs_blk_n1(i,1,gpn) = tns1 +
     &         cep(i,1,1)*ds1 + cep(i,1,2)*ds2 + cep(i,1,3)*dn
           local_work%urcs_blk_n1(i,2,gpn) = tns2 +
     &         cep(i,2,1)*ds1 + cep(i,2,2)*ds2 + cep(i,2,3)*dn
           local_work%urcs_blk_n1(i,3,gpn) = tnn +
     &         cep(i,3,1)*ds1 + cep(i,3,2)*ds2 + cep(i,3,3)*dn
          end do
      end if
c
c              store the lower-triangle of the 3x3 symmetric [D] each
c              element. subroutine should be inlined
c
      call rstgp1_store_cep( span, mxvl, gpn,
     &                       gbl_cep_blocks(now_blk)%vector,
     &                       6, cep )
c
      if( local_debug ) then
         write(iout,9000)
         write(iout,*) " .... extrapolated_du: ",extrapolated_du
         do i = 1, span
           write(iout,9100) felem+i-1, gpn
           write(iout,9110) cep(i,1:3,1:3)
         end do
      end if
c
      return
c
9000  format(5x,'drive_04_update. returned [D]s')
9100  format(10x,'...element, gpn: ',i7,i3)
9110  format(3(15x,3f10.3,/))
c
       contains
c      ========
c     ****************************************************************
c     *                                                              *
c     *         subroutine drive_04_update_iter_0_no_extrapolate     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 10/26/2015 rhd             *
c     *                                                              *
c     *     drive model 04 interface-cohesive to                     *
c     *     update stresses and history for all elements in the      *
c     *     for gauss point gpn                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_04_update_iter_0_no_extrapolate
c
c                linear and non-linear cohesive model
c
      call lcnst4(
     &   step, local_work%cohes_type, span, mxvl, nstr, gpn, felem,
     &   iout, time_n, dtime, cep, local_work%intf_prp_block )
c
      return
      end subroutine drive_04_update_iter_0_no_extrapolate

c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_04_update_std               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 10/26/2015 rhd             *
c     *                                                              *
c     *     drive model 04 interface-cohesive to                     *
c     *     update stresses and history for all elements in the      *
c     *     for gauss point gpn                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_04_update_std
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
        call cnst4(
     &   step, iter, felem, gpn, iout, span, imxvl, time_n, dtime,
     &   nonlocal, knumthreads, kthread, local_work%cohes_type,
     &   local_work%intf_prp_block,
     &   local_work%ddtse(1,1,gpn),
     &   local_work%elem_hist(1,1,gpn),
     &   local_work%elem_hist1(1,1,gpn),
     &   cep,
     &   local_work%cohes_temp_ref(1),
     &   local_work%cohes_dtemp(1),
     &   local_work%cohes_temp_n(1),
     &   local_work%top_surf_solid_elements(1),
     &   local_work%bott_surf_solid_elements(1),
     &   local_work%top_surf_solid_stresses_n(1,1),
     &   local_work%bott_surf_solid_stresses_n(1,1),
     &   local_work%top_surf_solid_eps_n(1,1),
     &   local_work%bott_surf_solid_eps_n(1,1),
     &   local_work%nonlocal_stvals_top_n(1,1),
     &   local_work%nonlocal_stvals_bott_n(1,1),
     &   local_work%top_solid_matl(1),
     &   local_work%bott_solid_matl(1) )
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
        call cnst4(
     &   step, iter, felem, gpn, iout, span, imxvl, time_n, dtime,
     &   nonlocal, knumthreads, kthread, local_work%cohes_type,
     &   local_work%intf_prp_block,
     &   local_work%ddtse(1,1,gpn),
     &   local_work%elem_hist(1,1,gpn),
     &   local_work%elem_hist1(1,1,gpn),
     &   cep,
     &   local_work%cohes_temp_ref(1),
     &   local_work%cohes_dtemp(1),
     &   local_work%cohes_temp_n(1),
     &   idummy(1),
     &   idummy(1),
     &   ddummy(1),
     &   ddummy(1),
     &   ddummy(1),
     &   ddummy(1),
     &   ddummy(1),
     &   ddummy(1),
     &   idummy(1),
     &   idummy(1) )
      end if
c
      return
      end subroutine drive_04_update_std
      end subroutine drive_04_update
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_05_update                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/26/2017 rhd              *
c     *                                                              *
c     *     drives material model 05 (cyclic plastcity) to           *
c     *     update stresses and history for all elements in the      *
c     *     for gauss point gpn                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_05_update( gpn, props, lprops, iprops,
     &                            local_work, uddt_displ, iout )
      use main_data, only : matprp, lmtprp
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks,
     &                            nonlocal_flags, nonlocal_data_n1
      use main_data, only : extrapolated_du, non_zero_imposed_du
c
      implicit none
      include 'param_def'
c
c                      parameter declarations
c
      real    ::  props(mxelpr,*)   ! all same but read only
      logical :: lprops(mxelpr,*)
      integer :: iprops(mxelpr,*)
      integer :: gpn, iout
      double precision ::  uddt_displ(mxvl,nstr)
      include 'include_sig_up'
c
c
c                       locally defined variables
c
      integer :: span, felem, type, order, ngp, nnode, ndof, step,
     &           iter, now_blk, mat_type, number_points, curve_set,
     &           hist_size_for_blk, curve_type, elem_type, i,
     &           matnum
c
      double precision ::
     &  gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero, gp_alpha, dtime, sig_tol, ddummy,
     &  nh_sigma_0_vec(mxvl), nh_q_u_vec(mxvl), nh_b_u_vec(mxvl),
     &  nh_h_u_vec(mxvl), nh_gamma_u_vec(mxvl), gp_tau_vec(mxvl),
     &  uddt_temps(mxvl,nstr), uddt(mxvl,nstr), cep(mxvl,6,6)
c
      logical :: signal_flag, local_debug, temperatures,
     &           temperatures_ref, adaptive_possible, geonl,
     &           cut_step_size_now, fgm_enode_props,
     &           segmental, nonlin_hard, generalized_pl
      data local_debug, zero / .false., 0.0d00 /
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
      ndof              = local_work%num_enode_dof
      geonl             = local_work%geo_non_flg
      order             = local_work%int_order
      mat_type          = local_work%mat_type
      nnode             = local_work%num_enodes
      signal_flag       = local_work%signal_flag
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      segmental         = local_work%segmental
      number_points     = local_work%number_points
      curve_set         = local_work%curve_set_number
      fgm_enode_props   = local_work%fgm_enode_props
      adaptive_possible = local_work%adaptive_flag .and.
     &                    step .gt. 1
      cut_step_size_now = .false.
      hist_size_for_blk = local_work%hist_size_for_blk
c
      local_debug       = .false. ! felem .eq. 1 .and. gpn .eq. 3
      if( local_debug ) then
        write(iout,9000) felem, gpn, span
        write(iout,9010) dtime, type, order, nnode, ndof, geonl, step,
     &                   iter, now_blk, mat_type,
     &                   temperatures, temperatures_ref, segmental,
     &                   number_points, curve_set,
     &                   fgm_enode_props, hist_size_for_blk
      end if

c
c          determine if the material elastic and properties are
c          temperature dependent. The temperature dependent props
c          are e, nu, alpha, gp_sigma_0, gp_h_u, gp_beta_u,
c          gp_delta_u. Points on a curve are defined in to
c          data to satisfy all error checks but are not used here.
c
      curve_type = -1
      if( segmental ) call set_segmental_type( curve_set, curve_type,
     &                                         local_work%eps_curve )
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
      if( local_debug ) write(iout,9030)
c
c          get temperature dependent young's modulus, poisson's
c          ratio, alpha, gp_sigma_0, gp_h_u, gp_beta_u, gp_delta_u
c          for the temperature at n and n+1
c
      if( curve_type .eq. 1 ) then
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
         if( local_debug ) write(iout,9040)
      end if
c
c            get the thermal strain increment (actually negative
c            of increment)
c
!DIR$ VECTOR ALIGNED
      uddt_temps = zero
      if( temperatures ) then
        call gp_temp_eps( span, uddt_temps, local_work%alpha_vec,
     &                    gp_dtemps, gp_temps, gp_rtemps,
     &                    local_work%alpha_vec_n, type )
        if( local_debug ) write(iout,9050)
      end if
c
c            uddt_displ - strain increment due to displacement increment
c            uddt_temps - (negative) of strain increment just due
c                         to temperature change
c
!DIR$ VECTOR ALIGNED
      uddt = uddt_displ + uddt_temps
!DIR$ VECTOR ALIGNED
      cep  = zero
c
      do i = 1, span
       if( local_work%killed_status_vec(i) ) uddt(i,1:nstr) = zero
      end do
c
c            now standard update process. use nonlinear update and [D]
c            computation for iter = 0 and extrapolation or iter > 1
c            for iter = 0 and no extrapolation, use linear-elastic [D]
c            with props at n+1.

      if( iter >= 1 .or. extrapolated_du ) then !nonlinear update
       if( local_debug ) write(iout,9060)
       call drive_05_update_c
       if ( adaptive_possible .and. cut_step_size_now ) then
          local_work%material_cut_step = .true.
          return
       end if
       if( local_work%block_has_nonlocal_solids )
     &           call drive_05_update_nonlocal
       call cnst5( span, felem, gpn, iter, iout, mxvl, nstr,
     &            local_work%e_vec, local_work%nu_vec,
     &            local_work%mm05_props,
     &            local_work%rtse(1,1,gpn),
     &            local_work%elem_hist(1,1,gpn),
     &            local_work%elem_hist1(1,1,gpn),
     &            local_work%urcs_blk_n1(1,1,gpn),
     &            cep,
     &            local_work%gp_h_u_vec, local_work%gp_beta_u_vec,
     &            local_work%gp_delta_u_vec,
     &            gp_tau_vec )
      else  ! linear-elastic update
        if( local_debug ) write(iout,9070)
        call drive_05_update_a
      end if
c
c          save the [D] matrices (lower-triangle)
c
      call rstgp1_store_cep( span, mxvl, gpn,
     &         gbl_cep_blocks(now_blk)%vector, 21, cep )
      if( local_debug ) write(iout,9080)
c
      return
c
 9000 format(1x,'.... debug mm05. felem, gpn, span: ',i7,i3,i3)
 9010 format(10x,'...dtime, type, order, nnode, ndof:',e14.6,4i5,
     &     /,10x,'...geonl, step, iter, now_blk, mat_type: ',l2,4i5,
     &     /,10x,'...temperatures, temperatures_ref: ',
     &               2l2,
     &     /,10x,'...segmental, number_points, curve_set: ',l2,i3,i3,
     &     /,10x,'...fgm_enode_props, hist_size_for_blk: ',
     &    l3,i4 )
 9030 format(10x,'...temperatures computed at integration point...')
 9040 format(10x,'...temperatures dependent properties computed...')
 9050 format(10x,'...thermal strains computed...')
 9060 format(10x,'...update stresses nonlinear procedure...' )
 9070 format(10x,'...update stresses use linear [D]...' )
 9080 format(10x,'...[D]s saved to global structure...')
c
      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_05_update_c                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/31/2015 rhd             *
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_05_update_c
      implicit none
c
c            build local data vectors for nonlinear_hardening option
c            to maintain consistency with generalized_plasticity option.
c            the gp option has possibly temperature dependent properties.
c            the nh option could have e, nu and alpha defined as temperature
c            dependent by the user thru curves (but we really don't show that
c            option in the manual.
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
        nh_sigma_0_vec(i) = local_work%sigyld_vec(i)
        nh_q_u_vec(i)     = matprp(55,matnum)
        nh_b_u_vec(i)     = matprp(56,matnum)
        nh_h_u_vec(i)     = matprp(57,matnum)
        nh_gamma_u_vec(i) = matprp(59,matnum)
        gp_tau_vec(i)     = matprp(56,matnum)
      end do
c
      nonlin_hard      = matprp(58,matnum) .gt. zero
      generalized_pl   = matprp(58,matnum) .lt. zero
      sig_tol          = matprp(60,matnum)
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
c     In local_work: same for both nonlinear_hardening (FA) and generalized_plas
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
     &           local_work%elem_hist1(1,1,gpn),
     &           local_work%block_has_nonlocal_solids,
     &           local_work%nonlocal_state_blk(1,1),
     &           nonlocal_shared_state_size ) ! value in param_def
c
      return
      end subroutine drive_05_update_c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_05_update_nonlocal          *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/26/2017 rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_05_update_nonlocal
      implicit none
c
      integer :: i, n, elem_num
      double precision :: real_npts

      n = nonlocal_shared_state_size ! for convenience from param_def
      if( local_debug ) write(iout,9010) n
c
      if( gpn .eq. 1 ) then  ! zero global values for elements
        do i = 1, span
          elem_num = felem + i - 1
          if( nonlocal_flags(elem_num) )
     &         nonlocal_data_n1(elem_num)%state_values(1:n) = zero
        end do
      end if
c
      do i = 1, span ! add in this gpn nonlocal values
       elem_num = felem + i - 1
       if( nonlocal_flags(elem_num) )
     &       nonlocal_data_n1(elem_num)%state_values(1:n) =
     &       nonlocal_data_n1(elem_num)%state_values(1:n) +
     &       local_work%nonlocal_state_blk(i,1:n)
      end do
c
      if( gpn .eq. local_work%num_int_points ) then
         real_npts = dble( local_work%num_int_points )
         do i = 1, span
          elem_num = felem + i - 1
          if( nonlocal_flags(elem_num) )  then
            nonlocal_data_n1(elem_num)%state_values(1:n) =
     &      nonlocal_data_n1(elem_num)%state_values(1:n) / real_npts
          end if
         end do
      end if
c
      return
c
 9010 format(/,'      processing nonlocal values. # values: ',i2 )
c
      end subroutine drive_05_update_nonlocal


c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_05_update_a                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/31/2015 rhd             *
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_05_update_a
      implicit none
c
      integer :: k, m, i
      double precision :: one, two, e, nu, c1, c2, c3, c4
      data one, two / 1.0d00, 2.0d00 /
c
c              get linear-elastic [D] with potentially temperature
c              dependent properties
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
         if( local_work%killed_status_vec(i) ) cycle
         e  = local_work%e_vec(i)
         nu = local_work%nu_vec(i)
         c1 = (e/((one+nu)*(one-two*nu)))
         c2 = (one-nu)*c1
         c3 = ((one-two*nu)/two)*c1
         c4 = nu*c1
         cep(i,1,1)= c2
         cep(i,2,2)= c2
         cep(i,3,3)= c2
         cep(i,4,4)= c3
         cep(i,5,5)= c3
         cep(i,6,6)= c3
         cep(i,1,2)= c4
         cep(i,1,3)= c4
         cep(i,2,1)= c4
         cep(i,3,1)= c4
         cep(i,2,3)= c4
         cep(i,3,2)= c4
      end do
c
c              stresses at n+1 using linear-elastic [D]
c
       call drive_05_update_b( span, mxvl, uddt, cep,
     &                        local_work%urcs_blk_n(1,1,gpn),
     &                        local_work%urcs_blk_n1(1,1,gpn),
     &                        local_work%killed_status_vec )
c
      return
      end subroutine drive_05_update_a
c
      end subroutine drive_05_update
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_05_update_b                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/21/2015 rhd             *
c     *                                                              *
c     *     support routine for mm05 material driver.                *
c     *     should be inlined                                        *
c     *                                                              *
c     ****************************************************************

      subroutine drive_05_update_b( span, mxvl, uddt,
     &                              local_cep, stress_n, stress_np1,
     &                              killed_status )
      implicit none
c
      integer :: span, mxvl
      logical :: killed_status(*)
      double precision :: local_cep(mxvl,6,6), stress_n(mxvl,6),
     &                    stress_np1(mxvl,6), uddt(mxvl,6)
      double precision, parameter :: zero = 0.0d00
c
      integer i, k, m
c
c              for each element in block, update stresses by
c              [D-elastic] * uddt. uddt contains thermal increment +
c              increment from imposed nodal displacements
c
!DIR$ VECTOR ALIGNED
      stress_np1 = stress_n
c
      do k = 1, 6
       do m = 1, 6
!DIR$ VECTOR ALIGNED
         do i = 1, span
           stress_np1(i,k) = stress_np1(i,k) +
     &                       local_cep(i,m,k) * uddt(i,m)
         end do
       end do
      end do
c
      do i = 1, span
        if( killed_status(i) ) stress_np1(i,1:6) = zero
      end do
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *            subroutine drive_06_update  (creep)               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 10/18/2015 rhd             *
c     *                                                              *
c     *     drives material model 06 to update stresses and history  *
c     *     for all elements in the block for gauss point gpn        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_06_update( gpn, props, lprops, iprops,
     &                            local_work, uddt_displ, iout )
      use main_data, only : matprp, lmtprp, extrapolated_du,
     &                      non_zero_imposed_du
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks,
     &                            nonlocal_flags, nonlocal_data_n1
c
      implicit integer (a-z)
      include 'param_def'
c
c                      parameter declarations
c
      real     :: props(mxelpr,*)   ! all 3 are same but read-only
      logical  :: lprops(mxelpr,*)
      integer  :: iprops(mxelpr,*)
      double precision :: uddt_displ(mxvl,nstr)
      include 'include_sig_up'
c
c                       locally defined variables
c
      double precision ::
     &  gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero, gp_alpha, dtime, real_npts, uddt_temps(mxvl,nstr),
     &  uddt(mxvl,nstr), cep(mxvl,6,6)
c
      logical :: signal_flag, local_debug, temperatures,
     &           temperatures_ref, process, compute_creep_strains
      data zero /  0.0d00 /
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
      local_debug = .false.
      if( local_debug ) then
        write(iout,*) '... entered drive_06_update'
        write(iout,9000) now_blk, felem, span, gpn
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
c            compute (negative) of thermal strain increment
c
!DIR$ VECTOR ALIGNED
      uddt_temps = zero
      if( temperatures ) then
        call gp_temp_eps( span, uddt_temps, local_work%alpha_vec,
     &                    gp_dtemps, gp_temps, gp_rtemps,
     &                    local_work%alpha_vec, type )
      end if
c
c            init block of nonlocal state variables. values array always
c            allocated but with size (1,1) for std. local analyses.
c            just makes passing args simpler.
c
      if( local_work%block_has_nonlocal_solids )
     &    local_work%nonlocal_state_blk = zero ! array
c
c            uddt_displ - strain increment due to displacement increment
c            uddt_temps - (negative) of strain increment just due
c                         to temperature change
c
c            iterno > 0 : add the two parts. give to mm06
c            iterno = 0 :  and global extrapolate ...
c              on: - add, use updated stresses computed by mm06
c                  - return updated [D]s
c             off: - pass only uddt_temps to mm06
c                  - mm06 will include creep strain increment in uddt
c                  - return linear-elastic [D]s
c
      compute_creep_strains = .false.
      if( iter .ge. 1 ) uddt = uddt_displ + uddt_temps
      if( iter .eq. 0 ) then
        if( extrapolated_du ) then
!DIR$ VECTOR ALIGNED
           uddt = uddt_displ + uddt_temps
        else !  iter = 0, no extrapolate
!DIR$ VECTOR ALIGNED
           uddt = uddt_temps
           compute_creep_strains = .true.
        end if
      end if
c
      call mm06( step, iter, felem, gpn, mxvl, hist_size_for_blk,
     &           nstrs, nstr, span, iout, dtime,
     &           local_work%mm06_props,
     &           local_work%e_vec,
     &           local_work%nu_vec,
     &           local_work%n_power_vec, local_work%rtse(1,1,gpn),
     &           local_work%urcs_blk_n(1,1,gpn),
     &           local_work%urcs_blk_n1(1,1,gpn),
     &           uddt, local_work%elem_hist(1,1,gpn),
     &           local_work%elem_hist1(1,1,gpn),
     &           local_work%killed_status_vec,
     &           local_work%block_has_nonlocal_solids,
     &           local_work%nonlocal_state_blk(1,1),
     &           nonlocal_shared_state_size, ! value in param_def
     &           cep, compute_creep_strains )
c
c          save the [D] matrices (lower-triangle)
c          computed by mm06. see comments above for special
c          computations with iter=0 and no extrapolation.
c
      call rstgp1_store_cep( span, mxvl, gpn,
     &         gbl_cep_blocks(now_blk)%vector, 21, cep )
c
      if( iter == 0 ) then
        if( extrapolated_du ) go to 9999 ! all updated
!DIR$ VECTOR ALIGNED
        uddt = uddt + uddt_displ ! temps + creep strain incr + imposed du
        call drive_06_update_a( span, mxvl, uddt, cep,
     &                          local_work%urcs_blk_n1(1,1,gpn) )
        go to 9999
      end if
c
      if( .not. local_work%block_has_nonlocal_solids ) go to 9999
c
      n = nonlocal_shared_state_size ! for convenience from param_def
      if( local_debug ) write(iout,9010) n
c
      if( gpn .eq. 1 ) then  ! zero global values for elements
        do i = 1, span
          elem_num = felem + i - 1
          if( nonlocal_flags(elem_num) )
     &         nonlocal_data_n1(elem_num)%state_values(1:n) = zero
        end do
      end if
c
      do i = 1, span ! add in this gpn nonlocal values
       elem_num = felem + i - 1
       if( nonlocal_flags(elem_num) )
     &       nonlocal_data_n1(elem_num)%state_values(1:n) =
     &       nonlocal_data_n1(elem_num)%state_values(1:n) +
     &       local_work%nonlocal_state_blk(i,1:n)
      end do
c
      if( gpn .eq. local_work%num_int_points ) then
         real_npts = dble( local_work%num_int_points )
         do i = 1, span
          elem_num = felem + i - 1
          if( nonlocal_flags(elem_num) )  then
            nonlocal_data_n1(elem_num)%state_values(1:n) =
     &      nonlocal_data_n1(elem_num)%state_values(1:n) / real_npts
          end if
         end do
      end if

 9999 continue
      if( local_debug ) write(iout,*) '... leave drive_06_update'
      return
c
 9000 format(/,'.... calling mm06. blk, felem, span, gpn: ',4i7)
 9010 format(/,'      processing nonlocal values. # values: ',i2 )
c
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_06_update_a                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 10/18/2015                 *
c     *                                                              *
c     *     support routine for creep material driver.               *
c     *     should be inlined                                        *
c     *                                                              *
c     ****************************************************************

      subroutine drive_06_update_a( span, mxvl, uddt,
     &                              local_cep, stress_np1 )
      implicit none
c
      integer :: span, mxvl
      double precision ::
     &  local_cep(mxvl,6,6), stress_np1(mxvl,6),
     &  uddt(mxvl,6)
c
      integer i, k, m
c
c              for each element in block, update stresses by
c              [D-elastic] * uddt. uddt contains creep increment +
c              thermal increment + increment from imposed
c              nodal displacements
c
      do k = 1, 6
       do m = 1, 6
!DIR$ VECTOR ALIGNED
         do i = 1, span
           stress_np1(i,k) = stress_np1(i,k) +
     &                       local_cep(i,m,k) * uddt(i,m)
         end do
       end do
      end do
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_07_update                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/2/2016 rhd               *
c     *                                                              *
c     *     drives material model 07 (mises + hydrogen) to           *
c     *     update stresses and history for all elements in the      *
c     *     for gauss point gpn                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_07_update( gpn, props, lprops, iprops,
     &                            local_work, uddt_displ, iout )
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
      use main_data, only : extrapolated_du, non_zero_imposed_du
c
      implicit none
      include 'param_def'
c
c                      parameter declarations
c
      real ::    props(mxelpr,*)   ! all same but read only
      logical :: lprops(mxelpr,*)
      integer :: iprops(mxelpr,*)
      integer :: gpn, iout
      double precision ::  uddt_displ(mxvl,nstr)
      include 'include_sig_up'
c
c                       locally defined variables
c
      integer :: span, felem, type, order, ngp, nnode, ndof, step,
     &           iter, now_blk, mat_type, number_points, curve_set,
     &           hist_size_for_blk, curve_type, elem_type, i
c
      double precision ::
     &  dtime, gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero, gp_alpha,  uddt_temps(mxvl,nstr),
     &  uddt(mxvl,nstr), cep(mxvl,6,6)
c
      logical :: geonl, local_debug, temperatures, segmental,
     &           temperatures_ref, fgm_enode_props, signal_flag,
     &           adaptive_possible, cut_step_size_now
c
      data zero / 0.0d0 /
c
      dtime             = local_work%dt
      span              = local_work%span
      felem             = local_work%felem
      step              = local_work%step
      iter              = local_work%iter
      now_blk           = local_work%blk
      type              = local_work%elem_type
      order             = local_work%int_order
      ndof              = local_work%num_enode_dof
      geonl             = local_work%geo_non_flg
      nnode             = local_work%num_enodes
      mat_type          = local_work%mat_type
      signal_flag       = local_work%signal_flag
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      hist_size_for_blk = local_work%hist_size_for_blk
      fgm_enode_props   = local_work%fgm_enode_props
      adaptive_possible = local_work%allow_cut
      local_debug       = .false.
c
      if( local_debug ) then
        write(iout,9000) felem, gpn, span
        write(iout,9010) dtime, type, order, nnode, ndof, geonl, step,
     &                   iter, now_blk, mat_type,
     &                   temperatures, temperatures_ref,
     &                   fgm_enode_props, hist_size_for_blk
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
c            subtract out the thermal strain increment from uddt (the
c            strain increment for step)
c
!DIR$ VECTOR ALIGNED
      uddt_temps = zero
      if ( temperatures ) then
        call gp_temp_eps( span, uddt_temps, local_work%alpha_vec,
     &                    gp_dtemps, gp_temps, gp_rtemps,
     &                    local_work%alpha_vec, type )
      end if
c
!DIR$ VECTOR ALIGNED
      uddt = uddt_displ + uddt_temps
!DIR$ VECTOR ALIGNED
      cep  = zero
c
      do i = 1, span
       if( local_work%killed_status_vec(i) ) uddt(i,1:nstr) = zero
      end do
c
      if( iter >= 1 .or. extrapolated_du ) then !nonlinear update
       if( local_debug ) write(iout,9060)
       cut_step_size_now = .false.
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
       local_work%material_cut_step = cut_step_size_now
       if( cut_step_size_now ) return
       call cnst7( span, felem, gpn, iter, iout, mxvl, nstr,
     &            local_work%e_vec, local_work%nu_vec,
     &            local_work%sigyld_vec, local_work%n_power_vec,
     &            local_work%mm07_props,
     &            local_work%rtse(1,1,gpn),
     &            local_work%elem_hist(1,1,gpn),
     &            local_work%elem_hist1(1,1,gpn),
     &            local_work%urcs_blk_n1(1,1,gpn), cep )
      else  ! linear-elastic update
        if( local_debug ) write(iout,9070)
        call drive_07_update_a
      end if
c
c          save the [D] matrices (lower-triangle)
c
      call rstgp1_store_cep( span, mxvl, gpn,
     &         gbl_cep_blocks(now_blk)%vector, 21, cep )
      if( local_debug ) write(iout,9080)
c
      return
c
 9000 format(1x,'.... debug mm07. felem, gpn, span: ',i7,i3,i3)
 9010 format(10x,'...dtime, type, order, nnode, ndof:',e14.6,4i5,
     &     /,10x,'...geonl, step, iter, now_blk, mat_type: ',l2,4i5,
     &     /,10x,'...temperatures, temperatures_ref: ',
     &               2l2,
     &     /,10x,'...segmental, number_points, curve_set: ',l2,i3,i3,
     &     /,10x,'...fgm_enode_props, hist_size_for_blk: ',
     &    l3,i4 )
 9610 format(' >> rate iterations to converge: ',i3 )
 9020 format(10x,'...fgm properties determined...')
 9030 format(10x,'...temperatures computed at integration point...')
 9040 format(10x,'...temperatures dependent properties computed...')
 9050 format(10x,'...thermal strains computed...')
 9060 format(10x,'...update stresses nonlinear procedure...' )
 9070 format(10x,'...update stresses use linear [D]...' )
 9080 format(10x,'...[D]s saved to global structure...')
c
      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_07_update_a                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/2/2016 rhd               *
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_07_update_a
      implicit none
c
      integer :: k, m, i
      double precision :: one, two, e, nu, c1, c2, c3, c4
      data one, two / 1.0d00, 2.0d00 /
c
c              get linear-elastic [D] with potentially temperature
c              dependent properties
c
!DIR$ VECTOR ALIGNED
      do i = 1, span
         if( local_work%killed_status_vec(i) ) cycle
         e  = local_work%e_vec(i)
         nu = local_work%nu_vec(i)
         c1 = (e/((one+nu)*(one-two*nu)))
         c2 = (one-nu)*c1
         c3 = ((one-two*nu)/two)*c1
         c4 = nu*c1
         cep(i,1,1)= c2
         cep(i,2,2)= c2
         cep(i,3,3)= c2
         cep(i,4,4)= c3
         cep(i,5,5)= c3
         cep(i,6,6)= c3
         cep(i,1,2)= c4
         cep(i,1,3)= c4
         cep(i,2,1)= c4
         cep(i,3,1)= c4
         cep(i,2,3)= c4
         cep(i,3,2)= c4
      end do
c
c              stresses at n+1 using linear-elastic [D]
c
       call drive_07_update_b( span, mxvl, uddt, cep,
     &                        local_work%urcs_blk_n(1,1,gpn),
     &                        local_work%urcs_blk_n1(1,1,gpn),
     &                        local_work%killed_status_vec )
c
      return
      end subroutine drive_07_update_a
c
      end subroutine drive_07_update
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_07_update_b                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/2/2016 rhd               *
c     *                                                              *
c     *     support routine for mm07 material driver.                *
c     *     should be inlined                                        *
c     *                                                              *
c     ****************************************************************

      subroutine drive_07_update_b( span, mxvl, uddt,
     &                              local_cep, stress_n, stress_np1,
     &                              killed_status )
      implicit none
c
      integer :: span, mxvl
      logical :: killed_status(*)
      double precision ::
     &  local_cep(mxvl,6,6), stress_n(mxvl,6), stress_np1(mxvl,6),
     &  uddt(mxvl,6), zero
      data zero / 0.0d00 /
c
      integer i, k, m

c
c              for each element in block, update stresses by
c              [D-elastic] * uddt. uddt contains thermal increment +
c              increment from imposed nodal displacements
c
!DIR$ VECTOR ALIGNED
      stress_np1 = stress_n
c
      do k = 1, 6
       do m = 1, 6
!DIR$ VECTOR ALIGNED
         do i = 1, span
           stress_np1(i,k) = stress_np1(i,k) +
     &                       local_cep(i,m,k) * uddt(i,m)
         end do
       end do
      end do
c
      do i = 1, span
        if( killed_status(i) ) stress_np1(i,1:6) = zero
      end do
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *            subroutine drive_umat_update  (umat)              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/19/2017 rhd             *
c     *                                                              *
c     *     drives material model 08 (Abaqus umat) to update         *
c     *     stresses and history for all elements in block at        *
c     *     integration point gpn                                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_umat_update( gpn, local_work, uddt, qn1, iout )
      use main_data, only : matprp, lmtprp, nonlocal_analysis,
     &                      extrapolated_du, non_zero_imposed_du
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : nonlocal_flags, nonlocal_data_n1,
     &                            gbl_cep_blocks => cep_blocks
c
      implicit none
      include 'param_def'
c
c                      parameter declarations
c
      integer :: gpn, iout
      double precision :: uddt(mxvl,nstr), qn1(mxvl,nstr,nstr)
      include 'include_sig_up'
c
c                       locally defined variables
c
      integer :: span, felem, type, order, ngp, nnode, step, iter,
     &           now_blk, hist_size_for_blk, knumthreads, kthread,
     &           max_nstatv, ndi, nshr, ntens, npt, layer, kspt,
     &           kstep, kinc, kiter, kout,nstatv, k, ielem, noel,
     &           nprops, j, nj, start_loc, n
      double precision ::
     &  gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero, one, gp_alpha, ddsddt(6), drplde(6), drpldt,
     &  big, pnewdt, predef(1), dpred(1), time(2), dtime,
     &  stress(6), stran(6), dstran(6),
     &  abq_props(50), temp_n, dtemp,
     &  statev(500), sse, spd, scd, coords(3), celent, rpl,
     &  dfgrd0(9), dfgrd1(9), dfgrd0_array(3,3), dfgrd1_array(3,3),
     &  drot(9), ddsdde(6,6),
     &  gp_coords(mxvl,3), symm_part_ddsdde(21), total_work_np1,
     &  plastic_work_np1, identity(9), check_key, t5, t6,
     &  temp_0, temp_n_0, s1, s2, ps(3), an(3,3),
     &  unrotated_cauchy(mxvl,6), real_npts,
     &  nonloc_ele_values(nonlocal_shared_state_size),
     &  sys_vals(nonlocal_shared_state_size), dstran_temps_only(6),
     &  dstran_displ_only(6), uddt_temps(mxvl,nstr)
c
      equivalence (dfgrd0, dfgrd0_array),  (dfgrd1, dfgrd1_array)
c
      logical :: signal_flag, local_debug, debug_now, temperatures,
     &        temperatures_ref, init_sig_eps, init_history,
     &        chk_umat_support, chk, chk2, do_nonlocal,
     &        process_flag
      integer :: map(6)
      character(len=8) :: cmname
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
      kthread           = local_work%now_thread
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
!DIR$ VECTOR ALIGNED
       local_work%urcs_blk_n = zero
!DIR$ VECTOR ALIGNED
       local_work%strain_n = zero
      end if
      if( init_history ) then
!DIR$ VECTOR ALIGNED
         local_work%elem_hist = zero
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
c           3. get negative of temp increments for step
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
!DIR$ VECTOR ALIGNED
      uddt_temps = zero ! replaced by -alpha * delta T
      if( temperatures ) ! global flag set by loads processor
c                           to indicate user-specified temp changes over
c                           load step
     &    call gp_temp_eps( span, uddt_temps, local_work%alpha_vec,
     &                      gp_dtemps, gp_temps, gp_rtemps,
     &                      local_work%alpha_vec, type )
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

!DIR$ VECTOR ALIGNED
      do k = 1, 9
       dfgrd0(k) = identity(k)
       dfgrd1(k) = identity(k)
       drot(k)   = identity(k)
      end do

c
      rpl         = zero
!DIR$ VECTOR ALIGNED
      ddsddt(1:6) = zero
!DIR$ VECTOR ALIGNED
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
c                   increment. adjust for thermal components and
c                   swap order of xz, yz shear terms to match Abaqus.
c                   thermal part of incremental strain handled above.
c                   WARP3D processing here for strain at n requires
c                   temperature invariant CTEs (subtract reference
c                   temperature). User should not define "alpha" in
c                   WARP3D input if umat handles temperature effects.
c                   alpha below will then be zero.
c
      do j = 1, 6
       nj = map(j)
       stress(nj) = local_work%urcs_blk_n(ielem,j,gpn)
       stran(nj)  = local_work%strain_n(ielem,j,gpn) -
     &                  local_work%alpha_vec(ielem,j) * temp_n_0
       dstran_temps_only(nj) = uddt_temps(ielem,j)
       dstran_displ_only(nj) = uddt(ielem,j)
       dstran(nj) = uddt_temps(ielem,j) + uddt(ielem,j)
      end do
c
      process_flag = .false.
      if( kiter .eq. 0 .and. .not. extrapolated_du ) then
        dstran(1:6) =  dstran_temps_only(1:6)
        process_flag = .true.
      end if

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
!DIR$ VECTOR ALIGNED
      ddsdde = zero
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
     k   noel, npt, layer,
     l   kspt, kstep, kinc, kiter, kout, kthread, knumthreads,
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
c                   imposed displacements, temperature increments,
c                   and creep strain increments (+ other "initial"
c                   strain effects from the umat.
c                   here we store the [D] returned from the umat
c                   (could be linear [D]). add stress increment from
c                   imposed non-zero displacement increment for the
c                   load step - we do not want umat to see those values
c                   which are often acting on just a few elements
c                   no need to store histories, nonlocal etc since this
c                   just to estimate load increment for step for
c                   other than directly applied forces/pressures
c
      if( kiter .eq. 0 ) then
        call rstgp1_make_symmetric_store( ddsdde, symm_part_ddsdde )
        start_loc = ( 21 * span * (gpn-1) ) + 21 * (ielem-1)
        do k = 1, 21
         gbl_cep_blocks(now_blk)%vector(start_loc+k) =
     &           symm_part_ddsdde(k)
        end do
c
        if( process_flag .and. non_zero_imposed_du )
     &     stress = stress + matmul( ddsdde, dstran_displ_only )
c
        local_work%urcs_blk_n1(ielem,1:4,gpn) = stress(1:4)
        local_work%urcs_blk_n1(ielem,5,gpn)   = stress(6)
        local_work%urcs_blk_n1(ielem,6,gpn)   = stress(5)
        if( debug_now ) write(iout,9106) stress(1:6)
        cycle ! to process next element in block
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
      if( local_work%allow_cut .and.  (pnewdt .lt. one) )
     &    local_work%material_cut_step = .true.
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
c     *                 subroutine drive_09_update                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/2/2016 rhd               *
c     *                                                              *
c     *     drives material model 09:; <available> to                *
c     *     update stresses and history for all elements in the      *
c     *     for gauss point gpn                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_09_update( gpn, props, lprops, iprops,
     &                            local_work, uddt_displ, iout )
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks
      use main_data, only : extrapolated_du, non_zero_imposed_du
c
      implicit none
      include 'param_def'
c
c                      parameter declarations
c
      real ::    props(mxelpr,*)   ! all same but read only
      logical :: lprops(mxelpr,*)
      integer :: iprops(mxelpr,*)
      integer :: gpn, iout
      double precision ::  uddt_displ(mxvl,nstr)
      include 'include_sig_up'
c
c                       locally defined variables
c
      integer :: span, felem, type, order, ngp, nnode, ndof, step,
     &           iter, now_blk, mat_type, number_points, curve_set,
     &           hist_size_for_blk, curve_type, elem_type, i
c
      double precision ::
     &  dtime, gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero, gp_alpha,  uddt_temps(mxvl,nstr),
     &  uddt(mxvl,nstr), cep(mxvl,6,6), weight, dj(128)
c
      logical :: geonl, local_debug, temperatures, segmental,
     &           temperatures_ref, fgm_enode_props, signal_flag,
     &           adaptive_possible, cut_step_size_now
c
      data zero / 0.0d0 /
c
      dtime             = local_work%dt
      span              = local_work%span
      felem             = local_work%felem
      step              = local_work%step
      iter              = local_work%iter
      now_blk           = local_work%blk
      type              = local_work%elem_type
      order             = local_work%int_order
      ndof              = local_work%num_enode_dof
      geonl             = local_work%geo_non_flg
      nnode             = local_work%num_enodes
      mat_type          = local_work%mat_type
      signal_flag       = local_work%signal_flag
      temperatures      = local_work%temperatures
      temperatures_ref  = local_work%temperatures_ref
      hist_size_for_blk = local_work%hist_size_for_blk
      fgm_enode_props   = local_work%fgm_enode_props
      adaptive_possible = local_work%allow_cut
c
      if( local_debug ) then
        write(iout,9000) felem, gpn, span
        write(iout,9010) dtime, type, order, nnode, ndof, geonl, step,
     &                   iter, now_blk, mat_type,
     &                   temperatures, temperatures_ref,
     &                   fgm_enode_props, hist_size_for_blk
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
c            subtract out the thermal strain increment from uddt (the
c            strain increment for step)
c
      uddt_temps = zero
      if ( temperatures ) then
        call gp_temp_eps( span, uddt_temps, local_work%alpha_vec,
     &                    gp_dtemps, gp_temps, gp_rtemps,
     &                    local_work%alpha_vec, type )
      end if
c

      uddt = uddt_displ + uddt_temps
      cep  = zero
c
      do i = 1, span
       if( local_work%killed_status_vec(i) ) uddt(i,1:nstr) = zero
      end do
c
      if( iter >= 1 .or. extrapolated_du ) then !nonlinear update
       if( local_debug ) write(iout,9060)
       cut_step_size_now = .false.
c       call mm09( step, iter, felem, gpn, mxvl,  hist_size_for_blk,
c     &           nstrs, nstr, span, iout,
c     &           signal_flag, adaptive_possible, cut_step_size_now,
c     &           local_work%mm07_props,
c     &           local_work%e_vec, local_work%tan_e_vec,
c     &           local_work%nu_vec, local_work%sigyld_vec,
c     &           local_work%n_power_vec, local_work%rtse(1,1,gpn),
c     &           local_work%urcs_blk_n(1,1,gpn),
c     &           local_work%urcs_blk_n1(1,1,gpn),
c     &           uddt, local_work%elem_hist(1,1,gpn),
c     &           local_work%elem_hist1(1,1,gpn) )
       local_work%material_cut_step = cut_step_size_now
       if( cut_step_size_now ) return
c       call cnst9( span, felem, gpn, iter, iout, mxvl, nstr,
c     &            local_work%e_vec, local_work%nu_vec,
c     &            local_work%sigyld_vec, local_work%n_power_vec,
c     &            local_work%mm07_props,
c     &            local_work%rtse(1,1,gpn),
c     &            local_work%elem_hist(1,1,gpn),
c     &            local_work%elem_hist1(1,1,gpn),
c     &            local_work%urcs_blk_n1(1,1,gpn), cep )
      else  ! linear-elastic update
        if( local_debug ) write(iout,9070)
        call drive_09_update_a
      end if
c
c          save the [D] matrices (lower-triangle)
c
      call rstgp1_store_cep( span, mxvl, gpn,
     &         gbl_cep_blocks(now_blk)%vector, 21, cep )
      if( local_debug ) write(iout,9080)
c
      return
c
 9000 format(1x,'.... debug mm09. felem, gpn, span: ',i7,i3,i3)
 9010 format(10x,'...dtime, type, order, nnode, ndof:',e14.6,4i5,
     &     /,10x,'...geonl, step, iter, now_blk, mat_type: ',l2,4i5,
     &     /,10x,'...temperatures, temperatures_ref: ',
     &               2l2,
     &     /,10x,'...segmental, number_points, curve_set: ',l2,i3,i3,
     &     /,10x,'...fgm_enode_props, hist_size_for_blk: ',
     &    l3,i4 )
 9610 format(' >> rate iterations to converge: ',i3 )
 9020 format(10x,'...fgm properties determined...')
 9030 format(10x,'...temperatures computed at integration point...')
 9040 format(10x,'...temperatures dependent properties computed...')
 9050 format(10x,'...thermal strains computed...')
 9060 format(10x,'...update stresses nonlinear procedure...' )
 9070 format(10x,'...update stresses use linear [D]...' )
 9080 format(10x,'...[D]s saved to global structure...')
c
      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_09_update_a                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/2/2016 rhd               *
c     *                                                              *
c     ****************************************************************
c
      subroutine drive_09_update_a
      implicit none
c
      integer :: k, m, i
      double precision :: one, two, e, nu, c1, c2, c3, c4
      data one, two / 1.0d00, 2.0d00 /
c
c              get linear-elastic [D] with potentially temperature
c              dependent properties
c
      do i = 1, span
         if( local_work%killed_status_vec(i) ) cycle
         e  = local_work%e_vec(i)
         nu = local_work%nu_vec(i)
         c1 = (e/((one+nu)*(one-two*nu)))
         c2 = (one-nu)*c1
         c3 = ((one-two*nu)/two)*c1
         c4 = nu*c1
         cep(i,1,1)= c2
         cep(i,2,2)= c2
         cep(i,3,3)= c2
         cep(i,4,4)= c3
         cep(i,5,5)= c3
         cep(i,6,6)= c3
         cep(i,1,2)= c4
         cep(i,1,3)= c4
         cep(i,2,1)= c4
         cep(i,3,1)= c4
         cep(i,2,3)= c4
         cep(i,3,2)= c4
      end do
c
c              stresses at n+1 using linear-elastic [D]
c
       call drive_09_update_b( span, mxvl, uddt, cep,
     &                        local_work%urcs_blk_n(1,1,gpn),
     &                        local_work%urcs_blk_n1(1,1,gpn),
     &                        local_work%killed_status_vec )
c
      return
      end subroutine drive_09_update_a
c
      end subroutine drive_09_update
c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_09_update_b                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/2/2016 rhd               *
c     *                                                              *
c     *     support routine for mm07 material driver.                *
c     *     should be inlined                                        *
c     *                                                              *
c     ****************************************************************

      subroutine drive_09_update_b( span, mxvl, uddt,
     &                              local_cep, stress_n, stress_np1,
     &                              killed_status )
      implicit none
c
      integer :: span, mxvl
      logical :: killed_status(*)
      double precision ::
     &  local_cep(mxvl,6,6), stress_n(mxvl,6), stress_np1(mxvl,6),
     &  uddt(mxvl,6), zero
      data zero / 0.0d00 /
c
      integer i, k, m
c
c              for each element in block, update stresses by
c              [D-elastic] * uddt. uddt contains thermal increment +
c              increment from imposed nodal displacements
c
      stress_np1 = stress_n
c
      do k = 1, 6
       do m = 1, 6
         do i = 1, span
           stress_np1(i,k) = stress_np1(i,k) +
     &                       local_cep(i,m,k) * uddt(i,m)
         end do
       end do
      end do
c
      do i = 1, span
        if( killed_status(i) ) stress_np1(i,1:6) = zero
      end do
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *     subroutine drive_10_update  (crystal plasticity)         *
c     *                                                              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/8/15                    *
c     *                                                              *
c     *     this subroutine drives material model 10 to              *
c     *     update stresses and history for all elements in the      *
c     *     for gauss point gpn                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_10_update( gpn, props, lprops, iprops,
     &                            local_work, uddt_displ, iout )
      use main_data, only : matprp, lmtprp, imatprp, dmatprp, smatprp,
     &                      extrapolated_du, non_zero_imposed_du
      use segmental_curves, only : max_seg_points
      use elem_block_data, only : gbl_cep_blocks => cep_blocks,
     &                            nonlocal_flags, nonlocal_data_n1
c
      implicit none
      include 'param_def'
c
c                      parameter declarations
c
      integer :: gpn, iout
      real    :: props(mxelpr,*)   ! all 3 are same but read-only here
      logical :: lprops(mxelpr,*)
      integer :: iprops(mxelpr,*)
      double precision :: uddt_displ(mxvl,nstr)
      include 'include_sig_up'
c
c
c                       locally defined variables
c
      integer :: ncrystals, iter, span, felem, step, type, order,
     &           nnode, hist_size_for_blk, now_blk,
     &           i, j, matnum, k, start_loc, m, n, igp
      double precision ::
     &  gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero, one, gp_alpha, dtime, uddt_temps(mxvl,nstr),
     &  uddt(mxvl,nstr), cep(mxvl,6,6), cep_vec(36), tol
c
      logical :: signal_flag, local_debug, temperatures,
     &           temperatures_ref, check_D, iter_0_extrapolate_off
      data zero, one / 0.0d0, 1.0d0 /

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
      matnum            = local_work%matnum
      local_debug       = .false. ! now_blk == 1  .and. gpn .eq. 1
      check_D           = .false.
      if( local_debug ) then
        write(iout,9000) felem, gpn, span
        write(iout,9010) dtime, type, order, nnode, step,
     &                   iter, now_blk,
     &                   temperatures, temperatures_ref,
     &                   hist_size_for_blk
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
c            subtract out the thermal strain increment from uddt (the
c            strain increment for step) (We're actually going to do this
c            internally)
c
!DIR$ VECTOR ALIGNED
      uddt_temps = zero
      if ( temperatures ) then
        call gp_temp_eps( span, uddt_temps, local_work%alpha_vec,
     &                    gp_dtemps, gp_temps, gp_rtemps,
     &                    local_work%alpha_vec, type )
      end if
c
c            init block of nonlocal state variables. values array always
c            allocated but with size (1,1) for std. local analyses.
c            just makes passing args simpler.
c
!DIR$ VECTOR ALIGNED
      if( local_work%block_has_nonlocal_solids )
     &    local_work%nonlocal_state_blk = zero ! array
c
c            for small displacement analysis, set integration point
c            rotations to identity.
c
      if( .not. local_work%geo_non_flg ) then ! set to identity
       if( gpn .eq. 1 ) then
!DIR$ VECTOR ALIGNED
         local_work%rot_blk_n1 = zero ! full array
         do igp = 1, local_work%num_int_points
!DIR$ VECTOR ALIGNED
           do i = 1, mxvl
             local_work%rot_blk_n1(i,1,igp) = one
             local_work%rot_blk_n1(i,5,igp) = one
             local_work%rot_blk_n1(i,9,igp) = one
           end do
         end do
       end if
      end if
c
c            uddt_displ - strain increment due to displacement
c                         increment
c            uddt_temps - (negative) of strain increment just due
c                         to temperature change
c            for iter > 1, do a usual nonlinear stress update.
c                          consistent tangent is in terms 1-36
c                          of history @ n+1 for the integration point
c            for iter = 0 and extrapolated, usual nonlinear update
c            for iter = 0 and extrapolate off, mm10 computes
c            sigma_n+1 = sigma_n + D_E * ( uddt - delta eps creep)
c            and puts D_E into history 1-36. D_E is linear-elastic
c            constitutive matrix.
c
!DIR$ VECTOR ALIGNED
      uddt = uddt_displ + uddt_temps
      iter_0_extrapolate_off = .false.
      if( iter .eq. 0 ) then
        iter_0_extrapolate_off = .not. extrapolated_du
        if( step .eq. 1 ) iter_0_extrapolate_off = .true.
      end if
      if( local_debug )  write(iout,9110) felem, gpn, span
      call mm10( gpn, local_work%span, local_work%ncrystals,
     &           hist_size_for_blk,
     &           local_work%elem_hist(1,1,gpn),
     &           local_work%elem_hist1(1,1,gpn),
     &           local_work, uddt, gp_temps,
     &           gp_dtemps, iout, signal_flag,
     &           local_work%block_has_nonlocal_solids,
     &           local_work%nonlocal_state_blk(1,1),
     &           nonlocal_shared_state_size,  ! value in param_def
     &           iter_0_extrapolate_off )
c
      if( local_work%block_has_nonlocal_solids )
     &    call drive_10_non_local ! finish nonlocal shared

      if( local_debug )  write(iout,9120) felem, gpn, span
c
      if( check_D ) then
        tol = 0.01d00
        do i = 1, span
          do j = 1, 6
            if( cep(i,j,j) .lt. tol ) then
               write(iout,*) ' .. fatal @ 1 in drive_10_update'
               call die_abort
            end if
          end do ! on j
          do j = 1, 6
            do k = 1, 6
             if( abs( cep(i,j,k) - cep(i,k,j) )
     &           .gt. 1.0d-8 ) then
               write(iout,*) ' .. fatal @ 2 drive_10_update'
               call die_abort
             end if
            end do ! on k
          end do ! on j
        end do  ! on i
      end if ! on cehck_D

      if( local_debug ) then
       write(iout,*) ".... linear elastic [D] for",
     &                 " CP in drive_10_update "
       do i = 1, span
          write(iout,*) '    ... element: ', felem+i-1
            write(iout,9100) cep(i,1:6,1:6)
       end do
      end if
c
      return

c
 9000 format(1x,'.... debug mm10. felem, gpn, span: ',i7,i3,i3)
 9010 format(10x,'... dtime, type, order, nnode:     ',e14.6,3i5,
     &     /,10x,'... step, iter, now_blk:           ',3i5,
     &     /,10x,'... temperatures, temperatures_ref: ',
     &               2l2,
     &     /,10x,'... hist_size_for_blk: ',i4 )
 9100 format(10x,6e14.5)
 9110 format(10x,'... call mm10 for felem, gpn, span: ',3i10)
 9120 format(10x,'... returned from mm10 for felem, gpn, span: ',3i10)
      contains
c     ========

c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_10_non_local                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/5/2016 rhd               *
c     *                                                              *
c     *     support routine for mm10 material driver.                *
c     *     should be inlined                                        *
c     *                                                              *
c     ****************************************************************

      subroutine drive_10_non_local
      implicit none

      integer :: n, i, elem_num
      double precision :: real_npts
c
      n = nonlocal_shared_state_size ! for convenience from param_def
      if( local_debug ) write(iout,9010) n
c
      if( gpn .eq. 1 ) then  ! zero global values for elements
        do i = 1, span
          elem_num = felem + i - 1
          if( nonlocal_flags(elem_num) )
     &         nonlocal_data_n1(elem_num)%state_values(1:n) = zero
        end do
      end if
c
      do i = 1, span ! add in this gpn nonlocal values
       elem_num = felem + i - 1
       if( nonlocal_flags(elem_num) )
     &       nonlocal_data_n1(elem_num)%state_values(1:n) =
     &       nonlocal_data_n1(elem_num)%state_values(1:n) +
     &       local_work%nonlocal_state_blk(i,1:n)
      end do
c
      if( gpn .eq. local_work%num_int_points ) then
         real_npts = dble( local_work%num_int_points )
         do i = 1, span
          elem_num = felem + i - 1
          if( nonlocal_flags(elem_num) )
     &      nonlocal_data_n1(elem_num)%state_values(1:n) =
     &      nonlocal_data_n1(elem_num)%state_values(1:n) / real_npts
         end do
c         write(iout,*) '.. drive_mm10. avg nonlocal. blk: ',now_blk
c         do i = 1, span
c         write(iout,9100) felem+i-1,
c     &      nonlocal_data_n1(elem_num)%state_values(1:n)
c         end do
      end if
c
 9010 format(/,'      processing nonlocal values. # values: ',i2 )
 9100 format(10x,i10,5e14.6)
      return
      end subroutine  drive_10_non_local



c     ****************************************************************
c     *                                                              *
c     *   ==> no longer called:  subroutine drive_10_update_b                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/9/2016 rhd               *
c     *                                                              *
c     *     support routine for mm10 material driver.                *
c     *     should be inlined                                        *
c     *                                                              *
c     ****************************************************************

      subroutine drive_10_update_b
      use crystal_data, only : c_array, angle_input, crystal_input,
     &                         data_offset
      implicit none
c
      integer :: i, elnum, ci, osn, cnum, ati, aci, tc, a, b
      double precision, allocatable :: cp_stiff(:,:,:,:),
     &                                 cp_g_rot(:,:,:,:)
      double precision :: angles(3), totalC(6,6), Cci(6,6), Srot(6,6),
     &                    Ct(6,6), local_rmat(3,3),
     &                    trans_local_rmat(3,3), trans_Srot(6,6)
      integer, allocatable :: ncrystals(:)
      character :: aconv*5, atype*7
c
      allocate( cp_stiff(mxvl,6,6,max_crystals) )
      allocate( cp_g_rot(mxvl,3,3,max_crystals) )
      allocate( ncrystals(mxvl) )
c
       do i = 1, span
          ncrystals(i) = imatprp(101,matnum)
          elnum = felem+i-1
          do ci = 1, ncrystals(i)
            if( imatprp(104,matnum) .eq. 1 ) then
                  cnum = imatprp(105,matnum)
            elseif( imatprp(104,matnum) .eq. 2 ) then
                  osn = data_offset(elnum)
                  cnum = crystal_input(osn,ci)
                  if( (cnum .gt. max_crystals) .or.
     &                  (cnum .lt. 0) ) then
                   write (iout,9501) cnum
                   call die_gracefully
                  end if
            else
                  write(iout,9502)
                  call die_gracefully
            end if
            cp_stiff(i,1:6,1:6,ci) = c_array(cnum)%elast_stiff
c
            if( imatprp(107,matnum) .eq. 1 ) then
                  angles(1) = dmatprp(108,matnum)
                  angles(2) = dmatprp(109,matnum)
                  angles(3) = dmatprp(110,matnum)
            elseif( imatprp(107,matnum) .eq. 2 ) then
                  osn = data_offset(elnum)
                  angles(1:3) = angle_input(osn,ci,1:3)
            else
                  write(iout,9502)
                  call die_gracefully
            end if
            aci = imatprp(102,matnum)
            ati = imatprp(103,matnum)
c              use helper to get the crystal -> reference rotation
            if( ati .eq. 1 ) then
                  atype = "degrees"
            elseif( ati .eq. 2) then
                  atype = "radians"
            else
                  write(iout,9503)
                  call die_gracefully
            end if
c
            if( aci .eq. 1 ) then
                  aconv = "kocks"
            else
                  write(iout,9504)
                  call die_gracefully
            end if
            call mm10_rotation_matrix( angles, aconv, atype,
     &                                 local_rmat, iout )
            cp_g_rot(i,1:3,1:3,ci) = local_rmat
          end do ! over ncrystals
      end do   !   over span
c
      cep = zero ! this is local to drive_10_update_b
c
       do i = 1, span
        tc = 0
        totalC = zero
        do ci = 1, ncrystals(i)
            local_rmat = cp_g_rot(i,1:3,1:3,ci)
            trans_local_rmat = transpose( local_rmat )
            Ct = cp_stiff(i,1:6,1:6,ci)
            Srot = zero
            call mm10_RT2RVE( trans_local_rmat, Srot)
            trans_Srot = transpose( Srot )
            Cci = matmul( Ct, trans_Srot )
            Cci = matmul( Srot, Cci )
            totalC = totalC + Cci
            tc = tc + 1
         end do
c
         totalC = totalC / dble(tc) ! average over all crystals
c
         cep(i,1:3,1:3) = totalC(1:3,1:3)
         cep(i,4,4) = totalC(4,4)
         cep(i,5,5) = totalC(5,5)
         cep(i,6,6) = totalC(6,6)
c
      end do  ! over span
c
      deallocate( cp_stiff, cp_g_rot, ncrystals )
      return
c
9501  format(/1x,
     &'>>>>> FATAL ERROR: detected in drive_10_update_b',
     & /,16x,'invalid crystal number detected: ',i3,
     & /,16x,'job aborted' )
9502  format(/1x,
     &'>>>>> FATAL ERROR: detected in drive_10_update_b',
     & /,16x,'invalid/inconsistent crystal input',
     & /,16x,'job aborted' )
9503  format(/1x,
     &'>>>>> FATAL ERROR: detected in drive_10_update_b',
     & /,16x,'unexpected/unknown angle measure',
     & /,16x,'job aborted' )
9504  format(/1x,
     &'>>>>> FATAL ERROR: detected in drive_10_update_b',
     & /,16x,'unexpected/unknown angle convention',
     & /,16x,'job aborted' )
c
      end subroutine drive_10_update_b
      end subroutine drive_10_update

c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_10_update_c                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/8/2015 rhd              *
c     *                                                              *
c     *     support routine for mm10 material driver.                *
c     *     should be inlined                                        *
c     *                                                              *
c     ****************************************************************

      subroutine drive_10_update_c( source, nrows, row, dest, nterms )
      implicit none
c
      integer :: nrows, row, nterms
      double precision ::
     &  source(nrows,nterms,nterms), dest(nterms)
c
      integer :: i, j, k
c
      k = 0
c
      do i = 1, 6
        do j = 1, 6
          k = k + 1
          dest(k) = source(row,j,i)
        end do
      end do
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine drive_10_update_a                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/8/2015 rhd              *
c     *                                                              *
c     *     support routine for mm10 material driver.                *
c     *     should be inlined                                        *
c     *                                                              *
c     ****************************************************************

      subroutine drive_10_update_a( span, mxvl, uddt,
     &                              local_cep, stress_n, stress_np1,
     &                              killed_status )
      implicit none
c
      integer :: span, mxvl
      logical :: killed_status(*)
      double precision ::
     &  local_cep(mxvl,6,6), stress_n(mxvl,6), stress_np1(mxvl,6),
     &  uddt(mxvl,6), zero
      data zero / 0.0d00 /
c
      integer i, k, m
c
c              for each element in block, update stresses by
c              [D-elastic] * uddt. uddt contains thermal increment +
c              increment from imposed nodal displacements.
c
!DIR$ VECTOR ALIGNED
      stress_np1 = stress_n
c
      do k = 1, 6
       do m = 1, 6
!DIR$ VECTOR ALIGNED
         do i = 1, span
           stress_np1(i,k) = stress_np1(i,k) +
     &                       local_cep(i,m,k) * uddt(i,m)
         end do
       end do
      end do
c
      do i = 1, span
        if( killed_status(i) ) stress_np1(i,1:6) = zero
      end do
c
      return
      end
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
      include 'param_def'
c
c                      parameter declarations
c
      real    props(mxelpr,*)
      logical lprops(mxelpr,*)
      integer iprops(mxelpr,*)
      double precision
     &  uddt(mxvl,nstr)
      include 'include_sig_up'
c
c
c                       locally defined variables
c
      double precision
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
     &                    gp_temps, gp_rtemps, local_work%alpha_vec,
     &                    type )
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
      use global_data ! old common.main
      implicit integer (a-z)
c
c                      local data
c
      integer :: info_vector(10)
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
         write(out,9000) 1
         call die_gracefully
      end if
c
      local_element_no = element_no
      if( element_no .le. 0 ) then
        if( block_no .gt. nelblk .or. block_no .le. 0 ) then
          write(out,9000) 2
          call die_gracefully
        else
          local_element_no = elblks(1,block_no)
        end if
      end if
c
      if( local_element_no .le. 0 .or.
     &    local_element_no .gt. noelem ) then
             write(out,9000) 3
             call die_gracefully
      end if
c
      mat_type = iprops(25,local_element_no)
c
c              See if we're actually a interface-damaged model
c
      if( iprops(42, local_element_no) .ne. -1 ) then
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
        write(out,9000) 4
        call die_gracefully
      end select
c
c              change history length if we are actually an
c              interface damaged material
c
      if( is_inter_dmg )
     &   call mm11_set_sizes_special(inter_mat,info_vector,
     &                               local_element_no)
c
      if( info_type .gt. 0 .and. info_type .le. 4 ) then
         value = info_vector(info_type)
      else
         write(out,9000) 5
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
c     *                   last modified : 9/28/2015 rhd              *
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
      double precision ::
     & deps(mxvl,*), strain_np1(mxvl,*)
c
!DIR$ VECTOR ALIGNED
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
c     *                   last modified : 09/28/2015 rhd             *
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
      include 'param_def'
c
c                      parameter declarations
c
      double precision ::
     &  gp_coords(mxvl,3), node_coords(mxvl,*)
c
c                     locally defined arrays-variables
c
      double precision ::
     &  sf(mxndel), xi, eta, zeta, weight, zero
      logical :: local_debug
      data zero, local_debug / 0.0d00, .false. /
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
!DIR$ VECTOR ALIGNED
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
!DIR$ VECTOR ALIGNED
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
      if( .not. local_debug ) return
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
c     *             subroutine rstgp1_make_symmetric_store           *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 10/27/20015 rhd            *
c     *                                                              *
c     *                   support for umat processing                *
c     *                                                              *
c     ****************************************************************
c
      subroutine rstgp1_make_symmetric_store( matrix, symm_vector )
      implicit none

      double precision :: matrix(6,6), symm_vector(21)
c
      double precision :: tp(6,6), symm_version(6,6), half
      integer :: i, j, k, map(6)
      data half / 0.5d00 /
      data map / 1,2,3,4,6,5 /
c
c         1. compute transpose of 6 x 6 matrix
c         2. compute symmetrized version
c         3. swap rows, cols 5 & 6 to make shear ordering
c            compatible with WARP3D
c         4. store 21 terms in lower triangle by row
c
      tp = transpose( matrix )
c
      do j = 1, 6
        do i = 1, 6
          symm_version(map(i),map(j)) = half * ( matrix(i,j) +
     &                                  tp(i,j) )
        end do
      end do
c
c      k = 1
c      do i = 1, 6
c        do j = 1, i
c          symm_vector(k) = symm_version(i,j)
c          k = k + 1
c        end do
c      end do

      symm_vector(1)  = symm_version(1,1)
      symm_vector(2)  = symm_version(2,1)
      symm_vector(3)  = symm_version(2,2)
      symm_vector(4)  = symm_version(3,1)
      symm_vector(5)  = symm_version(3,2)
      symm_vector(6)  = symm_version(3,3)
      symm_vector(7)  = symm_version(4,1)
      symm_vector(8)  = symm_version(4,2)
      symm_vector(9)  = symm_version(4,3)
      symm_vector(10) = symm_version(4,4)
      symm_vector(11) = symm_version(5,1)
      symm_vector(12) = symm_version(5,2)
      symm_vector(13) = symm_version(5,3)
      symm_vector(14) = symm_version(5,4)
      symm_vector(15) = symm_version(5,5)
      symm_vector(16) = symm_version(6,1)
      symm_vector(17) = symm_version(6,2)
      symm_vector(18) = symm_version(6,3)
      symm_vector(19) = symm_version(6,4)
      symm_vector(20) = symm_version(6,5)
      symm_vector(21) = symm_version(6,6)
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *    supporting routines for rstgp1 (to be inlined)            *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 10/18/15                       *
c     *                                                              *
c     *  support routines for rstgp1. include here so they can be    *
c     *  inlined.                                                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rstgp1_a( ndof, nnode, span, ue, due, uenh, uen1,
     &                     mxvl )
      integer :: span
      double precision ::
     &  ue(mxvl,*), due(mxvl,*), uenh(mxvl,*), uen1(mxvl,*),
     &  half
      data half / 0.5d00 /
c
      do j = 1, ndof*nnode
!DIR$ VECTOR ALIGNED
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
     &                     dfn1, itype, is_bar_elem, is_link_elem,
     &                     bar_vols )
      implicit none
c
      integer :: span, itype, i
      logical :: is_bar_elem, is_link_elem
      double precision ::
     &  internal_energy, plastic_work, gp_energies(*),
     &  det_j(*), dfn1(*), gp_plast_work(*), bar_vols(*)
c
      if( is_bar_elem ) then
!DIR$ VECTOR ALIGNED
        do i = 1, span
          internal_energy = internal_energy + gp_energies(i) *
     &                      bar_vols(i)
          plastic_work    = plastic_work +
     &                     gp_plast_work(i) * bar_vols(i)
        end do
        return
      end if
c
      if( is_link_elem ) then
!DIR$ VECTOR ALIGNED
        do i = 1, span    ! no delta-pls work. link is linear
          internal_energy = internal_energy + gp_energies(i)
        end do
        return
      end if

      if( itype .eq. 1 ) then   ! large strain
!DIR$ VECTOR ALIGNED
        do i = 1, span
          internal_energy = internal_energy + gp_energies(i) *
     &                      dfn1(i) * det_j(i)
          plastic_work    = plastic_work +
     &                      gp_plast_work(i) * dfn1(i) * det_j(i)
        end do
        return
      end if
c
      if( itype .eq. 2 ) then   ! small strain
!DIR$ VECTOR ALIGNED
        do i = 1, span
         internal_energy = internal_energy + gp_energies(i) *
     &                     det_j(i)
         plastic_work    = plastic_work + gp_plast_work(i) * det_j(i)
        end do
        return
      end if
c
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
      integer :: mxvl, span, type, nrow_ceps_blk, gpn
      double precision ::
     &  ceps_blk(nrow_ceps_blk,span,*), deps_blk(mxvl,*),
     &  stress_n(mxvl,6), stress_np1(mxvl,6)
c
c                      locals

      integer :: ielem, i, j, k
      double precision ::
     & full_cep(mxvl,6,6), zero
      data zero  / 0.0d00 /
c
c              handle solid elements (type = 1) and cohesive elements
c              (type = 2 ) to let compiler optimize loops.
c
      if( type .eq. 2 ) go to 1000
c
c              expand compressed (symmetric) [Dts] to full 6x6
c              for simplicity in coding next loop.
c
      k = 1
      do i = 1, 6
       do j = 1, i
!DIR$ VECTOR ALIGNED
         do ielem = 1, span
          full_cep(ielem,i,j) = ceps_blk(k,ielem,gpn)
          full_cep(ielem,j,i) = full_cep(ielem,i,j)
         end do
        k = k + 1
       end do
      end do
c
c              compute stress @ n+1 = stress @ n + [Dt]* deps for
c              each element in block at this integration point
c
!DIR$ VECTOR ALIGNED
      stress_np1(1:mxvl,1:6) = stress_n(1:mxvl,1:6)
c
      do i = 1, 6
       do k = 1, 6
!DIR$ VECTOR ALIGNED
         do ielem = 1, span
           stress_np1(ielem,i) = stress_np1(ielem,i) +
     &         full_cep(ielem,i,k) * deps_blk(ielem,k)
         end do
       end do
      end do
c
      return
c
c              cohesive elements have 3x3 [Dt]
c
 1000 continue
      k = 1
      do i = 1, 3
       do j = 1, i
!DIR$ VECTOR ALIGNED
         do ielem = 1, span
          full_cep(ielem,i,j) = ceps_blk(k,ielem,gpn)
          full_cep(ielem,j,i) = full_cep(ielem,i,j)
         end do
        k = k + 1
       end do
      end do
c
c              compute stress @ n+1 = stress @ n + [Dt]* deps for each
c              element in block at this integration point
c
      stress_np1(1:mxvl,1:3) = stress_n(1:mxvl,1:3)
c
      do i = 1, 3
       do k = 1, 3
!DIR$ VECTOR ALIGNED
         do ielem = 1, span
           stress_np1(ielem,i) =  stress_np1(ielem,i) +
     &         full_cep(ielem,i,k) * deps_blk(ielem,k)
         end do
       end do
      end do
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine rstgp1_store_cep                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 10/25/2015 rhd             *
c     *                                                              *
c     *  store full cep from mm.. routine into symmetric global data *
c     *  structure for elements in block for gauss point gpn         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rstgp1_store_cep( span, mxvl, gpn, gbl_ceps_blk,
     &                             nrow_ceps_blk, local_cep )
      implicit none
c
      integer :: span, mxvl, gpn, nrow_ceps_blk
      double precision ::
     &  gbl_ceps_blk(nrow_ceps_blk,span,*), local_cep(mxvl,6,6)
c
      integer i, k, ii, jj
c
      if( nrow_ceps_blk .eq. 21 ) then ! symmetric [D] 6x6
!DIR$ VECTOR ALIGNED
        do i = 1, span
          gbl_ceps_blk(1,i,gpn)  = local_cep(i,1,1)
          gbl_ceps_blk(2,i,gpn)  = local_cep(i,2,1)
          gbl_ceps_blk(3,i,gpn)  = local_cep(i,2,2)
          gbl_ceps_blk(4,i,gpn)  = local_cep(i,3,1)
          gbl_ceps_blk(5,i,gpn)  = local_cep(i,3,2)
          gbl_ceps_blk(6,i,gpn)  = local_cep(i,3,3)
          gbl_ceps_blk(7,i,gpn)  = local_cep(i,4,1)
          gbl_ceps_blk(8,i,gpn)  = local_cep(i,4,2)
          gbl_ceps_blk(9,i,gpn)  = local_cep(i,4,3)
          gbl_ceps_blk(10,i,gpn) = local_cep(i,4,4)
       end do
!DIR$ VECTOR ALIGNED
       do i = 1, span
          gbl_ceps_blk(11,i,gpn) = local_cep(i,5,1)
          gbl_ceps_blk(12,i,gpn) = local_cep(i,5,2)
          gbl_ceps_blk(13,i,gpn) = local_cep(i,5,3)
          gbl_ceps_blk(14,i,gpn) = local_cep(i,5,4)
          gbl_ceps_blk(15,i,gpn) = local_cep(i,5,5)
          gbl_ceps_blk(16,i,gpn) = local_cep(i,6,1)
          gbl_ceps_blk(17,i,gpn) = local_cep(i,6,2)
          gbl_ceps_blk(18,i,gpn) = local_cep(i,6,3)
          gbl_ceps_blk(19,i,gpn) = local_cep(i,6,4)
          gbl_ceps_blk(20,i,gpn) = local_cep(i,6,5)
          gbl_ceps_blk(21,i,gpn) = local_cep(i,6,6)
        end do
      elseif( nrow_ceps_blk .eq. 6 ) then ! symmetric [D] 3x3
!DIR$ VECTOR ALIGNED
        do i = 1, span
          gbl_ceps_blk(1,i,gpn)  = local_cep(i,1,1)
          gbl_ceps_blk(2,i,gpn)  = local_cep(i,2,1)
          gbl_ceps_blk(3,i,gpn)  = local_cep(i,2,2)
          gbl_ceps_blk(4,i,gpn)  = local_cep(i,3,1)
          gbl_ceps_blk(5,i,gpn)  = local_cep(i,3,2)
          gbl_ceps_blk(6,i,gpn)  = local_cep(i,3,3)
        end do
      elseif( nrow_ceps_blk .eq. 36 ) then ! non-symmetric [D] 6x6
          k = 1
          do i = 1, span
            do ii = 1, 6
!DIR$ VECTOR ALIGNED
              do jj = 1, 6
                gbl_ceps_blk(k,i,gpn)  = local_cep(i,ii,jj)
                k = k +1
              end do
            end do
          end do
      else
       write(*,9000)
       call die_abort
      end if
c
      return
 9000 format(/,3x,">>> FATAL ERROR: wrong. size. rstgp1_store_cep",
     &       /,3x,"                 job aborted" )
      end
