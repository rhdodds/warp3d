c     ****************************************************************
c     *                                                              *
c     *                      subroutine rknstr                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     drive updating of strains/stresses for a block of        *
c     *     elements                                                 *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rknstr( props, lprops, iprops, local_work )
      use segmental_curves
      use main_data, only : matprp, lmtprp, imatprp, dmatprp,
     &                      initial_stresses
      use mm10_defs, only : indexes_common, index_crys_hist
c
      implicit none
      include 'param_def'
c
      real    :: props(mxelpr,mxvl)   ! all 3 the same. read only
      logical :: lprops(mxelpr,mxvl)  ! props(1,1) -> 1st element of blk
      integer :: iprops(mxelpr,mxvl)
      include 'include_sig_up'
c
c                    locals
c
      integer :: span, felem, elem_type, order, ngp, nnode, step,
     &           iter, mat_type, totdof, cohesive_type, iout,
     &           i, gpn, ndof, cohes_type, ielem

      double precision, parameter :: zero = 0.0d0
      logical :: geonl, bbar, adaptive_flag, segmental, cohesive_elem,
     &           fgm_enode_props, compute_shape
      logical, parameter :: local_debug = .false.
c
c           pull values from the local block definition
c
      span            = local_work%span
      felem           = local_work%felem
      elem_type       = local_work%elem_type
      order           = local_work%int_order
      ngp             = local_work%num_int_points
      nnode           = local_work%num_enodes
      ndof            = local_work%num_enode_dof
      geonl           = local_work%geo_non_flg
      step            = local_work%step
      iter            = local_work%iter
      bbar            = local_work%bbar_flg
      mat_type        = local_work%mat_type
      totdof          = local_work%totdof
      adaptive_flag   = local_work%adaptive_flag
      cohes_type      = local_work%cohes_type
      cohesive_elem   = local_work%is_cohes_elem
      fgm_enode_props = local_work%fgm_enode_props
      compute_shape   = cohesive_elem .or. fgm_enode_props
      iout            = local_work%iout
c
c            set up to compute element volumes for bbar option and for
c            [F] bar options
c
      if( bbar .and. elem_type .eq. 2 )
     &  call rknstr_zero_vol( local_work%vol_block,
     &                        local_work%volume_block, span, mxvl )
      if( local_work%compute_f_bar ) then
!DIR$ VECTOR ALIGNED
         local_work%volume_block_0(1:span)  = zero
!DIR$ VECTOR ALIGNED
         local_work%volume_block_n(1:span)  = zero
!DIR$ VECTOR ALIGNED
         local_work%volume_block_n1(1:span) = zero
      end if
c
c           for interface/cohesive elements the global coordinates
c           and global displacements are rotated to a coordinate
c           system in which the normal axis (Z rotated) is
c           perpendicular to the surface to the cohesive element. for a
c           block of geometric nonlinear elements additional processing
c           maybe needed for reference surfaces other than the
c           mid-surface.
c
c           note: for small displacements, ce_mid, ce_n1 = ce_0
c
       if( cohesive_elem ) then
           call cohes_rot_mat( span, felem, nnode, elem_type,
     &                         local_work%ce_n1,
     &                         local_work%cohes_rot_block )
           if( geonl )
     &       call cohes_mirror_refsurf( span, mxvl, totdof, nnode,
     &                                  local_work%ce_mid )
       end if
c
c           for geometrically linear elements, compute all the shape
c           function derivatives and inverse coordinate jacobians (use
c           undeformed coordinates).
c           compute volume terms for b-bar option and finish up b-bar
c           when done with integration point loop.
c
      if( geonl ) then
         if( local_work%compute_f_bar ) then
            call rknstr_geonl_f_bar
         else
            call rknstr_geonl
         end if
      else
         call rknstr_sm_displ
      end if
c
c           all done with loop over gauss points for geometrically
c           linear and nonlinear options.
c
c           average temperature data and fgm nodal property values for
c           linear elements (or cohesive elements).
c
c           for linear elements, the reference temp, change in temp for
c           the step and temp at end of step for element nodes are
c           averaged so that the element undergoes uniform change over
c           the volume (prevents thermal shear and volumetric locking).
c
c           cohesive elements are treated as having constant
c           temperature. for later convenience, build element values
c           for cohesive temps.
c
      call rknstr_fix_temps
c
c           *** set up material data for block ***
c
c           we need vectors of material properties for elements in the
c           block. makes calling models easier and makes vectorization
c           possible for simple models. note that some properties can
c           be different for each element in the block (e.g., Young's
c           modulus). other data must be the same and we only need to
c           look at the first element in the block to set all element
c           values (e.g., if one element in the block uses a segmental
c           sig-eps curve, then all elements in block must use the same
c           curve).
c
c           for new models, all elements in the block must have the
c           same property values.
c
c           if global initial stresses exist, load array for block
c           as required based on material model.
c           (original allocate set zero)
c
      call rknstr_set_up_materials
c
      call rknstr_initial_stresses
      if( allocated( initial_stresses ) ) then
!DIR$ VECTOR ALIGNED
         do i = 1, span
           ielem = felem + i - 1
           local_work%initial_stresses(1:6,i) =
     &                     initial_stresses(1:6,ielem)
c           if( ielem .eq. 2) then
c             write(*,*) '... loading initial stresses elem 2 '
c             write(*,*) local_work%initial_stresses(1:6,i)
c           end if
         end do
      end if

c
c           compute the updated strains and stresses for all elements
c           in the block. outer loop is over integration points, inner
c           loop at lower levels is over elements in block.
c
c
      do gpn = 1, ngp
        local_work%gpn = gpn
        if ( geonl ) call rstgp1( props, lprops, iprops,
     &                            local_work )
        if ( .not. geonl ) call rstgp2( props, lprops, iprops,
     &                                  local_work )
        if ( local_work%material_cut_step ) return
      end do

c
c           For CP model, calculate the gradient of the elastic
c           rotations at the element level by linear curve fit.
c
c           For linear models this will just be based on the plastic rotations
c           which may or may not be a realistic assumption
c
      if( local_work%mat_type .eq. 10 ) call rknstr_finish_cp
c
      return
c
      contains
c     ========
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine rknstr_initial_stresses         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rknstr_initial_stresses
      implicit none
c
      integer :: i, k, ielem
      logical :: stresses_exist, value_found
      double precision, parameter :: zero = 0.0d0, tol = 1.0d-08
c
      stresses_exist = allocated( initial_stresses )
      if( .not. stresses_exist ) return
c
      value_found = .false.
      do k = 1, 6
!DIR$ VECTOR ALIGNED
         do i = 1, span
           ielem = felem + i - 1
           if( abs( initial_stresses(k,ielem) ) > tol ) then
               value_found = .true.
               go to 100
           end if
         end do
      end do
 100  continue
      if( .not. value_found ) return

      select case ( mat_type )
c
        case ( 1, 3, 5, 6, 7, 8, 10 ) ! bilinear, mises, cyclic,
!                                       creep, H_2, umat, CP
!DIR$ VECTOR ALIGNED
         do i = 1, span
           ielem = felem + i - 1
           local_work%initial_stresses(1:6,i) =
     &                     initial_stresses(1:6,ielem)
c           if( ielem .eq. 2) then
c             write(*,*) '... loading initial stresses elem 2 '
c             write(*,*) local_work%initial_stresses(1:6,i)
c           end if
         end do
        case( 2 )
           write(iout,9000) 'deformation'
           call die_gracefully
        case( 4 ) ! cohesive
           write(iout,9000) 'cohesive'
           call die_gracefully
        case default
         write(iout,9020)
         call die_gracefully
      end select
c
      return
c
 9000 format(//,'>>> Error: initial stresses not supported',
     &  ' material model: ',a,
     &  /,      '           job terminated....', //)
 9020 format(//,'>>> Fatal Error: routine rknstr_initial_stresses',
     &  /,      '                 job terminated....', //)
c
      end subroutine rknstr_initial_stresses
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine rknstr_finish_cp                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     make calls to the specific material model for block      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rknstr_finish_cp
      implicit none
c
      integer :: sh, eh, i, iblkrow, rs, re
c
c              calculate the gradient of the elastic
c              rotations at the element level by linear curve fit.
c
      sh  = indexes_common(2,1) ! first index of grad_fe
      eh  = indexes_common(2,2) ! last index of grad_fe
c
      do i = 1, local_work%span
        iblkrow = i
        if( local_work%ncrystals(i) .gt. 1 )  then
             local_work%elem_hist1(i,sh:eh,1:ngp) = zero
        else
         rs = index_crys_hist(1,3,1) ! first index of Rp, 1st crystal
         re = index_crys_hist(1,3,2) ! last index of Rp, 1st crystal
         call mm10_calc_grads( ngp, elem_type, order, geonl,
     &       local_work%rot_blk_n1(1,1,1),
     &       local_work%jac(1,1,1),
     &       local_work%elem_hist(1,rs,1),
     &       local_work%elem_hist1(1,sh,1), iout, iblkrow, span,
     &       local_work%hist_size_for_blk )
        end if
      end do
c
      return
      end subroutine rknstr_finish_cp



c     ****************************************************************
c     *                                                              *
c     *                   subroutine rknstr_set_up_materials         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     make calls to the specific material model for block      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rknstr_set_up_materials
      implicit none
c
      logical :: adaptive
c
      adaptive = adaptive_flag .and. step .gt. 1
c
      select case ( mat_type )
c
        case ( 1 )
          call setup_mm01_rknstr( span, props, lprops, iprops,
     &                            local_work )
          call setup_segmental(  span, props, lprops, iprops,
     &                           local_work )
c
        case ( 2 )
          call setup_mm02_rknstr( span, props, lprops, iprops,
     &                            local_work )
c
        case ( 3 )
          call setup_mm03_rknstr( span, props, lprops, iprops,
     &                            adaptive, local_work )
          call setup_segmental(  span, props, lprops, iprops,
     &                           local_work )
c
        case ( 4 )
          call mm04_init( iout, span, felem, props, lprops, iprops,
     &                    local_work%cohes_type,
     &                    local_work%intf_prp_block,
     &                    matprp(1,local_work%matnum),
     &                    local_work%cohes_rot_block )
c
        case ( 5 )
          call setup_mm05_rknstr( span, props, lprops, iprops,
     &                            adaptive, local_work )
          call setup_segmental(  span, props, lprops, iprops,
     &                           local_work )
c
        case ( 6 )
          call setup_mm06_rknstr( span, props, lprops, iprops,
     &                            adaptive, local_work )
c
        case ( 7 )
          call setup_mm07_rknstr( span, props, lprops, iprops,
     &                            adaptive, local_work )
c
        case ( 8 )
          call setup_umat_rknstr( span, props, lprops, iprops,
     &                            adaptive, local_work )
c
        case ( 10 )
          call setup_mm10_rknstr( span, props, lprops, iprops,
     &                            adaptive, local_work )
c
        case default
          write(iout,*) '>>> invalid material model number'
          write(iout,*) '    in rknstr. abort execution'
          call die_abort
          stop
c
      end select
c
c           If this material is going to have interface damage applied to it,
c           call the setup function for mm11
c
      if( local_work%is_inter_dmg ) then
        call setup_mm11_rknstr( span, props, lprops, iprops,
     &                            adaptive, local_work )
      end if
c
      return
      end subroutine rknstr_set_up_materials



c     ****************************************************************
c     *                                                              *
c     *                   subroutine rknstr_fix_temps                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     make temp for higher-order elements vary linearly        *
c     *     between nodes to stop locking from temp loading          *
c     *     intended for this routine to be inlined                  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rknstr_fix_temps
      implicit none
c
      integer :: i
      logical :: average
      double precision :: temp_ref, d_temp, temp_np1
c
      average =  local_work%linear_displ_elem .or.
     &             local_work%is_cohes_elem
      if( average ) then
         call average_nodal_temps( local_work%temperatures,
     &      local_work%temps_node_to_process,
     &      local_work%temperatures_ref, nnode,
     &      local_work%dtemps_node_blk, local_work%temps_node_blk,
     &      local_work%temps_ref_node_blk, mxvl, span, 1 )
      end if
c
      if( local_work%linear_displ_elem .and. fgm_enode_props ) then
         call average_fgm_properties( local_work%enode_mat_props,
     &                           mxndel, mxvl, mxndpr, nnode, span )
      end if
c
      if( local_work%is_cohes_elem ) then  ! all zero if no temps
!DIR$ VECTOR ALIGNED
         do i = 1, span
           temp_ref = local_work%temps_ref_node_blk(i,1)
           d_temp   = local_work%dtemps_node_blk(i,1)
           temp_np1 = local_work%temps_node_blk(i,1)
           local_work%cohes_temp_ref(i) = temp_ref
           local_work%cohes_dtemp(i)    = d_temp
           local_work%cohes_temp_n(i)   = temp_np1 - d_temp
         end do
      end if
c
      return
c
      end subroutine rknstr_fix_temps
c     ****************************************************************
c     *                                                              *
c     *                   subroutine rknstr_geonl                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     driver all computations to set up subsequent strain      *
c     *     computations for large displacement elements where F is  *
c     *     not needed                                               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rknstr_geonl
      implicit none
c
      integer :: gpn, enode
      double precision :: xi, eta, zeta

c
c           for geometrically nonlinear elements, the deformation
c           jacobians at (n + 1/2) and (n + 1) relative to n = 0 are
c           computed at a lower level. here we compute coordinate
c           jacobians and their inverses for the nodal coordinates
c           updated to the mid-step configuration. these will be used
c           to construct the linear-form for [B] at mid-step to multiply
c           into displacement increment to define the strain increment.
c           b-bar formulation is included as in the small displacement
c           model but using mid-step geometry.
c
c           for cohesive elements, the mid-increment formulation
c           to compute "strains" (displacement jumps) is not used.
c           we formulate the [B] using coordinates at n+1.
c
c           for some material models, we need to also compute [F] @ n.
c           A logical flag has been set by drive_eps_sig_internal_forces
c           for this case.
c
c           for the hex-8 elements, we need to perform the equivalent
c           of a b-bar adjustment on each [F] (but only needed at n
c           and n+1). drive_eps_sig_internal_forces set a logical flag
c           for this case.
c
c           the additional computations required for bar-[F] makes it
c           simpler to use a separate loop structure.
c
c
      do gpn = 1, ngp
       if( local_debug ) write(*,9050)  gpn, elem_type
       local_work%gpn = gpn
       call getgpts( elem_type, order, gpn, xi, eta, zeta,
     &               local_work%weights(gpn) )
       call derivs( elem_type, xi, eta, zeta, local_work%nxi(1,gpn),
     &              local_work%neta(1,gpn), local_work%nzeta(1,gpn) )
       if( local_debug .and. gpn .eq. 1 ) then
          write(iout,9000) gpn, elem_type, xi, eta, zeta
          write(iout,9005)
          do enode = 1, nnode
            write(iout,9010) enode, local_work%nxi(enode,gpn),
     &        local_work%neta(enode,gpn),local_work% nzeta(enode,gpn)
          end do
       end if
       if( cohesive_elem ) then
          call jacob1( elem_type, span, felem, gpn, local_work%jac,
     &              local_work%det_j_mid(1,gpn),
     &              local_work%gama_mid(1,1,1,gpn),
     &              local_work%cohes_rot_block,
     &              local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &              local_work%nzeta(1,gpn),
     &              local_work%ce_n1, nnode )
        else
          call jacob1( elem_type, span, felem, gpn, local_work%jac,
     &              local_work%det_j_mid(1,gpn),
     &              local_work%gama_mid(1,1,1,gpn),
     &              local_work%cohes_rot_block,
     &              local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &              local_work%nzeta(1,gpn),
     &              local_work%ce_mid, nnode )
        end if
c
       if( compute_shape )
     &     call shapef( elem_type, xi, eta, zeta,
     &                  local_work%shape(1,gpn) )
c
       if( bbar ) then
         call vol_terms( local_work%gama_mid(1,1,1,gpn),
     &                   local_work%det_j_mid(1,gpn),
     &                   local_work%vol_block, local_work%nxi(1,gpn),
     &                   local_work%neta(1,gpn),
     &                   local_work%nzeta(1,gpn),
     &                   local_work%volume_block, span, mxvl )
       end if
      end do
c
      if( bbar .and. elem_type .eq. 2 )
     &  call vol_avg( local_work%vol_block, local_work%volume_block,
     &                span, mxvl )

      return
 9000 format(5x,"... gpn, elem_type,  xi, eta, zeta: ",
     &  2i4, 3f10.4 )
 9005 format(10x,"... shape function derivatives ..." )
 9010 format(10x,i4,3f15.6)
 9050 format ( '>>> ready to calculate deformation gradient gpn,',
     &         ' etype: ',2i3)
      end subroutine rknstr_geonl

c     ****************************************************************
c     *                                                              *
c     *                   subroutine rknstr_sm_displ                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     driver all computations to set up subsequent strain      *
c     *     computations for small displacement elements             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rknstr_sm_displ
      implicit none
c
      integer :: gpn, enode
      double precision :: xi, eta, zeta
c
      do gpn = 1, ngp
        local_work%gpn = gpn
        call getgpts( elem_type, order, gpn, xi, eta, zeta,
     &                local_work%weights(gpn) )
        call derivs( elem_type, xi, eta, zeta, local_work%nxi(1,gpn),
     &               local_work%neta(1,gpn),local_work% nzeta(1,gpn) )
        if( local_debug .and. gpn .eq. 1 ) then
          write(iout,9000) gpn, elem_type, xi, eta, zeta
          write(iout,9005)
          do enode = 1, nnode
            write(iout,9010) enode, local_work%nxi(enode,gpn),
     &        local_work%neta(enode,gpn),local_work% nzeta(enode,gpn)
          end do
        end if
        call jacob1( elem_type, span, felem, gpn, local_work%jac,
     &    local_work%det_j(1,gpn), local_work%gama(1,1,1,gpn),
     &    local_work%cohes_rot_block,
     &    local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &    local_work%nzeta(1,gpn), local_work%ce_0, nnode )
c
        if( compute_shape )
     &     call shapef( elem_type, xi, eta, zeta,
     &                  local_work%shape(1,gpn) )
c
        if( bbar ) then
          call vol_terms( local_work%gama(1,1,1,gpn),
     &                    local_work%det_j(1,gpn),
     &                    local_work%vol_block,
     &                    local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &                    local_work%nzeta(1,gpn),
     &                    local_work%volume_block, span, mxvl )
        end if
      end do
c
      if( bbar .and. elem_type .eq. 2 )
     &  call vol_avg( local_work%vol_block, local_work%volume_block,
     &                span, mxvl )
      return
 9000 format(5x,"... gpn, elem_type,  xi, eta, zeta: ",
     &  2i4, 3f10.4 )
 9005 format(10x,"... shape function derivatives ..." )
 9010 format(10x,i4,3f15.6)
 9050 format ( '>>> ready to calculate deformation gradient gpn,',
     &         ' etype: ',2i3)

      end subroutine rknstr_sm_displ

c     ****************************************************************
c     *                                                              *
c     *                   subroutine rknstr_geonl_f_bar              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     driver all computations to set up subsequent strain      *
c     *     computations for large displacement elements             *
c     *     where F-bar is needed                                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rknstr_geonl_f_bar
      implicit none
c
      integer :: gpn, enode, i
      double precision :: xi, eta, zeta
c
c           loop for geometric nonlinear when we also need to compute
c           terms needed for [F] bar at same time. mainly we're just
c           computing volume of deformed element at  n=0, n, n+1.
c           regular coordinate jacobian is based on n+1/2 element
c           deformed shape.
c
      do gpn = 1, ngp
       if( local_debug ) write(iout,9050)  gpn, elem_type
       local_work%gpn = gpn
       call getgpts( elem_type, order, gpn, xi, eta, zeta,
     &               local_work%weights(gpn) )
       call derivs( elem_type, xi, eta, zeta, local_work%nxi(1,gpn),
     &              local_work%neta(1,gpn), local_work%nzeta(1,gpn) )
       if( local_debug .and. gpn .eq. 1 ) then
          write(iout,9000) gpn, elem_type, xi, eta, zeta
          write(iout,9005)
          do enode = 1, nnode
            write(iout,9010) enode, local_work%nxi(enode,gpn),
     &        local_work%neta(enode,gpn),local_work% nzeta(enode,gpn)
          end do
       end if
       call jacob1( elem_type, span, felem, gpn, local_work%jac,
     &              local_work%det_j_mid(1,gpn),
     &              local_work%gama_mid(1,1,1,gpn),
     &              local_work%cohes_rot_block,
     &              local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &              local_work%nzeta(1,gpn),
     &              local_work%ce_mid, nnode )
       call jacob1( elem_type, span, felem, gpn, local_work%jac,
     &              local_work%det_j(1,gpn),
     &              local_work%gama(1,1,1,gpn),
     &              local_work%cohes_rot_block,
     &              local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &              local_work%nzeta(1,gpn),
     &              local_work%ce_0, nnode )
!DIR$ VECTOR ALIGNED
       do i = 1, span
         local_work%volume_block_0(i) = local_work%volume_block_0(i)
     &                + local_work%det_j(i,gpn)
       end do
       call jacob1( elem_type, span, felem, gpn, local_work%jac,
     &              local_work%det_j(1,gpn),
     &              local_work%gama(1,1,1,gpn),
     &              local_work%cohes_rot_block,
     &              local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &              local_work%nzeta(1,gpn),
     &              local_work%ce_n, nnode )
!DIR$ VECTOR ALIGNED
        do i = 1, span
           local_work%volume_block_n(i) = local_work%volume_block_n(i)
     &                + local_work%det_j(i,gpn)
        end do
        call jacob1( elem_type, span, felem, gpn, local_work%jac,
     &              local_work%det_j(1,gpn),
     &              local_work%gama(1,1,1,gpn),
     &              local_work%cohes_rot_block,
     &              local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &              local_work%nzeta(1,gpn),
     &              local_work%ce_n1, nnode )
!DIR$ VECTOR ALIGNED
        do i = 1, span
          local_work%volume_block_n1(i) = local_work%volume_block_n1(i)
     &              + local_work%det_j(i,gpn)
        end do
c
       if( compute_shape )
     &     call shapef( elem_type, xi, eta, zeta,
     &                  local_work%shape(1,gpn) )
c
       if( bbar ) then
         call vol_terms( local_work%gama_mid(1,1,1,gpn),
     &                   local_work%det_j_mid(1,gpn),
     &                   local_work%vol_block, local_work%nxi(1,gpn),
     &                   local_work%neta(1,gpn),
     &                   local_work%nzeta(1,gpn),
     &                   local_work%volume_block, span, mxvl )
       end if
      end do
c
      if( bbar .and. elem_type .eq. 2 )
     &  call vol_avg( local_work%vol_block, local_work%volume_block,
     &                span, mxvl )
c
      return
c
 9000 format(5x,"... gpn, elem_type,  xi, eta, zeta: ",
     &  2i4, 3f10.4 )
 9005 format(10x,"... shape function derivatives ..." )
 9010 format(10x,i4,3f15.6)
 9050 format ( '>>> ready to calculate deformation gradient gpn,',
     &         ' etype: ',2i3)
c
      end subroutine rknstr_geonl_f_bar

      end subroutine rknstr

c     ****************************************************************
c     *                                                              *
c     *                   subroutine setup_segmental                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/13/2010 add only        *
c     *                                                              *
c     *     set up segemetnal stress-strain curves for use in        *
c     *     stress updating                                          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine setup_segmental( span, props, lprops, iprops,
     &                            local_work )
      use segmental_curves, only: seg_curve_table,  num_seg_points,
     &                            seg_curves, max_seg_points
c
c
      implicit none
      include 'param_def'
c
c                    parameter declarations
c
      real    props(mxelpr,mxvl)
      logical lprops(mxelpr,mxvl)
      integer iprops(mxelpr,mxvl), span
      include 'include_sig_up'
c
c                    local declarations
c
      integer :: bit_flags, curve_set_number, first_curve, no_pts
      integer, intrinsic :: iand
      real :: dumr
      double precision :: dumd
c
c                  determine if the material stress-strain properties
c                  are from a set of curves defined for the material.
c                  and set-up curve definition. all elements in
c                  the block use same curve set since they all must have
c                  the same material number. load the plastic strain values
c                  common to all curves in the set into a local work space.
c                  material properties cannot be both segmental and
c                  fgm defined at model nodes.
c
      bit_flags            = iprops(24,1)
      local_work%segmental = .true.
      if( iand(bit_flags,4) .eq. 0 ) local_work%segmental = .false.
      local_work%power_law = .not.  local_work%segmental
c
      if( local_work%segmental ) then
        if( local_work%fgm_enode_props )
     &       call errmsg2( 31, dumr, '  ', dumr, dumd  )
        curve_set_number  = iprops(21,1)
        first_curve       = seg_curve_table(2,curve_set_number)
        no_pts            = num_seg_points(first_curve)
        local_work%number_points    = no_pts
        local_work%curve_set_number = curve_set_number
        local_work%eps_curve(1:no_pts) =
     &          seg_curves(1:no_pts,1,first_curve)
      end if
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                subroutine average_nodal_temps                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 04/26/01 mcw               *
c     *                                                              *
c     *     for a block of "linear" displacement field elements      *
c     *     make nodal temperatures the average value for element    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine average_nodal_temps( temperatures,
     &      temps_node_to_process, temperatures_ref, nnodel,
     &      dtemps_node_blk, temps_node_blk, temps_ref_node_blk,
     &      mxvl, span, average_case )
      implicit none
c
c                    parameter declarations
c
      integer nnodel, mxvl, span, average_case
      logical  temperatures, temps_node_to_process,
     &         temperatures_ref
c
      double precision
     &  dtemps_node_blk(mxvl,*), temps_node_blk(mxvl,*),
     &  temps_ref_node_blk(mxvl,*)
c
c                    local
c
      logical do_average
      integer i, enode
      double precision
     & sum1(span), sum2(span), sum3(span), avg1, avg2, avg3, zero,
     & fnnodel
       data zero / 0.0d0 /


c
c                    elements with linear displacement fields
c                    can shear lock when the temperature
c                    change in the elements is not uniform.
c                    here we have a block of such elements
c                    with total, incremental and reference
c                    temperatures available at the element nodes.
c                    for each set of temperatures, redefine the
c                    nodal values to be the simple average
c                    of all the nodal values.
c
      do_average = temperatures .or.  temps_node_to_process .or.
     &             temperatures_ref
      if( .not. do_average ) return
c
      fnnodel = nnodel
c
      if( average_case .eq. 1) then
!DIR$ VECTOR ALIGNED
         sum1(1:span) = zero
!DIR$ VECTOR ALIGNED
         sum2(1:span) = zero
!DIR$ VECTOR ALIGNED
         sum3(1:span) = zero
         do enode = 1, nnodel
!DIR$ VECTOR ALIGNED
            do i = 1, span
               sum1(i) = sum1(i) + dtemps_node_blk(i,enode)
               sum2(i) = sum2(i) + temps_node_blk(i,enode)
               sum3(i) = sum3(i) + temps_ref_node_blk(i,enode)
            end do
         end do
         do enode = 1, nnodel
!DIR$ VECTOR ALIGNED
            do i = 1, span
               avg1 = sum1(i) / fnnodel
               avg2 = sum2(i) / fnnodel
               avg3 = sum3(i) / fnnodel
               dtemps_node_blk(i,enode)    = avg1
               temps_node_blk(i,enode)     = avg2
               temps_ref_node_blk(i,enode) = avg3
            end do
         end do
      end if
c
      if( average_case .eq. 2) then
 !DIR$ VECTOR ALIGNED
        sum1(1:span) = zero
         do enode = 1, nnodel
!DIR$ VECTOR ALIGNED
            do i = 1, span
               sum1(i) = sum1(i) + temps_node_blk(i,enode)
            end do
         end do
         do enode = 1, nnodel
!DIR$ VECTOR ALIGNED
            do i = 1, span
               avg1 = sum1(i) / fnnodel
               temps_node_blk(i,enode) = avg1
            end do
         end do
      end if
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                subroutine average_fgm_properties             *
c     *                                                              *
c     *                       written by : mcw                       *
c     *                                                              *
c     *                   last modified : 04/26/01 mcw               *
c     *                                                              *
c     *     for a block of "linear" displacement field elements,     *
c     *     average nodal fgm properties for elements in block.      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine average_fgm_properties( enode_mat_props, mxndel,
     &      mxvl, num_props, nnode, span)
      use main_data, only: fgm_node_values_cols
      implicit none
c
c                    parameter declarations
c
      integer nnode, mxndel, mxvl, num_props, span
      double precision
     &  enode_mat_props(mxndel,mxvl,num_props)
c
c                    local
c
      integer enode, elem, prop
c
      double precision
     & sum(span), avg, zero
c
      data zero / 0.d0 /
c
c                    elements with linear displacement fields
c                    can shear lock when properties within the
c                    the elements are not uniform. here we have
c                    a block of such elements with several fgm
c                    properties available at the element nodes.
c                    for each property, redefine the nodal values
c                    to be the simple average of all the nodal values.
c
      do prop = 1, num_props
!DIR$ VECTOR ALIGNED
         sum(1:span) = zero
         do enode = 1, nnode
!DIR$ VECTOR ALIGNED
            do elem = 1, span
               sum(elem) = sum(elem) + enode_mat_props(enode,elem,prop)
            end do
         end do
         do enode = 1, nnode
!DIR$ VECTOR ALIGNED
            do elem = 1, span
               avg = sum(elem) / nnode
               enode_mat_props(enode,elem,prop) = avg
            end do
         end do
      end do
c
      return
      end


c     ****************************************************************
c     *                                                              *
c     *                   subroutine setup_mm01_rknstr               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *              last modified : 9/26/2017 rhd                   *
c     *                                                              *
c     *     set up material model #1 (bilinear mises) for stress     *
c     *     updating                                                 *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine setup_mm01_rknstr( span, props, lprops, iprops,
     &                              local_work )
      use segmental_curves
c
      implicit none
      include 'param_def'
c
c                    parameter declarations
c
      real :: props(mxelpr,mxvl)
      logical :: lprops(mxelpr,mxvl)
      integer :: iprops(mxelpr,mxvl), span
      include 'include_sig_up'
c
      integer :: i
c
!DIR$ IVDEP
      do i = 1, span
         local_work%e_vec(i)         = props(7,i)
         local_work%e_vec_n(i)       = props(7,i)
         local_work%fgm_flags(i,1)   = props(7,i)
c
         local_work%nu_vec(i)        = props(8,i)
         local_work%nu_vec_n(i)      = props(8,i)
         local_work%fgm_flags(i,2)   = props(8,i)
c
         local_work%fgm_flags(i,3)   = props(9,i)
         local_work%alpha_vec(i,1)   = props(9,i)
         local_work%alpha_vec(i,2)   = props(13,i)
         local_work%alpha_vec(i,3)   = props(34,i)
         local_work%alpha_vec(i,4)   = props(35,i)
         local_work%alpha_vec(i,5)   = props(36,i)
         local_work%alpha_vec(i,6)   = props(37,i)
c
         local_work%alpha_vec_n(i,1) = props(9,i)
         local_work%alpha_vec_n(i,2) = props(13,i)
         local_work%alpha_vec_n(i,3) = props(34,i)
         local_work%alpha_vec_n(i,4) = props(35,i)
         local_work%alpha_vec_n(i,5) = props(36,i)
         local_work%alpha_vec_n(i,6) = props(37,i)
c
         local_work%beta_vec(i)      = props(14,i)
         local_work%h_vec(i)         = props(15,i)
         local_work%fgm_flags(i,6)   = local_work%tan_e_vec(i)
         local_work%lnelas_vec(i)    = lprops(17,i)
         local_work%sigyld_vec(i)    = props(23,i)
         local_work%fgm_flags(i,7)   = props(23,i)
c
      end do
c
      return
c
      end
c     ****************************************************************
c     *                                                              *
c     *                   subroutine setup_mm02_rknstr               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     set up material model #2 (deformation plasticity)        *
c     *     for stress updating                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine setup_mm02_rknstr( span, props, lprops, iprops,
     &                              local_work )
      use segmental_curves
c
c
      implicit none
      include 'param_def'
c
c                    parameter declarations
c
      real :: props(mxelpr,mxvl)
      logical :: lprops(mxelpr,mxvl)
      integer :: iprops(mxelpr,mxvl), span
      include 'include_sig_up'
c
      integer :: i
c
!DIR$ IVDEP
      do i = 1, span
         local_work%e_vec(i)       = props(7,i)
         local_work%nu_vec(i)      = props(8,i)
         local_work%alpha_vec(i,1) = props(9,i)
         local_work%alpha_vec(i,2) = props(13,i)
         local_work%alpha_vec(i,3) = props(34,i)
         local_work%alpha_vec(i,4) = props(35,i)
         local_work%alpha_vec(i,5) = props(36,i)
         local_work%alpha_vec(i,6) = props(37,i)
         local_work%lnelas_vec(i)  = lprops(17,i)
         local_work%sigyld_vec(i)  = props(23,i)
         local_work%n_power_vec(i) = props(21,i)
c
         local_work%fgm_flags(i,1) = props(7,i)
         local_work%fgm_flags(i,2) = props(8,i)
         local_work%fgm_flags(i,3) = props(9,i)
         local_work%fgm_flags(i,7) = props(23,i)
         local_work%fgm_flags(i,8) = props(21,i)
       end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                   subroutine setup_mm03_rknstr               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     set up material model #3 (general mises and gurson)      *
c     *     for stress updating                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine setup_mm03_rknstr( span, props, lprops, iprops,
     &                              adaptive, local_work )
      use segmental_curves
c
      implicit none
      include 'param_def'
c
c                    parameter declarations
c
      real ::   props(mxelpr,mxvl)
      logical :: lprops(mxelpr,mxvl), adaptive
      integer :: iprops(mxelpr,mxvl), span
      include 'include_sig_up'
c
c                    local
c
      integer :: i, bit_flags
      integer, intrinsic :: iand
c
!DIR$ IVDEP
      do i = 1, span
           local_work%e_vec(i)       = props(7,i)
           local_work%e_vec_n(i)     = props(7,i)
           local_work%fgm_flags(i,1) = props(7,i)
c
           local_work%nu_vec(i)      = props(8,i)
           local_work%nu_vec_n(i)    = props(8,i)
           local_work%fgm_flags(i,2) = props(8,i)
c
           local_work%fgm_flags(i,3) = props(9,i)
           local_work%alpha_vec(i,1) = props(9,i)
           local_work%alpha_vec(i,2) = props(13,i)
           local_work%alpha_vec(i,3) = props(34,i)
           local_work%alpha_vec(i,4) = props(35,i)
           local_work%alpha_vec(i,5) = props(36,i)
           local_work%alpha_vec(i,6) = props(37,i)
c
           local_work%alpha_vec_n(i,1) = props(9,i)
           local_work%alpha_vec_n(i,2) = props(13,i)
           local_work%alpha_vec_n(i,3) = props(34,i)
           local_work%alpha_vec_n(i,4) = props(35,i)
           local_work%alpha_vec_n(i,5) = props(36,i)
           local_work%alpha_vec_n(i,6) = props(37,i)
c
           local_work%f0_vec(i)      = props(26,i)
c
           local_work%h_vec(i)       = props(15,i)
           local_work%fgm_flags(i,6) = local_work%tan_e_vec(i)
           local_work%sigyld_vec(i)  = props(23,i)
           local_work%fgm_flags(i,7) = props(23,i)
           local_work%n_power_vec(i) = props(21,i)
           local_work%fgm_flags(i,8) = props(21,i)
c
           local_work%eps_ref_vec(i) = props(22,i)
           local_work%m_power_vec(i) = props(20,i)
           local_work%q1_vec(i)      = props(27,i)
           local_work%q2_vec(i)      = props(28,i)
           local_work%q3_vec(i)      = props(29,i)
           local_work%nuc_vec(i)     = .true.
           local_work%nuc_s_n_vec(i) = props(31,i)
           local_work%nuc_e_n_vec(i) = props(32,i)
           local_work%nuc_f_n_vec(i) = props(33,i)
c
      end do
c
      bit_flags = iprops(24,1)
c
!DIR$ IVDEP
      do i = 1, span
        if ( iand(iprops(30,i),1) .eq. 0 )
     &         local_work%nuc_vec(i) = .false.
      end do
c
c                   determine if material model can call for a
c                   reduction in the adaptive step size
c
      local_work%allow_cut = adaptive
      if ( iand(bit_flags,2) .eq. 0 ) local_work%allow_cut = .false.
c
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                   subroutine setup_mm05_rknstr               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     set up material model #5 (cyclic plasticity)             *
c     *     for stress updating: values constant across all g. pts.  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine setup_mm05_rknstr( span, props, lprops, iprops,
     &                              adaptive, local_work )
      use segmental_curves
      use main_data, only : matprp, lmtprp
c
      implicit none
      include 'param_def'
c
c                    parameter declarations
c
      real  ::  props(mxelpr,mxvl)
      logical :: lprops(mxelpr,mxvl), adaptive
      integer :: iprops(mxelpr,mxvl), span
      include 'include_sig_up'
c
c                    local
c
      integer :: matnum, i, bit_flags
      integer, intrinsic :: iand
c
c                  NOTE:  at present, all elements in the block must be
c                         same cyclic material defined by the user.
c                         this restriction can be removed by tracking
c                         thru and allowing matnum to vary by element in
c                         the block
c
c
      matnum = local_work%matnum
c                  storage layout:
c
c                              FA option      GP option
c                     55        q_u             gp_h_u
c                     56        b_u             gp_tau
c                     57        h_u             gp_beta_u
c                     58        1.0              -1.0
c                     59        gamma_u         gp_delta_u
c                     60 sig_tol
c                     61-64 <available>
c
!DIR$ IVDEP
      do i = 1, span
c
           local_work%e_vec(i)       = props(7,i)
           local_work%e_vec_n(i)     = props(7,i)
           local_work%nu_vec(i)      = props(8,i)
           local_work%nu_vec_n(i)    = props(8,i)
c
           local_work%alpha_vec(i,1) = props(9,i)
           local_work%alpha_vec(i,2) = props(13,i)
           local_work%alpha_vec(i,3) = props(34,i)
           local_work%alpha_vec(i,4) = props(35,i)
           local_work%alpha_vec(i,5) = props(36,i)
           local_work%alpha_vec(i,6) = props(37,i)
c
           local_work%alpha_vec_n(i,1) = props(9,i)
           local_work%alpha_vec_n(i,2) = props(13,i)
           local_work%alpha_vec_n(i,3) = props(34,i)
           local_work%alpha_vec_n(i,4) = props(35,i)
           local_work%alpha_vec_n(i,5) = props(36,i)
           local_work%alpha_vec_n(i,6) = props(37,i)
c
           local_work%gp_sig_0_vec(i)     = props(23,i)
           local_work%gp_sig_0_vec_n(i)   = props(23,i)
           local_work%gp_h_u_vec(i)       = matprp(55,matnum)
           local_work%gp_h_u_vec_n(i)     = matprp(55,matnum)
           local_work%gp_beta_u_vec(i)    = matprp(57,matnum)
           local_work%gp_beta_u_vec_n(i)  = matprp(57,matnum)
           local_work%gp_delta_u_vec(i)   = matprp(59,matnum)
           local_work%gp_delta_u_vec_n(i) = matprp(59,matnum)
c
           local_work%sigyld_vec(i)  = props(23,i)
c  ????           local_work%n_power_vec(i) = props(21,i)
c
c  ????           local_work%tan_e_vec(i)    = matprp(4,matnum)
           local_work%mm05_props(i,1:10) = matprp(55:64,matnum)
c
      end do
c
      bit_flags = iprops(24,1)
c
c                   determine if material model can call for a
c                   reduction in the adaptive step size
c
      local_work%allow_cut = adaptive
      if ( iand(bit_flags,2) .eq. 0 ) local_work%allow_cut = .false.
c
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                   subroutine setup_mm06_rknstr               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     set up material model #6 (creep)                         *
c     *     for stress updating: values constant across all g. pts.  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine setup_mm06_rknstr( span, props, lprops, iprops,
     &                              adaptive, local_work )
      use main_data, only : matprp, lmtprp
c
      implicit none
      include 'param_def'
c
c                    parameter declarations
c
      real    ::  props(mxelpr,mxvl)
      logical ::  lprops(mxelpr,mxvl), adaptive
      integer ::  iprops(mxelpr,mxvl), span
      include 'include_sig_up'
c
c                    local
c
      integer :: i, jout, felem, bit_flags, matnum
      logical :: local_debug
      integer, intrinsic :: iand
c
      local_debug = .false.
      matnum = local_work%matnum
c
!DIR$ IVDEP
      do i = 1, span
c
           local_work%e_vec(i)       = props(7,i)
           local_work%nu_vec(i)      = props(8,i)
           local_work%alpha_vec(i,1) = props(9,i)
           local_work%alpha_vec(i,2) = props(13,i)
           local_work%alpha_vec(i,3) = props(34,i)
           local_work%alpha_vec(i,4) = props(35,i)
           local_work%alpha_vec(i,5) = props(36,i)
           local_work%alpha_vec(i,6) = props(37,i)
           local_work%n_power_vec(i) = props(21,i)
           local_work%mm06_props(i,1) = matprp(80,matnum)
c
      end do

      if( local_debug ) then
         jout = local_work%iout
         write(jout,9000)
         write(jout,9005) span, matnum
         felem = local_work%felem - 1
         do i = 1, span
           write(jout,9010) felem+i, local_work%e_vec(i),
     &                      local_work%nu_vec(i),
     &                      local_work%n_power_vec(i),
     &                      local_work%mm06_props(i,1)
           write(jout,9015) local_work%alpha_vec(i,1:6)
         end do
      end if
c
c
c                   determine if material model can call for a
c                   reduction in the adaptive step size
c
      bit_flags = iprops(24,1)
      local_work%allow_cut = adaptive
      if ( iand(bit_flags,2) .eq. 0 ) local_work%allow_cut = .false.
c
      if( local_debug ) write(jout,9002)
      return
c
 9000 format(//,3x,'.... setup stress-strain update props for creep')
 9002 format(//,3x,'.... leaving setup props for creep')
 9005 format(//,3x,'.... span, matnum: ',2i6 )
 9010 format(8x,i8, 4e18.6 )
 9015 format(8x,8x,6e18.6 )
      end

c     ****************************************************************
c     *                                                              *
c     *                   subroutine setup_mm07_rknstr               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     set up material model #7 (mises + hydrogen)              *
c     *     for stress updating: values constant across all g. pts.  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine setup_mm07_rknstr( span, props, lprops, iprops,
     &                              adaptive, local_work )
      use segmental_curves
      use main_data, only : matprp, lmtprp
c
      implicit none
      include 'param_def'
c
c                    parameter declarations
c
      real  ::  props(mxelpr,mxvl)
      logical :: lprops(mxelpr,mxvl), adaptive
      integer :: iprops(mxelpr,mxvl), span
      include 'include_sig_up'
c
c                    local
c
      integer :: matnum, i, bit_flags
      integer, intrinsic :: iand
c
      matnum = local_work%matnum
c
!DIR$ IVDEP
      do i = 1, span
c
           local_work%e_vec(i)       = props(7,i)
           local_work%nu_vec(i)      = props(8,i)
           local_work%alpha_vec(i,1) = props(9,i)
           local_work%alpha_vec(i,2) = props(13,i)
           local_work%alpha_vec(i,3) = props(34,i)
           local_work%alpha_vec(i,4) = props(35,i)
           local_work%alpha_vec(i,5) = props(36,i)
           local_work%alpha_vec(i,6) = props(37,i)
           local_work%sigyld_vec(i)  = props(23,i)
           local_work%n_power_vec(i) = props(21,i)
c
           local_work%tan_e_vec(i)       = matprp(4,matnum)
           local_work%mm07_props(i,1:10) = matprp(70:79,matnum)
      end do
c
      bit_flags = iprops(24,1)
c
c                   determine if material model can call for a
c                   reduction in the adaptive step size
c
      local_work%allow_cut = adaptive
      if ( iand(bit_flags,2) .eq. 0 ) local_work%allow_cut = .false.
c
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                   subroutine setup_umat_rknstr               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *   set up material model #8 (Abaqus compatible UMAT)          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine setup_umat_rknstr( span, props, lprops, iprops,
     &                              adaptive, local_work )
      use segmental_curves
      use main_data, only : matprp, lmtprp, dmatprp
c
      implicit none
      include 'param_def'
c
c                    parameter declarations
c
      real ::   props(mxelpr,mxvl)
      logical :: lprops(mxelpr,mxvl), adaptive
      integer :: iprops(mxelpr,mxvl), span
      include 'include_sig_up'
c
c                    local
c
      integer :: i, matnum, bit_flags
      integer, intrinsic :: iand
c
      matnum = local_work%matnum
c
!DIR$ IVDEP
      do i = 1, span
        local_work%e_vec(i)       = props(7,i)
        local_work%nu_vec(i)      = props(8,i)
        local_work%alpha_vec(i,1) = props(9,i)
        local_work%alpha_vec(i,2) = props(13,i)
        local_work%alpha_vec(i,3) = props(34,i)
        local_work%alpha_vec(i,4) = props(35,i)
        local_work%alpha_vec(i,5) = props(36,i)
        local_work%alpha_vec(i,6) = props(37,i)
        local_work%umat_props(i,1:50) = dmatprp(151:200,matnum)
      end do
c
c                   for each element in block, compute a
c                   "characteristic" length per Abaqus specification.
c
       call characteristic_elem_length(
     &    local_work%elem_type, span,
     &    local_work%num_enodes,
     &    local_work%ce_n1(1,1),
     &    local_work%characteristic_length(1),
     &    local_work%iout  )
c
      bit_flags = iprops(24,1)
c
c                   determine if material model can call for a
c                   reduction in the adaptive step size
c
      local_work%allow_cut = adaptive
      if ( iand(bit_flags,2) .eq. 0 ) local_work%allow_cut = .false.
c
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine setup_mm10_rknstr               *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 1/9/2017 rhd               *
c     *                                                              *
c     *     set up material model #10 (crystal plasticity)           *
c     *     for stress updating: values constant across all g. pts.  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine setup_mm10_rknstr( span, props, lprops, iprops,
     &                              adaptive, local_work )
      use segmental_curves
      use main_data, only : matprp, lmtprp, imatprp, dmatprp, smatprp
      use crystal_data, only : c_array, angle_input, crystal_input,
     &                              data_offset
c
      implicit none
      include 'param_def'
c
c                    parameter declarations
c
      integer :: iprops(mxelpr,mxvl), span
      real    :: props(mxelpr,mxvl)
      logical :: lprops(mxelpr,mxvl), adaptive
      include 'include_sig_up'
c
c                    local
c
      integer :: i, out, matnum, ctotal, c, cnum, s, elnum, osn
      double precision :: angles(3), bs(3), ns(3), temp66(6,6),
     &                    rot(3,3), trans_rot(3,3),
     &                    temp_vec3(3),
     &                    A(3,3), trans_A(3,3), A_symm(3,3),
     &                    A_asymm(3,3),
     &                    Rstiff(6,6), trans_Rstiff(6,6)
      double precision, parameter :: half = 0.5d00, zero = 0.0d00
      character :: aconv*5, atype*7
c
      ctotal = 0
      out = local_work%iout
      matnum = local_work%matnum
c
      do i = 1, span
c
       call setup_mm10_rknstr_a( 1 )  ! selected props -> local_work
c
c              will need to (in the near future) extract crystal
c              and orientation information.  Also possibly change
c              to allow for variable number of crystals in each
c              element block.
c
       elnum = local_work%felem+i-1
c
       do c = 1, local_work%ncrystals(i)
c
         call setup_mm10_rknstr_a( 2 ) ! get data for this crystal
         call setup_mm10_rknstr_a( 3 ) ! put props into local_work
         call setup_mm10_rknstr_a( 4 ) ! get crystal->reference rot
c
c                     set up and rotate our orientation tensors
c
         rot       = local_work%c_props(i,c)%rotation_g
         call set_up_mm10_rknstr_c( trans_rot, rot, 3 )
         call setup_mm10_rknstr_a( 5 ) ! get crystal->reference rots
c
c                      rotate forward our stiffness tensor
c
         call mm10_RT2RVE( trans_rot, Rstiff ) ! makes Rstiff
         call set_up_mm10_rknstr_c( trans_Rstiff, Rstiff, 6)
         temp66 = matmul( local_work%c_props(i,c)%init_elast_stiff,
     &                       trans_Rstiff )
         local_work%c_props(i,c)%init_elast_stiff = matmul(
     &               Rstiff, temp66 )   ! 6x6 * 6xx
       end do ! on crystals
c
           ctotal = ctotal + local_work%ncrystals(i)
c
      end do ! on span

c
c                   determine if material model can call for a
c                   reduction in the adaptive step size
c
      local_work%allow_cut = adaptive .and. lmtprp(22,matnum)
c
c
      return

 9501 format(/,1x,'>>>> Not implemented yet in rknstr setup!'/)
 9503 format(/,1x,'>>>> System error: unexpected angle type in rknstr!',
     &            ' Aborting.'/)
 9504 format(/,1x,'>>>> System error: unexpected angle conv in rknstr!',
     &            ' Aborting.'/)

      contains
c     ========
      subroutine set_up_mm10_rknstr_b( a, b, c )
      implicit none
      double precision :: a(3), b(3,3), c(3)
c
c                     a = [b] * c
c
      a(1) = b(1,1)*c(1) + b(1,2)*c(2) + b(1,3)*c(3)
      a(2) = b(2,1)*c(1) + b(2,2)*c(2) + b(2,3)*c(3)
      a(3) = b(3,1)*c(1) + b(3,2)*c(2) + b(3,3)*c(3)
c
      return
      end subroutine set_up_mm10_rknstr_b

      subroutine set_up_mm10_rknstr_c( a, b, n )
      implicit none
      integer :: n, i, j
      double precision :: a(n,n), b(n,n)
c
c                     a = trans[ b ]
c
      if( n .eq. 3 ) then
         do i = 1, 3
!DIR$ VECTOR ALIGNED
            do j = 1, 3
               a(i,j) = b(j,i)
            end do
         end do
         return
      end if
c
       do i = 1, 6
!DIR$ VECTOR ALIGNED
         do j = 1, 6
            a(i,j) = b(j,i)
         end do
      end do
c
      return
      end subroutine set_up_mm10_rknstr_c

      subroutine set_up_mm10_rknstr_d( trans_A, A, bs, ns )
      implicit none
c
      integer :: i, j
      double precision :: trans_A(3,3), A(3,3), bs(3), ns(3)
c
      double precision :: spreadbs(3,3), spreadns(3,3)
c
      spreadbs(1,1) = bs(1)
      spreadbs(2,1) = bs(2)
      spreadbs(3,1) = bs(3)
      spreadbs(1,2) = bs(1)
      spreadbs(2,2) = bs(2)
      spreadbs(3,2) = bs(3)
      spreadbs(1,3) = bs(1)
      spreadbs(2,3) = bs(2)
      spreadbs(3,3) = bs(3)
c
      spreadns(1,1) = ns(1)
      spreadns(2,1) = ns(1)
      spreadns(3,1) = ns(1)
      spreadns(1,2) = ns(2)
      spreadns(2,2) = ns(2)
      spreadns(3,2) = ns(2)
      spreadns(1,3) = ns(3)
      spreadns(2,3) = ns(3)
      spreadns(3,3) = ns(3)
c
      do i = 1, 3
!DIR$ VECTOR ALIGNED
        do j = 1, 3
         A(i,j) = spreadbs(i,j) * spreadns(i,j)
         trans_A(j,i) = A(i,j)
        end do
      end do
c
      return
c
      end  subroutine set_up_mm10_rknstr_d




      subroutine setup_mm10_rknstr_a( dowhat )
      implicit none
c
      integer :: dowhat

      select case( dowhat )

      case( 1 )
c
        local_work%alpha_vec(i,1) = props(9,i)
        local_work%alpha_vec(i,2) = props(13,i)
        local_work%alpha_vec(i,3) = props(34,i)
        local_work%alpha_vec(i,4) = props(35,i)
        local_work%alpha_vec(i,5) = props(36,i)
        local_work%alpha_vec(i,6) = props(37,i)
c
        local_work%alpha_vec_n(i,1) = props(9,i)
        local_work%alpha_vec_n(i,2) = props(13,i)
        local_work%alpha_vec_n(i,3) = props(34,i)
        local_work%alpha_vec_n(i,4) = props(35,i)
        local_work%alpha_vec_n(i,5) = props(36,i)
        local_work%alpha_vec_n(i,6) = props(37,i)
c
        local_work%debug_flag(i) = lmtprp(13,matnum)
        local_work%local_tol(i) = dmatprp(100,matnum)
c
c              may eventually change this to allow for
c              different # of crystals in block
c
        local_work%ncrystals(i) = imatprp(101,matnum)
c
        local_work%angle_convention(i) = imatprp(102,matnum)
        local_work%angle_type(i) = imatprp(103,matnum)
c
      case( 2 )
c
c              get the local crystal number
c
        if( imatprp(104,matnum) .eq. 1 ) then
              cnum = imatprp(105,matnum)
        elseif( imatprp(104,matnum) .eq. 2 ) then
              osn = data_offset(elnum)
              cnum = crystal_input(osn,c)
c
c             couldn't do this earlier, so check here
c
              if( (cnum .gt. max_crystals) .or.
     &              (cnum .lt. 0) ) then
               write (out,'("Crystal ", i3, " not valid")')
     &              cnum
                    call die_gracefully
              end if
        else
              write(out,9502)
              call die_gracefully
        end if
c
c             get the local orientation
c
        if (imatprp(107,matnum ) .eq. 1) then
             angles(1) = dmatprp(108,matnum)
             angles(2) = dmatprp(109,matnum)
             angles(3) = dmatprp(110,matnum)
        elseif( imatprp(107,matnum) .eq. 2 ) then
             osn = data_offset(elnum)
             angles(1:3) = angle_input(osn,c,1:3)
        else
             write (out,9502)
             call die_gracefully
        end if
c
      case( 3 )
c
c             we have the properties, we just need to extract
c             into our local structure
c
        local_work%c_props(i,c)%init_elast_stiff =
     &        c_array(cnum)%elast_stiff
        local_work%c_props(i,c)%init_angles = angles
        local_work%c_props(i,c)%nslip = c_array(cnum)%nslip
        local_work%c_props(i,c)%rateN = c_array(cnum)%harden_n
        local_work%c_props(i,c)%tauHat_y = c_array(cnum)%tau_hat_y
        local_work%c_props(i,c)%Go_y = c_array(cnum)%g_o_y
        local_work%c_props(i,c)%tauHat_v = c_array(cnum)%tau_hat_v
        local_work%c_props(i,c)%Go_v = c_array(cnum)%g_o_v
        local_work%c_props(i,c)%tau_a = c_array(cnum)%tau_a
        local_work%c_props(i,c)%burgers = c_array(cnum)%b
        local_work%c_props(i,c)%p_v = c_array(cnum)%p_v
        local_work%c_props(i,c)%q_v = c_array(cnum)%q_v
        local_work%c_props(i,c)%p_y = c_array(cnum)%p_y
        local_work%c_props(i,c)%q_y = c_array(cnum)%q_y
        local_work%c_props(i,c)%boltzman = c_array(cnum)%boltz
        local_work%c_props(i,c)%theta_o = c_array(cnum)%theta_o
        local_work%c_props(i,c)%eps_dot_o_v = c_array(cnum)%eps_dot_o_v
        local_work%c_props(i,c)%eps_dot_o_y = c_array(cnum)%eps_dot_o_y
        local_work%c_props(i,c)%mu_o = c_array(cnum)%mu_o
        local_work%c_props(i,c)%D_o  = c_array(cnum)%D_o
c
        local_work%c_props(i,c)%t_o = c_array(cnum)%t_o
        local_work%c_props(i,c)%tau_a = c_array(cnum)%tau_a
        local_work%c_props(i,c)%k_o = c_array(cnum)%k_o
        local_work%c_props(i,c)%h_type = c_array(cnum)%h_type
        local_work%c_props(i,c)%s_type = c_array(cnum)%slip_type
        local_work%c_props(i,c)%cnum = cnum
        local_work%c_props(i,c)%num_hard = c_array(cnum)%num_hard
        local_work%c_props(i,c)%tang_calc = c_array(cnum)%tang_calc
        local_work%c_props(i,c)%u1 = c_array(cnum)%u1
        local_work%c_props(i,c)%u2 = c_array(cnum)%u2
        local_work%c_props(i,c)%u3 = c_array(cnum)%u3
        local_work%c_props(i,c)%u4 = c_array(cnum)%u4
        local_work%c_props(i,c)%u5 = c_array(cnum)%u5
        local_work%c_props(i,c)%u6 = c_array(cnum)%u6
        local_work%c_props(i,c)%u7 = c_array(cnum)%u7
        local_work%c_props(i,c)%u8 = c_array(cnum)%u8
        local_work%c_props(i,c)%u9 = c_array(cnum)%u9
        local_work%c_props(i,c)%u10 = c_array(cnum)%u10
        local_work%c_props(i,c)%tau_y = c_array(cnum)%tau_y
        local_work%c_props(i,c)%tau_v = c_array(cnum)%tau_v
        local_work%c_props(i,c)%voche_m = c_array(cnum)%voche_m
        local_work%c_props(i,c)%iD_v = c_array(cnum)%iD_v
c
        local_work%c_props(i,c)%cp_001 = c_array(cnum)%cp_001
        local_work%c_props(i,c)%cp_002 = c_array(cnum)%cp_002
        local_work%c_props(i,c)%cp_003 = c_array(cnum)%cp_003
        local_work%c_props(i,c)%cp_004 = c_array(cnum)%cp_004
        local_work%c_props(i,c)%cp_005 = c_array(cnum)%cp_005
        local_work%c_props(i,c)%cp_006 = c_array(cnum)%cp_006
        local_work%c_props(i,c)%cp_007 = c_array(cnum)%cp_007
        local_work%c_props(i,c)%cp_008 = c_array(cnum)%cp_008
        local_work%c_props(i,c)%cp_009 = c_array(cnum)%cp_009
        local_work%c_props(i,c)%cp_010 = c_array(cnum)%cp_010
        local_work%c_props(i,c)%cp_011 = c_array(cnum)%cp_011
        local_work%c_props(i,c)%cp_012 = c_array(cnum)%cp_012
        local_work%c_props(i,c)%cp_013 = c_array(cnum)%cp_013
        local_work%c_props(i,c)%cp_014 = c_array(cnum)%cp_014
        local_work%c_props(i,c)%cp_015 = c_array(cnum)%cp_015
        local_work%c_props(i,c)%cp_016 = c_array(cnum)%cp_016
        local_work%c_props(i,c)%cp_017 = c_array(cnum)%cp_017
        local_work%c_props(i,c)%cp_018 = c_array(cnum)%cp_018
        local_work%c_props(i,c)%cp_019 = c_array(cnum)%cp_019
        local_work%c_props(i,c)%cp_020 = c_array(cnum)%cp_020
        local_work%c_props(i,c)%cp_021 = c_array(cnum)%cp_021
        local_work%c_props(i,c)%cp_022 = c_array(cnum)%cp_022
        local_work%c_props(i,c)%cp_023 = c_array(cnum)%cp_023
        local_work%c_props(i,c)%cp_024 = c_array(cnum)%cp_024
        local_work%c_props(i,c)%cp_025 = c_array(cnum)%cp_025
        local_work%c_props(i,c)%cp_026 = c_array(cnum)%cp_026
        local_work%c_props(i,c)%cp_027 = c_array(cnum)%cp_027
        local_work%c_props(i,c)%cp_028 = c_array(cnum)%cp_028
        local_work%c_props(i,c)%cp_029 = c_array(cnum)%cp_029
        local_work%c_props(i,c)%cp_030 = c_array(cnum)%cp_030
        local_work%c_props(i,c)%cp_031 = c_array(cnum)%cp_031
        local_work%c_props(i,c)%cp_032 = c_array(cnum)%cp_032
        local_work%c_props(i,c)%cp_033 = c_array(cnum)%cp_033
        local_work%c_props(i,c)%cp_034 = c_array(cnum)%cp_034
        local_work%c_props(i,c)%cp_035 = c_array(cnum)%cp_035
        local_work%c_props(i,c)%cp_036 = c_array(cnum)%cp_036
        local_work%c_props(i,c)%cp_037 = c_array(cnum)%cp_037
        local_work%c_props(i,c)%cp_038 = c_array(cnum)%cp_038
        local_work%c_props(i,c)%cp_039 = c_array(cnum)%cp_039
        local_work%c_props(i,c)%cp_040 = c_array(cnum)%cp_040
        local_work%c_props(i,c)%cp_041 = c_array(cnum)%cp_041
        local_work%c_props(i,c)%cp_042 = c_array(cnum)%cp_042
        local_work%c_props(i,c)%cp_043 = c_array(cnum)%cp_043
        local_work%c_props(i,c)%cp_044 = c_array(cnum)%cp_044
        local_work%c_props(i,c)%cp_045 = c_array(cnum)%cp_045
        local_work%c_props(i,c)%cp_046 = c_array(cnum)%cp_046
        local_work%c_props(i,c)%cp_047 = c_array(cnum)%cp_047
        local_work%c_props(i,c)%cp_048 = c_array(cnum)%cp_048
        local_work%c_props(i,c)%cp_049 = c_array(cnum)%cp_049
        local_work%c_props(i,c)%cp_050 = c_array(cnum)%cp_050
        local_work%c_props(i,c)%cp_051 = c_array(cnum)%cp_051
        local_work%c_props(i,c)%cp_052 = c_array(cnum)%cp_052
        local_work%c_props(i,c)%cp_053 = c_array(cnum)%cp_053
        local_work%c_props(i,c)%cp_054 = c_array(cnum)%cp_054
        local_work%c_props(i,c)%cp_055 = c_array(cnum)%cp_055
        local_work%c_props(i,c)%cp_056 = c_array(cnum)%cp_056
        local_work%c_props(i,c)%cp_057 = c_array(cnum)%cp_057
        local_work%c_props(i,c)%cp_058 = c_array(cnum)%cp_058
        local_work%c_props(i,c)%cp_059 = c_array(cnum)%cp_059
        local_work%c_props(i,c)%cp_060 = c_array(cnum)%cp_060
        local_work%c_props(i,c)%cp_061 = c_array(cnum)%cp_061
        local_work%c_props(i,c)%cp_062 = c_array(cnum)%cp_062
        local_work%c_props(i,c)%cp_063 = c_array(cnum)%cp_063
        local_work%c_props(i,c)%cp_064 = c_array(cnum)%cp_064
        local_work%c_props(i,c)%cp_065 = c_array(cnum)%cp_065
        local_work%c_props(i,c)%cp_066 = c_array(cnum)%cp_066
        local_work%c_props(i,c)%cp_067 = c_array(cnum)%cp_067
        local_work%c_props(i,c)%cp_068 = c_array(cnum)%cp_068
        local_work%c_props(i,c)%cp_069 = c_array(cnum)%cp_069
        local_work%c_props(i,c)%cp_070 = c_array(cnum)%cp_070
        local_work%c_props(i,c)%cp_071 = c_array(cnum)%cp_071
        local_work%c_props(i,c)%cp_072 = c_array(cnum)%cp_072
        local_work%c_props(i,c)%cp_073 = c_array(cnum)%cp_073
        local_work%c_props(i,c)%cp_074 = c_array(cnum)%cp_074
        local_work%c_props(i,c)%cp_075 = c_array(cnum)%cp_075
        local_work%c_props(i,c)%cp_076 = c_array(cnum)%cp_076
        local_work%c_props(i,c)%cp_077 = c_array(cnum)%cp_077
        local_work%c_props(i,c)%cp_078 = c_array(cnum)%cp_078
        local_work%c_props(i,c)%cp_079 = c_array(cnum)%cp_079
        local_work%c_props(i,c)%cp_080 = c_array(cnum)%cp_080
        local_work%c_props(i,c)%cp_081 = c_array(cnum)%cp_081
        local_work%c_props(i,c)%cp_082 = c_array(cnum)%cp_082
        local_work%c_props(i,c)%cp_083 = c_array(cnum)%cp_083
        local_work%c_props(i,c)%cp_084 = c_array(cnum)%cp_084
        local_work%c_props(i,c)%cp_085 = c_array(cnum)%cp_085
        local_work%c_props(i,c)%cp_086 = c_array(cnum)%cp_086
        local_work%c_props(i,c)%cp_087 = c_array(cnum)%cp_087
        local_work%c_props(i,c)%cp_088 = c_array(cnum)%cp_088
        local_work%c_props(i,c)%cp_089 = c_array(cnum)%cp_089
        local_work%c_props(i,c)%cp_090 = c_array(cnum)%cp_090
        local_work%c_props(i,c)%cp_091 = c_array(cnum)%cp_091
        local_work%c_props(i,c)%cp_092 = c_array(cnum)%cp_092
        local_work%c_props(i,c)%cp_093 = c_array(cnum)%cp_093
        local_work%c_props(i,c)%cp_094 = c_array(cnum)%cp_094
        local_work%c_props(i,c)%cp_095 = c_array(cnum)%cp_095
        local_work%c_props(i,c)%cp_096 = c_array(cnum)%cp_096
        local_work%c_props(i,c)%cp_097 = c_array(cnum)%cp_097
        local_work%c_props(i,c)%cp_098 = c_array(cnum)%cp_098
        local_work%c_props(i,c)%cp_099 = c_array(cnum)%cp_099
        local_work%c_props(i,c)%cp_100 = c_array(cnum)%cp_100
c
c          flags to control CP internal solvers
c
        local_work%c_props(i,c)%solver = c_array(cnum)%solver
        local_work%c_props(i,c)%strategy = c_array(cnum)%strategy
        local_work%c_props(i,c)%gpall = c_array(cnum)%gpall
        local_work%c_props(i,c)%gpp = c_array(cnum)%gpp
        local_work%c_props(i,c)%st_it(1:3) = c_array(cnum)%st_it(1:3)
        local_work%c_props(i,c)%method = c_array(cnum)%method
        local_work%c_props(i,c)%miter = c_array(cnum)%miter
        local_work%c_props(i,c)%atol = c_array(cnum)%atol
        local_work%c_props(i,c)%atol1 = c_array(cnum)%atol1
        local_work%c_props(i,c)%rtol = c_array(cnum)%rtol
        local_work%c_props(i,c)%rtol1 = c_array(cnum)%rtol1
        local_work%c_props(i,c)%xtol = c_array(cnum)%xtol
        local_work%c_props(i,c)%xtol1 = c_array(cnum)%xtol1
c
c          alternative model flag
c
        local_work%c_props(i,c)%alter_mode = c_array(cnum)%alter_mode
c
      case( 4 )
c
c          reference rotation
c
        if (local_work%angle_type(i) .eq. 1) then
              atype = "degrees"
        elseif (local_work%angle_type(i) .eq. 2) then
              atype = "radians"
        else
              write(out,9503)
              call die_gracefully
        end if
        if (local_work%angle_convention(i) .eq. 1) then
              aconv="kocks"
        else
              write(out,9504)
              call die_gracefully
        end if
        call mm10_rotation_matrix( local_work%c_props(i,c)%init_angles,
     &        aconv, atype, local_work%c_props(i,c)%rotation_g, out )
c
      case( 5 )
c
c          rotation on each slip system for this crystal
c
        do s = 1, local_work%c_props(i,c)%nslip
            temp_vec3 = c_array(cnum)%bi(s,1:3)
            call set_up_mm10_rknstr_b( bs, trans_rot, temp_vec3 )
            temp_vec3 =  c_array(cnum)%ni(s,1:3)
            call set_up_mm10_rknstr_b( ns, trans_rot, temp_vec3 )
            local_work%c_props(i,c)%ns(1:3,s) = ns
            call set_up_mm10_rknstr_d( trans_A, A, bs, ns )
            A_symm  = half * ( A + trans_A )
            A_asymm = half * ( A - trans_A )
            call mm10_ET2EV( A_symm,  local_work%c_props(i,c)%ms(1,s) )
            call mm10_WT2WV( A_asymm, local_work%c_props(i,c)%qs(1,s) )
        end do ! on slip systems

      end select
c
      return
 9502 format(/,1x,'>>>> System error: unexpected input type in rknstr!',
     &            ' Aborting.'/)
 9503 format(/,1x,'>>>> System error: unexpected angle type in rknstr!',
     &            ' Aborting.'/)
 9504 format(/,1x,'>>>> System error: unexpected angle conv in rknstr!',
     &            ' Aborting.'/)

c
      end subroutine setup_mm10_rknstr_a

c
      end subroutine setup_mm10_rknstr

c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine setup_mm11_rknstr               *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 03/10/2014                 *
c     *                                                              *
c     *     set up material model #11 (crystal interface dmg)        *
c     *     for stress updating: values constant across all g. pts.  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine setup_mm11_rknstr( span, props, lprops, iprops,
     &                              adaptive, local_work )
      use segmental_curves
      use main_data, only : matprp, lmtprp, imatprp, dmatprp, smatprp
      use crystal_data, only : c_array, angle_input, crystal_input,
     &                              data_offset, simple_angles, nangles,
     &                              mc_array
c
      implicit integer (a-z)
      include 'param_def'
c
c                    parameter declarations
c
      real    props(mxelpr,mxvl)
      logical lprops(mxelpr,mxvl)
      integer iprops(mxelpr,mxvl)
      include 'include_sig_up'
c
c                    local
c
      logical adaptive
      integer :: ctotal, c, cnum, s, e
      double precision, dimension(3) :: angles
      character :: aconv*5
      character :: atype*7
      double precision, dimension(3) :: bs, ns
      double precision, dimension(3,3) :: A
      double precision, dimension(6,6) :: Rstiff, temp
      integer :: elnum, osn
c
c
      ctotal = 0
c
      matnum = local_work%inter_mat
c
      local_work%sv(1) = dmatprp(116, matnum)
      local_work%sv(2) = dmatprp(117, matnum)
      local_work%sv(3) = dmatprp(118, matnum)
c
      local_work%lv(1) = dmatprp(119, matnum)
      local_work%lv(2) = dmatprp(120, matnum)
      local_work%lv(3) = dmatprp(121, matnum)
c
      local_work%tv(1) = dmatprp(122, matnum)
      local_work%tv(2) = dmatprp(123, matnum)
      local_work%tv(3) = dmatprp(124, matnum)
c
      local_work%ls = dmatprp(126, matnum)
      local_work%ll = dmatprp(127, matnum)
      local_work%lt = dmatprp(128, matnum)
c
      local_work%alpha_dmg = dmatprp(129, matnum)
c
      do i = 1, span
           e = i + local_work%felem
c
c            Calculate the number of crystals per GP based on the initial
c            element size
c
           call mm11_elem_size(local_work%elem_type,
     &          local_work%num_enodes, local_work%ce_0(i,1:3*
     &          local_work%num_enodes),
     &          local_work%sv, local_work%lv, local_work%tv,
     &          local_work%ls, local_work%ll, local_work%lt,
     &          local_work%nstacks(i), local_work%nper(i))
c
           local_work%debug_flag(i) = lmtprp(13,matnum)
           local_work%local_tol(i) = dmatprp(100,matnum)
           local_work%ncrystals(i) = local_work%nstacks(i)*
     &            local_work%nper(i)*2
c
           local_work%angle_convention(i) = imatprp(102,matnum)
           local_work%angle_type(i) = imatprp(103,matnum)
c
c           Make sure we aren't asking for more orientations than we have
c
      if ((local_work%ncrystals(i) .gt. nangles) .or.
     &      (local_work%ncrystals(i) .gt. imatprp(101,matnum))) then
            write(*,*) "Too many angles required!"
            write(*,*) i+local_work%felem-1
            write(*,*) local_work%ncrystals(i)
            write(*,*) local_work%nstacks(i)
            write(*,*) local_work%nper(i)
            call die_abort
      end if
c
c          Extract props for each crystal
c
           do c= 1, local_work%ncrystals(i)
c
c                       Get the local crystal number
c
                  cnum = imatprp(105,matnum)
c
c                       Get the local orientation
c
                  angles(1:3) = simple_angles(mc_array(e,c),1:3)
c
c                       Now we have the properties, we just need
c                        to extract into our local structure
c
                  local_work%c_props(i,c)%init_elast_stiff =
     &                  c_array(cnum)%elast_stiff
                  local_work%c_props(i,c)%init_angles = angles
                  local_work%c_props(i,c)%nslip = c_array(cnum)%nslip
                  local_work%c_props(i,c)%rateN = c_array(cnum)%harden_n
                  local_work%c_props(i,c)%tauHat_y
     &                  = c_array(cnum)%tau_hat_y
                  local_work%c_props(i,c)%Go_y = c_array(cnum)%g_o_y
                  local_work%c_props(i,c)%tauHat_v
     &                  = c_array(cnum)%tau_hat_v
                  local_work%c_props(i,c)%Go_v = c_array(cnum)%g_o_v
                  local_work%c_props(i,c)%tau_a = c_array(cnum)%tau_a
                  local_work%c_props(i,c)%burgers = c_array(cnum)%b
                  local_work%c_props(i,c)%p_v = c_array(cnum)%p_v
                  local_work%c_props(i,c)%q_v = c_array(cnum)%q_v
                  local_work%c_props(i,c)%p_y = c_array(cnum)%p_y
                  local_work%c_props(i,c)%q_y = c_array(cnum)%q_y
                  local_work%c_props(i,c)%boltzman = c_array(cnum)%boltz
                  local_work%c_props(i,c)%theta_o =
     &                  c_array(cnum)%theta_o
                  local_work%c_props(i,c)%eps_dot_o_v =
     &                  c_array(cnum)%eps_dot_o_v
                  local_work%c_props(i,c)%eps_dot_o_y =
     &                  c_array(cnum)%eps_dot_o_y
                  local_work%c_props(i,c)%mu_o =
     &                  c_array(cnum)%mu_o
                  local_work%c_props(i,c)%D_o  =
     &                  c_array(cnum)%D_o
                  local_work%c_props(i,c)%t_o = c_array(cnum)%t_o
                  local_work%c_props(i,c)%tau_a = c_array(cnum)%tau_a
                  local_work%c_props(i,c)%k_o = c_array(cnum)%k_o
                  local_work%c_props(i,c)%h_type =
     &                  c_array(cnum)%h_type
                  local_work%c_props(i,c)%u1 = c_array(cnum)%u1
                  local_work%c_props(i,c)%u2 = c_array(cnum)%u2
                  local_work%c_props(i,c)%u3 = c_array(cnum)%u3
                  local_work%c_props(i,c)%u4 = c_array(cnum)%u4
                  local_work%c_props(i,c)%u5 = c_array(cnum)%u5
                  local_work%c_props(i,c)%u6 = c_array(cnum)%u6
                  local_work%c_props(i,c)%tau_y = c_array(cnum)%tau_y
                  local_work%c_props(i,c)%tau_v = c_array(cnum)%tau_v
                  local_work%c_props(i,c)%voche_m =
     &                  c_array(cnum)%voche_m
c
c                 Call a helper to get the crystal -> reference rotation
c
                  if (local_work%angle_type(i) .eq. 1) then
                        atype = "degrees"
                  elseif (local_work%angle_type(i) .eq. 2) then
                        atype = "radians"
                  else
                        write(out,9503)
                        call die_gracefully
                  end if

                  if (local_work%angle_convention(i) .eq. 1) then
                        aconv="kocks"
                  else
                        write(out,9504)
                        call die_gracefully
                  end if
c
                  call mm10_rotation_matrix(
     &                  local_work%c_props(i,c)%init_angles,
     &                  aconv, atype,
     &                  local_work%c_props(i,c)%rotation_g, out)
c
c                 Now that we have that, we can set up and rotate our
c                 orientation tensors
c
                  do s = 1, local_work%c_props(i,c)%nslip
                        bs = matmul(transpose(
     &                   local_work%c_props(i,c)%rotation_g),
     &                   c_array(cnum)%bi(s,:))
                        ns = matmul(transpose(
     &                   local_work%c_props(i,c)%rotation_g),
     &                   c_array(cnum)%ni(s,:))
                        local_work%c_props(i,c)%ns(1:3,s) =
     &                        ns
                        A = spread(bs,dim=2,ncopies=size(ns))*spread(
     &                        ns,dim=1,ncopies=size(bs))
                        call mm10_ET2EV(0.5*(A+transpose(A)),
     &                        local_work%c_props(i,c)%ms(:,s))
                        call mm10_WT2WV(0.5*(A-transpose(A)),
     &                        local_work%c_props(i,c)%qs(:,s))
                  end do
c
c           We can also rotate forward our stiffness tensor
c
            call mm10_RT2RVE( transpose(
     &            local_work%c_props(i,c)%rotation_g), Rstiff)
            local_work%c_props(i,c)%init_elast_stiff = matmul(
     &            Rstiff, matmul(
     &            local_work%c_props(i,c)%init_elast_stiff,
     &            transpose(Rstiff)))

           end do
c
           ctotal = ctotal + local_work%ncrystals(i)
c
      end do
c
c                   determine if material model can call for a
c                   reduction in the adaptive step size
c
      local_work%allow_cut = adaptive .and. lmtprp(22,matnum)
c
      return

 9501 format(/,1x,'>>>> Not implemented yet in rknstr setup!'/)
 9502 format(/,1x,'>>>> System error: unexpected input type in rknstr!',
     &            ' Aborting.'/)
 9503 format(/,1x,'>>>> System error: unexpected angle type in rknstr!',
     &            ' Aborting.'/)
 9504 format(/,1x,'>>>> System error: unexpected angle conv in rknstr!',
     &            ' Aborting.'/)
      end


c
c     ****************************************************************
c     *                                                              *
c     *            subroutine characteristic_elem_length             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/26/2017 rhd              *
c     *                                                              *
c     *     compute characteristic length for elements in block      *
c     *     for possible use by material update routines             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine characteristic_elem_length(
     &       etype, span, nnodel, node_coords, lengths, iout  )
      implicit none
      include 'param_def'
c
c                      parameter declarations
c
      integer :: etype, span, nnodel, iout
      double precision :: node_coords(mxvl,*), lengths(*)
c
c                     locally defined arrays-variables
c
      double precision :: rnlengths, xa, xb,
     &  ya, yb, za, zb, local_sums(mxvl), scale_factor
      double precision, parameter :: half=0.5d0, one=1.0d0, zero=0.0d0
      logical :: brick, tet, wedge, linear
      logical, parameter :: local_debug = .false.
      integer :: j, nodea, nodeb, i, node_pairs(4,2), nlengths

c
      brick = etype .ge. 1  .and. etype .le. 5
      tet   = etype .eq. 6 .or. etype .eq. 13
      wedge = etype .eq. 7
      linear = etype .eq. 2  .or.  etype .eq. 13
      if( .not. ( brick .or. tet ) ) then
          write(iout,9100)
          call die_abort
      end if
c
c                     set pairs of element nodes to use for computing
c                     a characteristic length. for hex we use the average
c                     of cross diagonals. for tets just average of edge
c                     lengths. take half of size for quadratic elements.
c
      if( brick ) then
        nlengths = 4
        rnlengths = one / real( nlengths )
        node_pairs(1,1) = 1; node_pairs(1,2) = 7
        node_pairs(2,1) = 2; node_pairs(2,2) = 8
        node_pairs(3,1) = 4; node_pairs(3,2) = 6
        node_pairs(4,1) = 3; node_pairs(4,2) = 5
      end if
      if( tet ) then
        nlengths = 4
        rnlengths =  one / real( nlengths )
        node_pairs(1,1) = 1; node_pairs(1,2) = 4
        node_pairs(2,1) = 1; node_pairs(2,2) = 2
        node_pairs(3,1) = 1; node_pairs(3,2) = 3
        node_pairs(4,1) = 2; node_pairs(4,2) = 3
      end if
c
      local_sums(1:span) = zero
c
c                     compute length for each pair of nodes for elements and
c                     sum. then take average and apply linear vs. quadratic
c                     scale factor.
c
      do j = 1, nlengths
        nodea = node_pairs(j,1)
        nodeb = node_pairs(j,2)
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
        do i = 1, span
          xa = node_coords(i,nodea)
          ya = node_coords(i,nodea+nnodel)
          za = node_coords(i,nodea+nnodel+nnodel)
          xb = node_coords(i,nodeb)
          yb = node_coords(i,nodeb+nnodel)
          zb = node_coords(i,nodeb+nnodel+nnodel)
          local_sums(i) = local_sums(i) + sqrt( (xa-xb)**2 + (ya-yb)**2
     &                                    + (za-zb)**2 )
         end do
      end do
c
      scale_factor = half
      if( linear ) scale_factor = one
c
!DIR$ IVDEP
!DIR$ VECTOR ALIGNED
      do i = 1, span
       lengths(i) = local_sums(i) * rnlengths * scale_factor
      end do

c
      if ( .not. local_debug ) return
         write(iout,*) '>> in characteristic_elem_length'
         write(iout,*) '.. element lengths:'
         write(iout,9000) ( i, lengths(i), i=1, span )
      return
c
 9000 format(1x,i4,3x,f15.6 )
 9100 format('>>>> Fatal Error: in characteristic_elem_length. invalid',
     &   /,  '                  element type. Terminate execution.' )
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine zero_vol                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 09/5/2017 rhd              *
c     *                                   meant to be inlined        *
c     *                                                              *
c     ****************************************************************
c
      subroutine rknstr_zero_vol ( vol, volume, span, mxvl )
      implicit none
      integer :: span, mxvl
      double precision  :: vol(mxvl,8,3), volume(span)
      double precision, parameter :: zero = 0.0d0
      vol = zero
      volume = zero
      return
      end

