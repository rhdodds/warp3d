c     ****************************************************************
c     *                                                              *
c     *                      subroutine rklstf                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 09/16/2015 rhd             *
c     *                                                              *
c     *     supervises the computation of linear stiffness matrices  *
c     *     for a block of similar elements.                         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rklstf( local_work, props, iprops, ek, nsz )
      use main_data, only : matprp, lmtprp, asymmetric_assembly, dmatprp
      implicit integer (a-z)
$add param_def
c
c                parameter declarations
c
$add include_lin_ek
      real props(mxelpr,*)
#dbl      double precision
#sgl      real
     & ek(nsz, *)
      dimension iprops(mxelpr,*)
c
c               locally defined
c
#dbl      double precision
#sgl      real
     & ce(mxvl,mxecor), sum_temps(mxvl), sum_rtemps(mxvl)
c
      logical local_debug, bbar, compute_shape, set_up_e_nu,
     &        cohesive, umat_solid, regular, crystal, average
c
#dbl      double precision
#sgl      real
     &   zero, one, beta_fact, eps_bbar, xi, eta, zeta, fnnode, avg
c
      data local_debug, zero, one
     &     / .false., 0.0d00, 1.0d00 /
      real :: cp_e, cp_nu
c
c          set local variables from the data structure
c          passed into routine. just saves lots of
c          long structure references
c
      felem           = local_work%felem
      elem_type       = local_work%elem_type
      int_order       = local_work%int_order
      mat_type        = local_work%mat_type
      matnum          = local_work%matnum
      nnode           = local_work%num_enodes
      ndof            = local_work%num_enode_dof
      totdof          = local_work%totdof
      bbar            = local_work%bbar_flg
      ngp             = local_work%num_int_points
      span            = local_work%span
      utsz            = local_work%utsz
      beta_fact       = local_work%beta_fact
      eps_bbar        = local_work%eps_bbar
      ce              = local_work%ce
      compute_shape   = local_work%is_cohes_elem .or.
     &                  local_work%is_axisymm_elem .or.
     &                  local_work%fgm_enode_props
      iout            = local_work%iout_local
c
c          set up material and other data needed for
c          linear stiffness computation.
c
c          we have regular material models builtin to
c          WARP3D, a  cohesive model builtin and
c          an Abaqus compatible umat.
c
c          for regular models, the material
c          is simple homogenous, isotropic but can be
c          temperature dependent. Element [K]s
c          compute using properties at end of step (n+1).
c          if modulus and poisson ratio are
c          temperature dependent, these values are
c          overwritten later with gauss point dependent
c          values. set flags for later use about
c          segmental stress-strain curves which could
c          be temperature dependent.
c
c          for cohesive elements pull their material props.
c
c          for umat, just pull user defined properties.
c          compute the characteristic element length
c          needed by the umat. the umat must take care of
c          temperature dependent material props (we give
c          gauss point temperature to umat).
c
c          local_work%fgm_flags(i,j) = -99 indicates that
c          element in block has an fgm material property
c          requiring interpolation.
c
      cohesive    = mat_type .eq. 4
      umat_solid  = mat_type .eq. 8
      crystal     = mat_type .eq. 10
      regular     = .not. ( cohesive .or. umat_solid )
      local_work%segmental = .false.
c      
      if( cohesive )
     &    call mm04_init( iout, span, felem, props, lprops, iprops,
     &                    local_work%cohes_type,
     &                    local_work%intf_prp_block,
     &                    matprp(1,matnum),
     &                    local_work%cohes_rot_block  )
c
      if( umat_solid ) then    ! Abaqus umat
@!DIR$ LOOP COUNT MAX=###      
         do i = 1, span
           local_work%umat_props(i,1:50) = dmatprp(151:200,matnum)
         end do
         call characteristic_elem_length(
     &      elem_type, span, nnode,
     &      local_work%ce(1,1),
     &      local_work%characteristic_length(1),
     &      local_work%iout_local  )
      end if
c
c          need to set the average properties differently for cp
c
      if( regular ) then
c
         if ( crystal ) then         
@!DIR$ LOOP COUNT MAX=###      
           do i = 1, span
              call mm10_set_e_nu(matnum,felem+i-1,cp_e,cp_nu)
              local_work%e_block(i)     = cp_e
              local_work%nu_block(i)    = cp_nu
              local_work%fgm_flags(i,1) = props(7,i)
              local_work%fgm_flags(i,2) = props(8,i)
           end do
         else
@!DIR$ LOOP COUNT MAX=###      
           do i = 1, span
              local_work%e_block(i)     = props(7,i)
              local_work%nu_block(i)    = props(8,i)
              local_work%fgm_flags(i,1) = props(7,i)
              local_work%fgm_flags(i,2) = props(8,i)
           end do
         end if
c
         bit_flags                 = iprops(24,1)
         local_work%segmental      = .true.
         if( iand(bit_flags,4) .eq. 0 ) local_work%segmental = .false.
c         
      end if
c
c           zero element stiffnesses in local work block.
c           we compute locally then copy to global.
c
      ntermsk = utsz*span
      if( asymmetric_assembly ) ntermsk = totdof*totdof*span
      call rklstf_zero_vector( ek, ntermsk )
c
c                       compute all the shape function derivatives, and
c                       inverse jacobians. also calculate volume
c                       terms if using bbar. ce is element
c                       nodal coordinates at time = 0.
c
      if( bbar .and. (elem_type .eq. 2) )
     &   call zero_vol( local_work%vol_block, local_work%volume_block,
     &                  span, mxvl )
c
c           for cohesive-interface element block,
c           element coordinates are rotated to a coordinate
c           system in which the normal axis (Z rotated) is
c           perpendicular to the surface ot the cohesive element
c
      if( local_work%is_cohes_elem )
     &     call cohes_rot_mat( span, felem, nnode, elem_type, ce,
     &                         local_work%cohes_rot_block )
c
      do gpn = 1, ngp
         call getgpts( elem_type, int_order, gpn, xi, eta, zeta,
     &                 local_work%weights(gpn) )
         call derivs( elem_type, xi, eta, zeta, local_work%nxi(1,gpn),
     &             local_work%neta(1,gpn), local_work%nzeta(1,gpn) )
         if( compute_shape )
     &     call shapef( elem_type, xi, eta, zeta,
     &                  local_work%shape(1,gpn) )
        if( local_debug .and. gpn .eq. 1 ) then
          write(iout,9000) gpn, elem_type, xi, eta, zeta
          write(iout,9005)
          do enode = 1, nnode
            write(iout,9010) enode, local_work%nxi(enode,gpn),
     &        local_work%neta(enode,gpn),local_work% nzeta(enode,gpn)
          end do
         end if
         call jacob1( elem_type,span, felem, gpn, local_work%jac_block,
     &                local_work%det_jac_block(1,gpn),
     &                local_work%gama(1,1,1,gpn),
     &                local_work%cohes_rot_block,
     &                local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &                local_work%nzeta(1,gpn), ce, nnode )
c
        if( bbar ) then
          call vol_terms(
     &    local_work%gama(1,1,1,gpn), local_work%det_jac_block(1,gpn),
     &    local_work%vol_block, local_work%nxi(1,gpn),
     &    local_work%neta(1,gpn), local_work%nzeta(1,gpn),
     &    local_work%volume_block, span, mxvl )
        end if
      end do
c
      if( bbar ) call vol_avg( local_work%vol_block,
     &                          local_work%volume_block, span, mxvl )
c
c          for linear displacement elements and
c          cohesive-interface elements, set the average
c          (single) value of temperature for
c          element nodes at n+1. same for reference temp
c          this prevents shear and volumetric locking under
c          temperature gradients over elements.
c
      average = ( local_work%linear_displ_elem .or.
     &            local_work%is_cohes_elem ) .and.
     &            local_work%temps_node_to_process
      if( average )  call rklstf_avg_temps( span, mxvl, nnode,
     &                local_work%temps_node_blk(1,1),
     &                local_work%temps_ref_node_blk(1,1) )
c
c          average fgm nodal property values of
c          e and nu for linear displacment elements.
c          (prevents shear and volumetric locking).
c
      if( local_work%linear_displ_elem .and.
     &    local_work%fgm_enode_props ) then
         call average_fgm_properties( local_work%enode_mat_props,
     &                                mxndel, mxvl, 2, nnode, span)
      end if
c
c
c          compute linear stiffnesses, one
c          gauss point at a time. zero local
c          [b] matrices. lower level routines
c          assume zero. also the elastic [d]
c          matrices. evaluate e and nu at
c          temperature for gauss point if properties
c          are temperature dependent.
c
c          for functionally graded materials, we interpolate
c          some material props at the gauss points from
c          model node values. at present only e and nu for
c          solid elements and volume fraction for
c          ductile component of fgm cohesive are
c          interpolated.
c
      nterms =  nstr * mxedof * mxvl
      call rklstf_zero_vector( local_work%b_block, nterms )
      call rklstf_zero_vector( local_work%bt_block, nterms )
      call rklstf_zero_vector( local_work%bd_block, nterms )
      call rklstf_zero_vector( local_work%cep, nstr*nstr*mxvl )
c
      set_up_e_nu = local_work%segmental .and.
     &              (mat_type .eq. 1) .or. (mat_type .eq. 3) .or.
     &              (mat_type .eq. 5) .or. (mat_type .eq. 6)
c
      do  gpn = 1, ngp
c
       if ( set_up_e_nu )
     &       call set_e_nu_for_block(
     &          span, local_work%nu_block, local_work%e_block,
     &          local_work%segmental, felem, elem_type,
     &          int_order, gpn, nnode,
     &          local_work%temps_node_to_process,
     &          local_work%temps_node_blk )
c
       if ( local_work%fgm_enode_props ) then
             call set_fgm_solid_props_for_block(
     &                span, felem, elem_type, gpn, nnode,
     &                local_work%e_block, local_work%shape(1,gpn),
     &                local_work%enode_mat_props, 1,
     &                local_work%fgm_flags(1,1) )
             call set_fgm_solid_props_for_block(
     &                span, felem, elem_type, gpn, nnode,
     &                local_work%nu_block, local_work%shape(1,gpn),
     &                local_work%enode_mat_props, 2,
     &                local_work%fgm_flags(1,2) )
       end if
c
       igpn = gpn
       if( .not. asymmetric_assembly ) then
         call gplns1(
     &    igpn, local_work%icp,
     &    local_work%e_block, local_work%nu_block,
     &    local_work%intf_prp_block, ek,
     &    local_work%gama(1,1,1,gpn), local_work%det_jac_block(1,gpn),
     &    local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &    local_work%nzeta(1,gpn), local_work%shape(1,gpn),
     &    local_work%ce, local_work%weights(gpn),
     &    local_work%vol_block, local_work%b_block,
     &    local_work%bt_block, local_work%bd_block, local_work%cep,
     &    local_work%cohes_rot_block,
     &    utsz, local_work, utsz)
       else
         call gplns1(
     &    igpn, local_work%icp,
     &    local_work%e_block, local_work%nu_block,
     &    local_work%intf_prp_block, ek,
     &    local_work%gama(1,1,1,gpn), local_work%det_jac_block(1,gpn),
     &    local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &    local_work%nzeta(1,gpn), local_work%shape(1,gpn),
     &    local_work%ce, local_work%weights(gpn),
     &    local_work%vol_block, local_work%b_block,
     &    local_work%bt_block, local_work%bd_block, local_work%cep,
     &    local_work%cohes_rot_block,
     &    utsz, local_work, totdof*totdof)
       end if
c
      end do  ! on gpn 
c
c
c		       modify element stiffness matrix by thickness
c         factor for a plane strain analysis
c
      if( beta_fact .ne. one ) then
            if( .not. asymmetric_assembly ) then
              ek(1:utsz,1:span) = beta_fact * ek(1:utsz,1:span)
            else
              ek(1:totdof*totdof,1:span) = beta_fact * 
     &            ek(1:totdof*totdof,1:span)
            end if
      end if
c
c          transform element stiffnesses to constraint
c          compatible coordinates if required. skip
c          whole block or specfic element to reduce work.
c
      if( local_work%trn_e_block ) then
@!DIR$ LOOP COUNT MAX=###      
        do i = 1, span
         if( local_work%trn_e_flags(i) ) then
           if( .not. asymmetric_assembly ) then
             call trnmtx( ek, local_work%cp,
     &                 local_work%trnmte, local_work%trne, ndof, nnode,
     &                 totdof, i, utsz )
           else
             call trnmtx( ek, local_work%cp,
     &                 local_work%trnmte, local_work%trne, ndof, nnode,
     &                 totdof, i, totdof*totdof ) 
           end if
         end if
        end do ! on i over span
      end if
c
      return
c
 9000 format(5x,"... gpn, elem_type,  xi, eta, zeta: ",
     &  2i4, 3f10.4 )
 9005 format(10x,"... shape function derivatives ..." )
 9010 format(10x,i4,3f15.6)
 9100 format(8x,/,'>> entered rklstf...' )
 9200 format(12x,//,'>> upper-triangular stiffnesses for elements',
     &   ' in block...')
 9210 format(3x,i8,1000(6e15.6,/11x))
 9300 format(1x,////,
     & '   ELEMENT TESTING, values in rklstf.f',//,
     & '           span = ',i6,'  (number of elements in block)',/,
     & '          nnode = ',i6,/,
     & '         totdof = ',i6,/,
     & '      elem_type = ',i6,/,
     & '          Felem = ',i6,'  (first element in block)',/,
     & '           utsz = ',i6,'  (size of ek storage vector)',/,
     & '   num Gauss pt = ',i6,/)
 9400 format(1x,//,'>> Upper triangluar matrix:')
 9405 format(1x,'   element #',i6)
 9410 format(1x,'j=',i3,2x,'row=',i3,2x,'col=',i3,2x,'ek=',f25.8)
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine rklstf_zero_vector           *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 04/8/12                    *
c     *                                                              *
c     *     zero a vector of specified length w/ floating zero       *
c     *                                                              *
c     ****************************************************************
c
      subroutine rklstf_zero_vector( vec, n )
      implicit none
c      
#dbl      double precision ::
#sgl      real ::
     &  vec(*)
      integer :: n
c
      vec(1:n) = 0.0d00
c
      return
      end


c     ****************************************************************
c     *                                                              *
c     *              subroutine rklstf_avg_temps                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/7/13 rhd                *
c     *                                                              *
c     *   for linear-displacement elements, replace nodal temps      *
c     *   with average for element - prevents shear-volumetric       *
c     *   locking for thermal gradients                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine rklstf_avg_temps( span, mxvl, nnode, temps_node_blk,
     &                             temps_ref_node_blk )
      implicit none
c
c               parameters
c
      integer  span, mxvl, nnode
#dbl      double precision
#sgl      real
     &  temps_node_blk(mxvl,*), temps_ref_node_blk(mxvl,*)
c
c               locally defined
c
      integer enode, i
#dbl      double precision
#sgl      real
     & sum_temps(mxvl), sum_rtemps(mxvl)
c
#dbl      double precision
#sgl      real
     &   zero, fnnode
c
      data  zero /  0.0d00 /
c
c
@!DIR$ LOOP COUNT MAX=###      
      sum_temps(1:span)  = zero
@!DIR$ LOOP COUNT MAX=###      
      sum_rtemps(1:span) = zero
c
      do enode = 1, nnode
@!DIR$ LOOP COUNT MAX=###      
       do i = 1, span
         sum_temps(i)  = sum_temps(i) + temps_node_blk(i,enode)
         sum_rtemps(i) = sum_rtemps(i) + temps_ref_node_blk(i,enode)
       end do
      end do
c
      fnnode =  real( nnode )
      do enode = 1, nnode
@!DIR$ LOOP COUNT MAX=###      
        do i = 1, span
          temps_node_blk(i,enode) = sum_temps(i) / fnnode
          temps_ref_node_blk(i,enode) = sum_rtemps(i) / fnnode
        end do
      end do
c
      return
      end





