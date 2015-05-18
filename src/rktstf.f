c     ****************************************************************
c     *                                                              *
c     *                      subroutine rktstf                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 12/09/13 mcm               *
c     *                                                              *
c     *     drive computation of tangent stiffness matrices for a    *
c     *     block of similar elements.                               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rktstf( props, iprops, lprops, ek, nrowek, ispan,
     &                   local_work )
      use main_data, only : matprp, lmtprp, asymmetric_assembly
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                 parameter declarations
c
#dbl      double precision
#sgl      real
     &   ek(nrowek,ispan)
      real    props(mxelpr,*)
      integer iprops(mxelpr,*)
      logical lprops(mxelpr,*)
c
c                 local
c
      logical geonl, local_debug, bbar, first, qbar_flag,
     &        compute_shape, cohes_mirror, dummy_logic,
     &        average
c
#dbl      double precision
#sgl      real
     &  xi, eta, zeta, beta_fact, eps_bbar, zero, one, dummy,
     &  temp_ref, d_temp, temp_np1
c
      data local_debug, dummy_logic / .false., .true. /
      data zero, one, dummy / 0.d0, 1.d0, 1.d0 /
c
      if ( local_debug ) then
         call iodevn( innum, out, dummy, 1 )
         write(out,9100)
      end if
c
c              set common properties for all elements in block
c
      span      = local_work%span
      felem     = local_work%felem
      utsz      = local_work%utsz
      beta_fact = local_work%beta_fact
      eps_bbar  = local_work%eps_bbar
      first     = local_work%first
      iter      = local_work%iter
      qbar_flag = local_work%qbar_flag
c
      type      = local_work%elem_type
      order     = local_work%int_order
      ngp       = local_work%num_int_points
      nnode     = local_work%num_enodes
      ndof      = local_work%num_enode_dof
      geonl     = local_work%geo_non_flg
      totdof    = nnode * ndof
      bbar      = local_work%bbar_flg
      mat_type  = local_work%mat_type
      surf       = local_work%surface
      compute_shape =  local_work%is_cohes_elem .or.
     &                 local_work%is_axisymm_elem .or.
     &                 local_work%fgm_enode_props
c
      if (.not. asymmetric_assembly) then
         call rktstf_zero_vec( ek, utsz*span )
      else
         call rktstf_zero_vec( ek, totdof*totdof*span)
      end if
c
c               compute all the shape function derivatives,
c               inverse coordinate jacobians and determinants.
c               get weight function values for gauss integration.
c               for bbar elements (8-node bricks) we also
c               the additional volume terms. for geonl,
c               we use the nodal coordinates updated to
c               n+1 to compute jacobians, bbar terms,etc.
c
      if( bbar .and. type .eq. 2 ) then
         call zero_vol( local_work%vol_block,
     &                  local_work%volume_block, span, mxvl )
      end if
c
c             for cohesive elements the element coordinates and
c             the global displacements are rotated to a coordinate
c             system in which the normal axis (Z rotated) is
c             perpendicular to the referernce surface of the
c             cohesive element
c
      if( local_work%is_cohes_elem )  then
            call cohes_rot_mat( span, felem, nnode, type,
     &                          local_work%ce,
     &                          local_work%cohes_rot_block )
          if ( geonl )
     &      call cohes_mirror_refsurf( span, mxvl, totdof, nnode,
     &                                 local_work%ce )
      end if
c
      do gpn = 1, ngp
         call getgpts( type, order, gpn, xi, eta, zeta,
     &                 local_work%weights(gpn) )
         call derivs( type, xi, eta, zeta, local_work%nxi(1,gpn),
     &                local_work%neta(1,gpn), local_work%nzeta(1,gpn) )
         if( compute_shape )
     &       call shapef( type, xi, eta, zeta, local_work%shape(1,gpn) )
         call jacob1( type, span, felem, gpn, local_work%jac_block,
     &                local_work%det_jac_block(1,gpn),
     &                local_work%gama_block(1,1,1,gpn),
     &                local_work%cohes_rot_block,
     &                local_work%nxi(1,gpn), local_work%neta(1,gpn),
     &                local_work%nzeta(1,gpn), local_work%ce,
     &                nnode )
         if ( bbar ) then
           call vol_terms( local_work%gama_block(1,1,1,gpn),
     &                     local_work%det_jac_block(1,gpn),
     &                     local_work%vol_block, local_work%nxi(1,gpn),
     &                     local_work%neta(1,gpn),
     &                     local_work%nzeta(1,gpn),
     &                     local_work%volume_block,
     &                     span, mxvl )
         end if
      end do
c
      if( bbar ) call vol_avg( local_work%vol_block,
     &                          local_work%volume_block, span, mxvl )
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
       average =  local_work%linear_displ_elem .or.
     &            local_work%is_cohes_elem
      if( average ) then
         call average_nodal_temps( local_work%temperatures,
     &      local_work%temps_node_to_process,
     &      local_work%temperatures_ref, nnode,
     &      local_work%dtemps_node_blk, local_work%temps_node_blk,
     &      local_work%temps_ref_node_blk, mxvl, span, 1 )
      end if
c
c               average fgm nodal property values of e and nu for
c               linear displacement elements to prevent shear and
c               volumetric locking. e and nu are the first 2 fgm
c               property values.
c
      if( local_work%linear_displ_elem .and.
     &    local_work%fgm_enode_props ) then
         call average_fgm_properties( local_work%enode_mat_props,
     &                                mxndel, mxvl, 2, nnode, span )
      end if
c
      if( local_work%is_cohes_elem ) then  ! all zero if no temps
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
c               generate tangent stiffnesses for elements
c               in block, one gauss point at a time. zero
c               [b] matrices used in computations for block -
c               lower routines fill in only non-zero terms.
c
      nlength = mxvl * mxedof * nstr
      call rktstf_zero_vec( local_work%b_block,  nlength )
      call rktstf_zero_vec( local_work%bd_block, nlength )
      call rktstf_zero_vec( local_work%bt_block, nlength )
c
c               for material models, we define block size
c               vectors of key material properies. otherwise,
c               the models could not be vectorized.
c
      select case ( mat_type )
         case( 1 )
            call rktstf_set_01( props, iprops, lprops, span,
     &                          local_work )
            bit_flags            = iprops(24,1)
            local_work%segmental = iand(bit_flags,4) .ne. 0
         case( 2 )
            call rktstf_set_02( props, iprops, lprops, span,
     &                          local_work )
            local_work%segmental = .false.
         case( 3 )
            call rktstf_set_03( props, iprops, lprops, span,
     &                          local_work )
            bit_flags            = iprops(24,1)
            local_work%segmental = iand(bit_flags,4) .ne. 0
         case ( 4 )
            call mm04_init( local_work%iout, span, felem, props,
     &                      lprops,
     &                      iprops, local_work%cohes_type,
     &                      local_work%intf_prp_block,
     &                      matprp(1,local_work%matnum) )
         case( 5 )
            call rktstf_set_05( props, iprops, lprops, span,
     &                          local_work )
            bit_flags            = iprops(24,1)
            local_work%segmental =  iand(bit_flags,4) .ne. 0
         case( 6 )
            call rktstf_set_06( props, iprops, lprops, span,
     &                          local_work )
            local_work%segmental = .false.
         case( 7 )
            call rktstf_set_07( props, iprops, lprops, span,
     &                          local_work )
            local_work%segmental = .false.
         case( 8 )
            call rktstf_set_08( props, iprops, lprops, span,
     &                          local_work )
            local_work%segmental = .false.
         case (10)
            call rktstf_set_10( props, iprops, lprops, span,
     &                          local_work )
            local_work%segmental = .false.
         case default
          write(*,9500)
          call die_abort
      end select
c
c               loop to process all gauss points. the gptns1 routine
c               computes stiffness contribution for all elements in
c               the block at a gauss point.
c
      do gpn = 1, ngp
        if ( .not. asymmetric_assembly) then
        call gptns1( local_work%cp, local_work%icp, gpn, props,
     &               iprops, ek, local_work, utsz )
        else
        call gptns1( local_work%cp, local_work%icp, gpn, props,
     &               iprops, ek, local_work, totdof*totdof )
        end if
      end do
c
c		             modify element stiffness matrix by thickness factor for
c             		plane strain analysis
c
      if ( beta_fact .ne. one ) ek = beta_fact * ek
c
c               transform element stiffnesses to constraint
c               compatible coordinates if required. skip
c               whole block or specfic element to reduce work.
c
      if ( local_work%trn_e_block ) then
        do i = 1, span
           if ( local_work%trn_e_flags(i) ) then
             if (.not. asymmetric_assembly) then
               call trnmtx( ek, local_work%cp,
     &                      local_work%trnmte, local_work%trne, ndof,
     &                      nnode, totdof, i, utsz )
             else
               call trnmtx( ek, local_work%cp,
     &                      local_work%trnmte, local_work%trne, ndof,
     &                      nnode, totdof, i, totdof*totdof )
             end if
           end if
         end do
      end if
c
      go to 9999
c
 9999 continue
c
 9100 format(8x,'>> entered rktstf...' )
 9200 format(12x,'>> upper-triangular stiffnesses for elements',
     &   ' in block...')
 9210 format(3x,i8,1000(6e15.6,/11x))
 9400 format(1x,//,'>> Upper triangluar matrix:')
 9405 format(1x,'   element #',i6)
 9410 format(1x,'j=',i3,2x,'row=',i3,2x,'col=',i3,2x,'ek=',f25.8)
 9300 format(1x,////,
     & '  ELEMENT TESTING, values in rktstf.f',//,
     & '       span = ',i6,/,
     & '      nnode = ',i6,/,
     & '     totdof = ',i6,/,
     & '  elem_type = ',i6,/,
     & '      felem = ',i6,/,
     & '       utsz = ',i6,/)
 9500 format(1x,'>> Fatal Error: rktstf. invalid material type..',
     &    /, 1x,'                job terminated' )
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine rktstf_set_01                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/28/02 rhd               *
c     *                                                              *
c     *          set up material model #1 (bilinear mises) for       *
c     *          tangent stiffness computations                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rktstf_set_01( props, iprops, lprops, span,
     &                          local_work )
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                 parameter declarations
c
      real    props(mxelpr,*)
      integer iprops(mxelpr,*)
      logical lprops(mxelpr,*)
c
      do i = 1, span
        local_work%e_v(i)         = props(7,i)
        local_work%fgm_flags(i,1) = props(7,i)
        local_work%nu_v(i)        = props(8,i)
        local_work%fgm_flags(i,2) = props(8,i)
        local_work%beta_v(i)      = props(14,i)
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine rktstf_set_02                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/28/02 rhd               *
c     *                                                              *
c     *      set up material model #2 (deformation plasticity)       *
c     *      for tangent stiffness computations                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rktstf_set_02( props, iprops, lprops, span,
     &                          local_work )
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                 parameter declarations
c
      real    props(mxelpr,*)
      integer iprops(mxelpr,*)
      logical lprops(mxelpr,*)
c
      do i = 1, span
        local_work%nu_v(i)        = props(8,i)
        local_work%e_v(i)         = props(7,i)
        local_work%sigyld_v(i)    = props(23,i)
        local_work%n_power_v(i)   = props(21,i)
        local_work%fgm_flags(i,1) = props(7,i)
        local_work%fgm_flags(i,2) = props(8,i)
        local_work%fgm_flags(i,7) = props(23,i)
        local_work%fgm_flags(i,8) = props(21,i)
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                   subroutine rktstf_set_03                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/28/02 rhd               *
c     *                                                              *
c     *        set up material model #3 (general mises and gurson)   *
c     *        bilinear mises) for tangent stiffness computations    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rktstf_set_03( props, iprops, lprops, span,
     &                          local_work )
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                 parameter declarations
c
      real    props(mxelpr,*)
      integer iprops(mxelpr,*)
      logical lprops(mxelpr,*)
c
      do i = 1, span
         local_work%e_v(i)         = props(7,i)
         local_work%fgm_flags(i,1) = props(7,i)
         local_work%nu_v(i)        = props(8,i)
         local_work%fgm_flags(i,2) = props(8,i)
         local_work%q1_v(i)        = props(27,i)
         local_work%q2_v(i)        = props(28,i)
         local_work%q3_v(i)        = props(29,i)
         local_work%nuc_v(i)       = .true.
         local_work%nuc_s_n_v(i)   = props(31,i)
         local_work%nuc_e_n_v(i)   = props(32,i)
         local_work%nuc_f_n_v(i)   = props(33,i)
         if ( iand (iprops(30,i),1) .eq. 0 )
     &       local_work%nuc_v(i) = .false.
      end do
c
      return
      end



c     ****************************************************************
c     *                                                              *
c     *                   subroutine rktstf_set_05                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/18/02                   *
c     *                                                              *
c     *        set up material model #5 (adv. cyclic plasticity)     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rktstf_set_05( props, iprops, lprops, span,
     &                          local_work )
      use main_data, only : matprp, lmtprp
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                 parameter declarations
c

      real    props(mxelpr,*)
      integer iprops(mxelpr,*)
      logical lprops(mxelpr,*)
c
      matnum = local_work%matnum
c
      do i = 1, span
c
         local_work%e_v(i)   = props(7,i)
         local_work%nu_v(i)  = props(8,i)
         local_work%mm05_props(i,1:10) = matprp(55:64,matnum)
c
      end do
c
      return
      end


c     ****************************************************************
c     *                                                              *
c     *                   subroutine rktstf_set_06                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/18/02                   *
c     *                                                              *
c     *        set up material model #6 (adv. gurson)                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rktstf_set_06( props, iprops, lprops, span,
     &                          local_work )
      use main_data, only : matprp, lmtprp
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                 parameter declarations
c

      real    props(mxelpr,*)
      integer iprops(mxelpr,*)
      logical lprops(mxelpr,*)
c
      matnum = local_work%matnum
c
      do i = 1, span
         local_work%e_v(i)         = props(7,i)
         local_work%nu_v(i)        = props(8,i)
         local_work%sigyld_v(i)    = props(23,i)
         local_work%n_power_v(i)   = props(21,i)
         local_work%f0_v(i)        = props(26,i)
         local_work%q1_v(i)        = props(27,i)
         local_work%q2_v(i)        = props(28,i)
         local_work%q3_v(i)        = props(29,i)
         local_work%nuc_v(i)       = .true.
         local_work%nuc_s_n_v(i)   = props(31,i)
         local_work%nuc_e_n_v(i)   = props(32,i)
         local_work%nuc_f_n_v(i)   = props(33,i)
         if ( iand (iprops(30,i),1) .eq. 0 )
     &             local_work%nuc_v(i) = .false.
         local_work%mm06_props(i,1:5) = matprp(65:69,matnum)

      end do
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                   subroutine rktstf_set_07                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/18/02                   *
c     *                                                              *
c     *        set up material model #7 (mises + hydrogen)           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rktstf_set_07( props, iprops, lprops, span,
     &                          local_work )
      use main_data, only : matprp, lmtprp
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                 parameter declarations
c

      real    props(mxelpr,*)
      integer iprops(mxelpr,*)
      logical lprops(mxelpr,*)
c
      matnum = local_work%matnum
c
      do i = 1, span
c
         local_work%e_v(i)             = props(7,i)
         local_work%nu_v(i)            = props(8,i)
         local_work%sigyld_v(i)        = props(23,i)
         local_work%n_power_v(i)       = props(21,i)
         local_work%mm07_props(i,1:10) = matprp(70:79,matnum)
c
      end do
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                   subroutine rktstf_set_08                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/18/02                   *
c     *                                                              *
c     *        set up material model #8 (warp3d UMAT)                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rktstf_set_08( props, iprops, lprops, span,
     &                          local_work )
      use main_data, only : matprp, lmtprp
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                 parameter declarations
c

      real    props(mxelpr,*)
      integer iprops(mxelpr,*)
      logical lprops(mxelpr,*)
c
c                 no props are needed
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine rktstf_set_10                *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 03/29/12 mcm               *
c     *                                                              *
c     *          set up material model #10 (crystal plasticity) for  *
c     *          tangent stiffness computations                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rktstf_set_10( props, iprops, lprops, span,
     &                          local_work )
      use segmental_curves
      use main_data, only: matprp, lmtprp, imatprp, dmatprp, smatprp
      use crystal_data, only : c_array
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                 parameter declarations
c
      real    props(mxelpr,*)
      integer iprops(mxelpr,*)
      logical lprops(mxelpr,*)
c
      integer :: cnum, matnum,i
c
c           Copy over whatever properties we need
c
      matnum = local_work%matnum
c
      do i=1, span
            local_work%ncrystals(i) = imatprp(101,matnum)
      end do

      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine rktstf_zero_vec                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/18/02                   *
c     *                                                              *
c     ****************************************************************
c
      subroutine rktstf_zero_vec( vec, nterms )
      implicit none
c
#dbl      double precision
#sgl      real
     &   vec(*)
      integer i, nterms
c
      do i = 1, nterms
#sgl         vec(i) = 0.0
#dbl         vec(i) = 0.0d00
      end do
c
      return
      end
