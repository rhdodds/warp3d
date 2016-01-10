c     ****************************************************************
c     *                                                              *
c     *                      subroutine rktstf                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 1/10/2016 rhd              *
c     *                                                              *
c     *     drive computation of tangent stiffness matrices for a    *
c     *     block of similar elements.                               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rktstf( props, iprops, lprops, ek, nrow_ek, ispan,
     &                   local_work )
      use main_data, only : matprp, lmtprp
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c                 parameter declarations
c
#dbl      double precision ::
#sgl      real ::
     &   ek(nrow_ek,ispan)
      real    :: props(mxelpr,*)   !  all 3 are same. read-only here
      integer :: iprops(mxelpr,*)
      logical :: lprops(mxelpr,*)
c
c                 local
c
      logical :: geonl, local_debug, bbar, first, qbar_flag,
     &           compute_shape, cohes_mirror, dummy_logic,
     &           average
c
#dbl      double precision ::
#sgl      real ::
     &  xi, eta, zeta, beta_fact, eps_bbar, zero, one, dummy,
     &  temp_ref, d_temp, temp_np1
c
      data local_debug, dummy_logic / .false., .true. /
      data zero, one, dummy / 0.d0, 1.d0, 1.d0 /
c
      local_iout = local_work%iout
      if( local_debug ) write(out,9100)
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
      call rktstf_zero_vec( ek, nrow_ek*span )
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
      do gpn = 1, ngp
        call gptns1( local_work%cp, local_work%icp, gpn, props,
     &               iprops, ek, local_work, nrow_ek )
      end do
c
c		             modify element stiffness matrix by thickness factor for
c             	plane strain analysis
c
      if( beta_fact .ne. one ) ek = beta_fact * ek
c
c               transform element stiffnesses to constraint
c               compatible coordinates if required. skip
c               whole block or specfic element to reduce work.
c
      if( local_work%trn_e_block ) then
        do i = 1, span
           i_local = i
           if( local_work%trn_e_flags(i) )
     &         call trnmtx( ek, local_work%cp,
     &                      local_work%trnmte, local_work%trne, ndof,
     &                      nnode, totdof, i_local, nrow_ek )
         end do
      end if
c
      return
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
