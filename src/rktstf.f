c     ****************************************************************
c     *                                                              *
c     *                      subroutine rktstf                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 9/5/2017  rhd              *
c     *                                                              *
c     *     drive computation of tangent stiffness matrices for a    *
c     *     block of similar elements.                               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine rktstf( props, iprops, lprops, glb_ek_blk, nrow_ek,
     &                   ispan, local_work )
      use main_data, only : matprp, lmtprp, asymmetric_assembly
      implicit none
      include 'param_def'
      include 'include_tan_ek'
c
c                 parameter declarations
c
      integer :: nrow_ek  ! set for symmetric or asymmetric
      integer :: ispan
      double precision ::  glb_ek_blk(nrow_ek,ispan)
      real    :: props(mxelpr,*)   !  all 3 are same. read-only here
      integer :: iprops(mxelpr,*)  !  1st col is for felem
      logical :: lprops(mxelpr,*)
c
c                 local
c
      integer :: span, felem, utsz, iter, type, order, ngp, nnode,
     &           ndof, totdof, mat_type, surf, local_iout, gpn,
     &           nlength, i, i_local

      logical :: geonl, local_debug, bbar, first, qbar_flag,
     &           compute_shape, cohes_mirror, dummy_logic,
     &           average, symmetric_assembly
c
      double precision ::
     &  xi, eta, zeta, beta_fact, eps_bbar, zero, one, dummy,
     &  temp_ref, d_temp, temp_np1
c
      data local_debug, dummy_logic / .false., .true. /
      data zero, one, dummy / 0.d0, 1.d0, 1.d0 /
c
c
      local_iout = local_work%iout
      if( local_debug ) write(local_iout,9100)
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
      symmetric_assembly = .not. asymmetric_assembly
c
      if( symmetric_assembly ) then
            allocate( local_work%ek_symm(span,nrow_ek) )
            local_work%ek_symm = zero
            allocate( local_work%ek_full(1,1) )  ! for safe passing
      else
            allocate( local_work%ek_full(span,nrow_ek) )
            local_work%ek_full = zero
            allocate( local_work%ek_symm(1,1) )
      end if
c
c      if( asymmetric_assembly ) ek = zero
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
         call rktstf_zero_vol( local_work%vol_block,
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
     &      call shapef( type, xi, eta, zeta, local_work%shape(1,gpn) )
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
c
      do gpn = 1, ngp
         call gptns1( local_work%cp, local_work%icp, gpn, props,
     &                iprops, glb_ek_blk, local_work )
       end do
c
      if( symmetric_assembly ) then
          call rktstf_do_transpose( local_work%ek_symm, glb_ek_blk,
     &                              span, nrow_ek )
      else
        call rktstf_do_transpose(  local_work%ek_full, glb_ek_blk,
     &                            span, nrow_ek )
      end if
      deallocate( local_work%ek_symm, local_work%ek_full )
c
c               modify element stiffness matrix by thickness
c               factor for plane strain analysis
c
      if( beta_fact .ne. one ) glb_ek_blk = beta_fact * glb_ek_blk
c
c               transform element stiffnesses to constraint
c               compatible coordinates if required. skip
c               whole block or specfic element to reduce work.
c
      if( local_work%trn_e_block ) then
        do i = 1, span
           i_local = i
           if( local_work%trn_e_flags(i) )
     &         call trnmtx( glb_ek_blk, local_work%cp,
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
 9405 format(1x,'   element #',i7)
 9410 format(1x,'j=',i3,2x,'row=',i3,2x,'col=',i3,2x,'ek=',f25.8)
 9300 format(1x,////,
     & '  ELEMENT TESTING, values in rktstf.f',//,
     & '       span = ',i7,/,
     & '      nnode = ',i7,/,
     & '     totdof = ',i7,/,
     & '  elem_type = ',i7,/,
     & '      felem = ',i7,/,
     & '       utsz = ',i7,/)
 9500 format(1x,'>> Fatal Error: rktstf. invalid material type..',
     &    /, 1x,'                job terminated' )
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine rktstf_do_transpose             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 09/19/2016rhd              *
c     *                                                              *
c     ****************************************************************
c

      subroutine rktstf_do_transpose( local_ek, ek, span, nrow_ek )
      implicit none
c
      integer :: span, nrow_ek
      integer :: i, j
      double precision ::
     & local_ek(span,nrow_ek), ek(nrow_ek,span)
c
c                 time 4.59 secs, tanstf time 12.64
c
c       ek = transpose( local_ek )
c       return
c
c
c           time 2.22 secs,  for 1830  tanstf.f 12.73395100
c             when nrow_ek is used, 2.78 secs

c
      select case ( nrow_ek ) ! careful. type could be integer constant
c
      case( 300 )
        do j = 1, span
!DIR$ VECTOR ALIGNED
              do i = 1, 300
                    ek(i,j) = local_ek(j,i)
              end do
        end do
c
      case( 465 )
        do j = 1, span
!DIR$ VECTOR ALIGNED
              do i = 1, 465
                    ek(i,j) = local_ek(j,i)
              end do
        end do
c
      case( 666 )
        do j = 1, span
!DIR$ VECTOR ALIGNED
              do i = 1, 666
                    ek(i,j) = local_ek(j,i)
              end do
        end do
c
      case( 1830 )
        do j = 1, span
!DIR$ VECTOR ALIGNED
              do i = 1, 1830
                    ek(i,j) = local_ek(j,i)
              end do
        end do
c
      case default
        do j = 1, span
!DIR$ VECTOR ALIGNED
              do i = 1, nrow_ek
                    ek(i,j) = local_ek(j,i)
              end do
        end do
      end select
c
c
c
c                  time 5.16 secs
c
c      do i = 1, 1830
c            do j = 1, span
c                  ek(i,j) = local_ek(j,i)
c            end do
c      end do


      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine rktstf_zero_vec                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/8/2016 rhd               *
c     *                                                              *
c     ****************************************************************
c
      subroutine rktstf_zero_vec( vec, nterms )
      implicit none
c
      integer :: nterms
      double precision :: vec(nterms), zero
      data zero / 0.0d00 /
c
!DIR$ VECTOR ALIGNED
      vec = zero
c
      return
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
      subroutine rktstf_zero_vol ( vol, volume, span, mxvl )
      implicit none
      integer :: span, mxvl
      double precision  :: vol(mxvl,8,3), volume(span)
      double precision, parameter :: zero = 0.0d0
      vol = zero
      volume = zero
      return
      end

