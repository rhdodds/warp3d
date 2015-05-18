c     ****************************************************************
c     *                                                              *
c     *                      subroutine gplns1                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 1/24/13 rhd                *
c     *                                                              *
c     *     process elements in the block for this integration       *
c     *     point. add to element stiffness for this integration pt  *
c     *                                                              *
c     *     Added the shape functions, n() array, added the node     *
c     *     coordinates array, ce().  See the supporting subroutines *
c     *     at the bottom of the file.  (gvt)                        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine gplns1( gpn, icp,
     &                   e_block, nu_block, intf_prp_block, ek,
     &                   gama_block, det_jac_block, nxi, neta, nzeta,
     &                   shape, ce, weight, vol_block, b_block,
     &                   bt_block, bd_block, local_cep_block,
     &                   cohes_rot_block,
     &                   utsz, local_work, nsz)
      use elem_block_data, only : global_cep_blocks => cep_blocks
      use main_data, only: asymmetric_assembly
      implicit integer (a-z)
$add param_def
$add include_lin_ek
c
c                       parameter declaration
c
#dbl      double precision
#sgl      real
     &   ek(nsz,*), e_block(*), nu_block(*),
     &   local_cep_block(mxvl,nstr,*),
     &   vol_block(*),
     &   gama_block(*), det_jac_block(*), nxi(*), neta(*), nzeta(*),
     &   shape(*), ce(mxvl,*),
     &   weight, b_block(*), bt_block(*), bd_block(*),
     &   cohes_rot_block(*), intf_prp_block(*)
c
      integer icp(mxutsz,*)
c
c                       locals
#dbl      double precision
#sgl      real
     &   rad(mxvl), eps_bbar, factors(mxvl), one
c
      logical bbar, local_debug, cohesive_elem, need_adjust,
     &        axisymm_elem, umat_solid, crystal
      data local_debug / .false. /
      data one
#sgl     &  / 1.0 /
#dbl     &  / 1.0d00 /
c
c
      nnode           = local_work%num_enodes
      span            = local_work%span
      totdof          = local_work%totdof
      elem_type       = local_work%elem_type
      eps_bbar        = local_work%eps_bbar
      bbar            = local_work%bbar_flg
      cohesive_elem   = local_work%is_cohes_elem
      need_adjust     = local_work%adjust_const_elem
      axisymm_elem    = local_work%is_axisymm_elem
      mat_type        = local_work%mat_type
      cohes_type      = local_work%cohes_type
      umat_solid      = mat_type .eq. 8
      crystal         = mat_type .eq. 10
c
c                       for axisymmetric elements, compute the radius to
c                       the current gauss point for each element
c                       in the block.
c
      if( axisymm_elem ) call get_radius( rad, nnode, span, shape,
     &                                    ce, mxvl )
c
c                       compute the strain-displacement matrices
c                       for the given gauss point for all elements
c                       in the block.
c
      if ( cohesive_elem ) then
        call blcmp_cohes( span, b_block, cohes_rot_block, shape,
     &                    elem_type, nnode )
      else
        call blcmp1_srt
     &     ( span, b_block, gama_block, nxi, neta, nzeta, shape,
     &             ce, rad, elem_type, nnode )
      end if
      if ( bbar .and. elem_type .eq. 2 )
     &   call bmod ( b_block, vol_block, span, mxvl, eps_bbar, mxedof )
c
c                       set the constitutive matrices for this
c                       gauss point for all elements in block.
c                       each element can have different elastic
c                       properties as specified by the user or
c                       indirectly by temperature dependent
c                       values. we call a separate driver for Abaqus
c                       compatible UMATs to set up all the local data
c                       requried by the UMAT.
c
      if ( cohesive_elem ) then
        ispan  = span
        imxvl  = mxvl
        instr  = nstr
        igpn   = gpn
        ifelem = local_work%felem
        iout   = local_work%iout_local
        istep  = local_work%now_step
        call lcnst4( istep, cohes_type, ispan, imxvl, instr, igpn,
     &               ifelem, iout, local_work%time_n, local_work%dt,
     &               local_cep_block, intf_prp_block,
     &               det_jac_block, weight,
     &               local_work%temps_node_blk(1,1),
     &               local_work%temps_ref_node_blk(1,1) )
        go to 100
      end if
c
      if( crystal ) then
        call lnstff10( gpn, span, local_work )
      else if( umat_solid ) then
        call drive_umat_lnstf( gpn, local_work )
      else
        call lcnst1( span, local_cep_block, nu_block, e_block,
     &               det_jac_block, weight )
      end if                  
c
c                       store cep (i.e. [D]) for all elements in
c                       block at this integration point. scale out
c                       det[J] * weight included above. put in
c                       global blocked data structure
c
 100  continue
      do i = 1, span
         factors(i) = one / ( det_jac_block(i) * weight )
      end do
c
      now_blk = local_work%blk
      call tanstf_store_cep( span, mxvl, gpn, local_work%cep_sym_size,
     &       local_cep_block, global_cep_blocks(now_blk)%vector,
     &       factors )
c
c                       include other required scalars with the
c                       material matrix for axisymmetric (triangle
c                       or quadrillateral) or planar triangle
c                       elements.  A scalar to adjust the Jacobian
c                       determinant is also needed for tetrahedral elements.
c
      if( need_adjust ) call adjust_cnst( elem_type, nstr, mxvl,
     &                                    span, rad, local_cep_block )
c
c                       compute each part of the element linear
c                       stiffness matrices and add it in to the
c                       total. assume the [D] matrix is full.
c
      if (.not. asymmetric_assembly) then
      if ( totdof .eq. 24 ) then
        call bdbt1( span, icp, b_block, bt_block, bd_block,
     &              local_cep_block, ek, mxvl, mxedof, utsz, nstr,
     &              mxutsz )
      else
        call bdbtgen1( span, icp, b_block, bt_block, bd_block,
     &                 local_cep_block, ek, mxvl, mxedof, utsz, nstr,
     &                 totdof, mxutsz )
      end if

      else
        call bdbt_asym( span, icp, local_work%b_block, 
     &                  local_work%bt_block, local_work%bd_block,
     &                  local_work%cep, ek, mxvl, mxedof, utsz,
     &                  nstr, totdof, mxutsz)
      end if
c
      return
c
400   format(1x,//,
     &  '>>>>>>>>>> gplns1 (called each Gauss point) <<<<<<<<<<',/,
     &               '        span = ',i6,/,
     &               '       nnode = ',i6,/,
     &               '      totdof = ',i6,/,
     &               '   elem_type = ',i6,//,
     & ' [n,nxi,neta,nzeta] = [shape function, 3 derivatives]')
410   format(1x,4e14.6)
500   format(1x,//,'>>>>> gplns1.f     axisymmetric element',/,
     &               '      radius to Gauss point:')
510   format(1x,'    elem #',i6,'  radius =',e14.6)
520   format(1x,'    Gauss integration weight = ',f14.8)
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *                     subroutine get_radius                    *
c     *                                                              *
c     *                       written by: gvt                        *
c     *                   last modified : 08/25/98                   *
c     *                                                              *
c     ****************************************************************
c
c         for axisymmetric elements, compute the radius to the current
c         gauss point for each element in the block.  The radius is
c         used for the hoop strain = u/r, and in the summation of the
c         element stiffness matrix. Use the first nnode values in ce()
c         for x=radial coordinates.
c
c         Variables:
c
c          nnode = number of nodes
c          span  = number of elements in the block
c          mxvl  = maximum number of elements in a block
c          rad() = radius to the current Gauss point for each element
c                  in the current block. Using the shape functions
c                  and node coordinates, the radius is given by:
c                  r = Sum[n(i)*x(i)], i=1,nnode where n(i) =
c                  shape functions, x(i) = x=radial coordinate
c          n()   = element shape functions, already evaluated at the curent
c                  gauss integration point
c          ce()  = matrix of node coordinates, the first nnode columns
c                  contain the x=radial coordinates
c
      subroutine get_radius( rad, nnode, span, n, ce, mxvl )
      implicit none
      integer nnode, span, mxvl
#dbl      double precision
#sgl      real
     &         rad(*), n(*), ce(mxvl,*)
c
      integer i,j
#dbl      double precision
#sgl      real
     &         zero
      data zero / 0.0 /
c
      do i = 1, span
         rad(i) = zero
      end do
c
      do j = 1, nnode
          do i = 1, span
            rad(i) = rad(i) + n(j)*ce(i,j)
          end do
      end do
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine bdbt1                        *
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
      subroutine bdbt1( span, icp, b, bt, bd, d, ek,
     &                  mxvl, mxedof, utsz, nstr, mxutsz )
      implicit integer (a-z)
c
c                       parameter declarations
c
#dbl      double precision
#sgl      real
     &   b(mxvl,mxedof,*), ek(utsz,*), d(mxvl,nstr,*),
     &   bd(mxvl,mxedof,*), bt(mxvl,nstr,*)
      integer icp(mxutsz,*)
c
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
c                       full matrix.
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
c     ****************************************************************
c     *                                                              *
c     *                      subroutine bdbtgen1                     *
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
      subroutine bdbtgen1( span, icp, b, bt, bd, d, ek,
     &                     mxvl, mxedof, utsz, nstr, totdof, mxutsz )
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
c                       full matrix.
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
c     ****************************************************************
c     *                                                              *
c     *            subroutine drive_umat_lnstf  (umat)               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 04/8/12                    *
c     *                                                              *
c     *     this routine drives material model 08 to                 *
c     *     compute the linear-elastic [D] for all element in the    *
c     *     block at this integration point                          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine drive_umat_lnstf( gpn, local_work )
c
      implicit integer (a-z)
$add param_def
$add include_lin_ek  ! defines local_work
c
c                      locally defined variables
c
#dbl      double precision
#sgl      real
     &  gp_temps(mxvl), gp_rtemps(mxvl), gp_dtemps(mxvl),
     &  zero, one, ddsddt(6), drplde(6), drpldt,
     &  big, pnewdt, predef(1), dpred(1), time(2), dtime,
     &  stress(6), stran(6), dstran(6), abq_props(50), temp_n, dtemp,
     &  statev(500), sse, spd, scd, coords(3), celent,
     &  dfgrd0(9), dfgrd1(9), drot(9), ddsdde(6,6),
     &  gp_coords(mxvl,3), identity(9), local_cep(6,6),
     &  nonloc_ele_values(nonlocal_shared_state_size)
c
      logical local_debug, debug_now, temperatures_ref, temperatures
      character * 8  cmname
      integer  map(6)
      data map / 1,2,3,4,6,5 /
#sgl      data zero, one, big / 0.0, 1.0, 1.0 /
#dbl      data zero, one, big / 0.0d00, 1.0d00, 1.0d06 /
#sgl      data identity /1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,0.0, 1.0 /
#dbl      data identity /1.0d00, 0.0d00, 0.0d00, 0.0d00, 1.0d00,
#dbl     &   0.0d00, 0.0d00, 0.0d00, 1.0d00 /
c
c           1a. pull a few values from local_work for block
c
      step            = local_work%now_step
      iter            = local_work%now_iter
c
      iout            = local_work%iout_local
      felem           = local_work%felem
      elem_type       = local_work%elem_type
      int_order       = local_work%int_order
      nnode           = local_work%num_enodes
      ndof            = local_work%num_enode_dof
      totdof          = local_work%totdof
      ngp             = local_work%num_int_points
      span            = local_work%span
c
      knumthreads     = local_work%num_threads
      kthread         = omp_get_thread_num() + 1
c
      temperatures_ref  =  local_work%temperatures_ref
      temperatures      =  local_work%temperatures

      local_debug = .false.
c
c           1b. values that remain constant over each element in
c               block. including times.
c
      call material_model_info( felem, 0, 1, hist_size_for_blk )
      ndi    = 3
      nshr   = 3
      ntens  = 6
      nstatv = hist_size_for_blk
      npt    = gpn
      layer  = 1
      kspt   = 1
      kstep  = 1          !  WARP3D has no concept of Abaqus "step"
      kinc   = step       !  Abaqus increments are really WARP3D steps
      kiter  = 1          ! so UMAT will not ignore request
      kout   = local_work%iout_local
      pnewdt = big
c
      dtime   = local_work%dt
      time(1) = zero
      time(2) = dtime
      local_debug = .false.

      if( local_debug ) then
         write(iout,*) ' '
         write(iout,9000)
         write(iout,9001) span, felem, gpn, hist_size_for_blk
         write(iout,9002) step, iter, dtime
         write(iout,9004) temperatures_ref
         write(iout,9008) time(1), time(2)
      end if
c
c           2. interpolate temperatures at material point
c              for elements in block. for step 1, iter 0 this defaults
c              down to the user-specified reference temperatures.
c
      call gauss_pt_temps(
     &   local_work%dtemps_node_blk, gpn, elem_type, span, int_order,
     &   nnode, gp_dtemps, local_work%temps_node_blk,
     &   gp_temps, temperatures, local_work%temps_node_to_process,
     &   temperatures_ref, local_work%temps_ref_node_blk,
     &   gp_rtemps )
c
c
c           3. (x,y,z) undeformed coordinates at this integration
c              point for all elements in the block.
c
      call gauss_pt_coords( gpn, elem_type, span, int_order, nnode,
     &                      gp_coords, local_work%ce, iout )
c
c
c           4. deformation gradients are identity for linear stiffness.
c              zero other items as well.
c
      dfgrd0 = identity
      dfgrd1 = identity
      drot   = identity
      rpl    = zero ! scalar
      ddsddt = zero
      drplde = zero
      drpldt = zero ! scalar
c
      predef(1) = zero; dpred(1) = zero
c
      cmname(1:) = "UMAT-WRP"
c
c               6. drive update over all elements in the block at this
c                  integration point. the Abaqus umat processes 1
c                  integration point per call.
c
      do ielem = 1, span
c
      noel = felem + ielem - 1
      debug_now = local_debug .and. ielem .eq. 1
      if( debug_now )  write(iout,9100) ielem, noel, gpn
c
c               6.1 temperature at n and increment over n -> n+1
c
      dtemp   =  gp_dtemps(ielem)  ! change over step
      temp_n  =  gp_temps(ielem) - gp_dtemps(ielem)
c
c               6.2 set properties and state(history) vectors.
c                   set statev = {0} so umat will return linear [D].
c
      nprops = 50
      abq_props(1:nprops) = local_work%umat_props(ielem,1:nprops)
      statev(1:nstatv) = zero
c
c               6.3 set zero vectors for old stress, old strains and
c                   strain increment. we just want the linear [D]
c
      do j = 1, 6
        stress(j) = zero
        stran(j)  = zero
        dstran(j) = zero
      end do
c
c               6.4 global coordinates of integration point and
c                   characteristic element length per Abaqus spec
c
      coords(1:3) = gp_coords(ielem,1:3)
      celent = local_work%characteristic_length(ielem)
c
c               6.5 umats expect material stiffness to
c                   be initialized zero
c
      ddsdde(1:6,1:6) = zero
c
c               6.6 zero starting values of specific energy,
c                   dissipation.
c
      sse = zero; spd = zero; scd = zero
c
      if( debug_now ) then
       write(iout,9125) coords(1:3), celent
       write(iout,9115) dtemp, temp_n
       write(iout,9110) abq_props(1:10)
      end if
c
c               6.65 zero vector for umat to put nonlocal material
c                   values if it wants
c
      nonloc_ele_values(1:nonlocal_shared_state_size) = zero
c
c               6.7 call the umat for material point
c
      call umat(
     1   stress,
     2   statev,
     3   ddsdde,
     4   sse,
     5   spd,
     6   scd,
     7   rpl, ddsddt, drplde, drpldt,
     8   stran, dstran,
     9   time, dtime,
     a   temp_n, dtemp,
     b   predef, dpred,
     c   cmname,
     d   ndi, nshr, ntens, nstatv,
     e   abq_props, nprops,
     f   coords,
     g   drot,
     h   pnewdt,
     i   celent,
     j   dfgrd0, dfgrd1,
     k   noel, npt,
     l   kslay, kspt, kstep, kinc, kiter, kout, kthread, knumthreads,
     m   nonloc_ele_values, nonlocal_shared_state_size )
c
c               6.8  put linear-elastic [D] into WARP3D local data
c                    structure for the block. swap rows 5, 6 from umat (yz and xz)
c                    which are reversed in WARP3D. make the [D] symmetric
c                    in case UMAT returned non-symemtric version.
c                    scale by integration  weight factor and det of coordinate Jacobian.
c
      if( debug_now ) then
         write(iout,*) '..[D] from umat ... '
         do i = 1, 6
           write(iout,9200) ddsdde(i,1:6)
         end do
      end if
c
      call gplns1_make_symmetric( ddsdde, local_cep )
c
      do j = 1, 6
       do i = 1, 6
         local_work%cep(ielem,i,j) = local_cep(i,j) *
     &    local_work%weights(gpn) * local_work%det_jac_block(ielem,gpn)
        end do
      end do
c
c               6.10 end of lop over all element in block at this
c                    material point.
c
      if( debug_now ) write(iout,*) ' '
      end do
c
      if( debug_now ) write(iout,9900)
c
 9000 format('>> Enter UMAT linear [k] driver...')
 9001 format(5x,'span, felem, gpn, hist_size: ',i4,i10,i3,i5 )
 9002 format(5x,'step, iter, dt: ',i8, i4,e14.6 )
 9004 format(5x,'ref temps defined: ', 2l3 )
 9006 format(5x,'num history terms: ',i4 )
 9008 format(5x,'time_n, delta_t: ', 2e20.9 )
 9100 format(5x,'... processing i, elem, gpn: ',i4,i10, i3)
 9110 format(5x,'props: ',5e14.6,/,12x,5e14.6)
 9115 format(5x,'delta-temperature: ',f10.3,
     &    /, 5x,'temperature @ n:   ',f10.3 )
 9125 format(5x,'coords: ', 3e14.6,/,
     &       5x,'characteristic length: ',e14.6 )
 9200 format(5x,6e14.6)
 9900 format('>> Leave UMAT linear [k] driver...')
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine gplns1_make_symmetric                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 05/14/12                   *
c     *                                                              *
c     *                   support for umat processing                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine gplns1_make_symmetric( ddsdde, cep )
      implicit none
#dbl      double precision
#sgl      real
     & ddsdde(6,6), cep(6,6)
c
#dbl      double precision
#sgl      real
     & transpose(6,6), half
      integer i, j, k, map(6)
#sgl      data half / 0.5 /
#dbl      data half / 0.5d00 /
      data map / 1,2,3,4,6,5 /
c
c         1. compute transpose of 6 x 6 ddsdde from UMAT
c         2. compute symmetrized version
c         3. swap rows, cols 5 & 6 to make shear ordering
c            compatible with WARP3D
c
      do i = 1, 6
        do j = 1, 6
          transpose(i,j) = ddsdde(j,i)
        end do
      end do
c
      do j = 1, 6
        do i = 1, 6
          cep(map(i),map(j)) = half * ( ddsdde(i,j) +
     &                                  transpose(i,j) )
        end do
      end do
c
      return
      end

