c     ****************************************************************
c     *                                                              *
c     *                      subroutine lnstff                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 9/28/2015 rhd              *
c     *                                                              *
c     *     drive linear stiffness computation for all elements.     *
c     *     assemble diagonal stiffness vector for strucure or       *
c     *     domain (MPI)                                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine lnstff( now_step, now_iter )
c
      use elem_block_data, only : estiff_blocks, edest_blocks
      use main_data,       only : umat_serial, 
     &                            asymmetric_assembly
      implicit integer (a-z)
$add common.main
c
c             local declarations
c
#dbl      double precision ::
#sgl      real ::
     &  zero, mag
      logical :: local_debug, blks_reqd_serial, umat_matl
      allocatable :: blks_reqd_serial(:)      
      data local_debug, zero / .false., 0.0d00 /
c
c             if MPI:
c                tell the worker processors to join us in this
c                routine. Then compute number of degrees of
c                freedom access by the processor.
c
      call wmpi_alert_slaves ( 3 )
      call wmpi_bcast_int ( now_step )
      call wmpi_bcast_int ( now_iter )
c
      call thyme( 3, 1 )
c
c             compute linear stiffnesses for each
c             block of elements in structure or domain.
c             Stiffness is based on elastic properties
c             for temperature at end of step for
c             zero nodal displacements.
c
c             we compute the element stiffnesses in a local block
c             array then copy to globally allocated array. reduces
c             access to globals from threads.
c
      call estiff_allocate ( 1 )
c
c             if MPI:
c                elblks(2,blk) holds which processor owns the
c                block.  If we don't own the block, then skip its
c                computation. otherwise elblks(2,blk) = 0
c
c                set logical flag vector to indicate which if any
c                element blocks must be run in a serial version of
c                the block loop. E.g. element blocks using a
c                umat where the user declared the umat must run
c                serial.
c
      allocate( blks_reqd_serial(nelblk) )
c
      do blk = 1, nelblk
        blks_reqd_serial(blk)   = .false.
        felem                   = elblks(1,blk)
        mat_type                = iprops(25,felem)
        umat_matl               = mat_type .eq. 8
        if( umat_matl .and. umat_serial )
     &       blks_reqd_serial(blk)   = .true.
      end do
c
      call omp_set_dynamic( .false. )
c$OMP PARALLEL DO PRIVATE( blk, now_thread )
c$OMP&            SHARED( now_step, now_iter, nelblk, elblks )
      do blk = 1, nelblk
         if( elblks(2,blk) .ne. myid ) cycle
         if( blks_reqd_serial(blk) ) cycle
         now_thread = omp_get_thread_num() + 1
         call do_lnek_block( now_step, now_iter, blk )
      end do
c$OMP END PARALLEL DO
c
c             serial version of the block loop to
c             catch non-thread safe code (possibly umats)
c
      do blk = 1, nelblk
         if( elblks(2,blk) .ne. myid ) cycle
         if( .not. blks_reqd_serial(blk) ) cycle
         now_thread = omp_get_thread_num() + 1
         call do_lnek_block( now_step, now_iter, blk )
      end do
c
c             flags indicating that the linear stiffness
c             matrix has been calculated.
c
      lkcomp = .true.
      tkcomp = .false.
      call thyme( 3, 2 )
c
      return
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine do_lnek_block                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/28/2015 rhd              *
c     *                                                              *
c     *     drive computation of linear stiffness matrices for a     *
c     *     block of elements                                        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine do_lnek_block( now_step, now_iter, blk )
c
      use elem_block_data,    only : estiff_blocks, cdest_blocks,
     &                               edest_blocks
      use elem_extinct_data,  only : dam_blk_killed, dam_state
      use main_data,          only : incmap, incid, 
     &                               temper_nodes, temper_elems,
     &                               temperatures_ref,
     &                               fgm_node_values_defined,
     &                               cohesive_ele_types,
     &                               linear_displ_ele_types,
     &                               adjust_constants_ele_types,
     &                               axisymm_ele_types,
     &                               asymmetric_assembly
c
      use damage_data, only : dam_ptr, growth_by_kill
      use contact, only : use_contact     
c
      implicit integer (a-z)
$add common.main
c
c             local declarations
c
$add include_lin_ek
#dbl      double precision ::
#sgl      real ::
     &  zero
      logical :: local_debug, geo_non_flg, bbar_flg
      data local_debug, zero / .false., 0.0d00 /
c
      felem          = elblks(1,blk)
      elem_type      = iprops(1,felem)
      int_order      = iprops(5,felem)
      mat_type       = iprops(25,felem)
      num_enodes     = iprops(2,felem)
      num_enode_dof  = iprops(4,felem)
      totdof         = num_enodes * num_enode_dof
      geo_non_flg    = lprops(18,felem)
      bbar_flg       = lprops(19,felem)
      num_int_points = iprops(6,felem)
      span           = elblks(0,blk)
      utsz           = ((totdof*totdof)-totdof)/2 + totdof
      cohes_type     = iprops(27,felem)
      matnum         = iprops(38,felem)
c
      local_work%felem          = felem
      local_work%blk            = blk
      local_work%num_threads    = num_threads
      local_work%elem_type      = elem_type
      local_work%int_order      = int_order
      local_work%mat_type       = mat_type
      local_work%num_enodes     = num_enodes
      local_work%num_enode_dof  = num_enode_dof
      local_work%totdof         = totdof
      local_work%geo_non_flg    = geo_non_flg
      local_work%bbar_flg       = bbar_flg
      local_work%num_int_points = num_int_points
      local_work%span           = span
      local_work%now_step       = now_step
      local_work%now_iter       = now_iter
      local_work%dt             = dt
      local_work%time_n         = total_model_time
      local_work%temperatures_ref = temperatures_ref
      local_work%temperatures     = temperatures ! from common.main
      local_work%matnum         = matnum
      local_work%utsz           = utsz
      local_work%beta_fact      = beta_fact
      local_work%eps_bbar       = eps_bbar
      local_work%cohes_type     = cohes_type
      local_work%fgm_enode_props = fgm_node_values_defined
      local_work%is_cohes_elem  = cohesive_ele_types(elem_type)
      local_work%linear_displ_elem = linear_displ_ele_types(elem_type)
      local_work%adjust_const_elem =
     &                             adjust_constants_ele_types(elem_type)
      local_work%is_axisymm_elem = axisymm_ele_types(elem_type)
      local_work%iout_local     = out
      local_work%cep_sym_size   = 21
      if( local_work%is_cohes_elem ) local_work%cep_sym_size = 6
c
c             if all elements in block killed skip
c             all calculations. just zero the stiffnesses for block.
c
      nrow_ek = utsz
      if( asymmetric_assembly ) nrow_ek = totdof**2
      if( growth_by_kill ) then
        if( dam_blk_killed(blk) ) then
          if( local_debug ) write(*,*) 'blk ',blk,' killed - skip.'
          call lnstff_zero_vector( estiff_blocks(blk)%ptr(1,1),
     &                             nrow_ek*span )
          return
        end if
      end if
c
c             build data structures for elements in this block.
c             this is a gather operation on nodal coordinates,
c             constraint transformation, etc. the gather is into
c             local (stack) data structures.
c
      if( local_debug )
     &  write(out,9100) blk, span, felem, mat_type, num_enodes,
     &                  num_enode_dof, totdof, num_int_points
c
      call dplstf( cdest_blocks(blk)%ptr(1,1), totdof,
     &             incid(incmap(felem)), num_enodes, local_work )
c
c             compute linear stiffness for each element
c             in the block. symmetric element stiffnesses are stored
c             in upper triangular form.
c
      if( local_debug )
     &          write(out,9200) blk, span, felem, elem_type, int_order,
     &                               geo_non_flg, bbar_flg
c
      call rklstf( local_work, props(1,felem), props(1,felem),
     &             estiff_blocks(blk)%ptr(1,1), nrow_ek )
c
c             check if this block has any killed elements -- if so,
c             zero computed linear stifffness matrices for killed
c             elements.
c
      if( growth_by_kill ) then
        block_is_killable = iand( iprops(30,felem),2 ) .ne. 0
        if( block_is_killable ) then
@!DIR$ LOOP COUNT MAX=###  
          do relem = 1, span
           element = felem + relem - 1
           if( dam_ptr(element) .eq. 0 ) cycle
           if( dam_state(dam_ptr(element) ) .ne. 0 )
     &     call lnstff_zero_vector( estiff_blocks(blk)%ptr(1,relem),
     &                              nrow_ek )
         end do
        end if 
      end if 
c
c              if contact is included in this analysis, then add
c              the contact spring stiffnesses into the corresponding
c              element stiffness matricies.
c
      if( use_contact )
     &   call contact_stfadd (span, felem, totdof,
     &     edest_blocks(blk)%ptr(1,1),
     &     estiff_blocks(blk)%ptr(1,1), nrow_ek, num_enodes,
     &     incid(incmap(felem)) )
c
c              release all the allocated data structures private to
c              this block
c
      deallocate( local_work%ce, local_work%trnmte,
     &  local_work%det_jac_block,
     &  local_work%shape,
     &  local_work%nxi,
     &  local_work%neta,
     &  local_work%nzeta,
     &  local_work%gama )
      deallocate(
     &  local_work%vol_block,
     &  local_work%jac_block,
     &  local_work%b_block,
     &  local_work%bt_block,
     &  local_work%bd_block,
     &  local_work%cep,
     &  local_work%cohes_rot_block,
     &  local_work%intf_prp_block,
     &  local_work%temps_node_blk,
     &  local_work%dtemps_node_blk,
     &  local_work%temps_ref_node_blk )
      deallocate(
     &  local_work%enode_mat_props,
     &  local_work%fgm_flags,
     &  local_work%umat_props,
     &  local_work%trne )
c
      if( local_work%mat_type .eq. 10 ) then
            deallocate(local_work%cp_stiff)
            deallocate(local_work%ncrystals)
            deallocate(local_work%cp_g_rot)
      end if

      return
c
 9100 format(5x,'>>> ready to call dplstf:',
     &     /,10x,'blk, span, felem, mat_model:      ',4i10,
     &     /,10x,'num_enodes, num_enode_dof, totdof:',3i10,
     &     /,10x,'num_int_points:                   ',i10 )
 9200 format(5x,'>>> ready to call rklstf:',
     &     /,10x,'blk, span, felem                 :',3i10,
     &     /,10x,'elem_type, int_order, geo_non_flg:',2i10,l10,
     &     /,10x,'bbar_flg:                        :',l10 )
c
      end


c     ****************************************************************
c     *                                                              *
c     *                      subroutine dplstf                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 4/21/2014 rhd              *
c     *                                                              *
c     *     set up block local data structure for element linear     *
c     *     stiffness computations.                                  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine dplstf( bcdst, totdof, belinc, num_enodes, local_work )
c
      use main_data, only: fgm_node_values, fgm_node_values_defined,
     &                     fgm_node_values_cols, temperatures_ref,
     &                     temper_nodes, temper_elems, trn,
     &                     temper_nodes_ref, dtemp_nodes, dtemp_elems
c
      use elem_block_data, only : history_blk_list, history_blocks
      use main_data, only : matprp, lmtprp, imatprp, dmatprp, smatprp
      use crystal_data, only : c_array, angle_input, crystal_input,
     &                              data_offset
c
      implicit integer (a-z)
$add common.main
$add include_lin_ek
c
c             parameter declarations (global, read-only arrays)
c
      dimension bcdst(totdof,*), belinc(num_enodes,*)
c
c             local declarations
c
      logical :: local_debug
#dbl      double precision ::
#sgl      real ::
     &  zero, one
      data local_debug, zero, one / .false., 0.0d00, 1.0d00 /
      integer :: matnum, ci, cnum, elnum, osn, ati, aci, ngp, co, blk,
     &            hist_sz
      character :: aconv*5
      character :: atype*7
      double precision :: angles(3), eye(3,3)
c       
      eye = zero
      eye(1,1) = one; eye(2,2) = one; eye(3,3) = one
c
      now_step = local_work%now_step
      now_iter = local_work%now_iter
      span     = local_work%span
      felem    = local_work%felem
      matnum   = local_work%matnum
      ngp      = local_work%num_int_points
      blk      = local_work%blk
c
c             allocate large data structures that hang from local_work.
c             these previously were not allocated and thus went on the
c             runtime stack - making it very large.
c
      if( local_debug ) write(out,9000)
     &          now_step, now_iter, span, felem, num_enodes
      allocate( local_work%ce(mxvl,mxecor) )
      local_work%ce = zero
      allocate( local_work%trnmte(mxvl,mxedof,mxndof) )
      allocate( local_work%det_jac_block(mxvl,mxgp) )
      allocate( local_work%shape(mxndel,mxgp) )
      allocate( local_work%nxi(mxndel,mxgp) )
      allocate( local_work%neta(mxndel,mxgp) )
      allocate( local_work%nzeta(mxndel,mxgp) )
      allocate( local_work%gama(mxvl,3,3,mxgp) )
      allocate( local_work%vol_block(mxvl,8,3) )
      allocate( local_work%jac_block(mxvl,3,3) )
      allocate( local_work%b_block(mxvl,mxedof,nstr) )
      allocate( local_work%bt_block(mxvl,nstr,mxedof) )
      allocate( local_work%bd_block(mxvl,mxedof,nstr) )
      allocate( local_work%cep(mxvl,6,6) )  ! even for cohesive elems
      allocate( local_work%cohes_rot_block(mxvl,3,3) )
      allocate( local_work%intf_prp_block(mxvl,max_interface_props) )
      allocate( local_work%temps_node_blk(mxvl,mxndel) )
      allocate( local_work%dtemps_node_blk(mxvl,mxndel) )
      allocate( local_work%temps_ref_node_blk(mxvl,mxndel) )
      allocate( local_work%enode_mat_props(mxndel,mxvl,mxndpr) )
      allocate( local_work%fgm_flags(mxvl,mxndpr) )
      allocate( local_work%umat_props(mxvl,50) )
      allocate( local_work%trne(mxvl,mxndel) )
c
c             transformation matrices are used to define a "local"
c             coordinate systen for imposition of constraints, e.g.,
c             skewed. if transformations are present, they are used
c             to rotate the computed element stiffness terms for those
c             nodes with constraint coordinate systems. all element
c             computations are performed in global coordinates.
c
c             build table of element nodal coordinates and
c             transformation flags for block.
c             this is done for each node of each element in the block.
c             we then set a unique flag for each element
c             and for the whole block (which can be skipped most times).
c
      local_work%trn_e_flags(1:span) = .false.
      local_work%trn_e_block = .false.
      k = 1
      do j = 1, num_enodes
@!DIR$ LOOP COUNT MAX=###  
        do i = 1, span
          local_work%ce(i,k)        = c(bcdst(k,i))
          local_work%ce(i,k+1)      = c(bcdst(k+1,i))
          local_work%ce(i,k+2)      = c(bcdst(k+2,i))
          local_work%trne(i,j)      = trn(belinc(j,i))
          local_work%trn_e_flags(i) = local_work%trn_e_flags(i) .or.
     &                                local_work%trne(i,j)
          local_work%trn_e_block    = local_work%trn_e_block .or.
     &                                local_work%trne(i,j)
        end do
        k = k + 3
      end do
c
      if( local_work%is_cohes_elem ) then
        if( local_work%elem_type .eq. 14 ) then
          call check_cohes_quad( local_work%elem_type, span, mxvl,
     &                        local_work%felem, local_work%num_enodes,
     &                        local_work%iout_local, local_work%ce )
        else
          call check_cohes_tri( local_work%elem_type, span, mxvl,
     &                        local_work%felem, local_work%num_enodes,
     &                        local_work%iout_local, local_work%ce )
        end if
       end if
c
c             gather element-level transformation matrices
c
      if( local_work%trn_e_block )
     &        call duptrans( span, felem, local_work%trnmte )
c
c             make copies of vectors of subscripts
c             used in matrix multiplies for element
c             stiffnesses. this reduces access to shared
c             global variables during threaded processing
c             of element blocks.
c
      local_work%cp(1:mxedof) = cp(1:mxedof)
      do i = 1, mxutsz
       local_work%icp(i,1) = icp(i,1)
       local_work%icp(i,2) = icp(i,2)
      end do
c
c             build nodal temperatures for elements in the block
c             at start of step and increment over step. temperatures
c             are needed to support temperature dependent, linear
c             material properties.
c
c             we can have reference temps, temps at n and dtemps
c             over n -> n+1.
c
c             temper_nodes, temper_elems: global variables for
c             temperatures @ n. temper_nodes are loaded with
c             user-specified reference values at n = 0
c             dtemp_nodes, dtemp_elems: global changes over step
c
c             temperatures flag set by loads processor to indicate
c             a specified temperature change over the step.
c
      do j = 1, num_enodes
@!DIR$ LOOP COUNT MAX=###  
        do i = 1, span
          local_work%dtemps_node_blk(i,j)    = zero
          local_work%temps_ref_node_blk(i,j) = zero
        end do
      end do

      if( temperatures ) then ! global flag
        if( local_debug )  write(out,9622)
        call gadtemps( dtemp_nodes, dtemp_elems(felem), belinc,
     &                 num_enodes, span, felem,
     &                 local_work%dtemps_node_blk, mxvl )
      end if
c
c             gather reference temperatures for element nodes
c             from the global vector of reference values
c             (if they are defined). we construct a set of reference
c             nodal temperatures for each element in block.
c
      if( temperatures_ref ) then ! global flag
        if( local_debug )  write(out,9610)
        call gartemps( temper_nodes_ref, belinc, num_enodes, span,
     &                 felem, local_work%temps_ref_node_blk, mxvl )
      end if
c
c             build nodal temperatures for elements in the block
c             at end of step (includes both imposed nodal and element
c             temperatures). set logical flag is lower routines
c             need to handles temperatures.
c
      if( local_debug )  write(out,9630)
      call gatemps( temper_nodes, temper_elems(felem), belinc,
     &              num_enodes, span, felem,
     &              local_work%temps_node_blk, mxvl,
     &              local_work%dtemps_node_blk,
     &              local_work%temps_node_to_process )
c
      local_work%temps_node_to_process = .false.
      enode_loop: do j = 1, num_enodes
@!DIR$ LOOP COUNT MAX=###  
       do i = 1, span
        if( abs( local_work%temps_node_blk(i,j) ) .gt. zero ) then
           local_work%temps_node_to_process = .true.
           exit enode_loop
        end if
        end do
       end do enode_loop

c
      if( local_debug ) then
        write(out,*) "...  temps_node_to_process: ",
     &                 local_work%temps_node_to_process
        write(out,*) "... element block dtemps temps:"
        do i = 1, span
          write(out,9100) i, local_work%dtemps_node_blk(i,1:num_enodes)
        end do
        write(out,*) " "
        write(out,*) "... element block node ref temps:"
        do i = 1, span
         write(out,9100) i,
     &          local_work%temps_ref_node_blk(i,1:num_enodes)
        end do
        write(out,*) " "
        write(out,*) "... element block node temps:"
        do i = 1, span
          write(out,9100) i, local_work%temps_node_blk(i,1:num_enodes)
        end do
      end if
c
c             if the model has fgm properties at the model nodes,
c             build a table of values for nodes of elements in the
c             block
c
      if ( fgm_node_values_defined ) then  ! global flag
        do j = 1, fgm_node_values_cols
         do i = 1, num_enodes
@!DIR$ LOOP COUNT MAX=###  
          do k = 1, span
           local_work%enode_mat_props(i,k,j) =
     &               fgm_node_values(belinc(i,k),j)
          end do
         end do
        end do
      end if
c
c
c             if the model has fgm properties at the model nodes,
c             allocate and set up some data structures specific to 
c             the CP materials need the initial, anisotropic stiffness,
c             the crystallographic rotations, and the elastic 
c             rotations (a history variable...).
c
      if( local_work%mat_type .ne. 10) return
        allocate(local_work%cp_stiff(mxvl, 6, 6, max_crystals))
        allocate(local_work%cp_g_rot(mxvl, 3, 3, max_crystals))
        allocate(local_work%ncrystals(mxvl))
        do i = 1, span
          local_work%ncrystals(i) = imatprp(101,matnum)
          elnum = felem+i-1
          do ci =1, local_work%ncrystals(i)
c                 Get the local crystal number
                  if (imatprp(104,matnum) .eq. 1) then
                        cnum = imatprp(105,matnum)
                  elseif (imatprp(104,matnum) .eq. 2) then
                        osn = data_offset(elnum)
                        cnum = crystal_input(osn,ci)
c                       Couldn't do this earlier, so check here
                        if ((cnum .gt. max_crystals) .or.
     &                        (cnum .lt. 0)) then
                         write (out,'("Crystal ", i3, " not valid")')
     &                        cnum
                              call die_gracefully
                        end if
                  else
                        write(out,9502)
                        call die_gracefully
                  end if
                  local_work%cp_stiff(i,1:6,1:6,ci) = 
     &                  c_array(cnum)%elast_stiff
c
c                       Get the local orientation
                  if (imatprp(107,matnum) .eq. 1) then
                        angles(1) = dmatprp(108,matnum)
                        angles(2) = dmatprp(109,matnum)
                        angles(3) = dmatprp(110,matnum)
                  elseif (imatprp(107,matnum) .eq. 2) then
                        osn = data_offset(elnum)
                        angles(1:3) = angle_input(osn,ci,1:3)
                  else
                        write (out,9502)
                        call die_gracefully
                  end if
                  aci = imatprp(102,matnum)
                  ati = imatprp(103,matnum)
c                 Call a helper to get the crystal -> reference rotation
                  if (ati .eq. 1) then
                        atype = "degrees"
                  elseif (ati .eq. 2) then
                        atype = "radians"
                  else
                        write(out,9503)
                        call die_gracefully
                  end if

                  if (aci .eq. 1) then
                        aconv="kocks"
                  else
                        write(out,9504)
                        call die_gracefully
                  end if
                  call mm10_rotation_matrix(angles,
     &                  aconv, atype,
     &                  local_work%cp_g_rot(i,1:3,1:3,ci), 
     &                  local_work%iout_local)
            end do ! over ncrystals
      end do   !   over span
      return
c
 9000 format(/,'>> running dplstf ...',
     &  /,5x,'step, iter, span, felem, num_enodes: ',i8,i4,i4,i8,i3 )
 9100 format(5x,i3,1x,10(f7.2,1x))
 9610 format(12x,'>> dplstf: set reference temps...' )
 9620 format(12x,'>> dplstf: set dtemps = 0 for step 1...' )
 9622 format(12x,'>> dplstf: set dtemps for step ne 1...' )
 9630 format(12x,'>> dplstf: set final element temps...' )
 9502 format(12x,'>> dplstf: bad crystal input...')
 9503 format(/,1x,'>>>> System error: unexpected angle type in dplstf!',
     &            ' Aborting.'/)
 9504 format(/,1x,'>>>> System error: unexpected angle conv in dplstf!',
     &            ' Aborting.'/)
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *               subroutine return_rotation                     *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 10/4/13 mcm                *
c     *                                                              *
c     *     Return the rotation from a history vector...             *
c     *                                                              *
c     ****************************************************************
c
      subroutine return_rotation(history_in, ngp, hs, span, co, gp, i,
     &      R)
      implicit none
c      
      integer :: ngp, hs, span, co, gp , i
      double precision, dimension(3,3) :: R
#dbl      double precision
#sgl      real
     & history_in(hs,ngp,span)
c
      R = reshape(history_in((co+12):(co+21), gp, i), (/3,3/))
c
      end 

c     ****************************************************************
c     *                                                              *
c     *               subroutine lnstff_zero_vector                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/23/15 rhd               *
c     *                                                              *
c     *     zero a vector of specified length w/ floating zero       *
c     *     signle or double based on this port                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine lnstff_zero_vector( vec, n )
#dbl      double precision
#sgl      real
     &  vec(*), zero
      data zero / 0.0d00 /
      vec(1:n) = zero
      return
      end

