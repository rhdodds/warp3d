c     *                                                              *
c     *                drive_eps_sig_internal_forces                 *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 12/21/2017 rhd             *
c     *                                                              *
c     *      recovers all the strains, stresses                      *
c     *      and internal forces (integral B-transpose * sigma)      *
c     *      for all the elements in the structure at state (n+1).   *
c     *                                                              *
c     *      The element internal force vectors are assembled        *
c     *      into the structure level internal force vector.         *
c     *      Blocks are processed in parallei; using threads (OMP),  *
c     *      domains in parallel using MPI                           *
c     *                                                              *
c     ****************************************************************
c
c
c
c              files containing code in this tree:
c                rstgp1.f, dupstr.f rknstr.f gtlsn1.f rplstr.f
c
c
c
      subroutine drive_eps_sig_internal_forces( step, iter,
     &                                          material_cut_step )
      use global_data ! old common.main
c
      use elem_block_data,   only : einfvec_blocks, edest_blocks
      use elem_extinct_data, only : dam_blk_killed, dam_ifv, dam_state
      use damage_data, only : growth_by_kill
      use main_data, only: umat_serial
c
      implicit none
c
      integer :: step, iter
      logical :: material_cut_step
c
c             locals
c
      integer :: blk, felem, mat_type, now_thread, num_enodes,
     &           num_enode_dof, span, i, j, kkk, totdof, ksize
      integer :: num_term_ifv_threads(max_threads),
     &           idummy1(1), idummy2(1)
      integer, external :: omp_get_thread_num
c
      double precision, allocatable :: block_energies(:),
     &                                 block_plastic_work(:),
     &                                 block_stress_norm2s(:),
     &                                 block_strain_norm2s(:)
      double precision :: mag, start_time, sum_stress_norm2s,
     &                    sum_strain_norm2s,
     &                    end_time, sum_ifv_threads(max_threads)
      double precision, parameter :: zero = 0.0d0, one = 1.0d0
      double precision, external :: omp_get_wtime
      logical, allocatable, dimension(:)  :: step_cut_flags,
     &                                       blks_reqd_serial
c
      logical :: umat_matl, run_serial_loop
      logical, parameter :: local_debug = .false.,
     &                      local_debug_sums = .false.
c
      integer :: init_file_no, ielem
      logical :: initsig
      integer, external :: warp3d_get_device_number
c
c             MPI:
c               alert workers we are in the strain-
c               stress-internal force routines.
c
      call wmpi_alert_slaves( 14 )
      call wmpi_bcast_int( step )
      call wmpi_bcast_int( iter )
c
      call thyme( 6, 1 )
c
c             recover and replace the stresss for each
c             block of elements. compute internal forces
c             for all elements using updated stresses.
c
c             features for elements in a block are obtained by
c             examining the first element in the block.
c
c             note: iter = 0 doing the psuedo-stress update at the
c             start of a step to support non-zero imposed displacements
c             and/or extrapolated displacement increment, and/or
c             imposed temperature increments. do not scatter state
c             variables back to structure data for this situation.
c
      if( local_debug ) then
         write(out,*) '>>>> entered drive_eps_sig_internal_forces: ',
     &            ' strain/stress update'
         write(out,*) '     step, iter: ',step, iter
      end if
c
c             allocate and initialize vectors of
c             internal energy values and logicals
c             for step cut parameters. Each block
c             sets its value, then we build structure
c             value when blocks are done. This approach
c             allows (threaded) parallel processing of blocks
c             without conflicts. allocate blocks of
c             element internal force vectors. we compute
c             the blocks of vectors in parallel then
c             scatter in serial to eliminate conflicts.
c
c             some element blocks may need to be done in
c             serial version of the block loop. for now that
c             is blocks with umat material where the user specified
c             the umat must run serial.
c
      allocate ( block_energies(nelblk), block_plastic_work(nelblk),
     &           step_cut_flags(nelblk), blks_reqd_serial(nelblk) )
      ksize = 1
      if( local_debug_sums ) ksize = nelblk
      allocate( block_stress_norm2s(ksize),
     &          block_strain_norm2s(ksize) ) !  for checking MPI
c
      do blk = 1, nelblk
        block_energies(blk)     = zero
        block_plastic_work(blk) = zero
        step_cut_flags(blk)     = .false.
        blks_reqd_serial(blk)   = .false.
        felem                   = elblks(1,blk)
        mat_type                = iprops(25,felem)
        umat_matl               = mat_type .eq. 8
        if( umat_matl .and. umat_serial ) blks_reqd_serial(blk) = .true.
      end do

c
c             allocate blocks of element internal force vectors.
c             we compute the blocks of vectors in parallel then
c             scatter in serial to eliminate conflicts.
c
      call allocate_ifv( 1 )
c
c             update element strains, stresses, internal forces.
c             MPI:
c               elblks(2,blk) holds which processor owns the
c               block.  If we don't own the block, then skip its
c               computation.
c             Else:
c               elblks(2,blk) is all equal to 0, so all blocks
c               are processed.
c
c
c             Process blocks in parallel with threads. the block
c             data structures are all designed to support this
c             high-level parallel operation.
c
      call omp_set_dynamic( .false. )
      if( local_debug ) start_time = omp_get_wtime()
      run_serial_loop = .false.
c
c$OMP PARALLEL DO PRIVATE( blk, now_thread )
c$OMP&            SHARED( nelblk, elblks, myid, iter, step,
c$OMP&                    step_cut_flags, block_energies,
c$OMP&                    block_plastic_work )
      do blk = 1, nelblk
         if( elblks(2,blk) .ne. myid ) cycle
         if( blks_reqd_serial(blk) ) then
c$OMP ATOMIC WRITE
            run_serial_loop = .true.
            cycle
         end if
         now_thread = omp_get_thread_num() + 1
         call do_nleps_block( blk, iter, step, step_cut_flags(blk),
     &                        block_energies(blk),
     &                        block_plastic_work(blk),
     &                        block_stress_norm2s(1),
     &                        block_strain_norm2s(1),
     &                        local_debug_sums )
      end do
c$OMP END PARALLEL DO
c
c             now run a serial version of the block loop to
c             catch left over blocks that must be run
c             serial.
c
      if( run_serial_loop ) then
       do blk = 1, nelblk
         if( elblks(2,blk) .ne. myid ) cycle
         if( .not. blks_reqd_serial(blk) ) cycle
         now_thread = omp_get_thread_num() + 1
         call do_nleps_block( blk, iter, step, step_cut_flags(blk),
     &                        block_energies(blk),
     &                        block_plastic_work(blk),
     &                        block_stress_norm2s(1),
     &                        block_strain_norm2s(1),
     &                        local_debug_sums )
       end do
      end if
c
      if( local_debug ) then
         end_time = omp_get_wtime()
         write(out,*) '>> threaded eps-sig-ifv: ', end_time - start_time
      end if
c
c             if any of the element blocks requested a step
c             size reduction during the stress update,
c             set the controlling flag passed in. also
c             build total internal energy. deallocate space.
c
      material_cut_step = .false.
      internal_energy   = zero
      plastic_work      = zero
      sum_stress_norm2s = zero      ! for MPI checking
      sum_strain_norm2s = zero      ! for MPI checking
      do blk = 1, nelblk
         if( elblks(2,blk) .ne. myid ) cycle
         if( step_cut_flags(blk) )  material_cut_step = .true.
         internal_energy = internal_energy + block_energies(blk)
         plastic_work    = plastic_work + block_plastic_work(blk)
         if( local_debug_sums ) then
            sum_stress_norm2s = sum_stress_norm2s +
     &                           block_stress_norm2s(blk)
            sum_strain_norm2s = sum_strain_norm2s +
     &                           block_strain_norm2s(blk)
         end if
      end do

      if( myid .eq. -1 ) then   !  save for future MPI debug

        write(out,*) '....  myid = 5, block norm2 stresses'
        do blk = 1, nelblk
           if( elblks(2,blk) .ne. 5 ) cycle
           write(out,8900)  blk, block_strain_norm2s(blk),
     &           block_stress_norm2s(blk)
        end do
        end if
c
      deallocate( step_cut_flags, block_energies, block_plastic_work,
     &            blks_reqd_serial, block_stress_norm2s,
     &            block_strain_norm2s )
c
c             For MPI:
c               reduce the logical flag step_cut_flags to all processors
c               so that if any element requested a step cutback, then
c               all processors are aware of it.  Also reduce back the
c               (scalar) internal energy to the root procaessor.
c
      call wmpi_redlog( material_cut_step )
      call wmpi_reduce_vec( internal_energy, 1 )
      call wmpi_reduce_vec( plastic_work, 1 )
c
c             scatter the element internal forces (stored in blocks) into
c             the global vector (unless we are cutting step)
c             initialize the global internal force vector.
c
      ifv(1:nodof) = zero    ! for this rank
      sum_ifv      = zero    ! for this rank
      num_term_ifv = 0       ! for this rank
      if( material_cut_step ) then
        call allocate_ifv( 2 )
        return
      end if
c
c             assemble internal force vectors for elements into
c             structure vector, and store ifv contributions from
c             killable elements. if block has all killed elements
c             then skip then entire block.  Each processor only computes
c             the blocks which they specifically own. Use threads to
c             process blocks in parallel on this domain.
c
c             timing studies show that use of ATOMIC UPDATE in
c             addifv is 10 x faster than allocating an ifv vector per
c             thread, then gathering vectors. even then, this whole
c             gather process is quite fast
c
c             we use local vectors to store sums, counts for each thread
c             to elminate atomic update on those variable (may still be
c             faster to use atomic).
c
      sum_ifv_threads(1:num_threads) = zero
      num_term_ifv_threads(1:num_threads) = 0
c
      call omp_set_dynamic( .false. )
      if( local_debug ) start_time = omp_get_wtime()
c
c$OMP PARALLEL DO PRIVATE( blk, now_thread, felem, num_enodes,
c$OMP&                     num_enode_dof, totdof, span )
c
      do blk = 1, nelblk
         if( elblks(2,blk) .ne. myid ) cycle
         if( growth_by_kill ) then
           if( dam_blk_killed(blk) ) cycle
         end if
         now_thread    = omp_get_thread_num() + 1
         felem         = elblks(1,blk)
         num_enodes    = iprops(2,felem)
         num_enode_dof = iprops(4,felem)
         totdof        = num_enodes * num_enode_dof
         span          = elblks(0,blk)
         if( local_debug ) write(out,9300) myid, blk, span, felem
         if( allocated( dam_ifv ) ) then
           call addifv( span, edest_blocks(blk)%ptr(1,1), totdof,
     &                ifv(1),
     &                iprops, felem, sum_ifv_threads(now_thread),
     &                num_term_ifv_threads(now_thread),
     &                einfvec_blocks(blk)%ptr(1,1), dam_ifv,
     &                dam_state, out )
         else
           call addifv( span, edest_blocks(blk)%ptr(1,1), totdof,
     &                ifv(1),
     &                iprops, felem, sum_ifv_threads(now_thread),
     &                num_term_ifv_threads(now_thread),
     &                einfvec_blocks(blk)%ptr(1,1), idummy1,
     &                idummy2, out )
         end if
c
      end do
c
c$OMP END PARALLEL DO
c
      if( local_debug ) then
         end_time = omp_get_wtime()
         write(out,*) '>> atomic scatter ifv: ', end_time - start_time
      end if
c
c             reduction of scalars sum_ifv and
c             num_term_ifv for all threads used to process blocks on
c             this rank.
c
      do j = 1, num_threads
        sum_ifv      = sum_ifv + sum_ifv_threads(j)
        num_term_ifv = num_term_ifv + num_term_ifv_threads(j)
      end do

      if( local_debug_sums ) then
         write (out,9000) myid, num_term_ifv, norm2( ifv ),
     &                    sum_stress_norm2s, sum_strain_norm2s
      end if
c
      if( iter .eq. 0 .and. local_debug ) then
       write(out,*) '... drive sig-eps ifv for iter 0'
       do kkk = 1, min(100,nodof)
         write(out,fmt="(i8,f15.8)")  kkk, ifv(kkk)
       end do
      end if
c
c            MPI:
c              sum the worker contributions to the ifv terms,
c              the number of ifv terms, and the whole ifv vector back
c              on the root processor.
c
      call wmpi_reduce_vec( sum_ifv, 1 )
      call wmpi_redint( num_term_ifv )
      call wmpi_reduce_vec( ifv(1), nodof )
      if( local_debug_sums ) then
         if( root_processor ) write(out,9010) norm2( ifv )
      end if

c
c            deallocate space for blocks of element internal force
c            vectors.
c
      call allocate_ifv( 2 )
c
c            all done
c
      if( slave_processor ) return
c
c            modify ifv by beta_factor for thickness other than 1.
c            only used in 2d analysis
c
      if( beta_fact .ne. one ) ifv(1:nodof) =
     &                            beta_fact * ifv(1:nodof)
c
c            set flag indicating that the internal force
c            vector has been calculated.
c
      ifvcmp = .true.
      call thyme( 6, 2 )
c
      return
c
c
8900  format(1x,'drive_eps_sig. rank=5, blk norms strains, stresses:',
     &       i10,e14.6,2x,e14.6 )
 9000 format(1x,'drive_eps_sig. rank, terms_ifv, norm ifv,',
     & ' sum_sig_norm sum_eps_norm:',
     &       i4,i12,e14.6,2x,e14.6,2x,e14.6)
 9010 format(/,1x,
     &  'drive_eps_sig. norm2 ifv after reduce to root: ',e14.6)
 9300 format(i5,':>>> ready to call addifv:',
     &     /,10x,'blk, span, felem                 :',3i10)
 9400 format(5x,'>>> dump of internal force vector:')
 9410 format(i10,2x,f20.8)
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine do_nleps_block               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/21/2017 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine do_nleps_block( blk, iter, step, material_cut_step,
     &                           block_energy, block_plastic_work,
     &                           block_stress_sums, block_strain_sums,
     &                           local_debug_sums )
      use global_data ! old common.main
c
      use elem_extinct_data, only : dam_blk_killed, dam_state
      use elem_block_data,   only : einfvec_blocks, cdest_blocks,
     &                              edest_blocks, element_vol_blocks,
     &                              nonlocal_flags, nonlocal_data_n1
      use main_data,         only : trn, incid, incmap,
     &                              temperatures_ref,
     &                              fgm_node_values_defined,
     &                              cohesive_ele_types, bar_types,
     &                              linear_displ_ele_types,
     &                              adjust_constants_ele_types,
     &                              axisymm_ele_types, link_types,
     &                              nonlocal_analysis, imatprp,
     &                              initial_stresses
      use segmental_curves, only : max_seg_points, max_seg_curves
      use damage_data, only : dam_ptr, growth_by_kill
c
      implicit none
      include 'include_sig_up'
c
      integer :: blk, iter, step
      double precision :: block_energy, block_plastic_work,
     &                    block_stress_sums(1), block_strain_sums(1)
      logical :: material_cut_step, local_debug_sums
c
      integer :: elem, ii, jj, span, felem, matnum, mat_type,
     &           num_enodes, num_enode_dof, totdof, num_int_points,
     &           elem_type, int_order, cohes_type, surface
      double precision, parameter :: zero=0.0d0
      integer, external :: omp_get_thread_num
      logical :: geo_non_flg, bbar_flg, tet_elem, tri_elem,
     &           axisymm_elem, cohesive_elem, used_flg
      logical, parameter :: local_debug = .false.
c
      span           = elblks(0,blk)
      felem          = elblks(1,blk)
      matnum         = iprops(38,felem)
      mat_type       = iprops(25,felem)
      num_enodes     = iprops(2,felem)
      num_enode_dof  = iprops(4,felem)
      totdof         = num_enodes * num_enode_dof
      num_int_points = iprops(6,felem)
      geo_non_flg    = lprops(18,felem)
      elem_type      = iprops(1,felem)
      int_order      = iprops(5,felem)
      bbar_flg       = lprops(19,felem)
      cohes_type     = iprops(27,felem)
      surface        = iprops(26,felem)
c
      call chk_killed_blk( blk, local_work%killed_status_vec,
     &                      local_work%block_killed )
      if( local_work%block_killed ) then
        einfvec_blocks(blk)%ptr(1:span,1:totdof) = zero
        return
      end if
c
      local_work%now_thread         = omp_get_thread_num()
      local_work%span               = span
      local_work%felem              = felem
      local_work%blk                = blk
      local_work%iout               = out
      local_work%num_threads        = num_threads
      local_work%matnum             = matnum
      local_work%mat_type           = mat_type
      local_work%num_enodes         = num_enodes
      local_work%num_enode_dof      = num_enode_dof
      local_work%totdof             = totdof
      local_work%num_int_points     = num_int_points
      local_work%geo_non_flg        = geo_non_flg
      local_work%elem_type          = elem_type
      local_work%int_order          = int_order
      local_work%bbar_flg           = bbar_flg
      local_work%iter               = iter
      local_work%step               = step
      local_work%material_cut_step  = .false.
      local_work%block_energy       = zero
      local_work%block_plastic_work = zero
      local_work%eps_bbar           = eps_bbar
      local_work%dt                 = dt
      local_work%time_n             = total_model_time
      local_work%beta_fact          = beta_fact
      local_work%signal_flag        = signal_flag
      local_work%adaptive_flag      = adaptive_flag
      local_work%temperatures       = temperatures
      local_work%temperatures_ref   = temperatures_ref
      local_work%step_scale_fact    = scaling_factor
      local_work%cohes_type         = cohes_type
      local_work%surface            = surface
      local_work%fgm_enode_props    = fgm_node_values_defined
      local_work%is_cohes_elem      = cohesive_ele_types(elem_type)
      local_work%is_cohes_nonlocal  = nonlocal_analysis .and.
     &                                local_work%is_cohes_elem
      local_work%is_solid_matl      = .not. local_work%is_cohes_elem
      local_work%is_umat            = mat_type .eq. 8
      local_work%is_crys_pls        = mat_type .eq. 10
      local_work%linear_displ_elem  = linear_displ_ele_types(elem_type)
      local_work%adjust_const_elem  =
     &                             adjust_constants_ele_types(elem_type)
      local_work%is_axisymm_elem    = axisymm_ele_types(elem_type)
      local_work%cep_sym_size = 21
      if( local_work%is_cohes_elem ) local_work%cep_sym_size = 6
      if( local_work%is_umat ) call material_model_info( felem, 0, 3,
     &                                 local_work%umat_stress_type )
      local_work%compute_f_bar =  bbar_flg .and.
     &            ( elem_type  .eq. 2 ) .and.
     &           ( local_work%is_umat .or. local_work%is_crys_pls )
      local_work%compute_f_n = geo_non_flg .and.
     &        ( local_work%is_umat .or. local_work%is_crys_pls )
      tet_elem = elem_type .eq. 6 .or. elem_type .eq. 13
      tri_elem = elem_type .eq. 8
      axisymm_elem = axisymm_ele_types(elem_type)
      local_work%is_bar_elem  = bar_types(elem_type)
      local_work%is_link_elem = link_types(elem_type)
      local_work%process_initial_stresses =
     &     allocated( initial_stresses ) .and. step == 1

c
c             See if we're actually an interface damaged material
c
      local_work%is_inter_dmg = .false.
      if( iprops(42,felem) .ne. -1 ) then
            local_work%is_inter_dmg = .true.
            local_work%inter_mat = iprops(42,felem)
            local_work%macro_sz = imatprp(132, local_work%inter_mat)
            local_work%cp_sz = imatprp(133, local_work%inter_mat)
      end if
c
c             build data structures for elements in this block.
c             this is a gather operation on nodal coordinates,
c             nodal displacements, stresses at time n, current
c             estimate of strain increment over step n->n+1,
c             various material state and history data at step n
c             for the specific material model
c             associated with elements of the block. we might be able
c             to skip processing of this block.
c
c             we do this in 3 steps: (1) allocate all data structures in
c             local_work (local_work is on stack but components are
c             allocatable), (2) data that is stored globally
c             in simple vectors and (3) data that is stored globally in
c             blocked structures. a mixure of arguments are passed
c             vs. data in modules to optimize indexing into blocked
c             data structures.
c
      if( local_debug ) then
            write(out,9100) blk, span, felem, mat_type, num_enodes,
     &                      num_enode_dof, totdof, num_int_points,
     &                      elem_type, int_order, bbar_flg
      end if
c
c             are there solid elements in this block involved in a
c             nonlocal analysis? set flags.
c
      call do_nleps_block_a
c
      call recstr_allocate( local_work )
c
      call dupstr( span, edest_blocks(blk)%ptr(1,1),
     &             incid(incmap(felem)), felem,
     &             num_enodes, num_enode_dof, totdof,
     &             local_work%trn_e_flags,
     &             local_work%trn_e_block,
     &             local_work%trne,
     &             local_work%trnmte,
     &             local_work%ue,
     &             local_work%due, trn )
c
      call dupstr_blocked( blk, span, felem, num_int_points,
     &   num_enodes, num_enode_dof, totdof, mat_type,
     &   geo_non_flg, step, iter, incid(incmap(felem)),
     &   cdest_blocks(blk)%ptr(1,1), local_work%ce_0,
     &   local_work%ce_n, local_work%ce_mid,
     &   local_work%ce_n1, local_work%ue, local_work%due, local_work )
c
c             compute updated strain increment for specified
c             displacement increment. update stresses at all
c             strain points of elements in block.
c
      if( local_debug ) then
         write(out,9200) blk, span, felem, geo_non_flg
      end if
c
      call rknstr( props(1,felem), lprops(1,felem),
     &             iprops(1,felem), local_work )
c
c            get the cut step size flag for the block based
c            on material model computations. set passed in value
c            for the block. just leave now if we need a cut.
c
      material_cut_step = local_work%material_cut_step
      if( material_cut_step ) return
c
c            save the summed internal energy for elements
c            in this block into the block value passed in.
c
      block_energy       = local_work%block_energy
      block_plastic_work = local_work%block_plastic_work
c
c            perform scatter operation to update global
c            data structures at n+1 from element block
c            data structures. blocks being processed in parallel can
c            access these data structures w/o conflict.
c
      if( local_debug ) write(out,9520) blk, span, felem
      call rplstr( span, felem, num_int_points, mat_type, iter,
     &             geo_non_flg, local_work, blk, block_stress_sums,
     &             block_strain_sums, local_debug_sums )
      if( local_debug ) write(out,9522) blk, span, felem
c
c            compute updated internal force vectors for each element
c            in the block.
c
      if( local_debug ) then
        write(out,9500) blk, span, felem
        write(out,*) ">>> size of element_vol_blocks: ",
     &         sizeof(    element_vol_blocks(blk)%ptr )
      end if
      call rknifv( einfvec_blocks(blk)%ptr(1,1),
     &             element_vol_blocks(blk)%ptr(1), props, span,
     &             local_work )
      if( local_debug ) write(out,9505) blk, span, felem
c
c             release all allocated data for block
c
      if( local_debug ) write(out,9510) blk, span, felem
      call recstr_deallocate( local_work )
      if( local_debug ) write(out,9515) blk, span, felem

      if( local_debug ) write(out,9525) blk, span, felem
      return
c
 9100 format(5x,'>>> ready to call dupstr:',
     &     /,10x,'blk, span, felem, mat_model:      ',4i10,
     &     /,10x,'num_enodes, num_enode_dof, totdof:',3i10,
     &     /,10x,'num_int_points:                   ',i10,
     &     /,10x,'elem_type, int_order, bbar_flg lg:',2i10,l10)
 9200 format(5x,'>>> ready to call rknstr:',
     &     /,10x,'blk, span, felem                 :',3i10,
     &     /,10x,'geo_non_flg                      :',l10 )
 9300 format(5x,'>>> ready to call rplstr:',
     &     /,10x,'blk, span, felem                 :',3i10)
 9400 format(/,5x,'>>> nodal displacments at n and increments: '
     &      // )
 9410 format(8x,i9,1x,6f15.6)
 9500 format(5x,'>>> ready to call rknifv. blk, span, felem: ',3i6)
 9505 format(5x,'>>> back from rknifv. blk, span, felem: ',3i6)
 9510 format(5x,'>>> ready to call recstr_deallocate. ',
     &     'blk, span, felem:',3i6)
 9515 format(5x,'>>> back from recstr_deallocate. ',
     &     'blk, span, felem:',3i6)
 9520 format(5x,'>>> ready to call rplstr. blk, span, felem:',3i6)
 9522 format(5x,'>>> back from rplstr. blk, span, felem:',3i6)
 9525 format(5x,'>>> leaving do_nleps_block. blk, span, felem:',3i6)
c
      contains
c     ========
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine do_nleps_block_a                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/30/2015 rhd              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine do_nleps_block_a
      implicit none
c
      integer :: i, elem_num
c
      local_work%block_has_nonlocal_solids = .false.
      if( .not. nonlocal_analysis ) return ! regular local analysis
      if( local_work%is_cohes_elem ) return
c
c              globally we have a nonlocal analysis.
c              does this block have any solid elements
c              that need nonlocal, shared variables -- means
c              they have an allocated nonlocal state vector.
c
      do i = 1, span ! check each element in block
        elem_num = felem + i - 1
        if( .not. nonlocal_flags(elem_num) ) cycle ! local element
        if( allocated( nonlocal_data_n1(elem_num)%state_values ) ) then
           local_work%block_has_nonlocal_solids = .true.
        else ! we have a bad error before here
           write(out,9000) elem_num, blk
           call die_abort
        end if
      end do
c
      return
c
 9000 format('>> FATAL ERROR: do_nleps_block_a'
     &  /,   '  unallocated nonlocal state vector for element: ',i7,
     &  /,   '  in block: ',i7,
     &  /,   '  job terminated' )
c
      end subroutine do_nleps_block_a
      end subroutine do_nleps_block
c     ****************************************************************
c     *                                                              *
c     *                   subroutine recstr_allocate                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified :  9/25,2017 rhd             *
c     *                                                              *
c     *     allocate data structure in local_work for updating       *
c     *     strains-stresses-internal forces.                        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine recstr_allocate( local_work )
      use global_data ! old common.main
      use segmental_curves, only : max_seg_points, max_seg_curves
      use elem_block_data, only: history_blk_list
      implicit none

      include 'include_sig_up'
c
      integer :: local_mt, error, span, blk, ngp, hist_size, nlsize
      double precision :: zero
      data zero / 0.0d00 /
c
      local_mt = local_work%mat_type
c
      allocate(
     &   local_work%ce_0(mxvl,mxecor),
     &   local_work%ce_n(mxvl,mxecor),
     &   local_work%ce_mid(mxvl,mxecor),
     &   local_work%ce_n1(mxvl,mxecor), stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 1
         call die_abort
      end if
c
      allocate( local_work%trnmte(mxvl,mxedof,mxndof) )
c
      allocate(
     1 local_work%det_j(mxvl,mxgp),
     2 local_work%det_j_mid(mxvl,mxgp),
     3 local_work%nxi(mxndel,mxgp),
     4 local_work%neta(mxndel,mxgp),
     5 local_work%nzeta(mxndel,mxgp),
     6 local_work%gama(mxvl,3,3,mxgp),
     7 local_work%gama_mid(mxvl,3,3,mxgp), stat=error  )
      if( error .ne. 0 ) then
         write(out,9000) 2
         call die_abort
      end if
!DIR$ VECTOR ALIGNED
       local_work%det_j = zero
!DIR$ VECTOR ALIGNED
       local_work%det_j_mid = zero
!DIR$ VECTOR ALIGNED
       local_work%gama = zero
!DIR$ VECTOR ALIGNED
       local_work%gama_mid = zero
!DIR$ VECTOR ALIGNED
       local_work%neta = zero
!DIR$ VECTOR ALIGNED
       local_work%nxi = zero
!DIR$ VECTOR ALIGNED
       local_work%nzeta = zero
c
      if( local_work%geo_non_flg ) then
        allocate( local_work%fn(mxvl,3,3),
     &           local_work%fn1(mxvl,3,3),
     &           local_work%dfn1(mxvl) )
!DIR$ VECTOR ALIGNED
        local_work%fn  = zero
!DIR$ VECTOR ALIGNED
        local_work%fn1 = zero
!DIR$ VECTOR ALIGNED
        local_work%dfn1 = zero
      end if


      allocate(
     &  local_work%vol_block(mxvl,8,3),
     &  local_work%volume_block(mxvl),
     &  local_work%volume_block_0(mxvl),
     &  local_work%volume_block_n(mxvl),
     &  local_work%volume_block_n1(mxvl),
     &  local_work%jac(mxvl,3,3),
     &  local_work%b(mxvl,mxedof,nstr),
     &  local_work%ue(mxvl,mxedof),
     &  local_work%due(mxvl,mxedof),
     &  local_work%uenh(mxvl,mxedof), stat=error  )
      if( error .ne. 0 ) then
         write(out,9000) 3
         call die_abort
      end if

!DIR$ VECTOR ALIGNED
      local_work%b = zero
c
      allocate( local_work%uen1(mxvl,mxedof),
     &  local_work%urcs_blk_n(mxvl,nstrs,mxgp),
     &  local_work%urcs_blk_n1(mxvl,nstrs,mxgp),
     &  local_work%initial_stresses(6,mxvl),
     3  local_work%rot_blk_n1(mxvl,9,mxgp),
     4  local_work%rtse(mxvl,nstr,mxgp), stat=error  )
      if( error .ne. 0 ) then
         write(out,9000) 4
         call die_abort
      end if
!DIR$ VECTOR ALIGNED
      local_work%det_j       = zero
!DIR$ VECTOR ALIGNED
      local_work%rot_blk_n1  = zero
!DIR$ VECTOR ALIGNED
      local_work%urcs_blk_n1 = zero
!DIR$ VECTOR ALIGNED
      local_work%urcs_blk_n  = zero
!DIR$ VECTOR ALIGNED
      local_work%initial_stresses = zero

c
      allocate( local_work%ddtse(mxvl,nstr,mxgp),
     1   local_work%strain_n(mxvl,nstr,mxgp),
     2   local_work%dtemps_node_blk(mxvl,mxndel),
     3   local_work%temps_ref_node_blk(mxvl,mxndel),
     4   local_work%temps_node_blk(mxvl,mxndel),
     5   local_work%temps_node_ref_blk(mxvl,mxndel),
     6   local_work%nu_vec(mxvl),
     7   local_work%beta_vec(mxvl),
     8   local_work%h_vec(mxvl),
     9   local_work%tan_e_vec(mxvl),
     a   local_work%e_vec(mxvl), stat=error  )
      if( error .ne. 0 ) then
         write(out,9000) 5
         call die_abort
      end if
c
      allocate( local_work%sigyld_vec(mxvl),
     1   local_work%alpha_vec(mxvl,6),
     2   local_work%e_vec_n(mxvl),
     3   local_work%nu_vec_n(mxvl),
     4   local_work%gp_sig_0_vec(mxvl),
     5   local_work%gp_sig_0_vec_n(mxvl),
     6   local_work%gp_h_u_vec(mxvl),
     7   local_work%gp_h_u_vec_n(mxvl),
     8   local_work%gp_beta_u_vec(mxvl),
     9   local_work%gp_beta_u_vec_n(mxvl), stat=error  )
      if( error .ne. 0 ) then
         write(out,9000) 6
         call die_abort
      end if
c
      allocate( local_work%gp_delta_u_vec(mxvl),
     1   local_work%gp_delta_u_vec_n(mxvl),
     2   local_work%alpha_vec_n(mxvl,6),
     3   local_work%h_vec_n(mxvl),
     4   local_work%n_power_vec(mxvl) )
c
      local_mt = local_work%mat_type
c
      if( local_mt .eq. 3 ) then
        allocate(
     1   local_work%f0_vec(mxvl),
     2   local_work%eps_ref_vec(mxvl),
     3   local_work%m_power_vec(mxvl),
     4   local_work%q1_vec(mxvl),
     5   local_work%q2_vec(mxvl),
     6   local_work%q3_vec(mxvl),
     7   local_work%nuc_s_n_vec(mxvl),
     8   local_work%nuc_e_n_vec(mxvl),
     9   local_work%nuc_f_n_vec(mxvl), stat=error  )
        if( error .ne. 0 ) then
         write(out,9000) 7
         call die_abort
        end if
      end if
c
      allocate( local_work%eps_curve(max_seg_points),
     1    local_work%shape(mxndel,mxgp),
     5    local_work%enode_mat_props(mxndel,mxvl,mxndpr),
     6    local_work%fgm_flags(mxvl,mxndpr), stat=error  )
      if( error .ne. 0 ) then
         write(out,9000) 9
         call die_abort
      end if
c
      if( local_mt .eq. 5 ) allocate( local_work%mm05_props(mxvl,10) )
      if( local_mt .eq. 6 ) allocate( local_work%mm06_props(mxvl,10) )
      if( local_mt .eq. 7 ) allocate( local_work%mm07_props(mxvl,10) )
      if( local_work%is_umat ) then
        allocate(local_work%umat_props(mxvl,50),
     1           local_work%characteristic_length(mxvl) )
      end if

c
      allocate( local_work%trne(mxvl,mxndel), stat=error  )
      if( error .ne. 0 ) then
         write(out,9000) 11
         call die_abort
      end if
      local_work%trne = .false.
c
      if( local_mt .eq. 10 .or. local_mt .eq. 11 ) then
        allocate( local_work%debug_flag(mxvl),
     1    local_work%local_tol(mxvl),
     2    local_work%ncrystals(mxvl),
     3    local_work%angle_type(mxvl),
     4    local_work%angle_convention(mxvl),
     5    local_work%c_props(mxvl,max_crystals),
     6    local_work%nstacks(mxvl),
     7    local_work%nper(mxvl))
      end if
c
      span                         = local_work%span
      blk                          = local_work%blk
      ngp                          = local_work%num_int_points
      hist_size                    = history_blk_list(blk)
      local_work%hist_size_for_blk = hist_size
c
      allocate( local_work%elem_hist1(span,hist_size,ngp),
     &          local_work%elem_hist(span,hist_size,ngp), stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 12
         call die_abort
      end if
c
      if( local_work%is_cohes_nonlocal ) then
         nlsize = nonlocal_shared_state_size
         allocate( local_work%top_surf_solid_stresses_n(mxvl,nstrs),
     &      local_work%bott_surf_solid_stresses_n(mxvl,nstrs),
     &      local_work%top_surf_solid_eps_n(mxvl,nstr),
     &      local_work%bott_surf_solid_eps_n(mxvl,nstr),
     &      local_work%top_surf_solid_elements(mxvl),
     &      local_work%bott_surf_solid_elements(mxvl),
     &      local_work%nonlocal_stvals_bott_n(mxvl,nlsize),
     &      local_work%nonlocal_stvals_top_n(mxvl,nlsize),
     &      stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 13
           call die_abort
         end if
         allocate( local_work%top_solid_matl(mxvl),
     &      local_work%bott_solid_matl(mxvl), stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 14
           call die_abort
         end if
      end if
c
c             always allocate cohes_rot_block even if not
c             a block of cohesive elements. the
c             array is passed in lots of places
c             that process both solid and interface elements.
c             So it needs to exist !

      allocate( local_work%cohes_rot_block(mxvl,3,3), stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 15
         call die_abort
      end if
c
      if( local_work%is_cohes_elem ) then
         allocate( local_work%cohes_temp_ref(mxvl),
     1      local_work%cohes_dtemp(mxvl),
     2      local_work%cohes_temp_n(mxvl),
     3      local_work%intf_prp_block(mxvl,max_interface_props),
     4      stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 16
           call die_abort
         end if
      end if
c
c             when nonlocal analysis in effect, allocate space
c             for block nonlocal state values. we need to store span x
c             num local values (fixed global in param_def). Space
c             will be passed to material model code to insert values
c             for elements in the block at the integration point
c             being processed.
c             allocate dummy 1x1 for std. local analysis to
c             simplify calls later.
c
      if( local_work%block_has_nonlocal_solids ) then
         allocate( local_work%nonlocal_state_blk(mxvl,
     &             nonlocal_shared_state_size), stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 17
           call die_abort
         end if
      else
         allocate( local_work%nonlocal_state_blk(1,1), stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 18
           call die_abort
         end if
      end if
c
      allocate( local_work%weights(mxgp), stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 19
         call die_abort
      end if
c
      allocate( local_work%sv(3), local_work%lv(3),
     &          local_work%tv(3), stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 20
         call die_abort
      end if
c
      return
c
 9000 format('>> FATAL ERROR: recstr_allocate'
     &  /,   '                failure status= ',i5,
     &  /,   '                job terminated' )
c
      end

c     ****************************************************************
c     *                                                              *
c     *                   subroutine recstr_deallocate               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified :  3/15/2017 rhd             *
c     *                                                              *
c     *     release data structure in local_work for updating        *
c     *     strains-stresses-internal forces.                        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine recstr_deallocate( local_work )
      use global_data ! old common.main
      implicit none
      include 'include_sig_up'
c
      integer :: local_mt, error
      logical :: local_debug
c
      local_debug = .false.
      if( local_debug ) write(out,*) "..recstr_dell @ 1"
      local_mt = local_work%mat_type
c

      deallocate(
     & local_work%ce_0,
     & local_work%ce_n,
     & local_work%ce_mid,
     & local_work%ce_n1, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 1
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 2"
      deallocate( local_work%trnmte, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 100
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 3"
      deallocate(
     1 local_work%det_j,
     2 local_work%det_j_mid,
     3 local_work%nxi,
     4 local_work%neta,
     5 local_work%nzeta,
     6 local_work%gama,
     7 local_work%gama_mid, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 2
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 4"
      if( local_work%geo_non_flg ) then
        deallocate( local_work%fn,
     &              local_work%fn1,
     &              local_work%dfn1, stat=error )
        if( error .ne. 0 ) then
         write(out,9000) 102
         call die_abort
        end if
      end if

      if( local_debug ) write(out,*) "..recstr_dell @ 5"
      deallocate(
     &  local_work%vol_block,
     &  local_work%volume_block, local_work%volume_block_0,
     &  local_work%volume_block_n, local_work%volume_block_n1,
     &  local_work%jac,
     &  local_work%b,
     &  local_work%ue,
     &  local_work%due,
     &  local_work%uenh,
     &  local_work%uen1,
     &  local_work%urcs_blk_n,
     &  local_work%urcs_blk_n1,
     &  local_work%initial_stresses,
     &  local_work%rot_blk_n1,
     &  local_work%rtse, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 3
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 6"
      deallocate( local_work%ddtse,
     1   local_work%strain_n,
     2   local_work%dtemps_node_blk,
     3   local_work%temps_ref_node_blk,
     4   local_work%temps_node_blk,
     5   local_work%temps_node_ref_blk,
     6   local_work%nu_vec,
     7   local_work%beta_vec,
     8   local_work%h_vec,
     9   local_work%tan_e_vec,
     a   local_work%e_vec, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 4
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 7"
      deallocate( local_work%sigyld_vec,
     1   local_work%alpha_vec,
     2   local_work%e_vec_n,
     3   local_work%nu_vec_n,
     4   local_work%gp_sig_0_vec,
     5   local_work%gp_sig_0_vec_n,
     6   local_work%gp_h_u_vec,
     7   local_work%gp_h_u_vec_n,
     8   local_work%gp_beta_u_vec,
     9   local_work%gp_beta_u_vec_n, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 5
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 8"
      deallocate( local_work%gp_delta_u_vec,
     1   local_work%gp_delta_u_vec_n,
     2   local_work%alpha_vec_n,
     3   local_work%h_vec_n,
     4   local_work%n_power_vec, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 104
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 9"
      if( local_mt .eq. 3 ) then
        deallocate(
     1   local_work%f0_vec,
     2   local_work%eps_ref_vec,
     3   local_work%m_power_vec,
     4   local_work%q1_vec,
     5   local_work%q2_vec,
     6   local_work%q3_vec,
     7   local_work%nuc_s_n_vec,
     8   local_work%nuc_e_n_vec,
     9   local_work%nuc_f_n_vec, stat=error )
        if( error .ne. 0 ) then
         write(out,9000) 6
         call die_abort
        end if
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 10"
      deallocate( local_work%eps_curve,
     1    local_work%shape,
     5    local_work%enode_mat_props,
     6    local_work%fgm_flags, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 7
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 11"
      if( local_mt .eq. 5 ) then
        deallocate( local_work%mm05_props, stat=error )
        if( error .ne. 0 ) then
         write(out,9000) 106
         call die_abort
        end if
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 12"
      if( local_mt .eq. 6 ) then
        deallocate( local_work%mm06_props, stat=error )
        if( error .ne. 0 ) then
         write(out,9000) 108
         call die_abort
        end if
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 13"
      if( local_mt .eq. 7 ) then
        deallocate( local_work%mm07_props, stat=error )
        if( error .ne. 0 ) then
         write(out,9000) 110
         call die_abort
        end if
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 14"
      if( local_work%is_umat ) then
        deallocate( local_work%umat_props,
     &              local_work%characteristic_length, stat=error )
        if( error .ne. 0 ) then
         write(out,9000) 112
         call die_abort
        end if
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 15"
      deallocate( local_work%trne, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 9
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 16"
      deallocate( local_work%elem_hist, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 10
         call die_abort
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 16.5"
      deallocate( local_work%elem_hist1, stat=error )
      if( error .ne. 0 ) then
         write(out,9000) 10
         call die_abort
      end if

c
      if( local_debug ) write(out,*) "..recstr_dell @ 17"
      if( local_mt .eq. 10 .or. local_mt .eq. 11 ) then
        deallocate( local_work%debug_flag,
     1    local_work%local_tol,
     2    local_work%ncrystals,
     3    local_work%angle_type,
     4    local_work%angle_convention,
     5    local_work%c_props,
     6    local_work%nstacks,
     7    local_work%nper, stat=error)
        if( error .ne. 0 ) then
         write(out,9000) 114
         call die_abort
        end if
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 18"
      if( local_work%is_cohes_nonlocal ) then
         deallocate( local_work%top_surf_solid_stresses_n,
     1      local_work%bott_surf_solid_stresses_n,
     2      local_work%top_surf_solid_eps_n,
     3      local_work%bott_surf_solid_eps_n,
     4      local_work%top_surf_solid_elements,
     5      local_work%bott_surf_solid_elements,
     6      local_work%nonlocal_stvals_bott_n,
     7      local_work%nonlocal_stvals_top_n,
     8      stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 12
           call die_abort
         end if
         if( local_debug ) write(out,*) "..recstr_dell @ 19"
         deallocate( local_work%top_solid_matl,
     &      local_work%bott_solid_matl, stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 13
           call die_abort
         end if
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 20"
      deallocate( local_work%cohes_rot_block, stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 14
           call die_abort
      end if
      if( local_work%is_cohes_elem ) then
         deallocate( local_work%cohes_temp_ref,
     1      local_work%cohes_dtemp,
     2      local_work%cohes_temp_n,
     3      local_work%intf_prp_block, stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 14
           call die_abort
         end if
      end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 21"
      deallocate( local_work%nonlocal_state_blk, stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 15
           call die_abort
         end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 22"
      deallocate( local_work%weights, stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 16
           call die_abort
         end if
c
      if( local_debug ) write(out,*) "..recstr_dell @ 23"
      deallocate( local_work%sv, local_work%lv, local_work%tv,
     &            stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 17
           call die_abort
         end if
      return
c
 9000 format('>> FATAL ERROR: recstr_deallocate'
     &  /,   '                failure status= ',i5,
     &  /,   '                job terminated' )
c
      end


c     ****************************************************************
c     *                                                              *
c     *                   subroutine dupstr_blocked                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/9/2016 rhd               *
c     *                                                              *
c     *     creates a separate copy of element                       *
c     *     blocked data necessary for global stress vector recovery *
c     *     each element in  block                                   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine dupstr_blocked(
     & blk, span, felem, ngp, nnode, ndof, totdof, mat_type,
     & geonl, step, iter, belinc, bcdst, ce_0, ce_n, ce_mid, ce_n1,
     & ue, due, local_work )
      use global_data ! old common.main
c
      use elem_block_data, only:  history_blocks,
     &                            eps_n_blocks, urcs_n_blocks,
     &                            history_blk_list
      use main_data,       only:  dtemp_nodes, dtemp_elems,
     &                            temper_nodes, temper_elems,
     &                            temper_nodes_ref, temperatures_ref,
     &                            fgm_node_values,
     &                            fgm_node_values_defined,
     &                            fgm_node_values_cols, matprp,
     &                            lmtprp
c
      use segmental_curves, only : max_seg_points, max_seg_curves
c
      implicit none
      include 'include_sig_up'
c
c           parameter declarations
c
      integer :: blk, span, felem, ngp, nnode, ndof, totdof, mat_type,
     &           step, iter, hist_size
      logical :: geonl
      integer :: belinc(nnode,*), bcdst(totdof,*)
      double precision ::
     & ce_0(mxvl,*), ce_mid(mxvl,*), ue(mxvl,*), due(mxvl,*),
     & ce_n(mxvl,*), ce_n1(mxvl,*)
c
c           local declarations
c
      integer :: elem_type, surf, k, j, i, matl_no
      logical ::local_debug, update, update_coords, middle_surface
      double precision ::
     &   half, zero, one, mag, mags(3), djcoh(mxvl)
      data local_debug, half, zero, one
     &  / .false., 0.5d00, 0.0d00, 1.0d00 /

c
      if( local_debug ) write(out,9100)
c
      elem_type      = local_work%elem_type
      surf           = local_work%surface
      middle_surface = surf .eq. 2
      djcoh(1:span)  = zero
c

c           pull coordinates at t=0 from global input vector.
c
      k = 1
      do j = 1, nnode
!DIR$ IVDEP
         do i = 1, span
            ce_0(i,k)   = c(bcdst(k,i))
            ce_0(i,k+1) = c(bcdst(k+1,i))
            ce_0(i,k+2) = c(bcdst(k+2,i))
         end do
         k = k + 3
      end do
c
c           for geonl, create a set of nodal coordinates
c           at mid-step and at the end of step. for step 1,
c           iter=0 we don't update for imposed displacements
c           since we are using linear stiffness.
c           if any elements are killed, their displacements
c           were set to zero above so we are just
c           setting the initial node coords. for iter=0 in other
c           steps, just make the mid and n1 coordinates the
c           coordinates at start of step. iter=0 means we are
c           computing equiv nodal forces for imposed/extrapolated
c           displacements to start a step.
c
      update = .true.
      if( step .eq. 1 .and. iter .eq. 0 ) update = .false.
      update_coords = geonl .and. update
c
      if( update_coords ) then
       do  j = 1, totdof
!DIR$ IVDEP
          do i = 1, span
            ce_n(i,j)   = ce_0(i,j) + ue(i,j)
            ce_mid(i,j) = ce_0(i,j) + ue(i,j) + half*due(i,j)
            ce_n1(i,j)  = ce_0(i,j) + ue(i,j) + due(i,j)
          end do
       end do
c
         if( local_work%is_cohes_elem ) then ! geonl interface elem
          if( middle_surface ) call chk_cohes_penetrate( span, mxvl,
     &               felem, mxndel, nnode, elem_type, ce_n1, djcoh )
          call cohes_ref_surface( span, mxvl, mxecor,
     &                            local_work%surface, nnode,
     &                            totdof, ce_0, ce_n1, djcoh )
          end if
      end if   !  update_coords
c
      if( .not. update_coords ) then
          do  j = 1, totdof
!DIR$ IVDEP
            do i = 1, span
              ce_n(i,j)   = ce_0(i,j)
              ce_mid(i,j) = ce_0(i,j)
              ce_n1(i,j)  = ce_0(i,j)
            end do
          end do
      end if
c
c           gather nodal and element temperature change over load step
c           (if they are defined). construct a set of incremental nodal
c           temperatures for each element in block.
c
      if( temperatures ) then
        if( local_debug )  write(out,9610)
        call gadtemps( dtemp_nodes, dtemp_elems(felem), belinc,
     &                 nnode, span, felem, local_work%dtemps_node_blk,
     &                 mxvl )
      else
        local_work%dtemps_node_blk(1:span,1:nnode) = zero
      end if
c
c           gather reference temperatures for element nodes from the
c           global vector of reference values at(if they are defined).
c           construct a set of reference nodal temperatures for each
c           element in block.
c
      if( temperatures_ref ) then
        if( local_debug )  write(out,9610)
        call gartemps( temper_nodes_ref, belinc, nnode, span,
     &                 felem, local_work%temps_ref_node_blk, mxvl )
      else
        local_work%temps_ref_node_blk(1:span,1:nnode) = zero
      end if
c
c           build nodal temperatures for elements in the block
c           at end of step (includes both imposed nodal and element
c           temperatures)
c
      if( local_debug )  write(out,9620)
      call gatemps( temper_nodes, temper_elems(felem), belinc,
     &              nnode, span, felem, local_work%temps_node_blk,
     &              mxvl, local_work%dtemps_node_blk,
     &              local_work%temps_node_to_process )

c
c           if the model has fgm properties at the model nodes, build a
c           table of values for nodes of elements in the block
c
      if( fgm_node_values_defined ) then
        do j = 1,  fgm_node_values_cols
          do i = 1, nnode
!DIR$ IVDEP
            do k = 1, span
              local_work%enode_mat_props(i,k,j) =
     &                     fgm_node_values(belinc(i,k),j)
            end do
          end do
        end do
      end if
c
c
c           gather material specific data for elements
c           in the block. we split operations based on
c           the material model associated with block
c
      if( local_debug )  write(out,9600)
c
c           gather element data at n from global blocks:
c            a) stresses -  unrotated cauchy stresses for geonl
c            b) strain_n -  strains at n
c            b) ddtse -     strains at start of step n, subsequently
c                           updated to strains at n+1 during strain-
c                           stress updating
c            c) elem_hist - integration point history data for material
c                           models
c           History data:
c            o The global blocks are sized(hist_size,ngp,span)
c            o The local block is sized (span,hist_size,ngp).
c              This makes it possible to pass a 2-D array slice for
c              all elements of the block for a single integration pt.
c
      hist_size = local_work%hist_size_for_blk
c
      call dptstf_copy_history(
     &  local_work%elem_hist(1,1,1), history_blocks(blk)%ptr(1),
     &            ngp, hist_size, span )
c
      call recstr_gastr( local_work%ddtse, eps_n_blocks(blk)%ptr(1),
     &                   ngp, nstr, span )
      call recstr_gastr( local_work%strain_n, eps_n_blocks(blk)%ptr(1),
     &                   ngp, nstr, span )
c
      call recstr_gastr( local_work%urcs_blk_n,
     &                   urcs_n_blocks(blk)%ptr(1),
     &                   ngp, nstrs, span )
c
c
      select case ( mat_type )
c
      case ( 1 )
c
c           vectorized mises plasticty model.
c
!DIR$ IVDEP
        do i = 1, span
           matl_no = iprops(38,felem+i-1)
           local_work%tan_e_vec(i) = matprp(4,matl_no)
        end do
c
      case ( 2 )
c
c           nonlinear elastic material model (deformation plasticity).
c
         if ( local_debug ) write(out,9950)
c
c
      case ( 3 )
c
c           general mises/gurson model.
c
        if ( local_debug ) write(out,9950)
!DIR$ IVDEP
        do i = 1, span
           matl_no = iprops(38,felem+i-1)
           local_work%tan_e_vec(i) = matprp(4,matl_no)
        end do
c
      case ( 4 )
c
c           cohesive zone models
c
         if ( local_debug ) write(out,9950)
         if( local_work%is_cohes_nonlocal ) then
           call recstr_build_cohes_nonlocal( local_work, iprops )
         end if
c
      case ( 5 )
c
c           advanced cyclic plasticity model
c
        if ( local_debug ) write(out,9950)
c
      case ( 6 )
c
c           creep model
c
        if ( local_debug ) write(out,9960)

      case ( 7 )
c
c           mises model + hydrogen
c
        if ( local_debug ) write(out,9970)

      case ( 8 )
c
c           Abaqus compatible UMAT
c
        if ( local_debug ) write(out,9980)
c
      case ( 10 )
c
c           CP model
c
        if ( local_debug) write(out,9990)
c
c
      case default
          write(out,*) '>>> invalid material model number'
          write(out,*) '    in dupstr_blocked'
          call die_abort
          stop
c
      end select
c
      if ( local_debug ) write(out,9150)
c
      return
c
 9100 format(8x,'>> entered dupstr_blocked...' )
 9150 format(8x,'>> leaving dupstr_blocked...' )
 9600 format(12x,'>> gather element stresses/strains at step n...' )
 9610 format(12x,'>> gather element incremental temperatures...' )
 9620 format(12x,'>> gather element total temperatures...' )
 9800 format(12x,'>> gather material data for type 1...' )
 9900 format(15x,'>> gather plast. parms, back stress, state...' )
 9950 format(12x,'>> gather data for model type 5...' )
 9960 format(12x,'>> gather data for model type creep...' )
 9970 format(12x,'>> gather data for model type 7...' )
 9980 format(12x,'>> gather data for model type 8...' )
 9990 format(12x,'>> gather data for model type 10...' )
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine recstr_gastr                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/19/2013 rhd             *
c     *                                                              *
c     *     gathers element stresses from the global stress data     *
c     *     structure to local block array                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine recstr_gastr( mlocal, mglobal, ngp, nprm, span )
      implicit none
      include 'param_def'
c
c               parameter declarations
c
      integer :: ngp, nprm, span
      double precision ::
     & mlocal(mxvl,nprm,*), mglobal(nprm,ngp,*)
c
      integer :: k, j, i
c
c           unroll inner loop for most common number of integration
c           points (ngp).
c
c
      if( ngp .ne. 8 ) then
        do k = 1, ngp
         do  j = 1, nprm
!DIR$ VECTOR ALIGNED
            do  i = 1, span
               mlocal(i,j,k) = mglobal(j,k,i)
            end do
         end do
        end do
        return
      end if
c
c                number of integration points = 8, unroll.
c
      do  j = 1, nprm
!DIR$ VECTOR ALIGNED
        do  i = 1, span
            mlocal(i,j,1) = mglobal(j,1,i)
            mlocal(i,j,2) = mglobal(j,2,i)
            mlocal(i,j,3) = mglobal(j,3,i)
            mlocal(i,j,4) = mglobal(j,4,i)
            mlocal(i,j,5) = mglobal(j,5,i)
            mlocal(i,j,6) = mglobal(j,6,i)
            mlocal(i,j,7) = mglobal(j,7,i)
            mlocal(i,j,8) = mglobal(j,8,i)
        end do
      end do
c
      return
      end


c     ****************************************************************
c     *                                                              *
c     *                      subroutine allocate_ifv                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/1/2013 rhd               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine allocate_ifv( action  )
      use global_data ! old common.main
      use elem_block_data, only:  einfvec_blocks
      implicit none
c
      integer :: action
c
      integer :: iok, blk, felem, num_enodes, num_enode_dof, totdof,
     &           span
c
      select case ( action )
      case( 1 )
c
c           allocate blks of elem internal force vectors
c
      allocate ( einfvec_blocks(nelblk), stat=iok )
      if ( iok .ne. 0 ) then
          write(out,9100) iok
          call die_abort
          stop
      end if
c
      do blk = 1, nelblk
c
c           MPI:
c             elblks(2,blk) holds which processor owns block. If we
c             don't own the block, skip allocation.
c           Threads-only:
c             elblks(2,blk) all = 0, allocate all blocks
c
         if( myid .ne. elblks(2,blk) ) cycle
c
         felem         = elblks(1,blk)
         num_enodes    = iprops(2,felem)
         num_enode_dof = iprops(4,felem)
         totdof        = num_enodes * num_enode_dof
         span          = elblks(0,blk)
         allocate( einfvec_blocks(blk)%ptr(span,totdof),stat=iok )
         if( iok .ne. 0 ) then
           write(out,9100) iok
           call die_abort
         end if
      end do
      return
c
      case( 2 )
c
c           deallocate blks of elem internal force vectors
c           see comments about MPI, threads in action type 1
c
      do blk = 1, nelblk
         if( myid .ne. elblks(2,blk) ) cycle
         deallocate( einfvec_blocks(blk)%ptr,stat=iok )
         if ( iok .ne. 0 ) then
           write(out,9200) iok
           call die_abort
         end if
      end do
c
      deallocate( einfvec_blocks, stat=iok )
      if( iok .ne. 0 ) then
          write(out,9200) iok
          call die_abort
      end if
      return
c
      case default
c
        write(out,9100) 1
        call die_abort
c
      end select
c
 9999 return
 9100 format('>> FATAL ERROR: einfvec_allocate, memory deallocate'
     &  /,   '                failure status= ',i5,
     &  /,   '                job terminated' )
 9200 format('>> FATAL ERROR: einfvec_allocate, memory deallocate',
     &  /,   '                failure status= ',i5,
     &  /,   '                job terminated' )
 9300 format('>> FATAL ERROR: allocate_ifv reports error: ',i2,
     &     /,'                job aborted' )

      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine recstr_build_cohes_nonlocal       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 03/1/2013 rhd              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine recstr_build_cohes_nonlocal( local_work, iprops  )

      use elem_block_data, only:  solid_interface_lists
      use main_data, only:  nonlocal_analysis
c
      implicit none
      include 'param_def'
      include 'include_sig_up'
c
      integer :: iprops(mxelpr,*)  ! global element props array
c
c           local declarations.
c
      integer :: span, felem, iout, iter, blk, rel_elem, j, n, step,
     &           elem_top, elem_bott, matl_top, ngp_top, ngp_bott,
     &           matl_bott, top_mat_model, bott_mat_model
      logical :: local_debug
      double precision ::
     &   zero,
     &   top_stress_n_avg(nstrs), bott_stress_n_avg(nstrs),
     &   top_eps_n_avg(nstr), bott_eps_n_avg(nstr),
     &   top_stress_n(nstrs,mxgp), bott_stress_n(nstrs,mxgp),
     &   top_eps_n(nstr,mxgp), bott_eps_n(nstr,mxgp),
     &   top_local_vals(nonlocal_shared_state_size),
     &   bott_local_vals(nonlocal_shared_state_size)
c
      data  zero  / 0.0d00 /
c
c           build averages of stresses and strains at start of load
c           step for the two solid elements connected to each cohesive
c           element in the block (top surface solid element, bottom
c           surface solid element).
c
c           this requires accessing blocks of stress-strain results
c           for the solid elements stored across the blocked
c           data structure for the full model.
c
c           call lower-level routines to hide and simplify extracting
c           stress-strain results from 3D arrays for the solid
c           elements.
c
c           return averages of integration point values for
c           the two solid elements. these are in model (global)
c           coordinates.
c
c           nonlocal data. the umat (and other solid matl models in
c           future) may have supplied a vector of material state
c           variables for use in nonlocal cohesive computations.
c           stress updating created an average value for solid elements
c           element for each state variable. pull those averaged values
c           into block data for top and bottom surface elements.
c           we give cohesive material values for converged solution
c           at n.
c           info for current block of cohesive elements

      span   = local_work%span
      felem  = local_work%felem
      iout   = local_work%iout
      iter   = local_work%iter
      step   = local_work%step
      blk    = local_work%blk
c
      local_debug = .false.
c
      if( local_debug ) then
        write(iout,9000)
        write(iout,9020) blk, span
      end if

      do rel_elem = 1, span
c
c           1. get two solid elements attached to this cohesive elem.
c              save in block local work to send to cohesive model.
c
        elem_top  = solid_interface_lists(blk)%list(rel_elem,1)
        elem_bott = solid_interface_lists(blk)%list(rel_elem,2)
        local_work%top_surf_solid_elements(rel_elem) = elem_top
        local_work%bott_surf_solid_elements(rel_elem) = elem_bott
        matl_top  = iprops(25,elem_top)
        matl_bott = iprops(25,elem_bott)
        if( local_debug ) write(iout,9030) rel_elem, felem+rel_elem-1,
     &                              elem_top, elem_bott, matl_top,
     &                              matl_bott
c
c           2. zero vectors to store integration point average for the
c              two solid elements
c
        do j = 1, nstrs
           top_stress_n_avg(j)  = zero
           bott_stress_n_avg(j) = zero
        end do
        do j = 1, nstr
           top_eps_n_avg(j)  = zero
           bott_eps_n_avg(j) = zero
        end do
c
c           3. get table of integration point values for strains,
c              stresses in solid element attached to top surface.
c              also get number of integration points for the solid
c              element. compute averages of integration point values.
c              these are passed on later to the cohesive element.
c              two solid elements
c
c              repeat for solid element on bottom surface.
c
c              either solid may be zero - cohesive element is attached
c              to a symmetry plane.
c
        if( elem_top .gt. 0 ) then
           call recstr_get_solid_results( elem_top, top_stress_n,
     &                                    top_eps_n, ngp_top )
           call recstr_make_avg( nstrs, ngp_top, top_stress_n,
     &                           top_stress_n_avg )
           call recstr_make_avg( nstr, ngp_top, top_eps_n,
     &                           top_eps_n_avg )
         end if
c
        if( elem_bott .gt. 0 ) then
           call recstr_get_solid_results( elem_bott, bott_stress_n,
     &                                    bott_eps_n, ngp_bott )
           call recstr_make_avg( nstrs, ngp_bott, bott_stress_n,
     &                           bott_stress_n_avg )
           call recstr_make_avg( nstr, ngp_bott, bott_eps_n,
     &                           bott_eps_n_avg  )
        end if
c
c           4. make average strains. stresses the same for two solid
c              elements when cohesive element is on symmetry plane
c
        if( elem_bott .eq. 0 ) then
           bott_stress_n_avg(1:nstrs) = top_stress_n_avg(1:nstrs)
           bott_eps_n_avg(1:nstr)     = top_eps_n_avg(1:nstr)
        end if
        if( elem_top .eq. 0 ) then
           top_stress_n_avg(1:nstrs)  = bott_stress_n_avg(1:nstrs)
           top_eps_n_avg(1:nstr)      = bott_eps_n_avg(1:nstr)
        end if
c
c           5. put average vectors for two solid elements into local
c              structure for this block of cohesive elements

        do j = 1, nstrs
          local_work%top_surf_solid_stresses_n(rel_elem,j) =
     &                   top_stress_n_avg(j)
          local_work%bott_surf_solid_stresses_n(rel_elem,j) =
     &                   bott_stress_n_avg(j)
        end do
        do j = 1, nstr
          local_work%top_surf_solid_eps_n(rel_elem,j)  =
     &                    top_eps_n_avg(j)
          local_work%bott_surf_solid_eps_n(rel_elem,j) =
     &                    bott_eps_n_avg(j)
        end do
c

c           6. nonlocal shared state data for solid element material
c              models. pull already averaged values from global data
c              structures for top & bottom surface solids. do the
c              same trick above when either top or bottom elements
c              are on symmetry plane.
c
c              get converged values at end of step n.
c
        n = nonlocal_shared_state_size
        call recstr_get_solid_nonlocal( elem_top, elem_bott,
     &               top_local_vals, bott_local_vals, iout, n )
        do j = 1, n
         local_work%nonlocal_stvals_bott_n(rel_elem,j) =
     &              bott_local_vals(j)
         local_work%nonlocal_stvals_top_n(rel_elem,j) =
     &              top_local_vals(j)
        end do
c
c           7. get the WARP3D material type id for the top
c              and bottom solid elements
c
        call recstr_get_solid_matl( elem_top, elem_bott,
     &                  top_mat_model, bott_mat_model )
        local_work%top_solid_matl(rel_elem) = top_mat_model
        local_work%bott_solid_matl(rel_elem) = bott_mat_model
c
c           6. debug output
c
        if( local_debug ) then
          write(iout,9070) top_stress_n_avg(1:nstrs)
          write(iout,9075) top_eps_n_avg(1:nstr)
          write(iout,9080) bott_stress_n_avg(1:nstrs)
          write(iout,9085) bott_eps_n_avg(1:nstr)
        end if
c
      end do  ! rel_elem loop
c
      if ( local_debug ) write(iout,9010)
      return
c
c
 9000 format(" .... entered recstr_build_cohes_nonlocal ....")
 9010 format(" .... leaving recstr_build_cohes_nonlocal ....")
 9020 format(" ....    processing cohesive block, span: ",2i5 )
 9030 format(10x,"rel_elem, cohes elem, solid top, solid bott,',
     & ' matl top, matl_bott: ",i3, 3i8,2i8)
 9050 format(/,".... summary of nonlocal cohesive setup. block: ",i5)
 9060 format(10x,"rel_elem, cohes elem, ele top, ele bott:",i4,3i8)
 9070 format(12x,'sig top: ',9f10.3)
 9075 format(12x,'eps top: ',6e14.6)
 9080 format(12x,'sig bot: ',9f10.3)
 9085 format(12x,'eps bot: ',6e14.6)
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine recstr_get_solid_matl             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 03/1/2013 rhd              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine recstr_get_solid_matl( elem_top, elem_bott,
     &              top_model, bott_model )
      use global_data ! old common.main
c
      implicit none
c
      integer :: elem_top, elem_bott, top_model, bott_model
c
c           extract WARP3D material type for top & bottom
c           solid element
c
      if( elem_top .ne. 0 ) top_model = iprops(25,elem_top)
      if( elem_bott .ne. 0 ) bott_model = iprops(25,elem_bott)
c
c           for symmetry case, make top and bottom surface the same
c
      if( elem_top .eq. 0 ) top_model = bott_model
      if( elem_bott .eq. 0 ) bott_model = top_model
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine recstr_get_solid_nonlocal         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 03/4/2013 rhd              *
c     *                                                              *
c     *    shared nonlocal state values from connected solid elems   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine recstr_get_solid_nonlocal( elem_top, elem_bott,
     &                top_local_vals, bott_local_vals, iout, nsize )
      use elem_block_data, only:  nonlocal_flags, nonlocal_data_n
c
      implicit none
c
c           parameter declarations
c
      integer :: elem_top, elem_bott, iout, nsize
      double precision ::
     &  top_local_vals(nsize), bott_local_vals(nsize)
c
c           local declarations
c
      integer :: n
      logical :: chk1, chk2
c
      n = nsize
c
c           extract nonlocal state values for top solid element
c           converged values at end of last step (n)

      if( elem_top .ne. 0 ) then
        chk1 = nonlocal_flags(elem_top)
        chk2 = allocated( nonlocal_data_n(elem_top)%state_values )
        if( chk1 .and. chk2 ) then
           top_local_vals(1:n) =
     &     nonlocal_data_n(elem_top)%state_values(1:n)
        else
            write(iout,9000) elem_top
            call die_abort
        end if
      end if
c
c           extract nonlocal state values for bottom solid element
c
      if( elem_bott .ne. 0 ) then
        chk1 = nonlocal_flags(elem_bott)
        chk2 = allocated( nonlocal_data_n(elem_bott)%state_values )
        if( chk1 .and. chk2 ) then
           bott_local_vals(1:n) =
     &     nonlocal_data_n(elem_bott)%state_values(1:n)
        else
            write(iout,9000) elem_bott
            call die_abort
        end if
      end if
c
c           for symmetry case, make top and bottom surface nonlocal state
c           values the same
c
      if( elem_top .eq. 0 )  top_local_vals(1:n) = bott_local_vals(1:n)
      if( elem_bott .eq. 0 ) bott_local_vals(1:n) = top_local_vals(1:n)
c
      return
c
 9000 format(">>>> FATAL ERROR: recstr_get_solid_nonlocal. elem: ",i8,
     &   /,  "                  job terminated..." )
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine recstr_get_solid_results          *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/22/2013 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine recstr_get_solid_results( solid_elem, stress_n,
     &                                     eps_n, ngp_solid )
      use global_data ! old common.main
c
      use main_data, only: elems_to_blocks
      use elem_block_data, only: urcs_n_blocks, eps_n_blocks
      implicit none
c
c           global declarations
c
      integer :: solid_elem, ngp_solid
      double precision ::
     &    stress_n(nstrs,mxgp), eps_n(nstr,mxgp)
c
c           local declarations
c
      integer :: blk, span, felem, rel_elem, i
      logical :: local_debug
      double precision :: zero

c
c           given solid element, build 2D arrays of strain-stress
c           values at integration points. use service routine
c           to let compiler work out indexing.
c
      local_debug = .false.
c
c           info for block that contains the solid element
c
      blk       = elems_to_blocks(solid_elem,1)
      span      = elblks(0,blk)
      felem     = elblks(1,blk)
      rel_elem  = solid_elem - felem + 1
      ngp_solid = iprops(6,solid_elem)  ! note -- returned
c
      call recstr_copy_results( rel_elem, stress_n(1,1),
     &          urcs_n_blocks(blk)%ptr(1), nstrs, ngp_solid, span )
c
      call recstr_copy_results( rel_elem, eps_n(1,1),
     &          eps_n_blocks(blk)%ptr(1), nstr, ngp_solid, span )
c
      if( local_debug ) then
        write(out,9000) solid_elem, blk, span, felem, rel_elem,
     &                  ngp_solid
        write(out,9005)
        do i = 1, ngp_solid
          write(out,9070) i, stress_n(1:6,i)
        end do
        write(out,9010)
        do i = 1, ngp_solid
          write(out,9075) i, eps_n(1:6,i)
        end do
        write(out,9100)
      end if
c
      return
c
 9000 format(" .... recstr_get_solid_results ....",
     & /, 10x,"solid elem, blk, span, felem, rel_elem, ngp: ",i8,
     &            5i6)
 9100 format(" .... leaving recstr_get_solid_results ....")
 9005 format(12x,"stresses at integration points:")
 9010 format(12x,"strains at integration points:")
 9070 format(15x,i2,9f10.3)
 9075 format(15x,i2,6e14.6)
c
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine recstr_copy_results               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/22/2013 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine recstr_copy_results( kindex_to_copy,
     &                                outmat,
     &                                in3dmat, nrow, ncol, nz )
      implicit  none
      integer  :: kindex_to_copy, nrow, ncol, nz, i, j
      double precision ::
     &  outmat(nrow,ncol), in3dmat(nrow,ncol,nz)
c
c           pull results from k-plane of 3D array into 2D array.
c           used as it exposes structure of 3D array. compiler
c           should inline this routine.
c
      do j = 1, ncol
        do i = 1, nrow
         outmat(i,j) = in3dmat(i,j,kindex_to_copy)
        end do
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine recstr_make_avg                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 01/22/2013 rhd                 *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine recstr_make_avg( nrows, ncols, matrix, averages )
      implicit  none
      integer :: nrows, ncols, i, j
      double precision ::
     &  averages(nrows), matrix(nrows,ncols)
c
c           compute the average of each row in matrix. averages was
c           zeroed before entry. compiler should inline this routine.
c
      do j = 1, ncols
         do i = 1, nrows
           averages(i) = averages(i) + matrix(i,j)
         end do
      end do
c
      do i = 1, nrows
         averages(i) = averages(i) / real(ncols)
      end do
c
      return
      end
