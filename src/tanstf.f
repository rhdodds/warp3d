c     ****************************************************************
c     *                                                              *
c     *                      subroutine tanstf                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 01/23/12 rhd               *
c     *                                                              *
c     *     drive computation of all element [K]s. can be symmetric  *
c     *     (store upper-triangle) or asymmetric (store full [K])    *
c     *     get [K] rotated into constraint compatible coords as     *
c     *     needed                                                   *                              
c     *                                                              *
c     ****************************************************************
c
c
      subroutine tanstf( first, now_step, now_iter ) 
c
      use elem_block_data, only : estiff_blocks, edest_blocks
      use main_data,       only : asymmetric_assembly
c
      implicit integer (a-z)
$add common.main
c
c                       parameter dclarations
c
      logical ::  first 
c                       local declarations
c
#dbl      double precision ::
#sgl      real ::
     &  zero, start_estiff, end_estiff, omp_get_wtime
      logical :: local_debug
      data local_debug, zero / .false., 0.0d00 /
c
c
c             For MPI:
c               alert MPI worker processors that we are in the tanstf routine.
c               also send workers the step and iteration number
c
      if( local_debug ) write(*,*) '... start of tanstf...'
      call wmpi_alert_slaves ( 4 )
      call wmpi_bcast_int ( now_step )
      call wmpi_bcast_int ( now_iter )
c
      call thyme( 2, 1 )
c
c             allocate block data arrays for element
c             stiffness matrices. we compute the element
c             stiffnesses in a local block array then
c             copy to globally allocated array. reduces
c             access to globals in parallel and gets better cache use.

      call estiff_allocate ( 1 )
      if( local_debug ) write(*,*) ' @ 1 tanstf'
c
c             compute element nonlinear [k] matrices. data structures
c             are set up so this can be done in parallel over blocks
c             using threads, with vectorization of loops inside of blocks.
c
c             if we are using MPI:
c               elblks(2,blk) holds which processor owns the block. If
c               we don't own the block, then skip its computation.
c             if we are using the non-MPI version:
c               elblks(2,blk) is all equal to 0, so all blocks
c               are processed.
c
c             this code runs serial, mpi, omp, or omp under mpi.
c
      if( local_debug ) write(*,*) ' @ 2 tanstf'
c
      call omp_set_dynamic( .false. )
      if( local_debug ) then
         start_estiff = omp_get_wtime()
         write(out,*) '... num_threads: ',num_threads
      end if
c
c$OMP PARALLEL DO  PRIVATE( blk, now_thread )
c$OMP&            SHARED( nelblk, elblks, first, now_iter,
c$OMP&                    now_step )
       do blk = 1, nelblk
         if( elblks(2,blk) .ne. myid ) cycle
         now_thread = omp_get_thread_num() + 1
         call do_nlek_block( blk, first, now_iter, now_step )
      end do
c$OMP END PARALLEL DO
c
      if( local_debug ) then
         end_estiff = omp_get_wtime()
         write(out,*) '>> threaded estiff: ', end_estiff - start_estiff
      end if
c
c                       set flags indicating that the tangent stiffness
c                       matrix has been calculated.
c
      tkcomp = .true.
      lkcomp = .false.
      call thyme(2,2)
c      
      return
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine do_nlek_block                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/27/2015 rhd              *
c     *                                                              *
c     *     computes the global nonlinear stiffness                  *
c     *     matrices for a block of elements. the data structures    *
c     *     enable this routine to run in thread parallel with       *
c     *     local vectorization inside a block. the matrix           *
c     *     coordinates. common.main and module damage_data are      *
c     *     read-only in this process                                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine do_nlek_block( blk, first, now_iter, now_step )
c
      use elem_block_data,   only : estiff_blocks, cdest_blocks,
     &                              edest_blocks
      use elem_extinct_data, only : dam_blk_killed, dam_state
c
      use main_data,         only : trn, incid, incmap,
     &                              fgm_node_values_defined,
     &                              cohesive_ele_types,
     &                              linear_displ_ele_types,
     &                              adjust_constants_ele_types,
     &                              axisymm_ele_types,
     &                              temperatures_ref,
     &                              nonlocal_analysis,
     &                              asymmetric_assembly,
     &                              dmatprp, imatprp
c
      use damage_data, only : dam_ptr, growth_by_kill
c 
      use contact, only : use_contact     
c
      implicit integer (a-z)
$add common.main
c
c                       parameter declarations
c
      logical :: first
      integer :: blk, now_iter, now_step
c
c                       local declarations 
c
$add include_tan_ek
#dbl      double precision ::
#sgl      real ::
     &  zero, lambda(mxvl,3,3) ! on stack
      logical :: local_debug, geo_non_flg, bbar_flg, 
     &           symmetric_assembly, block_is_killable
      data local_debug, zero / .false., 0.0d0 /
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
      surface        = iprops(26,felem)
      matnum         = iprops(38,felem)
c
      local_work%felem          = felem
      local_work%blk            = blk
      local_work%num_threads    = num_threads
      local_work%elem_type      = elem_type
      local_work%int_order      = int_order
      local_work%mat_type       = mat_type
      local_work%matnum         = matnum
      local_work%iout           = out
      local_work%num_enodes     = num_enodes
      local_work%num_enode_dof  = num_enode_dof
      local_work%totdof         = totdof
      local_work%geo_non_flg    = geo_non_flg
      local_work%bbar_flg       = bbar_flg
      local_work%num_int_points = num_int_points
      local_work%span           = span
      local_work%utsz           = utsz
      local_work%beta_fact      = beta_fact
      local_work%eps_bbar       = eps_bbar
      local_work%dt             = dt
      local_work%time_n         = total_model_time
      local_work%first          = first
      local_work%iter           = now_iter
      local_work%step           = now_step
      local_work%temperatures       = temperatures
      local_work%temperatures_ref   = temperatures_ref
      local_work%qbar_flag      = qbar_flag
      local_work%cohes_type     = cohes_type
      local_work%surface        = surface
      local_work%fgm_enode_props = fgm_node_values_defined
      local_work%is_cohes_elem  = cohesive_ele_types(elem_type)
      local_work%is_cohes_nonlocal  = nonlocal_analysis .and.
     &                                local_work%is_cohes_elem
      local_work%linear_displ_elem = linear_displ_ele_types(elem_type)
      local_work%adjust_const_elem =
     &                             adjust_constants_ele_types(elem_type)
      local_work%is_axisymm_elem = axisymm_ele_types(elem_type)
      local_work%is_solid_matl  = .not. local_work%is_cohes_elem
      local_work%is_umat        = mat_type .eq. 8
      local_work%is_deform_plas = mat_type .eq. 2
      local_work%is_crys_pls    = mat_type .eq. 10
      local_work%cep_sym_size       = 21
c
      if( local_work%is_umat ) call material_model_info( felem, 0, 3,
     &                                 local_work%umat_stress_type )
      if( local_work%is_cohes_elem ) local_work%cep_sym_size = 6
      symmetric_assembly = .not. asymmetric_assembly
c
c             See if we're actually an interface damaged material.
c        
c
      local_work%is_inter_dmg = .false.
      if( iprops(42,felem) .ne. -1 ) then
        local_work%is_inter_dmg = .true.
        local_work%inter_mat = iprops(42,felem)
        local_work%macro_sz = imatprp(132, local_work%inter_mat)
        local_work%cp_sz = imatprp(133, local_work%inter_mat)
        tm = local_work%inter_mat
c
        local_work%sv(1) = dmatprp(116, tm)
        local_work%sv(2) = dmatprp(117, tm)
        local_work%sv(3) = dmatprp(118, tm)
c
        local_work%lv(1) = dmatprp(119, tm)
        local_work%lv(2) = dmatprp(120, tm)
        local_work%lv(3) = dmatprp(121, tm)
c
        local_work%tv(1) = dmatprp(122, tm)
        local_work%tv(2) = dmatprp(123, tm)
        local_work%tv(3) = dmatprp(124, tm)
      end if
c
      call chk_killed_blk( blk, local_work%killed_status_vec,
     &                     local_work%block_killed )
c
c             check if blk has all killed elements -- if so skip
c             all calculations. just zero element [k]s
c
      nrow_ek = utsz
      if( asymmetric_assembly ) nrow_ek = totdof**2
c      
      if( growth_by_kill ) then  ! note return inside here
        if( local_work%block_killed ) then
          if( local_debug ) write (*,*)'blk ',blk,' killed, skip.'
          estiff_blocks(blk)%ptr(1:nrow_ek,1:span) = zero
          return  ! no tanstf_deallocate  needed
        end if
      end if
c
c             build data structures for elements in this block.
c             this is a gather operation on nodal coordinates,
c             nodal displacements, stresses at time n+1, material
c             state/history data needed to form consistent
c             tangent matrices. the gathered data for this block
c             is stored in the "local_work" data structure with base
c             definition on the stack (each thread thus has a 
c             private copy). allocatables inside local_work are
c             also unique to the thread.
c
      if( local_debug ) write(out,9100) blk, span, felem, mat_type, 
     &            num_enodes, num_enode_dof, totdof, num_int_points
c
      call tanstf_allocate( local_work )
c
      call dptstf( span,
     &             edest_blocks(blk)%ptr(1,1),
     &             cdest_blocks(blk)%ptr(1,1),
     &             incid(incmap(felem)),
     &             felem,
     &             num_int_points,
     &             num_enodes,
     &             num_enode_dof,
     &             geo_non_flg,
     &             totdof,
     &             mat_type,
     &             local_work%trn_e_flags,
     &             local_work%trn_e_block,
     &             local_work%ce,
     &             local_work%trne,
     &             local_work%trnmte,
     &             local_work%ue,
     &             local_work%due,
     &             local_work%cp,
     &             local_work%icp, trn,
     &             elem_type, local_work%surface,
     &             local_work%is_cohes_elem )
c
      if( local_debug ) write(*,*) '.. calling dptstf_blocks..'
      call dptstf_blocks( blk, span, incid(incmap(felem)), felem,
     &                    num_int_points, num_enodes, num_enode_dof,
     &                    geo_non_flg, totdof, mat_type, local_work )
c
c             compute updated tangent stiffness for each element
c             in the block. element stiffnesses are stored
c             in upper triangular form.
c
      if( local_debug ) write(out,9200) blk, span, felem, elem_type,
     &                 int_order, geo_non_flg, bbar_flg
c
c             compute element stiffness for the block. note we
c             pass first element in block of props table.
c
      ispan  = span   ! just protects span value
      call rktstf( props(1,felem), iprops(1,felem),
     &             lprops(1,felem), estiff_blocks(blk)%ptr(1,1),
     &             nrow_ek, ispan, local_work )
c
c             check if this block has any killed elements -- if so,
c             zero computed nonlinear stifffness matrices for killed
c             elements (the [D] may not have been zeroed)
c
      if( growth_by_kill ) then
        block_is_killable = iand( iprops(30,felem),2 ) .ne. 0
        if( block_is_killable ) then
@!DIR$ LOOP COUNT MAX=###  
         do relem = 1, span
           element = felem + relem - 1
           if( dam_ptr(element) .eq. 0 ) cycle
           if( dam_state(dam_ptr(element)) .ne. 0 )
     &       call tanstf_zero_vector( estiff_blocks(blk)%ptr(1,relem),
     &                                nrow_ek )
         end do
        end if 
      end if 
c
c              if contact, add penalty stiffnesses
c
      if( use_contact ) then
        call contact_stfadd( span, felem, totdof,
     &     edest_blocks(blk)%ptr(1,1),
     &     estiff_blocks(blk)%ptr(1,1), nrow_ek, num_enodes,
     &     incid(incmap(felem)) )
      end if
c
c             release all allocated data for block. data never 
c             allocated above for a killed block of elements
c
      call tanstf_deallocate( local_work )
      return
c
 9100 format(5x,'>>> ready to call dptstf:',
     &     /,10x,'blk, span, felem, mat_model:      ',4i10,
     &     /,10x,'num_enodes, num_enode_dof, totdof:',3i10,
     &     /,10x,'num_int_points:                   ',i10 )
 9200 format(5x,'>>> ready to call rktstf:',
     &     /,10x,'blk, span, felem                 :',3i10,
     &     /,10x,'elem_type, int_order, geo_non_flg:',2i10,l10,
     &     /,10x,'bbar_flg:                        :',l10 )
c
      end
c     ****************************************************************
c     *                                                              *
c     *                subroutine estiff_allocate                    *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *                last modified : 9/27/2015 rhd                 *
c     *                                                              *
c     *     create the blocked data structure for storage of element *
c     *     stiffness matrices                                       *
c     *        symmetric - upper triangle including diagonal         *
c     *       asymmetric - lower and upper triangle                  *
c     *     runs outside any threaded region.                        *
c     *                                                              *
c     ****************************************************************
c
      subroutine estiff_allocate( type )
c
      use elem_block_data, only:  estiff_blocks, pcm_blocks 
      use main_data, only: asymmetric_assembly
c
      implicit integer (a-z)
$add common.main
      logical :: myblk
c
c
c            the data structure is a 2-D array for each element
c            block (dynamically allocated). the arrays are
c            hung from a dynamically allocated pointer vector
c            of length nelblk.
c
c            type  = 1, 4 allocate stiffness blocks.
c            type  = 2 no longer used. old ebe
c            type =  3 not used
c            type =  5 deallocate blocks
c
      select case ( type ) ! careful. type could be integer constant
      case ( 1, 4 ) ! symmetric or asymmetric.
c
         if( .not. allocated( estiff_blocks ) ) then
           allocate( estiff_blocks(nelblk), stat=iok )
           if( iok .ne. 0 ) then
              call iodevn( idummy, iout, dummy, 1 )
              write(iout,9100) iok
              call die_abort
           end if
           do blk = 1, nelblk
            nullify( estiff_blocks(blk)%ptr )
           end do
         end if
c
c             MPI:
c               elblks(2,blk) holds which rank owns the
c               block.  For worker ranks, if we don't own
c               the block, don't allocate it.  For root,
c               allocate everything for type = 4.
c             for threads only, allocate all blocks
c             skip blocks if the pointer to them is associated.
c
         do blk = 1, nelblk
            myblk = myid .eq. elblks(2,blk)
            if( myid .eq. 0 .and. type .eq. 4 ) myblk = .true.
            if( .not. myblk ) cycle
            if( associated(estiff_blocks(blk)%ptr) ) cycle
            felem         = elblks(1,blk)
            num_enodes    = iprops(2,felem)
            num_enode_dof = iprops(4,felem)
            totdof        = num_enodes * num_enode_dof
            span          = elblks(0,blk)
            utsz          = ((totdof*totdof)-totdof)/2 + totdof
            nterms        = utsz
            if( asymmetric_assembly ) nterms = totdof * totdof
            allocate( estiff_blocks(blk)%ptr(nterms,span),stat=iok )
            if( iok .ne. 0 ) then
               call iodevn( idummy, iout, dummy, 1 )
               write(iout,9100) iok
            end if
         end do
c
      case( 5 ) ! deallocate estiff_blocks
c
         if( .not. allocated( estiff_blocks ) ) return
         if( myid .ne. 0 ) return ! only do this on root
         do blk = 1, nelblk
            myblk = myid .eq. elblks(2,blk)
            if( .not. myblk ) cycle
            deallocate( estiff_blocks(blk)%ptr, stat=iok )
            if( iok .ne. 0 ) then
               call iodevn( idummy, iout, dummy, 1 )
               write(iout,9200) iok
            end if
            nullify( estiff_blocks(blk)%ptr )
         end do
c
      case default 
         call iodevn( idummy, iout, dummy, 1 )
         write(iout,9300) 
         call die_abort
      end select
c
      return      
c
 9100 format('>> FATAL ERROR: estiff_allocate, memory allocate failure',
     &  /,   '                status= ',i5,
     &  /,   '                job terminated' )
 9200 format('>> FATAL ERROR: estiff_allocate, memory deallocate',
     &  /,   '                failure. status= ',i5,
     &  /,   '                job terminated' )
 9300 format('>> FATAL ERROR: estiff_allocate, unknown state',
     &  /,   '                job terminated' )
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine dptstf_blocks                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/27/2015 rhd              *
c     *                                                              *
c     *     this subroutine creates a separate copy of element       *
c     *     data necessary for the tangent stiffness computation of  *
c     *     each element in a block of similar                       *
c     *     elements. processes only data stored globally in         *
c     *     blocked data structures.                                 *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine dptstf_blocks(
     &   blk, span, belinc, felem, ngp, nnode, ndof, geonl, totdof,
     &   mat_type, local_work )
c
      use main_data, only:        temper_nodes, dtemp_nodes,
     &                            temper_elems, dtemp_elems,
     &                            fgm_node_values,
     &                            fgm_node_values_defined,
     &                            fgm_node_values_cols,
     &                            temperatures_ref,
     &                            temper_nodes_ref
c
      use elem_block_data, only:  history_blocks, rts_blocks,
     &                            rot_n1_blocks, history1_blocks,
     &                            eps_n1_blocks, urcs_n1_blocks,
     &                            history_blk_list
c
      implicit integer (a-z)
$add common.main
$add include_tan_ek
c
c           parameter declarations
c
      dimension :: belinc(nnode,*)
      logical ::   geonl
c
c           local declarations
c
      logical :: local_debug
#dbl      double precision ::
#sgl      real ::
     & zero
      data zero, local_debug / 0.0d00, .false. /
c
c
c           get data not dependent on material model type first. need
c           values for current estimate of end-of-step conditions
c               1) rotation matrices at integration points for geonl
c                  to unrotate cuchy stresses into global
c               2) element histories (allocate the local storage block)
c               3) unrotated cauchy stresses
c               4) nodal temperatures for elements in the block to
c                  support temperature dependent material properties
c
c           History data:
c            o The global blocks are sized(hist_size,ngp,span)
c            o The local block is sized (span,hist_size,ngp).
c              This makes it possible to pass a 2-D array slice for
c              all elements of the block for a single gauss point.
c            o History data is not needed for warp3d_umat since
c              we'll just extract already computed [Dt]s in global
c              blocks
c
      if( geonl )  call gastr( local_work%rot_blk_n1,
     &                          rot_n1_blocks(blk)%ptr(1), ngp,
     &                          9, span )
c
      hist_size = history_blk_list(blk)
      local_work%hist_size_for_blk = hist_size
c
      if( .not. local_work%is_umat ) then
         allocate( local_work%elem_hist1(span,hist_size,ngp),
     &             local_work%elem_hist(span,hist_size,ngp) )
         call dptstf_copy_history(
     &       local_work%elem_hist1(1,1,1), history1_blocks(blk)%ptr(1),
     &       ngp, hist_size, span )
      end if
c
      call gastr( local_work%urcs_blk_n1, urcs_n1_blocks(blk)%ptr(1),
     &            ngp, nstrs, span )
c
c                 gather nodal and element temperature change
c                 over load step (if they are defined). we
c                 construct a set of incremental nodal temperatures
c                 for each element in block.
c
      if( temperatures ) then
        if( local_debug )  write(*,9610)
        call gadtemps( dtemp_nodes, dtemp_elems(felem), belinc,
     &                 nnode, span, felem, local_work%dtemps_node_blk,
     &                 mxvl )
      else
        call tanstf_zero_vector( local_work%dtemps_node_blk,
     &                 mxvl*mxndel )
      end if
c
c           gather reference temperatures for element nodes from the
c           global vector of reference values at(if they are defined).
c           construct a set of reference nodal temperatures for each
c           element in block.
c
      if( temperatures_ref ) then
        call gartemps( temper_nodes_ref, belinc, nnode, span,
     &                 felem, local_work%temps_ref_node_blk, mxvl )
      else
        call tanstf_zero_vector( local_work%temps_ref_node_blk,
     &                           mxvl*mxndel )
      end if
c
c                 build nodal temperatures for elements in the block
c                 at end of step (includes both imposed nodal and element
c                 temperatures)
c
      if( local_debug )  write(*,9620)
      call gatemps( temper_nodes, temper_elems(felem), belinc,
     &              nnode, span, felem, local_work%temps_node_blk,
     &              mxvl, local_work%dtemps_node_blk,
     &              local_work%temps_node_to_process )
c
c
c
      select case( mat_type )
c      ------------------------
c
      case( 1 )  ! bilinear
c     =========
c
c           gather the deviatoric components of the trial elastic
c           stress increment for the step. this used by material
c           model 1 for consistent tangent stiffness evaluation.
c           gather plasticity parameters and state variables for
c           the vectorized plasticity model
c
        call gastr( local_work%rtse, rts_blocks(blk)%ptr(1), ngp, nstr,
     &               span )

c
      case( 2 )   ! deformation plasticity
c     =========
c
        call gastr( local_work%ddtse, eps_n1_blocks(blk)%ptr(1),
     &              ngp, nstr, span )
        call dptstf_copy_history(
     &           local_work%elem_hist(1,1,1), 
     &           history_blocks(blk)%ptr(1), ngp, hist_size, span )
c
      case( 4 )   ! cohesive
c     =========
c
         call gastr( local_work%ddtse, eps_n1_blocks(blk)%ptr(1),
     &               ngp, nstr, span )
         call dptstf_copy_history(
     &         local_work%elem_hist(1,1,1), history_blocks(blk)%ptr(1),
     &         ngp, hist_size, span )
         if( local_work%is_cohes_nonlocal ) 
     &      call tanstf_build_cohes_nonlocal( local_work )
c
      case( 3, 5, 6, 7, 10 )
c     ==================
c
c           gather history data at n, n+1, trial elastic
c           stresses at n+1 for these models
c            (3) -- general mises plasticity with optional Gurson  model
c            (5) -- cyclic plasticity
c            (6) -- creep
c            (7) -- mises + hyrdogen effects
c           (10) -- crystal plasticity
c
         call dptstf_copy_history(
     &         local_work%elem_hist(1,1,1), history_blocks(blk)%ptr(1),
     &         ngp, hist_size, span )
         call gastr( local_work%rtse, rts_blocks(blk)%ptr(1),
     &               ngp, nstr, span )
c
      case( 8 )   ! umat. nothing to do for now.
c     =========
c
      continue  ! just so it is clear that 8 is umat
c
      case default
          write(local_work%iout,*) '>>> invalid material model number'
          write(local_work%iout,*) '    in dptstf_blocks'
          call die_abort
      end select
c
c                 if the model has fgm properties at the model nodes,
c                 build a table of values for nodes of elements in the
c                 block
c
      if( fgm_node_values_defined ) then
        do j = 1,  fgm_node_values_cols
          do i = 1, nnode
@!DIR$ LOOP COUNT MAX=###  
            do k = 1, span
              local_work%enode_mat_props(i,k,j) =
     &                     fgm_node_values(belinc(i,k),j)
            end do
          end do
        end do
      end if
c
      return
c
 9610 format(12x,'>> gather element incremental temperatures...' )
 9620 format(12x,'>> gather element total temperatures...' )
c
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine duptrans                         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 09/27/2015 rhd             *
c     *                                                              *
c     *     this subroutine creates a separate copy of element       *
c     *     data necessary to transform displacements from global to *
c     *     the constraint coordinate system at the nodes            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine duptrans( span, felem, trnmte )
c
      use main_data, only : trnmat, trn, incid, incmap
c
      implicit integer (a-z)
$add common.main
c
c           parameters
c
#dbl      double precision
#sgl      real
     &  trnmte(mxvl,mxedof,3)
c
c           for this block of elements, gather the transformation
c           matrices (3x3) used to rotate between global and constraint
c           coordinate systems. if there is no rotation matrix
c           for that node, then trn(node) = .false. and the
c           corresponding matrix trnmat(node)%mat is unallocated.
c           if there is a rotation matrix for the node, then the
c           matrix is stored in trnmat(node)%mat.
c
      nnode  = iprops(2,felem)
      ndof   = iprops(4,felem)
      totdof = nnode * ndof
      call tanstf_zero_vector( trnmte, mxvl*mxedof*3 ) 
c
c           this code below depends on ndof per node = 3
c
      do j = 1, nnode
@!DIR$ LOOP COUNT MAX=###  
         do k = 1, span
            node = incid(incmap(felem+k-1) + j-1)
            jj = (j-1)*3
            if ( trn(node) ) then
               trnmte(k,jj+1,1) = trnmat(node)%mat(1,1)
               trnmte(k,jj+1,2) = trnmat(node)%mat(1,2)
               trnmte(k,jj+1,3) = trnmat(node)%mat(1,3)
               trnmte(k,jj+2,1) = trnmat(node)%mat(2,1)
               trnmte(k,jj+2,2) = trnmat(node)%mat(2,2)
               trnmte(k,jj+2,3) = trnmat(node)%mat(2,3)
               trnmte(k,jj+3,1) = trnmat(node)%mat(3,1)
               trnmte(k,jj+3,2) = trnmat(node)%mat(3,2)
               trnmte(k,jj+3,3) = trnmat(node)%mat(3,3)
            endif
         end do
      end do
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *               subroutine dptstf_copy_history                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 09/27/2015 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine dptstf_copy_history( local_hist, global_hist,
     &                                ngp, hist_size, span )
      implicit integer (a-z)
c
c               parameter declarations
c
#dbl      double precision
#sgl      real
     & local_hist(span,hist_size,ngp),
     & global_hist(hist_size,ngp,span)
c
      if ( ngp .ne. 8 ) then
        do k = 1, ngp
         do  j = 1, hist_size
@!DIR$ LOOP COUNT MAX=###  
            do  i = 1, span
               local_hist(i,j,k) = global_hist(j,k,i)
            end do
         end do
        end do
        return
      end if
c
c                number of gauss points = 8, unroll.
c
      do  j = 1, hist_size
@!DIR$ LOOP COUNT MAX=###  
        do  i = 1, span
            local_hist(i,j,1) = global_hist(j,1,i)
            local_hist(i,j,2) = global_hist(j,2,i)
            local_hist(i,j,3) = global_hist(j,3,i)
            local_hist(i,j,4) = global_hist(j,4,i)
            local_hist(i,j,5) = global_hist(j,5,i)
            local_hist(i,j,6) = global_hist(j,6,i)
            local_hist(i,j,7) = global_hist(j,7,i)
            local_hist(i,j,8) = global_hist(j,8,i)
        end do
      end do
c
      return
      end


c     ****************************************************************
c     *                                                              *
c     *                   subroutine tanstf_allocate                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/24/2015 rhd              *
c     *                                                              *
c     *     allocate data structure in local_work for updating       *
c     *     element stiffnesses                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine tanstf_allocate( local_work )
      use segmental_curves, only : max_seg_points, max_seg_curves
      implicit integer (a-z)

$add common.main
$add include_tan_ek
#dbl      double precision
#sgl      real
     &  zero
      data zero / 0.0d00 /
c
c
c               history data for block allocated in dptstf_blocks
c
      allocate( local_work%ce(mxvl,mxecor) )
c
      allocate( 
     1 local_work%det_jac_block(mxvl,mxgp),
     2 local_work%shape(mxndel,mxgp),
     3 local_work%nxi(mxndel,mxgp),
     4 local_work%neta(mxndel,mxgp),
     5 local_work%nzeta(mxndel,mxgp),
     6 local_work%gama_block(mxvl,3,3,mxgp), stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 1
           call die_abort
      end if
c
      allocate( local_work%vol_block(mxvl,8,3),
     1  local_work%volume_block(mxvl),
     2  local_work%jac_block(mxvl,3,3),
     3  local_work%b_block(mxvl,mxedof,nstr),
     4  local_work%bt_block(mxvl,nstr,mxedof),
     5  local_work%bd_block(mxvl,mxedof,nstr), stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 2
           call die_abort
      end if
c
      allocate( local_work%ue(mxvl,mxedof),
     1  local_work%due(mxvl,mxedof),
     2  local_work%urcs_blk_n1(mxvl,nstrs,mxgp),
     3  local_work%rot_blk_n1(mxvl,9,mxgp),
     4  local_work%rtse(mxvl,nstr,mxgp),
     5  local_work%ddtse(mxvl,nstr,mxgp), stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 3
           call die_abort
      end if
c
      allocate( 
     1  local_work%temps_node_blk(mxvl,mxndel),
     2  local_work%dtemps_node_blk(mxvl,mxndel),
     3  local_work%temps_ref_node_blk(mxvl,mxndel), stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 4
           call die_abort
      end if
c    
      allocate( local_work%fgm_flags(mxvl,mxndpr) )
      if( local_work%fgm_enode_props ) 
     &   allocate(local_work%enode_mat_props(mxndel,mxvl,mxndpr) )
c
c             local cep, i.e, [Dt] must be 6x6 for all
c             trans[B] [Dt] [B] to work correctly.
c             global cep's can be 6x6 or 3x3
c
      allocate( local_work%cep(mxvl,6,6),
     1  local_work%qn1(mxvl,nstr,nstr),
     2  local_work%cs_blk_n1(mxvl,nstr), stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 5
           call die_abort
      end if
c
      local_work%cep = zero
      local_work%qn1 = zero
      local_work%cs_blk_n1 = zero
c
      if( local_work%mat_type .eq. 5 )
     &  allocate( local_work%mm05_props(mxvl,10) )
      if( local_work%mat_type .eq. 6 )
     &  allocate( local_work%mm06_props(mxvl,10) )
      if( local_work%mat_type .eq. 7 )
     &  allocate( local_work%mm07_props(mxvl,10) )
c
      allocate( 
     a    local_work%e_v(mxvl), 
     1    local_work%beta_v(mxvl),
     2    local_work%nu_v(mxvl), 
     3    local_work%sigyld_v(mxvl),
     4    local_work%n_power_v(mxvl), 
     5    local_work%e_block(mxvl),
     6    local_work%nu_block(mxvl),
     7    local_work%weights(mxgp), stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 6
           call die_abort
      end if
c      
      if( local_work%mat_type .eq. 3 )  then
        allocate( local_work%f0_v(mxvl),
     1            local_work%q1_v(mxvl), 
     2            local_work%q2_v(mxvl),
     3            local_work%q3_v(mxvl), 
     4            local_work%nuc_s_n_v(mxvl),
     5            local_work%nuc_e_n_v(mxvl),
     6            local_work%nuc_f_n_v(mxvl) )
      end if
c
      allocate( local_work%cp(mxedof), local_work%icp(mxutsz,2),
     &          stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 7
           call die_abort
      end if
c
      allocate( 
     1  local_work%trn_e_flags(mxvl),
     2  local_work%trne(mxvl,mxndel),
     3  local_work%trnmte(mxvl,mxedof,mxndof), stat=error ) 
      if( error .ne. 0 ) then
           write(out,9000) 8
           call die_abort
      end if
c
      allocate( local_work%ncrystals(mxvl) )
c
c                always allocate cohes_rot_block. it gets passed as
c                parameter even when block is not interface elements

      allocate( local_work%cohes_rot_block(mxvl,3,3), stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 9
           call die_abort
      end if
      if( local_work%is_cohes_elem ) then
         allocate( local_work%cohes_temp_ref(mxvl),
     1      local_work%cohes_dtemp(mxvl),
     2      local_work%cohes_temp_n(mxvl),
     3      local_work%intf_prp_block(mxvl,max_interface_props),
     4      stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 9
           call die_abort
         end if
      end if
c
      if( local_work%is_cohes_nonlocal ) then
         nlsize = nonlocal_shared_state_size
         allocate( local_work%top_surf_solid_stresses_n1(mxvl,nstrs),
     &      local_work%bott_surf_solid_stresses_n1(mxvl,nstrs),
     &      local_work%top_surf_solid_eps_n1(mxvl,nstr),
     &      local_work%bott_surf_solid_eps_n1(mxvl,nstr),
     &      local_work%top_surf_solid_elements(mxvl),
     &      local_work%bott_surf_solid_elements(mxvl),
     &      local_work%nonlocal_stvals_bott_n1(mxvl,nlsize),
     &      local_work%nonlocal_stvals_top_n1(mxvl,nlsize),
     &      stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 13
           call die_abort
         end if
         allocate( local_work%top_solid_matl(mxvl),
     &         local_work%bott_solid_matl(mxvl), stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 14
           call die_abort
         end if
      end if
c
      return
 9000 format('>> FATAL ERROR: tanstf_allocate'
     &  /,   '                failure status= ',i5,
     &  /,   '                job terminated' )
c
      end
c     ****************************************************************
c     *                                                              *
c     *                   subroutine tanstf_deallocate               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/27/2015 rhd              *
c     *                                                              *
c     *     release data structure in local_work for updating        *
c     *     strains-stresses-internal forces.                        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine tanstf_deallocate( local_work )
      implicit integer (a-z)
$add common.main
$add include_tan_ek
c
      logical :: local_debug
c
      local_debug = .false.
      if( local_debug ) write(out,*) "..tanstf_deall @ 1"
      local_mt = local_work%mat_type
c
      deallocate( local_work%ce,
     1 local_work%det_jac_block,
     2 local_work%shape,
     3 local_work%nxi,
     4 local_work%neta,
     5 local_work%nzeta,
     6 local_work%gama_block, stat=error )
       if( error .ne. 0 ) then
           write(out,9000) 1
           call die_abort
       end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 5"
c
      deallocate( local_work%vol_block,
     1  local_work%volume_block,
     2  local_work%jac_block,
     3  local_work%b_block,
     4  local_work%bt_block,
     5  local_work%bd_block, stat=error )
       if( error .ne. 0 ) then
           write(out,9000) 2
           call die_abort
       end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 10"
c
      deallocate( local_work%ue,
     1  local_work%due,
     2  local_work%urcs_blk_n1,
     3  local_work%rot_blk_n1,
     4  local_work%rtse,
     5  local_work%ddtse, stat=error )
       if( error .ne. 0 ) then
           write(out,9000) 3
           call die_abort
       end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 15"
c
      deallocate(
     1  local_work%temps_node_blk,
     2  local_work%dtemps_node_blk,
     3  local_work%temps_ref_node_blk, stat=error )
       if( error .ne. 0 ) then
           write(out,9000) 4
           call die_abort
       end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 20"
c    
      deallocate( local_work%fgm_flags, stat=error )
       if( error .ne. 0 ) then
           write(out,9000) 200
           call die_abort
       end if
      
      if( local_work%fgm_enode_props ) then
        deallocate( local_work%enode_mat_props, stat=error )
        if( error .ne. 0 ) then
           write(out,9000) 205
           call die_abort
        end if
      end if  
     
      if( local_debug ) write(out,*) "..tanstf_deall @ 30"
c
      deallocate( local_work%cep,
     1  local_work%qn1,
     2  local_work%cs_blk_n1, stat=error )
       if( error .ne. 0 ) then
           write(out,9000) 5
           call die_abort
       end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 40"
c       
      if( local_mt .eq. 5 ) then
        deallocate( local_work%mm05_props, stat=error  )
        if( error .ne. 0 ) then
           write(out,9000) 210
           call die_abort
        end if
      end if
        
      if( local_mt .eq. 6 ) then 
         deallocate( local_work%mm06_props, stat=error )
        if( error .ne. 0 ) then
           write(out,9000) 215
           call die_abort
        end if
      end if
c
      if( local_mt .eq. 7 ) then 
         deallocate( local_work%mm07_props, stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 220
           call die_abort
        end if
      end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 50"
c
      deallocate( 
     a    local_work%e_v, 
     1    local_work%beta_v,
     2    local_work%nu_v, 
     3    local_work%sigyld_v,
     4    local_work%n_power_v, 
     5    local_work%e_block,
     6    local_work%nu_block,
     7    local_work%weights, stat=error )
       if( error .ne. 0 ) then
           write(out,9000) 6
           call die_abort
       end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 55"
c      
      if( local_work%mat_type .eq. 3 )  then
        deallocate( local_work%f0_v,
     1              local_work%q1_v, 
     2              local_work%q2_v,
     3              local_work%q3_v, 
     4              local_work%nuc_s_n_v,
     5              local_work%nuc_e_n_v,
     6              local_work%nuc_f_n_v )
      end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 60"

      deallocate( local_work%cp, local_work%icp, stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 7
           call die_abort
      end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 65"
c
      deallocate( local_work%trn_e_flags, local_work%trne,
     &            local_work%trnmte, stat=error )
       if( error .ne. 0 ) then
           write(out,9000) 8
           call die_abort
       end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 70"
c
      deallocate( local_work%ncrystals, stat=error)
      if( error .ne. 0 ) then
           write(out,9000) 225
           call die_abort
      end if
c
      if( allocated(local_work%elem_hist1) )
     &    deallocate(local_work%elem_hist1, stat=error )
       if( error .ne. 0 ) then
           write(out,9000) 9
           call die_abort
       end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 75"
      if( allocated(local_work%elem_hist) )
     &      deallocate(local_work%elem_hist, stat=error )
       if( error .ne. 0 ) then
           write(out,9000) 10
           call die_abort
       end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 80"
c
      if( local_work%is_cohes_nonlocal ) then
         deallocate( local_work%top_surf_solid_stresses_n1,
     1      local_work%bott_surf_solid_stresses_n1,
     2      local_work%top_surf_solid_eps_n1,
     3      local_work%bott_surf_solid_eps_n1,
     4      local_work%top_surf_solid_elements,
     5      local_work%bott_surf_solid_elements,
     6      local_work%nonlocal_stvals_bott_n1,
     7      local_work%nonlocal_stvals_top_n1,
     8      stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 11
           call die_abort
         end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 85"
         deallocate( local_work%top_solid_matl,
     &         local_work%bott_solid_matl, stat=error )
         if( error .ne. 0 ) then
           write(out,9000) 12
           call die_abort
         end if
      end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 90"
c
c                cohes_rot_block is always allocated. tt gets passed as
c                parameter even when block is not interface elements
c
      deallocate( local_work%cohes_rot_block, stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 13
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
      if( local_debug ) write(out,*) "..tanstf_deall @ 95"
c
      return
c
 9000 format('>> FATAL ERROR: tanstf_deallocate'
     &  /,   '                failure status= ',i5,
     &  /,   '                job terminated' )
c
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine tanstf_zero_vector                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 09/27/2015 rhd             *
c     *                                                              *
c     *     zero a vector of specified length w/ floating zero       *
c     *                                                              *
c     ****************************************************************
c
      subroutine tanstf_zero_vector( vec, n )
#dbl      double precision
#sgl      real
     &  vec(n), zero
      data zero / 0.0d00 /
c
      vec(1:n) = zero
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine tanstf_build_cohes_nonlocal       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 03/4/2013 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine tanstf_build_cohes_nonlocal( local_work  )

      use elem_block_data, only: solid_interface_lists
      use main_data, only:  nonlocal_analysis
c
      implicit integer (a-z)
$add param_def
$add include_tan_ek
c
c           local declarations. arrays are on stack.
c
      logical local_debug
#dbl      double precision
#sgl      real
     &   zero,
     &   top_stress_n1_avg(nstrs), bott_stress_n1_avg(nstrs),
     &   top_eps_n1_avg(nstr), bott_eps_n1_avg(nstr),
     &   top_stress_n1(nstrs,mxgp), bott_stress_n1(nstrs,mxgp),
     &   top_eps_n1(nstr,mxgp), bott_eps_n1(nstr,mxgp),
     &   top_local_vals(nonlocal_shared_state_size),
     &   bott_local_vals(nonlocal_shared_state_size)
c
      data  zero /  0.0d00 /
c
c           build averages of stresses and strains for current estimate
c           of solution at end of step for the two solid elements
c           connected to each cohesive
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
c           element for each state variable.. pull those averaged values
c           into block data for top and bottom surface elements.
c
c           we give cohesive material the nonlocal shared values for
c           current estimate of solution at n+1.
c
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
        if( local_debug ) write(iout,9030) rel_elem, felem+rel_elem-1,
     &                              elem_top, elem_bott
c
c           2. zero vectors to store integration point average for the
c              two solid elements
c
        do j = 1, nstrs
           top_stress_n1_avg(j)  = zero
           bott_stress_n1_avg(j) = zero
        end do
        do j = 1, nstr
           top_eps_n1_avg(j)  = zero
           bott_eps_n1_avg(j) = zero
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
           call tanstf_get_solid_results( elem_top, top_stress_n1,
     &                                    top_eps_n1, ngp_top )
           call tanstf_make_avg( nstrs, ngp_top, top_stress_n1,
     &                           top_stress_n1_avg )
           call tanstf_make_avg( nstr, ngp_top, top_eps_n1,
     &                           top_eps_n1_avg )
         end if
c
        if( elem_bott .gt. 0 ) then
           call tanstf_get_solid_results( elem_bott, bott_stress_n1,
     &                                    bott_eps_n1, ngp_bott )
           call tanstf_make_avg( nstrs, ngp_bott, bott_stress_n1,
     &                           bott_stress_n1_avg )
           call tanstf_make_avg( nstr, ngp_bott, bott_eps_n1,
     &                           bott_eps_n1_avg  )
        end if
c
c           4. make average strains. stresses the same for two solid
c              elements when cohesive element is on symmetry plane
c
        if( elem_bott .eq. 0 ) then
           bott_stress_n1_avg(1:nstrs) = top_stress_n1_avg(1:nstrs)
           bott_eps_n1_avg(1:nstr)     = top_eps_n1_avg(1:nstr)
        end if
        if( elem_top .eq. 0 ) then
           top_stress_n1_avg(1:nstrs)  = bott_stress_n1_avg(1:nstrs)
           top_eps_n1_avg(1:nstr)      = bott_eps_n1_avg(1:nstr)
        end if
c
c           5. put average vectors for two solid elements into local
c              structure for this block of cohesive elements

        do j = 1, nstrs
          local_work%top_surf_solid_stresses_n1(rel_elem,j) =
     &                   top_stress_n1_avg(j)
          local_work%bott_surf_solid_stresses_n1(rel_elem,j) =
     &                   bott_stress_n1_avg(j)
        end do
        do j = 1, nstr
          local_work%top_surf_solid_eps_n1(rel_elem,j)  =
     &                    top_eps_n1_avg(j)
          local_work%bott_surf_solid_eps_n1(rel_elem,j) =
     &                    bott_eps_n1_avg(j)
        end do
c
c
c           6. nonlocal shared state data for solid element material
c              models. pull already averaged values from global data
c              structures for top & bottom surface solids. do the
c              same trick above when either top or bottom elements
c              are on symmetry plane.
c
        n = nonlocal_shared_state_size
        call tanstf_get_solid_nonlocal( elem_top, elem_bott,
     &                    top_local_vals, bott_local_vals, iout, n )
        do j = 1, n
         local_work%nonlocal_stvals_bott_n1(rel_elem,j) =
     &              bott_local_vals(j)
         local_work%nonlocal_stvals_top_n1(rel_elem,j) =
     &              top_local_vals(j)
        end do
c
c           7. get the WARP3D material type id for the top
c              and bottom solid elements
c
        call tanstf_get_solid_matl( elem_top, elem_bott,
     &                  top_mat_model, bott_mat_model, iout )
        local_work%top_solid_matl(rel_elem) = top_mat_model
        local_work%bott_solid_matl(rel_elem) = bott_mat_model
c
c           8. debug output
c
        if( local_debug ) then
          write(iout,9070) top_stress_n1_avg(1:nstrs)
          write(iout,9075) top_eps_n1_avg(1:nstr)
          write(iout,9080) bott_stress_n1_avg(1:nstrs)
          write(iout,9085) bott_eps_n1_avg(1:nstr)
        end if
c
      end do  ! rel_elem loop
c
      if ( local_debug ) write(iout,9010)
      return
c
c
 9000 format(" .... entered tanstf_build_cohes_nonlocal ....")
 9010 format(" .... leaving tanstf_build_cohes_nonlocal ....")
 9020 format(" ....    processing cohesive block, span: ",2i5 )
 9030 format(10x,"rel_elem, cohes elem, solid top, solid bott: ",i3,
     &            3i8)
 9050 format(/,".... summary of nonlocal cohesive setup. block: ",i5)
 9060 format(10x,"rel_elem, cohes elem, ele top, ele bott:",i4,3i8)
 9070 format(12x,'sig n1 top: ',9f10.3)
 9075 format(12x,'eps n1 top: ',6e14.6)
 9080 format(12x,'sig n1 bot: ',9f10.3)
 9085 format(12x,'eps n1 bot: ',6e14.6)
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine tanstf_get_solid_matl             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/27/2013 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine tanstf_get_solid_matl( elem_top, elem_bott,
     &              top_model, bott_model )
c
      implicit integer (a-z)
$add common.main
c
c           parameter declarations
c
c
c           local declarations
c
      logical local_debug
c
      local_debug = .false.
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
c     *                 subroutine tanstf_get_solid_nonlocal         *
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
      subroutine tanstf_get_solid_nonlocal( elem_top, elem_bott,
     &                    top_local_vals, bott_local_vals, iout,
     &                    nsize )
      use elem_block_data, only: nonlocal_flags, nonlocal_data_n1
c
      implicit integer (a-z)
c
c           parameter declarations
c
#dbl      double precision
#sgl      real
     &  top_local_vals(nsize), bott_local_vals(nsize)
c
c           local declarations
c
      logical chk1
c
c           extract nonlocal state values for top solid element
c           use estimate of solution at n+1
c
      if( elem_top .ne. 0 ) then
        chk1 = nonlocal_flags(elem_top)
        if( chk1 ) then
           top_local_vals(1:nsize) =
     &     nonlocal_data_n1(elem_top)%state_values(1:nsize)
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
        if( chk1 ) then
           bott_local_vals(1:nsize) =
     &     nonlocal_data_n1(elem_bott)%state_values(1:nsize)
        else
            write(iout,9000) elem_bott
            call die_abort
        end if
      end if
c
c           for symmetry case, make top and bottom surface nonlocal state
c           values the same
c
      n = nsize
      if(  elem_top .eq. 0 ) top_local_vals(1:n) = bott_local_vals(1:n)
      if(  elem_bott .eq. 0 ) bott_local_vals(1:n) = top_local_vals(1:n)
c
      return
c
 9000 format(">>>> FATAL ERROR: tanstf_get_solid_nonlocal. elem: ",i8,
     &   /,  "                  job terminated..." )
c
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine tanstf_get_solid_results          *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/31/2013 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine tanstf_get_solid_results( solid_elem, stress_n1,
     &                                     eps_n1, ngp_solid )
c
      use main_data, only: elems_to_blocks
      use elem_block_data, only: urcs_n1_blocks, eps_n1_blocks
      implicit integer (a-z)
$add common.main
c
c           local declarations
c
      logical local_debug
#dbl      double precision
#sgl      real
     &    zero, stress_n1(nstrs,mxgp), eps_n1(nstr,mxgp)
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
      call tanstf_copy_results( rel_elem, stress_n1(1,1),
     &          urcs_n1_blocks(blk)%ptr(1), nstrs, ngp_solid, span )
c
      call tanstf_copy_results( rel_elem, eps_n1(1,1),
     &          eps_n1_blocks(blk)%ptr(1), nstr, ngp_solid, span )
c
      if( local_debug ) then
        write(out,9000) solid_elem, blk, span, felem, rel_elem,
     &                  ngp_solid
        write(out,9005)
        do i = 1, ngp_solid
          write(out,9070) i, stress_n1(1:6,i)
        end do
        write(out,9010)
        do i = 1, ngp_solid
          write(out,9075) i, eps_n1(1:6,i)
        end do
        write(out,9100)
      end if
c
      return
c
 9000 format(" .... tanstf_get_solid_results ....",
     & /, 10x,"solid elem, blk, span, felem, rel_elem, ngp: ",i8,
     &            5i6)
 9100 format(" .... leaving tanstf_get_solid_results ....")
 9005 format(12x,"stresses n1 at integration points:")
 9010 format(12x,"strains n1 at integration points:")
 9070 format(15x,i2,9f10.3)
 9075 format(15x,i2,6e14.6)
c
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine tanstf_copy_results               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/31/2013 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine tanstf_copy_results( kindex_to_copy,
     &                                outmat,
     &                                in3dmat, nrow, ncol, nz )
      implicit  none
      integer  kindex_to_copy, nrow, ncol, nz, i, j
#dbl      double precision
#sgl      real
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
c     *                 subroutine tanstf_make_avg                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 01/31/2013 rhd                 *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine tanstf_make_avg( nrows, ncols, matrix, averages )
      implicit  none
      integer nrows, ncols, i, j
#dbl      double precision
#sgl      real
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
