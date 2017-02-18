c     ****************************************************************
c     *                                                              *
c     *                      subroutine tanstf                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 02/17/2017 rhd             *
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
      implicit none
      include 'common.main'
c
c                       parameter dclarations
c
      logical :: first 
      integer :: now_step, now_iter
c      
c                       local declarations
c
      double precision ::
     &  zero, start_estiff, end_estiff
      double precision, external :: omp_get_wtime
      logical :: local_debug
      integer :: blk, now_thread
      integer, external :: omp_get_thread_num
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
c     *                   last modified : 8/16/2016 rhd              *
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
     &                              cohesive_ele_types,
     &                              linear_displ_ele_types,
     &                              adjust_constants_ele_types,
     &                              axisymm_ele_types,
     &                              nonlocal_analysis,
     &                              asymmetric_assembly,
     &                              dmatprp, imatprp,
     &                              temperatures_ref, 
     &                              fgm_node_values_defined
c
      use damage_data, only : dam_ptr, growth_by_kill
c 
      use contact, only : use_contact     
c
      implicit none
      include 'common.main'
c
c                       parameter declarations
c
      logical :: first
      integer :: blk, now_iter, now_step
c
c                       local declarations 
c
      include 'include_tan_ek'
      double precision ::
     &  zero, lambda(mxvl,3,3) ! on stack
      logical :: local_debug, geo_non_flg, bbar_flg, 
     &           symmetric_assembly, block_is_killable 
      integer :: felem, elem_type, int_order, mat_type, num_enodes,
     &           num_enode_dof, totdof, num_int_points, span, utsz,
     &           cohes_type, surface, matnum, nrow_ek, ispan, relem,
     &           element
c     
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
c      
      symmetric_assembly = .not. asymmetric_assembly
      nrow_ek = utsz
      if( asymmetric_assembly ) nrow_ek = totdof**2
      
c
c             See if we're actually an interface damaged CP material.
c             code commented until Mark resumes work on this model
c        
      local_work%is_inter_dmg = .false.
c     if( iprops(42,felem) .ne. -1 ) then
c        local_work%is_inter_dmg = .true.
c        local_work%inter_mat = iprops(42,felem)
c        local_work%macro_sz = imatprp(132, local_work%inter_mat)
c        local_work%cp_sz = imatprp(133, local_work%inter_mat)
c        tm = local_work%inter_mat
c
c        local_work%sv(1) = dmatprp(116, tm)
c        local_work%sv(2) = dmatprp(117, tm)
c        local_work%sv(3) = dmatprp(118, tm)
c
c        local_work%lv(1) = dmatprp(119, tm)
c        local_work%lv(2) = dmatprp(120, tm)
c       local_work%lv(3) = dmatprp(121, tm)
c
c        local_work%tv(1) = dmatprp(122, tm)
c        local_work%tv(2) = dmatprp(123, tm)
c        local_work%tv(3) = dmatprp(124, tm)
c      end if
c
      call chk_killed_blk( blk, local_work%killed_status_vec,
     &                     local_work%block_killed )
c
c             check if blk has all killed elements -- if so skip
c             all calculations. just zero element [k]s
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
!DIR$ LOOP COUNT MAX=128  
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
c     *                last modified : 1/9/2016 rhd                  *
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
      use elem_block_data, only:  estiff_blocks 
      use main_data, only: asymmetric_assembly
c
      implicit none
      include 'common.main'
c
      integer :: type
      
      integer :: iok, idummy, iout, dummy, blk, felem, num_enodes, 
     &           num_enode_dof, totdof, span, utsz, nterms


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
c     *                   last modified : 1/10/2016 rhd              *
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
      use elem_block_data, only:  history_blocks, 
     &                            rot_n1_blocks, history1_blocks,
     &                            eps_n1_blocks, urcs_n1_blocks,
     &                            history_blk_list
c
      implicit none
      include 'common.main'
      include 'include_tan_ek'
c
c           parameter declarations
c
      integer :: blk, span, felem, ngp, nnode, ndof, totdof, mat_type
      integer :: belinc(nnode,*)
      logical :: geonl
c
c           local declarations
c
      integer :: hist_size 
      logical :: local_debug
      double precision :: zero
      data zero, local_debug / 0.0d00, .false. /
c
c               1) rotation matrices at integration points for geonl
c                  to unrotate cuchy stresses into global
c               2) element histories (allocate the local storage block)
c                  only for certain material models (eg crystal plasticity)
c               3) unrotated cauchy stresses
c
c               History data (if needed):
c                 o The global blocks are sized(hist_size,ngp,span)
c                 o The local block is sized (span,hist_size,ngp).
c                This makes it possible to pass a 2-D array slice for
c                all elements of the block for a single gauss point.
c
      if( geonl )  call gastr( local_work%rot_blk_n1,
     &                         rot_n1_blocks(blk)%ptr(1),
     &                         ngp, 9, span )
c
      hist_size = history_blk_list(blk)
      local_work%hist_size_for_blk = hist_size
c
      call gastr( local_work%urcs_blk_n1, urcs_n1_blocks(blk)%ptr(1),
     &            ngp, nstrs, span )
c
c
c
      select case( mat_type )
c      ------------------------
c
      case( 1,2,3,4,5,6,7,8,9 )  ! nothing to do
c     ==========================
c
        continue
c
      case( 10, 11 )
c     ==================
c
c           gather history data at n+1. has the [D] matrices
c
         allocate( local_work%elem_hist1(span,hist_size,ngp),
     &             local_work%elem_hist(span,hist_size,ngp) )
         call dptstf_copy_history(
     &     local_work%elem_hist1(1,1,1), history1_blocks(blk)%ptr(1),
     &     ngp, hist_size, span )
c
      case default
          write(local_work%iout,*) '>>> invalid material model number'
          write(local_work%iout,*) '    in dptstf_blocks'
          call die_abort
      end select
c
      return
c
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine duptrans                         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/10/2016 rhd              *
c     *                                                              *
c     *     creates a separate copy of element                       *
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
      implicit none
      include 'common.main'
c
c           parameters
c
       integer :: span, felem
      double precision ::  trnmte(mxvl,mxedof,3)
c
c
      integer :: nnode, ndof, totdof, j, k, node, jj
c!DIR$ ASSUME_ALIGNED trnmte:64
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
!DIR$ LOOP COUNT MAX=128  
!DIR$ IVDEP
       do k = 1, span
         node = incid(incmap(felem+k-1) + j-1)
         if ( .not. trn(node) ) cycle
         jj = (j-1)*3
         trnmte(k,jj+1,1) = trnmat(node)%mat(1,1)
         trnmte(k,jj+1,2) = trnmat(node)%mat(1,2)
         trnmte(k,jj+1,3) = trnmat(node)%mat(1,3)
         trnmte(k,jj+2,1) = trnmat(node)%mat(2,1)
         trnmte(k,jj+2,2) = trnmat(node)%mat(2,2)
         trnmte(k,jj+2,3) = trnmat(node)%mat(2,3)
         trnmte(k,jj+3,1) = trnmat(node)%mat(3,1)
         trnmte(k,jj+3,2) = trnmat(node)%mat(3,2)
         trnmte(k,jj+3,3) = trnmat(node)%mat(3,3)
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
      implicit none
c
c               parameter declarations
c
      integer :: ngp, hist_size, span
      double precision ::
     & local_hist(span,hist_size,ngp),
     & global_hist(hist_size,ngp,span)
c
      integer :: k, j, i     
c!DIR$ ASSUME_ALIGNED local_hist:64, global_hist:64
c
c      if( ngp .ne. 8 ) then
        do k = 1, ngp
         do  j = 1, hist_size
!DIR$ LOOP COUNT MAX=128  
!DIR$ IVDEP
            do  i = 1, span
               local_hist(i,j,k) = global_hist(j,k,i)
            end do
         end do
        end do
        return
c      end if
c
c                number of gauss points = 8, unroll.
c
c     do  j = 1, hist_size
c@!DIR$ LOOP COUNT MAX=128  
c@!DIR$ IVDEP
c        do  i = 1, span
c            local_hist(i,j,1) = global_hist(j,1,i)
c            local_hist(i,j,2) = global_hist(j,2,i)
c            local_hist(i,j,3) = global_hist(j,3,i)
c            local_hist(i,j,4) = global_hist(j,4,i)
c            local_hist(i,j,5) = global_hist(j,5,i)
c            local_hist(i,j,6) = global_hist(j,6,i)
c            local_hist(i,j,7) = global_hist(j,7,i)
c            local_hist(i,j,8) = global_hist(j,8,i)
c        end do
c      end do
c
      return
      end


c     ****************************************************************
c     *                                                              *
c     *                   subroutine tanstf_allocate                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/10/2016 rhd              *
c     *                                                              *
c     *     allocate data structure in local_work for updating       *
c     *     element stiffnesses                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine tanstf_allocate( local_work )
      use segmental_curves, only : max_seg_points, max_seg_curves
      implicit none

      include 'common.main'
      include 'include_tan_ek'
c
      integer :: error
      double precision :: zero
      data zero / 0.0d00 /
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
      local_work%b_block = zero
c
      allocate( local_work%ue(mxvl,mxedof),
     1  local_work%due(mxvl,mxedof),
     2  local_work%urcs_blk_n1(mxvl,nstrs,mxgp),
     3  local_work%rot_blk_n1(mxvl,9,mxgp), stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 3
           call die_abort
      end if
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
      allocate( local_work%weights(mxgp), stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 6
           call die_abort
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
c                always allocate cohes_rot_block. it gets passed as
c                parameter even when block is not interface elements
c
      allocate( local_work%cohes_rot_block(mxvl,3,3), stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 9
           call die_abort
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
c     *                   last modified : 1/10/2016 rhd              *
c     *                                                              *
c     *     release data structure in local_work for updating        *
c     *     strains-stresses-internal forces.                        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine tanstf_deallocate( local_work )
      implicit none
      include 'common.main'
      include 'include_tan_ek'
c
      integer :: local_mt, error
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
     3  local_work%rot_blk_n1, stat=error )
       if( error .ne. 0 ) then
           write(out,9000) 3
           call die_abort
       end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 15"
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
      deallocate( local_work%weights, stat=error )
       if( error .ne. 0 ) then
           write(out,9000) 6
           call die_abort
       end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 55"
c      
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
      if( allocated(local_work%elem_hist1) )
     &    deallocate(local_work%elem_hist1, stat=error )
        if( error .ne. 0 ) then
           write(out,9000) 9
           call die_abort
        end if
c       
      if( local_debug ) write(out,*) "..tanstf_deall @ 75"
      if( allocated(local_work%elem_hist) )
     &      deallocate(local_work%elem_hist, stat=error )
       if( error .ne. 0 ) then
           write(out,9000) 10
           call die_abort
       end if
      if( local_debug ) write(out,*) "..tanstf_deall @ 80"
c
c                cohes_rot_block is always allocated. it gets passed as
c                parameter even when block is not interface elements
c
      deallocate( local_work%cohes_rot_block, stat=error )
      if( error .ne. 0 ) then
           write(out,9000) 13
           call die_abort
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
      implicit none
c
      integer :: n      
      double precision :: vec(n), zero
      data zero / 0.0d00 /
c!DIR$ ASSUME_ALIGNED vec:64
c
!DIR$ IVDEP
      vec = zero
c
      return
      end
