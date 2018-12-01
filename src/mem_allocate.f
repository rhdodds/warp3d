c     ****************************************************************
c     *                                                              *
c     *               subroutine mem_allocate                        *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *                last modified : 11/26/2018 rhd                *
c     *                                                              *
c     *     provides the general allocation/deallocation of arrays   *
c     *     during problem solution.                                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mem_allocate( itype )
      use global_data ! old common.main
c
      use main_data
      use contact, only : contact_cause, maxcontact, contact_force
c
      implicit none
c
      include 'mkl.fi'
      integer :: itype, i, alloc_stat, k
      double precision, parameter :: zero=0.0d0
      logical, parameter :: local_debug=.false.
c
      if( local_debug ) write(out,*) myid,
     &        ': =>in mem_allocate for number', itype
c
      select case( itype )
c
c              allocate and zero node temperature vectors.
c              totals and increments over a load step.
c
      case( 1 )
c
         if( allocated( temper_nodes ) ) then
           write(out,9900)
           write(out,9930) 1
           call die_abort
         end if
         if( allocated( dtemp_nodes ) ) then
           write(out,9900)
           write(out,9930) 2
           call die_abort
         end if
         if( allocated( temper_nodes_ref ) ) then
           write(out,9900)
           write(out,9930) 0
           call die_abort
         end if
c
         allocate( temper_nodes(nonode), dtemp_nodes(nonode),
     &             temper_nodes_ref(nonode), stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910)
           call die_abort
         end if
         do i = 1, nonode
          temper_nodes(i)     = zero
          temper_nodes_ref(i) = zero
         end do
c
c              allocate and zero element temperature vectors.
c              totals and increments over a load step.
c
      case( 2 )
c
         if( allocated( temper_elems ) ) then
           write(out,9900)
           write(out,9930) 3
           call die_abort
         end if
         if( allocated( dtemp_elems ) ) then
           write(out,9900)
           write(out,9930) 4
           call die_abort
         end if
c
         allocate( temper_elems(noelem), dtemp_elems(noelem),
     &              stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9920)
           call die_abort
         endif
         do i = 1, noelem
          temper_elems(i) = zero
         end do
c
c              vectors of integer and logical data to support
c              constraints in non-global coordinates. vectors
c              are length number of structure nodes. allocatables
c              in derived types deleted automatically
c
      case( 3 )
c
         if( allocated( trn ) ) deallocate( trn )
         if( allocated( trnmat ) ) deallocate( trnmat )
         allocate( trnmat(nonode), trn(nonode), stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9940)
           call die_abort
         end if
         trn = .false. ! all terms
c
c              double vectors based on number of nodes.
c
      case( 4 )
         k = 3 * nonode
         allocate( u(k), v(k), a(k), du(k), idu(k), load(k),
     &             res(k), ifv(k), c(k) )
         allocate( dstmap(nonode), cstmap(k) )
         do i = 1, k
          u(i)  = zero
          v(i) = zero
          a(i) = zero
          du(i) = zero
         end do
         max_mpc = 500 ! resized as needed in incon.f
         max_mpc_tied = 500 ! resized as needed in tied_mesh.f
c
c              integer arrays to store part of nodal load definitons
c              and the real blocked array for loading values (loddat)
c
      case( 5 )
c
         if( allocated(node_load_defs) ) then
           deallocate( node_load_defs )
           write(out,9900)
           write(out,9930) 20
         end if
         allocate( node_load_defs(mxlc), stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9962)
           call die_abort
         end if
         do i = 1, mxlc
          node_load_defs(i)%node_count = 0
          node_load_defs(i)%how_defined = 0
          node_load_defs(i)%user_file_name = ' '
          nullify( node_load_defs(i)%nodal_loads )
         end do
c
         num_loddat_blks    = 1
         sizeof_loddat_blks = 5000
         next_loddat_col    = 0
         max_loddat_blks    = 100000  !  rhd, 9/14/12. was 3000
         if( allocated(loddat_blocks) ) then
           deallocate( loddat_blocks )
           write(out,9900)
           write(out,9930) 21
         end if
         allocate( loddat_blocks(max_loddat_blks), stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9963)
           call die_abort
         end if
         allocate(
     &       loddat_blocks(1)%block(1:mxndldcm,1:sizeof_loddat_blks) )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9963)
           call die_abort
         end if
         loddat_blocks(1)%block(1:mxndldcm,1:sizeof_loddat_blks) = 0.0
c
c              elstor table used during input only
c
      case( 6 )
c
         if( allocated( elstor ) ) then
           deallocate( elstor )
           write(out,9900)
           write(out,9930) 9
         end if
c
         allocate( elstor(mxsepr,noelem), stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9970)
           call die_abort
         end if
c
c              deallocate the elstor table once not needed.
c
      case( 7 )
         if( allocated( elstor ) ) then
           deallocate( elstor )
         else
           write(out,9900)
           write(out,9980)
         endif
c
c              incmap and incid vector to start reading of input.
c
      case( 8 )
c
         if( allocated( incmap ) ) then
           write(out,9900)
           write(out,9930) 10
           call die_abort
         end if
         if( allocated( incid ) ) then
           write(out,9900)
           write(out,9930) 11
           call die_abort
         end if
c
         allocate( incmap(noelem), incid(mxndel*noelem),
     &             stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9990)
           call die_abort
         end if
c
c              incmap and incid vector for restart or to resize.
c
      case( 9 )
c
         if( allocated( incmap ) ) then
           write(out,9900)
           write(out,9930) 10
           call die_abort
         end if
         if( allocated( incid ) ) then
           write(out,9900)
           write(out,9930) 11
           call die_abort
         end if
c
         allocate( incmap(noelem), incid(inctop),
     &             stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9992)
           call die_abort
         end if
c
c              resize the incid vector.
c
      case( 10 )
c
         deallocate( incid )
         allocate( incid(inctop), stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9994)
           call die_abort
         end if
c
c              global elements table (props, iprops, lprops )
c              we use mkl_malloc to obtain a better aligned
c              array. WARP3D must have mkl working so we can
c              use mkl_malloc
c
      case( 11 )
         if( ptr_iprops .ne. 0 ) call mkl_free( ptr_iprops )
         ptr_iprops = mkl_malloc( int8(4*noelem*mxelpr), 64 )
         ptr_props  = ptr_iprops
         ptr_lprops = ptr_iprops
         mxel = noelem
c
c              allocate the diagonal mass vectors.
c
      case( 12 )
c
         if( allocated( mdiag ) ) then
           write(out,9900)
           write(out,9930) 13
           call die_abort
         end if
c
         allocate( mdiag(3*nonode), stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           call die_abort
         end if
c
c              allocate the effective dynamic load, pbar
c
      case( 13 )
c
         if( allocated( pbar ) ) then
           write(out,9900)
           write(out,9930) 14
           call die_abort
         end if
c
         allocate( pbar(3*nonode), stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9820)
           call die_abort
         end if
c
c              allocate the vector for coordinate indexes
c
      case( 14 )
c
         if( allocated( crdmap ) ) then
           write(out,9900)
           write(out,9930) 15
           call die_abort
         end if
c
         allocate( crdmap(nonode), stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9830)
           call die_abort
         end if
c
c              allocate the vectors for constraint definitions
c
      case( 15 )
c
         if( allocated( cnstrn ) ) then
           write(out,9900)
           write(out,9930) 16
           call die_abort
         end if
c
         if( allocated( cnstrn_in ) ) then
           write(out,9900)
           write(out,9930) 17
           call die_abort
         end if
c
         allocate( cnstrn(3*nonode), cnstrn_in(3*nonode),
     &             stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9840)
           call die_abort
         end if
c
c              allocate the vectors for the step incremental nodal loads
c              and the vectors of total applied nodal forces at n
c              and n-1 (these include only user applied nodal forces and
c              equivalent nodal force from specified element loads).
c              Allocate and initialize the total (col 1) and incremental
c              load pattern multipliers. Col 1 stores the accumulated,
c              actual multipliers for each pattern thru the current load
c              step (including any step size adjustments made by the
c              crack growth algorithms). Col 2 stores the pattern
c              increments for the current step as specified by the user
c              (no reduction factors requested by crack growth
c              algorithms).
c
      case( 16 )
c
         if( allocated( dload ) ) then
           write(out,9900)
           write(out,9930) 18
           call die_abort
         end if
         if( allocated( rload ) ) then
           write(out,9900)
           write(out,9930) 19
           call die_abort
         end if
         if( allocated( rload_nm1 ) ) then
           write(out,9900)
           write(out,9930) 22
           call die_abort
         end if
         if( allocated(load_pattern_factors) ) then
           write(out,9900)
           write(out,9930) 23
           call die_abort
         end if
c
         allocate( dload(3*nonode), rload(3*nonode),
     &             rload_nm1(3*nonode), stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9850)
           call die_abort
         end if
         allocate( load_pattern_factors(mxlc,2), stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9850)
           call die_abort
         end if
c
         do i = 1, 3*nonode
           rload_nm1(i) = zero
           dload(i) = zero
         end do
         do i = 1, mxlc
           load_pattern_factors(i,1) = zero
           load_pattern_factors(i,2) = zero
         end do
c
c              element blocking table.
c
      case( 17 )
       if( allocated( elblks ) ) deallocate( elblks )
       allocate( elblks(0:3,1:mxnmbl), stat = alloc_stat )
       if( alloc_stat .ne. 0 ) then
          write(out,9900)
          write(out,9908)
          call die_abort
       end if
c
c              available
c
      case( 18 )
c
c              allocate the data structure to store the definition
c              for nonlinear loading steps and initialize
c
      case( 19 )
c
       allocate( step_load_data(max_step_limit), stat = alloc_stat )
       if( alloc_stat .ne. 0 ) then
          write(out,9900)
          write(out,9996)
          call die_abort
       end if
c
       step_load_data(1:max_step_limit)%num_load_patterns = 0
c
c              allocate the array to store nodal values of material
c              properties for functionally graded materials
c
      case( 20 )
c
       allocate( fgm_node_values(nonode,fgm_node_values_cols),
     &           stat = alloc_stat )
       fgm_node_values_defined = .true.
       if( alloc_stat .ne. 0 ) then
          write(out,9900)
          write(out,9998)
          call die_abort
       end if
       fgm_node_values(1:nonode,1:fgm_node_values_cols) = 0.0  ! single
c
c              allocate the pointer-vector data structures to store
c              the equivalent nodal element forces for the load step.
c
      case( 21 )
c
       if( allocated(elem_eq_loads) ) then
         do i = 1, noelem
          if( associated( elem_eq_loads(i)%forces ) ) then
             deallocate( elem_eq_loads(i)%forces )
             nullify( elem_eq_loads(i)%forces )
             elem_eq_loads(i)%ncols = 0
          end if
         end do
         return
       end if
c
       allocate( elem_eq_loads(1:noelem), stat = alloc_stat )
       if( alloc_stat .ne. 0 ) then
          write(out,9900)
          write(out,9902)
          call die_abort
       end if
       do i = 1, noelem
         nullify( elem_eq_loads(i)%forces )
         elem_eq_loads(i)%ncols = 0
       end do
c
c              allocate the packed vector data structure to store
c              the equivalent nodal element forces for the load step.
c              compute required sizes based on existing pointer-vector
c              data structure
c
      case( 22 )
c
       if( .not. allocated( eq_node_force_indexes ) )
     &       allocate( eq_node_force_indexes(1:noelem) )
       eq_node_force_indexes(1:noelem) = 0
       eq_node_force_len = 0
       if( allocated( eq_node_forces ) ) deallocate( eq_node_forces )
       do i = 1, noelem
          eq_node_force_len = eq_node_force_len +
     &                        elem_eq_loads(i)%ncols
       end do
       eq_node_force_len = eq_node_force_len * 3
       if( eq_node_force_len .gt. 0 )
     &   allocate( eq_node_forces(1:eq_node_force_len) )
c
c              allocate the packed vector data structure to store
c              the equivalent nodal element forces for the load step.
c              sizes are already known.
c
      case( 23 )
c
       if( .not. allocated( eq_node_force_indexes ) )
     &       allocate( eq_node_force_indexes(1:noelem) )
       eq_node_force_indexes(1:noelem) = 0
       if( allocated( eq_node_forces ) ) deallocate( eq_node_forces )
       if( eq_node_force_len .gt. 0 )
     &   allocate( eq_node_forces(1:eq_node_force_len) )
c
c              allocate the one of two data structures for storage of
c              rigid contact information during solution
c
      case( 24 )
       if( .not. allocated( contact_force ) ) then
         allocate( contact_force(nonode*mxndof), stat=alloc_stat )
       end if
       if( alloc_stat .ne. 0 ) then
          write(out,9900)
          write(out,9904)
          call die_abort
       end if
       contact_force(1:nonode*mxndof) = zero
c
c              allocate one of two data structures for storage of
c              rigid contact information during solution
c
      case( 25 )
c
       if( .not. allocated( contact_cause ) ) then
         allocate( contact_cause(maxcontact,nonode), stat=alloc_stat )
       end if
       if( alloc_stat .ne. 0 ) then
          write(out,9900)
          write(out,9906)
          call die_abort
       end if
       contact_cause(1:maxcontact,1:nonode) = 0
c
      case default
c
       write(out,9900)
       write(out,9800)
       call die_abort

      end select
c
 9800 format('                 bad case select')
 9810 format('                 mdiag allocate failure')
 9820 format('                 pbar allocate failure')
 9830 format('                 crdmap allocate failure')
 9840 format('                 cnstrn allocate failure')
 9850 format('                 dload, rload, rload_nm1',
     & /,    '                 allocate failure')
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 nodal temperatures')
 9920 format('                 element temperatures')
 9930 format('                 array already allocated. type: ',
     &      i3 )
 9940 format('                 node constraint transformations')
 9960 format('                 nodlod tables')
 9962 format('                 node_load_defs table')
 9963 format('                 loddat_blocks table')
 9970 format('                 elstor tables')
 9980 format('                 elstor deallocate failure')
 9990 format('                 incid, incmap allocate failure')
 9992 format('                 incid, incmap re-allocate failure')
 9994 format('                 incid re-allocate failure')
 9996 format('                 step_load_data allocate failure')
 9998 format('                 fgm_node_values allocate failure')
 9902 format('                 element equiv. forces')
 9904 format('                 contact forces')
 9906 format('                 contact causes')
 9908 format('                 element blocking/domains table')
c
      end
