c     ****************************************************************
c     *                                                              *
c     *               subroutine mem_allocate                        *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *                last modified : 1/23/2015 rhd                 *
c     *                              : 02/07/17  mcm                 *
c     *                                                              *
c     *     this subroutine provides the general allocation/         *
c     *     deallocation of large arrays during problem solution.    * 
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine mem_allocate( itype )
c
      use main_data
      use contact, only : contact_cause, maxcontact, contact_force
c
      implicit integer (a-z)
      include 'common.main'
      double precision
     &  zero
      data zero / 0.0 /
c
c      write (*,*) myid,': =>in mem_allocate for number',itype
c
      go to ( 100, 200, 300, 400, 500, 600, 700, 800, 900,
     &        1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700,
     &        1800, 1900, 2000, 2100, 2200, 2300, 2400,
     &        2500 ), itype 
c
c              allocate and zero node temperature vectors.
c              totals and increments over a load step.
c
 100  continue
c
         if ( allocated( temper_nodes ) ) then
           write(out,9900)
           write(out,9930) 1
           call die_abort
           stop
         end if
         if ( allocated( dtemp_nodes ) ) then
           write(out,9900)
           write(out,9930) 2
           call die_abort
           stop
         end if
         if ( allocated( temper_nodes_ref ) ) then
           write(out,9900)
           write(out,9930) 0
           call die_abort
           stop
         end if
c
         allocate( temper_nodes(nonode), dtemp_nodes(nonode),
     &             temper_nodes_ref(nonode), stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910)
           call die_abort
           stop
         end if 
         do i = 1, nonode
          temper_nodes(i) = zero
          temper_nodes_ref(i) = zero
         end do
         return
c
c              allocate and zero element temperature vectors.
c              totals and increments over a load step.
c
 200  continue
c
         if ( allocated( temper_elems ) ) then
           write(out,9900)
           write(out,9930) 3
           call die_abort
           stop
         end if
         if ( allocated( dtemp_elems ) ) then
           write(out,9900)
           write(out,9930) 4
           call die_abort
           stop
         end if
c
         allocate(  temper_elems(noelem), dtemp_elems(noelem),
     &              stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9920)
           call die_abort
           stop
         endif 
         do i = 1, noelem
          temper_elems(i) = zero
         end do
         return
c
c              vectors of integer and logical data to support
c              constraints in non-global coordinates. vectors
c              are length number of structure nodes.
c
 300  continue
c
         if ( allocated( trnmat ) ) then
           write(out,9900)
           write(out,9930) 5
           call die_abort
           stop
         end if
         if ( allocated( trn ) ) then
           write(out,9900)
           write(out,9930) 6
           call die_abort
           stop
         end if
c
         allocate( trnmat(nonode), trn(nonode), stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9940)
           call die_abort
           stop
         end if 
         do i=1, nonode
            trn(i) = .false.
         enddo
         return
c
c
c              not used.
c
c
 400  continue
c
         return
c
c              integer arrays to store part of nodal load definitons
c              and the real blocked array for loading values (loddat)
c
 500  continue
c
         if ( allocated(node_load_defs) ) then
           deallocate( node_load_defs )
           write(out,9900)
           write(out,9930) 20
         end if
         allocate( node_load_defs(mxlc), stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9962)
           call die_abort
           stop
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
         if ( allocated(loddat_blocks) ) then
           deallocate( loddat_blocks )
           write(out,9900)
           write(out,9930) 21
         end if
         allocate( loddat_blocks(max_loddat_blks), stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9963)
           call die_abort
           stop
         end if 
         allocate(
     &       loddat_blocks(1)%block(1:mxndldcm,1:sizeof_loddat_blks) )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9963)
           call die_abort
           stop
         end if 
         loddat_blocks(1)%block(1:mxndldcm,1:sizeof_loddat_blks) = 0.0
         return
c
c              elstor table used during input only
c
 600  continue
c
         if ( allocated( elstor ) ) then
           deallocate( elstor )
           write(out,9900)
           write(out,9930) 9
         end if
c
         allocate( elstor(mxsepr,noelem), stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9970)
           call die_abort
           stop
         end if 

         if ( allocated( relstor ) ) then
           deallocate( relstor )
           write(out,9900)
           write(out,9930) 9
         end if
c
         allocate( relstor(mxsepr,noelem), stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9970)
           call die_abort
           stop
         end if 

         return
c
c              deallocate the elstor table once not needed.
c
 700  continue
         if ( allocated( elstor ) ) then
           deallocate( elstor )
         else
           write(out,9900)
           write(out,9980) 
         endif

         if ( allocated( relstor ) ) then
           deallocate( relstor )
         else
           write(out,9900)
           write(out,9980) 
         endif

         return
c
c              incmap and incid vector to start reading of input.
c
 800  continue
c
         if ( allocated( incmap ) ) then
           write(out,9900)
           write(out,9930) 10
           call die_abort
           stop
         end if
         if ( allocated( incid ) ) then
           write(out,9900)
           write(out,9930) 11
           call die_abort
           stop
         end if
c
         allocate( incmap(noelem), incid(mxndel*noelem),
     &             stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9990)
           call die_abort
           stop
         end if 
         return
c
c              incmap and incid vector for restart or to resize.
c
 900  continue
c
         if ( allocated( incmap ) ) then
           write(out,9900)
           write(out,9930) 10
           call die_abort
           stop
         end if
         if ( allocated( incid ) ) then
           write(out,9900)
           write(out,9930) 11
           call die_abort
           stop
         end if
c
         allocate( incmap(noelem), incid(inctop),
     &             stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9992)
           call die_abort
           stop
         end if 
         return
c
c              resize the incid vector.
c
 1000  continue
c
         deallocate( incid )
         allocate( incid(inctop), stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9994)
           call die_abort
           stop
         end if 
         return
c
c              available. deprecated allocate of diagonal stiffness
c
 1100  continue
       return
c
c              allocate the diagonal mass vectors.
c
 1200  continue
         if ( allocated( mdiag ) ) then
           write(out,9900)
           write(out,9930) 13
           call die_abort
           stop
         end if
c
         allocate( mdiag(3*nonode), stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           call die_abort
           stop
         end if 
         return
c
c              allocate the effective dynamic load, pbar
c
 1300  continue
         if ( allocated( pbar ) ) then
           write(out,9900)
           write(out,9930) 14
           call die_abort
           stop
         end if
c
         allocate( pbar(3*nonode), stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9820)
           call die_abort
           stop
         end if 
         return
c
c              allocate the vector for coordinate indexes
c
 1400  continue
         if ( allocated( crdmap ) ) then
           write(out,9900)
           write(out,9930) 15
           call die_abort
           stop
         end if
c
         allocate( crdmap(nonode), stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9830)
           call die_abort
           stop
         end if 
         return
c
c              allocate the vectors for constraint definitions
c
 1500  continue
         if ( allocated( cnstrn ) ) then
           write(out,9900)
           write(out,9930) 16
           call die_abort
           stop
         end if
         if ( allocated( cnstrn_in ) ) then
           write(out,9900)
           write(out,9930) 17
           call die_abort
           stop
         end if
c
         allocate( cnstrn(3*nonode), cnstrn_in(3*nonode),
     &             stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9840)
           call die_abort
           stop
         end if 
         return
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
 1600  continue
         if ( allocated( dload ) ) then
           write(out,9900)
           write(out,9930) 18
           call die_abort
           stop
         end if
         if ( allocated( rload ) ) then
           write(out,9900)
           write(out,9930) 19
           call die_abort
           stop
         end if
         if ( allocated( rload_nm1 ) ) then
           write(out,9900)
           write(out,9930) 22
           call die_abort
           stop
         end if
         if ( allocated(load_pattern_factors) ) then
           write(out,9900)
           write(out,9930) 23
           call die_abort
           stop
         end if
c
         allocate( dload(3*nonode), rload(3*nonode), 
     &             rload_nm1(3*nonode),
     &             stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9850)
           call die_abort
           stop
         end if 
         allocate( load_pattern_factors(mxlc,2), stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9850)
           call die_abort
           stop
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
         return
c
c              available. deprecated allocate the vectors used 
c                         only by the pcg solver.
c
 1700  continue
       return
c
c              deprecated allocate the local diagonal mass vector.
c     
 1800  continue
       return
c
c              allocate the data structure to store the definition
c              for nonlinear loading steps and initialize
c     
 1900  continue
       allocate( step_load_data(mxstep), stat = alloc_stat )
       if ( alloc_stat .ne. 0 ) then
          write(out,9900)
          write(out,9996)
          call die_abort
          stop
       end if
c
       step_load_data(1:mxstep)%num_load_patterns = 0
       return
c
c              allocate the array to store nodal values of material
c              properties for functionally graded materials
c     
 2000  continue
       allocate( fgm_node_values(nonode,fgm_node_values_cols),
     &           stat = alloc_stat )
       fgm_node_values_defined = .true.
       if ( alloc_stat .ne. 0 ) then
          write(out,9900)
          write(out,9998)
          call die_abort
          stop
       end if
       fgm_node_values(1:nonode,1:fgm_node_values_cols) = 0.0
       return
c
c              allocate the pointer-vector data structures to store
c              the equivalent nodal element forces for the load step.
c     
 2100  continue
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
       if ( alloc_stat .ne. 0 ) then
          write(out,9900)
          write(out,9902)
          call die_abort
          stop
       end if
       do i = 1, noelem
         nullify( elem_eq_loads(i)%forces )
         elem_eq_loads(i)%ncols = 0
       end do
       return
c
c              allocate the packed vector data structure to store
c              the equivalent nodal element forces for the load step.
c              compute required sizes based on existing pointer-vector
c              data structure
c     
 2200  continue
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
       return
c
c              allocate the packed vector data structure to store
c              the equivalent nodal element forces for the load step.
c              sizes are already known.
c     
 2300  continue
       if( .not. allocated( eq_node_force_indexes ) ) 
     &       allocate( eq_node_force_indexes(1:noelem) )
       eq_node_force_indexes(1:noelem) = 0
       if( allocated( eq_node_forces ) ) deallocate( eq_node_forces )
       if( eq_node_force_len .gt. 0 ) 
     &   allocate( eq_node_forces(1:eq_node_force_len) )
       return
c
c              allocate the one of two data structures for storage of 
c              rigid contact information during solution
c     
 2400  continue
       if( .not. allocated( contact_force ) ) then
         allocate( contact_force(nonode*mxndof), stat=alloc_stat )
       end if
       if ( alloc_stat .ne. 0 ) then
          write(out,9900)
          write(out,9904)
          call die_abort
          stop
       end if
       contact_force(1:nonode*mxndof) = zero
       return
c
c              allocate one of two data structures for storage of 
c              rigid contact information during solution
c     
 2500  continue
       if( .not. allocated( contact_cause ) ) then
         allocate( contact_cause(maxcontact,nonode), stat=alloc_stat )
       end if       
       if ( alloc_stat .ne. 0 ) then
          write(out,9900)
          write(out,9906)
          call die_abort
          stop
       end if
       contact_cause(1:maxcontact,1:nonode) = 0
       return
c
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
c
      end
