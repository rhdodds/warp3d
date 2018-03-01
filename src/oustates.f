
c     ****************************************************************
c     *                                                              *
c     *                   subroutine oustates_files                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 2/2/2017 rhd               *
c     *                                                              *
c     *    writes patran or flat files of element state results      *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_files( ouflg, oubin, ouasc, ounod,
     &                           flat_file, stream_file, text_file,
     &                           compressed )
      use global_data ! old common.main
c
      use elem_block_data, only: history_blk_list ! read only
      use main_data, only: material_model_names,    ! all read only
     &                     output_states_type_opt1,
     &                     output_states_type_opt2
      implicit none
      logical :: ouflg, oubin, ouasc, ounod, flat_file, stream_file,
     &           text_file, compressed
c
c           locals
c
      integer, parameter ::  max_comment_lines=200
      logical :: local_debug, patran_file, loop_debug, found_cp_matl,
     &           last_blk
      integer, external :: warp3d_get_device_number, omp_get_thread_num,
     &                     warp3d_matl_num
      double precision, allocatable, dimension(:,:) :: st_values
      double precision ::
     &  start_time, end_time, zero, start_time2, end_time2
      double precision, external :: omp_get_wtime
      integer :: num_values, num_states, num_matl_models, ncrystals,
     &           matl_model_list(mxmat), mcount, num_comment_lines,
     &           cp_model_type, matnum, max_cp_states_values,
     &           max_states_values, nstates,
     &           allocate_error, count_num_cp_matls, max_num_crystals,
     &           blk, felem, elem_type, nterms_crystal_list,
     &           mat_type, i, k, nchar_matl_id, matl_count,
     &           warp3d_matl_model, last, header_file,
     &           now_thread, int_points, span, hist_size, size_state,
     &           data_type, patran_file_number, flat_file_number
      integer :: crystal_list(max_crystals) !in param_def
      integer, allocatable :: states_cp_matl_info(:,:),
     &                        elements_in_this_domain(:)
      real, allocatable :: sngl_values(:)
c
      character(len=:), allocatable :: matl_name_id
      character(len=80), dimension (:), allocatable :: comment_lines
      character(len=1)  :: dummy_string = " "
      character(len=8), allocatable :: state_labels(:)
      character(len=60), allocatable :: state_descriptors(:)
      character(len=100) :: header_name
      character(len=20) :: mnames(mxmat)
c
      data zero / 0.0d00 /
c
      call wmpi_alert_slaves ( 32 )
      call wmpi_bcast_log( ouflg )
      call wmpi_bcast_log( oubin )
      call wmpi_bcast_log( ouasc )
      call wmpi_bcast_log( ounod )
      call wmpi_bcast_log( flat_file )
      call wmpi_bcast_log( stream_file )
      call wmpi_bcast_log( text_file )
      call wmpi_bcast_log( compressed )
c
      call wmpi_bcast_int( output_states_type_opt1 )
      call wmpi_bcast_int( output_states_type_opt2 )
      call wmpi_bcast_int ( ltmstp ) ! load step number
c
      call wmpi_bcast_string( stname, 8 )
      call wmpi_bcast_string( lsldnm, 8 )


      local_debug = .false.
      if( local_debug ) then
          write(out,9000) myid, ouflg, oubin, ouasc, ounod,
     &    flat_file, stream_file, text_file, compressed
      end if
c
c      0.   check for not-implemented cases: nodal output
c
      if( ounod )  then
          write(out,9010)
          return
      end if
c
c      1.   count the number of different WARP3D (built-in & umat)
c           material models that will need a state
c           variable result file. compute other useful counts
c           for crystal plasticity as it has a more complex
c           states output design.
c
c           Ex. if 5 materials defined in FE in model have type mises,
c           and 3 materials are type CP, the mises and CP
c           material types appear only once in the material list.
c
c           for CP model, make an integer list of which
c           crystals will have results included in the states file

c
      max_cp_states_values     = -1
      max_states_values        = -1
      max_num_crystals         = -1
      count_num_cp_matls       = 0
      cp_model_type            = warp3d_matl_num( "CP",2)
      loop_debug               = .false.
      found_cp_matl            = .false.
      num_matl_models          = 0
      matl_model_list(1:mxmat) = 0
c
      call oustates_make_matl_lists ! sets size_state
      call oustates_make_crystal_list

c
c      2.   write a header file for each WARP3D material that appears
c           in state variable output. this file defines number of
c           output state variables and a label/descriptor for each
c           state variable. this will delete any existing header file
c           for the same model.
c
      nchar_matl_id = len( material_model_names(1) )
      allocate( character(len=nchar_matl_id) :: matl_name_id )
      allocate( comment_lines(max_comment_lines) )
      allocate( sngl_values(size_state), state_labels(size_state),
     &          state_descriptors(size_state) )
c
      call oustates_write_states_headers ! only on rank 0
c
c      3.   loop over the WARP3D (built-in/umat) material models
c           in FE model. output a material states file for
c           each type of material that appears in the model.
c           the material model name (e.g. mises) is included in the
c           file name. Omit models that do not support states output.
c
c           we're processing 1 load(time) step of results here. for mpi,
c           keep a list of elements in this domain just for convenience
c           at time to write states file
c
      start_time2 = omp_get_wtime()
      if( use_mpi ) then
         allocate( elements_in_this_domain(noelem) )
         elements_in_this_domain = 0
      end if
c
      do matl_count = 1, num_matl_models ! serial loop
c
c         3a.  set up values for this material defined in FE model
c
        warp3d_matl_model = matl_model_list(matl_count)
        matl_name_id = material_model_names(warp3d_matl_model)
        if( warp3d_matl_model .eq. cp_model_type ) then
           num_values = max_cp_states_values ! consistency w/ old code
        else
           call oustates_number_values( warp3d_matl_model, num_values,
     &                                  out )
        end if
        num_states = num_values ! consistency with old code
        if( num_states .eq. 0 ) cycle  ! material type no state output
        matl_name_id = material_model_names(warp3d_matl_model)
c
c         3b.  make array of states results w/ 1 column/element.
c              process all element blocks in parallel. handle allocate
c              here for clarity
c
        if( allocated( st_values ) ) then
         write(out,9110) 1
         call die_gracefully
        end if
        allocate( st_values(num_states,noelem), stat=allocate_error )
        if( allocate_error .ne. 0 ) then
          write(out,9110)
          call die_gracefully
        end if
        st_values = zero ! simplifies computational code
c
        call oustates_process_blks !
c
c         3c.  create the states result file, write results in requested
c              form, close file. release array of states values.
c
        call oustates_make_states_file
        deallocate( st_values )
c
      end do  ! over material count.
c
c      4.  all done with work. deallocates, return
c
      deallocate(  matl_name_id, comment_lines )
      deallocate( sngl_values, state_labels, state_descriptors )
      if( allocated( states_cp_matl_info ) )
     &   deallocate( states_cp_matl_info )
      if( use_mpi ) deallocate( elements_in_this_domain )
c
      end_time2 = omp_get_wtime()

      if( local_debug ) then
          write(out,*) '... wall time for states output: ',
     &       end_time2 - start_time2
          write(out,9060)
      end if
c
      return
c
 9000 format(//,2x,"... Entered oustates_files for states output",
     &  /,10x,"myid: ",i5,
     &  /,10x,"ouflg, oubin, ouasc, ounod:                    ",4l5,
     &  /,10x,"flat_file, stream_file, text_file, compressed: ",4l5  )
 9010 format(/1x,
     &'>>>>> Warning: material state variable output at nodes ',
     & 'not supported',
     & /,16x,'command ignored....'/)
 9020 format(/1x,
     &'>>>>> Error: material state variable output not yet available ',
     & 'for MPI-based execution',
     & /,14x,'job terminated....'/)
 9060 format(2x,"... Leaving oustates_files" )
 9110 format(/1x,
     &'>>>>> FATAL ERROR: routine oustates_files @ ', i2,
     & /,16x,'Job terminated....'/)
 9120 format(/1x,
     &'>>>>> Error: allocate failed @ 1 in oustates_files',2i7,i20,
     & /,16x,'Job terminated....'/)

c
      contains
c     ========
c

c     ****************************************************************
c     *                                                              *
c     *   (in contains) subroutine oustates_make_matl_lists          *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 5/28/2016  rhd             *
c     *                                                              *
c     *    make list of materials, material types, some counts, etc. *
c     *    note of CP materials are present                          *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_make_matl_lists
      implicit none
c
      do blk = 1, nelblk
        last_blk    = blk .eq. nelblk
        felem       = elblks(1,blk)
        elem_type   = iprops(1,felem)
        mat_type    = iprops(25,felem)
        matnum      = iprops(38,felem)
        if( loop_debug ) write(out,9070) blk, felem, elem_type, mat_type
        if( mat_type .eq. cp_model_type ) then
           found_cp_matl = .true.
           call oustates_set_cp_info
        end if
        call oustates_count_matl  !  add to list
        if( mat_type .ne. cp_model_type ) then
          call oustates_number_values( mat_type, nstates, out )
          max_states_values = max( max_states_values, nstates )
        end if
      end do
c
      max_states_values = max( max_states_values, max_cp_states_values )
      size_state = max_states_values
c
      if( local_debug ) then
         write(out,9080) num_matl_models
         do i = 1, num_matl_models
           k = matl_model_list(i)
           write(out,9090)  k, material_model_names(k)
         end do
         write(out,9240) max_states_values, max_cp_states_values
      end if
c
      return
c
 9070 format(10x,"block, felem, etype, mtype:  ",4i7 )
 9080 format(/,10x,"no. WARP3D material models used:  ",i7 )
 9090 format(15x,"model #, name: ",i3,2x,a20)
 9240 format(10x,"max_states_values, max_cp_states_values:  ",2i7)
c
      end subroutine oustates_make_matl_lists
c

c     ****************************************************************
c     *                                                              *
c     * (in contains)  subroutine oustates_count_matl                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/13/2014 (rhd)           *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_count_matl
c
c              locals
c
      integer :: i
c
c              add this WARP3D (built-in) material model
c              number to list of those used in this FE model
c
      if( num_matl_models .eq. 0 ) then
        num_matl_models = 1
        matl_model_list(1) = mat_type
        return
      end if
c
      do i = 1, num_matl_models
        if( mat_type .eq. matl_model_list(i) )  return
      end do
      num_matl_models = num_matl_models + 1
      if( num_matl_models .gt. mxmat ) then
         write(out,9000) mxmat
         call die_abort
      endif
      matl_model_list(num_matl_models) = mat_type
c
      return
c
 9000 format(/1x,
     &'>>>>> Error: too many WARP3D material models for state output ',
     & /,14x,'maximum supported: ',i5," at ou_get_states_count_matl",
     & /,14x,'job terminated....'/)
c
      end subroutine oustates_count_matl
c     ****************************************************************
c     *                                                              *
c     * (in contains)  subroutine oustates_set_cp_info               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                last modified :  6/7/2016 (rhd)               *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_set_cp_info
c
      use main_data, only: matprp, imatprp
      implicit none
c
c              we can have multiple CP materials defined in FE model and
c              they can have a varying number of crystals.
c
c              build an info array for materials that have the
c              CP model: material number, number of crystals.
c
c              this array may or may not prove useful.
c
c              update max number of crystals over all materials in
c              FE model for later use.
c
      integer :: ncrystals, info_vector(10)
      logical :: here_debug, found_in_list
c
      here_debug = .false.
c
      if( .not. allocated( states_cp_matl_info ) ) then
          allocate( states_cp_matl_info(mxmat,2) ) ! in param_def
          states_cp_matl_info = 0
          count_num_cp_matls  = 0
      end if
c
      ncrystals = imatprp(101,matnum) ! CP for an element blk
      max_num_crystals = max( max_num_crystals, ncrystals )
c
      if( count_num_cp_matls == 0 ) then ! first in list
        count_num_cp_matls = count_num_cp_matls + 1
        states_cp_matl_info(count_num_cp_matls,1) = matnum
        states_cp_matl_info(count_num_cp_matls,2) = ncrystals
      end if
c
      found_in_list = .false.
      do i = 1, count_num_cp_matls ! aready in list ?
       if( states_cp_matl_info(i,1) == matnum ) found_in_list = .true.
      end do
      if( .not. found_in_list ) then
        count_num_cp_matls = count_num_cp_matls + 1 ! add to list
        states_cp_matl_info(count_num_cp_matls,1) = matnum
        states_cp_matl_info(count_num_cp_matls,2) = ncrystals
      end if
c
      call mm10_set_state_sizes( info_vector, ncrystals,
     &                           output_states_type_opt1,
     &                           output_states_type_opt2, out   )
c
      max_cp_states_values = max( info_vector(1), max_cp_states_values )
c
      if( .not. last_blk ) return
c
      if( count_num_cp_matls .eq. 0 ) then
         write(out,9000) 1
         call die_abort
      end if
c
      if( here_debug ) then
         write(out,9010) count_num_cp_matls
         do i = 1, count_num_cp_matls
           write(out,9020) i, states_cp_matl_info(i,1),
     &                     states_cp_matl_info(i,2)
        end do
        write(out,9030) max_num_crystals
      end if
c
      return
c
 9000 format('>> FATAL ERROR: routine oustates_set_cp_info @ ',i2,
     & /,    '                job aborted',//)
 9010 format(/,5x,'... routine oustates_set_cp_info ...',
     & /, 12x,'count_num_cp_matls: ',i3,
     & /, 12x,'states_cp_matl_info array: ')
 9020 format(12x,3i6)
 9030 format(12x,'max number of crystals in FE model: ',i4 )
 9040 format('>> FATAL ERROR: routine oustates_set_cp_info @ ',i2,
     & /,    '                crystal plasticity materials must',
     & /,    '                all have the same states_output_level.',
     & /,    '                material # detected in conflict: ',i3,
     & /,    '                job terminated',//)
c
      end subroutine oustates_set_cp_info
c
c
c     ****************************************************************
c     *                                                              *
c     *   (in contains) subroutine oustates_make_crystal_list        *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/2/2016  rhd              *
c     *                                                              *
c     *    make list of crystals to be included in results for CP    *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_make_crystal_list
      implicit none

      integer :: i
c
      crystal_list(1:)    = 0
      nterms_crystal_list = 0
      if( .not. found_cp_matl ) return
c
c           max_num_crystals is largest number of crystals present for
c           all CP materials in this FE model. Must be <= max_crystal
c
      if( output_states_type_opt2 .eq. 0 ) then  ! all crystals
        do i = 1, max_num_crystals
          crystal_list(i) = i
        end do
        nterms_crystal_list = max_num_crystals
        return
      end if
c
c           there is a crystal list given in the output states
c           command. for now only 1 crystal number is allowed
c           but this could be expanded easily in the future with
c           an <integerlist> of crystal numbers.
c
      crystal_list(1) = output_states_type_opt2
      nterms_crystal_list = 1
c
      return
      end subroutine oustates_make_crystal_list

c     ****************************************************************
c     *                                                              *
c     * (in contains) subroutine oustates_count_warp3d_matls         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/13/2014 (rhd)           *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_count_warp3d_matls
      implicit none

      integer i
      mcount = 0
      do i = 1, mxmat
        if( material_model_names(i)(1:) .eq. "not_used" ) cycle
        mcount = mcount + 1
        mnames(mcount)(1:) =   material_model_names(i)(1:)
      end do
      return
      end subroutine oustates_count_warp3d_matls


c     ****************************************************************
c     *                                                              *
c     *  (in contains) subroutine oustates_write_states_headers      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 2/9/2017 rhd            *
c     *                                                              *
c     *    write the states header file for each type WARP3D         *
c     *    material appearing in current FE model                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_write_states_headers
      implicit none
c
c           write a header file for each WARP3D material that appears
c           in state variable output. this file defines number of
c           output state variables and a label/descriptor for each
c           state variable. this will delete any existing header file
c           for the same model.
c
c           crystal plasticity model processed as exception as it has
c           a more complex states output design and needs additional
c           data to set the number of states values in output and
c           which/how many labels to write into header.
c
      do matl_count = 1, num_matl_models
c
        warp3d_matl_model = matl_model_list(matl_count)
        num_comment_lines = 0
        if( warp3d_matl_model .eq. cp_model_type ) then
          call mm10_states_labels( max_num_crystals,
     &           output_states_type_opt1, output_states_type_opt2,
     &           max_cp_states_values, size_state, state_labels,
     &           state_descriptors, out, comment_lines,
     &           max_comment_lines, num_comment_lines,
     &           nterms_crystal_list, crystal_list  )
          num_states = max_cp_states_values ! consistency
        else
           call oustates_get_state_labels( warp3d_matl_model,
     &           size_state, num_states, state_labels,
     &           state_descriptors, out, comment_lines,
     &           max_comment_lines, num_comment_lines )
        end if
        if( num_states .eq. 0 ) cycle  ! material type no state output
c
        matl_name_id = adjustl(material_model_names(warp3d_matl_model))
        last = len_trim( matl_name_id )
        header_name(1:) = "states_header_" // matl_name_id(1:last)
        header_file = warp3d_get_device_number()
        if( header_file .eq. -1 ) then
             write(out,9110) 1
             call die_gracefully
        end if
c

        if( myid .ne. 0 ) cycle
        open( unit=header_file,file=header_name,access="sequential",
     &       form="formatted" )
        write(header_file,9100) matl_name_id
        if( num_comment_lines .gt. 0 ) then
           write(header_file,9210) " "
           do k = 1, num_comment_lines
             write(header_file,9210) trim(comment_lines(k))
           end do
           write(header_file,9210) " "
        end if
        write(header_file,9200) num_states
        do k = 1, num_states
           write(header_file,9220) k, state_labels(k),
     &                             state_descriptors(k)
        end do
        close(unit=header_file,status="keep")
      end do  ! matl_count
c
      return
c
 9100 format("#",/,"#  header file for state variable output",
     & /, "#  WARP3D material: ",a )
 9110 format(/1x,
     &'>>>>> FATAL ERROR: routine oustates_write_states_headers @ ', i2,
     & /,16x,'Job terminated....'/)
 9200 format("#",/,
     & "#  8 character state labels and longer descriptors",
     & /, "#  material model number, number of state variables ",
     & /,"#",/,2i6 )
 9210 format("#  ",a)
 9220 format(i8, 2x, a8, 2x, a )
c
      end subroutine oustates_write_states_headers

c
c     ****************************************************************
c     *                                                              *
c     *  (in contains) subroutine oustates_process_blks              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 8/22/2017 rhd              *
c     *                                                              *
c     *    drive over all element blocks in parallel. build          *
c     *    an array of states output for entire model.               *
c     *    this acts on all blocks that are a specified WARP3D       *
c     *    material type, e.g. mises, CP, ...                        *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_process_blks
      use main_data, only: imatprp, bar_types, link_types
      implicit none
c
      integer :: ncrystals_blk, i, elem
      logical :: is_bar_elem, is_link_elem
c
c           process all blocks. call a material routine to fill state
c           variables to output for a block. write out zero values for
c           elements not using current WARP3D material type.
c
c           all elements in the block have the same material type, e.g.,
c           bilinear, mises, CP, etc.
c
c           we are processing all blocks that have warp3d_matl_model
c
      call omp_set_dynamic( .false. )
      if( local_debug ) start_time = omp_get_wtime()
c
c$OMP PARALLEL DO  PRIVATE( blk, now_thread, felem, elem_type,
c$OMP&                      mat_type, int_points, span, hist_size,
c$OMP&                      ncrystals_blk, matnum )
c
      do blk = 1, nelblk
          if( elblks(2,blk) .ne. myid ) cycle
          now_thread  = omp_get_thread_num() + 1
          felem       = elblks(1,blk)
          elem_type   = iprops(1,felem)
          mat_type    = iprops(25,felem)
          int_points  = iprops(6,felem)
          span        = elblks(0,blk)
          hist_size   = history_blk_list(blk)
          ncrystals_blk = 0
          matnum      = iprops(38,felem)
c
          if( local_debug ) write(out,9050) myid, blk, felem, elem_type,
     &                      mat_type, int_points, span, hist_size
          is_bar_elem  = bar_types(elem_type)
          is_link_elem = link_types(elem_type)
          if( is_bar_elem .or. is_link_elem ) cycle
c
          if( mat_type .eq. warp3d_matl_model ) then
             if( warp3d_matl_model .eq. cp_model_type )
     &           ncrystals_blk = imatprp(101,matnum) ! this block
             call oustates_drive_mm( blk, st_values(1,felem),
     &                            num_states, mat_type, out,
     &                            iprops(1,felem), felem, span,
     &                            ncrystals_blk, max_num_crystals,
     &                            output_states_type_opt1,
     &                            output_states_type_opt2,
     &                            nterms_crystal_list, crystal_list )
             if( use_mpi ) then ! add blk elems to list for domain
              do i = 1, span
                elem = felem + i - 1
                elements_in_this_domain(elem) = 1
              end do
             end if
          end if
      end do ! over blks
c$OMP END PARALLEL DO
      if( local_debug ) end_time = omp_get_wtime()
c
      return
c
 9050 format(10x,"myid, block, felem, etype, mtype:  ",5i7,
     &  /,10x,   "int_pts, span, hist_size:    ",3i7 )
c
      end subroutine oustates_process_blks

c
c
c     ****************************************************************
c     *                                                              *
c     *  (in contains)  subroutine oustates_make_states_file         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/7/2016 rhd               *
c     *                                                              *
c     *    create, open, write close states file for a single        *
c     *    WARP3D material code type for 1 load (time) step          *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_make_states_file
      implicit none
c
c           set up the Patran file or the flat file.
c
c           Patran: open either binary or formatted files
c                   for output. assign names to
c                   the files based on type of output,
c                   step number (ltmstp) and optionally mpi rank
c
c           Flat: open either text or stream file
c                 assign name to the file based on type
c                 of output, step number (ltmstp) and optionally mpi rank
c
c           ouocst_elem creates the right file name and opens the
c           file with correct attributes. for Patran output, we send
c           the file number to it. for flat files, it gets a file
c           number, opens file and sends back the file number.
c           differences caused by a mix of old-new code.
c
      patran_file = oubin .or. ouasc
      data_type = 3   !  material states in elements output here
c
      patran_file_number = -1
      if( patran_file ) then
       patran_file_number = warp3d_get_device_number()
       if( patran_file_number .eq. -1 ) then
           write(out,9030)
           call die_abort
       end if
      end if
c
      call ouocst_elem( data_type, ltmstp, oubin, ouasc,
     &                  patran_file_number, 1, use_mpi, myid,
     &                  flat_file, stream_file, text_file, compressed,
     &                  flat_file_number, matl_name_id )
c
c           write the header info for each file type. the flat
c           stream file has only number of elements and
c           number of state values.
c
c           the write the entire array of results (all elements)
c
      if( patran_file ) call oustates_patran_header
      if( flat_file ) call oustates_flat_header( flat_file_number,
     &    text_file, stream_file, nonode, noelem, ltmstp, num_states,
     &    stname  )
c
      call oustates_wrt
c
c           close the states file
c
      call ouocst_elem( data_type, ltmstp, oubin, ouasc,
     &                  patran_file_number, 2, use_mpi, myid,
     &                  flat_file, stream_file, text_file, compressed,
     &                  flat_file_number, dummy_string )
c
      return
c
9030  format('>> WARNING: could not find a unit number to write',
     & /     '            patran or flat file. action skipped...',/)

      end subroutine oustates_make_states_file
c
c
c     ****************************************************************
c     *                                                              *
c     * (in contains) subroutine oustates_patran_header              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *              last modified : 12/15/2014 (rhd)                *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_patran_header
c
c                       local declarations
c
      character(len=4) :: title(80), title1(80)
      character(len=80) :: string, strng1, stepstring*7
      integer :: titl(80), titl1(80), i
      equivalence (title,titl), (title1,titl1)
c
      string = ' '
      strng1 = ' '
      string = 'element (material) state values '
      strng1 = 'loading '
c
      string(38:45) = stname
      strng1(9:16)  = lsldnm
      write(stepstring,'(i7)') ltmstp
      string(46:) = ', step '//stepstring
c
      do i = 1, 80
         title(i)  = string(i:i)
         title1(i) = strng1(i:i)
      end do
c
c                     write header records for each Patran
c                     file type.
c
      if( patran_file .and. oubin ) then
       write(patran_file_number) titl, num_states
       write(patran_file_number) titl1
       write(patran_file_number) titl1
       write(patran_file_number) noelem
      end if
c
      if( patran_file .and. ouasc ) then
       write(patran_file_number,900) titl
       write(patran_file_number,915) num_states
       write(patran_file_number,900) titl1
       write(patran_file_number,900) titl1
       write(patran_file_number,*) noelem
      end if
c
      return
c
 900  format(80a1)
 915  format(i5)
 916  format(i10)
c
      end subroutine oustates_patran_header

c     ****************************************************************
c     *                                                              *
c     * (in contains)  subroutine oustates_wrt                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *              last modified : 2/9/2017 rhd                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_wrt
c
c                       local declarations
c
      integer :: elem, k
      double precision, parameter :: small_tol=1.0d-30, zero=0.0d0
c
c           sanity check
c
      if( num_states .le. 0 ) then
        write(out,9000) 1
        call die_gracefully
      endif
c
c           zero v. small values to prevent 3-digit exponents
c           in formatted output.
c
c           Patran binary files are single precision
c
      where( abs( st_values ) .lt. small_tol ) st_values = zero
c
      if( patran_file .and. oubin ) then
         do elem = 1, noelem
           if( use_mpi ) then
              if( elements_in_this_domain(elem) .eq. 0 ) cycle
           end if
           do k = 1, num_states
             sngl_values(k) = sngl( st_values(k,elem) )
           end do
           write(patran_file_number) elem, 8,
     &                                sngl_values(1:num_states)
         end do
      end if
c
      if( patran_file .and. ouasc ) then
         do elem = 1, noelem
           if( use_mpi ) then
              if( elements_in_this_domain(elem) .eq. 0 ) cycle
           end if
           do k = 1, num_states
             sngl_values(k) = sngl( st_values(k,elem) )
           end do
           write(patran_file_number,920) elem, 8,
     &                                   sngl_values(1:num_states)
         end do
      end if
c
      if( flat_file .and. stream_file  ) then
         if( use_mpi ) then
           write(flat_file_number) num_states
           do elem = 1, noelem
            if( elements_in_this_domain(elem) .eq. 0 ) cycle
            write(flat_file_number) elem, st_values(1:num_states,elem)
           end do
         else
           write(flat_file_number) st_values ! big 2d array
         end if
      end if
c
      if( flat_file .and. text_file ) then
         if( num_states .gt. 25000 ) then
           write(out,9010)  num_states, 25000
           call die_gracefully
         end if
         if( use_mpi )write(flat_file_number,*) num_states
         do elem = 1, noelem
           if( use_mpi ) then
            if( elements_in_this_domain(elem) .eq. 0 ) cycle
            write(flat_file_number,9200) elem,
     &                  st_values(1:num_states,elem)
           else
            write(flat_file_number,9100) st_values(1:num_states,elem)
           end if
         end do
      end if
c
      return
c
 920  format(2i8,/,(6e13.6))
 9000 format('>> FATAL ERROR: routine oustates_wrt @ ',i2,
     &   /,  '                job terminated.',//)
 9010 format('>> FATAL ERROR: the number of states values per ',
     &  'element: ',i8,
     & /,    '                exceeds limit of: ',i8,
     &   /,  '                job terminated.',//)
 9100 format(25000e15.6)
 9200 format(i8,25000d15.6)
c
      end subroutine  oustates_wrt

      end subroutine oustates_files

c
c     ****************************************************************
c     *                                                              *
c     *             subroutine oustates_get_state_labels             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/3/2015 (rhd)             *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_get_state_labels( warp3d_matl, size_state,
     &             num_states, state_labels, state_descriptors, out,
     &             comment_lines, max_comment_lines, num_comment_lines )
      implicit none
c
c                       parameters
c
      integer :: warp3d_matl, size_state, num_states, out,
     &           max_comment_lines, num_comment_lines
      character(len=8)  :: state_labels(size_state)
      character(len=60) :: state_descriptors(size_state)
      character(len=80) :: comment_lines(max_comment_lines)
c
c                       locals
c
c
      select case ( warp3d_matl )
c
      case ( 1 )
        call mm01_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
      case ( 2 )
        call mm02_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
      case ( 3 )
        call mm03_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
      case ( 4 )
        call mm04_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
      case ( 5 )
        call mm05_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
      case ( 6 )
        call mm06_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
      case ( 7 )
        call mm07_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
      case ( 8 )  !  the UMAT in WARP3D
        call mm08_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
      case ( 9 )  !
        call mm09_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
      case ( 10 )   ! handled directly in other oustates routine
         write(out,9010) 1
         call die_gracefully
      case ( 11 )
        call mm11_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
      case default
           write(out,9000) warp3d_matl
           call die_abort
      end select
c
      return
c
 9000 format(/1x,
     &'>>>>> Error: material state output not supported for WARP3D ',
     & /,14x,'built-in model: ',i7,
     & /,14x,'job terminated' )
 9010 format('>> FATAL ERROR: routine oustates_get_state_labels @ ',i2,
     & /,    '                job aborted',//)
c
      end

c
c     ****************************************************************
c     *                                                              *
c     *             subroutine oustates_number_values                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 5/26/2016 rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_number_values( warp3d_matl, nstates, out )
      implicit none
c
c                       parameters
c
      integer :: warp3d_matl, nstates, out
c
c                       locals
c
      integer :: info_vector(10)

      select case ( warp3d_matl )

      case ( 1 )
         call mm01_set_sizes( info_vector )
         nstates = info_vector(4)
      case ( 2 )
         call mm02_set_sizes( info_vector )
         nstates = info_vector(4)
      case ( 3 )
         call mm03_set_sizes( info_vector )
         nstates = info_vector(4)
      case ( 4 )
         call mm04_set_sizes( info_vector )
         nstates = info_vector(4)
      case ( 5 )
         call mm05_set_sizes( info_vector )
         nstates = info_vector(4)
      case ( 6 )
         call mm06_set_sizes( info_vector )
         nstates = info_vector(4)
      case ( 7 )
         call mm07_set_sizes( info_vector )
         nstates = info_vector(4)
      case ( 8 )  !  abaqus umat
         call umat_set_features( info_vector )
         nstates = info_vector(4)
      case ( 9 )
         call mm09_set_sizes( info_vector )
         nstates = info_vector(4)
      case ( 10 ) ! handled directly by outher oustates routine
         write(out,9010) 1
         call die_gracefully
      case ( 11 )
         call mm11_set_sizes( info_vector )
         nstates = info_vector(4)
      case default
           write(out,9000) warp3d_matl
           call die_abort
      end select
c
      return
c
 9000 format(/1x,
     &'>>>>> Error: unknown WARP3D material model in ',
     & "ou_get_number_states: ", i4,
     & /,14x,'job terminated....'/)
 9010 format('>> FATAL ERROR: routine oustates_number_values @ ',i2,
     & /,    '                job aborted',//)
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine oustates_flat_header                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *           last modified : 12/16/2014 rhd                     *
c     *                                                              *
c     *     Write header lines for flat text files of state values   *
c     *     & sizes for stream files                                 *
c     *                                                              *
c     ****************************************************************

      subroutine oustates_flat_header(
     &  flat_file_number, text_file, stream_file, nonode, noelem,
     &  ltmstp, num_states, stname  )
      implicit none
c
c                       parameter declarations
c
      integer :: flat_file_number, nonode, noelem, ltmstp, num_states
      logical :: text_file, stream_file
      character(len=8) :: stname
c
c                       local declarations
c
      character(len=24) :: sdate_time_tmp
      integer :: step_num
c
c                       for now, we are trying results files
c                       without the header line
c
      if( stream_file ) then
c         write(flat_file_number) noelem, num_states
         return
      end if
c
      step_num = ltmstp
c
      write(flat_file_number,9000)
      write(flat_file_number,9014)
      write(flat_file_number,9020) stname
      write(flat_file_number,9030) nonode, noelem
      call fdate( sdate_time_tmp )
      write(flat_file_number,9040) sdate_time_tmp
      write(flat_file_number,9050) step_num
      write(flat_file_number,9000)
c
      return
c
 9000 format('#')
 9014 format('#  WARP3D element results: states')
 9020 format('#  Structure name: ',a8 )
 9030 format('#  Model nodes, elements: ',2i8)
 9040 format('#  ',a24)
 9050 format('#  Load(time) step: ',i8 )
 9060 format(2i8)
c
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine oustates_drive_mm                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 5/30/2016 rhd              *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_drive_mm( blk, elem_states_output,
     &                              num_states, warp3d_matl, iout,
     &                              iprops, felem, span,
     &                              ncrystals_blk, max_num_crystals,
     &                              output_states_type_opt1,
     &                              output_states_type_opt2,
     &                              nterms_crystal_list, crystal_list )
c
      use elem_block_data, only: history_blocks, history_blk_list
      use main_data, only: elems_to_blocks
      implicit none
c
c                       parameters
c
      integer :: blk, warp3d_matl, iout, num_states, iprops(*),
     &           felem, span, ncrystals_blk, max_num_crystals,
     &           output_states_type_opt1, output_states_type_opt2,
     &           nterms_crystal_list, crystal_list(*)
      double precision :: elem_states_output(num_states,*)
c
c                       locals
c
      integer :: nrow_states, blockno, hist_size,
     &           int_points, elnum
      logical do_block, do_a_block
c
      nrow_states = num_states

c              blk > 0 => this is the block number. do all elements
c                         in the block
c
c              blk < 0 => this is an element number. put state
c                           values into column 1 of results
      do_block = .true.
      if( blk > 0 ) then
         do_a_block = .true.
         blockno = blk
      else
         do_a_block = .false.
         elnum = -blk
         blockno = elems_to_blocks(elnum,1)
      end if
c
      int_points  = iprops(6)
      hist_size   = history_blk_list(blockno)
c
      select case ( warp3d_matl )
      case ( 1 )
        call mm01_states_values( blk, elem_states_output,
     &                           nrow_states, num_states )
      case ( 2 )
        call mm02_states_values( blk, elem_states_output,
     &                           nrow_states, num_states )
      case ( 3 )
        call mm03_states_values( blk, elem_states_output,
     &                           nrow_states, num_states )
      case ( 4 )
        call mm04_states_values( blk, elem_states_output,
     &                           nrow_states, num_states )
      case ( 5 )
        call mm05_states_values( blk, elem_states_output,
     &                           nrow_states, num_states )
      case ( 6 )
        call mm06_states_values( blk, elem_states_output,
     &                           nrow_states, num_states )
      case ( 7 )
        call mm07_states_values( blk, elem_states_output,
     &                           nrow_states, num_states )
      case ( 8 )
        call mm08_states_values( blk, elem_states_output,
     &                           nrow_states, num_states )
      case ( 9 )
        call mm09_states_values( blk, elem_states_output,
     &                           nrow_states, num_states )
      case ( 10 )
        call mm10_states_values( blk, elem_states_output,
     &                           num_states, iout,
     &                           iprops, felem, hist_size, int_points,
     &                           span, history_blocks(blockno)%ptr(1),
     &                           ncrystals_blk, max_num_crystals,
     &                           output_states_type_opt1,
     &                           output_states_type_opt2,
     &                           nterms_crystal_list, crystal_list  )
      case default
           write(iout,9000) warp3d_matl
           call die_abort
      end select
c
      return
c
 9000 format(/1x,
     &'>>>>> Error: WARP3D material model: ',i4,
     & /,14x,'does not support states output. should not'
     & /,14x,'invoked in oustates_drive_mm...',
     & /,14x,'job terminated....'/)
c
      end
