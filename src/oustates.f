
c     ****************************************************************
c     *                                                              *
c     *                   subroutine oustates_files                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/28/15  rhd               *
c     *                                                              *
c     *    writes patran or flat files of element state results      *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_files( ouflg, oubin, ouasc, ounod,
     &    flat_file, stream_file, text_file, compressed )
c     
      use elem_block_data, only: history_blk_list ! read only
      use main_data, only: material_model_names ! read only
      implicit integer (a-z)
$add common.main
      logical :: ouflg, outbin, ouasc, ounod, flat_file, stream_file,
     &           text_file, compressed 
c
c           locals
c
      integer, parameter :: size_state=200, max_comment_lines=200
      logical :: local_debug, patran_file, loop_debug
      integer, external :: warp3d_get_device_number
#sgl      real, allocatable, dimension(:) ::
#dbl      double precision, allocatable, dimension(:,:) ::
     &  st_values
#dbl      double precision ::
#sgl      real ::
     &  start_time, end_time, zero
      integer :: num_values, num_states, num_matl_models, 
     &           matl_model_list(mxmat), mcount, num_comment_lines,
     &           allocate_error
      real sngl_values(size_state)
c
      character(len=:), allocatable :: matl_name_id
      character(len=80), dimension (:), allocatable :: comment_lines
      character(len=1)  :: dummy_string = " "
      character(len=8)  :: state_labels(size_state)
      character(len=60) :: state_descriptors(size_state)
      character(len=100) :: header_name
      character(len=20) :: mnames(mxmat)
      data zero / 0.0d00 /
c
      local_debug = .false.
      if( local_debug ) then
          write(out,9000) ouflg, oubin, ouasc, ounod,
     &    flat_file, stream_file, text_file, compressed 
      end if
c
c           Check for not-implemented cases: nodal output, mpi 
c
      if( ounod )  then
          write(out,9010)
          return
      end if
c
      if( use_mpi ) then
          write(out,9020)
          call die_gracefully
      end if
c
c           count the number of different WARP3D (built-in)
c           material models that will need a state 
c           variable result file. 
c
      loop_debug = .false.
      num_matl_models = 0; matl_model_list(1:mxmat) = 0
      do blk = 1,nelblk
        felem       = elblks(1,blk)
        elem_type   = iprops(1,felem)
        mat_type    = iprops(25,felem)
        if( loop_debug ) write(out,9070) blk, felem, elem_type,         
     &         mat_type
        call oustates_count_matl  !  add to list
      end do
c      
      if( local_debug ) then
         write(out,9080) num_matl_models
         do i = 1, num_matl_models
           k = matl_model_list(i)
           write(out,9090)  k, material_model_names(k)
         end do
      end if 
c
c           write a header file for each WARP3D material that appears
c           in state variable output. this file defines number of
c           output state variables and a label/descriptor for each
c           state variable. this will delete any existing header file
c           for the same model.
c
c           commented code in loop writes a list of all WARP3D material
c           models to each states header file. the current scheme
c           does not require this data.
c 
      nchar_matl_id = len( material_model_names(1) )
      allocate( character(len=nchar_matl_id) :: matl_name_id )   
      allocate( comment_lines(max_comment_lines) )   
c      
      do matl_count = 1, num_matl_models 
        warp3d_matl_model = matl_model_list(matl_count)
        num_comment_lines = 0
        call oustates_get_state_labels( warp3d_matl_model, size_state,
     &           num_states, state_labels, state_descriptors, out,
     &           comment_lines, max_comment_lines, num_comment_lines )  
        if( num_states .eq. 0 ) cycle  ! material type no state output
        matl_name_id = adjustl(material_model_names(warp3d_matl_model))
        last = len_trim( matl_name_id )
        header_name(1:) = "states_header_" // matl_name_id(1:last)
        header_file = warp3d_get_device_number()
        if( header_file .eq. -1 ) then
             write(out,*) " abort @ 1"
             call die_abort
        end if     
        open( unit=header_file,file=header_name,access="sequential",
     &       form="formatted" )
        write(header_file,9100) matl_name_id
c        call oustates_count_warp3d_matls
c        write(header_file,9190) mcount
c        do k = 1, mcount
c          write(header_file,9192) k, mnames(k)(1:)
c        end do  
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
c           loop over the WARP3D (built-in) material models appearing
c           in the FE model. output a material states file for
c           each type of material model the model type appears
c           in the file name. Omit models that do not support
c           state variable output.
c
      do matl_count = 1, num_matl_models 
c      
        warp3d_matl_model = matl_model_list(matl_count)
        call oustates_number_values( warp3d_matl_model, 
     &                               num_values, out )
        num_states = num_values
        if( num_states .eq. 0 ) cycle  ! material type no state output
        matl_name_id = material_model_names(warp3d_matl_model)
c
c           process all blocks. call a material
c           routine to send back a state variables to output for
c           a block. write out zero values for elements not using
c           current WARP3D material type.
c
        if( allocated( st_values ) ) then
          write(out,*) "<<< fatal. st_values already allocated"
        end if
        allocate( st_values(num_states,noelem), stat=allocate_error )
        if( allocate_err .ne. 0 ) then
          write(out,9110) num_states, noelem, allocate_err
          call die_gracefully
        end if
        st_values = zero
c
        call omp_set_dynamic( .false. )
        if( local_debug ) start_time = omp_get_wtime()
c
c$OMP PARALLEL DO  PRIVATE( blk, now_thread, felem, elem_type,
c$OMP&                     mat_type, int_points, span, hist_size )
c$OMP&            SHARED( nelblk, elblks,iprops, history_blk_list, l
c$OMP&                     local_debug, warp3d_matl_model, st_values,
c$OMP&                     num_states, out )
c
        do blk = 1,nelblk
          now_thread  = omp_get_thread_num() + 1
          felem       = elblks(1,blk)
          elem_type   = iprops(1,felem)
          mat_type    = iprops(25,felem)
          int_points  = iprops(6,felem)
          span        = elblks(0,blk)
          hist_size   = history_blk_list(blk)
          if( local_debug ) write(out,9050) blk, felem, elem_type,         
     &         mat_type, int_points, span, hist_size
          if( mat_type .eq. warp3d_matl_model ) 
     &       call oustates_drive_mm( blk, st_values(1,felem),
     &                            num_states, mat_type, out )
        end do ! over blks
c$OMP END PARALLEL DO
        if( local_debug ) end_time = omp_get_wtime()
c
c           set up the Patran file or the flat file.
c
c           Patran: open either binary or formatted files
c                   for output. assign names to 
c                   the files based on type of output,
c                   step number (ltmstp)
c
c           Flat: open either text or stream file
c                 assign name to the file based on type
c                 of output, step number (ltmstp).
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
        deallocate( st_values )
c      
c           close the states file      
c
        call ouocst_elem( data_type, ltmstp, oubin, ouasc, 
     &                  patran_file_number, 2, use_mpi, myid, 
     &                  flat_file, stream_file, text_file, compressed,
     &                  flat_file_number, dummy_string )
c

      end do  ! over material count. 
      deallocate(  matl_name_id, comment_lines )
      
      if( local_debug ) write(out,9060)
c      
      return
c
 9000 format(2x,"... Entered oustates_files for states output",
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
 9030 format('>> WARNING: could not find a unit number to write',
     & /     '            patran or flat file. action skipped...',/)
 9040 format(
     &'>>>>> Error: state variable output supported only ',
     & 'for crystal plasticity model',
     & /,14x,'Element does not use CP model: ',i7,
     & /,14x,'job terminated....'/)
 9050 format(10x,"block, felem, etype, mtype:  ",4i7,
     &  /,10x,   "int_pts, span, hist_size:    ",3i7 )
 9060 format(2x,"... Leaving oustates_files" )
 9070 format(10x,"block, felem, etype, mtype:  ",4i7 )
 9080 format(10x,"no. WARP3D material models used:  ",i7 )
 9090 format(15x,"model #, name: ",i3,2x,a20)
 9100 format("#",/,"#  header file for state variable output",
     & /, "#  WARP3D material: ",a )
 9110 format(/1x,
     &'>>>>> Error: allocate failed @ 1 in oustates_files',2i7,i20,
     & /,16x,'Job terminated....'/)
  9190 format("#",/,"# number of WARP3D model names and list",
     & /,"#",/,   i6 )
 9192 format(i6,2x,a20)
 9200 format("#",/,
     & "#  8 character state labels and longer descriptors",
     & /, "#  material model number, number of state variables ",
     & /,"#",/,2i6 )
 9210 format("#  ",a)    
 9220 format(i6, 2x, a8, 2x, a )
c
      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *             subroutine oustates_count_warp3d_matls           *
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
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine oustates_count_matl                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/13/2014 (rhd)           *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_count_matl
c
c                       locals
c
      integer :: i
c      
c                       add this WARP3D (built-in) material model 
c                       number to list of those used in this FE model
c
      if( num_matl_models .eq. 0 ) then
        num_matl_models = 1
        matl_model_list(1) = mat_type
        return
      end if
c      
      do i = 1, num_matl_models 
        if( mat_type .eq. matl_model_list(i) ) 	return
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
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine oustates_patran_header                *
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
#dbl      character(len=4) :: title(80), title1(80)
#sgl      character(len=8) :: title(80), title1(80)
      character(len=80) :: string, strng1, stepstring*6
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
      write(stepstring,'(i6)') ltmstp
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
      end if
c      
      if( patran_file .and. ouasc ) then
       write(patran_file_number,900) titl
       write(patran_file_number,915) num_states
       write(patran_file_number,900) titl1
       write(patran_file_number,900) titl1
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
c     *                   subroutine oustates_wrt                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *              last modified : 1/27/2015  (rhd)                *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_wrt
c
c                       local declarations
c
      integer :: elem, k  
#dbl      double precision :: small_tol, zero
#sgl      real :: small_tol, zero
      data small_tol, zero / 1.0d-80, 0.0d00 / 
c  
c           zero v. small values t prevent 3-digit exponents
c           in formatted output.
c
c           Patran binary files are single precision
c    
      where( abs( st_values ) .lt. small_tol ) st_values = zero
c      
      if( patran_file .and. oubin ) then
         do elem = 1, noelem
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
           do k = 1, num_states 
             sngl_values(k) = sngl( st_values(k,elem) )
           end do 
           write(patran_file_number,920) elem, 8,
     &                                   sngl_values(1:num_states) 
         end do
      end if
c
      if( flat_file .and. stream_file  ) then
         write(flat_file_number) st_values ! big 2d array
      end if
c      
      if( flat_file .and. text_file ) then
         do k = 1, noelem
           write(flat_file_number,9100) st_values(1:num_states,k)
         end do
      end if
c
      return
c
 920  format(2i8,/,(6e13.6))
 9100 format(200e15.6)
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
      case ( 10 )
        call mm10_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
      case ( 11 )  !  the UMAT in WARP3D
        call mm11_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
c      case ( 12 )  !  the UMAT in WARP3D
c        call mm12_states_labels( size_state,
c     &      num_states, state_labels, state_descriptors, out )
c      case ( 13 )  !  the UMAT in WARP3D
c        call mm13_states_labels( size_state,
c     &      num_states, state_labels, state_descriptors, out )
c      case ( 14 )  !  the UMAT in WARP3D
c        call mm14_states_labels( size_state,
c     &      num_states, state_labels, state_descriptors, out )
c      case ( 15 )  !  the UMAT in WARP3D
c        call mm15_states_labels( size_state,
c     &      num_states, state_labels, state_descriptors, out )
c      case default
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
c     
      end      

c
c     ****************************************************************
c     *                                                              *
c     *             subroutine oustates_number_values                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/14/2014 (rhd)           *
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
      case ( 10 )
         call mm10_set_sizes( info_vector )
         nstates = info_vector(4)
      case ( 11 )
         call mm11_set_sizes( info_vector )
         nstates = info_vector(4)
c      case ( 12 )
c         call mm12_set_sizes( info_vector )
c         nstates = info_vector(4)
c      case ( 13 )
c         call mm13_set_sizes( info_vector )
c         nstates = info_vector(4)
c      case ( 14 )
c         call mm14_set_sizes( info_vector )
c         nstates = info_vector(4)
c      case ( 15 )
c         call mm15_set_sizes( info_vector )
c         nstates = info_vector(4)
c      case default
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
c      write(flat_file_number,9060) noelem, num_states
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
c     *                   last modified : 1/27/2015 (rhd)            *
c     *                                                              *
c     ****************************************************************
c
      subroutine oustates_drive_mm( blk, elem_states_output, 
     &                              num_states, warp3d_matl, out )
      implicit none
c
c                       parameters
c
      integer :: blk, warp3d_matl, out, num_states
#dbl      double precision :: elem_states_output(num_states,*)
#sgl      real  :: elem_states_output(num_states,*)
c
c                       locals
c
      integer :: nrow_states
c
      nrow_states = num_states
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
     &                           nrow_states, num_states )
      case ( 11 )
        call mm11_states_values( blk, elem_states_output, 
     &                           nrow_states, num_states )
c      case ( 12 )
c        call mm12_states_values( blk, elem_states_output, 
c     &                           nrow_states, num_states )
c      case ( 13 )
c        call mm13_states_values( blk, elem_states_output, 
c     &                           nrow_states, num_states )
c      case ( 14 )
c        call mm14_states_values( blk, elem_states_output, 
c     &                           nrow_states, num_states )
c      case ( 15 )
c        call mm15_states_values( blk, elem_states_output, 
c     &                           nrow_states, num_states )
c      case default
           write(out,9000) warp3d_matl
           call die_abort
      end select 
c
      return
c
 9000 format(/1x,
     &'>>>>> Error: WARP3D material model: ',i4
     & /,14x,'does not support states output. should not'
     & /,14x,'invoked in oustates_drive_mm...'
     & /,14x,'job terminated....'/)
c     
      end subroutine    
