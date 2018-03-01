c     ****************************************************************
c     *                                                              *
c     *                      subroutine oudrive                      *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 6/3/2016 rhd               *
c     *                                                              *
c     *     drive output of any and all quantities requested by the  *
c     *     user. all phases of output except output of residual     *
c     *     loads are driven from this section of code.              *
c     *                                                              *
c     *     also handle the output commands file ...                 *
c     *                     output patran neutral  ...               *
c     *                     output model ...                         *
c     *                                                              *
c     *     this routine executes on root only                       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine oudrive( sbflg1, sbflg2, stname, ltmstp )
      use main_data, only: output_packets, output_states_type_opt1,
     &                                     output_states_type_opt2,
     &                                     windows_os
      implicit none
c
c                       parameter declarations
c
      character :: stname*8
      logical :: sbflg1, sbflg2
c
c                       local declarations
c
      real :: dumr
      double precision ::
     &   dumd
      character :: dums
      logical :: outyp, ouflg, ounod, wide, eform,
     &           oubin, ouasc, oupat, noheader, react_totals_only,
     &           compressed, ok, found, stream_file, text_file,
     &           flat_file
c
      logical, external :: matchs, endcrd, true, label, string, numi
      logical :: stress, prec, out_packet_now, neutral
      integer :: indev, outdev, idummy, jdummy, output_states_level,
     &           dum, osn, ltmstp, itype
c
      call iodevn( indev, outdev, idummy, jdummy )
c
c                       output commands ..... return
c
      if( matchs("commands",4) ) then
         call oudriv_cmds
         return
      end if
c
c                       output model in a flat file. return
c
      if( matchs( 'model', 4 ) ) then
        call oudriv_model_flat
        return
      end if
c
c                       output patran neutral file. return
c
      if( matchs( 'patran', 4 ) ) then
        neutral = matchs('neutral',4)
        if( .not. neutral ) then
          call backsp(2) ! must be patran results
          if( true(dumr) ) call splunj  ! a dummy routine
        else   !   outptut neutral file, leave
          call oudriv_patran_neutral
          return
        end if
      end if
c
c                       run possibly needed setups for special cases,
c                       i.e., output command just after restart
c                       w/o running thru lots of code that
c                       initizalizes variables.

      call compute_checks
c
c                       set defaults for regular output options
c
      wide              = .false.
      eform             = .false.
      oubin             = .false.
      ouasc             = .false.
      oupat             = .false.
      ounod             = .true.
      prec              = .false.
      noheader          = .false.
      react_totals_only = .false.
      out_packet_now    = .false.
      flat_file         = .false.
      stream_file       = .false.
      text_file         = .false.
      compressed        = .false.
c
      output_states_level = 1
c
c                       loop get all output qualifiers/options. can
c                       request only one of patran, flat or packets
c                       in an output command. fall out this loop
c                       when word is not an option/qualifier.
c
      found = .true.
c
      do while ( found )     !   over options before quantities
c
      if( matchs('wide',4) ) then
         wide = .true.
         cycle
      end if
c
      if( matchs('noheader',5) ) then
         noheader = .true.
         cycle
      end if
c
      if( matchs('totals',5) ) then
         if( matchs('only',4) ) call splunj
         react_totals_only = .true.
         noheader = .true.
         cycle
      end if
c
      if( matchs('eformat',5) ) then
         eform = .true.
         cycle
      end if
c
      if( matchs('precision',4) ) then
         prec = .true.
         cycle
      end if
c
      if( matchs('nodal',5) ) then
         ounod = .true.
         cycle
      end if
      if( matchs('node',4) ) then
         ounod = .true.
         cycle
      end if
c
      if( matchs('element',7) ) then
         ounod = .false.
         cycle
      end if
c
      if( matchs('flat',4) ) then
         flat_file = .true.
         if( matchs('stream',4) ) stream_file = .true.
         if( matchs('text',4) ) text_file = .true.
         if( matchs('compressed',4) ) compressed = .true.
         if( out_packet_now .or. oupat ) then
            call errmsg2( 32, dum, dums, dumr, dumd )
            go to 9999
         end if
         if( stream_file .and. text_file ) then
            call errmsg2( 82, dum, dums, dumr, dumd )
            go to 9999
         end if
         ok = stream_file .or. text_file
         if( .not. ok ) then
            call errmsg2( 83, dum, dums, dumr, dumd )
            go to 9999
         end if
         if( stream_file .and. compressed ) then
            call errmsg2( 84, dum, dums, dumr, dumd )
            compressed = .false.
         end if
         if( compressed .and. windows_os ) then
            call errmsg2( 86, dum, dums, dumr, dumd )
            compressed = .false.
         end if
         cycle
      end if
c
      if( matchs('packet',4) ) then
         if( output_packets ) then
            out_packet_now = .true.
         else
            call errmsg2( 28, dum, dums, dumr, dumd )
         end if
         if( out_packet_now .and. (oupat .or. flat_file) ) then
            call errmsg2( 32, dum, dums, dumr, dumd )
            go to 9999
         end if
         cycle
      end if
c
      if( matchs('patran',6) ) then
         if( matchs('binary',6) ) then
            oubin = .true.
            oupat = .true.
            if( matchs('format',6) ) ouasc = .true.
         else if( matchs('format',6) ) then
            ouasc = .true.
            oupat = .true.
            if( matchs('binary',6) ) oubin = .true.
         else
            ouasc = .true.
            oupat = .true.
         end if
         if( (out_packet_now .or.flat_file ) .and. oupat ) then
            call errmsg2( 32, dum, dums, dumr, dumd )
            go to 9999
         end if
         cycle
      end if
c
      found = .false.
      end do   ! do while over options at start of command
c
c
c    ===============  end of option/qualifier loop ==================
c
c                       all recognized 1st level options/qualifiers
c                       resolved. patran neutral file treated at
c                       very top
c
c                       get next output type on line, process and come
c                       back here, get another output quantity and
c                       continue until eol or error. on error,
c                       flush remainder of line and leave.
c
c                       ltmstp -> last completed load/time step
c                                 0 -> no computational results yet
c
c
      do ! over list of output quantities on the line

      osn = 0
      if( endcrd(dum) ) exit
c
      if( matchs('displacements',5) ) then
         osn   = 1
      else if( matchs('velocities',3) ) then
         osn   = 2
      else if( matchs('accelerations',3) ) then
         osn   = 3
      else if( matchs('stresses',6) ) then
         osn   = 4
      else if( matchs('strains',6) ) then
         osn   = 5
      else if( matchs('reactions',5) ) then
         osn   = 6
         if( matchs('totals',5) ) then
           if( matchs('only',4) ) call splunj
           react_totals_only = .true.
           noheader = .true.
         end if
      else if( matchs('temperatures',5) ) then
        osn = 8
      else if( matchs('states',6) ) then
        output_states_type_opt1 = 1
        output_states_type_opt2 = 0
        osn = 9
        if( matchs( 'level',5 ) ) call splunj
        if( matchs( 'type', 4 ) ) call splunj
        if( numi( output_states_type_opt1 ) ) call splunj
        if( matchs( 'crystal',4 ) ) call splunj
        if( matchs( 'option', 4 ) ) call splunj
        if( numi( output_states_type_opt2 ) ) call splunj
      end if
c
      if( ltmstp .eq. 0 ) then ! no results yet
        call scan_flushline
        call errmsg3( outdev, 14 )
        return
      end if
c
      select case ( osn )
c
c                       no recognizable output quantity found.
c                       error message, flush line, leave.
c
      case( 0 ) !  unrecognized quantity
        call errmsg( 138, dum, dums, dumr, dumd )
        call scan_flushline
        exit !  loop over quantities
c
c                       output nodal displacements, velocities,
c                       accelerations, reactions, temperatures
c
      case( 1, 2, 3, 6, 8 )
         itype = osn
         if( osn .eq. 6 ) itype = 4
         if( osn .eq. 8 ) itype = 5
         call oudva( itype, ouflg, oupat, oubin, ouasc, wide, eform,
     &               prec, noheader, react_totals_only,
     &               out_packet_now, flat_file, stream_file,
     &               text_file, compressed )
         if( ouflg ) then  ! parsing or list inside oustr failed
             call errmsg2( 87, dum, dums, dumr, dumd )
             call scan_flushline
             exit  !  loop over quantities
         end if
c
c                       stresses (=5), strains (=6)
c
      case( 4, 5 )
         if( osn .eq. 4 ) stress = .true.
         if( osn .eq. 5 ) stress = .false.
         if ( oupat .or. flat_file)  then
           call oustr_pat_flat_file( stress, ouflg, oubin, ouasc,
     &        ounod, flat_file, stream_file,
     &        text_file, compressed  )
         else
           call oustr( stress, ouflg, oupat, oubin, ouasc, ounod,
     &                 wide, eform, prec, noheader, out_packet_now )
           if( ouflg ) then  ! parsing or list inside oustr failed
             call errmsg2( 87, dum, dums, dumr, dumd )
             call scan_flushline
             exit  !  from loop over quantities
           end if
      end if
c
c                      -- case 7 not used. available --
c
c                      material states to a patran or flat file
c
      case( 9 )
         if ( oupat .or. flat_file )  then
           call  oustates_files( ouflg, oubin, ouasc, ounod,
     &        flat_file, stream_file, text_file, compressed )
         else
           call errmsg2( 85, dum, dums, dumr, dumd )
         end if

      end select
c
      end do   ! next type of quantity to output
c
c
 9999 sbflg1 = .false.
      sbflg2 = .false.
c
      return
c
      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *                 subroutine oudriv_patran_neutral             *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *                   last modified : 6/29/2014 rhd              *
c     *                                                              *
c     *     scan remainder of command and drive writing a Patran     *
c     *     neutral file for the model                               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine oudriv_patran_neutral
      implicit none

      integer :: nc, dummy
      character :: neut_name*80
      logical ::
     &        local_debug
c
c                       output patran neutral ... treat here for
c                       simplicity.
c
      local_debug = .false.
      if( local_debug ) write(outdev,*)
     &       '.., entered oudriv_patran_neutral'

      if ( matchs('file',4) ) call splunj
      neut_name = ' '    ! important to init string blank
      if ( label(dummy) ) then
          call entits( neut_name, nc )
      else if( string(dummy) ) then
          call entits( neut_name, nc )
      else
          nc = len(stname)
          call name_strip( stname, nc )
          neut_name(1:) = stname(1:nc) // '.neutral'
      end if
      call ouneut( neut_name )
      return
c
      end subroutine   oudriv_patran_neutral
c     ****************************************************************
c     *                                                              *
c     *                 subroutine oudriv_model_file                 *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *                   last modified : 9/13/2014 rhd              *
c     *                                                              *
c     *     scan remainder of command and drive writing a simple     *
c     *     flat file for model description                          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine oudriv_model_flat
      implicit none
c
      integer :: nc, dummy
      character :: flat_name*80
      logical ::
     &        local_debug, flat_file, stream_file, text_file,
     &        warp3d_convention, patran_convention
c
c                       output model flat ...
c
      local_debug = .false.
      if( local_debug ) write(outdev,*)
     &       '.., entered oudriv_model_flat'
c
      text_file         = .false.
      stream_file       = .false.
      compressed        = .false.
      patran_convention = .false.
      warp3d_convention = .false.
c
c                       must specify Patran or WARP3D element
c                       convention
c
      if( matchs('flat',4) ) call splunj ! dummy routine
      if( matchs('patran',4) ) then
        patran_convention = .true.
        if( matchs('convention',4) ) call splunj
      else if( matchs('warp3d',4) ) then
        warp3d_convention = .true.
        if( matchs('convention',4) ) call splunj
      else
        call errmsg2( 88, dum, dums, dumr, dumd )
        call scan_flushline
        return
      end if
c                       must specify text or stream. compressed is
c                       option for text
c
      if( matchs('text',4 ) ) then
        text_file = .true.
        if( matchs('compressed',4) ) compressed = .true.
      else if( matchs('stream', 4) ) then
        stream_file = .true.
      else
        call errmsg2( 89, dum, dums, dumr, dumd )
        call scan_flushline
        return
      end if
c                       must specify keyword file.
c                       if no file name, make one
c
      flat_name = " "
      if( matchs('file',4) ) then
        if ( label(dummy) ) then
          call entits( flat_name, nc )
        else if( string(dummy) ) then
          call entits( flat_name, nc )
        else
          nc = len(stname)
          call name_strip( stname, nc )
          flat_name(1:) = stname(1:nc) // '_model_flat'
        end if
      else
        call errmsg2( 90, dum, dums, dumr, dumd )
        call scan_flushline
        return
      end if

      call oumodel_flat( flat_name, text_file, stream_file,
     &                   compressed, warp3d_convention,
     &                   patran_convention )
c
      return
c
      end subroutine   oudriv_model_flat

      end subroutine   oudrive

c     ****************************************************************
c     *                                                              *
c     *                 subroutine oudriv_cmds                       *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *                   last modified : 2/10/2018 rhd              *
c     *                                                              *
c     *     scan store the file name for output commands file ...    *
c     *     get the integerlist of load steps and convert to a bit   *
c     *     map for very simple look up if output should be done     *
c     *     after a load step                                        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine oudriv_cmds
      use global_data ! old common.main
      use main_data, only: output_command_file, output_step_bitmap_list
      implicit none
c
c                      locals
c
      integer :: dummy, char, errnum, icn, iplist, nchar, istep, lenlst
      integer :: intlst(mxlsz)
      logical, external :: matchs, label, string
      logical :: ok
      character(len=80) :: bad_input
c
      output_command_file(1:) = " "
      if( allocated( output_step_bitmap_list ) )
     &      deallocate( output_step_bitmap_list )
c
c                      get the name of file that contains only
c                      output commands (and comments)
c                      splunj is a dummy routine so
c                      optimizers don't delete the call to matchs
c
c                      file must exit now or error
c
      if( matchs( "use", 3 ) ) call splunj
      if( matchs( "file", 4 ) ) call splunj
c
      if( label(dummy) ) then
        call entits( output_command_file(1:), nchar )
      elseif( string(dummy) ) then
        call entits( output_command_file(1:), nchar )
      else
        call entits( bad_input, nchar )
        write(out,9310) bad_input(1:nchar)
        num_error = num_error + 1
        call scan_flushline
        return
      end if
c
      inquire( file = output_command_file, exist = ok )
      if( .not. ok ) then
         write(out,9300) output_command_file
         output_command_file(1:) = " "
         num_error = num_error + 1
         call scan_flushline
         return
      end if
c
c                      must have keyword "steps" so we can detect
c                      the integer list
c
      if( matchs( "after", 3 ) ) call splunj
      if( matchs( "for",3 ) ) call splunj
      if( .not. matchs( "steps", 4 ) ) then
         write(out,9320)
         output_command_file(1:) = " "
         num_error = num_error + 1
         call scan_flushline
         return
      end if
c
c                      get the list of step numbers after which
c                      the file of output commands will be executed.
c                      a list of just keyword "all" is not ok since
c                      user may nothave given any step definitions yet
c
      call scan
      call trlist( intlst, mxlsz, 0, lenlst, errnum )
c
c                       branch on the return code from trlist. a
c                       value of 1 indicates no error. a value of
c                       2 indicates that the parse rules failed in
c                       the list. a value of 3 indicates that the
c                       list overflowed its maximum length of mxlsz.
c                       a value of 4 indicates that no list was found.
c                       in these last three cases, the illegal list
c                       will be ignored and a new compute command will
c                       be sought.
c
      bad_input(1:) = " "
      if( errnum .ne. 1 ) then
        if( errnum  == 2) then
          call entits( bad_input, nchar )
          write(out,9200) bad_input(1:)
        else if( errnum == 3 ) then
          write(out,9210)
        else if( errnum == 4) then
          call entits( bad_input, nchar )
          write(out,9220) bad_input(1:)
        end if
        num_error = num_error + 1
        output_command_file(1:) = " "
        call scan_flushline
        return
      end if
c
c                       list of steps found. store step numbers
c                       in a bitmap. use 30 bits/word. traverse
c                       the list to extract each step number.
c
      icn = 0; iplist = 1; ok = .true.
c
      do
         if( iplist .eq. 0 ) exit
         call trxlst( intlst, lenlst, iplist, icn, istep )
         if( istep .le. 0) then
            write(out,9100) istep
            ok = .false.
            exit
         end if
         call ouset_map_entry( istep )
      end do
c
      if( .not. ok ) then
         write(out,9330)
         deallocate( output_step_bitmap_list )
         output_command_file(1:) = " "
         num_error = num_error + 1
      end if
c
      return
c
 9100 format(/1x,'>>>>> error: step number in list is not valid: ',i7 )
 9200 format(/1x,'>>>>> error: cannot recognize required ',
     & 'list of step numbers',
     & /14x,'scanning: ',a,//)
 9210 format(/1x,'>>>>> error: list of step numbers contains too ',
     & 'many entries.',//)
 9220 format(/1x,'>>>>> error: no list of load steps found',
     & /14x,'scanning: ',a,//)
 9300 format(
     & /1x, '>>>>> error: the specified output commands',
     & /14x,'file does not exist. file: ',a80,//)
 9310 format(
     & /1x, '>>>>> error: no name for file of commands found.',
     & /14x,'scanning: ',a,//)
 9320 format(
     & /1x, '>>>>> error: required keyword steps not found.',//)
 9330 format(
     & /1x, '>>>>> error: this command contained errors',//)
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine ouset_map_entry                   *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *             last modified : 2/10/2018 rhd                    *
c     *                                                              *
c     *              set bit map entry for step                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouset_map_entry( step_no )
      use main_data, only: output_step_bitmap_list
      implicit none
c
      integer :: step_no
c
c                      locals
c
      integer :: word, bit
c
      call ouresize_step_bitmap( step_no )
      word = (step_no - 1 ) / 30 + 1
      bit  = step_no - ( word - 1 ) * 30 - 1
      output_step_bitmap_list(word) =
     &        ibset( output_step_bitmap_list(word), bit )
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine ouresize_step_bitmap              *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *                   last modified : 2/10/2018 rhd              *
c     *                                                              *
c     *     resize as needed the bitmap vector storing the list of   *
c     *     steps for the automatic output command feature.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouresize_step_bitmap( step_no )
      use main_data, only: output_step_bitmap_list
      implicit none
c
      integer :: step_no
c
c                      locals
c
      integer :: new_length, old_length, word
      integer, allocatable, dimension(:) :: new_bitmap
c
      if( .not. allocated( output_step_bitmap_list ) ) then
       allocate( output_step_bitmap_list(1) )
       output_step_bitmap_list(1) = 0
      end if
c
      old_length = size( output_step_bitmap_list )
      word = (step_no - 1 ) / 30 + 1
      if( word <= old_length ) return
c
      new_length = 2 * word ! double current req'd size
      allocate( new_bitmap(new_length) )
      new_bitmap = 0
      new_bitmap(1:old_length) =
     &                output_step_bitmap_list(1:old_length)
      call move_alloc( new_bitmap, output_step_bitmap_list )
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 function ouchk_map_entry                     *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *                   last modified : 2/10/2018 rhd              *
c     *                                                              *
c     *              check step output map for an entry              *
c     *                                                              *
c     ****************************************************************
c
c
      function ouchk_map_entry( step_no ) result( test )
      use main_data, only: output_step_bitmap_list
      implicit none
c
      integer :: step_no
      logical :: test
c
c                      locals
c
      integer :: map_length, word, bit
c
      test = .false.
      if( .not. allocated( output_step_bitmap_list ) ) return
c
      map_length = size( output_step_bitmap_list )
      word = (step_no - 1 ) / 30 + 1
      if( word > map_length ) return
c
      bit  = step_no - ( word - 1 ) * 30 - 1
      test = btest( output_step_bitmap_list(word), bit )
c
      return
      end


