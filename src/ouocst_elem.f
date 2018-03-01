c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouocst_elem                  *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 1/11/2017 rhd              *
c     *                                                              *
c     *     this subroutine opens or closes files for (1) Patran     *
c     *     binary or formatted output, or (2) flat file with        *
c     *     stream or text form                                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouocst_elem( data_type, stepno, oubin, ouasc, fileno,
     &                        opt, use_mpi, myid, flat_file,
     &                        stream_file, text_file, compressed,
     &                        flat_file_number, matl_name_id  )
      implicit none
c
c                       parameters
c
      integer :: data_type, stepno, fileno, opt,  myid,
     &           flat_file_number, dot_pos
      logical :: oubin, ouasc, use_mpi,
     &           flat_file, stream_file, text_file, compressed
      character(len=*) :: matl_name_id
c
c                       locals
c
      integer :: now_len
      logical :: ok, patran_file
      character(len=80) :: command
      character(len=80), save :: flat_name
c
c                       branch on whether files are to be opened or
c                       closed.
c
c                       data_type: 1 = strains, 2 = stresses,
c                                  3 = material states
c
c                       file_no is only for Patran results. calling
c                       routine must pass value.
c
c                       flat_file_number. routine here sets value
c                       on file open. calling routine must send
c                       value on close
c
      patran_file = oubin .or. ouasc
c
      select case( opt )
c
      case( 1 )   ! open files
        if( patran_file ) call ouocst_elem_pat_file
        if( flat_file )   call ouocst_elem_flat_file
        ok = patran_file .or. flat_file
        if( .not. ok ) then
          write(*,9000) opt
          call die_abort
        end if
c
      case( 2 )  !  close files
        ok = patran_file .or. flat_file
        if( .not. ok ) then
          write(*,9000) opt
          call die_abort
        end if
c
        if( patran_file ) then
          close( unit=fileno, status='keep' )
          return
        end if
c
        close( unit=flat_file_number, status='keep' )
        if( stream_file ) return
        if( compressed ) then
          command(1:) = ' '
          now_len = len_trim( flat_name )
          command(1:) = 'gzip ' // flat_name(1:now_len)
!win          result = system( command )
        end if
c
      end select

      return
 9000 format('>> FATAL ERROR: ouocst_elem. opt: ',i2,
     & /,    '                job aborted',// )
c
      contains
c     ========
c
c     ****************************************************************
c     *  (in contains)  subroutine ouocst_elem_pat_file              *
c     ****************************************************************
c
      subroutine ouocst_elem_pat_file
      implicit none
c
      integer :: now_len
      character :: patran_file_name*80, strtnm*4, form_type*20
c
c                       attach patran binary or ascii file
c
      if( oubin ) then
         strtnm = 'webe'
         if( data_type .eq. 2 ) strtnm = 'webs'
         if( data_type .eq. 3 ) strtnm = 'webm'
         form_type = "unformatted"
      else ! ascii
         strtnm = 'wefe'
         if( data_type .eq. 2 ) strtnm = 'wefs'
         if( data_type .eq. 3 ) strtnm = 'wefm'
         form_type = "formatted"
      end if
c
      patran_file_name = " "
      patran_file_name(1:4) = strtnm
      write(patran_file_name(5:),9000) stepno
c
      if( data_type .eq. 3 .and. matl_name_id(1:1) .ne. " " )
     &     patran_file_name(12:) = "_" // matl_name_id(1:)
c
      if( use_mpi ) then
         dot_pos = len_trim( patran_file_name ) + 1
         patran_file_name(dot_pos:dot_pos) = "."
         write(patran_file_name(dot_pos+1:),9010) myid
      end if
c
      open(unit=fileno,file=patran_file_name,status='unknown',
     &     access='sequential', form=form_type )
c
      return
c
 9000 format( i7.7 )
 9010 format( i4.4 )
c
      end subroutine ouocst_elem_pat_file
c
c
c     ****************************************************************
c     *  (in contains)  subroutine ouocst_elem_flat_file             *
c     ****************************************************************
c
      subroutine ouocst_elem_flat_file
      implicit none
c
      integer :: now_len, dot_pos
      integer, external :: warp3d_get_device_number
!win      integer, external  :: system
      character :: strtnm*4, form_type*20, access_type*20
c
      flat_file_number = warp3d_get_device_number()
      if( flat_file_number .eq. -1 ) then
        write(*,9110)
        call die_gracefully
      end if
c
      flat_name(1:) = ' '
      flat_name(1:3) = 'wee'
      if( data_type .eq. 2 ) flat_name(1:3) ='wes'
      if( data_type .eq. 3 ) flat_name(1:3) ='wem'
      write(flat_name(4:),9000) stepno
c
      if( flat_file .and. text_file ) then
         flat_name(11:) = '_text'
         if( data_type .eq. 3 .and. matl_name_id(1:1) .ne. " " ) then
            now_len = len_trim( flat_name )
            flat_name(now_len+1:) = "_" // matl_name_id(1:)
         end if
         if( use_mpi ) then
           dot_pos = len_trim( flat_name ) + 1
           flat_name(dot_pos:dot_pos) = "."
           write(flat_name(dot_pos+1:),9010) myid
         end if
         access_type = 'sequential'
         form_type   = 'formatted'
      end if
c
      if( flat_file .and. stream_file ) then
         flat_name(11:) = '_stream'
         if( data_type .eq. 3 .and. matl_name_id(1:1) .ne. " " ) then
            now_len = len_trim( flat_name )
            flat_name(now_len+1:) = "_" // matl_name_id(1:)
         end if
         if( use_mpi ) then
           dot_pos = len_trim( flat_name ) + 1
           flat_name(dot_pos:dot_pos) = "."
           write(flat_name(dot_pos+1:),9010) myid
         end if
         access_type = 'stream'
         form_type   = 'unformatted'
      end if

      open( unit=flat_file_number,file=flat_name, status='unknown',
     &     access=access_type, form=form_type )
c
      return
c
 9000 format( i7.7 )
 9010 format( i4.4 )
 9110 format(/1x,
     &'>>>>> FATAL ERROR: routine ouocst_elem_flat_file',
     & /,16x,'Job terminated....'/)
c
      end subroutine ouocst_elem_flat_file
      end subroutine ouocst_elem
