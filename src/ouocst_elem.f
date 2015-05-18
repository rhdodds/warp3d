c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouocst_elem                  *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 12/2/2014 rhd              *
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
c
      integer :: data_type, stepno, fileno, opt,  myid,
     &           flat_file_number 
      logical :: oubin, ouasc, use_mpi,
     &        flat_file, stream_file, text_file, compressed
      character (len=*) :: matl_name_id
c   
c                       locals
c
      character (len=80) :: bflnam, fflnam
      character (len=80), save :: flat_name
      character (len=4)  :: strtnm
      character (len=80) :: command
      integer :: result, now_len
      integer, external :: warp3d_get_device_number
!win      integer, external  :: system
      logical ::  patran_file
c
c                       branch on whether files are to be opened or  
c                       closed.
c
c                       data_type: 1 = strains, 2 = stresses, 
c                       3 = material states
c
      patran_file = oubin .or. ouasc
      select case( opt )
c
      case( 1 )   ! open files
c                       
c                       attach patran binary file
c                       
      if( patran_file .and. oubin ) then
         strtnm = 'webe'
         if( data_type .eq. 2 ) strtnm = 'webs'
         if( data_type .eq. 3 ) strtnm = 'webm'
         call ouflnm( strtnm, bflnam, stepno, use_mpi, myid )
         if( data_type .eq. 3 .and. matl_name_id(1:1) .ne. " " ) then
            bflnam(10:) = "_" // matl_name_id(1:)
         end if   
         open(fileno,file=bflnam,status='unknown',access='sequential',
     &               form='unformatted' )
      end if
c
c                       attach patran formatted file.
c                       
      if( patran_file .and. ouasc ) then
         strtnm = 'wefe'
         if( data_type .eq. 2 ) strtnm = 'wefs'
         if( data_type .eq. 3 ) strtnm = 'wefm'
         call ouflnm( strtnm, fflnam, stepno, use_mpi, myid )
         if( data_type .eq. 3 .and. matl_name_id(1:1) .ne. " " ) then
            fflnam(10:) = "_" // matl_name_id(1:)
         end if   
         open(fileno,file=fflnam,status='unknown',access='sequential',
     &               form='formatted' )
      end if
c
c                       attach flat text or stream file
c     
      if( .not. flat_file ) return
c      
      flat_file_number = warp3d_get_device_number()
      flat_name(1:) = ' '
      flat_name(1:3) = 'wee' 
      if( data_type .eq. 2 ) flat_name(1:3) ='wes'
      if( data_type .eq. 3 ) flat_name(1:3) ='wem'
      write(flat_name(4:),9000) stepno
c      
      if( flat_file .and. text_file ) then
         flat_name(9:) = '_text'
         if( use_mpi ) write(flat_name(14:),9010) myid
         if( data_type .eq. 3 .and. matl_name_id(1:1) .ne. " " ) then
            now_len = len_trim( flat_name )
            flat_name(now_len+1:) = "_" // matl_name_id(1:)
         end if   
         open(unit=flat_file_number,file=flat_name,
     &        status='unknown',access='sequential',
     &        form='formatted') 
      end if
c
      if( flat_file .and. stream_file ) then
         flat_name(9:) = '_stream'
         if( use_mpi ) write(flat_name(16:),9010) myid
         if( data_type .eq. 3 .and. matl_name_id(1:1) .ne. " " ) then
            now_len = len_trim( flat_name )
            flat_name(now_len+1:) = "_" // matl_name_id(1:)
         end if   
         open(unit=flat_file_number,file=flat_name,
     &        status='unknown',access='stream',
     &        form='unformatted') 
      end if
c
      case( 2 )  !  close files
c
c
      patran_file = oubin .or. ouasc
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
c
 9000 format( i5.5 )
 9010 format( i4.4 )
c 
      end

