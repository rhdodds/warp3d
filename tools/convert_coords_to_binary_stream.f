c     ****************************************************************
c     *                                                              *
c     *   support program to make a stream or binary file of nodal   *
c     *   coordinates. used to speed up input processing for very    *
c     *   large models.                                              *
c     *                                                              *
c     *   a text file of coordinates for all nodes is read and       *
c     *   written in required format. file name must be coords.txt   *
c     *                                                              *
c     *   the binary format uses a segmented record type to support  *
c     *   extremely long record lengths for very large models.       *
c     *   file name will be: coords.stream or coords.binary
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/1/2013                  *
c     *                                                              *
c     ****************************************************************
c
      program make_coords_file
      implicit none
c
      integer :: filenum, node, inode, cfile, nonode, termout
      logical :: stream, binary
      character(len=100) :: filename
      character(len=20) :: answer
      double precision, allocatable :: xcoor(:), ycoor(:), zcoor(:)
      double precision :: x, y, z
c
      cfile   = 10
      termout = 6
c
c              coords.txt file format:
c
c              <number of nodes>
c              <binary or stream>  keyword
c              <node #> <x> <y> <z>
c              <node #> <x> <y> <z>
c               ...
c
c              nodes not required to be in order. data for all nodes must
c              be given.
c
c              coords.txt may contain comment lines anywhere that
c              start with c, C, !, # in column 1 with blanks in cols 2,3,4
c
      open( unit=cfile, file='coords.txt', status='old')
      call skip_comments( cfile, 6 )
      read(cfile,*) nonode
      call skip_comments( cfile, 6 )
      read(cfile,fmt="(a)") answer
      stream = index( answer, "stream" ) > 0
      binary = .not. stream
c
      allocate( xcoor(nonode), ycoor(nonode), zcoor(nonode) )
      do node = 1, nonode
       read(cfile,*) inode, x, y, z
       xcoor(inode) = x
       ycoor(inode) = y
       zcoor(inode) = z
      end do
c
      close(cfile)
      write(termout,*) '... coords.txt read. file closed ...'
c
      if( stream ) then
         open( unit=cfile, file='coords.stream', status='unknown',
     &         access="stream", form="unformatted" )
      else
         open( unit=cfile, file='coords.binary', status='unknown',
     &         access='sequential', form='unformatted',
     &         recordtype='segmented' )
      end if
c
      write(cfile) xcoor
      write(cfile) ycoor
      write(cfile) zcoor
c
      close(cfile)
c
      if( stream ) write(termout,*) '... coords.stream file written ...'
      if( binary ) write(termout,*) '... coords.binary file written ...'
      write(termout,*) " "
      write(termout,*) "... normal termination ..."
      write(termout,*) " "
      stop
      end










c     ****************************************************************
c     *                                                              *
c     *                 subroutine skip_comments                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/1/2013                  *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine skip_comments( fileno, kout )
      implicit none
c
      integer  :: fileno, kout
c
c               local variables
c               ---------------
c
      logical :: blank_line
      integer :: read_flg, i
      character :: line*40, first_char*1
c
c                skip any comment and blank lines
c                (c, C, #, ! in col 1 or line wit cols 1-4 blank)
c
      do
        read(fileno,9100,iostat=read_flg) line
        if( read_flg .ne. 0 ) then
           write(kout,9220)
           stop
        end if
        first_char = line(1:1)
        blank_line = .true.
        do i = 1, 40
          if( line(i:i) .ne. ' ' ) then
              blank_line = .false.
              exit
          end if
        end do
        if( blank_line ) cycle
        if( first_char .eq. 'c' .or. first_char .eq. 'C' .or.
     &      first_char .eq. '#' .or. first_char .eq. '!' ) cycle
         backspace fileno
         return
      end do
c
      return
c
 9100 format(a40)
 9220 format(/,'>>>> FATAL ERROR: skip_comments routine',
     & /,18x,'trying to read line. encountered end of file',
     & ' or read error',
     & /,18x,'before normal eof expected',
     & /,18x,'job aborted.'//)
c
      end subroutine
