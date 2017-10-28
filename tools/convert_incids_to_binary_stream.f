c     ****************************************************************
c     *                                                              *
c     *   support program to make a stream or binary file of element *
c     *   incidences. used to speed up input processing for very     *
c     *   large models.                                              *
c     *                                                              *
c     *   a text file of incidences for all elements is read and     *
c     *   written in required format. the existing file name         *
c     *   must be incids.txt                                         *
c     *                                                              *
c     *   the binary format uses a segmented record type to support  *
c     *   extremely long record lengths for very large models.       *
c     *   file name will be: incids.stream or incids.binary          *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 10/28/2017                 *
c     *                                                              *
c     ****************************************************************
c
      program make_incids_file
      implicit none
c
      integer :: filenum, elem, ifile, noelem, termout, loc, i, nnode,
     &           ielem, einc(30)
      integer, allocatable :: incids(:)
      logical :: stream, binary
      character(len=100) :: filename
      character(len=20) :: answer
c
      ifile   = 10
      termout = 6
c
c          incids.txt file format:
c          ----------------------
c
c          <number of elements>
c          <binary or stream>  keyword
c          <element #> <number of nodes> <list of nodes>   ***
c          <element #> <number of nodes> <list of nodes>   ***
c           ...
c
c          elements MUST be in numerical order. data for all elements
c          must be given.
c
c          incids.txt may contain comment lines any where that
c          start with c, C, !, # in column 1 with blanks in cols 2,3,4
c
c         *** modify lines below marked with *** if you don't want to
c             include teh number of element nodes on each line
c
      open( unit=ifile, file='incids.txt', status='old')
      call skip_comments( ifile, termout )
      read(ifile,*) noelem
      call skip_comments( ifile, termout )
      read(ifile,fmt="(a)") answer
      stream = index( answer, "stream" ) > 0
      binary = .not. stream
c
      allocate( incids(noelem*27) )
      loc = 0
c
      do elem = 1, noelem
       call skip_comments( ifile, termout )
       read(ifile,*) ielem, nnode, einc(1:nnode)   !   ****
       do i = 1, nnode
          incids(loc+i) = einc(i)
       end do
       loc = loc + nnode
      end do
c
      close(ifile)
      write(termout,*) '... incids.txt read. file closed ...'
c
      if( stream ) then
         open( unit=ifile, file='incids.stream', status='unknown',
     &         access="stream", form="unformatted" )
      else
         open( unit=ifile, file='incids.binary', status='unknown',
     &         access='sequential', form='unformatted',
     &         recordtype='segmented' )
      end if
c
      write(ifile) loc
      write(ifile) incids(1:loc)
!      write(termout,*) incids(1:loc)
      write(termout,*) '... incid values writte to file: ', loc
      close(ifile)
c
      if( stream ) write(termout,*) '... incids.stream file written ...'
      if( binary ) write(termout,*) '... incids.binary file written ...'
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
