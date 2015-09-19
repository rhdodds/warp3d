c ************************************************************************
c *                                                                      *
c *   WARP3D source file builder - Mac OS X                              *
c *                                                                      *
c ************************************************************************
c
c
c
c     This program accepts Fortran source code files that contain
c     additional markup lines that indicate various dependencies and
c     size limitations.  It emits a new Fortran source code file
c     for compilation on a specific computer platform/operating system.
c
c     usage:   filter.exe < input_source  > output_source
c
c     Features
c
c     Program expands any $ADD <file> lines that maybe present. These
c     lines provides a platform independent method to "include" other
c     source code files. The $ADDs may be stacked 10 levels deep.
c
c     Filters based on certain characters appearing in column 1:
c
c        $add < file>
c
c        #<platform>  where <platform> is 3 characters, i.e.,
c            #mac,  #win,  #l64 or #lnx (last two for linux)
c          deprecated:  #dec,  #cry,  #sga, #r60, #hpi, #sun, #sgi
c
c          if <platform> matches variable machin, the #... is stripped
c          and the line emitted to filtered source file
c
c        #dbl or #sgl  lines to be included for double or single precision
c                      versions
c          if dbl or sgl matches variable precis, the #... is stripped
c          and the line emitted to filtered source file
c
c        !<platform>  means include this line (stripped)
c          if NOT building source for <platform> 
c
c        @...  strip @ and emit remainder of line to the filtered 
c              source file. @ often used ahead of !DIR$ ... compiler
c              directives. present version replaces appearance of ###
c              on such lines with the string in max_span
c
c        %...  not used at present. can be used in future. see routine
c              process_line
c
c     Above capabilities are readily extended with this updated version
c     of the filter program.
c
      call do_filter
      stop
      end
            
      subroutine do_filter
      implicit none
c      
      integer :: outfil, out, stklev, nowfil, nowlin, io
      character (len=256) :: upchar,lochar
      character (len=3)   :: machin, part1, precis, max_span
      character (len=16)  :: part3
      character (len=59)  :: part4
      character (len=16)  :: sfname
      character (len=16)  :: filtbl(10)
      character (len=80)  :: line
      character (len=150) :: error_filename
      logical :: debug, quotes
c
      machin = 'win'  !   Windows
      precis = 'DBL'
      quotes = .false.
      max_span = "128" !  for compiler directives
c
      outfil = 6
      out    = 10
      stklev = 1
      nowfil = 5
      nowlin = 0
      sfname = ' '
      debug = .false.
      filtbl(1:10) = ' '
c
c              generate upper to lower case/lower to upper case
c              mappings.
c
      call trmup      
c
c              open the debug output file
c
      if( debug ) open( unit = out,file = 
     &                  'expand.debug', status = 'unknown' )
c
c              loop to process each line of the source file. use
c              a stack to handle multiple level $add commands.
c      
      do 
       read(unit=nowfil,fmt=9006,iostat=io) line
       if( io == 0 ) then  ! read is ok
         call process_line
         cycle
       end if
       if( io > 0 ) then ! unknown file read error
         inquire(unit=nowfil,name=error_filename)
         write(out,9100) io, error_filename
         stop
       end if
       if( stklev == 1 ) then  ! eof on read. could be end of source
         if ( debug ) write(out,9011)
         return
       end if
       close(unit=nowfil)
       if( debug ) write(out,9010) stklev
       stklev = stklev - 1
       nowfil = stklev + 10
       if( stklev .eq. 1 ) nowfil = 5 ! back reading top level source
      end do
c       
 9006 format( a)
 9010 format(/,3x,'>>>>> Pop stack. File closed. Level = ',i3 )
 9011 format(/,3x,'>>>>> End of source file. Job done.' )
 9100 format(/,'....... FATAL ERROR .......',/,
     &       /,'iostat error: ',i6,' reading file: ',a)
 
      contains
      
c **********************************************************************
c *                                                                    *
c * process_line. handle the line and if and how to include it in the  *
c * filtered source code file                                          *
c *                                                                    *
c **********************************************************************
c
c
      subroutine process_line
      
      implicit none
      
      character (len=1) :: firstc, dollar, pound, exclam, part2, 
     &                     percent, atsign 
      character (len=3) :: addcrd
      integer :: linel, i, istats, pos
      data  addcrd, dollar, pound, exclam / 'ADD', '$', '#', '!' /
      data  percent, atsign / '%', '@' /

      firstc = line(1:1)
      part1  = line(2:4)
      part2  = line(5:5)
      part3  = line(6:21)
      part4  = line(22:80)
c      
      if ( debug ) write(out,9007) firstc, part1, part2, part3, part4
c
c              right strip the line.
c              replace double w/ single quotes if req'd
c
      linel = len_trim( line )
c
      if( quotes ) then
         do i = 1, linel
           if( line(i:i) .eq.'"' ) line(i:i) = ''''
         end do
      end if
c
c              look for a #, $, !, or %. if none just emit
c              the line.
c      
      select case( firstc )
c      
c              % line found. retained as future optional
c              filter sentinel. strip % and write line for now
c
      case( '%' )
        write(outfil,9006) line(2:linel)
        nowlin = nowlin + 1
c        
c                       $add line found.  pop down in stack to
c                       start reading the include file.
      case( '$' )
        stklev = stklev + 1
        nowfil = stklev + 10
        call convlc( part3 )
        sfname(1:16) = part3
        if( debug ) write(out,9008) stklev, part3
        filtbl(stklev) = part3
        open( unit=nowfil, file=sfname, iostat=istats, status = 'old')
        if( istats .ne. 0 ) then
           write(outfil,9100) sfname
           stop
        end if        
        if( debug ) write(out,9013)
c
c                       # card found with machine dependent or
c                       size dependent line.
c                       output line only if for current machine or
c                       size.
c
      case( '#' )
        if( match(machin,part1)  .or. match(precis,part1) ) 
     &      write(outfil,9012) part2, part3, part4
c
c                       ! line. include line unless on specified
c                       machine type. 
c                       
      case( '!' )
        if( .not. match(machin,part1) )
     &      write(outfil,9012) part2,part3,part4
c
c              @ sign starting line. a line with compiler directive
c              setting max value of span for key do loops. can
c              be used by compiler to assess approaches for 
c              optimization. replace ### on line with max allowed block
c              size. If no ###, just emit line w/o @ sign
c
      case( '@' )
        pos = index( line(1:), "###" ) 
        if( pos .gt. 0 ) line(pos:) = max_span
        write(outfil,9006) line(2:linel)
        nowlin = nowlin + 1
c
c              just an ordinary line
c
      case default
        write(outfil,9006) line(1:linel)
        nowlin = nowlin + 1
c
      end select 
c
      return
c
 9006 format( a)
 9007 format(3x,'>>> Line Read - ',a1,a3,a1,a16,a59)
 9008 format(/,3x,'>>>>> Pop Down. Level = ',i3, 'File = ',a16 )
 9009 format(/,3x,'>>>>> ADD file opened' )
 9012 format( a1, a16, a59 )
 9013 format(/,3x,'>>>>> Just opened a new input file' )
 9014 format(/,3x,'>>>>> Could not open the source file' )
 9016 format(6x,a)
 9100 format(/,'....... FATAL ERROR .......',/,
     &       /,'>>> could not open $ADD file: ',a )
c
      end subroutine   process_line    
      
                  
c **********************************************************************
c *                                                                    *
c * match - match two strings.  Case is converted to lower case for    *
c *         the compare.  Original strings are untouched               *
c *                                                                    *
c **********************************************************************
c
c
      logical function match( texta, textb )
c
      implicit none
      character (len=*) :: texta, textb
c            
      character (len=1) :: tex1, tex2
      integer :: len1, len2, nchar, i
c
      match = .false.
      len1 = len(texta)
      len2 = len(textb)
      nchar = min(len1,len2)
c   
      do i = 1, nchar
        tex1 = lochar( (ichar(texta(i:i)) +1): )
        tex2 = lochar( (ichar(textb(i:i)) +1): )
        if( tex1 .ne. tex2 ) return
      end do
c      
      match = .true.
      return
      end function match
c
c**********************************************************************
c                                                                     *
c      set up a table to map lower case letters to upper case         *
c      independent of character set.                                  *
c      output:     upchar = mapping vector for lower to upper         *
c                                                                     *
c      loones is the set to be mapped to uppercase                    *
c      upones is the uppercase mapping for loones                     *
c                                                                     *
c**********************************************************************
c
c      
      subroutine trmup
      implicit none
c      
      character (len=26)  ::  loones, upones
      integer :: subscr, maxchr, i, length
      data loones/'abcdefghijklmnopqrstuvwxyz'/
      data upones/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data maxchr/256/
c
c              first map all characters to themselves
c                               
      do i = 0, maxchr-1    
        lochar(i+1:) = char(i)
        upchar(i+1:) = char(i)
      end do 
c
c              now re-map certain lowercase ones to upper case
c              and vice versa
c
      length = len(loones)
c
      do i = 1, length
        subscr = ichar(upones(i:)) + 1
        lochar(subscr:subscr) = loones(i:i)
        subscr = ichar(loones(i:)) +1
        upchar(subscr:subscr) = upones(i:i)
      end do
c        
      return
      end subroutine trmup
c
c**********************************************************************
c                                                                     *
c        translate any upper case characters to lower case            *
c        equivalents defined  via trmup.                              *
c                                                                     *
c**********************************************************************
c
      subroutine convlc( text )
c
      implicit none
      character (len=*) :: text
      integer :: nchar, i, subscr
c
      nchar = len(text)                       
      do i = 1, nchar
        subscr    = ichar(text(i:i)) + 1
        text(i:i) = lochar( subscr : subscr )
      end do 
c      
      return
      end subroutine convlc
c      
      end subroutine do_filter

