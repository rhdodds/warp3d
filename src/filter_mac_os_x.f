c ************************************************************************
c *                                                                      *
c *   WARP3D source file builder - Mac OS X on Intel processors          *
c *                                                                      *
c ************************************************************************
c
c     see file filter_dvf.f for introductory comments
c
      implicit integer ( a-z)
      character * 1  firstc, dollar, pound, exclam, part2, percent
      character * 3 machin, pgmsiz, part1, addcrd,precis,fmt,form,ptr
c     &              ,exmach
      character part3*16, part4*59, sfname*16
      character * 16 filtbl(10)
      character * 80 line
c      character * 80 outnam
      logical    debug, quotes, match, cray_ptr
      equivalence (firstc,line( 1:1 ))
      equivalence (part1 ,line( 2:4 ))
      equivalence (part2 ,line( 5:5 ))
      equivalence (part3 ,line( 6:21))
      equivalence (part4 ,line(22:80))
      data  addcrd, dollar, pound, exclam / 'ADD', '$', '#', '!' /
      data  percent, ptr / '%', 'PTR' /
c
c             Mac OS X is same as Linux 64 bits.Code mark ups have been
c             added for "mac" in each place there was a l64 markup.
c
      machin = 'mac'  
      pgmsiz = 'SZ2'
      precis = 'DBL'
      fmt    = 'FMT'
      cray_ptr = .false.
      form   = ' z8'
      quotes = .false.
c
      outfil = 6
      out    = 10
      stklev = 1
      nowfil = 5
      nowlin = 0
      sfname = ' '
      debug = .false.
      do 10  i = 1, 10
  10   filtbl(i) = ' '
c
C
C                       generate upper to lower case/lower to upper case
C                       mappings.
c
      call trmup      
c
c                       open the debug output file
C
      if ( debug ) then
            open( unit = out,file = 'expand.debug',
     &      status = 'unknown')
      endif
c
c
c                       top of reading loop.  read next line from
c                       unexpanded source file.  if first character
c                       is not a special character, just write out the
c                       line.
c
c
 9999 continue
      read(unit=nowfil,fmt=9006,end=700) line
      if ( debug ) write(out,9007) firstc, part1, part2, part3, part4
c
c
c                       right strip the line
c
c
      do 15  i=1,80
         j = 81 - i
         if ( line(j:j) .ne. ' ' ) go to 20
   15 continue
   20 linel = j
c
c
c                       if there are double quotes get rid of them 
c                       unless its cdc or harris
c
c
      if ( quotes ) then
         do 50  i=1,linel
           if( line(i:i) .eq.'"' ) line(i:i) = ''''
  50     continue
      endif
c
c
c                       look for a #, $, !, or %. if none just emit
c                       the line.
c
c
c
      if ( firstc .eq. percent ) go to 100
      if ( firstc .eq. dollar  ) go to 200
      if ( firstc .eq. pound   ) go to 300
      if ( firstc .eq. exclam  ) go to 400
c
c
c                       just an ordinary line!
c
c
      write(outfil,9006) line(1:linel)
      nowlin = nowlin + 1
      go to 9999
c
c                       % card found. this would be a %include
c                       line for apollo computers. convert
c                       the line to lower case then emit.
c			for decs and suns, remove %.
c
 100  continue
      call convlc( line )
      if (machin .eq. 'HPU') write(outfil,9006) line(1:linel)
      if (machin .eq. 'SUN') write(outfil,9016) line(2:linel)
      if (machin .eq. 'DEC') write(outfil,9016) line(2:linel)
      nowlin = nowlin + 1
      go to 9999
c
c                       $add card found.  pop down in stack.
 200  continue
      stklev = stklev + 1
      nowfil = stklev + 10
      call convlc( part3 )
      sfname(1:16) = part3
      if ( debug ) write(out,9008) stklev, part3
      filtbl(stklev) = part3
      open( unit=nowfil, file=sfname, iostat=istats, status = 'old')
      if ( istats .ne. 0 ) go to 210
      if ( debug ) write(out,9013)
      go to 9999
 210  write(outfil,*) '>>> could not open $ADD file ',sfname
      stop
c
c                       # card found with machine dependent or
c                       size dependent line.
c                       output line only if for current machine or
c                       size.
c
 300  continue
c
      if ( match(machin,part1)  .or.  
     &     match(pgmsiz,part1)  .or.
     &     match(precis,part1) ) then
             write(outfil,9012) part2, part3, part4
c
c
c                       check for a #fmt card. if so insert the 
c                       correct format string for as many iterations
c                       as there are ### groups.
c                       then write the line.
c
c
      else if ( match(fmt,part1) ) then
 305     k = index(line,'###')
         if ( k .ne. 0 ) then
            line(k:k+2) = form
            go to 305
         else
            write(outfil,9012) part2,part3,part4
         endif
c
c                       if the string found is 'ptr' and we are using
c                       cray style pointers, then write line
c
      else if (match (ptr,part1) .and. cray_ptr) then
             write(outfil,9012) part2, part3, part4
      endif
      go to 9999
c
c
c                       got an exclamation line.   first check if it
c                       is a 'PTR' or a machine type. see if we should 
c                       emit it or skip it.
c
c                       
 400  if (match (ptr,part1) ) then
            if (.not. cray_ptr) then
                 write(outfil,9012) part2, part3, part4
            endif
      else if( .not.match(machin,part1) ) then 
            write(outfil,9012) part2,part3,part4
      endif
      go to 9999
c
c                       end of file encountered.  close current file
c                       pop file stack.
c
 700  continue
      if ( stklev .eq. 1 ) go to 1000
      close(unit=nowfil)
      if ( debug )write(out,9010) stklev
      stklev = stklev - 1
      nowfil = stklev + 10
      if (stklev .eq. 1) nowfil = 5
      go to 9999
c
c                       end of file on the unexpanded source.
c
 1000 continue
      if ( debug ) write(out,9011)
 9006 format( a)
 9007 format(3x,'>>> Line Read - ',a1,a3,a1,a16,a59)
 9008 format(/,3x,'>>>>> Pop Down. Level = ',i3, 'File = ',a16 )
 9009 format(/,3x,'>>>>> ADD file opened' )
 9010 format(/,3x,'>>>>> Pop stack. File closed. Level = ',i3 )
 9011 format(/,3x,'>>>>> End of source file. Job done.' )
 9012 format( a1, a16, a59 )
 9013 format(/,3x,'>>>>> Just opened a new input file' )
 9014 format(/,3x,'>>>>> Could not open the source file' )
 9016 format(6x,a)
      end
c **********************************************************************
c *                                                                    *
c * match - match two strings.  Case is converted to lower case for the*
c *         compare.  Original strings are untouched                   *
c *                                                                    *
c **********************************************************************
      logical function match(texta,textb)
      character *256 upchar,lochar
      common/tranup/ upchar,lochar
c
c
      character *(*) texta, textb
      character * 1  tex1, tex2
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
      match = .true.
      return
      end
      subroutine trmup
c
c**********************************************************************
c                       subprogram to set up a table to map lower     *
c                               case letters to upper case            *
c                               independent of character set.         *
c                                                                     *
c           input:      none                                          *
c           output:     upchar = mapping vector for lower to upper    *
c                                                                     *
c                                                                     *
c                       loones is the set to be mapped to uppercase   *
c                       upones is the uppercase mapping for loones    *
c                                                                     *
c**********************************************************************
C
      character *256 upchar,lochar
      character *26  loones, upones
      integer subscr
      common/tranup/ upchar,lochar
      data loones/'abcdefghijklmnopqrstuvwxyz'/
      data upones/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data maxchr/256/
c
c
c                       first map all characters to themselves
c
c                               
      do 10 i = 0,maxchr-1    
       lochar(i+1:) = char(i)
   10  upchar(i+1:) = char(i)
c
c
c                       now re-map certain lowercase ones to upper case
c                       and vice versa
c
      length = len(loones)
c
      do 20 i = 1,length
        subscr = ichar(upones(i:)) + 1
        lochar(subscr:subscr) = loones(i:i)
        subscr = ichar(loones(i:)) +1
   20   upchar(subscr:subscr) = upones(i:i)
      return
      end
c
      subroutine convlc(text)
c
      implicit integer(a-z)
c
c**********************************************************************
c                                                                     *
c                       routine to translate any upper case           *
c                       characters to lower case equivalents defined  *
c                       via trmup.                                    *
c**********************************************************************
c
      character  text * (*)                                           
      character *256 upchar,lochar
      common/tranup/ upchar,lochar
c
      nchar = len(text)                       
      do 10 i = 1,nchar
        subscr    =  ichar(text(i:i)) + 1
  10    text(i:i) = lochar( subscr : subscr )
      return
      end
