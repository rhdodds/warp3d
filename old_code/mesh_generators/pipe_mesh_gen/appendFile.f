c----67--1---------2---------3---------4---------5---------6---------712
c
c File: appendFile.f    Written By Greg Thorwald 1997
c                       Modified on July 10 1998  /JF
c
c----67--1---------2---------3---------4---------5---------6---------712
c
c  Append a given filename string with the desired extension.
c  Check if an extension is already present, then remove and
c  append new extension. Search from right to left, since the
c  filname might have blanks in the filename string.
c
      subroutine appendFile(filename,extension,appendName)
      implicit none
      character*(*) filename,extension,appendName
      character*1   cha(63)
      integer       lastDot,lastSlash,lastBackslash,lastchar,trial,i
c
      data  cha /'0','1','2','3','4','5','6','7','8','9',
     &           'a','b','c','d','e','f','g','h','i','j',
     &           'k','l','m','n','o','p','q','r','s','t',
     &           'u','v','w','x','y','z',
     &           'A','B','C','D','E','F','G','H','I','J',
     &           'K','L','M','N','O','P','Q','R','S','T',
     &           'U','V','W','X','Y','Z','$'/
c
c
c  Look for the end of the string (a blank) or an extension
c  (a "." then the extension).
c
      lastDot       = index(filename,'.',.true.)
      lastSlash     = index(filename,'\',.true.)
      lastBackslash = index(filename,'/',.true.)
c
c If lastDot > 0, lastdot > lastslash and lastdot > lastbackslash,
c then an extension exists, replace that extension with the given
c new extension string. If the above conditions are not fulfilled
c no extension is present. Then search, starting from the "back" of
c filename, the first character in the character set defined by 
c {0-9,a-z,A-Z,$}. This is assumed to mark the end of the filename
c string.
c
      if ( (lastDot.gt.0) .and. (lastDot.gt.lastSlash) .and.
     &     (lastDot.gt.lastBackslash) ) then
        appendName = filename(1:(lastDot-1)) // extension
      else
        lastchar = 0
        do i=1, 63
           trial    = index(filename,cha(i),.true.)
           lastchar = max(lastchar,trial)
        enddo
        appendName = filename(1:lastchar) //  extension
      end if
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
