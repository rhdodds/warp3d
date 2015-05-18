c
c     This program reads standard input and writes it back out again to
c     standard output.  If an input line has the character "\" as the last
c     non-balnk character, then it is assumed to be a continuation 
c     character.  In this case the "\" is omitted and the next line is 
c     appeneded to the current line in the output.
c
c
      program strip_make
      implicit integer (a-z)
      character *200 line
c
 10   continue
      read (*,9000, end = 1000) line
      call last_char (line, last)
c
      if ( line(last:last) .eq. "\" ) then
c
         write (*,9100) line(1:last-1)
c
 20      continue
         read (*,9000, end = 1000) line
         call first_char (line, first)
         call last_char (line, last)
         if ( line(last:last) .eq. "\" ) then
            write (*,9100) line(first:last-1)
            goto 20
         else
            write (*,9200) line(first:last)
         endif
c
      else
         write (*,9200) line(1:last)
      endif
c
      goto 10
c
c         
 1000     continue
 9000     format (a200)
 9100     format (a,$)
 9200     format (a)
          end
c
c ================================================================
c
      subroutine first_char (string, first)
      implicit integer (a-z)
      character *200 string 
c
      first = 1
      do i=1, 200
         if ( string (i:i) .ne. " ") then
            first = i
            exit
         endif
      enddo
c
      return
      end
c
c ================================================================
c
      subroutine last_char (string, last)
      implicit integer (a-z)
      character *200 string 
c
      last = 1
      do i=200,1,-1
         if ( string (i:i) .ne. " ") then
            last = i
            exit
         endif
      enddo
c
      return
      end
         

