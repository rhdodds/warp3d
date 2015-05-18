      program patint                                                    
c                                                                       
c **********************************************************************
c *            Convert Patran6 hardcopy file to Interleaf ASCII file   *
c *                                                                    *
c *                   By: M. Keppel 2/27/91  - patran5 version         *
c *                       B. Dodds  12/1/91  - added color support &   *
c *                                            formats for ileaf5      *
c *                       M. Keppel 12/16/91 - fixed bug for closing   *
c *                                            last polygon of row &   *
c *                                            models with hoels       *
c *                       R. Dodds   2/14/92 - fixes problems with     *
c *                                            text lines              *
c *                       R. Dodds   1/4/95  - fixes problems found in *
c *                                            port to IBM RISC 6000   *
c *                                                                    *
c *                  ------------- conversion to patran6 ------------  *
c *                                                                    *
c *                       A. Gullerud 7/9/96 - updated to use new      *
c *                                            patran6 hexadecimal     *
c *                                            format                  *
c *                       A. Gullerud 10/29/96 - updated to use patran6*
c *                                            polygon fills instead of*
c *                                            scanlines, use v7 for   *
c *                                            lines instead of v6, and*
c *                                            fixed color             *
c *                       A. Gullerud 2/3/97 - add .doc to end of file *
c *                                            name (ileaf6 needs it). *
c *                                            also remove all refs.   *
c *                                            to TPS-4 (ileaf4).      *
c *                                                                    *
c **********************************************************************
c *                                                                    *
c *     Procedure for Handling Scan Line Data                          *
c *                                                                    *
c *     Patran handles fringe plots in several ways.  In Patran 5,     *
c *     fringe plots are either output as a series of polygons and     *
c *     lines plotted from the point farthest from the viewer to       *
c *     the front (preferred), or as a series of thin rectangles which *
c *     sweep across the picture, defining a bitmap-like picture.      *
c *     These methods do not require the use of scanlines, which is    *
c *     the method used by Patran 2.5  However, the code for handling  *
c *     scanlines remains within this program in case it is needed in  *
c *     the future.                                                    *
c *                                                                    *
c *     Scan lines are translated into polygons in the following manner.
c *     Each scan line is pixinc pixels in height (pixinc is given in  *
c *     a type record 7). The polygon is started at the y coordinate   *
c *     y + pixinc / 2 (The top of the scan line) and extended straight*
c *     down to y (The middle of the scan line).  The polygon then     *
c *     connects to the middle of the next scan line and continues to  *
c *     connect to the middle of the next scan line until the number of*
c *     colors in the scan line changes.  When this occurs the polygon *
c *     is closed by extending straight down from y (The middle of the *
c *     scan line) to y - pixinc / 2 (The bottom of the scan line).    *
c *     The closed polygons are written to the file and new polygons   *
c *     are started.  The points in the polygon are stored in a linked *
c *                                                                    *
c *     Note: modifications for multiple, consecutive scan lines with  *
c *           with same Y coordinate have been added.                  *
c *                                                                    *
c **********************************************************************
c
      character*256 fname, string, string2, homed, timest*8, line,
     &     linetmp, dirname
      integer iin, iout, infil, outfil, intl_version, trgevn, str_len,
     &     num_char, str_len_tmp
      integer inrec(50000)
      integer head(200), tail(200), space(5000), ghead
      integer ncolor(200), ocolor(200), npoly, opoly, tot_len
      logical dopoly, tilde, ferr
      real x(5000), y(5000), xl(201)
      logical first, error
      real xy(4), xy1(2), xy2(2), xyc(2) 
      real magenta
      common /device/ xr, yr, ratior, pixinc, devx, devy
      integer counter 
c
c                 Initialize
c
      iin    = 5
      iout   = 6
      infil  = 20
      outfil = 21
      devx   = 11.0
      devy   = 8.5
      first  = .true.
      opoly  = 0
      maxlen2 = 256
      counter = 1
c
c                 Initialize garbage list
c
      do i = 1, 5000
          space(i) = i+1
      end do
      space(5000) = 0
      ghead = 1
c
c                 Print out header.
c
      write(iout,9100)
c
c                 Get Patran hardcopy file name. Open for input.
c
      ferr = .false.
 10   continue
      if (ferr) then
         write (iout,'(/," >> ERROR: file cannot be opened:",a,/)')
     &      fname
      endif
      ferr = .true.
      write(iout,9000)
     &  ' >> Patran Hardcopy File (default file is patran.hrd.01): '
      fname(1:) = ' '
      read(iin,9010) fname
      call stripf( fname, str_len, error )
      if ( error ) fname(1:) = 'patran.hrd.01'
      open(infil, file=fname, status='old', iostat=ierr )
      if ( ierr .ne. 0 ) then
          write(iout,*) '    >  Invalid Selection, Try Again...'
          go to 10
      end if
      write(iout,*) '>> Patran Hardcopy File Opened Ok'
c
c                 Get Version of Interleaf for Document Generation
c                 1) Ileaf 5.2 (B&W) or 2) Ileaf 5.2 (Color)
c
      write(iout,*) '>> Interleaf Document Type:'
      write(iout,*) '      1 -- Ileaf 5 or 6 (Black & White) <- default'
      write(iout,*) '      2 -- Ileaf 5 or 6 (Color)'
      write(iout,9000) ' >> Selection: '
      read(iin,*) intl_version  
      if ( intl_version .le. 0 .or. intl_version .gt. 2 )
     &                                   intl_version = 1
c
c                  Get name of directory to store the interleaf file
c
 20   continue
      tilde = .false.
      write(iout,9000)
     &  ' >> Directory for Ileaf file (default is ~/desktop5.0/ ): '
      dirname(1:) = ' '
      read(iin,9010) dirname
      call stripf( dirname, str_len, error )
      if ( error ) then
         dirname(1:) = '~/desktop5.0/'
         str_len = 13
      endif
c
c                       if no slash at end of name, put in a slash
c
      if (dirname(str_len:str_len).ne.'/') then 
	 dirname(str_len+1:str_len+1) = '/'
	 str_len = str_len + 1
      endif
c
c                       if ~ is specified, resolve home directory
c
      if (dirname(1:1) .eq.'~') then
	 tilde = .true.
         nchome = trgevn( 'HOME'//char(0), homed(1:) )
         if ( nchome .eq. 0 )  then
            write(iout,9022)
            stop
         end if                
      endif
c
c                 Form filename Patran_xx.xx.xx.doc where xx.xx.xx is the
c                 current time.  Open file for output on desktop.
c
      call time( timest )
      timest(3:3) = '.'
      timest(6:6) = '.'
      string(1:) = ' '
      write(string(1:),9020) timest(1:8)
      if (tilde) then
         fname(1:) = homed(1:nchome) // dirname(2:str_len) // 
     &         string(1:15) // ".doc"
         tot_len = nchome + (str_len - 1) + 15 + 4
      else
         fname(1:) = dirname(1:str_len) // string(1:15) // ".doc"
         tot_len = str_len + 15 + 4
      endif
      open(outfil, file=fname, status='new',
     &     recl=maxlen2, err=10,
     &     access='sequential', form='formatted', iostat=ierr)
      write(iout,9140) fname(1:tot_len)
c
c                 Initialize output file (Page & Color Definitions).
c                 Note: Color definitions account for Patran's weird
c                 ordering of colors.
c
      if ( intl_version .eq. 1 ) write(outfil,9111)
      if ( intl_version .eq. 2 ) write(outfil,9112)
c
c
c         ----------------
c         NOW CONVERT FILE
c         ----------------
c
c                 Discard first line
c
      read(infil, end=1100, fmt=9010) line
c
c                 Read a line from Patran hardcopy file; get record type
c
   25 continue
c
c        
c
      jptr = 1
      read(infil, end=1100, fmt=9010) line
      call stripf( line, str_len, error)
      if ( error ) then
         write (iout,*) ' >> ERROR: line is blank. Skip to next line.'
         go to 25
      endif
      call hex_to_int( line, jptr, 2, irec, str_len)
c
 50   continue
c
c                 Check for valid record type.
c
      if (irec .le. 0 .or. (irec .gt. 22 .and. .not. irec.eq.255)) 
     &     go to 25
c
      go to (1010, 1020, 1030, 1040, 1050, 1060, 1070, 1080, 1090,
     &       1100, 1110, 1120, 1130, 1140, 1150, 1160, 1170, 1180,
     &       1190, 1200, 1210, 1220  ), irec
c
c       if it didn't match goto above, must be text (irec = 255)
c
      goto 1900
c
c                 Record Type 1: Short Vector
c                     Generate black line from current position to xy1.
c                     Reset current position to xy1.
c
 1010 continue
      call hex_to_int( line, jptr, 4, ix, str_len)
      call hex_to_int( line, jptr, 4, iy, str_len)
      xy(1) = real(ix)
      xy(2) = real(iy)
      call map(xy, xy1)
      write(outfil,9030) counter, xyc, xy1, 0, 0
      counter = counter + 1
      xyc(1) = xy1(1)
      xyc(2) = xy1(2)
      go to 2000
c
c                 Record Type 2: Long Vector
c                     Generate black line from xy1 to xy2.
c                     Reset current position to xy2.
c
 1020 continue
      call hex_to_int( line, jptr, 4, ix, str_len)
      call hex_to_int( line, jptr, 4, iy, str_len)
      xy(1) = real(ix)
      xy(2) = real(iy)
      call map(xy, xy1)
      call hex_to_int( line, jptr, 4, ix, str_len)
      call hex_to_int( line, jptr, 4, iy, str_len)
      xy(1) = real(ix)
      xy(2) = real(iy)
      call map(xy, xy2)
      write(outfil,9030) counter, xy1, xy2, 0, 0
      counter = counter + 1
      xyc(1) = xy2(1)
      xyc(2) = xy2(2)
      go to 2000
c
c                 Record Type 3: Dot
c                     Generate zero length black line at xy1.
c                     Reset current position to xy1.
c
 1030 continue
      call hex_to_int( line, jptr, 4, ix, str_len)
      call hex_to_int( line, jptr, 4, iy, str_len)
      xy(1) = real(ix)
      xy(2) = real(iy)
      call map(xy, xy1)
      write(outfil,9030) counter, xy1, xy1, 0, 0
      counter = counter + 1
      xyc(1) = xy1(1)
      xyc(2) = xy1(2)
      go to 2000
c
c                 Record Type 4: Fill Region (Not used)
c
 1040 continue
      go to 2000
c
c                 Record Type 5: Color
c                     Get foreground color. Ignore background color.
c
 1050 continue
      call hex_to_int( line, jptr, 4, icolor, str_len)
      go to 2000
c
c                 Record Type 6: Dashed Vector (Not used)
c
 1060 continue
      go to 2000
c
c                 Record Type 7: Raster Specification
c                     Save raster information
c                     xr     - Width of picture in pixels
c                     yr     - Height of picture in pixels
c                     ratior - Ratio of width of pixel to height of pixel
c                     pixinc - Height of scan line in pixels.
c
 1070 continue
      call hex_to_int (line, jptr, 4, ix, str_len)
      call hex_to_int (line, jptr, 4, iy, str_len)
c
c                   NOTE: in patran5, xr and yr are switched in the 
c                          hardcopy file for the case of landscape
c                          page orientation.  Since all we will 
c                          output from patran is landscape, we switch
c                          it back here.
c
      yr = real(ix)
      xr = real(iy)
c
      call hex_to_real (line, jptr, ratior, str_len)
      call hex_to_real (line, jptr, pixinc, str_len)
c      
      realx = devy * xr * ratior / yr
      if (realx .gt. devx) then
          realy = devx * yr / xr / ratior
          devy = realy
      else
          devx = realx
      endif
      go to 2000
c
c                 Record Type 8: Frame Advance
c                     Close last scan line polygons and write them out.
c                     Put out end of frame information to close
c                     previous frame. Put out start of frame
c                     information for next frame. Skip close of 
c                     frame information for first frame.
c                     Reset number of polygons to 0.
c
 1080 continue
      yl = 2 * ya - yh
      do i = 1, opoly
c
c                 Close current polygons.
c
          j = ghead
          ghead = space(j)
          space(j) = head(i)
          head(i) = j
          x(j) = xl(i)
          y(j) = yl
          j = ghead
          ghead = space(j)
          space(tail(i)) = j
          space(j) = 0
          tail(i) = j
          x(j) = xl(i+1)
          y(j) = yl
c
c                 Write out polygon except for blanks (color = 17)
c
          if (ocolor(i) .ne. 17) then
              write(outfil,9040) counter, ocolor(i)
              write(outfil,9050) counter, counter
              n1 = head(i)
 1084         n2 = space(n1)
              write(outfil,9060) counter, x(n1), y(n1), x(n2), y(n2), 
     &             0, 127
              counter = counter + 1
              if (n2 .ne. tail(i)) then
                  n1 = n2
                  go to 1084
              endif
              n1 = tail(i)
              n2 = head(i)
              write(outfil,9060) counter, x(n1), y(n1), x(n2), y(n2), 
     &             0, 127
              counter = counter + 1
              write(outfil,9070)
          endif
c
c                 Free up polygon list
c
          space(tail(i)) = ghead
          ghead = head(i)
      end do
c
      if (.not. first) then
          write(outfil,9130) 
      endif
      if ( intl_version .eq. 1 ) then
         write(outfil,9121) 
       else
         if ( .not. first ) write(outfil,9121) 
      end if 
      first = .false.
      opoly = 0
      go to 2000
c
c                 Record Type 9: End of Buffer
c
 1090 continue
      go to 25
c
c                 Record Type 10: End of Data
c                     Put out end of frame information. Close files.
c
 1100 continue
      write(outfil,9130)
      close(infil)
      close(outfil)
      write(iout,*) '>> Execution Complete'
      stop
c
c                 Record Type 11: Scale Specification (Not used)
c
 1110 continue
      go to 2000
c
c                 Record Type 12: Polygon Fill
c                     Create Interleaf polygon with appropriate
c                     fill color and invisible edges.
c
 1120 continue
      call hex_to_int (line, jptr, 2, npts, str_len)
      call varlen( line, jptr, npts*2, 4, inrec, str_len, infil, iout)
      write(outfil,9040) counter, icolor
      write(outfil,9050) counter , counter 
      xy(1) = real(inrec(1))
      xy(2) = real(inrec(npts+1))
      call map(xy, xy1)
      do i = 2, npts
         xy(1) = real(inrec(i))
         xy(2) = real(inrec(npts+i))
         call map(xy, xy2) 
         write(outfil,9060) counter, xy1, xy2, 0, 127
         counter = counter + 1
         xy1(1) = xy2(1)
         xy1(2) = xy2(2)
      enddo
      xy(1) = real(inrec(1))
      xy(2) = real(inrec(npts+1))
      call map(xy, xy2)
      write(outfil,9060) counter, xy1, xy2, 0, 127
      counter = counter + 1
      write(outfil,9070)
      go to 2000
c
c                 Record Type 13: Fill Box
c                     Create Interleaf polygon with appropriate
c                     fill color and invisible edges.
c
 1130 continue
c
c                 Skip this if the second record after this is 
c                 a text record.  This avoids the black box used
c                 to erase the picture before the labels are put
c                 on the contours.
c
      read(infil, end=1100, fmt=9010) linetmp
      call stripf( linetmp, str_len_tmp, error)
      if (linetmp(1:2) .eq. 'FF') then
         line(1:) = linetmp(1:)
         jptr = 1
         str_len = str_len_tmp
         call hex_to_int( line, jptr, 2, irec, str_len)
         go to 2010
      endif
c
      call hex_to_int( line, jptr, 4, ix, str_len)
      call hex_to_int( line, jptr, 4, iy, str_len)
      xy(1) = real(ix)
      xy(2) = real(iy)
      call map(xy, xy1)
      call hex_to_int( line, jptr, 4, ix, str_len)
      call hex_to_int( line, jptr, 4, iy, str_len)
      xy(1) = real(ix)
      xy(2) = real(iy)
      call map(xy, xy2)
      write(outfil,9040) counter, icolor
      write(outfil,9050) counter, counter
      write(outfil,9060) counter,xy1(1), xy1(2), xy2(1), xy1(2), 0, 127
      counter = counter + 1
      write(outfil,9060) counter,xy2(1), xy1(2), xy2(1), xy2(2), 0, 127
      counter = counter + 1
      write(outfil,9060) counter,xy2(1), xy2(2), xy1(1), xy2(2), 0, 127
      counter = counter + 1
      write(outfil,9060) counter,xy1(1), xy2(2), xy1(1), xy1(2), 0, 127
      counter = counter + 1
      write(outfil,9070)
c
c         update new line
c
      line(1:) = linetmp(1:)
      jptr = 1
      str_len = str_len_tmp
      call hex_to_int( line, jptr, 2, irec, str_len)
      go to 2010
c
c                 Record Type 14: Scan Data Flag
c
 1140 go to 2000
c
c                 Record Type 15: Scan Line Data
c                     Create Interleaf polygons with appropriate
c                     fill color and invisible edges.
c
 1150 continue
      call scanline(npoly, ncolor, pixinc, ya, yh, line, 
     &     jptr, irec, infil, iout)
      if (npoly .gt. 200) then
        write(*,*) '*** Error: Insufficient number of polygons ***'
        write(*,*) '    Recompile with at least', npoly, ' polygons'
        write(*,*) '    Deleting output file'
        close(infil)
        close(outfil,status='delete')
        stop
      endif
c
c                 Check if number of colors is changed from 
c                 previous scan line.  If so, write out polygons.
c
      dopoly = npoly .ne. opoly
      do i = 1, min(opoly, npoly)
          dopoly = dopoly .or. (ncolor(i) .ne. ocolor(i))
      end do
      if ( dopoly ) then
          do i = 1, opoly
c
c                 Close current polygons.
c
              j = ghead
              ghead = space(j)
              space(j) = head(i)
              head(i) = j
              x(j) = xl(i)
              y(j) = yh
              j = ghead
              ghead = space(j)
              space(tail(i)) = j
              space(j) = 0
              tail(i) = j
              x(j) = xl(i+1)
              y(j) = yh
c
c                 Write out polygon except for blanks (color = 17)
c
              if (ocolor(i) .ne. 17) then
                  write(outfil,9040) counter, ocolor(i)
                  write(outfil,9050) counter, counter
                  n1 = head(i)
 1154             n2 = space(n1)
                  write(outfil,9060) counter, x(n1), y(n1), x(n2), 
     &                 y(n2), 0, 127
                  counter = counter + 1
                  if (n2 .ne. tail(i)) then
                      n1 = n2
                      go to 1154
                  endif
                  n1 = tail(i)
                  n2 = head(i)
                  write(outfil,9060) counter, x(n1), y(n1), x(n2), 
     &                 y(n2), 0, 127
                  counter = counter + 1
                  write(outfil,9070)
              endif
c
c                 Free up polygon list
c
              space(tail(i)) = ghead
              ghead = head(i)
          end do
      endif
c
c                 Get x positions on new scan line
c
      call getxpos(xl, npoly)
c
c                 Start new polygons
c
      if (dopoly) then
          opoly = npoly
          do i = 1, npoly
              ocolor(i) = ncolor(i)
              j = ghead
              ghead = space(j)
              head(i) = j
              x(j) = xl(i)
              y(j) = yh
              j = ghead
              ghead = space(j)
              tail(i) = j
              space(head(i)) = j
              space(j) = 0
              x(j) = xl(i+1)
              y(j) = yh
          end do
      endif
c
c                 Add to existing polygons
c
      do i = 1, npoly
          j = ghead
          ghead = space(j)
          space(j) = head(i)
          head(i) = j
          x(j) = xl(i)
          y(j) = ya
          j = ghead
          ghead = space(j)
          space(tail(i)) = j
          space(j) = 0
          tail(i) = j
          x(j) = xl(i+1)
          y(j) = ya
      end do
      go to 2010
c
c                 Record Type 16: Look-Up Table
c
c                     Black & White LUT were written by the
c                     document header format.        
c
c                     For color doc, we read the Patran file LUT of
c                     RGB here, convert to percentages of CMYK
c                     & write to doc. We also force color zero
c                     to be black so that element lines on color fringe
c                     plots are black.
c
c                     also issue frame advance because Patran puts
c                     out a frame advance before the LUT is given.
c                     the first frame advance instruction for
c                     color is ignored.
c
 1160 continue
      call hex_to_int( line, jptr, 4, nstart, str_len )
      call hex_to_int( line, jptr, 4, ncol, str_len )
c      if (ncol .gt. 16) then
c          close(outfil,status='delete')
c          write(iout,9010) ' >> Number of colors greater than 16.'
c          write(iout,9010) '    Abandoning processing.'
c          stop
c      endif
      call varlen( line, jptr, ncol*3, 2, inrec, str_len, infil, iout)
c
      if ( intl_version .eq. 2 ) then
          write(outfil,9113)
          write(outfil,9114) 0, 0.0, 0.0, 0.0, 100.0
          do i = 1, ncol
             if (nstart+i-1 .eq. 0) cycle
             red      = 100.0 * real(inrec(3*i - 2)) / 255.0
             green    = 100.0 * real(inrec(3*i - 1)) / 255.0
             blue     = 100.0 * real(inrec(3*i    )) / 255.0 
             cyan     = max(red,green,blue) - red
             magenta  = max(red,green,blue) - green
             yellow   = max(red,green,blue) - blue
             black    = 100.0 - max(red,green,blue)
             write(outfil,9114) nstart+i-1, cyan, magenta, yellow,
     &            black
          end do        
c
          write(outfil,9115)
          write(outfil,9121)
      end if                                          
      go to 2000
c
c                 Record Type 17: Continuation (ignore this)
c
 1170 continue
      go to 2000
c
c                 Record Type 18: Version 
c                     Print out version number.
c
 1180 continue
      write(iout,9150) line(ptr:str_len)
      go to 2000
c
c                 Record Type 19: Unknown (ignore this)
c
 1190 continue
      go to 2000
c
c                 Record Type 20: Text Style (ignore this)
c
 1200 continue
      go to 2000
c
c                 Record Type 21: Text Scale (ignore this)
c
 1210 continue
      go to 2000
c
c                 Record Type 22: Text Rotation
c                      Find and store the text rotation
c
 1220 continue
      call hex_to_int (line, jptr, 4, idum, str_len)
      call hex_to_real (line, jptr, rotation, str_len)
      rotation = 360.0 - rotation
      go to 2000
c
c
c                 Record Type 255: Text
c
 1900 continue
c
c                    get start position for texxt
c      
      call hex_to_int( line, jptr, 4, ix, str_len)
      call hex_to_int( line, jptr, 4, iy, str_len)
      xy(1) = real(ix)
      xy(2) = real(iy)
      call map(xy, xy1)
c
c                 Quote will escape certain characters in the text
c                 string so Interleaf can understand them. Output 
c                 text in black.  Add 0.1 to y coordinate since
c                 Patran gives upper left coordinate of string.
c                 and Interleaf wants middle left coordinate.
c
      call hex_to_int( line, jptr, 2, num_char, str_len)
      string(1:) = line(jptr:jptr + num_char)
      call quote(string, string2, num_char, len2, maxlen2, iout)
      xy1(2) = xy1(2) + 0.1
      write(outfil,9090) counter, xy1, rotation, string2(1:len2)
      counter = counter + 1
      go to 2000
c
c
c                 Move pointer to next record
c
 2000 continue
      go to 25
 2010 continue
      go to 50
c
c                 end of interpreter loop for hardcopy file
c
c                 Formats
c
 9000 format(a,$)
 9010 format(a)
 9020 format('Patran_',a8)
 9022 format('>> $HOME environment variable not found.',/,
     &       '    Job terminated.' )
c 9030 format('(v6,',i6,',0,',4(f8.5,','),i3,',',i3,',1,0)')
 9030 format(' (v7,',i6,',0,,',4(f8.5,','),i3,',',i3,',0.125,0)')
 9040 format('(p8,',i6,',0,5,',i3,',0')
 9050 format(' (g9,',i6,',0',/,'  (g9,',i6,',0')
c 9060 format('   (v6,',i6,',0,',4(f8.5,','),i3,',',i3,',1,0)')
 9060 format('   (v7,',i6,',0,,',4(f8.5,','),i3,',',i3,',0.125,0)')
 9070 format('  )',/, ' )', /, ')')
 9080 format(20a4)
 9090 format('(t14,',i6,',0,',2(f8.5,','),'0,0,0,',f11.7,
     &       ',,wst:helvps8,',a,')')
 9100 format(' ***************************************************',/,
     &       ' *                                                 *',/,
     &       ' *     Patran to InterLeaf Hardcopy Translator     *',/,
     &       ' *                                                 *',/,
     &       ' *           For Interleaf Versions 5 or 6         *',/,
     &       ' *              For Release 6 of Patran            *',/,
     &       ' *               (last update: 2-3-97)             *',/,
     &       ' *                                                 *',/,
     &       ' ***************************************************',///)
c
 9111 format( '<!OPS, Version= 8.0>',/,
     & '<!Font Definitions,',/,
     & ' F2= Times 10,F3= Times 12 Bold,F4= Times 12,',/,
     & ' F5= Helvetica 14 Bold,F6= Symbols 12>',/,
     & '<!Page,',/,
     & ' Height=8.50 Inches,Width=11 Inches,Bottom Margin=1 Inches,',/,
     & ' Left Margin=1 Inches,Right Margin=1 Inches,Hyphenation=on,',/,
     & ' Consecutive Hyphens= 3,Revision Bar Placement= Left,',/,
     & ' Feathering= off>',/,
     & '<!Color Definitions,',/,
     & '   C7 =  0.00,  0.00,  0.00,',/,
     & '   C6 =  6.67,  6.67,  6.67,',/,
     & '  C13 = 13.33, 13.33, 13.33,',/,
     & '   C4 = 20.00, 20.00, 20.00,',/,
     & '  C12 = 26.67, 26.67, 26.67,',/,
     & '   C2 = 33.33, 33.33, 33.33,',/,
     & '  C11 = 40.00, 40.00, 40.00,',/,
     & '  C10 = 46.67, 46.67, 46.67,',/,
     & '  C15 = 53.33, 53.33, 53.33,',/,
     & '  C14 = 60.00, 60.00, 60.00,',/,
     & '   C5 = 66.67, 66.67, 66.67,',/,
     & '   C3 = 73.33, 73.33, 73.33,',/,
     & '   C9 = 80.00, 80.00, 80.00,',/,
     & '   C8 = 86.67, 86.67, 86.67,',/,
     & '   C1 = 93.33, 93.33, 93.33,',/,
     & '   C0 =100.00,100.00,100.00>',/,
     & '<!Class, "para",',/,
     & '  Top Margin=0.04 Inches,Bottom Margin=0.04 Inches,',/,
     & '  Font=F4@Z7@Lam,Line Spacing=1.1606 lines,',/,
     & '  Left Tab= 0/1*3 Inches,Composition=Optimum>',/,
     & '<!Master Frame,',/,
     & '  Name="page",Placement=At Anchor,Vertical Alignment=Top,',/,
     & '  Width=2 Inches,Width=Page Without Margins,',/,
     & '  Height=1 Inches,Height=Page Without Margins,',/,
     & '  Diagram=',/,
     & 'V11,',/,
     & '(g9,0,0,',/,
     & ' (E16,0,0,,5,1,1,0.0533333,1,15,0,0,1,0,0,0,',/,
     & '  1,5,127,7,0,0,7,0,1,1,0.0666667,0.06',/,
     & '  66667,6,6,0,0.0666667,6))>',/,
     & '<!End Declarations>' )
c
 9112 format( '<!OPS, Version = 8.0>',
     & '<!Document,',/,
     & '        Header Page = no,',/,
     & '        Final Output Device =   "ileaf",',/,
     & '        Default Printer =       "p3129c",',/,
     & '        Default Page Stream Name = "page">',/,
     & '<!Font Definitions,',/,
     & '        F2 = Times 10,',/,
     & '        F4 = Thames 10>',/,
     & '<!Page,',/,
     & '        Height =                8.50 Inches,',/,
     & '        Width =                 11 Inches,',/,
     & '        Bottom Margin =         1 Inches,',/,
     & '        Left Margin =           0.75 Inches,',/,
     & '        Right Margin =          1.25 Inches,',/,
     & '        Hyphenation =           on,',/,
     & '        Consecutive Hyphens =   3,',/,
     & '        Revision Bar Placement = Left,',/,
     & '        Feathering =            off>',/,
     & '<!Class, "para",',/,
     & '        Top Margin =            0.07 Inches,',/,
     & '        Bottom Margin =         0.07 Inches,',/,
     & '        Font =                  F4@Z7@Lam,',/,
     & '        Line Spacing =          1.3132 lines>',/,
     & '<!Master Frame,',/,
     & '        Name =                  "page",',/,
     & '        Placement =             At Anchor,',/,
     & '        Vertical Alignment =    Top,',/,
     & '        Width =                 0.41 Inches,',/,
     & '        Width =                 Page Without Margins,',/,
     & '        Height =                0.1366667 Inches,',/,
     & '        Height =                Page Without Margins,',/,
     & '        Diagram =',/,
     & 'V11,',/,
     & '(g9,0,0,',/,
     & ' (E16,0,0,,5,1,1,0.0533333,1,15,0,0,1,0,0,0,1,',/,
     & '  5,127,7,0,0,7,0,1,1,0.0666667,0.06',/,
     & '  66667,6,6,0,0.0666667,6))>')
c
 9113 format( '<!Color Definitions,' )
 9114 format( '  C',i3.3,' =',3(f6.2,','),f6.2,',')
 9115 format( '>',/,'<!End Declarations>' )
 9116 format(
     & '  C0 =    0.0,   0.0,   0.0, 100.0,',/,
     & '  C16 =   0.0, 100.0, 100.0, 0.0,',/,
     & '  C17 = 100.0,   0.0, 100.0, 0.0,',/, 
     & '  C18 = 100.0, 100.0,   0.0, 0.0,',/,
     & '  C19 = 100.0,   0.0,   0.0, 0.0,',/,
     & '  C20 =   0.0, 100.0,   0.0, 0.0,',/,
     & '  C21 =   0.0,   0.0, 100.0, 0.0,',/,
     & '  C22 =   0.0,   0.0,   0.0, 0.0>'  )
c      
 9121 format( '<para>', /,
     & '<|,1>', /,
     & '<Frame,', /,
     & '      Name =page, Placement =At Anchor,', /,
     & '      Vertical Alignment =Top,', /,
     & '      Width =Page Without Margins,',/,
     & '      Height =Page Without Margins,',/,
     & '      Diagram =',/,
     & '  V11,',/,
     & '(g9,1,0') 
c
 9130 format( ') >')
 9140 format(' >> Interleaf Document: ', a)
 9150 format(' >> Patran Version: ', a)
      end
c
c *****************************************************************
c *                                                               *
c *      s u b r o u t i n e  -- v a r l e n                      *
c *                                                               *
c *****************************************************************
c
      subroutine varlen( line, ptr, num_vals, size, inrec, str_len,
     &     infil, iout) 
      implicit integer (a-z)
      integer inrec(*)
      character *(*) line
      logical error
c
c     
      counter = 1
      do while (.true.)
         call hex_to_int( line, ptr, size, inrec(counter), str_len)
         counter = counter + 1
         if (counter .gt. num_vals) exit
         if (ptr .gt. str_len) then
 10         continue
            read(infil, end=1000, fmt=9020) line
            call stripf( line, str_len, error )
            if ( error ) then
               write (iout,*) ' >> ERROR: line is blank. Skipping to'
               write (iout,*) ' >> next line.'
               go to 10
            endif
            ptr = 1
         endif
      enddo
c
      return
c
c         reached end of file
c
 1000 continue
      write (iout, *) ' >> fatal error: end of file reached before'
      write (iout, *) ' >>    end of record found.  Aborting.'
      stop
c
 9020 format(a)
      end      
c
c *****************************************************************
c *                                                               *
c *      s u b r o u t i n e  -- m a p                            *
c *                                                               *
c *****************************************************************
c
      subroutine map(xy, xym)
c
c                 Map coordinates in hardcopy file system (xy)
c                 to coordinates of output system (xym)
c 
      real xy(2), xym(2)
      common /device/ xr, yr, ratior, pixinc, devx, devy
      xym(1) = xy(1) * devx / xr
      xym(2) = devy - xy(2) * devy / yr
      return
      end
c *****************************************************************
c *                                                               *
c *      s u b r o u t i n e  -- q u o t e                        *
c *                                                               *
c *****************************************************************
c
      subroutine quote(strng1, strng2, len, len2, maxlen2, iout )
c
c                 Quote special characters ( < > ( ) , \  sp ) 
c                 with a backslash (\) so Interleaf can interpret
c                 them correctly 
c 
      character*(*) strng1, strng2
      len2 = 0
      do i = 1, len
          if ((strng1(i:i) .eq. '<') .or.
     &        (strng1(i:i) .eq. '>') .or.
     &        (strng1(i:i) .eq. '(') .or.
     &        (strng1(i:i) .eq. ')') .or.
     &        (strng1(i:i) .eq. ',') .or.
     &        (strng1(i:i) .eq. '\') .or.
     &        (strng1(i:i) .eq. ' ')) then
              len2 = len2 + 1                            
              if ( len2 .ge. maxlen2 ) then
                write(iout,9000) maxlen2
                stop
              end if
              strng2(len2:len2) = '\'
          endif
          len2 = len2 + 1
          if ( len2 .ge. maxlen2 ) then
            write(iout,9000) maxlen2
            stop
          end if
          strng2(len2:len2) = strng1(i:i)
      end do
      return          
c
 9000 format(/,'>> Fatal error: the string length generated',
     &       /,'                for Interleaf file exceeds the', 
     &       /,'                limit of: ',i5,' characters.',
     &       /,'                Shorten Patran results titles',
     &       /,'                and re-try! Job terminated.' )
c
      end
c
c *****************************************************************
c *                                                               *
c *      s u b r o u t i n e  -- s t r i p f                      *
c *                                                               *
c *****************************************************************
c
      subroutine stripf( string, str_len, error )  
      implicit integer (a-z) 
      character *(*) string, copy*256
      logical error
      save line
      data line /1/
c
c             strip leading and ending blanks from string and return it, 
c             along with length of string
c
      line = line + 1
      error = .false.
      last = len(string)
      copy(1:last) = string(1:last)
      start = 1
      do while (string(start:start) .eq. ' ')
         start = start + 1
         if (start .gt. last) go to 1000
      enddo
      do while (string(last:last) .eq. ' ')
         last = last - 1
         if (last .lt. start) go to 1000
      enddo
c
      string(1:) = copy(start:last)
      str_len = last - start + 1
c
      return
c
 1000 continue
      error = .true.
      return
      end
c
c *****************************************************************
c *                                                               *
c *      s u b r o u t i n e  -- s c a n l i n e                  *
c *                                                               *
c *****************************************************************
c
      subroutine scanline(npoly, ncolor, pixinc, ya, yh, line, 
     &     jptr, irec, infil, iout)
c 
c   Read in scan lines. Combine scan lines with same y level
c   into a single scan line. Return list of colors and
c   x positions for polygons and number of polygons. 
c 
      character*256 line
      common /xlpcom/ xlp(201)
      dimension ncolor(200), inrec(50000), xy(2), xy1(2), xy2(2)
      logical first, error
      integer str_len
c
c            Initialize
c
      npoly = 0
      first = .true.
c
  100 continue 
c
c            Map x and y positions from scan line data, and get number of
c            pixels
c
      call hex_to_int( line, jptr, 4, ix, str_len )
      call hex_to_int( line, jptr, 4, iy, str_len )
      xy(2) = real(iy) 
      xy(1) = real(ix) 
      call map(xy, xy1)
      call hex_to_int( line, jptr, 4, ix, str_len )
      xy(1) = real(ix) 
      xy(2) = xy(2) + pixinc / 2
      call map(xy, xy2)
      call hex_to_int( line, jptr, 4, npix, str_len )
c
c            Save y position of first scan line.
c
      if (first) ypos = xy1(2)
c
c            Check y position of new scan line against first scan line.
c            Add blank polygon between scan lines on the same y level.
c
      if (ypos .ne. xy1(2)) then
        return
      else if (first) then
        first = .false.
      else
        npoly = npoly + 1
        ncolor(npoly) = 17
      endif
c
c            Find number of colors in scan line and x positions of
c            color changes.
c
      ya = xy1(2)
      yh = xy2(2)
      xs = xy1(1)
      xf = xy2(1)
      call varlen( line, jptr, npix, 4, inrec, str_len, infil, iout)
      npoly = npoly + 1
      ncolor(npoly) = inrec(1)
      xlp(npoly) = xs
      do i = 1, npix
        if (inrec(i) .ne. ncolor(npoly)) then
          npoly = npoly + 1
          ncolor(npoly) = inrec(i)
          xlp(npoly) = xs + (xf-xs) * float(i-1) / npix
        endif
      end do
      xlp(npoly+1) = xf
c
c            Go get next line 
c
 200  continue
      jptr = 1
      read(infil, end=1000, fmt=9020) line
      call stripf( line, str_len, error)
      if ( error ) then
         write (iout,*) ' >> ERROR: line is blank. Skip to next line.'
         go to 200
      endif
      call hex_to_int( line, jptr, 2, irec, str_len)
c
c            Return if new line is not a scan line record.
c
      if (irec .ne. 15) return
      go to 100
c
c         reached end of file
c
 1000 continue
      write (iout, *) ' >> fatal error: end of file reached before'
      write (iout, *) ' >>    end of record found.  Aborting.'
      stop
c
 9020 format(a)
      end
c
c *****************************************************************
c *                                                               *
c *      s u b r o u t i n e  -- g e t x p o s                    *
c *                                                               *
c *****************************************************************
c
      subroutine getxpos(xl, npoly)
c
c             Return x positions of polygons
c 
      common /xlpcom/ xlp(201)
      dimension xl(201)
      do i = 1, npoly + 1
        xl(i) = xlp(i)
      end do
      return
      end
c
c *****************************************************************
c *                                                               *
c *      s u b r o u t i n e  -- h e x _ t o _ i n t              *
c *	                                                          *
c *          takes hex numbers as character data and converts     *
c *	     then into an integer.                                *
c *                                                               *
c *****************************************************************
c
      subroutine hex_to_int( string, ptr, num_hexchars, inum, len )  
      implicit integer (a-z) 
      character *(*) string, char*1
      data intzero, intnine, inta, intf / 48, 57, 65, 70 /
c	
c       string       = current line of characters to be treated as hex
c       ptr          = pointer to current place in string to process
c       num_hexchars = number of chars of hex to read
c       inum         = returned integer 
c       len          = length of initial string
c
      inum = 0
      do i = ptr, num_hexchars+ptr-1
c
         if (i .gt. len) then
            write (*,9010)
            stop
         endif
c
         char = string(i:i)
         char_num = ichar(char)
         if (char_num .ge. intzero .and. char_num .le. intnine) then
            inum = inum * 16 + char_num - intzero
         else if (char_num .ge. inta .and. char_num .le. intf) then
            inum = inum * 16 + char_num - inta + 10
         else
            write (*,9000) char
            stop
         endif
c
      enddo
      ptr = ptr + num_hexchars
c
c
      return
 9000 format(/,' >>> ERROR: character ',a1,' not recognized for',
     &       /,'            integer. Aborting Program.',/)
 9010 format(/,' >>> ERROR: attempt to access character beyond end',
     &       /,'            of line. Aborting Program.',/)
      end
c
c *****************************************************************
c *                                                               *
c *      s u b r o u t i n e  -- h e x _ t o _ r e a l            *
c *	                                                          *
c *          takes hex numbers as character data and converts     *
c *	     then into a real number.                             *
c *                                                               *
c *****************************************************************
c
      subroutine hex_to_real( string, ptr, rnum, len )  
      implicit integer (a-z) 
      character *(*) string
      real rnum
c
      call hex_to_int( string, ptr,   4, numerator, len)
      call hex_to_int( string, ptr, 4, denominator, len)
c
      rnum = real (numerator) / real (denominator)
c
c
      return
      end
