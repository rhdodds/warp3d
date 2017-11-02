c     ****************************************************************
c     *                                                              *
c     *                      subroutine ininc                        *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 10/27/2017 rhd             *
c     *                                                              *
c     *     read element incidences                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine ininc( sbflg1, sbflg2 )
c
      use global_data ! old common.main
      use main_data, only : incmap, incid
      implicit none
c
      real :: dumr
      double precision :: dumd
      character :: dums
      logical :: sbflg1, sbflg2
      logical, external :: integr, label, matchs, string, endcrd
      integer :: intlst(mxlsz), elem, param, errnum, nnode, i, count,
     &           icn, iplist, incpos, stnod, lenlst, dum
c
c                       if the subroutine has been previously
c                       entered and exited, then there was an
c                       error in the data of the last processed
c                       card. specifically, that error was the
c                       absence of the element number for incidence
c                       input. under this circumstance, print an
c                       error message and continue with input.
c
      if(sbflg1) call errmsg(39,dum,dums,dumr,dumd)
      if( matchs('from',4) ) call splunj
      if( matchs('file',4) ) call ininc_file
c
 605  call readsc
c
c
c **********************************************************************
c *                                                                    *
c *                     read in element. if no element is given, then  *
c *                     either the element has been forgotten or       *
c *                     incidence input has ended. in either           *
c *                     case, set sbflg1 to true and exit the          *
c *                     subroutine. if the element has been forgotten, *
c *                     then the rest of the line will be skipped.     *
c *                                                                    *
c **********************************************************************
c
c
      if(.not.integr(elem)) then
         go to 9999
      end if
c
c                       check that the element input does not exceed
c                       the number of elements in the structure.
c
      if(elem.gt.noelem) then
         param= elem
         call errmsg(35,param,dums,dumr,dumd)
         go to 605
      end if
c
c                       check that the element input is not negative.
c
      if(elem.lt.0) then
         param= elem
         call errmsg(86,param,dums,dumr,dumd)
         go to 605
      end if
c
c
      call scan
      call trlist(intlst,mxlsz,mxlsz,lenlst,errnum)
c
c                       branch on the return code from trlist. a
c                       value of 1 indicates no error. a value of
c                       2 indicates that the parse rules failed in
c                       the list. a value of 3 indicates that the
c                       list overflowed its maximum length of mxlsz.
c                       a value of 4 indicates that no list was found.
c                       in these last 3 cases, the rest of the card
c                       will be ignored and a new card will be sought.
c
      if(errnum.eq.2) then
         param= 1
         call errmsg(24,param,dums,dumr,dumd)
         go to 605
      else if(errnum.eq.3) then
         param= 2
         call errmsg(24,param,dums,dumr,dumd)
         go to 605
      else if(errnum.eq.4) then
         param=4
         call errmsg(24,param,dums,dumr,dumd)
         go to 605
      else
         if(errnum.eq.1) then
            call backsp(1)
            go to 610
         end if
         param= 3
         call errmsg(24,param,dums,dumr,dumd)
         go to 605
      end if
c
c
c **********************************************************************
c *                                                                    *
c *                     there is a valid list of incidences. parse     *
c *                     the list.                                      *
c *                                                                    *
c **********************************************************************
c
c
 610  nnode = iprops(2,elem)
c
c                       find the starting point for incidences in
c                       the incidences vector incid. if overwrite,
c                       take old st. point, if new, set it.
c
      if( incmap(elem) .eq. 0 ) then
         incpos = inctop
      else
         incpos = incmap(elem)-1
      end if
c
c                       initialize the incidence vector for error
c                       checking purposes.
c
      do i = 1, nnode
         incid(incpos+i) = 0
      end do
c
      count  = 0
      icn    = 0
      iplist = 1
 615  call trxlst(intlst,lenlst,iplist,icn,stnod)
c
      count= count+1
c
      if( stnod .gt. nonode ) then
         param = stnod
         call errmsg(16,param,dums,dumr,dumd)
         count = count-1
         go to 620
      end if
c
c                       check that the list node is not negative.
c
      if( stnod .lt. 0 ) then
         param = stnod
         call errmsg(58,param,dums,dumr,dumd)
         count = count-1
         go to 620
      end if
c
c                       check that the number of incidences input so
c                       far does not exceed the number of nodes for
c                       this element.
c
      if( count .gt. nnode ) then
         param = count*two16+nnode
         call errmsg(37,param,dums,dumr,dumd)
         iplist = 0
         go to 620
      end if
c
      incpos        = incpos+1
      incid(incpos) = stnod
c
 620  if( iplist .eq. 0 ) then
c
c                       check that as many incidences have been input
c                       as there are nodes for this element.
c
         if( count .lt. nnode ) then
            param = count*two16+nnode
            call errmsg(38,param,dums,dumr,dumd)
         else
c
c                       have the incidence mapping vector point to
c                       the current element's incidences.
c
            if( incmap(elem) .eq. 0 ) then
               incmap(elem) = inctop+1
               inctop       = incpos
            end if
         end if
      else
         go to 615
      end if
c
c                       the incidence list has been processed. return
c                       to input another element's incidences.
c
      go to 605
c
c
 9999 sbflg1 = .true.
      sbflg2 = .true.
c
c
      return
c
      contains
c     ========
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ininc_file                   *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 10/28/2017 rhd             *
c     *                                                              *
c     *     process reading element incidences from a stream or      *
c     *     binary unformatted file                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine ininc_file
      implicit none
c
      logical :: stream, nameok
      integer :: nchar, ierror, ifile, status, numentries
      integer, external :: warp3d_get_device_number
      integer, allocatable :: buffer(:)
      character(len=200) :: filename, tfilename
c
c              enter with word "file" having just been recognized by
c              scanner. allow label or string for file name. if string,
c              process ~ in name as required to generate absolute path
c
      filename = " "
      tfilename = " "
c
      if( endcrd(dum) ) return
c
      if( label(dum) ) then
          call entits ( filename, nchar )
          stream = .false.
          if( matchs('stream',5) ) stream = .true.
      elseif( string(dum) ) then
         call entits ( tfilename, nchar )
         call tilde( tfilename, filename, nameok )
         if( .not. nameok ) then
            call errmsg(189,dum,tfilename,dumr,dumd)
            return
         end if
         stream = .false.
         if( matchs('stream',5) ) stream = .true.
      else
         call errmsg( 175, dum, dums, dumr, dumd )
         return
      end if
c
c              confirm file exists.
c
      inquire( file = filename, iostat = ierror, exist = nameok )
      if( ierror .gt. 0 ) then
         call errmsg(176,dum,filename,dumr,dumd)
         return
      elseif( .not. nameok ) then
         call errmsg(177,dum,filename,dumr,dumd)
         return
      end if
c
c              open file. note use of "segmented" record type for
c              binary to support extremely long records.
c
      ifile = warp3d_get_device_number()
      if( stream ) then
         open( unit=ifile, file=filename, status='old', access="stream",
     &         form="unformatted" )
      else
         open( unit=ifile, file=filename, status='old',
     &         access='sequential',form='unformatted',
     &         recordtype='segmented' )
      end if
c
      write(out,9000)
c
c              allocate buffers and read file. note ordering,
c              all x, then all y, then all z.
c
      read(ifile,iostat=status) numentries
      if( status .ne. 0 ) then
         num_error = num_error + 1
         if( status < 0 ) write(out,9020)
         if( status > 0 ) write(out,9030)
         return
      end if
      if( numentries > noelem*mxndel ) then
           num_error = num_error + 1
           return
      end if
c
      allocate( buffer(numentries) )
c
      read(ifile,iostat=status) buffer
      if( status .ne. 0 ) then
         num_error = num_error + 1
         if( status < 0 ) write(out,9020)
         if( status > 0 ) write(out,9030)
         return
      end if
c
      close( unit=ifile )
      write(out,9010)
c
c              create the incmap for element values in order. save
c              structure node for each element inregular data structure
c
      inctop = 1
      do elem = 1, noelem
       nnode = iprops(2,elem)
       incmap(elem) = inctop
       do i = 1, nnode
         incid(inctop+i-1) = buffer(inctop+i-1)
       end do
       inctop = inctop + nnode
      end do
c
      deallocate( buffer )
c
      return
c
 9000 format(5x,'... file opened')
 9010 format(5x,'... incidences read. file closed')
 9020 format(/1x,'>>>>> error: end-of-file reading incidences' )
 9030 format(/1x,'>>>>> error: unknown error reading incidences' )
c
      end subroutine ininc_file
      end subroutine ininc




