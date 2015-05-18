c     ****************************************************************
c     *                                                              *
c     *                      subroutine ininc                        *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 07/27/89                   *
c     *                                                              *
c     *     this subroutine supervises and conducts the input of     *
c     *     element incidences.                                      * 
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine ininc( sbflg1, sbflg2 )
      use main_data, only : incmap, incid
      implicit integer (a-z)
$add common.main
      real dumr
#dbl      double precision
#sgl      real
     &   dumd
      character dums
      logical sbflg1,sbflg2
      logical integr
      dimension intlst(mxlsz)
c
c                       if the subroutine has been previously
c                       entered and exited, then there was an
c                       error in the data of the last processed
c                       card. specifically, that error was the
c                       absence of the element number for incidence
c                       input. under this circumstance, print an
c                       error message and continue with input.
c                       
      if(sbflg1) then
         call errmsg(39,dum,dums,dumr,dumd)
      end if
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
      end




