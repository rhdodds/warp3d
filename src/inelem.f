c     ****************************************************************
c     *                                                              *
c     *                      subroutine inelem                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 08/20/2017 rhd             *
c     *                                                              *
c     *     this subroutine supervises and conducts the input of     *
c     *     element type and properties.                             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine inelem( sbflg1, sbflg2 )
c
      use global_data ! old common.main
      use main_data, only: matprp, imatprp
c
      implicit none
      logical :: sbflg1, sbflg2
c
      integer :: intlst(mxlsz), minum, lenlst, param, dum, nc, i, type,
     &           outloc, strcon, intord, outfmt, geonl, bbar, surf,
     &           matn, matnum, errnum
      logical :: found, defmat, deftyp, macrointer
      logical, external :: matchs, endcrd, true, label, scanms, numr
      double precision :: dumd
      real :: dumr, area
      character :: name*80, dums*1
      character(len=8) :: tname
      character(len=24) :: mname
c
      macrointer = .false.
      minum = -1
      tname(1:8) = 'nodefalt'
      mname(1:8) = 'nodefalt'
c
c                       if sbflg1 is set, then there is re-entry
c                       into this subroutine. in this case, assign
c                       the default value of all to the element
c                       list.
c
      if(sbflg1) then
         intlst(1)= 1
         intlst(2)= -noelem
         intlst(3)= 1
         lenlst= 3
         go to 511
      end if
c
c
 505  call readsc
c
c                       translate the list of elements input on
c                       this line.
c
      call scan
      call trlist(intlst,mxlsz,noelem,lenlst,errnum)
c
c                       branch on the return code from trlist. a
c                       value of 1 indicates no error. a value of
c                       2 indicates that the parse rules failed in
c                       the list. a value of 3 indicates that the
c                       list overflowed its maximum length of mxlsz.
c                       in these last two cases, the rest of the card
c                       will be ignored and a new card will be sought.
c                       a value of 4 indicates that no list was found.
c                       in this case, either elements input has ceased
c                       or the default value of all is assigned to the
c                       list.
c
      if(errnum.eq.2) then
         param= 1
         call errmsg(24,param,dums,dumr,dumd)
         go to 505
      else if(errnum.eq.3) then
         param= 2
         call errmsg(24,param,dums,dumr,dumd)
         go to 505
      else if(errnum.eq.4) then
         call backsp(1)
         go to 9999
      else
         if(errnum.eq.1) then
            call backsp(1)
            go to 511
         end if
         param= 3
         call errmsg(24,param,dums,dumr,dumd)
         go to 505
      end if
c
c
c **********************************************************************
c *                                                                    *
c *                     an element list exists. read an element type   *
c *                     for these elements. an element type must be    *
c *                     given the first time through. on subsequent    *
c *                     lists, the previous lists' type will be used   *
c *                     as a default. if on the initial pass a list    *
c *                     is not given, or if the type given is not      *
c *                     found, ignore the current element list and     *
c *                     seek a new card.                               *
c *                                                                    *
c **********************************************************************
c
c
 511  if(matchs('type',4)) then
         if(label(dum)) then
            name= ' '
            tname= ' '
            call entits(name,nc)
            if(nc.gt.8) nc=8
            tname(1:nc)= name(1:nc)
            found= .false.
            do i = 1, nlibel
               if(scanms(tname,elelib(i),8)) then
                  found= .true.
                  type = i
                  go to 513
               end if
            end do
c
 513        if(.not.found) then
               call errmsg(26,dum,dums,dumr,dumd)
               tname= 'nodefalt'
               go to 505
            end if
            deftyp= .false.
         else
            call errmsg(27,dum,dums,dumr,dumd)
            tname= 'nodefalt'
            go to 505
         end if
      else
         if(scanms(tname,'nodefalt',8)) then
            call errmsg(32,dum,dums,dumr,dumd)
            go to 505
         else
            deftyp= .true.
         end if
      end if
c
c
c **********************************************************************
c *                                                                    *
c *                     a valid list and an element type for that list *
c *                     have been input. input properties.             *
c *                                                                    *
c **********************************************************************
c
c
c                       initialize the output location of stresses,
c                       strains, etc. and stress type to be output
c                       to default values.
c
      outloc= 4HDEFA
      strcon= 4HCURR
      intord= 4HDEFA
      outfmt= 4HDEFA
      geonl = 4HFLSE
      bbar  = 4hTRUE
      surf  = 4HDEFA
      defmat= .true.
      area = 0.0e0
c
      if(endcrd(dum)) call readsc
c
c                       properties on next card(s)
c
      if(matchs('properties',10)) then
         if(endcrd(dum)) go to 515
         go to 520
      else
         go to 520
      end if
c
c                       read the name of the property and its value.
c                       if there is no match between what is input
c                       and the property list, then print an error
c                       message and continue scanning for a property
c                       name. if an end of card is encountered, print
c                       another error message and return to read a
c                       high level command.
c
 515  call readsc
c
 520  if(matchs('order',5)) go to 530
      if(matchs('gausspts',8)) go to 540
      if(matchs('nodpts',6)) go to 550
      if(matchs('center_output',6)) go to 555
      if(matchs('current',4)) go to 560
      if(matchs('reference',3)) go to 565
      if(matchs('short',5)) go to 570
      if(matchs('long',4)) go to 575
      if(matchs('material',8)) go to 580
      if(matchs('linear',6)) go to 585
      if(matchs('nonlinear',6)) go to 586
      if(matchs('nlgeom',5)) go to 586
      if(matchs('geonl',5)) got o 586
      if(matchs('bbar',4)) go to 587
      if(matchs('no_bbar',6)) go to 588
      if(matchs('nobbar',6)) go to 588
      if(matchs('surface',4)) go to 595
      if(matchs('area',4)) go to 600
c
      if(matchs(',',1)) go to 525
c
c                       there is no match with existing properties.
c                       check for end of card. if not, print error
c                       message.
c
      if(endcrd(dum)) then
         go to 900
      else
         call errmsg(4,dum,dums,dumr,dumd)
         if(true(dum)) go to 520
      end if
c
c                       comma separator. if at the end of a line,
c                       denotes continuation. otherwise, ignore.
c
 525  continue
      if(endcrd(dum)) then
         go to 515
      else
         go to 520
      end if
c
c
c **********************************************************************
c *                                                                    *
c *                     input the values of the properties.            *
c *                                                                    *
c *                     input the order of integration for the         *
c *                     elements in the list. if the order is not      *
c *                     input, use the default value for the element.  *
c *                     as this default value is not known at the      *
c *                     moment, stuff an entity indicative of default  *
c *                     status in intord.                              *
c *                                                                    *
c **********************************************************************
c
c
 530  continue
      if(matchs('3x3x3',5)) then
         intord= 4HO333
      else if(matchs('2x2x2',5)) then
         intord= 4HO222
      else if(matchs('2x2ndl',5)) then
         intord= 4HO22N
      else if(matchs('2x2gs',5)) then
         intord= 4HO22G
      else if(matchs('14pt_rule',4)) then
         intord= 4HO14P
      else if(matchs('9pt_rule',3)) then
         intord= 4HO09P
      else if(matchs('1pt_rule',3)) then
         intord= 4HO01P
      else if(matchs('4pt_rule',3)) then
         intord= 4HO04P
      else if(matchs('3pt_rule',3)) then
         intord= 4HO03P
      else if(matchs('3mpt_rule',4)) then
         intord= 4HO3MP
      else if(matchs('5pt_rule',3)) then
         intord= 4HO05P
      else if(matchs('6pt_rule',3)) then
         intord= 4HO06P
      else if(matchs('7pt_rule',3)) then
         intord= 4HO07P
      else
         call errmsg(28,dum,'ordr',dumr,dumd)
      end if
      go to 520
c
c
c **********************************************************************
c *                                                                    *
c *                     input the output location for stresses,        *
c *                     strains, etc. if it is not input, assume       *
c *                     a default value for the elements in the list.  *
c *                     as this default value is not known at the      *
c *                     moment, an entity indicative of default        *
c *                     status has been previously stuffed in outloc.  *
c *                                                                    *
c **********************************************************************
c
c
c                       output location is gauss points
c
 540  continue
      outloc= 4HGAUS
      go to 520
c
c                       output location is node points
c
 550  continue
      outloc= 4HNODE
      go to 520
c
c                       output location is element center. only
c                       one set of values per element
c
 555  continue
      outloc= 4HCENT
      go to 520
c
c
c **********************************************************************
c *                                                                    *
c *                     input the configuration for output of stresses *
c *                     or strains. if the config. is not input, assume*
c *                     a default value for the elements in the list.  *
c *                     as this default value is not known at the      *
c *                     moment, an entity indicative of default        *
c *                     status has been previously stuffed in strcon.  *
c *                                                                    *
c **********************************************************************
c
c
c                       output stress type is cauchy stresses.
c                       this is the default for geonl.
c
 560  continue
      strcon= 4HCURR
      go to 520
c
c                       output stress type is 2nd_pk stresses.
c                       ** now obsolete feature unavailable **
c
 565  continue
      write(*,*) '>> * reference * is obsolete property. ignored...'
      go to 520
c
c
c **********************************************************************
c *                                                                    *
c *                     input output format. if the output format      *
c *                     is not input, assume a default value for the   *
c *                     elements in the list. as this default value is *
c *                     not known at the moment, an entity indicative  *
c *                     of default status has been previously stuffed  *
c *                     in outfmt.                                     *
c *                                                                    *
c **********************************************************************
c
c
c                       output format is short.
c
 570  continue
      outfmt= 4HSHRT
      go to 520
c
c                       output format is long.
c
 575  continue
      outfmt= 4HLONG
      go to 520
c
c
c **********************************************************************
c *                                                                    *
c *                     input the material for the elements in         *
c *                     the element list. a material type must be      *
c *                     given the first time through. on subsequent    *
c *                     lists, the previous lists' type will be used   *
c *                     as a default. if on the initial pass a list    *
c *                     is not given, ignore the current element list  *
c *                     and seek a new card. if the material specified *
c *                     does not exist in the material library, then   *
c *                     again ignore the current element list and      *
c *                     seek a new card.                               *
c *                                                                    *
c **********************************************************************
c
c
 580  continue
      if( label(dum) ) then
         name  = ' '
         mname = ' '
         call entits(name,nc)
         if( nc .gt. 24 ) nc = 24
         mname(1:nc) = name(1:nc)
         found = .false.
         matn = mathed/two16
 581     if( matn .eq. 32460 ) go to 582
            if( scanms(matnam(matn),mname,24) ) then
               matnum = matn
               found  = .true.
               go to 582
            end if
            matn = matlst(matn)/two16
         go to 581
 582     if( .not. found ) then
            call errmsg(29,dum,dums,dumr,dumd)
            mname(1:8) = 'nodefalt'
            go to 505
         end if
         defmat = .false.
c
c           Mark hack -- check to see if the type of this material is
c           11 (the interface material).  If it is, fake out the
c           element routines by reporting the macroscale material model
c           number here and setting a flag.
c
         if (matprp(9, matnum) .eq. 11) then
               macrointer = .true.
               minum = matnum
               matnum = imatprp(115, matnum)
         else
               macrointer = .false.
         end if
         go to 520
      else
         call errmsg(30,dum,dums,dumr,dumd)
         mname(1:8) = 'nodefalt'
         go to 505
      end if
c
c
c **********************************************************************
c *                                                                    *
c *                     element is geometrically linear.               *
c *                                                                    *
c **********************************************************************
c
c
 585  continue
      geonl= 4HFLSE
      go to 520
c
c
c **********************************************************************
c *                                                                    *
c *                     element is geometrically nonlinear.            *
c *                                                                    *
c **********************************************************************
c
c
 586  continue
      geonl= 4HTRUE
      go to 520
c
c
c **********************************************************************
c *                                                                    *
c *                     element uses/does not use b-bar                *
c *                     ** B-bar is the default formulation **         *
c *                                                                    *
c **********************************************************************
c
c
 587  continue
      bbar = 4HTRUE
      go to 520
 588  continue
      bbar = 4HFLSE
      go to 520
c
c
c **********************************************************************
c *                                                                    *
c *             for geometrically nonlinear  cohesive elements         *
c *             define the reference surface                           *
c *                                                                    *
c **********************************************************************
c
 595  continue
      if(matchs('top',3)) then
         surf = 4HO111
      else if(matchs('middle',3)) then
         surf = 4HO222
      else if(matchs('bottom',3)) then
         surf = 4HO333
      else
         call errmsg(318,dum,'surface',dumr,dumd)
      end if
      go to 520
c
c **********************************************************************
c *                                                                    *
c *             bar element area                                       *
c *                                                                    *
c **********************************************************************
c
 600  continue
      if( numr(area) ) then
         if( area <= 0.0e0 ) then
           call errmsg( 75, dum, dums, dumr, d umd )
           go to 520
         end if
         go to 520
      end if
      call errmsg( 85, dum, dums, dumr, dumd )
      go to 520
c
c **********************************************************************
c *                                                                    *
c *                     place element in temporary storage             *
c *                                                                    *
c **********************************************************************
c
c
c                       make sure that the material for the element
c                       list has not been omitted from input.
c
 900  continue
      if( scanms(mname,'nodefalt',8) ) then
         call errmsg(31,dum,dums,dumr,dumd)
         go to 505
      end if
c
c                       there is a valid list, an element type, a
c                       valid material, and either input or default
c                       other properties. place them in temporary
c                       element storage.
c
      call temsto( intlst, lenlst, intord, outloc, strcon, outfmt,
     &             defmat, deftyp, geonl, matnum, bbar, surf, type,
     &             macrointer, minum, area )
c
c                       all processed. examine another line for ele-
c                       ment data.
c
      go to 505
c
c
c **********************************************************************
c **********************************************************************
c
c
 9999 sbflg1= .true.
      sbflg2= .true.
c
c
      return
      end





