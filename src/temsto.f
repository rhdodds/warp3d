c     ****************************************************************
c     *                                                              *
c     *                      subroutine temsto                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 05/25/00                   *
c     *                                                              *
c     *     this subroutine sets up the elements in the given in-    *
c     *     teger list to have their properties placed in global     * 
c     *     storage, and calls the appropriate processing subroutine *
c     *     based on element type.                                   * 
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine temsto( intlst, lenlst, intord, outloc, strcon, outfmt,
     &                   defmat, deftyp, geonl, matnum, bbar, surf,
     &                   type, macrointer, minum )
      use main_data, only : elstor
      implicit integer (a-z)     
      include 'common.main'
      logical defmat, deftyp, macrointer
      real dumr
      double precision
     &     dumd
      character dums
      dimension intlst(*)               
c
c                       for each element in the list, set the
c                       element temporary storage array.
c                       
      icn    = 0
      iplist = 1
 20   call trxlst(intlst,lenlst,iplist,icn,elem)
c
c                       check that the list element is not negative.
c
      if( elem .lt. 0 ) then
         call errmsg(86,elem,dums,dumr,dumd)
         go to 30
      end if                              
c
c                       check that the list element does not exceed
c                       the number of elements in the structure.
c
      if( elem .gt. noelem ) then 
         call errmsg(35,elem,dums,dumr,dumd)
         go to 30
      end if
c
c                       store the material type. if the material for 
c                       this element has already been entered and the
c                       new material is a default value, then ignore
c                       the new material. default's cannot override
c                       previous data.
c                       
      if( elstor(1,elem) .ne. 0 ) then
         if( defmat ) go to 21
      end if
      elstor(2,elem) = matnum
c
c                       If we are actually faking out the material
c                       number for the interface damage model, store
c                       the "real" material number.  Else store -1
c                       (which will mean "normal material")
      if ( macrointer) then
            elstor(11,elem) = minum
      else
            elstor(11,elem) = -1
      end if
c
c
c                       store the integration order. if the order for
c                       this element has already been entered and the
c                       new order is a default value, then ignore
c                       the new order. default's cannot override
c                       previous data.
c                       
 21   if( elstor(1,elem) .ne. 0 ) then
         if( intord .eq. 4HDEFA ) go to 22
      end if
      elstor(3,elem) = intord
c
c                       store the output location. if the location for
c                       this element has already been entered and the
c                       new location is a default value, then ignore
c                       the new location. default's cannot override
c                       previous data.
c                       
 22   if( elstor(1,elem) .ne. 0 ) then
         if( outloc .eq. 4HDEFA ) go to 23
      end if
      elstor(4,elem) = outloc
c
c                       store the stress type. if the stress type for
c                       this element has already been entered and the
c                       new stress type is a default value, then ignore
c                       the new stress type. default's cannot override
c                       previous data.
c                       
 23   if( elstor(1,elem) .ne. 0 ) then
         if( strcon .eq. 4HDEFA ) go to 24
      end if
      elstor(5,elem) = strcon
c
c                       store the output format. if the output format for
c                       this element has already been entered and the
c                       new output format is a default value, then 
c                       ignore the new output format. default's cannot 
c                       override previous data.
c                       
 24   if( elstor(1,elem) .ne. 0 ) then
         if( outfmt .eq. 4HDEFA ) go to 25
      end if
      elstor(6,elem) = outfmt
c
c                       store the geometric nonlinearity flag. if this
c                       flag has already been stored for this element and
c                       the new flag is the default, then ignore the new flag.
c                       default's cannot override previous data.
c
 25   if( elstor(1,elem) .ne. 0 ) then
         if( geonl .eq. 4HDEFA ) go to 26
      end if
      elstor(7,elem) = geonl
c
c                       store the element type. if the type for this
c                       element has already been entered and the
c                       new type is a default value, then ignore
c                       the new type. default's cannot override
c                       previous data.
c                       
 26   if( elstor(1,elem) .ne. 0 ) then
         if(deftyp) go to 27
      end if
      elstor(1,elem) = type 
c
c                       store the bbar flag. if this
c                       flag has already been stored for this element and
c                       the new flag is the default, then ignore the new flag.
c                       default's cannot override previous data.
c
 27   if( elstor(1,elem) .ne. 0 ) then
         if( bbar .eq. 4HDEFA ) go to 28 
      end if
      elstor(8,elem) = bbar
c
c                       store the surface flag for nonlinear cohesive elements.
c                       if this flag has already been stored for this element
c                       and the new flag is the default, then ignore the new
c                       flag. default's cannot override previous data.
c
 28   if( elstor(1,elem) .ne. 0 ) then
         if( surf .eq. 4HDEFA ) go to 30
      end if
      elstor(10,elem) = surf
c
c
 30   if( iplist .ne. 0 ) go to 20
c
c
      return
      end

