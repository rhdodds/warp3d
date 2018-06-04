c     ****************************************************************
c     *                                                              *
c     *                      subroutine inalpha                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/06/98                   *
c     *                                                              *
c     *     this subroutine supervises and conducts the input of     *
c     *     element thermal expansion coefficients for anisotopic    *
c     *     materials. values are stored directly into props         *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine inalpha( sbflg1, sbflg2 )
      use global_data ! old common.main
      implicit integer (a-z)
      double precision
     &   dumd
      real dumr, alphax, alphay, alphaz, alphaxy, alphayz, alphaxz,
     &     zero
      character(len=1) :: dums
      logical sbflg1, sbflg2
      logical matchs, matchs_exact, endcrd, true, numr
      dimension intlst(mxlsz)
      data zero / 0.0 /
c
c
c                       if we just left with a list not starting with an
c                       integerlist or user list and command was
c                       not recognized in driver, we're back here with
c                       sbflg1 = .true. - error message. get next line
c
c
      if( sbflg1 ) call errmsg(24,1,dums,dumr,dumd)
c
 505  continue
      alphax  = zero
      alphay  = zero
      alphaz  = zero
      alphaxy = zero
      alphayz = zero
      alphaxz = zero
      call readsc
      if ( matchs('dump',4) ) then
       call inalpha_dump
       go to 505
      end if

c
c                       translate the list of elements input on
c                       this line.
c
      call trlist( intlst, mxlsz, noelem, lenlst, errnum )
c
c                       branch on the return code from trlist. a
c                       value of 1 indicates no error. a value of
c                       2 indicates that the parse rules failed in
c                       the list. a value of 3 indicates that the
c                       list overflowed its maximum length of mxlsz.
c                       in these last two cases, the rest of the card
c                       will be ignored and a new card will be sought.
c                       a value of 4 indicates that no list was found.
c                       in this case, either thrm. expansion input has
c                       ceased.
c
      if( errnum .eq. 1 ) go to 511
      if( errnum .eq. 2 ) then
         param = 1
         call errmsg(24,param,dums,dumr,dumd)
         call scan_flushline; go to 9999
      else if( errnum .eq. 3 ) then
         param = 2
         call errmsg(24,param,dums,dumr,dumd)
         call scan_flushline; go to 9999
      else ! = 4 no list found
         call reset
         if( true(idum) ) call splunj
         go to 9999
      end if
c
c
c **********************************************************************
c *                                                                    *
c *           an element list exists. read the list of alpha values    *
c *           for these elements. read values until end of line        *
c *                                                                    *
c **********************************************************************
c
c
 511  continue
      if ( matchs_exact('alphaxy') ) then
        if ( .not. numr(alphaxy) ) then
         call errmsg( 5,dum,'a_xy',dumr,dumd )
        else
         go to 511
        end if
      end if
      if ( matchs_exact('alphaxz') ) then
        if ( .not. numr(alphaxz) ) then
         call errmsg( 5,dum,'a_xz',dumr,dumd )
        else
         go to 511
        end if
      end if
      if( matchs_exact('alphayz') ) then
        if ( .not. numr(alphayz) ) then
         call errmsg( 5,dum,'a_yz',dumr,dumd )
        else
         go to 511
        end if
      end if
      if( matchs_exact('alphaz') ) then
        if ( .not. numr(alphaz) ) then
         call errmsg( 5,dum,'a_z ',dumr,dumd )
        else
         go to 511
        end if
      end if
      if( matchs_exact('alphay') ) then
        if ( .not. numr(alphay) ) then
         call errmsg( 5,dum,'a_y ',dumr,dumd )
        else
         go to 511
        end if
      end if
      if( matchs_exact('alphax') ) then
        if ( .not. numr(alphax) ) then
         call errmsg( 5,dum,'a_x ',dumr,dumd )
        else
         go to 511
        end if
      end if
      if ( matchs_exact('xy') ) then
        if ( .not. numr(alphaxy) ) then
         call errmsg( 5,dum,'a_xy',dumr,dumd )
        else
         go to 511
        end if
      end if
      if ( matchs_exact('xz') ) then
        if ( .not. numr(alphaxz) ) then
         call errmsg( 5,dum,'a_xz',dumr,dumd )
        else
         go to 511
        end if
      end if
      if( matchs_exact('yz') ) then
        if ( .not. numr(alphayz) ) then
         call errmsg( 5,dum,'a_yz',dumr,dumd )
        else
         go to 511
        end if
      end if
      if( matchs_exact('z') ) then
        if ( .not. numr(alphaz) ) then
         call errmsg( 5,dum,'a_z ',dumr,dumd )
        else
         go to 511
        end if
      end if
      if( matchs_exact('y') ) then
        if ( .not. numr(alphay) ) then
         call errmsg( 5,dum,'a_y ',dumr,dumd )
        else
         go to 511
        end if
      end if
      if( matchs_exact('x') ) then
        if ( .not. numr(alphax) ) then
         call errmsg( 5,dum,'a_x ',dumr,dumd )
        else
         go to 511
        end if
      end if
      if( matchs(',',1) ) go to 511
c
c                       there is no match for expansion coeff.
c                       check for end of card. if not, print error
c                       message.
c
      if( endcrd(dum) ) go to 590
      call errmsg(211,dum,dums,dumr,dumd)
      call scan_flushline; go to 9999
c
c **********************************************************************
c *                                                                    *
c *                     a valid list and alphas for that list          *
c *                     have been input. store values for elements     *
c *                                                                    *
c **********************************************************************
c
 590  continue
      call instore_alpha( intlst, lenlst, alphax, alphay, alphaz,
     &                    alphaxy, alphayz, alphaxz )
c
c                       all processed. examine another card for ele-
c                       ment data.
c
      go to 505
c
c
 9999 continue
      sbflg1 = .true.
      sbflg2 = .true.
c
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine instore_alpha                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/10/98                   *
c     *                                                              *
c     *     this subroutine stores the anisotropic thermal expansion *
c     *     coefficients for a list of elements                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine instore_alpha( intlst, lenlst, alphax, alphay, alphaz,
     &                    alphaxy, alphayz, alphaxz )
      use global_data ! old common.main
      implicit integer (a-z)
      real dumr, alphax, alphay, alphaz, alphaxy, alphayz, alphaxz
      double precision
     &     dumd
      character :: dums
      dimension intlst(*)
c
c                       for each element in the list, set the
c                       element temporary storage array.
c
      icn    = 0
      iplist = 1
 20   call trxlst( intlst, lenlst, iplist, icn, elem )
c
c                       check that the list element is not negative.
c
      if( elem .lt. 0 ) then
         call errmsg( 86, elem, dums, dumr, dumd )
         go to 30
      end if
c
c                       check that the list element does not exceed
c                       the number of elements in the structure.
c
      if( elem .gt. noelem ) then
         call errmsg( 35, elem, dums, dumr, dumd )
         go to 30
      end if
c
c                       store the alpha, ij values, get next element
c                       in the list
c                       note order of terms !!
c
      props(9,elem)  = alphax
      props(13,elem) = alphay
      props(34,elem) = alphaz
      props(35,elem) = alphaxy
      props(36,elem) = alphayz
      props(37,elem) = alphaxz
c
 30   if( iplist .ne. 0 ) go to 20
c
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine inalpha_dumpa                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/10/98                   *
c     *                                                              *
c     *     this subroutine prints the anisotropic thermal expansion *
c     *     coefficients for a list of elements                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine inalpha_dump
      use global_data ! old common.main
      implicit integer (a-z)
      real  alphax, alphay, alphaz, alphaxy, alphayz, alphaxz
c
c                       for each element
c                       element temporary storage array.
c
      write(out,9100)
      do elem = 1, noelem
        alphax  = props(9,elem)
        alphay  = props(13,elem)
        alphaz  = props(34,elem)
        alphaxy = props(35,elem)
        alphayz = props(36,elem)
        alphaxz = props(37,elem)
        write(out,9000) elem, alphax, alphay, alphaz, alphaxy, alphayz,
     &                        alphaxz
      end do
c
 9000 format(2x,i7,6e14.4)
 9100 format(//,'Dump of thermal expansion coefficients...',
     & //,4x,'elem',8x,'alphax',7x,'alphay',7x,'alphaz',
     &     7x,'alphaxy',7x, 'alphayz',7x,'alphaxz')
c
      return
      end
