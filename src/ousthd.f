c     ****************************************************************
c     *                                                              *
c     *                      subroutine ousthd                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 03/11/13 rhd               *
c     *                                                              *
c     *     prints the headers on a new page of output               *
c     *     put for output of stress or strain and their related     *
c     *     quantities.                                              *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine ousthd( pgnum, lnum, hedtyp, loctyp, strlbl, newhed,
     &                   nitm, nl, nstrou, fmtyp, noheader )
      use global_data ! old common.main
      implicit integer (a-z)
      logical newhed, noheader
      character(len=8) :: strlbl(*)
      character(len=*) :: hedtyp
      character(len=4) :: loctyp
      character(len=1) :: formfeed
c
      if ( noheader ) then
        newhed = .false.
        return
      end if
      formfeed = char(12)
c
c                       check to see if the page is full.
c
      if( lnum .gt. 55 ) then
c
c                       the page is full. reset line number and page
c                       number.
c
         newhed = .true.
         lnum   = 0
         pgnum  = pgnum+1
c
c                       print page headers.
c
         if(fmtyp.eq.1.or.fmtyp.eq.3) then
            write(out,900) formfeed,stname,pgnum,ltmstp,lsldnm,hedtyp
         else if(fmtyp.eq.2.or.fmtyp.eq.4) then
            write(out,905) formfeed,stname,pgnum,ltmstp,lsldnm,hedtyp
         end if
c
c
         if(nitm.gt.nstrou) then
            fnsh= nstrou
         else
            fnsh= nitm
         end if
c
         if(fmtyp.eq.1) then
            write(out,910) loctyp,(strlbl(i),i=1,fnsh)
         else if(fmtyp.eq.2) then
            write(out,915) loctyp,(strlbl(i),i=1,fnsh)
         else if(fmtyp.eq.3) then
            write(out,920) loctyp,(strlbl(i),i=1,fnsh)
         else if(fmtyp.eq.4) then
            write(out,925) loctyp,(strlbl(i),i=1,fnsh)
         end if
c
         cl= 2
 10      if(cl.gt.nl) go to 15
c
         strt= (cl-1)*nitm+1
         if(cl*nitm.gt.nstrou) then
            fnsh= nstrou
         else
            fnsh= cl*nitm
         end if
c
         if(fmtyp.eq.1) then
            write(out,911) (strlbl(i),i=strt,fnsh)
         else if(fmtyp.eq.2) then
            write(out,916) (strlbl(i),i=strt,fnsh)
         else if(fmtyp.eq.3) then
            write(out,921) (strlbl(i),i=strt,fnsh)
         else if(fmtyp.eq.4) then
            write(out,926) (strlbl(i),i=strt,fnsh)
         end if
c
         cl= cl+1
         go to 10
c
 15      lnum= lnum+10+(nl-1)
c
c
      else
c
c                       no new headers needed.
c
         newhed = .false.
c
      end if
c
c
 900  format(a1,///8x,'structure ',a8,51x,'page no. ',i4/8x,
     &       'step no. ',i7,24x,'loading ',a8//8x,a)
c
 905  format(a1,///8x,'structure ',a8,23x,'page no. ',i4/8x,
     &       'step no. ',i7,24x,'loading ',a8//8x,a)
c
 910  format(///6x,'elem',2x,a4,16x,a8,3(23x,a8))
 911  format(9x,4(23x,a8))
c
 915  format(///6x,'elem',2x,a4,16x,a8,23x,a8)
 916  format(9x,2(23x,a8))
c
 920  format(///6x,'elem',2x,a4,9x,a8,5(12x,a8))
 921  format(13x,6(12x,a8))
c
 925  format(///6x,'elem',2x,a4,9x,a8,2(12x,a8))
 926  format(13x,3(12x,a8))
c
      return
      end






