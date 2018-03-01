c     ****************************************************************
c     *                                                              *
c     *                      subroutine ourlhd                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 07/01//90                  *
c     *                                                              *
c     *     this subroutine prints the headers on a new page of out- *
c     *     put for output of residual loads.                        *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine ourlhd(pgnum,lnum,hedtyp,doflbl,ndof,newhed,ldnum,
     &                                           nitm,nl,fmtyp,iter)
      use global_data ! old common.main
      implicit integer (a-z)
      logical newhed
      character(len=8) :: doflbl(*),ldnam
      character(len=20) :: hedtyp
      character(len=1) :: formfeed
c
      formfeed = char(12)
c
c                       check to see if page is full.
c
      if(lnum.gt.55) then
c
c                       the page is full. reset line number and page
c                       number.
c
         newhed= .true.
         lnum= 0
         pgnum= pgnum+1
c
c                       set the current step number and its loading.
c
         stepno= ltmstp+1
         ldnam= lodnam(ldnum)
c
c                       print page headers.
c
         if(fmtyp.eq.1.or.fmtyp.eq.3) then
            write(out,900) formfeed,stname,pgnum,stepno,iter,ldnam,
     &                     hedtyp
         else if(fmtyp.eq.2.or.fmtyp.eq.4) then
            write(out,905) formfeed,stname,pgnum,stepno,iter,ldnam,
     &                     hedtyp
         end if
c
c
         if(nitm.gt.ndof) then
            fnsh= ndof
         else
            fnsh= nitm
         end if
c
         if(fmtyp.eq.1) then
            write(out,910) (doflbl(i),i=1,fnsh)
         else if(fmtyp.eq.2) then
            write(out,915) (doflbl(i),i=1,fnsh)
         else if(fmtyp.eq.3) then
            write(out,920) (doflbl(i),i=1,fnsh)
         else if(fmtyp.eq.4) then
            write(out,925) (doflbl(i),i=1,fnsh)
         end if
c
         cl= 2
 10      if(cl.gt.nl) go to 15
c
         strt= (cl-1)*nitm+1
         if(cl*nitm.gt.ndof) then
            fnsh= ndof
         else
            fnsh= cl*nitm
         end if
c
         if(fmtyp.eq.1) then
            write(out,911) (doflbl(i),i=strt,fnsh)
         else if(fmtyp.eq.2) then
            write(out,916) (doflbl(i),i=strt,fnsh)
         else if(fmtyp.eq.3) then
            write(out,921) (doflbl(i),i=strt,fnsh)
         else if(fmtyp.eq.4) then
            write(out,926) (doflbl(i),i=strt,fnsh)
         end if
c
         cl= cl+1
         go to 10
c
 15      lnum= lnum+10+(nl-1)
c
      else
c
c                       no new headers needed.
c
         newhed= .false.
c
      end if
c
c
 900  format(a1,///8x,'structure ',a8,51x,'page no. ',i4/8x,
     &       'step no. ',i7,4x,'iter. no. ',i7,4x,'loading ',a8//9x,
     &       'nodal ',a20)
c
 905  format(a1,///8x,'structure ',a8,23x,'page no. ',i4/8x,
     &       'step no. ',i7,4x,'iter. no. ',i7,4x,'loading ',a8//9x,
     &       'nodal ',a20)
c
 910  format(///8x,'node',12x,a8,3(21x,a8))
 911  format(24x,a8,3(21x,a8))
c
 915  format(///8x,'node',12x,a8,21x,a8)
 916  format(24x,a8,21x,a8)
c
 920  format(///8x,'node',6(12x,a8))
 921  format(12x,6(12x,a8))
c
 925  format(///8x,'node',3(12x,a8))
 926  format(12x,3(12x,a8))
c
      return
      end
