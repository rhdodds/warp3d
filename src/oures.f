c     ****************************************************************
c     *                                                              *
c     *                      subroutine oures                        *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/1/2017 rhd              *
c     *                                                              *
c     *     this subroutine performs the output of the residual load *
c     *     vector for the current iteration in the solution of the  *
c     *     (n+1) configuration.                                     *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine oures( ldnum, iter )
      use global_data ! old common.main
c
      use main_data, only : trn, trnmat, inverse_incidences,
     &                      inverse_dof_map
c
      implicit integer (a-z)
      double precision
     &    eres(mxvl,mxndof), zero, trnmte(mxvl,mxedof,mxndof)
      logical newhed, trne(mxvl,mxndel)
      character(len=8) :: doflbl(mxndof)
      character(len=20) :: hedtyp
      data zero /0.0/
c
c                       check that this iteration is to be output.
c
      icn    = 0
      iplist = 1
 5    call trxlst(prslst,nprs,iplist,icn,outitr)
      if( iter .eq. outitr ) go to 6
      if( iplist .ne. 0 ) go to 5
      go to 9999
c
c                       initialize parameters controlling output.
c
 6    lnum   = 56
      pgnum  = 0
      lbltyp = 0
      do i = 1, mxndof
         doflbl(i) = ' '
      end do
c
c                       loop over the number of nodes
c
      do 15 nod = 1, nonode
c
         elem  = inverse_incidences(nod)%element_list(1)
         elnod = inverse_dof_map(nod)%edof_table(1,1)
         type  = iprops(1,elem)
         ndof  = iprops(4,elem)
c
c                       set the number of lines to be output per node.
c                       assume wide output, and precision e format.
c
c
          nitm   = 4
          hfmtyp = 1
          fmtyp  = 1
          nl     = ndof/nitm
          if( (ndof-nl*nitm) .gt. 0 ) nl = nl+1
c
c                       set the dof labels for the structural nodes,
c                       if they have not already been established. if
c                       they have been established, make sure that
c                       the current node is compatable. if not, print
c                       a warning and change the dof labels.
c
         call oulbir( 2, hedtyp, lbltyp, type, elem, doflbl )
c
c                       check for new page. if so, print page headers.
c
         call ourlhd(pgnum,lnum,hedtyp,doflbl,ndof,newhed,ldnum,nitm,
     &               nl,hfmtyp,iter)
c
c                       set the temporary vector of residual loads to
c                       be printed.
c
         do dof = 1, ndof
            if( cstmap(dstmap(nod)+dof-1) .eq. 0 ) then
               eres(1,dof) = res(dstmap(nod)+dof-1)
            else
               eres(1,dof) = zero
            end if
         end do
c
c                       transform temporary vector at the current
c                       node to uniform global coordinates.
c
         trne(1,1) = trn(nod)
         if( trne(1,1) ) then
            trnmte(1,1,1) = trnmat(nod)%mat(1,1)
            trnmte(1,1,2) = trnmat(nod)%mat(1,2)
            trnmte(1,1,3) = trnmat(nod)%mat(1,3)
            trnmte(1,2,1) = trnmat(nod)%mat(2,1)
            trnmte(1,2,2) = trnmat(nod)%mat(2,2)
            trnmte(1,2,3) = trnmat(nod)%mat(2,3)
            trnmte(1,3,1) = trnmat(nod)%mat(3,1)
            trnmte(1,3,2) = trnmat(nod)%mat(3,2)
            trnmte(1,3,3) = trnmat(nod)%mat(3,3)
            call trnvec( eres, trnmte, trne, ndof, 1, 1, 2 )
c
         end if
c
c                       print the first line of nodal res. ld. output.
c
         if( nitm.gt. ndof ) then
            fnsh = ndof
         else
            fnsh = nitm
         end if
c
         if( fmtyp .eq. 1 ) then
            write(out,915) nod,(eres(1,i),i=1,fnsh)
         else if( fmtyp .eq. 2 ) then
            write(out,917) nod,(eres(1,i),i=1,fnsh)
         else if( fmtyp .eq. 3 ) then
            write(out,920) nod,(eres(1,i),i=1,fnsh)
         else if(fmtyp.eq.4) then
            write(out,922) nod,(eres(1,i),i=1,fnsh)
         else if(fmtyp.eq.5) then
            write(out,925) nod,(eres(1,i),i=1,fnsh)
         else if(fmtyp.eq.6) then
            write(out,927) nod,(eres(1,i),i=1,fnsh)
         else if(fmtyp.eq.7) then
            write(out,930) nod,(eres(1,i),i=1,fnsh)
         else if(fmtyp.eq.8) then
            write(out,932) nod,(eres(1,i),i=1,fnsh)
         end if
c
         lnum = lnum + 1
c
c                       print subsequent lines of nodal r.l. output,
c                       if necessary.
c
         cl = 2
 30      if( cl .gt. nl ) go to 15
c
         lnum = lnum + 1
         call ourlhd(pgnum,lnum,hedtyp,doflbl,ndof,newhed,ldnum,nitm,
     &                                                nl,hfmtyp,iter)
         if( newhed ) write(out,945)
c
         strt= (cl-1)*nitm+1
         if( cl*nitm .gt. ndof ) then
            fnsh = ndof
         else
            fnsh = cl*nitm
         end if
c
         if(fmtyp.eq.1) then
            write(out,916) (eres(1,i),i=strt,fnsh)
         else if(fmtyp.eq.2) then
            write(out,918) (eres(1,i),i=strt,fnsh)
         else if(fmtyp.eq.3) then
            write(out,921) (eres(1,i),i=strt,fnsh)
         else if(fmtyp.eq.4) then
            write(out,923) (eres(1,i),i=strt,fnsh)
         else if(fmtyp.eq.5) then
            write(out,926) (eres(1,i),i=strt,fnsh)
         else if(fmtyp.eq.6) then
            write(out,928) (eres(1,i),i=strt,fnsh)
         else if(fmtyp.eq.7) then
            write(out,931) (eres(1,i),i=strt,fnsh)
         else if(fmtyp.eq.8) then
            write(out,933) (eres(1,i),i=strt,fnsh)
         end if
c
         cl = cl + 1
         go to 30
c
 15   continue
c
c
 915  format(/6x,i7,4(2x,e26.19))
 916  format(13x,4(3x,e26.19))
 917  format(/6x,i7,4(2x,f26.16))
 918  format(13x,4(3x,f26.16))
c
 920  format(/6x,i7,2(2x,e26.19))
 921  format(13x,2(3x,e26.19))
 922  format(/6x,i7,2(2x,f26.16))
 923  format(13x,2(3x,f26.16))
c
 925  format(/6x,i7,6(78x,e12.5))
 926  format(13x,6(8x,e12.5))
 927  format(/6x,i7,6(7x,f12.6))
 928  format(13x,6(8x,f12.6))
c
 930  format(/6x,i7,3(7x,e12.5))
 931  format(13x,3(8x,e12.5))
 932  format(/6x,i7,3(7x,f12.6))
 933  format(13x,3(8x,f12.6))
c
 945  format(1x,' ')
c
 9999 return
      end
