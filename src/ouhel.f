c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouhel                        *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 12/1/2017 rhd              *
c     *                                                              *
c     *     hard copy output and/or packet output of strain or       *
c     *     stress values for a single element                       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouhel( bele, elem, hedtyp, nnode, ngp, nodpts, nstrou,
     &                  wide, eform, prec, strlbl, pgnum, lnum, newel,
     &                  center_output, noheader, out_packet_now )
      use global_data ! old common.main
      use main_data, only : incmap, incid, packet_file_no
      use elblk_data, only : elestr
      implicit integer (a-z)
      logical :: wide, eform, prec, newhed, nodpts, newel,
     &           center_output, noheader, out_packet_now
      character(len=8) :: strlbl(*)
      character(len=*) :: hedtyp
      double precision :: small_tol, zero
      data small_tol, zero / 1.0d-50, 0.0d00 /
c
c                       local declarations
c
      character(len=4) :: loctyp
c
c                       set up paramters for format of printed output.
c                       compute number f physical lines of output
c
      if( wide .and. prec ) then
         nitm   = 4
         hfmtyp = 1
         if( eform ) then
            fmtyp = 1
         else
            fmtyp = 2
         end if
      else if( .not. wide .and. prec ) then
         nitm   = 2
         hfmtyp = 2
         if( eform ) then
            fmtyp = 3
         else
            fmtyp = 4
         end if
      else if( wide .and. .not. prec ) then
         nitm   = 6
         hfmtyp = 3
         if( eform ) then
            fmtyp = 5
         else
            fmtyp = 6
         end if
      else if( .not. wide .and. .not. prec ) then
         nitm   = 3
         hfmtyp = 4
         if( eform ) then
            fmtyp = 7
         else
            fmtyp = 8
         end if
      end if

      nl = nstrou/nitm
      if( (nstrou-nl*nitm) .gt. 0 ) nl = nl+1
c
c                       set location type for the header for this
c                       element, and the number of stress/strain
c                       points in this element.
c
      if( nodpts ) then
         loctyp = 'node'
         npts   = nnode
      else if( center_output ) then
         loctyp = 'cent'
         npts   = 1
      else
         loctyp  = 'g.p.'
         npts    = ngp
      end if
c
c                       zero out small values to prevent 3 digit
c                       exponents
c
      do strpt = 1, npts
       do j = 1, nstrou
        if( abs(elestr(bele,j,strpt)) .le. small_tol )
     &    elestr(bele,j,strpt) = zero
       end do
      end do

c
c                      for packet output write entries & return.
c                      for output at element nodes, use the structure
c                      node number rather than element node number.
c
      if( out_packet_now ) then
        do strpt = 1, npts
          optno = strpt
          if( nodpts ) optno = incid(incmap(elem)+strpt-1)
          if( center_output ) optno = 0
          write(packet_file_no) elem, optno, elestr(bele,1:nstrou,strpt)
        end do
        return
      end if
c
c                       print the stresses or strains, etc, for this
c                       element. for output at element nodes,
c                       print the structure node number rather than
c                       the element node number.
c
      do strpt = 1, npts
c
      optno = strpt
      if( nodpts ) optno = incid(incmap(elem)+strpt-1)
      if( center_output ) optno = 1
c
c                       start new page if needed
c
        call ousthd( pgnum, lnum, hedtyp, loctyp, strlbl, newhed,
     &               nitm, nl, nstrou, hfmtyp, noheader )
c
c                       print the first line of stress/strain output.
c
        fnsh = nitm
        if( nitm .gt. nstrou ) fnsh = nstrou
c
        if( fmtyp .eq. 1 ) then
           if( newhed .or. newel ) then
              write(out,910) elem,optno,elestr(bele,1:fnsh,strpt)
           else
              write(out,911) optno,elestr(bele,1:fnsh,strpt)
           end if
        else if( fmtyp .eq. 2 ) then
           if( newhed .or. newel ) then
              write(out,912) elem,optno,elestr(bele,1:fnsh,strpt)
           else
              write(out,913) optno,elestr(bele,1:fnsh,strpt)
           end if
        else if( fmtyp .eq. 3 ) then
           if( newhed .or. newel ) then
              write(out,914) elem,optno,elestr(bele,1:fnsh,strpt)
           else
              write(out,915) optno,elestr(bele,1:fnsh,strpt)
           end if
        else if( fmtyp .eq. 4 ) then
           if(newhed.or.newel) then
              write(out,916) elem,optno,elestr(bele,1:fnsh,strpt)
           else
              write(out,917) optno,elestr(bele,1:fnsh,strpt)
           end if
        else if( fmtyp .eq. 5 ) then
           if( newhed .or. newel ) then
              write(out,918) elem,optno,elestr(bele,1:fnsh,strpt)
           else
              write(out,919) optno,elestr(bele,1:fnsh,strpt)
           end if
        else if( fmtyp .eq. 6 ) then
           if( newhed .or. newel ) then
              write(out,920) elem,optno,elestr(bele,1:fnsh,strpt)
           else
              write(out,921) optno,elestr(bele,1:fnsh,strpt)
           end if
        else if( fmtyp .eq. 7 ) then
           if( newhed .or. newel ) then
              write(out,922) elem,optno,elestr(bele,1:fnsh,strpt)
           else
              write(out,923) optno,elestr(bele,1:fnsh,strpt)
           end if
        else if( fmtyp .eq. 8 ) then
           if( newhed .or. newel ) then
              write(out,924) elem,optno,elestr(bele,1:fnsh,strpt)
           else
              write(out,925) optno,elestr(bele,1:fnsh,strpt)
           end if
        end if
c
        if( newel ) newel = .false.
        lnum = lnum + 1
c
c                    final lines as needed
c
        cl = 2
        if( cl .gt. nl ) cycle
        do cl = 2, nl
c
         lnum = lnum + 1
         call ousthd( pgnum, lnum, hedtyp, loctyp, strlbl, newhed,
     &                nitm, nl, nstrou, hfmtyp, noheader )
c
         strt = (cl-1)*nitm+1
         fnsh = cl*nitm
         if( cl*nitm .gt. nstrou ) fnsh = nstrou
c
         if( fmtyp .eq. 1 ) then
            if( newhed .or. newel ) then
               write(out,930) elem,elestr(bele,strt:fnsh,strpt)
            else
               write(out,931) elestr(bele,strt:fnsh,strpt)
            end if
         else if( fmtyp .eq. 2 ) then
            if( newhed .or. newel ) then
               write(out,932) elem,elestr(bele,strt:fnsh,strpt)
            else
               write(out,933) elestr(bele,strt:fnsh,strpt)
            end if
         else if( fmtyp .eq. 3 ) then
            if( newhed .or. newel ) then
               write(out,934) elem,elestr(bele,strt:fnsh,strpt)
            else
               write(out,935) elestr(bele,strt:fnsh,strpt)
            end if
         else if( fmtyp .eq. 4 ) then
            if( newhed .or. newel ) then
               write(out,936) elem,elestr(bele,strt:fnsh,strpt)
            else
               write(out,937) elestr(bele,strt:fnsh,strpt)
            end if
         else if( fmtyp .eq. 5 ) then
            if( newhed .or. newel ) then
               write(out,938) elem,elestr(bele,strt:fnsh,strpt)
            else
               write(out,939) elestr(bele,strt:fnsh,strpt)
            end if
         else if( fmtyp .eq. 6 ) then
            if( newhed .or. newel ) then
               write(out,940) elem,elestr(bele,strt:fnsh,strpt)
            else
               write(out,941) elestr(bele,strt:fnsh,strpt)
            end if
         else if( fmtyp .eq. 7 ) then
            if( newhed .or. newel ) then
               write(out,942) elem,elestr(bele,strt:fnsh,strpt)
            else
               write(out,943) elestr(bele,strt:fnsh,strpt)
            end if
         else if( fmtyp .eq. 8 ) then
            if( newhed .or. newel ) then
               write(out,944) elem,elestr(bele,strt:fnsh,strpt)
            else
               write(out,945) elestr(bele,strt:fnsh,strpt)
            end if
         end if
c
        end do
c
      end do
c
      return
c
 910  format(/4x,i7,1x,i6,4(3x,e26.19))
 911  format(/12x,i6,4(3x,e26.19))
 930  format(/4x,i7,7x,4(3x,e26.19))
 931  format(18x,4(3x,e26.19))
c
 912  format(/4x,i7,1x,i6,4(3x,f26.16))
 913  format(/12x,i6,4(3x,f26.16))
 932  format(/4x,i7,7x,4(3x,f26.16))
 933  format(18x,4(3x,f26.16))
c
 914  format(/4x,i7,1x,i6,2(3x,e26.19))
 915  format(/12x,i6,2(3x,e26.19))
 934  format(/4x,i7,7x,2(3x,e26.19))
 935  format(18x,2(3x,e26.19))
c
 916  format(/4x,i7,1x,i6,2(3x,f26.16))
 917  format(/12x,i6,2(3x,f26.16))
 936  format(/4x,i7,7x,2(3x,f26.16))
 937  format(18x,2(3x,f26.16))
c
 918  format(/4x,i7,1x,i6,4x,e12.5,5(8x,e12.5))
 919  format(/12x,i6,4x,e12.5,5(8x,e12.5))
 938  format(/4x,i7,3x,6(8x,e12.5))
 939  format(14x,6(8x,e12.5))
c
 920  format(/4x,i7,1x,i6,4x,f12.6,5(8x,f12.6))
 921  format(/12x,i6,4x,f12.6,5(8x,f12.6))
 940  format(/4x,i7,3x,6(8x,f12.6))
 941  format(14x,6(8x,f12.6))
c
 922  format(/4x,i7,1x,i6,4x,e12.5,2(8x,e12.5))
 923  format(/12x,i6,4x,e12.5,2(8x,e12.5))
 942  format(/4x,i7,3x,3(8x,e12.5))
 943  format(14x,3(8x,e12.5))
c
 924  format(/4x,i7,1x,i6,4x,f12.6,2(8x,f12.6))
 925  format(/12x,i6,4x,f12.6,2(8x,f12.6))
 944  format(/4x,i7,3x,3(8x,f12.6))
 945  format(14x,3(8x,f12.6))
c
      end












