      subroutine sparse_factor( neqns, xlnz, coefs,
     &                          diag, iloc, inum, itemp, jtemp,
     &                          link, jloc, c1 )
#dbl      implicit double precision (a-h,o-z)
#sgl      implicit real (a-h,o-z)
      integer    xlnz(*)
      dimension  coefs(*), diag(*), inum(*), iloc(*),
     &           itemp(*), jtemp(*), link(*), jloc(*), c1(*)
      data zero / 0.0 /
c
      neq = neqns
      loopl = 4
      loop  = 6
      do i = 1 , neq
        link(i) = 0
        jtemp(i) = 0
        c1(i) = zero
        itemp(i) = 0
      end do
c 
      do i = 1, neq
        ii = inum(i)
          do j = xlnz(i), xlnz(i+1)-1
            jtemp(iloc(ii)) = jtemp(iloc(ii)) + 1
            ii = ii + 1
          end do
      end do
c 
c            ... loop over all rows
c 
      do i = 1, neq
      diagj = zero
          istop = jtemp(i)
          istart = xlnz(i)
          iend   = xlnz(i+1) - 1
          mtimes = iend - istart + 1
          newk   = link(i)
          me     = istop
          ms     = 0
          do nn = 1, istop
            kim = newk
            if( kim .ne. 0 )  newk = link(kim)
            lll = kim
            kk = xlnz(lll) + itemp(lll)
            ibegin = kk + 1
            ifinsh = xlnz(lll+1) - 1
            ntimes = ifinsh - ibegin + 1
            if( ntimes .eq. mtimes ) then
              ms = ms + 1
              jloc(ms) = lll
            else
              jloc(me) = lll
              me = me - 1
            endif
            if( ifinsh .lt. ibegin ) go to  894
            if( kim .ne. 0 ) then
                isub = iloc(inum(lll)+itemp(lll)+1)
                link(kim) = link(isub)
                link(isub) = kim
            end if
894         continue
          end do
c --------------------------------
          mss = (ms/loopl) * loopl
          do jj = 1, mss, loopl
           ll1 = jloc(jj)
           ll2 = jloc(jj+1)
           ll3 = jloc(jj+2)
           ll4 = jloc(jj+3)
           kk  = xlnz(ll1) + itemp(ll1)
           kk2 = xlnz(ll2) + itemp(ll2)
           kk3 = xlnz(ll3) + itemp(ll3)
           kk4 = xlnz(ll4) + itemp(ll4)
           ibegin = kk + 1
           ifinsh = xlnz(ll1+1) - 1
           ibegin2 = kk2 + 1
           ibegin3 = kk3 + 1
           ibegin4 = kk4 + 1
           div  =  coefs(kk)
           div2 =  coefs(kk2)
           div3 =  coefs(kk3)
           div4 =  coefs(kk4)
           div11=  coefs(kk)  * diag(ll1)
           div22=  coefs(kk2) * diag(ll2)
           div33=  coefs(kk3) * diag(ll3)
           div44=  coefs(kk4) * diag(ll4)
           diagj = diagj + div * div11   + div2 * div22
     &                                   + div3 * div33
     &                                   + div4 * div44
           if( ifinsh .lt. ibegin ) go to  997
                kkk = istart
                l2  = ibegin2
                l3  = ibegin3
                l4  = ibegin4
                do ll = ibegin, ifinsh
                  coefs(kkk) = coefs(kkk) -  coefs(ll) * div11
     &                                    -  coefs(l2) * div22
     &                                    -  coefs(l3) * div33
     &                                    -  coefs(l4) * div44
                  kkk = kkk + 1
                  l2  = l2  + 1
                  l3  = l3  + 1
                  l4  = l4  + 1
                end do
997        continue
           itemp(ll1) = itemp(ll1) + 1
           itemp(ll2) = itemp(ll2) + 1
           itemp(ll3) = itemp(ll3) + 1
           itemp(ll4) = itemp(ll4) + 1
          end do
c ---------------------------------
          do jj = mss+1, ms
           lll = jloc(jj)
           kk = xlnz(lll) + itemp(lll)
           ibegin = kk + 1
           ifinsh = xlnz(lll+1) - 1
           div =  coefs(kk)
           div11 =  coefs(kk) * diag(lll)
           diagj = diagj + div * div11
           if( ifinsh .lt. ibegin ) go to  998
           mm = inum(lll) + itemp(lll) + 1
           kkk = istart
           do ll = ibegin, ifinsh
              coefs(kkk) = coefs(kkk) -  coefs(ll) * div11
              kkk = kkk + 1
           end do
998        continue
           itemp(lll) = itemp(lll) + 1
          end do
c ----------------------------------
        jj = ms
        if( istop .eq. ms ) go to 800
 700    continue
        jj = jj + 1
        lll = jloc(jj)
        kk = xlnz(lll) + itemp(lll)
        ibegin = kk + 1
        ifinsh = xlnz(lll+1) - 1
        ntimes1 = ifinsh - ibegin
        div =  coefs(kk)
        div11 =  coefs(kk) * diag(lll)
        mm = inum(lll) + itemp(lll) + 1
c
c               loop level 6
c
        nt1 = mm
        nt2 = 2000000
        nt3 = 3000000
        nt4 = 4000000
        nt5 = 5000000
        nt6 = 6000000
        if( jj .lt. istop-(loop-2) ) then
            ll2 = jloc(jj+1)
            kk2 = xlnz(ll2) + itemp(ll2)
            ibegin2 = kk2 + 1
            ifinsh2 = xlnz(ll2+1) - 1
            ntimes2 = ifinsh2 - ibegin2
            if( ntimes2 .ne. ntimes1 ) go to 750
            nt2 = inum(ll2) + itemp(ll2) + 1
            ll3 = jloc(jj+2)
            kk3 = xlnz(ll3) + itemp(ll3)
            ibegin3 = kk3 + 1
            ifinsh3 = xlnz(ll3+1) - 1
            ntimes3 = ifinsh3 - ibegin3
            if( ntimes3 .ne. ntimes1 ) go to 750
            nt3 = inum(ll3) + itemp(ll3) + 1
            ll4 = jloc(jj+3)
            kk4 = xlnz(ll4) + itemp(ll4)
            ibegin4 = kk4 + 1
            ifinsh4 = xlnz(ll4+1) - 1
            ntimes4 = ifinsh4 - ibegin4
            if( ntimes4 .ne. ntimes1 ) go to 750
            nt4 = inum(ll4) + itemp(ll4) + 1
            ll5 = jloc(jj+4)
            kk5 = xlnz(ll5) + itemp(ll5)
            ibegin5 = kk5 + 1
            ifinsh5 = xlnz(ll5+1) - 1
            ntimes5 = ifinsh5 - ibegin5
            if( ntimes5 .ne. ntimes1 ) go to 750
            nt5 = inum(ll5) + itemp(ll5) + 1
            ll6 = jloc(jj+5)
            kk6 = xlnz(ll6) + itemp(ll6)
            ibegin6 = kk6 + 1
            ifinsh6 = xlnz(ll6+1) - 1
            ntimes6 = ifinsh6 - ibegin6
            if( ntimes6 .ne. ntimes1 ) go to 750
            nt6 = inum(ll6) + itemp(ll6) + 1
        end if
750     continue
c
        if( nt1 .eq. nt2 .and .nt2. eq. nt3 .and. nt3 .eq. nt4
     &      .and. nt4. eq. nt5 .and. nt5.eq. nt6 ) then
c 
            kk2 = xlnz(ll2) + itemp(ll2)
            div2 =  coefs(kk2)
            div22 =  coefs(kk2) * diag(ll2)
            ibegin2 = kk2 + 1
c       	  
            kk3 = xlnz(ll3) + itemp(ll3)
            div3 =  coefs(kk3)
            div33 =  coefs(kk3) * diag(ll3)
            ibegin3 = kk3 + 1
c       	  
            kk4 = xlnz(ll4) + itemp(ll4)
            div4 =  coefs(kk4)
            div44 =  coefs(kk4) * diag(ll4)
            ibegin4 = kk4 + 1
c       	  
            kk5 = xlnz(ll5) + itemp(ll5)
            div5 =  coefs(kk5)
            div55 =  coefs(kk5) * diag(ll5)
            ibegin5 = kk5 + 1
c       	  
            kk6 = xlnz(ll6) + itemp(ll6)
            div6 =  coefs(kk6)
            div66 =  coefs(kk6) * diag(ll6)
            ibegin6 = kk6 + 1
c       	  
            diagj = diagj + div * div11 + div2 * div22 + div3 * div33
     &                    + div4 * div44
     &                    + div5 * div55
     &                    + div6 * div66
            jj = jj+5
            if( ifinsh .lt. ibegin ) go to  996
            l2 = ibegin2
            l3 = ibegin3
            l4 = ibegin4
            l5 = ibegin5
            l6 = ibegin6
            do ll = ibegin, ifinsh
              kkk = iloc(mm)
              c1(kkk) = c1(kkk) +  coefs(ll) * div11
     &                          +  coefs(l2) * div22
     &                          +  coefs(l3) * div33
     &                          +  coefs(l4) * div44
     &                          +  coefs(l5) * div55
     &                          +  coefs(l6) * div66
              mm = mm + 1
              l2 = l2 + 1
              l3 = l3 + 1
              l4 = l4 + 1
              l5 = l5 + 1
              l6 = l6 + 1
            end do
c 
996         itemp(ll6) = itemp(ll6) + 1
            itemp(ll5) = itemp(ll5) + 1
            itemp(ll4) = itemp(ll4) + 1
            itemp(ll3) = itemp(ll3) + 1
            itemp(ll2) = itemp(ll2) + 1
            itemp(lll) = itemp(lll) + 1
c
         elseif( nt1 .eq. nt2 .and .nt2. eq. nt3 .and. nt3 .eq. nt4
     &          .and. nt4. eq. nt5 ) then
c 
            kk2 = xlnz(ll2) + itemp(ll2)
            div2 =  coefs(kk2)
            div22 =  coefs(kk2) * diag(ll2)
            ibegin2 = kk2 + 1
c 	 
            kk3 = xlnz(ll3) + itemp(ll3)
            div3 =  coefs(kk3)
            div33 =  coefs(kk3) * diag(ll3)
            ibegin3 = kk3 + 1
c 	 
            kk4 = xlnz(ll4) + itemp(ll4)
            div4 =  coefs(kk4)
            div44 =  coefs(kk4) * diag(ll4)
            ibegin4 = kk4 + 1
c 	 
            kk5 = xlnz(ll5) + itemp(ll5)
            div5 =  coefs(kk5)
            div55 =  coefs(kk5) * diag(ll5)
            ibegin5 = kk5 + 1
c 	 
            diagj = diagj + div * div11 + div2 * div22 + div3 * div33
     &                    + div4 * div44
     &                    + div5 * div55
            jj = jj + 4
            if( ifinsh .lt. ibegin ) go to  995
            l2 = ibegin2
            l3 = ibegin3
            l4 = ibegin4
            l5 = ibegin5
            do ll = ibegin, ifinsh
              kkk = iloc(mm)
              c1(kkk) = c1(kkk) +  coefs(ll) * div11
     &                          +  coefs(l2) * div22
     &                          +  coefs(l3) * div33
     &                          +  coefs(l4) * div44
     &                          +  coefs(l5) * div55
              mm = mm + 1
              l2 = l2 + 1
              l3 = l3 + 1
              l4 = l4 + 1
              l5 = l5 + 1
            end do
995         itemp(ll5) = itemp(ll5) + 1
            itemp(ll4) = itemp(ll4) + 1
            itemp(ll3) = itemp(ll3) + 1
            itemp(ll2) = itemp(ll2) + 1
            itemp(lll) = itemp(lll) + 1
c 
        elseif( nt1 .eq. nt2 .and .nt2. eq. nt3 .and.
     &          nt3 .eq. nt4 ) then
c 
            kk2 = xlnz(ll2) + itemp(ll2)
            div2 =  coefs(kk2)
            div22 =  coefs(kk2) * diag(ll2)
            ibegin2 = kk2 + 1
c 	  
            kk3 = xlnz(ll3) + itemp(ll3)
            div3 =  coefs(kk3)
            div33 =  coefs(kk3) * diag(ll3)
            ibegin3 = kk3 + 1
c 	  
            kk4 = xlnz(ll4) + itemp(ll4)
            div4 =  coefs(kk4)
            div44=  coefs(kk4) * diag(ll4)
            ibegin4 = kk4 + 1
c 	  
            diagj = diagj + div * div11 + div2 * div22 + div3 * div33
     &                    + div4 * div44
            jj = jj+3
            if( ifinsh .lt. ibegin ) go to  994
            l2 = ibegin2
            l3 = ibegin3
            l4 = ibegin4
            do ll = ibegin, ifinsh
              kkk = iloc(mm)
              c1(kkk) = c1(kkk) +  coefs(ll) * div11
     &                          +  coefs(l2) * div22
     &                          +  coefs(l3) * div33
     &                          +  coefs(l4) * div44
              mm = mm + 1
              l2 = l2 + 1
              l3 = l3 + 1
              l4 = l4 + 1
            end do
994         itemp(ll4) = itemp(ll4) + 1
            itemp(ll3) = itemp(ll3) + 1
            itemp(ll2) = itemp(ll2) + 1
            itemp(lll) = itemp(lll) + 1
c 
        elseif( nt1 .eq. nt2 .and .nt2. eq. nt3 ) then
c 
            kk2 = xlnz(ll2) + itemp(ll2)
            div2 =  coefs(kk2)
            div22 =  coefs(kk2) * diag(ll2)
            ibegin2 = kk2 + 1
c 	
            kk3 = xlnz(ll3) + itemp(ll3)
            div3 =  coefs(kk3)
            div33 =  coefs(kk3) * diag(ll3)
            ibegin3 = kk3 + 1
c 	
            diagj = diagj + div * div11 + div22 * div2 + div3 * div33
            jj = jj + 2
            if( ifinsh .lt. ibegin ) go to  993
            l2 = ibegin2
            l3 = ibegin3
            do ll = ibegin, ifinsh
              kkk = iloc(mm)
              c1(kkk) = c1(kkk) +  coefs(ll) * div11
     &                          +  coefs(l2) * div22
     &                          +  coefs(l3) * div33
              mm = mm + 1
              l2 = l2 + 1
              l3 = l3 + 1
            end do
 993        itemp(ll3) = itemp(ll3) + 1
            itemp(ll2) = itemp(ll2) + 1
            itemp(lll) = itemp(lll) + 1
c
          elseif( nt1 .eq. nt2 ) then
c 
            kk2 = xlnz(ll2) + itemp(ll2)
            div2 =  coefs(kk2)
            div22 =  coefs(kk2) * diag(ll2)
            ibegin2 = kk2 + 1
            diagj = diagj + div * div11 + div2 * div22
            jj = jj + 1
            if( ifinsh .lt. ibegin ) go to  992
            l2 = ibegin2
            do ll = ibegin, ifinsh
              kkk = iloc(mm)
              c1(kkk) = c1(kkk) +  coefs(ll) * div11
     &                          +  coefs(l2) * div22
              mm = mm + 1
              l2 = l2 + 1
            end do
992         itemp(ll2) = itemp(ll2) + 1
            itemp(lll) = itemp(lll) + 1
c 
          else
            diagj = diagj + div * div11
            do ll = ibegin, ifinsh
              kkk = iloc(mm)
              c1(kkk) = c1(kkk) +  coefs(ll) * div11
              mm = mm + 1
            end do
            if( ifinsh .lt. ibegin ) go to  991
991         itemp(lll) = itemp(lll) + 1
          end if
c 
        if( jj. ne. istop ) go to 700
800     continue
c
        diagj = diag(i) - diagj
        diagal = diagj
        diag(i) = diagal
        mm = inum(i)
        if( iend .ge. istart ) then
          isub = iloc(inum(i))
          link(i) = link(isub)
          link(isub) = i
        end if
        do kik = istart, iend
          isub = iloc(mm)
          coefs(kik) = ( coefs(kik) - c1(isub) )/diagal
          c1(isub) = zero
          mm = mm + 1
        end do
c 
      end do
 
c-------------------  end of factorization ---------------
c
      return
      end
 









