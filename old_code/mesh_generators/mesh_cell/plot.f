c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine create_plotfil(nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3, 
     &       noda,na1,na2,na3,  nodb,nb1,nb2,nb3, job,jobh)
C-------------------------------------------------------C
C   Rutinen skapar 5 st plotfiler for HPPLOTERN .       C
C   1) Plotfil av sprickplanet.                         C
C   2) Plotfil av symmetri snittet.                     C
C   3) Plotfil av det bakre planet.                     C
C   4) Plotfil av ovansidan .                           C
C   5) Plotfil av undersidan.                           C
C-------------------------------------------------------C
      implicit none
c
      include 'common_nod.f'
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  nm
      parameter (nm=2000)
      integer jobh,n(nm),ia1,ia2,i,j,k,l,i1,di,dj,kend,io,ksr2
c
      double precision sf,z1,z2,x1,x2,da
c
      character job*40, fildat*40,filinf*40
c
      logical  left
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
c---------------------------------
c... (1) Plot file 1 ( crack plane, Y = 0)
c---------------------------------
c
c        write(*,'(t1,a,3i5)') ' ima, jma, kma : ',ima,jma,kma
c        write(*,'(t1,a,3i5)') ' na1, na2, na3 : ',na1,na2,na3
c
	do i=1, nm
	   n(i)=0
	enddo
	fildat = job(1:jobh)//'_1.dat'
        io = 31
	open(unit=io,file=fildat,status='unknown')
C . . . a) RAMEN !
        l=1
	n(1)=5
	write(io,*) 0.0, 0.0
	write(io,*) 0.0, t
	write(io,*) w  , t
	write(io,*) w  , 0.0
	write(io,*) 0.0, 0.0
C...rho-direction; Zone A
	ia1=2*m1+1
	ia2=2*(m1+mh+mh)+1
	do j=3, jma-2, 2
	   l=l+1
	   if (mod(j,4).eq.1) then
	      i1=1
	   else
	      i1=3
	   endif
	   do i=i1, ia1, 2
	      if (noda(i,j,1).gt.0) then
	         n(l)=n(l)+1
	         write(io,*) npos(noda(i,j,1),1), npos(noda(i,j,1),3)
	      endif
	   enddo
	   l=l+1
	   do i=ia2, ima, 2
	      n(l)=n(l)+1
	      write(io,*) npos(noda(i,j,1),1), npos(noda(i,j,1),3)
	   enddo
	enddo
C...rho-direction; Zone C & S
        do j=3, jms-2, 2
           if (sjred_type.eq.1) then
              kend = kms
           elseif (sjred_type.eq.2) then
              if (mod(j,4).eq.1) then
                 kend = kms
              else
                 kend = kms-2
              endif
           elseif (sjred_type.eq.3) then
              if (mod(j,6).eq.1) then
                 kend = kms
              else
                 kend = kms-2
              endif
           endif
	   l=l+1
           i = ims
           do k=kend, 1, -2
              n(l)=n(l)+1
              write(io,*) npos(nods(i,j,k),1),npos(nods(i,j,k),3)
           enddo
           l = l + 1
           n(l)=n(l)+2
           write(io,*) npos(nodc(1,j,3),1),npos(nodc(1,j,3),3)
           write(io,*) npos(nodc(3,j,3),1),npos(nodc(3,j,3),3)
           i = 1
           k = 1
           do i=5, imc, 2
              n(l)=n(l)+1
              write(io,*) npos(nodc(i,j,k),1),npos(nodc(i,j,k),3)
           enddo
           i = 1
           l = l + 1
           do k=1, kend, 2
              n(l)=n(l)+1
              write(io,*) npos(nods(i,j,k),1),npos(nods(i,j,k),3)
           enddo
        enddo
c
c...phi-direction; Zone A
        l=l+1
	if (mod(na,4).ne.0) then
	   left=.true.
	   write(io,*) npos(noda(3,1,3),1), npos(noda(3,1,3),3)
	else
	   left=.false.
	   write(io,*) npos(noda(1,1,3),1), npos(noda(1,1,3),3)
	endif
        n(l)=n(l)+1
	do j=3, jma-4, 4
	   if (left) then
	      write(io,*) npos(noda(3,j,3),1), npos(noda(3,j,3),3)
	      write(io,*) npos(noda(1,j+2,3),1), npos(noda(1,j+2,3),3)
              n(l)=n(l)+2
	      left=.false.
	   else
	      write(io,*) npos(noda(3,j,3),1), npos(noda(3,j,3),3)
	      write(io,*) npos(noda(3,j+2,3),1), npos(noda(3,j+2,3),3)
              n(l)=n(l)+2
	      left=.true.
	   endif
	enddo
	do i=5, ima-2, 2
	  if (noda(i,1,1).gt.0) then
	     l=l+1
	     do j=1, jma, 2
	        if (noda(i,j,1).gt.0) then
	          write(io,*) npos(noda(i,j,1),1), npos(noda(i,j,1),3)
		  n(l)=n(l)+1
	        endif
	     enddo
	  endif
	enddo
c
c...phi-direction; Zone S
        do k=1, kms-4, 2
           i = ims
	   l=l+1
           do j=1, jms, 2
              n(l)=n(l)+1
              write(io,*) npos(nods(i,j,k),1), npos(nods(i,j,k),3)
           enddo
           i = 1
	   l=l+1
           do j=1, jms, 2
              n(l)=n(l)+1
              write(io,*) npos(nods(i,j,k),1), npos(nods(i,j,k),3)
           enddo
        enddo
c
        if (sjred_type.eq.1) then
           k = kms - 2
           i = ims
           l=l+1
           do j=1, jms, 2
              n(l)=n(l)+1
              write(io,*) npos(nods(i,j,k),1), npos(nods(i,j,k),3)
           enddo
           i = 1
           l=l+1
           do j=1, jms, 2
              n(l)=n(l)+1
              write(io,*) npos(nods(i,j,k),1), npos(nods(i,j,k),3)
           enddo
        elseif (sjred_type.eq.2) then
           k = kms - 2
           i = ims
           l=l+1
           n(l)=n(l)+1
           write(io,*) npos(nods(i,1,k),1), npos(nods(i,1,k),3)
           do j=5, jms, 8
              n(l)=n(l)+2
              write(io,*) npos(nods(i,j-2,k),1), npos(nods(i,j-2,k),3)
              write(io,*) npos(nods(i,j,k+2),1), npos(nods(i,j,k+2),3)
              if (j.le.(jms-4)) then
                 n(l)=n(l)+2
                 write(io,*) npos(nods(i,j+2,k),1),npos(nods(i,j+2,k),3)
                 write(io,*) npos(nods(i,j+4,k),1),npos(nods(i,j+4,k),3)
              endif
           enddo
c
           i = 1
           l=l+1
           n(l)=n(l)+1
           write(io,*) npos(nods(i,1,k),1), npos(nods(i,1,k),3)
           do j=5, jms, 8
              n(l)=n(l)+2
              write(io,*) npos(nods(i,j-2,k),1), npos(nods(i,j-2,k),3)
              write(io,*) npos(nods(i,j,k+2),1), npos(nods(i,j,k+2),3)
              if (j.le.(jms-4)) then
                 n(l)=n(l)+2
                 write(io,*) npos(nods(i,j+2,k),1),npos(nods(i,j+2,k),3)
                 write(io,*) npos(nods(i,j+4,k),1),npos(nods(i,j+4,k),3)
              endif
           enddo
        elseif (sjred_type.eq.3) then
           l=l+1
           k = kms - 2
           i = ims
           n(l)=n(l)+1
           write(io,*) npos(nods(i,1,k+2),1), npos(nods(i,1,k+2),3)
           do j=7, jms, 6
              n(l)=n(l)+3
              write(io,*) npos(nods(i,j-4,k),1), npos(nods(i,j-4,k),3)
              write(io,*) npos(nods(i,j-2,k),1), npos(nods(i,j-2,k),3)
              write(io,*) npos(nods(i,j,k+2),1), npos(nods(i,j,k+2),3)
           enddo
c
           l=l+1
           k = kms - 2
           i = 1
           n(l)=n(l)+1
           write(io,*) npos(nods(i,1,k+2),1), npos(nods(i,1,k+2),3)
           do j=7, jms, 6
              n(l)=n(l)+3
              write(io,*) npos(nods(i,j-4,k),1), npos(nods(i,j-4,k),3)
              write(io,*) npos(nods(i,j-2,k),1), npos(nods(i,j-2,k),3)
              write(io,*) npos(nods(i,j,k+2),1), npos(nods(i,j,k+2),3)
           enddo
        endif
c
c...phi-direction; Zone C
        l = l + 1
        do j=1, jmc, 2
           n(l)=n(l)+1
           write(io,*) npos(nodc(3,j,3),1), npos(nodc(3,j,3),3)
        enddo
        k = 1
        do i=5, imc-2, 2
           l = l + 1
           do j=1, jmc, 2
              n(l)=n(l)+1
              write(io,*) npos(nodc(i,j,k),1), npos(nodc(i,j,k),3)
           enddo
        enddo
c
c...Z-direction; Zone B
	do i=1, imb-2, 2
	   l=l+1
	   do j=1, jmb, 2
	      write(io,*) npos(nodb(i,j,1),1), npos(nodb(i,j,1),3)
	      n(l)=n(l)+1
	   enddo
	enddo
c...X-direction; Zone B
	do j=3, jmb-2, 2
	   l=l+1
	   do i=1, imb, 2
	      write(io,*) npos(nodb(i,j,1),1), npos(nodb(i,j,1),3)
	      n(l)=n(l)+1
	   enddo
	enddo
	close(io)
c
c...Generate the information file
c
	filinf = job(1:jobh)//'_1.inf'
	open(unit=32,file=filinf,status='unknown')
        write(32,'(t1,a)') '*CURVE1'
 	if (t.gt.w) then
	   sf=t/20.d0
	else
	   sf=w/20.d0
	endif
	write(32,'(t1,2f5.1,tr2,6(tr1,g10.4))')
     &     3.0,3.0, -sf,21.*sf,-sf,21.*sf, 18.0,18.0
	write(32,'(t1,2(a,g10.4,tr2,i1),tr2,f4.2)')
     &     '0.0 ', w/2., 1, ' 0.0 ', t, 1, 0.35
	write(32,'(t1,i3)') l
	do i=1, l
	   write(32,'(t1,i4,a)') n(i),' 3 0 401 0 1'
	enddo
	write(32,'(t1,a)') '*FINISH'
	close(32)
c
c  Zoom into the crack front region
	filinf = job(1:jobh)//'_1z.inf'
	open(unit=32,file=filinf,status='unknown')
        write(32,'(t1,a)') '*CURVE1'
        sf = max( t, npos(nods(1,jms,kms),1) ) / 20.0
	write(32,'(t1,2f5.1,tr2,6(tr1,g10.4))')
     &    3.0,3.0, -sf,21.*sf,-sf,21.*sf, 18.0,18.0
	write(32,'(t1,2(a,g10.4,tr2,i1),tr2,f4.2)')
     &   '0.0 ', (t+2.*sf)/4, 1, ' 0.0 ', (t+2.*sf)/4, 1, 0.35
	write(32,'(t1,i3)') l
	do i=1, l
	   write(32,'(t1,i4,a)') n(i),' 3 0 401 0 1'
	enddo
	write(32,'(t1,a)') '*FINISH'
	close(32)
c
c-----------------------------------
c... (2) Plot file 2 (symmetry plane, X=0)
c-----------------------------------
c
	do i=1, nm
	   n(i)=0
	enddo
	l=0
	fildat = job(1:jobh)//'_2.dat'
	filinf = job(1:jobh)//'_2.inf'
        io =31
	open(unit=io,file=fildat,status='unknown')
c
c...Zone C:
c
        j = 1
        l=l+1
        do i=5, imc, 2
	   n(l)=n(l)+1
           write(io,*) -npos(nodc(i,j,1),2), npos(nodc(i,j,1),3)
        enddo
        do k=3, kmc-4, 2
           l=l+1
           do i=1, imc, 2
              n(l)=n(l)+1
              write(io,*) -npos(nodc(i,j,k),2), npos(nodc(i,j,k),3)
           enddo
        enddo
        k = kmc-2
        left = .true.
        l=l+1
        do i=3, imc-3, 4
         if (left) then
            n(l)=n(l)+3
            write(io,*) -npos(nodc(i-2,j,k),2), npos(nodc(i-2,j,k),3)
            write(io,*) -npos(nodc(i,j,k),2), npos(nodc(i,j,k),3)
            write(io,*) -npos(nodc(i+2,j,k+2),2),npos(nodc(i+2,j,k+2),3)
            left = .false.
         else
            n(l)=n(l)+1
            write(io,*) -npos(nodc(i,j,k),2),npos(nodc(i,j,k),3)
            left = .true.
         endif
      enddo
      n(l)=n(l)+1
      write(io,*) -npos(nodc(imc,j,kmc),2),npos(nodc(imc,j,kmc),3)
c
      do i=3, imc-2, 2
         l=l+1
         if (mod(i,4).eq.1) then
            do k=1, kmc, 2
               if (nodc(i,j,k).gt.0) then
                  n(l)=n(l)+1
                  write(io,*) -npos(nodc(i,j,k),2), npos(nodc(i,j,k),3)
               endif
            enddo
         else
            do k=1, kmc-2, 2
               if (nodc(i,j,k).gt.0) then
                  n(l)=n(l)+1
                  write(io,*) -npos(nodc(i,j,k),2), npos(nodc(i,j,k),3)
               endif
            enddo
         endif
      enddo
c
c...Zone S:
c
      ksr2 = ksr1 + 2
      j = 1
c
      if (sfred_type.eq.1) then
c
         do i=1, ims, 2
            l=l+1
            do k=1, kms, 2
               n(l)=n(l)+1
               write(io,*) -npos(nods(i,j,k),2), npos(nods(i,j,k),3)
            enddo
         enddo
c
         do k=1, kms, 2
            l=l+1
            do i=1, ims, 2
               n(l)=n(l)+1
               write(io,*) -npos(nods(i,j,k),2), npos(nods(i,j,k),3)
            enddo
         enddo
c
      elseif (sfred_type.eq.2) then
c
         do i=1, ims, 2
            l=l+1
            if (mod(i,4).eq.1) then
               do k=1, kms, 2
                  n(l)=n(l)+1
                  write(io,*) -npos(nods(i,j,k),2), npos(nods(i,j,k),3)
               enddo
            else
               do k=1, ksr1, 2
                  n(l)=n(l)+1
                  write(io,*) -npos(nods(i,j,k),2), npos(nods(i,j,k),3)
               enddo
            endif
         enddo
c
         do k=1, ksr1-2, 2
            l=l+1
            do i=1, ims, 2
               n(l)=n(l)+1
               write(io,*) -npos(nods(i,j,k),2), npos(nods(i,j,k),3)
            enddo
         enddo
         k = ksr1
         l=l+1
         n(l)=n(l)+1
         write(io,*) -npos(nods(1,j,k),2), npos(nods(1,j,k),3)
         do i=5, ims-4, 8
            n(l)=n(l)+4
            write(io,*) -npos(nods(i-2,j,k),2),npos(nods(i-2,j,k),3)
            write(io,*) -npos(nods(i,j,k+2),2), npos(nods(i,j,k+2),3)
            write(io,*) -npos(nods(i+2,j,k),2), npos(nods(i+2,j,k),3)
            write(io,*) -npos(nods(i+4,j,k),2), npos(nods(i+4,j,k),3)
         enddo
         do k=ksr2, kms, 2
            l=l+1
            do i=1, ims, 4
               n(l)=n(l)+1
               write(io,*) -npos(nods(i,j,k),2), npos(nods(i,j,k),3)
            enddo
         enddo
c
      elseif (sfred_type.eq.3) then
c
         do i=1, ims, 2
            l=l+1
            if (mod(i,6).eq.1) then
               do k=1, kms, 2
                  n(l)=n(l)+1
                  write(io,*) -npos(nods(i,j,k),2), npos(nods(i,j,k),3)
               enddo
            else
               do k=1, ksr1, 2
                  n(l)=n(l)+1
                  write(io,*) -npos(nods(i,j,k),2), npos(nods(i,j,k),3)
               enddo
            endif
         enddo
c
         do k=1, ksr1-2, 2
            l=l+1
            do i=1, ims, 2
               n(l)=n(l)+1
               write(io,*) -npos(nods(i,j,k),2), npos(nods(i,j,k),3)
            enddo
         enddo
c
         k = ksr1
         l=l+1
         n(l)=n(l)+1
         write(io,*) -npos(nods(1,j,k+2),2), npos(nods(1,j,k+2),3)
         do i=4, ims, 6
            n(l)=n(l)+3
            write(io,*) -npos(nods(i-1,j,k),2), npos(nods(i-1,j,k),3)
            write(io,*) -npos(nods(i+1,j,k),2), npos(nods(i+1,j,k),3)
            write(io,*) -npos(nods(i+3,j,k+2),2),npos(nods(i+3,j,k+2),3)
         enddo
         do k=ksr2, kms, 2
            l=l+1
            do i=1, ims, 3
               n(l)=n(l)+1
               write(io,*) -npos(nods(i,j,k),2), npos(nods(i,j,k),3)
            enddo
         enddo
      endif
c
c...Zone A:
c
C . Z-direction
	do k=1, 2*mv+1, 2
	   l=l+1
	   do i=1, ia1, 2
	      if (noda(i,1,k).gt.0) then
	         write(io,*) -npos(noda(i,1,k),2), npos(noda(i,1,k),3)
	         n(l)=n(l)+1
	      endif
	   enddo
	   l=l+1
	   do i=ia2, ima, 2
	      if (noda(i,1,k).gt.0) then
	         write(io,*) -npos(noda(i,1,k),2), npos(noda(i,1,k),3)
	         n(l)=n(l)+1
	      endif
	   enddo
	enddo
	do k=2*mv+3, kar1-2, 2
	   l=l+1
	   do i=1, ima, 2
 	     if (noda(i,1,k).gt.0) then
	       write(io,*) -npos(noda(i,1,k),2), npos(noda(i,1,k),3)
	       n(l)=n(l)+1
	     endif
	   enddo
	enddo
	if (rtype.eq.0) then
	   do k=kar1, kma, 2
	      l=l+1
	      do i=1, ima, 2
 	         if (noda(i,1,k).gt.0) then
	           write(io,*) -npos(noda(i,1,k),2),npos(noda(i,1,k),3)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	else
	   l=l+1
	   if (mod(ma,4).eq.0) then
	      write(io,*) -npos(noda(1,1,kar1),2),npos(noda(1,1,kar1),3)
	      write(io,*) -npos(noda(5,1,kar1),2),npos(noda(5,1,kar1),3)
	      n(l)=n(l)+2
	      do i=9, ima, 8
	         write(io,*) -npos(noda(i-2,1,kar1),2),
     &                        npos(noda(i-2,1,kar1),3)
	         write(io,*) -npos(noda(i,1,kar1+2),2),
     &                        npos(noda(i,1,kar1+2),3)
	         n(l)=n(l)+2
                 if ((i+4).lt.na1) then
	           if (noda(i+2,1,kar1).gt.0) then
	             write(io,*) -npos(noda(i+2,1,kar1),2),
     &                            npos(noda(i+2,1,kar1),3)
	             write(io,*) -npos(noda(i+4,1,kar1),2),
     &                            npos(noda(i+4,1,kar1),3)
	             n(l)=n(l)+2
	           endif
	         endif
	      enddo
	   else
	      write(io,*) -npos(noda(5,1,kar1+2),2),
     &                     npos(noda(5,1,kar1+2),3)
	      n(l)=n(l)+1
	      do i=9, ima, 8
	         write(io,*) -npos(noda(i-2,1,kar1),2),
     &                        npos(noda(i-2,1,kar1),3)
	         write(io,*) -npos(noda(i,1,kar1),2),
     &                        npos(noda(i,1,kar1),3)
	         write(io,*) -npos(noda(i+2,1,kar1),2),
     &                        npos(noda(i+2,1,kar1),3)
	         write(io,*) -npos(noda(i+4,1,kar1+2),2),
     &                        npos(noda(i+4,1,kar1+2),3)
	         n(l)=n(l)+4
	      enddo
	   endif
	   do k=kar1+2, kma, 2
	      l=l+1
	      if((rtype.eq.2).and.(mod(na,4).ne.0).and.(k.eq.kar2))then
                i1=3
	      elseif ( (rtype.eq.2).and.(k.eq.kar2) ) then
	        i1=ima
	        l=l-1
	      else
	        i1=1
	      endif
	      do i=i1, ima, 2
 	         if (noda(i,1,k).gt.0) then
	           write(io,*) -npos(noda(i,1,k),2),npos(noda(i,1,k),3)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	endif
 
C . Y-direction
	l=l+1
	do k=1, kma, 2
	   if (noda(1,1,k).gt.0) then
              write(io,*)  -npos(noda(1,1,k),2), npos(noda(1,1,k),3)
	      n(l)=n(l)+1
	   endif
	enddo
	if (rtype.lt.2) then
	   if (mod(na,4).ne.0) then
	      l=l+1
	      do k=1, kma, 2
	         if (noda(3,1,k).gt.0) then
                    write(io,*)  -npos(noda(3,1,k),2),
     &                            npos(noda(3,1,k),3)
	            n(l)=n(l)+1
	         endif
	      enddo
	   endif
	else
	   if (mod(na,4).ne.0) then
	      l=l+1
	      do k=1, kar2, 2
	         if (noda(3,1,k).gt.0) then
                    write(io,*)  -npos(noda(3,1,k),2),
     &                            npos(noda(3,1,k),3)
	            n(l)=n(l)+1
 	         endif
	      enddo
	      write(io,*) -npos(noda(1,1,kar2+2),2),
     &                     npos(noda(1,1,kar2+2),3)
	      n(l)=n(l)+1
	   endif
	endif
	if (rtype.eq.0) then
	   do i=5, ima, 2
	      l=l+1
 	      do k=1, kma, 2
	         if (noda(i,1,k).gt.0) then
	           write(io,*) -npos(noda(i,1,k),2),npos(noda(i,1,k),3)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	else
	   do i=5, ima, 2
              if (mod(i,4).eq.1) then
	         kend=kma
	      else
	         kend=kar1
	      endif
	      l=l+1
 	      do k=1, kend, 2
	         if (noda(i,1,k).gt.0) then
	           write(io,*) -npos(noda(i,1,k),2),npos(noda(i,1,k),3)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	endif
	close(io)
c
c...Generate the information file
c
	z2=t
	z1=0.
	do k=1, kma, 2
          if (npos(noda(ima,1,k),3).gt.z2) z2=npos(noda(ima,1,k),3)
          if (npos(noda(1,1,k),3).lt.z1) z1=npos(noda(1,1,k),3)
	enddo
	open(unit=32,file=filinf,status='unknown')
 	write(32,'(t1,a)') '*CURVE1'
	if ((z2-z1).gt.npos(noda(1,1,kma),2) ) then
	   sf=(z2-z1)/20.
	else
	   sf=npos(noda(1,1,kma),2)/20.
	endif
	write(32,'(t1,2f5.1,tr2,6(tr1,g10.4))')
     &  3.0,3.0,-21.*sf,sf,-sf+z1,21.*sf+z1,18.0 , 18.0
	write(32,'(t1,2(g12.6,tr2),i1,a,g10.4,tr2,i1,tr2,f4.2)')
     &  -1.*npos(noda(1,1,kma),2), npos(noda(1,1,kma),2)/2., 1,
     &  ' 0.0 ', (z2-z1), 1, 0.35
	write(32,'(t1,i3)') l
	do i=1, l
	   write(32,'(t1,i4,a)') n(i),' 3 0 401 0 1'
	enddo
	write(32,'(t1,a)') '*FINISH'
	close(32)
c
c  Zoom into the crack front region
	filinf = job(1:jobh)//'_2z.inf'
	open(unit=32,file=filinf,status='unknown')
        write(32,'(t1,a)') '*CURVE1'
        z1 = npos(nods(ims,1,kms),3)
        z2 = npos(nods(1,1,kms),3)
        da = npos(nodc(imc-2,1,1),3)-npos(nodc(5,1,1),3)
        sf = (z2-z1)/20.0
	write(32,'(t1,2f5.1,tr2,6(tr1,g10.4))')
     &    3.0,2.0, -21.*sf,sf, z1-sf, z1+21.*sf, 18.0,18.0
	write(32,'(t1,a,g14.4,tr2,i1,2g14.4,tr2,i1,tr2,f4.2)')
     &   '0.0 ', (z2-z1+2.*sf)/4, 1, a, da/2., 1, 0.35
	write(32,'(t1,i3)') l
	do i=1, l
	   write(32,'(t1,i4,a)') n(i),' 3 0 401 0 1'
	enddo
	write(32,'(t1,a)') '*FINISH'
	close(32)
c
c---------------------------------
c..... Plot file 3 ( the rear plane, Y=Ymax )
c---------------------------------
c
	do i=1, nm
	   n(i)=0
	enddo
	fildat = job(1:jobh)//'_3.dat'
	filinf = job(1:jobh)//'_3.inf'
        io = 31
	open(unit=io,file=fildat,status='unknown')
C . . . a) RAMEN !
	z2=npos(noda(ima,1,kma),3)
	z1=npos(noda(1,1,kma),3)
        l=1
	n(1)=5
	write(io,*) 0.0, z1
	write(io,*) 0.0, z2
	write(io,*) w  , z2
	write(io,*) w  , z1
	write(io,*) 0.0, z1
C . . . b) RHO-LED - ZON A !!
	if (rtype.le.1) then
	   dj=2
	else
	   dj=4
	endif
	do j=1+dj, jma, dj
	   l=l+1
	   if ( (rtype.le.1).and.(mod(j,4).eq.3) ) then
	      i1=3
	   else
	      i1=1
	   endif
	   do i=i1, ima
	      if (noda(i,j,kma).gt.0) then
	        n(l)=n(l)+1
	        write(io,*) npos(noda(i,j,kma),1),npos(noda(i,j,kma),3)
	      endif
	   enddo
	enddo
 
C . . . c) FI-LED - ZON A !!
	if (rtype.le.1) then
	   l=l+1
  	   if (mod(na,4).ne.0) then
	      left=.true.
	      write(io,*) npos(noda(3,1,kma),1), npos(noda(3,1,kma),3)
	   else
	      left=.false.
	      write(io,*) npos(noda(1,1,kma),1), npos(noda(1,1,kma),3)
	   endif
           n(l)=n(l)+1
	   do j=3, jma-4, 4
	      if (left) then
	       write(io,*) npos(noda(3,j,kma),1), npos(noda(3,j,kma),3)
	   write(io,*) npos(noda(1,j+2,kma),1), npos(noda(1,j+2,kma),3)
                n(l)=n(l)+2
	        left=.false.
	      else
	        write(io,*) npos(noda(3,j,kma),1), npos(noda(3,j,kma),3)
	   write(io,*) npos(noda(3,j+2,kma),1), npos(noda(3,j+2,kma),3)
                n(l)=n(l)+2
	        left=.true.
	      endif
	   enddo
	endif
	if (rtype.eq.0) then
	   di=2
	else
	   di=4
	endif
	do i=5, ima, di
	   l=l+1
	   do j=1, jma
	      if (noda(i,j,kma).gt.0) then
	         write(io,*)npos(noda(i,j,kma),1),npos(noda(i,j,kma),3)
	         n(l)=n(l)+1
	      endif
	   enddo
	enddo
 
C . . . d) ZON B !!
C . . . Z-LED . . . . . .
	do i=3, imb-2, 2
	   l=l+1
	   do j=1, jmb
	      if (nodb(i,j,kma).gt.0) then
	        write(io,*) npos(nodb(i,j,kma),1),npos(nodb(i,j,kma),3)
	        n(l)=n(l)+1
	      endif
	   enddo
	enddo
C . . . X-LED . . . . .
	do j=dj+1, jmb-dj, dj
	   l=l+1
	   do i=1, imb
	      write(io,*) npos(nodb(i,j,kma),1),npos(nodb(i,j,kma),3)
	      n(l)=n(l)+1
	   enddo
	enddo
	close(io)
c
c...Generate the information file
c
	open(unit=32,file=filinf,status='unknown')
        write(32,'(t1,a)') '*CURVE1'
	if ((z2-z1).gt.w) then
	   sf=(z2-z1)/20.
	else
	   sf=w/20.
	endif
	write(32,'(t1,2f5.1,tr2,6(tr1,g10.4))')
     &    3.0,3.0,-sf,21.*sf,-sf+z1,21.*sf+z1,18.,18.
	write(32,'(t1,2(a,g10.4,tr2,i1),tr2,f4.2)')
     &  '0.0 ', w/2., 1, ' 0.0 ', (z2-z1), 1, 0.35
	write(32,'(t1,i3)') l
	do i=1, l
	   write(32,'(t1,i4,a)') n(i),' 3 0 401 0 1'
	enddo
	write(32,'(t1,a)') '*FINISH'
	close(32)
 
c---------------------------------
c..... Plot file 4 ( the top side )
c---------------------------------
	do i=1, nm
	   n(i)=0
	enddo
	l=0
	fildat = job(1:jobh)//'_4.dat'
	filinf = job(1:jobh)//'_4.inf'
        io = io
	open(unit=io,file=fildat,status='unknown')
C Zon A :
C . .  Y-direction
       	if (rtype.le.1) then
	   do j=1, jma-jmb+1, 2
	      l=l+1
	      do k=1, kma
	         if (noda(ima,j,k).gt.0) then
	            write(io,*) npos(noda(ima,j,k),1),
     &                          npos(noda(ima,j,k),2)
	            n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	else
	   do j=1, jma-jmb+1, 2
	      l=l+1
	      if (mod(j,4).eq.3) then
	         kend=kar2
	      else
	         kend=kma
	      endif
	      do k=1, kend
	         if (noda(ima,j,k).gt.0) then
	            write(io,*) npos(noda(ima,j,k),1),
     &                          npos(noda(ima,j,k),2)
	            n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	endif
C . .  X-direction
       	if (rtype.le.1) then
	   do k=1, kma, 2
	      l=l+1
	      do j=1, jma-jmb+1
	         if (noda(ima,j,k).gt.0) then
	            write(io,*) npos(noda(ima,j,k),1),
     &                          npos(noda(ima,j,k),2)
	            n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	else
	   do k=1, kar2-2, 2
	      l=l+1
	      do j=1, jma-jmb+1
	         if (noda(ima,j,k).gt.0) then
	            write(io,*) npos(noda(ima,j,k),1),
     &                          npos(noda(ima,j,k),2)
	            n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	   if (mod(na,4).ne.0) then
	      l=l+1
              write(io,*) npos(noda(ima,1,kar2),1),
     &                    npos(noda(ima,1,kar2),2)
              write(io,*) npos(noda(ima,3,kar2),1),
     &                    npos(noda(ima,3,kar2),2)
              write(io,*) npos(noda(ima,5,kar2+2),1),
     &                    npos(noda(ima,5,kar2+2),2)
	      n(l)=n(l)+3
	      do j=9, jma-jmb+1, 8
	         write(io,*) npos(noda(ima,j-2,kar2),1),
     &                       npos(noda(ima,j-2,kar2),2)
	         write(io,*) npos(noda(ima,j,kar2),1),
     &                       npos(noda(ima,j,kar2),2)
	         n(l)=n(l)+2
                 if ((jma-jmb+1).gt.9) then
	            write(io,*) npos(noda(ima,j+2,kar2),1),
     &                          npos(noda(ima,j+2,kar2),2)
	            write(io,*) npos(noda(ima,j+4,kar2+2),1),
     &                          npos(noda(ima,j+4,kar2+2),2)
	            n(l)=n(l)+2
                 endif
	      enddo
	   else
	      l=l+1
              write(io,*) npos(noda(ima,1,kar2+2),1),
     &                    npos(noda(ima,1,kar2+2),2)
	      n(l)=n(l)+1
	      do j=5, jma-jmb+1, 8
	         write(io,*) npos(noda(ima,j-2,kar2),1),
     &                       npos(noda(ima,j-2,kar2),2)
	         write(io,*) npos(noda(ima,j,kar2),1),
     &                       npos(noda(ima,j,kar2),2)
	         write(io,*) npos(noda(ima,j+2,kar2),1),
     &                       npos(noda(ima,j+2,kar2),2)
	         write(io,*) npos(noda(ima,j+4,kar2+2),1),
     &                       npos(noda(ima,j+4,kar2+2),2)
	         n(l)=n(l)+4
	      enddo
	   endif
	   do k=kar2+2, kma, 2
	      l=l+1
	      do j=1, jma-jmb+1
	         if (noda(ima,j,k).gt.0) then
	            write(io,*) npos(noda(ima,j,k),1),
     &                          npos(noda(ima,j,k),2)
	            n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	endif
C Zon B :
C . .  Y-direction
        do i=3, imb, 2
	   l=l+1
	   do k=1, kma
	      if (nodb(i,1,k).gt.0) then
	         write(io,*) npos(nodb(i,1,k),1), npos(nodb(i,1,k),2)
	         n(l)=n(l)+1
	      endif
	   enddo
	enddo
C . .  X-direction
        do k=1, kma, 2
	   if (.not.( (rtype.eq.2).and.(mod(nb,4).eq.0).and.
     &                (k.eq.kar2) ) ) then
	      l=l+1
	      do i=1, imb
	         if (nodb(i,1,k).gt.0) then
	            write(io,*) npos(nodb(i,1,k),1), npos(nodb(i,1,k),2)
	            n(l)=n(l)+1
	         endif
	      enddo
	   endif
	enddo
	close(io)
c
c...Generate the information file
c
	open(unit=32,file=filinf,status='unknown')
        write(32,'(t1,a)') '*CURVE1'
	if (npos(noda(ima,1,kma),2).gt.w) then
	   sf=npos(noda(ima,1,kma),2)/20.
	else
	   sf=w/20.
	endif
	write(32,'(t1,2f5.1,tr2,6(tr1,g10.4))')
     &   3.0,3.0,-sf,21.*sf,-sf,21.*sf,18.,18.
	write(32,'(t1,2(a,g10.4,tr2,i1),tr2,f4.2)')
     &   '0.0 ', w/4., 1, ' 0.0 ', npos(noda(ima,1,kma),2)/4., 1, 0.35
	write(32,'(t1,i3)') l
	do i=1, l
	   write(32,'(t1,i4,a)') n(i),' 3 0 401 0 1'
	enddo
	write(32,'(t1,a)') '*FINISH'
	close(32)
 
c---------------------------------
c.... Plot file 5 ( the side/surface containing the surface, Z=0 )
c---------------------------------
	do i=1, nm
	   n(i)=0
	enddo
	l=0
c
	fildat = job(1:jobh)//'_5.dat'
	filinf = job(1:jobh)//'_5.inf'
        io=35
	open(unit=io,file=fildat,status='unknown')
c
c...Zone C:
c
        j = jmc
        l=l+1
        do i=5, imc, 2
	   n(l)=n(l)+1
           write(io,*) npos(nodc(i,j,1),1), npos(nodc(i,j,1),2)
        enddo
        do k=3, kmc-4, 2
           l=l+1
           do i=1, imc, 2
              n(l)=n(l)+1
              write(io,*) npos(nodc(i,j,k),1), npos(nodc(i,j,k),2)
           enddo
        enddo
        k = kmc-2
        left = .true.
        l=l+1
        do i=3, imc-3, 4
         if (left) then
            n(l)=n(l)+3
            write(io,*) npos(nodc(i-2,j,k),1), npos(nodc(i-2,j,k),2)
            write(io,*) npos(nodc(i,j,k),1), npos(nodc(i,j,k),2)
            write(io,*) npos(nodc(i+2,j,k+2),1),npos(nodc(i+2,j,k+2),2)
            left = .false.
         else
            n(l)=n(l)+1
            write(io,*) npos(nodc(i,j,k),1),npos(nodc(i,j,k),2)
            left = .true.
         endif
      enddo
      n(l)=n(l)+1
      write(io,*) npos(nodc(imc,j,kmc),1),npos(nodc(imc,j,kmc),2)
c
      do i=3, imc-2, 2
         l=l+1
         if (mod(i,4).eq.1) then
            do k=1, kmc, 2
               if (nodc(i,j,k).gt.0) then
                  n(l)=n(l)+1
                  write(io,*) npos(nodc(i,j,k),1), npos(nodc(i,j,k),2)
               endif
            enddo
         else
            do k=1, kmc-2, 2
               if (nodc(i,j,k).gt.0) then
                  n(l)=n(l)+1
                  write(io,*) npos(nodc(i,j,k),1), npos(nodc(i,j,k),2)
               endif
            enddo
         endif
      enddo
c
c...Zone S:
c
      ksr2 = ksr1 + 2
      j = jms
c
      if (sfred_type.eq.1) then
c
         do i=1, ims, 2
            l=l+1
            do k=1, kms, 2
               n(l)=n(l)+1
               write(io,*) npos(nods(i,j,k),1), npos(nods(i,j,k),2)
            enddo
         enddo
c
         do k=1,  kms, 2
            l=l+1
            do i=1, ims, 2
               n(l)=n(l)+1
               write(io,*) npos(nods(i,j,k),1), npos(nods(i,j,k),2)
            enddo
         enddo
c
      elseif (sfred_type.eq.2) then
c
         do i=1, ims, 2
            l=l+1
            if (mod(i,4).eq.1) then
               do k=1, kms, 2
                  n(l)=n(l)+1
                  write(io,*) npos(nods(i,j,k),1), npos(nods(i,j,k),2)
               enddo
            else
               do k=1, ksr1, 2
                  n(l)=n(l)+1
                  write(io,*) npos(nods(i,j,k),1), npos(nods(i,j,k),2)
               enddo
            endif
         enddo
c
         do k=1, ksr1-2, 2
            l=l+1
            do i=1, ims, 2
               n(l)=n(l)+1
               write(io,*) npos(nods(i,j,k),1), npos(nods(i,j,k),2)
            enddo
         enddo
         k = ksr1
         l=l+1
         n(l)=n(l)+1
         write(io,*)  npos(nods(1,j,k),1), npos(nods(1,j,k),2)
         do i=5, ims-4, 8
            n(l)=n(l)+4
            write(io,*) npos(nods(i-2,j,k),1), npos(nods(i-2,j,k),2)
            write(io,*) npos(nods(i,j,k+2),1), npos(nods(i,j,k+2),2)
            write(io,*) npos(nods(i+2,j,k),1), npos(nods(i+2,j,k),2)
            write(io,*) npos(nods(i+4,j,k),1), npos(nods(i+4,j,k),2)
         enddo
         do k=ksr2, kms, 2
            l=l+1
            do i=1, ims, 4
               n(l)=n(l)+1
               write(io,*) npos(nods(i,j,k),1), npos(nods(i,j,k),2)
            enddo
         enddo
c
      elseif (sfred_type.eq.3) then
c
         do i=1, ims, 2
            l=l+1
            if (mod(i,6).eq.1) then
               do k=1, kms, 2
                  n(l)=n(l)+1
                  write(io,*)  npos(nods(i,j,k),1), npos(nods(i,j,k),2)
               enddo
            else
               do k=1, ksr1, 2
                  n(l)=n(l)+1
                  write(io,*) npos(nods(i,j,k),1), npos(nods(i,j,k),2)
               enddo
            endif
         enddo
c
         do k=1, ksr1-2, 2
            l=l+1
            do i=1, ims, 2
               n(l)=n(l)+1
               write(io,*)  npos(nods(i,j,k),1), npos(nods(i,j,k),2)
            enddo
         enddo
c
         k = ksr1
         l=l+1
         n(l)=n(l)+1
         write(io,*)  npos(nods(1,j,k+2),1), npos(nods(1,j,k+2),2)
         do i=4, ims, 6
            n(l)=n(l)+3
            write(io,*)  npos(nods(i-1,j,k),2),  npos(nods(i-1,j,k),2)
            write(io,*)  npos(nods(i+1,j,k),2),  npos(nods(i+1,j,k),2)
            write(io,*)  npos(nods(i+3,j,k+2),2),npos(nods(i+3,j,k+2),2)
         enddo
         do k=ksr2, kms, 2
            l=l+1
            do i=1, ims, 3
               n(l)=n(l)+1
               write(io,*)  npos(nods(i,j,k),1), npos(nods(i,j,k),2)
            enddo
         enddo
      endif
c
c...Zone A:
c . .  X-direction, I=1
	do k=1, kma, 2
	   if (.not.( (rtype.ge.1).and.(mod(ma,4).ne.0).and.
     &           (k.eq.kar1) ) ) then
	      l=l+1
	      do j=1, jma-4, 2
	         if (noda(1,j,k).gt.0) then
	            write(io,*) npos(noda(1,j,k),1),npos(noda(1,j,k),2)
	            n(l)=n(l)+1
	         endif
	      enddo
	   endif
	enddo
C . .  Y-direction, I=1
	do j=1, jma-4, 4
	   l=l+1
	   do k=1, kma, 2
	      if (noda(1,j,k).gt.0) then
	         write(io,*) npos(noda(1,j,k),1),npos(noda(1,j,k),2)
	         n(l)=n(l)+1
	      endif
	   enddo
	enddo
c
c . . X-direction, I=IMA
	do k=1, 2*mv+1, 2
	   if (ia1.gt.5) then
	      l=l+1
	      do i=5, ia1, 2
	         if (noda(i,jma,k).gt.0) then
	            write(io,*) npos(noda(i,jma,k),1),
     &                           npos(noda(i,jma,k),2)
	            n(l)=n(l)+1
	         endif
	      enddo
	   endif
	   l=l+1
	   do i=ia2, ima, 2
	      if (noda(i,jma,k).gt.0) then
	         write(io,*) npos(noda(i,jma,k),1),npos(noda(i,jma,k),2)
	         n(l)=n(l)+1
	      endif
	   enddo
	enddo
	do k=2*mv+3, kar1-2, 2
	   l=l+1
	   do i=5, ima, 2
 	     if (noda(i,jma,k).gt.0) then
	       write(io,*) npos(noda(i,jma,k),1), npos(noda(i,jma,k),2)
	       n(l)=n(l)+1
	     endif
	   enddo
	enddo
	if (rtype.eq.0) then
	   do k=kar1, kma, 2
	      l=l+1
	      do i=5, ima, 2
 	         if (noda(i,jma,k).gt.0) then
	        write(io,*) npos(noda(i,jma,k),1),npos(noda(i,jma,k),2)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	else
	   l=l+1
	   if (mod(ma,4).eq.0) then
	      write(io,*) npos(noda(5,jma,kar1),1),
     &                    npos(noda(5,jma,kar1),2)
	      n(l)=n(l)+1
	      do i=9, ima, 8
	         write(io,*) npos(noda(i-2,jma,kar1),1),
     &                        npos(noda(i-2,jma,kar1),2)
	         write(io,*) npos(noda(i,jma,kar1+2),1),
     &                        npos(noda(i,jma,kar1+2),2)
	         n(l)=n(l)+2
                 if (i.lt.ima) then
	            write(io,*) npos(noda(i+2,jma,kar1),1),
     &                           npos(noda(i+2,jma,kar1),2)
	            write(io,*) npos(noda(i+4,jma,kar1),1),
     &                           npos(noda(i+4,jma,kar1),2)
	            n(l)=n(l)+2
	         endif
	      enddo
	   else
	      write(io,*) npos(noda(5,jma,kar1+2),1),
     &                    npos(noda(5,jma,kar1+2),2)
	      n(l)=n(l)+1
	      do i=9, ima, 8
	         write(io,*) npos(noda(i-2,jma,kar1),1),
     &                        npos(noda(i-2,jma,kar1),2)
	         write(io,*) npos(noda(i,jma,kar1),1),
     &                        npos(noda(i,jma,kar1),2)
	         write(io,*) npos(noda(i+2,jma,kar1),1),
     &                        npos(noda(i+2,jma,kar1),2)
	         write(io,*) npos(noda(i+4,jma,kar1+2),1),
     &                        npos(noda(i+4,jma,kar1+2),2)
	         n(l)=n(l)+4
	      enddo
	   endif
	   do k=kar1+2, kma, 2
	      l=l+1
	      do i=5, ima, 2
 	         if (noda(i,jma,k).gt.0) then
	        write(io,*) npos(noda(i,jma,k),1),npos(noda(i,jma,k),2)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	endif
C . .  Y-direction, I=1
	if (rtype.eq.0) then
	   do i=7, ima, 2
	      l=l+1
 	      do k=1, kma, 2
	         if (noda(i,jma,k).gt.0) then
	        write(io,*) npos(noda(i,jma,k),1),npos(noda(i,jma,k),2)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	else
	   do i=7, ima, 2
              if (mod(i,4).eq.1) then
	         kend=kma
	      else
	         kend=kar1
	      endif
	      l=l+1
 	      do k=1, kend, 2
	         if (noda(i,jma,k).gt.0) then
	       write(io,*) npos(noda(i,jma,k),1),npos(noda(i,jma,k),2)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	endif
 
C...ZON B
C . .  X-direction
	do k=1, kma, 2
	   if (.not.((rtype.ge.1).and.(k.eq.kar1)) ) then
	      l=l+1
	      do i=1, imb, 2
	         if (nodb(i,jmb,k).gt.0) then
	        write(io,*) npos(nodb(i,jmb,k),1),npos(nodb(i,jmb,k),2)
	            n(l)=n(l)+1
	         endif
	      enddo
	   endif
	enddo
C . .  Y-direction, I=1
	do i=3, imb, 2
	   l=l+1
	   do k=1, kma, 2
	      if (nodb(i,jmb,k).gt.0) then
	         write(io,*) npos(nodb(i,jmb,k),1),npos(nodb(i,jmb,k),2)
	         n(l)=n(l)+1
	      endif
	   enddo
	enddo
	close(io)
c
c...Generate the information file
c
	open(unit=32,file=filinf,status='unknown')
        write(32,'(t1,a)') '*CURVE1'
	if (npos(noda(1,1,kma),2).gt.w) then
	   sf=npos(noda(1,1,kma),2)/20.
	else
	   sf=w/20.
	endif
	write(32,'(t1,2f5.1,tr2,6(tr1,g10.4))')
     &     3.0,3.0,-sf,21.*sf,-sf,21.*sf,18.,18.
	write(32,'(t1,2(a,g10.4,tr2,i1),tr2,f4.2)')
     &    '0.0 ', w/4., 1, ' 0.0 ', npos(noda(1,1,kma),2)/4., 1, 0.35
	write(32,'(t1,i3)') l
	do i=1, l
	   write(32,'(t1,i4,a)') n(i),' 3 0 401 0 1'
	enddo
	write(32,'(t1,a)') '*FINISH'
	close(32)
c
c  Zoom into the crack front region
	filinf = job(1:jobh)//'_5z.inf'
	open(unit=32,file=filinf,status='unknown')
        write(32,'(t1,a)') '*CURVE1'
        x1 = npos(nods(ims,jms,kms),1)
        x2 = npos(nods(1,jms,kms),1)
        da = npos(nodc(imc-2,jmc,1),1)-npos(nodc(5,jmc,1),1)
        sf = (x2-x1)/20.0
	write(32,'(t1,2f5.1,tr2,6(tr1,g10.4))')
     &    3.0,3.0, x1-sf, x1+21.*sf, -sf,21.*sf,  18.0,18.0
	write(32,'(t1,2g14.4,tr2,i1,a,g14.4,tr2,i1,tr2,f4.2)')
     &     c, da/2., 1, ' 0.0 ', da/4, 1, 0.35
	write(32,'(t1,i3)') l
	do i=1, l
	   write(32,'(t1,i4,a)') n(i),' 3 0 401 0 1'
	enddo
	write(32,'(t1,a)') '*FINISH'
	close(32)
c
	return
	end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      subroutine plot_zone_c(nodc,nc1,nc2,nc3)
c
      implicit none
      include 'common_nod.f'
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3)
      integer  i,j,k,jj, io
      double precision x,y
      logical left
c
      integer      imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
      io = 25
      open(unit=io,file='cell_1.plt',status='unknown')
      k = 1
      do j=1, jmc, 4
         write(io,'(t1,a)') 'ZONE'
         do i=5, imc, 2
            x = npos(nodc(i,j,k),1)
            y = npos(nodc(i,j,k),3)
            write(io,'(2(g16.8))') x,y
         enddo
         jj = j + 2
         if (jj.le.jmc) then
            do i=imc, 5, -2
               x = npos(nodc(i,jj,k),1)
               y = npos(nodc(i,jj,k),3)
               write(io,'(2(g16.8))') x,y
            enddo
         endif
      enddo
c
      k = 1
      do i=5, imc, 2
         write(io,'(t1,a)') 'ZONE'
         do j=1, jmc, 2
            x = npos(nodc(i,j,k),1)
            y = npos(nodc(i,j,k),3)
            write(io,'(2(g16.8))') x,y
         enddo
      enddo
c
      close(io)
c
      io = 26
      open(unit=io,file='cell_2a.plt',status='unknown')
      j = 1
c
      k = 1
      write(io,'(t1,a)') 'ZONE'
      do i=5, imc, 2
         x = npos(nodc(i,j,k),3)
         y = npos(nodc(i,j,k),2)
         write(io,'(2(g16.8))') x,y
      enddo
c
      do k=3, 5, 2
         write(io,'(t1,a)') 'ZONE'
         do i=1, imc, 2
            x = npos(nodc(i,j,k),3)
            y = npos(nodc(i,j,k),2)
            write(io,'(2(g16.8))') x,y
         enddo
      enddo
c
      k = 7
      left = .true.
      write(io,'(t1,a)') 'ZONE'
      do i=3, imc-3, 4
         if (left) then
            x = npos(nodc(i-2,j,k),3)
            y = npos(nodc(i-2,j,k),2)
            write(io,'(2(g16.8))') x,y
            x = npos(nodc(i,j,k),3)
            y = npos(nodc(i,j,k),2)
            write(io,'(2(g16.8))') x,y
            x = npos(nodc(i+2,j,k+2),3)
            y = npos(nodc(i+2,j,k+2),2)
            write(io,'(2(g16.8))') x,y
            left = .false.
         else
            x = npos(nodc(i,j,k),3)
            y = npos(nodc(i,j,k),2)
            write(io,'(2(g16.8))') x,y
            left = .true.
         endif
      enddo
      x = npos(nodc(imc,j,kmc),3)
      y = npos(nodc(imc,j,kmc),2)
      write(io,'(2(g16.8))') x,y
c
      k = 9
      write(io,'(t1,a)') 'ZONE'
      do i=1, imc, 2
         x = npos(nodc(i,j,k),3)
         y = npos(nodc(i,j,k),2)
         write(io,'(2(g16.8))') x,y
      enddo
c
      do i=1, imc-2, 2
         write(io,'(t1,a)') 'ZONE'
         if (mod(i,4).eq.1) then
            do k=1, kmc, 2
               if (nodc(i,j,k).gt.0) then
                  x = npos(nodc(i,j,k),3)
                  y = npos(nodc(i,j,k),2)
                  write(io,'(2(g16.8))') x,y
               endif
            enddo
         else
            do k=1, kmc-2, 2
               if (nodc(i,j,k).gt.0) then
                  x = npos(nodc(i,j,k),3)
                  y = npos(nodc(i,j,k),2)
                  write(io,'(2(g16.8))') x,y
               endif
            enddo
         endif
      enddo
      write(io,'(t1,a)') 'ZONE'
      do k=1, kmc, 2
         x = npos(nodc(imc,j,k),3)
         y = npos(nodc(imc,j,k),2)
         write(io,'(2(g16.8))') x,y
      enddo
c
      close(io)
c
      io = 26
      open(unit=io,file='cell_2b.plt',status='unknown')
      j = jmc
c
      k = 1
      write(io,'(t1,a)') 'ZONE'
      do i=5, imc, 2
         x = npos(nodc(i,j,k),1)
         y = npos(nodc(i,j,k),2)
         write(io,'(2(g16.8))') x,y
      enddo
c
      do k=3, 5, 2
         write(io,'(t1,a)') 'ZONE'
         do i=1, imc, 2
            x = npos(nodc(i,j,k),1)
            y = npos(nodc(i,j,k),2)
            write(io,'(2(g16.8))') x,y
         enddo
      enddo
c
      k = 7
      left = .true.
      write(io,'(t1,a)') 'ZONE'
      do i=3, imc-3, 4
         if (left) then
            x = npos(nodc(i-2,j,k),1)
            y = npos(nodc(i-2,j,k),2)
            write(io,'(2(g16.8))') x,y
            x = npos(nodc(i,j,k),1)
            y = npos(nodc(i,j,k),2)
            write(io,'(2(g16.8))') x,y
            x = npos(nodc(i+2,j,k+2),1)
            y = npos(nodc(i+2,j,k+2),2)
            write(io,'(2(g16.8))') x,y
            left = .false.
         else
            x = npos(nodc(i,j,k),1)
            y = npos(nodc(i,j,k),2)
            write(io,'(2(g16.8))') x,y
            left = .true.
         endif
      enddo
      x = npos(nodc(imc,j,kmc),1)
      y = npos(nodc(imc,j,kmc),2)
      write(io,'(2(g16.8))') x,y
c
      k = 9
      write(io,'(t1,a)') 'ZONE'
      do i=1, imc, 2
         x = npos(nodc(i,j,k),1)
         y = npos(nodc(i,j,k),2)
         write(io,'(2(g16.8))') x,y
      enddo
c
      do i=1, imc-2, 2
         write(io,'(t1,a)') 'ZONE'
         if (mod(i,4).eq.1) then
            do k=1, kmc, 2
               if (nodc(i,j,k).gt.0) then
                  x = npos(nodc(i,j,k),1)
                  y = npos(nodc(i,j,k),2)
                  write(io,'(2(g16.8))') x,y
               endif
            enddo
         else
            do k=1, kmc-2, 2
               if (nodc(i,j,k).gt.0) then
                  x = npos(nodc(i,j,k),1)
                  y = npos(nodc(i,j,k),2)
                  write(io,'(2(g16.8))') x,y
               endif
            enddo
         endif
      enddo
      write(io,'(t1,a)') 'ZONE'
      do k=1, kmc, 2
         x = npos(nodc(imc,j,k),1)
         y = npos(nodc(imc,j,k),2)
         write(io,'(2(g16.8))') x,y
      enddo
c
      close(io)
c
      io = 26
      open(unit=io,file='cell_3.plt',status='unknown')
      j = 1
      do k=1, kmc
         do i=1, imc
            if (nodc(i,j,k).gt.0) then
               y = npos(nodc(i,j,k),2)
               x = npos(nodc(i,j,k),3)
               write(io,'(2(g16.8))') x,y
            endif
         enddo
      enddo
c
      close(io)
c
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      subroutine plot_zone_s(nods,ns1,ns2,ns3,rxy1,rxy2,fi,alfa_y)
c
      implicit none
      include 'common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3)
      integer  i,j,k, io,ksr2
      double precision rxy1(*),rxy2,fi(200,*),alfa_y
      double precision x,y,rho,fi1,psi
      double precision xaxel,yaxel
c
      integer      imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
      io = 26
c      open(unit=io,file='cell_4.plt',status='unknown')
      open(unit=io,file='cell_4.plt',status='scratch')
      do i=1, ims, 2
         write(io,'(t1,a)') 'ZONE'
         do k=1, 11
            psi = real(k-1)/real(10)
            rho = (1.-psi)*rxy1(i) + psi*rxy2
            fi1  = (1.-psi)*fi(i,1) + psi*fi(i,ksr1)
            x = xaxel(rho, fi1, alfa_y)
            y = yaxel(rho, fi1, alfa_y)
            write(io,*) x,y
         enddo
      enddo
      do i=1, 7, 2
         write(io,'(t1,a)') 'ZONE'
         x = xaxel(rxy1(i), fi(i,1), alfa_y)
         y = yaxel(rxy1(i), fi(i,1), alfa_y)
         write(io,*) x,y
         x = xaxel(rxy2, fi(i,ksr1), alfa_y)
         y = yaxel(rxy2, fi(i,ksr1), alfa_y)
         write(io,*) x,y
      enddo
      do i=ims-6, ims, 2
         write(io,'(t1,a)') 'ZONE'
         x = xaxel(rxy1(i), fi(i,1), alfa_y)
         y = yaxel(rxy1(i), fi(i,1), alfa_y)
         write(io,*) x,y
         x = xaxel(rxy2, fi(i,ksr1), alfa_y)
         y = yaxel(rxy2, fi(i,ksr1), alfa_y)
         write(io,*) x,y
      enddo
      close(io)
c
      io = 27
      if (sfred_type.eq.1) then
         ksr2 = ksr1
      elseif ((sfred_type.eq.2).or.(sfred_type.eq.3)) then
         ksr2 = ksr1 + 2
      endif
      open(unit=io,file='cell_5a.plt',status='unknown')
c
      if (sfred_type.eq.1) then
         j = 1
         do i=1, ims, 2
            write(io,'(t1,a)') 'ZONE'
            do k=1, kms, 2
               if (nods(i,j,k).gt.0) then 
                  x = npos(nods(i,j,k),3)
                  y = npos(nods(i,j,k),2)
                  write(io,*) real(x),real(y)
               endif
            enddo
         enddo
c
         do k=1, kms, 2
            write(io,'(t1,a)') 'ZONE'
            do i=1, ims, 2
               if (nods(i,j,k).gt.0) then 
                  x = npos(nods(i,j,k),3)
                  y = npos(nods(i,j,k),2)
                  write(io,*) real(x),real(y)
               endif
            enddo
         enddo
c
      elseif (sfred_type.eq.2) then
         j = 1
         do i=1, ims, 2
            write(io,'(t1,a)') 'ZONE'
            if (mod(i,4).eq.1) then
               do k=1, kms, 2
                  if (nods(i,j,k).gt.0) then 
                     x = npos(nods(i,j,k),3)
                     y = npos(nods(i,j,k),2)
                     write(io,*) real(x),real(y)
                  endif
               enddo
            else
               do k=1, ksr1, 2
                  if (nods(i,j,k).gt.0) then 
                     x = npos(nods(i,j,k),3)
                     y = npos(nods(i,j,k),2)
                     write(io,*) real(x),real(y)
                  endif
               enddo
            endif
         enddo
c
         do k=1, ksr1-2, 2
            write(io,'(t1,a)') 'ZONE'
            do i=1, ims, 2
               if (nods(i,j,k).gt.0) then 
                  x = npos(nods(i,j,k),3)
                  y = npos(nods(i,j,k),2)
                  write(io,*) real(x),real(y)
               endif
            enddo
         enddo
         write(io,'(t1,a)') 'ZONE'
         x = npos(nods(1,j,ksr1),3)
         y = npos(nods(1,j,ksr1),2)
         write(io,*) real(x),real(y)
         do i=5, ims-4, 8
            x = npos(nods(i-2,j,ksr1),3)
            y = npos(nods(i-2,j,ksr1),2)
            write(io,*) real(x),real(y)
            x = npos(nods(i,j,ksr2),3)
            y = npos(nods(i,j,ksr2),2)
            write(io,*) real(x),real(y)
            x = npos(nods(i+2,j,ksr1),3)
            y = npos(nods(i+2,j,ksr1),2)
            write(io,*) real(x),real(y)
            x = npos(nods(i+4,j,ksr1),3)
            y = npos(nods(i+4,j,ksr1),2)
            write(io,*) real(x),real(y)
         enddo
         do k=ksr2, kms, 2
            write(io,'(t1,a)') 'ZONE'
            do i=1, ims, 4
               x = npos(nods(i,j,k),3)
               y = npos(nods(i,j,k),2)
               write(io,*) real(x),real(y)
            enddo
         enddo
c
      elseif (sfred_type.eq.3) then
c
         j = 1
c
         do i=1, ims, 2
            write(io,'(t1,a)') 'ZONE'
            if (mod(i,6).eq.1) then
               do k=1, kms, 2
                  if (nods(i,1,k).gt.0) then 
                     x = npos(nods(i,j,k),3)
                     y = npos(nods(i,j,k),2)
                     write(io,*) real(x),real(y)
                  endif
               enddo
            else
               do k=1, ksr1, 2
                  if (nods(i,1,k).gt.0) then 
                     x = npos(nods(i,j,k),3)
                     y = npos(nods(i,j,k),2)
                     write(io,*) real(x),real(y)
                  endif
               enddo
            endif
         enddo
c
         do k=1, ksr1-2, 2
            write(io,'(t1,a)') 'ZONE'
            do i=1, ims, 2
               if (nods(i,1,k).gt.0) then 
                  x = npos(nods(i,j,k),3)
                  y = npos(nods(i,j,k),2)
                  write(io,*) real(x),real(y)
               endif
            enddo
         enddo
         write(io,'(t1,a)') 'ZONE'
         x = npos(nods(1,j,ksr2),3)
         y = npos(nods(1,j,ksr2),2)
         write(io,*) real(x),real(y)
         do i=4, ims, 6
            x = npos(nods(i-1,j,ksr1),3)
            y = npos(nods(i-1,j,ksr1),2)
            write(io,*) real(x),real(y)
            x = npos(nods(i+1,j,ksr1),3)
            y = npos(nods(i+1,j,ksr1),2)
            write(io,*) real(x),real(y)
            x = npos(nods(i+3,j,ksr2),3)
            y = npos(nods(i+3,j,ksr2),2)
            write(io,*) real(x),real(y)
         enddo
         do k=ksr2, kms, 2
            write(io,'(t1,a)') 'ZONE'
            do i=1, ims, 3
               x = npos(nods(i,j,k),3)
               y = npos(nods(i,j,k),2)
               write(io,*) real(x),real(y)
            enddo
         enddo
         close(io)
         return
c
      endif
c
      close(io)
c
      io = 27
      ksr2 = ksr1 + 2
      open(unit=io,file='cell_5b.plt',status='unknown')
c
      j = jms
c
      do i=1, ims, 2
         write(io,'(t1,a)') 'ZONE'
         if (mod(i,4).eq.1) then
            do k=1, kms, 2
               if (nods(i,1,k).gt.0) then 
                  x = npos(nods(i,j,k),1)
                  y = npos(nods(i,j,k),2)
                  write(io,*) real(x),real(y)
               endif
            enddo
         else
            do k=1, ksr1, 2
               if (nods(i,1,k).gt.0) then 
                  x = npos(nods(i,j,k),1)
                  y = npos(nods(i,j,k),2)
                  write(io,*) real(x),real(y)
               endif
            enddo
         endif
      enddo
c
      do k=1, ksr1-2, 2
         write(io,'(t1,a)') 'ZONE'
         do i=1, ims, 2
            if (nods(i,1,k).gt.0) then 
               x = npos(nods(i,j,k),1)
               y = npos(nods(i,j,k),2)
               write(io,*) real(x),real(y)
            endif
         enddo
      enddo
      write(io,'(t1,a)') 'ZONE'
      x = npos(nods(1,j,ksr1),1)
      y = npos(nods(1,j,ksr1),2)
      write(io,*) real(x),real(y)
      do i=5, ims-4, 8
         x = npos(nods(i-2,j,ksr1),1)
         y = npos(nods(i-2,j,ksr1),2)
         write(io,*) real(x),real(y)
         x = npos(nods(i,j,ksr2),1)
         y = npos(nods(i,j,ksr2),2)
         write(io,*) real(x),real(y)
         x = npos(nods(i+2,j,ksr1),1)
         y = npos(nods(i+2,j,ksr1),2)
         write(io,*) real(x),real(y)
         x = npos(nods(i+4,j,ksr1),1)
         y = npos(nods(i+4,j,ksr1),2)
         write(io,*) real(x),real(y)
      enddo
      do k=ksr2, kms, 2
         write(io,'(t1,a)') 'ZONE'
         do i=1, ims, 4
            x = npos(nods(i,j,k),1)
            y = npos(nods(i,j,k),2)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      close(io)
c
      io = 27
      open(unit=io,file='cell_5c.plt',status='unknown')
c
      do k=1, kms, 2
         write(io,'(t1,a)') 'ZONE'
         do j=1, jms, 2
            if (nods(1,j,k).gt.0) then 
               x = npos(nods(1,j,k),1)
               y = npos(nods(1,j,k),3)
               write(io,*) real(x),real(y)
            endif
         enddo
      enddo
c
      do k=1, kms, 2
         write(io,'(t1,a)') 'ZONE'
         do j=1, jms, 2
            if (nods(ims,j,k).gt.0) then 
               x = npos(nods(ims,j,k),1)
               y = npos(nods(ims,j,k),3)
               write(io,*) real(x),real(y)
            endif
         enddo
      enddo
c
      do j=1, jms, 2
         write(io,'(t1,a)') 'ZONE'
         do k=kms, 1, -2
            if (nods(ims,j,k).gt.0) then 
               x = npos(nods(ims,j,k),1)
               y = npos(nods(ims,j,k),3)
               write(io,*) real(x),real(y)
            endif
         enddo
         do k=1, kms, 2
            if (nods(1,j,k).gt.0) then 
               x = npos(nods(1,j,k),1)
               y = npos(nods(1,j,k),3)
               write(io,*) real(x),real(y)
            endif
         enddo
      enddo
c
      close(io)
c
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      subroutine plot_zone_a(noda,na1,na2,na3)
c
      implicit none
      include 'common_nod.f'
c
      integer  na1,na2,na3,noda(na1,na2,na3)
      integer  i,j,k,ia1,ia2,ka1,ka2,io
      double precision x,y
c
      integer      imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
      ia1 = 2*m1 + 1
      ia2 = 2*(m1+mh+mh) + 1
      ka1 = 2*mv + 1
      if (rtype.eq.0) then
         ka2 = kma
      else
         ka2 = kar1 - 2
      endif
c
      io = 27
      open(unit=io,file='cell_6a.plt',status='unknown')
c
      k=1
c
      do i=1, 3, 2
         write(io,'(t1,a)') 'ZONE'
         do j=1, jma-4, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
      do i=5, ia1, 2
         write(io,'(t1,a)') 'ZONE'
         do j=1, jma, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
      do i=ia2, ima, 2
         write(io,'(t1,a)') 'ZONE'
         do j=1, jma, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      do j=1, jma-4, 2
         write(io,'(t1,a)') 'ZONE'
         do i=1, ia1, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
         do i=ia2, ima, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
      do j=jma-2, jma, 2
         write(io,'(t1,a)') 'ZONE'
         do i=5, ia1, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
         do i=ia2, ima, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      close(io)
c
      io = 28
      open(unit=io,file='cell_6b.plt',status='unknown')
c
      k = ka1
c
      do i=1, 3, 2
         write(io,'(t1,a)') 'ZONE'
         do j=1, jma-4, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
      do i=5, ima, 2
         write(io,'(t1,a)') 'ZONE'
         do j=1, jma, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      do j=1, jma-4, 2
         write(io,'(t1,a)') 'ZONE'
         do i=1, ima, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
      do j=jma-2, jma, 2
         write(io,'(t1,a)') 'ZONE'
         do i=5, ima, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      close(io)
c
      io = 27
      open(unit=io,file='cell_6c.plt',status='unknown')
c
      k = ka2
c
      do i=1, 3, 2
         write(io,'(t1,a)') 'ZONE'
         do j=1, jma-4, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
      do i=5, ima, 2
         write(io,'(t1,a)') 'ZONE'
         do j=1, jma, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      do j=1, jma-4, 2
         write(io,'(t1,a)') 'ZONE'
         do i=1, ima, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
      do j=jma-2, jma, 2
         write(io,'(t1,a)') 'ZONE'
         do i=5, ima, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      close(io)
c
      io = 27
      open(unit=io,file='cell_6d.plt',status='unknown')
c
      k = kma
c
      i = 1
      write(io,'(t1,a)') 'ZONE'
      do j=1, jma-4, 2
         x = npos(noda(i,j,k),1)
         y = npos(noda(i,j,k),3)
         write(io,*) real(x),real(y)
      enddo
      do i=5, ima, 4
         write(io,'(t1,a)') 'ZONE'
         do j=1, jma, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      do j=1, jma-6, 4
         write(io,'(t1,a)') 'ZONE'
         do i=1, ima, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
      do j=jma-4, jma, 4
         write(io,'(t1,a)') 'ZONE'
         do i=5, ima, 2
            x = npos(noda(i,j,k),1)
            y = npos(noda(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      close(io)
c
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      subroutine plot_zone_b(nodb,nb1,nb2,nb3)
c
      implicit none
      include 'common_nod.f'
c
      integer  nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      integer  i,j,k,ka1,ka2,io
      double precision x,y
c
      integer      imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
      ka1 = 2*mv + 1
      if (rtype.eq.0) then
         ka2 = kma
      else
         ka2 = kar1 - 2
      endif
c
      io = 27
      open(unit=io,file='cell_7a.plt',status='unknown')
c
      k=1
c
      do i=1, imb, 2
         write(io,'(t1,a)') 'ZONE'
         do j=1, jmb, 2
            x = npos(nodb(i,j,k),1)
            y = npos(nodb(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      do j=1, jmb, 2
         write(io,'(t1,a)') 'ZONE'
         do i=1, imb, 2
            x = npos(nodb(i,j,k),1)
            y = npos(nodb(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      close(io)
c
      io = 27
      open(unit=io,file='cell_7b.plt',status='unknown')
c
      k=ka1
c
      do i=1, imb, 2
         write(io,'(t1,a)') 'ZONE'
         do j=1, jmb, 2
            x = npos(nodb(i,j,k),1)
            y = npos(nodb(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      do j=1, jmb, 2
         write(io,'(t1,a)') 'ZONE'
         do i=1, imb, 2
            x = npos(nodb(i,j,k),1)
            y = npos(nodb(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      close(io)
c
      io = 27
      open(unit=io,file='cell_7c.plt',status='unknown')
c
      k=ka2
c
      do i=1, imb, 2
         write(io,'(t1,a)') 'ZONE'
         do j=1, jmb, 2
            x = npos(nodb(i,j,k),1)
            y = npos(nodb(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      do j=1, jmb, 2
         write(io,'(t1,a)') 'ZONE'
         do i=1, imb, 2
            x = npos(nodb(i,j,k),1)
            y = npos(nodb(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      close(io)
c
      io = 27
      open(unit=io,file='cell_7d.plt',status='unknown')
c
      k=kma
c
      do i=1, imb, 2
         write(io,'(t1,a)') 'ZONE'
         do j=1, jmb, 2
            x = npos(nodb(i,j,k),1)
            y = npos(nodb(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      do j=1, jmb, 4
         write(io,'(t1,a)') 'ZONE'
         do i=1, imb, 2
            x = npos(nodb(i,j,k),1)
            y = npos(nodb(i,j,k),3)
            write(io,*) real(x),real(y)
         enddo
      enddo
c
      close(io)
c
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
