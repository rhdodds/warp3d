c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine create_plotfil(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                          nodb,nb1,nb2,nb3, job,jobh)
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
      include 'mesh3d_scp_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     &         sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma,
     &         kstart,kstep,n(1000),ia1,ia2,i,j,k,l,
     &         i1,di,j1,dj,kend,ksr1,kar1,kar2,jobh
 
      real*8   r,sf,zmax,zmin,t,w,c,a,
     &         kappa,alfa,r1,r2,rn,eta,my,x1,x2,x9,y1,y2,y9,
     &         z1,z2,z9,g
C     REAL*8 FUNCTIONS !
      real*8   xaxel,zaxel,rhoz,shapefcn
 
      character job*40,fil1p*40,fil1i*40,fil2p*40,fil2i*40,
     &          fil3p*40,fil3i*40,fil4p*40,fil4i*40,fil5p*40,fil5i*40
 
	logical  left
 
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
      common /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
      common /geom/t,w,c,a,kappa,alfa,r1,r2,rn,eta,my
      common /nblock/ kstart,kstep
      common /reduce/ ksr1,kar1,kar2
 
	fil1p=job(1:jobh)//'_1.dat'
	fil1i=job(1:jobh)//'_1.inf'
	fil2p=job(1:jobh)//'_2.dat'
	fil2i=job(1:jobh)//'_2.inf'
	fil3p=job(1:jobh)//'_3.dat'
	fil3i=job(1:jobh)//'_3.inf'
	fil4p=job(1:jobh)//'_4.dat'
	fil4i=job(1:jobh)//'_4.inf'
	fil5p=job(1:jobh)//'_5.dat'
        fil5i=job(1:jobh)//'_5.inf'
 
C---------------------------------
C.......1) Plotfil 1 ( sprickplanet )
C---------------------------------
	do i=1, 1000
	   n(i)=0
	enddo
	open(unit=31,file=fil1p,status='unknown')
C . . . a) RAMEN !
        l=1
	n(1)=5
	write(31,*) 0.0, 0.0
	write(31,*) 0.0, t
	write(31,*) w  , t
	write(31,*) w  , 0.0
	write(31,*) 0.0, 0.0
C . . . b) RHO-LED - ZON S & ZON A !!
	ia1=2*(m1-mh)+1
	ia2=2*(m1+mh)+1
	do j=3, jma-2, 2
	   l=l+1
	   if (mod(j,4).eq.1) then
	      i1=1
	   else
	      i1=3
	   endif
	   do i=i1, 5
	      if (noda(i,j,1).gt.0) then
	         n(l)=n(l)+1
	         write(31,*) npos(noda(i,j,1),1), npos(noda(i,j,1),3)
	      endif
	   enddo
	   do i=6, ia1
	      if (noda(i,j,1).gt.0) then
	         n(l)=n(l)+1
	         write(31,*) npos(noda(i,j,1),1), npos(noda(i,j,1),3)
	      endif
	   enddo
	   do k=kms-1, 3, -1
	      n(l)=n(l)+1
	      write(31,*) npos(nods(ims,j,k),1),npos(nods(ims,j,k),3)
	   enddo
	   do k=2, kms-1
	      n(l)=n(l)+1
	      write(31,*) npos(nods(1,j,k),1),npos(nods(1,j,k),3)
	   enddo
	   do i=ia2, ima
	      n(l)=n(l)+1
	      write(31,*) npos(noda(i,j,1),1), npos(noda(i,j,1),3)
	   enddo
	enddo
 
C . . . c) FI-LED - ZON A !!
	l=l+1
	if (mod(na,4).ne.0) then
	   left=.true.
	   write(31,*) npos(noda(3,1,3),1), npos(noda(3,1,3),3)
	else
	   left=.false.
	   write(31,*) npos(noda(1,1,3),1), npos(noda(1,1,3),3)
	endif
        n(l)=n(l)+1
	do j=3, jma-4, 4
	   if (left) then
	      write(31,*) npos(noda(3,j,3),1), npos(noda(3,j,3),3)
	      write(31,*) npos(noda(1,j+2,3),1), npos(noda(1,j+2,3),3)
              n(l)=n(l)+2
	      left=.false.
	   else
	      write(31,*) npos(noda(3,j,3),1), npos(noda(3,j,3),3)
	      write(31,*) npos(noda(3,j+2,3),1), npos(noda(3,j+2,3),3)
              n(l)=n(l)+2
	      left=.true.
	   endif
	enddo
	do i=5, ima-2, 2
	  if (noda(i,1,1).gt.0) then
	     l=l+1
	     do j=1, jma
	        if (noda(i,j,1).gt.0) then
	          write(31,*) npos(noda(i,j,1),1), npos(noda(i,j,1),3)
		  n(l)=n(l)+1
	        endif
	     enddo
	  endif
	enddo
C . . . d) FI-LED - ZON S !!
	do k=kms-2, 3, -2
	   l=l+1
	   n(l)=91
	   r=rhoz(npos(nods(ims,1,k),3),90.d0,alfa)
	   do j=0, 90
	      write(31,*) xaxel(r,dble(j),alfa),zaxel(r,dble(j),alfa)
	   enddo
	enddo
	do i=1, ims-4, 2
	   l=l+1
	   n(l)=91
	   r=rhoz(npos(nods(i,1,1),3),90.d0,alfa)
	   do j=0, 90
	      write(31,*) xaxel(r,dble(j),alfa),zaxel(r,dble(j),alfa)
	   enddo
	enddo
	do k=3, kms-2, 2
	   l=l+1
	   n(l)=91
	   r=rhoz(npos(nods(1,1,k),3),90.d0,alfa)
	   do j=0, 90
	      write(31,*) xaxel(r,dble(j),alfa),
     &                    zaxel(r,dble(j),alfa)
	   enddo
	enddo
C . . . e) ZON B !!
C . . . Z-LED . . . . . .
	do i=1, imb-2, 2
	   l=l+1
	   do j=1, jmb
	      write(31,*) npos(nodb(i,j,1),1), npos(nodb(i,j,1),3)
	      n(l)=n(l)+1
	   enddo
	enddo
C . . . X-LED . . . . .
	do j=3, jmb-2, 2
	   l=l+1
	   do i=1, imb
	      write(31,*) npos(nodb(i,j,1),1), npos(nodb(i,j,1),3)
	      n(l)=n(l)+1
	   enddo
	enddo
	close(31)
C.......Skapa "info" filen FIL1I .
	open(unit=32,file=fil1i,status='unknown')
        write(32,'(t1,a)') '*CURVE1'
 	if (t.gt.w) then
	   sf=t/20.d0
	else
	   sf=w/20.d0
	endif
	write(32,'(t1,2f5.1,tr2,6(tr1,g10.4))')
     &  1.5,1.5,-sf,21.*sf,-sf,21.*sf,18.0,18.0
	write(32,'(t1,2(a,g10.4,tr2,i1),tr2,f4.2)')
     &  '0.0 ', w/2., 1, ' 0.0 ', t, 1, 0.35
	write(32,'(t1,i3)') l
	do i=1, l
	   write(32,'(t1,i4,a)') n(i),' 3 0 401 0 1'
	enddo
 	write(32,'(7(t1,a/),t1,a)') '*SCALE ', '0 2700 0 1900 0',
     &  '*CHRDEF', '40 10 0.05', '*TEXT  ','400 1850 0.3 0. 1',
     &   job(1:jobh), '*FINISH'
	close(32)
 
C---------------------------------
C PLOTTFIL 2 ( symmetriplanet )
C---------------------------------
	do i=1, 1000
	   n(i)=0
	enddo
	open(unit=41,file=fil2p,status='unknown')
C . . . a) ZON S !
C . . . CIRKLAR !
	l=1
	write(41,*) -npos(nods(1,1,1),2), npos(nods(1,1,1),3)
	n(l)=n(l)+1
	do i=1, ims-4, 2
	   x1=npos(nods(i,1,1),2)
	   x2=npos(nods(i+2,1,1),2)
	   x9=npos(nods(i+1,1,1),2)
	   z1=npos(nods(i,1,1),3)
	   z2=npos(nods(i+2,1,1),3)
	   z9=npos(nods(i+1,1,1),3)
	   do j1=1, 10
	      g=0.2d0*dble(j1-5)
	      write(41,*) -shapefcn(g,x1,x2,x9),shapefcn(g,z1,z2,z9)
	      n(l)=n(l)+1
	   enddo
	enddo
	if (sfred.eq.0) then
C . . . . "CIRKLAR-ELLIPSER" !
	   do k=3, kms, 2
	      l=l+1
	      do i=1, ims
	         write(41,*) -npos(nods(i,1,k),2),npos(nods(i,1,k),3)
	         n(l)=n(l)+1
	      enddo
	   enddo
C . . .    RADIER !
	   do i=1, ims, 2
	      l=l+1
	      do k=1, kms
	         if (nods(i,1,k).gt.0) then
	           write(41,*) -npos(nods(i,1,k),2),npos(nods(i,1,k),3)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	else
C . . . . "CIRKLAR-ELLIPSER" !
	   do k=3, ksr1-2, 2
	      l=l+1
	      do i=1, ims
	         write(41,*) -npos(nods(i,1,k),2), npos(nods(i,1,k),3)
	         n(l)=n(l)+1
	      enddo
	   enddo
	   l=l+1
	   write(41,*) -npos(nods(1,1,ksr1),2),npos(nods(1,1,ksr1),3)
	   n(l)=1
           do i=2, ims, 8
	    write(41,*) -npos(nods(i,1,ksr1),2),npos(nods(i,1,ksr1),3)
	     write(41,*) -npos(nods(i+1,1,ksr1),2),
     &                       npos(nods(i+1,1,ksr1),3)
	     write(41,*) -npos(nods(i+2,1,ksr1+1),2),
     &                       npos(nods(i+2,1,ksr1+1),3)
             write(41,*) -npos(nods(i+3,1,ksr1+2),2),
     &                       npos(nods(i+3,1,ksr1+2),3)
	     write(41,*) -npos(nods(i+4,1,ksr1+1),2),
     &                       npos(nods(i+4,1,ksr1+1),3)
	     write(41,*) -npos(nods(i+5,1,ksr1),2),
     &                       npos(nods(i+5,1,ksr1),3)
	     write(41,*) -npos(nods(i+6,1,ksr1),2),
     &                       npos(nods(i+6,1,ksr1),3)
	     write(41,*) -npos(nods(i+7,1,ksr1),2),
     &                       npos(nods(i+7,1,ksr1),3)
	     n(l)=n(l)+8
	   enddo
	   do k=ksr1+2, kms, 2
	      l=l+1
	      do i=1, ims, 2
	         write(41,*) -npos(nods(i,1,k),2), npos(nods(i,1,k),3)
	         n(l)=n(l)+1
	      enddo
	   enddo
C . . .    RADIER !
	   do i=1, ims, 2
	      l=l+1
	      if (mod(i,4).eq.1) then
	         do k=1, kms
	            if (nods(i,1,k).gt.0) then
	               write(41,*) -npos(nods(i,1,k),2),
     &                              npos(nods(i,1,k),3)
	               n(l)=n(l)+1
	            endif
	         enddo
	      else
	         do k=1, ksr1
	            if (nods(i,1,k).gt.0) then
	               write(41,*) -npos(nods(i,1,k),2),
     &                              npos(nods(i,1,k),3)
	               n(l)=n(l)+1
	            endif
	         enddo
	      endif
	   enddo
	endif
C . ZON A
C . Z-direction
	do k=1, 2*mv+1, 2
	   l=l+1
	   do i=1, ia1
	      if (noda(i,1,k).gt.0) then
	         write(41,*) -npos(noda(i,1,k),2), npos(noda(i,1,k),3)
	         n(l)=n(l)+1
	      endif
	   enddo
	   l=l+1
	   do i=ia2, ima
	      if (noda(i,1,k).gt.0) then
	         write(41,*) -npos(noda(i,1,k),2), npos(noda(i,1,k),3)
	         n(l)=n(l)+1
	      endif
	   enddo
	enddo
	do k=2*mv+3, kar1-2, 2
	   l=l+1
	   do i=1, ima
 	     if (noda(i,1,k).gt.0) then
	       write(41,*) -npos(noda(i,1,k),2), npos(noda(i,1,k),3)
	       n(l)=n(l)+1
	     endif
	   enddo
	enddo
	if (rtype.eq.0) then
	   do k=kar1, kma, 2
	      l=l+1
	      do i=1, ima
 	         if (noda(i,1,k).gt.0) then
	           write(41,*) -npos(noda(i,1,k),2),npos(noda(i,1,k),3)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	else
	   l=l+1
	   if (mod(ma,4).eq.0) then
	      write(41,*) -npos(noda(1,1,kar1),2),npos(noda(1,1,kar1),3)
	      write(41,*) -npos(noda(5,1,kar1),2),npos(noda(5,1,kar1),3)
	      n(l)=n(l)+2
	      do i=9, ima, 8
	         write(41,*) -npos(noda(i-2,1,kar1),2),
     &                        npos(noda(i-2,1,kar1),3)
	         write(41,*) -npos(noda(i,1,kar1+2),2),
     &                        npos(noda(i,1,kar1+2),3)
	         n(l)=n(l)+2
	         if (noda(i+2,1,kar1).gt.0) then
	            write(41,*) -npos(noda(i+2,1,kar1),2),
     &                           npos(noda(i+2,1,kar1),3)
	            write(41,*) -npos(noda(i+4,1,kar1),2),
     &                           npos(noda(i+4,1,kar1),3)
	            n(l)=n(l)+2
	         endif
	      enddo
	   else
	      write(41,*) -npos(noda(5,1,kar1+2),2),
     &                     npos(noda(5,1,kar1+2),3)
	      n(l)=n(l)+1
	      do i=9, ima, 8
	         write(41,*) -npos(noda(i-2,1,kar1),2),
     &                        npos(noda(i-2,1,kar1),3)
	         write(41,*) -npos(noda(i,1,kar1),2),
     &                        npos(noda(i,1,kar1),3)
	         write(41,*) -npos(noda(i+2,1,kar1),2),
     &                        npos(noda(i+2,1,kar1),3)
	         write(41,*) -npos(noda(i+4,1,kar1+2),2),
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
	      do i=i1, ima
 	         if (noda(i,1,k).gt.0) then
	           write(41,*) -npos(noda(i,1,k),2),npos(noda(i,1,k),3)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	endif
 
C . Y-direction
	l=l+1
	do k=1, kma
	   if (noda(1,1,k).gt.0) then
              write(41,*)  -npos(noda(1,1,k),2), npos(noda(1,1,k),3)
	      n(l)=n(l)+1
	   endif
	enddo
	if (rtype.lt.2) then
	   if (mod(na,4).ne.0) then
	      l=l+1
	      do k=1, kma
	         if (noda(3,1,k).gt.0) then
                    write(41,*)  -npos(noda(3,1,k),2),
     &                            npos(noda(3,1,k),3)
	            n(l)=n(l)+1
	         endif
	      enddo
	   endif
	else
	   if (mod(na,4).ne.0) then
	      l=l+1
	      do k=1, kar2
	         if (noda(3,1,k).gt.0) then
                    write(41,*)  -npos(noda(3,1,k),2),
     &                            npos(noda(3,1,k),3)
	            n(l)=n(l)+1
 	         endif
	      enddo
	      write(41,*) -npos(noda(1,1,kar2+2),2),
     &                     npos(noda(1,1,kar2+2),3)
	      n(l)=n(l)+1
	   endif
	endif
	if (rtype.eq.0) then
	   do i=5, ima, 2
	      l=l+1
 	      do k=1, kma
	         if (noda(i,1,k).gt.0) then
	           write(41,*) -npos(noda(i,1,k),2),npos(noda(i,1,k),3)
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
 	      do k=1, kend
	         if (noda(i,1,k).gt.0) then
	           write(41,*) -npos(noda(i,1,k),2),npos(noda(i,1,k),3)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	endif
	close(41)
C.......Skapa "info" filen FIL2I .
	zmax=t
	zmin=0.
	do k=1, kma, 2
          if (npos(noda(ima,1,k),3).gt.zmax) zmax=npos(noda(ima,1,k),3)
          if (npos(noda(1,1,k),3).lt.zmin) zmin=npos(noda(1,1,k),3)
	enddo
	open(unit=42,file=fil2i,status='unknown')
 	write(42,'(t1,a)') '*CURVE1'
	if ((zmax-zmin).gt.npos(noda(1,1,kma),2) ) then
	   sf=(zmax-zmin)/20.
	else
	   sf=npos(noda(1,1,kma),2)/20.
	endif
	write(42,'(t1,2f5.1,tr2,6(tr1,g10.4))')
     &  1.5,1.5,-21.*sf,sf,-sf+zmin,21.*sf+zmin,18.0 , 18.0
	write(42,'(t1,2(g12.6,tr2),i1,a,g10.4,tr2,i1,tr2,f4.2)')
     &  -1.*npos(noda(1,1,kma),2), npos(noda(1,1,kma),2)/2., 1,
     &  ' 0.0 ', (zmax-zmin), 1, 0.35
	write(42,'(t1,i3)') l
	do i=1, l
	   write(42,'(t1,i4,a)') n(i),' 3 0 401 0 1'
	enddo
 	write(42,'(7(t1,a/),t1,a)') '*SCALE ', '0 2700 0 1900 0',
     &  '*CHRDEF', '40 10 0.05', '*TEXT  ','400 1850 0.3 0. 1',
     &   job(1:jobh), '*FINISH'
	close(42)
 
C---------------------------------
C Plotfil 3 ( det bakre planet )
C---------------------------------
	do i=1, 1000
	   n(i)=0
	enddo
	open(unit=31,file=fil3p,status='unknown')
C . . . a) RAMEN !
	zmax=npos(noda(ima,1,kma),3)
	zmin=npos(noda(1,1,kma),3)
        l=1
	n(1)=5
	write(31,*) 0.0, zmin
	write(31,*) 0.0, zmax
	write(31,*) w  , zmax
	write(31,*) w  , zmin
	write(31,*) 0.0, zmin
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
	        write(31,*) npos(noda(i,j,kma),1),npos(noda(i,j,kma),3)
	      endif
	   enddo
	enddo
 
C . . . c) FI-LED - ZON A !!
	if (rtype.le.1) then
	   l=l+1
  	   if (mod(na,4).ne.0) then
	      left=.true.
	      write(31,*) npos(noda(3,1,kma),1), npos(noda(3,1,kma),3)
	   else
	      left=.false.
	      write(31,*) npos(noda(1,1,kma),1), npos(noda(1,1,kma),3)
	   endif
           n(l)=n(l)+1
	   do j=3, jma-4, 4
	      if (left) then
	       write(31,*) npos(noda(3,j,kma),1), npos(noda(3,j,kma),3)
	   write(31,*) npos(noda(1,j+2,kma),1), npos(noda(1,j+2,kma),3)
                n(l)=n(l)+2
	        left=.false.
	      else
	        write(31,*) npos(noda(3,j,kma),1), npos(noda(3,j,kma),3)
	   write(31,*) npos(noda(3,j+2,kma),1), npos(noda(3,j+2,kma),3)
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
	         write(31,*)npos(noda(i,j,kma),1),npos(noda(i,j,kma),3)
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
	        write(31,*) npos(nodb(i,j,kma),1),npos(nodb(i,j,kma),3)
	        n(l)=n(l)+1
	      endif
	   enddo
	enddo
C . . . X-LED . . . . .
	do j=dj+1, jmb-dj, dj
	   l=l+1
	   do i=1, imb
	      write(31,*) npos(nodb(i,j,kma),1),npos(nodb(i,j,kma),3)
	      n(l)=n(l)+1
	   enddo
	enddo
	close(31)
C.......Skapa "info" filen FIL3I .
	open(unit=32,file=fil3i,status='unknown')
        write(32,'(t1,a)') '*CURVE1'
	if ((zmax-zmin).gt.w) then
	   sf=(zmax-zmin)/20.
	else
	   sf=w/20.
	endif
	write(32,'(t1,2f5.1,tr2,6(tr1,g10.4))')
     &  1.5,1.5,-sf,21.*sf,-sf+zmin,21.*sf+zmin,18.,18.
	write(32,'(t1,2(a,g10.4,tr2,i1),tr2,f4.2)')
     &  '0.0 ', w/2., 1, ' 0.0 ', (zmax-zmin), 1, 0.35
	write(32,'(t1,i3)') l
	do i=1, l
	   write(32,'(t1,i4,a)') n(i),' 3 0 401 0 1'
	enddo
        write(32,'(7(t1,a/),t1,a)') '*SCALE ', '0 2700 0 1900 0',
     &  '*CHRDEF', '40 10 0.05', '*TEXT  ','400 1850 0.3 0. 1',
     &   job(1:jobh), '*FINISH'
	close(32)
 
C---------------------------------
C Plotfil 4 ( ovansidan )
C---------------------------------
	do i=1, 1000
	   n(i)=0
	enddo
	l=0
	open(unit=33,file=fil4p,status='unknown')
C Zon A :
C . .  Y-direction
       	if (rtype.le.1) then
	   do j=1, jma-jmb+1, 2
	      l=l+1
	      do k=1, kma
	         if (noda(ima,j,k).gt.0) then
	            write(33,*) npos(noda(ima,j,k),1),
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
	            write(33,*) npos(noda(ima,j,k),1),
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
	            write(33,*) npos(noda(ima,j,k),1),
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
	            write(33,*) npos(noda(ima,j,k),1),
     &                          npos(noda(ima,j,k),2)
	            n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	   if (mod(na,4).ne.0) then
	      l=l+1
              write(33,*) npos(noda(ima,1,kar2),1),
     &                    npos(noda(ima,1,kar2),2)
              write(33,*) npos(noda(ima,3,kar2),1),
     &                    npos(noda(ima,3,kar2),2)
              write(33,*) npos(noda(ima,5,kar2+2),1),
     &                    npos(noda(ima,5,kar2+2),2)
	      n(l)=n(l)+3
	      do j=9, jma-jmb+1, 8
	         write(33,*) npos(noda(ima,j-2,kar2),1),
     &                       npos(noda(ima,j-2,kar2),2)
	         write(33,*) npos(noda(ima,j,kar2),1),
     &                       npos(noda(ima,j,kar2),2)
	         n(l)=n(l)+2
                 if ((jma-jmb+1).gt.9) then
	            write(33,*) npos(noda(ima,j+2,kar2),1),
     &                          npos(noda(ima,j+2,kar2),2)
	            write(33,*) npos(noda(ima,j+4,kar2+2),1),
     &                          npos(noda(ima,j+4,kar2+2),2)
	            n(l)=n(l)+2
                 endif
	      enddo
	   else
	      l=l+1
              write(33,*) npos(noda(ima,1,kar2+2),1),
     &                    npos(noda(ima,1,kar2+2),2)
	      n(l)=n(l)+1
	      do j=5, jma-jmb+1, 8
	         write(33,*) npos(noda(ima,j-2,kar2),1),
     &                       npos(noda(ima,j-2,kar2),2)
	         write(33,*) npos(noda(ima,j,kar2),1),
     &                       npos(noda(ima,j,kar2),2)
	         write(33,*) npos(noda(ima,j+2,kar2),1),
     &                       npos(noda(ima,j+2,kar2),2)
	         write(33,*) npos(noda(ima,j+4,kar2+2),1),
     &                       npos(noda(ima,j+4,kar2+2),2)
	         n(l)=n(l)+4
	      enddo
	   endif
	   do k=kar2+2, kma, 2
	      l=l+1
	      do j=1, jma-jmb+1
	         if (noda(ima,j,k).gt.0) then
	            write(33,*) npos(noda(ima,j,k),1),
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
	         write(33,*) npos(nodb(i,1,k),1), npos(nodb(i,1,k),2)
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
	            write(33,*) npos(nodb(i,1,k),1), npos(nodb(i,1,k),2)
	            n(l)=n(l)+1
	         endif
	      enddo
	   endif
	enddo
	close(33)
C.......Skapa "info" filen FIL4I .
	open(unit=34,file=fil4i,status='unknown')
        write(34,'(t1,a)') '*CURVE1'
	if (npos(noda(ima,1,kma),2).gt.w) then
	   sf=npos(noda(ima,1,kma),2)/20.
	else
	   sf=w/20.
	endif
	write(34,'(t1,2f5.1,tr2,6(tr1,g10.4))')
     &  1.5,1.5,-sf,21.*sf,-sf,21.*sf,18.,18.
	write(34,'(t1,2(a,g10.4,tr2,i1),tr2,f4.2)')
     &  '0.0 ', w/4., 1, ' 0.0 ', npos(noda(ima,1,kma),2)/4., 1, 0.35
	write(34,'(t1,i3)') l
	do i=1, l
	   write(34,'(t1,i4,a)') n(i),' 3 0 401 0 1'
	enddo
         write(34,'(7(t1,a/),t1,a)') '*SCALE ', '0 2700 0 1900 0',
     &  '*CHRDEF', '40 10 0.05', '*TEXT  ','400 1850 0.3 0. 1',
     &   job(1:jobh), '*FINISH'
	close(34)
 
C---------------------------------
C Plotfil 5 ( undersidan )
C---------------------------------
	do i=1, 1000
	   n(i)=0
	enddo
	open(unit=35,file=fil5p,status='unknown')
C . . . a) ZON S !
C . . . CIRKLAR !
	l=1
	write(35,*) npos(nods(1,jms,1),1), npos(nods(1,jms,1),2)
	n(l)=n(l)+1
	do i=1, ims-4, 2
	   x1=npos(nods(i,jms,1),1)
	   x2=npos(nods(i+2,jms,1),1)
	   x9=npos(nods(i+1,jms,1),1)
	   y1=npos(nods(i,jms,1),2)
	   y2=npos(nods(i+2,jms,1),2)
	   y9=npos(nods(i+1,jms,1),2)
	   do j1=1, 10
	      g=0.2d0*dble(j1-5)
	      write(35,*) shapefcn(g,x1,x2,x9),shapefcn(g,y1,y2,y9)
	      n(l)=n(l)+1
	   enddo
	enddo
	if (sfred.eq.0) then
C . . . . "CIRKLAR-ELLIPSER" !
	   do k=3, kms, 2
	      l=l+1
	      do i=1, ims
	         write(35,*) npos(nods(i,jms,k),1),npos(nods(i,jms,k),2)
	         n(l)=n(l)+1
	      enddo
	   enddo
C . . .    RADIER !
	   do i=1, ims, 2
	      l=l+1
	      do k=1, kms
	         if (nods(i,jms,k).gt.0) then
	         write(35,*) npos(nods(i,jms,k),1),npos(nods(i,jms,k),2)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	else
C . . . . "CIRKLAR-ELLIPSER" !
	   do k=3, ksr1-2, 2
	      l=l+1
	      do i=1, ims
	         write(35,*) npos(nods(i,jms,k),1),npos(nods(i,jms,k),2)
	         n(l)=n(l)+1
	      enddo
	   enddo
	   l=l+1
	   write(35,*) npos(nods(1,jms,ksr1),1),npos(nods(1,jms,ksr1),2)
	   n(l)=1
           do i=2, ims, 8
	     write(35,*) npos(nods(i,jms,ksr1),1),
     &                       npos(nods(i,jms,ksr1),2)
	     write(35,*) npos(nods(i+1,jms,ksr1),1),
     &                       npos(nods(i+1,jms,ksr1),2)
	     write(35,*) npos(nods(i+2,jms,ksr1+1),1),
     &                       npos(nods(i+2,jms,ksr1+1),2)
             write(35,*) npos(nods(i+3,jms,ksr1+2),1),
     &                       npos(nods(i+3,jms,ksr1+2),2)
	     write(35,*) npos(nods(i+4,jms,ksr1+1),1),
     &                       npos(nods(i+4,jms,ksr1+1),2)
	     write(35,*) npos(nods(i+5,jms,ksr1),1),
     &                       npos(nods(i+5,jms,ksr1),2)
	     write(35,*) npos(nods(i+6,jms,ksr1),1),
     &                       npos(nods(i+6,jms,ksr1),2)
	     write(35,*) npos(nods(i+7,jms,ksr1),1),
     &                       npos(nods(i+7,jms,ksr1),2)
	     n(l)=n(l)+8
	   enddo
	   do k=ksr1+2, kms, 2
	      l=l+1
	      do i=1, ims, 2
	         write(35,*) npos(nods(i,jms,k),1),npos(nods(i,jms,k),2)
	         n(l)=n(l)+1
	      enddo
	   enddo
C . . .    RADIER !
	   do i=1, ims, 2
	      l=l+1
	      if (mod(i,4).eq.1) then
	         do k=1, kms
	            if (nods(i,jms,k).gt.0) then
	               write(35,*) npos(nods(i,jms,k),1),
     &                              npos(nods(i,jms,k),2)
	               n(l)=n(l)+1
	            endif
	         enddo
	      else
	         do k=1, ksr1
	            if (nods(i,jms,k).gt.0) then
	               write(35,*) npos(nods(i,jms,k),1),
     &                              npos(nods(i,jms,k),2)
	               n(l)=n(l)+1
	            endif
	         enddo
	      endif
	   enddo
	endif
 
 
C Zon A :
C . .  X-direction, I=1
	do k=1, kma, 2
	   if (.not.( (rtype.ge.1).and.(mod(ma,4).ne.0).and.
     &           (k.eq.kar1) ) ) then
	      l=l+1
	      do j=1, jma-4
	         if (noda(1,j,k).gt.0) then
	            write(35,*) npos(noda(1,j,k),1),npos(noda(1,j,k),2)
	            n(l)=n(l)+1
	         endif
	      enddo
	   endif
	enddo
C . .  Y-direction, I=1
	do j=1, jma-4, 4
	   l=l+1
	   do k=1, kma
	      if (noda(1,j,k).gt.0) then
	         write(35,*) npos(noda(1,j,k),1),npos(noda(1,j,k),2)
	         n(l)=n(l)+1
	      endif
	   enddo
	enddo
 
C . . X-direction, I=IMA
	do k=1, 2*mv+1, 2
	   if (ia1.gt.5) then
	      l=l+1
	      do i=5, ia1
	         if (noda(i,jma,k).gt.0) then
	            write(35,*) npos(noda(i,jma,k),1),
     &                           npos(noda(i,jma,k),2)
	            n(l)=n(l)+1
	         endif
	      enddo
	   endif
	   l=l+1
	   do i=ia2, ima
	      if (noda(i,jma,k).gt.0) then
	         write(35,*) npos(noda(i,jma,k),1),npos(noda(i,jma,k),2)
	         n(l)=n(l)+1
	      endif
	   enddo
	enddo
	do k=2*mv+3, kar1-2, 2
	   l=l+1
	   do i=5, ima
 	     if (noda(i,jma,k).gt.0) then
	       write(35,*) npos(noda(i,jma,k),1), npos(noda(i,jma,k),2)
	       n(l)=n(l)+1
	     endif
	   enddo
	enddo
	if (rtype.eq.0) then
	   do k=kar1, kma, 2
	      l=l+1
	      do i=5, ima
 	         if (noda(i,jma,k).gt.0) then
	        write(35,*) npos(noda(i,jma,k),1),npos(noda(i,jma,k),2)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	else
	   l=l+1
	   if (mod(ma,4).eq.0) then
	      write(35,*) npos(noda(5,jma,kar1),1),
     &                    npos(noda(5,jma,kar1),2)
	      n(l)=n(l)+1
	      do i=9, ima, 8
	         write(35,*) npos(noda(i-2,jma,kar1),1),
     &                        npos(noda(i-2,jma,kar1),2)
	         write(35,*) npos(noda(i,jma,kar1+2),1),
     &                        npos(noda(i,jma,kar1+2),2)
	         n(l)=n(l)+2
C	         IF (NODA(I+2,JMA,KAR1).GT.0) THEN
                 if (i.lt.ima) then
	            write(35,*) npos(noda(i+2,jma,kar1),1),
     &                           npos(noda(i+2,jma,kar1),2)
	            write(35,*) npos(noda(i+4,jma,kar1),1),
     &                           npos(noda(i+4,jma,kar1),2)
	            n(l)=n(l)+2
	         endif
	      enddo
	   else
	      write(35,*) npos(noda(5,jma,kar1+2),1),
     &                    npos(noda(5,jma,kar1+2),2)
	      n(l)=n(l)+1
	      do i=9, ima, 8
	         write(35,*) npos(noda(i-2,jma,kar1),1),
     &                        npos(noda(i-2,jma,kar1),2)
	         write(35,*) npos(noda(i,jma,kar1),1),
     &                        npos(noda(i,jma,kar1),2)
	         write(35,*) npos(noda(i+2,jma,kar1),1),
     &                        npos(noda(i+2,jma,kar1),2)
	         write(35,*) npos(noda(i+4,jma,kar1+2),1),
     &                        npos(noda(i+4,jma,kar1+2),2)
	         n(l)=n(l)+4
	      enddo
	   endif
	   do k=kar1+2, kma, 2
	      l=l+1
	      do i=5, ima
 	         if (noda(i,jma,k).gt.0) then
	        write(35,*) npos(noda(i,jma,k),1),npos(noda(i,jma,k),2)
	           n(l)=n(l)+1
	         endif
	      enddo
	   enddo
	endif
C . .  Y-direction, I=1
	if (rtype.eq.0) then
	   do i=7, ima, 2
	      l=l+1
 	      do k=1, kma
	         if (noda(i,jma,k).gt.0) then
	        write(35,*) npos(noda(i,jma,k),1),npos(noda(i,jma,k),2)
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
 	      do k=1, kend
	         if (noda(i,jma,k).gt.0) then
	       write(35,*) npos(noda(i,jma,k),1),npos(noda(i,jma,k),2)
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
	      do i=1, imb
	         if (nodb(i,jmb,k).gt.0) then
	        write(35,*) npos(nodb(i,jmb,k),1),npos(nodb(i,jmb,k),2)
	            n(l)=n(l)+1
	         endif
	      enddo
	   endif
	enddo
C . .  Y-direction, I=1
	do i=3, imb, 2
	   l=l+1
	   do k=1, kma
	      if (nodb(i,jmb,k).gt.0) then
	         write(35,*) npos(nodb(i,jmb,k),1),npos(nodb(i,jmb,k),2)
	         n(l)=n(l)+1
	      endif
	   enddo
	enddo
	close(35)
C.......Skapa "info" filen FIL5I .
	open(unit=36,file=fil5i,status='unknown')
        write(36,'(t1,a)') '*CURVE1'
	if (npos(noda(1,1,kma),2).gt.w) then
	   sf=npos(noda(1,1,kma),2)/20.
	else
	   sf=w/20.
	endif
	write(36,'(t1,2f5.1,tr2,6(tr1,g10.4))')
     &  1.5,1.5,-sf,21.*sf,-sf,21.*sf,18.,18.
	write(36,'(t1,2(a,g10.4,tr2,i1),tr2,f4.2)')
     &  '0.0 ', w/4., 1, ' 0.0 ', npos(noda(1,1,kma),2)/4., 1, 0.35
	write(36,'(t1,i3)') l
	do i=1, l
	   write(36,'(t1,i4,a)') n(i),' 3 0 401 0 1'
	enddo
        write(36,'(7(t1,a/),t1,a)') '*SCALE ', '0 2700 0 1900 0',
     &  '*CHRDEF', '40 10 0.05', '*TEXT  ','400 1850 0.3 0. 1',
     &   job(1:jobh), '*FINISH'
	close(36)
 
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	real*8 function shapefcn(g,u1,u2,u9)
C----------------------------------------------------------------C
C This routine calaculates the physical coordinate of the iso-   C
C parametric variable G. U1, U2 and U9 are the nodal values in   C
C a solid 20-noded brick-element.                                C
C----------------------------------------------------------------C
	implicit none
	real*8 g,u1,u2,u9
	shapefcn=-0.5d0*g*(1.d0-g)*u1+0.5d0*g*(1.d0+g)*u2+
     &              (1.d0-g)*(1.d0+g)*u9
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
C----------------------------------------------------------------C
C        NEDANSAENDE SUBROUTINER ANVANDS ENBART VID              C
C               TEST AV SJALVA PROGRAMMET  !!!!!                 C
C----------------------------------------------------------------C
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine plot_nod_zons(nods,ns1,ns2,ns3,fi)
C----------------------------------------------------------C
C  Rutinen skapar en plotfil f|r noderna i Zon S.          C
C----------------------------------------------------------C
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     &         sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma,
     &         ns1,ns2,ns3,nods(ns1,ns2,ns3),ant1,ant2,n,i,j,k
      real*8   fi(50)
 
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
      common /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
 
	open(unit=14,file='zons_90c.plo',status='unknown')
	open(unit=15,file='zons_90m.plo',status='unknown')
	ant1=0
	ant2=0
	do i=1,ims
	   do k=1,kms
	      if (nods(i,1,k).gt.0) then
	        write (14,*) npos(nods(i,1,k),3),npos(nods(i,1,k),2)
	        ant1=ant1+1
	      endif
	      if (nods(i,2,k).gt.0) then
	        write(15,*) npos(nods(i,2,k),3),npos(nods(i,2,k),2)
	        ant2=ant2+1
	      endif
	   enddo
	enddo
	close(14)
	close(15)
	write(*,'(t10,a,i3/t10,a,i3)')
     &   '* ant. koordinater i filen "zons_90c.plo" : ', ant1,
     &   '* ant. koordinater i filen "zons_90m.plo" : ', ant2
 
        open(unit=14,file='zons_30zc.plo',status='unknown')
	open(unit=15,file='zons_30xc.plo',status='unknown')
	open(unit=16,file='zons_30xm.plo',status='unknown')
        ant1=0
        ant2=0
	j=1
10 	continue
	if (fi(j).gt.30.) then
	   j=j+2
	   goto 10
	endif
	n=j
	do i=1,ims
	   do k=1,kms
	      if (nods(i,n,k).gt.0) then
	         write(14,*) npos(nods(i,n,k),3),npos(nods(i,n,k),2)
	         write(15,*) npos(nods(i,n,k),1),npos(nods(i,n,k),2)
	         ant1=ant1+1
	      endif
	      if (nods(i,n-1,k).gt.0) then
	        write(16,*) npos(nods(i,n-1,k),1),npos(nods(i,n-1,k),2)
	        ant2=ant2+1
	      endif
	   enddo
	enddo
	close(14)
	close(15)
	close(16)
 
	write(*,'(t10,a,a,i3/t12,a,i3)') '* ant. koordinater i filen',
     &  ' "zons_30zc.plo" resp "zons_30xc.plo" : ', ant1,
     &  'ant. koord. ifilen "zons_30xm.plo" : ', ant2
 
	open(unit=14,file='zons_00c.plo',status='unknown')
	open(unit=15,file='zons_00m.plo',status='unknown')
	ant1=0
	ant2=0
	do i=1,ims
	  do k=1,kms
	    if (nods(i,jma,k).gt.0) then
	       write(14,*) npos(nods(i,jma,k),1),npos(nods(i,jma,k),2)
	       ant1=ant1+1
	    endif
	    if (nods(i,jma-1,k).gt.0) then
	     write(15,*)npos(nods(i,jma-1,k),1),npos(nods(i,jma-1,k),2)
	       ant2=ant2+1
	    endif
	  enddo
	enddo
	close(14)
	close(15)
	write(*,'(t10,a,i3/t10,a,i3)')
     &  '* ant. koordinater i filen "zons_00c.plo" : ', ant1,
     &  '* ant. koordinater i filen "zons_00m.plo" : ', ant2
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine plot_nod_zona(noda,na1,na2,na3)
C----------------------------------------------------------C
C  Rutinen skapar en plotfil f|r noderna i Zon A.          C
C----------------------------------------------------------C
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     &         sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma,
     &         na1,na2,na3,noda(na1,na2,na3),n,i,j,k
      character  fil*13
 
      common  /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
      common  /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
 
	write(*,*) ' kma = ',kma
 
	open(unit=14,file='zona_yz.plo',status='unknown')
	n=0
	do k=1, kma
	   do i=1, ima
	      if (noda(i,1,k).gt.0) then
                 write(14,*) npos(noda(i,1,k),2),npos(noda(i,1,k),3)
	         n=n+1
	      endif
	   enddo
	enddo
	write(*,'(t4,a,i3,a)') '* zona_yz.plo inneh. : ',n,' koord.'
	close(14)
 
	open(unit=14,file='zona_xyt.plo',status='unknown')
	n=0
	do k=1, kma
	   do j=1, 2*(na-nb)+1
	      if (noda(ima,j,k).gt.0) then
                write(14,*) npos(noda(ima,j,k),1),npos(noda(ima,j,k),2)
	        n=n+1
	      endif
	   enddo
	enddo
	write(*,'(t4,a,i3,a)') '* zona_xyt.plo inneh. : ',n,' koord.'
	close(14)
 
	open(unit=14,file='zona_xy0.plo',status='unknown')
	n=0
	do k=1, kma
	   do j=1, jma-4
	      if (noda(1,j,k).gt.0) then
                write(14,*) npos(noda(1,j,k),1),npos(noda(1,j,k),2)
	        n=n+1
	      endif
	   enddo
	   do i=1, ima
	      if (noda(i,jma,k).gt.0) then
                write(14,*) npos(noda(i,jma,k),1),npos(noda(i,jma,k),2)
	        n=n+1
	      endif
	   enddo
	enddo
	write(*,'(t4,a,i3,a)') '* zona_xy0.plo inneh. : ',n,' koord.'
	close(14)
 
	do k=1, kma
	   if (k.lt.10) then
	      write(fil,'(a,i1,a)') 'zona_xz0',k,'.plo'
	   else
	      write(fil,'(a,i2,a)') 'zona_xz',k,'.plo'
	   endif
	   open(unit=14,file=fil,status='unknown')
	   n=0
	   do i=1, ima
	      do j=1, jma
	         if (noda(i,j,k).gt.0) then
                   write(14,*) npos(noda(i,j,k),1),npos(noda(i,j,k),3)
		   n=n+1
	         endif
	      enddo
	   enddo
	   write(*,'(t4,a,a,a,i3,a)') '* ',fil,' inneh. : ',n,' koord.'
	   close(14)
	enddo
 
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine plot_nod_zonb(nodb,nb1,nb2,nb3)
C----------------------------------------------------------C
C  Rutinen skapar en plotfil f|r noderna i Zon A.          C
C----------------------------------------------------------C
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     &         sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma,
     &         nb1,nb2,nb3,nodb(nb1,nb2,nb3),n,i,j,k
c
      character  fil*13
 
      common  /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
      common  /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
	open(unit=14,file='zonb_yz.plo',status='unknown')
	n=0
	do k=1, kma
	   do j=1, jmb
	     if (nodb(imb,j,k).gt.0) then
               write(14,*) npos(nodb(imb,j,k),2),npos(nodb(imb,j,k),3)
	       n=n+1
	     endif
	   enddo
	enddo
	write(*,'(t4,a,i3,a)') '* zonb_yz.plo inneh. : ',n,' koord.'
	close(14)
 
	do k=1, kma
	   if (k.lt.10) then
	      write(fil,'(a,i1,a)') 'zonb_xz0',k,'.plo'
	   else
	      write(fil,'(a,i2,a)') 'zonb_xz',k,'.plo'
	   endif
	   open(unit=14,file=fil,status='unknown')
	   n=0
	   do i=1, imb
	      do j=1, jmb
	         if (nodb(i,j,k).gt.0) then
                   write(14,*) npos(nodb(i,j,k),1),npos(nodb(i,j,k),3)
		   n=n+1
	         endif
	      enddo
	   enddo
	   write(*,'(t4,a,a,a,i3,a)') '* ',fil,' inneh. : ',n,' koord.'
	   close(14)
	enddo
 
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
