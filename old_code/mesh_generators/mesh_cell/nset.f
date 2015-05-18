c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine make_nodset(uni,nset, nodc,nc1,nc2,nc3,
     &       nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3 )
c
C  * Rutinen skapar nodset for modellen, dessa ar :            *
C  *  (1) NFRONT  (noder ovanfor sprickspetsen, enbart den     *
C  *               "oversta sprickspetsnoden" ingar ! )        *
C  *  (2) NTIP    (enbart en den oversta sprickspetsnoden )    *
C  *  (3) NCRSUR  (noder under sprickfronten,                  *
C  *               "sprickspetsnoderna" ingar dock ej ! )      *
C  *  (4) CR01, ... ,CRNN (NN=JMA) (noder for resp koordnat)   *
C  *  (5) NSYM    (noder pa symmetrisnittet )                  *
C  *  (6) NBAK    (noder pa modellens baksida)                 *
C  *  (7) NUND    (noder pa modellens undersida )              *
C  *  (8) NSIDA   (noder pa sidan mittemot symmetrisnittet )   *
C  *  (9) NTOP    (noder pa modellens ovansida )               *
c
      implicit none
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  uni, nset(*),  n
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
      write(uni,'(t1,a)') '*nset,nset=nfront'
      call nodset_nfront_(nset,n, nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3,
     &                            noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)') '**  nfront contains ',n,' nodes'
c
      write(uni,'(t1,a)') '*nset,nset=ntip'
      call nodset_ntip_(nset,n, nodc,nc1,nc2,nc3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)')'** ntip contains ',n,' nodes'
c
      write(uni,'(t1,a)') '*nset,nset=ncrsur'
      call nodset_ncrsur_(nset,n, nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3, 
     &                            noda,na1,na2,na3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)')'** ncrsur contains ',n,' nodes'
c
      call nodset_cr1_cr1jms_(uni, nodc,nc1,nc2,nc3)
c
      write(uni,'(t1,a)') '*nset,nset=nsym'
      call nodset_nsym_(nset,n, nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3,
     &                          noda,na1,na2,na3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)')'** nsym contains ',n,' nodes'
c
      write(uni,'(t1,a)') '*nset,nset=nbak'
      call nodset_nbak_(nset,n,noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)')'** nbak contains ',n,' nodes'
c
      write(uni,'(t1,a)') '*nset,nset=nund'
      call nodset_nund_(nset,n, nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3,
     &                          noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)')'** nund contains ',n,' nodes'
c
      write(uni,'(t1,a)') '*nset,nset=nsida'
      call nodset_nsida_(nset,n, nodb,nb1,nb2,nb3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)')'** nsida contains ',n,' nodes'
c
      write(uni,'(t1,a)') '*nset,nset=ntop'
      call nodset_ntop_(nset,n,noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)')'** ntop contains ',n,' nodes'
c
      return
      end
c
C----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine nodset_nfront_(nset,n, nodc,nc1,nc2,nc3,
     &     nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c--(1) NFRONT  (nodes ahead of the crack front at the plane Y=0)
c
      implicit none
c
      include 'common_nod.f'
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer nset(*),n, i,j,k, ia2
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
      n = 0
c
c-Zone C:
c
      k = 1
      do j=1, jmc
         do i=5, imc
            if (nnr(nodc(i,j,k)).gt.0) then
               n=n+1
               nset(n) = nnr(nodc(i,j,k))
            endif
         enddo
      enddo
c
c-Zone S:
c
      i = 1
      do j=1, jms
         do k=2, kms
	      if (nnr(nods(i,j,k)).gt.0) then
                 n=n+1
                 nset(n) = nnr(nods(i,j,k))
	      endif
	   enddo
	enddo
c
c-Zone A:
c
      ia2 = 2*(m1+mh+mh)+1
      k = 1
      do j=1, jma
         do i=ia2+1, ima
            if (nnr(noda(i,j,k)).gt.0) then
               n=n+1
               nset(n) = nnr(noda(i,j,k))
            endif
         enddo
      enddo
c
c-Zone B:
c
      do j=1, jmb
         do i=2, imb
            if (nnr(nodb(i,j,1)).gt.0) then
               n=n+1
               nset(n) = nnr(nodb(i,j,1))
            endif
         enddo
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nodset_ntip_(nset,n, nodc,nc1,nc2,nc3)
c
c--(2) NTIP   (initial crack tip nodes)
c
      implicit none
c
      include 'common_nod.f'
c
      integer nset(*),n, nc1,nc2,nc3,nodc(nc1,nc2,nc3), i,j,k
c
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
      n = 0
c-Zone C:
      i = 5
      k = 1
      do j=1, jmc
         if (nnr(nodc(i,j,k)).gt.0) then
            n=n+1
            nset(n) = nnr(nodc(i,j,k))
         endif
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nodset_ncrsur_(nset,n, nodc,nc1,nc2,nc3,
     &                nods,ns1,ns2,ns3, noda,na1,na2,na3)
c
c--(3) NCRSUR  (nodes located on the initial crack surface)
c
      implicit none
c
      include 'common_nod.f'
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3)
c
      integer nset(*),n, i,j,k, ia1
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
        n = 0
c
c-Zone C:
c
      k = 3
      do j=1, jmc
         do i=1, 5
            if (nnr(nodc(i,j,k)).gt.0) then
               n=n+1
               nset(n) = nnr(nodc(i,j,k))
            endif
         enddo
      enddo
c
c-Zone S:
c
      i = ims
      do j=1, jms
         do k=2, kms
            if (nnr(nods(i,j,k)).gt.0) then
               n=n+1
               nset(n) = nnr(nods(i,j,k))
            endif
         enddo
      enddo
c
c-Zone A:
c
      ia1 = 2*m1 + 1
      k = 1
      do j=1, jma-4
         do i=1, ia1-1
            if (nnr(noda(i,j,k)).gt.0) then
               n=n+1
               nset(n) = nnr(noda(i,j,k))
            endif
         enddo
      enddo
      do j=jma-3, jma
         do i=6, ia1-1
            if (nnr(noda(i,j,k)).gt.0) then
               n=n+1
               nset(n) = nnr(noda(i,j,k))
            endif
         enddo
      enddo
c
      return
      end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nodset_cr1_cr1jms_(uni, nodc,nc1,nc2,nc3)
c
c--(4) cr01, ... ,crnn (nn=jmc) (nodes on the crack plane)
c
      implicit none
c
      include 'common_nod.f'
c
      integer uni, nc1,nc2,nc3, nodc(nc1,nc2,nc3),
     &        nv(12),i,j,k,l
c
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
c-Zone C:
c
      k = 1
      do j=1, jmc
         l=0
         if (j.lt.10) then
            write(uni,'(a,i1)') '*nset,nset=cr0',j
         else
            write(uni,'(a,i2)') '*nset,nset=cr',j
         endif
c
         nv(l)=nnr(nodc(1,j,k))
         do i=5, imc-6
            if (nnr(nodc(i,j,k)).gt.0) then
               if (l.lt.12) then
                  l=l+1
                  nv(l)=nnr(nodc(i,j,1))
               else
                  call write_in_nod_or_el_set_(uni,nv,12)
                  nv(1)=nnr(nodc(i,j,1))
                 l=1
               endif
            endif
         enddo
         call write_in_nod_or_el_set_(uni,nv,l)
      enddo
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nodset_nsym_(nset,n,nodc,nc1,nc2,nc3,
     &             nods,ns1,ns2,ns3, noda,na1,na2,na3)
c
c--(5) NSYM  (nodes on the symmetry plane, X=0)
c
      implicit none
c
      include 'common_nod.f'
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3)
c
      integer nset(*),n,i,j,k
c
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer         kstart,kstep
      common /nblock/ kstart,kstep
c
      n = 0
c
c-Zone C:
c
      j = 1
      do k=1, kmc
         do i=1, imc
            if (nnr(nodc(i,j,k)).gt.0) then
               n = n + 1
               nset(n) = nnr(nodc(i,j,k))
            endif
         enddo
      enddo
c
c-Zone S:
c
      j = 1
      do k=2, kms
         do i=1, ims
            if (nnr(nods(i,j,k)).gt.0) then
               n = n + 1
               nset(n) = nnr(nods(i,j,k))
            endif
         enddo
      enddo
c-Zone A:
      j = 1
      do k=1, kma
         do i=1, ima
            if ((nnr(noda(i,j,k)).gt.0).and.
     &          (noda(i,j,k).gt.kstart) ) then
               n = n + 1
               nset(n) = nnr(noda(i,j,k))
            endif
         enddo
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nodset_nbak_(nset,n,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
c
c-- (6) NBAK    (noder pa modellens baksida)
c
      implicit none
c
      include 'common_nod.f'
c
      integer  na1,na2,na3,noda(na1,na2,na3),
     &         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer nset(*),n, i,j,k
c
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
        n = 0
c-Zone A:
        k = kma
	do j=1, jma-4
	   do i=1, ima
	      if (nnr(noda(i,j,k)).gt.0) then
                 n = n + 1
                 nset(n) = nnr(noda(i,j,k))
	      endif
	   enddo
	enddo
	do j=jma-3, jma
	   do i=6, ima
	      if (nnr(noda(i,j,k)).gt.0) then
                 n = n + 1
                 nset(n) = nnr(noda(i,j,k))
	      endif
	   enddo
	enddo
c-Zone B:
        k = kma
	do j=1, jmb
	   do i=2, imb
	      if (nnr(nodb(i,j,k)).gt.0) then
                 n = n + 1
                 nset(n) = nnr(nodb(i,j,k))
	      endif
	   enddo
	enddo
c
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nodset_nund_(nset,n, nodc,nc1,nc2,nc3,
     &      nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c---- (7) NUND    (nodes on the plane Z=0)
c
      implicit none
c
      include 'common_nod.f'
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer nset(*),n, i,j,k
c
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer         kstart,kstep
      common /nblock/ kstart,kstep
c
        n = 0
c
c-Zone C:
c
      j = jmc
      do k=1, kmc
         do i=1, imc
            if (nnr(nodc(i,j,k)).gt.0) then
               n = n + 1
               nset(n) = nnr(nodc(i,j,k))
            endif
         enddo
      enddo
c
c-Zone S:
c
      j = jms
      do k=2, kms
         do i=1, ims
            if (nnr(nods(i,j,k)).gt.0) then
               n = n + 1
               nset(n) = nnr(nods(i,j,k))
            endif
         enddo
      enddo
c
c-Zone A:
c
      do k=1, kma
         i = 1
         do j=1, jma-4
            if ((nnr(noda(i,j,k)).gt.0).and.
     &          (noda(i,j,k).gt.kstart)) then
               n = n + 1
               nset(n) = nnr(noda(i,j,k))
            endif
         enddo
         j = jma
         do i=6, ima
            if ((nnr(noda(i,j,k)).gt.0).and.
     &          (noda(i,j,k).gt.kstart) ) then
               n = n + 1
               nset(n) = nnr(noda(i,j,k))
            endif
         enddo
      enddo
c
c-Zone B:
c
      j = jmb
      do k=1, kma
         do i=2, imb
            if (nnr(nodb(i,j,k)).gt.0) then
               n = n + 1
               nset(n) = nnr(nodb(i,j,k))
            endif
         enddo
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine nodset_nsida_(nset,n, nodb,nb1,nb2,nb3)
c
c-- (8) NSIDA   (noder pa sidan mittemot symmetrisnittet )
c
      implicit none
c
      include 'common_nod.f'
c
      integer  nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      integer  nset(*),n, i,j,k
c
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
      n = 0
c-Zone B:
c
      i = imb
c
      do k=1, kma
         do j=1, jmb
            if (nnr(nodb(i,j,k)).gt.0) then
               n = n + 1
               nset(n) = nnr(nodb(i,j,k))
            endif
         enddo
      enddo
c
      return
      end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nodset_ntop_(nset,n,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
c
c-- (9) NTOP    (noder pa modellens ovansida )
c
      implicit none
c
      include 'common_nod.f'
c
      integer  na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  nset(*),n, i,j,k
c
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
      n = 0
c
c-Zone A:
c
      i = ima
c
      do k=1, kma
         do j=1, jma-jmb+1
            if (nnr(noda(i,j,k)).gt.0) then
               n = n + 1
               nset(n) = nnr(noda(i,j,k))
            endif
         enddo
      enddo
c
c-Zone B:
c
      j = 1
c
      do k=1, kma
         do i=2, imb
            if (nnr(nodb(i,j,k)).gt.0) then
               n = n + 1
               nset(n) = nnr(nodb(i,j,k))
            endif
         enddo
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine write_out_set(uni,set,n)
      implicit none
      integer  uni,set(*),n, nv(12),i,j
C
      j = 0
      do i=1, n
         j = j + 1
         nv(j) = set(i)
         if (j.eq.12) then
            call write_in_nod_or_el_set_(uni,nv,j)
            j = 0
         endif
      enddo
      if (j.gt.0) then
         call write_in_nod_or_el_set_(uni,nv,j)
      endif
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine write_in_nod_or_el_set_(uni,nv,l)
c --- Rutinen skriver ut noder p} enheten UNIT=UNI
	implicit none
	integer uni,l,nv(*),i
	if (l.eq.1) then
	   write(uni,'(t1,i5)') nv(1)
	else
	   write(uni,'(t1,i5,11(a,i5))') nv(1),(',',nv(i), i=2,l)
	endif
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine elemset_(uni,el,iem,text)
C----- Rutinen skapar ett elementset av ytelement.
	implicit none
	integer uni,i,l,iem,el(iem),ev(12)
	character text*20
C.......Elementen skrivs ut p} filen.
	write(uni,'(t1,a)') text
	l=0
	do i=1, iem
	   if (el(i).gt.0) then
	      if (l.lt.12) then
                 l=l+1
	         ev(l)=i
	      else
	         call write_in_nod_or_el_set_(uni,ev,12)
	         ev(1)=i
		 l=1
	      endif
	   endif
	enddo
	call write_in_nod_or_el_set_(uni,ev,l)
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine define_normal_(uni,elset1,nset1,elset2,nset2)
c ---- Rutinen definierar normalriktningar f|r noder p}
c ---- plan som genomsk{rs av sprickspetsen.
	integer uni
	real*8 x1(3),x2(3)
	character elset1*7, elset2*7, nset1*4, nset2*4,text*13
	x1(1)=-1.
	x1(2)=0.
	x1(3)=0.
	x2(1)=0.
	x2(2)=0.
	x2(3)=-1.
	write(uni,'(t1,a)') '*normal'
	text=elset1//','//nset1//','
        write(uni,101) text, x1(1),',',x1(2),',',x1(3)
	text=elset2//','//nset2//','
        write(uni,101) text ,x2(1),',',x2(2),',',x2(3)
101	format(t3,a,f5.2,a,f5.2,a,f5.2)
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
