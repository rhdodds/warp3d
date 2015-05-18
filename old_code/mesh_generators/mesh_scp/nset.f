c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine make_nodset(uni,nset,elast, nods,ns1,ns2,ns3,
     &                     noda,na1,na2,na3, nodb,nb1,nb2,nb3 )
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
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  uni, nset(1), elast, n,
     1         mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     2         sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma
c
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &       /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
      write(uni,'(t1,a)') '*NSET,NSET=NFRONT'
      call nodset_nfront_(nset,n, nods,ns1,ns2,ns3,
     &          noda,na1,na2,na3, nodb,nb1,nb2,nb3 )
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)') '**  nfront contains ',n,' nodes'
c
      write(uni,'(t1,a)') '*NSET,NSET=NTIP'
      call nodset_ntip_(nset,n, nods,ns1,ns2,ns3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)')'** ntip contains ',n,' nodes'
c
      write(uni,'(t1,a)') '*NSET,NSET=NCRSUR'
      call nodset_ncrsur_(nset,n,elast, nods,ns1,ns2,ns3,
     &                                  noda,na1,na2,na3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)')'** ncrsur contains ',n,' nodes'
c
      call nodset_cr1_cr1jms_(uni,elast, nods,ns1,ns2,ns3)
c
      write(uni,'(t1,a)') '*NSET,NSET=NSYM'
      call nodset_nsym_(nset,n,elast, nods,ns1,ns2,ns3,
     &                                noda,na1,na2,na3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)')'** nsym contains ',n,' nodes'
c
      write(uni,'(t1,a)') '*NSET,NSET=NBAK'
      call nodset_nbak_(nset,n,noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)')'** nbak contains ',n,' nodes'
c
      write(uni,'(t1,a)') '*NSET,NSET=NUND'
      call nodset_nund_(nset,n,elast, nods,ns1,ns2,ns3,
     &              noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)')'** nund contains ',n,' nodes'
c
      write(uni,'(t1,a)') '*NSET,NSET=NSIDA'
      call nodset_nsida_(nset,n, nodb,nb1,nb2,nb3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)')'** nsida contains ',n,' nodes'
c
      write(uni,'(t1,a)') '*NSET,NSET=NTOP'
      call nodset_ntop_(nset,n,noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      call write_out_set(uni,nset,n)
      write(uni,'(t1,a,i5,a)')'** ntop contains ',n,' nodes'
c
      return
      end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nodset_nfront_(nset,n, nods,ns1,ns2,ns3,
     &                noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c-- (1) NFRONT  (noder ovanfor sprickspetsen, enbart den
c--              "oversta sprickspetsnoden" ingar ! )
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer nset(1),n, i,j,k,
     &        mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     &        sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma
 
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &       /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
        n = 0
c
c-Zone S:
c
	do j=1, jms
	   do k=1, kms
	      if (nnr(nods(1,j,k)).gt.0) then
                 n=n+1
                 nset(n) = nnr(nods(1,j,k))
	      endif
	   enddo
	enddo
c
c-Zone A:
c
	do j=1, jma
	   do i=2*(m1+mh)+2, ima
	      if (nnr(noda(i,j,1)).gt.0) then
                 n=n+1
                 nset(n) = nnr(noda(i,j,1))
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
      subroutine nodset_ntip_(nset,n, nods,ns1,ns2,ns3)
c
c-- (2) NTIP    (enbart en den "oversta" sprickspetsnoden )
c
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer nset(1),n, ns1,ns2,ns3,nods(ns1,ns2,ns3), i,j,
     &        mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     &        sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma
c
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &       /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
      n = 0
c-Zone S:
      do j=1, jms
         if (nnr(nods(1,j,1)).gt.0) then
            n=n+1
            nset(n) = nnr(nods(1,j,1))
         endif
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nodset_ncrsur_(nset,n,elast, nods,ns1,ns2,ns3,
     &                                        noda,na1,na2,na3)
c
c-- (3) NCRSUR  (noder under sprickfronten,
c--             "sprickspetsnoderna" ingar dock ej ! )
c
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     &         na1,na2,na3,noda(na1,na2,na3)
c
      integer nset(1),n, i,j,k, elast, ist,
     1        mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     2        sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma
 
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &       /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
	if (elast.eq.1) then
	   ist=(ims+1)/2+1
	else
	   ist=2
	endif
        n = 0
c
c-Zone S:
c
	do j=1, jms
	   do i=ist, ims-3
	      if (nnr(nods(i,j,1)).gt.0) then
                 n=n+1
                 nset(n) = nnr(nods(i,j,1))
	      endif
	   enddo
	   do k=3, kms
	      if (nnr(nods(ims,j,k)).gt.0) then
                 n=n+1
                 nset(n) = nnr(nods(ims,j,k))
	      endif
	   enddo
	enddo
c-Zone A:
	do j=1, jma-4
	   do i=1, 2*(m1-mh)
	      if (nnr(noda(i,j,1)).gt.0) then
                 n=n+1
                 nset(n) = nnr(noda(i,j,1))
	      endif
	   enddo
	enddo
	do j=jma-3, jma
	   do i=6, 2*(m1-mh)
	      if (nnr(noda(i,j,1)).gt.0) then
                 n=n+1
                 nset(n) = nnr(noda(i,j,1))
	      endif
	   enddo
	enddo
c
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nodset_cr1_cr1jms_(uni,elast, nods,ns1,ns2,ns3)
c
c-- (4) CR01, ... ,CRNN (NN=JMA) (noder for resp koordnat)
c
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer uni, ns1,ns2,ns3, nods(ns1,ns2,ns3),
     &        mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     &        sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma,
     &        nv(12),i,j,k,l, elast, ist, jj
c
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &       /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
	if (elast.eq.1) then
	   ist=(ims+1)/2+1
	else
	   ist=2
	endif
C.......Skriv ut noderna pa filen tillhorande nodsetet !
C . . . ! ZON S !
        jj = 0
	do j=1, jms
           if (nnr(nods(1,j,1)).gt.0) then
              jj = jj + 1
              l=0
              if (jj.lt.10) then
                 write(uni,'(a,i1)') '*nset,nset=cr0',jj
              else
                 write(uni,'(a,i2)') '*nset,nset=cr',jj
              endif
C
              l=l+1
              nv(l)=nnr(nods(1,j,1))
              do i=ist, ims-2
                 if (nnr(nods(i,j,1)).gt.0) then
                    if (l.lt.12) then
                       l=l+1
                       nv(l)=nnr(nods(i,j,1))
                    else
                       call write_in_nod_or_el_set_(uni,nv,12)
                       nv(1)=nnr(nods(i,j,1))
                       l=1
                    endif
                 endif
              enddo
C
              do k=2, 3
                 do i=1, ims-2
                    if (nnr(nods(i,j,k)).gt.0) then
                       if (l.lt.12) then
                          l=l+1
                          nv(l)=nnr(nods(i,j,k))
                       else
                          call write_in_nod_or_el_set_(uni,nv,12)
                          nv(1)=nnr(nods(i,j,k))
                          l=1
                       endif
                    endif
                 enddo
              enddo
              call write_in_nod_or_el_set_(uni,nv,l)
           endif
        enddo
c
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nodset_nsym_(nset,n,elast, nods,ns1,ns2,ns3,
     &                                      noda,na1,na2,na3 )
c
c-- (5) NSYM    (noder pa symmetrisnittet )
c
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3)
c
      integer nset(1),n,i,k, kstart,kstep, elast,ist,
     &        mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     &        sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma
 
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &       /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma,
     &       /nblock/ kstart,kstep
 
      if (elast.eq.1) then
         ist=(ims+1)/2+1
      else
         ist=2
      endif
      n = 0
c-Zone S:
      n = n + 1
      nset(n) = nnr(nods(1,1,1))
      do i=ist, ims-2
         if (nnr(nods(i,1,1)).gt.0) then
            n = n + 1
            nset(n) = nnr(nods(i,1,1))
         endif
      enddo
c
      do k=2, kms
         do i=1, ims
            if ((k.gt.3).or.(i.le.(ims-2))) then
               if (nnr(nods(i,1,k)).gt.0) then
                  n = n + 1
                  nset(n) = nnr(nods(i,1,k))
               endif
            endif
         enddo
      enddo
c-Zone A:
      do k=1, kma
         do i=1, ima
            if ((nnr(noda(i,1,k)).gt.0).and.
     &          (noda(i,1,k).gt.kstart) ) then
               n = n + 1
               nset(n) = nnr(noda(i,1,k))
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
      include 'mesh3d_scp_common_nod.f'
c
      integer  na1,na2,na3,noda(na1,na2,na3),
     &         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer nset(1),n, i,j,
     &        mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     &        sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma
c
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &       /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
        n = 0
c-Zone A:
	do j=1, jma-4
	   do i=1, ima
	      if (nnr(noda(i,j,kma)).gt.0) then
                 n = n + 1
                 nset(n) = nnr(noda(i,j,kma))
	      endif
	   enddo
	enddo
	do j=jma-3, jma
	   do i=6, ima
	      if (nnr(noda(i,j,kma)).gt.0) then
                 n = n + 1
                 nset(n) = nnr(noda(i,j,kma))
	      endif
	   enddo
	enddo
c-Zone B:
	do j=1, jmb
	   do i=2, imb
	      if (nnr(nodb(i,j,kma)).gt.0) then
                 n = n + 1
                 nset(n) = nnr(nodb(i,j,kma))
	      endif
	   enddo
	enddo
c
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nodset_nund_(nset,n,elast, nods,ns1,ns2,ns3,
     &                 noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c---- (7) NUND    (noder pa modellens undersida )
c
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer nset(1),n, i,j,k, kstart,kstep, elast,ist,
     1        mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     2        sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma
c
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &       /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma,
     &       /nblock/ kstart,kstep
c
	if (elast.eq.1) then
	   ist=(ims+1)/2+1
	else
	   ist=2
	endif
        n = 0
c
c-Zone S:
c
        n = n + 1
        nset(n) = nnr(nods(1,jms,1))
	do i=ist, ims
	   if (nnr(nods(i,jms,1)).gt.0) then
              n = n + 1
              nset(n) = nnr(nods(i,jms,1))
	   endif
	enddo
c
	do k=2, kms
	   do i=1, ims
	      if ( (k.gt.3).or.(i.le.(ims-2)) ) then
                 if (nnr(nods(i,jms,k)).gt.0) then
                    n = n + 1
                    nset(n) = nnr(nods(i,jms,k))
                 endif
              endif
           enddo
        enddo
c-Zone A:
	do k=1, kma
	   do j=1, jma-4
	      if ( (nnr(noda(1,j,k)).gt.0).and.
     &             (noda(1,j,k).gt.kstart) ) then
                 n = n + 1
                 nset(n) = nnr(noda(1,j,k))
	      endif
	   enddo
	   do i=6, ima
	      if ( (nnr(noda(i,jma,k)).gt.0).and.
     &             (noda(i,jma,k).gt.kstart) ) then
                 n = n + 1
                 nset(n) = nnr(noda(i,jma,k))
	      endif
	   enddo
	enddo
c-Zone B:
	do k=1, kma
	   do i=2, imb
	      if (nnr(nodb(i,jmb,k)).gt.0) then
                 n = n + 1
                 nset(n) = nnr(nodb(i,jmb,k))
	      endif
	   enddo
	enddo
c
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nodset_nsida_(nset,n, nodb,nb1,nb2,nb3)
c
c-- (8) NSIDA   (noder pa sidan mittemot symmetrisnittet )
c
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer nset(1),n, j,k,
     &        mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     &        sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma
 
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &       /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
        n = 0
c-Zone B:
	do k=1, kma
	   do j=1, jmb
	      if (nnr(nodb(imb,j,k)).gt.0) then
                 n = n + 1
                 nset(n) = nnr(nodb(imb,j,k))
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
      include 'mesh3d_scp_common_nod.f'
c
      integer  na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer nset(1),n, i,j,k,
     1        mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     2        sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma
c
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
      common /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
	n = 0
c-Zone A:
	do k=1, kma
	   do j=1, jma-jmb+1
	      if (nnr(noda(ima,j,k)).gt.0) then
                 n = n + 1
                 nset(n) = nnr(noda(ima,j,k))
	      endif
	   enddo
	enddo
c-Zone B:
	do k=1, kma
	   do i=2, imb
	      if (nnr(nodb(i,1,k)).gt.0) then
                 n = n + 1
                 nset(n) = nnr(nodb(i,1,k))
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
      integer  uni,set(1),n, nv(12),i,j
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
	integer uni,l,nv(1),i
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
	subroutine make_elemset(uni,estk,lt,elsym,elund,elsua)
C--- Rutinen laser de redan skapade elementseten ifran UNI1
C--- och skriver ut dem pa enhet UNI  , dvs input filen.
C--- Dessutom skapas  elementset for elementlagren (LT st.)
c
      include 'mesh3d_scp_common_eln.f'
      integer  eln_d(0:iem)
      integer  elsym(iem),elund(iem),elsua(iem)
c
      integer  uni,estk(0:100,5),lt,l
c
	write(uni,'(t1,a,i1,a)')'*elset,elset=layercr,generate'
	write(uni,'(t3,i5,a,i5)') estk(0,1), ',' ,estk(0,5)
	do l=1, lt
	  if (l.lt.10) then
	    write(uni,'(t1,a,i1,a)')'*elset,elset=layer0',l,',generate'
	  else
	    write(uni,'(t1,a,i2,a)')'*elset,elset=layer',l,',generate'
	  endif
	  write(uni,'(t3,i5,a,i5)') estk(l,1), ',' ,estk(l,4)
	enddo
C.......Elementset pa symmetriplan resp undersida av modellen.
C       Definiera {ven ytnormalen f|r noder i "ytorna" av dessa elem.set
C       som ligger p} planen vinkelr{t mot sprickplanet.
C.......Detta skrivs ocks} p} filen : NODELEM.LIS
	call elemset_(uni,elsym,iem,'*elset,elset=elemsym')
	call elemset_(uni,elund,iem,'*elset,elset=elemund')
	call elemset_(uni,elsua,iem,'*elset,elset=elsurfa')
	call define_normal_(uni,'elemsym','nsym',
     &                         'elemund','nund')
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
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
