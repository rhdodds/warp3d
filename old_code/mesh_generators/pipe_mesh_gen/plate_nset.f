c----67--1---------2---------3---------4---------5---------6---------712
c
c    File: plate_nset.f   latest modification April 23, 1998
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine nset_nfront(nset,n, nods,ns1,ns2,ns3,
     &                       noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c NFRONT: nodes ahead of the crack on the uncracked ligament on
c         the plane of the crack.
c
      implicit none
c
      include 'plate_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer nset(*),n, i,j,k, ia2
c
      integer      mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
      n = 0
c
c Zone S
c
      i = 1
      do j=1, jms
         do k=1, kms
            if (nnr(nods(i,j,k)).gt.0) then
                 n=n+1
                 nset(n) = nnr(nods(i,j,k))
            endif
         enddo
      enddo
c
c Zone A
c
      ia2 = 2*(m1+mh+mh)+1
      k = 1
      do j=1, jma
         do i=ia2, ima
            if (nnr(noda(i,j,k)).gt.0) then
               n=n+1
               nset(n) = nnr(noda(i,j,k))
            endif
         enddo
      enddo
c
c Zone B
c
      do j=1, jmb
         do i=1, imb
            if (nnr(nodb(i,j,1)).gt.0) then
               n=n+1
               nset(n) = nnr(nodb(i,j,1))
            endif
         enddo
      enddo
c
      call nset_unique(nset,n)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nset_ntip(nset,n, nods,ns1,ns2,ns3)
c
c NTIP: nodes at the crack tip along the crack front
c
      implicit none
c
      include 'plate_common_nod.f'
c
      integer nset(*),n, ns1,ns2,ns3,nods(ns1,ns2,ns3), i,j,k
c
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
      n = 0
c
c Zone S
c
      i = 1
      k = 1
      do j=1, jms
         if (nnr(nods(i,j,k)).gt.0) then
            n=n+1
            nset(n) = nnr(nods(i,j,k))
         endif
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine ncrd_ntip(nset,cset,n, nods,ns1,ns2,ns3)
c
c NTIP: nodes at the crack tip along the crack front
c
      implicit none
      include 'plate_common_nod.f'
c
      integer  nset(*),n, ns1,ns2,ns3,nods(ns1,ns2,ns3), i,j,k,np
      double precision cset(3,*) 
c
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
c Zone S
c
      n = 0
      i = 1
      k = 1
      do j=1, jms
         np = nods(i,j,k)
         if (nnr(np).gt.0) then
            n=n+1
            nset(n) = nnr(np)
            cset(1,n) = npos(np,1)
            cset(2,n) = npos(np,2)
            cset(3,n) = npos(np,3)
         endif
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine ncrd_ntipset(nset,cset,nct,ns, nods,ns1,ns2,ns3)
c
c NTIP: nodes at the crack tip along the crack front
c
      implicit none
      include 'plate_common_nod.f'
c
      integer  nset(*),nct,ns(2,*),ns1,ns2,ns3,nods(ns1,ns2,ns3),
     &         i,j,k,n,ntmp,np
      double precision cset(3,*) 
c
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
c Zone S
c
      n = 0
      k = 1
      do j=1, jms
         ntmp = 0
         do i=1, ims
            np = nods(i,j,k)
            if (nnr(np).gt.0) then
               n = n + 1
               ntmp = ntmp + 1
               nset(n) = nnr(np)
               cset(1,n) = npos(np,1)
               cset(2,n) = npos(np,2)
               cset(3,n) = npos(np,3)
            endif
         enddo
         ns(1,j) = n + 1 - ntmp
         ns(2,j) = n
      enddo
c
c Remove gaps in ns(i,j)
c
      ntmp = 0
      do j=1, jms
         if (ns(1,j).lt.ns(2,j)) then
            ntmp = ntmp + 1
            ns(1,ntmp) = ns(1,j)
            ns(2,ntmp) = ns(2,j)
         endif
      enddo
      do j=ntmp+1, jms
         ns(1,j) = 0
         ns(2,j) = 0
      enddo
      nct = ntmp
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nset_ntip_coincident(nset,n, nods,ns1,ns2,ns3)
c
c NTIPCOIN: nodes coincident with the single nodes at the crack tip
c           along the crack front. The single crack tip nodes are
c           not included.
c
      implicit none
c
      include 'plate_common_nod.f'
c
      integer nset(*),n, ns1,ns2,ns3,nods(ns1,ns2,ns3), i,j,k
c
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
      n = 0
c
c Zone S
c
      k = 1
      do j=1, jms
         do i=2, ims
            if ( nnr(nods(i,j,k)) .gt. 0 ) then
               n = n + 1
               nset(n) = nnr(nods(i,j,k))
            endif
         enddo
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nset_ntipall(nset,n,nods,ns1,ns2,ns3)
c
c NTIP: all nodes at in the focused mesh at the crack tip all
c       along the crack front
c
      implicit none
c
      include 'plate_common_nod.f'
c
      integer nset(*),n, ns1,ns2,ns3,nods(ns1,ns2,ns3), i,j,k
c
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
      n = 0
c
c Zone S
c
      k = 1
      do j=1, jms
         do i=1, ims
            if (nnr(nods(i,j,k)).gt.0) then
               n=n+1
               nset(n) = nnr(nods(i,j,k))
            endif
         enddo
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nset_cr1_cr1jms_(uni, nods,ns1,ns2,ns3,nc,
     &                            analysis,keyhole)
c
c Define crack tip node sets to be used in def. of domain integration
c        cr01, ... ,crnn (nn=jmc)
c
      implicit none
      include 'plate_common_nod.f'
c
      integer uni, ns1,ns2,ns3, nods(ns1,ns2,ns3), analysis,keyhole,
     &        ntip,nc,nv(12),i,j,k,l,iend
c
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
c Zone S
c
      if ( (analysis.eq.1).and.(keyhole.eq.0) ) then
         iend = 1
      else
         iend = ims
      endif
      k = 1
      nc = 0
      do j=1, jms
         ntip = nnr(nods(1,j,k))
         if (ntip.gt.0) then
            nc = nc + 1
            if (nc.lt.10) then
               write(uni,'(a,i1)') '*NSET,NSET=CR0',nc
            else
               write(uni,'(a,i2)') '*NSET,NSET=CR',nc
            endif
            l = 1
            nv(l) = ntip
            do i=2, iend
               if (nnr(nods(i,j,k)).gt.0) then
                  if (l.lt.12) then
                     l=l+1
                     nv(l)=nnr(nods(i,j,1))
                  else
                     call write_in_nod_or_el_set(uni,nv,12)
                     nv(1)=nnr(nods(i,j,1))
                     l=1
                  endif
               endif
            enddo
            call write_in_nod_or_el_set(uni,nv,l)
         endif
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nset_ncrsur(nset,n,nods,ns1,ns2,ns3,noda,na1,na2,na3)
c
c NCRSUR: nodes on the crack face excluding crack tip nodes and nodes
c         coincident with the crack tip nodes with the exception of
c         the "180 degree" node in the focused crack tip mesh
c 
      implicit none
c
      include 'plate_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     &         na1,na2,na3,noda(na1,na2,na3)
c
      integer nset(*),n, i,j,k, ia1
c
      integer      mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
        n = 0
c
c Zone S
c
      i = ims
      do j=1, jms
         do k=1, kms
            if (nnr(nods(i,j,k)).gt.0) then
               n=n+1
               nset(n) = nnr(nods(i,j,k))
            endif
         enddo
      enddo
c
c Zone A
c
      ia1 = 2*m1 + 1
      k = 1
      do j=1, jma
         do i=1, ia1
            if ( nnr(noda(i,j,k)) .gt. 0 ) then
               n=n+1
               nset(n) = nnr(noda(i,j,k))
            endif
         enddo
      enddo
c
      call nset_unique(nset,n)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nset_nasym(nset,n,nods,ns1,ns2,ns3, noda,na1,na2,na3)
c
c NASYM: nodes on the symmetry plane, X=0
c
      implicit none
c
      include 'plate_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     &         na1,na2,na3,noda(na1,na2,na3)
c
      integer nset(*),n,i,j,k
c
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
      n = 0
c
c Zone S
c
      j = 1
      do k=1, kms
         do i=1, ims
            if (nnr(nods(i,j,k)).gt.0) then
               n = n + 1
               nset(n) = nnr(nods(i,j,k))
            endif
         enddo
      enddo
c
c Zone A
c
      j = 1
      do k=1, kma
         do i=1, ima
            if ( nnr(noda(i,j,k)) .gt. 0) then
               n = n + 1
               nset(n) = nnr(noda(i,j,k))
            endif
         enddo
      enddo
c
      call nset_unique(nset,n)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine nset_ncsym(nset,n,nods,ns1,ns2,ns3,noda,na1,na2,na3,
     &                             nodb,nb1,nb2,nb3)
c
c NCSYM: nodes on the plane Z = 0 (adjacent and perpendicular to the
c                                   plane of the crack)
c
      implicit none
c
      include 'plate_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer nset(*),n, i,j,k
c
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
      integer         kstart,kstep
      common /nblock/ kstart,kstep
c
        n = 0
c
c Zone S
c
      j = jms
      do k=1, kms
         do i=1, ims
            if (nnr(nods(i,j,k)).gt.0) then
               n = n + 1
               nset(n) = nnr(nods(i,j,k))
            endif
         enddo
      enddo
c
c Zone A
c
      do k=1, kma
         i = 1
c         do j=1, jma-4
         do j=1, jma
            if ( nnr(noda(i,j,k)) .gt. 0 ) then
               n = n + 1
               nset(n) = nnr(noda(i,j,k))
            endif
         enddo
         j = jma
c         do i=6, ima
         do i=1, ima
            if ( nnr(noda(i,j,k)) .gt. 0 ) then
               n = n + 1
               nset(n) = nnr(noda(i,j,k))
            endif
         enddo
      enddo
c
c Zone B
c
      j = jmb
      do k=1, kma
         do i=1, imb
            if ( nnr(nodb(i,j,k)) .gt. 0) then
               n = n + 1
               nset(n) = nnr(nodb(i,j,k))
            endif
         enddo
      enddo
c
      call nset_unique(nset,n)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nset_nrem(nset,n,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
c
c NREM: nodes on the remote plane opposite the plane of the crack
c       ( Y = L)
c
      implicit none
c
      include 'plate_common_nod.f'
c
      integer  na1,na2,na3,noda(na1,na2,na3),
     &         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer nset(*),n, i,j,k
c
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
        n = 0
c
c Zone A
c
      k = kma
      do j=1, jma
         do i=1, ima
            if ( nnr(noda(i,j,k)) .gt. 0 ) then
               n = n + 1
               nset(n) = nnr(noda(i,j,k))
            endif
         enddo
      enddo
c
c Zone B
c
      k = kma
      do j=1, jmb
         do i=1, imb
            if ( nnr(nodb(i,j,k)) .gt. 0 ) then
                 n = n + 1
                 nset(n) = nnr(nodb(i,j,k))
            endif
         enddo
      enddo
c
      call nset_unique(nset,n)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine ncrd_nrem(nset,cset,n, noda,na1,na2,na3,
     &                                  nodb,nb1,nb2,nb3)
c
c NREM: nodes on the remote plane opposite the plane of the crack
c       ( Y = L)
c
      implicit none
c
      include 'plate_common_nod.f'
c
      integer  na1,na2,na3,noda(na1,na2,na3),
     &         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      integer nset(*),n, i,j,k,np
      double precision cset(3,*) 
c
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
        n = 0
c
c Zone A
c
      k = kma
      do j=1, jma
         do i=1, ima
            np = noda(i,j,k)
            if ( nnr(np) .gt. 0 ) then
               n = n + 1
               nset(n)   = nnr(np)
               cset(1,n) = npos(np,1)
               cset(2,n) = npos(np,2)
               cset(3,n) = npos(np,3)
            endif
         enddo
      enddo
c
c Zone B
c
      k = kma
      do j=1, jmb
         do i=1, imb
            np = nodb(i,j,k)
            if ( nnr(np) .gt. 0 ) then
                 n = n + 1
                 nset(n)   = nnr(np)
                 cset(1,n) = npos(np,1)
                 cset(2,n) = npos(np,2)
                 cset(3,n) = npos(np,3)
            endif
         enddo
      enddo
c
      call cset_unique(nset,cset,n)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nset_nsid(nset,n, nodb,nb1,nb2,nb3)
c
c NSID: nodes on the plane X = W (adjacent and perpendicular to the
c                                  plane of the crack)
c
      implicit none
c
      include 'plate_common_nod.f'
c
      integer  nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      integer  nset(*),n, i,j,k
c
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
      n = 0
c
c Zone B
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
      call nset_unique(nset,n)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine ncrd_nsid(nset,cset,n, nodb,nb1,nb2,nb3)
c
c NSID: nodes on the plane X = W (adjacent and perpendicular to the
c                                  plane of the crack)
c
      implicit none
c
      include 'plate_common_nod.f'
c
      integer  nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      integer  nset(*),n, i,j,k,np
      double precision cset(3,*) 
c
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
      n = 0
c
c Zone B
c
      i = imb
c
      do k=1, kma
         do j=1, jmb
            np = nodb(i,j,k)
            if ( nnr(np) .gt. 0 ) then
               n = n + 1
               nset(n)   = nnr(np)
               cset(1,n) = npos(np,1)
               cset(2,n) = npos(np,2)
               cset(3,n) = npos(np,3)
            endif
         enddo
      enddo
c
      call cset_unique(nset,cset,n)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nset_ntop(nset,n,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
c
c NTOP: nodes on the plane Z = T (adjacent and perpendicular to the
c                                 plane of the crack)
c
      implicit none
c
      include 'plate_common_nod.f'
c
      integer  na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  nset(*),n, i,j,k
c
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
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
         do i=1, imb
            if (nnr(nodb(i,j,k)).gt.0) then
               n = n + 1
               nset(n) = nnr(nodb(i,j,k))
            endif
         enddo
      enddo
c
      call nset_unique(nset,n)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nset_unique(nset,n)
c
c Removes node numbers if given more than once in the same node set
c
      implicit none
      integer nset(*),n, i,j,m
      logical unique
c
      m = 0
      do i=1, n
         unique = .true.
         do j=1, m
            if ( nset(i) .eq. nset(j) ) unique = .false.
         enddo
         if (unique) then
            m = m + 1
            nset(m) = nset(i)
         endif
      enddo
      n = m
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine cset_unique(nset,cset,n)
c
c Removes node numbers if given more than once in the same node set
c
      implicit none
      integer nset(*),n, i,j,m
      double precision cset(3,*)
      logical unique
c
      m = 0
      do i=1, n
         unique = .true.
         do j=1, m
            if ( nset(i) .eq. nset(j) ) unique = .false.
         enddo
         if (unique) then
            m = m + 1
            nset(m)   = nset(i)
            cset(1,m) = cset(1,i)
            cset(2,m) = cset(2,i)
            cset(3,m) = cset(3,i)
         endif
      enddo
      n = m
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine write_out_set(uni,set,n)
      implicit none
      integer  uni,set(*),n, nv(12),i,j
c
      j = 0
      do i=1, n
         j = j + 1
         nv(j) = set(i)
         if (j.eq.10) then
            call write_in_nod_or_el_set(uni,nv,j)
            j = 0
         endif
      enddo
      if (j.gt.0) then
         call write_in_nod_or_el_set(uni,nv,j)
      endif
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine write_in_nod_or_el_set(uni,nv,l)
c
c Write out the nodes on unit=unv
c
      implicit none
      integer uni,l,nv(*),i
      double precision r
      character form*7
      if (l.eq.1) then
         r = nv(1)
         r = log10(r)
         i = 1 + int(r)
         write(form,'(t1,a,i1,a)') '(t1,i',i,')'
         write(uni,form) nv(1)
c         write(uni,'(t1,i6)') nv(1)
      else
         write(uni,'(t1,i6,11(a,i6))') nv(1),(',',nv(i), i=2,l)
      endif
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine elemset_(uni,el,iem,text)
c
c write out a set of elemets on unit=uni
c
      implicit none
      integer uni,i,l,iem,el(iem),ev(12)
      character text*20
c
      if (text(1:2) .ne. '  ') write(uni,'(t1,a)') text
      l=0
      do i=1, iem
         if (el(i).gt.0) then
            if (l.lt.12) then
                 l=l+1
               ev(l)=i
            else
               call write_in_nod_or_el_set(uni,ev,12)
               ev(1)=i
               l=1
            endif
         endif
      enddo
      call write_in_nod_or_el_set(uni,ev,l)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
