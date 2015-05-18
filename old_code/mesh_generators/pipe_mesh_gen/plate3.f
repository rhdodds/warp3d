c----67--1---------2---------3---------4---------5---------6---------712
c
c    File: plate3.f   latest modification July 11, 1997
c    "write(*,...)" changed to "write(iws,...)", 7/29/97
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine node_coordinates(iws,nods,ns1,ns2,ns3,
     1           noda,na1,na2,na3,nodb,nb1,nb2,nb3,
     2           yg,xag,xbg,rfm,tcr, keyhole,errorFlag)
c
c Coordinates are generated in all zones
c
      implicit none
c
      include 'plate_common_nod.f'
c
      double precision   hyperbel_l
      external           hyperbel_l
c
      integer  nt1,nt2, iw, errorFlag
      parameter ( nt1=201, nt2=5, iw=21 )
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      integer  tcr,keyhole
c
      integer  iws,i,j,k,is_mid,is1,is2,isd,ist,it,ia1,ia2,j1,j2,
     2         ja,ja1,ja2, k1,ksr2,ka1,ka2,ilocal,icase,n,nk
c
      double precision xag(0:1000),xbg(0:1000),yg(0:1000),
     1                 rfm(100)
c
      double precision fi(201),rf(201),g_xy(2,201,201),
     1                 phi(500),phi4(500), rxy1(500),rxy2(500),
     2                 rxy3(500),rxy4(500),rho0,rho1,rcur,phicur,
     3                 x2star,yp(200),ypp(200),s3(1000),
     4                 psij(200),psijma(200),
     5                 xi,dx,x,x1,x2,y,y1,y2,z,zm,za,z1,z2,
     6                 psi,eta,chi,eps,rs1,rs2,l1,lambda,sl
c
c  double precision functions !
      double precision xaxel,yaxel,rhoy,calc_lambda,
     &                 calc_amap_length,calc_amap_yp,rhoxy
c
      logical          left,x2_lt_x1
c
      integer      mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      integer         kstart,kstep
      common /nblock/ kstart,kstep
c
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
c   =========================================================
c                         Z O N E    S
c   =========================================================
c
      write(iws,'(t10,a)') 'Zone S'
c
c  Determine the angular function
c
      call angular_distribution(phi,a,c,alfa,jms,tcr, iws,errorFlag)
      if(errorFlag.gt.0)then
        write(iws,*)'Stopped after error in angular_distribution subr.'
        return
      end if
      rho0 = sqrt( (c+a) / (c-a) )
c
c TEST
c
c      write(iws,'(t5,a)') ' PHI:'
c      do i=1, jms
c         write(iws,'(t3,a,i3,g12.4)') ' i=',i,phi(i)
c      enddo
c      write(iws,'(t5,a)') ' RFM:'
c      do i=1, mr
c         write(iws,'(t3,a,i3,g12.4)') ' i=',i,rfm(i)
c      enddo
c
c Define the focused mesh perpendicular to the front.
c Use a local coordinate system; origin: x=y=0.
c Refered to global coord. origin X=a, Y=Z=0
c
c Define the angles in the focused mesh
      fi(1) = 0.0
      fi(ims) = 180.0
      do i=1, ims 
         psi = real(i-1)/real(ims-1)
         fi(i) = (1.-psi)*fi(1) + psi*fi(ims)
      enddo
c Check if mesh coarsening is used in Zone S
      if (sfred_type.eq.1) then
         ksr2 = ksr1
      else
         ksr2 = ksr1 + 2
         do k=mr, sfred, -1
            rfm(k+1) = rfm(k)
         enddo
         rfm(sfred) = 0.5*( rfm(sfred-1) + rfm(sfred+1) )
      endif
      if (sjred_type.gt.1) then
         if (sfred_type.eq.1) then
            rfm(mr+1) = rfm(mr)
            rfm(mr) = 0.5*(rfm(mr-1)+rfm(mr+1))
         else
            rfm(mr+2) = rfm(mr+1)
            rfm(mr+1) = 0.5*(rfm(mr)+rfm(mr+2))
         endif
      endif
c Retrieve the radii in the focused mesh 
c  if keyhole = 0 a sharp crack tip is modeled
c  if keyhole = 1 the crack tip is modeled as a keyhole
c
c Note: in case of a quadratic element midnodes at 1/4-points
c
      if (keyhole.ne.0) then
c         rf(1) = 0.40*rfm(1)
         rf(1) = 0.20*rfm(1)
         rf(2) = 0.25*( rf(1) + rfm(1) )
         write(*,'(t1,a,g12.4/t1,a,f12.8)')
     &        '>> Maximun size of keyhole at cracktip R =',rf(1),
     &        '                                     R/a =',rf(1)/a
      else
         rf(1) = 0.0
         rf(2) = 0.25*rfm(1)
      endif
c
      do k=3, kms, 2
         j = (k-1)/2
         rf(k) = rfm(j)
         rf(k-1) = 0.5*(rf(k)+rf(k-2))
      enddo
c
      if (keyhole.ne.0) then
         k1 = 1
      else
         do i=1, ims
            g_xy(1,i,1) = 0.0
            g_xy(2,i,1) = 0.0
         enddo
         k1 = 2
      endif
      do k=k1, kms
         do i=1, ims
            g_xy(1,i,k) = rf(k)*cosd(fi(i))
            g_xy(2,i,k) = rf(k)*sind(fi(i))
         enddo
      enddo
c
c Remote region of Zone S
c
c  Define the remote boundary
c
      is_mid = (ims+1)/2
      is1 = 2*mv*sfred_type + 1
      is2 = ims - is1 + 1
      isd = sfred_type
      it = 0
      do i=1,is1
         g_xy(1,i,kms) = g_xy(1,1,kms)
      enddo
      do i=1+2*isd, is1, 2*isd
         it = it + 1
         g_xy(2,i,kms) = yg(it)
         j1 = i-2*isd
         j2 = i
         do j = j1+1, j2-1
            psi =	real(j-j1)/real(j2-j1)
            g_xy(2,j,kms)=(1.0-psi)*g_xy(2,j1,kms)+psi*g_xy(2,j2,kms)
         enddo
      enddo
c
      do i=ims, is2, -1
         j = ims+1-i
         g_xy(2,i,kms) = g_xy(2,j,kms)
         g_xy(1,i,kms) = g_xy(1,ims,kms)
      enddo
c
      do i=is1+1, is2-1
         psi = real(i-is1) / real(is2-is1)
         g_xy(1,i,kms) = (1.0-psi)*g_xy(1,is1,kms) + psi*g_xy(1,is2,kms)
         g_xy(2,i,kms) = g_xy(2,is1,kms)
      enddo
c
c Define the coordinates in the remote region
c
      do i=2, ims
         do k=ksr2+1, kms-1
            x1 = g_xy(1,1,ksr2) 
            x2 = g_xy(1,1,kms) 
            y1 = g_xy(2,1,ksr2) 
            y2 = g_xy(2,1,kms)
            psi = ( g_xy(1,1,k) - x1 ) / ( x2 - x1 ) 
            g_xy(1,i,k) = (1.0-psi)*g_xy(1,i,ksr2) + psi*g_xy(1,i,kms)
            g_xy(2,i,k) = (1.0-psi)*g_xy(2,i,ksr2) + psi*g_xy(2,i,kms)
         enddo
      enddo
c
c  Determine the coordinates for 1=<j=<jms
c
      zm = a
ccc      if (keyhole.ne.0) then
ccc         zm = a - rf(1)
ccc      else
ccc         zm = a
ccc      endif
      j = 1
      do k=1, kms, 2
         do i=1, ims, 2
            npos(nods(i,j,k),1) = 0.0
            npos(nods(i,j,k),2) = g_xy(2,i,k)
            npos(nods(i,j,k),3) = zm + g_xy(1,i,k)
         enddo
      enddo
c
c      open( unit = iw, file = 'zones.plt', status = 'unknown' )
c      write(iw,'(t1,a)') 'ZONE'
c      do k=1, kms, 2
c         do i=1, ims, 2
c            write(iw,'(t1,2g15.6)') g_xy(1,i,k),g_xy(2,i,k)
c         enddo
c      enddo
c      close(iw)
c
      z1 = npos(nods(ims,1,kms),3)
      za = a
ccc      if (keyhole.ne.0) then
ccc         za = a - rf(1)
ccc      else
ccc         za = a
ccc      endif
      rs1 = rhoy(z1,phi(1),alfa)
      rs2 = rho0
      eps = 1.0e-3
      do j=3, jms, 2
        do k=1, kms, 2
          do i=1, ims, 2
            z = npos(nods(i,1,k),3)
            psi = (z-z1) / (za-z1)
            sl = z - za
            if     ( abs((z-z1)/z1) .le. eps ) then
              rcur = rhoy(z,phi(1),alfa)
            elseif ( abs((z-za)/za) .le. eps ) then
              rcur = rhoy(z,phi(1),alfa)
            elseif ( z .gt. za ) then
              icase = 1
              call mapping_s(alfa,phi(j),icase,rcur,rs1,rs2,psi,sl,
     &             eps,iws)
            else
              icase = 2
              call mapping_s(alfa,phi(j),icase,rcur,rs1,rs2,psi,sl,
     &             eps,iws)
            endif
            npos(nods(i,j,k),1) = xaxel(rcur,phi(j),alfa)
            npos(nods(i,j,k),2) = g_xy(2,i,k)
            npos(nods(i,j,k),3) = yaxel(rcur,phi(j),alfa)
          enddo
        enddo
      enddo
c
c  Adjust the node coordinates on the crack plane
c
      do j=1, jms
         if (keyhole.ne.0) then
            npos(nods(1,j,1),2)   = 0.0
            npos(nods(ims,j,1),2)   = 0.0
            do k=2, kms
               npos(nods(1,j,k),2)   = 0.0
               npos(nods(ims,j,k),2) = 0.0
            enddo
         else
            do k=1, kms
               npos(nods(1,j,k),2)   = 0.0
               npos(nods(ims,j,k),2) = 0.0
            enddo
         endif
      enddo 
c
c Mid and Surface nodes
c
      do j=1, jms, 2
c
         do k=1, kms, 2
            do i=2, ims, 2
               npos(nods(i,j,k),1) = ( npos(nods(i-1,j,k),1) +
     &                                 npos(nods(i+1,j,k),1) ) / 2.
               npos(nods(i,j,k),2) = ( npos(nods(i-1,j,k),2) +
     &                                 npos(nods(i+1,j,k),2) ) / 2.
               npos(nods(i,j,k),3) = ( npos(nods(i-1,j,k),3) +
     &                                 npos(nods(i+1,j,k),3) ) / 2.
            enddo
         enddo
         do i=1, ims, 2
            do k=2, kms, 2
               npos(nods(i,j,k),1) = ( npos(nods(i,j,k-1),1) +
     &                                 npos(nods(i,j,k+1),1) ) / 2.
               npos(nods(i,j,k),2) = ( npos(nods(i,j,k-1),2) +
     &                                 npos(nods(i,j,k+1),2) ) / 2.
               npos(nods(i,j,k),3) = ( npos(nods(i,j,k-1),3) +
     &                                 npos(nods(i,j,k+1),3) ) / 2.
            enddo
         enddo
c
c Centroid or Mid-surface nodes:
c
         do k=2, kms, 2
            do i=2, ims, 2
               npos(nods(i,j,k),1) =
     &          0.25*( npos(nods(i-1,j,k),1) + npos(nods(i+1,j,k),1) +
     &                 npos(nods(i,j,k-1),1) + npos(nods(i,j,k+1),1) )
               npos(nods(i,j,k),2) =
     &          0.25*( npos(nods(i-1,j,k),2) + npos(nods(i+1,j,k),2) +
     &                 npos(nods(i,j,k-1),2) + npos(nods(i,j,k+1),2) )
               npos(nods(i,j,k),3) =
     &          0.25*( npos(nods(i-1,j,k),3) + npos(nods(i+1,j,k),3) +
     &                 npos(nods(i,j,k-1),3) + npos(nods(i,j,k+1),3) )
            enddo
         enddo
c
      enddo
c
      do j=2, jms, 2
         do k=2, kms
            do i=1, ims
               npos(nods(i,j,k),1) = 0.5*( npos(nods(i,j-1,k),1) +
     &                                     npos(nods(i,j+1,k),1) )
               npos(nods(i,j,k),2) = 0.5*( npos(nods(i,j-1,k),2) +
     &                                     npos(nods(i,j+1,k),2) )
               npos(nods(i,j,k),3) = 0.5*( npos(nods(i,j-1,k),3) +
     &                                     npos(nods(i,j+1,k),3) )
            enddo
         enddo
      enddo
c
c  Modify the coordinates in the mesh coarsening region
c
      if (sfred_type.eq.2) then
         k = ksr1
         do j=1, jms
            do i=3, ims-2, 4
               psi = abs(real(i-is_mid)) / real(ims-1)
               npos(nods(i,j,k),1) = (1.-psi)*npos(nods(i,j,k),1) +
     &                                    psi*npos(nods(i,j,k-2),1)
               npos(nods(i,j,k),2) = (1.-psi)*npos(nods(i,j,k),2) +
     &                                    psi*npos(nods(i,j,k-2),2)
               npos(nods(i,j,k),3) = (1.-psi)*npos(nods(i,j,k),3) +
     &                                    psi*npos(nods(i,j,k-2),3)
            enddo
         enddo
         left=.true.
         ilocal=1
         call reduce_coord_2_to_1(nods,ns1,ns2,ns3,ilocal,3,ims-2,4,
     &                            1,jms,1,ksr1,left)
         ist = 2
      elseif (sfred_type.eq.3) then
         ilocal=1
         call reduce_coord_plane_3_to_1(nods,ns1,ns2,ns3,ilocal,
     &                            4,ims,6,1,jms,1,ksr1)
         ist = 3
      endif
      if ((sfred_type.eq.2).or.(sfred_type.eq.3)) then
         do j=1, jms
c
            do k=ksr2, kms, 2
               do i=1+ist, ims, 2*ist
                  npos(nods(i,j,k),1) = ( npos(nods(i-ist,j,k),1) +
     &                                    npos(nods(i+ist,j,k),1) ) / 2.
                  npos(nods(i,j,k),2) = ( npos(nods(i-ist,j,k),2) +
     &                                    npos(nods(i+ist,j,k),2) ) / 2.
                  npos(nods(i,j,k),3) = ( npos(nods(i-ist,j,k),3) +
     &                                    npos(nods(i+ist,j,k),3) ) / 2.
               enddo
            enddo
            do k=ksr2+1, kms, 2
               do i=1+ist, ims, 2*ist
                  npos(nods(i,j,k),1) =  0.25 *
     &            ( npos(nods(i-ist,j,k),1) + npos(nods(i+ist,j,k),1) +
     &              npos(nods(i,j,k-1),1) + npos(nods(i,j,k+1),1) )
                  npos(nods(i,j,k),2) =  0.25 *
     &            ( npos(nods(i-ist,j,k),2) + npos(nods(i+ist,j,k),2) +
     &              npos(nods(i,j,k-1),2) + npos(nods(i,j,k+1),2) )
                  npos(nods(i,j,k),3) =  0.25 *
     &            ( npos(nods(i-ist,j,k),3) + npos(nods(i+ist,j,k),3) +
     &              npos(nods(i,j,k-1),3) + npos(nods(i,j,k+1),3) )
               enddo
            enddo
c
         enddo
      endif
c
      if (sjred_type.eq.2) then
         left=.true.
         ilocal=2
         call reduce_coord_2_to_1(nods,ns1,ns2,ns3,ilocal,1,ims,1,
     &                            3,jms-2,4,kms-2,left)
      elseif (sjred_type.eq.3) then
         ilocal=2
         call reduce_coord_plane_3_to_1(nods,ns1,ns2,ns3,ilocal,
     &                                  1,ims,1, 4,jms,6,kms-2)
      endif
c
c Make sure that coordinates on planes of symmetry
c  have appropriate values
c   Plane of Symmetry (X=0)
      do k=1, kms
         do i=1, ims
            npos(nods(i,1,k),1)=0.0
         enddo
      enddo
c . Free surface (Z=0)
      do k=1, kms
         do i=1, ims
            npos(nods(i,jms,k),3)=0.0
         enddo
      enddo
c
c      call plot_zone_s(nods,ns1,ns2,ns3) 
c      if (ims.gt.0) then
c          errorFlag = 1
c          return
c       endif
c
c   =========================================================
c                         Z O N E    A
c   =========================================================
c
      write(iws,'(t15,a)') 'Zone A'
c
c Indicies:
c
      ia1 = 2*m1 + 1
      ia2 = 2*(m1+mh+mh) + 1
      ka1 = 2*mv + 1
c
c  Transfer the angular function from Zone S to Zone A
c
      ja = 0
      if (sjred_type.eq.2) then
         do j=1, jms, 2
            ja = ja + 1
            phi(ja) = phi(j)
         enddo
         if (ja.ne.jma) then
           write(iws,'(t1,a,a)') '>> ERROR in transfer of angular fcn.',
     &        ' [ subroutine node_coordinates(...) ]'
            errorFlag = 1
            return
         endif
      elseif (sjred_type.eq.3) then
         do j=1, jms, 3
            ja = ja + 1
            phi(ja) = phi(j)
         enddo
         if (ja.ne.jma) then
           write(iws,'(t1,a,a)') '>> ERROR in transfer of angular fcn.',
     &        ' [ subroutine node_coordinates(...) ]'
            errorFlag = 1
            return
         endif
      endif
c
c  Determine the coordinates at X = Y = 0 ( j=1  &  k=1 )
c
      i = 1
      g_xy(1,i,1) = 0.0
      g_xy(2,i,1) = 0.0
      i = ia1
      g_xy(1,i,1) = 0.0
      g_xy(2,i,1) = npos(noda(ia1,1,1),3)
      do i=ia1, ia2
         g_xy(1,i,1) = 0.0
         g_xy(2,i,1) = npos(noda(i,1,ka1),3)
      enddo
      i = ia2
      g_xy(1,i,1) = 0.0
      g_xy(2,i,1) = npos(noda(ia2,1,1),3)
      i = ima
      g_xy(1,i,1) = 0.0
      g_xy(2,i,1) = t
c
      j = 1
      if (sjred_type.eq.1) then
         dx = npos(nods(ims,j,kms-2),3) - npos(nods(ims,j,kms),3)
         x1 = npos(nods(ims,j,kms-2),3)
         x2 = t - npos(nods(1,j,kms-2),3)
      else
         dx = npos(nods(ims,j,kms-4),3) - npos(nods(ims,j,kms),3)
         x1 = npos(nods(ims,j,kms-4),3)
         x2 = t - npos(nods(1,j,kms-4),3)
      endif
c
c   Mesh coarsening - elements below the crack front
c
      j = 1
      n = m1 + 1
      lambda = exp( calc_lambda(dx,x1,n) )
      if ( (x1/real(n)) .gt. (1.0001*dx) ) then
         it=0
         k = 0.8*n
         if (k.ge.n) k = n - 1          
         call calc_bias_int1(dx,x1,n,k,lambda)
         do i=ia1-2, 3, -2
            it = it + 1
            xi = it
            if ( it .ge. k ) xi = k - 1
            g_xy(1,i,j) = 0.0
            g_xy(2,i,j) = g_xy(2,i+2,j) - dx*lambda**xi
        enddo
      else
         do i=3, ia1, 2
            psi = real(i-1) / real(ia1-1)
            g_xy(1,i,j) = 0.0
            g_xy(2,i,j) = (1.-psi)*g_xy(2,1,j) + psi*g_xy(2,ia1,j)
         enddo
      endif
c
c   Mesh coarsening -  elements above the crack front
c
      j = 1
      n = m2 + 1
      nk = 0.8*n
      if (nk.ge.n) nk = n - 1          
      if ( (x2/real(n)) .gt. (1.0001*dx) ) then
         call calc_bias_int1(dx,x2,n,nk,lambda)
c         write(*,'(t2,a,2i3,a,g11.4)') 'n,nk:',n,nk,' x2=',x2
c         write(*,'(t4,2(a,g11.4))') ' lambda=',lambda
         it=0
         do i=ia2+2, ima-2, 2
            it = it + 1
            xi = it
            if ( it .ge. nk ) xi = nk - 1
            g_xy(1,i,j) = 0.0
            g_xy(2,i,j) = g_xy(2,i-2,j) + dx*lambda**xi
         enddo
      else
         do i=ia2+2, ima-2, 2
            psi = real(i-ia2) / real(ima-ia2)
            g_xy(1,i,j) = 0.0
            g_xy(2,i,j) = (1.-psi)*g_xy(2,ia2,j) + psi*g_xy(2,ima,j)
         enddo
      endif
c
c  Determine the boundaries of Zone A
c    => express in terms of radius rho and angle phi
c
      rxy1(1) = rhoy( g_xy(2,5,1) , phi(1), alfa)
      do j=2, jma
         rxy1(j) = rxy1(1)
      enddo
      rxy2(1) = rhoy( npos(noda(ia1,1,1),3) , phi(1), alfa )
      do j=1, jma,2
         rxy2(j) = rxy2(1)
         x1 = npos(noda(ia2,j,1),1)
         y1 = npos(noda(ia2,j,1),3)
         rxy3(j) = rhoxy( x1, y1, phi(j), alfa )
      enddo
c
c  Region above the front
c
      ja1 = 2*(na-nb)+1
      phi4(1) = phi(1)
      rxy4(1) = rhoy( t, phi4(1), alfa ) 
c    initial guess for rxy4(ja1) and phi4(ja1)
      phi4(ja1) = phi(ja1)
      rxy4(ja1) = rhoy( t, phi4(ja1), alfa ) 
      x2star = xaxel( rxy4(ja1), phi(ja1), alfa )
      x1 = xag(m2)
      if (x2star .lt. (x1-0.02*t) ) then
         x2 = x2star
         y1 = t / sqrt(1.0 - (x2/x1)**2.0 )
         x2_lt_x1 = .true.
      else
         x2 = x1
         x2_lt_x1 = .false.
      endif
c
c   Find rho(j) and phi(j) on  0 < X < X2, Y=0, Z=T
c
      if ( x2_lt_x1 ) then
         do j=1, ja1, 2
            phi4(j) = phi(j)
            rxy4(j) = rhoy( t, phi4(j), alfa )
         enddo
      else
         call calc_rf( x2, t, alfa, rxy4(ja1), phi4(ja1) )
         do j=3, ja1-2, 2
            psi = ( phi(1) - phi(j)) / ( phi(1) - phi(ja1) )
            phi4(j) = (1.0-psi)*phi4(1) + psi*phi4(ja1)
            rxy4(j) = rhoy(t,phi4(j),alfa)
         enddo
      endif
c
c   Find rho(j) and phi(j) on (x/x1)^2 + (z/z1)^2 = 1, Y=0
c
      dx = xag(m2) - xag(0)
      eta = dx / t
      chi = 1.25
      if (eta .gt. 1) eta = 1.0         
      do j=ja1, jma
         ypp(j) = t * real(jma-j)/real(jma-ja1)
      enddo
      s3(jma) = 0.0
      do j=jma-2, ja1, - 2
         x = npos(noda(ia2,j,1),1) - npos(noda(ia2,j+2,1),1)
         y = npos(noda(ia2,j,1),3) - npos(noda(ia2,j+2,1),3)
         s3(j) = s3(j+2) + sqrt( x*x + y*y )
      enddo
      if ( x2_lt_x1 ) then
         do j=ja1, jma-2, 2
            psi = s3(j) / s3(ja1)
            yp(j) =  calc_amap_yp(x1,y1,t,psi,eps,iws)
         enddo
         do j=ja1, jma-2, 2
            y = yp(j)*(1.0-eta**chi) + ypp(j)*eta**chi
            x = x1 * sqrt( 1.0 - (y/y1)**2.0 )
            call calc_rf( x, y, alfa, rxy4(j), phi4(j) )
         enddo
      else
         do j=ja1, jma, 2
            psi = s3(j) / s3(ja1)
            yp(j) = psi*t
         enddo
         do j=ja1, jma-2, 2
            y = yp(j)*(1.0-eta**chi) + ypp(j)*eta**chi
            x = x1
            call calc_rf( x, y, alfa, rxy4(j), phi4(j) )
         enddo
      endif
c
c  Determine the node coordinates in the front layer, Y=0
c
      do i=5, ia1-2, 2
         rho1 = rhoy( g_xy(2,i,1), phi(1), alfa )
         npos(noda(i,1,1),1) = 0.0
         npos(noda(i,1,1),2) = 0.0
         npos(noda(i,1,1),3) = g_xy(2,i,1)
         do j=3, jma, 2
            npos(noda(i,j,1),1) = xaxel( rho1, phi(j), alfa)
            npos(noda(i,j,1),2) = 0.0
            npos(noda(i,j,1),3) = yaxel( rho1, phi(j), alfa)
         enddo
      enddo
c
      ja2 = jma-4
      psi =  g_xy(2,3,1) / g_xy(2,5,1)
      do j=1, jma-6, 2
         npos(noda(1,j,1),1) = npos(noda(5,j,1),1)
         npos(noda(1,j,1),2) = 0.0
         npos(noda(1,j,1),3) = 0.0
c
         npos(noda(3,j,1),1) = npos(noda(5,j,1),1)
         npos(noda(3,j,1),2) = 0.0
         npos(noda(3,j,1),3) = npos(noda(5,j,1),3)*
     &          ( psi*real(ja2-j) + 0.50*real(j) ) / real(ja2)
      enddo
      if (noda(5,jma-4,k).gt.kstart) then
         npos(noda(1,ja2,1),1) = npos(noda(5,ja2,1),1)
         npos(noda(1,ja2,1),2) = 0.0
         npos(noda(1,ja2,1),3) = 0.0
c
         npos(noda(3,ja2,1),1) = npos(noda(5,ja2,1),1)
         npos(noda(3,ja2,1),2) = 0.0
         npos(noda(3,ja2,1),3) = npos(noda(5,j,1),3)*0.50
      endif
c
c Modify the geometry for the element i:3,5; j:jma-6,jma-4; k:1,3
c
      npos(noda(3,jma-6,1),1) = 0.5*( npos(noda(1,jma-8,1),1) +
     &                                npos(noda(1,jma-4,1),1) )
      npos(noda(3,jma-6,1),3) =  npos(noda(3,jma-4,1),3)
c
c Calculate the coordinates for Y=0 above the crack front.
c
      do i=ia2+2, ima, 2
         npos(noda(i,1,1),1) = 0.0
         npos(noda(i,1,1),2) = 0.0
         npos(noda(i,1,1),3) = g_xy(2,i,1)
         it = (i-ia2)/2
         npos(noda(i,jma,1),1) = xag(it)      
         npos(noda(i,jma,1),2) = 0.0
         npos(noda(i,jma,1),3) = 0.0
      enddo
      eps = 1.0e-04
      do j=3, ja1, 2
         psij(ia2) = 0.0
         dx = rfm(mr) - rfm(mr-1)
         x = calc_amap_length(alfa,phi(j),phi4(j),rxy3(j),rxy4(j),
     &       eps,iws)
         n = m2 + 1
         if ((x/real(n)) .gt. (1.0001*dx) ) then
            call calc_bias_int1(dx,x,n,nk,lambda)
            x = x - dx
            it = 0
            do i=ia2+2, ima-2, 2
               it = it + 1
               xi = it
               if (it .ge. nk) xi = nk - 1
               psij(i) = psij(i-2) + (dx/x)*lambda**xi
            enddo
         else
            do i=ia2+2, ima-2, 2
               psij(i) = real(i-ia2)/real(ima-ia2)
            enddo
         endif
         do i=ia2+2, ima-2, 2
            call mapping_a2(alfa,phi(j),phi4(j),rxy3(j),rxy4(j),
     &                      rcur,phicur,psij(i),eps,iws)
            npos(noda(i,j,1),1) = xaxel( rcur, phicur, alfa)
            npos(noda(i,j,1),2) = 0.0
            npos(noda(i,j,1),3) = yaxel( rcur, phicur, alfa)
         enddo
      enddo
      s3(1) = 0.0
      do j=3, jma, 2
         x = npos(noda(ia2,j,1),1) - npos(noda(ia2,j-2,1),1)
         y = npos(noda(ia2,j,1),3) - npos(noda(ia2,j-2,1),3)
         s3(j) = s3(j-2) + sqrt( x*x + y*y )
      enddo
      do i=ia2+2, ima,2
         it = (i-ia2)/2
         psijma(i) = (xag(it)-xag(0))/(xag(m2)-xag(0))
      enddo
      do j=ja1+2, jma-2, 2
         eta = ( s3(j) - s3(ja1) ) / ( s3(jma) - s3(ja1) )
         do i=ia2+2, ima-2, 2
            psi = (1.0-eta)*psij(i) + eta*psijma(i)
            call mapping_a2(alfa,phi(j),phi4(j),rxy3(j),rxy4(j),
     &                      rcur,phicur,psi,eps,iws)
            npos(noda(i,j,1),1) = xaxel( rcur, phicur, alfa)
            npos(noda(i,j,1),2) = 0.0
            npos(noda(i,j,1),3) = yaxel( rcur, phicur, alfa)
         enddo
      enddo
      do j=1, ja1, 2
         npos(noda(ima,j,1),1) = xaxel( rxy4(j), phi4(j), alfa)
         npos(noda(ima,j,1),2) = 0.0
         npos(noda(ima,j,1),3) = t
      enddo
      do j=ja1+2, jma-2, 2
         npos(noda(ima,j,1),1) = xaxel( rxy4(j), phi4(j), alfa)
         npos(noda(ima,j,1),2) = 0.0
         npos(noda(ima,j,1),3) = yaxel( rxy4(j), phi4(j), alfa)
      enddo
c
c  Node coordinates in the region embrazing Zone S (the crack zone)
c
      do k=3, ka1, 2
         do j=1, jma-4, 2
            do i=1, 3, 2
               npos(noda(i,j,k),1) = npos(noda(i,j,1),1)
               npos(noda(i,j,k),2) = yg((k-1)/2)
               npos(noda(i,j,k),3) = npos(noda(i,j,1),3)
            enddo
         enddo
c      ( Modify the element i:3,5; j:jma-6,jma-4; k:1,3)
         npos(noda(3,jma-6,k),1) = 0.5*( npos(noda(1,jma-8,1),1) +
     &                                   npos(noda(1,jma-4,1),1) )
         npos(noda(3,jma-6,k),3) =  npos(noda(3,jma-4,1),3)
         do j=1, jma, 2
            do i=5, ia1-2, 2
               npos(noda(i,j,k),1) = npos(noda(i,j,1),1)
               npos(noda(i,j,k),2) = yg((k-1)/2)
               npos(noda(i,j,k),3) = npos(noda(i,j,1),3)
            enddo
         enddo
         do j=1, jma, 2
            do i=ia2+2, ima, 2
               npos(noda(i,j,k),1) = npos(noda(i,j,1),1)
               npos(noda(i,j,k),2) = yg((k-1)/2)
               npos(noda(i,j,k),3) = npos(noda(i,j,1),3)
            enddo
         enddo
      enddo
c
c  Equal element size at k = kar1-2 or at k = kma
c
      if (rtype.eq.0) then
         ka2 = kma
      else
         ka2 = kar1 - 2
      endif
      x1 = t * real(5-1)/real(ima-1)
      y = yg((ka2-1)/2)
      rho1 = rhoy( x1, phi(1), alfa)
c
c  Region below the crack front
c
      do j=1, jma, 2
         npos(noda(5,j,ka2),1) = xaxel( rho1, phi(j), alfa)
         npos(noda(5,j,ka2),2) = y
         npos(noda(5,j,ka2),3) = yaxel( rho1, phi(j), alfa)
      enddo
c
      do j=1, jma-6, 2
         npos(noda(1,j,ka2),1) = npos(noda(5,j,ka2),1)
         npos(noda(1,j,ka2),2) = y
         npos(noda(1,j,ka2),3) = 0.0
c
         npos(noda(3,j,ka2),1) = npos(noda(5,j,ka2),1)
         npos(noda(3,j,ka2),2) = y
         npos(noda(3,j,ka2),3) = 0.50 * npos(noda(5,j,ka2),3)
      enddo
      if (noda(5,jma-4,k).gt.kstart) then
         npos(noda(1,ja2,ka2),1) = npos(noda(5,ja2,ka2),1)
         npos(noda(1,ja2,ka2),2) = y
         npos(noda(1,ja2,ka2),3) = 0.0
          npos(noda(3,ja2,ka2),1) = npos(noda(5,ja2,ka2),1)
         npos(noda(3,ja2,ka2),2) = y
         npos(noda(3,ja2,ka2),3) = 0.50 * npos(noda(5,j,ka2),3)
      endif
c
c Modify the geometry for the element i:3,5; j:jma-6,jma-4; k:1,3
c
      npos(noda(3,jma-6,ka2),1) = 0.5*( npos(noda(1,jma-8,ka2),1) +
     &                                  npos(noda(1,jma-4,ka2),1) )
      npos(noda(3,jma-6,ka2),3) =  npos(noda(3,jma-4,ka2),3)
c
c  Outer bounderies and interface to Zone B:
c
      do j=1, ja1, 2
         npos(noda(ima,j,ka2),1) = npos(noda(ima,j,1),1)
         npos(noda(ima,j,ka2),2) = y
         npos(noda(ima,j,ka2),3) = t
      enddo
      do j=ja1+2, jma, 2
         npos(noda(ima,j,ka2),1) = npos(noda(ima,j,1),1)
         npos(noda(ima,j,ka2),2) = y
         npos(noda(ima,j,ka2),3) = npos(noda(ima,j,1),3)
      enddo
c
c Region for i > 5; here mapping and linear scaling is used.
c
      eps = 1.0e-04
      do i=7, ima-2, 2
         psi = real(i-1)/real(ima-1)
         npos(noda(i,1,ka2),1) = 0.0
         npos(noda(i,1,ka2),2) = y
         npos(noda(i,1,ka2),3) = psi * t
         psi = real(i-5)/real(ima-5)
         do j=3, jma-4, 2
            call mapping_a2(alfa,phi(j),phi4(j),rho1,rxy4(j),
     &                      rcur,phicur,psi,eps,iws)
            npos(noda(i,j,ka2),1) = xaxel( rcur, phicur, alfa)
            npos(noda(i,j,ka2),2) = y
            npos(noda(i,j,ka2),3) = yaxel( rcur, phicur, alfa)
         enddo
c         do j=3, jma-4, 2
c            call mapping_a1(alfa,phi(j),rcur,rho1,rxy4(j),psi,eps,iws)
c            npos(noda(i,j,ka2),1) = xaxel( rcur, phi(j), alfa)
c            npos(noda(i,j,ka2),2) = y
c            npos(noda(i,j,ka2),3) = yaxel( rcur, phi(j), alfa)
c         enddo
c
         npos(noda(i,jma,ka2),1) = (1.-psi)*npos(noda(1,jma-4,ka2),1) +
     &                               psi *npos(noda(ima,jma,ka2),1)
         npos(noda(i,jma,ka2),2) = y
         npos(noda(i,jma,ka2),3) = 0.0
c
         npos(noda(i,jma-2,ka2),1) = 0.5*( npos(noda(i,jma-4,ka2),1) +
     &                                     npos(noda(i,jma,ka2),1) )
         npos(noda(i,jma-2,ka2),2) = y
         npos(noda(i,jma-2,ka2),3) = 0.5*( npos(noda(i,jma-4,ka2),3) +
     &                                     npos(noda(i,jma,ka2),3) )
c
      enddo                                
c
c   Node coordinates: ka1  <  k  <  ka2; linear interpolation
c
      y1 = yg((ka1-1)/2)
      y2 = yg((ka2-1)/2)
      do k=ka1+2, ka2-2, 2
         y = yg((k-1)/2)
         psi = (y-y1) / (y2-y1)
         do j=1, jma, 2
            do i=1, ima, 2
               npos(noda(i,j,k),1) = (1.-psi)*npos(noda(i,j,ka1),1) +
     &                                   psi *npos(noda(i,j,ka2),1)
               npos(noda(i,j,k),2) = y
               npos(noda(i,j,k),3) = (1.-psi)*npos(noda(i,j,ka1),3) +
     &                                   psi *npos(noda(i,j,ka2),3)
            enddo
         enddo
      enddo
c
c  Adjust the node coordinates on the crack plane
c
      do j=1, jma
         do i=1, ia1
            npos(noda(i,j,1),2) = npos(nods(ims,1,1),2)
         enddo
      enddo
c
c  Adjust coordinates if mesh coarsening, k = 2*(lred-1)+1
c
      if (rtype.eq.1) then
         do i=5, ima, 4
            npos(noda(i,jma-2,ka2),1) = 0.5*( npos(noda(i,jma-4,ka2),1)
     &                                      + npos(noda(i,jma,ka2),1) )
            npos(noda(i,jma-2,ka2),3) = 0.5*( npos(noda(i,jma-4,ka2),3)
     &                                      + npos(noda(i,jma,ka2),3) )
         enddo
         do j=1, jma, 2
            do i=7, ima, 4
               npos(noda(i,j,ka2),1) =  0.5*( npos(noda(i-2,j,ka2),1)
     &                                      + npos(noda(i+2,j,ka2),1) )
	         npos(noda(i,j,ka2),3) = 0.5*( npos(noda(i-2,j,ka2),3)
     &                                       + npos(noda(i+2,j,ka2),3) )
            enddo
         enddo
      endif
      if (rtype.eq.2) then
         do j=3, jma, 4
            do i=1, ima ,2
               npos(noda(i,j,ka2),1) = 0.5*( npos(noda(i,j-2,ka2),1)
     &                                     + npos(noda(i,j+2,ka2),1) )
               npos(noda(i,j,ka2),3) = 0.5*( npos(noda(i,j-2,ka2),3)
     &                                     + npos(noda(i,j+2,ka2),3) )
            enddo
         enddo
         do j=1, jma, 2
            do i=3, ima, 4
               npos(noda(i,j,ka2),1) = 0.5*( npos(noda(i-2,j,ka2),1)
     &                                     + npos(noda(i+2,j,ka2),1) )
               npos(noda(i,j,ka2),3) = 0.5*( npos(noda(i-2,j,ka2),3)
     &                                     + npos(noda(i+2,j,ka2),3) )
            enddo
         enddo
      endif
C..MID-, SURFACE- and CETROID-NODES !!!
C . Midnodes for even k:s
      do k=2, ka2, 2
         do i=1, ima, 2
            do j=1, jma, 2
               npos(noda(i,j,k),1) = 0.5*( npos(noda(i,j,k-1),1)
     &                                   + npos(noda(i,j,k+1),1) )
               npos(noda(i,j,k),2) = 0.5*( npos(noda(i,j,k-1),2)
     &                                   + npos(noda(i,j,k+1),2) )
               npos(noda(i,j,k),3) = 0.5*( npos(noda(i,j,k-1),3)
     &                                   + npos(noda(i,j,k+1),3) )
            enddo
         enddo
      enddo
C . Midnodes for even i:s
      do k=1, ka2
         do i=2, ima, 2
            do j=1, jma, 2
               npos(noda(i,j,k),1) = 0.5*( npos(noda(i-1,j,k),1)
     &                                   + npos(noda(i+1,j,k),1) )
               npos(noda(i,j,k),2) = 0.5*( npos(noda(i-1,j,k),2)
     &                                   + npos(noda(i+1,j,k),2) )
               npos(noda(i,j,k),3) = 0.5*( npos(noda(i-1,j,k),3)
     &                                   + npos(noda(i+1,j,k),3) )
            enddo
         enddo
      enddo
C . Midnodes for even j:s
      do k=1, ka2
         do i=1, ima, 2
            do j=2, jma,2
               npos(noda(i,j,k),1) = 0.5*( npos(noda(i,j-1,k),1)
     &                                   + npos(noda(i,j+1,k),1) )
               npos(noda(i,j,k),2) = 0.5*( npos(noda(i,j-1,k),2)
     &                                   + npos(noda(i,j+1,k),2) )
               npos(noda(i,j,k),3) = 0.5*( npos(noda(i,j-1,k),3)
     &                                   + npos(noda(i,j+1,k),3) )
            enddo
         enddo
      enddo
C . Surface and centroid nodes for even i:s and j:s
      do k=1, ka2
         do i=2, ima, 2
            do j=2, jma, 2
               npos(noda(i,j,k),1) = 0.25 *
     &               ( npos(noda(i-1,j,k),1) + npos(noda(i+1,j,k),1) +
     &                 npos(noda(i,j-1,k),1) + npos(noda(i,j+1,k),1) )
               npos(noda(i,j,k),2) = 0.25 *
     &               ( npos(noda(i-1,j,k),2) + npos(noda(i+1,j,k),2) +
     &                 npos(noda(i,j-1,k),2) + npos(noda(i,j+1,k),2) )
               npos(noda(i,j,k),3) = 0.25 *
     &               ( npos(noda(i-1,j,k),3) + npos(noda(i+1,j,k),3) +
     &                 npos(noda(i,j-1,k),3) + npos(noda(i,j+1,k),3) )
            enddo
         enddo
      enddo
C . The domain under the crack front: i =< 5 and 1 =< j =< jma-4
      if (mod(na,4).ne.0) then
         left=.true.
      else
         left=.false.
      endif
      call reduce_coord_under(noda,na1,na2,na3,3,jma-4,4,1,ka2,1,left)
 
C . The domain in the possibly reduced region for k > ka2
      do k=ka2+2, kma, 2
         do j=1, jma
            do i=1, ima
               npos(noda(i,j,k),1)=npos(noda(i,j,ka2),1)
               npos(noda(i,j,k),2)=yg((k-1)/2)
               npos(noda(i,j,k),3)=npos(noda(i,j,ka2),3)
c 
               npos(noda(i,j,k-1),1)=npos(noda(i,j,ka2),1)
               npos(noda(i,j,k-1),2)=(yg((k-1)/2)+yg((k-3)/2))/2.0
               npos(noda(i,j,k-1),3)=npos(noda(i,j,ka2),3)
            enddo
         enddo
      enddo
c
      if (rtype.ge.1) then
         if (mod(ma,4).eq.0) then
            left=.true.
         else
            left=.false.
         endif
         ilocal=1
         call reduce_coord_2_to_1(noda,na1,na2,na3,ilocal,7,ima,4,
     &                           1,jma,1,kar1,left)
      endif
      if (rtype.eq.2) then
         if (mod(na,4).ne.0) then
            left=.true.
         else
            left=.false.
         endif
         ilocal=2
         call reduce_coord_2_to_1(noda,na1,na2,na3,ilocal,5,ima,1,
     &                        3,jma,4,kar2,left)
C . . . The small limited area under the crack front
         if (mod(na,4).ne.0) then
            left=.true.
         else
            left=.false.
         endif
         call reduce_coord_3_to_1(noda,na1,na2,na3,3,jma,4,kar2,left)
      endif
c
c  Make sure that the outer boundaries to zone A has the correct coord.
c  X = 0:
      do k=1, kma
         do i=1, ima
            npos(noda(i,1,k),1)=0.0
         enddo
      enddo
c
      do k=1, kma
c  Z = t:
         do j=1, 2*(na-nb)+1
            npos(noda(ima,j,k),3)=t
         enddo
c  Z = 0:
         do j=1, jma-4
            npos(noda(1,j,k),3)=0.0
         enddo
         do i=5, ima
            npos(noda(i,jma,k),3)=0.0
         enddo
      enddo
c
c      call plot_zone_a(noda,na1,na2,na3)
c      if (ima.gt.0) then
c          errorFlag = 1
c          return
c      endif
c
c   =========================================================
c                         Z O N E    B
c   =========================================================
c
      if (imb .le. 1) goto 50
c
      write(iws,'(t20,a)') 'Zone B'
c
c Mesh coarsening - element sizes in the global X-direction
c
      it = 0
      do i=3, imb, 2
         it=it+1
         xi = it
         npos(nodb(i,jmb,1),1) = xbg(it)
         npos(nodb(i,jmb,1),2)=0.0
         npos(nodb(i,jmb,1),3)=0.0
      enddo
      npos(nodb(imb,jmb,1),1) = w
c
c Generate node coordinates in region  1 =< k =< ka2 :
c
      chi = (w - npos(nodb(1,jmb,1),1)) / t
      if (chi.gt.1) chi = 1.0
      do k=1, ka2, 2
         i = imb
         do j=1, jmb
            psi = real(j-1)/real(jmb-1)
            npos(nodb(i,j,k),1) = w
            npos(nodb(i,j,k),2) = yg((k-1)/2)
            npos(nodb(i,j,k),3) = chi*(1.-psi)*t
     &                          + (1.0-chi)*npos(nodb(1,j,k),3)
         enddo
         npos(nodb(i,jmb,k),3) = 0.0
      enddo
c
      chi = 1.5
      do k=1, ka2, 2
         l1 = npos(nodb(1,jmb,k),1)
         y = yg( (k-1)/2 )
         do j=1, jmb, 2
            z1 = npos(nodb(  1,j,k),3)
            z2 = npos(nodb(imb,j,k),3)
            x1 = npos(nodb(  1,j,k),1)
            x2 = npos(nodb(imb,j,k),1)
            do i=3, imb-2, 2
               x = npos(nodb(i,jmb,1),1)
               psi = (x-l1) / (w-l1)
               npos(nodb(i,j,k),1) = (1.-psi)*x1 + psi*x2
c               npos(nodb(i,j,k),1) = npos(nodb(i,j,k),1)*(1.-psi)**chi
c     &                             + x * (1.-(1.-psi)**chi)
               npos(nodb(i,j,k),2) = y
               npos(nodb(i,j,k),3) = z1*(1.-psi)**chi + 
     &                               z2*(1.-(1.-psi)**chi)
            enddo
         enddo
      enddo
c
c...Adjust coordinates in case of mesh coarsening
c
      if (rtype.eq.2) then
         do i=3, imb, 2
            do j=3, jmb-2, 4
               npos(nodb(i,j,ka2),1) = 0.5*( npos(nodb(i,j-2,ka2),1)
     &                                     + npos(nodb(i,j+2,ka2),1) )
               npos(nodb(i,j,ka2),3) = 0.5*( npos(nodb(i,j-2,ka2),3)
     &                                     + npos(nodb(i,j+2,ka2),3) )
            enddo
         enddo
      endif
c
c...MID-, SURFACE- and CETROID-NODES !!!
c . . Midnodes for even k:s
      do k=2, ka2, 2
         do i=3, imb, 2
            do j=1, jmb, 2
               npos(nodb(i,j,k),1) = 0.5*( npos(nodb(i,j,k-1),1)
     &                                   + npos(nodb(i,j,k+1),1) )
               npos(nodb(i,j,k),2) = 0.5*( npos(nodb(i,j,k-1),2)
     &                                   + npos(nodb(i,j,k+1),2) )
               npos(nodb(i,j,k),3) = 0.5*( npos(nodb(i,j,k-1),3)
     &                                   + npos(nodb(i,j,k+1),3) )
            enddo
         enddo
      enddo
C . Midnodes for even i:s
      do k=1, ka2
         do i=2, imb, 2
            do j=1, jmb, 2
               npos(nodb(i,j,k),1) = 0.5*( npos(nodb(i-1,j,k),1)
     &                                   + npos(nodb(i+1,j,k),1) )
               npos(nodb(i,j,k),2) = 0.5*( npos(nodb(i-1,j,k),2)
     &                                   + npos(nodb(i+1,j,k),2) )
               npos(nodb(i,j,k),3) = 0.5*( npos(nodb(i-1,j,k),3)
     &                                   + npos(nodb(i+1,j,k),3) )
            enddo
         enddo
      enddo
C . Midnodes for even j:s
      do k=1, ka2
         do i=3, imb, 2
            do j=2, jmb,2
               npos(nodb(i,j,k),1) = 0.5*( npos(nodb(i,j-1,k),1)
     &                                   + npos(nodb(i,j+1,k),1) ) 
               npos(nodb(i,j,k),2) = 0.5*( npos(nodb(i,j-1,k),2)
     &                                   + npos(nodb(i,j+1,k),2) )
               npos(nodb(i,j,k),3) = 0.5*( npos(nodb(i,j-1,k),3)
     &                                   + npos(nodb(i,j+1,k),3) )
            enddo
         enddo
      enddo
C . Surface and centroid nodes for even i:s and j:s
      do k=1, ka2
         do i=2, imb, 2
            do j=2, jmb, 2
               npos(nodb(i,j,k),1) = 0.25 *
     &            ( npos(nodb(i-1,j,k),1) + npos(nodb(i+1,j,k),1) +
     &              npos(nodb(i,j-1,k),1) + npos(nodb(i,j+1,k),1) )
               npos(nodb(i,j,k),2) = 0.25 *
     &            ( npos(nodb(i-1,j,k),2) + npos(nodb(i+1,j,k),2) +
     &              npos(nodb(i,j-1,k),2) + npos(nodb(i,j+1,k),2) )
               npos(nodb(i,j,k),3) = 0.25 *
     &            ( npos(nodb(i-1,j,k),3) + npos(nodb(i+1,j,k),3) +
     &              npos(nodb(i,j-1,k),3) + npos(nodb(i,j+1,k),3) )
            enddo
         enddo
      enddo
C . The domain in the possibly reduced region for k > ka2
      do k=ka2+2, kma, 2
         do j=1, jmb
            do i=2, imb
               npos(nodb(i,j,k),1) = npos(nodb(i,j,ka2),1)
               npos(nodb(i,j,k),2) = yg((k-1)/2)
               npos(nodb(i,j,k),3) = npos(nodb(i,j,ka2),3)
c
               npos(nodb(i,j,k-1),1) = npos(nodb(i,j,ka2),1)
               npos(nodb(i,j,k-1),2) = (yg((k-1)/2)+yg((k-3)/2))/2.0
               npos(nodb(i,j,k-1),3) = npos(nodb(i,j,ka2),3)
            enddo
         enddo
      enddo
C . If reduction of elements in Zon B:
      if (rtype.eq.2) then
         if (mod(nb,4).ne.0) then
            left=.true.
         else
            left=.false.
         endif
         ilocal=2
         call reduce_coord_2_to_1(nodb,nb1,nb2,nb3,ilocal,2,imb,1,
     &                           3,jmb,4,kar2,left)
      endif
c
c      call plot_zone_b(nodb,nb1,nb2,nb3)
c      if (imb.gt.0) then
c          errorFlag = 1
c          return
c      endif
c
c      call plot_zone_s(nods,ns1,ns2,ns3)
c      call plot_zone_a(noda,na1,na2,na3)
c      call plot_zone_b(nodb,nb1,nb2,nb3)
c      if (ims.gt.0) then
c          errorFlag = 1
c          return
c      endif
c 
 50   continue
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
c  Subroutines below is only used for testin and debugging the program
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      subroutine plot_zone_s(nods,ns1,ns2,ns3)
c
      implicit none
      include 'plate_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3)
      integer  i,j,k, io,ksr2
      double precision x,y
c
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
      io = 27
      if (sfred_type.eq.1) then
         ksr2 = ksr1
      elseif ((sfred_type.eq.2).or.(sfred_type.eq.3)) then
         ksr2 = ksr1 + 2
      endif
      open(unit=io,file='zoneS1.plt',status='unknown')
c
      if (sfred_type.eq.1) then
c
         j = 1
c
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
c
         j = 1
c
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
c
      endif
c
      close(io)
c
      io = 27
      open(unit=io,file='zoneS2.plt',status='unknown')
c
      if (sfred_type.eq.1) then
c
         j = jms
c
         do i=1, ims, 2
            write(io,'(t1,a)') 'ZONE'
            do k=1, kms, 2
               if (nods(i,j,k).gt.0) then 
                  x = npos(nods(i,j,k),1)
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
                  x = npos(nods(i,j,k),1)
                  y = npos(nods(i,j,k),2)
                  write(io,*) real(x),real(y)
               endif
            enddo
         enddo
c
      elseif (sfred_type.eq.2) then
c
         j = jms
c
         do i=1, ims, 2
            write(io,'(t1,a)') 'ZONE'
            if (mod(i,4).eq.1) then
               do k=1, kms, 2
                  if (nods(i,j,k).gt.0) then 
                     x = npos(nods(i,j,k),1)
                     y = npos(nods(i,j,k),2)
                     write(io,*) real(x),real(y)
                  endif
               enddo
            else
               do k=1, ksr1, 2
                  if (nods(i,j,k).gt.0) then 
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
               if (nods(i,j,k).gt.0) then 
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
      elseif (sfred_type.eq.3) then
c
         j = jms
c
         do i=1, ims, 2
            write(io,'(t1,a)') 'ZONE'
            if (mod(i,6).eq.1) then
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
         x = npos(nods(1,j,ksr2),1)
         y = npos(nods(1,j,ksr2),2)
         write(io,*) real(x),real(y)
         do i=4, ims, 6
            x = npos(nods(i-1,j,ksr1),1)
            y = npos(nods(i-1,j,ksr1),2)
            write(io,*) real(x),real(y)
            x = npos(nods(i+1,j,ksr1),1)
            y = npos(nods(i+1,j,ksr1),2)
            write(io,*) real(x),real(y)
            x = npos(nods(i+3,j,ksr2),1)
            y = npos(nods(i+3,j,ksr2),2)
            write(io,*) real(x),real(y)
         enddo
         do k=ksr2, kms, 2
            write(io,'(t1,a)') 'ZONE'
            do i=1, ims, 3
               x = npos(nods(i,j,k),1)
               y = npos(nods(i,j,k),2)
               write(io,*) real(x),real(y)
            enddo
         enddo
         close(io)
c
      endif
c
      close(io)
c
      io = 27
      open(unit=io,file='zoneS3.plt',status='unknown')
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
      include 'plate_common_nod.f'
c
      integer  na1,na2,na3,noda(na1,na2,na3)
      integer  i,j,k,ia1,ia2,ka1,ka2,io
      double precision x,y
c
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
      integer      mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
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
      open(unit=io,file='zoneA1.plt',status='unknown')
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
      open(unit=io,file='zoneA2.plt',status='unknown')
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
      open(unit=io,file='zoneA3.plt',status='unknown')
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
      open(unit=io,file='zoneA4.plt',status='unknown')
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
      include 'plate_common_nod.f'
c
      integer  nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      integer  i,j,k,ka1,ka2,io
      double precision x,y
c
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
      integer      mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
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
      open(unit=io,file='zoneB1.plt',status='unknown')
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
      open(unit=io,file='zoneB2.plt',status='unknown')
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
      open(unit=io,file='zoneB3.plt',status='unknown')
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
      open(unit=io,file='zoneB4.plt',status='unknown')
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
