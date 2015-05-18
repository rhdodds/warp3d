c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine node_coordinates(nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3,
     1           noda,na1,na2,na3, nodb,nb1,nb2,nb3,
     2           yl,z,job,jobh,pl,etyp,necb,mbtype,mb_bias)
C-----------------------------------------------------------------C
C   Rutinen skapar koordinater f|r noderna i alla zoner .         C
C                                                                 C
C-----------------------------------------------------------------C
      implicit none
c
      include 'common_nod.f'
c
      double precision   hyperbel_l
      external           hyperbel_l
c
      double precision pi
C     parameter ( pi=4.d0*atan(1.) )
      parameter ( pi = 3.141592654)
c
      integer  nt1,nt2
      parameter ( nt1=201, nt2=5 )
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  etyp,jobh,mbtype,pl,nodnum,necb,
     1         i,j,k,is_mid,is1,is2,isd,ist,it,ia1,ia2,
     2         ja,ja1,ja2, kend,ksr2,ka1,ka2,ilocal,icase,n, i2
c
      double precision yl(0:200),z(0:200,2),mb_bias,trc(50,5)
      double precision phi(500),rho0,rho1,rho2,rhotmp(25)
      double precision alfa_y,fi(200,50),g_xy(2,201,61),xy(2,100),
     1                 rxy1(200),rxy2(200),rxy3(200),rxy4(200),
     2                 ds,dx,dxi,x,x1,x2,xm,y,y1,y2,z1,z2,xdummy,
     3                 arg,rf(2,201),lsg,l1,ln,lambda, xcr,sl,
     4                 rcr,rs1,rs2,rcur,psi,eps,chi
c
c  double precision functions !
      double precision xaxel,yaxel,rhox,rhoy,calc_lambda,rho_inc,rhoxy
c
      character  job*40
      logical    left
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      integer         kstart,kstep
      common /nblock/ kstart,kstep
c
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
      double precision      dcell,eta_n,eta_t1,eta_t2
      common /geom_cell/    dcell,eta_n,eta_t1,eta_t2
      double precision      p_lcx1,p_lcx2,p_lcy2,p_alfac
      common /geom_zone_s/  p_lcx1,p_lcx2,p_lcy2,p_alfac
c
	write(*,'(t5,a/)') '* generating the coordinates . . . .'
c
c                   ==========================
c                         Z O N E    C
c                   ==========================
c
      write(*,'(t5,a)') 'Zone C'
c
c  Determine the angular function
c
      call angular_distribution(dcell,eta_t1,eta_t2,phi)
c
c  Determine the coordinates in the comp. cell part. This mesh is
c  aligned along segments of hyperbolas.
c
      rho0 = sqrt( (c+a) / (c-a) )
      do j=1, jmc
c
         k=1
c
         npos(nodc(5,j,k),1) = xaxel(rho0,phi(j),alfa)
         npos(nodc(5,j,k),2) = 0.0
         npos(nodc(5,j,k),3) = yaxel(rho0,phi(j),alfa)
c
         ds = 0.5*(eta_n*dcell)
         rho1 = rho0
         rho2 = rho1 + rho_inc(alfa,phi(j),ds,rho1)
         npos(nodc(7,j,k),1) = xaxel(rho2,phi(j),alfa)
         npos(nodc(7,j,k),2) = 0.0
         npos(nodc(7,j,k),3) = yaxel(rho2,phi(j),alfa)
c
         rho1 = rho2
         do i=9, imc-2, 2
            ds = (eta_n*dcell)
            rho2 = rho1 + rho_inc(alfa,phi(j),ds,rho1)
            npos(nodc(i,j,k),1) = xaxel(rho2,phi(j),alfa)
            npos(nodc(i,j,k),2) = 0.0
            npos(nodc(i,j,k),3) = yaxel(rho2,phi(j),alfa)
            rho1 = rho2
         enddo
         ds = 1.35*eta_n*dcell
         rho2 = rho1 + rho_inc(alfa,phi(j),ds,rho1)
         npos(nodc(imc,j,k),1) = xaxel(rho2,phi(j),alfa)
         npos(nodc(imc,j,k),2) = 0.0
         npos(nodc(imc,j,k),3) = yaxel(rho2,phi(j),alfa)
c
         k = 3
         do i=7, imc, 2
            npos(nodc(i,j,3),1) = npos(nodc(i,j,1),1)
            npos(nodc(i,j,3),2) = dcell/2.
            npos(nodc(i,j,3),3) = npos(nodc(i,j,1),3)
         enddo
c
         x = npos(nodc(7,j,1),1)
         y = npos(nodc(7,j,1),3)
         rho1 = rhoxy(x,y,phi(j),alfa)
         k = 3
         do i=5, 1, -2
            ds = -(eta_n*dcell)
            rho2 = rho1 + rho_inc(alfa,phi(j),ds,rho1)
            npos(nodc(i,j,k),1) = xaxel(rho2,phi(j),alfa)
            npos(nodc(i,j,k),2) = 0.25*dcell
            npos(nodc(i,j,k),3) = yaxel(rho2,phi(j),alfa)
            rho1 = rho2
         enddo
c
c Adjust the hight of the first nceb number of elements
c
         i2 = 5+2*necb
         k = 3
         y1 = npos(nodc(5,j,k),2)
         y2 = npos(nodc(i2,j,k),2)
         do i=5, i2, 2
            dxi = sqrt( (npos(nodc(i,j,k),1)-npos(nodc(5,j,k),1))**2.0
     &            + (npos(nodc(i,j,k),3)-npos(nodc(5,j,k),3))**2.0 )
            dx  = sqrt( (npos(nodc(i2,j,k),1)-npos(nodc(5,j,k),1))**2.0
     &            + (npos(nodc(i2,j,k),3)-npos(nodc(5,j,k),3))**2.0 )
            psi = dxi / dx
            npos(nodc(i,j,k),2) = (1.-psi)*y1 +  psi*y2
         enddo
c
         do k=5, kmc, 2
            do i=1, imc, 2
               npos(nodc(i,j,k),1) = npos(nodc(i,j,3),1)
               npos(nodc(i,j,k),2) = 0.5*real(k-1)*(dcell/2.)
               npos(nodc(i,j,k),3) = npos(nodc(i,j,3),3)
            enddo
         enddo
c
c Adjust the the coordinates behind the crack using conformal mapping.
         if (j.eq.1) then
            i = 1
            do k=3, kmc, 2
               y = npos(nodc(i,j,k),3)
               rhotmp(k) = rhoy(y,phi(j),alfa)
            enddo
         else
            do k=3, kmc, 2
               i = 1
               npos(nodc(i,j,k),1) = xaxel(rhotmp(k),phi(j),alfa)
               npos(nodc(i,j,k),3) = yaxel(rhotmp(k),phi(j),alfa)
               i = 3
               npos(nodc(i,j,k),1) = 0.5*( npos(nodc(i-2,j,k),1) +
     &                                     npos(nodc(i+2,j,k),1) )
               npos(nodc(i,j,k),3) = 0.5*( npos(nodc(i-2,j,k),3) +
     &                                     npos(nodc(i+2,j,k),3) )
            enddo
         endif
c
c Adjust the y-coordinates at end zones of the comp. cell el. layer
         psi = 0.85
         npos(nodc(1,j,kmc),2) = psi*npos(nodc(1,j,kmc),2)
         npos(nodc(3,j,kmc),2) = 0.5*(1.+psi)*npos(nodc(3,j,kmc),2)
         npos(nodc(imc,j,kmc),2) = psi*npos(nodc(imc,j,kmc),2)
         do k=5, kmc-2, 2
            psi = real(k-3)/real(kmc-3)
            npos(nodc(1,j,k),2) = (1.-psi)*npos(nodc(1,j,3),2) +
     &                                 psi*npos(nodc(1,j,kmc),2)
            npos(nodc(3,j,k),2) = ( npos(nodc(1,j,k),2) +
     &                              npos(nodc(3,j,k),2) ) / 2.
         enddo
         do k=3, kmc-4, 2
            psi = real(k-1)/real(kmc-3)
            npos(nodc(imc,j,k),2) = (1.-psi)*npos(nodc(imc,j,1),2) +
     &                                 psi*npos(nodc(imc,j,kmc),2)
         enddo
c
      enddo
c
c  test
c      do i=7, imc, 2
c         dx = npos(nodc(i,1,1),1) - npos(nodc(i-2,1,1),1)
c         dy = npos(nodc(i,1,1),3) - npos(nodc(i-2,1,1),3)
c         ds = sqrt(dx*dx+dy*dy)
c         write(28,*) ' i=',i,'  ds/(0.5*dcell) = ',ds/(0.5*dcell)
c      enddo
c
c Mid and Surface nodes
c
      do j=1, jmc
c
         do k=1, kmc, 2
            do i=2, imc, 2
               npos(nodc(i,j,k),1) = ( npos(nodc(i-1,j,k),1) +
     &                                 npos(nodc(i+1,j,k),1) ) / 2.
               npos(nodc(i,j,k),2) = ( npos(nodc(i-1,j,k),2) +
     &                                 npos(nodc(i+1,j,k),2) ) / 2.
               npos(nodc(i,j,k),3) = ( npos(nodc(i-1,j,k),3) +
     &                                 npos(nodc(i+1,j,k),3) ) / 2.
            enddo
         enddo
c
         npos(nodc(imc,j,kmc-2),1) = ( npos(nodc(imc,j,kmc-4),1) +
     &                                 npos(nodc(imc,j,kmc),1) ) / 2.
         npos(nodc(imc,j,kmc-2),2) = ( npos(nodc(imc,j,kmc-4),2) +
     &                                 npos(nodc(imc,j,kmc),2) ) / 2.
         npos(nodc(imc,j,kmc-2),3) = ( npos(nodc(imc,j,kmc-4),3) +
     &                                 npos(nodc(imc,j,kmc),3) ) / 2.
c
         do i=1, imc, 2
            do k=2, kmc, 2
               npos(nodc(i,j,k),1) = ( npos(nodc(i,j,k-1),1) +
     &                                 npos(nodc(i,j,k+1),1) ) / 2.
               npos(nodc(i,j,k),2) = ( npos(nodc(i,j,k-1),2) +
     &                                 npos(nodc(i,j,k+1),2) ) / 2.
               npos(nodc(i,j,k),3) = ( npos(nodc(i,j,k-1),3) +
     &                                 npos(nodc(i,j,k+1),3) ) / 2.
            enddo
         enddo
c Centroid or Mid-surface nodes:
         do k=2, kmc, 2
            do i=2, imc, 2
               npos(nodc(i,j,k),1) =
     &          0.25*( npos(nodc(i-1,j,k),1) + npos(nodc(i+1,j,k),1) +
     &                 npos(nodc(i,j,k-1),1) + npos(nodc(i,j,k+1),1) )
               npos(nodc(i,j,k),2) =
     &          0.25*( npos(nodc(i-1,j,k),2) + npos(nodc(i+1,j,k),2) +
     &                 npos(nodc(i,j,k-1),2) + npos(nodc(i,j,k+1),2) )
               npos(nodc(i,j,k),3) =
     &          0.25*( npos(nodc(i-1,j,k),3) + npos(nodc(i+1,j,k),3) +
     &                 npos(nodc(i,j,k-1),3) + npos(nodc(i,j,k+1),3) )
            enddo
         enddo
         i = imc-1
         k = kmc-2
         npos(nodc(i,j,k),1) =
     &      0.25*( npos(nodc(i-1,j,k),1) + npos(nodc(i+1,j,k),1) +
     &             npos(nodc(i,j,k-2),1) + npos(nodc(i,j,k+2),1) )
         npos(nodc(i,j,k),2) =
     &      0.25*( npos(nodc(i-1,j,k),2) + npos(nodc(i+1,j,k),2) +
     &             npos(nodc(i,j,k-2),2) + npos(nodc(i,j,k+2),2) )
         npos(nodc(i,j,k),3) =
     &      0.25*( npos(nodc(i-1,j,k),3) + npos(nodc(i+1,j,k),3) +
     &             npos(nodc(i,j,k-2),3) + npos(nodc(i,j,k+2),3) )
c
      enddo
c
c  Modify the coordinates in the mesh coarsening region
c
      left=.true.
      ilocal=1
      k = kmc-2
      call reduce_coord_2_to_1(nodc,nc1,nc2,nc3,ilocal,3,imc-4,4,
     &                         1,jmc,1,k,left)
c
c  Test of prog
c
c      call plot_zone_c(nodc,nc1,nc2,nc3)
c
c                   ==========================
c                         Z O N E    S
c                   ==========================
c
      write(*,'(t10,a)') 'Zone S'
c
c  Define the elliptical sections/planes perpendicular to the front.
c
c   Suitable paramater values are:  p_alfac = 0.45
c                                   p_lcx1 = 1.6
c                                   p_lcx2 = 2.8
c                                   p_lcy2 = 2.6
c
      alfa_y = p_alfac*(npos(nods(1,1,1),3)-npos(nods(ims,1,1),3))/2.
      is_mid = (ims+1)/2
      xm = npos(nods( is_mid, 1, 1 ),3)
      lsg = ( npos(nods(1,1,1),3) - npos(nods(ims,1,1),3) ) / 2.
      x1 = p_lcx1 * lsg
      call calc_rf( x1, 0.0001*x1, alfa_y, rxy2(1), xdummy )
      if ((xdummy-90.).gt.0.001) write(*,*)'>> ERR: Radius in Zone S'
c
      do i=1, ims, 2
         x = npos( nods(i,1,1),3 ) - xm
         y = npos( nods(i,1,1),2 )
         call calc_rf( x, y, alfa_y, rxy1(i), fi(i,1) )
c         write(28,112) ' is=',i,'  x=',x,'  y =',y,
c     &               '  rho=',rxy1(i),'  phi =',fi(i,1)
      enddo
 112  format(t1,a,i3,4(a,g16.8))
c
c  Use the same angular spacing where mesh-coarsening takes place.
c
      fi(1,1) = 0.0
      fi(1,ksr1) = 0.0
      fi(is_mid,ksr1) = 90.
      arg = npos(nods(ims,1,1),2)/( alfa_y*(rxy2(1)+1./rxy2(1)))
      fi(ims,ksr1) = 180. - dasind(arg)
      do i=1, is_mid
         psi = real(i-1) / real(is_mid-1)
         fi(i,ksr1) = (1.-psi)*fi(1,ksr1) + psi*fi(is_mid,ksr1)
      enddo
      do i=is_mid+1, ims-1
         psi = real(i-is_mid) / real(ims-is_mid)
         fi(i,ksr1) = (1.-psi)*fi(is_mid,ksr1) + psi*fi(ims,ksr1)
      enddo
c
      n = sfred
      l1 = p_lcx1 * eta_n*dcell
c
c  On the symmetry plane ahead of the crack
c
      g_xy(1,1,1) = npos(nods(1,1,1),3) - xm
      g_xy(2,1,1) = npos(nods(1,1,1),2)
      g_xy(1,1,ksr1) = xaxel(rxy2(1),fi(1,ksr1),alfa_y)
      g_xy(2,1,ksr1) = yaxel(rxy2(1),fi(1,ksr1),alfa_y)
      ln = g_xy(1,1,ksr1) - g_xy(1,1,1)
      lambda = calc_lambda(l1,ln,n)
      j = 0
      do k=3, ksr1-2, 2
         j = j + 1
         if (exp(lambda).ge.1) then
            g_xy(1,1,k) = g_xy(1,1,k-2) + l1*exp(lambda*real(j-1))
         else
            g_xy(1,1,k) = g_xy(1,1,k-2) + ln*real(j)/real(n)
         endif
         g_xy(2,1,k) = 0.
      enddo
c
c  On the symmetry plane behind the crack
c
      g_xy(1,ims,1) = npos(nods(ims,1,1),3) - xm
      g_xy(2,ims,1) = npos(nods(ims,1,1),2)
      g_xy(1,ims,ksr1) = xaxel(rxy2(1),fi(ims,ksr1),alfa_y)
      g_xy(2,ims,ksr1) = npos(nods(ims,1,1),2)
      ln = abs( g_xy(1,ims,ksr1) - g_xy(1,ims,1) )
      lambda = calc_lambda(l1,ln,n)
      j = 0
      do k=3, ksr1-2, 2
         j = j + 1
         if (exp(lambda).ge.1) then
            g_xy(1,ims,k) = g_xy(1,ims,k-2) - l1*exp(lambda*real(j-1))
         else
            g_xy(1,ims,k) = g_xy(1,ims,k-2) - ln*real(j)/real(n)
         endif
         g_xy(2,ims,k) = g_xy(2,ims,1)
      enddo
c
c  The region above the symmetry plane
c
      do i=3, ims-2, 2
         g_xy(1,i,1) = npos(nods(i,1,1),3) - xm
         g_xy(2,i,1) = npos(nods(i,1,1),2)
         g_xy(1,i,ksr1) = xaxel(rxy2(1),fi(i,ksr1),alfa_y)
         g_xy(2,i,ksr1) = yaxel(rxy2(1),fi(i,ksr1),alfa_y)
         call div_hyperbel(rxy1(i),rxy2(1), fi(i,1),fi(i,ksr1),
     &                     n ,l1,alfa_y,xy,rf)
         j = 0
         do k=3, ksr1-2, 2
            j = j + 1
            g_xy(1,i,k) = xy(1,j)
            g_xy(2,i,k) = xy(2,j)
         enddo
      enddo
c
c      do i=1, ims, 2
c         write(44,'(t1,a,i4)') '  i=',i
c         write(44,114) ' x:',(g_xy(1,i,k), k=1, ksr1, 2)
c         write(44,114) ' y:',(g_xy(2,i,k), k=1, ksr1, 2)
c      enddo
 114  format(t1,a,8g12.5)
c
c  Mesh coarsening
c
      if (sfred_type.eq.1) then
         ksr2 = ksr1
      elseif ((sfred_type.eq.2).or.(sfred_type.eq.3)) then
         ksr2 = ksr1 + 2
         do i=1, ims, 2
            g_xy(1,i,ksr2) = g_xy(1,i,ksr1)
            g_xy(2,i,ksr2) = g_xy(2,i,ksr1)
c
            g_xy(1,i,ksr1) = ( g_xy(1,i,ksr2) + g_xy(1,i,ksr1-2) ) / 2.
            g_xy(2,i,ksr1) = ( g_xy(2,i,ksr2) + g_xy(2,i,ksr1-2) ) / 2.
         enddo
      endif
c
c  Define the outer boundary to Zone S:
c
      is1 = 2*sfred_type*mv + 1
      is2 = ims - 2*sfred_type*mv
      isd = sfred_type
      x2 = p_lcx2 * lsg
c      y2 = p_lcy2 * lsg
c      do i=1, is1
c         psi = real(i-1)/real(is1-1)
c         g_xy(1,i,kms) = x2
c         g_xy(2,i,kms) = psi*y2
c      enddo
      it = -1
      do i=1, is1, 2*isd
         it = it + 1
         g_xy(1,i,kms) = x2
         g_xy(2,i,kms) = yl(it)
      enddo
      do i=1+isd, is1,2*isd
         g_xy(1,i,kms) = x2
         g_xy(2,i,kms) = 0.5*(g_xy(2,i-isd,kms) + g_xy(2,i+isd,kms))
      enddo
c      do i=1, is1, isd
c         write(*,*) ' i=',i,'   y=',real(g_xy(2,i,kms))
c      enddo
c
      do i=is1, is2
         psi = real(i-is1)/real(is2-is1)
         g_xy(1,i,kms) = (1.-psi)*x2 - psi*x2
         g_xy(2,i,kms) = g_xy(2,is1,kms)
c         g_xy(2,i,kms) = y2
      enddo
c
      it = it + 1
      do i=is2, ims, 2*isd
         it = it - 1
         g_xy(1,i,kms) = -x2
         g_xy(2,i,kms) = yl(it)
      enddo
      g_xy(2,ims,kms) = g_xy(2,ims,1)
      do i=is2+isd, ims, 2*isd
         g_xy(1,i,kms) = -x2
         g_xy(2,i,kms) = 0.5*(g_xy(2,i-isd,kms) + g_xy(2,i+isd,kms))
      enddo
c      do i=is2, ims, isd
c         write(*,*) ' i=',i,'   y=',real(g_xy(2,i,kms))
c      enddo
c
c      y1 = g_xy(2,ims,1)
c      do i=is2, ims
c         psi = real(i-is2)/real(ims-is2)
c         g_xy(1,i,kms) = -x2
c         g_xy(2,i,kms) = (1.-psi)*y2 + psi*y1
c      enddo
c
c  Check if y = yl(mv) is large enough!
c
      if (yl(mv).lt.g_xy(2,is_mid,ksr2)) then
          write(*,'(/t5,a,i3/t10,a,g16.8/t10,a,g16.8/)')
     &   'ERROR: increase the Y-coordinate for layer ',mv,
     &   'must be larger than ',g_xy(2,is_mid,ksr2),
     &   'current value = ',yl(mv)
         stop
      endif
c
c  Find the polar coordinates at k=kms
      do i=1, ims
         x = g_xy(1,i,kms)
         y = g_xy(2,i,kms)
         call calc_rf( x, y, alfa_y, rxy3(i), fi(i,kms) )
         fi(i,ksr2) = fi(i,ksr1)
      enddo
c
c  Determine the coordinates in the outer region in Zone S (j=1).
c
      if (sfred_type.eq.1) then
         n = mr - sfred
      elseif ((sfred_type.eq.2).or.(sfred_type.eq.3)) then
         n = mr - sfred - 1
      endif
      if ((sjred_type.eq.2).or.(sjred_type.eq.3)) then
         kend = kms-4
      else
         kend = kms-2
      endif
      do i=1, ims, 2
         l1 = 0.8*sqrt( ( g_xy(1,i,ksr1-2)-g_xy(1,i,ksr2) )**2. +
     &                  ( g_xy(2,i,ksr1-2)-g_xy(2,i,ksr2) )**2. )
c         call div_hyperbel(rxy2(1),rxy3(i), fi(i,ksr2),fi(i,kms),
c     &                        n,l1,alfa_y,xy,rf)
         call div_linear(rxy2(1),rxy3(i),fi(i,ksr2),fi(i,kms),
     &                  n,l1,alfa_y,xy)
         j = 0
         do k=ksr2+2, kend, 2
            j = j + 1
            g_xy(1,i,k) = xy(1,j)
            g_xy(2,i,k) = xy(2,j)
         enddo
      enddo
      if ((sjred_type.eq.2).or.(sjred_type.eq.3)) then
         do i=1, ims, 2
            g_xy(1,i,kms-2) = 0.5*( g_xy(1,i,kms-4) + g_xy(1,i,kms) )
            g_xy(2,i,kms-2) = 0.5*( g_xy(2,i,kms-4) + g_xy(2,i,kms) )
         enddo
      endif
c
c  Determine the coordinates for 1=<j=<jms
c
      j = 1
      do k=3, kms, 2
         do i=1, ims, 2
            npos(nods(i,j,k),1) = 0.0
            npos(nods(i,j,k),2) = g_xy(2,i,k)
            npos(nods(i,j,k),3) = xm + g_xy(1,i,k)
         enddo
      enddo
c
      xcr = abs( npos(nods( is_mid, 1, 1 ),3)-npos(nodc(5,1,1),3) )
      rcr = rho0
      ds = abs( npos(nodc(5,1,1),3) - npos(nodc(5,1,3),3) )
      y1 = xm + g_xy(1,ims,kms)
      y2 = xm - xcr - ds
      rs1 = rhoy(y1,phi(1),alfa)
      eps = 1.0d-4
      do j=3, jms, 2
        rs2 = rcr + rho_inc( alfa, phi(j), -ds, rcr )
        if (rs2.gt.rcr) write(*,*) '>> rs2 > rcr !'
        do k=3, kms, 2
           do i=1, ims, 2
             y = xm + g_xy(1,i,k)
             psi = (y-y1) / (y2-y1)
             sl = g_xy(1,i,k) + xcr + ds
             if (abs((y-y2)/y2).le.eps) then
               rcur = rhoy(y,phi(1),alfa)
             elseif (abs((y-y1)/y1).le.eps) then
               rcur = rhoy(y,phi(1),alfa)
             elseif (y.gt.y2) then
               icase = 1
               call mapping_s(alfa,phi(j),icase,rcur,rs1,rs2,psi,sl,eps)
             else
               icase = 2
               call mapping_s(alfa,phi(j),icase,rcur,rs1,rs2,psi,sl,eps)
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
         do k=1, kms
            npos(nods(ims,j,k),2) = npos(nodc(1,1,3),2)
         enddo
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
c...Make sure that coordinates on planes of symmetry have
c   appropriate values
c . Symmetrisnittet (x=0) !
      do k=1, kmc
         do i=1, imc
            npos(nodc(i,1,k),1)=0.d0
         enddo
      enddo
c . Undersidan (z=0) !
      do k=1, kms
         do i=1, ims
            npos(nods(i,jms,k),3)=0.d0
         enddo
      enddo
c
c      call plot_zone_s(nods,ns1,ns2,ns3,rxy1,rxy2(1),fi,alfa_y)
c      if (ims.gt.0) stop
c
c                   ==========================
c                         Z O N E    A
c                   ==========================
c
      write(*,'(t15,a)') 'Zone A'
c
c Indicies:
c
      ia1 = 2*m1 + 1
      ia2 = 2*(m1+mh+mh) + 1
      ka1 = 2*mv + 1
c
c  Transfer the angular function from Zone S to Zone A
      ja = 0
      if (sjred_type.eq.2) then
         do j=1, jms, 2
            ja = ja + 1
            phi(ja) = phi(j)
         enddo
         if (ja.ne.jma) then
            write(*,'(t1,a,a)') '>> ERROR in transfer of angular fcn.',
     &        ' [ subroutine node_coordinates(...) ]'
            stop
         endif
      elseif (sjred_type.eq.3) then
         do j=1, jms, 3
            ja = ja + 1
            phi(ja) = phi(j)
         enddo
         if (ja.ne.jma) then
            write(*,'(t1,a,a)') '>> ERROR in transfer of angular fcn.',
     &        ' [ subroutine node_coordinates(...) ]'
            stop
         endif
      endif
c
c  Determine the coordinates at Y = 0 ( j=1  &  k=1 )
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
c  Element grading below the crack front
      j = 1
      n = m1 + 1
      lambda = calc_lambda(dx,x1,n)
      if (exp(lambda).gt.(1.0)) then
         it=0
         do i=ia1-2, 3, -2
            it = it + 1
            g_xy(1,i,j) = 0.0
            g_xy(2,i,j) = g_xy(2,i+2,j) - dx*exp(lambda*real(it))
        enddo
      else
         do i=3, ia1, 2
            psi = real(i-1) / real(ia1-1)
            g_xy(1,i,j) = 0.0
            g_xy(2,i,j) = (1.-psi)*g_xy(2,1,j) + psi*g_xy(2,ia1,j)
         enddo
      endif
c
c  Element grading above the crack front
      j = 1
      n = m2 + 1
      lambda = calc_lambda(dx,x2,n)
      if (exp(lambda).gt.(1.0)) then
         it=0
         do i=ia2+2, ima-2, 2
            it = it + 1
            g_xy(1,i,j) = 0.0
            g_xy(2,i,j) = g_xy(2,i-2,j) + dx*exp(lambda*real(it))
         enddo
      else
         do i=ia2+2, ima-2, 2
            psi = real(i-ia2) / real(ima-ia2)
            g_xy(1,i,j) = 0.0
            g_xy(2,i,j) = (1.-psi)*g_xy(2,ia2,j) + psi*g_xy(2,ima,j)
         enddo
      endif
c
c  Determine the boundaries of Zone A (express in terms of radius rho)
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
      ja1 = 2*(na-nb)+1
      do j=1, ja1, 2
         rxy4(j) = rhoy( t, phi(j), alfa )
      enddo
      x1 = xaxel( rxy4(ja1), phi(ja1), alfa )
      dx = t - npos(noda(ia2,1,1),3)
      x2 = dx + npos(noda(ia2,jma,1),1)
      if (x1.lt.x2) then
         rxy4(jma) = rhox( x2, phi(jma), alfa)
         do j=ja1+2, jma-2, 2
            psi = ( phi(ja1)-phi(j) ) / phi(ja1)
            rxy4(j) = (1.0-psi)*rxy4(ja1) + psi*rxy4(jma)
        enddo
      else
         do j=ja1+2, jma, 2
            rxy4(j) = rhox( x2, phi(j), alfa)
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
     &          ( psi*dble(ja2-j) + 0.50*dble(j) ) / dble(ja2)
      enddo
      if (noda(5,jma-4,k).gt.kstart) then
         npos(noda(1,ja2,1),1) = npos(noda(5,ja2,1),1)
         npos(noda(1,ja2,1),2) = 0.0
         npos(noda(1,ja2,1),3) = 0.0
c
         npos(noda(3,ja2,1),1) = npos(noda(5,ja2,1),1)
         npos(noda(3,ja2,1),2) = 0.0
         npos(noda(3,ja2,1),3) = npos(noda(5,j,1),3)*0.5d0
      endif
c
c Modify the geometry for the element i:3,5; j:jma-6,jma-4; k:1,3
c
      npos(noda(3,jma-6,1),1) = 0.5*( npos(noda(1,jma-8,1),1) +
     &                                npos(noda(1,jma-4,1),1) )
      npos(noda(3,jma-6,1),3) =  npos(noda(3,jma-4,1),3)
c
      eps = 1.0d-04
      do i=ia2+2, ima-2, 2
         psi = (g_xy(2,i,1)-g_xy(2,ia2,1))/(g_xy(2,ima,1)-g_xy(2,ia2,1))
         npos(noda(i,1,1),1) = 0.0
         npos(noda(i,1,1),2) = 0.0
         npos(noda(i,1,1),3) = g_xy(2,i,1)
         do j=3, jma, 2
            call mapping_a(alfa,phi(j),rcur,rxy3(j),rxy4(j),psi,eps)
            npos(noda(i,j,1),1) = xaxel( rcur, phi(j), alfa)
            npos(noda(i,j,1),2) = 0.0
            npos(noda(i,j,1),3) = yaxel( rcur, phi(j), alfa)
         enddo
      enddo
      do j=1, ja1, 2
         npos(noda(ima,j,1),1) = xaxel( rxy4(j), phi(j), alfa)
         npos(noda(ima,j,1),2) = 0.0
         npos(noda(ima,j,1),3) = t
      enddo
      do j=ja1+2, jma, 2
         npos(noda(ima,j,1),1) = xaxel( rxy4(j), phi(j), alfa)
         npos(noda(ima,j,1),2) = 0.0
         npos(noda(ima,j,1),3) = yaxel( rxy4(j), phi(j), alfa)
      enddo
c
c  Node coordinates in the region embrazing Zone S (the crack zone)
c
      do k=3, ka1, 2
         do j=1, jma-4, 2
            do i=1, 3, 2
               npos(noda(i,j,k),1) = npos(noda(i,j,1),1)
               npos(noda(i,j,k),2) = yl((k-1)/2)
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
               npos(noda(i,j,k),2) = yl((k-1)/2)
               npos(noda(i,j,k),3) = npos(noda(i,j,1),3)
            enddo
         enddo
         do j=1, jma, 2
            do i=ia2+2, ima, 2
               npos(noda(i,j,k),1) = npos(noda(i,j,1),1)
               npos(noda(i,j,k),2) = yl((k-1)/2)
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
      y = yl((ka2-1)/2)
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
c
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
      eps = 1.0d-04
      do i=7, ima-2, 2
         psi = real(i-1)/real(ima-1)
         npos(noda(i,1,ka2),1) = 0.0
         npos(noda(i,1,ka2),2) = y
         npos(noda(i,1,ka2),3) = psi * t
         psi = real(i-5)/real(ima-5)
         do j=3, jma-4, 2
            call mapping_a(alfa,phi(j),rcur,rho1,rxy4(j),psi,eps)
            npos(noda(i,j,ka2),1) = xaxel( rcur, phi(j), alfa)
            npos(noda(i,j,ka2),2) = y
            npos(noda(i,j,ka2),3) = yaxel( rcur, phi(j), alfa)
         enddo
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
      y1 = yl((ka1-1)/2)
      y2 = yl((ka2-1)/2)
      do k=ka1+2, ka2-2, 2
         y = yl((k-1)/2)
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
            npos(noda(i,j,1),2) = npos(nodc(1,1,3),2)
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
               npos(noda(i,j,k),1) = 0.25d0 *
     &               ( npos(noda(i-1,j,k),1) + npos(noda(i+1,j,k),1) +
     &                 npos(noda(i,j-1,k),1) + npos(noda(i,j+1,k),1) )
               npos(noda(i,j,k),2) = 0.25d0 *
     &               ( npos(noda(i-1,j,k),2) + npos(noda(i+1,j,k),2) +
     &                 npos(noda(i,j-1,k),2) + npos(noda(i,j+1,k),2) )
               npos(noda(i,j,k),3) = 0.25d0 *
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
               npos(noda(i,j,k),2)=yl((k-1)/2)
               npos(noda(i,j,k),3)=npos(noda(i,j,ka2),3)
c 
               npos(noda(i,j,k-1),1)=npos(noda(i,j,ka2),1)
               npos(noda(i,j,k-1),2)=(yl((k-1)/2)+yl((k-3)/2))/2.d0
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
            npos(noda(i,1,k),1)=0.d0
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
c      if (ima.gt.0) stop

c
c                   ==========================
c                         Z O N E    B
c                   ==========================
c
      write(*,'(t20,a)') 'Zone B'
c
c  Element grading: Z = 0
c
      if (mbtype.eq.0) then
         dx = npos( noda(ima,jma,1),1 ) - npos( noda(ima-2,jma,1),1 )
         x = w - npos(noda(ima-2,jma,1),1)
         lambda = calc_lambda(dx,x,mb+1)
         if (exp(lambda).lt.1) then
            dx = (w - npos(noda(ima,jma,1),1)) / real(mb)
            lambda = 0.0
         endif
      else
         dx = ( w - npos(noda(ima,jma,1),1) ) * (1.d0-mb_bias) /
     &        ( mb_bias*( 1.d0 - mb_bias**dble(mb) ) )
         lambda=dlog(mb_bias)
      endif
      it=0
      do i=3, imb, 2
         it=it+1
         npos(nodb(i,jmb,1),1) = npos(nodb(i-2,jmb,1),1)+
     &              dx * dexp(lambda*dble(it))
         npos(nodb(i,jmb,1),2)=0.0
         npos(nodb(i,jmb,1),3)=0.0
      enddo
      npos(nodb(imb,jmb,1),1) = w
c
c Generate node coordinates in region  1 =< k =< ka2 :
c
      do k=1, ka2, 2
         i = imb
         do j=1, jmb
            psi = real(j-1)/real(jmb-1)
            npos(nodb(i,j,k),1) = w
            npos(nodb(i,j,k),2) = yl((k-1)/2)
            npos(nodb(i,j,k),3) = (1.-psi)*t
         enddo
         npos(nodb(i,jmb,k),3) = 0.0
      enddo
c
      chi = 2.0
      do k=1, ka2, 2
         l1 = npos(nodb(1,jmb,k),1)
         y = yl( (k-1)/2 )
         do j=1, jmb, 2
            z1 = npos(nodb(  1,j,k),3)
            z2 = npos(nodb(imb,j,k),3)
            x1 = npos(nodb(  1,j,k),1)
            x2 = npos(nodb(imb,j,k),1)
            do i=3, imb-2, 2
               x = npos(nodb(i,jmb,1),1)
               psi = (x-l1) / (w-l1)
               npos(nodb(i,j,k),1) = (1.-psi)*x1 + psi*x2
               npos(nodb(i,j,k),1) = npos(nodb(i,j,k),1)*(1.-psi)**chi
     &                             + x * (1.-(1.-psi)**chi)
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
               npos(nodb(i,j,k),1) = 0.25d0 *
     &            ( npos(nodb(i-1,j,k),1) + npos(nodb(i+1,j,k),1) +
     &              npos(nodb(i,j-1,k),1) + npos(nodb(i,j+1,k),1) )
               npos(nodb(i,j,k),2) = 0.25d0 *
     &            ( npos(nodb(i-1,j,k),2) + npos(nodb(i+1,j,k),2) +
     &              npos(nodb(i,j-1,k),2) + npos(nodb(i,j+1,k),2) )
               npos(nodb(i,j,k),3) = 0.25d0 *
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
               npos(nodb(i,j,k),2) = yl((k-1)/2)
               npos(nodb(i,j,k),3) = npos(nodb(i,j,ka2),3)
c
               npos(nodb(i,j,k-1),1) = npos(nodb(i,j,ka2),1)
               npos(nodb(i,j,k-1),2) = (yl((k-1)/2)+yl((k-3)/2))/2.d0
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
c      if (imb.gt.0) stop
c
c      call plot_zone_c(nodc,nc1,nc2,nc3)
c      call plot_zone_s(nods,ns1,ns2,ns3,rxy1,rxy2(1),fi,alfa_y)
c      call plot_zone_a(noda,na1,na2,na3)
c      call plot_zone_b(nodb,nb1,nb2,nb3)
c 
c-----------------------------------------------------------!
c  Coordinate transformation: z(y)
c-----------------------------------------------------------!
c
c.... Determine the transformation coefficient. Note tha the trans-
c     formation only will occur in the Z-direction (3-direction).
c
      call transformation_coefficient(trc,nt1,nt2,t,yl,z,lt)
c
      do i=0, inm
         nnr(i)=0
      enddo
c
c...Coordinate transformation in Zone C:
      do j=1, jmc
         do k=1, kmc
            do i=1, imc
               if ((nodc(i,j,k).gt.0).and.
     &            (nnr(nodc(i,j,k)).eq.0)) then
                  call coord_trans(nodc,nc1,nc2,nc3,i,j,k,trc,nt1,nt2)
                  nnr(nodc(i,j,k))=1
               endif
            enddo
         enddo
      enddo
c
c...Coordinate transformation in Zone S:
      do j=1, jms
         do k=1, kms
            do i=1, ims
               if ((nods(i,j,k).gt.0).and.
     &            (nnr(nods(i,j,k)).eq.0)) then
                  call coord_trans(nods,ns1,ns2,ns3,i,j,k,trc,nt1,nt2)
                  nnr(nods(i,j,k))=1
               endif
            enddo
         enddo
      enddo
c
C...Coordinate transformation in Zone A & B:
      do k=1, kma
         do j=1, jma
            do i=1, ima
               if ((noda(i,j,k).gt.0).and.(nnr(noda(i,j,k)).eq.0)) then
                  call coord_trans(noda,na1,na2,na3,i,j,k,trc,nt1,nt2)
                  nnr(noda(i,j,k))=1
               endif
            enddo
         enddo
         do j=1, jmb
            do i=1, imb
               if ((nodb(i,j,k).gt.0).and.(nnr(nodb(i,j,k)).eq.0)) then
                  call coord_trans(nodb,nb1,nb2,nb3,i,j,k,trc,nt1,nt2)
                  nnr(nodb(i,j,k))=1
               endif
            enddo
         enddo
      enddo
      nodnum=0
      do i=1, inm
         if(nnr(i).gt.0) nodnum=nodnum+1
      enddo
      write(*,'(t15,a,i6,a/)')  '=> ', nodnum,
     &                ' number of nodes has been defined !'
 
C...Generate plot files
      if (pl.eq.1) then
        call create_plotfil(nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3,
     &    noda,na1,na2,na3, nodb,nb1,nb2,nb3, job,jobh)
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
