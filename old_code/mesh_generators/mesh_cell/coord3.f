c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function xaxel(r,fi,alfa)
      double precision r,fi,alfa
      xaxel=(r+1.d0/r)*alfa*dcosd(fi)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function yaxel(r,fi,alfa)
      double precision r,fi,alfa
      yaxel=(r-1.d0/r)*alfa*dsind(fi)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision  function rhox(x,fi,alfa)
      double precision  x,fi,alfa
      rhox=(x+dsqrt(x**2.d0-4.d0*(alfa*dcosd(fi))**2.d0))/
     &          (2.d0*alfa*dcosd(fi))
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function rhoy(y,fi,alfa)
      double precision   y,fi,alfa
      rhoy=(y+dsqrt(y**2.d0+4.d0*(alfa*dsind(fi))**2.d0))/
     &           (2.d0*alfa*dsind(fi))
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function rhoxy(x,y,fi,alfa)
      double precision   x,y,fi,alfa,fic
      double precision   rhox,rhoy
      parameter (fic = 45.0)
      if (fi.le.fic) then
         rhoxy = rhox(x,fi,alfa)
      else
         rhoxy = rhoy(y,fi,alfa)
      endif
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
	double precision function rho_ellips(a,b,fi,alfa)
	   double precision fi,a,b,alfa,s1,c1,k
	   c1=dcosd(fi)
	   s1=dsind(fi)
	   k=( 2.d0*( (b*c1)**2.d0-(a*s1)**2.d0 )-(a*b/alfa)**2.d0 )/
     &         ((b*c1)**2.d0+(a*s1)**2.d0)
	   rho_ellips=dsqrt( (-k+dsqrt(k**2.d0-4.d0))/2.d0 )
	   return
	end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function calc_lambda(l1,ln,n)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  Rutinen ber{knar lambda. Parametern lambda anv{nds vid
c  gradering av n{tet i olika omr}den.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
      double precision eps
      parameter (eps=1.e-06)
      integer  n,i
      double precision l1,ln,lambda,xtot,dlambda,sum,sumold
      xtot=ln/l1
      lambda=0.015d0
      dlambda=0.005d0
      sum=0
 10   continue
        lambda=lambda+dlambda
        sumold=sum
        sum=0
        do i=1, n
           sum=sum+exp(lambda*dble(i-1))
        enddo
      if (abs(sum/xtot-1.).gt.eps) then
         if (((xtot-sum)*(xtot-sumold)).lt.0)  dlambda=-dlambda/2.
         goto 10
      endif
      calc_lambda=lambda
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      double precision function hyperbel_length(r1,r2,fi,alfa)
C----------------------------------------------------------------C
C  Rutinen ber{knar en str{cka utmed en hyperbel.                C
C   => Hyperbelns ekv (x/a)**2 - (z/b)**2 = 1                    C
C----------------------------------------------------------------C
	double precision r1,r2,fi,alfa,a,b,dr
	a=alfa*dcosd(fi)
	b=alfa*dsind(fi)
	dr=r2-r1
	sum=dsqrt((a*(1.d0-1.d0/r1**2.d0))**2.d0+
     &             (b*(1.d0+1.d0/r1**2.d0))**2.d0)
C . . . Dela in intervallet i 50 st lika delar. Anvand Simpsons rule.
	do i=2, 50, 2
	   r=r1+dr*dble(i-1)/50.d0
	   sum=sum+4.d0*dsqrt((a*(1.d0-1.d0/r**2.d0))**2.d0+
     &            (b*(1.d0+1.d0/r**2.d0))**2.d0)
	enddo
	do i=3, 49, 2
	   r=r1+dr*dble(i-1)/50.d0
	   sum=sum+2.d0*dsqrt((a*(1.d0-1.d0/r**2.d0))**2.d0+
     &               (b*(1.d0+1.d0/r**2.d0))**2.d0)
	enddo
	sum=sum+dsqrt((a*(1.d0-1.d0/r2**2.d0))**2.d0+
     &           (b*(1.d0+1.d0/r2**2.d0))**2.d0)
	hyperbel_length=sum*dr/(150.d0)
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine transformation_coefficient(trc,n1,n2,t,y,z,lt)
C--- The routine determines the transformation coefficients
      implicit none
      integer  n1,n2,lt,i
      double precision trc(n1,n2),t,y(0:200),z(0:200,2),a(0:1),b(0:1)
      do i=1, lt
         a(1) = ( z(i,1)-z(i-1,1) ) / ( y(i)-y(i-1) )
         a(0) = z(i,1) - a(1)*y(i)
         b(1) = ( z(i,2)-z(i-1,2) ) / ( y(i)-y(i-1) )
         b(0) = z(i,2) - b(1)*y(i)
         trc(i,1) = a(0)
         trc(i,2) = a(1)
         trc(i,3) = ( b(0)-a(0) ) / t
         trc(i,4) = ( b(1)-a(1) ) / t
         trc(i,5) = y(i) + ( y(i)-y(i-1) ) * 1.d-5
      enddo
      trc(lt+1,5) = 2.d0*y(lt)
      return
      end
c
C----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine coord_trans(nod,n1,n2,n3,i,j,k,trc,n4,n5)
C--The routine performs cordinate transformation
      implicit none
c
      include 'common_nod.f'
c
      integer  n1,n2,n3,n4,n5,nod(n1,n2,n3),i,j,k,m
      double precision   trc(n4,n5)
      m=0
10    m=m+1
      if (npos(nod(i,j,k),2).gt.trc(m,5)) goto 10
      npos(nod(i,j,k),3) = trc(m,1) + trc(m,2)*npos(nod(i,j,k),2) +
     &           trc(m,3)*npos(nod(i,j,k),3) +
     &           trc(m,4)*npos(nod(i,j,k),2)*npos(nod(i,j,k),3)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine angular_distribution(dcell,eta_t1,eta_t2,phi)
c
      implicit none
c
      double precision  ellipse_l
      external          ellipse_l
c
      double precision pi,eps
C     parameter ( pi = 4.*atan(1.), eps=1.e-7 )
      parameter ( pi = 3.141592654, eps=1.e-7 )
c
      integer i,j
      double precision dcell,eta_t1,eta_t2,phi(*), phi_tmp(1000),
     &       rel,rel_old,psi,dphi,phi1,phi2,dstar,dcur
c
      integer       ncd1,ncd2,ncd3
      common /cell/ ncd1,ncd2,ncd3
      integer      imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
      phi1 = 0.0
      do i=1, ncd1
         dcur = dcell*eta_t1
         rel_old = 1.
         dphi = 0.002*pi
         phi2 = phi1
 10      phi2 = phi2 + dphi
         call qromb(ellipse_l,phi1,phi2,dstar,1.0d-6)
         rel = (dcur-dstar)/dcur
         if (abs(rel).gt.eps) then
            if ( (rel_old*rel).lt.0.) dphi = -dphi/2.
            rel_old = rel
            goto 10
         endif
         phi_tmp(i) = (phi2/pi)*180.0
         phi1=phi2
      enddo     
c
      do i=ncd1+1, ncd1+ncd2
         psi = real(i-ncd1)/real(ncd2+1)
         dcur = dcell*( (1.-psi)*eta_t1 + psi*eta_t2 )
         rel_old = 1.
         dphi = 0.002*pi
         phi2 = phi1
 20      phi2 = phi2 + dphi
         call qromb(ellipse_l,phi1,phi2,dstar,1.0d-6)
         rel = (dcur-dstar)/dcur
         if (abs(rel).gt.eps) then
            if ( (rel_old*rel).lt.0.) dphi = -dphi/2.
            rel_old = rel
            goto 20
         endif
         phi_tmp(i) = (phi2/pi)*180.
         phi1=phi2
      enddo     
c
      do i=ncd1+ncd2+1, ncd1+ncd2+ncd3
         dcur = dcell*eta_t2
         rel_old = 1.
         dphi = 0.002*pi
         phi2 = phi1
 30      phi2 = phi2 + dphi
         call qromb(ellipse_l,phi1,phi2,dstar,1.0d-6)
         rel = (dcur-dstar) / dcur
         if (abs(rel).gt.eps) then
            if ( (rel_old*rel).lt.0.) dphi = -dphi/2.
            rel_old = rel
            goto 30
         endif
         phi_tmp(i) = (phi2/pi)*180.
         phi1=phi2
      enddo
c
      j = 1
      phi(j) = 90.0
      do i = ncd1+ncd2+ncd3-1, 1, -1
         j = j + 2
         phi(j) = phi_tmp(i)
         phi(j-1) = ( phi(j)+phi(j-2) ) / 2.
      enddo
      j = j + 2
      phi(j) = 0.0
      phi(j-1) = ( phi(j)+phi(j-2) ) / 2.
      if (j.ne.jmc) then
         write(*,*) '>> Angular distribution erroneous!'
         stop
      endif
c
      return
      end

c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine mapping_s(alfa,phi,icase,r,r1,r2,psi,sl,eps)
c
      implicit none
c
      double precision  hyperbel_l
      external          hyperbel_l
c
c      double precision eps
c      parameter ( eps=1.e-3 )
c
      integer     icase,it,itmax
      parameter ( itmax = 40 )
      logical map
      double precision alfa,phi,r,r1,r2,psi,sl,eps
      double precision s,ds,dstar,drho,rho2,rel,rel_old,eps1
c
      double precision  a1,b1
      common /hyperbel/ a1,b1
c
      a1 = alfa*cosd(phi)
      b1 = alfa*sind(phi)
      it = 0
      eps1 = eps/10.0
c
      if (icase.eq.1) then
c
c Keep the physical length ahead of the crack
         map=.false.
         ds = sl
         drho = (r2-r1)*(psi-1.)
         rho2 = r2 + drho
         drho = 0.10*drho
c         call qromb(hyperbel_l,r2,rho2,dstar,eps1)
         rel_old = 1.0
c         if ((dstar-ds).lt.0) rel_old = -rel_old
c
c         drho = 0.50*(r2-r1)*(psi-1.)
c         rho2 = r2 - drho/2.
c         call qromb(hyperbel_l,r2,rho2,dstar,eps1)
c         rel_old = 1.0
c         if ((ds-dstar).lt.0) rel_old = -rel_old
 10      continue
            it = it + 1
            rho2 = rho2 + drho
            call qromb(hyperbel_l,r2,rho2,dstar,eps1)
            rel = (ds-dstar)/ds
            
         if (it.gt.itmax) then
            write(*,110) '* it >',itmax,' in subr. mapping_s(...)'
            goto 30
         elseif (abs(rel).gt.eps) then
            if ( (rel_old*rel).lt.0.d0) drho = -drho/2.
            rel_old = rel
            goto 10
         endif
         r = rho2
c
      else
c
c Conformal mapping together with linear scaling
         map=.true.
         call qromb(hyperbel_l,r1,r2,s,eps1/10.)
         ds = (1.0-psi)*s
c
c
         drho = (1.-psi)*(r2-r1)
         rho2 = r2 - drho
         drho = -0.10*drho
         rel_old = 1.0
c         call qromb(hyperbel_l,rho2,r2,dstar,eps1)
c         if (dstar.lt.ds) rel_old = -rel_old
c
c       Find a new rho
c         drho = 0.50*(1.0-psi)*(r2-r1)
c         rho2 = r2 - drho
c         rel_old = -1.0
  20     continue
            it = it + 1
            rho2 = rho2 + drho
            call qromb(hyperbel_l,rho2,r2,dstar,eps1)
            rel = (ds-dstar)/ds
         if (it.gt.itmax) then
            write(*,110) '* it >',itmax,' in subr. mapping_s(...)'
            goto 30
         elseif (abs(rel).gt.eps) then
            if ( (rel_old*rel).lt.0.d0) drho = -drho/2.
            rel_old = rel
            goto 20
         endif
         r = rho2
c
      endif
c
  30  continue
c      if (map) then
c        write(26,111) ' Con mapping it = ',it,' psi=',psi
c      else
c        write(26,111) ' Phys length it = ',it,' psi=',psi
c      endif
 111  format(t1,a,i4,a,g15.6)
 110  format(a,i3,a)
c
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      subroutine mapping_a(alfa,phi,r,r1,r2,psi,eps)
c
      implicit none
c
      double precision  hyperbel_l
      external          hyperbel_l
c
c      double precision eps
c      parameter ( eps=1.e-3 )
c
      integer     it,itmax
      parameter ( itmax = 40 )
      double precision alfa,phi,r,r1,r2,psi,eps,eps1
      double precision s,ds,dstar,drho,rho2,rel,rel_old
c
      double precision  a1,b1
      common /hyperbel/ a1,b1
c
      a1 = alfa*cosd(phi)
      b1 = alfa*sind(phi)
      it = 0
c
      eps1 = eps/10.0
      call qromb(hyperbel_l,r1,r2,s,eps1)
      ds = psi*s
c
c  Find a new rho
c      drho = 0.50*psi*(r2-r1)
c      rho2 = r1 - drho/2.
c      call qromb(hyperbel_l,r1,rho2,dstar,eps1)
c
c      rho2 = r1 - drho/2.
c      rel_old = 1.
c
      drho = psi*(r2-r1)
      rho2 = r1 + drho
      drho = 0.10*drho
      rel_old = 1.0
c      call qromb(hyperbel_l,rho2,r2,dstar,eps1)
c      if (dstar.gt.ds) rel_old = -rel_old
c
  10  continue
         it = it + 1
         rho2 = rho2 + drho
         call qromb(hyperbel_l,r1,rho2,dstar,eps1)
         rel = (ds-dstar)/ds
      if (it.gt.itmax) then
         write(*,110) '* it >',itmax,' in subr. mapping_a(...)'
         goto 20
      elseif (abs(rel).gt.eps) then
         if ( (rel_old*rel).lt.0.) drho = -drho/2.
         rel_old = rel
         goto 10
      endif
  20  continue
      r = rho2
c      write(27,'(t1,a,i4,a,g15.6)') ' it = ',it,' psi=',psi
 110  format(a,i3,a)
c
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      subroutine div_hyperbel(r1,r2,f1,f2,n,l1,alfa,xy,rf)
c
      implicit none
      integer  n,i,j,m
      double precision r1,r2,f1,f2,l1,alfa,xy(2,*),rf(2,*),
     &           psi,r,f,dx,x1,x2,dy,y1,y2,ds,ln,lambda,xl(100),s
      double precision xaxel,yaxel,calc_lambda
c
c Compute the total length of the line
c
      m=200
      ln = 0.
      x1 = xaxel(r1,f1,alfa)
      y1 = yaxel(r1,f1,alfa)
      do i=1, m
         psi=real(i-1)/real(m-1)
         r = (1.-psi)*r1 + psi*r2
         f = (1.-psi)*f1 + psi*f2
         x2 = xaxel(r,f,alfa)
         y2 = yaxel(r,f,alfa)
         dx = x2-x1
         dy = y2-y1
         ds = sqrt(dx*dx+dy*dy)
         ln = ln + ds
         x1 = x2
         y1 = y2
      enddo
c
      lambda = calc_lambda(l1,ln,n)
c         write(31,*) '>> Ln = ',real(ln),' l1=',real(l1),
c     &      ' lambda=',real(lambda),' exp(lambda)=',real(exp(lambda))
      if (exp(lambda).le.1.) then
         do i=1, n
            xl(i) = ln*real(i)/real(n)
         enddo
      else
         xl(1) = l1
         do i=2, n
            xl(i) = xl(i-1) + l1*exp(lambda*real(i-1))
         enddo
      endif
c
      m=501
      s = 0
      j = 1
      x1 = xaxel(r1,f1,alfa)
      y1 = yaxel(r1,f1,alfa)
      do i=1, m
         psi=real(i-1)/real(m-1)
         r = (1.-psi)*r1 + psi*r2
         f = (1.-psi)*f1 + psi*f2
         x2 = xaxel(r,f,alfa)
         y2 = yaxel(r,f,alfa)
         dx = x2-x1
         dy = y2-y1
         ds = sqrt(dx*dx+dy*dy)
         s = s + ds
         x1 = x2
         y1 = y2
         if ( s.ge.xl(j) ) then
            xy(1,j) = xaxel(r,f,alfa)
            xy(2,j) = yaxel(r,f,alfa)
            rf(1,j) = r
            rf(2,j) = f
            j = j + 1
         endif
      enddo
c
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      subroutine div_linear(r1,r2,f1,f2,n,l1,alfa,xy)
c
      implicit none
      integer  n,i
      double precision r1,r2,f1,f2,l1,alfa,xy(2,*),
     &                 x1,x2,y1,y2,ln,lambda,xl(100),norm_x,norm_y
      double precision xaxel,yaxel,calc_lambda
c
c Compute the total length of the line
c
      x1 = xaxel(r1,f1,alfa)
      y1 = yaxel(r1,f1,alfa)
      x2 = xaxel(r2,f2,alfa)
      y2 = yaxel(r2,f2,alfa)
      ln = sqrt( (x2-x1)**2. + (y2-y1)**2. )
c
c calc. norm vector
c
      norm_x = (x2-x1)/ln
      norm_y = (y2-y1)/ln
c
      lambda = calc_lambda(l1,ln,n)
      if (exp(lambda).le.1.) then
         do i=1, n
            xl(i) = ln*real(i)/real(n)
         enddo
      else
         xl(1) = l1
         do i=2, n
            xl(i) = xl(i-1) + l1*exp(lambda*real(i-1))
         enddo
      endif
      do i=1, n
         xy(1,i) = x1 + xl(i)*norm_x
         xy(2,i) = y1 + xl(i)*norm_y
      enddo
c
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      double precision function rho_inc(alfa,phi,ds,rho1)
c
      implicit none
c
      double precision  eps
      parameter ( eps=1.e-6 )
c
      double precision alfa,phi,ds,rho1
      double precision rel_old,rel,drho,rho2,dx,dy,ds_new,xsign
      double precision xaxel,yaxel
c----67
      rel_old = 1.0
      if (ds.gt.0) then
         drho =  0.001
         xsign = 1.0
      else
         drho = -0.001
         xsign = -1.0
      endif
      rho2 = rho1
 10   continue
         rho2 = rho2 + drho
         dx = xaxel(rho2,phi,alfa) - xaxel(rho1,phi,alfa)
         dy = yaxel(rho2,phi,alfa) - yaxel(rho1,phi,alfa)
         ds_new = xsign*sqrt( dx*dx + dy*dy )
         rel = ( ds - ds_new ) / ds
         if (abs(rel).gt.eps) then
            if ( (rel_old*rel).lt.0.d0) drho = -drho/2.
            rel_old = rel
            goto 10
         endif
      continue
      rho_inc = rho2 - rho1
c
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      subroutine calc_rf(x,y,alfa,rho,phi)
c-----------------------------------------------------------------
c Rutinen räknar ut Rho och PHI för X och Y är givna
c-----------------------------------------------------------------
      implicit none
      double precision  x,y,alfa,rho,phi,xn,yn,xy_p,xy_m,t1,t2,
     &                  f,fprim,eps
      eps=1.d-10
c... FIND RHO FIRST !
      xn=x/alfa
      yn=y/alfa
      xy_p = xn*xn+yn*yn
      xy_m = xn*xn-yn*yn
c... Find the root (Newton Raphson's method):
      t1 =  xy_p + 2.d0
 10   continue
         f=1.d0+t1*( -xy_p + t1*( 2.d0*(xy_m-1.d0) + t1*(t1-xy_p) ) )
         fprim=(-xy_p+t1*(4.d0*(xy_m-1.d0)+t1*(4.d0*t1-3.d0*xy_p) ))
         t2 = t1 - f / fprim
         if ( abs((t2-t1)/t2).lt.eps) then
            goto 20
         else
           t1 = t2
           goto 10
         endif
20    continue
      if (t2.lt.1.d0) then
         t1=1.d0/t2
         goto 10
      endif
      rho=dsqrt(t2)
c     phi = 0.5d0*dacosd((xyn-t2-1.d0/t2)*0.5d0)
      phi=dacosd(xn/(rho+1.d0/rho))
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c

c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function ellipse_l(phi)
c
      implicit none
      double precision phi,a,b,f
      common /ellipse/ a,b
c
      f = sqrt( (a*sin(phi))**2.0 + (b*cos(phi))**2.0 )
      ellipse_l = f
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      double precision function hyperbel_l(rho)
c
      implicit none
      double precision rho,a,b,f
      common /hyperbel/ a,b
c
      f = sqrt((a*(1.-1./(rho*rho)))**2.0+(b*(1.+1./(rho*rho)))**2.0)
      hyperbel_l = f
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine qromb(func,a,b,ss,eps)
      implicit none
c
      double precision a,b,ss,func
      external func
c
      integer  jmax,jmaxp,k,km
      double precision eps
c      parameter (eps=1.d-5, jmax=20, jmaxp=jmax+1, k=5, km=k-1)
      parameter (jmax=20, jmaxp=jmax+1, k=5, km=k-1)
c
      integer  j
      double precision dss, h(jmaxp), s(jmaxp)
c
      h(1) = 1.d0
      do j=1, jmax
         call trapzd(func,a,b,s(j),j)
         if (j.ge.k) then
            call polint(h(j-km), s(j-km), k, 0.d0, ss, dss)
            if (abs(dss).le.eps*abs(ss)) return
         endif
         s(j+1)=s(j)
         h(j+1)=0.25d0*h(j)
      enddo
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      subroutine trapzd(func, a, b, s, n)
      implicit none
c
      integer  n
      double precision a,b,s,func
      external func
c
      integer it,j
      double precision del, sum, tnm, x
c
      if (n.eq.1) then
         s = 0.5d0*(b-a)*(func(a)+func(b))
      else
         it = 2**(n-2)
         tnm=dble(it)
         del=(b-a)/tnm
         x = a + 0.5d0*del
         sum = 0.d0
         do j=1, it
            sum = sum + func(x)
            x = x + del
         enddo
         s = 0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      subroutine polint(xa, ya, n, x, y, dy)
      implicit none
c
      integer  n,nmax
      double precision dy,x,y,xa(n),ya(n)
      parameter  (nmax=10)
c
      integer  i,m,ns
      double precision den,dif,dift,ho,hp,w,c(nmax),d(nmax)
c
      ns = 1
      dif = abs(x-xa(1))
      do i=1, n
         dift=abs(x-xa(i))
         if (dift.lt.dif) then
            ns=i
            dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
      enddo
      y = ya(ns)
      ns = ns - 1
      do m=1, n-1
         do i=1, n-m
            ho = xa(i)-x
            hp = xa(i+m)-x
            w  = c(i+1)-d(i)
            den= ho - hp
            if (den.eq.0.d0) write(*,*) '>> Failure in subr. polint!'
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
         enddo
         if (2*ns.lt.(n-m)) then
            dy = c(ns+1)
         else
            dy = d(ns)
            ns = ns-1
         endif
         y = y+dy
      enddo
      return
      end

c
c----67--1---------2---------3---------4---------5---------6---------712
c
