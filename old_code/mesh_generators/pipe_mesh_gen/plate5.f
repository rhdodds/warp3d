c----67--1---------2---------3---------4---------5---------6---------712
c
c    File: plate5.f   latest modification July 11, 1997
c    "write(*,...)" changed to "write(iws,...)", 7/29/97
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function xaxel(r,fi,alfa)
      double precision r,fi,alfa
      xaxel=(r+1.0/r)*alfa*cosd(fi)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function yaxel(r,fi,alfa)
      double precision r,fi,alfa
      yaxel=(r-1.0/r)*alfa*sind(fi)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision  function rhox(x,fi,alfa)
      double precision  x,fi,alfa
      rhox=(x+sqrt(x**2.0 - 4.0*(alfa*cosd(fi))**2.0))/
     &     (2.0*alfa*cosd(fi))
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function rhoy(y,fi,alfa)
      double precision   y,fi,alfa
      rhoy=(y+sqrt(y**2.0 + 4.0*(alfa*sind(fi))**2.0))/
     &     (2.0*alfa*sind(fi))
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
      c1=cosd(fi)
      s1=sind(fi)
      k = ( 2.0*( (b*c1)**2.0 - (a*s1)**2.0 ) - (a*b/alfa)**2.0 )/
     &        ( (b*c1)**2.0 + (a*s1)**2.0 )
      rho_ellips = sqrt( (-k + sqrt(k**2.0 - 4.0) ) / 2.0 )
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function calc_lambda(l1,ln,n)
c
c This routine computes the bias parameter which is used in
c mesh refinement/coarsening.
c
      implicit none
      double precision eps
      parameter (eps=1.e-06)
      integer  n,i
      double precision l1,ln,lambda,xtot,dlambda,sum,sumold,xi
      xtot    = ln / l1
      lambda  = 0.0150
      dlambda = 0.0050
      sum=0
 10   continue
        lambda = lambda + dlambda
        sumold = sum
        sum = 0.0
        do i=1, n
           xi = i-1
           sum = sum + exp( lambda*xi )
        enddo
      if ( abs( sum/xtot - 1.0) .gt. eps) then
         if ( ((xtot-sum)*(xtot-sumold) ) .lt. 0.0) dlambda=-dlambda/2.
         goto 10
      endif
      calc_lambda=lambda
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      double precision function hyperbel_length(r1,r2,fi,alfa)
c
c Computes the arc-length of a segment of an hyperbolic curve
c 
c   Hyperble Eq.   (x/a)**2 - (z/b)**2 = 1
c
      double precision r1,r2,fi,alfa,a,b,dr,xi
      a = alfa*cosd(fi)
      b = alfa*sind(fi)
      dr = r2-r1
      sum = sqrt( (a*(1.0-1.0/r1**2.0))**2.0 +
     &            (b*(1.0+1.0/r1**2.0))**2.0 )
c
c Divide the interval in 50 equidistant segments, use Simpson's rule
c
      do i=2, 50, 2
         xi = i-1
         r  = r1 + dr*xi / 50.0
         sum = sum + 4.0*sqrt( (a*(1.0-1.0/r**2.0))**2.0 +
     &                         (b*(1.0+1.0/r**2.0))**2.0 )
      enddo
      do i=3, 49, 2
         xi  = i-1
         r   = r1 + dr*xi / 50.0
         sum = sum + 2.0*sqrt( (a*(1.0-1.0/r**2.0))**2.0 +
     &                         (b*(1.0+1.0/r**2.0))**2.0 )
      enddo
      sum = sum + sqrt( (a*(1.0-1.0/r2**2.0))**2.0 +
     &                  (b*(1.0+1.0/r2**2.0))**2.0 )
      hyperbel_length = sum*dr / 150.0
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine angular_distribution(phi,a,c,alfa,jms,tcr,
     &           iws,errorFlag)
c
      implicit none
c
      double precision  ellipse_l
      external          ellipse_l
c
      double precision pi,eps,unit
      parameter ( eps=1.e-7, unit=1.0 )
c
      integer jms,tcr,i,j,ns,n1,n2,k1,k2,dk,dk1,iws,errorFlag
c
      double precision eta
c      parameter (eta=0.40)
      double precision beta1,beta2
      parameter (beta1=1.8, beta2=0.5)
      double precision phi(200),a,c,alfa,rho0,tmp(200),rel,rel_old,psi,
     &                 dphi,phi1,phi2,s,ds(100),dstar,xi,bias1,bias2
c
      double precision a1,b1
      common /ellipse/ a1,b1
c
      pi = 4.*atan(unit)
      rho0 = sqrt( (c+a) / (c-a) )
      a1 = alfa*( rho0 + 1.0/rho0 )
      b1 = alfa*( rho0 - 1.0/rho0 )
c
c Determine the length of the segments along the crack front. 
c
      phi1 = 0.0
      phi2 = pi/2.
      call qromb(ellipse_l,phi1,phi2,s,eps, iws)
c
      ns = (jms-1)/2
c
c Embeddied crack or Semi elliptical crack - major axis free surface 
c
      if ((tcr .eq. 1) .or. (tcr.eq.2)) then
         eta = 0.50
         k1 = int(eta*ns+0.5)
         n1 = k1
         psi = 0.5**(a/c)
         ds(1)  = psi*(pi/2.)*(a/real(ns)) * (a/c)
         call calc_bias_int1(ds(1),s,ns,k1,bias1)
         do i=1, n1
            xi = i - 1
            ds(i) = ds(1)*bias1**xi
         enddo
         do i=n1+1, ns
            ds(i) = ds(n1)
         enddo 
      else
c
c Semi elliptical crack (minor axis free surface) or
c  1/4 elliptical crack (two free surfaces) 
c
         dk  = ns/4
         dk1 = dk*0.5*(a/c)**0.4
         k2  = ns/2 - (dk-dk1)
         k1  = ns/2 -  dk1
         n1 = k1
         n2 = ns - k2 + 1
         psi = 0.5**(a/c)
         ds(1)  = psi*(pi/2.)*(a/real(ns)) * (a/c)
         ds(ns) = 0.5*s/real(ns)
         call calc_bias_int2(ds(1),ds(ns),s,ns,k1,k2,bias1,bias2)
         do i=1, n1
            xi = i - 1
            ds(i) = ds(1)*bias1**xi
         enddo
         do i=ns, n2, -1
            xi = ns-i
            ds(i) = ds(ns)*bias2**xi
         enddo 
         do i=n1+1, n2-1
            ds(i) = ds(n1)
         enddo 
c
      endif
c
c Determine the anglular function
c
      phi1 = 0.0
      do i=1, ns
         rel_old = 1.
         dphi = 0.002*pi
         phi2 = phi1
 10      phi2 = phi2 + dphi
         call qromb(ellipse_l,phi1,phi2,dstar,eps, iws)
         rel = (ds(i)-dstar)/ds(i)
         if (abs(rel).gt.eps) then
            if ( (rel_old*rel).lt.0.) dphi = -dphi/2.
            rel_old = rel
            goto 10
         endif
         tmp(i) = (phi2/pi)*180.
         phi1=phi2
      enddo     
c
      j = 1
      phi(j) = 90.0
      do i = ns-1, 1, -1
         j = j + 2
         phi(j) = tmp(i)
         phi(j-1) = ( phi(j)+phi(j-2) ) / 2.
      enddo
      j = j + 2
      phi(j) = 0.0
      phi(j-1) = ( phi(j)+phi(j-2) ) / 2.
      if (j.ne.jms) then
         write(iws,*) '>> Angular distribution erroneous!'
         errorFlag = 1
         return
      endif
c
      return
      end
c

c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine calc_bias_int1(ds,s,n,k,bias)
c
c Computes bias to be used for k out of n elements 
c
      implicit none
      double precision eps
      parameter (eps=1.e-06)
      integer  n,k,i
      double precision ds,s,bias,dbias,sratio,sum,sumold,xi
      sratio = s/ds
      bias  = 0.8
      dbias = 0.1
      sum   = 0.0
 10   continue
        bias = bias + dbias
        sumold = sum
        sum = 1.0
        do i=2, k
           xi = i-1
           sum = sum + bias**xi
        enddo
        xi = k-1
        sum = sum + (n-k)*bias**xi
      if ( abs( sum/sratio - 1.0) .gt. eps) then
         if ( ((sratio-sum)*(sratio-sumold)) .lt. 0.0) dbias=-dbias/2.
         goto 10
      endif
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine calc_bias_int2(ds1,ds2,s,n,k1,k2,bias1,bias2)
c
c Computes bias to be used for k1, k2, elements respectivly.
c The total number of elements is n ( n > k1 + k2) 
c
      implicit none
      double precision eps
      parameter (eps=1.e-06)
      integer  n,k1,k2,i
      double precision ds1,ds2,s,bias1,bias2,b,db,sum,sumold,xi,xk1,xk2
      b   = 0.8
      db  = 0.1
      sum = 0.0
 10   continue
        b = b + db
        sumold = sum
        sum = ds1
        bias1 = b
        do i=2, k1
           xi = i-1
           sum = sum + ds1*bias1**xi
        enddo
        sum = sum + ds2
        xk1 = k1 - 1
        xk2 = k2 - 1
        bias2 = ((ds1/ds2)*bias1**xk1)**(1.0/xk2)
        do i=2, k2
           xi = i-1
           sum = sum + ds2*bias2**xi
        enddo
        xi = k1-1
        sum = sum + (n-k1-k2)*ds1*bias1**xi
      if ( abs( sum/s - 1.0) .gt. eps) then
         if ( ((s-sum)*(s-sumold)) .lt. 0.0) db = - db / 2.
         goto 10
      endif
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine mapping_s(alfa,phi,icase,r,r1,r2,psi,sl,eps,iws)
c
c Coordinates originally located on a plane are mapped on a single
c curved surface.
c The curvature of the surface are defined by a hyperble
c
      implicit none
c
      double precision  hyperbel_l
      external          hyperbel_l
c
      integer     icase,it,itmax,iws
      parameter ( itmax = 1000 )
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
         ds = sl
         drho = (r2-r1)*(psi-1.)
         rho2 = r2 + drho
         drho = 0.80*drho
         rel_old = 1.0
 10      continue
            it = it + 1
            rho2 = rho2 + drho
            call qromb(hyperbel_l,r2,rho2,dstar,eps1, iws)
            rel = (ds-dstar)/ds            
         if (it.gt.itmax) then
ccc            write(*,110) '* it >',itmax,'(1) in subr. mapping_s'
            write(iws,110) '* it >',itmax,'(1) in subr. mapping_s'
            goto 30
         elseif (abs(rel).gt.eps) then
            if ( (rel_old*rel).lt.0.d0) drho = -drho/2.
            rel_old = rel
            goto 10
         endif
         r = rho2
      else
c
c Conformal mapping in combined with linear scaling
c
         call qromb(hyperbel_l,r1,r2,s,eps1/10., iws)
         ds = (1.0-psi)*s
c
         drho = (1.-psi)*(r2-r1)
         rho2 = r2 - drho
         drho = -0.10*drho
         rel_old = 1.0
  20     continue
            it = it + 1
            rho2 = rho2 + drho
            call qromb(hyperbel_l,rho2,r2,dstar,eps1, iws)
            rel = (ds-dstar)/ds
         if (it.gt.itmax) then
ccc            write(*,110) '* it >',itmax,'(2) in subr. mapping_s'
            write(iws,110) '* it >',itmax,'(2) in subr. mapping_s'
            goto 30
         elseif (abs(rel).gt.eps) then
            if ( (rel_old*rel).lt.0.d0) drho = -drho/2.0
            rel_old = rel
            goto 20
         endif
         r = rho2
c
      endif
c
  30  continue
c
 110  format(a,i3,a)
c
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      subroutine mapping_a1(alfa,phi,r,r1,r2,psi,eps,iws)
c
c Coordinates originally located on a plane are mapped on
c a single curved surface. The curvature of the surface are
c defined by a hyperble
c
      implicit none
c
      double precision  hyperbel_l
      external          hyperbel_l
c
      integer     it,itmax,iws
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
      eps1 = 0.1*eps
      call qromb(hyperbel_l,r1,r2,s,eps1, iws)
      ds = psi*s
c
      drho = psi*(r2-r1)
      rho2 = r1 + drho
      drho = 0.10*drho
      rel_old = 1.0
c
  10  continue
         it = it + 1
         rho2 = rho2 + drho
         call qromb(hyperbel_l,r1,rho2,dstar,eps1, iws)
         rel = (ds-dstar)/ds
      if (it.gt.itmax) then
ccc         write(*,110) '* it >',itmax,' in subr. mapping_a(...)'
         write(iws,110) '* it >',itmax,' in subr. mapping_a(...)'
         goto 20
      elseif (abs(rel).gt.eps) then
         if ( (rel_old*rel).lt.0.) drho = -drho/2.
         rel_old = rel
         goto 10
      endif
  20  continue
      r = rho2
c
 110  format(a,i3,a)
c
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      subroutine mapping_a2(alfa,p1,p2,r1,r2,r,phi,psi,eps,iws)
c
c Coordinates originally located on a plane are mapped on
c a single curved surface. The curvature of the surface are
c defined in a elliptic-hyperbolic coordinate system.
c
      implicit none
c
      double precision  hyper_l
      external          hyper_l
c
      integer     it,itmax,iws
      parameter ( itmax = 100 )
      double precision alfa,p1,p2,r1,r2,r,phi,psi,eps,eps1
      double precision s,ds,dstar,drho,rho2,rel,rel_old
      double precision pi,phip,phi_fcn4
c
      double precision alpha,rr1,rr2,pp1,pp2
      common   /hyper/ alpha,rr1,rr2,pp1,pp2
c
      pi = 4.*atan(1.0)
      alpha = alfa
      rr1 = r1
      rr2 = r2
      pp1 = p1*pi/180.0
      pp2 = p2*pi/180.0
c
      eps1 = 0.1*eps
      call qromb(hyper_l,r1,r2,s,eps1, iws)
      ds = psi*s
c
      drho = psi*(r2-r1)
      rho2 = r1 + drho
      drho = 0.10*drho
      rel_old = 1.0
      it = 0
c
  10  continue
         it = it + 1
         rho2 = rho2 + drho
         call qromb(hyper_l,r1,rho2,dstar,eps1, iws)
         rel = (ds-dstar)/ds
      if (it.gt.itmax) then
ccc         write(*,110) '* it >',itmax,' in subr. mapping_a2(...)'
         write(iws,110) '* it >',itmax,' in subr. mapping_a2(...)'
         goto 20
      elseif (abs(rel).gt.eps) then
         if ( (rel_old*rel).lt.0.) drho = -drho/2.
         rel_old = rel
         goto 10
      endif
  20  continue
      r = rho2
      phi = phi_fcn4(phip,r,r1,r2,p1,p2)
c
 110  format(a,i3,a)
c
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      double precision function calc_amap_bias(eta,u1,um,v1,vm,m,
     &                          alfa,p1,p2,r1,r2,eps, iws)
      implicit none
      double precision hyper_l
      external         hyper_l
      integer          m, iws
      double precision eta,u1,um,v1,vm,alfa,p1,p2,r1,r2,eps,
     &                 t1,tn,s1,sn,beta,calc_lambda
c
      double precision alpha,rr1,rr2,pp1,pp2
      common   /hyper/ alpha,rr1,rr2,pp1,pp2
c
      alpha = alfa
      rr1 = r1
      rr2 = r2
      pp1 = p1
      pp2 = p2
      call qromb( hyper_l, r1, r2, sn, 0.1*eps, iws)
      t1 = (1.0-eta)*u1 + eta*v1
      tn = (1.0-eta)*um + eta*vm
      beta = 0.4
      s1 = t1*(sn/tn)**beta
      calc_amap_bias = exp( calc_lambda(s1,sn,m) )
c
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      double precision function calc_amap_length(alfa,p1,p2,r1,r2,
     &                          eps,iws)
      implicit none
      double precision hyper_l
      external         hyper_l
      double precision alfa,p1,p2,r1,r2,eps,sn,pi
      integer iws
c
      double precision alpha,rr1,rr2,pp1,pp2
      common   /hyper/ alpha,rr1,rr2,pp1,pp2
c
      pi = 4.*atan(1.0)
      alpha = alfa
      rr1 = r1
      rr2 = r2
      pp1 = p1*pi/180.0
      pp2 = p2*pi/180.0
      call qromb( hyper_l, r1, r2, sn, 0.1*eps, iws)
      calc_amap_length = sn
c
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      double precision function calc_amap_yp(a,b,y,psi,eps,iws)
      implicit none
      double precision ellipse_l
      external         ellipse_l
      integer     it,itmax,iws
      parameter ( itmax = 100 )
      double precision a,b,y,psi,eps,
     &                 pi,eps1,p1,p2,sn,ds,dstar,dp,rel,rel_old
c
      double precision a1,b1
      common /ellipse/ a1,b1
c
      pi = 4.*atan(1.0)
      a1 = a
      b1 = b
      p1 = 0.0
      p2 = asin(y/b)
      eps1 = 0.1*eps
      call qromb( ellipse_l, p1, p2, sn, eps1, iws)
c
      ds = psi*sn
c
      dp = psi*(p2-p1)
      p2 = p1 + dp
      dp = 0.20*dp
      rel_old = 1.0
      it = 0
c
  10  continue
         it = it + 1
         p2 = p2 + dp
         call qromb( ellipse_l, p1, p2, dstar, eps1, iws)
         rel = ( ds - dstar ) / ds
      if ( it .gt. itmax ) then
ccc         write(*,110) '* it >',itmax,' in fcn. calc_amap_yp( )'
         write(iws,110) '* it >',itmax,' in fcn. calc_amap_yp( )'
         goto 20
      elseif( abs(rel) .gt. eps ) then
         if ( (rel_old*rel) .lt. 0.0 ) dp = - dp / 2.0
         rel_old = rel
         goto 10
      endif
  20  continue
      calc_amap_yp = b*sin(p2)
c
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
      subroutine calc_rf(x,y,alfa,rho,phi)
c
c Coordinate transformation - calculates rho and phi
c for a given set of x and y
c
      implicit none
      double precision  x,y,alfa,rho,phi,xn,yn,xy_p,xy_m,t1,t2,
     &                  f,fprim,eps,one,four
      parameter (one=1.0, four=4.0)
      eps = 1.e-10
c Compute rho
      xn = x/alfa
      yn = y/alfa
      xy_p = xn*xn + yn*yn
      xy_m = xn*xn - yn*yn
c Find the root (Newton Raphson's method):
      t1 =  xy_p + 2.0
 10   continue
         f = one + t1*( -xy_p + t1*( 2.0*(xy_m - one) + t1*(t1-xy_p)))
         fprim=(-xy_p + t1*(four*(xy_m-one)+t1*(four*t1-3.0*xy_p) ))
         t2 = t1 - f / fprim
         if ( abs( (t2-t1) / t2 ) .lt. eps) then
            goto 20
         else
           t1 = t2
           goto 10
         endif
20    continue
      if (t2.lt.one) then
         t1=1.d0/t2
         goto 10
      endif
      rho=dsqrt(t2)
c phi = 0.5 * acosd((xyn-t2-1.d0/t2)*0.5)
      phi = acosd( xn / (rho+1.0/rho) )
      return
      end
c
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
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
c----67--1----*----2----*----3----*----4----*----5----*----6----*----712
c
      double precision function hyper_l(r)
c
      implicit none
      double precision r,psi,phi,phip,c,s,f
      double precision phi_fcn4
c
      double precision alpha,r1,r2,p1,p2
      common   /hyper/ alpha,r1,r2,p1,p2
c
      psi = (r-r1) / (r2-r1)
      phi = phi_fcn4(phip,r,r1,r2,p1,p2)
      c = cos(phi)
      s = sin(phi)
c
      f =   ( (1.0 - 1.0/(r*r))*c - (r + 1.0/r )*phip*s )**2.0
     &    + ( (1.0 + 1.0/(r*r))*s + (r - 1.0/r )*phip*c )**2.0
      f = alpha*sqrt(f)
      hyper_l = f
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function phi_fcn0(phip,r,r1,r2,p1,p2)
      implicit none
      double precision phi,phip,r,r1,r2,p1,p2,psi
      psi  = sqrt( (r-r1) / (r2-r1) )
      phi  =  p1*(1.0-psi) + p2*psi
c the derivative is singular at psi=0
      phip = (p2-p1) / (2.*sqrt(psi)*(r2-r1))
      phi_fcn0 = phi
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function phi_fcn1(phip,r,r1,r2,p1,p2)
      implicit none
      double precision phi,phip,r,r1,r2,p1,p2,psi
      psi  = (r-r1) / (r2-r1)
      phi  =  p1*(1.0-psi) + p2*psi
      phip = (p2-p1) / (r2-r1)
      phi_fcn1 = phi
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function phi_fcn2(phip,r,r1,r2,p1,p2)
      implicit none
      double precision phi,phip,r,r1,r2,p1,p2,psi
      psi  = (r-r1) / (r2-r1)
      phi  = p1*(1.0-psi*psi) + p2*psi*psi
      phip = ( 2.0*(p2-p1) / (r2-r1) )*psi
      phi_fcn2 = phi
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function phi_fcn3(phip,r,r1,r2,p1,p2)
      implicit none
      double precision phi,phip,r,r1,r2,p1,p2,psi,eta
      psi  = (r-r1) / (r2-r1)
      eta  = 0.5
      phi  = p1*( eta*(1.0-psi) + (1.0-eta)*(1.0-psi*psi) )
     &     + p2*( eta*(1.0-psi) + (1.0-eta)*(1.0-psi*psi) )
      phip = eta * (p2-p1)/(r2-r1)
     &    + (1.0-eta) * (2.*(p2-p1)/(r2-r1))*psi
      phi_fcn3 = phi
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function phi_fcn4(phip,r,r1,r2,p1,p2)
      implicit none
      double precision phi,phip,r,r1,r2,p1,p2,psi,pi2
      pi2 = 2.0*atan(1.0)
      psi  = sin( pi2 * (r-r1)/(r2-r1) )
      phi  = p1*(1.0-psi) + p2*psi
      phip = ( pi2*(p2-p1) / (r2-r1) )*cos(pi2*psi)
      phi_fcn4 = phi
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function phi_fcn5(phip,r,r1,r2,p1,p2)
      implicit none
      double precision phi,phip,r,r1,r2,p1,p2,psi,pi2
      pi2 = 2.0*atan(1.0)
      psi  = cos( pi2 * (r-r1)/(r2-r1) )
      phi  = p1*(1.0-psi) + p2*psi
      phip = ( pi2*(p2-p1) / (r2-r1) )*sin(pi2*psi)
      phi_fcn5 = phi
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
c Below, routines "inspired" from Numerical Recipies 
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine qromb(func,a,b,ss,eps,iws)
      implicit none
c
      double precision a,b,ss,func
      external func
c
      integer  jmax,jmaxp,k,km, iws
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
            call polint(h(j-km), s(j-km), k, 0.d0, ss, dss, iws)
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
      subroutine polint(xa, ya, n, x, y, dy, iws)
      implicit none
c
      integer  n,nmax,iws
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
ccc            if (den.eq.0.d0) write(*,*) '>> Failure in subr. polint!'
            if (den.eq.0.d0) write(iws,*) '>> Failure in subr. polint!'
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
