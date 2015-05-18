c----67--1---------2---------3---------4---------5---------6---------712
c
c    File: plate_wrptrc.f   Written July 2 1998
c                           Modified on July 6 1998
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_pre_traction(tge,tcr,tbc,rgeo,rload,ip,id,rd,
     &           ecrsur,ecrtip,erem,esid,no_of_nodes,iws,errorFlag)
c     &           ecrsur,ecsym,erem,esid,etop,no_of_nodes,iws)
      implicit none
      integer  tge,tcr,tbc,no_of_nodes,id(3,*),iws,errorFlag
      integer  ecrsur(2,*),ecrtip(2,*),erem(2,*),esid(2,*)
      double precision rgeo(*),rload(*),rd(3,*),ptol
      logical  ip
      parameter (ptol=1.e-10)
c
c  TBC = 2: Prescribed endtractions
c
      if (tbc.eq.2) then
         if (abs(rload(7)).gt.ptol) then
            ip =.true.
         else
            ip = .false.
         endif
         if (tge.eq.101) then
            call wrp_trac_tge101_2(tcr,rgeo,rload,id,rd,erem,
     &                            no_of_nodes,iws,errorFlag)
         elseif ( (tge.ge.201) .and. (tge.le.304) ) then
            call wrp_trac_tge200_2(tge,rgeo,rload,ip,id,rd,erem,esid,
     &                            no_of_nodes,iws,errorFlag)
         elseif ( (tge.ge.301) .and. (tge.le.308) ) then
            call wrp_trac_tge300_2(tge,rgeo,rload,ip,id,rd,erem,esid,
     &                             no_of_nodes,iws,errorFlag)
         endif
c
c  TBC = 3: Prescribed tractions on crack face
c
      elseif (tbc.eq.3) then
         ip = .false.
         if (tge.eq.101) then
            call wrp_trac_tge101_3(rload,id,rd,ecrsur,ecrtip,
     &                             no_of_nodes,iws,errorFlag)
         elseif ( (tge.ge.201) .and. (tge.le.204) ) then
            call wrp_trac_tge200_3(tge,rgeo,rload,id,rd,ecrsur,
     &                             no_of_nodes,iws,errorFlag)
         elseif ( (tge.ge.301) .and. (tge.le.308) ) then
            call wrp_trac_tge300_3(tge,rgeo,rload,id,rd,ecrsur,
     &                             no_of_nodes,iws,errorFlag)
         endif
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_trac_tge101_2(tcr,rgeo,rload,id,rd,
     &                            erem,no_of_nodes,iws,errorFlag)
c
c  Defines Traction Boundary Condition on geometry TGE = 101
c          TBC = 2: Prescribed endtractions  
c
c  Define the subroutines "ftrac101_2" and "fjac1" as external so that
c  they can be passed as arguments when "wrp_cnf_driver" is called.
c
      implicit none
      external ftrac101_2,fjac1
      integer  tcr,id(3,*),erem(2,*),no_of_nodes,iws,errorFlag,idof(3),i
      double precision rgeo(*),rload(*),rd(3,*)
      double precision w0,t0,fz,mx,my,area,ixx,iyy,ct(6),cj(6)
c
      do i=1,6
         ct(i) = 0.
         cj(i) = 0.
      enddo
c
      w0 = rgeo(1)
      t0 = rgeo(3)
c
      fz = rload(3)          
      mx = rload(4)          
      my = rload(5)
c
c  Define non-zero constants to be used to calculate the traction
c  at each Gauss point
c
      if ( tcr .eq. 1 ) then
         area = 4.*w0*t0
         ct(1) = fz / area
c
      elseif (tcr.eq.2) then
         area = 2.*w0*t0
         ixx  = 2.*w0*t0**3. / 12.
         ct(1)= - ( -fz/area + (mx/ixx)*(t0/2.) )
         ct(3)=   mx/ixx
c
      elseif (tcr.eq.3) then
         area = 2.*w0*t0
         iyy  = 2.*t0*w0**3. / 12.
         ct(1)=  ( fz/area + (my/iyy)*(w0/2.) )
         ct(2)= - my/iyy
c
      elseif (tcr.eq.4) then
         area = w0*t0
         ixx  = w0*t0**3. / 12.
         iyy  = t0*w0**3. / 12.
         ct(1)= - ( - fz/area + (mx/ixx)*(t0/2.) - (my/iyy)*(w0/2.) )
         ct(3)=   mx/ixx
         ct(2)= - my/iyy
c
      endif
c
c  Define constants used to calculate the partial derivatives of the
c  LOCAL surface coordinates with respect to the GLOBAL ccordinates
c  at each Gauss point.
c
      cj(1) = -1.0
      cj(5) =  1.0
c
c  Compute the consistent nodal force vector;
c
      idof(1) = 0
      idof(2) = 0
      idof(3) = 2
      call wrp_cnf_driver(ftrac101_2,fjac1,ct,cj,
     &     erem,idof,id,rd,no_of_nodes,iws,errorFlag)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine ftrac101_2(x,ct,traction)
      implicit none
      double precision x(*),ct(*),traction(3)
c
c Computes traction according to tx = ty = 0 tz = C0 + C1x*X + C1y*Y
c
      traction(1) = 0.
      traction(2) = 0.
      traction(3) = ct(1) + ct(2)*x(1) + ct(3)*x(2)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine fjac1(x,cj,dudx,dvdx)
      implicit none
      double precision x(*),cj(*),dudx(3),dvdx(3)
      dudx(1) = cj(1)
      dudx(2) = cj(2)
      dudx(3) = cj(3)
      dvdx(1) = cj(4)
      dvdx(2) = cj(5)
      dvdx(3) = cj(6)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_trac_tge200_2(tge,rgeo,rload,ip,id,rd,
     &           erem,esid,no_of_nodes,iws,errorFlag)
c
c  Defines Traction Boundary Condition on geometry TGE = 201 - 204
c          TBC = 2: Prescribed endtractions  
c
c  Define the subroutines "ftrac200_2" and "fjac1" as external so that
c  they can be passed as arguments when "wrp_cnf_driver" is called.
c
      implicit none
      external ftrac200_2,fjac1
      integer  tge,id(3,*),erem(2,*),esid(2,*),no_of_nodes,
     &         iws,errorFlag,idof(3),i
      double precision rgeo(*),rload(*),rd(3,*)
      double precision t,r0,nx,mz,area,izz,pi,p0,ct(6),cj(6)
      logical          ip
c
      do i=1,6
         ct(i) = 0.
         cj(i) = 0.
      enddo
c
      t = rgeo(3)
      r0 = rgeo(4)
      pi = 4.*atan(1.)
c
      nx = rload(1)          
      mz = rload(6)          
      if (ip) then
         p0 = rload(7)
      else
         p0 = 0.0
      endif          
      area = 2.*pi*( r0 + 0.5*t )*t
      izz  = 0.25*pi*( (r0+t)**4. - r0**4. )
      ct(1) = nx/area + p0*r0 / ( (2.+t/r0)*t )
      ct(2) = -mz/izz
c
c  Define constants used to calculate the partial derivatives of the
c  LOCAL surface coordinates with respect to the GLOBAL ccordinates
c  at each Gauss point.
c
      cj(3) = -1.0
      cj(5) =  1.0
c
c  Compute the consistent nodal force vector;
c
      idof(1) = 2
      idof(2) = 0
      idof(3) = 0
      if ( (tge.eq.201).or.(tge.eq.202) ) then
         call wrp_cnf_driver(ftrac200_2,fjac1,ct,cj,
     &        esid,idof,id,rd,no_of_nodes,iws,errorFlag)
      elseif ( (tge.eq.203).or.(tge.eq.204) ) then
         call wrp_cnf_driver(ftrac200_2,fjac1,ct,cj,
     &        erem,idof,id,rd,no_of_nodes,iws,errorFlag)
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine ftrac200_2(x,ct,traction)
      implicit none
      double precision x(*),ct(*),traction(3)
c
c Computes traction according to tx = C0 + C1y*Y  ty = tz = 0
c
      traction(1) = ct(1) + ct(2)*x(2)
      traction(2) = 0.
      traction(3) = 0.
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_trac_tge300_2(tge,rgeo,rload,ip,id,rd,
     &           erem,esid,no_of_nodes,iws,errorFlag)
c
c  Defines Traction Boundary Condition on geometry TGE = 301 - 308
c          TBC = 2: Prescribed endtractions  
c
c  Define the subroutines "ftrac300_2" and "fjac1" as external so that
c  they can be passed as arguments when "wrp_cnf_driver" is called.
c
      implicit none
      external ftrac300_2,fjac1
      integer  tge,id(3,*),erem(2,*),esid(2,*),no_of_nodes,
     &         iws,errorFlag,idof(3),i
      double precision rgeo(*),rload(*),rd(3,*)
      double precision t,rc,rp,n,mz,p0,area,jzz,pi,phi0,ct(10),cj(10)
      logical          ip
c
      do i=1,10
         ct(i) = 0.
         cj(i) = 0.
      enddo
c
      pi = 4.*atan(1.)
      t    = rgeo(3)
      rc   = rgeo(4)
      rp   = rgeo(5)
      phi0 = rgeo(6)
c
      if (ip) then
         p0 = rload(7)
      else
         p0 = 0.0
      endif          
      n    = p0*pi*rc*rc          
      mz   = rload(6)          
      area = 2.*pi*( rc + 0.5*t )*t
      jzz  = 2.*pi*rp*rp*( sqrt(1.-(rc/rp)**2.) -
     &        sqrt(1.-((rc+t)/rp)**2.) - (rc+0.5*t)*t/rp/rp )
      ct(1) = (n + mz/rp ) / area
      ct(2) = (mz/rp) / jzz     
      ct(3) = sin(phi0) / rp
      ct(4) = cos(phi0) / rp
      ct(5) = -(1.-cos(phi0))
      ct(6) =  cos(phi0)
      ct(7) = -sin(phi0)
c
c  Define constants used to calculate the partial derivatives of the
c  LOCAL surface coordinates with respect to the GLOBAL ccordinates
c  at each Gauss point.
c
      cj(3) = -1.0
      cj(4) = sin(phi0) 
      cj(5) = cos(phi0) 
c
c  Compute the consistent nodal force vector;
c
      idof(1) = 2
      idof(2) = 2
      idof(3) = 0
      if ( (tge.eq.301).or.(tge.eq.302).or.(tge.eq.305).or.
     &     (tge.eq.306) ) then
         call wrp_cnf_driver(ftrac300_2,fjac1,ct,cj,
     &        esid,idof,id,rd,no_of_nodes,iws,errorFlag)
      elseif ( (tge.eq.303).or.(tge.eq.304).or.(tge.eq.307).or.
     &         (tge.eq.308) ) then
         call wrp_cnf_driver(ftrac300_2,fjac1,ct,cj,
     &        erem,idof,id,rd,no_of_nodes,iws,errorFlag)
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine ftrac300_2(x,ct,traction)
      implicit none
      double precision x(*),ct(*),traction(3),snn,r
c
c Computes traction according to
c    tx = Snn*C6  ty = Snn*C7   tz = 0  where
c    Snn = c1 + r/(1+r)*c2   r = X*C3 + Y*C4 + C5
c
      r = x(1)*ct(3) + x(2)*ct(4) + ct(5)
      snn = ct(1) + r/(1.-r)*ct(2)
      traction(1) = snn*ct(6)
      traction(2) = snn*ct(7)
      traction(3) = 0.
      return
      end

c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_trac_tge101_3(rload,id,rd,ecrsur,ecrtip,
     &                             no_of_nodes,iws,errorFlag)
c														  
c  Defines Traction Boundary Condition on geometry TGE = 101
c          TBC = 3: Prescribed tractions on crack face
c
c  Define the subroutines "ftrac101_3" and "fjac1" as external so
c  that they can be passed as arguments when "wrp_cnf_driver" is called.
c
      implicit none
      external ftrac101_3,fjac1
      integer  id(3,*),ecrsur(2,*),ecrtip(2,*),no_of_nodes,
     &         iws,errorFlag,i,idof(3)
      double precision rload(*),rd(3,*),ct(10),cj(10)
c
      double precision t,w,c,a,alfa
      common /geom/    t,w,c,a,alfa
c
      do i=1,10
         ct(i) = 0.
         cj(i) = 0.
      enddo
c
      ct(1) = rload(1)
      ct(2) = rload(2)
      ct(3) = rload(3)
      ct(4) = rload(4)
      ct(5) = rload(5)
      ct(6) = rload(6)
      ct(7) = rload(7)
      ct(8) = 1./a
      ct(9) = 1./c
c
c  Define constants used to calculate the partial derivatives of the
c  LOCAL surface coordinates with respect to the GLOBAL ccordinates
c  at each Gauss point.
c
      cj(1) = 1.0
      cj(5) = 1.0
c
c  Compute the consistent nodal force vector;
c
      idof(1) = 0
      idof(2) = 0
      idof(3) = 3
      call wrp_cnf_driver(ftrac101_3,fjac1,ct,cj,
     &     ecrsur,idof,id,rd,no_of_nodes,iws,errorFlag)
c      call wrp_cnf_driver(ftrac101_3,fjac1,ct,cj,
c     &     ecrtip,idof,id,rd,no_of_nodes,iws,errorFlag)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine ftrac101_3(x,ct,traction)
      implicit none
      double precision x(*),ct(*),traction(3),u,v
c
c Computes traction according to tx = ty = 0  and
c    tz = -( c0 + u*(c1 + u*(c2 + u*c3)) + v*(c1 + v*(c2 + v*c3))) 
c    u = Y*c8   v = X*c9
c
      u = x(2)*ct(8)
      v = x(1)*ct(9)
      traction(1) = 0.
      traction(2) = 0.
      traction(3) = - ( ct(1) + u*(ct(2) + u*(ct(3) + u*ct(4) ) )
     &                        + v*(ct(5) + v*(ct(6) + v*ct(7) ) ) )
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_trac_tge200_3(tge,rgeo,rload,id,rd,ecrsur,
     &           no_of_nodes,iws,errorFlag)
c
c  Defines Traction Boundary Condition on geometry TGE = 201-204
c          TBC = 3: Prescribed tractions on crack face
c
c  Define the subroutines "ftrac201_3", "ftrac203_3" and "fjac1"
c  as external so that they can be passed as arguments when
c  "wrp_cnf_driver" is called.
c
      implicit none
      external ftrac201_3,ftrac203_3,fjac1
      integer  tge,id(3,*),ecrsur(2,*),no_of_nodes,iws,errorFlag,
     &         i,idof(3)
      double precision rgeo(*),rload(*),rd(3,*),ct(10),cj(10),rc
c
      double precision t,w,c,a,alfa
      common /geom/    t,w,c,a,alfa
c
      do i=1, 10
         ct(i) = 0.
         cj(i) = 0.
      enddo
      idof(1) = 0
      idof(2) = 0
      idof(3) = 0
      rc = rgeo(4)
c
      ct(1) = rload(1)
      ct(2) = rload(2)
      ct(3) = rload(3)
      ct(4) = rload(4)
      ct(5) = rload(5)
      ct(6) = rload(6)
      ct(7) = rload(7)
c
      if (tge.eq.201) then
         ct(8) = -rc/a
         ct(9) =  1./a
         ct(10)=  1./c
         cj(1) = 1.0
         cj(5) = 1.0
         idof(3) = 3
         call wrp_cnf_driver(ftrac201_3,fjac1,ct,cj,
     &        ecrsur,idof,id,rd,no_of_nodes,iws,errorFlag)
      elseif (tge.eq.202) then
         ct(8) = (rc+t)/a
         ct(9) = -1./a
         ct(10)=  1./c
         cj(1) = 1.0
         cj(5) = 1.0
         idof(3) = 3
         call wrp_cnf_driver(ftrac201_3,fjac1,ct,cj,
     &        ecrsur,idof,id,rd,no_of_nodes,iws,errorFlag)
      elseif (tge.eq.203) then
         ct(8) = -rc/a
         ct(9) =  1./a
         ct(10)=  rc/c
         cj(3) = 1.0
         cj(5) = 1.0
         idof(3) = 3
         call wrp_cnf_driver(ftrac203_3,fjac1,ct,cj,
     &        ecrsur,idof,id,rd,no_of_nodes,iws,errorFlag)
      elseif (tge.eq.204) then
         ct(8) = (rc+t)/a
         ct(9) = -1./a
         ct(10)= (rc+t)/c
         cj(3) = 1.0
         cj(5) = 1.0
         idof(3) = 3
         call wrp_cnf_driver(ftrac203_3,fjac1,ct,cj,
     &        ecrsur,idof,id,rd,no_of_nodes,iws,errorFlag)
      endif
c
c  Note that cj(i) is used to define constants to be used in the
c  calculation of the partial derivatives of the LOCAL surface
c  coordinates with respect to the GLOBAL ccordinates at each
c  Gauss point.
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine ftrac201_3(x,ct,traction)
      implicit none
      double precision x(*),ct(*),traction(3),u,v
c
c Computes traction according to tx = ty = 0  and
c        tz = c0 + u*(c1 + u*(c2 + u*c3)) + v*(c1 + v*(c2 + v*c3)) 
c        u = c8 + Y*c9   v = X*c10
c
      u = ct(8) + x(2)*ct(9)
      v = x(1)*ct(10)
      traction(1) = 0.
      traction(2) = 0.
      traction(3) = - ( ct(1) + u*(ct(2) + u*(ct(3) + u*ct(4) ) )
     &                        + v*(ct(5) + v*(ct(6) + v*ct(7) ) ) )
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine ftrac203_3(x,ct,traction)
      implicit none
      double precision x(*),ct(*),traction(3),r,phi,u,v
c
c Computes traction according to ty = tz = 0  and
c        tx = C0 + u*(c1 + u*(c2 + u*c3)) + v*(c1 + v*(c2 + v*c3)) 
c        u = c8 + Y*c9   v = X*c10
c
      r   = sqrt( x(2)*x(2) + x(3)*x(3) )
      phi = asin(abs(x(3))/r)
      u = ct(8) + r*ct(9)
      v = phi*ct(10)
      traction(1) = ct(1) + u*(ct(2) + u*(ct(3) + u*ct(4) ) )
     &                    + v*(ct(5) + v*(ct(6) + v*ct(7) ) )
      traction(2) = 0.
      traction(3) = 0.
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_trac_tge300_3(tge,rgeo,rload,id,rd,ecrsur,
     &                             no_of_nodes,iws,errorFlag)
c
c  Defines Traction Boundary Condition on geometry TGE = 301-308
c          TBC = 3: Prescribed tractions on crack face
c
c  Define the subroutines "ftrac301_3", "ftrac303_3" and "fjac1"
c  as external so that they can be passed as arguments when
c  "wrp_cnf_driver" is called.
c
      implicit none
      external ftrac301_3,ftrac303_3,fjac1
      integer  tge,id(3,*),ecrsur(2,*),no_of_nodes,iws,errorFlag,
     &         i,idof(3)
      double precision rgeo(*),rload(*),rd(3,*),ct(11),cj(11),rc,rp
c
      double precision t,w,c,a,alfa
      common /geom/    t,w,c,a,alfa
c
      do i=1, 11
         ct(i) = 0.
         cj(i) = 0.
      enddo
      idof(1) = 0
      idof(2) = 0
      idof(3) = 0
      rc = rgeo(4)
      rp = rgeo(5)
c
      ct(1) = rload(1)
      ct(2) = rload(2)
      ct(3) = rload(3)
      ct(4) = rload(4)
      ct(5) = rload(5)
      ct(6) = rload(6)
      ct(7) = rload(7)
c
      if (tge.eq.301) then
         ct(8) = -(rp+rc)/a
         ct(9) =   1./a
         ct(10)=  (rp+rc)/c
         ct(11)=   rp
         cj(1) = 1.0
         cj(5) = 1.0
         idof(3) = 3
         call wrp_cnf_driver(ftrac301_3,fjac1,ct,cj,
     &        ecrsur,idof,id,rd,no_of_nodes,iws,errorFlag)
      elseif (tge.eq.302) then
         ct(8) = (rp+rc+t)/a
         ct(9) =  -1./a
         ct(10)= (rp+rc+t)/c
         ct(11)=  rp
         cj(1) = 1.0
         cj(5) = 1.0
         idof(3) = 3
         call wrp_cnf_driver(ftrac301_3,fjac1,ct,cj,
     &        ecrsur,idof,id,rd,no_of_nodes,iws,errorFlag)
      elseif (tge.eq.305) then
         ct(8) = (rp-rc)/a
         ct(9) =  -1./a
         ct(10)= (rp-rc)/c
         ct(11)=  rp
         cj(1) = 1.0
         cj(5) = 1.0
         idof(3) = 3
         call wrp_cnf_driver(ftrac301_3,fjac1,ct,cj,
     &        ecrsur,idof,id,rd,no_of_nodes,iws,errorFlag)
      elseif (tge.eq.306) then
         ct(8) = -(rp-rc-t)/a
         ct(9) =   1./a
         ct(10)=  (rp-rc-t)/c
         ct(11)=   rp
         cj(1) = 1.0
         cj(5) = 1.0
         idof(3) = 3
         call wrp_cnf_driver(ftrac301_3,fjac1,ct,cj,
     &        ecrsur,idof,id,rd,no_of_nodes,iws,errorFlag)
      elseif (tge.eq.303) then
         ct(8) = -rc/a
         ct(9) =  1./a
         ct(10)=  rc/c
         cj(3) = 1.0
         cj(5) = 1.0
         idof(3) = 3
         call wrp_cnf_driver(ftrac303_3,fjac1,ct,cj,
     &        ecrsur,idof,id,rd,no_of_nodes,iws,errorFlag)
      elseif (tge.eq.304) then
         ct(8) = (rc+t)/a
         ct(9) = -1./a
         ct(10)= (rc+t)/c
         cj(3) = 1.0
         cj(5) = 1.0
         idof(3) = 3
         call wrp_cnf_driver(ftrac303_3,fjac1,ct,cj,
     &        ecrsur,idof,id,rd,no_of_nodes,iws,errorFlag)
      elseif (tge.eq.307) then
         ct(8) = -rc/a
         ct(9) =  1./a
         ct(10)=  rc/c
         cj(3) = 1.0
         cj(5) = 1.0
         idof(3) = 3
         call wrp_cnf_driver(ftrac303_3,fjac1,ct,cj,
     &        ecrsur,idof,id,rd,no_of_nodes,iws,errorFlag)
      elseif (tge.eq.308) then
         ct(8) = (rc+t)/a
         ct(9) = -1./a
         ct(10)= (rc+t)/c
         cj(3) = 1.0
         cj(5) = 1.0
         idof(3) = 3
         call wrp_cnf_driver(ftrac303_3,fjac1,ct,cj,
     &        ecrsur,idof,id,rd,no_of_nodes,iws,errorFlag)
      endif
c
c  Note that cj(i) is used to define constants to be used in the
c  calculation of the partial derivatives of the LOCAL surface
c  coordinates with respect to the GLOBAL ccordinates at each
c  Gauss point.
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine ftrac301_3(x,ct,traction)
      implicit none
      double precision x(*),ct(*),traction(3),u,v,r,phi
c
c Computes traction according to tx = ty = 0  and
c        tz = c0 + u*(c1 + u*(c2 + u*c3)) + v*(c1 + v*(c2 + v*c3)) 
c        u = c8 + r*c9   v = phi*c10
c
      r   = sqrt( x(1)*x(1) + (x(2)+ct(11))*(x(2)+ct(11)) )
      phi = asin( abs( x(1)/r ) ) 
      u = ct(8) + r*ct(9)
      v = phi*ct(10)
      traction(1) = 0.
      traction(2) = 0.
      traction(3) = - ( ct(1) + u*(ct(2) + u*(ct(3) + u*ct(4) ) )
     &                        + v*(ct(5) + v*(ct(6) + v*ct(7) ) ) )
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine ftrac303_3(x,ct,traction)
      implicit none
      double precision x(*),ct(*),traction(3),r,phi,u,v
c
c Computes traction according to ty = tz = 0  and
c        tx = C0 + u*(c1 + u*(c2 + u*c3)) + v*(c1 + v*(c2 + v*c3)) 
c        u = c8 + Y*c9   v = X*c10
c
      r   = sqrt( x(2)*x(2) + x(3)*x(3) )
      phi = asin(abs(x(3))/r)
      u = ct(8) + r*ct(9)
      v = phi*ct(10)
      traction(1) = ct(1) + u*(ct(2) + u*(ct(3) + u*ct(4) ) )
     &                    + v*(ct(5) + v*(ct(6) + v*ct(7) ) )
      traction(2) = 0.
      traction(3) = 0.
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_cnf_driver(ftrac,fjac,ct,cj,eset,idof,id,rd,
     &           no_of_nodes,iws,errorFlag)
c
c  Calculation of the consistent nodal force vector;
c  Repeat the calculation until the maximum norm of the relative
c  difference between two consequtive calculations is less than ETOL.
c
c  Note that the id(k,node) {k=1,2,3} will be set at the first call
c  of SUBR. "wrp_calc_cnf()"
c
c  id(3,*) = type of B.C.
c  rd(3,*)    stores the values of prescribed displacements/tractions
c  rdtmp(3,*) local version of rd(3,*) during the CNF calculation
c
      implicit none
      external  ftrac,fjac
      integer   eset(2,*),idof(*),id(3,*),no_of_nodes,iws,errorFlag
      double precision ct(*),cj(*),rd(3,*)
      integer    ngmax, ns,ng, nmax, n
      parameter (ngmax=20, ns=8*ngmax*ngmax, nmax=50000)
      double precision x1,x2,rnorm1,rnorm2,relerr,etol,rdtmp(3,nmax),
     &                 rg(ngmax),wg(ngmax),wij(ns),
     &                 shp(ns),dsdg(ns),dsdh(ns), area,
     &                 compute_cnf_norm
      parameter (etol=1.e-10)
c
      if (no_of_nodes.gt.nmax) then
         write(iws,'(t1,a,a,i6)') '>> increase nmax in array ',
     &                 ' rdtmp(3,namx);  nmax > ',no_of_nodes
         errorFlag = 1
         return
      endif
c 
c The numerical integration is carried out using a Gauss-Legendre
c scheme. The accuracy of the integration is checked by repeating
c the numerical integration by increasing the number of integration
c points by one each time. The maximum number of Gauss points used
c set set by NGMAX.
c
c Compute the Abscissas and weights in the gauss-Legendre scheme
c used for numerical integration (use subr. from "Numerical Recipes").
c The current number of Gauss points is given by NG
c
      x1 = -1.0
      x2 =  1.0
      ng = 3
      call gauleg(x1,x2,rg,wg,ng)
      call rdtmp_zero(rdtmp,nmax)
      call wrp_calc_cnf(ftrac,fjac,ct,cj,eset,ng,idof,id,rdtmp,
     &     rg,wg,wij,shp,dsdg,dsdh,area)
      call write_gauss_points(rg,wg,ng,area,iws)
      rnorm1 = compute_cnf_norm(rdtmp,no_of_nodes)
c      do ng=40, ngmax
      do ng=4, ngmax
         call gauleg(x1,x2,rg,wg,ng)
         call rdtmp_zero(rdtmp,nmax)
         call wrp_calc_cnf(ftrac,fjac,ct,cj,eset,ng,idof,id,rdtmp,
     &        rg,wg,wij,shp,dsdg,dsdh,area)
         call write_gauss_points(rg,wg,ng,area,iws)
c
         rnorm2 = compute_cnf_norm(rdtmp,no_of_nodes)
         relerr = abs( (rnorm1-rnorm2) / rnorm1 ) 
         write(iws,101) ' CNF: iter =',ng-3,'  RelError = ',relerr
 101     format(t1,a,i3,a,g12.4)
         if (relerr.lt.etol) goto 20
      enddo
c
 20   continue
c
c Add the local force vector, rdtmp(), to the global force vector rd()
c
      do n=1, no_of_nodes
       if((id(1,n).eq.2).or.(id(1,n).eq.3)) rd(1,n)=rd(1,n)+rdtmp(1,n)
       if((id(2,n).eq.2).or.(id(2,n).eq.3)) rd(2,n)=rd(2,n)+rdtmp(2,n)
       if((id(3,n).eq.2).or.(id(3,n).eq.3)) rd(3,n)=rd(3,n)+rdtmp(3,n)
      enddo
c
      write(iws,102) 'Consistent Node Force, iter =',ng,
     &               ' ; MaxIter =',ngmax
 102  format(t1,a,i3,a,i3)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine write_gauss_points(rg,wg,ng,area,iws)
      implicit none
      integer          iws,ng,i
      double precision rg(*),wg(*),area
      write(iws,'(t1,a,i4)') '>> # of Gauss points = ',ng
      do i=1, ng
         write(iws,'(t1,2(a,f16.12))') 'rg=',rg(i),'  wg=',wg(i)
      enddo
      write(iws,'(t1,a,g12.4)') ' Total Surface Area =',area
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function compute_cnf_norm(rd,no_of_nodes)
      implicit none
      integer          no_of_nodes,i
      double precision rd(3,*),rnorm
c
      rnorm = 0.
      do i=1, no_of_nodes
         rnorm = rnorm + rd(1,i)**2. + rd(2,i)**2. + rd(3,i)**2.
      enddo
      compute_cnf_norm = sqrt(rnorm)     
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine rdtmp_zero(rdtmp,nmax)
      implicit none
      integer          nmax,n
      double precision rdtmp(3,*)
c
      do n=1, nmax
         rdtmp(1,n) = 0.
         rdtmp(2,n) = 0.
         rdtmp(3,n) = 0.
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine  wrp_calc_cnf(ftrac,fjac,ct,cj,eset,ng,idof,id,rd,
     &                         rg,wg,wij,shp,dsdg,dsdh,area)
c
c The routine calculates the consistant nodal forces for a given
c stress distribution on an element face. The local node numbering
c for an element face follow the definition in WARP3D.
c
c The Current version ONLY works for an 8-NODED BRICK ELEMENT
c
      implicit none
      external ftrac,fjac
      include 'plate_common_nod.f'
      include 'plate_common_eln.f'
c
      integer          eset(2,*),ng,idof(3),id(3,*)
      double precision ct(*),cj(*),rd(3,*)
      integer    nn
      parameter (nn=4)
c
      integer    i,j,k,nloc,ie,ne,face_no,fnode(8,6),np(8)
      double precision xn(3,8),cnf(3,8),g,h,rg(*),wg(*),wij(ng,*),
     &                 shp(4,ng,*),shp4(4),dsdg(4,ng,*),
     &                 dsdh(4,ng,*),  dsdg4(4),dsdh4(4)
      double precision area,elarea
c     
      data  (fnode(k,1),k=1,8)/ 1, 4, 3, 2,  12,11,10, 9 /
      data  (fnode(k,2),k=1,8)/ 5, 6, 7, 8,  13,14,15,16 /
      data  (fnode(k,3),k=1,8)/ 1, 2, 6, 5,   9,18,13,17 /
      data  (fnode(k,4),k=1,8)/ 3, 4, 8, 7,  11,20,15,19 /
      data  (fnode(k,5),k=1,8)/ 2, 3, 7, 6,  10,19,14,18 /
      data  (fnode(k,6),k=1,8)/ 1, 5, 8, 4,  17,16,20,12 /
c
c  Calculate values of the shape functions and it's partial derivatives
c  at all Gauss points.
c
CTEST      write(10,'(t1,a,8f8.5)') ' g: ',(rg(i),i=1,ng)     
c      do j=1,ng
c         do i=1,ng
c            g = rg(i)
c            h = rg(j)
c           call shape2d4(g,h,shp4)
c           write(10,'(t1,4g12.4)') (shp4(k),k=1,nn)
c         enddo
c      enddo
c
      area = 0.
      do j=1, ng
         do i=1, ng
            wij(i,j) = wg(i)*wg(j)
            g = rg(i)
            h = rg(j)
            call shape2d4(g,h,shp4)
            call dshape2d4(g,h,dsdg4,dsdh4)
            do k=1, nn
               shp(k,i,j)  = shp4(k)
               dsdg(k,i,j) = dsdg4(k)
               dsdh(k,i,j) = dsdh4(k)
            enddo
         enddo
      enddo
CTEST
c      do j=1, ng
c         write(10,'(t1,10f12.6)') (wij(i,j),i=1,ng)
c      enddo
c      do j=1, ng
c        do i=1, ng
c          write(10,'(t1,4f10.4)') (shp(k,i,j),k=1,nn)
c        enddo
c      enddo
c
c  Calculate the contribution to the consistent nodal forces from
c  all element faces located on the surface. 
c
      ne = eset(1,1)
      do i=1, ne
         ie = eset(1,i+1)
         face_no = eset(2,i+1)
c
c      Retrieve local node numbers and coordinats in the element ie
c
         do k=1, nn
            nloc = fnode(k,face_no)
            xn(1,k) = npos( eln(ie,nloc), 1)
            xn(2,k) = npos( eln(ie,nloc), 2)
            xn(3,k) = npos( eln(ie,nloc), 3)
            np(k)   = nnr(  eln(ie,nloc) )
         enddo 
c
c      Calculate the consistent node forces 
c
CTEST
c         do j=1, nn
c            write(10,*) 'node=',np(j),' XYZ:',(real(xn(k,j)),k=1,3)
c         enddo
         call  wrp_elem_cnf(ftrac,fjac,ct,cj,xn,nn,wij,ng,
     &                      shp,dsdg,dsdh,cnf,elarea)
c
c      Add the the consistent node forces to the GLOBAL array
c
         do k=1, nn
            if (idof(1).gt.0) then
               rd(1,np(k)) = rd(1,np(k)) + cnf(1,k)
               id(1,np(k)) = idof(1)
            endif
            if (idof(2).gt.0) then
               rd(2,np(k)) = rd(2,np(k)) + cnf(2,k)
               id(2,np(k)) = idof(2)
            endif
            if (idof(3).gt.0) then
               rd(3,np(k)) = rd(3,np(k)) + cnf(3,k)
               id(3,np(k)) = idof(3)
            endif
         enddo
         area = area + elarea
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine  wrp_elem_cnf(ftrac,fjac,ct,cj,xn,nn,wg,ng,
     &                         shp,dsdg,dsdh,cnf,elarea)
      implicit none
c
c  Calculates the consistent nodal forces associated with one
c  element. The local element face is defined by the the local
c  node coordinates as given in xn(i)
c
c  xn(3,i) node coordinates;   nn # of node coordinates in the element
c  wg(ng,ng)     weights in the Gauss quadrature scheme
c  shp(i,ng,ng)  shape function i at all Gauss points  
c  dsdg(i,ng,ng) partial derivatives of shape function i
c  dsdh(i,ng,ng) partial derivatives of shape function i
c  ftrac(x,ct, ...) external function for traction, t(i), calculations
c          uses constants in array ct(i) 
c  fjac(x,cj, ...)  external function for calculating the jacobian
c          uses constants in array cj(i)
c  dudx() and dvdx() are partial derivatives of the local surface
c      coordinates w.r.t. the global coordinates X,Y,Z - neded
c      in the calculation of the determinant of the Jacobian 
c  
c  cnf(3,n) the computed consistent nodal force components (X,Y,Z)
c           at all the node points 
c
      external ftrac,fjac
      integer          ng,nn
      double precision cnf(3,*),xn(3,*),ct(*),cj(*),wg(ng,*),
     &                 shp(nn,ng,*),dsdg(nn,ng,*),dsdh(nn,ng,*)
      integer          i,j,k,n
      double precision x(3),t(3),jac(2,2),jdet,
     &                 dxdg(3),dxdh(3),dudx(3),dvdx(3),elarea
c
      do n=1, nn
         cnf(1,n) = 0.
         cnf(2,n) = 0.
         cnf(3,n) = 0.
      enddo
      elarea = 0.
CTEST
c      do j=1, ng
c         write(10,'(t1,10f12.6)') (wg(i,j),i=1,ng)
c      enddo
c      do j=1, ng
c        do i=1, ng
c          write(10,'(t1,4f10.4)') (shp(n,i,j),n=1,nn)
c        enddo
c      enddo
c      do n=1, nn
c         write(10,*) ' X:',(real(xn(k,n)),k=1,3)
c      enddo
c      if (ng.gt.2) then
c         close(10)
c         stop
c      endif
c
c Integration by Gauss-Legendre quadrature using a ng x ng points.
c The double summation is caried out over i and j.
c
      do j=1, ng
         do i=1, ng
c
c  Set to zero
c
            jac(1,1) = 0.
            jac(1,2) = 0.
            jac(2,1) = 0.
            jac(2,2) = 0.
            do k=1, 3
               x(k) = 0.
               t(k) = 0.
               dxdg(k) = 0.
               dxdh(k) = 0.
            enddo
c
c  Calculate the global coordinates (X,Y,Z) at the current Gauss point
c
            do n=1, nn
               x(1) = x(1) +  xn(1,n)*shp(n,i,j)
               x(2) = x(2) +  xn(2,n)*shp(n,i,j)
               x(3) = x(3) +  xn(3,n)*shp(n,i,j)
               dxdg(1) = dxdg(1) + xn(1,n)*dsdg(n,i,j)
               dxdg(2) = dxdg(2) + xn(2,n)*dsdg(n,i,j)
               dxdg(3) = dxdg(3) + xn(3,n)*dsdg(n,i,j)
               dxdh(1) = dxdh(1) + xn(1,n)*dsdh(n,i,j)
               dxdh(2) = dxdh(2) + xn(2,n)*dsdh(n,i,j)
               dxdh(3) = dxdh(3) + xn(3,n)*dsdh(n,i,j)
            enddo
cTEST            write(10,*) ' XYZg: ',(real(x(k)),k=1,3)
c
c  Calculate the traction, du/dx and dv/dx at the Gauss point using
c  a geometry specific routines.
c     
            call ftrac(x,ct,t)
            call fjac(x,cj,dudx,dvdx)
c            write(10,*) ' Du/DX :',(real(dudx(k)),k=1,3)
c            write(10,*) ' Dv/DX :',(real(dvdx(k)),k=1,3)
c
c  Transformation from local surface coordinates (u,v) to global
c  coordinates (X,Y,Z).
c
            do k=1, 3
               jac(1,1) = jac(1,1) + dudx(k)*dxdg(k)
               jac(1,2) = jac(1,2) + dudx(k)*dxdh(k)
               jac(2,1) = jac(2,1) + dvdx(k)*dxdg(k)
               jac(2,2) = jac(2,2) + dvdx(k)*dxdh(k)
            enddo
            jdet = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
c
c  Integration by Gauss quadrature
c
            do n=1, nn
               cnf(1,n) = cnf(1,n) + wg(i,j)*shp(n,i,j)*t(1)*jdet
               cnf(2,n) = cnf(2,n) + wg(i,j)*shp(n,i,j)*t(2)*jdet
               cnf(3,n) = cnf(3,n) + wg(i,j)*shp(n,i,j)*t(3)*jdet
            enddo
c
c  Check if element area is correct!
c
            elarea = elarea + wg(i,j)*jdet
c
         enddo
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine shape2d4(g,h,shp4)
c
c  This routine computes the shape functions for a 4 node 2D element at
c  the point (g,h). The shape functions can for instance be found in
c  the ABAQUS theory manual Ver. 5.7 pp. 3.2.3-3 - 3.2.3-4
c
      double precision g,h, shp4(8), one,two,four
      parameter (one=1.0, two=2.0, four=4.0 )
c
      shp4(1) = (one-g)*(one-h) / four
      shp4(2) = (one+g)*(one-h) / four
      shp4(3) = (one+g)*(one+h) / four
      shp4(4) = (one-g)*(one+h) / four
c     
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine dshape2d4(g,h,dsdg4,dsdh4)
c
c  This routine computes the derivatives of the shape functions for a
c  4 node 2D element at the point (g,h). The shape functions can for
c  instance be found in the ABAQUS theory manual Ver. 5.7
c  pp. 3.2.3-3 - 3.2.3-4
c
      double precision g,h, dsdg4(4),dsdh4(4), one,two,four
      parameter (one=1.0, two=2.0, four=4.0 )
c
      dsdg4(1) = -(one-h) / four
      dsdg4(2) =  (one-h) / four
      dsdg4(3) =  (one+h) / four
      dsdg4(4) = -(one+h) / four
c
      dsdh4(1) = -(one-g) / four
      dsdh4(2) = -(one+g) / four
      dsdh4(3) =  (one+g) / four
      dsdh4(4) =  (one-g) / four
c     
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine shape2d8(g,h,shp8)
c
c  This routine computes the shape functions for a 8 node 2D element
c  at the point (g,h).  The shape functions can for instance be found
c  in the ABAQUS theory manual Ver. 5.7 pp. 3.2.3-3 - 3.2.3-4
c
      double precision g,h, shp8(8), one,two,four
      parameter (one=1.0, two=2.0, four=4.0 )
c
      shp8(1) = -(one-g)*(one-h)*(one+g+h)/four
      shp8(2) = -(one+g)*(one-h)*(one-g+h)/four
      shp8(3) = -(one+g)*(one+h)*(one-g-h)/four
      shp8(4) = -(one-g)*(one+h)*(one+g-h)/four
      shp8(5) =  (one-g)*(one+g)*(one-h)/two
      shp8(6) =  (one-h)*(one+h)*(one+g)/two
      shp8(7) =  (one-g)*(one+g)*(one+h)/two
      shp8(8) =  (one-h)*(one+h)*(one-g)/two
c     
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine dshape2d8(g,h,dsdg8,dsdh8)
c
c  This routine computes the derivatives of the shape functions for a 8
c  node 2D element at the point (g,h). The shape functions can for
c  instance be found in the ABAQUS theory manual Ver. 5.7
c  pp. 3.2.3-3 - 3.2.3-4
c
      double precision g,h, dsdg8(8),dsdh8(8), one,two,four
      parameter (one=1.0, two=2.0, four=4.0 )
c
      dsdg8(1) = ( (one-h)*(one+g+h) - (one-g)*(one-h) ) / four
      dsdg8(2) = (-(one-h)*(one-g+h) + (one+g)*(one-h) ) / four
      dsdg8(3) = (-(one+h)*(one-g-h) + (one+g)*(one+h) ) / four
      dsdg8(4) = ( (one+h)*(one+g-h) - (one-g)*(one+h) ) / four
      dsdg8(5) = (-(one+g)*(one-h)   + (one-g)*(one-h) ) / two
      dsdg8(6) =   (one-h)*(one+h) / two
      dsdg8(7) = (-(one+g)*(one+h)   + (one-g)*(one+h) ) / two
      dsdg8(8) =  -(one-h)*(one+h) / two
c
      dsdh8(1) = ( (one-g)*(one+g+h) - (one-g)*(one-h) ) /four	
      dsdh8(2) = ( (one+g)*(one-g+h) - (one+g)*(one-h) ) /four
      dsdh8(3) = (-(one+g)*(one-g-h) + (one+g)*(one+h) ) /four
      dsdh8(4) = (-(one-g)*(one+g-h) + (one-g)*(one+h) ) /four
      dsdh8(5) =  -(one-g)*(one+g) / two
      dsdh8(6) = (-(one+g)*(one+h)   + (one+g)*(one-h) ) / two
      dsdh8(7) =   (one-g)*(one+g) / two
      dsdh8(8) = (-(one-g)*(one+h)   + (one-g)*(one-h) ) / two
c     
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine gauleg(x1,x2,x,w,n)
      implicit none
      integer n
      double precision x1,x2,x(n),w(n)
c
c This routine is taken rom "Numerical Recipes in FORTRAN" 2:nd Ed.
c (1992) pp. 140 - 149. 
c
c Given the lower and upper limits of integration x1 and x2, and n,
c this routine returns arrays x(1:n) and w(1:n) of length n,
c containing the abscissas and weights of the Gauss-Legendre n-point
c quadrature formula. High precision is a good idea for this routine.
c EPS is the relative precision.
c 
c Note that the roots are symmetric in the interval, so that only
c "half" of them has to be determined.
c 
c
      integer i,j,m
      double precision p1,p2,p3,pp,xl,xm,z,z1,eps
      parameter (eps=3.d-14)
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do i=1, m
         z=cos( 3.141592654d0 * (i-.25d0) / (n+0.5d0) )
10       continue
            p1=1.d0
            p2=0.d0
            do j=1, n
               p3=p2
               p2=p1
               p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
            enddo
            pp=n*(z*p1-p2)/(z*z-1.d0)
            z1=z
            z=z1-p1/pp
         if (abs(z-z1).gt.eps)goto 10
         x(i)     = xm - xl*z
         x(n+1-i) = xm + xl*z
         w(i)     = 2.d0*xl / ((1.d0-z*z)*pp*pp)
         w(n+1-i) = w(i)
      enddo
      return
      END
c
c----67--1---------2---------3---------4---------5---------6---------712
c
