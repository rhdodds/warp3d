c----67--1---------2---------3---------4---------5---------6---------712
c
c File: plate_crdtransf.f   latest revision 24 April 1998
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine coord_transf(tge,elnum,etyp,rgeo,crfnvec,cqvec,ncq,
     &                        ns1,ns2,ns3,nods)
c
c  Coordinate transformation: the basic plate is transformed
c  into the desired geometry as given by TGE
c
      implicit none
c
      integer  tge,elnum,etyp, ncq,
     &         ns1,ns2,ns3,nods(ns1,ns2,ns3)
      double precision rgeo(*),crfnvec(3),cqvec(3,*)
c
      if     (tge .eq. 101) then
         call  geo_pla_101(crfnvec,cqvec,ncq,ns1,ns2,ns3,nods)
c
      elseif (tge .eq. 201) then
         call  geo_cyl_201(rgeo,crfnvec,cqvec,ncq,ns1,ns2,ns3,nods)
c
      elseif (tge .eq. 202) then
         call  geo_cyl_202(rgeo,crfnvec,cqvec,ncq,elnum,etyp,
     &                     ns1,ns2,ns3,nods)
c
      elseif (tge .eq. 203) then
         call  geo_cyl_203(rgeo,crfnvec,cqvec,ncq,elnum,etyp,
     &                     ns1,ns2,ns3,nods)
c
      elseif (tge .eq. 204) then
         call  geo_cyl_204(rgeo,crfnvec,cqvec,ncq,ns1,ns2,ns3,nods)
c
      elseif (tge .eq. 301) then
         call  geo_elb_301(rgeo,crfnvec,cqvec,ncq,ns1,ns2,ns3,nods)
c
      elseif (tge .eq. 302) then
         call  geo_elb_302(rgeo,crfnvec,cqvec,ncq,elnum,etyp,
     &                     ns1,ns2,ns3,nods)
c
      elseif (tge .eq. 303) then
         call  geo_elb_303(rgeo,crfnvec,cqvec,ncq,elnum,etyp,
     &                     ns1,ns2,ns3,nods)
c
      elseif (tge .eq. 304) then
         call  geo_elb_304(rgeo,crfnvec,cqvec,ncq,ns1,ns2,ns3,nods)
c
      elseif (tge .eq. 305) then
         call  geo_elb_305(rgeo,crfnvec,cqvec,ncq,elnum,etyp,
     &                     ns1,ns2,ns3,nods)
c
      elseif (tge .eq. 306) then
         call  geo_elb_306(rgeo,crfnvec,cqvec,ncq,ns1,ns2,ns3,nods)
c
      elseif (tge .eq. 307) then
         call  geo_elb_307(rgeo,crfnvec,cqvec,ncq,ns1,ns2,ns3,nods)
c
      elseif (tge .eq. 308) then
         call  geo_elb_308(rgeo,crfnvec,cqvec,ncq,elnum,etyp,
     &                     ns1,ns2,ns3,nods)
c
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine geo_pla_101(crfnvec,cqvec,ncq, ns1,ns2,ns3,nods)
c
c Basic Plate: change the reference coordinate system
c
      implicit none
      include 'plate_common_nod.f'
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3), ncq
      integer  nset(200),i
      double precision crfnvec(3),cqvec(3,*)
      double precision x0,y0,z0, cset(3,200),x,y,qx,qy,qz,q
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      do i=1, inm
         if (nnr(i).gt.0) then
            x0 = npos(i,1)
            y0 = npos(i,2)
            z0 = npos(i,3)
c
            npos(i,1) =  x0
            npos(i,2) =  z0
            npos(i,3) = -y0
         endif
      enddo
c
c Define normal vector of crack face
c
      crfnvec(1) =  0.0 
      crfnvec(2) =  0.0
      crfnvec(3) =  1.0
c
c Calculate the normal vector of crack front
c
      call ncrd_ntip(nset,cset,ncq, nods,ns1,ns2,ns3)
cc      call calc_crack_front_normal(cqvec,cset,ncq,1,1,0)
c      write(44,'(i10)') ncq
      do i=1, ncq
         x = cset(1,i)
         y = cset(2,i)
         qx = 2.*x/(c*c)
         qy = 2.*y/(a*a)
         qz = 0.
         q  = sqrt(qx*qx + qy*qy + qz*qz)
         cqvec(1,i) = qx / q
         cqvec(2,i) = qy / q
         cqvec(3,i) = qz / q
c
c         q = acosd(x/c)
c         write(44,'(i5,4g15.6)') i, x,y,0., q
c
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine  geo_cyl_201(rgeo,crfnvec,cqvec,ncq,ns1,ns2,ns3,nods)
c
c Cylinder: internal axial crack 
c
      implicit none
      include 'plate_common_nod.f'
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3), ncq
      integer  nset(200),i
      double precision rgeo(*),crfnvec(3),cqvec(3,*)
      double precision rc,r0,r,theta,x0,y0,z0,
     &                 cset(3,200),x,y,qx,qy,qz,q
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      rc = rgeo(4)
      r0 = rc
c
      do i=1, inm
         if (nnr(i).gt.0) then
            x0 = npos(i,1)
            y0 = npos(i,2)
            z0 = npos(i,3)
c
            r = rc + z0
            theta = y0 / r0
            npos(i,1) =  x0
            npos(i,2) =  r*cos(theta)
            npos(i,3) = -r*sin(theta)
         endif
      enddo
c
c Define normal vector of crack face
c
      crfnvec(1) =  0.0
      crfnvec(2) =  0.0
      crfnvec(3) =  1.0
c
c Calculate the normal vector of crack front
c
      call ncrd_ntip(nset,cset,ncq, nods,ns1,ns2,ns3)
cc      call calc_crack_front_normal(cqvec,cset,ncq,1,1,0)
      do i=1, ncq
         x = cset(1,i)
         y = cset(2,i)
         qx = 2.*x/(c*c)
         qy = 2.*(y-rc)/(a*a)
         qz = 0.
         q  = sqrt(qx*qx + qy*qy + qz*qz)
         cqvec(1,i) = qx / q
         cqvec(2,i) = qy / q
         cqvec(3,i) = qz / q
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine geo_cyl_202(rgeo,crfnvec,cqvec,ncq,elnum,etyp,
     &                       ns1,ns2,ns3,nods)
c
c Cylinder: external axial crack
c
c - Reorder local node numbering order in elements
c
      implicit none
      include 'plate_common_nod.f'
      integer  elnum,etyp, ns1,ns2,ns3,nods(ns1,ns2,ns3), ncq
      integer  nset(200),i
      double precision rgeo(*),crfnvec(3),cqvec(3,*),rc,r0,r,
     &                 theta,x0,y0,z0,cset(3,200),x,y,qx,qy,qz,q
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      t  = rgeo(3)
      rc = rgeo(4)
      r0 = rc + t
c
      do i=1, inm
         if (nnr(i).gt.0) then
            x0 = npos(i,1)
            y0 = npos(i,2)
            z0 = npos(i,3)
c
            r = rc + t - z0
            theta = y0 / r0
            npos(i,1) =  x0
            npos(i,2) =  r*cos(theta)
            npos(i,3) = -r*sin(theta)
         endif
      enddo
c
c Define normal vector of crack face
c
      crfnvec(1) =  0.0
      crfnvec(2) =  0.0
      crfnvec(3) =  1.0
c
c Calculate the normal vector of crack front
c
      call ncrd_ntip(nset,cset,ncq, nods,ns1,ns2,ns3)
cc      call calc_crack_front_normal(cqvec,cset,ncq,-1,-1,0)
      do i=1, ncq
         x = cset(1,i)
         y = cset(2,i)
         qx = 2.*x/(c*c)
         qy = 2.*(y-rc-t)/(a*a)
         qz = 0.
         q  = sqrt(qx*qx + qy*qy + qz*qz)
         cqvec(1,i) = qx / q
         cqvec(2,i) = qy / q
         cqvec(3,i) = qz / q
      enddo
c
c Change local node numbering order in elements
c
      call reorder_elnodes(elnum,etyp)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine geo_cyl_203(rgeo,crfnvec,cqvec,ncq,elnum,etyp,
     &                       ns1,ns2,ns3,nods)
c
c Cylinder: internal circumferential crack
c
c - Reorder local node numbering order in elements
c
      implicit none
      include 'plate_common_nod.f'
      integer  elnum,etyp, ns1,ns2,ns3,nods(ns1,ns2,ns3), ncq
      integer  nset(200),i
      double precision rgeo(*),crfnvec(3),cqvec(3,*),rc,r0,r,theta,
     1                 x0,y0,z0,cset(3,200),y,z,qx,qy,qz,q,r2,
     2                 dfdy1,dfdz1,dy1dy2,dy1dz2,dz1dy2,dz1dz2,
     3                 calc_arctan
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      rc = rgeo(4)
      r0 = rc
c
      do i=1, inm
         if (nnr(i).gt.0) then
            x0 = npos(i,1)
            y0 = npos(i,2)
            z0 = npos(i,3)
c
            r = rc + z0
            theta = x0 / r0
            npos(i,1) =  y0
            npos(i,2) =  r*cos(theta)
            npos(i,3) = -r*sin(theta)
         endif
      enddo
c
c Define normal vector of crack face
c
      crfnvec(1) = -1.0
      crfnvec(2) =  0.0
      crfnvec(3) =  0.0
c
c Calculate the normal vector of crack front
c
      call ncrd_ntip(nset,cset,ncq, nods,ns1,ns2,ns3)
cc      call calc_crack_front_normal(cqvec,cset,ncq,0,1,1)
      do i=1, ncq
         y = cset(2,i)
         z = cset(3,i)
         r2 = y*y + z*z
         dfdy1 = 2.*(sqrt(r2)-rc) / (a*a)
         dfdz1 = 2.*rc*calc_arctan(2,z,y) / (c*c)
         dy1dy2 = y/sqrt(r2)
         dy1dz2 = z/sqrt(r2)
         dz1dy2 = -rc*z/r2
         dz1dz2 =  rc*y/r2
         qx = 0.
         qy = dfdy1*dy1dy2 + dfdz1*dz1dy2
         qz = dfdy1*dy1dz2 + dfdz1*dz1dz2
         q  = sqrt(qx*qx + qy*qy + qz*qz)
         cqvec(1,i) = qx / q
         cqvec(2,i) = qy / q
         cqvec(3,i) = qz / q
      enddo
c
c Change local node numbering order in elements
c
      call reorder_elnodes(elnum,etyp)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine geo_cyl_204(rgeo,crfnvec,cqvec,ncq,ns1,ns2,ns3,nods)
c
c Cylinder: external circumferential crack
c
c - Reorder local node numbering order in elements
c
      implicit none
      include 'plate_common_nod.f'
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3), ncq
      integer  nset(200),i
      double precision rgeo(*),crfnvec(3),cqvec(3,*),rc,r0,r,theta,
     1                 x0,y0,z0,cset(3,200),y,z,qx,qy,qz,q,r2,
     2                 dfdy1,dfdz1,dy1dy2,dy1dz2,dz1dy2,dz1dz2,
     3                 calc_arctan
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      t  = rgeo(3)
      rc = rgeo(4)
      r0 = rc + t
c
      do i=1, inm
         if (nnr(i).gt.0) then
            x0 = npos(i,1)
            y0 = npos(i,2)
            z0 = npos(i,3)
c
            r = rc + t - z0
            theta = x0 / r0
            npos(i,1) =  y0
            npos(i,2) =  r*cos(theta)
            npos(i,3) = -r*sin(theta)
         endif
      enddo
c
c Define normal vector of crack face
c
      crfnvec(1) = -1.0
      crfnvec(2) =  0.0
      crfnvec(3) =  0.0
c
c Calculate the normal vector of crack front
c
      call ncrd_ntip(nset,cset,ncq, nods,ns1,ns2,ns3)
cc      call calc_crack_front_normal(cqvec,cset,ncq,0,-1,-1)
      do i=1, ncq
         y = cset(2,i)
         z = cset(3,i)
         r2 = y*y + z*z
         dfdy1 = 2.*(sqrt(r2)-(rc+t))/(a*a)
         dfdz1 = 2.*(rc+t)*calc_arctan(2,z,y) / (c*c)
         dy1dy2 = y/sqrt(r2)
         dy1dz2 = z/sqrt(r2)
         dz1dy2 = -(rc+t)*z/r2
         dz1dz2 =  (rc+t)*y/r2
         qx = 0.
         qy = dfdy1*dy1dy2 + dfdz1*dz1dy2
         qz = dfdy1*dy1dz2 + dfdz1*dz1dz2
         q  = sqrt(qx*qx + qy*qy + qz*qz)
         cqvec(1,i) = qx / q
         cqvec(2,i) = qy / q
         cqvec(3,i) = qz / q
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine geo_elb_301(rgeo,crfnvec,cqvec,ncq,ns1,ns2,ns3,nods)
c
c Pipe-elbow: internal axial crack, location extrados
c
      implicit none
      include 'plate_common_nod.f'
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3), ncq
      integer  nset(200),i
      double precision rgeo(*),crfnvec(3),cqvec(3,*),rc,rp,r0c,r0p,r,
     1                 theta,phi,x0,y0,z0,x2,y2,z2,
     2                 cset(3,200),x,y,qx,qy,qz,q,r2,dfdx1,dfdy1,
     3                 dx1dx2,dx1dy2,dy1dx2,dy1dy2,calc_arctan
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      rc = rgeo(4)
      rp = rgeo(5)
      r0c = rc
      r0p = rp + rc
c
      do i=1, inm
         if (nnr(i).gt.0) then
            x0 = npos(i,1)
            y0 = npos(i,2)
            z0 = npos(i,3)
c
            r = rc + z0
            theta = y0 / r0c
            x2 =  x0
            y2 =  r*cos(theta)
            z2 = -r*sin(theta)
c
            r  =  rp + y2
            phi=  x2 / r0p
            npos(i,1) = r*sin(phi)
            npos(i,2) = r*cos(phi) - rp
            npos(i,3) = z2
         endif
      enddo
c
c Define normal vector of crack face
c
      crfnvec(1) = 0.0
      crfnvec(2) = 0.0
      crfnvec(3) = 1.0
c
c Calculate the normal vector of crack front
c
      call ncrd_ntip(nset,cset,ncq, nods,ns1,ns2,ns3)
cc      call calc_crack_front_normal(cqvec,cset,ncq,1,1,0)
      do i=1, ncq
         x = cset(1,i)
         y = cset(2,i)
         r2 = x*x + (y+rp)*(y+rp)
         dfdx1 = 2.*(rp+rc)*calc_arctan(1,x,y+rp) / (c*c)
         dfdy1 = 2.*( sqrt(r2) - (rp+rc) ) / (a*a)
         dx1dx2 =  (rp+rc)*(y+rp) / r2
         dx1dy2 = -(rp+rc)*x / r2
         dy1dx2 = x / sqrt(r2)
         dy1dy2 = (y+rp) / sqrt(r2)
         qx = dfdx1*dx1dx2 + dfdy1*dy1dx2
         qy = dfdx1*dx1dy2 + dfdy1*dy1dy2
         qz = 0.
         q  = sqrt(qx*qx + qy*qy + qz*qz)
         cqvec(1,i) = qx / q
         cqvec(2,i) = qy / q
         cqvec(3,i) = qz / q
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine geo_elb_302(rgeo,crfnvec,cqvec,ncq,elnum,etyp,
     &                       ns1,ns2,ns3,nods)
c
c Pipe-elbow: external axial crack, location extrados
c
c - Reorder local node numbering order in elements
c
      implicit none
      include 'plate_common_nod.f'
      integer  elnum,etyp, ns1,ns2,ns3,nods(ns1,ns2,ns3), ncq
      integer  nset(200),i
      double precision rgeo(*),crfnvec(3),cqvec(3,*),rc,rp,r0c,r0p,r,
     1                 theta,phi,x0,y0,z0,x2,y2,z2,
     2                 cset(3,200),x,y,qx,qy,qz,q,r2,dfdx1,dfdy1,
     3                 dx1dx2,dx1dy2,dy1dx2,dy1dy2,calc_arctan
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      t  = rgeo(3)
      rc = rgeo(4)
      rp = rgeo(5)
      r0c = rc + t
      r0p = rp + rc + t
c
      do i=1, inm
         if (nnr(i).gt.0) then
            x0 = npos(i,1)
            y0 = npos(i,2)
            z0 = npos(i,3)
c
            r = rc + t - z0
            theta = y0 / r0c
            x2 =  x0
            y2 =  r*cos(theta)
            z2 = -r*sin(theta)
c
            r  =  rp + y2
            phi=  x2 / r0p
            npos(i,1) = r*sin(phi)
            npos(i,2) = r*cos(phi) - rp
            npos(i,3) = z2
         endif
      enddo
c
c Define normal vector of crack face
c
      crfnvec(1) = 0.0
      crfnvec(2) = 0.0
      crfnvec(3) = 1.0
c
c Calculate the normal vector of crack front
c
      call ncrd_ntip(nset,cset,ncq, nods,ns1,ns2,ns3)
cc      call calc_crack_front_normal(cqvec,cset,ncq,-1,-1,0)
      do i=1, ncq
         x = cset(1,i)
         y = cset(2,i)
         r2 = x*x + (y+rp)*(y+rp)
         dfdx1 = 2.*(rp+rc+t)*calc_arctan(1,x,y+rp) / (c*c)
         dfdy1 = 2.*( sqrt(r2) - (rp+rc+t) ) / (a*a)
         dx1dx2 =  (rp+rc+t)*(y+rp) / r2
         dx1dy2 = -(rp+rc+t)*x / r2
         dy1dx2 = x / sqrt(r2)
         dy1dy2 = (y+rp) / sqrt(r2)
         qx = dfdx1*dx1dx2 + dfdy1*dy1dx2
         qy = dfdx1*dx1dy2 + dfdy1*dy1dy2
         qz = 0.
         q  = sqrt(qx*qx + qy*qy + qz*qz)
         cqvec(1,i) = qx / q
         cqvec(2,i) = qy / q
         cqvec(3,i) = qz / q
      enddo
c
c Change local node numbering order in elements
c
      call reorder_elnodes(elnum,etyp)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine geo_elb_303(rgeo,crfnvec,cqvec,ncq,elnum,etyp,
     &                       ns1,ns2,ns3,nods)
c
c Pipe-elbow: internal "circumferential" crack, location extrados
c
c - Reorder local node numbering order in elements
c
      implicit none
      include 'plate_common_nod.f'
      integer  elnum,etyp, ns1,ns2,ns3,nods(ns1,ns2,ns3), ncq
      integer  nset(200),i
      double precision rgeo(*),crfnvec(3),cqvec(3,*),rc,rp,r0c,r0p,
     1                 r,theta,phi,x0,y0,z0,x2,y2,z2,
     2                 cset(3,200),y,z,qx,qy,qz,q,r2,dfdy1,dfdz1,
     3                 dy1dy2,dy1dz2,dz1dy2,dz1dz2,calc_arctan
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      rc = rgeo(4)
      rp = rgeo(5)
      r0c = rc
      r0p = rp + rc
c
      do i=1, inm
         if (nnr(i).gt.0) then
            x0 = npos(i,1)
            y0 = npos(i,2)
            z0 = npos(i,3)
c
            r = rc + z0
            theta = x0 / r0c
            x2 =  y0
            y2 =  r*cos(theta)
            z2 = -r*sin(theta)
c
            r  =  rp + y2
            phi=  x2 / r0p 
            npos(i,1) = r*sin(phi)
            npos(i,2) = r*cos(phi) - rp
            npos(i,3) = z2
         endif
      enddo
c
c Calculate the normal vector of crack front
c
      call ncrd_ntip(nset,cset,ncq, nods,ns1,ns2,ns3)
cc      call calc_crack_front_normal(cqvec,cset,ncq,0,1,1)
      do i=1, ncq
         y = cset(2,i)
         z = cset(3,i)
         r2 = y*y + z*z
         dfdy1 = 2.*(sqrt(r2)-rc)/(a*a)
         dfdz1 = 2.*rc*calc_arctan(2,z,y) / (c*c)
         dy1dy2 = y/sqrt(r2)
         dy1dz2 = z/sqrt(r2)
         dz1dy2 = -rc*z/r2
         dz1dz2 =  rc*y/r2
         qx = 0.
         qy = dfdy1*dy1dy2 + dfdz1*dz1dy2
         qz = dfdy1*dy1dz2 + dfdz1*dz1dz2
         q  = sqrt(qx*qx + qy*qy + qz*qz)
         cqvec(1,i) = qx / q
         cqvec(2,i) = qy / q
         cqvec(3,i) = qz / q
      enddo
c
c Define normal vector of crack face
c
      crfnvec(1) = -1.0
      crfnvec(2) =  0.0
      crfnvec(3) =  0.0
c
c Change local node numbering order in elements
c
      call reorder_elnodes(elnum,etyp)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine geo_elb_304(rgeo,crfnvec,cqvec,ncq,ns1,ns2,ns3,nods)
c
c Pipe-elbow: external circumferential crack, location extrados 
c
c - Reorder local node numbering order in elements
c
      implicit none
      include 'plate_common_nod.f'
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3), ncq
      integer  nset(200),i
      double precision rgeo(*),crfnvec(3),cqvec(3,*),rc,rp,r0c,r0p,r,
     1                 theta,phi,x0,y0,z0,x2,y2,z2,
     2                 cset(3,200),y,z,qx,qy,qz,q,r2,dfdy1,dfdz1,
     3                 dy1dy2,dy1dz2,dz1dy2,dz1dz2,calc_arctan
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      t  = rgeo(3)
      rc = rgeo(4)
      rp = rgeo(5)
      r0c = rc + t
      r0p = rp + rc + t
c
      do i=1, inm
         if (nnr(i).gt.0) then
            x0 = npos(i,1)
            y0 = npos(i,2)
            z0 = npos(i,3)
c
            r = rc + t - z0
            theta = x0 / r0c
            x2 =  y0
            y2 =  r*cos(theta)
            z2 = -r*sin(theta)
c
            r = rp + y2
            phi = x2 / r0p
            npos(i,1) = r*sin(phi)
            npos(i,2) = r*cos(phi) - rp
            npos(i,3) = z2
         endif
      enddo
c
c Define normal vector of crack face
c
      crfnvec(1) = -1.0
      crfnvec(2) =  0.0
      crfnvec(3) =  0.0
c
c Calculate the normal vector of crack front
c
      call ncrd_ntip(nset,cset,ncq, nods,ns1,ns2,ns3)
cc      call calc_crack_front_normal(cqvec,cset,ncq,0,-1,-1)
      do i=1, ncq
         y = cset(2,i)
         z = cset(3,i)
         r2 = y*y + z*z
         dfdy1 = 2.*(sqrt(r2)-(rc+t))/(a*a)
         dfdz1 = 2.*(rc+t)*calc_arctan(2,z,y) / (c*c)
         dy1dy2 = y/sqrt(r2)
         dy1dz2 = z/sqrt(r2)
         dz1dy2 = -(rc+t)*z/r2
         dz1dz2 =  (rc+t)*y/r2
         qx = 0.
         qy = dfdy1*dy1dy2 + dfdz1*dz1dy2
         qz = dfdy1*dy1dz2 + dfdz1*dz1dz2
         q  = sqrt(qx*qx + qy*qy + qz*qz)
         cqvec(1,i) = qx / q
         cqvec(2,i) = qy / q
         cqvec(3,i) = qz / q
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine geo_elb_305(rgeo,crfnvec,cqvec,ncq,elnum,etyp,
     &                       ns1,ns2,ns3,nods)
c
c Pipe-elbow: internal axial crack, location intrados
c
      implicit none
      include 'plate_common_nod.f'
      integer  elnum,etyp, ns1,ns2,ns3,nods(ns1,ns2,ns3), ncq
      integer  nset(200),i
      double precision rgeo(*),crfnvec(3),cqvec(3,*),rc,rp,r0c,r0p,r,
     1                 theta,phi,x0,y0,z0,x2,y2,z2,
     2                 cset(3,200),x,y,qx,qy,qz,q,r2,dfdx1,dfdy1,
     3                 dx1dx2,dx1dy2,dy1dx2,dy1dy2,calc_arctan
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      rc = rgeo(4)
      rp = rgeo(5)
      r0c = rc
      r0p = rp - rc
c
      do i=1, inm
         if (nnr(i).gt.0) then
            x0 = npos(i,1)
            y0 = npos(i,2)
            z0 = npos(i,3)
c
            r = rc + z0
            theta = y0 / r0c
            x2 =  x0
            y2 =  r*cos(theta)
            z2 = -r*sin(theta)
c
            y2 = - y2
c
            r  =  rp + y2
            phi=  x2 / r0p
            npos(i,1) = r*sin(phi)
            npos(i,2) = r*cos(phi) - rp
            npos(i,3) = z2
         endif
      enddo
c
c Define normal vector of crack face
c
      crfnvec(1) = 0.0
      crfnvec(2) = 0.0
      crfnvec(3) = 1.0
c
c Calculate the normal vector of crack front
c
      call ncrd_ntip(nset,cset,ncq, nods,ns1,ns2,ns3)
cc      call calc_crack_front_normal(cqvec,cset,ncq,-1,-1,0)
      do i=1, ncq
         x = cset(1,i)
         y = cset(2,i)
         r2 = x*x + (y+rp)*(y+rp)
         dfdx1 = 2.*(rp-rc)*calc_arctan(1,x,y+rp) / (c*c)
         dfdy1 = 2.*(sqrt(r2) - (rp-rc)) / (a*a)
         dx1dx2 =  (rp-rc)*(y+rp) / r2
         dx1dy2 = -(rp-rc)*x      / r2
         dy1dx2 =      x / sqrt(r2)
         dy1dy2 = (y+rp) / sqrt(r2)
         qx = dfdx1*dx1dx2 + dfdy1*dy1dx2
         qy = dfdx1*dx1dy2 + dfdy1*dy1dy2
         qz = 0.
         q  = sqrt(qx*qx + qy*qy + qz*qz)
         cqvec(1,i) = qx / q
         cqvec(2,i) = qy / q
         cqvec(3,i) = qz / q
      enddo
c
c Change local node numbering order in elements
c
      call reorder_elnodes(elnum,etyp)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine geo_elb_306(rgeo,crfnvec,cqvec,ncq,ns1,ns2,ns3,nods)
c
c Pipe-elbow: external axial crack, location intrados
c
c - Reorder local node numbering order in elements
c
      implicit none
      include 'plate_common_nod.f'
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3), ncq
      integer  nset(200),i
      double precision rgeo(*),crfnvec(3),cqvec(3,*),rc,rp,r0c,r0p,r,
     1                 theta,phi,x0,y0,z0,x2,y2,z2,
     2                 cset(3,200),x,y,qx,qy,qz,q,r2,dfdx1,dfdy1,
     3                 dx1dx2,dx1dy2,dy1dx2,dy1dy2,calc_arctan
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      t  = rgeo(3)
      rc = rgeo(4)
      rp = rgeo(5)
      r0c = rc + t
      r0p = rp - rc - t
c
      do i=1, inm
         if (nnr(i).gt.0) then
            x0 = npos(i,1)
            y0 = npos(i,2)
            z0 = npos(i,3)
c
            r = rc + t - z0
            theta = y0 / r0c
            x2 =  x0
            y2 =  r*cos(theta)
            z2 = -r*sin(theta)
c
            y2 = - y2
c
            r  =  rp + y2
            phi=  x2 / r0p
            npos(i,1) = r*sin(phi)
            npos(i,2) = r*cos(phi) - rp
            npos(i,3) = z2
         endif
      enddo
c
c Define normal vector of crack face
c
      crfnvec(1) = 0.0
      crfnvec(2) = 0.0
      crfnvec(3) = 1.0
c
c Calculate the normal vector of crack front
c
      call ncrd_ntip(nset,cset,ncq, nods,ns1,ns2,ns3)
cc      call calc_crack_front_normal(cqvec,cset,ncq,1,1,0)
      do i=1, ncq
         x = cset(1,i)
         y = cset(2,i)
         r2 = x*x + (y+rp)*(y+rp)
         dfdx1 = 2.*(rp-rc-t)*calc_arctan(1,x,y+rp) / (c*c)
         dfdy1 = 2.*(sqrt(r2) - (rp-rc-t)) / (a*a)
         dx1dx2 =  (rp-rc-t)*(y+rp) / r2
         dx1dy2 = -(rp-rc-t)*x      / r2
         dy1dx2 =      x / sqrt(r2)
         dy1dy2 = (y+rp) / sqrt(r2)
         qx = dfdx1*dx1dx2 + dfdy1*dy1dx2
         qy = dfdx1*dx1dy2 + dfdy1*dy1dy2
         qz = 0.
         q  = sqrt(qx*qx + qy*qy + qz*qz)
         cqvec(1,i) = qx / q
         cqvec(2,i) = qy / q
         cqvec(3,i) = qz / q
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine geo_elb_307(rgeo,crfnvec,cqvec,ncq,ns1,ns2,ns3,nods)
c
c Pipe-elbow: internal "circumferential" crack, location intrados
c
c - Reorder local node numbering order in elements
c
      implicit none
      include 'plate_common_nod.f'
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3), ncq
      integer  nset(200),i
      double precision rgeo(*),crfnvec(3),cqvec(3,*),rc,rp,r0c,r0p,r,
     1                 theta,phi,x0,y0,z0,x2,y2,z2,
     2                 cset(3,200),y,z,qx,qy,qz,q,r2,dfdy1,dfdz1,
     3                 dy1dy2,dy1dz2,dz1dy2,dz1dz2,calc_arctan
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      rc = rgeo(4)
      rp = rgeo(5)
      r0c = rc
      r0p = rp - rc
c
      do i=1, inm
         if (nnr(i).gt.0) then
            x0 = npos(i,1)
            y0 = npos(i,2)
            z0 = npos(i,3)
c
            r = rc + z0
            theta = x0 / r0c
            x2 =  y0
            y2 =  r*cos(theta)
            z2 = -r*sin(theta)
c
            y2 = - y2
c
            r  =  rp + y2
            phi=  x2 / r0p 
            npos(i,1) = r*sin(phi)
            npos(i,2) = r*cos(phi) - rp
            npos(i,3) = z2
         endif
      enddo
c
c Define normal vector of crack face
c
      crfnvec(1) = -1.0
      crfnvec(2) =  0.0
      crfnvec(3) =  0.0
c
c Calculate the normal vector of crack front
c
      call ncrd_ntip(nset,cset,ncq, nods,ns1,ns2,ns3)
cc      call calc_crack_front_normal(cqvec,cset,ncq,0,-1,-1)
      do i=1, ncq
         y = cset(2,i)
         z = cset(3,i)
         r2 = y*y + z*z
         dfdy1 =   2.*( rc - sqrt(r2) ) / (a*a)
         dfdz1 = - 2.*rc*calc_arctan(1,-z,-y) / (c*c)
         dy1dy2 = - y/sqrt(r2)
         dy1dz2 = - z/sqrt(r2)
         dz1dy2 =   rc*z/r2
         dz1dz2 = - rc*y/r2
         qx = 0.
         qy = dfdy1*dy1dy2 + dfdz1*dz1dy2
         qz = dfdy1*dy1dz2 + dfdz1*dz1dz2
         q  = sqrt(qx*qx + qy*qy + qz*qz)
         cqvec(1,i) = qx / q
         cqvec(2,i) = qy / q
         cqvec(3,i) = qz / q
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine geo_elb_308(rgeo,crfnvec,cqvec,ncq,elnum,etyp,
     &                       ns1,ns2,ns3,nods)
c
c Pipe-elbow: external circumferential crack, location intrados 
c
c - Reorder local node numbering order in elements
c
      implicit none
      include 'plate_common_nod.f'
      integer  elnum,etyp, ns1,ns2,ns3,nods(ns1,ns2,ns3), ncq
      integer  nset(200),i
      double precision rgeo(*),crfnvec(3),cqvec(3,*),rc,rp,r0c,r0p,r,
     1                 theta,phi,x0,y0,z0,x2,y2,z2,
     2                 cset(3,200),y,z,qx,qy,qz,q,r2,dfdy1,dfdz1,
     3                 dy1dy2,dy1dz2,dz1dy2,dz1dz2,calc_arctan
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      t  = rgeo(3)
      rc = rgeo(4)
      rp = rgeo(5)
      r0c = rc + t
      r0p = rp - rc - t
c
      do i=1, inm
         if (nnr(i).gt.0) then
            x0 = npos(i,1)
            y0 = npos(i,2)
            z0 = npos(i,3)
c
            r = rc + t - z0
            theta = x0 / r0c
            x2 =  y0
            y2 =  r*cos(theta)
            z2 = -r*sin(theta)
c
            y2 = - y2
c
            r = rp + y2
            phi = x2 / r0p
            npos(i,1) = r*sin(phi)
            npos(i,2) = r*cos(phi) - rp
            npos(i,3) = z2
         endif
      enddo
c
c Define normal vector of crack face
c
      crfnvec(1) = -1.0
      crfnvec(2) =  0.0
      crfnvec(3) =  0.0
c
c Calculate the normal vector of crack front
c
      call ncrd_ntip(nset,cset,ncq, nods,ns1,ns2,ns3)
cc      call calc_crack_front_normal(cqvec,cset,ncq,0,1,1)
      do i=1, ncq
         y = cset(2,i)
         z = cset(3,i)
         r2 = y*y + z*z
         dfdy1 =  2.*(rc+t - sqrt(r2))/(a*a)
         dfdz1 = -2.*(rc+t)*calc_arctan(1,-z,-y) / (c*c)
         dy1dy2 = - y/sqrt(r2)
         dy1dz2 = - z/sqrt(r2)
         dz1dy2 =   (rc+t)*z/r2
         dz1dz2 = - (rc+t)*y/r2
         qx = 0.
         qy = dfdy1*dy1dy2 + dfdz1*dz1dy2
         qz = dfdy1*dy1dz2 + dfdz1*dz1dz2
         q  = sqrt(qx*qx + qy*qy + qz*qz)
         cqvec(1,i) = qx / q
         cqvec(2,i) = qy / q
         cqvec(3,i) = qz / q
      enddo
c
c Change local node numbering order in elements
c
      call reorder_elnodes(elnum,etyp)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function calc_arctan(case,alpha,beta)
      implicit none
      integer  case
      double precision alpha,beta,rho,one,theta,pi
      parameter (one=1.)
c
c  Calcultates theta = arctan(alpha/beta) according to a left hand
c  oriented cartesian coordinate system.
c
c    case = 1   0 < theta < pi;    case = 2  -pi < theta < 0
c
      pi = 4.*atan(one)
      rho = sqrt(alpha*alpha + beta*beta)
      if (rho.eq.0) then
         calc_arctan = 0.
         return
      endif
      if (case.eq.1) then
         if (alpha.lt.0) then
            if (beta.gt.0)  theta = 0.
            if (beta.lt.0)  theta = pi
         elseif (beta.gt.0) then
            theta = acos(beta/rho)
         elseif (beta.lt.0) then
            theta = pi - acos(beta/rho)
         else
            theta = pi/2.
         endif
      elseif (case.eq.2) then
         if (alpha.gt.0) then
            if (beta.gt.0)  theta =   0.
            if (beta.lt.0)  theta = - pi
         elseif (beta.gt.0) then
            theta = - acos(beta/rho)
         elseif (beta.lt.0) then
            theta = - ( pi - acos(beta/rho) )
         else
            theta = - pi/2.
         endif
      else
         theta = 0.
      endif
      calc_arctan = theta
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine reorder_elnodes(elnum,etyp)
c
c Reorder the local node number order. Necessary if coordinate
c transformation involves mirroring (to avoid a neg. Jacobian).
c
c Currently only 8-node bricks is supported.
c
      implicit none
      include 'plate_common_eln.f'
      integer  elnum,etyp, i,j,n(27)
c
      do i=1, elnum
         do j=1, etyp
            n(j) = eln(i,j)
         enddo
c         eln(i,1) = n(1)         
         eln(i,2) = n(4)         
c         eln(i,3) = n(3)         
         eln(i,4) = n(2)         
c         eln(i,5) = n(5)         
         eln(i,6) = n(8)         
c         eln(i,7) = n(7)         
         eln(i,8) = n(6)         
      enddo
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine calc_crack_front_normal(qnorm,coord,nc,sx,sy,sz)
      implicit none
      integer  nc,i,sx,sy,sz
      double precision qnorm(3,*),coord(3,*),
     &                 tmp(3,200),dx,dy,dz,dl,w1,w2, sn(200),s
      logical  wfront
c
      if (sz.eq.0) then
         do i=2, nc
            dx = coord(1,i) - coord(1,i-1)
            dy = coord(2,i) - coord(2,i-1)
            dl = sqrt(dx*dx+dy*dy)
            tmp(1,i) = -dy / dl
            tmp(2,i) =  dx / dl
            tmp(3,i) =  dl
         enddo
         qnorm(1,1) = tmp(1,2)
         qnorm(2,1) = tmp(2,2)
         qnorm(3,1) = 0.
         qnorm(1,nc) = tmp(1,nc)
         qnorm(2,nc) = tmp(2,nc)
         qnorm(3,nc) = 0.
         do i=2, nc-1
            w1 = tmp(3,i+1) / ( tmp(3,i)+tmp(3,i+1) )
            w2 = tmp(3,i)   / ( tmp(3,i)+tmp(3,i+1) )
            qnorm(1,i) = tmp(1,i)*w1 + tmp(1,i+1)*w2
            qnorm(2,i) = tmp(2,i)*w1 + tmp(2,i+1)*w2
            qnorm(3,i) = 0.
         enddo
      elseif (sy.eq.0) then
         write(*,*) '>> Not yet implemented !!'
         stop
      elseif (sx.eq.0) then
         do i=2, nc
            dy = coord(2,i) - coord(2,i-1)
            dz = coord(3,i) - coord(3,i-1)
            dl = sqrt(dy*dy+dz*dz)
            tmp(1,i) = -dz / dl
            tmp(2,i) =  dy / dl
            tmp(3,i) =  dl
         enddo
         qnorm(1,1) = 0.
         qnorm(2,1) = tmp(1,2)
         qnorm(3,1) = tmp(2,2)
         qnorm(1,nc) = 0.
         qnorm(2,nc) = tmp(1,nc)
         qnorm(3,nc) = tmp(2,nc)
         do i=2, nc-1
            w1 = tmp(3,i+1) / ( tmp(3,i)+tmp(3,i+1) )
            w2 = tmp(3,i)   / ( tmp(3,i)+tmp(3,i+1) )
            qnorm(1,i) = 0.
            qnorm(2,i) = tmp(1,i)*w1 + tmp(1,i+1)*w2
            qnorm(3,i) = tmp(2,i)*w1 + tmp(2,i+1)*w2
         enddo
      endif
c
      do i=1, nc
         qnorm(1,i) = sign( qnorm(1,i), sx*qnorm(1,i) )
         qnorm(2,i) = sign( qnorm(2,i), sy*qnorm(2,i) )
         qnorm(3,i) = sign( qnorm(3,i), sz*qnorm(3,i) )
ccccc         write(55,'(i3,3g15.6)') i,(real(qnorm(j,i)),j=1,3)
      enddo
c
c  Calculate a normalized crack front coordinate
c
c      wfront = .true.
      wfront = .false.
      if (wfront) then
         open(unit=49,file='plate3d_cfront.dat',status='unknown')
         s = 0
         sn(1) = s
         do i=2, nc
            dx = coord(1,i) - coord(1,i-1)
            dy = coord(2,i) - coord(2,i-1)
            dz = coord(3,i) - coord(3,i-1)
            s = s + sqrt(dx*dx+dy*dy+dz*dz)
            sn(i) = s
         enddo
         write(21,'(t1,i5)') nc
         do i=1, nc
            sn(i) = sn(i) / sn(nc)
            write(21,'(t1,i5,4g15.6)') i, coord(1,i), coord(2,i),
     &                                    coord(3,i), sn(i)
         enddo
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
