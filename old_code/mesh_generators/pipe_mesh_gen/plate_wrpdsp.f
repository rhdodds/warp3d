c----67--1---------2---------3---------4---------5---------6---------712
c
c    File: plate_wrpdsp.f   Written July 2 1998
c                           Modified on July 10 1998
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_fix_disp(analysis,tcr,tge,id,rd, nods,ns1,ns2,ns3,
     &           noda,na1,na2,na3, nodb,nb1,nb2,nb3, iws,errorFlag)
c
c Fixed Displacement Boundary Conditions
c
      implicit none
c
      integer  analysis,tcr,tge,id(3,*),iws,errorFlag
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      double precision rd(3,*)
c
c TGE = 101  (plate)
      if ( tge .eq. 101 ) then
         call wrp_fix_disp_tge101(analysis,tcr,id,rd, iws,errorFlag,
     &        nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c TGE = 201 - 204 (cylinder)
      elseif ( (tge.ge.201) .and. (tge.le.204) ) then
         call wrp_fix_disp_tge200(analysis,tge,id,rd, iws,errorFlag,
     &        nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c TGE = 301 - 308 (pipe-elbow)
      elseif ( (tge.ge.301) .and. (tge.le.308) ) then
         call wrp_fix_disp_tge300(analysis,tge,id,rd, iws,errorFlag,
     &        nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      else
         write(iws,*)'>>ERROR: in "SUBR. wrp_pre_disp"'
         write(iws,*)'>>while writing fixed displacement B.C.'
         write(iws,*)'>>Unknown value for TGE geometry parameter:'
         write(iws,*)'>>   TGE =',tge
         errorFlag = 1
         return
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_fix_disp_tge101(analysis,tcr,id,rd,iws,errorFlag,
     &           nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c Fixed Displacement Boundary Conditions TGE = 101 (basic plate)
c
      implicit none
      include 'plate_common_nod.f'
      integer  analysis,tcr,id(3,*), iws,errorFlag
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      double precision rd(3,*)
c
      integer    nsize
      parameter (nsize=15000)
      integer  nfront(nsize),nf,nasym(nsize),na,
     &         ncsym(nsize),nc, ntcoi(nsize),ntc,
     &         nset(1),n,dofx(3),dofy(3),dofz(3)
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
      data   dofx(1),dofx(2),dofx(3) / 1, 0, 0 /
      data   dofy(1),dofy(2),dofy(3) / 0, 1, 0 /
      data   dofz(1),dofz(2),dofz(3) / 0, 0, 1 /
c
c TCR = 1 (embeddied); = 2,3 (semielliptical); = 4 (edge)
c
c Note that NSET=ntcoi should only be constrained in the
c Linear Elastic case  - NO LONGER USED, SEE BELOW
c
c Same crack tip mesh in WARP3D as in ABAQUS 1 April 1999
c
ccc      if (analysis.eq.1) 
ccc     &   call nset_ntip_coincident(ntcoi,ntc,nods,ns1,ns2,ns3)
      call nset_nfront(nfront,nf,nods,ns1,ns2,ns3,noda,na1,na2,na3,
     &                 nodb,nb1,nb2,nb3)
      call nset_nasym(nasym,na,nods,ns1,ns2,ns3,noda,na1,na2,na3)
      call nset_ncsym(ncsym,nc,nods,ns1,ns2,ns3,noda,na1,na2,na3,
     &                nodb,nb1,nb2,nb3)
c
      if     (tcr.eq.1) then
ccc         if (analysis.eq.1) call wrp_zero_disp(ntcoi,ntc,dofz,id,rd)
         call wrp_zero_disp(nfront,nf,dofz,id,rd)
         call wrp_zero_disp(nasym, na,dofx,id,rd)
         call wrp_zero_disp(ncsym, nc,dofy,id,rd)
      elseif (tcr.eq.2) then
ccc         if (analysis.eq.1) call wrp_zero_disp(ntcoi,ntc,dofz,id,rd)
         call wrp_zero_disp(nfront,nf,dofz,id,rd)
         call wrp_zero_disp(nasym, na,dofx,id,rd)
         nset(1) = nnr(noda(ima,1,1))
         n = 1
         call wrp_zero_disp(nset,n,dofy,id,rd)
      elseif (tcr.eq.3) then
ccc         if (analysis.eq.1) call wrp_zero_disp(ntcoi,ntc,dofz,id,rd)
         call wrp_zero_disp(nfront,nf,dofz,id,rd)
         call wrp_zero_disp(ncsym, nc,dofy,id,rd)
         nset(1) = nnr(noda(1,1,1))
         n = 1
         call wrp_zero_disp(nset,n,dofx,id,rd)
      elseif (tcr.eq.4) then
ccc         if (analysis.eq.1) call wrp_zero_disp(ntcoi,ntc,dofz,id,rd)
         call wrp_zero_disp(nfront,nf,dofz,id,rd)
         nset(1) = nnr(noda(1,1,1))
         n = 1
         call wrp_zero_disp(nset,n,dofx,id,rd)
         call wrp_zero_disp(nset,n,dofy,id,rd)
         nset(1) = nnr(noda(ima,1,1))
         n = 1
         call wrp_zero_disp(nset,n,dofx,id,rd)
      else
         write(iws,*) '>> ERROR: in SUBR. "wrp_fix_disp_tge101"'
         errorFlag = 1
         return
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_fix_disp_tge200(analysis,tge,id,rd, iws,errorFlag,
     &           nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c Fixed Displacement Boundary Conditions TGE = 201-204 (Cylider)
c
      implicit none
      include 'plate_common_nod.f'
      integer  analysis,tge,id(3,*), iws,errorFlag
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      double precision rd(3,*)
c
      integer    nsize
      parameter (nsize=15000)
      integer  nfront(nsize),nf, nasym(nsize),na,  nrem(nsize),nr,
     &         nsid(nsize),ns,   ntcoi(nsize),ntc, nset(1),n,
     &         dofx(3),dofy(3),dofz(3)
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
      data   dofx(1),dofx(2),dofx(3) / 1, 0, 0 /
      data   dofy(1),dofy(2),dofy(3) / 0, 1, 0 /
      data   dofz(1),dofz(2),dofz(3) / 0, 0, 1 /
c
c Same crack tip mesh in WARP3D as in ABAQUS 1 April 1999
c
ccc      if (analysis.eq.1) 
ccc     &   call nset_ntip_coincident(ntcoi,ntc,nods,ns1,ns2,ns3)
      call nset_nfront(nfront,nf,nods,ns1,ns2,ns3,noda,na1,na2,na3,
     &                 nodb,nb1,nb2,nb3)
      call nset_nasym(nasym,na,nods,ns1,ns2,ns3,noda,na1,na2,na3)
      call nset_nrem(nrem,nr,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
      call nset_nsid(nsid,ns,nodb,nb1,nb2,nb3)
      n = 1
c
      if     (tge.eq.201) then
ccc         if (analysis.eq.1) call wrp_zero_disp(ntcoi,ntc,dofz,id,rd)
         call wrp_zero_disp(nfront,nf,dofz,id,rd)
         call wrp_zero_disp(nasym, na,dofx,id,rd)
         call wrp_zero_disp(nrem,  nr,dofz,id,rd)
         nset(1) = nnr(noda(ima,1,kma))
         call wrp_zero_disp(nset,  n, dofy,id,rd)
      elseif (tge.eq.202) then
ccc         if (analysis.eq.1) call wrp_zero_disp(ntcoi,ntc,dofz,id,rd)
         call wrp_zero_disp(nfront,nf,dofz,id,rd)
         call wrp_zero_disp(nasym, na,dofx,id,rd)
         call wrp_zero_disp(nrem,  nr,dofz,id,rd)
         nset(1) = nnr(noda(1,1,kma))
         call wrp_zero_disp(nset,  n, dofy,id,rd)
      elseif (tge.eq.203) then
ccc         if (analysis.eq.1) call wrp_zero_disp(ntcoi,ntc,dofx,id,rd)
         call wrp_zero_disp(nfront,nf,dofx,id,rd)
         call wrp_zero_disp(nasym, na,dofz,id,rd)
         call wrp_zero_disp(nsid,  ns,dofz,id,rd)
         nset(1) = nnr(nodb(imb,1,1))
         call wrp_zero_disp(nset,  n, dofy,id,rd)
      elseif (tge.eq.204) then
ccc         if (analysis.eq.1) call wrp_zero_disp(ntcoi,ntc,dofx,id,rd)
         call wrp_zero_disp(nfront,nf,dofx,id,rd)
         call wrp_zero_disp(nasym, na,dofz,id,rd)
         call wrp_zero_disp(nsid,  ns,dofz,id,rd)
         nset(1) = nnr(nodb(imb,jmb,1))
         call wrp_zero_disp(nset,  n, dofy,id,rd)
      else
         write(iws,*) '>> ERROR: in SUBR. "wrp_fix_disp_tge101"'
         errorFlag = 1
         return
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_fix_disp_tge300(analysis,tge,id,rd, iws,errorFlag,
     &           nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c Fixed Displacement Boundary Conditions TGE = 301-308 (Pipe elbow)
c
      implicit none
      include 'plate_common_nod.f'
      integer  analysis,tge,id(3,*), iws,errorFlag
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      double precision rd(3,*)
c
      integer    nsize
      parameter (nsize=15000)
      integer  nfront(nsize),nf,nasym(nsize),na,nrem(nsize),nr,
     &         nsid(nsize),ns,nset(1),n,ntcoi(nsize),ntc,
     &         dofx(3),dofy(3),dofz(3)
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
      data   dofx(1),dofx(2),dofx(3) / 1, 0, 0 /
      data   dofy(1),dofy(2),dofy(3) / 0, 1, 0 /
      data   dofz(1),dofz(2),dofz(3) / 0, 0, 1 /
c
c Same crack tip mesh in WARP3D as in ABAQUS 1 April 1999
c
ccc      if (analysis.eq.1) 
ccc     &   call nset_ntip_coincident(ntcoi,ntc,nods,ns1,ns2,ns3)
      call nset_nfront(nfront,nf,nods,ns1,ns2,ns3,noda,na1,na2,na3,
     &                 nodb,nb1,nb2,nb3)
      call nset_nasym(nasym,na,nods,ns1,ns2,ns3,noda,na1,na2,na3)
      call nset_nrem(nrem,nr,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
      call nset_nsid(nsid,ns,nodb,nb1,nb2,nb3)
      n = 1
c
      if     ( (tge.eq.301) .or. (tge.eq.305) ) then
ccc         if (analysis.eq.1) call wrp_zero_disp(ntcoi,ntc,dofz,id,rd)
         call wrp_zero_disp(nfront,nf,dofz,id,rd)
         call wrp_zero_disp(nasym, na,dofx,id,rd)
         call wrp_zero_disp(nrem,  nr,dofz,id,rd)
         nset(1) = nnr(noda(ima,1,kma))
         call wrp_zero_disp(nset,  n, dofy,id,rd)
      elseif ( (tge.eq.302) .or. (tge.eq.306) ) then
ccc         if (analysis.eq.1) call wrp_zero_disp(ntcoi,ntc,dofz,id,rd)
         call wrp_zero_disp(nfront,nf,dofz,id,rd)
         call wrp_zero_disp(nasym, na,dofx,id,rd)
         call wrp_zero_disp(nrem,  nr,dofz,id,rd)
         nset(1) = nnr(noda(1,1,kma))
         call wrp_zero_disp(nset,  n, dofy,id,rd)
      elseif ( (tge.eq.303) .or. (tge.eq.307) ) then
ccc         if (analysis.eq.1) call wrp_zero_disp(ntcoi,ntc,dofx,id,rd)
         call wrp_zero_disp(nfront,nf,dofx,id,rd)
         call wrp_zero_disp(nasym, na,dofz,id,rd)
         call wrp_zero_disp(nsid,  ns,dofz,id,rd)
         nset(1) = nnr(nodb(imb,1,1))
         call wrp_zero_disp(nset,  n, dofy,id,rd)
      elseif ( (tge.eq.304) .or. (tge.eq.308) ) then
ccc         if (analysis.eq.1) call wrp_zero_disp(ntcoi,ntc,dofx,id,rd)
         call wrp_zero_disp(nfront,nf,dofx,id,rd)
         call wrp_zero_disp(nasym, na,dofz,id,rd)
         call wrp_zero_disp(nsid,  ns,dofz,id,rd)
         nset(1) = nnr(nodb(imb,jmb,1))
         call wrp_zero_disp(nset,  n, dofy,id,rd)
      else
         write(iws,*) '>> ERROR: in SUBR. "wrp_fix_disp_tge101"'
         errorFlag = 1
         return
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_zero_disp(nset,n,idof,id,rd)
c
c Define zero displacement conditions at the DOF's listed in nset,idof
c
      implicit none
      integer  nset(*),n,idof(3),id(3,*),i,node
      double precision rd(3,*)
c
      do i=1, n
         node = nset(i)
         if (idof(1).eq.1) then
            id(1,node) = idof(1)
            rd(1,node) = 0
         endif
         if (idof(2).eq.1) then
            id(2,node) = idof(2)
            rd(2,node) = 0
         endif
         if (idof(3).eq.1) then
            id(3,node) = idof(3)
            rd(3,node) = 0
         endif
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_pre_disp(tcr,tge,rgeo,rload,id,rd, 
     &           noda,na1,na2,na3, nodb,nb1,nb2,nb3, iws,errorFlag)
c
c TBC = 1 : Prescribed Displacement Boundary Conditions
c
      implicit none
c
      integer  tcr,tge,id(3,*),iws,errorFlag
      integer  na1,na2,na3,noda(na1,na2,na3),
     &         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      double precision rgeo(*),rload(*),rd(3,*)
c
c TGE = 101  (plate)
      if ( tge .eq. 101 ) then
         call wrp_pre_disp_tge101(tcr,rgeo,rload,id,rd,
     &        noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c TGE = 201 - 204 (cylinder)
      elseif ( (tge.ge.201) .and. (tge.le.204) ) then
         call wrp_pre_disp_tge200(tge,rgeo,rload,id,rd,
     &        noda,na1,na2,na3,nodb,nb1,nb2,nb3)
      else
         write(iws,*)'>>ERROR: in "SUBR. wrp_pre_disp"'
         errorFlag = 1
         return
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_pre_disp_tge101(tcr,rgeo,rload,id,rd,
     &           noda,na1,na2,na3,nodb,nb1,nb2,nb3)
c
      implicit none
      include 'plate_common_nod.f'
      integer  tcr,id(3,*),na1,na2,na3,noda(na1,na2,na3),
     &                     nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      double precision rgeo(*),rload(*),rd(3,*),dsp(10)
c
      integer    nsize
      parameter (nsize=15000)
      integer	       nset(nsize),n,idof(3),i
      double precision cset(3,nsize),xref(3),t,w,length
c
      do i=1, 6
         dsp(i) = rload(i)
      enddo
      call ncrd_nrem(nset,cset,n,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
      idof(1) = 0
      idof(2) = 0
      idof(3) = 1
      t = rgeo(3)
      w = rgeo(1)
      length = rgeo(2)
c			
c     TBC = 1 (embeddied); = 2,3 (semielliptical); = 4 (edge)
c
      if     (tcr.eq.1) then
         xref(1) = 0. 
         xref(2) = 0.
         xref(3) = -length
      elseif (tcr.eq.2) then
         xref(1) = 0. 
         xref(2) = t/2.
         xref(3) = -length
      elseif (tcr.eq.3) then
         xref(1) = w/2.
         xref(2) = 0.
         xref(3) = -length
      elseif (tcr.eq.4) then
         xref(1) = w/2.
         xref(2) = t/2.
         xref(3) = -length
      endif
c
      call wrp_calc_disp(nset,cset,n,idof,id,rd,xref,dsp)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_pre_disp_tge200(tge,rgeo,rload,id,rd,
     &           noda,na1,na2,na3,nodb,nb1,nb2,nb3)
c
c Prescribed Displacement Boundary Conditions TGE = 201-204 (cylinder)
c
      implicit none
      include 'plate_common_nod.f'
      integer  tge,id(3,*),na1,na2,na3,noda(na1,na2,na3),
     &                     nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      double precision rgeo(*),rload(*),rd(3,*),dsp(10)
c
      integer    nsize
      parameter (nsize=15000)
      integer	   nset(nsize),n,idof(3),i
      double precision cset(3,nsize),xref(3)
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
      do i=1, 6
         dsp(i) = rload(i)
      enddo
      call ncrd_nrem(nset,cset,n,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
      idof(1) = 1
      idof(2) = 0
      idof(3) = 0
c
c  The length of the cylinder is stored in rgeo(5)
c
      xref(1) = rgeo(5) 
      xref(2) = 0.
      xref(3) = 0.
c			
      if     ( (tge.eq.201).and.(tge.eq.202) ) then
         call ncrd_nsid(nset,cset,n,nodb,nb1,nb2,nb3)
      elseif ( (tge.eq.203).and.(tge.eq.204) ) then
         call ncrd_nrem(nset,cset,n,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
      endif
c
      call wrp_calc_disp(nset,cset,n,idof,id,rd,xref,dsp)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_calc_disp(nset,cset,n,idof,id,rd,xref,dsp)
c
c The routine calculates the prescribed displacements acording to
c
c   Ux = dx - (Y-Yr)*wz + (Z-Zr)*wy
c   Uy = dy - (Z-Zr)*wx + (X-Xr)*wz
c   Uz = dz - (X-Xr)*wy + (Y-Yr)*wx
c
c where   dx = dsp(1)   dy = dsp(2)    dz = dsp(3)
c         wx = dsp(4)   wy = dsp(5)    wz = dsp(6)
c
      implicit none
      integer          nset(*),n,idof(3),id(3,*),i,node
      double precision cset(3,*),rd(3,*),xref(*),dsp(*),
     &                 dx,dy,dz,wx,wy,wz,xr,yr,zr,x,y,z
c
      xr = xref(1)
      yr = xref(2)
      zr = xref(3)
      dx = dsp(1)
      dy = dsp(2)
      dz = dsp(3)
      wx = dsp(4)
      wy = dsp(5)
      wz = dsp(6)
      do i=1, n
         node = nset(i)
         x = cset(1,i)
         y = cset(2,i)
         z = cset(3,i)
         if (idof(1).eq.1) then
            id(1,node) = idof(1)
            rd(1,node) = dx - (y-yr)*wz + (z-zr)*wy
         endif
         if (idof(2).eq.1) then
            id(2,node) = idof(2)
            rd(2,node) = dy - (z-zr)*wx + (x-xr)*wz
         endif
         if (idof(3).eq.1) then
            id(3,node) = idof(3)
            rd(3,node) = dz - (x-xr)*wy + (y-yr)*wx
         endif
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
