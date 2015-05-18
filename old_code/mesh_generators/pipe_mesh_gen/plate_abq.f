c----67--1---------2---------3---------4---------5---------6---------712
c
c    File: plate_abq.f   Written July 1997
c                        Modified on July 1, 1998
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abaqus_inp(etyp,elnum,analysis,numLoadStep,keyhole, 
     &                      tge,tcr,tbc,rgeo,rload,
     &                      dmat,ecrsur,ecsym,erem,esid,etop,
     &                      crfnvec,cqvec,ncq,
     &                      nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                      nodb,nb1,nb2,nb3, jobname,iws,errorFlag)
c
c  Generate in-put deck to ABAQUS
c
      implicit none
      include 'plate_common_nod.f'
      include 'plate_common_eln.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      integer  etyp,elnum, analysis,numLoadStep,keyhole,
     &         tcr,tge,tbc,ncq,iws,errorFlag
      integer  ecrsur(2,*),ecsym(2,*),erem(2,*),esid(2,*),etop(2,*)
c
      double precision rgeo(*),dmat(*),rload(*),crfnvec(3),cqvec(3,*)
c
      integer   io,io2,iload(10),idof(10),nc, i,j,nstresspoints
      double precision dt,dtmin,dtmax
      character jobname*200, file0*200,file1*200,file2*200,file3*200
      character abqel*5
      logical   singlefile,ip
c
      integer      mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
      singlefile = .true.
c      singlefile = .false.
      call appendFile(jobname,'_abq.inp',file0)
      call appendFile(jobname,'_abq_crd.inp',file1)
      call appendFile(jobname,'_abq_elm.inp',file2)
      call appendFile(jobname,'_abq_equ.inp',file3)
c
      io = 14
      open( unit=io, file=file0, status='unknown' )
      write(io,101) '*HEADING'
      write(io,101) '>> 3D-FEA Surface Crack Mesh <<'
      write(io,101) '**'
      write(io,101) '*PREPRINT,MODEL=NO,ECHO=NO'
101   format(t1,a,a)
c
c Write out coordinates on file = "plate_crd.inp"
c
      if (singlefile) then
         write(io,101) '*NODE,NSET=NALL'
         io2 = io
      else
         write(io,101) '*NODE,NSET=NALL,INPUT=',file1
         io2 = io + 1
         open( unit=io2, file=file1, status='unknown' )
      endif
      do i=1, inm
         if (nnr(i).gt.0) write(io2,102) nnr(i),',',
     &      npos(i,1),',',npos(i,2),',',npos(i,3)
      enddo
      if (.not.singlefile) close(io2)
102   format(t3,i7,a,2(g18.10,a),g18.10)
c
c Write out connectivity on file = "plate_elm.inp"
c
      if (etyp.eq.8)  abqel = 'C3D8'
      if (etyp.eq.20) abqel = 'C3D20'
      if (singlefile) then
         write(io,103)'*ELEMENT,ELSET=ELALL,TYPE=',abqel
         io2 = io
      else
         write(io,103)'*ELEMENT,ELSET=ELALL,TYPE=',abqel,',INPUT=',file2
         io2 = io + 1
         open( unit=io2, file=file2, status='unknown' )
      endif
      if (etyp.eq.8) then
         do i = 1, elnum
            write(io2,104) eln(i,0),(nnr(eln(i,j)), j=1, etyp)
         enddo
      elseif (etyp.eq.20) then
         do i = 1, elnum
            write(io2,105) eln(i,0),(nnr(eln(i,j)), j=1, etyp)
         enddo
      else
         do i = 1, elnum
            write(io2,106) eln(i,0),(nnr(eln(i,j)), j=1, etyp)
         enddo
      endif
      if (.not.singlefile) close(io2)
103   format(t1,a,a,a,a)
104   format( t1,i5,',',7(i6,','),i6)
105   format( t1,i5,',',10(i6,',')/t2,9(i6,','),i6)
106   format( t1,i5,',',10(i6,',')/t2,10(i6,',')/t2,9(i6,','),i6)
c    
c Define Material properties
c
      write(io,101) '*SOLID SECTION,ELSET=ELALL,MATERIAL=MAT1'
      write(io,101) '*MATERIAL,NAME=MAT1'
c
c  analysis = 1   static, linear elastic
c        dmat(1) = no of material constants defined = 2
c        dmat(2) = Young's modulii
c        dmat(3) = Poisson's ratio
      if (analysis.eq.1) then
         write(io,101) '*ELASTIC'
         write(io,107) dmat(2),',',dmat(3)
107      format(t1,g15.8,a,f8.4)
c
c  analysis = 2  static, nonlinear elastic (deformation Plasticity)
c        dmat(1) = no of material constants defined = 4 (5 if abaqus)
c        dmat(2) = Young's modulii
c        dmat(3) = Poisson's ratio
c        dmat(4) = sigY
c        dmat(5) = n = hardening
c        dmat(6) = alpha  = offset strain (Ramberg-Osgood model)
      elseif (analysis.eq.2) then
         write(io,101) '*DEFORMATION PLASTICITY'
         write(io,108) dmat(2),',',dmat(3),',',dmat(4),',',
     &                 dmat(5),',',dmat(6)
108      format(t1,g15.8,a,f8.4,a,g15.8,a,g12.5,a,g15.8)
c
c  analysis = 3   static,  incremental elastic-plastic
c        dmat(1) = no of material constants defined >= 4
c        dmat(2) = Young's modulii
c        dmat(3) = Poisson's ratio
c        dmat(4) = sig1
c        dmat(5) = eps1
c        dmat(6) = sig2
c        dmat(7) = eps3
      elseif (analysis.eq.3) then
         write(io,101) '*ELASTIC'
         write(io,107) dmat(2),',',dmat(3)
         write(io,101) '*PLASTIC'
         nstresspoints = int( dmat(1)-2 + 0.5 )
         do i=1, nstresspoints, 2
            write(io,109) dmat(i+3),',',dmat(i+4)
         enddo
 109     format(t1,g15.8,a,g15.8)
      endif
c
c Define node sets for Domain integral definition
c
      call nset_cr1_cr1jms_(io, nods,ns1,ns2,ns3, nc,analysis,keyhole)
c
c Define Prescribed Boundary Conditions.
c
      if (tbc.eq.1) then
         if (singlefile) then
            write(io,'(t1,a,a)') '*EQUATION'
            io2 = io
         else
            write(io,'(t1,a,a)') '*EQUATION,INPUT=',file3
            io2 = io + 1
            open( unit=io2, file=file3, status='unknown' )
         endif

         call abq_pre_displacement(io2, tcr,tge, rload,iload,idof,
     &            noda,na1,na2,na3, nodb,nb1,nb2,nb3, iws,errorFlag)
         if (.not.singlefile) close(io2)
      elseif ((tbc.eq.2).or.(tbc.eq.3)) then
         call abq_pre_traction(io,tge,tcr,tbc,rgeo,rload,ip)
      endif
c
c Fixed Displacement Boundary Conditions (symmetry cond. etc.)
c
      call abq_fixed_disp(io, tcr,tge, nods,ns1,ns2,ns3, 
     &               noda,na1,na2,na3, nodb,nb1,nb2,nb3 )
c
c Node and element sets used to define out-put results
c
      write(io,101) '**'
      write(io,101) '** result - information **'
      write(io,101) '**'
      write(io,101) '*RESTART,WRITE,FREQ=999'
      write(io,101) '**'
      write(io,101) '*NSET,NSET=NODPRD'
      write(io,120)  nnr(noda(1,1,kma)),',',nnr(noda(ima,1,kma))
      write(io,101) '*ELSET,ELSET=ELPRD'
      write(io,101) ' 1,2'
120   format(t1,10(i6,a))
c
c Define load step
c
      write(io,101) '**'
      write(io,101) '** ****** s t e p s  ************'
      write(io,101) '**'
c
      dt = 1.0/real(numLoadStep)
      dtmax = 0.1*dt
      dtmin = 2.0*dt
121   format(t1,f8.6,a,f4.2,a,f8.6,a,f8.6)
c
c  analysis = 1   static, linear elastic
      if (analysis.eq.1) then
         write(io,101) '*STEP'
         write(io,101) '*STATIC'
         call abq_step_load(io,tge,tbc,rload,iload,idof,
     &                      ecrsur,ecsym,erem,esid,etop,ip)
c
c  analysis = 2  static, nonlinear elastic (deformation Plasticity)
      elseif (analysis.eq.2) then
         write(io,101) '*STEP,INC=100'
         write(io,101) '*STATIC'
c         write(io,101) '0.1, 1.0, 0.01, 0.2'
         write(io,121) dt,',',1.0,',',dtmin,',',dtmax
         call abq_step_load(io,tge,tbc,rload,iload,idof,
     &                      ecrsur,ecsym,erem,esid,etop,ip)
c
c  analysis = 3   static,  incremental elastic-plastic
      elseif (analysis.eq.3) then
         write(io,101) '*STEP,INC=100'
         write(io,101) '*STATIC'
c         write(io,101) '0.1, 1.0, 0.01, 0.2'
         write(io,121) dt,',',1.0,',',dtmin,',',dtmax
         call abq_step_load(io,tge,tbc,rload,iload,idof,
     &                      ecrsur,ecsym,erem,esid,etop,ip)
      endif
c
      write(io,101) '*NODE PRINT,NSET=NODPRD'
      write(io,101) ' U,RF'
      write(io,101) '*EL PRINT,ELSET=ELPRD'
      write(io,101) ' S,MISES'
      write(io,101) ' E'
C
      write(io,101) '*NODE FILE,NSET=NALL'
      write(io,101) ' U,RF'
      write(io,101) '*EL FILE,ELSET=ELALL'
      write(io,101) ' S'
      write(io,101) ' E'
C FEMAP requires the ASCII format to read Abaqus results
      write(io,101) '*FILE FORMAT, ASCII'
c
      call abq_def_contour(io,mr,m1,m2,crfnvec,cqvec,ncq)
c
      write(io,101) '*ENDSTEP'
      close(io)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abaqus_post(tge,tcr,tbc,dmat,nods,ns1,ns2,ns3,jobname)
c
c The routine generates an info file for
c    Post-Processing of WARP3D results
c
      implicit none
      integer          tge,tcr,tbc,ns1,ns2,ns3,nods(ns1,ns2,ns3),
     &                 io,ntip(200),nt,i,k
      parameter    (io=24)
      double precision dmat(*),xtip(3,200)
      character        jobname*200,infofile*200,space*80
c
c  Retrieve tip coordinates.
c
      call ncrd_ntip(ntip,xtip,nt, nods,ns1,ns2,ns3)
c
c Open file = "plate_crk.inf"  and
c write out the coordinates of the crack front.
c
      call appendFile(jobname,'_crf.inf',infofile)
      open( unit=io, file=infofile, status='unknown' )
      write(io,'(t1,3i5)') tge,tcr,tbc
      write(io,'(t1,2g15.8)') dmat(2),dmat(3)
      write(io,'(t1,i5)') nt
      do i=1, nt
         write(io,'(t1,i5,3g16.8)') ntip(i),(xtip(k,i),k=1,3)
      enddo
      write(io,'(a)') jobname
c
c Add a few lines with just spaces - seems to avoid problems
c when reading "job_crf.inf" file 
c
      do i=1, 80
         space(i:i) = ' '
      enddo
      write(io,'(a80)') space
      write(io,'(a80)') space
      write(io,'(a80)') space
      close(io)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_pre_displacement(io, tcr,tge, rload,iload,idof,
     &           noda,na1,na2,na3, nodb,nb1,nb2,nb3, iws,errorFlag)
c
c TBC = 1 : Prescribed Displacement Boundary Conditions
c
      implicit none
c
      integer  io, tcr,tge,  iload(*),idof(*),iws,errorFlag
      integer  na1,na2,na3,noda(na1,na2,na3),
     &         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      double precision rload(*)
c
      if ( tge .eq. 101 ) then
c
c TGE = 101  (plate)
c
         call abq_pre_disp_tge101(io, tcr, rload,iload,idof,
     &          noda,na1,na2,na3, nodb,nb1,nb2,nb3, iws,errorFlag)
      elseif ( (tge.ge.201) .and. (tge.le.204) ) then
c
c TGE = 201 - 204 (cylinder)
c
         call abq_pre_disp_tge200(io,tge,rload,iload,idof,
     &             noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      else
         write(iws,*)'>>ERROR: in "subroutine abq_bc_equation"'
         errorFlag = 1
         return
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_pre_disp_tge101(io, tcr, rload,iload,idof,
     &           noda,na1,na2,na3, nodb,nb1,nb2,nb3, iws,errorFlag)
c
c Prescribed Displacement Boundary Conditions
c
      implicit none
      include 'plate_common_nod.f'
      integer  io,tcr,iload(*),idof(*),iws,errorFlag
      integer  na1,na2,na3,noda(na1,na2,na3),
     &         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      double precision rload(*)
c
      integer    nsize
      parameter (nsize=15000)
      integer	       rdof,sdof,nset(nsize),n,n1,n2
      double precision h,r1,a1,cset(3,nsize),delta,theta
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
c     TBC = 1 (embeddied); = 2,3 (semielliptical); = 4 (edge)
c
      if ( (tcr.eq.1).or.(tcr.eq.2).or.(tcr.eq.4) ) then
         n1 = noda(1,1,kma)
         n2 = noda(ima,1,kma)
         rdof = 2
      elseif ( tcr .eq. 3 ) then
         n1 = noda(1,1,kma)
         n2 = nodb(imb,jmb,kma)
         rdof = 1
      else
         write(iws,*)'>>ERROR: TGE=101 TCR must equal 1,2,3 or 4'
         errorFlag = 1
         return		    
      endif
      iload(1) = nnr( n1 )
      iload(2) = nnr( n2 )
      r1 = npos( n1, rdof)
      h = npos( n2, rdof ) - npos( n1, rdof )
      call ncrd_nrem(nset,cset,n,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
      a1   = 1.0
      sdof = 3
      idof(1) = 3
      idof(2) = 3
      if (tcr.eq.1)                 theta =  0.
      if ((tcr.eq.2).or.(tcr.eq.4)) theta =  rload(4)
      if (tcr.eq.3)                 theta = -rload(5)
      delta = rload(3)
c
      call abq_write_equ(io,sdof,a1,r1,h,iload,idof,rdof,nset,cset,n)
c
c     Convert delta and theta to generelized displacements u1 and u2
c
c     Let   rload(1) = u1  and  rload(2) = u2
c
      rload(1) =  delta - theta*h/2.0
      rload(2) =  delta + theta*h/2.0
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_pre_disp_tge200(io, tge, rload,iload,idof,
     &                 noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c Prescribed Displacement Boundary Conditions
c TGE = 201 - 204  (cylinder)
c
      implicit none
      include 'plate_common_nod.f'
      integer  io, tge, iload(*),idof(*)
      integer  na1,na2,na3,noda(na1,na2,na3),
     &         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      double precision rload(*)
c
      integer    nsize
      parameter (nsize=15000)
      integer	         rdof,sdof,nset(nsize),n,n1,n2
      double precision h,r1,a1,cset(3,nsize),delta,theta
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
c TGE = 201
c
      if ( tge .eq. 201 ) then
         n1 = nodb(imb,1,kma)
         n2 = nodb(imb,1,1)
         call ncrd_nsid( nset,cset,n, nodb,nb1,nb2,nb3 )
c
c TGE = 202
c
      elseif ( tge .eq. 202 ) then
         n1 = nodb(imb,jmb,kma)
         n2 = nodb(imb,jmb,1)
         call ncrd_nsid( nset,cset,n, nodb,nb1,nb2,nb3 )
c
c TGE = 203  (cylinder)
c
      elseif ( tge .eq. 203 ) then
         n1 = nodb(imb,1,kma)
         n2 = noda(ima,1,kma)
         call ncrd_nrem(nset,cset,n, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c TGE = 204  (cylinder)
c
      elseif ( tge .eq. 204 ) then
         n1 = nodb(imb,jmb,kma)
         n2 = noda(1,1,kma)
         call ncrd_nrem(nset,cset,n, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      endif
c
      iload(1) = nnr( n1 )
      iload(2) = nnr( n2 )
      rdof = 2
      r1 = npos( n1, rdof )
      h = npos( n2, rdof ) - npos( n1, rdof )
      a1   = 1.0
      sdof = 1
      idof(1) = 1
      idof(2) = 1
c
      call abq_write_equ(io,sdof,a1,r1,h,iload,idof,rdof,nset,cset,n)
c
c     Convert delta and theta to generelized displacements u1 and u2
c
c     Let   rload(1) = u1  and  rload(2) = u2
c
      delta = rload(1)
      theta = rload(6)
      rload(1) =  delta + theta*h/2.0
      rload(2) =  delta - theta*h/2.0
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_write_equ(io,sdof,a1,r1,h,iload,idof,rdof,
     &                         nset,cset,n)
      implicit none
      integer          io,sdof,rdof,iload(*),idof(*),nset(*),n,snode,i
      double precision cset(3,*),r1,h,a1, r,a2,a3
c
      do i=1, n
         snode = nset(i)
         if ((snode.ne.iload(1)).and.(snode.ne.iload(2))) then
            r  = cset( rdof, i) - r1
            a2 = - ( h - r ) / h
            a3 = - r / h
            write(io,101) '3'
            write(io,102) snode,sdof,a1,
     &                    iload(1),idof(1),a2, iload(2),idof(2),a3
         endif
      enddo
101   format(t1,a)
102   format(t1,i6,',',i2,',',f6.2,2(',',i6,',',i2,',',d19.12))
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_pre_traction(io,tge,tcr,tbc,rgeo,rload,ip)
      implicit none
      integer  io,tge,tcr,tbc
      double precision rgeo(*),rload(*),ptol
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
            call abq_traction2_sub_tge101(io,tcr,rgeo,rload)
         elseif ( (tge.ge.201) .and. (tge.le.304) ) then
            call abq_traction2_sub_tge200(io,rgeo,rload,ip)
         elseif ( (tge.ge.301) .and. (tge.le.308) ) then
            call abq_traction2_sub_tge300(io,rgeo,rload,ip)
         endif
c
c  TBC = 3: Prescribed tractions on crack face
c
      elseif (tbc.eq.3) then
         ip = .false.
         if (tge.eq.101) then
            call abq_traction3_sub_tge101(io,tcr,rload)
         elseif ( (tge.ge.201) .and. (tge.le.304) ) then
            call abq_traction3_sub_tge200(io,tge,rgeo,rload)
         elseif ( (tge.ge.301) .and. (tge.le.308) ) then
            call abq_traction3_sub_tge300(io,tge,rgeo,rload)
         endif
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_traction2_sub_tge101(io,tcr,rgeo,rload)
      implicit none
c
c  Defines Traction Boundary Condition on geometry TGE = 101
c          TBC = 2: Prescribed endtractions  
c
      integer          io,tcr
      double precision rgeo(*),rload(*)
      double precision w0,t0,fz,mx,my,area,ixx,iyy,c0,c1,c2
c
      call abq_dload_subrhead(io)
      w0 = rgeo(1)
      t0 = rgeo(3)
c
      fz = rload(3)          
      mx = rload(4)          
      my = rload(5)
C
      if ( tcr .eq. 1 ) then
         area = 4.*w0*t0
         c0 =  fz / area
         write(io,101)'      F =',c0
c
      elseif (tcr.eq.2) then
         area = 2.*w0*t0
         ixx  = 2.*w0*t0**3. / 12.
         c0   = - ( -fz/area + (mx/ixx)*(t0/2.) )
         c1   =   mx/ixx
         write(io,101)'      C0 =',c0
         write(io,101)'      C1 =',c1 
         write(io,102)'      F  = C0 + C1*COORDS(2)' 
c
      elseif (tcr.eq.3) then
         area = 2.*w0*t0
         iyy  = 2.*t0*w0**3. / 12.
         c0   =  ( fz/area + (my/iyy)*(w0/2.) )
         c1   = - my/iyy
         write(io,101)'      C0 =',c0 
         write(io,101)'      C1 =',c1 
         write(io,102)'      F  = C0 + C1*COORDS(1)' 
c
      elseif (tcr.eq.4) then
         area = w0*t0
         ixx  = w0*t0**3. / 12.
         iyy  = t0*w0**3. / 12.
         c0   = - ( - fz/area + (mx/ixx)*(t0/2.) - (my/iyy)*(w0/2.) )
         c1   =   mx/ixx
         c2   = - my/iyy
         write(io,101)'      C0 =',c0 
         write(io,101)'      C1 =',c1 
         write(io,101)'      C2 =',c2 
         write(io,102)'      F  = C0 + C1*COORDS(1) + C2*COORDS(2)' 
c
      endif
c
      write(io,102)'      RETURN'
      write(io,102)'      END'
c
 101  format(t1,a,g17.10)
 102  format(t1,a,a)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_traction2_sub_tge200(io,rgeo,rload,ip)
c
c  Defines Traction Boundary Condition on geometry TGE = 201 - 204
c          TBC = 2: Prescribed endtractions  
c
      implicit none
      integer          io
      double precision rgeo(*),rload(*)
      double precision t,r0,nx,mz,area,izz,c0,c1,pi,p0
      logical          ip
c
      call abq_dload_subrhead(io)
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
      c0 = - nx/area - p0*r0 / ( (2.+t/r0)*t )
      c1 = mz/izz
c
      write(io,101)'      C0 =',c0
      write(io,101)'      C1 =',c1
      write(io,102)'      F  = C0 + C1*COORDS(2)' 
      write(io,102)'      RETURN'
      write(io,102)'      END'
c
 101  format(t1,a,g17.10)
 102  format(t1,a,a)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_traction2_sub_tge300(io,rgeo,rload,ip)
c
c  Defines Traction Boundary Condition on geometry TGE = 301 - 308
c          TBC = 2: Prescribed endtractions  
c
      implicit none
      integer          io
      double precision rgeo(*),rload(*)
      double precision t,rc,rp,n,mz,p0,area,jzz,c0,c1,pi,phi0
      logical          ip
c
      call abq_dload_subrhead(io)
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
      c0 = - (n + mz/rp ) / area
      c1 = - (mz/rp) / jzz     
c
      write(io,101)'      C  =',cos(phi0)
      write(io,101)'      S  =',sin(phi0)
      write(io,101)'      RP =',rp
      write(io,102)'      X = COORDS(1)/RP'
      write(io,102)'      Y = COORDS(2)/RP'
      write(io,102)'      RHO = X*S + Y*C - 1.0 + C '
      write(io,101)'      C0 =',c0
      write(io,101)'      C1 =',c1
      write(io,102)'      F = C0 + C1*RHO/(1.+RHO)'
      write(io,102)'      RETURN'
      write(io,102)'      END'
c
 101  format(t1,a,g17.10)
 102  format(t1,a,a)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_traction3_sub_tge101(io,tcr,rload)
      implicit none
c
c  Defines Traction Boundary Condition on geometry TGE = 101
c          TBC = 3: Prescribed tractions on crack face
c
      integer          io,tcr
      double precision rload(*)
c
      double precision t,w,c,a,alfa
      common /geom/    t,w,c,a,alfa
c
      call abq_dload_subrhead(io)
      write(io,101)'      P0 =',rload(1)
      write(io,101)'      P1 =',rload(2)
      write(io,101)'      P2 =',rload(3)
      write(io,101)'      P3 =',rload(4)
      write(io,101)'      Q1 =',rload(5)
      write(io,101)'      Q2 =',rload(6)
      write(io,101)'      Q3 =',rload(7)
      write(io,101)'      A  =',a
      write(io,101)'      C  =',c
c
      if ( (tcr.ge.1) .and. (tcr.le.4) ) then
         write(io,102)'      U = COORDS(2) / A'
         write(io,102)'      V = COORDS(1) / C'
         write(io,102)'      F = P0 + U*(P1+U*(P2+U*P3))',
     &                            ' + V*(Q1+V*(Q2+V*Q3))'
      endif
      write(io,102)'      RETURN'
      write(io,102)'      END'
c
 101  format(t1,a,g17.10)
 102  format(t1,a,a)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_traction3_sub_tge200(io,tge,rgeo,rload)
      implicit none
c
c  Defines Traction Boundary Condition on geometry TGE = 201 - 204
c          TBC = 3: Prescribed tractions on crack face
c
      integer          io,tge
      double precision rgeo(*),rload(*)
c
      double precision t,w,c,a,alfa
      common /geom/    t,w,c,a,alfa
c
      call abq_dload_subrhead(io)
      write(io,101)'      P0 =',rload(1)
      write(io,101)'      P1 =',rload(2)
      write(io,101)'      P2 =',rload(3)
      write(io,101)'      P3 =',rload(4)
      write(io,101)'      Q1 =',rload(5)
      write(io,101)'      Q2 =',rload(6)
      write(io,101)'      Q3 =',rload(7)
      write(io,101)'      A  =',a
      write(io,101)'      C  =',c
      write(io,101)'      RC =',rgeo(4)
c
      if (tge .ge. 201 ) then
         write(io,102)'      U = ( COORDS(2) - RC ) / A'
         write(io,102)'      V = COORDS(1) / C'
      elseif (tge .ge. 202 ) then
         write(io,101)'      T  =',rgeo(3)
         write(io,102)'      U = ( RC + T - COORDS(2) ) / A'
         write(io,102)'      V = COORDS(1) / C'
      elseif (tge .ge. 203 ) then
         write(io,102)'      R = SQRT( COORDS(2)*COORDS(2) + ',
     &                                 'COORDS(3)*COORDS(3) )'
         write(io,102)'      THETA = ASIN( ABS(COORDS(3)) / R ) '
         write(io,102)'      U = ( R - RC ) / A'
         write(io,102)'      V = THETA*RC / C'
      elseif (tge .ge. 204 ) then
         write(io,101)'      T  =',rgeo(3)
         write(io,102)'      R = SQRT( COORDS(2)*COORDS(2) + ',
     &                                 'COORDS(3)*COORDS(3) )'
         write(io,102)'      THETA = ASIN( ABS(COORDS(3)) / R ) '
         write(io,102)'      U = ( RC + T - RC ) / A'
         write(io,102)'      V = THETA*(RC+T) / C'
      endif
c
      write(io,102)'      F = P0 + U*(P1+U*(P2+U*P3))',
     &                         ' + V*(Q1+V*(Q2+V*Q3))'
      write(io,102)'      RETURN'
      write(io,102)'      END'
c
 101  format(t1,a,g17.10)
 102  format(t1,a,a)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_traction3_sub_tge300(io,tge,rgeo,rload)
      implicit none
c
c  Defines Traction Boundary Condition on geometry TGE = 301 - 308
c          TBC = 3: Prescribed tractions on crack face
c
      integer          io,tge
      double precision rgeo(*),rload(*)
c
      double precision t,w,c,a,alfa
      common /geom/    t,w,c,a,alfa
c
      call abq_dload_subrhead(io)
      write(io,101)'      P0 =',rload(1)
      write(io,101)'      P1 =',rload(2)
      write(io,101)'      P2 =',rload(3)
      write(io,101)'      P3 =',rload(4)
      write(io,101)'      Q1 =',rload(5)
      write(io,101)'      Q2 =',rload(6)
      write(io,101)'      Q3 =',rload(7)
      write(io,101)'      A  =',a
      write(io,101)'      C  =',c
      write(io,101)'      RC =',rgeo(4)
c
      if     (tge .ge. 301) then
         write(io,101)'      RP =',rgeo(5)
         write(io,102)'      X = COORDS(1)'
         write(io,102)'      Y = COORDS(2)'
         write(io,102)'      R  = SQRT( X*X + (Y+RP)*(Y+RP) )'
         write(io,102)'      THETA = ATAN( X / (Y+RP) )'
         write(io,102)'      CT = (Y+RP) / R'
         write(io,102)'      ST =  X / R'
         write(io,102)'      U = ( X*ST + (Y+RP)*CT - (RP+RC) ) / A'
         write(io,102)'      V = THETA*(RP+RC) / C'
      elseif (tge .ge. 302) then
         write(io,101)'      RP =',rgeo(5)
         write(io,101)'      T  =',rgeo(3)
         write(io,102)'      X = COORDS(1)'
         write(io,102)'      Y = COORDS(2)'
         write(io,102)'      R  = SQRT( X*X + (Y+RP)*(Y+RP) )'
         write(io,102)'      THETA = ATAN( X / (Y+RP) )'
         write(io,102)'      CT = (Y+RP) / R'
         write(io,102)'      ST =  X / R'
         write(io,102)'      U = -( X*ST + (Y+RP)*CT - (RP+RC+T) ) / A'
         write(io,102)'      V = THETA*(RP+RC+T) / C'
      elseif ( (tge.ge.303) .or. (tge.ge.307)) then
         write(io,102)'      R  = SQRT(COORDS(2)**2 + COORDS(3)**2)'
         write(io,102)'      THETA = ATAN(ABS(COORDS(3)/COORDS(2)))'
         write(io,102)'      U = ( R - RC ) / A'
         write(io,102)'      V = THETA*RC / C'
      elseif ( (tge.ge.304) .or. (tge.ge.308)) then
         write(io,101)'      T  =',rgeo(3)
         write(io,102)'      R  = SQRT(COORDS(2)**2 + COORDS(3)**2)'
         write(io,102)'      THETA = ATAN(ABS(COORDS(3)/COORDS(2)))'
         write(io,102)'      U = ( RC + T - R ) / A'
         write(io,102)'      V = THETA*(RC+T) / C'
      elseif (tge .ge. 305) then
         write(io,101)'      RP =',rgeo(5)
         write(io,102)'      X = COORDS(1)'
         write(io,102)'      Y = COORDS(2)'
         write(io,102)'      R  = SQRT( X*X + (Y+RP)*(Y+RP) )'
         write(io,102)'      THETA = ATAN( X / (Y+RP) )'
         write(io,102)'      CT = (Y+RP) / R'
         write(io,102)'      ST =  X / R'
         write(io,102)'      U = - ( X*ST + (Y+RP)*CT - (RP-RC) ) / A'
         write(io,102)'      V = THETA*(RP-RC) / C'
      elseif (tge .ge. 306) then
         write(io,101)'      RP =',rgeo(5)
         write(io,101)'      T  =',rgeo(3)
         write(io,102)'      X = COORDS(1)'
         write(io,102)'      Y = COORDS(2)'
         write(io,102)'      R  = SQRT( X*X + (Y+RP)*(Y+RP) )'
         write(io,102)'      THETA = ATAN( X / (Y+RP) )'
         write(io,102)'      CT = (Y+RP) / R'
         write(io,102)'      ST =  X / R'
         write(io,102)'      U = -( X*ST + (Y+RP)*CT - (RP-RC-T) ) / A'
         write(io,102)'      V = THETA*(RP-RC-T) / C'
      endif
c
      write(io,102)'      F = P0 + U*(P1+U*(P2+U*P3))',
     &                         ' + V*(Q1+V*(Q2+V*Q3))'
      write(io,102)'      RETURN'
      write(io,102)'      END'
c
 101  format(t1,a,g17.10)
 102  format(t1,a,a)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_dload_subrhead(io)
      integer   io
c
      write(io,101)'**'
      write(io,101)'*USER SUBROUTINES'
      write(io,101)'      SUBROUTINE DLOAD(F,KSTEP,KINC,TIME,NOEL,',
     &                                    'NPT,LAYER,KSPT,'
      write(io,101)'     &                 COORDS,JLTYP)'
      write(io,101)'C'
      write(io,101)'C Traction Boundary Conditions'
      write(io,101)'C'
c
      write(io,101)'      INCLUDE ''ABA_PARAM.INC'''
      write(io,101)'C'
      write(io,101)'      DIMENSION TIME(2), COORDS(3)'
      write(io,101)'C'
c
 101  format(t1,a,a,a,a)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_fixed_disp(io, tcr,tge, nods,ns1,ns2,ns3, 
     &                     noda,na1,na2,na3, nodb,nb1,nb2,nb3 )
c
c Generate node sets used to define fixed displacement B.C.
c
      implicit none
c
      integer  io, tcr,tge
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
c  TGE = 101
c
      if (tge.eq.101) then
         call abq_fixed_disp101(io, tcr, nods,ns1,ns2,ns3,
     &                 noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c  TGE = 201 - 204
c
      elseif ( (tge.ge.201) .and. (tge.le.204) ) then
         call abq_fixed_disp200(io, tge, nods,ns1,ns2,ns3,
     &                 noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c  TGE = 301 - 308
c
      elseif ( (tge.ge.301) .and. (tge.le.308) ) then
         call abq_fixed_disp300(io,tge,nods,ns1,ns2,ns3,
     &               noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_fixed_disp101(io, tcr, nods,ns1,ns2,ns3,
     &                    noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c  Fixed displacement Boundary Conditions
c  TGE = 101  (basic plate geometry)
c
      implicit none
      include 'plate_common_nod.f'
c
      integer  io, tcr, ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      integer    nsize
      parameter (nsize=15000)
      integer	   nset(nsize),n,node1,node2
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
101   format(t1,a)
102   format(t1,i6,a)
c
      write(io,101) '**'
      write(io,101) '*NSET,NSET=NFRONT'
      call nset_nfront(nset,n, nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                         nodb,nb1,nb2,nb3)
      call write_out_set(io,nset,n)
      write(io,101) '**'
      node1 = nnr( noda(1,1,1)	)
      node2 = nnr( noda(ima,1,1) )
c
c     TRC = 1: embeddied elliptical crack
c
      if (tcr.eq.1) then
         write(io,101) '*NSET,NSET=NASYM'
         call nset_nasym(nset,n,nods,ns1,ns2,ns3,noda,na1,na2,na3)
         call write_out_set(io,nset,n)
         write(io,101) '**'
         write(io,101) '*NSET,NSET=NCSYM'
         call nset_ncsym(nset,n,nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                          nodb,nb1,nb2,nb3)
         call write_out_set(io,nset,n)
         write(io,101) '**'
c
         write(io,101) '*BOUNDARY'
         write(io,101) ' NFRONT,3,3'
         write(io,101) ' NASYM,1,1'
         write(io,101) ' NCSYM,2,2'
c
c     TRC = 2: semi elliptical surface crack; X = 0 symmetry plane
c
      elseif (tcr.eq.2) then
         write(io,101) '*NSET,NSET=NASYM'
         call nset_nasym(nset,n, nods,ns1,ns2,ns3, noda,na1,na2,na3)
         call write_out_set(io,nset,n)
         write(io,101) '**'
c
         write(io,101) '*BOUNDARY'
         write(io,101) ' NFRONT,3,3'
         write(io,101) ' NASYM,1,1'
         write(io,102)  node2,',2,2'
c
c     TRC = 3: semi elliptical surface crack; Z = 0 symmetry plane
c
      elseif (tcr.eq.3) then
         write(io,101) '*NSET,NSET=NCSYM'
         call nset_ncsym(nset,n,nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                          nodb,nb1,nb2,nb3)
         call write_out_set(io,nset,n)
         write(io,101) '**'
c
         write(io,101) '*BOUNDARY'
         write(io,101) ' NFRONT,3,3'
         write(io,101) ' NCSYM,2,2'
         write(io,102)  node1,',1,1'
c
c     TRC = 4: quarter (edge) elliptical surface crack
c
      elseif (tcr.eq.4) then
         write(io,101) '*BOUNDARY'
         write(io,101) ' NFRONT,3,3'
         write(io,102)  node1,',1,1'
         write(io,102)  node2,',1,1'
         write(io,102)  node2,',2,2'
      endif
      write(io,101) '**'
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_fixed_disp200(io, tge, nods,ns1,ns2,ns3,
     &                    noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c  Fixed displacement Boundary Conditions
c  TGE = 201 - 204  ( Cylider)
c
      implicit none
      include 'plate_common_nod.f'
c
      integer  io, tge,ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      integer    nsize
      parameter (nsize=15000)
      integer	   nset(nsize),n,node1
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
101   format(t1,a)
102   format(t1,i6,a)
c
      write(io,101) '**'
      write(io,101) '*NSET,NSET=NFRONT'
      call nset_nfront(nset,n, nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                         nodb,nb1,nb2,nb3)
      call write_out_set(io,nset,n)
      write(io,101) '**'
      write(io,101) '*NSET,NSET=NASYM'
      call nset_nasym(nset,n,nods,ns1,ns2,ns3,noda,na1,na2,na3)
      call write_out_set(io,nset,n)
      write(io,101) '**'
c
c TGE = 201,202 (cylinder)
      if ( (tge.eq.201) .or. (tge.eq.202) ) then
         write(io,101) '*NSET,NSET=NREM'
         call nset_nrem(nset,n,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
         call write_out_set(io,nset,n)
         write(io,101) '**'
         write(io,101) '*BOUNDARY'
         write(io,101) ' NFRONT,3,3'
         write(io,101) ' NASYM,1,1'
         write(io,101) ' NREM,3,3'
         node1 = nnr( noda(ima,1,kma) )
         if (tge.eq.201) write(io,102) node1,',2,2'
         node1 = nnr( noda(1,1,kma) )
         if (tge.eq.202) write(io,102) node1,',2,2'
c
c TGE = 203,204 (cylinder)
      elseif ( (tge.eq.203) .or. (tge.eq.204) ) then
         write(io,101) '*NSET,NSET=NSID'
         call nset_nsid(nset,n, nodb,nb1,nb2,nb3)
         call write_out_set(io,nset,n)
         write(io,101) '**'
         write(io,101) '*BOUNDARY'
         write(io,101) ' NFRONT,1,1'
         write(io,101) ' NASYM,3,3'
         write(io,101) ' NSID,3,3'
         node1 = nnr( nodb(imb,1,1) )
         if (tge.eq.203) write(io,102) node1,',2,2'
         node1 = nnr( nodb(imb,jmb,1) )
         if (tge.eq.204) write(io,102) node1,',2,2'
      endif
      write(io,101) '**'
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_fixed_disp300(io, tge, nods,ns1,ns2,ns3,
     &                    noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c  Fixed displacement Boundary Conditions
c  TGE = 301 - 308  ( pipe-elbow)
c
      implicit none
      include 'plate_common_nod.f'
c
      integer  io, tge, ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      integer    nsize
      parameter (nsize=15000)
      integer	   nset(nsize),n,node1
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
101   format(t1,a)
102   format(t1,i6,a)
c
      write(io,101) '**'
      write(io,101) '*NSET,NSET=NFRONT'
      call nset_nfront(nset,n, nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                         nodb,nb1,nb2,nb3)
      call write_out_set(io,nset,n)
      write(io,101) '**'
      write(io,101) '*NSET,NSET=NASYM'
      call nset_nasym(nset,n,nods,ns1,ns2,ns3,noda,na1,na2,na3)
      call write_out_set(io,nset,n)
      write(io,101) '**'
c
c TGE = 301, 305 (pipe-elbow)
      if ( (tge.eq.301) .or. (tge.eq.305) ) then
         write(io,101) '*NSET,NSET=NREM'
         call nset_nrem(nset,n,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
         call write_out_set(io,nset,n)
         write(io,101) '**'
         write(io,101) '*BOUNDARY'
         write(io,101) ' NFRONT,3,3'
         write(io,101) ' NASYM,1,1'
         write(io,101) ' NREM,3,3'
         node1 = nnr( noda(ima,1,kma) )
         write(io,102)  node1,',2,2'
c
c TGE = 302, 306 (pipe-elbow)
      elseif ( (tge.eq.302) .or. (tge.eq.306) ) then
         write(io,101) '*NSET,NSET=NREM'
         call nset_nrem(nset,n,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
         call write_out_set(io,nset,n)
         write(io,101) '**'
         write(io,101) '*BOUNDARY'
         write(io,101) ' NFRONT,3,3'
         write(io,101) ' NASYM,1,1'
         write(io,101) ' NREM,3,3'
         node1 = nnr( noda(1,1,kma) )
         write(io,102)  node1,',2,2'
c
c TGE = 303, 307 (pipe-elbow)
      elseif ( (tge.eq.303) .or. (tge.eq.307) ) then
         write(io,101) '*NSET,NSET=NSID'
         call nset_nsid(nset,n, nodb,nb1,nb2,nb3)
         call write_out_set(io,nset,n)
         write(io,101) '**'
         write(io,101) '*BOUNDARY'
         write(io,101) ' NFRONT,1,1'
         write(io,101) ' NASYM,3,3'
         write(io,101) ' NSID,3,3'
         node1 = nnr( nodb(imb,1,1) )
         write(io,102)  node1,',2,2'
c
c TGE = 304, 308 (pipe-elbow)
      elseif ( (tge.eq.304) .or. (tge.eq.308) ) then
         write(io,101) '*NSET,NSET=NSID'
         call nset_nsid(nset,n, nodb,nb1,nb2,nb3)
         call write_out_set(io,nset,n)
         write(io,101) '**'
         write(io,101) '*BOUNDARY'
         write(io,101) ' NFRONT,1,1'
         write(io,101) ' NASYM,3,3'
         write(io,101) ' NSID,3,3'
         node1 = nnr( nodb(imb,jmb,1) )
         write(io,102)  node1,',2,2'
      endif	   
      write(io,101) '**'
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_step_load(io,tge,tbc,rload,iload,idof,
     &                         ecrsur,ecsym,erem,esid,etop,ip)
c
c Definition of prescribed load
c
      implicit none
      include 'plate_common_nod.f'
c
      integer  io,tge,tbc,iload(*),idof(*)
      integer  ecrsur(2,*),ecsym(2,*),erem(2,*),esid(2,*),etop(2,*)
      double precision rload(*),pressure
      logical  ip
c
c TBC = 1:  Prescribed Displacements (delta, theta)
c           (include geometries TGE: 101,201,202,203,204)
c
      if ( tbc .eq. 1 ) then
         write(io,101) '*BOUNDARY'
         write(io,102) iload(1),idof(1),idof(1),rload(1)
         write(io,102) iload(2),idof(2),idof(2),rload(2)
101      format(t1,a)
102      format(t1,i6,',',i2,',',i2,',',g16.8)
c
c TBC = 2:  Prescribed Tractions (N, M and internal pressure)
c           (include geometries TGE: 101, 201-204, 301-308 )
c
      elseif (tbc.eq.2) then
         write(io,101) '*DLOAD'
c
         if (tge.eq.101) call abq_write_endtraction(io,erem)
c
         if (tge.eq.201) call abq_write_endtraction(io,esid)
         if (tge.eq.202) call abq_write_endtraction(io,esid)
         if (tge.eq.203) call abq_write_endtraction(io,erem)
         if (tge.eq.204) call abq_write_endtraction(io,erem)
c
         if (tge.eq.301) call abq_write_endtraction(io,esid)
         if (tge.eq.302) call abq_write_endtraction(io,esid)
         if (tge.eq.303) call abq_write_endtraction(io,erem)
         if (tge.eq.304) call abq_write_endtraction(io,erem)
         if (tge.eq.305) call abq_write_endtraction(io,esid)
         if (tge.eq.306) call abq_write_endtraction(io,esid)
         if (tge.eq.307) call abq_write_endtraction(io,erem)
         if (tge.eq.308) call abq_write_endtraction(io,erem)
c
         pressure = rload(7)
         if     ( ip .and. ( (tge.eq.201).or.(tge.eq.203) ) ) then
            call abq_write_pressure(io,ecsym,pressure)
            call abq_write_pressure(io,ecrsur,pressure)
         elseif ( ip .and. ( (tge.eq.202).or.(tge.eq.204) ) ) then
            call abq_write_pressure(io,etop,pressure)
c
         elseif ( ip .and. ( (tge.eq.301).or.(tge.eq.303).or.
     &                       (tge.eq.305).or.(tge.eq.307) ) ) then
            call abq_write_pressure(io,ecsym,pressure)
            call abq_write_pressure(io,ecrsur,pressure)
         elseif ( ip .and. ( (tge.eq.302).or.(tge.eq.304).or.
     &                       (tge.eq.306).or.(tge.eq.308) ) ) then
            call abq_write_pressure(io,etop,pressure)
         endif
c
c TBC = 3:  Prescribed Tractions on the Crack face
c           (include geometries TGE: 101, 201-204, 301-308 )
c
      elseif (tbc.eq.3) then
         write(io,101) '*DLOAD'
         call abq_write_endtraction(io,ecrsur)
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_write_endtraction(io,eset)
      implicit none
      integer   io,eset(2,*),ne,i,ie,face_no
      character*4 cface(6)
      data (cface(i),i=1,6)/'P1NU','P2NU','P3NU','P5NU','P4NU','P6NU'/
c
c  Note that default element faces are according to the WARP3D scheme 
c
      ne = eset(1,1)
c
      do i=1, ne
         ie      = eset(1,i+1)
         face_no = eset(2,i+1)
         write(io,101) ie,',',cface(face_no)
      enddo
 101  format(t1,i5,a,a)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_write_pressure(io,eset,p)
      implicit none
      integer          io,eset(2,*),ne,i,ie,face_no
      double precision p
      character*2 cface(6)
      data (cface(i),i=1,6)/'P1','P2','P3','P5','P4','P6'/
c
c  Note that default element faces are according to the WARP3D scheme 
c
      ne = eset(1,1)
c
      do i=1, ne
         ie      = eset(1,i+1)
         face_no = eset(2,i+1)
         write(io,101) ie,',',cface(face_no),',',p
      enddo
 101  format(t1,i5,a,a,a,g17.10)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abq_def_contour(io,mr,m1,m2,crfnvec,cqvec,ncq)
c
c Generate constraint on nodes coincident with the crack tip nodes
c
      implicit none
      include 'plate_common_nod.f'
c
      integer          io,mr,m1,m2,ncq
      double precision crfnvec(3),cqvec(3,*)
      integer     i,i1,i2,jdom
      character*5 jco(100)
      logical     give_qvec
c
c Define names of the crack tip node sets
c
      do i=1, 99
         i1 = i/10
         i2 = mod(i,10)
         jco(i) = ' CR'//char(48+i1)//char(48+i2)
      enddo
      jdom = mr + min(m1-2,m2-1)
c
      give_qvec = .true.
c      give_qvec = .false.
c
c Note that the normal vector of the crack plane is defined as positive
c in the outward direction.
c In ABAQUS this vector should be defined as positive in the inward
c direction.
c
      if (give_qvec) then
         write(io,101) '*CONTOUR INTEGRAL,CONTOURS=',jdom,
     &              ',TYPE=J,OUTPUT=BOTH,SYMM,FREQ=1'
         do i=1, ncq
            write(io,102)  jco(i),',',
     &                     cqvec(1,i),',',cqvec(2,i),',',cqvec(3,i)
         enddo
      else
         write(io,101) '*CONTOUR INTEGRAL,CONTOURS=',jdom,
     &              ',TYPE=J,OUTPUT=BOTH,NORMAL,SYMM,FREQ=1'
         write(io,103) -crfnvec(1),',',-crfnvec(2),',',-crfnvec(3)
         write(io,104) (jco(i),i=1,ncq)
      endif
 101  format(t1,a,i2,a)
 102  format(t1,a,3(a,g16.8))
 103  format(t1,f6.2,a,f6.2,a,f6.2)
 104  format(t1,15a/t1,15a/t1,15a/t1,15a/t1,15a/t1,15a/t1,15a)
c
      return
      end	  
c
c----67--1---------2---------3---------4---------5---------6---------712
c
