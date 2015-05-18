c----67--1---------2---------3---------4---------5---------6---------712
c
c    File: plate_wrp.f   Written July 2 1998
c                        Modified on July 3 1998
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine warp3d_inp(etyp,elnum,no_of_nodes,dmat,rgeo,rload,
     &           analysis,keyhole,numLoadStep,tcr,tge,tbc,eqsolver,
     &           crfnvec,ecrsur,ecrtip,ecsym,etop, erem,esid,
     &           nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3,
     &           jobname,iws,errorFlag)
c
c The routine generates an input deck for a FEM-analysis using WARP3D
c
      implicit none
      include 'plate_common_nod.f'
      include 'plate_common_eln.f'
c
      integer  analysis,keyhole,numLoadStep,tcr,tge,tbc,eqsolver,etyp,
     &         elnum,no_of_nodes,
     &         ecrsur(2,*),ecrtip(2,*),ecsym(2,*),etop(2,*),
     &         erem(2,*),esid(2,*),iws,errorFlag
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     &         na1,na2,na3,noda(na1,na2,na3),
     &         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      double precision  dmat(*),rgeo(*),rload(*),crfnvec(*)
      character jobname*200, file0*200,file1*200,file2*200,
     &                       file3*200,file4*200,file5*200
c
      integer    ndoftot,iodef, iotmp
      parameter (ndoftot = 100000, iodef=11, iotmp=21)
      integer    id(3,ndoftot),io, i,j,iblock,ipart,blocksize,
     &           nstresspoints,step,totstep
      double precision  rd(3,ndoftot),emod,xnu,sigy,r,hard
      character  gcase*24,form*9,elist*11,loadpatt*12,loadname*9
      logical    singlefile,ip
c
      integer      mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
c Define file names
c
      singlefile = .true.
c      singlefile = .false.
      call appendFile(jobname,'_wrp.inp',file0)
      call appendFile(jobname,'_wrp_crd.inp',file1)
      call appendFile(jobname,'_wrp_elm.inp',file2)
      call appendFile(jobname,'_wrp_dsp.inp',file3)
      call appendFile(jobname,'_wrp_trc.inp',file4)
      call appendFile(jobname,'_wrp_pre.inp',file5)
c
c  Set load arrays id(3,*), rd(3,i)to zero.
c
      do i=1, ndoftot
         id(1,i) = 0
         id(2,i) = 0
         id(3,i) = 0
         rd(1,i) = 0.
         rd(2,i) = 0.
         rd(3,i) = 0.
      enddo
c
c Start Generate the in-put deck to WARP3D
c
      io = iodef
      open( unit=io, file=file0, status='unknown' )
c
c  write out a general discription of the problem as a comment
c
      write(gcase,'(t1,3(a,i3))')' TGE=',tge,' TCR=',tcr,' TBC=',tbc
      write(io,101) 'c '
      write(io,101) 'c Geometry/Crack/BoundaryCond: ',gcase
      write(io,101) 'c '
      write(io,101) 'structure surface_crack'
      write(io,101) 'c'
101   format(t1,10a)
c
c Material definition
c
c  Analysis = 1;  Static Linear Elastic
c        dmat(1) = no of material constants defined = 2
c        dmat(2) = Young's modulii
c        dmat(3) = Poisson's ratio
      if (analysis .eq. 1) then
         emod = dmat(2)
         xnu  = dmat(3)
         sigy = 5.*emod
         write(io,101) 'material elastic'
         write(io,102) '  properties deformation  e ',emod,' nu ',xnu,
     &             ',','     yld_pt ',sigy,' n_power 1.001  rho 0.0'
102      format(t1,a,g15.8,a,f8.5,a/t1,a,g15.8,a)
c
c  Analysis = 2;  Deformation Plasticity (Non-Linear Elastic)
c        dmat(1) = no of material constants defined = 4 (5 if abaqus)
c        dmat(2) = Young's modulii
c        dmat(3) = Poisson's ratio
c        dmat(4) = sigY
c        dmat(5) = n = hardening
c        dmat(6) = alpha  = offset strain (Ramberg-Osgood model)
      elseif (analysis .eq. 2) then
         emod = dmat(2)
         xnu  = dmat(3)
         sigy = dmat(4)
         hard = dmat(5)
         write(io,101) 'material def_plast'
         write(io,103) '  properties deformation  e ',emod,' nu ',xnu,
     &         ',','     yld_pt ',sigy,' n_power ', hard,'  rho 0.0'
103     format(t1,a,g15.8,a,f8.5,a/t1,2(a,g15.8),a)
c
c  Analysis = 3;  Incremental Plasticity
c        dmat(1) = no of material constants defined >= 4
c        dmat(2) = Young's modulii
c        dmat(3) = Poisson's ratio
c        dmat(4) = sig1
c        dmat(5) = eps1
c        dmat(6) = sig2
c        dmat(7) = eps2
c        dmat(8) = ...
      elseif (analysis .eq. 3) then
         emod = dmat(2)
         xnu  = dmat(3)
         write(io,101) 'material inc_plast'
         write(io,104) '  properties deformation  e ',emod,' nu ',xnu,
     &             ',','     stress-strain curve 1  rho 0.0'
c
c   Define Stress-Strain Curve
         write(io,101) 'c'
         write(io,101) 'stress-strain curve 1'
         nstresspoints = int( dmat(1)-2 + 0.5 )
         do i=1, nstresspoints-2, 2
            write(io,105) dmat(i+4),dmat(i+3),','
         enddo
         write(io,105) dmat(nstresspoints+3),dmat(nstresspoints+2)
104      format(t1,a,g15.8,a,f8.5,a/t1,a)
105      format(t1,g15.8,g15.8,a)
      endif
c
c Define the size of the Model
c
      write(io,101) 'c'
      write(io,106) 'number of nodes ',no_of_nodes,'  elements ',elnum
106   format(t1,a,i6,a,i6)
c
c Coordinates
c
      write(io,101) 'c'
      write(io,101) '*echo off'
      write(io,101) 'coordinates'
      if (.not.singlefile) call open_2wrpfile(io,iotmp,file1)
      do i=1, inm
         if (nnr(i).gt.0) write(io,107) nnr(i),
     &                    npos(i,1),npos(i,2),npos(i,3)
      enddo
      if (.not.singlefile) call close_2wrpfile(io,iodef)
      write(io,101) 'c'
107   format(t1,i6,3(g20.10))
c
c Connectivity
c
      write(io,101) 'elements'
c
c  Analysis = 1;  Static Linear Elastic
      if (analysis .eq. 1) then
         r = elnum
         r = log10(r)
         i = 1 + int(r)
         write(form,'(t1,a,i1,a)') '(t1,a,i', i, ')'
         elist = '           '
         write(elist,form) '  1-',elnum
         write(io,108) elist,' type l3disop  linear ',
     &                       ' material elastic ,',
     &     '        ',' order 2x2x2    bbar    center_output  short'
c
c  Analysis = 2;  Deformation Plasticity (Non-Linear Elastic)
      elseif (analysis .eq. 2) then
         r = elnum
         r = log10(r)
         i = 1 + int(r)
         write(form,'(t1,a,i1,a)') '(t1,a,i', i, ')'
         elist = '           '
         write(elist,form) '  1-',elnum
         write(io,108) elist,' type l3disop  linear ',
     &                       ' material def_plast ,',
     &     '        ',' order 2x2x2    bbar    center_output  short'
c
c  Analysis = 3;  Incremental Plasticity
      elseif (analysis .eq. 3) then
         r = elnum
         r = log10(r)
         i = 1 + int(r)
         write(form,'(t1,a,i1,a)') '(t1,a,i', i, ')'
         elist = '           '
         write(elist,form) '  1-',elnum
         write(io,108) elist,' type l3disop  linear ',
     &                       ' material inc_plast ,',
     &     '        ',' order 2x2x2    bbar    center_output  short'
      endif
108   format(t1,a,a,a/t1,a,a)
c
      write(io,101) 'c'
      write(io,101) 'incidences'
      if (.not.singlefile) call open_2wrpfile(io,iotmp,file2)
      do i=1, elnum
         write(io,109) eln(i,0),( nnr(eln(i,j)), j=1, etyp)
      enddo
      if (.not.singlefile) call close_2wrpfile(io,iodef)
      write(io,101) 'c'
109   format(t1,9(i7) )
c
c Element Blocking
c
      write(io,101) 'blocking   $    sequential ordering'
      blocksize = 128
      iblock    = 0
      ipart     = mod(elnum,blocksize)
      do i=1, elnum-ipart, blocksize
         iblock = iblock + 1
         write(io,110) iblock, blocksize, i
      enddo
      if ( ipart .gt. 0 ) then
         iblock = iblock + 1
         write(io,110) iblock, ipart, elnum-ipart+1
      endif
110   format(t1,i10,i10,i10)
c
c Define displacement constraints (write out on file)
c
c Nodal Constraints and Equivalent Nodal forces are stored in rd(3,i)
c
      call wrp_fix_disp(analysis,tcr,tge,id,rd, nods,ns1,ns2,ns3,
     &         noda,na1,na2,na3, nodb,nb1,nb2,nb3,iws,errorFlag)
      if (tbc.eq.1) call wrp_pre_disp(tcr,tge,rgeo,rload, 
     &       id,rd, noda,na1,na2,na3, nodb,nb1,nb2,nb3, iws,errorFlag)
c
      write(io,101) 'constraints'
      write(io,101) '*echo off'
      if (.not.singlefile) call open_2wrpfile(io,iotmp,file3)
      call wrp_write_constraint(io,id,rd,no_of_nodes,iws) 
      if (.not.singlefile) call close_2wrpfile(io,iodef)
c
c Define Applied Load
c
c  TBC = 1  Displacements on Remote Boundary
c  TBC = 2  Tractions on Remote Boundary (includes Internal Pressure)
c  TBC = 3  Tractions on Crack Face
c
      if (tbc .eq. 1) then
         loadpatt = 'disp_pattern'
         loadname = 'pre_disp '
         write(io,101) 'c'
         write(io,101) 'loading ',loadpatt
         write(io,101) '  nodal loads '
         write(io,101) '    1 force_x 0.0'
         write(io,101) 'c'
      elseif (tbc.eq.2) then
         call wrp_pre_traction(tge,tcr,tbc,rgeo,rload,ip,id,rd,
     &            ecrsur,ecrtip,erem,esid,no_of_nodes,iws,errorFlag)
         loadpatt = 'trac_pattern'
         loadname = 'rem_trac '
         write(io,101) 'c'
         write(io,101) 'loading ',loadpatt
         write(io,101) '  nodal loads '
         if (.not.singlefile) call open_2wrpfile(io,iotmp,file4)
         call wrp_write_traction(io,id,rd,no_of_nodes,iws)
         if (.not.singlefile) call close_2wrpfile(io,iodef)
         if (ip) then
c    "element loads" also written in wrp_elem_pressure, skip it here
c            write(io,101) '  element loads '
            if (.not.singlefile) call open_2wrpfile(io,iotmp,file5)
            call wrp_elem_pressure(io,tge,rload,ecrsur,ecsym,etop)
         endif
         write(io,101) 'c'
      elseif (tbc.eq.3) then
         call wrp_pre_traction(tge,tcr,tbc,rgeo,rload,ip,id,rd,
     &            ecrsur,ecrtip,erem,esid,no_of_nodes,iws,errorFlag)
         loadpatt = 'trac_pattern'
         loadname = 'face_trac'
         write(io,101) 'c'
         write(io,101) 'loading ',loadpatt
         write(io,101) '  nodal loads '
         if (.not.singlefile) call open_2wrpfile(io,iotmp,file4)
         call wrp_write_traction(io,id,rd,no_of_nodes,iws)
         if (.not.singlefile) call close_2wrpfile(io,iodef)
         write(io,101) 'c'
      endif
      write(io,101) '*echo on'
c
c Define load steps
c
c  Define the number of load steps here to set the loading
c  history and loading increment.
c  One load step for elastic analysis, use numLoadStep for
c  elastic-plastic multi step analysis.
c  (up to 99 load steps for i2 format in line 151)
c
      if (analysis .eq.1) then
         totstep = 1
      elseif (analysis.eq.2) then
         totstep = numLoadStep
      elseif (analysis.eq.3) then
         totstep = numLoadStep
      endif
c check for a bad number of load steps and reset
      if (totstep.le.0) then
         totstep = 1
      end if
c
      write(io,101) 'c'
      write(io,101) 'c     step cases:'
      write(io,101) 'c'
      write(io,101) 'loading ',loadname
      write(io,101) '  nonlinear'
      write(io,101) '    steps  1-100 ',loadpatt,'  1.0'
c
c  A single load step, apply the full given load pattern,
c  use a 1.0 multiplier
cc      if (totstep.eq.1) then
cc         write(io,121) totstep,loadpatt
cc121      format(1x,'    steps ',i1,1x,a12,'  1.0')
c  For several load steps use an incremental load pattern multiplier,
c  use 1.0/totstep as a uniform incremental load multiplier,
c  the incremental load is additive each load step so that the full
c  load pattern is applied at the final load step.
cc      else
cc         write(io,122) totstep,loadpatt,(1.0/totstep)
cc122      format(1x,'    steps 1-',i2,1x,a12,1x,f6.4)
cc      end if
c
c Define Solution parameters
c
c   EQSOLVER = 1   Direct Sparse Equation Solver
c   EQSOLVER = 2   Conjugate Gradient Solver (diagonal precond.)
c
      write(io,101) 'c'
      write(io,101) 'c    solution paramters'
      write(io,101) 'c'
      write(io,101) 'nonlinear analysis parameters'
      if (eqsolver.eq.1) then
         write(io,101) '  solution technique direct sparse'
      else
         write(io,101) '  solution technique lnpcg'
         write(io,101) '  preconditioner type diagonal'
c
c    ?????? preconditioner type Hughes-Winget = ebe
c         write(io,101) '  preconditioner type Hughes-Winget'
         write(io,101) '  lnpcg conv test res tol 0.01'
         write(io,101) '  maximum linear iterations 5000'
      endif
      write(io,101) 'c  solution technique direct sparse'
      write(io,101) 'c  solution technique direct sparse hp'
      write(io,101) 'c  solution technique direct sparse sgi'
      write(io,101) 'c  solution technique direct sparse bcs'
      write(io,101) 'c  solution technique lnpcg'
      write(io,101) 'c  preconditioner type diagonal'
      write(io,101) 'c  preconditioner type ebe'
      write(io,101) 'c  lnpcg conv test res tol 0.01'
      write(io,101) 'c  maximum linear iterations 5000'
c
c Newton-Raphson parameters (necessary even in the case of lin. elast.)
c
      write(io,101) '  maximum iterations 5  $ default is 10'
      write(io,101) '  minimum iterations 1'
      write(io,101) '  convergence test norm res tol 0.001'

      write(io,101) '  nonconvergent solutions stop'
      write(io,101) '  adaptive solution on'
      write(io,101) '  linear stiffness for iteration one off'
      write(io,101) '  batch messages off'
      write(io,101) '  cpu time limit off'
      write(io,101) '  material messages off'
      write(io,101) '  bbar stabilization factor 0.0'
      write(io,101) '  consistent q-matrix on'
      write(io,101) '  time step 100000'
      write(io,101) '  trace solution on lpcg_solution off'
      write(io,101) '  extrapolate on'
c
c Define Compute steps
c
      write(io,101) 'c'
      write(io,101) 'c    start the analysis'
      write(io,101) 'c'
c
      do step=1, totstep
         write(io,111) 'compute displacements for load ',
     &                  loadname,' steps ',step
         if (no_of_nodes.lt.1000) then
           write(io,131) 'output eformat displ for nodes 1-',no_of_nodes
           write(io,131) 'output eformat react for nodes 1-',no_of_nodes
         elseif (no_of_nodes.lt.10000) then
           write(io,132) 'output eformat displ for nodes 1-',no_of_nodes
           write(io,132) 'output eformat react for nodes 1-',no_of_nodes
         elseif (no_of_nodes.lt.100000) then
           write(io,133) 'output eformat displ for nodes 1-',no_of_nodes
           write(io,133) 'output eformat react for nodes 1-',no_of_nodes
         else
           write(io,134) 'output eformat displ for nodes 1-',no_of_nodes
           write(io,134) 'output eformat react for nodes 1-',no_of_nodes
         endif
 131     format(t1,a,i3)
 132     format(t1,a,i4)
 133     format(t1,a,i5)
 134     format(t1,a,i6)
c
c Define Domain Integral Evaluation
c
         write(io,101) 'c'
         write(io,101) 'c     J-integral evaluations'
         write(io,101) 'c'
         if (keyhole.eq.0) then
            call wrp_write_domain1(io,mr,m1,m2,crfnvec,ns1,ns2,ns3,nods)
         else
            call wrp_write_domain2(io,mr,m1,m2,crfnvec,ns1,ns2,ns3,nods)
         endif
         write(io,101) 'c'
      enddo
 111  format(t1,a,a,a,i3)
      write(io,101) 'stop'
c
      close(io)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine open_2wrpfile(io,iotmp,file2nd)
      integer   io,iotmp
      character file2nd*17
      write(io,'(t1,a,a,a)') '*input from ''',file2nd,''''
      io = iotmp
      open( unit=iotmp, file=file2nd, status='unknown' )
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine close_2wrpfile(io,iodef)
      integer   io,iodef
      close(io)
      io = iodef
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine warp3d_post(tge,tcr,tbc,dmat,nods,ns1,ns2,ns3,jobname)
c
c The routine generates an info file for
c    Post-Processing of WARP3D results
c
      implicit none
      integer    tge,tcr,tbc,ns1,ns2,ns3,nods(ns1,ns2,ns3),
     &           io,ntip(200),nt,i,k
      parameter (io=24)
      double precision dmat(*),xtip(3,200)
      character        jobname*200, infofile*200,space*80
c
c  Retrieve tip coordinates.
c
      call ncrd_ntip(ntip,xtip,nt, nods,ns1,ns2,ns3)
c
c Open file = "plate_crk.inf"  and
c write out the coordinates of the crack front.
c
      call appendFile( jobname, '_crf.inf', infofile )
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
      subroutine wrp_write_constraint(io,id,rd,no_of_nodes,iws) 
c
c Write out constraint as defined in the global arrays id(3,i), rd(3,i)
c
      implicit none
      integer          io,id(3,*),no_of_nodes,iws,n,ix,iy,iz
      double precision rd(3,*), u,v,w
c
      do n=1, no_of_nodes
         ix = 0
         iy = 0
         iz = 0
         u  = rd(1,n)
         v  = rd(2,n)
         w  = rd(3,n)
         if ( id(1,n) .eq. 1 ) ix = 1
         if ( id(2,n) .eq. 1 ) iy = 1
         if ( id(3,n) .eq. 1 ) iz = 1
ctest         write(iws,'(t1,i6,3i2,3g15.6)') node,ix,iy,iz,u,v,w 
c
         if     ( (ix.eq.1) .and. (iy.ne.1) .and. (iz.ne.1) ) then
            write(io,101) n,'  u ',u
         elseif ( (ix.ne.1) .and. (iy.eq.1) .and. (iz.ne.1) ) then
            write(io,101) n,'  v ',v
         elseif ( (ix.ne.1) .and. (iy.ne.1) .and. (iz.eq.1) ) then
            write(io,101) n,'  w ',w
         elseif ( (ix.eq.1) .and. (iy.eq.1) .and. (iz.ne.1) ) then
            write(io,102) n,'  u ',u,' v ',v
         elseif ( (ix.eq.1) .and. (iy.ne.1) .and. (iz.eq.1) ) then
            write(io,102) n,'  u ',u,' w ',w
         elseif ( (ix.ne.1) .and. (iy.eq.1) .and. (iz.eq.1) ) then
            write(io,102) n,'  v ',v,' w ',w
         elseif ( (ix.eq.1) .and. (iy.eq.1) .and. (iz.eq.1) ) then
            write(io,103) n,'  u ',u,' v ',v,' w ',w
         endif
      enddo
c
 101  format(i6,a,g16.8)
 102  format(i6,2(a,g16.8))
 103  format(i6,3(a,g16.8))
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_write_traction(io,id,rd,no_of_nodes,iws)
c
c  Write out the equivalent nodal loads due to prescribed tractions
c
      implicit none
      integer          io,id(3,*),no_of_nodes,iws,ix,iy,iz,n
      double precision rd(3,*),fx,fy,fz
      character        cx*8,cy*8,cz*8
      data             cx,cy,cz /' force_x',' force_y',' force_z'/
c
c
      do n=1, no_of_nodes
         ix = 0
         iy = 0
         iz = 0
         fx = rd(1,n)
         fy = rd(2,n)
         fz = rd(3,n)
         if ((id(1,n).eq.2).or.(id(1,n).eq.3)) ix = 1
         if ((id(2,n).eq.2).or.(id(2,n).eq.3)) iy = 1
         if ((id(3,n).eq.2).or.(id(3,n).eq.3)) iz = 1
c         write(iws,201) n, id(1,n),id(2,n),id(3,n), fx,fy,fz 
c 201     format(t1,i6,3i2,3g15.6)
c
         if     ((ix.eq.1).and.(iy.ne.1).and.(iz.ne.1)) then
            write(io,102) n,cx,fx
         elseif ((ix.ne.1).and.(iy.eq.1).and.(iz.ne.1)) then
            write(io,102) n,cy,fy
         elseif ((ix.ne.1).and.(iy.ne.1).and.(iz.eq.1)) then
            write(io,102) n,cz,fz
         elseif ((ix.eq.1).and.(iy.eq.1).and.(iz.ne.1)) then
            write(io,103) n,cx,fx,cy,fy
         elseif ((ix.eq.1).and.(iy.ne.1).and.(iz.eq.1)) then
            write(io,103) n,cx,fx,cz,fz
         elseif ((ix.ne.1).and.(iy.eq.1).and.(iz.eq.1)) then
            write(io,103) n,cy,fy,cz,fz
         elseif ((ix.eq.1).and.(iy.eq.1).and.(iz.eq.1)) then
            write(io,104) n,cx,fx,cy,fy,cz,fz
        endif
      enddo
      write(io,101)'c'
 101  format(t1,a)
 102  format(t1,i6,a,g14.6)
 103  format(t1,i6,2(a,g14.6))
 104  format(t1,i6,3(a,g14.6))
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_elem_pressure(io,tge,rload,ecrsur,ecsym,etop)
c
c Define Uniform pressure on element faces
c
      implicit none
      include 'plate_common_nod.f'
c
      integer  io,tge,ecrsur(2,*),ecsym(2,*),etop(2,*)
      double precision rload(*),press
c
      press = rload(7)
      if     ( (tge.eq.201).or.(tge.eq.203) ) then
         call wrp_write_pressure(io,ecrsur,press)
         call wrp_write_pressure(io,ecsym,press)
      elseif ( (tge.eq.202).or.(tge.eq.204) ) then
         call wrp_write_pressure(io,etop,press)
c
      elseif ( (tge.eq.301).or.(tge.eq.303).or.
     &         (tge.eq.305).or.(tge.eq.307) ) then
         call wrp_write_pressure(io,ecrsur,press)
         call wrp_write_pressure(io,ecsym,press)
      elseif ( (tge.eq.302).or.(tge.eq.304).or.
     &         (tge.eq.306).or.(tge.eq.308) ) then
         call wrp_write_pressure(io,etop,press)
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_write_pressure(io,eset,p)
      implicit none
      integer          io,eset(2,*),ne,i,ie,face_no
      double precision p
      character*7 cface(6)
      data (cface(i),i=1,6) / ' face 1', ' face 2', ' face 3',
     &                        ' face 4', ' face 5', ' face 6' /
c
c  The default element faces are defined according to the WARP3D scheme 
c
      ne = eset(1,1)
c
      write(io,101)'c'
      write(io,101)'c Define Uniform Pressure on Element Faces'
      write(io,101)'c'
      write(io,101)'   element loads'
      do i=1, ne
         ie      = eset(1,i+1)
         face_no = eset(2,i+1)
         write(io,102) '      ',ie,cface(face_no),' pressure ',p
      enddo
      write(io,101)'c'
 101  format(t1,a)
 102  format(t1,a,i6,1x,a,a,g17.10)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_write_domain1(io,mr,m1,m2,crfnvec,ns1,ns2,ns3,nods)
c
c Writes out definitions for domain integral evaluations
c
c Current implementation (July -98) l3disop element only
c
      implicit none
      integer  io,mr,m1,m2,ns1,ns2,ns3,nods(ns1,ns2,ns3)
      integer          ntip(200),nt,nc,i
      double precision crfnvec(*),xtip(3,200),s(200),ds
      character*6      dname(200)
      character*18     dnode(200)
      character*6      dtype(200)
      character*4      dring
c
c  Retrieve tip coordinates and determine maximum # of domains
c
      call ncrd_ntip(ntip,xtip,nt, nods,ns1,ns2,ns3)
c      nc = min(m1-3,m2-2)
      nc = min(m1-2,m2-1)
      if (nc.lt.0) nc = 0
      nc = nc + mr
      if (nc.lt.9) then
         write(dring,'(t1,a,i1)') ' 3-',nc 
      else
         write(dring,'(t1,a,i2)') ' 3-',nc 
      endif
c
c  Compute a natural coordinate along the Crack Front, define the 
c  Domain names, define the front nodes and function type.
c
      dname(1) = 'dnr001'
      write(dnode(1),'(t1,i6,i6,a6)') ntip(1),ntip(2),'      '
      dtype(1) = 'type a'
      s(1) = 0.
      do i=2, nt
         ds = sqrt( (xtip(1,i)-xtip(1,i-1))**2.
     &            + (xtip(2,i)-xtip(2,i-1))**2.
     &            + (xtip(3,i)-xtip(3,i-1))**2. )
         s(i) = s(i-1) + ds
         if (i.lt.10) then
            write(dname(i),'(t1,a,i1)') 'dnr00',i
         elseif (i.lt.100) then
            write(dname(i),'(t1,a,i2)') 'dnr0',i
         else
            write(dname(i),'(t1,a,i3)') 'dnr',i
         endif
         write(dnode(i),'(t1,i6,i6,i6)') ntip(i-1),ntip(i),ntip(i+1)
         dtype(i) = 'type b'
      enddo
      write(dnode(nt),'(t1,i6,i6,a6)') ntip(nt-1),ntip(nt),'      '
      dtype(nt) = 'type c'
c
c Write out the domains on unit=io
c
      do i=1, nt
         write(io,101) 'domain ',dname(i),
     &                     '  $ Crack Front Pos. s =',s(i)/s(nt)
         write(io,102) '  normal plane  nx',crfnvec(1),
     &                                ' ny',crfnvec(2),' nz',crfnvec(3)
         write(io,103) '  symmetric'
         write(io,103) '  front nodes',dnode(i),' linear  verify'
         write(io,103) '  q-values automatic rings ',dring
         write(io,103) '  function ',dtype(i)
         write(io,103) '  print totals'
c         write(io,103) '  use 1 point rule'
         write(io,103) 'compute domain integral'
         write(io,'(t1,a)') 'c'
      enddo
c
 101  format(t1,a,a,a,f7.3) 
 102  format(t1,3(a,f10.6))
 103  format(t1,a,a,a,a)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine wrp_write_domain2(io,mr,m1,m2,crfnvec,
     &                                ns1,ns2,ns3,nods)
c
c Writes out definitions for domain integral evaluations
c
c Current implementation (July -98) l3disop element only
c
      implicit none
      integer  io,mr,m1,m2,ns1,ns2,ns3,nods(ns1,ns2,ns3)
      integer          ntip(200),nt,ns(2,200),nc,i,i1,i2
      double precision crfnvec(*),xtip(3,200),s(200),ds
      character*6      dname(200)
      character*9     dnode(200)
      character*6      dtype(200)
      character*4      dring
c
c  Retrieve tip coordinates and determine maximum # of domains
c
      call ncrd_ntipset(ntip,xtip,nt,ns, nods,ns1,ns2,ns3)
c
c      nc = min(m1-3,m2-2)
      nc = min(m1-2,m2-1)
      if (nc.lt.0) nc = 0
      nc = nc + mr
      if (nc.lt.9) then
         write(dring,'(t1,a,i1)') ' 3-',nc 
      else
         write(dring,'(t1,a,i2)') ' 3-',nc 
      endif
c
c  Compute a natural coordinate along the Crack Front, define the 
c  Domain names, define the front nodes and function type.
c
      dname(1) = 'dnr001'
      write(dnode(1),'(t1,i3,i3,a3)') 1,2,'   '
      dtype(1) = 'type a'
      s(1) = 0.
      do i=2, nt
         i1 = ns(1,i-1)
         i2 = ns(1,i)
         ds = sqrt( (xtip(1,i2)-xtip(1,i1))**2.
     &            + (xtip(2,i2)-xtip(2,i1))**2.
     &            + (xtip(3,i2)-xtip(3,i1))**2. )
         s(i) = s(i-1) + ds
         if (i.lt.10) then
            write(dname(i),'(t1,a,i1)') 'dnr00',i
         elseif (i.lt.100) then
            write(dname(i),'(t1,a,i2)') 'dnr0',i
         else
            write(dname(i),'(t1,a,i3)') 'dnr',i
         endif
         write(dnode(i),'(t1,i3,i3,i3)') i-1,i,i+1
         dtype(i) = 'type b'
      enddo
      write(dnode(nt),'(t1,i3,i3,a3)') nt-1,nt,'   '
      dtype(nt) = 'type c'
c
c Write out the domains on unit=io
c
      do i=1, nt
         write(io,101) 'domain ',dname(i),
     &                     '  $ Crack Front Pos. s =',s(i)/s(nt)
         write(io,102) '  normal plane  nx',crfnvec(1),
     &                                ' ny',crfnvec(2),' nz',crfnvec(3)
         write(io,103) '  symmetric'
c
c Write out node sets
c
         if (i.eq.1) then
            call write_wrp_tipsets(io,i,ns,ntip)
            call write_wrp_tipsets(io,i+1,ns,ntip)
         elseif (i.eq.nt) then
            call write_wrp_tipsets(io,i-1,ns,ntip)
            call write_wrp_tipsets(io,i,ns,ntip)
         else
            call write_wrp_tipsets(io,i-1,ns,ntip)
            call write_wrp_tipsets(io,i,ns,ntip)
            call write_wrp_tipsets(io,i+1,ns,ntip)
         endif
         write(io,103) '  front node sets',dnode(i),' linear'
         write(io,103) '  q-values automatic rings ',dring
         write(io,103) '  function ',dtype(i)
         write(io,103) '  print totals'
c         write(io,103) '  use 1 point rule'
         write(io,103) 'compute domain integral'
         write(io,'(t1,a)') 'c'
      enddo
c
 101  format(t1,a,a,a,f7.3) 
 102  format(t1,3(a,f10.6))
 103  format(t1,a,a,a,a)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine write_wrp_tipsets(io,i,ns,ntip)
      implicit none
      integer  io,i,ns(2,*),ntip(*), j,n,iv(100)
      n = 0
      do j=ns(1,i), ns(2,i)
         n = n + 1
         iv(n) = ntip(j)
      enddo
      if (n.le.11) then
         write(io,101) '  node set ',i,(iv(j),j=1,n)
      else
         write(io,102) '  node set ',i,(iv(j),j=1,11),',',
     &                '               ',(iv(j),j=12,n)
      endif
 101  format(t1,a,i2,'  ',11i5)
 102  format(t1,a,i2,'  ',11i5,a/a,11i5)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
