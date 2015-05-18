c----67--1---------2---------3---------4---------5---------6---------712
c
c    File: plate1.f   Revised April 23, 1998
c    plate1.f is the primary subroutine for the FEACode mesh generation
c
c    "write(*,...)" changed to "write(iws,...)", 7/29/97
c
c----67--1---------2---------3---------4---------5---------6---------712
c
ccc
c  Convert the main program "plate3d" to the main subroutine called
c  from the shell/DLL accessible subroutine.
c  (GVT 7/29/97, note: look for the "ccc" triple comment for new
c                      additions )
CJF      program plate3d
      subroutine plate1(errorFlag,filename,jobname)
c
c  filename = input file name ("file.in") for the input data file
c             giving the geometric parameters for the mesh generation,
c             use a long character string to allow for a lengthy file
c             name, filename should have no leading blanks
c
c  errorFlag = integer error flag to indicate if the program has had
c              an error and the user should check the *.sta status
c              file in the same directory as the input file
c              (4 byte integer)
c           = 0 = no error
c           = 1 = error, check the *.sta status file
c
c
c----67--1---------2---------3---------4---------5---------6---------712
c
c  This program generates a 3-D FE-model of a plate containing a
c  surface crack of elliptical shape. The crack front is embeddied
c  in a tube following the ellips. A focused mesh is used around
c  the crack front within the tube. In-put deck are generated for
c  the FE-programs WARP3D and ABAQUS.
c
c          VERSION July 11 1997   / Jonas Faleskog /
c
c          Modified  Nov. 22 1999 /JF/
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      implicit none
c
c Array bounds   Zone S: nods(i,j,k) nsm = ns1*ns2*ns3
c                Zone A: noda(i,j,k) nam = na1*na2*na3
c                Zone B: nodb(i,j,k) nbm = nb1*nb2*nb3
c
      include 'plate_common_eln.f'
      include 'plate_common_nod.f'
ccc  Declare the input arguments
      integer errorFlag
      character filename*200,jobname*200
ccc
      integer    nsm,nam,nbm, nesm,nesm2
       parameter (nsm=60000,nam=900000,nbm=1700000)
ccc      parameter (nsm=20000,nam=300000,nbm=600000)
ccc      parameter (nsm=20000,nam=140000,nbm=140000)
ccc      parameter (nsm=20000,nam=80000,nbm=80000)
      parameter (nesm=200, nesm2=nesm*nesm)
c
      integer  iws, ns1,ns2,ns3, na1,na2,na3, nb1,nb2,nb3
c
      integer  nea2,neb2, estk_s(2,nesm), estk_a(2,nesm2),
     &         estk_b(2,nesm2),nstk_s(3,401),nstk_a(3,401),
     &         nstk_b(3,401), tcr,tge,tbc,keyhole,analysis,
     &         numLoadStep,eqsolver,ncq, i
c
      integer  nods(nsm),noda(nam),nodb(nbm),
     1         etyp,elnum,no_of_nodes
c
      integer  eset_size
      parameter (eset_size = 20000)
      integer  efront(2,eset_size),ecrsur(2,eset_size),
     &         ecrtip(2,eset_size),easym(2,eset_size),
     &         ecsym(2,eset_size),erem(2,eset_size),
     &         esid(2,eset_size),etop(2,eset_size)
c
      double precision dmat(2000), rgeo(20),rload(10),
     1                 yg(0:200),xag(0:200),xbg(0:200),zg(2,0:200),
     2                 rfm(50),bias(10),crfnvec(3),cqvec(3,200)
c Make the status filename "stfile" the same size as the input filename
      character    prog*20,NeutralFormat*20,stfile*200
      logical      demoFlag, statusFlag
ccc
c
      integer      mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
      integer         kstart,kstep
      common /nblock/ kstart,kstep
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
      double precision t,w,c,a,alfa
      common /geom/    t,w,c,a,alfa
ccc
      errorFlag = 0
c Check for leading blank spaces in the filename string and remove
cgt      filename = ADJUSTL(filename)
ccc
c
c
c..Input is read from file mesh3d_scp.in
c
      call read_input( prog,etyp,analysis,numLoadStep,eqsolver,
     &          tcr,tge,tbc,rload, dmat,keyhole,
     &          nods,ns1,ns2,ns3,nsm, noda,na1,na2,na3,nam,
     &          nodb,nb1,nb2,nb3,nbm, nesm,nesm2,nea2,neb2,
     &          estk_s,estk_a,estk_b, rgeo,yg,xag,xbg,zg,rfm,bias,
     &          NeutralFormat, demoFlag,statusFlag, iws,errorFlag,
     &          filename,jobname,stfile)
      if ( errorFlag .gt. 0 ) goto 999
c
c TEST
c
c Note!  "write(iws,...)" to write to the status *.sta file
c
      if (demoFlag) write(iws,'(/t1,a/)') '>>> Demo Version <<<'
      write(iws,'(t5,a,i2,a,i4,a,i2)')  'TCR =',tcr,
     &                    '   TGE =',tge,'   TBC =',tbc
      write(iws,'(t5,3(a,g12.4))') 'W0 =',w,'  L0 =',yg(lt),'  t0 =',t
      write(iws,'(t5,2(a,g12.4)/)') 'Rc =',rgeo(4),'  Rp =',rgeo(5)
      write(iws,*) ' dmat:'
      do i=1, (int(dmat(1)+0.5) + 1)
         write(iws,*) real(dmat(i))
      enddo
c
c      write(iws,*) ' After read_input'
c      write(iws,*) ' indicies'
c      write(iws,'(t1,a,8i4)') ' mf,mr: ',mf,mr
c      write(iws,'(t1,a,8i4)') ' mv,mh,m1,m2,ma,na: ',mv,mh,m1,m2,ma,na
c      write(iws,'(t1,a,8i4)') ' mb,nb: ',mb,nb
c      write(iws,'(t1,a,8i4)') ' lt,lred: ',lt,lred
c      write(iws,'(t1,a,8i4)') ' ims,jms,kms: ',ims,jms,kms
c      write(iws,'(t1,a,8i4)') ' ima,jma,kma: ',ima,jma,kma
c      write(iws,'(t1,a,8i4)') ' imb,jmb: ',imb,jmb
c      write(iws,'(t1,a,8i4)') ' ksr1,kar1,kar2: ',ksr1,kar1,kar2
c      write(iws,'(t1,a,8i4)') ' rtype,sfred,sfred_type,sjred_type: ',
c     1                        rtype,sfred,sfred_type,sjred_type
c   write(iws,'(t1,a,5g12.4)') 't,w:',t,w
c   write(iws,'(t1,a,5g12.4)') 'c,a,alfa:',c,a,alfa
c      write(iws,'(t3,a)') 'BIAS'
c      do i=1, 8
c         write(iws,'(g15.6)') bias(i)
c      enddo
c      write(iws,'(t3,a)') 'rfm'
c      do i=1, mr
c         write(iws,'(g15.6)') rfm(i)
c      enddo
c      write(iws,'(t3,a)') 'yg, zg'
c      do i=0, lt
c         write(iws,'(3g15.6)') yg(i),zg(1,i), zg(2,i)
c      enddo
c      write(iws,'(t3,a,i2)') 'xag;  m2 =',m2
c      do i=0, m2
c         if (i.le.1) then
c            write(iws,'(g15.6)') xag(i)
c         else
c            write(iws,'(2g15.6)') xag(i),
c     &                          (xag(i)-xag(i-1))/(xag(i-1)-xag(i-2))
c         endif
c      enddo
c      write(iws,'(t3,a,i2)') 'xbg;  mb =',mb
c      do i=0, mb
c         if (i.le.1) then
c            write(iws,'(g15.6)') xbg(i)
c         else
c            write(iws,'(2g15.6)') xbg(i),
c     &                          (xbg(i)-xbg(i-1))/(xbg(i-1)-xbg(i-2))
c         endif
c      enddo
c
c Generate node numbers, store in the arrays NODS,NODA & NODB
c
      call nodnumber(nods,ns1,ns2,ns3,noda,na1,na2,na3,
     &               nodb,nb1,nb2,nb3, iws,errorFlag)
      if ( errorFlag .gt. 0 ) then
         write(iws,*)'Stopped after error in nodnumber subroutine'
         goto 999
      endif
c
c Generate node coordinates
c
      if (statusFlag) then
c         close(unit=iws)
c         open(unit=iws, file=stfile, position='append')
      endif
      write(iws,'(/t1,a)') '* generating coordinates'
      write(*,'(/t1,a)') '* generating coordinates'
      call node_coordinates(iws,nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                      nodb,nb1,nb2,nb3, yg,xag,xbg,rfm, tcr,
     &                      keyhole,errorFlag)                     
      if ( errorFlag .gt. 0 ) then
        write(iws,*)'Stopped after error in node_coordinates subroutine'
         goto 999
      endif
c
c..Generate elements
c
      if (statusFlag) then
         close(unit=iws)
         open(unit=iws, file=stfile, position='append')
      endif
      write(iws,'(/t1,a)') '* generating elements'
      write(*,'(/t1,a)') '* generating elements'
      call make_element(nods,ns1,ns2,ns3, noda,na1,na2,na3, 
     1                  nodb,nb1,nb2,nb3, elnum, analysis,keyhole,
     2                  nea2,neb2, estk_s,estk_a,estk_b, prog,
     3                  iws,errorFlag)
      if ( errorFlag .gt. 0 ) then
         write(iws,*)'Stopped after error in make_element subroutine'
         goto 999
      end if
      write(iws,'(t15,a,i5,a)') '=> # of elements = ',elnum
      write(*,'(t15,a,i5,a)') '=> # of elements = ',elnum
c
c..Sort out the nodes in the model
c
      if (statusFlag) then
         close(unit=iws)
         open(unit=iws, file=stfile, position='append')
      endif
      write(iws,'(/t1,a)') '* sort out the unique nodes'
      write(*,'(/t1,a)') '* sort out the unique nodes'
      call sort_out_nodes(elnum,no_of_nodes,etyp,prog, 
     1                    nstk_s,nstk_a,nstk_b, nods,ns1,ns2,ns3,
     2                    noda,na1,na2,na3, nodb,nb1,nb2,nb3, iws)
c
c Coordinate transformation - type defined by parameter TGE
c
      if (statusFlag) then
         close(unit=iws)
         open(unit=iws, file=stfile, position='append')
      endif
      write(iws,'(/t1,a)') '* Coordinate Transformation'
      write(*,'(/t1,a)') '* Coordinate Transformation'
      call coord_transf(tge,elnum,etyp,rgeo,crfnvec,cqvec,ncq,
     &                  ns1,ns2,ns3,nods)
c
c Generate element sets containing element on all the surfacres
c of the model.
c
      if (statusFlag) then
         close(unit=iws)
         open(unit=iws, file=stfile, position='append')
      endif
      write(iws,'(/t1,a)') '* Retrieve element surfaces'
      write(*,'(/t1,a)')'* Retrieve element surfaces - may take a while'
      call gen_surface(efront,ecrsur,ecrtip,easym,ecsym,erem,esid,etop,
     1                 nods,ns1,ns2,ns3, noda,na1,na2,na3,
     2                 nodb,nb1,nb2,nb3, elnum,iws,errorFlag)
c
c Generate plot files containing wire frames of key surfaces
c
ccc      if (statusFlag) then
ccc         close(unit=iws)
ccc         open(unit=iws, file=stfile, position='append')
ccc      endif
ccc      write(iws,'(/t1,a)') '* Generate Graphics files'
ccc      write(*,'(/t1,a)') '* Generate Graphics files'
ccc      call gen_graph(efront,ecrsur,easym,ecsym,erem,esid,etop,
ccc     &              tge,no_of_nodes,elnum,iws,jobname)
c
c Generate FEA in-put deck and an info-file for Post-Processing
c
      if ( demoFlag ) goto 999
      if (statusFlag) then
         close(unit=iws)
         open(unit=iws, file=stfile, position='append')
      endif
      if (NeutralFormat.eq.'PATRAN') then
         call patran_neu(etyp,elnum,no_of_nodes,dmat,rgeo,rload,
     &                   analysis,tcr,tge,tbc,ecrsur,ecrtip,ecsym,etop,
     &                   erem,esid, nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                   nodb,nb1,nb2,nb3, jobname,iws,errorFlag)
      elseif (NeutralFormat.eq.'FEMAP') then
         call femap_neu(etyp,elnum,dmat,jobname)
      elseif (prog.eq.'ABAQUS') then
         write(iws,'(/t1,a,a)') '* Generate in-put deck for ',prog
         write(*,'(/t1,a)') '* Generate ABAQUS in-put deck'
         call abaqus_inp(etyp,elnum,analysis,numLoadStep,keyhole,
     &                   tge,tcr,tbc,rgeo,rload,
     &                   dmat,ecrsur,ecsym,erem,esid,etop,
     &                   crfnvec,cqvec,ncq,
     &                   nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                   nodb,nb1,nb2,nb3, jobname,iws,errorFlag)
         call abaqus_post(tge,tcr,tbc,dmat,nods,ns1,ns2,ns3,jobname)
c
      elseif (prog .eq. 'WARP3D') then
         write(*,'(/t1,a)') '* Generate WARP3D in-put deck'
         call warp3d_inp(etyp,elnum,no_of_nodes,dmat,rgeo,rload,
     &                   analysis,keyhole,numLoadStep,
     &                   tcr,tge,tbc,eqsolver,crfnvec,
     &                   ecrsur,ecrtip,ecsym,etop, erem,esid,
     &                   nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                   nodb,nb1,nb2,nb3,
     &                   jobname,iws,errorFlag)
         call warp3d_post(tge,tcr,tbc,dmat,nods,ns1,ns2,ns3,jobname)
c
      endif
c
      if (statusFlag) then
         close(unit=iws)
         open(unit=iws, file=stfile, position='append')
      endif
c
      write(iws,'(/t1,a)') '*  Done with mesh generation'
c
999   continue
      close(iws)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine read_input( prog,etyp,analysis,numLoadStep,eqsolver,
     &           tcr,tge,tbc,rload,dmat,keyhole,
     &           nods,ns1,ns2,ns3,nsm, noda,na1,na2,na3,nam,
     &           nodb,nb1,nb2,nb3,nbm, nesm,nesm2,nea2,neb2,
     &           estk_s,estk_a,estk_b, rgeo,yg,xag,xbg,zg,rfm,bias,
     &           NeutralFormat, demoFlag,statusFlag, iws,errorFlag,
     &           filename,jobname,stfile)
c
c The subroutine reads in-put parameters from the file "plate3d.in"
c mesh3d_nsct.in. Some of the parameters are checked so that
c they are located in a reasonable interval.
c
      implicit none
      include 'plate_common_nod.f'
c
      integer ns1,ns2,ns3,nsm,nods(nsm), na1,na2,na3,nam,noda(nam),
     1        nb1,nb2,nb3,nbm,nodb(nbm)
      integer analysis,numLoadStep,eqsolver,tcr,tge,tbc,keyhole
      integer iws,nesm,nesm2,nes2,nea2,neb2,estk_s(2,nesm),
     1        estk_a(2,nesm2),estk_b(2,nesm2)
      double precision  rload(*),dmat(*),rgeo(*),yg(0:200),xag(0:200),
     1                  xbg(0:200),zg(2,0:200),rfm(50),bias(10),b
      integer errorFlag
      character filename*200,jobname*200,stfile*200
c
      integer i,j,k,n,num,etyp,max_node_no,io,ivek(20),
     1        m1min,m2min,m12m,m12e,m12,mb1,nmin
      double precision  length,rvek(15),dx0,xmax,dxmax,
     3                  x(200),xi,beta,beta0,betam,c1,c2,psi,tmp,
     4                  lambda,lambda0,lambdam,lambda1
      double precision  bias_crack
ccc  Make the inpfile string the same size as filename
      character prog*20,NeutralFormat*20,row*200,inpfile*200
ccc  Old string sizes
ccc      character prog*20,row*200,inpfile*10
      logical   ok, demoFlag,statusFlag
      double precision  pi,one
      parameter (one=1.0)
c
      integer      mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer      ims,jms,kms,ima,jma,kma,imb,jmb
      common /max/ ims,jms,kms,ima,jma,kma,imb,jmb
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      double precision t,w,c,a,alfa
      common    /geom/ t,w,c,a,alfa
c
      pi = 4.*atan(one)
c
ccc
c  Use the filename given by the calling routine
c      inpfile = 'plate3d.in'
      inpfile = filename
      inquire( file=inpfile, exist=ok )
      iws = 44
c      open( unit=iws, file='tmp.sta',  status='scratch' )
      open( unit=iws, status='scratch' )
      if (.not.ok) then
         write(iws,'(t1,3a)')'>> The File: ',inpfile(1:20),
     &                       ' does not exist!'
         errorFlag = 1
         return
      endif
      write(iws,201) '* reading in-put-data from infile:',inpfile
      close(iws)
201   format(/t1,a,/,a/)
c
      io = 12
      open(unit=io,file=inpfile,status='old')
c
c DemoFlag & statusFlag:
c
c     If DemoFlag set this program version is only for demonstration
c     purposes. Graphics files will still be generated for all
c     geometries but FE in-put deck will only be generated for the
c     basic PLATE model.
c
c     If statusFlag is set a "job.sta" will be generated otherwise not.
c     The status file is opened on UNIT=10
c
      read(io,'(a)') row
      call char_to_int(row,ivek,num)
cc      if (num.ne.2) goto 99
c Remove Demoflag-option as in input option /JF Nov. 22 1999/
c
      if (num.ne.1) goto 99
      ivek(2) = ivek(1)
      ivek(1) = 0
      if (ivek(1).eq.0) then
         demoFlag = .false.
      else
         demoFlag = .true.
      endif
c
c Open a mesh-generating status file on UNIT=10
c
ccc  Append the given filename with the "*.sta" extension
c      i = index(jobname,' ',.false.)
c      stfile = jobname(1:i)//'.sta'
      call appendFile(jobname,'.sta',stfile)
      iws = 10
      if (ivek(2).eq.0) then
         statusFlag = .false.
         open(unit=iws, status='scratch')
c         open(unit=iws, file=stfile, status='scratch')
      else
         statusFlag = .true.
         open(unit=iws, file=stfile, status='replace')
      endif
      call writeDate(iws)
c
c FE-code: prog (read and left adjust)
c
      read(io,'(a)') row
      i = 0
 10   continue
         i = i + 1
      if (row(i:i).eq.' ') goto 10
         j = i
 15   continue
         j = j + 1
      if (row(j:j).ne.' ') goto 15
      j = j - 1
      do k=1, 20
         prog(k:k) = ' '
      enddo
      prog(1:(j-i+1)) = row(i:j)
c
c Make sure prog is all upper case
c
      call ToUpper( prog )
c
c If in-put should be generated on neutral format use
c WARP3D format as default until the specific neutral format
c is generated.
c
      if     ( prog .eq. 'PATRAN' ) then
         NeutralFormat = prog
         prog = 'WARP3D    '
      elseif ( prog .eq. 'FEMAP' ) then
         NeutralFormat = prog
         prog = 'WARP3D    '
      else
         NeutralFormat = '          '
      endif
c
c Type of analysis:
c
c      analysis = 1   static,  linear elastic
c      analysis = 2   static,  nonlinear elastic
c      analysis = 3   static,  incremental elastic-plastic
c      analysis = 4 static,  incremental finite strain el.-pl.
c      analysis = 5   dynamic, linear elastic
c      analysis = 6   dynamic, nonlinear elastic
c      analysis = 7   dynamic, incremental elastic-plastic
c      analysis = 8   ...   
c      analysis = 9   ...   
c
c   Equation Solver   Direct Sparse Solver      eqsolver=1 
c                     Conjugate Gradient Solver eqsolver=2
c
      read(io,'(a)') row
      call char_to_int(row,ivek,num)
      if ((num.ne.2).and.(num.ne.3)) goto 99
      analysis = ivek(1)
      numLoadStep = ivek(2)
      if ( (num.eq.3) .and. (ivek(3).eq.1)) then
         eqsolver = 1
      else
         eqsolver = 2
      endif
c
c Material parameters; elastic: E, ny
c
c   Material parameters are stored in array dmat(*)
c
c      analysis = 1   static,  linear elastic
c        dmat(1) = no of material constants defined = 2
c        dmat(2) = Young's modulii
c        dmat(3) = Poisson's ratio
c
c      analysis = 2  static, nonlinear elastic (deformation Plasticity)
c        dmat(1) = no of material constants defined = 4 (5 if abaqus)
c        dmat(2) = Young's modulii
c        dmat(3) = Poisson's ratio
c        dmat(4) = sigY
c        dmat(5) = n = hardening
c        dmat(6) = alpha  = offset strain (Ramberg-Osgood model)
c
c      analysis = 3   static,  incremental elastic-plastic
c        dmat(1) = no of material constants defined >= 4
c        dmat(2) = Young's modulii
c        dmat(3) = Poisson's ratio
c        dmat(4) = sig1
c        dmat(5) = eps1
c        dmat(6) = sig2
c        dmat(7) = eps3
c        dmat(6) = ..
c        dmat(7) = ..
c
      read(io,'(a)') row
      call char_to_real(row,rvek,num)
      if (num.ne.1) goto 99
      dmat(1) = rvek(1)
      n = int( dmat(1) + 0.5 )
      j = 0
c
c Material data constants must be written in the in-put file
c using columns 1 to 200.
c
 20   continue
         read(io,'(a)') row
         call char_to_real(row,rvek,num)
         do i=1, num
            j = j + 1
            dmat(j+1) = rvek(i)
         enddo
      if (j.lt.n) then
         goto 20
      elseif (j.gt.n) then
         write(iws,*) 'ERROR: no. of material constants does not match'
         goto 99
      endif
c
c Geometry parameters: L, W, T, Rc, Rp, Phi0
c
      read(io,'(a)') row
      call char_to_real(row,rvek,num)
      call determ_geom(tge,w,length,t,rgeo,rvek,num, iws,errorFlag)
      if ( errorFlag .gt. 0 ) then
         write(iws,*)'Stop after error in determ_geom subroutine'
         return
      end if
c
c Crack geomtry: a, c
c
c  TCR = 1  embedded crack
c      = 2  semi elliptical crack (c - major axis = free surface)
c      = 3  semi elliptical crack (a - minor axis = free surface)
c      = 4  quarter (corner) crack (a,c - both free surfaces)
c
      read(io,'(a)') row
      call char_to_real(row,rvek,num)
      if (num.ne.4) goto 99
      tcr = int(rvek(1) + 0.5)
      a    = rvek(2)
      c    = rvek(3)
      keyhole = int( rvek(4) + 0.5 )
c
c Test overall dimensions of the body in relation to crack dimensions.
c
      if (a .gt. (0.99*t)) then
         write(iws,'(t2,a)') '>> a > 0.99*t => decrease a'
         errorFlag = 1
         return
      endif     
      if (w .lt. (c + 0.01*t)) then
         write(iws,'(t2,a)') '>> W < c + 0.01*t => increase W'
         errorFlag = 1
         return
      endif     
c
c Load definition
c
c  TBC = 1   prescribed end displacements
c  TBC = 2   prescribed end tractions and   internal pressure
c                           (cylinder, sphere, pipe)
c  TBC = 3   prescribed tractions on the crack face
c
      do i=1, 10
         rload(i) = 0
      enddo
      read(io,'(a)') row
      call char_to_real(row,rvek,num)
      tbc = int(rvek(1) + 0.5)
      do i=2, num
         rload(i-1) = rvek(i)
      enddo
c
c Mesh parameters: elemtype, mf
c
      read(io,'(a)') row
      call char_to_real(row,rvek,num)
      if (num.eq.1) then
         etyp = int(rvek(1)+0.5)
c . . . Default values
         mf = 10
         if (keyhole.ne.0) then
            mr = 7
         else
            mr = 5
         endif
      elseif (num.eq.2) then
         etyp = int(rvek(1)+0.5)
         mr   = int(rvek(2)+0.5)
      elseif (num.eq.3) then
         etyp = int(rvek(1)+0.5)
         mr   = int(rvek(2)+0.5)
         mf   = int(rvek(3)+0.5)
      else
         goto 99
      endif
c
      close(io)
c
c Check input data, so that parameters are within permitted limits
c
c    =>   To be implemented
c
c Determine model parameters
c
c  # of elements
c
c  Zone S
c
      sfred_type = 1
      sjred_type = 1
      na = 12
      nb = 8
c     na = 12*(c/a)**0.15
      na = nb + (na-nb)*(c/a)**0.5       
      if (mod(na,2) .ne. 0) na = na+1
c
c if s = 0 is a free surface => increase  na
c
      if ((tcr.eq.3) .or. (tcr.eq.4)) na = na + 2      
      sfred = mr - 2
      if (sfred.lt.2) sfred = mr - 1
      call determ_mv_mh(mf,mv,mh,sfred_type)
c
c  Zone A
c
      m1min = 2
      m2min = 3
      m12m  = 4
      psi = a/t
      if (psi .le. 0.5) then
      m12e  = 5
         m1  = (1.-2.*psi)*real(m1min) + 2.*psi*(real(m12m)+0.99)
         m12 = 2*int((1.-2.*psi)*real(m12e) + 2.*psi*(real(m12m)+0.99))
         m2  = m12 - m1
         if (mod((m1+m2),2).ne.0) m2 = m2
      else
      m12e  = 6
         psi = (1.-a/t)
         m2  = (1.-2.*psi)*real(m2min) + 2.*psi*(real(m12m)+0.99)
         m12 = 2*int((1.-2.*psi)*real(m12e) + 2.*psi*(real(m12m)+0.99))
         m1  = m12 - m2
         if (mod((m1+m2),2).ne.0) m1 = m1
      endif
      ma = m1 + m2 + 2*mh
c
c Geometry
c
c  Zone S
c
      alfa = sqrt(c*c-a*a)/2.
      beta0 = 0.35
      betam = 0.50
      psi  = abs(2.*(a/t)-1.)
      beta = (1.0-psi)*betam + psi*beta0
      c1 = 1.0 / ( 1.0 + beta*real(m1min) )
      c2 = 4.0 / ( 2.0 + beta*real(2*m12m) ) - 2.0*c1
      psi = a/t
      if (psi .gt. 0.5 ) psi = 1.0 - psi
      b = t*psi*(c1 + c2*psi)
c
 30   continue
c
c  Determine the bias in the focused mesh
      bias(1) =  bias_crack(b,b*beta,mr)
      xi = 1 - mr
      rfm(1) = beta*b*bias(1)**xi
      do i=2, mr
         xi = i - 1
         rfm(i) = rfm(i-1) + rfm(1)*bias(1)**xi
      enddo
c
      bias(2) = 1.2
      xi = mv
      yg(0) = 0.0
      yg(1) = b * ( bias(2)-1. ) / ( bias(2)**xi-1. )
      do i=2, mv
         xi = i-1
         yg(i) = yg(i-1) + yg(1)*bias(2)**xi
      enddo
c
c  Zone A
c
      bias(4) = 1.5
      dx0 = yg(mv)-yg(mv-1)
      dxmax = 7.5*t / real(ma/2)
      if (dxmax.gt.t) then
         if ( (tge.eq.201).or.(tge.eq.202).or.(tge.eq.301).or.
     &        (tge.eq.302).or.(tge.eq.305).or.(tge.eq.306) ) then
            dxmax = t
         endif
      endif
      xmax = length - yg(mv-1)
      nmin = 1
      call determ_remote_mesh(nmin,dx0,dxmax,xmax,bias(4),x,n)
c      lt = n + mv-1
      lt = n + mv
      do i=1, n
         j = i+mv
         yg(j) = yg(mv) + x(i)
      enddo
c
      do i=0, lt
         zg(1,i) = 0.0
         zg(2,i) = t
      enddo
      if ( length .lt. t) then
         rtype = 0
      elseif (length .lt. (1.6*t)) then
         rtype = 1
         lred = 0
         do i=mv, lt
            if ((yg(i).gt.(0.5*t)).and.(lred.eq.0)) lred = i
         enddo
      else
         rtype = 2
         lred = 0
         do i=mv, lt
            if ((yg(i).gt.(0.6*t)).and.(lred.eq.0)) lred = i
         enddo
      endif
      if (tge .gt. 200) then
         rtype = 2
         lred = 0
         do i=mv, lt
            if ((yg(i).gt.(0.5*t)).and.(lred.eq.0)) lred = i
         enddo
      endif
c
c If rtype = 0 no mesh coarsening
c    rtype = 1 mesh coarsening 2 to 1 in Zone A
c              (one extra layer is included at the red. site)
c    rtype = 2 mesh coarsening 4 to 1 in Zone A and 2 to 1 in Zone B
c              (currently not used)
c
      if (rtype.eq.0) then
         kma=2*lt+1
      elseif (rtype.ge.1) then
         lt=lt+1
         kma=2*lt+1
         do i=lt, lred+1, -1
            yg(i)   = yg(i-1)
            zg(1,i) = zg(1,i-1)
            zg(2,i) = zg(2,i-1)
         enddo
         yg(lred)   = 0.5*( yg(lred-1)   + yg(lred+1) )
         zg(1,lred) = 0.5*( zg(1,lred-1) + zg(1,lred+1) )
         zg(2,lred) = 0.5*( zg(2,lred-1) + zg(2,lred+1) )
      endif
      if (rtype.eq.2) then
         lt=lt+1
         kma=2*lt+1
         do i=lt,lred+3, -1
            yg(i)   = yg(i-1)
            zg(1,i) = zg(1,i-1)
            zg(2,i) = zg(2,i-1)
         enddo
         yg(lred+2)   = 0.5*( yg(lred+1)   + yg(lred+3) )
         zg(1,lred+2) = 0.5*( zg(1,lred+1) + zg(1,lred+3) )
         zg(2,lred+2) = 0.5*( zg(2,lred+1) + zg(2,lred+3) )
      endif
c
c  Mesh parameters in region X > c+b; Zone A and Zone B
c
      psi = (a+b)/t 
      lambda0 = 1.8
      lambdam = 1.1
      lambda1 = 1.8
      if (psi .lt. 0.5) then
         lambda = (1.-2.*psi)*lambda0 + 2.*psi*lambdam
      else
         lambda = 2.*(1.-psi)*lambdam + (2.*psi-1.)*lambda1
      endif
      bias(5) = lambda
c
c  define element lengths in X-dir. in zone A & B 
c
      if ( ((a+b)/t) .lt. 0.45 ) then
         nmin = m2
      else
         mb1 = 2
         mb  = int( 2.*(1.-psi) + (2.*psi-1.)*real(mb1) + 0.5 )
         nmin = m2 + mb
      endif
c   decrease b if W < tmp
      tmp = c + b + beta*b*real(nmin)*0.8
      if ( w .lt. tmp ) then
ccc         write(*,'(t2,a,g11.4)') 'old b =',b
         write(iws,'(t2,a,g11.4)') 'old b =',b
         b = 0.99*(w-c) / ( 1.0 + beta*real(nmin)*0.8 )
ccc         write(*,'(t2,a,g11.4)') 'new b =',b
         write(iws,'(t2,a,g11.4)') 'new b =',b
         goto 30 
      endif
c
c  If TGE > 200  => Relate dxmax to the inner radius of the Cylinder
c
      dx0 = rfm(mr) - rfm(mr-1)
      if ( rtype .lt. 2) then
         dxmax = 7.5*t / real(nb)
      else
         dxmax = 7.5*t / real(nb/2)
      endif
      if (dxmax.gt.t) then
         if ( (tge.eq.203).or.(tge.eq.204).or.(tge.eq.303).or.
     &        (tge.eq.304).or.(tge.eq.307).or.(tge.eq.308) ) then
            dxmax = t
         endif
      endif
      xmax = w - (c + b)
      call determ_remote_mesh(nmin,dx0,dxmax,xmax,bias(5),x,n)
      xag(0) = c + b
      do i=1, n
         xag(i) = xag(0) + x(i)
      enddo
      mb = n - m2
      do i=m2, n
         xbg(i-m2) = xag(i)
      enddo
c
c  Indicies in node number arrays
c
c  Zone S:
c
      ims = 2*mf+1
      jms = 2*na*sjred_type+1
      kms = 2*mr + 1
      if ( (sfred_type.eq.1) .and. (sjred_type.gt.1) ) then
         kms = kms + 2
      elseif ( (sfred_type.gt.1) .and. (sjred_type.eq.1) ) then
         kms = kms + 2
      elseif ( (sfred_type.gt.1) .and. (sjred_type.gt.1) ) then
         kms = kms + 4
      endif
      ksr1 = 2*sfred+1
c
c  Zone A:
      ima = 2*ma+1
      jma = 2*na+1
      kar1 = 2*lred+1
      kar2 = 2*(lred+2)+1
c
c  Zone B:
      if (mb.gt.0) then
         imb=2*mb+1
      else
         imb = 1
      endif
      jmb=2*nb+1
c
c TEST
c
c      write(iws,'(t3,a)') 'Inside read_input'
c      write(iws,'(t3,a)') 'BIAS'
c      do i=1, 8
c         write(iws,'(g15.6)') bias(i)
c      enddo
c      write(iws,'(t3,a)') 'rfm'
c      do i=1, mr
c         write(iws,'(g15.6)') rfm(i)
c      enddo
c      write(iws,'(t3,a)') 'yg, zg'
c      do i=0, lt
c         write(iws,'(3g15.6)') yg(i),zg(1,i), zg(2,i)
c      enddo
c      write(iws,'(t3,a)') 'xag'
c      do i=1, m2
c         write(iws,'(g15.6)') xag(i)
c      enddo
c      write(iws,'(t3,a)') 'xbg'
c      do i=1, mb
c         write(iws,'(g15.6)') xbg(i)
c      enddo
c
c  Set size of matrices 
c
      ns1 = ims
      ns2 = jms
      ns3 = kms
      if ((ns1*ns2*ns3).gt.nsm) then
         write(iws,'(t1,a)') '>> increase the size of array nods(i,j,k)'
         write(iws,'(t4,a,i8)') 'current  size = ', nsm
         write(iws,'(t4,a,i8)') 'required size = ', ns1*ns2*ns3
         errorFlag = 1
         return
      endif
      do i=1, nsm
         nods(i) = 0
      enddo
c
      na1 = ima
      na2 = jma
      na3 = kma
      if ((na1*na2*na3).gt.nam) then
         write(iws,'(t1,a)') '>> increase the size of array noda(i,j,k)'
         write(iws,'(t4,a,i8)') 'current  size = ', nam
         write(iws,'(t4,a,i8)') 'required size = ', na1*na2*na3
         errorFlag = 1
         return
      endif
      do i=1, nam
         noda(i) = 0
      enddo
c
      nb1 = imb
      nb2 = jmb
      nb3 = kma
      if ((nb1*nb2*nb3).gt.nbm) then
         write(iws,'(t1,a)') '>> increase the size of array nodb(i,j,k)'
         write(iws,'(t4,a,i8)') 'current  size = ', nbm
         write(iws,'(t4,a,i8)') 'required size = ', nb1*nb2*nb3
         errorFlag = 1
         return
      endif
      do i=1, nbm
         nodb(i) = 0
      enddo
c
c  Check if maximum size of npos (coordinate array) is sufficient
c
      max_node_no = jms*ims*kms + (jma*ima+jmb*imb)*kma
      if (max_node_no.gt.inm) then
        write(iws,'(t1,a)') '>> increase the size of array npos(inm,3)'
        write(iws,'(t4,a,i8)') 'current  size = ', inm
        write(iws,'(t4,a,i8)') 'required size = ', max_node_no
        errorFlag = 1
        return
      endif
c
      nes2 = (kms+1)/2
      if (nes2.gt.nesm) then
        write(iws,'(t1,a)') '>> increase the array size parameter nesm'
        write(iws,'(t4,a,i8)') 'current  size = ', nesm
        write(iws,'(t4,a,i8)') 'required size = ', nes2
        errorFlag = 1
        return
      endif
      nea2 = (kma+1)/2
      if ( (nea2*(jma+1)/2) .gt. (nesm*nesm) ) then
        write(iws,'(t1,a)') '>> increase the array size parameter nesm'
        write(iws,'(t4,a,i8)') 'current  size = ', nesm
        tmp = real(nea2*(jma+1)/2)
        write(iws,'(t4,a,i8)') 'required size = ', int( sqrt(tmp) )+1
        errorFlag = 1
        return
      endif
      neb2 = (kma+1)/2
      if ( (neb2*(imb+1)/2) .gt. (nesm*nesm) ) then
        write(iws,'(t1,a)') '>> increase the array size parameter nesm'
        write(iws,'(t4,a,i8)') 'current  size = ', nesm
        tmp = real(neb2*(imb+1)/2)
        write(iws,'(t4,a,i8)') 'required size = ', int( sqrt(tmp) )+1
        errorFlag = 1
        return
      endif
      if ( (neb2) .gt. (nesm) ) then
        write(iws,'(t1,a)') '>> increase the array size parameter nesm'
        write(iws,'(t4,a,i8)') 'current  size = ', nesm
        write(iws,'(t4,a,i8)') 'required size = ', neb2
        errorFlag = 1
        return
      endif
      do i=1, nesm
         do j=1,2
            estk_s(j,i) = 0
         enddo
      enddo
      do i=1, nesm2
         do j=1,2
            estk_a(j,i) = 0
            estk_b(j,i) = 0
         enddo
      enddo
c
      return
c
99    write(iws,'(t10,a/t10,a)')
     1  '>> the number of data given does not match the number of',
     2  '   data required at each line, => correct the in-put file!!'
      close(io)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine determ_geom(tge,w0,l0,t0,rgeo,rvek,num, iws,errorFlag)
c
c Determine the basic plate geometric parameters (W,L,t), given TGE
c
      implicit none
      integer          tge,num,iws,errorFlag
      double precision w0,l0,t0,rgeo(*),phi0,rvek(*), l,rc,rp,pi,one
      parameter (one = 1.0)
c
      tge = int(rvek(1) + 0.5)
      if ( (tge.lt.299) .and. (num.ne.4) ) then
         write(iws,'(t5,a)') '>> IN-PUT ERROR: geometric parameters'
         errorFlag = 1
         return
      elseif ( (tge.gt.299) .and. (num.ne.5) ) then
         write(iws,'(t5,a)') '>> IN-PUT ERROR: geometric parameters'
         errorFlag = 1
         return
      endif
c
      pi = 4.0*atan(one)
c
c  Basic Plate
c
      if (tge.eq.101) then
         w0 = rvek(2)
         l0 = rvek(3)
         t0 = rvek(4)
c
c  Cylinder
c
      elseif ((tge.gt.200).and.(tge.lt.300)) then
         rc = rvek(2)
         l  = rvek(3)
         t0 = rvek(4)
         if (tge.eq.201) then
            l0 = pi * rc
            w0 = l
         elseif (tge.eq.202) then
            l0 = pi * ( rc + t0 )
            w0 = l
         elseif (tge.eq.203) then
            l0 = l
            w0 = pi * rc
         elseif (tge.eq.204) then
            l0 = l
            w0 = pi * ( rc + t0 )
         else
            write(iws,'(t5,a)') '>> ERROR: cylider input'
            errorFlag = 1
            return
         endif
         rgeo(4) = rc
         rgeo(5) = l
c
c  Pipe - Elbow
c
      elseif ((tge.gt.300).and.(tge.lt.400)) then
         rc   = rvek(2)
         t0   = rvek(3)
         rp   = rvek(4)
         phi0 = rvek(5)
         if (tge.eq.301) then
            l0 = pi * rc
            w0 = phi0 * ( rp + rc )
         elseif (tge.eq.302) then
            l0 = pi * ( rc + t0 )
            w0 = phi0 * ( rp + rc + t0 )
         elseif (tge.eq.303) then
            l0 = phi0 * ( rp + rc )
            w0 = pi * rc
         elseif (tge.eq.304) then
            l0 = phi0 * ( rp + rc + t0 )
            w0 = pi * ( rc + t0 )
         elseif (tge.eq.305) then
            l0 = pi * rc
            w0 = phi0 * ( rp - rc )
         elseif (tge.eq.306) then
            l0 = pi * ( rc + t0 )
            w0 = phi0 * ( rp - rc - t0 )
         elseif (tge.eq.307) then
            l0 = phi0 * ( rp - rc )
            w0 = pi * rc
         elseif (tge.eq.308) then
            l0 = phi0 * ( rp - rc - t0 )
            w0 = pi * ( rc + t0 )
         else
            write(iws,'(t5,a)') '>> ERROR: pipe-elbow input'
            errorFlag = 1
            return
         endif
         rgeo(4) = rc
         rgeo(5) = rp
         rgeo(6) = phi0
      endif
c
      rgeo(1) = w0
      rgeo(2) = l0
      rgeo(3) = t0
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine determ_mv_mh(mf,mv,mh,sfred_type)
      implicit none
      integer  mf,mv,mh,sfred_type, mfa
c
      if (sfred_type.eq.1) then
         if (mod(mf,2).ne.0) mf = mf + mod(mf,2)
         mfa = mf
      elseif (sfred_type.eq.2) then
         if (mod(mf,4).ne.0) mf = mf + mod(mf,4)
         mfa = mf/2
      elseif (sfred_type.eq.3) then
         if (mod(mf,6).ne.0) mf = mf + mod(mf,6)
         mfa = mf/3
      endif
      if (mfa.lt.8)  mfa=8
      if (mfa.gt.20) mfa=20
      if     (mfa .le. 8) then
         mv = 2
      elseif (mfa .le. 12) then
         mv = 3 
      elseif (mfa .lt. 16) then
         mv = 4
      else
         mv = 5
      endif
      mh = (mfa-mv-mv)/2
      mf = sfred_type*mfa
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      double precision function bias_crack(xtot,xn,n)
c
c Given: Xtot, xn, and n;  Xtot = x1+x2+...+xn
c Find : bias (=b below) in a geometric series  xn > x1
c 
      implicit none
      integer  n,i,it
      double precision xtot,xn,rtot,b,db,tol,xi,r,rold
      parameter( tol = 1.e-6)
c
      rtot = xtot/xn
      b  = 0.8
      db = 0.2
      r  = 0.0
      it = 0
10    continue
         it = it + 1
         b = b + db
         rold = r
         r = 1.0
         do i=2, n
            xi = 1 - i
            r = r + b**xi
         enddo
         if (it.eq.1) goto 10
         if (abs((rtot-r)/rtot) .gt. tol) then
            if ( ((rtot-r)*(rtot-rold)) .lt. 0) db = - db/2.0
            goto 10
         endif
      bias_crack = b
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine determ_remote_mesh(nmin,dx0,dxmax,xm,bias,x,n)
      implicit none
      integer  nmin,n, i,k
      double precision dx0,dxmax,xm,bias,x(*),xnmin,xstar,dx,dx1,xi
      double precision calc_lambda
c
c nmin = minimum elements allowed
c nmax = maximum elements allowed (currently not implemented)
c dx0 = initial minimum reference length 
c dxmax = maximum element length allowed
c xm = total length to be divided into segments
c bias = initial ratio between adjacent elements 
c
      if ((xm/real(nmin)) .lt. dx0) then
         n = nmin
         do i=1, n
            x(i) = xm*real(i)/real(n)
         enddo
         return
      endif
c
      dx1 = bias*dx0
      xi = nmin
      xnmin = dx1*(bias**xi-1.0)/(bias-1.0)
c
      if (xnmin .ge. xm) then
         n = nmin
         bias = exp( calc_lambda( dx0, xm+dx0, n+1) )
         dx1 = dx0*bias
         x(1) = dx1
         do i=2, n
            xi = i - 1
            x(i) = x(i-1) + dx1*bias**xi
         enddo
         return
      endif
c
      k = 1
      x(k) = dx0*bias
10    continue
         k = k + 1
         xi = k
         dx = dx0*bias**xi
         x(k) = x(k-1) + dx
      if (dx .lt. dxmax) goto 10
      xstar = x(k-1)
c
      if (xstar .gt. xm) then
         n = 0
20       continue
            n = n + 1
         if (x(n) .lt. xm) goto 20
         bias = exp( calc_lambda( dx0, xm+dx0, n+1) )
         dx1 = dx0*bias
         x(1) = dx1
         do i=2, n
            xi = i - 1
            x(i) = x(i-1) + dx1*bias**xi
         enddo
         return
      endif
c
      n = k
30    continue
         n = n + 1
         xstar = xstar + dxmax
      if (xstar .le. (xm+0.2*dxmax)) goto 30
      call calc_bias_int1(dx0,xm,n+1,k+1,bias)
      dx1 = dx0*bias
      x(1) = dx1
      do i=2, n
         xi = i-1
         if ( i .ge. k ) xi = k - 1
         x(i) = x(i-1) + dx1*bias**xi
      enddo
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine Xdeterm_remote_mesh(bias,dx1,dxmax,xmax,x,n)
      implicit none
      integer  i,k,n
      double precision bias,dx1,dxmax,xmax,x(*),dx(200),xi
      double precision calc_lambda
      i = 1
      x(i) = dx1
10    continue
         i = i + 1
         xi = i - 1 
         dx(i) = dx1*bias**xi
         x(i) = x(i-1) + dx(i)
      if (dx(i) .lt. dxmax) goto 10
      k = i - 1
c
c Case 1: Xgradation < Xmax
c
      if (x(k) .lt. xmax ) then
         i = k
20       continue
            i = i + 1
            x(i) = x(i-1) + dxmax
         if (x(i) .lt. xmax) goto 20
         n = i
30       continue
            xi = n - k
            dx(k+1) = (xmax-x(k)) / xi
         if (dx(k+1) .lt. dx(k)) then
            k = k - 1
            goto 30
         endif
         do i=k+1, n
            x(i) = x(i-1) + dx(k+1)
         enddo 
      else
c
c Case 2: Xgradation > Xmax
c
         i = k
40       continue         
            i = i - 1
         if (x(i) .gt. xmax) goto 40
         if ((xmax-x(i)) .le. (0.1+0.025*i*dx(i)) ) then
            n = i
         else
            n = i + 1
         endif
         bias = exp(calc_lambda(dx1,xmax,n))
         do i=2, n
            xi = i-1
            x(i) = x(i-1) + dx1*bias**xi
         enddo
      endif
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine char_to_real(row,rvek,num)
c  The routine transform a character to a real number
      implicit   none
      integer    num,n,k,i,exponent,es, nr
      double precision  rvek(*),sign,xi,xe
      character  row*200
c
c  nr = to maximum position of characters in row
c
      nr = 200
      num=0
      n=1
 5    continue
         i=0
         k=0
         exponent=0
         es=1
         sign=1.0
 10      continue
c    Searcing for the 1:st digit !
         if ( ( (row(n:n).lt.'0').or.(row(n:n).gt.'9') ).and.
     &          (row(n:n).ne.'.') )  then
            n=n+1
            if (n.gt.nr) goto 99
            goto 10
         endif
c    Negative sign !
         if (row((n-1):(n-1)).eq.'-') sign=-1.0
 20      continue
c    Integer part
         if ( (row(n:n).ge.'0').and.(row(n:n).le.'9') ) then
            i=10*i+ichar(row(n:n))-48
            n=n+1
            if (n.gt.nr) goto 99
            goto 20
         endif
c
         if (row(n:n).ne.'.') goto 40
         n=n+1
         if (n.gt.nr) goto 99
 30      continue
c    Decimal part !
         if ( (row(n:n).ge.'0').and.(row(n:n).le.'9') ) then
            i=10*i+(ichar(row(n:n))-48)
            n=n+1
            if (n.gt.nr) goto 99
            k=k-1
            goto 30
         endif
 40      continue
         if ( (row(n:n).eq.'E').or.(row(n:n).eq.'e').or.
     &          (row(n:n).eq.'D').or.(row(n:n).eq.'d') ) then
c    Exponent
            n=n+1
            if (n.gt.nr) goto 99
            if (row(n:n).eq.'+') then
c    Positive or negative sign !
               es=1
               n=n+1
               if (n.gt.nr) goto 99
            elseif (row(n:n).eq.'-') then
               es=-1
               n=n+1
               if (n.gt.nr) goto 99
            endif
 50         if ( (row(n:n).ge.'0').and.(row(n:n).le.'9') ) then
               exponent=10*exponent+ichar(row(n:n))-48
               n=n+1
               if (n.gt.nr) goto 99
               goto 50
            endif
         endif
         xi = 1000*i
         k = k - 3
c         xi = i
c         k = k
         if (es.gt.0) then
            xe = k + exponent
         else
            xe = k - exponent
         endif
         num=num+1
      if (sign.gt.0) then
         rvek(num) =  xi*10.0**xe
      else
         rvek(num) = -xi*10.0**xe
      endif
      if (num.lt.nr) goto 5
 99   return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine char_to_int(row,ivek,num)
c-- The routine transform a character to an integer number
      integer    ivek(*),num,n,i
      character  row*200
      logical    minus
c
c  nr = to maximum position of characters in row
c
      nr = 200
      num=0
      n=1
5     continue
         i=0
         minus = .false.
10        continue
c... Search for the 1:st digit
         if ( (row(n:n).lt.'0') .or. (row(n:n).gt.'9') ) then
            n=n+1
            if (n.gt.nr) goto 99 
            goto 10
         endif
c... Negative sign !
         if (row((n-1):(n-1)).eq.'-') minus=.true.
20       continue
         if ( (row(n:n).ge.'0').and.(row(n:n).le.'9') ) then
            i=10*i+ichar(row(n:n))-48
            n=n+1
            if (n.gt.nr) goto 99 
            goto 20
         endif
         num=num+1
         ivek(num)=i
         if (minus) ivek(num)= - ivek(num)
      if (num.lt.nr) goto 5
99    return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
