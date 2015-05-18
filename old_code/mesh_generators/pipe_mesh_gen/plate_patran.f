c----67--1---------2---------3---------4---------5---------6---------712
c
c    File: plate_patran.f   Written on  July 7 1998
c                           Modified on July 8 1998
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine patran_neu(etyp,elnum,no_of_nodes,dmat,rgeo,rload,
     &           analysis,tcr,tge,tbc,ecrsur,ecrtip,ecsym,etop,
     &           erem,esid, nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &           nodb,nb1,nb2,nb3, jobname,iws,errorFlag)
c
c The routine generates an input deck for in the Neutral Patran format
c
      implicit none
c
      integer  analysis,tcr,tge,tbc,etyp,elnum,no_of_nodes,
     &         ecrsur(2,*),ecrtip(2,*),ecsym(2,*),etop(2,*),
     &         erem(2,*),esid(2,*),iws,errorFlag
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     &         na1,na2,na3,noda(na1,na2,na3),
     &         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      double precision  dmat(*),rgeo(*),rload(*)
      character  jobname*200, patfile*200
c
      integer    ndoftot,io, eset(2,2),no_eset
      parameter (ndoftot = 50000, io=11)
      integer           id(3,ndoftot),i
      double precision  rd(3,ndoftot),press
      logical    ip
c
      integer      mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
c  Set load arrays id(3,*), rd(3,i) to zero.
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
c Define Boundary Conditions
c
c Nodal Constraints and Equivalent Nodal forces are stored in	rd(3,i)
c
c  TBC = 1  Displacements on Remote Boundary
c  TBC = 2  Tractions on Remote Boundary (includes Internal Pressure)
c  TBC = 3  Tractions on Crack Face
c
      if (tbc.eq.1) then
         call wrp_pre_disp(tcr,tge,rgeo,rload,id,rd,
     &            noda,na1,na2,na3, nodb,nb1,nb2,nb3, iws,errorFlag)
      elseif ( (tbc.eq.2) .or. (tbc.eq.3) ) then
         call wrp_pre_traction(tge,tcr,tbc,rgeo,rload,ip,id,rd,
     &            ecrsur,ecrtip,erem,esid,no_of_nodes,iws,errorFlag)
      else
         write(iws,*)'>>ERROR: defining TBC in "SUBR. warp3d_inp"'
         errorFlag = 1
         return
      endif
      call wrp_fix_disp(analysis,tcr,tge,id,rd, nods,ns1,ns2,ns3,
     &         noda,na1,na2,na3, nodb,nb1,nb2,nb3,iws,errorFlag)
c
c Start Generate the in-put deck on NEUTRAL PATRAN FORMAT
c
c PATRAN 2.5 Neutral File Packet Header Contains
c
c  IT   ID   IV   KC   N1   N2   N3   N4   N5   format(I8,8I8)
c
c  IT = Packet Type
c  ID = Identification number. A "0" ID value means not applicable (n/a)
c  IV = Additional ID. A "0" value means not applicable (n/a)
c  KC = Card count (number of data cards after the header)
c  N1 to N5 = Supplemental integer values
c
c  The Neutral File Packet Header should be given in format(I8,8I8)
c
      call appendFile( jobname, '_pat.neu', patfile )
      open( unit=io, file=patfile, status='unknown' )
c
c PACKET Type 25: Title Card
c
      call patran_25(io)
c
c PACKET Type 26: Summary Data
c
      call patran_26(io,no_of_nodes,elnum)
c
c PACKET Type 01: Node Data
c 
      call patran_01(io)
c
c PACKET Type 02: Element Data
c
      call patran_02(io,etyp,elnum)
c
c PACKET Type 03: Material Properties
c
      call patran_03(io,analysis,dmat,iws)
c
c PACKET Type 04: Element Properties
c
      call patran_04(io)
c
c PACKET Type 05: Coordinate frames
c

c
c PACKET Type 06: Distributed Loads
c
      if (ip) then
         press = rload(7)
         if    ((tge.eq.201).or.(tge.eq.203)) then
            no_eset = 2
            call patran_06(io,no_eset,ecrsur,ecsym,press)
         elseif ( (tge.eq.202).or.(tge.eq.204) ) then
            no_eset = 1
            call patran_06(io,no_eset,etop,eset,press)
c
         elseif ( (tge.eq.301).or.(tge.eq.303).or.
     &         (tge.eq.305).or.(tge.eq.307) ) then
            no_eset = 2
            call patran_06(io,no_eset,ecrsur,ecsym,press)
         elseif ( (tge.eq.302).or.(tge.eq.304).or.
     &         (tge.eq.306).or.(tge.eq.308) ) then
            no_eset = 1
            call patran_06(io,no_eset,etop,eset,press)
         endif
      endif
c
c PACKET Type 07: Node Forces
c
      call patran_07(io,no_of_nodes,id,rd)
c
c PACKET Type 08: Node Displacements
c
      call patran_08(io,no_of_nodes,id,rd)
c
c PACKET Type 10: Node Temperatures
c

c
c PACKET Type 11: Element Temperatures
c

c
c PACKET Type 99: Termination Card
c
      write(io,'(i2,8i8)') 99, 0,0,0,0,0,0,0,0
c
      close(io)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine patran_25(io)
c
c PACKET Type 25: Title Card
c
c  Header: IT,ID,IV,KC, 0,0,0,0,0
c
c   IT = 1
c   ID = 0 n/a
c   IV = 0 n/a
c   KC = Card count (number of data cards after the header)
c
c  Card 1: TITLE
c
c   TITLE = Identifying title may contain up to 80 characters
c
      implicit none
      integer     io,it,id,iv,kc
      character*4 head(20)
      data  head/'3-D ','Surf','ace ','crac','k FE',
     &           '-mes','h   ','    ','    ','    ',
     &           '    ','    ','    ','    ','    ',
     &           '    ','    ','    ','    ','    '/
c
      it = 25
      id = 0
      iv = 0
      kc = 1
      write(io,101) it, id,iv,kc, 0,0,0,0,0
      write(io,102) head
101   format(i2,8i8)
102   format(20a4)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine patran_26(io,no_of_nodes,elnum)
c
c PACKET Type 26: Summary Data
c
c  Header: IT,ID,IV,KC,N1,N2,N3,N4,N5
c
c   IT = 1
c   ID = 0 n/a
c   IV = 0 n/a
c   KC = Card count (number of data cards after the header)
c   N1 = Number of Nodes,
c   N2 = Number of Elements
c   N3 = Number of Materials,
c   N4 = Number of Element Properties
c   N5 = Number of Coordinate Frames (use 1)
c        (the default Coordinate Frame = 0 - Cartesian Coord. Frame)
c
c  Card 1: DATE  TIME  VERSION
c
c    DATE    = Date neutral file was created
c    TIME    = Time neutral file was created
c    VERSION = PATRAN release number
c
      implicit none
      integer     io,no_of_nodes,elnum,it,id,iv,kc,n(5),i
      character   day*12,tim*8,versn*12
c
      it = 26
      id = 0
      iv = 0
      kc = 1
      n(1) = no_of_nodes
      n(2) = elnum
      n(3) = 1
      n(4) = 1
      n(5) = 1
      versn= '    3.0     '
      call get_time_date(day,tim)
      write(io,101) it, id,iv,kc,(n(i),i=1,5)
      write(io,102) day,tim,versn
101   format(i2,8i8)
102   format(a12,a8,a12)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine get_time_date(day,tim)
c
c  Find the current time and date
c
      character     day*12,tim*8, a*20,b*20,c*20
      character*4   month(12)
      integer       value(10), i,j,k
c
      call date_and_time(a,b,c,value)
c
c  Date
c
      i = value(2)
      j = value(3)
      k = value(1)
      if (j.lt.10) then
         write(day,'(a4,a1,i1,a1,i4,a1)') month(i),' ',j,' ',k,' '
      else
         write(day,'(a4,i2,a1,i4,a1)')    month(i),j,' ',k,' '
      endif
      data month /'Jan ','Feb ','Mar ','Apr ',
     &            'May ','Jun ','Jul ','Aug ',
     &            'Sep ','Oct ','Nov ','Dec '/
c
c  Time
c
      i = value(5)
      j = value(6)
      k = value(7)
      tim = '00:00:00'
      if (i.lt.10) then
         write( tim(2:2), '(i1)' ) i
      else
         write( tim(1:2), '(i2)' ) i
      endif
      if (j.lt.10) then
         write( tim(5:5), '(i1)' ) j
      else
         write( tim(4:5), '(i2)' ) j
      endif
      if (k.lt.10) then
         write( tim(8:8), '(i1)' ) k
      else
         write( tim(7:8), '(i2)' ) k
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine patran_01(io)
c
c PACKET Type 01: Node Data
c 
c  Header: IT,ID,IV,KC, 0,0,0,0,0
c
c   IT = 1
c   ID = Node ID (=node number)
c   IV = 0 n/a
c   KC = Card count (number of data cards after the header)
c
c  Card 1: X,Y,Z   (Cartesian Coordinates of Node)
c
c  Card 2: ICF, GTYPE, NDF, CONFIG, CID, PSPC
c
c   ICF*   =  Condensation flag (use 1)
c   GTYPE  =  Node type
c   NDF*   =  Number of degrees-of-freedom
c   CONFIG =  Node configuration       (use 0)
c   CID    =  Coordinate for analysis  (use 0)
c   PSPC*  =  6 permanent single point constraint flags 0 or 1
c    * These parameters are not currently used
c        ( Allways give three coordinates, i.e. x,y,z, )
c        ( e.g. a 2D problem  x, y, 0.0                )
c
      implicit none
      include 'plate_common_nod.f'
      integer          io, it,id,iv,kc,icf,ndf,config,cid, i
      double precision x,y,z
      character        gtype*1
      it = 1
      iv = 0
      kc = 2
      icf = 1
      gtype  = 'G'
      ndf    = 3
      config = 0
      cid    = 0
      do i=1, inm
         if (nnr(i).gt.0) then
            id = nnr(i)
            x = npos(i,1)
            y = npos(i,2)
            z = npos(i,3)
            write(io,101) it, id,iv,kc,0,0,0,0,0
            write(io,102) x, y, z
            write(io,103) icf,gtype,ndf,config,cid,  0,0,0,0,0,0
         endif
      enddo
101   format(i2,8i8)
102   format( 3e16.9 )
103   format( i1, a1, 3i8, 2x, 6i1 )
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine patran_02(io,etyp,elnum)
c
c PACKET Type 02: Element Data
c
c  Header: IT,ID,IV,KC, 0,0,0,0,0
c
c   IT = 2
c   ID = Element ID (=element number)
c   IV = Shape/Type of element
c      =2 Bar;         =3 Triangle; =4 Quadrilateral;
c      =5 Tetrahedral; =6 Pyramid;  =7 Wedge
c      =8 Hexahedral (=current application)
c   KC = Card count (number of data cards after the header)
c
c  Card 1: NODES, CONFIG, PID,  CEID,  q1,  q2,  q3
c
c   NODES  = total number of nodes
c   CONFIG = Element configuration (in case of several materials)
c   PID    = Property ID (+) or Material ID (-)
c   CEID   = Congruent element ID
c   q1,q2,q3 = 0.0  (currently not used only for bar elements)
c         (  PID  - element property set (we don't use)  )
c         (  CEID - element coordinate system (we don't use)  )
c
c
      implicit none
      include 'plate_common_nod.f'
      include 'plate_common_eln.f'
      integer  io,etyp,elnum
      integer  it,id,iv,kc,nodes,config,pid,ceid, i,j
c
      it = 2
      iv = 8
      if (etyp.le.10) then
         kc = 2
      elseif (etyp.le.20) then
         kc = 3
      else
         kc = 4
      endif
      nodes  = etyp
      config = 0
      pid    = 1
      ceid   = 0
      do i=1, elnum
         id = i
         write(io,101) it, id,iv,kc, 0,0,0,0,0
         write(io,102) etyp, config, pid, ceid, 0.0, 0.0, 0.0
         write(io,103) ( nnr(eln(i,j)), j=1, etyp )
      enddo
101   format(i2,8i8)
102   format( 4i8, 3e16.9 )
103   format( 10i8 )
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine patran_03(io,analysis,dmat,iws)
c
c PACKET Type 03: Material Properties
c 
c  Header: IT,ID,IV,KC, 0,0,0,0,0
c
c   IT = 3
c   ID = Material ID
c   IV = Material Type (1,2, ..., 13) (=1 Isotropic)
c   KC = Card count (number of data cards after the header)
c
c  Card 1: DATA
c
c   DATA(96)  96 Material properties
      implicit none
      integer          io,iws, analysis,it,id,iv,kc, i
      double precision dmat(*),mat(96),
     &                 elastmod,poisson,shearmod,density,reftemp
      it = 3
      id = 1
      iv = 1
      kc = 20
      elastmod = dmat(1)
      poisson  = dmat(2)
      shearmod = elastmod / (2.*(1.+poisson))
      density = 1.0
      reftemp = 0.0
      do i=1, 96
         mat(i) = 0
      enddo
      if (analysis.ne.1) then
         write(iws,'(t1,a,a)') 'PATRAN Neutral Format: ',
     &  'Only Isotropic Linear Elastic Material Behaviour is generated'
      endif 
      mat(1)  = density
      mat(2)  = reftemp
      mat(27) = elastmod
      mat(28) = elastmod
      mat(29) = elastmod
      mat(30) = poisson
      mat(31) = poisson
      mat(32) = poisson
      mat(33) = shearmod
      mat(34) = shearmod
      mat(35) = shearmod
      write(io,101) it, id,iv,kc,0,0,0,0,0
      write(io,102) (mat(i),i=1, 96)
101   format( i2,8i8)
102   format( 5e16.9 )
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine patran_04(io)
c
c PACKET Type 04: Element Properties
c 
c  Header: IT,ID,IV,KC, N1,N2,N3,N4,0
c
c   IT = 4
c   ID = Property ID
c   IV = Material ID
c   KC = Card count (number of data cards after the header)
c   N1* = Shape
c   N2* = Nodes
c   N3* = Configuration
c   N4* = Number of data fields
c    (* not needed in the current application)
c  Card 1: DATA
c
c   DATA = Property data for defined type as required by the
c          analysis program. (1 to 5 property fields per record
c          in 16 character fields.)
c
c modification in this subroutine to be able to connect the
c property of the element with the material (elastic-linear)
c   by christophe
c
      implicit none
      integer          io,it,id,iv,kc
      double precision prop(15) ! limited by 15 properties
      it = 4
      id = 1
      iv = 1
      kc = 1    ! put kc=1 and not 0
c you have to write shape and node in this case shape is define
c  in element packet =8
c  node is the number of node per element =8
c  this value below =1 number of data field
c  it seems we need to put this value equal 1
c  to connect the property to the material and write prop(1)=1
c
      prop(1)=1.0
      write(io,101) it,id,iv,kc,8,8,0,1,0
      write(io,102) prop(1) 
c      write(io,102) (prop(i),i=1, 5)
 101  format( i2,8i8)
 102  format( 5e16.9 )
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine patran_06(io,no_eset,eset1,eset2,press)
c
c PACKET Type 06: Distributed Loads
c
c  Header: IT,ID,IV,KC, 0,0,0,0,0
c
c   IT = 6
c   ID = Element ID (=element number)
c   IV = Load set ID  (use IV = 2, since IV=1 for nodal loads below)
c   KC = Card count (number of data cards after the header)
c
c  Card 1: LTYPE, EFLAG, CFLAG, ICOMP(6), NODE(8), NFE
c
c   LTYPE = Load Type (0=line or 1=surface)
c   EFLAG = Element Flag (0 or 1)
c   CFLAG = Node Flag (0 or 1)
c   ICOMP = 6 load component flags (0 or 1)
c   NODE  = 8 element node flags (0 or 1)
c   NFE   = Edge number (1-12) or FACE number (1-6)
c
c  Card 2: PDATA
c
c   PDATA  = non-x=zero load components (defined by the number
c           of non-zero components in ICOMP(i) and NODE(i) )
c
      ! modification to be able to read the pressure value by christophe
      implicit none
      integer  io,no_eset,eset1(2,*),eset2(2,*),i,j
      double precision press
      integer  it,id,iv,kc,nfe,ltype,eflag,cflag,
     &         ne,face_no(6),icomp(8),node(8),fnode(4,6)
c
c  Default face numbering according to WARP3D - different from PATRAN
c
c  WARP3D                      PATRAN
c  Face #  Local Vertex Nodes  Face #
c     1       1-4-3-2            3
c     2       5-6-7-8            5
c     3       1-2-6-5            1
c     4       3-4-8-7	           2
c     5       2-3-7-6            4
c     6       1-5-8-4            6
c
      data (face_no(i),i=1,6) / 3, 5, 1, 2, 4, 6 /
c
      data  (fnode(j,1),j=1,4)/ 1, 4, 3, 2 /
      data  (fnode(j,2),j=1,4)/ 5, 6, 7, 8 /
      data  (fnode(j,3),j=1,4)/ 1, 2, 6, 5 /
      data  (fnode(j,4),j=1,4)/ 3, 4, 8, 7 /
      data  (fnode(j,5),j=1,4)/ 2, 3, 7, 6 /
      data  (fnode(j,6),j=1,4)/ 1, 5, 8, 4 /
c
c      data (icomp(i),i=1, 6)  /1,0,0,0,0,0/
c      data (node(i),i=1, 8)   /0,0,0,0,0,0,0,0/
c
      it = 6
      iv = 2
      kc = 2
c
      ne = eset1(1,1)
      do i=1, ne
         do j=1, 8
            icomp(j) = 0
            node(j)  = 0
         enddo
         id = eset1(1,i+1)
         nfe = face_no( eset1(2,i+1) )
         do j=1, 4
            node( fnode( j, eset1(2,i+1) ) ) = 1
         enddo
         !icomp(nfe) = 1
         icomp(1)=1
         write(io,101) it, id,iv,kc, 0,0,0,0,0
         ltype = 1
         eflag = 1
         cflag = 0
         write(io,102) ltype,eflag,cflag,(icomp(j),j=1,6),
     &                  (node(j),j=1,8), 0	 ! modif nfe to 0
         write(io,103) press
      enddo
c
      if (no_eset.eq.1) return
c
      ne = eset2(1,1)
      do i=1, ne
         do j=1, 8
            icomp(j) = 0
            node(j)  = 0
         enddo
         id = eset2(1,i+1)
         nfe = face_no( eset2(2,i+1) )
         do j=1, 4
            node( fnode( j, eset2(2,i+1) ) ) = 1
         enddo
c           icomp(nfe) = 1
         icomp(1)=1
         write(io,101) it, id,iv,kc, 0,0,0,0,0
         ltype = 1
         eflag = 1
         cflag = 0
         write(io,102) ltype,eflag,cflag,(icomp(j),j=1,6),
     &                  (node(j),j=1,8), 0	! modif nfe to 0
         write(io,103) press
      enddo
c
101   format(i2,8i8)
102   format(i1,i1,i1,6i1,8i1,i2)
103   format(5e16.9)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine patran_07(io,no_of_nodes,idbc,rdbc)
c
c PACKET Type 07: Node Forces
c
c  Header: IT,ID,IV,KC, 0,0,0,0,0
c
c   IT = 7
c   ID = Node ID (=node number)
c   IV = Load set ID
c   KC = Card count (number of data cards after the header)
c
c  Card 1: CID, ICOMP
c
c   CID Coordinate frame ID (use default Frame =0)
c   ICOMP(i) = 6 force component flags (0 or 1)
c
c  Card 2: FDATA
c
c   FDATA  = non-x=zero load components
c   (defined by the number of non-zero components in ICOMP(i))
c
      implicit none
      integer          io,no_of_nodes,idbc(3,*)
      double precision rdbc(3,*),fdata(3)
      integer          it,id,iv,kc,cid,icomp(6),i,j,k
      logical          ok
      it = 7
      iv = 1
      kc = 2
      cid = 0
      do i=1, no_of_nodes
         id = i
         icomp(1) = 0
         icomp(2) = 0
         icomp(3) = 0
         k = 0
         ok = .false.
         if ( (idbc(1,i).eq.2) .or. (idbc(1,i).eq.3) ) then
            k = k + 1
            ok = .true.
            icomp(1) = 1
            fdata(k) = rdbc(1,i)
         endif
         if ( (idbc(2,i).eq.2) .or. (idbc(2,i).eq.3) ) then
            k = k + 1
            ok = .true.
            icomp(2) = 1
            fdata(k) = rdbc(2,i)
         endif
         if ( (idbc(3,i).eq.2) .or. (idbc(3,i).eq.3) ) then
            k = k + 1
            ok = .true.
            icomp(3) = 1
            fdata(k) = rdbc(3,i)
         endif
         if (ok) then
            write(io,101) it, id,iv,kc, 0,0,0,0,0
            write(io,102) cid, (icomp(j),j=1, 3),  0,0,0
            write(io,103) ( fdata(j), j=1, k )
         endif
      enddo
101   format(i2,8i8)
102   format(i8,6i1)
103   format(5e16.9)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine patran_08(io,no_of_nodes,idbc,rdbc)
c
c PACKET Type 08: Node Displacements
c
c Header: IT,ID,IV,KC, 0,0,0,0,0
c
c   IT = 8
c   ID = Node ID (=node number)
c   IV = Constraint set ID
c   KC = Card count (number of data cards after the header)
c
c Card 1: CID, ICOMP
c
c   CID Coordinate frame ID
c   ICOMP(i) = 6 displacement component flags (0 or 1)
c
c Card 2: FDATA
c
c   FDATA  = Non-blank displacement components (may be 0.0)
c   (# is defined by the number of non-zero components in ICOMP(i))
c
      implicit none
      integer          io,no_of_nodes,idbc(3,*)
      double precision rdbc(3,*),fdata(3)
      integer          it,id,iv,kc,cid,icomp(6),i,j,k
      logical          ok
      it = 8
      iv = 1
      kc = 2
      cid = 0
      do i=1, no_of_nodes
         id = i
         icomp(1) = 0
         icomp(2) = 0
         icomp(3) = 0
         k = 0
         ok = .false.
         if (idbc(1,i).eq.1) then
            k = k + 1
            ok = .true.
            icomp(1) = 1
            fdata(k) = rdbc(1,i)
         endif
         if (idbc(2,i).eq.1) then
            k = k + 1
            ok = .true.
            icomp(2) = 1
            fdata(k) = rdbc(2,i)
         endif
         if (idbc(3,i).eq.1) then
            k = k + 1
            ok = .true.
            icomp(3) = 1
            fdata(k) = rdbc(3,i)
         endif
         if (ok) then
            write(io,101) it, id,iv,kc, 0,0,0,0,0
            write(io,102) cid, (icomp(j),j=1, 3),  0,0,0
            write(io,103) ( fdata(j), j=1, k )
         endif
      enddo
101   format(i2,8i8)
102   format(i8,6i1)
103   format(5e16.9)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
