c----67--1---------2---------3---------4---------5---------6---------712
c
	program mesh3d_scpcell
c
c----67--1---------2---------3---------4---------5---------6---------712
c
C  The program creates a 3-dimensional finite element model of
C  a specimen with a semi-elliptical surface crack. The cracktip
C  is modelled as a semi-cirkular notch or as being initially
C  sharp. Nodes and elements are created the finite element
C  ABAQUS and WARP3D.
C
C    The following files is created :
C
C    =>  job.inp/in    ( Input information ABAQUS/WARP3D)
C    =>  job.##        ( node, elemeent, etc Input information )
C    =>  job_mesh.sta  ( Model-statistics  )
C    =>  job_set.sta   ( Node & element Info for postprocessing
C                        purposes )
C
C    =>  job_1.dat     ( Coordinates for postscript graphics)
C    =>  job_1.inf     ( plot-information )
C        . . . .         . . . . . . . . .
C    =>  job_5.dat
C    =>  job_5.inf
C
C     ( The plot files are optionel, will be created if PL=1 )
C
C          VERSION Oct. 18, 1993 (modified Nov, 1996)
C                        by
C                  JONAS  FALESKOG
C
c----67--1---------2---------3---------4---------5---------6---------712
c
      implicit none
c
c Array bounds   Zone C: nodc(i,j,k) ncm = nc1*nc2*nc3
c                Zone S: nods(i,j,k) nsm = ns1*ns2*ns3
c                Zone A: noda(i,j,k) nam = na1*na2*na3
c                Zone B: nodb(i,j,k) nbm = nb1*nb2*nb3
c
      include 'common_eln.f'
c
      include 'common_nod.f'
c
      integer    ncm,nsm,nam,nbm, nesm,nesm2
      parameter (ncm=160000,nsm=220000,nam=100000,nbm=60000)
      parameter (nesm=240, nesm2=nesm*nesm)
c
c....INTEGER VARIABLES
      integer  nc1,nc2,nc3, ns1,ns2,ns3, na1,na2,na3, nb1,nb2,nb3
c
      integer  nea2,neb2, estk_gc(2,nesm),estk_c(2,nesm),
     1         estk_s(2,nesm), estk_a(2,nesm2), estk_b(2,nesm2),
     2         nstk_c(3,401),nstk_s(3,401),nstk_a(3,201),nstk_b(3,201)
c
      integer  nodc(ncm),nods(nsm),noda(nam),nodb(nbm),
     1         pl,jobh,etyp,elnum,elnum_cell,no_of_nodes,
     2         mbtype, ncell,necb
c
c..REAL VARIABLES
      double precision   y(0:200),z(0:200,2),mb_bias
c
      character job*40,prog*20,datum*20
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer      imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer         kstart,kstep
      common /nblock/ kstart,kstep
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
      double precision t,w,c,a,alfa
      common /geom/    t,w,c,a,alfa
c
c..Input is read from file mesh3d_scp.in
c
      call indata_read(job,jobh,prog,datum,etyp,y,z,
     1                 mbtype,mb_bias,pl,
     2                 nodc,nc1,nc2,nc3,ncm,  nods,ns1,ns2,ns3,nsm,
     3                 noda,na1,na2,na3,nam,  nodb,nb1,nb2,nb3,nbm,
     4                 nesm,nesm2,nea2,neb2,estk_gc,estk_c,estk_s,
     5                 estk_a,estk_b, ncell,necb)
c
c..Generate node numbers and store those in the arrays NODS,NODA & NODB
      call nodnumber(nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3,
     &               noda,na1,na2,na3, nodb,nb1,nb2,nb3 )
c
c..Generate node coordinates
      call node_coordinates(nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3,
     1                      noda,na1,na2,na3, nodb,nb1,nb2,nb3,
     2                      y,z,job,jobh,pl,etyp,necb,mbtype,mb_bias)
c
c..Generate elements
      call make_element(nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3, 
     1      noda,na1,na2,na3, nodb,nb1,nb2,nb3, elnum,elnum_cell,
     2      nea2,neb2,estk_gc,estk_c,estk_s,estk_a,estk_b)
c
c..Sort out the nodes in the model
      call sort_out_nodes(job,jobh, elnum, no_of_nodes, etyp,prog, 
     1       nstk_c,nstk_s,nstk_a,nstk_b, nodc,nc1,nc2,nc3, 
     2     nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c..Write out node coord. and connectivity on files
      call write_out_nod_el(job,jobh,elnum,elnum_cell,etyp,prog)
c
c..SKAPA INDATA FILEN
c..Generate the input decks to the FEM-program in question
      if (prog.eq.'ABAQUS') then
         call abaqus_inp(job,jobh,etyp,  nea2,neb2, 
     1               estk_s,estk_a,estk_b,  nodc,nc1,nc2,nc3,
     2        nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3 )
      endif
c
      if ((prog.eq.'WARP3D') .or. (prog.eq.'PATRAN')) then
         call warp3d_inp(job,jobh,etyp,elnum_cell,elnum,no_of_nodes,
     1       estk_gc,estk_c,estk_s,nea2,estk_a, nodc,nc1,nc2,nc3,
     2       nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      endif
c
      if (prog.eq.'PATRAN') then
         call patran_neu(job,jobh,etyp,elnum_cell,elnum,
     1        no_of_nodes, nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3, 
     2        noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      endif
c
c..Generate bookkeeping files
c
      call mesh_statestik(job,jobh,datum,y,z,etyp,mbtype,mb_bias,ncell,
     &     elnum,elnum_cell,no_of_nodes,prog,  nea2,neb2, estk_gc,
     &     estk_c,estk_s,estk_a,estk_b, nstk_c,nstk_s,nstk_a,nstk_b)
c
      stop
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine indata_read(job,jobh,prog,datum,etyp,y,z,
     1           mbtype,mb_bias,pl,
     2           nodc,nc1,nc2,nc3,ncm,  nods,ns1,ns2,ns3,nsm,
     3           noda,na1,na2,na3,nam,  nodb,nb1,nb2,nb3,nbm,
     4           nesm,nesm2,nea2,neb2,estk_gc,estk_c,estk_s,
     5           estk_a,estk_b, ncell,necb )
c
c---  The subroutine reads indata parameters from the file
c---- mesh3d_nsct.in. Some of the parameters are checked so that
c---- they are located in a reasonable interval.
c
      implicit none
c
      external          ellipse_l
      double precision  ellipse_l
c
      include 'common_nod.f'
c
      double precision  pi
c      parameter ( pi = 4.*atan(1.))
      parameter ( pi = 3.141592654)
c
      integer nc1,nc2,nc3,ncm,nodc(ncm), ns1,ns2,ns3,nsm,nods(nsm),
     &        na1,na2,na3,nam,noda(nam), nb1,nb2,nb3,nbm,nodb(nbm)
      integer nesm,nesm2,nes2,nea2,neb2,estk_gc(2,nesm),estk_c(2,nesm),
     &        estk_s(2,nesm), estk_a(2,nesm2),estk_b(2,nesm2)
      integer i,j,k,num,pl,jobh,etyp,mbtype, max_node_no, io
      integer ncell,necb,ncf,ncz
c
      double precision  y(0:200),z(0:200,2),rvek(15),mb_bias,x,
     &                  rho0,phi1,phi2,cl
c
      character ch*3,job*40,prog*20,row*80,datum*20,infile*80
c
      logical   ok
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer      imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      integer       ncd1,ncd2,ncd3
      common /cell/ ncd1,ncd2,ncd3
c
      double precision t,w,c,a,alfa
      common    /geom/ t,w,c,a,alfa
      double precision      dcell,eta_n,eta_t1,eta_t2
      common    /geom_cell/ dcell,eta_n,eta_t1,eta_t2
      double precision      p_lcx1,p_lcx2,p_lcy2,p_alfac
      common /geom_zone_s/  p_lcx1,p_lcx2,p_lcy2,p_alfac
c
      double precision  a1,b1
      common /ellipse/  a1,b1
c
      write(*,8001) ' '
      write(*,8001) 
     & '>> Surface crack mesh generator (computational cells)'
      write(*,8001) '  '
      write(*,'(t1,a,$)') '>> Input file name: '
      read(*,'(a)') infile
      inquire(file=infile,exist=ok)
      if (.not.ok) then
            write(*,'(t4,3a)')'>>The file: ',infile,'doesn''t exist!'
            write(*,*) '>> Program aborted...'
         stop
      endif
      write(*,'(/t5,a,a/)') '* reading in-put-data from the file:',
     &      infile
      io = 25
      open(unit=25,file=infile,status='old')
 5    read(io,'(a)') ch
      if ((ch.ne.'*in').and.(ch.ne.'*IN')) goto 5
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
      jobh = j-i+1
      job(1:jobh) = row(i:j)
      do i=jobh+1, 40
         job(i:i)=' '
      enddo
c
c...Program
      read(io,'(a)') row
      i = 0
 20   continue
         i = i + 1
      if (row(i:i).eq.' ') goto 20
         j = i
 25   continue
         j = j + 1
      if (row(j:j).ne.' ') goto 25
      j = j - 1
      do k=1, 20
         prog(k:k) = ' '
      enddo
      prog(1:(j-i+1)) = row(i:j)
c
c... Date
      read(io,'(a)') row
      i=0
 26   continue
         i = i + 1
      if (row(i:i).eq.' ') goto 26
      datum = row(i:(i+19))
c
c... ncell, dcell, eta_n, necb:
c
        read(io,'(a)') row
	call char_to_real(row,rvek,num)
	if (num.ne.4) goto 99
	ncell = int(rvek(1)+0.5)
	dcell = rvek(2)
	eta_n = rvek(3)
	necb = int(rvek(4)+0.5)
c
c... eta_t1, eta_t2, ncd1, ncd2:
c
        read(io,'(a)') row
	call char_to_real(row,rvek,num)
	if (num.ne.4) goto 99
	eta_t1 = rvek(1)
	eta_t2 = rvek(2)
	ncd1 = int(rvek(3)+0.5)
	ncd2 = int(rvek(4)+0.5)
c
c... p_lcx1, p_lcx2, p_lcy2, p_alfac:
c
	read(io,'(a)') row
	call char_to_real(row,rvek,num)
	if (num.ne.3) goto 99
c	if (num.ne.4) goto 99
	p_lcx1 = rvek(1)
	p_lcx2 = rvek(2)
c	p_lcy2 = rvek(3)
	p_alfac = rvek(3)
c
c... mr, sfred, sfred_type, mv, sjred_type
c
        read(io,'(a)') row
        call char_to_real(row,rvek,num)
        if (num.ne.5) goto 99
        mr = int(rvek(1)+0.5)
        sfred = int(rvek(2)+0.5)
        sfred_type = int(rvek(3)+0.5)
        mv = int(rvek(4)+0.5)
        sjred_type  = int(rvek(5)+0.5)
c
c... m1, m2:  (nb is determined by ncd1, ncd2, ncd3, sjred_type)
c
        read(io,'(a)') row
        call char_to_real(row,rvek,num)
        if (num.ne.2) goto 99
        m1 = int(rvek(1)+0.5)
        m2 = int(rvek(2)+0.5)
c
c... mb, mbtype, mb_bias:
c
        read(io,'(a)') row
        call char_to_real(row,rvek,num)
        if (num.ne.4) goto 99
        mb = int(rvek(1)+0.5)
        nb = int(rvek(2)+0.5)
        mbtype = int(rvek(3)+0.5)
        mb_bias = rvek(4)
c
c... lt, lred, rtype:
c
        read(io,'(a)') row
        call char_to_real(row,rvek,num)
        if (num.ne.3) goto 99
        lt = int(rvek(1)+0.5)
        lred = int(rvek(2)+0.5)
        rtype = int(rvek(3)+0.5)
c
c... etyp:
c
        read(io,'(a)') row
        call char_to_real(row,rvek,num)
        if (num.ne.1) goto 99
        etyp = int(rvek(1)+0.5)
c
c... t, w, c, a:
c
c        read(io,*) t,w,c,a
        read(io,'(a)') row
	call char_to_real(row,rvek,num)
	if (num.ne.4) goto 99
	t=rvek(1)
	w=rvek(2)
	c=rvek(3)
	a=rvek(4)
c
c        write(*,'(t1,g30.20)') ' t = ',t
c        write(*,'(t1,g30.20)') ' w = ',w
c        write(*,'(t1,g30.20)') ' c = ',c
c        write(*,'(t1,g30.20)') ' a = ',a
c
c... pl plotting parameter pl=1 allways (modification 941005):
c
	pl=1
c
c... Yi, Z1i, Z2i :
c
	do i=0, lt
           read(io,'(a)') row
	   call char_to_real(row,rvek,num)
	   if (num.ne.3) goto 99
	   y(i)=rvek(1)
	   z(i,1)=rvek(2)
	   z(i,2)=rvek(3)
	   row='***empty**'
	enddo
	close(io)
c
c Determine model parameters:
c
      ncf = (ncell+2)/2+1 + 6
      if (sfred_type.eq.1) then
         if (mod(ncf,2).ne.0) then
            ncf = ncf + ( 2 - mod(ncf,2) )
            ncell = 2*(ncf-7) - 2
         endif
         mh = ( ncf - 2*mv ) / 2
      elseif (sfred_type.eq.2) then
         if (mod(ncf,4).ne.0) then
            ncf = ncf + ( 4 - mod(ncf,4) )
            ncell = 2*(ncf-7) - 2
         endif
         mh = ( ncf/2 - 2*mv ) / 2
      elseif (sfred_type.eq.3) then
         if (mod(ncf,6).ne.0) then
            ncf = ncf + ( 6 - mod(ncf,6) )
            ncell = 2*(ncf-7) - 2
         endif
         mh = ( ncf/3 - 2*mv ) / 2
      else
         write(*,*) '>> sfred_type must be equal to 1 or 2 or 3'
         stop
      endif
c
c   Angular grading.
c
      alfa = sqrt(c*c-a*a)/2.
      rho0 = dsqrt( (c+a) / (c-a) )
      a1 = alfa * ( rho0 + 1.0/rho0 )
      b1 = alfa * ( rho0 - 1.0/rho0 )
      phi1 = 0.0
      phi2 = pi/2.
      call qromb(ellipse_l,phi1,phi2,cl, 1.0d-6)
c
c   Adjust ncd3 and eta_t2
c      ncd3 = ( cl - dcell*( ncd1-ncd2*(1.+eta_t)*0.5 ) ) /
c     &                   (dcell*eta_t)
      ncd3 = ( cl/dcell - (ncd1*eta_t1 + ncd2*(eta_t1+eta_t2)*0.5) ) /
     &                   eta_t2
      ncz = ncd1+ncd2+ncd3
      if (sjred_type.eq.2) then
         if (mod(ncz,4).ne.0) then
            ncd3 = ncd3 + ( 4 - mod(ncz,4) )
            ncz = ncd1+ncd2+ncd3
         endif
         na = ncz / 2
      elseif (sjred_type.eq.3) then
         if (mod(ncz,6).ne.0) then
            ncd3 = ncd3 + ( 6 - mod(ncz,6) )
            ncz = ncd1+ncd2+ncd3
         endif
         na = ncz / 3
      else
         na = ncz
      endif
      eta_t2 = ( cl/dcell - eta_t1*(ncd1+0.5*ncd2) ) / (0.5*ncd2+ncd3)
c
      write(*,'(t10,a,i4,a/)')'>>',na,' elements along the crack front'
c
c  Indicies in node arrays
c
c  Zone C:
      imc = 2*ncell + 5 + 2
      jmc = 2*ncz + 1
      kmc = 9
c
c  Zone S:
      ims = 2*ncf+1
      jms = jmc
      if (sjred_type.eq.1) then
         kms = 2*mr + 1
      else
         kms = 2*(mr+1) + 1
      endif
      ksr1 = 2*sfred+1
c
c  Zone A:
      ma = m1 + m2 + mh + mh
      ima = 2*ma+1
      jma = 2*na+1
c
c  Zone B:
      imb=2*mb+1
      jmb=2*nb+1
c
c  If rtype=1, one extra element row is included. If rtype=2, two
c    extra element rows is included.
c
      if (rtype.eq.0) then
         kma=2*lt+1
      elseif (rtype.ge.1) then
         lt=lt+1
         kma=2*lt+1
         do i=lt, lred+1, -1
            y(i)=y(i-1)
            z(i,1)=z(i-1,1)
            z(i,2)=z(i-1,2)
         enddo
         y(lred)=(y(lred-1)+y(lred+1))/2.d0
         z(lred,1) =(z(lred-1,1)+z(lred+1,1))/2.d0
         z(lred,2) =(z(lred-1,2)+z(lred+1,2))/2.d0
      endif
      if (rtype.eq.2) then
         lt=lt+1
         kma=2*lt+1
         do i=lt,lred+3, -1
            y(i)=y(i-1)
            z(i,1)=z(i-1,1)
            z(i,2)=z(i-1,2)
         enddo
         y(lred+2)=(y(lred+1)+y(lred+3))/2.d0
         z(lred+2,1) =(z(lred+1,1)+z(lred+3,1))/2.d0
         z(lred+2,2) =(z(lred+1,2)+z(lred+3,2))/2.d0
      endif
      ksr1=2*sfred+1
      kar1=2*lred+1
      kar2=2*(lred+2)+1
c
c  Set matrix bounds
c
      nc1 = imc
      nc2 = jmc
      nc3 = kmc
      if ((nc1*nc2*nc3).gt.ncm) then
         write(*,'(t1,a)') '>> increase the size of array nodc(i,j,k)'
         write(*,'(t4,a,i8)') 'current  size = ', ncm
         write(*,'(t4,a,i8)') 'required size = ', nc1*nc2*nc3
         stop
      endif
      do i=1, ncm
         nodc(i) = 0
      enddo
c
      ns1 = ims
      ns2 = jms
      ns3 = kms
      if ((ns1*ns2*ns3).gt.nsm) then
         write(*,'(t1,a)') '>> increase the size of array nods(i,j,k)'
         write(*,'(t4,a,i8)') 'current  size = ', nsm
         write(*,'(t4,a,i8)') 'required size = ', ns1*ns2*ns3
         stop
      endif
      do i=1, nsm
         nods(i) = 0
      enddo
c
      na1 = ima
      na2 = jma
      na3 = kma
      if ((na1*na2*na3).gt.nam) then
         write(*,'(t1,a)') '>> increase the size of array noda(i,j,k)'
         write(*,'(t4,a,i8)') 'current  size = ', nam
         write(*,'(t4,a,i8)') 'required size = ', na1*na2*na3
         stop
      endif
      do i=1, nam
         noda(i) = 0
      enddo
c
      nb1 = imb
      nb2 = jmb
      nb3 = kma
      if ((nb1*nb2*nb3).gt.nbm) then
         write(*,'(t1,a)') '>> increase the size of array nodb(i,j,k)'
         write(*,'(t4,a,i8)') 'current  size = ', nbm
         write(*,'(t4,a,i8)') 'required size = ', nb1*nb2*nb3
         stop
      endif
      do i=1, nbm
         nodb(i) = 0
      enddo
c
c  Check array bounds of npos (coordinate vector)
c
      max_node_no = jmc*imc*kmc + jms*ims*kms + (jma*ima+jmb*imb)*kma
      if (max_node_no.gt.inm) then
         write(*,'(t1,a)') '>> increase the size of array npos(inm,3)'
         write(*,'(t4,a,i8)') 'current  size = ', inm
         write(*,'(t4,a,i8)') 'required size = ', max_node_no
         stop
      endif
c
      nes2 = (kms+1)/2
      if (nes2.gt.nesm) then
         write(*,'(t1,a)') '>> increase the array size parameter nesm'
         write(*,'(t4,a,i8)') 'current  size = ', nesm
         write(*,'(t4,a,i8)') 'required size = ', nes2
      endif
      nea2 = (kma+1)/2
      if ( (nea2*(jma+1)/2) .gt. (nesm*nesm) ) then
         write(*,'(t1,a)') '>> increase the array size parameter nesm'
         write(*,'(t4,a,i8)') 'current  size = ', nesm
         x = real(nea2*(jma+1)/2)
         write(*,'(t4,a,i8)') 'required size = ', int( sqrt(x) )+1
      endif
      neb2 = (kma+1)/2
      if ( (neb2*(imb+1)/2) .gt. (nesm*nesm) ) then
         write(*,'(t1,a)') '>> increase the array size parameter nesm'
         write(*,'(t4,a,i8)') 'current  size = ', nesm
         x = real(neb2*(imb+1)/2)
         write(*,'(t4,a,i8)') 'required size = ', int( sqrt(x) )+1
      endif
      if ( (neb2) .gt. (nesm) ) then
         write(*,'(t1,a)') '>> increase the array size parameter nesm'
         write(*,'(t4,a,i8)') 'current  size = ', nesm
         write(*,'(t4,a,i8)') 'required size = ', neb2
      endif
      do i=1, nesm
         do j=1,2
            estk_gc(j,i) = 0
            estk_c(j,i) = 0
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
c TEST
c
      write(*,'(5(a,i5))') ' ims=',ims,' kms=',kms,' mv=',mv,' mh=',mh
c
      call indata_test(prog,mr,m1,m2,nb,na,rtype,sfred,sfred_type,
     &   sjred_type,etyp,lt,p_lcx1,p_lcx2,eta_n,eta_t2,a,c,t,y,z)
c
      write(*,'(t10,a,a/)') '=> Generating in-put-data for ',prog
c
      return
99    write(*,'(t10,a/t10,a)')
     1  '* the number of data given does not match the number of',
     2  '  data required at each line, => correct the in-put file!!'
 8001 format(t1,a)
c
      stop
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine indata_test(prog,mr,m1,m2,nb,na,rtype,sfred,sfred_type,
     1           sjred_type,etyp,lt,p_lcx1,p_lcx2,eta_n,eta_t2,
     2           a,c,t,y,z )
c
c                         Check of in-put-data
c
      implicit none
c
      integer  mr,m1,m2,nb,na,rtype,sfred,sfred_type,sjred_type,
     &         etyp,lt,i
c
      double precision  p_lcx1,p_lcx2,eta_n,eta_t2,a,c,t,
     &                  y(0:200),z(0:200,2)
c
      character prog*20
c
      if (m1.lt.2) then
         write(*,'(/t4,a)') 'ERROR: m1 >= 2 not ok!' 
         goto 99
      endif
      if (m2.lt.2) then
         write(*,'(/t4,a)') 'ERROR: m2 >= 2 not ok!' 
         goto 99
      endif
      if ( (rtype.ne.0).and.(rtype.ne.1).and.(rtype.ne.2) ) then
         write(*,'(t4,a)') 'ERROR: rtype must be equal to 0, 1 or 2'
         goto 99
      endif
      if ( (rtype.eq.2).and.(mod(nb,2).ne.0) ) then
         write(*,'(t4,a)') 'ERROR: Elem. reduction; nb must be even'
         goto 99
      endif
      if ( (rtype.eq.2).and.(mod(na,2).ne.0) ) then
         write(*,'(t4,a)')
     &   'ERROR: Elem. reduction; na must be even'
         goto 99
      endif
      if ( ( rtype.gt.2 ).and.( mod((m1+m2),2).ne.0 ) ) then
         write(*,'(t4,a)') 'ERROR: if rtype > 0; m1+m2 must be even!'
         goto 99
      endif
      if (p_lcx1.gt.p_lcx2) then
         write(*,'(/t4,a)') 'ERROR: p_lcx2 > p_lcx1 not fullfilled' 
         goto 99
      endif
      if ( (sfred_type.ne.1).and.(sfred_type.ne.2).and.
     &     (sfred_type.ne.3) ) then
         write(*,'(/t4,a)')
     &              'ERROR: sfred_type must be equal to 1 or 2 or 3'
         goto 99
      endif
      if (sfred.ge.(mr-1)) then
         write(*,'(/t4,a)') 'ERROR: sfred < mr-1; not fullfilled! '
         goto 99
      endif
      if ( (sjred_type.ne.1).and.(sjred_type.ne.2).and.
     &     (sjred_type.ne.3) ) then
         write(*,'(/t4,a)')'ERROR: sfred_type must be equal to 1 2 or 3'
         goto 99
      endif
      if ((etyp.ne.8).and.(etyp.ne.20).and.(etyp.ne.27)) then
         write(*,'(/t4,a)') 'ERROR: etyp must be 8, 20 or 27'
         goto 99
      endif
      if ((prog.eq.'WARP3D').and.(etyp.ne.8)) then
         write(*,'(t4,a)')'if warp3d-input 8-node element must be used'
         goto 99
      endif
      if ( (prog.ne.'ABAQUS').and.(prog.ne.'WARP3D').and.
     &     (prog.ne.'PATRAN') ) then
         write(*,'(t1,a,a)') 'program name must be either',
     &                       ' ABAQUS, WARP3D or PATRAN.'
         stop
      endif
      if (a.gt.c) then
         write(*,'(/t4,a)')  ' c > a , not fulfilled !'
         goto 99
      endif
      if (a.gt.t) then
	 write(*,'(/t10,a)') ' t > a , not fulfilled !'
	 goto 99
      endif
      do i=1, lt
         if ( ( z(i,2).lt.z(i,1) ).or.( y(i).lt.y(i-1) ) ) then
            write(*,'(/t10,a,i2,a,i2,a)')
     &      ' z(' ,i,',2) > z(' ,i, ',1) , not fulfilled!'
            goto 99
         endif
      enddo
c
      if ( (eta_n.le.0.4999).or.(eta_n.ge.1.00001) ) then
         write(*,'(/t4,a,a,g15.7,a)') 'WARNING: ',
     &   ' eta_n should be between 0.5 and 1.0 (eta_n=',eta_n,')' 
      endif
c      if (eta_t2.lt.1.0) then
c         write(*,'(/t4,a,g15.7,a)')
c     &   'WARNING: eta_t2 should be  >  1.0 (eta_t2=',eta_t2,')'
c      endif
c
      return
 99   continue
      write(*,'(/t10,a/t12,a)')
     &     '* the modellparameters is/are incorrect !',
     &     '=> correct the modell.'
      stop
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
	subroutine char_to_real(row,rvek,num)
C	--- The routine transform a character to a real number !
	implicit   none
	integer    num,n,k,i,exponent,es
	double precision  rvek(*),sign,xi,xe
	character  row*80
 
	num=0
	n=1
5	continue
	   i=0
	   k=0
	   exponent=0
	   es=1
	   sign=1.0
10	   continue
C             ! Searcing for the 1:st digit !
	   if ( ( (row(n:n).lt.'0').or.(row(n:n).gt.'9') ).and.
     &          (row(n:n).ne.'.') )  then
	      n=n+1
	      if (n.gt.80) goto 99
              goto 10
	   endif
C              ! Negative sign !
	   if (row (n-1:n-1) .eq.'-') then
	   sign=-1.0
	   end if
20	   continue
C             !  Integer part
	   if ( (row(n:n).ge.'0').and.(row(n:n).le.'9') ) then
	      i=10*i+ichar(row(n:n))-48
	      n=n+1
	      if (n.gt.80) goto 99
	      goto 20
	   endif
 
	   if (row(n:n).ne.'.') goto 40
	   n=n+1
	   if (n.gt.80) goto 99
30	   continue
C             !   Decimal part !
	   if ( (row(n:n).ge.'0').and.(row(n:n).le.'9') ) then
	      i=10*i+(ichar(row(n:n))-48)
	      n=n+1
	      if (n.gt.80) goto 99
	      k=k-1
	      goto 30
	   endif
40	   continue
	   if ( (row(n:n).eq.'E').or.(row(n:n).eq.'e').or.
     &          (row(n:n).eq.'D').or.(row(n:n).eq.'d') ) then
C             ! Exponent
	      n=n+1
	      if (n.gt.80) goto 99
	      if (row(n:n).eq.'+') then
C               ! Positive or negative sign !
	         es=1
	         n=n+1
	         if (n.gt.80) goto 99
	      elseif (row(n:n).eq.'-') then
	         es=-1
	         n=n+1
	         if (n.gt.80) goto 99
	      endif
50	      if ( (row(n:n).ge.'0').and.(row(n:n).le.'9') ) then
	         exponent=10*exponent+ichar(row(n:n))-48
	         n=n+1
	         if (n.gt.80) goto 99
	         goto 50
	      endif
	   endif
           xi = 1000*i
           k = k - 3
c           xi = i
c           k = k
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
	if (num.lt.80) goto 5
99	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine char_to_int(row,ivek,num)
c-- The routine transform a character to an integer number
      integer    ivek(*),num,n,i
      character  row*100
      logical    minus
      num=0
      n=1
5     continue
         i=0
         minus = .false.
10        continue
c... Search for the 1:st digit
         if ( (row(n:n).lt.'0') .or. (row(n:n).gt.'9') ) then
            n=n+1
            if (n.gt.100) goto 99 
            goto 10
         endif
c... Negative sign !
         if (row((n-1):(n-1)).eq.'-') minus=.true.
20       continue
         if ( (row(n:n).ge.'0').and.(row(n:n).le.'9') ) then
            i=10*i+ichar(row(n:n))-48
            n=n+1
            if (n.gt.100) goto 99 
            goto 20
         endif
         num=num+1
         ivek(num)=i
         if (minus) ivek(num)= - ivek(num)
      if (num.lt.100) goto 5
99    return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
