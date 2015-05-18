c----67--1---------2---------3---------4---------5---------6---------7-!
	program mesh3d_scp
c----67--1---------2---------3---------4---------5---------6---------7-!
c
C  The program creates a 3-dimensional finite element model of
C  a specimen with a semi-elliptical surface crack. The cracktip
C  is modelled as a semi-cirkular notch or as being initially
C  sharp. Nodes and elements are created the finite element
C  ABAQUS and ADINA.
C
C    The following files are created :
C
C    =>  job.inp/in    ( Input information ABAQUS/ADINA)
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
C          VERSION Oct. 18, 1993 (modified Jan 12, 1996)
C                        by
C                  JONAS  FALESKOG
C
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      implicit none
c
c Array bounds   Zone S: NODS(i,j,k) NSM = NS1*NS2*NS3
c                Zone A: NODA(i,j,k) NAM = NA1*NA2*NA3
c                Zone B: NODB(i,j,k) NBM = NB1*NB2*NB3
c
      include 'mesh3d_scp_common_eln.f'
      integer  eln_d(0:iem)
      integer  elsym(iem),elund(iem),elsua(iem)
c
      include 'mesh3d_scp_common_nod.f'
c
      integer    nsm, nam,nbm, nesm,nesm2, njintm
      parameter (nsm=200000,nam=200000,nbm=200000)
      parameter (nesm=100, nesm2=nesm*nesm)
      parameter (njintm=10000)
c
c....INTEGER VARIABLES
      integer  ns1,ns2,ns3, na1,na2,na3, nb1,nb2,nb3
c
      integer  nea2,neb2,  estk_s(2,nesm),estk_s0(2,nesm),
     &         estk_a(2,nesm2), estk_b(2,nesm2),estk_b0(2,nesm)
c
      integer  el_jco(2,njintm), nj1,nj2,nj3
c
      integer  nods(nsm), noda(nam), nodb(nbm),
     1         nstk(0:500,7),pl,jobh,etyp,elnum,no_of_nodes,
     2         el_ch_zs,slice, no_of_jcont,cr_front(500,3),zb_xw,
     3         dsl(200),grp, rzero,elast
c COMMON INTEGER
      integer  mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,natyp,rtype,sfred,
     &         ims,jms,kms,is,ima,jma,imb,jmb,kma,
     &         kstart,kstep,ksr1,kar1,kar2
c..REAL*8 VARIABLES
      real*8   y(0:200),z(0:200,2),rs(200),omega,vsh_vec(200,3),rate_b
c COMMON REAL*8
      real*8   t,w,c,a,kappa,alfa,r1,r2,rn,eta,my
 
      character head*45,job*40,prog*20
 
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
      common /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
      common /geom/ t,w,c,a,kappa,alfa,r1,r2,rn,eta,my
      common /nblock/ kstart,kstep
      common /reduce/ ksr1,kar1,kar2
c
c..Input is read from file mesh3d_scp.in
c
      call indata_read(natyp,omega,y,z,etyp,el_ch_zs,
     1     slice,head,job,jobh,pl,prog,zb_xw,rate_b,rzero,elast,
     2     ns1,ns2,ns3,nsm,na1,na2,na3,nam,nb1,nb2,nb3,nbm,
     3     nesm,nesm2,nea2,neb2,estk_s,estk_s0,estk_a,estk_b,estk_b0,
     4     njintm, nj1,nj2,nj3 )
c
c..Generate node numbers and store those in the arrays NODS,NODA & NODB
      call nodnumber(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &               nodb,nb1,nb2,nb3, prog,rzero,elast)
c
c..Generate node coordinates
      call nodposition(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     1                 nodb,nb1,nb2,nb3,
     2     natyp,y,z,omega,job,jobh,pl,rs,etyp,zb_xw,rate_b,rzero)
c
c..Generate elements
      call make_element(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     1     nodb,nb1,nb2,nb3, elnum,nea2,neb2,estk_s,estk_a,estk_b)
c
c..Sort out the nodes in the model
      call sort_out_nodes(elnum,no_of_nodes,etyp,prog,
     1     el_ch_zs,slice, eln_d,nstk,dsl,grp,
     2     nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3,
     3     nea2,neb2, estk_s,estk_s0,estk_a,estk_b,estk_b0)
c
c..Write out node coord. and connectivity on files
      call write_out_nod_el(job,jobh,elnum,etyp,eln_d,
     &                      el_ch_zs,slice,prog)
c
c..SKAPA INDATA FILEN
c..Generate the input decks to the FEM-program in question
      if (prog.eq.'ABAQUS') then
	 call abaqus_inp(job,jobh,head,etyp,
     1        nea2,neb2,estk_s,estk_s0,estk_a,estk_b,estk_b0,
     2        el_ch_zs,no_of_jcont,elast,
     2        nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      endif
c
      if (prog.eq.'ADINA') then
c..Generate J-contours for the VCE-method (only in case of ADINA)
         call create_j_contours(estk_a,nea2,el_ch_zs,nods,ns1,ns2,ns3,
     &      eln_d,no_of_jcont,vsh_vec,cr_front,el_jco,nj1,nj2,nj3)
	 call adina_in(job,jobh,el_ch_zs,slice,el_jco,nj1,nj2,nj3,
     1        no_of_jcont,vsh_vec,cr_front,dsl,grp,etyp,elast,
     2        nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      endif
c
      if ((prog.eq.'WARP3D') .or. (prog.eq.'PATRAN')) then
         call warp3d_inp(job,jobh,etyp,elast,elnum,no_of_nodes,
     1        nea2,neb2, estk_s,estk_s0, estk_a, estk_b,estk_b0,
     2        nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      endif
c
      if (prog.eq.'PATRAN') then
         call patran_neu(job,jobh,etyp,elast,elnum,no_of_nodes,
     1        nea2,neb2, estk_s,estk_s0, estk_a, estk_b,estk_b0,
     2        nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      endif
c
c..Generate bookkeeping files
c
      call mesh_statestik(job,jobh,nstk,y,z,rs,natyp,omega,
     &        etyp,zb_xw,rate_b,head,el_ch_zs,slice,prog,
     &        nea2,neb2, estk_s,estk_s0,estk_a,estk_b,estk_b0,
     &        eln_d,elnum,dsl,grp, no_of_jcont )
c
      call model_info(eln_d,no_of_jcont,job,jobh, estk_a,nea2,
     1     nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3 )
c
      write(*,8001) ' '
      write(*,8001) '>> All done...'
      stop
 8001 format(t1,a)
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine indata_read(natyp,omega,y,z,etyp,el_ch_zs,slice,
     1     head,job,jobh,pl,prog,zb_xw,rate_b,rzero,elast,
     2     ns1,ns2,ns3,nsm, na1,na2,na3,nam, nb1,nb2,nb3,nbm,
     3     nesm,nesm2,nea2,neb2,estk_s,estk_s0,estk_a,estk_b,estk_b0,
     4     njintm, nj1,nj2,nj3 )
c
C---  The subroutine reads indata parameters from the file
C---- mesh3d_nsct.in. Some of the parameters are checked so that
C---- they are located in a reasonable interval.
c
      implicit none
c
      integer ns1,ns2,ns3,nsm,na1,na2,na3,nam,nb1,nb2,nb3,nbm
c
      integer njintm, nj1,nj2,nj3
c
      logical ok
c
      integer nesm,nesm2,nes2,nea2,neb2,estk_s(2,nesm),estk_s0(2,nesm),
     &        estk_a(2,nesm2), estk_b(2,nesm2),estk_b0(2,nesm)
c
      integer mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     1        natyp,ims,jms,kms,is,ima,jma,imb,jmb,kma,datum,
     2        i,j,k,num,pl,jobh,ksr1,kar1,kar2,etyp,
     3        el_ch_zs,zb_xw,slice, rzero,elast
c
      real*8  t,w,c,a,kappa,alfa,r1,r2,rn,eta,my,y(0:200),z(0:200,2),
     &        rvek(15),omega,rate_b, x
c
      character ch*3,job*40,prog*20,head*45,row*80,infile*80
c
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &       /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma,
     &       /geom/ t,w,c,a,kappa,alfa,r1,r2,rn,eta,my,
     &       /reduce/ ksr1,kar1,kar2
c
c
      write(*,8001) ' '
      write(*,8001) '>> Surface crack mesh generator (focused meshes)'
      write(*,8001) '  '
      write(*,'(t1,a,$)') '>> Input file name: '
      read(*,'(a)') infile
      inquire(file=infile,exist=ok)
      if ( .not. ok ) then
            write(*,'(t4,3a)')'>>The file: ',infile,'doesn''t exist!'
            write(*,*) '>> Program aborted...'
            stop
      endif
      write(*,'(/t5,a,a/)') '> reading input data from the file: ',
     &        infile
      open(unit=25,file=infile,status='old')
 5    read(25,'(a)') ch
      if ((ch.ne.'*in').and.(ch.ne.'*IN')) goto 5
      read(25,'(a)') row
      i = 0
 10   continue
         i = i + 1
      if (row(i:i).eq.' ') goto 10
         j = i
 20   continue
         j = j + 1
      if (row(j:j).ne.' ') goto 20
      j = j - 1
      jobh = j-i+1
      job(1:jobh) = row(i:j)
      do i=jobh+1, 40
         job(i:i)=' '
      enddo
C... Program
	read(25,'(a)') prog
        if ( (prog.ne.'ADINA').and.(prog.ne.'ABAQUS') .and.
     &       (prog.ne.'WARP3D').and.(prog.ne.'PATRAN') ) then
           write(*,'(t1,a,a)') 'program name must be either',
     &      ' ADINA, ABAQUS, WARP3D or PATRAN.'
           stop
        endif
c
c... Date
	read(25,*) datum
c... MR,MF,MV,SFRED:
	read(25,'(a)') row
	call char_to_real(row,rvek,num)
	if (num.ne.4) goto 99
	mr=int(rvek(1)+0.5)
	mf=int(rvek(2)+0.5)
	mv=int(rvek(3)+0.5)
	sfred=int(rvek(4)+0.5)
c... M1,M2,NA,NATYP,OMEGA:
	read(25,'(a)') row
	call char_to_real(row,rvek,num)
	if (num.ne.5) goto 99
	m1=int(rvek(1)+0.5)
	m2=int(rvek(2)+0.5)
	na=int(rvek(3)+0.5)
	natyp=int(rvek(4)+0.5)
	omega=rvek(5)
C... MB,NB:
	read(25,'(a)') row
	call char_to_real(row,rvek,num)
	if (num.ne.4) goto 99
	mb=int(rvek(1)+0.5)
	nb=int(rvek(2)+0.5)
	zb_xw = int(rvek(3)+0.5)
 	rate_b =rvek(4)
C... LT,LRED,RTYPE:
	read(25,'(a)') row
	call char_to_real(row,rvek,num)
	if (num.ne.3) goto 99
	lt=int(rvek(1)+0.5)
	lred=int(rvek(2)+0.5)
	rtype=int(rvek(3)+0.5)
C... ETYP,EL_CH_ZS,SLICE:
        read(25,'(a)') row
	call char_to_real(row,rvek,num)
	if (num.ne.3) goto 99
	etyp=int(rvek(1)+0.5)
	el_ch_zs=int(rvek(2)+0.5)
	slice=int(rvek(3)+0.5)
C
C... R1,R2,ETA,MY:
        read(25,'(a)') row
	call char_to_real(row,rvek,num)
	if (num.ne.4) goto 99
	r1=rvek(1)
	r2=rvek(2)
 	eta=rvek(3)
	my=rvek(4)
C alternatively Read L1 and compute ETA
C        ETA = R1/L1
C
C... RN,RZERO,ELAST:
        read(25,'(a)') row
	call char_to_real(row,rvek,num)
	if (num.ne.3) goto 99
	rn=rvek(1)
	rzero = int( rvek(2)+0.5 )
	elast = int( rvek(3)+0.5 )
	if (rzero.eq.1) then
           rn=0.0001*r1/eta
	elseif (rzero.ne.0) then
	   write(*,'(t5,a)') 'rzero must be equal to 0 or 1'
	   stop
	endif
	if ((elast.ne.0).and.(elast.ne.1)) then
	   write(*,'(t5,a)') 'elast must be equal to 0 or 1'
	   stop
	endif
C
	eta=eta-rn/r1
C... T,W,C,A:
        read(25,'(a)') row
	call char_to_real(row,rvek,num)
	if (num.ne.4) goto 99
	t=rvek(1)
	w=rvek(2)
	c=rvek(3)
	a=rvek(4)
C... ZB_XW,RATE_B:
C        READ(25,'(A)') ROW
C	CALL CHAR_TO_REAL(ROW,RVEK,NUM)
C	IF (NUM.NE.2) GOTO 99
C	ZB_XW=RVEK(1)
C	RATE_B=RVEK(2)
C... ETYP,EL_CH_ZS,SLICE:
C        READ(25,'(A)') ROW
C	CALL CHAR_TO_REAL(ROW,RVEK,NUM)
C	IF (NUM.NE.3) GOTO 99
C	ETYP=INT(RVEK(1)+0.5)
C	EL_CH_ZS=INT(RVEK(2)+0.5)
C	SLICE=INT(RVEK(3)+0.5)
C
	if (.not.( (etyp.eq.8).or.(etyp.eq.20).or.(etyp.eq.27) ) ) then
	  write(*,'(/t10,a)') '=> etyp should be 8, 20 or 27 !'
	  stop
	endif
	if ((prog.eq.'ABAQUS').and.(slice.gt.0)) then
	  write(*,'(/t10,a,a)') 'if abaqus-analysis slice must be zero!'
	  stop
	endif
	if ( ((sfred.gt.0).and.(el_ch_zs.gt.0)) .and.
     &        (el_ch_zs.gt.sfred) ) then
	  write(*,'(/t10,a)') 'note that  el_ch_zs =< sfred '
	  stop
	endif
C
C... PL plotting parameter PL=1 allways (modification 941005):
C        READ(25,'(A)') ROW
C	CALL CHAR_TO_REAL(ROW,RVEK,NUM)
C	IF (NUM.NE.1) GOTO 99
C	PL=INT(RVEK(1)+0.5)
	pl=1
C... Yi, Z1i, Z2i :
	do i=0, lt
           read(25,'(a)') row
	   call char_to_real(row,rvek,num)
	   if (num.ne.3) goto 99
	   y(i)=rvek(1)
	   z(i,1)=rvek(2)
	   z(i,2)=rvek(3)
	   row='***empty**'
	enddo
	close(25)
	write(head,'(a,a,i6,a)')
     &     job(1:10),' -  ',datum,'  - semi ellip surf crack'
C.......BER[KNA KAPPA OCH ALFA :
	kappa=c/a
	alfa=dsqrt(c**2.d0-a**2.d0)/2.d0
C.......Kolla indata och Ber{kning av ett antal konstanter :
	if (sfred.eq.0) then
	   mh=mf/2-mv
	else
	   mh=mf/4-mv
	endif
	jms=2*na+1
        ims=2*mf+1
	kms=2*mr+1
        is=2*m1+1
	ma=m1+m2
	ima=2*ma+1
        jma=2*na+1
	imb=2*mb+1
	jmb=2*nb+1
C .... If RTYPE=1, one extra ("help") element row is included or if
C .... RTYPE=2, two extra ("help") element rows is included.
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
	call indata_test(mr,mf,mv,mh,m1,m2,na,mb,nb,lred,rtype,sfred,
     &       t,a,c,lt,r1,r2,rn,rzero,eta,y,z,kappa,etyp,prog)
c
c  Set matrix limits
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
        na1 = ima
        na2 = jma
        na3 = kma
        if ((na1*na2*na3).gt.nam) then
          write(*,'(t1,a)') '>> increase the size of array noda(i,j,k)'
          write(*,'(t4,a,i8)') 'current  size = ', nam
          write(*,'(t4,a,i8)') 'required size = ', na1*na2*na3
          stop
        endif
        nb1 = imb
        nb2 = jmb
        nb3 = kma
        if ((nb1*nb2*nb3).gt.nbm) then
          write(*,'(t1,a)') '>> increase the size of array nodb(i,j,k)'
          write(*,'(t4,a,i8)') 'current  size = ', nbm
          write(*,'(t4,a,i8)') 'required size = ', nb1*nb2*nb3
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
              estk_s(j,i) = 0
              estk_s0(j,i) = 0
              estk_b0(j,i) = 0
           enddo
        enddo
        do i=1, nesm2
           do j=1,2
              estk_a(j,i) = 0
              estk_b(j,i) = 0
           enddo
        enddo
c
        nj1 = (jms+1)/2
        nj2 = (kms+1)/2 + max(m1,m2)
        nj3 =  ims + 1
        if ((nj1*nj2*nj3).gt.njintm) then
          write(*,'(t1,a)') '>> increase the array size par. njintm'
          write(*,'(t4,a,i8)') 'current  size = ', njintm
          write(*,'(t4,a,i8)') 'required size = ', nj1*nj2*nj3
        endif
c
	return
99	write(*,'(t10,a/t10,a)')
     1  '* the number of data given does not match the number of',
     2  '  data required at each line, => correct the in-put file!!'
	stop
c
 8001	format(t1,a)
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine indata_test(mr,mf,mv,mh,m1,m2,na,mb,nb,lred,rtype,
     &    sfred,t,a,c,lt,r1,r2,rn,rzero,eta,y,z,kappa,etyp,prog)
C	***********************************************************
C	*             TESTAR INDATA                               *
C	***********************************************************
	integer mr,mf,mv,mh,m1,m2,na,mb,nb,lred,rtype,sfred,
     &          lt,i, rzero, etyp
	real*8  t,a,c,r1,r2,rn,eta,y(0:200),z(0:200,2),kappa,chi
        character prog*20
C
	if (mf.lt.4) then
	  write(*,'(/t10,a/t10,a)')
     1    '* the no. of elements in the circumferential direction',
     2    '  in the tubular region must be > 4'
	  goto 99
	endif
	if (mr.lt.4) then
	  write(*,'(/t10,a/t10,a)')
     1    '* the no. of elements in the radial direction',
     2    '  in the tubular region must >  4'
	  goto 99
	endif
	if (sfred.gt.0) then
	   if (mod(mf,4).ne.0) then
	      write(*,'(/t10,a)')'* if sfred > 0, then mod(mf,4)=0'
	      goto 99
	   endif
           if (mf.lt.8) then
	      write(*,'(/t10,a)') '* if mfred > 0, then mf > = 8 !'
	      goto 99
	   endif
	   if (mv.ge.(mf/4.))then
	      write(*,'(/t10,a)')
     &        'the condition: mv < mf/4, not fulfilled !'
	      goto 99
	   endif
           if (sfred.ge.mr) then
	      write(*,'(/t10,a)')
     &        'the condition:  mfred < mr, not fulfilled !'
	      goto 99
	   endif
           if (sfred.le.2) then
	      write(*,'(/t10,a)')
     &        'the condition:  mfred > 2, not fulfilled !'
	      goto 99
	   endif
	else
	   if (mv.ge.(mf/2))then
	      write(*,'(/t10,a)')
     &        'the condition:  mv < mf/2, not fulfilled !'
	      goto 99
	   endif
	endif
        if (rzero.eq.0) then
	   if ( (rn.lt.(0.02d0*r1/eta)).or.
     &          (rn.gt.(10.d0*r1/eta)) ) then
	      write(*,'(/t10,a/t10,a/t10,a/t10,a)')
     1        'in case of a notch at the crack front (rzero=0)',
     2        '  the follwing condition for rn must be fulfilled:',
     3        '      0.02*r1 < eta*rn < 10*r1',
     4        '       => choose a new rn !!'
	      goto 99
	   endif
	endif
	if ( (na.lt.6).or.(na.gt.24) ) then
	   write(*,'(/t10,a,t10,a)')
     1     'the number of element along the crack front (=na)',
     2     '  must be between 6 and 24 !'
	   goto 99
	endif
	if ((mod((na-nb),2).ne.0).and.(rtype.eq.1)) then
	   write(*,'(/t10,a,i3,a)')
     1         'if element reduction:  rtype=1, then na-nb (=',
     2          na-nb,') must an even number, not fulfilled!'
 	   goto 99
	endif
	if ( (m1-mh).lt.2) then
	 write(*,'(/t10,a)')'the condition: m1-mh >= 2, not fulfilled!'
	  goto 99
	endif
	if ( (m2-mh).lt.2) then
	 write(*,'(/t10,a)')'the condition: m2-mh >= 2, not fulfilled!'
	  goto 99
	endif
	if (mod(m1+m2,2).ne.0) then
	   write(*,'(/t10,a)') 'm1+m2 have to be an even number !'
	   goto 99
	endif
 	if ( (mod(nb,2).ne.0).and.(rtype.eq.2) ) then
	   write(*,'(/t10,a)')
     &     'nb have to be an even number, when bred .ne. 0 ! '
	   goto 99
	endif
 	if (lred.le.mv) then
	  write(*,'(/t10,a)')
     &               'the condition: lred > mv , not fulfilled !'
	  goto 99
	endif
	if (a.gt.c) then
	  write(*,'(/t10,a)')  'the condition: c > a , not fulfilled !'
	  goto 99
	endif
	if (a.gt.t) then
	  write(*,'(/t10,a)') 'the condition: t > a , not fulfilled !'
	  goto 99
	endif
C	CHI=( DBLE(MV)*DSIND( 180.D0/DBLE(MF) ) )/
C     & (((1.D0+KAPPA)/(2.D0*KAPPA))*DSIND(180.D0*DBLE(MV)/DBLE(MF) ) )
C	IF (R2.LT.(CHI*R1)) THEN
C	   WRITE(*,'(/T10,A,G10.3,A)') 'Note that R2 > ',CHI,' * R1 !'
C	   GOTO 99
C	ENDIF
	do i=1, lt
	   if ( ( z(i,2).lt.z(i,1) ).or.( y(i).lt.y(i-1) ) ) then
	      write(*,'(/t10,a,i2,a,i2,a)') 'villkoret : z(' ,i,
     &        ',2) > z(' ,i, ',1) , ej uppfyllt !'
	      goto 99
	   endif
	enddo
c
        if ((prog.eq.'WARP3D').and.(etyp.ne.8)) then
           write(*,'(a)') 'if warp3d-input 8-node element must be used'
           goto 99
        endif
c
	return
99	continue
	   write(*,'(/t10,a/t12,a)')
     &     '* the modellparameters is/are incorrect !',
     &     '=> correct the modell.'
	   stop
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine char_to_real(row,rvek,num)
C	--- The routine transform a character to a real number !
	implicit   none
	integer    num,n,k,i,exponent,es
	real*8     rvek(15),sign
	character  row*80
 
	num=0
	n=1
5	continue
	   i=0
	   k=0
	   exponent=0
	   es=1
	   sign=1.d0
10	   continue
C             ! Searcing for the 1:st digit !
	   if ( ( (row(n:n).lt.'0').or.(row(n:n).gt.'9') ).and.
     &          (row(n:n).ne.'.') )  then
	      n=n+1
	      if (n.gt.80) goto 99
              goto 10
	   endif
C              ! Negative sign !
	   if (row((n-1):(n-1)).eq.'-') sign=-1.d0
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
	   num=num+1
	   rvek(num)=sign*dble(i)*10.d0**dble(k+es*exponent)
	if (num.lt.80) goto 5
99	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine char_to_int(row,ivek,num)
c-- The routine transform a character to an integer number
      integer    ivek(1),num,n,i
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
