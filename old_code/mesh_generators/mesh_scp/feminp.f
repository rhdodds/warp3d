c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      SUBROUTINE ADINA_IN(JOB,JOBH,EL_CH_ZS,SLICE,EL_JCO,NJ1,NJ2,NJ3,
     1           NO_OF_JCONT,VSH_VEC,CR_FRONT,DSL,GRP,ETYP,ELAST,
     2           NODS,NS1,NS2,NS3, NODA,NA1,NA2,NA3, NODB,NB1,NB2,NB3)
c
c--- The routine creates an indata file (job.in) to an ADINA analysis.
c
      IMPLICIT NONE
c
      INCLUDE 'mesh3d_scp_common_nod.f'
c
      INTEGER  NS1,NS2,NS3,NODS(NS1,NS2,NS3),
     1         NA1,NA2,NA3,NODA(NA1,NA2,NA3),
     2         NB1,NB2,NB3,NODB(NB1,NB2,NB3)
c
      INTEGER  NJ1,NJ2,NJ3, EL_JCO(2,NJ1,NJ2,0:NJ3)
c
      INTEGER  EL_CH_ZS,SLICE,JOBH,UNI,DSL(200),GRP,ETYP, NN,I,J,
     1         NO_OF_JCONT,CR_FRONT(500,3), ELAST, SLAVE,MAST,
     2         IMS,JMS,KMS,IS,IMA,JMA,IMB,JMB,KMA,
     3         NSET(20000)
c
      REAL*8   VSH_VEC(200,3)
c
      CHARACTER JOB*40,FILIN*40
      COMMON   /MAX/ IMS,JMS,KMS,IS,IMA,JMA,IMB,JMB,KMA

C.....File names
	FILIN=JOB(1:JOBH)//'.in'
	OPEN(UNIT=14,FILE=FILIN,STATUS='UNKNOWN')

	UNI=14
	WRITE(UNI,'(T1,A)') '** MODEL - information **'
	WRITE(UNI,'(T1,A,A)') '** job : ',JOB(1:JOBH)
C
	WRITE(UNI,'(T1,A)') 'DATABASE SCRATCH'
	WRITE(UNI,'(T1,A)') 'CONTROL ECHO=YES'
	WRITE(UNI,'(T1,A)') 'FILE ECHO=7 LOG=8'
	WRITE(UNI,'(T1,A,A)') 'MASTER IDOF=000111 REACT=YES',
     &  ' MODEX=EXECUTE N=400 DT=0.25E-02 IR=400'
	WRITE(UNI,'(T1,A)') 'ANALYSIS TYPE=STATIC'
	WRITE(UNI,'(T1,A)') 'KINEMATICS DISPL=L STRAIN=L'
	WRITE(UNI,'(T1,A)') '*'
	WRITE(UNI,'(T1,A)') '*    Solution Method'
	WRITE(UNI,'(T1,A)') 'SOLUTION METHOD=1 SHIFT=1.02'
	WRITE(UNI,'(T1,A)') '*'
	WRITE(UNI,'(T1,A)') 'TOLERANCE ETOL=1.E-05'
	WRITE(UNI,'(T1,A)') '*'
	WRITE(UNI,'(T1,A)') '*    Save Specification'
	WRITE(UNI,'(T1,A)') 'ESAVESTEPS FIRST1=50 LAST1=200 INCR1=50'
	WRITE(UNI,'(T1,A)') 'NSAVESTEPS FIRST1=50 LAST1=200 INCR1=50'
	WRITE(UNI,'(T1,A)') '*'
	WRITE(UNI,'(T1,A)') '*    Node Definition    *'
	WRITE(UNI,'(T1,A)') 'COORDINATES / ENTRIES NODE X Y Z'
	WRITE(UNI,'(T1,A)') 'READ 15'
	WRITE(UNI,'(T1,A)') '* SAVENODES'
	WRITE(UNI,'(T1,A)') '* READ 16'
	WRITE(UNI,'(T1,A)') '*'
	WRITE(UNI,'(T1,A)') '*    Element Definition    *'
C . . Elset no. 1
        WRITE(UNI,'(T1,A)') 
     &   'EGROUP N=1 THREEDS  MAT=1 RSINT=3 TINT=3 FORMULATION=1'
	IF (ETYP.EQ.27) THEN
	   WRITE(UNI,'(T1,A)') 
     &    'ENODES / ENTRIES EL N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12,',
     &    'N13 N14 N15 N16 N17 N18 N19 N20 N21 N22 N23 N24 N25 N26 N27'
	ELSEIF (ETYP.EQ.20) THEN
	   WRITE(UNI,'(T1,A)') 
     &    'ENODES / ENTRIES EL N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12,',
     &    'N13 N14 N15 N16 N17 N18 N19 N20'
	ELSE
	   WRITE(UNI,'(T1,A)') 
     &    'ENODES / ENTRIES EL N1 N2 N3 N4 N5 N6 N7 N8'
	ENDIF
	WRITE(UNI,'(T1,A)') 'READ 51'
	WRITE(UNI,'(T1,A)') 'EDATA / ENTRIES EL SAVE'
	WRITE(UNI,'(T1,A,I5,A)')' 1 yes / STEP 1 TO / ',DSL(1),' yes'
	WRITE(UNI,'(T1,A)') '*'
	IF ((SLICE.GT.0).OR.(EL_CH_ZS.GT.0)) THEN
C . . Elset no. 2 -
	  DO I=2, GRP
           WRITE(UNI,'(T1,A,I2,A)') 
     &     'EGROUP N=',I,' THREEDS  MAT=1 RSINT=2 TINT=2'
	  IF (ETYP.EQ.27) THEN
	   WRITE(UNI,'(T1,A)') 
     &    'ENODES / ENTRIES EL N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12,',
     &    'N13 N14 N15 N16 N17 N18 N19 N20 N21 N22 N23 N24 N25 N26 N27'
	  ELSEIF (ETYP.EQ.20) THEN
	   WRITE(UNI,'(T1,A)') 
     &    'ENODES / ENTRIES EL N1 N2 N3 N4 N5 N6 N7 N8 N9 N10 N11 N12,',
     &    'N13 N14 N15 N16 N17 N18 N19 N20'
	  ELSE
	   WRITE(UNI,'(T1,A)') 
     &    'ENODES / ENTRIES EL N1 N2 N3 N4 N5 N6 N7 N8'
	  ENDIF
	   WRITE(UNI,'(T1,A,I2)') 'READ ',50+I
	   WRITE(UNI,'(T1,A)') 'EDATA / ENTRIES EL SAVE'
	   WRITE(UNI,'(T1,A,I5,A)')' 1 yes / STEP 1 TO / ',DSL(I),' yes'
	   WRITE(UNI,'(T1,A)') '*'
	  ENDDO
	ENDIF
	WRITE(UNI,'(T1,A)') '*'
	WRITE(UNI,'(T1,A)') '*    Material Law   *'
	WRITE(UNI,'(T1,A)') '* MATERIAL N=1 ELASTIC E=206.E+09 NU=0.300'
	WRITE(UNI,'(T1,A)')
     &  'MATERIAL N=1 PLASTIC-M HARDENING=ISOTR E=206.E+09 NU=0.300'
	WRITE(UNI,'(T1,A)') 'READ 41'
	WRITE(UNI,'(T1,A)') '*'
	WRITE(UNI,'(T1,A)') '*    Constraint Equations   *'
	WRITE(UNI,'(T1,A)') 'CONSTRAINTS / READ 26'
	CALL CONSTRAINT_EQU_ADINA(NODA,NA1,NA2,NA3,
     &                            IMA,KMA,JOBH,JOB,26)
	IF (ELAST.EQ.1) THEN
	   WRITE(UNI,'(T1,A)') '*'
	   WRITE(UNI,'(T1,A/T1,A)')
     1     '**  Constraint Equations for degenerated elements along *',
     2     '**  the crack front in the case of an elastic analysis. *'
	   WRITE(UNI,'(T1,A)') 'CONSTRAINTS'
	   DO J=1, JMS
	      MAST=NNR(NODS(1,J,1))
	      DO I=2, (IMS+1)/2
                 IF ( NNR(NODS(I,J,1)).GT.0 ) THEN
	           SLAVE=NNR(NODS(I,J,1))
	           WRITE(UNI,'(T1,2(I5,A))') SLAVE,',1,',MAST,',1, 1.0'
	           WRITE(UNI,'(T1,2(I5,A))') SLAVE,',2,',MAST,',2, 1.0'
	           WRITE(UNI,'(T1,2(I5,A))') SLAVE,',3,',MAST,',3, 1.0'
	         ENDIF
	      ENDDO
	   ENDDO
	ENDIF
	WRITE(UNI,'(T1,A)') '*'
	WRITE(UNI,'(T1,A)') '*    Boundary Condition   *'
	WRITE(UNI,'(T1,A)') '* crackfront plane NFRONT'
	WRITE(UNI,'(T1,A)') 'FIXBOUNDARIES DIREC=2 TYPE=NODES'
C
        CALL NODSET_NFRONT_(NSET,NN, NODS,NS1,NS2,NS3, 
     &         NODA,NA1,NA2,NA3, NODB,NB1,NB2,NB3 )
        CALL WRITE_OUT_SET(UNI,NSET,NN)
        WRITE(UNI,'(T1,A,I5,A)') '**  NFRONT contains ',NN,' nodes'
C
	WRITE(UNI,'(T1,A)') '* symmetry plane NSYM'
	WRITE(UNI,'(T1,A)') 'FIXBOUNDARIES DIREC=1 TYPE=NODES'
        CALL NODSET_NSYM_(NSET,NN,ELAST, NODS,NS1,NS2,NS3, 
     &                                   NODA,NA1,NA2,NA3)
        CALL WRITE_OUT_SET(UNI,NSET,NN)
        WRITE(UNI,'(T1,A,I5,A)')'** NSYM contains ',NN,' nodes'
C
	WRITE(UNI,'(T1,A)') '* back plane NBAK'
	WRITE(UNI,'(T1,A)') 'FIXBOUNDARIES DIREC=13 TYPE=NODES'
        CALL NODSET_NBAK_(NSET,NN, NODA,NA1,NA2,NA3, NODB,NB1,NB2,NB3)
        CALL WRITE_OUT_SET(UNI,NSET,NN)
        WRITE(UNI,'(T1,A,I5,A)')'** NBAK contains ',NN,' nodes'
C
C	WRITE(*,'(T5,A,A,$)') 'Shall the crack surface be ',
C     &             'constrained in the 2-direction ? (Y/N): '
C	READ(*,'(A)') ANSWER
C	IF ((ANSWER.EQ.'Y').OR.(ANSWER.EQ.'y')) THEN
C	   WRITE(UNI,'(T1,A)') '* crack face NCRSUR'
C	   WRITE(UNI,'(T1,A)') 'FIXBOUNDARIES DIREC=2 TYPE=NODES'
C	   CALL NODSET_NCRSUR_(NSET,NN,ELAST, NODS,NS1,NS2,NS3,
c      &                                      NODA,NA1,NA2,NA3)
c        WRITE(UNI,'(T1,A,I5,A)')'** NCRSUR contains ',NN,' nodes'
C	ENDIF
	WRITE(UNI,'(T1,A)') '*'
	WRITE(UNI,'(T1,A)') '*    Loading Conditions  *'
	WRITE(UNI,'(T1,A)') 'LOADS DISPLACEMENTS TYPE=NODES'
	WRITE(UNI,'(T1,I5,A,I2,A)') NNR(NODA(1,1,KMA)),  ',2,??,1, 0.0'
	WRITE(UNI,'(T1,I5,A,I2,A)') NNR(NODA(IMA,1,KMA)),',2,??,1, 0.0'
	WRITE(UNI,'(T1,A)') '*'
	WRITE(UNI,'(T1,A)') 'TIMEFUNCTION N=1'
	WRITE(UNI,'(T1,A)') '0.0  0.0'
	WRITE(UNI,'(T1,A)') '1.0  1.0'
C... Include J-integral contours and the virtual shift vector
	WRITE(UNI,'(T1,A)') '* For J-integral computations'
	WRITE(UNI,'(T1,A)') 'FRACTURE METH=VIRTUAL IFDIM=3'
	WRITE(UNI,'(T1,A,I2)') 'CRACK-PROP NCRACK=',2
	DO J=1, JMS, 2
	   WRITE(UNI,'(T1,2I6)') ( CR_FRONT(J,I),I=1, 2 )
	ENDDO
	NN=0
	DO J=1, (JMS-1)/2
	   WRITE(UNI,'(T1,A)') '**'
	   WRITE(UNI,'(T1,A,I2,A)') '* ELEMENT SLICE NO.=',J,' ****'
	   WRITE(UNI,'(T1,A)') '**'
	   DO I=1, NO_OF_JCONT
	      NN=NN+1
	      WRITE(UNI,'(T1,A,I3,A,2(A,D17.10))') 
     &             'VSHIFT ',NN,' TYPE=ELEM VECTOR=INPUT',
     &             ' VX=',VSH_VEC(J,1),' VZ=',VSH_VEC(J,3)
	      CALL  WRITE_OUT_JCONT(UNI,EL_JCO,NJ1,NJ2,NJ3,J,I)
	   ENDDO
	ENDDO
	WRITE(UNI,'(T1,A)') 'CONTROL ECHO=YES'
	WRITE(UNI,'(T1,A)') 'FILE LOG=6 ECHO=6'
	WRITE(UNI,'(T1,A)') 'FILE 5'
	RETURN
	END
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine write_out_jcont(uni,el_jco,nj1,nj2,nj3,j,i)
	implicit none
	integer  uni, nj1,nj2,nj3,el_jco(2,nj1,nj2,0:nj3),
     &           j,i, ii,k,kk,n,ejc(1000)
	character form*29
	n=0
	do ii=1, i
	   do k=1, el_jco(1,j,ii,0)
	      n=n+2
	      ejc(n-1)=el_jco(1,j,ii,k)
	      ejc(n)  =el_jco(2,j,ii,k)
	   enddo
	enddo
	do k=1, n, 20
	   if (((k+19)/20).le.int(n/20)) then
	      write(uni,101) ( ejc(kk), kk=k, k+19 )
	   else
	      if (mod(n,20).eq.2) then
	         write(uni,102) ( ejc(kk), kk=k, n )
	      elseif (mod(n,20).eq.4) then
	         write(uni,103) ( ejc(kk), kk=k, n )
	      else
	         write(form,'(a,i1,a)')   '(i4,i3,',mod(n,20)/2-2,
     &                      '(''/'',i4,i3),''/''i4,i3)'
	         write(uni,form) ( ejc(kk), kk=k, n )
	      endif
	   endif
	enddo
101     format(i4,i3,9('/',i4,i3))
102     format(i4,i3)
103     format(i4,i3,'/',i4,i3)
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine constraint_equ_adina(noda,na1,na2,na3,
     &                       ima,kma,jobh,job,uni)
C--- The routine creates constraint equations
      implicit  none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  na1,na2,na3,noda(na1,na2,na3), ima,kma,jobh,m1,m2,uni,i
      real*8   l0,z1,z,b1,b2
      character job*40,fil26*40
c
	fil26=job(1:jobh)//'.26'
	open(uni,file=fil26,status='unknown')
	m1=nnr(noda(1,1,kma))
	m2=nnr(noda(ima,1,kma))
	z1=npos(noda(1,1,kma),3)
	l0=npos(noda(ima,1,kma),3)-npos(noda(1,1,kma),3)
	do i=1, inm
	   if ((nnr(i).gt.m1).and.(nnr(i).ne.m2)) then
	      z=npos(i,3)-z1
	      b1=(l0-z)/l0
	      b2=z/l0
	      if (b1.lt.(l0*1.d-7)) then
	         write(uni,122) nnr(i), 2, m2, 2, b2
	      elseif (b2.lt.(l0*1.d-7)) then
	         write(uni,122) nnr(i), 2, m1, 2, b1
	      else
	         write(uni,121) nnr(i), 2, m1, 2, b1, m2, 2, b2
	      endif
	   endif
	enddo
121     format(t1,i5,',',i2,2(',',i5,',',i2,',',d19.12))
122     format(t1,i5,',',i2,',',i5,',',i2,',',d19.12)
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine create_j_contours(estk_a,nea2,el_ch_zs,nods,
     1           ns1,ns2,ns3,eln_d,no_of_jcont,vsh_vec,cr_front,
     2           el_jco,nj1,nj2,nj3)
c
c-- The subroutine creates contours for J-integral evaluation with
c-- the so called virtual crack-extension-method.
c
c   EL_JCO(2,NJ1,NJ2,0:NJ3)
c   NO_OF_JCONT
c   VSH_VEC(200,3)
c
      implicit none
c
      include 'mesh3d_scp_common_eln.f'
      integer eln_d(0:iem)
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3)
c
      integer  nj1,nj2,nj3, el_jco(2,nj1,nj2,0:nj3)
c
      integer  nea2,estk_a(2,nea2,1),el_ch_zs, no_of_jcont,
     1         cr_front(500,3),ejc(100,100,2), ea(30,50,50),
     2         i,j,k, ii,kk,e1,i2,ee, t,t1,t2
c
C COMMON INTEGER
      integer  mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &           ims,jms,kms,is,ima,jma,imb,jmb,kma
c
      real*8   vsh_vec(200,3),norm(3),dx,dz,l0
C TEST VARIABLES
C     REAL*8   A0,B0,X0,Z0,TN(3),TNL
 
	common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
        common /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
 
C... Determine the  number of contours
	if (min((m1-mh-2),(m2-mh-2)).gt.0) then
	   no_of_jcont=mr+min((m1-mh-2),(m2-mh-2))
	else
	   no_of_jcont=mr
	endif
	if (sfred.gt.0) no_of_jcont=no_of_jcont-1
C... DETERMINE THE CONTOURS
C . . Zon S
C	E1=(MR*MF-1)*(J-1)
	e1=0
	if (sfred.eq.0) then
	   i2=mr
	else
	   i2=mr-1
	endif
	do j=1, (jms-1)/2
	   e1=e1+1
	   ejc(j,1,1)=e1
	   e1=e1+mf-2
	   ejc(j,1,2)=e1
	   do i=2, i2
	      e1=e1+1
	      ejc(j,i,1)=e1
	      if ((sfred.gt.0).and.(i.eq.sfred)) then
	         e1=e1+mf+mf/2-1
	      elseif ((sfred.gt.0).and.(i.gt.sfred)) then
	         e1=e1+mf/2-1
	      else
	         e1=e1+mf-1
	      endif
	      ejc(j,i,2)=e1
	   enddo
	enddo
	do j=1, (jms-1)/2
	   do i=1, i2
	      ee=0
	      do k=ejc(j,i,1), ejc(j,i,2)
	         ee=ee+1
                 el_jco(1,j,i,ee)=eln(k,0)
                 el_jco(2,j,i,ee)=eln_d(k)
	      enddo
	      el_jco(1,j,i,0)=ee
	   enddo
	enddo
C . . Zon A
	if (min((m1-mh-2),(m2-mh-2)).gt.0) then
C...      Create an element matrix in Zone A in order to make the gen-
	   do k=1, mv
	      do j=1, (jma-1)/2
	         ii=-1
	         do i=(ima-1)/2, 3, -1
	            if ( (i.gt.(m1+mh)).or.(i.lt.(m1+1-mh)) ) then
                       ii=ii+1
                       ea(k,j,i)=estk_a(2,k,j)-ii
	            endif
	         enddo
	      enddo
	   enddo
	   do k=mv+1, mv+min( (m1-mh-2),(m2-mh-2) )
	      do j=1, (jma-1)/2
	         do i=(ima-1)/2, 3, -1
	            ii=(ima-1)/2-i
	            ea(k,j,i)=estk_a(2,k,j)-ii
	         enddo
	      enddo
	   enddo
	   do j=1, (jma-1)/2
	      do i=i2+1, no_of_jcont
	         ee=0
	         kk=i-i2-1
	         t1=m1+mh+kk+1
	         t2=m1-mh-kk
	         do k=1, mv+kk
	            ee=ee+1
	            ii=ea(k,j,t1)
                    el_jco(1,j,i,ee)=eln(ii,0)
                    el_jco(2,j,i,ee)=eln_d(ii)
	         enddo
	         k=mv+kk+1
	         do t=t1, t2, -1
	            ee=ee+1
	            ii=ea(k,j,t)
                    el_jco(1,j,i,ee)=eln(ii,0)
                    el_jco(2,j,i,ee)=eln_d(ii)
	         enddo
	         do k=mv+kk, 1, -1
	            ee=ee+1
	            ii=ea(k,j,t2)
                    el_jco(1,j,i,ee)=eln(ii,0)
                    el_jco(2,j,i,ee)=eln_d(ii)
	         enddo
	         el_jco(1,j,i,0)=ee
	      enddo
	   enddo
	endif
C...... Determine the crack fronts nodes
	do j=1, jms, 2
	   cr_front(j,1)=nnr(nods(1,j,1))
	   cr_front(j,2)=nnr(nods(1,j,3))
	   cr_front(j,3)=nnr(nods(1,j,5))
	enddo
C...... Determine the virtual shift ( VSH_VEC(200,3) )
C	A0=NPOS(NODS(1,JMS,1),1)
C	B0=NPOS(NODS(1,1,1),3)
C	WRITE(49,*) ' A0=',A0,' B0=',B0
	do j=2, jms, 2
	   dx=dabs(npos(nods(1,j-1,1),1)-npos(nods(1,j+1,1),1))
	   dz=dabs(npos(nods(1,j-1,1),3)-npos(nods(1,j+1,1),3))
	   norm(1)= dz / dsqrt(dx*dx+dz*dz)
	   norm(2)= 0.d0
	   norm(3)= dx / dsqrt(dx*dx+dz*dz)
	   l0=dsqrt((npos(nods(1,j+1,1),1)-npos(nods(1,j+1,3),1))**2.d0
     &         + (npos(nods(1,j+1,1),3)-npos(nods(1,j+1,3),3))**2.d0 )
	   do i=1, 3
	      vsh_vec(j/2,i)=l0*norm(i)*0.01
	   enddo
C	   WRITE(49,*) 'J=',J,' norm: ',(REAL(NORM(I)), I=1, 3)
C TEST
C	   X0=NPOS(NODS(1,J,1),1)
C	   Z0=NPOS(NODS(1,J,1),3)
C	   TN(1)= (X0/(A0*A0))
C	   TN(3)= (Z0/(B0*B0))
C	   TNL= DSQRT( TN(1)*TN(1)+TN(3)*TN(3) )
C	   TN(1)= TN(1) / TNL
C	   TN(2)= 0.D0
C	   TN(3)= TN(3) / TNL
C	   WRITE(49,*) '       TN : ',(REAL(TN(I)), I=1, 3)
C	   WRITE(49,*) '  '
	enddo
 
C......TEST
C	DO I=(IMA-1)/2, 3, -1
C	   WRITE(49,'(10I5)')
C     &     ( EA(K,1,I), K=MV+MIN((M1-MH-2),  (M2-MH-2)), 1, -1 )
C	ENDDO
C	WRITE(49,*) ' J-CONTORS'
C	DO J=1, (JMS-1)/2
C	   WRITE(49,*) ' ---- J=',J,' -------'
C	   DO I=1, NO_OF_JCONT
C	      WRITE(49,'(T1,A,I2,20I6)') 'I=',I,
C     &          ( EL_JCO(1,J,I,K), K=1, EL_JCO(1,J,I,0) )
C	      WRITE(49,'(T1,A,I2,20I3)') 'I=',I,
C     &            ( EL_JCO(2,J,I,K), K=1, EL_JCO(1,J,I,0) )
C	   ENDDO
C	ENDDO
C	DO J=1, (JMS-1)/2
C	   WRITE(49,*) 'J=',J,' VSH-VEC: ',(REAL(VSH_VEC(J,I)), I=1, 3)
C	ENDDO
 
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      SUBROUTINE ABAQUS_INP(JOB,JOBH,HEAD,ETYP,
     1           NEA2,NEB2, ESTK_S,ESTK_S0, ESTK_A, ESTK_B,ESTK_B0,
     2           EL_CH_ZS, NO_OF_JCONT,ELAST,
     3           NODS,NS1,NS2,NS3, NODA,NA1,NA2,NA3, NODB,NB1,NB2,NB3)
c
C--- The routine creates an indata file (job.inp) to an ABAQUS analysis.
c
      IMPLICIT NONE
c
      INCLUDE 'mesh3d_scp_common_nod.f'
c
c
      INTEGER  NS1,NS2,NS3,NODS(NS1,NS2,NS3),
     1         NA1,NA2,NA3,NODA(NA1,NA2,NA3),
     2         NB1,NB2,NB3,NODB(NB1,NB2,NB3)
c
      INTEGER  NEA2,NEB2, ESTK_S(2,1),ESTK_S0(2,1), ESTK_A(2,NEA2,1),
     &         ESTK_B(2,NEB2,1),ESTK_B0(2,1)
c
      INTEGER  JOBH,ETYP, ELAST, NSET(20000),JCMAX,
     1         EL_CH_ZS,UNI,I,J,K, NO_OF_JCONT, KEND
C
      CHARACTER HEAD*45,JOB*40,FILINP*40,FIL1*40,FIL2*40,FIL3*40,
     &          FIL4*40,FIL5*40,FIL6*40
      CHARACTER*4 JCO(49)
      INTEGER          IMS,JMS,KMS,IS,IMA,JMA,IMB,JMB,KMA
      COMMON     /MAX/ IMS,JMS,KMS,IS,IMA,JMA,IMB,JMB,KMA
      integer      mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
c
      FIL1=JOB(1:JOBH)//'.015'
      FIL2=JOB(1:JOBH)//'.016'
      FIL3=JOB(1:JOBH)//'.017'
      FIL4=JOB(1:JOBH)//'.018'
      FIL5=JOB(1:JOBH)//'.019'
      FIL6=JOB(1:JOBH)//'.020'
      FILINP=JOB(1:JOBH)//'.inp'
      OPEN(UNIT=14,FILE=FILINP,STATUS='UNKNOWN')
	UNI=14
	WRITE(UNI,'(T1,A)') '*HEADING'
	WRITE(UNI,'(T1,A)') HEAD
	WRITE(UNI,'(T1,A)') '**'
	WRITE(UNI,'(T1,A)') '*PREPRINT,MODEL=NO,ECHO=NO'
	WRITE(UNI,'(T1,A,A)') '*NODE,NSET=NALL,INPUT=',FIL1
        WRITE(UNI,'(T1,A,I2,A,A)') '*ELEMENT,ELSET=ELEMALL,TYPE=',ETYP,
     &    ',INPUT=',FIL2
	IF (EL_CH_ZS.GT.0) THEN
           WRITE(UNI,'(T1,A,I2,A,A)') 
     &      '*ELEMENT,ELSET=ELEMALL,TYPE=',ETYP,',INPUT=',FIL3
	ENDIF
	CALL MAKE_NODSET(UNI,NSET,ELAST, NODS,NS1,NS2,NS3, 
     &                 NODA,NA1,NA2,NA3, NODB,NB1,NB2,NB3 )
	WRITE(UNI,'(T1,A,I1,A)')'*ELSET,ELSET=LAYERCR,GENERATE'
	IF (EL_CH_ZS.EQ.0) THEN
	   I=ESTK_S(1,1)
	ELSE
	   I=ESTK_S0(1,1)
	ENDIF
	WRITE(UNI,'(T1,I4,A,I4)') ESTK_S(1,1), ',' ,ESTK_S(2,(JMS-1)/2)
c
        IF ((LRED.LT.LT).OR.ETYP.EQ.8) THEN
           KEND = 2*LRED-3
        ELSE
           KEND = KMA
        ENDIF
	DO K=1, KEND-2, 2
	   IF (ESTK_A(1,(K+1)/2,1).GT.0) THEN
	      IF ( ((K+1)/2).LT.10) THEN
	         WRITE(UNI,'(T1,A,I1,A)')
     &              '*ELSET,ELSET=LAYER0',(K+1)/2,',GENERATE'
	      ELSE
	         WRITE(UNI,'(T1,A,I2,A)')
     &              '*ELSET,ELSET=LAYER',(K+1)/2,',GENERATE'
	      ENDIF
	      IF (ESTK_B((K+1)/2,(IMB-1)/2,1).GT.0) THEN
	         WRITE(UNI,'(T1,I4,A,I4)') ESTK_A(1,(K+1)/2,1),
     &               ',' ,ESTK_B(2,(K+1)/2,(IMB-1)/2)
	      ELSE
	         WRITE(UNI,'(T1,I4,A,I4)') ESTK_A(1, (K+1)/2, 1),
     &              ',' ,ESTK_A(2,(K+1)/2,(JMA-1)/2)
	      ENDIF
	   ENDIF
	ENDDO
c
	WRITE(UNI,'(T1,A)')
     &   '*SOLID SECTION,ELSET=ELEMALL,MATERIAL=MAT1'
	WRITE(UNI,'(T1,A)') '*MATERIAL,NAME=MAT1'
	WRITE(UNI,'(T1,A)') '*ELASTIC'
	WRITE(UNI,'(T1,A)') ' 200.0E+09, 0.3'
	WRITE(UNI,'(T1,A)') '**PLASTIC'
	WRITE(UNI,'(T1,A)') '*RESTART,WRITE,FREQ=999'
	WRITE(UNI,'(T1,A)') '**'
	WRITE(UNI,'(T1,A)') '** RESULT - INFORMATION **'
	WRITE(UNI,'(T1,A)') '**'
	WRITE(UNI,'(T1,A)') '*NSET,NSET=NODPRD'
	WRITE(UNI,'(T1,A,I5)') '1,2,',NNR(NODA(1,1,1))
	WRITE(UNI,'(T1,A)') '*NSET,NSET=NODPRF,GENERATE'
	WRITE(UNI,'(T1,I5,A,I5)')  NNR(NODS(1,1,1)),',',
     &                             NNR(NODA(IMA,JMA,13))
	WRITE(UNI,'(T1,A)') '*NSET,NSET=NODPRF'
	WRITE(UNI,'(T1,A)') 'NODPRF,NFRONT,NBAK,NCRSUR,NUND,NTOP'
	WRITE(UNI,'(T1,A)') '*ELSET,ELSET=ELPRD'
	WRITE(UNI,'(T1,A)') '1,2'
	WRITE(UNI,'(T1,A)') '*ELSET,ELSET=ELPRF'
	WRITE(UNI,'(T1,A,A)') 'LAYERCR,LAYER01,LAYER02,LAYER03,',
     &                        'LAYER04,LAYER05,LAYER06'
	WRITE(UNI,'(T1,A)') '*EQUATION,INPUT=FIL6'
	CALL CONSTRAINT_EQU_ABAQUS(IMA,KMA,FIL6,26,
     &                             NODA,NA1,NA2,NA3)
	IF (ELAST.EQ.1) THEN
 	   WRITE(UNI,'(T1,A/T1,A/T1,A)')'***',
     1      '** "Rigid links" between nodes in degenerated elements',
     2      '**         at the crack front (elastic analysis)'
           WRITE(uni,'(t1,a)') '*MPC'
	   DO J=1, JMS
              DO I=2, (IMS+1)/2
                 IF ( NNR(NODS(I,J,1)).GT.0 ) THEN
	            WRITE(UNI,'(T1,2(A,I5))') 'TIE, ',
     &                   NNR(NODS(I,J,1)),', ',NNR(NODS(1,J,1))
	         ENDIF
              ENDDO
	   ENDDO
	ENDIF
	WRITE(UNI,'(T1,A)') '*BOUNDARY'
	WRITE(UNI,'(T1,A)') ' NSYM,1'
	WRITE(UNI,'(T1,A)') ' NFRONT,2'
	WRITE(UNI,'(T1,A)') '**'
	WRITE(UNI,'(T1,A)') '** ****** S T E P = 1  ************'
	WRITE(UNI,'(T1,A)') '**'
	WRITE(UNI,'(T1,A)') '*STEP'
	WRITE(UNI,'(T1,A)') '*STATIC'
	WRITE(UNI,'(T1,A)') '*BOUNDARY'
	WRITE(UNI,'(T1,I5,A,I2,A)') NNR(NODA(1,1,KMA)),  ',2,2, 1.E-3??'
	WRITE(UNI,'(T1,I5,A,I2,A)') NNR(NODA(IMA,1,KMA)),',2,2, 1.E-3??'
	WRITE(UNI,'(T1,A)') '*NODE PRINT,NSET=NODPRD'
	WRITE(UNI,'(T1,A)') ' U,RF'
	WRITE(UNI,'(T1,A)') '*EL PRINT,ELSET=ELPRD'
	WRITE(UNI,'(T1,A)') ' S,MISES'
	WRITE(UNI,'(T1,A)') ' E'
	WRITE(UNI,'(T1,A)') '***NODE FILE,NSET=NODPRF,FREQ=10'
	WRITE(UNI,'(T1,A)') '** U,RF'
	WRITE(UNI,'(T1,A)') '***EL FILE,ELSET=ELPRF,FREQ=10'
	WRITE(UNI,'(T1,A)') '** S,SINV,ENER'
	WRITE(UNI,'(T1,A)') '** E,PE'
	WRITE(UNI,'(T1,A)') '*ENERGY FILE,FREQ=10'
	NO_OF_JCONT=MR+M1-MH-2
	WRITE(UNI,'(T1,A,I2,A)') '*CONTOUR INTEGRAL,CONTOURS=',
     &    MR+M1-MH-2,',TYPE=J,NORMAL,SYMM,OUTPUT=BOTH,FREQ=10'
	WRITE(UNI,'(T1,A)') ' 0.,-1.,0.'
	DO I=1, 9
	   JCO(I)='CR0'//CHAR(48+I)
	ENDDO
	JCO(I)='CR10'
	DO I=11, 19
	   JCO(I)='CR1'//CHAR(38+I)
	ENDDO
	JCO(I)='CR20'
	DO I=21, 29
	   JCO(I)='CR2'//CHAR(28+I)
	ENDDO
	JCO(I)='CR30'
	DO I=31, 39
	   JCO(I)='CR3'//CHAR(18+I)
	ENDDO
	JCO(I)='CR40'
	DO I=41, 49
	   JCO(I)='CR4'//CHAR(8+I)
	ENDDO
        IF (ETYP.LE.8) THEN
           JCMAX = (JMS+1)/2
        ELSE
           JCMAX = JMS
        ENDIF
	IF (JCMAX.LE.16) THEN
	   WRITE(UNI,101) ( JCO(J), J=1, JCMAX )
	ELSEIF (JCMAX.LE.32) THEN
	   WRITE(UNI,101) ( JCO(J), J=1, 16 )
	   WRITE(UNI,101) ( JCO(J), J=17, JCMAX )
	ELSE
	   WRITE(UNI,101) ( JCO(J), J=1, 16 )
	   WRITE(UNI,101) ( JCO(J), J=17, 32 )
	   WRITE(UNI,101) ( JCO(J), J=33, JCMAX )
	ENDIF
	WRITE(UNI,'(T1,A)') '*ENDSTEP'
	CLOSE(UNI)
c
        call write_out_tipnodes(nods,ns1,ns2,ns3,job,jobh)
c
101	FORMAT(T1,' ',A,15(' ',A))
	RETURN
	END
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine constraint_equ_abaqus(ima,kma,fil6,uni,
     &                                 noda,na1,na2,na3)
C--- The routine creates constraint equations for ABAQUS
c
      implicit  none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  na1,na2,na3,noda(na1,na2,na3)
c
      integer   ima,kma,uni,m1,m2,i
      real*8    l0,z1,z,b1,b2
      character fil6*14
c
	open(uni,file=fil6,status='unknown')
	m1=nnr(noda(1,1,kma))
	m2=nnr(noda(ima,1,kma))
	z1=npos(noda(1,1,kma),3)
	l0=npos(noda(ima,1,kma),3)-npos(noda(1,1,kma),3)
	do i=1, inm
	   if ((nnr(i).gt.m1).and.(nnr(i).ne.m2)) then
	      z=npos(i,3)-z1
	      b1=(l0-z)/l0
	      b2=z/l0
	      write(uni,'(t1,i1)') 3
	      write(uni,121) nnr(i), 2, -1.00, m1, 2, b1, m2, 2, b2
	   endif
	enddo
121     format(t1,i5,',',i2,',',f6.2,2(',',i5,',',i2,',',d19.12))
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine write_out_tipnodes(nods,ns1,ns2,ns3,job,jobh)
c
      implicit none
      include 'mesh3d_scp_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3)
      integer  jobh,ntip(401,2),i,j,n,io
      double precision  x,y,z,phi
      character job*40,tipfile*40
c
      integer       ims,jms,kms,is,ima,jma,imb,jmb,kma
      common  /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
      double precision  t,w,c,a,kappa,alfa,r1,r2,rn,eta,my
      common     /geom/ t,w,c,a,kappa,alfa,r1,r2,rn,eta,my
c
      n = 0
c-Zone S:
      do j=1, jms
         if (nnr(nods(1,j,1)).gt.0) then
            n=n+1
            ntip(n,1) = nnr(nods(1,j,1))
            ntip(n,2) = nods(1,j,1)
         endif
      enddo
c
      tipfile = job(1:jobh)//'_tip.inf'
      io = 19
      open(unit=io,file=tipfile,status='unknown')
      write(io,*) n
      do i=1, n
         if (i.eq.1) then
            x = npos(ntip(i,2),1)
            y = npos(ntip(i,2),2)
            z = npos(ntip(i,2),3)
            phi = acosd(x/c)
            write(io,111) ntip(i,1), x,y,z, 90.0
        elseif (i.eq.n) then
            x = npos(ntip(i,2),1)
            y = npos(ntip(i,2),2)
            z = npos(ntip(i,2),3)
            phi = acosd(x/c)
            write(io,111) ntip(i,1), x,y,z, 0.0
         else
            x = npos(ntip(i,2),1)
            y = npos(ntip(i,2),2)
            z = npos(ntip(i,2),3)
            phi = acosd(x/c)
            write(io,111) ntip(i,1), x,y,z, phi
         endif
      enddo
      close(io)
 111  format(t1,i7,4g15.7)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine warp3d_inp(job,jobh,etyp,elast,elnum,no_of_nodes,
     1           nea2,neb2, estk_s,estk_s0, estk_a, estk_b,estk_b0,
     2           nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c-- The routine generates an input deck for a FEM-analysis using WARP3D
c
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
      include 'mesh3d_scp_common_eln.f'
      integer  eln2o(iem),elo2n(iem)
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3), wblock_size
c
      integer  nea2,neb2, estk_s(2,1),estk_s0(2,1), estk_a(2,nea2,1),
     &         estk_b(2,neb2,1),estk_b0(2,1)
c
      integer  jobh,etyp, elnum,no_of_nodes, elast, iblock,ipart,
     2         uni,i,j,ii, ntip(400),njc,ialt
c
      integer  ndomj(101,1000),edomj(2,101,500),nnd,ned,
     &         jdomaut(0:48,201)
      integer  ivek1(20),ivek2(20),num1,num2,iend
c
      double precision  youngs,poisson,hard
c
      character job*40,filinp*40, el_list*8, answ*3,
     1          domain*10,dname*7, file_n2o*30, row*100,
     2          coordfile*40,elemfile*40,constfile*50,prdspfile*50,
     3          blckfile*40, ablck*60
c
      logical   ren_el,ok
c
      integer          ims,jms,kms,is,ima,jma,imb,jmb,kma
      common     /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
c  Check if the elements have been re-numbered
c
      write(*,'(t1,a,$)')
     &  '>> Have the elements been re-numbered using PATWARP ? (y/n): '
      read(*,'(a)') answ
      if ( (answ.eq.'Y').or.(answ.eq.'y') ) then
         ren_el = .true.
         write(*,'(t1,a,$)')
     &       '>> New to old element correspondence file name = '
         read(*,'(a)') file_n2o
         open(unit=13,file=file_n2o,status='unknown')
 10      continue
            read(13,'(a)',end=20) row
            call char_to_int(row,ivek1,num1)
            if (num1.ne.2) goto 21
            iend = ivek1(2)
c
            read(13,'(a)',end=21) row
            call char_to_int(row,ivek1,num1)
            read(13,'(a)',end=21) row
            call char_to_int(row,ivek2,num2)
            if (num1.ne.num2) goto 21
            do i=1, num1
               elo2n(ivek2(i)) = ivek1(i)
               eln2o(ivek1(i)) = ivek2(i)
            enddo
            write(77,'(t2,12i6)') ( ivek1(i),i=1,num1 )
            write(77,'(t2,12i6)') ( ivek2(i),i=1,num1 )
            write(77,*) ' *   *  *  *  *  *'
            goto 10
 20      continue
         close(13)
         if (iend.ne.elnum) then
            write(*,'(a,a,a)') 'No. of elements in file =',file_n2o,
     &        ' different from model!'
            stop
         endif
         goto 22
 21      write(*,'(a,a)') 'ERROR: in the file =',file_n2o
         close(13)
         stop
      else
         ren_el = .false.
      endif
 22   continue
c
c
c
      filinp=job(1:jobh)//'.warp_input'
      open(unit=14,file=filinp,status='unknown')
      uni=14
c
c  write out a general discription of the problem as a comment
c
	write(uni,'(4(t1,a/),t1,a)')
     & 'c ',
     & 'c  surface cracked plate:',
     & 'c ',
     & 'c  one quarter of the plate is modelled',
     & 'c '
c
	write(uni,8001) 'structure surface_crack'
        write(uni,8001) 'c'
c
        youngs = 30000.0
        poisson = 0.30
        hard = 10.0
c
        write(uni,8001) 'material steel'
        write(uni,9001) ' properties mises e ',youngs,' nu ',poisson,
     &   '  yld_pt 1.0e10  n_power ',hard,' ,',  'rho 0.0'
 8001	format(t1,a)
 8002	format(t1,a,i7)
 8003	format(t1,5a)
 9001	format(t1,a,f10.0,a,f5.3,a,f6.0,a/t20,a)
        write(uni,8001) 'c'
        write(uni,8002) 'number of nodes ',no_of_nodes
        write(uni,8002) 'number of elements',elnum
c
        write(uni,8001) 'c'
        write(uni,8001) 'coordinates'
        write(uni,8001) '*echo off'
        coordfile = job(1:jobh)//'.crd'
        write(uni,8003) '*input from ''',coordfile(1:(jobh+4)), ''''
c
        open(unit=48,file=coordfile(1:(jobh+4)),status='unknown')
        do i=1, inm
           if (nnr(i).gt.0) write(48,201) nnr(i),
     &                      npos(i,1),npos(i,2),npos(i,3)
        enddo
        close(48)
 201    format(t3,i6, 3(g20.10) )
c
        write(uni,8001) 'c'
        write(uni,8001) 'elements'
        if (elnum.lt.100) then
           write(el_list,'(t1,a,i2)') '1-', elnum
        elseif (elnum.lt.1000) then
           write(el_list,'(t1,a,i3)') '1-', elnum
        elseif (elnum.lt.10000) then
           write(el_list,'(t1,a,i4)') '1-', elnum
	endif
        write(uni,'(t4,a,t15,a)') el_list,
     &         'type l3disop   nonlinear  material steel  ,'
        write(uni,'(t15,a)')
     &         'order 2x2x2    bbar    center_output  short'
        write(uni,8001) 'c'
        write(uni,8001) 'incidences'
        elemfile = job(1:jobh)//'.elm'
        write(uni,8003) '*input from ''',elemfile(1:(jobh+4)),''''
c
        open(unit=49,file=elemfile(1:(jobh+4)),status='unknown')
        if (ren_el) then
           do i=1, elnum
              write(49,203) i, ( eln( eln2o(i) ,j), j=1, etyp)
           enddo
        else
           do i=1, elnum
              write(49,203) ( eln(i,j), j=0, etyp)
           enddo
        endif
        close(49)
 203    format( t1,9(i7) )
c
        write(uni,8001) 'c'
c
        if (ren_el) then
 40        continue
           write(*,'(t2,a,$)')
     &            '*Give input file with new element blocking : '
           read(*,'(a)') blckfile
           inquire(file=blckfile,exist=ok)
           if (.not.ok) then
            write(*,'(t4,3a)')'>>The file: ',blckfile,'doesn''t exist!'
              goto 40
           endif
           open(unit=46,file=blckfile,status='old')
 45        continue
           read(46,'(a)',end=50) ablck
           write(uni,'(a)') ablck
           goto 45
 50        continue
           close(46)
        else
          write(*,'(t1,a,$)')
     &     '>> Element block size (only scalar blocking supported): '
           read(*,*) wblock_size
           write(uni,8001) 'blocking $ scalar'
 	   iblock =0
	   ipart = mod(elnum,wblock_size)
	   do i= 1, elnum-ipart, wblock_size
	      iblock = iblock + 1
	      write(uni,'(t5,i4,t15,i3,t25,i6)') iblock, wblock_size, i
	   end do
           if (ipart.gt.0) then
              iblock = iblock + 1
              write(uni,'(t5,i4,t15,i3,t25,i6)')
     &                   iblock,ipart,elnum-ipart+1
         endif
      endif
c
      write(uni,8001) '*echo on'
      write(uni,8001) 'c'
      write(uni,8001) 'c        displacement increments'
      write(uni,8001) 'c'
 9010 format(t1,a,t37,a)
 9011 format(t1,a,t3,a,t18,a,t35,a,t55,a)
 9012 format(t1,a,t5,i2,t15,i4,a,i4,t30,g15.8,t50,g15.8  )
c
      write(uni,8001) 'c'
      write(uni,8001) 'c  define a dummy prescribed load step'
      write(uni,'(t1,a,a)') 'c  (necessary even though the ',
     &                        'loading is prescribed displacements)'
      write(uni,8001) 'c'
      write(uni,8001) 'loading dummy '
      write(uni,8001) '   nodal loads '
      write(uni,8001) '   2 force_x 0.0'
      write(uni,8001) 'c'
      write(uni,8001) 'loading predisp'
      write(uni,8001) '    nonlinear'
      write(uni,'(t1,a,i)') '    steps  1-100  dummy 1.0'
c
      write(uni,8001) 'c'
      write(uni,8001) 'constraints'
      write(uni,8001) '*echo off'
      write(uni,8001) 'c'
      write(uni,8001) 'c       symmetry plane and crack plane'
      write(uni,8001) 'c'
      constfile = job(1:jobh)//'_001.const'
      write(uni,8003) '*input from ''',constfile(1:(jobh+10)),''''
      prdspfile = job(1:jobh)//'_001.prdsp'
      write(uni,8001) 'c'
      write(uni,8001) 'c      imposed axial extension, remote end'
      write(uni,8001) 'c'
      write(uni,8003) '*input from ''',prdspfile(1:(jobh+10)),''''
      prdspfile = job(1:jobh)//'_002.prdsp'
      write(uni,8003) 'c  *input from ''',prdspfile(1:(jobh+10)),''''
      write(uni,8001) '*echo on'
      write(uni,8001) 'c'
      write(uni,8001) 'c    solution paramters'
      write(uni,8001) 'c'
      write(uni,8001) ' nonlinear analysis parameters'
      write(uni,8001) 'c   solution technique direct sparse '
      write(uni,8001) '   solution technique lnpcg'
      write(uni,8001) 'c   solution technique direct sparse '
      write(uni,8001) 'c   solution technique direct sparse hp'
      write(uni,8001) 'c   solution technique direct sparse bcs'
      write(uni,8001) 'c   solution technique direct sparse sgi'
      write(uni,8001) '   preconditioner type diagonal'
      write(uni,8001) '   lnr_pcg conv test res tol 0.01'
      write(uni,8001) '   maximum linear iterations 20000'
      write(uni,8001) '   maximum iterations 5'
      write(uni,8001) '   minimum iterations 2'
      write(uni,8001) '   convergence test norm res tol 0.01'
      write(uni,8001) '   nonconvergent solutions stop'
      write(uni,8001) '   adaptive on'
      write(uni,8001) '   linear stiffness for iteration one off'
      write(uni,8001) '   batch messages off'
      write(uni,8001) '   cpu time limit off'
      write(uni,8001) '   material messages off'
      write(uni,8001) '   bbar stabilization factor 0.0'
      write(uni,8001) '   consistent q-matrix on'
      write(uni,8001) '   extrapolation on'
      write(uni,8001) '   time step 1.0e06'
      write(uni,8001) 'c'
      write(uni,8001) 'c    start the analysis'
      write(uni,8001) 'c'
      write(uni,8001)
     &   ' compute displacements for load predisp steps 1'
c     &   ' compute displacements for load predisp steps 1-10'
      write(uni,8001) ' output patran binary displacements'
      write(uni,8001) ' output patran binary stresses'
      write(uni,8001) ' output patran binary strains'
      write(uni,8001) 'c'
      write(uni,8001) 'c     j-integral calculations'
      write(uni,8001) 'c'
      write(uni,8001) 'c   domain 001 is at point of max crack depth'
      write(uni,8001) 'c   last domain is at free surface'
      write(uni,8001) 'c'
c
c Call routine that define domains directly
c
      call warp3d_make_domain(ns1,ns2,ns3,nods,na1,na2,na3,noda,
     1                  nb1,nb2,nb3,nodb,etyp,elnum,no_of_nodes,
     2                  nnd,ned,ndomj,edomj )
c
      call def_aut_jdomain(ns1,ns2,ns3,nods,elast,jdomaut)
c
c
c      write(*,'(t1,a/t10,a/t10,a/t10,a,tr5,a,$)')
c     1     '* Domain definition alternative;',
c     2     '(1) automatic ','(2) two-el.-layer direct def.',
c     3     '(3) four-el.-layer direct def.',' >> Give alt. : '
c      read(*,*) ialt
       ialt = 1
c      if ((ialt.eq.3).and.(mod(jms-1,8).ne.0)) then
c         write(*,*) '>>If alt. = 3, then NA must be a multiple of 4'
c         stop
c      endif
c
      call nodset_ntip_(ntip,njc, nods,ns1,ns2,ns3)
      if (ialt.le.2) then
         do i=1, njc
            if (i.lt.10) then
               write(domain,'(t1,a,i1)') 'domain_00',i
               write(dname,'(t1,a,i1)') 'dnr_00',i
            elseif (i.lt.100) then
               write(domain,'(t1,a,i2)') 'domain_0',i
               write(dname,'(t1,a,i2)') 'dnr_0',i
            else
               write(domain,'(t1,a,i3)') 'domain_',i
               write(dname,'(t1,a,i3)') 'dnr_',i
            endif
            write(uni,8003) '*input from file ''',domain ,''''
            if (i.eq.1) then
               if (ialt.eq.1) then
                   call warp3d_domain(dname,domain,jdomaut,1,2,0,'a')
               else
                   call warp3d_q2_domain(dname,domain,
     &                  ntip(1),ntip(2),0,'a',ndomj,edomj,nnd,ned,1)
               endif
            elseif (i.eq.njc) then
               if (ialt.eq.1) then
                  call warp3d_domain(dname,domain,jdomaut,0,i-1,i,'c')
               else
                  call warp3d_q2_domain(dname,domain,
     &            0,ntip(njc-1),ntip(njc),'c',ndomj,edomj,nnd,ned,jms)
               endif
            else
               if (ialt.eq.1) then
                 call warp3d_domain(dname,domain,jdomaut,i-1,i,i+1,'b')
               else
                  call warp3d_q2_domain(dname,domain,ntip(i-1),
     &                 ntip(i),ntip(i+1),'b',ndomj,edomj,nnd,ned,2*i-1)
               endif
            endif
         enddo
      else
         ii = 0
         do i=3, njc-2, 4
            ii = ii + 1
            if (ii.lt.10) then
               write(domain,'(t1,a,i1)') 'domain_00',ii
               write(dname,'(t1,a,i1)') 'dnr_00',ii
            elseif (ii.lt.100) then
               write(domain,'(t1,a,i2)') 'domain_0',ii
               write(dname,'(t1,a,i2)') 'dnr_0',ii
            else
               write(domain,'(t1,a,i3)') 'domain_',ii
               write(dname,'(t1,a,i3)') 'dnr_',ii
            endif
            write(uni,8003) '*input from file ''',domain ,''''
            call warp3d_q4_domain(dname,domain,ntip(i-2),ntip(i-1),
     &      ntip(i),ntip(i+1),ntip(i+2),'b',ndomj,edomj,nnd,ned,2*i-1)
         enddo
      endif
      write(uni,8001) 'c'
      write(uni,8001) 'stop'
c
      close(uni)
c
c  Generate the constraint file
c
      call warp3d_const(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                  nodb,nb1,nb2,nb3, job,jobh,elast)
c
c  Generate the files containing the domain integral computations
c
      call warp3d_prdsp(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                  nodb,nb1,nb2,nb3, job,jobh,elast)
c
        call write_out_tipnodes(nods,ns1,ns2,ns3,job,jobh)
c 
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine warp3d_const(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &  nodb,nb1,nb2,nb3, job,jobh,elast)
c
c Routine determines the fixed BC.
c
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
c
      integer    nstm
      parameter (nstm=20000)
      integer    nfront(nstm),nsym(nstm),nconst(0:3,nstm)
c
      integer jobh, elast,i,j,cuni, n,n1,n2, node, iu,iv,iw
c
      double precision u,v,w
c
      character job*40,constfile*50
      logical   exist
c
      integer mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype
      integer sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
      common /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
c Node Set = NFRONT (Y=0)
c
      call nodset_nfront_(nfront,n1, nods,ns1,ns2,ns3,
     &             noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c Node Set = NSYM (X=0)
c
      call nodset_nsym_(nsym,n2,elast, nods,ns1,ns2,ns3,
     &                                 noda,na1,na2,na3)
c
c Sort out the unique nodes
c
      do i=1, n1
         nconst(0,i) = nfront(i)
         nconst(1,i) = 0
         nconst(2,i) = 1
         nconst(3,i) = 0
      enddo
      n = n1
      do j=1, n2
         exist = .false.
         do i=1, n1
            if ( nsym(j).eq.nconst(0,i) ) then
               nconst(1,i) = 1
               exist = .true.
            endif
         enddo
         if (.not.exist) then
            n = n + 1
            nconst(0,n) = nsym(j)
            nconst(1,n) = 1
            nconst(2,n) = 0
            nconst(3,n) = 0
         endif
      enddo
c
c  Write out the constraints on file:
c
      constfile = job(1:jobh)//'_001.const'
      cuni = 55
      open(unit=cuni,file=constfile(1:jobh+10),status='unknown')
      do i=1, n
         node = nconst(0,i)
         iu = nconst(1,i)
         iv = nconst(2,i)
         iw = nconst(3,i)
         u  = 0.0
         v  = 0.0
         w  = 0.0
         call fixed_constraint_warp3d(cuni,node,iu,iv,iw,u,v,w)
      enddo
      close(cuni)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine warp3d_prdsp(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &           nodb,nb1,nb2,nb3, job,jobh,elast)
c
c Routine determines the fixed BC.
c
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer    nstm
      parameter (nstm=20000)
      integer    nbak(nstm),nsym(nstm)
c
      integer jobh, elast,i,j,cuni, n,n2, node, iu,iv,iw
c
      double precision u,v,w, coord(3,nstm),zmin,zmax,zm,h,delta,theta
c
      character job*40,prdspfile*50
C
      logical   unique
c
      integer mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype
      integer sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
      common /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
c Node Set = NBAK ( Y = Lt )
c
        n = 0
c-Zone A:
	do j=1, jma-4
	   do i=1, ima
	      if (nnr(noda(i,j,kma)).gt.0) then
                 n = n + 1
                 nbak(n) = nnr(noda(i,j,kma))
                 coord(1,n) = npos(noda(i,j,kma),1)
                 coord(2,n) = npos(noda(i,j,kma),2)
                 coord(3,n) = npos(noda(i,j,kma),3)
	      endif
	   enddo
	enddo
	do j=jma-3, jma
	   do i=6, ima
	      if (nnr(noda(i,j,kma)).gt.0) then
                 n = n + 1
                 nbak(n) = nnr(noda(i,j,kma))
                 coord(1,n) = npos(noda(i,j,kma),1)
                 coord(2,n) = npos(noda(i,j,kma),2)
                 coord(3,n) = npos(noda(i,j,kma),3)
	      endif
	   enddo
	enddo
c-Zone B:
	do j=1, jmb
	   do i=2, imb
	      if (nnr(nodb(i,j,kma)).gt.0) then
                 n = n + 1
                 nbak(n) = nnr(nodb(i,j,kma))
                 coord(1,n) = npos(nodb(i,j,kma),1)
                 coord(2,n) = npos(nodb(i,j,kma),2)
                 coord(3,n) = npos(nodb(i,j,kma),3)
	      endif
	   enddo
	enddo
c
c Node Set = NSYM (X=0) (Needed to avoid specifying BC. twice)
c
      call nodset_nsym_(nsym,n2,elast, nods,ns1,ns2,ns3,
     &                                 noda,na1,na2,na3)
c
c  Write out the constraints on file:
c
      delta = 0.0001*coord(2,1)
c
c  Find Zmin and Zmax    =>    h = Zmax-Zmin
c
      zmin =  1.e+10
      zmax = -1.e+10
      do i=1, n
         zmax = max(zmax,coord(3,i))
         zmin = min(zmin,coord(3,i))
      enddo
c      write(*,*) '  ZMAX, ZMIN : ',ZMAX, ZMIN
      zm = (zmax-zmin)/2.
      h = abs(2*zm)
      theta = 2.*delta/h
c
c  Uniform displacement applied on the rear surface
c
      prdspfile = job(1:jobh)//'_001.prdsp'
      cuni = 55
      open(unit=cuni,file=prdspfile(1:jobh+10),status='unknown')
      do i=1, n
         node = nbak(i)
         unique = .true.
         do j=1, n2
            if (node.eq.nsym(j)) unique = .false.
         enddo
         if (unique) then
            iu = 1
            iv = 1
            iw = 1
         else
            iu = 0
            iv = 1
            iw = 1
         endif
         u  = 0.0
         v  = delta
         w  = 0.0
         call fixed_constraint_warp3d(cuni,node,iu,iv,iw,u,v,w)
      enddo
      close(cuni)
c
c  Linear varying displacement applied on the rear surface
c  (Correspondent to an angle theta)
c
      prdspfile = job(1:jobh)//'_002.prdsp'
      cuni = 55
      open(unit=cuni,file=prdspfile(1:jobh+10),status='unknown')
      do i=1, n
         node = nbak(i)
         unique = .true.
         do j=1, n2
            if (node.eq.nsym(j)) unique = .false.
         enddo
         if (unique) then
            iu = 1
            iv = 1
            iw = 1
         else
            iu = 0
            iv = 1
            iw = 1
         endif
         u  = 0.0
         v  = -theta * ( coord(3,i) - zm )
         w  = 0.0
         call fixed_constraint_warp3d(cuni,node,iu,iv,iw,u,v,w)
      enddo
      close(cuni)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine fixed_constraint_warp3d(cuni,node,iu,iv,iw,u,v,w)
c
c  Write out constraint on WARP3D-format on unit=CUNI
c
      implicit none
      integer cuni,node,iu,iv,iw
      double precision u,v,w
c
      if ( (iu.eq.1) .and. (iv.eq.0) .and. (iw.eq.0) ) then
         write(cuni,101) node,'  u ',u
      elseif ( (iu.eq.1) .and. (iv.eq.1) .and. (iw.eq.0) ) then
         write(cuni,102) node,'  u ',u,' v  ',v
      elseif ( (iu.eq.1) .and. (iv.eq.1) .and. (iw.eq.1) ) then
         write(cuni,103) node,'  u ',u,' v ',v,' w ',w
      elseif ( (iu.eq.0) .and. (iv.eq.1) .and. (iw.eq.0) ) then
         write(cuni,101) node,'  v ',v
      elseif ( (iu.eq.0) .and. (iv.eq.1) .and. (iw.eq.1) ) then
         write(cuni,102) node,'  v ',v,' w ',w
      elseif ( (iu.eq.0) .and. (iv.eq.0) .and. (iw.eq.1) ) then
         write(cuni,101) node,'  w ',w
      else
         write(*,*) '>> error in the constraint equations for warp3d'
         stop
      endif
 101  format(i6,a,g15.8)
 102  format(i6,2(a,g15.8))
 103  format(i6,3(a,g15.8))
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine def_aut_jdomain(ns1,ns2,ns3,nods,elast,jdomaut)
c
      implicit none
      include 'mesh3d_scp_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),elast,jdomaut(0:48,1)
      integer  n,i,j,k,ist,jj
c
      integer          ims,jms,kms,is,ima,jma,imb,jmb,kma
      common     /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
      if (elast.eq.1) then
         ist=(ims+1)/2+2
      else
         ist=3
      endif
c
      jj = 0
      do j=1, jms, 2
         jj = jj + 1
         k=1
         n=1
         jdomaut(n,jj) = nnr(nods(1,j,k))
         do i=ist, ims-4, 2
            n = n + 1
            jdomaut(n,jj) = nnr(nods(i,j,k))
         enddo
         jdomaut(0,jj) = n
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine warp3d_domain(dname,domain,jdomaut,n1,n2,n3,dtype)
c
c  Generates a file with domain integral commands
c
      implicit none
      integer  jdomaut(0:48,1),n1,n2,n3, io,nc
      character domain*10,dname*7,dtype*1,contours*5
c
      integer      mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
c
      nc = mr + min(m1,m2) - 2
      if (nc.lt.10) then
         write(contours,'(t1,a,i1)') '2-',nc
      elseif (nc.lt.100) then
         write(contours,'(t1,a,i2)') '2-',nc
      else
         write(contours,'(t1,a,i3)') '2-',nc
      endif
c
      io=49
      open(unit=io,file=domain,status='unknown')
c
      write(io,'(t1,a,a)') 'domain ',dname
      write(io,'(t4,a)')   'symmetric'
      write(io,'(t4,a)')   'normal plane nx 0  ny 0  ny -1.0'
c
      if (dtype.eq.'a') then
         call write_out_jnset(n1,jdomaut,io)
         call write_out_jnset(n2,jdomaut,io)
         write(io,'(t4,a,2i4,a)') 'front node sets',n1,n2,' linear'
         write(io,'(t4,a,a)')     'q-values automatic rings ',contours
         write(io,'(t4,a)')       'function type a'
      elseif (dtype.eq.'b') then
         call write_out_jnset(n1,jdomaut,io)
         call write_out_jnset(n2,jdomaut,io)
         call write_out_jnset(n3,jdomaut,io)
         write(io,'(t4,a,3i4,a)')'front node sets',n1,n2,n3,' linear'
         write(io,'(t4,a,a)')     'q-values automatic rings ',contours
         write(io,'(t4,a)')       'function type b'
      elseif (dtype.eq.'c') then
         call write_out_jnset(n2,jdomaut,io)
         call write_out_jnset(n3,jdomaut,io)
         write(io,'(t4,a,2i4,a)') 'front node sets',n2,n3,' linear'
         write(io,'(t4,a,a)')     'q-values automatic rings ',contours
         write(io,'(t4,a)')       'function type c'
      endif
      write(io,'(t4,a)')     'print totals'
      write(io,'(t4,a)')     'compute domain integral'
c
      close(io)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine write_out_jnset(n,jdomaut,io)
      integer n,jdomaut(0:48,1),io
      if (n.le.8) then
         write(io,108) 'node set ',n, (jdomaut(i,n),i=1,jdomaut(0,n))
      elseif (n.le.18) then
         write(io,118) 'node set ',n, (jdomaut(i,n),i=1,jdomaut(0,n))
      elseif (n.le.28) then
         write(io,128) 'node set ',n, (jdomaut(i,n),i=1,jdomaut(0,n))
      elseif (n.le.38) then
         write(io,138) 'node set ',n, (jdomaut(i,n),i=1,jdomaut(0,n))
      elseif (n.le.48) then
         write(io,148) 'node set ',n, (jdomaut(i,n),i=1,jdomaut(0,n))
      else
         write(io,'(t1,a,a)') 'ERROR: max circumferential nodes ',
     &    'in automatic J-domain routines = 48'
         write(io,'(t1,a,a)') 'Increase size of array jdomaut(0:48,*)'
      endif
c
 108  format(t4,a,i2,tr1,8i6)
 118  format(t4,a,i2,tr1,8i6,a/t4,10i6)
 128  format(t4,a,i2,tr1,8i6,a/t4,10i6,a/t4,10i6)
 138  format(t4,a,i2,tr1,8i6,a/t4,10i6,a/t4,10i6,a/t4,10i6)
 148  format(t4,a,i2,tr1,8i6,a/t4,10i6,a/t4,10i6,a/t4,10i6,a/t4,10i6)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine warp3d_q2_domain(dname,domain,n1,n2,n3,dtype,
     &                            ndomj,edomj,nnd,ned,jst)
c
c  Generates a file with domain integral commands
c
      implicit none
      integer  n1,n2,n3,ndomj(101,1),edomj(2,101,1),nnd,ned,jst
      integer  io,nv(1000),ev(1000),i,j,k,ipart,nednew
      character domain*10,dname*7,dtype*1, form*14
c
      integer mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
c
      io=49
      open(unit=io,file=domain,status='unknown')
c
      write(io,'(t1,a,a)') 'domain ',dname
      write(io,'(t2,a)')   'symmetric'
      write(io,'(t2,a)')   'normal plane nx 0  ny 0  ny -1.0'
c
      k=0
      if (dtype.eq.'a') then
         write(io,'(t2,a,2i7,a)') 'front nodes',n1,n2,' linear verify'
         write(io,'(t2,a)')       'function type a'
         do j=1, ned
            k=k+1
            ev(k) = edomj(2,jst,j)
         enddo
      elseif (dtype.eq.'b') then
         write(io,'(t2,a,3i7,a)')'front nodes',n1,n2,n3,' linear verify'
         write(io,'(t2,a)')       'function type b'
         do j=1, ned
            k=k+1
            ev(k) = edomj(1,jst,j)
         enddo
         do j=1, ned
            k=k+1
            ev(k) = edomj(2,jst,j)
         enddo
         nednew = 2*ned
      elseif (dtype.eq.'c') then
         write(io,'(t2,a,2i7,a)') 'front nodes',n2,n3,' linear verify'
         write(io,'(t2,a)')       'function type c'
         do j=1, ned
            k=k+1
            ev(k) = edomj(1,jst,j)
         enddo
      endif
c
      do i=1, nnd
         nv(i) = ndomj(jst,i)
      enddo
      ipart = mod(nnd,8)
      if (mod(nnd,8).eq.0) ipart = 8
      write(io,'(t2,a,t12,8i7,'','')') 'q-values', (nv(i),i=1,8)
      write(io,'(t12,8i7,'','')') (nv(i),i=9, nnd-ipart)
      write(form,'(a,i1,a)') '(t12,',ipart,'i7,f6.3)'
      write(io,form) (nv(i),i=nnd-ipart+1, nnd), 1.00
c
      ipart = mod(nednew,10)
      if (mod(nednew,10).eq.0) ipart = 10
      write(io,'(t2,a,t12,10i6,'','')') 'elements', (ev(i),i=1,10)
      write(io,'(t12,10i6,'','')') (ev(i),i=11, nednew-ipart)
      write(io,'(t12,10i6)')     (ev(i),i=nednew-ipart+1, nednew)
c
      write(io,'(t2,a)')     'print totals'
      write(io,'(t2,a)')     'compute domain integral'
c
      close(io)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine warp3d_q4_domain(dname,domain,n1,n2,n3,n4,n5,dtype,
     &                   ndomj,edomj,nnd,ned,jst)
c
c  Generates a file with domain integral commands
c
      implicit none
      integer  n1,n2,n3,n4,n5,ndomj(101,1),edomj(2,101,1),nnd,ned,jst
      integer  io,nv(1000),ev(1000),i,j,k,ipart,nednew
      character domain*10,dname*7,dtype*1, form*14
c
      integer mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
c
      io=49
      open(unit=io,file=domain,status='unknown')
c
      write(io,'(t1,a,a)') 'domain ',dname
      write(io,'(t2,a)')   'symmetric'
      write(io,'(t2,a)')   'normal plane nx 0  ny 0  ny -1.0'
c
      if (dtype.ne.'b') then
         write(*,'(t1,a,a)') 'The trial function should be',
     &                 ' equal to b! check the program'
         stop
      endif
      write(io,'(t2,a,5i7,a)') 'front nodes',n1,n2,n3,n4,n5,
     &                                     ' linear verify'
      write(io,'(t2,a)')       'function type b'
c
      k=0
      do j=1, ned
         k=k+1
         ev(k) = edomj(1,jst-2,j)
      enddo
      do j=1, ned
         k=k+1
         ev(k) = edomj(1,jst,j)
      enddo
      do j=1, ned
         k=k+1
         ev(k) = edomj(2,jst,j)
      enddo
      do j=1, ned
         k=k+1
         ev(k) = edomj(2,jst+2,j)
      enddo
      nednew = k
c
      do i=1, nnd
         nv(i) = ndomj(jst-2,i)
      enddo
      ipart = mod(nnd,8)
      if (mod(nnd,8).eq.0) ipart = 8
      write(io,'(t2,a,t12,8i7,'','')') 'q-values', (nv(i),i=1,8)
      write(io,'(t12,8i7,'','')') (nv(i),i=9, nnd-ipart)
      write(form,'(a,i1,a)') '(t12,',ipart,'i7,f6.3)'
      write(io,form) (nv(i),i=nnd-ipart+1, nnd), 0.5
c
      do i=1, nnd
         nv(i) = ndomj(jst,i)
      enddo
      ipart = mod(nnd,8)
      if (mod(nnd,8).eq.0) ipart = 8
      write(io,'(t2,a,t12,8i7,'','')') 'q-values', (nv(i),i=1,8)
      write(io,'(t12,8i7,'','')') (nv(i),i=9, nnd-ipart)
      write(form,'(a,i1,a)') '(t12,',ipart,'i7,f6.3)'
      write(io,form) (nv(i),i=nnd-ipart+1, nnd), 1.00
c
      do i=1, nnd
         nv(i) = ndomj(jst+2,i)
      enddo
      ipart = mod(nnd,8)
      if (mod(nnd,8).eq.0) ipart = 8
      write(io,'(t2,a,t12,8i7,'','')') 'q-values', (nv(i),i=1,8)
      write(io,'(t12,8i7,'','')') (nv(i),i=9, nnd-ipart)
      write(form,'(a,i1,a)') '(t12,',ipart,'i7,f6.3)'
      write(io,form) (nv(i),i=nnd-ipart+1, nnd), 0.50
c
      ipart = mod(nednew,10)
      if (mod(nednew,10).eq.0) ipart = 10
      write(io,'(t2,a,t12,10i6,'','')') 'elements', (ev(i),i=1,10)
      write(io,'(t12,10i6,'','')') (ev(i),i=11, nednew-ipart)
      write(io,'(t12,10i6)')     (ev(i),i=nednew-ipart+1, nednew)
c
      write(io,'(t2,a)')     'print totals'
      write(io,'(t2,a)')     'compute domain integral'
c
      close(io)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine warp3d_make_domain(ns1,ns2,ns3,nods,na1,na2,na3,noda,
     1           nb1,nb2,nb3,nodb,etyp,elnum,no_of_nodes,
     2           nnd,ned,ndomj,edomj )
c
c  Generates a file with domain integral commands
c
      implicit none
c
      include 'mesh3d_scp_common_eln.f'
      include 'mesh3d_scp_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  etyp,elnum,no_of_nodes, i,j,jj,k,kk,ie,in,
     1         nnd,ned,nno,eno,n1,n2,ndomj(101,1), edomj(2,101,1)
c     2         ,ndomj_coord(101,1)
c
      logical  unique,ok1,ok2
c
      integer          ims,jms,kms,is,ima,jma,imb,jmb,kma
      common     /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
      integer mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
c

c
c  (1) find all the nodes in each node layer of interest
c      ( Only one domain in each layer)
c
      do j=1, jms,2
         jj = 0
         do k=1, kms-2
            do i=1, ims
               unique = .true.
               nno = nnr(nods(i,j,k))
               if (nno.gt.0) then
                  do kk=1, jj
                     if (nno.eq.ndomj(j,kk)) unique = .false.
                  enddo
                  if (unique) then
                     jj = jj + 1
                     ndomj(j,jj) = nno
c                     ndomj_coord(j,jj) = nods(i,j,k)
                  endif
               endif
            enddo
         enddo
c         write(*,'(t1,a,i3,a,i5)') '>> j=',j,' nnd = ',jj
c         write(56,'(t1,a,i3,a,i5)') '>> j=',j,' nnd = ',jj
c         write(56,'(t5,10i6)') (ndomj(j,kk), kk=1, jj)
      enddo
      nnd = jj
c
c  (2) find the q-values for a four-element-layer domain
c      (As a fisrt approximation use q: 0.5  1.0  0.5 in each domain)
c
c      do j=5, jms-4, 8
c         do i=1, nnd
c            do k=1, 3
c               x1(k) = npos(ndomj_coord(j-4,i),k)
c               x2(k) = npos(ndomj_coord(j-2,i),k)
c               x3(k) = npos(ndomj_coord(j,i),k)
c               x4(k) = npos(ndomj_coord(j+2,i),k)
c               x5(k) = npos(ndomj_coord(j+4,i),k)
c            enddo
c            xl1 = sqrt( (x1(1)-x3(1))**2. + (x1(2)-x3(2))**2. +
c     &                  (x1(3)-x3(3))**2. )
c            xl2 = sqrt( (x2(1)-x3(1))**2. + (x2(2)-x3(2))**2. +
c     &                  (x2(3)-x3(3))**2. )
c            xl4 = sqrt( (x4(1)-x3(1))**2. + (x4(2)-x3(2))**2. +
c     &                  (x4(3)-x3(3))**2. )
c            xl5 = sqrt( (x5(1)-x3(1))**2. + (x5(2)-x3(2))**2. +
c     &                  (x5(3)-x3(3))**2. )
c            q1 = xl2/xl1
c            q2 = 1.0
c            q3 = xl4/xl5
c
c
c  (3) find all elements involved in the domains.
c
      do j=1, jms-2, 2
         jj = 0
         do i=1, nnd
            n1 = ndomj(j,i)
            n2 = ndomj(j+2,i)
            do ie=1, elnum
               eno = eln(ie,0)
               ok1 = .false.
               ok2 = .false.
               do in=1, etyp
                  if (eln(ie,in).eq.n1) then
                     ok1 = .true.
                  elseif (eln(ie,in).eq.n2) then
                     ok2 = .true.
                  endif
               enddo
               if (ok1.and.ok2) then
                  unique = .true.
                  do kk=1, jj
                     if (eno.eq.edomj(2,j,kk)) unique = .false.
                  enddo
                  if (unique) then
                     jj = jj + 1
                     edomj(2,j,jj) = eno
                  endif
c                  goto 10
               endif
            enddo
c
c  10        continue
c
         enddo
c         write(*,'(t1,a,i3,a,i5)') '>> j=',j,' ned = ',jj
c         write(56,'(t1,a,i3,a,i5)') '>> j=',j,' ned = ',jj
c         write(56,'(t5,10i6)') (edomj(2,j,kk), kk=1, jj)
      enddo
      ned = jj
c
      do i=1, ned
         edomj(1,1,i) = 0
         edomj(2,jms,i) = 0
      enddo
      do j=3, jms, 2
         do i=1, ned
            edomj(1,j,i) = edomj(2,j-2,i)
         enddo
      enddo
c
c      write(*,*) '>> nnd =',nnd,' ned =',ned
c
c      do j=1, jms, 2
c         write(56,'(t1,a,i3,a,i5)') '>> j=',j,' ned = ',jj
c         write(56,'(t5,10i6)') (edomj(1,j,kk), kk=1, ned)
c         write(56,'(t5,10i6)') (edomj(2,j,kk), kk=1, ned)
c      enddo
c
c----67      
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine patran_neu(job,jobh,etyp,elast,elnum,no_of_nodes,
     1           nea2,neb2, estk_s,estk_s0, estk_a, estk_b,estk_b0,
     2           nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
      implicit none
c
      include 'mesh3d_scp_common_eln.f'
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  nea2,neb2, estk_s(2,1),estk_s0(2,1), estk_a(2,nea2,1),
     &         estk_b(2,neb2,1),estk_b0(2,1)
c
      integer  jobh,etyp,elast,elnum,no_of_nodes, uni,i,j
      integer  node,ndof,el_no,packet,inf_lines
      integer  igtype,tyel, head(15),versn(3)
c
      double precision x,y,z
c
      character job*40
c
c... Generate an input file on PATRAN format
c
      open(unit=14,file='patran.out.1',status='unknown')
      uni = 14
c
c  PACKET 25   =>  Head
c
      data  head /4htest,4h mod,4hel f,4hor w,4harp3,
     &            4hsurf,4hace ,4hcrac,4hked ,4hplat,
     &            4he   ,4h    ,4h    ,4h    ,4h    /
      packet = 25
      inf_lines = 1
      write(uni,9251) packet, 0, 0, inf_lines, 0, 0, 0, 0, 0
      write(uni,9252) head
 9251 format(i2,8i8)
 9252 format(15a4)
c
c  PACKET 26   =>  Global model info
c
      data  versn /4h    ,4h 2.2,4h-2  /
      packet = 26
      inf_lines = 1
      write(uni,9261) packet,0,0,inf_lines, no_of_nodes,elnum, 0,0,0
      write(uni,9262) '01-jan-96',' 00:00pm',versn
 9261 format(i2,8i8)
 9262 format(a9,3x,a8,3a4)
c
c  PACKET  1   =>  Node coordinstes
c
c       Loop over all nodes:
c            ( ndof = number of degrees of freedom;        )
c            ( Allways give three coordinates, i.e. x,y,z, )
c            ( e.g. a 2D problem  x, y, 0.0                )
c
      data igtype /1hg/
      packet = 1
      ndof=3
      inf_lines = 2
      do i=1, inm
            if (nnr(i).gt.0) then
            node = nnr(i)
            x = npos(i,1)
            y = npos(i,2)
            z = npos(i,3)
            write(uni,9011) packet, node, 0, inf_lines, 0,0,0,0,0
            write(uni,9012) x, y, z
            write(uni,9013) 1,igtype,ndof, 0,0,  0,0,0,0,0,0
         endif
      enddo
 9011 format( i2, 8i8 )
 9012 format( 3e16.9 )
 9013 format( i1, a1, 3i8, 2x, 6i1 )
c
c  PACKET  8   =>  Node constraints
c
c  Define the type of constrait for each node
c
c      Loop over all nodes:
c
c      packet = 8
c      inf_lines = 2
c      do i=1, inn
c         if ( (icn(1,i).eq.1).or.(icn(2,i).eq.1).or.
c     &          (icn(3,i).eq.1) ) then
c            node = nodn(i)
c            cflag(1) = xcn(1,i)
c            cflag(2) = xcn(2,i)
c            cflag(3) = xcn(3,i)
c            write(uni,9081) packet, node, 1, inf_lines, 0,0,0,0,0
c            write(uni,9082) 0, (cflag(j),j=1, 3),  0, 0, 0
c            write(uni,9083) xcn(1,i), xcn(2,i), xcn(3,i)
c         endif
c      enddo
c 9081 format( i2, 8i8 )
c 9082 format( i8, 6i1 )
c 9083 format( 3e16.9 )
c
c  PACKET  2   =>  Element incidences
c
c   etyp=nodes per elemnet;
c
c  Type of element is given by tyel according to the following rule:\
c     tyel:  2 = Bar          3 = Triangle    4 = Quadrilateral
c            5 = tetrahedral  6 = pyramid     7 = wedge
c            8 = hexahedral
c
      packet = 2
      if (etyp.le.10) then
         inf_lines = 2
      elseif (etyp.le.20) then
         inf_lines = 3
      else
         inf_lines = 4
      endif
      tyel = 8
      etyp = 8
      do i=1, elnum
         el_no = eln(i,0)
         write(uni,9021) packet, el_no, tyel, inf_lines, 0,0,0,0,0
         write(uni,9022) etyp, 0,0,0, 0.0, 0.0, 0.0
         write(uni,9023) ( eln(i,j), j=1, etyp )
      enddo
 9021 format( i2, 8i8 )
 9022 format( 4i8, 3e16.9 )
 9023 format( 10i8 )
c
c  PACKET 99   =>  Termination card
c
      packet = 99
      inf_lines = 0
      write(uni,9991) packet, 0,0, inf_lines, 0,0,0,0,0
 9991 format( i2, 8i8 )
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine mesh_statestik(job,jobh,nstk,y,z,rs,natyp,omega,
     &        etyp,zb_xw,rate_b,head,el_ch_zs,slice,prog,
     &        nea2,neb2, estk_s,estk_s0,estk_a,estk_b,estk_b0,
     &        eln_d,elnum,dsl,grp, no_of_jcont )
C--------------------------------------------------------------*
C       Rutinen skapar en statestikfil innehallande            *
C       information om modellen.                               *
C--------------------------------------------------------------*
	implicit none
c
      include 'mesh3d_scp_common_eln.f'
      integer  eln_d(0:iem)
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  nea2,neb2, estk_s(2,1),estk_s0(2,1), estk_a(2,nea2,1),
     &         estk_b(2,neb2,1),estk_b0(2,1)
c
      integer  jobh,nstk(0:500,7),natyp,etyp,el_ch_zs,
     &         elnum,i,k,sum,uni,zb_xw,slice,dsl(200),grp, no_of_jcont
C COMMON INTEGERS!
      integer   mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     &          sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma
 
	real*8    y(0:200),z(0:200,2),rs(200),omega,rate_b
C COMMON REAL*8!
	real*8    t,w,c,a,kappa,alfa,r1,r2,rn,eta,my
 
	character head*45,job*40,prog*20,fil4*40
 
	common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &         /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma,
     &         /geom/t,w,c,a,kappa,alfa,r1,r2,rn,eta,my
 
	fil4=job(1:jobh)//'_mesh.sta'
	open(unit=21,file=fil4,status='unknown')
	uni=21
 
	write(uni,'(/t10,a//t10,a)')
     &  '* * * *  m o d e l - i n f o r m a t i o n   * * * *',
     &  'the model is generated with the program mesh3d_scp'
	write(uni,'(/t5,a,a)') 'program-system : ',prog
	write(uni,'(/t5,a,a)') 'job       : ',job(1:jobh)
	write(uni,'(t5,a,a)')  'date      : ',head(15:20)
	write(uni,'(t5,a,f10.5,3(a,f10.5))')
     &  'geometri  : t=' ,t, '  w=' ,w, '  c=' ,c, '  a=' ,a
	write(uni,'(t5,a,f10.5,4(a,f12.7)/t17,5(a,i2) )')
     &  'zone s    : r1=',r1,' r2=',r2,' rn=',rn,' eta=',eta,
     &  ' my=',my,'mr=',mr,' mf=',mf,' mv=',mv,' mh=',mh,
     &  ' sfred=',sfred
        write(uni,'(t5,a,t15,a,t25,a,t35,a,t45,a,t55,a,t65,a)')
     &  'radius    :','r1','r2','r3','r4','r5','r6'
	if ((mr-1).le.6) then
	   write(uni,'(t10,6(tr1,g9.3))') ( rs(i),i=1, mr-1 )
	elseif ((mr-1).le.12) then
	   write(uni,'(t10,6(tr1,g9.3))') ( rs(i),i=1, 6 )
           write(uni,'(t15,a,t25,a,t35,a,t45,a,t55,a,t65,a)')
     &     'r7','r8','r9','r10','r11','r12'
	   write(uni,'(t10,6(tr1,g9.3))') ( rs(i),i=7, mr-1 )
	elseif ((mr-1).le.18) then
	   write(uni,'(t10,6(tr1,g9.3))') ( rs(i),i=1, 6 )
           write(uni,'(t15,a,t25,a,t35,a,t45,a,t55,a,t65,a)')
     &     'r7','r8','r9','r10','r11','r12'
	   write(uni,'(t10,6(tr1,g9.3))') ( rs(i),i=7, 12 )
           write(uni,'(t15,a,t25,a,t35,a,t45,a,t55,a,t65,a)')
     &     'r13','r14','r15','r16','r17','r18'
	   write(uni,'(t10,6(tr1,g9.3))') ( rs(i),i=13, mr-1 )
	endif
	write(uni,'(t5,4(a,i2),a,f8.4)')  'zone a    : m1=',m1,
     &  '  m2=',m2,'  na=',na,'  natyp=',natyp,' ,omega=',omega
	write(uni,'(t13,3(a,i2))')
     &          '    rtype=',rtype,'  lt=',lt,'  lred=',lred
	write(uni,'(t5,2(a,i2))')'zone b    : mb=',mb,'  nb=',nb
        write(uni,'(t5,a,i2,a,f5.3))')
     &       ': zb_xw=',zb_xw,'  rate_b=',rate_b
	write(uni,'(t5,3(a,i2))')  'mesh refinements : el_ch_zs=',
     &            el_ch_zs, ' slice=',slice
	write(uni,'(t5,a,i2,a)') 'elementype: ',etyp,' -noded'
	write(uni,'(t5,a,i2)') 'no_of_jcont= ',no_of_jcont
 
	write(uni,'(//t5,a/)')
     &  '-------- n o d - i n f o r m a t i o n -----------'
        write(uni,'(t5,a,i5,a,i5/t12,a,i5/t12,a,i5/t12,a,i5/)')
     &  'zone s : the range of nodenumber : ',nstk(0,3),' - ',
     &   nstk(0,4),': number of nodes (j=odd  no) =',nstk(0,1),
     &  ': number of nodes (j=even no) =',nstk(0,2),
     &  ': total number of nodes in s  =',nstk(0,7)
	write(uni,'(t5,a,t15,a,t36,a/t6,a,t14,a,t35,a,t60,a/)')
     &  'nodeset','zone a :','zone b :','k  (l)',
     &  '1''st  -  last  suma ','1''st  -  last  sumb ','totsum'
	do k=1, kma
	   if (mod(k,2).eq.1) then
              write(uni,'(t5,i2,t14,2(i5,a,i5,i6,a),t60,i6)') k,
     &        nstk(k,1),' - ',nstk(k,2),nstk(k,3),', ',
     &        nstk(k,4),' - ',nstk(k,5),nstk(k,6),', ',nstk(k,7)
	   else
              write(uni,'(t5,i2,t9,i2,t14,2(i5,a,i5,i6,a),t60,i6)') k,
     &        k/2,nstk(k,1),' - ',nstk(k,2),nstk(k,3),', ',
     &        nstk(k,4),' - ',nstk(k,5),nstk(k,6),', ',nstk(k,7)
	   endif
	enddo
	sum=0
	do k=0,kma
	   sum=sum+nstk(k,7)
	enddo
	write(uni,'(/t5,a,i5)')'the total number of nodes    = ',sum
 
	write(uni,'(/t5,a)')
     &  '-------- e l e m e n t - i n f o r m a t i o n -----------'
	write(uni,'(/t5,a,i5)') 'total number of elements = ',elnum
	write(uni,'(/t5,a)') 'element layers/thickness :'
	write(uni,'(t5,a,t20,a,t40,a)')
     &          'elem. layer','y','dy (thickness)'
	write(uni,'(t10,i2,t17,e16.8,t37,a)') 0, y(0), '  -'
	do i=1, lt
	   write(uni,'(t10,i2,t17,e16.8,t37,e16.8)') i,y(i),y(i)-y(i-1)
	enddo
	write(uni,'(/t5,a,i2)') 'number of element groups = ',grp
	write(uni,'(/t5,a,t17,a)') 'group no.',' no. of el. in group'
	do i=1, grp
	   write(uni,'(t8,i2,t22,i4)') i,dsl(i)
	enddo
	call write_el_stat(uni, el_ch_zs,eln_d,elnum,
     &       nea2,neb2,estk_s,estk_s0,estk_a,estk_b,estk_b0)
	close(uni)
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine model_info(eln_d,no_of_jcont,job,jobh, estk_a,nea2,
     &           nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c-- The routine creates node and element sets for post processing
c-- of the results.
c
      implicit none
c
      include 'mesh3d_scp_common_eln.f'
      integer eln_d(0:iem)
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  nea2,estk_a(2,nea2,1),jobh, no_of_jcont,
     &         i,i1,i2,j,k,kk,n,el_ntf(20,20,20,2),iv(0:200),uni
c
c COMMON INTEGER
      integer  mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &         ims,jms,kms,is,ima,jma,imb,jmb,kma
c
      character job*40,info_file*40
 
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
      common /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
      info_file=job(1:jobh)//'_set.sta'
      uni=33
      open(unit=uni,file=info_file,status='unknown')
 
C... Reference node numbers:
	write(uni,'(/t2,a)') '* reference nodes / node no: x y z'
	write(uni,'(t5,i6,3(g14.6))')  nnr(nodb(imb,jmb,1)),
     &      ( npos(nodb(imb,jmb,1),k), k=1, 3 )
	write(uni,'(t5,i6,3(g14.6))')  nnr(nodb(imb,1,1)),
     &      ( npos(nodb(imb,1,1),k), k=1, 3 )
	write(uni,'(t5,i6,3(g14.6))') nnr(nodb(imb,jmb,kma)),
     &      ( npos(nodb(imb,jmb,kma),k), k=1, 3 )
	write(uni,'(t5,i6,3(g14.6))') nnr(nodb(imb,1,kma)),
     &      ( npos(nodb(imb,1,kma),k), k=1, 3 )
	write(uni,'(t10,4i6)')  nnr(nodb(imb,jmb,1)),
     &   nnr(nodb(imb,1,1)),nnr(nodb(imb,jmb,kma)),nnr(nodb(imb,1,kma))
 
C... CMOD & Delta node numbers:
	write(uni,'(/t2,a)') '* cmod nodes / node no: x y z'
	do j=1, jma-6, 2
	   write(uni,'(t5,i6,3(g14.6))') nnr(noda(1,j,1)),
     &      ( npos(noda(1,j,1),k), k=1, 3 )
	enddo
 
	write(uni,'(/t4,a)')'* delta_f & delta_b nodes/node no: x y z'
	do k=kma-8, kma
	 write(uni,'(/t5,a,i2,a,g14.6)')'k=',k,' y=',npos(noda(1,1,k),2)
           write(uni,'(t10,a,g14.6)')  ' z=',npos(noda(1,1,k),3)
	   do j=1, 2*(na-nb)+1
	      if (nnr(noda(1,j,k)).gt.0) write(uni,'(t5,i6,a,g14.6)')
     &              nnr(noda(1,j,k)), ' x=',npos(noda(1,j,k),1)
	   enddo
           write(uni,'(t10,a,g14.6)')  ' z=',npos(noda(ima,1,k),3)
	   do j=1, 2*(na-nb)+1
	      if (nnr(noda(ima,j,k)).gt.0) write(uni,'(t5,i6,a,g14.6)')
     &            nnr(noda(ima,j,k)), ' x=',npos(noda(ima,j,k),1)
	   enddo
	enddo
 
C... Crack front nodes:
	write(uni,'(/t2,a)') '* crack front nodes / node no:'
	do j=1, jms
	   iv(j)=nnr(nods(1,j,1))
	enddo
	iv(0)=jms
	call write_no_el_info(uni,iv)
 
C... Nodes for calculation of radial length to a virtual shift or
C    to a domain ABAQUS.
C    ( Note that this works only for SFRED=0)
	write(uni,'(/t2,a)') '* vshift nodes (radial length) / node no:'
	do j=1, jms
	   write(uni,'(t10,a,i2)') 'j=',j
	   n=0
	   do k=1, kms, 2
	      n=n+1
	      iv(n)=nnr(nods(1,j,k))
	   enddo
	   i1=2*(m1+mh+1)+1
	   i2=i1+2*(no_of_jcont-n)
	   do i=i1, i2, 2
	      n=n+1
	      iv(n)=nnr(noda(i,j,1))
	   enddo
	   iv(0)=n
	   write(uni,'(t10,i3,i4)') j, n
	   call write_no_el_info(uni,iv)
	enddo
 
C... Nodes for computing the CTOD:
	write(uni,'(/t2,a)')'* crack surface nodes (ctod-computations)'
	do j=1, jms
	   write(uni,'(t10,a,i2)') 'j=',j
	   n=0
	   do i=1, ims-2
	      if (nnr(nods(i,j,1)).gt.0) then
	         n=n+1
	         iv(n)=nnr(nods(i,j,1))
	      endif
	   enddo
	   do k=4, kms
	      if (nnr(nods(ims,j,k)).gt.0) then
	         n=n+1
	         iv(n)=nnr(nods(ims,j,k))
	      endif
	   enddo
	   iv(0)=n
	   write(uni,'(t10,i2,i4)') j,n
	   call write_no_el_info(uni,iv)
	enddo
 
C... Elements for studying the near-tip fields:
C    ( Note that this works only for SFRED=0)
C . . ZON S.
	write(uni,'(/t2,a,i3,a,i3)')
     &   '* near-tip elements, no. of slices=',jms,' mr =',mr
	n=0
	do j=1, na
	   do i=1, mf-1
	      n=n+1
	      el_ntf(i,j,1,1)=eln(n,0)
	      el_ntf(i,j,1,2)=eln_d(n)
	   enddo
	   do k=2, mr
	      do i=1, mf
	         n=n+1
	         el_ntf(i,j,k,1)=eln(n,0)
	         el_ntf(i,j,k,2)=eln_d(n)
	      enddo
	   enddo
	enddo
C . . Xtra elements in zon A ( I=1 )
	do j=1, na
	   do k=1, m2-mh
	      i=estk_a(2,1,j)-(m2-mh)+k
	      el_ntf(1,j,mr+k,1)=eln(i,0)
	      el_ntf(1,j,mr+k,2)=eln_d(i)
	   enddo
	enddo
C... Write out the elements on the info file:
	do j=1, na
	   write(uni,'(t1,a,i3)') '=> slice nr. ',j
	   do k=1, mr+m2-mh
	      iv(2*k-1)=el_ntf(1,j,k,1)
	      iv(2*k)=el_ntf(1,j,k,2)
	   enddo
	   iv(0)=2*(mr+m2-mh)
	   write(uni,'(t5,3i5)') j, 1, iv(0)
	   call write_no_el_info(uni,iv)
	   do i=2, mf-1
	      do k=1, mr
	         iv(2*k-1)=el_ntf(i,j,k,1)
	         iv(2*k)=el_ntf(i,j,k,2)
	      enddo
	      iv(0)=2*mr
	      write(uni,'(t5,3i5)') j, i, iv(0)
	      call write_no_el_info(uni,iv)
	   enddo
	   do k=2, mr
	      kk=k-1
	      iv(2*kk-1)=el_ntf(mf,j,k,1)
	      iv(2*kk)=el_ntf(mf,j,k,2)
	   enddo
	   iv(0)=2*(mr-1)
	   write(uni,'(t5,3i5)') j, mf, iv(0)
	   call write_no_el_info(uni,iv)
	enddo
	close(uni)
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine write_no_el_info(uni,iv)
	implicit none
	integer uni,iv(0:200),i,i2,i3,j,num_per_row,row
	num_per_row=10
	row=int(iv(0)/num_per_row)
	i2=row*num_per_row
	i3=mod(iv(0),num_per_row)
	do i=1, i2, num_per_row
	   write(uni,'(t5,10i6)') ( iv(j), j=i, i+num_per_row-1 )
	enddo
	if (i3.gt.0) write(uni,'(t5,10i6)') ( iv(j), j=i2+1, i2+i3 )
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
