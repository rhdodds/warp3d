c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	integer function nodn(zon,i,j,k)
C-------------------------------------------------------
C  Funktionen r{knar ut ett nodnummer f|r en NOD
C-------------------------------------------------------
	implicit none
	integer mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     &          sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma,
     &          i,j,k,kstart,kstep
	character zon*1
 
	common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &         /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma,
     &         /nblock/ kstart,kstep
 
	if (zon.eq.'s') then
	   nodn=(j-1)*ims*kms+(k-1)*ims+i
	elseif (zon.eq.'a') then
	   nodn=kstart+(k-1)*kstep+(j-1)*ima+i
	elseif (zon.eq.'b') then
	   nodn=kstart+(k-1)*kstep+ima*jma+(j-1)*imb+i
	endif
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine nodnumber(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                 nodb,nb1,nb2,nb3, prog,rzero,elast)
C-----------------------------------------------C
C  Rutinen skapar matriserna NODS,NODA & NODB   C
C-----------------------------------------------C
	implicit none
        integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1           na1,na2,na3,noda(na1,na2,na3),
     2           nb1,nb2,nb3,nodb(nb1,nb2,nb3)
        integer  i,j,k,ilocal,ktmp,rzero,elast
C.......INTEGER FUNCTIONS
        integer  nodn
	logical   left
        character prog*20
c test variables
c	CHARACTER SVAR*1
c
	integer       mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
	common  /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
        integer       ims,jms,kms,is,ima,jma,imb,jmb,kma
        common  /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
        integer          kstart,kstep
        common  /nblock/ kstart,kstep
        integer          ksr1,kar1,kar2
        common  /reduce/ ksr1,kar1,kar2
c
c...... Nollst{llning av alla matriser, samt ber. av konstanter :
	do k=1, ns3
	  do j=1, ns2
	     do i=1, ns1
	       nods(i,j,k)=0
	    enddo
	  enddo
	enddo
	do k=1, na3
	  do j=1, na2
	     do i=1, na1
	       noda(i,j,k)=0
	    enddo
	  enddo
	enddo
	do k=1, nb3
	  do j=1, nb2
	     do i=1, nb1
	       nodb(i,j,k)=0
	    enddo
	  enddo
	enddo
 	kstart=int(0.001*real(ims*jms*kms))*1000+1000
       	kstep=int(0.01*real(ima*jma+imb*jmb))*100+100
	write(*,'(t15,2(a,i6),a/)')'(kstart=',kstart,' kstep=',kstep,')'
 
C ZON S (embrazing the crack front) :
C====================================
	do j=1, jms
	   do k=1, kms
	      do i=1, ims
                 nods(i,j,k)=nodn('s',i,j,k)
	      enddo
	   enddo
	enddo
C . . . Ev. nodreducering i zon S ( om SFRED > 0 ) :
	if (sfred.gt.0) then
	   ilocal=1
	   left=.true.
	   call reduce_mesh_2_to_1(nods,ns1,ns2,ns3,ilocal,3,ims,4,
     &                    1,jms,1,ksr1,left)
	   do k=ksr1+2, kms
	      do i=2, ims, 2
	         do j=1, jms
	            nods(i,j,k)=0
	         enddo
	      enddo
	   enddo
	endif
 
C ZON A :
C========
C       Node numbers in Zone A before the reduction of elements.
	do k=1, kma
	   do j=1, jma
	      do i=1, ima
                 noda(i,j,k)=nodn('a',i,j,k)	
	      enddo
	   enddo
	enddo
C . . .The small limited area under the crack front adjacent to the
C      free surface on the cracked side of the specimen.
	if (rtype.le.1) then
	   ktmp=kma
	else
	   ktmp=kar1+2
	endif
	do j=jma-8, 1, -8
	   do k=1, ktmp
	      noda(1,j+3,k)=0
	      noda(2,j+3,k)=0
	      noda(1,j+1,k)=0
	      noda(4,j+1,k)=0
	      noda(2,j,k)=0
	      noda(4,j,k)=0
	   enddo
	   if (j.gt.1) then
	      do k=1, ktmp
	         noda(1,j-1,k)=0
	         noda(4,j-1,k)=0
	         noda(1,j-3,k)=0
	         noda(2,j-3,k)=0
	      enddo
	   endif
	enddo
C. . . Node numbers in Zone A after the reduction of elements.
	if (rtype.ge.1) then
	   ilocal=1
	   if (mod(ma,4).eq.0) then
              left=.true.
	   else
              left=.false.
	   endif
	   call reduce_mesh_2_to_1(noda,na1,na2,na3,ilocal,7,ima,4,
     &                    1,jma,1,kar1,left)
	   do k=kar1+2, kma
	      do i=6, ima, 2
	         do j=1, jma
	            noda(i,j,k)=0
	         enddo
	      enddo
	   enddo
C . . . . .The small limited area under the crack front
	   if (mod(ma,4).ne.0) then
	      do i=1, 5
	         do j=1, jma
	            noda(i,j,kar1-1)=0
	            noda(i,j,kar1+1)=0
	         enddo
	      enddo
	   endif
	endif
 
	if (rtype.eq.2) then
C . . . . .The small limited area under the crack front
	   if (mod(na,4).ne.0) then
              left=.true.
	   else
              left=.false.
	   endif
	   call reduce_mesh_3_to_1(noda,na1,na2,na3,3,jma,4,kar2,left)
	   do k=kar2+2, kma
	      do j=1, jma, 2
	         noda(2,j,k)=0
	         noda(4,j,k)=0
	      enddo
	   enddo
	   ilocal=2
	   if (mod(na,4).ne.0) then
              left=.true.
	   else
              left=.false.
	   endif
	   call reduce_mesh_2_to_1(noda,na1,na2,na3,ilocal,5,ima,1,
     &                    3,jma,4,kar2,left)
	   do k=kar2+2, kma
	      do j=2, jma, 2
	         do i=1, ima
	            noda(i,j,k)=0
	         enddo
	      enddo
	   enddo
	endif
 
C . .  Delete the node numbers in crack Zone area. !
C . . . ZON S :
	do j=1, jms
	   nods(ims-1,j,1)=0
	   nods(ims-1,j,2)=0
	   nods(ims,j,1)=0
	   nods(ims,j,2)=0
	enddo
C . . . ZON A :
	do k=1, 2*mv
	   do j=1, jma
	      do i=2*(m1-mh)+2, 2*(m1+mh)
	         noda(i,j,k)=0
	      enddo
	   enddo
	enddo
C . . . Delete the node numbers for j > jma-4 & i < 5 ( ZON A )
	do k=1, kma
	   do j=jma-3, jma
	      do i=1,4
	         noda(i,j,k)=0
	      enddo
	   enddo
	enddo
 
C ZON B :
C========
C       Node numbers in Zone B before the reduction of elements.
	do k=1, kma
	   do j=1, jmb
	      do i=1, imb
                 nodb(i,j,k)=nodn('b',i,j,k)	
	      enddo
	   enddo
	enddo
C. . . Node numbers in Zone B after the reduction of elements.
	if (rtype.ge.1) then
	   do i=1, imb
	      do j=1, jmb
	         nodb(i,j,kar1-1)=0
	         nodb(i,j,kar1+1)=0
	      enddo
	   enddo
	endif
	if (rtype.eq.2) then
	   ilocal=2
	   if (mod(nb,4).ne.0) then
	      left=.true.
	   else
	      left=.false.
	   endif
	   call reduce_mesh_2_to_1(nodb,nb1,nb2,nb3,ilocal,1,imb,1,
     &                    3,jmb,4,kar2,left)
	   do k=kar2+2, kma
	      do j=2, jmb, 2
	         do i=1, imb
	            nodb(i,j,k)=0
	         enddo
	      enddo
	   enddo
	endif
 
C.......6) Fixa till de olika zonernas gemensamma nodnummer :
	call common_nodes(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                    nodb,nb1,nb2,nb3, prog,rzero,elast)
C
C--------- For Test of the Program ! -------------------------
C.......7) Eventuell utskrift av nodpositinerna i matriserna :
C	WRITE(*,'(T10,A,$)')
C     &   '* Utskrift av  matriserna NODS,NODA & NODB (J/N) ? :'
C	READ(*,'(A)') SVAR
C	IF ((SVAR.EQ.'J').OR.(SVAR.EQ.'j')) THEN
C	   CALL PRINT_NODMAT(NODS,NS1,NS2,NS3, NODA,NA1,NA2,NA3,
c      &                       NODB,NB1,NB2,NB3)
C       ENDIF
 
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine reduce_mesh_2_to_1(nod,n1,n2,n3,ilocal,
     &           i1,i2,ist,j1,j2,jst,k,left)
C  The routnine removes nodes in reduced mesh areas.
	implicit none
	integer  n1,n2,n3,nod(n1,n2,n3),ilocal,i1,i2,ist,
     &           j1,j2,jst,k,i,j
	logical  left
	if (ilocal.eq.1) then
	   do i=i1, i2, ist
	      if (left) then
	         do j=j1, j2, jst
	            nod(i-1,j,k+1)=0
	            nod(i+1,j,k-1)=0
	            nod(i+2,j,k-1)=0
	            nod(i+2,j,k+1)=0
	         enddo
	         left=.false.
	      else
	         do j=j1, j2, jst
	            nod(i-2,j,k-1)=0
	            nod(i-2,j,k+1)=0
	            nod(i-1,j,k-1)=0
	            nod(i+1,j,k+1)=0
	         enddo
	         left=.true.
	      endif
	   enddo
	else
	   do j=j1, j2, jst
	      if (left) then
	         do i=i1, i2, ist
	            nod(i,j-1,k+1)=0
	            nod(i,j+1,k-1)=0
	            nod(i,j+2,k-1)=0
	            nod(i,j+2,k+1)=0
	         enddo
	         left=.false.
	      else
	         do i=i1, i2, ist
	            nod(i,j-2,k-1)=0
	            nod(i,j-2,k+1)=0
	            nod(i,j-1,k-1)=0
	            nod(i,j+1,k+1)=0
	         enddo
	         left=.true.
	      endif
	   enddo
	endif
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine reduce_mesh_3_to_1(nod,n1,n2,n3,j1,j2,jst,kr,left)
C  The routnine removes nodes in reduced mesh areas. Three elements
C  is reduced to one element.
	implicit none
	integer  n1,n2,n3,nod(n1,n2,n3),j1,j2,jst,kr,k,j,jj
	logical  left
	do j=j1, j2, jst
	   if (left) then
C . . . . . I=1 :
	      do k=kr-2, kr+2
	         nod(1,j-1,k)=0
	         nod(1,j+1,k)=0
	      enddo
	      do jj=j-2, j+2
	         nod(1,jj,kr-1)=0
	         nod(1,jj,kr+1)=0
	      enddo
C . . . . . I=2 :
	      do k=kr-2, kr+2
	         nod(2,j+2,k)=0
	         nod(2,j-1,k)=0
	      enddo
	      do jj=j-2, j+2
	         nod(2,jj,kr-1)=0
	         nod(2,jj,kr+2)=0
	      enddo
C . . . . . I=3 :
	      nod(3,j+2,kr-1)=0
	      nod(3,j+2,kr+1)=0
	      nod(3,j+1,kr-1)=0
	      nod(3,j+1,kr+2)=0
	      nod(3,j-1,kr+1)=0
	      nod(3,j-1,kr+2)=0
C . . . . . I=4 :
	      do k=kr-2, kr+2
	         nod(4,j+2,k)=0
	         nod(4,j+1,k)=0
	      enddo
	      do jj=j-2, j+2
	         nod(4,jj,kr+1)=0
	         nod(4,jj,kr+2)=0
	      enddo
C . . . . . I=5 :
	      nod(5,j+2,kr-1)=0
	      nod(5,j+2,kr+1)=0
	      nod(5,j+1,kr-1)=0
	      nod(5,j+1,kr+2)=0
	      nod(5,j-1,kr+1)=0
	      nod(5,j-1,kr+2)=0
	      left=.false.
	   else
C . . . . . I=1 :
	      do k=kr-2, kr+2
	         nod(1,j-1,k)=0
	         nod(1,j+1,k)=0
	      enddo
	      do jj=j-2, j+2
	         nod(1,jj,kr-1)=0
	         nod(1,jj,kr+1)=0
	      enddo
C . . . . . I=2 :
	      do k=kr-2, kr+2
	         nod(2,j-2,k)=0
	         nod(2,j+1,k)=0
	      enddo
	      do jj=j-2, j+2
	         nod(2,jj,kr-1)=0
	         nod(2,jj,kr+2)=0
	      enddo
C . . . . . I=3 :
	      nod(3,j-2,kr-1)=0
	      nod(3,j-2,kr+1)=0
	      nod(3,j-1,kr-1)=0
	      nod(3,j-1,kr+2)=0
	      nod(3,j+1,kr+1)=0
	      nod(3,j+1,kr+2)=0
C . . . . . I=4 :
	      do k=kr-2, kr+2
	         nod(4,j-2,k)=0
	         nod(4,j-1,k)=0
	      enddo
	      do jj=j-2, j+2
	         nod(4,jj,kr+1)=0
	         nod(4,jj,kr+2)=0
	      enddo
C . . . . . I=5 :
	      nod(5,j-2,kr-1)=0
	      nod(5,j-2,kr+1)=0
	      nod(5,j-1,kr-1)=0
	      nod(5,j-1,kr+2)=0
	      nod(5,j+1,kr+1)=0
	      nod(5,j+1,kr+2)=0
	      left=.true.
	   endif
	enddo
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine common_nodes(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                          nodb,nb1,nb2,nb3, prog,rzero,elast)
c
c-----------------------------------------------------------------C
c   Rutinen ser till att angr{nsande omr}den till Zon S ( Zon A)  C
c   resp Zon A (Zon B och ZonC) f}r Zon S resp Zon A :s nodnr     C
c   i de gemensamma noderna                                       C
c-----------------------------------------------------------------C
	implicit none
c
        integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1           na1,na2,na3,noda(na1,na2,na3),
     2           nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
	integer  i,j,t,ia1,ia2,k,t1,ts,ts1,ts2,tsd,rzero,elast
c
        character prog*20
c
        integer      mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
	common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
        integer      ims,jms,kms,is,ima,jma,imb,jmb,kma
        common /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
c.......1) Common nodes within  Zone S :
 	do j=1, jms
	   nods(ims-1,j,3)=nods(ims-2,j,2)
	   nods(ims,j,3)=nods(ims-2,j,1)
	enddo
        if ( ((prog.eq.'WARP3D').or.(prog.eq.'PATRAN')).and.
     &       (rzero.eq.1).and.(elast.eq.1) ) then
           do j=1, jms
              do i=2, (ims+1)/2
                 nods(i,j,1)=nods(1,j,1)
              enddo
           enddo
        endif
c
c.......2) Common nodes in Zon S and Zon A  :
	ia1=2*(m1+mh)+1
	ia2=2*(m1-mh)+1
	t1=2*mv+1
	if (sfred.eq.0) then
	   tsd=1
	   ts1=t1
	   ts2=2*(mf-mv)+1
	else
	   tsd=2
	   ts1=4*mv+1
	   ts2=2*(mf-2*mv)+1
	endif
	do j=1, jma
	   do t=1, t1
	      ts=tsd*(t-1)+1
	      noda(ia1,j,t)=nods(ts,j,kms)
	      noda(ia2,j,t)=nods(ims-(ts-1),j,kms)
	   enddo
	   do t=1,2*mh
	      ts=tsd*t
	      noda((ia1-t),j,t1)=nods(ts1+ts,j,kms)
	      noda((ia2+t),j,t1)=nods(ts2-ts,j,kms)
	   enddo
	enddo
 
C.......3) Common nodes within  Zone A  :
	if ((m1-mh).eq.2) then
	   do k=1, 2*mv+1
	      do t=4, 1, -1
	         noda(t,jma-4,k)=noda(5,jma-(t-1),k)
	      enddo
	   enddo
	   do k=2*mv+2, kma
	      do t=4, 1, -1
	         noda(5,jma-(t-1),k)=noda(t,jma-4,k)
	      enddo
	   enddo
	else
	   do k=1, kma
	      do t=4, 1, -1
	         noda(5,jma-(t-1),k)=noda(t,jma-4,k)
	      enddo
	   enddo
	endif
	
C.......4) Common nodes in Zon A and Zon B  :
	t1=2*(na-nb)
 
	do k=1, kma
	   do t=1, jmb
	      nodb(1,t,k)=noda(ima,t+t1,k)
	   enddo
	enddo
 
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
C----------------------------------------------------------------C
C        NEDANSTAENDE SUBROUTIN/ER ANVANDS ENBART VID            C
C               TEST AV SJALVA PROGRAMMET  !!!!!                 C
C----------------------------------------------------------------C
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine print_nodmat(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                          nodb,nb1,nb2,nb3)
C----------------------------------------------------------------C
C  Rutinen skriver ut matriserna NODA - NODS p} filen NODMAT.DAT C
C----------------------------------------------------------------C
      implicit none
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     &         sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma,i,j,k,
     &         vek(100),   ksr1,kar1,kar2
c
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     &       /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma,
     &       /reduce/ ksr1,kar1,kar2
c
	   open(unit=11,file='nodmat.dat',status='unknown')
	   open(unit=12,file='nodmatpr.dat',status='unknown')
c
	   write(11,'(//t10,a)') '  *** modell information *** '
	   write(11,'(t10,9(a,i2,tr1))')' mr=',mr,' mf=',mf,
     &     ' sfred=',sfred,' mv=',mv,
     &     ' mh=',mh,' m1=',m1,' m2=',m2,' ma=',ma,' na=',na
           write(11,'(t10,7(a,i2,tr1))') ' mb=',mb,' nb=',nb,' lt=',lt,
     &     ' lred=',lred,' is=',is,' ima= ',ima,' jma=',jma
           write(11,'(t10,6(a,i2,tr1))')' imb=',imb,' jmb= ',jmb,
     &     ' ims=',ims,'  jms=',jms,' kms=',kms,' kma=',kma
           write(11,'(t10,4(a,i2,tr1))')' rtype=',rtype,' ksr1= ',ksr1,
     &     ' kar1=',kar1,'  kar2=',kar2
 
	   write(12,'(//t10,a)') '  *** modell information *** '
	   write(12,'(t10,9(a,i2,tr1))')' mr=',mr,' mf=',mf,
     &     ' sfred=',sfred,' mv=',mv,
     &     ' mh=',mh,' m1=',m1,' m2=',m2,' ma=',ma,' na=',na
           write(12,'(t10,7(a,i2,tr1))') ' mb=',mb,' nb=',nb,' lt=',lt,
     &     ' lred=',lred,' is=',is,' ima= ',ima,' jma=',jma
           write(12,'(t10,6(a,i2,tr1))')' imb=',imb,' jmb= ',jmb,
     &     ' ims=',ims,'  jms=',jms,' kms=',kms,' kma=',kma
           write(12,'(t10,4(a,i2,tr1))')' rtype=',rtype,' ksr1= ',ksr1,
     &     ' kar1=',kar1,'  kar2=',kar2
 
	   write(11,'(//t10,a)') '=> matris nods(i,j,k) :'
	   write(11,'(t10,a/)') '   ===================='
	   write(12,'(//t10,a)') '=> matris nods(i,j,k) :'
	   write(12,'(t10,a/)') '   ===================='
	   do j=1, jms
C	      WRITE(11,'(/T2,A,I2,A)') '* J=',J,'  ====== I ========='
	      write(12,'(/t2,a,i2,a)') '* j=',j,'  ====== i ========='
	      do k=1,kms
C	         WRITE(11,'(T5,A,I2,TR3,9I6,3(/T14,9I6))')
C     &           'K=',K, (NODS(I,J,K), I=1, IMS)
	         write(12,'(t3,a,i2,tr2,17i6,3(/t15,17i6))')
     &           'k=',k, (nods(i,j,k), i=1, ims)
	      enddo
	   enddo
	   do j=1, jms
	      write(11,'(/t4,a,i2,a)') '* j=',j,'  ====== i ========='
	      do k=1, kms
	        do i=1, ims
	           if (nods(i,j,k).ne.0) then
	              vek(i)=1
	           else
	              vek(i)=0
	           endif
	        enddo
	        write(11,'(t10,a,i2,tr5,100i1)')'k=',k,(vek(i),i=1,ims)
	      enddo
	   enddo
 
	   write(11,'(//t10,a)') '=> matris noda(i,j,k) :'
	   write(11,'(t10,a/)') '   ===================='
	   write(12,'(//t10,a)') '=> matris noda(i,j,k) :'
	   write(12,'(t10,a/)') '   ===================='
	   do k=1, kma
C	      WRITE(11,'(/T2,A,I2,A)') '* K=',K,'  ====== I ========='
	      write(12,'(/t2,a,i2,a)') '* k=',k,'  ====== i ========='
	      do j=1,jma
C	         WRITE(11,'(T5,A,I2,TR3,10I6/T14,9I6/T14,9I6)')
C     &           'J=',J, (NODA(I,J,K), I=1, IMA)
	         write(12,'(t5,a,i2,tr3,19i6/t14,9i6)')
     &           'j=',j, (noda(i,j,k), i=1, ima)
 	      enddo
	   enddo
	   do k=1, kma
	      write(11,'(/t4,a,i2,a)') '* k=',k,'   ====== i  ======='
	      do j=1, jma
	       do i=1, ima
	          if (noda(i,j,k).ne.0) then
	            vek(i)=1
	          else
	            vek(i)=0
		  endif
	       enddo
	       write(11,'(t10,a,i2,tr5,100i1)')'j=',j,(vek(i),i=1,ima)
	      enddo
	   enddo
 
	   write(11,'(//t10,a)') '=> matris nodb(i,j,k) :'
	   write(11,'(t10,a/)') '   ===================='
	   write(12,'(//t10,a)') '=> matris nodb(i,j,k) :'
	   write(12,'(t10,a/)') '   ===================='
	   do k=1, kma
C	      WRITE(11,'(/T2,A,I2,A)') '* K=',K,'  ====== I ========='
	      write(12,'(/t2,a,i2,a)') '* k=',k,'  ====== i ========='
	      do j=1,jmb
C	         WRITE(11,'(T5,A,I2,TR3,10I6/T14,9I6/T14,9I6)')
C     &           'J=',J, (NODB(I,J,K), I=1, IMB)
	         write(12,'(t5,a,i2,tr3,19i6/t14,9i6)')
     &           'j=',j, (nodb(i,j,k), i=1, imb)
 	      enddo
	   enddo
	   do k=1, kma
	      write(11,'(/t4,a,i2,a)') '* k=',k,'  ====== i  ======='
	      do j=1, jmb
	       do i=1, imb
	          if (nodb(i,j,k).ne.0) then
	            vek(i)=1
	          else
	            vek(i)=0
		  endif
	       enddo
	       write(11,'(t10,a,i2,tr5,100i1)')'j=',j,(vek(i),i=1,imb)
	      enddo
	   enddo
 
	   close(11)
	   close(12)
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
