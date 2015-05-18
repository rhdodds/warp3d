c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine sort_out_nodes( elnum, no_of_nodes,
     1     etyp,prog, el_ch_zs, slice, eln_d,nstk,dsl,grp,
     2     nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3,
     3     nea2,neb2, estk_s,estk_s0,estk_a,estk_b,estk_b0)
c
c--- The routine sort out the uniqe nodes in the models
c
      implicit none
c
      include 'mesh3d_scp_common_eln.f'
      integer  eln_d(0:iem),e(iem,27)
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  nea2,neb2, estk_s(2,1),estk_s0(2,1), estk_a(2,nea2,1),
     1         estk_b(2,neb2,1),estk_b0(2,1),
     2         es_a(25,25,2),es_b(25,25,2)
c
      integer  elnum,no_of_nodes, etyp, el_ch_zs,slice,
     1         i,ii,i1,i2,j,jj,j1,k,nn,enum,el_nodes,d1,d2,d3,d4,
     2         nstk(0:500,7),n,in1,in2,kst,dsl(200),grp
CCOMMON INTEGER
      integer  ims,jms,kms,is,ima,jma,imb,jmb,kma,
     &         mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
 
C  INTEGER FUNTIONS
      integer  numb_of_el_nodes
      character prog*20
      common /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
      common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred
 
C.... If elements with variable nodes will be used
      eln_d(0)=0
      do i=1, elnum
         eln_d(i)=1
      enddo
      grp=1
      dsl(1)=elnum
C... Change local node numbering order in elements if
C    ABAQUS-analysis or WARP3D-analysis
       if (prog.eq.'ABAQUS') then
	   do i=1, elnum
	      call change_elnum_to_abaqus(i)
	   enddo
        elseif ( (prog.eq.'WARP3D').or.(prog.eq.'PATRAN') ) then
	   do i=1, elnum
	      call change_elnum_to_abaqus(i)
	   enddo
	   do i=1, elnum
	      call change_elnum_to_warp3d(i)
	   enddo
        endif
C.... Delete unnessasary nodes in ELN
	do i=1, elnum
	   do j=etyp+1, 27
	      eln(i,j)=0
	   enddo
	enddo
C.. If  EL_CH_ZS > 0 (inner element rings will full integration)
	if (el_ch_zs.gt.0) then
	   do i=1, elnum
	      eln_d(i)=2
	   enddo
	   do j=1, (jma-1)/2
	      i1=estk_s(1,j)
	      i2=estk_s(1,j) + (ims-1)/2-2 + (el_ch_zs-1)*(ims-1)/2
	      estk_s(1,j)=i2+1
	      estk_s0(1,j)=i1
	      estk_s0(2,j)=i2
	      do i=i1, i2
	         eln_d(i)=1
	      enddo
	   enddo
	endif
C.. Change the element numbering order if ADINA-analysis & EL_CH_ZS > 0
        if ((el_ch_zs.gt.0).and.(prog.eq.'ADINA')) then
	   d1=0
	   d2=0
	   do i=1, elnum
	      if (eln_d(i).eq.1) then
	         d1=d1+1
	         eln(i,0)=d1
	      else
	         d2=d2+1
	         eln(i,0)=d2
	      endif
	   enddo
	   dsl(1)=d1
	   dsl(2)=d2
	   grp=2
	endif
 
C.. Change the element numbering order if ADINA-analysis & SLICE > 0
        if ((slice.gt.0).and.(prog.eq.'ADINA')) then
	   if (el_ch_zs.gt.0) then
	      dsl(1)=d1
	      grp=1
	   else
	      grp=0
	   endif
	   do i=grp+1, 200
	      dsl(i)=0
	   enddo
C... Slices around the crack front in Zon S and Zon A :
	   do j=1, na
	      grp=grp+1
	      i1=estk_s(1,j)
	      i2=estk_s(2,j)
	      do i=i1, i2
	         eln_d(i)=grp
	         dsl(grp)=dsl(grp)+1
	         eln(i,0)=dsl(grp)
	      enddo
	      do k=1, slice
	         i1=estk_a(1,k,j)
	         i2=estk_a(2,k,j)
	         do i=i1, i2
	            eln_d(i)=grp
	            dsl(grp)=dsl(grp)+1
	            eln(i,0)=dsl(grp)
	         enddo
	      enddo
	   enddo
C... First part in Zon B :
	   grp=grp+1
	   do k=1, slice
	      do i=estk_b(1,k,1), estk_b(2,k,mb)
	         eln_d(i)=grp
	         dsl(grp)=dsl(grp)+1
	         eln(i,0)=dsl(grp)
	      enddo
	   enddo
C... The remaining elements in Zon A and Zon B:
	   grp=grp+1
	   do i=estk_a(1,slice+1,1),  estk_b(2,lt,mb)
	      eln_d(i)=grp
	      dsl(grp)=dsl(grp)+1
	      eln(i,0)=dsl(grp)
	   enddo
	endif
C.... Correct the element statistics
C	IF ( (PROG.EQ.'ADINA').AND.
C     &       ( (SLICE.GT.0).OR.(EL_CH_ZS.GT.0) ) ) THEN
C	   DO I=1, 25
C	    DO K=1, 2
C	       IF (ESTK_S(K,I).GT.0)  ESTK_S(K,I)=ELN(ESTK_S(K,I),0)
C	       IF (ESTK_S0(K,I).GT.0) ESTK_S0(K,I)=ELN(ESTK_S0(K,I),0)
C	    ENDDO
C	    DO J=1, 25
C	     DO K=1, 2
C	      IF (ESTK_A(K,I,J).GT.0) ESTK_A(K,I,J)=ELN(ESTK_A(K,I,J),0)
C	      IF (ESTK_B(K,I,J).GT.0) ESTK_B(K,I,J)=ELN(ESTK_B(K,I,J),0)
C	     ENDDO
C	    ENDDO
C	   ENDDO
C	ENDIF
CCCCCCCC
C.... Determine the unique nodes
	do i=0, inm
	   nnr(i)=0
	enddo
	do i=1, elnum
	   el_nodes = numb_of_el_nodes(i)
	   do j=1, el_nodes
	      nnr(eln(i,j))=1
	   enddo
	enddo
C... Take away jumps in node-numbers if ADINA analysis
        if ((prog.eq.'ADINA').or.(prog.eq.'WARP3D').or.
     &      (prog.eq.'PATRAN').or.(prog.eq.'ABAQUS')) then
	   nn=0
C . . . . Change the node numbers
	   do i=1, inm
	      if (nnr(i).gt.0) then
	         nn=nn+1
	         nnr(i)=nn
	      endif
	   enddo
C . . . . Change the node numbers in the elements
	   do i=1, elnum
	      do j=1, 27
	         if (eln(i,j).gt.0) then
	            eln(i,j)=nnr(eln(i,j))
	         endif
	      enddo
	   enddo
	else
	   nn=0
	   do i=1, inm
	      if (nnr(i).gt.0) then
	         nn=nn+1
	         nnr(i)=i
	      endif
	   enddo
	endif
C... Create node statistics
C  ZON S
	nstk(0,3)=nnr(nods(1,1,1))
	nstk(0,4)=nnr(nods(ims,jms,kms))
	n=0
	do i=nods(1,1,1), nods(ims,1,kms)
	   if (nnr(i).gt.0) n=n+1
	enddo
	nstk(0,1)=n
	n=0
	do i=nods(ims,1,kms)+1, nods(1,3,1)-1
	   if (nnr(i).gt.0) n=n+1
	enddo
	nstk(0,2)=n
	n=0
	do i=nods(1,1,1), nods(ims,jms,kms)
	   if (nnr(i).gt.0) n=n+1
	enddo
	nstk(0,7)=n
C  ZON A & ZON B
	kst=nnr(noda(1,1,1))-1
	do k=1, kma
C . . . zon A:
	   do j=1, jma
	      do i=1, ima
	         if (nnr(noda(i,j,k)).gt.kst) then
	            nstk(k,1)=nnr(noda(i,j,k))
	            in1=noda(i,j,k)
	            goto 10
	         endif
	      enddo
	   enddo
	   nstk(k,1)=0
10	   continue
	   if (nstk(k,1).gt.kst) then
	      do j=jma, 1, -1
	         do i=ima, 1, -1
	            if (nnr(noda(i,j,k)).gt.kst) then
	               nstk(k,2)=nnr(noda(i,j,k))
	               in2=noda(i,j,k)
	               goto 15
	            endif
	         enddo
	      enddo
	   endif
	   nstk(k,2)=0
15	   continue
	   if (nstk(k,1).gt.kst) then
	      n=0
	      do i=in1, in2
	         if (nnr(i).gt.kst) n=n+1
	      enddo
	      nstk(k,3)=n
	   else
	      nstk(k,3)=0
	   endif
C . . . zon B:
	   do j=1, jmb
	      do i=2, imb
	         if (nnr(nodb(i,j,k)).gt.kst) then
	            nstk(k,4)=nnr(nodb(i,j,k))
	            in1=nodb(i,j,k)
	            goto 20
	         endif
	      enddo
	   enddo
	   nstk(k,4)=0
20	   continue
	   if (nstk(k,4).gt.kst) then
	      do j=jmb, 1, -1
	         do i=imb, 2, -1
	            if (nnr(nodb(i,j,k)).gt.kst) then
	               nstk(k,5)=nnr(nodb(i,j,k))
	               in2=nodb(i,j,k)
	               goto 25
	            endif
	         enddo
	      enddo
	   endif
	   nstk(k,5)=0
25	   continue
	   if (nstk(k,4).gt.kst) then
	      n=0
	      do i=in1, in2
	         if (nnr(i).gt.kst) n=n+1
	      enddo
	      nstk(k,6)=n
	   else
	      nstk(k,6)=0
	   endif
	   nstk(k,7)=nstk(k,3)+nstk(k,6)
	enddo
	n=0
	do k=0, kma
	   n=n+nstk(k,7)
	enddo
	write(*,'(t15,a,i5,a,i5,a)') '=> the model contains ',
     &               nn,' unique nodes (nstk=',n,')'
c
        no_of_nodes = nn
c
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine change_elnum_to_abaqus(i)
c
c-- The routine change the local element numbering from ADINA to ABAQUS
c
	implicit none
c
        include 'mesh3d_scp_common_eln.f'
c
	integer  i,j,n(27)
	do j=1, 27
	   n(j)=eln(i,j)
	enddo
	do j=1, 4
	   eln(i,j)=n(j+4)
	   eln(i,j+8)=n(j+12)
	   eln(i,j+4)=n(j)
	   eln(i,j+12)=n(j+8)
	   eln(i,j+23)=n(j+21)
	enddo
	eln(i,22)=n(27)
	eln(i,23)=n(26)
C . . . The nodes 17, 18, 19, 20 and 21 will be unchanged
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine change_elnum_to_warp3d(i)
c
c--The routine change the local element numbering from ABAQUS to WARP3D
c
	implicit none
c
        include 'mesh3d_scp_common_eln.f'
c
	integer i,j,n(27)
	do j=1, 20
	   n(j)=eln(i,j)
	enddo
	do j=1, 12
	   eln(i,j)=n(j)
	enddo
        eln(i,13) = n(17)
        eln(i,14) = n(18)
        eln(i,15) = n(19)
        eln(i,16) = n(20)
        eln(i,17) = n(13)
        eln(i,18) = n(14)
        eln(i,19) = n(15)
        eln(i,20) = n(16)
C
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine write_out_nod_el(job,jobh,elnum,etyp,eln_d,
     &                            el_ch_zs,slice,prog)
c
c-- The routine writes out nod coordinates and element definitions
c-- on file.
c
      implicit none
c
      include 'mesh3d_scp_common_eln.f'
      integer  eln_d(0:iem)
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  elnum,etyp, el_ch_zs, slice,
     &         jobh,i,j,el_nodes,uni
c
c  INTEGER FUNTIONS
      integer   numb_of_el_nodes
      character job*40,fil1*40,fil2*40,fil3*40,fil4*40,fil5*40,prog*20
c
      write(*,'(/t5,a,a/)') '* writing node coord. and element ',
     &                      'def. on output files . . . .'
c
c.. Write out the nod coordinates
      if (prog.eq.'ABAQUS') then
         fil1=job(1:jobh)//'.015'
      elseif (prog.eq.'ADINA') then
         fil1=job(1:jobh)//'.15'
      elseif ((prog.eq.'WARP3D').or.(prog.eq.'PATRAN')) then
         fil1=job(1:jobh)//'.crd'
      endif
      open(unit=15,file=fil1,status='unknown')
      if ((prog.eq.'ABAQUS').or.(prog.eq.'ADINA')) then
         do i=1, inm
            if (nnr(i).gt.0) write(15,101) nnr(i),',',
     &                       npos(i,1),',',npos(i,2),',',npos(i,3)
         enddo
c      elseif ((prog.eq.'WARP3D').or.(prog.eq.'PATRAN')) then
c         do i=1, inm
c            if (nnr(i).gt.0) write(15,201) nnr(i),
c     &                       npos(i,1),npos(i,2),npos(i,3)
c         enddo
      endif
      close(15)
c
c.. Write out the element definitions
c..ABAQUS:
      if (prog.eq.'ABAQUS') then
         fil2=job(1:jobh)//'.016'
         open(unit=16,file=fil2,status='unknown')
         if (el_ch_zs.gt.0) then
	      fil3=job(1:jobh)//'.017'
	      open(unit=17,file=fil3,status='unknown')
	      do i=1, elnum
		 el_nodes=numb_of_el_nodes(i)
	         if (eln_d(i).eq.1) then
	            write(16,104) (eln(i,j), j=0, el_nodes)
	         else
	            write(17,104) (eln(i,j), j=0, el_nodes)
	         endif
	      enddo
	      close(17)
         elseif (etyp.eq.8) then
	      do i=1, elnum
	         write(16,102) (eln(i,j), j=0, etyp)
	      enddo
         elseif (etyp.eq.20) then
	      do i=1, elnum
	         write(16,103) (eln(i,j), j=0, etyp)
	      enddo
         else
	      do i=1, elnum
	         write(16,104) (eln(i,j), j=0, etyp)
	      enddo
         endif
         close(16)
      elseif (prog.eq.'ADINA') then
C..ADINA:
         if ((el_ch_zs.gt.0).or.(slice.gt.0)) then
	    do i=1, elnum
	       uni=50+eln_d(i)
	       el_nodes=numb_of_el_nodes(i)
	       if (el_nodes.eq.27) then
	          write(uni,114) (eln(i,j), j=0, 27)
	       elseif (el_nodes.eq.20) then
	          write(uni,113) (eln(i,j), j=0, 20)
	       else
	          write(uni,112) (eln(i,j), j=0, 8)
	       endif
	    enddo
         elseif (etyp.eq.8) then
	    fil2=job(1:jobh)//'.21'
	    open(unit=16,file=fil2,status='unknown')
	    do i=1, elnum
	       write(16,112) (eln(i,j), j=0, etyp)
	    enddo
	    close(16)
         elseif (etyp.eq.20) then
	    fil2=job(1:jobh)//'.21'
	    open(unit=16,file=fil2,status='unknown')
	    do i=1, elnum
	       write(16,113) (eln(i,j), j=0, etyp)
	    enddo
	    close(16)
         else
	    fil2=job(1:jobh)//'.21'
	    open(unit=16,file=fil2,status='unknown')
	    do i=1, elnum
	       write(16,114) (eln(i,j), j=0, etyp)
	    enddo
	    close(16)
         endif
c      elseif ((prog.eq.'WARP3D').or.(prog.eq.'PATRAN')) then
C..WARP3D & PATRAN:
c         fil2=job(1:jobh)//'.elm'
c         open(unit=16,file=fil2,status='unknown')
c         if (etyp.eq.8) then
c            do i=1, elnum
c               write(16,203) (eln(i,j), j=0, etyp)
c            enddo
c         elseif (etyp.eq.20) then
c            do i=1, elnum
c               write(16,204) (eln(i,j), j=0, etyp)
c            enddo
c         else
c            do i=1, elnum
c               write(16,205) (eln(i,j), j=0, etyp)
c            enddo
c         endif
c         close(16)
      endif
C
 101  format(t3,i5,a,2(g16.9,a),g16.9)
 102  format( t1,9(i5) )
 103  format( t1,16(i5)/t1,5(i5) )
 104  format( t1,16(i5)/t1,12(i5) )
C
 201  format(t3,i6, 3(g20.10) )
 203  format( t1,9(i7) )
 204  format( t1,9(i7),','/t8,8(i7),','/t8,4(i7))
 205  format( t1,9(i7),','/t8,8(i7),','/t8,8(i7),','/t8,8(i7))
c
 112  format( t1,i5,',',4(i5,','),tr2,3(i5,','),i5 )
 113  format( t1,i5,',',4(i5,','),4(i5,','),4(i5,',')/
     &                 t7,4(i5,','),3(i5,','),i5  )
 114  format( t1,i5,',',4(i5,','),4(i5,','),4(i5,',')/
     &                 t7,4(i5,','),4(i5,',')/t7,6(i5,','),i5  )
 115  format( t1,i5,',',4(i5,','),4(i5,','),4(i5,',')/
     &                 t7,4(i5,','),4(i5,','),i5,',',i5  )
      return
      end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	integer function numb_of_el_nodes(i)
	implicit none
c
        include 'mesh3d_scp_common_eln.f'
c
	integer  i,j
	j=27
10      if (eln(i,j).eq.0) then
	   j=j-1
	   goto 10
	endif
	numb_of_el_nodes = j
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine write_el_stat(uni, el_ch_zs,eln_d,elnum,
     &       nea2,neb2,estk_s,estk_s0,estk_a,estk_b,estk_b0)
C--- The routine writes out element statistics
      implicit none
c
      include 'mesh3d_scp_common_eln.f'
      integer  eln_d(0:iem)
c
      integer  nea2,neb2, estk_s(2,1),estk_s0(2,1), estk_a(2,nea2,1),
     &         estk_b(2,neb2,1),estk_b0(2,1)
c
      integer  el_ch_zs, uni,elnum,sum,e1,e2,ed,
     &         ims,jms,kms,is,ima,jma,imb,jmb,kma, i,j,k
C COMMON INTEGER
 
        common   /maxa/ ims,jms,kms,is,ima,jma,imb,jmb,kma
c
C	FIL=JOB(1:JOBH)//'_el_sta.dat'
C	OPEN(UNIT=29,FILE=FIL,STATUS='UNKNOWN')
C	UNI=29
	write(uni,'(/t5,a/)') ' *** element-statisics ***'
	write(uni,'(t5,a,i2/)') 'el_ch_zs  = ',el_ch_zs
C... ZON S
	write(uni,'(t5,a)') 'zon s :'
	write(uni,'(t5,a)') '======='
	write(uni,'(t3,a,t12,a,t22,a,t30,a,t38,a)')
     &             'j','i1','i2','number','el. group'
	if (el_ch_zs.gt.0) then
	   do j=1, jms-2, 2
	      e1=estk_s0(1,(j+1)/2)
	      e2=estk_s0(2,(j+1)/2)
	      ed=eln_d(e1)
	      if (e1.gt.0) then
	         e1=eln(e1,0)
	      else
	         e1=0
	      endif
	      if (e2.gt.0) then
	         e2=eln(e2,0)
	      else
	         e2=0
	      endif
	      sum=e2-e1+1
	      if (sum.eq.1) sum=0
              write(uni,'(t3,i2,t10,i5,t20,i5,t30,i5,t40,i2)')
     &                j, e1, e2, sum, ed
	   enddo
	endif
	write(uni,'(t2,a)') '- - - - - - - - - - - - - -'
	do j=1, jms-2, 2
	   e1=estk_s(1,(j+1)/2)
	   e2=estk_s(2,(j+1)/2)
	   ed=eln_d(e1)
	   if (e1.gt.0) then
	      e1=eln(e1,0)
	   else
	      e1=0
	   endif
	   if (e2.gt.0) then
	      e2=eln(e2,0)
	   else
	      e2=0
	   endif
	   sum=e2-e1+1
	   if (sum.eq.1) sum=0
	   write(uni,'(t3,i2,t10,i5,t20,i5,t30,i5,t40,i2)')
     &              j, e1, e2, sum, ed
	enddo
C... ZON A
	write(uni,'(/t5,a)') 'zon a :'
	write(uni,'(t5,a)') '======='
	do k=1, kma-2, 2
	   write(uni,'(t10,a,i2)') 'k = ',k
	   write(uni,'(t3,a,t12,a,t22,a,t30,a,t38,a)')
     &             'j','i1','i2','number','el. group'
	   do j=1, jma-2, 2
	      e1=estk_a(1,(k+1)/2,(j+1)/2)
	      e2=estk_a(2,(k+1)/2,(j+1)/2)
	      ed=eln_d(e1)
	      if (e1.gt.0) then
	         e1=eln(e1,0)
	      else
	         e1=0
	      endif
	      if (e2.gt.0) then
	         e2=eln(e2,0)
	      else
	         e2=0
	      endif
	      sum=e2-e1+1
	      if (sum.eq.1) sum=0
	      write(uni,'(t3,i2,t10,i5,t20,i5,t30,i5,t40,i2)')
     &               j, e1, e2, sum, ed
	   enddo
	enddo
C... ZON B
	write(uni,'(/t5,a)') 'zon b :'
	write(uni,'(t5,a)') '======='
	write(uni,'(t2,a)') '- - - - - - - - - - - - - -'
	do k=1, kma-2, 2
	   write(uni,'(t10,a,i2)') 'k = ',k
	   write(uni,'(t3,a,t12,a,t22,a,t30,a,t38,a)')
     &                 'ib','j1','j2','number','el. group'
	   do i=1, imb-2, 2
	      e1=estk_b(1, (k+1)/2, (i+1)/2)
	      e2=estk_b(2, (k+1)/2, (i+1)/2)
	      ed=eln_d(e1)
	      if (e1.gt.0) then
	         e1=eln(e1,0)
	      else
	         e1=0
	      endif
	      if (e2.gt.0) then
	         e2=eln(e2,0)
	      else
	         e2=0
	      endif
	      sum=e2-e1+1
	      if (sum.eq.1) sum=0
	      write(uni,'(t3,i2,t10,i5,t20,i5,t30,i5,t40,i2)')
     &                 i, e1, e2, sum, ed
	   enddo
	enddo
C	CLOSE(UNI)
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c

