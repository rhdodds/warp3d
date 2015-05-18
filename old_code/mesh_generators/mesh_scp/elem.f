c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
        subroutine make_element(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &       nodb,nb1,nb2,nb3, elnum, nea2,neb2,estk_s,estk_a,estk_b)
C-- Rutinen skapar elementen i alla elementlager.
C-- Elementen, 27-nodiga solidelement tilldelas
C-- noder i ADINA foreskriven ordning.
c
      implicit none
c
      include 'mesh3d_scp_common_eln.f'
c
      integer ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1        na1,na2,na3,noda(na1,na2,na3),
     2        nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer nea2,neb2,estk_s(2,1),estk_a(2,nea2,1),estk_b(2,neb2,1)
c
      integer mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,
     1        sfred,ims,jms,kms,is,ima,jma,imb,jmb,kma,
     2        elnum,ii,it,it0,jt,jt0,kend,dk,i,j,k,
     3        ksr1,kar1,kar2
 
      common /nod/  mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     1       /max/  ims,jms,kms,is,ima,jma,imb,jmb,kma,
     2       /reduce/ ksr1,kar1,kar2
 
      write(*,'(t5,a)') '* generatating the elements . . . .'
c
      elnum=0
c
c... ZON S :
c===========
	if (sfred.eq.0) then
 	   do j=1, jms-2, 2
	      k=1
	      estk_s(1,(j+1)/2)=elnum+1
	      do i=1, ims-4, 2
	         elnum=elnum+1
	         call element_ik(elnum, nods,ns1,ns2,ns3,
     1                i,i+1,i+2, i,i+1,i+2,  i,i+1,i+2, j,j+1,j+2,
     2                k,k,k, k+1,k+1,k+1, k+2,k+2,k+2)
	      enddo
	      do k=3, kms-2, 2
	         do i=1, ims-2, 2
	            elnum=elnum+1
	            call element_ik(elnum, nods,ns1,ns2,ns3,
     1                   i,i+1,i+2, i,i+1,i+2, i,i+1,i+2, j,j+1,j+2,
     2                   k,k,k, k+1,k+1,k+1, k+2,k+2,k+2)
	         enddo
	      enddo
	      estk_s(2,(j+1)/2)=elnum
	   enddo
	else
 	   do j=1, jms-2, 2
	      k=1
	      estk_s(1,(j+1)/2)=elnum+1
	      do i=1, ims-4, 2
	         elnum=elnum+1
	         call element_ik(elnum,  nods,ns1,ns2,ns3,
     1                i,i+1,i+2, i,i+1,i+2,  i,i+1,i+2, j,j+1,j+2,
     2                k,k,k, k+1,k+1,k+1, k+2,k+2,k+2)
	      enddo
	      do k=3, ksr1-4, 2
	         do i=1, ims-2, 2
	            elnum=elnum+1
	            call element_ik(elnum,  nods,ns1,ns2,ns3,
     1                   i,i+1,i+2, i,i+1,i+2, i,i+1,i+2, j,j+1,j+2,
     2                   k,k,k, k+1,k+1,k+1, k+2,k+2,k+2)
	         enddo
	      enddo
C . . . . . Reduction of elements
	      k=ksr1-2
	      do ii=5, ims, 8
	         i=ii-4
	         elnum=elnum+1
	         call element_ik(elnum, nods,ns1,ns2,ns3,
     1                i,i+1,i+2, i,i+1,i+2, i,i+1,i+2, j,j+1,j+2,
     2                k,k,k, k+1,k+1,k+1, k+2,k+2,k+2)
	         i=ii-2
	         elnum=elnum+1
	         call element_ik(elnum, nods,ns1,ns2,ns3,
     1                i,i+1,i+2, i,i+1,i+2,  i,i+1,i+2,  j,j+1,j+2,
     2                k,k,k, k+1,k+2,k+2, k+2,k+3,k+4)
	         i=ii
	         elnum=elnum+1
	         call element_ik(elnum, nods,ns1,ns2,ns3,
     1                i,i+1,i+2, i,i+1,i+2, i,i+1,i+2,  j,j+1,j+2,
     2                k,k,k, k+2,k+2,k+1, k+4,k+3,k+2)
	         i=ii+2
	         elnum=elnum+1
	         call element_ik(elnum, nods,ns1,ns2,ns3,
     1                i,i+1,i+2, i,i+1,i+2, i,i+1,i+2, j,j+1,j+2,
     2                k,k,k, k+1,k+1,k+1, k+2,k+2,k+2)
	      enddo
	      k=ksr1
	      do ii=5, ims, 8
	         i=ii-4
	         elnum=elnum+1
	         call element_ik(elnum, nods,ns1,ns2,ns3,
     1                i,i+1,i+2, i,i+2,i+3,  i,i+2,i+4, j,j+1,j+2,
     2                k,k,k, k+1,k+1,k+1, k+2,k+2,k+2)
	         i=ii
	         elnum=elnum+1
	         call element_ik(elnum, nods,ns1,ns2,ns3,
     1                i+2,i+3,i+4, i+1,i+2,i+4, i,i+2,i+4, j,j+1,j+2,
     2                k,k,k, k+1,k+1,k+1, k+2,k+2,k+2)
	      enddo
	      do k=ksr1+2, kms-2, 2
	         do i=1, ims-4, 4
	            elnum=elnum+1
	            call element_ik(elnum, nods,ns1,ns2,ns3,
     1                   i,i+2,i+4, i,i+2,i+4,  i,i+2,i+4, j,j+1,j+2,
     2                   k,k,k, k+1,k+1,k+1, k+2,k+2,k+2)
	         enddo
	      enddo
	      estk_s(2,(j+1)/2)=elnum
	   enddo
	endif
 
C... ZON A & ZON B :
C===================
	if (mod(na,4).ne.0) then
	   jt0=0
	else
	   jt0=2
	endif
	do k=1, 2*mv-1, 2
C . . ZON A
	   jt=jt0
	   do j=1, jma-2, 2
	      estk_a(1, (k+1)/2, (j+1)/2)=elnum+1
	      if (j.le.jma-6) call element_au(elnum, noda,na1,na2,na3,
     &                             j,k,1,jt)
	      do i=5, 2*(m1-mh)-1, 2
	         elnum=elnum+1
	         call element_ij(elnum, noda,na1,na2,na3,
     1                i,i,i, i+1,i+1,i+1,  i+2,i+2,i+2,
     2                j+2,j+1,j, j+2,j+1,j, j+2,j+1,j, k,k+1,k+2)
	      enddo
	      do i=2*(m1+mh)+1, ima-2, 2
	         elnum=elnum+1
	         call element_ij(elnum, noda,na1,na2,na3,
     1                i,i,i, i+1,i+1,i+1,  i+2,i+2,i+2,
     2                j+2,j+1,j, j+2,j+1,j, j+2,j+1,j, k,k+1,k+2)
	      enddo
	      estk_a(2,(k+1)/2,(j+1)/2)=elnum
	   enddo
C . . ZON B
	   do i=1, imb-2, 2
	      estk_b(1,(k+1)/2,(i+1)/2)=elnum+1
	      do j=1, jmb-2, 2
	         elnum=elnum+1
	         call element_ij(elnum, nodb,nb1,nb2,nb3,
     1                i+2,i+1,i, i+2,i+1,i,  i+2,i+1,i,
     2                j+2,j+2,j+2, j+1,j+1,j+1, j,j,j, k,k+1,k+2)
	      enddo
	      estk_b(2,(k+1)/2,(i+1)/2)=elnum
	   enddo
	enddo
	if (rtype.eq.0) then
	   kend=kma-2
	else
	   kend=kar1-4
	endif
	do k=2*mv+1, kend, 2
C . . ZON A
	   jt=jt0
	   do j=1, jma-2, 2
	      estk_a(1, (k+1)/2, (j+1)/2)=elnum+1
	      if (j.le.jma-6) call element_au(elnum, noda,na1,na2,na3,
     &                             j,k,1,jt)
	      do i=5, ima-2, 2
	         elnum=elnum+1
	         call element_ij(elnum, noda,na1,na2,na3,
     1                i,i,i, i+1,i+1,i+1, i+2,i+2,i+2,
     2                j+2,j+1,j, j+2,j+1,j, j+2,j+1,j, k,k+1,k+2)
	      enddo
	      estk_a(2, (k+1)/2, (j+1)/2)=elnum
	   enddo
C . . ZON B
	   do i=1, imb-2, 2
	      estk_b(1,(k+1)/2,(i+1)/2)=elnum+1
	      do j=1, jmb-2, 2
	         elnum=elnum+1
	         call element_ij(elnum, nodb,nb1,nb2,nb3,
     1                i+2,i+1,i, i+2,i+1,i,  i+2,i+1,i,
     2                j+2,j+2,j+2, j+1,j+1,j+1, j,j,j, k,k+1,k+2)
	      enddo
	      estk_b(2,(k+1)/2,(i+1)/2)=elnum
	   enddo
	enddo
C..MESH REDUCTION >= 1 :
	if (rtype.ge.1) then
	   if (mod(ma,4).eq.0) then
	      dk=1
	      it0=0
	   else
	      dk=2
	      it0=2
	   endif
	   jt=jt0
	   k=kar1-2
C . . ZON A
	   do j=1, jma-2, 2
	      estk_a(1,(k+1)/2,(j+1)/2)=elnum+1
	      if (j.le.jma-6) call element_au(elnum, noda,na1,na2,na3,
     &                                        j,k,dk,jt)
	      it=it0
	      do i=5, ima-2, 2
	         elnum=elnum+1
	         call element_red1_ai(elnum, noda,na1,na2,na3,
     &                                i,j,k,it)
	      enddo
	      estk_a(2,(k+1)/2,(j+1)/2)=elnum
	   enddo
C . . ZON B
	   do i=1, imb-2, 2
	      estk_b(1,(k+1)/2,(i+1)/2)=elnum+1
	      do j=1, jmb-2, 2
	         elnum=elnum+1
	         call element_ij(elnum, nodb,nb1,nb2,nb3,
     1                i+2,i+1,i, i+2,i+1,i,  i+2,i+1,i,
     2                j+2,j+2,j+2, j+1,j+1,j+1, j,j,j, k,k+2,k+4)
	      enddo
	      estk_b(2,(k+1)/2,(i+1)/2)=elnum
	   enddo
C . . ZON A
	   jt=jt0
	   k=kar1
	   do j=1, jma-2, 2
	      estk_a(1,(k+1)/2,(j+1)/2)=elnum+1
	      if ( (j.le.jma-6).and.(dk.eq.1) )
     &             call element_au(elnum,noda,na1,na2,na3,j,k,dk,jt)
	      it=it0
	      do i=5, ima-2, 4
	         elnum=elnum+1
	         call element_red1_aii(elnum,noda,na1,na2,na3,i,j,k,it)
	      enddo
	      estk_a(2,(k+1)/2,(j+1)/2)=elnum
	   enddo
	endif
C..MESH REDUCTION = 1 :
	if (rtype.eq.1) then
	   do k=kar1+2, kma-2, 2
C . . . . ZON A
	      jt=jt0
	      do j=1, jma-2, 2
	         estk_a(1,(k+1)/2,(j+1)/2)=elnum+1
	         if (j.le.jma-6) call element_au(elnum,
     &                                noda,na1,na2,na3, j,k,1,jt)
	         do i=5, ima-2, 4
	            elnum=elnum+1
	            call element_ij(elnum, noda,na1,na2,na3,
     1                   i,i,i, i+2,i+2,i+2,  i+4,i+4,i+4,
     2                   j+2,j+1,j, j+2,j+1,j, j+2,j+1,j, k,k+1,k+2)
	         enddo
	         estk_a(2, (k+1)/2, (j+1)/2)=elnum
	      enddo
C . . . . ZON B
	      do i=1, imb-2, 2
	         estk_b(1, (k+1)/2, (i+1)/2)=elnum+1
	         do j=1, jmb-2, 2
	            elnum=elnum+1
	            call element_ij(elnum, nodb,nb1,nb2,nb3,
     1                   i+2,i+1,i, i+2,i+1,i, i+2,i+1,i,
     2                   j+2,j+2,j+2, j+1,j+1,j+1, j,j,j, k,k+1,k+2)
	         enddo
	         estk_b(2, (k+1)/2, (i+1)/2)=elnum
	      enddo
	   enddo
	endif
C..MESH REDUCTION = 2 :
	if (rtype.eq.2) then
	   k=kar2
	   jt=jt0
C . . . . ZON A
	   do j=3, jma, 4
	      estk_a(1, (k+1)/2, (j+1)/4)=elnum+1
	      if (j.le.jma-6)
     &           call element_red2_au(elnum, noda,na1,na2,na3,3,j,k,jt)
	      do i=7, ima, 4
	         if (mod(jt,4).eq.0) then
	            elnum=elnum+1
	            call element_jk(elnum, noda,na1,na2,na3,
     1                   i-2,i,i+2, j,j-1,j-2, j,j-1,j-2, j,j-1,j-2,
     2                   k-2,k-2,k-2,  k-1,k-1,k-1,  k,k,k)
	            elnum=elnum+1
	            call element_jk(elnum, noda,na1,na2,na3,
     1                   i-2,i,i+2, j+2,j+1,j, j+2,j+1,j, j+2,j+1,j,
     2                   k-2,k-2,k-2,  k,k,k-1,  k+2,k+1,k)
	            elnum=elnum+1
	            call element_jk(elnum, noda,na1,na2,na3,
     1                   i-2,i,i+2, j,j-1,j-2, j+1,j,j-2, j+2,j,j-2,
     2                   k,k,k,  k+1,k+1,k+1,  k+2,k+2,k+2)
	         else
	            elnum=elnum+1
	            call element_jk(elnum, noda,na1,na2,na3,
     1                   i-2,i,i+2, j+2,j+1,j, j+2,j+1,j, j+2,j+1,j,
     2                   k-2,k-2,k-2,  k-1,k-1,k-1,  k,k,k)
	            elnum=elnum+1
	            call element_jk(elnum, noda,na1,na2,na3,
     1                   i-2,i,i+2, j,j-1,j-2,  j,j-1,j-2, j,j-1,j-2,
     2                   k-2,k-2,k-2,  k-1,k,k,  k,k+1,k+2)
	            elnum=elnum+1
	            call element_jk(elnum, noda,na1,na2,na3,
     1                   i-2,i,i+2, j+2,j+1,j,  j+2,j,j-1, j+2,j,j-2,
     2                   k,k,k,  k+1,k+1,k+1,  k+2,k+2,k+2)
	         endif
	      enddo
	      jt=jt+2
	      estk_a(2, (k+1)/2, (j+1)/4)=elnum
	   enddo
C . . . . ZON B
	   if (mod(nb,4).ne.0) then
	      jt0=0
	   else
	      jt0=2
	   endif
	   do i=1, imb-2, 2
	      jt=jt0
	      estk_b(1, (k+1)/2, (i+1)/2)=elnum+1
	      do j=3, jmb-2, 4
	         if (mod(jt,4).eq.0) then
	            elnum=elnum+1
	            call element_jk2(elnum, nodb,nb1,nb2,nb3,
     1                   i,i+1,i+2, j,j,j,  j-1,j-1,j-1, j-2,j-2,j-2,
     2                   k-2,k-1,k, k-2,k-1,k, k-2,k-1,k)
	            elnum=elnum+1
	            call element_jk2(elnum, nodb,nb1,nb2,nb3,
     1                   i,i+1,i+2, j+2,j+2,j+2,  j+1,j+1,j+1, j,j,j,
     2                   k-2,k,k+2, k-2,k,k+1, k-2,k-1,k)
	            elnum=elnum+1
	            call element_jk2(elnum, nodb,nb1,nb2,nb3,
     1                   i,i+1,i+2, j,j+1,j+2, j-1,j,j, j-2,j-2,j-2,
     2                   k,k+1,k+2, k,k+1,k+2, k,k+1,k+2)
	            jt=jt+2
	         else
	            elnum=elnum+1
	            call element_jk2(elnum, nodb,nb1,nb2,nb3,
     1                   i,i+1,i+2, j+2,j+2,j+2, j+1,j+1,j+1, j,j,j,
     2                   k-2,k-1,k, k-2,k-1,k, k-2,k-1,k)
	            elnum=elnum+1
	            call element_jk2(elnum, nodb,nb1,nb2,nb3,
     1                   i,i+1,i+2, j,j,j,  j-1,j-1,j-1, j-2,j-2,j-2,
     2                   k-2,k-1,k, k-2,k,k+1, k-2,k,k+2)
	            elnum=elnum+1
	            call element_jk2(elnum, nodb,nb1,nb2,nb3,
     1                   i,i+1,i+2, j+2,j+2,j+2,  j+1,j,j, j,j-1,j-2,
     2                   k,k+1,k+2, k,k+1,k+2, k,k+1,k+2)
	            jt=jt+2
	         endif
	      enddo
	      estk_b(2, (k+1)/2, (i+1)/2)=elnum
	   enddo
C . . . . ZON A
	   do k=kar2+2, kma-2, 2
	      do j=1, jma-6, 4
	         estk_a(1, (k+1)/2, (j+3)/4)=elnum+1
	         do i=1, ima-2, 4
	            elnum=elnum+1
	            call element_ij(elnum, noda,na1,na2,na3,
     1                   i,i,i, i+2,i+2,i+2,  i+4,i+4,i+4,
     2                   j+4,j+2,j, j+4,j+2,j, j+4,j+2,j, k,k+1,k+2)
	         enddo
	         estk_a(2, (k+1)/2, (j+3)/4)=elnum
	      enddo
	      j=jma-4
	      estk_a(1, (k+1)/2, (j+3)/4)=elnum+1
	      do i=5, ima-2, 4
	        elnum=elnum+1
	        call element_ij(elnum, noda,na1,na2,na3,
     1               i,i,i, i+2,i+2,i+2,  i+4,i+4,i+4,
     2               j+4,j+2,j, j+4,j+2,j, j+4,j+2,j, k,k+1,k+2)
	      enddo
	      estk_a(2, (k+1)/2, (j+3)/4)=elnum
C . . . . ZON B
	      do i=1, imb-2, 2
	         estk_b(1, (k+1)/2, (i+1)/2)=elnum+1
	         do j=1, jmb-2, 4
	            elnum=elnum+1
	            call element_ij(elnum, nodb,nb1,nb2,nb3,
     1                   i+2,i+1,i, i+2,i+1,i, i+2,i+1,i,
     2                   j+4,j+4,j+4, j+2,j+2,j+2, j,j,j, k,k+1,k+2)
	         enddo
	         estk_b(2, (k+1)/2, (i+1)/2)=elnum
	      enddo
	   enddo
	endif
 
	write(*,'(t15,a,i5,a/)')    '=> ', elnum,
     &            ' number of elements has been defined !'
 
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine element_ij(elnum,nod,n1,n2,n3,i1,i2,i3,i4,i5,i6,
     &             i7,i8,i9, j1,j2,j3,j4,j5,j6,j7,j8,j9, k1,k2,k3)
c-- The routine assigns nodes to an element in Zon S.
	implicit none
c
        include 'mesh3d_scp_common_eln.f'
c
	integer n1,n2,n3,nod(n1,n2,n3),elnum,
     1          i1,i2,i3,i4,i5,i6,i7,i8,i9,
     2          j1,j2,j3,j4,j5,j6,j7,j8,j9, k1,k2,k3
c
	eln(elnum,0)=elnum
	eln(elnum,1) = nod(i7,j7,k1)
	eln(elnum,2) = nod(i7,j7,k3)
	eln(elnum,3) = nod(i9,j9,k3)
	eln(elnum,4) = nod(i9,j9,k1)
 
	eln(elnum,5) = nod(i1,j1,k1)
	eln(elnum,6) = nod(i1,j1,k3)
	eln(elnum,7) = nod(i3,j3,k3)
	eln(elnum,8) = nod(i3,j3,k1)
 
	eln(elnum,9) =  nod(i7,j7,k2)
	eln(elnum,10) = nod(i8,j8,k3)
	eln(elnum,11) = nod(i9,j9,k2)
	eln(elnum,12) = nod(i8,j8,k1)
 
	eln(elnum,13) = nod(i1,j1,k2)
	eln(elnum,14) = nod(i2,j2,k3)
	eln(elnum,15) = nod(i3,j3,k2)
	eln(elnum,16) = nod(i2,j2,k1)
 
	eln(elnum,17) = nod(i4,j4,k1)
	eln(elnum,18) = nod(i4,j4,k3)
	eln(elnum,19) = nod(i6,j6,k3)
	eln(elnum,20) = nod(i6,j6,k1)
 
	eln(elnum,21) = nod(i5,j5,k2)
 
	eln(elnum,22) = nod(i4,j4,k2)
	eln(elnum,23) = nod(i5,j5,k3)
	eln(elnum,24) = nod(i6,j6,k2)
	eln(elnum,25) = nod(i5,j5,k1)
 
	eln(elnum,26) = nod(i8,j8,k2)
	eln(elnum,27) = nod(i2,j2,k2)
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine element_ik(elnum,nod,n1,n2,n3,i1,i2,i3,i4,i5,i6,
     &             i7,i8,i9,j1,j2,j3, k1,k2,k3,k4,k5,k6,k7,k8,k9)
C--- The routine assigns nodes to an element in Zon S.
        implicit none
c
      include 'mesh3d_scp_common_eln.f'
c
	integer n1,n2,n3,nod(n1,n2,n3),elnum,
     1          i1,i2,i3,i4,i5, i6,i7,i8,i9,
     2          j1,j2,j3, k1,k2,k3,k4,k5,k6,k7,k8,k9
c
	eln(elnum,0)=elnum
	eln(elnum,1)=nod(i7,j3,k7)
	eln(elnum,2)=nod(i9,j3,k9)
	eln(elnum,3)=nod(i9,j1,k9)
	eln(elnum,4)=nod(i7,j1,k7)
 
	eln(elnum,5)=nod(i1,j3,k1)
	eln(elnum,6)=nod(i3,j3,k3)
	eln(elnum,7)=nod(i3,j1,k3)
	eln(elnum,8)=nod(i1,j1,k1)
 
	eln(elnum,9) =nod(i8,j3,k8)
	eln(elnum,10)=nod(i9,j2,k9)
	eln(elnum,11)=nod(i8,j1,k8)
	eln(elnum,12)=nod(i7,j2,k7)
 
	eln(elnum,13)=nod(i2,j3,k2)
	eln(elnum,14)=nod(i3,j2,k3)
	eln(elnum,15)=nod(i2,j1,k2)
	eln(elnum,16)=nod(i1,j2,k1)
 
	eln(elnum,17)=nod(i4,j3,k4)
	eln(elnum,18)=nod(i6,j3,k6)
	eln(elnum,19)=nod(i6,j1,k6)
	eln(elnum,20)=nod(i4,j1,k4)
 
	eln(elnum,21)=nod(i5,j2,k5)
 
	eln(elnum,22)=nod(i5,j3,k5)
	eln(elnum,23)=nod(i6,j2,k6)
	eln(elnum,24)=nod(i5,j1,k5)
	eln(elnum,25)=nod(i4,j2,k4)
 
	eln(elnum,26)=nod(i8,j2,k8)
	eln(elnum,27)=nod(i2,j2,k2)
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine element_jk(elnum,nod,n1,n2,n3,i1,i2,i3, j1,j2,j3,
     &             j4,j5,j6, j7,j8,j9, k1,k2,k3,k4,k5,k6,k7,k8,k9)
C--- The routine assigns nodes to an element in Zon S.
	implicit none
c
        include 'mesh3d_scp_common_eln.f'
c
	integer n1,n2,n3,nod(n1,n2,n3),elnum, i1,i2,i3,
     &     j1,j2,j3,j4,j5,j6,j7,j8,j9,  k1,k2,k3,k4,k5,k6,k7,k8,k9
	eln(elnum,0)=elnum
	eln(elnum,1)=nod(i3,j1,k1)
	eln(elnum,2)=nod(i3,j7,k7)
	eln(elnum,3)=nod(i3,j9,k9)
	eln(elnum,4)=nod(i3,j3,k3)
 
	eln(elnum,5)=nod(i1,j1,k1)
	eln(elnum,6)=nod(i1,j7,k7)
	eln(elnum,7)=nod(i1,j9,k9)
	eln(elnum,8)=nod(i1,j3,k3)
 
	eln(elnum,9) =nod(i3,j4,k4)
	eln(elnum,10)=nod(i3,j8,k8)
	eln(elnum,11)=nod(i3,j6,k6)
	eln(elnum,12)=nod(i3,j2,k2)
 
	eln(elnum,13)=nod(i1,j4,k4)
	eln(elnum,14)=nod(i1,j8,k8)
	eln(elnum,15)=nod(i1,j6,k6)
	eln(elnum,16)=nod(i1,j2,k2)
 
	eln(elnum,17)=nod(i2,j1,k1)
	eln(elnum,18)=nod(i2,j7,k7)
	eln(elnum,19)=nod(i2,j9,k9)
	eln(elnum,20)=nod(i2,j3,k3)
 
	eln(elnum,21)=nod(i2,j5,k5)
 
	eln(elnum,22)=nod(i2,j4,k4)
	eln(elnum,23)=nod(i2,j8,k8)
	eln(elnum,24)=nod(i2,j6,k6)
	eln(elnum,25)=nod(i2,j2,k2)
 
	eln(elnum,26)=nod(i3,j5,k5)
	eln(elnum,27)=nod(i1,j5,k5)
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine element_jk2(elnum,nod,n1,n2,n3,i1,i2,i3, j1,j2,j3,
     &             j4,j5,j6, j7,j8,j9, k1,k2,k3,k4,k5,k6,k7,k8,k9)
C--- The routine assigns nodes to an element in Zon S.
	implicit none
c
        include 'mesh3d_scp_common_eln.f'
c
	integer n1,n2,n3,nod(n1,n2,n3),elnum, i1,i2,i3,
     &      j1,j2,j3,j4,j5,j6,j7,j8,j9,  k1,k2,k3,k4,k5,k6,k7,k8,k9
	eln(elnum,0)=elnum
	eln(elnum,1)=nod(i3,j7,k7)
	eln(elnum,2)=nod(i3,j9,k9)
	eln(elnum,3)=nod(i1,j9,k9)
	eln(elnum,4)=nod(i1,j7,k7)
 
	eln(elnum,5)=nod(i3,j1,k1)
	eln(elnum,6)=nod(i3,j3,k3)
	eln(elnum,7)=nod(i1,j3,k3)
	eln(elnum,8)=nod(i1,j1,k1)
 
	eln(elnum,9) =nod(i3,j8,k8)
	eln(elnum,10)=nod(i2,j9,k9)
	eln(elnum,11)=nod(i1,j8,k8)
	eln(elnum,12)=nod(i2,j7,k7)
 
	eln(elnum,13)=nod(i3,j2,k2)
	eln(elnum,14)=nod(i2,j3,k3)
	eln(elnum,15)=nod(i1,j2,k2)
	eln(elnum,16)=nod(i2,j1,k1)
 
	eln(elnum,17)=nod(i3,j4,k4)
	eln(elnum,18)=nod(i3,j6,k6)
	eln(elnum,19)=nod(i1,j6,k6)
	eln(elnum,20)=nod(i1,j4,k4)
 
	eln(elnum,21)=nod(i2,j5,k5)
 
	eln(elnum,22)=nod(i3,j5,k5)
	eln(elnum,23)=nod(i2,j6,k6)
	eln(elnum,24)=nod(i1,j5,k5)
	eln(elnum,25)=nod(i2,j4,k4)
 
	eln(elnum,26)=nod(i2,j8,k8)
	eln(elnum,27)=nod(i2,j2,k2)
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine element_ijk(elnum,nod,n1,n2,n3,
     &          i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,
     &          i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,
     &          j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15,j16,
     &          j17,j18,j19,j20,j21,j22,j23,j24,j25,j26,j27,
     &          k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,
     &          k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27)
C--- The routine assigns nodes to an element in Zon S.
	implicit none
c
        include 'mesh3d_scp_common_eln.f'
c
	integer n1,n2,n3,nod(n1,n2,n3),elnum,
     &          i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,
     &          i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27,
     &          j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15,j16,
     &          j17,j18,j19,j20,j21,j22,j23,j24,j25,j26,j27,
     &          k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,
     &          k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27
	eln(elnum,0)=elnum
	eln(elnum,1)=nod(i7,j7,k7)
	eln(elnum,2)=nod(i25,j25,k25)
	eln(elnum,3)=nod(i27,j27,k27)
	eln(elnum,4)=nod(i9,j9,k9)
 
	eln(elnum,5)=nod(i1,j1,k1)
	eln(elnum,6)=nod(i19,j19,k19)
	eln(elnum,7)=nod(i21,j21,k21)
	eln(elnum,8)=nod(i3,j3,k3)
 
	eln(elnum,9) =nod(i16,j16,k16)
	eln(elnum,10)=nod(i26,j26,k26)
	eln(elnum,11)=nod(i18,j18,k18)
	eln(elnum,12)=nod(i8,j8,k8)
 
	eln(elnum,13)=nod(i10,j10,k10)
	eln(elnum,14)=nod(i20,j20,k20)
	eln(elnum,15)=nod(i12,j12,k12)
	eln(elnum,16)=nod(i2,j2,k2)
 
	eln(elnum,17)=nod(i4,j4,k4)
	eln(elnum,18)=nod(i22,j22,k22)
	eln(elnum,19)=nod(i24,j24,k24)
	eln(elnum,20)=nod(i6,j6,k6)
 
	eln(elnum,21)=nod(i14,j14,k14)
 
	eln(elnum,22)=nod(i13,j13,k13)
	eln(elnum,23)=nod(i23,j23,k23)
	eln(elnum,24)=nod(i15,j15,k15)
	eln(elnum,25)=nod(i5,j5,k5)
 
	eln(elnum,26)=nod(i17,j17,k17)
	eln(elnum,27)=nod(i11,j11,k11)
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine element_au(elnum,nod,n1,n2,n3,j,k,dk,jt)
C--- The routine assigns nodes to elements in Zon A under the front.
	implicit none
	integer  n1,n2,n3,nod(n1,n2,n3),elnum,j,k,dk,jt,i
	i=1
	if (mod(jt,4).eq.0) then
	   elnum=elnum+1
	   call element_ij(elnum, nod,n1,n2,n3, i,i,i, i+1,i+1,i+1,
     &     i+2,i+2,i+2, j+4,j+2,j, j+3,j+2,j, j+2,j+1,j, k,k+dk,k+2*dk)
	   elnum=elnum+1
	   call element_ij(elnum, nod,n1,n2,n3,i+2,i+2,i+2, i+3,i+3,i+3,
     &     i+4,i+4,i+4, j+2,j+1,j, j+2,j+1,j, j+2,j+1,j, k,k+dk,k+2*dk)
	   jt=jt+1
	elseif (mod(jt,4).eq.1) then
	   elnum=elnum+1
	   call element_ij(elnum, nod,n1,n2,n3, i,i+1,i+2, i+2,i+2,i+3,
     &     i+4,i+4,i+4, j+2,j+1,j, j+2,j+1,j, j+2,j+1,j, k,k+dk,k+2*dk)
	   jt=jt+1
	elseif (mod(jt,4).eq.2) then
	   elnum=elnum+1
	   call element_ij(elnum, nod,n1,n2,n3, i+2,i+1,i, i+3,i+2,i+2,
     &     i+4,i+4,i+4, j+2,j+1,j, j+2,j+1,j, j+2,j+1,j, k,k+dk,k+2*dk)
	   jt=jt+1
	elseif (mod(jt,4).eq.3) then
	   elnum=elnum+1
	   call element_ij(elnum, nod,n1,n2,n3, i,i,i, i+1,i+1,i+1,
     &     i+2,i+2,i+2, j+2,j,j-2, j+2,j,j-1, j+2,j+1,j, k,k+dk,k+2*dk)
	   elnum=elnum+1
	   call element_ij(elnum,nod,n1,n2,n3,i+2,i+2,i+2, i+3,i+3,i+3,
     &     i+4,i+4,i+4, j+2,j+1,j, j+2,j+1,j, j+2,j+1,j, k,k+dk,k+2*dk)
	   jt=jt+1
	endif
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine element_red1_ai(elnum,nod,n1,n2,n3,i,j,k,it)
C--- The routine assigns nodes to elements in Zon A under the front.
	implicit none
	integer n1,n2,n3,nod(n1,n2,n3),elnum,i,j,k,it
	if (mod(it,4).eq.0) then
           call element_ik(elnum, nod,n1,n2,n3,
     1          i,i,i, i+1,i+1,i+1, i+2,i+2,i+2, j,j+1,j+2,
     2          k,k+1,k+2, k,k+1,k+2, k,k+1,k+2)
	   it=it+1
	elseif (mod(it,4).eq.1) then
           call element_ik(elnum, nod,n1,n2,n3,
     1          i,i,i, i+1,i+1,i+1, i+2,i+2,i+2, j,j+1,j+2,
     2          k,k+1,k+2, k,k+2,k+3, k,k+2,k+4)
	   it=it+1
	elseif (mod(it,4).eq.2) then
           call element_ik(elnum, nod,n1,n2,n3,
     1          i,i,i, i+1,i+1,i+1, i+2,i+2,i+2, j,j+1,j+2,
     2          k,k+2,k+4, k,k+2,k+3, k,k+1,k+2)
	   it=it+1
	elseif (mod(it,4).eq.3) then
           call element_ik(elnum, nod,n1,n2,n3,
     1          i,i,i, i+1,i+1,i+1, i+2,i+2,i+2, j,j+1,j+2,
     2          k,k+1,k+2, k,k+1,k+2, k,k+1,k+2)
	   it=it+1
	endif
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine element_red1_aii(elnum,nod,n1,n2,n3,i,j,k,it)
C--- The routine assigns nodes to elements in Zon A under the front.
	implicit none
	integer  n1,n2,n3,nod(n1,n2,n3),elnum,i,j,k,it
	if (mod(it,4).eq.0) then
           call element_ik(elnum, nod,n1,n2,n3,
     1          i,i,i,  i+1,i+2,i+2, i+2,i+3,i+4, j,j+1,j+2,
     2          k,k+1,k+2, k,k+1,k+2, k,k+1,k+2)
	   it=it+2
	elseif (mod(it,4).eq.2) then
           call element_ik(elnum, nod,n1,n2,n3,
     1          i+2,i+1,i,  i+3,i+2,i+2,  i+4,i+4,i+4, j,j+1,j+2,
     2          k,k+1,k+2, k,k+1,k+2, k,k+1,k+2)
	   it=it+2
	endif
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine element_red2_au(elnum,nod,n1,n2,n3,i,j,k,jt)
C--- The routine assigns nodes to elements in Zon A under the front
C--- for the mesh reduction 2 case.
	implicit none
	integer  n1,n2,n3,nod(n1,n2,n3),elnum,i,j,k,jt
	if (mod(jt,4).eq.0) then
	   elnum=elnum+1
	   call element_jk(elnum, nod,n1,n2,n3, i,i+1,i+2, j,j-1,j-2,
     &     j,j-1,j-2, j,j-1,j-2, k-2,k-2,k-2,  k-1,k-1,k-1,  k,k,k)
	   elnum=elnum+1
	   call element_ijk(elnum, nod,n1,n2,n3,
     &i-2,i-2,i-2, i-1,i-1,i-1, i,i,i, i-2,i-2,i-2, i-1,i-1,i-1, i,i,i,
     &i-2,i-2,i-2, i-1,i-1,i-1, i,i,i, j+2,j,j-2, j+1,j,j-2, j,j-1,j-2,
     &j+2,j,j-2, j+1,j,j-2, j,j-1,j-2, j+2,j,j-2, j+1,j,j-2, j,j-1,j-2,
     &k-2,k-2,k-2, k-2,k-2,k-2, k-2,k-2,k-2, k,k,k, k,k,k, k-1,k-1,k-1,
     &k+2,k+2,k+2, k+1,k+1,k+1, k,k,k)
	   elnum=elnum+1
	   call element_ijk(elnum, nod,n1,n2,n3,
     &i-2,i-1,i, i,i,i+1, i+2,i+2,i+2, i-2,i-1,i, i,i,i+1, i+2,i+2,i+2,
     &i-2,i-1,i, i,i,i+1, i+2,i+2,i+2, j+2,j+1,j, j+2,j+1,j, j+2,j+1,j,
     &j+2,j+1,j, j+2,j+1,j, j+2,j+1,j, j+2,j+1,j, j+2,j+1,j, j+2,j+1,j,
     &k-2,k-2,k-2, k-2,k-2,k-2, k-2,k-2,k-2, k,k,k-1, k,k,k-1, k,k,k-1,
     &k+2,k+1,k, k+2,k+1,k, k+2,k+1,k)
	   elnum=elnum+1
	   call element_ijk(elnum, nod,n1,n2,n3,
     &i,i,i, i+1,i+1,i+1, i+2,i+2,i+2, i-1,i-1,i-1, i,i,i, i+2,i+2,i+2,
     &i-2,i-2,i-2, i,i,i, i+2,i+2,i+2, j,j-1,j-2, j,j-1,j-2, j,j-1,j-2,
     &j+1,j,j-2, j+1,j,j-2, j+1,j,j-2, j+2,j,j-2, j+2,j,j-2, j+2,j,j-2,
     &k,k,k, k,k,k, k,k,k, k+1,k+1,k+1, k+1,k+1,k+1, k+1,k+1,k+1,
     &k+2,k+2,k+2, k+2,k+2,k+2, k+2,k+2,k+2)
	else
	   elnum=elnum+1
	   call element_jk(elnum, nod,n1,n2,n3, i,i+1,i+2, j+2,j+1,j,
     &     j+2,j+1,j, j+2,j+1,j, k-2,k-2,k-2,  k-1,k-1,k-1,  k,k,k)
	   elnum=elnum+1
	   call element_ijk(elnum, nod,n1,n2,n3,
     &i-2,i-2,i-2, i-1,i-1,i-1, i,i,i, i-2,i-2,i-2, i-1,i-1,i-1, i,i,i,
     &i-2,i-2,i-2, i-1,i-1,i-1, i,i,i, j+2,j,j-2, j+2,j,j-1, j+2,j+1,j,
     &j+2,j,j-2, j+2,j,j-1, j+2,j+1,j, j+2,j,j-2, j+2,j,j-1, j+2,j+1,j,
     &k-2,k-2,k-2, k-2,k-2,k-2, k-2,k-2,k-2, k,k,k, k,k,k, k-1,k-1,k-1,
     &k+2,k+2,k+2, k+1,k+1,k+1, k,k,k)
	   elnum=elnum+1
	   call element_ijk(elnum, nod,n1,n2,n3,
     &i,i-1,i-2, i+1,i,i, i+2,i+2,i+2, i,i-1,i-2, i+1,i,i, i+2,i+2,i+2,
     &i,i-1,i-2, i+1,i,i, i+2,i+2,i+2, j,j-1,j-2, j,j-1,j-2, j,j-1,j-2,
     &j,j-1,j-2, j,j-1,j-2, j,j-1,j-2, j,j-1,j-2, j,j-1,j-2, j,j-1,j-2,
     &k-2,k-2,k-2, k-2,k-2,k-2, k-2,k-2,k-2, k-1,k,k, k-1,k,k, k-1,k,k,
     &k,k+1,k+2, k,k+1,k+2, k,k+1,k+2)
	   elnum=elnum+1
	   call element_ijk(elnum, nod,n1,n2,n3,
     &i,i,i, i+1,i+1,i+1, i+2,i+2,i+2, i-1,i-1,i-1, i,i,i, i+2,i+2,i+2,
     &i-2,i-2,i-2, i,i,i, i+2,i+2,i+2, j+2,j+1,j, j+2,j+1,j, j+2,j+1,j,
     &j+2,j,j-1, j+2,j,j-1, j+2,j,j-1, j+2,j,j-2, j+2,j,j-2, j+2,j,j-2,
     &k,k,k, k,k,k, k,k,k, k+1,k+1,k+1, k+1,k+1,k+1, k+1,k+1,k+1,
     &k+2,k+2,k+2, k+2,k+2,k+2, k+2,k+2,k+2)
	endif
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
