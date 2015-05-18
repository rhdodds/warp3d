c
c----67--1---------2---------3---------4---------5---------6---------712
c
        integer function nodn(zon,i,j,k)
c-------------------------------------------------------
c  Funktionen r{knar ut ett nodnummer f|r en NOD
c-------------------------------------------------------
      implicit none
      integer  i,j,k, ns
      character zon*1
c
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer          kstart,kstep
      common  /nblock/ kstart,kstep
c
      if (zon.eq.'c') then
         nodn = (j-1)*imc*kmc + (k-1)*imc + i
      elseif (zon.eq.'s') then
         ns = jmc*imc*kmc
         nodn = ns + (j-1)*ims*kms + (k-1)*ims + i
      elseif (zon.eq.'a') then
         nodn = kstart + (k-1)*kstep + (j-1)*ima + i
      elseif (zon.eq.'b') then
         nodn = kstart + (k-1)*kstep + ima*jma + (j-1)*imb + i
      endif
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine nodnumber( nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3,
     &                      noda,na1,na2,na3, nodb,nb1,nb2,nb3 )
c-------------------------------------------------------
c  Rutinen skapar matriserna nodc, nods, noda & nodb
c-------------------------------------------------------
c
      implicit none
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      integer  i,j,k,ilocal,ktmp
c...Integer functions
      integer  nodn
      logical   left
c test variables
c     character svar*1
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer          kstart,kstep
      common  /nblock/ kstart,kstep
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
      kstart = jmc*imc*kmc + jms*ims*kms
      kstep  = ima*jma + imb*jmb
c                  ---------------------------------
c                             Zone C:
c                  ---------------------------------
      do j=1, jmc
         do k=1, kmc
            do i=1, imc
               nodc(i,j,k) = nodn('c',i,j,k)
            enddo
         enddo
      enddo
c
c Mesh coarsening, remove nodes
c
      do j=1, jmc
         do i=1, 4
            nodc(i,j,1) = 0
            nodc(i,j,2) = 0
         enddo
         nodc(imc-1,j,kmc-3) = 0
         nodc(imc,  j,kmc-3) = 0
         nodc(imc-1,j,kmc-1) = 0
         nodc(imc,  j,kmc-1) = 0
      enddo
c
      ilocal=1
      left=.true.
      call reduce_mesh_2_to_1(nodc,nc1,nc2,nc3,ilocal,3,imc-3,4,
     &                        1,jmc,1, kmc-2  ,left)
      do i=2, imc-3, 2
         do j=1, jmc
            nodc(i,j,kmc) = 0
         enddo
      enddo
c                  ---------------------------------
c                             Zone S:
c                  ---------------------------------
      do j=1, jms
         do k=1, kms
            do i=1, ims
               nods(i,j,k) = nodn('s',i,j,k)
            enddo
         enddo
      enddo
      if (sfred_type.eq.2) then
         ilocal=1
         left = .true.
         call reduce_mesh_2_to_1(nods,ns1,ns2,ns3,ilocal,3,ims,4,
     &                           1,jms,1,ksr1,left)
         do k=ksr1+2, kms
            do i=2, ims, 2
               do j=1, jms
                  nods(i,j,k)=0
               enddo
            enddo
         enddo
      elseif (sfred_type.eq.3) then
         ilocal=1
         call reduce_plane_3_to_1(nods,ns1,ns2,ns3,ilocal,4,ims,6,
     &                            1,jms,1,ksr1)
         do k=ksr1+2, kms
            do i=4, ims, 3
               do j=1, jms
                  nods(i-1,j,k)=0
                  nods(i-2,j,k)=0
               enddo
            enddo
         enddo
      endif
c
      if (sjred_type.eq.2) then
         ilocal = 2
         left = .true.
         call reduce_mesh_2_to_1(nods,ns1,ns2,ns3,ilocal,1,ims,1,
     &                           3,jms,4,kms-2,left)
      elseif (sjred_type.eq.3) then
         ilocal = 2
         call reduce_plane_3_to_1(nods,ns1,ns2,ns3,ilocal,1,ims,1,
     &                            4,jms,6,kms-2)
      endif
c                  ---------------------------------
c                             Zone A:
c                  ---------------------------------
c...Node numbers in Zone A before the reduction of elements.
	do k=1, kma
	   do j=1, jma
	      do i=1, ima
                 noda(i,j,k)=nodn('a',i,j,k)	
	      enddo
	   enddo
	enddo
C...The small limited area under the crack front adjacent to the
C   free surface on the cracked side of the specimen.
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
C...Node numbers in Zone A after the reduction of elements.
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
C...The small limited area under the crack front
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
 
C...Delete the node numbers in crack Zone area. !
C . . Zone A:
      do k=1, 2*mv
         do j=1, jma
            do i=2*m1+2, 2*m1+4*mh
               noda(i,j,k)=0
            enddo
         enddo
      enddo
C . . Delete the node numbers for j > jma-4 & i < 5 ( ZON A )
      do k=1, kma
         do j=jma-3, jma
            do i=1,4
               noda(i,j,k)=0
            enddo
         enddo
      enddo
c                  ---------------------------------
c                             Zone B:
c                  ---------------------------------
c...Node numbers in Zone B before the reduction of elements.
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
     &                           3,jmb,4,kar2,left)
         do k=kar2+2, kma
            do j=2, jmb, 2
               do i=1, imb
                  nodb(i,j,k)=0
               enddo
            enddo
         enddo
      endif
c
c.......6) Fixa till de olika zonernas gemensamma nodnummer :
      call common_nodes(nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3, 
     &                  noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c--Test of Program ! -------------------------
c      call print_node_arrays( nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3,
c     &                        noda,na1,na2,na3, nodb,nb1,nb2,nb3 )
c 
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
	subroutine common_nodes(nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3, 
     &                          noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c-----------------------------------------------------------------C
c   Rutinen ser till att angr{nsande omr}den till Zon S ( Zon A)  C
c   resp Zon A (Zon B och ZonC) f}r Zon S resp Zon A :s nodnr     C
c   i de gemensamma noderna                                       C
c-----------------------------------------------------------------C
	implicit none
c
        integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1           ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2           na1,na2,na3,noda(na1,na2,na3),
     3           nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  i,j,k,is,ks,ia,ka,t,t1,n1,n2,ia1,ia2,n,ja,jst
c
      logical ok
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
c  Common nodes between Zone C and Zone S:
c
      do j=1, jmc
         is = 0
         ks = 1
         do k=1, kmc-1
            if (nodc(imc,j,k).gt.0) then
               is = is + 1
               nods(is,j,ks) = nodc(imc,j,k)
            endif
         enddo
c
         do i=imc, 1, -1
            if (nodc(i,j,kmc).gt.0) then
               is = is + 1
               nods(is,j,ks) = nodc(i,j,kmc)
            endif
         enddo
c
         do k=kmc-1, 3, -1
            if (nodc(1,j,k).gt.0) then
               is = is + 1
               nods(is,j,ks) = nodc(1,j,k)
            endif
         enddo
         if (is.ne.ims) then
            write(*,'(/t5,a,a/)') 'ERROR: Common Node # between',
     1         ' Zone C and Zone S are incompatible, ',
     2          'see subroutine common_nodes(...)'
            write(*,*)'   is =' ,is,'  ims = ',ims
            stop
         endif
      enddo
c
c  Common nodes between Zone S and Zone A:
c
      ia1  = 2*m1 + 1
      ia2  = 2*(m1+mh+mh) + 1
      n1 = 2*mv + 1
      n2 = 2*mv + 1 + (ia2-ia1)
      if (sjred_type.eq.1) then
         jst = 1
      elseif (sjred_type.eq.2) then
         jst = 2
      elseif (sjred_type.eq.3) then
         jst = 3
      else
         write(*,*) ' ERROR: sjred_type should be equal to 1, 2 or 3'
         stop
      endif
      ja = 0
      do j=1, jms, jst
         ja = ja + 1
         n = 0
         ka = 0
         ia = ia2
         ok = .false.
         do i=1, ims
            if (nods(i,j,kms).gt.0) then
               ok = .true.
               n = n + 1
               if (n.le.n1) then
                  ka = ka + 1
                  noda(ia,ja,ka) = nods(i,j,kms)
               elseif (n.le.n2) then
                  ia = ia - 1
                  noda(ia,ja,ka) = nods(i,j,kms)
               else
                  ka = ka - 1
                  noda(ia,ja,ka) = nods(i,j,kms)
               endif
            endif
         enddo
         if ( ok .and. (ka.ne.1) ) then
            write(*,'(/t5,a,a/)') 'ERROR: Common Node # between',
     1          ' Zone S and Zone A are incompatible, ',
     2          'see subroutine common_nodes(...)'
            write(*,'(t2,2(a,i4))') ' ka =',ka,'  n =',n
            stop
         endif
      enddo
      if (ja.ne.jma) then
         write(*,'(/t5,a/t10,a/)') 'ERROR: Common Node # between',
     1             ' Zone S and Zone A are incompatible, ',
     2             'see subroutine common_nodes(...)'
         write(*,'(t2,2(a,i4))') ' ja =',ja,'  jma =',jma
         stop
      endif
c
c  Common nodes within Zone A:
c
      if (m1.eq.2) then
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
c
c  Common nodes between Zone A and Zone B:
c
	t1=2*(na-nb)
 	do k=1, kma
	   do t=1, jmb
	      nodb(1,t,k)=noda(ima,t+t1,k)
	   enddo
	enddo
c
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine reduce_mesh_2_to_1(nod,n1,n2,n3,ilocal,
     &           i1,i2,ist,j1,j2,jst,k,left)
C  The routnine removes nodes in reduced mesh areas.
      implicit none
      integer  n1,n2,n3,nod(n1,n2,n3),ilocal,i1,i2,ist,j1,j2,jst,k,i,j
      logical  left
      if (ilocal.eq.1) then
         do i=i1, i2, ist
            if (left) then
               do j=j1, j2, jst
                  nod(i-1,j,k+1)=0
                  nod(i+1,j,k-1)=0
                  nod(i+2,j,k-1)=0
                  nod(i+2,j,k+1)=0
                  nod(i-1,j,k+2)=0
                  nod(i+1,j,k+2)=0
               enddo
               left=.false.
            else
               do j=j1, j2, jst
                  nod(i-2,j,k-1)=0
                  nod(i-2,j,k+1)=0
                  nod(i-1,j,k-1)=0
                  nod(i+1,j,k+1)=0
                  nod(i-1,j,k+2)=0
                  nod(i+1,j,k+2)=0
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
                  nod(i,j-1,k+2)=0
                  nod(i,j+1,k+2)=0
               enddo
               left=.false.
            else
               do i=i1, i2, ist
                  nod(i,j-2,k-1)=0
                  nod(i,j-2,k+1)=0
                  nod(i,j-1,k-1)=0
                  nod(i,j+1,k+1)=0
                  nod(i,j-1,k+2)=0
                  nod(i,j+1,k+2)=0
               enddo
               left=.true.
            endif
         enddo
      endif
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine reduce_plane_3_to_1(nod,n1,n2,n3,ilocal,
     &           i1,i2,ist,j1,j2,jst,k)
C  The routnine removes nodes in reduced mesh areas.
      implicit none
      integer  n1,n2,n3,nod(n1,n2,n3),ilocal,i1,i2,ist,j1,j2,jst,k,i,j
      if (ilocal.eq.1) then
         do i=i1, i2, ist
	    do j=j1, j2, jst
               nod(i-3,j,k-1)=0
               nod(i-2,j,k-1)=0
               nod(i+2,j,k-1)=0
               nod(i+3,j,k-1)=0
c
               nod(i-3,j,k+1)=0
               nod(i-1,j,k+1)=0
               nod(i+1,j,k+1)=0
               nod(i+3,j,k+1)=0
c
               nod(i-2,j,k+2)=0
               nod(i-1,j,k+2)=0
               nod(i+1,j,k+2)=0
               nod(i+2,j,k+2)=0
            enddo
         enddo
      else
	 do j=j1, j2, jst
            do i=i1, i2, ist
               nod(i,j-3,k-1)=0
               nod(i,j-2,k-1)=0
               nod(i,j+2,k-1)=0
               nod(i,j+3,k-1)=0
c
               nod(i,j-3,k+1)=0
               nod(i,j-1,k+1)=0
               nod(i,j+1,k+1)=0
               nod(i,j+3,k+1)=0
c
               nod(i,j-2,k+2)=0
               nod(i,j-1,k+2)=0
               nod(i,j+1,k+2)=0
               nod(i,j+2,k+2)=0
            enddo
         enddo
      endif
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
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
c----67--1---------2---------3---------4---------5---------6---------712
c
c
c----67--1---------2---------3---------4---------5---------6---------712
c
c----------------------------------------------------------------
c        NEDANSTAENDE SUBROUTIN/ER ANVANDS ENBART VID
c               TEST AV SJALVA PROGRAMMET  !!!!!
c----------------------------------------------------------------
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine print_node_arrays( nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3,
     &                        noda,na1,na2,na3, nodb,nb1,nb2,nb3 )
c
      implicit none
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  i,j,k,jend, itmp(200)
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer          kstart,kstep
      common  /nblock/ kstart,kstep
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
      open(unit=11,file='nodmat.dat',status='unknown')
      open(unit=12,file='nodmatpr.dat',status='unknown')
c
      write(11,'(//t10,a)') '  *** modell information *** '
      write(11,'(t2,11(a,i2,tr1))') 'mr=',mr,' mv=',mv,' mh=',mh,
     &  ' m1=',m1,' m2=',m2,' ma=',ma,' na=',na,' mb=',mb,
     &  ' nb=',nb,' lt=',lt,' lred=',lred
c
      write(11,'(/t2,3(a,i4))') ' imc=',imc,' jmc=',jmc,' kmc=',kmc
      write(11,'(/t2,3(a,i4))') ' ims=',ims,' jms=',jms,' kms=',kms
      write(11,'(/t2,3(a,i4))') ' ima=',ima,' jma=',jma,' kma=',kma
      write(11,'(/t2,3(a,i4))') ' imb=',imb,' jmb=',jmb,' kma=',kma
c
      write(11,'(/t2,7(a,i4))') 'ksr1=',ksr1,' kar1=',kar1,
     &            ' kar2=',kar2,' rtype=',rtype,' sfred=',sfred,
     &            ' sfred_type=',sfred_type,' sjred_type=',sjred_type
c
      write(12,'(//t10,a)') '  *** modell information *** '
      write(12,'(t2,11(a,i2,tr1))') 'mr=',mr,' mv=',mv,' mh=',mh,
     &  ' m1=',m1,' m2=',m2,' ma=',ma,' na=',na,' mb=',mb,
     &  ' nb=',nb,' lt=',lt,' lred=',lred
c
      write(12,'(/t2,3(a,i4))') ' imc=',imc,' jmc=',jmc,' kmc=',kmc
      write(12,'(/t2,3(a,i4))') ' ims=',ims,' jms=',jms,' kms=',kms
      write(12,'(/t2,3(a,i4))') ' ima=',ima,' jma=',jma,' kma=',kma
      write(12,'(/t2,3(a,i4))') ' imb=',imb,' jmb=',jmb,' kma=',kma
c
      write(12,'(/t2,7(a,i4))') 'ksr1=',ksr1,' kar1=',kar1,
     &            ' kar2=',kar2,' rtype=',rtype,' sfred=',sfred,
     &            ' sfred_type=',sfred_type,' sjred_type=',sjred_type
c
      write(11,'(//t10,a)') '=> matris nodc(i,j,k) :'
      write(11,'(t10,a/)')  '   ===================='
      write(12,'(//t10,a)') '=> matris nodc(i,j,k) :'
      write(12,'(t10,a/)')  '   ===================='
      do j=1, 11
         write(11,'(t1,a,i5)') '*** jc=',j
         write(12,'(t1,a,i5)') '*** jc=',j
         do k=kmc, 1, -1
            do i=1, imc
               if (nodc(i,j,k).gt.0) then
                  itmp(i)=1
               else
                  itmp(i)=0
               endif
            enddo
            write(11,'(t1,120i1)') ( itmp(i),i=1, imc )
            write(12,'(t1,a,i5)')  '    k=',k
            write(12,'(t1,10i7)')  ( nodc(i,j,k), i=1, imc )
         enddo
      enddo
c
      write(11,'(//t10,a)') '=> matris nods(i,j,k) :'
      write(11,'(t10,a/)') '   ===================='
      write(12,'(//t10,a)') '=> matris nods(i,j,k) :'
      write(12,'(t10,a/)') '   ===================='
      write(11,'(//t10,a//)') '***********  i-k plane  *************'
      write(12,'(//t10,a//)') '***********  i-k plane  *************'
      do j=1, 11
         write(11,'(t1,a,i5)') '*** js=',j
         write(12,'(t1,a,i5)') '*** js=',j
         do k=kms, 1, -1
            do i=1, ims
               if (nods(i,j,k).gt.0) then
                  itmp(i)=1
               else
                  itmp(i)=0
               endif
            enddo
            write(11,'(t1,120i1)') ( itmp(i), i=1, ims )
            write(12,'(t1,a,i5)')  '  ks=',k
            write(12,'(t1,10i7)')  ( nods(i,1,k), i=1, ims )
         enddo
      enddo
c
      write(11,'(//t10,a//)') '***********  k-j plane  *************'
      write(12,'(//t10,a//)') '***********  k-j plane  *************'
      jend = 72
      if (jend.gt.jms) jend = jms
      do i=1, 15
         write(11,'(t1,a,i5)') '*** is=',i
         write(12,'(t1,a,i5)') '*** is=',i
         do k=kms, 1, -1
            do j=1, 72
               if (nods(i,j,k).gt.0) then
                  itmp(j)=1
               else
                  itmp(j)=0
               endif
            enddo
            write(11,'(t1,120i1)') (itmp(j),j=1, jend)
            write(12,'(t1,a,i5)')  '  ks=',k
            write(12,'(t1,10i7)')  (nods(i,j,k) ,j=1, jend)
         enddo
      enddo
c
      write(11,'(//t10,a)') '=> matris noda(i,j,k) :'
      write(11,'(t10,a/)') '   ===================='
      write(12,'(//t10,a)') '=> matris noda(i,j,k) :'
      write(12,'(t10,a/)') '   ===================='
      write(11,'(//t10,a//)') '***********  i-j plane  *************'
      do k=1, kma
	 write(11,'(/t2,a,i2,a)') '* ka=',k,'  ====== i ========='
	 write(12,'(/t2,a,i2,a)') '* ka=',k,'  ====== i ========='
	 do j=1,jma
            do i=1, ima
               if (noda(i,j,k).ne.0) then
                  itmp(i)=1
               else
                  itmp(i)=0
               endif
            enddo
            write(11,'(t1,120i1)') ( itmp(i), i=1, ima )
            write(12,'(t1,a,i5)')  '            ja=',j
            write(12,'(t1,10i7)')  (noda(i,j,k) ,i=1, ima)
 	 enddo
      enddo
c
      write(11,'(//t10,a)') '=> matris nodb(i,j,k) :'
      write(11,'(t10,a/)') '   ===================='
      write(12,'(//t10,a)') '=> matris nodb(i,j,k) :'
      write(12,'(t10,a/)') '   ===================='
      write(11,'(//t10,a//)') '***********  i-j plane  *************'
      do k=1, kma
         write(11,'(/t2,a,i2,a)') '* kb=',k,'  ====== i ========='
         write(12,'(/t2,a,i2,a)') '* kb=',k,'  ====== i ========='
         do j=1,jmb
            do i=1, imb
               if (nodb(i,j,k).ne.0) then
                  itmp(i) = 1
               else
                  itmp(i) = 0
               endif
            enddo
            write(11,'(t1,120i1)') ( itmp(i), i=1, imb )
            write(12,'(t1,a,i5)')  '            jb=',j
            write(12,'(t1,10i7)')  (nodb(i,j,k) ,i=1, imb)
         enddo
      enddo
c
      close(11)
      close(12)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
