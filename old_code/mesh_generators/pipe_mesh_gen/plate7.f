c----67--1---------2---------3---------4---------5---------6---------712
c
c    File: plate7.f   latest modification Nov 29, 1997
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine sort_out_nodes(elnum,no_of_nodes,etyp,prog,
     1           nstk_s,nstk_a,nstk_b, nods,ns1,ns2,ns3,
     2           noda,na1,na2,na3, nodb,nb1,nb2,nb3, iws)
c
c The routine sorts out the uniqe nodes and remove the excessive nodes.
c
      implicit none
c
      include 'plate_common_eln.f'
      include 'plate_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3),iws
c
      integer  elnum,no_of_nodes, etyp,
     1         i,j,k,nn, el_nodes,n, nstart,
     2         nstk_s(3,*),nstk_a(3,*),nstk_b(3,*)
c integer fcn
      integer   numb_of_el_nodes
c
      character prog*20
      logical   basic_plate
c
      integer      mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
c
c Change the local node numbering order in elements if in-put
c deck is generated for ABAQUS or WARP3D (local node numbering
c scheme in default element follows that of ADINA).  
c
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
c
c If etyp < 27 remove excessive nodes 
c
      do i=1, elnum
         do j=etyp+1, 27
            eln(i,j)=0
         enddo
      enddo
c
c Sort out the unique nodes
c
      do i=0, inm
         nnr(i)=0
      enddo
      do i=1, elnum
         el_nodes = numb_of_el_nodes(i)
         do j=1, el_nodes
            nnr(eln(i,j))=1
         enddo
      enddo
c
c Remove jumps in node-numbers - ALL FE CODES
c
c (Note that the pointer to the nnr(i) and npos(i,k) is lost
c
c      if ( (prog.eq.'WARP3D') .or. (prog.eq.'PATRAN') ) then
         nn=0
c   Change the node numbers
         do i=1, inm
            if (nnr(i).gt.0) then
               nn=nn+1
               nnr(i)=nn
            endif
         enddo
c
c modified 11/29/97 eln(i,j) is pointing to nnr(k) containing
c                            the node number.
c
c   Change the node numbers in the elements
c         do i=1, elnum
c            do j=1, 27
c               if (eln(i,j).gt.0) then
c                  eln(i,j)=nnr(eln(i,j))
c               endif
c            enddo
c         enddo
c      else
c         nn=0
c         do i=1, inm
c            if (nnr(i).gt.0) then
c               nn=nn+1
c               nnr(i)=i
c            endif
c         enddo
c      endif
c
      write(*,'(t15,a,i5,a)') '=> # of nodes    = ',nn
      write(iws,'(t15,a,i5,a)') '=> # of nodes    = ',nn
      no_of_nodes = nn
c
c  Generate graphics file to view the basic plate model (X0,Y0,Z0)
c
c      basic_plate = .true.
      basic_plate = .false.
      if (basic_plate) then
         call plot_xz_y0(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                   nodb,nb1,nb2,nb3)
         call plot_yz_x0(nods,ns1,ns2,ns3, noda,na1,na2,na3)
         call plot_xz_ylt(noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      endif
c
c Node statistics
c
c  Zone S:
c
      do j=1, jms
         if (nnr(nods(1,j,2)).gt.0) then
            nstk_s(1,j) = nnr(nods(1,j,2))
         else
            nstk_s(1,j) = nnr(nods(1,j,3))
         endif
         nstk_s(2,j) = nnr(nods(ims,j,kms))
         n = 0
         do k=2, kms
            do i=1, ims
               if (nnr(nods(i,j,k)).gt.0) n=n+1
            enddo
         enddo
         nstk_s(3,j) = n
c         write(34,'(a,i3,3i7)') ' C: j=',j, (nstk_s(i,j),i=1,3)
      enddo
c
c  Zone A:
c
      nstart = nnr(nods(ims,jms,kms))
      do k=1, kma
         nstk_a(1,k) = nnr(noda(1,1,k))
         nstk_a(2,k) = nnr(noda(ima,jma,k))
         n = 0
         do j=1, jma
            do i=1, ima
               if (nnr(noda(i,j,k)).gt.nstart) n=n+1
            enddo
         enddo
         nstk_a(3,k) = n
c         write(34,'(a,i3,3i7)') ' A: k=',k, (nstk_a(i,k),i=1,3)
      enddo
c
c  Zone B:
c
      do k=1, kma
         if (nnr(noda(2,1,k)).gt.nstart) then
            nstk_b(1,k) = nnr(nodb(2,1,k))
         else
            nstk_b(1,k) = nnr(nodb(3,1,k))
         endif
         nstk_b(2,k) = nnr(nodb(imb,jmb,k))
         n = 0
         do j=1, jmb
            do i=1, imb
               if (nnr(nodb(i,j,k)).gt.nstart) n=n+1
            enddo
         enddo
         nstk_b(3,k) = n
c         write(34,'(a,i3,3i7)') ' B: k=',k, (nstk_b(i,k),i=1,3)
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine change_elnum_to_abaqus(i)
c
c The routine change the local element numbering from ADINA to ABAQUS
c
      implicit none
c
      include 'plate_common_eln.f'
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
c  The local nodes 17, 18, 19, 20 and 21 remain unchanged
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine change_elnum_to_warp3d(i)
c
c The routine change the local element numbering from ABAQUS to WARP3D
c
      implicit none
c
      include 'plate_common_eln.f'
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
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      integer function numb_of_el_nodes(i)
      implicit none
c
      include 'plate_common_eln.f'
c
      integer  i,j
      j=27
10    if (eln(i,j).eq.0) then
         j=j-1
         goto 10
      endif
      numb_of_el_nodes = j
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine plot_xz_y0(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                      nodb,nb1,nb2,nb3)
c
c Generate a input file for tecplot; plane of the crack
c
      implicit none
      include 'plate_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3),
     2         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  io
      parameter (io=31)
      integer  i,j,k, js,jb, ia1,ia2,i1
c
      character tecfile*40
c
      logical  left
c
      integer      mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      tecfile = 'plate_crp.plt'
      open( unit=io, file=tecfile, status = 'unknown' )
 101  format(t1,2g16.8)
 102  format(t1,a)
c
c Boundary
c
      write(io,102) 'ZONE'
      write(io,101) 0.0, 0.0
      write(io,101) 0.0, t
      write(io,101) w  , t
      write(io,101) w  , 0.0
      write(io,101) 0.0, 0.0
c
c Zone A - Zone S - Zone B; rho-direction
c
      ia1=2*m1+1
      ia2=2*(m1+mh+mh)+1
      do j=3, jma-2, 2
         write(io,102) 'ZONE'
c    Zone A
         if (mod(j,4).eq.1) then
            i1=1
         else
            i1=3
         endif
         do i=i1, ia1, 2
            if (noda(i,j,1).gt.0) then
               write(io,101) npos(noda(i,j,1),1),npos(noda(i,j,1),3)
            endif
         enddo
c    Zone S
         js = j*sjred_type - 1
         do k=kms-2, 1, -2
            write(io,101) npos(nods(ims,j,k),1),npos(nods(ims,j,k),3)
         enddo
         do k=3, kms-2, 2
            write(io,101) npos(nods(1,j,k),1),npos(nods(1,j,k),3)
         enddo
c    Zone A
         do i=ia2, ima, 2
            write(io,101) npos(noda(i,j,1),1), npos(noda(i,j,1),3)
         enddo
c    Zone B
         jb = j - (jma-jmb)
         if (jb .gt. 1) then
            do i=3, imb, 2
               write(io,101) npos(nodb(i,jb,1),1), npos(nodb(i,jb,1),3)
            enddo
         endif
      enddo
c
c Zone A; phi-direction
c
      write(io,102) 'ZONE'
      if (mod(na,4).ne.0) then
         left=.true.
         write(io,101) npos(noda(3,1,3),1), npos(noda(3,1,3),3)
      else
         left=.false.
         write(io,101) npos(noda(1,1,3),1), npos(noda(1,1,3),3)
      endif
      do j=3, jma-4, 4
         if (left) then
            write(io,101) npos(noda(3,j,3),1), npos(noda(3,j,3),3)
            write(io,101) npos(noda(1,j+2,3),1), npos(noda(1,j+2,3),3)
            left=.false.
         else
            write(io,101) npos(noda(3,j,3),1), npos(noda(3,j,3),3)
            write(io,101) npos(noda(3,j+2,3),1), npos(noda(3,j+2,3),3)
            left=.true.
         endif
      enddo
      do i=5, ima-2, 2
         if ( noda(i,1,1) .gt. 0 ) then
            write(io,102) 'ZONE'
            do j=1, jma, 2
               if ( noda(i,j,1) .gt. 0 ) then
                  write(io,101) npos(noda(i,j,1),1), npos(noda(i,j,1),3)
               endif
            enddo
         endif
      enddo
c
c Zone S; phi-direction
c
      do k=1, kms-4, 2
         write(io,102) 'ZONE'
         i = ims
         do j=1, jms, 2
            write(io,101) npos(nods(i,j,k),1), npos(nods(i,j,k),3)
         enddo
         i = 1
         write(io,102) 'ZONE'
         do j=1, jms, 2
            write(io,101) npos(nods(i,j,k),1), npos(nods(i,j,k),3)
         enddo
      enddo
c
      if (sjred_type.eq.1) then
         k = kms - 2
         i = ims
         write(io,102) 'ZONE'
         do j=1, jms, 2
            write(io,101) npos(nods(i,j,k),1), npos(nods(i,j,k),3)
         enddo
         i = 1
         write(io,102) 'ZONE'
         do j=1, jms, 2
            write(io,101) npos(nods(i,j,k),1), npos(nods(i,j,k),3)
         enddo
      elseif (sjred_type.eq.2) then
         k = kms - 2
         i = ims
         write(io,102) 'ZONE'
         write(io,101) npos(nods(i,1,k),1), npos(nods(i,1,k),3)
         do j=5, jms, 8
            write(io,101) npos(nods(i,j-2,k),1), npos(nods(i,j-2,k),3)
            write(io,101) npos(nods(i,j,k+2),1), npos(nods(i,j,k+2),3)
            if (j.le.(jms-4)) then
               write(io,101) npos(nods(i,j+2,k),1),npos(nods(i,j+2,k),3)
               write(io,101) npos(nods(i,j+4,k),1),npos(nods(i,j+4,k),3)
            endif
         enddo
c
         i = 1
         write(io,102) 'ZONE'
         write(io,101) npos(nods(i,1,k),1), npos(nods(i,1,k),3)
         do j=5, jms, 8
            write(io,101) npos(nods(i,j-2,k),1), npos(nods(i,j-2,k),3)
            write(io,101) npos(nods(i,j,k+2),1), npos(nods(i,j,k+2),3)
            if (j.le.(jms-4)) then
               write(io,101) npos(nods(i,j+2,k),1),npos(nods(i,j+2,k),3)
               write(io,101) npos(nods(i,j+4,k),1),npos(nods(i,j+4,k),3)
            endif
         enddo
      elseif (sjred_type.eq.3) then
         k = kms - 2
         i = ims
         write(io,102) 'ZONE'
         write(io,101) npos(nods(i,1,k+2),1), npos(nods(i,1,k+2),3)
         do j=7, jms, 6
            write(io,101) npos(nods(i,j-4,k),1), npos(nods(i,j-4,k),3)
            write(io,101) npos(nods(i,j-2,k),1), npos(nods(i,j-2,k),3)
            write(io,101) npos(nods(i,j,k+2),1), npos(nods(i,j,k+2),3)
         enddo
c
         k = kms - 2
         i = 1
         write(io,102) 'ZONE'
         write(io,101) npos(nods(i,1,k+2),1), npos(nods(i,1,k+2),3)
         do j=7, jms, 6
            write(io,101) npos(nods(i,j-4,k),1), npos(nods(i,j-4,k),3)
            write(io,101) npos(nods(i,j-2,k),1), npos(nods(i,j-2,k),3)
            write(io,101) npos(nods(i,j,k+2),1), npos(nods(i,j,k+2),3)
         enddo
      endif
c
c Zone B; Z-direction
c
      do i=1, imb-2, 2
         write(io,102) 'ZONE'
         do j=1, jmb, 2
            write(io,101) npos(nodb(i,j,1),1), npos(nodb(i,j,1),3)
         enddo
      enddo
c
      close(io)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine plot_yz_x0(nods,ns1,ns2,ns3, noda,na1,na2,na3)
c
c Generate a input file for tecplot; symmetry plane 
c
      implicit none
      include 'plate_common_nod.f'
c
      integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1         na1,na2,na3,noda(na1,na2,na3)
c
      integer  io
      parameter (io=31)
      integer  i,j,k, ia1,ia2,i1,ksr2,kend
c
      character tecfile*40
c
      integer      mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      tecfile = 'plate_sym.plt'
      open( unit=io, file=tecfile, status = 'unknown' )
 101  format(t1,2g16.8)
 102  format(t1,a)
c
c...Zone S:
c
      ksr2 = ksr1 + 2
      j = 1
c
      if (sfred_type.eq.1) then
c R-direction
         do i=1, ims, 2
            write(io,102) 'ZONE'
            do k=1, kms, 2
               write(io,101) npos(nods(i,j,k),2),npos(nods(i,j,k),3)
            enddo
         enddo
c Phi-direction
         do k=1, kms, 2
            write(io,102) 'ZONE'
            do i=1, ims, 2
               write(io,101) npos(nods(i,j,k),2),npos(nods(i,j,k),3)
            enddo
         enddo
c
      elseif (sfred_type.eq.2) then
c
         do i=3, ims-2, 2
            write(io,102) 'ZONE'
            if (mod(i,4).eq.1) then
               do k=1, kms, 2
                  write(io,101) npos(nods(i,j,k),2),npos(nods(i,j,k),3)
               enddo
            else
               do k=1, ksr1, 2
                  write(io,101) npos(nods(i,j,k),2),npos(nods(i,j,k),3)
               enddo
            endif
         enddo
c
         do k=1, ksr1-2, 2
            write(io,102) 'ZONE'
            do i=1, ims, 2
               write(io,101) npos(nods(i,j,k),2),npos(nods(i,j,k),3)
            enddo
         enddo
         k = ksr1
         write(io,102) 'ZONE'
         write(io,101) npos(nods(1,j,k),2),npos(nods(1,j,k),3)
         do i=5, ims-4, 8
            write(io,101) npos(nods(i-2,j,k),2),npos(nods(i-2,j,k),3)
            write(io,101) npos(nods(i,j,k+2),2),npos(nods(i,j,k+2),3)
            write(io,101) npos(nods(i+2,j,k),2),npos(nods(i+2,j,k),3)
            write(io,101) npos(nods(i+4,j,k),2),npos(nods(i+4,j,k),3)
         enddo
         do k=ksr2, kms, 2
            write(io,102) 'ZONE'
            do i=1, ims, 4
               write(io,101) npos(nods(i,j,k),2),npos(nods(i,j,k),3)
            enddo
         enddo
c
      elseif (sfred_type.eq.3) then
c
         do i=1, ims, 2
            write(io,102) 'ZONE'
            if (mod(i,6).eq.1) then
               do k=1, kms, 2
                  write(io,101) npos(nods(i,j,k),2),npos(nods(i,j,k),3)
               enddo
            else
               do k=1, ksr1, 2
                  write(io,101) npos(nods(i,j,k),2),npos(nods(i,j,k),3)
               enddo
            endif
         enddo
c
         do k=1, ksr1-2, 2
            write(io,102) 'ZONE'
            do i=1, ims, 2
               write(io,101) npos(nods(i,j,k),2),npos(nods(i,j,k),3)
            enddo
         enddo
c
         k = ksr1
         write(io,102) 'ZONE'
         write(io,101) npos(nods(1,j,k+2),2),npos(nods(1,j,k+2),3)
         do i=4, ims, 6
          write(io,101) npos(nods(i-1,j,k),2),npos(nods(i-1,j,k),3)
          write(io,101) npos(nods(i+1,j,k),2),npos(nods(i+1,j,k),3)
          write(io,101) npos(nods(i+3,j,k+2),2),npos(nods(i+3,j,k+2),3)
         enddo
         do k=ksr2, kms, 2
            write(io,102) 'ZONE'
            do i=1, ims, 3
               write(io,101) npos(nods(i,j,k),2),npos(nods(i,j,k),3)
            enddo
         enddo
      endif
c
c Zone A
c
      ia1=2*m1+1
      ia2=2*(m1+mh+mh)+1
c    Z-direction
      do k=1, 2*mv+1, 2
         write(io,102) 'ZONE'
         do i=1, ia1, 2
            if (noda(i,1,k).gt.0) then
               write(io,101) npos(noda(i,1,k),2),npos(noda(i,1,k),3)
            endif
         enddo
         write(io,102) 'ZONE'
         do i=ia2, ima, 2
            if (noda(i,1,k).gt.0) then
               write(io,101) npos(noda(i,1,k),2),npos(noda(i,1,k),3)
            endif
         enddo
      enddo
      do k=2*mv+3, kar1-2, 2
         write(io,102) 'ZONE'
         do i=1, ima, 2
             if (noda(i,1,k).gt.0) then
               write(io,101) npos(noda(i,1,k),2),npos(noda(i,1,k),3)
            endif
         enddo
      enddo
      if (rtype.eq.0) then
         do k=kar1, kma, 2
            write(io,102) 'ZONE'
            do i=1, ima, 2
               if (noda(i,1,k).gt.0) then
                 write(io,101) npos(noda(i,1,k),2),npos(noda(i,1,k),3)
               endif
            enddo
         enddo
      else
         if (mod(ma,4).eq.0) then
            write(io,102) 'ZONE'
          write(io,101) npos(noda(1,1,kar1),2),npos(noda(1,1,kar1),3)
          write(io,101) npos(noda(5,1,kar1),2),npos(noda(5,1,kar1),3)
          do i=9, ima, 8
               write(io,101) npos(noda(i-2,1,kar1),2),
     &                       npos(noda(i-2,1,kar1),3)
               write(io,101) npos(noda(i,1,kar1+2),2),
     &                       npos(noda(i,1,kar1+2),3)
               if ((i+4).lt.na1) then
                  if (noda(i+2,1,kar1).gt.0) then
                      write(io,101) npos(noda(i+2,1,kar1),2),
     &                              npos(noda(i+2,1,kar1),3)
                      write(io,101) npos(noda(i+4,1,kar1),2),
     &                              npos(noda(i+4,1,kar1),3)
                  endif
               endif
            enddo
         else
            write(io,102) 'ZONE'
            write(io,101) npos(noda(5,1,kar1+2),2),
     &                    npos(noda(5,1,kar1+2),3)
            do i=9, ima, 8
               write(io,101) npos(noda(i-2,1,kar1),2),
     &                       npos(noda(i-2,1,kar1),3)
               write(io,101) npos(noda(i,1,kar1),2),
     &                       npos(noda(i,1,kar1),3)
               write(io,101) npos(noda(i+2,1,kar1),2),
     &                       npos(noda(i+2,1,kar1),3)
               write(io,101) npos(noda(i+4,1,kar1+2),2),
     &                       npos(noda(i+4,1,kar1+2),3)
            enddo
         endif
         do k=kar1+2, kma, 2
            if((rtype.eq.2).and.(mod(na,4).ne.0).and.(k.eq.kar2))then
                i1=3
            elseif ( (rtype.eq.2).and.(k.eq.kar2) ) then
              i1=ima
            else
              i1=1
            endif
            write(io,102) 'ZONE'
            do i=i1, ima, 2
                if (noda(i,1,k).gt.0) then
                 write(io,101) npos(noda(i,1,k),2),npos(noda(i,1,k),3)
               endif
            enddo
         enddo
      endif
c    Y-direction
      write(io,102) 'ZONE'
      do k=1, kma, 2
         if (noda(1,1,k).gt.0) then
            write(io,101)  npos(noda(1,1,k),2),npos(noda(1,1,k),3)
         endif
      enddo
      if (rtype.lt.2) then
         if (mod(na,4).ne.0) then
            write(io,102) 'ZONE'
            do k=1, kma, 2
               if (noda(3,1,k).gt.0) then
                  write(io,101) npos(noda(3,1,k),2),npos(noda(3,1,k),3)
               endif
            enddo
         endif
      else
         if (mod(na,4).ne.0) then
            write(io,102) 'ZONE'
            do k=1, kar2, 2
               if (noda(3,1,k).gt.0) then
                  write(io,101) npos(noda(3,1,k),2),npos(noda(3,1,k),3)
                endif
            enddo
            write(io,101) npos(noda(1,1,kar2+2),2),
     &                    npos(noda(1,1,kar2+2),3)
         endif
      endif
      if (rtype.eq.0) then
         do i=5, ima, 2
            write(io,102) 'ZONE'
            do k=1, kma, 2
               if (noda(i,1,k).gt.0) then
                 write(io,101) npos(noda(i,1,k),2),npos(noda(i,1,k),3)
               endif
            enddo
         enddo
      else
         do i=5, ima, 2
            write(io,102) 'ZONE'
            if (mod(i,4).eq.1) then
               kend=kma
            else
               kend=kar1
            endif
             do k=1, kend, 2
               if (noda(i,1,k).gt.0) then
                 write(io,101) npos(noda(i,1,k),2),npos(noda(i,1,k),3)
               endif
            enddo
         enddo
      endif
c
      close(io)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine plot_xz_ylt(noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c Generate a input file for tecplot; remote plane opposing the crack
c
      implicit none
      include 'plate_common_nod.f'
c
      integer  na1,na2,na3,noda(na1,na2,na3),
     1         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  io
      parameter (io=31)
      integer  i,j, jb, di,dj,i1
      double precision  z1,z2
c
      character tecfile*40
c
      logical  left
c
      integer      mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mf,mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       ims,jms,kms,ima,jma,kma,imb,jmb
      common  /max/ ims,jms,kms,ima,jma,kma,imb,jmb
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
      double precision  t,w,c,a,alfa
      common /geom/     t,w,c,a,alfa
c
      tecfile = 'plate_bak.plt'
      open( unit=io, file=tecfile, status = 'unknown' )
 101  format(t1,2g16.8)
 102  format(t1,a)
c
c Boundary
c
      z2 = npos(noda(ima,1,kma),3)
      z1 = npos(noda(1,1,kma),3)
      write(io,102) 'ZONE'
      write(io,101) 0.0, z1
      write(io,101) 0.0, z2
      write(io,101) w  , z2
      write(io,101) w  , z1
      write(io,101) 0.0, z1
c
c Zone A - Zone S - Zone B; rho-direction
c
      if (rtype.le.1) then
         dj=2
      else
         dj=4
      endif
      do j=1+dj, jma, dj
         write(io,102) 'ZONE'
c    Zone A
         if ( (rtype.le.1).and.(mod(j,4).eq.3) ) then
            i1=3
         else
            i1=1
         endif
         do i=i1, ima
            if (noda(i,j,kma).gt.0) then
             write(io,101) npos(noda(i,j,kma),1),npos(noda(i,j,kma),3)
            endif
         enddo
c    Zone B
         jb = j - (jma-jmb)
         if (jb .gt. 1) then
            do i=3, imb, 2
               write(io,101) npos(nodb(i,jb,1),1), npos(nodb(i,jb,1),3)
            enddo
         endif
      enddo
c
c Zone A; phi-direction
c
      if (rtype.le.1) then
         write(io,102) 'ZONE'
           if (mod(na,4).ne.0) then
            left=.true.
            write(io,101) npos(noda(3,1,kma),1), npos(noda(3,1,kma),3)
         else
            left=.false.
            write(io,101) npos(noda(1,1,kma),1), npos(noda(1,1,kma),3)
         endif
         do j=3, jma-4, 4
            if (left) then
             write(io,101) npos(noda(3,j,kma),1),npos(noda(3,j,kma),3)
             write(io,101) npos(noda(1,j+2,kma),1),
     &                       npos(noda(1,j+2,kma),3)
               left=.false.
            else
             write(io,101) npos(noda(3,j,kma),1),npos(noda(3,j,kma),3)
              write(io,101) npos(noda(3,j+2,kma),1),
     &                      npos(noda(3,j+2,kma),3)
              left=.true.
            endif
         enddo
      endif
      if (rtype.eq.0) then
         di=2
      else
         di=4
      endif
      do i=5, ima, di
         write(io,102) 'ZONE'
         do j=1, jma
            if (noda(i,j,kma).gt.0) write(io,101) 
     &         npos(noda(i,j,kma),1),npos(noda(i,j,kma),3)
         enddo
      enddo
c
c Zone B; Z-direction
c
      do i=3, imb-2, 2
         write(io,102) 'ZONE'
            do j=1, jmb
               if (nodb(i,j,kma).gt.0) write(io,101)
     &            npos(nodb(i,j,kma),1), npos(nodb(i,j,kma),3)
         enddo
      enddo
c
      close(io)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
