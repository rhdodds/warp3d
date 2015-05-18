c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine sort_out_nodes(job,jobh,elnum,no_of_nodes,etyp,prog, 
     1       nstk_c,nstk_s,nstk_a,nstk_b, nodc,nc1,nc2,nc3, 
     2     nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c--- The routine sort out the uniqe nodes in the models
c
      implicit none
c
      include 'common_eln.f'
c
      include 'common_nod.f'
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  jobh,elnum,no_of_nodes, etyp,
     1         i,j,k,nn, el_nodes,n, io, nstart,
     2         nstk_c(3,*),nstk_s(3,*),nstk_a(3,*),nstk_b(3,*)

c
C  INTEGER FUNCTIONS
      integer  numb_of_el_nodes
      character job*40,prog*20,form*16,tec_file*30
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
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
c
c.... Determine the unique nodes
	do i=0, inm
	   nnr(i)=0
	enddo
	do i=1, elnum
	   el_nodes = numb_of_el_nodes(i)
	   do j=1, el_nodes
	      nnr(eln(i,j))=1
	   enddo
	enddo
C... Take away jumps in node-numbers if WARP3D analysis
        if ( (prog.eq.'WARP3D') .or. (prog.eq.'PATRAN') ) then
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
c
	write(*,'(t15,a,i5,a)') '=> the model contains ',
     &                           nn,' unique nodes'
        no_of_nodes = nn
c
c  Create a tec-plot file
c
      tec_file = job(1:jobh)//'.tec'
      io = 31
      i = log10(real(no_of_nodes)) + 1
      j = log10(real(elnum)) + 1
      write(form,'(a,i1,a,i1,a)') '(t1,a,i',i,',a,i',j,',a)'
      open(unit=io,file=tec_file,status='unknown')
      write(io,form) 'ZONE T="SCP", N=',no_of_nodes,
     &               ', E=',elnum,', ET=BRICK, F=FEPOINT'
      do i=1, inm
         if (nnr(i).gt.0) then
            write(io,'(t1,2(g16.8,tr2),g16.8)') ( npos(i,k), k=1, 3 )
         endif
      enddo
      do i=1, elnum
         write(io,'(t1,8i8)') ( eln(i,k), k=1, 8 )
      enddo
      close(31)
c
c... Create node statistics
c
c... Zone C:
      do j=1, jmc
         nstk_c(1,j) = nnr(nodc(5,j,1))
         nstk_c(2,j) = nnr(nodc(imc,j,kmc))
         n = 0
         do k=1, kmc
            do i=1, imc
               if (nnr(nodc(i,j,k)).gt.0) n=n+1
            enddo
         enddo
         nstk_c(3,j) = n
c         write(34,'(a,i3,3i7)') ' C: j=',j, (nstk_c(i,j),i=1,3)
      enddo
c
c... Zone S:
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
c... Zone A:
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
c... Zone B:
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
c-- The routine change the local element numbering from ADINA to ABAQUS
c
	implicit none
c
        include 'common_eln.f'
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
        include 'common_eln.f'
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
      subroutine write_out_nod_el(job,jobh,elnum,elnum_cell,etyp,prog)
c
c---- The routine writes out nod coordinates and element definitions
c---- on file.
c
      implicit none
c
      include 'common_eln.f'
c
      include 'common_nod.f'
c
      integer  elnum,elnum_cell,etyp, jobh,i,j
c
c  INTEGER FUNTIONS
      character job*40,fil1*40,fil2*40,fil3*40,prog*20
c
      write(*,'(/t5,a,a/)') '* writing node coord. and element ',
     &                      'def. on out-put files . . . .'
c
c.. Write out the nod coordinates
      if (prog.eq.'ABAQUS') then
         fil1=job(1:jobh)//'.015'
      elseif ((prog.eq.'WARP3D').or.(prog.eq.'PATRAN')) then
         fil1=job(1:jobh)//'.crd'
      endif
      open(unit=15,file=fil1,status='unknown')
      if (prog.eq.'ABAQUS') then
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
         fil3=job(1:jobh)//'.017'
         open(unit=17,file=fil3,status='unknown')
         if (etyp.eq.8) then
            do i=1, elnum_cell
               write(16,102) (eln(i,j), j=0, etyp)
            enddo
            do i=elnum_cell+1, elnum
               write(17,102) (eln(i,j), j=0, etyp)
            enddo
         elseif (etyp.eq.20) then
            do i=1, elnum_cell
               write(16,103) (eln(i,j), j=0, etyp)
            enddo
            do i=elnum_cell+1, elnum
               write(17,103) (eln(i,j), j=0, etyp)
            enddo
         else
            do i=1, elnum_cell
               write(16,104) (eln(i,j), j=0, etyp)
            enddo
            do i=elnum_cell+1, elnum
               write(17,104) (eln(i,j), j=0, etyp)
            enddo
         endif
         close(16)
         close(17)
      endif
C
 101  format(t3,i7,a,2(g18.10,a),g18.10)
 102  format( t1,i5,',',7(i6,','),i6)
 103  format( t1,i5,',',10(i6,',')/t2,9(i6,','),i6)
 104  format( t1,i5,',',10(i6,',')/t2,10(i6,',')/t2,9(i6,','),i6)
C
 201  format(t3,i6, 3(g20.10) )
 203  format( t1,9(i7) )
 204  format( t1,9(i7),','/t8,8(i7),','/t8,4(i7))
 205  format( t1,9(i7),','/t8,8(i7),','/t8,8(i7),','/t8,8(i7))
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      integer function numb_of_el_nodes(i)
      implicit none
c
      include 'common_eln.f'
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
