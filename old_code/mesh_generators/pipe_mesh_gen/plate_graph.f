c----67--1---------2---------3---------4---------5---------6---------7-2
c
c  File: plate_graph2.f  Written  December 9, 1997 by Jonas F.
c                        Modified July 6 1998	 J.F.
c
c----67--1---------2---------3---------4---------5---------6---------7-2
c
      subroutine gen_graph(efront,ecrsur,easym,ecsym,erem,esid,etop,
     &           tge,no_of_nodes,elnum,iws,jobname)
c
c Generate graphic files to be used in FECRACK user interface
c Six files are generated containing meshes of key model surfaces
c Face numbers according to definitions in WARP3D 
c
c Graphic file names:  graph1.plt  graph2.plt   graph3.plt
c                      graph4.plt  graph5.plt   graph6.plt
c                      graph7.plt
c
c A 3D-mesh Tecplot file is also generated
c
      implicit none
c
      integer efront(2,*),ecrsur(2,*),easym(2,*),ecsym(2,*),
     &        erem(2,*),esid(2,*),etop(2,*)
      integer tge, no_of_nodes,elnum,iws,noel
c
      integer    io
      parameter (io=21)
      integer    surface_no,face
      integer    select_face
c
ccc      character grpfile*10
      character grpfile*200,jobname*200
      logical   tecplot
c
c Generate the graphics files
c
c 1. Find out the total number of elements on all surfaces
c
c 2. Write out all the element surfaces
c
      noel = 0
      do surface_no = 1, 6
         face = select_face(surface_no,tge)
c face_front
         if ( face .eq. 1) noel = noel + efront(1,1)
c face_asym
         if ( face .eq. 2) noel = noel + easym(1,1)
c face_rem
         if ( face .eq. 3) noel = noel + erem(1,1)
c face_sid
         if ( face .eq. 4) noel = noel + esid(1,1)
c face_csym
         if ( face .eq. 5) noel = noel + ecsym(1,1)
c face_top
         if ( face .eq. 6) noel = noel + etop(1,1)
      enddo
c
      call appendFile(jobname,'_grp1.plt',grpfile)
      open(unit=io, file=grpfile, status='unknown')
      write(io,'(i5)') noel
c
      do surface_no = 1, 6
         face = select_face(surface_no,tge)
c
c face_front
         if ( face .eq. 1) then
            write(iws,101) ' Graph no.=',surface_no,' face = front'
            call write_out_surface(io, efront)
c
c face_asym
         elseif ( face .eq. 2) then
            write(iws,101) ' Graph no.=',surface_no,' face = asym'
            call write_out_surface(io, easym)
c
c face_rem
         elseif ( face .eq. 3) then
            write(iws,101) ' Graph no.=',surface_no,' face = rem'
            call write_out_surface(io, erem)
c
c face_sid
         elseif ( face .eq. 4) then
            write(iws,101) ' Graph no.=',surface_no,' face = sid'
            call write_out_surface(io, esid)
c
c face_csym
         elseif ( face .eq. 5) then
            write(iws,101) ' Graph no.=',surface_no,' face = csym'
            call write_out_surface(io, ecsym)
c
c face_top
         elseif ( face .eq. 6) then
            write(iws,101) ' Graph no.=',surface_no,' face = top'
            call write_out_surface(io, etop)
c
         endif
      enddo
      write(io,'(a3)')'END'
      close(io)
c
 101  format(t2,a,i2,a)
c
c   Put the crack surface in a separate graphics file.
c
      call appendFile(jobname,'_grp2.plt',grpfile)
      open(unit=io, file=grpfile, status='unknown')
      write(iws,101) ' Graph no.=',2,' face = crsur'
      write(io,'(i5)') ecrsur(1,1)
      call write_out_surface(io, ecrsur)
      write(io,'(a3)')'END'
      close(io)
c
      tecplot = .false.
c      tecplot = .true.
      if (tecplot) call gen_3dtecplot(no_of_nodes,elnum)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      integer function select_face(graph_no,tge)
c
c Selects which surface should be drawn depending on the iface & tge
c
c Choose face according to:
c
c              face(i,1,j) => face_front
c              face(i,2,j) => face_asym
c              face(i,3,j) => face_rem
c              face(i,4,j) => face_sid
c              face(i,5,j) => face_csym
c              face(i,6,j) => face_top
c
      implicit none
      integer graph_no,tge,face(6,6,10),select,i,j
c
      data (face(1,1,j),j=1,10) /101,201,202,301,302,0,0,0,0,0/
      data (face(1,2,j),j=1,10) /203,204,303,304,0,0,0,0,0,0/
      data (face(1,3,j),j=1,10) /305,306,0,0,0,0,0,0,0,0/
      data (face(1,4,j),j=1,10) /307,308,0,0,0,0,0,0,0,0/
      data (face(1,5,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(1,6,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/

c
      data (face(2,1,j),j=1,10) /203,204,303,304,307,308,0,0,0,0/
      data (face(2,2,j),j=1,10) /101,201,202,301,302,305,306,0,0,0/
      data (face(2,3,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(2,4,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(2,5,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(2,6,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
c
      data (face(3,1,j),j=1,10) /305,306,0,0,0,0,0,0,0,0/
      data (face(3,2,j),j=1,10) /307,308,0,0,0,0,0,0,0,0/
      data (face(3,3,j),j=1,10) /101,201,202,301,302,0,0,0,0,0/
      data (face(3,4,j),j=1,10) /203,204,303,304,0,0,0,0,0,0/
      data (face(3,5,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(3,6,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
c
      data (face(4,1,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(4,2,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(4,3,j),j=1,10) /203,204,303,304,307,308,0,0,0,0 /
      data (face(4,4,j),j=1,10) /101,201,202,301,302,305,306,0,0,0/
      data (face(4,5,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(4,6,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
c
      data (face(5,1,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(5,2,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(5,3,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(5,4,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(5,5,j),j=1,10) /101,201,203,301,303,305,307,0,0,0/
      data (face(5,6,j),j=1,10) /202,204,302,304,306,308,0,0,0,0/
c
      data (face(6,1,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(6,2,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(6,3,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(6,4,j),j=1,10) /0,0,0,0,0,0,0,0,0,0/
      data (face(6,5,j),j=1,10) /202,204,302,304,306,308,0,0,0,0/
      data (face(6,6,j),j=1,10) /101,201,203,301,303,305,307,0,0,0/
c
      select = 0
      do i=1, 6
         do j=1, 10
            if ( face(graph_no,i,j) .eq. tge ) select = i
         enddo
      enddo
c
      select_face = select
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine write_out_surface(io, eset)
c
c  Write out coordinates on graphics file unit = io
c
      implicit none
      include 'plate_common_nod.f'
      include 'plate_common_eln.f'
c
      integer  io, eset(2,*)
      integer  i,j,fnode(4,6),ne,ie,face_no,nodeptr
      double precision x,y,z
c
      data  (fnode(j,1),j=1,4)/ 1, 4, 3, 2 /
      data  (fnode(j,2),j=1,4)/ 5, 6, 7, 8 /
      data  (fnode(j,3),j=1,4)/ 1, 2, 6, 5 /
      data  (fnode(j,4),j=1,4)/ 3, 4, 8, 7 /
      data  (fnode(j,5),j=1,4)/ 2, 3, 7, 6 /
      data  (fnode(j,6),j=1,4)/ 1, 5, 8, 4 /
c
      ne = eset(1,1)
c
      do i=1, ne
         ie      = eset(1,i+1)
         face_no = eset(2,i+1)
cc         write(io,'(t1,a)') 'ZONE'
         do j=1, 4
            nodeptr = eln( ie, fnode(j,face_no) )
            x = npos( nodeptr, 1 )
            y = npos( nodeptr, 2 )
            z = npos( nodeptr, 3 )
            write(io,101) x, y, z
         enddo
      enddo
c
 101  format(t1,3g16.8)
c 102  format(t1,g16.8,', ',g16.8,', ',g16.8)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine gen_3dtecplot(no_of_nodes,elnum)
c
c Generate a 3D-mesh input deck to Tecplot
c
      implicit none
      include 'plate_common_eln.f'
      include 'plate_common_nod.f'
c
      integer no_of_nodes,elnum, io,i,j
      parameter (io = 22)
      character tecfile*11,form*16
c
      tecfile = 'plate3d.tec'
      open( unit=io, file=tecfile, status='unknown' )
c
      i = log10( real( no_of_nodes ) ) + 1
      j = log10( real( elnum ) ) + 1
      write(form,'(a,i1,a,i1,a)') '(t1,a,i',i,',a,i',j,',a)'
      write(io,form) 'ZONE T="SCP", N=',no_of_nodes,', E=',elnum,
     &               ', ET=BRICK, F=FEPOINT'
      do i=1, inm
         if (nnr(i).gt.0) write(io,101) ( npos(i,j), j=1, 3 )
      enddo
      do i=1, elnum
         write(io,102) ( nnr( eln(i,j) ), j=1, 8 )
      enddo
c
      close(io)
101   format( t1,2(g16.8,tr2),g16.8 )
102   format( t1,8i8 )
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
