c----67--1---------2---------3---------4---------5---------6---------712
c
c    File: plate_surface.f   eilatest modification April 22, 1998
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine gen_surface(efront,ecrsur,ecrtip,easym,ecsym,erem,
     &           esid,etop,nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &           nodb,nb1,nb2,nb3, elnum, iws,errorFlag)
c
c  This routine define all the element sets that is located on one of
c  the boundaries of the model.
c
c    1. find nodes on the boundary in question
c
c    2. locate the elements containing these nodes and determine the
c       specific element face.
c
c    Note that faces are defined according to scheme in WARP3D
c
c    Hexahedron (brick) element faces
c
c      WARP3D    ABAQUS    local nodes
c
c      face 1    face 1     1-4-3-2
c      face 2    face 2     5-6-7-8
c      face 3    face 3     1-2-6-5
c      face 4    face 5     3-4-8-7
c      face 5    face 4     3-7-6-2
c      face 6    face 6     1-5-8-4
c
      implicit none
c
      integer efront(2,*),ecrsur(2,*),ecrtip(2,*),easym(2,*),
     &        ecsym(2,*),erem(2,*),esid(2,*),etop(2,*),elnum,
     &        iws,errorFlag
      integer ns1,ns2,ns3,nods(ns1,ns2,ns3),
     &        na1,na2,na3,noda(na1,na2,na3),
     &        nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      integer nsize
      parameter (nsize = 15000)
      integer nset(nsize),nn
c      integer io
c
c  Open temporary file to write out the node sets
c
c      io = 23
c      open(unit=io,file='nset.dat',status='unknown')
c
c  face_front
c
ccc      write(io,'(t1,a/t10,a/t1,a)') '**','nfront','**'
ccc      write(iws,'(t2,a)') ' nfront '
      call nset_nfront(nset,nn, nods,ns1,ns2,ns3,noda,na1,na2,na3,
     &                          nodb,nb1,nb2,nb3)
      write(*,'(t2,a,i6)') ' nfront n=',nn
ccc      call write_set(nn,nset,io)
      call find_elem_face(elnum,efront,nset,nn,iws,errorFlag)
      write(*,'(t2,a,i6)') ' efront n=',efront(1,1)
c
c  face_crsur
c
ccc      write(io,'(t1,a/t10,a/t1,a)') '**','ncrsur','**'
ccc      write(iws,'(t2,a)') ' ncrsur '
      call nset_ncrsur(nset,nn,nods,ns1,ns2,ns3,noda,na1,na2,na3)
      write(*,'(t2,a,i6)') ' ncrsur n=',nn
ccc      call write_set(nn,nset,io)
      call find_elem_face(elnum,ecrsur,nset,nn,iws,errorFlag)
      write(*,'(t2,a,i6)') ' ecrsur n=',ecrsur(1,1)
c
c  face_crtip  (element surfaces in the focused mesh at the crack tip)
c
ccc      write(io,'(t1,a/t10,a/t1,a)') '**','ncrtip','**'
ccc      write(iws,'(t2,a)') ' ncrtip '
      call nset_ntipall(nset,nn,nods,ns1,ns2,ns3)
      write(*,'(t2,a,i6)') ' ncrtip n=',nn
ccc      call write_set(nn,nset,io)
      call find_elem_face(elnum,ecrtip,nset,nn,iws,errorFlag)
      write(*,'(t2,a,i6)') ' ecrtip n=',ecrtip(1,1)
c
c  face_asym
c
ccc      write(io,'(t1,a/t10,a/t1,a)') '**','nasym','**'
ccc      write(iws,'(t2,a)') ' nasym '
      call nset_nasym(nset,nn,nods,ns1,ns2,ns3, noda,na1,na2,na3)
      write(*,'(t2,a,i6)') ' nasym n=',nn
ccc      call write_set(nn,nset,io)
      call find_elem_face(elnum,easym,nset,nn,iws,errorFlag)
      write(*,'(t2,a,i6)') ' easym n=',easym(1,1)
c
c  face_rem
c
ccc      write(io,'(t1,a/t10,a/t1,a)') '**','nrem','**'
ccc      write(iws,'(t2,a,i6)') ' nrem '
      call nset_nrem(nset,nn,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
      write(*,'(t2,a,i6)') ' nrem n=',nn
ccc      call write_set(nn,nset,io)
      call find_elem_face(elnum,erem,nset,nn,iws,errorFlag)
      write(*,'(t2,a,i6)') ' erem n=',erem(1,1)
c
c  face_sid
c
ccc      write(io,'(t1,a/t10,a/t1,a)') '**','nsid','**'
ccc      write(iws,'(t2,a)') ' nsid '
      call nset_nsid(nset,nn, nodb,nb1,nb2,nb3)
      write(*,'(t2,a,i6)') ' nsid n=',nn
ccc      call write_set(nn,nset,io)
      call find_elem_face(elnum,esid,nset,nn,iws,errorFlag)
      write(*,'(t2,a,i6)') ' esid n=',esid(1,1)
c
c  face_csym
c
ccc      write(io,'(t1,a/t10,a/t1,a)') '**','ncsym','**'
ccc      write(iws,'(t2,a)') ' ncsym '
      call nset_ncsym(nset,nn,nods,ns1,ns2,ns3,noda,na1,na2,na3,
     &                        nodb,nb1,nb2,nb3)
      write(*,'(t2,a,i6)') ' ncsym n=',nn
ccc      call write_set(nn,nset,io)
      call find_elem_face(elnum,ecsym,nset,nn,iws,errorFlag)
      write(*,'(t2,a,i6)') ' ecsym n=',ecsym(1,1)
c
c  face_top
c
ccc      write(io,'(t1,a/t10,a/t1,a)') '**','ntop','**'
ccc      write(iws,'(t2,a)') ' ntop '
      call nset_ntop(nset,nn,noda,na1,na2,na3,nodb,nb1,nb2,nb3)
      write(*,'(t2,a,i6)') ' ntop n=',nn
ccc      call write_set(nn,nset,io)
      call find_elem_face(elnum,etop,nset,nn,iws,errorFlag)
      write(*,'(t2,a,i6)') ' etop n=',etop(1,1)
c
ccc      close(io)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine find_elem_face(elnum,eset,nset,nn,iws,errorFlag)
c
c  Given a node set find the element faces that contain the nodes
c
      implicit none
      include 'plate_common_nod.f'
      include 'plate_common_eln.f'
      integer  elnum,eset(2,*),nset(*),nn,iws,errorFlag
      integer  ne,i,j,k,l,jv(3),eloc,ev(8),enp,nnp,face
      logical  ok
c
      data jv(1),jv(2),jv(3) / 1,3,5 /
c
      ne = 0
      do i=1, elnum
         do j=1,3
            eloc = jv(j)
            enp = nnr( eln( i, eloc ) )
            do k=1, nn
               nnp = nset(k)
               if ( enp .eq. nnp ) then
                  do l=1, 8
                     ev(l) = nnr( eln( i, l ) )
                  enddo
                  call elface(i,ev,nset,nn,eloc,ok,face,iws,errorFlag)
                  if (ok) then
                     ne = ne + 1
                     eset(1,ne) = i
                     eset(2,ne) = face
                  endif
                  goto 10
               endif
            enddo
         enddo
 10      continue
      enddo
c
c  Rearrange the eset.
c
      do i=ne+1, 2, -1
         eset(1,i) = eset(1,i-1)
         eset(2,i) = eset(2,i-1)
      enddo
      eset(1,1) = ne
      eset(2,1) = ne
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine elface(ei,ev,nset,nn,eloc,ok,face,iws,errorFlag)
c
c Given a node set find the element faces that contain the nodes.
c Note that in case an element only has two nodes on the surface
c it should not belong to the surface element set. In this case
c the variable face = 0.
c
      implicit none
      integer  ei,ev(8),nset(*),nn,eloc,iws,errorFlag
      integer  enp3,enp6,enp7,enp8,i,nnp,face
      logical  ok
c
      face = 0
c
c Local elem. node = 1 => face = 1, 3, or 6
c
      if ( eloc .eq. 1 ) then
         enp3 = ev(3) 
         enp6 = ev(6) 
         enp7 = ev(7) 
         enp8 = ev(8)
         do i=1, nn
            nnp = nset(i)
            if ( enp3 .eq. nnp ) then
               face = 1
               goto 10
            elseif ( enp6 .eq. nnp ) then 
               face = 3
               goto 10
            elseif ( enp8 .eq. nnp ) then 
               face = 6
               goto 10
            endif
         enddo
c
c Local elem. node = 3 => face = 4 or 5
c
      elseif ( eloc .eq. 3 ) then
         enp6 = ev(6) 
         enp8 = ev(8) 
         do i=1, nn
            nnp = nset(i)
            if ( enp6 .eq. nnp ) then
               face = 5
               goto 10
            elseif ( enp8 .eq. nnp ) then 
               face = 4
               goto 10
            endif
         enddo
c
c Local elem. node = 5 => face = 2
c
      elseif ( eloc .eq. 5 ) then
         enp7 = ev(7)
         do i=1, nn
            nnp = nset(i)
            if ( enp7 .eq. nnp ) then
               face = 2
               goto 10
            endif
         enddo
      endif
c
 10   continue
      call elface_check(face,ev,nset,nn,iws,errorFlag)
c
      if ( face .eq. 0 ) then
         ok = .false.
      else
         ok = .true.
      endif
ccc      if (ok) write(iws,'(t1,i5,a,i2,a,8i5)') ei,' face = ',face,
ccc     & ' ',(ev(i),i=1,8)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine elface_check(face,ev,nset,nn,iws,errorFlag)
c
c  This routine check if the surface is collapsed and if so if there
c  is another non-collapsed surface in the nset. 
c
      implicit none
      integer  face,ev(*),nset(*),nn,iws,errorFlag
      integer  nloc(8,6),etmp(8),newface,i,j,k
      logical  ok,new
c
      data  (nloc(j,1),j=1,8)/ 1, 4, 3, 2,   5, 6, 7, 8 /
      data  (nloc(j,2),j=1,8)/ 5, 6, 7, 8,   1, 2, 3, 4 /
      data  (nloc(j,3),j=1,8)/ 1, 2, 6, 5,   3, 4, 7, 8 /
      data  (nloc(j,4),j=1,8)/ 3, 4, 8, 7,   1, 2, 5, 6 /
      data  (nloc(j,5),j=1,8)/ 2, 3, 7, 6,   1, 4, 5, 8 /
      data  (nloc(j,6),j=1,8)/ 1, 5, 8, 4,   2, 3, 6, 7 /
c
      if ( face .eq. 0 ) return
      k = 0
      do i=1, 3
         do j=i+1, 4
            if ( ev(nloc(i,face)) .eq. ev(nloc(j,face)) ) k = k + 1
         enddo
      enddo
c
      if ( k .lt. 2 ) return
c
      do i=1, 4
         etmp( nloc(   i, face ) ) = 1
         etmp( nloc( i+4, face ) ) = 0
      enddo
c
      do i=1, nn
         do j=5, 8
            if ( nset(i) .eq. ev( nloc(j,face) ) ) then
               etmp( nloc( j, face ) ) = 2
            endif
         enddo
      enddo
      newface = 0
      do i=1, 6
         ok  = .true.
         new = .false.
         do j=1, 4
            if ( etmp( nloc(j,i) ) .lt. 1 ) ok   = .false.
            if ( etmp( nloc(j,i) ) .eq. 2 ) new  = .true.
         enddo
         if ( ok .and. new ) newface = i
      enddo
      if ( newface .eq. face ) then
         write(iws,101) '>> ERROR in the element face search routines'
         errorFlag = 1
      else
         face = newface
      endif
 101  format(t1,a)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine write_set(nn,nset,io)
      implicit none
      integer nset(*),nn,io,i
      write(io,'(10i6)') (nset(i),i=1,nn)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
