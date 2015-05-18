c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine abaqus_inp(job,jobh,etyp,   nea2,neb2, 
     1           estk_s,estk_a,estk_b,  nodc,nc1,nc2,nc3,
     2           nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c--The routine creates an indata file (job.inp) to an ABAQUS analysis.
c
      implicit none
c
      include 'common_nod.f'
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  nea2,neb2,estk_s(2,*),estk_a(2,nea2,*),estk_b(2,neb2,*)
c
      integer  jobh,etyp, io,i,j,k, no_of_jcont, e1,e2
c
      character job*40,filinp*40,fil1*40,fil2*40,fil3*40,
     &          fil4*40,fil5*40,fil6*40
      character*4 jco(100)
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer      imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
      fil1=job(1:jobh)//'.015'
      fil2=job(1:jobh)//'.016'
      fil3=job(1:jobh)//'.017'
      fil4=job(1:jobh)//'.018'
      fil5=job(1:jobh)//'.019'
      fil6=job(1:jobh)//'.020'
      filinp=job(1:jobh)//'.inp'
      open(unit=14,file=filinp,status='unknown')
      io=14
      write(io,'(t1,a)') '*HEADING'
      write(io,'(t1,a)') '>> Computational Cell Element Model <<'
      write(io,'(t1,a)') '**'
      write(io,'(t1,a)') '*PREPRINT,MODEL=NO,ECHO=NO'
      write(io,'(t1,a,a)') '*NODE,NSET=NALL,INPUT=',fil1
      write(io,'(t1,a,i2,a,a)') '*ELEMENT,ELSET=ELGURS,TYPE="',etyp,
     &                           '",INPUT=',fil2
      write(io,'(t1,a,i2,a,a)') '*ELEMENT,ELSET=ELREGL,TYPE="',etyp,
     &                           '",INPUT=',fil3
      write(io,'(t1,a,i1,a)') '*ELSET,ELSET=LAYER_S,GENERATE'
      write(io,'(t1,i5,a,i5)') estk_s(1,1), ',' ,estk_s(2,(jms-1)/2)
c
      do k=1, kma-2, 2
         if (estk_a(1,(k+1)/2,1).gt.0) then
            if ( ((k+1)/2).lt.10) then
               write(io,'(t1,a,i1,a)')
     &              '*ELSET,ELSET=LAYER0',(k+1)/2,',GENERATE'
            else
               write(io,'(t1,a,i2,a)')
     &              '*ELSET,ELSET=LAYER',(k+1)/2,',GENERATE'
            endif
            e1 = estk_a(1,(k+1)/2,1)
            i = imb-1
 10         continue
               i = i - 1
            if (i.gt.0) then
               if ( estk_b(2,(k+1)/2,(i+1)/2 ).le.0 ) goto 10
               e2 = estk_b(2,(k+1)/2,(i+1)/2 )
               goto 20
            else
               j=jma-1
 15            continue
                  j = j - 1
               if (estk_a(2,(k+1)/2,(j+1)/2).le.0) goto 15
               e2 = estk_a(2,(k+1)/2,(j+1)/2)
               goto 20
            endif
 20         continue
            write(io,'(t1,i5,a,i5)') e1,',',e2
         endif
      enddo
c
      write(io,'(t1,a)') '*SOLID SECTION,ELSET=ELREGL,MATERIAL=MAT1'
      write(io,'(t1,a)') '*MATERIAL,NAME=MAT1'
      write(io,'(t1,a)') '*ELASTIC'
      write(io,'(t1,a)') ' 200.0E+09, 0.3'
      write(io,'(t1,a)') '**PLASTIC'
c
      write(io,'(t1,a)') '*SOLID SECTION,ELSET=ELGURS,MATERIAL=MAT2'
      write(io,'(t1,a)') '*MATERIAL,NAME=MAT2'
      write(io,'(t1,a)') '*ELASTIC'
      write(io,'(t1,a)') ' 200.0E+09, 0.3'
      write(io,'(t1,a)') '**PLASTIC'
      write(io,'(t1,a,a)') '*POROUS METAL PLASTICITY,',
     &                      'RELATIVE DENSITY=0.001'
      write(io,'(t1,a)') ' 1.50, 1.0, 2.25'
c
      write(io,'(t1,a)') '*RESTART,WRITE,FREQ=999'
      write(io,'(t1,a)') '**'
      write(io,'(t1,a)') '** result - information **'
      write(io,'(t1,a)') '**'
c      write(io,'(t1,a)') '*NSET,NSET=NODPRD'
c      write(io,'(t1,a,i5)') '1,2,',noda(1,1,1)
c      write(io,'(t1,a)') '*NSET,NSET=NODPRF,GENERATE'
c      write(io,'(t1,i5,a,i5)')  nods(1,1,1),',',noda(ima,jma,13)
c      write(io,'(t1,a)') '*NSET,NSET=NODPRF'
c      write(io,'(t1,a)') 'NODPRF,NFRONT,NBAK,NCRSUR,NUND,NTOP'
c      write(io,'(t1,a)') '*ELSET,ELSET=ELPRD'
c      write(io,'(t1,a)') '1,2'
c      write(io,'(t1,a)') '*ELSET,ELSET=ELPRF'
c      write(io,'(t1,a,a)') 'LAYERCR,LAYER01,LAYER02,LAYER03,',
c     &                      'LAYER04,LAYER05,LAYER06'
c      write(io,'(t1,a)') '*EQUATION,INPUT=FIL6'
c      call constraint_equ_abaqus(ima,kma,fil6,26,noda,na1,na2,na3)
c      write(io,'(t1,a)') '*BOUNDARY'
c      write(io,'(t1,a)') ' NSYM,1'
c      write(io,'(t1,a)') ' NFRONT,2'
c
       if (ima.gt.0) then
          close(io)
          return
       endif
c
      write(io,'(t1,a)') '**'
      write(io,'(t1,a)') '** ****** s t e p = 1  ************'
      write(io,'(t1,a)') '**'
      write(io,'(t1,a)') '*STEP'
      write(io,'(t1,a)') '*STATIC'
      write(io,'(t1,a)') '*BOUNDARY'
      write(io,'(t1,i5,a,i2,a)') nnr(noda(1,1,kma)),  ',2,2, 1.e-3??'
      write(io,'(t1,i5,a,i2,a)') nnr(noda(ima,1,kma)),',2,2, 1.e-3??'
      write(io,'(t1,a)') '*NODE PRINT,NSET=NODPRD'
      write(io,'(t1,a)') ' U,RF'
      write(io,'(t1,a)') '*EL PRINT,ELSET=ELPRD'
      write(io,'(t1,a)') ' S,MISES'
      write(io,'(t1,a)') ' E'
      write(io,'(t1,a)') '***NODE FILE,NSET=NODPRF,FREQ=10'
      write(io,'(t1,a)') '** U,RF'
      write(io,'(t1,a)') '***EL FILE,ELSET=ELPRF,FREQ=10'
      write(io,'(t1,a)') '** S,SINV,ENER'
      write(io,'(t1,a)') '** E,PE'
      write(io,'(t1,a)') '*ENERGY FILE,FREQ=10'
      no_of_jcont=mr+m1-2
      write(io,'(t1,a,i2,a)') '*CONTOUR INTEGRAL,CONTOURS=',
     &  no_of_jcont,',TYPE=J,NORMAL,SYMM,OUTPUT=BOTH,FREQ=10'
      write(io,'(t1,a)') '0.0, -1.0, 0.0'
      do i=1, 9
         jco(i)='CR0'//char(48+i)
      enddo
      jco(10)='CR10'
      do i=11, 19
         jco(i)='CR1'//char(38+i)
      enddo
      jco(20)='CR20'
      do i=21, 29
         jco(i)='CR2'//char(28+i)
      enddo
      jco(30)='CR30'
      do i=31, 39
         jco(i)='CR3'//char(18+i)
      enddo
      jco(40)='CR40'
      do i=41, 49
         jco(i)='CR4'//char(8+i)
      enddo
      jco(50)='CR50'
      do i=51, 59
         jco(i)='CR5'//char(-2+i)
      enddo
      jco(60)='CR60'
      do i=61, 69
         jco(i)='CR6'//char(-12+i)
      enddo
      jco(70)='CR70'
      do i=71, 79
         jco(i)='CR7'//char(-22+i)
      enddo
      jco(80)='CR80'
      do i=81, 89
         jco(i)='CR8'//char(-32+i)
      enddo
      jco(90)='CR90'
      do i=91, 99
         jco(i)='CR9'//char(-42+i)
      enddo
      if (jms.le.16) then
         write(io,101) ( jco(j), j=1, jms )
      elseif (jms.le.32) then
         write(io,101) ( jco(j), j=1, 16 )
         write(io,101) ( jco(j), j=17, jms )
      elseif (jms.le.48) then
         write(io,101) ( jco(j), j=1, 16 )
         write(io,101) ( jco(j), j=17, 32 )
         write(io,101) ( jco(j), j=33, jms )
      elseif (jms.le.64) then
         write(io,101) ( jco(j), j=1, 16 )
         write(io,101) ( jco(j), j=17, 32 )
         write(io,101) ( jco(j), j=33, 48 )
         write(io,101) ( jco(j), j=49, jms )
      elseif (jms.le.80) then
         write(io,101) ( jco(j), j=1, 16 )
         write(io,101) ( jco(j), j=17, 32 )
         write(io,101) ( jco(j), j=33, 48 )
         write(io,101) ( jco(j), j=49, 64 )
         write(io,101) ( jco(j), j=65, jms )
      elseif (jms.le.80) then
         write(io,101) ( jco(j), j=1, 16 )
         write(io,101) ( jco(j), j=17, 32 )
         write(io,101) ( jco(j), j=33, 48 )
         write(io,101) ( jco(j), j=49, 64 )
         write(io,101) ( jco(j), j=65, 80 )
         write(io,101) ( jco(j), j=80, jms )
      else
         write(*,'(t1,a)')
     &   'Max crack front nodes for J-calc., currently equal to 96'
      endif
      write(io,'(t1,a)') '*endstep'
      close(io)
 101  format(t1,' ',a,15(' ',a))
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine constraint_equ_abaqus(ima,kma,fil6,uni,
     &                                 noda,na1,na2,na3)
C--- The routine creates constraint equations for ABAQUS
c
      implicit  none
c
      include 'common_nod.f'
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
      subroutine warp3d_inp(job,jobh,etyp,elnum_cell,elnum,no_of_nodes,
     1           estk_gc,estk_c,estk_s,nea2,estk_a, nodc,nc1,nc2,nc3,
     2           nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c-- The routine generates an input deck for a FEM-analysis using WARP3D
c
      implicit none
c
      include 'common_nod.f'
      include 'common_eln.f'
      integer  elo2n(iem),eln2o(iem)
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer  estk_gc(2,*),estk_c(2,*),estk_s(2,*),
     &         nea2,estk_a(2,nea2,*)
c
      integer  jobh,etyp,elnum_cell,elnum,no_of_nodes, ntip(401),njc,
     1         elno,io,i,j,jj,jn(3),je(3),id,ioj,
     2         iblock,ipart,blocksize
      integer  ndomj(5000), edomj(5000),nn,ne
      integer  ivek1(20),ivek2(20),num1,num2,iend
c
      double precision  youngs,poisson,sigy,hard
c
      character job*40,filinp*40,   el_list*10, answ*3, row*100,
     1          domain*10,  crdfile*40, file_n2o*30,
     2          elfile1*40,elfile2*40,constfile*50,prdspfile*50,
     3          blckfile*40, ablck*60, dname1*7,dname2*7,dname3*7
c
      logical ren_el,ok
c
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      double precision      dcell,eta_n,eta_t1,eta_t2
      common    /geom_cell/ dcell,eta_n,eta_t1,eta_t2
c
c  Check if the elements has been re-numbered
c
      write(*,'(t1,a,$)')
     &  '>> Have the elements been re-numbered using PATWARP ? (y/n): '
      read(*,'(a)') answ
      if ( (answ.eq.'Y').or.(answ.eq.'y') ) then
         ren_el = .true.
 7       continue
         write(*,'(t1,a,$)')
     &       '>> New to old element correspondence file name = '
         read(*,'(a)') file_n2o
         inquire(file=file_n2o,exist=ok)
         if (.not.ok) then
            write(*,'(t4,3a)')'>>The file: ',file_n2o,'doesn''t exist!'
            goto 7
         endif
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
      filinp=job(1:jobh)//'.warp_input'
      open(unit=14,file=filinp,status='unknown')
      io=14
c
c  write out a general discription of the problem as a comment
c
	write(io,'(4(t1,a/),t1,a)')
     & 'c ',
     & 'c   ***  CRACK GROWTH ANALYSIS  ***',
     & 'c ',
     & 'c  surface cracked plate:',
     & 'c  (one quarter of the plate is modelled)',
     & 'c '
c
      write(io,8001) 'structure surface_crack'
      write(io,8001) 'c'
c
      youngs = 30000.0
      poisson = 0.30
      sigy = 60.0
      hard = 10.0
c
      write(io,8001) 'material gurson_layer'
      write(io,8910) ' properties gurson  e ',youngs,' nu ',poisson,
     &                ' yld_pt ',sigy,' n_power ',hard,','
      write(io,8920) ' killable  f_0 0.001  q1 1.25 q2 1.0 ',
     &                'q3 1.5625  rho 0.0,'
      write(io,8001) ' nucleation  e_n 0.40  s_n 0.05  f_n 0.50'
c
 8001 format(t1,a)
 8002 format(t1,a,i7)
 8003 format(t1,5a)
 8004 format(t1,a,g15.8)
 8910 format(t1,a,f8.0,a,f5.3,a,f3.0,a,f3.0,a)
 8920 format(t1,a,a)
 9001 format(t1,a,f8.0,a,f5.3,a,f3.0,a/t20,a)
c
      
      write(io,*) 'c'
      write(io,8001) 'material steel'
      write(io,9001) ' properties mises e ',youngs,' nu ',poisson,
     &               '  yld_pt 60.0  n_power ',hard,' ,',  'rho 0.0'
c
      write(io,8001) 'c'
      write(io,8002) 'number of nodes ',no_of_nodes
      write(io,8002) 'number of elements',elnum
c
c Coordinates
c
      write(io,8001) 'c'
      write(io,8001) 'coordinates'
      crdfile = job(1:jobh)//'.crd'
      write(io,8003) '*input from ''',crdfile(1:(jobh+4)), ''''
      write(io,8001) 'c'
c
        open(unit=47,file=crdfile(1:(jobh+4)),status='unknown')
        do i=1, inm
           if (nnr(i).gt.0) write(47,201) nnr(i),
     &                      npos(i,1),npos(i,2),npos(i,3)
        enddo
        close(47)
 201    format(t3,i6, 3(g20.10) )
c
c Element incidences
c
      write(io,8001) 'elements'
      el_list = '          '
      if (elnum_cell.lt.100) then
         write(el_list,'(t1,a,i2)') '1-',elnum_cell
      elseif (elnum_cell.lt.1000) then
         write(el_list,'(t1,a,i3)') '1-',elnum_cell
      elseif (elnum_cell.lt.10000) then
         write(el_list,'(t1,a,i4)') '1-',elnum_cell
      endif
      write(io,'(t4,3a)') el_list,' type l3disop   nonlinear  ',
     &                                ' material gurson_layer   ,'
      write(io,'(t15,a)')
     &         'order 2x2x2  bbar center_output  short'
c
      write(io,8001) 'c'
      if (elnum.lt.100) then
         write(el_list,'(t1,i4,a,i2)') elnum_cell+1,'-',elnum
      elseif (elnum.lt.1000) then
         write(el_list,'(t1,i4,a,i3)') elnum_cell+1,'-',elnum
      elseif (elnum.lt.10000) then
         write(el_list,'(t1,i4,a,i4)') elnum_cell+1,'-',elnum
      elseif (elnum.lt.100000) then
         write(el_list,'(t1,i4,a,i5)') elnum_cell+1,'-',elnum
      endif
      write(io,'(t4,3a)') el_list,' type l3disop   nonlinear  ',
     &                                ' material steel ,'
      write(io,'(t15,a)')
     &         'order 2x2x2  bbar  center_output short'
      write(io,8001) 'c'
      write(io,8001) 'incidences'
      elfile1=job(1:jobh)//'_cl.elm'
      elfile2=job(1:jobh)//'_rg.elm'
      write(io,8003) '*input from ''',elfile1(1:(jobh+7)),''''
      write(io,8003) '*input from ''',elfile2(1:(jobh+7)),''''
      write(io,8001) 'c'
c
      open(unit=48,file=elfile1(1:(jobh+7)),status='unknown')
      open(unit=49,file=elfile2(1:(jobh+7)),status='unknown')
      if (ren_el) then
         do i=1, elnum_cell
            write(48,203) i, ( eln(eln2o(i),j), j=1, etyp)
         enddo
         do i=elnum_cell+1, elnum
            write(49,203) i, ( eln(eln2o(i),j), j=1, etyp)
         enddo
      else
         do i=1, elnum_cell
            write(48,203) ( eln(i,j), j=0, etyp)
         enddo
         do i=elnum_cell+1, elnum
            write(49,203) ( eln(i,j), j=0, etyp)
         enddo
      endif
      close(48)
      close(49)
 203  format( t1,9(i7) )
c
      if (ren_el) then
 40      continue
         write(*,'(t2,a,$)')
     &          '*Give input file with new element blocking : '
         read(*,'(a)') blckfile
         inquire(file=blckfile,exist=ok)
         if (.not.ok) then
            write(*,'(t4,3a)')'>>The file: ',blckfile,'doesn''t exist!'
            goto 40
         endif
         open(unit=46,file=blckfile,status='old')
 45      continue
         read(46,'(a)',end=50) ablck
         write(io,'(a)') ablck
         goto 45
 50      continue
         close(46)
      else
         write(io,8001) 'blocking'
          write(*,'(t1,a,$)')
     &     '>> WARP3D block size (only scalar blocking supported): '
           read(*,*) blocksize
         iblock = 0
         ipart = mod(elnum_cell,blocksize)
         do i = 1, elnum_cell-ipart, blocksize
            iblock = iblock + 1
            write(io,'(t5,i4,t15,i3,t25,i6)') iblock, blocksize, i
         enddo
         if (ipart.gt.0) then
            iblock = iblock + 1
            write(io,'(t5,i4,t15,i3,t25,i6)') iblock,ipart,
     &                elnum_cell-ipart+1
         endif
c
         elno = elnum-elnum_cell
         ipart = mod(elno,blocksize)
         do i=1, elno-ipart, blocksize
            iblock = iblock + 1
            write(io,'(t5,i4,t15,i3,t25,i6)')
     &           iblock,blocksize,i+elnum_cell
         enddo
         if (ipart.gt.0) then
            iblock = iblock + 1
            write(io,'(t5,i4,t15,i3,t25,i6)') iblock,ipart,
     &                 elno-ipart+1 + elnum_cell
         endif
c
      endif
c
      write(io,8001) 'c'
      write(io,8001) 'c     step cases:'
      write(io,8001) 'c'
      write(io,9010) 'c','displacement increments'
      write(io,8001) 'c'
      write(io,8001) 'c'
 9010 format(t1,a,t37,a)
 9011 format(t1,a,t3,a,t18,a,t35,a,t55,a)
 9012 format(t1,a,t5,i2,t15,i4,a,i4,t30,g15.8,t50,g15.8  )
c
      write(io,8001) 'c'
      write(io,8001) 'c  define a "fake" prescribed load step'
      write(io,'(t1,a,a)') 'c  (necessary even though the ',
     &                        'loading is prescribed displacements'
      write(io,8001) 'c'
      write(io,8001) 'loading fake '
      write(io,8001) '   nodal loads '
      write(io,8001) '   2 force_x 0.0'
      write(io,8001) 'c'
      write(io,8001) 'loading predisp'
      write(io,8001) '    nonlinear'
      write(io,'(t1,a,i)') '    steps  1-1000  fake  1.0'
c
      write(io,8001) 'c'
      write(io,8001) 'c    solution paramters'
      write(io,8001) 'c'
      write(io,8001) ' nonlinear analysis parameters'
      write(io,8001) 'c   solution technique direct sparse '
      write(io,8001) '   solution technique lnpcg'
      write(io,8001) 'c   solution technique direct sparse '
      write(io,8001) 'c   solution technique direct sparse hp'
      write(io,8001) 'c   solution technique direct sparse bcs'
      write(io,8001) 'c   solution technique direct sparse sgi'
      write(io,8001) '   preconditioner type diagonal'
      write(io,8001) '   lnr_pcg conv test res tol 0.01'
      write(io,8001) '   maximum linear iterations 20000'
      write(io,8001) '   maximum iterations 5'
      write(io,8001) '   minimum iterations 2'
      write(io,8001) '   convergence test norm res tol 0.01'
      write(io,8001) '   nonconvergent solutions stop'
      write(io,8001) '   adaptive on'
      write(io,8001) '   linear stiffness for iteration one off'
      write(io,8001) '   batch messages off'
      write(io,8001) '   cpu time limit off'
      write(io,8001) '   material messages off'
      write(io,8001) '   bbar stabilization factor 0.0'
      write(io,8001) '   consistent q-matrix on'
      write(io,8001) '   extrapolation on'
      write(io,8001) '   time step 1.0e06'
      write(io,8001) 'c'
      write(io,8001) 'crack growth parameters'
      write(io,8001) '   type of crack growth element_extinction'
      write(io,8001) '   critical porosity 0.2'
      write(io,8001) '   force release type traction-separation'
      write(io,8004) '   cell height',dcell*0.5
      write(io,8001) '   release fraction 0.10'
      write(io,8001) '   crack plane normal y coordinate 0.0'
      if (elnum_cell.lt.100) then
         write(io,'(t1,a,i2)') '   print status on order 1-',elnum_cell
      elseif (elnum_cell.lt.1000) then
         write(io,'(t1,a,i3)') '   print status on order 1-',elnum_cell
      elseif (elnum_cell.lt.10000) then
         write(io,'(t1,a,i4)') '   print status on order 1-',elnum_cell
      elseif (elnum_cell.lt.100000) then
         write(io,'(t1,a,i5)') '   print status on order 1-',elnum_cell
      endif
c
      write(io,8001) 'c'
      write(io,8001) 'c    start the analysis'
c
      write(io,8001) 'c'
      write(io,8001) 'constraints'
      write(io,8001) '*echo off'
      constfile = job(1:jobh)//'_001.const'
      write(io,8003) '*input from ''',constfile(1:(jobh+10)),''''
      prdspfile = job(1:jobh)//'_001.prdsp'
      write(io,8003) '*input from ''',prdspfile(1:(jobh+10)),''''
      prdspfile = job(1:jobh)//'_002.prdsp'
      write(io,8003) 'c  *input from ''',prdspfile(1:(jobh+10)),''''
      write(io,8001) '*echo on'
      write(io,8001)
     &    'compute displacements for load predisp steps 1'
c     &   'compute displacements for load predisp steps 1-10'
      write(io,8001) 'output patran binary displacement'
      write(io,8001) 'output patran binary stresses'
      write(io,8001) 'output patran binary strains'
      write(io,8001) 'c'
      write(io,8001) 'c     j-integral calculations'
      write(io,8001) 'c'
c
c  Currently only one domain is choosen which is located inside Zone S
c
      call nodset_ntip_(ntip,njc, nodc,nc1,nc2,nc3)
c
      if (etyp.ne.8) write(*,'(/t1,a,a/)') '*WARNING! if not ',
     &   '8-noded elements, the J-domains will be wrong!'
c
c . . find nodes where q = 1 and the corresponding elements.
      id = 0
      do j=1, jmc, 2
         jj = (j+1)/2
         call j_domain_lists(j,nn,ne,ndomj,edomj,jn,je,etyp, estk_gc,
     &        estk_c,estk_s, nc1,nc2,nc3,nodc, ns1,ns2,ns3,nods)
         if (jj.lt.10) then
            write(domain,'(t1,a,i1)') 'domain_00',jj
         elseif (jj.lt.100) then
            write(domain,'(t1,a,i2)') 'domain_0',jj
         else
            write(domain,'(t1,a,i3)') 'domain_',jj
         endif
c
         id = id + 1
         if (id.lt.10) then
            write(dname1,'(t1,a,i1)')  'dnr_00',id
         elseif (id.lt.100) then
            write(dname1,'(t1,a,i2)')  'dnr_0',id
         else
            write(dname1,'(t1,a,i3)')  'dnr_',id
         endif
c
         id = id + 1
         if (id.lt.10) then
            write(dname2,'(t1,a,i1)')  'dnr_00',id
         elseif (id.lt.100) then
            write(dname2,'(t1,a,i2)')  'dnr_0',id
         else
            write(dname2,'(t1,a,i3)')  'dnr_',id
         endif
c
         id = id + 1
         if (id.lt.10) then
            write(dname3,'(t1,a,i1)')  'dnr_00',id
         elseif (id.lt.100) then
            write(dname3,'(t1,a,i2)')  'dnr_0',id
         else
            write(dname3,'(t1,a,i3)')  'dnr_',id
         endif
c
         write(io,8003) '*input from file ''',domain ,''''
c
         ioj=49
         open(unit=ioj,file=domain,status='unknown')
         if (jj.eq.1) then
            call warp3d_q2_domain(dname1,ioj,ntip(1),ntip(2),0,'a',
     &          jn(1),1,je(1),ndomj,edomj,ren_el,elo2n)
            call warp3d_q2_domain(dname2,ioj,ntip(1),ntip(2),0,'a',
     &          jn(2),je(1)+1,je(2),ndomj,edomj,ren_el,elo2n)
            call warp3d_q2_domain(dname3,ioj,ntip(1),ntip(2),0,'a',
     &          jn(3),je(2)+1,je(3),ndomj,edomj,ren_el,elo2n)
         elseif (jj.eq.njc) then
            call warp3d_q2_domain(dname1,ioj,0,ntip(njc-1),ntip(njc),
     &           'c',jn(1),1,je(1),ndomj,edomj,ren_el,elo2n)
            call warp3d_q2_domain(dname2,ioj,0,ntip(njc-1),ntip(njc),
     &           'c',jn(2),je(1)+1,je(2),ndomj,edomj,ren_el,elo2n)
            call warp3d_q2_domain(dname3,ioj,0,ntip(njc-1),ntip(njc),
     &           'c',jn(3),je(2)+1,je(3),ndomj,edomj,ren_el,elo2n)
         else
            call warp3d_q2_domain(dname1,ioj,ntip(jj-1),ntip(jj),
     &       ntip(jj+1),'b',jn(1),1,je(1),ndomj,edomj,ren_el,elo2n)
            call warp3d_q2_domain(dname2,ioj,ntip(jj-1),ntip(jj),
     &      ntip(jj+1),'b',jn(2),je(1)+1,je(2),ndomj,edomj,ren_el,elo2n)
            call warp3d_q2_domain(dname3,ioj,ntip(jj-1),ntip(jj),
     &      ntip(jj+1),'b',jn(3),je(2)+1,je(3),ndomj,edomj,ren_el,elo2n)
         endif
         close(ioj)
      enddo
      write(io,8001) 'c'
      write(io,8001) 'stop'
c
      close(io)
c
c  Generate the constraint file
c
      call warp3d_const(nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3, 
     &                  noda,na1,na2,na3, nodb,nb1,nb2,nb3, job,jobh)
c
c  Generate the files containing the domain integral computations
c
      call warp3d_prdsp(nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3, 
     &                  noda,na1,na2,na3, nodb,nb1,nb2,nb3, job,jobh)
c
      call write_out_tipnodes(nodc,nc1,nc2,nc3,job,jobh)
c
      call write_out_crackplane(nodc,nc1,nc2,nc3,job,jobh,
     &               estk_gc,ren_el,elo2n)
c
      call  write_out_stress_layers(nodc,nc1,nc2,nc3,
     &      nods,ns1,ns2,ns3, noda,na1,na2,na3, job,jobh,
     &      estk_gc,estk_c,estk_s,nea2,estk_a,ren_el,elo2n)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine write_out_stress_layers(nodc,nc1,nc2,nc3,
     &      nods,ns1,ns2,ns3, noda,na1,na2,na3, job,jobh,
     &      estk_gc,estk_c,estk_s,nea2,estk_a,ren_el,elo2n)
c
      implicit none
c
      include 'common_nod.f'
      include 'common_eln.f'
      integer  elo2n(iem)
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3)
c
      integer  estk_gc(2,*),estk_c(2,*),estk_s(2,*),
     &         nea2,estk_a(2,nea2,*),jobh
c
      integer  io,elay(2000,50),ie(100),i,j,k,i1,i2,j1,j2,nl,nl1,
     &        jt0,jt
c
      character job*40,layfile*40
c
      logical ren_el
c
      integer       mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/  mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
      layfile = job(1:jobh)//'_lay.inf'
c
c  find the elements on the different layers perpendicular to the
c  crack front
c
      nl = 0
      do j=1, jmc-2, 2*sjred_type
         nl = nl + 1
         ie(nl) = 0
c
         j1 = (j+1)/2
         j2 = (j+1)/2 + (sjred_type - 1)
c
c Zone C (Gurson elements):
         i1 = estk_gc(1,j1)
         i2 = estk_gc(2,j2)
         do i=i1, i2
            ie(nl) = ie(nl) + 1
            if (ren_el) then
               elay( ie(nl), nl) = eln(elo2n(i),0)
            else
               elay( ie(nl), nl) = eln(i,0)
            endif
         enddo
c
c Zone C:
         i1 = estk_c(1,j1)
         i2 = estk_c(2,j2)
         do i=i1, i2
            ie(nl) = ie(nl) + 1
            if (ren_el) then
               elay( ie(nl), nl) = eln(elo2n(i),0)
            else
               elay( ie(nl), nl) = eln(i,0)
            endif
         enddo
c
c Zone S:
         i1 = estk_s(1,j1)
         i2 = estk_s(2,j2)
         do i=i1, i2
            ie(nl) = ie(nl) + 1
            if (ren_el) then
               elay( ie(nl), nl) = eln(elo2n(i),0)
            else
               elay( ie(nl), nl) = eln(i,0)
            endif
         enddo
c
      enddo
      nl1 = nl
c
c Zone A:
      nl = 0
      if (mod(na,4).ne.0) then
         jt0=0
      else
         jt0=2
      endif
      jt = jt0
      do j=1, jma-2, 2
         nl = nl + 1
         j1 = (j+1)/2
         j2 = (j+1)/2
         do k=1, kar1-2, 2
            if (nl.ge.(na-1)) then
               i1 = estk_a( 1, (k+1)/2, j1)
            elseif (mod(jt,4).eq.0) then
               i1 = estk_a( 1, (k+1)/2, j1) + 2
            elseif (mod(jt,4).eq.1) then
               i1 = estk_a( 1, (k+1)/2, j1) + 1
            elseif (mod(jt,4).eq.2) then
               i1 = estk_a( 1, (k+1)/2, j1) + 1
            elseif (mod(jt,4).eq.3) then
               i1 = estk_a( 1, (k+1)/2, j1) + 2
            endif
            i2 = estk_a( 2, (k+1)/2, j2)
            do i=i1, i2
               ie(nl) = ie(nl) + 1
               if (ren_el) then
                  elay( ie(nl), nl) = eln(elo2n(i),0)
               else
                  elay( ie(nl), nl) = eln(i,0)
               endif
            enddo
         enddo
         jt = jt + 1
c
      enddo
c
      if (nl1.ne.nl) then
         write(*,'(t1,a,a)') '>> WARNING  nl1 ne nl :in subroutine ',
     &         'write_out_stress_layers CHECK THE CODE!'
      endif
c
      layfile = job(1:jobh)//'_lay.inf'
      io = 19
      open(unit=io,file=layfile,status='unknown')
      write(io,'(a,a)')'*Element layers along the crack front ',
     &   '<no-elem-layers>  /    <layer-no> <no-elem-per-layer>'
      write(io,'(t1,2i7)')  nl
      do j=1, nl
         write(io,'(t3,10i7)')   j, ie(j)
         write(io,'(t3,10i7)') ( elay(i,j), i=1, ie(j) )
      enddo
      close(io)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine write_out_crackplane(nodc,nc1,nc2,nc3,job,jobh,
     &           estk_gc,ren_el,elo2n)
c
      implicit none
      include 'common_nod.f'
      include 'common_eln.f'
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),jobh,estk_gc(2,*),elo2n(*)
      integer  i,j,k,ii,jj,i1,i2,nni,nnj,nei,nej,ncrk(100,200),
     &         ecrk(100,200),io
      character job*40,crkfile*40
      logical ren_el
c
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
      crkfile = job(1:jobh)//'_crk.inf'
c
c     find nodes on the crack plane
c
      k = 1
      jj = 0
      do j=1, jmc, 2
         jj = jj + 1
         ii = 0
         do i=1, imc-2, 2
            if (nnr(nodc(i,j,k)).gt.0) then
               ii = ii + 1
               ncrk(ii,jj) = nnr(nodc(i,j,k))
            endif
         enddo
      enddo
      nni = ii
      nnj = jj
c
c     find the elements on the crack plane
c
      do j=1, jmc-2, 2
         jj = (j+1)/2
         i1 = estk_gc(1,jj)
         i2 = estk_gc(2,jj)
         ii = 0
         do i=i1, i2
            ii = ii + 1
            if (ren_el) then
               ecrk(ii,jj) = eln(elo2n(i),0)
            else
               ecrk(ii,jj) = eln(i,0)
            endif
         enddo
      enddo
      nei = ii
      nej = jj
c
      crkfile = job(1:jobh)//'_crk.inf'
      io = 19
      open(unit=io,file=crkfile,status='unknown')
      write(io,'(a,a)')'*crack plane nodes    ',
     &      '<no-nodes-phi-dir>     <no-nodes-growth-dir>'
      write(io,'(t1,2i7)')  nnj, nni
      do j=1, nnj
         write(io,'(t3,10i7)') (ncrk(i,j),i=1, nni)
      enddo
      write(io,'(a,a)')'*crack plane elemnts    ',
     &      '<no-elem-phi-dir>     <no-elem-growth-dir>'
      write(io,'(t1,2i7)')  nej, nei
      do j=1, nej
         write(io,'(t3,10i7)') (ecrk(i,j),i=1, nei)
      enddo
      close(io)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine write_out_tipnodes(nodc,nc1,nc2,nc3,job,jobh)
c
      implicit none
      include 'common_nod.f'
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3)
      integer  jobh,ntip(401,2),i,j,k,n,io
      double precision  x,y,z,phi
      character job*40,tipfile*40
c
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      double precision  t,w,c,a,alfa
      common  /geom/    t,w,c,a,alfa
c
      n = 0
c-Zone C:
      i = 5
      k = 1
      do j=1, jmc
         if (nnr(nodc(i,j,k)).gt.0) then
            n=n+1
            ntip(n,1) = nnr(nodc(i,j,k))
            ntip(n,2) = nodc(i,j,k)
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
 111  format(t1,i7,4g16.8)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine warp3d_const(nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3,
     &      noda,na1,na2,na3, nodb,nb1,nb2,nb3, job,jobh)
c
c Routine determines the fixed BC.
c
      implicit none
c
      include 'common_nod.f'
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
c
      integer    nstm
      parameter (nstm=20000)
      integer    nfront(nstm),nsym(nstm),nconst(0:3,nstm)
c
      integer jobh, i,j,cuni, n,n1,n2, node, iu,iv,iw
c
      double precision u,v,w
c
      character job*40,constfile*50
      logical   exist
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer      imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
c Node Set = NFRONT (Y=0)
c
      call nodset_nfront_(nfront,n1, nodc,nc1,nc2,nc3,
     &     nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
c
c Node Set = NSYM (X=0)
c
      call nodset_nsym_(nsym,n2, nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3, 
     &                           noda,na1,na2,na3)
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
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine warp3d_prdsp(nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3, 
     &      noda,na1,na2,na3, nodb,nb1,nb2,nb3, job,jobh)
c
c Routine determines the fixed BC.
c
      implicit none
c
      include 'common_nod.f'
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
      integer    nstm
      parameter (nstm=20000)
      integer    nbak(nstm),nsym(nstm),nline(1000)
c
      integer jobh, i,j,k,jj,jend,cuni, n,n2,nl1,nl2, node,iu,iv,iw
c
      double precision u,v,w, coord(3,nstm),zmin,zmax,zm,h,delta,theta
c
      character job*40,prdspfile*50
C
      logical   unique
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer      imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
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
c Node Set = NSYM (X=0) (Avoid specifying BC. twice)
c
      call nodset_nsym_(nsym,n2, nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3,
     &                           noda,na1,na2,na3)
c
c Find nodes for a four point bend loading case
c
      k = kma
      i = 1
      jj = 0
      do j=1, jma-4
         if (nnr(noda(i,j,k)).gt.0) then
            jj = jj + 1
            nline(jj) = nnr(noda(i,j,k))
         endif
      enddo
      j = jma
      do i=6, ima
         if (nnr(noda(i,j,k)).gt.0) then
            jj = jj + 1
            nline(jj) = nnr(noda(i,j,k))
         endif
      enddo
      j = jmb
      do i=2, imb
         if (nnr(nodb(i,j,k)).gt.0) then
            jj = jj + 1
            nline(jj) = nnr(nodb(i,j,k))
         endif
      enddo
      nl1 = jj
c
      nl2 = 0
      k = kma - 4
      i = ima
      jend = jma - jmb + 1
      do j=1, jend
         if (nnr(noda(i,j,k)).gt.0) then
            jj = jj + 1
            nline(jj) = nnr(noda(i,j,k))
         endif
      enddo
      j = 1
      do i=2, imb
         if (nnr(nodb(i,j,k)).gt.0) then
            jj = jj + 1
            nline(jj) = nnr(nodb(i,j,k))
         endif
      enddo
      nl2 = jj
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
c  Four point bending
c
      prdspfile = job(1:jobh)//'_003.prdsp'
      cuni = 55
      open(unit=cuni,file=prdspfile(1:jobh+10),status='unknown')
      do i=1, nl2
         node = nline(i)
         if (i.le.nl1) then
            iu = 0
            iv = 0
            iw = 1
            u  = 0.0
            v  = 0.0
            w  = 1.0
         else
            iu = 0
            iv = 0
            iw = 1
            u  = 0.0
            v  = 0.0
            w  = 0.0
         endif
         call fixed_constraint_warp3d(cuni,node,iu,iv,iw,u,v,w)
      enddo
      close(cuni)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine fixed_constraint_warp3d(cuni,node,iu,iv,iw,u,v,w)
c
c  Write out constraint on WARP3D-format on unit=CUNI
c
      implicit none
      integer cuni,node,iu,iv,iw
      double precision u,v,w
c
      if     ( (iu.eq.1) .and. (iv.eq.0) .and. (iw.eq.0) ) then
         write(cuni,101) node,'  u ',u
      elseif ( (iu.eq.1) .and. (iv.eq.1) .and. (iw.eq.0) ) then
         write(cuni,102) node,'  u ',u,' v ',v
      elseif ( (iu.eq.1) .and. (iv.eq.0) .and. (iw.eq.1) ) then
         write(cuni,102) node,'  u ',u,' w ',w
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
c
 101  format(i6,a,g15.8)
 102  format(i6,2(a,g15.8))
 103  format(i6,3(a,g15.8))
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine j_domain_lists(j,nn,ne,ndomj,edomj,jn,je,etyp,
     &    estk_gc,estk_c,estk_s, nc1,nc2,nc3,nodc, ns1,ns2,ns3,nods)
c
      implicit none
c
      include 'common_eln.f'
      include 'common_nod.f'
c
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     &         ns1,ns2,ns3,nods(ns1,ns2,ns3)
c
      integer  ndomj(*),edomj(*),etyp,etmp(4000),iok(4000)
      integer  estk_gc(2,*),estk_c(2,*),estk_s(2,*)
c
      integer  nn,ne,i,j,k,j1,j2,jj,nod,eno,jn(3),je(3),
     &         in,ie,i1gc,i2gc,i1c,i2c,i1s,i2s,iemax
c
c      integer  it
c      character tecfile*11
c
c      logical  ok1,unique
c
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
c
c  TEST open tec files
c
c      if (j.le.9) then
c         write(tecfile,'(t1,a,i1,a)') 'jtec_0',j,'.plt'
c      else
c         write(tecfile,'(t1,a,i2,a)') 'jtec_',j,'.plt'
c      endif
c      it = 88
c      open(unit=it,file=tecfile,status='unknown')
c      write(it,'(t1,a)') 'ZONE'
c
c  Find the nodes where q = 1
c
      nn = 0
      do k=1, kmc
         do i=1, imc
            nod = nnr(nodc(i,j,k))
            if (nod.gt.0) then
               nn = nn + 1
               ndomj(nn) = nod
c             write(it,*) real( npos(nodc(i,j,k),1)),
c     &             real(npos(nodc(i,j,k),2)), real(npos(nodc(i,j,k),3))
            endif
         enddo
      enddo
c
      do k=2, kms-4
         if (k.eq.(kms-6)) then
            jn(1) = nn
         elseif (k.eq.(kms-4)) then
            jn(2) = nn
         endif
         do i=1, ims
            nod = nnr(nods(i,j,k))
            if (nod.gt.0) then
               nn = nn + 1
               ndomj(nn) = nod
c             write(it,*) real(npos(nods(i,j,k),1)),
c     &             real(npos(nods(i,j,k),2)), real(npos(nods(i,j,k),3))
            endif
         enddo
       enddo
       jn(3) = nn
c
c      close(it)
c
c  Find the elements containing the nodes where q = 1
c
      j2 = (j+1)/2
      j1 = j2 - 1
      jj = ((jmc+1)/2)
      if (j1.lt.1) then
         i1gc = estk_gc(1,1)
         i2gc = estk_gc(2,1)
         i1c = estk_c(1,1)
         i2c = estk_c(2,1)
         i1s = estk_s(1,1)
         i2s = estk_s(2,1)
      elseif (j2.ge.jj) then
         i1gc = estk_gc(1,j1)
         i2gc = estk_gc(2,j1)
         i1c = estk_c(1,j1)
         i2c = estk_c(2,j1)
         i1s = estk_s(1,j1)
         i2s = estk_s(2,j1)
      else
         i1gc = estk_gc(1,j1)
         i2gc = estk_gc(2,j2)
         i1c = estk_c(1,j1)
         i2c = estk_c(2,j2)
         i1s = estk_s(1,j1)
         i2s = estk_s(2,j2)
      endif
      k = 0
      do i=i1gc, i2gc
         k = k + 1
         etmp(k) = i
      enddo
      do i=i1c, i2c
         k = k + 1
         etmp(k) = i
      enddo
      ne = k
      do i=i1s, i2s
         k = k + 1
         etmp(k) = i
      enddo
      iemax = k
c
c  Find the elements in zone S that belongs to the j-domains
c
      do i=1, iemax
         iok(i) = 0
      enddo
c
c  Elements in ring 1:
c
      do i=1, jn(1)
         nod = ndomj(i)
         do ie=ne+1, iemax
            eno = etmp(ie)
            do in=1, etyp
               if (eln(eno,in).eq.nod) iok(ie) = 1
            enddo
         enddo
      enddo
c
      do i=1, ne
         edomj(i) = etmp(i)
      enddo
      je(1) = ne
c
      do ie=ne+1, iemax
         if ( (iok(ie).gt.0) .and. (iok(ie).le.1) ) then
            je(1) = je(1) + 1
            edomj(je(1)) = etmp(ie)
         endif
      enddo
c
c  Elements in ring 2:
c
      do i=1, jn(2)
         nod = ndomj(i)
         do ie=ne+1, iemax
            eno = etmp(ie)
            do in=1, etyp
               if (eln(eno,in).eq.nod) iok(ie) = 2
            enddo
         enddo
      enddo
c
      je(2) = je(1)
      do i=1, ne
         je(2) = je(2) + 1
         edomj(je(2)) = etmp(i)
      enddo
c
      do ie=ne+1, iemax
         if ( (iok(ie).gt.0) .and. (iok(ie).le.2) ) then
            je(2) = je(2) + 1
            edomj(je(2)) = etmp(ie)
         endif
      enddo
c
c  Elements in ring 3:
c
      do i=1, jn(3)
         nod = ndomj(i)
         do ie=ne+1, iemax
            eno = etmp(ie)
            do in=1, etyp
               if (eln(eno,in).eq.nod) iok(ie) = 3
            enddo
         enddo
      enddo
c
      je(3) = je(2)
      do i=1, ne
         je(3) = je(3) + 1
         edomj(je(3)) = etmp(i)
      enddo
c
      do ie=ne+1, iemax
         if ( (iok(ie).gt.0) .and. (iok(ie).le.3) ) then
            je(3) = je(3) + 1
            edomj(je(3)) = etmp(ie)
         endif
      enddo
c
      return
      end
c
c      return
c      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine warp3d_q2_domain(dname,io,n1,n2,n3,dtype,
     &                     nn,ne1,ne2,ndomj,edomj,ren_el,elo2n)
c
c  Generates a file with domain integral commands
c
      implicit none
c
      include 'common_eln.f'
      integer  elo2n(iem)
c
      integer  n1,n2,n3,ndomj(*),edomj(*),nn,ne,ne1,ne2,
     &         io,i,ipart,etmp(3000)
      character dname*7,dtype*1, form*14
c      character domain*10
c
      logical ren_el
c
c      io=49
c      open(unit=io,file=domain,status='unknown')
c
      write(io,'(t1,a)') 'c'
      write(io,'(t1,a,a)') 'domain ',dname
      write(io,'(t1,a)')   '*echo off'
      write(io,'(t2,a)')   'symmetric'
      write(io,'(t2,a)')   'normal plane nx 0  ny -1.0  nz 0.0'
c
      if (dtype.eq.'a') then
         write(io,'(t2,a,2i7,a)') 'front nodes',n1,n2,' linear verify'
         write(io,'(t2,a)')       'function type a'
      elseif (dtype.eq.'b') then
         write(io,'(t2,a,3i7,a)')'front nodes',n1,n2,n3,' linear verify'
         write(io,'(t2,a)')       'function type b'
      elseif (dtype.eq.'c') then
         write(io,'(t2,a,2i7,a)') 'front nodes',n2,n3,' linear verify'
         write(io,'(t2,a)')       'function type c'
      endif
      write(io,'(t2,a)')          'ignore crack face loading'
c
      ipart = mod(nn,8)
      if (mod(nn,8).eq.0) ipart = 8
      write(io,'(t2,a,t12,8i7,'','')') 'q-values', (ndomj(i),i=1,8)
      write(io,'(t12,8i7,'','')')       (ndomj(i),i=9, nn-ipart)
      write(form,'(a,i1,a)') '(t12,',ipart,'i7,f6.3)'
      write(io,form)               (ndomj(i),i=nn-ipart+1, nn), 1.00
c
c  Element renumbered?
c
      ne =0
      if (ren_el) then
         do i=ne1, ne2
            ne = ne + 1
            etmp(ne) = elo2n(edomj(i))
         enddo
      else
         do i=ne1, ne2
            ne = ne + 1
            etmp(ne) = edomj(i)
         enddo
      endif
c
      ipart = mod(ne,10)
      if (ipart.eq.0) ipart = 10
      write(io,'(t2,a,t12,10i6,'','')') 'elements', (etmp(i),i=1,10)
      write(io,'(t12,10i6,'','')')      (etmp(i),i=11, ne-ipart)
      write(io,'(t12,10i6)')            (etmp(i),i=ne-ipart+1, ne)
c
      write(io,'(t1,a)')   '*echo on'
      write(io,'(t2,a)')     'print totals'
      write(io,'(t2,a)')     'compute domain integral'
c
c      close(io)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine patran_neu(job,jobh,etyp,elnum_cell,elnum,
     1           no_of_nodes, nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3, 
     2           noda,na1,na2,na3, nodb,nb1,nb2,nb3)
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
      integer  jobh,etyp,elnum,elnum_cell,no_of_nodes, uni,i,j
      integer  node,ndof,cflag(3),el_no,packet,inf_lines
      integer  igtype,tyel,config,pid,cid,ceid, head(15),versn(3)
c
      integer    nnc
      parameter (nnc=20000)
      integer    nconst(0:3,nnc), incst
      double precision uconst(3,nnc)
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
c      uni = 51
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
      call patran_const(nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3, 
     &     noda,na1,na2,na3, nodb,nb1,nb2,nb3, nnc,nconst,uconst,incst)
c
c      uni = 52
c
      packet = 8
      inf_lines = 2
      cid = 0
      do i=1, incst
         node     = nconst(0,i)
         cflag(1) = nconst(1,i)
         cflag(2) = nconst(2,i)
         cflag(3) = nconst(3,i)
         write(uni,9081) packet, node, 1, inf_lines, 0,0,0,0,0
         write(uni,9082) cid, (cflag(j),j=1, 3),  0, 0, 0
         write(uni,9083) uconst(1,i), uconst(2,i), uconst(3,i)
      enddo
 9081 format( i2, 8i8 )
 9082 format( i8, 6i1 )
 9083 format( 3e16.9 )
c
c  PACKET  2   =>  Element incidences
c
c   Gurson elements   config = 1
c   J2     elements   config = 2
c   etyp=nodes per elemnet;
c    (  pid - element property set (we don't use)  )
c    (  ceid - element coordinate system (we don't use)  )
c    (  <angles> - we don't use....  )
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
      pid  = 0
      ceid = 0
c
      config = 1
      do i=1, elnum_cell
         el_no = eln(i,0)
         write(uni,9021) packet, el_no, tyel, inf_lines, 0,0,0,0,0
         write(uni,9022) etyp, config, pid, ceid, 0.0, 0.0, 0.0
         write(uni,9023) ( eln(i,j), j=1, etyp )
      enddo
c
      config = 2
      do i=elnum_cell+1, elnum
         el_no = eln(i,0)
         write(uni,9021) packet, el_no, tyel, inf_lines, 0,0,0,0,0
         write(uni,9022) etyp, config, pid, ceid, 0.0, 0.0, 0.0
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
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine patran_const(nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3, 
     &   noda,na1,na2,na3, nodb,nb1,nb2,nb3, nnc,nconst,uconst,incst)
c
      implicit none
      include 'common_nod.f'
      integer  nc1,nc2,nc3,nodc(nc1,nc2,nc3),
     1         ns1,ns2,ns3,nods(ns1,ns2,ns3),
     2         na1,na2,na3,noda(na1,na2,na3),
     3         nb1,nb2,nb3,nodb(nb1,nb2,nb3)
      integer  nnc,nconst(0:3,*),incst 
      double precision uconst(3,*)
c
      integer    nstm
      parameter (nstm=20000)
      integer    nfront(nstm),nsym(nstm)
c
      integer   i,j, n,n1,n2
      logical   exist
c
c Node Set = NFRONT (Y=0)
c
      call nodset_nfront_(nfront,n1, nodc,nc1,nc2,nc3,
     &             nods,ns1,ns2,ns3, noda,na1,na2,na3, nodb,nb1,nb2,nb3)
      if (n1.gt.nstm) then
         write(*,'(/t1,a/t1,a,i6/t1,a,i6)')
     1    '>>Increase the size of nstm in subroutine patran_const(..)',
     2    '  Current size nstm = ',nstm,'  Needed size > ',n1
      endif
c
c Node Set = NSYM (X=0)
c
      call nodset_nsym_(nsym,n2, nodc,nc1,nc2,nc3, nods,ns1,ns2,ns3,
     &                           noda,na1,na2,na3)
      if (n2.gt.nstm) then
         write(*,'(/t1,a/t1,a,i6/t1,a,i6)')
     1    '>>Increase the size of nstm in subroutine patran_const(..)',
     2    '  Current size nstm = ',nstm,'  Needed size > ',n2
      elseif ((n1+n2).gt.nnc) then
         write(*,'(/t1,a/t1,a,i6/t1,a,i6)')
     1    '>>Increase the size of nnc in subroutine patran_neu(..)',
     2    '  Current size nnc = ',nnc,'  Needed size > ',(n1+n2)
      endif
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
      incst = n
      do i=1, incst
         uconst(1,i) = 0.0
         uconst(2,i) = 0.0
         uconst(3,i) = 0.0
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine mesh_statestik(job,jobh,datum,y,z,etyp,mbtype,mb_bias,
     &  ncell,elnum,elnum_cell,no_of_nodes,prog,  nea2,neb2, estk_gc,
     &   estk_c,estk_s,estk_a,estk_b, nstk_c,nstk_s,nstk_a,nstk_b)
C--------------------------------------------------------------*
C       Rutinen skapar en statestikfil innehallande            *
C       information om modellen.                               *
C--------------------------------------------------------------*
	implicit none
c
      include 'common_eln.f'
c
      include 'common_nod.f'
c
      integer  ncell,nea2,neb2, estk_gc(2,*),estk_c(2,*), estk_s(2,*),
     1         estk_a(2,nea2,*),estk_b(2,neb2,*),
     2         nstk_c(3,*),nstk_s(3,*),nstk_a(3,*),nstk_b(3,*)
c
      integer  jobh,etyp,elnum,elnum_cell,no_of_nodes,mbtype,
     &         uni,i,j,k,nnc,nns,nna,nnb,nez
c
      double precision  y(0:200),z(0:200,2),mb_bias
c 
      character job*40,prog*20,fil4*40,datum*20
c
      integer      mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      common /nod/ mr,mv,mh,m1,m2,ma,na,mb,nb,lt,lred
      integer       imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      common  /maxa/ imc,kmc,jmc,ims,jms,kms,ima,jma,kma,imb,jmb
      integer         ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
      common /reduce/ ksr1,kar1,kar2,rtype,sfred,sfred_type,sjred_type
c
      double precision  t,w,c,a,alfa
      common  /geom/    t,w,c,a,alfa
      double precision      dcell,eta_n,eta_t1,eta_t2
      common    /geom_cell/ dcell,eta_n,eta_t1,eta_t2
c
      fil4=job(1:jobh)//'_mesh.sta'
      open(unit=21,file=fil4,status='unknown')
      uni=21
c
      write(uni,'(/t10,a//t10,a)')
     &  '* * * *  m o d e l - i n f o r m a t i o n   * * * *',
     &  'the model is generated with the program mesh3d_scpcell'
      write(uni,'(/t5,a,a)') 'program-system : ',prog
      write(uni,'(/t5,a,a)') 'job       : ',job(1:jobh)
      write(uni,'(t5,a,a)')  'date      : ',datum
c
      write(uni,'(t5,a,f10.5,3(a,f10.5))')
     &  'geometri  : t=' ,t, '  w=' ,w, '  c=' ,c, '  a=' ,a
      write(uni,'(t5,a/t10,a,i3/t10,4(a,g14.6))')
     1  'zone c; Containing Computational Cells',
     2  'No. of cells ahead of initial crack front =',ncell,
     3  'D = ',dcell,'  eta_n=',eta_n,'  eta_t1=',eta_t1,
     4  '  eta_t2=',eta_t2
c
      write(uni,'(t5,6(a,i2))') 'zone s    : mr=',mr,
     1    ' mv=',mv,' mh=',mh,' sfred=',sfred,
     2    ' sfred_type=',sfred_type,' sjred_type=',sjred_type
c
      write(uni,'(t5,3(a,i2))')  'zone a    : m1=',m1,
     &   '  m2=',m2,'  na=',na
      write(uni,'(t13,3(a,i2))')
     &          '    rtype=',rtype,'  lt=',lt,'  lred=',lred
c
      write(uni,'(t5,2(a,i2),a,i2,a,f6.3)') 'zone b    : mb=',mb,
     &    '  nb=',nb,'  mbtype=',mbtype,'  mb_bias=',mb_bias
c
      write(uni,'(t5,a,i2,a//)') 'elementype: ',etyp,' -noded'
c
      write(uni,'(t5,a,i6)') 'Total number of nodes     = ',no_of_nodes
      write(uni,'(t5,a,i6)') 'Total number of elements  = ',elnum
      write(uni,'(t5,a,i6)') 'Number of Gurson elements = ',elnum_cell
c
      nnc = 0
      do j=1, jmc
         nnc = nnc + nstk_c(3,j)
      enddo
      nns = 0
      doj=1, jms
         nns = nns + nstk_s(3,j)
      enddo
      nnb = 0
      do k=1, kma
         nnb = nnb + nstk_b(3,k)
      enddo
      nna = no_of_nodes - nnc - nns - nnb
c
      write(uni,'(//t5,a/)')
     & '-------- n o d - i n f o r m a t i o n -----------'
c
      write(uni,1001) 'zone', 'range', 'tot. no. of nodes'
      write(uni,1002) 'c', nstk_c(1,1), '-', nstk_c(2,jmc), nnc
      write(uni,1002) 's', nstk_s(1,1), '-', nstk_s(2,jms), nns
      write(uni,1002) 'a', nstk_a(1,1), '-', nstk_a(2,kma), nna
      write(uni,1002) 'b', nstk_b(1,1), '-', nstk_b(2,kma), nnb
1001  format(t4,a,t14,a,t25,a)
1002  format(t6,a,t10,i6,a,i6,t25,i6)
c
      write(uni,'(t4,a,i6,a,i6)') 'Node interval on the rear surface: ',
     &           nstk_a(1,kma), ' - ', nstk_b(2,kma)
c
      write(uni,'(/t5,a/)')
     &  '-------- e l e m e n t - i n f o r m a t i o n -----------'
c
      write(uni,1011) 'zone', 'No. of elements in zone'
      do j=jmc, 1, -1
         if (estk_gc(2,(j+1)/2).gt.0) then
            nez = estk_gc(2,(j+1)/2)
            goto 22
         endif
      enddo
22    write(uni,1012) 'Gurson', nez - estk_gc(1,1) + 1
      do j=jmc, 1, -1
         if (estk_c(2,(j+1)/2).gt.0) then
            nez = estk_c(2,(j+1)/2)
            goto 23
         endif
      enddo
23    write(uni,1012) 'c', nez - estk_c(1,1) + 1
      do j=jms, 1, -1
         if (estk_s(2,(j+1)/2).gt.0) then
            nez = estk_s(2,(j+1)/2)
            goto 24
         endif
      enddo
24    write(uni,1012) 's', nez - estk_s(1,1) + 1
      do i=imb, 1, -1
         if (estk_b(2,(kma-1)/2,(i+1)/2).gt.0) then
            nez = estk_b(2,(kma-1)/2,(i+1)/2)
            goto 25
         endif
      enddo
25      write(uni,1012) 'a+b',  nez - estk_a(1,1,1) + 1
1011  format(t4,a,t14,a)
1012  format(t6,a,t16,i6)
c
      write(uni,'(/t5,a)') 'element layers/thickness :'
      write(uni,'(t5,a,t20,a,t40,a)')
     &          'elem. layer','y','dy (thickness)'
      write(uni,'(t10,i2,t17,e16.8,t37,a)') 0, y(0), '  -'
      do i=1, lt
	   write(uni,'(t10,i2,t17,e16.8,t37,e16.8)') i,y(i),y(i)-y(i-1)
      enddo
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
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
C----67--1---------2---------3---------4---------5---------6---------712
c
