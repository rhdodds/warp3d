c
c  this program filters converts the patran compatible WARP result file
c  to TECPLOT loadable data file
c
c
      program tecplot_patran
      implicit integer (a-z)
      parameter(n_max=500)
      parameter(num_node=95000)
      parameter(num_elem=50000)
      parameter(maxcols=26)
c
      common   /text/ textstring
c
      double precision ux(num_node), uy(num_node), uz(num_node)
      double precision x(3, num_node), elval(maxcols, num_elem)
      real xval, yval, zval, pvals(maxcols), rtemp
      character*9 ch_disp(n_max), ch_strn(n_max), ch_strs(n_max)
      character*12  dispfile(n_max), stressfile(n_max),
     &     strainfile(n_max)
      character binnam*80, ok*1, cstep*6, yes*6, strsname*80, struc*30,
     &          strnname*80
      character textstring*80, coordinates*30, incidences*30, dispname*12
      dimension title(80), ele(20, num_elem)
      dimension nstep(n_max), ino(num_node), iel(num_elem),
     &          lex(num_elem)
c
c
c
c--------------------------------------
c
      write(*,*)
      write(*,*)' --------------------------------------------------'
      write(*,*)'| Program prepares data for Tecplot 360 program    |'
      write(*,*)'|                                                  |'
      write(*,*)'| The input files is: web* wnb*                    |'
      write(*,*)'| The output files are: disp*.dat strs*.dat        |'
      write(*,*)'|                       and strn*.dat              |'
      write(*,*)'|                                                  |'
      write(*,*)'|                                                  |'
      write(*,*)'| Computational Fracture Mechanics Research Group  |'
      write(*,*)'| University of Illinois at Urbana-Champaign       |'
      write(*,*)'|                                                  |'
      write(*,*)'| Last modification: June  18, 2007                |'
      write(*,*)' --------------------------------------------------'
      write(*,*)
c
c--------------------------------------
c
c
      write(*,*) 'How many steps at which warp result files exist?'
      read(*,*) nn
      write(*,*) 'For every ? steps the results were output?'
      read(*,*) nk
3     write(*,*) 'nodal stress or strain (s for stress, n for strain)?'
      read(*,'(a)') ok
         if(ok.ne.'s'.and.ok.ne.'S'.and.ok.ne.'n'
     &   .and.ok.ne.'N') then
          write(*,*) 'please enter "s" for stress or "n" for strain'
          go to 3
         endif
         write(*,*) 'output for nodal displacement, (yes or no)?'
         read(*,'(a)') yes
      write(*,*) 'Total element number in the model?'
      read(*,*) nelem
      write(*,*) 'Total node number in the model?'
      read(*,*) n_nod
         write(*,*) 'prefix of the input file?'
         read(*,'(a)') struc
5        write(*,*) '8 or 20 node elements?'
         read(*,*) ncon
c
         if(ncon.ne.8.and.ncon.ne.20)then
          write(*,*) 'please enter 8 or 20!'
          go to 5
         endif
c
c
      do i=1,nn
      nstep(i) = nk*i
      write(cstep, 100) nstep(i)
      ch_disp(i)='wnbd'//cstep(2:6)
      dispfile(i)='disp'//cstep(3:6)//'.dat'
      stressfile(i)='strs'//cstep(3:6)//'.dat'
      strainfile(i)='strn'//cstep(3:6)//'.dat'
       ch_strn(i)='wnbe'//cstep(2:6)
       ch_strs(i)='wnbs'//cstep(2:6)
      enddo
100   format('c',i5.5)
c
c   Read in the nodal coordinates and element connectivity
c
      textstring=struc//'   '
      call left_centred(30,irp)
      struc=textstring(1:30)
      il=1
10       il=il+1
      if (il.lt.irp) then
      if (struc(il:il).ne.'.') goto 10
      il=il-1
      endif
c
      coordinates = struc(1:il)//'.crd'
      incidences = struc(1:il)//'.elm'
c
      open(unit=1,file=coordinates,status='old')
      do i=1,n_nod
          read(1,*) ii, x(1,ii), x(2,ii), x(3,ii)
         enddo
         close(unit=1)
      open(unit=1,file=incidences,status='old')
      do i=1,nelem
          if(ncon.eq.8)then
          read(1,*) ii, ele(1,ii), ele(2,ii), ele(3,ii), ele(4, ii),
     &  ele(5, ii), ele(6, ii), ele(7, ii), ele(8,ii)
          else
          read(1,*) ii, ele(1,ii), ele(2,ii), ele(3,ii), ele(4, ii),
     &  ele(5, ii), ele(6, ii), ele(7, ii), ele(8,ii), ele(9,ii),
     &  ele(10,ii), ele(11,ii), ele(12,ii), ele(13,ii), ele(14,ii),
     &  ele(15,ii), ele(16,ii), ele(17,ii), ele(18,ii), ele(19,ii),
     &  ele(20, ii)
          endif
         enddo
         close(unit=1)
c
c *********************************************************
c *                                                       *
c *    read displacement binary patran file               *
c *                                                       *
c *********************************************************
c
      if((yes(1:1).eq.'y').or.(yes(1:1).eq.'Y')) then
c
      termin=5
      termot=6
      binfil=10
         tec=15
c
      write(termot,*) ' '
      write(termot,*) ' '
      write(termot,*) '>> binary node value processing program'
      write(termot,*) ' '
c
      do 200 i=1,nn
      open(unit=tec,file=dispfile(i),status='unknown')
         write(tec,*) 'Variables= "X", "Y","Z", "Node UserID",
     &"U", "V", "W","Step ID"'
         write(tec,*) 'Zone T="nodaldisp", N=', n_nod,',','E=',nelem,','
     &  ,'datapacking=point, zonetype=FEBRICK'
       binnam = ch_disp(i)
       open(unit=binfil,file=binnam,status='old',recl=3000,
     &      form='unformatted')
      write(termot,*) ' '
      write(termot,'(a,a,a)') ' > file ', ch_disp(i), ' open ok'
c
c           read the binary results file of nodal values.
c           read x, y, z components. patran results are single
c           precision. read as single and store as double.
c
      read(binfil) title, nnode, ii, rtemp, ii, nvals
      read(binfil) title
      read(binfil) title
c
c            read values for each node into a double array
c
      write(termot,*) '> reading nodal results file..'
      do node = 1,nnode
        read(binfil) ii, xval, yval, zval
        ux(node) = xval
       uy(node) = yval
       uz(node) = zval
       ino(node) = ii
c
      if(ii.ne.node) then
       write(*,*) 'error during reading nodal results !!!'
       stop
      endif
c
c => n_step  node-number  u1   u2   u3
      write(tec,150)  x(1,node), x(2,node),x(3,node),
     & ino(node), ux(node), uy(node), uz(node),nstep(i)
      enddo
150   format(3(2x,g15.6),(2x,i6),3(2x,g15.6),(2x,i6))
         do elem=1, nelem
         write(tec, 180) (ele(ii,elem), ii=1,8)
         enddo
180   format(8(2x,i6))
      close(unit=tec)
      close(unit=binfil)
      write(termot,*) '> writing nodal results to file..'
200   continue
      endif
c
c ************************************************************
c *                                                          *
c *  Read in stress strain  binary patran files              *
c *                                                          *
c ************************************************************
c
      write(termot,*) ' '
      write(termot,*) ' '
      write(termot,*) '>> binary stress/strain processing program'
      write(termot,*) ' '
c
      if((ok(1:1).eq.'s') .or. (ok(1:1).eq.'S')) then
      do 300 i=1,nn
      binnam = ch_strs(i)
      strsname=stressfile(i)
      open(unit=tec,file=strsname,status='unknown')
      write(tec,*) 'Variables="X","Y","Z", "NodeID"'
      write(tec,*) '"Sxx","Syy","Szz","Sxy","Syz","Sxz",'
      write(tec,*) '"U0","mises", "c1", "c2", "c3", "I1", "I2","I3",'
      write(tec,*) ' "S1","S2","S3","l1","m1","n1","l2","m2","n2", '
          write(tec,*) '"l3","m3","n3","Step ID"'
          write(tec,*) 'Zone T="nodalstress", N=', n_nod, 'E=',nelem,
     & 'datapacking=point, zonetype=FEBRICK'
       write(*,*)'Processing file: ',binnam
       open(unit=binfil,file=binnam,status='old',recl=3000,
     &     form='unformatted')
      write(termot,*) ' '
      write(termot,'(a,a,a)') ' > file ', ch_strs(i), ' open ok'
c
c           read the binary results file of nodal strains/stresses.
c           patran results are single precision. read as single and
c           store as double.
c
      read(binfil) title, nnode, ii, rtemp, ii, nvals
      read(binfil) title
      read(binfil) title
c
c            read values for each node into a double array
c
      write(termot,*) '> reading nodal stress/strain file..'
      write(termot,*) '> writing nodal stress/strain to file..'
c
c
      do node = 1,nnode
       read(binfil) ii, (pvals(jj),jj=1,nvals)
       ino(node) = ii
      if(ii.ne.node) then
       write(*,*) 'error during reading nodal results !!!'
       stop
      endif
c
      write(tec,250) x(1,node), x(2,node),x(3,node),
     &  ino(node),(pvals(jj),jj=1,nvals),nstep(i)
c
      enddo 
250   format(3(2x,g15.6)(2x,i6),26(2x,g15.6),(2x,i6))
      close(unit=binfil)
         do elem=1, nelem
         write(tec, 280) (ele(ii,elem), ii=1,8)
         enddo
280   format(8(2x,i6))
      close(unit=tec)
c
300   continue
c
      else
c
      do 400 i=1,nn
      binnam = ch_strn(i)
      strnname=strainfile(i)
      open(unit=tec,file=strnname,status='unknown')
      write(tec,*) 'Variables="X","Y","Z","NodeID",'
      write(tec,*) '"Exx","Eyy","Ezz","Exy","Eyz","Exz","Eeff",'
      write(tec,*) '"I1", "I2", "I3", "e1", "e2", "e3","l1","m1",'
      write(tec,*) '"n1","l2","m2","n2","l3","m3","n3", "STEP ID"'
          write(tec,*) 'Zone T="nodestrain", N=',n_nod,'E=',nelem,
     & 'datapacking=point, zonetype=FEBRICK'
      write(*,*)'Processing file: ',binnam
      open(unit=binfil,file=binnam,status='old',recl=3000,
     &     form='unformatted')
      write(termot,*) ' '
      write(termot,'(a,a,a)') ' > file ', ch_strn(i), ' open ok'
c
c           read the binary results file of nodal strains/stresses.
c           patran results are single precision. read as single and
c           store as double.
c
      read(binfil) title, nnode, ii, rtemp, ii, nvals
      read(binfil) title
      read(binfil) title
c
c            read values for each node into a double array
c
      write(termot,*) '> reading nodal stress/strain file..'
      write(termot,*) '> writing nodal stress/strain to file..'
c
c
      do node = 1,nnode
        read(binfil) ii, (pvals(jj),jj=1,nvals)
        ino(node) = ii
      if(ii.ne.node) then
       write(*,*) 'error during reading nodal results !!!'
       stop
      endif
c
      write(tec,350) x(1,node), x(2,node),x(3,node),
     &  ino(node),(pvals(jj),jj=1,nvals),nstep(i)
c
      enddo
350   format(3(2x,g15.6),(2x,i6),22(2x,g15.6),(2x,i6))
      close(unit=binfil)
         do elem=1, nelem
         write(tec, 380) (ele(ii,elem), ii=1,8)
         enddo
380   format(8(2x,i6))
      close(unit=tec)
c
400   continue
      endif
c
         stop
         end
c
c ********************************************************************
c
c     this subroutine rewrites the string from the first non-empty
c     letter to the last letter, and also gives the total number of
c     letters in the string
c
c ********************************************************************
c----67--1---------2---------3---------4---------5---------6---------7-!
      subroutine left_centred(nch,rpo)
      implicit none
c
      integer   nch,rpo,l,r,i
      character textstring*80
c
      common   /text/ textstring
c
      l=0
10    l=l+1
      if ( textstring(l:l).eq.' ') goto 10
      r=l
20    r=r+1
      if (r.lt.nch) then
        if (textstring(r:r).ne.' ') goto 20
        r=r-1
      endif
      rpo=0
      do i=l, r
        rpo=rpo+1
        textstring(rpo:rpo)=textstring(i:i)
      enddo
      return
      end
c
