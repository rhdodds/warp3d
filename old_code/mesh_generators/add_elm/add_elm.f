	program mesh_scp_add
C----------------------------------------------------------------------
C          A program to add layers of elements into a presmatic type   |
C	       finite element model.                                   |
C	Written by:   Rami M. Haj-Ali and  Robert H. Dodds, Jr.        |
C         date:        8/8/1996                                        |
C       Last modified:    ............     by: ...............         |
C                                                                      |
C----------------------------------------------------------------------
	implicit real*8(a-h,o-z)
	parameter (node_max=30000, ielm_max=30000, node_per_el=8)
	parameter (node_surf_max=1000,layers_max=20)
	parameter (num_face=6,num_el_block=64)
	parameter (Tol=0.001)
	dimension nodes(node_max),coords(3,node_max)
	dimension id_el(node_per_el,ielm_max)
	dimension iels(2,ielm_max),iface(4,num_face)
	dimension nodef(4)
C
	dimension nodes_new(node_max),coords_new(3,node_max)
	dimension id_el_new(node_per_el,ielm_max)
	dimension iels_new(2,ielm_max)
	dimension nodesl(node_max)
	dimension node_surface(layers_max,node_surf_max)
C
	logical debug
	character inp_file*80,inp_crd*80,inp_els*80,out_crd*80,out_els*80
	character in_abaqus*80,out_nsets*80
	common/files/inp_file,inp_crd,inp_els,out_crd,out_els
	common/files2/in_abaqus,out_nsets
	common/check1/debug
	data iface/1,2,3,4,5,6,7,8,1,2,6,5,4,3,7,8,2,3,7,6,1,4,8,5/
C
	debug=.false.
C
	write(*,*) 'Input file name ?'
	read(*,'(a)') inp_file
	inpf=1
	open(inpf,file=inp_file)
C
	call input(inpf,ixyz,xyzi,xyzf,biase,layers,nodal_sets)
	close(1)
C
	if(debug) then
	 write(*,*) 'DEBUG>>>>> (1)  End of input reading'
	 write(*,*)
	end if
C
	inp_c=1
	inp_e=2
	iout_c=3
	iout_e=4
C
	open(inp_c,file=inp_crd)
	open(inp_e,file=inp_els)
	open(iout_c,file=out_crd)
	open(iout_e,file=out_els)
C
C
	if(debug) then
	 write(*,*) 'DEBUG>>>>> (2) End of open/create input/output files'
	 write(*,*)
	end if
C
	num_nod=0
	num_elm=0
C
	do i=1,node_max
	read(inp_c,*,end=10) nodes(i),coords(1,i),coords(2,i),coords(3,i)
	num_nod=num_nod+1
	end do
C
 10	if(num_nod.ne.nodes(num_nod)) then
	 write(*,*) '>>>>>ERROR-1 -  Number of nodes is not sequential',
     &   num_nod,nodes(num_nod)
	 stop
	end if
C
	do i=1,ielm_max
	read(inp_e,*,end=20) iels(1,i),(id_el(k,i),k=1,node_per_el)
	num_elm=num_elm+1
	end do
C
 20	if(num_elm.ne.iels(1,num_elm)) then
	 write(*,*) '>>>>>ERROR-2 -Number of elements is not sequential',
     &   num_elm,iels(1,num_elm)
	 stop
	end if
C
	close(inp_c)
	close(inp_e)
C
	do i=1,num_elm
	   iels(2,i)=0
	end do
C
C Start searching elements faces that belong to the specified surface
C
	do 200 i=1,num_elm
	 do 100 j= 1,num_face
	   nodef(1) = id_el( iface(1,j) , i )
	   nodef(2) = id_el( iface(2,j) , i )
	   nodef(3) = id_el( iface(3,j) , i )
	   nodef(4) = id_el( iface(4,j) , i )
C
	   iflag=1
	      do 50 k=1,4
	       if( coords(ixyz, nodef(k)).GE.xyzi-TOL.AND.
     &            coords(ixyz,nodef(k)).LE.xyzi+TOL) then
	        goto 50
	       end if
	       iflag=0
	       goto 100
  50	      continue
	   iels(2,i)=j
C
	if(debug) then
	   write(*,*) iels(1,i),iels(2,i)
	end if
C
	   goto 200
 100	 continue
 200     continue
C
C
	if(debug) then
	 write(*,*) 'DEBUG>>>>> Calling subroutine add_new_nodes'
	 write(*,*)
	end if
C
	call   add_new_nodes(nodes,coords,id_el,iels,iface,
     &                    nodes_new,coords_new,id_el_new,iels_new,
     &                    ixyz,xyzi,xyzf,biase,layers,
     &                    node_surface,nodesl,
     &                    num_nod,num_elm,
     &                    num_nod_new,num_elm_new,
     &                    node_max,ielm_max,node_per_el,num_face,tol,
     &                    node_surf_max,layers_max,num_nod_surf)
C
	write(*,1000)
	write(*,1010) num_nod
	write(*,1020) num_elm
	write(*,1021) ixyz,xyzi
	write(*,1022) layers,ixyz,xyzi,xyzf,biase
	write(*,1030) num_nod_surf
	write(*,1040) num_nod_new
	write(*,1050) num_elm_new
	write(*,1060) num_nod+num_nod_new
	write(*,1070) num_elm+num_elm_new
1000    format(//)
1010    format(5x,'Number of nodes in the model...........',I6,/)
1020    format(5x,'Number of elements in the model........',I6,/)
1021    format(5x,'Identified edge surface with  X(',I1,')= ',F8.2,/)
1022    format(5x,'New ',I2,' layers of elements are added where:',//,
     & 5x,'X(',I1,')= from ',F8.2,' to ',F8.2,' with biase =', F4.2,/)
1030    format(5x,'Number of nodes on the given edge...',I6,/)
1040    format(5x,'Number of added new nodes...........',I6,/)
1050    format(5x,'Number of added new elements........',I6,/)
1060    format(5x,'Total number of nodes...............',I6,/)
1070    format(5x,'Total number of elements............',I6,/)
C
C
	call  output_mesh(nodes,coords,id_el,iels,
     &                    nodes_new,coords_new,id_el_new,iels_new,
     &	                  iout_c,iout_e,
     &                    num_nod,num_elm,
     &                    num_nod_new,num_elm_new,
     &                    node_max,ielm_max,node_per_el,num_el_block)
C
	if(nodal_sets.ne.0) then
C
	call  output_nsets(nodes,coords,id_el,iels,
     &                     nodes_new,coords_new,id_el_new,iels_new,
     &	                   iout_c,iout_e,
     &                     num_nod,num_elm,
     &                     num_nod_new,num_elm_new,
     &                     node_max,ielm_max,node_per_el,num_face,
     &                     tol,
     &                     in_abaqus,out_nsets)
C
	end if
	STOP
	END
	
	


	subroutine add_new_nodes(nodes,coords,id_el,iels,iface,
     &                      nodes_new,coords_new,id_el_new,iels_new,
     &                      ixyz,xyzi,xyzf,biase,layers,
     &                      node_surface,nodesl,
     &                      num_nod,num_elm,
     &                      num_nod_new,num_elm_new,
     &                      node_max,ielm_max,node_per_el,num_face,tol,
     &                      node_surf_max,layers_max,num_nod_surf)
	implicit real*8(a-h,o-z)
C
	dimension nodes(node_max),coords(3,node_max)
	dimension nodes_new(node_max),coords_new(3,node_max)
	dimension id_el(node_per_el,ielm_max)
	dimension id_el_new(node_per_el,ielm_max)
	dimension iels_new(2,ielm_max)
	dimension iels(2,ielm_max),iface(4,num_face)
	dimension nodesl(node_max)
	dimension node_surface(layers_max,node_surf_max)
C
	logical debug
	common/check1/debug
C
C
C	Local arrays
C
	dimension nodef(4),node8(8),iface2(6)
	data iface2/2,1,4,3,6,5/
C
C iface2 - denote the oposite face number definition for hexagonal element
C	   Note: Currently this defintion corresponds to WARP3D element
C
	do i=1, node_max
	    nodesl(i)=0
	end do
C
	do i=1,node_surf_max
	  do l=1,layers_max
	    node_surface(l,i)=0
	  end do
	end do
C
	num_nod_surf=0
C
C	Collect the surface nodes
C
C
	if(debug) then
	 write(*,*) 'DEBUG/add_new_nodes/ start collecting nodes '
	 write(*,*)
	end if
C
	do 30 i=1,num_elm
	 j=iels(2,i)
C
	 if(j.ne.0) then
	   nodef(1) = id_el( iface(1,j) , i )
	   nodef(2) = id_el( iface(2,j) , i )
	   nodef(3) = id_el( iface(3,j) , i )
	   nodef(4) = id_el( iface(4,j) , i )
	   do 20 i1=1,4
	     iadd=1
	     do 10 i2=1,num_nod_surf
	     if(nodef(i1).eq.node_surface(1,i2)) then
	      iadd=0
	      goto 20
	     end if
 10	    continue
	    num_nod_surf=num_nod_surf + 1
	    node_surface(1,num_nod_surf) = nodef(i1)
	    nodesl( nodef(i1) ) = num_nod_surf
 20	   continue
	 end if
C
 30	continue
C
C        Create new nodes by mapping the specified surface nodes
C
C
	if(debug) then
	write(*,*) 'num of surface nodes =', num_nod_surf
	 write(*,*) 'DEBUG/add_new_nodes/ start creating nodes/els '
	 write(*,*)
	end if
C
	num_nod_new=0
	num_elm_new=0
C
	Sum=1.0
	do i=2,layers
	 Sum = Sum + biase ** ( float(i-1))
	end do
	a1 = ( xyzf - xyzi ) / ( Sum )
C
	do 40  i = 1, layers
C
	 if(i.eq.1) then
	  an = a1
	  Sumi=a1
	 else
	  an = a1 * biase**( float(i-1) )
	  Sumi = Sumi + an
	 end if
C
	  do j =1, num_nod_surf
	   num_nod_new=num_nod_new + 1
	   nodes_new(num_nod_new) = num_nod + num_nod_new
	   coords_new(1,num_nod_new) = coords(1, node_surface(1,j) )
	   coords_new(2,num_nod_new) = coords(2, node_surface(1,j) )
	   coords_new(3,num_nod_new) = coords(3, node_surface(1,j) )
	   coords_new(ixyz,num_nod_new) = coords(ixyz, node_surface(1,j) )
     &                                    + Sumi
	   node_surface(i+1,j)= nodes_new( num_nod_new )
	  end do
 40     continue
C
C	Start definig the new elements
C
	do 100 i=1,layers
	do 50 ii = 1, num_elm
	 ifacel1=iels(2,ii)
C
	  if(ifacel1.ne.0) then
	    ns1 = id_el( iface (1, ifacel1), ii)
	    ns2 = id_el( iface (2, ifacel1), ii)
	    ns3 = id_el( iface (3, ifacel1), ii)
	    ns4 = id_el( iface (4, ifacel1), ii)
C
	    nsl1 = nodesl( ns1 )
	    nsl2 = nodesl( ns2 )
	    nsl3 = nodesl( ns3 )
	    nsl4 = nodesl( ns4 )
C
	    node8(1) = node_surface(i, nsl1)
	    node8(2) = node_surface(i, nsl2)
	    node8(3) = node_surface(i, nsl3)
	    node8(4) = node_surface(i, nsl4)
C
	    node8(5) = node_surface(i+1, nsl1)
	    node8(6) = node_surface(i+1, nsl2)
	    node8(7) = node_surface(i+1, nsl3)
	    node8(8) = node_surface(i+1, nsl4)
C
	    num_elm_new = num_elm_new +1
	    iels_new(1,num_elm_new) = num_elm + num_elm_new
	    do k=1,8
	    id_el_new(k,num_elm_new) = node8(k)
	    end do
	  end if
C
 50     continue
C
	if(debug) then
	 write(*,*) 'DEBUG/add_new_nodes/ finish layer nmber',i
	 write(*,*)
	end if
C
100    continue
	return
	end
	
	  
	
	
	
	


	subroutine input(inpf,ixyz,xyzi,xyzf,biase,layers,nodal_sets)
	implicit real*8(a-h,o-z)
	character inp_file*80,inp_crd*80,inp_els*80,out_crd*80,out_els*80
	character in_abaqus*80,out_nsets*80
	character dumc1*80,dumc2*80
	logical debug
	common/check1/debug
	common/files/inp_file,inp_crd,inp_els,out_crd,out_els
	common/files2/in_abaqus,out_nsets
C
C	1)The FE model must be presmatic in shape in order to properly
C	  locate the elements and nodes on the desired edge.
C
C	2) The elements must be 8-node brick (3D) elements
C
C
C	ixyz= (1,2 or 3 ) the coord direction needed for extending the mesh
C
	call inp_back(inpf,ifirst,ilast,ierr)
	read(inpf,*) ixyz
C
Cinput xyzi and xyzf -  The initial existing coordinate and the new final coord
C
	call inp_back(inpf,ifirst,ilast,ierr)
	read(inpf,*) xyzi,xyzf
C
Cinput the biase-mesh grading and the number of layers or elements to be added
C
	call inp_back(inpf,ifirst,ilast,ierr)
	read(inpf,*) biase,layers
C
C	input 'xxx' string for file name xxx.crd and xxx.elm
C
	call inp_back(inpf,ifirst,ilast,ierr)
	read(inpf,'(a)') dumc1
        inp_crd = dumc1 (ifirst:ilast) // '.crd'
        inp_els = dumc1 (ifirst:ilast) // '.elm'
C
        out_crd = dumc1 (ifirst:ilast) // '_mod.crd'
        out_els = dumc1 (ifirst:ilast) // '_mod.elm'
C
	inpo=11
	inpn=12
        dumc2  = dumc1 (ifirst:ilast) // '.inp'
	open(inpo,file= dumc2 )
        dumc2  = dumc1 (ifirst:ilast) // '_mod.inp'
	open(inpn,file= dumc2 )
C
C	input flag for reading nodal sets and outputing new lists
C		nodal_sets = .EQ. 0  - no reading/generating new sets
C			     ELSE - reading and generating nodal_sets
C
	call inp_back(inpf,ifirst,ilast,ierr)
	read(inpf,*) nodal_sets
C
	if(nodal_sets.ne.0) then
C
C	input abaqus input file generated from mesh3d_scp.f program
C
	  call inp_back(inpf,ifirst,ilast,ierr)
	  read(inpf,'(a)') in_abaqus
	  call inp_back(inpf,ifirst,ilast,ierr)
	  read(inpf,'(a)') out_nsets
	end if
C
	if(debug) then
	 write(*,*) 'DEBUG>>>>> (0)  End of reading input data '
	 write(*,*)
	end if
C
	return
	end



	SUBROUTINE INP_BACK(INP,IFIRST,ILAST,IERR)
	IMPLICIT REAL*8(A-H,O-Z)
	PARAMETER (MAXC=20)
C
C	Skip input-file lines that:
C			1) begin with 'C ' or 'c '
C			2) first 20 characters are blank
C
	CHARACTER DUMC*80
C
	IFIRST=0
	ILAST=0
C
	IERR=1
10	READ(INP,'(A)',END=999) DUMC
	IF( DUMC (1:2) .EQ.'C '.OR. DUMC (1:2) .EQ.'c ') THEN
	GOTO 10
	END IF
C
	IF(INDEX ( DUMC (1:20), '                    ' ) .NE.0) THEN
	GOTO 10
	END IF
C
	DO I = 1,MAXC
	  IF(DUMC (I:I) .NE. ' ') THEN
	  IFIRST=I
	  GOTO 20
	  END IF
	END DO
C
 20	CONTINUE
	ILAST = IFIRST + INDEX(  DUMC (IFIRST:80), ' ') -2
 	BACKSPACE (INP)
	IERR=0
999	IF(IERR.EQ.1) THEN
	 WRITE(*,*) '>>>> ERROR(0)  Premature End of Input File '
	END IF
  	RETURN
	END



	SUBROUTINE INP_BACK2(INP,CINP,LENGTH,IERR)
	IMPLICIT REAL*8(A-H,O-Z)
C
C	Skip input-file lines until CINP is found
C
	CHARACTER DUMC*80,CINP*80
C
	IERR=1
	IFLAG=0
10	READ(INP,'(A)',END=999) DUMC
	IFLAG = INDEX (DUMC, CINP (1:LENGTH) )
	IF(IFLAG.EQ.0) THEN
	GOTO 10
	ELSE 
	IERR=0
C	BACKSPACE (INP)
	END IF
C
999	IF(IERR.EQ.1) THEN
	 WRITE(*,*) '>>>> ERROR(0)  Premature End of Input File '
	END IF
  	RETURN
	END



	subroutine output_mesh(nodes,coords,id_el,iels,
     &                    nodes_new,coords_new,id_el_new,iels_new,
     &	                  iout_c,iout_e,
     &                    num_nod,num_elm,
     &                    num_nod_new,num_elm_new,
     &                    node_max,ielm_max,node_per_el,num_el_block)
	implicit real*8(a-h,o-z)
C
	dimension nodes(node_max),coords(3,node_max)
	dimension nodes_new(node_max),coords_new(3,node_max)
	dimension id_el(node_per_el,ielm_max)
	dimension id_el_new(node_per_el,ielm_max)
	dimension iels_new(2,ielm_max)
	dimension iels(2,ielm_max)
	dimension nodesl(node_max)
	character dumc1*80
C
	logical debug
	common/check1/debug
C
C
C	Local arrays
C
C		Output modefied nodes and elements
	do i=1,num_nod
	write(iout_c,100) nodes(i),(coords(k,i),k=1,3)
	end do
C
	do i=1,num_nod_new
	write(iout_c,100) nodes_new(i),(coords_new(k,i),k=1,3)
	end do
C
	do i=1,num_elm
	write(iout_e,200) iels(1,i),(id_el(k,i),k=1,8)
	end do
C
	do i=1,num_elm_new
	write(iout_e,200) iels_new(1,i),(id_el_new(k,i),k=1,8)
	end do
C
100	format(i6,3(1x,e16.9))
200	format(10i6)
C
	close(iout_c)
	close(iout_e)
C
	inpo=11
	inpn=12
C
400     continue
	read(inpo,'(a)') dumc1
	if( dumc1 (1:8) .NE. 'blocking' ) then
	  if( dumc1 (1:15) .EQ. 'number of nodes' ) then
	  write(inpn,*) 'number of nodes  ',
     &    num_nod + num_nod_new
	  goto 400
	  end if
	  if( dumc1 (1:18) .EQ. 'number of elements' ) then
	  write(inpn,*) 'number of elements  ',
     &    num_elm + num_elm_new
	  goto 400
	  end if
	write(inpn,'(a)') dumc1
	goto 400
	else 
	write(inpn,'(a)') dumc1
	goto 600
	end if
C
600	continue
c
	nblck=num_el_block
	do i=1,9999999
	  nc=(i*nblck - (nblck -1) )
	  if( (nc+nblck).GT.(num_elm+num_elm_new) ) then
	  nblckf = num_elm+num_elm_new - ((i-1)*nblck )
	  write(inpn,1000) i,nblckf,nc
	  goto 650
	  end if
	  write(inpn,1000) i,nblck,nc
	end do
650	write(inpn,1010)
C
800	continue
	read(inpo,'(a)') dumc1
	if( index ( dumc1, 'step cases:' ) .EQ. 0 ) then
	goto 800
	else
	write(inpn,1010)
	write(inpn,'(a)') dumc1
	end if
C
900     continue
	read(inpo,'(a)',end=9999) dumc1
	write(inpn,'(a)') dumc1
	goto 900
C	
1000    format(i8,x,i8,2x,i10)
1010    format('c')
9999	close(inpo)
	close(inpn)
	return
	end



	subroutine  output_nsets(nodes,coords,id_el,iels,
     &                     nodes_new,coords_new,id_el_new,iels_new,
     &	                   iout_c,iout_e,
     &                     num_nod,num_elm,
     &                     num_nod_new,num_elm_new,
     &                     node_max,ielm_max,node_per_el,num_face,
     &                     tol,
     &                     in_abaqus,out_nsets)
C
	implicit real*8(a-h,o-z)
	parameter( max_ndset=5000)
C
	dimension nodes(node_max),coords(3,node_max)
	dimension nodes_new(node_max),coords_new(3,node_max)
	dimension id_el(node_per_el,ielm_max)
	dimension id_el_new(node_per_el,ielm_max)
	dimension iels_new(2,ielm_max)
	dimension iels(2,ielm_max)
	dimension nodesl(node_max)
	character dumc1*80
	character in_abaqus*80,out_nsets*80
C
C	local arrays
C
	dimension nfront(max_ndset),ntip(max_ndset),ncrsur(max_ndset)
	dimension nsym(max_ndset),nbak(max_ndset),nund(max_ndset),
     &            nsida(max_ndset),ntop(max_ndset)
C
	dimension nfront_new(max_ndset)
	dimension nsym_new(max_ndset),nbak_new(max_ndset),
     &            nsida_new(max_ndset),ntop_new(max_ndset)
C
	inp=22
	iout=23
	open(inp,file=in_abaqus)
	open(iout,file=out_nsets)
C
C       Locate NSET=NFRONT
C
	call INP_BACK2(inp,'NSET=NFRONT',11,IERR)
	read(inp,*) nfront_nd
	ipr= float( nfront_nd ) / 12.
	iprl = nfront_nd - ipr*12
	do i=1,ipr
	i1=(i-1)*12+1
	i2=i*12
	read(inp,*) (nfront(k),k=i1,i2)
	end do
	if(iprl.ne.0) then
	read(inp,*) (nfront(k),k=i2+1,i2+iprl)
	end if
C
C	Locate NSET=NTIP
C
	call INP_BACK2(inp,'NSET=NTIP',9,IERR)
	read(inp,*) ntip_nd
	ipr=float(ntip_nd ) / 12.
	iprl = ntip_nd - ipr*12
	do i=1,ipr
	i1=(i-1)*12+1
	i2=i*12
	read(inp,*) (ntip(k),k=i1,i2)
	end do
	if(iprl.ne.0) then
	read(inp,*) (ntip(k),k=i2+1,i2+iprl)
	end if
C
C	Locate NSET=NCRSUR
C
	call INP_BACK2(inp,'NSET=NCRSUR',11,IERR)
	read(inp,*) ncrsur_nd
	ipr= float( ncrsur_nd ) / 12.
	iprl = ncrsur_nd - ipr*12
	do i=1,ipr
	i1=(i-1)*12+1
	i2=i*12
	read(inp,*) (ncrsur(k),k=i1,i2)
	end do
	if(iprl.ne.0) then
	read(inp,*) (ncrsur(k),k=i2+1,i2+iprl)
	end if
C
C	Locate NSET=NSYM
C
	call INP_BACK2(inp,'NSET=NSYM',9,IERR)
	read(inp,*) nsym_nd
	ipr=float(nsym_nd) / 12.
	iprl = nsym_nd - ipr*12
	do i=1,ipr
	i1=(i-1)*12+1
	i2=i*12
	read(inp,*) (nsym(k),k=i1,i2)
	end do
	if(iprl.ne.0) then
	read(inp,*) (nsym(k),k=i2+1,i2+iprl)
	end if
C
C
C	Locate NSET=NBAK
C
	call INP_BACK2(inp,'NSET=NBAK',9,IERR)
	read(inp,*) nbak_nd
	ipr= float( nbak_nd ) / 12.
	iprl = nbak_nd - ipr*12
	do i=1,ipr
	i1=(i-1)*12+1
	i2=i*12
	read(inp,*) (nbak(k),k=i1,i2)
	end do
	if(iprl.ne.0) then
	read(inp,*) (nbak(k),k=i2+1,i2+iprl)
	end if
C
C
C	Locate NSET=NUND
C
	call INP_BACK2(inp,'NSET=NUND',9,IERR)
	read(inp,*) nund_nd
	ipr= float(nund_nd) / 12.
	iprl = nund_nd - ipr*12
	do i=1,ipr
	i1=(i-1)*12+1
	i2=i*12
	read(inp,*) (nund(k),k=i1,i2)
	end do
	if(iprl.ne.0) then
	read(inp,*) (nund(k),k=i2+1,i2+iprl)
	end if
C
C
C	Locate NSET=NSIDA
C
	call INP_BACK2(inp,'NSET=NSIDA',10,IERR)
	read(inp,*) nsida_nd
	ipr= float(nsida_nd) / 12.0
	iprl = nsida_nd - ipr*12
	do i=1,ipr
	i1=(i-1)*12+1
	i2=i*12
	read(inp,*) (nsida(k),k=i1,i2)
	end do
	if(iprl.ne.0) then
	read(inp,*) (nsida(k),k=i2+1,i2+iprl)
	end if
C
C	Locate NSET=NTOP
C
	call INP_BACK2(inp,'NSET=NTOP',9,IERR)
	read(inp,*) ntop_nd
	ipr= float(ntop_nd) / 12.
	iprl = ntop_nd - ipr*12
	do i=1,ipr
	i1=(i-1)*12+1
	i2=i*12
	read(inp,*) (ntop(k),k=i1,i2)
	end do
	if(iprl.ne.0) then
	read(inp,*) (ntop(k),k=i2+1,i2+iprl)
	end if
C
C      Locate max/min coordinates
C
	Xmax=-99999999.0
	Ymax=-99999999.0
	Zmax=-99999999.0
	Xmin=0.0
	Ymin=0.0
	Zmin=0.0
C
	do i=1,num_nod_new
	  Xc=coords_new(1,i)
	  Yc=coords_new(2,i)
	  Zc=coords_new(3,i)
	  if(Xmax.LT.Xc) then
	   Xmax=Xc
	  end if
	  if(Ymax.LT.Yc) then
	   Ymax=Yc
	  end if
	  if(Zmax.LT.Zc) then
	   Zmax=Zc
	  end if
	  if(Xmin.GT.Xc) then
	   Xmin=Xc
	  end if
	  if(Ymin.GT.Yc) then
	   Ymin=Yc
	  end if
	  if(Zmin.GT.Zc) then
	   Zmin=Zc
	  end if
	end do
C
	 ntop_nd=0
C
	do i=1,num_nod_new
	  Xc=coords_new(1,i)
	  Yc=coords_new(2,i)
	  Zc=coords_new(3,i)
C
	  if(Zc.GT.(Zmax-tol).AND.Zc.Lt.(Zmax+tol) ) then
	    ntop_nd = ntop_nd +1
	    ntop( ntop_nd ) = nodes_new( i )
	  end if
C
	  if(Yc.GT.(Ymin-tol).AND.Yc.Lt.(Ymin+tol) ) then
	    nfront_nd = nfront_nd +1
	    nfront( nfront_nd ) = nodes_new( i )
	  end if
C
	  if(Xc.GT.(Xmin-tol).AND.Xc.Lt.(Xmin+tol) ) then
	    nsym_nd = nsym_nd +1
	    nsym( nsym_nd ) = nodes_new( i )
	  end if
C
	  if(Xc.GT.(Xmax-tol).AND.Xc.Lt.(Xmax+tol) ) then
	    nsida_nd = nsida_nd +1
	    nsida( nsida_nd ) = nodes_new( i )
	  end if
C
	  if(Yc.GT.(Ymax-tol).AND.Yc.Lt.(Ymax+tol) ) then
	    nbak_nd = nbak_nd +1
	    nbak( nbak_nd ) = nodes_new( i )
	  end if
C
	end do
C
C
C
C	Output nsets
C
	write(iout,*) 'NSET=NTOP'
	write(iout,*) '**  number of ntop nodes =', ntop_nd
	ipr= float(ntop_nd) / 12.0
	iprl = ntop_nd - ipr*12
	do i=1,ipr
	i1=(i-1)*12+1
	i2=i*12
	write(iout,1000) (ntop(k),k=i1,i2)
	end do
	if(iprl.ne.0) then
	write(iout,1000) (ntop(k),k=i2+1,i2+iprl)
	end if
C
	write(iout,*) 'NSET=NFRONT'
	write(iout,*) '**  number of nfront nodes =', nfront_nd
	ipr= float(nfront_nd)/ 12.0
	iprl = nfront_nd - ipr*12
	do i=1,ipr
	i1=(i-1)*12+1
	i2=i*12
	write(iout,1000) (nfront(k),k=i1,i2)
	end do
	if(iprl.ne.0) then
	write(iout,1000) (nfront(k),k=i2+1,i2+iprl)
	end if
C
	write(iout,*) 'NSET=NSYM'
	write(iout,*) '**  number of nsym nodes =', nsym_nd
	ipr= float(nsym_nd) / 12.0
	iprl = nsym_nd - ipr*12
	do i=1,ipr
	i1=(i-1)*12+1
	i2=i*12
	write(iout,1000) (nsym(k),k=i1,i2)
	end do
	if(iprl.ne.0) then
	write(iout,1000) (nsym(k),k=i2+1,i2+iprl)
	end if
C
	write(iout,*) 'NSET=NSIDA'
	write(iout,*) '**  number of nsida nodes =', nsida_nd
	ipr= float(nsida_nd) / 12.0
	iprl = nsida_nd - ipr*12
	do i=1,ipr
	i1=(i-1)*12+1
	i2=i*12
	write(iout,1000) (nsida(k),k=i1,i2)
	end do
	if(iprl.ne.0) then
	write(iout,1000) (nsida(k),k=i2+1,i2+iprl)
	end if
C
	write(iout,*) 'NSET=NBAK'
	write(iout,*) '**  number of nbak nodes =', nbak_nd
	ipr= float(nbak_nd) / 12.0
	iprl = nbak_nd - ipr*12
	do i=1,ipr
	i1=(i-1)*12+1
	i2=i*12
	write(iout,1000) (nbak(k),k=i1,i2)
	end do
	if(iprl.ne.0) then
	write(iout,1000) (nbak(k),k=i2+1,i2+iprl)
	end if
C
C	output symmetry boundary conditions in for warp3d
C
C
	dum0=0.0
C
	write(iout,*) 
	write(iout,*) 'WARP3D    NFRONT V=0'
	write(iout,*) 
	do i=1,nfront_nd
	write(iout,1100) nfront(i),dum0
	end do
C
C
	write(iout,*) 
	write(iout,*) 'WARP3D    NSYM  U=0'
	write(iout,*) 
	do i=1,nsym_nd
	write(iout,1200) nsym(i),dum0
	end do
C
C
	dum0=5.0
	write(iout,*) 
	write(iout,*) 'WARP3D   NBACK  V=c'
	write(iout,*) 
	do i=1,nbak_nd
	write(iout,1100) nbak(i),dum0
	end do
C
1000    format(i5,11(',',i5))
1100    format(i5,2x,'v',2x,E13.8)
1200    format(i5,2x,'u',2x,E13.8)
	return
	end
