	program mesh_ssy
c
c     *********************************************************************
c     *                                                                   *
c     *        Small Scale Yielding (SSY) with a Blunt Notch              *
c     *                    2D/3D Mesh Generator                           *
c     *                                                                   *
c     *    This program generates a 3D/2D FE model for SSY with           *
c     *    a blunt notch.  The output of this prog. is a full set         *
c     *    of WARP3D input files that define the FE geometry,             *
c     *    applied symmetry and elastic K or J field displacement         *
c     *    constraints along with the necessary J-domain definitions.     *
c     *                                                                   *
c     *    The user can specify the number of layers through thickness    *
c     *    in order to generate a full 3D model.  See example input       *
c     *    file for the definitions of the input and options.             *
c     *                                                                   *
c     *       written by  : Rami Haj-ALi (rh), Arne Gullerud (AG), and    *
c     *                     Robert Dodds, Jr. (rhd)                       *
c     *                                                                   *
c     *           Dept. of Civil Engineering                              *
c     *           University of Illinois at Urbana-Champaign              *
c     *                                                                   *
c     *                 (c)     All Rights Reserved                       *
c     *                                                                   *
c     *                  modified by :                                    *
c     *                  last modified : 04/03/97                         *
c     *                                                                   *
c     *********************************************************************
c
	implicit none
	integer i, j, k, max_nodes, max_elm, max_layers
	integer max_surfn, max_surfe, elem_type
	integer num_nod, num_elem, num_quads, num_nod_quad
	parameter (max_nodes=25000, max_surfn=25000,  max_layers=20)
	parameter (max_elm=25000, max_surfe=25000, elem_type=8 )
	integer nodes(max_nodes), elem(elem_type, max_elm)
	integer quad_surfe( 4, max_surfe, max_layers+1)
	integer quad_surfn(max_surfn, max_layers+1)
	real*8  quad_coord(3,max_surfn,max_layers+1)
	real*8 pi, zero, tol
	parameter (pi=3.14159265358979)
	real*8 coords(3,max_nodes)
	real*8 x0,y0, z
	data zero, tol /0.0, 1.0E-16/
c
c	input vars
c
        integer inpf
	integer mteta, num_layers, num_r, max_r, num_dr
	parameter(max_r=25)
	integer mr(max_r)
	real*8 Rn, l1, l2, R1(max_r), t(max_layers+1), biase(max_r)
	real*8 Rt
	common /input1/ num_layers, num_r, mr, mteta, num_dr
	common /input2/ x0, y0, z, Rn, Rt, l1, l2, R1, t, biase
        character*80 inp_file
c
c		WARP3D parameters
c              Warp input ( material and constraints definitions )
c
	real*8 appl_J, appl_K, t_ratio, E, nu, sig0
	integer num_steps, num_output, iw_face1, iw_face2,ple
	character*80 warp_file, warp_struct, warp_crd,
     &               warp_elem, warp_const, warp_out, warp_domain
	common /K_field/ appl_J, appl_K, t_ratio
	common /warp1/ E, nu, sig0
	common /warp2/ num_steps, num_output, ple
	common /warp3/ iw_face1, iw_face2
	common /warp4/ warp_file, warp_struct, warp_crd,
     &                 warp_elem, warp_const, warp_out, warp_domain
c
c	Input 
c
      write(*,*) 'Input file name ?'
      read(*,'(a)') inp_file
      inpf=1
      open(inpf,file=inp_file)
c
      call input_vars(inpf)
      close(1)
c
c	Loop over layers/surfaces and generate surface elements/nodes
c
        num_nod = 0
	do i=1,num_layers+1
c		  Location of root notch
	 z = t(i)
c
      call gen_surface( i, quad_surfn, quad_surfe, quad_coord,
     &                  max_surfn, max_surfe, max_layers, 
     &                  num_nod, num_quads, num_nod_quad )
c
        num_nod = num_nod + num_nod_quad
c
	end do
c
c	generate HEX8 brick elements
c
      num_elem = 0      
      do k=1,num_layers
       do i=1,num_quads
	num_elem = num_elem + 1
	do j=1,4
	elem(j,num_elem)   = quad_surfe(j,i,k+1)
	elem(j+4,num_elem) = quad_surfe(j,i,k)
	end do
       end do
      end do
c
      num_nod = 0
      do k=1, num_layers+1
       do i=1,num_nod_quad
	num_nod = num_nod + 1
	nodes( num_nod ) = quad_surfn(i,k)
	coords(1,num_nod) = quad_coord(1,i,k)
	coords(2,num_nod) = quad_coord(2,i,k)
	coords(3,num_nod) = quad_coord(3,i,k)
       end do
      end do
c
      call  output_warp(nodes, coords, elem, quad_surfn,
     &                  quad_coord,
     &                  num_nod, num_elem, max_surfn,
     &                  max_layers, num_nod_quad)
c
      stop
      end



      subroutine output_warp(nodes, coords, elem, quad_surfn,
     &                       quad_coord,
     &                       num_nod, num_elem, max_surfn,
     &                       max_layers1, num_nod_quad)
      implicit real*8(a-h,o-z)
      integer quad_surfn(max_surfn, max_layers1+1)
      real*8  quad_coord(3, max_surfn, max_layers1+1)
      dimension nodes(*), coords(3,*)
      integer elem(8,*)
c
c	input vars
c
	integer mteta, num_layers, num_r, max_r, num_dr
	parameter(max_r=25, max_layers=20)
	integer mr(max_r)
	real*8 Rn, l1, l2, R1(max_r), t(max_layers+1), biase(max_r)
	real*8 Rt
	common /input1/ num_layers, num_r, mr, mteta, num_dr
	common /input2/ x0, y0, z, Rn, Rt, l1, l2, R1, t, biase
c
c		WARP3D parameters
c              Warp input ( material and constraints definitions )
c
	real*8 appl_J, appl_K, t_ratio, E, nu, sig0
	integer num_steps, num_output, iw_face1, iw_face2,ple
	character*80 warp_file, warp_struct, warp_crd,
     &               warp_elem, warp_const, warp_out, warp_domain
	common /K_field/ appl_J, appl_K, t_ratio
	common /warp1/ E, nu, sig0
	common /warp2/ num_steps, num_output, ple
	common /warp3/ iw_face1, iw_face2
	common /warp4/ warp_file, warp_struct, warp_crd,
     &                 warp_elem, warp_const, warp_out, warp_domain
c
c	Loacal Variables
c
	parameter(max_bc=80000, max_nd_out=10000)
	dimension xyrt(4,max_bc), uv_k(2,max_bc)
	dimension node_uvc(max_bc), node_v0(max_bc), node_w0(max_bc)
	dimension nd_out(max_nd_out)
	character*80 dumc1,dumc2,dumc3
c
	data zero /0.0/
c

      open(1,file=warp_crd)
      do i=1,num_nod
      write(1,1000) nodes(i),(coords(k,i),k=1,3)
      end do
      close(1)
c
      open(1,file=warp_elem)
      do i=1,num_elem
      write(1,1100) i,(elem(k,i),k=1,8)
      end do
      close(1)
c
c
c    Locate nodes for boundary condtions
c
c	V=0  bc
c
      num_v0 = 0
      do k=1,num_layers+1
      num_v0 = num_v0 + 1
      node_v0( num_v0 ) = quad_surfn(1,k)
      num_v0 = num_v0 + 1
      node_v0( num_v0 ) = quad_surfn(1+mteta,k)
      num_v0 = num_v0 + 1
      node_v0( num_v0 ) = quad_surfn(1+2*mteta,k)
       do i=4, num_dr
        num_v0 = num_v0 + 1
        nadd = 1+ 2*mteta + (i-3)*(mteta+1)
        node_v0( num_v0 ) = quad_surfn(nadd,k)
       end do
      end do
c
c      W=0 bc
c
      num_w0 = 0
c
c	face 1  ( Z=t(1) )
c
      if(iw_face1.eq.0) then
       do i=1,num_nod_quad
       num_w0 = num_w0 + 1
       node_w0( num_w0 ) = quad_surfn(i,1)
       end do
      end if 
c
c	face 2  ( Z=t( num_layers + 1 ) )
c
      if(iw_face2.eq.0) then
       do i=1,num_nod_quad
       num_w0 = num_w0 + 1
       node_w0( num_w0 ) = quad_surfn(i,num_layers+1)
       end do
      end if 
c
c    U and V from applied_K field
c
c   1) Prepare the nodes and their polar coordinates for calculations
c
        num_uvc = 0
       do k = 1, num_layers+1
        nadd = 1 + 2*mteta + (num_dr-3)*(mteta+1) - 1
        do i=1, mteta+1
        num_uvc = num_uvc + 1
        node_uvc( num_uvc) = quad_surfn( nadd+i, k)
        xr = quad_coord(1, nadd+i, k)
        yr = quad_coord(2, nadd+i, k)
        xyrt(1, num_uvc) = xr
        xyrt(2, num_uvc) = yr
        xyrt(3, num_uvc) = dsqrt( xr*xr + yr*yr )
        call get_teta(xr, yr, tr_rad, tr_deg )
        xyrt(4, num_uvc) = tr_rad
        end do
       end do
c
c
c       Output constrain file
c
      open(1,file=warp_const)
      write(1,10) 
      write(1,12) 
      write(1,10) 
      do i=1, num_v0
       write(1,22) node_v0(i), zero
      end do
c
      write(1,10) 
      write(1,13) 
      write(1,10) 
      do i=1, num_w0
       write(1,23) node_w0(i), zero
      end do
c
      call get_k_field(xyrt, uv_k, num_uvc, E, nu,
     &                 appl_J, appl_K, ple)
c
      write(1,10) 
      write(1,14) 
      write(1,15) appl_J, appl_k
      write(1,10) 
      do i=1, num_uvc
       if(uv_k(2,i).ne.0.0) then
       write(1,21) node_uvc(i),uv_k(1,i),uv_k(2,i)
       else
       write(1,24) node_uvc(i),uv_k(1,i)
       end if
      end do
c
      close(1)
c
c	Domain Integral definition
c
      iring1 = 0.75*float(num_dr)
      iring2 = 0.9*float(num_dr)
c
      open(1,file=warp_domain)
      if( num_layers.eq.1) then
c
c		2D Model
c
        dumc3= 'domain dssy'
        call make_name(dumc3,1,dumc1)
	write(1,'(a)') dumc1
	write(1,32)
	write(1,33)
	write(1,34)  1
c	write(1,30) ( quad_surfn(i,1), i=1,mteta-1)
        do i=1,mteta-1
	nd_out(i) = quad_surfn(i,1)
	end do
	call warp_node_out(nd_out,mteta-1,1)
	write(1,34)  2
c	write(1,30) ( quad_surfn(i,2), i=1,mteta-1)
        do i=1,mteta-1
	nd_out(i) = quad_surfn(i,2)
	end do
	call warp_node_out(nd_out,mteta-1,1)
	write(1,37) 1, 2
c	write(1,38) iring1, iring2
	write(1,38) 
	write(1,39)
	write(1,40)
	write(1,42)
c
c		3D Model
c
	else
c
c
        dumc3= 'domain dssy'
        call make_name(dumc3,1,dumc1)
	write(1,'(a)') dumc1
	write(1,32)
	write(1,33)
	write(1,34) 1
c	write(1,30) ( quad_surfn(i,1), i=1,mteta-1)
        do i=1,mteta-1
	nd_out(i) = quad_surfn(i,1)
	end do
	call warp_node_out(nd_out,mteta-1,1)
	write(1,34) 2
c	write(1,30) ( quad_surfn(i,2), i=1,mteta-1)
        do i=1,mteta-1
	nd_out(i) = quad_surfn(i,2)
	end do
	call warp_node_out(nd_out,mteta-1,1)
	write(1,37) 1,2
c	write(1,38) iring1, iring2
	write(1,38) 
	write(1,39)
	write(1,40)
	write(1,42)
c
        do k=2, num_layers
        call make_name(dumc3,k,dumc1)
	write(1,'(a)') dumc1
	write(1,32)
	write(1,33)
	write(1,34) k-1
c	write(1,30) ( quad_surfn(i,k-1), i=1,mteta-1)
        do i=1,mteta-1
	nd_out(i) = quad_surfn(i,k-1)
	end do
	call warp_node_out(nd_out,mteta-1,1)
	write(1,34) k
c	write(1,30) ( quad_surfn(i,k), i=1,mteta-1)
        do i=1,mteta-1
	nd_out(i) = quad_surfn(i,k)
	end do
	call warp_node_out(nd_out,mteta-1,1)
	write(1,34) k+1
c	write(1,30) ( quad_surfn(i,k+1), i=1,mteta-1)
        do i=1,mteta-1
	nd_out(i) = quad_surfn(i,k+1)
	end do
	call warp_node_out(nd_out,mteta-1,1)
	write(1,35) k-1, k, k+1
c	write(1,38) iring1, iring2
	write(1,38) 
	write(1,44)
	write(1,40)
	write(1,42)
	end do
c
c
        k=num_layers+1
        call make_name(dumc3,k,dumc1)
	write(1,'(a)') dumc1
	write(1,32)
	write(1,33)
	write(1,34) k-1
c	write(1,30) ( quad_surfn(i,k-1), i=1,mteta-1)
        do i=1,mteta-1
	nd_out(i) = quad_surfn(i,k-1)
	end do
	call warp_node_out(nd_out,mteta-1,1)
	write(1,34) k
c	write(1,30) ( quad_surfn(i,k), i=1,mteta-1)
        do i=1,mteta-1
	nd_out(i) = quad_surfn(i,k)
	end do
	call warp_node_out(nd_out,mteta-1,1)
	write(1,37) k-1, k
c	write(1,38) iring1, iring2
	write(1,38) 
	write(1,45)
	write(1,40)
	write(1,42)
c
       end if
c
       close(1)
c
c
c	Warp main input file
c
        open(1,file=warp_file)
c
c
	istring = index(warp_struct,' ')
	if(istring.eq.0) then
	write(*,*) '>>>> Fatal Error'
	stop
	end if
c
	istring = istring - 1
        dumc2 =  "*input from file '" // warp_struct (1:istring) 
       write(1,10)
       write(1,101)  Rn, Rt
       write(1,10)
       write(1,10)
       dumc1 =  'structure ' // warp_struct (1:40)
       write(1,'(a)') dumc1
       write(1,102)
       write(1,103) E, nu, sig0
       write(1,104)
       write(1,10)
       write(1,10)
       write(1,105) num_nod
       write(1,106) num_elem
       write(1,10)
       write(1,107)
        dumc1 =  dumc2 (1:18+istring) // ".crd' " 
       write(1,10)
       write(1,8)
       write(1,'(a)') dumc1
       write(1,10)
       write(1,9)
       write(1,108)
       write(1,10)
       write(1,109)
       write(1,110) 1,num_elem
       write(1,111) 
       write(1,10)
       write(1,112) 
        dumc1 =  dumc2 (1:18+istring) // ".elm' " 
       write(1,10)
       write(1,8)
       write(1,'(a)') dumc1
       write(1,10)
       write(1,9)
       write(1,10)
       write(1,113)
c
        nblock=128
	iblock = dfloat(num_elem)/nblock
	do i=1,iblock
	write(1,1200) i,nblock,1+(i-1)*nblock
	end do
	ilast=nblock*iblock
	if(ilast.lt.num_elem) then
	write(1,1200) iblock+1,num_elem-ilast,1+ilast
	end if
c
        write(1,10)
        write(1,114)
       write(1,10)
       write(1,8)
        dumc1 =  dumc2 (1:18+istring) // ".const' " 
        write(1,'(a)') dumc1
       write(1,10)
       write(1,9)
        write(1,10)
        write(1,10)
c
        write(1,10) 
        write(1,14) 
        write(1,15) appl_J, appl_k
        write(1,10) 
        write(1,10)
c
        write(1,201)
        write(1,202)
        write(1,203)
        write(1,10)
        write(1,10)
        write(1,204)
        write(1,205)
	write(1,206) num_steps
        write(1,10)
        write(1,10)
        write(1,501)
        write(1,502)
        write(1,10)
        write(1,10)
        write(1,207)
        write(1,208)
        write(1,209)
        write(1,210)
        write(1,211)
c        write(1,212)
        write(1,213)
        write(1,214)
        write(1,215)
        write(1,216)
        write(1,217)
        write(1,218)
        write(1,220)
c
	write(1,10)
	write(1,10)
c
        dumc1 =  dumc2 (1:18+istring) // ".output' " 
        dumc3 =  dumc2 (1:18+istring) // ".domain' " 
	iold = 0
       do i=num_output, num_steps, num_output
	write(1,10)
        write(1,301) iold+1,i
        write(1,'(a)') dumc1
        write(1,'(a)') dumc3
	write(1,10)
	iold = i
       end do
c
        write(1,302)
	close(1)
c
c
c	Generate Output file
c
	open(1,file=warp_out)
	write(1,10)
	write(1,401)
	write(1,402)
	write(1,403)
	write(1,10)
	close(1)
c
8	format('*echo off')
9	format('*echo on')
10	format('c')
12	format('c',8x,'v=0   constraints')
13	format('c',8x,'w=0   constraints')
14	format('c',8x,'u and v K_field constraints')
15	format('c',15x,'Applied J increment is J_inc = ',e17.9,/,
     &         'c',15x,'Applied K increment is K_inc = ',e17.9)
21	format(i7,2x,'u',2x,e17.9,4x,'v',2x,e17.9)
24	format(i7,2x,'u',2x,e17.9)
22	format(i7,2x,'v',2x,e17.9)
23	format(i7,2x,'w',2x,e17.9)
c
30	  format(2x,i4,20(',',i4))
32        format('    symmetric')
33        format('    normal plane  ny  1.0')
34        format('    node set ',i3,',')
37        format('    front node sets ',i2,x,i2,x,'linear ')
35        format('    front node sets ',i2,x,i2,x,i2,x,'linear ')
c 38        format('    q-values automatic rings ',i3,' - ',i3)
38        format('    q-values automatic rings 50-70 ')
39        format('    function type a')
40        format('    print totals')
42        format('    compute domain integral')
44        format('    function type b')
45        format('    function type c')
c
101    format('c    SSY 3D Model with Blunt Notch  Rn = ',
     &           f9.6,'     R= ',f15.7)
102    format('material steel')
103    format('    properties  mises  e ',f6.0,'  nu ',f3.1,
     &        '  yld_pt ',f9.4,'  n_power 6.,')
104    format('                rho 0.0')
105    format('number of nodes ',i8)
106    format('number of elements ',i8)
107    format('coordinates')
108    format('elements')
109    format('c   for config number   0')
110    format(14x,i1,' - ',i5,' type l3disop nonlinear material steel,')
111    format(23x,'order 2x2x2 bbar center_output short')
112    format('incidences')
113    format('blocking')
114    format('constraints')
c
c
201     format ('loading unit_j')
202     format (' nodal_loads')
203     format ('    1  force_y 0.0')
c
204     format(' loading load1')
205     format('  nonlinear')
206     format('    step 1 -',i3,' unit_j 1.0')
501     format('c output patran neutral')
502     format('c stop')
c
c
207     format(' dynamic analysis parameters')
208     format('   solution technique direct sparse')
209     format('c   preconditioner type diagonal')
210     format('c   lnr_pcg conv test res tol 0.01')
211     format('c   maximum linear iterations 2000')
212     format('   stiffness updates before iterations all')
213     format('   maximum iterations 20')
214     format('   convergence test norm res tol .005 max res tol 0.05')  
215     format('   adaptive on')
216     format('   time step 100000')
217     format('   nonconvergent solutions stop')
218     format('   trace solution on lpcg_solution off')
219     format('   accelerate off')
220     format('   extrapolate on')
c
301     format('compute displacements for loading load1 steps ',
     &          i3,' - ',i3)
302     format('stop')
c
c
401     format(' output patran binary displacements')
402     format(' output patran binary element strains')
403     format(' output patran binary element stresses')
c
c
1000	format(i9,3(2x,e17.9))
1100    format(9(1x,i6))
1200	format(i7,3x,i3,3x,i9)
c
	return
	end


	subroutine warp_node_out(nd,num_nod,out)
	implicit real*8(a-h,o-z)
	dimension nd(*)
	parameter(nd_line=12)
	integer out
c
        if(num_nod.lt.nd_line) then
	write(out,100) (nd(k),k=1,num_nod)
	return
100	format(3x,i5,11(',',i5))
	end if
c
	nline1 = (dfloat(num_nod)-0.1)/dfloat(nd_line)
	nrest = num_nod - nline1*nd_line
	if(nrest.eq.0) then
	goto 1
	end if
	 do i=1,nline1
	  ncount = 1 + (i-1)*nd_line
	  write(out,200) (nd(k),k=ncount,ncount+nd_line-1)
	 end do
	 ncount = 1 + (nline1)*nd_line
	 if(nrest.eq.1) then
	 write(out,301) (nd(k),k=ncount, num_nod)
	 end if
	 if(nrest.eq.2) then
	 write(out,302) (nd(k),k=ncount, num_nod)
	 end if
	 if(nrest.eq.3) then
	 write(out,303) (nd(k),k=ncount, num_nod)
	 end if
	 if(nrest.eq.4) then
	 write(out,304) (nd(k),k=ncount, num_nod)
	 end if
	 if(nrest.eq.5) then
	 write(out,305) (nd(k),k=ncount, num_nod)
	 end if
	 if(nrest.eq.6) then
	 write(out,306) (nd(k),k=ncount, num_nod)
	 end if
	 if(nrest.eq.7) then
	 write(out,307) (nd(k),k=ncount, num_nod)
	 end if
	 if(nrest.eq.8) then
	 write(out,308) (nd(k),k=ncount, num_nod)
	 end if
	 if(nrest.eq.9) then
	 write(out,309) (nd(k),k=ncount, num_nod)
	 end if
	 if(nrest.eq.10) then
	 write(out,310) (nd(k),k=ncount, num_nod)
	 end if
	 if(nrest.eq.11) then
	 write(out,311) (nd(k),k=ncount, num_nod)
	 end if
	 if(nrest.eq.12) then
	 write(out,312) (nd(k),k=ncount, num_nod)
	 end if
	 return
c
1	continue
c
	 do i=1,nline1-1
	  ncount = 1 + (i-1)*nd_line
	  write(out,200) (nd(k),k=ncount,ncount+nd_line-1)
	 end do
	 ncount = 1 + (nline1)*nd_line
	 write(out,210) (nd(k),k=ncount, num_nod)
c
200	format(3x,12(i5,','))
210	format(3x,i5,11(',',i5))
c
301	format(3x,i5)
302	format(3x,i5,(',',i5))
303	format(3x,i5,2(',',i5))
304	format(3x,3(i5,','),i5)
305	format(3x,i5,4(',',i5))
306	format(3x,5(i5,','),i5)
307	format(3x,i5,6(',',i5))
308	format(3x,7(i5,','),i5)
309	format(3x,i5,8(',',i5))
310	format(3x,9(i5,','),i5)
311	format(3x,i5,10(',',i5))
312	format(3x,11(i5,','),i5)
c

	return
	end


       subroutine get_k_field(xyrt, uv_k, num_uvc, E, nu,
     &                        applied_J, applied_K, ple)
       implicit real*8(a-h,o-z)
       real*8 nu, kappa
       integer ple
       parameter (pi=3.14159265358979)
       dimension xyrt(4,*), uv_k(2,*)
c
        G = E/(2.0*(1.0+nu))
c
	if(ple.eq.1) then
c
c		Plane Strain
c 
         H = E / ( 1.0 -nu*nu )
	 applied_k = dsqrt ( H * applied_J )
	 kappa = 3.0 - 4.0*nu
c
	else
c
c	       Plane Stress
c
         H = E 
	 applied_k = dsqrt ( H * applied_J )
	 kappa = (3.0 - nu) / ( 1.0 + nu )
c
	end if
c
	term1 = applied_K / (2.0*G)
c
	do i=1, num_uvc
	      r = xyrt(3,i)
	 teta_r = xyrt(4,i)
	 ct = dcos(teta_r)
	 st = dsin(teta_r)
	 cto2 = dcos(teta_r/2.0)
	 sto2 = dsin(teta_r/2.0)
c
	term2 = dsqrt(r/(2.0*pi))
c
c
         u = term1*term2*cto2*( kappa -1.0 + 2.0*sto2*sto2)   
         v = term1*term2*sto2*( kappa +1.0 - 2.0*cto2*cto2)   
c
	 uv_k(1,i) = u
	 uv_k(2,i) = v
c
         end do
	 return
	 end


      subroutine gen_hex8_layer( quad1, quad2, elem, num_quads)
      implicit real*8(a-h,o-z)
      integer quad1(4,*), quad2(4,*),elem(8,*)
       do i=1,num_quads
	do k=1,4
	elem(k,i)   = quad2(k,i)
	elem(k+4,i) = quad1(k,i)
	end do
       end do
c
       return
       end

      subroutine gen_surface( ilr, quad_surfn, quad_surfe, quad_coord,
     &                 max_surfn, max_surfe, max_layersl, 
     &                 num_nodg, num_quads, num_nod )
c
	implicit real*8(a-h,o-z)
	integer quad_surfe(4, max_surfe, max_layersl+1)
	integer quad_surfn(max_surfn, max_layersl+1)
	real*8  quad_coord(3,max_surfn,max_layersl+1)
	parameter (pi=3.14159265358979)
	data zero, tol /0.0, 1.0E-16/
c
c	input vars
c
	integer mteta, num_layers, num_r, max_r, num_dr
	parameter(max_r=25, max_layers=20)
	integer mr(max_r)
	real*8 Rn, l1, l2, R1(max_r), t(max_layers+1), biase(max_r)
	real*8 Rt
	common /input1/ num_layers, num_r, mr, mteta, num_dr
	common /input2/ x0, y0, z, Rn, Rt, l1, l2, R1, t, biase
c
c	Loacal Variables
c
	parameter(max_dr=1500, max_nteta=1500)
	dimension dr(max_dr)
	dimension iquadl(4,max_nteta)
	dimension nd_crv1(max_nteta), nd_crv2(max_nteta)
c
c	Notch nodes:  first two lines of nodes
c
	num_dr = 3
	do i=1,num_r
	num_dr = num_dr + mr(i)
	end do
c
c
	dr(1) = Rn
	dr(2) = dr(1) + l1
	dr(3) = dr(2) + l2
	ii=3
	do k=1,num_r
	 Sum=1.0
	 do i=2,mr(k)
	  Sum = Sum + biase(k)**( dfloat(i-1) )
         end do
	 a1 = R1(k) / Sum
	 do i=1,mr(k)
	  ii= ii + 1
	  an = a1 * biase(k)**( dfloat(i-1))
	  dr(ii) = dr(ii-1) + an
	 end do
	end do
c
	Rt = dr(ii)
c
	dteta_r=(pi/2.0)/( dfloat(mteta)/2.0 )
	num_nod = 0
	teta_r=0.0
	do i=1,mteta/2 + 1
	 teta_r = (dfloat(i)-1.0)*dteta_r
	 num_nod = num_nod + 1
	 quad_surfn(num_nod,ilr) = num_nod + num_nodg
	 quad_coord(1,num_nod,ilr) = x0 + dr(1)*dcos(teta_r)
	 quad_coord(2,num_nod,ilr) = y0 + dr(1)*dsin(teta_r)
	 quad_coord(3,num_nod,ilr) = z
	end do
c
	dl= ( dr(2)+0.5*l2 )/(dfloat(mteta)/2.0 - 1.0 )
	do i=1,mteta/2 -1
	 num_nod = num_nod + 1
	 quad_surfn(num_nod,ilr) = num_nod + num_nodg
	 quad_coord(1,num_nod,ilr) = x0 - dl*dfloat(i)
	 quad_coord(2,num_nod,ilr) = y0 + dr(1)
	 quad_coord(3,num_nod,ilr) = z
	end do
c
c	Second circum. nodes
c
	teta_r=0.0
	do i=1,mteta/2 + 1
	 teta_r = (dfloat(i)-1.0)*dteta_r
	 num_nod = num_nod + 1
	 quad_surfn(num_nod,ilr) = num_nod + num_nodg
	 quad_coord(1,num_nod,ilr) = x0 + dr(2)*dcos(teta_r)
	 quad_coord(2,num_nod,ilr) = y0 + dr(2)*dsin(teta_r)
	 quad_coord(3,num_nod,ilr) = z
	end do
c
c
	do i=1,mteta/2 -2
	 num_nod = num_nod + 1
	 quad_surfn(num_nod,ilr) = num_nod + num_nodg
	 quad_coord(1,num_nod,ilr) = x0 - dl*dfloat(i)
	                        xl = x0 - dl*dfloat(i)
	 quad_coord(2,num_nod,ilr) = y0 +
     &           ellipse_y( xl, dl*(dfloat(mteta)/2.0-1.0), dr(2))
	   yl =  ellipse_y( xl, dl*(dfloat(mteta)/2.0-1.0), dr(2))
	 quad_coord(3,num_nod,ilr) = z
	end do
c
c
	 i=mteta/2 -1
	 num_nod = num_nod + 1
	 quad_surfn(num_nod,ilr) = num_nod + num_nodg
	 quad_coord(1,num_nod,ilr) = x0 - dl*dfloat(i)
	                        xl = x0 - dl*dfloat(i)
	 quad_coord(2,num_nod,ilr) = y0 + dr(1) +  0.4*(yl-dr(1))
	 quad_coord(3,num_nod,ilr) = z
c
c	Loop over the rest of the arcs
c
	do k = 3, num_dr
c
c	    Right Side
c
c
	dteta_r=(pi/2.0)/( dfloat(mteta)/2.0 )
	teta_r=0.0
	do i=1,mteta/2 + 1
	 teta_r = dfloat(i-1)*dteta_r
	 num_nod = num_nod + 1
	 quad_surfn(num_nod,ilr) = num_nod + num_nodg
	 quad_coord(1,num_nod,ilr) = x0 + dr(k)*dcos(teta_r)
	 quad_coord(2,num_nod,ilr) = y0 + dr(k)*dsin(teta_r)
	 quad_coord(3,num_nod,ilr) = z
	end do
c
c	   Left Side
c
	xl = dsqrt( dr(k)*dr(k) - dr(1)*dr(1) )
	alpha = dAtan( dr(1)/xl )
	teta_l = (PI/2.0) - alpha
	dteta_l = teta_l / dfloat( mteta/2 )
c
	do i=1,mteta/2 -1
	 num_nod = num_nod + 1
	 teta_l = dfloat(i)*dteta_l
	 quad_surfn(num_nod,ilr) = num_nod + num_nodg
	 quad_coord(1,num_nod,ilr) = x0 - dr(k)*dsin(teta_l)
	 quad_coord(2,num_nod,ilr) = y0 + dr(k)*dcos(teta_l)
	 quad_coord(3,num_nod,ilr) = z
	end do
c
c	Last point on the straight line
c
	 num_nod = num_nod + 1
	 quad_surfn(num_nod,ilr) = num_nod + num_nodg
	 quad_coord(1,num_nod,ilr) = x0 - xl
	 quad_coord(2,num_nod,ilr) = y0 + dr(1)
	 quad_coord(3,num_nod,ilr) = z
c
	end do
c
c
c		End of Generating nodes
c
c	Begin generating surface elements (quads )
c
c	   First layer of elements
c
      num_quads = 0
      do i = 1, mteta
       nd_crv1(i) = quad_surfn(i,ilr)
       nd_crv2(i) = quad_surfn(i+mteta,ilr)
      end do
c
      call quads_2crvs(nd_crv1, nd_crv2, iquadl, mteta, num_quadi)
      do i=1, num_quadi
       ii = num_quads + i
       do j=1,4
	 quad_surfe(j,ii, ilr) = iquadl(j,i)
       end do
      end do
      num_quads = num_quads + num_quadi
c
c	second layer of elements
c
	ndum1 = nd_crv1(mteta)
	do i=1,mteta
	nd_crv1(i) = nd_crv2(i)
	end do
	nd_crv1(mteta+1) = ndum1
c
      do k = 1, num_dr-2
       nc = 2*mteta + (mteta+1)*(k-1)
       do i = 1, mteta + 1
        nd_crv2(i) = quad_surfn(i+nc,ilr)
       end do
c
      call quads_2crvs(nd_crv1, nd_crv2, iquadl, mteta+1, num_quadi)
      do i=1, num_quadi
       ii = num_quads + i
       do j=1,4
	 quad_surfe(j,ii, ilr) = iquadl(j,i)
       end do
      end do
      num_quads = num_quads + num_quadi
c
      do i=1,mteta+1
       nd_crv1(i) = nd_crv2(i)
      end do
c
      end do
c
      return
      end


      subroutine quads_2crvs(nd_crv1, nd_crv2, ielm, ndt, num_elm)
       implicit real*8(a-h,o-z)
       dimension nd_crv1(*), nd_crv2(*), ielm(4,*)
       num_elm = ndt - 1
	do i= 1, num_elm
	 ielm(1,i) = nd_crv2(i)
	 ielm(2,i) = nd_crv2(i+1)
	 ielm(3,i) = nd_crv1(i+1)
	 ielm(4,i) = nd_crv1(i)
	end do
c
	return
	end


	real*8 function ellipse_y(x,a,b)
	    real*8 x, a, b
	    ellipse_y = dsqrt( b*b -(b*b/(a*a))*x*x )
	    return
	end




      subroutine input_vars(inpf)
       implicit none
       integer i
       character*80 dumc1
       character overwrt*1
       integer inpf
       logical file_exist
c
c	input vars
c
	integer mteta, num_layers, num_r, max_r, max_layers, num_dr
        integer ifirst, ilast, ierr
	parameter(max_r=25, max_layers=20)
	integer mr(max_r)
	real*8 Rn, l1, l2, R1(max_r), t(max_layers+1), biase(max_r)
	real*8 x0, y0, z
	real*8 Rt
	common /input1/ num_layers, num_r, mr, mteta, num_dr
	common /input2/ x0, y0, z, Rn, Rt, l1, l2, R1, t, biase
c
c		WARP3D parameters
c              Warp input ( material and constraints definitions )
c
	real*8 appl_J, appl_K, t_ratio, E, nu, sig0
	integer num_steps, num_output, iw_face1, iw_face2, ple
	character*80 warp_file, warp_struct, warp_crd,
     &               warp_elem, warp_const, warp_out, warp_domain
	common /K_field/ appl_J, appl_K, t_ratio
	common /warp1/ E, nu, sig0
	common /warp2/ num_steps, num_output, ple
	common /warp3/ iw_face1, iw_face2
	common /warp4/ warp_file, warp_struct, warp_crd,
     &                 warp_elem, warp_const, warp_out, warp_domain
c
c
c	Start Reading Input Parameters
c
c     Rn - root notch radius    l1, l2 - length of elemnts 1 and 2 in r-dir
c
      call inp_back(inpf,ifirst,ilast,ierr)
      read(inpf,*) Rn, l1, l2
c
c     num_r - number of segments in the r-dir that will
c             be divided (aside from (l1+l2)
c
      call inp_back(inpf,ifirst,ilast,ierr)
      read(inpf,*) num_r
c
c     R(k) (k=1,num_r) -- line segments' length in the r-dir
c     mr(k) (k=1,num_r) -- Number of elements for line segments
c     biase(k)  --  geometric series mesh grading biase :
c                           a(n+1)/a(n) = biase**(n-1)
      do i = 1, num_r
       call inp_back(inpf,ifirst,ilast,ierr)
       read(inpf,*) R1(i), mr(i), biase(i)
      end do
c
c	ADJUST R1(1) --->   R1(1) - ( Rn + l1 + l2 )
c       This is in order to have a final R_total that is "nice number"
c
	R1(1) = R1(1) - ( Rn + l1 + l2 )
c
c   mteta - number of elements circumfrentially ( Must be even !!!! )
c
      call inp_back(inpf,ifirst,ilast,ierr)
      read(inpf,*) mteta
c
c    num_layers - number of elements through the thickness ( z-direction )
c
      call inp_back(inpf,ifirst,ilast,ierr)
      read(inpf,*) num_layers
c
c     t(k) - Z-ccord of surface k ; where k=1,num_layers+1
c            it is recomended that Z coordinate be negative
c
c    Location of root notch ( Note: this can be changed in the future )
c
	 x0 = 0.0
	 y0 = 0.0
c
      do i=1, num_layers+1
       call inp_back(inpf,ifirst,ilast,ierr)
       read(inpf,*) t(i)
      end do
c
c                              WARP3D
c
c     Input warp structure nmae (xxx)
c     This will generate the following files:
c       xxx.inp  - warp input file
c       xxx.crd  - coordinates file
c       xxx.elm  - elements id file
c       xxx.const -  applied displacement constriants
c       xxx.output - output requests at specified intervals
c       xxx.domain - J integral domain and output definitions 
c
       call inp_back(inpf,ifirst,ilast,ierr)
       read(inpf,'(a)')  dumc1
       warp_struct = dumc1 (ifirst:ilast)
       warp_file   = dumc1 (ifirst:ilast) // '.inp'
       warp_crd    = dumc1 (ifirst:ilast) // '.crd'
       warp_elem   = dumc1 (ifirst:ilast) // '.elm'
       warp_const  = dumc1 (ifirst:ilast) // '.const'
       warp_out    = dumc1 (ifirst:ilast) // '.output'
       warp_domain = dumc1 (ifirst:ilast) // '.domain'
c
c           Check if output files exist
c
       file_exist = .false.
       inquire (file = warp_file  , exist = file_exist )
       inquire (file = warp_crd   , exist = file_exist )
       inquire (file = warp_elem  , exist = file_exist )
       inquire (file = warp_const , exist = file_exist )
       inquire (file = warp_out   , exist = file_exist )
       inquire (file = warp_domain, exist = file_exist )
       if(file_exist) then
       write(*,*)
       write(*,*) '>>> Some or all of warp output file already exist '
88     write(*,*)
       write(*,*) '>>> Overwrite these files? [y]es / [n]o ? '
       write(*,*)
       read(*,'(a)') overwrt
        if(overwrt .eq. 'n') then
          write(*,*)
          write(*,*) '>>> Program terminated'
          write(*,*)
	  stop
        elseif (overwrt .ne. 'y') then
	goto 88
	end if
       end if
c
c
c      Input w boundary conditions on the two surfaces/faces (Z=0 and Z=B or B/2)
c             Boundary_flag = 0  --- Apply w=0 constrain
c             Boundary_flag = 1  --- Traction Free Surface ( NO W constraints are applied )
c
c      Note: for plane-strain case w=0 constraint should be applied on booth faces
c
c        Input:  iw_face1   iw_face2
c
       call inp_back(inpf,ifirst,ilast,ierr)
       read(inpf,*)  iw_face1, iw_face2
c
c	Input:  ple, num_steps, num_output
c                ple     -- (Integer) flag:  =1 - plane strain .else.  plane stress
c                num_steps- Number of solution steps
c                num_output-Interval of load steps Number for warp files output and calc. J
c
       call inp_back(inpf,ifirst,ilast,ierr)
       read(inpf,*) ple, num_steps, num_output
c
c	Input: E, nu, sig0
c	Where:
c                E       -- Young's Modulus
c                nu      -- Piosson's Ratio
c                sig0    -- Yield stress
c
       call inp_back(inpf,ifirst,ilast,ierr)
       read(inpf,*)  E, nu, sig0
c
c	Input: appl_J
c	Where:
c                appl_J -- Increment of J if the material were elastic
c                appl_K -- will be calculated from appl_J for pl-e or pl-s
c
       call inp_back(inpf,ifirst,ilast,ierr)
       read(inpf,*)  appl_J
c
      return
      end




	SUBROUTINE GET_TETA(COORD1,COORD2,TETAR,TETAD)
	IMPLICIT REAL*8 (A-H,O-Z)
C
	pi=3.14159265358979
	TETAR=0.
	XN=DABS ( COORD1 )
	YN=DABS ( COORD2 )
	IF(COORD1.GT.0..AND.COORD2.GT.0.) THEN
	TETAR=ATAN( YN/XN )
	END IF
C
	IF(COORD1.LT.0..AND.COORD2.GT.0.) THEN
	TETAR=ATAN( YN/XN )
	TETAR= PI - TETAR
	END IF
C
	IF(COORD1.LT.0..AND.COORD2.LT.0.) THEN
	TETAR=ATAN( YN/XN )
	TETAR= PI + TETAR
	END IF
C
	IF(COORD1.GT.0..AND.COORD2.LT.0.) THEN
	TETAR=ATAN( YN/XN )
	TETAR= 2.*PI - TETAR
	END IF
C
	IF(COORD1.EQ.0.) THEN
	  IF(COORD2.GT.0. ) THEN
	  TETAR=PI/2.
	  ELSE
	  TETAR=3.*PI/2.
	  END IF
	END IF
C
	IF(COORD2.EQ.0.) THEN
	  IF(COORD1.GT.0. ) THEN
	  TETAR=0.
	  ELSE
	  TETAR=PI
	  END IF
	END IF
C
	TETAD=TETAR*180/PI
	RETURN
	END



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
	IF( DUMC (1:1) .EQ.'C '.OR. DUMC (1:1) .EQ.'c ') THEN
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




	subroutine make_name(part1,n,fname)
C
C	THIS SUB. ATTACHES AN INTEGER NUMBER (3 Digits) TO A CHARCTER
C
	implicit real*8 (a-h,o-z)
	character part1*80,fname*80
	character num(10)*1
c
c
	num(1)="0"
	num(2)="1"
	num(3)="2"
	num(4)="3"
	num(5)="4"
	num(6)="5"
	num(7)="6"
	num(8)="7"
	num(9)="8"
	num(10)="9"
c
c	
	a= dfloat(n)/100.
	ia=int(a)
	n2 = n - ia * 100
	b = dfloat(n2) / 10.
	ib = int (b)
	n3 = n2 - ib * 10
	ic = n3
c
	j=index (part1,"   ")
c	fname = part1 (1:j-1)//' '//num(ia+1)//num(ib+1)//num(ic+1)
	fname = part1 (1:j-1)//'_'//num(ia+1)//num(ib+1)//num(ic+1)
c	fname = part1 (1:j-1)//num(ia+1)//num(ib+1)//num(ic+1)
c
	return
	end

