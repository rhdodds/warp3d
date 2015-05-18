c----67--1---------2---------3---------4---------5---------6---------712
c
c    File: plate_femap.f   Written July 7 1998   /Jonas Faleskog
c                          Modified on  8 1998   /JF
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine femap_neu(etyp,elnum,dmat,jobname)
c
c The routine generates an input deck for in the Neutral Patran format
c
      implicit none
c
      integer           etyp,elnum
      double precision  dmat(*)
      character  jobname*200, femapfile*200
c
      integer    io,pid,mid,etypid
      parameter (io=11)
      real       version
      character  mark*5
c
c Open FEMAP Neutral Format file.
c      Each each data block is separated by a mark = -1
c
      call appendFile( jobname, '_fmp.neu', femapfile )
      open( unit=io, file=femapfile, status='unknown' )
c
c Separator mark for each data block in the FEMAP neutral
c file, "-1" must be in columns 4 and 5
c
c Recommend using FEMAP version 4.1 data block formats
c
      mark = '   -1'
      version = 4.1
c
c Write the Neutral File Header (title and version number)
c
      call femap_header(io,mark,version)
c
c Write the node data
c
      call femap_nodes(io,mark)
c
c Write the element connectivity data
c
      pid = 1
      call femap_elements(io,mark,pid,etypid,etyp,elnum)
c
c Write material data.
c
      mid = 1
      call femap_material(io,mark,version,mid,dmat)
c
c Write element property data
c
      call femap_elemproperty(io,mark,pid,mid,etypid)
c
      close(io)
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine femap_header(io,mark,version)
c
c Write the Neutral File Header (title and version number)
c
      implicit none
      integer   io,id
      real      version
      character mark*5,title*78
c
      id = 100
      title = '3-D Surface Crack FE-mesh' 
      write(io,'(a5)')   mark
      write(io,'(i3)')   id
      write(io,'(a78)')  title
      write(io,'(f3.1)') version
      write(io,'(a5)')   mark
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine femap_nodes(io,mark)
c
c Write the node data
c
      implicit none
      include 'plate_common_nod.f'
      integer          io,id,color,typ,nodnr,i
      double precision x,y,z
      character mark*5
c
c the node data line contains:
c  node number, coord sys ID, output coord sys ID, layer ID, color ID, 
c  6 permanent constraint values, x, y, z(=0) coord, node type(=0)
c (see Appendix G, p. G-12 of the FEMAP manual)
c
      id    = 403            ! ID number for the data block
      color = 46             ! green
      typ   = 0              ! node type
      write(io,'(a5)') mark
      write(io,'(i3)') id
      do i=1, inm
         if (nnr(i).gt.0) then
            nodnr = nnr(i)
            x = npos(i,1)
            y = npos(i,2)
            z = npos(i,3)
            write(io,101) nodnr,0,0,1,color,0,0,0,0,0,0,x,y,z,typ
         endif
      enddo
      write(io,'(a5)') mark
101   format(i6,3i4,i4,6i3,3e16.8,i3)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine femap_elements(io,mark,pid,etypid,etyp,elnum)
c
c Write the element connectivity data
c
      implicit none
      include 'plate_common_nod.f'
      include 'plate_common_eln.f'
      integer   io,etyp,elnum,id,color,pid,etypid,shape,ie,j,k
      character mark*5
c
      id    = 404    ! ID number for the data block
      color = 124         ! white
      etypid = 25   ! linear solid element type number
      shape = 8      ! 8 node brick
c
      write(io,'(a5)') mark
      write(io,'(i3)') id
c
c Set property type=1, only one property data set given below
c
      do ie=1, elnum
c List first 8 nodes for solid elem., fill remainder of line with zeros
c fill in zeros for 10 values on remainding line
         write(io,101) ie,color,pid,etypid,shape,0,0,0
         write(io,102) ( nnr(eln(ie,j)), j=1, etyp ), 0,0
         write(io,102) 0,0,0,0,0,0,0,0,0,0
c
         write(io,103) (0.0,k=1,3)
         write(io,103) (0.0,k=1,3)
         write(io,103) (0.0,k=1,3)
         write(io,104) (0,k=1,16)
      enddo
      write(io,'(a5)') mark
c
101   format(8i8)
102   format(10i8)
103   format(3f12.6)
104   format(12i3,4i8)
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine femap_material(io,mark,version,mid,dmat)
c
c Write material data.
c Write a data set with mostly zero values to fill the needed
c data values to avoid reading errors in FEMAP.
c Note that there are 18 expected records if the FEMAP version
c number is 4.1, and 33 records if the version is 4.5
c
      implicit none
      integer   io,id,mid,color,lid,matype
      real      version
      double precision dmat(*),Emod,Gmod,poisson
      character        mark*5,title*78
c
      id     = 401    ! ID number for the data block
      color  = 55     ! "color"
      mid    = 1      ! Material ID number
      matype = 0      ! Material Type Isotropic
      lid    = 1      ! Layer ID
      title  = 'Isotrophic Linear Elastic Material'
      Emod    = dmat(1)
      poisson = dmat(2)
      Gmod    = Emod / (2.*(1.+poisson))
c
      write(io,'(a5)')  mark
      write(io,'(i3)')  id
      write(io,'(5i6)') mid,color,matype,lid, 0
      write(io,'(a78)')   title
      write(io,'(3g16.9)')Emod,Emod,Emod
      write(io,'(3g16.9)')Gmod,Gmod,Gmod
      write(io,'(3g16.9)')poisson,poisson,poisson
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(f14.2)')0.0
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(f14.2)')0.0
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(f14.2)')0.0
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(f14.2)')0.0
c
c specific heat, density, damping, ref. temperature
c  (arbitrary at the moment)
c
      write(io,'(4f14.4)')0.1,0.283,0.0,70.0
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      if (version.ge.4.5)then
         write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
         write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
         write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
         write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
         write(io,'(f14.2)')0.0
         write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
         write(io,'(3f14.2)')0.0,0.0,0.0
         write(io,'(3f14.2,2i8)')0.0,0.0,0.0,0,0
         write(io,'(5i8)')0,0,0,0,0
         write(io,'(6i8)')0,0,0,0,0,0
         write(io,'(3f14.2)')0.0,0.0,0.0
         write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
         write(io,'(2f14.2)')0.0,0.0
         write(io,'(4i8)')0,0,0,0
         write(io,'(4i8)')0,0,0,0
      end if
      write(io,'(a5)')mark
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      subroutine femap_elemproperty(io,mark,pid,mid,etypid)
c
c Write element property data
c
      implicit none
      integer   io,id,color,pid,mid,lid,etypid
      character mark*5,title*78
c
c   pid = Property ID ( as defined in for the element  )
c   mid = Material ID ( as defined in for the Material )
c
      id    = 402   ! ID number for the data block
      color = 110   ! "color"
      lid    = 1      ! Layer ID
      title = 'Basic Element Properties'
      write(io,'(a5)') mark
      write(io,'(i3)') id
c
c property ID, color, material ID (see above),
c      property type (same as the element type),
c          layer ID, reference coord. system
c
      write(io,'(6i6)') pid,color,mid,etypid,lid,0
      write(io,'(a78)')   title
      write(io,'(4i8)') 0,0,0,0
      write(io,'(i8)')  0
c
c maximum number of property values to follow
c
      write(io,'(i8)')42
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(5f14.2)')0.0,0.0,0.0,0.0,0.0
      write(io,'(2f14.2)')0.0,0.0
      write(io,'(a5)')mark
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
