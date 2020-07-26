c
c2345678901234567890
c              given a mesh of 20-node elements. make all mid-side
c              nodes lie on straight edges between the corner nodes
c
      implicit none
c
      integer :: numnodes, numelems, i, node, elem, edge, ic1, ic2,
     &           im, snodec1, snodec2, snodem
      integer :: inc(20),  elem_edge_nodes(3,12)
      integer, allocatable :: incid(:,:)
   
!123456
      double precision :: x, y, z, cnew(3), half
      double precision, allocatable :: coords_old(:,:),
     &                                 coords_new(:,:)

      data elem_edge_nodes / 1, 2, 9,
     &                       2, 3, 10,
     &                       3, 4, 11,
     &                       4, 1, 12,
     &                       5, 6, 13,
     &                       6, 7, 14,
     &                       7, 8, 15,
     &                       8, 5, 16,
     &                       1, 5, 17,
     &                       2, 6, 18,
     &                       3, 7, 19,
     &                       4, 8, 20 /


      numnodes = 105905
      numelems = 22930
      allocate( incid(20,numelems), coords_old(3,numnodes),
     &          coords_new(3,numnodes) )

      half = 0.5d0
c
c234567890123456
c              coordinates:  coords.inp
c              incidences:   incid.inp
c              both must be bare files with no other lines
c              new coordinates written: new_coords.inp
c

      open(unit=10,file='coords.inp')
      open(unit=20,file='incid.inp')
c
      do i = 1, numnodes
        read(10,*) node, x, y, z
        coords_old(1,node) = x
        coords_old(2,node) = y
        coords_old(3,node) = z
      end do
      write(*,*) '...  coords read'
      do i = 1, numelems
        read(20,*) elem, inc(1:20)
        incid(1:20,elem) = inc
      end do
      write(*,*) '...  incidences read'
      close(10)
      close(20)
c
      coords_new = coords_old
c
      do elem = 1, numelems
        do edge = 1, 12
          ic1 = elem_edge_nodes(1,edge)
          ic2 = elem_edge_nodes(2,edge)
          im =  elem_edge_nodes(3,edge)
          snodec1 = incid(ic1,elem)
          snodec2 = incid(ic2,elem)
          snodem  = incid(im,elem)
          cnew(1:3) = half * ( coords_old(1:3,snodec1) + 
     &                         coords_old(1:3,snodec2) )
          coords_new(1:3,snodem) = cnew(1:3)
       end do  ! element edge number
      end do   ! elements
      write(*,*) '...  mid-edge nodes adjusted'
c
      open(unit=10,file='new_coords.inp')
      write(10,*) "*echo off"
      write(10,*) "!"
      write(10,*) " coordinates"
      do node = 1, numnodes
        x = coords_new(1,node)
        y = coords_new(2,node)
        z = coords_new(3,node)
        write(10,9000) node, x, y, z
      end do
      write(*,*) '...  new coords written'
c           
      stop
 9000 format(i8,3e17.9)
      end

  


