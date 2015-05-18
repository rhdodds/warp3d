c ****************************************************************
c
c             8to20.f
c
c             for sharp crack-tip meshes modeled with collapsed
c             element faces:
c             1. convert 8-noded bricks to 20-noded bricks 
c             2. condense numbering of crack-front nodes
c             3. create 1/4-point nodes.
c
c ***************************************************************
c ***************************************************************
c
c             module: declaration of global variables
c
c ***************************************************************
c
      module main_data
c
      integer, parameter :: keyboard=5, screen=6, input_file=3,
     &                output_file=4, max_blocks=800, block_size=128
      integer :: num_elems, num_nodes, condense, max_nodes,
     &           cf_node_counter, last_node, num_blocks, alloc_stat
      integer, dimension(:,:), save :: elem_edge_nodes(3,12)
      integer, allocatable, dimension(:), save :: node_nums, elem_nums
      integer, allocatable, dimension(:,:), save :: incids, blocking,
     &                                              coincident_nodes
      integer, allocatable, dimension(:,:,:), save :: node_conn
c
      double precision, allocatable, dimension(:,:), save :: coords
c
      logical :: condense_nodes, qp_nodes
      logical, allocatable, dimension(:), save :: cf_node_flags
c
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
c
      end module
c
c
c******************************************************************
c
c             main program
c
c******************************************************************
      program convert_elems
      use main_data
      implicit none
c
c             read coordinates & incidences of 8-noded mesh in
c             warp3d format. read flags for condensing crack-
c             front nodes and generating 1/4-point elements.
c
      call read_data
c
c             locate crack-front nodes.they are on collapsed
c             element faces.
c
      call find_cf_nodes
c
c             if requested by user, condense node numbering of
c             crack-front nodes
c
      if( condense_nodes ) call condense_cf_nodes
c
c             insert node on the midside of every element edge.
c             if user has requested 1/4-point nodes, generate
c             them here.
c
      call create_midside_nodes
c
c             output coordinates and incidences arrays in warp3d
c             format. blocking output here is identical to blocking
c             read from original incidences file.
c
      call write_data(3)
c
c   call deallocate_arrays
c
c
      end
c
c****************************************************************
c
c             read user data and generate coordinates and
c             incidences arrays. note: this program reads
c             element blocking based on continuous blocks of
c             size 128.
c
c****************************************************************
      subroutine read_data
      use main_data
      implicit none
c
c             local variables
c
      integer :: ios, i, j, elemno, blocking_flag, condense_flag,
     &           qp_flag, elem_count, len_coordfile, len_incidfile
      double precision :: ratio
      character(len=80) :: coordfile, incidfile, scratch
      logical :: debug, file_exist
c
      debug = .false.
c
c             read coordinates, incidences, and user-defined input.
c
      num_elems     = 0
      num_nodes     = 0
      blocking_flag = 0
      num_blocks    = 0
      condense_flag = 0
      qp_flag       = 0
c
      write( screen,100,advance="no" ) "> number of elements: "
      read( keyboard,* ) num_elems
c
      write( screen,100,advance="no" ) "> number of nodes: "
      read( keyboard,* ) num_nodes
c
      file_exist = .false.
      do while ( .not.file_exist )
         write( screen,100,advance="no" ) "> name of coordinates file: "
         read( keyboard,100 ) coordfile
         coordfile     = ADJUSTL( coordfile )
         len_coordfile = LEN_TRIM( coordfile )
         inquire( exist=file_exist,file=coordfile )
         if( .not.file_exist ) then
            write( screen,200 ) coordfile(1:len_coordfile),
     &                          " does not exist."
            if( coordfile.eq."stop" ) stop
            cycle
         end if
         file_exist = .true.
      end do
c
      file_exist = .false.
      do while ( .not.file_exist )
         write( screen,100,advance="no" ) "> name of incidences file: "
         read( keyboard,100 ) incidfile
         incidfile     = ADJUSTL( incidfile )
         len_incidfile = LEN_TRIM( incidfile )
         inquire( exist=file_exist,file=incidfile )
         if( .not.file_exist ) then
            write( screen,200 ) incidfile(1:len_incidfile),
     &                          " does not exist."
            if( incidfile.eq."stop" ) stop
            cycle
         end if
         file_exist = .true.
      end do
c
      write( screen,100,advance="no" ) "> condense numbering of",
     & " coincident nodes? ( 1=yes, 0=no ): "
      read( keyboard,* ) condense_flag
      condense_nodes = .false.
      if( condense_flag.eq.1 ) condense_nodes = .true.
c     
      write( screen,100,advance="no" ) "> create 1/4-point nodes?",
     & " ( 1=yes, 0=no ): "
      read( keyboard,* ) qp_flag
      qp_nodes = .false.
      if( qp_flag.eq.1 ) qp_nodes = .true.
c
c             allocate and initialize coordinate & incidences arrays.
c
      max_nodes = num_elems * 20
      allocate( coords(max_nodes,3) )
      allocate( node_nums(max_nodes) )
      allocate( incids(num_elems,20) )
      allocate( elem_nums(num_elems) )
      allocate( blocking(max_blocks,3) )
      coords(1:max_nodes,1:3)    = 0.0
      node_nums(1:max_nodes)     = 0
      incids(1:num_elems,1:20)   = 0
      elem_nums(1:num_elems)     = 0
      blocking(1:max_blocks,1:3) = 0
c
c             read and store info from coordinates file.
c
      ios = 0  
      open( input_file,file=coordfile,status="old",iostat=ios )
      if( ios.ne.0 ) then
         write( screen,100 ) "ios open coordfile .ne. zero"
         stop
      end if
c
      read( input_file,100 ) scratch
      read( input_file,100 ) scratch
      read( input_file,100 ) scratch
      do i=1,num_nodes
         read( input_file,* ) node_nums(i),
     &                       (coords(node_nums(i),j),j=1,3)
      end do
      close( unit=input_file,status='keep' )
c
c             read & store info from incidences file.
c
      ios = 0
      open( input_file,file=incidfile,status="old",iostat=ios )
      if( ios.ne.0 ) then
         write( screen,100 ) "ios open incidfile .ne. zero"
         stop
      end if
c
      read( input_file,100 ) scratch
      read( input_file,100 ) scratch
      read( input_file,100 ) scratch
      do i=1,num_elems
         read( input_file,* ) elem_nums(i),
     &                       (incids(elem_nums(i),j),j=1,8)
      end do
      do while( scratch(1:5).ne."block" )
         read( input_file,100 ) scratch
         scratch = ADJUSTL( scratch )
      end do
      i=0
      do while( elem_count.lt.num_elems )
         i = i + 1
         read( input_file,* ) (blocking(i,j),j=1,3)
         elem_count = blocking(i,3) + blocking(i,2) - 1
      end do
      num_blocks = i
      close( unit=input_file,status='keep' )
c
c             write files for debugging
c
      if( debug ) call write_data(1)
c
      write( screen,300 ) "finished reading mesh data."
c
      return
c
 100  format(a)
 200  format(2a)
 300  format(/a)
c
      end subroutine read_data
c
c****************************************************************
c
c             find crack-front nodes by locating coincident nodes.
c
c****************************************************************
      subroutine find_cf_nodes
      use main_data
      implicit none
c
c             local variables
c
      integer :: node1, node2, row, ios, i, j, node_count, elem,
     &     counter
      double precision, parameter :: toler1=0.1e-3, toler=0.1e-6
      double precision :: x1, x2, x3, y1, y2, y3, z1, z2, z3,
     &     dx, dy, dz, dist, axis1, axis2
      logical :: debug
c
      debug = .true.
c
c             allocate arrays for tracking crack-front nodes.
c
      allocate( coincident_nodes(0:100,max_nodes) )
      allocate( cf_node_flags(1:max_nodes) )
      coincident_nodes(0:100,1:num_nodes) = 0
      cf_node_flags(1:max_nodes)          = .false.
c     
c             loop through coordinates. identify & record location
c             of coincident nodes. this will identify collapsed-face
c             crack-front nodes with unique node numbers.
c
      cf_node_counter = 0
      do node1=1,num_nodes-1
         do node2=node1+1,num_nodes
            x1   = coords(node1,1)
            y1   = coords(node1,2)
            z1   = coords(node1,3)
            x2   = coords(node2,1)
            y2   = coords(node2,2)
            z2   = coords(node2,3)
            dx   = abs(x1 - x2)
            dy   = abs(y1 - y2)
            dz   = abs(z1 - z2)
            dist = (dx**2 + dy**2 + dz**2)**0.5
c
c             insert lines here to omit nodes on an axis of symmetry
c             or on closed crack faces
c
c         axis1 = (x1**2 + z1**2)**0.5
c         axis2 = (x2**2 + z2**2)**0.5
c         if( axis1.le.toler1.or.axis2.le.toler1 ) cycle
c
            if( dist.le.toler ) then
               if( .not.cf_node_flags(node2) ) then
                  cf_node_flags(node1)        = .true.
                  cf_node_flags(node2)        = .true.
                  cf_node_counter             = cf_node_counter + 1
                  coincident_nodes(0,node1)= coincident_nodes(0,node1)+1
                  row                      = coincident_nodes(0,node1)
                  coincident_nodes(row,node1) = node2
               end if
            end if
         end do
      end do
c
c             loop through incidences. identify node numbers that are
c             listed twice on the same element. this indentifies
c             collapsed-face crack-front nodes with repeated node
c             numbering.
c
      elemloop: do elem=1,num_elems
      do i=1,7
         node1 = incids(elem,i)
         do j=i+1,8
            node2 = incids(elem,j)
            if( node1.eq.node2 ) then
c
c             insert lines here to omit nodes on an axis of symmetry
c             or on closed crack faces
c
c               x1 = coords(node1,1)
c               y1 = coords(node1,2)
c               z1 = coords(node1,3)
c               x2 = coords(node2,1)
c               y2 = coords(node2,2)
c               z2 = coords(node2,3)
c               axis1 = (x1**2 + z1**2)**0.5
c               axis2 = (x2**2 + z2**2)**0.5
c               if( axis1.le.toler1.or.axis2.le.toler1 ) cycle
c
               cf_node_flags(node1) = .true.
            end if
         end do
      end do
      end do elemloop
      counter = COUNT( cf_node_flags.eq..true. )
c
      if( debug ) call write_data(2)
c
      write( screen,100 ) "finished search for coincident nodes."
      write( screen,90 ) "total coincident nodes found: ",counter   
      write( screen,90 ) "crack-front locations: ",
     &                    counter - cf_node_counter
      write( screen,90 ) "duplicate coincident nodes: ",cf_node_counter
c
      return
c
 90   format(a,i6)
 100  format(/a)
c
      end subroutine find_cf_nodes
c
c****************************************************************
c
c             renumber crack-front nodes so all coincident nodes
c             will share one node number.
c
c****************************************************************
      subroutine condense_cf_nodes
      use main_data
      implicit none
c
c             local variables
c
      integer :: i, j, node, node1, node2, node_count,
     &           delete_node, elem, ios, deleted_node_count, row
      integer, allocatable, dimension(:) :: temp_node_nums
      integer, allocatable, dimension(:,:) :: temp_incids
      double precision, allocatable, dimension(:,:) :: temp_coords
      logical :: debug
c
      debug = .false.
c
c             allocate space to store temporary arrays
c
      allocate( temp_incids(num_elems,20) )
      allocate( temp_coords(max_nodes,3) )
      allocate( temp_node_nums(max_nodes) )
      temp_incids(1:num_elems,1:20) = incids(1:num_elems,1:20)
      temp_coords(1:num_nodes,1:3)  = 0.0
      temp_node_nums(1:max_nodes)   = 0
c
c             renumber nodes to account for removal of repeated nodes.
c             the original node number = the corresponding index of
c             'node_nums.' the new node number will be the number
c             actually stored in 'node_nums.'
c
      deleted_node_count = 0
      do node=1,num_nodes
         node_count = coincident_nodes(0,node)
         if( node_count.ne.0 ) then
            do i=1,node_count
               delete_node                = coincident_nodes(i,node)
               deleted_node_count         = deleted_node_count + 1
               cf_node_flags(delete_node) = .false.
               do j=1,num_nodes
                  if( node_nums(j) .le. node_nums(delete_node)) cycle
                  if( j .eq. delete_node ) cycle
                  node_nums(j) = node_nums(j) - 1
               end do
               node_nums(delete_node) = 0            
            end do
         end if
      end do
c
c             replace entries of coincident nodes listed in 'incids'
c             with lowest-numbered coincident node.
c
      do node=1,num_nodes
         node_count = coincident_nodes(0,node)
         if( node_count.ne.0 ) then
            do i=1,node_count
               node1                    = coincident_nodes(i,node)
               coincident_nodes(i,node) = 0
               coincident_nodes(0,node) = coincident_nodes(0,node) - 1
               do elem=1,num_elems
                  do j=1,8
                     node2 = incids(elem,j)
                     if( node1.eq.node2 ) temp_incids(elem,j) = node
                  end do
               end do
            end do
         end if
      end do
      incids(1:num_elems,1:20) = temp_incids(1:num_elems,1:20)
c
c             renumber entries in incids according to data in 'node_nums.'
c
      do elem=1,num_elems
         do i=1,8
            node           = incids(elem,i)
            incids(elem,i) = node_nums(node)
         end do
      end do
c
c             update cf_node_flags to reflect changes in node numbering.
c
      do i=1,num_nodes
         row = node_nums(i)
         if( row.eq.0 .or. row.eq.i ) cycle
         if( i.gt.max_nodes .or. row.gt.max_nodes .or. row.lt.0 ) then
            write(*,*) "error in node numbering"
            stop
         end if
         if( cf_node_flags(i) .eq. .true. ) then
            cf_node_flags(i)   = .false.
            cf_node_flags(row) = .true.
         end if
      end do
c
c             create new coordinates array.
c
      do i=1,num_nodes
         row = node_nums(i)
         if( row.eq.0 ) cycle
         temp_coords(row,1:3) = coords(i,1:3)
      end do
      coords(1:max_nodes,1:3) = temp_coords(1:max_nodes,1:3)
      num_nodes               = num_nodes - deleted_node_count
c
      if( debug ) then
         call write_data(2)
         call write_data(1)
      end if
c
      deallocate( temp_incids )
      deallocate( temp_coords )
      deallocate( temp_node_nums )
c
      write( screen,100 )"finished condensing crack-front nodes."
      write( screen,90 )"total number of nodes deleted:",
     &                   deleted_node_count
      write( screen,90 )"total remaining coincident nodes",
     &                  cf_node_counter - deleted_node_count
      write( screen,90 )"total nodes remaining in model:",num_nodes
c     
      return
c
 90   format(a,i6)
 100  format(/a)
c
      end subroutine condense_cf_nodes
c
c****************************************************************
c
c             create mid-side nodes.
c
c****************************************************************
      subroutine create_midside_nodes
      use main_data
      implicit none
c
c             local variables
c
      integer :: next_node, elem, pair, node1, node2, node_count,
     &           i, node, row1, row2, mid_node
      double precision :: x1, y1, z1, x2, y2, z2, xmid, ymid, zmid,
     &                    xq, yq, zq
      logical :: debug, qp_edge, node1_cf, node2_cf
c     
      debug = .false.
c
c             allocate array 'node_conn' for edge-node information.
c             node_conn(0,i,1) contains the number of nodes that
c             share an edge with node i. node_conn(1:,i,1) lists
c             the node numbers that share an edge with node i.
c             node_conn(0,:,2) is not used. node_conn(1:,i,2) lists
c             node numbers of mid-side nodes on the edge connecting
c             node i and node(1:,i,1)
c
      allocate( node_conn(0:50,max_nodes,2) )
      node_conn(0:50,0:max_nodes,1:2) = 0
c
c             loop through element incidences. generate node_conn.
c             create midside nodes where they doesn't exist. place
c             mid-side node number in incids and node_conn. place
c             mid-side node coordinates in coords.
c
      next_node = num_nodes + 1
      elemloop: do elem=1,num_elems
      pairloop: do pair=1,12
c
c             obtain nodes on next element edge.
c
      node1 = incids(elem,elem_edge_nodes(1,pair))
      node2 = incids(elem,elem_edge_nodes(2,pair))
      if( ANY( node_conn(1:,node1,1) .eq. node2 ) ) then
c
c             midside node already exists. retrieve and place in incids.
c
         node_count = node_conn(0,node1,1)
         do i=1,node_count
            node = node_conn(i,node1,1)
            if( node.eq.node2 ) then
               row1 = i
               exit
            end if
         end do
         mid_node = node_conn(row1,node1,2)
         incids(elem,elem_edge_nodes(3,pair)) = mid_node
         cycle pairloop
      end if
c
c             midside node doesn't exist. assign it a number.
c
      if( node1.ne.node2 ) then
         mid_node  = next_node
         num_nodes = mid_node
         next_node = next_node + 1
      else
         mid_node  = node1
      end if
c
c             update incids
c
      incids(elem,elem_edge_nodes(3,pair)) = mid_node
c
c             update coords. if user has requested generation
c             of quarter-point nodes, generate them here. a quarter-
c             point node is identified as follows: if node1 and node2
c             are coincident with each other, the mid_node is on a
c             collapsed edge of an element. if node1 and node2 are
c             each coincident with other nodes (and not each other),
c             mid_node is on a crack front edge. therefore, if only
c             one of node1 and node2 is coincident with other nodes,
c             a quarter-point node may be generated between them.
c
      x1   = coords(node1,1)
      y1   = coords(node1,2)
      z1   = coords(node1,3)
      x2   = coords(node2,1)
      y2   = coords(node2,2)
      z2   = coords(node2,3)
      xmid = (x1 + x2)/2
      ymid = (y1 + y2)/2
      zmid = (z1 + z2)/2
      if( qp_nodes ) then
         qp_edge  = .false.
         node1_cf = .false.
         node2_cf = .false.
         if( cf_node_flags(node1) ) node1_cf = .true.
         if( cf_node_flags(node2) ) node2_cf = .true.
         if( node1_cf.or.node2_cf ) qp_edge  = .true.
         if( node1_cf.and.node2_cf ) qp_edge = .false.
         if( qp_edge ) then
            if( node1_cf ) then
               xq = (x1 + xmid)/2
               yq = (y1 + ymid)/2
               zq = (z1 + zmid)/2
            end if
            if( node2_cf ) then
               xq = (x2 + xmid)/2
               yq = (y2 + ymid)/2
               zq = (z2 + zmid)/2
            end if
            xmid = xq
            ymid = yq
            zmid = zq
         end if
      end if
      coords(mid_node,1) = xmid
      coords(mid_node,2) = ymid
      coords(mid_node,3) = zmid
c
c             update node_conn
c
      node_conn(0,node1,1)    = node_conn(0,node1,1) + 1
      node_conn(0,node2,1)    = node_conn(0,node2,1) + 1
      row1                    = node_conn(0,node1,1)
      row2                    = node_conn(0,node2,1)         
      node_conn(row1,node1,1) = node2
      node_conn(row2,node2,1) = node1
      node_conn(row1,node1,2) = mid_node
      node_conn(row2,node2,2) = mid_node
c
      end do pairloop
      end do elemloop
c
      if( debug ) call write_data(1)
c
      write( screen,100 ) "finished creating midside nodes."
      write( screen,90 ) "new total nodes in model:",num_nodes
c
      return
c
 90   format(a,i6)
 100  format(/a)
c
      end subroutine create_midside_nodes
c
c****************************************************************
c
c             write coordinates and incidences in warp3d format.
c             blocking must be appended to incidences.
c
c****************************************************************
      subroutine write_data(flag)
      use main_data
      implicit none
c     
c             dummy variable
c
      integer, intent(in) :: flag
c
c             local variables
c
      integer :: ios, i, j, node_count
      character(len=80) :: coord, incid, stat, posit
c
c             output will be for debugging purposes or
c             for writing new warp3d files.
c
      select case(flag)
      case(1)
         coord = "debugfile"
         incid = "debugfile"
         stat  = "unknown"
         posit = "append"
      case(2)
         coord = "debugfile"
         incid = "debugfile"
         stat  = "unknown"
         posit = "append"
         goto 1000
      case(3)
         coord = "new.coordinates"
         incid = "new.incid_and_blocking"
         stat  = "replace"
         posit = "rewind"
      end select
c
c             write coordinates.
c
      ios = 0
      open(output_file,file=coord,status=stat,position=posit,iostat=ios)
      if( ios.ne.0 ) then
         write( screen,100 ) "ios output coord .ne. zero"
         stop
      end if
      write( output_file,200 )
      do i=1,num_nodes
         write( output_file,300 ) i,(coords(i,j),j=1,3)      
      end do
      write( output_file,400 )
      close( unit=output_file,status='keep' )
c
c             write incidences.
c
      ios = 0
      open(output_file,file=incid,status=stat,position=posit,iostat=ios)
      if( ios.ne.0 ) then
         write( screen,100 ) "ios output incid .ne. zero"
         stop
      end if
      write( output_file,500 )
      do i=1,num_elems
         write(output_file,600) elem_nums(i), (incids(i,j),j=1,8)
         write(output_file,600) (incids(i,j),j=9,17)
         write(output_file,700) (incids(i,j),j=18,20)
      end do
      write( output_file,400 )
      write( output_file,800 )
      do i=1,num_blocks
         write( output_file,700 ) (blocking(i,j),j=1,3)
      end do
      write( output_file,900 )
      close( unit=output_file,status='keep' )
c
      return
c
c             write crack-front node data
c
 1000 continue
      ios = 0
      open( output_file,file="debugfile",status="unknown",
     &      position="append",iostat=ios )
      if( ios.ne.0 ) then
         write( screen,100 ) "open failed in find_cf_nodes"
         stop
      end if
      write( output_file,100 ) "all crack-front nodes"
      do i=1,num_nodes
         if( cf_node_flags(i) ) write( output_file,950 ) i
      end do
      write( output_file,100 ) "crack-front node groups"
      do i=1,num_nodes
         if( cf_node_flags(i) ) then
            write( output_file,950,advance="no" ) i
            node_count = coincident_nodes(0,i)
            if( node_count.ne.0 ) then
               do j=1,node_count
                  write( output_file,960,advance="no" )
     &                  coincident_nodes(j,i)
               end do
               write( output_file,970 )
            else
               write( output_file,970 )
            end if
         end if
      end do
      close( unit=output_file,status='keep' )
c
      return
c
 100  format(/a)
 200  format('c',/,'coordinates',/,'*echo off')
 300  format(i7,3(2x,e16.9))
 400  format('*echo on',/,'c')
 500  format('c',/,'incidences',/,'*echo off')
 600  format(9i8,',')
 700  format(3i8)
 800  format('c',/,'blocking',5x,'$  scalar')
 900  format('c',/,'c')
 950  format(i7)
 960  format(',',1x,i6)
 970  format(/)
c
      end subroutine write_data
c
c****************************************************************
c
c             deallocate arrays
c
c****************************************************************
      subroutine deallocate_arrays
      use main_data
      implicit none
c
      if( allocated(coords) ) deallocate( coords )
      if( allocated(node_nums) ) deallocate( node_nums )
      if( allocated(incids) ) deallocate( incids )
      if( allocated(elem_nums) ) deallocate( elem_nums )
      if( allocated(blocking) ) deallocate( blocking )
      if( allocated(coincident_nodes) ) deallocate( coincident_nodes )
      if( allocated(cf_node_flags) ) deallocate( cf_node_flags )
c
      return
c
 100  format(/a)
c
      end subroutine deallocate_arrays
c


