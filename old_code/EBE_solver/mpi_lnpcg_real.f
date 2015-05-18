c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_init_owner              *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/04/98                   *
c     *                                                              *
c     *     this subroutine calculates the node owner structure      *
c     *     for domain decomposition portions of the MPI code        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_init_owner
      use main_data, only: elems_to_blocks, inverse_incidences
      use mpi_lnpcg
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      character *80 line
      integer, allocatable, dimension(:,:) ::  node2proc
      dimension num_private_nodes(0:max_procs-1), 
     &     num_own_shared_nodes(0:max_procs-1), 
     &     num_shared_nodes(0:max_procs-1),
     &     num_local_nodes(0:max_procs-1),
     &     count_priv(0:max_procs-1), count_own(0:max_procs-1), 
     &     count_shar(0:max_procs-1), 
     &     conn(max_procs,0:max_procs-1), conn_num(0:max_procs-1),
     &     ordering(max_procs,max_procs), order_num(max_procs)
      logical debug, been_here, blk_on_boundary(mxnmbl)
      real wcputime
      external wcputime
      data debug, been_here /.false., .false./
      save been_here
c
c
      if ( debug ) write (out,*) myid,': ->> in wmpi_init_owner'
c 
c            if we have executed this routine previously,
c            then skip it.
c
      if ( been_here ) return
      been_here = .true.
c
c            init variables
c
      num_private_nodes(0:numprocs-1) = 0 
      num_own_shared_nodes(0:numprocs-1) = 0 
      num_shared_nodes(0:numprocs-1) = 0
      num_local_nodes(0:numprocs-1) = 0
      count_priv(0:numprocs-1) = 0 
      count_own(0:numprocs-1) = 0 
      count_shar(0:numprocs-1) = 0
c
c            allocate structure which, when given a specific node, 
c            tells what processors access it.
c
      allocate (node2proc(0:numprocs,nonode), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
      node2proc(0:numprocs,1:nonode) = 0
c
c            set up data structures needed to calculate
c            the owner data structures
c
      call wmpi_owner_setup ( node2proc, num_private_nodes, 
     &     num_own_shared_nodes, num_shared_nodes, num_local_nodes,
     &     blk_on_boundary)
c
      if (debug) write (out,*) '...allocating'
c
c            allocate space on root proc for all the info. We will
c            deallocate this space once all processors have been
c            sent their information.
c
      if ( allocated(proc_nodes)) deallocate (proc_nodes)
      allocate (proc_nodes(0:numprocs-1), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
      do proc = 0, numprocs - 1
c
c              allocate data structure for each processor
c
         call allocate_owner(num_private_nodes(proc), 
     &        num_own_shared_nodes(proc), num_shared_nodes(proc),
     &        nelblk, proc_nodes(proc))
c
c              allocate vectors to hold lists of shared nodes for
c              each processor.  Here we allocate these vectors to
c              be larger than we need; since this is a temporary
c              data structure, this is not a problem.
c
         do i = 0, numprocs -1
c
            allocate (proc_nodes(proc)%sharing_procs(i)%ptr(
     &           num_own_shared_nodes(proc)), stat = alloc_stat)
            if ( alloc_stat .ne. 0 ) then
                write(out,9900)
               call die_abort
               stop
            end if
c
            allocate (proc_nodes(proc)%shared_owner(i)%ptr(
     &           num_shared_nodes(proc)), stat = alloc_stat)
            if ( alloc_stat .ne. 0 ) then
                write(out,9900)
               call die_abort
               stop
            end if
c
         enddo
c
      enddo
c
c            now fill the data structures
c     
      if (debug) write (out,*) '...filling structure'
      do node = 1, nonode
         own_proc = node2proc(1,node)
c
c               if there is only one processor referencing this node, then
c               this node is private to that processor.  increment the number
c               of private nodes by one for that processor
c
         if (node2proc(0,node).eq. 1) then
            count_priv(own_proc) = count_priv(own_proc) + 1
            proc_nodes(own_proc)%private(count_priv(own_proc)) = node
         else
c
c               there are more than 1 processors which reference the node.
c               for the owning processor, increase the number of own_shared
c               nodes by one. Also, for all other processors that reference
c               the node but don't own it, increase their number of shared
c               nodes.
c
            count_own(own_proc) = count_own(own_proc) + 1
            proc_nodes(own_proc)%own_shared(count_own(own_proc)) = node
            do i= 2, node2proc(0,node)
               new_proc = node2proc(i,node)
c
c                   put node in owner's list of nodes shared with
c                   sharing processor new_proc
c
               proc_nodes(own_proc)%sharing_count(new_proc) = 
     &              proc_nodes(own_proc)%sharing_count(new_proc) +1
               proc_nodes(own_proc)%sharing_procs(new_proc)%ptr(
     &              proc_nodes(own_proc)%sharing_count(new_proc)) =
     &              count_own(own_proc)
c
c                   put node in sharer's list of nodes shared with
c                   the owning processor own_proc.
c     
               count_shar(new_proc) = count_shar(new_proc) + 1
               proc_nodes(new_proc)%shared(count_shar(new_proc)) = node
               proc_nodes(new_proc)%shared_owner(own_proc)%ptr(
     &              proc_nodes(own_proc)%sharing_count(new_proc)) =
     &              count_shar(new_proc)
               proc_nodes(new_proc)%shared_count(own_proc) =
     &              proc_nodes(own_proc)%sharing_count(new_proc) 
c     
            enddo
         endif
      enddo
c
c            for ebe preconditioner, we want to know which
c            blocks are completely internal (have no shared
c            nodes).  
c
c               First find out which blocks
c               have elements on the boundary between two
c               processors, and which do not.
c
      do blk = 1, nelblk
         blk_on_boundary(blk) = .false.
      enddo
c
      do node = 1, nonode
         if ( node2proc(0,node) .eq. 1 ) cycle
         do j = 1, inverse_incidences(node)%element_count
            elem = inverse_incidences(node)%element_list(j)
            blk  = elems_to_blocks(elem, 1)
            blk_on_boundary(blk) = .true.
         enddo
      enddo
c
c               fill list of internal blocks for each processor
c
      do proc = 0, numprocs - 1 
         proc_nodes(proc)%num_int_blks = 0 
      enddo
      do blk = 1, nelblk
         if ( blk_on_boundary(blk) ) cycle
         proc = elblks(2,blk)
         proc_nodes(proc)%num_int_blks = 
     &        proc_nodes(proc)%num_int_blks + 1
         proc_nodes(proc)%internal_blks(proc_nodes(proc)%num_int_blks) =
     &        blk
      enddo
c
c            create local dof to global dof pointer array -- note this
c            assumes each node has three degrees of freedom.
c            also create global dof to local dof pointer array, pointing
c            the opposite direction.
c
c               loop over the processors
c
      do proc = 0, numprocs-1
c
c                  process the private nodes of the current processor
c
         do i=1, proc_nodes(proc)%num_private
c
            proc_nodes(proc)%local2global((i-1)*3 + 1) = 
     &           dstmap(proc_nodes(proc)%private(i))
            proc_nodes(proc)%local2global((i-1)*3 + 2) = 
     &           dstmap(proc_nodes(proc)%private(i))+1
            proc_nodes(proc)%local2global(i*3) = 
     &           dstmap(proc_nodes(proc)%private(i))+2
c
            proc_nodes(proc)%global2local(dstmap(
     &           proc_nodes(proc)%private(i))) = (i-1)*3 + 1
            proc_nodes(proc)%global2local( dstmap(
     &           proc_nodes(proc)%private(i)) + 1 ) = (i-1)*3 + 2
            proc_nodes(proc)%global2local(dstmap(
     &           proc_nodes(proc)%private(i)) + 2 ) = i*3 
c
         enddo
c
c                  process the nodes shared and owned by the current processor
c
         start = proc_nodes(proc)%num_private * 3
         do i=1, proc_nodes(proc)%num_own_shared
c
            proc_nodes(proc)%local2global((i-1)*3 + 1 + start)= 
     &           dstmap(proc_nodes(proc)%own_shared(i))
            proc_nodes(proc)%local2global((i-1)*3 + 2 + start)= 
     &           dstmap(proc_nodes(proc)%own_shared(i)) + 1
            proc_nodes(proc)%local2global(i*3 + start)= 
     &           dstmap(proc_nodes(proc)%own_shared(i)) + 2
c
            proc_nodes(proc)%global2local(dstmap(
     &           proc_nodes(proc)%own_shared(i))) = (i-1)*3+1+start
            proc_nodes(proc)%global2local(dstmap(
     &           proc_nodes(proc)%own_shared(i))+1) = (i-1)*3+2+start
            proc_nodes(proc)%global2local(dstmap(
     &           proc_nodes(proc)%own_shared(i))+2) = i*3+start
c
         enddo
c
c                  process the nodes shared but not owned by the processor
c
         start = ( proc_nodes(proc)%num_private +
     &        proc_nodes(proc)%num_own_shared ) * 3
         do i=1, proc_nodes(proc)%num_shared
c
            proc_nodes(proc)%local2global((i-1)*3 + 1 + start)= 
     &           dstmap(proc_nodes(proc)%shared(i))
            proc_nodes(proc)%local2global((i-1)*3 + 2 + start)= 
     &           dstmap(proc_nodes(proc)%shared(i)) + 1
            proc_nodes(proc)%local2global(i*3 + start) =
     &           dstmap(proc_nodes(proc)%shared(i)) + 2
c
            proc_nodes(proc)%global2local(dstmap(
     &           proc_nodes(proc)%shared(i))) = (i-1)*3 + 1 + start
            proc_nodes(proc)%global2local(dstmap(
     &           proc_nodes(proc)%shared(i))+1) = (i-1)*3 + 2 + start
            proc_nodes(proc)%global2local(dstmap(
     &           proc_nodes(proc)%shared(i))+2) = i*3 + start
c
         enddo
      enddo
c
c            print out the count for private and shared nodes, and the
c            number shared with each other processor.
c
      write (out,'("  -------------------------")')
      write (out,'("  Node Ownership Statistics")')
      write (out,'("  -------------------------")')
      do i=0, numprocs-1
         call print_node_own (proc_nodes(i),i)
      enddo
      write (out,'("  -------------------------")')
c
c            fill processor zero's own local_nodes structure -- this
c            is the permanent data stucture for the node ownership
c            data for the root processor
c
      if (debug) write (out,*) '...do proc 0s local_nodes'
      call allocate_owner(num_private_nodes(0), 
     &     num_own_shared_nodes(0), num_shared_nodes(0),
     &     proc_nodes(0)%num_int_blks, local_nodes)
c     
      call copy_own (proc_nodes(0), local_nodes)
c
c            for the root processor, store all the mappings from
c            processor local to global degrees of freedom
c
c               first allocate the base data structure
c
      allocate (procdof2glob(numprocs-1), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
c               now allocate and fill the mappings for each slave proc
c
      do proc = 1, numprocs - 1
	 procdof2glob(proc)%num_dof = num_local_nodes(proc)*mxndof
	 allocate (procdof2glob(proc)%dof(num_local_nodes(proc)*mxndof),
     &       stat = alloc_stat)
         if ( alloc_stat .ne. 0 ) then
            write(out,9900)
            call die_abort
            stop
         end if
	 do i = 1, procdof2glob(proc)%num_dof 
            procdof2glob(proc)%dof(i) = proc_nodes(proc)%local2global(i)
         enddo
      enddo
c
c            for each processor, have the processor allocate
c            the appropriate data structure, then have the
c            root processor send the correct information to 
c            the processor.
c
      if (debug) write (out,*) '...doing multi-processor thang'
      call wmpi_owner_send(num_private_nodes, 
     &     num_own_shared_nodes, num_shared_nodes, num_local_nodes)
c
c
c            do splitting thing for parallel ebe.
c     
      call find_blk_groups ( node2proc )
c
c            deallocate all the structures that are no longer needed
c
      if (debug) write (out,*) '...deallocate'
      deallocate (node2proc)
c
      do proc = 0, numprocs - 1  
         do i = 0, numprocs - 1
            deallocate(proc_nodes(proc)%sharing_procs(i)%ptr)
            deallocate(proc_nodes(proc)%shared_owner(i)%ptr)
         enddo
      enddo
c
      do proc = 0, numprocs - 1  
         call deallocate_owner (proc_nodes(proc))
      enddo
      deallocate(proc_nodes)
c
      write (out,'(" >> Time required for input (sec):",f10.4)') 
     &    wcputime(1)
c
      return
 9900 format('>>> FATAL ERROR: memory allocate failure...')
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_owner_setup             *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/04/98                   *
c     *                                                              *
c     *     this subroutine calculates the node owner structure      *
c     *     for domain decomposition portions of the MPI code        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_owner_setup (node2proc, num_private_nodes, 
     &     num_own_shared_nodes, num_shared_nodes, num_local_nodes,
     &     blk_on_boundary)
c
      use main_data, only: inverse_incidences
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      character *80 line
      dimension num_private_nodes(0:max_procs-1), 
     &     num_own_shared_nodes(0:max_procs-1), 
     &     num_shared_nodes(0:max_procs-1),
     &     num_local_nodes(0:max_procs-1),
     &     elem2blk(mxel), node2proc(0:numprocs,1:nonode)
      logical its_new_proc, debug, blk_on_boundary(*)
      data debug /.false./
c
c                       build element to block list
c
      do blk = 1, nelblk
         blk_on_boundary(blk) = .false.
         felem          = elblks(1,blk)
         span           = elblks(0,blk)
         do elem = felem, felem + span - 1
            elem2blk(elem) = blk
         enddo
      enddo
c
c                       build structure which, when given a specific node, 
c                       tells what processors access it.
c
      do node = 1, nonode
         do elem_loop = 1, inverse_incidences(node)%element_count
            elem = inverse_incidences(node)%element_list(elem_loop)
            blk = elem2blk(elem)
            proc = elblks(2,blk) 
            if (node2proc(0,node) .eq. 0) then
               node2proc(0,node) = node2proc(0,node) + 1
               node2proc(1,node) = proc 
               cycle
            else
               its_new_proc = .true.
               do i=1, node2proc(0,node)
                  if ( node2proc(i,node) .eq. proc) its_new_proc=.false.
               enddo
               if (.not. its_new_proc) cycle
               node2proc(0,node) = node2proc(0,node) + 1
               node2proc(node2proc(0,node),node) = proc
            endif
c
c                           if node is even, make lowest number processor
c                           its owner.  The owner is the first processor
c                           in the list.  Thus make sure lowest number
c                           processor is first in list
c
c
            if ( mod(node,2) .eq. 0) then
               if (node2proc(node2proc(0,node),node) .lt. 
     &              node2proc(1,node)) then
                  tmp = node2proc(node2proc(0,node),node) 
                  node2proc(node2proc(0,node),node) = node2proc(1,node)
                  node2proc(1,node) = tmp
               endif
            else
c
c                           if node is odd, make highest number processor
c                           its owner.  
c
               if (node2proc(node2proc(0,node),node) .gt. 
     &              node2proc(1,node)) then
                  tmp = node2proc(node2proc(0,node),node) 
                  node2proc(node2proc(0,node),node) = node2proc(1,node)
                  node2proc(1,node) = tmp
               endif
            endif
c
         enddo
      enddo
c
c                       write out structure for checking
c
      if (debug) then
         write (out,*) '>>> proc accesses of nodes:'
         do i=1, nonode
            write (out,'(4x,i7,3x,10i3)') i, (node2proc(j,i),j=1,
     &           node2proc(0,i) )
         enddo
      endif
c
c                       count number of nodes for each processor of
c                       the following types:
c                        1:  private -- only one processor references the node
c                        2:  own_shared -- owned by processor, but shared with
c                                        others
c                        3:  shared -- owned by other processor, but referenced
c                                        by current processor
c         
      if (debug) write (out,*) '...counting'
      do node = 1, nonode
         proc = node2proc(1,node)
c
c             if there is only one processor referencing this node, then
c             this node is private to that processor.  increment the number
c             of private nodes by one for that processor
c
         if (node2proc(0,node).eq. 1) then
            num_private_nodes(proc) = num_private_nodes(proc) + 1
         else
c
c             there are more than 1 processors which reference the node.
c             for the owning processor, increase the number of own_shared
c             nodes by one. Also, for all other processors that reference
c             the node but don't own it, increase their number of shared
c             nodes.
c
            num_own_shared_nodes(proc) = num_own_shared_nodes(proc) + 1
            do i= 2, node2proc(0,node)
               num_shared_nodes(node2proc(i,node)) = 
     &              num_shared_nodes(node2proc(i,node)) + 1
            enddo
         endif
      enddo
c
      do proc = 0, numprocs - 1
         num_local_nodes(proc) =  num_private_nodes(proc) +  
     &        num_own_shared_nodes(proc) +  num_shared_nodes(proc)
      enddo
c
c
c     
      if (debug) then
         write (out,*) ' >> proc:   private  own_shared   shared'
         do proc = 0, numprocs - 1
            write (out,'(2x,i5,3(3x,i7))')
     &           proc, num_private_nodes(proc), 
     &           num_own_shared_nodes(proc), num_shared_nodes(proc)
         enddo
      endif
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_owner_send              *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/04/98                   *
c     *                                                              *
c     *     this subroutine sneds the node owner structures to       *
c     *     the slave processors                                     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_owner_send(num_private_nodes, 
     &     num_own_shared_nodes, num_shared_nodes, num_local_nodes)
      use mpi_lnpcg, only: local_nodes, proc_nodes
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      character *80 line
      dimension num_private_nodes(0:max_procs-1), 
     &     num_own_shared_nodes(0:max_procs-1), 
     &     num_shared_nodes(0:max_procs-1),
     &     num_local_nodes(0:max_procs-1),
     &     disp(mxdof), blksz(mxdof)
      logical debug
      data debug /.false./
c
      dimension status(MPI_STATUS_SIZE)
c
c
      if ( debug ) write (out,*) myid,': ->> in wmpi_owner_send'
c 
c            notify slaves
c
      call wmpi_alert_slaves ( 16 )
c
c            send the data
c
      if ( root_processor ) then
c
         do proc = 1, numprocs - 1
c     
c              send out the constants first
c
            call MPI_SEND (num_private_nodes(proc),1,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
            call MPI_SEND (num_own_shared_nodes(proc),1,MPI_INTEGER,
     &           proc,proc,MPI_COMM_WORLD,ierr)
            call MPI_SEND (num_shared_nodes(proc),1,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
            call MPI_SEND (proc_nodes(proc)%num_int_blks,1,MPI_INTEGER,
     &           proc,proc,MPI_COMM_WORLD,ierr)
c     
c              now send out the arrays
c
c                 send out information about local to/from global mapping
c
            call MPI_SEND (proc_nodes(proc)%local2global,
     &           num_local_nodes(proc)*mxndof,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
            call MPI_SEND (proc_nodes(proc)%global2local,
     &           nodof,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
c
c                 send out information about private nodes
c
            call MPI_SEND (proc_nodes(proc)%private,
     &           num_private_nodes(proc),MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
c
c                 send out information about shared but owned nodes
c
            call MPI_SEND (proc_nodes(proc)%own_shared,
     &           num_own_shared_nodes(proc),MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
            call MPI_SEND (proc_nodes(proc)%sharing_count,
     &           numprocs,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
c
c                    send out the pointers which indicate which shared and
c                    owned node are shared with each processor
c
            do i=0, numprocs-1
               if ( proc_nodes(proc)%sharing_count(i) .eq. 0) cycle
               call MPI_SEND (proc_nodes(proc)%sharing_procs(i)%ptr,
     &              proc_nodes(proc)%sharing_count(i),MPI_INTEGER,proc,
     &              proc,MPI_COMM_WORLD,ierr)
            enddo
c
c                 send out information about shared and not owned nodes
c
            call MPI_SEND (proc_nodes(proc)%shared,
     &           num_shared_nodes(proc),MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
            call MPI_SEND (proc_nodes(proc)%shared_count,
     &           numprocs,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
c
c                    send out the pointers which indicate which shared but
c                    not owned node are owned by each processor
c
            do i=0, numprocs-1
               if ( proc_nodes(proc)%shared_count(i) .eq. 0) cycle
               call MPI_SEND (proc_nodes(proc)%shared_owner(i)%ptr,
     &              proc_nodes(proc)%shared_count(i),MPI_INTEGER,proc,
     &              proc,MPI_COMM_WORLD,ierr)
            enddo
c
c                 send out information about internal blocks
c
            call MPI_SEND (proc_nodes(proc)%internal_blks,
     &           proc_nodes(proc)%num_int_blks,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
c
         enddo
c
      else
c     
c              get the constants
c
         call MPI_RECV (num_private,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,
     &        status, ierr)
         call MPI_RECV (num_own_shared,1,MPI_INTEGER,0,myid,
     &        MPI_COMM_WORLD, status, ierr)
         call MPI_RECV (num_shared,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,
     &        status, ierr)
         call MPI_RECV (local_nodes%num_int_blks,1,MPI_INTEGER,0,myid,
     &        MPI_COMM_WORLD, status, ierr)
c
c              allocate the space
c     
         call allocate_owner(num_private, num_own_shared, num_shared,
     &        local_nodes%num_int_blks, local_nodes)
c     
c              now fill the arrays
c
c                 recv info about local to/from global mapping
c
         call MPI_RECV (local_nodes%local2global,
     &        local_nodes%num_local_nodes*mxndof,MPI_INTEGER,0,
     &        myid,MPI_COMM_WORLD,status,ierr)
         call MPI_RECV (local_nodes%global2local,
     &        nodof,MPI_INTEGER,0,
     &        myid,MPI_COMM_WORLD,status,ierr)
c     
c                 recv info about private nodes
c
         call MPI_RECV (local_nodes%private,
     &        num_private,MPI_INTEGER,0,
     &        myid,MPI_COMM_WORLD,status,ierr)
c     
c                 recv info about shared but owned nodes
c
         call MPI_RECV (local_nodes%own_shared,
     &        num_own_shared,MPI_INTEGER,0,
     &        myid,MPI_COMM_WORLD,status,ierr)
         call MPI_RECV (local_nodes%sharing_count,numprocs,
     &        MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
c
c                    recv the pointers which indicate which shared and
c                    owned node are shared by each processor.  Allocate
c                    each pointer vector when we know how big it should be.
c
         do i=0, numprocs - 1
            if ( local_nodes%sharing_count(i) .eq. 0) cycle
            allocate(local_nodes%sharing_procs(i)%ptr(
     &           local_nodes%sharing_count(i)))
            call MPI_RECV (local_nodes%sharing_procs(i)%ptr,
     &           local_nodes%sharing_count(i),MPI_INTEGER,0,
     &           myid,MPI_COMM_WORLD,status,ierr)
         enddo
c     
c                 recv info about shared but not owned nodes
c
         call MPI_RECV (local_nodes%shared,
     &        num_shared,MPI_INTEGER,0,
     &        myid,MPI_COMM_WORLD,status,ierr)
         call MPI_RECV (local_nodes%shared_count,numprocs,
     &        MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
c
c                    recv the pointers which indicate which shared but
c                    not owned node are owned by each processor.  Allocate
c                    each pointer vector when we know how big it should be.
c
         do i=0, numprocs-1
            if ( local_nodes%shared_count(i) .eq. 0) cycle
            allocate(local_nodes%shared_owner(i)%ptr(
     &           local_nodes%shared_count(i)))
            call MPI_RECV (local_nodes%shared_owner(i)%ptr,
     &           local_nodes%shared_count(i),MPI_INTEGER,0,
     &           myid,MPI_COMM_WORLD,status,ierr)
         enddo
c     
c                 recv info about internal blocks
c
         call MPI_RECV (local_nodes%internal_blks,
     &        local_nodes%num_int_blks,
     &        MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
c
      endif
c
c           print out the data structures if debug is on.
c
      if (debug) then
         do proc = 0, numprocs - 1
            if (myid .eq. proc) then
               write (out,*) '>=>=>=>=>=> This is proc', myid
               call print_node_own_all(local_nodes)
            endif
            call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         enddo
      endif
c
c           create data types for transmission of the extra nodes
c
      do proc = 0, numprocs - 1
c
         if ( local_nodes%sharing_count(proc) .gt. 0) then
c
            do i=1, local_nodes%sharing_count(proc)
               disp(i) = (local_nodes%sharing_procs(proc)%ptr(i)-1) * 3
               blksz(i) = 3
            enddo
            call MPI_TYPE_INDEXED (local_nodes%sharing_count(proc),
     &           blksz,disp,MPI_VAL, local_nodes%MPI_sharing_type(proc), 
     &           ierr)
            call MPI_TYPE_COMMIT (local_nodes%MPI_sharing_type(proc), 
     &           ierr)
c
         endif
c
         if ( local_nodes%shared_count(proc) .gt. 0) then
c
            do i=1, local_nodes%shared_count(proc)
               disp(i) = (local_nodes%shared_owner(proc)%ptr(i)-1) * 3
               blksz(i) = 3
            enddo
            call MPI_TYPE_INDEXED (local_nodes%shared_count(proc),
     &           blksz,disp,MPI_VAL, local_nodes%MPI_shared_type(proc),
     &           ierr)
            call MPI_TYPE_COMMIT (local_nodes%MPI_shared_type(proc),
     &           ierr)
c
         endif
c
      enddo
c
c                  create a new local version of edest_blocks which points
c                  directly to the locally-renumbers dofs instead of the 
c                  global dofs.
c
      call ledest_init
c
c                  allocate local_mdiag and kdiag
c
      call mem_allocate ( 11 )
      call mem_allocate ( 18 )
c
c
      if ( debug ) write (out,*) myid,': <<- leaving wmpi_owner_send'
      return
 9900 format('>>> FATAL ERROR: memory allocate failure...')
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine allocate_owner               *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/04/98                   *
c     *                                                              *
c     *     this subroutine calculates the node owner structure      *
c     *     for domain decomposition portions of the MPI code        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine allocate_owner(num_priv, num_own_shar, num_shar,
     &     num_int_blks, own)
      use mpi_lnpcg
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      type (node_owner_type) :: own
c
      dimension status(MPI_STATUS_SIZE)
c     
c         allocate pointer from local dofs to global dofs
c
      own%num_local_nodes = num_priv + num_own_shar + num_shar
      allocate(own%local2global(own%num_local_nodes*mxndof), 
     &     stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
      own%local2global(1:own%num_local_nodes*mxndof) = 0
      allocate(own%global2local(nodof), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
      own%global2local(1:nodof) = 0
c     
c         allocate list of private nodes -- nodes only referenced by
c         the owning processor
c
      own%num_private = num_priv
      allocate (own%private(num_priv), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
      own%private(1:num_priv) = 0
c
c         allocate list of own_shared (sharing) nodes -- nodes owned 
c         by this processor, but referenced by others
c
      own%num_own_shared = num_own_shar
      allocate (own%own_shared(num_own_shar), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
      own%own_shared(1:num_own_shar) = 0
c
      allocate (own%sharing_count(0:numprocs-1), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
      own%sharing_count(0:numprocs-1) = 0
c
c            allocate ptr vector which indicates which localled owned but
c            shared nodes are shared with each processor
c
      allocate (own%sharing_procs(0:numprocs-1), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
      allocate (own%MPI_sharing_type(0:numprocs-1), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
c         allocate list of shared nodes -- nodes referenced by this
c         processor , but owned by another
c
      own%num_shared = num_shar
      allocate (own%shared(num_shar), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
      own%shared(1:num_shar) = 0
c
      allocate (own%shared_count(0:numprocs-1), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
      own%shared_count(0:numprocs-1) = 0
c
c            allocate ptr vector which indicates which shared but not 
c            shared nodes are shared with each processor
c
      allocate (own%shared_owner(0:numprocs-1), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
      allocate (own%MPI_shared_type(0:numprocs-1), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
c         allocate list of blocks with no elements on boundary.
c
      own%num_int_blks = num_int_blks
      allocate (own%internal_blks(num_int_blks), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
c
      return
 9900 format('>>> FATAL ERROR: memory allocate failure...')
      end
c
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine deallocate_owner             *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/04/98                   *
c     *                                                              *
c     *     this subroutine calculates the node owner structure      *
c     *     for domain decomposition portions of the MPI code        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine deallocate_owner(own)
      use mpi_lnpcg
      implicit integer (a-z)
      include "mpif.h"
      type (node_owner_type) :: own
c
c                       local declarations
c
      dimension status(MPI_STATUS_SIZE)
c     
c         allocate pointer from local dofs to global dofs
c
      deallocate (own%local2global)
      deallocate (own%global2local)
      deallocate (own%private)
      deallocate (own%own_shared)
      deallocate (own%sharing_count)
      deallocate (own%sharing_procs)
      deallocate (own%MPI_sharing_type)
      deallocate (own%shared)
      deallocate (own%shared_count)
      deallocate (own%shared_owner)
      deallocate (own%MPI_shared_type)
      deallocate (own%internal_blks)
c
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine print_node_own_all           *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/04/98                   *
c     *                                                              *
c     *     this subroutine prints the entire node owner structure   *
c     *     for domain decomposition portions of the MPI code        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine print_node_own_all(own)
      use mpi_lnpcg
      implicit integer (a-z)
$add common.main
      type (node_owner_type) :: own
c
c
      write (out,*) '>>> # of nodes:',own%num_local_nodes
      write (out,*) '>>> local2global:'
      write (out,'(5x,12i6)') (own%local2global(i),i=1,
     &     own%num_local_nodes*mxndof)
c
      write (out,*) '>>> global2local:'
      write (out,'(17x,"1",17x,"2",17x,"3",17x,"4")')
      write (out,'(5x,12i6)') (own%global2local(i),i=1,nodof)
c
      write (out,*) '>>> private nodes:',own%num_private
      write (out,'(5x,10i6)') (own%private(i),i=1,own%num_private)
c
      write (out,*) '>>> owned but shared nodes:',own%num_own_shared
      write (out,'(5x,10i6)') (own%own_shared(i),i=1,own%num_own_shared)
      do proc=0, numprocs-1
         write (out,*) '   ==> nodes shared w/ proc',proc,':',
     &        own%sharing_count(proc)
         if (own%sharing_count(proc) .eq. 0) cycle
         write (out,'(8x,12i6)') 
     &        (own%sharing_procs(proc)%ptr(i),
     &        i=1,own%sharing_count(proc))
      enddo
c
      write (out,*) '>>> shared nodes:',own%num_shared
      write (out,'(5x,10i6)') (own%shared(i),i=1,own%num_shared)
      do proc=0, numprocs-1
         write (out,*) '   ==> nodes owned by proc',proc,':',
     &        own%shared_count(proc)
         if ( own%shared_count(proc) .eq. 0 ) cycle
         write (out,'(8x,12i6)')
     &        (own%shared_owner(proc)%ptr(i),i=1,own%shared_count(proc))
      enddo
c
      write (out,*) '>>> list of internal blks'
      write (out,'("      internal:",20i3)') 
     &     (own%internal_blks(j),j=1,own%num_int_blks)
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine print_node_own               *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/04/98                   *
c     *                                                              *
c     *     this subroutine prints selected information about node   *
c     *     ownership for domain decomposition portions of the MPI   *
c     *     code.                                                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine print_node_own(own, procnum)
      use mpi_lnpcg
      implicit integer (a-z)
$add common.main
      type (node_owner_type) :: own
c
c
      sharing_procs = 0
      owning_procs = 0
      do proc = 0, numprocs - 1
	 if ( proc .eq. procnum) cycle
	 if ( own%sharing_count(proc) .gt. 0 ) 
     &        sharing_procs = sharing_procs  + 1
	 if ( own%shared_count(proc) .gt. 0 ) 
     &        owning_procs = owning_procs + 1
      enddo
c
      write (out,'(5x,"Processor ",i3,":")') procnum
      write (out,'(10x,"total # of nodes     :",i8)')own%num_local_nodes
      write (out,'(10x,"private nodes        :",i8," (",f5.1,"%)")')
     &       own%num_private, 100.0 *
     &       float(own%num_private)/float(own%num_local_nodes)
c
      write (out,'(10x,"nodes owned & shared :",i8," (",f5.1,"%)")')
     &       own%num_own_shared, 100.0 *
     &       float(own%num_own_shared)/float(own%num_local_nodes)
      write (out,'(10x,"# of sharing procs   :",i8)') sharing_procs
c
      write (out,'(10x,"nodes not owned      :",i8," (",f5.1,"%)")')
     &       own%num_shared, 100.0 *
     &       float(own%num_shared)/float(own%num_local_nodes)
      write (out,'(10x,"# of owning procs    :",i8,/)') owning_procs
c
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine copy_own                     *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/05/98                   *
c     *                                                              *
c     *     this subroutine copies the node owner structure          *
c     *     from one structure to another, allocated to the same     *
c     *     size. Assumes that the constants have already been       *
c     *     set in the new copy.                                     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine copy_own(old, new)
      use mpi_lnpcg
      implicit integer (a-z)
$add common.main
      type (node_owner_type) :: old, new
c
c
      do i=1, old%num_local_nodes * mxndof
         new%local2global(i) = old%local2global(i)
      enddo
      do i=1, nodof
         new%global2local(i) = old%global2local(i)
      enddo
c
      do i=1, old%num_private
         new%private(i) = old%private(i)
      enddo
c
c
      do i=1, old%num_own_shared
         new%own_shared(i) = old%own_shared(i)
      enddo
c
      do i=0, numprocs-1
         new%sharing_count(i) = old%sharing_count(i) 
         if (new%sharing_count(i) .eq. 0 ) cycle
c
	 allocate (new%sharing_procs(i)%ptr(new%sharing_count(i)),
     &        stat = allo_stat )
	 if ( allo_stat .ne. 0) then
	     write(out,*) '>>>> allocate error in copy_own'
	     call die_abort
	     stop
         endif
c
         do j=1, new%sharing_count(i)
            new%sharing_procs(i)%ptr(j) = old%sharing_procs(i)%ptr(j)
         enddo
c 
      enddo
c
c
      do i=1, old%num_shared
         new%shared(i) = old%shared(i)
      enddo
c
      do i=0, numprocs-1
         new%shared_count(i) = old%shared_count(i) 
         if (new%shared_count(i) .eq. 0 ) cycle
c
	 allocate (new%shared_owner(i)%ptr(new%shared_count(i)),
     &        stat = allo_stat )
	 if ( allo_stat .ne. 0) then
	     write(out,*) '>>>> allocate error in copy_own'
	     call die_abort
	     stop
         endif
c
         do j=1, new%shared_count(i)
            new%shared_owner(i)%ptr(j) = old%shared_owner(i)%ptr(j)
         enddo
      enddo
c
      new%num_int_blks = old%num_int_blks
      do i=1, new%num_int_blks
         new%internal_blks(i) = old%internal_blks(i) 
      enddo
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine dd_pcg_init                  *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/06/98                   *
c     *                                                              *
c     *     this subroutine sends the initial copies of diag and     *
c     *     lres to lnpcg.  The sent copies have their dofs          *
c     *     rearranged so that all the private nodes are first,      *
c     *     the shared-but-owned nodes are next, then the shared     *
c     *     only nodes are last.                                     *
c     *                                                              *
c     ****************************************************************
c
      subroutine dd_pcg_init( lres_local, svec_local, svec, 
     &     mydof, refdof)
      use mpi_lnpcg
      use pcg_data, only: lres
      implicit integer (a-z)
      include "mpif.h"
$add common.main
#dbl      double precision
#sgl      real
     &     lres_local(*), svec_local(*), svec(*), zero
      data zero /0.0/
c
c
c
      mydof = (local_nodes%num_private + local_nodes%num_own_shared)*3
      refdof = local_nodes%num_local_nodes * 3
c
c            allocate pcg data structures if not done already
c
      if (.not. allocated (lres) ) call mem_allocate ( 17 )
c
c            send lres to all processors
c
      if ( slave_processor ) then
         lres(1:nodof) = zero
      endif
      call MPI_BCAST (lres, nodof, MPI_VAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST (svec, nodof, MPI_VAL, 0, MPI_COMM_WORLD, ierr)
c      
c            convert diag and lres to local number sequence
c     
      do i=1, local_nodes%num_local_nodes*mxndof
         lres_local(i)=lres(local_nodes%local2global(i))
         svec_local(i)=svec(local_nodes%local2global(i))
      enddo
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine dd_combine                   *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/06/98                   *
c     *                                                              *
c     *     this subroutine takes an array divided between           *
c     *     processors and combines them in the correct global       *
c     *     ordering back on processor 0.                            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine dd_combine(array_in, array_out)
      use mpi_lnpcg
      implicit integer (a-z)
      include "mpif.h"
$add common.main
#dbl      double precision
#sgl      real
     &     array_in(*), array_tmp(mxdof), zero, array_out(*)
      data zero /0.0/
c
c
c
      mydof = (local_nodes%num_private + local_nodes%num_own_shared)*3
c
c          copy values related to owned nodes to global data strucutre
c
      array_tmp(1:nodof) = zero
      do i = 1, mydof
         array_tmp(local_nodes%local2global(i)) = array_in(i)
      enddo
c
c          now do a reduction to get all the copies of the global
c          data structure to root proc
c
      call MPI_REDUCE (array_tmp, array_out, nodof, 
#dbl     &     MPI_DOUBLE_PRECISION,
#sgl     &     MPI_REAL,
     &     MPI_SUM,0,MPI_COMM_WORLD, ierr)
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine dd_print                     *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/06/98                   *
c     *                                                              *
c     *     this subroutine takes an array divided between           *
c     *     processors and combines them in the correct global       *
c     *     ordering back on processor 0, then prints the array.     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine dd_print(array)
      use mpi_lnpcg
      implicit integer (a-z)
      include "mpif.h"
$add common.main
#dbl      double precision
#sgl      real
     &     array(*), local(mxdof), zero
      data zero /0.0/
c
c
c
      call dd_combine ( array, local ) 
c     
      if ( root_processor ) then
         write (out,*) '    =>   array is:'
         do i=1, nodof
            write (out,'(5x,i5,2x,e20.12)') i, local(i)
         enddo
      endif
c     
      call MPI_BARRIER ( MPI_COMM_WORLD, ierr )
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine lfinda_comm                  *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/26/98                   *
c     *                                                              *
c     *     this subroutine updates values in a distributed vector   *
c     *     that lie on the boundary between processors.  For a node *
c     *     on a boundary, the sharing processors each have a        *
c     *     contribution to the final term; these contributions      *
c     *    must be gathered together, then the final result is sent  *
c     *    to each of the sharing processors. Combination of terms   *
c     *     can either be additive or multiplicative.                *
c     *
c     *     communication occurs using nonblocking sends and         *
c     *     recieves to improve parallel efficiency                  *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine lnpcg_comm(vec, mydof, refdof, oper )
      use mpi_lnpcg
      implicit integer (a-z)
      include "mpif.h"
$add common.main
#dbl      double precision
#sgl      real
     &     vec_tmp(mxdof), vec(*), zero
      dimension status (MPI_STATUS_SIZE), request(2)
      logical debug, found_proc
      data zero, debug / 0.0, .false. /
c
c
c           oper specifies if combination of contributions is additive of
c           mulitpolictaive.
c               oper = 1 -- add 
c               oper = 2 -- mult
c
c
      if (numprocs .eq. 1) return
c
c            exchange only the critical info that the other
c            processors need to know.
c           
c                first gather shared terms into the owning
c                processors for full calculation of vec.
c           
c           
c                   start by initiating a send request and a 
c                   receive request -- find first processor
c                   to send to and to receive from, and initiate
c                   those processes.
c
      if (debug) write (out,*) myid,':>>> start vec gather.'
      send_proc = -1
      do proc = 0, numprocs -1 
         if ( local_nodes%shared_count(proc) .eq. 0 ) then
            cycle
         else
            send_proc = proc
            exit
         endif
      enddo
      if (debug) write (out,*) myid,':   => send_proc =',send_proc
      if ( send_proc .ne. -1) then 
         call MPI_ISEND(vec(mydof+1),1,
     &        local_nodes%MPI_shared_type(send_proc),
     &        send_proc, send_proc, MPI_COMM_WORLD, request(1), ierr)
      endif
c
      recv_proc = -1
      do proc = 0, numprocs -1 
         if ( local_nodes%sharing_count(proc) .eq. 0 ) then
            cycle
         else
            recv_proc = proc
            exit
         endif
      enddo
      if (debug) write (out,*) myid,':   => recv_proc =',recv_proc
      if ( recv_proc .ne. -1) then 
         vec_tmp(1:refdof) = zero
         call MPI_IRECV(vec_tmp,1,
     &        local_nodes%MPI_sharing_type(recv_proc),
     &        recv_proc, myid, MPI_COMM_WORLD, request(2), ierr)
      endif
c
c                   now run a loop, waiting for the next
c                   completed communication operation. Do
c                   the send if the send completes, or 
c                   process the recv and do the next recv
c                   if the recv completes.
c
      do while (.true.)
c
c                       if both send_proc and recv_proc are -1, then
c                       we have completed all the sends and receives.
c                       exit this loop.
c     
         if ( send_proc .eq. -1 .and. recv_proc .eq. -1 ) then
            if (debug) 
     &           write(out,*)myid,':>>> all requests are done. Move on.'
            exit
         endif
c
c                       wait for next completed communication request
c         
         if (debug)write (out,*) myid,':  >> waiting for a request'
         call MPI_WAITANY ( 2, request, do_req, status, ierr)
         if (debug)write (out,*) myid,':  >> got a request =',do_req
c     
c                       process completed request.
c                           1 = do next send
c                           2 = do next recv
c
         if ( do_req .eq. 1) then
c
c                           doing send.  find next proc to send. If we are
c                           out of sends to do, set send_proc to -1 and cycle.
c            
            found_proc = .false.
            do proc = send_proc + 1, numprocs - 1 
               if ( local_nodes%shared_count(proc) .eq. 0 ) then
                  cycle
               else
                  send_proc = proc
                  found_proc = .true.
                  exit
               endif
            enddo
            if ( .not. found_proc ) then
               send_proc = -1
               if (debug)write(out,*)myid,':   => send_proc =',send_proc
               cycle
            else
               if (debug)write(out,*)myid,':   => send_proc =',send_proc
c
c                           we have found a send to do ... do it.
c
               call MPI_ISEND(vec(mydof+1),1,
     &              local_nodes%MPI_shared_type(send_proc),
     &              send_proc, send_proc, MPI_COMM_WORLD, request(1), 
     &              ierr)
            endif
c
         else
c
c                           do the next recv.  Process the recv that was
c                           just completed, and then set up another.
c
            do i=1, local_nodes%sharing_count(recv_proc)
               lindx =( local_nodes%sharing_procs(recv_proc)%ptr(i) 
     &                  - 1 ) * 3
               gindx =( local_nodes%num_private +
     &              local_nodes%sharing_procs(recv_proc)%ptr(i) 
     &                  - 1 ) * 3 
               if ( oper .eq. 1 ) then
                  vec(gindx+1) = vec(gindx+1) + vec_tmp(lindx+1)
                  vec(gindx+2) = vec(gindx+2) + vec_tmp(lindx+2)
                  vec(gindx+3) = vec(gindx+3) + vec_tmp(lindx+3)
               else if ( oper .eq. 2 ) then
                  vec(gindx+1) = vec(gindx+1) * vec_tmp(lindx+1)
                  vec(gindx+2) = vec(gindx+2) * vec_tmp(lindx+2)
                  vec(gindx+3) = vec(gindx+3) * vec_tmp(lindx+3)
               endif
            enddo
c
c                           now intiate next recv request. find next proc 
c                           to recv. If we are out of recvs to do, set 
c                           recv_proc to -1 and cycle.
c            
            found_proc = .false.
            do proc = recv_proc + 1, numprocs - 1 
               if ( local_nodes%sharing_count(proc) .eq. 0 ) then
                  cycle
               else
                  recv_proc = proc
                  found_proc = .true.
                  exit
               endif
            enddo
            if ( .not. found_proc ) then
               recv_proc = -1
               if (debug)write(out,*)myid,':   => recv_proc =',recv_proc
               cycle
            else
c
c                           we have found a send to do ... do it.
c
               if (debug)write(out,*)myid,':   => recv_proc =',recv_proc
               vec_tmp(1:refdof) = zero
               call MPI_IRECV(vec_tmp,1,
     &              local_nodes%MPI_sharing_type(recv_proc),
     &              recv_proc, myid, MPI_COMM_WORLD, request(2), ierr)
            endif
c
         endif
c
      enddo
c
c
      if (debug) then
         if (myid .eq. 0) write (out,*) '>>>>> this is vec:'
         call dd_print(vec)
      endif
c
c                Now we have processed all of the contibutions to the
c                owned but shared dofs, and each processor has a 
c                correct version of vec for the nodes it owns.  Now
c                send the correct vec terms to the processors which
c                reference that information.
c           
c                   start by initiating a send request and a 
c                   receive request -- find first processor
c                   to send to and to receive from, and initiate
c                   those processes.
c
      if (debug) write (out,*) myid,':>>> start good vec scatter.'
      send_proc = -1
      do proc = 0, numprocs -1 
         if ( local_nodes%sharing_count(proc) .eq. 0 ) then
            cycle
         else
            send_proc = proc
            exit
         endif
      enddo
      if (debug)write (out,*) myid,':   => send_proc =',send_proc
      if ( send_proc .ne. -1) then 
         call MPI_ISEND(vec(local_nodes%num_private*3+1),1,
     &        local_nodes%MPI_sharing_type(send_proc),
     &        send_proc, send_proc, MPI_COMM_WORLD, request(1), ierr)
      endif
c
      recv_proc = -1
      do proc = 0, numprocs -1 
         if ( local_nodes%shared_count(proc) .eq. 0 ) then
            cycle
         else
            recv_proc = proc
            exit
         endif
      enddo
      if (debug)write (out,*) myid,':   => recv_proc =',recv_proc
      if ( recv_proc .ne. -1) then 
         vec_tmp(1:refdof) = zero
         call MPI_IRECV(vec(mydof + 1),1,
     &        local_nodes%MPI_shared_type(recv_proc),
     &        recv_proc, myid, MPI_COMM_WORLD, request(2), ierr)
      endif
c
c                   now run a loop, waiting for the next
c                   completed communication operation. Do
c                   the send if the send completes, or 
c                   process the recv and do the next recv
c                   if the recv completes.
c
      do while (.true.)
c
c                       if both send_proc and recv_proc are -1, then
c                       we have completed all the sends and receives.
c                       exit this loop.
c     
         if ( send_proc .eq. -1 .and. recv_proc .eq. -1 ) then
            if (debug)
     &           write(out,*)myid,':>>> all requests are done. Move on.'
            exit
         endif
c
c                       wait for next completed communication request
c         
         if (debug)write (out,*) myid,':  >> waiting for a request'
         call MPI_WAITANY ( 2, request, do_req, status, ierr)
         if (debug)write (out,*) myid,':  >> got a request =',do_req
c     
c                       process completed request.
c                           1 = do next send
c                           2 = do next recv
c
         if ( do_req .eq. 1) then
c
c                           doing send.  find next proc to send. If we are
c                           out of sends to do, set send_proc to -1 and cycle.
c            
            found_proc = .false.
            do proc = send_proc + 1, numprocs - 1 
               if ( local_nodes%sharing_count(proc) .eq. 0 ) then
                  cycle
               else
                  send_proc = proc
                  found_proc = .true.
                  exit
               endif
            enddo
            if ( .not. found_proc ) then
               send_proc = -1
               if (debug)write(out,*)myid,':   => send_proc =',send_proc
               cycle
            else
               if (debug)write(out,*)myid,':   => send_proc =',send_proc
c
c                           we have found a send to do ... do it.
c
               call MPI_ISEND(vec(local_nodes%num_private*3+1),1,
     &              local_nodes%MPI_sharing_type(send_proc),
     &              send_proc, send_proc, MPI_COMM_WORLD, request(1), 
     &              ierr)
            endif
c
         else
c
c                           do the next recv. find next proc 
c                           to recv. If we are out of recvs to do, set 
c                           recv_proc to -1 and cycle.
c            
            found_proc = .false.
            do proc = recv_proc + 1, numprocs - 1 
               if ( local_nodes%shared_count(proc) .eq. 0 ) then
                  cycle
               else
                  recv_proc = proc
                  found_proc = .true.
                  exit
               endif
            enddo
            if ( .not. found_proc ) then
               recv_proc = -1
               if (debug)write(out,*)myid,':   => recv_proc =',recv_proc
               cycle
            else
c
c                           we have found a send to do ... do it.
c
               if (debug)write(out,*)myid,':   => recv_proc =',recv_proc
               call MPI_IRECV(vec(mydof + 1),1,
     &              local_nodes%MPI_shared_type(recv_proc),
     &              recv_proc, myid, MPI_COMM_WORLD, request(2), ierr)
            endif
c
         endif
c
      enddo
c
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine ledest_init                  *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/27/98                   *
c     *                                                              *
c     *     create the blocked data structure for storage of element *
c     *     displacement indexes (cdest). fill integer indexes       *
c     *                                                              *
c     ****************************************************************
c
      subroutine ledest_init
c
      use elem_block_data, only : edest_blocks, edest_blk_list
      use mpi_lnpcg, only : ledest_blocks, ledest_blk_list, local_nodes
      use main_data, only : incid, incmap, elems_to_blocks
c
      implicit integer (a-z)     
$add common.main
c
      integer, dimension (:,:), pointer :: ledest
      logical myblk
c
c            the data structure is a vector for each element
c            block (dynamically allocated). the vectors are
c            hung from a dynamically allocated pointer vector
c            of length nelblk.
c
      allocate( ledest_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 1
           call die_abort
           stop
      end if
c
      allocate( ledest_blk_list(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 5
           call die_abort
           stop
      end if
      ledest_blk_list(1:nelblk) = 0
c
c            loop over all element blocks. allocate blocks.
c
      do blk = 1, nelblk
c
         myblk = myid .eq. elblks(2,blk)
         if ( root_processor ) myblk = .true.
         if (.not. myblk) cycle
c
         felem         = elblks(1,blk)
         span          = elblks(0,blk)
         nnode         = iprops(2,felem)
         num_enode_dof = iprops(4,felem)
         totdof        = nnode * num_enode_dof
c     
         allocate( ledest_blocks(blk)%ptr(totdof,span),
     &             stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
            write(out,9900)
            write(out,9910) 3
            call die_abort
            stop
         end if
         ledest_blk_list(blk) = 1
      end do
c
c            loop over element blocks.  For all of the blocks that this
c            processor owns, copy over the edest_blocks structure, and
c            modify it to reflect the local renumbering.
c
      do blk = 1, nelblk
c
c         myblk = myid .eq. elblks(2,blk)
c         if ( root_processor ) myblk = .true.
c         if (.not. myblk) cycle
c
         if ( myid .ne. elblks(2,blk) ) cycle
c
         felem         = elblks(1,blk)
         span          = elblks(0,blk)
         nnode         = iprops(2,felem)
         num_enode_dof = iprops(4,felem)
         totdof        = nnode * num_enode_dof
c
         do j=1, span
            do i=1, totdof
               ledest_blocks(blk)%ptr(i,j) = 
     &              local_nodes%global2local(edest_blocks(blk)%ptr(i,j))
            enddo
         enddo
c
      end do
c
      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 ledest_init: @',i2,/,
     &       '>>> Job terminated....')
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine wmpi_chknode                     *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 04/06/98                   *
c     *                                                              *
c     *           The processor assesses whether the node number     *
c     *           passed in as an argument is referenced by this     *
c     *           processor.  If it is referenced, then it also      *
c     *           checks to see if this processor owns the node.     *
c     *                                                              *
c     ****************************************************************
c    
c
      subroutine wmpi_chknode ( node, referenced, owned )
      use mpi_lnpcg, only : local_nodes
      implicit integer (a-z)
$add common.main
      logical referenced, owned
      include "mpif.h"
c
c
      ptr = local_nodes%global2local(dstmap(node)) 
c
      if ( ptr .eq. 0 ) then
c
         referenced = .false.
         owned = .false.
c
      else
c
         referenced = .true.
c
         if ( (ptr+2) / 3 .le. (local_nodes%num_private + 
     &        local_nodes%num_own_shared) ) then
            owned = .true.
         else
            owned = .false.
         endif
c
      endif
c
c
      return 
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine find_blk_groups              *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 10/16/98                   *
c     *                                                              *
c     *      This subroutine drives the process of scheduling the    *
c     *      blocks of elements for the ebe preconditioner so that   *
c     *      a maximum amount of parallelism can be achieved.        *
c     *      blocks of elements with elements that lie on the        *
c     *      boundary between processors must be scheduled carefully *
c     *      so that processors do not simultaneously compute the    *
c     *      ebe contribution for elements which share nodes.        *
c     *                                                              *
c     **************************************************************** 
c 
c
      subroutine find_blk_groups (node2proc)
      use mpi_lnpcg
      use main_data, only: incid, incmap
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      dimension 
     &     ordering(max_procs * 5,mxnmbl), order_num(mxnmbl),
     &     node2proc(0:numprocs,nonode),
     &     max_blks(mxnmbl), num_blks(0:max_procs-1),
     &     max_elems(mxnmbl), num_elem(0:max_procs-1)
      logical debug
      data debug /.false./
c
      type (graph_type), dimension(:), allocatable  :: graph(:)
c
      if ( debug ) write (out,*) myid,': ->> in find_blk_groups'
c 
c           in order to overlap computation and communication for the
c           ebe preconditioner, first find all of the blocks which have a
c           element on the boundary.  Then break the boundary blocks into 
c           groups in which no two blocks share a node.  calculation of
c           boundary terms is thus parallel within a group. The following
c           code computes these groups and sends the info to the processors.
c
c              allocate data structure to hold graph data, init vars
c
      allocate (graph(nelblk), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
c             build graph structure which will describe the dependencies
c             between bounary blocks
c
      call build_graph ( node2proc, graph, tot_num_graph, debug )
c     
c	       call subroutine which, based on the dependency data, 
c              breaks the processors into groups, in which no two
c	       processors share a degree of freedom.
c
      call find_par_groups (graph, ordering, order_num, tot_num_graph )
c
c              create data structures in graph for link information
c
      call fill_grp_links ( node2proc, graph, tot_num_graph, debug )
c
c              if debug is on, print out the full graph
c
      if (debug) then
         write (*,*) '>>> Alright, folks, here is the graph structure:'
         do i=1, tot_num_graph
            write (*,*) ' -> for graph ',i
            call print_grnode(graph(i))
         enddo
      endif
c
c              calculate the maximum number of external blocks 
c              that are calculated during each parallel group.  This
c              is used to determine the ordering of operations for
c              the ebe calculations
c
      tot_blks = 0
      tot_elem = 0
      do grp = 1, num_groups
         max_blks(grp) = 0
         max_elems(grp) = 0
         num_blks(0:numprocs-1) = 0
         num_elem(0:numprocs-1) = 0
         do i = 1, order_num(grp)
            grnode = ordering(i,grp)
            owner = graph(grnode)%owner
            num_blks(owner) = num_blks(owner) + 
     &           graph(grnode)%num_blks
            num_elem(owner) = num_elem(owner) + 
     &           graph(grnode)%num_elem
            max_blks(grp) = max( max_blks(grp), num_blks(owner) )
            max_elems(grp) = max( max_elems(grp), num_elem(owner) )
         enddo
         tot_blks = tot_blks + max_blks(grp)
c
         min_elem = mxel
         do proc = 0, numprocs -1  
            min_elem = min(min_elem, num_elem(proc))
         enddo
c
         tot_elem = tot_elem + max_elems(grp) - min_elem
c
      enddo
c
c              Broadcast all the graph info to the various processors,
c              and determine the ordering of processes for the ebe
c              preconditioner application.
c
      call wmpi_graph_send (tot_num_graph, graph, max_elems)
c
c              deallocate global graph structure
c
      if (debug) write (out,*) '...deallocate'
      deallocate(graph)
c
      return
 9900 format('>>> FATAL ERROR: memory allocate failure...')
      end
c

c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine build_graph                  *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 10/16/98                   *
c     *                                                              *
c     *       This subroutine builds the dependency graph which      *
c     *       specifies which boundary blocks are linked to which    *
c     *       other boundary blocks. boundary blocks that have the   *
c     *       same processor owner and the same processor accesses   *
c     *       are lumped together and treated as a boundary group    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine build_graph (node2proc, graph, tot_num_graph, debug)
      use mpi_lnpcg
      use main_data, only: incid, incmap
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      dimension node2proc(0:numprocs,nonode)
      logical debug, found, same, found_proc, keep_on
c
      dimension blkprocs(max_procs), elem2grnode(mxel)
      type (graph_type) :: graph(*)
c
      if ( debug ) write (out,*) myid,': ->> in find_blk_groups'
c 
c           in order to overlap computation and communication for the
c           ebe preconditioner, first find all of the blocks which have a
c           element on the boundary.  Then break the boundary blocks into 
c           groups in which no two blocks share a node.  calculation of
c           boundary terms is thus parallel within a group. The following
c           code computes these groups and sends the info to the processors.
c
c             build graph structure; loop over all blocks
c
      tot_num_graph = 0
      do blk = 1, nelblk
         span   = elblks(0,blk)
         felem  = elblks(1,blk)
         owner_proc  = elblks(2,blk)
         nnode = iprops(2,felem)
c
c                for this block, find all procs which it references
c
         num_blkprocs = 0
         do elem = felem, felem + span - 1
            do i = 1, nnode
               node = incid(incmap(elem)+i-1)
               do j= 1, node2proc(0,node)
                  proc = node2proc(j,node)
                  if (proc .eq. owner_proc) cycle 
                  if ( num_blkprocs .eq. 0 ) then
                     num_blkprocs = 1
                     blkprocs(1) = proc
                  else
                     found = .false.
                     do k = 1,  num_blkprocs
                        if (proc .eq. blkprocs(k)) then
                           found = .true.
                           exit
                        endif
                     enddo
                     if ( .not. found) then
                        num_blkprocs = num_blkprocs + 1
                        blkprocs(num_blkprocs) = proc
                     endif
                  endif
               enddo
            enddo
         enddo
c
c                if block accesses no processors other than itself,
c                then it only contains internal elements (not on the 
c                boundary).  We do not include internal blocks in
c                the graph, so go to the next block.
c
         if ( num_blkprocs .eq. 0 ) cycle
c
c                now loop over previous graph nodes to see if this blk
c                fits into a previous graph node.  If so, add it, otherwise,
c                create a new graph node
c     
         found = .false.
         do i=1, tot_num_graph
c
c                   if blk is owned by different processors than the graph
c                   node, or if the number of processors accessed by the 
c                   graph node is different than the current block, then
c                   skip to the next block
c
            if ( graph(i)%owner .ne. elblks(2,blk)) cycle
            if ( graph(i)%num_proc_share .ne. num_blkprocs) cycle
c
c                   check if the processors accessed by this block
c                   are the same as those accessed by the others in
c                   the current graph node.
c
            same = .true.
            do j = 1, num_blkprocs
               proc = blkprocs(j)
               found_proc = .false.
               do k = 1, graph(i)%num_proc_share
                  if (proc .eq. graph(i)%proc_share(k)) then
                     found_proc = .true.
                     exit
                  endif
               enddo
               if (.not. found_proc) then
                  same = .false.
                  exit
               endif
            enddo
c
            if ( .not. same) cycle
c
c                   graph node accesses match current block.  add block 
c                   to list of blocks for this graph node.
c
            graph(i)%num_blks = graph(i)%num_blks + 1
            graph(i)%blks(graph(i)%num_blks) = blk
            graph(i)%num_elem = graph(i)%num_elem + elblks(0,blk)
            found = .true.
            exit
c                     
         enddo
c
         if ( .not. found ) then
c
c                   this block has different access requirements than
c                   any currently existing graph nodes.  create a new
c                   graph node for this block.
c
            tot_num_graph = tot_num_graph + 1
            call allo_grnode (mxnmbl, num_blkprocs, mxnmbl, 
     &           mxnmbl, graph(tot_num_graph))
            graph(tot_num_graph)%id = tot_num_graph
            graph(tot_num_graph)%num_blks = 1
            graph(tot_num_graph)%num_ext_links = 0
            graph(tot_num_graph)%num_int_links = 0
            graph(tot_num_graph)%blks(1) = blk
            graph(tot_num_graph)%num_elem = elblks(0,blk)
            graph(tot_num_graph)%owner = owner_proc
            graph(tot_num_graph)%num_proc_share = num_blkprocs
            do j=1, num_blkprocs
               graph(tot_num_graph)%proc_share(j) = blkprocs(j)
            enddo
c
c                   since we created a new graph node, go through all
c                   the other graph nodes and update the graph links.
c     
            found = .false.
            do i = 1, tot_num_graph - 1
c
c                      if graph node is owned by same processor as the 
c                      new graph node, register a link between them.  Thus
c                      all graph nodes owned by the same processor are
c                      joined by links.
c
               if ( graph(tot_num_graph)%owner .eq. 
     &              graph(i)%owner) then
                  graph(i)%num_int_links = 
     &                 graph(i)%num_int_links + 1
                  graph(tot_num_graph)%num_int_links = 
     &                 graph(tot_num_graph)%num_int_links + 1
                  graph(i)%int_links(graph(i)%num_int_links) = 
     &                 tot_num_graph
                  graph(tot_num_graph)%int_links(
     &                 graph(tot_num_graph)%num_int_links) = i
c
c                      if graph node is owned by different processor than 
c                      the new graph node, register a link between them if
c                      they reference each other's processor. 
c
               else
                  keep_on = .true.
                  do j = 1, graph(tot_num_graph)%num_proc_share
                     proc1 = graph(tot_num_graph)%proc_share(j)
                     if ( graph(i)%owner .ne. proc1) cycle
                     do k = 1, graph(i)%num_proc_share
                        proc2 = graph(i)%proc_share(k)
                        if ( proc2 .eq. graph(tot_num_graph)%owner) then
                           graph(i)%num_ext_links = 
     &                          graph(i)%num_ext_links + 1
                           graph(tot_num_graph)%num_ext_links = 
     &                          graph(tot_num_graph)%num_ext_links + 1
                           graph(i)%ext_links(graph(i)%num_ext_links) = 
     &                          tot_num_graph
                           graph(tot_num_graph)%ext_links(
     &                          graph(tot_num_graph)%num_ext_links) = i
                           keep_on = .false.
                           exit
                        endif
                     enddo
                     if ( .not. keep_on) exit
                  enddo
c
               endif
c
            enddo
c     
         endif
c
      enddo
c     
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine find_par_groups                  *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 10/13/98                   *
c     *                                                              *
c     *           This subroutine breaks the processors into         *
c     *           groups in which no processors are adjacent.        *
c     *           This is essentially a graph coloring algorithm     *
c     *           but which is balanced so that each group contains  *
c     *           roughly the same amount of work.                   *
c     *                                                              *
c     ****************************************************************
c    
c      
      subroutine find_par_groups (graph, ordering, order_num,
     &     graph_size )
      use mpi_lnpcg
      implicit integer (a-z)
      type (graph_type) :: graph(graph_size)
c
$add common.main
c
c           argument variables
c
      dimension ordering(max_procs * 5,*), order_num(*)
c 
c           local variables -- notice here that we assume that the
c           largest graph size is smaller than the number of 
c           blocks in the structure.
c
      dimension list(mxnmbl), conn(max_procs * 5,mxnmbl), 
     &     conn_num(mxnmbl), 
     &     grp_w_grph(mxnmbl), proc_elems(0:max_procs - 1 ),
     &     proc_list(0:max_procs - 1 )
      logical found, assigned(mxnmbl), conflict, debug
      data debug /.false./
c     
c          copy graph into single data structure for graph coloring 
c          routine
c
      do grnode = 1, graph_size
         conn_num(grnode) = graph(grnode)%num_int_links +
     &        graph(grnode)%num_ext_links
         do i = 1, graph(grnode)%num_int_links
            conn(i,grnode) = graph(grnode)%int_links(i)
         enddo
         do i = 1, graph(grnode)%num_ext_links
            conn(i+graph(grnode)%num_int_links,grnode) = 
     &           graph(grnode)%ext_links(i)
         enddo
      enddo
c
c
c         Write out data structure conn, which lists for each graph node
c         the other graph nodes on which it is dependent.
c
      if ( debug ) then
         write (*,*) '>> Dependency data structure: grnode, #dep, dep:'
         do i = 1, graph_size
            write (*,'(2x,i2,":",i2,"=",50(1x,i2))') i, conn_num(i),
     &           (conn(j,i),j=1,conn_num(i))
         enddo
      endif
c
c         Create a sorted list of the graph nodes in order of decreasing
c         number of dependent graph nodes -- thus the first element has
c         the largest number of graph nodes on which it is dependent, etc.
c
c         New try -- order based on number of elems instead.
c
c
      do i=1, graph_size
         list(i) = 0
      enddo
      list(1) = 1
      do i = 2, graph_size 
         found = .false.
         do j = 1, i-1
            if ( graph(i)%num_elem.gt.graph(list(j))%num_elem ) then
               do k = i, j+1, -1
                  list(k) = list(k-1)
               enddo
               list(j) = i
               found = .true.
               exit
            endif
         enddo
         if ( .not. found) list(i) = i
      enddo
c
c         Write out graph nodes in order of decreasing number of dependent
c         graph nodes
c
      if ( debug) then
         write (*,*) 'Sorted order: (graph, #blks, #elems)'
         do i=1, graph_size
            write (*,*) list(i), graph(list(i))%num_blks, 
     &           graph(list(i))%num_elem
         enddo
      endif
c     
c         now, based on the new ordering, do graph coloring to break the
c         graph nodes into non-conflicting groups.  Within each group,
c         no two graph nodes are dependent on each other, thus they can
c         caluclate border terms in parallel.
c
c            init vars
c
      do i=1, graph_size
         order_num(i) = 0
         assigned(i) = .false.
	 grp_w_grph(i) = 0
      enddo
      grp = 0
c
c            loop over groups until we have taken care of all graph nodes
c
      do while (.true.)
c
c               initialize vars
c
         found = .false.
         grp = grp + 1
         do i = 0, numprocs -1 
            proc_elems(i) = 0
            proc_list(i) = 0
         enddo
c
c               loop over graph nodes to find elements to put into group
c
         do k = 1, graph_size
            grnode = list(k)
            if (assigned(grnode)) cycle
c
c                   if first unassigned graph node we have found, start group
c
            if ( .not. found) then
               order_num(grp) = 1
               ordering(1,grp) = grnode
	       grp_w_grph(grnode) = grp
               assigned(grnode) = .true.
	       proc_elems(graph(grnode)%owner) = graph(grnode)%num_elem
               found = .true.
               cycle
            endif
c
c                   if not unassigned graph node we have found, check conflict
c                   and, if not found, add to group
c
            conflict = .false.
            do i = 1, order_num(grp)
               do j = 1, conn_num(grnode)
                  if (conn(j,grnode) .eq. ordering(i,grp)) then
                     conflict = .true.
                     exit
                  endif
               enddo
               if ( conflict ) exit
            enddo
            if ( .not. conflict) then
               order_num(grp) = order_num(grp) + 1
               ordering(order_num(grp),grp) = grnode
               assigned(grnode) = .true.
	       grp_w_grph(grnode) = grp
	       proc_elems(graph(grnode)%owner) = graph(grnode)%num_elem
            endif
c     
         enddo
c
c               if no unassigned graph node found, stop it
c
         if ( .not. found ) then
            num_groups = grp - 1
            exit
         endif
c
c               now check if there are any additional graph nodes that
c               can be used for load balance
c
c                  first sort the procs by # of elements in list
c
         proc_list(1) = 0
         do proc = 1, numprocs - 1
            found = .false.
            do chk = 1, proc 
	       if  (proc_elems(proc) .gt. proc_elems(proc_list(chk)))
     &              then
		  do k = proc+1, chk+1, -1
		     proc_list(k) = proc_list(k-1)
                  enddo
                  proc_list(chk) = proc
		  found = .true.
		  exit
               endif
            enddo
	    if ( .not. found) proc_list(proc+1) = proc
         enddo
c
c                 If only one proc has all the remaining groups, then
c                 put all the remaining graph nodes into this group.
c
         if ( proc_elems(proc_list(2)) .eq. 0 ) then
            proc = proc_list(1)
            do grnode = 1, graph_size
               if ( graph(grnode)%owner .ne. proc ) cycle
               if ( assigned(grnode) ) cycle
               order_num(grp) = order_num(grp) + 1
               ordering(order_num(grp),grp) = grnode
               assigned(grnode) = .true.
               grp_w_grph(grnode) = grp
               proc_elems(proc) = proc_elems(proc) 
     &              + graph(grnode)%num_elem
            enddo
            cycle
         endif
c
c                 go thru the list backwards now.  For each processor with 
c                 fewer elements than the largest, find any additional 
c                 non-conflicting graph nodes that can be processed 
c                 during this grouping
c
	 do indx = numprocs, 2, -1
c
c                    compute how many fewer elements this processor
c                    has in the group versus the maximum.
c
	    proc = proc_list(indx)
	    deficit = proc_elems(proc_list(1)) - proc_elems(proc)
c
c                    loop over all the graph nodes owned by this processor
c                    that have fewer of equal numbers of elements as the 
c                    current deficit
c
	    do grnode = 1, graph_size
	       if ( graph(grnode)%owner .ne. proc ) cycle
	       if ( graph(grnode)%num_elem .gt. deficit) cycle
	       if ( assigned(grnode) ) cycle
c
c 	                loop over external links for this graph node;
c                       skip if it has links to any node owned by another
c                       processor that is in this group
c
	       conflict = .false.
	       do link = 1, graph(grnode)%num_ext_links
		  link_node = graph(grnode)%ext_links(link)
		  if ( grp_w_grph(link_node) .eq. grp ) then
		     conflict = .true.
		     exit
		  endif
               enddo
c
c                       if this graph node does not conflict with any others
c                       in this group, then add it to the group
c
	       if ( .not. conflict ) then
		  order_num(grp) = order_num(grp) + 1
                  ordering(order_num(grp),grp) = grnode
                  assigned(grnode) = .true.
	          grp_w_grph(grnode) = grp
       	          proc_elems(proc) = proc_elems(proc) 
     &                 + graph(grnode)%num_elem
	          deficit = proc_elems(proc_list(1)) - proc_elems(proc)
	       endif
c
	    enddo
	 enddo
c
      enddo
c     
c         copy group data into graph. 
c     
      do i = 1, num_groups
         do j = 1, order_num(i)
            graph(ordering(j,i))%order = i
         enddo
      enddo
c
c         Write out grouping
c
      if ( debug ) then
         write (*,*) 'Here are the ',num_groups,' groups'
         do i=1, num_groups
            write (*,'(2x,i2,":",i2,"=>",50(1x,i2))') i, order_num(i),
     &           (ordering(j,i),j=1,order_num(i))
         enddo
         write (*,*) 'Here are the groups for the graph nodes:'
         do i=1, graph_size
            write (*,'(2x,i5,2x,i5)') i,grp_w_grph(i)
         enddo
      endif
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine fill_grp_links               *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 10/16/98                   *
c     *                                                              *
c     *      This subroutine takes the graph information and         *
c     *      constructs a data strucutre describing all the links    *
c     *      and the ordering in which the boundary groups should    *
c     *      occur for parallel execution.                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine fill_grp_links ( node2proc, graph, tot_num_graph,
     &     debug )
      use mpi_lnpcg
      use main_data, only: incid, incmap, inverse_incidences 
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      dimension 
     &    grnd(mxconn), node2proc(0:numprocs,nonode), num_dofs(mxconn)
      logical debug, found, same, found_proc, keep_on, repeat_grnd
c
      dimension elem2grnode(mxel)
      type (graph_type) :: graph(*)
c
      if ( debug ) write (out,*) myid,': ->> in find_blk_groups'
c 
c              allocate data structures for link information
c
      do i = 1, tot_num_graph
         do j=1, graph(i)%num_ext_links
            num_dofs(j) = mxedof * mxvl * graph(i)%num_blks
         enddo
         call allo_link (graph(i)%num_ext_links, num_dofs, graph(i))
         do j=1, graph(i)%num_ext_links
            graph(i)%link(j)%owner = graph(graph(i)%ext_links(j))%owner
            graph(i)%link(j)%order = graph(graph(i)%ext_links(j))%order
            graph(i)%link(j)%id = graph(graph(i)%ext_links(j))%id
            graph(i)%link(j)%num_dofs = 0
         enddo
      enddo
c
c              for each graph node calculate which dofs are shared with
c              the link
c
c                 first create an element to graph node pointer vector.
c
      do elem = 1, noelem
         elem2grnode(elem) = 0
      enddo
c
      do grnode = 1, tot_num_graph
         do blk_loop = 1, graph(grnode)%num_blks
            blk = graph(grnode)%blks(blk_loop)
            span   = elblks(0,blk)
            felem  = elblks(1,blk)
            do elem = felem, felem + span - 1
               elem2grnode(elem) = grnode
            enddo
         enddo
      enddo
c
c                 now loop over all the nodes and, if shared by several
c                 groups, add the dofs to the corresponding link locations
c     
      do node = 1, nonode
         if (node2proc(0,node) .eq. 1 ) cycle
         do i=1, mxconn 
            grnd(i) = 0
         enddo
         num_grnd = 0
c
c                   find the groups which access this node. loop over the
c                   elements connected to the node and makes a list in
c                   grnd of the groups the elements belong to.
c         
         do elem_loop = 1, inverse_incidences(node)%element_count
            elem = inverse_incidences(node)%element_list(elem_loop)
            grnode = elem2grnode(elem)
            if ( grnode .lt. 1) cycle
            repeat_grnd = .false.
            do i = 1, num_grnd
               if (grnd(i) .eq. grnode) then
                  repeat_grnd = .true.
                  exit
               endif
            enddo
            if (repeat_grnd) cycle
            num_grnd = num_grnd + 1
            grnd(num_grnd) = grnode
         enddo
c
c                   for each group which accesses this node, put the 
c                   corresponding dofs into the link data strcuture
c
c                   NOTE: this assumes 3 dofs per node!!!!
c
         do i = 1, num_grnd
            grnode = grnd(i)
            do j = 1, num_grnd
               link = grnd(j)
               if (link .eq. grnode) cycle
               do k = 1, graph(grnode)%num_ext_links
                  if ( graph(grnode)%ext_links(k) .ne. link ) cycle
                  start = graph(grnode)%link(k)%num_dofs
                  if ( start .lt. 0 ) then
                     write (*,*) '>> fatal error in fill grp links',start
                     call die_abort
                  endif
                  do dof = 1, 3
                     graph(grnode)%link(k)%dofs(
     &                    graph(grnode)%link(k)%num_dofs + dof) = 
     &                    dstmap(node) + dof - 1
                  enddo
                  graph(grnode)%link(k)%num_dofs = 
     &                 graph(grnode)%link(k)%num_dofs + 3
               enddo
            enddo
         enddo
c     
c
      enddo
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine allo_grnode                  *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 10/19/98                   *
c     *                                                              *
c     *     this subroutine allocates the data structure for a       *
c     *     graph node.                                              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine allo_grnode (num_blks, num_proc_share, num_int_links,
     &     num_ext_links, grnode)
      use mpi_lnpcg
      implicit integer (a-z)
$add common.main
      type (graph_type) :: grnode
c
c         allocate list of blocks in graph node
c     
      allocate (grnode%blks(num_blks), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
c         allocate list of procs which share graph node
c     
      allocate (grnode%proc_share(num_proc_share), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
c         allocate list of links between graph node and the other graph nodes
c     
      allocate (grnode%int_links(num_int_links), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
      allocate (grnode%ext_links(num_ext_links), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
c         indicate that the data structure for each link has not been 
c         allocated yet
c
      grnode%allo_link = .false.
c
      return
 9900 format('>>> FATAL ERROR: memory allocate failure...')
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine allo_link                    *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 10/19/98                   *
c     *                                                              *
c     *     this subroutine allocates the data structure for         *
c     *     each link to a graph node.  The data structure contains  *
c     *     info about the degrees of freedom shared between the     *
c     *     current graph node and the linked graph node.            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine allo_link (num_links, num_dofs, grnode)
      use mpi_lnpcg
      implicit integer (a-z)
$add common.main
      type (graph_type) :: grnode
      dimension num_dofs(*)
c
c         allocate structure of links in the graph node structure
c     
      allocate (grnode%link(num_links), stat = alloc_stat)
      if ( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
c         for each link, allocate arrays to hold the dofs shared between 
c         the local graph node and its links.
c
      do lnk = 1, num_links
         grnode%link(lnk)%num_dofs = num_dofs(lnk) 
         allocate (grnode%link(lnk)%dofs(num_dofs(lnk)), 
     &        stat = alloc_stat)
         if ( alloc_stat .ne. 0 ) then
            write(out,9900)
            call die_abort
            stop
         end if
      enddo
c
c         indicate that the data structure for each link has now been 
c         allocated
c
      grnode%allo_link = .true.
c
c
      return
 9900 format('>>> FATAL ERROR: memory allocate failure...')
      end
c
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine deallocate_owner             *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/04/98                   *
c     *                                                              *
c     *     this subroutine calculates the node owner structure      *
c     *     for domain decomposition portions of the MPI code        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine deallo_grnode(grnode)
      use mpi_lnpcg
      implicit integer (a-z)
      type (graph_type) :: grnode
c
c         allocate pointer from local dofs to global dofs
c
      do i = 1, grnode%num_ext_links
         deallocate (grnode%link(i)%dofs)
      enddo
      deallocate (grnode%link)
      deallocate (grnode%int_links)
      deallocate (grnode%ext_links)
      deallocate (grnode%proc_share)
      deallocate (grnode%blks)
c
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine copy_grnode                  *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 10/19/98                   *
c     *                                                              *
c     *     this subroutine copies the data structure for a          *
c     *     graph node.                                              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine copy_grnode (old, new)
      use mpi_lnpcg
      implicit integer (a-z)
      type (graph_type) :: old, new
c
c
      new%id = old%id
      new%owner = old%owner
      new%order = old%order
      new%num_blks = old%num_blks
      do i = 1, new%num_blks
         new%blks(i) = old%blks(i)
      enddo
      new%num_elem = old%num_elem
      new%num_proc_share = old%num_proc_share
      do i = 1, new%num_proc_share
         new%proc_share(i) = old%proc_share(i)
      enddo
      new%num_int_links = old%num_int_links
      do i = 1, new%num_int_links
         new%int_links(i) = old%int_links(i)
      enddo
      new%num_ext_links = old%num_ext_links
      do i = 1, new%num_ext_links
         new%ext_links(i) = old%ext_links(i)
      enddo
c     
      new%allo_link = old%allo_link
      if (new%allo_link) then
         do lnk = 1, new%num_ext_links
            new%link(lnk)%id = old%link(lnk)%id
            new%link(lnk)%owner = old%link(lnk)%owner
            new%link(lnk)%order = old%link(lnk)%order
            new%link(lnk)%num_dofs = old%link(lnk)%num_dofs
            do i = 1, new%link(lnk)%num_dofs
               new%link(lnk)%dofs(i) = old%link(lnk)%dofs(i)
            enddo
         enddo
      endif
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine print_grnode                 *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 10/19/98                   *
c     *                                                              *
c     *     this subroutine prints the data structure for a          *
c     *     graph node.                                              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine print_grnode (grnode)
      use mpi_lnpcg
      implicit integer (a-z)
$add common.main
      type (graph_type) :: grnode
c
c
      write (out,*) '   -> output of graph node:'
      write (out,*) '    == id:',grnode%id
      write (out,*) '    == owner:',grnode%owner
      write (out,*) '    == order:',grnode%order
      write (out,*) '    == num blocks:',grnode%num_blks
      write (out,*) '    == blocks:'
      write (out,'(10x,10i6)') (grnode%blks(j),j=1,grnode%num_blks)
      write (out,*) '    == num elem:',grnode%num_elem
      write (out,*) '    == num shared procs:',grnode%num_proc_share
      write (out,*) '    == shared procs:'
      write (out,'(10x,10i6)') (grnode%proc_share(j),j=1,
     &     grnode%num_proc_share)
      write (out,*) '    == num internal links:',grnode%num_int_links
      write (out,*) '    == internal links:'
      write (out,'(10x,10i6)') (grnode%int_links(j),
     &     j=1,grnode%num_int_links)
      write (out,*) '    == num ext links:',grnode%num_ext_links
      if ( grnode%allo_link) then
         write (out,*) '    == ext link data:'
         do i=1, grnode%num_ext_links
            write (out,*) '      -- ext link:',grnode%ext_links(i)
            write (out,*) '        -- id:',grnode%link(i)%id
            write (out,*) '        -- owner:',grnode%link(i)%owner
            write (out,*) '        -- order:',grnode%link(i)%order
            write (out,*) '        -- num_dofs:',grnode%link(i)%num_dofs
            if ( grnode%link(i)%num_dofs .gt. 0) then
               write (out,*) '        -- dofs:'
               write (out,'(15x,10i6)') (grnode%link(i)%dofs(j),j=1,
     &              grnode%link(i)%num_dofs)
            endif
         enddo
      else
         write (out,*) '    == external links:'
         write (out,'(10x,10i6)') (grnode%ext_links(j),
     &        j=1,grnode%num_ext_links)
      endif
      write (out,*) '   -> done with graph node:'
c         
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine dd_send                      *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 11/20/98                   *
c     *                                                              *
c     *     this subroutine controls the sending of boundary data    *
c     *     between processors for the ebe preconditioner            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine dd_send ( vec, send_done, num_sends, grnode, oper )
      use mpi_lnpcg
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      logical debug, send_done
      data debug /.false./
#sgl      real
#dbl      double precision 
     &     vec(*), tmp(mxdof)
      dimension send_req(mxnmbl), send_ind(mxnmbl)
      save send_req, tmp, send_ind
c
      dimension status_array(MPI_STATUS_SIZE,mxnmbl)
c
      if (debug)write (*,*) myid,':>>>>> in dd_send, operation:',oper
      goto (100, 200, 300, 400) oper
c
c           set up all the sends for this step
c
  100 continue 
c
c              for this graph node, loop over all the external links and
c              set up sends for the corresponding data.  allocates a
c              temporary array for all the sends
c   
      if ( num_sends .eq. 0) send_ind(1) = 1
      do lnk = 1, local_graph(grnode)%num_ext_links
         num_sends = num_sends + 1
         num_vals = local_graph(grnode)%link(lnk)%num_dofs
	 send_ind(num_sends+1) = send_ind(num_sends) + num_vals
c
c
         do i = 1, num_vals
            tmp(send_ind(num_sends)+i-1) = 
     &           vec(local_graph(grnode)%link(lnk)%dofs(i))
c
         enddo
c  
         call MPI_ISEND (tmp(send_ind(num_sends)),num_vals,MPI_VAL,
     &        local_graph(grnode)%link(lnk)%owner,
     &        local_graph(grnode)%id,
     &        MPI_COMM_WORLD, send_req(num_sends),
     &        ierr)
c
      enddo
c
      goto 9999
c
c         OPERATION = 2 -- test for send
c
  200 continue
      call MPI_testall ( num_sends, send_req, send_done, status_array,
     &      ierr )
      goto 9999
c
c         OPERATION = 3 -- wait for sends
c
  300 continue
      call MPI_waitall ( num_sends, send_req, status_array, ierr)
      goto 9999
c
c         OPERATION = 4 -- clean up after sends (currnetly nothing to do)
c
  400 continue
c
      goto 9999
c
c
 9999 continue
      if (debug)write (*,*) myid,':<<<< leaving dd_send'
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine dd_recv                      *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 11/20/98                   *
c     *                                                              *
c     *     this subroutine controls the recving of boundary data    *
c     *     between processors for the ebe preconditioner            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine dd_recv ( group, vec, recv_done, oper )
      use mpi_lnpcg
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      logical debug, recv_done
      data debug /.false./
#sgl      real
#dbl      double precision 
     &     vec(*), tmp(mxdof), last_val
      dimension recv_req(mxnmbl), recv_ind(mxnmbl)
      integer num_recvs, oper, group
      save num_recvs, recv_req, recv_ind, tmp
c
      dimension status_array(MPI_STATUS_SIZE,mxnmbl)
c
c
      if (debug)write (*,*) myid,':>>>>> in dd_recv, operation:',oper
      goto (100, 200, 300, 400) oper
c
c           set up all the receives for this step
c
  100 continue 
c
c              loop over all the local graph nodes -- for all those that
c              depend on data from the nodes processed during this step,
c              set up a recv for that data. allocate a temporary vectore
c	       to get that information
c   
      num_recvs = 0
      recv_ind(1) = 1
      do grnode = 1, num_graph
c
c                       loop over links for this graph node -- make
c                       a receive for each link
c
         do lnk = 1, local_graph(grnode)%num_ext_links
c 
           if (local_graph(grnode)%link(lnk)%order .ne. group)
     &           cycle
            num_recvs = num_recvs + 1
            num_vals = local_graph(grnode)%link(lnk)%num_dofs
            recv_ind(num_recvs+1) = recv_ind(num_recvs) + num_vals 
c
            call MPI_IRECV (tmp(recv_ind(num_recvs)), num_vals,
     &           MPI_VAL, local_graph(grnode)%link(lnk)%owner,
     &           local_graph(grnode)%link(lnk)%id, MPI_COMM_WORLD,
     &           recv_req(num_recvs), ierr)
         enddo
      enddo
c
      goto 9999
c
c         OPERATION = 2 -- test for recv
c
  200 continue
      call MPI_testall ( num_recvs, recv_req, recv_done, status_array,
     &      ierr )
      goto 9999
c
c         OPERATION = 3 -- wait for recvs
c
  300 continue
      call MPI_waitall ( num_recvs, recv_req, status_array, ierr)
      goto 9999
c
c         OPERATION = 4 -- process recieved data
c
  400 continue
c
c            loop over all the received information -- put into the
c            correct place in the actual vector
c
      num = 0 
      do grnode = 1, num_graph
c
c                       loop over links for this graph node -- make
c                       a receive for each link
c
         do lnk = 1, local_graph(grnode)%num_ext_links
            if (local_graph(grnode)%link(lnk)%order .ne. group)
     &           cycle
            num = num + 1
            num_vals = local_graph(grnode)%link(lnk)%num_dofs
            do i = 1, num_vals
               vec(local_graph(grnode)%link(lnk)%dofs(i)) =
     &              tmp(recv_ind(num)+i-1)
            enddo
         enddo
      enddo
c 
      goto 9999
c
c
 9999 continue
      if (debug)write (*,*) myid,':<<<< leaving dd_recv'
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_graph_send              *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/04/98                   *
c     *                                                              *
c     *     this subroutine sends the graph node owner structures    *
c     *     used in the ebe preconditioner to the processors         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_graph_send(tot_num_graph, graph, max_elems)
      use mpi_lnpcg
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      dimension num_grnodes(0:max_procs-1), num_dofs(mxconn),
     &    blklen(mxdof), disp(mxdof), max_elems(mxnmbl)
      logical debug, found_break
      data debug /.false./
      double precision test(mxdof)
c
      dimension status(MPI_STATUS_SIZE)
c
      type (graph_type) :: graph(*)
c
      if ( debug ) write (out,*) myid,': ->> in wmpi_graph_send'
c 
c            notify slaves
c
      call wmpi_alert_slaves ( 31 )
c
c            For the root processor:
c
      if ( root_processor ) then
c
c              first find out how many graph structures to send to 
c              each proc
c
         do i= 0, numprocs - 1 
            num_grnodes(i) = 0
         enddo
         do grnode = 1, tot_num_graph
            num_grnodes(graph(grnode)%owner) =  
     &           num_grnodes(graph(grnode)%owner) + 1
         enddo
         num_graph = num_grnodes(0)
c
c              allocate the local graph data structures for root
c
         allocate (local_graph(num_graph), stat = alloc_stat)
         if ( alloc_stat .ne. 0 ) then
            write(out,9900)
            call die_abort
            stop
         end if
c
c              copy root processor info into local_graph structure
c
         ptr = 0
         do i = 1, tot_num_graph
            if ( graph(i)%owner .ne. 0 ) cycle
            ptr = ptr + 1
            if ( ptr .gt. num_graph) then
               write (out,*) '>> Fatal error during sending of graphs'
               call die_abort
            endif
            call allo_grnode (graph(i)%num_blks, 
     &           graph(i)%num_proc_share, graph(i)%num_int_links,
     &           graph(i)%num_ext_links, local_graph(ptr))
            do j = 1, graph(i)%num_ext_links
               num_dofs(j) = graph(i)%link(j)%num_dofs
            enddo
            call allo_link (graph(i)%num_ext_links, num_dofs, 
     &           local_graph(ptr))
            call copy_grnode (graph(i),local_graph(ptr))
         enddo
c
c              now loop over the processors and send each proc the information
c              they need. 
c     
         do proc = 1, numprocs - 1
c
c              send out # of graph nodes to allocate on processor, and
c              number of groups found in graph coloring scheme
c
            call MPI_SEND (num_grnodes(proc),1,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
            call MPI_SEND (num_groups,1,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
c
c              loop over graph nodes, skip if not owned by recieving proc
c
            do grnode = 1, tot_num_graph
               if ( graph(grnode)%owner .ne. proc ) cycle
c
c                 send out contstants
c
               call MPI_SEND (graph(grnode)%id,1,MPI_INTEGER,
     &              proc,proc,MPI_COMM_WORLD,ierr)
               call MPI_SEND (graph(grnode)%owner,1,MPI_INTEGER,
     &              proc,proc,MPI_COMM_WORLD,ierr)
               call MPI_SEND (graph(grnode)%order,1,MPI_INTEGER,proc,
     &              proc,MPI_COMM_WORLD,ierr)
               call MPI_SEND (graph(grnode)%num_blks,1,MPI_INTEGER,proc,
     &              proc,MPI_COMM_WORLD,ierr)
               call MPI_SEND (graph(grnode)%num_elem,1,MPI_INTEGER,proc,
     &              proc,MPI_COMM_WORLD,ierr)
               call MPI_SEND (graph(grnode)%num_proc_share,1,
     &              MPI_INTEGER,proc,proc,MPI_COMM_WORLD,ierr)
               call MPI_SEND (graph(grnode)%num_int_links,1,MPI_INTEGER,
     &              proc,proc,MPI_COMM_WORLD,ierr)
               call MPI_SEND (graph(grnode)%num_ext_links,1,MPI_INTEGER,
     &              proc,proc,MPI_COMM_WORLD,ierr)
c     
c                 now send out the arrays
c
               call MPI_SEND (graph(grnode)%blks,
     &              graph(grnode)%num_blks,
     &              MPI_INTEGER,proc,proc,MPI_COMM_WORLD,ierr)
               call MPI_SEND (graph(grnode)%proc_share,  
     &              graph(grnode)%num_proc_share,        
     &              MPI_INTEGER,proc,proc,MPI_COMM_WORLD,ierr)
               call MPI_SEND (graph(grnode)%ext_links,   
     &              graph(grnode)%num_ext_links,         
     &              MPI_INTEGER,proc,proc,MPI_COMM_WORLD,ierr)
               call MPI_SEND (graph(grnode)%int_links,   
     &              graph(grnode)%num_int_links,         
     &              MPI_INTEGER,proc,proc,MPI_COMM_WORLD,ierr)
c
c                 send out links alloc data for this node
c
               do lnk = 1, graph(grnode)%num_ext_links
                  num_dofs(lnk) = graph(grnode)%link(lnk)%num_dofs
               enddo
               call MPI_SEND (num_dofs,graph(grnode)%num_ext_links,
     &              MPI_INTEGER,proc,proc,MPI_COMM_WORLD,ierr)
c
c                 send out links info for this node
c
               do lnk = 1, graph(grnode)%num_ext_links
                  call MPI_SEND (graph(grnode)%link(lnk)%id,1,
     &                 MPI_INTEGER,proc,proc,MPI_COMM_WORLD,ierr)
                  call MPI_SEND (graph(grnode)%link(lnk)%owner,1,
     &                 MPI_INTEGER,proc,proc,MPI_COMM_WORLD,ierr)
                  call MPI_SEND (graph(grnode)%link(lnk)%order,1,
     &                 MPI_INTEGER,proc,proc,MPI_COMM_WORLD,ierr)
                  call MPI_SEND (graph(grnode)%link(lnk)%dofs,
     &                 graph(grnode)%link(lnk)%num_dofs,
     &                 MPI_INTEGER,proc,proc,MPI_COMM_WORLD,ierr)
               enddo
c     
            enddo
c
c              send out the info about ebe ordering so the procs can
c              do stuff
c
            call MPI_SEND (max_elems,num_groups,MPI_INTEGER,proc,
     &           45,MPI_COMM_WORLD,ierr)
c
         enddo
c
      else
c
c           this is for the slave processors
c
c              get the number of graph nodes it needs to allocate, and
c              the number of coloring groups
c
         call MPI_RECV (num_graph,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,
     &        status, ierr)
         call MPI_RECV (num_groups,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,
     &        status, ierr)
c
c              allocate the graph data structures
c
         allocate (local_graph(num_graph), stat = alloc_stat)
         if ( alloc_stat .ne. 0 ) then
            write(out,9900)
            call die_abort
            stop
         end if
c
c              loop over graph nodes
c     
         do i=1, num_graph
c
c                 get constants needed for allocation
c
            call MPI_RECV (id,1,MPI_INTEGER,
     &           0,myid,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV (owner,1,MPI_INTEGER,
     &           0,myid,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV (order,1,MPI_INTEGER,
     &           0,myid,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV (num_blks,1,MPI_INTEGER,
     &           0,myid,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV (num_elem,1,MPI_INTEGER,
     &           0,myid,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV (num_proc_share,1,
     &           MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV (num_int_links,1,
     &           MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV (num_ext_links,1,
     &           MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
c
c                 allocate graph structures
c
            call allo_grnode (num_blks, num_proc_share, num_int_links,
     &           num_ext_links, local_graph(i))
c
c                 fill graph structures
c
            local_graph(i)%id = id
            local_graph(i)%owner = owner
            local_graph(i)%order = order
            local_graph(i)%num_blks = num_blks
            local_graph(i)%num_elem = num_elem
            local_graph(i)%num_proc_share = num_proc_share 
            local_graph(i)%num_int_links  = num_int_links  
            local_graph(i)%num_ext_links  = num_ext_links  
            call MPI_RECV (local_graph(i)%blks,num_blks,MPI_INTEGER,0,
     &           myid,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV (local_graph(i)%proc_share,num_proc_share,
     &           MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV (local_graph(i)%ext_links,num_ext_links,
     &           MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV (local_graph(i)%int_links,num_int_links,
     &           MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
c
c                 loop over external links -- get num dofs for each link
c
            call MPI_RECV (num_dofs,num_ext_links,MPI_INTEGER,0,
     &           myid,MPI_COMM_WORLD,status,ierr)
            call allo_link (num_ext_links, num_dofs, local_graph(i))
c
c                 now fill link data structures
c
            do lnk = 1, num_ext_links
               call MPI_RECV (local_graph(i)%link(lnk)%id,1,
     &              MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
               call MPI_RECV (local_graph(i)%link(lnk)%owner,1,
     &              MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
               call MPI_RECV (local_graph(i)%link(lnk)%order,1,
     &              MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
               local_graph(i)%link(lnk)%num_dofs = num_dofs(lnk)
               call MPI_RECV (local_graph(i)%link(lnk)%dofs,
     &              num_dofs(lnk),MPI_INTEGER,0,
     &              myid,MPI_COMM_WORLD,status,ierr)
            enddo
c
         enddo
c
c              send out the info about ebe ordering so the procs can
c              do stuff
c
         call MPI_RECV (max_elems,num_groups,MPI_INTEGER,0,45,
     &        MPI_COMM_WORLD,status,ierr)
c
c
      endif
c
c           print out the data structures if debug is on.
c
      if (debug) then
         write (*,*) '>>>> graphs in global dof ordering'
         do proc = 0, numprocs - 1
            if (myid .eq. proc) then
               write (out,*) '>=>=>=>=>=> This is proc', myid
               do i = 1, num_graph
                  write (*,*) '   -> for graph ',i
                  call print_grnode(local_graph(i))
               enddo
            endif
            call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         enddo
      endif
c
c
c           now have everyone change dof numberings from global to local
c
      do grnode = 1, num_graph
         do lnk = 1, local_graph(grnode)%num_ext_links
            do dof = 1, local_graph(grnode)%link(lnk)%num_dofs
               local_graph(grnode)%link(lnk)%dofs(dof) = 
     &              local_nodes%global2local(
     &              local_graph(grnode)%link(lnk)%dofs(dof))
            enddo
         enddo
      enddo
c
c           print out the data structures if debug is on.
c
      if (debug) then
         write (*,*) '>>>> graphs in local dof ordering'
         do proc = 0, numprocs - 1
            if (myid .ne. proc) then
	       call wmpi_wait
	       cycle
            endif
            write (out,*) '>=>=>=>=>=> This is proc', myid
            do i = 1, num_graph
               write (*,*) '   -> for graph ',i
               call print_grnode(local_graph(i))
            enddo
            call wmpi_wait
         enddo
      endif
c
c
c         have each processor calculcate the ordering for application of
c         the ebe preconditioner
c
      call find_ebe_order ( max_elems )
c
      if ( debug ) write (out,*) myid,': <<- leaving wmpi_graph_send'
c
      call wmpi_wait
c
      return
 9900 format('>>> FATAL ERROR: memory allocate failure...')
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine find_ebe_order               *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 12/19/98                   *
c     *                                                              *
c     *     this subroutine determines an ordering for the steps in  *
c     *     the ebe preconditioner which will operate in as close    *
c     *     to full parallel as possible. The ordering is stored     *
c     *     in the array order_list.                                 *
c     *                                                              *
c     ****************************************************************
c
c
c          format for ordering list:
c             type  ( order_list(1,#) ) 	value ( order_list(2,#) )
c             ----   				-----
c               1 : compute internal blocks      # of internal blocks
c		2 : compute external group       # of group of external
c		    of blocks (graph node)	    blocks to compute
c		3 : new parallel group 		 parallel group #
c		    indicator for ebksub 
c		4 : new parallel group 		 parallel group #
c		    indicator for efrwrd
c
c          typical ordering for 2 parallel groups:
c
c               3   1 -- seperator for parallel group 1
c               1   3 --   do 3 internal blocks
c               4   1 -- seperator for parallel group 1
c               2   1 --   do boundary group #1
c               2   3 --   do boundary group #3
c               3   2 -- seperator for parallel group 2
c               1   2 --   do 2 internal blocks
c               4   2 -- seperator for parallel group 2
c               2   2 --   do boundary group #2
c               1   3 --   do 1 internal block to balance work
c               3   3 -- seperator for parallel group 3  
c               1   4 --   do 4 internal blocks
c               4   3 -- seperator for parallel group 3
c               
c
      subroutine find_ebe_order ( max_elems )
      use mpi_lnpcg
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      dimension max_elems(mxnmbl), num_ext(mxnmbl), num_need(mxnmbl),
     &     num_max (mxnmbl)
      real fraction
      logical debug, need_filled
      data debug /.false./
c
      call wmpi_wait
      if ( debug ) write (out,*) myid,':   ->> in find_ebe_order'
c
c                 if there is only one processor, build simple scheduling
c                 list   
c 
      if ( numprocs .eq. 1 ) then
c
         num_order = 3
         if ( allocated (order_list)) then
            write (out,9000) 1
            call die_abort
         endif
c     
         allocate (order_list(2,num_order), stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
            write(out,9000) 2
            call die_abort
            stop
         end if 
c
         order_list(1,1) = 3
         order_list(2,1) = 1
         order_list(1,2) = 1
         order_list(2,2) = nelblk
         order_list(1,3) = 4
         order_list(2,3) = 1
         order_ptr = 3
         goto 1000
      endif
c
c           count the number of boundary blocks that are to be computed
c           during each parallel group.  Store this in num_ext.
c
      tot_int = 0
      tot_ext = 0
      do grp = 1, num_groups
         num_ext(grp) = 0
      enddo
c     
      do grnode = 1, num_graph
         grp = local_graph(grnode)%order
         num_ext(grp) = num_ext(grp) + 
     &        local_graph(grnode)%num_elem
         tot_ext = tot_ext + local_graph(grnode)%num_elem
      enddo
c     
      do i = 1, local_nodes%num_int_blks
         blk = local_nodes%internal_blks(i)
         tot_int = tot_int + elblks(0,blk)
      enddo
c     
c                compute the total number of internal elements needed to
c                achieve load balance.  
c
      num_need(1) = max_elems(1) - num_ext(1) 
      num_need(num_groups + 1) = max_elems(num_groups) - 
     &     num_ext(num_groups) 
      tot_need = num_need(1) + num_need( num_groups + 1 )
      num_max(1) = max_elems(1)
      num_max(num_groups + 1) = max_elems(num_groups)
      tot_max = num_max(1) + num_max(num_groups +1)
      do grp = 2, num_groups
         old_need = max_elems(grp-1) - num_ext(grp-1)
         new_need = max_elems(grp) - num_ext(grp)
         if ( old_need .gt. new_need ) then          
            num_need(grp) = old_need
            num_max(grp) = max_elems(grp - 1)
         else
            num_need(grp) = new_need
            num_max(grp) = max_elems(grp)
         endif
         tot_need = tot_need + num_need(grp)
         tot_max = tot_max + num_max(grp)
      enddo
c     
c           allocate the order_list to hold the scheduling list for 
c           processes during application of the ebe preconditioner
c
      num_order = (num_groups + 1)*3 + num_graph
      if ( allocated (order_list)) then
         write (out,9000) 1
         call die_abort
      endif
c     
      allocate (order_list(2,num_order), stat = alloc_stat )
      if ( alloc_stat .ne. 0 ) then
         write(out,9000) 2
         call die_abort
         stop
      end if 
c     
c           **----------------**
c           Create ordering list
c           **----------------**
c
c              Case 1) -- we don't have enough internal blocks to achieve 
c              -------    good load balance.  
c
      if ( tot_int .lt. tot_need ) then
         if ( debug) write (out,*) myid,':>>>>>>> lousy load balance'
c
c                 build scheduling list. Loop
c                 over the groups, assigning internal blocks to aleviate 
c                 load balancing problems.  Since we don't have enough
c                 internal elements, we will not be able to fully
c                 optimize the load balance.
c     
         ptr = 0
         entry = 0
c
c                    calculate the fraction of internal elements
c                    to number of elements needed for load balance.
c                    we will aim for that fraction for this group
c
         fraction = float( tot_int + tot_ext ) / float(tot_max)
c
c                    loop over groups
c
         do grp = 1, num_groups
c
c                    indicate start of internal block list to 
c                    compute
c
            entry = entry + 1
            order_list(1,entry) = 3
            order_list(2,entry) = grp
c
c                    init vars for this group
c
            need_filled = .false.
            need_num = 0
            extra_num = 0
            use_int = 0
            num_int = 0
c
c                    count the number of internal blocks required
c                    to get within "fraction" of the need
c
            do while (.true.)
c
c                       check to see if we are within 90% of "fraction" of the 
c                       need -- if so, exit loop 
c     
               if ( num_int + num_ext(grp) .ge. 
     &              .90 * fraction * num_max(grp) ) then
                  exit
               endif
c
c                       increment pointer to internal block to use. If
c                       we have exceeded the number of internal blocks
c                       available, then skip the sceduling of internal
c                       blocks.
c
               ptr = ptr + 1
               if ( ptr .gt. local_nodes%num_int_blks ) then
                  exit
               endif
               blk = local_nodes%internal_blks(ptr)
               use_int = use_int + 1
               num_int = num_int + elblks(0,blk)
c     
            enddo
c
c                    write scheduling entry for the computed internal
c                    blocks
c
            if ( use_int .gt. 0 ) then
               entry = entry + 1
               order_list(1,entry) = 1
               order_list(2,entry) = use_int
            endif
c
c                    indicate beginning of list of external graph nodes
c                    to compute
c
            entry = entry + 1
            order_list(1,entry) = 4
            order_list(2,entry) = grp
c
c                   loop over graph nodes -- find the ones which are 
c                   processed at this point and make entries in scheduling
c                   list
c
            do grnode = 1, num_graph
               if ( local_graph(grnode)%order .eq. grp) then
                  entry = entry + 1
                  order_list(1,entry) = 2
                  order_list(2,entry) = grnode
               endif   
            enddo
c     
         enddo                                     
c
c                 now finish off any remaining internal blocks
c
         entry = entry + 1
         order_list(1,entry) = 3
         order_list(2,entry) = num_groups + 1
c     
         entry = entry + 1
         use_int = local_nodes%num_int_blks - ptr
         if ( use_int .gt. 0 ) then
            order_list(1,entry) = 1
            order_list(2,entry) = use_int
         endif
c     
         entry = entry + 1
         order_list(1,entry) = 4
         order_list(2,entry) = num_groups + 1
c     
         order_ptr = entry
c     
c
      else
c
c              Case 2) -- we have at least enough internal elements to achieve 
c              -------    good load balance.  Take remainder and distribute 
c                         evenly over the internal block calculations
c
c                 compute the number of extra blocks to add to each one
c
         if ( debug) write (out,*) myid,':>>>>>>> great load balance!'
c
c                 build scheduling list. Loop
c                 over the groups, assigning internal blocks to aleviate 
c                 load balancing problems.
c
         ptr = 0
         entry = 0
         do grp = 1, num_groups
c
c                    indicate start of internal block list to 
c                    compute
c
            entry = entry + 1
            order_list(1,entry) = 3
            order_list(2,entry) = grp
c
c                    init vars for this group
c
            need_filled = .false.
            need_num = 0
            extra_num = 0
            use_int = 0
c
c                    calc the number of extra internal elements
c                    per remaining group.  We will try to 
c                    schedule this many extra elements for this 
c                    internal group. 
c
            extra = ( tot_int - tot_need ) / (num_groups - grp + 2)
            if ( extra .le. 0 .and. num_need(grp) .eq. 0) goto 100
c
c                    count the number of internal blocks requiured
c                    to satisfy the deficit requirements and the
c                    extra elements
c
            do while (.true.)
c
c                       increment pointer to internal block to use. If
c                       we have exceeded the number of internal blocks
c                       available, then skip the sceduling of internal
c                       blocks.
c
               ptr = ptr + 1
               if ( ptr .gt. local_nodes%num_int_blks ) then
                  exit
               endif
               blk = local_nodes%internal_blks(ptr)
               use_int = use_int + 1
               tot_int = tot_int - elblks(0,blk)
c    
c                       schedule the internal blocks required to 
c                       fill the needed amount
c
               if ( .not. need_filled) then                    
                  need_num = need_num + elblks(0,blk)
                  if ( need_num .ge. num_need(grp) ) then
                     need_filled = .true.
                     extra_num = need_num - num_need(grp)
                     tot_need = tot_need - num_need(grp) 
                  endif  
               else
c
c                       scedule internal blocks for extra elements
c
                  extra_num = extra_num + elblks(0,blk)
               endif
c
c                       if the number of extra assigned elements is
c                       within 90% of our traget, then stop asssigning
c                       internal blocks for this group
c
               if ( extra_num .ge. .90 * extra ) then
                  exit
               endif
            enddo
c
c                    write scheduling entry for the computed internal
c                    blocks
c
            if ( use_int .gt. 0 ) then
               entry = entry + 1
               order_list(1,entry) = 1
               order_list(2,entry) = use_int
            endif
c
c                    indicate beginning of list of external graph nodes
c                    to compute
c
 100        continue
            entry = entry + 1
            order_list(1,entry) = 4
            order_list(2,entry) = grp
c     
c                   loop over graph nodes -- find the ones which are 
c                   processed at this point and make entries in scheduling
c                   list
c
            do grnode = 1, num_graph
               if ( local_graph(grnode)%order .eq. grp) then
                  entry = entry + 1
                  order_list(1,entry) = 2
                  order_list(2,entry) = grnode
               endif   
            enddo
c     
         enddo                                     
c
c                 now finish off any remaining internal blocks
c
         entry = entry + 1
         order_list(1,entry) = 3
         order_list(2,entry) = num_groups + 1
c     
         entry = entry + 1
         use_int = local_nodes%num_int_blks - ptr
         if ( use_int .gt. 0 ) then
            order_list(1,entry) = 1
            order_list(2,entry) = use_int
         endif
c     
         entry = entry + 1
         order_list(1,entry) = 4
         order_list(2,entry) = num_groups + 1
c     
         order_ptr = entry
c
      endif
c     
c     
 1000 continue
c
c           print out the data structures if debug is on.
c
      if (debug) then
	 call wmpi_wait
         write (out,*) '>>>> ebe procedure ordering'
         do proc = 0, numprocs - 1
            if (myid .ne. proc) then
	       call wmpi_wait
	       cycle
            endif
            write (out,*) '>=>=>=>=>=> This is proc', myid
            tot = 0
            do i = 1, order_ptr
	       type = order_list(1,i)
	       if (type .eq. 3 ) then
		  write (out,'("-> 3 : ", i7)') order_list(2,i)
	       else if (type .eq. 4) then
		  write (out,'("<- 4 : ", i7)') order_list(2,i)
               else if (type .eq. 1) then
		  write (out,'(" * 1 : ", i7)') order_list(2,i)
                  tot = tot + order_list(2,i)
               else if (type .eq. 2) then
		  write (out,'("   2 : ", i7," -> #bl:",i7)') 
     &                 order_list(2,i), 
     &                 local_graph(order_list(2,i))%num_blks
                  tot = tot + local_graph(order_list(2,i))%num_blks
	       endif	     
            enddo
            write (out,*) '   -------> total blocks:', tot
            call wmpi_wait
         enddo
      endif
c
      call wmpi_wait
c
      if ( debug ) write (out,*) myid,':   <<- leave find_ebe_order'
c
 9999 continue
 9000 format (">>> Fatal Allocation error in find_ebe_order, type:",i2)
c
      return
      end

