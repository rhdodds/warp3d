c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_init                    *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                             *
c     *                   last modified : 12/17/2017 rhd             *
c     *                                                              *
c     *     This routine initializes MPI on all processors at the    *
c     *     beginning of a MPI warp run.  It also creates the        *
c     *     MPI_VAL datatype which is equivalent to either a real or *
c     *     double precision variable depending on whether the       *
c     *     platform is double precision or single precision. This   *
c     *     greatly simplifies later MPI calls without requiring     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_init
      use global_data ! old common.main
      use hypre_parameters
      implicit none
      include "mpif.h"
c
      integer :: ierr, iallowed, idum
      logical :: ldum1, ldum2
      logical, parameter :: local_debug = .false.
      character(len=1) :: dums
      real :: dumr
      double precision :: dumd
c
c              initialize MPI, find the total number of processors
c              for this run, and get the id number of the local
c              processor in the total processor numbering. Processor 0
c              will be the root(manager) processor, while all others
c              be workers.
c
c      call MPI_INIT_THREAD (MPI_THREAD_FUNNELED, iallowed, ierr)
c      if( iallowed .ne. MPI_THREAD_FUNNELED ) then
c        write(out,*) '>>> fatal error on MPI start up.'
c        write(out,*) '    MPI_THREAD_FUNNELED, iallowed:',
c     &                  MPI_THREAD_FUNNELED, iallowed
c        call die_abort
c      end if
      call MPI_INIT_THREAD (MPI_THREAD_MULTIPLE, iallowed, ierr)
      if( iallowed .ne. MPI_THREAD_MULTIPLE ) then
        write(out,*) '>>> fatal error on MPI start up.'
        write(out,*) '    MPI_THREAD_MULTIPLE, iallowed:',
     &                  MPI_THREAD_FUNNELED, iallowed
        call die_abort
      end if
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      if( local_debug )
     &  write(out,'(" >>> numprocs:",i4," myid:", i4)') numprocs, myid
      if( numprocs .gt. max_procs ) then
         call errmsg (312, idum, dums, dumr, dumd)
         call MPI_FINALIZE (ierr)
         call die_abort
         stop
      endif
c
c             create MPI_VAL datatype to be equal to either MPI_REAL
c             MPI_VAL is stored in common.main, while MPI_REAL and
c             MPI_DOUBLE_PRECISION is in mpif.h.
c
c             MPI_TYPE_CONTIGUOUS defines MPI_VAL to be equivalent to
c             either MPI_REAL or MPI_DOUBLE_PRECISION,
c             while MPI_TYPE_COMMIT officially registers the datatype
c             with the MPI routines.
c
      call MPI_TYPE_CONTIGUOUS(1,MPI_DOUBLE_PRECISION,MPI_VAL,ierr)
      call MPI_TYPE_COMMIT (MPI_VAL,ierr)
c
c             set the logical flags root_processor and worker_processor
c             for simple identification of whether we are the worker
c             or the master.
c
      if( myid .eq. 0 ) then
         root_processor   = .true.
         slave_processor  = .false.
         worker_processor = .false.
      else
         root_processor   = .false.
         slave_processor  = .true.
         worker_processor = .true.
      endif
c
c             set use_mpi flag to inidcate that this is the MPI version
c
      use_mpi = .true.
c
c             if we are the root processor, return to driver, and
c             continue executing the code as normal, reading input, etc.
c
      if( root_processor ) return
c
c             processors executing code past this point are worker
c             processors. Call initst to initialize variables they
c             may need, then throw them in the worker handler.  In this
c             routine they wait for instructions from the root processor
c             about what job to do. The worker processors never return
c             from the worker handler; they die when told to do so by
c             the root processor.
c
      call initst( ldum1, ldum2 )
c
      call wmpi_handle_slaves
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_procids                 *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 2/17/2017 rhd              *
c     *                                                              *
c     *     This routine has all the worker processors find their    *
c     *     process ids, then send these Ids back to the root.       *
c     *     This allows root to suspend the MPI threads if a         *
c     *     threaded sparse solver is called.                        *
c     *                                                              *
c     *       ==> deprecated by RHD on 6/17/2017 <==                 *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_procids
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: proc, ierr, status(MPI_STATUS_SIZE), procid
      integer, external :: getpid
      logical, parameter :: debug = .false.

c
c          if this is the the root processor, then recieve data from each
c          processor.
c
      if( root_processor ) then
c
        do proc = 1, numprocs -1
        call MPI_RECV ( proc_pids(proc), 1, MPI_INTEGER,
     &          proc, 1, MPI_COMM_WORLD, status, ierr)
        end do
c
         if( debug ) then
           write (out,*) '>>> Here are proc ids gathered on root:'
           do proc = 1, numprocs -1
             write (out,*) '    proc:', proc, ' id:',proc_pids(proc)
           end do
           procid = getpid ()
           if( ierr .ne. 0 ) then
             write (out,*) '>>> FATAL ERROR in getting process id.'
             call die_abort
           endif
         write(out,*) '>> the pid of root is ',procid
         endif
c
      else    ! worker
c
        procid = getpid ()
        call MPI_SEND( procid, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD,
     &                 ierr)
c
      endif
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_suspend                 *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 2/17/2017 rhd              *
c     *                                                              *
c     *     This routine allows the root to suspend the worker       *
c     *     processes when a threaded sparse solver is called.       *
c     *     Once the solver completes, this routine can also         *
c     *     restart the mpi processes.                               *
c     *                                                              *
c     *       ==> deprecated by RHD on 6/17/2017 <==                 *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_suspend( option )
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: option
c
      integer :: proc, ierr
      integer, external :: kill
      logical, parameter :: debug = .false.
      logical, save :: suspended
      data suspended /.false./
c
c	             suspend the slave mpi processes using system 'kill'
c              commands and hardware-specific signal numbers.
c
      if( option .eq. 1 ) then
         do proc = 1, numprocs - 1
          if( debug ) write (out,*) '>>> root suspend ',proc_pids(proc)
          ierr = kill (proc_pids(proc), 19)
          if( ierr .ne. 0 ) then
            write (out,*) '>>> FATAL ERROR: root unable to suspend',
     &                  ' pid ',proc_pids(proc)
            call die_abort
           endif
         end do
         suspended = .true.
c
      else if( option .eq. 2 .and. suspended ) then
c
c              awaken the suspended mpi processes
c
         do proc = 1, numprocs - 1
          if( debug ) write (out,*) '>>> root reviving ',proc_pids(proc)
          ierr = kill (proc_pids(proc), 18)
          if( ierr .ne. 0 ) then
           write(out,*) '>>> FATAL ERROR: root unable to awaken',
     &                ' pid ',proc_pids(proc)
               call die_abort
          end if
         end do
         suspended = .false.
      end if
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_alert_slaves            *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/17/2017 rhd             *
c     *                                                              *
c     *  The MPI implementation of warp3d follows a manager-worker   *
c     *  approach for most of the code.  The manager processor       *
c     *  (processor 0) runs warp3d in full, notifying the worker     *
c     *  processors when they have work to do.                       *
c     *                                                              *
c     *  This routine allows the root processor to contact the       *
c     *  worker processes and tell them what to do next.             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_alert_slaves( do_it )
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: do_it
c
      integer :: ierr
      logical, parameter :: debug = .false.
c
c              return for workers - workers are not
c              allowed to tell other workers what to do.
c
      if( worker_processor ) return
c
c              do_it is the code used by handle_mpi_workers to
c              tell the workers what to do.
c
      if( debug ) write (out,'("====> Root is sending alert ",i3)')
     &     do_it
      call MPI_BCAST( do_it, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_reduce_vec              *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 12/20/2017 rhd             *
c     *                                                              *
c     *     This routine takes processor local vectors which hold    *
c     *     contributions to a global vector, and sums them on       *
c     *     the root processor.  Thus once this routine is complete, *
c     *     the root processor has a full copy of the vector, as if  *
c     *     it had calculated the whole vector on its own.           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_reduce_vec( vec, size )
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: size
      double precision :: vec(size)

      call wmpi_reduce_vec_inplace( vec, size )
      return
c
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine wmpi_allreduce_dble              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/29/2018 rhd              *
c     *                                                              *
c     *     takes processor local vectors which hold                 *
c     *     contributions to a global vector, and sums them on       *
c     *     the root processor - then broadcasts to all ranks        *
c     *                                                              *
c     *     this version uses the MPI_IN_PLACE option to eliminate   *
c     *     the usual allocation/copy of temporary on root           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_allreduce_dble( vec, size )
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: size
      double precision :: vec(size)
c
      integer ::  ierr
c
      call MPI_ALLREDUCE( MPI_IN_PLACE, vec, size, MPI_DOUBLE,
     &                    MPI_SUM, MPI_COMM_WORLD, ierr )
      return
c
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine wmpi_allreduce_real              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/29/2018 rhd              *
c     *                                                              *
c     *     takes processor local vectors which hold                 *
c     *     contributions to a global vector, and sums them on       *
c     *     the root processor - then broadcasts to all ranks        *
c     *                                                              *
c     *     this version uses the MPI_IN_PLACE option to eliminate   *
c     *     the usual allocation/copy of temporary on root           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_allreduce_real( vec, size )
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: size
      real :: vec(size)
c
      integer ::  ierr

      call MPI_ALLREDUCE( MPI_IN_PLACE, vec, size, MPI_REAL,
     &                    MPI_SUM, MPI_COMM_WORLD, ierr )
      return
c
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine wmpi_allreduce_int               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/29/2018 rhd              *
c     *                                                              *
c     *     takes processor local vectors which hold                 *
c     *     contributions to a global vector, and sums them on       *
c     *     the root processor - then broadcasts to all ranks        *
c     *                                                              *
c     *     this version uses the MPI_IN_PLACE option to eliminate   *
c     *     the usual allocation/copy of temporary on root           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_allreduce_int( vec, size )
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: size, vec(size)
      integer ::  ierr
c
      call MPI_ALLREDUCE( MPI_IN_PLACE, vec, size, MPI_INTEGER,
     &                    MPI_SUM, MPI_COMM_WORLD, ierr )
      return
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_reduce_vec_inplace      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/13/2017 rhd             *
c     *                                                              *
c     *     takes processor local vectors which hold                 *
c     *     contributions to a global vector, and sums them on       *
c     *     the root processor.  Thus once this routine is complete, *
c     *     the root processor has a full copy of the vector, as if  *
c     *     it had calculated the whole vector on its own.           *
c     *                                                              *
c     *     this version uses the MPI_IN_PLACE option to eliminate   *
c     *     the usual allocation/copy of temporary on root           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_reduce_vec_inplace( vec, size )
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: size
      double precision :: vec(size)
c
      integer ::  ierr
c
      if( root_processor ) then
        call MPI_REDUCE( MPI_IN_PLACE, vec, size, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr )
      else
        call MPI_REDUCE( vec, vec, size, MPI_DOUBLE_PRECISION,
     &                   MPI_SUM, 0, MPI_COMM_WORLD, ierr )
      end if
c
      return
c
      end

c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_reduce_vec_std          *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/17/2017 rhd             *
c     *                                                              *
c     *     This routine takes processor local vectors which hold    *
c     *     contributions to a global vector, and sums them on       *
c     *     the root processor.  Thus once this routine is complete, *
c     *     the root processor has a full copy of the vector, as if  *
c     *     it had calculated the whole vector on its own.           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_reduce_vec_std( vec, size )
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: size
      double precision :: vec(size)
c
      integer :: alloc_stat, ierr
      double precision, allocatable :: vec_tmp(:)
      double precision, parameter :: zero = 0.0d00
c
c              we cannot reduce into the same vector we are sending, so
c              allocate a temporary array to hold the data
c
      allocate ( vec_tmp(size), stat = alloc_stat )
      if ( alloc_stat .ne. 0) then
         write(out,9000)
         call die_abort
      endif
c
c              copy processor local vector into temporary vector, zero
c              local
c
      vec_tmp = vec
      vec = zero
      call MPI_REDUCE( vec_tmp, vec, size, MPI_DOUBLE_PRECISION,
     &                 MPI_SUM, 0, MPI_COMM_WORLD, ierr )
      deallocate ( vec_tmp )
c
      return
 9000 format ('>>>  FATAL ERROR:  could not allocate additional',/,
     &        '>>>       memory in wmpi_reduce_vec.')

      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine wmpi_allreduce_vec_log          *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 4/21/2018 rhd              *
c     *                                                              *
c     *     reduce LOR worker vectors to root then bcast to all      *
c     *                                                              *
c     ****************************************************************
      subroutine wmpi_allreduce_vec_log( vec, size )
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: size
      logical :: vec(size)
c
      integer :: alloc_stat, ierr
c
      call MPI_ALLREDUCE( MPI_IN_PLACE, vec, size, MPI_LOGICAL,
     &                    MPI_LOR, MPI_COMM_WORLD, ierr )
c
      return
c
      end
c
c
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine wmpi_reduce_vec_new              *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 6/17/2017 rhd              *
c     *                                                              *
c     *     This routine takes processor local vectors which hold    *
c     *     contributions to a global vector, and sums them on       *
c     *     the root processor.  Thus once this routine is complete, *
c     *     the root processor has a full copy of the vector, as if  *
c     *     it had calculated the whole vector on its own.           *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_reduce_vec_new( vec )
      use global_data ! old common.main
      use mpi_lnpcg
      implicit none
      include "mpif.h"
c
      double precision :: vec(*)
c
      integer :: status(MPI_STATUS_SIZE), proc, ierr, i, num_dof
      double precision, allocatable :: vec_tmp(:)
      double precision, parameter :: zero = 0.0d00


      allocate( vec_tmp(nodof) )
      vec_tmp = zero
c
c          if this is the the root processor, then recieve data
c          from each processor.
c
      if( root_processor ) then
c
        do proc = 1, numprocs -1
         call MPI_RECV( vec_tmp, procdof2glob(proc)%num_dof,
     &                  MPI_VAL, proc, 1, MPI_COMM_WORLD, status, ierr)
            do i = 1, procdof2glob(proc)%num_dof
               vec(procdof2glob(proc)%dof(i)) =
     &             vec(procdof2glob(proc)%dof(i)) + vec_tmp(i)
            end do
        end do
c
      else  ! worker process
c
        num_dof = local_nodes%num_local_nodes * mxndof
        do i = 1, num_dof
            vec_tmp(i) = vec(local_nodes%local2global(i))
        end do
        call MPI_SEND( vec_tmp, num_dof, MPI_VAL, 0, 1,
     &        MPI_COMM_WORLD, ierr )
c
      endif
c
      return
 9000 format ('>>>  FATAL ERROR:  could not allocate additional',/,
     &        '>>>       memory in wmpi_reduce_vec.')
c
      end
c
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_red_intvec              *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/17/2017 rhd             *
c     *                                                              *
c     *     This routine takes processor local vectors which hold    *
c     *     contributions to a global vector, and sums them on       *
c     *     the root processor.  Thus once this routine is complete, *
c     *     the root processor has a full copy of the vector, as if  *
c     *     it had calculated the whole vector on its own.           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_red_intvec( vec, size )
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: size, vec(size)
      integer, allocatable ::  vec_tmp(:)
c
      integer :: alloc_stat, ierr
c
c              we cannot reduce into the same vector we are sending, so
c              allocate a temporary array to hold the data
c
      allocate( vec_tmp(size), stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
         write (out,9000)
         call die_abort
      endif
c
      vec_tmp = vec
      vec     = 0
c
c              reduce the vector using the MPI_REDUCE command.  This
c              reduces the result only to the root processor.
c
      call MPI_REDUCE( vec_tmp, vec, size,MPI_INTEGER, MPI_SUM, 0,
     &                 MPI_COMM_WORLD, ierr )
      deallocate( vec_tmp )
c
      return
 9000 format ('>>>  FATAL ERROR:  could not allocate additional',/,
     &        '>>>       memory in wmpi_red_intvec.')

      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_wait                    *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/05/98                   *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_wait
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: ierr
c
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_dotprod                 *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/18/2017 rhd             *
c     *                                                              *
c     *     This routine drives the reduction via mpi of a           *
c     *     dot-product distributed among processors. Each           *
c     *     processor computes their contribution to the dot product *
c     *     as a scalar, then the scalars are summed together and    *
c     *     the result given to all processors.                      *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_dotprod( scalar )
      implicit none
      include "mpif.h"
c
      double precision :: scalar, scalar_tmp, ierr
c
      scalar_tmp = scalar
      call MPI_ALLREDUCE( scalar_tmp, scalar, 1, MPI_DOUBLE_PRECISION,
     &                    MPI_SUM, MPI_COMM_WORLD, ierr )
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_redint                  *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/18/2017 rhd             *
c     *                                                              *
c     *     This routine drives the reduction via mpi of a           *
c     *     integer which has a contribution on each processor.      *
c     *     The total is only available on the root processor.       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_redint( scalar )
      implicit none
      include "mpif.h"
c
      integer :: scalar_tmp, scalar, ierr

      scalar_tmp = scalar
      call MPI_REDUCE( scalar_tmp, scalar, 1, MPI_INTEGER,
     &                 MPI_SUM, 0, MPI_COMM_WORLD, ierr )
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_redlog                  *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/18/2017 rhd             *
c     *                                                              *
c     *     This routine drives the reduction via mpi of a           *
c     *     logical which has a contribution on each processor.      *
c     *     The result of an OR on all the logicals is only         *
c     *     available on the root processor.                         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_redlog( log_var )
      implicit none
      include "mpif.h"
c
      integer :: ierr
      logical :: log_var, log_var_tmp
c
      log_var_tmp = log_var
      call MPI_ALLREDUCE( log_var_tmp, log_var, 1, MPI_LOGICAL,
     &                    MPI_LOR, MPI_COMM_WORLD, ierr )
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_bcast_int               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/18/2017 rhd             *
c     *                                                              *
c     *      Broadcast an integer from root to all the worker        *
c     *      processors.                                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_bcast_int ( int_var )
      implicit none
      include "mpif.h"
c
      integer :: int_var
      integer :: ierr
c
      call MPI_BCAST( int_var, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_bcast_real              *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/19/2017 rhd             *
c     *                                                              *
c     *      Broadcast a real from root to all the worker            *
c     *      processors.                                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_bcast_real ( real_var )
      implicit none
      include "mpif.h"
c
      real :: real_var
c
      integer :: ierr
c
      call MPI_BCAST(real_var,1,MPI_REAL,0,MPI_COMM_WORLD, ierr)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_bcast_log               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 2/19/2017 rhd              *
c     *                                                              *
c     *      Broadcast a logical variable from root to all the       *
c     *      workers.                                                *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_bcast_log ( log_var )
      implicit none
      include "mpif.h"
c
      logical :: log_var
c
      integer :: ierr
c
      call MPI_BCAST( log_var, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_bcast_string            *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/19/2017 rhd             *
c     *                                                              *
c     *      Broadcast a charater string to workers                  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_bcast_string ( string, nchars )
      implicit none
      include "mpif.h"
c
      integer :: nchars
      character(len=nchars) :: string
c
      integer :: ierr
c
      call MPI_BCAST( string, nchars, MPI_CHARACTER, 0, MPI_COMM_WORLD,
     &               ierr )
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_do_external_db          *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 2/19/2017 rhd              *
c     *                                                              *
c     *          invoke  UEXTERNALDB Abaqus support routine and      *
c     *          other specific UEXTERNADB routines                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_do_uexternaldb
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: dimension status (MPI_STATUS_SIZE), ierr
      integer :: aba_lop, aba_lrestart, aba_kstep, aba_kinc
      double precision :: aba_time(2), aba_dtime
      double precision, parameter :: zero = 0.0d0
      logical, parameter :: local_debug = .false.
c
c         tell workers we are about to run uexternadb
c
      call wmpi_alert_slaves ( 51 )
c
      if( local_debug )
     &  write (out,'("=> proc ",i3," run externaldb")')myid
c
c         broadcast several solution parameters to be passed to
c         uexternaldb on each worker
c
      call wmpi_wait
      call MPI_BCAST(total_model_time,1,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dt,1,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ltmstp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(douextdb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
      select case ( douextdb )
       case ( 1 )   ! from incomp start of step 1
         aba_lop = 0
         aba_lrestart = 0
         aba_time(1) = zero; aba_time(2) = zero
         aba_dtime = dt
         aba_kstep = 1
         aba_kinc = 1
       case ( 2 )  ! from driver (normal termination)
        aba_lop = 3
        aba_lrestart = 0
        aba_time(1) = total_model_time
        aba_time(2) = total_model_time
        aba_dtime = dt
        aba_kstep = 1
        aba_kinc  = ltmstp   ! current WARP3D load step
       case ( 3 )  ! from mnralg (start of step, restart adaptive step)
        aba_lop = 1
        aba_lrestart = 0
        aba_time(1) = total_model_time
        aba_time(2) = total_model_time
        aba_dtime = dt
        aba_kstep = 1
        aba_kinc  = ltmstp + 1  ! current WARP3D load step
       case ( 4 ) ! from mnralg (end step or adaptive substep)
        aba_lop = 2
        aba_lrestart = 0
        aba_time(1) = total_model_time
        aba_time(2) = total_model_time
        aba_dtime = dt
        aba_kstep = 1
        aba_kinc  = ltmstp + 1  ! current WARP3D load step
       case ( 5 )   ! from store (make restart file)
        aba_lop = 2
        aba_lrestart = 1
        aba_time(1) = total_model_time
        aba_time(2) = total_model_time
        aba_dtime = dt
        aba_kstep = 1
        aba_kinc  = ltmstp   ! current WARP3D load step
       case ( 6 )   ! from reopen (restarting analysis)
        aba_lop = 4
        aba_lrestart = 0
        aba_time(1) = total_model_time
        aba_time(2) = total_model_time
        aba_dtime = dt
        aba_kstep = 1
        aba_kinc  = ltmstp + 1  ! next WARP3D load step to be solved
       case default
           write(out,9000)
           call die_abort
       end select
c
      call uexternaldb( aba_lop, aba_lrestart, aba_time,
     &                  aba_dtime, aba_kstep, aba_kinc )
c
      call uexternaldb_mm04_cavity( aba_lop, aba_lrestart, aba_time,
     &                  aba_dtime, aba_kstep, aba_kinc )
      call wmpi_wait
      return
c
 9000 format(/,2x,"FATAL ERROR: invalid douextdb in ",
     & "wmpi_do_uexternaldb. ",/,2x,"Analysis terminated...",//)
c
       end

c
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_basic              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/15/2018 rhd              *
c     *                                                              *
c     *       This subroutine allows the root processor to send      *
c     *       basic model data, such as the coordinates, incidences, *
c     *       etc., to all of workers.  This data is                 *
c     *       only sent once during the analysis.                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_basic
      use global_data ! old common.main
c
      use main_data, only: incmap, incid, crdmap, dtemp_nodes,
     &      dtemp_elems, temper_nodes, temper_elems, invdst,
     &      temper_nodes_ref, temperatures_ref, fgm_node_values,
     &      fgm_node_values_defined, fgm_node_values_cols,
     &      fgm_node_values_used,
     &      matprp, lmtprp, imatprp, dmatprp, smatprp,
     &      nonlocal_analysis, asymmetric_assembly,
     &      initial_stresses_user_routine, initial_stresses_file,
     &      initial_stresses, initial_state_option, initial_state_step,
     &      initial_stresses_input

      use elem_block_data, only: edest_blocks, cdest_blocks,
     &                           edest_blk_list, cdest_blk_list
      use segmental_curves
      use contact, only : use_contact
      use mm10_defs, only : one_crystal_hist_size, common_hist_size

c
      implicit none
      include "mpif.h"
c
      integer :: status(MPI_STATUS_SIZE), ierr
      integer, allocatable :: element_node_counts(:)
      double precision, parameter :: zero = 0.0d0
      logical :: initial_stresses_exist
c
c         tell workers we are about to send them basic data
c
      call wmpi_alert_slaves ( 2 )
c
c      write (out,'("=> proc ",i3," is doing basic data transfer")')myid
c
c         broadcast element properties and element blocking info. These
c         values are fixed once the analysis starts.
c
c             constants:
c
      call MPI_BCAST(noelem,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nonode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nelblk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nodof,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(beta_fact,1,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(inctop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(max_current_pts,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_BCAST(max_current_curves,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_BCAST(num_seg_curve_sets,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_BCAST(lgnmcn,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(temperatures_ref,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_BCAST(fgm_node_values_defined,1,MPI_LOGICAL,0,
     &               MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fgm_node_values_used,1,MPI_LOGICAL,0,
     &               MPI_COMM_WORLD,ierr)
      call MPI_BCAST (use_contact,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nonlocal_analysis,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_BCAST(asymmetric_assembly,1,MPI_LOGICAL,0,
     &      MPI_COMM_WORLD, ierr)
      call MPI_BCAST(max_mpc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(max_mpc_tied,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      if( nonlocal_analysis .and. myid .eq. 0 ) then
         write(out,9100)
         call die_abort
      end if
      call MPI_BCAST( one_crystal_hist_size, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( common_hist_size, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( common_hist_size, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( initial_stresses_user_routine, 1, MPI_LOGICAL,
     &                0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast( initial_stresses_file, 100, MPI_CHARACTER, 0,
     &                MPI_COMM_WORLD,ierr)
      call MPI_BCAST( initial_state_option, 1, MPI_LOGICAL,
     &                0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST( initial_state_step, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( initial_stresses_input, 1, MPI_LOGICAL, 0,
     &                MPI_COMM_WORLD, ierr )
c
c             static arrays:
c
      call MPI_BCAST(props,mxelpr*noelem,MPI_REAL,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_BCAST(cp,mxedof,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(icp,mxutsz*2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dcp,mxedof,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(elblks,4*nelblk,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(num_seg_points,max_seg_curves,MPI_INTEGER,0,
     &               MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curves,max_seg_curves*max_seg_points*2,
     &               MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curves_min_stress,max_seg_curves,
     &               MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curve_def,max_seg_curves,MPI_LOGICAL,0,
     &               MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curves_type,max_seg_curves,MPI_INTEGER,0,
     &               MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curve_table,
     &            (max_seg_curves+1)*max_seg_curve_sets, MPI_INTEGER,0,
     &            MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curves_value,max_seg_curves,
     &               MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curves_ym,max_seg_curves,
     &               MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curves_nu,max_seg_curves,
     &               MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curves_alpha,max_seg_curves,
     &               MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curves_gp_sigma_0,max_seg_curves,
     &               MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curves_gp_h_u,max_seg_curves,
     &               MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curves_gp_beta_u,max_seg_curves,
     &               MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curves_gp_delta_u,max_seg_curves,
     &               MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(matprp,mxmtpr*mxmat,MPI_REAL,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(imatprp,mxmtpr*mxmat,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_Bcast(dmatprp,mxmtpr*mxmat,MPI_DOUBLE_PRECISION,0,
     &               MPI_COMM_WORLD,ierr)
      call MPI_Bcast(smatprp,mxmtpr*mxmat*24,MPI_CHARACTER,0,
     &               MPI_COMM_WORLD,ierr)
c
c                  allocated arrays:
c
c                  incidence data structures, inverse incidence
c                  dta structures, inverse dof maps
c
      if( worker_processor ) then
         call mem_allocate( 4 )  ! u, v, a, ... zeroed
         call mem_allocate(9)
         call init_maps ( 0, 1 )
         call init_maps ( 0, 2 )
      endif
      call MPI_BCAST(c,nodof,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dstmap,nonode,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST( incid, inctop, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( incmap, noelem, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( invdst, nodof, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &                ierr )
c
c            build a vector of no. of nodes on all elements and broadcast.
c            this makes the process of building inverse incidences and
c            inverse dof mappings no dependent on slaves having the
c            full properties table (for the future). setup_slave routine
c            builds these data structures for each slave rather than
c            sending complex, pointer-based array structures.
c
      allocate( element_node_counts(noelem) )
      if( root_processor ) element_node_counts(1:noelem) =
     &                     iprops(2,1:noelem)
      call MPI_BCAST( element_node_counts, noelem, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      if( worker_processor ) call setup_slave( element_node_counts )
      deallocate( element_node_counts )
c
c                  functionally graded material properties at
c                  model nodes
c
      if( fgm_node_values_defined ) then
        if(  worker_processor ) call mem_allocate(20)
        call MPI_BCAST( fgm_node_values, nonode*fgm_node_values_cols,
     &                 MPI_REAL, 0, MPI_COMM_WORLD, ierr )
      end if
c
c                  crdmap
c
      if( worker_processor ) call mem_allocate(14)
      call MPI_BCAST( crdmap, nonode, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &                ierr)
c
c                  temperature loadings
c
      if( worker_processor ) then
         call mem_allocate(1)
         call mem_allocate(2)
         call mem_allocate(16)
      endif
      call MPI_BCAST( dtemp_nodes, nonode, MPI_VAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( dtemp_elems, noelem, MPI_VAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( temper_nodes, nonode, MPI_VAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( temper_nodes_ref, nonode, MPI_VAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( temper_elems, noelem, MPI_VAL, 0, MPI_COMM_WORLD,
     &                ierr)
c
c                  elems_to_blocks structure
c
      if( worker_processor ) call init_eblock_map
c
c                  allocate space for mdiag
c
      if( worker_processor ) call mem_allocate ( 12 )
c
c             set up data distribution information so all processors
c             know what blocks of elements they own.
c
      call wmpi_calc_dist
c
c             initial stresses if they exist
c
      if( root_processor ) initial_stresses_exist =
     &                   allocated( initial_stresses )
      call MPI_BCAST( initial_stresses_exist, 1, MPI_LOGICAL,
     &                0, MPI_COMM_WORLD, ierr)
      if( initial_stresses_exist ) then
        if( worker_processor ) allocate( initial_stresses(6,noelem) )
        call MPI_Bcast( initial_stresses, 6*noelem,
     &                  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,
     &                  ierr )
      end if

c
c             now communicate allocated blocked structures:
c
c                call the standard initialization routines for these
c                structures, so that each will contain the proper starting
c                information.  Note that the root processor has full
c                copys of each of these arrays; however, the workers
c                only store the information pertaining to the elements they
c                own.  The routines mostly allocate blocks.
c
      if( worker_processor ) then
         call cdest_init
         call edest_init
         call history_cep_init( 0, 1 )
         call rotation_init( 0, 1 )
         call strains_init( 0, 1 )
         call stresses_init( 0, 1 )
         call element_volumes_init( 0, 1 )
         call mem_allocate( 24 )
         if( use_contact ) call mem_allocate( 25 )
         call allocate_damage( 12 )
      end if
c
      return
c
 9100 format(">>>> Fatal Error: nonlocal cohesive modeling not",
     & /,    "                  supported with MPI. ",
     & "Analysis terminated.",//)
c
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_analysis           *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 6/29/20198 rhd             *
c     *                                                              *
c     *       send data from anaylsis parameters to all the MPI      *
c     *       processors                                             *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_analysis
      use global_data ! old common.main
      use main_data,only : umat_serial, cp_matls_present,
     &                     cp_unloading, creep_model_used,
     &                     initial_state_option, initial_state_step
      implicit none
      include "mpif.h"
c
      integer :: ierr
c
c         tell slaves we are about to send them data about
c         the analysis parameters
c
      call wmpi_alert_slaves ( 9 )
c
c      write (out,'("=> proc ",i3," is doing analysis data transfer")')
c     &         myid
c
c         broadcast element properties and element blocking info
c
c         values from analysis parameters
c
      call MPI_BCAST( eps_bbar, 1, MPI_VAL, 0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( qbar_flag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( mxlitr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( dt, 1, MPI_VAL, 0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( nbeta, 1, MPI_VAL, 0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( convrg, mxcvtests, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( signal_flag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( adaptive_flag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( umat_serial, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( cp_matls_present, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( cp_unloading, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( creep_model_used, 1, MPI_LOGICAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( initial_state_option, 1, MPI_LOGICAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( initial_state_step, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_const              *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 4/10/2018 rhd              *
c     *                                                              *
c     *       send data about constraints to all the MPI workers     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_const
      use global_data ! old common.main
      use main_data, only: trn, trnmat
      implicit none
      include "mpif.h"
c
      integer :: ierr, node, dum, count
      double precision, allocatable :: buffer(:,:,:)
c
c              Basic data must be sent before we can call this routine.
c              Send simple constraint map.
c              Send logical vector for possible local rotation of
c                  constraints at nodes.
c              If node rotations for constraints exist, pack 3x3 node
c                  matrices into a buffer, bcast the buffer, then
c                  unpack back into 3x3 rotation matrices on each
c                  node as required.
c                  Packing just saves possibly many thousands of
c                  bcasts for 3x3 materices.
c
      call wmpi_alert_slaves ( 10 )
c
c      write (out,'("=> proc ",i3," is doing const data transfer")')myid
c
c                   broadcast values from constraint input
c
      call MPI_BCAST( csthed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
      call MPI_BCAST( cstmap,nodof,MPI_INTEGER,0,MPI_COMM_WORLD,ierr )
      if( worker_processor ) call mem_allocate( 3 ) ! trn, trnmat
      call MPI_BCAST( trn,nonode,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr )
c
      count = 0
      do node = 1, nonode
         if( trn(node) ) count = count + 1
      end do
      if( count .eq. 0 ) return ! no 3x3 rotation mats to send
c
      allocate( buffer(3,3,count) )
      if( root_processor ) then
        count = 0
        do node = 1, nonode
          if( trn(node) ) then
            count = count + 1
            buffer(1:3,1:3,count) = trnmat(node)%mat
          end if
        end do
      end if
c
      call MPI_BCAST( buffer,9*count, MPI_VAL, 0,
     &                MPI_COMM_WORLD, ierr )

      if( worker_processor ) then
       count = 0
       do node = 1, nonode
         if( trn(node) ) then
          count = count + 1
          allocate( trnmat(node)%mat(3,3) )
          trnmat(node)%mat = buffer(1:3,1:3,count)
         end if
       end do
      end if
c
      deallocate( buffer )
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_itern              *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/20/2017 rhd             *
c     *                                                              *
c     *   send workers the information needed for Newton iteration.  *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_itern
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: ierr
      logical, parameter :: local_debug = .false.
c
c                   incremental displacements as determined
c                   from the equation solve.
c
      call wmpi_alert_slaves( 24 )
      call MPI_BCAST( du,nodof,MPI_VAL,0,MPI_COMM_WORLD,ierr )
      call wmpi_wait
c
      if( local_debug ) write(out,9000) myid, norm2( du )
c
      return
c
 9000 format(1x,'>>> L2 norm du after bcast update: ',i4,e14.6)
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_step               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 7/1/2018 rhd               *
c     *                                                              *
c     *       send workers the information needed for                *
c     *       each load step. This includes u, v, a                  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_step
      use global_data ! old common.main
      use main_data, only: asymmetric_assembly, extrapolated_du,
     &                     non_zero_imposed_du, initial_state_option,
     &                     initial_state_step
      use adaptive_steps, only : adapt_disp_fact, adapt_temper_fact,
     &                           adapt_load_fact
      implicit none
      include "mpif.h"
c
      integer :: ierr
      logical, parameter :: local_debug = .false.
c
      call wmpi_alert_slaves ( 12 )
c
      call MPI_BCAST( total_model_time, 1, MPI_VAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( dt, 1, MPI_VAL, 0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( scaling_factor, 1, MPI_VAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( adapt_temper_fact, 1, MPI_VAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( adapt_disp_fact, 1, MPI_VAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( adapt_load_fact, 1, MPI_VAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( asymmetric_assembly, 1, MPI_LOGICAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( extrapolated_du, 1, MPI_LOGICAL, 0,
     &                MPI_COMM_WORLD, ierr)
      call MPI_BCAST( non_zero_imposed_du, 1, MPI_LOGICAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( initial_state_option, 1, MPI_LOGICAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( initial_state_step, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
c
      call MPI_BCAST( u, nodof, MPI_VAL, 0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( v, nodof, MPI_VAL, 0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( a, nodof, MPI_VAL, 0, MPI_COMM_WORLD, ierr )
c
      if( local_debug ) write(out,9000) myid, norm2( u )
c
      return
 9000 format(1x,'.. wmpi_send_step. rank, norm2(u): ', i4, e14.6 )
      end
c
c     ****************************************************************
c     *                                                              *
c     *                subbroutine wmpi_send_temp_eqloads            *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/19/2017 rhd             *
c     *                                                              *
c     *       send workers the temperature data for                  *
c     *       current step. send element equivalent nodal loads      *
c     *       if they exist.                                         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_temp_eqloads
      use global_data ! old common.main
      use main_data, only: dtemp_nodes, dtemp_elems, eq_node_force_len,
     &                     eq_node_force_indexes, eq_node_forces
      implicit none
      include "mpif.h"
c
      integer :: ierr
      double precision :: mag
      double precision, parameter :: zero = 0.0d0
      logical, parameter :: local_debug = .false.
c
      call wmpi_alert_slaves ( 22 )
c
c      write (out,'("=> proc ",i3," is doing step data transfer")')myid
c
c            send logical flag indicating if there are any temperatures
c            to be considered in the current loading.
c
      call MPI_BCAST( temperatures, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     &                ierr )
c
c            send the incremental element and nodal temperatures,
c            if there are any.
c
      if( temperatures ) then
         call MPI_BCAST( dtemp_nodes, nonode, MPI_VAL, 0,
     &                   MPI_COMM_WORLD, ierr )
         call MPI_BCAST( dtemp_elems, noelem, MPI_VAL, 0,
     &                   MPI_COMM_WORLD, ierr )
         if( local_debug ) write(out,9000) myid, norm2(dtemp_nodes),
     &                                     norm2(dtemp_elems)
      end if
c
c            send the global integer flag which sets the length of the
c            packed vectors of element equivalent nodal forces for the
c            step. then send the indexing vector and the packed
c            vector of equiv. force values.
c
      call MPI_BCAST( eq_node_force_len,  1,  MPI_INTEGER,  0,
     &                MPI_COMM_WORLD,  ierr )
      if( eq_node_force_len > 0 ) then
        if( worker_processor ) call mem_allocate( 23 )
        call MPI_BCAST( eq_node_force_indexes,  noelem, MPI_INTEGER, 0,
     &                  MPI_COMM_WORLD,  ierr )
        call MPI_BCAST( eq_node_forces, eq_node_force_len, MPI_VAL, 0,
     &                  MPI_COMM_WORLD, ierr )
        if( local_debug ) write(out,9010) myid,
     &                    sum(eq_node_force_indexes),
     &                    norm2(eq_node_forces(1:eq_node_force_len))
       end if
c
      return
c
 9000 format(1x,'..wmpi_send_temp_eqloads. myid, norm2a, norm2b: ',i4,
     &     e14.6, 2x, e14.6 )
 9010 format(1x,'..wmpi_send_temp_eqloads. myid, sum-c, norm2d: ',i4,
     &     e14.6, 2x, e14.6 )
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine wmpi_send_contact               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 07/5/2017 rhd              *
c     *                                                              *
c     *                 send contact info to all processors          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_contact( restart )
      use global_data ! old common.main
      use contact
      implicit none
c
      include "mpif.h"
c
      integer :: ierr, proc, shape, j, i
      integer :: status(MPI_STATUS_SIZE)
      double precision, parameter :: zero = 0.0d0
      logical :: restart, referenced, ldum
      logical, parameter :: debug = .false.
c
      call wmpi_alert_slaves ( 28 )
c
      if( debug )
     &   write(out,'("=> proc ",i3," do contact data trans")') myid
c
c         broadcast contact information:
c
c             constants:
c
      call MPI_BCAST (use_contact,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c
      if( .not. use_contact) then
         if( debug) write (*,*) myid,': no use contact, baby.'
         return
      endif
c
      call MPI_BCAST(num_contact,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c
      if( debug ) write (*,*) myid,': -> num_contact = ', num_contact
c
c             static arrays:
c
      call MPI_BCAST( cplane_vec, 3*2*num_contact, MPI_VAL, 0,
     &                MPI_COMM_WORLD, ierr)
      call MPI_BCAST( cshape_norm, 3*num_contact, MPI_VAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( cshape_pnt, 3*num_contact, MPI_VAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( cshape_rate, 3*num_contact, MPI_VAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( cshape_param, 10*num_contact, MPI_VAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( contact_depth, num_contact, MPI_VAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( contact_stiff, num_contact, MPI_VAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( contact_fric, num_contact, MPI_VAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( contact_shape, num_contact, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( contact_outside, num_contact, MPI_LOGICAL, 0,
     &                MPI_COMM_WORLD, ierr )
c
c             if we just came back from a restart, we must also
c             broadcast contact cause.  We just broadcast the entire
c             array to all processors, then let them clear out the
c             entries that they don't reference.  This is an extremely
c             clumsy and inefficient broadcast, but it works...
c             little time of basic thought could make this better.
c
      if( restart ) then
         call MPI_BCAST( contact_cause, maxcontact*nonode, MPI_INTEGER,
     &                   0, MPI_COMM_WORLD, ierr )
c
      end if
c
      if ( .not. debug ) return
c
      write (*,*) '>>> Here are the contact surfaces:'
      do proc = 0, numprocs -1
         if (myid .eq. proc ) then
            write (*,*) '**** proc = ',proc,' of ',numprocs,' ****'
            do shape = 1, num_contact
               if ( contact_shape(shape) .eq. 1) then
                  write (*,*) '>>>> Contact Params for shape ',shape
                  write (*,*) '   type: ', contact_shape(shape),
     &                 '   stiff:', contact_stiff(shape),
     &                 '   fric: ', contact_fric(shape),
     &                 '   depth:', contact_depth(shape),
     &                 '   rate: ', (cshape_rate(j,shape),j=1,3)
                  write (*,*) '>  input vecs:'
                  do j = 1, 2
                     write (*, *) '     ', cplane_vec (1,j,shape),
     &                    cplane_vec (2,j,shape),
     &                    cplane_vec (3,j,shape)
                  end do
                  write (*,*) '>  norm vec:'
                  write (*, *) '     ', (cshape_norm (i,shape),i=1,3)
               else if( contact_shape(shape) .eq. 2) then
                  write (*,*) '>>>> Contact Params for shape ',shape
                  write (*,*) '   type: ', contact_shape(shape),
     &                 '   stiff:', contact_stiff(shape),
     &                 '   fric: ', contact_fric(shape),
     &                 '   depth:', contact_depth(shape),
     &                 '   rate: ', (cshape_rate(j,shape),j=1,3)
                  write (*,*) '>  base point:'
                  write (*,*) '     ', (cshape_pnt (i,shape),i=1,3)
                  write (*,*) '>  direction:'
                  write (*,*) '     ', (cshape_norm (i,shape),i=1,3)
                  write (*,*) '>  radius:', cshape_param (1,shape),
     &                 '  length:',cshape_param (2,shape)
               else
                  write (*,*) '>>>> No Contact def for shape ',shape
               endif
            end do
            call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         else
            write (*,*) '> --> proc ',myid,' is waiting. '
            call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         endif
      end do
      write (*,*) '>>> proc ',myid,' is done. '
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_contact_gthr            *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 06/17/2017 rhd             *
c     *                                                              *
c     *        This routine gathers all the contact information      *
c     *        calculated by each processor back to the root         *
c     *        processor. This includes the contact force vector,    *
c     *        the contact-induced transformation matricies, and     *
c     *        the contact cause vector.                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_contact_gthr
      use global_data ! old common.main
      use main_data, only: trn, trnmat
      use mpi_lnpcg, only: local_nodes
      use contact
      implicit none
      include "mpif.h"
c
      integer :: ierr, status(MPI_STATUS_SIZE), i, j, proc, idum, node
      double precision :: mag
      double precision, parameter :: zero = 0.0d0
      logical :: keep_going
      logical, parameter ::  debug = .false.
c
      if( debug )
     &  write (out,'("=> proc ",i3," is doing contact gather")')  myid
c
c               reduce the contact_force vector back to the
c               root processor.
c
      call wmpi_reduce_vec_std( contact_force, nodof )
c
      if( root_processor .and. debug ) then
         write (out,*) '>> NONZERO CONTACT FORCE TERMS ON ROOT'
         do i=1, nodof
            if ( contact_force(i) .ne. zero ) then
               write (*,*) '     ',i,'   ',contact_force(i)
            endif
         end do
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      else
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      endif
c
c               gather transformation matrices
c
      if( root_processor ) then
c
c               the root processor loops over the processors.  For each
c	 	            processor, it recieves the number of a node undergoing
c		            contact, then the transformation matrix and cause of
c		            contact for that node.  This continues until all
c		            contacting nodes are processed, then the root processor
c		            moves on to the next processor.
c
         do proc = 1, numprocs - 1
            if( debug ) write (*,*) '======> root get trn from proc ',
     &           proc
c
            keep_going = .true.
            do while( keep_going )
c
               call MPI_RECV( node, 1, MPI_INTEGER, proc, 1,
     &              MPI_COMM_WORLD, status, ierr )
               if( node .le. 0 ) then
                  keep_going = .false.
               else
                  if( .not. trn (node) ) then
                     call allo_trnmat( node, 1, idum )
                     trn(node) = .true.
                  endif
                  call MPI_RECV( trnmat(node)%mat, 9, MPI_VAL, proc, 2,
     &                           MPI_COMM_WORLD, status, ierr )
                  call MPI_RECV( contact_cause(1,node), num_contact,
     &                           MPI_INTEGER, proc, 3, MPI_COMM_WORLD,
     &                           status, ierr )
c
                  if( debug ) then
                     write (*,*) '    > trnmat:'
                     write (*,'(10x,3e13.6)')
     &                    ((trnmat(node)%mat(i,j),j=1,3),i=1,3)
                  endif
c
                  call trn2all( node, 1 )
               endif
c
            end do ! on do while
c
         end do ! on proc
c
      else
c            Each worker sends the number of a node which is
c	           undergoing contact, then follows this with the
c		         transformation matrix and cause of contact. Continue
c		         until worker sends all the information about contacting
c		         nodes which it owns, then it sends a 0 as a node number
c		         to indicate completion.
c
         if( debug )
     &      write (*,*) '======>',myid,' is sending trn to root'
         do i = 1, local_nodes%num_private
            node = local_nodes%private(i)
            if( contact_cause(1,node) .eq. 0 ) cycle
            if( .not. trn(node) ) cycle
            call MPI_SEND( node, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD,
     &                     ierr )
            call MPI_SEND( trnmat(node)%mat, 9, MPI_VAL, 0, 2,
     &                     MPI_COMM_WORLD, ierr )
            call MPI_SEND( contact_cause(1,node), num_contact,
     &                     MPI_INTEGER, 0, 3, MPI_COMM_WORLD, ierr )
         end do ! on i
         call MPI_SEND( -1, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ierr )
c
      endif
c
      if( debug ) write (*,*) '<<=====',myid,' is done w/ trn transfer'
c
      if( debug ) then
         do proc = 0, numprocs - 1
            if( myid .eq. proc ) then
               write(out,*) '********* contact cause for proc:',proc
               do i = 1, nonode
                  if( contact_cause(1,i) .ne. 0 ) then
                     write(out,*) i,'    ',(contact_cause(j,i),j=1,
     &                    num_contact)
                  endif
               end do ! on i
               call MPI_BARRIER( MPI_COMM_WORLD, ierr )
            else
               call MPI_BARRIER( MPI_COMM_WORLD, ierr )
            endif
         end do ! on proc
      endif
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_growth_init             *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/21/2017 rhd             *
c     *                                                              *
c     *       send initial crack growth data to all the MPI          *
c     *       processors                                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_growth_init
      use global_data ! old common.main
c
      use elem_extinct_data, only: dam_blk_killed
      use damage_data, only : dam_ptr, num_kill_elem, growth_by_kill
c
      implicit none
      include "mpif.h"
c
      integer :: ierr
      logical, save :: been_here = .false.
      logical, parameter :: debug = .false.
c
      if( debug ) write (out,'("=> proc ",i3," do init grow trans")')
     &            myid
c
c           return if we have already done the initialization.
c
      if( been_here ) then
         if( debug ) write (out,'("=> proc ",i3," already init grow")')
     &        myid
         return
      end if
      been_here = .true.
c
      call wmpi_alert_slaves ( 25 )
c
c           send the number of killable elements and the index vector
c           which goes from element number to the entry in the crack
c           growth data strucutres.
c
      call MPI_BCAST( num_kill_elem,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr )
      call MPI_BCAST( dam_ptr,noelem,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr )
      if( worker_processor ) growth_by_kill = .true.
c
c           workers just allocate
c           the other crack growth data structures we need
c
      if( worker_processor ) call allocate_damage ( 1 )
c
      if( debug )
     &   write (out,'("=> proc ",i3," done w/ init grow data")')
     &     myid
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine wmpi_send_growth                *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/21/2017 rhd             *
c     *                                                              *
c     *           send crack growth data to all the MPI worker       *
c     *           processors after root has completed processing     *
c     *           at start of the step                               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_growth ( killed_this_time )
      use global_data ! old common.main
c
      use elem_extinct_data, only: dam_blk_killed, dam_state
      use main_data, only: elems_to_blocks
      use damage_data, only : dam_ptr, growth_by_kill, num_kill_elem
c
      implicit none
      include "mpif.h"
c
      integer :: ierr, elem
      logical :: killed_this_time
      logical, parameter :: debug = .false.
      double precision, parameter :: zero = 0.0d0
c
      call wmpi_alert_slaves ( 13 )
c
      if( debug ) write (out,'("=> proc ",i3," do grwth data trans")')
     &     myid
c
c           send them the state vectors describing what elements
c           and blocks have been killed.
c
      call MPI_BCAST( dam_blk_killed, nelblk, MPI_LOGICAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( dam_state, num_kill_elem, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
c
c           let MPI processors know if an element has been killed
c           this time
c
      call MPI_BCAST( killed_this_time, 1, MPI_LOGICAL, 0,
     &                MPI_COMM_WORLD, ierr )
c
c           if we are the root processor, skip the rest of the routine.
c           also, if no elements were killed this step, return.
c
      if( root_processor ) return
      if( .not. killed_this_time ) return
c
c           if any elements have just been killed, then figure
c           out which ones from dam_state and dam_ptr.  If the
c           killed element is owned by the local processor, then
c           kill the element (set its stresses, strains, poisson's
c           ratio and youngs modulus to zero).
c
      do elem = 1, noelem
         if( dam_ptr(elem) .eq. 0 ) cycle
         if( elblks(2,(elems_to_blocks(elem,1))) .ne. myid ) cycle
         if( dam_state(dam_ptr(elem)) .eq. 1 ) then
            if( debug )
     &          write(out,'("=> proc ",i3,": elem",i5," is dead")')
     &          myid, elem
            call kill_element ( elem, debug )
         end if
      end do
c
      if( debug )
     &     write(out,'("=> proc ",i3," done w/ grth data trans ")')
     &              myid
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                    subroutine wmpi_get_grow                  *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 12/24/99                   *
c     *                                                              *
c     *       Gather data about killable elements that the           *
c     *       workers own and send it back to root so                *
c     *       that root can conduct the growth calculations.         *
c     *       The needed data includes dam_ifv (ifv terms for        *
c     *       killable elements), element stresses, strains, volumes *
c     *       and material history. Only send the information        *
c     *       about killable elements - this dramtically reduces the *
c     *       amount of data that has to be broadcast.               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_get_grow
      use global_data ! old common.main
c
      use main_data, only: elems_to_blocks
      use elem_extinct_data, only: dam_ifv, dam_state, dam_blk_killed
      use elem_block_data, only: history_blocks, urcs_n_blocks,
     &                           eps_n_blocks, urcs_n1_blocks,
     &                           element_vol_blocks, history_blk_list
      use damage_data, only : dam_ptr, growth_by_kill
c
      implicit none
      include "mpif.h"
c
      integer :: ierr, status(MPI_STATUS_SIZE), blk, felem, span, owner,
     &           ngp, histsize, elem, elem_ptr, block_size
      double precision, parameter :: zero = 0.0d0
c
c                   if we aren't using crack growth by element
c                   extinction, then return
c
      if( .not. growth_by_kill ) return
c
c                   tell slaves we need to gather the element information
c                   for killable elements
c
      call wmpi_alert_slaves ( 26 )
c
c                   Loop over the blocks.  For blocks whose
c                   elements are killable and which are owned by a
c                   non-root processor, send the proper element
c                   information back to the root processor.
c
      do blk = 1, nelblk
         felem         = elblks(1,blk)
         span          = elblks(0,blk)
         owner         = elblks(2,blk)
         ngp           = iprops(6,felem)
         histsize      = history_blk_list(blk)
c
c                      this is for root processor.  skip block if
c                      owned by root, or if block has no killable
c                      elements, or if the entire block has been
c                      killed. Otherwise, receive the element
c                      data from the owning slave processor.
c
         if( root_processor ) then
c
            if( owner .eq. 0 ) cycle
            if( iand( iprops(30,felem),2 ) .eq. 0 ) cycle
            if( dam_blk_killed(blk) ) cycle
c
c                         get internal forces for killable elements
c                         not yet killed (in this block)
c
            do elem = felem, felem + span - 1
               elem_ptr = dam_ptr(elem)
               if( dam_state(elem_ptr) .gt. 0 ) cycle
               call MPI_RECV( dam_ifv(1,elem_ptr), mxedof,
     &                        MPI_VAL, owner, owner, MPI_COMM_WORLD,
     &                        status, ierr )
            end do
c
c                         receive stresses, strains, histories
c                         and current volumes for all elements in block
c
            block_size = span * ngp * histsize
            call MPI_RECV( history_blocks(blk)%ptr(1), block_size,
     &                     MPI_VAL, elblks(2,blk), 14, MPI_COMM_WORLD,
     &                     status, ierr )
c
            block_size = span * ngp * nstrs
            call MPI_RECV( urcs_n_blocks(blk)%ptr(1), block_size,
     &                     MPI_VAL, elblks(2,blk), 14, MPI_COMM_WORLD,
     &                     status, ierr )
c
            block_size = span * ngp * nstr
            call MPI_RECV( eps_n_blocks(blk)%ptr(1), block_size,
     &                     MPI_VAL, elblks(2,blk), 14, MPI_COMM_WORLD,
     &                     status, ierr )
c
            block_size = span
            call MPI_RECV( element_vol_blocks(blk)%ptr(1), block_size,
     &                     MPI_VAL, elblks(2,blk), 14, MPI_COMM_WORLD,
     &                     status, ierr )

c
c                      this is for slave processors.  skip block if
c                      not owned by processor, or if block has no killable
c                      elements, or if the entire block has been
c                      killed. Otherwise, send the element
c                      data to the root processor.
c
         else
c
            if( owner .ne. myid ) cycle
            if( iand( iprops(30,felem),2 ) .eq. 0 ) cycle
            if( dam_blk_killed(blk) ) cycle
c
c                         get internal forces for killable elements
c                         not yet killed (in this block)
c
            do elem = felem, felem + span - 1
               elem_ptr = dam_ptr(elem)
               if( dam_state(elem_ptr) .gt. 0 ) cycle
               call MPI_SEND( dam_ifv(1,elem_ptr), mxedof,
     &                        MPI_VAL, 0, myid, MPI_COMM_WORLD, ierr )
            end do
c
c                         receive stresses, strains, histories
c                         and current volumes for all elements in block
c
            block_size = span * ngp * histsize
            call MPI_SEND( history_blocks(blk)%ptr(1), block_size,
     &                     MPI_VAL, 0, 14, MPI_COMM_WORLD, ierr )
c
            block_size    = span * ngp * nstrs
            call MPI_SEND( urcs_n_blocks(blk)%ptr(1), block_size,
     &                     MPI_VAL, 0, 14, MPI_COMM_WORLD, ierr )
c
            block_size = span * ngp * nstr
            call MPI_SEND( eps_n_blocks(blk)%ptr(1), block_size,
     &                     MPI_VAL, 0, 14, MPI_COMM_WORLD, ierr )
c
            block_size = span
            call MPI_SEND( element_vol_blocks(blk)%ptr(1), block_size,
     &                     MPI_VAL, 0, 14, MPI_COMM_WORLD, ierr )
c
         end if
c
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_calc_dist               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 6/5/2017 rhd               *
c     *                                                              *
c     *       compute information about the distribution of          *
c     *       elements, nodes, and degrees of freedom across         *
c     *       the MPI processors.                                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_calc_dist
      use global_data ! old common.main
      use elem_block_data, only: edest_blocks
      implicit none
      include "mpif.h"
c
      integer :: ierr, k,
     &           status(MPI_STATUS_SIZE), blk_ptr(0:max_procs-1),
     &           i, j, blk, felem, num_enodes, num_enode_dof, span,
     &           elem, dof, block, tot_size, num_dof_chunks, proc, idum,
     &           owner_proc
      integer, allocatable :: size_dof_chunks(:), start_dof_chunks(:)
      logical :: open, assign, fail
      logical, allocatable :: access_dof(:)
      logical, parameter :: debug = .false.
      logical, save :: sent = .false.
      real :: dumr
      double precision ::dumd
      character(len=5) :: name, dums
c
      if( debug )
     & write (out,'("=> proc ",i3," is doing dist calcs")')myid
c
      k = 3 * nonode
      allocate( size_dof_chunks(k), start_dof_chunks(k), access_dof(k) )
c
c               make sure we only execute this code once per run.
c
      if( sent ) return
      sent = .true.
      fail = .false.
c
c              first check if the user specified
c              the processor assignment of the blocks in the blocking
c              part of the input file.  If not, do a round robin
c              assignment of the blocks to the different processors.
c
      assign = .false.
      do i = 1, nelblk
         if( elblks(2,i) .lt. 0 ) then
            assign = .true.
            exit
         endif
      enddo
c
      if( assign ) then
c
         if( root_processor ) call errmsg( 313, idum, dums, dumr, dumd )
         do i = 1, nelblk
            elblks(2,i) = mod(i,numprocs)
         end do
c
      else
c
c                 if processor assignment was given, then we must check
c                 to make sure it is valid -- doesnt refer to more
c                 processors than we have, etc.
c
         do i = 1, nelblk
            if( elblks(2,i) .gt. numprocs - 1 ) then
               if( .not. fail .and. root_processor ) then
                  call errmsg( 314, idum, dums, dumr, dumd )
               endif
               fail = .true.
            endif
         end do
c
      endif ! on assign
c
      call MPI_BARRIER( MPI_COMM_WORLD, ierr ) ! wait for everyone
      if( fail ) then
         call MPI_FINALIZE ( ierr )
         stop
      endif
c
c              build index structure onto elblks which points to the
c              next block in sequence which belongs to the processor
c              which owns the current block. Thus if blocks 5 and 8 were
c              owned by processor 2, then elblks(3,5) = 8, because 8 is
c              the next block which belongs to processor 2.
c
      blk_ptr_head(0:numprocs-1) = -1
      blk_ptr(0:numprocs-1) = -1
c
      do blk = 1, nelblk
         owner_proc = elblks(2,blk)
         if( blk_ptr_head(owner_proc) .eq. -1 ) then
            blk_ptr_head(owner_proc) = blk
         else
            elblks(3,blk_ptr(owner_proc)) = blk
         endif
         elblks(3,blk) = -1
         blk_ptr(owner_proc) = blk
      end do
c
c              print out elblks structure if debug
c
      if( debug .and. root_processor ) then
         write (out,*) '>>>>>>  elblks: (span,felem,proc,nextblk)'
         do i = 1, nelblk
            write (out,'(8x,i5,":",4(1x,i7))') i,(elblks(j,i),j=0,3)
         end do
         write (out,*) '>>>>>>  blk_ptr_head: (per proc)'
         write (out,'(8x,12i6)') (blk_ptr_head(i),i=0,numprocs-1)
      endif
c
c              The rest of this subroutine creates datatypes which
c              automatically send the data pertaining to the degrees
c              of freedom effected by the elements owned by the local
c              processor.
c
c             =====> NOTE:  since these data structures are not
c                    directly needed right now, this code is skipped.
c                    However, it may prove to be very useful in the
c                    future, so we have left it in the source.
c
      return
c
c
c              now compute which degrees of freedom are local
c              to the processor.
c
      if( root_processor ) goto 1000
c
      access_dof(1:nodof) = .false.
      do blk = 1, nelblk
c
         if( myid .ne. elblks(2,blk) ) cycle
c
         felem         = elblks(1,blk)
         num_enodes    = iprops(2,felem)
         num_enode_dof = iprops(4,felem)
         span          = elblks(0,blk)
c
         do elem = 1, span
            do i = 1, num_enodes * num_enode_dof
               access_dof ( edest_blocks(blk)%ptr(i,elem) ) = .true.
            end do
         end do
      end do
c
c              Construct a datatype which includes
c              only these local data points. Using this datatype
c              causes broadcasts to only broadcast the actual dof
c              information a processor needs.
c
c                 loop over the degrees of freedom
c
      dof = 0
      block = 0
      tot_size = 0
      open = .false.
c
      do while (.true.)
c
         dof = dof + 1
c
c                    if we are at last dof, close any blocks left
c                    open
c
         if( dof .gt. nodof ) then
            if( open ) then
              size_dof_chunks(block) =
     &              nodof - start_dof_chunks ( block )
              tot_size = tot_size + size_dof_chunks(block)
            endif
            exit
         endif
c
c                    if this dof is needed by the processor, and the
c                    last dof was not, then this is the start of
c                    a chunk.
c
         if( access_dof(dof) ) then
            if( .not. open ) then
               block = block + 1
               start_dof_chunks ( block ) = dof - 1
               open = .true.
            endif
         else
c
c                    this dof is not needed by processor. If last one
c                    was, then set the size of the chunk.
c
            if( open ) then
               size_dof_chunks(block) = dof - start_dof_chunks(block)+1
               tot_size = tot_size + size_dof_chunks(block)
               open = .false.
            endif
         endif
c
      end do
      num_dof_chunks = block
      num_dof_local(myid) = tot_size
c
c       output dof info just so we can check
c
c      write (name,'("dofs",i1)') myid
c      my_unit = 68 + myid
c      open (unit = my_unit, file = name, status = 'unknown')
c      write (my_unit, *) '>>>> dof need list:'
c      do i=1, nodof
c         if (access_dof(i)) then
c            write (my_unit,'(3x,i6,"*")')i
c         else
c            write (my_unit,'(3x,i6)')i
c         endif
c      enddo
c      write (my_unit, *) '>>>> dofs needed for proc:', myid
c      write (my_unit, *) 'start   size'
c      do i=1, num_dof_chunks
c         write (my_unit, *) start_dof_chunks(i), size_dof_chunks(i)
c      enddo
c      close (my_unit)
c
c                    Now create an MPI datatype that corresponds to
c                    this chunked structure.  Using this data type
c                    in a broadcast enables broadcast of only
c                    the data the processor needs.
c
c                    The data type is only used between the root processor
c                    and the local processor, so both the root and local
c                    processor must create and commit the datatype at the
c                    same time.
c
 1000 continue
      if(  root_processor  ) then
         do proc = 1, numprocs - 1
            call MPI_RECV( num_dof_chunks, 1, MPI_INTEGER, proc,
     &           proc, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(num_dof_local(proc), 1, MPI_INTEGER, proc,
     &           proc, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV( start_dof_chunks, num_dof_chunks,
     &           MPI_INTEGER, proc, proc, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV( size_dof_chunks, num_dof_chunks,
     &           MPI_INTEGER, proc, proc, MPI_COMM_WORLD, status, ierr)
            call MPI_TYPE_INDEXED ( num_dof_chunks, size_dof_chunks,
     &           start_dof_chunks, MPI_VAL, MPI_DOF_LOCAL(proc), ierr)
            call MPI_TYPE_COMMIT ( MPI_DOF_LOCAL(proc), ierr )
         end do
      else
         call MPI_SEND( num_dof_chunks, 1, MPI_INTEGER, 0, myid,
     &        MPI_COMM_WORLD, ierr)
         call MPI_SEND( num_dof_local(myid), 1, MPI_INTEGER, 0, myid,
     &        MPI_COMM_WORLD, ierr)
         call MPI_SEND( start_dof_chunks, num_dof_chunks,
     &        MPI_INTEGER, 0, myid, MPI_COMM_WORLD, ierr)
         call MPI_SEND( size_dof_chunks, num_dof_chunks,
     &        MPI_INTEGER, 0, myid, MPI_COMM_WORLD, ierr)
         call MPI_TYPE_INDEXED ( num_dof_chunks, size_dof_chunks,
     &        start_dof_chunks, MPI_VAL, MPI_DOF_LOCAL(myid), ierr)
         call MPI_TYPE_COMMIT ( MPI_DOF_LOCAL(myid), ierr )
      endif
c
      if( debug )
     &  write (out,*) myid,':>> hey, weza levin wmpi_calc_dist'
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine die_gracefully               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 6/17/2017 rhd              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine die_gracefully
      implicit none
      include "mpif.h"
c
      integer :: ierr
c
c              tell workers to quit
c
      call wmpi_alert_slaves ( 1 )
c
c              tell MPI that we are ending now
c
      call MPI_FINALIZE( ierr )
c
      stop
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine die_abort                    *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/21/2017 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine die_abort
      implicit none
      include "mpif.h"
c
      integer :: idum, ierr
c
c              force all processes to quit
c
      call MPI_ABORT( MPI_COMM_WORLD, idum, ierr )
c
      stop
c
      return
      end
c
c
c         we have 3 versions of the routine to bring element stiffness
c         blocks back to rank 0. All 3 work correctly.
c
c             1. loop over all blocks sequentially. use send/recv to
c                request and retrieve each block.
c
c             2. OpenMP parallel over all ranks. each thread loops
c                over all blocks for that rank and brings them back.
c                uses send/recv
c
c             3. Use isend/irecv to post requests for all blocks
c                from workers and handle as they arrive on rank 0.
c
c          Benchmarking as of 12/21/2017 at OSC using DAPL fabric
c          shows that overall method (1) is fastest. Surprising.
c
c          Methods 2,3 can be 50% slower wall time.
c
c          We use method 1 for now.
c
c
c     ****************************************************************
c     *                                                              *
c     *               subroutine wmpi_check_error                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/21/2017  rhd             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine  wmpi_check_error( warp3d_error, myid,
     &                              mpi_error_code, out )
      implicit none
      include "mpif.h"
c
      integer :: warp3d_error, myid, mpi_error_code, out, nchars,
     &           now_mpi_error
      character(len=1024) :: message
c
      if(  mpi_error_code .eq. 0 ) return
c
      call MPI_ERROR_STRING( mpi_error_code,  message, nchars,
     &                       now_mpi_error )

      write(out,9000) warp3d_error, myid, mpi_error_code,
     &                message(1:nchars)
      call die_abort
c
 9000 format(
     &   15x, 'FATAL ERRROR: check on MPI call @ ',i4,' rank: ',i4,
     & /,15x, '              MPI error message string:',
     & /,5x,a,
     & /,15x  '              Job aborted.')
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *           subroutine wmpi_combine_stf_non_blocking           *
c     *                                                              *
c     *                   last modified : 12/21/2017 rhd             *
c     *                                                              *
c     *     gathers the block element stiffnesses from workers.      *
c     *     when done, root has copy of all element stiff blocks.    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_combine_stf
      use global_data ! old common.main
      use elem_block_data, only: estiff_blocks
      use main_data, only: asymmetric_assembly
      implicit none
      include "mpif.h"
c
      integer :: ierr, status(MPI_STATUS_SIZE), blk, felem, num_enodes,
     &           num_enode_dof, totdof, span, utsz, blk_owner,
     &           inblk, i, count, count_to_receive, count_to_send,
     &           nrow_block, received_count, send_size_data(2), source,
     &           tag, request, num_requests, idx
      integer, save :: local_count = 0
      integer, allocatable :: requests(:)
      logical :: ok
      logical, parameter :: debug = .false., check_1 = .false.,
     &                      check_2 = .false., check_3 = .false.,
     &                      local_debug = .false., check_4 = .false.
      double precision :: start_move_time, end_move_time
c
      if( numprocs .eq. 1 ) return
c
c              alert workers to sync here
c
      call wmpi_alert_slaves ( 7 )
c
c              create blocks on root for element stiffness owned
c              by workers. optionally check blocks owned by this
c              rank for <nans>
c
c
      if( local_debug ) call wmpi_combine_stf_mess( 0 )
      if( root_processor ) call estiff_allocate( 2 )
      if( check_1 ) call estiff_allocate( 3 )
c
c              root drives process to bring blocks back in order
c              1, 2, 3, ... nelblk. Sends message to worker who owns
c              the block to transmit it. workers run a loop blocking on
c              a receive waiting to get a block number.
c
c              when done, root sends workers a 0 block number.
c
      start_move_time = MPI_WTIME()
c
      if( root_processor ) then
c
         num_requests = 0
         do blk = 1, nelblk ! count blocks not owned by root
           blk_owner      = elblks(2,blk)
           if( blk_owner .eq. 0 ) cycle ! skip blocks owned by root
           num_requests = num_requests + 1
         end do  ! over all blocks

         allocate ( requests( num_requests ) )

         idx = 1
         do blk = 1, nelblk ! bring blocks to root in order
           blk_owner      = elblks(2,blk)
           if( blk_owner .eq. 0 ) cycle ! skip blocks owned by root
           felem          = elblks(1,blk)
           num_enodes     = iprops(2,felem)
           num_enode_dof  = iprops(4,felem)
           totdof         = num_enodes * num_enode_dof
           span           = elblks(0,blk)
           utsz           = ((totdof*totdof)-totdof)/2 + totdof
           nrow_block     = utsz
           if( asymmetric_assembly ) nrow_block = totdof**2
           count_to_receive = nrow_block * span
           if( local_debug ) call wmpi_combine_stf_mess( 4 )
c           call MPI_SEND( blk, 1, MPI_INTEGER, blk_owner, 21,
c     &                    MPI_COMM_WORLD, ierr)
c           call wmpi_check_error( 1, myid, ierr, out )
c           if( check_4 ) then
c             call MPI_RECV( send_size_data, 2, MPI_INTEGER, blk_owner,
c     &                      13, MPI_COMM_WORLD, status, ierr )
c             call wmpi_check_error( 5, myid, ierr, out )
c             ok = ( send_size_data(1) .eq. nrow_block ) .and.
c     &          ( send_size_data(2) .eq. span )
c             if( .not. ok ) call wmpi_combine_stf_mess( 10 )
c           end if
c           call MPI_RECV( estiff_blocks(blk)%ptr(1,1),
c     &                    count_to_receive, MPI_VAL, blk_owner,
c     &                    14, MPI_COMM_WORLD, status, ierr )
           call MPI_IRECV( estiff_blocks(blk)%ptr(1,1),
     &                     count_to_receive, MPI_VAL, blk_owner,
     &                     14, MPI_COMM_WORLD, request, ierr )
           requests(idx) = request
           idx = idx + 1
           call wmpi_check_error( 2, myid, ierr, out )
c           call MPI_GET_COUNT( status, MPI_VAL, received_count, ierr )
c           if( received_count .ne. count_to_receive )
c     &          call wmpi_combine_stf_mess( 1 )
           if( check_3 ) call wmpi_combine_stf_check( 0, blk )
         end do  ! over all blocks

         call MPI_WAITALL( num_requests, requests, MPI_STATUSES_IGNORE,
     &                     ierr )
         deallocate ( requests )
c
         if( local_debug ) call wmpi_combine_stf_mess( 2 )
c         do i = 1, numprocs - 1
c           call MPI_SEND( 0, 1, MPI_INTEGER, i, 21, MPI_COMM_WORLD,
c     &                    ierr )
c         end do
         if( local_debug ) call wmpi_combine_stf_mess( 3 )
c
      else ! worker code
c
         !do
         do blk = 1, nelblk
c           call MPI_RECV( inblk, 1, MPI_INTEGER, 0, 21, MPI_COMM_WORLD,
c     &                    status, ierr )
c           call wmpi_check_error( 3, myid, ierr, out )
c           if( local_debug ) call wmpi_combine_stf_mess( 7 )
c           if( inblk .eq. 0 ) then
c             if( local_debug ) call wmpi_combine_stf_mess( 5 )
c             exit  ! done with loop
c           end if
           blk_owner      = elblks(2,blk)
           if( myid .ne. blk_owner )
     &         cycle
c     &        call wmpi_combine_stf_mess( 6 ) ! root sent wrong block #
           felem          = elblks(1,blk)
           num_enodes     = iprops(2,felem)
           num_enode_dof  = iprops(4,felem)
           totdof         = num_enodes * num_enode_dof
           span           = elblks(0,blk)
           utsz           = ((totdof*totdof)-totdof)/2 + totdof
           nrow_block     = utsz
           if( asymmetric_assembly ) nrow_block = totdof**2
           count_to_send = nrow_block * span
           if( check_3 ) call wmpi_combine_stf_check( 1, blk )
c           if( check_4 ) then
c              send_size_data(1) = nrow_block
c              send_size_data(2) = span
c              call MPI_SEND( send_size_data, 2, MPI_INTEGER, 0, 13,
c     &                       MPI_COMM_WORLD, ierr )
c           end if
c           call MPI_SEND( estiff_blocks(blk)%ptr(1,1),
c     &                       count_to_send, MPI_VAL, 0, 14,
c     &                       MPI_COMM_WORLD, ierr )
           call MPI_ISEND( estiff_blocks(blk)%ptr(1,1),
     &                     count_to_send, MPI_VAL, 0, 14,
     &                     MPI_COMM_WORLD, request, ierr )
c           call MPI_WAIT(request, status, ierr)
           call wmpi_check_error( 4, myid, ierr, out )
          end do ! worker wait loop
      end if
c
      end_move_time = MPI_WTIME()
c
      if( root_processor .and. check_2 ) then
          if( local_debug ) call wmpi_combine_stf_mess( 9 )
          call estiff_allocate( 4 ) ! check all blks root
      end if
c
      if( root_processor ) then ! .and. local_count .eq. 0 ) then
         write(out,9010) end_move_time - start_move_time
         local_count = 1
      end if
c
      if( local_debug ) call  wmpi_combine_stf_mess( 8 )
      return
c
 9010 format(10x,'>>> walltime to move elem stiffs:',f14.3 )
 9020 format(
     &   15x, 'FATAL ERRROR: wmpi_combine_stf. inconsistent',
     & /,15x, '              block size. blk: ',i8,
     & /,15x  '              Job aborted.')
c
      contains
c     ========

      subroutine wmpi_combine_stf_mess( messno )
      implicit none
c
      integer :: messno
c
      select case( messno )
c
        case( 0 )
          write(out,*) '... entered wmpi_combine_stf_mess on rank: ',
     &                   myid
        case( 1 )
          write(out,9010) blk, blk_owner, count_to_receive,
     &                    received_count
          call die_abort
        case( 2 )
          write(out,*) '.... sending termination 0 to each worker'
        case( 3 )
          write(out,*) '.... all workers told done'
        case( 4 )
          write(out,*) '.... root requests block, size : ', blk,
     &                count_to_receive
        case( 5 )
          write(out,*) '.... worker recd quit from root: ',myid
        case( 6 )
          write(out,9030) myid, inblk
          call die_abort
        case( 7 )
          write(out,*) '.... worker recd block number. myid, block: ',
     &                 myid, inblk
        case( 8 )
          write(out,*) '... leaving wmpi_combine_stf_mess on rank: ',
     &                   myid
        case( 9 )
          write(out,*) '... all blocks on root. check again for nans'
        case( 10 )
          write(out,9040) blk, blk_owner, nrow_block, span,
     &         send_size_data(1), send_size_data(2)
          call die_abort
        case default
          write(out,9020) myid, messno
c
      end select
c
      return
c
 9010 format( '>>>  FATAL ERROR: wmpi_combine_stf',
     & /,     '     block recd from worker has wrong size',
     & /,10x,'blk, blk_owner, count_to_receive, received_count: ',
     & 4i6,
     & /,     '     job terminated')
 9020 format( '>>>  FATAL ERROR: wmpi_combine_stf',
     & /,     '     invalid message number. myid, error:',2i5,
     & /,     '     job terminated')
 9030 format( '>>>  FATAL ERROR: wmpi_combine_stf',
     & /,     '     root sent wrong block number to worker',
     & /,10x,'myid, inblk: ', 2i6,
     & /,     '     job terminated')
 9040 format( '>>>  FATAL ERROR: wmpi_combine_stf',
     & /,     '     mismatched sizes on send to root',
     & /,10x,'blk, blk_owner, nrow_block, span, send_size_data(1,2): ',
     & 6i6,
     & /,     '     job terminated')
c
      end subroutine wmpi_combine_stf_mess

      subroutine wmpi_combine_stf_check( type, blk_to_chk )
      implicit none
c
      integer :: type, blk_to_chk, i, j, mcount
      logical :: header
c
      header = .true.
      mcount = 0
c
      do i = 1, span
        do j = 1, nrow_block ! look for a nan value
         if( isnan( estiff_blocks(blk_to_chk)%ptr(1,1) ) ) then
             if( header ) then
               write(out,9000) type
               header = .false.
             end if
             write(out,9010) myid, blk_to_chk, span, nrow_block,
     &                       felem+i-1, j
             mcount = mcount + 1
             if( mcount > 10 ) call die_abort
         end if
        end do
      end do
c
      return

 9000 format( '>> Internal Errors: wmpi_combine_stf_check. type: ',i3)
 9010 format(10x,' myid, blk, span, nterms, element, stiff term:',6i7)
c
      end subroutine wmpi_combine_stf_check
c
      end subroutine wmpi_combine_stf

c     ****************************************************************
c     *                                                              *
c     *              subroutine wmpi_combine_stf_parallel            *
c     *                                                              *
c     *                   last modified : 12/20/2017 rhd             *
c     *                                                              *
c     *     gathers the block element stiffnesses from workers.      *
c     *     when done, root has copy of all element stiff blocks.    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_combine_stf_parallel
      use global_data ! old common.main
      use elem_block_data, only: estiff_blocks
      use main_data, only: asymmetric_assembly
      implicit none
      include "mpif.h"
c
      integer :: ierr, status(MPI_STATUS_SIZE), blk, felem, num_enodes,
     &           num_enode_dof, totdof, span, utsz, blk_owner,
     &           inblk, i, count, count_to_receive, count_to_send,
     &           nrow_block, received_count, send_size_data(2),
     &           rank
      integer, save :: local_count = 0
      logical :: ok
      logical, parameter :: debug = .false., check_1 = .false.,
     &                      check_2 = .false., check_3 = .false.,
     &                      local_debug = .false., check_4 = .false.
      double precision :: start_move_time, end_move_time
c
      if( numprocs .eq. 1 ) return
c
c              alert workers to sync here
c
      call wmpi_alert_slaves ( 7 )
c
c              create blocks on root for element stiffness owned
c              by workers. optionally check blocks owned by this
c              rank for <nans>
c
c
      if( local_debug ) call wmpi_combine_stf_mess( 0 )
      if( root_processor ) call estiff_allocate( 2 )
      if( check_1 ) call estiff_allocate( 3 )
c
c              root drives process to bring blocks back in order
c              1, 2, 3, ... nelblk. Sends message to worker who owns
c              the block to transmit it. workers run a loop blocking on
c              a receive waiting to get a block number.
c
c              when done, root sends workers a 0 block number.
c
      start_move_time = MPI_WTIME()
c
      if( root_processor ) then
c
c
c$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)
c$OMP&   PRIVATE( rank, blk, blk_owner, felem, num_enodes,
c$OMP&   num_enode_dof, totdof, span, utsz, nrow_block,
c$OMP&   count_to_receive, ierr, myid, send_size_data,
c$OMP&   status, received_count, out)
       do rank = 1, numprocs - 1
         do blk = 1, nelblk ! bring blocks to root in order
           blk_owner      = elblks(2,blk)
           if( blk_owner .ne. rank ) cycle ! skip blocks owned by root
           felem          = elblks(1,blk)
           num_enodes     = iprops(2,felem)
           num_enode_dof  = iprops(4,felem)
           totdof         = num_enodes * num_enode_dof
           span           = elblks(0,blk)
           utsz           = ((totdof*totdof)-totdof)/2 + totdof
           nrow_block     = utsz
           if( asymmetric_assembly ) nrow_block = totdof**2
           count_to_receive = nrow_block * span
           if( local_debug ) call wmpi_combine_stf_mess( 4 )
           call MPI_SEND( blk, 1, MPI_INTEGER, blk_owner, 21,
     &                    MPI_COMM_WORLD, ierr)
           call wmpi_check_error( 1, 0, ierr, out )
           if( check_4 ) then
             call MPI_RECV( send_size_data, 2, MPI_INTEGER, blk_owner,
     &                      13, MPI_COMM_WORLD, status, ierr )
             call wmpi_check_error( 5, myid, ierr, out )
             ok = ( send_size_data(1) .eq. nrow_block ) .and.
     &          ( send_size_data(2) .eq. span )
             if( .not. ok ) call wmpi_combine_stf_mess( 10 )
           end if
           call MPI_RECV( estiff_blocks(blk)%ptr(1,1),
     &                    count_to_receive, MPI_VAL, blk_owner,
     &                    14, MPI_COMM_WORLD, status, ierr )
           call wmpi_check_error( 2, 0, ierr, out )
           call MPI_GET_COUNT( status, MPI_VAL, received_count, ierr )
           if( received_count .ne. count_to_receive )
     &          call wmpi_combine_stf_mess( 1 )
           if( check_3 ) call wmpi_combine_stf_check( 0, blk )
         end do  ! over all blocks
      end do ! over ranks
c$OMP END PARALLEL DO

         if( local_debug ) call wmpi_combine_stf_mess( 2 )
         do i = 1, numprocs - 1
           call MPI_SEND( 0, 1, MPI_INTEGER, i, 21, MPI_COMM_WORLD,
     &                    ierr )
         end do
         if( local_debug ) call wmpi_combine_stf_mess( 3 )
c
      else ! worker code
c
         do
           call MPI_RECV( inblk, 1, MPI_INTEGER, 0, 21, MPI_COMM_WORLD,
     &                    status, ierr )
           call wmpi_check_error( 3, myid, ierr, out )
           if( local_debug ) call wmpi_combine_stf_mess( 7 )
           if( inblk .eq. 0 ) then
             if( local_debug ) call wmpi_combine_stf_mess( 5 )
             exit  ! done with loop
           end if
           if( myid .ne. elblks(2,inblk) )
     &        call wmpi_combine_stf_mess( 6 ) ! root sent wrong block #
           felem          = elblks(1,inblk)
           num_enodes     = iprops(2,felem)
           num_enode_dof  = iprops(4,felem)
           totdof         = num_enodes * num_enode_dof
           span           = elblks(0,inblk)
           utsz           = ((totdof*totdof)-totdof)/2 + totdof
           nrow_block     = utsz
           if( asymmetric_assembly ) nrow_block = totdof**2
           count_to_send = nrow_block * span
           if( check_3 ) call wmpi_combine_stf_check( 1, inblk )
           if( check_4 ) then
              send_size_data(1) = nrow_block
              send_size_data(2) = span
              call MPI_SEND( send_size_data, 2, MPI_INTEGER, 0, 13,
     &                       MPI_COMM_WORLD, ierr )
           end if
           call MPI_SEND( estiff_blocks(inblk)%ptr(1,1),
     &                       count_to_send, MPI_VAL, 0, 14,
     &                       MPI_COMM_WORLD, ierr )
           call wmpi_check_error( 4, myid, ierr, out )
          end do ! worker wait loop
      end if
c
      call wmpi_wait
      end_move_time = MPI_WTIME()
c
      if( root_processor .and. check_2 ) then
          if( local_debug ) call wmpi_combine_stf_mess( 9 )
          call estiff_allocate( 4 ) ! check all blks root
      end if
c
      if( root_processor .and. local_count .eq. 0 ) then
         write(out,9010) end_move_time - start_move_time
         local_count = 1
      end if
c
      if( local_debug ) call  wmpi_combine_stf_mess( 8 )
      return
c
 9010 format(10x,'>>> walltime to move elem stiffs:',f14.3 )
 9020 format(
     &   15x, 'FATAL ERRROR: wmpi_combine_stf. inconsistent',
     & /,15x, '              block size. blk: ',i8,
     & /,15x  '              Job aborted.')
c
      contains
c     ========

      subroutine wmpi_combine_stf_mess( messno )
      implicit none
c
      integer :: messno
c
      select case( messno )
c
        case( 0 )
          write(out,*) '... entered wmpi_combine_stf_mess on rank: ',
     &                   myid
        case( 1 )
          write(out,9010) blk, blk_owner, count_to_receive,
     &                    received_count
          call die_abort
        case( 2 )
          write(out,*) '.... sending termination 0 to each worker'
        case( 3 )
          write(out,*) '.... all workers told done'
        case( 4 )
          write(out,*) '.... root requests block, size : ', blk,
     &                count_to_receive
        case( 5 )
          write(out,*) '.... worker recd quit from root: ',myid
        case( 6 )
          write(out,9030) myid, inblk
          call die_abort
        case( 7 )
          write(out,*) '.... worker recd block number. myid, block: ',
     &                 myid, inblk
        case( 8 )
          write(out,*) '... leaving wmpi_combine_stf_mess on rank: ',
     &                   myid
        case( 9 )
          write(out,*) '... all blocks on root. check again for nans'
        case( 10 )
          write(out,9040) blk, blk_owner, nrow_block, span,
     &         send_size_data(1), send_size_data(2)
          call die_abort
        case default
          write(out,9020) myid, messno
c
      end select
c
      return
c
 9010 format( '>>>  FATAL ERROR: wmpi_combine_stf',
     & /,     '     block recd from worker has wrong size',
     & /,10x,'blk, blk_owner, count_to_receive, received_count: ',
     & 4i6,
     & /,     '     job terminated')
 9020 format( '>>>  FATAL ERROR: wmpi_combine_stf',
     & /,     '     invalid message number. myid, error:',2i5,
     & /,     '     job terminated')
 9030 format( '>>>  FATAL ERROR: wmpi_combine_stf',
     & /,     '     root sent wrong block number to worker',
     & /,10x,'myid, inblk: ', 2i6,
     & /,     '     job terminated')
 9040 format( '>>>  FATAL ERROR: wmpi_combine_stf',
     & /,     '     mismatched sizes on send to root',
     & /,10x,'blk, blk_owner, nrow_block, span, send_size_data(1,2): ',
     & 6i6,
     & /,     '     job terminated')
c
      end subroutine wmpi_combine_stf_mess

      subroutine wmpi_combine_stf_check( type, blk_to_chk )
      implicit none
c
      integer :: type, blk_to_chk, i, j, mcount
      logical :: header
c
      header = .true.
      mcount = 0
c
      do i = 1, span
        do j = 1, nrow_block ! look for a nan value
         if( isnan( estiff_blocks(blk_to_chk)%ptr(1,1) ) ) then
             if( header ) then
               write(out,9000) type
               header = .false.
             end if
             write(out,9010) myid, blk_to_chk, span, nrow_block,
     &                       felem+i-1, j
             mcount = mcount + 1
             if( mcount > 10 ) call die_abort
         end if
        end do
      end do
c
      return

 9000 format( '>> Internal Errors: wmpi_combine_stf_check. type: ',i3)
 9010 format(10x,' myid, blk, span, nterms, element, stiff term:',6i7)
c
      end subroutine wmpi_combine_stf_check
c
      end subroutine wmpi_combine_stf_parallel
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_combine_stf             *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 6/23/2017 rhd              *
c     *                                                              *
c     *     gathers the block element stiffnesses from workers.      *
c     *     when done, root has copy of all element stiff blocks.    *                     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_combine_stf_original
      use global_data ! old common.main
      use elem_block_data, only: estiff_blocks
      use main_data, only: asymmetric_assembly
      implicit none
      include "mpif.h"
c
      integer :: ierr, status(MPI_STATUS_SIZE), blk, felem, num_enodes,
     &           num_enode_dof, totdof, span, utsz, blk_owner,
     &           inblk, i, count, count_to_receive, count_to_send,
     &           nrow_block, received_count, send_size_data(2)
      integer, save :: local_count = 0
      logical :: ok
      logical, parameter :: debug = .false., check_1 = .false.,
     &                      check_2 = .false., check_3 = .false.,
     &                      local_debug = .false., check_4 = .false.
      double precision :: start_move_time, end_move_time
c
      if( numprocs .eq. 1 ) return
c
c              alert workers to sync here
c
      call wmpi_alert_slaves ( 7 )
c
c              create blocks on root for element stiffness owned
c              by workers. optionally check blocks owned by this
c              rank for <nans>
c
c
      if( local_debug ) call wmpi_combine_stf_mess( 0 )
      if( root_processor ) call estiff_allocate( 2 )
      if( check_1 ) call estiff_allocate( 3 )
c
c              root drives process to bring blocks back in order
c              1, 2, 3, ... nelblk. Sends message to worker who owns
c              the block to transmit it. workers run a loop blocking on
c              a receive waiting to get a block number.
c
c              when done, root sends workers a 0 block number.
c
      start_move_time = MPI_WTIME()
c
      if( root_processor ) then
c
         do blk = 1, nelblk ! bring blocks to root in order
           blk_owner      = elblks(2,blk)
           if( blk_owner .eq. 0 ) cycle ! skip blocks owned by root
           felem          = elblks(1,blk)
           num_enodes     = iprops(2,felem)
           num_enode_dof  = iprops(4,felem)
           totdof         = num_enodes * num_enode_dof
           span           = elblks(0,blk)
           utsz           = ((totdof*totdof)-totdof)/2 + totdof
           nrow_block     = utsz
           if( asymmetric_assembly ) nrow_block = totdof**2
           count_to_receive = nrow_block * span
           if( local_debug ) call wmpi_combine_stf_mess( 4 )
           call MPI_SEND( blk, 1, MPI_INTEGER, blk_owner, 21,
     &                    MPI_COMM_WORLD, ierr)
           call wmpi_check_error( 1, myid, ierr, out )
           if( check_4 ) then
             call MPI_RECV( send_size_data, 2, MPI_INTEGER, blk_owner,
     &                      13, MPI_COMM_WORLD, status, ierr )
             call wmpi_check_error( 5, myid, ierr, out )
             ok = ( send_size_data(1) .eq. nrow_block ) .and.
     &          ( send_size_data(2) .eq. span )
             if( .not. ok ) call wmpi_combine_stf_mess( 10 )
           end if
           call MPI_RECV( estiff_blocks(blk)%ptr(1,1),
     &                    count_to_receive, MPI_VAL, blk_owner,
     &                    14, MPI_COMM_WORLD, status, ierr )
           call wmpi_check_error( 2, myid, ierr, out )
           call MPI_GET_COUNT( status, MPI_VAL, received_count, ierr )
           if( received_count .ne. count_to_receive )
     &          call wmpi_combine_stf_mess( 1 )
           if( check_3 ) call wmpi_combine_stf_check( 0, blk )
         end do  ! over all blocks
c
         if( local_debug ) call wmpi_combine_stf_mess( 2 )
         do i = 1, numprocs - 1
           call MPI_SEND( 0, 1, MPI_INTEGER, i, 21, MPI_COMM_WORLD,
     &                    ierr )
         end do
         if( local_debug ) call wmpi_combine_stf_mess( 3 )
c
      else ! worker code
c
         do
           call MPI_RECV( inblk, 1, MPI_INTEGER, 0, 21, MPI_COMM_WORLD,
     &                    status, ierr )
           call wmpi_check_error( 3, myid, ierr, out )
           if( local_debug ) call wmpi_combine_stf_mess( 7 )
           if( inblk .eq. 0 ) then
             if( local_debug ) call wmpi_combine_stf_mess( 5 )
             exit  ! done with loop
           end if
           if( myid .ne. elblks(2,inblk) )
     &        call wmpi_combine_stf_mess( 6 ) ! root sent wrong block #
           felem          = elblks(1,inblk)
           num_enodes     = iprops(2,felem)
           num_enode_dof  = iprops(4,felem)
           totdof         = num_enodes * num_enode_dof
           span           = elblks(0,inblk)
           utsz           = ((totdof*totdof)-totdof)/2 + totdof
           nrow_block     = utsz
           if( asymmetric_assembly ) nrow_block = totdof**2
           count_to_send = nrow_block * span
           if( check_3 ) call wmpi_combine_stf_check( 1, inblk )
           if( check_4 ) then
              send_size_data(1) = nrow_block
              send_size_data(2) = span
              call MPI_SEND( send_size_data, 2, MPI_INTEGER, 0, 13,
     &                       MPI_COMM_WORLD, ierr )
           end if
           call MPI_SEND( estiff_blocks(inblk)%ptr(1,1),
     &                       count_to_send, MPI_VAL, 0, 14,
     &                       MPI_COMM_WORLD, ierr )
           call wmpi_check_error( 4, myid, ierr, out )
          end do ! worker wait loop
      end if
c
      end_move_time = MPI_WTIME()
c
      if( root_processor .and. check_2 ) then
          if( local_debug ) call wmpi_combine_stf_mess( 9 )
          call estiff_allocate( 4 ) ! check all blks root
      end if
c
      if( root_processor .and. local_count .eq. 0 ) then
         write(out,9010) end_move_time - start_move_time
         local_count = 1
      end if
c
      if( local_debug ) call  wmpi_combine_stf_mess( 8 )
      return
c
 9010 format(10x,'>>> walltime to move elem stiffs:',f14.3 )
 9020 format(
     &   15x, 'FATAL ERRROR: wmpi_combine_stf. inconsistent',
     & /,15x, '              block size. blk: ',i8,
     & /,15x  '              Job aborted.')
c
      contains
c     ========

      subroutine wmpi_combine_stf_mess( messno )
      implicit none
c
      integer :: messno
c
      select case( messno )
c
        case( 0 )
          write(out,*) '... entered wmpi_combine_stf_mess on rank: ',
     &                   myid
        case( 1 )
          write(out,9010) blk, blk_owner, count_to_receive,
     &                    received_count
          call die_abort
        case( 2 )
          write(out,*) '.... sending termination 0 to each worker'
        case( 3 )
          write(out,*) '.... all workers told done'
        case( 4 )
          write(out,*) '.... root requests block, size : ', blk,
     &                count_to_receive
        case( 5 )
          write(out,*) '.... worker recd quit from root: ',myid
        case( 6 )
          write(out,9030) myid, inblk
          call die_abort
        case( 7 )
          write(out,*) '.... worker recd block number. myid, block: ',
     &                 myid, inblk
        case( 8 )
          write(out,*) '... leaving wmpi_combine_stf_mess on rank: ',
     &                   myid
        case( 9 )
          write(out,*) '... all blocks on root. check again for nans'
        case( 10 )
          write(out,9040) blk, blk_owner, nrow_block, span,
     &         send_size_data(1), send_size_data(2)
          call die_abort
        case default
          write(out,9020) myid, messno
c
      end select
c
      return
c
 9010 format( '>>>  FATAL ERROR: wmpi_combine_stf',
     & /,     '     block recd from worker has wrong size',
     & /,10x,'blk, blk_owner, count_to_receive, received_count: ',
     & 4i6,
     & /,     '     job terminated')
 9020 format( '>>>  FATAL ERROR: wmpi_combine_stf',
     & /,     '     invalid message number. myid, error:',2i5,
     & /,     '     job terminated')
 9030 format( '>>>  FATAL ERROR: wmpi_combine_stf',
     & /,     '     root sent wrong block number to worker',
     & /,10x,'myid, inblk: ', 2i6,
     & /,     '     job terminated')
 9040 format( '>>>  FATAL ERROR: wmpi_combine_stf',
     & /,     '     mismatched sizes on send to root',
     & /,10x,'blk, blk_owner, nrow_block, span, send_size_data(1,2): ',
     & 6i6,
     & /,     '     job terminated')
c
      end subroutine wmpi_combine_stf_mess

      subroutine wmpi_combine_stf_check( type, blk_to_chk )
      implicit none
c
      integer :: type, blk_to_chk, i, j, mcount
      logical :: header
c
      header = .true.
      mcount = 0
c
      do i = 1, span
        do j = 1, nrow_block ! look for a nan value
         if( isnan( estiff_blocks(blk_to_chk)%ptr(1,1) ) ) then
             if( header ) then
               write(out,9000) type
               header = .false.
             end if
             write(out,9010) myid, blk_to_chk, span, nrow_block,
     &                       felem+i-1, j
             mcount = mcount + 1
             if( mcount > 10 ) call die_abort
         end if
        end do
      end do
c
      return

 9000 format( '>> Internal Errors: wmpi_combine_stf_check. type: ',i3)
 9010 format(10x,' myid, blk, span, nterms, element, stiff term:',6i7)
c
      end subroutine wmpi_combine_stf_check
c
      end subroutine wmpi_combine_stf_original
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_get_str                 *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 6/28/2018 rhd              *
c     *                                                              *
c     *     this subroutine gathers element data (stresses,strains,  *
c     *     element volumes)                                         *
c     *     distributed among processors back to the root processor. *
c     *     we keep track if we have gathered the data since         *
c     *     the last time step so that we won't do all the           *
c     *     communication twice.                                     *
c     *                                                              *
c     *     do specifies what action to take in relation to the      *
c     *     blocked data structures:                                 *
c     *        1 = tell root its element data for non-owned elems    *
c     *            is wrong                                          *
c     *        2 = gather element stresses, material history and     *
c     *            [D]s on root                                      *
c     *        3 = <available >                                      *
c     *        4 = gather element strains on root                    *
c     *        5 = gather element volumes on root                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_get_str( do )
      use global_data ! old common.main
      use elem_block_data
      implicit none
      include "mpif.h"
c
      integer :: do
c
c                    local declarations
c
      integer :: ierr, status(MPI_STATUS_SIZE), blk, felem, num_enodes,
     &           num_enode_dof, totdof, span, utsz, mat_type, histsize,
     &           cepsize, ngp, dowhat
      character(len=80) :: line
      logical, parameter :: local_debug = .false.
      logical, save :: root_has_stress, root_has_strain
      double precision, parameter :: zero = 0.0d0
c
      if( numprocs .eq. 1 ) return
c
c                    root does some checks. may be no work to do
c
      if( root_processor ) then
         call wmpi_get_str_a( dowhat )
         if( dowhat == 1 ) return
      end if
c
c                   we need to get all strains,
c                   stresses, volume data from workers. let workers
c                   the suncronization here works via root and
c                   workers processing all blocks. only 1 worker owns a
c                   block so the sync to-from root to the worker is
c                   simple but likely not very effificent.
c
      call wmpi_alert_slaves ( 20 )
      call wmpi_bcast_int ( do )
c
c                   get all of the element block info back to
c                   the root processor. This code synchronizes
c                   the transfers. root and all workers are executing
c                   this code.
c                   if root owns block, skip.
c                   if worker does not own block, skip.
c                   after this logic, block "blk" is on the worker
c                   with data for sending to root and root posts
c                   receives to that rank.
c
      do blk = 1, nelblk
c
        if( worker_processor ) then
            if( elblks(2,blk) .ne. myid ) cycle
         end if
         if( root_processor ) then
            if( elblks(2,blk) .eq. myid ) cycle
         end if
c
         felem          = elblks(1,blk)
         num_enodes     = iprops(2,felem)
         num_enode_dof  = iprops(4,felem)
         totdof         = num_enodes * num_enode_dof
         span           = elblks(0,blk)
         utsz           = ((totdof*totdof)-totdof)/2 + totdof
         ngp            = iprops(6,felem)
         mat_type       = iprops(25,felem)
         histsize       = history_blk_list(blk)
         cepsize        = cep_blk_list(blk)
c
c                         if do = 2, get history, cep, unrotated
c                         stresses, and rotation data
c
         select case( do )
         case( 2 )
            call wmpi_get_str_2 ! history, cep [D], unrotated sig, rot mats
c         case( 3 )    !! available
c            continue
         case( 4 )
           call wmpi_get_str_4  ! strains
         case( 5 )
           call wmpi_get_str_5  ! element volumnes
         case default
           write(out,9000)
           call die_abort
         end select
c
      end do
c
      return
c
 9000 format ('>>>  FATAL ERROR:   wmpi_get_str',/,
     &        '                    do: ',i10,
     &    /,  '                    job terminated')
c

      contains
c     ========
c
      subroutine wmpi_get_str_a( action )
      implicit none
c
      integer :: action
c
c               action = 1    caller exectues a return
c
      select case( do )
      case( 1 )
c
c               mark stresses and strains that the root processor has
c               for blocks it does not own as out of date. Thus next
c               time we return to this routine looking
c               to get the full stresses, we will have to do communication.
c
         root_has_stress = .false.
         root_has_strain = .false.
         action = 1
c
      case( 2 )
c
c               getting the major stress data from workers to root.
c               This data is needed
c               by several different algorithms, including stress
c               output, and restart file construction. if root already
c               has the current stresses for this step, return.
c               otherwise go get the stresses.
c
        if( root_has_stress ) then
          action = 1
        else
          action = 0
          root_has_stress = .true.
        end if
c
c               available
c
      case( 3 )
        action = 1
c
      case( 4 )
c
c               getting the strain data. This data is needed by
c               several different algorithms, including strain output,
c               crack growth calculation, and database construction.
c               if root already has the current strains for this
c               step, otherwise go get the strains.
c
        if( root_has_strain ) then
          action = 1
        else
          action = 0
          root_has_strain = .true.
        end if
c
      case( 5 )
c
c               element volumnes. alwasy get them
c
        action = 0
c
      case default
        write(out,9000)
        call die_abort
c
      end select
      return
c
 9000 format ('>>>  FATAL ERROR:   wmpi_get_str_a',/,
     &        '                    do: ',i10,
     &    /,  '                    job terminated')
c
      end subroutine wmpi_get_str_a
c
      subroutine wmpi_get_str_2
      implicit none
c
      integer :: block_size
c
c              element [D] tangents
c
      block_size = span * ngp * cepsize
      if( root_processor ) then
         cep_blocks(blk)%vector(1:block_size) = zero
         call MPI_RECV( cep_blocks(blk)%vector(1), block_size,
     &        MPI_VAL, elblks(2,blk), 14, MPI_COMM_WORLD, status,
     &        ierr )
      else
         call MPI_SEND( cep_blocks(blk)%vector(1), block_size,
     &        MPI_VAL, 0, 14, MPI_COMM_WORLD, ierr )
      end if
c
c              element histories
c
      block_size = span * ngp * histsize
      if( root_processor ) then
         history_blocks(blk)%ptr(1:block_size) = zero
         call MPI_RECV( history_blocks(blk)%ptr(1), block_size,
     &        MPI_VAL, elblks(2,blk), 14, MPI_COMM_WORLD, status,
     &        ierr )
      else
         call MPI_SEND( history_blocks(blk)%ptr(1), block_size,
     &        MPI_VAL, 0, 14, MPI_COMM_WORLD, ierr )
      end if
C
c              element stresses
c
      block_size    = span * ngp * nstrs
      if( root_processor ) then
         urcs_n_blocks(blk)%ptr(1:block_size) = zero
         call MPI_RECV( urcs_n_blocks(blk)%ptr(1), block_size,
     &        MPI_VAL, elblks(2,blk), 14, MPI_COMM_WORLD, status,
     &        ierr )
      else
         call MPI_SEND( urcs_n_blocks(blk)%ptr(1), block_size,
     &        MPI_VAL, 0, 14, MPI_COMM_WORLD, ierr )
      end if
c
c              element rotations
c
      if( rot_blk_list(blk) .eq. 0) return
      block_size = span * ngp * 9
      if( root_processor ) then
         if( allocated(rot_n1_blocks) ) then
            rot_n1_blocks(blk)%ptr(1:block_size) = zero
            call MPI_RECV( rot_n1_blocks(blk)%ptr(1), block_size,
     &           MPI_VAL, elblks(2,blk), 15, MPI_COMM_WORLD,
     &           status, ierr )
         end if
      else
         if( allocated(rot_n1_blocks) ) then
            call MPI_SEND( rot_n1_blocks(blk)%ptr(1), block_size,
     &           MPI_VAL, 0, 15, MPI_COMM_WORLD, ierr )
         end if
      end if
c
      return

      end subroutine wmpi_get_str_2
c
      subroutine wmpi_get_str_4
      implicit none
c
      integer :: block_size

      block_size = span * ngp * nstr
      if( root_processor ) then
        eps_n_blocks(blk)%ptr(1:block_size) = zero
        call MPI_RECV( eps_n_blocks(blk)%ptr(1), block_size,
     &                 MPI_VAL, elblks(2,blk), 14, MPI_COMM_WORLD,
     &                 status, ierr )
      else
        call MPI_SEND( eps_n_blocks(blk)%ptr(1), block_size,
     &                 MPI_VAL, 0, 14, MPI_COMM_WORLD, ierr )
      end if
c
      return
c
      end subroutine wmpi_get_str_4
c
      subroutine wmpi_get_str_5
      implicit none
c
      integer :: block_size
c
      block_size = span
      if( root_processor ) then
        element_vol_blocks(blk)%ptr(1:block_size) = zero
        call MPI_RECV( element_vol_blocks(blk)%ptr(1),
     &                 block_size, MPI_VAL, elblks(2,blk), 14,
     &                 MPI_COMM_WORLD, status, ierr )
      else
        call MPI_SEND( element_vol_blocks(blk)%ptr(1),
     &                 block_size, MPI_VAL, 0, 14, MPI_COMM_WORLD,
     &                 ierr )
      end if
c
      return
c
      end subroutine wmpi_get_str_5
      end subroutine wmpi_get_str

c     ****************************************************************
c     *                                                              *
c     *                subroutine wmpi_get_initial_state             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/16/2018 rhd              *
c     *                                                              *
c     *     makes all intiial state arrays present on rank 0         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_get_initial_state
      use global_data ! old common.main
      use elem_block_data, only: estiff_blocks, initial_state_data
      implicit none
      include "mpif.h"
c
      integer :: ierr, status(MPI_STATUS_SIZE), blk, felem,
     &           span, blk_owner, ngp, inblk, i, count_to_receive,
     &           count_to_send, received_count
     &
      logical, parameter :: local_debug = .false.
c
      if( numprocs .eq. 1 ) return
c
c              alert workers to sync here
c
      call wmpi_alert_slaves ( 33 )
c
      if( local_debug ) call wmpi_get_initial_state_mess( 0 )
c
c              root drives process to bring blocks back in order
c              1, 2, 3, ... nelblk. Sends message to worker who owns
c              the block to transmit it. workers run a loop blocking on
c              a receive waiting to get a block number.
c              allocate space on root for worker owned blocks
c
c              when done, root sends workers a 0 block number.
c
      if( root_processor ) then
c
         associate( x => initial_state_data )
c
         do blk = 1, nelblk ! bring blocks to root in order
           blk_owner      = elblks(2,blk)
           if( blk_owner .eq. 0 ) cycle ! skip blocks owned by root
           felem          = elblks(1,blk)
           span           = elblks(0,blk)
           ngp            = iprops(6,felem)
           count_to_receive = ngp * span
           if( .not. allocated( x(blk)%W_plastic_nis_block ) )
     &          allocate( x(blk)%W_plastic_nis_block(span,ngp) )
           call MPI_SEND( blk, 1, MPI_INTEGER, blk_owner, 21,
     &                    MPI_COMM_WORLD, ierr)
           call wmpi_check_error( 1, myid, ierr, out )
           call MPI_RECV( x(blk)%W_plastic_nis_block(1,1),
     &                    count_to_receive, MPI_VAL, blk_owner,
     &                    14, MPI_COMM_WORLD, status, ierr )
           call wmpi_check_error( 2, myid, ierr, out )
           call MPI_GET_COUNT( status, MPI_VAL, received_count, ierr )
           if( received_count .ne. count_to_receive )
     &          call wmpi_get_initial_state_mess( 1 )
         end do  ! over all blocks
c
         end associate
c
         do i = 1, numprocs - 1   !   workers we're done
           call MPI_SEND( 0, 1, MPI_INTEGER, i, 21, MPI_COMM_WORLD,
     &                    ierr )
         end do
c
      else ! worker code
c
         associate( y => initial_state_data )
c
         do
           call MPI_RECV( inblk, 1, MPI_INTEGER, 0, 21, MPI_COMM_WORLD,
     &                    status, ierr )
           call wmpi_check_error( 5, myid, ierr, out )
           if( inblk .eq. 0 ) exit  ! done with loop
           if( myid .ne. elblks(2,inblk) )
     &        call  wmpi_get_initial_state_mess( 6 ) ! wrong blk # sent
           felem          = elblks(1,inblk)
           span           = elblks(0,inblk)
           ngp            = iprops(6,felem)
           count_to_send  = ngp * span
           call MPI_SEND( y(inblk)%W_plastic_nis_block(1,1),
     &                       count_to_send, MPI_VAL, 0, 14,
     &                       MPI_COMM_WORLD, ierr )
           call wmpi_check_error( 6, myid, ierr, out )
          end do ! worker wait loop
c
          end associate
      end if
c
      return
c
      contains
c     ========
c
      subroutine wmpi_get_initial_state_mess( messno )
      implicit none
c
      integer :: messno
c
      select case( messno )
c
        case( 0 )
          write(out,*)
     &     '... entered wmpi_get_initial_state_mess on rank: ',
     &                   myid
        case( 1 )
          write(out,9010) blk, blk_owner, count_to_receive,
     &                    received_count
          call die_abort
        case( 3 )
          write(out,*) '.... send worker blk number'
        case( 4 )
        case( 5 )
          write(out,*) '.... worker recd quit from root: ',myid
        case( 6 )
          write(out,9030) myid, inblk
          call die_abort
        case( 7 )
          write(out,*) '.... worker recd block number. myid, block: ',
     &                 myid, inblk
        case( 8 )
          write(out,*) '... leaving  wmpi_get_initial_state_mess: ',
     &                   myid
        case( 10 )
      end select
c
      return
c
 9010 format( '>>>  FATAL ERROR: wmpi_get_initial_state_mess',
     & /,     '     block recd from worker has wrong size',
     & /,10x,'blk, blk_owner, count_to_receive, received_count: ',
     & 4i6,
     & /,     '     job terminated')
 9020 format( '>>>  FATAL ERROR: wmpi_get_initial_state_mess',
     & /,     '     invalid message number. myid, error:',2i5,
     & /,     '     job terminated')
 9030 format( '>>>  FATAL ERROR: wmpi_get_initial_state_mess',
     & /,     '     root sent wrong block number to worker',
     & /,10x,'myid, inblk: ', 2i6,
     & /,     '     job terminated')
c
      end subroutine wmpi_get_initial_state_mess
      end subroutine wmpi_get_initial_state


c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine wmpi_send_reopen                  *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 9/16/2018 rhd              *
c     *                                                              *
c     *     this subroutine sends information that the               *
c     *     workers need after a restart.                            *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_reopen
      use global_data ! old common.main
      use elem_block_data
      use main_data, only: eq_node_force_len, eq_node_force_indexes,
     &                     eq_node_forces, initial_state_option
      implicit none
      include "mpif.h"
c
c                    local declarations
c
      integer :: ierr, status(MPI_STATUS_SIZE), felem, num_enodes,
     &           num_enode_dof, totdof, span, utsz, ngp, mat_type,
     &           block_size, blk, idummy, wrkr_rank
      logical, save :: local_debug, root_has_stress, root_has_strain
      character(len=80) :: line
c
      if( numprocs .eq. 1 ) return
c
      call wmpi_alert_slaves ( 23 )
c
c                   displacements, velocities, accelerations.
c
      call MPI_BCAST( u, nodof, MPI_VAL, 0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( v, nodof, MPI_VAL, 0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( a, nodof, MPI_VAL, 0, MPI_COMM_WORLD, ierr )
c
c                   blocked data structure. synching
c                   again thru blocks - a block is owned by only 1
c                   worker.
c
      root_has_stress  = .false.
      root_has_strain = .false.
c
      do blk = 1, nelblk
c
         if( worker_processor ) then
            if( elblks(2,blk) .ne. myid ) cycle
         end if
         if( root_processor ) then
            if( elblks(2,blk) .eq. myid ) cycle
         end if
         if( root_processor ) wrkr_rank = elblks(2,blk)

c
         felem          = elblks(1,blk)
         num_enodes     = iprops(2,felem)
         num_enode_dof  = iprops(4,felem)
         totdof         = num_enodes * num_enode_dof
         span           = elblks(0,blk)
         utsz           = ((totdof*totdof)-totdof)/2 + totdof
         ngp            = iprops(6,felem)
         mat_type       = iprops(25,felem)
c
c                               element histories
c
         block_size = span * ngp * history_blk_list(blk)
         if( worker_processor ) then
            call MPI_RECV(history_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,0,13,MPI_COMM_WORLD,status,ierr)
            call vec_ops( history1_blocks(blk)%ptr(1),
     &           history_blocks(blk)%ptr(1), idummy, block_size, 5 )
         else
            call MPI_SEND(history_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,wrkr_rank,13,MPI_COMM_WORLD,ierr)
         end if
c
c                               element [D]s
c
         block_size = span * ngp * cep_blk_list(blk)
         if( worker_processor ) then
            call MPI_RECV(cep_blocks(blk)%vector(1),block_size,
     &           MPI_VAL,0,14,MPI_COMM_WORLD,status,ierr)
         else
            call MPI_SEND(cep_blocks(blk)%vector(1),block_size,
     &           MPI_VAL,wrkr_rank,14,MPI_COMM_WORLD,ierr)
         end if
c
c                               element stresses
c
         block_size = span * ngp * nstrs
         if( worker_processor ) then
            call MPI_RECV(urcs_n_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,0,15,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(urcs_n1_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,0,15,MPI_COMM_WORLD,status,ierr)
         else
            call MPI_SEND(urcs_n_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,wrkr_rank,15,MPI_COMM_WORLD,ierr)
            call MPI_SEND(urcs_n1_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,wrkr_rank,15,MPI_COMM_WORLD,ierr)
         end if
c
c
c                               element integration point rotations
c
         if( rot_blk_list(blk) .ne. 0 ) then
            block_size = span * ngp * 9
            if( worker_processor ) then
               if( allocated(rot_n1_blocks) ) then
                  call MPI_RECV(rot_n1_blocks(blk)%ptr(1),block_size,
     &                 MPI_VAL,0,16,MPI_COMM_WORLD,status,
     &                 ierr)
               end if
            else
               if( allocated(rot_n1_blocks) ) then
                  call MPI_SEND(rot_n1_blocks(blk)%ptr(1),block_size,
     &                 MPI_VAL,wrkr_rank,16,MPI_COMM_WORLD,ierr)
               end if
            end if
         end if
c
c                               element strains
c
         block_size = span * ngp * nstr
         if( worker_processor ) then
            call MPI_RECV(eps_n_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,0,17,MPI_COMM_WORLD,status,ierr)
         else
            call MPI_SEND(eps_n_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,wrkr_rank,17,MPI_COMM_WORLD,ierr)
         end if
c
c                               element volumes
c
         block_size = span
         if( worker_processor ) then
            call MPI_RECV(element_vol_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,0,18,MPI_COMM_WORLD,status,ierr)
         else
            call MPI_SEND(element_vol_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,wrkr_rank,18,MPI_COMM_WORLD,ierr)
         end if
c
c                              send initial state data arrays. we can
c                              delete the blocks for workers not on root
c                              that were needed to support restart.
c
         if( .not. initial_state_option ) cycle
         if( worker_processor ) then
           if( .not. allocated( initial_state_data ) )
     &           allocate( initial_state_data(nelblk) )
             associate( x => initial_state_data )
             allocate( x(blk)%W_plastic_nis_block(span,ngp) )
             call MPI_RECV( x(blk)%W_plastic_nis_block(1,1), span*ngp,
     &           MPI_VAL, 0, 19, MPI_COMM_WORLD, status, ierr )
             end associate
         else
             associate( x => initial_state_data )
             call MPI_send( x(blk)%W_plastic_nis_block(1,1), span*ngp,
     &           MPI_VAL, wrkr_rank, 19,
     &           MPI_COMM_WORLD, status, ierr )
             deallocate( x(blk)%W_plastic_nis_block )
             end associate
         end if
c
      end do ! on blk
c
c                   send the packed vector form of the element
c                   equivalent nodal forces if they exist.
c
      call MPI_BCAST( eq_node_force_len, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      if( eq_node_force_len .gt. 0 ) then
        if( worker_processor ) call mem_allocate( 23 )
        call MPI_BCAST( eq_node_force_indexes, noelem,
     &           MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
        call MPI_BCAST( eq_node_forces, eq_node_force_len,
     &           MPI_VAL, 0, MPI_COMM_WORLD, ierr )
      end if
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine wmpi_send_jint                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/28/2018 rhd              *
c     *                                                              *
c     *           send/receive data to processes to set up           *
c     *           j-integral and i-integral computation for          *
c     *           this single domain                                 *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_jint
      use global_data, only : numprocs, nonode, noelem, out, myid,
     &                 root_processor,  worker_processor, MPI_VAL
      use j_data, only : num_front_nodes, max_exp_front, q_values,
     &    q_element_maps, front_element_list, expanded_front_nodes,
     &    first_domain, front_nodes, front_order, front_element_list,
     &    domain_type, omit_crack_front_elems, face_loading,
     &    symmetric_domain, q_vals_linear, one_point_rule, static_j,
     &    ignore_face_loads, front_list_length,
     &    debug_driver, debug_elements, print_elem_values, out_pstress,
     &    out_pstrain, cf_traction_flags, domain_rot, domain_origin,
     &    cf_tractions, front_coords, e_front, nu_front,
     &    crack_curvature, max_front_nodes, j_geonl_formulation,
     &    j_linear_formulation, temperatures_on_model,
     &    process_initial_state, process_temperatures
c
      implicit none
      include "mpif.h"
c
      integer :: ierr, flags(4), num
      logical :: logical_vec(17)
      logical, parameter :: debug = .false.
c
      if( numprocs .eq. 1 ) return
c
c         send /receive shared data needed to process this domain.
c         we have two major cases: (1) this is the first domain
c         to be processed and most all data needs to be broadcast,
c         (2) on subsequent domains, we only need to broadcast
c         domain dependent information.
c
c         broadcast domain dependent information. these 2 bcasts just make
c         logic simpler
c
      call MPI_BCAST( num_front_nodes, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( max_exp_front, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( front_list_length, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD,  ierr )
c
      flags = 0
      if( worker_processor ) then
       if( allocated( q_values ) ) deallocate( q_values )
       allocate( q_values(nonode), stat=flags(1))
       if( allocated( q_element_maps ) ) deallocate( q_element_maps )
       allocate( q_element_maps(noelem/30+1),stat=flags(2) )
       if( allocated( front_element_list) )
     &     deallocate( front_element_list )
       allocate( front_element_list(front_list_length), stat=flags(3) )
       if( allocated( expanded_front_nodes ) )
     &     deallocate( expanded_front_nodes )
       allocate( expanded_front_nodes(0:max_exp_front,num_front_nodes),
     &           stat=flags(4) )
       if( any( flags .ne. 0 ) ) then
          write(out,9100); call die_abort
       end if
      end if
c
      call MPI_BCAST( q_values, nonode, MPI_REAL, 0,
     &                MPI_COMM_WORLD, ierr )
c
      call MPI_BCAST( q_element_maps, noelem/30+1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
c
      call MPI_BCAST( front_element_list, front_list_length,
     &                 MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
c
      call MPI_BCAST( first_domain, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     &                ierr )
c
      call MPI_BCAST( front_nodes, max_front_nodes, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
c
      call MPI_BCAST( num_front_nodes, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
c
      call MPI_BCAST( front_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &                ierr )
c
      call MPI_BCAST( domain_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &                ierr )
c
      call MPI_BCAST( omit_crack_front_elems, 1, MPI_LOGICAL, 0,
     &                MPI_COMM_WORLD, ierr )
c
      call MPI_BCAST( max_exp_front, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &                ierr )
c
      call MPI_BCAST( face_loading, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     &                ierr )
c
      if( .not. first_domain ) return
c
c         processing first domain. broadcast data used for all
c         domains defined at this crack front position.
c
c         a set of logical flags. for workers turn off flags
c         that would otherwise let them print intermediate results
c         during j computation.
c
      if( root_processor ) then
         logical_vec(1)  = symmetric_domain
         logical_vec(2)  = q_vals_linear
         logical_vec(3)  = one_point_rule
         logical_vec(4)  = static_j
         logical_vec(5)  = ignore_face_loads
         logical_vec(6)  = j_linear_formulation
         logical_vec(7)  = j_geonl_formulation
         logical_vec(8)  = debug_driver
         logical_vec(9)  = debug_elements
         logical_vec(10) = print_elem_values
         logical_vec(11) = out_pstress
         logical_vec(12) = out_pstrain
         logical_vec(13) = cf_traction_flags(1)
         logical_vec(14) = cf_traction_flags(2)
         logical_vec(15) = cf_traction_flags(3)
         logical_vec(16) = temperatures_on_model
         logical_vec(17) = process_initial_state
         logical_vec(18) = process_temperatures
      end if
c
      call MPI_BCAST( logical_vec, 18, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     &                ierr )
c
      if( worker_processor ) then
         symmetric_domain     = logical_vec(1)
         q_vals_linear        = logical_vec(2)
         one_point_rule       = logical_vec(3)
         static_j             = logical_vec(4)
         ignore_face_loads    = logical_vec(5)
         j_linear_formulation = logical_vec(6)
         j_geonl_formulation  = logical_vec(7)
         debug_driver         = logical_vec(8)
         debug_elements       = logical_vec(9)
         print_elem_values    = logical_vec(10)
         out_pstress          = logical_vec(11)
         out_pstrain          = logical_vec(12)
         cf_traction_flags(1) = logical_vec(13)
         cf_traction_flags(2) = logical_vec(14)
         cf_traction_flags(3) = logical_vec(15)
         temperatures_on_model = logical_vec(16)
         process_initial_state = logical_vec(17)
         process_temperatures  = logical_vec(18)
      end if
c
c             1) 3x3 rotation matrix for this crack front poisiton
c
      call MPI_BCAST( domain_rot, 9, MPI_VAL, 0, MPI_COMM_WORLD,
     &                ierr )
c
      call MPI_BCAST( domain_origin, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
c
      call MPI_BCAST( cf_tractions, 3, MPI_VAL, 0,
     &                MPI_COMM_WORLD, ierr )
c
      call MPI_BCAST( front_coords, 3*max_front_nodes, MPI_VAL, 0,
     &                MPI_COMM_WORLD,ierr)
c
      num = (max_exp_front + 1) * num_front_nodes
      call MPI_BCAST( expanded_front_nodes, num, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
c
      call MPI_BCAST( e_front, 1, MPI_VAL, 0, MPI_COMM_WORLD, ierr )
c
      call MPI_BCAST( nu_front, 1, MPI_VAL, 0, MPI_COMM_WORLD, ierr )
c
      call MPI_BCAST( crack_curvature, 7, MPI_VAL, 0,
     &                MPI_COMM_WORLD, ierr )
c
      if( debug )  then
        write (*,*) '<<< just leaving wmpi_send_jint,', myid
      end if
c
      return
c
 9100 format('>> FATAL ERROR: allocate error in wmpi_send_jint: ')
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *                subroutine wmpi_reduce_real_max               *
c     *                                                              *
c     *                       written by :rhd                        *
c     *                                                              *
c     *                   last modified : 6/29/2018 rhd              *
c     *                                                              *
c     *     reduction operator of single precision to root with      *
c     *     MPI_MAX rather than MPI_SUM. each final entry on root    *
c     *     is the max value across workers                          *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_reduce_real_max( vector, length )
      use global_data ! old common.main
      implicit none
      include "mpif.h"
c
      integer :: length
      real :: vector(length)
c                 Local
      real, dimension(:), allocatable :: temp_vec
      integer :: ierr
c
c                 Setup receive buffer, reduce via *MAX*, copy final
c                 vector to actual root vector
c
      allocate( temp_vec(length) )
      call MPI_REDUCE( vector, temp_vec, length, MPI_REAL,
     &                 MPI_MAX, 0, MPI_COMM_WORLD, ierr )
      if( root_processor ) vector(1:length) = temp_vec(1:length)
c
      return
c
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *              subroutine wmpi_send_simple_angles              *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 4/26/2017 rhd              *
c     *                                                              *
c     *           Send the simple angle properties, if required      *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_simple_angles
      use global_data ! old common.main
      use crystal_data, only: srequired, simple_angles, nangles,
     &                        mc_array
      implicit none
      include 'mpif.h'
      integer :: ierr
c
c              get the workers synced to this point
c
      if( root_processor ) call wmpi_alert_slaves(52)
c
      call MPI_Bcast( srequired, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     &                ierr )
c
      if( .not. srequired ) return
c
      call MPI_Bcast( nangles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &                ierr )
c
      if( worker_processor ) then
        if( .not. allocated(simple_angles) )
     &             allocate(simple_angles(nangles,3))
        if( .not. allocated(mc_array) )
     &             allocate(mc_array(noelem,max_crystals))
      end if
c
      call MPI_Bcast( simple_angles, 3*nangles, MPI_DOUBLE_PRECISION,
     &                0, MPI_COMM_WORLD, ierr )
      call MPI_Bcast( mc_array, noelem*max_crystals, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
c
      return
c
      end subroutine

c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_crystals           *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 5/12/2017 rhd              *
c     *                                                              *
c     *           Send all the crystal properties to the workers     *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_crystals
      use global_data ! old common.main
      use crystal_data, only: c_array, angle_input,
     &                        crystal_input, data_offset, print_crystal
      implicit none
      include 'mpif.h'
      integer, parameter :: count_struct=9
      integer ::  nelem, mxcry, ierr, i, blocklens(0:8), indices(0:8),
     &            otypes(0:8), size_int, size_dp, size_log, ntype
      logical, parameter :: local_debug = .false.
c
      if( myid .eq. 0 ) then
          call wmpi_alert_slaves( 49 )
          nelem = size(crystal_input,1) ! num rows
          mxcry = size(crystal_input,2) ! num cols
      end if

      call MPI_Bcast( nelem, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr)
      call MPI_Bcast( mxcry, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr)
      call MPI_Bcast( noelem, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr)

      if( nelem > 0 ) then
c
        if( worker_processor ) then
          if( .not. allocated(data_offset) )
     &          allocate( data_offset(noelem) )
          if( .not. allocated(angle_input) )
     &          allocate( angle_input(nelem,mxcry,3) )
          if( .not. allocated(crystal_input) )
     &          allocate( crystal_input(nelem,mxcry) )
        end if
        call MPI_Bcast( data_offset, noelem, MPI_INTEGER, 0,
     &                  MPI_COMM_WORLD, ierr )
        call MPI_Bcast( angle_input, nelem*mxcry*3,
     &                  MPI_DOUBLE_PRECISION, 0,
     &                  MPI_COMM_WORLD, ierr )
        call MPI_Bcast( crystal_input, nelem*mxcry,
     &                  MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
c
      end if
c
      call MPI_Type_extent( MPI_INTEGER, size_int, ierr )
      call MPI_Type_extent( MPI_DOUBLE_PRECISION, size_dp, ierr )
      call MPI_Type_extent( MPI_LOGICAL, size_log, ierr )
c
      otypes(0) = MPI_DOUBLE_PRECISION
      otypes(1) = MPI_DOUBLE_PRECISION
      otypes(2) = MPI_DOUBLE_PRECISION
      otypes(3) = MPI_DOUBLE_PRECISION
      otypes(4) = MPI_DOUBLE_PRECISION
      otypes(5) = MPI_DOUBLE_PRECISION
      otypes(6) = MPI_INTEGER
      otypes(7) = MPI_DOUBLE_PRECISION
      otypes(8) = MPI_LOGICAL
c
      blocklens(0) = 2*6*6
      blocklens(1) = 2*max_slip_sys*3
      blocklens(2) = 27
      blocklens(3) = 10
      blocklens(4) = 6
      blocklens(5) = 6
      blocklens(6) = 12
      blocklens(7) = 100
      blocklens(8) = 5
c
      indices(0) = 0
      indices(1) = 2*6*6*size_dp
      indices(2) = indices(1) + 2*max_slip_sys*3*size_dp
      indices(3) = indices(2) + 27*size_dp
      indices(4) = indices(3) + 10*size_dp
      indices(5) = indices(4) + 6*size_dp
      indices(6) = indices(5) + 6*size_dp
      indices(7) = indices(6) + 12*size_int
      indices(8) = indices(7) + 100*size_dp
c
      call MPI_Type_struct( count_struct,  blocklens,  indices,
     &                      otypes,  ntype,  ierr )
      call MPI_Type_commit( ntype,  ierr )
      call MPI_Barrier( MPI_COMM_WORLD,  ierr )
c
      if( local_debug ) then
       do i=0,3
         if (i .eq. myid) then
           write (*,*) "..............................."
           write (*,*) myid
           write (*,*)
           call print_crystal(1)
           write (*,*)
         end if
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
       end do
      end if
c
      if( local_debug )
     &       write(*,*) '...   transmitting c_array...',myid
c
       call MPI_Bcast( c_array, max_crystals, ntype,
     &                 0, MPI_COMM_WORLD, ierr )
c
      if( local_debug )
     &  write(*,*) '...   done transmitting c_array...',myid
      call MPI_Type_free( ntype, ierr )
      call MPI_Barrier( MPI_COMM_WORLD, ierr )
c
      if( local_debug ) then
           do i=0,3
           if (i .eq. myid) then
           write (*,*) "..............................."
           write (*,*) myid
           write (*,*)
           call print_crystal(1)
           write (*,*)
           end if
           call MPI_Barrier(MPI_COMM_WORLD,ierr)
           end do
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
          write(*,*) '......   at die gracefully debug'
          call die_gracefully
      end if
c
      return
c
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *               subroutine wmpi_dealloc_crystals               *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 4/26/2017 rhd              *
c     *                                                              *
c     *                Dealloc crystal data structures               *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_dealloc_crystals
      use global_data ! old common.main
      use crystal_data, only: c_array, angle_input, crystal_input,
     &                        data_offset, mc_array, simple_angles
      implicit none
      include 'mpif.h'
c
      if( root_processor ) then
         call wmpi_alert_slaves( 50 )
      else
         if( allocated(angle_input) )   deallocate( angle_input )
         if( allocated(crystal_input) ) deallocate( crystal_input )
         if( allocated(data_offset) )   deallocate( data_offset )
         if( allocated(simple_angles) ) deallocate( simple_angles )
         if( allocated(mc_array) )      deallocate( mc_array )
      end if
      return
c
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *            subroutine wmpi_compute_set_history_locs          *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified : 5/12/2017 tjt, rhd         *
c     *                                                              *
c     *           Send CP history sizing data to workers             *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_compute_set_history_locs
      use global_data ! old common.main
      use mm10_defs, only : indexes_common, index_crys_hist,
     &                      num_common_indexes, num_crystal_terms,
     &                      length_crys_hist, one_crystal_hist_size,
     &                      common_hist_size, length_comm_hist
c
      implicit none
      include 'mpif.h'
      integer :: ierr, k
      logical :: local_debug
c
      local_debug = .false.
c
c              get the workers synced to this point
c
      if( root_processor ) call wmpi_alert_slaves( 53 )
c
      call MPI_Bcast( num_common_indexes, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_Bcast( num_crystal_terms, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_Bcast( common_hist_size, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_Bcast( one_crystal_hist_size, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
c
      if( worker_processor ) then
c
        if( .not. allocated( indexes_common ) )
     &     allocate( indexes_common(num_common_indexes,2) )
        if( .not. allocated( index_crys_hist ) )
     &     allocate( index_crys_hist(max_crystals,
     &               num_crystal_terms,2) )
        if( .not. allocated( length_comm_hist ) )
     &     allocate( length_comm_hist(num_common_indexes) )
        if( .not. allocated( length_crys_hist ) )
     &     allocate( length_crys_hist(num_crystal_terms) )
      end if
c
      call MPI_Bcast( indexes_common, num_common_indexes*2,
     &                MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
      k = max_crystals*num_crystal_terms*2
      call MPI_Bcast( index_crys_hist, k, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_Bcast( length_comm_hist, num_common_indexes,
     &                MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
      call MPI_Bcast( length_crys_hist, num_crystal_terms,
     &                MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
c
      return
c
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_init_owner              *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 2/25/2017 rhd              *
c     *                                                              *
c     *                     generates node owner structure           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_init_owner
      use global_data ! old common.main
      use main_data, only: elems_to_blocks, inverse_incidences
      use mpi_lnpcg
      implicit none
      include "mpif.h"
c
      integer :: ierr, alloc_stat, proc, i, j, node, own_proc,
     &           blk,  new_proc, start, elem
      integer, allocatable, dimension(:,:) ::  node2proc
      integer :: num_private_nodes(0:max_procs-1),
     &     num_own_shared_nodes(0:max_procs-1),
     &     num_shared_nodes(0:max_procs-1),
     &     num_local_nodes(0:max_procs-1),
     &     count_priv(0:max_procs-1), count_own(0:max_procs-1),
     &     count_shar(0:max_procs-1),
     &     conn(max_procs,0:max_procs-1), conn_num(0:max_procs-1),
     &     ordering(max_procs,max_procs), order_num(max_procs)
      logical, save :: been_here = .false.
      logical, parameter :: debug = .false.
      logical :: blk_on_boundary(mxnmbl)
      real, external ::  wcputime
      character(len=80) :: line
c
      if ( debug ) write (out,*) myid,': ->> in wmpi_init_owner'
c
      if( been_here ) return ! do this only once per run
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
      allocate( node2proc(0:numprocs,nonode), stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
      node2proc(0:numprocs,1:nonode) = 0
c
c            set up data structures needed to calculate
c            the owner data structures
c
      call wmpi_owner_setup( node2proc, num_private_nodes,
     &     num_own_shared_nodes, num_shared_nodes, num_local_nodes,
     &     blk_on_boundary )
c
      if( debug ) write (out,*) '...allocating'
c
c            allocate space on root proc for all the info. We will
c            deallocate this space once all processors have been
c            sent their information.
c
      if( allocated(proc_nodes) ) deallocate( proc_nodes )
      allocate(proc_nodes(0:numprocs-1), stat = alloc_stat)
      if( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
      do proc = 0, numprocs - 1
c
c              allocate data structure for each processor
c
         call wmpi_allocate_owner( num_private_nodes(proc),
     &        num_own_shared_nodes(proc), num_shared_nodes(proc),
     &        nelblk, proc_nodes(proc) )
c
c              allocate vectors to hold lists of shared nodes for
c              each processor.  Here we allocate these vectors to
c              be larger than we need; since this is a temporary
c              data structure, this is not a problem.
c
         do i = 0, numprocs -1
c
            allocate( proc_nodes(proc)%sharing_procs(i)%ptr(
     &           num_own_shared_nodes(proc)), stat = alloc_stat )
            if( alloc_stat .ne. 0 ) then
                write(out,9900)
               call die_abort
               stop
            end if
c
            allocate( proc_nodes(proc)%shared_owner(i)%ptr(
     &           num_shared_nodes(proc)), stat = alloc_stat)
            if( alloc_stat .ne. 0 ) then
                write(out,9900)
               call die_abort
               stop
            end if
c
         end do  ! on i procs
c
      end do ! on proc procs
c
c            now fill the data structures
c
      if( debug ) write (out,*) '...filling structure'
      do node = 1, nonode
         own_proc = node2proc(1,node)
c
c               if there is only one processor referencing this node, then
c               this node is private to that processor.  increment the number
c               of private nodes by one for that processor
c
         if( node2proc(0,node) == 1 ) then
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
            do i = 2, node2proc(0,node)
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
            end do  ! on i
         endif
      end do ! on node
c
c            in the future we may want to know which
c            blocks are completely internal (have no shared
c            nodes). this was part of the ebe pcg solver in the
c            past.
c
c               First find out which blocks
c               have elements on the boundary between two
c               processors, and which do not.
c
      blk_on_boundary(1:nelblk) = .false.
c
      do node = 1, nonode
         if( node2proc(0,node) == 1 ) cycle
         do j = 1, inverse_incidences(node)%element_count
            elem = inverse_incidences(node)%element_list(j)
            blk  = elems_to_blocks(elem, 1)
            blk_on_boundary(blk) = .true.
         end do
      end do
c
c               fill list of internal blocks for each processor
c
      do proc = 0, numprocs - 1
         proc_nodes(proc)%num_int_blks = 0
      end do
      do blk = 1, nelblk
         if( blk_on_boundary(blk) ) cycle
         proc = elblks(2,blk)
         proc_nodes(proc)%num_int_blks =
     &        proc_nodes(proc)%num_int_blks + 1
         proc_nodes(proc)%internal_blks(proc_nodes(proc)%num_int_blks) =
     &        blk
      end do
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
         do i = 1, proc_nodes(proc)%num_private
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
         end do
c
c                  process the nodes shared and owned by the current processor
c
         start = proc_nodes(proc)%num_private * 3
         do i = 1, proc_nodes(proc)%num_own_shared
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
         end do
c
c                  process the nodes shared but not owned by the processor
c
         start = ( proc_nodes(proc)%num_private +
     &        proc_nodes(proc)%num_own_shared ) * 3
         do i = 1, proc_nodes(proc)%num_shared
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
         end do
      end do
c
c            print out the count for private and shared nodes, and the
c            number shared with each other processor.
c
      write(out,'("  -------------------------")')
      write(out,'("  Node Ownership Statistics")')
      write(out,'("  -------------------------")')
      do i = 0, numprocs-1
         call wmpi_print_node_own( proc_nodes(i), i )
      end do
      write(out,'("  -------------------------")')
c
c            fill processor zero's own local_nodes structure -- this
c            is the permanent data stucture for the node ownership
c            data for the root processor
c
      if( debug ) write(out,*) '...do proc 0s local_nodes'
      call wmpi_allocate_owner( num_private_nodes(0),
     &     num_own_shared_nodes(0), num_shared_nodes(0),
     &     proc_nodes(0)%num_int_blks, local_nodes )
c
      call wmpi_copy_own( proc_nodes(0), local_nodes )
c
c            for the root processor, store all the mappings from
c            processor local to global degrees of freedom
c
c               first allocate the base data structure
c
      allocate( procdof2glob(numprocs-1), stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
c               now allocate and fill the mappings for each slave proc
c
      do proc = 1, numprocs - 1
          procdof2glob(proc)%num_dof = num_local_nodes(proc)*mxndof
          allocate(
     &       procdof2glob(proc)%dof(num_local_nodes(proc)*mxndof),
     &       stat = alloc_stat )
         if( alloc_stat .ne. 0 ) then
            write(out,9900)
            call die_abort
            stop
         end if
         do i = 1, procdof2glob(proc)%num_dof
            procdof2glob(proc)%dof(i) = proc_nodes(proc)%local2global(i)
         end do
      end do
c
c            for each processor, have the processor allocate
c            the appropriate data structure, then have the
c            root processor send the correct information to
c            the processor.
c
      if( debug ) write (out,*) '...doing multi-processor thang'
      call wmpi_owner_send( num_private_nodes,
     &     num_own_shared_nodes, num_shared_nodes, num_local_nodes )
c
c            deallocate all the structures that are no longer needed
c
      if( debug ) write(out,*) '...deallocate'
      deallocate( node2proc )
c
      do proc = 0, numprocs - 1
         do i = 0, numprocs - 1
            deallocate(proc_nodes(proc)%sharing_procs(i)%ptr)
            deallocate(proc_nodes(proc)%shared_owner(i)%ptr)
         end do
      end do
c
      do proc = 0, numprocs - 1
         call wmpi_deallocate_owner( proc_nodes(proc) )
      end do
      deallocate( proc_nodes )
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
c     *                   last modified : 02/26/2017 rhd             *
c     *                                                              *
c     *                        generates node owner structure        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_owner_setup(node2proc, num_private_nodes,
     &     num_own_shared_nodes, num_shared_nodes, num_local_nodes,
     &     blk_on_boundary )
      use global_data ! old common.main
c
      use main_data, only: inverse_incidences
      implicit none
      include "mpif.h"
c
      integer :: ierr, blk, felem, span, elem, node, elem_loop, proc,
     &           i, j, tmp
      integer :: num_private_nodes(0:max_procs-1),
     &     num_own_shared_nodes(0:max_procs-1),
     &     num_shared_nodes(0:max_procs-1),
     &     num_local_nodes(0:max_procs-1),
     &     elem2blk(mxel), node2proc(0:numprocs,1:nonode)
      logical :: its_new_proc, blk_on_boundary(*)
      logical, parameter :: debug = .false.
      character(len=80) :: line
c
c          build element to block list
c
      do blk = 1, nelblk
         blk_on_boundary(blk) = .false.
         felem          = elblks(1,blk)
         span           = elblks(0,blk)
         do elem = felem, felem + span - 1
            elem2blk(elem) = blk
         end do
      end do
c
c          build structure which, when given a specific node,
c          tells what processors access it.
c
      do node = 1, nonode
         do elem_loop = 1, inverse_incidences(node)%element_count
            elem = inverse_incidences(node)%element_list(elem_loop)
            blk = elem2blk(elem)
            proc = elblks(2,blk)
            if( node2proc(0,node) == 0 ) then
               node2proc(0,node) = node2proc(0,node) + 1
               node2proc(1,node) = proc
               cycle
            else
               its_new_proc = .true.
               do i = 1, node2proc(0,node)
                  if( node2proc(i,node) == proc ) its_new_proc=.false.
               end do
               if( .not. its_new_proc ) cycle
               node2proc(0,node) = node2proc(0,node) + 1
               node2proc(node2proc(0,node),node) = proc
            endif
c
c          if node is even, make lowest number processor
c          its owner.  The owner is the first processor
c          in the list.  Thus make sure lowest number
c          processor is first in list
c
            if( mod(node,2) == 0 ) then
               if( node2proc(node2proc(0,node),node) <
     &              node2proc(1,node) ) then
                  tmp = node2proc(node2proc(0,node),node)
                  node2proc(node2proc(0,node),node) = node2proc(1,node)
                  node2proc(1,node) = tmp
               endif
            else
c
c          if node is odd, make highest number processor its owner.
c
               if( node2proc(node2proc(0,node),node) >
     &              node2proc(1,node) ) then
                  tmp = node2proc(node2proc(0,node),node)
                  node2proc(node2proc(0,node),node) = node2proc(1,node)
                  node2proc(1,node) = tmp
               endif
            endif
c
         end do ! on elem_loop
      end do ! on node
c
c          write out structure for checking
c
      if( debug ) then
         write(out,*) '>>> proc accesses of nodes:'
         do i = 1, nonode
            write(out,'(4x,i7,3x,10i3)') i, (node2proc(j,i),j=1,
     &           node2proc(0,i) )
         end do
      endif
c
c          count number of nodes for each processor of
c          the following types:
c           1:  private -- only one processor references the node
c           2:  own_shared -- owned by processor, but shared with
c                           others
c           3:  shared -- owned by other processor, but referenced
c                            by current processor
c
      if( debug ) write (out,*) '...counting'
      do node = 1, nonode
         proc = node2proc(1,node)
c
c          if there is only one processor referencing this node, then
c          this node is private to that processor.  increment the number
c          of private nodes by one for that processor
c
         if( node2proc(0,node) == 1 ) then
            num_private_nodes(proc) = num_private_nodes(proc) + 1
         else
c
c          there are more than 1 processors which reference the node.
c          for the owning processor, increase the number of own_shared
c          nodes by one. Also, for all other processors that reference
c          the node but don't own it, increase their number of shared
c          nodes.
c
            num_own_shared_nodes(proc) = num_own_shared_nodes(proc) + 1
            do i = 2, node2proc(0,node)
               num_shared_nodes(node2proc(i,node)) =
     &              num_shared_nodes(node2proc(i,node)) + 1
            end do
         endif
      end do
c
      do proc = 0, numprocs - 1
         num_local_nodes(proc) = num_private_nodes(proc) +
     &        num_own_shared_nodes(proc) +  num_shared_nodes(proc)
      end do
c
      if( debug ) then
         write(out,*) ' >> proc:   private  own_shared   shared'
         do proc = 0, numprocs - 1
            write(out,'(2x,i5,3(3x,i7))')
     &           proc, num_private_nodes(proc),
     &           num_own_shared_nodes(proc), num_shared_nodes(proc)
         end do
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
c     *                   last modified : 02/27/2017 rhd             *
c     *                                                              *
c     *     send node owner structures to workers                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_owner_send( num_private_nodes,
     &     num_own_shared_nodes, num_shared_nodes, num_local_nodes )
      use global_data ! old common.main
      use mpi_lnpcg, only: local_nodes, proc_nodes
      implicit none
      include "mpif.h"
c
      integer :: ierr, status(MPI_STATUS_SIZE), proc, i, num_private,
     &           num_own_shared, num_shared
      integer :: num_private_nodes(0:max_procs-1),
     &     num_own_shared_nodes(0:max_procs-1),
     &     num_shared_nodes(0:max_procs-1),
     &     num_local_nodes(0:max_procs-1)
      integer, allocatable :: disp(:), blksz(:)
      logical, parameter ::  debug = .false.
      character(len=80) :: line
c
      if( debug ) write(out,*) myid,': ->> in wmpi_owner_send'
c
      allocate( disp(3*nonode), blksz(3*nonode) )
c
c            notify workers
c
      call wmpi_alert_slaves ( 16 )
c
c            send the data
c
      if( root_processor ) then
c
         do proc = 1, numprocs - 1
c
c              send out the constants first
c
            call MPI_SEND(num_private_nodes(proc),1,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(num_own_shared_nodes(proc),1,MPI_INTEGER,
     &           proc,proc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(num_shared_nodes(proc),1,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(proc_nodes(proc)%num_int_blks,1,MPI_INTEGER,
     &           proc,proc,MPI_COMM_WORLD,ierr)
c
c              now send out the arrays
c
c                 send out information about local to/from global mapping
c
            call MPI_SEND(proc_nodes(proc)%local2global,
     &           num_local_nodes(proc)*mxndof,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(proc_nodes(proc)%global2local,
     &           nodof,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
c
c                 send out information about private nodes
c
            call MPI_SEND(proc_nodes(proc)%private,
     &           num_private_nodes(proc),MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
c
c                 send out information about shared but owned nodes
c
            call MPI_SEND(proc_nodes(proc)%own_shared,
     &           num_own_shared_nodes(proc),MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(proc_nodes(proc)%sharing_count,
     &           numprocs,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
c
c                    send out the pointers which indicate which shared and
c                    owned node are shared with each processor
c
            do i = 0, numprocs-1
               if( proc_nodes(proc)%sharing_count(i) == 0 ) cycle
               call MPI_SEND(proc_nodes(proc)%sharing_procs(i)%ptr,
     &              proc_nodes(proc)%sharing_count(i),MPI_INTEGER,proc,
     &              proc,MPI_COMM_WORLD,ierr)
            end do
c
c                 send out information about shared and not owned nodes
c
            call MPI_SEND(proc_nodes(proc)%shared,
     &           num_shared_nodes(proc),MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
            call MPI_SEND(proc_nodes(proc)%shared_count,
     &           numprocs,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
c
c                    send out the pointers which indicate which shared but
c                    not owned node are owned by, num_shared each processor
c
            do i = 0, numprocs-1
               if( proc_nodes(proc)%shared_count(i) == 0 ) cycle
               call MPI_SEND(proc_nodes(proc)%shared_owner(i)%ptr,
     &              proc_nodes(proc)%shared_count(i),MPI_INTEGER,proc,
     &              proc,MPI_COMM_WORLD,ierr)
            end do
c
c                 send out information about internal blocks
c
            call MPI_SEND(proc_nodes(proc)%internal_blks,
     &           proc_nodes(proc)%num_int_blks,MPI_INTEGER,proc,
     &           proc,MPI_COMM_WORLD,ierr)
c
         end do ! over procs
c
      else
c
c              get the constants
c
         call MPI_RECV(num_private,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,
     &        status, ierr)
         call MPI_RECV(num_own_shared,1,MPI_INTEGER,0,myid,
     &        MPI_COMM_WORLD, status, ierr)
         call MPI_RECV(num_shared,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,
     &        status, ierr)
         call MPI_RECV(local_nodes%num_int_blks,1,MPI_INTEGER,0,myid,
     &        MPI_COMM_WORLD, status, ierr)
c
c              allocate the space
c
         call wmpi_allocate_owner(num_private, num_own_shared,
     &              num_shared, local_nodes%num_int_blks, local_nodes)
c
c              now fill the arrays
c
c                 recv info about local to/from global mapping
c
         call MPI_RECV(local_nodes%local2global,
     &        local_nodes%num_local_nodes*mxndof,MPI_INTEGER,0,
     &        myid,MPI_COMM_WORLD,status,ierr)
         call MPI_RECV(local_nodes%global2local,
     &        nodof,MPI_INTEGER,0,
     &        myid,MPI_COMM_WORLD,status,ierr)
c
c                 recv info about private nodes
c
         call MPI_RECV(local_nodes%private,
     &        num_private,MPI_INTEGER,0,
     &        myid,MPI_COMM_WORLD,status,ierr)
c
c                 recv info about shared but owned nodes
c
         call MPI_RECV(local_nodes%own_shared,
     &        num_own_shared,MPI_INTEGER,0,
     &        myid,MPI_COMM_WORLD,status,ierr)
         call MPI_RECV(local_nodes%sharing_count,numprocs,
     &        MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
c
c                    recv the pointers which indicate which shared and
c                    owned node are shared by each processor.  Allocate
c                    each pointer vector when we know how big it should be.
c
         do i = 0, numprocs - 1
            if( local_nodes%sharing_count(i) == 0 ) cycle
            allocate(local_nodes%sharing_procs(i)%ptr(
     &           local_nodes%sharing_count(i)))
            call MPI_RECV(local_nodes%sharing_procs(i)%ptr,
     &           local_nodes%sharing_count(i),MPI_INTEGER,0,
     &           myid,MPI_COMM_WORLD,status,ierr)
         end do
c
c                 recv info about shared but not owned nodes
c
         call MPI_RECV(local_nodes%shared,
     &        num_shared,MPI_INTEGER,0,
     &        myid,MPI_COMM_WORLD,status,ierr)
         call MPI_RECV(local_nodes%shared_count,numprocs,
     &        MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
c
c                    recv the pointers which indicate which shared but
c                    not owned node are owned by each processor.  Allocate
c                    each pointer vector when we know how big it should be.
c
         do i = 0, numprocs-1
            if( local_nodes%shared_count(i) == 0 ) cycle
            allocate(local_nodes%shared_owner(i)%ptr(
     &           local_nodes%shared_count(i)))
            call MPI_RECV(local_nodes%shared_owner(i)%ptr,
     &           local_nodes%shared_count(i),MPI_INTEGER,0,
     &           myid,MPI_COMM_WORLD,status,ierr)
         end do
c
c                 recv info about internal blocks
c
         call MPI_RECV(local_nodes%internal_blks,
     &        local_nodes%num_int_blks,
     &        MPI_INTEGER,0,myid,MPI_COMM_WORLD,status,ierr)
c
      endif
c
c           print out the data structures if debug is on.
c
      if( debug ) then
         do proc = 0, numprocs - 1
            if( myid == proc ) then
               write(out,*) '>=>=>=>=>=> This is proc', myid
               call wmpi_print_node_own_all(local_nodes)
            endif
            call MPI_BARRIER (MPI_COMM_WORLD, ierr)
         end do
      endif
c
c           create data types for transmission of the extra nodes
c
      do proc = 0, numprocs - 1
c
         if( local_nodes%sharing_count(proc) > 0 ) then
c
            do i = 1, local_nodes%sharing_count(proc)
               disp(i) = (local_nodes%sharing_procs(proc)%ptr(i)-1) * 3
               blksz(i) = 3
            end do
            call MPI_TYPE_INDEXED( local_nodes%sharing_count(proc),
     &           blksz,disp,MPI_VAL, local_nodes%MPI_sharing_type(proc),
     &           ierr )
            call MPI_TYPE_COMMIT( local_nodes%MPI_sharing_type(proc),
     &           ierr )
c
         endif
c
         if( local_nodes%shared_count(proc) > 0 ) then
c
            do i = 1, local_nodes%shared_count(proc)
               disp(i) = (local_nodes%shared_owner(proc)%ptr(i)-1) * 3
               blksz(i) = 3
            end do
            call MPI_TYPE_INDEXED( local_nodes%shared_count(proc),
     &           blksz,disp,MPI_VAL, local_nodes%MPI_shared_type(proc),
     &           ierr )
            call MPI_TYPE_COMMIT( local_nodes%MPI_shared_type(proc),
     &           ierr )
c
         endif
c
      end do
c
c                  create a new local version of edest_blocks which points
c                  directly to the locally-renumbers dofs instead of the
c                  global dofs.
c
c      call ledest_init
c
c                  deprecated allocate kdiag, local_mdiag
c
      if( debug ) write(out,*) myid,': <<- leaving wmpi_owner_send'
      return
 9900 format('>>> FATAL ERROR: memory allocate failure...')
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine wmpi_allocate_owner               *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/27/2017 rhd             *
c     *                                                              *
c     *                   creates  node owner structure              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_allocate_owner( num_priv, num_own_shar,
     &                                num_shar, num_int_blks, own )
      use global_data ! old common.main
      use mpi_lnpcg
      implicit none
      include "mpif.h"
c
      integer :: num_priv, num_own_shar,  num_shar, num_int_blks
      type (node_owner_type) :: own
c
      integer :: ierr, ptr, status(MPI_STATUS_SIZE), alloc_stat
c
c         allocate pointer from local dofs to global dofs
c
      own%num_local_nodes = num_priv + num_own_shar + num_shar
      allocate( own%local2global(own%num_local_nodes*mxndof),
     &     stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
      own%local2global(1:own%num_local_nodes*mxndof) = 0
      allocate( own%global2local(nodof), stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
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
      allocate( own%private(num_priv), stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
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
      allocate( own%own_shared(num_own_shar), stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
      own%own_shared(1:num_own_shar) = 0
c
      allocate( own%sharing_count(0:numprocs-1), stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
      own%sharing_count(0:numprocs-1) = 0
c
c            allocate ptr vector which indicates which localled owned but
c            shared nodes are shared with each processor
c
      allocate( own%sharing_procs(0:numprocs-1), stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
      allocate( own%MPI_sharing_type(0:numprocs-1), stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
c         allocate list of shared nodes -- nodes referenced by this
c         processor , but owned by another
c
      own%num_shared = num_shar
      allocate( own%shared(num_shar), stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
      own%shared(1:num_shar) = 0
c
      allocate( own%shared_count(0:numprocs-1), stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
      own%shared_count(0:numprocs-1) = 0
c
c            allocate ptr vector which indicates which shared but not
c            shared nodes are shared with each processor
c
      allocate( own%shared_owner(0:numprocs-1), stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
      allocate( own%MPI_shared_type(0:numprocs-1), stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
         write(out,9900)
         call die_abort
         stop
      end if
c
c         allocate list of blocks with no elements on boundary.
c
      own%num_int_blks = num_int_blks
      allocate( own%internal_blks(num_int_blks), stat = alloc_stat )
      if( alloc_stat .ne. 0 ) then
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
c     *                 subroutine wmpi_deallocate_owner             *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 2/23/2017 rhd              *
c     *                                                              *
c     *     support for dallocation of node owner setup              *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_deallocate_owner(own)
      use mpi_lnpcg
      implicit none
      include "mpif.h"
      type (node_owner_type) :: own
c
c                       local declarations
c
      integer :: ierr, status(MPI_STATUS_SIZE)
c
c         allocate pointer from local dofs to global dofs
c
      deallocate( own%local2global )
      deallocate( own%global2local )
      deallocate( own%private )
      deallocate( own%own_shared )
      deallocate( own%sharing_count )
      deallocate( own%sharing_procs )
      deallocate( own%MPI_sharing_type )
      deallocate( own%shared )
      deallocate( own%shared_count )
      deallocate( own%shared_owner )
      deallocate( own%MPI_shared_type )
      deallocate( own%internal_blks )
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine wmpi_print_node_own_all           *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 1/24/2015 rhd              *
c     *                                                              *
c     *     this subroutine prints the entire node owner structure   *
c     *     for domain decomposition portions of the MPI code        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_print_node_own_all( own )
      use global_data ! old common.main
      use mpi_lnpcg
      implicit none
c
      type (node_owner_type) :: own
c
      integer :: proc, i, j
c
      write(out,*) '>>> # of nodes:',own%num_local_nodes
      write(out,*) '>>> local2global:'
      write(out,'(5x,12i7)') (own%local2global(i),i=1,
     &     own%num_local_nodes*mxndof)
c
      write(out,*) '>>> global2local:'
      write(out,'(17x,"1",17x,"2",17x,"3",17x,"4")')
      write(out,'(5x,12i7)') (own%global2local(i),i=1,nodof)
c
      write(out,*) '>>> private nodes:',own%num_private
      write(out,'(5x,10i7)') (own%private(i),i=1,own%num_private)
c
      write(out,*) '>>> owned but shared nodes:',own%num_own_shared
      write(out,'(5x,10i7)') (own%own_shared(i),i=1,own%num_own_shared)
      do proc = 0, numprocs-1
         write(out,*) '   ==> nodes shared w/ proc',proc,':',
     &        own%sharing_count(proc)
         if( own%sharing_count(proc) .eq. 0 ) cycle
         write(out,'(8x,12i7)')
     &        (own%sharing_procs(proc)%ptr(i),
     &        i=1,own%sharing_count(proc))
      end do
c
      write(out,*) '>>> shared nodes:',own%num_shared
      write(out,'(5x,10i7)') (own%shared(i),i=1,own%num_shared)
      do proc = 0, numprocs-1
         write(out,*) '   ==> nodes owned by proc',proc,':',
     &        own%shared_count(proc)
         if( own%shared_count(proc) .eq. 0 ) cycle
         write(out,'(8x,12i7)')
     &        (own%shared_owner(proc)%ptr(i),i=1,own%shared_count(proc))
      end do
c
      write(out,*) '>>> list of internal blks'
      write(out,'("      internal:",20i4)')
     &     (own%internal_blks(j),j=1,own%num_int_blks)
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine wmpi_print_node_own               *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 2/24/2017 rhd              *
c     *                                                              *
c     *     prints selected information about node ownership         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_print_node_own( own, procnum )
      use global_data ! old common.main
      use mpi_lnpcg
      implicit none
c
      integer :: procnum
      type (node_owner_type) :: own
c
      integer :: sharing_procs, owning_procs, proc
c
      sharing_procs = 0
      owning_procs  = 0
c
      do proc = 0, numprocs - 1
        if( proc .eq. procnum ) cycle
        if( own%sharing_count(proc) .gt. 0 )
     &        sharing_procs = sharing_procs  + 1
        if( own%shared_count(proc) .gt. 0 )
     &        owning_procs = owning_procs + 1
      end do
c
      write(out,'(5x,"Processor ",i3,":")') procnum
      write(out,'(10x,"total # of nodes     :",i8)')own%num_local_nodes
      write(out,'(10x,"private nodes        :",i8," (",f5.1,"%)")')
     &       own%num_private, 100.0 *
     &       float(own%num_private)/float(own%num_local_nodes)
c
      write(out,'(10x,"nodes owned & shared :",i8," (",f5.1,"%)")')
     &       own%num_own_shared, 100.0 *
     &       float(own%num_own_shared)/float(own%num_local_nodes)
      write(out,'(10x,"# of sharing procs   :",i8)') sharing_procs
c
      write(out,'(10x,"nodes not owned      :",i8," (",f5.1,"%)")')
     &       own%num_shared, 100.0 *
     &       float(own%num_shared)/float(own%num_local_nodes)
      write(out,'(10x,"# of owning procs    :",i8,/)') owning_procs
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_copy_own                *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/24/2017 rhd             *
c     *                                                              *
c     *     this subroutine copies the node owner structure          *
c     *     from one structure to another, allocated to the same     *
c     *     size. Assumes that the constants have already been       *
c     *     set in the new copy.                                     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_copy_own( old, new )
      use global_data ! old common.main
      use mpi_lnpcg
      implicit none
c
      type (node_owner_type) :: old, new
c
      integer :: ierr, i, j, allo_stat
c
      do i = 1, old%num_local_nodes * mxndof
         new%local2global(i) = old%local2global(i)
      end do
      do i = 1, nodof
         new%global2local(i) = old%global2local(i)
      end do
c
      do i = 1, old%num_private
         new%private(i) = old%private(i)
      end do
c
c
      do i = 1, old%num_own_shared
         new%own_shared(i) = old%own_shared(i)
      end do
c
      do i = 0, numprocs-1
         new%sharing_count(i) = old%sharing_count(i)
         if( new%sharing_count(i) .eq. 0 ) cycle
c
         allocate( new%sharing_procs(i)%ptr(new%sharing_count(i)),
     &        stat = allo_stat )
         if(  allo_stat .ne. 0) then
            write(out,*) '>>>> allocate error in copy_own'
            call die_abort
         endif
c
         do j = 1, new%sharing_count(i)
            new%sharing_procs(i)%ptr(j) = old%sharing_procs(i)%ptr(j)
         end do
c
      end do
c
      do i = 1, old%num_shared
         new%shared(i) = old%shared(i)
      end do
c
      do i = 0, numprocs-1
         new%shared_count(i) = old%shared_count(i)
         if( new%shared_count(i) .eq. 0 ) cycle
c
         allocate( new%shared_owner(i)%ptr(new%shared_count(i)),
     &        stat = allo_stat )
         if(  allo_stat .ne. 0 ) then
            write(out,*) '>>>> allocate error in copy_own'
            call die_abort
         endif
c
         do j = 1, new%shared_count(i)
            new%shared_owner(i)%ptr(j) = old%shared_owner(i)%ptr(j)
         end do
      end do
c
      new%num_int_blks = old%num_int_blks
      do i = 1, new%num_int_blks
         new%internal_blks(i) = old%internal_blks(i)
      end do
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine wmpi_chknode                     *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/24/2017 rhd             *
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
      use global_data ! old common.main
      use mpi_lnpcg, only : local_nodes
      implicit none
      include "mpif.h"
c
      integer :: node, ptr
      logical :: referenced, owned
c
      ptr = local_nodes%global2local(dstmap(node))
c
      if( ptr .eq. 0 ) then
         referenced = .false.
         owned = .false.
      else
         referenced = .true.
         if( (ptr+2) / 3 .le. (local_nodes%num_private +
     &        local_nodes%num_own_shared) ) then
            owned = .true.
         else
            owned = .false.
         endif
      endif
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *              drive the unsymmetric CPardiso solution         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                  last modified : 1/12/20167 rhd              *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine cpardiso_unsymmetric( neq, ncoeff, k_pointers,
     &            k_indexes, eqn_coeffs, rhs, solution_vec,
     &            print_cpu_stats, itype, out, myrank )
      implicit none
      include 'mpif.h'
c
c                parameter declarations
c
      integer :: neq, ncoeff, k_pointers(neq+1), k_indexes(ncoeff),
     &           itype, out, myrank
      double precision :: eqn_coeffs(ncoeff), rhs(neq),
     &                    solution_vec(neq)
      logical :: print_cpu_stats
c
c                locals for driver
c
      logical :: use_iterative, use_direct, time_stats, master, worker,
     &           debug_worker, debug_master
      logical, save :: cpardiso_mat_defined
      logical, parameter :: debug_master_flag = .false.,
     &                      debug_worker_flag = .false.
      integer, save :: num_calls
      real, external :: wwalltime
c
c                locald for cpardiso
c
      integer(kind=8), save ::  pt(64)
      integer :: maxfct, mnum, phase, nrhs, error, idum, len_iparm,
     &           ierror
      integer, save :: iparm(64), mtype, msglvl
      double precision :: ddum
c
      data num_calls, cpardiso_mat_defined / 0, .false. /
c
c
c       Compressed Sparse Row format (CSR) for unsymmetric equations
c
c        neq           --  number of equations
c        ncoeff        --  number of non-zero terms
c        rhs           --  right-hand side
c        solution_vec  --  solution vector
c        eqn_coeffs    --  non-zero values on each row
c        k_pointers    --  dimension (n+1). k_pointers(i) points to the
c                          first column index of row i in the array
c                          k_indexes in compressed sparse row format.
c                          k_pointers(i) gives the index of the element
c                          in vector eqn_coeffs that contains the first
c                          non-zero element from row i of eqn_coeffs.
c                          The last element k_pointers(n+1) is taken
c                          to be equal to the number of non-zero
c                          elements in eqn_coeffs, plus one.
c        k_indexes     --  contains column indexes of the sparse
c                          matrix. The indices in each row must
c                          be sorted in increasing order.
c
c              Vectors describing the equations are dummies of length 1
c              on all workers
c
      master        = myrank .eq. 0
      worker        = myrank > 0
      debug_master  = master .and. debug_master_flag
      debug_worker  = worker .and. debug_worker_flag
      nrhs          = 1
      maxfct        = 1
      mnum          = 1
      len_iparm     = 64      ! CPardiso required
      use_iterative = .false. ! CPardiso does not support
      use_direct    = .true.  ! direct option
c
      call wmpi_bcast_int( itype ) ! workers need this
c
      if( master )
     &  call cpardiso_messages( 12, out, error, print_cpu_stats,
     &                          iparm, master )
c
c              solution types (itype):
c                 1 - first time solution for a matrix:
c                               setup ordering method and perform
c                               pre-processing steps.
c                 2 - Solution of above same matrix equations
c                     with a new set of coefficients but same
c                     sparsity
c                 3 - no solution. just release data.
c
      select case( itype )
c
      case( 1 )
        call cpardiso_unsymmetric_setup
        call cpardiso_unsymmetric_solve
      case( 2 )
        call cpardiso_unsymmetric_solve
      case( 3 ) ! release all CPardiso data and return
        if( worker ) return
        call cpardiso_unsymmetric_release
      case default ! then die
         call cpardiso_messages( 8, out, error, print_cpu_stats,
     &                           iparm, master )
c
      end select
      return
c
      contains
c     ========
c
c     ******************************************************************
c     *       contains:   cpardiso_symmetric_release                   *
c     ******************************************************************
c

      subroutine cpardiso_unsymmetric_release
      implicit none
c
      if( .not. cpardiso_mat_defined ) then ! job will be aborted
         call cpardiso_messages( 7, out, error, print_cpu_stats,
     &                           iparm, master )
      end if
c
      phase = -1 ! release internal memory
      call cluster_sparse_solver( pt, maxfct, mnum, mtype,
     &                phase, neq, ddum, idum, idum, idum, nrhs, iparm,
     &                msglvl, ddum, ddum, MPI_COMM_WORLD, error )
      cpardiso_mat_defined = .false.
      num_calls = 0
c
      return
      end  subroutine cpardiso_unsymmetric_release
c
c     ******************************************************************
c     *       contains:   cpardiso_unsymmetric_setup                    *
c     ******************************************************************
c
      subroutine cpardiso_unsymmetric_setup
      implicit none
c
      if( master ) then
          if( cpardiso_mat_defined ) then
             phase = -1 ! release internal memory
             call cluster_sparse_solver( pt, maxfct, mnum, mtype,
     &                phase, neq, ddum, idum, idum, idum, nrhs, iparm,
     &                msglvl, ddum, ddum, MPI_COMM_WORLD, error )
             cpardiso_mat_defined = .false.
             call cpardiso_messages( 6, out, error, print_cpu_stats,
     &                               iparm, master )
          end if
      end if
c
      call thyme( 23, 1 )
      pt(1:64) = 0 ! CPardiso required input (note I*8)
c
      call cpardiso_unsymmetric_set_iparm( iparm, len_iparm,
     &         use_iterative, error, msglvl, mtype, out  )
c
      call cpardiso_messages( 11, out, error, print_cpu_stats, iparm,
     &                        master )
c
c                  input sparsity structure, reorder, symbolic
c                  factorization. CPardiso sends data to workers.
c
      phase = 11 ! reordering and symbolic factorization
      call cpardiso_messages( 10, out, error, print_cpu_stats, iparm,
     &                        master )
      call cluster_sparse_solver( pt, maxfct, mnum, mtype, phase, neq,
     &              eqn_coeffs, k_pointers, k_indexes, idum, nrhs,
     &              iparm, msglvl, ddum, ddum, MPI_COMM_WORLD, error )
      call cpardiso_messages( 2, out, error, print_cpu_stats, iparm,
     &                        master)
      call thyme( 23, 2 )
c
      return
      end  subroutine cpardiso_unsymmetric_setup

c
c     ******************************************************************
c     *       contains:   cpardiso_unsymmetric_solve                   *
c     ******************************************************************
c

      subroutine cpardiso_unsymmetric_solve
      implicit none
c
c              barrier may not be needed
c
      num_calls = num_calls + 1
      call thyme( 25, 1)
      phase = 22 ! only factorization
      call cpardiso_messages( 4, out, error, print_cpu_stats, iparm,
     &                        master )
      call MPI_BARRIER( MPI_COMM_WORLD, ierror )
      call cluster_sparse_solver( pt, maxfct, mnum, mtype, phase, neq,
     &                eqn_coeffs, k_pointers, k_indexes, idum, nrhs,
     &                iparm, msglvl, ddum, ddum, MPI_COMM_WORLD, error )
      call cpardiso_messages( 3, out, error, print_cpu_stats, iparm,
     &                        master )
      call thyme( 25, 2 )
c
      call thyme( 26, 1)
      phase = 33   ! forward/backward solve
      call cluster_sparse_solver( pt, maxfct, mnum, mtype, phase, neq,
     &                eqn_coeffs, k_pointers, k_indexes, idum, nrhs,
     &                iparm, msglvl, rhs, solution_vec, MPI_COMM_WORLD,
     &                error )
      call cpardiso_messages( 5, out, error, print_cpu_stats, iparm,
     &                        master )
      call thyme( 26, 2 )
c
      return
      end  subroutine cpardiso_unsymmetric_solve
      end  subroutine cpardiso_unsymmetric


c ********************************************************************
c *                                                                  *
c *            setup iparm() for CPardiso unsymmetric                *
c *                                                                  *
c *                  last updated:  1/12/2016 rhd                    *
c *                                                                  *
c ********************************************************************
c
      subroutine cpardiso_unsymmetric_set_iparm( iparm, len_iparm,
     &    use_iterative, error, msglvl, mtype, iout  )
      implicit none
c
      integer :: len_iparm, iparm(len_iparm), error, msglvl, mtype,
     &           iout
      logical :: use_iterative
c
      error  = 0
      msglvl = 0    !  request no messages from CPardiso
      mtype  = 1 ! structurally symmetric, values not symmetric
c
      iparm(1) = 1 ! we set all iparm values here
      iparm(2) = 3 ! threaded, non-MPI reordering
      iparm(3) = 0 ! CPardiso reserved
      iparm(4) = 0 ! not used by CPardiso on input
      iparm(5) = 0 ! not used by CPardiso on input
      iparm(6) = 0 ! rhs NOT overwritten by solution
      iparm(7) = 0 ! output - actual # iterative refinement steps
      iparm(8) = 30! lots of iterative refinement steps
      iparm(9) = 0 ! not used by CPardiso on input
      iparm(10) = 13 ! perturb the pivot elements with 1E-13
                     ! default value is 1E-8
      iparm(11) = 1 ! scaling.
      iparm(12) = 0 ! not used by CPardiso on input
      iparm(13) = 1 ! Add some extra permutation for
c                     non-symmetric matrices
      iparm(14) = 0 ! not used by CPardiso on input
      iparm(15) = 0 ! not used by CPardiso on input
      iparm(16) = 0 ! not used by CPardiso on input
      iparm(17) = 0 ! not used by CPardiso on input
      iparm(18) = -1 ! return: number of nonzeros in the factor LU
      iparm(19) = -1 ! return: Mflops for LU factorization
      iparm(20) = -1 ! return: Numbers of CG Iterations
      iparm(21) = 0 ! 1x1 pivoting. try = 1 if issues
      iparm(22) = 0 ! not used by CPardiso on input
      iparm(23) = 0 ! not used by CPardiso on input
      iparm(24) = 0 ! not used by CPardiso on input
      iparm(25) = 0 ! not used by CPardiso on input
      iparm(26) = 0 ! not used by CPardiso on input
      iparm(27) = 1 ! CPardiso checks matrix for indexing errors
c                     set = 0 for production
      iparm(28) = 0 ! input arrays are double precision
      iparm(29) = 0 ! not used by CPardiso on input
      iparm(30) = 0 ! not used by CPardiso on input
      iparm(31) = 0 ! not used by CPardiso on input
      iparm(32) = 0 ! not used by CPardiso on input
      iparm(33) = 0 ! not used by CPardiso on input
      iparm(34) = 0 ! not used by CPardiso on input
      iparm(35) = 0 ! input arryas are 1-based (Fortran)
      iparm(36) = 0 ! not used by CPardiso on input
      iparm(37) = 0 ! CSR arrya format. we convert to CSR from VSS
      iparm(38) = 0 ! CPardiso reserved
      iparm(39) = 0 ! CPardiso reserved
      iparm(40) = 0 ! only rank 0 get A matrix. workers do not need A.
c                     solution vector only on rank 0
      iparm(41) = 0 ! not used for iparm(40) = 0.
      iparm(42) = 0 ! not used for iparm(40) = 0.
c
      iparm(43:64) = 0 ! CPardiso reserved
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *              drive the symmetric CPardiso solution           *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/7/20167 rhd              *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine cpardiso_symmetric( neq, ncoeff, k_diag, rhs,
     &   solution_vec, eqn_coeffs, k_pointers, k_indices,
     &   print_cpu_stats, itype, out, myrank )
      implicit none
      include 'mpif.h'
c
c              parameters
c
      integer :: neq, ncoeff, k_pointers(*), k_indices(*), out,
     &           myrank, itype
      logical :: print_cpu_stats
      double precision :: solution_vec(*), eqn_coeffs(*), k_diag(*),
     &                    rhs(*)
c
c              locals for driver
c
      logical :: use_iterative, use_direct, time_stats, master, worker,
     &           debug_worker, debug_master
      logical, save :: cpardiso_mat_defined
      logical, parameter :: debug_master_flag = .false.,
     &                      debug_worker_flag = .false.
      integer :: i, nterms
      integer, save :: num_calls
      real, external :: wwalltime
c
c              locals for CPardiso
c
      integer(kind=8), save ::  pt(64)
      integer :: maxfct, mnum, phase, nrhs, error,
     &           idum, len_iparm, ierror
      integer, save :: iparm(64), mtype, msglvl
      double precision :: ddum
c
      data num_calls, cpardiso_mat_defined / 0, .false. /
c
c
c       Input description of sparse equations in GPS/VSS format
c
c        neq           --  number of equations
c        ncoeff        --  number of off-diagonal coefficients
c        k_diag        --  diagonal values
c        rhs           --  right-hand side
c        solution_vec  --  solution vector
c        eqn_coeffs    --  off-diagonal values
c        k_pointers    --  number of non-zero values on each row
c                          excluding the diagonal
c       k_indices      --  column numbers off diagonal terms
c
c
c       Compressed Sparse Row format (CSR) for symmetric equations
c       after mapping. These vectors describe equations for CPardiso.
c
c        neq           --  number of equations
c        ncoeff        --  number of non-zero terms in upper-triangle
c                          A (including k_diag)
c        rhs           --  right-hand side
c        solution_vec  --  solution vector
c        eqn_coeffs    --  non-zero values on each row starting with
c                          diagonal (was off_off_diagonals)
c        k_pointers    --  see diagram example in drive_pardiso.f
c        k_indices     --  contains column indices of the sparse
c                          matrix. The indices in each row must
c                          be sorted in increasing order.
c
c              Vectors describing the equations are dummies of length 1
c              on all workers
c
      master        = myrank .eq. 0
      worker        = myrank > 0
      debug_master  = master .and. debug_master_flag
      debug_worker  = worker .and. debug_worker_flag
      nrhs          = 1
      maxfct        = 1
      mnum          = 1
      len_iparm     = 64      ! CPardiso required
      use_iterative = .false. ! CPardiso does not support
      use_direct    = .true.  ! direct option
c
      call wmpi_bcast_int( itype ) ! workers need this
c
      if( master )
     &  call cpardiso_messages( 1, out, error, print_cpu_stats,
     &                          iparm, master )
c
c              solution types (itype):
c                 1 - first time solution for a matrix:
c                               setup ordering method and perform
c                               pre-processing steps.
c                 2 - Solution of above same matrix equations
c                     with a new set of coefficients but same
c                     sparsity
c                 3 - no solution. just release data.
c
      select case( itype )
c
      case( 1 )
        call cpardiso_symmetric_setup
        call cpardiso_symmetric_solve
      case( 2 )
        if( master ) call pardiso_symmetric_map( neq, ncoeff, k_diag,
     &                       rhs, eqn_coeffs, k_pointers, k_indices )
        call cpardiso_symmetric_solve
      case( 3 ) ! release all CPardiso data and return
        if( worker ) return
        call cpardiso_symmetric_release
      case default ! then die
         call cpardiso_messages( 8, out, error, print_cpu_stats,
     &                           iparm, master )
c
      end select
      return
c
      contains
c     ========
c
c     ******************************************************************
c     *       contains:   cpardiso_symmetric_release                   *
c     ******************************************************************
c

      subroutine cpardiso_symmetric_release
      implicit none
c
      if( .not. cpardiso_mat_defined ) then ! job will be aborted
         call cpardiso_messages( 7, out, error, print_cpu_stats,
     &                           iparm, master )
      end if
c
      phase = -1 ! release internal memory
      call cluster_sparse_solver( pt, maxfct, mnum, mtype,
     &                phase, neq, ddum, idum, idum, idum, nrhs, iparm,
     &                msglvl, ddum, ddum, MPI_COMM_WORLD, error )
      cpardiso_mat_defined = .false.
      num_calls = 0
c
      return
      end  subroutine cpardiso_symmetric_release
c
c     ******************************************************************
c     *       contains:   cpardiso_symmetric_setup                     *
c     ******************************************************************
c
      subroutine cpardiso_symmetric_setup
      implicit none
c
      if( master ) then
          call pardiso_symmetric_map( neq, ncoeff, k_diag,
     &             rhs, eqn_coeffs, k_pointers, k_indices )
          if( cpardiso_mat_defined ) then
             phase = -1 ! release internal memory
             call cluster_sparse_solver( pt, maxfct, mnum, mtype,
     &                phase, neq, ddum, idum, idum, idum, nrhs, iparm,
     &                msglvl,ddum, ddum, MPI_COMM_WORLD, error )
             cpardiso_mat_defined = .false.
             call cpardiso_messages( 6, out, error, print_cpu_stats,
     &                               iparm, master )
          end if
      end if
c
      call thyme( 23, 1 )
      pt(1:64) = 0 ! CPardiso required input (note I*8)
c
      call cpardiso_symmetric_set_iparm( iparm, len_iparm,
     &         use_iterative, error, msglvl, mtype, out  )
c
      call cpardiso_messages( 11, out, error, print_cpu_stats, iparm,
     &                        master )
c
c                  input sparsity structure, reorder, symbolic
c                  factorization. CPardiso sends data to workers.
c
      phase = 11 ! reordering and symbolic factorization
      call cpardiso_messages( 10, out, error, print_cpu_stats, iparm,
     &                        master )
      call cluster_sparse_solver( pt, maxfct, mnum, mtype, phase, neq,
     &              eqn_coeffs, k_pointers, k_indices, idum, nrhs,
     &              iparm, msglvl, ddum, ddum, MPI_COMM_WORLD, error )
      call cpardiso_messages( 2, out, error, print_cpu_stats, iparm,
     &                        master)
      call thyme( 23, 2 )
c
      return
      end  subroutine cpardiso_symmetric_setup

c
c     ******************************************************************
c     *       contains:   cpardiso_symmetric_solve                     *
c     ******************************************************************
c

      subroutine cpardiso_symmetric_solve
      implicit none
c
      num_calls = num_calls + 1
      call thyme( 25, 1)
      phase = 22 ! only factorization
      call cpardiso_messages( 4, out, error, print_cpu_stats, iparm,
     &                        master )
      call MPI_BARRIER( MPI_COMM_WORLD, ierror )
      call cluster_sparse_solver( pt, maxfct, mnum, mtype, phase, neq,
     &                eqn_coeffs, k_pointers, k_indices, idum, nrhs,
     &                iparm, msglvl, ddum, ddum, MPI_COMM_WORLD, error )
      call cpardiso_messages( 3, out, error, print_cpu_stats, iparm,
     &                        master )
      call thyme( 25, 2 )
c
      call thyme( 26, 1)
      phase = 33   ! forward/backward solve
      call cluster_sparse_solver( pt, maxfct, mnum, mtype, phase, neq,
     &                eqn_coeffs, k_pointers, k_indices, idum, nrhs,
     &                iparm, msglvl, rhs, solution_vec, MPI_COMM_WORLD,
     &                error )
      call cpardiso_messages( 5, out, error, print_cpu_stats, iparm,
     &                        master )
      call thyme( 26, 2 )
c
      return
      end  subroutine cpardiso_symmetric_solve
      end  subroutine cpardiso_symmetric


c ********************************************************************
c *                                                                  *
c *                setup iparm() for CPardiso                        *
c *                                                                  *
c *                  last updated:  1/5/2016 rhd d                   *
c *                                                                  *
c ********************************************************************
c
      subroutine cpardiso_symmetric_set_iparm( iparm, len_iparm,
     &    use_iterative, error, msglvl, mtype, iout  )
      implicit none
c
      integer :: len_iparm, iparm(len_iparm), error, msglvl, mtype,
     &           iout
      logical :: use_iterative
c
      error  = 0
      msglvl = 0    !  request no messages from CPardiso
      mtype  = -2   !  real, symmetric, indefinite
c
      iparm(1) = 1 ! we set all iparm values here
      iparm(2) = 3 ! threaded, non-MPI reordering
      iparm(3) = 0 ! CPardiso reserved
      iparm(4) = 0 ! not used by CPardiso on input
      iparm(5) = 0 ! not used by CPardiso on input
      iparm(6) = 0 ! rhs NOT overwritten by solution
      iparm(7) = 0 ! output - actual # iterative refinement steps
      iparm(8) = 0 ! use default # iterative refinement steps
      iparm(9) = 0 ! not used by CPardiso on input
      iparm(10) = 13 ! perturb the pivot elements with 1E-13
                     ! default value is 1E-8
      iparm(11) = 0 ! scaling.none for real, symm, indef
      iparm(12) = 0 ! not used by CPardiso on input
      iparm(13) = 0 ! no symmetric weighting. default for symmetric.
c                     Try iparm(11,13) = 1 in case of inappropriate
c                     accuracy
      iparm(14) = 0 ! not used by CPardiso on input
      iparm(15) = 0 ! not used by CPardiso on input
      iparm(16) = 0 ! not used by CPardiso on input
      iparm(17) = 0 ! not used by CPardiso on input
      iparm(18) = -1 ! output number of non-zero terms in factor
      iparm(19) = -1 ! output MFLOP for factorization & solution
      iparm(20) = 0 ! not used by CPardiso on input
      iparm(21) = 0 ! 1x1 pivoting. try = 1 if issues
      iparm(22) = 0 ! not used by CPardiso on input
      iparm(23) = 0 ! not used by CPardiso on input
      iparm(24) = 0 ! not used by CPardiso on input
      iparm(25) = 0 ! not used by CPardiso on input
      iparm(26) = 0 ! not used by CPardiso on input
      iparm(27) = 1 ! CPardiso checks matrix for indexing errors
c                     set = 0 for production
      iparm(28) = 0 ! input arrays are double precision
      iparm(29) = 0 ! not used by CPardiso on input
      iparm(30) = 0 ! not used by CPardiso on input
      iparm(31) = 0 ! not used by CPardiso on input
      iparm(32) = 0 ! not used by CPardiso on input
      iparm(33) = 0 ! not used by CPardiso on input
      iparm(34) = 0 ! not used by CPardiso on input
      iparm(35) = 0 ! input arryas are 1-based (Fortran)
      iparm(36) = 0 ! not used by CPardiso on input
      iparm(37) = 0 ! CSR arrya format. we convert to CSR from VSS
      iparm(38) = 0 ! CPardiso reserved
      iparm(39) = 0 ! CPardiso reserved
      iparm(40) = 0 ! only rank 0 get A matrix. workers do not need A.
c                     solution vector only on rank 0
      iparm(41) = 0 ! not used for iparm(40) = 0.
      iparm(42) = 0 ! not used for iparm(40) = 0.
c
      iparm(43:64) = 0 ! CPardiso reserved
c
      return
      end
c
c ********************************************************************
c *                                                                  *
c *                message output from driver program                *
c *                                                                  *
c *                  last updated:  1/7/2017 rhd                     *
c *                                                                  *
c ********************************************************************
c
      subroutine cpardiso_messages( mess_no, iout, ier,
     &                              print_time_stats, iparm, master )
      implicit none
c
c             parameters
c
      integer :: mess_no, iout, ier
      logical :: print_time_stats, master
      integer, dimension (*) :: iparm
c
c             locals
c
      real, save :: start_factor_wtime
      real, external :: wwalltime
c
      if( .not. master ) return
c
      select case ( mess_no )
c
        case( 1 )  ! used
         if( print_time_stats ) write(iout,9480) wwalltime(1)
c
        case( 2 ) ! used
         if( ier .ne. 0 ) then
           write(iout,9470)
           write(iout,*) '       >> @ phase 11, ier: ',ier
           stop
         end if
         if( print_time_stats )  then
           write(iout,9482) wwalltime(1)
           write(iout,2010) dble(iparm(18))/1000.d0/1000.d0/1000.d0
           write(iout,2020) dble(iparm(19))/1000.d0
           write(iout,2022) dble(iparm(15))/1000.d0/1000.d0
         end if
c
        case( 3 )  ! used
         if( ier .ne. 0 ) then
            write(iout,9471) ier
            call die_abort
         end if
         if( print_time_stats ) then
           write(iout,9490) wwalltime(1)
           write(iout,2030) dble(iparm(17))/1000.d0/1000.d0
         end if

        case( 4 ) ! used
         if( print_time_stats ) then
           write(iout,9184)
           if( master ) start_factor_wtime = wwalltime(1)
         end if
c
        case( 5 )   ! used
         if( ier .ne. 0 ) then
           write(iout,9470)
           write(iout,*) '       >> @ phase 33, ier: ',ier
           call die_abort
         end if
         if( print_time_stats ) then
           write(iout,9492) wwalltime(1)
           write(iout,9494) wwalltime(1) - start_factor_wtime
         end if
c
        case( 6 )  !  used
         if( ier .ne. 0 ) then
           write(iout,9470)
           write(iout,*) '       >>  phase -1, ier: ',ier
           write(iout,*) '       >>  releasing solver data'
           call die_abort
         end if
c
        case( 7 )    ! used
           write(iout,9470)
           call die_abort
c
        case( 8 )
           write(iout,9485)
           call die_abort
c
        case( 9 )
c
        case( 10 ) !used
          if( print_time_stats ) write(iout,9300) wwalltime(1)
c
        case( 11 )  !  used
          if( print_time_stats ) write(iout,9305) wwalltime(1)
c
        case( 12 )  !  used
          if( print_time_stats ) write(iout,9500) wwalltime(1)
c
      end select
      return
c
 2010 format(
     &  15x,'no. terms in [L][U] factor:      ', f9.3,' x 10**9' )
 2020 format(
     &  15x,'factor + sol op count (GFlop):   ', f9.3)
 2022 format(
     &  15x,'reordering memory (GB):          ', f9.3)
 2030 format(
     &  15x,'factorization memory (GB):       ', f9.3)
 9184  format(
     &  15x, 'start cluster factorization    ')
 9300  format(
     &  15x, 'start reorder, symb. factor.  @ ',f10.2 )
 9305  format(
     &  15x, 'done mapping vss -> csr       @ ',f10.2 )
 9470  format(
     &   15x, 'FATAL ERRROR: WARP3D has requested deletion of all',
     & /,15x, '              CPardiso data but no equations have been',
     & /,15x  '              solved. Job aborted.')
 9471  format(
     &   15x, 'FATAL ERRROR: CPardiso failed during factorization with',
     & /,15x, '              error code: ',i5,
     & /,15x  '              Job aborted.')
 9480  format(
     &  15x, 'start cluster symmetric solve @ ',f10.2 )
 9482  format(
     &  15x, 'reorder-symbolic factor. done @ ',f10.2 )
 9485  format(
     &   15x, 'FATAL ERRROR: WARP3D has requested an invalid solution',
     & /,15x, '              type in cpardiso_symmetric',
     & /,15x  '              Job aborted.')
 9490  format(
     &  15x, 'numeric factorization done    @ ',f10.2 )
 9492  format(
     &  15x, 'numeric loadpass done         @ ',f10.2 )
 9494  format(
     &  15x, 'factor + solve wall time:       ',f10.2 )
 9500  format(
     &  15x, 'start cluster unsymmet  solve @ ',f10.2 )
c
      end
c
