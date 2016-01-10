c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_init                    *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : rhd 1/9/2016 rhd           *
c     *                                                              *
c     *     This routine initializes MPI on all processors at the    *
c     *     beginning of a MPI warp run.  It also creates the        *
c     *     MPI_VAL datatype which is equivalent to either a real or *
c     *     double precision variable depending on whether the       *
c     *     platform is double precision or single precision. This   *
c     *     greatly simplifies later MPI calls without requiring     *
c     *     excessive #sgl,#dbl pairs.                               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_init
      use hypre_parameters
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      logical ldum1,ldum2
      character *1 dums
      real dumr
#dbl      double precision
#sgl      real
     &     dumd
c
c             initialize MPI, find the total number of processors involved
c             with this run, and get the id number of the local processor in
c             the total processor numbering. Processor 0 will be the root
c             processor, while all others will be slaves.
c
      call MPI_INIT_THREAD (MPI_THREAD_FUNNELED, iallowed, ierr)
      if( iallowed .ne. MPI_THREAD_FUNNELED ) then
        write(*,*) '>>> fatal error on MPI start up.'
        write(*,*) '    MPI_THREAD_FUNNELED, iallowed:',
     &                  MPI_THREAD_FUNNELED, iallowed
        call die_abort
      end if
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      write (*,'(" >>> numprocs:",i4," myid:", i4)') numprocs, myid
      if (numprocs .gt. max_procs) then
         call errmsg (312, dum, dums, dumr, dumd)
         call MPI_FINALIZE (ierr)
         call die_abort
         stop
      endif
c
c             create MPI_VAL datatype to be equal to either MPI_REAL
c             or MPI_DOUBLE_PRECISION; simplifies the #dbl #sgl specification.
c             MPI_VAL is stored in common.main, while MPI_REAL and
c             MPI_DOUBLE_PRECISION is in mpif.h.
c
c             MPI_TYPE_CONTIGUOUS defines MPI_VAL to be equivalent to
c             either MPI_REAL or MPI_DOUBLE_PRECISION, while MPI_TYPE_COMMIT
c             officially registers the datatype with the MPI routines.
c
#sgl      call MPI_TYPE_CONTIGUOUS(1,MPI_REAL,MPI_VAL,ierr)
#dbl      call MPI_TYPE_CONTIGUOUS(1,MPI_DOUBLE_PRECISION,MPI_VAL,ierr)
      call MPI_TYPE_COMMIT (MPI_VAL,ierr)
c
c             set the logical flags root_processor and slave_processor
c             for simple identification of whether we are the slave
c             or the master.
c
      if (myid .eq. 0) then
         root_processor = .true.
         slave_processor = .false.
      else
         root_processor = .false.
         slave_processor = .true.
      endif
c
c             set use_mpi flag to inidcate that this is the MPI version
c
      use_mpi = .true.
c
c	      have all the slave provessors find out their process IDs
c	      and then send these numbers back to the root process.  These
c	      will be used to suspend the mpi procs if a threaded
c	      sparse solver is called.
c
#l64      call wmpi_procids
c
c             if we are the root processor, return to driver, and
c             continue executing the code as normal, reading input, etc.
c
      if (root_processor) return
c
c             processors executing code past this point are slave
c             processors. Call initst to initialize variables they
c             may need, then throw them in the slave handler.  In this
c             routine they wait for instructions from the root processor
c             about what job to do. The slave processors never return
c             from the slave handler; they die when told to do so by
c             the root processor.
c
      call initst(ldum1,ldum2)
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
c     *                   last modified : 1/24/2015 rhd              *
c     *                                                              *
c     *     This routine has all the slave processors find their     *
c     *     process ids, then send these Ids back to the root.       *
c     *     This allows root to suspend the MPI threads if a         *
c     *     threaded sparse solver is called.                        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_procids
      use mpi_lnpcg
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      logical debug
      data debug /.false./
      integer :: status(MPI_STATUS_SIZE)
c
c          if this is the the root processor, then recieve data from each
c          processor.
c
      if ( root_processor ) then
c
        do proc = 1, numprocs -1
        call MPI_RECV ( proc_pids(proc), 1, MPI_INTEGER,
     &          proc, 1, MPI_COMM_WORLD, status, ierr)
        enddo
c
      if ( debug ) then
        write (out,*) '>>> Here are proc ids gathered on root:'
        do proc = 1, numprocs -1
          write (out,*) '    proc:', proc, ' id:',proc_pids(proc)
        end do
        procid = getpid ()
        if ( ierr .ne. 0) then
          write (out,*) '>>> FATAL ERROR in getting process id.'
          call die_abort
        endif
      write (out,*) '>> the pid of root is ',procid
      endif
c
c
      else
c
c           if this is a slave processor, find proc id and send
c           to root
c
      ierr = 0
      procid = getpid ()
      if ( ierr .ne. 0) then
        write (out,*) '>>> FATAL ERROR in retrieving process id.'
        call die_abort
      endif
c
      call MPI_SEND (procid, 1, MPI_INTEGER, 0, 1,
     &        MPI_COMM_WORLD, ierr)
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
c     *                   last modified : 1/24/2015 rhd              *
c     *                                                              *
c     *     This routine allows the root to suspend the slave        *
c     *     processes when a threaded sparse solver is called.       *
c     *     Once the solver completes, this routine can also         *
c     *     restart the mpi processes.                               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_suspend (option)
      use mpi_lnpcg
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      logical debug, suspended
      save suspended
      data debug, suspended /.false., .false. /
c
c	  if option = 1, suspend the slave mpi processes using system 'kill'
c         commands and hardware-specific signal numbers.
c
      if ( option .eq. 1 ) then
c
         do proc = 1, numprocs - 1
          if( debug ) write (out,*) '>>> root suspend ',proc_pids(proc)
#l64            ierr = kill (proc_pids(proc), 19)
          if ( ierr .ne. 0) then
	       write (out,*) '>>> FATAL ERROR: root unable to suspend',
     &                ' pid ',proc_pids(proc)
               call die_abort
            endif
         end do
         suspended = .true.
c
      else if ( option .eq. 2 .and. suspended ) then
c
c	  if option = 2, awaken the suspended slave mpi processes
c
         do proc = 1, numprocs - 1
          if( debug )write (out,*) '>>> root reviving ',proc_pids(proc)
#l64            ierr = kill (proc_pids(proc), 18)
          if ( ierr .ne. 0) then
	       write (out,*) '>>> FATAL ERROR: root unable to awaken',
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
c     *                   last modified : 03/05/98                   *
c     *                                                              *
c     *     The MPI implementation of warp3d follows a master-slave  *
c     *     approach for most of the code.  The root processor       *
c     *     (processor 0) runs warp3d in full, notifying the slave   *
c     *     processors when they have work to do.                    *
c     *                                                              *
c     *     This routine allows the root processor to contact the    *
c     *     slave processes and tell them what to do next.           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_alert_slaves ( do_it )
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      logical debug
      data debug / .false. /
c
c              for all slave processors, return; slaves are not
c              allowed to tell other slaves what to do.
c
      if ( slave_processor ) return
c
c              do_it is the code used by handle_mpi_slaves to
c              tell the slaves what to do.
c
      if(debug)write (out,'("====> Root is sending alert ",i3)')
     &     do_it
      call MPI_BCAST(do_it,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
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
c     *                   last modified : 03/05/98                   *
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
      subroutine wmpi_reduce_vec ( vec, size )
      implicit integer (a-z)
      include "mpif.h"
$add common.main
c
#dbl      double precision
#sgl      real
     &     vec(*)
c
c
      if ( size .eq. nodof ) then
         call wmpi_reduce_vec_new ( vec )
      else
         call wmpi_reduce_vec_std ( vec, size )
      endif
c
c      call wmpi_reduce_vec_std ( vec, size )
c      call wmpi_reduce_vec_new ( vec )
c      call wmpi_reduce_vec_log ( vec, size )
c
      return
c
      end
c
c
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_reduce_vec_std          *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/05/98                   *
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
      subroutine wmpi_reduce_vec_std ( vec, size )
      implicit integer (a-z)
      include "mpif.h"
$add common.main
c
      allocatable ::  vec_tmp(:)
c
#dbl      double precision
#sgl      real
     &     vec(*), vec_tmp, zero
      data zero /0.0/
c
c              we cannot reduce into the same vector we are sending, so
c              allocate a temporary array to hold the data
c
      allocate ( vec_tmp(size), stat = alloc_stat )
      if ( alloc_stat .ne. 0) then
         write (out,9000)
         call die_abort
      endif
c
c              copy processor local vector into temporary vector
c
      do i = 1, size
         vec_tmp(i) = vec(i)
      enddo
c
c              zero out processor local vector
c
      do i = 1, size
         vec(i) = zero
      enddo
c
c              reduce the vector using the MPI_REDUCE command.  This
c              reduces the result only to the root processor.
c
      call MPI_REDUCE (vec_tmp, vec, size,
#dbl     &     MPI_DOUBLE_PRECISION,
#sgl     &     MPI_REAL,
     &    MPI_SUM,0,MPI_COMM_WORLD, ierr)
c
c              deallocate temporary vector
c
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
c     *                      subroutine wmpi_reduce_vec_log          *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/05/98                   *
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
      subroutine wmpi_reduce_vec_log ( vec, size )
      implicit integer (a-z)
      include "mpif.h"
$add common.main
c
      allocatable ::  vec_tmp(:)
c
#dbl      double precision
#sgl      real
     &     vec(*), vec_tmp, zero
      data zero /0.0/
c
c              we cannot reduce into the same vector we are sending, so
c              allocate a temporary array to hold the data
c
      allocate ( vec_tmp(size), stat = alloc_stat )
      if ( alloc_stat .ne. 0) then
         write (out,9000)
         call die_abort
      endif
c
c              copy processor local vector into temporary vector
c
      do i = 1, size
         vec_tmp(i) = vec(i)
      enddo
c
c              zero out processor local vector
c
      do i = 1, size
         vec(i) = zero
      enddo
c
c              reduce the vector using the MPI_REDUCE command.  This
c              reduces the result only to the root processor.
c
      call MPI_REDUCE (vec_tmp, vec, size,
#dbl     &     MPI_DOUBLE_PRECISION,
#sgl     &     MPI_REAL,
     &    MPI_SUM,0,MPI_COMM_WORLD, ierr)
c
c              deallocate temporary vector
c
      deallocate ( vec_tmp )
c
      return
 9000 format ('>>>  FATAL ERROR:  could not allocate additional',/,
     &        '>>>       memory in wmpi_reduce_vec.')

      end
c
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_reduce_vec_new          *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/05/98                   *
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
      subroutine wmpi_reduce_vec_new ( vec )
      use mpi_lnpcg
      implicit integer (a-z)
      include "mpif.h"
$add common.main
c
c
#dbl      double precision
#sgl      real
     &     vec(*), vec_tmp(mxdof), zero
      dimension status (MPI_STATUS_SIZE)
      data zero /0.0/
c
c          if this is the the root processor, then recieve data from each
c          processor.
c
      if ( root_processor ) then
        do proc = 1, numprocs -1
         call MPI_RECV ( vec_tmp, procdof2glob(proc)%num_dof,
     &          MPI_VAL, proc, 1, MPI_COMM_WORLD, status, ierr)
            do i = 1, procdof2glob(proc)%num_dof
               vec(procdof2glob(proc)%dof(i)) =
     &             vec(procdof2glob(proc)%dof(i)) + vec_tmp(i)
            enddo
        enddo
      else
c
c           if this is a slave processor, send local info to root
c
      num_dof = local_nodes%num_local_nodes * mxndof
         do i = 1, num_dof
            vec_tmp(i) = vec(local_nodes%local2global(i))
         enddo
         call MPI_SEND (vec_tmp, num_dof, MPI_VAL, 0, 1,
     &        MPI_COMM_WORLD, ierr)
c
      endif
c
      return
 9000 format ('>>>  FATAL ERROR:  could not allocate additional',/,
     &        '>>>       memory in wmpi_reduce_vec.')

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
c     *                   last modified : 03/05/98                   *
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
      subroutine wmpi_red_intvec ( vec, size )
      implicit integer (a-z)
      include "mpif.h"
$add common.main
c
      integer, allocatable ::  vec_tmp(:)
c
      dimension vec(*)
c
c              we cannot reduce into the same vector we are sending, so
c              allocate a temporary array to hold the data
c
      allocate ( vec_tmp(size), stat = alloc_stat )
      if ( alloc_stat .ne. 0) then
         write (out,9000)
         call die_abort
      endif
c
c              copy processor local vector into temporary vector
c
      do i = 1, size
         vec_tmp(i) = vec(i)
      enddo
c
c              zero out processor local vector
c
      do i = 1, size
         vec(i) = 0
      enddo
c
c              reduce the vector using the MPI_REDUCE command.  This
c              reduces the result only to the root processor.
c
      call MPI_REDUCE (vec_tmp, vec, size,MPI_INTEGER,
     &     MPI_SUM,0,MPI_COMM_WORLD, ierr)
c
c              deallocate temporary vector
c
      deallocate ( vec_tmp )
c
      return
 9000 format ('>>>  FATAL ERROR:  could not allocate additional',/,
     &        '>>>       memory in wmpi_red_intvec.')

      end
c
c
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
      implicit integer (a-z)
      include "mpif.h"
$add common.main
c
      call MPI_BARRIER (MPI_COMM_WORLD, ierr)
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
c     *                   last modified : 03/05/98                   *
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
      subroutine wmpi_dotprod ( scalar )
      implicit integer (a-z)
      include "mpif.h"
c
#dbl      double precision
#sgl      real
     &     scalar, scalar_tmp
c
c
      scalar_tmp = scalar
c
      call MPI_ALLREDUCE (scalar_tmp, scalar, 1,
#dbl     &     MPI_DOUBLE_PRECISION,
#sgl     &     MPI_REAL,
     &    MPI_SUM, MPI_COMM_WORLD, ierr)
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_redint                  *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/05/98                   *
c     *                                                              *
c     *     This routine drives the reduction via mpi of a           *
c     *     integer which has a contribution on each processor.      *
c     *     The total is only available on the root processor.       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_redint ( scalar )
      implicit integer (a-z)
      include "mpif.h"
c
      scalar_tmp = scalar
c
      call MPI_REDUCE (scalar_tmp, scalar, 1, MPI_INTEGER,
     &    MPI_SUM, 0, MPI_COMM_WORLD, ierr)
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
c     *                   last modified : 03/05/98                   *
c     *                                                              *
c     *     This routine drives the reduction via mpi of a           *
c     *     logical which has a contribution on each processor.      *
c     *     The result of an OR on all the logicals is only         *
c     *     available on the root processor.                         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_redlog ( log_var )
      implicit integer (a-z)
      include "mpif.h"
c
      logical log_var, log_var_tmp
c
      log_var_tmp = log_var
c
      call MPI_ALLREDUCE (log_var_tmp, log_var, 1, MPI_LOGICAL,
     &    MPI_LOR, MPI_COMM_WORLD, ierr)
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
c     *                   last modified : 03/05/98                   *
c     *                                                              *
c     *      Broadcast an integer from root to all the slave         *
c     *      processors.                                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_bcast_int ( int_var )
      implicit integer (a-z)
      include "mpif.h"
c
      call MPI_BCAST(int_var,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
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
c     *                   last modified : 03/05/98                   *
c     *                                                              *
c     *      Broadcast a real from root to all the slave             *
c     *      processors.                                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_bcast_real ( real_var )
      implicit integer (a-z)
      real real_var
      include "mpif.h"
c
      call MPI_BCAST(real_var,1,MPI_REAL,0,MPI_COMM_WORLD, ierr)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_real               *
c     *                                                              *
c     *                       written by : mcw                       *
c     *                                                              *
c     *                   last modified : 02/17/05                   *
c     *                                                              *
c     *      Broadcast a real from slaves to root. root keeps the    *
c     *      non-zero value(s)                                       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_real ( real_vec, size )
      implicit integer (a-z)
      include "mpif.h"
$add common.main
c
      real real_vec(*), zero
      real, allocatable :: temp_vec(:,:)
      logical debug
      data debug, zero / .false., 0.0 /
c
      if( debug ) write(out,*) ' myid = ', myid, ' size = ', size,
     &                         ' real_vec() = ', (real_vec(i),i=1,size)
c
      if( root_processor ) then
c
         allocate ( temp_vec(size,numprocs), stat = alloc_stat )
         if ( alloc_stat .ne. 0) then
            write (out,9000)
            call die_abort
         endif
c
         do proc = 1, numprocs - 1
c
            call MPI_RECV(real_vec, size, MPI_REAL, proc, 1,
     &                    MPI_COMM_WORLD, status, ierr)
c
            do i = 1, size
               temp_vec(i,proc) = real_vec(i)
               if( debug ) then
                  write(out,*) 'real_vec(', i,') = ', real_vec(i)
                  write(out,*) 'temp_vec(', i,',',proc,') = ',
     &                          temp_vec(i,proc)
               end if
            end do
c
         end do
c
c             assign non-zero values of recieved data to real_vec
c             on root.
c
         do proc = 1, numprocs - 1
            do i = 1, size
               if( temp_vec(i,proc) .ne. zero ) then
                  real_vec(i) = temp_vec(i,proc)
               end if
            end do
         end do
         if( debug ) write(out,*) 'real_vec on root = ',
     &                            (real_vec(i),i=1,size)
c
         deallocate ( temp_vec )
c
      else
         call MPI_SEND(real_vec, size, MPI_REAL, 0, 1,
     &                 MPI_COMM_WORLD, ierr)
      end if
c
      return
 9000 format ('>>>  FATAL ERROR:  could not allocate additional',/,
     &        '>>>       memory in wmpi_reduce_vec.')
c
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_bcast_log               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/05/98                   *
c     *                                                              *
c     *      Broadcast a logical variable from root to all the slave *
c     *      processors.                                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_bcast_log ( log_var )
      implicit integer (a-z)
      include "mpif.h"
      logical log_var
c
      call MPI_BCAST(log_var,1,MPI_LOGICAL,0,MPI_COMM_WORLD, ierr)
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
c     *                   last modified : 05/10/99                   *
c     *                                                              *
c     *      Broadcast a logical variable from root to all the slave *
c     *      processors.                                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_bcast_string ( string, nchars )
      implicit integer (a-z)
      include "mpif.h"
      character *(*) string
c
      call MPI_BCAST(string,nchars,MPI_CHARACTER,0,MPI_COMM_WORLD,
     &               ierr)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_do_external_db          *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 615/2013                   *
c     *                                                              *
c     *         call for UEXTERNALDB Abaqus support routine          *
c     *                                                              *
c     ****************************************************************
c      
      subroutine wmpi_do_uexternaldb
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      dimension status (MPI_STATUS_SIZE)
#dbl      double precision
#sgl      real
     &     zero, aba_time(2), aba_dtime
      logical local_debug
      data zero / 0.0d00 /
c
      local_debug = .false.
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
      return
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
c     *                   last modified : 06/15/01                   *
c     *                                                              *
c     *       This subroutine allows the root processor to send      *
c     *       basic model data, such as the coordinates, incidences, *
c     *       etc., to all of the slave processors.  This data is    *
c     *       only sent once during the analysis.                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_basic
c
      use main_data, only: incmap, incid, crdmap, dtemp_nodes,
     &      dtemp_elems, temper_nodes, temper_elems, invdst,
     &      temper_nodes_ref, temperatures_ref, fgm_node_values,
     &      fgm_node_values_defined, fgm_node_values_cols,
     &      matprp, lmtprp, imatprp, dmatprp, smatprp,
     &      nonlocal_analysis, asymmetric_assembly

      use elem_block_data, only: edest_blocks, cdest_blocks,
     &      edest_blk_list, cdest_blk_list
      use segmental_curves
      use contact, only : use_contact
c
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      dimension status (MPI_STATUS_SIZE)
#dbl      double precision
#sgl      real
     &     zero
      allocatable element_node_counts(:)
      data zero / 0.0d0 /
c
c         tell slaves we are about to send them basic data
c
      call wmpi_alert_slaves ( 2 )
c
      write (out,'("=> proc ",i3," is doing basic data transfer")')myid
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
      call MPI_BCAST(geonl,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
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
      call MPI_BCAST (use_contact,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nonlocal_analysis,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_BCAST(asymmetric_assembly,1,MPI_LOGICAL,0,
     &      MPI_COMM_WORLD, ierr)
      if( nonlocal_analysis .and. myid .eq. 0 ) then
         write(out,9100)
         call die_abort
      end if
c
c             static arrays:
c
      call MPI_BCAST(props,mxelpr*noelem,MPI_REAL,0,MPI_COMM_WORLD,
     &     ierr)
      call MPI_BCAST(cp,mxedof,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(icp,mxutsz*2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dcp,mxedof,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(elblks,4*nelblk,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(c,nodof,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dstmap,nonode,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(num_seg_points,max_seg_curves,MPI_INTEGER,0,
     &     MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curves,max_seg_curves*max_seg_points*2,
     &     MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curves_min_stress,max_seg_curves,
     &     MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curve_def,max_seg_curves,MPI_LOGICAL,0,
     &     MPI_COMM_WORLD,ierr)
      call MPI_BCAST(seg_curves_type,max_seg_curves,MPI_INTEGER,0,
     &     MPI_COMM_WORLD,ierr)
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
c                  if slave, zero out displacements
c
      if ( slave_processor ) then
         do i = 1, nodof
            u(i) = zero
         end do
      end if
c
c             allocated arrays:
c
c                  incidence data structures, inverse incidence
c                  dta structures, inverse dof maps
c
      if (  slave_processor ) then
         call mem_allocate(9)
         call init_maps ( 0, 1 )
         call init_maps ( 0, 2 )
      endif
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
      if ( root_processor ) then
        do elem = 1, noelem
           element_node_counts(elem) = iprops(2,elem)
        end do
      end if
      call MPI_BCAST( element_node_counts, noelem, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr)
      if ( slave_processor ) call setup_slave( element_node_counts )
      deallocate( element_node_counts )
c
c                  functionally graded material properties at
c                  model nodes
c
      if ( fgm_node_values_defined ) then
        if (  slave_processor ) call mem_allocate(20)
        call MPI_BCAST(fgm_node_values, nonode*fgm_node_values_cols,
     &                 MPI_REAL,0, MPI_COMM_WORLD, ierr)
      end if
c
c                  crdmap
c
      if (  slave_processor ) call mem_allocate(14)
      call MPI_BCAST(crdmap,nonode,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
c                  temperature loadings
c
      if (  slave_processor ) then
         call mem_allocate(1)
         call mem_allocate(2)
         call mem_allocate(16)
      endif
      call MPI_BCAST(dtemp_nodes,nonode,MPI_VAL,0,MPI_COMM_WORLD,
     &     ierr)
      call MPI_BCAST(dtemp_elems,noelem,MPI_VAL,0,MPI_COMM_WORLD,
     &     ierr)
      call MPI_BCAST(temper_nodes,nonode,MPI_VAL,0,MPI_COMM_WORLD,
     &     ierr)
      call MPI_BCAST(temper_nodes_ref,nonode,MPI_VAL,0,MPI_COMM_WORLD,
     &     ierr)
      call MPI_BCAST(temper_elems,noelem,MPI_VAL,0,MPI_COMM_WORLD,
     &     ierr)
c
c                  elems_to_blocks structure
c
      if (  slave_processor ) call init_eblock_map
c
c                  allocate space for mdiag
c
      if (  slave_processor  ) then
         call mem_allocate ( 12 )
      endif
c
c             set up data distribution information so all processors
c             know what blocks of elements they own.
c
      call wmpi_calc_dist
c
c             now communicate allocated blocked structures:
c
c                call the standard initialization routines for these
c                structures, so that each will contain the proper starting
c                information.  Note that the root processor has full
c                copys of each of these arrays; however, the slave processors
c                only store the information pertaining to the elements they
c                own.
c
      if ( slave_processor ) then
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
c     *                   last modified : 06/2/2015 rhd              *
c     *                                                              *
c     *       send data from anaylsis parameters to all the MPI      *
c     *       processors                                             *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_analysis
      use main_data,       only : umat_serial
      implicit integer (a-z)
      include "mpif.h"
$add common.main
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
c            values from analysis parameters
c
      call MPI_BCAST(eps_bbar,1,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(qbar_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mxlitr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restrt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dt,1,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nbeta,1,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(convrg,mxcvtests,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(signal_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(adaptive_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(umat_serial,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c
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
c     *                   last modified : 02/27/98                   *
c     *                                                              *
c     *       send data about constraints to all the MPI slave       *
c     *       processors.                                            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_const
      use main_data, only: trn, trnmat
      implicit integer (a-z)
      include "mpif.h"
$add common.main
c
c                   tell slaves we are about to send them data about
c                   the constraints.  Note that basic data must be sent
c                   before we can call this routine -- each processor
c                   must have the number of elements in the structure,
c                   number of nodes, etc.
c
      call wmpi_alert_slaves ( 10 )
c
c      write (out,'("=> proc ",i3," is doing const data transfer")')myid
c
c         broadcast values from constraint input
c
c            constants:
c
      call MPI_BCAST(csthed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
c            static arrays:
c
      call MPI_BCAST(cstmap,nodof,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
c            allocated arrays:
c
      if (.not. allocated (trn)) call mem_allocate(3)
      call MPI_BCAST(trn,nonode,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &     ierr)
c
      do node = 1, nonode
         if ( trn(node) ) then
            if ( slave_processor ) call allo_trnmat (node,1,dum)
            call MPI_BCAST(trnmat(node)%mat,9,MPI_VAL,0,
     &           MPI_COMM_WORLD,ierr)
         endif
      enddo
c
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
c     *                   last modified : 02/27/98                   *
c     *                                                              *
c     *       send to slave processors the information needed for    *
c     *       each newton iteration.                                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_itern
      implicit integer (a-z)
      include "mpif.h"
$add common.main
c
c                   tell slaves we are about to send them data which
c                   they need during the newtons iterations.  This
c                   includes the incremental displacements as determined
c                   from the equation solve.
c
      call wmpi_alert_slaves ( 24 )
c
c         broadcast values
c
      call MPI_BCAST(du,nodof,MPI_VAL,0,MPI_COMM_WORLD,ierr)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_step               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/27/98                   *
c     *                                                              *
c     *       send to slave processors the information needed for    *
c     *       each load step. This includes displacements            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_step
      use main_data, only: asymmetric_assembly
      use adaptive_steps, only : adapt_disp_fact, adapt_temper_fact,
     &                           adapt_load_fact
      implicit integer (a-z)
      include "mpif.h"
$add common.main
#dbl      double precision
#sgl      real
     &     mag, zero
      data zero /0.0/
c
c                   tell slaves we are about to send them data they
c                   need for this newton's iteration.
c
      call wmpi_alert_slaves ( 12 )
c
c         broadcast values from constraint input
c
c            constants:
c
      call MPI_BCAST(tot_anal_time,1,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dt,1,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(scaling_factor,1,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(adapt_temper_fact,1,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(adapt_disp_fact,1,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(adapt_load_fact,1,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(asymmetric_assembly,1,MPI_LOGICAL,0,
     &      MPI_COMM_WORLD,ierr)
c
c            static arrays:
c
      call MPI_BCAST(u,nodof,MPI_VAL,0,MPI_COMM_WORLD,ierr)
c
c            allocated arrays:
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                subbroutine wmpi_send_temp_eqloads            *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 06/15/02                   *
c     *                                                              *
c     *       send to slave processors the temperature data for      *
c     *       current step. send element equivalent nodal loads      *
c     *       if they exist.                                         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_temp_eqloads
      use main_data, only: dtemp_nodes, dtemp_elems, eq_node_force_len,
     &                     eq_node_force_indexes, eq_node_forces
      implicit integer (a-z)
      include "mpif.h"
$add common.main
#dbl      double precision
#sgl      real
     &     mag, zero
      data zero /0.0/
c
c                   tell slaves we are about to send them data they
c                   need for this step
c
      call wmpi_alert_slaves ( 22 )
c
c      write (out,'("=> proc ",i3," is doing step data transfer")')myid
c
c            send logical flag indicating if there are any temperatures
c            to be considered in the current loading.
c
      call MPI_BCAST(temperatures,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c
c            send the incremental element and nodal temperatures,
c            if there are any.
c
      if ( temperatures ) then
c
         call MPI_BCAST(dtemp_nodes,nonode,MPI_VAL,0,MPI_COMM_WORLD,
     &        ierr)
         call MPI_BCAST(dtemp_elems,noelem,MPI_VAL,0,MPI_COMM_WORLD,
     &        ierr)
c
      end if
c
c            send the global integer flag which sets the length of the
c            packed vectors of element equivalent nodal forces for the
c            step. then send the indexing vector and the packed
c            vector of equiv. force values.
c
      call MPI_BCAST( eq_node_force_len, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      if( eq_node_force_len .gt. 0 ) then
        if( slave_processor ) call mem_allocate( 23 )
        call MPI_BCAST( eq_node_force_indexes, noelem,
     &           MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
        call MPI_BCAST( eq_node_forces, eq_node_force_len,
     &           MPI_VAL, 0, MPI_COMM_WORLD, ierr )
      end if
c
c
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine wmpi_send_contact               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 05/10/04                   *
c     *                                                              *
c     *       This subroutine sends the contact information to all   *
c     *       processors.                                            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_contact( restart )
      use contact
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      dimension status (MPI_STATUS_SIZE)
      logical debug, restart, referenced, ldum
#dbl      double precision
#sgl      real
     &     zero
      data debug, zero /.false., 0.0/
c
c         tell slaves we are about to send them contact data
c
      call wmpi_alert_slaves ( 28 )
c
      if(debug)write(out,'("=> proc ",i3," do contact data trans")')myid
c
c         broadcast contact information:
c
c             constants:
c
      call MPI_BCAST (use_contact,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c
      if ( .not. use_contact) then
         if ( debug) write (*,*) myid,': no use contact, baby.'
         return
      endif
c
      call MPI_BCAST(num_contact,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(restart,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c
      if ( debug) write (*,*) myid,': -> num_contact = ', num_contact
c
c             static arrays:
c
      call MPI_BCAST(cplane_vec,3*2*num_contact,MPI_VAL,0,
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(cshape_norm,3*num_contact,MPI_VAL,0,
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(cshape_pnt,3*num_contact,MPI_VAL,0,
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(cshape_rate,3*num_contact,MPI_VAL,0,
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(cshape_param,10*num_contact,MPI_VAL,0,
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(contact_depth,num_contact,MPI_VAL,0,
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(contact_stiff,num_contact,MPI_VAL,0,
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(contact_fric,num_contact,MPI_VAL,0,
     &     MPI_COMM_WORLD, ierr)
      call MPI_BCAST(contact_shape,num_contact,MPI_INTEGER,0,
     &     MPI_COMM_WORLD, ierr)
c
c             if we just came back from a restart, we must also
c             broadcast contact cause.  We just broadcast the entire
c             array to all processors, then let them clear out the
c             entries that they don't reference.  This is an extremely
c             clumsy and inefficient broadcast, but it works...
c             little time of basic thought could make this better.
c
      if( restart ) then
         call MPI_BCAST(contact_cause,maxcontact*nonode,MPI_INTEGER,0,
     &        MPI_COMM_WORLD, ierr)
c
      end if
c
      if ( debug ) then
         write (*,*) '>>> Here are the contact surfaces:'
         do proc = 0, numprocs -1
            if (myid .eq. proc ) then
               write (*,*) '**** proc = ',proc,' of ',numprocs,' ****'
               do shape = 1, num_contact
                  if ( contact_shape(shape) .eq. 1) then
                     write (*,*) '>>>> Contact Params for shape ',shape
                     write (*,*) '   type: ', contact_shape(shape),
     &                    '   stiff:', contact_stiff(shape),
     &                    '   fric: ', contact_fric(shape),
     &                    '   depth:', contact_depth(shape),
     &                    '   rate: ', (cshape_rate(j,shape),j=1,3)
                     write (*,*) '>  input vecs:'
                     do j = 1, 2
                        write (*, *) '     ', cplane_vec (1,j,shape),
     &                       cplane_vec (2,j,shape),
     &                       cplane_vec (3,j,shape)
                     enddo
                     write (*,*) '>  norm vec:'
                     write (*, *) '     ', (cshape_norm (i,shape),i=1,3)
                  else if ( contact_shape(shape) .eq. 2) then
                     write (*,*) '>>>> Contact Params for shape ',shape
                     write (*,*) '   type: ', contact_shape(shape),
     &                    '   stiff:', contact_stiff(shape),
     &                    '   fric: ', contact_fric(shape),
     &                    '   depth:', contact_depth(shape),
     &                    '   rate: ', (cshape_rate(j,shape),j=1,3)
                     write (*,*) '>  base point:'
                     write (*,*) '     ', (cshape_pnt (i,shape),i=1,3)
                     write (*,*) '>  direction:'
                     write (*,*) '     ', (cshape_norm (i,shape),i=1,3)
                     write (*,*) '>  radius:', cshape_param (1,shape),
     &                    '  length:',cshape_param (2,shape)
                  else
                     write (*,*) '>>>> No Contact def for shape ',shape
                  endif
               enddo
               call MPI_BARRIER (MPI_COMM_WORLD, ierr)
            else
               write (*,*) '> --> proc ',myid,' is waiting. '
               call MPI_BARRIER (MPI_COMM_WORLD, ierr)
            endif
         enddo
         write (*,*) '>>> proc ',myid,' is done. '
      endif
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
c     *                   last modified : 02/27/98                   *
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
      use main_data, only: trn, trnmat
      use mpi_lnpcg, only: local_nodes
      use contact
      implicit integer (a-z)
      include "mpif.h"
$add common.main
#dbl      double precision
#sgl      real
     &     zero, mag
      dimension status (MPI_STATUS_SIZE)
      logical debug
      data debug, zero /.false., 0.0/
c
      logical keep_going
c
c
      if(debug)write (out,'("=> proc ",i3," is doing contact gather")')
     &     myid
c
c               reduce the contact_force vector back to the
c               root processor.
c
      call wmpi_reduce_vec (contact_force, nodof)
c
      if (root_processor .and. debug) then
         write (*,*) '>> NONZERO CONTACT FORCE TERMS ON ROOT'
         do i=1, nodof
            if ( contact_force(i) .ne. zero ) then
               write (*,*) '     ',i,'   ',contact_force(i)
            endif
         enddo
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      else
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      endif
c
c               gather transformation matrices
c
      if ( root_processor ) then
c
c                  the root processor loops over the processors.  For each
c	 	   processor, it recieves the number of a node undergoing
c		   contact, then the transformation matrix and cause of
c		   contact for that node.  This continues until all
c		   contacting nodes are processed, then the root processor
c		   moves on to the next processor.
c
         do proc = 1, numprocs - 1
            if ( debug) write (*,*) '======> root get trn from proc ',
     &           proc
c
            keep_going = .true.
            do while ( keep_going )
c
               call MPI_RECV (node, 1, MPI_INTEGER, proc, 1,
     &              MPI_COMM_WORLD, status, ierr)
               if ( node .le. 0 ) then
                  keep_going = .false.
               else
                  if ( .not. trn (node) ) then
                     call allo_trnmat (node, 1, dum )
                     trn ( node ) = .true.
                  endif
                  call MPI_RECV (trnmat(node)%mat, 9, MPI_VAL, proc, 2,
     &                 MPI_COMM_WORLD, status, ierr)
                  call MPI_RECV (contact_cause(1,node), num_contact,
     &                 MPI_INTEGER, proc, 3, MPI_COMM_WORLD, status,
     &                 ierr)
c
                  if ( debug ) then
                     write (*,*) '    > trnmat:'
                     write (*,'(10x,3e13.6)')
     &                    ((trnmat(node)%mat(i,j),j=1,3),i=1,3)
                  endif
c
                  call trn2all (node, 1)
               endif
c
            enddo
c
         enddo
c
      else
c
c                  Each slave processor sends the number of a node which is
c	           undergoing contact, then follows this with the
c		   transformation matric and cause of contact.  It continus
c		   until it has sent all the information about contacting
c		   nodes which it owns, then it sends a 0 asa node number
c		   to indicate that it is finished.
c
         if (debug) write (*,*) '======>',myid,' is sending trn to root'
         do i = 1, local_nodes%num_private
            node = local_nodes%private(i)
            if ( contact_cause(1,node) .eq. 0 ) cycle
            if ( .not. trn(node)) cycle
            call MPI_SEND (node, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD,
     &           ierr)
            call MPI_SEND (trnmat(node)%mat, 9, MPI_VAL, 0, 2,
     &           MPI_COMM_WORLD, ierr)
            call MPI_SEND (contact_cause(1,node), num_contact,
     &           MPI_INTEGER, 0, 3, MPI_COMM_WORLD, ierr)
         enddo
         call MPI_SEND (-1, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ierr)
c
      endif
      if ( debug) write (*,*) '<<=====',myid,' is done w/ trn transfer'
c
      if ( debug ) then
         do proc = 0, numprocs - 1
            if ( myid .eq. proc ) then
               write (*,*) '********* contact cause for proc:',proc
               do i=1, nonode
                  if ( contact_cause(1,i) .ne. 0) then
                     write (*,*) i,'    ',(contact_cause(j,i),j=1,
     &                    num_contact)
                  endif
               enddo
               call MPI_BARRIER (MPI_COMM_WORLD,ierr)
            else
               call MPI_BARRIER (MPI_COMM_WORLD,ierr)
            endif
         enddo
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
c     *                   last modified : 01/07/98                   *
c     *                                                              *
c     *       send initial crack growth data to all the MPI          *
c     *       processors                                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_growth_init
c
      use elem_extinct_data, only: dam_blk_killed
      use damage_data, only : dam_ptr, num_kill_elem, growth_by_kill
c
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      logical been_here, debug
      save been_here
      data been_here, debug /.false., .false./
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
c           tell the slaves to initialize the crack growth data structures.
c
      call wmpi_alert_slaves ( 25 )
c
c           send the number of killable elements and the index vector
c           which goes from element number to the entry in the crack
c           growth data strucutres.
c
      call MPI_BCAST(num_kill_elem,1,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_BCAST(dam_ptr,noelem,MPI_INTEGER,0,MPI_COMM_WORLD,
     &               ierr)
      if( slave_processor ) growth_by_kill = .true.
c
c           if we are one of the MPI slave processors, go allocate
c           the other crack growth data structures we need
c
      if( slave_processor ) call allocate_damage ( 1 )
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
c     *                   last modified : 03/25/98                   *
c     *                                                              *
c     *           send crack growth data to all the MPI slave        *
c     *           processors after root has completed processing     *
c     *           at start of the step                               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_growth ( killed_this_time )
c
      use elem_extinct_data, only: dam_blk_killed, dam_state
      use main_data, only: elems_to_blocks
      use damage_data, only : dam_ptr, growth_by_kill, num_kill_elem
c
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      logical killed_this_time, debug
#sgl      real
#dbl      double precision
     &      zero
      data zero, debug /0.0, .false./
c
c           tell the slaves that we are about to send the the current
c           crack growth information.
c
      call wmpi_alert_slaves ( 13 )
c
      if( debug ) write (out,'("=> proc ",i3," do grwth data trans")')
     &     myid
c
c           send them the state vectors describing what elements
c           and blocks have been killed.
c
      call MPI_BCAST(dam_blk_killed,nelblk,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &               ierr)
      call MPI_BCAST(dam_state,num_kill_elem,MPI_INTEGER,0,
     &               MPI_COMM_WORLD,ierr)
c
c           let MPI processors know if an element has been killed this time
c
      call MPI_BCAST(killed_this_time,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &               ierr)
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
c     *       slave processors own and send it back to root so       *
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
c
      use main_data, only: elems_to_blocks
      use elem_extinct_data, only: dam_ifv, dam_state, dam_blk_killed
      use elem_block_data, only: history_blocks, urcs_n_blocks,
     &                           eps_n_blocks, urcs_n1_blocks,
     &                           element_vol_blocks, history_blk_list
      use damage_data, only : dam_ptr, growth_by_kill
c
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      dimension status ( MPI_STATUS_SIZE )
#dbl      double precision
#sgl      real
     &     mag, zero
      data zero /0.0/
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
c     *                   last modified : 01/23/98                   *
c     *                                                              *
c     *       compute information about the distribution of          *
c     *       elements, nodes, and degrees of freedom across         *
c     *       the MPI processors.                                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_calc_dist
      use elem_block_data, only: edest_blocks
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      logical sent, access_dof(mxdof), open, assign, debug, fail
      dimension size_dof_chunks(mxdof), start_dof_chunks(mxdof),
     &    status(MPI_STATUS_SIZE), blk_ptr(0:max_procs-1)
      character *5 name, dums
      data sent, debug /.false., .false./
      save sent
      real dumr
#dbl      double precision
#sgl      real
     &   dumd
c
      if (debug)write (out,'("=> proc ",i3," is doing dist calcs")')myid
c
c               make sure we only execute this code once per run.
c
      if (sent) return
      sent = .true.
      fail = .false.
c
c              first check if the user specified
c              the processor assignment of the blocks in the blocking
c              part of the input file.  If not, do a round robin
c              assignment of the blocks to the different processors.
c
      assign = .false.
      do i=1, nelblk
         if ( elblks(2,i) .lt. 0 ) then
            assign = .true.
            exit
         endif
      enddo
c
      if ( assign ) then
c
         call errmsg ( 313, dum, dums, dumr, dumd )
         do i=1, nelblk
            elblks(2,i) = mod(i,numprocs)
         enddo
c
      else
c
c                 if processor assignment was given, then we must check
c                 to make sure it is valid -- doesnt refer to more
c                 processors than we have, etc.
c
         do i=1, nelblk
c
            if (elblks(2,i) .gt. numprocs - 1 ) then
               if ( .not. fail .and. root_processor ) then
                  call errmsg ( 314, dum, dums, dumr, dumd )
               endif
               fail = .true.
            endif
c
         enddo
c
      endif
c
      call MPI_BARRIER ( MPI_COMM_WORLD, ierr )
      if ( fail ) then
         call MPI_FINALIZE ( ierr )
         stop
      endif
c
c              build pointer structure onto elblks which points to the
c              next block in sequence which belongs to the processor
c              which owns the current block.  Thus if blocks 5 and 8 were
c              owned by processor 2, then elblks(3,5) = 8, because 8 is
c              the next block which belongs to processor 2.
c
      blk_ptr_head(0:numprocs-1) = -1
      blk_ptr(0:numprocs-1) = -1
c
      do blk = 1, nelblk
         owner_proc = elblks(2,blk)
         if ( blk_ptr_head(owner_proc) .eq. -1) then
            blk_ptr_head(owner_proc) = blk
         else
            elblks(3,blk_ptr(owner_proc)) = blk
         endif
         elblks(3,blk) = -1
         blk_ptr(owner_proc) = blk
      enddo
c
c              print out elblks structure if debug
c
      if ( debug .and. root_processor ) then
         write (out,*) '>>>>>>  elblks: (span,felem,proc,nextblk)'
         do i=1, nelblk
            write (out,'(8x,i5,":",4(1x,i7))') i,(elblks(j,i),j=0,3)
         enddo
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
      if ( root_processor ) goto 1000
c
      access_dof(1:nodof) = .false.
      do blk = 1, nelblk
c
         if ( myid .ne. elblks(2,blk) ) cycle
c
         felem         = elblks(1,blk)
         num_enodes    = iprops(2,felem)
         num_enode_dof = iprops(4,felem)
         span          = elblks(0,blk)
c
         do elem = 1, span
c            write (out,*) myid,':   elem:',elem
            do i = 1, num_enodes * num_enode_dof
c               write (out,*) myid,':       dof:',
c     &           edest_blocks(blk)%ptr(i,elem)
               access_dof ( edest_blocks(blk)%ptr(i,elem) ) = .true.
            enddo
c
         enddo
      enddo
c
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
         if ( dof .gt. nodof ) then
            if ( open ) then
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
         if (access_dof(dof)) then
            if (.not. open) then
               block = block + 1
               start_dof_chunks ( block ) = dof - 1
               open = .true.
            endif
         else
c
c                    this dof is not needed by processor. If last one
c                    was, then set the size of the chunk.
c
            if ( open ) then
               size_dof_chunks(block) = dof - start_dof_chunks(block)+1
               tot_size = tot_size + size_dof_chunks(block)
               open = .false.
            endif
         endif
c
      enddo
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
      if (  root_processor  ) then
         do proc = 1, numprocs - 1
            call MPI_RECV ( num_dof_chunks, 1, MPI_INTEGER, proc,
     &           proc, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV (num_dof_local(proc), 1, MPI_INTEGER, proc,
     &           proc, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV ( start_dof_chunks, num_dof_chunks,
     &           MPI_INTEGER, proc, proc, MPI_COMM_WORLD, status, ierr)
            call MPI_RECV ( size_dof_chunks, num_dof_chunks,
     &           MPI_INTEGER, proc, proc, MPI_COMM_WORLD, status, ierr)
            call MPI_TYPE_INDEXED ( num_dof_chunks, size_dof_chunks,
     &           start_dof_chunks, MPI_VAL, MPI_DOF_LOCAL(proc), ierr)
            call MPI_TYPE_COMMIT ( MPI_DOF_LOCAL(proc), ierr )
         enddo
      else
c
         call MPI_SEND ( num_dof_chunks, 1, MPI_INTEGER, 0, myid,
     &        MPI_COMM_WORLD, ierr)
         call MPI_SEND ( num_dof_local(myid), 1, MPI_INTEGER, 0, myid,
     &        MPI_COMM_WORLD, ierr)
         call MPI_SEND ( start_dof_chunks, num_dof_chunks,
     &        MPI_INTEGER, 0, myid, MPI_COMM_WORLD, ierr)
         call MPI_SEND ( size_dof_chunks, num_dof_chunks,
     &        MPI_INTEGER, 0, myid, MPI_COMM_WORLD, ierr)
         call MPI_TYPE_INDEXED ( num_dof_chunks, size_dof_chunks,
     &        start_dof_chunks, MPI_VAL, MPI_DOF_LOCAL(myid), ierr)
         call MPI_TYPE_COMMIT ( MPI_DOF_LOCAL(myid), ierr )
      endif
c
      if (debug) write (out,*) myid,':>> hey, weza levin wmpi_calc_dist'
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
c     *                   last modified : 03/02/98                   *
c     *                                                              *
c     *     In cases where slave processors are waiting in the slave *
c     *     processor and the root processor determines that warp3d  *
c     *     should end, this routine allows a graceful exit of MPI.  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine die_gracefully
      implicit integer (a-z)
      include "mpif.h"
#dbl      double precision
#sgl      real
     &   aba_time(2), aba_dtime, dbl_val, dt, total_model_time
       real sgl_val
c
c         wake up slaves if they are asleep
c
      call wmpi_suspend(2)
c
c         tell slave processors to quit
c
      call wmpi_alert_slaves ( 1 )
c
c         tell MPI that we are ending now
c
      call MPI_FINALIZE(ierr)
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
c     *                   last modified : 03/02/98                   *
c     *                                                              *
c     *     This subroutine tries to abort all MPI processes in an   *
c     *     emergency situation where slave processors may not be    *
c     *     able to return to the slave handler before aborting.     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine die_abort
      implicit integer (a-z)
      include "mpif.h"
c
c         force all processes to quit
c
      call MPI_ABORT(MPI_COMM_WORLD, dum, ierr)
c
      stop
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_combine_stf             *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 12/10/13 mcm               *
c     *                                                              *
c     *     this subroutine gathers the block element stiffnesses    *
c     *     from the various MPI slave processes and combines them   *
c     *     on the root processor.                                   *
c     *                                                              *
c     *     Modified for parallel assembly                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_combine_stf
      use elem_block_data, only: estiff_blocks
      use main_data, only: asymmetric_assembly
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      character *80 line
c
c                       local declarations
c
      logical local_debug, ldummy
      dimension status(MPI_STATUS_SIZE)
#dbl      double precision
#sgl      real
     &    start_move_time, end_move_time
      data local_count / 0 /
c
c                       if there is only one processor, then skip
c
      if ( numprocs .eq. 1 ) return
c
c          on root, we must create blocks for element stiffness owned
c          by the slaves.
c
      if ( root_processor ) call estiff_allocate( 4 )
c
c          let MPI slaves know we are gathering element stiffnesses
c
      call wmpi_alert_slaves ( 7 )
c
c                       get all of the element block info back to
c                       the root processor
c
      start_move_time = MPI_WTIME()
      do blk = 1, nelblk
c
         if( slave_processor ) then
            if( elblks(2,blk) .ne. myid ) cycle
         else
            if( elblks(2,blk) .eq. myid ) cycle
         endif
c
         felem          = elblks(1,blk)
         num_enodes     = iprops(2,felem)
         num_enode_dof  = iprops(4,felem)
         totdof         = num_enodes * num_enode_dof
         span           = elblks(0,blk)
         utsz           = ((totdof*totdof)-totdof)/2 + totdof
c
         if (.not. asymmetric_assembly) then
         if ( root_processor ) then
            call MPI_RECV( estiff_blocks(blk)%ptr(1,1),
     &                     span*utsz, MPI_VAL, elblks(2,blk),
     &                     14, MPI_COMM_WORLD, status, ierr )
         else
            call MPI_SEND( estiff_blocks(blk)%ptr(1,1), span*utsz,
     &                     MPI_VAL, 0, 14, MPI_COMM_WORLD, ierr )
         end if
         else
         if ( root_processor ) then
            call MPI_RECV( estiff_blocks(blk)%ptr(1,1),
     &                     span*totdof*totdof, MPI_VAL, elblks(2,blk),
     &                     14, MPI_COMM_WORLD, status, ierr )
         else
            call MPI_SEND( estiff_blocks(blk)%ptr(1,1), 
     &                     span*totdof*totdof,
     &                     MPI_VAL, 0, 14, MPI_COMM_WORLD, ierr )
         end if
         end if
      end do
c
      end_move_time = MPI_WTIME()
      if ( root_processor .and. local_count .eq. 0 ) then
         call iodevn( idummy, iout, ldummy, 1 )
         write(iout,9010) end_move_time - start_move_time
         local_count = 1
      end if
c
      return
 9010 format(10x,'>>> walltime to move elem stiffs:',f14.3 )
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_get_str                 *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 12/25/99                   *
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
c     *        3 = gather trial elastic stress vectors on root       *
c     *        4 = gather element strains on root                    *
c     *        5 = gather element volumes on root                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_get_str ( do )
      use elem_block_data
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      character *80 line
c
c                    local declarations
c
      logical local_debug, root_has_stress, root_has_strain
      save root_has_stress, root_has_strain
      dimension status(MPI_STATUS_SIZE)
#dbl      double precision
#sgl      real
     &     zero
#sgl      data zero / 0.0 /
#dbl      data zero / 0.0d00 /
c
c                    if there is only one processor, then skip
c
      if( numprocs .eq. 1 ) return
c
c                    if we are the root processor, then do some checks.
c
      if ( root_processor ) then
c
c                       if do = 1, then the stresses and strains
c                       that the root processor
c                       has for blocks it does not own are out of date.
c                       Thus next time we return to this routine looking
c                       to get the full stresses, we will have to do
c                       communication.
c
         if( do .eq. 1 ) then
            root_has_stress = .false.
            root_has_strain = .false.
            return
c
c                       if do = 2 , then we are getting the major stress
c                       data. This data is needed by several different
c                       algorithms, including stress output, domain
c                       integral calculations, and database construction.
c                       if root already has the current stresses for this
c                       step, then return.  otherwise go get the stresses.
c
         else if( do .eq. 2 ) then
            if( root_has_stress ) then
               return
            end if
            root_has_stress = .true.
c
c                       if do = 3, then we need to get the elastic stress
c                       vectors.  This is only needed for database
c                       construction, so we don't have to worry about
c                       getting these terms twice.
c
         else if( do .eq. 3 ) then
c
c                       if do = 4 , then we are getting the strain
c                       data. This data is needed by several different
c                       algorithms, including strain output, crack
c                       growth calculation, and database construction.
c                       if root already has the current strains for this
c                       step, then return.  otherwise go get the strains.
c
         else if( do .eq. 4 ) then
            if( root_has_strain ) then
               return
            end if
            root_has_strain = .true.
c
         end if
c
      end if
c
c                   if we are here, then we need to get the strains,
c                   stresses, volume data
c                   from the other processors. let MPI slaves know
c                   we are gathering stresses
c
      call wmpi_alert_slaves ( 20 )
      call wmpi_bcast_int ( do )
c
c                   get all of the element block info back to
c                   the root processor
c
      do blk = 1, nelblk
c
         if( slave_processor ) then
            if( elblks(2,blk) .ne. myid ) cycle
         else
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
         if( do .eq. 2 ) then
c
c                               element [D] tangents
c
            block_size = span * ngp * cepsize
            if ( root_processor ) then
               cep_blocks(blk)%vector(1:block_size) = zero
               call MPI_RECV(cep_blocks(blk)%vector(1),block_size,
     &              MPI_VAL,elblks(2,blk),14,MPI_COMM_WORLD,status,ierr)
            else
               call MPI_SEND(cep_blocks(blk)%vector(1),block_size,
     &              MPI_VAL,0,14,MPI_COMM_WORLD,ierr)
            endif
c
c                               element histories
c
            block_size = span * ngp * histsize
            if ( root_processor ) then
               history_blocks(blk)%ptr(1:block_size) = zero
               call MPI_RECV(history_blocks(blk)%ptr(1),block_size,
     &              MPI_VAL,elblks(2,blk),14,MPI_COMM_WORLD,status,ierr)
            else
               call MPI_SEND(history_blocks(blk)%ptr(1),block_size,
     &              MPI_VAL,0,14,MPI_COMM_WORLD,ierr)
            endif
C
c                               element stresses
c
            block_size    = span * ngp * nstrs
            if (root_processor) then
               urcs_n_blocks(blk)%ptr(1:block_size) = zero
               call MPI_RECV(urcs_n_blocks(blk)%ptr(1),block_size,
     &              MPI_VAL,elblks(2,blk),14,MPI_COMM_WORLD,status,ierr)
            else
               call MPI_SEND(urcs_n_blocks(blk)%ptr(1),block_size,
     &              MPI_VAL,0,14,MPI_COMM_WORLD,ierr)
            endif
c
c                               element rotations
c
            if ( rot_blk_list(blk) .eq. 0) cycle
            block_size       = span * ngp * 9
            if (root_processor) then
               if ( allocated(rot_n1_blocks)) then
                  rot_n1_blocks(blk)%ptr(1:block_size) = zero
                  call MPI_RECV(rot_n1_blocks(blk)%ptr(1),block_size,
     &                 MPI_VAL,elblks(2,blk),15,MPI_COMM_WORLD,status,
     &                 ierr)
               endif
            else
               if ( allocated(rot_n1_blocks)) then
                  call MPI_SEND(rot_n1_blocks(blk)%ptr(1),block_size,
     &                 MPI_VAL,0,15,MPI_COMM_WORLD,ierr)
               endif
            endif
c
c                         if do = 3, <available>
c
         else if ( do .eq. 3 ) then
c
c                         if do = 4, get strains
c
         else if ( do .eq. 4 ) then
c
            block_size = span * ngp * nstr
            if (root_processor) then
               eps_n_blocks(blk)%ptr(1:block_size) = zero
               call MPI_RECV(eps_n_blocks(blk)%ptr(1),block_size,
     &              MPI_VAL,elblks(2,blk),14,MPI_COMM_WORLD,status,ierr)
            else
               call MPI_SEND(eps_n_blocks(blk)%ptr(1),block_size,
     &              MPI_VAL,0,14,MPI_COMM_WORLD,ierr)
            endif
c
c                         if do = 5, get element volumes
c
         else if ( do .eq. 5 ) then
c
            block_size = span
            if( root_processor ) then
               element_vol_blocks(blk)%ptr(1:block_size) = zero
               call MPI_RECV(element_vol_blocks(blk)%ptr(1),block_size,
     &              MPI_VAL,elblks(2,blk),14,MPI_COMM_WORLD,status,ierr)
            else
               call MPI_SEND(element_vol_blocks(blk)%ptr(1),block_size,
     &              MPI_VAL,0,14,MPI_COMM_WORLD,ierr)
            end if
c
         end if
c
      end do
c
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine wmpi_send_reopen                  *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 06/15/02                   *
c     *                                                              *
c     *     this subroutine sends all the information that the       *
c     *     slaves need after a restart.                             *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_reopen
      use elem_block_data
      use main_data, only: eq_node_force_len, eq_node_force_indexes,
     &                     eq_node_forces
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      character *80 line
c
c                    local declarations
c
      logical local_debug, root_has_stress, root_has_strain
      save root_has_stress, root_has_strain
      dimension status(MPI_STATUS_SIZE)
c
c                    if there is only one processor, then skip
c
      if( numprocs .eq. 1 ) return
c
c                   tell the slaves that we are going to send them
c                   all the information they need after a restart.
c
      call wmpi_alert_slaves ( 23 )
c
c                   First, send out the element blocks of stresses and
c                   strains to the processors which own them.
c
      do blk = 1, nelblk
c
         if ( slave_processor ) then
            if (elblks(2,blk).ne.myid) cycle
         else
            if (elblks(2,blk).eq.myid) cycle
         endif
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
         if ( slave_processor ) then
            call MPI_RECV(history_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,0,14,MPI_COMM_WORLD,status,ierr)
            call vec_ops( history1_blocks(blk)%ptr(1),
     &           history_blocks(blk)%ptr(1),
     &           dummy, block_size, 5 )
         else
            call MPI_SEND(history_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,elblks(2,blk),14,MPI_COMM_WORLD,ierr)
         endif
c
c                               element [D]s
c
         block_size = span * ngp * cep_blk_list(blk)
         if ( slave_processor ) then
            call MPI_RECV(cep_blocks(blk)%vector(1),block_size,
     &           MPI_VAL,0,14,MPI_COMM_WORLD,status,ierr)
         else
            call MPI_SEND(cep_blocks(blk)%vector(1),block_size,
     &           MPI_VAL,elblks(2,blk),14,MPI_COMM_WORLD,ierr)
         endif
c
c                               element stresses
c
         block_size    = span * ngp * nstrs
         if ( slave_processor ) then
            call MPI_RECV(urcs_n_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,0,15,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(urcs_n1_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,0,15,MPI_COMM_WORLD,status,ierr)
         else
            call MPI_SEND(urcs_n_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,elblks(2,blk),15,MPI_COMM_WORLD,ierr)
            call MPI_SEND(urcs_n1_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,elblks(2,blk),15,MPI_COMM_WORLD,ierr)
         endif
c
c
c                               element integration point rotations
c
         if ( rot_blk_list(blk) .ne. 0) then
            block_size       = span * ngp * 9
            if ( slave_processor ) then
               if ( allocated(rot_n1_blocks)) then
                  call MPI_RECV(rot_n1_blocks(blk)%ptr(1),block_size,
     &                 MPI_VAL,0,16,MPI_COMM_WORLD,status,
     &                 ierr)
               endif
            else
               if ( allocated(rot_n1_blocks)) then
                  call MPI_SEND(rot_n1_blocks(blk)%ptr(1),block_size,
     &                 MPI_VAL,elblks(2,blk),16,MPI_COMM_WORLD,ierr)
               endif
            endif
         endif
c
c                               element strains
c
         block_size = span * ngp * nstr
         if ( slave_processor ) then
            call MPI_RECV(eps_n_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,0,17,MPI_COMM_WORLD,status,ierr)
         else
            call MPI_SEND(eps_n_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,elblks(2,blk),17,MPI_COMM_WORLD,ierr)
         endif
c
c                               element volumes
c
         block_size = span
         if ( slave_processor ) then
            call MPI_RECV(element_vol_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,0,17,MPI_COMM_WORLD,status,ierr)
         else
            call MPI_SEND(element_vol_blocks(blk)%ptr(1),block_size,
     &           MPI_VAL,elblks(2,blk),17,MPI_COMM_WORLD,ierr)
         endif
c
      end do
c
c
c                   send the packed vector form of the element
c                   equivalent nodal forces if they exist.
c
c
      call MPI_BCAST( eq_node_force_len, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      if( eq_node_force_len .gt. 0 ) then
        if( slave_processor ) call mem_allocate( 23 )
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
c     *                   last modified : 02/05                      *
c     *                              by : mcw                        *
c     *                                                              *
c     *           send/receive data to processes to set up           *
c     *           j-integral and i-integral computation for          *
c     *           this single domain                                 *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_jint
      use main_data
      use elem_block_data, only : history_blocks, rot_n1_blocks,
     &                            urcs_n_blocks, cdest_blocks,
     &                            edest_blocks, eps_n_blocks
      use j_data
      implicit integer (a-z)
      include "mpif.h"
#dbl      double precision
#sgl      real
     &  zero
$add common.main
      logical logical_vec(15), handle_face_loads, debug
      data zero, debug / 0.0, .false. /
c
c          skip this routine if we only have one processor
c
          if ( numprocs .eq. 1 ) return
c
c         send /receive shared data needed to process this domain.
c         we have two major cases: (1) this is the first domain
c         to be processed and most all data needs to be broadcast,
c         (2) on subsequent domains, we only need to broadcast
c         domain dependent information.
c
c         broadcast domain dependent information.
c             1) q-values at nodes
c             2) element list to process (stored as a bitmap)
c
      if ( slave_processor ) then
          allocate( q_values(nonode), stat = alloc_stat )
          if ( alloc_stat .ne. 0 ) then
             write(out,9100)
             call die_abort
             stop
          end if
          allocate( q_element_maps(noelem/30+1), stat = alloc_stat )
          if ( alloc_stat .ne. 0 ) then
             write(out,9100)
             call die_abort
             stop
          end if
          allocate( crack_front_elem(noelem), stat = alloc_stat )
          if ( alloc_stat .ne. 0 ) then
             write(out,9100)
             call die_abort
             stop
          end if
      end if
c
      call MPI_BCAST( q_values, nonode, MPI_REAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( q_element_maps, noelem/30+1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( crack_front_elem, noelem, MPI_LOGICAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( first_domain, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( front_nodes, 30, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( num_front_nodes, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( front_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( domain_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( omit_crack_front_elems, 1, MPI_LOGICAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( max_exp_front, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &                ierr )
      call MPI_BCAST( face_loading, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     &                ierr )
      if( slave_processor ) then
         allocate(expanded_front_nodes(0:max_exp_front,num_front_nodes),
     &        stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
            write(out,9100)
            call die_abort
            stop
         end if
      end if
c
      if ( .not. first_domain ) goto 9999
c
c         processing first domain. broadcast data used for all
c         domains defined at this crack front position.
c
c         a set of logical flags. for slaves, turn off flags
c         that would otherwise let them print intermediate results
c         during j computation.
c
      if( root_processor ) then
         logical_vec(1)  = symmetric_domain
         logical_vec(2)  = q_vals_linear
         logical_vec(3)  = one_point_rule
         logical_vec(4)  = static_j
         logical_vec(5)  = ignore_face_loads
         logical_vec(6)  = process_velocities
         logical_vec(7)  = process_accels
         logical_vec(8) = debug_driver
         logical_vec(9) = debug_elements
         logical_vec(10) = print_elem_values
         logical_vec(11) = out_pstress
         logical_vec(12) = out_pstrain
         logical_vec(13) = cf_traction_flags(1)
         logical_vec(14) = cf_traction_flags(2)
         logical_vec(15) = cf_traction_flags(3)
      end if
c
      call MPI_BCAST( logical_vec, 15, MPI_LOGICAL, 0, MPI_COMM_WORLD,
     &                ierr )
c
      if ( slave_processor ) then
         symmetric_domain     = logical_vec(1)
         q_vals_linear        = logical_vec(2)
         one_point_rule       = logical_vec(3)
         static_j             = logical_vec(4)
         ignore_face_loads    = logical_vec(5)
         process_velocities   = logical_vec(6)
         process_accels       = logical_vec(7)
         debug_driver         = logical_vec(8)
         debug_elements       = logical_vec(9)
         print_elem_values    = logical_vec(10)
         out_pstress          = logical_vec(11)
         out_pstrain          = logical_vec(12)
         cf_traction_flags(1) = logical_vec(13)
         cf_traction_flags(2) = logical_vec(14)
         cf_traction_flags(3) = logical_vec(15)
      end if
c
c             1) 3x3 rotation matrix for this crack front poisiton
c             2) nodal displacements, velocities and accelerations
c
      call MPI_BCAST( domain_rot, 9, MPI_VAL, 0, MPI_COMM_WORLD,
     &                ierr)
      call MPI_BCAST( u, nodof, MPI_VAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST( v, nodof, MPI_VAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST( a, nodof, MPI_VAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST( domain_origin, 1, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( cf_tractions, 3, MPI_VAL, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( front_coords, 90, MPI_VAL, 0, MPI_COMM_WORLD,
     &                ierr)
      num = (max_exp_front + 1) * num_front_nodes
      call MPI_BCAST( expanded_front_nodes, num, MPI_INTEGER, 0,
     &                MPI_COMM_WORLD, ierr )
      call MPI_BCAST( e_front, 1, MPI_VAL, 0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( nu_front, 1, MPI_VAL, 0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( crack_curvature, 7, MPI_VAL, 0,
     &                MPI_COMM_WORLD, ierr )
c
 9999 continue
c
      if(debug)write (*,*) '<<< just leaving wmpi_send_jint,',myid
c
      return
c
 9100 format('>> FATAL ERROR: allocate error in wmpi_send_jint: ')
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_real_new           *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 06/11                      *
c     *                                                              *
c     *     Get the non-zero data in a vector of a given length      *
c     *     onto the root processor                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_real_new(vector,length)
            implicit integer (a-z)
            include "mpif.h"
$add common.main
c                 Dummy
            real :: vector(*)
            integer :: length
c                 Local
            real, dimension(:), allocatable :: temp_vec
            integer :: ierr
c                 Setup receive buffer, reduce via max, copy over
            allocate(temp_vec(length))
            call MPI_REDUCE(vector,temp_vec,length,MPI_REAL,
     &                      MPI_MAX,0,MPI_COMM_WORLD,ierr)
            if (myid .eq. 0) then
                  vector(1:length) = temp_vec(1:length)
            end if
c
            return
c
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_simple_angles      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 3/12/14                    *
c     *                                                              *
c     *           Send the simple angle properties, if required      *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_simple_angles
            use crystal_data, only: srequired, simple_angles, nangles,
     &            mc_array
            implicit integer(a-z)
            include 'mpif.h'
$add common.main
c     
            if (myid .eq. 0) then
                  call wmpi_alert_slaves(49)
            end if

            call MPI_Bcast(srequired, 1, MPI_LOGICAL, 0,
     &            MPI_COMM_WORLD, ierr)

            if (srequired) then
              call MPI_Bcast(nangles, 1, MPI_INTEGER, 0,
     &            MPI_COMM_WORLD, ierr)
              if (myid .ne. 0) then
                if (.not. allocated(simple_angles)) then
                  allocate(simple_angles(nangles,3))
                end if
                if (.not. allocated(mc_array)) then
                  allocate(mc_array(noelem,max_crystals))
                end if
              end if

              call MPI_Bcast(simple_angles, 3*nangles, 
     &            MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
              call MPI_Bcast(mc_array, noelem*max_crystals,
     &            MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            end if

            return

       end subroutine

c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_crystals           *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 8/16/12                    *
c     *                                                              *
c     *           Send all the crystal properties to the workers     *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_crystals
            use crystal_data, only: c_array, angle_input,
     &            crystal_input, data_offset, print_crystal
            implicit integer(a-z)
            include 'mpif.h'
$add common.main
            integer :: nelem, mxcry, ierr
            integer, parameter :: count_struct=5
            integer :: blocklens(0:4), indices(0:4), otypes(0:4)
            integer :: size_int, size_dp, size_log, ntype
c
            if (myid .eq. 0) then
                  call wmpi_alert_slaves(49)

                  nelem = size(crystal_input,1)
                  mxcry = size(crystal_input,2)
            end if

            call MPI_Bcast(nelem, 1, MPI_INTEGER, 0,
     &            MPI_COMM_WORLD, ierr)
            call MPI_Bcast(mxcry, 1, MPI_INTEGER, 0,
     &            MPI_COMM_WORLD, ierr)
            call MPI_Bcast(noelem, 1, MPI_INTEGER, 0,
     &            MPI_COMM_WORLD, ierr)

            if (nelem .eq. 0 ) go to 100
c
            if (myid .ne. 0) then
              if (.not. allocated(data_offset)) then
                allocate(data_offset(noelem))
              end if
              if (.not. allocated(angle_input)) then
                allocate(angle_input(nelem,mxcry,3))
              end if
              if (.not. allocated(crystal_input)) then
                allocate(crystal_input(nelem,mxcry))
              end if
            end if
c           Now get the data out there
c           First the easy ones
            call MPI_Bcast(data_offset,noelem,MPI_INTEGER,0,
     &            MPI_COMM_WORLD, ierr)
            call MPI_Bcast(angle_input,nelem*mxcry*3,
     &            MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,ierr)
            call MPI_Bcast(crystal_input,nelem*mxcry,
     &            MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
c
100         continue
c           Now the struct
c           3 integers
c           25 double precisions
c           2 dp 6x6
c           2 dp 12x3
c           2 logical
c           = 5 blocks
c
            call MPI_Type_extent(MPI_INTEGER,size_int,ierr)
            call MPI_Type_extent(MPI_DOUBLE_PRECISION,size_dp,ierr)
            call MPI_Type_extent(MPI_LOGICAL,size_log,ierr)
c           Set up all the wonderful structs
            otypes(0) = MPI_INTEGER
            otypes(1) = MPI_DOUBLE_PRECISION
            otypes(2) = MPI_DOUBLE_PRECISION
            otypes(3) = MPI_DOUBLE_PRECISION
            otypes(4) = MPI_LOGICAL
c
            blocklens(0) = 3
            blocklens(1) = 25
            blocklens(2) = 2*6*6
            blocklens(3) = 2*12*3
            blocklens(4) = 2
c
            indices(0) = 0
            indices(1) = 3*size_int
            indices(2) = indices(1) + 25*size_dp
            indices(3) = indices(2) + 2*6*6*size_dp
            indices(4) = indices(3) + 2*12*3*size_dp
c
            call MPI_Type_struct(count_struct, blocklens, indices,
     &            otypes, ntype, ierr)
            call MPI_Type_commit(ntype, ierr)
c
c           Actually transmit
            call MPI_Bcast(c_array, max_crystals, ntype,
     &            0, MPI_COMM_WORLD,ierr)


            call MPI_Type_free(ntype,ierr)


c            call MPI_Barrier(MPI_COMM_WORLD,ierr)
c            do i=0,3
c            if (i .eq. myid) then
c            write (*,*) "..............................."
c            write (*,*) myid
c            write (*,*)
c            call print_crystal(1)
c            write (*,*)
c            end if
c            call MPI_Barrier(MPI_COMM_WORLD,ierr)
c            end do
c            call MPI_Barrier(MPI_COMM_WORLD,ierr)
c            call die_gracefully

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_dealloc_crystals        *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 8/16/12                    *
c     *                                                              *
c     *           Dealloc crystal data structures                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_dealloc_crystals
            use crystal_data, only: c_array, angle_input,
     &            crystal_input, data_offset, mc_array,
     &            simple_angles
            implicit integer (a-z)
            include 'mpif.h'
$add common.main
c
            if (myid .eq. 0) then
             call wmpi_alert_slaves(50)
            else
             if (allocated(angle_input)) deallocate(angle_input)
             if (allocated(crystal_input)) deallocate(crystal_input)
             if (allocated(data_offset)) deallocate(data_offset)
             if (allocated(simple_angles)) deallocate(simple_angles)
             if (allocated(mc_array)) deallocate(mc_array)
            end if

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_init_owner              *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 1/24/2015 rhd              *
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
         call wmpi_allocate_owner(num_private_nodes(proc),
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
c            in the future we may want to know which
c            blocks are completely internal (have no shared
c            nodes). this was part of the ebe pcg solver in the 
c            past.
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
         call wmpi_print_node_own (proc_nodes(i),i)
      enddo
      write (out,'("  -------------------------")')
c
c            fill processor zero's own local_nodes structure -- this
c            is the permanent data stucture for the node ownership
c            data for the root processor
c
      if (debug) write (out,*) '...do proc 0s local_nodes'
      call wmpi_allocate_owner(num_private_nodes(0),
     &     num_own_shared_nodes(0), num_shared_nodes(0),
     &     proc_nodes(0)%num_int_blks, local_nodes)
c
      call wmpi_copy_own (proc_nodes(0), local_nodes)
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
          allocate( 
     &       procdof2glob(proc)%dof(num_local_nodes(proc)*mxndof),
     &       stat = alloc_stat)
         if ( alloc_stat .ne. 0 ) then
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
      if (debug) write (out,*) '...doing multi-processor thang'
      call wmpi_owner_send(num_private_nodes,
     &     num_own_shared_nodes, num_shared_nodes, num_local_nodes)
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
         call wmpi_deallocate_owner (proc_nodes(proc))
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
c     *                   last modified : 01/23/2015 rhd             *
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
         call wmpi_allocate_owner(num_private, num_own_shared,
     &              num_shared, local_nodes%num_int_blks, local_nodes)
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
               call wmpi_print_node_own_all(local_nodes)
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
c      call ledest_init
c
c                  deprecated allocate kdiag, local_mdiag
c
      if ( debug ) write (out,*) myid,': <<- leaving wmpi_owner_send'
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
c     *                   last modified : 01/24/2015 rhd             *
c     *                                                              *
c     *     this subroutine calculates the node owner structure      *
c     *     for domain decomposition portions of the MPI code        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_allocate_owner(num_priv, num_own_shar,
     &         num_shar, num_int_blks, own)
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
c     *                 subroutine wmpi_deallocate_owner             *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 1/24/2015                  *
c     *                                                              *
c     *     this subroutine calculates the node owner structure      *
c     *     for domain decomposition portions of the MPI code        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_deallocate_owner(own)
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
      subroutine wmpi_print_node_own_all(own)
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
c     *                 subroutine wmpi_print_node_own               *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 2/24/2015 rhd              *
c     *                                                              *
c     *     this subroutine prints selected information about node   *
c     *     ownership for domain decomposition portions of the MPI   *
c     *     code.                                                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_print_node_own(own, procnum)
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
c     *                      subroutine wmpi_copy_own                *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 01/24/2015                 *
c     *                                                              *
c     *     this subroutine copies the node owner structure          *
c     *     from one structure to another, allocated to the same     *
c     *     size. Assumes that the constants have already been       *
c     *     set in the new copy.                                     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_copy_own(old, new)
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
