c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_init                    *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : rhd 7/3/12                 *
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
#hpi      call wmpi_procids
#h11      call wmpi_procids
#sgi      call wmpi_procids
#sga      call wmpi_procids
#dec      call wmpi_procids
#r60      call wmpi_procids
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
c     *                   last modified : 10/01/01                   *
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
	    enddo
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
c     *                   last modified : 10/01/01                   *
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
#sgi            call pxfkill(proc_pids(proc), 23, ierr)
#sga            ierr = kill (proc_pids(proc), 19)
#l64            ierr = kill (proc_pids(proc), 19)
#hpi            ierr = kill (proc_pids(proc), 24)
#h11            ierr = kill (proc_pids(proc), 24)
#dec            ierr = kill (proc_pids(proc), 17)
#r60            ierr = kill (%val(proc_pids(proc)), %val(17))
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
#sgi            call pxfkill(proc_pids(proc), 25, ierr)
#sga            ierr = kill (proc_pids(proc), 18)
#l64            ierr = kill (proc_pids(proc), 18)
#hpi            ierr = kill (proc_pids(proc), 26)
#h11            ierr = kill (proc_pids(proc), 26)
#dec            ierr = kill (proc_pids(proc), 19)
#r60            ierr = kill (%val(proc_pids(proc)), %val(19))
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
c     *                      subroutine wmpi_handle_slaves           *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : mcm 05/11                  *
c     *                                                              *
c     *     The MPI implementation of warp3d follows a master-slave  *
c     *     approach for most of the code.  The root processor       *
c     *     (processor 0) runs warp3d in full, notifying the slave   *
c     *     processors when they have work to do.                    *
c     *                                                              *
c     *     This routine serves as the waiting point for all of the  *
c     *     slave processors, where they wait for a message from     *
c     *     the root telling them what to do next.                   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_handle_slaves
      implicit integer (a-z)
      include "mpif.h"
$add common.main
      logical duml, debug, ldum1, ldum2, ldum3, ldum4, ldum5
      data debug /.false./
#dbl      double precision
#sgl      real
     &   dumd
      dimension adum1(max_procs), adum2(max_procs), adum3(max_procs),
     &     adum4(max_procs), num_blks(mxnmbl)
c
c           stay in loop, processing the master requests
c
      if(debug)write (out,'("=> Proc",i3," in slave driver")') myid
      do while (.true.)
c
c               wait for a job instruction from the root processor.
c               In this case we call an MPI_BCAST, where the slave
c               processors are waiting for root to broadcast the
c               next instruction.  Note that all slaves are given
c               the same instruction.
c
         if(debug)write (out,'("====> Proc",i3," waiting for do ")')
     &        myid
c
         call MPI_BCAST(do,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
         if(debug)write (out,'("====> Proc",i3," just got ",i3)')
     &        myid, do
c
c               Now branch based on the instruction number sent.
c
         goto ( 100, 200, 300, 400, 500, 600, 700, 800, 900,1000,
     &         1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,
     &         2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,
     &         3100,3200,3300,3400,3500,3600,3700,3800) do
c
         goto 9000
c
c           do = 1: call a routine to kill the slaves (call stop).
c
 100     continue
         if(debug)write (out,'("=> Proc",i3," is exiting")') myid
         call die_gracefully
         goto 9000
c
c           do = 2: send basic model data (coordinates, incidences, etc.)
c
 200     continue
         if(debug)write (out,'("=> Proc",i3," get basic data")') myid
         call wmpi_send_basic
         goto 9000
c
c           do = 3: execute lnstff, the routine which generates the linear
c                   stiffness matrix for each element.  Element stiffnesses
c                   are only calculated for the elements that the processor
c                   owns.
c
 300     continue
         if(debug)write (out,'("=> Proc",i3," do lnstff")') myid
         dumi = 1
         dumj = 1
         call lnstff( dumi, dumj )
         goto 9000
c
c           do = 4: execute tanstf, the routine which generates the tangent
c                   stiffness matrix for each element. Element stiffnesses
c                   are only calculated for the elements that the processor
c                   owns.
c
 400     continue
         if(debug)write (out,'("=> Proc",i3," do tanstf")') myid
         duml = .false.
         dumi = 1
         dumj = 1
         call tanstf( duml, dumi, dumj)
         goto 9000
c
c           do = 5: enter into the linear preconditioned conjugate
c                   gradient solver.  Runs with full domain decomposition,
c                   where each processor owns not only a subset of the
c                   elements, but also a subset of the nodes.
c
 500     continue
         if(debug)write (out,'("=> Proc",i3," do lnpcg")') myid
         call lnpcg(idu,dum1,dum2,dumd)
         goto 9000
c
c           do = 6: call the routine which adds the lumped mass terms
c                   into the element stiffness matricies.  Modifys only
c                   the element stiffness matricies which this processor
c                   owns.
c
 600     continue
         if(debug)write (out,'("=> Proc",i3," go to inclmass")') myid
         call inclmass
         goto 9000
c
c           do = 7: gather the element stiffnesses to the root processor
c                   so that it can assemble a full stiffness matrix for
c                   the direct solvers.
c
 700     continue
         if(debug)write (out,'("=> Proc",i3," gather estf")') myid
         call wmpi_combine_stf
         goto 9000
c
c           do = 8: calculate the lumped mass matrix for all of the elements
c                   which the processor owns.
c
 800     continue
         if(debug)write (out,'("=> Proc",i3," go to cmpmas")') myid
         call cmpmas
         goto 9000
c
c           do = 9: send analysis parameter data from root to all the
c                   slave processors.
c
 900     continue
         if(debug)write (out,'("=> Proc",i3," snd analysis params")')
     &        myid
         call wmpi_send_analysis
         goto 9000
c
c           do = 10: send constraint data from root to all the slave
c                    processors.
c
 1000    continue
         if(debug)write (out,'("=> Proc",i3," send const data")') myid
         call wmpi_send_const
         goto 9000
c
c           do = 11: remove transform. mat. due to contact
c
 1100    continue
         if(debug)write (out,'("=> Proc",i3,"kill cont trnmats")') myid
         call contact_remove (.false.)
         goto 9000
c
c           do = 12: send data needed each load step to all
c                    the slave processors.
c
 1200    continue
         if(debug)write (out,'("=> Proc",i3," send step data")') myid
         call wmpi_send_step
         goto 9000
c
c           do = 13: send info from root to the slaves about the current
c                    amount of crack growth` and whether they should
c                    kill any of their owned elements.
c
 1300    continue
         if(debug)write (out,'("=> Proc",i3," get crk growth info")')
     &        myid
         call wmpi_send_growth ( duml )
         goto 9000
c
c           do = 14: call the routine which computes the stresses for the
c                    elements.  Each processor copmutes the stresses only
c                    for the elements it owns.
c
 1400    continue
         if(debug)write (out,'("=> Proc",i3," do stress calcs")') myid
         call drive_eps_sig_internal_forces(dumi, dumi2, duml)
         goto 9000
c
c           do = 15: call routines that calculate strains and stress work
c                    density at element nodes for J-computations
c
 1500    continue
         if(debug)write (out,'("=> Proc",i3,"do fgm setup")') myid
         call di_fgm_setup(dumi,dumi2,dumi3,duml)
	 goto 9000
c
c           do = 16: setup nodal ownership data structures for the
c                    domain decomposition version of lnpcg
c
 1600    continue
         if(debug)write (out,'("=> Proc",i3," set up dom. dec. mpi")')
     &        myid
         call wmpi_owner_send(adum1,adum2,adum3,adum4)
         goto 9000
c
c           do = 17: call the routine which updates the global solutions
c                    (stresses, strains, total displacements) following
c                    the successful completion of a load step.
c
 1700    continue
         if(debug)write (out,'("=> Proc",i3," do update")') myid
         call update
         goto 9000
c
c           do = 18: save data (stresses, strains, etc.) required for
c                    restarting a non-converged load increment using
c                    the adaptive load reduction mechanism.
c
 1800    continue
         if(debug)write (out,'("=> Proc",i3," adapt save")') myid
         call adaptive_save
         goto 9000
c
c           do = 19: this routine reloads the solutions for the previous
c                    load step when the current load step has failed to
c                    converge, and we are using the adaptive mechanism
c                    to try a smaller load step.
c
 1900    continue
         if(debug)write (out,'("=> Proc",i3," adapt reset")') myid
         call adaptive_reset
         goto 9000
c
c           do = 20: gather the element stresses from the slave processors
c                    and send them to the root processor.  This is necessary
c                    for stress output, restarts, etc.
c
 2000    continue
         if(debug)write (out,'("=> Proc",i3," gather stress")') myid
         call wmpi_get_str ( dum )
         goto 9000
c
c           do = 21: each processor owns a complete copy of the temperatures
c                    this routine scales the temperatures for a load
c                    step due to load multipliers or adaptive load
c                    reduction.
c
 2100    continue
         if(debug)write (out,'("=> Proc",i3," deal w/ temps")') myid
         call mnralg_scale_temps ( dum1, dum2 )
         goto 9000
c
c           do = 22: send temperature data and element equiv nodal forces
c                    from root to the slave processors.
c                    each processor has a complete copy of the
c                    temperature data structures.
 2200    continue
         if(debug)write (out,'("=> Proc",i3," send temp data")') myid
         call wmpi_send_temp_eqloads
         goto 9000
c
c           do = 23: this routines sends the data needed by the slave
c                    processors which were read in by root after a restart.
c
 2300    continue
         if(debug)write (out,'("=> Proc",i3," send reopen")') myid
         call wmpi_send_reopen
         goto 9000
c
c           do = 24: have root send data the slave processors need in the
c                    middle of each newtons iteration (du vector, etc.).
c
 2400    continue
         if(debug)write (out,'("=> Proc",i3," send du")') myid
         call wmpi_send_itern
         goto 9000
c
c           do = 25: have the slave processors initialize the crack
c                    growth data structures and get critical growth
c                    information from the root processor.
c
 2500    continue
         if(debug)write (out,'("=> Proc",i3," init crkgrowth")') myid
         call wmpi_growth_init
         goto 9000
c
c           do = 26: retrieve dam_ifv terms and other information
c                    owned by the slave processors but needed by
c                    root to correctly compute crack growth.
c
 2600    continue
         if(debug)write (out,'("=> Proc",i3," str crkgrowth")') myid
         call wmpi_get_grow
         goto 9000
c
c           do = 27: have each processor compute the parts of the
c                    J-integral domain for the elements which they
c                    own.
c
 2700    continue
         if(debug)write (out,'("=> Proc",i3,"compute J")') myid
         call dicmj
         goto 9000
c
c           do = 28: send each processor the contact information
c
 2800    continue
         if(debug)write (out,'("=> Proc",i3," get contact info")') myid
         call wmpi_send_contact (duml)
         goto 9000
c
c           do = 29: find contact for nodes referenced by each processor
c
 2900    continue
         if(debug)write (out,'("=> Proc",i3,"find contact")') myid
         call contact_find
         goto 9000
c
c           do = 30: allocate data structures for J calculations with
c                    temperature loadings
c
 3000    continue
         if(debug)write (out,'("=> Proc",i3,"all j temp stuff")') myid
         call di_node_props_setup (dum1,dum2)
         goto 9000
c
c           do = 31: distribute graph structures for ebe preconditioner
c
 3100    continue
         if(debug)write (out,'("=> Proc",i3," get ebe graph")') myid
         call wmpi_graph_send (dum1,dum2,num_blks)
         goto 9000
c
c           do = 32: init preconditioner for lnpcg
c
 3200    continue
         if(debug)write (out,'("=> Proc",i3," init prec")') myid
         call updpcm
         goto 9000
c
c           do = 33: allocate pcm matrices
c
 3300    continue
         if(debug)write (out,'("=> Proc",i3," allo pcm")') myid
         call estiff_allocate(2)
         goto 9000
c
c           do = 34: test stuff
c
 3400    continue
         goto 9000
c
c
c           do = 35: patran element stress-strain file
c
 3500    continue
         if(debug)write (out,'("=> Proc",i3," get pat ele")') myid
         call oustr_pat_file( ldum1,ldum2,ldum3,ldum4,ldum5 )
         goto 9000
c
c
c           do = 36: call hypre routine with dummy arguments
c
 3600    continue
            if (debug) write (out,'("=> Proc",i3," enter Hypre")') myid
            call warp_hypre( 0, 0, 0, 0,
     &                        0, 0, 0,
     &                        0, 0 )
            goto 9000
c
c           do =37: call wmpi_set_mkl_threads
c
 3700    continue
            if (debug) write (out,'("=> Proc",i3," enter MKL set")') myid
            call wmpi_set_mkl_threads(0)
            goto 9000
c
c           do =38: call wmpi_set_omp_threads
c
 3800    continue
            if (debug) write (out,'("=> Proc",i3," enter OMP set")') myid
            call wmpi_set_omp_threads(0)
            goto 9000
c
c
 9000 continue
      enddo
c
      return
c
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
     &      matprp, lmtprp, nonlocal_analysis  

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
      call MPI_BCAST(use_contact,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nonlocal_analysis,1,MPI_LOGICAL,0,MPI_COMM_WORLD,
     &               ierr)
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
         call rts_init( 0, 1 )
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
c     *                   last modified : 09/27/12 rhd               *
c     *                                                              *
c     *       send data from anaylsis parameters to all the MPI      *
c     *       processors                                             *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_analysis
      use pcg_data,        only : ebe_pre
      use thread_data,     only : warp_threads, umat_threads
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
      call MPI_BCAST(convrg,10,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lcnvrg,10,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ltol,10,MPI_VAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(signal_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(adaptive_flag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ebe_pre,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(umat_threads,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(warp_threads,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
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
c
      return
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
c     *                   last modified : 01/21/98                   *
c     *                                                              *
c     *     this subroutine gathers the block element stiffnesses    *
c     *     from the various MPI slave processes and combines them   *
c     *     on the root processor.                                   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_combine_stf
      use elem_block_data, only: estiff_blocks
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
         if ( root_processor ) then
            call MPI_RECV( estiff_blocks(blk)%ptr(1,1),
     &                     span*utsz, MPI_VAL, elblks(2,blk),
     &                     14, MPI_COMM_WORLD, status, ierr )
         else
            call MPI_SEND( estiff_blocks(blk)%ptr(1,1), span*utsz,
     &                     MPI_VAL, 0, 14, MPI_COMM_WORLD, ierr )
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
c                         if do = 3, get elastic stress vectors
c
         else if ( do .eq. 3 ) then
c
            if ( rts_blk_list(blk) .eq. 0) cycle
            block_size = span * ngp * nstr
            if (root_processor) then
               rts_blocks(blk)%ptr(1:block_size) = zero
               call MPI_RECV(rts_blocks(blk)%ptr(1),block_size,
     &              MPI_VAL,elblks(2,blk),14,MPI_COMM_WORLD,status,ierr)
            else
               call MPI_SEND(rts_blocks(blk)%ptr(1),block_size,
     &              MPI_VAL,0,14,MPI_COMM_WORLD,ierr)
            endif
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
c                               elements (trial) elastic stress vectors
c
c
         if ( rts_blk_list(blk) .ne. 0) then
            block_size = span * ngp * nstr
            if ( slave_processor ) then
               if ( allocated(rts_blocks)) then
                  call MPI_RECV(rts_blocks(blk)%ptr(1),block_size,
     &                 MPI_VAL,0,14,MPI_COMM_WORLD,status,ierr)
               endif
            else
               if ( allocated(rts_blocks)) then
                  call MPI_SEND(rts_blocks(blk)%ptr(1),block_size,
     &                 MPI_VAL,elblks(2,blk),14,MPI_COMM_WORLD,ierr)
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
c     *                      subroutine wmpi_set_mlk_threads         *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 06/23/11                   *
c     *                                                              *
c     *     Has all MPI ranks set the number of MKL threads          *
c     *     to the provided number.  Leaves you responsible for      *
c     *     tracking what that value should be.                      *
c     *                                                              *
c     ****************************************************************
      subroutine wmpi_set_mkl_threads(num)
            implicit integer (a-z)
            include 'mpif.h'
$add common.main
c           Dummy
            integer :: num
c           Declarations
            integer :: mpierr, numcpy
c           Sync up and set
            if (myid .eq. 0) then
                  numcpy=num
            end if
            call MPI_Bcast(numcpy,1,MPI_INTEGER
     &            ,0,MPI_COMM_WORLD,mpierr)
            call mkl_set_num_threads(numcpy)
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_set_omp_threads         *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 05/31/11                   *
c     *                                                              *
c     *     Has all MPI ranks set the number of OMP threads          *
c     *     We also need to set the WARP parameters num_threads and  *
c     *     use_omp.                                                 *
c     *                                                              *
c     ****************************************************************
      subroutine wmpi_set_omp_threads(num)
            use thread_data
            implicit integer (a-z)
            include 'mpif.h'
$add common.main
c           Dummy vars
            integer :: num
c           Local variable
            integer :: mpierr, numcpy
c           Sync and set
            if (myid .eq. 0) then
                  numcpy = num
            end if
            call MPI_Bcast(numcpy,1,MPI_INTEGER
     &            ,0,MPI_COMM_WORLD,mpierr)
            call omp_set_num_threads(numcpy)
            num_threads=numcpy
      end subroutine
c
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
