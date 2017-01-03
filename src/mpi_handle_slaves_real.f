c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_handle_slaves           *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : rhd 1/24/2015              *
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
      include 'common.main'
      logical duml, debug, ldum1, ldum2, ldum3, ldum4, ldum5, ldum6,
     &        ldum7, ldum8, ldum9, ldum10
      data debug /.false./
      double precision
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
     &         3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,
     &         4100,4200,4300,4400,4500,4600,4700,4800,4900,
     &         5000,5100,5200) do
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
c           do = 3: available. deprecated the call to lnstff.
c
 300     continue
         if(debug)write (out,'("=> Proc",i3," <available>")') myid
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
c           do = 5: obsoleted EBE solver
c
 500     continue
         write (out,'("=> Proc",i3," obsolete lnpcg")') myid
         call die_gracefully
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
c           do = 31: deprecated distribute graph structures
c                    for ebe preconditioner
c
 3100    continue
         goto 9000
c
c           do = 32: init preconditioner for lnpcg
c
 3200    continue
         write (out,'("=> Proc",i3," obsolete updpcm")') myid
         goto 9000
c
c           do = 33: no longer used
c
 3300    continue
         if(debug)write (out,'("=> Proc",i3," allo pcm")') myid
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
         call oustr_pat_flat_file( ldum1,ldum2,ldum3,ldum4,ldum5,
     &                             ldum6,ldum7,ldum8,ldum9 )
         goto 9000
c
c
c           do = 36: call hypre routine with dummy arguments
c
 3600    continue
            if (debug) write (out,'("=> Proc",i3," enter Hypre")') myid
            write(out,*) "Call no longer available.  Should not be here"
            call die_gracefully
            goto 9000
c
c           do =37: call wmpi_set_mkl_threads
c
 3700    continue
            if (debug) write (out,'("=> Proc",i3," enter MKL set")') myid
            write(out,*) "Call no longer available.  Should not be here"
c            call wmpi_set_mkl_threads(0)
            goto 9000
c
c           do =38: call wmpi_set_omp_threads
c
 3800    continue
            if (debug) write (out,'("=> Proc",i3," enter OMP set")') myid
            write(out,*) "Call no longer available.  Should not be here"
c            call wmpi_set_omp_threads(0)
            goto 9000
c
c           do = 39: call distribute_from_assembled
c
 3900    continue
            if (debug) write (out,'("=> Proc",i3," enter distr")') myid
            call distribute_from_assembled(0,0,0,0,0,0,0,0,0,0,0)
            goto 9000
c
c           do = 40: enter hypre solver
c
 4000    continue
            if (debug) write (out,'("=> Proc",i3," enter hypre")') myid
            call iterative_sparse_hypre(0,0,0)
            goto 9000
c
 4100    continue
            if (debug) write (out,'("=> Proc", i3,"local_sparse")') myid
            call determine_local_sparsity
            goto 9000
c
 4200    continue
            if (debug) write (out,'("=> Proc", i3,"initial_map")') myid
            call determine_initial_map
            goto 9000
c
 4300    continue
            if (debug) write(out,'("=> Proc", i3,"assem_sparse")') myid
            call assemble_sparsity
            goto 9000
c
 4400    continue
            if (debug) write(out,'("=> Proc", i3,"det_ord")') myid
            call determine_ordering
            goto 9000
c
 4500    continue
            if (debug) write(out,'("=> Proc", i3,"mov_sparse")') myid
            call move_sparsity
            goto 9000
c    
 4600     continue
            if (debug) write(out,'("=> Proc", i3,"assem_coefs")') myid
            call assemble_coefs
            goto 9000
c
 4700     continue
            if (debug) write(out,'("=> Proc", i3,"assem_load")') myid
            call assem_load_vec
            goto 9000
c
 4800     continue
            if (debug) write(out,'("=> Proc", i3,"final_set")') myid
            call dist_final_setup
            goto 9000
c
 4900     continue
            if (debug) write(out,'("=> Proc", i3,"send cryst.")') myid
            call wmpi_send_crystals
            goto 9000
c
 5000     continue
            if (debug) write(out,'("=> Proc", i3,"deal_crys.")') myid
            call wmpi_dealloc_crystals
            goto 9000
c
 5100     continue
            if (debug) write(out,'("=> Proc", i3,
     &       "run uexternaldb. ")') myid
            call wmpi_do_uexternaldb
            goto 9000
c
 5200     continue
            if (debug) write(out,'("=> Proc", i3,"send simple")') myid
            call wmpi_send_simple_angles
            goto 9000
c
 9000 continue
      end do
c
      return
c
      end

