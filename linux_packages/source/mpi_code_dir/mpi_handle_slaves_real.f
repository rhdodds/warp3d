
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_handle_slaves           *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : rhd 1/7/2017 rhd           *
c     *                                                              *
c     *     The MPI implementation of warp3d follows a master-worker *
c     *     approach for most of the code.  The root processor       *
c     *     (processor 0) runs warp3d in full, notifying the worker  *
c     *     processors when they have work to do.                    *
c     *                                                              *
c     *     This routine serves as the waiting point for all of the  *
c     *     worker processors, where they wait for a message from    *
c     *     the root telling them what to do next.                   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_handle_slaves
      implicit none              
      include "mpif.h"
      include 'common.main'
c
      integer :: idummy1, idummy2, idumm3, idummy4, idummy5,
     &           do, ierr, idummy3, dumi, dumi2, dumi3, dum, dum2,
     &           dumj, dum1
c     
      integer, dimension(1) :: ivec_dummy1, ivec_dummy2, ivec_dummy3,
     &                         ivec_dummy4, ivec_dummy5
      integer :: adum1(max_procs), adum2(max_procs), adum3(max_procs),
     &           adum4(max_procs), num_blks(mxnmbl)
c     
      logical :: ldummy1, ldummy2, ldummy3, ldummy4, ldummy5
      double precision, dimension(1) :: dvec_dummy1, dvec_dummy2,
     &                                  dvec_dummy3, dvec_dummy4,
     &                                  dvec_dummy5
c
      logical :: duml, ldum1, ldum2, ldum3, ldum4, ldum5, ldum6,
     &           ldum7, ldum8, ldum9, ldum10
      logical, parameter :: debug = .false.
      double precision :: dumd
c
c           stay in loop, processing the master requests
c
      if( debug ) write(out,'("=> Proc",i3," in worker driver")') myid
c
      do while (.true.)
c
c               wait for a job instruction from the root processor.
c               In this case we call an MPI_BCAST, where the worker
c               processors are waiting for root to broadcast the
c               next instruction. All workers are sent same instruction
c
      if( debug ) write(out,'("====> Proc",i3," waiting for do ")') 
     &        myid
c
      call MPI_BCAST( do, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
c
      if( debug ) write(out,'("====> Proc",i3," just got ",i3)') 
     &        myid, do
c
c          now branch based on the instruction number sent.
c
      select case ( do )
c
c           do = 1: call a routine to kill the slaves (call stop).
c
      case( 1 ) 
         if(debug) write(out,'("=> Proc",i3," is exiting")') myid
         call die_gracefully
c
c           do = 2: send basic model data (coordinates, incidences, etc.)
c
      case( 2 ) 
         if(debug) write(out,'("=> Proc",i3," get basic data")') myid
         call wmpi_send_basic
c
c           do = 3: Cpardiso symmetric
c
      case( 3 ) 
         if(debug) write(out,'("=> Proc",i3," cpardiso_symmetric")')
     &           myid
         call cpardiso_symmetric( idummy1, idummy2, dvec_dummy1, 
     &        dvec_dummy2, dvec_dummy3, dvec_dummy4, ivec_dummy1, 
     &        ivec_dummy2,ldummy1, idummy3, out, myid )
c
c           do = 4: execute tanstf, the routine which generates the tangent
c                   stiffness matrix for each element. Element stiffnesses
c                   are only calculated for the elements that the processor
c                   owns.
c
      case( 4 ) 
         if(debug) write(out,'("=> Proc",i3," do tanstf")') myid
         duml = .false.
         dumi = 1
         dumj = 1
         call tanstf( duml, dumi, dumj)
c
c           do = 5: obsoleted EBE solver
c
      case( 5 ) 
         write (out,'("=> Proc",i3," obsolete lnpcg")') myid
         call die_gracefully
c
c           do = 6: call the routine which adds the lumped mass terms 
c                   into the element stiffness matricies.  Modifys only
c                   the element stiffness matricies which this processor
c                   owns.
c
      case( 6 ) 
         if(debug) write(out,'("=> Proc",i3," go to inclmass")') myid
         call inclmass
c
c           do = 7: gather the element stiffnesses to the root processor
c                   so that it can assemble a full stiffness matrix for 
c                   the direct solvers.
c
      case( 7 ) 
         if(debug) write(out,'("=> Proc",i3," gather estf")') myid
         call wmpi_combine_stf
c
c           do = 8: calculate the lumped mass matrix for all of the elements
c                   which the processor owns.
c
      case( 8 )
         if(debug) write(out,'("=> Proc",i3," go to cmpmas")') myid
         call cmpmas
c
c           do = 9: send analysis parameter data from root to all the
c                   slave processors.
c
      case( 9 ) 
         if(debug) write(out,'("=> Proc",i3," snd analysis params")') 
     &        myid
         call wmpi_send_analysis
c
c           do = 10: send constraint data from root to all the slave
c                    processors.
c
      case( 10 ) 
         if(debug) write(out,'("=> Proc",i3," send const data")') myid
         call wmpi_send_const
c
c           do = 11: remove transform. mat. due to contact
c
      case( 11 ) 
         if(debug) write(out,'("=> Proc",i3,"kill cont trnmats")') myid
         call contact_remove (.false.)
c
c           do = 12: send data needed each load step to all
c                    the slave processors.
c
      case( 12 ) 
         if(debug) write(out,'("=> Proc",i3," send step data")') myid
         call wmpi_send_step
c
c           do = 13: send info from root to the slaves about the current
c                    amount of crack growth` and whether they should
c                    kill any of their owned elements.
c
      case( 13 ) 
         if(debug) write(out,'("=> Proc",i3," get crk growth info")') 
     &        myid
         call wmpi_send_growth ( duml )
c
c           do = 14: call the routine which computes the stresses for the
c                    elements.  Each processor copmutes the stresses only
c                    for the elements it owns. 
c
      case( 14 ) 
         if(debug) write(out,'("=> Proc",i3," do stress calcs")') myid
         call drive_eps_sig_internal_forces(dumi, dumi2, duml)
c
c           do = 15: call routines that calculate strains and stress work
c                    density at element nodes for J-computations
c
      case( 15 ) 
         if(debug) write(out,'("=> Proc",i3,"do fgm setup")') myid
         call di_fgm_setup(dumi,dumi2,dumi3,duml)
c
c           do = 16: setup nodal ownership data structures for the
c                    domain decomposition version of lnpcg
c
      case( 16 )
         if(debug) write(out,'("=> Proc",i3," set up dom. dec. mpi")')
     &        myid
         call wmpi_owner_send(adum1,adum2,adum3,adum4)
c
c           do = 17: call the routine which updates the global solutions
c                    (stresses, strains, total displacements) following
c                    the successful completion of a load step.
c
      case( 17 )
         if(debug) write(out,'("=> Proc",i3," do update")') myid
         call update
c
c           do = 18: save data (stresses, strains, etc.) required for
c                    restarting a non-converged load increment using
c                    the adaptive load reduction mechanism.
c
      case( 18 )
         if(debug) write(out,'("=> Proc",i3," adapt save")') myid
         call adaptive_save
c
c           do = 19: this routine reloads the solutions for the previous
c                    load step when the current load step has failed to 
c                    converge, and we are using the adaptive mechanism
c                    to try a smaller load step.
c
      case( 19 )
         if(debug) write(out,'("=> Proc",i3," adapt reset")') myid
         call adaptive_reset
c
c           do = 20: gather the element stresses from the slave processors
c                    and send them to the root processor.  This is necessary
c                    for stress output, restarts, etc.
c
      case( 20 ) 
         if(debug) write(out,'("=> Proc",i3," gather stress")') myid
         call wmpi_get_str ( dum )
c
c           do = 21: each processor owns a complete copy of the temperatures
c                    this routine scales the temperatures for a load 
c                    step due to load multipliers or adaptive load 
c                    reduction.
c
      case( 21 ) 
         if(debug) write(out,'("=> Proc",i3," deal w/ temps")') myid
         call mnralg_scale_temps ( dum1, dum2 )
c
c           do = 22: send temperature data and element equiv nodal forces
c                    from root to the slave processors.
c                    each processor has a complete copy of the
c                    temperature data structures.
      case( 22 ) 
         if(debug) write(out,'("=> Proc",i3," send temp data")') myid
         call wmpi_send_temp_eqloads
c
c           do = 23: this routines sends the data needed by the slave
c                    processors which were read in by root after a restart. 
c
      case( 23 )
         if(debug) write(out,'("=> Proc",i3," send reopen")') myid
         call wmpi_send_reopen
c
c           do = 24: have root send data the slave processors need in the 
c                    middle of each newtons iteration (du vector, etc.).
c
      case( 24 ) 
         if(debug) write(out,'("=> Proc",i3," send du")') myid
         call wmpi_send_itern
c
c           do = 25: have the slave processors initialize the crack 
c                    growth data structures and get critical growth
c                    information from the root processor.
c
      case( 25 ) 
         if(debug) write(out,'("=> Proc",i3," init crkgrowth")') myid
         call wmpi_growth_init
c
c           do = 26: retrieve dam_ifv terms and other information
c                    owned by the slave processors but needed by
c                    root to correctly compute crack growth. 
c
      case( 26 ) 
         if(debug) write(out,'("=> Proc",i3," str crkgrowth")') myid
         call wmpi_get_grow 
c
c           do = 27: have each processor compute the parts of the
c                    J-integral domain for the elements which they
c                    own.
c
      case( 27 )
         if(debug) write(out,'("=> Proc",i3,"compute J")') myid
         call dicmj
c
c           do = 28: send each processor the contact information
c
      case( 28 )
         if(debug) write(out,'("=> Proc",i3," get contact info")') myid
         call wmpi_send_contact (duml)
c
c           do = 29: find contact for nodes referenced by each processor
c
      case( 29 )
         if(debug) write(out,'("=> Proc",i3,"find contact")') myid
         call contact_find
c
c           do = 30: allocate data structures for J calculations with
c                    temperature loadings
c
      case( 30 ) 
         if(debug) write(out,'("=> Proc",i3,"all j temp stuff")') myid
         call di_node_props_setup (dum1,dum2)
c
c           do = 31: available
c
      case( 31 )     
c
c           do = 32:  available
c
      case( 32 ) 
         write (out,'("=> Proc",i3," obsolete updpcm")') myid
c
c           do = 33: available
c
      case( 33 ) 
         if(debug) write(out,'("=> Proc",i3," allo pcm")') myid
c
c           do = 34: available
c
      case( 34 ) 
c
c
c           do = 35: patran element stress-strain file
c
      case( 35 )
         if(debug) write(out,'("=> Proc",i3," get pat ele")') myid
         call oustr_pat_flat_file( ldum1,ldum2,ldum3,ldum4,ldum5,
     &                             ldum6,ldum7,ldum8,ldum9 )
c
c
c           do = 36: call hypre routine with dummy arguments
c
      case( 36 ) 
         if (debug) write (out,'("=> Proc",i3," enter Hypre")') myid
            write(out,*) "Call no longer available.  Should not be here"
         call die_gracefully
c
c           do =37: available
c
      case( 37 ) 
         if (debug) write (out,'("=> Proc",i3," enter MKL set")')
     &                  myid
         write(out,*) "Call no longer available.  Should not be here"
c
c           do =38: available
c
      case( 38 ) 
         if (debug) write (out,'("=> Proc",i3," enter OMP set")')
     &                               myid
            write(out,*) "Call no longer available.  Should not be here"
c
c           do = 39: call distribute_from_assembled
c
      case( 39 ) 
         if (debug) write (out,'("=> Proc",i3," enter distr")') myid
         call distribute_from_assembled(0,0,0,0,0,0,0,0,0,0,0)
c
c           do = 40: enter hypre solver
c
      case( 40 ) 
         if (debug) write (out,'("=> Proc",i3," enter hypre")') myid
         call iterative_sparse_hypre(0,0,0)
c
      case( 41 ) 
         if (debug) write (out,'("=> Proc", i3,"local_sparse")') myid
         call determine_local_sparsity
c
      case( 42 ) 
         if (debug) write (out,'("=> Proc", i3,"initial_map")') myid
         call determine_initial_map
c
      case( 43 ) 
         if (debug) write(out,'("=> Proc", i3,"assem_sparse")') myid
         call assemble_sparsity
c
      case( 44 ) 
         if (debug) write(out,'("=> Proc", i3,"det_ord")') myid
         call determine_ordering
c
      case( 45 ) 
         if (debug) write(out,'("=> Proc", i3,"mov_sparse")') myid
         call move_sparsity
c    
      case( 46 ) 
         if (debug) write(out,'("=> Proc", i3,"assem_coefs")') myid
         call assemble_coefs
c
      case( 47 ) 
         if (debug) write(out,'("=> Proc", i3,"assem_load")') myid
         call assem_load_vec
c
      case( 48 ) 
         if (debug) write(out,'("=> Proc", i3,"final_set")') myid
         call dist_final_setup
c
      case( 49 ) 
         if (debug) write(out,'("=> Proc", i3,"send cryst.")') myid
         call wmpi_send_crystals
c
      case( 50 ) 
         if (debug) write(out,'("=> Proc", i3,"deal_crys.")') myid
         call wmpi_dealloc_crystals
c
      case( 51 )
         if (debug) write(out,'("=> Proc", i3,
     &       "run uexternaldb. ")') myid
         call wmpi_do_uexternaldb
c
      case( 52 ) 
         if (debug) write(out,'("=> Proc", i3,"send simple")') myid
         call wmpi_send_simple_angles
c
      case default
         write(out,9000)
         call die_abort
      end select
c
      end do !   on main do while loop
c      
      return
c
 9000 format ('>>>  FATAL ERROR:  wmpi_handle_slaves. invalid op code',
     &  /,    '                   job aborted...')
      end

