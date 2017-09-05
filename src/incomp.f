c     ****************************************************************
c     *                                                              *
c     *                      subroutine incomp                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 6/15/2013 rhd              *
c     *                                                              *
c     *     performs the initial computations, data structure setup  *
c     *     necessary for first (time) step i                        *
c     *                                                              *
c     ****************************************************************
c
      subroutine incomp
      use global_data ! old common.main
      use main_data, only : incid
      use contact, only : use_contact
      implicit integer (a-z)
      double precision
     &   zero
      integer, dimension(:), allocatable ::  temp_vec
      logical local_debug
      real, external :: wcputime
      data local_debug / .false. /
      data zero / 0.0d00 /
c
c          this code has been re-written a number of
c          times as WARP3D evolved. At present, I believe
c          this routine is called only at the start of
c          solution for time step 1 and never again.
c
c          re-size the incid vector to its required size
c          rather than the maximum possible size.
c
      if( local_debug ) write(*,*) '.... inside incomp....'
c
      allocate( temp_vec(inctop),  stat = alloc_stat )
      if ( alloc_stat .ne. 0 ) then
         write(out,9900); write(out,9994); call die_abort
      end if
      temp_vec(1:inctop) = incid(1:inctop)
      call mem_allocate( 10 )
      incid(1:inctop) = temp_vec(1:inctop)
      deallocate( temp_vec )
c
c          allocate and initialize:
c            o   build an element to blocks mapping
c                data structure for general use.
c            o   deallocate the elstor table used only
c                during input.
c            o   history dependent data for elements to
c                support material nonlinearity models.
c            o   rotation matrices at gauss points to
c                support finite strain/large displacement
c                analyses.
c            o   strains for all elements
c            o   stresses for all elements
c            o   cdest indexes for all elements.
c            o   edest indexes for all elements.
c            o   volumes for all elements
c            o   vector of rigid contact forces (always)
c            o   array of contact "causes" -- only if
c                rigid contact surfaces defined
c            o   lists of solid elements connected to
c                cohesive-interface elements
c
      call init_eblock_map
      call mem_allocate( 7 )
      call history_cep_init( 0, 1 )
      call rotation_init( 0, 1 )
      call strains_init( 0, 1 )
      call stresses_init( 0, 1 )
      call cdest_init
      call edest_init
      call element_volumes_init( 0, 1 )
      call solid_cohes_init
c
c          MPI:
c             broadcast basic model, constraint, and analysis
c             parameter information to all the worker processors.
c
      if( use_mpi ) write(out,9000) wcputime(1)
      call wmpi_send_basic
      call wmpi_send_const
      call wmpi_send_analysis
      call wmpi_init_owner
      if( use_mpi ) write(out,9005) wcputime(1)
c
c          allocate space for rigid contact info
c
      call mem_allocate( 24 )
      call mem_allocate( 25 )
c
c          compute the mass matrix at time zero. as
c          required, set flag indicating that a new
c          mass for use in the iterative procedure
c          has been computed.
c
      write(out,9150)
      call cmpmas
      newmas = .true.
c
c          generate mpc equations from input tied mesh data.
c
      call tied_mesh_driver()
c
c          uexternaldb for Abaqus compatible support
c
      douextdb = 1  ! in common.main
      call wmpi_do_uexternaldb
c
c          compute the initial internal forces if not
c          done already. we now just zero them - applied
c          initial stresses at t=0 not supported.
c
      ifv(1:nodof) = zero
c
c          set flag (init computations flag) signaling
c          that initial computations-setup have been made
c
      incflg = .true.

      if( local_debug ) write(out,*) '... leaving incomp ...'
c
c
 9000 format(/,7x,
     & '>> Start MPI toplogy setup, worker init @',f10.2)
 9005 format(/,7x,
     & '>> Done MPI toplogy setup, worker init  @',f10.2)
 9150 format(/,7x,
     & '>> computing structural diagonal mass (@ t=0)')
 9800 format('>>> FATAL ERROR: routine incomp. send to WARP3D',
     & /,    '                 developers. Job terminated.',/)
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9994 format('                 incid re-allocate failure. incomp')
 9999 return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine tied_mesh_driver                  *
c     *                                                              *
c     *                       written by : bjb                       *
c     *                                                              *
c     *                   last modified : 9/3/2017 rhd               *
c     *                                                              *
c     *     this subroutine prepares necessary temporary data        *
c     *              to call the tied mesh processor                 *
c     *                                                              *
c     ****************************************************************

      subroutine tied_mesh_driver()
      use global_data ! old common.main
      use main_data, only : crdmap, incid, incmap
      use mod_mpc, only : tied_sets_exist, tied_con_mpcs_constructed,
     &                    surface_table, tied_contact_table,
     &                    tied_con_mpc_table
      integer  err, dimnode, dimcoord, dimelem, meshformat, nnodes,
     &         icol, dumi
      integer, allocatable, dimension (:,:) :: elemcon
      real  epsilon, dumr, wcputime
      double precision, allocatable, dimension (:,:) :: nodecoord
      character(len=1) :: dums
      double precision dumd
      external wcputime
c
      if (.not. tied_sets_exist)  return
      allocate (tied_con_mpc_table(max_mpc_tied),
     &          nodecoord(nonode,3),
     &          elemcon(noelem,mxndel), stat=err)
      if (err .ne. 0) then
         call errmsg2(58,dumi,dums,dumr,dumd)
         call die_abort
      end if
c
      do i = 1, nonode
         do j = 1, 3
            nodecoord(i,j) = 0
         end do
      end do
c
      do i = 1, nonode
         nodecoord(i,1) = c(crdmap(i))
         nodecoord(i,2) = c(crdmap(i)+1)
         nodecoord(i,3) = c(crdmap(i)+2)
      end do
c
      do i = 1, noelem
         do j = 1, mxndel
            elemcon(i,j) = 0
         end do
      end do
c
      do i = 1, noelem
         nnodes = iprops(2,i)
         icol = 0
         do j = incmap(i), incmap(i)+nnodes-1
            icol = icol + 1
            elemcon(i,icol) = incid(j)
         end do
      end do
c
      dimnode    = nonode
      dimcoord   = 3
      dimelem    = noelem
      meshformat = 2
      epsilon    = 0.000001
c
      call  tied_mainprocessor(nodecoord, dimnode, dimcoord, nonode,
     &                         elemcon, dimelem, mxndel, noelem,
     &                         cstmap, meshformat, props, out, epsilon)
      write(out,9000) wcputime(1)
c
      tied_con_mpcs_constructed = .true.
      deallocate (nodecoord,elemcon,surface_table,tied_contact_table)
c
      return
c
 9000 format(/,7x,
     & '>> tied mesh processor finished  @ ', f10.2)
c
      end
