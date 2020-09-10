c     ****************************************************************
c     *                                                              *
c     *                      subroutine incomp                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 9/9/2020 rhd               *
c     *                                                              *
c     *     performs the initial computations, data structure setup  *
c     *     necessary for first (time) step i                        *
c     *                                                              *
c     ****************************************************************
c
      subroutine incomp
      use global_data ! old common.main
      use main_data, only : incid, initial_stresses_user_routine
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
c            o   get initial stresses from a user-routine (if defined)
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
      if( initial_stresses_user_routine ) call initial_sig_user_routine
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
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9994 format('                 incid re-allocate failure. incomp')
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
c     ****************************************************************          
c     *                                                              *          
c     *          subroutine initial_sig_user_routine                 *          
c     *                                                              *          
c     *                  written by : rhd                            *          
c     *                                                              *          
c     *            last modified : 9/9/2020 rhd                      *          
c     *                                                              *          
c     *  Drive over all elements to call user routine to obtain      *
c     *  six components of initial stress at approximate center of   *
c     *  the element                                                 *
c     *                                                              *          
c     ****************************************************************          
c
      subroutine initial_sig_user_routine
      use global_data, only : elblks, iprops, out, nelblk
      use main_data, only : initial_stresses_user_routine,
     &                      initial_stresses_file,
     &                      initial_stresses
      use elem_block_data, only : cdest_blocks
      implicit none
c
      integer :: blk, felem, mat_type, now_thread, totdof
      integer :: ntens, ncrds, npt, layer, kspt, lrebar, local_out,
     &           noel
      double precision :: dummy_sig(6), dummy_coords(3)
      character(len=80) :: names(2)
      logical, parameter :: ldebug = .false.
      double precision, parameter :: zero = 0.0d0
      integer, external :: omp_get_thread_num
c
      if( ldebug ) then
       write(out,*) '.... inside initial_sig_user_routine ....'
       write(out,*) '          initial_stresses_user_routine: ',
     &             initial_stresses_user_routine
       write(out,*) '          initial_stresses_file: ',
     &             initial_stresses_file
      end if
c
      initial_stresses = zero  ! entire array
c
c              call the user routine with element 0 to allow for any
c              setup before the parallel region is entered.
c              there is a module setup for the storage of
c              sigini use. see code just above the example 
c              sigini routine.
c
      ntens    = 0
      ncrds    = 0
      npt      = 0
      layer    = 0
      kspt     = 0
      lrebar   = 0
      local_out = out
      noel = 0
      now_thread = 0
      dummy_coords = zero
      names(1) = " "
      names(2) = " "
c
      call sigini( dummy_sig, dummy_coords, ntens, ncrds, noel, npt,
     &             layer, kspt, lrebar, names, now_thread, local_out,
     &             initial_stresses_file )
c
c              process all element blocks. call user routine to
c              get initial stresses at element center.
c
c$OMP PARALLEL DO PRIVATE( blk, felem, mat_type, now_thread, 
c$OMP&                     totdof )
      do blk = 1, nelblk
        felem      = elblks(1,blk)
        mat_type   = iprops(25,felem)
        now_thread = omp_get_thread_num() + 1
        select case( mat_type ) ! some matl models no sig ini support
          case( 1, 3, 5, 6, 7, 8, 10 ) ! bilinear, mises, cyclic,
c                                        creep, H_2, umat, CP
            totdof = iprops(2,felem) * iprops(4,felem) 
c                    num_enodes * num_enode_dof
            call initial_sig_user_do_blk( blk, felem, now_thread,
     &                                    totdof, 
     &                                    cdest_blocks(blk)%ptr )
          case( 2 )
             write(out,9000) 'deformation'
             call die_gracefully
          case( 4 ) ! cohesive
             write(out,9000) 'cohesive'
             call die_gracefully
          case default
             write(out,9020)
             call die_gracefully
        end select
      end do
c
      return
 9000 format(//,'>>> Error: initial stresses not supported',
     &  ' for material model: ',a,
     &  /,      '           Please use equivalent eigenstrains imposed',
     &  /,      '           via temperatures in step 1.',
     &  /,      '           job terminated....', //)
 9020 format(//,'>>> Fatal Error: routine initial_sig_user_routine',
     &  /,      '                 job terminated....', //)
c
      end
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine initial_sig_user_do_blk           *          
c     *                                                              *          
c     *                  written by : rhd                            *          
c     *                                                              *          
c     *                last modified : 9/9/2020 rhd                  *          
c     *                                                              *          
c     *  Process element block to get user routine initial stresses  *
c     *                                                              *          
c     ****************************************************************          
c
      subroutine initial_sig_user_do_blk( blk, felem, now_thread,
     &                                    totdof, bcdst )
      use global_data, only : elblks, iprops, out, elelib,
     &                        mxvl, mxecor, scoords => c 
      use main_data, only : initial_stresses, initial_stresses_file
      implicit none
c
      integer, intent(in) :: blk, felem, now_thread, totdof
      integer, intent(in) :: bcdst(totdof,*)  ! totdof x span
c
c              local variables
c
      integer :: span, num_enodes, num_enode_dof, local_out,
     &           elem_type, relem, noel, node, k, j,
     &           ntens, ncrds, npt, layer, kspt, lrebar
c
      double precision :: xbar, ybar, zbar, x, y, z
      double precision :: enode_coords(mxvl,mxecor), 
     &                    e_cntr_coords(3,mxvl), eisig(6)
      double precision, parameter :: zero = 0.0d0
      character(len=80) :: names(2)
      logical, parameter :: ldebug = .false.
c
      if( ldebug ) then
        write(out,*) ' ...   blk, felem, now_thread: ',blk, felem,
     &        now_thread
      end if
c
      span           = elblks(0,blk)
      num_enodes     = iprops(2,felem)
      num_enode_dof  = iprops(4,felem)
      elem_type      = iprops(1,felem)
c
c           pull coordinates at t=0 from global input vector for all
c           elements in the block
c
      k = 1
      do j = 1, num_enodes
!DIR$ IVDEP
         do relem = 1, span
            enode_coords(relem,k)   = scoords(bcdst(k,relem))
            enode_coords(relem,k+1) = scoords(bcdst(k+1,relem))
            enode_coords(relem,k+2) = scoords(bcdst(k+2,relem))
         end do
         k = k + 3
      end do
c
c         compute coordinates of approximate center of elements
c
      do relem = 1, span
        xbar = zero; ybar = zero; zbar = zero
        do node = 1, num_enodes
          x = enode_coords(relem,node)
          y = enode_coords(relem,node+num_enodes)
          z = enode_coords(relem,node+num_enodes+num_enodes)
          xbar = xbar + x
          ybar = ybar + y
          zbar = zbar + z
       end do
       e_cntr_coords(1,relem) = xbar / dble(num_enodes)
       e_cntr_coords(2,relem) = ybar / dble(num_enodes)
       e_cntr_coords(3,relem) = zbar / dble(num_enodes)
      end do
c
      ntens    = 6
      ncrds    = 3
      npt      = 1
      layer    = 0
      kspt     = 0
      lrebar   = 0
      names(1) = " "
      names(2) = elelib(elem_type)
      local_out = out
c
c         call user initial stress routine once for all elements 
c         in block. we only request initial stresses at element center.
c
      do relem = 1, span
        noel = felem + relem - 1
        eisig = zero
        call sigini( eisig, e_cntr_coords(1,relem), 
     &               ntens, ncrds, noel, npt, layer, kspt, lrebar,
     &               names, now_thread, local_out, 
     &               initial_stresses_file )
        initial_stresses(1:4,noel) = eisig(1:4)
        initial_stresses(5,noel) = eisig(6) ! switch to WARP3D order
        initial_stresses(6,noel) = eisig(5)
      end do 
c 
      return
      end


