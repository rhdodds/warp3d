c     ****************************************************************
c     *                                                              *
c     *                     module mpc_init                          *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 11/26/2018 rhd            *
c     *                                                              *
c     *  routine to initialize the multipoint constraint (mpc)       *
c     *  data structure variables                                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine mpc_init()
c
c           declare modules
c
      use mod_mpc
      use global_data, only : max_mpc_tied
c
c           declare variables
c
      implicit none
      integer logflag,logunit,mflag
      integer, external :: warp3d_get_device_number
      character(len=80) :: logfile,message
c
c           local variables
c
      integer :: debug_level

      logflag = 0                  ! set = 0 for no debugging log file
      logunit = warp3d_get_device_number()
      logfile = 'echo_mpc_init.log'
      message = '[replace with an error or warning message]'
      mflag = 0
c
      debug_level = 5

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'----- begin mpc_init -----'
        write(logunit,*)'         max_mpc_tied =',max_mpc_tied
        write(logunit,*)'         max_surfaces =',max_surfaces
        write(logunit,*)'        max_tied_sets =',max_tied_sets
        close(unit=logunit)
      end if
c
c           error checking
c
      if(max_mpc_tied.le.0)then
        mflag = 2
        message = 'error: max_mpc_tied = 0 (mpc_init).'
      else if(max_surfaces.le.0)then
        mflag = 2
        message = 'error: max_surfaces = 0 (mpc_init).'
      else if(max_tied_sets.le.0)then
        mflag = 2
        message = 'error: max_tied_sets = 0 (mpc_init).'
      end if
      if (mflag .gt. 0)  goto 900
c
c           initialize arrays for the mpc equations
c
      num_user_mpc     = 0
      num_tied_con_mpc = 0
      nmpc             = 0
      mpcs_exist                = .false.
      tied_sets_exist           = .false.
      tied_con_mpcs_constructed = .false.
      display_mpcs              = .false.
c
      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'         num_user_mpc =',num_user_mpc
        write(logunit,*)'     num_tied_con_mpc =',num_tied_con_mpc
        close(unit=logunit)
      end if
c
c           initialize arrays for the tied contact data
c
      num_surfaces = 0
      num_tied_sets = 0
c
      dep_trms_len = 0
      if (allocated(dep_check))     deallocate(dep_check)
      if (allocated(dep_dof))       deallocate(dep_dof)
      if (allocated(dep_ptr))       deallocate(dep_ptr)
      if (allocated(dep_trms))      deallocate(dep_trms)
      if (allocated(ind_dof))       deallocate(ind_dof)
      if (allocated(num_terms))     deallocate(num_terms)
      if (allocated(num_dep_trms))  deallocate(num_dep_trms)
      if (allocated(abs_dep_ptr))   deallocate(abs_dep_ptr)
      if (allocated(multi_list))    deallocate(multi_list)
      if (allocated(dep_coef))      deallocate(dep_coef)
      if (allocated(dep_rhs))       deallocate(dep_rhs)
c
      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'         num_surfaces =',num_surfaces
        write(logunit,*)'        num_tied_sets =',num_tied_sets
        close(unit=logunit)
      end if
c
c           jump here on an error
c
900   continue
c
      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'----- end mpc_init -----'
        close(unit=logunit)
      end if
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine tied_mainprocessor                *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 6/5/2017 rhd              *
c     *                                                              *
c     *  process the tied contact surface data and compute the       *
c     *  multipoint constraint (mpc) equations for the slave nodes;  *
c     *  the mesh surfaces can be defined by lists of elements and   *
c     *  the element faces; the slave node's location on the master  *
c     *  element face is used to generate the mpc equations.         *
c     *                                                              *
c     *  the slave (dependent) mesh surface nodes are tied           *
c     *  (constrained) to the element faces on the matching master   *
c     *  mesh surface.  to generate the constraint equation          *
c     *  coefficient values, the isoparametric coordinates of the    *
c     *  slave node on the matching master element face are needed.  *
c     *  the mpc equations are stored in the tied_con_mpc_table()    *
c     *  data type; see the mod_mpc module.                          *
c     *                                                              *
c     ****************************************************************
c
      subroutine tied_mainprocessor(
     &  nodecoord,dimnode,dimcoord,numnode,
     &  elemcon,dimelem,maxnodeperelem,numelem,
     &  cstmap,meshformat,props,out,epsilon_mpc)
c
c           declare modules
c
      use mod_mpc
      implicit none
c
c           declare variables
c
      integer
     &  dimnode,dimcoord,numnode,dimelem,maxnodeperelem,numelem,
     &  meshformat,out,
     &  elemcon(dimelem,maxnodeperelem),
     &  cstmap(*)
c
      double precision
     &  nodecoord(dimnode,dimcoord)

      real props(mxelpr,*),epsilon_mpc
c
c           variables:
c
c           nodecoord() = global x,y,z node coordinates of all the
c             nodes in the mesh
c           dimnode = row dimension for the nodecoord() array
c           dimcoord = column dimension for the nodecoord() array;
c             expecting 3 for the x,y,z dof values per node
c           numnode = the given number of nodes in the mesh from the
c             fea input data
c           elemcon() = array of the element connectivity node list for
c             each element in the mesh from the fea input file
c           dimelem = row dimension for the elemcon() array; the
c             maximum number of elements
c           maxnodeperelem = column dimension for the elemcon() array;
c             the maximum number of nodes per element (usually 20 or 27
c             nodes per element for a higher order brick element)
c           numelem = the given number of elements in the mesh from the
c             fea input file
c           cstmap() = given constraint flag for each dof,
c             0=unconstrained and active dof
c           meshformat = format of the fea input data, affects the
c             element connectivity ordering and element face order;
c             1=abaqus input data; 2=warp3d input data
c           props() = element property array, use to get values like
c             the element type (brick or tetrahedron), and the number
c             of nodes per element
c           out = the output file unit number
c           epsilon_mpc = small value to avoid saving very small
c             mpc coefficient multiplier values, only save the mpc
c             coefficient value if it is larger than epsilon_mpc
c
c
c           local variables:
c
c           tiednodeid() = list the slave node id in column 1, list the
c             master element id that the slave node is connected to in
c             column 2, list the master element face number in column 3,
c             give the node tied flag 1=yes, 0=no in column 4, if the
c             node is within the gap tolerance distance then set the
c             flag = 1 = yes to be included in the tied contact
c           tiednodeisocoord()  = element local isoparametric
c             coordinates (xi,eta,zeta) of the slave nodes on the
c             master element face; the local isoparametric coordinates
c             are used to compute the constraint equation coefficient
c             values
c           tiednodeglobalcoord() = x,y,z global coordinates of the
c             slave nodes; the slave node coordinates can be updated to
c             be located directly on the master element face
c           tiednodeadjustflag() = true/false flag to adjust the initial
c             location of the slave node onto the master element face;
c             true = adjust the slave node's initial position
c           tiednodegaptolerance() = surface gap tolerance distance
c             between the master and slave surfaces; the slave node must
c             be within this distance to be considered for tied contact
c             with the master surface
c           maxtiednode = array row dimension; the maximum number of
c             tied slave nodes
c           maxtieddatacol = array column dimension for the
c             tiednodeid() array
c           maxisodof = array column dimension; expecting 3 for the
c             xi,eta,zeta isoparametric coordinate dof values
c           numtiednode = return the number of tied slave nodes
c           logflag,logunit,logfile = log file for debugging and testing;
c             the flag, file unit, and file name for the log file
c           message = error or warning message string returned from a
c             subroutine, check the value of mflag
c           mflag = 0=no message, 1=warning, 2=error; check the message
c             text string
c
c           note: see also the mpc module for the module variables
c
c
c           declare local variables
c
      integer
     &  maxtiednode,maxtieddatacol,maxisodof,numtiednode,
     &  allocate_status,slave_surface_count,
     &  tiedset,numsurfpairs,row,surfaceid,
     &  logflag,logunit,mflag
      integer, parameter :: numchar = 512
      integer, external :: warp3d_get_device_number
      logical
     &  export_mpc_abaqus
      character(len=numchar)
     &  logfile,message
      integer, allocatable :: tiednodeid(:,:)
      logical, allocatable :: tiednodeadjustflag(:)

      double precision
     &  zero,epsilon,tiednodeisocoord,tiednodeglobalcoord
      allocatable
     &  tiednodeisocoord(:,:),tiednodeglobalcoord(:,:)
      real, allocatable ::    tiednodegaptolerance(:)
      data zero
     &   / 0.0d00 /
c
c           set local debugging log file values
c
c           set logflag = 0 for no debugging information written to
c           the log file; set logflag from 1 (minimum debug) to 10
c           (maximum debug) for debugging information to be written
c           to the log file
c
      logflag = 0                  ! set = 0 for no debugging log file
      logunit = warp3d_get_device_number()
      logfile = 'echo_tiednodedebug.log'
      message = '[replace with an error or warning message]'
      mflag = 0

      if(logflag.ge.1)then
        open(unit=logunit,file=logfile,status='unknown')
        write(logunit,*)'----- log file tied_mainprocessor -----'
        write(logunit,*)'           meshformat =',meshformat
        write(logunit,*)'              dimnode =',dimnode
        write(logunit,*)'             dimcoord =',dimcoord
        write(logunit,*)'              numnode =',numnode
        write(logunit,*)'              dimelem =',dimelem
        write(logunit,*)'       maxnodeperelem =',maxnodeperelem
        write(logunit,*)'              numelem =',numelem
        write(logunit,*)'        file unit out =',out
        write(logunit,*)'          epsilon_mpc =',epsilon_mpc
        write(logunit,*)'              logflag =',logflag
        close(unit=logunit)
      end if
c
c           error checking
c
      if(numnode.le.0)then
        mflag = 2
        message = 'error: numnode = 0 (tied_mainprocessor).'
      else if(numnode.gt.dimnode)then
        mflag = 2
        message = 'error: numnode > dimnode (tied_mainprocessor).'
      else if(dimcoord.le.0)then
        mflag = 2
        message = 'error: dimcoord = 0 (tied_mainprocessor).'
      else if(numelem.le.0)then
        mflag = 2
        message = 'error: numelem = 0 (tied_mainprocessor).'
      else if(numelem.gt.dimelem)then
        mflag = 2
        message = 'error: numelem > dimelem (tied_mainprocessor).'
      else if(maxnodeperelem.le.0)then
        mflag = 2
        message = 'error: maxnodeperelem =0 (tied_mainprocessor).'
      else if(num_surfaces.le.0)then
        mflag = 2
        message = 'error: num_surfaces = 0 (tied_mainprocessor).'
      else if(epsilon_mpc.lt.zero)then
        mflag = 2
        message = 'error: the given epsilon_mpc < 0.0; expecting a'//
     &    ' value for epsilon_mpc >= 0.0 (tied_mainprocessor).'
      end if
      if(mflag.ge.2)goto 900
c
c           initialize values
c
      slave_surface_count = 0
c
c           count the number of elements in the slave surface lists;
c           the product of the number of slave elements with the
c           maximum number of nodes per element gives an estimate for
c           the maximum number of slave nodes and the array dimension
c
      do tiedset=1,num_tied_sets
        numsurfpairs = tied_contact_table(tiedset)%num_pairs

        if(numsurfpairs.gt.0)then
          do row=1,numsurfpairs
            surfaceid = tied_contact_table(tiedset)%slave_list(row)

            if(surfaceid.gt.0 .and. surfaceid.le.num_surfaces)then
              slave_surface_count = slave_surface_count +
     &                              surface_table(surfaceid)%num_elems
            else
              mflag = 2
              message = 'error: invalid mesh surface id number;'//
     &          ' check the tied_contact_table()%slave_list()'//
     &          ' (tied_mainprocessor).'
              goto 900
            end if

          end do
        end if

      end do

      if(slave_surface_count.gt.0)then
        maxtiednode = maxnodeperelem*slave_surface_count
      else
        maxtiednode = numnode            ! default array dimension
        mflag = 1
        message = 'warning: zero slave elements counted; setting the'//
     &            ' local array dimension maxtiednode to default'//
     &            ' size (tied_mainprocessor).'
        if(logflag.ge.1)
     &    call srt_writelogmessage(mflag,message,logfile,logunit)
      end if

      maxtieddatacol = 4
      maxisodof = 3
      numtiednode = 0
c
c           allocate local arrays and initialize their value
c
      allocate(tiednodeid(maxtiednode,maxtieddatacol),
     &  stat=allocate_status)
      if(allocate_status.gt.0)then
        message = 'memory error during allocation of the local'//
     &            ' tiednodeid array.'
        goto 900
      end if

      allocate(tiednodeisocoord(maxtiednode,maxisodof),
     &  stat=allocate_status)
      if(allocate_status.gt.0)then
        message = 'memory error during allocation of the local'//
     &            ' tiednodeisocoord array.'
        goto 900
      end if

      allocate(tiednodeglobalcoord(maxtiednode,dimcoord),
     &  stat=allocate_status)
      if(allocate_status.gt.0)then
        message = 'memory error during allocation of the local'//
     &            ' tiednodeglobalcoord array.'
        goto 900
      end if

      allocate(tiednodeadjustflag(maxtiednode),stat=allocate_status)
      if(allocate_status.gt.0)then
        message = 'memory error during allocation of the local'//
     &            ' tiednodeadjustflag array.'
        goto 900
      end if

      allocate(tiednodegaptolerance(maxtiednode),stat=allocate_status)
      if(allocate_status.gt.0)then
        message = 'memory error during allocation of the local'//
     &            ' tiednodegaptolerance array.'
        goto 900
      end if

      tiednodeid = 0
      tiednodeisocoord = zero
      tiednodeglobalcoord = zero
      tiednodeadjustflag = .false.
      tiednodegaptolerance = 0.0 ! single

      if(logflag.ge.1)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'initialize local values'
        write(logunit,*)'         max_surfaces =',max_surfaces
        write(logunit,*)'         num_surfaces =',num_surfaces
        write(logunit,*)'  slave_surface_count =',slave_surface_count
        write(logunit,*)'          maxtiednode =',maxtiednode
        write(logunit,*)'       maxtieddatacol =',maxtieddatacol
        write(logunit,*)'            maxisodof =',maxisodof
        close(unit=logunit)
      end if
c
c           compute the location of the slave nodes in isoparametric
c           coordinates on the master element faces
c
      call tied_nodelocate(tiednodeid,tiednodeisocoord,
     &  tiednodeglobalcoord,tiednodeadjustflag,tiednodegaptolerance,
     &  maxtiednode,maxtieddatacol,maxisodof,numtiednode,
     &  nodecoord,dimnode,dimcoord,numnode,
     &  elemcon,dimelem,maxnodeperelem,numelem,
     &  meshformat,props,out,
     &  logflag,logunit,logfile,mflag,message)
      if(mflag.ge.2)goto 900
c
c           compute the tied contact constraint equation coefficient
c           values for the slave nodes to tie them to the master
c           element node degrees of freedom
c
      call tied_nodecompute(tiednodeid,tiednodeisocoord,
     &  maxtiednode,maxtieddatacol,maxisodof,numtiednode,
     &  elemcon,dimelem,maxnodeperelem,numelem,
     &  numnode,dimcoord,meshformat,props,cstmap,epsilon_mpc,
     &  logflag,logunit,logfile,mflag,message)
      if(mflag.ge.2)goto 900
c
c           export the generated mpc equations to an output file
c           in warp3d user defined mpc format or in abaqus *equation
c           format
c
c           export the slave node and assigned master element to an
c           output file (6/9/04 gvt)
c
c
      export_mpc_abaqus = .false.

      if(display_mpcs)then
        call tied_exportmpc_nodelist(tiednodeid,
     &         maxtiednode,maxtieddatacol,numtiednode,out)

        epsilon = epsilon_mpc
        call tied_exportmpc_warp3dformat(epsilon,out)
      end if

      if(export_mpc_abaqus)then
        epsilon = epsilon_mpc
        call tied_exportmpc_abaqusformat(epsilon,out)
      end if
c----------------------------------------------------------------------
c           jump here on an error for exit procedures
c
900   continue
c
c           deallocate local arrays
c
      if(allocated(tiednodeid))deallocate(tiednodeid)
      if(allocated(tiednodeadjustflag))deallocate(tiednodeadjustflag)
      if(allocated(tiednodeisocoord))deallocate(tiednodeisocoord)
      if(allocated(tiednodeglobalcoord))deallocate(tiednodeglobalcoord)
      if(allocated(tiednodegaptolerance))
     &    deallocate(tiednodegaptolerance)
c
c           report error or warning messages
c
      if(logflag.ge.1)then
        open(unit=logunit,file=logfile,position='append')
        if(mflag.ge.1)then
          write(logunit,*)'    mflag =',mflag
          write(logunit,'(a)')trim(message)
        end if
        write(logunit,*)'----- end tied_mainprocessor -----'
        close(unit=logunit)
      end if

      if(mflag.ge.2)then
        write(out,9001) mflag,trim(message)
        call die_gracefully
        stop
      end if

9001  format(/,
     &1x,'>> fatal error: job aborted during tied contact
     &processing.',/,
     &1x,'>> error flag status: mflag =',i2,/,
     &1x,'>> additional error message = ',/,
     &1x,'>> ',a,/)

      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine tied_nodelist                   *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 11/01/02                  *
c     *                                    02/26/03 gvt              *
c     *                                                              *
c     *  build a list of the tied slave node id numbers from the     *
c     *  element faces listed on the slave mesh surface; use the     *
c     *  element connectivity to get the node id numbers on the      *
c     *  slave mesh surface element faces; also refer to             *
c     *  subroutine locatetiednode.f                                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine tied_nodelist(
     &  tiednodeid,tiednodeadjustflag,tiednodegaptolerance,
     &  maxtiednode,maxtieddatacol,numtiednode,
     &  elemcon,dimelem,maxnodeperelem,numelem,
     &  meshformat,props,
     &  logflag,logunit,logfile,mflag,message)
c
c           declare modules
c
      use mod_mpc
      implicit none
c
c           declare variables
c
      integer
     &  maxtiednode,maxtieddatacol,numtiednode,
     &  dimelem,maxnodeperelem,numelem,
     &  meshformat,
     &  logflag,logunit,mflag,
     &  tiednodeid(maxtiednode,maxtieddatacol),
     &  elemcon(dimelem,maxnodeperelem)

      real
     &  tiednodegaptolerance(maxtiednode),
     &  props(mxelpr,*)

      logical
     &  tiednodeadjustflag(maxtiednode)

      character(len=*) ::
     &  logfile,message
c
c           local variables
c
      integer, parameter :: slave_surface_flag = 2
      integer, parameter :: maxelemface = 6
      integer, parameter :: maxnodeperface = 9
      integer
     &  debug_level,allocate_error,
     &  surfaceid,row,element_id,face_id,node_per_elem,numnodeperface,
     &  maxnodelist,numnodelist,maxelempersurface,maxslavesurface,
     &  numelemlist,numelemface,nd,node_id,elemchoice,
     &  numcornerperface,tiedset,numsurfpairs,iword,
     &  localfacenode(maxelemface,maxnodeperface)
      integer, allocatable :: nodelist(:)
      integer, allocatable :: surfaceflagtable(:)
      logical, allocatable :: surfaceadjustgap(:)
      logical, allocatable :: nodeadjustgap(:)
      real, allocatable :: surfacegaptolerance(:)
      real, allocatable :: nodegaptolerance(:)
      real rword
      logical
     &  local_debug
      equivalence (rword,iword)
c
c           debugging information
c
      debug_level = 4
      local_debug = .false.

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'----- begin tied_nodelist -----'
        write(logunit,*)'         num_surfaces =',num_surfaces
        write(logunit,*)'         max_surfaces =',max_surfaces
        write(logunit,*)'        num_tied_sets =',num_tied_sets
        write(logunit,*)'        max_tied_sets =',max_tied_sets
        write(logunit,*)'          maxtiednode =',maxtiednode
        write(logunit,*)'       maxtieddatacol =',maxtieddatacol
        write(logunit,*)'              dimelem =',dimelem
        write(logunit,*)'       maxnodeperelem =',maxnodeperelem
        write(logunit,*)'              numelem =',numelem
        write(logunit,*)'           meshformat =',meshformat
        if(local_debug)then
          write(logunit,*)'          local_debug = ',local_debug
        end if
        close(unit=logunit)
      end if
c
c           error checking
c
      if(numelem.le.0)then
        mflag = 2
        message = 'error: numelem = 0 (tied_nodelist).'
      else if(numelem.gt.dimelem)then
        mflag = 2
        message = 'error: numelem > dimelem (tied_nodelist).'
      else if(maxtiednode.le.0)then
        mflag = 2
        message = 'error: maxtiednode = 0 (tied_nodelist).'
      else if(maxtieddatacol.le.0)then
        mflag = 2
        message = 'error: maxtieddatacol = 0 (tied_nodelist).'
      else if(num_surfaces.le.0)then
        mflag = 2
        message = 'error: num_surfaces = 0 (tied_nodelist).'
      else if(num_tied_sets.le.0)then
        mflag = 2
        message = 'error: num_tied_sets = 0 (tied_nodelist).'
      end if
      if(mflag.ge.2)goto 900
c
c           initialize values
c
c           get the maximum number of elements for a mesh surface;
c           note that the number of elements per mesh surface is
c           stored in the surface_table()%num_elems data structure;
c           use the maximum number of elements to allocate the local
c           nodelist() array;
c
      maxelempersurface = 0
      maxslavesurface = 0

      do surfaceid=1,num_surfaces
        if(surface_table(surfaceid)%num_elems.gt.maxelempersurface)then
          maxelempersurface = surface_table(surfaceid)%num_elems
        end if
      end do

      do tiedset=1,num_tied_sets
        numsurfpairs = tied_contact_table(tiedset)%num_pairs
        if(numsurfpairs.gt.0)then
          maxslavesurface = maxslavesurface + numsurfpairs
        end if
      end do

      numnodelist = 0
      maxnodelist = maxelempersurface*maxslavesurface*maxnodeperface

      if(maxnodelist.le.0)then
        mflag = 2
        message = 'internal error: invalid local array dimension;'//
     &    ' maxnodelist = 0; check the local "nodelist()" array'//
     &    ' (tied_nodelist).'
        goto 900
      end if

      allocate(nodelist(maxnodelist),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "nodelist()" array'//
     &    ' (tied_nodelist).'
        goto 900
      end if
      nodelist = 0

      allocate(nodeadjustgap(maxnodelist),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "nodeadjustgap()"'//
     &    ' array (tied_nodelist).'
        goto 900
      end if
      nodeadjustgap = .false.

      allocate(nodegaptolerance(maxnodelist),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "nodegaptolerance()"'//
     &    ' array (tied_nodelist).'
        goto 900
      end if
      nodegaptolerance = 0
c
c           set a flag for the slave mesh surfaces as identified in
c           the tied_contact_table()%slave_list() for each contact
c           pair; also save the adjust gap flag and surface gap
c           tolerance distance from the tied contact data to the
c           mesh surface list;
c           first allocate local arrays for each surface
c
      allocate(surfaceflagtable(num_surfaces),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "surfaceflagtable()"'//
     &    ' array (tied_nodelist).'
        goto 900
      end if
      surfaceflagtable = 0

      allocate(surfaceadjustgap(num_surfaces),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "surfaceadjustgap()"'//
     &    ' array (tied_nodelist).'
        goto 900
      end if
      surfaceadjustgap = .false.

      allocate(surfacegaptolerance(num_surfaces),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &   ' error status > 0; check the local "surfacegaptolerance()"'//
     &   ' array (tied_nodelist).'
        goto 900
      end if
      surfacegaptolerance = 0.0

      do tiedset=1,num_tied_sets
        numsurfpairs = tied_contact_table(tiedset)%num_pairs
        if(numsurfpairs.gt.0)then
          do row=1,numsurfpairs
            surfaceid = tied_contact_table(tiedset)%slave_list(row)
            if(surfaceid.gt.0 .and. surfaceid.le.num_surfaces)then
              surfaceflagtable(surfaceid) = slave_surface_flag
              surfaceadjustgap(surfaceid) =
     &            tied_contact_table(tiedset)%adjust_gap
              surfacegaptolerance(surfaceid) =
     &            tied_contact_table(tiedset)%tolerance
            else
              mflag = 2
              message = 'error: invalid mesh surface id number;'//
     &          ' check the tied_contact_table()%slave_list()'//
     &          ' (tied_nodelist).'
              goto 900
            end if
          end do
        end if
      end do

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'local values:'
        write(logunit,*)'    maxelempersurface =',maxelempersurface
        write(logunit,*)'      maxslavesurface =',maxslavesurface
        write(logunit,*)'          maxnodelist =',maxnodelist
        write(logunit,1)'surfaceid','surfaceflagtable',
     &    'surfaceadjustgap','surfacegaptolerance'
        do surfaceid=1,num_surfaces
          write(logunit,2)surfaceid,surfaceflagtable(surfaceid),
     &      surfaceadjustgap(surfaceid),surfacegaptolerance(surfaceid)
        end do
        close(unit=logunit)
      end if

1     format(a10,3a20)
2     format(i10,i20,l20,f20.8)
c
c           build the list of slave nodes from the given mesh surface
c           defined by element lists.
c
c           check the surfaceflagtable() for slave mesh surfaces with
c           flag = 2; get the slave surface node id numbers from the
c           element faces on the mesh surface
c
      do surfaceid=1,num_surfaces
        if(surfaceflagtable(surfaceid).eq.slave_surface_flag)then
          if(logflag.ge.debug_level)then
            open(unit=logunit,file=logfile,position='append')
            write(logunit,*)
            write(logunit,*)'               mesh surface =',surfaceid
            write(logunit,*)'  surface_table()%num_elems =',
     &                         surface_table(surfaceid)%num_elems
            write(logunit,*)'         surface_table()%id = ',
     &                         trim(surface_table(surfaceid)%id)
            write(logunit,*)'         surfaceflagtable() =',
     &                         surfaceflagtable(surfaceid)
            write(logunit,*)
            close(unit=logunit)
          end if

          numelemlist = surface_table(surfaceid)%num_elems
          if(numelemlist.le.0)then
            mflag = 2
            message = 'error: numelemlist = 0 for a slave mesh'//
     &                ' surface; check surface_table()%num_elems'//
     &                ' (tied_nodelist).'
          end if
          if(mflag.ge.2)goto 900

          do row=1,numelemlist
c
c               get the element id number, the element face number,
c               and the number of nodes per element for this
c               slave mesh surface element; use the nodes per element
c               to identify the element type
c
            element_id = surface_table(surfaceid)%elem_list(row)
            face_id = surface_table(surfaceid)%face_list(row)

            if(element_id.le.0)then
              mflag = 2
              message = 'error: element_id = 0 for a slave mesh'//
     &                  ' element; check surface_table()%elem_list()'//
     &                  ' (tied_nodelist).'
            else if(element_id.gt.numelem)then
              mflag = 2
              message = 'error: element_id > numelem for a slave'//
     &                  ' mesh element; check column 1 of the'//
     &                  ' surfaceelementlist() array (tied_nodelist).'
            else if(face_id.le.0)then
              mflag = 2
              message = 'error: face_id = 0 for a slave mesh'//
     &                  ' element; check surface_table()%elem_list()'//
     &                  ' (tied_nodelist).'
            end if
            if(mflag.ge.2)goto 900

            rword = props(1,element_id)
            elemchoice = iword

            rword = props(2,element_id)
            node_per_elem = iword
c
c               get the local element face node numbers (the columns
c               in the connectivity array) for this element face
c
            select case(meshformat)
            case(1)      ! abaqus
              call srt_setlocalfacenodes(localfacenode,maxelemface,
     &               maxnodeperface,numelemface,numnodeperface,
     &               numcornerperface,node_per_elem,elemchoice,
     &               logflag,logunit,logfile,mflag,message)
              if(mflag.ge.2)goto 900
            case(2)      ! warp3d
              call srt_setwarp3dlocalfacenodes(localfacenode,
     &               maxelemface,maxnodeperface,numelemface,
     &               numnodeperface,numcornerperface,node_per_elem,
     &               elemchoice,
     &               logflag,logunit,logfile,mflag,message)
              if(mflag.ge.2)goto 900
            case default
              mflag = 2
              message = 'error: unexpected meshformat choice;'//
     &                  ' cannot get local element face nodes'//
     &                  ' (tied_nodelist).'
            end select      ! meshformat
            if(mflag.ge.2)goto 900
c
c               save the element face nodes on the slave mesh surface
c               to the local nodelist() array
c
            if(face_id.gt.numelemface)then
              mflag = 2
              message = 'error: invalid face_id value;'//
     &                  ' face_id > numelemface; check the surface'//
     &                  ' data (tied_nodelist).'
            end if
            if(mflag.ge.2)goto 900

            do nd=1,numnodeperface
              node_id = elemcon(element_id,localfacenode(face_id,nd))
              if(node_id.le.0)then
                mflag = 1
                message = 'warning: element face node = 0;'//
     &                    ' skipping the slave surface face node'//
     &                    ' (tied_nodelist).'
                if(logflag.ge.1)
     &            call srt_writelogmessage(mflag,message,logfile,
     &                                     logunit)
              else      ! save the node_id to the local nodelist() array
                numnodelist = numnodelist + 1
                if(numnodelist.gt.maxnodelist)then
                  mflag = 2
                  message = 'error: numnodelist > maxnodelist;'//
     &                  ' increase the local array dimension'//
     &                  ' maxnodelist (tied_nodelist).'
                  goto 900
                end if
c
c                   save the slave node id number;
c                   also save the the adjust gap flag and gap tolerance
c                   distance from the contact pair to the node
c
                nodelist(numnodelist) = node_id
                nodeadjustgap(numnodelist) =
     &              surfaceadjustgap(surfaceid)
                nodegaptolerance(numnodelist) =
     &              surfacegaptolerance(surfaceid)
              end if
            end do
          end do
        end if
      end do
c
c               check for duplicate node id numbers in the local
c               nodelist() array; remove the duplicate id numbers
c
      call srt_checkduplicatemultiarray(nodelist,nodeadjustgap,
     &        nodegaptolerance,maxnodelist,numnodelist,
     &        mflag,message,logflag,logunit,logfile)
      if(mflag.ge.2)goto 900
c
c               sort the node list into ascending order
c
      call srt_sortintmultiarray(nodelist,maxnodelist,1,numnodelist,1,
     &             nodeadjustgap,nodegaptolerance,1)
c
c               save the slave node id numbers from the nodelist() array
c               to column 1 in the tiednodeid() array;
c               save the gap adjust flag and surface gap tolerance
c               distance for each node to the tied node arrays
c
      do nd=1,numnodelist
        numtiednode = numtiednode + 1
        if(numtiednode.gt.maxtiednode)then
          mflag = 2
          message = 'error: too many nodes in the tiednodeid() array'//
     &              ' numtiednode > maxtiednode;'//
     &              ' increase the maxtiednode array dimension'//
     &              ' (tied_nodelist).'
          goto 900
        end if
        tiednodeid(numtiednode,1) = nodelist(nd)
        tiednodeadjustflag(numtiednode) = nodeadjustgap(nd)
        tiednodegaptolerance(numtiednode) = nodegaptolerance(nd)
      end do

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'      numnodelist =',numnodelist
        write(logunit,*)'      numtiednode =',numtiednode
        write(logunit,*)'      maxtiednode =',maxtiednode
        write(logunit,*)' tiednodeid() array:'
        write(logunit,3)'row','node_id','adjust','gap_tolerance'
        do nd=1,numtiednode
          write(logunit,4)nd,tiednodeid(nd,1),nodeadjustgap(nd),
     &      nodegaptolerance(nd)
        end do
        close(unit=logunit)
      end if

3     format(a8,a10,a8,a16)
4     format(i8,i10,l8,es16.6)
c
c----------------------------------------------------------------------
c           jump here on an error for exit procedures
c
900   continue
c
c           deallocate local arrays
c
      if(allocated(nodelist))deallocate(nodelist)
      if(allocated(surfaceflagtable))deallocate(surfaceflagtable)
      if(allocated(surfaceadjustgap))deallocate(surfaceadjustgap)
      if(allocated(nodeadjustgap))deallocate(nodeadjustgap)
      if(allocated(surfacegaptolerance))deallocate(surfacegaptolerance)
      if(allocated(nodegaptolerance))deallocate(nodegaptolerance)

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'----- end tied_nodelist -----'
        close(unit=logunit)
      end if

      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine tied_nodelocate                  *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 11/01/02                  *
c     *                                    02/25/03 gvt              *
c     *                                                              *
c     *  compute the location of the slave nodes (the dependent,     *
c     *  constrained nodes) in isoparametric coordinates on the      *
c     *  master element faces (the independent nodes); assign each   *
c     *  slave node to a master element face; an option is available *
c     *  to update the slave node global x,y,z coordinates to lie    *
c     *  exactly on the master element face.                         *
c     *                                                              *
c     ****************************************************************
c
      subroutine tied_nodelocate(
     &  tiednodeid,tiednodeisocoord,tiednodeglobalcoord,
     &  tiednodeadjustflag,tiednodegaptolerance,
     &  maxtiednode,maxtieddatacol,maxisodof,numtiednode,
     &  nodecoord,dimnode,dimcoord,numnode,
     &  elemcon,dimelem,maxnodeperelem,numelem,
     &  meshformat,props,out,
     &  logflag,logunit,logfile,mflag,message)
c
c           declare modules
c
      use mod_mpc
      implicit none
c
c           declare variables
c
      integer
     &  maxtiednode,maxtieddatacol,maxisodof,numtiednode,
     &  dimnode,dimcoord,numnode,dimelem,maxnodeperelem,numelem,
     &  meshformat,out,
     &  logflag,logunit,mflag,
     &  tiednodeid(maxtiednode,maxtieddatacol),
     &  elemcon(dimelem,maxnodeperelem)

      double precision
     &  tiednodeisocoord(maxtiednode,maxisodof),
     &  tiednodeglobalcoord(maxtiednode,dimcoord),
     &  nodecoord(dimnode,dimcoord)

      real
     &  tiednodegaptolerance(maxtiednode),
     &  props(mxelpr,*)

      logical
     &  tiednodeadjustflag(maxtiednode)

      character(len=*) ::
     &  logfile,message
c
c           variables:
c
c           tiednodeid() = list the slave node id in column 1, list the
c             master element id that the slave node is connected to in
c             column 2, list the master element face number in column 3,
c             give the node tied flag 1=yes, 0=no in column 4, if the
c             node is within the gap tolerance distance then set the
c             flag = 1 = yes to be included in the tied contact
c           tiednodeisocoord()  = element local isoparametric
c             coordinates (xi,eta,zeta) of the slave nodes on the master
c             element face; the local isoparametric coordinates are used
c             to compute the constraint equation coefficient values
c           tiednodeglobalcoord() = x,y,z global coordinates of the slave
c             nodes; the slave node coordinates can be updated to be
c             located directly on the master element face
c           tiednodeadjustflag() = true/false flag to adjust the initial
c             location of the slave node onto the master element face;
c             true = adjust the slave node's initial position
c           tiednodegaptolerance() = surface gap tolerance distance
c             between the master and slave surfaces; the slave node must
c             be within this distance to be considered for tied contact
c             with the master surface
c           maxtiednode = array row dimension; the maximum number of
c             tied slave nodes
c           maxtieddatacol = array column dimension
c           maxisodof = array column dimension; expecting 3 for the
c             xi,eta,zeta isoparametric coordinate dof values
c           numtiednode = return the number of tied slave nodes
c           nodecoord() = global x,y,z node coordinates of all the nodes
c             in the mesh
c           dimnode = row dimension for the nodecoord() array
c           dimcoord = column dimension for the nodecoord() array;
c             expecting 3 for the x,y,z dof values per node
c           numnode = the given number of nodes in the mesh from the fea
c             input file
c           elemcon() = array of the element connectivity node list for
c             each element in the mesh from the fea input file
c           dimelem = row dimension for the elemcon() array; the maximum
c             number of elements
c           maxnodeperelem = column dimension for the elemcon() array;
c             the maximum number of nodes per element (usually 20 or 27
c             nodes per element for a higher order brick element)
c           numelem = the given number of elements in the mesh from the
c             fea input file
c           meshformat = format of the fea input file;
c             1=abaqus input data file; 2=warp3d input data file
c
c           note: see also the mpc module for the module variables
c
c           logflag,logunit,logfile = log file for debugging and testing;
c             the flag, file unit, and file name for the log file
c           mflag,message = error or warning message; integer flag
c             0=no message, 1=warning, 2=error; message text string
c
c
c           local variables
c
      integer debug_level,row
c
      debug_level = 2

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'----- begin tied_nodelocate -----'
        write(logunit,*)'         num_surfaces =',num_surfaces
        write(logunit,*)'          maxtiednode =',maxtiednode
        write(logunit,*)'       maxtieddatacol =',maxtieddatacol
        write(logunit,*)'            maxisodof =',maxisodof
        write(logunit,*)'              dimnode =',dimnode
        write(logunit,*)'             dimcoord =',dimcoord
        write(logunit,*)'              numnode =',numnode
        write(logunit,*)'              dimelem =',dimelem
        write(logunit,*)'       maxnodeperelem =',maxnodeperelem
        write(logunit,*)'              numelem =',numelem
        write(logunit,*)'           meshformat =',meshformat
        close(unit=logunit)
      end if
c
c           error checking
c
      if(numnode.le.0)then
        mflag = 2
        message = 'error: numnode = 0 (tied_nodelocate).'
      else if(numnode.gt.dimnode)then
        mflag = 2
        message = 'error: numnode > dimnode (tied_nodelocate).'
      else if(dimcoord.le.0)then
        mflag = 2
        message = 'error: dimcoord = 0 (tied_nodelocate).'
      else if(numelem.le.0)then
        mflag = 2
        message = 'error: numelem = 0 (tied_nodelocate).'
      else if(numelem.gt.dimelem)then
        mflag = 2
        message = 'error: numelem > dimelem (tied_nodelocate).'
      else if(num_surfaces.le.0)then
        mflag = 2
        message = 'error: num_surfaces = 0 (tied_nodelocate).'
      else if(maxtiednode.le.0)then
        mflag = 2
        message = 'error: maxtiednode = 0 (tied_nodelocate).'
      else if(maxtieddatacol.le.0)then
        mflag = 2
        message = 'error: maxtieddatacol = 0 (tied_nodelocate).'
      else if(maxisodof.le.0)then
        mflag = 2
        message = 'error: maxisodof = 0 (tied_nodelocate).'
      end if
      if(mflag.ge.2)goto 900
c
c           build a list of the tied slave node id numbers from the
c           element faces listed on the slave mesh surface
c
      call tied_nodelist(tiednodeid,tiednodeadjustflag,
     &  tiednodegaptolerance,maxtiednode,maxtieddatacol,numtiednode,
     &  elemcon,dimelem,maxnodeperelem,numelem,
     &  meshformat,props,
     &  logflag,logunit,logfile,mflag,message)
      if(mflag.ge.2)goto 900
c
c           assign each slave node to a master element face
c
      call tied_nodeassign(tiednodeid,tiednodeadjustflag,
     &  tiednodegaptolerance,
     &  maxtiednode,maxtieddatacol,numtiednode,
     &  elemcon,dimelem,maxnodeperelem,numelem,
     &  nodecoord,dimnode,dimcoord,numnode,meshformat,props,out,
     &  logflag,logunit,logfile,mflag,message)
      if(mflag.ge.2)goto 900
c
c           compute the isoparametric coordinates for each slave node
c           on its assigned master element face
c
      call tied_nodeiso(tiednodeid,tiednodeisocoord,
     &  tiednodeglobalcoord,tiednodeadjustflag,tiednodegaptolerance,
     &  maxtiednode,maxtieddatacol,maxisodof,
     &  numtiednode,elemcon,dimelem,maxnodeperelem,numelem,
     &  nodecoord,dimnode,dimcoord,numnode,
     &  meshformat,props,out,
     &  logflag,logunit,logfile,mflag,message)
      if(mflag.ge.2)goto 900
c
c           debugging output of the computed slave node coordinates
c
      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'      num_surfaces =',num_surfaces
        write(logunit,*)'       numtiednode =',numtiednode
        write(logunit,*)
        write(logunit,1)'row','node_id','elem_id','face_id','flag',
     &                  'xi_iso1','eta_iso2','zeta_iso3',
     &                  'x_dof1','y_dof2','z_dof3'
        do row=1,numtiednode
          write(logunit,4)row,tiednodeid(row,1:maxtieddatacol),
     &              tiednodeisocoord(row,1:maxisodof),
     &              tiednodeglobalcoord(row,1:dimcoord)
        end do
        close(unit=logunit)
      end if

1     format(a8,3a10,a5,3a16,3a16)
4     format(i8,3i10,i5,3f16.8,3es16.6)
c
c----------------------------------------------------------------------
c           jump here on an error for exit procedures
c
900   continue

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'----- end tied_nodelocate -----'
        close(unit=logunit)
      end if

      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine tied_nodecompute                  *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 9/3/2017 rhd              *
c     *                                                              *
c     *  compute the tied contact constraint equation coefficient    *
c     *  values for the slave nodes (the dependent, constrained      *
c     *  nodes) to tie the slave nodes to the master element face    *
c     *  node degrees of freedom (the independent nodes); evaluate   *
c     *  the master element shape functions at the slave node        *
c     *  isoparametric coordinate location on the master element     *
c     *  face to get the constraint equation coefficient values.     *
c     *                                                              *
c     *  12/22/03 added a minimum value check to only save the mpc   *
c     *  multiplier value if it is greater than epsilon_mpc; added   *
c     *  epsilon_mpc to the call statement                           *
c     *                                                              *
c     ****************************************************************
c
      subroutine tied_nodecompute(
     &  tiednodeid,tiednodeisocoord,
     &  maxtiednode,maxtieddatacol,maxisodof,numtiednode,
     &  elemcon,dimelem,maxnodeperelem,numelem,
     &  numnode,dimcoord,meshformat,props,cstmap,epsilon_mpc,
     &  logflag,logunit,logfile,mflag,message)
c
c           declare modules
c
      use mod_mpc
      use global_data, only : max_mpc_tied
      implicit none
c
c           declare variables
c
      integer
     &  maxtiednode,maxtieddatacol,maxisodof,numtiednode,
     &  dimelem,maxnodeperelem,numelem,
     &  numnode,dimcoord,meshformat,
     &  logflag,logunit,mflag,
     &  tiednodeid(maxtiednode,maxtieddatacol),
     &  elemcon(dimelem,maxnodeperelem),
     &  cstmap(*)

      double precision
     &  tiednodeisocoord(maxtiednode,maxisodof)

      real props(mxelpr,*),epsilon_mpc

      character(len=*) ::
     &  logfile,message
c
c           variables:
c
c           tiednodeid() = list the slave node id in column 1, list the
c             master element id that the slave node is connected to in
c             column 2, list the master element face number in column 3,
c             give the node tied flag 1=yes, 0=no in column 4, if the
c             node is within the gap tolerance distance then set the
c             flag = 1 = yes to be included in the tied contact
c           tiednodeisocoord()  = element local isoparametric
c             coordinates (xi,eta,zeta) of the slave nodes on the master
c             element face; the local isoparametric coordinates are used
c             to compute the constraint equation coefficient values
c           maxtiednode = array row dimension; the maximum number of
c             tied slave nodes
c           maxtieddatacol = array column dimension
c           maxisodof = array column dimension; expecting 3 for the
c             xi,eta,zeta isoparametric coordinate dof values
c           numtiednode = return the number of tied slave nodes
c           elemcon() = array of the element connectivity node list for
c             each element in the mesh from the fea input file
c           dimelem = row dimension for the elemcon() array; the maximum
c             number of elements
c           maxnodeperelem = column dimension for the elemcon() array;
c             the maximum number of nodes per element (usually 20 or 27
c             nodes per element for a higher order brick element)
c           numelem = the given number of elements in the mesh from the
c             fea input file
c           meshformat = format of the fea input file;
c             1=abaqus input data file; 2=warp3d input data file
c           numnode = the given number of nodes in the mesh from the fea
c             input file
c           dimcoord = the maximum number of degrees-of-freedom (dof)
c             per node, usually 3 for x,y,z dof per node
c           cstmap() = given constraint flag for each dof,
c             0=unconstrained and active dof
c           epsilon_mpc = small value to avoid saving very small
c             mpc coefficient multiplier values, only save the mpc
c             coefficient value if it is larger than epsilon_mpc
c
c           note: see also the mpc module for the module variables
c
c           logflag,logunit,logfile = log file for debugging and testing;
c             the flag, file unit, and file name for the log file
c           mflag,message = error or warning message; integer flag
c             0=no message, 1=warning, 2=error; message text string
c
c
c           local variables
c
      integer, parameter :: maxelemface = 6
      integer, parameter :: maxnodeperface = 9
      integer
     &  debug_level,allocate_error,row,total_dof,numequations,dof,i,
     &  eqn,node_id,tied_node_row,tied_node_flag,local_dof,elem_id,
     &  face_id,node_per_elem,elemchoice,last_node,face_node,column,
     &  master_node_id,master_node_dof,master_eqn,count,iword,
     &  numelemface,numnodeperface,numcornerperface,maxlist,numlist,
     &  localfacenode(maxelemface,maxnodeperface)
      integer, allocatable :: dof_eqn_map(:)
      integer, allocatable :: eqn_node_map(:)
      logical
     &  local_debug,warp3d_mpc_export,abaqus_mpc_export

      double precision
     &  row_sum,epsilon,zero,one,
     &  shapefunc,norm_mult
      real rword
      data zero, one
     &  / 0.0d00, 1.0d00 /
      allocatable
     &  shapefunc(:)
      character text*512,word*10,xyz_letter(3)*1
c
c           variables to sort the mpc node list for testing
c           when sort_list_flag = true
c
      integer, parameter :: maxrow = mxndel      ! max 20 nodes per element
      integer, parameter :: maxcol = 2
      integer numrow,intarray(maxrow,maxcol)
      logical sort_list_flag,logicallist(maxrow)
      real real4list(maxrow)

      equivalence (rword,iword)
c
      debug_level = 2
      local_debug = .false.
      sort_list_flag = .false.

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'----- begin tied_nodecompute -----'
        write(logunit,*)'          maxtiednode =',maxtiednode
        write(logunit,*)'       maxtieddatacol =',maxtieddatacol
        write(logunit,*)'            maxisodof =',maxisodof
        write(logunit,*)'          numtiednode =',numtiednode
        write(logunit,*)'              dimelem =',dimelem
        write(logunit,*)'       maxnodeperelem =',maxnodeperelem
        write(logunit,*)'              numelem =',numelem
        write(logunit,*)'           meshformat =',meshformat
        write(logunit,*)'              numnode =',numnode
        write(logunit,*)'             dimcoord =',dimcoord
        write(logunit,*)'          epsilon_mpc =',epsilon_mpc
        write(logunit,*)'       sort_list_flag = ',sort_list_flag
        if(local_debug)then
          write(logunit,*)'          local_debug = ',local_debug
          write(logunit,*)'tied node data and isop. position from the'
          write(logunit,*)'tiednodeid() and tiednodeisocoord() arrays'
          write(logunit,1)'row','node_id','elem_id','face_id','flag',
     &                    'xi_iso1','eta_iso2','zeta_iso3'
          do row=1,numtiednode
            write(logunit,2)row,tiednodeid(row,1:maxtieddatacol),
     &                      tiednodeisocoord(row,1:maxisodof)
          end do
        end if
        close(unit=logunit)
      end if
1     format(4a10,a5,6a16)
2     format(4i10,i5,3f16.8)
c
c           error checking
c
      if(numnode.le.0)then
        mflag = 2
        message = 'error: numnode = 0 (tied_nodecompute).'
      else if(dimcoord.le.0)then
        mflag = 2
        message = 'error: dimcoord = 0 (tied_nodecompute).'
      else if(numelem.le.0)then
        mflag = 2
        message = 'error: numelem = 0 (tied_nodecompute).'
      else if(numelem.gt.dimelem)then
        mflag = 2
        message = 'error: numelem > dimelem (tied_nodecompute).'
      else if(maxtiednode.le.0)then
        mflag = 2
        message = 'error: maxtiednode = 0 (tied_nodecompute).'
      else if(maxtieddatacol.le.0)then
        mflag = 2
        message = 'error: maxtieddatacol = 0 (tied_nodecompute).'
      else if(maxisodof.le.0)then
        mflag = 2
        message = 'error: maxisodof = 0 (tied_nodecompute).'
      else if(maxnodeperelem.le.0)then
        mflag = 2
        message = 'error: maxnodeperelem = 0 (tied_nodecompute).'
      end if
      if(mflag.ge.2)goto 900
c
c           initialize values
c
      allocate(shapefunc(maxnodeperelem),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "shapefunc()" array'//
     &    ' (tied_nodecompute).'
        goto 900
      end if
c
c           get arrays to map from local dof to global equation number
c           and to map from the global equation number to the global
c           node number; constrained nodes are removed from the active
c           model dof; use the total number of dof for the local array
c           dimension
c
      total_dof = numnode*dimcoord
c
      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'            total_dof =',total_dof
        close(unit=logunit)
      end if
c
      allocate(dof_eqn_map(total_dof),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "dof_eqn_map()" array'//
     &    ' (tied_nodecompute).'
        goto 900
      end if
c
      allocate(eqn_node_map(total_dof),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "eqn_node_map()" array'//
     &    ' (tied_nodecompute).'
        goto 900
      end if
c
      call dof_map(dof_eqn_map,cstmap,numnode,dimcoord,
     &             eqn_node_map,numequations)

      if(local_debug .and. logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'constraint and dof arrays'
        write(logunit,*)'              numnode =',numnode
        write(logunit,*)'             dimcoord =',dimcoord
        write(logunit,*)'            total_dof =',total_dof
        write(logunit,3)'dof','nodeid','nodedof','cstmap','dof_eqn_map'
        do dof=1,total_dof
          node_id = ((dof-1)/dimcoord) + 1
          local_dof = dof - (node_id-1)*dimcoord
          write(logunit,4)dof,node_id,local_dof,cstmap(dof),
     &                    dof_eqn_map(dof)
        end do
        write(logunit,*)'         numequations =',numequations
        write(logunit,5)'eq_num','eqn_node_map'
        do i=1,numequations
          write(logunit,6)i,eqn_node_map(i)
        end do
        close(unit=logunit)
      end if
3     format(3a8,2a14)
4     format(3i8,2i14)
5     format(a8,2a14)
6     format(i8,2i14)
c
c           compute the tied contact node constraint equation
c           coefficient values; check for active global dof (nodes
c           that are unconstrained) if the node id for the global
c           equation matches a dependent (slave) node id, then compute
c           the constraint coefficient values for that global equation;
c           will also need to check that the corresponding dof of the
c           master nodes are active equations
c
      last_node = 0

      loop_1: do dof=1,total_dof
        eqn = dof_eqn_map(dof)

        if(eqn.gt.0)then
          node_id = eqn_node_map(eqn)
        else
          cycle loop_1
        end if

        tied_node_row = 0
        tied_node_flag = 0

        loop_2: do row=1,numtiednode
          if(tiednodeid(row,1).eq.node_id)then
            tied_node_row = row
            tied_node_flag = tiednodeid(row,4)
            exit loop_2
          end if
        end do loop_2

        if(tied_node_row.gt.0 .and. tied_node_flag.gt.0)then
          local_dof = dof - (node_id-1)*dimcoord
          elem_id = tiednodeid(row,2)
          face_id = tiednodeid(row,3)

          if(local_dof.le.0 .or. local_dof.gt.dimcoord)then
            mflag = 2
            message = 'internal error: invalid local_dof for the'//
     &         ' slave node; expecting 1 <= local_dof <= dimcoord'//
     &         ' (tied_nodecompute).'
          else if(elem_id.le.0 .or. elem_id.gt.numelem)then
            mflag = 2
            message = 'internal error: invalid elem_id from the'//
     &         ' tiednodeid() array (tied_nodecompute).'
          else if(face_id.le.0)then
            mflag = 2
            message = 'internal error: invalid face_id from the'//
     &         ' tiednodeid() array (tied_nodecompute).'
          end if
          if(mflag.ge.2)goto 900

          rword = props(1,elem_id)
          elemchoice = iword

          rword = props(2,elem_id)
          node_per_elem = iword
c
c               evaluate the master element shape functions at the
c               slave node's isoparametric location; if the node id
c               has not changed then use the previous shape functions
c
          if(node_id.ne.last_node)then
            call shapef(elemchoice,tiednodeisocoord(row,1),
     &                tiednodeisocoord(row,2),tiednodeisocoord(row,3),
     &                shapefunc)
c
c               get the local element face node numbers (the columns
c               in the connectivity array) for the master element face
c
            select case(meshformat)
            case(1)      ! abaqus
               call srt_setlocalfacenodes(localfacenode,maxelemface,
     &                maxnodeperface,numelemface,numnodeperface,
     &                numcornerperface,node_per_elem,elemchoice,
     &                logflag,logunit,logfile,mflag,message)
               if(mflag.ge.2)goto 900
            case(2)      ! warp3d
              call srt_setwarp3dlocalfacenodes(localfacenode,
     &               maxelemface,maxnodeperface,numelemface,
     &               numnodeperface,numcornerperface,node_per_elem,
     &               elemchoice,
     &               logflag,logunit,logfile,mflag,message)
              if(mflag.ge.2)goto 900
            case default
              mflag = 2
              message = 'error: unexpected meshformat choice;'//
     &              ' cannot get local element face nodes'//
     &              ' (tied_nodecompute).'
            end select      ! meshformat
            if(mflag.ge.2)goto 900
          end if
c
c               allocate the array components for this mpc equation
c               in the tied_con_mpc_table() array; set the maximum
c               number of terms for the mpc equation, maxlist, to the
c               number of nodes on the master element face plus one
c               for the slave node -1.0 coefficient value
c
          num_tied_con_mpc = num_tied_con_mpc + 1
          maxlist = numnodeperface + 1

          if(num_tied_con_mpc.gt.max_mpc_tied)then
            call srt_mpcs_resize
c           mflag = 2
            message = 'error: too many mpc tied contact constraint'//
     &              ' equations; num_tied_con_mpc > max_mpc;'//
     &              ' increase the max_mpc array dimension'//
     &              ' (tied_nodecompute).'
          end if
          if(mflag.ge.2)goto 900

          allocate(
     &       tied_con_mpc_table(num_tied_con_mpc)%node_list(maxlist),
     &       stat=allocate_error)
          if(allocate_error.gt.0)then
            mflag = 2
            message = 'internal error: memory error; array allocation'//
     &        ' error status > 0; check the'//
     &        ' "tied_con_mpc_table()%node_list()" array component'//
     &        ' of the tied_con_mpc_table data type'//
     &        ' (tied_nodecompute).'
            goto 900
          end if

          allocate(
     &       tied_con_mpc_table(num_tied_con_mpc)%dof_list(maxlist),
     &       stat=allocate_error)
          if(allocate_error.gt.0)then
            mflag = 2
            message = 'internal error: memory error; array allocation'//
     &        ' error status > 0; check the'//
     &        ' "tied_con_mpc_table()%dof_list()" array component'//
     &        ' of the tied_con_mpc_table data type'//
     &        ' (tied_nodecompute).'
            goto 900
          end if

          allocate(
     &    tied_con_mpc_table(num_tied_con_mpc)%multiplier_list(maxlist),
     &    stat=allocate_error)
          if(allocate_error.gt.0)then
            mflag = 2
            message = 'internal error: memory error; array allocation'//
     &        ' error status > 0; check the'//
     &        ' "tied_con_mpc_table()%multiplier_list()" array'//
     &        ' component of the tied_con_mpc_table data type'//
     &        ' (tied_nodecompute).'
            goto 900
          end if
c
c               save the slave node dof -1.0 equation coefficient;
c               set the tied node mpc equation constant to zero
c               (right hand side constant terms is zero for the tied
c               node mpc equation);
c               save the slave node id and slave node local dof to
c               the tied_con_mpc_table() array
c
          numlist = 1
          tied_con_mpc_table(num_tied_con_mpc)%constant = 0.0 ! single
          tied_con_mpc_table(num_tied_con_mpc)%node_list(numlist) =
     &      node_id
          tied_con_mpc_table(num_tied_con_mpc)%dof_list(numlist) =
     &      local_dof
          tied_con_mpc_table(num_tied_con_mpc)%multiplier_list(numlist)=
     &      -1.0 ! single
c
c               check that the master element face nodes have active
c               global equations for the current local dof; note that
c               the local node dof is the same for the slave and master
c               nodes, for example, the x=u dof would be tied between
c               the slave and master nodes
c
          loop_3: do face_node=1,numnodeperface
            column = localfacenode(face_id,face_node)

            if(column.le.0 .or. column.gt.maxnodeperelem)then
              mflag = 2
              message = 'error: invalid local element face node;'//
     &              ' expecting 1 <= column <= maxnodeperelem'//
     &              ' (tied_nodecompute).'
              goto 900
            end if

            master_node_id = elemcon(elem_id,column)
            master_node_dof = (master_node_id-1)*dimcoord + local_dof

            if(master_node_id.le.0 .or. master_node_id.gt.numnode)then
              mflag = 2
              message = 'error: invalid node id from the master;'//
     &              ' element; expecting 1 <= master_node_id <='//
     &              ' numnode (tied_nodecompute).'
            else if(master_node_dof.le.0 .or.
     &              master_node_dof.gt.total_dof)then
              mflag = 2
              message = 'error: invalid global dof from a node on'//
     &              ' the master element; expecting 1 <='//
     &              ' master_node_dof <= total_dof (tied_nodecompute).'
            end if
            if(mflag.ge.2)goto 900
c
c               save the shape function value for each master element
c               face node as the tied contact constraint equation
c               coefficient for the master element face node;
c               the master and slave node dof are the same for each
c               tied node mpc equation, for example, tie the x dof
c               of the master and slave nodes;
c               added a check for the master node mpc coefficient value
c               to avoid saving very small coefficient values to reduce
c               the number of mpc coefficients per slave node
c               (12/22/03 gvt)
c
            master_eqn = dof_eqn_map(master_node_dof)

            if(master_eqn.gt.0 .and.
     &      ABS(shapefunc(column)).gt.epsilon_mpc)then
              numlist = numlist + 1
              if(numlist.gt.maxlist)then
                mflag = 2
                message = 'error: too many terms in the tied node'//
     &              ' mpc equation; numlist > maxlist;'//
     &              ' increase the local maxlist array dimension'//
     &              ' (tied_nodecompute).'
                goto 900
              end if

              tied_con_mpc_table(num_tied_con_mpc)%node_list(numlist) =
     &          master_node_id
              tied_con_mpc_table(num_tied_con_mpc)%dof_list(numlist) =
     &          local_dof
              tied_con_mpc_table(num_tied_con_mpc)%
     &          multiplier_list(numlist) = real( shapefunc(column) )

            else if(master_eqn.gt.0 .and.
     &      ABS(shapefunc(column)).le.epsilon_mpc)then
              if(logflag.ge.1)then
                open(unit=logunit,file=logfile,position='append')
                write(logunit,*)
                write(logunit,*)'info: mpc coefficient value too small;'
                write(logunit,*)'value not saved'
                write(logunit,*)'    master_node_id =',master_node_id
                write(logunit,*)'   master_node_dof =',master_node_dof
                write(logunit,*)'        master_eqn =',master_eqn
                write(logunit,*)'            column =',column
                write(logunit,*)'       shapefunc() =',shapefunc(column)
                write(logunit,*)'       epsilon_mpc =',epsilon_mpc
                write(logunit,*)
                close(unit=logunit)
              end if

            else
              mflag = 1
              message = 'warning: a global dof for a master node'//
     &              ' does not have an active global equation;'//
     &              ' the node may be constrained;'//
     &              ' the constraint coefficient is not saved;'//
     &              ' check the log file (tied_nodecompute).'
              if(logflag.ge.1)then
                open(unit=logunit,file=logfile,position='append')
                write(logunit,*)
                write(logunit,*)'    mflag =',mflag
                write(logunit,'(a)')trim(message)
                write(logunit,*)'    master_node_id =',master_node_id
                write(logunit,*)'   master_node_dof =',master_node_dof
                write(logunit,*)'        master_eqn =',master_eqn
                write(logunit,*)
                close(unit=logunit)
              end if
            end if
          end do loop_3
c
c           save the number of terms in the tied node mpc equation
c
          tied_con_mpc_table(num_tied_con_mpc)%num_terms = numlist
          last_node = node_id
c
c           normalize the mpc equation coefficient values so that the
c           sum of the master node terms equals 1.0; the slave node
c           coefficient is saved in location 1, see above for numlist = 1;
c           expecting the master dof coefficient values to be saved in
c           locations 2 to numlist; normalize the numlist-1 coefficient
c           values for this slave node (12/22/03 GVT)
c
          if(local_debug .and. logflag.ge.debug_level)then
            open(unit=logunit,file=logfile,position='append')
            write(logunit,*)
            write(logunit,*)'before normalize coefficient'
            write(logunit,*)'  num_tied_con_mpc =',num_tied_con_mpc
            write(logunit,*)'           numlist =',numlist
            write(logunit,11)'term','node_id','dof','coefficient'
            row_sum = zero
            do column=1,numlist
              write(logunit,12)column,
     &          tied_con_mpc_table(num_tied_con_mpc)%node_list(column),
     &          tied_con_mpc_table(num_tied_con_mpc)%dof_list(column),
     &          tied_con_mpc_table(num_tied_con_mpc)%
     &          multiplier_list(column)
                row_sum = row_sum +
     &          tied_con_mpc_table(num_tied_con_mpc)%
     &          multiplier_list(column)
            end do
            write(logunit,*)'  equation coefficient row_sum =',row_sum
            close(unit=logunit)
          end if

          row_sum = zero
          do i=2,numlist
            row_sum = row_sum +
     &        tied_con_mpc_table(num_tied_con_mpc)%multiplier_list(i)
          end do

          if(row_sum.gt.zero)then
            norm_mult = one/row_sum            ! normalize multiplier
            do i=2,numlist
              tied_con_mpc_table(num_tied_con_mpc)%multiplier_list(i) =
     &        tied_con_mpc_table(num_tied_con_mpc)%multiplier_list(i)*
     &        real( norm_mult )
            end do
          else
            norm_mult = zero
            mflag = 1
            message = 'warning: cannot normalize the mpc equation'//
     &              ' coefficient values; the coefficient row sum'//
     &              ' value is less than or equal to zero;'//
     &              ' the given epsilon_mpc value may be too large'//
     &              ' and may be causing valid mpc coefficient'//
     &              ' values to be eliminated;'//
     &              ' check the log file (tied_nodecompute).'
            if(logflag.ge.1)then
              open(unit=logunit,file=logfile,position='append')
              write(logunit,*)
              write(logunit,*)'    mflag =',mflag
              write(logunit,'(a)')trim(message)
              write(logunit,*)'        equation # =',num_tied_con_mpc
              write(logunit,*)'     slave node_id =',node_id
              write(logunit,*)'      coef row sum =',row_sum
              write(logunit,*)'  num master nodes =',numlist
              write(logunit,*)'       epsilon_mpc =',epsilon_mpc
              write(logunit,*)
              close(unit=logunit)
            end if
          end if

          if(local_debug .and. logflag.ge.debug_level)then
            open(unit=logunit,file=logfile,position='append')
            write(logunit,*)
            write(logunit,*)'after normalize coefficient'
            write(logunit,*)'      coef row_sum =',row_sum
            write(logunit,*)'         norm_mult =',norm_mult
            write(logunit,*)'  num_tied_con_mpc =',num_tied_con_mpc
            write(logunit,*)'           numlist =',numlist
            write(logunit,11)'term','node_id','dof','coefficient'
            row_sum = zero
            do column=1,numlist
              write(logunit,12)column,
     &          tied_con_mpc_table(num_tied_con_mpc)%node_list(column),
     &          tied_con_mpc_table(num_tied_con_mpc)%dof_list(column),
     &          tied_con_mpc_table(num_tied_con_mpc)%
     &          multiplier_list(column)
                row_sum = row_sum +
     &          tied_con_mpc_table(num_tied_con_mpc)%
     &          multiplier_list(column)
            end do
            write(logunit,*)'  equation coefficient row_sum =',row_sum
            close(unit=logunit)
          end if
        end if
c
c           for testing, sort the mpc list values by node id;
c           allows for easier comparison of mpc equations between
c           abaqus and warp3d input files that use different element
c           face numbering
c
        if(sort_list_flag)then
          intarray = 0
          logicallist = .false.
          real4list = 0.0

          numrow = tied_con_mpc_table(num_tied_con_mpc)%num_terms

          do i=1,numrow
            intarray(i,1) =
     &        tied_con_mpc_table(num_tied_con_mpc)%node_list(i)
            intarray(i,2) =
     &        tied_con_mpc_table(num_tied_con_mpc)%dof_list(i)
            real4list(i) =
     &        tied_con_mpc_table(num_tied_con_mpc)%multiplier_list(i)
          end do

          call srt_sortintmultiarray(intarray,maxrow,maxcol,numrow,1,
     &                   logicallist,real4list,1)

          do i=1,numrow
            tied_con_mpc_table(num_tied_con_mpc)%node_list(i) =
     &        intarray(i,1)
            tied_con_mpc_table(num_tied_con_mpc)%dof_list(i) =
     &        intarray(i,2)
            tied_con_mpc_table(num_tied_con_mpc)%multiplier_list(i) =
     &        real4list(i)
          end do
        end if
      end do loop_1
c
c           write the constraint coefficient matrix values and the
c           constraint equations with node id numbers to the log file
c
      if(logflag.ge.debug_level)then
        warp3d_mpc_export = .true.
        abaqus_mpc_export = .true.

        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'tied node mpc equation coefficient values:'
        write(logunit,*)'     num_tied_con_mpc =',num_tied_con_mpc
        do row=1,num_tied_con_mpc
          numlist = tied_con_mpc_table(row)%num_terms
          write(logunit,*)
          write(logunit,*)'     equation #',row
          write(logunit,*)'      numlist =',numlist
          write(logunit,11)'term','node_id','dof','coefficient'
          row_sum = zero
          do column=1,numlist
            write(logunit,12)column,
     &        tied_con_mpc_table(row)%node_list(column),
     &        tied_con_mpc_table(row)%dof_list(column),
     &        tied_con_mpc_table(row)%multiplier_list(column)
            row_sum = row_sum +
     &        tied_con_mpc_table(row)%multiplier_list(column)
          end do
          write(logunit,*)'  equation coefficient row_sum =',row_sum
        end do
11      format(a6,a10,a5,a17)
12      format(i7,i10,i5,es17.7)

        xyz_letter(1) = 'x'
        xyz_letter(2) = 'y'
        xyz_letter(3) = 'z'

        write(logunit,*)
        write(logunit,*)'tied node mpc equations'
        write(logunit,*)'stored in the tied_con_mpc_table() data type'
        write(logunit,*)'  number of mpc equations =',num_tied_con_mpc

        do row=1,num_tied_con_mpc
          eqn = 0

          row_sum = zero
          count = 0
          text = ' '
          numlist = tied_con_mpc_table(row)%num_terms
          do column=1,numlist
            if(tied_con_mpc_table(row)%multiplier_list(column)
     &      .ne.zero .or. .true.)then
              count = count + 1
              row_sum = row_sum +
     &          tied_con_mpc_table(row)%multiplier_list(column)
              write(word,'(f10.4)')
     &          tied_con_mpc_table(row)%multiplier_list(column)
              if(column.gt.1)then
                text = trim(text)//' + ('//trim(adjustl(word))//')*'
              else
                text = trim(text)//'('//trim(adjustl(word))//')*'
              end if
              node_id = tied_con_mpc_table(row)%node_list(column)
              local_dof = tied_con_mpc_table(row)%dof_list(column)
              write(word,'(i10)')node_id
              if(local_dof.ge.1 .and. local_dof.le.3)then
                text = trim(text)//'node '//trim(adjustl(word))//
     &                   xyz_letter(local_dof)
              else
                text = trim(text)//'node '//trim(adjustl(word))
              end if
            end if
          end do
          text = trim(text)//' = 0.0'
          write(logunit,8)'equation ',row,'   sum =',row_sum,
     &                      '  num =',count,':   ',trim(text)
        end do
c
c           write the tied contact mpc equations in warp3d format
c           to support testing of the tied contact and user defined
c           mpc in warp3d; the listed user mpc equations can replace
c           the corresponding tied contact data in the warp3d input
c           file
c
c           warp3d format is:
c           <nodes> <multiplier> <u,v,w> + (,) ... + (,) = <constant>
c
c           for example:
c           6 1.5 u + 6 3.2 v + 6 -3.1 w = 0.0
c
        if(warp3d_mpc_export)then
          epsilon = 1.0d-15
          call tied_exportmpc_warp3dformat(epsilon,logunit)
        end if      ! warp3d_mpc_export
c
c           write the tied contact mpc equations in abaqus format
c           to support testing of the tied contact and user defined
c           mpcs; the listed user mpc equations can replace
c           the corresponding tied contact data in the abaqus input
c           file; note that the first node is removed as a degree of
c           freedom to impose the constraint, don't apply any other
c           boundary conditions or constraints to that node;
c           the slave node should already be first in the equation list;
c           the data line can have a maximum of four node terms per
c           line, use extra lines for more terms in the equation
c
c           abaqus format is:
c           line 1: n (the number of nodes in the equation)
c           line 2: node, dof, mult, ... (terms for each node)
c
c           for example:
c           3
c           6, 1, 1.5, 6, 2, 3.2, 6, 3, -3.1
c           (node 6, x, 1.5, node 6, y, 3.2, node 6, z, -3.1 = 0.0)
c
        if(abaqus_mpc_export)then
          epsilon = 1.0d-5
          call tied_exportmpc_abaqusformat(epsilon,logunit)
        end if      ! abaqus_mpc_export

        close(unit=logunit)
      end if

8     format(a,i9,a,f12.6,a,i3,a,a)
c
c----------------------------------------------------------------------
c           jump here on an error for exit procedures
c
900   continue
c
c           deallocate local arrays
c
      if(allocated(dof_eqn_map))deallocate(dof_eqn_map)
      if(allocated(eqn_node_map))deallocate(eqn_node_map)
      if(allocated(shapefunc))deallocate(shapefunc)

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'----- end tied_nodecompute -----'
        close(unit=logunit)
      end if

      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine tied_nodeassign                  *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 11/04/02                  *
c     *                                    12/23/02 gvt              *
c     *                                    02/26/03 gvt              *
c     *                                    03/06/03 gvt              *
c     *                                    05/12/04 gvt              *
c     *                                                              *
c     *  assign each slave node to a master element face;            *
c     *  the slave (dependent) nodes are listed in the tiednodeid()  *
c     *  array; use the element connectivity and node coordinates    *
c     *  to compute an approximate master element face centroid      *
c     *  location; assign each slave node to the nearest master      *
c     *  element face centroid point                                 *
c     *                                                              *
c     *  added a distance check to the logic that checks for the     *
c     *  best master element (03/05/03 gvt)                          *
c     *                                                              *
c     *  added checkstatus = 0 initialization to avoid a warning     *
c     *  from the HP compiler (03/06/03 gvt)                         *
c     *                                                              *
c     *  added the tiednodegaptolerance() array to the call          *
c     *  statement (05/12/04 gvt)                                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine tied_nodeassign(tiednodeid,tiednodeadjustflag,
     &  tiednodegaptolerance,
     &  maxtiednode,maxtieddatacol,numtiednode,
     &  elemcon,dimelem,maxnodeperelem,numelem,
     &  nodecoord,dimnode,dimcoord,numnode,meshformat,props,out,
     &  logflag,logunit,logfile,mflag,message)
c
c           declare modules
c
      use mod_mpc
      implicit none
c
c           declare variables
c
      integer
     &  maxtiednode,maxtieddatacol,numtiednode,
     &  dimelem,maxnodeperelem,numelem,
     &  dimnode,dimcoord,numnode,meshformat,out,
     &  logflag,logunit,mflag,
     &  tiednodeid(maxtiednode,maxtieddatacol),
     &  elemcon(dimelem,maxnodeperelem)

      double precision
     &  nodecoord(dimnode,dimcoord)

      real props(mxelpr,*),
     &  tiednodegaptolerance(maxtiednode)

      logical
     &  tiednodeadjustflag(maxtiednode)

      character(len=*) ::
     &  logfile,message
c
c           local variables
c
      integer, parameter :: master_surface_flag = 1
      integer, parameter :: maxelemface = 6
      integer, parameter :: maxnodeperface = 9
      integer, parameter :: maxdof = 3
      integer
     &  debug_level,allocate_error,corner,nd1,dist_count,iword,
     &  surfaceid,row,element_id,face_id,node_per_elem,numnodeperface,
     &  maxmasterlist,maxmastercol,nummasterlist,
     &  maxelempersurface,maxmastersurface,
     &  numelemlist,numelemface,nd,dof,node_id,elemchoice,
     &  numcornerperface,checkstatus,best_status,
     &  best_elem_row,closest_elem_row,tiedset,numsurfpairs,
     &  max_masternodelist,num_masternodelist,face_node,
     &  localfacenode(maxelemface,maxnodeperface)
      integer, allocatable :: elemlist(:,:)
      integer, allocatable :: surfaceflagtable(:)
      integer, allocatable :: masternodelist(:)

      logical
     &  local_debug,inside_flag,dist_gap_flag,check1_flag,check2_flag,
     &  check3_flag,check_redundant_flag

      double precision
     &  zero,distance,given_epsilon,face_dist,
     &  max_elem_dist,sum_elem_dist,closest_distance,
     &  elemfacecentroid,normal_dist,gap_tolerance,
     &  best_distance
      real rword
      data zero
     &  / 0.0d00 /
      allocatable
     &  elemfacecentroid(:,:)
c
c           declare functions
c
      double precision
     &  srt_dist3d

      equivalence (rword,iword)
c
      debug_level = 7
      local_debug = .false.

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'----- begin tied_nodeassign -----'
        write(logunit,*)'         num_surfaces =',num_surfaces
        write(logunit,*)'          maxtiednode =',maxtiednode
        write(logunit,*)'       maxtieddatacol =',maxtieddatacol
        write(logunit,*)'          numtiednode =',numtiednode
        write(logunit,*)'              dimelem =',dimelem
        write(logunit,*)'       maxnodeperelem =',maxnodeperelem
        write(logunit,*)'              numelem =',numelem
        write(logunit,*)'              dimnode =',dimnode
        write(logunit,*)'             dimcoord =',dimcoord
        write(logunit,*)'              numnode =',numnode
        write(logunit,*)'           meshformat =',meshformat
        if(local_debug)then
          write(logunit,*)'          local_debug = ',local_debug
        end if
        close(unit=logunit)
      end if
c
c           error checking
c
      if(numnode.le.0)then
        mflag = 2
        message = 'error: numnode = 0 (tied_nodeassign).'
      else if(numnode.gt.dimnode)then
        mflag = 2
        message = 'error: numnode > dimnode (tied_nodeassign).'
      else if(dimcoord.le.0)then
        mflag = 2
        message = 'error: dimcoord = 0 (tied_nodeassign).'
      else if(dimcoord.lt.maxdof)then
        mflag = 2
        message = 'error: dimcoord < 3 (tied_nodeassign).'
      else if(numelem.le.0)then
        mflag = 2
        message = 'error: numelem = 0 (tied_nodeassign).'
      else if(numelem.gt.dimelem)then
        mflag = 2
        message = 'error: numelem > dimelem (tied_nodeassign).'
      else if(num_surfaces.le.0)then
        mflag = 2
        message = 'error: num_surfaces = 0 (tied_nodeassign).'
      else if(maxtiednode.le.0)then
        mflag = 2
        message = 'error: maxtiednode = 0 (tied_nodeassign).'
      else if(maxtieddatacol.le.0)then
        mflag = 2
        message = 'error: maxtieddatacol = 0 (tied_nodeassign).'
      else if(maxtieddatacol.lt.3)then
        mflag = 2
        message = 'error: maxtieddatacol < 3 (tied_nodeassign).'
      else if(numtiednode.le.0)then
        mflag = 2
        message = 'error: numtiednode = 0 (tied_nodeassign).'
      else if(numtiednode.gt.maxtiednode)then
        mflag = 2
        message = 'error: numtiednode > maxtiednode (tied_nodeassign).'
      end if
      if(mflag.ge.2)goto 900
c
c           initialize values
c
c           get the maximum number of elements for a mesh surface;
c           note that the number of elements per mesh surface is
c           stored in the surface_table()%num_elems data structure;
c           use the maximum number of elements to allocate the local
c           nodelist() array;
c
c           added checkstatus = 0 initialization to avoid a warning
c           from the HP compiler; note that checkstatus is also
c           initialized in the srt_insideelemface routine (03/06/03 gvt)
c
      check_redundant_flag = .true.
      checkstatus = 0
      maxelempersurface = 0
      maxmastersurface = 0

      do surfaceid=1,num_surfaces
        if(surface_table(surfaceid)%num_elems.gt.maxelempersurface)then
          maxelempersurface = surface_table(surfaceid)%num_elems
        end if
      end do

      do tiedset=1,num_tied_sets
        numsurfpairs = tied_contact_table(tiedset)%num_pairs
        if(numsurfpairs.gt.0)then
          maxmastersurface = maxmastersurface + numsurfpairs
        end if
      end do

      nummasterlist = 0
      maxmasterlist = maxelempersurface*maxmastersurface
      maxmastercol = 3
      max_masternodelist = maxmasterlist*maxnodeperelem
      num_masternodelist = 0

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'local values:'
        write(logunit,*)'     maxelempersurface =',maxelempersurface
        write(logunit,*)'      maxmastersurface =',maxmastersurface
        write(logunit,*)'         maxmasterlist =',maxmasterlist
        write(logunit,*)'          maxmastercol =',maxmastercol
        write(logunit,*)'    max_masternodelist =',max_masternodelist
        write(logunit,*)'    num_masternodelist =',num_masternodelist
        write(logunit,*)'  check_redundant_flag = ',check_redundant_flag
        close(unit=logunit)
      end if

      allocate(elemlist(maxmasterlist,maxmastercol),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "elemlist()" array'//
     &    ' (tied_nodeassign).'
        goto 900
      end if
      elemlist = 0

      allocate(elemfacecentroid(maxmasterlist,maxdof),
     &         stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "elemfacecentroid()"'//
     &    ' array (tied_nodeassign).'
        goto 900
      end if
      elemfacecentroid = zero

      if(check_redundant_flag)then
        allocate(masternodelist(max_masternodelist),stat=allocate_error)
        if(allocate_error.gt.0)then
          mflag = 2
          message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "masternodelist()"'//
     &    ' array (tied_nodeassign).'
          goto 900
        end if
        masternodelist = 0
      end if

      sum_elem_dist = zero
      dist_count = 0
c
c           set a flag for the master mesh surfaces as identified in
c           the tied_contact_table()%master_list() for each contact
c           pair
c
      allocate(surfaceflagtable(num_surfaces),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "surfaceflagtable()"'//
     &    ' array (tied_nodeassign).'
        goto 900
      end if
      surfaceflagtable = 0

      do tiedset=1,num_tied_sets
        numsurfpairs = tied_contact_table(tiedset)%num_pairs
        if(numsurfpairs.gt.0)then
          do row=1,numsurfpairs
            surfaceid = tied_contact_table(tiedset)%master_list(row)
            if(surfaceid.gt.0 .and. surfaceid.le.num_surfaces)then
              surfaceflagtable(surfaceid) = master_surface_flag
            else
              mflag = 2
              message = 'error: invalid mesh surface id number;'//
     &          ' check the tied_contact_table()%master_list()'//
     &          ' (tied_nodeassign).'
              goto 900
            end if
          end do
        end if
      end do

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'    surfaceid,    surfaceflagtable()'
        do surfaceid=1,num_surfaces
          write(logunit,*)surfaceid,surfaceflagtable(surfaceid)
        end do
        close(unit=logunit)
      end if
c
c           get a list of the elements on the master mesh surface;
c           check the surfaceflagtable() for master mesh surfaces with
c           flag = 1; get the approximate master element face centroid
c           location using the element face corner nodes
c
      do surfaceid=1,num_surfaces
        if(surfaceflagtable(surfaceid).eq.master_surface_flag)then
          if(logflag.ge.debug_level)then
            open(unit=logunit,file=logfile,position='append')
            write(logunit,*)
            write(logunit,*)'               mesh surface =',surfaceid
            write(logunit,*)'  surface_table()%num_elems =',
     &                         surface_table(surfaceid)%num_elems
            write(logunit,*)'         surface_table()%id = ',
     &                         trim(surface_table(surfaceid)%id)
            write(logunit,*)'         surfaceflagtable() =',
     &                         surfaceflagtable(surfaceid)
            write(logunit,*)
            close(unit=logunit)
          end if

          numelemlist = surface_table(surfaceid)%num_elems
          if(numelemlist.le.0)then
            mflag = 2
            message = 'error: numelemlist = 0 for a master mesh'//
     &                ' surface; check surface_table()%num_elems'//
     &                ' (tied_nodeassign).'
          end if
          if(mflag.ge.2)goto 900

          do row=1,numelemlist
c
c               get the element id number, the element face number,
c               and the number of nodes per element for this
c               master mesh surface element; use the nodes per element
c               to identify the element type
c
            element_id = surface_table(surfaceid)%elem_list(row)
            face_id = surface_table(surfaceid)%face_list(row)
            max_elem_dist = zero

            if(element_id.le.0)then
              mflag = 2
              message = 'error: element_id = 0 for a master mesh'//
     &                  ' element; check surface_table()%elem_list()'//
     &                  ' (tied_nodeassign).'
            else if(element_id.gt.numelem)then
              mflag = 2
              message = 'error: element_id > numelem for a master'//
     &                  ' element; check surface_table()%elem_list()'//
     &                  ' (tied_nodeassign).'
            else if(face_id.le.0)then
              mflag = 2
              message = 'error: face_id = 0 for a master mesh'//
     &                  ' element; check surface_table()%face_list()'//
     &                  ' (tied_nodeassign).'
            end if
            if(mflag.ge.2)goto 900

            rword = props(1,element_id)
            elemchoice = iword

            rword = props(2,element_id)
            node_per_elem = iword
c
c               get the local element face node numbers (the columns
c               in the connectivity array) for this element face
c
            select case(meshformat)
            case(1)      ! abaqus
              call srt_setlocalfacenodes(localfacenode,maxelemface,
     &               maxnodeperface,numelemface,numnodeperface,
     &               numcornerperface,node_per_elem,elemchoice,
     &               logflag,logunit,logfile,mflag,message)
              if(mflag.ge.2)goto 900
            case(2)      ! warp3d
              call srt_setwarp3dlocalfacenodes(localfacenode,
     &               maxelemface,maxnodeperface,numelemface,
     &               numnodeperface,numcornerperface,node_per_elem,
     &               elemchoice,
     &               logflag,logunit,logfile,mflag,message)
              if(mflag.ge.2)goto 900
            case default
              mflag = 2
              message = 'error: unexpected meshformat choice;'//
     &                  ' cannot get local element face nodes'//
     &                  ' (tied_nodeassign).'
            end select      ! meshformat
            if(mflag.ge.2)goto 900
c
c               save the master element id and element face id
c               to the local arrays
c
            nummasterlist = nummasterlist + 1
            if(nummasterlist.gt.maxmasterlist)then
              mflag = 2
              message = 'error: nummasterlist > maxmasterlist;'//
     &                  ' increase the local array dimension'//
     &                  ' maxmasterlist (tied_nodeassign).'
              goto 900
            end if
            elemlist(nummasterlist,1) = element_id
            elemlist(nummasterlist,2) = face_id
            elemlist(nummasterlist,3) = elemchoice
c
c               use the master element face corner nodes to compute an
c               approximate face centroid x,y,z location
c
            elemfacecentroid(nummasterlist,1:maxdof) = zero

            do corner=1,numcornerperface
              node_id =
     &            elemcon(element_id,localfacenode(face_id,corner))
              if(node_id.le.0)then
                mflag = 1
                message = 'warning: element face node id = 0;'//
     &                    ' skipping the master surface face node'//
     &                    ' when computing the element face'//
     &                    ' centroid location (tied_nodeassign).'
                if(logflag.ge.1)
     &            call srt_writelogmessage(mflag,message,logfile,
     &                                     logunit)
              else if(node_id.gt.numnode)then
                mflag = 1
                message = 'warning: element face node id > numnode;'//
     &                    ' skipping the master surface face node'//
     &                    ' when computing the element face'//
     &                    ' centroid location (tied_nodeassign).'
                if(logflag.ge.1)
     &            call srt_writelogmessage(mflag,message,logfile,
     &                                     logunit)
              end if
c
c               sum the master element approximate face centroid
c               location in the local array
c
              do dof=1,maxdof
                elemfacecentroid(nummasterlist,dof) =
     &              elemfacecentroid(nummasterlist,dof) +
     &              nodecoord(node_id,dof)
              end do
c
c               get the maximum element distance from corner 1 to the
c               other corners on the element face; use later to set an
c               epsilon distance
c
              if(corner.eq.1)then
                nd1 = node_id
              else
                face_dist = srt_dist3d(nodecoord(nd1,1:dimcoord),
     &                        nodecoord(node_id,1:dimcoord),dimcoord)
                if(face_dist.gt.max_elem_dist)then
                  max_elem_dist = face_dist
                end if
              end if
            end do
c
c               save the node id number on the master element face;
c               use later to check for redundant slave nodes
c               (6/10/04 gvt)
c
            if(check_redundant_flag)then
              do face_node=1,numnodeperface
                node_id =
     &            elemcon(element_id,localfacenode(face_id,face_node))
                num_masternodelist = num_masternodelist + 1

                if(num_masternodelist.gt.max_masternodelist)then
                  mflag = 1
                  message = 'error: too many nodes in the local'//
     &              ' master node list; increase the local array'//
     &              ' dimension; cannot check for redundant slave'//
     &              ' nodes (tied_nodeassign).'
                  goto 900
                end if

                masternodelist(num_masternodelist) = node_id
              end do
            end if
c
c               average the master element approximate face centroid
c               location in the local array
c
            do dof=1,maxdof
              elemfacecentroid(nummasterlist,dof) =
     &            elemfacecentroid(nummasterlist,dof)/
     &            dble(numcornerperface)
            end do

            dist_count = dist_count + 1
            sum_elem_dist = sum_elem_dist + max_elem_dist

          end do
        end if
      end do
c
c               remove duplicate master node id numbers from the master
c               node list before checking for redundant slave nodes
c
      if(check_redundant_flag)then
        call srt_checkduplicate(masternodelist,max_masternodelist,
     &    num_masternodelist,
     &    logflag,logunit,logfile,mflag,message)
        if(mflag.ge.2)goto 900
      end if

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'       dist_count =',dist_count
        write(logunit,*)'    sum_elem_dist =',sum_elem_dist
        write(logunit,*)'master surface, element centroid coord'
        write(logunit,*)'    nummasterlist =',nummasterlist
        write(logunit,*)'    maxmasterlist =',maxmasterlist
        write(logunit,8)'row','elem_id','face_id','type',
     &                  'x-coord','y-coord','z-coord'
        do row=1,nummasterlist
          write(logunit,9)row,elemlist(row,1),elemlist(row,2),
     &                elemlist(row,3),elemfacecentroid(row,1:maxdof)
        end do
        write(logunit,*)
        write(logunit,*)'    num_masternodelist =',num_masternodelist
        write(logunit,*)'    max_masternodelist =',max_masternodelist
        write(logunit,*)'  check_redundant_flag = ',check_redundant_flag
        if(check_redundant_flag)then
          write(logunit,*)'masternodelist():'
          write(logunit,12)'row','master_node_id'
          do row=1,num_masternodelist
            write(logunit,13)row,masternodelist(row)
          end do
        end if
        close(unit=logunit)
      end if
8     format(a8,3a8,3a16)
9     format(i8,3i8,3es16.6)
12    format(a8,a16)
13    format(i8,i16)
c
c               check for slave nodes in the master node list;
c               remove any redundant slave nodes to avoid generating
c               unneeded redundant mpc equations (6/10/04 gvt)
c
      if(check_redundant_flag)then
        call tied_nodecheckredundant(tiednodeid,tiednodeadjustflag,
     &    tiednodegaptolerance,maxtiednode,maxtieddatacol,numtiednode,
     &    masternodelist,max_masternodelist,num_masternodelist,out,
     &    logflag,logunit,logfile,mflag,message)
        if(mflag.ge.2)goto 900
      end if
c
c               for each slave node find the nearest master element face
c               centroid; assign the slave node to that master element
c
      do nd=1,numtiednode
        node_id = tiednodeid(nd,1)
        gap_tolerance = DBLE(tiednodegaptolerance(nd))
c
c               use the given gap tolerance distance to check if a node
c               is close to an element face
c
        given_epsilon = gap_tolerance

        if(local_debug .and.logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)
          write(logunit,*)'- - - - - - - - - - - - - - - - - - - - - -'
          write(logunit,*)'   node loop index, nd =',nd
          write(logunit,*)'         slave node_id =',node_id
          write(logunit,*)'         gap_tolerance =',gap_tolerance
          write(logunit,*)'         given_epsilon =',given_epsilon
          close(unit=logunit)
        end if

        if(node_id.le.0)then
          mflag = 2
          message = 'error: node_id = 0; invalid node id for the'//
     &              ' slave node; check the tiednodeid() array'//
     &              ' (tied_nodeassign).'
        else if(node_id.gt.numnode)then
          mflag = 2
          message = 'error: node_id > numnode; invalid node id for'//
     &              ' the slave node; check the tiednodeid() array'//
     &              ' (tied_nodeassign).'
        end if
        if(mflag.ge.2)goto 900
c
c               loop through the master element list to assign the best
c               master element to the current slave node
c
        do row=1,nummasterlist
c
c               compute the distance from this slave node to each master
c               element centroid point (for debugging); update the default
c               nearest master element when a nearer element is found;
c
          distance = srt_dist3d(nodecoord(node_id,1:maxdof),
     &                   elemfacecentroid(row,1:maxdof),maxdof)

          if(local_debug .and.logflag.ge.debug_level)then
            open(unit=logunit,file=logfile,position='append')
            write(logunit,*)
            write(logunit,*)'   node loop index, nd =',nd
            write(logunit,*)'         slave node_id =',node_id
            write(logunit,*)'         gap_tolerance =',gap_tolerance
            write(logunit,*)'  elem loop index, row =',row
            write(logunit,*)'        master elem id =',elemlist(row,1)
            write(logunit,*)'      master elem face =',elemlist(row,2)
            write(logunit,*)'             elem type =',elemlist(row,3)
            write(logunit,3)'       slave node coord =',
     &                               nodecoord(node_id,1:maxdof)
            write(logunit,3)'               distance =',distance
            close(unit=logunit)
          end if
3         format(a,3es16.6)
c
c               initialize values for the first trial master element
c
          if(row.eq.1)then
            best_status = 0
            best_elem_row = 0
            best_distance = zero
            closest_elem_row = 0
            closest_distance = zero
          end if

          if(.true.)then      ! check all master elements
c
c               check that the node is within the element face
c               boundaries
c
            element_id = elemlist(row,1)
            face_id = elemlist(row,2)
            elemchoice = elemlist(row,3)

            rword = props(2,element_id)
            node_per_elem = iword

            if(local_debug .and.logflag.ge.debug_level)then
              open(unit=logunit,file=logfile,position='append')
              write(logunit,*)
              write(logunit,*)'check if node within element face'
              write(logunit,*)'               row =',row
              write(logunit,*)'        element_id =',element_id
              write(logunit,*)'           face_id =',face_id
              write(logunit,*)'        elemchoice =',elemchoice
              write(logunit,*)'     node_per_elem =',node_per_elem
              close(unit=logunit)
            end if

            select case(meshformat)
            case(1)      ! abaqus
              call srt_setlocalfacenodes(localfacenode,maxelemface,
     &               maxnodeperface,numelemface,numnodeperface,
     &               numcornerperface,node_per_elem,elemchoice,
     &               logflag,logunit,logfile,mflag,message)
              if(mflag.ge.2)goto 900
            case(2)      ! warp3d
              call srt_setwarp3dlocalfacenodes(localfacenode,
     &               maxelemface,maxnodeperface,numelemface,
     &               numnodeperface,numcornerperface,node_per_elem,
     &               elemchoice,
     &               logflag,logunit,logfile,mflag,message)
              if(mflag.ge.2)goto 900
            case default
              mflag = 2
              message = 'error: unexpected meshformat choice;'//
     &                  ' cannot get local element face nodes'//
     &                  ' (tied_nodeassign).'
            end select      ! meshformat
            if(mflag.ge.2)goto 900

            call srt_insideelemface(inside_flag,checkstatus,element_id,
     &                     face_id,elemchoice,node_per_elem,
     &                     nodecoord(node_id,1:maxdof),maxdof,
     &                     localfacenode,maxelemface,maxnodeperface,
     &                     numelemface,numnodeperface,numcornerperface,
     &                     elemcon,dimelem,maxnodeperelem,numelem,
     &                     nodecoord,dimnode,dimcoord,numnode,
     &                     given_epsilon,normal_dist,
     &                     logflag,logunit,logfile,mflag,message)
            if(mflag.ge.2)goto 900

            if(local_debug .and.logflag.ge.debug_level)then
              open(unit=logunit,file=logfile,position='append')
              write(logunit,*)
              write(logunit,*)'check for nearest element update'
              write(logunit,*)'     slave node_id =',node_id
              write(logunit,*)'       inside_flag =',inside_flag
              write(logunit,*)'          elem row =',row
              write(logunit,*)'        element_id =',element_id
              write(logunit,*)'           face_id =',face_id
              write(logunit,*)'       checkstatus =',checkstatus
              write(logunit,*)'       best_status =',best_status
              write(logunit,*)'       normal_dist =',normal_dist
              close(unit=logunit)
            end if
c
c               check for the nearest element face centroid distance;
c               use as the default master element if needed
c
            if(distance.lt.closest_distance .or. row.eq.1)then
              closest_distance = distance
              closest_elem_row = row

              if(local_debug .and.logflag.ge.debug_level)then
                open(unit=logunit,file=logfile,position='append')
                write(logunit,*)
                write(logunit,*)'update closest element face centroid'
                write(logunit,*)'  closest_distance =',closest_distance
                write(logunit,*)'  closest_elem_row =',closest_elem_row
                close(unit=logunit)
              end if
            end if
c
c               1. normal distance check to the element face to
c               get the correct master element when there are nearby
c               perpendicular surfaces included in the tied contact
c               2. use the slave node gap_tolerance value to select only
c               master elements very near the slave node
c
            dist_gap_flag = .false.
            check1_flag = .false.
            check2_flag = .false.
            check3_flag = .false.
c
c               check that the slave node normal distance to the
c               master element is within the gap tolerance
c
            if(ABS(normal_dist).le.gap_tolerance)dist_gap_flag = .true.
c
c               check the slave node status score compared to
c               previous trial master elements; want a closer
c               master element but also a better status score;
c               check 1: keep a master element that is closer
c               with an equal or better status score;
c               check 2: keep a master element with a better
c               status score that is nearly as close
c
            if(checkstatus.ge.best_status .and.
     &          (ABS(normal_dist).lt.best_distance .or.
     &           best_elem_row.le.0))check1_flag = .true.

            if(checkstatus.gt.best_status .and.
     &          (ABS(normal_dist).le.(2.0*best_distance) .or.
     &           best_elem_row.le.0))check2_flag = .true.

            if(checkstatus.gt.best_status .and.
     &          dist_gap_flag .and. inside_flag)check3_flag = .true.
c
c               update the current best master element if the logic
c               conditions are met
c
            if(inside_flag .and. dist_gap_flag .and.
     &          (check1_flag .or. check2_flag .or. check3_flag))then

              best_elem_row = row
              best_distance = ABS(normal_dist)
              best_status = checkstatus

              if(local_debug .and.logflag.ge.debug_level)then
                open(unit=logunit,file=logfile,position='append')
                write(logunit,*)
                write(logunit,*)'update nearest element flags'
                write(logunit,*)'     slave node_id =',node_id
                write(logunit,*)'       inside_flag = ',inside_flag
                write(logunit,*)'     dist_gap_flag = ',dist_gap_flag
                write(logunit,*)'       check1_flag = ',check1_flag
                write(logunit,*)'       check2_flag = ',check2_flag
                write(logunit,*)'       check3_flag = ',check3_flag
                write(logunit,*)'     best_elem_row =',best_elem_row
                write(logunit,*)'     best_distance =',best_distance
                write(logunit,*)'       best_status =',best_status
                write(logunit,*)'  closest_elem_row =',closest_elem_row
                write(logunit,*)'  closest_distance =',closest_distance
                write(logunit,*)'       normal_dist =',normal_dist
                write(logunit,*)'     gap_tolerance =',gap_tolerance
                close(unit=logunit)
              end if
            else
              if(local_debug .and.logflag.ge.debug_level)then
                open(unit=logunit,file=logfile,position='append')
                write(logunit,*)
                write(logunit,*)'skip current trial element; no update'
                write(logunit,*)'     slave node_id =',node_id
                write(logunit,*)'       inside_flag = ',inside_flag
                write(logunit,*)'     dist_gap_flag = ',dist_gap_flag
                write(logunit,*)'       check1_flag = ',check1_flag
                write(logunit,*)'       check2_flag = ',check2_flag
                write(logunit,*)'  closest_elem_row =',closest_elem_row
                write(logunit,*)'  closest_distance =',closest_distance
                write(logunit,*)'     gap_tolerance =',gap_tolerance
                write(logunit,*)'       normal_dist =',normal_dist
                write(logunit,*)'       checkstatus =',checkstatus
                write(logunit,*)'       best_status =',best_status
                close(unit=logunit)
              end if
            end if

          else
            if(local_debug .and.logflag.ge.debug_level)then
              open(unit=logunit,file=logfile,position='append')
              write(logunit,*)
              write(logunit,*)'skip the element'
              write(logunit,*)'               row =',row
              write(logunit,*)'       normal_dist =',normal_dist
              close(unit=logunit)
            end if

          end if      ! .true.
        end do
c
c               save the nearest master element as the assigned element
c               for this slave node; the best master element is the
c               closest_elem_row, but if best_status = 0 the slave node
c               was not matched with any master elements, in that case
c               use the closest_elem_row as the default master element
c
        if(best_status.gt.0 .and. best_elem_row.gt.0)then
          tiednodeid(nd,2) = elemlist(best_elem_row,1)
          tiednodeid(nd,3) = elemlist(best_elem_row,2)

          if(logflag.ge.debug_level)then
            open(unit=logunit,file=logfile,position='append')
            write(logunit,*)
            write(logunit,*)'save slave node data for best master elem'
            write(logunit,*)'     row index, nd =',nd
            write(logunit,*)'     slave node id =',tiednodeid(nd,1)
            write(logunit,*)'     best_distance =',best_distance
            write(logunit,*)'     best_elem_row =',best_elem_row
            write(logunit,*)'     best elem ID  =',
     &                            elemlist(best_elem_row,1)
            write(logunit,*)'       best_status =',best_status
            write(logunit,*)'  closest_distance =',closest_distance
            write(logunit,*)'  closest_elem_row =',closest_elem_row
            write(logunit,*)'  closest elem ID  =',
     &                           elemlist(closest_elem_row,1)
            write(logunit,*)'    master elem id =',
     &                           tiednodeid(nd,2)
            write(logunit,*)'    master face id =',
     &                           tiednodeid(nd,3)
            close(unit=logunit)
          end if
c
c               save the nearest master element as the default
c               master element assigned to the slave node
c
        else
          if(closest_elem_row.gt.0)then
            tiednodeid(nd,2) = elemlist(closest_elem_row,1)
            tiednodeid(nd,3) = elemlist(closest_elem_row,2)
          else
            mflag = 2
            message = 'internal error: invalid array row index;'//
     &        ' closest_elem_row = 0; cannot save the default'//
     &        ' master element (tied_nodeassign).'
            goto 900
          end if

          mflag = 1
          message = 'warning: saving the nearest master element as'//
     &      ' the default for this slave node (tied_nodeassign).'
          if(logflag.ge.1)then
            open(unit=logunit,file=logfile,position='append')
            write(logunit,*)
            write(logunit,*)TRIM(message)
            write(logunit,*)'              mflag =',mflag
            write(logunit,*)'default master element for slave node'
            write(logunit,*)'      row index, nd =',nd
            write(logunit,*)'      slave node id =',tiednodeid(nd,1)
            write(logunit,*)'   closest_elem_row =',closest_elem_row
            write(logunit,*)'   closest_distance =',closest_distance
            write(logunit,*)'        best_status =',best_status
            write(logunit,*)'      gap_tolerance =',gap_tolerance
            write(logunit,*)'     master elem id =',
     &                              tiednodeid(nd,2)
            write(logunit,*)'     master face id =',
     &                              tiednodeid(nd,3)
            write(logunit,*)
            close(unit=logunit)
          end if
        end if
      end do

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'assigned tied slave node list'
        write(logunit,*)'      numtiednode =',numtiednode
        write(logunit,*)'   maxtieddatacol =',maxtieddatacol
        write(logunit,10)'row','node_id','elem_id','face_id','flag'
        do row=1,numtiednode
          write(logunit,11)row,tiednodeid(row,1:maxtieddatacol)
        end do
        close(unit=logunit)
      end if
10    format(a10,4(a10))
11    format(i10,4(i10))
c
c           jump here on an error for exit procedures
c
900   continue
c
c           deallocate local arrays
c
      if(allocated(elemlist))deallocate(elemlist)
      if(allocated(surfaceflagtable))deallocate(surfaceflagtable)
      if(allocated(masternodelist))deallocate(masternodelist)
      if(allocated(elemfacecentroid))deallocate(elemfacecentroid)

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'----- end tied_nodeassign -----'
        close(unit=logunit)
      end if

      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine tied_tetfaceisopoint              *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 11/13/02                  *
c     *                                    02/26/03 gvt              *
c     *                                                              *
c     *  compute the tetrahedron isoparametric coordinates for a     *
c     *  point on one of the tetrahedron element faces given in the  *
c     *  global x,y,z coordinate system                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine tied_tetfaceisopoint(isotetpoint,maxtetdof,
     &                  globalpoint,maxdof,constantisodof,
     &                  elemcoord,maxnodeperelem,maxisodof,
     &                  nodesperelem,faceid,
     &                  localfacenode,maxelemface,maxnodeperface,
     &                  numelemface,numnodeperface,numcornerperface,
     &                  mflag,message,logflag,logunit,logfile)
c
c           declare variables
c
      implicit none
      integer
     &  maxtetdof,maxdof,constantisodof,maxnodeperelem,maxisodof,
     &  nodesperelem,faceid,
     &  maxelemface,maxnodeperface,numelemface,numnodeperface,
     &  numcornerperface,
     &  logflag,logunit,mflag,
     &  localfacenode(maxelemface,maxnodeperface)
      double precision
     &  isotetpoint(maxtetdof),globalpoint(maxdof),
     &  elemcoord(maxnodeperelem,maxisodof)
      character(len=*) :: logfile,message
c
c           local variables
c
      integer, parameter :: dim3 = 3
      integer
     &  debug_level,corner1,corner2,corner3
      logical
     &  local_debug
      double precision
     &  s1,s2,s3,s4,face_length,point_dist,iso_tet_sum,
     &  edge_vec_1(dim3),edge_vec_2(dim3),position(dim3),
     &  normal(dim3),edgenormal(dim3),
     &  midedge1(dim3),midedge2(dim3),
     &  zero,one
      data zero, one
     & / 0.0d00, 1.0d00 /
c
c           declare functions
c
      double precision
     &  srt_dotprod
c
      debug_level = 10
      local_debug = .false.

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'----- begin tied_tetfaceisopoint -----'
        write(logunit,*)'         maxtetdof =',maxtetdof
        write(logunit,*)'            maxdof =',maxdof
        write(logunit,*)'    constantisodof =',constantisodof
        write(logunit,*)'      nodesperelem =',nodesperelem
        write(logunit,*)'            faceid =',faceid
        write(logunit,*)'       numelemface =',numelemface
        write(logunit,*)'    numnodeperface =',numnodeperface
        write(logunit,*)'  numcornerperface =',numcornerperface
        write(logunit,3)'      globalpoint() =',globalpoint(1:maxdof)
        if(local_debug)then
          write(logunit,*)'       local_debug = ',local_debug
        end if
        close(unit=logunit)
      end if
3     format(a,3es16.6)
c
c           error checking
c
      if(nodesperelem.le.0)then
        mflag = 2
        message = 'error: nodesperelem = 0 (tied_tetfaceisopoint).'
      else if(maxtetdof.ne.4)then
        mflag = 2
        message = 'error: maxtetdof /= 4 (tied_tetfaceisopoint).'
      else if(maxdof.ne.dim3)then
        mflag = 2
        message = 'error: maxdof /= 3 (tied_tetfaceisopoint).'
      else if(faceid.le.0)then
        mflag = 2
        message = 'error: faceid = 0 (tied_tetfaceisopoint).'
      else if(faceid.gt.numelemface)then
        mflag = 2
        message = 'error: faceid > numelemface (tied_tetfaceisopoint).'
      end if
      if(mflag.ge.2)goto 900
c
c           initialize values
c
      isotetpoint = zero

      corner1 = localfacenode(faceid,1)
      corner2 = localfacenode(faceid,2)
      corner3 = localfacenode(faceid,3)

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'local values:'
        write(logunit,*)'           corner1 =',corner1
        write(logunit,*)'           corner2 =',corner2
        write(logunit,*)'           corner3 =',corner3
        close(unit=logunit)
      end if
c
c           use the corner coordinates to get a vector along two edges
c           of the triangular face of the tetrahedron element;
c           normalize the edge and face normal vectors for unit vectors
c
      call srt_getunitvec(elemcoord(corner1,:),elemcoord(corner2,:),
     &                edge_vec_1,dim3,0)
      call srt_normalize(edge_vec_1,maxdof)

      call srt_getunitvec(elemcoord(corner1,:),elemcoord(corner3,:),
     &                edge_vec_2,dim3,0)
      call srt_normalize(edge_vec_2,maxdof)

      call srt_crossprod(edge_vec_1,edge_vec_2,normal,dim3)
      call srt_normalize(normal,maxdof)

      call srt_midpoint(elemcoord(corner1,:),elemcoord(corner2,:),
     &     midedge1,dim3)

      call srt_midpoint(elemcoord(corner1,:),elemcoord(corner3,:),
     &     midedge2,dim3)

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,3)'      edge_vec_1() =',edge_vec_1(1:dim3)
        write(logunit,3)'      edge_vec_2() =',edge_vec_2(1:dim3)
        write(logunit,3)'          normal() =',normal(1:dim3)
        write(logunit,3)'        midedge1() =',midedge1(1:dim3)
        write(logunit,3)'        midedge2() =',midedge2(1:dim3)
        close(unit=logunit)
      end if
c
c           use dot products to project the position vector from corner
c           point 1 to the face point along the edges of the triangular
c           face; use the edge length projection to get the tetrahedron
c           isoparametric coordinate
c
      select case(faceid)
      case(1)      ! face 1, corners 1-3-2, s4 = 0.0
c
c           get the projected length along the normal to face edge 1
c           (corner 1 to 3) for s2; use the cross product of the vector
c           along edge 1 and the face normal to get the edge normal
c           vector
c
        call srt_crossprod(normal,edge_vec_1,edgenormal,dim3)
        call srt_normalize(edgenormal,maxdof)

        call srt_getunitvec(midedge1,elemcoord(corner3,:),
     &                  position,dim3,0)
        face_length = srt_dotprod(position,edgenormal,dim3)

        call srt_getunitvec(midedge1,globalpoint,
     &                  position,dim3,0)
        point_dist = srt_dotprod(position,edgenormal,dim3)

        s2 = point_dist/face_length

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)'face 1, edge 1 (corner 1-3) and s2:'
          write(logunit,3)'       edgenormal() =',edgenormal(1:dim3)
          write(logunit,*)'        point_dist =',point_dist
          write(logunit,*)'       face_length =',face_length
          write(logunit,*)'                s2 =',s2
          close(unit=logunit)
        end if
c
c           get the projected length along the normal to face edge 2
c           (corner 1 to 2) for s3; use the cross product of the vector
c           along edge 2 and the face normal to get the edge normal
c           vector
c
        call srt_crossprod(edge_vec_2,normal,edgenormal,dim3)
        call srt_normalize(edgenormal,maxdof)

        call srt_getunitvec(midedge2,elemcoord(corner2,:),
     &                  position,dim3,0)
        face_length = srt_dotprod(position,edgenormal,dim3)

        call srt_getunitvec(midedge2,globalpoint,
     &                  position,dim3,0)
        point_dist = srt_dotprod(position,edgenormal,dim3)

        s3 = point_dist/face_length

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)'face 1, edge 2 (corner 1-2) and s3:'
          write(logunit,3)'       edgenormal() =',edgenormal(1:dim3)
          write(logunit,*)'        point_dist =',point_dist
          write(logunit,*)'       face_length =',face_length
          write(logunit,*)'                s3 =',s3
          close(unit=logunit)
        end if
c
c           update the dependent tetrahedron coordinate
c
        s1 = one - s2 - s3
        s4 = zero
        if(constantisodof.ne.4)then
          mflag = 2
          message = 'error: constantisodof /= 4 on tetrahedron'//
     &              ' face 1 (tied_tetfaceisopoint).'
        end if

      case(2)      ! face 2, corners 1-4-2, s3 = 0.0
c
c           get the projected length along the normal to face edge 1
c           (corner 1 to 2) for s4; use the cross product of the vector
c           along edge 1 and the face normal to get the edge normal
c           vector
c
        call srt_crossprod(normal,edge_vec_1,edgenormal,dim3)
        call srt_normalize(edgenormal,maxdof)

        call srt_getunitvec(midedge1,elemcoord(corner3,:),
     &                  position,dim3,0)
        face_length = srt_dotprod(position,edgenormal,dim3)

        call srt_getunitvec(midedge1,globalpoint,
     &                  position,dim3,0)
        point_dist = srt_dotprod(position,edgenormal,dim3)

        s4 = point_dist/face_length

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)'face 2, edge 1 (corner 1-2) and s4:'
          write(logunit,3)'       edgenormal() =',edgenormal(1:dim3)
          write(logunit,*)'        point_dist =',point_dist
          write(logunit,*)'       face_length =',face_length
          write(logunit,*)'                s4 =',s4
          close(unit=logunit)
        end if
c
c           get the projected length along the normal to face edge 2
c           (corner 1 to 3) for s2; use the cross product of the vector
c           along edge 2 and the face normal to get the edge normal
c           vector
c
        call srt_crossprod(edge_vec_2,normal,edgenormal,dim3)
        call srt_normalize(edgenormal,maxdof)

        call srt_getunitvec(midedge2,elemcoord(corner2,:),
     &                  position,dim3,0)
        face_length = srt_dotprod(position,edgenormal,dim3)

        call srt_getunitvec(midedge2,globalpoint,
     &                  position,dim3,0)
        point_dist = srt_dotprod(position,edgenormal,dim3)

        s2 = point_dist/face_length

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)'face 2, edge 2 (corner 1-3) and s2:'
          write(logunit,3)'       edgenormal() =',edgenormal(1:dim3)
          write(logunit,*)'        point_dist =',point_dist
          write(logunit,*)'       face_length =',face_length
          write(logunit,*)'                s2 =',s2
          close(unit=logunit)
        end if
c
c           update the dependent tetrahedron coordinate
c
        s1 = one - s2 - s4
        s3 = zero
        if(constantisodof.ne.3)then
          mflag = 2
          message = 'error: constantisodof /= 3 on tetrahedron'//
     &              ' face 2 (tied_tetfaceisopoint).'
        end if

      case(3)      ! face 3, corners 2-3-4, s1 = 0.0
c
c           get the projected length along the normal to face edge 1
c           (corner 2 to 3) for s4; use the cross product of the face
c           normal and vetor along edge 1 to get the edge normal
c           vector
c
        call srt_crossprod(normal,edge_vec_1,edgenormal,dim3)
        call srt_normalize(edgenormal,maxdof)

        call srt_getunitvec(midedge1,elemcoord(corner3,:),
     &                  position,dim3,0)
        face_length = srt_dotprod(position,edgenormal,dim3)

        call srt_getunitvec(midedge1,globalpoint,
     &                  position,dim3,0)
        point_dist = srt_dotprod(position,edgenormal,dim3)

        s4 = point_dist/face_length

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)'face 3, edge 1 (corner 2-3) and s4:'
          write(logunit,3)'       edgenormal() =',edgenormal(1:dim3)
          write(logunit,*)'        point_dist =',point_dist
          write(logunit,*)'       face_length =',face_length
          write(logunit,*)'                s4 =',s4
          close(unit=logunit)
        end if
c
c           get the projected length along the normal to face edge 2
c           (corner 2 to 4) for s3; use the cross product of the vector
c           along edge 2 and the face normal to get the edge normal
c           vector
c
        call srt_crossprod(edge_vec_2,normal,edgenormal,dim3)
        call srt_normalize(edgenormal,maxdof)

        call srt_getunitvec(midedge2,elemcoord(corner2,:),
     &                  position,dim3,0)
        face_length = srt_dotprod(position,edgenormal,dim3)

        call srt_getunitvec(midedge2,globalpoint,
     &                  position,dim3,0)
        point_dist = srt_dotprod(position,edgenormal,dim3)

        s3 = point_dist/face_length

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)'face 3, edge 2 (corner 2-4) and s3:'
          write(logunit,3)'       edgenormal() =',edgenormal(1:dim3)
          write(logunit,*)'        point_dist =',point_dist
          write(logunit,*)'       face_length =',face_length
          write(logunit,*)'                s3 =',s3
          close(unit=logunit)
        end if
c
c           update the dependent tetrahedron coordinate
c
        s2 = one - s3 - s4
        s1 = zero
        if(constantisodof.ne.1)then
          mflag = 2
          message = 'error: constantisodof /= 1 on tetrahedron'//
     &              ' face 3 (tied_tetfaceisopoint).'
        end if

      case(4)      ! face 4, corners 1-4-3, s2 = 0.0
c
c           get the projected length along the normal to face edge 1
c           (corner 1 to 4) for s3; use the cross product of the face
c           normal and vetor along edge 1 to get the edge normal
c           vector
c
        call srt_crossprod(normal,edge_vec_1,edgenormal,dim3)
        call srt_normalize(edgenormal,maxdof)

        call srt_getunitvec(midedge1,elemcoord(corner3,:),
     &                  position,dim3,0)
        face_length = srt_dotprod(position,edgenormal,dim3)

        call srt_getunitvec(midedge1,globalpoint,
     &                  position,dim3,0)
        point_dist = srt_dotprod(position,edgenormal,dim3)

        s3 = point_dist/face_length

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)'face 4, edge 1 (corner 1-4) and s3:'
          write(logunit,3)'       edgenormal() =',edgenormal(1:dim3)
          write(logunit,*)'        point_dist =',point_dist
          write(logunit,*)'       face_length =',face_length
          write(logunit,*)'                s3 =',s3
          close(unit=logunit)
        end if
c
c           get the projected length along the normal to face edge 2
c           (corner 1 to 3) for s4; use the cross product of the vector
c           along edge 2 and the face normal to get the edge normal
c           vector
c
        call srt_crossprod(edge_vec_2,normal,edgenormal,dim3)
        call srt_normalize(edgenormal,maxdof)

        call srt_getunitvec(midedge2,elemcoord(corner2,:),
     &                  position,dim3,0)
        face_length = srt_dotprod(position,edgenormal,dim3)

        call srt_getunitvec(midedge2,globalpoint,
     &                  position,dim3,0)
        point_dist = srt_dotprod(position,edgenormal,dim3)

        s4 = point_dist/face_length

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)'face 3, edge 2 (corner 1-3) and s4:'
          write(logunit,3)'       edgenormal() =',edgenormal(1:dim3)
          write(logunit,*)'        point_dist =',point_dist
          write(logunit,*)'       face_length =',face_length
          write(logunit,*)'                s4 =',s4
          close(unit=logunit)
        end if
c
c           update the dependent tetrahedron coordinate
c
        s1 = one - s3 - s4
        s2 = zero
        if(constantisodof.ne.2)then
          mflag = 2
          message = 'error: constantisodof /= 2 on tetrahedron'//
     &              ' face 4 (tied_tetfaceisopoint).'
        end if

      case default
        mflag = 2
        message = 'error: unexpected faceid (tied_tetfaceisopoint).'
      end select
      if(mflag.ge.2)goto 900

      isotetpoint(1) = s1
      isotetpoint(2) = s2
      isotetpoint(3) = s3
      isotetpoint(4) = s4

      if(local_debug)then
        iso_tet_sum = isotetpoint(1) + isotetpoint(2) +
     &                isotetpoint(3) + isotetpoint(4)

        if(logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)'                s1 =',s1
          write(logunit,*)'                s2 =',s2
          write(logunit,*)'                s3 =',s3
          write(logunit,*)'                s4 =',s4
          write(logunit,*)'       iso_tet_sum =',iso_tet_sum
          close(unit=logunit)
        end if
      end if
c
c           jump here on an error for exit procedures
c
900   continue

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'----- end tied_tetfaceisopoint -----'
        close(unit=logunit)
      end if

      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine tied_tetfacetangent               *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 11/13/02                  *
c     *                                    02/26/03 gvt              *
c     *                                                              *
c     *  use the triangular face corners to get orthogonal tangent   *
c     *  vectors on the tetrahedron element face                     *
c     *                                                              *
c     ****************************************************************
c
      subroutine tied_tetfacetangent(tangent_a,tangent_b,normal,maxdof,
     &                  elemcoord,maxnodeperelem,maxisodof,
     &                  nodesperelem,faceid,
     &                  localfacenode,maxelemface,maxnodeperface,
     &                  numelemface,numnodeperface,numcornerperface,
     &                  mflag,message,logflag,logunit,logfile)
c
c           declare variables
c
      implicit none
      integer
     &  maxdof,maxnodeperelem,maxisodof,nodesperelem,faceid,
     &  maxelemface,maxnodeperface,numelemface,numnodeperface,
     &  numcornerperface,
     &  logflag,logunit,mflag,
     &  localfacenode(maxelemface,maxnodeperface)
      double precision
     &  tangent_a(maxdof),tangent_b(maxdof),normal(maxdof),
     &  elemcoord(maxnodeperelem,maxisodof)
      character(len=*) :: logfile,message
c
c           local variables
c
      integer, parameter :: dim3 = 3
      integer
     &  debug_level,corner1,corner2,corner3
      logical
     &  local_debug
      double precision
     &  zero,edge_vec_1(dim3),edge_vec_2(dim3)
      data zero
     & / 0.0d00 /
c
      debug_level = 10
      local_debug = .false.

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'----- begin tied_tetfacetangent -----'
        write(logunit,*)'            maxdof =',maxdof
        write(logunit,*)'      nodesperelem =',nodesperelem
        write(logunit,*)'            faceid =',faceid
        write(logunit,*)'       numelemface =',numelemface
        write(logunit,*)'    numnodeperface =',numnodeperface
        write(logunit,*)'  numcornerperface =',numcornerperface
        if(local_debug)then
          write(logunit,*)'       local_debug = ',local_debug
        end if
        close(unit=logunit)
      end if
c
c           error checking
c
      if(nodesperelem.le.0)then
        mflag = 2
        message = 'error: nodesperelem = 0 (tied_tetfacetangent).'
      else if(maxdof.ne.dim3)then
        mflag = 2
        message = 'error: maxdof /= 3 (tied_tetfacetangent).'
      else if(faceid.le.0)then
        mflag = 2
        message = 'error: faceid = 0 (tied_tetfacetangent).'
      else if(faceid.gt.numelemface)then
        mflag = 2
        message = 'error: faceid > numelemface (tied_tetfacetangent).'
      end if
      if(mflag.ge.2)goto 900
c
c           initialize values
c
      tangent_a = zero
      tangent_b = zero
      normal = zero

      corner1 = localfacenode(faceid,1)
      corner2 = localfacenode(faceid,2)
      corner3 = localfacenode(faceid,3)

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'local values:'
        write(logunit,*)'           corner1 =',corner1
        write(logunit,*)'           corner2 =',corner2
        write(logunit,*)'           corner3 =',corner3
        close(unit=logunit)
      end if
c
c           use the corner coordinates to get a vector along two edges
c           of the triangular face of the tetrahedron element;
c           the face outward normal is the cross product of the edge
c           vectors; the first tangent vector is along edge 1 from
c           corners 1 to 2; the orthogonal tangent is the cross product
c           of the outward normal and the first tangent vector
c
      call srt_getunitvec(elemcoord(corner1,:),elemcoord(corner2,:),
     &                edge_vec_1,dim3,0)

      call srt_getunitvec(elemcoord(corner1,:),elemcoord(corner3,:),
     &                edge_vec_2,dim3,0)

      call srt_crossprod(edge_vec_1,edge_vec_2,normal,dim3)

      call srt_crossprod(normal,edge_vec_1,tangent_b,dim3)

      tangent_a(1:maxdof) = edge_vec_1(1:dim3)

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,3)'      edge_vec_1() =',edge_vec_1(1:dim3)
        write(logunit,3)'      edge_vec_2() =',edge_vec_2(1:dim3)
        write(logunit,3)'          normal() =',normal(1:maxdof)
        write(logunit,3)'       tangent_a() =',tangent_a(1:maxdof)
        write(logunit,3)'       tangent_b() =',tangent_b(1:maxdof)
        close(unit=logunit)
      end if
3     format(a,3es16.6)
c
c           jump here on an error for exit procedures
c
900   continue

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'----- end tied_tetfacetangent -----'
        close(unit=logunit)
      end if

      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine tied_isoconstdof                *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 11/06/02                  *
c     *                                    02/26/03 gvt              *
c     *                                                              *
c     *  set the isoparametric dof number and coordinate value that  *
c     *  is constant on the given element face                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine tied_isoconstdof(constantisodof,constantisocoord,
     &              faceid,nodesperelem,elemchoice,meshformat,
     &              logflag,logunit,logfile,mflag,message)
c
c           declare variables
c
      implicit none
      integer
     &  constantisodof,faceid,nodesperelem,elemchoice,meshformat,
     &  logflag,logunit,mflag
      double precision
     &  constantisocoord
      character(len=*) :: logfile,message
c
c           local variables
c
      integer, parameter :: maxelemface = 6
      integer, parameter :: maxdataformat = 2
      integer
     &  debug_level,i,brickisodof(maxelemface,maxdataformat),
     &  tetisodof(maxelemface,maxdataformat)
      double precision
     &  zero,one,brickisocoord(maxelemface,maxdataformat)
      data zero, one
     &   / 0.0d00, 1.0d00 /
c
c           set the isoparametric dof and coordinate value for each
c           face of a brick or tetrahedron element for either
c           abaqus (column 1) or warp3d (column 2) input data format
c
      data brickisodof(1:maxelemface,1) / 1,1,2,3,2,3 /
      data brickisodof(1:maxelemface,2) / 1,1,2,2,3,3 /
      data brickisocoord(1:maxelemface,1)
     &  /-1.0d0,1.0d0,-1.0d0,-1.0d0,1.0d0,1.0d0/
      data brickisocoord(1:maxelemface,2)
     &  /-1.0d0,1.0d0,-1.0d0,1.0d0,-1.0d0,1.0d0/

      data tetisodof(1:maxelemface,1) / 4,3,1,2,0,0 /
      data tetisodof(1:maxelemface,2) / 4,3,1,2,0,0 /
c
      debug_level = 10

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'----- begin tied_isoconstdof -----'
        write(logunit,*)'            faceid =',faceid
        write(logunit,*)'      nodesperelem =',nodesperelem
        write(logunit,*)'        elemchoice =',elemchoice
        write(logunit,*)'        meshformat =',meshformat
        write(logunit,*)'local values:'
        write(logunit,*)'       maxelemface =',maxelemface
        write(logunit,*)'     maxdataformat =',maxdataformat
        write(logunit,*)'  face,  brickisodof()'
        do i=1,maxelemface
          write(logunit,*)i,brickisodof(i,1:maxdataformat)
        end do
        write(logunit,*)'  face,  brickisocoord()'
        do i=1,maxelemface
          write(logunit,*)i,brickisocoord(i,1:maxdataformat)
        end do
        write(logunit,*)'  face,  tetisodof()'
        do i=1,maxelemface
          write(logunit,*)i,tetisodof(i,1:maxdataformat)
        end do
        close(unit=logunit)
      end if
c
c           error checking
c
      if(nodesperelem.le.0)then
        mflag = 2
        message = 'error: nodesperelem = 0 (tied_isoconstdof).'
      else if(meshformat.le.0 .or. meshformat.gt.maxdataformat)then
        mflag = 2
        message = 'error: invalid meshformat; expecting'//
     &    ' meshformat = 1,2 (tied_isoconstdof).'
      end if
      if(mflag.ge.2)goto 900
c
c           set the constant isoparametric dof value
c
c           for a brick element the isoparametric dof 1=xi, 2=eta, 3=zeta
c           and on the brick faces the constant value is +/- 1.0
c
c           abaqus or warp3d format for the input data
c
      select case(elemchoice)
      case(1,2)      ! 20 or 8 node brick
        select case(faceid)
        case(1:6)
          constantisodof = brickisodof(faceid,meshformat)
          constantisocoord = brickisocoord(faceid,meshformat)

        case default
          mflag = 2
          select case(meshformat)
          case(1)      ! 1=abaqus
            message = 'error: unexpected element face id for a'//
     &                ' brick element; cannot set the element face'//
     &                ' constant isoparametric values; abaqus format'//
     &                ' (tied_isoconstdof).'
          case(2)      ! 2=warp3d
            message = 'error: unexpected element face id for a'//
     &                ' brick element; cannot set the element face'//
     &                ' constant isoparametric values; warp3d format'//
     &                ' (tied_isoconstdof).'
          end select
          goto 900
        end select

      case(6,13)      ! 10 or 4 node tetrahedron
c
c           note that for a tetrahedron there are 4 isoparametric
c           coordinates, and s1+s2+s3+s4=1 or s1=1-s2-s3-s4, so only
c           3 of the 4 values are independent, one must be dependent;
c           the isoparametric coordinate on a tetrahedron face is
c           equal to zero
c
        select case(faceid)
        case(1:4)
          constantisodof = tetisodof(faceid,meshformat)
          constantisocoord = zero

        case default
          mflag = 2
          select case(meshformat)
          case(1)      ! 1=abaqus
            message = 'error: unexpected element face id for a'//
     &                ' tetrahedron element; cannot set the element'//
     &                ' face constant isoparametric values;'//
     &                ' abaqus input format (tied_isoconstdof).'
          case(2)      ! 2=warp3d
            message = 'error: unexpected element face id for a'//
     &                ' tetrahedron element; cannot set the element'//
     &                ' face constant isoparametric values;'//
     &                ' warp3d input format (tied_isoconstdof).'
          end select
          goto 900
        end select

      case default
        mflag = 2
          message = 'error: unexpected element choice; cannot set'//
     &              ' the element face constant isoparametric values;'//
     &              ' check "elemchoice"; abaqus element format'//
     &              ' (tied_isoconstdof).'
        goto 900
      end select
c
c           jump here on an error for exit procedures
c
900   continue

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'    constantisodof =',constantisodof
        write(logunit,*)'  constantisocoord =',constantisocoord
        write(logunit,*)'----- end tied_isoconstdof -----'
        close(unit=logunit)
      end if

      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *              subroutine tied_globalcoordfromiso              *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                   last modified : 11/06/02 gvt               *
c     *                                   02/26/03 gvt               *
c     *                                                              *
c     *  compute the face isoparametric point in the global x,y,z    *
c     *  coordinates; use the element shape functions evaluated      *
c     *  at the isoparametric point                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine tied_globalcoordfromiso(globalpoint,maxdof,
     &                shapefunc,maxnodeperelem,nodesperelem,
     &                elemcoord,maxelemcoordrow,maxelemcoordcol,
     &                logflag,logunit,logfile,mflag,message)
c
c           declare variables
c
      implicit none
      integer
     &  maxdof,maxnodeperelem,nodesperelem,maxelemcoordrow,
     &  maxelemcoordcol,
     &  logflag,logunit,mflag

      double precision
     &  globalpoint(maxdof),shapefunc(maxnodeperelem),
     &  elemcoord(maxelemcoordrow,maxelemcoordcol)

      character(len=*) ::
     &  logfile,message
c
c           local variables
c
      integer
     &  debug_level,i,j

      double precision
     &  zero
      data zero
     & / 0.0d00 /
      debug_level = 9

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'----- begin tied_globalcoordfromiso -----'
        write(logunit,*)'               maxdof =',maxdof
        write(logunit,*)'       maxnodeperelem =',maxnodeperelem
        write(logunit,*)'         nodesperelem =',nodesperelem
        write(logunit,*)'      maxelemcoordrow =',maxelemcoordrow
        write(logunit,*)'      maxelemcoordcol =',maxelemcoordcol
        write(logunit,*)'      maxelemcoordcol =',maxelemcoordcol
        close(unit=logunit)
      end if
c
c           error checking
c
      if(maxdof.le.0)then
        mflag = 2
        message = 'error: maxdof = 0'//
     &            ' (tied_globalcoordfromiso).'
      else if(nodesperelem.le.0)then
        mflag = 2
        message = 'error: nodesperelem = 0'//
     &            ' (tied_globalcoordfromiso).'
      end if
      if(mflag.ge.2)goto 900
c
c           initialize values
c
      globalpoint = zero
c
c           sum the shape functions at each element node to compute
c           the global x,y,z coordinates at the given isoparametric
c           point (the shapefunc() array should already be evaluated at
c           the isoparametric point)
c
      do i=1,nodesperelem
        do j=1,maxdof
          globalpoint(j) = globalpoint(j) +
     &                         shapefunc(i)*elemcoord(i,j)
        end do
      end do

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,1)'         globalpoint() =',globalpoint(1:maxdof)
        close(unit=logunit)
      end if
1     format(a,3es16.6)
c
c           jump here on an error for exit procedures
c
900   continue

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'----- end tied_globalcoordfromiso -----'
        close(unit=logunit)
      end if

      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *            subroutine tied_exportmpc_warp3dformat            *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 02/26/03                  *
c     *                                                              *
c     *  export the computed multi point constraint (mpc) equations  *
c     *  for tied contact in the warp3d mpc equation format; use the *
c     *  mpc equations in the input file for testing the tied        *
c     *  contact                                                     *
c     *                                                              *
c     ****************************************************************
c
      subroutine tied_exportmpc_warp3dformat(epsilon,io)
c
c           declare modules
c
      use mod_mpc
c
c           declare variables
c
      implicit none
      integer
     &  io

      double precision
     &  epsilon
c
c           variables:
c           epsilon = tolerance for the mpc equation coefficient
c             value; write the mpc coefficient if it is larger
c             than epsilon
c           io = file unit number to write the mpc export to;
c             the file should already be open and connected to
c             the file unit number
c
c
c           local variables
c
      integer
     &  i,row,numlist,column,node_id,local_dof
      double precision
     &  zero,mult
      data zero
     & / 0.0d00 /
      character xyz_letter(3)*1
c
c           initialize values
c
      xyz_letter(1) = 'u'
      xyz_letter(2) = 'v'
      xyz_letter(3) = 'w'
c
c           write the tied contact mpc equations in warp3d format
c           to support testing of the tied contact and user defined
c           mpc in warp3d; the listed user mpc equations can replace
c           the corresponding tied contact data in the warp3d input
c           file
c
c           warp3d format is:
c           <nodes> <multiplier> <u,v,w> + (,) ... + (,) = <constant>
c
c           for example:
c           6 1.5 u + 6 3.2 v + 6 -3.1 w = 0.0
c
      write(io,*)
      write(io,*)
      write(io,9)'c'
      write(io,9)'c  tied node contact processing'
      write(io,9)'c  export mpc equations in warp3d format'
      write(io,9)'c'
      write(io,12)'c          epsilon =',epsilon
      do i=1,3
        write(io,10)'c    xyz_letter(',i,') = ',xyz_letter(i)
      end do
      write(io,9)'c'
      write(io,11)'c  number of mpc equations =',num_tied_con_mpc
      write(io,9)'c'
      write(io,9)'constraints'

      do row=1,num_tied_con_mpc
        numlist = tied_con_mpc_table(row)%num_terms
        write(io,*)
        do column=1,numlist
          if(abs(tied_con_mpc_table(row)%multiplier_list(column))
     &    .gt.epsilon)then

            node_id = tied_con_mpc_table(row)%node_list(column)
            mult = tied_con_mpc_table(row)%multiplier_list(column)
            local_dof = tied_con_mpc_table(row)%dof_list(column)

            if(column.gt.1)then
              write(io,13)' + ',node_id,mult,
     &                          xyz_letter(local_dof),','
            else
              write(io,13)'   ',node_id,mult,
     &                          xyz_letter(local_dof),','
            end if
          end if
        end do

        write(io,14)' = ',zero
      end do

      write(io,9)'c'
      write(io,9)'c  end of tied node mpc equations'
      write(io,9)'c'
      write(io,*)
      write(io,*)

9     format(a)
10    format(a,i1,a,a)
11    format(a,i8)
12    format(a,es16.6)
13    format(8x,a,i7,es16.7,a2,a)
14    format(8x,a,f15.4,a2,a)

      return
      end




c
c     ****************************************************************
c     *                                                              *
c     *            subroutine tied_exportmpc_abaqusformat            *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 02/26/03                  *
c     *                                                              *
c     *  export the computed multi point constraint (mpc) equations  *
c     *  for tied contact in the abaqus *equation format; use the    *
c     *  mpc equations in the input file for testing the tied        *
c     *  contact                                                     *
c     *                                                              *
c     ****************************************************************
c
      subroutine tied_exportmpc_abaqusformat(epsilon,io)
c
c           declare modules
c
      use mod_mpc
c
c           declare variables
c
      implicit none
      integer
     &  io

      double precision
     &  epsilon
c
c           variables:
c           epsilon = tolerance for the mpc equation coefficient
c             value; write the mpc coefficient if it is larger
c             than epsilon
c           io = file unit number to write the mpc export to;
c             the file should already be open and connected to
c             the file unit number
c
c
c           local variables
c
      integer
     &  i,row,count,numlist,column,count_terms,node_id,local_dof
      double precision
     &  zero, mult
      data zero
     &       / 0.0d00 /
      character text*512,word*10,xyz_letter(3)*1
c
c           initialize values
c
      xyz_letter(1) = '1'
      xyz_letter(2) = '2'
      xyz_letter(3) = '3'
c
c           write the tied contact mpc equations in abaqus format
c           to support testing of the tied contact and user defined
c           mpcs; the listed user mpc equations can replace
c           the corresponding tied contact data in the abaqus input
c           file; note that the first node is removed as a degree of
c           freedom to impose the constraint, don't apply any other
c           boundary conditions or constraints to that node;
c           the slave node should already be first in the equation list;
c           the data line can have a maximum of four node terms per
c           line, use extra lines for more terms in the equation
c
c           abaqus format is:
c           line 1: n (the number of nodes in the equation)
c           line 2: node, dof, mult, ... (terms for each node)
c
c           for example:
c           3
c           6, 1, 1.5, 6, 2, 3.2, 6, 3, -3.1
c           (node 6, x, 1.5, node 6, y, 3.2, node 6, z, -3.1 = 0.0)
c
      write(io,*)
      write(io,*)
      write(io,9)'**'
      write(io,9)'**  tied node contact processing'
      write(io,9)'**  export mpc equations in abaqus format'
      write(io,9)'**'
      write(io,12)'**          epsilon =',epsilon
      do i=1,3
        write(io,10)'**    xyz_letter(',i,') = ',xyz_letter(i)
      end do
      write(io,9)'**'
      write(io,11)'**  number of mpc equations =',num_tied_con_mpc
      write(io,9)'**'

      write(io,9)'*equation'

      do row=1,num_tied_con_mpc
        count = 0
        numlist = tied_con_mpc_table(row)%num_terms
c
c           count the number of nodes in the mpc equation
c
        do column=1,numlist
          if(abs(tied_con_mpc_table(row)%multiplier_list(column))
     &    .gt.epsilon)then
            count = count + 1
          end if
        end do

        write(io,13)count
c
c           get the mpc equation terms
c
        count_terms = 0
        text = ' '

        do column=1,numlist
          if(abs(tied_con_mpc_table(row)%multiplier_list(column))
     &    .gt.epsilon)then

            count_terms = count_terms + 1
            node_id = tied_con_mpc_table(row)%node_list(column)
            local_dof = tied_con_mpc_table(row)%dof_list(column)
            mult = tied_con_mpc_table(row)%multiplier_list(column)

            write(word,31)node_id
            text = trim(text)//trim(adjustl(word))//', '
            text = trim(text)//xyz_letter(local_dof)//', '
            write(word,32)mult
            text = trim(text)//trim(adjustl(word))//', '
c
c                 maximum of 4 terms in the mpc equation per data line
c
            if(count_terms.eq.4)then
              write(io,9)trim(text)
              count_terms = 0
              text = ' '
            end if
          end if
        end do
c
c                 write remaining terms in the mpc equation to the data line
c
        if(count_terms.gt.0)then
          write(io,9)trim(text)
        end if
      end do

      write(io,9)'**'
      write(io,9)'**  end of tied node mpc equations'
      write(io,9)'**'
      write(io,*)
      write(io,*)

9     format(a)
10    format(a,i1,a,a)
11    format(a,i8)
12    format(a,es16.6)
13    format(i4)
31    format(i10)
32    format(f10.5)

      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *           subroutine tied_exportmpc_nodelist                 *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 06/09/04                  *
c     *                                                              *
c     *  write a table of the slave node and its assigned master     *
c     *  element for all the slave nodes in the tied contact         *
c     *                                                              *
c     ****************************************************************
c
      subroutine tied_exportmpc_nodelist(tiednodeid,
     &         maxtiednode,maxtieddatacol,numtiednode,io)
c
c           declare variables
c
      implicit none
      integer
     &  maxtiednode,maxtieddatacol,numtiednode,
     &  io,
     &  tiednodeid(maxtiednode,maxtieddatacol)
c
c           variables:
c
c           tiednodeid() = list the slave node id in column 1, list the
c             master element id that the slave node is connected to in
c             column 2, list the master element face number in column 3,
c             give the node tied flag 1=yes, 0=no in column 4, if the
c             node is within the gap tolerance distance then set the
c             flag = 1 = yes to be included in the tied contact
c           maxtiednode = array row dimension; the maximum number of
c             tied slave nodes
c           maxtieddatacol = array column dimension
c           numtiednode = the number of tied slave nodes
c
c
c           local variables
c
      integer row
c
c           write a table of the slave node and the assigned master
c           element, element face number, and the yes/no tied flag
c
      write(io,*)
      write(io,*)
      write(io,*)
      write(io,2)'      >> tied contact slave node table: list the'
      write(io,2)'      >> slave node and the master element it is'
      write(io,2)'      >> assigned to'
      write(io,*)
      write(io,3)'      total number of slave nodes =',numtiednode
      write(io,*)
      write(io,1)'row','slave_node_id','master_elem_id',
     &           'elem_face_id','tied_flag'
      do row=1,numtiednode
        if(tiednodeid(row,maxtieddatacol).eq.1)then
          write(io,4)row,tiednodeid(row,1:3),'yes'
        else
          write(io,4)row,tiednodeid(row,1:3),'no'
        end if
      end do
      write(io,*)

1     format(a8,3a16,a12)
2     format(a)
3     format(a,i10)
4     format(i8,3i16,a12)

      return
      end
