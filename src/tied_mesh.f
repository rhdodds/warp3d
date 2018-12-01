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
      integer i,debug_level

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
      tiednodegaptolerance = zero

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
c     *                    subroutine tied_nodeiso                   *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 11/04/02                  *
c     *                                    02/26/03 gvt              *
c     *                                                              *
c     *  compute the isoparametric coordinates of the slave nodes    *
c     *  (the dependent, constrained nodes) on its assigned master   *
c     *  element face (the independent nodes);                       *
c     *  the slave node global x,y,z coordinates can also be updated *
c     *  to lie exactly on the master element face.                  *
c     *                                                              *
c     *  added the warning message output to the default "out" file  *
c     *  when the slave node is outside the given gap distance       *
c     *  (3/03/03 gvt)                                               *
c     *                                                              *
c     ****************************************************************
c
      subroutine tied_nodeiso(tiednodeid,tiednodeisocoord,
     &  tiednodeglobalcoord,tiednodeadjustflag,tiednodegaptolerance,
     &  maxtiednode,maxtieddatacol,maxisodof,numtiednode,
     &  elemcon,dimelem,maxnodeperelem,numelem,
     &  nodecoord,dimnode,dimcoord,numnode,
     &  meshformat,props,out,
     &  logflag,logunit,logfile,mflag,message)
c
c           declare modules
c
      implicit none
      include 'param_def'
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
c           local variables
c
      integer, parameter :: maxdof = 3
      integer, parameter :: maxtetdof = 4
      integer, parameter :: maxelemface = 6
      integer, parameter :: maxnodeperface = 9
      integer, parameter :: maxiterloop = 1000
      integer
     &  debug_level,allocate_error,i,j,nd,row,
     &  numelemface,numnodeperface,numcornerperface,
     &  node_id,elem_id,face_id,node_per_elem,elemchoice,
     &  constantisodof,itercount,node_row,corner,
     &  min_iter_count,max_iter_count,min_iter_nodeid,max_iter_nodeid,
     &  xi_error_count,eta_error_count,zeta_error_count,
     &  check_iso_coord_flag,ierr,iword,err,gap_cnt,cnt,dumi,
     &  localfacenode(maxelemface,maxnodeperface)
      integer, allocatable, dimension (:) :: gap_nodes,gap_elems,
     &                                       gap_faces

      logical
     &  local_debug,twod_flag,iso_adjust_flag

      double precision
     &  zero,one,three,pi,
     &  constantisocoord,dist_point,epsilon,angle_point,
     &  epsilon_angle,delta1,delta2,delta_iso,scale1,scale2,
     &  iso_error_val,dist_gap,
     &  isopoint(maxdof),globalpoint(maxdof),position(maxdof),
     &  normal(maxdof),tangent1(maxdof),tangent2(maxdof),
     &  jacobian(maxdof,maxdof),invj(maxdof,maxdof),detj,
     &  isopoint_old(maxdof),slavenode(maxdof),vec(maxdof),
     &  isotetpoint(maxtetdof),vec1(maxdof),vec2(maxdof),
     &  dsf(32,3),coord(3,32),
     &  shapefunc,shapederivs,elemcoord
      double precision,
     &  allocatable, dimension (:) :: gap_dists, gap_toler
      real rword,dumr
      double precision  dumd
      character(len=1) :: dums
      data zero, one, three
     & / 0.0d00, 1.0d00, 3.0d00 /
      allocatable
     &  shapefunc(:),shapederivs(:,:),elemcoord(:,:)
c
c           declare functions
c
      double precision
     &  srt_dist3d,srt_getangle,srt_dotprod

      equivalence (rword,iword)
c
      debug_level = 7
      local_debug = .false.

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'----- begin tied_nodeiso -----'
        write(logunit,*)'          maxtiednode =',maxtiednode
        write(logunit,*)'       maxtieddatacol =',maxtieddatacol
        write(logunit,*)'            maxisodof =',maxisodof
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
        message = 'error: numnode = 0 (tied_nodeiso).'
      else if(numnode.gt.dimnode)then
        mflag = 2
        message = 'error: numnode > dimnode (tied_nodeiso).'
      else if(dimcoord.le.0)then
        mflag = 2
        message = 'error: dimcoord = 0 (tied_nodeiso).'
      else if(dimcoord.lt.maxdof)then
        mflag = 2
        message = 'error: dimcoord < 3 (tied_nodeiso).'
      else if(numelem.le.0)then
        mflag = 2
        message = 'error: numelem = 0 (tied_nodeiso).'
      else if(numelem.gt.dimelem)then
        mflag = 2
        message = 'error: numelem > dimelem (tied_nodeiso).'
      else if(maxtiednode.le.0)then
        mflag = 2
        message = 'error: maxtiednode = 0 (tied_nodeiso).'
      else if(numtiednode.le.0)then
        mflag = 2
        message = 'error: numtiednode = 0 (tied_nodeiso).'
      else if(maxtieddatacol.le.0)then
        mflag = 2
        message = 'error: maxtieddatacol = 0 (tied_nodeiso).'
      else if(maxtieddatacol.lt.4)then
        mflag = 2
        message = 'error: maxtieddatacol < 3 (tied_nodeiso).'
      else if(maxisodof.le.0)then
        mflag = 2
        message = 'error: maxisodof = 0 (tied_nodeiso).'
      else if(maxisodof.ne.maxdof)then
        mflag = 2
        message = 'error: maxisodof /= 3 (tied_nodeiso).'
      end if
      if(mflag.ge.2)goto 900
c
c           initialize values
c
      pi = 4.0d0*atan(1.0d0)
      epsilon = 1.0d-6
      epsilon_angle = 0.001d0*(pi/180.0d0)      ! 0.001 degrees
c
c           set the check_iso_coord_flag value:
c             0=no check or change to the isopoint() values
c             1=check that isopoint() values are within +/- 1.0
c             2=set maximum isopoint() values to be within +/- 1.0
c
      check_iso_coord_flag = 2
      select case(check_iso_coord_flag)
      case(1)      ! check xi,eta,zeta and give a warning
        iso_error_val = 1.001d0
      case(2)      ! set xi,eta,zeta to have a maximum of +/- 1.0
        iso_error_val = one
      case default
        iso_error_val = 1.001d0
      end select

      min_iter_count = maxiterloop
      max_iter_count = 0
      min_iter_nodeid = 0
      max_iter_nodeid = 0
      xi_error_count = 0
      eta_error_count = 0
      zeta_error_count = 0
      gap_cnt = 0

      allocate(shapefunc(maxnodeperelem),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "shapefunc()" array'//
     &    ' (tied_nodeiso).'
        goto 900
      end if

      allocate(shapederivs(maxnodeperelem,maxisodof),
     &         stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "shapederivs()" array'//
     &    ' (tied_nodeiso).'
        goto 900
      end if

      allocate(elemcoord(maxnodeperelem,maxisodof),
     &         stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "elemcoord()" array'//
     &    ' (tied_nodeiso).'
        goto 900
      end if

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'local values:'
        write(logunit,*)'                   pi =',pi
        write(logunit,*)'              epsilon =',epsilon
        write(logunit,*)'        epsilon_angle =',epsilon_angle
        write(logunit,*)'        iso_error_val =',iso_error_val
        write(logunit,*)' check_iso_coord_flag =',check_iso_coord_flag
        close(unit=logunit)
      end if
c
c           locate the slave node on its assigned master element face;
c           compute the isoparametric element coordinates at the slave
c           node location
c
      do nd=1,numtiednode
        node_id = tiednodeid(nd,1)      ! slave node id
        elem_id = tiednodeid(nd,2)      ! assigned master element id
        face_id = tiednodeid(nd,3)      ! assigned master element face

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)
          write(logunit,*)'  slave node row index, nd =',nd
          write(logunit,*)'             slave node_id =',node_id
          write(logunit,*)'            master elem_id =',elem_id
          write(logunit,*)'    master element face_id =',face_id
          close(unit=logunit)
        end if

        if(node_id.le.0)then
          mflag = 2
          message = 'error: node_id = 0; invalid node id for the'//
     &              ' slave node; check the tiednodeid() array'//
     &              ' (tied_nodeiso).'
        else if(node_id.gt.numnode)then
          mflag = 2
          message = 'error: node_id > numnode; invalid node id for'//
     &              ' the slave node; check the tiednodeid() array'//
     &              ' (tied_nodeiso).'
        else if(elem_id.le.0)then
          mflag = 2
          message = 'error: elem_id = 0; invalid element id for'//
     &              ' the master element; check the tiednodeid()'//
     &              ' array (tied_nodeiso).'
        else if(elem_id.gt.numelem)then
          mflag = 2
          message = 'error: elem_id > numelem; invalid element id for'//
     &              ' the master element; check the tiednodeid()'//
     &              ' array (tied_nodeiso).'
        else if(face_id.le.0)then
          mflag = 2
          message = 'error: face_id = 0; invalid element face id for'//
     &              ' the master element; check the tiednodeid()'//
     &              ' array (tied_nodeiso).'
        end if
        if(mflag.ge.2)goto 900

            rword = props(1,elem_id)
            elemchoice = iword

            rword = props(2,elem_id)
            node_per_elem = iword
c
c           get the local element face node numbers (the columns
c           in the connectivity array) for this element face
c
        select case(meshformat)
        case(1)      ! abaqus
          call srt_setlocalfacenodes(localfacenode,maxelemface,
     &           maxnodeperface,numelemface,numnodeperface,
     &           numcornerperface,node_per_elem,elemchoice,
     &           logflag,logunit,logfile,mflag,message)
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
     &              ' (tied_nodeiso).'
        end select      ! meshformat
        if(mflag.ge.2)goto 900
c
c               set the isoparametric dof that is constant on the master
c               element face; set the constant isoparametric value
c
        call tied_isoconstdof(constantisodof,constantisocoord,
     &          face_id,node_per_elem,elemchoice,meshformat,
     &          logflag,logunit,logfile,mflag,message)
        if(mflag.ge.2)goto 900
c
c               set the initial face isoparametric point location at the
c               center of the master element face; at the center of a
c               brick element two of the isoparametric coordinates are
c               zero and the constant coordinate is +/- 1.0;
c               at the center of a tetrahedron element face one of the
c               coordinates is zero and the other three are 1/3
c
        select case(elemchoice)
        case(1,2)      ! 20 or 8 node brick element
          isopoint = zero
          isopoint(constantisodof) = constantisocoord
          twod_flag = .false.

        case(6,13)      ! 10 or 4 node tetrahedron
          isotetpoint = one/three
          isotetpoint(constantisodof) = zero
          twod_flag = .false.
c
c               save s2,s3,s4 as the independent isoparametric points on
c               the tetrahedron face, use isopoint() when evaluating the
c               shape functions
c
          do i=1,maxdof
            isopoint(i) = isotetpoint(i+1)
          end do

        case default
          mflag = 2
          message = 'error: unexpected elemchoice; cannot set the'//
     &              ' initial isoparametric point location'//
     &              ' (tied_nodeiso).'
          goto 900
        end select
c
c               get the slave node global x,y,z coordinates
c
        slavenode(1:maxdof) = nodecoord(node_id,1:maxdof)
c
c               put the element node coordinates into the local
c               elemcoord() array; used to compute the jacobian matrix
c
        do i=1,node_per_elem
          node_row = elemcon(elem_id,i)
          if(node_row.le.0 .or. node_row.gt.numnode)then
            mflag = 2
            message = 'error: invalid node id from the element'//
     &                ' connectivity; cannot get the node '//
     &                ' coordinates of the element (tied_nodeiso).'
            goto 900
          end if
          elemcoord(i,1:maxdof) = nodecoord(node_row,1:maxdof)
        end do

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)'     node_per_elem =',node_per_elem
          write(logunit,*)'        elemchoice =',elemchoice
          write(logunit,*)'    constantisodof =',constantisodof
          write(logunit,3)'        isopoint() =',isopoint(1:maxdof)
          write(logunit,3)'       slavenode() =',slavenode(1:maxdof)
          close(unit=logunit)
        end if
c
c               iteration loop to project the slave node onto the master
c               element face and compute the local face isoparametric
c               coordinates of the slave node's position
c
        isopoint_old = zero

        loop_1: do itercount=1,maxiterloop
          if(local_debug .and. logflag.ge.debug_level)then
            open(unit=logunit,file=logfile,position='append')
            write(logunit,*)
            write(logunit,*)'    itercount =',itercount
            write(logunit,3)'     isopoint() =',isopoint(1:maxdof)
            write(logunit,3)'    slavenode() =',slavenode(1:maxdof)
            close(unit=logunit)
          end if
c
c               compute the face isoparametric point in the global
c               x,y,z coordinates; use the element shape functions
c               evaluated at the isoparametric point
c
          call shapef(elemchoice,isopoint(1),isopoint(2),isopoint(3),
     &                shapefunc)

          call tied_globalcoordfromiso(globalpoint,maxdof,
     &                shapefunc,maxnodeperelem,node_per_elem,
     &                elemcoord,maxnodeperelem,maxisodof,
     &                logflag,logunit,logfile,mflag,message)
          if(mflag.ge.2)goto 900
c
c               get the position vector from the face isoparametric
c               point to the slave node; use the global x,y,z position
c               of the face point to compute the position vector
c
          call srt_getunitvec(globalpoint,slavenode,position,maxdof,0)
c
c               get the unit normal and unit tangent vectors to the
c               master element face at the face isoparametric point;
c               evaluate the element shape functions at the current
c               isoparametric coordinates and compute the jacobian
c               matrix
c
          select case(elemchoice)
          case(1,2)      ! 20 or 8 node brick element
c
c               evaluate the element shape functions at the current
c               isoparametric coordinates and compute the jacobian
c               matrix
c
            call derivs(elemchoice,isopoint(1),isopoint(2),isopoint(3),
     &                  shapederivs(:,1),shapederivs(:,2),
     &                  shapederivs(:,3))

            ierr = 0
            do i=1,node_per_elem
              do j=1,3
                dsf(i,j) = shapederivs(i,j)
                coord(j,i) = elemcoord(i,j)
              end do
            end do

            call eqldjb(dsf,coord,node_per_elem,jacobian,invj,detj,ierr)

            if(ierr.gt.0)then
              mflag = 2
              message = 'internal error: error while computing the'//
     &                  ' jacobian matrix; ierr > 0 (tied_nodeiso).'
              goto 900
            end if
c
c               use the row of the jacobian() matrix that matches the
c               constant isoparametric dof for the face normal vector;
c               the other rows from the jacobian() matrix give the
c               master element face tangent vectors
c
            normal(1:maxdof) = jacobian(constantisodof,1:maxdof)

            select case(constantisodof)
            case(1)
              tangent1(1:maxdof) = jacobian(2,1:maxdof)
              tangent2(1:maxdof) = jacobian(3,1:maxdof)
            case(2)
              tangent1(1:maxdof) = jacobian(1,1:maxdof)
              tangent2(1:maxdof) = jacobian(3,1:maxdof)
            case(3)
              tangent1(1:maxdof) = jacobian(1,1:maxdof)
              tangent2(1:maxdof) = jacobian(2,1:maxdof)
            case default
              mflag = 2
              message = 'internal error: unexpected constantisodof'//
     &                ' value; cannot set the element face tangent '//
     &                ' vectors (tied_nodeiso).'
              goto 900
            end select

          case(6,13)      ! 10 or 4 node tetrahedron
c
c               use the triangular face corners to get orthogonal
c               tangent vectors on the tetrahedron element face
c
            call tied_tetfacetangent(tangent1,tangent2,normal,maxdof,
     &                  elemcoord,maxnodeperelem,maxisodof,
     &                  node_per_elem,face_id,
     &                  localfacenode,maxelemface,maxnodeperface,
     &                  numelemface,numnodeperface,numcornerperface,
     &                  mflag,message,logflag,logunit,logfile)
            if(mflag.ge.2)goto 900

          case default
            mflag = 2
            message = 'error: unexpected elemchoice; cannot set the'//
     &              ' initial isoparametric point location'//
     &              ' (tied_nodeiso).'
            goto 900
          end select
c
c               normalize to get unit vectors
c
          call srt_normalize(normal,maxdof)
          call srt_normalize(tangent1,maxdof)
          call srt_normalize(tangent2,maxdof)

          if(local_debug .and. logflag.ge.debug_level)then
            open(unit=logunit,file=logfile,position='append')
            write(logunit,*)
            write(logunit,*)'     itercount =',itercount
            write(logunit,*)'          detj =',detj
            write(logunit,3)'     isopoint() =',isopoint(1:maxdof)
            write(logunit,3)'  globalpoint() =',globalpoint(1:maxdof)
            write(logunit,3)'     position() =',position(1:maxdof)
            write(logunit,3)'       normal() =',normal(1:maxdof)
            write(logunit,3)'     tangent1() =',tangent1(1:maxdof)
            write(logunit,3)'     tangent2() =',tangent2(1:maxdof)
            close(unit=logunit)
          end if
3         format(a,3es16.6)
c
c               check the distance from the face point to the slave
c               node; if close enough then the face point is located
c               at the slave node on the master element face
c
          dist_point = srt_dist3d(globalpoint,slavenode,maxdof)

          if(dist_point.le.epsilon)then
            if(local_debug .and. logflag.ge.debug_level)then
              open(unit=logunit,file=logfile,position='append')
              write(logunit,*)
              write(logunit,*)'face point location converged (distance)'
              write(logunit,*)'     itercount =',itercount
              write(logunit,*)'    dist_point =',dist_point
              write(logunit,*)'       epsilon =',epsilon
              write(logunit,*)'    globalpoint(1) =',globalpoint(1)
              write(logunit,*)'    globalpoint(2) =',globalpoint(2)
              write(logunit,*)'    globalpoint(3) =',globalpoint(3)
              write(logunit,*)'  slave node row index, nd =',nd
              write(logunit,*)'             slave node_id =',node_id
              write(logunit,*)'            master elem_id =',elem_id
              write(logunit,*)'    master element face_id =',face_id
              write(logunit,*)
              close(unit=logunit)
            end if

            exit loop_1
          end if
c
c               check the position change of the element face point
c               from the previous iteration; if the change of point
c               position is small enough then the face point is
c               converged to the mapped slave node location
c
          if(itercount.ge.2)then
            delta_iso = srt_dist3d(isopoint,isopoint_old,maxdof)

            if(delta_iso.le.epsilon)then
              if(local_debug .and. logflag.ge.debug_level)then
                open(unit=logunit,file=logfile,position='append')
                write(logunit,*)
                write(logunit,*)'face point location converged'//
     &                          ' (small iso change)'
                write(logunit,*)'     itercount =',itercount
                write(logunit,*)'     delta_iso =',delta_iso
                write(logunit,*)'       epsilon =',epsilon
                write(logunit,*)'    globalpoint(1) =',globalpoint(1)
                write(logunit,*)'    globalpoint(2) =',globalpoint(2)
                write(logunit,*)'    globalpoint(3) =',globalpoint(3)
                write(logunit,*)'  slave node row index, nd =',nd
                write(logunit,*)'             slave node_id =',node_id
                write(logunit,*)'            master elem_id =',elem_id
                write(logunit,*)'    master element face_id =',face_id
                write(logunit,*)
                close(unit=logunit)
              end if

              exit loop_1
            end if
          end if
c
c               check the angle from the face point to the slave node;
c               the slave node may be away from the master element face
c               but if the angle between the face normal and the
c               position is small then the slave node is projected onto
c               the master element face at the face point location
c
          angle_point = srt_getangle(position,normal,maxdof)

          if(angle_point.le.epsilon_angle .or.
     &     dabs(angle_point-pi).le.epsilon_angle)then
            if(local_debug .and. logflag.ge.debug_level)then
              open(unit=logunit,file=logfile,position='append')
              write(logunit,*)
              write(logunit,*)'face point location converged (angle)'
              write(logunit,*)'     itercount =',itercount
              write(logunit,*)'   angle_point =',angle_point
              write(logunit,*)' epsilon_angle =',epsilon_angle
              write(logunit,*)'    globalpoint(1) =',globalpoint(1)
              write(logunit,*)'    globalpoint(2) =',globalpoint(2)
              write(logunit,*)'    globalpoint(3) =',globalpoint(3)
              write(logunit,*)'  slave node row index, nd =',nd
              write(logunit,*)'             slave node_id =',node_id
              write(logunit,*)'            master elem_id =',elem_id
              write(logunit,*)'    master element face_id =',face_id
              write(logunit,*)
              close(unit=logunit)
            end if

            exit loop_1
          end if
c
c               set the scaling factors used to update the face point
c               location
c
          if(itercount.eq.1)then
            scale1 = zero
            scale2 = zero

            select case(elemchoice)
            case(1,2)      ! 20 or 8 node brick element
c
c               compute the average projected distance of the vectors
c               from the master element face center (face point in the
c               first iteration) to the corner nodes onto the tangent
c               vectors; use the scale1,scale2 values to scale the
c               delta1,delta2 values (this is an approximate mapping of
c               the element size back to the isoparametric element size)
c
              do corner=1,numcornerperface
                node_row = localfacenode(face_id,corner)
                call srt_getunitvec(globalpoint,
     &                  elemcoord(node_row,1:maxdof),vec,maxdof,0)
                scale1 = scale1 +
     &                   abs(srt_dotprod(vec,tangent1,maxdof))
                scale2 = scale2 +
     &                   abs(srt_dotprod(vec,tangent2,maxdof))
              end do
              scale1 = scale1/dble(numcornerperface)
              scale2 = scale2/dble(numcornerperface)

            case(6,13)      ! 10 or 4 node tetrahedron
              scale1 = one
              scale2 = one
            case default
              scale1 = one
              scale2 = one
            end select

            if(local_debug .and. logflag.ge.debug_level)then
              open(unit=logunit,file=logfile,position='append')
              write(logunit,*)'  numcornerperface =',numcornerperface
              write(logunit,*)'            scale1 =',scale1
              write(logunit,*)'            scale2 =',scale2
              close(unit=logunit)
            end if

            if(dabs(scale1).lt.1.0d-8)then
              mflag = 2
              message = 'internal error: scale1 = 0.0;'//
     &                ' cannot scale the delta1 value '//
     &                ' (tied_nodeiso).'
            else if(dabs(scale2).lt.1.0d-8)then
              mflag = 2
              message = 'internal error: scale2 = 0.0;'//
     &                ' cannot scale the delta2 value '//
     &                ' (tied_nodeiso).'
            end if
            if(mflag.ge.2)goto 900
          end if
c
c               project the position vector onto the face tangent
c               vectors; use the local element face coordinates to
c               update the isoparametric point position
c
          delta1 = srt_dotprod(position,tangent1,maxdof)/scale1
          delta2 = srt_dotprod(position,tangent2,maxdof)/scale2
          isopoint_old(1:maxdof) = isopoint(1:maxdof)

          if(local_debug .and. logflag.ge.debug_level)then
            open(unit=logunit,file=logfile,position='append')
            write(logunit,*)'face point not converged; update location'
            write(logunit,*)'     itercount =',itercount
            write(logunit,*)'    dist_point =',dist_point
            write(logunit,*)'     delta_iso =',delta_iso
            write(logunit,*)'   angle_point =',angle_point
            write(logunit,*)'       epsilon =',epsilon
            write(logunit,*)' epsilon_angle =',epsilon_angle
            write(logunit,3)'  globalpoint() =',globalpoint(1:maxdof)
            write(logunit,3)'    slavenode() =',slavenode(1:maxdof)
            write(logunit,*)'        delta1 =',delta1
            write(logunit,*)'        delta2 =',delta2
            close(unit=logunit)
          end if

          select case(elemchoice)
          case(1,2)      ! 20 or 8 node brick element
            select case(constantisodof)
            case(1)
              isopoint(2) = isopoint(2) + delta1
              isopoint(3) = isopoint(3) + delta2
            case(2)
              isopoint(1) = isopoint(1) + delta1
              isopoint(3) = isopoint(3) + delta2
            case(3)
              isopoint(1) = isopoint(1) + delta1
              isopoint(2) = isopoint(2) + delta2
            case default
              mflag = 2
              message = 'internal error: unexpected constantisodof'//
     &                ' value; cannot update the element face point'//
     &                ' location (tied_nodeiso).'
              goto 900
            end select

          case(6,13)      ! 10 or 4 node tetrahedron
            do i=1,maxdof
              vec1(i) = delta1*tangent1(i)
              vec2(i) = delta2*tangent2(i)
            end do
            do i=1,maxdof
              globalpoint(i) = globalpoint(i) + vec1(i) + vec2(i)
            end do

            if(local_debug .and. logflag.ge.debug_level)then
              open(unit=logunit,file=logfile,position='append')
              write(logunit,*)'update location for tetrahedron'
              write(logunit,*)'        delta1 =',delta1
              write(logunit,*)'        delta2 =',delta2
              write(logunit,3)'         vec1() =',vec1(1:maxdof)
              write(logunit,3)'         vec2() =',vec2(1:maxdof)
              write(logunit,3)'  globalpoint() =',globalpoint(1:maxdof)
              close(unit=logunit)
            end if

            call tied_tetfaceisopoint(isotetpoint,maxtetdof,
     &                  globalpoint,maxdof,constantisodof,
     &                  elemcoord,maxnodeperelem,maxisodof,
     &                  node_per_elem,face_id,
     &                  localfacenode,maxelemface,maxnodeperface,
     &                  numelemface,numnodeperface,numcornerperface,
     &                  mflag,message,logflag,logunit,logfile)
            if(mflag.ge.2)goto 900

            do i=1,maxtetdof
              if(isotetpoint(i).lt.-0.01d0)then
                mflag = 2
                message = 'internal error: tetrahedron isoparametric'//
     &                ' coordinate < 0.0; check the log file'//
     &                ' (tied_nodeiso).'
              else if(isotetpoint(i).gt.1.01d0)then
                mflag = 2
                message = 'internal error: tetrahedron isoparametric'//
     &                ' coordinate > 1.01; check the log file'//
     &                ' (tied_nodeiso).'
              end if
              if(mflag.ge.2)then
                if(logflag.ge.1)then
                  open(unit=logunit,file=logfile,position='append')
                  write(logunit,*)
                  write(logunit,*)trim(message)
                  write(logunit,*)'            mflag =',mflag
                  write(logunit,*)'   isotetpoint(1) =',isotetpoint(1)
                  write(logunit,*)'   isotetpoint(2) =',isotetpoint(2)
                  write(logunit,*)'   isotetpoint(3) =',isotetpoint(3)
                  write(logunit,*)'   isotetpoint(4) =',isotetpoint(4)
                  write(logunit,*)'        itercount =',itercount
                  write(logunit,*)'          node_id =',node_id
                  write(logunit,*)'          elem_id =',elem_id
                  write(logunit,*)'          face_id =',face_id
                  write(logunit,*)'           delta1 =',delta1
                  write(logunit,*)'           delta2 =',delta2
                  write(logunit,3)'     globalpoint() =',
     &                                  globalpoint(1:maxdof)
                  write(logunit,3)'       slavenode() =',
     &                                    slavenode(1:maxdof)
                  write(logunit,*)
                  close(unit=logunit)
                end if
                goto 900
              end if
            end do
c
c               save s2,s3,s4 as the independent isoparametric points on
c               the tetrahedron face, use isopoint() when evaluating the
c               shape functions
c
            do i=1,maxdof
              isopoint(i) = isotetpoint(i+1)
            end do

          case default
            mflag = 2
            message = 'error: unexpected elemchoice; cannot update'//
     &              ' the isoparametric point location'//
     &              ' (tied_nodeiso).'
            goto 900
          end select

          if(local_debug .and. logflag.ge.debug_level)then
            open(unit=logunit,file=logfile,position='append')
            write(logunit,*)'        delta1 =',delta1
            write(logunit,*)'        delta2 =',delta2
            write(logunit,*)'updated face point position:'
            write(logunit,*)'   isopoint(1) =',isopoint(1)
            write(logunit,*)'   isopoint(2) =',isopoint(2)
            write(logunit,*)'   isopoint(3) =',isopoint(3)
            close(unit=logunit)
          end if
        end do loop_1

        if(itercount.lt.min_iter_count)then
          min_iter_count = itercount
          min_iter_nodeid = node_id
        end if
        if(itercount.gt.max_iter_count)then
          max_iter_count = itercount
          max_iter_nodeid = node_id
        end if

        if(itercount.ge.maxiterloop)then
          mflag = 1
          message = 'warning: maximum iterations used to locate the'//
     &              ' slave node on the master element face; check'//
     &              ' the isoparametric element face coordinates'//
     &              ' (tied_nodeiso).'
          if(logflag.ge.1)then
            call srt_writelogmessage(mflag,message,logfile,logunit)
          end if
        end if
c
c               check if the face point isoparametric coordinates are
c               outside the +/- 1.0 range; expecting xi, eta, zeta to
c               be within +/- iso_error_val to allow for some round off;
c               can also reset the isoparametric value to be within the
c               element face, set the xi,eta,zeta value to be +/- 1.0
c
        iso_adjust_flag = .false.

        if(check_iso_coord_flag.gt.0)then
          if(abs(isopoint(1)).gt.iso_error_val)then
            xi_error_count = xi_error_count + 1
            select case(check_iso_coord_flag)
            case(1)      ! give a warning
              mflag = 1
              message = 'warning: the "xi" isoparametric coordinate'//
     &                ' is outside the element; expecting'//
     &                ' -1.0 <= xi <= 1.0; check the log file'//
     &                ' (tied_nodeiso).'
              if(logflag.ge.1)then
                open(unit=logunit,file=logfile,position='append')
                write(logunit,*)
                write(logunit,*)'         mflag =',mflag
                write(logunit,'(a)')trim(message)
                write(logunit,*)'   xi = isopoint(1) =',isopoint(1)
                write(logunit,*)'      iso_error_val =',iso_error_val
                write(logunit,*)'      slave node_id =',node_id
                write(logunit,*)'     xi_error_count =',xi_error_count
                write(logunit,*)
                close(unit=logunit)
              end if
            case(2)      ! set the maximum isoparametric value
              if(logflag.ge.1)then
                open(unit=logunit,file=logfile,position='append')
                write(logunit,*)
                write(logunit,'(a)')'note: set xi = limit value'
                write(logunit,*)'   xi = isopoint(1) =',isopoint(1)
                write(logunit,*)'   xi = reset to    =',
     &                          sign(one,isopoint(1))
                write(logunit,*)'      iso_error_val =',iso_error_val
                write(logunit,*)'      slave node_id =',node_id
                write(logunit,*)'     xi_error_count =',xi_error_count
                write(logunit,*)
                close(unit=logunit)
              end if
              isopoint(1) = sign(one,isopoint(1))
              iso_adjust_flag = .true.
            end select
          end if

          if(abs(isopoint(2)).gt.iso_error_val)then
            eta_error_count = eta_error_count + 1
            select case(check_iso_coord_flag)
            case(1)      ! give a warning
              mflag = 1
              message = 'warning: the "eta" isoparametric coordinate'//
     &                ' is outside the element; expecting'//
     &                ' -1.0 <= eta <= 1.0; check the log file'//
     &                ' (tied_nodeiso).'
              if(logflag.ge.1)then
                open(unit=logunit,file=logfile,position='append')
                write(logunit,*)
                write(logunit,*)'         mflag =',mflag
                write(logunit,'(a)')trim(message)
                write(logunit,*)'  eta = isopoint(2) =',isopoint(2)
                write(logunit,*)'      iso_error_val =',iso_error_val
                write(logunit,*)'      slave node_id =',node_id
                write(logunit,*)'    eta_error_count =',eta_error_count
                write(logunit,*)
                close(unit=logunit)
              end if
            case(2)      ! set the maximum isoparametric value
              if(logflag.ge.1)then
                open(unit=logunit,file=logfile,position='append')
                write(logunit,*)
                write(logunit,'(a)')'note: set eta = limit value'
                write(logunit,*)'  eta = isopoint(2) =',isopoint(2)
                write(logunit,*)'  eta = reset to    =',
     &                          sign(one,isopoint(2))
                write(logunit,*)'      iso_error_val =',iso_error_val
                write(logunit,*)'      slave node_id =',node_id
                write(logunit,*)'    eta_error_count =',eta_error_count
                write(logunit,*)
                close(unit=logunit)
              end if
              isopoint(2) = sign(one,isopoint(2))
              iso_adjust_flag = .true.
            end select
          end if

          if(abs(isopoint(3)).gt.iso_error_val)then
            zeta_error_count = zeta_error_count + 1
            select case(check_iso_coord_flag)
            case(1)      ! give a warning
              mflag = 1
              message = 'warning: the "zeta" isoparametric coordinate'//
     &                ' is outside the element; expecting'//
     &                ' -1.0 <= zeta <= 1.0; check the log file'//
     &                ' (tied_nodeiso).'
              if(logflag.ge.1)then
                open(unit=logunit,file=logfile,position='append')
                write(logunit,*)
                write(logunit,*)'         mflag =',mflag
                write(logunit,'(a)')trim(message)
                write(logunit,*)' zeta = isopoint(3) =',isopoint(3)
                write(logunit,*)'      iso_error_val =',iso_error_val
                write(logunit,*)'      slave node_id =',node_id
                write(logunit,*)'   zeta_error_count =',zeta_error_count
                write(logunit,*)
                close(unit=logunit)
              end if
            case(2)      ! set the maximum isoparametric value
              if(logflag.ge.1)then
                open(unit=logunit,file=logfile,position='append')
                write(logunit,*)
                write(logunit,'(a)')'note: set zeta = limit value'
                write(logunit,*)' zeta = isopoint(3) =',isopoint(3)
                write(logunit,*)' zeta = reset to    =',
     &                          sign(one,isopoint(3))
                write(logunit,*)'      iso_error_val =',iso_error_val
                write(logunit,*)'      slave node_id =',node_id
                write(logunit,*)'   zeta_error_count =',zeta_error_count
                write(logunit,*)
                close(unit=logunit)
              end if
              isopoint(3) = sign(one,isopoint(3))
              iso_adjust_flag = .true.
            end select
          end if
        end if
c
c               save the isoparametric coordinates of the slave node
c               on the master element
c
        tiednodeisocoord(nd,1:maxisodof) = isopoint(1:maxdof)
c
c               check if the node's global x,y,z location needs to be
c               updated due to an adjustment of the master element face
c               point getting adjusted; reevaluate the shape functions
c               to get the node's updated global coordinates
c
        if((iso_adjust_flag .and. tiednodeadjustflag(nd)) .or.
     &     tiednodegaptolerance(nd).gt.zero)then

          call shapef(elemchoice,isopoint(1),isopoint(2),isopoint(3),
     &                shapefunc)

          call tied_globalcoordfromiso(globalpoint,maxdof,
     &                shapefunc,maxnodeperelem,node_per_elem,
     &                elemcoord,maxnodeperelem,maxisodof,
     &                logflag,logunit,logfile,mflag,message)
          if(mflag.ge.2)goto 900
        end if
c
c               check if the distance between the element face point and
c               the slave node is within the given gap tolerance, set the
c               tied node flag = 1, otherwise the node is too far from
c               the master surface, set the tied node flag = 0
c
        if(tiednodegaptolerance(nd).gt.zero)then
          dist_gap = srt_dist3d(globalpoint,slavenode,maxdof)
          if(dist_gap.le.tiednodegaptolerance(nd))then
            tiednodeid(nd,4) = 1
          else
            tiednodeid(nd,4) = 0
            mflag = 1
            message = 'warning: a dependent slave node was not'//
     &                ' within the given gap tolerance of its'//
     &                ' assigned master element; the node will'//
     &                ' not be tied to the master surface;'//
     &                ' check the surface definitions (tied_nodeiso).'
c
c               report the slave node warning to the default file
c               (3/03/03 gvt)
c
            if(logunit.ne.out)then
              if (gap_cnt .eq. 0) then
                 allocate( gap_nodes(maxtiednode),
     &                     gap_elems(maxtiednode),
     &                     gap_faces(maxtiednode),
     &                     gap_dists(maxtiednode),
     &                     gap_toler(maxtiednode), stat=err)
                 if (err .ne. 0) then
                    call errmsg2(48,dumi,dums,dumr,dumd)
                    call die_abort
                 end if
              end if
              gap_cnt = gap_cnt + 1
              gap_nodes(gap_cnt) = node_id
              gap_elems(gap_cnt) = elem_id
              gap_faces(gap_cnt) = face_id
              gap_dists(gap_cnt) = dist_gap
              gap_toler(gap_cnt) = tiednodegaptolerance(nd)
            end if
c
c               report the slave node warning to the debugging log file
c
            if(logflag.ge.1)then
              open(unit=logunit,file=logfile,position='append')
              write(logunit,*)
              write(logunit,'(a)')trim(message)
              write(logunit,*)'  message flag, mflag         =',mflag
              write(logunit,*)'  dependent slave node id     =',node_id
              write(logunit,*)'  assigned master element id  =',elem_id
              write(logunit,*)'  master element face id      =',face_id
              write(logunit,*)'  node gap distance           =',dist_gap
              write(logunit,*)'  given gap tolerance         =',
     &                           tiednodegaptolerance(nd)
              write(logunit,*)'  slave node row, nd          =',nd
              write(logunit,*)'  tied node flag (set = 0)    =',
     &                           tiednodeid(nd,4)
              write(logunit,*)'  numtiednode                 =',
     &                           numtiednode
              write(logunit,*)
              close(unit=logunit)
            end if
          end if
        end if

        if(tiednodeadjustflag(nd))then
          tiednodeglobalcoord(nd,1:maxisodof) = globalpoint(1:maxdof)
          nodecoord(node_id,1:maxdof) = globalpoint(1:maxdof)
        end if
      end do

      if (gap_cnt .gt. 0) then
         write(out,9002)
         do cnt = 1, gap_cnt
            write(out,9003) gap_nodes(cnt),gap_elems(cnt),
     &                      gap_faces(cnt),gap_dists(cnt),
     &                      gap_toler(cnt)
         end do
         deallocate(gap_nodes,gap_elems,gap_faces,gap_dists,gap_toler)
      end if
c
c           write slave node information to the debugging log file
c
      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'tied slave node isoparametric location'
        write(logunit,*)'      numtiednode =',numtiednode
        write(logunit,*)'   min_iter_count =',min_iter_count
        write(logunit,*)'  min_iter_nodeid =',min_iter_nodeid
        write(logunit,*)'   max_iter_count =',max_iter_count
        write(logunit,*)'  max_iter_nodeid =',max_iter_nodeid
        write(logunit,*)'      xi_error_count =',xi_error_count
        write(logunit,*)'     eta_error_count =',eta_error_count
        write(logunit,*)'    zeta_error_count =',zeta_error_count
        write(logunit,5)'row','node_id','elem_id','face_id','flag',
     &                   'xi_iso1','eta_iso2','zeta_iso3'
        do row=1,numtiednode
          write(logunit,4)row,tiednodeid(row,1:maxtieddatacol),
     &                     tiednodeisocoord(row,1:maxisodof)
        end do
        if(local_debug)then
          write(logunit,*)
          write(logunit,*)'tied slave node updated global location'
          write(logunit,*)'      numtiednode =',numtiednode
          write(logunit,1)'row','node_id','elem_id','face_id','flag',
     &                     'adjust','x_dof1','y_dof2','z_dof3',
     &                     'x_original','y_original','z_original',
     &                     'delta'
          do row=1,numtiednode
            node_id = tiednodeid(row,1)      ! slave node id
            dist_point = srt_dist3d(tiednodeglobalcoord(row,1:dimcoord),
     &                      nodecoord(node_id,1:dimcoord),dimcoord)

            write(logunit,2)row,tiednodeid(row,1:maxtieddatacol),
     &                       tiednodeadjustflag(row),
     &                       tiednodeglobalcoord(row,1:dimcoord),
     &                       nodecoord(node_id,1:dimcoord),dist_point
          end do
          write(logunit,*)
          write(logunit,*)'tied slave node updated global location'
          write(logunit,*)'femap netural file node format:'
          write(logunit,*)'(copy and paste to femap data block 403)'
          write(logunit,*)'(may need to offset the node id numbers)'
          do row=1,numtiednode
            write(logunit,7)row,tiednodeglobalcoord(row,1:dimcoord)
          end do
          write(logunit,*)
        end if
        close(unit=logunit)
      end if
1     format(4a10,a5,a8,6a16,a16)
2     format(4i10,i5,l8,3es16.6,3es16.6,es16.6)
4     format(4i10,i5,3f16.8)
5     format(4a10,a5,3a16)
c          femap netural file node format:
7     format(i7,'   0 0 1 46 0 0 0 0 0 0  ',3es16.8,' 0')
9002  format(//,10x,'>> warning: tied contact processing',/,
     &      10x,'>>   the following dependent slave nodes were not',/,
     &      10x,'>>   within the given gap tolerance of their',/,
     &      10x,'>>   assigned master elements',/,
     &      10x,'>>   the nodes will not be tied to their master',/,
     &      10x,'>>   surfaces; check the surface definitions',//,
     & 5x,'slave node',3x,'master elem',3x,'face',6x,'gap distance',
     & 7x,'gap tolerance',/,
     & 80('-'))
9003  format(7x,i7,8x,i5,8x,i1,4x,es15.6,5x,es15.6)
c
c----------------------------------------------------------------------
c           jump here on an error for exit procedures
c
900   continue
c
c           deallocate local arrays
c
      if(allocated(shapefunc))deallocate(shapefunc)
      if(allocated(shapederivs))deallocate(shapederivs)
      if(allocated(elemcoord))deallocate(elemcoord)

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'----- end tied_nodeiso -----'
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
          tied_con_mpc_table(num_tied_con_mpc)%constant = zero
          tied_con_mpc_table(num_tied_con_mpc)%node_list(numlist) =
     &      node_id
          tied_con_mpc_table(num_tied_con_mpc)%dof_list(numlist) =
     &      local_dof
          tied_con_mpc_table(num_tied_con_mpc)%multiplier_list(numlist)=
     &      -one
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
     &          multiplier_list(numlist) = shapefunc(column)

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
     &        norm_mult
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
        masternodelist = zero
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
c
c     ****************************************************************
c     *                                                              *
c     *                subroutine srt_checkduplicate                 *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 06/10/04 gvt              *
c     *                                                              *
c     *  check for duplicate integer values in the rows of the       *
c     *  nodelist() array; for example, node id numbers could be     *
c     *  listed in the nodelist() array; remove any duplicated node  *
c     *  numbers and reduce the count of numlist                     *
c     *                                                              *
c     ****************************************************************
c
c         variables:
c
c         nodelist() = the list of integer data values;
c           for example, a list of node id numbers
c         dimlist = row dimension for the nodelist() array
c         numlist = number of rows of values (nodes or other integer
c           values) listed in the nodelist() array;
c           require numlist <= dimlist
c
c         mflag = 0=no error, 1=warning message, 2=error message
c         message = text string with a warning or error message when
c           mflag >= 1
c         logflag = for debugging output, 0=no output, >=1=output write
c           to the log file
c         logunit = debugging log file unit number
c         logfile = debugging log file name
c
      subroutine srt_checkduplicate(nodelist,dimlist,numlist,
     &                        logflag,logunit,logfile,mflag,message)
c
c         declare variables
c
      implicit none
      integer
     &  dimlist,numlist,mflag,logflag,logunit,
     &  nodelist(dimlist)
      character(len=*) :: message,logfile
c
c         local variables
c
      integer
     &  nd1,nd2,i,k,count,numold,allocate_error
      integer, allocatable :: newid(:)
      integer, allocatable :: mergeflag(:)
      integer, allocatable :: renumberlist(:)
      integer, allocatable :: updateid(:)
      logical
     &  done
c
c         start the log file
c
      if(logflag.ge.8)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'---- begin srt_checkduplicate ----'
        write(logunit,*)'       numlist =',numlist
        write(logunit,*)'       dimlist =',dimlist
        if(logflag.ge.10)then
          write(logunit,*)'array before checking for duplicates:'
          write(logunit,*)'    i,  nodelist'
          do i=1,numlist
            write(logunit,*)i,nodelist(i)
          end do
        end if
        close(unit=logunit)
      end if
c
c         error checking
c
      if(numlist.gt.dimlist)then
        mflag = 2
        message = 'error: too many values in the list;'//
     &    ' numlist > dimlist; check values of "dimlist" and'//
     &    ' "numlist" (srt_checkduplicate).'
      end if
      if(mflag.ge.2)goto 900
c
c         initialize values
c
      allocate(newid(numlist),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "newid()" array'//
     &    ' (srt_checkduplicate).'
        goto 900
      end if

      allocate(mergeflag(numlist),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "mergeflag()" array'//
     &    ' (srt_checkduplicate).'
        goto 900
      end if

      allocate(renumberlist(numlist),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "renumberlist()" array'//
     &    ' (srt_checkduplicate).'
        goto 900
      end if

      allocate(updateid(numlist),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "updateid()" array'//
     &    ' (srt_checkduplicate).'
        goto 900
      end if

      do i=1,numlist
        newid(i) = i
      end do

      mergeflag = 0
      renumberlist = 0
      updateid = 0

      numold = numlist
      count = 0
c
c         loop through the integer value list and check for
c         duplicated integer values (duplicate node numbers)
c
      do nd1=1,numlist
        done = .false.
        nd2 = nd1

        do while(.not.done)
          nd2 = nd2 + 1

          if(nd2.gt.numlist)then
            done=.true.
          else if((mergeflag(nd1).eq.0) .and.
     &    (nodelist(nd1).eq.nodelist(nd2)))then
c
c                                    flag the duplicate integer value in the lis
c
            count = count + 1
            mergeflag(nd2) = 1      ! merge node j
            newid(nd2) = nd1            ! replace node id
c
c                                    found a duplicate node id
c
            if(.false. .and. logflag.ge.10)then
              open(unit=logunit,file=logfile,position='append')
              write(logunit,*)'duplicate node found: nd1 =',nd1,
     &            '  nd2 =',nd2
              write(logunit,*)'            count =',count
              write(logunit,*)'    nodelist(nd1) =',nodelist(nd1),
     &            '  nodelist(nd2) =',nodelist(nd2)
              close(unit=logunit)
            end if

          end if
        end do
      end do
c
c         build two more integer pointer arrays.
c
c         renumberlist() contains the sequential node id numbering for
c           the remaining nodes.
c         updateid() contains the new node id numbering for all the
c           remaining and old/merged nodes; the nodes being merged need
c           to have their node id set to the renumbered node id of the
c           coincident node.
c
      k = 0
      do i=1,numlist
        if(mergeflag(i).eq.0)then
          k = k + 1
          renumberlist(i) = k
          updateid(i) = k
        else
          renumberlist(i) = 0
          updateid(i) = renumberlist(newid(i))
        end if
      end do
c
      if(logflag.ge.10)then      ! log file output
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'after comparing the node list:'
        write(logunit,*)'       count =',count
        write(logunit,*)'      numold =',numold
        write(logunit,*)'integer pointer arrays:'
        write(logunit,'(5a10)')'i','mergeflag','newid','renumlist',
     &                         'updateid'
        do i=1,numlist
          write(logunit,'(5i10)')i,mergeflag(i),newid(i),
     &                           renumberlist(i),updateid(i)
        end do
        close(unit=logunit)
      end if
c
c         update the arrays to remove rows with duplicate values;
c         use the integer pointer arrays to update the arrays;
c         to remove the merged nodes shift the nodes lower in the array
c         up to the empty rows.
c
      do i=1,numlist
        if(mergeflag(i).eq.0)then            ! skip merged rows
          nd1 = newid(i)                        ! old node id of remaining node
          nd2 = updateid(i)                        ! shifted/new node id
          if(nd1.ne.nd2)then
            nodelist(nd2) = nodelist(nd1)
          end if
        end if
      end do
c
c         update the number of remaining nodes.
c
      numlist = numlist - count
c
c         erase values in the last "count" rows of the arrays
c
      do i=(numlist+1),numold
        nodelist(i) = 0
      end do
c
c         jump here on an error
c
900   continue
c
c         deallocate local arrays
c
      if(allocated(newid))deallocate(newid)
      if(allocated(mergeflag))deallocate(mergeflag)
      if(allocated(renumberlist))deallocate(renumberlist)
      if(allocated(updateid))deallocate(updateid)

      if(logflag.ge.8)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'after duplicate check:'
        write(logunit,*)'         count =',count
        write(logunit,*)'       numlist =',numlist
        if(logflag.ge.10)then
          write(logunit,*)'array after removing duplicates:'
          write(logunit,*)'    i,  nodelist'
          do i=1,numlist
            write(logunit,*)i,nodelist(i)
          end do
        end if
        write(logunit,*)'---- end srt_checkduplicate ----'
        close(unit=logunit)
      end if
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *           subroutine srt_checkduplicatemultiarray            *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 12/30/02                  *
c     *                                    02/26/03 gvt              *
c     *                                                              *
c     *  check for duplicate integer values in the rows of the       *
c     *  nodelist() array; for example, node id numbers could be     *
c     *  listed in the nodelist() array; remove any duplicated node  *
c     *  numbers and reduce the count of numlistnode;                *
c     *  keep the values on the same row of the various arrays       *
c     *  together on the same row during the update.                 *
c     *                                                              *
c     ****************************************************************
c
c         variables:
c
c         nodelist() = the primary list of integer data values,
c           for example, a list of node id numbers; keep the values on
c           the same row of the various arrays together on the same row
c           during the update
c         nodelogical() = a matching list of logical values for each
c           item in the list
c         nodereal() = a matching list of real values for each item in
c           the list
c         dimlist = row dimension for the nodelist() array
c         numlistnode = number of nodes (integer values) listed in the
c           nodelist() array, need numlistnode <= dim,
c
c         mflag = 0=no error, 1=warning message, 2=error message
c         message = text string with a warning or error message when
c           mflag >= 1
c         logflag = for debugging output, 0=no output, >=1=output write
c           to the log file
c         logunit = debugging log file unit number
c         logfile = debugging log file name
c
      subroutine srt_checkduplicatemultiarray(nodelist,nodelogical,
     &                        nodereal,dimlist,numlistnode,
     &                        mflag,message,logflag,logunit,logfile)
c
c         declare variables
c
      implicit none
      integer
     &  dimlist,numlistnode,mflag,logflag,logunit,
     &  nodelist(dimlist)
      logical
     &  nodelogical(dimlist)
      real
     &  nodereal(dimlist)
      character(len=*) :: message,logfile
c
c         local variables
c
      integer
     &  nd1,nd2,i,k,count,numold,allocate_error
      integer, allocatable :: newid(:)
      integer, allocatable :: mergeflag(:)
      integer, allocatable :: renumberlist(:)
      integer, allocatable :: updateid(:)
      logical
     &  done
c
c         start the log file
c
      if(logflag.ge.8)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'---- begin srt_checkduplicatemultiarray ----'
        write(logunit,*)'   numlistnode =',numlistnode
        write(logunit,*)'       dimlist =',dimlist
        if(logflag.ge.10)then
          write(logunit,*)'arrays before checking for duplicates:'
          write(logunit,*)'    i,  nodelist, nodelogical, nodereal'
          do i=1,numlistnode
            write(logunit,*)i,nodelist(i),nodelogical(i),nodereal(i)
          end do
        end if
        close(unit=logunit)
      end if
c
c         error checking
c
      if(numlistnode.gt.dimlist)then
        mflag = 2
        message = 'error: too many values in the list;'//
     &    ' numlistnode > dimlist; check values of "dimlist" and'//
     &    ' "numlistnode" (srt_checkduplicatemultiarray).'
      end if
      if(mflag.ge.2)goto 900
c
c         initialize values
c
      allocate(newid(numlistnode),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "newid()" array'//
     &    ' (srt_checkduplicatemultiarray).'
        goto 900
      end if

      allocate(mergeflag(numlistnode),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "mergeflag()" array'//
     &    ' (srt_checkduplicatemultiarray).'
        goto 900
      end if

      allocate(renumberlist(numlistnode),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "renumberlist()" array'//
     &    ' (srt_checkduplicatemultiarray).'
        goto 900
      end if

      allocate(updateid(numlistnode),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "updateid()" array'//
     &    ' (srt_checkduplicatemultiarray).'
        goto 900
      end if

      do i=1,numlistnode
        newid(i) = i
      end do

      mergeflag = 0
      renumberlist = 0
      updateid = 0

      numold = numlistnode
      count = 0
c
c         loop through the node list and check for duplicated
c         node numbers
c
      do nd1=1,numlistnode
        done = .false.
        nd2 = nd1

        do while(.not.done)
          nd2 = nd2 + 1

          if(nd2.gt.numlistnode)then
            done=.true.
          else if((mergeflag(nd1).eq.0) .and.
     &    (nodelist(nd1).eq.nodelist(nd2)))then
c
c                                    flag the duplicate integer value in the lis
c
            count = count + 1
            mergeflag(nd2) = 1      ! merge node j
            newid(nd2) = nd1            ! replace node id
c
c                                    found a duplicate node id
c
            if(.false. .and. logflag.ge.10)then
              open(unit=logunit,file=logfile,position='append')
              write(logunit,*)'duplicate node found: nd1 =',nd1,
     &            '  nd2 =',nd2
              write(logunit,*)'            count =',count
              write(logunit,*)'    nodelist(nd1) =',nodelist(nd1),
     &            '  nodelist(nd2) =',nodelist(nd2)
              close(unit=logunit)
            end if

          end if
        end do
      end do
c
c         build two more integer pointer arrays.
c
c         renumberlist() contains the sequential node id numbering for
c           the remaining nodes.
c         updateid() contains the new node id numbering for all the
c           remaining and old/merged nodes; the nodes being merged need
c           to have their node id set to the renumbered node id of the
c           coincident node.
c
      k = 0
      do i=1,numlistnode
        if(mergeflag(i).eq.0)then
          k = k + 1
          renumberlist(i) = k
          updateid(i) = k
        else
          renumberlist(i) = 0
          updateid(i) = renumberlist(newid(i))
        end if
      end do
c
      if(logflag.ge.10)then      ! log file output
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'after comparing the node list:'
        write(logunit,*)'       count =',count
        write(logunit,*)'      numold =',numold
        write(logunit,*)'integer pointer arrays:'
        write(logunit,'(5a10)')'i','mergeflag','newid','renumlist',
     &                         'updateid'
        do i=1,numlistnode
          write(logunit,'(5i10)')i,mergeflag(i),newid(i),
     &                           renumberlist(i),updateid(i)
        end do
        close(unit=logunit)
      end if
c
c         update the arrays to remove rows with duplicate values;
c         use the integer pointer arrays to update the arrays;
c         to remove the merged nodes shift the nodes lower in the array
c         up to the empty rows.
c
      do i=1,numlistnode
        if(mergeflag(i).eq.0)then            ! skip merged rows
          nd1 = newid(i)                        ! old node id of remaining node
          nd2 = updateid(i)                  ! shifted/new node id
          if(nd1.ne.nd2)then
            nodelist(nd2) = nodelist(nd1)
            nodelogical(nd2) = nodelogical(nd1)
            nodereal(nd2) = nodereal(nd1)
          end if
        end if
      end do
c
c         update the number of remaining nodes.
c
      numlistnode = numlistnode - count
c
c         erase values in the last "count" rows of the arrays
c
      do i=(numlistnode+1),numold
        nodelist(i) = 0
        nodelogical(i) = .false.
        nodereal(i) = 0.0
      end do
c
c         jump here on an error
c
900   continue
c
c         deallocate local arrays
c
      if(allocated(newid))deallocate(newid)
      if(allocated(mergeflag))deallocate(mergeflag)
      if(allocated(renumberlist))deallocate(renumberlist)
      if(allocated(updateid))deallocate(updateid)

      if(logflag.ge.8)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'after duplicate check:'
        write(logunit,*)'         count =',count
        write(logunit,*)'   numlistnode =',numlistnode
        if(logflag.ge.10)then
          write(logunit,*)'arrays after removing duplicates:'
          write(logunit,*)'    i,  nodelist, nodelogical, nodereal'
          do i=1,numlistnode
            write(logunit,*)i,nodelist(i),nodelogical(i),nodereal(i)
          end do
        end if
        write(logunit,*)'---- end srt_checkduplicatemultiarray ----'
        close(unit=logunit)
      end if
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *  compute the vector cross product:  {vec} = {a} x {b}        *
c     *                                                              *
c     ****************************************************************
c
c         variables:
c
c         a(),b() = 3x1 vectors
c         vec() = return the vector cross product
c         dim = vector dimensions, expecting dim=3 for x,y,z coordinates
c
      subroutine srt_crossprod(a,b,vec,dim)
c
c         declare variables
c
      implicit none
      integer dim
      double precision
     &  a(dim),b(dim),vec(dim)
c
c         compute the vector cross product
c
      vec(1) = a(2)*b(3) - a(3)*b(2)
      vec(2) = a(3)*b(1) - a(1)*b(3)
      vec(3) = a(1)*b(2) - a(2)*b(1)
c
      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *  compute the scalar dot product of two vectors:              *
c     *     projected length = {a}.{b}                               *
c     *                                                              *
c     ****************************************************************
c
      function srt_dotprod(a,b,dim)
c
c           declare variables
c

      implicit none
      integer dim
      double precision
     &  srt_dotprod,a(dim),b(dim)
c
c           local variables
c
      integer i
      double precision
     &  sum
c
      sum = 0.0d00

      if(dim.eq.3)then
        sum = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
      else
        do i=1,dim
          sum = sum + a(i)*b(i)
        end do
      end if

      srt_dotprod = sum
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *  compute the angle between two vectors.                      *
c     *  angle value returned is in radians.                         *
c     *  use the value of the dot product and an inverse cosine:     *
c     *        {a}.{b} = |a|*|b|*cos(theta)                          *
c     *        theta = acos({a}.{b} / |a|*|b|)                       *
c     *                                                              *
c     ****************************************************************
c
      function srt_getangle(a,b,dim)
c
c         declare variables
c
      implicit none
      integer dim
      double precision
     &  srt_getangle,a(dim),b(dim)
c
c         a(),b() = vectors
c         dim = vector dimension
c
c
c         local variables
c
      double precision
     &  dot,amag,bmag,x, one
           data one / 1.0d00 /
c
c         declare functions
c
      double precision
     &  srt_dotprod,srt_veclength
c
c         compute the dot product and vector magnitude/lengths
c
      dot  = srt_dotprod(a,b,dim)
      amag = srt_veclength(a,dim)
      bmag = srt_veclength(b,dim)
c
c         compute the angle between vectors a and b
c
      x = dot/(amag*bmag)

      if(abs(x).gt.one)then
        x = sign(one,x)
      end if
c
c         get the angle in the range 0 to pi, use the dacos() function
c
      srt_getangle = acos(x)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *  build a unit vector, nhat, along the line defined by 2      *
c     *  points: a, b.  the vector is oriented from point a to       *
c     *  point b                                                     *
c     *                                                              *
c     ****************************************************************
c
      subroutine srt_getunitvec(a,b,nhat,dim,normflag)
c
c         declare varibles
c
      implicit none
      integer
     &  dim,normflag
      double precision
     &  a(dim),b(dim),nhat(dim)
c
c         dim=vector dimensions, expecting dim=3 for x,y,z coordinates
c         a(),b() = coordinates of points on the line
c         nhat() = vector along the line, can be normalized for a unit vector
c         normflag = 0=don't normalize the vector from a to b
c            1=normalize to get a unit vector nhat
c
c
c         local varibles
c
      integer i
c
c         get components of direction vector from a to b
c
      do i =1, dim
        nhat(i) = b(i) - a(i)
      end do
c
c            normalize the vector to get a unit direction vector
c
      if( normflag .gt. 0 )then
        call srt_normalize( nhat, dim )
      end if
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine srt_insideelemface                *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 11/15/02                  *
c     *                                    02/26/03 gvt              *
c     *                                    03/04/03 gvt              *
c     *                                    05/05/05 gvt              *
c     *                                                              *
c     *  check that the point is within the element face boundaries; *
c     *  from each element face corner check that the point is       *
c     *  within the angle between the element faces edges; return a  *
c     *  logical flag if the point is inside the element face, and   *
c     *  return an integer flag that indicates the status of how     *
c     *  well the point passed the geometric checks                  *
c     *                                                              *
c     *  added checkstatus = 50 logic when the point is inside the   *
c     *  element face (03/04/03 gvt)                                 *
c     *                                                              *
c     *  added normal_dist to the call statement; return the normal  *
c     *  distance from the element face to the point (05/05/04 gvt)  *
c     *                                                              *
c     ****************************************************************
c
      subroutine srt_insideelemface(insideflag,checkstatus,elemid,
     &                     faceid,elemchoice,nodesperelem,
     &                     point,maxdof,
     &                     localfacenode,maxelemface,maxnodeperface,
     &                     numelemface,numnodeperface,numcornerperface,
     &                     elemcon,dimelem,maxnodeperelem,numelem,
     &                     nodecoord,dimnode,dimcoord,numnode,
     &                     given_epsilon,normal_dist,
     &                     logflag,logunit,logfile,mflag,message)
c
c           declare variables
c
c           insideflag = return true if the point is found to be within
c                        the element face boundaries, return false if
c                        the point is outside the element face
c           checkstatus = return 0 if the point is outside the element
c                        face, return 1 if the point is inside the
c                        element face but passed only secondary checks,
c                        the point may be near or just outside an edge,
c                        return 2 if the point passed all the primary
c                        checks and is within the element face
c           normal_dist = return the average normal distance from the
c                        element face to the pointc
      implicit none
      logical
     &  insideflag
      integer
     &  checkstatus,elemid,faceid,elemchoice,nodesperelem,maxdof,
     &  maxelemface,maxnodeperface,numelemface,numnodeperface,
     &  numcornerperface,dimelem,maxnodeperelem,numelem,
     &  dimnode,dimcoord,numnode,
     &  logflag,logunit,mflag,
     &  localfacenode(maxelemface,maxnodeperface),
     &  elemcon(dimelem,maxnodeperelem)
      double precision
     &  point(maxdof),nodecoord(dimnode,dimcoord),given_epsilon,
     &  normal_dist
      character(len=*) :: logfile,message
c
c           local variables
c
      integer, parameter :: dim3 = 3
      integer, parameter :: maxfacecorner = 4
      integer, parameter :: maxedgecol = 2
      integer
     &  debug_level,i,corner,nd,nd1,nd2,count_inside,count_best,
     &  cornerid(maxfacecorner),edgeid(maxfacecorner,maxedgecol)
      logical
     &  local_debug
      double precision
     &  edge_vec_1(dim3),edge_vec_2(dim3),normal(dim3),position(dim3),
     &  xvec(dim3),yvec(dim3),dist,epsilon,dx,dy,dz,theta_edge2,
     &  point_angle,angle_epsilon,pi,angle_pt_edge1,angle_pt_edge2,
     &  edge_rad_1,edge_rad_2,corner_dist,zero,four,one,tol,
     &  hundred,one_eighty,sum_norm_dist,dot_edge_1,dot_edge_2,
     &  len_edge_1,len_edge_2
       data zero, four, one, tol, hundred, one_eighty
     &   / 0.0d00, 4.0d00, 1.0d00, 1.0d-06, 100.0d00, 180.0d00 /
      character checklabel*80
c
c           declare functions
c
      double precision
     &  srt_dotprod,srt_dist3d,srt_getangle,srt_veclength
c
      debug_level = 9
      local_debug = .false.

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'----- begin srt_insideelemface -----'
        write(logunit,*)'            elemid =',elemid
        write(logunit,*)'            faceid =',faceid
        write(logunit,*)'        elemchoice =',elemchoice
        write(logunit,*)'      nodesperelem =',nodesperelem
        write(logunit,*)'            maxdof =',maxdof
        write(logunit,*)'       numelemface =',numelemface
        write(logunit,*)'    numnodeperface =',numnodeperface
        write(logunit,*)'  numcornerperface =',numcornerperface
        write(logunit,*)'           numelem =',numelem
        write(logunit,*)'           numnode =',numnode
        write(logunit,3)'            point() =',point(1:maxdof)
        write(logunit,3)'      given_epsilon =',given_epsilon
        if(local_debug)then
          write(logunit,*)'       local_debug = ',local_debug
          write(logunit,*)'              dim3 =',dim3
          write(logunit,*)'     maxfacecorner =',maxfacecorner
        end if
        close(unit=logunit)
      end if
3     format(a,3es16.6)
c
c           error checking
c
      if(elemid.le.0)then
        mflag = 2
        message = 'error: elemid = 0 (srt_insideelemface).'
      else if(elemid.gt.numelem)then
        mflag = 2
        message = 'error: elemid > numelem (srt_insideelemface).'
      else if(faceid.le.0)then
        mflag = 2
        message = 'faceid: elemid = 0 (srt_insideelemface).'
      else if(faceid.gt.numelemface)then
        mflag = 2
        message = 'faceid: elemid > numelemface (srt_insideelemface).'
      else if(nodesperelem.le.0)then
        mflag = 2
        message = 'error: nodesperelem = 0 (srt_insideelemface).'
      else if(maxdof.ne.dim3)then
        mflag = 2
        message = 'error: maxdof /= 3 (srt_insideelemface).'
      else if(dimcoord.ne.dim3)then
        mflag = 2
        message = 'error: dimcoord /= 3 (srt_insideelemface).'
      end if
      if(mflag.ge.2)goto 900
c
c           initialize values
c
      pi = four*atan(one)
      if(given_epsilon.gt.zero)then
        epsilon = given_epsilon/hundred
      else
        epsilon = tol
      end if
      angle_epsilon = one*(pi/one_eighty)
      count_inside = 0
      count_best = 0
      insideflag = .false.
      checkstatus = 0
      normal_dist = zero
      sum_norm_dist = zero

      cornerid = 0
      do i=1,numcornerperface
        cornerid(i) = localfacenode(faceid,i)
      end do

      edgeid = 0
      select case(numcornerperface)
      case(3)      ! 3 corners on a triangular face
c           neighbor corner node for edge 1
        edgeid(1,1) = localfacenode(faceid,2)
        edgeid(2,1) = localfacenode(faceid,3)
        edgeid(3,1) = localfacenode(faceid,1)
c           neighbor corner node for edge 2
        edgeid(1,2) = localfacenode(faceid,3)
        edgeid(2,2) = localfacenode(faceid,1)
        edgeid(3,2) = localfacenode(faceid,2)

        if(elemchoice.ne.6 .and. elemchoice.ne.13)then
          mflag = 2
          message = 'error: unexpected elemchoice for 3 corners'//
     &    ' per element face; check elemchoice and numcornerperface'//
     &    ' (srt_insideelemface).'
        end if
      case(4)      ! 4 corners on a quadrillateral face
c           neighbor corner node for edge 1
        edgeid(1,1) = localfacenode(faceid,2)
        edgeid(2,1) = localfacenode(faceid,3)
        edgeid(3,1) = localfacenode(faceid,4)
        edgeid(4,1) = localfacenode(faceid,1)
c           neighbor corner node for edge 2
        edgeid(1,2) = localfacenode(faceid,4)
        edgeid(2,2) = localfacenode(faceid,1)
        edgeid(3,2) = localfacenode(faceid,2)
        edgeid(4,2) = localfacenode(faceid,3)

        if(elemchoice.ne.1 .and. elemchoice.ne.2)then
          mflag = 2
          message = 'error: unexpected elemchoice for 4 corners'//
     &    ' per element face; check elemchoice and numcornerperface'//
     &    ' (srt_insideelemface).'
        end if
      case default
        mflag = 2
        message = 'error: unexpected numcornerperface; cannot set'//
     &    ' the corner points per edge (srt_insideelemface).'
      end select
      if(mflag.ge.2)goto 900

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'local values:'
        write(logunit,*)'           epsilon =',epsilon
        write(logunit,*)'     angle_epsilon =',angle_epsilon
        write(logunit,*)'    i,    cornerid()'
        do i=1,numcornerperface
          write(logunit,*)i,cornerid(i)
        end do
        write(logunit,*)'    i,    edgeid()'
        do i=1,numcornerperface
          write(logunit,*)i,edgeid(i,1:maxedgecol)
        end do
        close(unit=logunit)
      end if
c
c           at each element face corner get two edge vectors from the
c           corner node to the neighboring corner nodes; the cross
c           product gives the normal vector to the element face at that
c           corner; use edge 1 as a local x-axis and get the local
c           y-axis from the cross product of the face normal and edge 1
c           vectors; compute the angle between edge 1 (local x-axis)
c           and edge 2; compute the angle from edge 1 to the point and
c           compare the point angle to the edge 2 angle to determine if
c           the point is within the current edges; any corner check
c           that gives the point outside those element face edges means
c           that the point is outside the element face;
c           for the point to be within the element face it must pass the
c           check at each corner node
c
      loop_1: do corner=1,numcornerperface
        nd = elemcon(elemid,cornerid(corner))
        nd1 = elemcon(elemid,edgeid(corner,1))
        nd2 = elemcon(elemid,edgeid(corner,2))

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)
          write(logunit,*)'  loop index, corner =',corner
          write(logunit,*)'    corner node, nd  =',nd
          write(logunit,*)'    edge 1 node, nd1 =',nd1
          write(logunit,*)'    edge 2 node, nd2 =',nd2
          close(unit=logunit)
        end if

        if(nd.le.0 .or. nd.gt.numnode)then
          mflag = 2
          message = 'error: invalid corner node; cannot'//
     &      ' get the element face edge vectors (srt_insideelemface).'
        else if(nd1.le.0 .or. nd1.gt.numnode)then
          mflag = 2
          message = 'error: invalid corner node for edge 1; cannot'//
     &      ' get the element face edge 1 vector (srt_insideelemface).'
        else if(nd2.le.0 .or. nd2.gt.numnode)then
          mflag = 2
          message = 'error: invalid corner node for edge 2; cannot'//
     &      ' get the element face edge 2 vector (srt_insideelemface).'
        end if
        if(mflag.ge.2)goto 900
c
c           check the distance from the corner node to the point;
c           if the point is coincident with the corner node then it is
c           within the element face; also avoid a zero length position
c           vector below
c
        dist = srt_dist3d(nodecoord(nd,:),point,dimcoord)

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)'              dist =',dist
          write(logunit,*)'           epsilon =',epsilon
          close(unit=logunit)
        end if

        if(dist.le.epsilon)then
          insideflag = .true.
          checkstatus = 100
          exit loop_1
        end if
c
c           get the edge, normal, and local x and y vectors
c
        call srt_getunitvec(nodecoord(nd,:),nodecoord(nd1,:),
     &                edge_vec_1,dimcoord,0)
        len_edge_1 = srt_veclength(edge_vec_1,dim3)

        call srt_getunitvec(nodecoord(nd,:),nodecoord(nd2,:),
     &                edge_vec_2,dimcoord,0)
        len_edge_2 = srt_veclength(edge_vec_2,dim3)

        call srt_crossprod(edge_vec_1,edge_vec_2,normal,dim3)

        xvec(1:dim3) = edge_vec_1(1:dim3)

        call srt_crossprod(normal,xvec,yvec,dim3)

        call srt_normalize(xvec,dim3)
        call srt_normalize(yvec,dim3)
        call srt_normalize(normal,dim3)

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,3)'       edge_vec_1() =',edge_vec_1(1:dim3)
          write(logunit,3)'         len_edge_1 =',len_edge_1
          write(logunit,3)'       edge_vec_2() =',edge_vec_2(1:dim3)
          write(logunit,3)'         len_edge_2 =',len_edge_2
          write(logunit,3)'           normal() =',normal(1:dim3)
          write(logunit,3)'             xvec() =',xvec(1:dim3)
          write(logunit,3)'             yvec() =',yvec(1:dim3)
          close(unit=logunit)
        end if
c
c           get the angle between edge 1 (the local x-axis) and edge 2
c
        dx = srt_dotprod(edge_vec_2,xvec,dim3)
        dy = srt_dotprod(edge_vec_2,yvec,dim3)
        theta_edge2 = atan2(dy,dx)
        if(theta_edge2.le.zero)then
          mflag = 2
          message = 'internal error: theta_edge2 <= 0.0; expecting'//
     &      ' the angle between the element face edges to be'//
     &      ' positive (srt_insideelemface).'
          goto 900
        end if

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)'angle from edge 1 to edge 2:'
          write(logunit,*)'                dx =',dx
          write(logunit,*)'                dy =',dy
          write(logunit,*)'       theta_edge2 =',theta_edge2
          close(unit=logunit)
        end if
c
c           get the angle between the local x-axis and the point;
c           use a position vector from the corner node to the point;
c           compare the point angle to the angle between edges 1 and 2;
c           compute the in plane distance from the corner
c
        call srt_getunitvec(nodecoord(nd,:),point,position,dimcoord,0)

        dx = srt_dotprod(position,xvec,dim3)
        dy = srt_dotprod(position,yvec,dim3)
        dz = srt_dotprod(position,normal,dim3)

        dot_edge_1 = srt_dotprod(position,edge_vec_1,dim3)
        if(len_edge_1.gt.zero)then
          dot_edge_1 = dot_edge_1/len_edge_1
        end if
        dot_edge_2 = srt_dotprod(position,edge_vec_2,dim3)
        if(len_edge_2.gt.zero)then
          dot_edge_2 = dot_edge_2/len_edge_2
        end if

        point_angle = atan2(dy,dx)
        corner_dist = sqrt(dx**2 + dy**2)
c
c           get the angle between each edge and the point; use the angle
c           and distance from the corner to the point to compute the
c           radial distance from the edge to the point;
c           compare the radial distance to the gap tolerance epsilon to
c           check if the point is close to either edge
c
        angle_pt_edge1 = srt_getangle(edge_vec_1,position,dim3)
        edge_rad_1 = dist*sin(angle_pt_edge1)
        angle_pt_edge2 = srt_getangle(edge_vec_2,position,dim3)
        edge_rad_2 = dist*sin(angle_pt_edge2)
c
c           add to the sum of the normal distance from the element face
c           to the point; get the average normal distance after the loop
c           below (5/5/2004 gvt)
c
        sum_norm_dist = sum_norm_dist + dz

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,*)'angle to the point:'
          write(logunit,3)'        position() =',position(1:dim3)
          write(logunit,*)'                dx =',dx
          write(logunit,*)'                dy =',dy
          write(logunit,*)'                dz =',dz
          write(logunit,*)'       point_angle =',point_angle
          write(logunit,*)'     angle_epsilon =',angle_epsilon
          write(logunit,*)'    angle_pt_edge1 =',angle_pt_edge1
          write(logunit,*)'        edge_rad_1 =',edge_rad_1
          write(logunit,*)'        dot_edge_1 =',dot_edge_1
          write(logunit,*)'        len_edge_1 =',len_edge_1
          write(logunit,*)'    angle_pt_edge2 =',angle_pt_edge2
          write(logunit,*)'        edge_rad_2 =',edge_rad_2
          write(logunit,*)'        dot_edge_2 =',dot_edge_2
          write(logunit,*)'        len_edge_2 =',len_edge_2
          write(logunit,*)'       corner_dist =',corner_dist
          close(unit=logunit)
        end if
c
c           compare the angles and distances to check if the point is
c           within the element face boundaries; the primary check is
c           for the point to be within the edge 1 to 2 angle for each
c           corner; secondary checks permit the point to be near or
c           just outside an element edge; also check the in plane
c           distance from the corner to the point
c
        checklabel = ' '
        if(point_angle.le.theta_edge2 .and. point_angle.ge.zero)then
          count_inside = count_inside + 1
          count_best = count_best + 1
          checklabel = 'point within angle range; point_angle check'
          checkstatus = checkstatus + 10
        else if(corner_dist.le.given_epsilon)then
          count_inside = count_inside + 1
          checklabel = 'point near corner; corner_dist check'
          checkstatus = checkstatus + 8
        else if(point_angle.le.(theta_edge2 + angle_epsilon) .and.
     &          point_angle.ge.-angle_epsilon)then
          count_inside = count_inside + 1
          checklabel = 'point near angle range; 2nd point_angle check'
          checkstatus = checkstatus + 5
        else if(edge_rad_1.le.given_epsilon .and.
     &          (dot_edge_1.ge.-given_epsilon .and.
     &           dot_edge_1.le.(given_epsilon+len_edge_1)))then
          count_inside = count_inside + 1
          checklabel = 'point close to edge 1; edge_rad_1 check'
          checkstatus = checkstatus + 1
        else if(edge_rad_2.le.given_epsilon .and.
     &          (dot_edge_2.ge.-given_epsilon .and.
     &           dot_edge_2.le.(given_epsilon+len_edge_2)))then
          count_inside = count_inside + 1
          checklabel = 'point close to edge 2; edge_rad_2 check'
          checkstatus = checkstatus + 1
        else
          checklabel = 'point outside the element face'
          insideflag = .false.
          checkstatus = 0
          exit loop_1
        end if

        if(local_debug .and. logflag.ge.debug_level)then
          open(unit=logunit,file=logfile,position='append')
          write(logunit,'(a)')trim(checklabel)
          write(logunit,*)'       checkstatus =',checkstatus
          write(logunit,*)'      count_inside =',count_inside
          write(logunit,*)'        count_best =',count_best
          close(unit=logunit)
        end if

      end do loop_1
c
c compute the average normal distance from the element face to the point
c
      normal_dist = sum_norm_dist/DBLE(numcornerperface)

      if(local_debug .and. logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'average normal distance to the point'
        write(logunit,*)'     sum_norm_dist =',sum_norm_dist
        write(logunit,*)'  numcornerperface =',numcornerperface
        write(logunit,*)'       normal_dist =',normal_dist
        close(unit=logunit)
      end if

      if(count_inside.ge.numcornerperface)then
        insideflag = .true.
c
c           set checkstatus = 50 if insideflag is true, that is when
c           the point is inside the element face
c           (03/04/03 gvt)
c           require that the count_best equals the number of face
c           corner nodes to set this high value for the checkstatus
c           (5/4/2004 gvt)
c
        if(insideflag .and. count_best.ge.numcornerperface)then
          checkstatus = 50
        end if
      end if
c
c           jump here on an error for exit procedures
c
900   continue

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'       normal_dist =',normal_dist
        write(logunit,*)'        insideflag = ',insideflag
        write(logunit,*)'       checkstatus = ',checkstatus
        write(logunit,*)'----- end srt_insideelemface -----'
        close(unit=logunit)
      end if

      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *               subroutine srt_setlocalfacenodes               *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 11/06/02                  *
c     *                                    02/26/03 gvt              *
c     *                                                              *
c     *  provide the list of local face node numbers on each face of *
c     *  the element; get the local element face node numbers which  *
c     *  are the columns in the connectivity array;                  *
c     *  the local element face nodes are given for an abaqus format *
c     *  element.                                                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine srt_setlocalfacenodes(localfacenode,maxelemface,
     &               maxnodeperface,numelemface,numnodeperface,
     &               numcornerperface,nodesperelem,elemchoice,
     &               logflag,logunit,logfile,mflag,message)
c
c          declare variables
c
      implicit none
      integer
     & maxelemface,maxnodeperface,
     & numelemface,numnodeperface,numcornerperface,nodesperelem,
     & elemchoice,
     & localfacenode(maxelemface,maxnodeperface),
     & logflag,logunit,mflag
      character(len=*) :: logfile,message
c
c          variables:
c
c          localfacenode() = integer array to store the local element
c            face node list; for 3-d brick elements there are 6 element
c            faces with 4 or 8 nodes on each element face; return the
c            local face node list on each element face
c          maxelemface = the maximum number of faces per element; row
c            dimension for the localfacenode() array
c          maxnodeperface = the maximum number of nodes per face;
c            column dimension for the localfacenode() array
c          numelemface = return number of element faces, 6 element
c            faces for 3-d brick elements, 4 element faces for 3-d
c            tetrahedron elements
c          numnodeperface = number of nodes on each element face, an
c            8-node brick element has 4 nodes per face, a 20-node
c            brick has 8 nodes per face, a 4-node tetrahedron has 3
c            nodes per face, a 10-node tetrahedron has 6 nodes per face
c          numcornerperface = return the number of corner nodes per
c            element face; the first numcornerperface nodes per row
c            are the face corner nodes
c          nodesperelem = the number of nodes per element; this value
c            must be consistent with the elemchoice value
c          elemchoice = flag to choose the element type;
c            1 = 20-node brick element, also set numelemface = 6 and
c                numnodeperface = 8; expect nodesperelem = 20
c            2 = 8-node brick element, also set numelemface = 6 and
c                numnodeperface = 4; expect nodesperelem = 8
c            6 = 10-node tetrahedron element, also set numelemface = 4
c                and numnodeperface = 6; expect nodesperelem = 10
c           13 = 4-node tetrahedron element, also set numelemface = 4
c                and numnodeperface = 3; expect nodesperelem = 4
c
c          local variables
c
      integer
     & debug_level,i,j,
     & facenode8(6,4),facenode20(6,8),tet_10node(4,6),tet_4node(4,3)
c
c          set the face nodes in outward normal right-hand order for the
c          8-node brick; match the order of the 4-node quadrillateral
c
      data (facenode8(1,j),j=1,4) / 4,3,2,1 /
      data (facenode8(2,j),j=1,4) / 5,6,7,8 /
      data (facenode8(3,j),j=1,4) / 1,2,6,5 /
      data (facenode8(4,j),j=1,4) / 2,3,7,6 /
      data (facenode8(5,j),j=1,4) / 4,8,7,3 /
      data (facenode8(6,j),j=1,4) / 4,1,5,8 /
c
c          set the face nodes in outward normal right-hand order for the
c          20-node brick; corner nodes then edge nodes to match the order
c          of the 8-node quadrillateral; note that the corner numbers are
c          the same as the 8-node brick
c
      data (facenode20(1,j),j=1,8) / 4,3,2,1, 11,10, 9,12 /
      data (facenode20(2,j),j=1,8) / 5,6,7,8, 13,14,15,16 /
      data (facenode20(3,j),j=1,8) / 1,2,6,5,  9,18,13,17 /
      data (facenode20(4,j),j=1,8) / 2,3,7,6, 10,19,14,18 /
      data (facenode20(5,j),j=1,8) / 4,8,7,3, 20,15,19,11 /
      data (facenode20(6,j),j=1,8) / 4,1,5,8, 12,17,16,20 /
c
c          set the face nodes in outward normal right-hand order for the
c          4 and 10 node tetrahedrons; the corners listed first followed
c          by the edge nodes for the 10 node tetrahedron
c
      data (tet_4node(1,j),j=1,3) / 1,3,2 /
      data (tet_4node(2,j),j=1,3) / 1,2,4 /
      data (tet_4node(3,j),j=1,3) / 2,3,4 /
      data (tet_4node(4,j),j=1,3) / 1,4,3 /

      data (tet_10node(1,j),j=1,6) / 1,3,2, 7,6,5 /
      data (tet_10node(2,j),j=1,6) / 1,2,4, 5,9,8 /
      data (tet_10node(3,j),j=1,6) / 2,3,4, 6,10,9 /
      data (tet_10node(4,j),j=1,6) / 1,4,3, 8,10,7 /
c
      debug_level = 10

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'----- begin srt_setlocalfacenodes -----'
        write(logunit,*)'        maxelemface =',maxelemface
        write(logunit,*)'     maxnodeperface =',maxnodeperface
        write(logunit,*)'       nodesperelem =',nodesperelem
        write(logunit,*)'         elemchoice =',elemchoice
        close(unit=logunit)
      end if
c
c          set the element local face node numbers for the element
c          choice; check that the values for numelemface and
c          numnodeperface are correct for the element choice.
c
      select case(elemchoice)
      case(1)      ! 20-node brick
        numelemface = 6
        numnodeperface = 8
        numcornerperface = 4
        if(maxelemface.lt.numelemface)then
          mflag = 2
          message = 'error: maxelemface < 6; expecting 6 element'//
     &      ' faces for the 20-node brick element; cannot set the'//
     &      ' element local face node numbers'//
     &      ' (srt_setlocalfacenodes).'
        else if(numnodeperface.lt.numnodeperface)then
          mflag = 2
          message = 'error: maxnodeperface < 8; expecting 8 nodes'//
     &      ' per face for the 20-node brick element; cannot set'//
     &      ' the element local face node numbers'//
     &      ' (srt_setlocalfacenodes).'
        else if(nodesperelem.ne.20)then
          mflag = 2
          message = 'error: nodesperelem /= 20; expecting 20 nodes'//
     &      ' for the 20-node brick element; cannot set the'//
     &      ' element local face node numbers'//
     &      ' (srt_setlocalfacenodes).'
        end if
        if(mflag.ge.2)goto 900      ! exit on an error

        do i=1,numelemface
          do j=1,numnodeperface
            localfacenode(i,j) = facenode20(i,j)
          end do
        end do

      case(2)      ! 8-node brick
        numelemface = 6
        numnodeperface = 4
        numcornerperface = 4
        if(maxelemface.lt.numelemface)then
          mflag = 2
          message = 'error: maxelemface < 6; expecting 6 element'//
     &      ' faces for the 8-node brick element; cannot set the'//
     &      ' element local face node numbers (srt_setlocalfacenodes).'
        else if(maxnodeperface.lt.numnodeperface)then
          mflag = 2
          message = 'error: maxnodeperface < 4; expecting 4 nodes'//
     &      ' per face for the 8-node brick element; cannot set the'//
     &      ' element local face node numbers (srt_setlocalfacenodes).'
        else if(nodesperelem.ne.8)then
          mflag = 2
          message = 'error: nodesperelem /= 8; expecting 8 nodes'//
     &      ' for the 8-node brick element; cannot set the'//
     &      ' element local face node numbers (srt_setlocalfacenodes).'
        end if
        if(mflag.ge.2)goto 900      ! exit on an error

        do i=1,numelemface
          do j=1,numnodeperface
            localfacenode(i,j) = facenode8(i,j)
          end do
        end do

      case(6)      ! 10-node tetrahedron
        numelemface = 4
        numnodeperface = 6
        numcornerperface = 3
        if(maxelemface.lt.numelemface)then
          mflag = 2
          message = 'error: maxelemface < 4; expecting 4 element'//
     &      ' faces for the 10-node tetrahedron element; cannot set'//
     &      ' the  element local face node numbers'//
     &      ' (srt_setlocalfacenodes).'
        else if(maxnodeperface.lt.numnodeperface)then
          mflag = 2
          message = 'error: maxnodeperface < 6; expecting 6 nodes'//
     &      ' per face for the 10-node tetrahedron element; cannot'//
     &      ' set the  element local face node numbers'//
     &      ' (srt_setlocalfacenodes).'
        else if(nodesperelem.ne.10)then
          mflag = 2
          message = 'error: nodesperelem /= 10; expecting 10 nodes'//
     &      ' for the 10-node tetrahedron element; cannot set the'//
     &      ' element local face node numbers (srt_setlocalfacenodes).'
        end if
        if(mflag.ge.2)goto 900      ! exit on an error

        do i=1,numelemface
          do j=1,numnodeperface
            localfacenode(i,j) = tet_10node(i,j)
          end do
        end do

      case(13)      ! 4-node tetrahedron
        numelemface = 4
        numnodeperface = 3
        numcornerperface = 3
        if(maxelemface.lt.numelemface)then
          mflag = 2
          message = 'error: maxelemface < 4; expecting 4 element'//
     &      ' faces for the 4-node tetrahedron element; cannot set'//
     &      ' the  element local face node numbers'//
     &      ' (srt_setlocalfacenodes).'
        else if(maxnodeperface.lt.numnodeperface)then
          mflag = 2
          message = 'error: maxnodeperface < 3; expecting 3 nodes'//
     &      ' per face for the 4-node tetrahedron element; cannot'//
     &      ' set the  element local face node numbers'//
     &      ' (srt_setlocalfacenodes).'
        else if(nodesperelem.ne.10)then
          mflag = 2
          message = 'error: nodesperelem /= 4; expecting 4 nodes'//
     &      ' for the 4-node tetrahedron element; cannot set the'//
     &      ' element local face node numbers (srt_setlocalfacenodes).'
        end if
        if(mflag.ge.2)goto 900      ! exit on an error

        do i=1,numelemface
          do j=1,numnodeperface
            localfacenode(i,j) = tet_4node(i,j)
          end do
        end do

      case default
        mflag = 2
        message = 'error: unknown element type choice; cannot set'//
     &    ' the element local face node numbers; check the value of'//
     &    ' "elemchoice" (srt_setlocalfacenodes).'
        goto 900      ! exit on an error
      end select      ! elemchoice
c
c jump here on an error
c
900   continue

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'        numelemface =',numelemface
        write(logunit,*)'     numnodeperface =',numnodeperface
        write(logunit,*)'   numcornerperface =',numcornerperface
        write(logunit,*)'   localfacenode():'
        do i=1,numelemface
          write(logunit,'(8i4)')(localfacenode(i,j),j=1,numnodeperface)
        end do
        write(logunit,*)'----- end srt_setlocalfacenodes -----'
        close(unit=logunit)
      end if

      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *            subroutine srt_setwarp3dlocalfacenodes            *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 12/30/02                  *
c     *                                    02/26/03 gvt              *
c     *                                                              *
c     *  provide the list of local face node numbers on each face of *
c     *  the element; get the local element face node numbers which  *
c     *  are the columns in the connectivity array;                  *
c     *  the local element face nodes are given for a warp3d format  *
c     *  element.                                                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine srt_setwarp3dlocalfacenodes(localfacenode,maxelemface,
     &               maxnodeperface,numelemface,numnodeperface,
     &               numcornerperface,nodesperelem,elemchoice,
     &               logflag,logunit,logfile,mflag,message)
c
c          declare variables
c
      implicit none
      integer
     & maxelemface,maxnodeperface,
     & numelemface,numnodeperface,numcornerperface,nodesperelem,
     & elemchoice,
     & localfacenode(maxelemface,maxnodeperface),
     & logflag,logunit,mflag
      character(len=*) :: logfile,message
c
c          variables:
c
c          localfacenode() = integer array to store the local element
c            face node list; for 3-d brick elements there are 6 element
c            faces with 4 or 8 nodes on each element face; return the
c            local face node list on each element face
c          maxelemface = the maximum number of faces per element; row
c            dimension for the localfacenode() array
c          maxnodeperface = the maximum number of nodes per face;
c            column dimension for the localfacenode() array
c          numelemface = return number of element faces, 6 element
c            faces for 3-d brick elements, 4 element faces for 3-d
c            tetrahedron elements
c          numnodeperface = number of nodes on each element face, an
c            8-node brick element has 4 nodes per face, a 20-node
c            brick has 8 nodes per face, a 4-node tetrahedron has 3
c            nodes per face, a 10-node tetrahedron has 6 nodes per face
c          numcornerperface = return the number of corner nodes per
c            element face; the first numcornerperface nodes per row
c            are the face corner nodes
c          nodesperelem = the number of nodes per element; this value
c            must be consistent with the elemchoice value
c          elemchoice = flag to choose the element type;
c            1 = 20-node brick element, also set numelemface = 6 and
c                numnodeperface = 8; expect nodesperelem = 20
c            2 = 8-node brick element, also set numelemface = 6 and
c                numnodeperface = 4; expect nodesperelem = 8
c            6 = 10-node tetrahedron element, also set numelemface = 4
c                and numnodeperface = 6; expect nodesperelem = 10
c           13 = 4-node tetrahedron element, also set numelemface = 4
c                and numnodeperface = 3; expect nodesperelem = 4
c
c          local variables
c
      integer
     & debug_level,i,j,face,nfnode,
     & fnodes(10),tet_10_swap(6),num_face_nodes(13),
     & num_face_corners(13),num_elem_faces(13)
c
c          set the number of nodes per face and the number of corner
c          nodes per face for the various element types;
c          (maximum element type index is 13)
c          elemchoice = flag to choose the element type;
c            1 = 20-node brick element, also set numelemface = 6 and
c                numnodeperface = 8; expect nodesperelem = 20
c            2 = 8-node brick element, also set numelemface = 6 and
c                numnodeperface = 4; expect nodesperelem = 8
c            6 = 10-node tetrahedron element, also set numelemface = 4
c                and numnodeperface = 6; expect nodesperelem = 10
c           13 = 4-node tetrahedron element, also set numelemface = 4
c                and numnodeperface = 3; expect nodesperelem = 4
c
      data tet_10_swap(1:6)       / 1,3,5,2,4,6 /
      data num_face_nodes(1:13)   / 8,4,0,0,0,6,0,0,0,0,0,0,3 /
      data num_face_corners(1:13) / 4,4,0,0,0,3,0,0,0,0,0,0,3 /
      data num_elem_faces(1:13)   / 6,6,0,0,0,4,0,0,0,0,0,0,4 /
c
      debug_level = 10

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'----- begin srt_setwarp3dlocalfacenodes -----'
        write(logunit,*)'        maxelemface =',maxelemface
        write(logunit,*)'     maxnodeperface =',maxnodeperface
        write(logunit,*)'       nodesperelem =',nodesperelem
        write(logunit,*)'         elemchoice =',elemchoice
        close(unit=logunit)
      end if
c
c          set the element local face node numbers for the element
c          choice; check that the values for numelemface and
c          numnodeperface are correct for the element choice.
c
      numelemface      = num_elem_faces(elemchoice)
      numnodeperface   = num_face_nodes(elemchoice)
      numcornerperface = num_face_corners(elemchoice)

      select case(elemchoice)
      case(1:2)      ! 20 and 8 node brick elements
        do face=1,numelemface
          call eqelfn( fnodes, elemchoice, face, nfnode )

          do j=1,nfnode
            localfacenode(face,j) = fnodes(j)
          end do
        end do

      case(6,13)      ! 10 and 4 node tetrahedron elements
        do face=1,numelemface
          call tet_get_nodes( fnodes, elemchoice, face, nfnode )

          select case(elemchoice)
          case(6)      ! 10 node tetrahedron
c
c                   swap fnodes() for corners then mid side nodes;
c                   need face corners first for some routines
c
            do j=1,nfnode
              localfacenode(face,j) = fnodes(tet_10_swap(j))
            end do

          case default
            do j=1,nfnode
              localfacenode(face,j) = fnodes(j)
            end do
          end select
        end do

      case default
        mflag = 2
        message = 'error: unknown element type choice; cannot set'//
     &    ' the element local face node numbers; check the value of'//
     &    ' "elemchoice" (srt_setwarp3dlocalfacenodes).'
        goto 900      ! exit on an error
      end select      ! elemchoice
c
c jump here on an error
c
900   continue

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'        numelemface =',numelemface
        write(logunit,*)'     numnodeperface =',numnodeperface
        write(logunit,*)'   numcornerperface =',numcornerperface
        write(logunit,*)'   localfacenode():'
        do i=1,numelemface
          write(logunit,'(8i4)')(localfacenode(i,j),j=1,numnodeperface)
        end do
        write(logunit,*)'----- end srt_setwarp3dlocalfacenodes -----'
        close(unit=logunit)
      end if

      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                    subroutine srt_sortint                    *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 06/10/04 gvt              *
c     *                                                              *
c     *  sort an integer array in ascending or descending order      *
c     *                                                              *
c     ****************************************************************
c
      subroutine srt_sortint(intarray,dimrow,dimcol,num,col,sort_flag)
c
c         declare variables
c
      implicit none
      integer
     &  dimrow,dimcol,num,col,sort_flag,
     &  intarray(dimrow,dimcol)
c
c         variables:
c         intarray() = integer array of values for sorting
c         dimrow = array row dimension
c         dimcol = array column dimension
c         num = number of rows of values to sort
c         col = number of columns of values to sort
c         sort_flag = choose the sorting direction;
c           1 = sort values in ascending order
c          -1 = sort values in descending order
c
c
c         local variables
c
      integer
     &  i,j,k,swap
c
c         sort the integer values;
c         value tested moves down column one row after the
c         minimum (or maximum) value has been found; look for
c         the next lowest (or highest) value.
c
      if(sort_flag.le.-1)then      ! descending order
        do j=1,num-1
c
c               test each value below "top" of column
c
          do i=j+1,num
            if(intarray(i,col).gt.intarray(j,col))then
c
c               swap row values in every column of the row;
c               keep values on a row together
c
              do k=1,dimcol
                swap = intarray(j,k)
                intarray(j,k) = intarray(i,k)
                intarray(i,k) = swap
              end do
            end if
          end do
        end do

      else      ! ascending order
        do j=1,num-1
c
c               test each value below "top" of column
c
          do i=j+1,num
            if(intarray(i,col).lt.intarray(j,col))then
c
c               swap row values in every column of the row;
c               keep values on a row together
c
              do k=1,dimcol
                swap = intarray(j,k)
                intarray(j,k) = intarray(i,k)
                intarray(i,k) = swap
              end do
            end if
          end do
        end do
      end if
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *               subroutine srt_sortintmultiarray               *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 12/30/02                  *
c     *                                    02/26/03 gvt              *
c     *                                                              *
c     *  sort an integer array in ascending order.                   *
c     *                                                              *
c     *  sort the column, col, of the 2-d integer array; swap rows   *
c     *  in all the columns to keep values on a row together.  sort  *
c     *  rows 1 to num of the array.  sort the secondary logical and *
c     *  real arrays with the same row swapping to match the integer *
c     *  array to keep the rows of all arrays together.              *
c     *                                                              *
c     *  note that the real4list() is single precision               *
c     *                                                              *
c     ****************************************************************
c
      subroutine srt_sortintmultiarray(intarray,dimrow,dimcol,num,col,
     &                   logicallist,real4list,sort_flag)
c
c         declare variables
c
      implicit none
      integer
     &  dimrow,dimcol,num,col,sort_flag,intarray(dimrow,dimcol)
      logical
     &  logicallist(dimrow)
      real
     &  real4list(dimrow)
c
c         variables:
c         intarray() = integer array of values for sorting
c         dimrow = array row dimension
c         dimcol = array column dimension
c         num = number of rows of values to sort
c         col = number of columns of values to sort
c         logicallist() = array of logical values; keep the rows
c           together with any changes to the intarray() during sorting
c         real4list() = array of real values; keep the rows
c           together with any changes to the intarray() during sorting
c         sort_flag = choose the sorting direction;
c           1 = sort values in ascending order
c          -1 = sort values in descending order
c
c
c         local variables
c
      integer
     &  i,j,k,swap
      logical
     &  temp_logical
      real
     &  temp_real
c
c         value tested moves down column one row after the
c         minimum value has been found.  look for the next
c         lowest value.
c
c
c         sort the integer values;
c         value tested moves down column one row after the
c         minimum (or maximum) value has been found; look for
c         the next lowest (or highest) value.
c
      if(sort_flag.le.-1)then      ! descending order
        do j=1,num-1
c
c               test each value below "top" of column
c
          do i=j+1,num
            if(intarray(i,col).gt.intarray(j,col))then
c
c               swap row values in every column of the row;
c               keep values on a row together
c
              do k=1,dimcol
                swap = intarray(j,k)
                intarray(j,k) = intarray(i,k)
                intarray(i,k) = swap
              end do
c
c               swap rows of other secondary arrays
c
              temp_logical = logicallist(j)
              logicallist(j) = logicallist(i)
              logicallist(i) = temp_logical

              temp_real = real4list(j)
              real4list(j) = real4list(i)
              real4list(i) = temp_real
            end if
          end do
        end do

      else      ! ascending order
        do j=1,num-1
c
c               test each value below "top" of column
c
          do i=j+1,num
            if(intarray(i,col).lt.intarray(j,col))then
c
c               swap row values in every column of the row;
c               keep values on a row together
c
              do k=1,dimcol
                swap = intarray(j,k)
                intarray(j,k) = intarray(i,k)
                intarray(i,k) = swap
              end do
c
c               swap rows of other secondary arrays
c
              temp_logical = logicallist(j)
              logicallist(j) = logicallist(i)
              logicallist(i) = temp_logical

              temp_real = real4list(j)
              real4list(j) = real4list(i)
              real4list(i) = temp_real
            end if
          end do
        end do
      end if
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *  compute the length (the magnitude) of the vector, a(dim).   *
c     *                                                              *
c     ****************************************************************
c
      function srt_veclength(a,dim)
c
c         declare variables
c
      implicit none
      integer dim
      double precision
     &  srt_veclength,a(dim)
c
c         local variables
c
      integer i
      double precision
     &  sum, zero
        data zero / 0.0d00 /
c
c         compute the vector length
c
      sum = zero

      if( dim .eq. 3 )then
        sum = a(1)**2 + a(2)**2 + a(3)**2
      else
        do i= 1, dim
          sum = sum + a(i)**2
        end do
      end if

      srt_veclength = sqrt(sum)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *  compute the cartesian distance between two points:          *
c     *  a() and b(); the 3-d distance = sqrt(dx^2 + dy^2 + dz^2)    *
c     *  where dx = b(1) - a(1)                                      *
c     *        dy = b(2) - a(2)                                      *
c     *        dz = b(3) - a(3)                                      *
c     *                                                              *
c     ****************************************************************
c
      function srt_dist3d(a,b,dim)
c
c           declare variables
c
      implicit none
      integer dim
      double precision
     &  srt_dist3d,a(dim),b(dim)
c
c           local variables
c
      integer i
      double precision
     &  sum
c
c           initialize values
c
      sum = 0.0d00
c
c           compute the distance between the points
c
      if(dim.eq.3)then
        sum = (b(1)-a(1))**2 + (b(2)-a(2))**2 + (b(3)-a(3))**2
      else
        do i=1,dim
          sum = sum + (b(i) - a(i))**2
        end do
      end if

      srt_dist3d = sqrt(sum)

      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *  compute the location of the midpoint, mid, between the 2    *
c     *  given points, a, b.                                         *
c     *                                                              *
c     ****************************************************************
c
      subroutine srt_midpoint(a,b,mid,dim)
c
c         declare variables
c
      implicit none
      integer dim
      double precision
     &  a(dim), b(dim), mid(dim), two
          data two / 2.0d00 /
c
c         local variables
c
      integer i
c
c         dim = number of dimensions, expecting dim=3 for 3-d x,y,z
c               coordinates
c         a(),b() = given points
c         mid() = new point located at the midpoint between points
c               a() and b()
c
      do i=1,dim
        mid(i) = (a(i) + b(i)) / two
      end do
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *  normalize the given vector to get a unit vector, nhat.      *
c     *                                                              *
c     ****************************************************************
c
      subroutine srt_normalize(nhat,dim)
c
c         declare variables
c
      implicit none
      integer dim
      double precision
     &  nhat(dim)
c
c         nhat = vector to normalize
c         dim = vector dimensions, expecting dim=3 for x,y,z coordinates
c
c
c         local variables
c
      integer i
      double precision
     &  mag,zero
      data zero / 0.0d0 /
c
c         declare functions
c
      double precision
     &  srt_veclength
c
c            normalize the vector,
c            get the vector magnitude
c
      mag = srt_veclength(nhat,dim)
c
c            normalize by the magnitude
c
      if(mag.ne.zero)then
        do i=1,dim
          nhat(i) = nhat(i)/mag
        end do
      end if
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *  write the error or warning message to the debugging         *
c     *  log file.                                                   *
c     *                                                              *
c     ****************************************************************
c
      subroutine srt_writelogmessage(mflag,message,logfile,logunit)
c
c         declare variables
c
      implicit none
      integer mflag,logunit
      character(len=*) :: message,logfile
c
      open(file=logfile,unit=logunit,position='append')
      write(logunit,*)
      write(logunit,*)'  mflag =',mflag
      write(logunit,'(a)')trim(message)
      write(logunit,*)
      close(unit=logunit)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine tied_nodecheckredundant               *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 06/10/04                  *
c     *                                                              *
c     *  check for redundant slave nodes to avoid generating         *
c     *  redundant mpc equations; if a slave node is also used by    *
c     *  a master element then a mpc equation connecting the slave   *
c     *  node to itself on the master surface is created but is not  *
c     *  needed since the master and slave surfaces are already      *
c     *  connected by the common node; compare the slave node id     *
c     *  number to the node id numbers used by the master elements   *
c     *  on the master surface                                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine tied_nodecheckredundant(tiednodeid,tiednodeadjustflag,
     &    tiednodegaptolerance,maxtiednode,maxtieddatacol,numtiednode,
     &    masternodelist,max_masternodelist,num_masternodelist,out,
     &    logflag,logunit,logfile,mflag,message)
c
c           declare variables
c
      implicit none
      integer
     &  maxtiednode,maxtieddatacol,numtiednode,
     &  max_masternodelist,num_masternodelist,out,
     &  logflag,logunit,mflag,i,num,rem,
     &  tiednodeid(maxtiednode,maxtieddatacol),
     &  masternodelist(max_masternodelist),
     &  red_nod(maxtiednode)

      real
     &  tiednodegaptolerance(maxtiednode),mlt

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
c           numtiednode = return the number of tied slave nodes
c           masternodelist() = list of the node id numbers on the master
c             surface
c           max_masternodelist = array row dimension
c           num_masternodelist = number of master node id numbers
c           out = the output file unit number
c
c           logflag,logunit,logfile = log file for debugging and testing;
c             the flag, file unit, and file name for the log file
c           mflag,message = error or warning message; integer flag
c             0=no message, 1=warning, 2=error; message text string
c
c
c           local variables
c
      integer debug_level,row,sort_flag,node_id,master_row
      integer remove_count
c
      debug_level = 7

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'----- begin tied_nodecheckredundant -----'
        write(logunit,*)'          numtiednode =',numtiednode
        write(logunit,*)'          maxtiednode =',maxtiednode
        write(logunit,*)'       maxtieddatacol =',maxtieddatacol
        write(logunit,*)'   max_masternodelist =',max_masternodelist
        write(logunit,*)'   num_masternodelist =',num_masternodelist
        write(logunit,*)'       file unit, out =',out
        close(unit=logunit)
      end if
c
c           error checking
c
      if(numtiednode.le.0)then
        mflag = 2
        message = 'error: numtiednode = 0 (tied_nodecheckredundant).'
      else if(maxtiednode.le.0)then
        mflag = 2
        message = 'error: maxtiednode = 0 (tied_nodecheckredundant).'
      else if(maxtieddatacol.le.0)then
        mflag = 2
        message = 'error: maxtieddatacol = 0 (tied_nodecheckredundant).'
      else if(max_masternodelist.le.0)then
        mflag = 2
        message = 'error: max_masternodelist = 0'//
     &    ' (tied_nodecheckredundant).'
      else if(num_masternodelist.le.0)then
        mflag = 2
        message = 'error: num_masternodelist = 0'//
     &    ' (tied_nodecheckredundant).'
      end if
      if(mflag.ge.2)goto 900
c
c           initialize values
c
      remove_count = 0
c           sort the master node id list by ascending order
c
      sort_flag = 1      ! 1=ascending order, -1=descending order

      call srt_sortint(masternodelist,max_masternodelist,1,
     &  num_masternodelist,1,sort_flag)

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'sorted master node id list'
        write(logunit,*)'  num_masternodelist =',num_masternodelist
        write(logunit,*)'  max_masternodelist =',max_masternodelist
        write(logunit,*)'           sort_flag =',sort_flag
        write(logunit,*)
        write(logunit,3)'row','master_node_id'
        do master_row=1,num_masternodelist
          write(logunit,4)master_row,masternodelist(master_row)
        end do
        close(unit=logunit)
      end if
c
c           check for slave node id numbers also used by the master
c           elements; remove the redundant slave nodes from all of
c           the tied node arrays by setting the slave node id number
c           to zero; the array will be sorted later to remove the
c           zeroed out rows;
c           when the master node id is larger than the slave node id
c           the inner loop is exited since there can't be any matches
c           between the slave and master node id numbers;
c           write a warning message to the output file when a slave
c           node is removed
c
      loop_1: do row=1,numtiednode
        node_id = tiednodeid(row,1)

        loop_2: do master_row=1,num_masternodelist
          if(node_id.eq.masternodelist(master_row))then
            remove_count = remove_count + 1
            red_nod(remove_count) = node_id
            tiednodeid(row,1:maxtieddatacol) = 0
            tiednodeadjustflag(row) = .false.
            tiednodegaptolerance(row) = 0.0

            mflag = 1
            message = 'warning: a slave node id number is also used'//
     &        ' by the master surface elements; the redundant slave'//
     &        ' node has been removed from the tied node list;'//
     &        ' check the log file (tied_nodecheckredundant).'

            if(logflag.ge.1)then
              open(unit=logunit,file=logfile,position='append')
              write(logunit,*)
              write(logunit,*)TRIM(message)
              write(logunit,*)'                   mflag =',mflag
              write(logunit,*)'   slave node loop index =',row
              write(logunit,*)'    redundant node count =',remove_count
              write(logunit,*)'           slave node id =',node_id
              write(logunit,*)'  master node loop index =',master_row
              write(logunit,*)'          master node id =',
     &                                   masternodelist(master_row)
              write(logunit,*)
              close(unit=logunit)
            end if

          else if(masternodelist(master_row).gt.node_id)then
            exit loop_2
          end if
        end do loop_2
      end do loop_1

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'number of redundant slave node id numbers:'
        write(logunit,*)'        remove_count =',remove_count
        close(unit=logunit)
      end if
c
c           update the tied node arrays if any slave node id numbers
c           were removed
c
      if(remove_count.gt.0)then

         write(out,9001) remove_count
         mlt = real(remove_count)/8.0
         num = int(mlt)*8
         rem = remove_count - num
         do i = 1, num, 8
            write(out,9004) red_nod( i ),red_nod(i+1),red_nod(i+2),
     &                      red_nod(i+3),red_nod(i+4),red_nod(i+5),
     &                      red_nod(i+6),red_nod(i+7)
         end do
         if (rem .eq. 1)  write(out,9005) red_nod(num+1)
         if (rem .eq. 2)  write(out,9006) red_nod(num+1),red_nod(num+2)
         if (rem .eq. 3)  write(out,9007) red_nod(num+1),red_nod(num+2),
     &                                    red_nod(num+3)
         if (rem .eq. 4)  write(out,9008) red_nod(num+1),red_nod(num+2),
     &                                    red_nod(num+3),red_nod(num+4)
         if (rem .eq. 5)  write(out,9009) red_nod(num+1),red_nod(num+2),
     &                                    red_nod(num+3),red_nod(num+4),
     &                                    red_nod(num+5)
         if (rem .eq. 6)  write(out,9010) red_nod(num+1),red_nod(num+2),
     &                                    red_nod(num+3),red_nod(num+4),
     &                                    red_nod(num+5),red_nod(num+6)
         if (rem .eq. 7)  write(out,9011) red_nod(num+1),red_nod(num+2),
     &                                    red_nod(num+3),red_nod(num+4),
     &                                    red_nod(num+5),red_nod(num+6),
     &                                    red_nod(num+7)
         write(out,9003)

c
c           sort the tied node arrays by descending node id number to
c           get the zeroed rows at the bottom of the list
c
        sort_flag = -1    ! 1=ascending order, -1=descending order

        call srt_sortintmultiarray(tiednodeid,maxtiednode,
     &               maxtieddatacol,numtiednode,1,
     &               tiednodeadjustflag,tiednodegaptolerance,sort_flag)
c
c           update the reduced number of slave nodes;
c           sort the slave node id numbers by ascending order
c
        numtiednode = numtiednode - remove_count
        sort_flag = 1     ! 1=ascending order, -1=descending order

        call srt_sortintmultiarray(tiednodeid,maxtiednode,
     &               maxtieddatacol,numtiednode,1,
     &               tiednodeadjustflag,tiednodegaptolerance,sort_flag)
      end if
c
c           debugging output of the computed slave node coordinates
c
      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)
        write(logunit,*)'after checking for redundant slave nodes'
        write(logunit,*)'      remove_count =',remove_count
        write(logunit,*)'       numtiednode =',numtiednode
        write(logunit,*)
        write(logunit,1)'row','node_id','elem_id','face_id','flag',
     &                  'adjust_flag','gap_tol'
        do row=1,numtiednode
          write(logunit,2)row,tiednodeid(row,1:maxtieddatacol),
     &              tiednodeadjustflag(row),
     &              tiednodegaptolerance(row)
        end do
        close(unit=logunit)
      end if
1     format(a8,3a10,a5,a16,a16)
2     format(i8,3i10,i5,l16,es16.6)
3     format(a8,a18)
4     format(i8,i18)

9001  format(/,
     & 1x,'>> warning: the following',i9,' nodes have been removed',/,
     & 1x,'>>          from the tied mesh slave node list to prevent',/,
     & 1x,'>>          redundant mpc equations from being generated',/)
9003  format(/,
     & 1x,'>> info: redundant slave nodes are caused when elements',/,
     & 1x,'>>       in the slave and master surface element lists',/,
     & 1x,'>>       use the same nodes; for example, two surfaces',/,
     & 1x,'>>       may have some merged nodes',/)
 9004 format(8(i9))
 9005 format(i9)
 9006 format(2(i9))
 9007 format(3(i9))
 9008 format(4(i9))
 9009 format(5(i9))
 9010 format(6(i9))
 9011 format(7(i9))
c
c----------------------------------------------------------------------
c           jump here on an error for exit procedures
c
900   continue

      if(logflag.ge.debug_level)then
        open(unit=logunit,file=logfile,position='append')
        write(logunit,*)'----- end tied_nodecheckredundant -----'
        close(unit=logunit)
      end if

      return
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine srt_mpcs_resize                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 9/3/2017 rhd               *
c     *                                                              *
c     *     increase size of the tied_con_mpc_table.                 *
c     *     make larger one, deep copy old -> new, move allocation.  *
c     *     release old table                                        *
c     *                                                              *
c     ****************************************************************
c

      subroutine srt_mpcs_resize
      use global_data ! old common.main
      use mod_mpc, only : tied_con_mpc_table, mpc_eqn
      implicit none
c
      integer :: old_mpc_size, i, j, nt
      logical, parameter :: local_debug = .false.
      type (mpc_eqn), allocatable, dimension (:) :: new_mpc_table
c
      if( local_debug ) then
         write(*,*) '...  resizing mpc tied table.....'
         write(*,*) '       old_max_mpc_tied: ', max_mpc_tied
      end if
c
      old_mpc_size = max_mpc_tied
      max_mpc_tied = 2 * old_mpc_size
      allocate( new_mpc_table(max_mpc_tied) )
c
c              deep copy as required
c
      do i = 1, old_mpc_size
          nt = tied_con_mpc_table(i)%num_terms
          new_mpc_table(i)%num_terms = nt
          new_mpc_table(i)%constant = tied_con_mpc_table(i)%constant
          allocate( new_mpc_table(i)%node_list(nt),
     &              new_mpc_table(i)%dof_list(nt),
     &              new_mpc_table(i)%multiplier_list(nt)  )
          do j = 1, nt
           new_mpc_table(i)%node_list(j) =
     &                   tied_con_mpc_table(i)%node_list(j)
           new_mpc_table(i)%dof_list(j) =
     &                   tied_con_mpc_table(i)%dof_list(j)
           new_mpc_table(i)%multiplier_list(j) =
     &                   tied_con_mpc_table(i)%multiplier_list(j)
          end do ! on j
      end do
c
c              initialize remainder of new table. probably not req'd
c
      do i =  old_mpc_size+1, max_mpc_tied
          new_mpc_table(i)%num_terms = 0
          new_mpc_table(i)%constant = 0.0
          new_mpc_table(i)%node_list => null()
          new_mpc_table(i)%dof_list  => null()
          new_mpc_table(i)%multiplier_list  => null()
      end do ! on i
c
c              release old table and move allocation. F2003 &
c              later does a deep release on derived types. move_alloc
c              does a deep release on new_mpc_table
c
      deallocate( tied_con_mpc_table )
      call move_alloc( new_mpc_table, tied_con_mpc_table )
c
      return
      end
c
