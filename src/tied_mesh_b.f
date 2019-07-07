

c     ****************************************************************
c     *                                                              *
c     *                    subroutine tied_nodeiso                   *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 11/04/02                  *
c     *                                    7/5/2019 rhd              *
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
      integer ::
     &  maxtiednode,maxtieddatacol,maxisodof,numtiednode,
     &  dimnode,dimcoord,numnode,dimelem,maxnodeperelem,numelem,
     &  meshformat,out,
     &  logflag,logunit,mflag,
     &  tiednodeid(maxtiednode,maxtieddatacol),
     &  elemcon(dimelem,maxnodeperelem)
      double precision ::
     &  tiednodeisocoord(maxtiednode,maxisodof),
     &  tiednodeglobalcoord(maxtiednode,dimcoord),
     &  nodecoord(dimnode,dimcoord)
      real :: tiednodegaptolerance(maxtiednode), props(mxelpr,*)
      logical :: tiednodeadjustflag(maxtiednode)
      character(len=*) :: logfile,message
c
c           local variables
c
      integer, parameter :: maxdof = 3
      integer, parameter :: maxtetdof = 4
      integer, parameter :: maxelemface = 6
      integer, parameter :: maxnodeperface = 9
      integer, parameter :: maxiterloop = 1000
      integer ::
     &  debug_level,allocate_error,i,j,nd,row,enode,iout,dumi2,
     &  numelemface,numnodeperface,numcornerperface,
     &  node_id,elem_id,face_id,node_per_elem,elemchoice,
     &  constantisodof,itercount,node_row,corner,
     &  min_iter_count,max_iter_count,min_iter_nodeid,max_iter_nodeid,
     &  xi_error_count,eta_error_count,zeta_error_count,
     &  check_iso_coord_flag,ierr,iword,err,gap_cnt,cnt,dumi,
     &  localfacenode(maxelemface,maxnodeperface)
      integer, allocatable, dimension (:) :: gap_nodes,gap_elems,
     &                                       gap_faces
      logical :: local_debug, twod_flag, iso_adjust_flag, debug_now,
     &           tied_node_is_a_master_node
c
      double precision ::  
     &  constantisocoord,dist_point,epsilon,angle_point,
     &  epsilon_angle,delta1,delta2,delta_iso,scale1,scale2,
     &  iso_error_val,dist_gap,
     &  isopoint(maxdof),globalpoint(maxdof),position(maxdof),
     &  normal(maxdof),tangent1(maxdof),tangent2(maxdof),
     &  jacobian(maxdof,maxdof),invj(maxdof,maxdof),detj,
     &  isopoint_old(maxdof),slavenode(maxdof),vec(maxdof),
     &  isotetpoint(maxtetdof),vec1(maxdof),vec2(maxdof),
     &  dsf(32,3),coord(3,32), dumd, character_size
      double precision, allocatable :: shapefunc(:),
     &  shapederivs(:,:), elemcoord(:,:), gap_dists(:), gap_toler(:)
      double precision, parameter :: zero = 0.d0, one = 1.d0,
     &                               three = 3.d0,
     &                               pi = 4.0d0*atan(1.0d0)
      real :: rword, dumr
      character(len=1) :: dums
c
c           declare exgternal functions
c
      double precision, external :: srt_dist3d, srt_getangle, 
     &                              srt_dotprod
c
      equivalence (rword,iword)
c
      debug_level = 7
      local_debug = .false.
      debug_now   = local_debug .and. logflag >= debug_level
      call iodevn( dumi, iout, dumi2, 1 )
      if( logflag >= debug_level ) call tied_nodeiso_a
c
c           error checking
c
      call tied_nodeiso_b
      if( mflag >=2 ) goto 900
c
c           initialize values
c
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

      call tied_nodeiso_bb !  allocates w/checks
      if( mflag >= 2 ) go to 900
      if( logflag >= debug_level ) call tied_nodeiso_c
c
c           locate the slave node on its assigned master element face;
c           compute the isoparametric element coordinates at the slave
c           node location
c
      do nd = 1, numtiednode    !  *** loop over all slave nodes ****
        node_id = tiednodeid(nd,1)      ! slave node id
        elem_id = tiednodeid(nd,2)      ! assigned master element id
        face_id = tiednodeid(nd,3)      ! assigned master element face
        if( debug_now ) call tied_nodeiso_d
        call tied_nodeiso_e ! error chceking
        if( mflag >=2 ) goto 900
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
          if( mflag >=2 ) goto 900
        case(2)      ! warp3d
          call srt_setwarp3dlocalfacenodes(localfacenode,
     &               maxelemface,maxnodeperface,numelemface,
     &               numnodeperface,numcornerperface,node_per_elem,
     &               elemchoice,
     &               logflag,logunit,logfile,mflag,message)
              if( mflag >=2 ) goto 900
        case default
          mflag = 2
          message = 'error: unexpected meshformat choice;'//
     &              ' cannot get local element face nodes'//
     &              ' (tied_nodeiso).'
        end select      ! meshformat
        if( mflag >=2 ) goto 900
c
c               set the isoparametric dof that is constant on the master
c               element face; set the constant isoparametric value
c
        call tied_isoconstdof(constantisodof,constantisocoord,
     &          face_id,node_per_elem,elemchoice,meshformat,
     &          logflag,logunit,logfile,mflag,message)
        if( mflag >=2 ) goto 900
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
          do i = 1, maxdof
            isopoint(i) = isotetpoint(i+1)
          end do
c
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
c               check that element incidences are within bounds
c               at this point (sanity check).
c               put the element node coordinates into the local
c               elemcoord() array; used to compute the jacobian matrix
c
        call tied_nodeiso_dd
        if( mflag>= 2 ) go to 900
c
        do enode = 1, node_per_elem
          node_row = elemcon(elem_id,enode)
          elemcoord(enode,1:maxdof) = nodecoord(node_row,1:maxdof)
        end do
        if( debug_now ) call tied_nodeiso_f
c
c               get a characteristic element size to support tolerance
c               checking
c
        call tied_nodeiso_gg
c
c               is the slave node actually a corner node on the
c               master element. if so, we do not want to allow
c               the subsequent gap type adjustment of the master
c               node coordinates. the iterative procedure below
c               seems to work ok in finding that the slave node
c               actually coincides with a master element node. no
c               need for special cases.    
c
        call tied_nodeiso_ee ! set a flag
        if( debug_now ) call tied_nodeiso_ff
c
c               iteration loop to project the slave node onto the master
c               element face and compute the local face isoparametric
c               coordinates of the slave node's position
c
        if( .not. tied_node_is_a_master_node ) call tied_nodeiso_hh
        if( mflag >= 2 ) go to 900
c

        if( itercount < min_iter_count ) then
          min_iter_count = itercount
          min_iter_nodeid = node_id
        end if
        if( itercount > max_iter_coun t) then
          max_iter_count = itercount
          max_iter_nodeid = node_id
        end if

        if( itercount >= maxiterloop ) then
          mflag = 1
          message = 'warning: maximum iterations used to locate the'//
     &              ' slave node on the master element face; check'//
     &              ' the isoparametric element face coordinates'//
     &              ' (tied_nodeiso).'
          if( logflag >= 1 ) call srt_writelogmessage(mflag,
     &          message,logfile,logunit)
        end if
c
c 
c
c               check if the face point isoparametric coordinates are
c               outside the +/- 1.0 range; expecting xi, eta, zeta to
c               be within +/- iso_error_val to allow for some round off;
c               can also reset the isoparametric value to be within the
c               element face, set the xi,eta,zeta value to be +/- 1.0
c               tied_nodeiso_cc can change isopoint(1:3) values
c
        iso_adjust_flag = .false.
        if( check_iso_coord_flag > 0 ) call tied_nodeiso_cc
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
        if( (iso_adjust_flag .and. tiednodeadjustflag(nd)) .or.
     &     tiednodegaptolerance(nd) > zero ) then

          call shapef(elemchoice,isopoint(1),isopoint(2),isopoint(3),
     &                shapefunc)

          call tied_globalcoordfromiso(globalpoint,maxdof,
     &                shapefunc,maxnodeperelem,node_per_elem,
     &                elemcoord,maxnodeperelem,maxisodof,
     &                logflag,logunit,logfile,mflag,message)
          if( mflag >=2 ) goto 900
           end if
!        end if
c
c               check if the distance between the element face point and
c               the slave node is within the given gap tolerance, set the
c               tied node flag = 1, otherwise the node is too far from
c               the master surface, set the tied node flag = 0
c
        if( tiednodegaptolerance(nd) > zero )then
          dist_gap = srt_dist3d(globalpoint,slavenode,maxdof)
          if( dist_gap <= tiednodegaptolerance(nd) ) then
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
            if( logunit .ne. out ) call tied_nodeiso_w
c
c               report the slave node warning to the debugging log file
c
            if( logflag >= 1 ) call tied_nodeiso_x
          end if
        end if

        if( tiednodeadjustflag(nd) ) then 
          tiednodeglobalcoord(nd,1:maxisodof) = globalpoint(1:maxdof)
          nodecoord(node_id,1:maxdof) = globalpoint(1:maxdof)
        end if

      end do  !      end loop over slave nodes ****

      if( gap_cnt > 0 ) call tied_nodeiso_z ! write-deallocate
c
c           write slave node information to the debugging log file
c
      if( logflag >= debug_level ) call tied_nodeiso_y
c
c----------------------------------------------------------------------
c           jump here on an error for exit procedures
c           =========================================
c
900   continue
c
c           deallocate local arrays
c
      if( allocated(shapefunc) ) deallocate(shapefunc)
      if( allocated(shapederivs) ) deallocate(shapederivs)
      if( allocated(elemcoord) ) deallocate(elemcoord)
c
      if( logflag >= debug_level ) call tied_nodeiso_aa

      return
c
      contains
c     ========
c
c
c 
c           debuuging output, logout support to cleanup code
c           ================================================

      subroutine tied_nodeiso_a
c
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
c
      return
      end subroutine tied_nodeiso_a

      subroutine tied_nodeiso_b
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
c
      return
      end subroutine tied_nodeiso_b

      subroutine tied_nodeiso_c
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)'local values:'
      write(logunit,*)'                   pi =',pi
      write(logunit,*)'              epsilon =',epsilon
      write(logunit,*)'        epsilon_angle =',epsilon_angle
      write(logunit,*)'        iso_error_val =',iso_error_val
      write(logunit,*)' check_iso_coord_flag =',check_iso_coord_flag
      close(unit=logunit)
c
      return
      end subroutine tied_nodeiso_c


      subroutine tied_nodeiso_d
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)
      write(logunit,*)'  slave node row index, nd =',nd
      write(logunit,*)'             slave node_id =',node_id
      write(logunit,*)'            master elem_id =',elem_id
      write(logunit,*)'    master element face_id =',face_id
      close(unit=logunit)
c
      return
      end subroutine tied_nodeiso_d

      subroutine tied_nodeiso_e
c
      if(node_id.le.0)then
        mflag = 2
        message = 'error: node_id = 0; invalid node id for the'//
     &            ' slave node; check the tiednodeid() array'//
     &            ' (tied_nodeiso).'
      else if(node_id.gt.numnode)then
        mflag = 2
        message = 'error: node_id > numnode; invalid node id for'//
     &            ' the slave node; check the tiednodeid() array'//
     &            ' (tied_nodeiso).'
      else if(elem_id.le.0)then
        mflag = 2
        message = 'error: elem_id = 0; invalid element id for'//
     &            ' the master element; check the tiednodeid()'//
     &            ' array (tied_nodeiso).'
      else if(elem_id.gt.numelem)then
        mflag = 2
        message = 'error: elem_id > numelem; invalid element id for'//
     &            ' the master element; check the tiednodeid()'//
     &            ' array (tied_nodeiso).'
      else if(face_id.le.0)then
        mflag = 2
        message = 'error: face_id = 0; invalid element face id for'//
     &            ' the master element; check the tiednodeid()'//
     &            ' array (tied_nodeiso).'
      end if
c
      return
      end subroutine tied_nodeiso_e

      subroutine tied_nodeiso_f
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)'     node_per_elem =',node_per_elem
      write(logunit,*)'        elemchoice =',elemchoice
      write(logunit,*)'    constantisodof =',constantisodof
      write(logunit,3)'        isopoint() =',isopoint(1:maxdof)
      write(logunit,3)'       slavenode() =',slavenode(1:maxdof)
      close(unit=logunit)
c
      return
 3    format(a,3es16.6)
      end subroutine tied_nodeiso_f


      subroutine tied_nodeiso_g
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)
      write(logunit,*)'    itercount =',itercount
      write(logunit,3)'     isopoint() =',isopoint(1:maxdof)
      write(logunit,3)'    slavenode() =',slavenode(1:maxdof)
      close(unit=logunit)
c
      return
 3    format(a,3es16.6)
      end subroutine tied_nodeiso_g



      subroutine tied_nodeiso_h
c
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
c
      return
 3    format(a,3es16.6)
      end subroutine tied_nodeiso_h


      subroutine tied_nodeiso_i
c
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
c
      return
      end subroutine tied_nodeiso_i


      subroutine tied_nodeiso_j
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)
      write(logunit,*)'face point location converged'//
     &                ' (small iso change)'
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
c
      return
      end subroutine tied_nodeiso_j

      subroutine tied_nodeiso_k
c
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
c
      return
      end subroutine tied_nodeiso_k

      subroutine tied_nodeiso_l
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)'  numcornerperface =',numcornerperface
      write(logunit,*)'            scale1 =',scale1
      write(logunit,*)'            scale2 =',scale2
      close(unit=logunit)
c
      return
      end subroutine tied_nodeiso_l


      subroutine tied_nodeiso_m
c
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
c
      return
 3    format(a,3es16.6)
      end subroutine tied_nodeiso_m

      subroutine tied_nodeiso_n
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)'update location for tetrahedron'
      write(logunit,*)'        delta1 =',delta1
      write(logunit,*)'        delta2 =',delta2
      write(logunit,3)'         vec1() =',vec1(1:maxdof)
      write(logunit,3)'         vec2() =',vec2(1:maxdof)
      write(logunit,3)'  globalpoint() =',globalpoint(1:maxdof)
      close(unit=logunit)
c
      return
 3    format(a,3es16.6)
      end subroutine tied_nodeiso_n

      subroutine tied_nodeiso_o
c
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
     &                      globalpoint(1:maxdof)
      write(logunit,3)'       slavenode() =',
     &                        slavenode(1:maxdof)
      write(logunit,*)
      close(unit=logunit)
c
      return
 3    format(a,3es16.6)
      end subroutine tied_nodeiso_o


      subroutine tied_nodeiso_p
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)'        delta1 =',delta1
      write(logunit,*)'        delta2 =',delta2
      write(logunit,*)'updated face point position:'
      write(logunit,*)'   isopoint(1) =',isopoint(1)
      write(logunit,*)'   isopoint(2) =',isopoint(2)
      write(logunit,*)'   isopoint(3) =',isopoint(3)
      close(unit=logunit)
c
      return
      end subroutine tied_nodeiso_p
c


      subroutine tied_nodeiso_q
c
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
c
      return
      end subroutine tied_nodeiso_q
c

      subroutine tied_nodeiso_r
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)
      write(logunit,'(a)')'note: set xi = limit value'
      write(logunit,*)'   xi = isopoint(1) =',isopoint(1)
      write(logunit,*)'   xi = reset to    =',
     &                sign(one,isopoint(1))
      write(logunit,*)'      iso_error_val =',iso_error_val
      write(logunit,*)'      slave node_id =',node_id
      write(logunit,*)'     xi_error_count =',xi_error_count
      write(logunit,*)
      close(unit=logunit)
c
      return
      end subroutine tied_nodeiso_r


      subroutine tied_nodeiso_s
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)
      write(logunit,*)'         mflag =',mflag
      write(logunit,'(a)')trim(message)
      write(logunit,*)'  eta = isopoint(2) =',isopoint(2)
      write(logunit,*)'      iso_error_val =',iso_error_val
      write(logunit,*)'      slave node_id =',node_id
      write(logunit,*)'    eta_error_count =',eta_error_count
      write(logunit,*)
c
      end subroutine tied_nodeiso_s

      subroutine tied_nodeiso_t
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)
      write(logunit,'(a)')'note: set eta = limit value'
      write(logunit,*)'  eta = isopoint(2) =',isopoint(2)
      write(logunit,*)'  eta = reset to    =',
     &                sign(one,isopoint(2))
      write(logunit,*)'      iso_error_val =',iso_error_val
      write(logunit,*)'      slave node_id =',node_id
      write(logunit,*)'    eta_error_count =',eta_error_count
      write(logunit,*)
      close(unit=logunit)
c
      end subroutine tied_nodeiso_t


      subroutine tied_nodeiso_u
c
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
c
      end subroutine tied_nodeiso_u

      subroutine tied_nodeiso_v
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)
      write(logunit,'(a)')'note: set zeta = limit value'
      write(logunit,*)' zeta = isopoint(3) =',isopoint(3)
      write(logunit,*)' zeta = reset to    =',
     &                sign(one,isopoint(3))
      write(logunit,*)'      iso_error_val =',iso_error_val
      write(logunit,*)'      slave node_id =',node_id
      write(logunit,*)'   zeta_error_count =',zeta_error_count
      write(logunit,*)
      close(unit=logunit)
c
      end subroutine tied_nodeiso_v

      subroutine tied_nodeiso_w
c
      if( gap_cnt == 0 ) then 
          allocate( gap_nodes(maxtiednode),
     &              gap_elems(maxtiednode),
     &              gap_faces(maxtiednode),
     &              gap_dists(maxtiednode),
     &              gap_toler(maxtiednode), stat=err)
          if( err .ne. 0 ) then
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
c
       return
c
      end subroutine tied_nodeiso_w

      subroutine tied_nodeiso_x
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)
      write(logunit,'(a)')trim(message)
      write(logunit,*)'  message flag, mflag         =',mflag
      write(logunit,*)'  dependent slave node id     =',node_id
      write(logunit,*)'  assigned master element id  =',elem_id
      write(logunit,*)'  master element face id      =',face_id
      write(logunit,*)'  node gap distance           =',dist_gap
      write(logunit,*)'  given gap tolerance         =',
     &                   tiednodegaptolerance(nd)
      write(logunit,*)'  slave node row, nd          =',nd
      write(logunit,*)'  tied node flag (set = 0)    =',
     &                   tiednodeid(nd,4)
      write(logunit,*)'  numtiednode                 =',
     &                   numtiednode
      write(logunit,*)
      close(unit=logunit)
      return
c
      end subroutine tied_nodeiso_x

      subroutine tied_nodeiso_y
c
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
     &                 'xi_iso1','eta_iso2','zeta_iso3'
      do row = 1, numtiednode
        write(logunit,4)row,tiednodeid(row,1:maxtieddatacol),
     &                   tiednodeisocoord(row,1:maxisodof)
      end do
      if( local_debug )then
        write(logunit,*)
        write(logunit,*)'tied slave node updated global location'
        write(logunit,*)'      numtiednode =',numtiednode
        write(logunit,1)'row','node_id','elem_id','face_id','flag',
     &                   'adjust','x_dof1','y_dof2','z_dof3',
     &                   'x_original','y_original','z_original',
     &                   'delta'
        do row = 1, numtiednode
          node_id = tiednodeid(row,1)      ! slave node id
          dist_point = srt_dist3d(tiednodeglobalcoord(row,1:dimcoord),
     &                    nodecoord(node_id,1:dimcoord),dimcoord)

          write(logunit,2)row,tiednodeid(row,1:maxtieddatacol),
     &                     tiednodeadjustflag(row),
     &                     tiednodeglobalcoord(row,1:dimcoord),
     &                     nodecoord(node_id,1:dimcoord),dist_point
        end do
        write(logunit,*)
        write(logunit,*)'tied slave node updated global location'
        write(logunit,*)'femap netural file node format:'
        write(logunit,*)'(copy and paste to femap data block 403)'
        write(logunit,*)'(may need to offset the node id numbers)'
        do row = 1, numtiednode
          write(logunit,7)row,tiednodeglobalcoord(row,1:dimcoord)
        end do
        write(logunit,*)
      end if
      close(unit=logunit)
      return
c
 1    format(4a10,a5,a8,6a16,a16)
 2    format(4i10,i5,l8,3es16.6,3es16.6,es16.6)
 4    format(4i10,i5,3f16.8)
 5    format(4a10,a5,3a16)
 7    format(i7,'   0 0 1 46 0 0 0 0 0 0  ',3es16.8,' 0')
c
      end subroutine tied_nodeiso_y

      subroutine tied_nodeiso_z
c
      write(out,9002)
      do cnt = 1, gap_cnt
            write(out,9003) gap_nodes(cnt),gap_elems(cnt),
     &                      gap_faces(cnt),gap_dists(cnt),
     &                      gap_toler(cnt)
      end do
      deallocate(gap_nodes,gap_elems,gap_faces,gap_dists,gap_toler)
c
      return
 9002 format(//,10x,'>> warning: tied contact processing',/,
     &      10x,'>>   the following dependent slave nodes were not',/,
     &      10x,'>>   within the given gap tolerance of their',/,
     &      10x,'>>   assigned master elements',/,
     &      10x,'>>   the nodes will not be tied to their master',/,
     &      10x,'>>   surfaces; check the surface definitions',//,
     & 5x,'slave node',3x,'master elem',3x,'face',6x,'gap distance',
     & 7x,'gap tolerance',/,
     & 80('-'))
 9003 format(7x,i7,8x,i5,8x,i1,4x,es15.6,5x,es15.6)
c
      end subroutine tied_nodeiso_z


      subroutine tied_nodeiso_aa
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)'----- end tied_nodeiso -----'
      close(unit=logunit)
c
      return
c
      end subroutine tied_nodeiso_aa
c
      subroutine tied_nodeiso_bb
c
      allocate(shapefunc(maxnodeperelem),stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "shapefunc()" array'//
     &    ' (tied_nodeiso).'
        return
      end if

      allocate(shapederivs(maxnodeperelem,maxisodof),
     &         stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "shapederivs()" array'//
     &    ' (tied_nodeiso).'
        return
      end if

      allocate(elemcoord(maxnodeperelem,maxisodof),
     &         stat=allocate_error)
      if(allocate_error.gt.0)then
        mflag = 2
        message = 'internal error: memory error; array allocation'//
     &    ' error status > 0; check the local "elemcoord()" array'//
     &    ' (tied_nodeiso).'
        return
      end if

      return
c
      end subroutine tied_nodeiso_bb


      subroutine tied_nodeiso_cc
c
      if(abs(isopoint(1)).gt.iso_error_val)then
        xi_error_count = xi_error_count + 1
        select case(check_iso_coord_flag)
        case(1)      ! give a warning
          mflag = 1
          message = 'warning: the "xi" isoparametric coordinate'//
     &            ' is outside the element; expecting'//
     &            ' -1.0 <= xi <= 1.0; check the log file'//
     &            ' (tied_nodeiso).'
          if( logflag >= 1 ) call tied_nodeiso_q
        case(2)      ! set the maximum isoparametric value
          if( logflag >= 1 ) call tied_nodeiso_r
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
     &            ' is outside the element; expecting'//
     &            ' -1.0 <= eta <= 1.0; check the log file'//
     &            ' (tied_nodeiso).'
          if( logflag >= 1 ) call tied_nodeiso_s
        case(2)      ! set the maximum isoparametric value
          if( logflag >= 1 ) call tied_nodeiso_t
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
     &            ' is outside the element; expecting'//
     &            ' -1.0 <= zeta <= 1.0; check the log file'//
     &            ' (tied_nodeiso).'
          if( logflag >= 1 ) call tied_nodeiso_u
        case(2)      ! set the maximum isoparametric value
          if( logflag >= 1 ) call tied_nodeiso_v
          isopoint(3) = sign(one,isopoint(3))
          iso_adjust_flag = .true.
        end select
      end if
c
      return
c
      end subroutine tied_nodeiso_cc

      subroutine tied_nodeiso_dd
c
      do enode = 1, node_per_elem
         node_row = elemcon(elem_id,enode)
         if( node_row <= 0 .or. node_row > numnode) then
            mflag = 2
            message = 'error: invalid node id from the element'//
     &                ' connectivity; cannot get the node '//
     &                ' coordinates of the element (tied_nodeiso).'
            return
         end if
      end do

      return
c
      end subroutine tied_nodeiso_dd

      subroutine tied_nodeiso_ee
c
      integer :: ii, etype, idummy
      double precision :: enode_coords(3), ddummy
c        
c              is the slave node coincident on a node of master
c              element face?
c              if so, return iso coords of that node on the master
c              element.
c              the follow on search loop will converge on 1st
c              iteration.
c  
      tied_node_is_a_master_node = .false.
c
      do ii = 1, numnodeperface
       enode = localfacenode(face_id,ii)
       enode_coords(1:3) = elemcoord(enode,1:3) 
       if( srt_dist3d( slavenode, enode_coords, 3 ) <= 
     &     epsilon*character_size ) then
         tied_node_is_a_master_node = .true.
         etype = elemchoice  ! protect elemchoice value
         call ndpts1( idummy, numnodeperface, ddummy, etype, enode,                   
     &                isopoint(1), isopoint(2), isopoint(3) )
         return
       end if
      end do
c
      return
c
      end subroutine tied_nodeiso_ee

      subroutine tied_nodeiso_ff
c
      open(unit=logunit,file=logfile,position='append')
      write(logunit,*)
      write(logunit,*)  '                  slave node =',node_id
      if( tied_node_is_a_master_node ) then
       write(logunit,*) '  is a node on master element=',elem_id
      else
       write(logunit,*) '  is NOT node on master element=',elem_id
      end if
      write(logunit,*)
      close(unit=logunit)
c
      end subroutine tied_nodeiso_ff

      subroutine tied_nodeiso_gg
c
      integer :: enode, k, j
      double precision :: lengths(1),
     &                    enode_coords(mxvl,3*node_per_elem)

      j = node_per_elem + 1
      k = j + node_per_elem 
c     
      do enode = 1, node_per_elem 
        enode_coords(1,enode) = elemcoord(enode,1) 
        enode_coords(1,j) = elemcoord(enode,2) 
        enode_coords(1,k) = elemcoord(enode,3) 
        j = j + 1
        k = k + 1
      end do
c
      call characteristic_elem_length( elemchoice, 1, node_per_elem,
     &                                 enode_coords, lengths,  iout  )
      character_size = lengths(1)
c
      return
c
      end subroutine tied_nodeiso_gg

      subroutine tied_nodeiso_hh
c
      isopoint_old = zero
c
      loop_1: do itercount = 1, maxiterloop
        if( debug_now ) call tied_nodeiso_g
c
c             compute the face isoparametric point in the global
c             x,y,z coordinates; use the element shape functions
c             evaluated at the isoparametric point
c
        call shapef(elemchoice,isopoint(1),isopoint(2),isopoint(3),
     &              shapefunc)

        call tied_globalcoordfromiso(globalpoint,maxdof,
     &              shapefunc,maxnodeperelem,node_per_elem,
     &              elemcoord,maxnodeperelem,maxisodof,
     &              logflag,logunit,logfile,mflag,message)
        if( mflag >=2 ) return
c
c             get the position vector from the face isoparametric
c             point to the slave node; use the global x,y,z position
c             of the face point to compute the position vector
c
        call srt_getunitvec(globalpoint,slavenode,position,maxdof,0)
c
c             get the unit normal and unit tangent vectors to the
c             master element face at the face isoparametric point;
c             evaluate the element shape functions at the current
c             isoparametric coordinates and compute the jacobian
c             matrix
c
        select case(elemchoice)
        case(1,2)      ! 20 or 8 node brick element
c
c             evaluate the element shape functions at the current
c             isoparametric coordinates and compute the jacobian
c             matrix
c
          call derivs(elemchoice,isopoint(1),isopoint(2),isopoint(3),
     &                shapederivs(:,1),shapederivs(:,2),
     &                shapederivs(:,3))
          ierr = 0
          do i = 1, node_per_elem
              dsf(i,1:3) = shapederivs(i,1:3)
              coord(1:3,i) = elemcoord(i,1:3)
          end do
          call eqldjb(dsf,coord,node_per_elem,jacobian,invj,detj,ierr)
          if( ierr > 0 ) then
            mflag = 2
            message = 'internal error: error while computing the'//
     &                ' jacobian matrix; ierr > 0 (tied_nodeiso).'
            return
          end if
c
c             use the row of the jacobian() matrix that matches the
c             constant isoparametric dof for the face normal vector;
c             the other rows from the jacobian() matrix give the
c             master element face tangent vectors
c
          normal(1:maxdof) = jacobian(constantisodof,1:maxdof)
c
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
     &              ' value; cannot set the element face tangent '//
     &              ' vectors (tied_nodeiso).'
            return
          end select
c
        case(6,13)      ! 10 or 4 node tetrahedron
c
c             use the triangular face corners to get orthogonal
c             tangent vectors on the tetrahedron element face
c
          call tied_tetfacetangent(tangent1,tangent2,normal,maxdof,
     &                elemcoord,maxnodeperelem,maxisodof,
     &                node_per_elem,face_id,
     &                localfacenode,maxelemface,maxnodeperface,
     &                numelemface,numnodeperface,numcornerperface,
     &                mflag,message,logflag,logunit,logfile)
          if( mflag >=2 ) return
c
        case default
          mflag = 2
          message = 'error: unexpected elemchoice; cannot set the'//
     &            ' initial isoparametric point location'//
     &            ' (tied_nodeiso).'
          return
        end select
c
c             normalize to get unit vectors
c
        call srt_normalize(normal,maxdof)
        call srt_normalize(tangent1,maxdof)
        call srt_normalize(tangent2,maxdof)
c
        if( debug_now ) call tied_nodeiso_h
c
c             check the distance from the face point to the slave
c             node; if close enough then the face point is located
c             at the slave node on the master element face
c
        dist_point = srt_dist3d(globalpoint,slavenode,maxdof)
        if( dist_point <= epsilon*character_size) then
          if( debug_now ) call tied_nodeiso_i
          exit loop_1
        end if
c
c             check the position change of the element face point
c             from the previous iteration; if the change of point
c             position is small enough then the face point is
c             converged to the mapped slave node location
c
        if( itercount >= 2 ) then
          delta_iso = srt_dist3d(isopoint,isopoint_old,maxdof)
          if( delta_iso <= epsilon*character_size) then
            if( debug_now ) call tied_nodeiso_j
            exit loop_1
          end if
        end if
c
c             check the angle from the face point to the slave node;
c             the slave node may be away from the master element face
c             but if the angle between the face normal and the
c             position is small then the slave node is projected onto
c             the master element face at the face point location
c
        angle_point = srt_getangle(position,normal,maxdof)
        if( angle_point <= epsilon_angle .or.
     &   dabs(angle_point-pi) <= epsilon_angle) then
          if( debug_now ) call tied_nodeiso_k
          exit loop_1
        end if
c
c             set the scaling factors used to update the face point
c             location
c
        if( itercount == 1 ) then
          scale1 = zero
          scale2 = zero

          select case(elemchoice)
          case(1,2)      ! 20 or 8 node brick element
c
c             compute the average projected distance of the vectors
c             from the master element face center (face point in the
c             first iteration) to the corner nodes onto the tangent
c             vectors; use the scale1,scale2 values to scale the
c             delta1,delta2 values (this is an approximate mapping of
c             the element size back to the isoparametric element size)
c
            do corner = 1, numcornerperface
              node_row = localfacenode(face_id,corner)
              call srt_getunitvec(globalpoint,
     &                elemcoord(node_row,1:maxdof),vec,maxdof,0)
              scale1 = scale1 +
     &                 abs(srt_dotprod(vec,tangent1,maxdof))
              scale2 = scale2 +
     &                 abs(srt_dotprod(vec,tangent2,maxdof))
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
c
          if( debug_now )  call tied_nodeiso_l 
          if( dabs(scale1) < 1.0d-8 )then
            mflag = 2
            message = 'internal error: scale1 = 0.0;'//
     &              ' cannot scale the delta1 value '//
     &              ' (tied_nodeiso).'
          else if( dabs(scale2) < 1.0d-8 ) then
            mflag = 2
            message = 'internal error: scale2 = 0.0;'//
     &              ' cannot scale the delta2 value '//
     &              ' (tied_nodeiso).'
          end if
          if( mflag >=2 ) return
        end if
c
c             project the position vector onto the face tangent
c             vectors; use the local element face coordinates to
c             update the isoparametric point position
c
        delta1 = srt_dotprod(position,tangent1,maxdof)/scale1
        delta2 = srt_dotprod(position,tangent2,maxdof)/scale2
        isopoint_old(1:maxdof) = isopoint(1:maxdof)

        if( debug_now ) call tied_nodeiso_m
c
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
     &              ' value; cannot update the element face point'//
     &              ' location (tied_nodeiso).'
            return
          end select

        case(6,13)      ! 10 or 4 node tetrahedron
          do i = 1 , maxdof
            vec1(i) = delta1*tangent1(i)
            vec2(i) = delta2*tangent2(i)
          end do
          do i = 1, maxdof
            globalpoint(i) = globalpoint(i) + vec1(i) + vec2(i)
          end do

          if( debug_now ) call tied_nodeiso_n

          call tied_tetfaceisopoint(isotetpoint,maxtetdof,
     &                globalpoint,maxdof,constantisodof,
     &                elemcoord,maxnodeperelem,maxisodof,
     &                node_per_elem,face_id,
     &                localfacenode,maxelemface,maxnodeperface,
     &                numelemface,numnodeperface,numcornerperface,
     &                mflag,message,logflag,logunit,logfile)
          if( mflag >=2 ) return

          do i = 1, maxtetdof
            if( isotetpoint(i) < -0.01d0 )then
              mflag = 2
              message = 'internal error: tetrahedron isoparametric'//
     &              ' coordinate < 0.0; check the log file'//
     &              ' (tied_nodeiso).'
            else if(isotetpoint(i).gt.1.01d0)then
              mflag = 2
              message = 'internal error: tetrahedron isoparametric'//
     &              ' coordinate > 1.01; check the log file'//
     &              ' (tied_nodeiso).'
            end if
            if( mflag >= 2 ) then
              if( logflag >= 1 ) call tied_nodeiso_o
              return
            end if
          end do ! on i
c
c             save s2,s3,s4 as the independent isoparametric points on
c             the tetrahedron face, use isopoint() when evaluating the
c             shape functions
c
          do i = 1, maxdof
            isopoint(i) = isotetpoint(i+1)
          end do

        case default
          mflag = 2
          message = 'error: unexpected elemchoice; cannot update'//
     &            ' the isoparametric point location'//
     &            ' (tied_nodeiso).'
          return
        end select
c
        if( debug_now )  call tied_nodeiso_p
c
      end do loop_1   !      *** end loop_1  ****
c
      return

      end subroutine tied_nodeiso_hh



      end subroutine tied_nodeiso


