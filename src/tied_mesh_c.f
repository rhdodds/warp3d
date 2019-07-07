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
