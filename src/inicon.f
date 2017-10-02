c     ****************************************************************
c     *                                                              *
c     *                      subroutine inicon                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 09/21/2017 rhd             *
c     *                                                              *
c     *     supervises and conducts the input of the                 *
c     *     desired initial conditions for the structure at time 0.  *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine inicon( sbflg1, sbflg2 )
      use global_data ! old common.main
c
      use main_data, only : trn, trnmat, temper_nodes,
     &                      temper_nodes_ref, temperatures_ref,
     &                      inverse_incidences, initial_stresses,
     &                      initial_stresses_user_routine,
     &                      initial_stresses_file
      implicit none
c
      logical :: sbflg1, sbflg2
c
c                       local declarations
c
      integer :: i, j, cond, lenlst, errnum, param, dofn, icn,
     &           iplist, node, elem, type, nedof, idum, nc, scan_stat
      integer, allocatable :: intlst(:)
      real :: dumr
      double precision :: cval, mpfact, dumd
      double precision, allocatable :: edva(:,:), trnmte(:,:,:)
      double precision, parameter :: zero = 0.0d0
      character name*80, iclnam*8, dums, curtyp*1
      logical :: found, dvaflg(mxndof), dump
      logical, allocatable :: trne(:,:)
      logical, external :: matchs, endcrd, true, label, numd, scanms,
     &                     string
c
c                       if sub flag 1 is on, there is reentry into
c                       inicon after an error in the input of type
c                       of initial condition input.
c
      allocate( edva(mxvl,mxndof), trnmte(mxvl,mxedof,mxndof),
     &          trne(mxvl,mxndel), intlst(mxlsz) )
c
      if( sbflg1 ) then
         call errmsg(119,idum,dums,dumr,dumd)
      end if
c
      call readsc  ! line after line containing keyword: initial
c
c                       branch on the type of initial condition to
c                       be input. if key words not encountered,
c                       return to driver subroutine and look for a
c                       high level command.
c
      do
        if( matchs('nodal',4) ) call splunj ! skip optional word
        if( matchs('displacement',7) ) then
           cond = 1
           call inicon_node_values ! returns on intlst not found
        elseif( matchs('velocity',3) ) then
           cond = 2
           call inicon_node_values
        elseif( matchs('temperature',4)  ) then
           cond = 4
           call inicon_node_values
        elseif( matchs('stresses',4)  ) then
           cond = 5
           call inicon_initial_stresses
           if( scan_stat .eq. 1 ) call readsc
        else
           sbflg1 = .true.
           sbflg2 = .false.
           return
        end if
      end do

      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *     contains:   inicon_node_values                           *
c     *                                                              *
c     *                 last modified : 09/21/2017 rhd               *
c     *                                                              *
c     ****************************************************************


      subroutine inicon_node_values
      implicit none
c
c                       read the list of nodes whose initial conditions
c                       are to be input.
c
 1120 call readsc
      if( matchs('nodes',4) ) call scan
      call trlist(intlst,mxlsz,nonode,lenlst,errnum)
c
c                       branch on the return code from trlist. a
c                       value of 1 indicates no error. a value of
c                       2 indicates that the parse rules failed in
c                       the list. a value of 3 indicates that the
c                       list overflowed its maximum length of mxlsz.
c                       in these last two cases, the illegal list
c                       will be ignored and a new node list will
c                       be sought. a value of 4 indicates that no list
c                       was found. in this case, input of initial cond-
c                       itions has ceased.
c
      if( errnum .eq. 2 ) then
         param = 1
         call errmsg( 24, param, dums, dumr, dumd )
         go to 1120
      else if( errnum .eq. 3 ) then
         param = 2
         call errmsg( 24, param, dums, dumr, dumd )
         go to 1120
      else if( errnum .eq. 4 ) then ! no list found
         return
      else
         if( errnum .eq. 1 ) then
            call backsp( 1 )
            if( true(idum) ) go to 1125
         end if
         param = 3
         call errmsg( 24, param, dums, dumr, dumd )
         go to 1120
      end if
c
c                       there is a valid nodal list. read the dof
c                       to be assigned initial conditions.
c
c                       initialize the d/v/a/ and dof type arrays.
c                       use dofn = 0 for temperature
c
 1125 do i = 1, mxndof
         edva(1,i) = zero
         dvaflg(i) = .false.
      end do
c
c                       input the dof type. if found, branch to
c                       store its initial value.
c
 1127 if( matchs('u',1) ) then
         dofn   = 1
         curtyp = 'U'
         go to 1129
      end if
c
      if( matchs('v',1) ) then
         dofn   = 2
         curtyp = 'V'
         go to 1129
      end if
c
      if( matchs('w',1) ) then
         dofn   = 3
         curtyp = 'W'
         go to 1129
      end if
c
      if( matchs('temperature',4) ) then
         dofn   = 0
         curtyp = 'T'
         go to 1129
      end if
c
      if( matchs(',',1) ) go to 1128
c
c                       if there is an end of card, in. cond. input has
c                       ended. branch to store the temp. vec. globally.
c                       if not, ignore the current entity and search
c                       for another dof to input.
c
      if( endcrd(idum) ) then
         go to 1135
      else
         call errmsg( 72, idum, dums, dumr, dumd )
         if( true(idum) ) go to 1127
      end if
c
c                       if there is a comma at the end of a line, the
c                       input line is continued.
c
 1128 continue
      if( endcrd(idum) ) then
         call readsc
      end if
      go to 1127
c
c                       store the initial conditions.
c
 1129 if( matchs('=',1) ) call splunj
      if( .not. numd(cval) ) then
         call errmsg( 73, idum, dums, dumr, dumd)
      else
c
c                       make sure that the current dof has not
c                       already been input on the same nodal initial
c                       conditions command.
c
         if( dofn .eq. 0 ) go to 1127
         if( dvaflg(dofn) ) then
            call errmsg( 117, idum, curtyp, dumr, dumd )
            go to 1127
         end if
c
         dvaflg(dofn) = .true.
         edva(1,dofn) = cval
c
      end if
c
      go to 1127
c
c                       store initial conditions globally.
c
 1135 icn    = 0
      iplist = 1
 1136 call trxlst( intlst, lenlst, iplist, icn, node )
c
c                       check that the list node does not exceed
c                       the number of nodes in the structure.
c
         if( node .gt. nonode ) then
            param = node
            call errmsg( 16, param, dums, dumr, dumd )
            go to 1140
         end if
c
c                       check that the list node is not negative.
c
         if( node .lt. 0 ) then
            param = node
            call errmsg( 58, param, dums, dumr, dumd )
            go to 1140
         end if
c
         if ( dofn .eq. 0 ) then
            temper_nodes(node) = cval
            temper_nodes_ref(node) = cval
            if ( abs(cval) .ne. zero ) temperatures_ref = .true.
            go to 1140
         end if
c
         elem = inverse_incidences(node)%element_list(1)
         type = iprops(1,elem)
         nedof = iprops(4,elem)  !  should always = 3
c
c                       as the displacements, velocities, and
c                       accelerations are input in uniform global
c                       coordinates, they must be transformed to
c                       constraint compatable global coordinates
c                       before global storage.
c
c                       extract transformation matrix for this node
c                       and transform to ccg coordinates.
c
         trne(1,1) = trn(node)
         if( trne(1,1) ) then
            trnmte(1,1:nedof,1:nedof) =
     &              trnmat(node)%mat(1:nedof,1:nedof)
            call trnvec( edva, trnmte, trne, nedof, 1, 1, 1 )
         end if
c
c                       store the transformed temporary vector
c                       globally, which contains the dis/vel/acc com-
c                       patable with the current node in consecutive
c                       order.
c
         do i = 1, nedof
            if( cond .eq. 1 ) then
               u(dstmap(node)+i-1) = edva(1,i)
            else if( cond .eq. 2 ) then
               v(dstmap(node)+i-1) = edva(1,i)
            else
               a(dstmap(node)+i-1) = edva(1,i)
            end if
         end do
c
c
 1140 if( iplist .ne. 0 ) go to 1136
c
c                       return to process more lists.
c
      go to 1120 ! next list of nodes on new line
c
c
c
      return
      end subroutine inicon_node_values
c     ****************************************************************
c     *                                                              *
c     *     contains:   inicon_initial_stresses                      *
c     *                                                              *
c     *                 last modified : 09/22/2017 rhd               *
c     *                                                              *
c     ****************************************************************

      subroutine inicon_initial_stresses
      implicit none

      logical :: file_exists, found
      character(len=100) :: workstr
      integer :: result

      initial_stresses_user_routine = .false.
      initial_stresses_file = " "
      scan_stat = 0
c
      if( .not. allocated( initial_stresses ) )
     &    allocate( initial_stresses(6,noelem) )
      initial_stresses = zero
c
c              user routine with optional file name
c
      if( matchs('user_routine',8) ) then
        initial_stresses_user_routine = .true.
        if( matchs('file',4) ) call splunj
        workstr(1:100) = " "
        if( label(idum) ) then
           call entits( workstr, nc )
        elseif( string(idum) ) then
           call entits( workstr, nc )
        end if
        initial_stresses_file = workstr
        inquire( file=workstr, exist=file_exists )
        if( .not. file_exists ) write(out,9000)
     &                        workstr(1:len_trim(workstr))
        scan_stat = 1
        return
      end if
c
c              regular lists of elements and stress values on each
c              logical input lines.
c
      scan_stat = 0 ! no read on return
      do
c
        call readsc
        if( matchs('dump',4) ) then
           call inicon_dump_stresses
           cycle
        end if
        call inicon_get_list( result )
        if( result .eq. 4 ) then ! no list. do we have non-zero stresses
           call inicon_chk_stress_tbl( found )
           if( .not. found ) then
            deallocate( initial_stresses )
            write(out,9010)
           end if
           return ! no list found
        end if
        if( result .ne. 1 ) then  ! bad list syntax
           initial_stresses = zero
           cycle
        end if
        call inicon_get_store_stresses
c
      end do

      return
c
 9000 format(/1x,'>>>>> warning: file does not exist: ', a,/)
 9010 format(/1x,'>>>>> warning: no valid initial stress data',
     &  ' found',/)
c
      end subroutine inicon_initial_stresses
c
c     ****************************************************************
c     *                                                              *
c     *     contains:   inicon_chk_stress_tbl                        *
c     *                                                              *
c     *                 last modified : 09/22/2017 rhd               *
c     *                                                              *
c     ****************************************************************

      subroutine inicon_chk_stress_tbl( found )
      implicit none
c
      integer :: i, j
      logical :: found
      double precision :: sum

      do i = 1, noelem
        do j = 1, 6
           sum = sum + abs( initial_stresses(j,i) )
        end do
      end do
c
      found = sum > zero
c
      return
c
      end subroutine inicon_chk_stress_tbl
c     ****************************************************************
c     *                                                              *
c     *     contains:   inicon_get_store_stresses                    *
c     *                                                              *
c     *                 last modified : 09/22/2017 rhd               *
c     *                                                              *
c     ****************************************************************

      subroutine inicon_get_store_stresses
      implicit none
c
      integer :: i, nc, icn, iplist, elem
      double precision :: local_vals(6)
      character(len=80) bad_data
c
      local_vals = zero
c
c              get up to 6 values from remainder of line
c              ordering: xx, yy, zz, xy, yz, xz
c
      do i = 1, 6
       if( endcrd() ) exit
       if( numd( local_vals(i) ) ) cycle
       call entits( bad_data, nc )
       write(out,9000) bad_data(1:nc)
       num_fatal = num_fatal + 1
       initial_stresses = zero
       return
      end do
c
c              save all these values to each element in the list
c
      icn = 0; iplist = 1
c
      do
       if( iplist .eq. 0 ) exit
       call trxlst( intlst, lenlst, iplist, icn, elem )
       if( elem > 0 .and. elem <= noelem ) then
         initial_stresses(1:6,elem) = local_vals
       else
         write(out,9010) elem
         num_fatal = num_fatal + 1
         initial_stresses = zero
       end if
      end do


      return
c
 9000 format(/1x,'>>>>> error: expecting a stress value. ',
     &   /16x,'found: ',a,/ )
 9010 format(/1x,'>>>>> error: invalid element number: ',i10,/ )
c
      end subroutine inicon_get_store_stresses
c     ****************************************************************
c     *                                                              *
c     *     contains:   inicon_get_list                              *
c     *                                                              *
c     *                 last modified : 09/22/2017 rhd               *
c     *                                                              *
c     ****************************************************************
c
      subroutine inicon_get_list( result )
      implicit none
c
      integer :: result
c
      if( matchs('elements',4) ) call scan
      call trlist( intlst, mxlsz, noelem, lenlst, result )
c
      if( result .eq. 1 ) then ! list found set scan for stresses list
          call backsp( 1 )
          if( true(idum) ) call splunj
          return
      end if
c
      if( result .eq. 4 ) return ! no list found @ start of line
c
      if( result .eq. 2 ) then ! parse rules failed
         param = 1
         call errmsg( 24, param, dums, dumr, dumd )
         return
      else if( result .eq. 3 ) then ! list overflow
         param = 2
         call errmsg( 24, param, dums, dumr, dumd )
         return
      else
         param = 3
         call errmsg( 24, param, dums, dumr, dumd )
         return
      end if
c
      return
c
      end subroutine inicon_get_list
c     ****************************************************************
c     *                                                              *
c     *     contains:   inicon_dump_stresses                         *
c     *                                                              *
c     *                 last modified : 09/22/2017 rhd               *
c     *                                                              *
c     ****************************************************************
c
      subroutine inicon_dump_stresses
      implicit none
c
      integer :: i, j
      double precision :: sum
c
      write(out,9000)
      do i = 1, noelem
        sum = zero
        do j = 1, 6
         sum = sum + abs(initial_stresses(j,i))
        end do
        if( sum > zero ) write(out,9010) i, initial_stresses(1:6,i)
      end do
c
      return
c
 9000 format(5x,"..... initial stresses. only elements with non-zero",
     &   " values listed ....." )
 9010 format(2x,i8,6d16.6)
c
      end subroutine inicon_dump_stresses
c
      end subroutine inicon

