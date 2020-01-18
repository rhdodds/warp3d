c     ****************************************************************
c     *                                                              *
c     *                subroutine release_constraints                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 12/3/2019 rhd               *
c     *                                                              *
c     *     read/store data for the release constraints command      *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine release_constraints( sbflg1, sbflg2 )
      use global_data ! old common.main
c
      use main_data, only : cnstrn_in, release_cons_table,
     &                      release_cons_steps, mdiag, rload
      use mod_mpc, only:    num_user_mpc, user_mpc_table
      use allocated_integer_list
c
      implicit none
c
c              parameters
c
      logical :: sbflg1, sbflg2
c
c              locals
c
      integer :: list_size, list_size_mpc, lenlst, 
     &           mpc_lenlst, warn_mess_num_steps
      integer, save :: count_stored
      integer, allocatable :: user_mpc_col_list(:), intlst(:),
     &                        mpc_intlst(:)
      logical ::  found_list, bad_list, release_flags(3), debug,
     &            compact_mpc_table, warn
      logical, save :: delete_release_cons
c
c              if sub flag 1 is on, there is re-entry after an error
c              in input. if we have bad input somewhere here,
c              release the data structure on exit
c
      debug = .false.
      if( debug ) then
         write(out,*) '... entered release_constraints ...'
         write(out,*) '    sbflg1, sbflg2: ', sbflg1, sbflg2
      end if
c 
      allocate( intlst(10) )
      count_stored      = 0
      compact_mpc_table = .false.
      warn_mess_num_steps = 0
c
      if( num_user_mpc > 0 ) then
        allocate( user_mpc_col_list(max_mpc) )
      else
        allocate( user_mpc_col_list(1) )
      end if
      user_mpc_col_list = 0
c
      if( sbflg1 ) then
         write(out,9000) ; num_error = num_error + 1
      else
         delete_release_cons = .false.
      end if
c
c              get number of steps used to release reaction
c              forces to zero. default = 1 step
c
      if( .not. sbflg1 ) call release_cons_init 
c
      do !   outer read loop over lines of node lists and dofs
c
        call readsc
        call release_cons_scan
        if( .not. found_list ) exit  ! out of read loop. back to main
        if( bad_list ) then
          delete_release_cons = .true.
          cycle
        end if
c
c              input line looks ok. update data structures
c              for the releases. only 1 release step
c              allowed for MPCs
c
        if( mpc_lenlst == 0 ) then  ! cmd had only absolute constraints
          call release_cons_abs
        else
          if( num_user_mpc .eq. 0 ) then  ! cmd had MPC but none defined
             write(out,9100)
             delete_release_cons = .true.
             num_error = num_error + 1
          else
             warn = warn_mess_num_steps .eq. 0 .and. 
     &              release_cons_steps .ne. 1
             if( warn ) then
               warn_mess_num_steps = warn_mess_num_steps + 1
               write(out,9110)
             end if
             call release_cons_mpc
          end if
        end if
c
      end do  ! get next input line
c
c              all release commands read/processed. cleanup.
c
      sbflg1         = .true.
      sbflg2         = .true.
c
      if( debug ) then
         write(out,*) '... count_stored: ',count_stored
         write(out,*) '... delete_release_cons: ', delete_release_cons
         write(out,*) '... compact_mpc_table: ', compact_mpc_table
      end if
c
      if( delete_release_cons ) then
        if( allocated( release_cons_table ) )
     &      deallocate( release_cons_table )
      end if
c
      call release_empty_cons_table
c
      if( compact_mpc_table ) then
        if( debug ) then
           write(out,*) '.. ready to compact. num_user_mpc: ',
     &                  num_user_mpc
           write(out,*) '.. user_mpc_col_list: ',
     &                  user_mpc_col_list(1:num_user_mpc)
         end if
         call incon_mpcs_resize( 2, user_mpc_col_list )
      end if
      if( debug ) write(out,*) '... leaving release_constraints ...'
c
      return
c
 9000 format(/1x,'>>>>> error: syntax error in release constraints',
     &   /,14x,'input. read new line ...',/)
 9100 format(/1x,'>>>>> error: no user MPCs defined for model.',
     &   ' read new line ...',/)
 9110 format(/1x,'>>>>> Warning: only 1 release step allowed for ',
     &   'MPCs.',
     &  /,14x,'Using 1 step',/)
c
      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *              subroutine release_cons_init                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/20/2018 rhd             *
c     *                                                              *
c     *           initialize the release cons table                  *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine release_cons_init
      implicit none
c
      integer :: i, dummy
      logical, external :: matchs, endcrd, integr
      double precision, parameter :: zero = 0.d0
c
      release_cons_steps = 1
      if( matchs( 'constraints', 4 ) ) call splunj
      if( matchs( 'steps',4 ) ) call splunj
      if( matchs( '=',1 ) ) call splunj
      if( .not. integr( release_cons_steps ) ) then
         if( .not. endcrd(dummy) ) then
           write(out,9120)
           num_error = num_error + 1
         end if
      end if
c
c                       allocate release constraints table if needed
c
      if( .not. allocated( release_cons_table ) ) then
        if( debug ) write(out,*) '... allocating release_cons_table'
        allocate( release_cons_table(3,nonode) )
        do i = 1, nonode
           release_cons_table(1:3,i)%num_release_steps = 0
           release_cons_table(1:3,i)%remaining_steps_for_release = 0
           release_cons_table(1:3,i)%reaction_force = zero
        end do
      end if
c
      return
c
 9120 format(/1x,'>>>>> error: invalid number of release steps',
     &   /,14x,'input. read new line ...')
c
      end subroutine release_cons_init

c     ****************************************************************
c     *                                                              *
c     *              subroutine release_cons_scan                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/3/2019 rhd              *
c     *                                                              *
c     *     scan/store list of nodes and constraint directions to    *
c     *     be released                                              *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine release_cons_scan
c
c              locals
c
      integer :: strlng, errnum, param, dummy
      logical :: ok, done_input
      logical, external :: matchs, match_exact, endcrd, integr, true
      character :: dums*1, string*80
      real :: dumr
      double precision :: dumd
c
      if( debug ) write(out,*) ' '
      if( debug ) write(out,*) '    ... entered release_cons_scan ...'
c
c              command format:
c
c                <node list> [ dof ]   absolute cons
c                <node list a> [ dof] <node list b>        MPCs
c
c              for MPC, the 1st entry in list a must match 1st entry
c              in list b in the user-defined MPC equations for model.
c              same for 2nd, 3rd, 4th, .. entries in each list.
c              each such node pair order does not have to match the
c              user-defined MPC equation. For MPC
c
c                - 4823 1.0 u + 9954 1.0 u  = 0.0
c
c              these forms are ok for in the release input
c                  4823 u 9954
c                  9954 u 4823
c
c              the processors here *assume* the MPC is always of the
c              type   node_a/dof = node_b/dof  
c
      if( matchs( 'plane', 5 ) ) call splunj ! implement later
c
      mpc_lenlst = 0
c
      call trlist_allocated( intlst, list_size, nonode,
     &                       lenlst, errnum )
c
c              branch on the return code from trlist. a
c              value of 1 indicates no error. a value of
c              2 indicates that the parse rules failed in
c              the list. a value of 3 indicates that the
c              list overflowed its maximum length of mxlsz.
c              in these last two cases, the illegal list
c              will be ignored and a new node list will
c              be sought. a value of 4 indicates that no list
c              was found. in this case, release cons input
c              has ceased.
c
      if( debug ) write(out,*)
     & '    ... trlist done, errnum, lenlst: ', errnum, lenlst
c
      bad_list   = .false.
      found_list = .false.
c
      if( errnum .eq. 4 ) return
      if( errnum .eq. 2 ) then
         param = 1
         call errmsg(24,param,dums,dumr,dumd)
         bad_list = .true.
         call scan_flushline
         return
      end if
c
      if( errnum .eq. 3 ) then
         param = 2
         call errmsg(24,param,dums,dumr,dumd)
         bad_list = .true.
         call scan_flushline
         return
      end if
c
      if( errnum .eq. 1 ) then
          found_list = .true.
          call backsp(1)
          if( true(dummy) ) call splunj
      else
         write(out,9000)
         call die_abort
      end if
c
c                       get one or more directions to be released
c
      release_flags(1:3) = .false.
      do
         if( match_exact( ',',1 ) ) then
             call splunj ! do nothing forces compiler to execute
             cycle
         end if
         if( match_exact( 'u',1 ) ) then
             release_flags(1) = .true.
             cycle
         end if
         if( match_exact( 'v',1 ) ) then
             release_flags(2) = .true.
             cycle
         end if
         if( match_exact( 'w',1 ) ) then
             release_flags(3) = .true.
             cycle
         end if
         if( endcrd(dummy) ) then
             done_input = .true.
             exit
         end if
         call trlist_allocated( mpc_intlst, list_size_mpc, nonode,
     &                          mpc_lenlst, errnum )
         if( errnum .eq. 4 ) then ! no list found but not endcrd
           bad_list = .true.
           num_error = num_error + 1
           call entits( string, strlng )
           write(out,9010) string(1:strlng)
           call scan_flushline
           return
         else
           exit  ! list found but could have list syntax errors
           done_input = .false.
        end if
      end do
c
      ok = release_flags(1) .or. release_flags(2) .or. release_flags(3)
      if( .not. ok ) then
          write(out,9020)
          num_error = num_error + 1
          call scan_flushline
          return
      end if
c
c              done or is this an MPC release ?
c
      if( debug ) write(out,9030) release_flags, done_input
      if( done_input ) return
c
c              we found an integerlist following the list of u, v, w
c
      bad_list = .false.
      if( debug ) write(out,*)
     & '    ... trlist found for mpc nodes, errnum, mpc_lenlst: ',
     &      errnum,lenlst
      if( errnum .eq. 2 ) then
         param = 1
         call errmsg(24,param,dums,dumr,dumd)
         bad_list = .true.
         call scan_flushline
         return
      end if
c
      if( errnum .eq. 3 ) then
         param = 2
         call errmsg(24,param,dums,dumr,dumd)
         bad_list = .true.
         call scan_flushline
         return
      end if
c
      if( errnum .eq. 1 ) then
         call backsp(1)
         if( true(dummy) ) call splunj
      else
         write(out,9000)
         call die_abort
      end if
c
      return
c
 9000 format(/1x,'>>>>> error: invalid return on trlist in ',
     &  /14x,'release_cons_scan. system error. job aborted.',//)
 9010 format(/1x,'>>>>> error: unrecognized data in list of released',
     &  /14x,'dof. scanning: ', a )
 9020 format(/1x,'>>>>> error: no valid components to release found')
 9030 format(5x,"... release flags, done_input:",4l3)
c
      end subroutine release_cons_scan
c     ****************************************************************
c     *                                                              *
c     *                 subroutine release_cons_abs                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/20/2018 rhd             *
c     *                                                              *
c     *     process list of absolute constraints to release          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine release_cons_abs
      implicit none
c
      integer :: icn, iplist, node, param, idof, nsteps, sdof
      character :: dums*1, dof_names(3)*1
      real :: dumr
      double precision :: dumd, react
      double precision, parameter :: d32460 = 32460.0d00 
      logical :: bad
      data dof_names / 'u', 'v', 'w' /
c
c              process list of nodes in release command for 
c              absolute constraints. list entries must be valid
c              node number, not already in release and have
c              requested release dof  defined in global constraint table.
c
      icn    = 0
      iplist = 1
c
      do while ( iplist .ne. 0 )
c
         call trxlst( intlst, lenlst, iplist, icn, node )
         if( node .gt. nonode ) then
            param = node
            call errmsg(16,param,dums,dumr,dumd)
            cycle
         end if
         if( node .le. 0 )  then
            param = node
            call errmsg(58,param,dums,dumr,dumd)
            cycle
         end if
c
         bad = .false.
         do idof = 1, 3
c
           if( .not. release_flags(idof) ) cycle
           nsteps = release_cons_table(idof,node)%num_release_steps
           if( nsteps .ne. 0 ) then
              write(out,9100) node, dof_names(idof)
              bad = .true.
           end if
           sdof = dstmap(node) + idof - 1
           if( cnstrn_in(sdof) .eq. d32460 ) then ! no constraint exists
              write(out,9110) node,  dof_names(idof)
              bad = .true.
           end if
c
         end do ! over idof = 1,2,3
c
         if( bad ) then
           num_error = num_error + 1
           cycle
         end if
c
c              valid node and currently constrained dofs. pull reaction
c              values from global load vector and store for gradual
c              reduction to zero of specified release steps.
c              remove absolute constraint from constraints data
c              structure.
c
         do idof = 1, 3
c
           if( .not. release_flags(idof) ) cycle
           sdof = 3 * (node-1) + idof
           release_cons_table(idof,node)%num_release_steps =
     &            release_cons_steps
           release_cons_table(idof,node)%remaining_steps_for_release =
     &            release_cons_steps
           react = -rload(sdof) + ifv(sdof) + mdiag(sdof)*a(sdof)
           release_cons_table(idof,node)%reaction_force = react
           call release_cons_update_constraints( sdof )
           count_stored = count_stored + 1
c
         end do ! over idof  = 1,2,3
c
      end do  ! over loop to process nodes in list
c
      return
c
 9100 format(/1x,'>>>>> error: node: ',i0,' dof: ',a1,
     &  /14x,'already in release')
 9110 format(/1x,'>>>>> error: node: ',i0,' dof: ',a1,
     &  /14x,'has no constraint to release')
c
      end subroutine release_cons_abs

c     ****************************************************************
c     *                                                              *
c     *                 subroutine release_cons_mpc                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/23/2018 rhd             *
c     *                                                              *
c     *     process list of mpc constraints to release               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine release_cons_mpc
c
      use stiffness_data, only : total_lagrange_forces, 
     &                           d_lagrange_forces, i_lagrange_forces                        
c
      implicit none
c
      integer, external :: iszlst
      integer :: nodes_lst_a, nodes_lst_b, icn, iplist, node, count,
     &           node_a, node_b, i, j, idof, param, 
     &           nsteps_a, nsteps_b, node_1, node_2, mpc_dof_1,
     &           mpc_dof_2, minab, maxab, min12, max12, sdof
      integer, allocatable :: list_a(:), list_b(:)
      double precision, parameter :: d32460 = 32460.0d00, zero = 0.d0 
      double precision :: dumd, react
      logical :: badlist, bad, t1, t2, found_mpc_entry, ok,
     &           local_debug1 
      real :: dumr, mpc_constant
      character :: dof_names(3)*1, dums*1
      data dof_names / 'u', 'v', 'w' /
c
c              the two node lists must have the same number of terms.
c
      local_debug1 = .false.
      if( debug ) write(out,9000)
c
      nodes_lst_a = iszlst( intlst, lenlst )
      nodes_lst_b = iszlst( mpc_intlst, mpc_lenlst )
      if( nodes_lst_a .ne. nodes_lst_b ) then
           write(out,9130) nodes_lst_a, nodes_lst_b
           num_error = num_error + 1
           return
      end if
c
c              expand each list of nodes. check for invalid entries
c
      if( local_debug1 ) write(out,9010) lenlst
      allocate( list_a(lenlst), list_b(lenlst) )
c
      icn    = 0
      iplist = 1
      count  = 1
      badlist = .false.
      do while ( iplist .ne. 0 )
         call trxlst( intlst, lenlst, iplist, icn, node )
         if( node .gt. nonode ) then
            param = node
            badlist = .true.
            call errmsg(16,param,dums,dumr,dumd)
            cycle
         end if
         if( node .le. 0 )  then
            param = node
            badlist = .true.
            call errmsg(58,param,dums,dumr,dumd)
            cycle
         end if
         list_a(count) = node
         count = count + 1
      end do  ! while loop to process nodes in list
      if( badlist ) return
c
      icn    = 0
      iplist = 1
      count  = 1
      badlist = .false.
      do while ( iplist .ne. 0 )
         call trxlst( mpc_intlst, mpc_lenlst, iplist, icn, node )
         if( node .gt. nonode ) then
            param = node
            badlist = .true.
            call errmsg(16,param,dums,dumr,dumd)
            cycle
         end if
         if( node .le. 0 )  then
            param = node
            badlist = .true.
            call errmsg(58,param,dums,dumr,dumd)
            cycle
         end if
         list_b(count) = node
         count = count + 1
      end do  ! while loop to process nodes in list
      if( badlist ) return
c
      if( local_debug1 ) then
         write(out,9019)
         do i = 1, lenlst
           write(out,9020) list_a(i), list_b(i)
         end do
      end if
c
c              process each pair of MPC nodes for each released dof. 
c               1. loop over nodes in the list (both lists already 
c                  checked for same number entries)
c               2. loop over 3 dof (u,v,w). skip if dof not listed
c                  for release in the command
c               3. verify neither node & dof are currently being 
c                  released.
c               4. search the table of user-defined MPC equations
c                   - a valid entry can only have 2 nodes
c                   - node pairs of the MPC must match user defined
c                     pair. see comments in release_cons_scan
c                   - node pairs and dof must match with release
c                     input
c                   - add both nodes that dof to the lists of nodes
c                     and dof to be released. this is same as adding
c                     two absolute constraints to the table since the
c                     MPC "force" on each node must be relaxed to zero.
c               5. columns in the user_mpc_table are marked for
c                  deletion. Set num terms = 0. table will be 
c                  compacted when all done.
c
      if( .not. allocated ( total_lagrange_forces ) ) then
        write(out,9120)  ! we're in a really bad condition
        call die_abort
      end if                         
c
      do i = 1, nodes_lst_a  ! same # as in list b
c
       node_a = list_a(i)
       node_b = list_b(i)
       minab  = min( node_a, node_b )
       maxab  = max( node_a, node_b )
       bad    = .false.
       if( local_debug1 ) write(out,9040) minab, maxab
c
       do idof = 1, 3  ! over u, v, w
c
         if( .not. release_flags(idof) ) cycle ! not mentioned in cmd
         nsteps_a = release_cons_table(idof,node_a)%num_release_steps
         nsteps_b = release_cons_table(idof,node_b)%num_release_steps
         t1 = nsteps_a .ne. 0
         t2 = nsteps_b .ne. 0
         if( t1 .or. t2 ) then  ! already being released
              write(out,9100) node_a, node_b, dof_names(idof)
              bad = .true.
              cycle
         end if
c
         found_mpc_entry = .false.
         do j = 1, num_user_mpc
c
           if( user_mpc_table(j)%num_terms .ne. 2 ) cycle
           mpc_dof_1 = user_mpc_table(j)%dof_list(1)
           mpc_dof_2 = user_mpc_table(j)%dof_list(2)
           if( debug ) write(out,9050) j, mpc_dof_1, mpc_dof_2
           ok = mpc_dof_1 .eq. idof  .and.  mpc_dof_2 .eq. idof
           if( .not. ok ) cycle ! not same dof as release cmd
           node_1 = user_mpc_table(j)%node_list(1)
           node_2 = user_mpc_table(j)%node_list(2)
           mpc_constant = user_mpc_table(j)%constant
           min12 = min( node_1, node_2 ) ! use min,max for easy check
           max12 = max( node_1, node_2 ) ! if 2 release nodes are same
           ok =  minab .eq. min12  .and.  maxab .eq. max12
           if( local_debug1 ) write(out,9060) min12, max12
           if( .not. ok ) cycle ! no match w/ release node pair
c
c                        release cmd node pair and dof matches 
c
           sdof = 3 * (node_a-1) + idof
           release_cons_table(idof,node_a)%num_release_steps =
     &            release_cons_steps
           release_cons_table(idof,node_a)%remaining_steps_for_release =
     &            release_cons_steps
           react = -total_lagrange_forces(sdof) + ifv(sdof)
           total_lagrange_forces(sdof) = zero
           d_lagrange_forces(sdof) = zero
           i_lagrange_forces(sdof) = zero
           release_cons_table(idof,node_a)%reaction_force = react
           count_stored = count_stored + 1
c
           sdof = 3 * (node_b-1) + idof
           release_cons_table(idof,node_b)%num_release_steps =
     &            release_cons_steps
           release_cons_table(idof,node_b)%remaining_steps_for_release =
     &            release_cons_steps
           react = -total_lagrange_forces(sdof) + ifv(sdof)
           total_lagrange_forces(sdof) = zero
           d_lagrange_forces(sdof) = zero
           i_lagrange_forces(sdof) = zero
           release_cons_table(idof,node_b)%reaction_force = react
           count_stored = count_stored + 1
c
           user_mpc_col_list(j) = 1
           compact_mpc_table = .true.
           found_mpc_entry = .true.
           if( local_debug1 ) write(out,9030) node_a, node_b,
     &                        dof_names(idof)
c
           if( local_debug1 ) then
             if( idof .eq. 2 .and. 
     &          (node_a .eq. 1 .or. node_a .eq. 85) ) then
              write(out,*) '.... data for node 1 v release...'
              sdof = 3 * (1-1) + idof
              react = release_cons_table(idof,node_a)%reaction_force
              write(out,*) '.......sdof,rload,ifv, LF:',sdof,
     &             rload(sdof),ifv(sdof), react
              write(out,*) '.... data for node 85 v release...'
              sdof = 3 * (85-1) + idof
              react = release_cons_table(idof,node_b)%reaction_force
              write(out,*) '.......sdof,rload,ifv, LF:',sdof,
     &              rload(sdof),ifv(sdof), react
             end if
           end if ! on local_debug1
c                 
         end do ! over user_mpc_table entries  
c
         if( found_mpc_entry ) cycle
         write(out,9110) node_a, node_b, dof_names(idof)
         bad = .true.
c
       end do ! over idof = 1, 2, 3list of release MPC node pairs
c
       if( bad )  num_error = num_error + 1
c
      end do  ! over nodes in lists           
c
      return
c
 9000 format(2x,"... entered release_cons_mpc ...")
 9010 format(10x,"... num nodes in each list: ",i0)
 9019 format(10x,"... paired nodes in lists")
 9020 format(20x,2i8)
 9030 format(15x,"... releasing nodes: ",i0,1x,i0,' dof: ', a1)
 9040 format(15x,"... minab, maxab: ",i0,1x,i0)
 9050 format(15x,"j, mpc_dof_1, mpc_dof_2:",3(1x,i0))
 9060 format(15x,"... min12, max12: ",2(1x,i0))
 9100 format(/1x,'>>>>> error: nodes: ',2(1x,i0),'dof: ',a1,'.',
     & /14x,'already in release state' )
 9110 format(/1x,'>>>>> error: nodes: ',i0,1x,i0,' dof: ',a1,'.',
     & ' do not have a defined MPC in constraints input.',/ )
 9120 format(/1x,'>>>>> fatal error: inconsistent data structures',
     &  /14x,'in release_cons_mpc. job terminated.',// )
 9130 format(/1x,'>>>>> error: the number of nodes in each list is ',
     & 'not identical.....'
     & /14x,'1st and 2nd list terms: ',2i7,/)
c
      end subroutine release_cons_mpc
c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine release empty table               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/23/03                   *
c     *                                                              *
c     *     skips to the end of a logical line of the input file     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine release_empty_cons_table
c
c              delete an empty release constraints table
c
      logical :: found_entry
      integer :: i, n1, n2, n3
c
      if(  .not. allocated( release_cons_table ) ) return
c
      found_entry = .false.
      do i = 1, nonode
          n1 = release_cons_table(1,i)%remaining_steps_for_release
          n2 = release_cons_table(2,i)%remaining_steps_for_release
          n3 = release_cons_table(3,i)%remaining_steps_for_release
          if( n1 .ne. 0 ) found_entry = .true.
          if( n2 .ne. 0 ) found_entry = .true.
          if( n3 .ne. 0 ) found_entry = .true.
          if( found_entry ) exit
      end do
      if( .not. found_entry ) then
           deallocate( release_cons_table )
           if( debug )  write(out,*) '... deleted empty table ...'
      end if
c
      return
c
      end subroutine release_empty_cons_table

      end subroutine release_constraints
c     ****************************************************************
c     *                                                              *
c     *           subroutine release_cons_update_constraints         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/24/03                   *
c     *                                                              *
c     * remove absolute constraint on 1 dof from constraint data     *
c     * structure                                                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine release_cons_update_constraints( sdof )
      use global_data ! old common.main
c
      use main_data, only : cnstrn_in, cnstrn
      use damage_data, only : csttail
c
      implicit none
c
      integer :: sdof
c
      integer :: cst_ptr, above_ptr
      double precision :: d32460
      data d32460 / 32460.0d00 /
c
c               traverse the singly-linked list of constraints in the
c               model. If the user has input new constraints, then some
c               of the previously released nodes may be re-constrained.
c               If a released node has been re-constrained, remove the
c               constraint.
c
c               set up indexes for constraint linked list
c
      cst_ptr   = csthed
      above_ptr = -1
c
c               enter top of constraint linked list loop. run till
c               dof is found or end of constraints.
c
      do while ( cst_ptr .ne. -1 )
c
      if( cst_ptr .eq. sdof ) then
         cnstrn(cst_ptr)    = d32460
         cnstrn_in(cst_ptr) = d32460
c
         if( above_ptr .eq. -1 ) then
c
c                           at top of list.  move head pointer.
c
            csthed          = cstmap(cst_ptr)
            cstmap(cst_ptr) = 0
            cst_ptr         = csthed
         else
c
c                           in middle of list. correct link indexes
c
            cstmap(above_ptr) = cstmap(cst_ptr)
            cstmap(cst_ptr)   = 0
            if ( csttail .eq. cst_ptr ) csttail = above_ptr
            cst_ptr = cstmap(above_ptr)
         end if
         exit
      end if
c
c                    examine next constraint in list
c
      above_ptr = cst_ptr
      cst_ptr   = cstmap(cst_ptr)
c
      end do
c
      return
      end



c     ****************************************************************
c     *                                                              *
c     *                      subroutine incon                        *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 12/3/2019 rhd              *
c     *                                                              *
c     *     input of nodal displacement constraints                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine incon( sbflg1, sbflg2, olddof )
      use global_data ! old common.main
      use allocated_integer_list
c
      use main_data, only : trn, trnmat, cnstrn_in,
     &                      inverse_incidences, force_solver_rebuild
      use mod_mpc, only : mpcs_exist, num_user_mpc, user_mpc_table
      use damage_data, only : csttail
c
      implicit integer (a-z)
c
c                       locally allocated
c
      integer, allocatable :: intlst(:)
      integer :: list_size
      double precision :: 
     &  convec(mxndof), trans(mxndof,mxndof), tval, cval, dumd, 
     &  rlen1, rlen2, rlen3, t11, t12, t13,
     &  t21, t22, t23, t31, t32, t33
      real :: dumr, rnode
      character :: dums, curtyp *1
      logical :: sbflg1, sbflg2, skew, inpflg(mxndof), defcon(3),
     &           rflag1, rflag2
      logical, save :: cons_defined = .false. 
      logical, external :: matchs, endcrd, true, numd
      double precision, parameter :: zero = 0.d0, one = 1.d0,
     &                               d32460 = 32460.0d0,
     &                               rottol = 0.0001d0
c
c                       if sub flag 1 is on, there is reentry into
c                       incon after an error in constraints input.
c
      allocate( intlst(10) )
c
      if( sbflg1 ) then
         call errmsg(71,dum,dums,dumr,dumd)
         go to 710
      end if
c
c                       make sure that any previous transformation
c                       matrices created due to contact are re-applied
c                       to rotate the corresponding degrees of freedom
c                       to global coordinates.
c
      call contact_remove (.true.)
c
c                       warn user that all constraints are being
c                       destroyed and that all constraints must be
c                       be re-defined....  deallocate also deletes
c                       allocatable sub-objects (F2003)
c
      new_constraints = .true.
      force_solver_rebuild = .true.
c
      if ( cons_defined ) then
         call errmsg(224,dum,dums,dumr,dumd)
         call errmsg2(59,dum,dums,dumr,dumd)
         if (mpcs_exist) then
            if (allocated(user_mpc_table)) deallocate(user_mpc_table)
            num_user_mpc = 0
            mpcs_exist = .false.
         end if
      endif
      cons_defined = .true.
c
c                       initialize the constraints link list.
c
      csthed = -1
      olddof = 0
      do i = 1, nodof
         cstmap(i)    = 0
         cnstrn_in(i) = d32460
      end do
c
c                       initialize transformation matrix indexes
c                       and flags. deallocate any old transformation
c                       matrices
c
      do i = 1, nonode
         if ( trn(i) ) call allo_trnmat(i,2,dum)
         trn(i) = .false.
      end do
c
c                      initialize defined constraint flags to catch
c                      problems with one or more totally unconstrained
c                      direction
c
      defcon(1) = .false.
      defcon(2) = .false.
      defcon(3) = .false.
c
c

 710  call readsc
 711  continue
c
c                       look for transformation matrix input. set
c                       flag depending on whether or not it is found.
c
      if( matchs('transformation',5) ) then
         if( matchs('matrix',5) ) call splunj
         if( matchs('nodes',4) ) call splunj
         skew = .true.
      else
         skew = .false.
      end if
c
c                       look for dump command to dump out all the
c                       constraints
c
      if ( matchs('dump',4) ) then
         call con_dump(olddof)
         go to 710
      end if
c
c                       look for a command to impose constraints
c                       on all nodes on a plane
c
      if ( matchs('plane',4) ) then
         call inconplane( olddof, defcon )
         go to 710
      end if
c
c                       look for a command to start input of multi-point
c                       constraint equations. returns when line does not
c                       start with an integer. Did a reset, true before
c                       return.
c                       delete data if no mpc equations actually input
c                       correctly.
c
      if ( matchs('multipoint',5) ) then
         call incon_mpcs
         if( num_user_mpc == 0 ) then
           mpcs_exist = .false.
           deallocate( user_mpc_table )
         end if
         go to 711
      end if
c
c                       input node list
c
      deallocate( intlst )
      allocate( intlst(10) )
      call trlist_allocated( intlst, list_size, nonode, 
     &                       lenlst, errnum)
c
c                       branch on the return code from trlist. a
c                       value of 1 indicates no error. a value of
c                       2 indicates that the parse rules failed in
c                       the list. a value of 3 indicates that the
c                       list overflowed its maximum length of mxlsz.
c                       in these last two cases, the illegal list
c                       will be ignored and a new node list will
c                       be sought. a value of 4 indicates that no list
c                       was found. in this case, constraints input
c                       has ceased.
c
      if( errnum .eq. 2 ) then
         param= 1
         call errmsg(24,param,dums,dumr,dumd)
         go to 710
      else if( errnum .eq. 3 ) then
         param = 2
         call errmsg(24,param,dums,dumr,dumd)
         go to 710
      else if( errnum .eq. 4 ) then
         go to 9999
      else
         if( errnum .eq. 1 ) then
            call backsp(1)
            if( true(dummy) ) go to 715
         end if
         param = 3
         call errmsg(24,param,dums,dumr,dumd)
         go to 710
      end if
c
c
c **********************************************************************
c *                                                                    *
c *                     input of transformation matrix at list of      *
c *                     nodes which defines constraint compatable      *
c *                     global coordinates.                            *
c *                                                                    *
c **********************************************************************
c
c
 715  if( skew ) then
c
c                       initialize nodal input transformation array
c                       and the input flag vector.
c
         do i = 1, mxndof
            inpflg(i) = .false.
            do j = 1, mxndof
               if( i .eq. j ) then
                  trans(i,j) = one
               else
                  trans(i,j) = zero
               end if
            end do
         end do
c
c                       branch on row of matrix currently being input.
c
 717     if( matchs('row_1',5) ) then
            row = 1
            go to 719
         end if
c
         if( matchs('row_2',5) ) then
            row = 2
            go to 719
         end if
c
         if( matchs('row_3',5) ) then
            row = 3
            go to 719
         end if
c
         if( matchs(',',1) ) go to 718
c
c                       if there is an end of card, matrix input has
c                       ended. branch to store the matrix globally.
c                       if not, ignore the current entity and search
c                       for another row to input. check validity of
c                       3x3 rotation matrix.
c
         if( endcrd(dum) ) then
          rlen1 = sqrt( trans(1,1)**2 + trans(1,2)**2 +trans(1,3)**2 )
          rlen2 = sqrt( trans(2,1)**2 + trans(2,2)**2 +trans(2,3)**2 )
          rlen3 = sqrt( trans(3,1)**2 + trans(3,2)**2 +trans(3,3)**2 )
          rflag1 = abs(rlen1-one) .gt. rottol .or.
     &             abs(rlen2-one) .gt. rottol .or.
     &             abs(rlen3-one) .gt. rottol
          t11 = trans(1,1)**2 + trans(2,1)**2 + trans(3,1)**2
          t12 = trans(1,1)*trans(1,2) + trans(2,1)*trans(2,2) +
     &          trans(3,1)*trans(3,2)
          t13 = trans(1,1)*trans(1,3) + trans(2,1)*trans(2,3) +
     &          trans(3,1)*trans(3,3)
          t21 = trans(1,1)*trans(1,2) + trans(2,1)*trans(2,2) +
     &          trans(3,1)*trans(3,2)
          t22 = trans(1,2)**2 + trans(2,2)**2 + trans(3,2)**2
          t23 = trans(1,2)*trans(1,3) + trans(2,2)*trans(2,3) +
     &          trans(3,2)*trans(3,3)
          t31 = trans(1,1)*trans(1,3) + trans(2,1)*trans(2,3) +
     &          trans(3,1)*trans(3,3)
          t32 = trans(1,2)*trans(1,3) + trans(2,2)*trans(2,3) +
     &          trans(3,2)*trans(3,3)
          t33 = trans(1,3)**2 + trans(2,3)**2 + trans(3,3)**2
          rflag2 = abs(t11-one) .gt. rottol .or.
     &             abs(t22-one) .gt. rottol .or.
     &             abs(t33-one) .gt. rottol .or.
     &             abs(t12) .gt. rottol     .or.
     &             abs(t13) .gt. rottol     .or.
     &             abs(t21) .gt. rottol     .or.
     &             abs(t23) .gt. rottol     .or.
     &             abs(t31) .gt. rottol     .or.
     &             abs(t32) .gt. rottol
          if ( rflag1 .or. rflag2 ) then
                  call errmsg( 266, param, dums, dumr, dumd )
          end if
          go to 730
         else
            call errmsg(84,dum,dums,dumr,dumd)
            if( true(dum) ) go to 717
         end if
c
c                       if there is a comma at the end of a line, the
c                       input line is continued.
c
 718     continue
         if( endcrd(dum) ) then
            call readsc
         end if
         go to 717
c
c                       check to make sure that the current row has
c                       not been input twice. set input flag if not.
c
 719     if( inpflg(row) ) then
            param = row
            call errmsg(135,param,dums,dumr,dumd)
            go to 717
         else
            inpflg(row) = .true.
         end if
c
c                       store the row input.
c
         col= 0
 720     if( .not. numd(tval) ) then
            go to 717
         else
            col = col+1
            if( col .gt. mxndof ) then
               col = mxndof
               write(out,9073) 'columns ', node, 'columns '
               num_error = num_error + 1
               if( true(dum) ) go to 717
            end if
            trans(row,col) = tval
         end if
         go to 720
c
c
c **********************************************************************
c *                                                                    *
c *                     input of nodal constraint values, in con-      *
c *                     straint compatable global coordinates.         *
c *                                                                    *
c **********************************************************************
c
c
      else
c
c                       initialize the temporary constraint vector.
c
         do i = 1,mxndof
            convec(i) = d32460
            inpflg(i) = .false.
         end do
c
c                       branch on dof type input.
c
 726     if( matchs('u',1) ) then
            dofn= 1
            curtyp= 'u'
            go to 728
         end if
c
         if( matchs('v',1) ) then
            dofn= 2
            curtyp= 'v'
            go to 728
         end if
c
         if( matchs('w',1) ) then
            dofn= 3
            curtyp= 'w'
            go to 728
         end if
c
         if( matchs(',',1) ) go to 727
c
c                       if there is an end of card, constr. input has
c                       ended. branch to store the temp. vec. globally.
c                       if not, ignore the current entity and search
c                       for another dof to input.
c
         if(endcrd(dum)) then
            go to 730
         else
            call errmsg(72,dum,dums,dumr,dumd)
            if(true(dum)) go to 726
         end if
c
c                       if there is a comma at the end of a line, the
c                       input line is continued.
c
 727     continue
         if(endcrd(dum)) then
            call readsc
         end if
         go to 726
c
c                       store the constraints in the temporary vector
c
 728     if(matchs('=',1)) call splunj
         if(.not.numd(cval)) then
            call errmsg(73,dum,dums,dumr,dumd)
         else
c
c                       make sure that the current dof has not
c                       already been input on the same nodal con-
c                       straint command.
c
            if( inpflg(dofn) ) then
               call errmsg(109,dum,curtyp,dumr,dumd)
               go to 726
            end if
            defcon(dofn) = .true.
            inpflg(dofn) = .true.
            convec(dofn) = cval
         end if
         go to 726
      end if
c
c
c **********************************************************************
c *                                                                    *
c *                     store constraints or transformation matrices   *
c *                     globally.                                      *
c *                                                                    *
c **********************************************************************
c
c
 730  icn    = 0
      iplist = 1
 731  call trxlst( intlst, lenlst, iplist, icn, node )
c
c                       check that the list node does not exceed
c                       the number of nodes in the structure.
c
         if( node .gt. nonode ) then
            param= node
            call errmsg(16,param,dums,dumr,dumd)
            go to 734
         end if
c
c                       check that the list node is not negative.
c
         if( node .lt. 0 )  then
            param= node
            call errmsg(58,param,dums,dumr,dumd)
            go to 734
         end if
c
         felem = inverse_incidences(node)%element_list(1)
         type  = iprops(1,felem)
         ndof  = iprops(4,felem)
c
c
c **********************************************************************
c *                                                                    *
c *                     store transformation matrix globally           *
c *                                                                    *
c **********************************************************************
c
c
         if( skew ) then
c
c                       set flag that there is a trn. mat. for the node.
c                       if needed, allocate matrix for this node.
c                       if the node already has had a trn. mat. input,
c                       then overwrite the previous matrix.
c
            if( .not. trn(node) ) then
               trn(node) = .true.
               call allo_trnmat ( node, 1, dum )
            endif
c
c                       store the transformation matrix globally
c                       for the node.
c
            do row = 1, 3
               trnmat(node)%mat(row,1) = trans(row,1)
               trnmat(node)%mat(row,2) = trans(row,2)
               trnmat(node)%mat(row,3) = trans(row,3)
            end do
c
c
c
c **********************************************************************
c *                                                                    *
c *                     store the constraint vector globally           *
c *                                                                    *
c *                     ndof  -   number of dof per node               *
c *                                                                    *
c *                                                                    *
c *                                                                    *
c **********************************************************************
c
c
         else
c
            do i = 1, ndof
               if( convec(i).eq.d32460 ) cycle
               dof = dstmap(node)+i-1
c
c                       make sure this dof hasn't been previously
c                       constrained.
c
               if( cnstrn_in(dof) .ne. d32460 ) then
                  param= dof
                  rnode = node
                  call errmsg(136,param,dums,rnode,dumd)
                  cycle
               end if
c
c                       place the dof in the constraint mapping data
c                       structure.
c
               if( csthed .eq. -1 ) then
                  csthed = dof
               else
                  cstmap(olddof) = dof
               end if
               olddof         = dof
               cnstrn_in(dof) = convec(i)
            end do
c
         end if
c
 734  if( iplist .ne. 0 ) go to 731
c
      go to 710
c
c
c **********************************************************************
c **********************************************************************
c
c
 9999 sbflg1         = .true.
      sbflg2         = .true.
      cstmap(olddof) = -1
      csttail        = olddof
c
c                  user must specify at least 1 u, v, w
c                  constraint, else warning message
c
      do i = 1, 3
         if ( .not. defcon(i) ) then
            call errmsg(186,dum,dums,dumr,dumd)
            exit
         end if
      end do
c
      return
c
 9073 format(/1x,'>>>>> error: the number of ',a8,' input for a row of',
     &           ' the transformation'/14x,'matrix of node ',i7,
     &           ' exceeds the maximum number of'/14x,'degrees of',
     &           ' freedom allowed any node. the number of'/14x,a8,
     &           ' stored will be the above maximum.'/)
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine inconplane                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 4/12/2018 rhd              *
c     *                                                              *
c     *     process plane command option for constraints             *
c     *     constraints imposed on a specified plane                 K*
c     *                                                              *
c     ****************************************************************
c
c
      subroutine inconplane( olddof, defcon )
      use global_data ! old common.main
c
      use main_data, only : cnstrn_in, crdmap
c
      implicit integer (a-z)
c
c                       parameters
c
      logical defcon(*)
c
c
c                       locally allocated
c
      double precision
     & convec(mxndof), dumd, zero, proximity_distance,
     & d32460, xmin, xmax, ymin, ymax, zmin, zmax, x, y, z,
     & coordtol(3), plane_coord, ctol, coordvec(3)
      real dumr
      character(len=1) :: dums
      logical inpflg(mxndof), verify
      logical matchs, numd, local_debug, endcrd
      data zero, d32460 / 0.0d00, 32460.0 /
      data ctol, local_debug / 0.00001, .false. /
c
c                  find min and max coordinates for the body
c
      xmin =  1.0e30
      xmax = -1.0e30
      ymin =  1.0e30
      ymax = -1.0e30
      zmin =  1.0e30
      zmax = -1.0e30
c
      do node = 1, nonode
        x    = c(crdmap(node))
        y    = c(crdmap(node)+1)
        z    = c(crdmap(node)+2)
        xmax = max( xmax, x )
        xmin = min( xmin, x )
        ymax = max( ymax, y )
        ymin = min( ymin, y )
        zmax = max( zmax, z )
        zmin = min( zmin, z )
      end do
c
      coordtol(1) = abs( xmax-xmin ) * ctol
      coordtol(2) = abs( ymax-ymin ) * ctol
      coordtol(3) = abs( zmax-zmin ) * ctol
c
      if ( local_debug ) then
        write(*,*) '>> inside inconplane:'
        write(*,*) ' xmin,max: ',xmin,xmax
        write(*,*) ' ymin,max: ',ymin,ymax
        write(*,*) ' zmin,max: ',zmin,zmax
        write(*,*) ' xtol, ytol, ztol: ', coordtol
      end if
c
c
c                  which plane to constrain
c
      plane = 0
      if ( matchs('x',1) ) plane = 1
      if ( matchs('y',1) ) plane = 2
      if ( matchs('z',1) ) plane = 3
      if ( plane .eq. 0 ) then
        call errmsg2( 4, dum, dums, dumr, dumd )
        return
      end if
c
c                  coordinate of plane is zero by default
c
      plane_coord = zero
      if ( matchs('=',1) ) then
       if ( .not. numd(plane_coord) ) then
         call errmsg2( 5, dum, dums, dumr, dumd )
         return
       end if
      end if
c
      if ( local_debug ) then
        write(*,*) '>> plane id, coord: ',plane,plane_coord
      end if
c
      if( matchs('proximity',4) ) then
        if( numd( proximity_distance ) ) then
          coordtol(1:3) =  proximity_distance
        end if
      end if
c
c                  get the type of constraints to impose on nodes
c                  that lie on the plane, store values, return
c
c                  fixed, symmetry or list of u, v, w values
c
      convec(1:3) = d32460
      inpflg(1:3) = .false.
      verify = .false.
      if ( matchs('verify',3) ) verify = .true.

c
      if ( matchs('fixed',3) ) then
        if ( matchs('verify',3) ) verify = .true.
        convec(1:3) = zero
        inpflg(1:3) = .true.
        if( matchs('proximity',4) ) then
          if( numd( proximity_distance ) ) then
            coordtol(1:3) =  proximity_distance
          end if
        end if
        write(out,9100) coordtol(1:3)
        if ( matchs('verify',3) ) verify = .true.
        call inconplane_store
        return
      end if
c
      if ( matchs('symmetry',3) ) then
        if ( matchs('verify',3) ) verify = .true.
        convec(plane) = zero
        inpflg(plane) = .true.
        if( matchs('proximity',4) ) then
          if( numd( proximity_distance ) ) then
            coordtol(1:3) =  proximity_distance
          end if
        end if
        write(out,9100) coordtol(1:3)
        if ( matchs('verify',3) ) verify = .true.
        call inconplane_store
        return
      end if
c
      do   !  u, v, w values
       if( matchs('proximity',4) ) then
         if( numd( proximity_distance ) ) then
           coordtol(1:3) =  proximity_distance
         end if
         cycle
        end if
        if ( matchs('verify',3) ) then
           verify = .true.
           cycle
        end if
        if( endcrd( dummy ) ) then
          write(out,9100) coordtol(1:3)
          call inconplane_store
          return
        end if
        if ( matchs('u',1) ) then
          dofn = 1
        elseif ( matchs('v',1) ) then
          dofn = 2
        elseif ( matchs('w',1) ) then
          dofn = 3
        else
          call errmsg2( 6, dum, dums, dumr, dumd )
          return
        end if
c
        inpflg(dofn) = .true.
        if ( matchs('=',1) ) call splunj
        if ( .not. numd(convec(dofn)) ) then
          call errmsg2( 7, dum, dums, dumr, dumd )
          return
        end if
      end do
c
 9100 format(10x,">>>> Distances (x,y,z) for proximity test: ",
     &    3f15.6)

      contains  ! avoids go tos in above code.
c     ========
c     ****************************************************************
c     *                                                              *
c     *                subroutine inconplane_store                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 07/3/2014                  *
c     *                                                              *
c     *     store constraints imposed by plane command               *
c     *                                                              *
c     ****************************************************************
c
      subroutine inconplane_store
c
c                  search for nodes that lie on the specified plane.
c                  impose constraints on nodes. output verification list
c                  if requested. we don't warn here about constraints
c                  imposed on nodes that already have constraints.
c
c
      if ( local_debug ) then
        write(*,*) '>> constraints to impose on matching nodes'
        write(*,*) ' convec: ', convec
        write(*,*) ' inpflg: ', inpflg
      end if
      if ( verify ) write(out,9000)
      con_node_count = 0
c
      do node = 1, nonode
        coordvec(1) = c(crdmap(node))
        coordvec(2) = c(crdmap(node)+1)
        coordvec(3) = c(crdmap(node)+2)
        if ( abs( plane_coord -coordvec(plane) ) .gt.
     &        coordtol(plane) ) cycle
        con_node_count = con_node_count + 1
        do i = 1, 3
          if( convec(i) .eq. d32460 ) cycle
          dof = dstmap(node) + i -1
          if( csthed .eq. -1 ) then
             csthed = dof
          else
             cstmap(olddof) = dof
          end if
          olddof         = dof
          cnstrn_in(dof) = convec(i)
          defcon(i)      = .true.
          if ( verify ) then
            write(out,9100) node,i,convec(i)
          end if
        end do
      end do
      write(out,9200) con_node_count
c
      return
c
 9000 format(//,1x,
     &'**** Constraints imposed by just entered plane command',/,/,
     & 5x,' node   dof   constraint value')
 9100 format(5x,i5,3x,i1,8x,f12.7)
 9200 format(/,'>> Plane command applied constraints to: ',i5,
     & ' nodes...',//)
c
      end subroutine inconplane_store
c
      end subroutine inconplane


c     ****************************************************************
c     *                                                              *
c     *                    subroutine incon_mpcs                     *
c     *                                                              *
c     *                       written by : bjb                       *
c     *                                                              *
c     *                   last modified : 08/22/2017 rhd             *
c     *                                                              *
c     *     this subroutine supervises and conducts the input of     *
c     *     multi-point constraint equations                         *
c     *     - update to following fix so leading term can have a +,- *
c     *     - major error fix 7/29/11 to make eqn arrangement match  *
c     *       that for tied contact. Byron implemented enforcement   *
c     *       of mpcs assuming they were in form of tied contact     *
c     *       arrangement. see notes below                           *
c     *     - fixed several errors when MPCs are changed during      *
c     *       an analysis
c     *                                                              *
c     ****************************************************************
c
c
      subroutine incon_mpcs
      use global_data ! old common.main
      use mod_mpc, only : mpcs_exist, num_user_mpc, user_mpc_table
      use main_data, only : modified_mpcs
c
      implicit none
c
c              locals
c
      integer, parameter :: max_trm = 100 ! limit on terms per mpc eqn
c
c              scan support
c
      integer :: dumi
      real    :: dumr
      double precision ::  dumd
      character(len=1) :: dums
      logical, external :: matchs, realn, numr, integr, endcrd, true
      logical, parameter :: debug = .false.
c
c              local storage for translating, checking
c
      real :: pre_term_multiplier, now_term_coeff, eqn_rhs_constant,
     &        rtemp
      integer, dimension(:) :: eqn_nodes(max_trm), eqn_dofs(max_trm)
      integer :: now_dof, ierr, idummy(1), now_node, nterm, nmpc
      logical :: found_mpc, found_rigid
      real, dimension(:) :: eqn_coeffs(max_trm)
c
c              Set a flag so that we know later on the constraints have
c              been touched
c
      modified_mpcs = .true.
c
c              notify user we're discarding any existing mpcs
c              in the mpc module. allocate permanent data structure
c              for user mpcs in the module (see comments in mod_mpc)
c              F2003 - also deallocates allocated sub-objects.
c
      if( mpcs_exist ) then
         call errmsg2( 59, dumi, dums, dumr, dumd )
         if( allocated(user_mpc_table) ) deallocate( user_mpc_table )
         num_user_mpc = 0
      end if
c
      allocate( user_mpc_table(max_mpc), stat=ierr )
      if( ierr .ne. 0 ) then
         call errmsg2( 47, dumi, dums, dumr, dumd )
         call die_abort
      end if
c
      mpcs_exist = .true.
      nmpc  = 0
c
c               loop to process next mpc equation. line must start
c               with an +, - or integer. if not, treat as end of
c               mpc equations
c
      do    ! over mpc equations in inout
         call readsc
         found_mpc = .false.
         call incon_mpcs_rigid ! process rigid option separately
         if( found_rigid ) cycle
         call reset
         call scan
         if( matchs('+',1) ) then
            found_mpc = .true.
         elseif( matchs('-',1) ) then
            found_mpc = .true.
         elseif( integr(idummy) ) then
            found_mpc = .true.
         else
            found_mpc = .false.
         end if
c
         call reset
         if( true(idummy) )  call splunj
         if( .not. found_mpc ) return
c
c              loop across line to extract each node number, its
c              dof and multipler, the (+,-) between terms and the
c              rhs const. for now, the constant = 0.0. when done.
c              store the equation terms in the module data structure.
c              loop over multiple lines continued by commas.
c
c              for curremt mpc implementation, the rhs
c              constant must = 0.0
c
c              since we have done a line reset, the leading term may
c              now have a pre-multiplier ( '+' or '-' ) just like other
c              terms
c
         nterm          = 0
         now_node       = 0
         now_dof        = 0
         now_term_coeff = 0.0
         nmpc = nmpc + 1
         if( nmpc .gt. max_mpc ) call  incon_mpcs_resize( 1, idummy ) 
c
         do ! all terms for this mpc eqn
c
            if( matchs('=',1) ) then
              eqn_rhs_constant = 0.0
              if( numr(rtemp) ) then
                 eqn_rhs_constant = rtemp
                 if( rtemp .ne. 0.0 ) then
                  call errmsg( 320, 1, dums, dumr, dumd )
                  call incon_flushline
                  nmpc = nmpc - 1
                  exit ! from processing this mpc eqn
                 end if
              end if
              call incon_mpcs_store( nterm, eqn_rhs_constant,
     &                     eqn_nodes, eqn_dofs, eqn_coeffs )
               exit  !  start at loop for new mpc eqn
            end if
c
            nterm = nterm + 1
c
            if( nterm .gt. max_trm ) then
               call errmsg2( 43, max_trm, dums, dumr, dumd )
               call die_abort
            end if
c
            pre_term_multiplier = 1.0
            if( matchs('+',1) ) then
                  pre_term_multiplier = 1.0
               else if (matchs('-',1) ) then
                  pre_term_multiplier = -1.0
            end if
c
c                     node number
c
            if( integr(now_node) ) then
               if( now_node .gt. nonode .or. now_node .le. 0  ) then
                  call errmsg2( 34, now_node, dums, dumr, dumd )
                  call reset
                  if( true(dumi) ) call splunj
                  call incon_flushline
                  nmpc = nmpc - 1
                  exit ! from processing this mpc eqn
               end if
               eqn_nodes(nterm) = now_node
            end if
c
c                    coefficient after node number. must be
c                    present and a real number.
c
            if( realn(now_term_coeff) ) then
             eqn_coeffs(nterm) = pre_term_multiplier * now_term_coeff
            else
             call errmsg( 319, 1, dums, dumr, dumd )
             call reset
             if( true(dumi) ) call splunj
             call incon_flushline
             nmpc = nmpc - 1
             exit ! from processing this mpc eqn
            end if
c
c                    get dof: u, v, w
c
            now_dof = 0
            if( matchs('u',1) ) then
                  now_dof = 1
               else if( matchs('v',1) ) then
                  now_dof = 2
               else if( matchs('w',1) ) then
                  now_dof = 3
            end if
            if( now_dof. eq. 0 ) then
             call errmsg( 210, 1, dums, dumr, dumd )
             call reset
             if( true(dumi) ) call splunj
             call incon_flushline
             nmpc = nmpc - 1
             exit ! from processing this mpc eqn
            end if
            eqn_dofs(nterm) = now_dof
c
c                comma followed by eol is continuation of current mpc
c
            if( matchs(',',1) ) then
               if( endcrd(dumi) ) then
                  call readsc
                  cycle  ! top of loop looking for terms in this mpc
               end if
            end if
c
         end do ! loop looking for terms in an mpc
      end do  ! loop to process all mpc eqns
c
      return
c
c    Notes on storage and implementation of MPC equations.
c
c      (a) the equations must be homogeneous (rhs = 0)
c
c      (b) none of the mentioned dof can also have a zero
c          or non-zero (absolute) constraint applied (which
c          effectively makes the equation non-homogeneous.
c
c      (c) equations must be normalized so that the multiplier
c          stored for the leading term is -1.0
c
c      (d) a simple form is available to rigidly couple nodes
c          (make u, v, w identical).
c
c           <node list> rigid <node list>
c
c    The tied-contact capability was implemented before user
c    defined mpcs. The tied contact arrangement of equations has
c    the form (for example):
c
c       ux = 0.25 ub + 0.3 uc + 0.4uc + 0.05 ud
c
c    Here node 'x' is located on an element face which has corner nodes
c    a, b, c, d. The u-displ at 'x' is then a linear combination of
c    displacements of the 4 face nodes with coefficients determined by
c    the geometric location of node 'x' on the face.
c
c    The tied contact input system stores the above equation in the form:
c
c     -1.0 ux + 0.25 ub + 0.3 uc + 0.4uc + 0.05 ud= 0.0
c
c    Note the leading coefficient is -1.0. Byron implemented the above
c    mpc assuming the leading term is -1.0 and that ux is the dependent
c    dof
c
c    Similarly, for two coincident nodes connected by tied contact say
c    20 and 50, the tied contact input system will store:
c
c      -1.0 u20 + 1.0 u50 = 0.0
c
c    Our user-defined mpc's must be transformed to follow the exact
c    storage for tied contact since Byron's code to enforce the mpcs
c    is based on the tied contact data structure.
c
c    A user mpc of the form:
c
c        53 1.0 u - 43 1.0 u = 0
c
c    to make the two u displacements equal must be stored in the form:
c
c         -1.0 u53 + 1.0 u43 = 0
c
c    The user can write this same equation as
c
c      - 53 1.0 u + 43 1.0 u = 0
c      - 53 -1.0 u - 43 1.0 u = 0
c      - 53 1.0 u - 43 -1.0 u = 0
c       53 -42.5 + 43 u 42.5 = 0
c
c    All forms must be stored as
c
c        -1.0 u53 + 1.0 u43 = 0.0
c
c    More complex example. The user wants the v-dispacement on all 4
c    nodes of an element face to sum = 0. For face nodes 23,
c    99, 103, 42 the user might write as input:
c
c            23 1.0 v + 99 1.0 v + 103 1.0 v + 1.0 42 v = 0
c
c    We must store this equation as:
c
c            -1.0 u23  - 1.0 v99 - 1.0 v103 - 1.0 v42 = 0.0
c
c    The code scans for the pre-multiplier sign (+,-), the
c    node number, the coefficient and the dof.
c    The node number, dof and the multiplier = pre-multiplier *
c    the coefficient are stored in simple vectors while scanning.
c
c    While storing into the global MPC data structure that Brian's
c    code uses, scale the coefficients so that the leading term
c    multiplier is -1.0.
c
      contains
c     ========
c     ****************************************************************
c     *                                                              *
c     *                  subroutine incon_mpcs_rigid                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/3/2019 rhd              *
c     *                                                              *
c     *     scan and store rigid constraints between 2 node lists    *
c     *                                                              *
c     ****************************************************************
      subroutine incon_mpcs_rigid()
c
      use allocated_integer_list
c
      logical, external :: match_exact, true, iszlst
      integer :: errnum, list_size_a, lenlst_a, lenlst_b, param, 
     &           dumi, list_size_b, nodes_lst_a,
     &           nodes_lst_b, nterms, node_list(2), idof,
     &           icn_a, icn_b, iplist_a, iplist_b, node_a, node_b,
     &           dof_list(2)
      integer, allocatable :: intlst_a(:), intlst_b(:)
      real :: dumr, constant, multipliers(2)
      double precision :: dumd
      character * 1 :: dums
c
c              <node list> rigid <node list>
c
c              line must start with integerlist and keyword rigid
c              otherwise reset to start of line for likely
c              processing of regular MPC equation of node with MPCs
c
      if( debug ) write(out,*) '..... entered incon_mpcs_rigid'
      found_rigid = .false.
c
      allocate( intlst_a(10), intlst_b(10) ) 
      call scan
      call trlist_allocated( intlst_a, list_size_a, nonode, 
     &                       lenlst_a, errnum )
c
c              branch on the return code from trlist. a
c              value of 1 indicates no error. a value of
c              2 indicates that the parse rules failed in
c              the list. a value of 3 indicates that the
c              list overflowed its maximum length of mxlsz.
c              in these last two cases, the illegal list
c              will be ignored and a new node list will
c              be sought. a value of 4 indicates that no list
c              was found. in this case, release cons input
c              has ceased.
c
      if( debug ) write(out,9010) errnum, lenlst_a
      if( errnum .eq. 4 ) return
      if( errnum .eq. 2 ) then
         param = 1
         call errmsg(24,param,dums,dumr,dumd)
         call scan_flushline
         return
      end if
c
      if( errnum .eq. 3 ) then
         param = 2
         call errmsg(24,param,dums,dumr,dumd)
         call scan_flushline
         return
      end if
c
      if( errnum .eq. 1 ) then
          call backsp(1)
          if( true(dumi) ) call splunj
      else
         write(out,9000)
         call die_abort
      end if
c
      if( .not. match_exact( 'rigid', 5 ) ) return
c
c              we appear to have a rigid constraint line. get the
c              matching list of nodes. both lists must have
c              the same number of terms.
c
      found_rigid = .true.
c
      call scan
      call trlist_allocated( intlst_b, list_size_b, nonode, 
     &                       lenlst_b, errnum )
      if( debug ) write(out,9010) errnum, lenlst_b
      if( errnum .eq. 4 ) return
      if( errnum .eq. 2 ) then
         param = 1
         call errmsg(24,param,dums,dumr,dumd)
         call scan_flushline
         return
      end if
      if( errnum .eq. 3 ) then
         param = 2
         call errmsg(24,param,dums,dumr,dumd)
         call scan_flushline
         return
      end if
      if( errnum .eq. 1 ) then
          call backsp(1)
          if( true(dumi) ) call splunj
      else
         write(out,9000)
         call die_abort
      end if
c
      nodes_lst_a = iszlst( intlst_a, lenlst_a )
      nodes_lst_b = iszlst( intlst_b, lenlst_b )
      if( debug ) write(out,9020) nodes_lst_a, nodes_lst_b
      if( nodes_lst_a .ne. nodes_lst_b ) then
           write(out,9100)  nodes_lst_a, nodes_lst_a
           num_error = num_error + 1
           return
      end if
c
c              traverse each list to extract pairs of nodes
c              and insert rigid constraint on u, v, w
c
      if( debug ) write(out,9030) nmpc
      icn_a    = 0
      iplist_a = 1
      icn_b    = 0
      iplist_b = 1
c 
      do while ( iplist_a .ne. 0 )  ! both lists same # terms
c
        call trxlst(intlst_a,lenlst_a,iplist_a,icn_a,node_a)
        call trxlst(intlst_b,lenlst_b,iplist_b,icn_b,node_b)
        if( debug ) write(out,*) '... node_a, node_b', node_a, node_b
c
        if( node_a .gt. nonode ) then
            param = node_a
            call errmsg(16,param,dums,dumr,dumd)
            cycle
        end if
        if( node_a .le. 0 )  then
            param = node_a
            call errmsg(58,param,dums,dumr,dumd)
            cycle
        end if
        if( node_b .gt. nonode ) then
            param = node_b
            call errmsg(16,param,dums,dumr,dumd)
            cycle
        end if
        if( node_b .le. 0 )  then
            param = node_b
            call errmsg(58,param,dums,dumr,dumd)
            cycle
        end if
c
        nterms         = 2
        constant       = 0.0   !  real not double
        node_list(1)   = node_a
        node_list(2)   = node_b
        multipliers(1) = -1.0  !  real not double
        multipliers(2) = 1.0
c
        do idof = 1, 3   !  u, v, w
          nmpc = nmpc + 1          
          if( nmpc .gt. max_mpc ) call incon_mpcs_resize( 1, idummy ) 
          dof_list(1) = idof
          dof_list(2) = idof
          call incon_mpcs_store( nterms, constant, node_list, 
     &                           dof_list, multipliers )
        end do ! over idof
c
      end do ! over node pair extraction
c
      if( debug ) write(out,9040) nmpc
c
      return
c
 9000 format(/1x,'>>>>> error: invalid return on trlist in ',
     &  /14x,'release_cons_scan. system error. job aborted.',//)
 9010 format('.... trlist done, errnum, lenlst: ',2i4)
 9020 format('.... nodes_lst_a and _b: ',2i6)
 9030 format('.... insert rigid MPOCs, nmpc now: ',i6)
 9040 format('.... leaving... updated nmpc: ',i6)
 9100 format(/1x,'>>>>> error: the number of nodes in each list is ',
     & 'not the same.',
     &  /14x,'1st and 2nd lists: ',2i7)
c
      end subroutine incon_mpcs_rigid
c
      end subroutine incon_mpcs
c
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine incon_mpcs_resize                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/21/2018 rhd             *
c     *                                                              *
c     *     increase size of the user_mpc_table. make larger one,    *
c     *     deep copy old -> new, move allocation. release old table *
c     *                                                              *
c     ****************************************************************
c

      subroutine incon_mpcs_resize( type, user_mpc_col_list )
      use global_data ! old common.main
      use main_data, only : modified_mpcs, force_solver_rebuild
      use mod_mpc, only : mpcs_exist, user_mpc_table, mpc_eqn, 
     &                    num_user_mpc
      implicit none
c
c               type = 1 create new larger table. copy old into new
c               type = 2 make new table but omit zeroed entries
c                        in current table
c
      integer :: old_mpc_size, i, j, nt, type, inew,
     &           user_mpc_col_list(*) 
      logical, parameter :: local_debug = .false.
      type (mpc_eqn), allocatable, dimension (:) :: new_mpc_table
c
      if( local_debug ) then
         write(*,*) '...  resizing/compressing mpc user table.....'
         write(*,*) '       old_max_mpc: ', max_mpc
      end if
c
      old_mpc_size = max_mpc
      if( type == 1 ) max_mpc = 2 * old_mpc_size
      allocate( new_mpc_table(max_mpc) )
c
c              initialize new table. 
c
      do i =  1, max_mpc
          new_mpc_table(i)%num_terms = 0
          new_mpc_table(i)%constant = 0.0
          new_mpc_table(i)%node_list => null()
          new_mpc_table(i)%dof_list  => null()
          new_mpc_table(i)%multiplier_list  => null()
      end do ! on i
c
      if( type .eq. 1 ) call incon_mpcs_resize_1
      if( type .eq. 2 ) call incon_mpcs_resize_2
c
      return
c
      contains
c     ========
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine incon_mpcs_resize_1              *
c     *                                                              *
c     *                       written by :rhd                        *
c     *                                                              *
c     *                   last modified : 12/22/2018 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine incon_mpcs_resize_1
      implicit none
c
c              deep copy as required
c
      do i = 1, old_mpc_size
          nt = user_mpc_table(i)%num_terms
          new_mpc_table(i)%num_terms = nt
          new_mpc_table(i)%constant = user_mpc_table(i)%constant
          allocate( new_mpc_table(i)%node_list(nt),
     &              new_mpc_table(i)%dof_list(nt),
     &              new_mpc_table(i)%multiplier_list(nt)  )
          do j = 1, nt
           new_mpc_table(i)%node_list(j) =
     &                   user_mpc_table(i)%node_list(j)
           new_mpc_table(i)%dof_list(j) =
     &                   user_mpc_table(i)%dof_list(j)
           new_mpc_table(i)%multiplier_list(j) =
     &                   user_mpc_table(i)%multiplier_list(j)
          end do ! on j
      end do
c
c              initialize remainder of new table. probably not req'd
c
      do i =  old_mpc_size+1, max_mpc
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
      deallocate( user_mpc_table )
      call move_alloc( new_mpc_table, user_mpc_table )
c
      return
      end subroutine incon_mpcs_resize_1

c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine incon_mpcs_resize_2              *
c     *                                                              *
c     *                       written by :rhd                        *
c     *                                                              *
c     *                   last modified : 12/22/2018 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine incon_mpcs_resize_2
      implicit none
c
c              deep copy as required
c
      inew = 0
      do i = 1,  num_user_mpc
        if( user_mpc_col_list(i) .eq. 1 ) cycle
        modified_mpcs = .true.
        inew = inew + 1
        nt = user_mpc_table(i)%num_terms
        new_mpc_table(inew)%num_terms = nt
        new_mpc_table(inew)%constant = user_mpc_table(i)%constant
        allocate( new_mpc_table(inew)%node_list(nt),
     &              new_mpc_table(inew)%dof_list(nt),
     &              new_mpc_table(inew)%multiplier_list(nt)  )
        do j = 1, nt
           new_mpc_table(inew)%node_list(j) =
     &                   user_mpc_table(i)%node_list(j)
           new_mpc_table(inew)%dof_list(j) =
     &                   user_mpc_table(i)%dof_list(j)
           new_mpc_table(inew)%multiplier_list(j) =
     &                   user_mpc_table(i)%multiplier_list(j)
        end do ! on j
      end do ! on i over prior number user mpcs
c
      num_user_mpc = inew
c
c              we could have released all user-defined MPCs in this
c              process. delete prior user_mpc. A new multipoints
c              command will cause new one to be created.
c
c              if we still have user MOCS remaing, release old table 
c              and move allocation. F2003 &
c              later does a deep release on derived types. move_alloc
c              does a deep release on new_mpc_table
c
      if( num_user_mpc == 0 ) then
         deallocate( user_mpc_table )
         new_constraints = .true.
         mpcs_exist = .false. 
      else
         deallocate( user_mpc_table )
         call move_alloc( new_mpc_table, user_mpc_table )
         modified_mpcs = .true.
         mpcs_exist = .true.
         new_constraints = .true.
      end if
c
      force_solver_rebuild = .true.
c
      return 
      end subroutine incon_mpcs_resize_2
c
      end subroutine incon_mpcs_resize
c
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine incon_mpcs_store                 *
c     *                                                              *
c     *                       written by : bjb                       *
c     *                                                              *
c     *                   last modified : 12/22/2018 rhd             *
c     *                                                              *
c     *     stores a user-defined MPC equation into global data      *
c     *     data structure                                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine incon_mpcs_store( nterm, const, node, dof, multi )
      use mod_mpc, only : num_user_mpc, user_mpc_table
      use global_data, only : max_mpc
c
      implicit none
      include 'param_def'
c
c              parameters
c
      integer :: nterm, node(*), dof(*)
      real ::    const, multi(*)
c
c              locals
c
      integer ::  err, dumi, idummy(1), i
      real ::     dumr, factor
      double precision :: dumd
      character(len=1) :: dums
c
      num_user_mpc = num_user_mpc + 1
      if( num_user_mpc .gt. max_mpc ) 
     &       call incon_mpcs_resize( 1, idummy ) 
c
      user_mpc_table(num_user_mpc)%num_terms = nterm
      user_mpc_table(num_user_mpc)%constant  = const
c
      allocate( user_mpc_table(num_user_mpc)%node_list(nterm),
     &          user_mpc_table(num_user_mpc)%dof_list(nterm),
     &          user_mpc_table(num_user_mpc)%multiplier_list(nterm),
     &          stat=err )
c
      if( err .ne. 0 ) then
         call errmsg2( 47,dumi,dums,dumr,dumd )
         call die_abort
      end if
c
c              normalize the multiplier on each dof in the
c              equations such that the leading term is -1.0
c
      factor = 1.0 / abs( multi(1) )
      if( multi(1) .gt. 0.0 ) factor = -1.0 * factor
c
      do i = 1, nterm
         user_mpc_table(num_user_mpc)%node_list(i)       = node(i)
         user_mpc_table(num_user_mpc)%dof_list(i)        = dof(i)
         user_mpc_table(num_user_mpc)%multiplier_list(i) =
     &         multi(i)*factor
      end do
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                     subroutine incon_flushline               *
c     *                                                              *
c     *                       written by : bjb                       *
c     *                                                              *
c     *                   last modified : 02/24/03                   *
c     *                                                              *
c     *     this subroutine skips a logical line of the input file   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine incon_flushline
      integer  dum
      logical  matchs, endcrd, true
c
      call readsc
      do
         if (matchs(',',1)) then
            if (endcrd(dum)) then
               call readsc
               cycle
            end if
         end if
         if (.not.matchs(',',1)) then
            if (endcrd(dum)) return
            if (true(dum)) then
               call splunj
               cycle
            end if
         end if
      end do
c
      end
