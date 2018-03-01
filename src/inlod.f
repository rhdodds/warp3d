c     ****************************************************************
c     *                                                              *
c     *                      subroutine inlod                        *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 2/8/2018 rhd               *
c     *                                                              *
c     *              translate and store loading definitions         *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine inlod( sbflg1, sbflg2, lodnum, path )
      use global_data ! old common.main
      use elem_load_data, only : numfaces
      use main_data, only : step_load_data, temp_nodmap, temp_nodlod,
     &                      node_load_defs, tables, stpchk,
     &                      user_cnstrn_stp_factors, max_step_limit,
     &                      actual_cnstrn_stp_factors,
     &                      load_data_for_a_step
      implicit none
c
c                       parameter declarations
c
      integer :: lodnum, path
      logical :: sbflg1, sbflg2
c
c                       local declarations
c
      integer :: i, dummy, dum, nc, lodn, param, open, openxt, prelod,
     &           step_count, lenlst, errnum, icn, iplist, dofn, iii,
     &           num_patterns, step, stldnm, word1, bit, test, stnod,
     &           node, idummy, column, face, dof, pist_tabn
      integer, allocatable :: intlst(:), step_load_list(:),
     &                        list_of_steps(:)
      integer, save :: nodlim, nodcol
      double precision ::
     &  forval, mpfact, dumd, body_force(mxndof),
     &  face_force(numfaces,mxndof), press_force(numfaces),
     &  elem_temper, step_load_factors(mxlc)
      double precision, parameter :: zero = 0.0d0
      character :: name*80, lname*8, stlnam*8, dums*1, curtyp*4,
     &             pname*24, workstr*80
      logical :: found, elndld, inpflg(mxndldcm),
     &     face_inpflg(numfaces,mxndof), body_inpflg(mxndof),
     &     press_set(numfaces), constraint_factor,
     &     pist_set(numfaces), abaqus_face_flag, string
      logical :: matchs, endcrd, true, label, numd, scanms,
     &           numi, integr
      logical, parameter :: debug = .false.
      real :: dumr, force(mxndldcm)
c
c                       if sub flag 1 is set the subroutine is re-
c                       entered.
c
      if( debug ) write (*,*) '>>>>>>>>> in inlod'
      allocate( intlst(mxlsz), step_load_list(mxlc),
     &          list_of_steps(max_step_limit) ) ! auto deallocated
      if( sbflg1 ) then
c
c                       re-enter nodal or element load input
c
       elndld = .true.
       go to 810
      end if
c
c
c **********************************************************************
c *                                                                    *
c *                     get the name of the loading being input.       *
c *                     each loading must have a name for id           *
c *                     purposes. if there is no name, or if the vector*
c *                     containing the load names is full, print an    *
c *                     error message and read another high level      *
c *                     command.                                       *
c *                                                                    *
c **********************************************************************
c
c
      if( label(dummy) ) then ! start of processing loading cond name
c
         name  = ' '
         lname = ' '
         call entits(name,nc)
         if( nc .gt. 8 ) nc = 8
         lname(1:nc) = name(1:nc)
c
c                       search for the specified loading in the
c                       existing loading library. if it is found,
c                       use the give load number.
c
         lodn = lodhed/two16
c
c                       the first time through, lodhed =2127298561
c                       and statement 801 is true.
c                       when the list reaches the last value, lodn =
c                       32460.  so this loops through all the load names
c                       to be sure the name hasn't already appeared.
c
 801     if(lodn.eq.32460) go to 802
            if(scanms (lodnam(lodn),lname,8)) then
               call errmsg(209,param,dums,dumr,dumd)
               lodnum = lodn
               go to 804
            end if
            lodn= lodlst(lodn)/two16
         go to 801
c
c                       not in library. check to make sure library is
c                       not full.
c
 802     numlod= numlod+1
         if(numlod.gt.mxlc) then
            numlod= mxlc
            call errmsg(53,param,dums,dumr,dumd)
            path= 1
            go to 9997
         end if
c
c                       find the first open slot in the loading lib-
c                       rary vector lodnam and assign it to the spec-
c                       ified loading.  (integer arithmetic)
c
         lodn= lodhed/two16
         open= lodhed-lodn*two16
c
c                       make sure there isn't an error in the slot
c                       available.
c
         if(open.le.0.or.open.gt.mxlc) then
            param= 2
            call errmsg(114,param,dums,dumr,dumd)
            path= 1
            go to 9997
         end if
c
         lodnum= open
         lodnam(lodnum)= lname
c
c                       update the list of open and occupied slots
c                       in lodnam. note that if a slot is occupied,
c                       it cannot also be open. therefore, if the
c                       slot is occupied, set the open linked list
c                       to zero for that slot, and vice-versa.
c
c                       find the next open slot
c
         openxt= lodlst(open)-(lodlst(open)/two16)*two16
c
c                       check for the condition that the head of the
c                       occupied list is also the end. if so, modify
c                       both lists. in either case, modify the heads
c                       of both lists accordingly.
c
         if(lodn.eq.32460) then
            lodhed= open*two16+openxt
            lodlst(open)= lodn*two16
            go to 804
         else
            lodhed= lodn*two16+openxt
         end if
c
c                       the head of the occupied linked list is not
c                       the end. find the end and modify both lists.
c
 803     prelod = lodn
         lodn   = lodlst(prelod)/two16
         if( lodn .eq. 32460 ) then
            lodlst(prelod) = open*two16
            lodlst(open)   = lodn*two16
         else
            go to 803
         end if
c
c                       initialize the nodal load arrays
c                       and the node bit map
c
 804     nodcol = 0
         nodlim = nonode/31 + 1
         if( .not. allocated( temp_nodlod ) )
     &       allocate( temp_nodlod(1:nonode,1:2) )
         temp_nodlod(1:nonode,1:2) = 0
         if ( .not. allocated( temp_nodmap ) )
     &       allocate( temp_nodmap(1:nodlim) )
         temp_nodmap(1:nodlim) = 0
c
      else
c
         call errmsg(54,dum,dums,dumr,dumd)
         path = 1
         go to 9997
c
      end if  ! scan/setting up of loading cond. name
c
c
c **********************************************************************
c *                                                                    *
c *                     look for loading commands. the command         *
c *                     nonlinear must stand alone for                 *
c *                     each loading in which it appears. that is,     *
c *                     no other command can be used in conjunction    *
c *                     with nonlinear                                 *
c *                                                                    *
c **********************************************************************
c
c
 808  continue
      elndld = .false.
      call readsc
 810  continue
c
      if( matchs('nonlinear',9) .or. matchs('dynamic',7) ) then
         if( .not. elndld ) then
            go to 820
         else
            call errmsg(55,1,dums,dumr,dumd)
            call readsc
            go to 9998
         end if
      end if
c
      if( matchs('nodal',5) ) then
         elndld = .true.
         go to 860
      end if
c
      if( matchs('element',7) ) then
         elndld = .true.
         go to 900
      end if
c
c                       there is no match with existing commands.
c                       if there is an end of card, leave the sub-
c                       routine. if not, print an error message and
c                       scan for another command at this level.
c
      if( endcrd(dum) ) then
         call errmsg(65,dum,dums,dumr,dumd)
         path = 1
         go to 9997
      else
         call errmsg(56,dum,dums,dumr,dumd)
         if( true(dum) ) go to 810
      end if
c
c
c **********************************************************************
c *                                                                    *
c *                     nonlinear : loads for time steps               *
c *                     --------------------------------               *
c *                                                                    *
c **********************************************************************
c
c
 820  continue
c
c                       initialize variables which keep track of the
c                       lowest and highest step number given for the
c                       loading.
c
      histep = 0     ! highest step
      lowstp = 1000000
c
c                       step check vector initialized by initst.f
c                       set flag indicating that this loading is a
c                       time step definition.
c
      lodtyp(lodnum) = 'TIMESTEP'
c
c                       allocate the data structure to store the list
c                       of loading pattern numbers and multipliers for
c                       each load step.
c
      if ( .not. allocated( step_load_data ) ) call mem_allocate( 19 )
c
c                       look for the command step. if it is not pre-
c                       ent, skip out to look for high level commands,
c                       returning here if can't find one.
c
 822  call readsc
c
      if( matchs('steps',4) ) then
c
c                       get a list of step numbers and expand list into
c                       local vector for easy reference. "all" is not
c                       allowed here.
c
         list_of_steps(1:max_step_limit) = 0
         step_count = 0
         call scan
         call trlist(intlst,mxlsz,0,lenlst,errnum)
c
c                       branch on the return code from trlist. a
c                       value of 1 indicates no error. a value of
c                       2 indicates that the parse rules failed in
c                       the list. a value of 3 indicates that the
c                       list overflowed its maximum length of mxlsz.
c                       a value of 4 indicates that no list was found.
c                       in these last 3 cases, the rest of the card
c                       will be ignored and a new card will be sought.
c
         if( errnum .eq. 2 ) then
            param= 1
            call errmsg(24,param,dums,dumr,dumd)
            go to 822
         else if( errnum .eq. 3 ) then
            param= 2
            call errmsg(24,param,dums,dumr,dumd)
            go to 822
         else if( errnum .eq. 4 ) then
            param=4
            call errmsg(24,param,dums,dumr,dumd)
            go to 822
         else
            if( errnum .eq. 1 ) then
               call backsp(1)
               go to 823
            end if
            param = 3
            call errmsg(24,param,dums,dumr,dumd)
            go to 822
         end if
c
c                       parse the list of steps, initialize the loading
c                       patterns defined for each step
c
 823     icn    = 0
         iplist = 1
 824     call trxlst(intlst,lenlst,iplist,icn,step)
c
         if( step .gt. max_step_limit )  call inlod_resize_steps
c
         if( step .lt. 0 ) then
            param = step
            call errmsg(69,param,dums,dumr,dumd)
            go to 826
         end if
c
         step_count                = step_count + 1
         list_of_steps(step_count) = step
         if( step .gt. histep ) histep = step
         if( step .lt. lowstp ) lowstp = step
c
c                       initialize the load factors for this step
c                       for the given loading to zero. the default
c                       factor for incremental constraints for each
c                       step is 1.0
c
         user_cnstrn_stp_factors(step) = 1.0
         if ( step_load_data(step)%num_load_patterns .gt. 0 ) then
          deallocate( step_load_data(step)%load_patt_num )
          deallocate( step_load_data(step)%load_patt_factor )
         end if
c
c                       set the step check vector for this step to
c                       true, signifying that this step has been
c                       defined (at least the incremental constraints).
c
         stpchk(step) = .true.
c
 826     if( iplist .ne. 0 ) go to 824
c
c                       list of steps traversed, checked and existing
c                       list of loading paterns destroyed
c
         if ( step_count .eq. 0 ) go to 822
         step_load_list(1:mxlc) = 0
         step_load_factors(1:mxlc) = zero
         num_patterns = 0
c
c                       input the loadings and their weights associated
c                       with the steps in the list.
c
 827     continue
         if( label(dummy) ) go to 828
c
c                       if there is a comma, check for continuation.
c                       if there is an end of card, save the
c                       definitions for the current step list and
c                       return for another step list. otherwise,
c                       skip the current entity and try to input
c                       a loading name.
c
         if( endcrd(dummy) ) then
            call save_step_definitions( step_load_list, list_of_steps,
     &                                  step_load_factors, num_patterns,
     &                                  step_count )
            go to 822
         else if( matchs(',',1) ) then
            if( endcrd(dummy) ) call readsc
            go to 827
         else
            call errmsg(59,dum,dums,dumr,dumd)
            if( true(dummy) ) go to 827
         end if
c
 828     continue
c
c                       if the loading name input cannot be found,
c                       then skip it and look for a new loading name.
c
         found   = .false.
         name    = ' '
         stlnam  = ' '
         call entits(name,nc)
         if( nc .gt. 8 ) nc = 8
         stlnam(1:nc)      = name(1:nc)
         constraint_factor = .false.
         if( scanms('constrai', stlnam, 8) ) then
           constraint_factor = .true.
           found             = .true.
           go to 830
         end if
c
c                       search the loading library for the specified
c                       loading pattern used in step definition.
c
         lodn = lodhed/two16
 829     if( lodn.eq.32460 ) go to 830
            if( scanms(lodnam(lodn),stlnam,8) ) then
               if( .not. scanms(lodtyp(lodn),'REGULAR ',8) ) then
                  param = lodn
                  call errmsg(88,param,dums,dumr,dumd)
                  go to 827
               end if
c
               stldnm = lodn
               found  = .true.
               num_patterns = num_patterns + 1
               go to 830
            end if
c
            lodn = lodlst(lodn)/two16
c
         go to 829
c
 830     if( .not. found ) then
           call errmsg(60,dum,dums,dumr,dumd)
           go to 827
         end if
c
c                       get the multiplication factor for the loading
c                       pattern (or constraints)
c
         if( numd(mpfact) ) then
            if ( constraint_factor ) then
              do i = 1, step_count
                user_cnstrn_stp_factors(list_of_steps(i)) = mpfact
              end do
            else
              step_load_list(num_patterns)    = lodn
              step_load_factors(num_patterns) = mpfact
            end if
         else
            call errmsg(67,dum,dums,dumr,dumd)
         end if
c
c                       return to process more loadings
c
         go to 827
c
      else       !   end of the if ( matchs( 'step'...
c
         path = 0
         go to 9996
c
      end if
c
c **********************************************************************
c *                                                                    *
c *                     nodal forces and temperatures                  *
c *                     -----------------------------                  *
c *                                                                    *
c **********************************************************************
c
c
 860  continue
c
c                       set flag indicating that this loading is a
c                       regular definition of actual loads.
c
c                       how_defined = 0 lists of nodes and values
c                                       in text input (default)
c                                   = 1 a user nodal loads routine
c                                       to be called during solution.
c                                       store a user load file name if
c                                       one given
c
      lodtyp(lodnum)= 'REGULAR '
      node_load_defs(lodnum)%how_defined = 0
c
      if( matchs('loads',4) ) call splunj
 861  call readsc  ! loop to get line of node numbers & forces, temp

      if( matchs('user_routine',4) ) then
          node_load_defs(lodnum)%how_defined = 1
          if( matchs('file',4) ) call splunj
          workstr(1:80) = " "
          if( label(dummy) ) then
               call entits( workstr, nc )
          elseif( string(dummy) ) then
               call entits( workstr, nc )
          end if
          node_load_defs(lodnum)%user_file_name(1:) =  workstr
          call readsc
      else
          call reset
          if( true(idummy) ) call splunj
      end if
c
c                       translate list of nodes upon which loads occur.
c
      call scan
      call trlist(intlst,mxlsz,nonode,lenlst,errnum)
c
c                       branch on the return code from trlist. a
c                       value of 1 indicates line start with an
c                       <integerlist>  -- start processing.

c                       a value of 2 indicates that the parse rules
c                       failed in the list.

c                       a value of 3 indicates that the
c                       list overflowed its maximum length of mxlsz.
c                       in these last two cases, the rest of the line
c                       will be ignored and a new card will be sought.

c                       a value of 4 indicates that no list was found.
c                       in this case, either node loads input has ceased
c                       or the an error has been made by omitting the
c                       list. save the node loading in permanent data
c                       structure.
c
      if( errnum .eq. 2 ) then
         param = 1
         call errmsg(24,param,dums,dumr,dumd)
         go to 861
      else if( errnum .eq. 3 ) then
         param = 2
         call errmsg(24,param,dums,dumr,dumd)
         go to 861
      else if( errnum .eq. 4 ) then
         call backsp(1)
         if ( node_load_defs(lodnum)%node_count .ne. 0 ) then
           deallocate( node_load_defs(lodnum)%nodal_loads )
           nullify( node_load_defs(lodnum)%nodal_loads )
         end if
         if( nodcol .gt. 0 ) then ! zero for user routine, file, etc.
           allocate( node_load_defs(lodnum)%nodal_loads(1:nodcol,1:2) )
           node_load_defs(lodnum)%nodal_loads(1:nodcol,1:2) =
     &        temp_nodlod(1:nodcol,1:2)
           node_load_defs(lodnum)%node_count = nodcol
         end if
         go to 9995
      else   ! line starts with <integerlist>
         if( errnum .eq. 1 ) then
            call backsp(1)
            go to 862
         end if
         param = 3
         call errmsg(24,param,dums,dumr,dumd)
         go to 861
      end if
c
c                       there is a valid node list. input the loads for
c                       the nodes in the list. initialize the temporary
c                       storage vectors. note: a temperature is placed
c                       in the last row of the force vector and
c                       in the loads table.
c
 862  do i = 1, mxndldcm
         force(i)  = zero
         inpflg(i) = .false.
      end do
      go to 865
c
 864  call readsc   ! there was a comma at end of previous line
c
c                       branch on load type input.
c
 865  if( matchs('force_x',7) ) then
         dofn   = 1
         curtyp = 'forx'
         go to 867
      end if
      if( matchs('force_y',7) ) then
         dofn   = 2
         curtyp = 'fory'
         go to 867
      end if
      if( matchs('force_z',7) ) then
         dofn   = 3
         curtyp = 'forz'
         go to 867
      end if
      if( matchs('temperature',4) ) then
         dofn   = mxndof + 1
         curtyp = 'temp'
         go to 867
      end if
c
c                       match for node keyword failed. if end line,
c                       ok to start storing values. if comma, loop
c                       back to get new line and keep processing node
c                       keywords. otherwise ignore bad keyword,
c                       try to skip over and keep going.
c
      if( matchs(',',1) ) go to 866
      if( endcrd(dum) ) then
         go to 873
      else
         call errmsg(66,dum,dums,dumr,dumd)
         if( true(dum) ) go to 865
      end if
c
c                       if there is a comma at the end of a line, the
c                       input line is continued.
c
 866  continue
      if( endcrd(dum) ) then
         go to 864
      end if
      go to 865
c
c                       get the load value following keyword
c
 867  if( matchs('=',1) ) call splunj
      if( .not. numd(forval) ) then
         call errmsg(57,dum,dums,dumr,dumd)
      else
c
c                       make sure that the node force value has
c                       not already been input on the same command
c
         if( inpflg(dofn) ) then
            call errmsg(137,dum,curtyp,dumr,dumd)
            go to 865
         end if
         inpflg(dofn) = .true.
         force(dofn)  = forval
c
      end if
      go to 865
c
c                       store the vector of load magnitudes in the
c                       load data array and set the pointer, etc, in
c                       the nodal load array for each node in the
c                       current list.
c
 873  icn    = 0
      iplist = 1
c
 874  call trxlst(intlst,lenlst,iplist,icn,stnod)
c
c                       check that the list node is ok
c
      if( stnod .gt. nonode ) then
         param = stnod
         call errmsg(16,param,dums,dumr,dumd)
         go to 882
      end if
      if( stnod .lt .0 ) then
         param = stnod
         call errmsg(58,param,dums,dumr,dumd)
         go to 882
      end if
c
c                       check to see if the node has already been
c                       mentioned during input for this loading
c                       condition
c
      word1 = (stnod-1)/31 + 1
      bit   = stnod-31*(word1-1)
      test  = 0
      test  = iand(temp_nodmap(word1),bits(bit))
c
 875  if( test .eq. 0 ) then
c
c                       bit is off. the node has not been previously
c                       accessed. loddat_blocks data strucutre contains
c                       input nodal forces and temperatures for all
c                       loading conditions. the column for a specific
c                       loading condition & node number is contained in
c                       the node_load_defs dtaa structure.
c
         numcol = numcol + 1  ! this is a major global variable
c                               in common.main that should be better
c                               protected. it sets the last used column
c                               in the loddat_blocks structure
c
         call loddat_ops( 1, force, iii )
c
c                       set the index to the load magnitudes and the
c                       node number in the nodal load array. nodcol
c                       is initialized to zero above when the loading
c                       condition is defined. turn on the bit for the
c                       node in the temporary bit map list so we
c                       can know fast if the same node is mentioned
c                       again during the nodal load input.
c
         nodcol = nodcol + 1
         temp_nodlod(nodcol,1) = stnod
         temp_nodlod(nodcol,2) = numcol
         temp_nodmap(word1)    = ior(temp_nodmap(word1),bits(bit))
c
      else
c
c                       bit for the specified node is already on.
c                       a force/temperature has already been
c                       specified for this node in this loading.
c                       find which column in loddat has the node's
c                       load values.
c
         do i = 1, nodcol
            node  = temp_nodlod(i,1)
            found = .false.
            if( node .eq. stnod ) then
               column = temp_nodlod(i,2)
               found  = .true.
               go to 880
            end if
         end do
c
c                       make sure that the node was actually found.
c                       if not, then branch to poriton of subroutine
c                       where new columns in loddat are added for the
c                       nodal load option.
c
 880           if( .not. found ) then
                  test = 0
                  go to 875
               end if
c
c                       add the vector of load magnitudes to the
c                       existing data in the load data array.
c
         call loddat_ops( 2, force, column )
c
      end if
c
 882  if( iplist .ne. 0 ) go to 874   ! more nodes in current
c                                       <integerlist>
c
c                       return to possibly process more load input
c
      go to 861  ! read a new line and look for another node list
c
c **********************************************************************
c *                                                                    *
c *                     element loads                                  *
c *                     -------------                                  *
c *                                                                    *
c * body force, face tractions, face pressure, constant temperature,   *
c * and piston loads for aerodynamic pressure on deforming surfaces    *
c *                                                                    *
c **********************************************************************
c
c
 900  continue
c
c                       set flag indicating that this loading is a
c                       regular definition of actual loads.
c
      lodtyp(lodnum)= 'REGULAR '
      if(matchs('loads',4)) call splunj
 905  call readsc
c
c                       translate list of elements upon which loads occur.
c
      call scan
      call trlist(intlst,mxlsz,noelem,lenlst,errnum)
c
c                       branch on the return code from trlist. a
c                       value of 1 indicates no error. a value of
c                       2 indicates that the parse rules failed in
c                       the list. a value of 3 indicates that the
c                       list overflowed its maximum length of mxlsz.
c                       in these last two cases, the rest of the card
c                       will be ignored and a new card will be sought.
c                       a value of 4 indicates that no list was found.
c                       in this case, either element loads input has ceased
c                       or the an error has been made by omitting the
c                       list.
c
      if(errnum.eq.2) then
         param= 1
         call errmsg(24,param,dums,dumr,dumd)
         go to 905
      else if(errnum.eq.3) then
         param= 2
         call errmsg(24,param,dums,dumr,dumd)
         go to 905
      else if(errnum.eq.4) then
         call backsp(1)
         go to 9995
      else
         if(errnum.eq.1) then
            call backsp(1)
            go to 907
         end if
         param= 3
         call errmsg(24,param,dums,dumr,dumd)
         go to 905
      end if
c
c                       there is a valid list. input the loads for the
c                       elements in the list.
c                       initialize the temporary storage vectors.
c
 907  continue
      do dof = 1, mxndof
         body_force(dof)  = zero
         body_inpflg(dof) = .false.
         do face = 1, numfaces
            face_force(face,dof) = zero
            face_inpflg(face,dof) = .false.
         end do
      end do
      do face = 1, numfaces
         press_force(face) = zero
         press_set(face)   = .false.
         pist_set(face)    = .false.
      end do
      elem_temper = zero
      go to 920
c
 910  call readsc   ! top of loop to keep reading element loads
c
c                       branch on load type input.  If end of line,
c                       store the loadings in a temporary structure.
c
 920  continue
      if ( matchs('body',4) ) go to 1000
      if ( matchs('face',4) ) go to 1100
      if ( matchs('temperature',6) ) go to 1200
      if ( endcrd(dum) ) then
         call store_loadings(body_inpflg,face_inpflg,press_set,
     &        body_force,face_force, press_force, intlst, lenlst,
     &        elem_temper, pist_tabn, pist_set )
         go to 905
      else
         call errmsg(66,dum,dums,dumr,dumd)
       call scan
         go to 920
      end if
c
c                     ----------------
c                     | body loading |
c                     ----------------
 1000 continue
      if (debug) write (*,*) 'body loading'
      if (matchs('forces',5)) call splunj
c
c                       branch on body loading input
c
 1010 continue
c
      if(matchs('bx',2)) then
         dofn= 1
         curtyp= 'bx'
         go to 1020
      end if
c
      if(matchs('by',2)) then
         dofn= 2
         curtyp= 'by'
         go to 1020
      end if
c
      if(matchs('bz',2)) then
         dofn= 3
         curtyp= 'bz'
         go to 1020
      end if
c
      goto 920
c
c                       store the body forces in the temporary
c                       force vector.
c
 1020 if( matchs('=',1) ) call splunj
      if( .not.numd(forval) ) then
         call errmsg(57,dum,dums,dumr,dumd)
      else
c
c                       make sure that the body force for the current
c                       dof type has not already been input on the same
c                       initial element loading command line.
c
         if(body_inpflg(dofn)) then
            call errmsg(137,dum,curtyp,dumr,dumd)
            go to 1010
         end if
c
         body_inpflg(dofn)= .true.
         body_force(dofn)= forval
c
      end if
      go to 1010
c
c                     ----------------
c                     | face loading |
c                     ----------------
c
 1100 continue
c
c                      get face number and check it
c
      if (debug) write (*,*) 'face loading'
      abaqus_face_flag = .false.
      if( matchs('abaqus',4) ) abaqus_face_flag = .true.
      if (.not. integr(face)) then
         call errmsg(226,dum,dums,dumr,dumd)
         goto 905
      endif
c
      if (face.lt.0 .or. face.gt.numfaces) then
         call errmsg(227,face,dums,dumr,dumd)
         goto 905
      endif
c
c                      the abaqus flag applies only to hex
c                      elements, faces (4,5) in WARP3D are (5,4)
c                      in abaqus. this option is convenient to allow
c                      input of abaqus face numbers by user. should
c                      not be used for tet elements
c
      if( abaqus_face_flag ) then
        if( face .eq. 4 ) then
           face = 5
        elseif( face .eq. 5 ) then
           face = 4
        endif
       endif
       abaqus_face_flag = .false.
c
c                      get pressure load if specified.
c
      if (matchs('pressure',5)) then
         if( matchs('=',1) ) call splunj
         if (.not. numd(forval)) then
            call errmsg(57,dum,dums,dumr,dumd)
         else
c
c                       make sure that the pressure for the current
c                       face has not already been input on the same
c                       initial element loading command line.
c
            if (press_set(face)) then
               call errmsg(228,face,dums,dumr,dumd)
               goto 920
            endif
c
            press_force(face) = forval
            press_set(face) = .true.
         endif
         goto 920
      endif
c
c                      get piston load on face if specified
c
      if (matchs('piston',6)) then
c
c                       initiate loading flags and vectors
c
         if (debug) write (*,*) 'piston loading'
         pist_set(face) = .false.
c
c                       check that the next entry on the line
c                       is a string.
c
         if ( .not. label(dummy) ) then
c
            call errmsg(54,dum,dums,dumr,dumd)
            go to 920
c
         end if
c
         name  = ' '
         pname = ' '
         call entits(name,nc)
         if( nc .gt. 24 ) nc = 24
         pname(1:nc) = name(1:nc)
         if ( debug ) write(*,*)
     &        'Piston loading name = ', pname
c
c                       look for name of previously defnined
c                       table of type 'piston'. if table has
c                       not been defined previously, this is
c                       a fatal error
c
         pist_tabn = -1
         do i = 1,max_tables
            if ( scanms(tables(i)%table_name,pname,24) ) then
               if ( scanms(tables(i)%table_type,'PISTON  ',8) ) then
                  pist_tabn = i
                  if ( debug ) write(*,*)
     &                 'Piston loading # =', pist_tabn
                  pist_set(face) = .true.
               end if
            end if
         end do
c
         if ( pist_tabn .eq. -1 )
     &        call errmsg(334,dum,dums,dumr,dumd)
         go to 920
c
      end if      ! end of piston loading definition
c
c                      get traction values
c
      if (matchs('force',5)) call splunj
      if (matchs('tractions',5)) call splunj
      if (matchs('force',5)) call splunj
c
c                      parse through traction specifications
c
 1110 continue
      if(matchs('tx',2)) then
         dofn= 1
         curtyp= 'tx'
         go to 1120
      end if
      if(matchs('ty',2)) then
         dofn= 2
         curtyp= 'ty'
         go to 1120
      end if
      if(matchs('tz',2)) then
         dofn= 3
         curtyp= 'tz'
         go to 1120
      end if
      go to 920
c
c                       store the face loadings in the temporary
c                       force vector.
c
 1120 continue
      if( matchs('=',1) ) call splunj
      if( .not.numd(forval) ) then
         call errmsg(57,dum,dums,dumr,dumd)
      else
c
c                       make sure that the current face force has
c                       not already been input on the same initial
c                       element loading command line.
c
         if( face_inpflg(face,dofn) ) then
            call errmsg(137,dum,curtyp,dumr,dumd)
            go to 1110
         end if
         face_inpflg(face,dofn) = .true.
         face_force(face,dofn)  = forval
      end if
      go to 1110
c
c                     -------------------------------
c                     | element temperature loading |
c                     -------------------------------
 1200 continue
      if( debug ) write (*,*) 'element temperature loading'
      if( matchs('=',1) ) call splunj
      if( .not.numd(elem_temper) ) then
         call errmsg(57,dum,dums,dumr,dumd)
      end if
      go to 920
c
c
c
c **********************************************************************
c **********************************************************************
c
c
 9995 sbflg1= .true.
      sbflg2= .true.
      if (debug) write (*,*) '9995'
      go to 9999
c
 9996 sbflg1= .false.
      sbflg2= .true.
      if (debug) write (*,*) '9996'
      go to 9999
c
 9997 sbflg1= .false.
      sbflg2= .false.
      if (debug) write (*,*) '9997'
      go to 9999
c
 9998 sbflg1= .true.
      sbflg2= .false.
      if (debug) write (*,*) '9998'
c
c
 9999 continue
      return
c
      contains
c     ========
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine inlod_resize_steps               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 2/8/2018 rhd               *
c     *                                                              *
c     *     resize global data structures that store definition for  *
c     *     each load step to increase the allowable number of steps *
c     *                                                              *
c     ****************************************************************

      subroutine inlod_resize_steps
      implicit none
c
      integer :: old_max_steps, i, alloc_stat, np
      integer, allocatable, dimension(:) :: new_list
      real, allocatable, dimension(:) :: new_user, new_actual
      logical, allocatable, dimension(:) :: new_stpchk
      type(load_data_for_a_step), allocatable,
     &                            dimension(:) :: new_step_data
c
c              double the maxim number of load steps.
c
      old_max_steps = max_step_limit
      max_step_limit = max_step_limit * 2
c
c              resize the working step list for inload driver
c
      allocate( new_list(max_step_limit) )
      new_list(1:old_max_steps) = list_of_steps(1:old_max_steps)
      call move_alloc( new_list, list_of_steps )
c
c              resize simple global vectors. move_alloc deallocates
c              as needed.
c
      allocate( new_user(max_step_limit), new_actual(max_step_limit),
     &          new_stpchk(max_step_limit) )
c
      do i = 1, old_max_steps
        new_user(i)   = user_cnstrn_stp_factors(i)
        new_actual(i) = actual_cnstrn_stp_factors(i)
        new_stpchk(i) = stpchk(i)
      end do
      call move_alloc( new_user, user_cnstrn_stp_factors )
      call move_alloc( new_actual, actual_cnstrn_stp_factors )
      call move_alloc( new_stpchk, stpchk )
c
c              for step_load_data we have to do the deep copy of
c              allocated vectors of the derived type. move_alloc
c              does not do this.
c
      allocate( new_step_data(max_step_limit), stat = alloc_stat )
c
      do i = old_max_steps + 1, max_step_limit
        new_step_data(i)%num_load_patterns  = 0
      end do
c
      do i = 1, old_max_steps
       np = step_load_data(i)%num_load_patterns
       if( np == 0 ) cycle
       new_step_data(i)%num_load_patterns = np
       allocate( new_step_data(i)%load_patt_num(np) )
       new_step_data(i)%load_patt_num(1:np) =
     &         step_load_data(i)%load_patt_num(1:np)
       deallocate( step_load_data(i)%load_patt_num )
       allocate( new_step_data(i)%load_patt_factor(np) )
       new_step_data(i)%load_patt_factor(1:np) =
     &         step_load_data(i)%load_patt_factor(1:np)
       deallocate( step_load_data(i)%load_patt_factor )
      end do
c
      call move_alloc( new_step_data, step_load_data )
c
      return

      end subroutine inlod_resize_steps
      end subroutine inlod
c     ****************************************************************
c     *                                                              *
c     *                      subroutine store_loadings               *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 11/15/11 jcs               *
c     *                                                              *
c     *     this subroutine stores the face and body loading         *
c     *     information into the temporary storage arrays.           *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine store_loadings( body_inpflg, face_inpflg, press_set,
     &     body_force, face_force, press_force, intlst, lenlst,
     &     elem_temper, pist_tabn, pist_set )
      use global_data ! old common.main
      use elem_load_data
      implicit none
c
      integer :: intlst(*), lenlst, pist_tabn
      double precision ::
     &  body_force(*), face_force(numfaces,*),
     &  press_force(*), dumd, elem_temper
      logical :: face_inpflg(numfaces,*), body_inpflg(*),
     &           press_set(*), pist_set(*)
c
      integer :: icn, iplist, dof, face, stelem
      logical, parameter :: debug = .false.
      double precision, parameter :: one = 1.0d0
      real :: dumr
      character(len=1) :: dums
c
      if( debug ) write (*,*) ' >>> inside store_loadings'
c
c                 allocate the temporary vectors if it hasn't already
c                 been done
c
      call allocate_temp_load(1)
c
c                       store the vector of body force and face loading
c                       magnitudes in temporary force arrays.  Check all
c                       values for errors.
c
 10   icn    = 0
      iplist = 1
c
 20   continue
      call trxlst(intlst,lenlst,iplist,icn,stelem)
c
c                       check that the list element does not exceed
c                       the number of elements in the structure.
c
      if( stelem.gt.noelem ) then
         call errmsg(35,stelem,dums,dumr,dumd)
         go to 1000
      end if
c
c                       check that the list node is not negative.
c
      if( stelem.lt.0 ) then
         call errmsg(86,stelem,dums,dumr,dumd)
         go to 1000
      end if
c
c                       now branch on the loading types
c                       1) body forces
c
      do dof = 1, mxndof
       if ( body_inpflg(dof) ) temp_body(stelem,dof) = body_force(dof)
      end do
c
c                       2) face forces
c
      do face = 1, numfaces
         if ( press_set(face) ) then
            temp_press(stelem,face) = press_force(face)
         elseif ( pist_set(face) ) then
            temp_piston(stelem,face) = pist_tabn
         else
            do dof = 1, mxndof
               if ( face_inpflg(face,dof) ) then
                  temp_face(stelem,face,dof) = face_force(face,dof)
               end if
            end do
         end if
      end do
c
c                       3) element temperature change
c
      temp_temper(stelem) = elem_temper
c
c
 1000 if( iplist.ne.0 ) go to 20
c
c
      if (debug) write (*,*) ' <<< leaving store_loadings'
      return
      end


c     ****************************************************************
c     *                                                              *
c     *               subroutine save_step_definitions               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified :  1/16/2018 rhd             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine save_step_definitions( step_load_list, list_of_steps,
     &                                  step_load_factors, num_patterns,
     &                                  step_count )
      use global_data ! old common.main
      use main_data, only : step_load_data
      implicit none
c
      integer :: step_load_list(*), list_of_steps(*),
     &           num_patterns, step_count
      double precision :: step_load_factors(*)

      integer :: i, step
      logical, parameter :: local_debug = .false.
c
      if( local_debug ) write(*,9100) num_patterns
      if( num_patterns .eq. 0 ) return
      if( local_debug ) then
        write(*,9200)
        do i = 1,  num_patterns
         write(*,9300) i,  step_load_list(i),  step_load_factors(i)
        end do
        write(*,9400) list_of_steps(1:step_count)
      end if
c
      do i = 1, step_count
       step = list_of_steps(i)
       if( allocated( step_load_data(step)%load_patt_num ) .or.
     &     allocated( step_load_data(step)%load_patt_factor) ) then
         write(*,9000) step
         call die_abort
       end if
c
       allocate( step_load_data(step)%load_patt_num(num_patterns) )
       allocate( step_load_data(step)%load_patt_factor(num_patterns) )
c
       step_load_data(step)%num_load_patterns = num_patterns
       step_load_data(step)%load_patt_num(1:num_patterns) =
     &         step_load_list(1:num_patterns)
       step_load_data(step)%load_patt_factor(1:num_patterns) =
     &         step_load_factors(1:num_patterns)
      end do
c
      return
 9000 format('>> FATAL ERROR: routine save_step_definitions.',
     &   /,  '                definitions already present for ',
     &   /,  '                load step: ',i7,
     &   /,  '                job aborted....')
 9100 format('>> inside save_step_definitions.  num_patterns: ',i4)
 9200 format(' load patt. no.  load pattern id    load factor')
 9300 format(i8,5x,i8,3x,f10.5)
 9400 format(' for load steps: ',100(10i8,/))
      end

c     ****************************************************************
c     *                                                              *
c     *                   subroutine loddat_ops                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/16/2018 rhd              *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine loddat_ops( opcode, node_data, loddat_col )
      use global_data ! old common.main
      use main_data, only : num_loddat_blks, sizeof_loddat_blks,
     &                      next_loddat_col, max_loddat_blks,
     &                      loddat_blocks
      implicit none
c
      integer :: opcode, loddat_col
      real :: node_data(*)
c
      integer :: now_block, relcol, alloc_stat
c
      select case ( opcode )

      case ( 1 )
c
c               store node values in next available column of the blocked
c               loddat data array. open up new block if needed
c

      next_loddat_col = next_loddat_col + 1
      loddat_col      = next_loddat_col
c
      now_block  = ( loddat_col - 1 ) / sizeof_loddat_blks + 1
      if ( now_block .le. num_loddat_blks ) then
        relcol = loddat_col - (now_block - 1) * sizeof_loddat_blks
        loddat_blocks(now_block)%block(1:mxndldcm,relcol) =
     &        node_data(1:mxndldcm)
        return
      end if

      num_loddat_blks = num_loddat_blks + 1
      now_block = num_loddat_blks
      if ( now_block .gt. max_loddat_blks ) then
         write(*,9100) 1
         call die_abort
         stop
      end if
      allocate( loddat_blocks(now_block)%block(1:mxndldcm,
     &            1:sizeof_loddat_blks), stat = alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(*,9100) 2
           call die_abort
           stop
      end if
      loddat_blocks(now_block)%
     &         block(1:mxndldcm,1:sizeof_loddat_blks) = 0.0
      loddat_blocks(now_block)%block(1:mxndldcm,1) =
     &        node_data(1:mxndldcm)
      return


      case ( 2 )
c
c               add node data to an already existing column of the
c               blocked array
c
c
      now_block  = ( loddat_col - 1 ) / sizeof_loddat_blks + 1
      if ( now_block .ge. 1 .and. now_block .le. num_loddat_blks ) then
        relcol = loddat_col - (now_block - 1) * sizeof_loddat_blks
        loddat_blocks(now_block)%block(1:mxndldcm,relcol) =
     &      loddat_blocks(now_block)%block(1:mxndldcm,relcol) +
     &        node_data(1:mxndldcm)
        return
      end if
      write(*,9100) 4
      call die_abort
      stop

      case ( 3 )
c
c               return the nodal loads vector for requested column
c               of loads table
c
      now_block  = ( loddat_col - 1 ) / sizeof_loddat_blks + 1
      if ( now_block .ge. 1 .and. now_block .le. num_loddat_blks ) then
        relcol = loddat_col - (now_block - 1) * sizeof_loddat_blks
        node_data(1:mxndldcm) =
     &   loddat_blocks(now_block)%block(1:mxndldcm,relcol)
        return
      end if
      write(*,9100) 5
      call die_abort
      stop

      case default
           write(*,9100) 3
           call die_abort
           stop

      end select
c
 9100 format('>>>> FATAL ERROR: routine loddat_ops @ : ',i2,
     & /,    '                  job terminated.')
      return
      end
