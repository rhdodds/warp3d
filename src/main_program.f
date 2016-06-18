c     ****************************************************************
c     *                                                              *
c     *                      Main program for WARP3D                 *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 5/22/2016 rhd              *
c     *                                                              *
c     *                      main program for WARP3D                 *
c     *                                                              *
c     ****************************************************************
c
c
c
      program warp3d
      use file_info
      use main_data, only : output_packets
      use performance_data, only : time_assembly, assembly_total,
     &            ntimes_assembly
      implicit integer (a-z)
$add common.main
      real t1, wcputime, dumr
      external wcputime
      character*8 stcnam, dums,
     &     sdate_*24
      character*80 name, stflnm, rtflnm
      logical hilcmd,sbflg1,sbflg2
      logical endcrd,label,matchs,debug1,debug2,debug,endfil,
     &        string, matchs_exact
c
      common/errprm/ erprmd(10),erprmr(10),erprmi(10),erprms
#dbl      double precision
#sgl      real
     &  erprmd, dumd
      real erprmr
      character erprms *50
      real :: wall, real_start, real_end, real_rate
c
      common/erflgs/ numnod,numel,fatal,coor,elprop,elinc,constr,block
      logical fatal,coor,numnod,numel,elprop,elinc,constr,block
      integer return_type
c
c                       MPI: initialize all processors
c                       and make workers go into the worker handler. if
c                       else   wmpi_init is a dummy routine.
c
      call wmpi_init
c
c                       initialize the load step timing and debug
c
      debug1   = .false.
      debug2   = .false.
      call setstarttime
      t1 = wcputime ( 0 )
      call steptime ( dummy, 1 )
c
c                       initialize scan and variables necessary
c                       for program execution
c
      call initst(sbflg1,sbflg2)
c
c                       print warp3d header
c
#lnx      call fdate (sdate_)
#l64      call fdate (sdate_)
#mac      call fdate (sdate_)
#win      call fdate (sdate_)
c

      write (*,9000) sdate_ , mxnod, mxel
c
c                       read a high level command
c
c
      nsn = 0
 10   continue
      call readsc
c
c                       branch on high level command
c
 20   continue
      hilcmd= .false.
c
      if(matchs('structure',9)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 1
         go to 25
      end if
c
      if(matchs('material',8)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 2
         go to 25
      end if
c
      if(matchs('crystal',7)) then
         lsn = nsn
         hilcmd = .true.
         nsn = 32
         go to 25
      end if
c
      if(matchs('number',3)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 3
         go to 25
      end if
c
      if(matchs('coordinates',5)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 4
         go to 25
      end if
c
      if(matchs('elements',8)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 5
         go to 25
      end if
c
      if(matchs('incidences',5)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 6
         go to 25
      end if
c
      if(matchs('constraints',11)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 7
         go to 25
      end if
c
      if(matchs('release',7)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 33
         go to 25
      end if
c
      if(matchs('loading',7)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 8
         go to 25
      end if
c
      if( matchs('dynamic',7) .or. matchs('nonlinear',7) ) then
         lsn= nsn
         hilcmd= .true.
         nsn= 10
         go to 25
      end if
c
      if(matchs('solution',4)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 10
         go to 25
      end if
c
      if(matchs('initial',4)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 11
         go to 25
      end if
c
      if(matchs('compute',4)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 12
         go to 25
      end if
c
      if(matchs('output',6)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 13
         go to 25
      end if
c
      if(matchs('convert',4)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 14
         go to 25
      end if
c
      if(matchs('save',4)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 15
         go to 25
      end if
c
      if(matchs('retrieve',4)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 16
         go to 25
      end if
c
      if(matchs('restart',4)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 16
         go to 25
      end if
c
      if(matchs('blocking',5)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 18
         go to 25
      end if
c
      if(matchs('stop',4) .or. matchs('exit',4) .or. matchs('quit',4))
     &      then
         lsn= nsn
         hilcmd= .true.
         nsn= 19
         go to 25
      end if
c
      if(matchs('domain',5)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 20
         go to 25
      endif
c
      if(matchs('crack',5)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 21
         go to 25
      endif
c
      if(matchs('debug',5)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 22
         go to 25
      endif
c
      if(matchs('stress',6)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 23
         go to 25
      endif
c
      if(matchs('license',4)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 24
         go to 25
      endif
c
      if(matchs('contact',4)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 25
         go to 25
      endif
c
      if(matchs('thermal',4)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 26
         go to 25
      endif
c
      if(matchs('functionally',4)) then
         lsn    = nsn
         hilcmd = .true.
         nsn    = 27
         go to 25
      endif
c
      if(matchs('surface',7)) then
         lsn    = nsn
         hilcmd = .true.
         nsn    = 28
         go to 25
      endif
c
      if(matchs('tie',3)) then
         if(matchs('mesh',4)) then
            lsn    = nsn
            hilcmd = .true.
            nsn    = 29
            go to 25
         end if
      endif
c
      if(matchs('table',3)) then
         lsn    = nsn
         hilcmd = .true.
         nsn    = 30
         go to 25
      endif
c
      if(matchs('list',4)) then
         lsn= nsn
         hilcmd= .true.
         nsn= 31
         go to 25
      end if

c
      if( endfil(dum) ) goto 30
c
c                       if a high level command is encountered,
c                       first check to see if the previous sub-
c                       routine executed requires error checking.
c                       if there  are no fatal errors found, then
c                       clear the subroutine flags for the next
c                       subroutine to be accessed.
 25   continue
      if ( hilcmd ) then
         if ( sbflg2 ) call errchk(lsn,chkprm,debug1)
         if ( fatal ) go to 9999
         sbflg1 = .false.
         sbflg2 = .false.
      else
c
c                       check if a subroutine is to be reentered or
c                       there is an error.
c
         if ( .not. sbflg1 ) nsn=0
      end if
c

c
c                       branch on the new subroutine number to be
c                       accessed.
c
      go to (100,200,300,400,500,600,700,800,900,1000,1100,
     &       1200,1300,1400,1500,1600,1700,1800,1900,2000,
     &       2100,2200,2300,2400,2500,2600,2700,2800,2900,
     &       3000,3100,3200,3300), nsn
c
c                       if a high level command is not
c                       encountered, print an error message
c                       and scan another line in search for one
c
      if (endcrd(dum)) goto 10
c
      call errmsg(1,dum,dums,dumr,dumd)
      go to 10
c
c                       reached an end on input file. decrease
c                       file count by one -- set input device.
c
 30   continue
      close (inlun(filcnt))
      filcnt=filcnt-1
      if (filcnt.eq.0) then
         call errmsg(204,dum,dums,dumr,dumd)
         goto 9999
      end if
      call setin(inlun(filcnt))
      in = inlun(filcnt)
      call errmsg(179,dum,dums,dumr, dumd)
      go to 10
c
c                       structure command. the name of the structure
c                       is input here. the structure name is optional,
c                       the default being blank.
c
 100  continue
      if(matchs('name',4)) call splunj
      if(label(dummy)) then
         name= ' '
         call entits(name,nc)
         if(nc.gt.8) nc=8
         stname(1:nc)= name(1:nc)
      else
         call errmsg(44,dum,dums,dumr,dumd)
      end if
      go to 10
c
c                       material command. properties relted to
c                       materials to be used in the problem are
c                       input here.
c
 200  continue
      call inmat(sbflg1,sbflg2,chkprm)
      go to 10
c
c                       number command. the number of nodes
c                       and/or elements is input here.
c
 300  continue
c
c                       make sure there is an attempt to input
c                       something. if there isn't, return for a
c                       high level command. fatal = .true. ind-
c                       icates that a fatal error has been made
c                       and execution will be halted.
c
      if(endcrd(dum)) then
         call errmsg(22,dum,dums,dumr,dumd)
      else
         call innum(sbflg1,sbflg2)
         if(fatal) go to 9999
      end if
      go to 10
c
c                       coordinates command. the uniform global
c                       coordinates of each node are input here
c                       (reference configuration)
c
 400  continue
c
c                       check to see if the number of nodes in the
c                       strucure has been set. if not, print error
c                       message and ignore the current command.
c                       return to high order commands to receive
c                       number of nodes command
c
      if(.not.numnod) then
         param= 1
         call errmsg(11,param,dums,dumr,dumd)
         go to 10
      else if(coor) then
         call errmsg(113,dum,dums,dumr,dumd)
         go to 10
      else
         call incoor(sbflg1,sbflg2)
         go to 20
      end if
c
c                       elements command. all necessary element
c                       properties for each element are set here.
c
c
 500  continue
c
c                       check to see if the number of elements in the
c                       structure has been set. if not, print error
c                       message and ignore the current command.
c                       return to high order commands to receive
c                       number of elements command.
c
      if(.not.numel) then
         param= 2
         call errmsg(11,param,dums,dumr,dumd)
         go to 10
      else if(elprop) then
         call errmsg(108,dum,dums,dumr,dumd)
         go to 10
      else
         call inelem(sbflg1,sbflg2)
         go to 20
      end if
c
c                       incidences command. the incidences for each
c                       element in the structure are input here.
c
 600  continue
c
c                       check to see if the element properties for
c                       each element and the number of nodes have been
c                       set. if not, print an error message and ignore
c                       the current command. return to high order
c                       commands to receive the necessary commands.
c
      if(.not.elprop.or..not.numnod) then
         param= 3
         call errmsg(11,param,dums,dumr,dumd)
         go to 10
      else if(elinc) then
         call errmsg(112,dum,dums,dumr,dumd)
         go to 10
      else
         call ininc(sbflg1,sbflg2)
         go to 20
      end if
c
c                       constraints command. the constraints in
c                       constraint compatable coordinates are input
c                       here, along with the transformation matrices
c                       defining the constraint compatable directions
c                       at each affected node. constraints are appli-
c                       ed at each  step, the same until changed.
c
 700  continue
c
c                       check to see if the element properties and
c                       incidences for each element and the number of
c                       nodes have been set. if not, print an error
c                       message and ignore the current command. return
c                       to high order commands to receive the necess-
c                       ary commands.
c
      if(.not.elprop.or..not.numnod.or..not.elinc) then
         param= 4
         call errmsg(11,param,dums,dumr,dumd)
         go to 10
      else
         call incon(sbflg1,sbflg2,olddof)
         go to 20
      end if
c
c                       loading command. all loadings for the structure
c                       are defined here. there can be up to 10 load-
c                       ings per strucure.
c
 800  continue
c
c                       check to see if the element properties and
c                       incidences for each element and the number of
c                       nodes have been set. if not, print an error
c                       message and ignore the current command. return
c                       to high order commands to receive the necess-
c                       ary commands.
c
      if(.not.elprop.or..not.numnod.or..not.elinc) then
         param= 4
         call errmsg(11,param,dums,dumr,dumd)
         go to 10
      else
         call inlod(sbflg1,sbflg2,chkprm,path)
         if(sbflg1) then
            go to 20
         else
            if(path.eq.0) then
               go to 20
            else
               go to 10
            end if
         end if
      end if
c
c                       available command
c
 900  continue
      go to 10
c
c                       nonlinear & dynamic analysis parameters.
c                       all parameters associated with the
c                       solution algorithm are input here.
c
 1000 continue
c
c                       check for proper command syntax and branch to
c                       appropriate subroutine if found. print error
c                       message if not.
c
      if(matchs('analysis',5)) goto 1000
      if(matchs('parameters',5)) then
         call indypm(sbflg1,sbflg2)
         go to 20
      else
         if(sbflg1) then
            call indypm(sbflg1,sbflg2)
            go to 20
         end if
         call errmsg(91,dum,dums,dumr,dumd)
         go to 10
      end if
c
c                       initial conditions command. initial displace-
c                       ments, velocities, and loads can be input
c                       here.
 1100 continue
c
      if(matchs('conditions',4)) call splunj
c
c                       check to see if the element properties and
c                       incidences for each element and the number of
c                       nodes have been set. if not, print an error
c                       message and ignore the current command. return
c                       to high order commands to receive the necess-
c                       ary commands.
c
      if(.not.elprop.or..not.numnod.or..not.elinc) then
         param= 4
         call errmsg(11,param,dums,dumr,dumd)
         go to 10
      else
         call inicon(sbflg1,sbflg2)
         go to 20
      end if
c
c                       compute command. advancement of the sol-
c                       ution and various other computational re-
c                       quests can be made here.
c
 1200 continue
c
c                       before any computations can be made, the
c                       structure must be completely defined. check
c                       for this condition and allow the request
c                       only if it exists.
c
      if (.not.input_ok) then
         call errmsg(183,dum,dums,dumr,dumd)
         goto 9999
      endif
      if(.not.coor.or..not.elprop.or..not.elinc.or..not.block) then
         param= 5
         call errmsg(11,param,dums,dumr,dumd)
      else
c
c                       print out # of error messages encountered so far
c                       and invoke compute to driver solution
c
         call error_count(out,.false.)
         call compute
      end if
      go to 10
c
c                       output command. output of various quantities
c                       can be requested here.
c
 1300 continue
c
c                       before any output can be made, the
c                       structure must be completely defined. check
c                       for this condition and allow the request
c                       only if it exists.
c
      if(.not.coor.or..not.elprop.or..not.elinc) then
         param= 5
         call errmsg(11,param,dums,dumr,dumd)
      else
         call oudrive(sbflg1,sbflg2,stname,ltmstp)
      end if
      go to 10
c
c
c
 1400 continue
c
c                       Convert binary packet file to ascii.
c
      if (matchs('binary',3))   call splunj
      if (matchs('packets',4))  call splunj
c
c                       Before packets can be converted, the structure must
c                       be completely defined. check for this condition and
c                       allow the request only if it exists.
c
      if(.not.coor.or..not.elprop.or..not.elinc) then
         param=5
         call errmsg(11, param, dums, dumr, dumd)
         go to 10
      end if
c
c                       Check for to make sure binary packet output is
c                       turned on.
c
      if( .not. output_packets ) then
         call errmsg2(67, dum, dums, dumr, dumd)
         go to 10
      end if

      call pconvdriv
      go to 10
c
c
c
c                       save command. the current structure in its
c                       present state is stored on disk. default
c                       save file name is structure name_db
c
 1500 continue
      stflnm = " "
      stcnam = " "
      stflnm(1:) = 'namenone'
      stcnam = stname
      if( endcrd(dummy) ) then
        call store( stcnam, stflnm, sbflg1, sbflg2 )
        go to 10
      end if
c
      if( matchs('to',2) ) call splunj
      stcnam(1:) = ' '
      stflnm(1:) = ' '
      if( matchs_exact('file') ) call splunj
      if ( endcrd(dummy) ) then
        stflnm(1:)= 'namenone'
        stcnam = stname
        call store( stcnam, stflnm, sbflg1, sbflg2 )
        go to 10
      end if
      stcnam(1:) = ' '
      stflnm(1:) = ' '
      if( label(dummy) ) then
           stflnm(1:) = ' '
           call entits( stflnm, nc )
           call store( stcnam, stflnm, sbflg1, sbflg2 )
         go to 10
      end if
      if( string(dummy) ) then
           stflnm(1:) = ' '
           call entits( stflnm, nc )
           call store( stcnam, stflnm, sbflg1, sbflg2 )
           go to 10
       end if
       write(*,9310)
       stflnm = 'namenone'
       stcnam = stname
       call store( stcnam, stflnm, sbflg1, sbflg2 )
       go to 10

c
c                       restart command. a structure previously stored
c                       on disk is read back. the current structure is
c                       destroyed if not saved before this step.
c
 1600 continue
      stcnam = " "
      rtflnm = " "
      if ( endcrd(dummy) ) then
       write(out,9300)
       go to 10
      end if
      stcnam = ' '
      if( matchs('structure',9) ) then
        if( label(dummy) ) then
          name = ' '
          call entits( name,nc )
          if( nc .gt. 8 ) nc = 8
          stcnam(1:nc)  = name(1:nc)
          rtflnm(1:nc)  =  stcnam(1:nc)
          rtflnm(nc+1:) = '_db'
        end if
        call reopen( stcnam, rtflnm, sbflg1, sbflg2 )
        go to 10
      end if
c
      if( matchs_exact('from') ) call splunj
      if( matchs_exact('file') ) call splunj
      if( label(dummy) ) then
          rtflnm = ' '
          call entits( rtflnm, nc )
      else if( string(dummy) ) then
          rtflnm = ' '
          call entits( rtflnm, nc )
      end if
      call reopen(stcnam,rtflnm,sbflg1,sbflg2)
      go to 10
c
      write(*,9200)
      go to 10
c
c                       1700 is Not Used at this .
c
 1700 continue
      write(*,9100)
      go to 10
c
c                       element blocking command. the particulars of
c                       the element blocking data structures are
c                       input here.
c
 1800 continue
      if(.not.numel) then
         param= 2
         call errmsg(11,param,dums,dumr,dumd)
         go to 10
      else if(block) then
         call errmsg(154,dum,dums,dumr,dumd)
         go to 10
      else
         call inelbk(sbflg1,sbflg2)
c                       Set up the crystal element properties as well
         call read_crystal_data
         call read_simple_angles
         call avg_cry_elast_props
         go to 20
      end if
c
c                       stop command. session is completed.
c
 1900 continue
      goto 9999
c
c                       domain computation command
c
c                       before any domain computations can be made, the
c                       structure must be completely defined. check
c                       for this condition and allow the request
c                       only if it exists.
c
 2000 continue
      if(.not.coor.or..not.elprop.or..not.elinc) then
         param= 5
         call errmsg(11,param,dums,dumr,dumd)
      else
         call indom( sbflg1, sbflg2 )
         go to 20
      end if
      go to 10
c
c                       crack growth parameters
c
 2100 continue
      if (sbflg1) then
         call incrack(sbflg1,sbflg2)
         goto 20
      else
         if (matchs('growth',4)) call splunj
         if (matchs('parameters',4)) call splunj
         call incrack(sbflg1,sbflg2)
         goto 20
      endif
      goto 10
c
c                       debug command.
c
 2200 continue
      if (matchs('off',2)) then
         debug1 = .false.
         debug2 = .false.
         debug = .false.
      else
         if (matchs('on',2)) call splunj
         debug1 = .true.
         debug2 = .true.
         debug = .true.
      endif
 2210 continue
      sbflg1 = .true.
      sbflg2 = .true.
      goto 10
c
c                       stress-strain command
c
 2300 continue
      call incurv(sbflg1, sbflg2)
      go to 20
c
c                       license command -- print out license agreement
c
 2400 continue
      call license
      go to 20
c
c                       contact plane command -- defines contact planes
c
 2500 continue
      if (matchs('plane',4)) call splunj
      if (matchs('surfaces',4)) call splunj
      call incontact(sbflg1,sbflg2)
      go to 20
c
c                       anisotropic thermal expansion coefficients
c
 2600 continue
      if( .not.elprop.or..not.numnod ) then
         param = 4
         call errmsg(11,param,dums,dumr,dumd)
      else
         call inalpha( sbflg1, sbflg2 )
         go to 20
      end if
      go to 20
c
c                       fgm properties at model nodes
c
 2700 continue
      if( .not.numnod ) then
         param = 4
         call errmsg(11,param,dums,dumr,dumd)
      else
         call infgm( sbflg1, sbflg2 )
         go to 20
      end if
      go to 20
c
c                       surfaces for tied contact
c
 2800 continue
      if (.not. numel) then
         param = 2
         call errmsg(11,param,dums,dumr,dumd)
      else
         call insurf()
         go to 20
      end if
      go to 20
c
c                       tied mesh sets
c
 2900 continue
      if (.not. numel) then
         param = 2
         call errmsg(11,param,dums,dumr,dumd)
      else
         call intied()
         go to 20
      end if
      go to 20
c
c                       table definition
c
 3000 continue
      call intab(sbflg1,sbflg2,chkprm,path)
      if(sbflg1) then
         go to 20
      else
         go to 10
      end if
c
c                       user defined list of something
c
 3100 continue
      call trwlist(sbflg1,sbflg2)
      if(sbflg1) then
         go to 20
      else
         go to 10
      end if
c
c                       crystal definition
c
 3200 continue
      call incrystal( sbflg1, sbflg2, chkprm, out )
      sbflg2 = .true.
      go to 10
c
c                       release constraints command.
c
 3300 continue
c
c                       check to see if the element properties and
c                       incidences for each element and the number of
c                       nodes have been set. if not, print an error
c                       message and ignore the current command. return
c                       to high order commands to receive the necess-
c                       ary commands.
c
      if(.not.elprop.or..not.numnod.or..not.elinc) then
         param= 4
         call errmsg(11,param,dums,dumr,dumd)
         go to 10
      else
         call release_constraints( sbflg1, sbflg2 )
         go to 20
      end if
      
c
c                       output timings. end execution.
c
 9999 continue
      call warp3d_normal_stop
c
 9000 format('    ********************************************',
     &       '***********************',/,
     &     '    **                                             ',
     &     '                  **',/,
     &     '    **                                             ',
     &     '                  **',/,
     &     '    **   W       W    AAAAA    RRRRRR    PPPPPP    ',
     &     '   33333  DDDD    **',/,
     &     '    **   W       W   A     A   R     R   P     P   ',
     &     '  3     3 D   D   **',/,
     &     '    **   W       W   A     A   R     R   P     P   ',
     &     '        3 D    D  **',/,
     &     '    **   W       W   A     A   R     R   P     P   ',
     &     '        3 D    D  **',/,
     &     '    **   W   W   W   AAAAAAA   RRRRRR    PPPPPP  --',
     &     '--  3333  D    D  **',/,
     &     '    **   W   W   W   A     A   R RR      P         ',
     &     '        3 D    D  **',/,
     &     '    **   W   W   W   A     A   R   RR    P         ',
     &     '  3     3 D   D   **',/,
     &     '    **    WWW WWW    A     A   R     RR  P         ',
     &     '   33333  DDDD    **',/,
     &     '    **                                             ',
     &     '                  **',/,
#win     &     '    **     Windows (Intel)            -dev-    Release: ',
#lnx     &     '    **     Linux (Intel)              -dev-    Release: ',
#l64     &     '    **     Intel 64-bit on Linux      -dev-    Release: ',
#mac     &     '    **     Mac OS X (Intel)           -dev-    Release: ',
     &     ' 17.7.3      **',/,
#win     &     '    **     Code Build Number: 3204             ',
#win     &     '                     **',/,
!win     &     "    **     Built on: Fri Jun 17 14:56:17 EDT 2016 ",
!win     &     '                   **',/,
     &     '    **     University of Illinois @ U-C.',
     &     '    Civil & Env Engineering  **',/,
     &     '    **     Today: ',a24,27x,'**',/,
     &     '    **                                             ',
     &     '                  **',/,
     &     '    **     NOTICE:  Use of Program Implies Agreement with',
     &     ' Terms &    **',/,
     &     '    **              Conditions Set Forth in File',
     &     " 'license_agreement' **",/,
     &     "    **              Enter the Command 'license' to",
     &     ' Display Text      **',/,
     &     '    **                                             ',
     &     '                  **',/,
     &     '    **     Limits (nodes, elements): ',2i10,
     &     '            **',/,
     &     '    **                                             ',
     &     '                  **',/,
     &     '    ***********************************************',
     &     '********************')
c
 9100 format('>>>>> error: the structure library no longer exists,',
     &       ' and thus this command',/,7x,'is no longer valid.')
 9200 format(//'>>>>> error: not able to construct restart file',
     &  ' name...',//)
 9300 format(//'>>>>> error: no file or structure name given...')
 9310 format(//'>>>>> error: no save file name given. default used...')
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine error_count                  *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 05/29/2016 rhd             *
c     *                                                              *
c     *     write the summary of # of errors, warnings, and fatal    *
c     *     errors encountered during the run.                       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine error_count( outdev, clear )
      implicit integer (a-z)
$add common.main
      logical already_called, clear
      data already_called / .false. /
      save already_called
c
c            clear specifies if the message counts are to be cleared.
c            already_called is how we tell if we are outputing on the
c            first step or not.
c
c            if we jill execution, make sure to also write the
c            terminate message to stdout [ variable out in common.main]
c
     &  already_called
      if( .not. already_called ) then
       if( num_warn.gt.0 .or. num_error.gt.0 .or.
     &      num_fatal.gt.0 ) then
            write(outdev,9000) num_warn, num_error, num_fatal
       end if
       if( clear ) already_called = .true.
       if( .not. clear .and. num_fatal .gt. 0 .or.
     &       num_error .gt. 0 ) then
         write(outdev,9200)
         call die_abort
       end if
      else  !  already called
       if( num_warn.gt.0 .or. num_error.gt.0 .or.
     &      num_fatal.gt.0 ) then
           write(outdev,9100) num_warn, num_error, num_fatal
           if( outdev .ne. out ) write (out,9100) num_warn, 
     &                            num_error, num_fatal
       end if
       if( .not. clear .and. num_fatal .gt. 0 .or.
     &       num_error .gt. 0 ) then
         write(outdev,9200)
         if( outdev .ne. out ) write(out,9200)
         call die_abort
       end if
      end if
c
      if( clear ) then
         num_warn  = 0
         num_error = 0
         num_fatal = 0
      end if
c
 9000 format (/,' >>>>> Messages before first step:',/,
     &          ' >>         Warnings:      ',i4,/,
     &          ' >>         Errors:        ',i4,/,
     &          ' >>         Fatal Errors:  ',i4,/)
c
 9100 format (/,' >>>>> Messages since last step:',/,
     &          ' >>         Warnings:      ',i4,/,
     &          ' >>         Errors:        ',i4,/,
     &          ' >>         Fatal Errors:  ',i4,/)
c
 9200 format (/,' >>>>> Errors Prevent Analysis.....',/,
     &          ' >>    Job Terminated',/)
      return
c
      end

c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine warp3d_normal_stop           *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/20/2103                 *
c     *                                                              *
c     *     execute a normal shutdown with messages, etc.            *
c     *                                                              *
c     ****************************************************************
c
c

      subroutine warp3d_normal_stop
      
      use file_info
      use main_data, only : output_packets
      use performance_data, only : time_assembly, assembly_total,
     &            ntimes_assembly
      implicit integer (a-z)
$add common.main
      real t1, wcputime, dumr
      external wcputime
      character*8 stcnam, dums,
     &     sdate_*24
      character*80 name, stflnm, rtflnm
      logical hilcmd,sbflg1,sbflg2
      logical endcrd,label,matchs,debug1,debug2,debug,endfil,
     &        string, matchs_exact
c
      common/errprm/ erprmd(10),erprmr(10),erprmi(10),erprms
#dbl      double precision
#sgl      real
     &  erprmd, dumd
      real erprmr
      character erprms *50
      real :: wall, real_start, real_end, real_rate
c
      common/erflgs/ numnod,numel,fatal,coor,elprop,elinc,constr,block
      logical fatal,coor,numnod,numel,elprop,elinc,constr,block
c
c                       cleanup some allocs first
c
      call cleanup_crystal
c
      call outime
      if( time_assembly ) then
       write(out,*)
       write(out,*)
       write(out,'(">> avg. assembly wall time (secs): ", f12.6)')
     &                  assembly_total/ dble(ntimes_assembly)
      end if

      write(out,*)
      write(out,*)
      t1 = wcputime ( 1 )
      write(out,'(">> total job wall time (secs): ", f12.2)') t1
c
c                      close input and output files
c
      if( outing ) close (out)
      if( filcnt .gt. 1 ) then
         do i = filcnt, 2, -1
            close(i)
         end do
      end if
c
c          uexternaldb for Abaqus compatible support
c
      douextdb = 2   ! in common. tell uexternaldb to terminate
      call wmpi_do_uexternaldb
c
c         MPI:
c            tell workers we are ending now and then stop.
c         threads - Fortran stop
c
      call die_gracefully
c
      end




