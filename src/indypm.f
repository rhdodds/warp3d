c     ****************************************************************
c     *                                                              *
c     *                      subroutine indypm                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 1/6/2016 rhd               *
c     *                                                              *
c     *     input parameters controlling how the solution is         *
c     *     performed for analysis                                   *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine indypm( sbflg1, sbflg2 )
      use global_data ! old common.main
c
      use mod_mpc, only : display_mpcs
      use main_data, only : output_packets, packet_file_name,
     &                      run_user_solution_routine, cp_unloading,
     &                      divergence_check, diverge_check_strict,
     &                      asymmetric_assembly,
     &                      batch_mess_fname, extrapolate,
     &                      extrap_off_next_step, line_search,
     &                      ls_details, ls_min_step_length,
     &                      ls_max_step_length, ls_rho,
     &                      ls_slack_tol, umat_serial
      use hypre_parameters
      use performance_data
      use distributed_stiffness_data, only : parallel_assembly_allowed,
     &                                       initial_map_type,
     &                                       final_map_type
c
      implicit none
c
      logical :: sbflg1, sbflg2
c
      integer :: intlst(mxlsz), dum, lenlst, param, i,testyp, nc, idum,
     &           errnum, trctyp, ncerror
      real :: dumr
      double precision ::  dnum, dumd
      double precision, parameter :: zero = 0.0d00
      logical :: linlst, lstflg, msg_flag, local_direct_flag,
     &           local_direct
      logical, external :: matchs, integr, endcrd, true, numd, numr,
     &                     string, numi, label, matchs_exact, match
      character(len=1) :: dums
      character (len=50) :: error_string
c
      msg_flag = .true.
c
c                       if sub flag 1 is on, indypm has been re-
c                       entered and the last command was in error.
c
      if( sbflg1 ) call errmsg(92,dum,dums,dumr,dumd)
c
c                       set flag indicating that the dynamic analysis
c                       parameters have changed.
c
      new_analysis_param = .true.
c
c
 10   call readsc
c
c                       branch on type of dynamic analysis parameter
c                       command.
c
 20   if( matchs('print',5)         ) go to 100
      if( matchs('newmark',4)       ) go to 300
      if( matchs('time',4)          ) go to 400
      if( matchs('extrapolate',7)   ) go to 500
      if( matchs('convergence',5)   ) go to 700
      if( matchs('nonconvergent',6) ) go to 800
      if( matchs('maximum',3)       ) go to 900
      if( matchs('minimum',3)       ) go to 950
      if( matchs('linear',6)        ) go to 1000 ! deprecated. send warn
      if( matchs('solution',8)      ) go to 1100
      if( matchs('trace',5)         ) go to 1200
      if( matchs('adaptive',5)      ) go to 1500
      if( matchs('batch',5)         ) go to 1600
      if( matchs('cpu',3)           ) go to 1700
      if( matchs('wall',3)          ) go to 1700
      if( matchs('material',4)      ) go to 1800
      if( matchs('bbar',4)          ) go to 2100
      if( matchs('consistent',5)    ) go to 2200
      if( matchs('solver',6)        ) go to 2300
      if( matchs('show',4)          ) go to 2400
      if( matchs('reset',5)         ) go to 2500
      if( matchs('sparse',5)        ) go to 2600
      if( matchs('binary',5)        ) go to 2700
      if( matchs('display',4)       ) go to 2800
      if( matchs('hypre',5)         ) go to 2900
      if( matchs('umat',4)          ) go to 3000
      if( matchs('assembly',5)      ) go to 3100
      if( matchs('user_routine',6)  ) go to 3200
      if( matchs('unloading',6)     ) go to 3300
      if( matchs('divergence',5)    ) go to 3400
      if( matchs_exact('line')      ) go to 3500 ! line search
c
c                       no match with solutions parameters command.
c                       return to driver subroutine to look for high
c                       level command.
c
      go to 9999
c
c
c **********************************************************************
c *                                                                    *
c *                     print residual loads command.                  *
c *                                                                    *
c **********************************************************************
c
c
 100  continue
c
      if(matchs('linear',6)) then
         linlst= .true.
      else
         linlst= .false.
      end if
c
      if(matchs('residual',3)) then
c
 110     if(matchs('loads',4)) then
            go to 110
         else if(matchs('for',3)) then
            go to 110
         else if(matchs('iterations',4)) then
            go to 110
         else
            call backsp(1)
         end if
c
         call scan
         call trlist(intlst,mxlsz,mxiter,lenlst,errnum)
c
c                    branch on the return code from trlist. a
c                    value of 1 indicates no error. a value of
c                    2 indicates that the parse rules failed in
c                    the list. a value of 3 indicates that the
c                    list overflowed its maximum length of mxlsz.
c                    a value of 4 indicates that no list was found.
c                    for the first two error cases, set the print
c                    residuals flag to false. for the third, assume
c                    that no list means all
c
         if(errnum.eq.1) then
            lstflg= .true.
         else if(errnum.eq.2) then
            param= 1
            call errmsg(24,param,dums,dumr,dumd)
            lstflg= .false.
         else if(errnum.eq.3) then
            param= 2
            call errmsg(24,param,dums,dumr,dumd)
            lstflg= .false.
         else if(errnum.eq.4) then
            intlst(1)= 1
            intlst(2)= -mxiter
            intlst(3)= 1
            lenlst= 3
            lstflg= .true.
         else
            param= 3
            call errmsg(24,param,dums,dumr,dumd)
            lstflg= .false.
         end if
c
c
         if(linlst) then
c
            prlres= lstflg
            if(lstflg) then
               do i= 1,lenlst
                  plrlst(i:lenlst)= intlst(i:lenlst)
               end do
               nplrs= lenlst
            end if
c
         else
c
            prnres= lstflg
            if(lstflg) then
               do i= 1,lenlst
                  prslst(i:lenlst)= intlst(i:lenlst)
               end do
               nprs= lenlst
            end if
c
         end if
c
      else
c
         call errmsg(94,dum,dums,dumr,dumd)
c
      end if
c
      go to 10
c
c **********************************************************************
c *                                                                    *
c *                     newmark beta parameter                         *
c *                                                                    *
c **********************************************************************
c
c
 300  continue
c
      if(matchs('beta',4)) then
         if(numd(dnum)) then
            nbeta= dnum
         else
            call errmsg(97,dum,dums,dumr,dumd)
         end if
      else
         call errmsg(98,dum,dums,dumr,dumd)
      end if
c
      go to 10
c
c
c **********************************************************************
c *                                                                    *
c *                     time step increment parameter                  *
c *                                                                    *
c **********************************************************************
c
c
 400  continue
c
      if(matchs('step',4)) then
         if(matchs('del_t',5)) call splunj
         if(numd(dnum)) then
            dt= dnum
         else
            call errmsg(97,dum,dums,dumr,dumd)
         end if
      else
         call errmsg(99,dum,dums,dumr,dumd)
      end if
c
      go to 10
c
c
c **********************************************************************
c *                                                                    *
c *                     parameter controlling whether to extrapolate   *
c *                     the incremental displacement vector for        *
c *                     newton iterations                              *
c *                                                                    *
c **********************************************************************
c
c
 500  continue
c
      if( matchs('solution',4) ) call splunj
      if( matchs('on',2) ) then
         extrapolate = .true.
      else if( matchs('off',3) ) then
         extrapolate = .false.
      else
         call errmsg(174,dum,dums,dumr,dumd)
      end if
      if( endcrd() ) go to 10
      if( matchs('next',2) ) then
         extrap_off_next_step = .true.
      else
         call errmsg(170,dum,dums,dumr,dumd)
      end if
      go to 10
c
c **********************************************************************
c *                                                                    *
c *                     parameters controlling the type of converg-    *
c *                     ence tests to be used.                         *
c *                                                                    *
c **********************************************************************
c
c
 700  continue
c
c                       initialize the appropriate vectors so that
c                       only the requests currently being specified
c                       will be executed.
c
      do i = 1, mxcvtests
         convrg(i)= .false.
         tol(i)= zero
      end do
c
      if(matchs('tests',4)) call splunj
c
 710  if(matchs('norm',4)) then
c
         if(matchs('displacement',5)) then
            testyp= 1
            convrg(1)= .true.
            if(matchs('absolute',5)) then
              testyp= 5
              convrg(5)= .true.
              convrg(1)= .false.
            end if
            go to 720
         else if(matchs('residual',3)) then
            if(matchs('loads',4)) call splunj
            testyp= 2
            convrg(2)= .true.
            go to 720
         else
            call errmsg(105,dum,dums,dumr,dumd)
            if(true(dum)) go to 710
         end if
c
      else if(matchs('maximum',3)) then
c
         if(matchs('displacement',5)) then
            testyp= 3
            convrg(3)= .true.
            if(matchs('absolute',5)) then
              testyp= 6
              convrg(6)= .true.
              convrg(3)= .false.
            end if
            go to 720
         else if(matchs('residual',3)) then
            if(matchs('loads',4)) call splunj
            testyp= 4
            convrg(4)= .true.
            go to 720
         else
            call errmsg(105,dum,dums,dumr,dumd)
            if(true(dum)) go to 710
         end if
c
      else if(matchs(',',1)) then
c
         if(endcrd(dum)) call readsc
         go to 710
c
      else
c
         if(endcrd(dum)) then
            go to 10
         else
            call errmsg(106,dum,dums,dumr,dumd)
            if(true(dum)) go to 710
         end if
c
      end if
c
c                       input the tolerance by which the given con-
c                       vergence test is to be evaluated.
c
 720  if(matchs('tolerance',3)) call splunj
c
      if(numd(dnum)) then
         tol(testyp)= dnum
      else
         call errmsg(97,dum,dums,dumr,dumd)
         convrg(testyp)= .false.
      end if
      go to 710
c
c
c **********************************************************************
c *                                                                    *
c *                     parameter controlling solution if there is     *
c *                     not convergence after mxiter iterations.       *
c *                                                                    *
c **********************************************************************
c
c
 800  continue
c
      if(matchs('solutions',3)) then
         if(matchs('stop',4)) then
            halt= .true.
         else if(matchs('continue',4)) then
            halt= .false.
         else
            call errmsg(102,dum,dums,dumr,dumd)
         end if
      else
         call errmsg(107,dum,dums,dumr,dumd)
      end if
c
      go to 10
c
c
c **********************************************************************
c *                                                                    *
c *                     parameter defining the maximum number of       *
c *                     newton iterations                              *
c *                     to execute before possibly moving to the       *
c *                     next step if there is no convergence.          *
c *                                                                    *
c **********************************************************************
c
c
 900  continue
c
      if(matchs('iterations',4)) then
         if(.not.integr(mxiter)) then
            call errmsg(103,dum,dums,dumr,dumd)
         end if
      else
         call errmsg(104,dum,dums,dumr,dumd)
      end if
c
      go to 10
c **********************************************************************
c *                                                                    *
c *                     parameter defining the minimum number of       *
c *                     newton iterations                              *
c *                     to execute before moving to the next step      *
c *                                                                    *
c **********************************************************************
c
c
 950  continue
c
      if(matchs('iterations',4)) then
         if(.not.integr(mniter)) then
            call errmsg(103,dum,dums,dumr,dumd)
         end if
      else
         call errmsg(104,dum,dums,dumr,dumd)
      end if
c
      go to 10
c
c
c **********************************************************************
c *                                                                    *
c *                    *** deprecated ***                              *
c *             *** done now thru extrapolation command ***            *
c *                                                                    *
c *                     request use of linear stiffness                *
c *                     iteration 1 of each load step                  *
c *                                                                    *
c *  (a)  linear stiffness for iteration 1 on | off                    *
c *          (subsequent steps)                                        *
c *  (b)  linear stiffness for iteration one next step                 *
c *                                                                    *
c *    use these commands to replace                                   *
c *                                                                    *
c *       extrapolate off => same  as (a)                              *
c *       extrapolate off next step => same  as (b)                    *
c *                                                                    *
c **********************************************************************
c
c
1000  continue
      write(out,9000)
      write(out,*)
      go to 10
c
      if( matchs('stiffness',5) ) then
c
         if( matchs('for',3) ) call splunj
c
         if( matchs('iteration',4) ) then
            if( matchs('one',3) ) then
               if( matchs('on',2) ) then
                  extrapolate = .false. !lnkit1= .true.
                  extrap_off_next_step = .true.
               else if( matchs('off',3) ) then
                  extrapolate = .false. ! lnkit1= .false.
                  extrap_off_next_step = .true.
               else if( matchs('next',4) ) then ! step
                   extrap_off_next_step = .true.
               else
                  call errmsg(93,dum,dums,dumr,dumd)
               end if
            else
               call errmsg(153,dum,dums,dumr,dumd)
            end if
         else
            call errmsg(93,dum,dums,dumr,dumd)
         end if
      else if(matchs('mass',4)) then
         if(matchs('on',2)) then
            linmas= .true.
         else if(matchs('off',3)) then
            linmas= .false.
         else
            call errmsg(93,dum,dums,dumr,dumd)
         end if
      else
         call errmsg(96,dum,dums,dumr,dumd)
      end if
c
      go to 10

c
c **********************************************************************
c *                                                                    *
c *                     Selection of Solver Type                       *
c *                                                                    *
c *      See initst.f. change solver after step 1 not allowed          *
c *                                                                    *
c **********************************************************************
c
c
 1100 continue
      if( ltmstp .ne. 0 ) then
         call errmsg(207,dum,dums,dumr,dumd)
         go to 10
      end if
      solver_flag = 0
      solver_mkl_iterative = .false.
      local_direct_flag = .false.
      if ( matchs('technique',4) ) call splunj
      if ( matchs('type',4) ) call splunj
c
      if ( matchs('direct',6) ) then
        local_direct = .true.
        if ( endcrd(dum) ) then
c          if (use_mpi) then  ! removed for now
c            call errmsg(333,dum,dums,dumr,dumd)
c          end if
          solver_mkl_iterative = .false.
          solver_flag = 0
          go to 1150
        end if
      end if
c
      if ( matchs('pardiso',6) ) then
        local_direct_flag = .true.
        if ( endcrd(dum) ) then
c          if (use_mpi) then   ! removed for now
c            call errmsg(333,dum,dums,dumr,dumd)
c          end if
          solver_mkl_iterative = .false.
          solver_flag = 0
          go to 1150
        end if
      end if
c
c                Pardiso direct. try to allow most possible
c                combinations of words direct sparse
c
      if ( matchs('sparse',5) ) then
c        if (use_mpi) then
c            call errmsg(333,dum,dums,dumr,dumd)
c        end if
        if ( endcrd(dum) ) then
               solver_flag = 7
               solver_mkl_iterative = .false.
               go to 1150
        end if
        if (  matchs('direct',6) ) then
               solver_flag = 7
               solver_mkl_iterative = .false.
               go to 1150
        end if
        if (  matchs('windows',6) ) then
               solver_flag = 7
               solver_mkl_iterative = .false.
               go to 1150
        end if
        if ( local_direct_flag ) then
               solver_flag = 7
               solver_mkl_iterative = .false.
               go to 1150
        end if
      end if
c
c                 Pardiso sparse iterative
c
      if ( matchs('iterative',8)) then
c         if (use_mpi) then
c            call errmsg(333,dum,dums,dumr,dumd)
c         end if
         solver_mkl_iterative = .true.
         solver_flag = 7
         go to 1150
      end if
      if ( matchs('windows',6) ) solver_flag = 7
c
c                 Asymmetric pardiso (iterative or direct). Force
c                 asymmetric assembly.
c                 check code where assembly command is translated
c
      if ( matchs('asymmetric',4)) then
            solver_flag = 8
            asymmetric_assembly = .true.
      end if

      if ( matchs('iterative',4)) then
            solver_mkl_iterative = .true.
            go to 1150
      end if

      if ( matchs('direct',6)) then
            solver_mkl_iterative = .false.
            go to 1150
      end if
c
      if( matchs('hypre',5) ) solver_flag = 9
c
      if( matchs('cluster',4) ) then
        solver_flag = 10
        if( asymmetric_assembly ) solver_flag = 11
        solver_mkl_iterative = .false.
      end if
c
 1150 continue
      if( solver_flag .eq. 0 ) then
        write(*,*) ">>>>> invalid solver specified..."
        write(*,*) "......... solver: ", solver_flag
        write(*,*) "......... mkl_iter_flg: ", solver_mkl_iterative
        write(*,*) " "
      end if
      if( .not. use_mpi ) then
        if( solver_flag == 10 .or. solver_flag == 11 ) then
          num_error = num_error + 1
          write(out,9200)
          solver_flag = 7
        end if
        if( solver_flag == 9 ) then
          num_error = num_error + 1
          write(out,9210)
          solver_flag = 7
        end if
      end if
c      write(*,*) "......... solver: ", solver_flag
c      write(*,*) "......... mkl_iter_flg: ", solver_mkl_iterative
c      write(*,*) " "
      go to 10
c
c
c **********************************************************************
c *                                                                    *
c *                     trace solution convergence command             *
c *                                                                    *
c **********************************************************************
c
c
 1200 continue
c
c                       initialize the trace vector so that only the
c                       tracings currently being specified will be
c                       executed.
c
      trace(1) = .false.
      trace(2) = .false.
      trace(3) = .false.
      trace(4) = .false.
c
 1210 if(matchs('solution',3)) then
         trctyp= 1
         go to 1220
      else if(matchs('material',3)) then
         trctyp= 2
         go to 1220
      else if(matchs('acceleration',5)) then
         trctyp= 3
         go to 1220
      else if(matchs(',',1)) then
         if(endcrd(dum)) call readsc
         go to 1210
      else
         if(endcrd(dum)) then
            go to 10
         else
            call errmsg(172,dum,dums,dumr,dumd)
            if(true(dum)) go to 1210
         end if
      end if
c
 1220 if(matchs('on',2)) then
         trace(trctyp)= .true.
      else if(matchs('off',3)) then
         trace(trctyp)= .false.
      else if(matchs(',',1)) then
         if(endcrd(dum)) call readsc
         go to 1220
      else
         call errmsg(142,dum,dums,dumr,dumd)
      end if
c
      go to 1210
c
c **********************************************************************
c *                                                                    *
c *                     adaptive solution for static problems          *
c *                                                                    *
c **********************************************************************
c
c
 1500 continue
      if( matchs('solution',3) ) call splunj
      if( matchs('on',      2) ) adaptive_flag = .true.
      if( matchs('off',     2) ) adaptive_flag = .false.
      go to 10
c
c
c **********************************************************************
c *                                                                    *
c *        batch messages on/off and file to use                       *
c *        after each newton iteration to show convergence info        *
c *                                                                    *
c **********************************************************************
c
c
 1600 continue
      if( matchs('messages',4)  ) call splunj
      if( matchs('on',      2)  ) batch_messages = .true.
      if( matchs('off',     2)  ) batch_messages = .false.
      batch_mess_fname(1:) = " "
      if( endcrd( dum ) ) go to 10
      if ( matchs( 'file',4 ) ) call splunj
      if ( label( dumr ) ) then
          call entits( batch_mess_fname(1:), nc )
      else if ( string( dumr ) ) then
          call entits( batch_mess_fname(1:), nc )
      else
          call errmsg(322,dum,dums,dumr,dumd)
      end if
      if( nc. eq. 80 ) write(out,9100)
      go to 10
c
c
c **********************************************************************
c *                                                                    *
c *                     set wall time limit -- used in systems where    *
c *                     there is a job queueing system so that data    *
c *                     wont be lost if the queued time runs out.      *
c *                                                                    *
c **********************************************************************
c
c
 1700 continue
      if (matchs('time',4)) call splunj
      if (matchs('limit',5)) call splunj
      if (matchs('off',3)) then
         time_limit = 0.0
      else
         if (matchs('on',2)) call splunj
         if (numr(time_limit)) then
            call errmsg(197,dum,dums,time_limit,dumd)
         else
            call errmsg(196,dum,dums,dumr,dumd)
         endif
      endif
      call steptime ( dum, 4 )
      go to 10
c
c
c **********************************************************************
c *                                                                    *
c *                     specifies if message about material            *
c *                     behavior should be printed or not              *
c *                                                                    *
c **********************************************************************
c
c
 1800 continue
      if (matchs('messages',4)) call splunj
      if (matchs('on',2)) then
         signal_flag = .true.
      else if (matchs('off',3)) then
         signal_flag = .false.
      endif
      go to 10
c
c
c **********************************************************************
c *                                                                    *
c *      input the b-bar element stabilization factor                  *
c *                                                                    *
c **********************************************************************
c
c
 2100 continue
c
      if( matchs('stabilization',4) ) call splunj
      if( matchs('factor',4) ) call splunj
      if( numd(dnum) ) then
        eps_bbar = dnum
      else
        call errmsg(97,dum,dums,dumr,dumd)
      end if
      go to 10
c
c
c **********************************************************************
c *                                                                    *
c *      input the consistent [Q] matrix flag - on/off                 *
c *                                                                    *
c **********************************************************************
c
c
 2200 continue
c
      if( matchs('q',1) ) call splunj
      if( matchs('-',1) ) call splunj
      if( matchs('matrix',4) ) call splunj
      if( matchs('on',2) ) then
         qbar_flag = .true.
      elseif( matchs('off',3) ) then
         qbar_flag = .false.
      else
        call errmsg(271,dum,dums,dumr,dumd)
      end if
      go to 10
c
c **********************************************************************
c *                                                                    *
c *      input the various solver options                              *
c *            out-of-core on|off                                      *
c *            memory <integer>      (units of MB)                     *
c *            scratch directory     <string>                          *
c *                                                                    *
c *                                                                    *
c **********************************************************************
c
c
 2300 continue
c
      if( matchs('out',3) ) then
        if ( matchs('-',1) )      call splunj
        if ( matchs('of',2) )     call splunj
        if ( matchs('-',1) )      call splunj
        if ( matchs('memory',3) ) call splunj
        if ( matchs('core',3) )   call splunj
        solver_out_of_core = .false.
        if ( matchs('on',2) ) then
          solver_out_of_core = .true.
        elseif ( matchs('off',3) ) then
          solver_out_of_core = .false.
        else
          call errmsg(280,dum,dums,dumr,dumd)
        end if
      elseif ( matchs('memory',3) ) then
        solver_memory = 0
        if ( .not. numi( solver_memory ) ) then
          call errmsg(277,dum,dums,dumr,dumd)
        end if
      elseif ( matchs('scratch',3) ) then
        if ( matchs('directory',3) ) call splunj
        if ( string( dumr ) ) then
          solver_scr_dir(1:) = ' '
          call entits( solver_scr_dir, dum )
        else
          call errmsg(278,dum,dums,dumr,dumd)
        end if
      else
        call errmsg(279,dum,dums,dumr,dumd)
      end if
      go to 10
c **********************************************************************
c *                                                                    *
c *                     parameter controlling whether to show          *
c *                     detailed messages during nonlinear solution    *
c *                                                                    *
c **********************************************************************
c
c
 2400 continue
      if( matchs('details',4) ) call splunj
      if( matchs('on',2) ) then
         show_details = .true.
      else if( matchs('off',3) ) then
         show_details = .false.
      else
         call errmsg(287,dum,dums,dumr,dumd)
      end if
c
      go to 10
c
c **********************************************************************
c *                                                                    *
c *                     reset the current permanent load reduction     *
c *                     factor to the specified value                  *
c *                                                                    *
c **********************************************************************
c
c
 2500 continue
      if( matchs('load',4) ) then
         if( matchs('reduction',4) ) then
           if( matchs('factor',4) ) call splunj
           if ( numr(dumr) ) then
               call insave_value( dumr, dum, 1 )
           else
               call errmsg(5,dum,dums,dumr,dumd)
           end if
           go to 10
         else
           call errmsg(316,dum,dums,dumr,dumd)
           go to 10
         end if
      else
         call errmsg(316,dum,dums,dumr,dumd)
      end if
      go to 10
c
c **********************************************************************
c *                                                                    *
c *                     command to request output of the assembled     *
c *                     sparse matrix in a formatted or unformatted    *
c *                     file                                           *
c *                                                                    *
c **********************************************************************
c
c
 2600 continue
      if ( matchs( 'stiffness', 5 ) ) call splunj
      if ( matchs( 'output', 5 )    ) call splunj
      sparse_stiff_output = .false.
      sparse_stiff_binary = .true.
      sparse_research     = .false.
      sparse_stiff_file_name(1:) = ' '
      sparse_stiff_file_name(1:) = 'sparse_stiffness_output'
      if ( matchs( 'on', 2 ) ) then
        sparse_stiff_output = .true.
      else if ( matchs( 'yes', 3 ) ) then
        sparse_stiff_output = .true.
      else if ( matchs( 'off', 3 ) ) then
        sparse_stiff_output = .false.
      else if ( matchs( 'no', 3 ) ) then
        sparse_stiff_output = .false.
      else
           call errmsg(271,dum,dums,dumr,dumd)
           go to 10
      end if
      if ( .not. sparse_stiff_output ) go to 10
      if ( endcrd(dum) ) go to 10
      if ( matchs( 'research',5 )    ) sparse_research = .true.
      if ( matchs( 'formatted',5 ) ) sparse_stiff_binary = .false.
      if ( matchs( 'unformatted',5 ) ) sparse_stiff_binary = .true.
      if ( matchs( 'binary',5 )      ) sparse_stiff_binary = .true.
      if ( endcrd(dum) ) go to 10
      if ( matchs( 'file',4 ) ) call splunj
      sparse_stiff_file_name(1:) = ' '
      if ( label( dumr ) ) then
          call entits( sparse_stiff_file_name(1:), dum )
      else if ( string( dumr ) ) then
          call entits( sparse_stiff_file_name(1:), dum )
      else
           call errmsg(322,dum,dums,dumr,dumd)
      end if
      go to 10
c
c **********************************************************************
c *                                                                    *
c *                     command to request output of the binary        *
c *                     result packets                                 *
c *                                                                    *
c **********************************************************************
c
c
 2700 continue
      if ( .not. matchs( 'packets', 5 ) ) then
           call errmsg2(24,dum,dums,dumr,dumd)
           go to 10
      end if
      output_packets = .false.
      call close_packets_file(msg_flag)
      packet_file_name(1:) = ' '
      packet_file_name(1:) = 'warp_binary_packets'
      if ( matchs( 'on', 2 ) ) then
        output_packets = .true.
      else if ( matchs( 'yes', 3 ) ) then
        output_packets = .true.
      else if ( matchs( 'off', 3 ) ) then
        output_packets = .false.
      else if ( matchs( 'no', 3 ) ) then
        output_packets = .false.
      else
           call errmsg(24,dum,dums,dumr,dumd)
           output_packets = .false.
           go to 10
      end if
      if ( .not. output_packets ) go to 10
      if ( endcrd(dum) ) go to 10
      if ( .not. matchs( 'file',4 ) ) then
          call errmsg(24,dum,dums,dumr,dumd)
          output_packets = .false.
      end if
      if ( label( dumr ) ) then
          packet_file_name(1:) = ' '
          call entits( packet_file_name(1:), dum )
      else if ( string( dumr ) ) then
          packet_file_name(1:) = ' '
          call entits( packet_file_name(1:), dum )
      else
           call errmsg(25,dum,dums,dumr,dumd)
           output_packets = .false.
      end if
      call open_packets_file(msg_flag)
      go to 10
c **********************************************************************
c *                                                                    *
c *                     parameter controlling whether to display       *
c *                     the mpc eqns generated by the tied mesh        *
c *                     processor                                      *
c *                                                                    *
c **********************************************************************
c
c
 2800 continue
      if( matchs('tied',4) ) call splunj
      if( matchs('mesh',4) ) call splunj
      if( matchs('mpcs',4) ) then
         if( matchs('on',2) ) then
            display_mpcs = .true.
         else if( matchs('off',3) ) then
            display_mpcs = .false.
         else
            call errmsg(63,dum,dums,dumr,dumd)
         end if
      end if
c
      go to 10
c
c **********************************************************************
c *                                                                    *
c *                     all hypre parameters                           *
c *    Input parameters for hypre iterative solver and                 *
c *    parasails or boomeramg preconditioner.                          *
c *                                                                    *
c **********************************************************************
c
c
 2900 continue
      if (matchs('printlevel',5)) then
            if(.not.integr(idum)) then
               call errmsg(103,dum,dums,dumr,dumd)
            else
                precond_printlevel = idum
                solver_printlevel = idum
            end if
      else if (matchs('preconditioner',7) ) then
            if (matchs('parasails',4) ) then
                  precond_type = 1
            else if (matchs('boomeramg',4) ) then
                  precond_type = 2
            else
                  call errmsg(332,dum,dums,dumr,dumd)
            end if
      else if (matchs('coarsening',7)) then
            if (matchs('CLJP',4)) then
                  coarsening = 0
            elseif (matchs('falgout',7)) then
                  coarsening = 6
            elseif (matchs('PMIS',4)) then
                  coarsening = 8
            elseif (matchs('HMIS',4)) then
                  coarsening = 10
            else
                  call errmsg(332,dum,dums,dumr,dumd)
            end if
      else if (matchs('solver',6 ) ) then
            if (matchs('pcg',3) ) then
                  hsolver_type = 1
            else
                  call errmsg(332,dum,dums,dumr,dumd)
            end if
      else if (matchs('interpolation',6)) then
            if (matchs('classical',7)) then
                  interpolation = 0
            elseif (matchs('direct',6)) then
                  interpolation = 3
            elseif (matchs('standard',8)) then
                  interpolation = 9
            elseif (matchs('multipass',5)) then
                  interpolation = 5
            elseif (matchs('ext_classical',3)) then
                  interpolation = 6
            else
                  call errmsg(332,dum,dums,dumr,dumd)
            end if
      else if (matchs('relaxation',6)) then
            if (matchs('jacobi',6)) then
                  relaxation = 0
            elseif (matchs('gs',2)) then
                  relaxation = 6
            else
                  call errmsg(332,dum,dums,dumr,dumd)
            end if
      else if (matchs('max_levels',9) ) then
            if(.not.integr(max_levels)) then
               call errmsg(103,dum,dums,dumr,dumd)
            else
            end if
      else if (matchs('levels',3) ) then
            if(.not.integr(levels)) then
               call errmsg(103,dum,dums,dumr,dumd)
            else
            end if
      else if (matchs('threshold',6) ) then
            if( numd(dnum) ) then
                  threshold = dnum
            else
                  call errmsg(91,dum,dums,dumr,dumd)
            end if
      else if (matchs('mg_threshold',8) ) then
            if( numd(dnum) ) then
                  mg_threshold = dnum
            else
                  call errmsg(91,dum,dums,dumr,dumd)
            end if
      else if (matchs('filter',4) ) then
            if( numd(dnum) ) then
                  filter = dnum
            else
                  call errmsg(91,dum,dums,dumr,dumd)
            end if
      else if (matchs('symmetry',3) ) then
            if(.not.integr(symme)) then
               call errmsg(103,dum,dums,dumr,dumd)
            else
            end if
      else if (matchs('balance',3) ) then
            if( numd(dnum) ) then
                  loadbal = dnum
            else
                  call errmsg(91,dum,dums,dumr,dumd)
            end if
      else if (matchs('tolerance',4) ) then
            if (numd(dnum)) then
                  hypre_tol = dnum
            else
                  call errmsg(91,dum,dums,dumr,dumd)
            end if
      else if (matchs('iterations',5) ) then
            if (.not. integr(hypre_max) ) then
                  call errmsg(103,dum,dums,dumr,dumd)
            else
            end if
      else if (matchs('cycle_type',5))  then
            if (matchs('V',1)) then
                  cycle_type = 1
            elseif (matchs('W',1)) then
                  cycle_type = 2
            else
                  call errmsg(332,dum,dums,dumr,dumd)
            end if
      else if (matchs('sweeps',6)) then
            if (.not. integr(sweeps) ) then
                  call errmsg(103,dum,dums,dumr,dumd)
            end if
      else if (matchs('agg_levels',5) ) then
            if (.not. integr(agg_levels) ) then
                  call errmsg(103,dum,dums,dumr,dumd)
            end if
      else if (matchs('truncation',5) ) then
            if (numd(dnum) ) then
                  truncation = dnum
            else
                  call errmsg(103,dum,dums,dumr,dumd)
            end if
      else if (matchs('wt_relax',8) ) then
            if (numd(dnum) ) then
                  relax_wt = dnum
            else
                  call errmsg(103,dum,dums,dumr,dumd)
            end if
      else if (matchs('wt_outer',8) ) then
            if (numd(dnum) ) then
                  relax_outer_wt = dnum
            else
                  call errmsg(103,dum,dums,dumr,dumd)
            end if
      else if (matchs('cs',2) ) then
            if (matchs('T',1) ) then
                  cf = 1
            elseif (matchs('F',1)) then
                  cf = 0
            else
                  call errmsg(332,dum,dums,dumr,dumd)
            end if
      else
            call errmsg(332,dum,dums,dumr,dumd)
      end if
c
      go to 10
c
c **********************************************************************
c *                                                                    *
c *                     Parameters to set num threads                  *
c *                 for use in various parts of the code               *
c *                                                                    *
c **********************************************************************
c
c
 3000 continue
      if( matchs("serial",6) ) then
        if (matchs('on',2)) then
           umat_serial = .true.
        else if (matchs('off',3)) then
           umat_serial = .false.
        else
           call errmsg(343,dum,dums,dumr,dumd)
        end if
      else
        call errmsg(332,dum,dums,dumr,dumd)
      end if
c
      go to 10
c
c
c **********************************************************************
c *                                                                    *
c *                     Parameters dealing with                        *
c *                    parallel and asymmetric assembly routines       *
c *                                                                    *
c **********************************************************************
c
c
 3100 continue
      if (matchs('parallel',4)) then
            if (matchs('on',2)) then
                  parallel_assembly_allowed=.true.
            else if (matchs('off',3)) then
                  parallel_assembly_allowed=.false.
            else
                  call errmsg(343,dum,dums,dumr,dumd)
            end if
      else if (matchs('initial',4)) then
            if (matchs('srows',4)) then
                  initial_map_type = 1
            else if (matchs('brows',4)) then
                  initial_map_type = 2
            else
                  call errmsg(342,dum,dums,dumr,dumd)
            end if
      else if (matchs('final',5)) then
            if (matchs('initial',7)) then
                  final_map_type = 1
            else if (matchs('srows',4)) then
                  final_map_type = 2
            else if (matchs('mat-vec',7)) then
                  final_map_type = 3
            else
                  call errmsg(342,dum,dums,dumr,dumd)
            end if
      else if (matchs('time',4)) then
            if (matchs('on',2)) then
                  time_assembly = .true.
            else if (matchs('off',3)) then
                  time_assembly = .false.
            else
                  call errmsg(343,dum,dums,dumr,dumd)
            end if
      else if (matchs('asymmetric',5)) then
            if (matchs('on',2)) then
                  asymmetric_assembly = .true.
            else if (matchs('off',3)) then
                  asymmetric_assembly = .false.
            else
                  call errmsg(343,dum,dums,dumr,dumd)
            end if
      else
            call errmsg(340,dum,dums,dumr,dumd)
      end if

      go to 10
c
c **********************************************************************
c *                                                                    *
c *      run on/off of the user routine to set solution parameters     *
c *      and/or the definition of the next loading step.               *
c *                                                                    *
c **********************************************************************
c
c
 3200 continue
      if( matchs('solution',4) ) call splunj
      if( matchs('routine',4) ) call splunj
      if( matchs('parameters',3) ) call splunj
      if( matchs('nonlinear',3) ) call splunj
      if( matchs('on',2) ) then
          run_user_solution_routine = .true.
      else if( matchs('off',3) ) then
          run_user_solution_routine = .false.
      else
          call errmsg(326,dum,dums,dumr,dumd)
      end if
      go to 10
c
c
c
c **********************************************************************
c *                                                                    *
c *                     parameter controlling the cp unloading         *
c *                     behavior                                       *
c *                                                                    *
c **********************************************************************
c
c
 3300  continue
c
      if(matchs('on',2)) then
c
         cp_unloading = .true.
c
      else if(matchs('off',3)) then
         cp_unloading = .false.
      else
         call errmsg(350,dum,dums,dumr,dumd)
      end if
c
      go to 10
c
c
c **********************************************************************
c *                                                                    *
c *                     divergence check. default is on. mnralg will   *
c *                     trigger adaptive if solution seems to be non-  *
c *                     convergence                                    *
c *                                                                    *
c **********************************************************************
c
c
 3400  continue
c
      if( matchs('check',2) ) call splunj
      if( matchs('on',2) ) then
         divergence_check = .true.
         if( match('strict',4) ) diverge_check_strict = .true.
      else if( matchs('off',3) ) then
         divergence_check = .false.
         diverge_check_strict = .false.
      else
         call errmsg2(81,dum,dums,dumr,dumd)
      end if
      go to 10
c
c **********************************************************************
c *                                                                    *
c *                     line search and parameters                     *
c *                                                                    *
c **********************************************************************
c
c
 3500 continue  ! make sure consistent with initst
      line_search    =  .true.
      ls_details     =  .false.
      ls_min_step_length = 0.01d00
      ls_max_step_length = 1.0d00
      ls_rho         = 0.7d00
      ls_slack_tol   = 0.5d00
c
      if( matchs('search',5) ) call splunj
      if( matchs_exact('on') ) line_search = .true.
      if( matchs_exact('off') ) then
         line_search = .false.
         call scan_flushline
         go to 10
      end if
c
 3510 continue
       if( endcrd( ) ) go to 10
       if( matchs_exact(',') ) call splunj
       if( matchs_exact('rho') ) then
          if( numd(ls_rho) ) go to 3510
          call entits( error_string, ncerror )
          write(out,9514) error_string(1:ncerror)
          call scan_flushline
          go to 10
       end if
       if( matchs('details',5) ) then
          ls_details = .true.
          go to 3510
       end if
       if( matchs_exact('min_step') ) then
          if( numd(ls_min_step_length) ) go to 3510
          call entits( error_string, ncerror )
          write(out,9516)  error_string(1:ncerror)
          call scan_flushline
          go to 10
       end if
       if( matchs_exact('max_step') ) then
          if( numd(ls_max_step_length) ) go to 3510
          call entits( error_string, ncerror )
          write(out,9518)  error_string(1:ncerror)
          call scan_flushline
          go to 10
       end if
       if( matchs('slack_toler',5) ) then
          if( numd(ls_slack_tol) ) go to 3510
          call entits( error_string, ncerror )
          write(out,9520)  error_string(1:ncerror)
          call scan_flushline
          go to 10
       end if
       if( matchs('dump',4) ) then
         write(out,9510) line_search, ls_details,
     &          ls_min_step_length, ls_max_step_length, ls_rho,
     &          ls_slack_tol
         go to 3510
       end if
       call entits( error_string, ncerror )
       write(out,9512) error_string(1:ncerror)
       num_error = num_error + 1
       call scan_flushline
       go to 10
c
 9999 sbflg1 = .true.
      sbflg2 = .false.
c
c
      return
c
 9000 format(/1x'>>>>>> Warning: the option:  linear stiffness ... ',
     &   'has been deleted from WARP3D.',
     &/,1x,'                Use the option: extrapolate off   to '
     &      ' force use of linear [D] at start of subsequent steps.',
     & /,1x,'                Or use the option: extrapolate off  ',
     &   'next step'
     & /,1x,'                The: linear stiffness ... command will',
     & ' be deleted in a coming release. '  )
 9100 format(/1x,'>>>>> Warning:  the file name has been truncated ',
     &           'to 80 characters...',/)
 9200 format(/1x,'>>>>> Error: cluster solver not compatible ',
     &   /16x,'with threads-only execution.'/)
 9210 format(/1x,'>>>>> Error: hypre solver not compatible ',
     &   /16x,'with threads-only execution.'/)
 9500 format(/1x,'>>>>> Warning:  the predict/extrapolate nonlinear',
     &      /,1x,'                solution with a user specified',
     &  ' factor is',
     &      /,1x '                no longer supported. command',
     &     ' ignored...',/)
 9510 format(/1x,'.... dump of line search values ....',
     &      /,10x,'line_search:          ', l1,
     &      /,10x,'ls_details:           ', l1,
     &      /,10x,'ls_min_step_length:   ', f10.4,
     &      /,10x,'ls_max_step_length:   ', f10.4,
     &      /,10x,'ls_rho:               ', f10.4,
     &      /,10x,'ls_slack_tol:         ', f10.4, //)
 9512 format(/1x,'>>>>> error: unknown line search option: ',a,/)
 9514 format(/1x,'>>>>> error: expecting rho value: ',a,/)
 9516 format(/1x,'>>>>> error: expecting min step length value: ',a,/)
 9518 format(/1x,'>>>>> error: expecting max step length value: ',a,/)
 9520 format(/1x,'>>>>> error: expecting slack tolerance value: ',a,/)
c
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine insave_value                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/10/98                   *
c     *                                                              *
c     *     set a value in the crack growth module block             *
c     *                                                              *
c     ****************************************************************
c
      subroutine  insave_value( value, ivalue, action )
c
      use damage_data, only : perm_load_fact
c
      implicit integer (a-z)
      real value, zero, one
      character :: dums
      double precision
     &   dumd
      include 'param_def'
      data zero, one / 0.0, 1.0 /
c
c
      go to ( 100 ), action
c
 100  continue
      if ( value .le. zero ) then
        call errmsg(317,dum,dums,value,dumd)
        return
      end if
      perm_load_fact = 1.0 / value
      return
c
      end



