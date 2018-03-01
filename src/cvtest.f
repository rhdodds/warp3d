c     ****************************************************************
c     *                                                              *
c     *                      subroutine cvtest                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 1/11/2018 rhd              *
c     *                                                              *
c     *     perform convergence tests specified                      *
c     *     by the user. if any one of the convergence criteria is   *
c     *     met, the convergence flag is set. solution was not       *
c     *     allowed to proceed w/o at least one test defined         *
c     *                                                              *
c     *     set flag if solution appears to be diverging             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine cvtest( cnverg, step, iter, magdu1, mgload,
     &                   adapt_load_fact, diverging_flag )
      use global_data ! old common.main
      use mod_mpc, only : tied_con_mpcs_constructed, mpcs_exist
      use stiffness_data, only : d_lagrange_forces
      use main_data, only : diverge_check_strict
      implicit integer (a-z)
      double precision
     & magdu1, mgload, adapt_load_fact, res_max
      logical cnverg, mducom, trcsol, diverging_flag
c
c                       locals
c
      integer, parameter :: num_local_flags = 15
      double precision ::
     &  testvl(mxcvtests), ratio(mxcvtests), absval, zero,
     &  hundred, avg_force,
     &  term_count, sum, val_1, val_2, val_3, resforce
      double precision, save ::
     &  norm_resids_for_checking(50)
      logical local_flags(num_local_flags), have_mpc_equations
      data zero, hundred / 0.0d00, 100.0d00 /
c
      have_mpc_equations = tied_con_mpcs_constructed .or. mpcs_exist
c
      max_tests = mxcvtests  ! but only 6 tests current available
      avg_force = zero
c
c                       compute norm of residual forces. compare
c                       with values from prior iterations. examine user
c                       specified conditions to adapt no if solution
c                       appears to be diverging.
c                       note save attribute for vector of prior norms
c                       terms in res() for absolute constraints are
c                       zero at this point.
c
c                       do not include Lagrange MPC forces in making
c                       convergence tests
c
      if( iter .eq. 1 ) norm_resids_for_checking(1:50) = zero
      res_max   = -hundred
      sum = zero
      do j = 1, nodof
        resforce = res(j)
        if( have_mpc_equations ) then
          if( d_lagrange_forces(j) .ne. zero ) resforce = zero
        end if
        sum = sum + resforce**2
        if ( abs(resforce) .gt. res_max ) then
          res_max     = abs(resforce)
          res_max_dof = j
        end if
      end do
      sum = sqrt( sum )
      res_max      = res(res_max_dof)
      res_max_node = (res_max_dof-1)/3 + 1
      res_max_node_dof = res_max_dof - (res_max_node-1)*3
c
      norm_resids_for_checking(iter) = sum
      diverging_flag = .false.
      if( diverge_check_strict .and. iter .ge. 2 ) then
        val_2 = norm_resids_for_checking(iter)
        val_1 = norm_resids_for_checking(iter-1)
        diverging_flag = val_2 .ge. val_1
      elseif( iter .ge. 3 ) then
        val_3 = norm_resids_for_checking(iter)
        val_2 = norm_resids_for_checking(iter-1)
        val_1 = norm_resids_for_checking(iter-2)
        diverging_flag = (val_3 .ge. val_2) .and. (val_2 .ge. val_1)
      end if
c
c                       set flags.
c
      trcsol = trace(1)
      mducom = .false.
      local_flags(1:num_local_flags)  = .false.
c
c                       loop over all possible tests.
c
      do i = 1, max_tests
       if( .not. convrg(i) ) cycle
       select case ( i )
c
       case( 1 )
c      =========
c
c                  ****    test 1    ****
c
c                  test compares the ratio of the norm of the
c                  corrective displacement vector to the norm of
c                  the corrective displacement vector for the
c                  first iteration of the current step.
c
c                  the norm of the corrective displacement vector for
c                  the first iteration of the current step is needed.
c                  if the iteration number is one, compute it and
c                  set a flag indicating that it has been computed.
c
        if( iter .eq. 1 ) then
            magdu1  = zero
            do j = 1, nodof
               magdu1 = magdu1 + idu(j)**2
            end do
            magdu1 = sqrt(magdu1)
            mducom = .true.
c
c                  test to see if convergence has been met. since
c                  it is the first iteration, convergence will be
c                  met only if tol is greater than or equal to 100
c
            if( tol(1) .gt. hundred) local_flags(1) = .true.
c
         else
c
c                  subsequent iteration. compute the magnitude of
c                  the corrective displacement vector and compare
c                  it to magdu1.
c
            testvl(1) = zero
            do j = 1, nodof
               testvl(1) = testvl(1) + idu(j)**2
            end do
            testvl(1) = sqrt(testvl(1))
c
c                  perform convergence test - avoid zero divide
c
            if ( magdu1 .eq. zero ) then
               ratio(1) = zero
               local_flags(1) = .true.
            else
               ratio(1) = (testvl(1)/magdu1)*hundred
               if( ratio(1) .lt. tol(1) ) local_flags(1) = .true.
            end if
c
         end if
c
       case( 2 )
c      =========
c
c                ****    test 2    ****
c
c
c                test compares the ratio of the norm of the
c                residual load vector to the norm of the current
c                load vector for the step.
c
c                test to see if convergence has been met.
c                compute the magnitude of residual load vector
c                for unconstrained dof only and compare it to
c                mgload. if norm of the residual load vector
c                and the applied load is < 1.0e-08, we assume
c                convergence no matter  what the tolerance.
c
         testvl(2) = sum
         if ( mgload .eq. zero ) then
             ratio(2) = zero
             local_flags(2) = .true.
         else
             ratio(2) = (testvl(2)/mgload)*hundred
             if( ratio(2) .le. tol(2) ) local_flags(2) = .true.
             if ( mgload .le. 1.0e-08 .and.  testvl(2) .le.
     &            1.0e-08 ) then
                   local_flags(2) = .true.
                   local_flags(6) = .true.
             end if
         end if
c
       case( 3 )
c      =========
c
c                 ****    test 3    ****
c
c                 test compares the absolute value of the max-
c                 imum entry in the corrective displacement
c                 vector and the norm of the first iteration
c                 corrective displacement vector.
c
c                 the norm of the corrective displacement vector for
c                 the first iteration of the current step is needed.
c                 if the iteration number is one, and it has not
c                 already been computed, compute it and set a
c                 flag indicating that it has been computed.
c
         if( iter .eq. 1 ) then
            if( .not. mducom ) then
               magdu1 = zero
               do j = 1, nodof
                  magdu1 = magdu1 + idu(j)**2
               end do
               magdu1 = sqrt(magdu1)
               mducom = .true.
            end if
         end if
c
c                 test to see if convergence has been met.
c                 compute the absolute value of the maximum
c                 entry in the corrective displacement vector
c                 and compare it to magdu1.
c
         testvl(3) = zero
         do j = 1, nodof
            absval = abs( idu(j) )
            if( absval .gt. testvl(3) ) testvl(3) = absval
         end do
c
c                 perform convergence test.
c
         if ( magdu1 .eq. zero ) then
           ratio(3) = zero
           local_flags(3) = .true.
         else
           ratio(3) = (testvl(3)/magdu1)*hundred
           if( ratio(3) .lt. tol(3) ) local_flags(3) = .true.
         end if
c
       case( 4 )
c      =========
c
c                 ****    test 4    ****
c
c                 test compares the absolute value of the max-
c                 imum entry in the residual load vector and
c                 the average nodal force in the model
c
         testvl(4)  = abs( res_max )
         term_count = dble( num_term_ifv + num_term_loads )
         avg_force  = (sum_ifv + sum_loads) / term_count
c
c                 perform convergence test.
c
         if ( avg_force .eq. zero ) then
             ratio(4) = zero
             local_flags(4) = .true.
         else
             ratio(4) = (testvl(4)/avg_force)*hundred
             if( ratio(4) .le. tol(4) ) local_flags(4) = .true.
             if ( testvl(4) .le. 1.0e-08 .and. avg_force .le.
     &            1.0e-08 ) then
                 local_flags(4) = .true.
                 local_flags(8) = .true.
             end if
        end if
c
       case( 5 )
c      =========
c
c                 ****    test 5    ****
c
c                 test compares the norm of the
c                 corrective displacement vector to an
c                 absolute value
c
        testvl(5) = zero
        do j = 1, nodof
          testvl(5) = testvl(5) + idu(j)**2
        end do
        testvl(5) = sqrt(testvl(5))
        local_flags(10)  = testvl(5) .le. tol(5)
c
       case( 6 )
c      =========
c
c                 ****    test 6    ****
c
c                 test compares max entry in the
c                 corrective displacement vector to an
c                 absolute value
c
        testvl(6) = zero
        do j = 1, nodof
            absval = abs( idu(j) )
            if( absval .gt. testvl(6) ) testvl(6) = absval
        end do
        local_flags(11)  = testvl(6) .le. tol(6)
c
      end select

      end do  ! over i
c
c                output the results of the convergence tests
c                if desired by the user. all tests performed
c                must be satisfied for acceptance of the
c                solution.
c
      cnverg = .true.
      if ( convrg(1) .and. .not. local_flags(1) ) cnverg = .false.
      if ( convrg(2) .and. .not. local_flags(2) ) cnverg = .false.
      if ( convrg(3) .and. .not. local_flags(3) ) cnverg = .false.
      if ( convrg(4) .and. .not. local_flags(4) ) cnverg = .false.
      if ( convrg(5) .and. .not. local_flags(10) ) cnverg = .false.
      if ( convrg(6) .and. .not. local_flags(11) ) cnverg = .false.
      if( trcsol ) then
         call cv_outrac( step, iter, out, convrg, magdu1,
     &                   mgload, testvl, ratio, adapt_load_fact,
     &                   res_max, res_max_node, avg_force,
     &                   local_flags, stname, show_details,
     &                   res_max_node_dof, max_tests )
      end if
c
      return
c
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine cv_outrac                    *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 05/12/2015 rhd             *
c     *                                                              *
c     *     outputs status of the convergence                        *
c     *     tests for global Newton iterations                       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine cv_outrac( step, iter, out, convrg, magdu1,
     &                      magdp, testvl, ratio, adapt_load_fact,
     &                     res_max, res_max_node, avg_force,
     &                     local_flags, stname, show_details,
     &                     res_max_node_dof, max_tests )
      implicit none
c
c                        parameters
c
      integer :: step, iter, out,  res_max_node_dof, res_max_node,
     &           max_tests
      double precision ::
     & magdu1, magdp, testvl(*), ratio(*), adapt_load_fact,
     & res_max, avg_force
      logical :: convrg(*), local_flags(*),
     &           show_details, write_max_resid
      character (*) :: stname
c
c                        locals
c
      real :: time
      real, external :: wwalltime
      integer :: bmout, test_num, num_conv_tests
      integer, external :: warp3d_get_device_number
      logical, external :: warp3d_batch_message_file_info
      logical :: batch_messages
      character (8) ::  passed, failed, test_result
      character (1) ::  global_dir(3)
      data passed(1:), failed(1:) / '(passed)', '(failed)' /
      data global_dir / "X", "Y", "Z" /
c
c                        write convergence messages to a message
c                        file in addition to the default output
c                        device if user requests the option.
c
      batch_messages = warp3d_batch_message_file_info( bmout )
      time = wwalltime( 1 )
c
      if( show_details ) write(out,9000) step, iter, time
      if( show_details ) write(out,9005) adapt_load_fact
      if( batch_messages ) then
         call error_count( bmout, .true. )
         write(bmout,9100) step, iter, time, adapt_load_fact
         write(bmout,9015) res_max, res_max_node,
     &                     global_dir(res_max_node_dof)
      end if
c
c                       loop over all possible tests.
c
      write_max_resid = .true.
      num_conv_tests = max_tests
c
      do test_num = 1, num_conv_tests
c
      if( .not. convrg(test_num) ) cycle ! test not active
      select case( test_num )
c
c                       test comparing the ratio of the norm of the
c                       corrective displacement vector to the norm of
c                       the displacement vector for the first iteration
c                       of the current step.
c
       case( 1 )
         if( iter .eq. 1 ) then
            ratio(1) = 100.0
            testvl(1) = magdu1
         end if
         if( write_max_resid ) then
            if( show_details ) write(out,9015) res_max, res_max_node,
     &                         global_dir(res_max_node_dof)
            write_max_resid = .false.
         end if
         test_result(1:) = failed
         if( local_flags(1) ) test_result(1:) = passed
         if( show_details ) write(out,9010) testvl(1), magdu1,
     &                     ratio(1),test_result
         if( batch_messages ) then
             write(bmout,9011) testvl(1), magdu1, ratio(1),
     &                         test_result
         end if
         cycle
c
c                       test comparing the ratio of the norm of the
c                       residual load vector to the norm of the
c                       total load vector for the current step.
c
       case( 2 )
         if( write_max_resid ) then
           if( show_details ) write(out,9015) res_max, res_max_node,
     &                        global_dir(res_max_node_dof)
           write_max_resid = .false.
         end if
         test_result(1:) = failed
         if( local_flags(2) ) test_result(1:) = passed
         if( show_details ) write(out,9021) testvl(2), magdp,
     &                   ratio(2), test_result
         if( local_flags(6) .and. show_details ) write(out,9022)
         if( batch_messages ) then
           write(bmout,9020) testvl(2), magdp, ratio(2), test_result
           if( local_flags(6) ) write(bmout,9022)
         end if
         cycle
c
c                       test comparing the absolute value of the max-
c                       imum entry in the corrective displacement
c                       vector and the norm of the first iteration
c                       displacement vector.
c
       case( 3 )
        if( write_max_resid ) then
           if( show_details ) write(out,9015) res_max, res_max_node,
     &                        global_dir(res_max_node_dof)
           write_max_resid = .false.
        end if
        test_result(1:) = failed
        if( local_flags(3) ) test_result(1:) = passed
        if( show_details )  write(out,9030) testvl(3),
     &                  magdu1, ratio(3), test_result
        if( batch_messages ) then
          write(bmout,9031) testvl(3), magdu1, ratio(3), test_result
        end if
        cycle
c
c                       test comparing the absolute value of the max-
c                       imum entry in the residual load vector and
c                       the average of all forces exerted on nodes
c                       (internal forces, applied loads, reactions,
c                        inertia forces)
c
       case( 4 )
         if( write_max_resid ) then
            if( show_details ) write(out,9015) res_max, res_max_node,
     &                        global_dir(res_max_node_dof)
           write_max_resid = .false.
         end if
         test_result(1:) = failed
         if( local_flags(4) ) test_result(1:) = passed
         if( show_details ) write(out,9040) testvl(4),
     &                   avg_force,ratio(4),test_result
         if( local_flags(8) .and. show_details ) write(out,9022)
         if( batch_messages ) then
           write(bmout,9041) testvl(4),avg_force,ratio(4),test_result
           if( local_flags(8) ) write(bmout,9022)
         end if
         cycle
c
c                       test comparing the norm of corrective
c                       displacement vector an absolute value
c
       case( 5 )
         if( write_max_resid ) then
           if( show_details ) write(out,9015) res_max, res_max_node,
     &                        global_dir(res_max_node_dof)
           write_max_resid = .false.
         end if
         test_result(1:) = failed
         if( local_flags(10) ) test_result(1:) = passed
         if( show_details )  write(out,9050) testvl(5), test_result
         if( batch_messages ) write(bmout,9051) testvl(5), test_result
         cycle
c
c                       test comparing the max value corrective
c                       displacement vector an absolute value
c
       case( 6 )
         if( write_max_resid ) then
           if( show_details ) write(out,9015) res_max, res_max_node,
     &                        global_dir(res_max_node_dof)
           write_max_resid = .false.
         end if
         test_result(1:) = failed
         if( local_flags(11) ) test_result(1:) = passed
         if( show_details )  write(out,9060) testvl(6), test_result
         if( batch_messages ) write(bmout,9061) testvl(6), test_result
         cycle
c
       end select
c
      end do
c
      return
c
 9000 format(/4x,'newton convergence tests step: ',i7,1x,
     &       'iteration: ',i3,1x,'@ wall time: ',f8.1,
     7       /4x,75('-') )
 9005 format(4x,'completed fraction over step: ',f7.4)
c
 9010 format(
     & 4x,'test 1: norm of corrective displacement vector:  ',e12.5,
     &/4x,'        norm of displacement vector iteration 1: ',e12.5,
     &/4x,'        ratio*100: ',f15.5,3x,a8/ )

 9011 format(
     & 4x,'test 1: norm corrective displs, displs iter 1:  ',e12.5,
     &    e12.5,' ratio*100: ',f11.5,2x,a8 )
c
 9015 format(
     & 4x,'maximum residual force: ',e14.6,' @ (node,dir): ',i7,1x,a1)
c
 9020 format(
     & 4x,'test 2: norm resid, total loads:',e12.5,e12.5,
     & ' ratio*100: ',f11.5,2x,a8 )
c
 9021 format(
     & 4x,'test 2: norm of residual load vector:         ',e12.5,
     &/4x,'        norm of total load vector:            ',e12.5,
     &/4x,'        ratio*100: ',f15.5,3x,a8/ )

 9022 format(
     & 4x,'        *** convergence assumed due to very',
     & ' small values ***',/)
c
 9030 format(
     & 4x,'test 3: max entry in corrective displacement vector: ',e12.5,
     &/4x,'        norm of displacement vector iteration 1:     ',e12.5,
     &/4x,'        ratio*100: ',f15.5,3x,a8/ )
 9031 format(
     & 4x,'test 3: max corrective displ, norm displ iter 1 :',
     & e12.5,e12.5,' ratio*100: ',f11.5,2x,a8 )
c
 9040 format(
     & 4x,'test 4: max entry in residual load vector:    ',e12.5,
     &/4x,'        average of all nodal forces:          ',e12.5,
     &/4x,'        ratio*100: ',f15.5,3x,a8/ )
 9041 format(
     & 4x,'test 4: max residual, avg force:',
     & e12.5,e12.5,' ratio*100: ',f11.5,2x,a8 )
c
 9050 format(
     & 4x,'test 5: norm of corrective displacement vector:  ',e12.5,
     & 1x,a8/ )
 9051 format(
     & 4x,'test 5: norm of corrective displacement vector:',
     & e12.5,1x,a8 )
c
 9060 format(
     & 4x,'test 6: max term of corrective displacement vector:  ',e12.5,
     & 1x,a8/ )
 9061 format(
     & 4x,'test 6: max term corrective displacement vector:',
     & e12.5,1x,a8 )
c
 9100 format(/2x,'convergence tests step: ',i7,1x,
     &       'iter: ',i2,1x,'wall time: ',f6.0,1x,
     &       'step fraction: ',f7.4 )
c
      end

c     ****************************************************************
c     *                                                              *
c     *   logical function warp3d_batch_message_file_info            *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/5/2014                   *
c     *                                                              *
c     *    if batch messages are on, return file number to use       *
c     *                                                              *
c     ****************************************************************
c
      logical function warp3d_batch_message_file_info( bmout )
      use main_data, only :  batch_mess_fname
      use global_data, only : batch_messages, stname
      implicit none
c
      integer :: bmout
c
      integer, external :: warp3d_get_device_number
      integer :: lastchar, now_unit
      character(len=80) :: batch_file
      logical :: bm_file_exists, now_open
c
      bmout = -99
      warp3d_batch_message_file_info = .false.
      if( .not. batch_messages ) return
c
c             build file name for batch messages or user
c             input specified name
c
      warp3d_batch_message_file_info = .true.
      if( batch_mess_fname(1:1) .eq. " " ) then
         lastchar =  index( stname(1:8), ' ' ) - 1
         if ( lastchar .le. 0 ) lastchar = 8
         batch_file(1:) = stname(1:lastchar) //'.batch_messages'
      else
         batch_file(1:) = batch_mess_fname(1:)
      end if
c
c             if somebody left the message file open, close
c             and re-open to get in sequence
c
      inquire( file=batch_file, opened=now_open,
     &         exist=bm_file_exists, number=now_unit )
      if( now_open .and. bm_file_exists ) then
          close(unit=now_unit)
      end if
c
      bmout  = warp3d_get_device_number()
      if( bm_file_exists ) then
           open(unit=bmout,file=batch_file,
     &        form='formatted',status='old',position='append')
        else
           open(unit=bmout,file=batch_file,
     &        form='formatted',status='replace' )
      end if
c
      return
      end





