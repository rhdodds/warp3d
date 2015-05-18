c     ****************************************************************
c     *                                                              *
c     *                      subroutine outrac                       *
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
c
      subroutine outrac( step, iter, out, convrg, magdu1, magdp,
     &                   testvl, ratio, 
     &                   adapt_load_fact, res_max, res_max_node,
     &                   avg_force, local_flags, stname,
     &                   show_details, res_max_node_dof,
     &                   mxcvtests )
      implicit none
c
c                        parameters
c     
      integer :: step, iter, out,  res_max_node_dof, res_max_node,
     &           mxcvtests
#dbl      double precision ::
#sgl      real ::
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
      num_conv_tests = mxcvtests
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
 9000 format(/4x,'newton convergence tests step: ',i5,1x,
     &       'iteration: ',i3,1x,'@ wall time: ',f8.1,
     7       /4x,73('-') )
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
 9100 format(/2x,'convergence tests step: ',i6,1x,
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
      implicit integer (a-z)
$add common.main

c
      integer, external :: warp3d_get_device_number
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
 

