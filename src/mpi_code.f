c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_init                    *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/2/98                    *
c     *                                                              *
c     *         This version is for serial version of warp3d.        *
c     *         This sets the id of the only processor to 0, and     *
c     *         sets the number of processors to 1.                  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_init
      use global_data
      implicit integer (a-z)
c
      myid = 0
      numprocs = 1
      root_processor = .true.
      worker_processor = .false.
      use_mpi = .false.
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_wait                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 10/29/2017 rhd             *
c     *                                                              *
c     *         This is a dummy routine that just returns.           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_wait
      implicit integer (a-z)
c
      return
      end

c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_suspend                 *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/12/99                   *
c     *                                                              *
c     *         This is a dummy routine that just returns.           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_suspend(option)
      implicit integer (a-z)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_alert_slaves            *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/05/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_alert_slaves ( do_it )
      implicit integer (a-z)
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_reduce_vec              *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/05/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_reduce_vec ( in, size )
      implicit integer (a-z)
c
      double precision   in(*)
c
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine wmpi_reduce_vec_std              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 6/17/2017 rhd              *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_reduce_vec_std ( in, size )
      implicit integer (a-z)
c
      double precision   in(*)
c
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_dotprod                 *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/05/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_dotprod ( in )
      implicit integer (a-z)
c
      double precision    in
c
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_redint                  *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/05/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_redint ( int )
      implicit integer (a-z)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_redlog                  *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/05/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_redlog ( log_var )
      implicit integer (a-z)
c
      logical log_var
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_bcast_int               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/05/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_bcast_int (int_var)
      implicit integer (a-z)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_real               *
c     *                                                              *
c     *                       written by : mcw                       *
c     *                                                              *
c     *                   last modified : 02/17/05                   *
c     *                                                              *
c     *      Broadcast a real from slave 'source' to root.           *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_real ( real_vec, size )
      implicit integer (a-z)
      real real_vec(*)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_bcast_log               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/05/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_bcast_log (log_var)
      implicit integer (a-z)
c
      logical log_var
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_do_external_db          *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 4/16/2016 rhd              *
c     *                                                              *
c     *          invoke  UEXTERNALDB Abaqus support routine and      *
c     *          other specific UEXTERNADB routines                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_do_uexternaldb
      use global_data
      implicit integer (a-z)
      double precision ::  zero, aba_time(2), aba_dtime
      logical local_debug
      data zero / 0.0d00 /
c
      local_debug = .false.
c
      select case ( douextdb )
       case ( 1 )   ! from incomp start of step 1
         aba_lop = 0
         aba_lrestart = 0
         aba_time(1) = zero; aba_time(2) = zero
         aba_dtime = dt
         aba_kstep = 1
         aba_kinc = 1
       case ( 2 )  ! from driver (normal termination)
        aba_lop = 3
        aba_lrestart = 0
        aba_time(1) = total_model_time
        aba_time(2) = total_model_time
        aba_dtime = dt
        aba_kstep = 1
        aba_kinc  = ltmstp   ! current WARP3D load step
       case ( 3 )  ! from mnralg (start of step, restart adaptive step)
        aba_lop = 1
        aba_lrestart = 0
        aba_time(1) = total_model_time
        aba_time(2) = total_model_time
        aba_dtime = dt
        aba_kstep = 1
        aba_kinc  = ltmstp + 1  ! current WARP3D load step
       case ( 4 ) ! from mnralg (end step or adaptive substep)
        aba_lop = 2
        aba_lrestart = 0
        aba_time(1) = total_model_time
        aba_time(2) = total_model_time
        aba_dtime = dt
        aba_kstep = 1
        aba_kinc  = ltmstp + 1  ! current WARP3D load step
       case ( 5 )   ! from store (make restart file)
        aba_lop = 2
        aba_lrestart = 1
        aba_time(1) = total_model_time
        aba_time(2) = total_model_time
        aba_dtime = dt
        aba_kstep = 1
        aba_kinc  = ltmstp   ! current WARP3D load step
       case ( 6 )   ! from reopen (restarting analysis)
        aba_lop = 4
        aba_lrestart = 0
        aba_time(1) = total_model_time
        aba_time(2) = total_model_time
        aba_dtime = dt
        aba_kstep = 1
        aba_kinc  = ltmstp + 1  ! next WARP3D load step to be solved
       case default
           write(out,9000)
           call die_abort
       end select
c
      call uexternaldb( aba_lop, aba_lrestart, aba_time,
     &                  aba_dtime, aba_kstep, aba_kinc )
c
      call uexternaldb_mm04_cavity( aba_lop, aba_lrestart, aba_time,
     &                  aba_dtime, aba_kstep, aba_kinc )
      return
 9000 format(/,2x,"FATAL ERROR: invalid douextdb in ",
     & "wmpi_do_uexternaldb. ",/,2x,"Analysis terminated...",//)
c
       end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_basic              *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/27/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_basic
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_analysis           *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/27/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_analysis
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_const              *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/27/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_const
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_iter1              *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/27/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_iter1
      implicit integer (a-z)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_itern              *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/27/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_itern
      implicit integer (a-z)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_step               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 02/27/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_step
      implicit integer (a-z)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *            subroutine wmpi_send_temp_eqloads                 *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 06/15/02                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_temp_eqloads
      implicit integer (a-z)
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                   subroutine wmpi_send_contact               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 07/20/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_contact (ldum)
      logical ldum
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_growth             *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 01/07/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_growth ( killed_this_time )
      logical killed_this_time
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_growth_init             *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 01/07/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_growth_init
      logical killed_this_time
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_get_grow                *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 01/07/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_get_grow
      logical killed_this_time
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_calc_dist               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 01/23/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_calc_dist
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine die_gracefully               *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 03/02/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine die_gracefully
c
      stop
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine die_abort                    *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 10/8/2017 rhd              *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine die_abort
c
      use ifcore
      call tracebackqq('... warp3d stop ...',1)
      stop
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_combine_stf             *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 01/21/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_combine_stf
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_get_str                 *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 01/21/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_get_str (int_var)
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_reopen             *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 01/21/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_reopen
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_scatter_pcm             *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 01/30/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_scatter_pcm
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_jint               *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 01/21/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_send_jint
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_contact_gthr            *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 09/04/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_contact_gthr
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_bcast_string            *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 05/10/99                   *
c     *                                                              *
c     *      Broadcast a logical variable from root to all the slave *
c     *      processors.                                             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_bcast_string ( string, nchars )
      implicit integer (a-z)
      character(len=*) :: string
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_set_mlk_threads         *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 05/31/11                   *
c     *                                                              *
c     *     Has all MPI ranks set the number of MKL threads          *
c     *                                                              *
c     ****************************************************************
      subroutine wmpi_set_mkl_threads(num)
            implicit none
            integer :: num
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_set_omp_threads         *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 05/31/11                   *
c     *                                                              *
c     *     Has all MPI ranks set the number of OMP threads          *
c     *                                                              *
c     ****************************************************************
      subroutine wmpi_set_omp_threads(num)
            implicit none
            integer :: num
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_real_new           *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 06/11                      *
c     *                                                              *
c     *     Get the non-zero data in a vector of a given length      *
c     *     onto the root processor                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_real_new(vector,length)
      implicit none
      real :: vector(*)
      integer :: length
c
      return
      end subroutine
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_simple_angles      *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 3/12/14                    *
c     *                                                              *
c     *           Send the simple angle properties, if required      *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_simple_angles
      implicit none
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_send_crystals           *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 8/16/12                    *
c     *                                                              *
c     *           Send all the crystal properties to the workers     *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_send_crystals
      implicit none
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_dealloc_crystals        *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 8/16/12                    *
c     *                                                              *
c     *           Dealloc crystal data structures                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine wmpi_dealloc_crystals
      implicit none
      return
      end subroutine

c
c     ****************************************************************
c     *                                                              *
c     *                wmpi_compute_set_history_locs                 *
c     *                                                              *
c     *                       written by : tjt                       *
c     *                                                              *
c     *                   last modified : 5/12/2017 tjt              *
c     *                                                              *
c     ****************************************************************
c
      subroutine  wmpi_compute_set_history_locs
      implicit none
      return
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_init_owner              *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 5/2/2013 rhd moved here    *
c     *                                                              *
c     *       Dummy routine for serial version of warp3d.            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_init_owner
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine wmpi_chknode                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 5/2/2013 moved here        *
c     *                                                              *
c     *           The processor assesses whether the node number     *
c     *           passed in as an argument is accessed by this       *
c     *           processor. In threads only execution               *
c     *           it is always true.                                 *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_chknode ( node, referenced, owned )
      implicit integer (a-z)
      logical referenced, owned
c
      referenced = .true.
      owned = .true.
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *              drive the symmetric CPardiso solution           *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/7/20167 rhd              *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine cpardiso_symmetric( neq, ncoeff, k_diag, rhs,
     &   solution_vec, eqn_coeffs, k_pointers, k_indices,
     &   print_cpu_stats, itype, out, myrank )
      implicit integer (a-z)
      return
      end
c     ****************************************************************
c     *                                                              *
c     *              drive the unsymmetric CPardiso solution         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 1/7/20167 rhd              *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine cpardiso_unsymmetric( neqns, nnz, k_ptrs, k_indexes,
     &            k_coeffs,  p_vec, u_vec, cpu_stats, itype, out,
     &            myid )
      implicit integer (a-z)
      return
      end
