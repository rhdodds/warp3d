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
      implicit integer (a-z)
$add common.main
c      
      myid = 0
      numprocs = 1
      root_processor = .true.
      slave_processor = .false.
      use_mpi = .false.
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
$add common.main
c      
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_handle_slaves           *
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
      subroutine wmpi_handle_slaves
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
#dbl      double precision
#sgl      real
     &     in(*)
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
#dbl      double precision
#sgl      real
     &     in 
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
c     *                   last modified : 03/02/98                   *
c     *                                                              *
c     *         Dummy routine for serial version of warp3d.          *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine die_abort
c
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
      character *(*) string
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
c                 Dummy
            real :: vector(*)
            integer :: length
c
            return
      end subroutine
