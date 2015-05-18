c     ****************************************************************
c     *                                                              *
c     *                      module performance_data                 *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 06/10/2013 rhd             *
c     *                                                              *
c     *      stores profiling data (right now, just the wall time)   *
c     *                                                              *
c     ****************************************************************

      module performance_data
      implicit none
c
      double precision, save :: start_wall_time
      logical :: time_assembly
      double precision, save :: start_assembly_step, assembly_total
      integer, save :: ntimes_assembly
c
      contains
            subroutine t_start_assembly(tstart)
                  implicit none
                  double precision :: omp_get_wtime, tstart

                  tstart = omp_get_wtime()
c
            end subroutine
c
            subroutine t_end_assembly(tstore,tstart)
                  implicit none
                  double precision :: omp_get_wtime,tstore,tstart

                  tstore = tstore + (omp_get_wtime()-
     &                   tstart)
            end subroutine
c
      end module performance_data
