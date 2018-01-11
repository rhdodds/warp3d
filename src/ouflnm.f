c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouflnm                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 1/11/2018 rhd              *
c     *                                                              *
c     *     creates a file name for patran output                    *
c     *     based on the output type and the step number.            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouflnm( strtnm, flnam, stepno, use_mpi, myid )
      implicit none
c
      integer :: stepno, myid
      character(len=4) :: strtnm
      character(len=*) :: flnam
      logical :: use_mpi
c
      integer :: step_number
c
      step_number = stepno
c
c                check to make sure the step number is not
c                greater than 6 digits.
c
      if( step_number .gt. 999999 )
     &     step_number = step_number - 999999
c
      flnam(1:) = ' '
      flnam(1:4) = strtnm  ! first part of name
      write(flnam(5:),9100) step_number
c
c                for mpi, we add the process id as a suffix
c                to the file name.
c
      if ( .not. use_mpi ) return
      flnam(11:11) = '.'
      write(flnam(12:15),9000) myid
      return
c
 9000 format(i4.4)
 9100 format(i6.6)
c
      end







