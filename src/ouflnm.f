c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouflnm                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 2/27/2018 rhd              *
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
      integer :: step_number  ! can be 7 digits
c
      step_number = stepno
      if( len(flnam) < 16 ) then
        write(*,9200)
        call die_gracefully
      end if
c
c                check to make sure the step number is not
c                greater than 7 digits.
c
      if( step_number .gt. 9999999 )
     &     step_number = step_number - 9999999
c
      flnam(1:) = ' '
      flnam(1:4) = strtnm  ! first part of name
      write(flnam(5:),9100) step_number
c
c                for mpi, we add the process id as a suffix
c                to the file name.
c
      if ( .not. use_mpi ) return
      flnam(12:12) = '.'
      write(flnam(13:16),9000) myid
      return
c
 9000 format(i4.4)
 9100 format(i7.7)
 9200 format('>>>> FATAL ERROR: invlaid string length. ',
     &    'ouflnm.')
c
      end







