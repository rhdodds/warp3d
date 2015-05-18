c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouflnm                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 04/13/2014 rhd             *
c     *                                                              *
c     *     creates a file name for patran output                    *
c     *     based on the output type and the step number.            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouflnm( strtnm, flnam, stepno, use_mpi, myid )
      implicit integer (a-z)
c      
      character *4 strtnm
      character *(*) flnam
      logical use_mpi
c
      step_number = stepno 
c
c                check to make sure the step number is not
c                greater than 5 digits. 
c
      if( step_number .gt. 99999 ) step_number = step_number - 99999
c
      flnam(1:) = ' '             
      flnam(1:4) = strtnm  ! first part of name
      write(flnam(5:),9100) step_number
c
c                for mpi, we add the process id as a suffix
c                to the file name.
c
      if ( .not. use_mpi ) return
      flnam(10:10) = '.'
      write(flnam(11:14),9000) myid
      return
c
 9000 format(i4.4)
 9100 format(i5.5)
c
      end







