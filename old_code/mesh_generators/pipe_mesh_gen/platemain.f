c----67--1---------2---------3---------4---------5---------6---------712
c
c    File: platemain.f   Written July 7 1998
c
c    plate1.f is the primary subroutine for the FEACode mesh generation
c
c    Fortran 90 complie with f90
c
c----67--1---------2---------3---------4---------5---------6---------712
c
      program platemain
c
      implicit none
      integer    errorFlag, len_trim
      character  filename*200, jobname*200, adjustl
      logical    ok
c
c
      filename = 'plate3d.in' 
      jobname  = 'test'
c
      write(*,8001) ' '
      write(*,8001) '>> Surface crack mesh generator (focused meshes)'
      write(*,8001) '  '
      write(*,8002) '>> Input file name (no quotes): '
      read(*,8003) filename
      filename = adjustl(filename)
      inquire(file=filename,exist=ok)
      if ( .not. ok ) then
            write(*,8004) filename(1:len_trim(filename))
            stop
      endif
      write(*,8005) filename(1:len_trim(filename))
      write(*,8007)
c
c
      call plate1( errorFlag, filename, jobname )
      write(*,8006) 
c
      stop
c
 8001 format(t1,a)
 8002 format(t1,a,$)
 8003 format(a)
 8004 format('>> Fatal Error.....',
     & /,    '      the specified file: ',a,
     & /,    '      does not exist. program aborted.',//)
 8005 format(/,t5,'> reading input data from the file: ',a)
 8006 format(//,t1,'>> Program completed...',/)
 8007 format(/,t5,'> all generated files have the prefix test_')
c
      end
c
c----67--1---------2---------3---------4---------5---------6---------712
c
