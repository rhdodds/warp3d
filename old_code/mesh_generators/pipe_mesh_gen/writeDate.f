c
c----67--1---------2---------3---------4---------5---------6---------7-2
c Write the date and time to the given file unit number.
c Useful to mark the log files to check that the file is current.
c
      subroutine writeDate(io)
      implicit none
c
      character     day*12,tim*8, a*20,b*20,c*20
      character*4   month(12)
      integer       value(10), io,i,j,k
c
      call date_and_time(a,b,c,value)
c
c  Date
c
      data month /'Jan ','Feb ','Mar ','Apr ',
     &            'May ','Jun ','Jul ','Aug ',
     &            'Sep ','Oct ','Nov ','Dec '/
      i = value(2)
      j = value(3)
      k = value(1)
      if (j.lt.10) then
         write(day,'(a4,a1,i1,a1,i4,a1)') month(i),' ',j,' ',k,' '
      else
         write(day,'(a4,i2,a1,i4,a1)')    month(i),j,' ',k,' '
      endif
      write(io,'(t1,a,a)') 'current date = ',day
c
c  Time
c
      i = value(5)
      j = value(6)
      k = value(7)
      tim = '00:00:00'
      if (i.lt.10) then
         write( tim(2:2), '(i1)' ) i
      else
         write( tim(1:2), '(i2)' ) i
      endif
      if (j.lt.10) then
         write( tim(5:5), '(i1)' ) j
      else
         write( tim(4:5), '(i2)' ) j
      endif
      if (k.lt.10) then
         write( tim(8:8), '(i1)' ) k
      else
         write( tim(7:8), '(i2)' ) k
      endif
      write(io,'(t1,a,a)') 'current time = ',tim
c
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-2
c
