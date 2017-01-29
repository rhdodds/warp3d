      implicit integer (a-z)
      character * 50 file1_name, file2_name
      real results1(30,100000), results2(30,100000) 
      character * 1 title1_formatted(80), title2_formatted(80)
      real dum2
c
      write(*,*) '>> ascii file 1: '
      read(*,*) file1_name
      write(*,*) '>> ascii file 2: '
      read(*,*) file2_name
c
      file_no1 = 20
      file_no2 = 21
c
       open( unit=file_no1, file=file1_name, status='old', 
     &       access='sequential', form='formatted',
     &       recl=350 )
       open( unit=file_no2, file=file2_name, status='old', 
     &       access='sequential', form='formatted',
     &       recl=350 )
c
      write(*,*) '>> files opened'
c
      read(unit=file_no1,fmt=9300) title1_formatted
      read(unit=file_no1,fmt=9315) num_nodes, dum1, dum2, dum3, 
     &  num_values
      read(unit=file_no1,fmt=9300) title2_formatted
      read(unit=file_no1,fmt=9300) title2_formatted
      write(*,*) '> header records read, file 1...'
c
      read(unit=file_no2,fmt=9300) title1_formatted
      read(unit=file_no2,fmt=9315) num_nodes, dum1, dum2, dum3, 
     &         num_values
      read(unit=file_no2,fmt=9300) title2_formatted
      read(unit=file_no2,fmt=9300) title2_formatted
      write(*,*) '> header records read, file 2...'
c
      write(*,*) '> start node processing...'

      do node = 1, num_nodes
        read(file_no1,fmt=9330)  nod, results1(1:num_values,node)
        read(file_no2,fmt=9330)  nod, results2(1:num_values,node)
c        write(*,*) '> compare node: ',node
        do j = 1, num_values
         value1 = results1(j,node)
         value2 = results2(j,node)
         if ( abs(value1) .le. 1e-10 ) value1 = 0.0
         if ( abs(value2) .le. 1e-10 ) value2 = 0.0
         if ( value1 .ne. value2 ) then
           write(*,*) '>> diff values. node, comp. :',node,j
           write(*,9000) value1, value2
         end if
        end do
      end do
c
c
      write(*,*) '>> all done'
      stop
 9300 format(80a1)
 9315 format(2i9,e15.6,2i9)
 9000 format(20x,e13.6,2x,e13.6)
 9330 format(i8,(5e13.6))
      end
