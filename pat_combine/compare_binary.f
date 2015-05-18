      implicit integer (a-z)
      character * 100 file1_name, file2_name
      real results1(30,100000), results2(30,100000) 
      character * 4 title1_binary(80), title2_binary(80),
     &              titlet1_binary(80)
      real dum2
c
      write(*,*) '>> binary file 1: '
      read(*,*) file1_name
      write(*,*) '>> binary file 2: '
      read(*,*) file2_name
c
      file_no1 = 20
      file_no2 = 21
c
       open( unit=file_no1, file=file1_name, status='old', 
     &       access='sequential', form='unformatted')
c     &       recl=350 )
       open( unit=file_no2, file=file2_name, status='old', 
     &       access='sequential', form='unformatted')
c     &       recl=350 )
c
      write(*,*) '>> files opened'
c
      read(unit=file_no1) title1_binary,
     &                 num_nodes, dum1, dum2, dum3, num_values1 
      read(unit=file_no1) 
      read(unit=file_no1) 
      write(*,*) '> header records read, file 1...'
c
      read(unit=file_no2) title1_binary,
     &                 num_nodes, dum1, dum2, dum3, num_values2
      read(unit=file_no2) 
      read(unit=file_no2) 
      write(*,*) '> header records read, file 2...'
c
      if ( num_values1 .ne. num_values2 ) then
         write(*,*) '>> num_values do not agree'
         stop
      end if
      num_values = num_values1
       
      write(*,*) '> start node processing...'

      do node = 1, num_nodes
        read(file_no1) nod, results1(1:num_values,node)
        read(file_no2) nod, results2(1:num_values,node)
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
 9000 format(20x,e13.6,2x,e13.6)
      end
