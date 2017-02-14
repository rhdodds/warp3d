      program compare_patran_element 
c      
      implicit none
      integer :: out, num_chars_to_compare, file_no1, file_no2,
     &           header_written, stat1, stat2, num_elements,
     &           num_values1, num_values2, dum1, dum3, j,
     &           num_values, num_vals_to_compare, idummy,
     &           len1, len2, elem1, elem2, elem, last_slash
      logical :: bad1, bad2, binary, formatted, states1, states2,
     &           states 
      real, allocatable :: results1(:,:), results2(:,:)
      real :: value1, value2, dum2, zero_tol
      character(len=80) :: file1_name, file2_name, form, matl1, matl2
      character(len=1) :: title1_formatted(80), title2_formatted(80)
      character(len=4) :: title1_binary(80), title2_binary(80)
      character(len=9) :: root1, root2
      character(len=20) charval1, charval2
c
      out = 6
      num_chars_to_compare = 5
      num_vals_to_compare  = 6
      file_no1 = 20
      file_no2 = 21
c      
      write(out,*) " "
      write(out,*) ".. .Compare Patran Element Results"
      write(out,*) " "
      write(out,*) '>> Patran element results file 1: '
      read(*,*) file1_name
      write(out,*) '>> Patran element results file 2: '
      read(*,*) file2_name 
      write(out,*) 
     & '>> # digits to compare, # values to compare each node:'
      read(*,*) num_elements, num_chars_to_compare,
     &          num_vals_to_compare
      write(out,*) 
     & '>> set values smaller than this to zero:'
      read(*,*) zero_tol
c
      bad1 = num_elements <= 0 .or.
     &       num_chars_to_compare <= 0  .or.
     &       num_vals_to_compare <= 0   .or.
     &       zero_tol <= 0.0  
      if( bad1 ) then
        write(out,9120) num_elements, num_chars_to_compare,
     &                  num_vals_to_compare, zero_tol
        stop
      end if  
c
c              get file name w/o directories
c
      len1 = len_trim( file1_name )
      last_slash = index( file1_name, "/", .true. )
      if( last_slash == 0 ) then
         root1 = file1_name
      else
         root1 = file1_name(last_slash+1:)
      end if
      last_slash = index( file2_name, "/", .true. )
      if( last_slash == 0 ) then
         root2 = file2_name
      else
         root2 = file2_name(last_slash+1:)
      end if

      bad1 = (root1(2:2) .ne. 'e')  .or. (root2(2:2) .ne. 'e')
      if( bad1 ) then
        write(out,9130)
        stop
      end if
c      
      binary    = (root1(3:3) == 'b') .and. 
     &            (root2(3:3) == 'b')
      formatted = (root1(3:3) == 'f') .and.
     &            (root2(3:3) == 'f')
      if( .not. (binary .or. formatted ) ) then
        write(out,9140)
        stop
      end if

c
      states1 = root1(4:4) == 'm'
      states2 = root2(4:4) == 'm'
      states  = states1 .and. states2
      if( states ) then ! both files must be same material
       matl1 = root1(11:)
       matl2 = root2(11:)
       if( matl1 .ne. matl2 ) then
        write(out,9150)
        stop
       end if
      end if 

c    
      form = "formatted"
      if( binary ) form = "unformatted"
      open( unit=file_no1, file=file1_name, status='old', iostat=stat1,
     &      form =form   )
      open( unit=file_no2, file=file2_name, status='old', iostat=stat2,
     &      form=form   )
      if( stat1 .ne. 0 ) then 
        write(out,*) '>> FATAL. file not found: ',
     &      file1_name(1:len_trim(file1_name))
        stop
      end if
      if( stat2 .ne. 0 ) then
        write(out,*) '>> FATAL. file not found: ',
     &      file2_name(1:len_trim(file2_name))
        stop
      end if

c
      write(out,9160) file1_name(1:len_trim(file1_name)),
     &                file2_name(1:len_trim(file2_name))
      write(out,*) '>> files opened'
c
      if( formatted ) then
         read(unit=file_no1,fmt=9300) title1_formatted
         read(unit=file_no1,fmt="(i)") num_values1
         read(unit=file_no1,fmt=9300) title1_formatted
         read(unit=file_no1,fmt=9300) title1_formatted
         write(out,*) '> header records read, file 1...'
c   
         read(unit=file_no2,fmt=9300) title2_formatted
         read(unit=file_no2,fmt="(i)") num_values2
         read(unit=file_no2,fmt=9300) title2_formatted
         read(unit=file_no2,fmt=9300) title2_formatted
         write(out,*) '> header records read, file 2...'
      end if

      if( binary ) then 
         read(unit=file_no1) title1_binary, num_values1 
         read(unit=file_no1) 
         read(unit=file_no1) 
         write(out,*) '> header records read, file 1...'
         read(unit=file_no2) title2_binary, num_values2 
         read(unit=file_no2) 
         read(unit=file_no2) 
         write(out,*) '> header records read, file 2...'
      end if
c

      bad2 = num_values1 .ne. num_values2    
      if( bad2 ) then
        write(out,9100) num_values1, num_values2
        stop
      end if

c
      num_values = num_values1      
      allocate( results1(num_values,num_elements), 
     &          results2(num_values,num_elements) )
c
      write(out,*) '> start elements processing...'
      header_written = .false.
c
      do elem = 1, num_elements
       if( mod(elem,100) == 0 )
     &         write(out,*) "... comparing element: ", elem
        if( formatted ) then
          read(unit=file_no1,fmt="(2i)")  elem1, idummy
          read(unit=file_no1,fmt=9330) results1(1:num_values,elem)
          read(unit=file_no2,fmt="(2i)")  elem2, idummy
          read(unit=file_no2,fmt=9330) results2(1:num_values,elem)
        end if
        if( binary ) then
          read(file_no1)  elem1, idummy, results1(1:num_values,elem)
          read(file_no2)  elem2, idummy, results2(1:num_values,elem)
        end if
        if( elem1 .ne. elem2 ) then
          write(out,9110) elem1, elem2
          stop
        end if
        do j = 1, num_vals_to_compare
         value1 = results1(j,elem)
         value2 = results2(j,elem)
         if( abs(value1) .le. zero_tol ) value1 = 0.0
         if( abs(value2) .le. zero_tol ) value2 = 0.0
         write(charval1,9210) value1
         write(charval2,9210) value2
         if( charval1(1:num_chars_to_compare) .ne.
     &       charval2(1:num_chars_to_compare) ) then
           if( .not. header_written ) write(out,9200)
           header_written = .true.
           write(out,9010) elem, j, value1, value2
         end if
        end do
      end do
c
      write(*,*) '>> all done'
      stop
c
 9000 format(20x,e13.6,2x,e13.6)
 9010 format(20x,i8,i4,2x,e13.6,2x,e13.6)
 9100 format('>> FATAL: inconsistent file 1 and 2 data:',
     & /,10x,'...num_values1, num_values2: ',2i10)      
 9110 format('>> FATAL: inconsistent file 1 and 2 data:',
     & /,10x,'...elem1, elem2: ',2i10 )     
 9120 format('>> FATAL: inconsistent input:',
     & /,10x,'...num_elements: ',i10,
     & /,10x,'...num_chars_to_compare: ',i10,
     & /,10x,'...num_vals_to_compare:  ',i10,
     & /,10x,'...zero_tol:             ',e14.6 )
 9130 format('>> FATAL: not element result files')
 9140 format('>> FATAL: inconsistent input:',
     & /,10x,'...both files must be formatted or binary')
 9150 format('>> FATAL: inconsistent input:',
     & /,10x,'...both states files must be same material model')
 9160 format('...comparing 2 files: ',
     & /10x,a,/,10x,a)
 9200 format('... mismatch @ element    component, 2 values')
 9210 format(e14.6) 
 9300 format(80a1)
 9315 format(2i9,e15.6,2i9)
 9330 format(6e13.6)

c 
      end
