      program compare_patran_nodal
c      
      implicit none
      integer :: out, num_chars_to_compare, file_no1, file_no2,
     &           header_written, stat1, stat2, num_nodes1, num_nodes2,
     &           num_values1, num_values2, dum1, dum3, nod1, nod2, j,
     &           node, num_nodes, num_values, num_vals_to_compare,
     &           len1, len2
      integer, parameter :: max_nodes = 5000000
      logical :: bad1, bad2, binary, formatted
      real, allocatable :: results1(:,:), results2(:,:)
      real :: value1, value2, dum2, zero_tol
      character(len=80) :: file1_name, file2_name, form
      character(len=1) :: title1_formatted(80), title2_formatted(80)
      character(len=4) :: title1_binary(80), title2_binary(80)
      character(len=9) :: root1, root2
      character(len=20) charval1, charval2
c
      allocate( results1(50,max_nodes), results2(50,max_nodes) )
      out = 6
      write(out,*) " "
      write(out,*) ".. .Compare Patran Nodal Results"
      write(out,*) " "
      write(out,*) '>> Patran node results file 1: '
      read(*,*) file1_name
      write(out,*) '>> Patran node results file 2: '
      read(*,*) file2_name 
      write(out,*) 
     & '>> # digits to compare, # values to compare each node:'
      read(*,*) num_chars_to_compare, num_vals_to_compare
      write(out,*) 
     & '>> set values smaller than this to zero:'
      read(*,*) zero_tol
c
      bad1 = num_chars_to_compare <= 0  .or.
     &       num_vals_to_compare <= 0   .or.
     &       zero_tol <= 0.0  
      if( bad1 ) then
        write(out,9120) num_chars_to_compare, num_vals_to_compare,
     &                  zero_tol
        stop
      end if  
c
      num_chars_to_compare = 5
      num_vals_to_compare = 6
      file_no1 = 20
      file_no2 = 21
      len1 = len_trim( file1_name )
      len2 = len_trim( file2_name )
      root1 = file1_name(len1-8:len1)
      root2 = file2_name(len2-8:len2)
c      
      bad1 = (root1(2:2) .ne. 'n')  .or. (root2(2:2) .ne. 'n')
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
      write(out,9150) file1_name(1:len_trim(file1_name)),
     &                file2_name(1:len_trim(file2_name))
      write(out,*) '>> files opened'
c
      if( formatted ) then
         read(unit=file_no1,fmt=9300) title1_formatted
         read(unit=file_no1,fmt=9315) num_nodes1, dum1, dum2, dum3, 
     &                                num_values1
         read(unit=file_no1,fmt=9300) title2_formatted
         read(unit=file_no1,fmt=9300) title2_formatted
         write(out,*) '> header records read, file 1...'
c   
         read(unit=file_no2,fmt=9300) title1_formatted
         read(unit=file_no2,fmt=9315) num_nodes2, dum1, dum2, dum3, 
     &                                num_values2
         read(unit=file_no2,fmt=9300) title2_formatted
         read(unit=file_no2,fmt=9300) title2_formatted
         write(out,*) '> header records read, file 2...'
      end if

      if( binary ) then 
         read(unit=file_no1) title1_binary, num_nodes1, dum1, dum2,
     &                       dum3, num_values1 
         read(unit=file_no1) 
         read(unit=file_no1) 
         write(out,*) '> header records read, file 1...'
         read(unit=file_no2) title2_binary, num_nodes2, dum1, dum2,
     &                       dum3, num_values2
         read(unit=file_no2) 
         read(unit=file_no2) 
         write(out,*) '> header records read, file 2...'
      end if
c

      bad1 = num_nodes1 .ne. num_nodes2 
      bad2 = num_values1 .ne. num_values2    
      if( bad1 .or. bad2 ) then
        write(out,9100) num_nodes1, num_nodes2, 
     &                  num_values1, num_values2
        stop
      end if
c
      num_nodes  = num_nodes1
      num_values = num_values1      
c
      write(out,*) '> start node processing...'
      header_written = .false.

      do node = 1, num_nodes
       if( mod(node,100) == 0 )
     &         write(out,*) "... comparing node: ", node
        if( formatted ) then
          read(file_no1,fmt=9330)  nod1, results1(1:num_values,node)
          read(file_no2,fmt=9330)  nod2, results2(1:num_values,node)
        end if
        if( binary ) then
          read(file_no1)  nod1, results1(1:num_values,node)
          read(file_no2)  nod2, results2(1:num_values,node)
        end if
        if( nod1 .ne. nod2 ) then
          write(out,9110) nod1, nod2
          stop
        end if
        do j = 1, num_vals_to_compare
         value1 = results1(j,node)
         value2 = results2(j,node)
         if( abs(value1) .le. zero_tol ) value1 = 0.0
         if( abs(value2) .le. zero_tol ) value2 = 0.0
         write(charval1,9210) value1
         write(charval2,9210) value2
         if( charval1(1:num_chars_to_compare) .ne.
     &       charval2(1:num_chars_to_compare) ) then
           if( .not. header_written ) write(out,9200)
           header_written = .true.
           write(out,9010) node, j, value1, value2
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
     & /,10x,'...num_nodes1, num_nodes2:   ',2i10,
     & /,10x,'...num_values1, num_values2: ',2i10)      
 9110 format('>> FATAL: inconsistent file 1 and 2 data:',
     & /,10x,'...nod1, nod2: ',2i10 )     
 9120 format('>> FATAL: inconsistent input:',
     & /,10x,'...num_chars_to_compare: ',i10,
     & /,10x,'...num_vals_to_compare:  ',i10,
     & /,10x,'...zero_tol:             ',e14.6 )
 9130 format('>> FATAL: not nodal result files')
 9140 format('>> FATAL: inconsistent input:',
     & /,10x,'...both files must be formatted or binary')
 9150 format('...comparing 2 files: ',
     & /10x,a,/,10x,a)
 9200 format('... mismatch @ node     component, 2 values')
 9210 format(e14.6) 
 9300 format(80a1)
 9315 format(2i9,e15.6,2i9)
 9330 format(i8,(5e13.6))
c 
      end
