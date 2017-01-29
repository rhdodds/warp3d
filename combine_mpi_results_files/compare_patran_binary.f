      implicit integer (a-z)
      integer, parameter :: max = 3000000
      character(len=100) :: file1_name, file2_name
      real  :: results1(30,max), results2(30,max) 
      character(len=4) :: title1_binary(80), title2_binary(80),
     &              titlet1_binary(80)
      character(len=14) charval1, charval2
      real :: dum2, value1, value2
      logical :: header_written, compare_nodal_results,
     &           compare_element_results 
c

      compare_nodal_results    = .false.
      compare_element_results  = .not. compare_nodal_results
c      
      write(*,*) '>> binary file 1: '
c      read(*,*) file1_name
      file1_name = "./webs00001_mpi"
      write(*,*) '>> binary file 2: '
c      read(*,*) file2_name
      file2_name = "./webs00001_threads"
c
      file_no1 = 20
      file_no2 = 21
      out = 6
      num_chars_to_compare = 5
      num_elements = 109120   ! required to compare element results
c
      open( unit=file_no1, file=file1_name, status='old', 
     &       access='sequential', form='unformatted')
c     &       recl=350 )
       open( unit=file_no2, file=file2_name, status='old', 
     &       access='sequential', form='unformatted')
c     &       recl=350 )
c
      write(out,*) '>> files opened'
c
      if( compare_nodal_results ) then 
         read(unit=file_no1) title1_binary,
     &                 num_rows1, dum1, dum2, dum3, num_values1 
         read(unit=file_no1) 
         read(unit=file_no1) 
         write(out,*) '> header records read, file 1...'
      else
         read(unit=file_no1) title1_binary, num_values1 
         read(unit=file_no1) 
         read(unit=file_no1) 
         write(out,*) '> header records read, file 1...'
         num_rows1 = num_elements
      end if    
c
      if( compare_nodal_results ) then 
         read(unit=file_no2) title1_binary,
     &                 num_rows2, dum1, dum2, dum3, num_values2
         read(unit=file_no2) 
         read(unit=file_no2) 
         write(out,*) '> header records read, file 2...'
      else   
         read(unit=file_no2) title1_binary, num_values2 
         read(unit=file_no2) 
         read(unit=file_no2) 
         write(out,*) '> header records read, file 2...'
         num_rows2 = num_elements
      end if
c
      if ( num_rows1 .ne. num_rows2 ) then
         write(out,*) '>> num records on the two files do not agree'
         stop
      end if
      num_rows = num_rows1

      if( num_values1 .ne. num_values2 ) then 
         write(out,*)
     &       '>> num values per record on the two files do not agree'
         stop
      end if
      num_values = num_values1
       
      write(out,*) '> start processing file records...'

      header_written = .false.
      
      do row = 1, num_rows
        if( mod(row,50000) == 0 ) 
     &    write(out,*) "... comparing row: ", row
        if( compare_nodal_results ) then
           read(file_no1) irow, results1(1:num_values,row)
           read(file_no2) jrow, results2(1:num_values,row)
        else   
           read(file_no1) irow, kk, results1(1:num_values,row)
           read(file_no2) jrow, kk, results2(1:num_values,row)
        end if           
        if( irow .ne. jrow ) then
          write(out,9300) row
          stop
        end if  
        do j = 1, num_values
         value1 = results1(j,row)
         value2 = results2(j,row)
         if( abs(value1) .le. 1e-10 ) value1 = 0.0  ! note 
         if( abs(value2) .le. 1e-10 ) value2 = 0.0  ! note
         write(charval1,9100) value1
         write(charval2,9100) value2
         if( charval1(1:num_chars_to_compare) .ne.
     &       charval2(1:num_chars_to_compare) ) then
           if( .not. header_written ) write(out,9200)
           header_written = .true.
           write(*,9000) row, j, value1, value2
         end if
        end do
      end do
c
c
      write(*,*) '>> all done'
      stop
 9000 format(20x,i8,i4,2x,e13.6,2x,e13.6)
 9100 format(e14.6)  !     compare to 3 figures 
 9200 format('... mismatch @ node/element   component, 2 values')
 9300 format('... fatal error type 1 @ row: ',i8)
      end
