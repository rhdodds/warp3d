      implicit none
c      
      integer, parameter :: max = 3000000
      integer, parameter :: num_flat_nodal_strains = 22
      integer, parameter :: num_flat_nodal_stresses = 26


      integer :: num_values, file_no1, file_no2, out,
     &           num_digits_to_compare, num_nodes, row, j,
     &           num_components_to_compare
      character(len=100) :: file1_name, file2_name
      double precision :: results1(30,max), results2(30,max),
     &                    value1, value2, small_value
      character(len=14) charval1, charval2
      logical :: header_written, do_stresses, do_strains
c
      out                   = 6
      num_digits_to_compare = 6
      num_components_to_compare = 20
      num_nodes             = 1568 ! required to compare element results
      small_value           = 1.0d-4
c      
      do_stresses              = .true.
      do_strains               = .not. do_stresses
c      
      if( do_stresses ) num_values = num_flat_nodal_stresses
      if( do_strains )  num_values = num_flat_nodal_strains
c      
      write(*,*) '>> binary file 1: '
c      read(*,*) file1_name
      file1_name = "./wns00030_stream"
      write(*,*) '>> binary file 2: '
c      read(*,*) file2_name
      file2_name = "./mpi_wns00030_stream"
c
      file_no1 = 20
      file_no2 = 21
c
      open( unit=file_no1, file=file1_name, status='old', 
     &       access='stream', form='unformatted' )
      open( unit=file_no2, file=file2_name, status='old', 
     &       access='stream', form='unformatted' )
c
      write(out,*) '>> files opened'
      write(out,*) '> start processing file records...'
c
      header_written = .false.
c      
      do row = 1, num_nodes
        if( mod(row,50000) == 0 ) 
     &    write(out,*) "... comparing row: ", row
        read(file_no1) results1(1:num_values,row)
        read(file_no2) results2(1:num_values,row)
        do j = 1, num_components_to_compare
         value1 = results1(j,row)
         value2 = results2(j,row)
         if( abs(value1) .le. small_value ) value1 = 0.0d0  ! note 
         if( abs(value2) .le. small_value ) value2 = 0.0d0  ! note
         write(charval1,9100) value1
         write(charval2,9100) value2
         if( charval1(1:num_digits_to_compare) .ne.
     &       charval2(1:num_digits_to_compare) ) then
           if( .not. header_written ) write(out,9200)
           header_written = .true.
           write(*,9000) row, j, value1, value2
         end if
        end do
      end do
c
      write(out,*) '>> all done'
      stop
c      
 9000 format(20x,i8,i4,2x,e13.6,2x,e13.6)
 9100 format(e14.6) 
 9200 format('... mismatch @ node/element   component, 2 values')
 9300 format('... fatal error type 1 @ row: ',i8)
      end
