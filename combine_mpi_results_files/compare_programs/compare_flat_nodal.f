      program compare_flat_nodal 
c      
      implicit none
c      
      integer :: out, num_chars_to_compare, file_no1, file_no2,
     &           header_written, stat1, stat2, num_nodes,
     &           j, num_values, num_vals_to_compare,
     &           node, last_slash, num_flat_node_strains,
     &           num_flat_node_stresses
      logical :: bad1, bad2, stream, stream1, stream2, 
     &           do_stresses, do_strains, text, text1, text2 
      double precision, allocatable :: results1(:,:), results2(:,:)
      double precision :: value1, value2, zero_tol
      character(len=80) :: file1_name, file2_name
      character(len=50) :: root1, root2
      character(len=20) charval1, charval2, form, access
c
      out = 6
      num_chars_to_compare = 5
      num_vals_to_compare  = 6
      file_no1 = 20
      file_no2 = 21
      num_flat_node_strains  = 22
      num_flat_node_stresses = 26
c      
      write(out,9005)
      write(out,*) "... Compare Flat Nodal Results"
      write(out,*) "    =========================="
      write(out,*) " "
      write(out,*) '>> flat nodal results file 1: '
      read(*,*) file1_name
      write(out,*) '>> flat nodal results file 2: '
      read(*,*) file2_name 
      write(out,*) 
     & '>> # digits to compare, # values to compare each node:'
      read(*,*) num_nodes, num_chars_to_compare, num_vals_to_compare
      write(out,*) 
     & '>> set values smaller than this to zero:'
      read(*,*) zero_tol
c
      bad1 = num_nodes <= 0 .or.
     &       num_chars_to_compare <= 0  .or.
     &       num_vals_to_compare <= 0   .or.
     &       zero_tol <= 0.0  
      if( bad1 ) then
        write(out,9120) num_nodes, num_chars_to_compare,
     &                  num_vals_to_compare, zero_tol
        stop
      end if  
c
c              get file name w/o directories
c
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
c
      bad1 = (root1(2:2) .ne. 'n')  .or. (root2(2:2) .ne. 'n')
      if( bad1 ) then
        write(out,9130)
        stop
      end if
c      
      stream1 = index( root1, "stream" ) > 0
      text1   = index( root1, "text" )   > 0
      stream2 = index( root2, "stream" ) > 0
      text2   = index( root2, "text" )   > 0
      stream = stream1 .and. stream2
      text   = text1 .and. text2
      if( .not. (stream .or. text ) ) then
        write(out,9140)
        stop
      end if
c
c      write(out,*) '.. root1: ', root1
c      write(out,*) '.. root2: ', root2
c      write(out,*) '.. file 1: ', file1_name
c      write(out,*) '.. file 2: ', file2_name

      num_values = 0
c      
      do_stresses = root1(3:3) .eq. 's' .and.
     &              root2(3:3) .eq. 's'
      do_strains  = root1(3:3) .eq. 'e' .and.
     &              root2(3:3) .eq. 'e'
       
      if( do_stresses ) num_values = num_flat_node_stresses
      if( do_strains )  num_values = num_flat_node_strains

c
      form = "formatted"
      if( stream ) form = "unformatted"
      access = 'sequential'
      if( stream ) access = "stream"
      open( unit=file_no1, file=file1_name, status='old', iostat=stat1,
     &      form = form, access = access   )
      open( unit=file_no2, file=file2_name, status='old', iostat=stat2,
     &      form = form, access = access   )
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
      if( text ) then
         call skip_comment_lines( file_no1 )
         call skip_comment_lines( file_no2 )
         write(out,*) '> header records read...'
      end if
c

      bad2 = num_values .eq. 0   
      if( bad2 ) then
        write(out,9100) num_values
        stop
      end if

c
      allocate( results1(num_values,num_nodes), 
     &          results2(num_values,num_nodes) )
c
      write(out,*) '>> start elements processing...'
      header_written = .false.
c
      do node = 1, num_nodes
       if( mod(node,1000) == 0 )
     &         write(out,*) "... comparing node: ", node
        if( text ) then
          read(unit=file_no1,fmt=9330) results1(1:num_values,node)
          read(unit=file_no2,fmt=9330) results2(1:num_values,node)
        end if
        if( stream ) then
          read(file_no1)  results1(1:num_values,node)
          read(file_no2)  results2(1:num_values,node)
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
 9005 format(///)
 9010 format(20x,i8,i4,2x,e13.6,2x,e13.6)
 9100 format('>> FATAL: inconsistent file 1 and 2 data:',
     & /,10x,'...num_values: ',i10)      
 9110 format('>> FATAL: inconsistent file 1 and 2 data:',
     & /,10x,'...node1, node2: ',2i10 )     
 9120 format('>> FATAL: inconsistent input:',
     & /,10x,'...num_nodes: ',i10,
     & /,10x,'...num_chars_to_compare: ',i10,
     & /,10x,'...num_vals_to_compare:  ',i10,
     & /,10x,'...zero_tol:             ',e14.6 )
 9130 format('>> FATAL: not element result files')
 9140 format('>> FATAL: inconsistent input:',
     & /,10x,'...both files must be text or stream')
 9160 format(' >> comparing 2 files: ',
     & /10x,a,/,10x,a)
 9200 format('... mismatch @   node       component, 2 values')
 9210 format(e14.6) 
 9300 format(80a1)
 9315 format(2i9,e15.6,2i9)
 9330 format(10000e15.6)

c 
      end
c **********************************************************************
c *                                                                    *
c *          skip_comment_lines                                        *
c *                                                                    *
c **********************************************************************
c
c
      subroutine skip_comment_lines( file_no )
      implicit none
c
      integer :: file_no
c
      character :: first_char*1
c
      do 
        read(file_no,fmt="(a1)") first_char
        if( first_char .eq. "#"  .or.
     &      first_char .eq. "!"  .or.
     &      first_char .eq. "c"  .or.
     &      first_char .eq. "C"  ) cycle
        backspace(file_no)
        return
       end do
c
      return
      end subroutine skip_comment_lines