      implicit integer (a-z)
      character * 20 search_string
      character * 80 line
c
      open(unit=10,file='Makefile.hpux.8000',status='old')
c
 1    continue
      write(*,*) 'string to find: '
      read(*,*)  search_string
      search_string = adjustl( search_string)
      num_search_chars = len_trim(search_string)
c
      rewind(unit=10)
      do 
       read(10,1000,end=2000) line(1:)
       line = adjustl( line )
       if ( index( line(1:),
     &      search_string(1:num_search_chars) ) .gt. 0 ) then
         write(*,*) ' '
         write(*,*) trim(line)
         read(10,1000) line
         write(*,*) trim(line)
         read(10,1000) line
         write(*,*) trim(line)
         write(*,*) ' '
         cycle
       end if
      end do
c
 2000 continue
      write(*,*) '>> end of file...'
      go to 1
c
1000  format(a80)
      end

     
           
