c**********************************************************************         
c                                                                     *         
c                  subroutine tilde                                   *         
c                                                                     *         
c            written by;  rhd                                         *         
c            updated:  rhd 12/19/2016 rhd                             *         
c                                                                     *         
c            given a file name, replace ~/ if present with user's     *         
c            home directory (full path name)                          *         
c                                                                     *         
c**********************************************************************         
c                                                                               
      subroutine tilde( name_in, name_out, ok )                                 
      use main_data, only : windows_os                                          
      implicit none                                                             
c                                                                               
c              parameter declarations                                           
c                                                                               
      character (len=*) :: name_in, name_out  !  in -> out                      
      logical :: ok                                                             
c                                                                               
c              local declarations                                               
c                                                                               
      character (len=1000) :: work_name, home_name                              
      character (len=1) :: dums, file_separator                                 
      integer :: len_out, len_in, last_work_name, last_home_name,               
     &           reqd_length, dumi, i                                           
      logical :: debug                                                          
      double precision :: dumd                                                  
      real dumr                                                                 
                                                                                
                                                                                
      ok    = .false.                                                           
      debug = .true.                                                            
      name_out = " "                                                            
                                                                                
      len_out = len( name_out )                                                 
      len_in  = len( name_in )                                                  
      work_name = " "                                                           
      work_name(1:) = adjustl( name_in )                                        
      last_work_name = index( work_name, " " ) - 1                              
c                                                                               
      if( work_name(1:1) .eq. "/" ) then                                        
        call errmsg( 187, dumi, dums, dumr, dumd )                              
        return                                                                  
      end if                                                                    
c                                                                               
      if( work_name(1:2) .ne. "~/" ) then                                       
         if( len_out .lt. last_work_name ) return                               
         name_out(1:) = work_name(1:last_work_name)                             
         ok = .true.                                                            
         return                                                                 
      end if                                                                    
c                                                                               
      home_name = " "                                                           
      call getenv( "HOME", home_name )                                          
      home_name = adjustl( home_name )                                          
      last_home_name = index( home_name, " " ) - 1                              
c                                                                               
c               must find last non-blank character the hard way                 
c               since the $HOME directory may have blanks in the                
c               name (ie Windows)                                               
c                                                                               
      do i = 1000, 1, -1  ! handles when name contains blanks                   
        last_home_name = i                                                      
        if( home_name(i:i) .ne. " " ) exit                                      
      end do                                                                    
      if( last_home_name .le. 1 ) then                                          
        call errmsg( 189, dumi, dums, dumr, dumd )                              
        return                                                                  
      end if                                                                    
c                                                                               
      reqd_length = last_home_name + last_work_name                             
      if( reqd_length .gt. len_out ) then                                       
        ok = .false.                                                            
        return                                                                  
      end if                                                                    
c                                                                               
      file_separator = "/"                                                      
      if( windows_os ) file_separator = "\"                                     
      name_out(1:) = home_name(1:last_home_name) // file_separator              
     &                      // work_name(3:last_work_name)                      
      ok = .true.                                                               
      return                                                                    
c                                                                               
      end                                                                       
