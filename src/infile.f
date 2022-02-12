c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine infile                       *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified :  02/12/2022 rhd            *          
c     *                                                              *          
c     *     gets the name of a file for input as part of the         *
c     *     *input from <...>  command                               *          
c     *     opens and attaches file or window to the scanner         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                                                                                                             
      subroutine infile  
c                                                       
      use global_data ! old common.main
      use file_info
c                                                             
      implicit none  
c
      integer :: dum, filchr, ierror                                                  
      real :: dumr                                                              
      character(len=8) :: dums                                                  
      character(len=80) :: filnam, infil                                        
      logical :: ok, nameok, isatty
      logical, external :: endcrd, label, matchs, string                         
      double precision :: dumd                                                  
c                                                                               
c              if "to window", then set screen as output device.
c              allowed only if stdin *not* set to a file
c                                                                               
      if( matchs('from',4) ) call splunj   
      ok = .false.                                     
      if( matchs('terminal',5) ) ok = .true.
      if( matchs('window',5) )   ok = .true.
      if( matchs('screen',5) )   ok = .true.
      if( matchs('display',5) )  ok = .true.
      if( ok ) then                                     
         filcnt = 1                                                             
         call setin( inlun(filcnt) )                                             
         in = inlun(filcnt)   
         if( .not. isatty( in ) ) then
           write(out,*) " "
           write(out,9020)
         end if                                           
         return                                                                 
      end if                                                                     
c                                                                               
c              from a file. file name must be given. resolve full path
c              name if ~ is present
c                                                                               
      if( matchs('file',4) ) call splunj                                        
      if( endcrd(dum) ) then                                                     
         call errmsg( 175, dum, dums, dumr, dumd )                                    
         return                                                                 
      end if                                                                     
      filnam = ' '                                                              
      infil = ' '                                                               
      if( label(dum) ) then                                                     
         call entits( filnam, filchr )                                            
         infil = filnam(1:filchr)                                               
      else if( string(dum) ) then                                               
         call entits( filnam, filchr )                                            
         call tilde( filnam, infil, nameok )                                       
         if( .not. nameok ) then                                                 
            call errmsg( 189, dum, infil, dumr, dumd )                                
            return                                                              
         end if                                                                  
      else                                                                      
         call errmsg( 175, dum, dums, dumr, dumd)                                    
         return                                                                 
      end if                                                                     
c                                                                               
c              file must exist. add to input file stack. 
c              open file, attach to scan as current source
c                                                                                             
      inquire( file = infil, iostat = ierror, exist = ok )                      
      if( ierror > 0 ) then                                                  
         call errmsg( 176, dum, infil, dumr, dumd )                                   
         return                                                                 
      else if( .not. ok ) then                                                  
         call errmsg( 177, dum, infil, dumr, dumd )                                   
         return                                                                 
      end if                                                                     
      filcnt = filcnt + 1                                                       
      if( filcnt > max_opened_input_files ) then                             
         call errmsg( 178, dum, dums, dumr, dumd )                                    
         filcnt = filcnt - 1                                                      
         return                                                                 
      end if                                                                     
      open( unit=inlun(filcnt), file = infil, status = 'old' )                  
      in = inlun(filcnt)                                                        
      call setin( inlun(filcnt) )                                                 
      write(out,9000) trim( infil )                                             
c                                                                               
      return                                                                    
 9000 format (1x,'>>>>> input file is: ',a)                                     
 9020 format('>>>>> *input command ignored. stdin set to a file',
     & /,    '       on startup. stdin must be the interactive',
     & /,    '       window for this command to work....',//)                                      
      end                                                                       
                                                                                
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *               subroutine infile_stpdrv_open                  *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 02/12/2022 rhd             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine infile_stpdrv_open( file_name ) 
c                               
      use global_data, only : in, out 
      use file_info       
c                                                      
      implicit none
c                                                    
      character(len=80) :: file_name                                            
      logical :: ok                                                                
      integer :: ierror                                                            
c                                                                               
      inquire ( file = file_name, iostat = ierror, exist = ok )                 
c                                                                               
      if( .not. ok ) then                                                       
         write(out,9110) file_name                                              
         call die_gracefully                                                    
      end if                                                                    
c                                                                               
      filcnt = filcnt + 1                                                       
      if( filcnt .gt. max_opened_input_files ) then                             
         write(out,9120) max_opened_input_files                                 
         call die_gracefully                                                    
      end if                                                                    
c                                                                               
      open( unit = inlun(filcnt), file = file_name, status = 'old')             
      in = inlun(filcnt)                                                        
      call setin(inlun(filcnt))                                                 
c                                                                               
      return                                                                    
c                                                                               
 9110 format(                                                                   
     & /1x, '>>>>> error: the previously specified output commands',            
     & /14x,'file no longer exists. file: ',a80,                                
     & /14x,'job terminated...',//)                                             
 9120 format(                                                                   
     & /1x, '>>>>> error: trying to open the specified output commands',        
     & /14x,'file. too many files are open = ',i5,                              
     & /14x,'job terminated...',//)                                             
c                                                                               
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *               subroutine infile_stpdrv_close                 *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 02/12/2022 rhd             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine infile_stpdrv_close( file_name )     
c                          
      use global_data, only : out,in
      use file_info    
c                                                         
      implicit none
c
      character(len=80) :: file_name                                            
      logical :: now_open                                                      
      integer :: ierror, dum
      real :: dumr
      double precision :: dumd
      character(len=1) :: dums                                                            
                                                                                
c                                                                               
      inquire ( file = file_name, iostat = ierror, opened = now_open )                                            
c                                                                               
      if( .not. now_open ) then                                                 
         write(out,9110) file_name                                              
         call die_gracefully                                                    
      end if                                                                    
c                                                                               
      close( inlun(filcnt) )                                                    
      filcnt = filcnt-1                                                         
      if( filcnt .eq. 0 ) then                                                  
         call errmsg( 204, dum, dums, dumr, dumd )                                    
         call die_gracefully                                                    
      end if                                                                    
c                                                                               
      call setin( inlun(filcnt) )                                                 
      in = inlun(filcnt)                                                        
c                                                                               
      return                                                                    
c                                                                               
 9110 format(                                                                   
     & /1x, '>>>>> error: the previously specified output commands',            
     & /14x,'file is not open. an invalid close is being attempted in'          
     & /14x,'routine infile_stpdrv_close'                                       
     & /14x,'output commands file: ',a80,                                       
     & /14x,'job terminated...',//)                                             
c                                                                               
      end                                                                       
                                                                                
                                                                                
                                                                                
