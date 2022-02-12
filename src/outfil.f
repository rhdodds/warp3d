c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine outfil                       *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 02/12/2022 rhd             *          
c     *                                                              *          
c     *     this subroutine gets and sets the name of a file to      *          
c     *     dump all the output to.                                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine outfil                                                         
      use global_data ! old common.main
      use file_info                                                             
      implicit none                                                             
c                                                                               
      integer :: dum, filchr, ierror, iunit, istat                                            
      real :: dumr                                                              
      double precision :: dumd                                                  
      character(len=8) :: dums                                                  
      character(len=80) :: filnam, outfle                                       
      logical :: file_exists, nameok, ok                                                    
      logical, external :: endcrd, label, matchs, string
      logical :: isatty 
c                                                                             
c              if "to display", then set screen as output device                        
c              Only allowed if stdout was *not* set to a file
c                                                                                
      if( matchs('to',2) )   call splunj       
      ok = .false.                                  
      if( matchs('display',5) ) ok = .true.
      if( matchs('window',5) ) ok = .true.
      if( ok ) then                                             
         if( outing ) close( out )   ! if now output to a file                                             
         call setout( outlun(1) )                                                 
         out = outlun(1)                                                        
         outing = .false.            
         if( .not. isatty( out ) ) then
           write(out,*) " "
           write(out,9020)
         end if                                           
         return
      end if                                                                     
c                                                                               
c              output to a file:                                       
c                check to see if filename is given, then get name as
c                   a label or a string
c                if ~  precedes name, get it converted to the full
c                   path name                      
c                                                                               
      if( matchs('file',4) ) call splunj                                         
      if( endcrd(dum) ) then                                                     
         call errmsg( 175, dum, dums, dumr, dumd )                                    
         return
      end if                                                                     
      if( outing ) then ! already writing to file (not stdout)                                                         
         call errmsg( 180, dum, dums, dumr, dumd )                                    
         return
      end if                                                                     
      filnam = ' '                                                              
      outfle = ' '                                                              
      if( label(dum) ) then                                                      
         call entits( filnam, filchr )                                            
         outfle = filnam(1:filchr)                                              
      else if( string(dum) ) then                                                
         call entits( filnam, filchr )                                            
         call tilde( filnam, outfle, nameok )                                       
       if( .not. nameok ) then                                                    
         call errmsg( 189, dum, outfle, dumr, dumd )                                  
         return
       end if
      else ! no file name given                                                                     
         call errmsg( 175, dum, dums, dumr, dumd )                                    
         return
      end if                                                                     
c                                                                               
c              check if file exists, then increase file number 
c              and open file.       
c              if file exists, delete.                                 
c                                                                               
      inquire( file = outfle, iostat = ierror, exist = file_exists )                     
      if( ierror > 0 ) then                                                     
         call errmsg( 176, dum, outfle, dumr, dumd )                                  
         return
      endif   
      if( file_exists ) then
         open( newunit=iunit, file=outfle, status='OLD', 
     &         iostat=istat )
         close( iunit, status='DELETE', iostat=istat )  
         write(out,9010) outfle
      end if                                                                
      open( unit=outlun(2), file = outfle, status = 'unknown' )                   
      write(out,9000) outfle                                                    
      out = outlun(2)                                                           
      call setout(outlun(2))                                                    
      outing = .true.                                                           
      return
c                                                                    
 9000 format('>>>>> output file is:',a)    
 9010 format('>>>>> existing file for output deleted: 'a) 
 9020 format('>>>>> *output command ignored. stdout set to a file',
     & /,    '       on startup. stdout must be the interactive',
     & /,    '       window for this command to work....',//)                                      
      end                                                                       
