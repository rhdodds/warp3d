c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine infile                       *          
c     *                                                              *          
c     *                       written by : AG                        *          
c     *                                                              *          
c     *                   last modified :  01/20/2016 rhd            *          
c     *                                                              *          
c     *     this subroutine gets the name of a file for input        *          
c     *     as part of the *input from file command                  *          
c     *     opens and attaches file to the scanner                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine infile                                                         
      use global_data ! old common.main
      use file_info                                                             
      implicit integer (a-z)                                                    
      real :: dumr                                                              
      character(len=8) :: dums                                                  
      character(len=80) :: filnam, infil                                        
      logical :: ok, nameok, endcrd,label,matchs,string                         
      double precision :: dumd                                                  
c                                                                               
c                       if "from terminal" or "from display"                    
c                       then set keyboard as input  device                      
c                                                                               
      if( matchs('from',4) ) call splunj                                        
      if( matchs('terminal',5) ) then                                           
         filcnt = 1                                                             
         call setin (inlun(filcnt))                                             
         in = inlun(filcnt)                                                     
         return                                                                 
      endif                                                                     
      if( matchs('display',5) ) then                                            
         filcnt = 1                                                             
         call setin (inlun(filcnt))                                             
         in = inlun(filcnt)                                                     
         return                                                                 
      endif                                                                     
c                                                                               
c                       file is the input device:                               
c                         check to see if the file name is given,               
c                         then get name as label or string                      
c                                                                               
      if( matchs('file',4) ) call splunj                                        
      if (endcrd(dum)) then                                                     
         call errmsg(175,dum,dums,dumr,dumd)                                    
         return                                                                 
      endif                                                                     
      filnam = ' '                                                              
      infil = ' '                                                               
      if( label(dum) ) then                                                     
         call entits (filnam,filchr)                                            
         infil = filnam(1:filchr)                                               
      else if( string(dum) ) then                                               
         call entits (filnam,filchr)                                            
         call tilde (filnam,infil,nameok)                                       
         if( .not.nameok ) then                                                 
            call errmsg(189,dum,infil,dumr,dumd)                                
            return                                                              
         endif                                                                  
      else                                                                      
         call errmsg(175,dum,dums,dumr,dumd)                                    
         return                                                                 
      endif                                                                     
c                                                                               
c                  check if file is there, then increase                        
c                  file number and open                                         
c                                                                               
      inquire( file = infil, iostat = ierror, exist = ok )                      
      if( ierror .gt. 0 ) then                                                  
         call errmsg(176,dum,infil,dumr,dumd)                                   
         return                                                                 
      else if( .not. ok ) then                                                  
         call errmsg(177,dum,infil,dumr,dumd)                                   
         return                                                                 
      endif                                                                     
      filcnt = filcnt + 1                                                       
      if( filcnt .gt. max_opened_input_files ) then                             
         call errmsg(178,dum,dums,dumr,dumd)                                    
         filcnt = filcnt-1                                                      
         return                                                                 
      endif                                                                     
      open( unit=inlun(filcnt), file = infil, status = 'old' )                  
      in = inlun(filcnt)                                                        
      call setin(inlun(filcnt))                                                 
      write(out,9000) trim( infil )                                             
c                                                                               
      return                                                                    
 9000 format (1x,'>>>>> input file is: ',a)                                     
      end                                                                       
                                                                                
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *               subroutine infile_stpdrv_open                  *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 06/27/2014                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine infile_stpdrv_open( file_name )                                
      use global_data ! old common.main
      use file_info                                                             
      implicit integer (a-z)                                                    
      character(len=80) :: file_name                                            
      logical ok                                                                
      integer ierror                                                            
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
c     *                   last modified : 06/27/2014                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine infile_stpdrv_close( file_name )                               
      use global_data ! old common.main
      use file_info                                                             
      implicit integer (a-z)                                                    
      character(len=80) :: file_name                                            
      logical ok, now_open                                                      
      integer ierror                                                            
                                                                                
c                                                                               
      inquire ( file = file_name, iostat = ierror,                              
     &          opened = now_open  )                                            
c                                                                               
      if( .not. now_open ) then                                                 
         write(out,9110) file_name                                              
         call die_gracefully                                                    
      end if                                                                    
c                                                                               
      close( inlun(filcnt) )                                                    
      filcnt = filcnt-1                                                         
      if( filcnt .eq. 0 ) then                                                  
         call errmsg(204,dum,dums,dumr,dumd)                                    
         call die_gracefully                                                    
      end if                                                                    
c                                                                               
      call setin(inlun(filcnt))                                                 
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
                                                                                
                                                                                
                                                                                
