c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine outfil                       *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 01/20/2017 rhd             *          
c     *                                                              *          
c     *     this subroutine gets and sets the name of a file to      *          
c     *     dump all the output to.                                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine outfil                                                         
      use global_data ! old common.main
      use file_info                                                             
      implicit none                                                             
c                                                                               
      integer :: dum, filchr, ierror                                            
      real :: dumr                                                              
      double precision :: dumd                                                  
      character(len=8) :: dums                                                  
      character(len=80) :: filnam, outfle                                       
      logical :: ok, nameok                                                     
      logical, external :: endcrd, label, matchs, string                        
c                                                                               
c                       if "to terminal" or "to display"                        
c                       then set screen as output device                        
c                                                                               
      if (matchs('from',4)) call splunj                                         
      if (matchs('to',2))   call splunj                                         
      if (matchs('is',2))   call splunj                                         
      if (matchs('terminal',5)) then                                            
         if (outing) close (out)                                                
         call setout(outlun(1))                                                 
         out = outlun(1)                                                        
         outing = .false.                                                       
         go to 2110                                                             
      endif                                                                     
      if (matchs('display',5)) then                                             
         if (outing) close (out)                                                
         call setout(outlun(1))                                                 
         out = outlun(1)                                                        
         outing = .false.                                                       
         go to 2110                                                             
      endif                                                                     
c                                                                               
c                       output to a file:                                       
c                         check to see if filename is given, then               
c                         get name as a label or a string                       
c                                                                               
c                                                                               
      if (matchs('file',4)) call splunj                                         
      if (endcrd(dum)) then                                                     
         call errmsg(175,dum,dums,dumr,dumd)                                    
         goto 2110                                                              
      endif                                                                     
      if (outing) then                                                          
         call errmsg(180,dum,dums,dumr,dumd)                                    
         goto 2110                                                              
      endif                                                                     
      filnam = ' '                                                              
      outfle = ' '                                                              
      if (label(dum)) then                                                      
         call entits (filnam,filchr)                                            
         outfle = filnam(1:filchr)                                              
      else if (string(dum)) then                                                
         call entits (filnam,filchr)                                            
         call tilde(filnam,outfle,nameok)                                       
       if (.not.nameok) then                                                    
         call errmsg(189,dum,outfle,dumr,dumd)                                  
         goto 2110                                                              
         endif                                                                  
      else                                                                      
         call errmsg(175,dum,dums,dumr,dumd)                                    
       goto 2110                                                                
      endif                                                                     
c                                                                               
c                  check if file exists, then increase file                     
c                  number and open file.                                        
c                                                                               
      inquire ( file = outfle, iostat = ierror, exist = ok)                     
      if (ierror.gt.0) then                                                     
         call errmsg(176,dum,outfle,dumr,dumd)                                  
         goto 2110                                                              
      endif                                                                     
      open(unit=outlun(2), file = outfle, status = 'unknown')                   
      write(out,9000) outfle                                                    
      out = outlun(2)                                                           
      call setout(outlun(2))                                                    
      outing = .true.                                                           
c                                                                               
 2110 continue                                                                  
c                                                                               
      return                                                                    
 9000 format('>>>>> output file is:',a)                                         
      end                                                                       
