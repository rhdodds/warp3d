c  *********************************************************************        
c  *                                                                   *        
c  *       subroutine name_strip                                       *        
c  *                                                                   *        
c  *            this routine strips the trailing blanks from the       *        
c  *            input string                                           *        
c  *                                                                   *        
c  *********************************************************************        
      subroutine name_strip ( instring, length)                                 
c
      implicit none
c 
      integer :: length
      character(LEN=length) :: instring                                            
c 
      integer :: last
      logical :: keep_strip                                                        
c                                                                               
c        input parameters                                                       
c               instring -- the input string                                    
c               length   -- passed in:   initial length                         
c                           passed out:  stripped length                        
c                                                                               
      last = length                                                             
      keep_strip = .true.                                                       
      do while (keep_strip)                                                     
         if (instring(last:last).eq.' ') then                                   
            last = last -1                                                      
            if (last.eq.0 ) keep_strip = .false.                                
         else                                                                   
            keep_strip  = .false.                                               
         endif                                                                  
      end do                                                                     
      length = last 
c                                                                  
      return                                                                    
      end                                                                       
