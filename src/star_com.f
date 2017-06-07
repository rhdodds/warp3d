                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine star_com                     *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 06/27/2014 rhd             *          
c     *                                                              *          
c     *     interprets the special star commands:                    *          
c     *     commands that are preceeded by an astrick and are        *          
c     *     trapped directly by scan.                                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine star_com                                                       
      use global_data ! old common.main
      implicit integer (a-z)                                                    
      double precision                                                          
     &  dumd                                                                    
      real t1, wcputime, dumr                                                   
      external wcputime                                                         
      character(len=1) :: dums                                                  
      logical promsw,echosw,comsw,atrdsw,eolsw,eofsw,menusw,ptsw,               
     &     signsw                                                               
      logical matchs, debug                                                     
      data debug /.false./                                                      
      if (debug) write (out,*) '>>>>>>> inside star_com'                        
c                                                                               
c               command: input -- getting input from file                       
c                                                                               
      if(matchs('input',5)) then                                                
         call infile                                                            
c                                                                               
c               command: output -- writing output to a file                     
c                                                                               
      else if (matchs('output',6)) then                                         
         call outfil                                                            
c                                                                               
c               command: time -- writes out total wall time                     
c                                                                               
      else if(matchs('time',4)) then                                            
         t1 = wcputime(1)                                                       
         call errmsg (182,dum,dums,t1,dumd)                                     
c                                                                               
c               command: reset -- resets after a fatal error                    
c                                                                               
      else if(matchs('reset',5)) then                                           
         input_ok = .true.                                                      
         call errmsg(184,dum,dums,dumr,dumd)                                    
c                                                                               
c               command: echo -- sets the scan echo on or off                   
c                                                                               
      else if(matchs('echo',4)) then                                            
         nblank= 20                                                             
         reclen= 80                                                             
         endchr= 1h$                                                            
         promsw= .false.                                                        
         comsw= .false.                                                         
         atrdsw= .false.                                                        
         eolsw= .true.                                                          
         eofsw= .true.                                                          
         menusw= .false.                                                        
         ptsw= .false.                                                          
         signsw= .false.                                                        
         call scinit(nblank,reclen,endchr,promsw,echosw,comsw,atrdsw,           
     &        eolsw,eofsw,menusw,ptsw,signsw)                                   
c                                                                               
         if (matchs('off',3)) then                                              
            echosw= .false.                                                     
            call scinit(nblank,reclen,endchr,promsw,echosw,comsw,               
     &           atrdsw,eolsw,eofsw,menusw,ptsw,signsw)                         
         else                                                                   
            echosw= .true.                                                      
            call scinit(nblank,reclen,endchr,promsw,echosw,comsw,               
     &           atrdsw,eolsw,eofsw,menusw,ptsw,signsw)                         
         endif                                                                  
c                                                                               
      else                                                                      
         call errmsg (206,dum,dums,dumr,dumd)                                   
      endif                                                                     
c                                                                               
 9999 continue                                                                  
      if (debug) write (out,*) '<<<<<< leaving star_com'                        
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
