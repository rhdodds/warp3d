c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine iodevn                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 06/23/91                   *          
c     *                                                              *          
c     *     this subroutine provides a calling subprogram with the   *          
c     *     input and output device number, and whether or not the   *          
c     *     trace solution flag is on.                               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine iodevn(innum,outnum,trc,trctyp)                                
      use global_data ! old common.main
      implicit integer (a-z)                                                    
      logical trc                                                               
c                                                                               
      innum= in                                                                 
      outnum= out                                                               
c                                                                               
      go to (100,200,300,400,500) trctyp                                        
c                                                                               
 100  continue                                                                  
      trc= trace(1)                                                             
      go to 9999                                                                
c                                                                               
 200  continue                                                                  
      trc= trace(2)                                                             
      go to 9999                                                                
c                                                                               
 300  continue                                                                  
      trc= trace(3)                                                             
      go to 9999                                                                
c                                                                               
 400  continue                                                                  
      trc= trace(4)                                                             
      go to 9999                                                                
c                                                                               
 500  continue                                                                  
      trc= trace(5)                                                             
      go to 9999                                                                
c                                                                               
 9999 return                                                                    
      end                                                                       
