c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine iodevn                       *          
c     *                                                              *          
c     *                       written by : bh                        *          
c     *                                                              *          
c     *                   last modified : 4/28/26 rhd                *          
c     *                                                              *          
c     *     provides a calling subprogram with the                   *          
c     *     input and output device number, and whether or not the   *          
c     *     trace solution flag is on.                               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine iodevn( innum, outnum, trc, trctyp)                                 
c
      use global_data, only : in, out, trace
c      
      implicit none
c 
      integer :: innum, outnum, trctyp      
      logical :: trc                                                               
c                                                                               
      innum = in                                                                 
      outnum = out                                                               
c                                                                               
      select case( trctyp )
c                                                                               
      case( 1 ) 
      trc = trace(1)                                                             
c                                                                               
      case( 2 ) 
      trc = trace(2)                                                             
c                                                                               
      case( 3 ) 
      trc = trace(3)                                                             
c                                                                               
      case( 4 ) 
      trc = trace(4)                                                             
c                                                                               
      case( 5 ) 
      trc = trace(5)                                                             
c  
      end select
c
      return      
      end                                                                       
