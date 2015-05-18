c     ****************************************************************
c     *                                                              *
c     *                      subroutine aplacc                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 06/23/91                   *
c     *                                                              *
c     *     this subroutine applies the crisfield accelerating       *
c     *     algorithm to the current incremental displacement        *
c     *     vector.                                                  *
c     *           (this routine currently inactive)                  *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine aplacc
c      implicit integer (a-z)
c      logical trcacc
c$add common.main
c#dbl      double precision
c#sgl      real
c     &   zero, one
c      data zero, one / 0.0, 1.0 /
cc
c      trcacc = trace(3)
cc
cc                       compute the accelerating parameters.
cc      
c      if( iter .gt. 1 ) then
cc
c         aparm= zero
c         bparm= zero
c         cparm= zero
c         do i= 1, nodof                  
c            aparm= aparm + oldidu(i)*oldres(i)                  
c            bparm= bparm + oldidu(i)*(res(i)-oldres(i))
c            cparm= cparm + idu(i)*(res(i)-oldres(i))
c         end do
c         fparm= -(aparm/bparm)
c         eparm= fparm*(one-(cparm/bparm))-one
cc
cc                       check to make sure eparm and fparm are within
cc                       user specified bounds.        
cc                    
c         if( abs(eparm) .gt. emax ) eparm= emax
c         if( fparm .gt. fmax )      fparm= fmax
cc
cc                       apply the accelerators to the current    
cc                       incremental displacement vector, idu.
cc
c         do i = 1, nodof
c            idu(i)= eparm*oldidu(i) + fparm*idu(i)   
c         end do
cc                                                             
c      else                              
c         aparm= zero
c         bparm= zero   
c         cparm= zero   
c         eparm= zero   
c         fparm= one
c      end if
cc
cc                       set the previous incremental displacement
cc                       vector and the previous residual load vector
cc                       for use in the next iteration.
cc              
c      do i = 1, nodof
c         oldidu(i)= idu(i)
c         oldres(i)= res(i)
c      end do
cc
cc                       print out the values of the parameters if
cc                       the trace acceleration flag is on.
cc              
c      if(trcacc) call ouparm(aparm,bparm,cparm,eparm,fparm,step,iter)
cc
cc
c      return
      end
