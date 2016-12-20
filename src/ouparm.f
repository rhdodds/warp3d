c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouparm                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 02/17/88                   *
c     *                                                              *
c     *     this subroutine outputs the parameters of the crisfield  *
c     *     accelerating algorithm.                                  * 
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine ouparm(aparm,bparm,cparm,eparm,fparm,step,iter)
      implicit integer (a-z)
      double precision
     &     aparm,bparm,cparm,eparm,fparm                                 
      logical dummy
c
c                       set device numbers.
c
      call iodevn(in,out,dummy,1)
c
c                       write out header
c
      write(out,9000) step,iter
c
c                       write out parameters.
c                       
      write(out,9010) aparm,bparm,cparm,eparm,fparm                                 
c
c        
 9000 format(///20x,'crisfield accelerating parameters'//8x,
     &       'step number',6x,': ',i6/8x,'iteration number : ',i6//)
c
 9010 format(//11x,'parameter a : ',e12.5/11x,
     &             'parameter b : ',e12.5/11x,
     &             'parameter c : ',e12.5/11x,
     &             'parameter e : ',e12.5/11x,
     &             'parameter f : ',e12.5 )
c
c
      return
      end
