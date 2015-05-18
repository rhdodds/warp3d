c     ****************************************************************
c     *                                                              *
c     *                      subroutine oualph                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 03/15/94                   *
c     *                                                              *
c     *     this subroutine outputs the value of the step length     *
c     *     for an iteration of the linear pcg algorithm.            * 
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine oualph(step,stpitr,iter,alpha)
      implicit integer (a-z)
#dbl      double precision
#sgl      real
     &     alpha        
      logical dummy
c
c                       set device numbers.
c
      call iodevn(in,out,dummy,1)
c
c                       write out header
c
      write(out,9000) step,stpitr,iter
c
c                       write out parameters.
c                       
      write(out,9010) alpha
c
c        
 9000 format(///20x,'step length for linear pcg algorithm'//8x,
     &       'step number',6x,': ',i6/8x,'step iteration number : ',
     &       i6/8x,'iteration number : ',i6//)
c
 9010 format(//11x,'step length = ',e12.5/)
c
c
      return
      end
