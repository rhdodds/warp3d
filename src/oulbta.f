c     ****************************************************************
c     *                                                              *
c     *                      subroutine oulbta                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 08/16/89                   *
c     *                                                              *
c     *     this subroutine outputs the value of beta for an iter-   *
c     *     ation of the linear pcg algorithm                        *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine oulbta(step,stpitr,iter,beta)
      implicit integer (a-z)
      double precision
     &     beta
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
      write(out,9010) beta
c
c
 9000 format(///20x,'beta parameter for linear pcg algorithm'//8x,
     &       'step number',6x,': ',i7/8x,'step iteration number : ',
     &       i7/8x,'iteration number : ',i7//)
c
 9010 format(//11x,'beta parameter = ',e12.5/)
c
c
      return
      end
