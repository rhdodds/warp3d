c     ****************************************************************
c     *                                                              *
c     *                      subroutine zero_vol                     *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 02/15/92                   *
c     *                                 : 02/08/94 rhd               *
c     *                                                              *
c     *     initialize volumetric terms for bbar                     *
c     *                                                              *
c     ****************************************************************
c
      subroutine zero_vol ( vol, volume, span, mxvl )     
      implicit integer ( a-z )
#dbl      double precision
#sgl      real
     &  vol(mxvl,8,*), volume(*), zero
      data zero / 0.0 /
      
c
      do i = 1, span 
c           
c           zero out 1st term for each node
c
         vol ( i, 1, 1 ) = zero
         vol ( i, 2, 1 ) = zero
         vol ( i, 3, 1 ) = zero
         vol ( i, 4, 1 ) = zero
         vol ( i, 5, 1 ) = zero
         vol ( i, 6, 1 ) = zero
         vol ( i, 7, 1 ) = zero
         vol ( i, 8, 1 ) = zero
c           
c           zero out 2nd term for each node
c   
         vol ( i, 1, 2 ) = zero
         vol ( i, 2, 2 ) = zero
         vol ( i, 3, 2 ) = zero
         vol ( i, 4, 2 ) = zero
         vol ( i, 5, 2 ) = zero
         vol ( i, 6, 2 ) = zero
         vol ( i, 7, 2 ) = zero
         vol ( i, 8, 2 ) = zero 
c           
c           zero out 3rd term for each node
c 
         vol ( i, 1, 3 ) = zero
         vol ( i, 2, 3 ) = zero
         vol ( i, 3, 3 ) = zero
         vol ( i, 4, 3 ) = zero
         vol ( i, 5, 3 ) = zero
         vol ( i, 6, 3 ) = zero
         vol ( i, 7, 3 ) = zero
         vol ( i, 8, 3 ) = zero 
c           
c           zero out total volume term for element
c 
         volume ( i ) = zero
c
      end do     
      return
      end
     


        





