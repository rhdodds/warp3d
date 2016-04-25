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
     &  vol(mxvl,8,3), volume(*), zero
      data zero / 0.0d0 /
@!DIR$ ASSUME_ALIGNED vol:64, volume:64,      
c
      vol = zero
@!DIR$ LOOP COUNT MAX=###
      volume(1:span) = zero  
      return
      end
     


        





