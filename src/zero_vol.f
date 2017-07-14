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
      double precision
     &  vol(mxvl,8,3), volume(span), zero
         data zero / 0.0d0 /  
      vol = zero
      volume = zero
      return
      end









