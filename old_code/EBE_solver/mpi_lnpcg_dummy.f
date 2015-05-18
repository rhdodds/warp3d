c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine wmpi_init_owner              *
c     *                                                              *
c     *                       written by : ag                        *
c     *                                                              *
c     *                   last modified : 02/04/98                   *
c     *                                                              *
c     *       Dummy routine for serial version of warp3d.            *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine wmpi_init_owner
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine wmpi_chknode                     *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 04/06/98                   *
c     *                                                              *
c     *           The processor assesses whether the node number     *
c     *           passed in as an argument is accessed by this       *
c     *           processor. In serial case, it is always true.      *
c     *                                                              *
c     ****************************************************************
c     
c     
      subroutine wmpi_chknode ( node, referenced, owned )
      implicit integer (a-z)
      logical referenced, owned
c     
      referenced = .true.
      owned = .true.         
c     
      return 
      end
c
c     ****************************************************************
c     *                                                              *
c     *                  subroutine lnpcg_comm                       *
c     *                                                              *
c     *                       written by : asg                       *
c     *                                                              *
c     *                   last modified : 11/18/98                   *
c     *                                                              *
c     *           Dummy routine for serial version.                  *
c     *                                                              *
c     ****************************************************************
c     
c     
      subroutine lnpcg_comm ( ddum, idum1, idum2, idum3)
      implicit integer (a-z)
#dbl      double precision 
#sgl      real
     &     ddum
c
      return
      end     


