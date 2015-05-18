c     ****************************************************************
c     *                                                              *
c     *                      subroutine oupele                       *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 3/11/04 (rhd)              *
c     *                                                              *
c     *     this subroutine copies the strain/stress results at the  *
c     *     element center into the data array for all element       *
c     *     results.                                                 *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine oupele( span, num_strain, num_stress, stress,
     &                   elem_results, nrowd )
      use elblk_data, only : elestr
      implicit integer (a-z)
      logical stress
$add param_def
#dbl      double precision
#sgl      real 
     &  elem_results(nrowd,*)
c
c                       copy element results into global
c                       data structure for element center results.
c
      num_vals = num_strain
      if ( stress ) num_vals = num_stress
c
      do k = 1, num_vals
            do i = 1, span
               elem_results(i,k) =  elestr(i,k,1)
            end do
      end do
c
      return
      end







