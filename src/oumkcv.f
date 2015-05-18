c     ****************************************************************
c     *                                                              *
c     *                      subroutine oumkcv                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/10/04                    *
c     *                                                              *
c     *     this subroutine avarages element strains or stresses     *
c     *     at the gauss points to define a single set over the      *
c     *     element. this set is copied over all element gps.        *
c     *                                                              *
c     ****************************************************************
c
c           
      subroutine oumkcv( span, ngp, stress, num_short_stress,
     &                   num_short_strain )
      use elblk_data, only : elestr
      implicit integer (a-z)
$add param_def
      logical stress
c
c                       local declarations
c
#dbl      double precision
#sgl      real
     &     temvals(mxvl,mxoupr), zero, rngp
      data zero /0.0/
c
c                       set the actual number of primary output
c                       values. can be different for stresses
c                       and strains
c      
      num_vals = num_short_strain
      if ( stress )  num_vals = num_short_stress
c      
c                       1. zero accumulation array for average.
c          
      do j = 1, num_vals
        do i = 1, span
          temvals(i,j) = zero
        end do
      end do
c
c                       2. build summed values over gauss points
c                          for each stress/strain value. do all
c                          elements in block at same time
c          
      do k = 1, ngp
        do j = 1, num_vals
          do i = 1, span
             temvals(i,j) = temvals(i,j) + elestr(i,j,k)
          end do   
        end do
      end do   
c
c                       3. compute average value of each strain/stress
c                          value for each element.
c 
#sgl      rngp = real(ngp)         
#dbl      rngp = dble(ngp)         
      do j = 1, num_vals
        do i = 1, span
           temvals(i,j) = temvals(i,j) / rngp
        end do 
      end do
c  
c                       4. put average values within element at every
c                          gauss point of element
c          
      do k = 1, ngp
        do j = 1, num_vals
          do i = 1, span
             elestr(i,j,k) = temvals(i,j)
          end do   
        end do
      end do   
c
      return
      end
