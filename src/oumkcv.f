c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oumkcv                       *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 08/4/2025 rhd              *          
c     *                                                              *          
c     *     averages element strains or stresses                     *          
c     *     at the gauss points to define a single set over the      *          
c     *     element. this set is copied over all element gps.        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oumkcv( span, ngp, do_stresses, num_short_stress,              
     &                   num_short_strain )                                     
      use elblk_data, only : elestr  
      use constants                                           
      implicit none                                                             
c                                                                               
      include 'param_def'                                                       
      integer :: span, ngp, num_short_stress, num_short_strain                  
      logical :: do_stresses                                                    
c                                                                               
c                       local declarations                                      
c                                                                               
      integer :: i, j, k, num_vals                                              
      double precision :: temvals(mxvl,mxoupr), rngp                            
c                                                                               
c                       set the actual number of primary output                 
c                       values. can be different for stresses                   
c                       and strains                                             
c                                                                               
      num_vals = num_short_strain                                               
      if( do_stresses ) num_vals = num_short_stress                            
c                                                                               
c                       1. zero accumulation array for average.                 
c                                                                               
      temvals = zero                                                            
c                                                                               
c                       2. build summed values over gauss points                
c                          for each stress/strain value. do all                 
c                          elements in block at same time                       
c      
      if( ngp .eq. 8 ) then                                                     
c                                                                               
       do k = 1, 8                                                              
        do j = 1, num_vals                                                      
!$omp simd
          do i = 1, span                                                        
             temvals(i,j) = temvals(i,j) + elestr(i,j,k)                        
          end do                                                                
        end do                                                                  
       end do                                                                   
c                                                                               
      else  ! not 8 int points                                                  
c                                                                               
       do k = 1, ngp                                                            
        do j = 1, num_vals                                                      
!$omp simd
          do i = 1, span                                                        
             temvals(i,j) = temvals(i,j) + elestr(i,j,k)                        
          end do                                                                
        end do                                                                  
       end do                                                                   
c                                                                               
      end if                                                                    
c                                                                               
c                       3. compute average value of each strain/stress          
c                          value for each element.                              
c                                                                               
      rngp = one / dble(ngp)                                                          
      temvals = temvals * rngp   ! entire array                                    
c                                                                               
c                       4. put average values within element at every           
c                          gauss point of element                               
c                                                                               
      if( ngp .eq. 8 ) then                                                     
c                                                                               
         do k = 1, 8                                                            
           do j = 1, num_vals                                                   
!$omp simd
             do i = 1, span                                                     
                elestr(i,j,k) = temvals(i,j)                                    
             end do                                                             
           end do                                                               
         end do                                                                 
c                                                                               
      else                                                                      
c                                                                               
         do k = 1, ngp                                                          
           do j = 1, num_vals                                                   
!$omp simd
             do i = 1, span                                                     
                elestr(i,j,k) = temvals(i,j)                                    
             end do                                                             
           end do                                                               
         end do                                                                 
c                                                                               
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
