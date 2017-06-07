c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine oupele                       *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 1/19/2017 rhd              *          
c     *                                                              *          
c     *     copy strain/stress results at the                        *          
c     *     element center into the data array for all element       *          
c     *     results.                                                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine oupele( span, num_strain, num_stress, do_stresses,             
     &                   elem_results, nrowd )                                  
      use elblk_data, only : elestr                                             
      implicit none                                                             
c                                                                               
      include 'param_def'                                                       
c                                                                               
      integer :: span, num_strain, num_stress, nrowd                            
      logical :: do_stresses                                                    
      double precision :: elem_results(nrowd,*)                                 
c                                                                               
      integer :: num_vals, k, i                                                 
c                                                                               
c                       copy element results into global                        
c                       data structure for element center results.              
c                                                                               
      num_vals = num_strain                                                     
      if( do_stresses ) num_vals = num_stress                                   
c                                                                               
      do k = 1, num_vals                                                        
        do i = 1, span                                                          
          elem_results(i,k) =  elestr(i,k,1)                                    
        end do                                                                  
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
