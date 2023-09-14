c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine gurson_cut_step                 *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 9/10/2023 rhd               *          
c     *                                                              *          
c     *         Checks if the maximum porosity change in an element  *
c     *         is > target limit or < target limit.                 *          
c     *         set new global load factor if we need to decrease or *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine gurson_cut_step
c                                      
      use global_data, only : out, noelem
      use elem_extinct_data, only : gt_old_porosity, dam_state                     
      use damage_data, only : all_elems_killed, dam_ptr, mxstp_store,
     &                        del_poros, 
     &                        max_porosity_change, load_reduced,
     &                        perm_load_fact, num_steps_min, 
     &                        no_killed_elems   
      use constants    
c            
      implicit none
c                                 
c           local declarations                                                  
c          
      integer :: elem, elem_ptr, i, nn                                                                    
      double precision :: max_del_poros, values(5)    
      logical :: debug, not_cut                                  
c   
      debug = .false.                                                                          
      if( debug ) write(out,*) '>>>>>> in gurson_cut_step'                     
c                                                                               
c           if all elements in the model have been killed, we                   
c           can skip processing here                                            
c                                                                               
      if( all_elems_killed ) return                                             
c                                                                               
c           loop over all killable gurson elements to find maximum              
c           change in porosity over last step. skip elements                    
c           that (1) are not killable, (2) have already been killed.            
c                                                                               
      max_del_poros = zero                                                      
c 
c$OMP PARALLEL DO  PRIVATE( elem, elem_ptr, values ) 
c$OMP&             reduction(max: max_del_poros) 
c                                                                        
      do elem = 1, noelem    
         elem_ptr = dam_ptr( elem )                                             
         if( elem_ptr .eq. 0 ) cycle                                            
         if( dam_state(elem_ptr) .ne. 0 ) cycle                                 
         call mm_return_values( "avg_porosity", elem, values, 5 )     
c                                                                               
c              find the change in porosity over the last step for this          
c              element -- if the change is the largest so far, store it.        
c              also store the porosity for comparison next step.                
c                                                                               
         max_del_poros = max( max_del_poros,                                    
     &               values(1) - gt_old_porosity(elem_ptr ) )         
         gt_old_porosity(elem_ptr) =  values(1)                                  
c                                                                               
      end do ! elem     
c$OMP END PARALLEL DO       
c                                                       
      if( debug ) then
         write(out,9100)  max_del_poros 
         write(out,*) ' no_killed_elems: ',no_killed_elems
      end if       
c                                                                               
c           store the maximum porosity in a data structure which            
c           holds the porosity values last "mxstp_store" steps                  
c                                                                               
      do i = mxstp_store, 2, -1                                                 
         del_poros(i) = del_poros(i-1)                                          
      end do                                                                    
      del_poros(1) = max_del_poros                                              
c                                                                               
c           if we have not yet killed any elements, do not allow                
c           the load control mechanism to change the load step size.            
c                                                                               
      if( no_killed_elems ) go to 9999                                          
c                                                                               
c           based on past porosities, evaluate load control and return          
c           with a load multiplier for the next step.                           
c                                                                               
      call gurson_load_factor( del_poros, mxstp_store,                          
     &     max_porosity_change, load_reduced, perm_load_fact,                   
     &     num_steps_min, no_killed_elems, out )
c                                                                               
 9999 continue                                                                  
      if( debug ) write( out, * ) '<<<<< leaving gurson_cut_step'               
c                                                                               
      return  
c
 9100 format(/,">>>>> max porosity change: ",f10.5)
c                                                                  
      end subroutine gurson_cut_step                                                                
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine gurson_load_factor              *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 9/10/23 rhd                *          
c     *                                   tfactor set back to 2.0    *          
c     *                                   from 4.0 to comply with    *          
c     *                                   manual                     *          
c     *                                                              *          
c     *         This subroutine modifies the load factor for         *          
c     *         gurson crack growth analyses.  The user specifies    *          
c     *         a target for the maximum increase in porosity per    *          
c     *         load step. If the actual increase in porosity is     *          
c     *         20% larger than that, then this algorithm computes   *          
c     *         a new, smaller load factor to reduce the loading.    *          
c     *         If the change in porosity during each of the last    *          
c     *         "max_steps_min" steps is 20% less than the target    *          
c     *         value, then increase the loading.  This allows the   *          
c     *         analysis to automatically adjust the load step size  *          
c     *         to prevent convergence problems, etc.                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine gurson_load_factor( del_poros, mxstp_store,                    
     &     max_porosity_change, load_reduced, perm_load_fact,                   
     &     num_steps_min, no_killed_elems, out )
      use constants                      
      implicit none                                                    
c
      integer, intent(in) :: out, mxstp_store
      integer, intent(out) :: num_steps_min
      logical, intent(in) ::  no_killed_elems  
      logical, intent(out):: load_reduced 
      double precision, intent(in) :: del_poros(*), max_porosity_change     
      double precision, intent(out) :: perm_load_fact              
c                                                                               
c           local declarations                                                  
c         
      integer :: j                                                                      
      integer, parameter :: max_steps_min = 3    
      logical :: debug                                         
      double precision :: ave_del_poros, ratio, tfactor                                                         
      double precision, parameter ::  max_factor = 1.05d0,
     &                                min_factor = 0.8d0            
c       
      debug = .false.                                                                       
      if( debug ) write (out,*) '   >>>>>> in gurson_load_factor'            
c                                                                               
c         evaluate ratio of the max change in porosity for last step over       
c         the maximum requested change in porosity. If the ratio is             
c         greater than 1.0, then the porosity increased more than               
c         requested.  If the ratio is less than 1.0, then the ratio is          
c         less than requested.                                                  
c                                                                               
      ratio = del_poros(1) / max_porosity_change                                
      if( debug ) write (*,*) '      >> ratio:',ratio                           
c                                                                               
c         if the ratio is more than max_factor, then reduce the                 
c         loading.  Estimate the reduction within 20 percent                    
c         of the user-requested porosity change value.                          
c                                                                               
c         note: These expressions currently assume that the gurson growth       
c         does not have overshoot control. If it does, then that needs to       
c         be factored into the load factor equation.                            
c                                                                               
c         overshoot control has not been implemented for gurson growth.         
c                                                                               
c         load_reduced is a global factor initialized to false on analysis      
c         startup. once set .true. it remains unchanged.                        
c                                                                               
      if( ratio .gt. max_factor ) then                                         
         load_reduced   = .true.                                                
         write(out,9338) del_poros(1), max_factor,
     &      max_porosity_change * max_factor,         
     &                   perm_load_fact, point_eight / ratio,                   
     &                   (perm_load_fact / ratio) * point_eight 
         perm_load_fact = ( perm_load_fact / ratio ) * point_eight              
      end if                                                                    
c                                                                               
c         to keep the load control algorithm from affecting the initial         
c         stages of the solution, we prevent any load reductions from           
c         taking place until the load has been reduced at least once.           
c                                                                               
      if( .not. load_reduced ) go to 9999                                       
c                                                                               
c         if the ratio is less than min_factor, and the ratio for the           
c         previous two steps has also been less than min_factor, then           
c         increase the loading.  Use the average of the ratios for the          
c         last three steps to approximate what the new loading should be.       
c         Only increase the loading if it hasn't been increased during the      
c         last three steps. The max load increase is a factor of 2.             
c                                                                               
c         note: These expressions currently assume that the gurson model        
c         does not have overshoot control. If it does, then that needs to       
c         be factored into the load factor equation.                            
c                                                                               
      if( ratio .lt. min_factor ) then                                          
         num_steps_min = num_steps_min + 1                                      
         if( num_steps_min .eq. max_steps_min ) then                            
c                                                                               
c               compute average change in porosity for the last                 
c               three steps                                                     
c                                                                               
            ave_del_poros = zero                                                
            do j = 1, num_steps_min                                             
               ave_del_poros = ave_del_poros + del_poros(j)                     
            end do                                                              
            ave_del_poros = ave_del_poros / num_steps_min                       
c                                                                               
c               compute load factor                                             
c                                                                               
            if( abs(ave_del_poros) .le. 1.0e-10 ) then                         
              tfactor = two                                                     
            else                                                                
              tfactor = one / ( ave_del_poros / max_porosity_change )           
            end if                                                              
            if( tfactor .gt. two ) tfactor = two                               
            write(out,9339) ave_del_poros,                                      
     &                      max_factor * max_porosity_change,                   
     &                      tfactor, perm_load_fact,                            
     &                      perm_load_fact*tfactor                              
            perm_load_fact = perm_load_fact * tfactor                           
            num_steps_min = 0                                                   
         end if                                                                 
      else                                                                      
         num_steps_min = 0                                                      
      end if                                                                    
c                                                                                                                                                             
 9999 continue                                                                  
      if( debug ) write( out, * ) '   <<<<< leaving gurson_load_factor'         
c                                                                               
      return                                                                    
c  
 9338 format(                                                                   
     &/1x,'>>>>> Note: maximum delta-f:',f8.4,' over last step exceeds ',       
     & f5.2,' x allowed delta-f of',f7.3,
     &/12x,' Current step size:',          
     & f7.3, ' reduced by factor of:',f6.2,' to:',f6.2)                
                                                                            
 9339 format(                                                                   
     &/1x,'>>>>> Note: the average delta-f: ',f8.4,' over last 3',              
     &/,13x,'steps is smaller that target value: ',f8.4,'. Step size',          
     &/,13x,'increased by a factor of: ',f5.2,' x the',                         
     &/,13x,'current value:',f7.2,' to a new value:',f7.2,/)                    
                                                                                
      end                                                                       
