c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine cohes_cut_step                  *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 9/3/2010 RHD               *          
c     *                                                              *          
c     *         This routine checks if the load step size is too     *          
c     *         large for analyses using cohesive elements to model  *          
c     *         crack growth. If the criterion in any killable       *          
c     *         element has grown more than the maximum allowed      *          
c     *         percent of the critical value then adjust load step  *          
c     *         size down. if change remains very small, we can      *          
c     *         adjust up the load step size                         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine cohes_cut_step( debug )                                        
      use global_data ! old common.main
      use elem_extinct_data, only : old_deff, dam_state                         
      use damage_data                                                           
      implicit integer (a-z)                                                    
c                                                                               
      logical debug                                                             
c                                                                               
c           local declarations                                                  
c                                                                               
      double precision                                                          
     &     values(20), two, dumd1, dumd2, dumd3, dumd4, zero,                   
     &     max_del_deff, new_deff, deff, peak_displ,                            
     &     new_deff_normalized,  max_del_deff_normalized                        
      logical not_cut, duml, all_killed, option_exponential,                    
     &        option_ppr, local_debug                                           
      character(len=1) :: dums                                                  
      real dumr                                                                 
      data two, zero / 2.0, 0.0 /                                               
c                                                                               
      local_debug = .false.                                                     
      if( local_debug ) write(out,9800)                                         
c                                                                               
c           if all cohesive elements in the model have been killed, we          
c           can skip processing here                                            
c                                                                               
      if ( all_elems_killed ) return                                            
c                                                                               
c           loop over all killable interface-cohesive elements                  
c           to find the maximum change in adaptive criterion.                   
c           skip elements with no "damage" yet.                                 
c                                                                               
c           Adaptive criterion is usually the change                            
c           in an effective interface displacement relative to the              
c           value of that displacement at the peak traction * a user            
c           specified input fraction (default = 0.2).                           
c                                                                               
c           here, we find if new load step size                                 
c           needs to be increased or decreased. skip elements                   
c           that: (1) are not killable, (2) have already been killed.           
c                                                                               
      all_killed    = .true.                                                    
      max_del_deff  = zero                                                      
      max_del_deff_normalized = zero                                            
c                                                                               
      do elem = 1, noelem                                                       
         elem_ptr = dam_ptr( elem )                                             
         if( elem_ptr .eq. 0 ) cycle                                            
         if( dam_state(elem_ptr) .ne. 0 ) cycle                                 
         all_killed = .false.                                                   
         cohes_type  = iprops(27,elem)                                          
         option_exponential = cohes_type .eq. 4                                 
         option_ppr         = cohes_type .eq. 6                                 
c                                                                               
c              calculate current state of "damage" in the                       
c              interface-cohesive element                                       
c                                                                               
         call dam_param_cohes( elem, duml, debug, values, 1 )                   
c                                                                               
c              compute value(s) to determine if the just converged              
c              load step size was too large.                                    
c                                                                               
c              deff = an "effective" scalar measure of the interface            
c                     displacement for the cohesive element                     
c                     (see dam_param_cohes for conversion of multiple           
c                      sliding and normal displacements to a single             
c                      "effective" value)                                       
c                                                                               
c              peak_displ = the value of deff at some measure of peak           
c                     value on the traction separation curves.                  
c                                                                               
c              for either pure opening or pure shear loading, these             
c              quantities are very clear. for mixed mode loading,               
c              dam_parm_cohes generates single "effective" values for           
c              use here in the adaptive load stepping (to make this             
c              adaptive algorithm reasonably straightforward).                  
c                                                                               
c              the global data vector old_deff( ) stores normalized             
c              values as defined below.                                         
c                                                                               
         if( option_exponential ) then                                          
            deff         = values(6)                                            
            peak_displ   = values(7)                                            
         elseif( option_ppr ) then                                              
            deff         = values(11)                                           
            peak_displ   = values(12)                                           
            if( local_debug )                                                   
     &          write(out,*) '.. elem, deff, peak_displ: ',                     
     &                        elem, deff, peak_displ                            
         endif                                                                  
c                                                                               
         new_deff_normalized     = deff / peak_displ                            
         max_del_deff_normalized = max( max_del_deff_normalized,                
     &                  new_deff_normalized - old_deff(elem_ptr) )              
         old_deff(elem_ptr) = new_deff_normalized                               
c                                                                               
      end do                                                                    
      if( local_debug ) write(out,9900) max_del_deff_normalized*100.0           
c                                                                               
c           now store the maximum increment in the (effective) interface        
c           displacement for all cohesive elements over the last step --        
c           holds the deff values last "mxstp_store" steps. Just                
c           converged step becomes first entry in vector.                       
c                                                                               
      do i = mxstp_store, 2, -1                                                 
         del_deff(i) = del_deff(i-1)                                            
      end do                                                                    
      del_deff(1) = max_del_deff_normalized                                     
c                                                                               
c           decide if the load step size should remain the same,                
c           be increased or decreased. this can be done from the                
c           very beginning of the analysis or only after the first              
c           cohesive element has been killed. uncomment the appropriate         
c           if stm. below to implement desired behavior.                        
c                                                                               
c      if( no_killed_elems ) go to 9999                                         
      if( all_killed )      go to 9999                                          
c                                                                               
c           based on past deff values, evaluate load control and return         
c           with a load multiplier for the next step. max_deff_change           
c           was input by user (default=0.2). it is used to assess the           
c           normalized change in (effective) separation displacement            
c           [computed/stored above] for the just computed step.                 
c                                                                               
      call cohes_load_factor( del_deff, mxstp_store,                            
     &     max_deff_change, load_reduced, perm_load_fact,                       
     &     num_steps_min, no_killed_elems, out )                                
c                                                                               
 9999 continue                                                                  
      if( local_debug ) write(out,*) "<<<<< leaving cohes_cut_step"             
c                                                                               
      return                                                                    
c                                                                               
 9800 format(//,">>>>>> in cohes_cut_step")                                     
 9900 format(                                                                   
     &  2x,".... max change in (effective) displacement jump for"               
     &        " adaptive",                                                      
     &/,2x,".... as % of (effective) jump displacement",                        
     &/,2x,".... at peak stress: ", f5.1)                                       
c                                                                               
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine cohes_load_factor               *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 9/3/2010 RHD               *          
c     *                                                              *          
c     *         This subroutine modifies the load factor for         *          
c     *         cohesive crack growth analyses. The user specifies   *          
c     *         a target for the maximum increase in effective       *          
c     *         interface displacement per load step. If the actual  *          
c     *         increase in this value is 5% larger than that, then  *          
c     *         this algorithm computes                              *          
c     *         a new, smaller load factor to reduce the loading.    *          
c     *         If the change in effective interface displacement    *          
c     *         during each of the last                              *          
c     *         "max_steps_min" steps is 20% less than the target    *          
c     *         value, then increase the loading.  This allows the   *          
c     *         analysis to automatically adjust the load step size  *          
c     *         to prevent convergence problems, etc.                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine cohes_load_factor( del_deff, mxstp_store,                      
     &     max_deff_change, load_reduced, perm_load_fact,                       
     &     num_steps_min, no_killed_elems, out )                                
      implicit integer (a-z)                                                    
c                                                                               
      logical load_reduced, no_killed_elems                                     
      parameter (max_steps_min = 3)                                             
c                                                                               
c           local declarations                                                  
c                                                                               
      double precision                                                          
     &     two, dumd1, dumd2, dumd3, dumd4, max_factor, min_factor,             
     &     del_deff(*), max_deff_change, zero, one, point_eight,                
     &     ave_del_deff, ratio, perm_load_fact, tfactor, four                   
      logical duml, local_debug                                                 
      character(len=1) :: dums                                                  
      real dumr                                                                 
      data two, max_factor, min_factor, one, point_eight, zero, four            
     &     / 2.0, 1.01, 0.5, 1.0, 0.8, 0.0, 4.0 /                               
c                                                                               
c                                                                               
      local_debug = .false.                                                     
      if( local_debug ) write (out,*) '   >>>>>> in cohes_load_factor'          
c                                                                               
c         del_eff(1) -- change in (effective) jump displacement                 
c         for last step [max over all cohesive elements] normalized             
c         by value of (effective) jump displacement at some                     
c         measure of "peak" value on traction-separation curve.                 
c                                                                               
c         compare to maximum normalized change requested by user.               
c         If the ratio is >1.0, then the relative displacement                  
c         increased more than requested. If the ratio is < 1.0,                 
c         then the ratio is < requested.                                        
c                                                                               
      ratio      =  del_deff(1) / max_deff_change                               
      if( local_debug ) write (out,9900) del_deff(1),                           
     &                  max_deff_change, ratio                                  
c                                                                               
c         if the ratio is more than max_factor, then reduce the                 
c         loading. Estimate the reduction within 20 percent                     
c         of the user-requested change in relative displacement.                
c                                                                               
c         note: These expressions currently assume that cohesive growth         
c         does not have "overshoot" control (use projections of changes         
c         in next step to pre-adjust load size). If it does, then               
c         that needs to be factored into the load factor equation.              
c                                                                               
      if ( ratio .gt. max_factor ) then                                         
         perm_load_fact = perm_load_fact / ratio * point_eight                  
         load_reduced   = .true.                                                
         write(out,9983) del_deff(1), max_deff_change, perm_load_fact           
      end if                                                                    
c                                                                               
c         to keep the load control algorithm from affecting the initial         
c         stages of the solution, we prevent any load increases from            
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
c         note: These expressions currently assume that cohesive growth         
c         does not have overshoot control. If it does, then that needs to       
c         be factored into the load factor equation.                            
c                                                                               
      if( ratio .lt. min_factor ) then                                          
         num_steps_min = num_steps_min + 1                                      
         if( num_steps_min .eq. max_steps_min ) then                            
c                                                                               
c               compute average change in relative displacement                 
c               for the last three steps                                        
c                                                                               
            ave_del_deff = zero                                                 
            do j = 1, num_steps_min                                             
               ave_del_deff = ave_del_deff + del_deff(j)                        
            end do                                                              
            ave_del_deff = ave_del_deff / num_steps_min                         
c                                                                               
c               compute load factor                                             
c                                                                               
            if ( abs(ave_del_deff) .le. 1.0e-10 ) then                          
              tfactor = two                                                     
            else                                                                
              tfactor = one / ( ave_del_deff / max_deff_change )                
            end if                                                              
            if ( tfactor .gt. two ) tfactor = two                               
            perm_load_fact = perm_load_fact * tfactor                           
            write(out,9337) ave_del_deff, perm_load_fact                        
            num_steps_min = 0                                                   
         end if                                                                 
      else                                                                      
         num_steps_min = 0                                                      
      end if                                                                    
c                                                                               
c                                                                               
 9999 continue                                                                  
      if( local_debug )                                                         
     &   write(out,*) '      <<<<< leaving cohesive_load_factor'                
c                                                                               
      return                                                                    
c                                                                               
 9983 format(                                                                   
     & /1x,'>>>>> Adaptive solution for cohesive crack growth ',                
     &            'reducing load step size',/,                                  
     & 8x,'-> max increment of (effective) jump displacement ',                 
     &                  'over last step normalized by ',/                       
     & 14x,'(effective) displacement jump at peak stress: ',f11.8,/             
     & 8x,'-> exceeds user-specifed value: ',f5.2,                              
     &                       ' by more than 5%.',/,                             
     & 8x,'-> future load step sizes reduced by ',                              
     &                   'factor: ',f7.3,/)                                     
 9337        format(                                                            
     & /1x,'>>>>> Adaptive solution for cohesive crack growth ',                
     &            'increasing load step size',/,                                
     & 8x,'-> avg max increment of (effective) jump displacement ',             
     &                  'over last 3 steps normalized by ',/                    
     & 14x,'(effective) displacement jump at peak stress: ',f11.8,/,            
     & 14x,'is smaller than target value by more than 20%. ',/                  
     & 8x,'-> step size increased to ',f7.3,' x the original',                  
     &                   ' loading.',/)                                         
 9900 format(2x,".... del_eff(1):                ",e14.6,                       
     & /,    2x,".... max_deff_change:           ",e14.6,                       
     & /,    2x,".... ratio:                     ",e14.6)                       
c                                                                               
      end                                                                       
                                                                                
                                                                                
                                                                                
