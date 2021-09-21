c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine smcs_cut_step                   *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 4/13/21 rhd                *          
c     *                                                              *          
c     *         adaptively control the global load step size based   *          
c     *         on increments of plastic strain between load steps   *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine smcs_cut_step( debug )                                        
      use global_data, only : noelem, out   ! old common.main
      use elem_extinct_data, only : old_plast_strain, dam_state  
      use damage_data, only : dam_ptr,
     &           allowable_change => max_plast_strain_change, 
     &           perm_load_fact, num_elements_in_force_release,
     &           smcs_adapt_alpha_min, smcs_adapt_alpha_max,    
     &           smcs_allowable_in_release
      use dam_param_code, only :  dam_param_3_get_values
      use constants
c                                                 
      implicit none                                                    
      logical :: debug                                                                                                                                            
c                                                                               
c              local declarations                                                  
c                              
      integer :: elem, elem_ptr, count, elem_max_change, 
     &           allowed_elems_in_release                                              
      logical :: ldebug, increase_size, ldummy                                         
      double precision :: sig_mean, sig_mises, triaxiality, 
     &                    new_plast_strain, mean_zeta, mean_bar_theta,
     &                    max_change, old_perm_factor, eps_crit, factor,
     &                    default_cut_factor, mean_omega, sig_1,
     &                    tear_param, increase_multiplier, delta 
      double precision, parameter :: one_toler = 0.999d0
c                
      ldebug = debug
c
      if( ldebug ) write ( out, * ) '>>>>>> in smcs_cut_step'                   
      max_change = zero 
      old_perm_factor = perm_load_fact
      count = 0
c
      default_cut_factor       = smcs_adapt_alpha_min
      increase_multiplier      = smcs_adapt_alpha_max
      allowed_elems_in_release = smcs_allowable_in_release
      elem_max_change          = 0
c
c              get plastic eps in each killable element that has not
c              already been killed (=2 means no history update)                            
c                                                                               
      do elem = 1, noelem                                                       
        elem_ptr = dam_ptr( elem )                                             
        if( elem_ptr .eq. 0 ) cycle   ! element not killable 
        if( dam_state(elem_ptr) .gt. 0 ) cycle ! elem killed already     
c                                              or being killed                
        call dam_param_3_get_values(
     &      elem, debug, new_plast_strain, eps_crit, sig_mean, 
     &      sig_mises, triaxiality, mean_zeta, mean_omega,
     &      mean_bar_theta, 2, tear_param, ldummy, sig_1 )   
        delta = new_plast_strain - old_plast_strain(elem_ptr) 
        if( delta > max_change ) then
           max_change = delta
           elem_max_change = elem
        end if
        old_plast_strain(elem_ptr) = new_plast_strain  
        count = count + 1
      end do ! over elements 
c
c              is max change too large, too small, ok  ? 
c
      if( count == 0 ) return ! all elements killed
      if( max_change <= zero ) return ! no eps_plas 
      if( ldebug ) write(out,9050) max_change, allowable_change           
c
c              rule for plastic strain increment too large:
c              reduce global load factor by at
c              least * default_cut_factor,
c              or actual reduction to match allowable if smaller
c 
      if( max_change > allowable_change ) then   ! too big
        factor = min( default_cut_factor,  
     &                allowable_change / max_change )
        perm_load_fact = perm_load_fact * factor
        write(out,9000) max_change, elem_max_change,
     &                  allowable_change, factor, perm_load_fact 
        if( ldebug ) write(out,*) '<<<<< leaving smcs_cut_step'  
        return                 
      end if
c
c              load reduction not required based on plastic strain
c              increment.
c
c              may want to reduce future step sizes if the current
c              number of elements is release exceeds a threshold.
c
      if( num_elements_in_force_release >
     &    allowed_elems_in_release ) then
        factor = default_cut_factor
        perm_load_fact = perm_load_fact * factor
        write(out,9010) num_elements_in_force_release, 
     &                  allowed_elems_in_release,
     &                  factor, perm_load_fact 
        if( ldebug ) write(out,*) '<<<<< leaving smcs_cut_step'  
        return                 
      end if
 
c
c              may be able to increase global loading given
c              small plastic strain increments
c              rules:
c              - if global load multiplier = 1.0, no increase
c              - if elements are beiing released, no increase
c              - multiplier is smaller of 1.5 and actual value to
c                match allowable strain increment.
c              - in no case can perm_load_factor exceed 1.0
c              so we never more than double the global load factor
c              in a single step.
c
      increase_size = max_change < allowable_change
      if( perm_load_fact > one_toler ) increase_size = .false.
      if( num_elements_in_force_release > 0 ) increase_size = .false.
      if( increase_size ) then   ! too small
        factor = allowable_change / max_change 
        if( factor > increase_multiplier ) factor = increase_multiplier
        if( perm_load_fact * factor > one ) then
           perm_load_fact = one
        else
           perm_load_fact = perm_load_fact * factor
        end if
        write(out,9100) max_change, elem_max_change, allowable_change,
     &                  factor, perm_load_fact 
      end if
c                                                                               
      if( ldebug ) write(out,*) '<<<<< leaving smcs_cut_step'                   
c                                                                               
      return   
c
 9000 format(">>> crack growth processor decreasing global ",
     & "load step size:",
     & /,5x,"max change plastic strain over step:",f15.8,2x,
     &  "element: ",i8, /,                                
     & 5x,"user-specified allowable strain change:",f15.8,    
     & /,5x,"step-size (decrease) factor multiplier: ",f10.5,
     & /,5x,"updated global step multiplier: ",e13.6,//)
 9010 format(">>> crack growth processor decreasing global ",
     & "load step size:",
     & /,5x,"max # elements in release: ",i4,
     & 5x,"exceeds allowable of: ",i4,    
     & /,5x,"step-size (decrease) factor multiplier: ",f10.5,
     & /,5x,"updated global step multiplier: ",e13.6,//)
 9100 format(">>> crack growth processor increasing global ",
     & "load step size:",
     & /,5x,"max change plastic strain over step:    ",f15.8,2x, 
     &  "element: ",i8, /,                                
     & 5x,"user-specified allowable strain change: ",f15.8,    
     & /,5x,"step-size (increase) factor multiplier: ",f10.5,
     & /,5x,"updated global step multiplier: ",e13.6,//)
 9050 format(5x,"**** max change, allocable change:",f15.8)
c                                       
            end        
