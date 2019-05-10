c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine smcs_cut_step                   *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 4/28/2019 rhd              *          
c     *                                                              *          
c     *         This routine checks if the load step size is too     *          
c     *         large. If the plastic strain in any killable mises   *          
c     *         element has grown more than max_plast_strain_change, *          
c     *         then permanently cut the load step size by half.     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine smcs_cut_step( debug )                                        
      use global_data, only : noelem, out   ! old common.main
      use elem_extinct_data, only : old_plast_strain                            
      use damage_data, only : dam_ptr, max_plast_strain_change, 
     &                        perm_load_fact         
c                                                 
      implicit none                                                    
      logical :: debug                                                                                                                                            
c                                                                               
c           local declarations                                                  
c                              
      integer :: elem, elem_ptr, dum                                                 
      logical :: not_cut, kill_now                                                     
      double precision :: eps_plas, eps_crit, sig_mean, sig_mises,
     &                    triaxiality, new_plast_strain, mean_zeta   
      double precision, parameter :: two = 2.0d0                
      character(len=1) :: dums                                                  
      real :: dumr                                                                 
c                                                                               
      if ( debug ) write ( out, * ) '>>>>>> in smcs_cut_step'                   
c                                                                               
c              loop over all killable elements                               
c                                                                               
      not_cut = .true.   
                                                             
      do elem = 1, noelem                                                       
        elem_ptr = dam_ptr( elem )                                             
        if( elem_ptr .eq. 0 ) cycle   ! element not killable                                        
c                                                                               
c              get plastic eps in element (=2 means no history update)                            
c                                                                               
        call dam_param_3_get_values( elem, debug, eps_plas, eps_crit, 
     &                               sig_mean, sig_mises, triaxiality, 
     &                               mean_zeta, 2, kill_now )   
        new_plast_strain = eps_plas   
c                                                                               
c              compare old plast_strain with new plast_strain -- if             
c              change is larger than the acceptible max, cut load               
c              step size.                                                       
c                                                                               
        if( new_plast_strain - old_plast_strain(elem_ptr) .gt.                
     &        max_plast_strain_change ) then                                    
          if( not_cut ) then                                                   
               perm_load_fact = perm_load_fact / two                            
               call errmsg ( 273, dum, dums, dumr, perm_load_fact )             
             not_cut = .false.                                                  
            end if                                                              
         end if                                                                 
         old_plast_strain(elem_ptr) = new_plast_strain                            
      end do                                                                    
c                                                                               
      if( debug ) write(out,*) '<<<<< leaving smcs_cut_step'                   
c                                                                               
      return                                                                    
      end                                                                       
