c ********************************************************************          
c *                                                                  *          
c *  routines to support crack growth by cohesive criterion          *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                subroutine dam_param_cohes                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 7/31/2013 rhd              *          
c     *                                                              *          
c     *     for a killable cohesive element not yet killed,          *          
c     *     determine if the element should be killed now.           *          
c     *     return values of cohesive parameters to support optional *          
c     *     printing of values at start of each step                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine dam_param_cohes( elem, kill_now, debug, values,                
     &                            compute_type )                                
      use global_data ! old common.main
c                                                                               
      use main_data,       only : elems_to_blocks                               
      use elem_block_data, only : history_blocks, eps_n_blocks,                 
     &                            urcs_n_blocks, history_blk_list               
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c              parameter declarations                                           
c                                                                               
      logical kill_now, debug                                                   
      double precision                                                          
     &     values(*)  ! see caller for size                                     
c                                                                               
c              local declarations                                               
c                                                                               
      double precision                                                          
     &     zero, one, half, t1, t2, tn, d1, d2, dn, fpngp, beta,                
     &     deff, deff_peak, teff, normalized_teff, teff_peak,                   
     &     normalized_deff, peak_normal_stress,                                 
     &     peak_shear_stress, dn_max_hist, ds_max_hist,                         
     &     dn_limit, ds_limit, dn_at_peak, ds_at_peak,                          
     &     ratio_normal, ratio_shear, dn_half_peak, ds_half_peak,               
     &     kill_fract, deff_ratio, deff_gp, deff_peak_gp                        
      real chk_value                                                            
      double precision,                                                         
     &  dimension(:), pointer :: history, urcs_n, eps_n                         
      logical beyond_peak, fgm_cohes, option_ppr,                               
     &        option_exponential, kill_criterion_element,                       
     &        kill_gp, local_debug, option_cavit                                
      data zero, one, half                                                      
     &       / 0.0d00, 1.0d00, 0.5d00 /                                         
c                                                                               
c                                                                               
c               elem -- absolute number for interface element                   
c               kill_now -- output, logical                                     
c               debug -- input, logical                                         
c               values -- output vector, float                                  
c               compute_type -- input integer.                                  
c                    = 1 compute/return the "values" vector                     
c                        but not kill_now flag (for checking,                   
c                        printing, etc.)                                        
c                    = 2 compute/return the "values" vector                     
c                        and the kill_now flag                                  
c                                                                               
c                                                                               
c               get basic information about element and the                     
c               cohesive material model for the element. set pointers           
c               into block data structures to simplify access                   
c               to stress, strain, history values for this specific             
c               element.                                                        
c                                                                               
c               cohesive types:                                                 
c                  1 -- linear elastic; 2 -- bilinear; 3 -- ramp                
c                  4 -- exponential_1; 5 -- exponential_2                       
c                  6 -- PPR; 7 -- cavit                                         
c                  Only types 1, 4, 6, 7 implemented.                           
c                                                                               
      local_debug = .false.                                                     
      if( local_debug ) write(out,9030)                                         
      ngp         = iprops(6,elem)                                              
      blk         = elems_to_blocks(elem,1)                                     
      rel_elem    = elems_to_blocks(elem,2)                                     
      hist_size   = history_blk_list(blk)                                       
      hoffset     = (rel_elem-1)*hist_size*ngp                                  
      epsoffset   = (rel_elem-1)*nstr*ngp                                       
      sigoffset   = (rel_elem-1)*nstrs*ngp                                      
      urcs_n      => urcs_n_blocks(blk)%ptr                                     
      eps_n       => eps_n_blocks(blk)%ptr                                      
      history     => history_blocks(blk)%ptr                                    
      cohes_type  = iprops(27,elem)                                             
      option_exponential = cohes_type .eq. 4                                    
      option_ppr         = cohes_type .eq. 6                                    
      option_cavit       = cohes_type .eq. 7                                    
c                                                                               
c               loop over element integration pts. compute the simple           
c               average of the normal traction and (relative)                   
c               displacement. Also compute average of the two shear             
c               tractions and (total) sliding displacement                      
c                                                                               
      t1 = zero; t2 = zero; tn = zero; d1 = zero; d2 = zero; dn = zero          
c                                                                               
      do gp = 1, ngp                                                            
        t1 = t1 + urcs_n(sigoffset+1)                                           
        t2 = t2 + urcs_n(sigoffset+2)                                           
        tn = tn + urcs_n(sigoffset+3)                                           
        d1 = d1 + eps_n(epsoffset+1)                                            
        d2 = d2 + eps_n(epsoffset+2)                                            
        dn = dn + eps_n(epsoffset+3)                                            
        if( local_debug) then                                                   
         i = sigoffset; j = epsoffset                                           
          write(out,9020) gp, urcs_n(i+1), urcs_n(i+2), urcs_n(i+3),            
     &                    eps_n(j+1), eps_n(j+2), eps_n(j+3)                    
        end if                                                                  
        epsoffset = epsoffset + nstr                                            
        sigoffset = sigoffset + nstrs                                           
      end do                                                                    
c                                                                               
c    finish computing the average value for each of the                         
c    3 tractions and 3 interface displacements. compute                         
c    resultant shear traction and resultant sliding displacement.               
c                                                                               
c       values        content                                                   
c        1         averaged normal traction                                     
c        2         averaged total shear traction (1+2)                          
c        3         averaged normal opening                                      
c        4         averaged (total) sliding displacement (1+2)                  
c                                                                               
c                    for exponential option (=4)                                
c                    ===========================                                
c                                                                               
c        5         averaged single "effective" traction (reflects a             
c                  user-specified beta value)                                   
c        6         averaged single "effective" interface displacement (reflects 
c                  user-specified beta value)                                   
c        7         user-specified "equivalent" interface displacement           
c                  at peak stress                                               
c        8         averaged equivalent stress / user-specified peak stress      
c        9         averaged equivalent interface displacement / user-specified  
c                  equivalent displacement at peak stress                       
c                                                                               
c                    for PPR option (=6)                                        
c                    ===================                                        
c                                                                               
c        5         averaged normal traction / peak normal traction              
c                  on curve                                                     
c        6         averaged shear traction / peak shear traction                
c                  on curve                                                     
c        7         the number of Gauss points as follows                        
c                   for each GP: get maximum normal displacement                
c                   over loading history.                                       
c                   if it exceeds the normal displacement at peak               
c                   normal stress on curve, increment counter                   
c                   of GPs                                                      
c        8         the number of Gauss points as follows                        
c                   for each GP: get maximum sliding displacement               
c                   (1+2) over loading history.                                 
c                   if it exceeds the sliding displacement at peak shear        
c                   stress on curve, increment counter of GPs                   
c        9         the number of Gauss points as follows                        
c                   for each GP: get maximum normal displacement                
c                   over loading history.                                       
c                   if it exceeds the 50% of normal displacement at             
c                   peak normal stress on curve, increment GP counter           
c       10         the number of Gauss points as follows                        
c                   for each GP: get maximum shear displacement                 
c                   over loading history.                                       
c                   if it exceeds the sliding displacement at peak shear        
c                   stress on curve, increment counter                          
c       11          maximum value of an "effective" interface                   
c                   displacement at Gauss points                                
c                    = sqrt( dnormal**2 + dsliding**2 )                         
c       12          maximum value of effective peak displacement at Gauss points
c                   corresponding to above value. = sqrt( dn_peak**2 +          
c                   ds_peak**2 )                                                
c                                                                               
c                    for cavit option (=7)                                      
c                    ===================                                        
c                                                                               
c        5         averaged normal traction / peak normal traction              
c                  on curve                                                     
c        6         averaged shear traction / peak shear traction                
c                  on curve                                                     
c                                                                               
c                                                                               
      fpngp = one / dble(ngp)                                                   
      t1 = t1 * fpngp                                                           
      t2 = t2 * fpngp                                                           
      tn = tn * fpngp                                                           
      d1 = d1 * fpngp                                                           
      d2 = d2 * fpngp                                                           
      dn = dn * fpngp                                                           
      values(1) = tn                                                            
      values(2) = sqrt( t1*t1 + t2*t2 )                                         
      values(3) = dn                                                            
      values(4) = sqrt( d1*d1 + d2*d2 )                                         
c                                                                               
c              "exponential_1" option                                           
c                                                                               
      if( option_exponential ) then                                             
        call dam_param_cohes_exp                                                
        return                                                                  
      end if                                                                    
c                                                                               
c              "PPR" option                                                     
c                                                                               
      if( option_ppr ) then                                                     
        call dam_param_cohes_ppr                                                
        return                                                                  
      end if                                                                    
c                                                                               
c              "cavit" option                                                   
c                                                                               
      if( option_cavit ) then                                                   
        call dam_param_cohes_cavit                                              
        return                                                                  
      end if                                                                    
c                                                                               
c              unknown cohesive option                                          
c                                                                               
      write(out,*) ' '                                                          
      write(out,*) '>>> FATAL ERROR: dam_param_cohes (1)'                       
      write(out,*) '    unknown cohesive option. job terminated'                
      call die_abort                                                            
c                                                                               
 9020 format(                                                                   
     &   5x,'gp:                  ',i2,                                         
     & /,5x,'t1:                  ',f12.4,                                      
     & /,5x,'t2:                  ',f12.4,                                      
     & /,5x,'tn:                  ',f12.4,                                      
     & /,5x,'d1:                  ',f12.8,                                      
     & /,5x,'d2:                  ',f12.8,                                      
     & /,5x,'dn:                  ',f12.8 )                                     
 9030 format(/,"..... inside dam_param_cohes ....")                             
c                                                                               
c                                                                               
      contains                                                                  
c     ========                                                                  
c     ****************************************************************          
c     *                                                              *          
c     *                subroutine dam_param_cohes_exp                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 7/31/2013 rhd              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine dam_param_cohes_exp                                            
      implicit integer (a-z)                                                    
c                                                                               
      beta        = props(23,elem)                                              
      fgm_cohes   = abs(props(34,elem)) .gt. 0.1                                
      if( beta .eq. zero ) then                                                 
        values(4) = zero                                                        
        values(5) = tn                                                          
        values(6) = dn                                                          
      else                                                                      
        if( dn .lt. zero ) then                                                 
           tn = zero                                                            
           dn = zero                                                            
        end if                                                                  
        values(5) = sqrt( (t1*t1 + t2*t2)/beta/beta + tn*tn )                   
        values(6) = sqrt( beta*beta*(d1*d1+d2*d2) + dn*dn )                     
      end if                                                                    
      teff      = values(5)                                                     
      deff      = values(6)                                                     
c                                                                               
      values(7) = props(22,elem)  ! equivalent separation distance              
                                  ! at peak stress                              
      teff_peak = props(13,elem)  ! peak normal stress on expon curve           
      deff_peak = props(22,elem)                                                
      if( fgm_cohes ) then                                                      
         values(7) = props(35,elem) ! peak normal stress on curve (fgm)         
         teff_peak = props(37,elem) ! critical stress ductile (fgm)             
         deff_peak = props(35,elem)                                             
      end if                                                                    
c                                                                               
      normalized_teff = teff / teff_peak                                        
      normalized_deff = deff / deff_peak                                        
      values(8)       = normalized_teff                                         
      values(9)       = normalized_deff                                         
c                                                                               
c              does element need to be killed?                                  
c                                                                               
      kill_now    = .false.                                                     
      if( compute_type .ne. 2 ) return                                          
      deff_peak   = props(22,elem)   ! separation distance at peak on           
                                     ! exponential curve                        
      if( fgm_cohes )  deff_peak = props(35,elem)                               
      deff        = values(6)                                                   
      beyond_peak =  deff .gt. deff_peak                                        
      if ( .not. beyond_peak ) return                                           
      if ( normalized_deff .gt. critical_cohes_deff_fract )                     
     &      kill_now = .true.                                                   
c                                                                               
      return                                                                    
      end subroutine dam_param_cohes_exp                                        
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                subroutine dam_param_cohes_ppr                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 7/31/2013 rhd              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine dam_param_cohes_ppr                                            
      implicit integer (a-z)                                                    
c                                                                               
c              loop over integration pts of this interface element.             
c               o pull from history the maximum normal and sliding              
c                 displacements ever reached over loading history               
c               o pull from history the normal and sliding displacements        
c                 at which tractions degrade to zero                            
c               o pull from history the fractions of displacements              
c                 at peak values of normal and shear tractions                  
c               o count number of gauss points at which maximum                 
c                 interface displacements over history exceed 50% of            
c                 values at peak tractions                                      
c               o count number of gauss points at which maximum                 
c                 interface displacements over history exceed values            
c                 at peak tractions                                             
c                                                                               
c              Does the integration pt status meet criterion for element        
c              extinction? (All GP must meet criterion for extincting           
c              the element).                                                    
c                                                                               
c              After much discussion, Kyoungsoo recommends using the            
c              "or" condition on normal and shear loading. If either            
c              the normal opening or sliding displacement exceeds a             
c              fraction (e.g., 90%) of the corresponding displacement           
c              when the traction degrades fully to zero, kill the element       
c              (must be satisfied at all interface element Gauss points)        
c                                                                               
c              The reason is that the remaining shear capacity                  
c              for example, degrades severely as the normal                     
c              traction capacity degrades. Same for normal capacity             
c              under shear loading.                                             
c                                                                               
      epsoffset        = (rel_elem-1)*nstr*ngp                                  
      num_dn_past_peak = 0                                                      
      num_ds_past_peak = 0                                                      
      num_dn_half_peak = 0                                                      
      num_ds_half_peak = 0                                                      
      kill_fract = ppr_kill_displ_fraction                                      
      kill_criterion_element = .true.                                           
      deff_ratio = -1.0 ! initialize to small value                             
      values(11) = zero                                                         
      values(12) = one                                                          
      if( local_debug ) write(out,9000) elem                                    
c                                                                               
      do gp = 1, ngp                                                            
c                                                                               
c                  values for possible printing                                 
c                                                                               
        dn_limit    = history(hoffset+1)                                        
        ds_limit    = history(hoffset+2)                                        
        dn_max_hist = history(hoffset+3)                                        
        ds_max_hist = abs(history(hoffset+4))                                   
        ratio_normal = props(37,elem)                                           
        ratio_shear  = props(39,elem)                                           
        dn_at_peak = dn_limit * ratio_normal                                    
        ds_at_peak = ds_limit * ratio_shear                                     
        dn_half_peak = half * dn_at_peak                                        
        ds_half_peak = half * ds_at_peak                                        
        if( dn_max_hist .gt. dn_at_peak )                                       
     &                       num_dn_past_peak = num_dn_past_peak + 1            
        if( ds_max_hist .gt. ds_at_peak )                                       
     &                       num_ds_past_peak = num_ds_past_peak + 1            
        if( dn_max_hist .gt. dn_half_peak )                                     
     &                       num_dn_half_peak = num_dn_half_peak + 1            
        if( ds_max_hist .gt. ds_half_peak )                                     
     &                       num_ds_half_peak = num_ds_half_peak + 1            
c                                                                               
c                  does this Gauss point satisfy criterion to                   
c                  kill element? If not, element is not killable                
c                  at this load step                                            
c                                                                               
        kill_gp = ( dn_max_hist .gt. kill_fract * dn_limit )  .or.              
     &            ( ds_max_hist .gt. kill_fract * ds_limit )                    
        if( .not. kill_gp ) kill_criterion_element = .false.                    
c                                                                               
c                  values that may be used by adapative load step               
c                  processor                                                    
c                                                                               
        d1 = eps_n(epsoffset+1)                                                 
        d2 = eps_n(epsoffset+2)                                                 
        dn = eps_n(epsoffset+3)                                                 
        epsoffset = epsoffset + nstr                                            
        deff_gp = sqrt( dn**2 + d1**2 + d2**2 )                                 
        deff_peak_gp = sqrt( dn_at_peak**2 + ds_at_peak**2 )                    
        normalized_deff = deff_gp / deff_peak_gp                                
c        write(*,*) '...gp, d1, d2, dn, deff_gp, dn_at_peak,ds_at_peak,         
c     &     deff_peak_gp'                                                       
c        write(*,*) gp, d1, d2, dn, deff_gp, dn_at_peak,ds_at_peak,             
c     &       deff_peak_gp                                                      
        if( normalized_deff .gt. deff_ratio ) then                              
            values(11) = deff_gp                                                
            values(12) = deff_peak_gp                                           
            deff_ratio = normalized_deff                                        
        end if                                                                  
c                                                                               
        if( local_debug .and. gp .eq. 1 )                                       
     &       write(out,9010) gp, dn_limit, ds_limit,                            
     &       dn_max_hist, ds_max_hist, ratio_normal, ratio_shear,               
     &       dn_at_peak, ds_at_peak, dn_half_peak, ds_half_peak,                
     &       kill_gp, deff_gp, normalized_deff                                  
        hoffset   = hoffset + hist_size                                         
      end do                                                                    
c                                                                               
      peak_normal_stress = props(13,elem)                                       
      peak_shear_stress  = props(14,elem)                                       
      values(5) = zero; values(6) = zero                                        
      if( peak_normal_stress .gt. zero )                                        
     &    values(5) = values(1) / peak_normal_stress                            
      if( peak_shear_stress .gt. zero )                                         
     &    values(6) = values(2) / peak_shear_stress                             
c                                                                               
      values(7) = zero; values(8) = zero                                        
      values(9) = zero; values(10) = zero                                       
      if( num_dn_past_peak .gt. 0 ) then                                        
         values(7) = dble(num_dn_past_peak)                                     
      end if                                                                    
      if( num_ds_past_peak .gt. 0 ) then                                        
         values(8) = dble(num_ds_past_peak)                                     
      end if                                                                    
c                                                                               
      if( num_dn_half_peak .gt. 0 ) then                                        
         values(9) = dble(num_dn_half_peak)                                     
      end if                                                                    
      if( num_ds_half_peak .gt. 0 ) then                                        
         values(10) = dble(num_ds_half_peak)                                    
      end if                                                                    
c                                                                               
c              does element need to be killed? needs fixing ....!               
c                                                                               
      kill_now    = .false.                                                     
      if ( compute_type .ne. 2 ) return                                         
      kill_now =  kill_criterion_element                                        
      return                                                                    
c                                                                               
 9000 format('**** debug ppr in dam_param_cohes element: ',i6)                  
 9010 format(                                                                   
     &   5x,'gp:                  ',i2,                                         
     & /,5x,'dn_limit:            ',f12.8,                                      
     & /,5x,'ds_limit:            ',f12.8,                                      
     & /,5x,'dn_max_hist:         ',f12.8,                                      
     & /,5x,'ds_max_hist:         ',f12.8,                                      
     & /,5x,'ratio_normal:        ',f12.8,                                      
     & /,5x,'ratio_shear:         ',f12.8,                                      
     & /,5x,'dn_at_peak:          ',f12.8,                                      
     & /,5x,'ds_at_peak:          ',f12.8,                                      
     & /,5x,'dn_half_peak:        ',f12.8,                                      
     & /,5x,'ds_half_peak:        ',f12.8,                                      
     & /,5x,'kill_gp:             ',l1,                                         
     & /,5x,'deff_gp:             ',f12.8,                                      
     & /,5x,'normalized_deff:     ',f12.8)                                      
                                                                                
c                                                                               
      end subroutine dam_param_cohes_ppr                                        
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                subroutine dam_param_cohes_cavit              *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 7/31/2013 rhd              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine dam_param_cohes_cavit                                          
      implicit integer (a-z)                                                    
      double precision                                                          
     &     tn_max_over_history, dn_max_at_tn_max                                
c                                                                               
c              loop over integration pts of this interface element.             
c                                                                               
c               o pull from history the maximum normal traction                 
c                 reached over loading history                                  
c               o count number of integration pts at which the                  
c                 current normal tractions is less than 20% of the              
c                 peak value                                                    
c                                                                               
c              For element death, all integration points must have              
c              have a normal traction below 20% of peak value.                  
c                                                                               
      kill_fract = 0.20   !  make input parameter later                         
      kill_criterion_element = .true.                                           
      if( local_debug ) write(out,9100) elem                                    
      hoffset     = (rel_elem-1)*hist_size*ngp                                  
      epsoffset   = (rel_elem-1)*nstr*ngp                                       
      sigoffset   = (rel_elem-1)*nstrs*ngp                                      
c                                                                               
      count_past_peak = 0                                                       
c                                                                               
      do gp = 1, ngp                                                            
c                                                                               
c                  values for possible printing                                 
c                                                                               
        tn = urcs_n(sigoffset+3)                                                
        dn = eps_n(epsoffset+3)                                                 
        tn_max_over_history =  history(hoffset+12)                              
        dn_max_at_tn_max    =  history(hoffset+13)                              
        epsoffset = epsoffset + nstr                                            
        sigoffset = sigoffset + nstrs                                           
        hoffset   = hoffset + hist_size                                         
c                                                                               
c                  does this point satisfy criterion to                         
c                  kill element? If not, element is not killable                
c                  at this load step                                            
c                                                                               
        kill_gp =  dn .gt. dn_max_at_tn_max  .and.                              
     &             tn .le. kill_fract * tn_max_over_history                     
        if( dn .gt. dn_max_at_tn_max )                                          
     &      count_past_peak = count_past_peak + 1                               
        if( .not. kill_gp ) kill_criterion_element = .false.                    
c                                                                               
        if( local_debug .and. gp .eq. 1 )                                       
     &       write(out,9010) gp, tn, dn, tn_max_over_history,                   
     &                       dn_max_at_tn_max                                   
      end do                                                                    
c                                                                               
      values(5) = dble( count_past_peak )                                       
c                                                                               
c              does element need to be killed? needs fixing ....!               
c                                                                               
      kill_now    = .false.                                                     
      if ( compute_type .ne. 2 ) return                                         
      kill_now =  kill_criterion_element                                        
      return                                                                    
c                                                                               
 9100 format('**** debug cavit in dam_param_cohes element: ',i6)                
 9010 format(                                                                   
     &   5x,'gp:                  ',i2,                                         
     & /,5x,'tn:                  ',f12.8,                                      
     & /,5x,'dn:                  ',f12.8,                                      
     & /,5x,'tn_max_over_history: ',f12.8,                                      
     & /,5x,'dn_max_at_tn_max:    ',f12.8 )                                     
c                                                                               
      end subroutine dam_param_cohes_cavit                                      
c                                                                               
      end subroutine dam_param_cohes                                            
c     contains completed                                                        
c     ==================                                                        
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine dam_print_elem4                   *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 11/28/2010 rhd             *          
c     *                                                              *          
c     *     This routine prints out the status of cohesive elements  *          
c     *     marked as killable at the beginning of a load step       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine dam_print_elem4( step, iter )                                  
      use global_data ! old common.main
      use elem_extinct_data, only : dam_state, dam_print_list                   
      use main_data, only : output_packets, packet_file_no                      
      use damage_data                                                           
      implicit integer (a-z)                                                    
c                                                                               
c                    declare local variables                                    
c                                                                               
      double precision                                                          
     &    dummy, eps_plas, eps_crit, sig_mean, sig_mises, eps_plas_tol,         
     &    d_eps_plas, max_d_eps_plas, values(20),half                           
      character(len=10) :: special_char                                         
      logical ldummy, debug, all_killed, beyond_peak,                           
     & found_fgm_cohes, found_homog_cohes, found_ppr,                           
     & found_exponential, option_ppr, option_exponential,                       
     & write_header, show_results, option_cavit, found_cavit                    
      data zero, half, debug, eps_plas_tol / 0.0, 0.5, .false., 1.0e-9 /        
c                                                                               
c          1.0  print the header                                                
c                                                                               
      write(out,*) ' '                                                          
      write(out,*) ' ********************************************* '            
      write(out,*) ' ***       Cohesive element status         *** '            
      write(out,*) ' ********************************************* '            
      write(out,*) ' '                                                          
c                                                                               
c         2.0  check to see if all elements in print list have been killed.     
c              if so, print a message and return                                
c                                                                               
      all_killed = .true.                                                       
      do elem_loop = 1, num_print_list                                          
         element  = dam_print_list(elem_loop)                                   
         elem_ptr = dam_ptr(element)                                            
         if( dam_ptr(element) .eq. 0 ) cycle                                    
         if( dam_state(elem_ptr) .ne. 0 ) cycle                                 
         all_killed = .false.                                                   
      end do                                                                    
c                                                                               
      if( all_killed ) then                                                     
         write(out,9040)                                                        
         return                                                                 
      end if                                                                    
c                                                                               
c         3.0  process each element in the list of elements specifed            
c              by the user. skipped already killed eleemnts.                    
c                                                                               
      local_count       = 0                                                     
      found_homog_cohes = .false.                                               
      found_fgm_cohes   = .false.                                               
      found_ppr         = .false.                                               
      found_exponential = .false.                                               
      found_cavit       = .false.                                               
      write_header      = .true.                                                
c                                                                               
c          loop over all model elements                                         
c                                                                               
      do elem_loop = 1, num_print_list                                          
         element  = dam_print_list(elem_loop)                                   
         elem_ptr = dam_ptr(element)                                            
c                                                                               
c             check if element is a killable and/or if it has                   
c             already been killed.                                              
c                                                                               
         if ( dam_ptr(element) .eq. 0 ) cycle                                   
         if ( dam_state(elem_ptr) .ne. 0 ) cycle                                
c                                                                               
c             cohesive element is still active in the model. call a             
c             routine to compute a vector of values to print for                
c             this element. these are average values computed from              
c             the element gauss point values.                                   
c                                                                               
         call dam_param_cohes( element, ldummy, debug, values, 1 )              
         cohes_type         = iprops(27,element)                                
         option_exponential = cohes_type .eq. 4                                 
         option_ppr         = cohes_type .eq. 6                                 
         option_cavit       = cohes_type .eq. 7                                 
         if( option_exponential ) found_exponential = .true.                    
         if( option_ppr )         found_ppr = .true.                            
         if( option_cavit )       found_cavit = .true.                          
c                                                                               
c              write results for exponential model                              
c                                                                               
         if( option_exponential ) then                                          
           if( write_header ) then                                              
              write(out,9030)                                                   
              write_header = .false.                                            
           end if                                                               
           beyond_peak = values(6) .gt.  values(7)                              
           special_char(1:10) = ' '                                             
           if ( beyond_peak ) special_char(1:3) = '(*)'                         
           if ( values(6) .gt. half*values(7) ) then                            
             write(out,9004) element, (values(k),k=1,6),                        
     &                       special_char(1:3)                                  
             local_count = local_count + 1                                      
             if( abs(props(34,element)) .gt. 0.1 ) then                         
                 found_fgm_cohes = .true.                                       
             else                                                               
                 found_exponential = .true.                                     
             end if                                                             
           end if                                                               
           cycle                                                                
         end if                                                                 
c                                                                               
c              write results for ppr model                                      
c                                                                               
         if( option_ppr ) then                                                  
           if( write_header ) then                                              
              write(out,9032)                                                   
              write_header = .false.                                            
           end if                                                               
           num_gp_post_peak_normal = int(values(7))                             
           num_gp_post_peak_shear  = int(values(8))                             
           num_gp_half_peak_normal = int(values(9))                             
           num_gp_half_peak_shear  = int(values(10))                            
           show_results =                                                       
     &       (num_gp_half_peak_normal + num_gp_half_peak_shear) .gt. 0          
           if( show_results ) then                                              
              write(out,9005) element, (values(k),k=1,6),                       
     &                num_gp_post_peak_normal, num_gp_post_peak_shear           
              local_count = local_count + 1                                     
           end if                                                               
           cycle                                                                
         end if                                                                 
c                                                                               
c              write results for cavit model                                    
c                                                                               
         if( option_cavit ) then                                                
           if( write_header ) then                                              
              write(out,9034)                                                   
              write_header = .false.                                            
           end if                                                               
           write(out,9008) element, (values(k),k=1,4),                          
     &                     int(values(5))                                       
           cycle                                                                
         end if                                                                 
c                                                                               
         write(out,9100)                                                        
         call die_gracefully                                                    
      end do                                                                    
c                                                                               
c         4.0  write nice ending explanations for table of status               
c              values just printed                                              
c                                                                               
      write(out,*) ' '                                                          
      if( found_exponential ) then                                              
           if( local_count .gt. 0 ) write(out,9050)                             
           if( local_count .eq. 0 ) write(out,9052)                             
      end if                                                                    
      if( found_ppr ) then                                                      
         if( local_count .gt. 0 ) write(out,9055)                               
         if( local_count .eq. 0 ) write(out,9058)                               
      end if                                                                    
      if( found_cavit ) then                                                    
         write(out,9060)                                                        
      end if                                                                    
      write(out,*) ' '                                                          
c                                                                               
c         5.0 for packet output, we now know the number of data lines           
c             (elements) to be output. run thru the above loop again to         
c             write packets.                                                    
c                                                                               
      if ( local_count .eq. 0 ) go to 9000                                      
      if ( .not. output_packets ) go to 9000                                    
      write(packet_file_no) 7, local_count, step, iter                          
      do elem_loop = 1, num_print_list                                          
         element  = dam_print_list(elem_loop)                                   
         elem_ptr = dam_ptr(element)                                            
c                                                                               
c             check if element is a killable element and/or if it has           
c             already been killed.                                              
c                                                                               
         if ( dam_ptr(element) .eq. 0 ) cycle                                   
         if ( dam_state(elem_ptr) .ne. 0 ) cycle                                
c                                                                               
c             cohesive element is still active in the model. call a             
c             routine to compute a vector of values to issue values into        
c             binary file                                                       
c                                                                               
         call dam_param_cohes( element, ldummy, debug, values, 1 )              
         cohes_type  = iprops(27,element)                                       
         option_exponential = cohes_type .eq. 4                                 
         option_ppr         = cohes_type .eq. 6                                 
         if( option_exponential ) then                                          
          beyond_peak = values(6) .gt. values(7)                                
          special_char(1:3) = ' '                                               
          if ( beyond_peak ) special_char(1:3) = '(*)'                          
          if ( values(6) .gt. 0.5*values(7) )                                   
     &     write(packet_file_no) element, cohes_type, (values(k),k=1,6),        
     &                           special_char                                   
         end if                                                                 
         if( option_ppr ) then                                                  
           num_gp_post_peak_normal = int(values(7))                             
           num_gp_post_peak_shear  = int(values(8))                             
           num_gp_half_peak_normal = int(values(9))                             
           num_gp_half_peak_shear  = int(values(10))                            
           show_results =                                                       
     &       (num_gp_half_peak_normal + num_gp_half_peak_shear) .gt. 0          
           if( show_results )                                                   
     &        write(packet_file_no) element, cohes_type,                        
     &           (values(k),k=1,6),                                             
     &           num_gp_post_peak_normal, num_gp_post_peak_shear                
         end if                                                                 
      end do                                                                    
c                                                                               
 9000 continue                                                                  
      return                                                                    
c                                                                               
 9040 format(1x,'*** NOTE: all elements in the killable element',               
     & ' status list have been killed.',/)                                      
 9030 format(1x,'element       Tn             Ts            Dn    ',            
     &       '        Ds           Teff       ',                                
     &       '  Deff    ',/,1x,                                                 
     &       '-------   --------       ---------      --------',                
     &       '     ---------    ----------- ',                                  
     &       '  -----------')                                                   
 9032 format(1x,'element    T-normal       T-shear        D-normal ',           
     &       '    D-shear     Tn/Tn-peak   ',                                   
     &       ' Ts/Ts-peak   ',/,1x,                                             
     &       '-------   --------       ---------      --------',                
     &       '     ---------    ----------- ',                                  
     &       '  -----------')                                                   
 9034 format(1x,'element    T-normal       T-shear        D-normal ',           
     &       '    D-shear     (1)',/,1x,                                        
     &       '-------   --------       ---------      --------',                
     &       '     ---------    ---- ' )                                        
 9050 format(                                                                   
     &   1x,'Exponential option: **** only elements with Deff > 0.5 x ',        
     &      'Dpeak listed ****',                                                
     & /,1x,'(*) denotes deformation beyond peak stress',                       
     & /,1x,'**** T and D values denote independent averages',                  
     &      ' computed over element integration points ****',                   
     & /,6x,'(T & D values may thus not exactly satisfy cohesive',              
     &      ' constitutive model)' )                                            
 9052 format(                                                                   
     &   1x,'Exponential option: **** no elements with Deff > 0.5 x ',          
     &      'Dpeak ****' )                                                      
 9055 format(                                                                   
     &    1x,'PPR option: ',                                                    
     &       'T -> traction values, D -> relative displacement values',         
     &  /,1x,'**** only elements with at least 50% peak loading ',              
     &      'listed ****',                                                      
     &  /,1x,'(#,#) denotes number of Gauss points for (normal,shear)',         
     &       ' displacement that have exceeded the coresponding',               
     &  /,1x,'      values at peak (normal,shear) stress over the',             
     &       ' loading history',                                                
     &  /,1x,'**** T and D values shown denote independent averages',           
     &       ' computed over element integration points ****',                  
     &  /,6x,                                                                   
     & '(these T & D values may thus not exactly satisfy cohesive',             
     &       ' constitutive model)' )                                           
 9058 format(                                                                   
     &    1x,'PPR option: ',                                                    
     &       '**** no elements with at least 50% peak loading **** ')           
 9060 format(                                                                   
     &   1x,'Cavit option: integration point average tractions & ',             
     &      'displacement jumps',                                               
     &  /,1x,'(1) Number of integration points with post-peak state ',          
     &         'based on T-normal, D-normal')                                   
 9004 format(1x,i7,1x,e14.6,1x,e14.6,                                           
     &       e14.6, e14.6, e14.6, e14.6, 1x,a3)                                 
 9005 format(1x,i7,1x,e14.6,1x,e14.6,                                           
     &       e14.6, e14.6, f11.5, f13.5, 3x,'(',i1,',',i1,')')                  
 9008 format(1x,i7,1x,e14.6,1x,e14.6,                                           
     &       e14.6, e14.6, 3x, I1)                                              
 9100 format(//,'>>> FATAL ERROR: dam_print_elem4'                              
     &        /,'                 job terminated' )                             
c                                                                               
      end                                                                       
