c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine incurv                       *          
c     *                                                              *          
c     *                       written by : kck                       *          
c     *                                                              *          
c     *                   last modified : 7/29/2011 rhd              *          
c     *                                                              *          
c     *     this subroutine supervises and conducts the input of     *          
c     *     segmentally defined stress-strain curves including       *          
c     *     temperature dependent and strain-rate dependent cases    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine incurv( sbflg1, sbflg2 )                                       
      use segmental_curves                                                      
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
c                                                                               
c                local declarations                                             
c                                                                               
      logical matchs, numi, numd, sbflg1, sbflg2, new_line, local_debug,        
     &        endcrd, local_gp_model_flag, found_points                         
      double precision                                                          
     &   stress_pt, strain_pt, dumd, ym, zero                                   
      real dumr                                                                 
      character(len=1) :: dums                                                  
      data zero, local_debug / 0.0, .false. /                                   
c                                                                               
      if( sbflg1 ) then                                                         
         call errmsg(222,dum,dums,dumr,dumd)                                    
         go to 30                                                               
      end if                                                                    
c                                                                               
c             called by driver. num_points, max_current_pts,                    
c             max_current_curves, num_points,                                   
c             num_seg_curve_sets initialized in initst and then                 
c             after each curve definition by errck_23 called from driver        
c             immediately after this routine                                    
c                                                                               
c             continue reading the start up line                                
c                                                                               
      if( matchs('-',1) )      call splunj                                      
      if( matchs('strain',6) ) call splunj                                      
      if( matchs('curve',5) )  call splunj                                      
      if( .not. numi(num_curve) ) then                                          
         call errmsg(217,dum,'curve',dumr,dumd)                                 
         num_curve = 1                                                          
      end if                                                                    
      if( num_curve .le. 0 .or. num_curve .gt. max_seg_curves ) then            
         call errmsg(315,max_seg_curves,dums,dumr,dumd)                         
         call die_gracefully                                                    
         stop                                                                   
      end if                                                                    
c                                                                               
c            get temperature or strain rate value for this curve if             
c            specified. mark the 'type" of data provided for the                
c            curve. =0 (temp. & rate insensitive), =1 (temp. dependent)         
c            =2 (rate dependent). for temperature dependent, the user           
c            must also specify young's modulus, poisson's ratio and             
c            thermal expansion coefficient. for strain-rate dependent,          
c            the modulus, nu and alpha are all constant and specified           
c            in the material definition.                                        
c                                                                               
      local_gp_model_flag = .false.                                             
      seg_curves_value(num_curve) = zero                                        
      seg_curves_type(num_curve)  = 0                                           
      if ( endcrd(dumr) ) go to 30                                              
c                                                                               
      if ( matchs('temperature',4) ) then                                       
        call incurv_temperature( local_gp_model_flag )                          
        go to 30                                                                
      end if                                                                    
c                                                                               
      if ( matchs('plastic',4) ) then                                           
        call incurv_strain_rate                                                 
        go to 30                                                                
      end if                                                                    
c                                                                               
      call errmsg2( 3, dumr, dums, dumr, dumd )                                 
c                                                                               
c                                                                               
c             read in the strain,stress values to define the curve              
c                                                                               
c             going to line 10 means read a new line, go to 20                  
c             means look for another pair on the same line.                     
c                                                                               
c             the temperature dependent cyclic model ("gp" here) has only       
c             properties defined. no curve is defined. to make a lot            
c             error checking, save, reopen, etc. code work, we define           
c             two phoney points on a curve.                                     
c                                                                               
 30   continue                                                                  
      if ( local_debug ) then                                                   
        write(*,*) '>> curve type: ',seg_curves_type(num_curve)                 
        write(*,*) '>> curve type value: ',seg_curves_value(num_curve)          
      end if                                                                    
c                                                                               
      found_points = .false.                                                    
 10   continue                                                                  
      call readsc                                                               
      new_line = .true.                                                         
 20   continue                                                                  
      if( matchs(',',1) ) call splunj                                           
      if( numd(strain_pt) ) then                                                
         if( numd(stress_pt) ) then                                             
            num_points = num_points + 1                                         
            if( num_points .gt. max_seg_points ) then                           
               num_points = max_seg_points                                      
               call errmsg(220,max_seg_points,dums,dumr,dumd)                   
               go to 10                                                         
            endif                                                               
            seg_curves(num_points,1,num_curve) = strain_pt                      
            seg_curves(num_points,2,num_curve) = stress_pt                      
          new_line = .false.                                                    
          go to 20                                                              
         else                                                                   
            call errmsg(219,num_points,dums,dumr,dumd)                          
          go to 10                                                              
         end if                                                                 
      else                                                                      
         if( new_line ) then                                                    
          go to 9999                                                            
         else                                                                   
          go to 10                                                              
         end if                                                                 
      end if                                                                    
      go to 20                                                                  
c                                                                               
 9999 continue                                                                  
c                                                                               
c         if points are given for the gp_model, this is an error.               
c         if correct input for gp_model is given (no points),                   
c         create two dummy points so all subsequent data checks                 
c         will pass. The values are never used in real processing               
c         routines for the gp option of the cyclic model                        
c                                                                               
      if( local_gp_model_flag .and. (num_points .ne. 0) )                       
     &       call errmsg2( 80, dumr, dums, dumr, dumd )                         
      if( local_gp_model_flag ) then                                            
            num_points = 1                                                      
            seg_curves(num_points,1,num_curve) = 0.0                            
            seg_curves(num_points,2,num_curve) = 1.0                            
            num_points = 2                                                      
            seg_curves(num_points,1,num_curve) = 0.5                            
            seg_curves(num_points,2,num_curve) = 2.0                            
      end if                                                                    
c                                                                               
c         for type=0 curves, the first point cannot have zero strain            
c         (if it does, we cannot compute the modulus!)                          
c                                                                               
      if (  seg_curves_type(num_curve) .eq. 0 .and.                             
     &      seg_curves(1,1,num_curve) .eq. zero ) then                          
        call errmsg(233,num_points,dums,dumr,dumd)                              
        do i = 1, num_points                                                    
          seg_curves(i,1,num_curve) = zero                                      
          seg_curves(i,2,num_curve) = zero                                      
        end do                                                                  
        return                                                                  
      end if                                                                    
c                                                                               
c         for type=1,2 curves, the first point must have a zero                 
c         plastic strain so that we know the yield stress                       
c                                                                               
      if (  seg_curves_type(num_curve) .ne. 0 .and.                             
     &      seg_curves(1,1,num_curve) .ne. zero ) then                          
        call errmsg2(21,num_points,dums,dumr,dumd)                              
        do i = 1, num_points                                                    
          seg_curves(i,1,num_curve) = zero                                      
          seg_curves(i,2,num_curve) = zero                                      
        end do                                                                  
        return                                                                  
      end if                                                                    
c                                                                               
c         the strain values must be monotonically increasing                    
c                                                                               
      do i = 2, num_points                                                      
           if ( seg_curves(i,1,num_curve) .le.                                  
     &          seg_curves(i-1,1,num_curve) ) then                              
             call errmsg2( 15, num_points, dums, dumr, dumd )                   
             do j = 1, num_points                                               
                seg_curves(i,1,num_curve) = zero                                
                seg_curves(i,2,num_curve) = zero                                
             end do                                                             
             return                                                             
           end if                                                               
       end do                                                                   
c                                                                               
c         for type = 0 stress-strain curves (temperature and rate               
c         independent), convert the strain values to plastic strain             
c         by subtracting off the elastic contribution. the curve                
c         input to the material model is thus stress vs. plastic strain.        
c         for temperature or rate dependent curves, the strain                  
c         input values must be plastic strains.                                 
c                                                                               
       if ( seg_curves_type(num_curve) .eq. 0 ) then                            
         ym = seg_curves(1,2,num_curve) / seg_curves(1,1,num_curve)             
         do i = 1, num_points                                                   
           seg_curves(i,1,num_curve) = seg_curves(i,1,num_curve) -              
     &                                 seg_curves(i,2,num_curve) / ym           
         end do                                                                 
       end if                                                                   
c                                                                               
c         locate and save the minimum stress value on the curve.                
c                                                                               
       seg_curves_min_stress(num_curve) = 1.0e20                                
       do i = 1, num_points                                                     
           seg_curves_min_stress(num_curve) =                                   
     &       min( seg_curves_min_stress(num_curve),                             
     &          seg_curves(i,2,num_curve) )                                     
       end do                                                                   
                                                                                
c                                                                               
       if ( local_debug ) then                                                  
         write(*,*) '>> debug for incurv'                                       
         write(*,*) ' curve number: ',num_curve                                 
         write(*,*) ' plastic strain, stress values'                            
         do i = 1, num_points                                                   
           write(*,9000) i, seg_curves(i,1,num_curve),                          
     &                      seg_curves(i,2,num_curve)                           
         end do                                                                 
         write(*,9100) seg_curves_min_stress(num_curve)                         
       end if                                                                   
c                                                                               
       sbflg1 = .true.                                                          
       sbflg2 = .true.                                                          
       return                                                                   
c                                                                               
 9000  format(3x,i3,f14.6,f15.4)                                                
 9100  format(3x,'min stress value: ',f15.4)                                    
       end                                                                      
c     ****************************************************************          
c     *                                                              *          
c     *                subroutine incurv_temperature                 *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 7/28/2011 rhd              *          
c     *                                                              *          
c     *     this subroutine supervises and conducts the input of the *          
c     *     temperature dependent modulus, poisson ratio and         *          
c     *     expansion coefficient for a segmental stress-strain      *          
c     *     curve                                                    *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine incurv_temperature( gp_model_flag )                            
      use segmental_curves                                                      
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
c                                                                               
c                paramter declarations                                          
c                                                                               
      logical gp_model_flag                                                     
c                                                                               
c                local declarations                                             
c                                                                               
      logical matchs, matchs_exact, numd,  new_line, local_debug,               
     &        endcrd, ym_input, nu_input, alpha_input, ok,                      
     &        gp_sigma_0_input, gp_h_u_input,                                   
     &        gp_beta_u_input, gp_delta_u_input,                                
     &        local_match_flag                                                  
c                                                                               
      double precision                                                          
     &    dumd, value                                                           
      real dumr                                                                 
      character(len=1) :: dums                                                  
      data local_debug / .false. /                                              
c                                                                               
      gp_model_flag = .false.                                                   
c                                                                               
c                get the temperature for the curve. if not given, curve         
c                is stored as temperature independent                           
c                                                                               
      seg_curves_type(num_curve)  = 1                                           
      if ( .not. numd(seg_curves_value(num_curve)) ) then                       
         call errmsg2( 1, dumr, dums, dumr, dumd )                              
         seg_curves_type(num_curve)  = 0                                        
         return                                                                 
      end if                                                                    
c                                                                               
c                get young's modulus, poisson's ratio and thermal               
c                expansion coefficient for this temperature.                    
c                values must be given else error.                               
c                                                                               
      ym_input    = .false.                                                     
      nu_input    = .false.                                                     
      alpha_input = .false.                                                     
c                                                                               
c                generalized plasticity terms                                   
c                                                                               
      gp_sigma_0_input = .false.                                                
      gp_h_u_input     = .false.                                                
      gp_beta_u_input  = .false.                                                
      gp_delta_u_input = .false.                                                
c                                                                               
      seg_curves_ym(num_curve)         = -99.0                                  
      seg_curves_nu(num_curve)         = -99.0                                  
      seg_curves_alpha(num_curve)      = -99.0                                  
      seg_curves_gp_sigma_0(num_curve) = -99.0                                  
      seg_curves_gp_h_u(num_curve)     = -99.0                                  
      seg_curves_gp_beta_u(num_curve)  = -99.0                                  
      seg_curves_gp_delta_u(num_curve) = -99.0                                  
c                                                                               
c                loop till end of line to extract materia property              
c                values at this temperature                                     
c                                                                               
      do while ( .not. endcrd(dumr) )                                           
c                                                                               
c                e, nu and alpha                                                
c                                                                               
        if ( matchs_exact('e') ) then                                           
           if ( numd(value) ) then                                              
             ym_input = .true.                                                  
             seg_curves_ym(num_curve) = value                                   
             cycle                                                              
           else                                                                 
             call errmsg2( 16, dumr, dums, dumr, dumd )                         
           end if                                                               
        end if                                                                  
c                                                                               
        if ( matchs_exact('nu') ) then                                          
           if ( numd(value) ) then                                              
             nu_input = .true.                                                  
             seg_curves_nu(num_curve) = value                                   
             cycle                                                              
           else                                                                 
             call errmsg2( 16, dumr, dums, dumr, dumd )                         
           end if                                                               
        end if                                                                  
c                                                                               
        if ( matchs_exact('alpha') ) then                                       
           if ( numd(value) ) then                                              
             alpha_input = .true.                                               
             seg_curves_alpha(num_curve) = value                                
             cycle                                                              
           else                                                                 
             call errmsg2( 16, dumr, dums, dumr, dumd )                         
           end if                                                               
        end if                                                                  
c                                                                               
c                generalized plasticity terms                                   
c                                                                               
        local_match_flag = .false.                                              
        if( matchs_exact('gp_sigma_0') )  local_match_flag = .true.             
        if( matchs_exact('gp_sig_0') )    local_match_flag = .true.             
        if( matchs_exact('gp_yld_pt') )   local_match_flag = .true.             
        if( local_match_flag ) then                                             
           if ( numd(value) ) then                                              
             gp_sigma_0_input = .true.                                          
             seg_curves_gp_sigma_0(num_curve) = value                           
             gp_model_flag = .true.                                             
             cycle                                                              
           else                                                                 
             call errmsg2( 16, dumr, dums, dumr, dumd )                         
           end if                                                               
        end if                                                                  
c                                                                               
        if ( matchs_exact('gp_h_u') ) then                                      
           if ( numd(value) ) then                                              
             gp_h_u_input = .true.                                              
             seg_curves_gp_h_u(num_curve) = value                               
             gp_model_flag = .true.                                             
             cycle                                                              
           else                                                                 
             call errmsg2( 16, dumr, dums, dumr, dumd )                         
           end if                                                               
        end if                                                                  
c                                                                               
        if ( matchs_exact('gp_beta_u') ) then                                   
           if ( numd(value) ) then                                              
             gp_beta_u_input = .true.                                           
             seg_curves_gp_beta_u(num_curve) = value                            
             gp_model_flag = .true.                                             
             cycle                                                              
           else                                                                 
             call errmsg2( 16, dumr, dums, dumr, dumd )                         
           end if                                                               
        end if                                                                  
c                                                                               
        if ( matchs_exact('gp_delta_u') ) then                                  
           if ( numd(value) ) then                                              
             gp_delta_u_input = .true.                                          
             seg_curves_gp_delta_u(num_curve) = value                           
             gp_model_flag = .true.                                             
             cycle                                                              
           else                                                                 
             call errmsg2( 16, dumr, dums, dumr, dumd )                         
           end if                                                               
        end if                                                                  
c                                                                               
c                property value not recognized                                  
c                                                                               
        call errmsg2( 17, dumr, dums, dumr, dumd )                              
        seg_curves_type(num_curve)  = 0                                         
        return                                                                  
c                                                                               
      end do                                                                    
c                                                                               
c                                                                               
c                end of data line. make sure all required data                  
c                has been entered. only e, nu, and alpha must be                
c                given at this point for all models.                            
c                                                                               
c                                                                               
      ok = ym_input .and. nu_input .and. alpha_input                            
      if ( .not. ok ) then                                                      
        call errmsg2( 18, dumr, dums, dumr, dumd )                              
        seg_curves_type(num_curve)  = 0                                         
        return                                                                  
      end if                                                                    
c                                                                               
c                all 4 properties of the gp model must be given if              
c                any are given.                                                 
c                                                                               
      if(  gp_model_flag ) then                                                 
        ok = gp_sigma_0_input .and. gp_h_u_input .and. gp_beta_u_input          
     &       .and. gp_delta_u_input                                             
        if( .not. ok ) then                                                     
          call errmsg2( 79, dumr, dums, dumr, dumd )                            
          seg_curves_type(num_curve)  = 0                                       
        end if                                                                  
      end if                                                                    
                                                                                
c                                                                               
      if ( local_debug ) then                                                   
        write(*,9000) seg_curves_value(num_curve),                              
     &             seg_curves_ym(num_curve), seg_curves_nu(num_curve),          
     &             seg_curves_alpha(num_curve),                                 
     &             seg_curves_gp_sigma_0(num_curve),                            
     &             seg_curves_gp_h_u(num_curve),                                
     &             seg_curves_gp_beta_u(num_curve),                             
     &             seg_curves_gp_delta_u(num_curve)                             
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format('>> temperature:',f10.3,' ym: ',f10.1,' nu: ',f7.2,                
     &     ' alpha:', f12.8,/,                                                  
     &     ' gp_sigma_0: ',f7.3,' gp_h_u: ',f10.2,/,                            
     &     ' gp_beta_u: ',f10.3,' gp_delta_u: ',f10.4 )                         
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                subroutine incurv_strain_rate                 *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 4/15/00 rhd                *          
c     *                                                              *          
c     *     this subroutine supervises and conducts the input of the *          
c     *     strain-rate dependent plastic-strain vs. stress,         *          
c     *     segmental curve                                          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine incurv_strain_rate                                             
      use segmental_curves                                                      
      implicit integer (a-z)                                                    
      include 'param_def'                                                       
c                                                                               
c                local declarations                                             
c                                                                               
      logical matchs, numd, local_debug                                         
c                                                                               
      double precision                                                          
     &    dumd                                                                  
      real dumr                                                                 
      character(len=1) :: dums                                                  
      data local_debug / .false. /                                              
c                                                                               
c                get the plastic strain rate for the curve.                     
c                if not given, curve is stored as rate independent              
c                                                                               
      seg_curves_type(num_curve)  = 2                                           
      if ( matchs('strain',4) ) call splunj                                     
      if ( matchs('-',1) )      call splunj                                     
      if ( matchs('rate',4) )   call splunj                                     
      if ( .not. numd(seg_curves_value(num_curve)) ) then                       
         call errmsg2( 2, dumr, dums, dumr, dumd )                              
         seg_curves_type(num_curve)  = 0                                        
         return                                                                 
      end if                                                                    
c                                                                               
      if ( local_debug ) then                                                   
        write(*,9000) seg_curves_value(num_curve)                               
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format('>> strain-rate:',f10.3)                                           
c                                                                               
      end                                                                       
                                                                                
                                                                                
                                                                                
