c *******************************************************************           
c *                                                                 *           
c *        material model # 9 -- <available>                        *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine mm09(                                                          
     &  step, iter, felem, gpn, mxvl, hist_size, nstrs, nstrn, span,            
     &  iout, signal_flag, adaptive_possible, cut_step_size_now,                
     &  mm_props, e_vec, tan_e_vec, nu_vec, sigyld_vec,                         
     &  n_power_vec, trial_elas_stress_np1, stress_n, stress_np1,               
     &  deps, history_n, history_np1, killed_list_vec )                         
      implicit none                                                             
c                                                                               
c                   parameter declarations                                      
c                   ----------------------                                      
c                                                                               
      integer                                                                   
     &  step, iter, felem, gpn, mxvl, hist_size, span,                          
     &  iout, nstrs, nstrn                                                      
c                                                                               
      logical                                                                   
     &   signal_flag, adaptive_possible, cut_step_size_now,                     
     &   killed_list_vec(mxvl)                                                  
c                                                                               
c                                                                               
      double precision                                                          
     & mm_props(mxvl,10), e_vec(mxvl), tan_e_vec(mxvl), nu_vec(mxvl),           
     & sigyld_vec(mxvl), n_power_vec(mxvl), stress_n(mxvl,nstrs),               
     & stress_np1(mxvl,nstrs), deps(mxvl,nstrn),                                
     & trial_elas_stress_np1(mxvl,nstrn), history_n(span,hist_size),            
     & history_np1(span,hist_size)                                              
c                                                                               
c               description of parameters                                       
c               -------------------------                                       
c                                                                               
c     step              : current load step number                              
c     iter              : current newton iteration number                       
c     felem             : first element of the current block                    
c     gpn               : gauss point number being processed for block          
c     mxvl              : maximum no. elements per block                        
c     hist_size         : number of history words per gauss point               
c     nstrs             : number of stress terms (6 + extras)                   
c     nstrn             : number of incremental strain terms (6)                
c     span              : number of elements in current block                   
c     iout              : write messates to this device number                  
c     signal_flag       : user wants notification messages for key              
c                         events in material response                           
c     adaptive_possible : .true. if the material model may request              
c                         immediate reduction of global load step size.         
c                         no stress-histroy update required                     
c (*) cut_step_size_now : set .true. if material model wants immediate          
c                         reduction of global load step size.                   
c                         no stress-histroy update required                     
c     mm_props          : material parameter values input by user for           
c                         each element in block                                 
c     e_vec             : Young's modulus for each element in block             
c     tan_e_vec         : constant tangent modulus after yield (may not         
c                         be defined) for each element in block                 
c     nu_vec            : Poisson's ratio for each element in block             
c     sigyld_vec        : yield stress for each element in block                
c     n_power_vec       : power-law hardening exponent for each element         
c                         in block (may not be defined)                         
c (*) trial_elas_stress_np1 : trial elastic stress vector to be used later by   
c                         consistent tangent routine for model                  
c (**)stress_n          : stresses at start of load step (n) for all            
c                         elements in block for this gauss point                
c (*) stress_np1        : stresses at end of load step (n+1) for all            
c                         elements in block for this gauss point                
c     deps              : current estimate of strain increment over the         
c                         load step (minus increment of thermal strain)         
c (**)history_n         : history values at start of load step (n) for all      
c                         elements in block for this gauss point                
c (*) history_np1       : history values at end of load step (n+1) for all      
c                         elements in block for this gauss point                
c     killed_list_vec   : .true. if this element has been previously killed     
c                         by crack growth operations                            
c                                                                               
c    (*)  values to be updated by this material model                           
c    (**) needs to be initialized on step 1                                     
c                                                                               
c                                                                               
c                                                                               
c   strain ordering:                                                            
c     deps-xx, deps-yy, deps-zz, gamma-xy, gamma-yz, gamma-xz                   
c                                                                               
c   stress ordering (at n and n+1):                                             
c     (1) sig-xx                                                                
c     (2) sig-yy                                                                
c     (3) sig-zz                                                                
c     (4) tau-xy                                                                
c     (5) tau-yz                                                                
c     (6) tau-xz                                                                
c     (7) total work density                                                    
c     (8) total plastic work density                                            
c     (9) total plastic strain                                                  
c                                                                               
c                                                                               
c                                                                               
c                   local variables                                             
c                   ---------------                                             
c                                                                               
      integer i                                                                 
      logical local_debug                                                       
      double precision                                                          
     &     shear_mod, c, a, zero, one, two, b, delastic,                        
     &     half                                                                 
      data zero, one, two, half / 0.0d00, 1.0d00, 2.0d00, 0.5d00 /              
      data local_debug / .false. /                                              
c                                                                               
c                                                                               
      if ( local_debug ) then                                                   
         write(iout,*) '... inside mm09 ...'                                    
         write(iout,9000) step, iter, felem, gpn, mxvl, hist_size,              
     &                    nstrs, nstrn, span, signal_flag,                      
     &                    adaptive_possible, cut_step_size_now                  
         do i = 1, span                                                         
          write(iout,9010)  felem+i-1, deps(i,1:6)                              
         end do                                                                 
      end if                                                                    
c                                                                               
c                initialize stresses on step 1                                  
c                                                                               
      if ( step .eq. 1 ) then                                                   
        do i = 1, span                                                          
          stress_n(i,1) = zero                                                  
          stress_n(i,2) = zero                                                  
          stress_n(i,3) = zero                                                  
          stress_n(i,4) = zero                                                  
          stress_n(i,5) = zero                                                  
          stress_n(i,6) = zero                                                  
          stress_n(i,7) = zero                                                  
          stress_n(i,8) = zero                                                  
          stress_n(i,9) = zero                                                  
        end do                                                                  
      end if                                                                    
c                                                                               
c                update stresses as just isotropic, linear-elastic              
c                                                                               
      do i = 1, span                                                            
c                                                                               
        c = e_vec(i) / ( ( one + nu_vec(i) ) *                                  
     &      ( one - two * nu_vec(i) ) )                                         
        a = c * ( one - nu_vec(i) )                                             
        b = nu_vec(i) * c                                                       
        shear_mod = e_vec(i) / ( two*(one+nu_vec(i)))                           
c                                                                               
        stress_np1(i,1) = stress_n(i,1) + a * deps(i,1) +                       
     &                    b * deps(i,2) +  b * deps(i,3)                        
        stress_np1(i,2) = stress_n(i,2) + a * deps(i,2) +                       
     &                    b * deps(i,1) +  b * deps(i,3)                        
        stress_np1(i,3) = stress_n(i,3) + a * deps(i,3) +                       
     &                    b * deps(i,1) +  b * deps(i,2)                        
        stress_np1(i,4) = stress_n(i,4) + shear_mod * deps(i,4)                 
        stress_np1(i,5) = stress_n(i,5) + shear_mod * deps(i,5)                 
        stress_np1(i,6) = stress_n(i,6) + shear_mod * deps(i,6)                 
c                                                                               
      end do                                                                    
c                                                                               
c                update energy densities                                        
c                                                                               
      do i = 1, span                                                            
         delastic =                                                             
     &  half * ( stress_np1(i,1) + stress_n(i,1) ) *  deps(i,1) +               
     &  half * ( stress_np1(i,2) + stress_n(i,2) ) *  deps(i,2) +               
     &  half * ( stress_np1(i,3) + stress_n(i,3) ) *  deps(i,3) +               
     &  half * ( stress_np1(i,4) + stress_n(i,4) ) *  deps(i,4) +               
     &  half * ( stress_np1(i,5) + stress_n(i,5) ) *  deps(i,5) +               
     &  half * ( stress_np1(i,6) + stress_n(i,6) ) *  deps(i,6)                 
        stress_np1(i,7) = stress_n(i,7) +  delastic                             
        stress_np1(i,8) = zero                                                  
        stress_np1(i,9) = zero                                                  
      end do                                                                    
                                                                                
                                                                                
      if ( local_debug ) then                                                   
         write(iout,*) '... updated stresses'                                   
         do i = 1, span                                                         
          write(iout,9020) felem+i-1, stress_np1(i,1:7)                         
         end do                                                                 
      end if                                                                    
                                                                                
      return                                                                    
c                                                                               
 9000 format(/,10x,'step, iter, felem, gpn:         ',4i8,                      
     &       /,10x,'mxvl, hist_size, nstrs:         ',3i8,                      
     &       /,10x,'nstrn, span:                    ',2i8,                      
     &       /,10x,'signal_flag, adaptive_possible: ',2l3,                      
     &       /,10x,'cut_step_size_now:              ',l3,/)                     
 9010 format(5x,i5,6f10.6)                                                      
 9020 format(5x,i5,7f10.3)                                                      
                                                                                
                                                                                
       end                                                                      
                                                                                
                                                                                
c *******************************************************************           
c *                                                                 *           
c *        material model # 9 -- <available>                        *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine cnst9(                                                         
     &  span, felem, gpn, first, iter, iout, mxvl, nstrn,                       
     &  weight, e_vec, nu_vec, sigyld_vec, n_power_vec, mm_props,               
     &  trial_elas_stress, history_n,                                           
     &  history_np1, stress_np1, dmat, det_jac_block,                           
     &  killed_list_vec )                                                       
      implicit none                                                             
c                                                                               
c                   parameter declarations                                      
c                   ----------------------                                      
c                                                                               
      integer                                                                   
     &  span, felem, gpn, iter, iout, mxvl, nstrn                               
c                                                                               
      logical                                                                   
     &   first, killed_list_vec(mxvl)                                           
c                                                                               
      double precision                                                          
     & weight, mm_props(mxvl,5), e_vec(mxvl), nu_vec(mxvl),                     
     & trial_elas_stress(mxvl,nstrn), history_n(span,*),                        
     & history_np1(span,*), dmat(mxvl,nstrn,nstrn),                             
     & det_jac_block(mxvl), stress_np1(mxvl,nstrn),                             
     & sigyld_vec(mxvl), n_power_vec(mxvl)                                      
c                                                                               
c                                                                               
c                                                                               
c               description of parameters                                       
c               -------------------------                                       
c                                                                               
c     step              : current load step number                              
c     iter              : current newton iteration number. iter 1 is            
c                         for application of the "real" load increment.         
c     felem             : first element of the current block                    
c     gpn               : gauss point number being processed for block          
c     first             : logical flag indicating if this                       
c                         is the very first call to the routine                 
c                         the linear elastic [d] is returned                    
c                         whenever first = .true.                               
c     mxvl              : maximum no. elements per block                        
c     nstrn             : number of strain-stress components (=6)               
c     span              : number of elements in current block                   
c     iout              : write messages to this device number                  
c     mm_props          : material parameter values input by user for           
c                         each element in block                                 
c     e_vec             : Young's modulus for each element in block             
c     nu_vec            : Poisson's ratio for each element in block             
c     sigyld_vec        : yield stress for each element in block                
c     n_power_vec       : power-law hardening exponent                          
c (#) trial_elas_stress_np1 : trial elastic stress vector defined by stress     
c                             update routine                                    
c     history_n         : history values at start of load step (n) for all      
c                         elements in block for this gauss point                
c     history_np1       : history values at end of load step (n+1) for all      
c                         elements in block for this gauss point                
c     det_jac_block     : |J| at this gauss point for each element in           
c                         block                                                 
c (!) stress_np1        : current estimate of 6 stress components for           
c                         end of step (see ordering below)                      
c     weight            : integration point weight factor                       
c (*) dmat              : 6x6 (symmetric) tangent (consistent) for              
c                         this gauss point for each element of block            
c                         (see stress ordering below)                           
c                                                                               
c    (*)  values to be updated by this material model                           
c    (!)  for finite strain computations, these are unrotated                   
c         Cauchy stresses                                                       
c    (#)  used by constitutive update procedures based on some                  
c         for of elastic predictor - return mapping algorithm.                  
c         the contents of this array are set by the corresponding               
c         stress update routine for the model. for finite                       
c         strain computations, the contents will be unrotated                   
c         Cauchy stress terms of some form as set by the stress                 
c         update routine.                                                       
c                                                                               
c   Finite strain issues:                                                       
c     this constitutive model only sees stresses that have already              
c     been rotation "neutralized" by WARP3D. this routine can                   
c     thus operate simply as "small strain" theory. WARP3D will                 
c     "rotate" the compute [D] matrices on return for finite strain-            
c     rotation effects.                                                         
c                                                                               
c   strain ordering:                                                            
c     deps-xx, deps-yy, deps-zz, gamma-xy, gamma-yz, gamma-xz                   
c                                                                               
c   stress ordering (at n and n+1):                                             
c     sig-xx, sig-yy, sig-zz, tau-xy, tau-yz, tau-xz                            
c                                                                               
c                                                                               
c   Note: warp3d expects all dmat[] values to be multiplied by                  
c         weight * det_jac_block(i), where i = relative                         
c         element number of block                                               
c                                                                               
c                                                                               
c                   local variables                                             
c                   ---------------                                             
c                                                                               
      integer i                                                                 
      logical debug                                                             
      double precision                                                          
     &     c1, c2, c3, c4, fact, zero, one, two                                 
      data zero,one,two / 0.0d00, 1.0d00, 2.0d00 /                              
c                                                                               
                                                                                
      do i = 1, span                                                            
         dmat(i,1,4) = zero                                                     
         dmat(i,1,5) = zero                                                     
         dmat(i,1,6) = zero                                                     
         dmat(i,2,4) = zero                                                     
         dmat(i,2,5) = zero                                                     
         dmat(i,2,6) = zero                                                     
         dmat(i,3,4) = zero                                                     
         dmat(i,3,5) = zero                                                     
         dmat(i,3,6) = zero                                                     
         dmat(i,4,1) = zero                                                     
         dmat(i,4,2) = zero                                                     
         dmat(i,4,3) = zero                                                     
         dmat(i,4,5) = zero                                                     
         dmat(i,4,6) = zero                                                     
         dmat(i,5,1) = zero                                                     
         dmat(i,5,2) = zero                                                     
         dmat(i,5,3) = zero                                                     
         dmat(i,5,4) = zero                                                     
         dmat(i,5,6) = zero                                                     
         dmat(i,6,1) = zero                                                     
         dmat(i,6,2) = zero                                                     
         dmat(i,6,3) = zero                                                     
         dmat(i,6,4) = zero                                                     
         dmat(i,6,5) = zero                                                     
c                                                                               
         fact = weight * det_jac_block(i)                                       
         c1 = (e_vec(i)/((one+nu_vec(i))*(one-two*nu_vec(i))))*fact             
         c2 = (one-nu_vec(i))*c1                                                
         c3 = ((one-two*nu_vec(i))/two)*c1                                      
         c4 = nu_vec(i)*c1                                                      
c                                                                               
         dmat(i,1,1)= c2                                                        
         dmat(i,2,2)= c2                                                        
         dmat(i,3,3)= c2                                                        
         dmat(i,4,4)= c3                                                        
         dmat(i,5,5)= c3                                                        
         dmat(i,6,6)= c3                                                        
         dmat(i,1,2)= c4                                                        
         dmat(i,1,3)= c4                                                        
         dmat(i,2,1)= c4                                                        
         dmat(i,3,1)= c4                                                        
         dmat(i,2,3)= c4                                                        
         dmat(i,3,2)= c4                                                        
c                                                                               
      end do                                                                    
c                                                                               
      debug = .false.                                                           
      if( debug .and. felem .eq. 1 .and. gpn .eq. 5 ) then                      
       write(iout,*) '... stresses at n+1, element 1, gp 5...'                  
       write(iout,9100) stress_np1(1,1:6)                                       
      end if                                                                    
c                                                                               
      return                                                                    
 9100 format(3x,6f10.3)                                                         
      end                                                                       
c                                                                               
c *******************************************************************           
c *                                                                 *           
c *        material model # 9 -- <available>                        *           
c *                                                                 *           
c *           set 3 material model dependent output values          *           
c *                                                                 *           
c *******************************************************************           
c                                                                               
c                                                                               
      subroutine oumm09( gpn, mxvl, span, iout, elestr,                         
     &                   stress, history )                                      
      implicit none                                                             
c                                                                               
c                   parameter declarations                                      
c                   ----------------------                                      
c                                                                               
      integer                                                                   
     &  gpn, mxvl, span, iout                                                   
c                                                                               
c                                                                               
      double precision                                                          
     & stress(mxvl,*), elestr(mxvl,*), history(mxvl,*)                          
c                                                                               
c               description of parameters                                       
c               -------------------------                                       
c                                                                               
c     gpn               : gauss point number being processed for block          
c     mxvl              : maximum no. elements per block                        
c     span              : number of elements in current block                   
c     iout              : write messages to this device number                  
c     stress            : current stresses for all                              
c                         elements in block for this gauss point                
c (*) elestr            : stresses to be output for elements                    
c     history           : current history values for all                        
c                         elements in block for this gauss point                
c                                                                               
c    (*)  values to be updated by this material model                           
c                                                                               
c                                                                               
c   stress ordering                elestr ordering                              
c     (1) sig-xx                   (1) sig-xx                                   
c     (2) sig-yy                   (2) sig-yy                                   
c     (3) sig-zz                   (3) sig-zz                                   
c     (4) tau-xy                   (4) tau-xy                                   
c     (5) tau-yz                   (5) tau-yz                                   
c     (6) tau-xz                   (6) tau-xz                                   
c     (7) total work density       (7) total work density                       
c                                  (8) mises equiv. stress                      
c                             (*)  (9) mat_val1                                 
c                             (*) (10) mat_val2                                 
c                             (*) (11) mat_val3                                 
c                                                                               
c  NOTE:  do not modify "stress" array or columns 1-8 of the "elestr"           
c         array. only modify columns 9-11 of "elestr". These are the            
c         3 "material model" dependent values output in the stress              
c         values for a gauss point. See Section 2.12 of manual and              
c         description of each material model. The output labels for             
c         columns 9-11 are "c1", "c2", "c3".                                    
c                                                                               
c                                                                               
c                                                                               
c                                                                               
       return                                                                   
       end                                                                      
                                                                                
                                                                                
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm09_set_sizes                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified: 7/1/12  rhd                 *          
c     *                                                              *          
c     *    called by warp3d for each material model to obtain        *          
c     *    various sizes of data for the model                       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm09_set_sizes( info_vector )                                  
      dimension info_vector(*)                                                  
c                                                                               
c        set infor_data                                                         
c                                                                               
c         1        number of history values per integration                     
c                  point. Abaqus calls these "statev". Values                   
c                  double or single precsion based on hardware.                 
c                                                                               
c         2        number of values in the symmetric part of the                
c                  [D] for each integration point. for solid                    
c                  elements this is 21, for cohesive elements this 6.           
c                                                                               
c         3        = 0, the material model returns "unrotated"                  
c                       Cauchy stresses at n+1                                  
c                  = 1, the material model returns the standard                 
c                       Cauchy stresses at n+1                                  
c                                                                               
c         4        number of state variables per point to be output             
c                  when user requests this type of results                      
c                                                                               
      info_vector(1) = 1                                                        
      info_vector(2) = 21                                                       
      info_vector(3) = 0                                                        
      info_vector(4) = 0                                                        
c                                                                               
      return                                                                    
      end                                                                       
c            dummy routines for model not yet supporting                        
c            states output                                                      
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *             subroutine mm09_states_values                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *               last modified : 1/3/2015 (rhd))                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm09_states_values( itype, elem_states_output,                 
     &                                 nrow_states, num_states  )               
      use global_data ! old common.main
c                                                                               
c                       access some global data structures                      
c                                                                               
      use elem_block_data, only: history_blocks, history_blk_list               
      use main_data, only: elems_to_blocks, cohesive_ele_types                  
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c                       parameters                                              
c                                                                               
      integer :: nrow_states, itype, num_states                                 
      double precision :: elem_states_output(nrow_states,*)                     
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine mm09_states_labels                *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *               last modified : 1/11/2015 (rhd)                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine mm09_states_labels( size_state,                                
     &      num_states, state_labels, state_descriptors, out,                   
     &      comment_lines, max_comment_lines, num_comment_lines )               
      implicit none                                                             
c                                                                               
c                       parameters                                              
c                                                                               
      integer :: size_state, num_states, out, max_comment_lines,                
     &           num_comment_lines                                              
      character(len=8)  :: state_labels(size_state)                             
      character(len=60) :: state_descriptors(size_state)                        
      character(len=80) :: comment_lines(max_comment_lines)                     
c                                                                               
c                       locals                                                  
c                                                                               
      integer :: i                                                              
c                                                                               
      num_states = 0                                                            
      num_comment_lines = 0                                                     
      state_labels(1) = "..."                                                   
      state_descriptors(1) = "...."                                             
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
