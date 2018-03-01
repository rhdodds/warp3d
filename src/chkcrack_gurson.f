c ********************************************************************          
c *                                                                  *          
c *  routines to support crack growth by gurson criterion            *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine dam_param_gt                   *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 09/22/02                   *          
c     *                                                              *          
c     *    evaluate current damage state for an element associated   *          
c     *    with the standard GT material model (#3)                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine dam_param_gt( elem, kill_now, debug, porosity,                 
     &                         sig_mean, sig_mises )                            
      use global_data ! old common.main
c                                                                               
c        elem      -- (input)   element number to be checked for                
c                               killing                                         
c        kill_now  -- (output)  set .true. if element should be                 
c                               killed now                                      
c        debug     -- (input)   .true. if the calling routine wants             
c                               debugging info output                           
c        porosity  -- (output)  for Gurson-type material models,                
c                               this is the average element porosity.           
c                               Just used for output messages                   
c        sig_mean  -- (output)  all models. average mean stress                 
c                               over element. just used for output              
c                               messages                                        
c        sig_mises -- (ouput)   all models. average mises stress                
c                               over element. just used for output.             
c                                                                               
c                                                                               
      use main_data,       only : elems_to_blocks                               
      use elem_block_data, only : history_blocks, urcs_n_blocks,                
     &                            history_blk_list                              
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c              parameter declarations                                           
c                                                                               
      logical kill_now, debug                                                   
      double precision                                                          
     &     porosity, sig_mean, sig_mises                                        
      double precision,                                                         
     &  dimension(:), pointer :: history, urcs_n                                
c                                                                               
c              local declarations                                               
c                                                                               
      double precision                                                          
     &     zero,                                                                
     &     third, six, iroot2,                                                  
     &     sig_xx, sig_yy, sig_zz,                                              
     &     sig_xy, sig_xz, sig_yz,eps_xx, eps_yy, eps_zz,                       
     &     gam_xy, gam_xz, gam_yz, fpngp                                        
      data zero /0.0/                                                           
      data third / 0.333333333 /                                                
      data iroot2, six / 0.70711, 6.0 /                                         
c                                                                               
c          set up to process element. get number of gauss points.               
c          find the block that the element belongs to. get the                  
c          relative number of element within the block. used to                 
c          index into data structures stored by block (from pointers)           
c          set pointers to stress and history blocks that contain this          
c          element.                                                             
c                                                                               
      porosity    = zero                                                        
      sig_mean    = zero                                                        
      sig_mises   = zero                                                        
      ngp         = iprops(6,elem)                                              
      blk         = elems_to_blocks(elem,1)                                     
      rel_elem    = elems_to_blocks(elem,2)                                     
      hist_size   = history_blk_list(blk)                                       
      offset      = (rel_elem-1)*hist_size*ngp + 1                              
      urcs_n      => urcs_n_blocks(blk)%ptr                                     
      history     => history_blocks(blk)%ptr                                    
c                                                                               
c          loop over the element gauss points to compute the average            
c          value of porosity, mean stress and (macro) mises stress.             
c                                                                               
c          at each gauss point:                                                 
c              history(5) = current porosity                                    
c          urcs_n: current values of unrotated cauchy stress -note              
c                  component ordering                                           
c                  nstrs: number of stress values stored per gauss point        
c       hist_size: number of history values stored per gauss point              
c                                                                               
      do gp = 1, ngp                                                            
        sigoffset   = (rel_elem-1)*nstrs*ngp + (gp-1)*nstrs                     
        sig_xx = urcs_n(sigoffset + 1)                                          
        sig_yy = urcs_n(sigoffset + 2)                                          
        sig_zz = urcs_n(sigoffset + 3)                                          
        sig_xy = urcs_n(sigoffset + 4)                                          
        sig_yz = urcs_n(sigoffset + 5)                                          
        sig_xz = urcs_n(sigoffset + 6)                                          
        sig_mean  = sig_mean + (sig_xx + sig_yy + sig_zz )                      
        sig_mises = sig_mises +                                                 
     &           sqrt( (sig_xx-sig_yy)**2 + (sig_yy-sig_zz)**2 +                
     &           (sig_xx-sig_zz)**2 + six*( sig_xy**2 +                         
     &           sig_yz**2 + sig_xz**2 ) ) * iroot2                             
        hist_offset = offset+4+(gp-1)*hist_size                                 
        porosity = porosity + history(hist_offset)                              
      end do                                                                    
c                                                                               
c          final average values over element. check average porosity            
c          against user specified limit (in module damage_data)                 
c                                                                               
       fpngp = dble( ngp )                                                      
      porosity = porosity / fpngp                                               
      sig_mean  = sig_mean * third / fpngp                                      
      sig_mises = sig_mises / fpngp                                             
      kill_now  =  porosity .gt. porosity_limit                                 
c                                                                               
      if ( debug ) then                                                         
            write(out,*) 'element is:',elem                                     
            write(out,'("    porosity=",e14.6)') porosity                       
            write(out,'("    sig_mean=",e14.6)') sig_mean                       
            write(out,'("    sig_mises=",e14.6)') sig_mises                     
      end if                                                                    
      return                                                                    
c                                                                               
      end                                                                       
                                                                                
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                    subroutine dam_param_agt                  *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 06/12/03                   *          
c     *                                                              *          
c     *    evaluate current damage state for an element associated   *          
c     *    with the extended GT material model (#6)                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine dam_param_agt( elem, kill_now, debug, porosity,                
     &                          sig_mean, sig_mises, ext_spacing,               
     &                          ext_shape )                                     
      use global_data ! old common.main
c                                                                               
c   elem        -- (input)   element number to be checked for                   
c                            killing                                            
c   kill_now    -- (output)  set .true. if element should be                    
c                            killed now                                         
c   debug       -- (input)   .true. if the calling routine wants                
c                            debugging info output                              
c   porosity    -- (output)  for Gurson-type material models,                   
c                            this is the average element porosity.              
c                            Just used for output messages                      
c   sig_mean    -- (output)  all models. average mean stress                    
c                            over element. just used for output                 
c                            messages                                           
c   sig_mises   -- (ouput)   all models. average mises stress                   
c                            over element. just used for output.                
c   ext_shape   -- (output)  average value of W (shape) parameter               
c                            over element                                       
c   ext_spacing -- (output)  average value of X (spacing) parameter             
c                            over element                                       
c                                                                               
c                                                                               
      use main_data,       only : elems_to_blocks                               
      use elem_block_data, only : history_blocks, urcs_n_blocks,                
     &                            history_blk_list                              
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
c              parameter declarations                                           
c                                                                               
      logical kill_now, debug                                                   
      double precision                                                          
     &     porosity, sig_mean, sig_mises, ext_spacing, ext_shape                
      double precision,                                                         
     &  dimension(:), pointer :: history, urcs_n                                
c                                                                               
c              local declarations                                               
c                                                                               
      double precision                                                          
     &     zero, coalesce_count, one, half,                                     
     &     third, six, iroot2,                                                  
     &     sig_xx, sig_yy, sig_zz,                                              
     &     sig_xy, sig_xz, sig_yz,eps_xx, eps_yy, eps_zz,                       
     &     gam_xy, gam_xz, gam_yz, fpngp                                        
      data zero, one, half / 0.0, 1.0, 0.5 /                                    
      data third / 0.333333333 /                                                
      data iroot2, six / 0.70711, 6.0 /                                         
c                                                                               
c          set up to process element. get number of gauss points.               
c          find the block that the element belongs to. get the                  
c          relative number of element within the block. used to                 
c          index into data structures stored by block (from pointers)           
c          set pointers to stress and history blocks that contain this          
c          element.                                                             
c                                                                               
      porosity    = zero                                                        
      sig_mean    = zero                                                        
      sig_mises   = zero                                                        
      ext_shape   = zero                                                        
      ext_spacing = zero                                                        
      ngp         = iprops(6,elem)                                              
      blk         = elems_to_blocks(elem,1)                                     
      rel_elem    = elems_to_blocks(elem,2)                                     
      hist_size   = history_blk_list(blk)                                       
      offset      = (rel_elem-1)*hist_size*ngp + 1                              
      urcs_n      => urcs_n_blocks(blk)%ptr                                     
      history     => history_blocks(blk)%ptr                                    
c                                                                               
c          loop over the element gauss points to compute the average            
c          value of porosity, mean stress, (macro) mises stress,                
c          spacing parameter (X), shape parameter (W) and count of              
c          no. of points on coalesence surface.                                 
c                                                                               
c          at each gauss point:                                                 
c              history(1) = current porosity                                    
c              history(2) = matrix flow stress                                  
c              history(3) = void shape parameter (W)                            
c              history(4) = void spacing parameter (X)                          
c              history(5) = matrix plastic strain                               
c              history(6) = > 0.0 if point is plastic                           
c              history(7) = state flag (see mm06 for values)                    
c              history(8) = > 0.0 if point is on coalesence surface             
c                                                                               
c          urcs_n: current values of unrotated cauchy stress -note              
c                  component ordering                                           
c                  nstrs: number of stress values stored per gauss point        
c       hist_size: number of history values stored per gauss point              
c                  (same for all material models. value is set in               
c                  param_def)                                                   
c                                                                               
      coalesce_count = zero                                                     
      do gp = 1, ngp                                                            
        sigoffset   = (rel_elem-1)*nstrs*ngp + (gp-1)*nstrs                     
        sig_xx = urcs_n(sigoffset + 1)                                          
        sig_yy = urcs_n(sigoffset + 2)                                          
        sig_zz = urcs_n(sigoffset + 3)                                          
        sig_xy = urcs_n(sigoffset + 4)                                          
        sig_yz = urcs_n(sigoffset + 5)                                          
        sig_xz = urcs_n(sigoffset + 6)                                          
        sig_mean  = sig_mean + (sig_xx + sig_yy + sig_zz )                      
        sig_mises = sig_mises +                                                 
     &           sqrt( (sig_xx-sig_yy)**2 + (sig_yy-sig_zz)**2 +                
     &           (sig_xx-sig_zz)**2 + six*( sig_xy**2 +                         
     &           sig_yz**2 + sig_xz**2 ) ) * iroot2                             
        hist_offset = offset+(gp-1)*hist_size                                   
        porosity    = porosity + history(hist_offset+0)                         
        ext_spacing = ext_spacing + history(hist_offset+2)                      
        ext_shape   = ext_shape + history(hist_offset+3)                        
        if( history(hist_offset+7) .gt. zero ) coalesce_count =                 
     &       coalesce_count + one                                               
      end do                                                                    
c                                                                               
c          final average values over element. check number of points            
c          on coalesence surface to set kill flag                               
c                                                                               
       fpngp = dble( ngp )                                                      
      porosity    = porosity / fpngp                                            
      sig_mean    = sig_mean * third / fpngp                                    
      sig_mises   = sig_mises / fpngp                                           
      ext_spacing = ext_spacing / fpngp                                         
      ext_shape   = ext_shape  / fpngp                                          
      kill_now    = coalesce_count / fpngp .ge. half                            
c                                                                               
      if ( debug ) then                                                         
         write(out,*) 'element is:',elem                                        
         write(out,'("    porosity=",e14.6)') porosity                          
         write(out,'("    sig_mean=",e14.6)') sig_mean                          
         write(out,'("    sig_mises=",e14.6)') sig_mises                        
         write(out,'("    ext_shape=",e14.6)') ext_shape                        
         write(out,'("    ext_spacing=",e14.6)') ext_spacing                    
         write(out,'("    coalesce_count=",f5.1)') coalesce_count               
      end if                                                                    
      return                                                                    
c                                                                               
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_print_elem1              *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 6/11/03 rhd                *          
c     *                                                              *          
c     *     This routine prints out the status of killable Gurson    *          
c     *     and extended Gurson elements at the beginning of a       *          
c     *     load step                                                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine dam_print_elem1(step,iter)                                     
      use global_data ! old common.main
      use main_data                                                             
      use elem_extinct_data, only : dam_state, dam_print_list,                  
     &     old_porosity, old_mises, old_mean                                    
      use elem_block_data, only : history_blocks, history_blk_list              
      use damage_data                                                           
      implicit integer (a-z)                                                    
c                                                                               
c                                                                               
                                                                                
      double precision                                                          
     &     porosity, orig_porosity, zero, ebarp, sigma_bar, d_poros,            
     &     max_d_poros, fpngp, ddum1, ddum2, sig_mean, sig_mises,               
     &     ext_shape, ext_spacing                                               
      double precision,                                                         
     &  dimension(:), pointer :: history                                        
      logical ldummy, debug, all_killed, lmises(num_print_list),                
     &                                   lmean(num_print_list),                 
     &                                   ext_gurson                             
      data zero, debug / 0.0, .false. /                                         
      character(len=1) :: mises_char, mean_char                                 
c                                                                               
c           check to see if all elements in print list have been killed.        
c           if so, print a message and return                                   
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
      max_d_poros = zero                                                        
c                                                                               
c           print the header                                                    
c                                                                               
      write(out,*) ' '                                                          
      write(out,*) ' ********************************************* '            
      write(out,*) ' ***         Killable element status       *** '            
      write(out,*) ' ********************************************* '            
c                                                                               
      write(out,*)                                                              
      if( load_size_control_crk_grth ) then                                     
         write(out,9020)                                                        
      else                                                                      
         write(out,9030)                                                        
      end if                                                                    
c                                                                               
c          parse the order list                                                 
c                                                                               
      local_count = 0                                                           
      lmises = .false.                                                          
      lmean  = .false.                                                          
      do elem_loop = 1, num_print_list                                          
       element  = dam_print_list(elem_loop)                                     
         elem_ptr = dam_ptr(element)                                            
c                                                                               
c             check if element is a killable element and if it has              
c             been killed.                                                      
c                                                                               
         if( dam_ptr(element) .eq. 0 ) cycle                                    
         if( dam_state(elem_ptr) .ne. 0 ) cycle                                 
c                                                                               
c             average porosity values.  print if new porosity is                
c             different than the old porosity.                                  
c                                                                               
         porosity  = zero                                                       
         sig_mean  = zero                                                       
         sig_mises = zero                                                       
         ngp       = iprops(6,element)                                          
         fpngp     = ngp                                                        
         blk       = elems_to_blocks(element,1)                                 
         rel_elem  = elems_to_blocks(element,2)                                 
         hist_size = history_blk_list(blk)                                      
         offset    = (rel_elem-1)*hist_size*ngp + 1                             
         history   => history_blocks(blk)%ptr                                   
c                                                                               
         call dam_param( element, ldummy, debug, porosity,                      
     &                   ddum1, ddum2, sig_mean, sig_mises,                     
     &                   ext_gurson, ext_shape, ext_spacing )                   
c                                                                               
c             if the mises stress or mean stress from this step are             
c             lower than the previous step, then put a * by the                 
c             value to indicate this state.                                     
c                                                                               
       if( sig_mises .lt. old_mises(elem_loop) ) then                           
          mises_char = '*'                                                      
            lmises(elem_loop)=.true.                                            
         else                                                                   
          mises_char = ' '                                                      
         end if                                                                 
c                                                                               
       if( sig_mean .lt. old_mean(elem_loop) ) then                             
          mean_char = '*'                                                       
            lmean(elem_loop)=.true.                                             
         else                                                                   
          mean_char = ' '                                                       
         end if                                                                 
c                                                                               
c             if we have porosity greater than the initial porosity,            
c             and if we have not yet killed the element, output the             
c             element data -- porosity, mean stress, etc. Also output           
c             change in porosity if we are using automatic load control.        
c                                                                               
         ebarp_loc  = 0                                                         
         sigbar_loc = 1                                                         
         if( ext_gurson ) then                                                  
           ebarp_loc  = 4                                                       
           sigbar_loc = 1                                                       
         end if                                                                 
c                                                                               
         orig_porosity = props(26,element)                                      
         if( porosity .ne. orig_porosity ) then                                 
            ebarp     = zero                                                    
            sigma_bar = zero                                                    
            do gp = 1, ngp                                                      
              hist_offset = offset+0+(gp-1)*hist_size                           
              ebarp = ebarp +                                                   
     &                history(hist_offset+ebarp_loc) / fpngp                    
              sigma_bar = sigma_bar +                                           
     &                history(hist_offset+sigbar_loc) / fpngp                   
            end do                                                              
c                                                                               
            if( load_size_control_crk_grth ) then                               
               d_poros = porosity - old_porosity(dam_ptr(element))              
               max_d_poros = max(max_d_poros, d_poros)                          
               write(out,9005) element, orig_porosity, porosity, ebarp,         
     &              sigma_bar, sig_mean, mean_char, sig_mises,                  
     &              mises_char, d_poros                                         
            else                                                                
               write(out,9000) element, orig_porosity, porosity, ebarp,         
     &              sigma_bar, sig_mean, mean_char, sig_mises,                  
     &              mises_char                                                  
            end if                                                              
c                                                                               
         local_count = local_count + 1                                          
c                                                                               
         end if                                                                 
c                                                                               
c             store mises and mean stresses for this step for comparison        
c             after next step                                                   
c                                                                               
       old_mises(elem_loop) = sig_mises                                         
       old_mean(elem_loop)  = sig_mean                                          
c                                                                               
      end do                                                                    
c                                                                               
c           packet output for gurson cell status                                
c           ====================================                                
c                                                                               
c           repeat loop over num_print_list for packet output                   
c                                                                               
c                                                                               
c                                                                               
      if( local_count .ne. 0 .and. output_packets ) then                        
         if( load_size_control_crk_grth ) then                                  
           write( packet_file_no ) 6, local_count, step, iter                   
         else                                                                   
           write( packet_file_no ) 5, local_count, step, iter                   
         end if                                                                 
c                                                                               
c                                                                               
      do elem_loop = 1, num_print_list                                          
       element  = dam_print_list(elem_loop)                                     
         elem_ptr = dam_ptr(element)                                            
c                                                                               
c             check if element is a killable element and if it has              
c             been killed.                                                      
c                                                                               
         if( dam_ptr(element) .eq. 0 ) cycle                                    
         if( dam_state(elem_ptr) .ne. 0 ) cycle                                 
c                                                                               
c             average porosity values.  print if new porosity is                
c             different than the old porosity.                                  
c                                                                               
         porosity  = zero                                                       
         sig_mean  = zero                                                       
         sig_mises = zero                                                       
         ngp       = iprops(6,element)                                          
         fpngp     = ngp                                                        
         blk       = elems_to_blocks(element,1)                                 
         rel_elem  = elems_to_blocks(element,2)                                 
         hist_size = history_blk_list(blk)                                      
         offset    = (rel_elem-1)*hist_size*ngp + 1                             
         history   => history_blocks(blk)%ptr                                   
c                                                                               
         call dam_param( element, ldummy, debug, porosity,                      
     &                   ddum1, ddum2, sig_mean, sig_mises,                     
     &                   ext_gurson, ext_shape, ext_spacing )                   
c                                                                               
c             if the mises stress or mean stress from this step are             
c             lower than the previous step, then put a 0 or 1 flag              
c             value to indicate this state.                                     
c                                                                               
         if( lmises(elem_loop) )then                                            
            mises_flag = 0                                                      
         else                                                                   
            mises_flag = 1                                                      
         end if                                                                 
c                                                                               
         if( lmean(elem_loop) )then                                             
            mean_flag = 0                                                       
         else                                                                   
            mean_flag = 1                                                       
         end if                                                                 
c                                                                               
c                                                                               
c             if we have porosity greater than the initial porosity,            
c             and if we have not yet killed the element, output the             
c             element data -- porosity, mean stress, etc. Also output           
c             change in porosity if we are using automatic load control.        
c                                                                               
         ebarp_loc  = 0                                                         
         sigbar_loc = 1                                                         
         if( ext_gurson ) then                                                  
           ebarp_loc  = 4                                                       
           sigbar_loc = 1                                                       
         end if                                                                 
c                                                                               
         orig_porosity = props(26,element)                                      
         if( porosity .ne. orig_porosity ) then                                 
            ebarp     = zero                                                    
            sigma_bar = zero                                                    
            do gp = 1, ngp                                                      
              hist_offset = offset+0+(gp-1)*hist_size                           
              ebarp = ebarp +                                                   
     &                history(hist_offset+ebarp_loc) / fpngp                    
              sigma_bar = sigma_bar +                                           
     &                history(hist_offset+sigbar_loc) / fpngp                   
            end do                                                              
c                                                                               
            if( load_size_control_crk_grth ) then                               
               d_poros = porosity - old_porosity(dam_ptr(element))              
               max_d_poros = max(max_d_poros, d_poros)                          
               write(packet_file_no) element, orig_porosity,                    
     &              porosity, ebarp,                                            
     &              sigma_bar, sig_mean, mean_flag, sig_mises,                  
     &              mises_flag, d_poros                                         
            else                                                                
               write(packet_file_no) element, orig_porosity,                    
     &               porosity, ebarp,                                           
     &              sigma_bar, sig_mean, mean_flag, sig_mises,                  
     &              mises_flag                                                  
            end if                                                              
c                                                                               
         end if                                                                 
c                                                                               
c                                                                               
      end do                                                                    
c                                                                               
c                                                                               
      end if                                                                    
c                                                                               
c            end of packet output loop                                          
c            =========================                                          
c                                                                               
      write(out,9015) num_elements_killed                                       
c                                                                               
c         print out information about the largest increase in porosity          
c         over the last load step. get the original porosity for the            
c         first element in the print list                                       
c                                                                               
      element       = dam_print_list(1)                                         
      orig_porosity = props(26,element)                                         
c                                                                               
c           now write out info about change in porosity                         
c                                                                               
      if( load_size_control_crk_grth ) then                                     
         write(out,9010) max_d_poros, 100.0*max_d_poros/porosity_limit,         
     &        100.0*max_d_poros/orig_porosity                                   
      else                                                                      
         write(out,*) ' '                                                       
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format(1x,i7,1x,f12.5,3x,f12.5,3(1x,e14.6),a1,(1x,e14.6),a1)              
 9005 format(1x,i7,1x,f12.5,3x,f12.5,                                           
     &       e14.6,                                                             
     &       2x,e14.6,                                                          
     &       3x,e14.6,a1,                                                       
     &       1x,e14.6,a1,                                                       
     &       1x,f10.6 )                                                         
 9010 format(5x,'>> max. delta-f for step:',f8.5,                               
     &      '.  f-critical (%): ',f7.3,                                         
     &     '. f-initial (%): ',f7.3,/)                                          
 9015 format(/,5x,'>> * = value has decreased since last step.',                
     &       /,5x,'>> total number of elements killed: ',i7)                    
 9020 format(1x,'element     inital f     current f       Ep    ',              
     &       '        sigma bar      mean stress    ',                          
     &       'mises stress     delta f        ',/,1x,                           
     &       '-------     --------     ---------       --    ',                 
     &       '        ---------      -----------    ',                          
     &       '------------     -------        ')                                
 9030 format(1x,'element     inital f     current f       Ep    ',              
     &       '      sigma bar      mean stress    ',                            
     &       'mises stress',/,1x,                                               
     &       '-------     --------     ---------       --    ',                 
     &       '      ---------      -----------    ',                            
     &       '------------')                                                    
 9040 format(1x,'*** NOTE: all elements in the killable element',               
     & ' status list have been killed.',/)                                      
c                                                                               
      end                                                                       
