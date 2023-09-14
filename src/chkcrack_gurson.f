c ********************************************************************          
c *                                                                  *          
c *  routines to support crack growth by gurson criterion            *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *         subroutine dam_param_gt_get_values                   *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 4/21/23 rhd                *          
c     *                                                              *          
c     *    evaluate current damage state for an element associated   *          
c     *    with the standard GT material model (#3)                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine dam_param_gt_get_values( elem, kill_now, porosity,                 
     &                         sig_mean, sig_mises, eps_plastic,
     &                         ebarp, sigma_bar )                            
c                                                                               
c        elem      -- (input)   element number to be checked for                
c                               killing                                         
c        kill_now  -- (output)  set .true. if element should be                 
c                               killed now                                      
c        porosity  -- (output)  average element porosity.           
c        sig_mean  -- (output)  average mean stress over element
c        sig_mises -- (ouput)   average mises stress over element
c        eps_plastic -- (output) average (macroscopic) plastic
c                               strain over element
c        ebarp      -- (output)  average matrix plastic
c                               strain over element
c        sigma_bar  -- (output) average matrix stress over element                               over element
c                                                                                                                                                           
      use global_data,     only : iprops, nstrs, out
      use main_data,       only : elems_to_blocks                               
      use elem_block_data, only : history_blocks, urcs_n_blocks,                
     &                            history_blk_list                              
      use damage_data,     only : porosity_limit   
      use constants                                                        
c                                                                               
      implicit none                                                    
c                                                                               
c              parameter declarations                                           
c               
      integer, intent(in) :: elem                                                                
      logical, intent(out):: kill_now                                                   
      double precision, intent(out) :: porosity, sig_mean, sig_mises, 
     &                                 eps_plastic, ebarp, sigma_bar
c                                        
c              local declarations                                               
c      
      integer :: ngp, blk, rel_elem, hist_size, offset, gp, sigoffset,
     &           hist_offset, ebarp_loc, sigbar_loc, ho
      logical :: local_debug
      double precision, dimension(:), pointer :: history, urcs_n                                
      double precision :: sig_xx, sig_yy, sig_zz,                                              
     &     sig_xy, sig_xz, sig_yz,eps_xx, eps_yy, eps_zz,                       
     &     gam_xy, gam_xz, gam_yz, fpngp       
c                                                                               
c          set up to process element. get number of gauss points.               
c          find the block that the element belongs to. get the                  
c          relative number of element within the block. used to                 
c          index into data structures stored by block (from pointers)           
c          set pointers to stress and history blocks that contain this          
c          element.                                                             
c  
      local_debug = elem == -1
      porosity    = zero                                                        
      sig_mean    = zero                                                        
      sig_mises   = zero   
      eps_plastic = zero   
      ebarp       = zero ! matrix plastic strain
      sigma_bar   = zero ! matrix stress                                                  
      ngp         = iprops(6,elem)                                              
      blk         = elems_to_blocks(elem,1)                                     
      rel_elem    = elems_to_blocks(elem,2)                                     
      hist_size   = history_blk_list(blk)                                       
      offset      = (rel_elem-1)*hist_size*ngp + 1                              
      urcs_n      => urcs_n_blocks(blk)%ptr                                     
      history     => history_blocks(blk)%ptr                                    
      ebarp_loc   = 0       ! history offsets                                                  
      sigbar_loc  = 1                                                         
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
        sigoffset = (rel_elem-1)*nstrs*ngp + (gp-1)*nstrs                     
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
        eps_plastic = eps_plastic + urcs_n(sigoffset + 9) 
        ho = offset+0+(gp-1)*hist_size                           
        ebarp = ebarp + history(ho+ebarp_loc)
        sigma_bar = sigma_bar + history(ho+sigbar_loc)
      end do                                                                    
c                                                                               
c          final average values over element. check average porosity            
c          against user specified limit (in module damage_data)                 
c                                                                               
      fpngp = dble( ngp )                                                      
      porosity    = porosity / fpngp                                               
      sig_mean    = sig_mean * third / fpngp                                      
      sig_mises   = sig_mises / fpngp
      eps_plastic = eps_plastic / fpngp  
      ebarp       = ebarp / fpngp        
      sigma_bar   = sigma_bar / fpngp                              
      kill_now    =  porosity .gt. porosity_limit                                 
c                                                                               
      if( local_debug ) then                                                         
            write(out,*) 'element is:',elem                                     
            write(out,'("    porosity=",e14.6)') porosity                       
            write(out,'("    sig_mean=",e14.6)') sig_mean                       
            write(out,'("    sig_mises=",e14.6)') sig_mises                     
      end if 
c                                                                   
      return                                                                    
c                                                                               
      end  subroutine dam_param_gt_get_values                                                                
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine dam_print_gt                      *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 9/10/2023 rhd              *          
c     *                                                              *          
c     *     prints the status of killable Gurson and extended Gurson *          
c     *     elements at the beginning of a load step                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine dam_print_gt( step, iter )                                     
c
      use global_data, only : out, noelem, iprops, props, lprops
      use main_data,   only : output_packets, packet_file_no   
      use elem_extinct_data, only : dam_state, dam_print_list,                  
     &                              gt_old_porosity, Oddy_metrics                                    
      use damage_data, only : dam_ptr, load_size_control_crk_grth,
     &                        num_elements_killed, porosity_limit,
     &                        num_print_list, use_distortion_metric,
     &                        print_status, print_top_list,
     &                        num_top_list, num_kill_elem 
      use constants
c                                                        
      implicit none                                                                      
c                                                                               
c              parameter declarations                                           
c
      integer, intent(in) :: step, iter
      include 'include_damage_values'
      type(values_T) :: gt_elem_values
c                                        
c              local declarations                                               
c  
      integer :: elem_loop, element, elem_ptr, local_count,
     &           i, list_count, list_length, first_element     
      double precision :: porosity, orig_porosity, ebarp, sigma_bar,
     &     d_poros, max_d_poros, sig_mean, max_Oddy_ratio,             
     &     sig_mises, avg_eps_plas, max_oddy
      logical :: debug, write_packets, all_killed,  geo_non_flg,
     &           print_distort_metric   
      logical, save :: print_info = .true.  
c
      double precision, allocatable :: damage_factors(:)
      integer, allocatable :: element_list(:), sort_index_vec(:)
c                                                                               
      debug = .false.
c
c              are all elements in the user specified print list killed.        
c              if so, print a message and return                                   
c
      if( print_status ) then  ! user gave list of elements                                                                            
        all_killed = .true.                                                       
        do elem_loop = 1, num_print_list                                          
           element  = dam_print_list(elem_loop)                                   
           elem_ptr = dam_ptr(element)                                            
           if( elem_ptr .eq. 0 ) cycle    ! not killable                                     
           if( dam_state(elem_ptr) .ne. 0 ) cycle ! already killed                                 
           all_killed = .false.                                                   
        end do                                                                    
        if( all_killed ) then                                                     
           write(out,9040)                                                        
           return                                                                 
        end if      
      end if                                                                      
c
c              Support for printing elements with largest damage
c              factors. build sorted list of elements.              
c          
      if( print_top_list ) then                                                          
        allocate( damage_factors(num_kill_elem),
     &            element_list(num_kill_elem),
     &            sort_index_vec(num_kill_elem) ) 
        damage_factors = zero
        sort_index_vec = 0
        list_count     = 0
c
        do element = 1, noelem
          elem_ptr = dam_ptr(element) 
          if( elem_ptr == 0 ) cycle  ! element not killable
          if( dam_state(elem_ptr) .ne. 0 ) cycle ! already killed
          list_count = list_count + 1
          element_list(list_count) = element
        end do  
c$OMP PARALLEL DO PRIVATE( i, element, elem_ptr, gt_elem_values ) 
c$OMP&          IF( noelem > 10000)
        do i = 1, list_count
          element = element_list(i)
          call dam_param_gt( element, gt_elem_values )
          damage_factors(i) = gt_elem_values%porosity
        end do
c$OMP END PARALLEL DO
        call hpsort_eps_epw( list_count, damage_factors, 
     &                        sort_index_vec, 1.0d-07 )
      end if ! print_top_list
c
      write(out,*) ' '                                                          
      write(out,*) ' ********************************************* '            
      write(out,*) ' ***         Killable element status       *** '            
      write(out,*) ' ********************************************* '            
c                                                                               
      write(out,*) " "                                                             
      if( load_size_control_crk_grth ) write(out,9020)                                                        
      if( .not. load_size_control_crk_grth ) write(out,9030)                                                        
c                                                                               
c              print element status for GTN. follow user-specified
c              print order or the sorted list by descending porosity
c                                                                               
      local_count = 0                                                           
      max_d_poros = zero 
      max_oddy = minus_one 
      if( print_status ) then
          list_length = num_print_list
          first_element = dam_print_list(1)  ! to get initial porosity
      end if               
      if( print_top_list ) then
          list_length = min( num_top_list, list_count )
          first_element = element_list(1)
      end if          
c                                                                
      do elem_loop = 1, list_length
c      
         if( print_top_list ) element = 
     &              element_list(sort_index_vec(elem_loop))
         if( print_status ) then
            element = dam_print_list(elem_loop)
            elem_ptr = dam_ptr(element)    
            if( elem_ptr == 0 ) cycle   ! not killable 
            if( dam_state(elem_ptr) .ne. 0 ) cycle ! already killed
         end if  
         elem_ptr = dam_ptr(element)                                  
         call dam_param_gt( element, gt_elem_values )
         sig_mises    = gt_elem_values%sig_mises         
         sig_mean     = gt_elem_values%sig_mean
         ebarp        = gt_elem_values%gt_matrix_eps_plas
         sigma_bar    = gt_elem_values%gt_matrix_stress
         avg_eps_plas = gt_elem_values%eps_plas
         porosity     = gt_elem_values%porosity
         geo_non_flg  = lprops(18,element)    
c                                                                               
c             if we have porosity greater than the initial porosity,            
c             and if we have not yet killed the element, output the             
c             element data -- porosity, mean stress, etc. 
c         
         orig_porosity = props(26,element)                                      
         if( porosity .le. orig_porosity ) cycle
c         
         print_distort_metric = .false.          
         if( use_distortion_metric .and. geo_non_flg ) then
                 max_Oddy_ratio = Oddy_metrics(elem_ptr,2) /
     &                            Oddy_metrics(elem_ptr,1)     
                 max_oddy = max( max_oddy, max_Oddy_ratio )
                 print_distort_metric = .true.
         end if                                             
c                                                             
         if( load_size_control_crk_grth ) then                               
               d_poros = porosity - gt_old_porosity(dam_ptr(element))              
               max_d_poros = max(max_d_poros, d_poros)  
               if( .not.  print_distort_metric )                        
     &          write(out,9005) element, orig_porosity, porosity, ebarp,         
     &              sigma_bar, sig_mean, sig_mises,                  
     &              d_poros, avg_eps_plas                                       
               if( print_distort_metric )                        
     &          write(out,9005) element, orig_porosity, porosity, ebarp,         
     &              sigma_bar, sig_mean, sig_mises,                  
     &              d_poros, avg_eps_plas, max_Oddy_ratio                                       
         else 
              if( .not. print_distort_metric )                        
     &         write(out,9000) element, orig_porosity, porosity, ebarp,         
     &              sigma_bar, sig_mean, sig_mises,                  
     &              avg_eps_plas                                                  
             if( print_distort_metric )                        
     &         write(out,9000) element, orig_porosity, porosity, ebarp,         
     &              sigma_bar, sig_mean,sig_mises,                  
     &              avg_eps_plas, max_Oddy_ratio                                                  
         end if                                                              
c                                                                               
         local_count = local_count + 1                                          
c                                                                               
      end do     ! elem_loop
c   
      write_packets = local_count .ne. 0 .and. output_packets 
      if( write_packets ) call dam_print_gt_packets
c
      if( print_top_list ) write(out,9018)
      write(out,9015) num_elements_killed  
      if( print_info ) then
        write(out,9016)  
        print_info = .false.
      end if                                     
c                                                                               
c         information about the largest increase in porosity          
c         over the last load step. get the original porosity for the            
c         first element in the print list                                       
c                                                                               
      orig_porosity = props(26,first_element)    
c
      if( use_distortion_metric .and. geo_non_flg ) 
     &         write(out,9050) max_oddy                                     
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
 9000 format(1x,i7,1x,f12.5,3x,f11.5,1x,f10.4,4x,2(1x,g14.6),1x,
     & (1x,g14.6),1x,
     &        4x,f7.5,6x,f7.3)
 9005 format(1x,i7,1x,f12.5,3x,f12.5,                                           
     &       f10.4,4x,                                                             
     &       2x,g14.6,                                                          
     &       3x,g14.6,1x,                                                       
     &       1x,g14.6,1x,                                                       
     &       f10.5,8x,f7.5,4x,f7.3 )                                                         
 9010 format(5x,'>> max. delta-f for step:',f8.5,                               
     &      '.  f-critical (%): ',f7.3,                                         
     &     '. f-initial (%): ',f10.3,/)                                          
 9015 format(5x,'>> total number of elements killed: ',i7) 
 9016 format(5x,'f: porosity',
     &  /,5x,'Ep: matrix plastic strain',
     &  /,5x,'sigma_bar: matrix stress',
     &  /,5x,'mean stress, mises_stress, eps_plastic: ',
     &  'element (macroscale) averages')     
 9018 format(/,5x,'>> ordering by decreasing porosity')                    
 9020 format(1x,'element     inital f     current f       Ep    ',              
     &       '        sigma bar      mean stress    ',                          
     &       'mises stress     delta f    eps_plastic  Oddy ratio',                         
     &  /,1x '-------     --------     ---------       --    ',                 
     &       '        ---------      -----------    ',                          
     &       '------------     -------    -----------  ----------')                                
 9030 format(1x,'element     inital f     current f       Ep    ',              
     &       '      sigma bar      mean stress    ',                            
     &       'mises stress    eps_plastic  Oddy ratio',/,1x,                                               
     &       '-------     --------     ---------       --    ',                 
     &       '      ---------      -----------    ',                            
     &       '------------    -----------  ----------')    
 9040 format(1x,'*** NOTE: all elements in the killable element',               
     & ' status list have been killed.',/)      
 9050 format(5x,'>> maximum Oddy ratio: ',f6.3)
c 
      contains          
c     ========                                                                         
c                                     
      subroutine  dam_print_gt_packets   
c
      implicit none 
c
      integer :: packet_type, now_count   
      logical :: write_debug_text         
c
c               packet output for gurson cell status                                
c             
      write_debug_text = .false.
c
c              packet written only for elements that appeared
c              int the printed list just generated
c
      if( print_status )   list_length = num_print_list
      if( print_top_list ) list_length = min( num_top_list,
     &                                        list_count )
      
      
      packet_type = 5                                                                 
      if( load_size_control_crk_grth ) packet_type = 6 ! adaptive P control
      write( packet_file_no ) packet_type, local_count, step, iter  
      if( write_debug_text ) write(out,9010) packet_type, local_count, 
     &                       step, iter                                                            
c
      now_count = 0
c      
      do elem_loop = 1, list_length
c
         if( print_top_list ) element = 
     &              element_list(sort_index_vec(elem_loop))
         if( print_status ) then
            element = dam_print_list(elem_loop)
            elem_ptr = dam_ptr(element)    
            if( elem_ptr == 0 ) cycle   ! not killable 
            if( dam_state(elem_ptr) .ne. 0 ) cycle ! already killed
         end if     
c                                       
         call dam_param_gt( element, gt_elem_values )
         sig_mises    = gt_elem_values%sig_mises         
         sig_mean     = gt_elem_values%sig_mean
         ebarp        = gt_elem_values%gt_matrix_eps_plas
         sigma_bar    = gt_elem_values%gt_matrix_stress
         avg_eps_plas = gt_elem_values%eps_plas
         porosity     = gt_elem_values%porosity
c                                                                               
c             if we have porosity greater than the initial porosity,            
c             and if we have not yet killed the element, output the             
c             element data -- porosity, mean stress, etc. Also output           
c             change in porosity if we are using automatic load control.        
c                                                                                 
         orig_porosity = props(26,element)                                      
         if( porosity .le. orig_porosity ) cycle
         if( load_size_control_crk_grth ) then                               
           d_poros = porosity - gt_old_porosity(dam_ptr(element))              
           max_d_poros = max(max_d_poros, d_poros)                          
           write(packet_file_no) element, orig_porosity, porosity,
     &         ebarp, sigma_bar, sig_mean, sig_mises, d_poros  
           if( write_debug_text ) 
     &          write(out,9020) element, orig_porosity, porosity,
     &            ebarp, sigma_bar, sig_mean, sig_mises, d_poros                                                              
         else                                                                
           write(packet_file_no) element, orig_porosity,                    
     &           porosity, ebarp, sigma_bar, sig_mean, sig_mises
           if( write_debug_text ) 
     &          write(out,9030) element, orig_porosity, porosity,
     &            ebarp, sigma_bar, sig_mean, sig_mises                                                            
         end if         
c
         now_count = now_count + 1
c                                                                               
      end do   ! elem_loop    
      
      if( now_count .ne. local_count ) then
        write(out,9000) now_count, local_count
        call die_abort
      end if  
c
      return
c      
 9000 format(/,">>>> FATAL ERROR. routine dam_print_gt_packets",
     & /,15x,'now_count, local_count: ', 2i5 )     
 9010 format(/,".... binary packet file:",
     &  /,10x,"type, line count, step, iter: ",i2,i4,i8,i4) 
 9020 format(10x,"element, f_0, f, ebarp, sigbar, mean, mises, d(f): ",
     &           i7,f7.4,f7.4,f8.5,f8.2,f8.2,f8.2,f7.4 )
 9030 format(10x,"element, f_0, f, ebarp, sigbar, mean, mises: ",
     &    i7,f7.4,f7.4,f8.5, f8.2,f8.2,f8.2) 
c
      end subroutine dam_print_gt_packets

      end subroutine dam_print_gt
