c ********************************************************************          
c *                                                                  *          
c *  routines to support crack growth by smcs criterion              *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_print_elem3              *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 04/18/96                   *          
c     *                                                              *          
c     *     This routine prints out the status of SMCS elements      *          
c     *     marked as killable at the beginning of a load step       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dam_print_elem3( step, iter )                                  
      use global_data ! old common.main
      use elem_extinct_data, only : dam_state, dam_print_list,                  
     &     old_plast_strain                                                     
      use main_data, only : output_packets, packet_file_no                      
      use damage_data                                                           
      implicit integer (a-z)                                                    
c                                                                               
      double precision                                                          
     &    dummy, eps_plas, eps_crit, sig_mean, sig_mises, eps_plas_tol,         
     &    d_eps_plas, max_d_eps_plas, ddummy1, ddummy2, ddummy3                 
      logical ldummy, debug, ldummy2                                            
      data zero, debug, eps_plas_tol / 0.0, .false., 1.0e-9 /                   
c                                                                               
      max_d_eps_plas = zero                                                     
c                                                                               
c           print the header                                                    
c                                                                               
      write(out,*) ' '                                                          
      write(out,*) ' ********************************************* '            
      write(out,*) ' ***         Killable element status       *** '            
      write(out,*) ' ********************************************* '            
c                                                                               
      write(out,*)                                                              
      if ( load_size_control_crk_grth ) then                                    
         write(out,9020)                                                        
      else                                                                      
         write(out,9030)                                                        
      end if                                                                    
c                                                                               
c          parse the order list                                                 
c                                                                               
      local_count = 0                                                           
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
c             get current smcs parameters for element and print.                
c                                                                               
         ldummy2 = .false.                                                      
         call dam_param( element, ldummy, debug, dummy, eps_plas,               
     &                   eps_crit, sig_mean, sig_mises,                         
     &                   ldummy2, ddummy2, ddummy3 )                            
         if ( eps_plas .gt. eps_plas_tol ) then                                 
            if ( load_size_control_crk_grth ) then                              
               d_eps_plas = eps_plas -                                          
     &             old_plast_strain(dam_ptr(element))                           
               max_d_eps_plas = max(max_d_eps_plas, d_eps_plas)                 
               write(out,9000) element, eps_plas, eps_crit, sig_mean,           
     &              sig_mises, d_eps_plas                                       
            else                                                                
               write(out,9000) element, eps_plas, eps_crit, sig_mean,           
     &              sig_mises                                                   
            endif                                                               
            local_count = local_count + 1                                       
         end if                                                                 
      end do                                                                    
c                                                                               
c                                                                               
c                                                                               
c                                                                               
c          packet output section                                                
c                                                                               
c                                                                               
c                                                                               
      if ( output_packets .and. local_count .ne. 0 )then                        
c                                                                               
c                                                                               
         if (load_size_control_crk_grth) then                                   
           write( packet_file_no ) 23, local_count, step, iter                  
         else                                                                   
           write( packet_file_no ) 22, local_count, step, iter                  
         end if                                                                 
c                                                                               
c                                                                               
c                                                                               
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
c             get current smcs parameters for element and print.                
c                                                                               
         ldummy2 = .false.                                                      
         call dam_param( element, ldummy, debug, dummy, eps_plas,               
     &                   eps_crit, sig_mean, sig_mises,                         
     &                   ldummy2, ddummy2, ddummy3 )                            
         if ( eps_plas .gt. eps_plas_tol ) then                                 
            if ( load_size_control_crk_grth ) then                              
               d_eps_plas = eps_plas -                                          
     &             old_plast_strain(dam_ptr(element))                           
               max_d_eps_plas = max(max_d_eps_plas, d_eps_plas)                 
               write(packet_file_no) element, eps_plas, eps_crit,               
     &              sig_mean, sig_mises, d_eps_plas                             
            else                                                                
               write(packet_file_no) element, eps_plas, eps_crit,               
     &              sig_mean, sig_mises                                         
            endif                                                               
         end if                                                                 
      end do                                                                    
c                                                                               
c                                                                               
      end if                                                                    
c                                                                               
c                                                                               
c         end of packet output section                                          
c                                                                               
c                                                                               
c                                                                               
c                                                                               
      if ( load_size_control_crk_grth ) then                                    
         write (out,9010) max_d_eps_plas                                        
      else                                                                      
         write (out,*) ' '                                                      
      end if                                                                    
c                                                                               
      return                                                                    
c                                                                               
 9000 format(1x,i7,5(1x,e14.6))                                                 
 9010 format(/,5x,'>> Maximum change in plastic strain for this step:',         
     &     e13.6,/)                                                             
 9020 format(1x,'element   eps-pls        eps-crit      ',                      
     &       'sig-mean       sig-mises       d(eps-pls)',/,1x,                  
     &       '-------   -------        --------      ',                         
     &       '--------       ---------       ----------')                       
 9030 format(1x,'element   eps-pls        eps-crit      ',                      
     &       'sig-mean       sig-mises',/,1x,                                   
     &       '-------   -------        --------      ',                         
     &       '--------       ---------')                                        
c                                                                               
      end                                                                       
