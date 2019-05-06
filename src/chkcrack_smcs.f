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
c     *                   last modified : 04/28/2019 rhd             *          
c     *                                                              *          
c     *     This routine prints out the status of SMCS elements      *          
c     *     marked as killable at the beginning of a load step       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dam_print_elem3( step, iter )                                  
      use global_data, only : out  ! old common.main
      use elem_extinct_data, only : dam_state, dam_print_list,                  
     &                              old_plast_strain                                                     
      use main_data, only :   output_packets, packet_file_no                      
      use damage_data, only : dam_ptr, load_size_control_crk_grth,
     &                        num_print_list  
c                                                        
      implicit none                                                    
c
c              parameters 
c          
      integer :: step, iter
c
c              local declarations
c                                                                    
      integer :: elem_loop, elem_ptr, element, local_count, kill_count
      double precision :: dummy, eps_plas, eps_crit, sig_mean,
     &    sig_mises, d_eps_plas, max_d_eps_plas, ddummy1, ddummy2,
     &    ddummy3, triaxiality
      double precision, parameter :: zero = 0.0d0,
     &                               eps_plas_tol = 1.0d-06                 
      logical :: ldummy, debug, do_packets, skip_element
      character(len=1) :: cflag                                          
c                                                                               
      write(out,*) ' '                                                          
      write(out,*) ' ********************************************* '            
      write(out,*) ' ***         Killable element status       *** '            
      write(out,*) ' ********************************************* '            
      write(out,*)                                                              
      if( load_size_control_crk_grth ) then                                    
         write(out,9020)                                                        
      else                                                                      
         write(out,9030)                                                        
      end if                                                                    
c                                                                             
c              print element status for smcs. follow user-specified
c              print order
c                                                                               
      max_d_eps_plas = zero                                                     
      local_count = 0 
      kill_count = 0
c                                                          
      do elem_loop = 1, num_print_list                                          
         element  = dam_print_list(elem_loop)                                   
         elem_ptr = dam_ptr(element)                                            
c                                                                               
c              check if element is a killable element and/or if it has           
c              already been killed.                                              
c                                                                               
         if( dam_ptr(element) .eq. 0 ) cycle    
         if( dam_state(elem_ptr) .ne. 0 ) then
             kill_count = kill_count + 1
             cycle         
         end if    
c
c              get current smcs parameters for element and print.
c              skip printing elements with very small plastic strain                
c                                                                               
         call dam_param_3_get_values( element, debug, eps_plas,               
     &                   eps_crit, sig_mean, sig_mises,                         
     &                   triaxiality, 2, ldummy )  
         skip_element = eps_plas < eps_plas_tol
         if( skip_element ) cycle
         cflag = " "
         if( triaxiality < 1.0d0 ) cflag = "*"                  
         if( load_size_control_crk_grth ) then                              
            d_eps_plas = eps_plas - old_plast_strain(dam_ptr(element))    
            max_d_eps_plas = max(max_d_eps_plas, d_eps_plas)                 
            write(out,9110) element, eps_plas, eps_crit, sig_mean,           
     &              sig_mises, d_eps_plas, triaxiality, cflag                                       
         else 
            write(out,9100) element, eps_plas, eps_crit, sig_mean,           
     &              sig_mises, triaxiality, cflag 
         end if                                                 
         local_count = local_count + 1                                       
      end do  ! elem_loop
c
      if( kill_count > 0 ) write(out,9220) kill_count
      if( local_count > 0 ) then
        write(out,9200); write(out,9210)
      end if                                                              
      if( load_size_control_crk_grth ) then                                    
         write (out,9010) max_d_eps_plas                                        
      else                                                                      
         write (out,*) ' '                                                      
      end if                                                                    
c
c
c              packet output section                                                
c              ---------------------
c                                                                 
      do_packets = output_packets .and. local_count .ne. 0                                                                     
      if( .not. do_packets ) return
      if( load_size_control_crk_grth ) then                                   
           write( packet_file_no ) 23, local_count, step, iter                  
      else                                                                   
           write( packet_file_no ) 22, local_count, step, iter                  
      end if                                                                 
c                                                                               
      do elem_loop = 1, num_print_list                                          
         element  = dam_print_list(elem_loop)                                   
         elem_ptr = dam_ptr(element)                                            
c                                                                               
c             check if element is a killable element and/or if it has           
c             already been killed.                                              
c                                                                               
         if( dam_ptr(element) .eq. 0 ) cycle                                   
         if( dam_state(elem_ptr) .ne. 0 ) cycle                                
c                                                                               
c             get current smcs parameters for element and print.                
c                                                                               
         call dam_param_3_get_values( element, debug, eps_plas,               
     &                   eps_crit, sig_mean, sig_mises,                         
     &                   triaxiality, 2, ldummy )  
         skip_element = eps_plas < eps_plas_tol
         if( skip_element ) cycle
         if( load_size_control_crk_grth ) then                              
            d_eps_plas = eps_plas - old_plast_strain(dam_ptr(element))                           
            max_d_eps_plas = max(max_d_eps_plas, d_eps_plas)                 
            write(packet_file_no) element, eps_plas, eps_crit,               
     &              sig_mean, sig_mises, d_eps_plas                             
         else                                                                
            write(packet_file_no) element, eps_plas, eps_crit,               
     &              sig_mean, sig_mises                                         
         end if                                                               
      end do   ! elem_loop                                                                   
c                                                                               
      return                                                                    
c                                                                               
 9010 format(/,5x,'>> Maximum change in plastic strain for this step:',         
     &     f9.7,/)                                                             
 9020 format(2x,'element   eps-pls    eps-crit    ',                      
     &       'sig-mean  sig-mises  d(eps-pls)    triaxiality (T)',/,2x,                                   
     &          '-------   -------    --------    ',                        
     &       '--------  ---------  ----------    ---------------')                                        
 9030 format(2x,'element   eps-pls    eps-crit    ',                      
     &       'sig-mean  sig-mises  triaxiality (T)',/,2x,                                   
     &       '-------   -------    --------    ',                        
     &       '--------  ---------  ---------------')                                        
 9100 format(1x,i8,2(2x,f9.6),2x,f9.3,2x,f9.3,4x,f6.2,2x,a1) 
 9110 format(1x,i8,2(2x,f9.6),2x,f9.3,2x,f9.3,3x,f9.7,4x,f6.2,2x,a1) 
 9200 format(10x,"* -- Triaxiality (sig_mean/sig_mises) < 1.0")                                                
 9210 format(10x,"Triaxiality is mean over history weighted by ",
     & "plastic strain")   
 9220 format(10x,"Total elements now killed: ",i5)                                                
c                                             
      end                                                                       
