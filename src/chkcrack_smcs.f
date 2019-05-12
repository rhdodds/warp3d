c ********************************************************************          
c *                                                                  *          
c *  routines to support crack growth by smcs criterion              *          
c *                                                                  *          
c ********************************************************************          
c
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_print_elem3              *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 5/12/2019 rhd              *          
c     *                                                              *          
c     *     This routine prints out the status of SMCS elements      *          
c     *     marked as killable at the beginning of a load step       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dam_print_elem3( step, iter )                                  
      use global_data, only : out, noelem  ! old common.main
      use elem_extinct_data, only : dam_state, dam_print_list,                  
     &                              old_plast_strain
      use main_data, only :   output_packets, packet_file_no                      
      use damage_data, only : dam_ptr, load_size_control_crk_grth,
     &                        num_print_list, num_elements_killed,
     &                        num_kill_elem, print_status, 
     &                        print_top_list, num_top_list 
c                                                        
      implicit none                                                    
c
c              parameters 
c          
      integer :: step, iter
c
c              local declarations
c                                                                    
      integer :: elem_loop, elem_ptr, element, local_count,
     &           list_count, list_length
      integer, allocatable :: element_list(:), sort_index_vec(:)
      double precision :: dummy, eps_plas, eps_crit, sig_mean,
     &    sig_mises, d_eps_plas, max_d_eps_plas, ddummy1, ddummy2,
     &    ddummy3, triaxiality, mean_zeta, ratio
      double precision, allocatable :: strain_ratios(:)
      double precision, parameter :: zero = 0.0d0,
     &                               eps_plas_tol = 1.0d-06                 
      logical :: ldummy, debug, skip_element
      character(len=1) :: cflag                                          
c                                                                               
      write(out,*) ' '                                                          
      write(out,*) ' ********************************************* '            
      write(out,*) ' ***         Killable element status       *** '            
      write(out,*) ' ********************************************* '            
      write(out,*)                                                              
      if( load_size_control_crk_grth ) write(out,9020)                                    
      if( .not. load_size_control_crk_grth ) write(out,9030)                                    
c
c              Support for printing elements with largest ratio
c              of plastic strain / critical strain. build sorted
c              list of elements.              
c          
      if( print_top_list ) then                                                          
        allocate( strain_ratios(num_kill_elem),
     &            element_list(num_kill_elem),
     &            sort_index_vec(num_kill_elem) ) 
        strain_ratios  = zero
        sort_index_vec = 0
        list_count     = 0
        do element = 1, noelem
          elem_ptr = dam_ptr(element) 
          if( elem_ptr == 0 ) cycle  ! element not killable
          if( dam_state(elem_ptr) .ne. 0 ) cycle ! already killed
          list_count = list_count + 1
          element_list(list_count) = element
        end do
        call dam_print_elem3_fill
        call  hpsort_eps_epw( list_count, strain_ratios, 
     &                        sort_index_vec, 1.0d-07 )
      end if
c                                                                             
c              print element status for smcs. follow user-specified
c              print order or the sorted list by descending
c              critical strain ratio
c                                                                               
      max_d_eps_plas = zero                                                     
      local_count = 0
      if( print_status )   list_length = num_print_list
      if( print_top_list ) list_length = min( num_top_list,
     &                                        list_count )
c                                                          
      do elem_loop = 1, list_length
         if( print_top_list ) then
            element = element_list(sort_index_vec(elem_loop))
         end if
         if( print_status ) then
            element = dam_print_list(elem_loop)
            elem_ptr = dam_ptr(element)    
            if( elem_ptr == 0 ) cycle   ! not killable 
            if( dam_state(elem_ptr) .ne. 0 ) cycle ! already killed
         end if                                
c
c              get current smcs parameters for element and print.
c              skip printing elements with very small plastic strain                
c                                                                               
         call dam_param_3_get_values( element, debug, eps_plas,               
     &                   eps_crit, sig_mean, sig_mises,                         
     &                   triaxiality, mean_zeta, 2, ldummy )  
         skip_element = eps_plas < eps_plas_tol
         if( skip_element ) cycle
         cflag = " "
         if( triaxiality < 1.0d0 ) cflag = "*"                  
         if( load_size_control_crk_grth ) then                              
            d_eps_plas = eps_plas - old_plast_strain(dam_ptr(element))    
            max_d_eps_plas = max(max_d_eps_plas, d_eps_plas)  
            ratio = eps_plas / eps_crit
            if( ratio >= 0.001 )
     &       write(out,9110) element, eps_plas, eps_crit, 
     &              eps_plas / eps_crit, sig_mean,           
     &              sig_mises, d_eps_plas, triaxiality, mean_zeta,
     &              cflag                                       
         else 
            ratio = eps_plas / eps_crit
            if( ratio >= 0.001 )
     &         write(out,9100) element, eps_plas, eps_crit, 
     &              ratio, sig_mean,           
     &              sig_mises, triaxiality, mean_zeta, cflag 
         end if                                                 
         local_count = local_count + 1                                       
      end do  ! elem_loop
c
      write(out,9220) num_elements_killed
      if( local_count > 0 ) then
        write(out,9200); write(out,9210)
      end if                                                              
      if( load_size_control_crk_grth ) write(out,9010) max_d_eps_plas                                        
c
c
c              packet output section. output results for all elements
c              not yet killed.                                                
c              ------------------------------------------------------
c  
      if( .not. output_packets ) then
         call dam_print_elem3_rel
         return
      end if
      if( .not. allocated( element_list ) )
     &    allocate( element_list(num_kill_elem) ) 
      list_count     = 0
      do element = 1, noelem
        elem_ptr = dam_ptr(element) 
        if( elem_ptr == 0 ) cycle  ! element not killable
        if( dam_state(elem_ptr) .ne. 0 ) cycle ! already killed
        list_count = list_count + 1
        element_list(list_count) = element
      end do
      if( list_count == 0 ) then  ! all elements killed
         call dam_print_elem3_rel
         return
      end if 
c                        
      if( load_size_control_crk_grth )                                   
     &      write( packet_file_no ) 23, list_count, step, iter                  
      if( .not. load_size_control_crk_grth )                                   
     &        write( packet_file_no ) 22, list_count, step, iter                  
c                                                                               
      do elem_loop = 1, list_count                                          
         element  = element_list(elem_loop)    
         call dam_param_3_get_values( element, debug, eps_plas,               
     &                   eps_crit, sig_mean, sig_mises,                         
     &                   triaxiality, mean_zeta, 2, ldummy )  
         if( load_size_control_crk_grth ) then                              
            d_eps_plas = eps_plas - old_plast_strain(dam_ptr(element))                           
            max_d_eps_plas = max(max_d_eps_plas, d_eps_plas)                 
            write(packet_file_no) element, eps_plas, eps_crit,               
     &              sig_mean, sig_mises, d_eps_plas, triaxiality,
     &              mean_zeta                             
         else                                                                
            write(packet_file_no) element, eps_plas, eps_crit,               
     &              sig_mean, sig_mises, triaxiality, mean_zeta                                         
         end if                                                               
      end do   ! elem_loop    
c
      write(out,9300) list_count                                                               
c        
      call dam_print_elem3_rel
c                                                                       
      return                                                                    
c                                                                               
 9010 format(/,5x,'>> Maximum change in plastic strain for this step:',         
     &     f10.7,/)                                                             
 9020 format(2x,'element   eps-pls    eps-crit  ratio  ',                      
     & 'sig-mean  sig-mises  d(eps-pls)    triaxiality (T)   zeta',/,2x,                                   
     &          '-------   -------    --------  -----  ',                        
     & '--------  ---------  ----------    ---------------   ----')                                        
 9030 format(2x,'element   eps-pls    eps-crit  ratio  ',                      
     &       'sig-mean  sig-mises  triaxiality (T)   zeta',/,2x,                                   
     &          '-------   -------    --------  -----  ',                        
     &       '--------  ---------  ---------------   ----')                                        
 9100 format(1x,i8,2(2x,f9.6),2x,f5.3,1x,f9.3,2x,f9.3,4x,f6.2,8x,
     &       f6.2,2x,a1) 
 9110 format(1x,i8,2(2x,f9.6),2x,f5.3,1x,f9.3,2x,f9.3,3x,f9.7,
     &  7x,f6.2,7x,f6.2,2x,a1) 
 9200 format(10x,"* -- Triaxiality (sig_mean/sig_mises) < 1.0")                                                
 9210 format(10x,"-- Triaxiality is mean over history weighted by ",
     & "plastic strain",/,
     & 10x,"-- Zeta is mean of normalized Lode angle weighted by ",
     & "plastic strain",
     & /,10x,"-- Elements with ratio < 0.001 not shown")   
 9220 format(10x,"Total elements now killed: ",i5)       
 9300 format(10x,"-- Packet data written for element count: ",i7)                                         
c  
      contains
c     ========
      subroutine dam_print_elem3_rel     
      implicit none
c
      if( allocated( strain_ratios ) ) deallocate( strain_ratios )
      if( allocated( element_list ) ) deallocate( element_list )
      if( allocated( sort_index_vec ) ) deallocate( sort_index_vec )
c
      write(out,*) " " 
c
      return
c
      end subroutine dam_print_elem3_rel
c
      subroutine dam_print_elem3_fill     
      implicit none
      integer :: k, element
c
      do k = 1, list_count
         element  = element_list(k)
         call dam_param_3_get_values( element, debug, eps_plas,               
     &                   eps_crit, sig_mean, sig_mises,                         
     &                   triaxiality, mean_zeta, 2, ldummy )  
         strain_ratios(k) =  eps_plas / eps_crit
      end do  ! elem_loop
c
      return
c
      end subroutine dam_print_elem3_fill
      end subroutine dam_print_elem3
                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_print_elem3              *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 5/11/2019 rhd              *          
c     *                                                              *          
c     *     This routine prints out the status of SMCS elements      *          
c     *     marked as killable at the beginning of a load step       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dam_print_elem3x( step, iter )                                  
      use global_data, only : out, noelem  ! old common.main
      use elem_extinct_data, only : dam_state, dam_print_list,                  
     &                              old_plast_strain
      use main_data, only :   output_packets, packet_file_no                      
      use damage_data, only : dam_ptr, load_size_control_crk_grth,
     &                        num_print_list, num_elements_killed,
     &                        num_kill_elem, print_status, 
     &                        print_top_list, num_top_list 
c                                                        
      implicit none                                                    
c
c              parameters 
c          
      integer :: step, iter
c
c              local declarations
c                                                                    
      integer :: elem_loop, elem_ptr, element, local_count,
     &           list_count, list_length
      integer, allocatable :: element_list(:), sort_index_vec(:)
      double precision :: dummy, eps_plas, eps_crit, sig_mean,
     &    sig_mises, d_eps_plas, max_d_eps_plas, ddummy1, ddummy2,
     &    ddummy3, triaxiality, mean_zeta, ratio
      double precision, allocatable :: strain_ratios(:)
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
      if( load_size_control_crk_grth ) write(out,9020)                                    
      if( .not. load_size_control_crk_grth ) write(out,9030)                                    
c
c              Support for printing elements with largest ratio
c              of plastic strain / critical strain. build sorted
c              list of elements.              
c          
      if( print_top_list ) then                                                          
        allocate( strain_ratios(num_kill_elem),
     &            element_list(num_kill_elem),
     &            sort_index_vec(num_kill_elem) ) 
        strain_ratios  = zero
        sort_index_vec = 0
        list_count     = 0
        do element = 1, noelem
          elem_ptr = dam_ptr(element) 
          if( elem_ptr == 0 ) cycle  ! element not killable
          if( dam_state(elem_ptr) .ne. 0 ) cycle ! already killed
          list_count = list_count + 1
          element_list(list_count) = element
        end do
        call dam_print_elem3_fill
        call  hpsort_eps_epw( list_count, strain_ratios, 
     &                        sort_index_vec, 1.0d-07 )
      end if
c                                                                             
c              print element status for smcs. follow user-specified
c              print order or the sorted list by descending
c              critical strain ratio
c                                                                               
      max_d_eps_plas = zero                                                     
      local_count = 0
      if( print_status )   list_length = num_print_list
      if( print_top_list ) list_length = min( num_top_list,
     &                                        list_count )
c                                                          
      do elem_loop = 1, list_length
         if( print_top_list ) then
            element = element_list(sort_index_vec(elem_loop))
         end if
         if( print_status ) then
            element = dam_print_list(elem_loop)
            elem_ptr = dam_ptr(element)    
            if( elem_ptr == 0 ) cycle   ! not killable 
            if( dam_state(elem_ptr) .ne. 0 ) cycle ! already killed
         end if                                
c
c              get current smcs parameters for element and print.
c              skip printing elements with very small plastic strain                
c                                                                               
         call dam_param_3_get_values( element, debug, eps_plas,               
     &                   eps_crit, sig_mean, sig_mises,                         
     &                   triaxiality, mean_zeta, 2, ldummy )  
         skip_element = eps_plas < eps_plas_tol
         if( skip_element ) cycle
         cflag = " "
         if( triaxiality < 1.0d0 ) cflag = "*"                  
         if( load_size_control_crk_grth ) then                              
            d_eps_plas = eps_plas - old_plast_strain(dam_ptr(element))    
            max_d_eps_plas = max(max_d_eps_plas, d_eps_plas)  
            ratio = eps_plas / eps_crit
            if( ratio >= 0.001 )
     &       write(out,9110) element, eps_plas, eps_crit, 
     &              eps_plas / eps_crit, sig_mean,           
     &              sig_mises, d_eps_plas, triaxiality, mean_zeta,
     &              cflag                                       
         else 
            ratio = eps_plas / eps_crit
            if( ratio >= 0.001 )
     &         write(out,9100) element, eps_plas, eps_crit, 
     &              ratio, sig_mean,           
     &              sig_mises, triaxiality, mean_zeta, cflag 
         end if                                                 
         local_count = local_count + 1                                       
      end do  ! elem_loop
c
      write(out,9220) num_elements_killed
      if( local_count > 0 ) then
        write(out,9200); write(out,9210)
      end if                                                              
      if( load_size_control_crk_grth ) write (out,9010) max_d_eps_plas                                        
      if( .not. load_size_control_crk_grth ) write (out,*) ' '                                                      
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
     &                   triaxiality, mean_zeta, 2, ldummy )  
         skip_element = eps_plas < eps_plas_tol
         if( skip_element ) cycle
         if( load_size_control_crk_grth ) then                              
            d_eps_plas = eps_plas - old_plast_strain(dam_ptr(element))                           
            max_d_eps_plas = max(max_d_eps_plas, d_eps_plas)                 
            write(packet_file_no) element, eps_plas, eps_crit,               
     &              sig_mean, sig_mises, d_eps_plas, triaxiality,
     &              mean_zeta                             
         else                                                                
            write(packet_file_no) element, eps_plas, eps_crit,               
     &              sig_mean, sig_mises, triaxiality, mean_zeta                                         
         end if                                                               
      end do   ! elem_loop                                                                   
c                                                                               
      return                                                                    
c                                                                               
 9010 format(/,5x,'>> Maximum change in plastic strain for this step:',         
     &     f10.7,/)                                                             
 9020 format(2x,'element   eps-pls    eps-crit  ratio  ',                      
     & 'sig-mean  sig-mises  d(eps-pls)    triaxiality (T)   zeta',/,2x,                                   
     &          '-------   -------    --------  -----  ',                        
     & '--------  ---------  ----------    ---------------   ----')                                        
 9030 format(2x,'element   eps-pls    eps-crit  ratio  ',                      
     &       'sig-mean  sig-mises  triaxiality (T)   zeta',/,2x,                                   
     &          '-------   -------    --------  -----  ',                        
     &       '--------  ---------  ---------------   ----')                                        
 9100 format(1x,i8,2(2x,f9.6),2x,f5.3,1x,f9.3,2x,f9.3,4x,f6.2,8x,
     &       f6.2,2x,a1) 
 9110 format(1x,i8,2(2x,f9.6),2x,f5.3,1x,f9.3,2x,f9.3,3x,f9.7,
     &  7x,f6.2,7x,f6.2,2x,a1) 
 9200 format(10x,"* -- Triaxiality (sig_mean/sig_mises) < 1.0")                                                
 9210 format(10x,"-- Triaxiality is mean over history weighted by ",
     & "plastic strain",/,
     & 10x,"-- Zeta is mean of normalized Lode angle weighted by ",
     & "plastic strain",
     & /,10x,"-- Elements with ratio < 0.001 not shown")   
 9220 format(10x,"Total elements now killed: ",i5)                                                
c  
      contains
c     ========
c                                       
      subroutine dam_print_elem3_fill     
      implicit none
      integer :: k, element
c
      do k = 1, list_count
         element  = element_list(k)
         call dam_param_3_get_values( element, debug, eps_plas,               
     &                   eps_crit, sig_mean, sig_mises,                         
     &                   triaxiality, mean_zeta, 2, ldummy )  
         strain_ratios(k) =  eps_plas / eps_crit
      end do  ! elem_loop
c
      return
c
      end subroutine dam_print_elem3_fill
      end subroutine dam_print_elem3x
