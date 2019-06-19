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
c     *                   last modified : 6/18/2019 rhd              *          
c     *                                                              *          
c     *     This routine prints out the status of SMCS elements      *          
c     *     marked as killable at the beginning of a load step       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dam_print_elem3( step, iter )                                  
      use global_data, only : out, noelem, nonode, ltmstp, stname
      use elem_extinct_data, only : dam_state, dam_print_list,                  
     &                              old_plast_strain
      use main_data, only :   output_packets, packet_file_no                      
      use damage_data, only : dam_ptr, load_size_control_crk_grth,
     &                        num_print_list, num_elements_killed,
     &                        num_kill_elem, print_status, 
     &                        print_top_list, num_top_list, 
     &                        smcs_states, smcs_stream, smcs_text,
     &                        smcs_states_intlst

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
     &    ddummy3, triaxiality, mean_zeta, ratio, mean_omega
      double precision, allocatable :: strain_ratios(:)
      double precision, parameter :: zero = 0.0d0, one = 1.0d0,
     &                               eps_plas_tol = 1.0d-06,  
     &                               ratio_tol = 0.001d0              
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
c              get current smcs parameters for element and print.
c              skip printing elements with very small plastic strain                
c                                                                               
         call dam_param_3_get_values( element, debug, eps_plas,               
     &                   eps_crit, sig_mean, sig_mises,                         
     &                   triaxiality, mean_zeta, mean_omega, 
     &                   2, ldummy )  
         skip_element = eps_plas < eps_plas_tol
         if( skip_element ) cycle
         cflag = " "
         if( triaxiality < one ) cflag = "*"                  
         if( load_size_control_crk_grth ) then                              
            d_eps_plas = eps_plas - old_plast_strain(dam_ptr(element))    
            max_d_eps_plas = max(max_d_eps_plas, d_eps_plas)  
            ratio = eps_plas / eps_crit
            if( ratio >= ratio_tol )
     &       write(out,9110) element, eps_plas, eps_crit, 
     &              eps_plas / eps_crit, sig_mean,           
     &              sig_mises, d_eps_plas, triaxiality, mean_zeta,
     &              mean_omega, cflag                                       
         else 
            ratio = eps_plas / eps_crit
            if( ratio >= ratio_tol )
     &         write(out,9100) element, eps_plas, eps_crit, 
     &              ratio, sig_mean,           
     &              sig_mises, triaxiality, mean_zeta, mean_omega,
     &              cflag 
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
c              output an smcs states file for step if requested
c
      if( smcs_states ) call dam_print_elem3_states_file
c
c
c              packet output section. output results for all elements
c              not yet killed.                                                
c              ------------------------------------------------------
c  
      if( .not. output_packets ) then
         call dam_print_elem3_rel
         return
      else
         call dam_print_elem3_packets
         write(out,9300) list_count                                                               
         call dam_print_elem3_rel
      end if
c
c              flat file of state variables for SMCS
c
c
      return                                                                    
c                                                                               
 9010 format(/,5x,'>> Maximum change in plastic strain for this step:',         
     &     f10.7,/)                                                             
 9020 format(2x,'element   eps-pls    eps-crit  ratio  ',                      
     & 'sig-mean  sig-mises  d(eps-pls)    triaxiality (T)   zeta',
     & '   omega',/,2x,                                   
     &          '-------   -------    --------  -----  ',                        
     & '--------  ---------  ----------    ---------------   ----',
     & '   -----')                                        
 9030 format(2x,'element   eps-pls    eps-crit  ratio  ',                      
     &       'sig-mean  sig-mises  triaxiality (T)   zeta',
     &       '   omega',/,2x,                                   
     &          '-------   -------    --------  -----  ',                        
     &       '--------  ---------  ---------------   ----   -----')                                        
 9100 format(1x,i8,2(2x,f9.6),2x,f5.3,1x,f9.3,2x,f9.3,4x,f6.2,8x,
     &       f6.2,2x,f6.2,2x,a1) 
 9110 format(1x,i8,2(2x,f9.6),2x,f5.3,1x,f9.3,2x,f9.3,3x,f9.7,
     &  7x,f6.2,7x,f6.2,2x,f6.2,2x,a1) 
 9200 format(10x,"* -- Triaxiality (sig_mean/sig_mises) < 1.0")                                                
 9210 format(10x,"-- Triaxiality is mean over history weighted by ",
     & "plastic strain",/,
     & 10x,"-- Zeta is mean of normalized Lode angle weighted by ",
     & "plastic strain",/,
     & 10x,"-- Omega is mean of Lode angle parameter weighted by ",
     & "plastic strain",
     & /,10x,"-- Elements with ratio < 0.001 not shown")   
 9220 format(10x,"Total elements now killed: ",i5)       
 9300 format(10x,"-- Packet data written for element count: ",i7)                                         
c  
      contains
c     ========
c     ****************************************************************          
c     *                                                              *          
c     *               subroutine dam_print_elem3_rel                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dam_print_elem3_rel     
      implicit none
c
      if( allocated( strain_ratios ) )  deallocate( strain_ratios )
      if( allocated( element_list ) )   deallocate( element_list )
      if( allocated( sort_index_vec ) ) deallocate( sort_index_vec )
c
      write(out,*) " " 
c
      return
c
      end subroutine dam_print_elem3_rel
c
c     ****************************************************************          
c     *                                                              *          
c     *               subroutine dam_print_elem3_fill                *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine dam_print_elem3_fill     
      implicit none
      integer :: k, element
c
      do k = 1, list_count
         element  = element_list(k)
         call dam_param_3_get_values( element, debug, eps_plas,               
     &                   eps_crit, sig_mean, sig_mises,                         
     &                   triaxiality, mean_zeta, mean_omega,
     &                   2, ldummy )  
         strain_ratios(k) =  eps_plas / eps_crit
      end do  ! elem_loop
c
      return
c
      end subroutine dam_print_elem3_fill
c     ****************************************************************          
c     *                                                              *          
c     *               subroutine dam_print_elem3_packets             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine dam_print_elem3_packets
      implicit none
c
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
     &                   triaxiality, mean_zeta, mean_omega,
     &                   2, ldummy )  
         if( load_size_control_crk_grth ) then                              
            d_eps_plas = eps_plas - old_plast_strain(dam_ptr(element))                           
            max_d_eps_plas = max(max_d_eps_plas, d_eps_plas)                 
            write(packet_file_no) element, eps_plas, eps_crit,               
     &              sig_mean, sig_mises, d_eps_plas, triaxiality,
     &              mean_zeta, mean_omega                             
         else                                                                
            write(packet_file_no) element, eps_plas, eps_crit,               
     &              sig_mean, sig_mises, triaxiality, mean_zeta, 
     &              mean_omega                                         
         end if                                                               
      end do   ! elem_loop    
c
      return
      end subroutine dam_print_elem3_packets
c     ****************************************************************          
c     *                                                              *          
c     *          subroutine dam_print_elem3_states_file              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine dam_print_elem3_states_file
      implicit none
c
      integer :: header_file_number, flat_file_number, step_num,
     &           now_len, lenlst
      logical ::  get_values, output_this_step
      logical, external :: scan_entry_in_list
      double precision :: evalues(8)
      character(len=24) :: sdate_time_tmp
      character(len=20) :: form_type, access_type
      character(len=80) :: flat_name
c
c              user requested output of element states file with smcs
c              values. write file if this step in user list
c
      step_num = ltmstp  ! current step number
      lenlst = sizeof( smcs_states_intlst )
      output_this_step = scan_entry_in_list( step_num,
     &                      smcs_states_intlst, lenlst )
      if( .not. output_this_step ) return
c
c              open/write/close the states header file. just overwrites
c              an existing file
c
      open( newunit=header_file_number,file="states_header_smcs",
     &      access="sequential", form="formatted" )
      write(header_file_number,9100) 
      close(unit=header_file_number,status="keep")
c
c              make name of the smcs states results file for this
c              load(time) step
c                wem#######_stream_smcs   wem#######_text_smcs
c
      flat_name(1:3) ='wem'
      write(flat_name(4:),9200) step_num
c
      if( smcs_text ) then
         flat_name(11:) = '_text'
         now_len = len_trim( flat_name )
         flat_name(now_len+1:) = "_" // "smcs"
         access_type = 'sequential'
         form_type   = 'formatted'
      end if
c
      if( smcs_stream ) then
         flat_name(11:) = '_stream'
         now_len = len_trim( flat_name )
         flat_name(now_len+1:) = "_" // "smcs"
         access_type = 'stream'
         form_type   = 'unformatted'
      end if
c
      open( newunit=flat_file_number,file=flat_name, status='unknown',
     &     access=access_type, form=form_type )
c
      if( smcs_text ) then 
          write(flat_file_number,9000)
          write(flat_file_number,9014)
          write(flat_file_number,9020) stname
          write(flat_file_number,9030) nonode, noelem
          call fdate( sdate_time_tmp )
          write(flat_file_number,9040) sdate_time_tmp
          write(flat_file_number,9050) step_num
          write(flat_file_number,9000)
      end if
c
c              write a result record for every element. nonkillable
c              elements have a large critical plastic strain.
c   
      do element = 1, noelem
        elem_ptr = dam_ptr(element)    
        if( elem_ptr == 0 ) then ! element not killable
          get_values = .true.
        elseif( dam_state(elem_ptr) .ne. 0 ) then !already killed
           evalues = zero
           get_values = .false.
        else ! killable elem not yet killed
           get_values = .true.
        end if
        if( get_values ) then
           call dam_param_3_get_values( element, debug, eps_plas,               
     &           eps_crit, sig_mean, sig_mises, triaxiality, 
     &           mean_zeta, mean_omega, 3, ldummy )  
           evalues(1) = eps_plas
           evalues(2) = eps_crit
           evalues(3) = eps_plas / eps_crit
           evalues(4) = sig_mean           
           evalues(5) = sig_mises
           evalues(6) = triaxiality
           evalues(7) = mean_zeta
           evalues(8) = mean_omega
        end if  
        where( abs( evalues ) .lt. eps_plas_tol ) evalues = zero
       if( smcs_stream ) write(flat_file_number) evalues
       if( smcs_text )   write(flat_file_number,9110) evalues
c
      end do ! over all elements
c
      close(unit=flat_file_number,status="keep")

      return
c
 9000 format('#')
 9014 format('#  WARP3D element results: states-smcs')
 9020 format('#  Structure name: ',a8 )
 9030 format('#  Model nodes, elements: ',2i8)
 9040 format('#  ',a24)
 9050 format('#  Load(time) step: ',i8 )
 9060 format(2i8)
 9100 format(
     &"#",/,
     &"#  header file for smcs state output",/,
     &"#",/,
     &"#  8 character state labels and longer descriptors",/,
     &"#  material model number, number of state variables",/,
     &"#",/,
     &"     8",/,
     &"      1  epspls   Eq. plastic strain",/,
     &"      2  epscrit  Critical strain",/,
     &"      3  ratio    epspls/epscrit",/,
     &"      4  mean     mean stress",/,
     &"      5  mises    mean mises",/,
     &"      6  triax    mean triaxiality",/,
     &"      7  xi       mean xi",/,
     &"      8  omega    mean omega" )
 9110 format(8e15.6)
 9200 format( i7.7 )                              
c
      end subroutine dam_print_elem3_states_file
c
      end subroutine dam_print_elem3
