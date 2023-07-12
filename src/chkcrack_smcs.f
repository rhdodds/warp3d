!c ********************************************************************          
c *                                                                  *          
c *  routines to support crack growth by smcs criterion              *          
c *                                                                  *          
c ********************************************************************          
c
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_print_smcs               *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 7/11/23 rhd                *          
c     *                                                              *          
c     *     Prints out the status of SMCS elements                   *          
c     *     marked as killable at the beginning of a load step       *  
c     *     drives writing binary packets and states file results    *        
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dam_print_smcs( step, iter )    
c
      use global_data, only : out, noelem, nonode, ltmstp, stname
      use elem_extinct_data, only : dam_state, dam_print_list,                  
     &                              old_plast_strain
      use main_data, only :   output_packets, packet_file_no                      
      use damage_data, only : dam_ptr, load_size_control_crk_grth,
     &                        num_print_list, num_elements_killed,
     &                        num_kill_elem, print_status, 
     &                        print_top_list, num_top_list, 
     &                        smcs_states, smcs_stream, smcs_text,
     &                        smcs_states_intlst, smcs_type,
     &                        use_weighted
      use constants
c                                                        
      implicit none                                                    
c
c              parameters 
c          
      integer, intent(in) :: step, iter
c
c              local declarations
c                                                                    
      integer :: i, elem_loop, elem_ptr, element, local_count,
     &           list_count, list_length, dowhat 
      integer, allocatable :: element_list(:), sort_index_vec(:)
      double precision :: eps_plas, d_eps_plas, max_d_eps_plas, sig_1
      double precision, allocatable :: damage_factors(:)
      double precision, parameter :: damage_tol = 0.05d00
      logical :: debug, skip_element, get_princ, write_packets
      logical, save :: print_info = .true.
      character(len=1) :: cflag      
c
      include 'include_damage_values'
      type(values_T) :: values
c
      debug = .false.
c
      write(out,*) ' '                                                          
      write(out,*) ' ********************************************* '            
      write(out,*) ' ***         Killable element status       *** '            
      write(out,*) ' ********************************************* '            
      write(out,*)           
      if( load_size_control_crk_grth ) write(out,9020)                                    
      if( .not. load_size_control_crk_grth ) write(out,9030)                                    
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
        dowhat = 1  
        get_princ = .false.  ! princ stresses not needed
c
        do element = 1, noelem
          elem_ptr = dam_ptr(element) 
          if( elem_ptr == 0 ) cycle  ! element not killable
          if( dam_state(elem_ptr) .ne. 0 ) cycle ! already killed
          list_count = list_count + 1
          element_list(list_count) = element
        end do  
c$OMP PARALLEL DO PRIVATE( i, element, elem_ptr, values )
        do i = 1, list_count
          element = element_list(i)
          call dam_param_smcs_get_values( element, dowhat, values,
     &                                    get_princ )                                                                   
          damage_factors(i) = values%damage_factor
        end do
c$OMP END PARALLEL DO
        call hpsort_eps_epw( list_count, damage_factors, 
     &                        sort_index_vec, 1.0d-07 )
      end if ! print_top_list
c                                                                             
c              print element status for smcs. follow user-specified
c              print order or the sorted list by descending
c              damage factor
c                                                                               
      max_d_eps_plas = zero                                                     
      local_count = 0
      if( print_status )   list_length = num_print_list
      if( print_top_list ) list_length = min( num_top_list,
     &                                        list_count )
c     
      associate( v => values )
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
c              skip printing elements with very small damage. dowhat = 1
c              means don't update anything. just return values.
c      
         dowhat = 1 
         get_princ = .false.  ! no sig1,.. needed     
         call dam_param_smcs_get_values( element, dowhat, v, get_princ )                                                                   
         skip_element = v%damage_factor < damage_tol
         if( skip_element ) cycle
         eps_plas = v%eps_plas
         cflag = " "
         if( v%triaxiality < one ) cflag = "*"                  
         if( load_size_control_crk_grth ) then                              
            d_eps_plas = eps_plas - old_plast_strain(dam_ptr(element))    
            max_d_eps_plas = max(max_d_eps_plas, d_eps_plas)  
            write(out,9110) element, eps_plas, v%eps_critical, 
     &              v%damage_factor, v%sig_mean, v%sig_mises, 
     &              d_eps_plas, v%triaxiality, v%zeta,
     &              v%omega, v%bar_theta, v%tearing_param, cflag                                       
         else 
            write(out,9100) element, eps_plas, v%eps_critical, 
     &              v%damage_factor, v%sig_mean, v%sig_mises,
     &              v%triaxiality, v%zeta, v%omega, v%bar_theta,
     &              v%tearing_param, cflag    
         end if                                                 
         local_count = local_count + 1                                       
c
      end do  ! elem_loop
c
      end associate 
c      
      if( load_size_control_crk_grth ) write(out,9010) max_d_eps_plas 
      write(out,9220) num_elements_killed
      if( local_count > 0 ) then
        if( print_info ) then
           write(out,9200)
           if( use_weighted ) write(out,9205)
           if( .not. use_weighted ) write(out,9206)
           write(out,9210)
         end if   
        if( .not. print_info ) write(out,9202)
        print_info = .false.
        if( smcs_type == 5 ) write(out,9215)
      end if                                                              
c
c              output an smcs states file for step if requested
c
      if( smcs_states ) call dam_print_smcs_states
c
c              packet output section. output results for all elements
c              not yet killed.                                                
c  
      write_packets = local_count .ne. 0 .and. output_packets 
      if( write_packets ) then
         call dam_print_smcs_packets
         write(out,9300) local_count
      end if   
c
c              flat file of state variables for SMCS
c
      return                                                                    
c                                                                               
 9010 format(/,10x,'Max delta plastic strain over step ',
     &  '(listed elements):',f10.7)                                                             
 9020 format(2x,'element   eps-pls    eps-crit  damage  ',                      
     & 'sig-mean  sig-mises  d(eps-pls)     T       xi',
     & '       omega      theta    tear_param',/,2x, 
     &          '-------   -------    --------  ------  ',                        
     & '--------  ---------  ----------   -----   ------',
     & '   ---------  ---------  ----------')                                        
 9030 format(2x,'element   eps-pls    eps-crit   damage  ',                      
     & 'sig-mean  sig-mises    T       xi',
     & '       omega       theta   tear_param',/,2x,                                   
     & '-------   -------    --------  -------  ',                        
     & '--------  ---------  -----   ------   ---------  ',
     & '  -------  ----------')
 9100 format(1x,i8,2(2x,f9.6),4x,f5.3,1x,f9.3,2x,f9.3,1x,f6.3,3x,
     &       f6.3,6x,f6.3,5x,f6.3,5x,f7.3,2x,a1) 
 9110 format(1x,i8,2(2x,f9.6),2x,f5.3,1x,f9.3,2x,f9.3,3x,f9.7,
     &  3x,f6.3,3x,f6.3,6x,f6.3,5x,f6.3,5x,f7.3,2x,a1) 
 9200 format(8x,"* -- T (sig_mean/sig_mises) < 1.0")                                                
 9202 format(10x,"** See earlier messages w/ first element values ",
     & "about T, xi, ... definitions")                                                
 9205 format(10x,
     & "-- T, xi, omega, theta are *plastic weighted* average values")
 9206 format(10x,
     & "-- T, xi, omega, theta are *current* values")
 9210 format(
     & 10x,"-- T is element average triaxiality",/,
     & 10x,"-- xi is element average of normalized Lode angle",/,
     & 10x,"-- omega is element average of Lode angle parameter",/,
     & 10x,"-- theta is element average of Lode angle parameter",
     &       " - type 4 MMC",/,
     & 10x,"-- tear_param is element average Sandia/Wellman tearing",
     &     " parameter- type 5",
     & /,10x,"-- Elements with small damage factors not shown",/)  
 9215 format(10x,"** for Type 5 (tearing parameter), eps-crit set",
     & ' to large value for printing purposes. not used',
     &   ' in computations',/)
 9220 format(/,10x,"Total elements now killed: ",i7,/)       
 9300 format(10x,"-- Packet data written for element count: ",i7)                                         
c  
      contains
c
c     ****************************************************************          
c     *                                                              *          
c     *               subroutine dam_print_smcs_packets              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine dam_print_smcs_packets
c      
      implicit none
c
      integer :: packet_type 
      logical :: write_debug_text
c
      write_debug_text = .false.
c
c              packet written only for elements that appeared
c              int the printed list just generated
c
      if( print_status )   list_length = num_print_list
      if( print_top_list ) list_length = min( num_top_list,
     &                                        list_count )
c
      packet_type = 22                                                                 
      if( load_size_control_crk_grth ) packet_type = 23 ! adaptive P control
      write( packet_file_no ) packet_type, local_count, step, iter    
      if( write_debug_text ) write(out,9010 ) packet_type, local_count,
     &                        step, iter               
c
      associate( v => values )
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
         dowhat = 1
         get_princ = .false.                           
         call dam_param_smcs_get_values( element, dowhat, values,
     &                                   get_princ )
         skip_element = v%damage_factor < damage_tol
         if( skip_element ) cycle
         if( load_size_control_crk_grth ) then                              
            d_eps_plas = v%eps_plas -
     &                   old_plast_strain(dam_ptr(element))                           
            write(packet_file_no) element, v%eps_plas, v%eps_critical,               
     &              v%sig_mean, v%sig_mises, d_eps_plas, v%triaxiality,
     &              v%zeta, v%omega, v%bar_theta, v%tearing_param                            
            if( write_debug_text ) 
     &         write(out,9020) element, v%eps_plas, v%eps_critical,               
     &            v%sig_mean, v%sig_mises, d_eps_plas, v%triaxiality,
     &            v%zeta, v%omega, v%bar_theta, v%tearing_param                                                      
         else                                                                
            write(packet_file_no) element, v%eps_plas, v%eps_critical,               
     &              v%sig_mean, v%sig_mises, v%triaxiality,
     &              v%zeta, v%omega, v%bar_theta, v%tearing_param                            
            if( write_debug_text ) 
     &         write(out,9030) element, v%eps_plas, v%eps_critical,               
     &            v%sig_mean, v%sig_mises, v%triaxiality,
     &            v%zeta, v%omega, v%bar_theta, v%tearing_param                                                      
                                            
         end if                                                               
      end do   ! elem_loop    
c
      end associate 
c
      return
c
 9010 format(/,".... binary packet file:",
     &  /,10x,"type, line count, step, iter: ",i2,i4,i8,i4)
 9020 format(10x,"element, eps_plas, eps_crit, sig_mean: ",
     &           i7,f8.5,f8.5,f10.4,
     &     /,15x,"mises, deps_plas, triaxiality: ",
     &       f7.4,f8.6,f8.3,
     &     /15x,"zeta, omega, bar_theta, TP: ",
     &       f7.4,f7.4,f7.4,f7.4 ) 
9030  format(10x,"element, eps_plas, eps_crit, sig_mean: ",
     &           i7,f8.5,f8.5,f10.4,
     &     /,15x,"mises, triaxiality: ",
     &       f7.4,f8.3,
     &     /15x,"zeta, omega, bar_theta, TP: ",
     &       f8.4,f9.4,f9.4,f9.4 ) 
c      
      end subroutine dam_print_smcs_packets
c     ****************************************************************          
c     *                                                              *          
c     *          subroutine dam_print_smcs_states                    *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine dam_print_smcs_states  
c      
      implicit none
c
      integer :: header_file_number, flat_file_number, step_num,
     &           now_len, lenlst
      logical ::  get_values, output_this_step, ldummy
      logical, external :: scan_entry_in_list
      double precision :: evalues(10)
      character(len=24) :: sdate_time_tmp
      character(len=20) :: form_type, access_type
      character(len=80) :: flat_name
c
c              user requested output of element states file with smcs
c              values. write file if this step in user list
c
      step_num = ltmstp  ! current step number
      lenlst = sizeof( smcs_states_intlst ) / 4
      output_this_step = scan_entry_in_list( step_num,
     &                        smcs_states_intlst, lenlst )
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
c              write a result record for every element. non-killable
c              elements have all zero values
c     
      associate( v => values )
c  
      do element = 1, noelem
        elem_ptr = dam_ptr(element)    
        if( elem_ptr == 0 ) then ! element not killable
          get_values = .false.
        elseif( dam_state(elem_ptr) .ne. 0 ) then !already killed
           get_values = .false.
        else ! killable elem not yet killed
           get_values = .true.
        end if
        evalues = zero
        if( get_values ) then
           dowhat = 1
           get_princ = .false.
           call dam_param_smcs_get_values( element, dowhat, values,
     &                                     get_princ )                                                                   
           evalues(1) = v%eps_plas
           evalues(2) = v%eps_critical
           evalues(3) = v%damage_factor
           evalues(4) = v%sig_mean           
           evalues(5) = v%sig_mises
           evalues(6) = v%triaxiality
           evalues(7) = v%zeta
           evalues(8) = v%omega
           evalues(9) = v%bar_theta
           evalues(10) = v%tearing_param
        end if  
       if( smcs_stream ) write(flat_file_number) evalues
       if( smcs_text )   write(flat_file_number,9110) evalues
c
      end do ! over all elements
c      
      end associate 
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
     &"     10",/,
     &"      1  epspls   Eq. plastic strain",/,
     &"      2  epscrit  Critical strain",/,
     &"      3  damage   0 <= damage <=1",/,
     &"      4  mean     mean stress",/,
     &"      5  mises    mean mises",/,
     &"      6  triax    Triaxiality",/,
     &"      7  xi       xi",/,
     &"      8  omega    omega",/,
     &"      9  theta    hat theta",/
     &"     10  tear_p   tearing_parm")
 9110 format(10e15.6)
 9200 format( i7.7 )                              
c
      end subroutine dam_print_smcs_states
c
      end subroutine dam_print_smcs
c
c     ****************************************************************          
c     *                                                              *          
c     *           subroutine dam_param_smcs_get_values               *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 5/18/23 rhd                *          
c     *                                                              *          
c     *             SMCS compute values for this element w/ optional *
c     *             update of values                                 *
c     *                                                              *          
c     ****************************************************************          
c
                                                                     
      subroutine dam_param_smcs_get_values( elem, dowhat, values,
     &                                      get_princ )
c              
      use global_data, only : iprops, nstr, nstrs, out,
     &                        now_step => current_load_time_step
      use main_data,       only : elems_to_blocks                               
      use elem_block_data, only : history_blocks, eps_n_blocks,                 
     &                            urcs_n_blocks, history_blk_list               
      use damage_data, only     : dam_ptr, smcs_gamma, smcs_alpha, 
     &                            smcs_beta,  smcs_type, 
     &                            smcs_beta_1,smcs_beta_2, 
     &                            smcs_a_plus, smcs_a_minus,
     &                            smcs_kappa, smcs_alpha_1, 
     &                            smcs_alpha_2, max_eps_critical,
     &                            smcs_cutoff_triaxiality,
     &                            smcs_type_4_A, smcs_type_4_n,
     &                            smcs_type_4_c1, smcs_type_4_c2,
     &                            smcs_type_4_c3, smcs_type_5_power,
     &                            smcs_type_5_tp_critical,
     &                            num_kill_elem, use_weighted
      use elem_extinct_data, only : smcs_old_epsplas, 
     &                              smcs_tear_param,
     &                              smcs_weighted_T,
     &                              smcs_weighted_zeta,
     &                              smcs_weighted_bar_theta
c
      use constants
c                                                                               
c        elem      -- (input)   element number to be checked for                
c        dowhat    -- (input)   = 1 return damage variables - no update.
c                                    set kill now flag
c                                = 2 just return values - no update 
c                                = 3 no longer used
c                                = 4 compute/update damage variables  
c        values    -- (output) derived type to store values
c        get_princ -- (input)  .true. always return principal stresses
c                                                                                                                                                              
      implicit none                                                    
c                                                                               
c               parameter declarations                                           
c                             
      integer, intent(in) :: elem, dowhat       
      logical, intent(in) :: get_princ                                       
c
      include 'include_damage_values'
      type(values_T), intent(out) :: values
c                                                                               
c               local declarations                                               
c
      integer :: gp, blk, epsoffset, elem_type,
     &           ngp, rel_elem, sigoffset, j, elem_ptr
      double precision, dimension(:), pointer :: urcs_n, eps_n                         
      double precision ::  eps_plas, eps_crit, sig_mean, sig_mises,
     &                     triaxiality, tear_param   
      double precision :: 
     &   sig(6), eps(6), fpngp, epspls_vec(27), sig_1, sig_2, sig_3,
     &   tp_kernel, deps_plas, s11, s22, s33, s12, s13, s23, theta_mit,
     &   j2, j3, zeta, omega, t1, t2, t3, t4, bar_theta,
     &    princ_dcosines(3,3), princ_values(3), t
      logical :: update, local_debug, threed_solid_elem, l1, l2, l3,
     &           l4, l5, l6, l7, l8, l9, l10, kill_now_local
      double precision, save, allocatable :: weighted_T(:)
         
      double precision, parameter :: 
     &     small_plastic_strain = 1.0d-07, small_mises = 1.0d-07,
     &     type_2_tol = 1.0d-06, 
     &     type_4_toler = 0.0001d0
c                              
c
      update = dowhat .eq. 4
c
      local_debug = .false. !elem == 1 .or. elem == 31
      if( local_debug )
     &    write(out,*) '... dam_param_smcs_get_values, dowhat, elem: ',
     &                     dowhat, elem
c
      elem_type   = iprops(1,elem)
      elem_ptr    = dam_ptr(elem)
      ngp         = iprops(6,elem)                                           
      blk         = elems_to_blocks(elem,1)                                  
      rel_elem    = elems_to_blocks(elem,2)                                  
      epsoffset   = (rel_elem-1)*nstr*ngp                                    
      sigoffset   = (rel_elem-1)*nstrs*ngp                                   
      urcs_n      => urcs_n_blocks(blk)%ptr                                  
      eps_n       => eps_n_blocks(blk)%ptr    
c
      eps_plas    = zero
      eps_crit    = max_eps_critical
      sig_mean    = zero
      sig_mises   = zero
      triaxiality = zero
      zeta        = zero
      omega       = zero  
      bar_theta   = zero    
      tp_kernel   = zero   
      sig         = zero  ! 6x1
      eps         = zero  ! 6x1  strains not used current code       
c
      if( now_step < 2 ) then ! no need to check later
      call set_element_type( elem_type, threed_solid_elem, l1, l2, l3,
     &                       l4, l5, l6, l7, l8, l9, l10 )  
      if( .not. threed_solid_elem ) then
          write(out,9300) elem
          call die_abort
        end if
      end if
c
c               average plastic strain over all integration points
c
      call mm_return_values( "plastic_strain", elem, epspls_vec, ngp )
c
c               average stresses over all integration points
c               stress/strain vectors are (x,y,z,xy,yz,xz).                
c
      do gp = 1, ngp    
        do j = 1, 6                                                     
         sig(j) = sig(j) + urcs_n(sigoffset+j)  
         eps(j) = eps(j) + eps_n(epsoffset+j)
        end do                                      
        eps_plas  = eps_plas + epspls_vec(gp)
        epsoffset = epsoffset + nstr                                         
        sigoffset = sigoffset + nstrs                                        
      end do 
c
c               average mean stress, mises stress, plastic strain,
c               Lode parameters, principal stresses. 
c               * only the tearing parameter criterion (#5)
c                 need the principal stresses. 
c               or calling routine can request them as well
c               This is an option since principal computation is most
c               most expensive of present calculations
c            
c
      fpngp     = one / dble( ngp )
      eps_plas  = fpngp * eps_plas     !  averaged at element center
      sig = sig * fpngp
      princ_dcosines = zero
      princ_values   = zero
      if( smcs_type == 5 .or. get_princ )
     &   call principal_values( sig, princ_dcosines, princ_values )
      sig_1 = princ_values(1)
      sig_2 = princ_values(2)
      sig_3 = princ_values(3)
      sig_mean  = (sig(1) + sig(2) + sig(3)) * third                    
      s11 = sig(1) - sig_mean
      s22 = sig(2) - sig_mean
      s33 = sig(3) - sig_mean
      s12 = sig(4)
      s13 = sig(5)
      s23 = sig(6)
c    
      j3 = s11*s22*s33 + two*s12*s23*s13 - s11*s23*s23 -
     &     s22*s13*s13 - s33*s12*s12
      j2 = half * ( s11**2 + s22**2 + s33**2 + 
     &              two*(s12**2 + s23**2 + s13**2))
      sig_mises = sqrt( three * j2 )
      zeta      = zero  !  called \xi in manual
      omega     = zero
      bar_theta = zero
      tear_param = zero
      triaxiality = zero
c
      if( sig_mises > small_mises ) then
        triaxiality = sig_mean / sig_mises
        zeta  = onept5*root3 * j3 / j2**onept5
        omega = one - zeta*zeta 
        t2 = thirteenpt5 * j3 / (three*j2)**onept5
        if( abs(t2) > one+type_4_toler ) then
          write(out,9118) elem, t2
          call die_abort
        end if
        theta_mit = third * acos(t2)
        if( theta_mit < zero .or. theta_mit > pi/three ) then
          write(out,9119) elem, theta_mit
          call die_abort
        end if
        bar_theta = one - six * theta_mit / pi
        if( abs(bar_theta) > one+type_4_toler ) then
          write(out,9120) elem, bar_theta  ! improve message
          call die_abort
        endif
      end if
      if( local_debug ) write(out,9100) sig_mean, sig_mises, eps_plas,
     &                   princ_values, zeta, omega, bar_theta,
     &                   triaxiality
c
      if( dowhat == 3 ) then  ! no longer supported
         write(out,9200)
         call die_abort
      end if
c          
      t =  smcs_old_epsplas(elem_ptr)
      deps_plas = eps_plas - t
      if( update ) then
         if( local_debug ) write(out,9102) t, deps_plas, t + deps_plas
         smcs_old_epsplas(elem_ptr) = eps_plas 
      end if
c
      if( use_weighted .and. update )
     &    call dam_param_smcs_get_values_weighted
      if( use_weighted ) then ! history avg values
         triaxiality = smcs_weighted_T(elem_ptr) / eps_plas
         zeta        = smcs_weighted_zeta(elem_ptr) / eps_plas
         omega       = one - zeta*zeta 
         bar_theta   = smcs_weighted_bar_theta(elem_ptr) / eps_plas
      end if
c
c               if no plastic strain yet, set defaults and return.
c
      values%eps_plas         = eps_plas
      values%old_eps_plas     = smcs_old_epsplas(elem_ptr)
      values%deps_plas        = deps_plas
      values%sig_mean         = sig_mean
      values%sig_mises        = sig_mises
      values%triaxiality      = triaxiality
      values%sig_princ        = princ_values
      values%max_princ_stress = sig_1
      values%zeta             = zeta
      values%omega            = omega
      values%bar_theta        = bar_theta
      values%tearing_param    = zero
      values%damage_factor    = zero
      values%eps_critical     = max_eps_critical
      values%kill_now         = .false.
c      
      if( eps_plas < small_plastic_strain ) return
c
c               use current or weighted triaxiality, zeta, omega, 
c               bar_theta values to set critical plastic strain. or
c               check tearing parameter
c
      eps_crit = max_eps_critical
c
      select case( smcs_type )
c
       case( 1 )
        if( triaxiality <= smcs_cutoff_triaxiality ) return
        t2 = exp(-smcs_kappa*abs(zeta))
        t1 = smcs_alpha * exp(-smcs_beta * triaxiality)
        eps_crit = smcs_gamma + t1*t2
        kill_now_local = ( eps_plas - eps_crit ) >= zero .or. 
     &                   eps_plas >= max_eps_critical
        values%eps_plas     = eps_plas
        values%eps_critical = eps_crit
        values%kill_now     = kill_now_local
        values%triaxiality  = triaxiality
        values%damage_factor = eps_plas / eps_crit
        if( local_debug ) write(out,9110) t1, t2, eps_crit,
     &                    values%damage_factor
c
       case( 2 ) 
c
        if( triaxiality <= smcs_cutoff_triaxiality ) return
        t3 = exp(-smcs_kappa*abs(zeta))
        t2 = smcs_beta_2 * exp(-smcs_a_minus * triaxiality)
        t1 = smcs_beta_1 * exp(smcs_a_plus * triaxiality)
        t4 = t1 - t2
        if( abs(t4) < type_2_tol .or. t4 < zero ) then
          eps_crit = max_eps_critical
        else
          eps_crit = smcs_gamma + t3 / t4
        end if
        kill_now_local = ( eps_plas - eps_crit ) >= zero .or. 
     &                   eps_plas >= max_eps_critical
        values%eps_plas     = eps_plas
        values%eps_critical = eps_crit
        values%kill_now     = kill_now_local
        values%damage_factor = eps_plas / eps_crit
        if( local_debug ) write(out,9112) t1, t2, t3, t4, eps_crit,
     &                    values%damage_factor
c
       case( 3 ) 
c
        if( triaxiality <= smcs_cutoff_triaxiality ) return
        block
          double precision :: t2, t1, p_omega
          p_omega = (one - omega)**2
          t1 = smcs_alpha_1*exp(-smcs_beta_1 * triaxiality )
          t2 = smcs_alpha_2*exp(-smcs_beta_2 * triaxiality )
          eps_crit = (one-p_omega)*t1 + p_omega*t2
          kill_now_local = ( eps_plas - eps_crit ) >= zero .or. 
     &                    eps_plas >= max_eps_critical
          values%eps_plas     = eps_plas
          values%eps_critical = eps_crit
          values%kill_now     = kill_now_local
          values%damage_factor = eps_plas / eps_crit
        end block
c
       case( 4 )
c
        if( triaxiality <= smcs_cutoff_triaxiality ) return
        block
          double precision :: A, n, c1, c2, c3, eta, t1, t2, ta, tb,
     &                        t1c, t1s, t1sec
          A     = smcs_type_4_A
          n     = smcs_type_4_n
          c1    = smcs_type_4_c1
          c2    = smcs_type_4_c2
          c3    = smcs_type_4_c3
          eta   = triaxiality
          t1    = bar_theta * pi / six
          t1c   = cos( t1 )
          t1s   = sin( t1 )
          t1sec = one / t1c
          t2 = root3 / ( two - root3 )
          ta = (A/c2) * ( c3 + t2*(one-c3)*(t1sec-one) )
          tb = sqrt( third*(one+c1*c1) ) * t1c + 
     &         c1*(eta+third*t1s)
          eps_crit = ( ta * tb )**(-one/n)
          kill_now_local = ( eps_plas - eps_crit ) >= zero .or. 
     &                     eps_plas >= max_eps_critical
          values%eps_plas     = eps_plas
          values%eps_critical = eps_crit
          values%kill_now     = kill_now_local
          values%damage_factor = eps_plas / eps_crit
        end block
c
      case( 5 ) 
c
        if( update ) then
          tp_kernel = twothird * sig_1/(sig_1-sig_mean)
          if( tp_kernel < zero ) tp_kernel = zero
          tp_kernel = tp_kernel ** smcs_type_5_power
          t = smcs_tear_param(elem_ptr)
          if( local_debug ) write(out,9310) t, deps_plas * tp_kernel,
     &         t + deps_plas * tp_kernel
          smcs_tear_param(elem_ptr) = t + deps_plas * tp_kernel
        end if
c
        tear_param = smcs_tear_param(elem_ptr)
        eps_crit = max_eps_critical
        kill_now_local = tear_param >= smcs_type_5_tp_critical .or. 
     &                     eps_plas >= max_eps_critical
        values%eps_critical  = ten ! just for printing
        values%kill_now      = kill_now_local
        values%damage_factor = tear_param / smcs_type_5_tp_critical
        values%tearing_param = tear_param
        if( local_debug ) write(out,9315) tear_param,
     &     values%damage_factor, kill_now_local
c
       case default
          write(out,9130) elem, smcs_type
          call die_abort
      end select
c      
      return
c
 9100 format(10x,'sig-mean, sig-mises, eps_pls:',2f10.3,f12.6,
     & /,10x,'princ. sigs:',3f10.3,
     & /,10x,'zeta, omega, bar_theta:', 3f12.4,
     & /,10x,'triaxiality: ',f8.3)
 9102 format(10x,'updating damage quantities:',
     & /,15x,'smcs_old_epsplas:       ',f11.8, 
     & /,15x,'deps_plas:              ',f11.8,
     & /,15x,'updated smcs_old_epsplas',f11.8 )
 9110 format(10x,'updating SMCS type 1-  t1, t2',2f10.3,
     &   /,15x,'eps_crit, damage: ',f8.5, 2x, f6.3)
 9112 format(10x,'updating SMCS type 2-  t1-t4',4f8.3,
     &   /,15x,'eps_crit, damage: ',f8.5, 2x, f6.3)
 9118 format(/1x,'>>>>> error: inconsistency in bar_theta computation',
     & ' for smcs type 4 (MIT)',/,
     & 14x,'element: ',i8,' acos argument invalid: ',f6.2,/,
     & 14x,'value must lie in range (-1.0,1.0)',/,
     & 14x,'job aborted by dam_param_smcs_get_values',/)
 9119 format(/1x,'>>>>> error: inconsistency in bar_theta computation',
     & ' for smcs type 4 (MIT)',/,
     & 14x,'element: ',i8,' theta value: ',f6.2,/,
     & 14x,'value must lie in range (0.0, pi/3)',/,
     & 14x,'job aborted by dam_param_smcs_get_values',/)
 9120 format(/1x,'>>>>> error: inconsistency in bar_theta computation',
     & ' for smcs type 4 (MIT)',/,
     & 14x,'element: ',i8,' bar_theta value: ',f6.2,/,
     & 14x,'value must lie in range (-1.0,1.0)',/,
     & 14x,'job aborted by dam_param_smcs_get_values',/)
 9130 format(/1x,'>>>>> error: inconsistency in smcs_type',/,
     & 14x,'element: ',i8,' smcs_type: ',i5,/,
     & 14x,'job aborted by dam_param_smcs_get_values',/)
 9200 format(/1x,'>>>>> FATAL ERROR: code inconsistency in routine: ',
     &   /,      '                   dam_param_smcs_get_values',
     &   /,      '                   analysis terminated'//)
 9300 format('>> FATAL ERROR: routine ',
     & 'dam_param_smcs_get_values.',
     & 10x,'inconsistent condition. for element:',i8,
     & ' job terminated.'//)
 9310 format(10x,'updating Sandia TP:',
     & /,15x,'old TP: ',f6.3,'  incr TP: ',f6.3,'  updated TP:',f6.3 )
 9315 format(10x,'TP returned values:',
     & /,15x,'TP: ',f6.3,'  damage: ',f6.3,'  kill:',l2 )
c
      contains
c     ****************************************************************          
c     *                                                              *          
c     *       subroutine dam_param_smcs_get_values_weighted          *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 4/15/23 rhd                *          
c     *                                                              *          
c     *     For use of weighted by plastic strain values, accumulate *
c     *     current values * deps_plas into global vectors           *  
c     *     drives writing binary packets and states file results    *        
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c
      subroutine dam_param_smcs_get_values_weighted
      implicit none
c
      double precision :: tp_kernel, tear, x, now_triaxiality
c
      if( sig_mises <= small_mises ) return
c
      now_triaxiality = sig_mean / sig_mises 
      if( now_triaxiality < smcs_cutoff_triaxiality ) return
c
      x = smcs_weighted_T(elem_ptr)
      x = x + deps_plas * now_triaxiality
      smcs_weighted_T(elem_ptr) = x
c
      x = smcs_weighted_zeta(elem_ptr)
      x = x + deps_plas * zeta
      smcs_weighted_zeta(elem_ptr) = x
c
      x = smcs_weighted_bar_theta(elem_ptr)     
      x = x + deps_plas * bar_theta
      smcs_weighted_bar_theta(elem_ptr) = x
c
      return  
c
      end subroutine dam_param_smcs_get_values_weighted
      end subroutine dam_param_smcs_get_values

