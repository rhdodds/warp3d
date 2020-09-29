c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine chkcrack                     *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 04/28/2019 rhd             *          
c     *                                                              *          
c     *        This routine drives the check on conditions to tell   *          
c     *        if crack growth will occur.  This calls routines for  *          
c     *        both crack growth and discrete crack growth.          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine chkcrack( step, iter )                                         
      use global_data ! old common.main
      use damage_data                                                           
      implicit none                                                    
c                                                                               
      integer :: step, iter
      logical :: debug                                                             
c                                                                               
      debug = .false.                                                           
c                                                                               
      if ( debug ) write(out,*) '>>> enter chkcrack <<<'                        
c                                                                               
c              Return if no crack growth specified.                             
c                                                                               
c              crack_growth_type:                                               
c                                                                               
c                 (1)  element extinction by critical porosity                  
c                      in Gurson-Tvergaard material model                       
c                      or the extended Gurson-Tvergaard model                   
c                 (2)  node release on attainment of critical                   
c                      ctoa                                                     
c                 (3)  element extinction by stress-modified                    
c                      critical strain criterion. any material                  
c                      model.                                                   
c                 (4)  crack growth using cohesive elements              
c                 (5)  user requested by input command. element list
c                      provided in module damage_data       
c                                                                               
c              Lower level routines are separated into other .f files by        
c              the type of crack growth.                                        
c                                                                               
c                                                                               
      select case( crack_growth_type )                                           
      case(0)                                                                  
         if(debug) write (out,*) '>>>> no crack growth specified.'             
      case(1)                                                                  
         if(debug) write (out,*) '>>>> use gurson crack growth.'               
         call chk_elem_kill( debug, step, iter )                                
      case(2)                                                                  
         if(debug) write (out,*) '>>>> use node release crack growth.'         
         call chk_node_release( debug, step, iter )                             
      case(3)                                                                  
         if(debug) write (out,*) '>>>> use smcs crack growth.'                 
         call chk_elem_kill( debug, step, iter )                                
      case(4)    !  chohesive                                                               
         call chk_elem_kill( debug, step, iter )                                
         return                                                                 
      case(5)                                                                  
         if(debug) write (out,*) '>>>> user specified kill element.'      
         call chk_elem_kill( debug, step, iter )                                
      end select                                                                
c                                                                               
      if( debug ) then                                                         
         call dam_debug                                                         
         write(out,*) '>>>>>>>> leaving chkcrack <<<<<<<'                       
      end if                                                                    
c                                                                               
      return                                                                    
      end 

c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine chk_elem_kill                *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 09/7/2019 rhd              *          
c     *                                                              *          
c     *        checks conditions to see if crack                     *          
c     *        growth will occur by element extinction               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine chk_elem_kill( debug, step, iter )                             
      use global_data ! old common.main
      use elem_extinct_data, only : dam_blk_killed, dam_state,                  
     &                              kill_order_list, dam_print_list
      use main_data, only : output_packets, packet_file_no                      
      use damage_data                                                           
      implicit none
c
      integer :: step, iter
      logical :: debug                                                    
c                                                                               
      integer, parameter :: max_local_list=200 ! max elements that
c                                                can be killed this call    
c                                       
      integer :: local_packet(max_local_list), elem, elem_ptr,
     &           local_count, blk, chk_kill, cohes_type, count, 
     &           element, felem, i, span
      integer, parameter :: local_pkt_type(5) = (/ 20, 0, 21, 10, 0 / )                   
      logical :: blk_killed, killed_this_time, killed_found,                
     &           kill_the_elem, elems_left, ldummy1, ldummy2,
     &           local_debug                                                                                                                           
      double precision ::                                                        
     &   dummy_arg, porosity, eps_plas, eps_crit, dummy_arg2,                   
     &   values(20), local_status(max_local_list,3), orig_poros,                
     &   ext_shape, ext_spacing, triaxiality                                                                                                                                
      logical :: fgm_cohes, ext_gurson, option_exponential,                        
     &           option_ppr, local_write_packets, do_kill_in_order,                
     &           found_exponential, found_ppr, found_cavit,                        
     &           option_cavit                                                      
c        
      local_debug = .false. ! debug
      if( local_debug ) then      
         write(out,*) '... entered chk_elem_kill ...'                                                    
         write(out,*) '    ... no_killed_elems:',no_killed_elems
      end if                                                                    
c                                                                               
c              1) Loop for all elements. process killable elements.
c                 kill as required.                 
c                                                                               
      killed_this_time  = .false.                                                
      elems_left        = .false.                                                
      local_count       = 0                                                           
      found_exponential = .false.                                               
      found_ppr         = .false.                                               
      found_cavit       = .false.                                               
c     
                              
      do elem = 1, noelem
c          
         elem_ptr = dam_ptr(elem)                                               
         if( elem_ptr .eq. 0 ) cycle     ! element not killable                                      
         if( local_debug ) write(out,*) 'killable element is:',elem                           
c                                                                               
c                If element has been previously killed, update it's             
c                unloading state vector, and then skip the rest of the          
c                loop. we do this for traction-separation law as well           
c                since other routines look at dam_state.                        
c                                                                               
         if( dam_state(elem_ptr) .gt. 0 ) then ! element killed already    
            if( dam_state(elem_ptr) .le. max_dam_state )                       
     &           dam_state(elem_ptr) = dam_state(elem_ptr) + 1                  
            cycle                                                               
         end if                                                                 
c                                                                               
c                 calculate damage parameter(s) for element,                    
c                 then check it. if exceeded, set-up data structures for        
c                 internal force release using fixed no.                        
c                 of load steps or the traction-separation law.                 
c                                                                               
         select case( crack_growth_type )
         case( 1, 3 )  ! GTN or SMCS   
            call dam_param( elem, kill_the_elem, debug, porosity,              
     &                      eps_plas, eps_crit, dummy_arg, dummy_arg2,           
     &                      ext_gurson, ext_shape, ext_spacing )
         case( 4 )
            call dam_param_cohes( elem, kill_the_elem, debug,                  
     &                            values, 2 )  
         case( 5 )
            kill_the_elem = .false.
            do i = 1, num_user_kill_elems
             if( elem  == user_kill_list_now(i) ) then
                 kill_the_elem = .true.
                 exit
             end if
            end do
         end select 
c                              
         killed_this_time = killed_this_time .or. kill_the_elem                 
         if( .not. kill_the_elem ) then                                         
           elems_left = .true. ! killable elements remain to be killed                   
           cycle                                                                
         end if                                                                 
c                                                                               
c                  kill now and make output packet data as needed.
c
         call chk_elem_kill_now
         if( output_packets ) call chk_elem_kill_make_packets
c                                                                               
c                 initialize global variables for 1st killed element 
c                 null props and history for this newly killed element          
c                 so subsequent stiffness will be zero. zero stresses,          
c                 start the element force reduction process, constrain          
c                 completely uncoupled nodes if required.                       
c                                                                               
         if( no_killed_elems ) then                                             
            call allocate_damage( 4 )                                           
            if( release_type .eq. 2 ) call allocate_damage( 5 )                 
            call dam_init2( debug )                                             
         end if                                                                 
         num_elements_killed = num_elements_killed + 1                       
         call update_killed_energy( elem )                                   
         call kill_element( elem, debug )                                    
         call store_ifv( elem, elem_ptr, debug )                             
         call update_node_elecnt( elem, debug )                              
         if ( release_type .eq. 2 ) call growth_set_dbar( elem,              
     &             elem_ptr, debug, -1, dummy_arg ) 
c                         
      end do  ! over all model elements                                         
c                                                                               
c                                                                               
c         (2) packet output for newly killed elements at                        
c             for gurson, extended-gurson, smcs, cohesive                       
c                                                                               
      local_write_packets = output_packets .and. local_count .gt. 0             
      if( local_write_packets ) call chk_elem_kill_output_packets
c                                                                               
c         (3) are elements left to be killed in model?
c                                                                               
      if( .not. elems_left ) all_elems_killed = .true.                          
c                                                                               

c         (4) if a killing order has been specified then make sure that         
c             no element "holes" have been left.                                
c                                                                               
      do_kill_in_order = kill_order .and. killed_this_time                      
      if( do_kill_in_order ) call chk_elem_kill_in_order 
c                                                                               

c         (5) check the element blocks to see if all elements in                
c             them have been killed.  Set dam_blk_killed if so.                 
c                                                                               
      if( killed_this_time ) call chk_elem_kill_chk_blk
c                                                                               
c         (6) re-check for free nodes if any elements have ever                 
c             been killed. Just makes triple sure that there are                
c             no free nodes.                                                    
c                                                                               
      if( .not. no_killed_elems ) call chk_free_nodes( debug )                 
c                                                                               

c         (7) if MPI, then all processors need to know what blocks              
c             and elements have been killed.  Send dam_blk_killed               
c             and dam_state.  Also have processors kill any elements            
c             which they own which should be killed.                            
c             this is a dummy for non-MPI                                      
c                                                                               
      call wmpi_send_growth( killed_this_time )                                
c                                                                               
          
      if( num_elements_killed > killed_element_limit ) then
          write(out,9000) killed_element_limit
          call store( ' ','kill_limit_restart.db', ldummy1, ldummy2 )
          call warp3d_normal_stop
      end if
c      
      return    
c
 9000 format(/,1x,">>>>> Number of currently killed (deleted)",
     & " elements has reached the ",
     &       /,7x,"user-specified limit of: ",i6,".",
     &       /,7x,"Restart file: kill_limit_restart.db written.",
     &       /,7x,"Job ended normally."// )
c                                                                               
c
      contains
c     ========
c
      subroutine chk_elem_kill_now
      implicit none
c
c                 kill this element right now. record status for                
c                 subsequent packet output if neeeded. write usual              
c                 output about element death.                                   
c
      local_count = local_count + 1                                          
      select case( crack_growth_type )
      case( 1 )
          if( .not. ext_gurson ) write(out,9000) elem, porosity              
          if( ext_gurson ) write(out,9002) elem, porosity,                   
     &                     ext_spacing, ext_shape                            
      case( 3 )
          write(out,9010) elem, eps_plas, eps_crit                           
      case( 4 )
          cohes_type         = iprops(27,elem)                               
          option_exponential = cohes_type .eq. 4                             
          option_ppr         = cohes_type .eq. 6                             
          option_cavit       = cohes_type .eq. 7                             
          found_exponential  = found_exponential .or.                         
     &                         option_exponential                             
          found_ppr   = found_ppr .or. option_ppr                            
          found_cavit = found_cavit .or. option_cavit                        
          count = 0                                                          
          if( found_exponential ) count = count + 1                          
          if( found_ppr ) count = count + 1                                  
          if( found_cavit ) count = count + 1                                
          if( count .gt. 1 ) then                                            
            write(out,9400) elem                                            
            call die_gracefully                                              
          end if                                                             
          if( option_exponential )                                           
     &       write(out,9200) elem, values(6)/values(7), values(8)            
          if( option_ppr )   write(out,9220) elem                            
          if( option_cavit ) write(out,9230) elem   
      case( 5 )
          write(out,9240) elem                         
      end select
c
      return
c
 9000 format(/,'   >> element death invoked for element: ',i7,                  
     & '.   f: ',f6.3)                                                          
 9002 format(/,'   >> element death invoked for element: ',i7,                  
     & '.   f, W, X: ',3f8.3)                                                   
 9010 format(/,' >> element death option invoked before next step',             
     &       /,'    element: ',i7,' is now killed.',                            
     &       /,'    plastic strain: ',f8.5,' is > limit of: ',f8.5)             
 9200 format(/,'   >> element death invoked for element: ',i7,                  
     & '.   Deff/Dpeak: ',f5.2,' Teff/Tpeak: ',f5.2)                            
 9220 format(/,'   >> element death invoked for element: ',i7,                  
     & '.  PPR cohesive option')                                                
 9230 format(/,'   >> element death invoked for element: ',i7,                  
     & '.  cavit cohesive option')                                              
 9240 format(/,'   >> element death invoked for element: ',i7,                  
     & '.  user-directed removal')                                              
 9400 format(/,'FATAL ERROR: mixed cohesive options not allowed',               
     & /,      '             at present in WARP3D with crack growth.'           
     & /,      '             found while processing element: ',i7,              
     & /,      '             job aborted.' )  
c 
      end subroutine chk_elem_kill_now
c
      subroutine chk_elem_kill_output_packets
      implicit none
c
      write(packet_file_no) local_pkt_type(crack_growth_type),                  
     &                      local_count, step, iter                             
      do elem = 1, local_count                                                  
        element = local_packet(elem)                                            
        select case( crack_growth_type )
        case( 1 )
          orig_poros = props(26,element)                                        
          write(packet_file_no) element,orig_poros,                             
     &                  local_status(elem,1),                                   
     &                  local_status(elem,2),local_status(elem,3)               
        case( 3 )
          write(packet_file_no) element, local_status(elem,1),                  
     &                                  local_status(elem,2)                    
        case( 4 ) 
          write(packet_file_no) element, local_status(elem,1),                  
     &          local_status(elem,2), int(local_status(elem,3))                 
        end select
      end do   
c                                                                 
      return
      end subroutine chk_elem_kill_output_packets
c
      subroutine chk_elem_kill_in_order
      implicit none
c
      if ( debug ) write (out,*) '>>> checking for holes....'                   
      killed_found = .false.                                                    
      do chk_kill = num_kill_order_list,1, -1                                   
        elem     = kill_order_list(chk_kill)                                    
        elem_ptr = dam_ptr(elem)                                                
        if( .not. killed_found ) then                                          
          if( dam_state(elem_ptr) .ne. 0 ) killed_found = .true.               
        else                                                                    
          if( dam_state(elem_ptr) .eq. 0 ) then                                
             write (out,9100) elem                                              
             num_elements_killed = num_elements_killed + 1                      
             call update_killed_energy( elem )                                  
             call kill_element( elem, debug )                                   
             call store_ifv( elem, elem_ptr, debug )                            
             call update_node_elecnt( elem, debug )                             
             if ( release_type .eq. 2 )                                         
     &             call growth_set_dbar( elem, elem_ptr, debug,                 
     &                                        -1, dummy_arg )                   
          end if                                                                
        end if                                                                  
      end do 
c   
      return
c
 9100 format(/,' >> element death option invoked before next step',             
     &       /,'    element: ',i7,' is now killed to make crack',               
     &       /,'    face uniform.')                                             
c
      end subroutine   chk_elem_kill_in_order 
c                                                             
      subroutine chk_elem_kill_make_packets
      implicit none
c
      if( local_count .gt. max_local_list ) then                           
         write(out,9300)                                                   
         call die_gracefully                                               
      end if                                                               
      select case( crack_growth_type )
      case( 1 )
           local_packet(local_count) = elem                                
           local_status(local_count,1) = porosity                          
           local_status(local_count,2) = ext_spacing                       
           local_status(local_count,3) = ext_shape                         
      case( 3 ) 
           local_packet(local_count) = elem                                
           local_status(local_count,1) = eps_plas                          
           local_status(local_count,2) = eps_crit                          
      case( 4 )
           if( option_exponential ) then                                   
             local_packet(local_count) = elem                              
             local_status(local_count,1) = values(6)/values(7)             
             local_status(local_count,2) = values(8)                       
             local_status(local_count,3) = 4.0                             
           end if                                                          
           if( option_ppr ) then                                           
             local_packet(local_count) = elem                              
             local_status(local_count,1) = values(5)                       
             local_status(local_count,2) = values(6)                       
             local_status(local_count,3) = 6.0                             
           end if                                                          
           if( option_cavit) then                                          
             write(*,*) ' cavit death packets not implemented'             
             call die_abort                                                
           end if  
      end select 
c                                                       
      return
c
 9300 format(/,'FATAL ERROR: list length exceeded in chk_elem_kill',            
     & /,      '             job aborted.' )                                    
c
      end subroutine chk_elem_kill_make_packets
c
      subroutine chk_elem_kill_chk_blk
      implicit none
c
      do blk = 1, nelblk                                                        
        if( dam_blk_killed(blk) ) cycle                                        
        span  = elblks(0,blk)                                                   
        felem = elblks(1,blk)                                                   
        if( iand( iprops(30,felem), 2 ) .eq. 0 ) cycle                         
        blk_killed = .true.                                                     
        do i = 1, span                                                          
         if( dam_state(dam_ptr(felem+i-1)) .eq. 0 ) blk_killed =.false.        
        end do                                                                  
        if( blk_killed ) dam_blk_killed(blk) = .true.                          
      end do       
c                                                             
      return
      end subroutine chk_elem_kill_chk_blk
c
      end subroutine chk_elem_kill                                                                      
                                                                      
c                                          
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_param                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 04/29/2019 rhd             *          
c     *                                                              *          
c     *     for a killable element not yet killed, determine if the  *          
c     *     element should be killed now. isolating decision here    *          
c     *     makes possible different kinds of killing criteria       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dam_param( elem, kill_now, debug, porosity,                    
     &                      eps_plas, eps_crit, sig_mean, sig_mises,            
     &                      ext_gurson, ext_shape, ext_spacing )
      use global_data, only : iprops, out  ! old common.main
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
c        eps_plas  -- (output)  for SMCS type models. this the                  
c                               average plastic strain over element.            
c                               Just used for output                            
c        eps_crit  -- (output)  for SMCS type models. this the                  
c                               critical plastic strain value.                  
c                               Just used for output                            
c        sig_mean  -- (output)  all models. average mean stress                 
c                               over element. just used for output              
c                               messages                                        
c        sig_mises -- (ouput)   all models. average mises stress                
c                               over element. just used for output.             
c        ext_gurson -- (output) logical = .true. if element death using         
c                               the extended gurson model                       
c        ext_shape -- (output)  W (shape) parameter for extended gurson model   
c        ext_spacing -- (output) X (spacing) parameter for extended             
c                                gurson model    
c        triaxiality -- (output) plastic strain weighted average of
c                                triaxiality over loading history                               
c                                                                               
c                                                                               
      use main_data,       only : elems_to_blocks                               
      use elem_block_data, only : history_blocks, eps_n_blocks,                 
     &                            urcs_n_blocks, history_blk_list               
      use damage_data, only     : crack_growth_type 
c                                                         
      implicit none                                                    
c                                                                               
c              parameter declarations                                           
c                             
      integer :: elem                                                  
      logical :: kill_now, debug, ext_gurson                                       
      double precision :: porosity, eps_plas, eps_crit, sig_mean,
     &                    sig_mises, ext_shape, ext_spacing,
     &                    triaxiality  
      double precision, dimension(:), pointer :: history, urcs_n, eps_n                         
c                                                                               
c              local declarations                                               
c       
      integer :: mat_model
      double precision :: mean_zeta, mean_omega
      double precision, parameter :: zero = 0.d0, one = 1.0d0                         
c   
      ext_gurson  = .false.                                                     
      ext_shape   = one                                                         
      ext_spacing = zero                                                        
c                                                                               
      select case( crack_growth_type )                         
c                                                                               
c                                                                               
c           1 gurson model- the material model 3 is the standard         
c                           Gurson model, 6 is the extended                
c                           Gurson model. call model dependent routine     
c                           to assess element status for killing and       
c                           to compute parameters.                         
c                                                                               
      case( 1) 
       mat_model = iprops(25,elem)                                               
       if( mat_model .eq. 3 ) then                                               
         call dam_param_gt( elem, kill_now, debug, porosity,                    
     &                      sig_mean, sig_mises )                               
       else if( mat_model .eq. 6 ) then                                          
         ext_gurson = .true.                                                    
         call dam_param_agt( elem, kill_now, debug, porosity,                   
     &                       sig_mean, sig_mises, ext_spacing,                  
     &                       ext_shape )                                        
       else                                                                      
         write(out,9100) 1; write(out,9120)                                                        
         call die_gracefully                                                       
       end if                                                                    
c                                                                               
c           2  ctoa - not processed here                  
c                                                                               
      case(2) 
       write(out,9100) 2; write(out,9120)                                                      
       call die_gracefully                                                       
c                                                                               
c           3 smcs - failure based on a stress modified critical strain     
c
      case(3)   
       call dam_param_3_get_values( elem, debug, eps_plas, eps_crit, 
     &                              sig_mean, sig_mises, triaxiality, 
     &                              mean_zeta, mean_omega, 1,
     &                              kill_now )                
c                                                                               
c           4 cohesive - not processed here               
c                                                                                                                                                             
      case(4)
        write(out,9100) 4; write(out,9120)                                                          
        call die_gracefully                                                       
c
      end select  
c
      return                                                              
c                                                                               
 9100 format('>>>> FATAL ERROR: invalid model no. in dam_param @:',i1 )                
 9120 format('                  job terminated....')                
c    
      end
c     ****************************************************************          
c     *                                                              *          
c     *              subroutine dam_param_3_get_values               *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 09/7/2019 rhd              *          
c     *                                                              *          
c     *             SMCS compute values for this element             *
c     *                                                              *          
c     ****************************************************************          
c
                                                                     
      subroutine dam_param_3_get_values( elem, debug,                     
     &                      eps_plas, eps_crit, sig_mean, sig_mises,            
     &                      triaxiality, mean_zeta, mean_omega,
     &                      dowhat, kill_now )                
      use global_data, only : iprops, nstr, nstrs, out ! old common.main
      use main_data,       only : elems_to_blocks                               
      use elem_block_data, only : history_blocks, eps_n_blocks,                 
     &                            urcs_n_blocks, history_blk_list               
      use damage_data, only     : dam_ptr, smcs_gamma, smcs_alpha, 
     &                            smcs_beta,  smcs_type, 
     &                            smcs_beta_1,smcs_beta_2, 
     &                            smcs_a_plus, smcs_a_minus,
     &                            smcs_kappa, smcs_alpha_1, 
     &                            smcs_alpha_2, 
     &                            smcs_cutoff_triaxiality                                                           
      use elem_extinct_data, only : smcs_old_epsplas, smcs_weighted_T,
     &                              smcs_weighted_zeta
c                                                                               
c        elem      -- (input)   element number to be checked for                
c                               killing                                         
c        debug     -- (input)   .true. if the calling routine wants             
c                               debugging info output                           
c        eps_plas  -- (output)  for SMCS type models. this the                  
c                               average plastic strain over element.            
c                               Just used for output                            
c        eps_crit  -- (output)  for SMCS type models. this the                  
c                               critical plastic strain value.                  
c                               Just used for output                            
c        sig_mean  -- (output)  all models. average mean stress                 
c                               over element. just used for output              
c                               messages                                        
c        sig_mises -- (ouput)   all models. average mises stress                
c                               over element. just used for output.             
c        triaxiality -- (output) plastic strain weighted average of
c                                triaxiality over loading history 
c        mean_zeta -- (output)  Lode angle parameter average weighted
c                               by plastic strain
c        mean_omega -- (output)  Nahshon-Hutchinson Lode angle parameter
c                                average weighted by plastic strain
c        dowhat    --  (input)   = 1 return damage variables - no update.
c                                    set kill now flag
c                                = 2 just return values - no update 
c                                = 3 return values for non-killable elem
c                                = 4 compute/update damage variables  
c        kill_now  --  (output)  if dowhat = 1, set flag if critical
c                                condition reached                           
c                                                                                                                                                              
      implicit none                                                    
c                                                                               
c               parameter declarations                                           
c                             
      integer :: elem, dowhat                                                  
      logical :: debug, kill_now                                       
      double precision :: eps_crit, sig_mean, sig_mises, triaxiality,
     &                    eps_plas, mean_zeta, mean_omega
c                                                                               
c               local declarations                                               
c
      double precision, dimension(:), pointer :: urcs_n, eps_n                         
      integer :: gp, blk, epsoffset, mat_model, elem_type,
     &           ngp, rel_elem, sigoffset, j, elem_ptr
      double precision ::  sig(6), eps(6), fpngp,
     &                     mean_triaxiality,
     &                     epspls_vec(27),s11, s22, s33, s12, s13, s23,
     &                     j2, j3, zeta, omega, t1, t2, t3, t4,
     &                     p_omega, omega2
      logical :: update, local_debug, threed_solid_elem, l1, l2, l3,
     &           l4, l5, l6, l7, l8, l9, l10
     &          
      double precision, parameter :: zero = 0.d0,
     &     third = 1.d0/3.0d0, iroot2 = 1.0d0/dsqrt(2.0d0), six = 6.0d0,
     &     one = 1.0d0, two = 2.0d0, onept5 = 1.5d0, three = 3.d0,
     &     root3 = dsqrt(3.0d0), half = 0.5d0, thirteenpt5 = 13.5d0,
     &     small_plastic_strain = 1.0d-07, small_mises = 1.0d-07,
     &     type_2_tol = 1.0d-06, nine = 9.d0
c                              
c               stress/strain vectors are (x,y,z,xy,yz,xz).                
c
      update      = dowhat .eq. 4
      local_debug = .false.  !   elem == 16
      if( local_debug )then
         write(out,*) '... dam_param_3_get_values, dowhat: ...', dowhat
      end if
      kill_now = .false.
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
      eps_crit    = nine
      sig_mean    = zero
      sig_mises   = zero
      triaxiality = zero
      mean_zeta   = zero
      mean_omega  = zero                               
      sig         = zero
      eps         = zero  !   strains not used current code       
      call set_element_type( elem_type, threed_solid_elem, l1, l2, l3,
     &                       l4, l5, l6, l7, l8, l9, l10 )  
      if( .not. threed_solid_elem ) return                                 
c
      call mm_return_values( "plastic_strain", elem, epspls_vec, ngp )
c
c               average stresses over all integration points
c
      if( ngp == 8 ) then  ! for optimizer
       do gp = 1, 8    
        do j = 1, 6                                                     
         sig(j) = sig(j) + urcs_n(sigoffset+j)  
         eps(j) = eps(j) + eps_n(epsoffset+j)
        end do                                      
        eps_plas  = eps_plas + epspls_vec(gp)
        epsoffset = epsoffset + nstr                                         
        sigoffset = sigoffset + nstrs                                        
       end do   
      else  ! general case
       do gp = 1, ngp    
        do j = 1, 6                                                     
         sig(j) = sig(j) + urcs_n(sigoffset+j)  
         eps(j) = eps(j) + eps_n(epsoffset+j)
        end do                                      
        eps_plas  = eps_plas + epspls_vec(gp)
        epsoffset = epsoffset + nstr                                         
        sigoffset = sigoffset + nstrs                                        
       end do 
      end if  
c
c               average mean stress, mises stress, plastic strain,
c               Lode parameters
c
      fpngp     = one / dble( ngp )
      sig = sig * fpngp  !  averaged stresses at element center
      sig_mean  = (sig(1) + sig(2) + sig(3)) * third                    
      s11 = sig(1) - sig_mean
      s22 = sig(2) - sig_mean
      s33 = sig(3) - sig_mean
      s12 = sig(4)
      s13 = sig(5)
      s23 = sig(6)
      j2 = half * ( s11**2 + s22**2 + s33**2 + 
     &              two*(s12**2 + s23**2 + s13**2))
      j3 = s11*s22*s33 + two*s12*s23*s13 - s11*s23*s23 -
     &     s22*s13*s13 - s33*s12*s12
      sig_mises = sqrt( three * j2 )
      zeta = zero  !  called \xi in manual
      omega = zero
      if( sig_mises > small_mises ) then
        zeta = onept5*root3 * j3 / j2**onept5
        omega = one - zeta*zeta 
      end if
      eps_plas  = fpngp * eps_plas   
      if( local_debug ) write(out,9100) sig_mean,sig_mises, eps_plas
c
      if( dowhat == 3 ) then ! compute values for non-killable elem
         eps_crit = nine
         mean_zeta = zeta
         mean_omega = omega
         triaxiality = 100.0d0   ! in cases mises -> 0                                             
         if( sig_mises .gt. small_mises )
     &         triaxiality = sig_mean / sig_mises 
         return
       end if
c
c               update plastic strain weighted triaxiality, zeta,
c               and omega. update saved plastic strain.
c               only update if instantaneous triaxiality exceeds
c               user-specified cutoff value
c              
      if( update ) call dam_param_3_get_values_a
c
c               if no plastic strain yet, set defaults and return.
c               cannot compute weighed averages over history (divide
c               by zero plastic strain).
c
      if( eps_plas < small_plastic_strain ) then 
        triaxiality = zero
        eps_crit = nine  ! a large value for return
        mean_zeta  = zero
        mean_omega = zero
        return
      end if
c
c               use weighted triaxiality, zeta, omega values to 
c               set critical plastic strain
c
      mean_triaxiality =  smcs_weighted_T(elem_ptr) / eps_plas
      triaxiality = mean_triaxiality  ! to return
      mean_zeta = smcs_weighted_zeta(elem_ptr) / eps_plas
      mean_omega = one - mean_zeta**2
      eps_crit = nine
      if( triaxiality <= smcs_cutoff_triaxiality ) return
c
      if( smcs_type .eq. 1 ) then
          t2 = exp(-smcs_kappa*abs(mean_zeta))
          t1 = smcs_alpha * exp(-smcs_beta * mean_triaxiality)
          eps_crit = smcs_gamma + t1*t2
        if( local_debug ) write(out,9110) mean_triaxiality,eps_crit

      elseif( smcs_type .eq. 2 ) then
          t3 = exp(-smcs_kappa*abs(mean_zeta))
          t2 = smcs_beta_2 * exp(-smcs_a_minus * mean_triaxiality)
          t1 = smcs_beta_1 * exp(smcs_a_plus * mean_triaxiality)
          t4 = t1 - t2
          if( abs(t4) < type_2_tol .or. t4 < zero ) then
            eps_crit = nine ! a large value
          else
            eps_crit = smcs_gamma + t3 / t4
          end if
      else
          p_omega = (one - mean_omega)**2
          t1 = smcs_alpha_1*exp(-smcs_beta_1 * mean_triaxiality )
          t2 = smcs_alpha_2*exp(-smcs_beta_2 * mean_triaxiality )
          eps_crit = (one-p_omega)*t1 + p_omega*t2
      end if
c      
      kill_now = ( eps_plas - eps_crit ) >= zero
c
      return
c
 9100 format(10x,'mean, mises, eps_pls:',2f10.3,f12.6)
 9110  format(10x,'T_mean, eps crit: ',f10.3,f12.6) 
c
      contains
c     ========
c
      subroutine dam_param_3_get_values_a
      implicit none
c
      double precision :: deps_plas, now_triaxiality
c
      deps_plas = eps_plas - smcs_old_epsplas(elem_ptr)
      smcs_old_epsplas(elem_ptr) = eps_plas 
      if( sig_mises .gt. small_mises ) then
        now_triaxiality = sig_mean / sig_mises 
        if( now_triaxiality < smcs_cutoff_triaxiality ) return
        smcs_weighted_T(elem_ptr)    = smcs_weighted_T(elem_ptr) +
     &                                 deps_plas * now_triaxiality
        smcs_weighted_zeta(elem_ptr) = smcs_weighted_zeta(elem_ptr) +
     &                                 deps_plas * zeta
      end if
c
      return
      end subroutine dam_param_3_get_values_a

      end subroutine dam_param_3_get_values


c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_init2                    *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 04/19/2019 rhd             *          
c     *                                                              *          
c     *        This routine is called when damage to an element      *          
c     *        has occurred.  It initializes arrays to support       *          
c     *        releasing of the internal forces present on element   *          
c     *        at extinction.                                        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine dam_init2( debug )                                             
      use global_data ! old common.main
c                                                                               
      use elem_extinct_data, only : dam_node_elecnt, dam_face_nodes             
      use main_data, only : incmap, incid, elems_to_blocks,                     
     &                      inverse_incidences                                  
      use elem_block_data, only : cdest_blocks                                  
      use damage_data                                                           
c                                                                               
      implicit none                                                    
c                                                                               
      logical :: debug
c
      integer :: i, elem, elem_ptr, incptr, num_enodes, k1, k2,
     &           face_count, blk, rel_elem, enode, snode
      logical, parameter :: local_debug = .false.                                            
      double precision :: node_coord(3), coord
      double precision, parameter :: plane_tol = 0.01d0                                 
      integer, dimension(:,:), pointer :: cdest                                 
c                                                                               
      if( debug ) write(out,*) '>>>> in dam_init2'                             
c                                                                               
      no_killed_elems = .false.                                                 
      do i = 1, nonode                                                          
         dam_node_elecnt(i) = inverse_incidences(i)%element_count               
      end do                                                                    
      if( debug ) then                                                         
         do i = 1, nonode                                                       
            write (out,9030) i, dam_node_elecnt(i)                              
         end do                                                                 
      end if                                                                    
c                                                                               
c               set up the force release by traction-displacement               
c               separation law if that option is being used.                    
c               for each killable element on the crack plane,                   
c               find the 4 nodes on the face opposite the crack                 
c               plane face used to define the element deformation.              
c               we examine only the first 8 nodes of the element                
c               (corner nodes).                                                 
c                                                                               
      if( release_type .ne. 2 ) then
        if( debug ) write(out,*) '<<<< leaving dam_init2'                        
        return                                                                    
      end if
c
      do elem  = 1, noelem                                                      
        elem_ptr = dam_ptr(elem)                                                
        if( elem_ptr .eq. 0 ) cycle                                            
        incptr     = incmap(elem)-1                                             
        num_enodes = iprops(2,elem)                                             
        k1         = num_enodes                                                 
        k2         = 2*num_enodes                                               
        face_count = 0                                                          
        blk        = elems_to_blocks(elem,1)                                    
        rel_elem   = elems_to_blocks(elem,2)                                    
        cdest      => cdest_blocks(blk)%ptr                                     
        do enode = 1, 8                                                         
          snode         = incid(incptr+enode)                                   
          node_coord(1) = c( cdest(enode,rel_elem) )                            
          node_coord(2) = c( cdest(k1+enode,rel_elem) )                         
          node_coord(3) = c( cdest(k2+enode,rel_elem) )                         
          coord         = node_coord(crk_pln_normal_idx)                        
          if( abs(coord-crack_plane_coord) .le. plane_tol*                     
     &        gurson_cell_size ) cycle                                         
          face_count = face_count + 1                                           
          dam_face_nodes(face_count,elem_ptr) = snode                           
        end do   ! enode                                                               
        if( face_count .ne. 4 ) then                                           
             write(out,*) 'FATAL error 1 in dam_init2'                          
             write(out,*) 'invalid plane definition for growth'                 
             write(out,*)  elem, elem_ptr, num_enodes,face_count                
             call die_gracefully                                                
             stop                                                               
        end if                                                                  
      end do ! elem                                                                  
c                                                                              
      if( local_debug ) then                                                  
       write(out,*) '> element release type 2. face node table'                
       do elem  = 1, noelem                                                    
        elem_ptr = dam_ptr(elem)                                               
        if( elem_ptr .eq. 0 ) cycle                                           
        write(out,9050) elem, dam_face_nodes(1,elem_ptr),                      
     &       dam_face_nodes(2,elem_ptr), dam_face_nodes(3,elem_ptr),         
     &       dam_face_nodes(4,elem_ptr)                                      
       end do                                                                  
      end if                                                                   
c                                                                                                                                             
      if( debug ) write(out,*) '<<<< leaving dam_init2'                        
      return                                                                    
 9030 format(' dam_node_elecnt(',i4,')=',i2)                                    
 9050 format(' elem: ',i7, '  face nodes: ',8i8)                                
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine chk_free_nodes               *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 04/19/2019 rhd             *          
c     *                                                              *          
c     *        check each node in the structure to                   *          
c     *        see if they are no longer connected to any element.   *          
c     *        if so, it constrains them.                            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine chk_free_nodes( debug )                                        
      use global_data ! old common.main
c                                                                               
      use elem_extinct_data, only : dam_node_elecnt                             
      use main_data,         only : cnstrn, cnstrn_in                          
      use damage_data                                                           
c                                                                               
      implicit none                                                    
c                                                                               
      integer :: node, num_dof, j, glbal_dof
      logical :: debug, local_debug                                                                                                                                          
      double precision, parameter :: zero = 0.0d0                                                          
c
      local_debug = debug
      if( local_debug ) write(out,*) '>>> in chk_free_nodes <<<'                     
c                                                                               
c              Loop over each node -- if a node is now attatched to             
c              no elements, constrain it.                                       
c              Note: this could be more efficient -- we are looping             
c                over all the nodes each time, and we could just                
c                constrain the newly zeroed nodes.  This, however,              
c                will also take care of the problem of a user                   
c                putting in a new set of constraints.                           
c                                                                               
      do node = 1, nonode                                                       
       if( dam_node_elecnt(node) .ne. 0 ) cycle                             
       if( local_debug ) then                                                   
         write(out,*) 'free node ',node                                   
         write (out,9010) node,cnstrn(dstmap(node)),                      
     &          cnstrn(dstmap(node)+1),cnstrn(dstmap(node)+2)               
       end if                                                              
       do j = 1, 3    ! num_dof                                                   
         glbal_dof            = dstmap(node) + j-1                        
         cnstrn(glbal_dof)    = zero                                      
         cnstrn_in(glbal_dof) = zero                                      
         if( cstmap(glbal_dof) .eq. 0 ) then                             
           cstmap(csttail) = glbal_dof                                   
           csttail         = glbal_dof                                   
         end if                                                           
       end do                                                              
       cstmap(csttail) = -1                                                
      end do                                                                    
c                                                                               
 9999 continue                                                                  
      if ( local_debug ) write(out,*) '>>> in chk_free_nodes <<<'                     
      return                                                                    
 9010 format('  old constraints for node ',i7,' :',3(1x,e14.6))                 
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_print                    *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 9/29/2020 rhd              *          
c     *                                                              *          
c     *     print status of killable elements or released nodes at   *
c     *     the beginning of a load step.                            *       
c     *     optional smcs states output for smcs model
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine dam_print( step, iter )                                        
      use global_data ! old common.main
      use damage_data                                                           
      implicit none    
c
      integer :: step, iter  
      logical :: doprint                                              
c                                                                               
c           check to make sure crack growth is on and some type of
c           output is requested
c                                                                               
      select case( crack_growth_type )                                           
      case(0)       ! no action                                                             
      case(1)                                                                  
         doprint = print_status .or. print_top_list
         if( .not. doprint ) return
         call dam_print_elem1( step, iter )                                  
      case(2)
         doprint = print_status .or. print_top_list
         if( .not. doprint ) return
         if( const_front ) then                                                
           call dam_print_front( step, iter)                                
         else                                                                 
           call dam_print_node( step, iter )                                
         end if    
      case(3)                                                            
         doprint = print_status .or. print_top_list .or. smcs_states
         if( .not. doprint ) return
         call dam_print_elem3( step, iter )                                  
      case(4)
         doprint = print_status .or. print_top_list
         if( .not. doprint ) return
         call dam_print_elem4( step, iter )                                  
      end select                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_debug                    *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 06/23/94                   *          
c     *                                                              *          
c     *     This routine prints out all of the crack growth          *          
c     *     parameters.                                              *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               

      subroutine dam_debug                                                      
      use global_data ! old common.main
      use elem_extinct_data                                                     
      use node_release_data                                                     
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
      double precision                                                          
     &     zero                                                                 
      data zero / 0.0 /                                                         
                                                                                
c                                                                               
      write (out,*)                                                             
      write (out,*) '>>>>>>>> dumping all relevent crack growth info:'          
      write (out,*)                                                             
      write (out,*) '> crack_growth_type:',crack_growth_type                    
c                                                                               
c              -- Gurson crack growth                                           
c                                                                               
      if( crack_growth_type.eq.1 ) then                                         
         write (out,*) '> max_dam_state:', max_dam_state                        
         write (out,*) '> num_kill_elem:', num_kill_elem                        
         write (out,*) '> dam_ptr:'                                             
         do i=1, noelem                                                         
            write (out,'(" dam_ptr(",i4,")=",i4)') i,dam_ptr(i)                 
         end do                                                                 
         write (out,*) '> dam_state:'                                           
         do i = 1, num_kill_elem                                                
            write (out,'("dam_state(",i3,")=",i3)')i,dam_state(i)               
         end do                                                                 
         write (out,*) '> dam_ifv:'                                             
         do i=1, num_kill_elem                                                  
            write (out,*) 'ifv entries for ptr:',i                              
            write(out,9000) (dam_ifv(j,i), j=1, mxedof)                         
         end do                                                                 
         write (out,'(" > porosity_limit:",e14.6)') porosity_limit              
         write (out,*) '> dam_blk_killed:'                                      
         do blk = 1, nelblk                                                     
            if( dam_blk_killed(blk) ) then                                      
               write (out,'(" dam_blk_killed(",i3,") is true")') blk            
            else                                                                
               write (out,'(" dam_blk_killed(",i3,") is false")') blk           
            endif                                                               
         end do                                                                 
         if( kill_order ) then                                                  
            write (out,*) '> kill_order_list:'                                  
            do i = 1, num_kill_order_list                                       
               write (out,'("  entry:",i7," elem:",i6)')i,                      
     &              kill_order_list(i)                                          
            end do                                                              
         end if                                                                 
         if( .not.no_killed_elems ) then                                        
            write(out,*) '> elems had been killed'                              
            write(out,*) '> dam_node_elecnt:'                                   
            do i = 1, nonode                                                    
               write (out,'(" dam_node_elecnt(",i5,")=",i3)')i,                 
     &              dam_node_elecnt(i)                                          
            end do                                                              
            write (out,*) '> dam_dbar_elems(current height, fraction)'          
            if( release_type .eq. 2 ) then                                      
               do i = 1, num_kill_elem                                          
                  write (out,'(e13.6,2x,e13.6)') dam_dbar_elems(1,i),           
     &                 dam_dbar_elems(2,i)                                      
               end do                                                           
            end if                                                              
         else                                                                   
            write (out,*) '> no elems have been killed.'                        
         end if                                                                 
c                                                                               
c              -- Node release crack growth                                     
c                                                                               
      else if (crack_growth_type.eq.2) then                                     
c                                                                               
         write (out,*) '> Number of crack plane nodes:',                        
     &        num_crack_plane_nodes                                             
         write (out,'(" > half of CTOA for release:",e13.6)')                   
     &        critical_angle                                                    
         write (out,*) '> entry, node, state, ifv:'                             
         do i=1,num_crack_plane_nodes                                           
            write (out,'(2x,i7,i6,i6,3x,e13.6)')i,crack_plane_nodes(i),         
     &           crkpln_nodes_state(i), crkpln_nodes_react(i)                   
         enddo                                                                  
         write (out,*) '> global nodes => crack plane nodes'                    
         do i=1, nonode                                                         
            write (out,*) ' ',i,':',inv_crkpln_nodes(i)                         
         enddo                                                                  
         write (out,*) '> crack front list:'                                    
         call write_list                                                        
         write (out,*) '> neighbor_nodes:'                                      
         do i=1, num_crack_plane_nodes                                          
            write (out,*) ' ',i,':  for node:',crack_plane_nodes(i),            
     &           ' # neighbors=', num_neighbors(i),' which are:'                
            write (out,'(7x,12i5)') (neighbor_nodes(j,i),                       
     &           j=1,num_neighbors(i))                                          
         enddo                                                                  
c                                                                               
         if (release_type .eq. 2) then                                          
            write(out,'(" > charlen,rel.frac,rel.height:",3e13.6)')             
     &           char_length, release_fraction, release_height                  
            write (out,*) '> node_release_frac:'                                
            do i = 1, num_crack_plane_nodes                                     
               write (out,'(2x,i7,2x,e13.6)')i,node_release_frac(i)             
            enddo                                                               
         endif                                                                  
c                                                                               
         if ( overshoot_control_crk_grth ) then                                 
            write (out,*) '> overshoot control is on:'                          
            write (out,'("   -control load fact:",e13.6)')                      
     &           control_load_fact                                              
            write (out,'("   -old load fact:",e13.6)') old_load_fact            
            write (out,'("   -min load fact:",e13.6)') min_load_fact            
            write (out,'("   -overshoot_limit:",e13.6)')                        
     &           overshoot_limit                                                
            write (out,'("   -CTOA range:",e13.6)') CTOA_range                  
            if (overshoot_allocated) then                                       
               write (out,*) '   -overshoot stuff is allocated'                 
               write (out,'("   -old angles:")')                                
               write (out,'("     entry  angles")')                             
               do i=1, num_crack_plane_nodes                                    
                  write (out,'(6x,i4,1x,12e13.6)') i,                           
     &                 (old_angles_at_front(i,j),j=1,12)                        
               enddo                                                            
            else                                                                
               write (out,*) '   -overshoot stuff is not allocated'             
            endif                                                               
         else                                                                   
            write (out,*) '> no overshoot control'                              
         endif                                                                  
c                                                                               
         if ( load_size_control_crk_grth ) then                                 
            write (out,*) '> load size control is on'                           
            write (out,'("   -perm load fact:",e13.6)') perm_load_fact          
            write (out,'("   -min steps:",i3)') min_steps_for_release           
         else                                                                   
            write (out,*) '> no load size control'                              
         endif                                                                  
c                                                                               
      if (const_front) then                                                     
       write (out,*) '> const front is set.'                                    
       write (out,*) '  -master nodes:'                                         
       do i=1, num_crack_fronts                                                 
         write (out,'(5x,i7)') master_nodes(i)                                  
       enddo                                                                    
       write (out,*) '  -crack_front_list:'                                     
       do i=1, num_crack_fronts * num_nodes_grwinc                              
        write (out,'(5x,i3,":",40i7)') i,                                       
     &             (crack_front_list(i,j),j=1,num_nodes_thick)                  
       enddo                                                                    
c                                                                               
       write (out,*) '> measured CTOA'                                          
       write (out,'("   -init_ctoa_dist:",e13.6)') init_ctoa_dist               
       write (out,'("   -ctoa_dist:",e13.6)') ctoa_dist                         
       write (out,*) '  -num_nodes_back =',num_nodes_back                       
       write (out,*) '  -master list (first 20 nodes):'                         
       do i=1, num_crack_fronts                                                 
         write (out,'(5x,i7,":",20i6)') i, (master_lines(i,j),                  
     &              j = 1, 20)                                                  
       enddo                                                                    
      endif                                                                     
c                                                                               
         write (out,*)                                                          
c                                                                               
c              -- No crack growth                                               
c                                                                               
      else if (crack_growth_type.eq.0) then                                     
         write (out,*) '> no crack growth specified'                            
c                                                                               
c              -- Error                                                         
c                                                                               
      else                                                                      
         write (out,*) ' Incorrect crack_growth_type:',                         
     &        crack_growth_type                                                 
      endif                                                                     
c                                                                               
      write (out,*)                                                             
      write (out,*) '<<<<<<<< finished dumping crack growth info.'              
      write (out,*)                                                             
 9999 return                                                                    
c                                                                               
 9000 format(5(3x,4e14.6,/))                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                 subroutine growth_set_dbar                   *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 04/19/2019 rhd             *          
c     *                                                              *          
c     *  supports traction-separation law release option. compute    *          
c     *  elongation of element normal to crack plane when first      *          
c     *  killed (option=0) or at some load step later (option=1)     *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine growth_set_dbar( elem, elem_ptr, debug, action,                
     &                            dbar_now )                                    
      use global_data ! old common.main
      use elem_extinct_data, only : dam_face_nodes, dam_dbar_elems              
      use damage_data                                                           
c                                                                               
      implicit none
c                
      integer :: elem, elem_ptr, action                                    
      logical :: debug
c
      integer :: i, snode
      integer, parameter :: num_face_nodes = 4
      logical, parameter :: local_debug = .false.                                                
      double precision :: sum, node_displ(3), dbar_now          
      double precision, parameter :: zero = 0.d0                           
c                                                                               
      if( debug ) write (out,*) '>>>> in growth_set_dbar'                      
c                                                                               
      sum = zero       
c                                                         
      if( action .lt. 0 ) then                                                 
        do i = 1, num_face_nodes                                                
          snode         = dam_face_nodes(i,elem_ptr)                            
          node_displ(1) = u(dstmap(snode)+0)                                    
          node_displ(2) = u(dstmap(snode)+1)                                    
          node_displ(3) = u(dstmap(snode)+2)                                    
          sum           = sum + node_displ(crk_pln_normal_idx)                  
        end do                                                                  
      end if   
c                                                                 
      if( action .gt. 0 ) then                                                 
        do i = 1, num_face_nodes                                                
          snode         = dam_face_nodes(i,elem_ptr)                            
          node_displ(1) = u(dstmap(snode)+0) + du(dstmap(snode)+0)              
          node_displ(2) = u(dstmap(snode)+1) + du(dstmap(snode)+1)              
          node_displ(3) = u(dstmap(snode)+2) + du(dstmap(snode)+2)              
          sum           = sum + node_displ(crk_pln_normal_idx)                  
        end do                                                                  
      end if                                                                    
c                                                                               
      sum = sum / dble(num_face_nodes)                                          
      dbar_now = gurson_cell_size + sum                                         
      if( abs(action) .eq. 1 ) then                                            
       dam_dbar_elems(1,elem_ptr) = dbar_now                                    
       dam_dbar_elems(2,elem_ptr) = zero                                        
      end if                                                                    
c                                                                               
      if( local_debug ) write(out,9000)  action, elem, sum, dbar_now                            
c                                                                               
      return                                                                    
c                                                                               
 9000 format(' > element release type 2. Cell height computation.',             
     & /,5x,'action: ',i1,10x,'element: ',i7,' avg. face displ: ',f10.6,        
     & ' D-bar now: ',f10.6 )                                                   
c                                                                               
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      function chk_killed                     *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 07/18/2019 rhd             *          
c     *                                                              *          
c     *     given an element, returns true if the                    *          
c     *     element had been killed, or false otherwise.             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      function chk_killed( elem ) result( t_chk_killed )                                      
c                                                                               
      use elem_extinct_data, only : dam_state                                   
      use damage_data                                                           
c                                                                               
      implicit none
c
      integer :: elem
      logical :: t_chk_killed    
c
      integer :: elem_ptr                                         
c                                                                               
      t_chk_killed = .false.                                                      
      elem_ptr =  dam_ptr(elem)                                               
      if( elem_ptr .ne. 0 ) then  ! element is killable                                             
       if( dam_state( elem_ptr) .gt. 0 )
     &       t_chk_killed = .true.                                               
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c
c
c     ****************************************************************
c     *                                                              *
c     *                      function chk_killed_blk                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 7/21/2019                  *
c     *                                                              *
c     *     for the specified block of elements, return a logical    *
c     *     vector indicating if each element is killed or active.   *
c     *     also a flag if the whole block is killed                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine chk_killed_blk( block_no, status_vec, block_killed )
      use global_data, only : out, elblks, mxvl ! old common.main
c
      use elem_extinct_data, only : dam_state, dam_blk_killed
      use damage_data, only : growth_by_kill, dam_ptr
c
      implicit none
      integer :: block_no
      logical :: status_vec(mxvl), block_killed
c      
      integer :: span, felem, i, n, nn, elem, elem_ptr
      logical, parameter :: check_consistency = .false.
c
      span  = elblks(0,block_no)
      felem = elblks(1,block_no) - 1
c
      block_killed = .false.
      status_vec   = .false.
      if( .not. growth_by_kill ) return
c
      if( .not. allocated(dam_blk_killed) ) then                                
        write(out,9100) 1                                                         
        call die_abort                                                          
      end if                                                                    

      if( dam_blk_killed(block_no) ) then
        block_killed = .true.
        status_vec   = .true.
        return
      end if
c
      if( .not. allocated(dam_state) ) then                                     
        write(out,9110) 2                                                         
        call die_abort                                                          
      end if 
c
      if( check_consistency ) then      
        n  = sizeof( dam_state )
        nn = sizeof( dam_ptr )
        do i = 1, span
          elem = felem + i
          if( elem > nn ) then
            write(out,9120) elem, nn
            call die_abort
          end if   
          elem_ptr = dam_ptr( elem )
          if( elem_ptr .ne. 0 ) then
            if( elem_ptr > n .or. elem_ptr <=0 ) then
               write(out,9130) elem, n
               call die_abort
            end if    
            if( dam_state(elem_ptr) .gt. 0 ) status_vec(i) = .true.
          end if
        end do
        return
      end if
c
      do i = 1, span  ! no consistency check
        elem = felem + i
        elem_ptr = dam_ptr( elem )
        if( elem_ptr .ne. 0 ) then
          if( dam_state(elem_ptr) .gt. 0 ) status_vec(i) = .true.
        end if
      end do
      return
c      
 9100 format(/,'FATAL ERROR: chk_killed_nlk. Contact WARP3D group',             
     &       /,'             Job terminated at ',i1,//)                         
 9110 format(/,'FATAL ERROR: chk_killed_nlk. Contact WARP3D group',             
     &       /,'             Job terminated at ',i1,//)   
 9120 format(/,'FATAL ERROR: chk_killed_nlk at 3. elem, nn:',2i8,
     &       /,'             Contact WARP3D group',             
     &       /,'             Job terminated....'//)                                             
 9130 format(/,'FATAL ERROR: chk_killed_nlk at 4. elem, n:',2i8,
     &       /,'             Contact WARP3D group',             
     &       /,'             Job terminated....'//)                                             
c
      end
c     ****************************************************************          
c     *                                                              *          
c     *                     subroutine damage_update                 *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 09/7/2019 rhd              *          
c     *                                                              *          
c     *        certain damage models need variables updated outside  *
c     *        element stress/history data (e.g. smcs )              *

c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine damage_update( now_step, now_iter )
c
      use global_data ! old common.main
      use elem_extinct_data, only : dam_state
      use damage_data, only :  crack_growth_type, dam_ptr                                                      
      implicit none
c
      integer :: now_step, now_iter
c                                       
      integer :: elem, elem_ptr, dowhat
      logical :: debug, local_debug, l1, process                                                    
      double precision :: d1, d2, d3, d4, d5, d6, d7
c        
      debug = .false.
      local_debug = .false. ! debug
      if( local_debug ) then      
         write(out,*) '... entered damage_update ...'                                                    
      end if                                                                    
c 
c              only damage model using SMCS requires this update
c              using data outside the stress/history data.
c
c              the Gurson model, for example, only needs the instantaneous
c              values stored in the element stress/history data to
c              determine if element needs to be killed, a step
c              size reduction is needed, values for printing, etc.
c      
      process = crack_growth_type .eq. 3
      if( .not. process ) return 
c
c              Loop for all elements. process killable elements.
c                                                                               
      do elem = 1, noelem
         elem_ptr = dam_ptr(elem)                                               
         if( elem_ptr .eq. 0 ) cycle     ! element not killable                                      
         if( dam_state(elem_ptr) .gt. 0 ) cycle ! element killed already    
         if( local_debug ) write(out,*) 'killable element is:',elem                           
c                                                                               
c                 calculate/update damage parameter(s) for element
c                 if this damage model requires it. 
c
c                 only SMCS needs this now to support plastic strain
c                 weighted averages.                   
c                                                                               
         select case( crack_growth_type )
         case( 3 ) !  SMCS   
          dowhat = 4
          call dam_param_3_get_values( elem, debug, d1, d2, d3, d4,            
     &                                 d5, d6, d7, dowhat, l1 )                
         end select 
c                         
      end do  ! over all model elements                                         
c      
      return    
      end
