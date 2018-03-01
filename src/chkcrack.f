c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine chkcrack                     *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 10/21/09 RHD (messages)    *          
c     *                                   11/24/2010 RHD             *          
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
      implicit integer (a-z)                                                    
c                                                                               
      logical debug                                                             
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
c                                                                               
c              Lower level routines are separated into other .f files by        
c              the type of crack growth.                                        
c                                                                               
c                                                                               
      select case (crack_growth_type)                                           
      case (0)                                                                  
         if (debug) write (out,*) '>>>> no crack growth specified.'             
      case (1)                                                                  
         if (debug) write (out,*) '>>>> use gurson crack growth.'               
         call chk_elem_kill( debug, step, iter )                                
      case (2)                                                                  
         if (debug) write (out,*) '>>>> use node release crack growth.'         
         call chk_node_release( debug, step, iter )                             
      case (3)                                                                  
         if (debug) write (out,*) '>>>> use smcs crack growth.'                 
         call chk_elem_kill( debug, step, iter )                                
      case (4)                                                                  
         call chk_elem_kill( debug, step, iter )                                
         return                                                                 
      end select                                                                
c                                                                               
      if ( debug ) then                                                         
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
c     *                   last modified : 07/30/2013 rhd             *          
c     *                                                              *          
c     *        This routine checks conditions to see if crack        *          
c     *        growth will occur by element extinction               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine chk_elem_kill( debug, step, iter )                             
      use global_data ! old common.main
      use elem_extinct_data, only : dam_blk_killed, dam_state,                  
     &                              kill_order_list, dam_print_list             
      use main_data, only : output_packets, packet_file_no                      
      use damage_data                                                           
      implicit integer (a-z)                                                    
c                                                                               
      parameter (max_local_list=200)                                            
      logical debug, blk_killed, killed_this_time, killed_found,                
     &   kill_the_elem, elems_left                                              
c                                                                               
      double precision                                                          
     &   dummy_arg, porosity, eps_plas, eps_crit, dummy_arg2,                   
     &   values(20), local_status(max_local_list,3), orig_poros,                
     &   ext_shape, ext_spacing                                                 
c                                                                               
      integer local_packet(max_local_list), local_pkt_type(5)                   
      logical fgm_cohes, ext_gurson, option_exponential,                        
     &        option_ppr, local_write_packets, do_kill_in_order,                
     &        found_exponential, found_ppr, found_cavit,                        
     &        option_cavit                                                      
      data local_pkt_type / 20, 0, 21, 10, 0 /                                  
c                                                                               
      if ( debug ) then                                                         
         if ( .not. no_killed_elems )                                           
     &        write(out,*) ' an element is dead.'                               
      end if                                                                    
c                                                                               
c              1) Loop for all elements. process killable elems                 
c                 ---------------------------------------------                 
c                                                                               
      killed_this_time = .false.                                                
      elems_left       = .false.                                                
      local_count = 0                                                           
      found_exponential = .false.                                               
      found_ppr         = .false.                                               
      found_cavit       = .false.                                               
c                                                                               
      do elem = 1, noelem                                                       
         elem_ptr = dam_ptr(elem)                                               
         if ( elem_ptr .eq. 0 ) cycle                                           
         if ( debug ) write(out,*) 'element is:',elem                           
c                                                                               
c                If element has been previously killed, update it's             
c                unloading state vector, and then skip the rest of the          
c                loop. we do this for traction-separation law as well           
c                since other routines look at dam_state.                        
c                                                                               
         if ( dam_state(elem_ptr) .gt. 0 ) then                                 
            if ( dam_state(elem_ptr) .le. max_dam_state )                       
     &           dam_state(elem_ptr) = dam_state(elem_ptr) + 1                  
            cycle                                                               
         end if                                                                 
c                                                                               
c                 calculate damage parameter(s) for element,                    
c                 then check it. if exceeded, set-up data structures for        
c                 internal force release using fixed no.                        
c                 of load steps or the traction-separation law.                 
c                                                                               
         if ( crack_growth_type .eq. 1 .or.                                     
     &        crack_growth_type .eq. 3 )                                        
     &       call dam_param( elem, kill_the_elem, debug, porosity,              
     &                     eps_plas, eps_crit, dummy_arg, dummy_arg2,           
     &                     ext_gurson, ext_shape, ext_spacing )                 
         if ( crack_growth_type .eq. 4 )                                        
     &       call dam_param_cohes( elem, kill_the_elem, debug,                  
     &                             values, 2 )                                  
         killed_this_time = killed_this_time .or. kill_the_elem                 
         if( .not. kill_the_elem ) then                                         
           elems_left = .true. ! elements remain to be killed                   
           cycle                                                                
         end if                                                                 
c                                                                               
c                 kill this element right now. record status for                
c                 subsequent packet output if neeeded. write usual              
c                 output about element death.                                   
c                                                                               
         local_count = local_count + 1                                          
         if( crack_growth_type .eq. 1 ) then                                    
             if( .not. ext_gurson ) write(out,9000) elem, porosity              
             if( ext_gurson ) write(out,9002) elem, porosity,                   
     &                        ext_spacing, ext_shape                            
         else if( crack_growth_type .eq. 3 ) then                               
             write(out,9010) elem, eps_plas, eps_crit                           
         else if( crack_growth_type .eq. 4 ) then                               
             cohes_type         = iprops(27,elem)                               
             option_exponential = cohes_type .eq. 4                             
             option_ppr         = cohes_type .eq. 6                             
             option_cavit       = cohes_type .eq. 7                             
             found_exponential = found_exponential .or.                         
     &                           option_exponential                             
             found_ppr   = found_ppr .or. option_ppr                            
             found_cavit = found_cavit .or. option_cavit                        
             count = 0                                                          
             if( found_exponential ) count = count + 1                          
             if( found_ppr ) count = count + 1                                  
             if( found_cavit ) count = count + 1                                
             if( count .gt. 1 ) then                                            
               write(iout,9400) elem                                            
               call die_gracefully                                              
             end if                                                             
             if( option_exponential )                                           
     &          write(out,9200) elem, values(6)/values(7), values(8)            
             if( option_ppr )   write(out,9220) elem                            
             if( option_cavit ) write(out,9230) elem                            
         end if                                                                 
c                                                                               
         if( output_packets ) then                                              
           if( local_count .gt. max_local_list ) then                           
              write(out,9300)                                                   
              call die_gracefully                                               
           end if                                                               
           if ( crack_growth_type .eq. 1 ) then                                 
                local_packet(local_count) = elem                                
                local_status(local_count,1) = porosity                          
                local_status(local_count,2) = ext_spacing                       
                local_status(local_count,3) = ext_shape                         
           else if ( crack_growth_type .eq. 3 ) then                            
                local_packet(local_count) = elem                                
                local_status(local_count,1) = eps_plas                          
                local_status(local_count,2) = eps_crit                          
           else if ( crack_growth_type .eq. 4 ) then                            
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
           end if                                                               
         end if                                                                 
c                                                                               
c                 if this is the first element to be killed, then               
c                 initialize needed variables.                                  
c                                                                               
         if( no_killed_elems ) then                                             
            call allocate_damage( 4 )                                           
            if( release_type .eq. 2 ) call allocate_damage( 5 )                 
            call dam_init2( debug )                                             
         end if                                                                 
c                                                                               
c                 null props and history for this newly killed element          
c                 so subsequent stiffness will be zero. zero stresses,          
c                 start the element force reduction process, constrain          
c                 completely uncoupled nodes if required.                       
c                                                                               
            num_elements_killed = num_elements_killed + 1                       
            call update_killed_energy( elem )                                   
            call kill_element( elem, debug )                                    
            call store_ifv( elem, elem_ptr, debug )                             
            call update_node_elecnt( elem, debug )                              
            if ( release_type .eq. 2 ) call growth_set_dbar( elem,              
     &                elem_ptr, debug, -1, dummy_arg )                          
      end do  ! over all model elements                                         
c                                                                               
c         end (1) loop over elements to check for new those now                 
c         exceeding the growth criterion                                        
c         ----------------------------------------------------                  
c                                                                               
c =======================================================================       
c                                                                               
c         (2) packet output for newly killed elements at                        
c             for gurson, extended-gurson, smcs, cohesive                       
c             -------------------------------------------                       
c                                                                               
      local_write_packets = output_packets .and. local_count .gt. 0             
      if( .not. local_write_packets ) go to 500                                 
      write(packet_file_no) local_pkt_type(crack_growth_type),                  
     &                      local_count, step, iter                             
      do elem = 1, local_count                                                  
        element = local_packet(elem)                                            
        if( crack_growth_type .eq. 1 ) then                                     
          orig_poros = props(26,element)                                        
          write(packet_file_no) element,orig_poros,                             
     &                  local_status(elem,1),                                   
     &                  local_status(elem,2),local_status(elem,3)               
        else if ( crack_growth_type .eq. 3 ) then                               
          write(packet_file_no) element, local_status(elem,1),                  
     &                                  local_status(elem,2)                    
        else if ( crack_growth_type .eq. 4 ) then                               
          write(packet_file_no) element, local_status(elem,1),                  
     &          local_status(elem,2), int(local_status(elem,3))                 
        end if                                                                  
      end do                                                                    
c                                                                               
c =======================================================================       
c                                                                               
c         (3) if no elements remain to be killed in the model, set              
c             global flag to indicate this condition.                           
c                                                                               
 500  continue                                                                  
      if( .not. elems_left ) all_elems_killed = .true.                          
c                                                                               
c =======================================================================       
c                                                                               
c         (4) if a killing order has been specified then make sure that         
c             no element "holes" have been left.                                
c                                                                               
      do_kill_in_order = kill_order .and. killed_this_time                      
      if( .not.  do_kill_in_order ) go to 600                                   
      if ( debug ) write (out,*) '>>> checking for holes....'                   
      killed_found = .false.                                                    
      do chk_kill = num_kill_order_list,1, -1                                   
        elem     = kill_order_list(chk_kill)                                    
        elem_ptr = dam_ptr(elem)                                                
        if ( .not. killed_found ) then                                          
          if ( dam_state(elem_ptr) .ne. 0 ) killed_found = .true.               
        else                                                                    
          if ( dam_state(elem_ptr) .eq. 0 ) then                                
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
c =======================================================================       
c                                                                               
c         (5) check the element blocks to see if all elements in                
c             them have been killed.  Set dam_blk_killed if so.                 
c                                                                               
 600  continue                                                                  
      if ( .not. killed_this_time ) go to 700                                   
      do blk = 1, nelblk                                                        
        if ( dam_blk_killed(blk) ) cycle                                        
        span  = elblks(0,blk)                                                   
        felem = elblks(1,blk)                                                   
        if ( iand( iprops(30,felem), 2 ) .eq. 0 ) cycle                         
        blk_killed = .true.                                                     
        do i = 1, span                                                          
         if ( dam_state(dam_ptr(felem+i-1)) .eq. 0 ) blk_killed =.false.        
        end do                                                                  
        if ( blk_killed ) dam_blk_killed(blk) = .true.                          
      end do                                                                    
c                                                                               
c =======================================================================       
c                                                                               
c         (6) re-check for free nodes if any elements have ever                 
c             been killed. Just makes triple sure that there are                
c             no free nodes.                                                    
c                                                                               
 700  continue                                                                  
      if ( .not. no_killed_elems ) call chk_free_nodes( debug )                 
c                                                                               
c =======================================================================       
c                                                                               
c         (7) if MPI, then all processors need to know what blocks              
c             and elements have been killed.  Send dam_blk_killed               
c             and dam_state.  Also have processors kill any elements            
c             which they own which should be killed.                            
c             this is a dummey for non-MPI                                      
c                                                                               
      call wmpi_send_growth ( killed_this_time )                                
c                                                                               
 9999 continue                                                                  
      return                                                                    
c                                                                               
 9000 format(/,'   >> element death invoked for element: ',i7,                  
     & '.   f: ',f6.3)                                                          
 9002 format(/,'   >> element death invoked for element: ',i7,                  
     & '.   f, W, X: ',3f8.3)                                                   
 9010 format(/,' >> element death option invoked before next step',             
     &       /,'    element: ',i7,' is now killed.',                            
     &       /,'    plastic strain: ',f8.5,' is > limit of: ',f8.5)             
 9100 format(/,' >> element death option invoked before next step',             
     &       /,'    element: ',i7,' is now killed to make crack',               
     &       /,'    face uniform.')                                             
 9200 format(/,'   >> element death invoked for element: ',i7,                  
     & '.   Deff/Dpeak: ',f5.2,' Teff/Tpeak: ',f5.2)                            
 9220 format(/,'   >> element death invoked for element: ',i7,                  
     & '.  PPR cohesive option')                                                
 9230 format(/,'   >> element death invoked for element: ',i7,                  
     & '.  cavit cohesive option')                                              
 9300 format(/,'FATAL ERROR: list length exceeded in chk_elem_kill',            
     & /,      '             job aborted.' )                                    
 9400 format(/,'FATAL ERROR: mixed cohesive options not allowed',               
     & /,      '             at present in WARP3D with crack growth.'           
     & /,      '             found while processing element: ',i7,              
     & /,      '             job aborted.' )                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_param                    *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 06/11/03  rhd              *          
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
c                                                                               
c                                                                               
      use main_data,       only : elems_to_blocks                               
      use elem_block_data, only : history_blocks, eps_n_blocks,                 
     &                            urcs_n_blocks, history_blk_list               
c                                                                               
      use damage_data                                                           
      implicit integer (a-z)                                                    
c                                                                               
c              parameter declarations                                           
c                                                                               
      logical kill_now, debug, ext_gurson                                       
      double precision                                                          
     &     porosity, eps_plas, eps_crit, sig_mean, sig_mises,                   
     &     ext_shape, ext_spacing                                               
      double precision,                                                         
     &  dimension(:), pointer :: history, urcs_n, eps_n                         
c                                                                               
c              local declarations                                               
c                                                                               
      double precision                                                          
     &     zero, one,                                                           
     &     third, six, iroot2,                                                  
     &     sig_xx, sig_yy, sig_zz,                                              
     &     sig_xy, sig_xz, sig_yz,eps_xx, eps_yy, eps_zz,                       
     &     gam_xy, gam_xz, gam_yz, smcs,  fpngp                                 
      data zero, one / 0.0, 1.0 /                                               
      data third / 0.333333333 /                                                
      data iroot2, six / 0.70711, 6.0 /                                         
c                                                                               
c                                                                               
      ext_gurson  = .false.                                                     
      ext_shape   = one                                                         
      ext_spacing = zero                                                        
c                                                                               
      go to ( 100, 200, 300, 400 ), crack_growth_type                           
c                                                                               
c                                                                               
c           1.0 gurson criterion - the material model 3 is the standard         
c                                Gurson model, 6 is the extended                
c                                Gurson model. call model dependent routine     
c                                to assess element status for killing and       
c                                to compute parameters.                         
c                                                                               
 100  continue                                                                  
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
         write(out,9100)                                                        
         call die_abort                                                         
      end if                                                                    
      return                                                                    
c                                                                               
c           2.0  growth type 2 is by ctoa - not processed here                  
c                                                                               
 200  continue                                                                  
      write(out,*) ' '                                                          
      write(out,*) '>>> FATAL ERROR: dam_param'                                 
      write(out,*) '    job terminated'                                         
      call die_gracefully                                                       
      stop                                                                      
c                                                                               
c           3.0   smcs - failure based on a stress modified critical strain     
c                    stress/strain vectors are (x,y,z,xy,yz,xz).                
c                                                                               
c             note:  we compute the mises and mean stress at each gp            
c                    then average the gp values. we cannot "average"            
c                    gp stresses in a geo. nonlin. analysis since the           
c                    stresses available now are those in "unrotated"            
c                    coordinates. but since mean stress and mises               
c                    are invariants, we can compute and average them.           
c                                                                               
c                    get plastic strain from history data at gp                 
c                                                                               
 300  continue                                                                  
         sig_mean    = zero                                                     
         sig_mises   = zero                                                     
         eps_plas    = zero                                                     
         ngp         = iprops(6,elem)                                           
         blk         = elems_to_blocks(elem,1)                                  
         rel_elem    = elems_to_blocks(elem,2)                                  
         hist_size   = history_blk_list(blk)                                    
         hoffset     = (rel_elem-1)*hist_size*ngp + 1                           
         epsoffset   = (rel_elem-1)*nstr*ngp                                    
         sigoffset   = (rel_elem-1)*nstrs*ngp                                   
         urcs_n      => urcs_n_blocks(blk)%ptr                                  
         eps_n       => eps_n_blocks(blk)%ptr                                   
         history     => history_blocks(blk)%ptr                                 
c                                                                               
         do gp = 1, ngp                                                         
           sig_xx = urcs_n(sigoffset+1)                                         
           sig_yy = urcs_n(sigoffset+2)                                         
           sig_zz = urcs_n(sigoffset+3)                                         
           sig_xy = urcs_n(sigoffset+4)                                         
           sig_yz = urcs_n(sigoffset+5)                                         
           sig_xz = urcs_n(sigoffset+6)                                         
           eps_xx = eps_n(epsoffset+1)                                          
           eps_yy = eps_n(epsoffset+2)                                          
           eps_zz = eps_n(epsoffset+3)                                          
           gam_xy = eps_n(epsoffset+4)                                          
           gam_yz = eps_n(epsoffset+5)                                          
           gam_xz = eps_n(epsoffset+6)                                          
           sig_mean  = sig_mean + (sig_xx + sig_yy + sig_zz )                   
           sig_mises = sig_mises +                                              
     &                 sqrt( (sig_xx-sig_yy)**2 + (sig_yy-sig_zz)**2 +          
     &                       (sig_xx-sig_zz)**2 + six*( sig_xy**2 +             
     &                       sig_yz**2 + sig_xz**2 ) ) * iroot2                 
           eps_plas  = eps_plas + history(hoffset+0+(gp-1)*hist_size)           
           epsoffset = epsoffset + nstr                                         
           sigoffset = sigoffset + nstrs                                        
         end do                                                                 
c                                                                               
         fpngp     = ngp                                                        
         sig_mean  = sig_mean * third / fpngp                                   
         sig_mises = sig_mises / fpngp                                          
         eps_plas  = eps_plas / fpngp                                           
         eps_crit  = 100.0                                                      
         if ( sig_mises .gt. 1.0e-06 )                                          
     &       eps_crit = smcs_alpha * exp(-smcs_beta*sig_mean/sig_mises)         
         smcs     = eps_plas - eps_crit                                         
         kill_now = smcs .ge. zero                                              
         return                                                                 
c                                                                               
c           4.0 growth type 4 is by cohesive - not processed here               
c                                                                               
c                                                                               
 400  continue                                                                  
      write(out,*) ' '                                                          
      write(out,*) '>>> FATAL ERROR: dam_param (4)'                             
      write(out,*) '    job terminated'                                         
      call die_gracefully                                                       
      stop                                                                      
c                                                                               
 9000 format(1x,'> smcs failure criterion check. element: ',i7,                 
     & /,5x,'ngp, ym, locsig, loceps: ',i2,f12.0,2i8,                           
     & /,5x,'smcs, kill_now: ',                                                 
     & /,5x,f15.6,5x,l1,                                                        
     & /,5x,'mises, mean stress: ',2f15.6,                                      
     & /,5x,'eps_eff, eps_plastic: ',2f15.9 )                                   
 9100 format('>>>> FATAL ERROR: invalid model no. in dam_param')                
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_init2                    *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 06/24/94                   *          
c     *                   last modified : 07/28/95                   *          
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
      implicit integer (a-z)                                                    
c                                                                               
      logical debug, local_debug                                                
      double precision                                                          
     &          node_coord(3), coord, plane_tol                                 
      integer, dimension(:,:), pointer :: cdest                                 
      data plane_tol, local_debug  / 0.01, .false. /                            
c                                                                               
      if ( debug ) write(out,*) '>>>> in dam_init2'                             
c                                                                               
      no_killed_elems = .false.                                                 
      do i = 1, nonode                                                          
         dam_node_elecnt(i) = inverse_incidences(i)%element_count               
      end do                                                                    
      if ( debug ) then                                                         
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
      if ( release_type .ne. 2 ) go to 999                                      
      do elem  = 1, noelem                                                      
        elem_ptr = dam_ptr(elem)                                                
        if ( elem_ptr .eq. 0 ) cycle                                            
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
          if ( abs(coord-crack_plane_coord) .le. plane_tol*                     
     &         gurson_cell_size ) cycle                                         
          face_count = face_count + 1                                           
          dam_face_nodes(face_count,elem_ptr) = snode                           
        end do                                                                  
        if ( face_count .ne. 4 ) then                                           
             write(out,*) 'FATAL error 1 in dam_init2'                          
             write(out,*) 'invalid plane definition for growth'                 
             write(out,*)  elem, elem_ptr, num_enodes,face_count                
             call die_gracefully                                                
             stop                                                               
        end if                                                                  
       end do                                                                   
c                                                                               
       if ( local_debug ) then                                                  
        write(out,*) '> element release type 2. face node table'                
        do elem  = 1, noelem                                                    
         elem_ptr = dam_ptr(elem)                                               
         if ( elem_ptr .eq. 0 ) cycle                                           
         write(out,9050) elem, dam_face_nodes(1,elem_ptr),                      
     &          dam_face_nodes(2,elem_ptr), dam_face_nodes(3,elem_ptr),         
     &          dam_face_nodes(4,elem_ptr)                                      
        end do                                                                  
       end if                                                                   
c                                                                               
 999  continue                                                                  
      if ( debug ) write(out,*) '<<<< leaving dam_init2'                        
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
c     *                   last modified : 06/24/94                   *          
c     *                                                              *          
c     *        This routine checks each node in the structure to     *          
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
      use main_data,         only : cnstrn, cnstrn_in,                          
     &                              inverse_incidences                          
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      logical debug                                                             
c                                                                               
      double precision                                                          
     &   zero                                                                   
      data zero /0.0/                                                           
      if ( debug ) write(out,*) '>>> in chk_free_nodes <<<'                     
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
         if ( dam_node_elecnt(node) .eq. 0 ) then                               
            num_dof = iprops(4,inverse_incidences(node)%element_list(1))        
            if ( debug ) then                                                   
               write(out,*) 'free node ',node                                   
               write (out,9010) node,cnstrn(dstmap(node)),                      
     &              cnstrn(dstmap(node)+1),cnstrn(dstmap(node)+2)               
            end if                                                              
            do j = 1, num_dof                                                   
               glbal_dof            = dstmap(node) + j-1                        
               cnstrn(glbal_dof)    = zero                                      
               cnstrn_in(glbal_dof) = zero                                      
               if ( cstmap(glbal_dof) .eq. 0 ) then                             
                  cstmap(csttail) = glbal_dof                                   
                  csttail         = glbal_dof                                   
               end if                                                           
            end do                                                              
            cstmap(csttail) = -1                                                
         end if                                                                 
      end do                                                                    
c                                                                               
 9999 continue                                                                  
      if ( debug ) write(out,*) '>>> in chk_free_nodes <<<'                     
      return                                                                    
 9010 format('  old constraints for node ',i7,' :',3(1x,e14.6))                 
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dam_print                    *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 10/14/94                   *          
c     *                                                              *          
c     *     This routine calls routines to print out the status of   *          
c     *     killed elements or released nodes at the beginning of a  *          
c     *     load step.                                               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine dam_print( step, iter )                                        
      use global_data ! old common.main
      use damage_data                                                           
      implicit integer (a-z)                                                    
c                                                                               
c           check to make sure crack growth is on and printing is               
c           specified                                                           
c                                                                               
      if ( print_status ) then                                                  
         if ( crack_growth_type .eq. 1 ) then                                   
            call dam_print_elem1( step, iter )                                  
         else if ( crack_growth_type .eq. 2 ) then                              
           if (const_front) then                                                
               call dam_print_front( step, iter)                                
           else                                                                 
               call dam_print_node( step, iter )                                
           endif                                                                
         else if ( crack_growth_type .eq. 3 ) then                              
            call dam_print_elem3( step, iter )                                  
         else if ( crack_growth_type .eq. 4 ) then                              
            call dam_print_elem4( step, iter )                                  
         end if                                                                 
      end if                                                                    
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
c                                                                               
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
c     *                   last modified : 07/28/95                   *          
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
      implicit integer (a-z)                                                    
      logical debug, local_debug                                                
      double precision                                                          
     &   zero, sum, node_displ(3), dbar_now                                     
      data local_debug, zero, num_face_nodes / .false., 0.0, 4 /                
c                                                                               
      if ( debug ) write (out,*) '>>>> in growth_set_dbar'                      
c                                                                               
      sum = zero                                                                
      if ( action .lt. 0 ) then                                                 
        do i = 1, num_face_nodes                                                
          snode         = dam_face_nodes(i,elem_ptr)                            
          node_displ(1) = u(dstmap(snode)+0)                                    
          node_displ(2) = u(dstmap(snode)+1)                                    
          node_displ(3) = u(dstmap(snode)+2)                                    
          sum           = sum + node_displ(crk_pln_normal_idx)                  
        end do                                                                  
      end if                                                                    
      if ( action .gt. 0 ) then                                                 
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
      if ( abs(action) .eq. 1 ) then                                            
       dam_dbar_elems(1,elem_ptr) = dbar_now                                    
       dam_dbar_elems(2,elem_ptr) = zero                                        
      end if                                                                    
c                                                                               
      if ( local_debug ) then                                                   
        write(out,9000)  action, elem, sum, dbar_now                            
      end if                                                                    
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
c     *                   last modified : 04/23/96                   *          
c     *                                                              *          
c     *     This function, given an element, returns true if the     *          
c     *     element had been killed, or false otherwise.             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      logical function chk_killed( elem )                                       
c                                                                               
      use elem_extinct_data, only : dam_state                                   
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      chk_killed = .false.                                                      
      elem_ptr =  dam_ptr( elem )                                               
      if ( elem_ptr .ne. 0 ) then                                               
         if ( dam_state( dam_ptr( elem_ptr ) ) .gt. 0 )                         
     &        chk_killed = .true.                                               
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      function chk_killed_blk                 *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 06/14/03                   *          
c     *                                                              *          
c     *     for the specified block of elements, return a logical    *          
c     *     vector indicating if each element is killed or active.   *          
c     *     also a flag if the whole block is killed                 *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine chk_killed_blk( block_no, status_vec, block_killed )           
      use global_data ! old common.main
c                                                                               
      use elem_extinct_data, only : dam_state, dam_blk_killed                   
      use damage_data                                                           
c                                                                               
      implicit integer (a-z)                                                    
c                                                                               
      logical status_vec(*), block_killed                                       
c                                                                               
      span  = elblks(0,block_no)                                                
      felem = elblks(1,block_no) - 1                                            
c                                                                               
      block_killed       = .false.                                              
      status_vec(1:span) = .false.                                              
      if( .not. growth_by_kill ) return                                         
c                                                                               
      if( .not. allocated(dam_blk_killed) ) then                                
        write(*,9100) 1                                                         
        call die_abort                                                          
      end if                                                                    
c                                                                               
      if( dam_blk_killed(block_no) ) then                                       
        block_killed       = .true.                                             
        status_vec(1:span) = .true.                                             
        return                                                                  
      end if                                                                    
c                                                                               
      if( .not. allocated(dam_state) ) then                                     
        write(*,9110) 2                                                         
        call die_abort                                                          
      end if                                                                    
                                                                                
      do i = 1, span                                                            
        elem = felem + i                                                        
        elem_ptr = dam_ptr( elem )                                              
        if( elem_ptr .ne. 0 ) then                                              
           if( dam_state( dam_ptr( elem_ptr ) ) .gt. 0 )                        
     &         status_vec(i) = .true.                                           
        end if                                                                  
      end do                                                                    
c                                                                               
      return                                                                    
 9100 format(/,'FATAL ERROR: chk_killed_nlk. Contact WARP3D group',             
     &       /,'             Job terminated at ',i1,//)                         
 9110 format(/,'FATAL ERROR: chk_killed_nlk. Contact WARP3D group',             
     &       /,'             Job terminated at ',i1,//)                         
      end                                                                       
