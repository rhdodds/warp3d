c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine chkcrack                     *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 12/6/20 rhd                *          
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
c                                      
      use global_data, only : out
      use damage_data, only : crack_growth_type
c                                                          
      implicit none                                                    
c                                                                               
      integer, intent(in) :: step, iter
      logical :: debug, ldebug                                                            
c                                                                               
      debug  = .false.    
      ldebug = .false.                                                       
c                                                                               
      if( ldebug ) write(out,*) '>>> entered chkcrack ...'                        
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
         if(ldebug) write (out,*) '>>>> use smcs crack growth.'                 
         call chk_elem_kill( debug, step, iter )                                
      case(4)    !  cohesive                                                               
         call chk_elem_kill( debug, step, iter )                                
         return                                                                 
      case(5)                                                                  
         if(debug) write (out,*) '>>>> user specified kill element.'      
         call chk_elem_kill( debug, step, iter )                                
      end select                                                                
c                                                                               
      if( debug ) call dam_debug    
      if( ldebug ) then 
        write(out,*) '>>> leaving chkcrack ...'                       
        write(out,*) ' '
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
c     *                   last modified : 9/15/21 rhd                *          
c     *                                                              *          
c     *        checks conditions to see if crack growth will occur   *          
c     *        by element extinction. if so start the deletion.      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine chk_elem_kill( debug, step, iter )
c                             
      use global_data, only : out, noelem, current_load_time_step, 
     &                        iprops, props, nelblk, elblks, num_error,
     &                        mxgp
      use dam_param_code, only : dam_param, growth_set_dbar
      use elem_extinct_data, only : dam_state, smcs_start_kill_step,                 
     &                              kill_order_list, dam_print_list,
     &                              smcs_d_values,
     &                              smcs_eps_plas_at_death,
     &                              smcs_stress_at_death,
     &                              Oddy_metrics
      use main_data, only : output_packets, packet_file_no
      use damage_data, only : no_killed_elems, dam_ptr, max_dam_state,
     &     crack_growth_type, num_user_kill_elems, release_type, 
     &     user_kill_list_now, all_elems_killed, killed_element_limit,
     &     num_elements_killed, kill_order, smcs_deleted_list_file_flag, 
     &     smcs_deleted_list_file_name, num_kill_order_list, smcs_type,
     &     stop_killed_elist_length, deleted_elist_to_stop, 
     &     distortion_plastic_limit, smcs_type_5_tp_critical,
     &     use_mesh_regularization,
     &     tol_regular =>
     &     tolerance_mesh_regularization, use_distortion_metric, 
     &     Oddy_critical_ratio    
      use constants                   
c                                                         
      implicit none
c
      integer, intent(in) :: step, iter
      logical, intent(in) :: debug                                                                                                                                   
c                                                                               
c               local declarations                                               
c
      integer, parameter :: max_local_list=500 ! max elements that
c                                                can be killed this call    
c                                       
      integer :: local_packet(max_local_list), elem, elem_ptr,
     &           local_count, blk,  cohes_type, count, 
     &           element, felem, i, span
      integer, parameter :: local_pkt_type(5) = (/ 20, 0, 21, 10, 0 / )                   
      logical :: blk_killed, killed_this_time, distortion_killed,           
     &           kill_the_elem, ldummy1, ldummy2,
     &           local_debug, stop_solution, standard_kill, 
     &           elements_have_been_killed,  
     &           fgm_cohes, ext_gurson, option_exponential,                        
     &           option_ppr, local_write_packets, do_kill_in_order,                
     &           found_exponential, found_ppr, found_cavit,                        
     &           option_cavit, elem_deleted                                                     
      logical, external :: scan_entry_in_list, chkcrack_is_elem_deleted                                                                                                                       
      double precision ::                                                        
     &   porosity, eps_plas, eps_crit, d_old, d_new,           
     &   values(20), local_status(max_local_list,3), orig_poros,                
     &   ext_shape, ext_spacing, avg_triaxiality, avg_zeta,
     &   avg_bar_theta, tear_param, sig_mises, mises_at_death,
     &   sig_1, sig_1_at_death, sig_mean
c        
      local_debug = .false.
      if( local_debug ) then      
         write(out,*) '... entered chk_elem_kill ...'                                                    
         write(out,*) '    ... no_killed_elems:',no_killed_elems
      end if                                                                    
c                                                                               
c              1) Loop for all elements. process killable elements.
c                 start kill as required.           
c
c                 standard kill: immediately delete element. no
c                                [Ke], stresses, properties, histories.
c                                element forces save for prescribed
c                                reduction to zero in future steps.
c                 use_mesh_regularization: element stiffness, stresses
c                                gradually released to zero following
c                                a user-defined damage function.      
c                                                                               
      killed_this_time  = .false.                                                
      local_count       = 0                                                           
      found_exponential = .false.                                               
      found_ppr         = .false.                                               
      found_cavit       = .false.    
      stop_solution     = .false.   
      standard_kill     = .not. use_mesh_regularization                                       
c     
      do elem = 1, noelem
c          
         elem_ptr = dam_ptr(elem)                                               
         if( elem_ptr == 0 ) cycle     ! element not killable                                      
         if( local_debug ) write(out,*) 'killable element, dam state:',
     &                                  elem, dam_state(elem_ptr)
c                                                                               
c                If we previously started killing element, update it's             
c                unloading state vector or damage parameter, 
c                and then skip the rest of the          
c                loop. we do this for traction-separation law as well           
c                since other routines look at dam_state.                        
c              
         if( dam_state(elem_ptr) > 0 ) then ! already being killed
            if( use_mesh_regularization ) then
               call chk_elem_kill_upmr
               cycle 
            end if
            if( standard_kill ) then
               if( dam_state(elem_ptr) <=  max_dam_state )                        
     &             dam_state(elem_ptr) = dam_state(elem_ptr) + 1                  
               cycle   
            end if
            write(out,9100) elem
            call die_abort
         end if                                                                 
c                                                                               
c                 calculate damage parameter(s) for element,                    
c                 then check it. if exceeded, set-up data structures for        
c                 internal force release using fixed no.                        
c                 of load steps or the traction-separation law.                 
c                                                                               
         select case( crack_growth_type )
         case( 1 )  ! GTN
            call dam_param(
     &         elem=elem, kill_now=kill_the_elem, debug=debug,
     &         porosity=porosity, eps_plas=eps_plas, eps_crit=eps_crit,
     &         sig_mean=sig_mean, sig_mises=sig_mises, 
     &         ext_gurson=ext_gurson, ext_shape=ext_shape,
     &         ext_spacing=ext_spacing )
            distortion_killed = .false.
            if( .not. kill_the_elem ) then ! distortion may kill it
               if( use_distortion_metric ) call chk_gt_kill_distortion
            end if      
         case( 3 )  ! SMCS   
            call dam_param(
     &         elem=elem, kill_now=kill_the_elem, debug=debug,              
     &         eps_plas=eps_plas, eps_crit=eps_crit, 
     &         sig_mises=sig_mises,          
     &         avg_triaxiality=avg_triaxiality, avg_zeta=avg_zeta, 
     &         avg_bar_theta=avg_bar_theta, tear_param=tear_param,
     &         max_princ_stress=sig_1 )
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
         if( .not. kill_the_elem ) cycle                                                                
c                                                                               
c                  for standard element deletion, the element is 
c                  immediately removed from the mesh by zeroing 
c                  strains, stresses, history, stiffness. forces 
c                  exerted by element on model nodes are gradually
c                  relaxed to zero using one of several methods.
c
c                  for mesh regularization, the element stiffness and
c                  internal forces are degraded towards zero with 
c                  increasing deformation based on a damage 
c                  parameter (d). When d > (1.0-tol), the element
c                  stiffness, strains, stresses, .. are zeroed and no
c                  further element computations occur
c
         local_count = local_count + 1                                          
         if( crack_growth_type == 1 ) call chk_output_kill_messages_1
         if( crack_growth_type == 3 ) call chk_output_kill_messages_3
         if( crack_growth_type == 4 ) call chk_output_kill_messages_4
         if( crack_growth_type == 5 ) call chk_output_kill_messages_5
c
         if( output_packets ) call chk_elem_kill_make_packets
c                                                                               
c                 initialize global variables for 1st killed element.
c                 for standard kill process (no damage tracking),  
c                 null props and history for this newly killed element          
c                 so subsequent stiffness will be zero. zero stresses,          
c                 start the element force reduction process, constrain          
c                 completely uncoupled nodes if required.                       
c     
         if( no_killed_elems ) then     
            call allocate_damage( 4 ) ! setup to find newly free nodes                                         
            if( release_type .eq. 2 ) call allocate_damage( 5 )                 
            call chk_dam_init2( debug )                                             
         end if                                                                 
         num_elements_killed = num_elements_killed + 1                       
         if( standard_kill ) then
           call chk_update_killed_energy( elem )                                   
           call chk_kill_element_now( elem, debug )                                    
           call chk_store_ifv( elem, elem_ptr, debug )                             
           call chk_update_node_elecnt( elem, debug )                              
           if( release_type .eq. 2 ) call growth_set_dbar( elem,              
     &             elem_ptr, debug, -1 ) 
         end if
         if( use_mesh_regularization ) then
            mises_at_death = sig_mises
            sig_1_at_death = sig_1
            smcs_eps_plas_at_death(elem_ptr) = eps_plas
            smcs_stress_at_death(elem_ptr) = sig_1_at_death
            smcs_start_kill_step(elem_ptr) = current_load_time_step
            dam_state(elem_ptr) = 1
            call chk_get_d( elem, zero, sig_1_at_death, d_new, out )
            d_old = smcs_d_values(elem_ptr)
            if( d_new < d_old .and. iter > 1) then
              write(out,9020) elem, d_old, d_new
              num_error = num_error + 1 ! stop at next compute cmd
            end if
            smcs_d_values(elem_ptr) = d_new
         end if
c
c                 the user may have provided a list if elements.
c                 if this newly killed element is is the list,
c                 stop solution at cleanup in this routine.
c
         if( stop_killed_elist_length > 0 ) then
           if( scan_entry_in_list( elem, deleted_elist_to_stop, 
     &         stop_killed_elist_length ) ) stop_solution = .true.
         end if

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
c         (3) are all elements in model now killed ?
c             std kill means all killable elements have dam state >0
c             regularization kill means damage > one
c          
      all_elems_killed = .true.
      do elem = 1, noelem
         elem_deleted = chkcrack_is_elem_deleted( elem )
         if( elem_deleted ) cycle
         all_elems_killed = .false.
         exit
      end do
c                                                                               
c         (4) if a killing order has been specified then make sure that         
c             no element "holes" have been left.                                
c 
      if( local_debug )write(out,*) '.. kill_order, killed this time:',
     &      kill_order, killed_this_time                                                     
      do_kill_in_order = kill_order .and. killed_this_time                      
      if( do_kill_in_order ) call chk_elem_kill_in_order 
c                                                                               
c                                                                               
c         (5) deprecated code
c                                                                               
c         (6) re-check for free nodes if any elements have ever                 
c             been killed. Just makes triple sure that there are                
c             no free nodes.                                                    
c              
      elements_have_been_killed =  .not. no_killed_elems                                                                 
      if( elements_have_been_killed ) call chk_free_nodes( debug )                 
c                                                                               
c         (7) if MPI, then all processors need to know what blocks              
c             and elements have been killed.  send dam_state,...          
c             Also have processors kill any elements            
c             which they own which should be killed.                            
c             this is a dummy for non-MPI                                      
c                                                                               
      call wmpi_send_growth( killed_this_time )                                
c                                                                               
c         (8) if limit on number of killed elements is exceeded,
c             stop the solution.
c
      if( num_elements_killed >= killed_element_limit ) then
          write(out,9000) killed_element_limit
          call store( ' ','kill_limit_restart.db', ldummy1, ldummy2 )
          call warp3d_normal_stop
      end if
c                                                                               
c         (9) the user may have provided a list of elements
c             to stop the analysis if <any> of them are killed
c             stop the solution.
c
      if( stop_solution ) then
          write(out,9010) 
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
 9010 format(/,1x,">>>>> An element just killed appears in the user",
     & " defined list to stop the solution.",
     &       /,7x,"Restart file: kill_limit_restart.db written.",
     &       /,7x,"Job ended normally."// )
 9020 format(/,1x,">>>>> Inconsistent condition for element: ",i8,
     &       /,1x,"      new damage paramter < old damage parameter",
     &       /,1x,"      values: ", 2f8.4,
     &       /,1x,"      in routine chk_elem_kill", 
     &       /,1x,"      Execution will stop at next Compute cmd",
     &       /)
 9100 format('>> FATAL ERROR: routine ',
     & 'chk_elem_kill. @ 1',
     & 10x,'inconsistent condition. for element:',i8,
     & ' job terminated.'//)
c                                                                               
c
      contains
c     ========
c
c
c     ****************************************************************
c     *                                                              *
c     *          internal subroutine chk_gt_kill_distortion          *
c     *                                                              *
c     *                   last modified : 8/29/21 rhd                *
c     *                                                              *
c     *     Gurson element did not fail by critical porosity.        *
c     *     Check if it fails by distortion metrics                  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine chk_gt_kill_distortion
c
      implicit none
c
      logical :: local_debug, kill_by_Oddy, kill_immediately,
     &           kill_by_plastic_limit
      double precision :: total_plastic_strain, max_Oddy
c
      local_debug = elem == -1
      if( local_debug ) write(out,9000) elem
c
c              current plastic strain
c
      total_plastic_strain = eps_plas
      if( eps_plas < zero ) total_plastic_strain = zero ! roundoff
c
c               could also delete by distortion
c               limit or plastic strain limit.
c
      max_Oddy = Oddy_metrics(elem_ptr,2) / Oddy_metrics(elem_ptr,1)
      kill_by_Oddy = max_Oddy >= Oddy_critical_ratio 
      kill_by_plastic_limit = total_plastic_strain >= 
     &                        distortion_plastic_limit
      kill_immediately = kill_by_Oddy  .or. kill_by_plastic_limit
      if( local_debug ) then
        write(out,9020) total_plastic_strain, max_Oddy
      end if  
      if( .not. kill_immediately ) return
c      
      kill_the_elem = .true.
      distortion_killed = .true. 
      return
c
 9000 format(5x,"... entered chk_xxx: ",i8)
 9010 format(10x,'eps_plas ',f12.8)
 9020 format(10x,'eps_plas, gt_max_Oddy: ',f8.5,f7.2)
 9200 format(/,' @2  >> element removed by distortion metrics: ',i7,
     & '  max Oddy ratio: ',f7.2,2x,'plastic strain: ',
     & f4.2,2x,a)
c
      end subroutine chk_gt_kill_distortion
cc     ****************************************************************
c     *                                                              *
c     *          internal subroutine chk_elem_kill_upmr              *
c     *                                                              *
c     *                   last modified : 8/28/21 rhd                *
c     *                                                              *
c     *     element in deletion via SMCS using mesh regularization.  *
c     *     update damage parameter. if fully deleted, updated       *
c     *     other data structures                                    *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine chk_elem_kill_upmr
c
      implicit none
c
c              element is in progress of being deleted. mesh
c              regularization is being used. Update the damage
c              parameter "d" for the element. It was zeroed
c              upon triggering deletion. When d >=1, the
c              element is to be fully deleted as in the standard
c              killing process.
c
      logical :: local_debug, kill_by_Oddy, kill_by_d_value, 
     &           kill_immediately, kill_by_plastic_limit
      double precision ::  avg_eps_plas, avg_sig_mises, d_now, d_old,
     &                     deps_plas, stress_at_death, avg_sig_1,
     &                     max_Oddy, total_plastic_strain, max_eps_plas
      character(len=50) :: distortion_msg
c
      local_debug = elem == -1
      if( local_debug ) write(out,9000) elem
c
      if( smcs_d_values(elem_ptr) > tol_regular ) return ! already fully killed
c
c              current plastic strain and mises stress averages
c              for element. avg_sig_1 is *max* principal stress
c
      call chk_get_ele_vals( elem, avg_eps_plas, avg_sig_mises, 
     &                       avg_sig_1, max_eps_plas )
      deps_plas = avg_eps_plas - smcs_eps_plas_at_death(elem_ptr)
      if( local_debug ) then
        write(out,9010) avg_eps_plas, deps_plas, 
     &                  smcs_eps_plas_at_death(elem_ptr)
      end if
      if( deps_plas < zero ) deps_plas = zero ! round-off @ just killed
c
c               get current damage parameter value "d". decide to
c               fully delete element. could also delete by distortion
c               limit or plastic strain limit.
c
      stress_at_death = smcs_stress_at_death(elem_ptr)
      call chk_get_d( elem, deps_plas, stress_at_death, d_now, out )
      d_old = smcs_d_values(elem_ptr)
      if( d_now < d_old .and. iter > 1 ) then
         write(out,9040) elem, d_old, d_now
         num_error = num_error + 1 ! stop at next compute cmd
      end if
      smcs_d_values(elem_ptr) = d_now
      if( local_debug ) write(out,9020) d_old, d_now, stress_at_death
c
      if( use_distortion_metric ) then
         max_Oddy = Oddy_metrics(elem_ptr,2) / Oddy_metrics(elem_ptr,1)
         total_plastic_strain = max_eps_plas
      else
         max_Oddy = zero
         total_plastic_strain = -one ! actual value ust be larger
      end if 
c
      kill_by_Oddy = max_Oddy >= Oddy_critical_ratio 
      kill_by_d_value = d_now > tol_regular
      kill_by_plastic_limit = total_plastic_strain >= 
     &                        distortion_plastic_limit
c      
      kill_immediately = kill_by_Oddy  .or. kill_by_d_value .or. 
     &                   kill_by_plastic_limit
      if( .not. kill_immediately ) return
      smcs_d_values(elem_ptr) = hundred
c
c               fully kill the element. 
c               data structures like a standard (immediate) deletion
c               process.
c                   - energy ??
c                   - zero properties, stresses, history
c                   - decreases by 1 the count of elements attached to
c                     the now fully deleted element
c                   - if a node now has no elements attached, add
c                     constraints to prevent rbm
c
      if( local_debug )  write(out,9030)
      call chk_update_killed_energy( elem )     
      call chk_kill_element_now( elem, debug )                                    
      call chk_update_node_elecnt( elem, debug )  
      call chk_free_nodes( .false. ) ! debug flag  
      if( use_distortion_metric ) then
        distortion_msg = ' '
        if( kill_by_Oddy ) distortion_msg = '** Oddy exceeds limit'
        if( kill_by_plastic_limit ) distortion_msg =
     &                '** eps plastic exceeds limit'
        if( kill_by_Oddy .and. kill_by_plastic_limit )
     &     distortion_msg = '** Oddy & eps plastic exceed limits'
        write(out,9200) elem, d_now, max_Oddy, total_plastic_strain, 
     &                  distortion_msg(1:len(distortion_msg))
      else
        write(out,9210) elem, d_now
      end if  
c               
      return
c
 9000 format(5x"... entered chk_elem_kill_upmr. elem: ",i8)
 9010 format(10x,'avg eps pls, deps_pls, eps_pls @ death: ',3f12.8)
 
 9020 format(10x,'d_old: ',f6.3,'d_now: ',f6.3,' sig_1 @ death:',f8.3)
 9030 format(10x,'ready to leave. fully deleting element')
 9040 format(/,1x,">>>>> Inconsistent condition for element: ",i8,
     &       /,1x,"      new damage paramter < old damage parameter",
     &       /,1x,"      values: ", 2f8.4,
     &       /,1x,"      in routine chk_elem_kill_upmr", 
     &       /,1x,"      Execution will stop at next Compute cmd",
     &       /)
 9100 format('>> FATAL ERROR: routine ',
     & 'chk_elem_kill_update_mesh_regularization.',
     & 10x,'inconsistent condition. for element:',i8,
     & ' job terminated.'//)
 9200 format(/,' >> element removed by regularization: ',i7,
     & '  d: ', f5.2,'  max Oddy ratio: ',f7.2,2x,'plastic strain: ',
     & f4.2,2x,a)
 9210 format(/,' >> element removed by regularization: ',i7,
     & '  d: ', f5.2 )
c
      end subroutine chk_elem_kill_upmr
c     ****************************************************************
c     *                                                              *
c     *          internal subroutine chk_get_ele_vals                *
c     *                                                              *
c     *                   last modified : 8/26/21 rhd                *
c     *                                                              *
c     *         service routine to get strain-stress values for a    *
c     *         single element                                       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine chk_get_ele_vals( elem, avg_eps_plas, avg_sig_mises,
     &                             avg_sig_1, max_eps_plas  )
c
      implicit none
c
c              return current values of average plastic strain and 
c              average mises stress for element
c
      integer, intent(in) :: elem
      double precision, intent(out) :: avg_eps_plas, avg_sig_mises,
     &                                 avg_sig_1, max_eps_plas
c
      integer :: elem_type, ngp
      logical :: threed_solid_elem, l1, l2, l3, l4, l5, l6, l7, l8,
     &           l9, l10
      logical :: ldebug 
      double precision :: workvec(27)
c
      ldebug = .false. ! elem .eq. 4 .and. step > 57
      if( ldebug ) 
     &    write(out,9000) elem
c
      elem_type = iprops(1,elem)
      ngp       = iprops(6,elem)    
c
      if( step < 2 ) then  ! sanity check. step passed to chkcrack
        call set_element_type( elem_type, threed_solid_elem, l1, l2, l3,
     &                       l4, l5, l6, l7, l8, l9, l10 )  
        if( .not. threed_solid_elem ) then
          write(out,9100) elem
          call die_abort
        end if
      end if
c 
      call mm_return_values( "plastic_strain", elem, workvec, ngp )
      avg_eps_plas = sum( workvec(1:ngp) ) * ( one / dble( ngp ) )
      max_eps_plas = maxval( workvec(1:ngp) )
c
      call mm_return_values( "avg_mises", elem, workvec, ngp )
      avg_sig_mises = workvec(1)
      call mm_return_values( "avg_princ_stress", elem, workvec, ngp )
      avg_sig_1 = workvec(1) ! maximum value
      if( ldebug ) then
        write(out,9010) avg_eps_plas, avg_sig_mises
      end if
      return
c
 9000 format(5x"... enter chk_get_ele_vals. elem: ",i8)
 9010 format(10x,'... leave chk_get_ele_vals with:',
     &   /,10x,'avg eps pls:',f12.8,' avg_sig_mises:',f8.3)
 9100 format('>> FATAL ERROR: routine ',
     & 'chk_get_ele_vals.',
     & 10x,'inconsistent condition. for element:',i8,
     & ' job terminated.'//)
c
      end subroutine chk_get_ele_vals
c
c
c     ****************************************************************
c     *                                                              *
c     *        internal subroutine chk_output_kill_messages_1        *
c     *                                                              *
c     *                   last modified : 9/19/21 rhd                *
c     *                                                              *
c     *                     Gurson crack growth                      *
c     *                                                              *
c     ****************************************************************
c
      subroutine chk_output_kill_messages_1
c
      implicit none
c
      integer :: device
      logical :: fexists, kill_by_Oddy, kill_by_plastic_limit
      double precision ::  total_plastic_strain, max_Oddy
      character(len=50) :: distortion_msg      
c
      if( ext_gurson ) then
         write(out,9002) elem, porosity, ext_spacing, ext_shape    
         return
      endif
c     
      if( distortion_killed ) then
        total_plastic_strain = eps_plas
        if( eps_plas < zero ) total_plastic_strain = zero ! roundoff
        max_Oddy = Oddy_metrics(elem_ptr,2) / Oddy_metrics(elem_ptr,1)
        kill_by_Oddy = max_Oddy >= Oddy_critical_ratio 
        kill_by_plastic_limit = total_plastic_strain >= 
     &                          distortion_plastic_limit
        distortion_msg = ' '
        if( kill_by_Oddy ) distortion_msg = '** Oddy exceeds limit'
        if( kill_by_plastic_limit ) distortion_msg =
     &                '** eps plastic exceeds limit'
        if( kill_by_Oddy .and. kill_by_plastic_limit )
     &     distortion_msg = '** Oddy & eps plastic exceed limits'
        write(out,9200) elem, porosity, max_Oddy, total_plastic_strain, 
     &                  distortion_msg(1:len(distortion_msg))
      else  !  killed by porosity
      if( .not. use_distortion_metric ) write(out,9000) elem, porosity
      if( use_distortion_metric )  write(out,9010) elem, porosity,
     &              Oddy_metrics(elem_ptr,2) / Oddy_metrics(elem_ptr,1)   
      end if 
c   
c              use same file name for list of deleted elements
c              as used for scs since both cannot be used in
c              same model 
      
      if( .not. smcs_deleted_list_file_flag ) return ! delete elements file
      inquire( file=smcs_deleted_list_file_name, exist=fexists )
      open( newunit=device,file=smcs_deleted_list_file_name,
     &      form='formatted', status='unknown', 
     &      access='sequential', position='append' )   
      if( .not. fexists ) then
        if( .not. use_distortion_metric ) write(device,9018)
        if( use_distortion_metric ) write(device,9019)
      end if   
      if( .not. use_distortion_metric )
     &  write(device,9020) current_load_time_step,
     &                   elem, porosity, sig_mean, sig_mises, eps_plas
      if( use_distortion_metric ) then
         max_Oddy = Oddy_metrics(elem_ptr,2) / Oddy_metrics(elem_ptr,1)
         write(device,9020) current_load_time_step,
     &                   elem, porosity, sig_mean, sig_mises, eps_plas,
     &                   max_Oddy
      end if         
      close(unit=device)    
                 

 9000 format(/,' @1  >> element death invoked for element: ',i7,                  
     & '.   f: ',f6.3)                                                          
 9002 format(/,'   >> element death invoked for element: ',i7,                  
     & '.   f, W, X: ',3f8.3)                                                   
 9010 format(/,' @1  >> element death invoked for element: ',i7,                  
     & '.   f: ',f6.3,2x,'Oddy metric: ',f8.3)                                                          
 9018 format('!',/,'!   Gurson quantities when elements deleted',/,
     & '!   (element deletion begins at listed step number)',
     & /,'!',/,
     & 6x,'step     elem   porosity     mean stress    mises stress',
     & '  plastic strain')
 9019 format('!',/,'!  Gurson quantities when elements deleted',/,
     & '!   (element deletion begins at listed step number)',
     & /,'!',/,
     & 6x,'step     elem   porosity     mean stress    mises stress',
     & '  plastic strain  Oddy metric')
 9020 format(2x,i8,i9,3x,f8.5,2x,e14.6,2x,e14.6,2x,e14.6,2x,f7.2)    
 9200 format(/,' @2  >> element removed by distortion metrics: ',i7,
     & '  f: ', f5.3,' max Oddy ratio: ',f7.2,2x,'plastic strain: ',
     & f4.2,2x,a)
c
      end subroutine chk_output_kill_messages_1
c
c
c     ****************************************************************
c     *                                                              *
c     *        internal subroutine chk_output_kill_messages_3        *
c     *                                                              *
c     *                   last modified : 8/29/21 rhd                *
c     *                                                              *
c     *                     SMCS crack growth                        *
c     *                                                              *
c     ****************************************************************
c
      subroutine chk_output_kill_messages_3
c
      implicit none
c
      integer :: device
      logical :: fexists
c
      if( smcs_type <= 4 ) then
        write(out,9010) elem, eps_plas, eps_crit, avg_triaxiality,
     &                avg_zeta, avg_bar_theta
      end if
      if( smcs_type == 5 ) then
        write(out,9025) elem, tear_param, smcs_type_5_tp_critical
      end if
      if( .not. smcs_deleted_list_file_flag ) return ! delete elements file
      inquire( file=smcs_deleted_list_file_name, exist=fexists )
      open( newunit=device,file=smcs_deleted_list_file_name,
     &      form='formatted', status='unknown', 
     &      access='sequential', position='append' )   
      if( .not. fexists ) write(device,9018)
      write(device,9020) current_load_time_step,
     &                   elem, eps_plas, eps_crit,
     &                   avg_triaxiality, avg_zeta, avg_bar_theta,
     &                   tear_param 
      close(unit=device)    
c
      return
c
 9010 format(/,' >> element death option invoked before next step',             
     &       /,'    element: ',i7,' is now killed.',                            
     &       /,'    plastic strain: ',f8.5,' is > limit of: ',f8.5,
     &       /,'    plastic strain weighted values: bar_T, ',
     &         'bar_xi, bar_theta:',f6.3,2x,f6.3,2x,f6.3)    
 9018 format('!',/,'!  SMCS quantities when critical plastic ',
     & 'strain reached for elements', /,
     & '!   (element deletion begins at listed step number)',
     & /,'!',/,
     & 3x,'step     elem   ebarp  ebarp_crit bar_T  bar_xi',
     &  '   bar_theta   tear_param')  
 9020 format(2x,2i8,f8.5,1x,f8.5,1x,f6.3,2x,f6.3,2x,f6.3,4x,f8.3)      
 9025 format(/,' >> element death option invoked before next step',             
     &       /,'    element: ',i7,' is now killed.',                            
     &       /,'    tearing parameter: ',f9.3,' is > limit of: ',
     &        f9.3 )
c 
      end subroutine chk_output_kill_messages_3
c
c
c     ****************************************************************
c     *                                                              *
c     *        internal subroutine chk_output_kill_messages_4        *
c     *                                                              *
c     *                   last modified : 8/29/21 rhd                *
c     *                                                              *
c     *                    cohesive element/material growth          *
c     *                                                              *
c     ****************************************************************
c
      subroutine chk_output_kill_messages_4
c
      implicit none
c
c                 start killing this element right now. record status               
c                 for subsequent packet output if neeeded. write usual              
c                 output about element death.                                   
c
      cohes_type         = iprops(27,elem)                               
      option_exponential = cohes_type .eq. 4                             
      option_ppr         = cohes_type .eq. 6                             
      option_cavit       = cohes_type .eq. 7                             
      found_exponential  = found_exponential .or.                         
     &                     option_exponential                             
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
     &   write(out,9200) elem, values(6)/values(7), values(8)            
      if( option_ppr )   write(out,9220) elem                            
      if( option_cavit ) write(out,9230) elem  
c
      return
c
 9200 format(/,'   >> element death invoked for element: ',i7,                  
     & '.   Deff/Dpeak: ',f5.2,' Teff/Tpeak: ',f5.2)                            
 9220 format(/,'   >> element death invoked for element: ',i7,                  
     & '.  PPR cohesive option')                                                
 9230 format(/,'   >> element death invoked for element: ',i7,                  
     & '.  cavit cohesive option')                                              
 9400 format(/,'FATAL ERROR: mixed cohesive options not allowed',               
     & /,      '             at present in WARP3D with crack growth.'           
     & /,      '             found while processing element: ',i7,              
     & /,      '             job aborted.' )  
c 
      end subroutine chk_output_kill_messages_4
c     ****************************************************************
c     *                                                              *
c     *        internal subroutine chk_output_kill_messages_5        *
c     *                                                              *
c     *                   last modified : 8/29/21 rhd                *
c     *                                                              *
c     *                 user directed element deletion               *
c     *                                                              *
c     ****************************************************************
c
      subroutine chk_output_kill_messages_5
c
      implicit none
c
      write(out,9240) elem                         
c
      return
c
 9240 format(/,'   >> element death invoked for element: ',i7,                  
     & '.  user-directed removal')                                              
c
      end subroutine chk_output_kill_messages_5
c     ****************************************************************
c     *                                                              *
c     *        internal subroutine chk_elem_kill_output_packets      *
c     *                                                              *
c     *                   last modified : 3/7/21 rhd                 *
c     *                                                              *
c     *     packet output for new element starting deletion          *
c     *                                                              *
c     ****************************************************************
c
      subroutine chk_elem_kill_output_packets
c      implicit none
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
c     ****************************************************************
c     *                                                              *
c     *        internal subroutine chk_elem_kill_in_order            *
c     *                                                              *
c     *                   last modified : 4/13/21 rhd                *
c     *                                                              *
c     *     after starting deletion of new element, see if other     *
c     *     elements must also be killed to satisfy the user input   *
c     *     command to kill in order                                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine chk_elem_kill_in_order
c
      use constants
c
      implicit none
c
      integer :: i, index_killed, elem, elem_ptr
      logical :: killed_found_in_order_list
      logical, parameter :: ldebug = .false.
      double precision :: eps_plas_at_death, mises_at_death_local,
     &                    sig_1_at_death_local, max_eps_plas
c
      if( ldebug ) write (out,*) '>>> chk_elem_kill_in_order....'       
c
c              have started deaths for 1 or more elements just now.
c              user has provided a kill-in-order list.
c              start at end of user supplied order list. work backwards.
c              find first killed element if any. kill all elements from 
c              that point to beginning of list unless already 
c              being killed.
c            
      killed_found_in_order_list = .false.
c
      do i = num_kill_order_list, 1, -1                                   
        elem     = kill_order_list(i)                                    
        elem_ptr = dam_ptr(elem)                                                
        if( dam_state(elem_ptr) == 0 ) cycle ! not yet killed
        killed_found_in_order_list = .true.
        index_killed = i
        exit
      end do ! i
c
      if( .not. killed_found_in_order_list ) return 
c
c              start killing elements from beginning of in-order list up to
c              first element in list already being killed.
c              
      do i = 1, index_killed-1
        elem     = kill_order_list(i)                                    
        elem_ptr = dam_ptr(elem)                                                
        if( dam_state(elem_ptr) .ne. 0 ) cycle ! already being killed
        write(out,9100) elem                                              
        num_elements_killed = num_elements_killed + 1  
        if( standard_kill ) then                    
          call chk_update_killed_energy( elem )                                  
          call chk_kill_element_now( elem, ldebug )                                   
          call chk_store_ifv( elem, elem_ptr, ldebug )                            
          call chk_update_node_elecnt( elem, ldebug )                             
          if ( release_type .eq. 2 )                                         
     &       call growth_set_dbar( elem, elem_ptr, ldebug, -1 ) 
        elseif( use_mesh_regularization ) then
          call chk_get_ele_vals( elem, eps_plas_at_death,
     &                           mises_at_death_local, 
     &                           sig_1_at_death_local, ! max prin value
     &                           max_eps_plas )
          smcs_eps_plas_at_death(elem_ptr) = eps_plas_at_death
          smcs_stress_at_death(elem_ptr) = sig_1_at_death_local
          dam_state(elem_ptr) = 1
          call chk_get_d( elem, zero, sig_1_at_death_local, 
     &                    smcs_d_values(elem_ptr), out )
        else
          write(out,9200) elem
          call die_abort
        end if                                                                  
      end do  ! i 
c   
      return
c
 9100 format(/,' >> element death option invoked before next step',             
     &       /,'    element: ',i7,' is now killed per user request',               
     &       /,'    via kill in order input command')    
 9200 format('>>>> FATAL ERROR: chk_elem_kill_in_order. elem',i8,
     &     /,10x,'Job terminated...',//)             
                                         
c
      end subroutine   chk_elem_kill_in_order 
c                                                             
c     ****************************************************************
c     *                                                              *
c     *        internal subroutine chk_elem_kill_make_packets        *
c     *                                                              *
c     *                   last modified : 3/7/21 rhd                 *
c     *                                                              *
c     *     new element starting deletion. output packets            *
c     *                                                              *
c     ****************************************************************
c
      subroutine chk_elem_kill_make_packets
c
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
             write(out,*) ' cavit death packets not implemented'             
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
      end subroutine chk_elem_kill                                                                      
c
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine chk_dam_init2                *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 12/6/20 rhd                *          
c     *                                                              *          
c     *        called when 1st element is deleted to setup.          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine chk_dam_init2( debug ) 
c                                            
      use global_data, only : out, nonode, noelem, iprops, c ! coords
c                                                                               
      use elem_extinct_data, only : dam_node_elecnt, dam_face_nodes             
      use main_data, only : incmap, incid, elems_to_blocks,                     
     &                      inverse_incidences                                  
      use elem_block_data, only : cdest_blocks                                  
      use damage_data, only : no_killed_elems, release_type, dam_ptr, 
     &                        crk_pln_normal_idx, crack_plane_coord,
     &                        gurson_cell_size                                                           
c                                                                               
      implicit none                                                    
c                                                                               
      logical, intent(in) :: debug
c
      integer :: i, elem, elem_ptr, incptr, num_enodes, k1, k2,
     &           face_count, blk, rel_elem, enode, snode
      logical, parameter :: local_debug = .false.                                            
      double precision :: node_coord(3), coord
      double precision, parameter :: plane_tol = 0.01d0                                 
      integer, dimension(:,:), pointer :: cdest                                 
c                                                                               
      if( debug ) write(out,*) '>>>> in chk_dam_init2'                             
c                                                                               
      no_killed_elems = .false.                                                 
      do i = 1, nonode                                                          
         dam_node_elecnt(i) = inverse_incidences(i)%element_count               
      end do                                                                    
      if( debug ) then                                                         
         do i = 1, nonode                                                       
            write(out,9030) i, dam_node_elecnt(i)                              
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
        if( debug ) write(out,*) '<<<< leaving chkdam_init2'                        
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
      if( debug ) write(out,*) '<<<< leaving chk_dam_init2'                        
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
c     *                   last modified : 12/6/20 rhd                *          
c     *                                                              *          
c     *        check each node in the structure to                   *          
c     *        see if they are no longer connected to any element.   *          
c     *        if so, it constrains them.                            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine chk_free_nodes( debug )                                        
      use global_data ! old common.main
c                                                                               
      use elem_extinct_data, only : dam_node_elecnt                             
      use main_data,         only : cnstrn, cnstrn_in                          
      use damage_data, only : csttail  
      use constants                                                    
c                                                                               
      implicit none                                                    
c              
      logical, intent(in) :: debug
c
      integer :: node, num_dof, j, glbal_dof
c
      if( debug ) write(out,*) '>>> in chk_free_nodes <<<'                     
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
       if( debug ) then                                                   
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
c     *                   last modified : 12/6/20 rhd                *          
c     *                                                              *          
c     *     print status of killable elements or released nodes at   *
c     *     the beginning of a load step.                            *       
c     *     optional smcs states output for smcs model
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine dam_print( step, iter ) 
c                                       
      use damage_data, only : crack_growth_type, print_status, 
     &                        print_top_list, const_front, smcs_states                                                           
c
      implicit none    
c
      integer, intent(in) :: step, iter  
c
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
c                                                                           
c     ****************************************************************          
c     *                                                              *          
c     *               service function chk_killed                    *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 7/7/21 rhd                 *          
c     *                                                              *          
c     *     given an element, returns true if the                    *          
c     *     element has been or is in-progress of being killed       *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      logical function chk_killed( elem ) result( t_chk_killed )                                      
c                                                                               
      use elem_extinct_data, only : dam_state                                   
      use damage_data, only : dam_ptr                                                          
c                                                                               
      implicit none
c
      integer, intent(in) :: elem
c
      integer :: elem_ptr                                         
c                                                                               
      t_chk_killed = .false.  
      if( .not. allocated( dam_ptr ) ) return                                                    
      elem_ptr =  dam_ptr(elem)    ! assume checking elsewhere                                            
      if( elem_ptr == 0 ) return ! element is not killable                                             
      t_chk_killed = dam_state( elem_ptr) > 0
c                                                                               
      return                                                                    
      end                                                                       
c
c     ****************************************************************
c     *                                                              *
c     *                    subroutine chk_killed_blk                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/28/21 rhd                *
c     *                                                              *
c     *     for the specified block of elements, return a logical    *
c     *     vector indicating if each element is killed or active.   *
c     *     "killed" can have two meanings based on type             *
c     *     also a flag if the whole block is killed                 *
c     *                                                              *
c     *     service routine for tanstf and do_nleps_block            *
c     *                                                              *
c     ****************************************************************
c
      subroutine chk_killed_blk( span, felem, step, block_no, 
     &                           status_vec, block_killed, type )
c
      use global_data, only : out, mxvl 
      use elem_extinct_data, only : dam_state, smcs_d_values
      use damage_data, only : growth_by_kill, dam_ptr, 
     &                        use_mesh_regularization,
     &                        tol_regular =>
     &                        tolerance_mesh_regularization     
      use constants
c
      implicit none
c
c              parameters
c
      integer, intent(in)  :: span, block_no, felem, type, step
      logical, intent(out) :: status_vec(mxvl), block_killed
c      
c              locals
c
      integer :: elem, elem_ptr, i
      logical :: standard_kill_method
      logical, parameter :: ldebug = .false.
      double precision :: chk_value
c
      block_killed = .false.
      status_vec   = .false. ! all entries

      if( .not. growth_by_kill ) return
c 
      if( step <= 2 ) then
        if( .not. allocated(dam_state) ) then         ! sanity check                                
          write(out,9100) 1                                                         
          call die_abort                                                          
        end if 
      end if
c
      if( dam_ptr(felem) == 0 ) return ! no killable elems in blk
c
c               standard_kill_method = no mesh regularization
c               for standard, element is killed immediately. state > 0
c
c               for mesh regularization, have to check current value
c               of the damage parameter 'd'. if > 0, element is 
c               being killed (and d<1.0) or has been 
c               fully killed (d>=1.0)
c
c               type = 1 used in tanstf. mark element killed if
c                      d > 0.0. tanstf will use either the [Ke] saved
c                      at death degraded by 1.0 - d or set [Ke] =0
c                      if d >= 1.0
c               type = 2 used in drive_eps_sig_internal_forces.
c                      mark element killed only if d > 1.0.
c                      for d < 1.0, the element strains, stresses
c                      and internal forces are computed as usual.
c                      routine addifv will degrade them by 1 - d before
c                      assembly into structure. if d >= 1.0 we don't
c                      want to process strains, stresses, internal 
c                      forces for the element.
c 
c
c               run sanity check. every element in blk must be 
c               killable. 
c
      if( step <= 2 ) then  
         do  i = 1, span  
           elem = felem + i - 1
           elem_ptr = dam_ptr( elem )
           if( elem_ptr == 0 ) then ! bad data structures
             write(out,9100) 3                                                         
             call die_abort                                                          
           end if 
         end do 
      end if 
c
      standard_kill_method = .not. use_mesh_regularization
c
      if( standard_kill_method ) then
        do i = 1, span  
          elem = felem + i -1
          elem_ptr = dam_ptr( elem )
          status_vec(i) = dam_state(elem_ptr) > 0
        end do
        if( all( status_vec(1:span) ) ) block_killed = .true.
        return
      end if
c
      if( .not. use_mesh_regularization ) then
          write(out,9100) 4
          call die_abort
      end if
c
      if( .not. allocated( smcs_d_values ) ) then ! really bad
          write(out,9100) 5
          call die_abort
      end if
c
      if( type < 1 .or. type > 2 ) then
          write(out,9100) 6
          call die_abort
      end if
c 
      chk_value = zero
      if( type == 2 ) chk_value = tol_regular
c
      do i = 1, span  
        elem = felem + i - 1
        elem_ptr = dam_ptr( elem )
        status_vec(i) = smcs_d_values(elem_ptr) > chk_value
      end do
c
      if( all( status_vec(1:span) ) ) block_killed = .true.
c
      return       
c      
 9100 format(/,'FATAL ERROR: chk_killed_blk. Contact WARP3D group',             
     &       /,'             Job terminated at ',i1,//)                         
c
      end
c     ****************************************************************          
c     *                                                              *          
c     *                     subroutine damage_update                 *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 12/18/2020 rhd             *          
c     *                                                              *          
c     *        certain damage models need variables updated outside  *
c     *        element stress/history data (e.g. smcs )              *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine damage_update( now_step, now_iter )
c
      use global_data, only : out, noelem
      use elem_extinct_data, only : dam_state
      use damage_data, only :  crack_growth_type, dam_ptr        
      use dam_param_code, only : dam_param_3_get_values 
c                                             
      implicit none
c
      integer :: now_step, now_iter
c                                       
      integer :: elem, elem_ptr, dowhat
      logical :: debug, local_debug, process, ld1                                                   
      double precision :: d1, d2, d3, d4, d5, d6, d7, d8, d9, d10
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
c$OMP PARALLEL DO PRIVATE( elem, elem_ptr, dowhat, d1, d2, d3, d4, 
c$OMP&                     d5, d6, d7, d8, d9, ld1, d10 )
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
     &                                 d5, d6, d7, d8, dowhat, d9,
     &                                 ld1, d10 )                
         end select 
c                         
      end do  ! over all model elements                                         
c$OMP END PARALLEL DO
c      
      return    
      end
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine chk_store_ifv                *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 03/9/21 rhd                *          
c     *                                                              *          
c     *        setup to gradually reduce newly killed element        *
c     *        internal forces -- for standard death process         *
c     *        use_mesh_regularization -- nothing to do.             *
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine chk_store_ifv( elem, elem_ptr, gdebug )  
c                           
      use global_data, only : out, iprops, load                                                                               
      use elem_extinct_data, only : dam_ifv, dam_state                          
      use elem_block_data,   only : edest_blocks                                
      use main_data,         only : elems_to_blocks                             
      use damage_data, only : release_type, use_mesh_regularization  
      use constants                                    
c                                                                               
      implicit none
c
      integer, intent(in) :: elem, elem_ptr
      logical, intent(in) :: gdebug
c
      integer :: num_edof, i, blk, rel_elem                                                     
      integer, dimension(:,:), pointer :: edest
      logical :: ldebug                                 
c     
      ldebug = .false.
c                                                                          
      if( ldebug ) write (*,*) '>>>>> Entering chk_store_ifv'      
c
      if( use_mesh_regularization ) return                 
c                                                                               
c            initialize state of release for element                            
c                                                                               
      dam_state(elem_ptr) = 1                                                   
c                                                                               
c            subtract the element internal forces from the                      
c            structure load vector. for all killable elements                   
c            which are not in the process of releasing, we                      
c            grab their internal forces at start of each load                   
c            step in another routine so they are available                      
c            here.                                                              
c                                                                               
      num_edof = iprops(2,elem) * iprops(4,elem)                                
c                                                                               
c            for the traction separation release, provide                       
c            an option to immediately                                           
c            drop the internal forces by a fraction, then release the           
c            remaining forces linearly with opening displacement.               
c            current reduction is zero.                                         
c                                                                               
      if( release_type .eq. 2 ) then                                           
        do i = 1, num_edof                                                      
            dam_ifv(i,elem_ptr) = one * dam_ifv(i,elem_ptr)                     
        end do                                                                  
      end if                                                                    
c                                                                               
      blk         = elems_to_blocks(elem,1)                                     
      rel_elem    = elems_to_blocks(elem,2)                                     
      edest       => edest_blocks(blk)%ptr                                      
      do i = 1, num_edof                                                        
        load(edest(i,rel_elem)) = load(edest(i,rel_elem)) -                     
     &                            dam_ifv(i,elem_ptr)                           
      end do                                                                    
c                                                                               
      if( ldebug ) then                                                         
         write (*,*) '>>>>> load in storeifv for elem:',elem                    
         do i = 1, num_edof                                                     
            write (*,'("new load:",i5,1x,e14.6)') edest(i,rel_elem),            
     &        load(edest(i,rel_elem))                                           
         end do                                                                 
         write (*,*) '>>>>>> Leaving chk_store_ifv'                                  
      end if                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                  subroutine chk_kill_element_now             *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 9/19/21 rhd                *          
c     *                                                              *          
c     *        This routine actually kills an element in the model.  *          
c     *        Material properties and history data for the element  *          
c     *        are modified so that subsequent element [k]s will be  *          
c     *        zero. the existing stresses are zeroed. the props and *          
c     *        history changes are such that stresses will remain    *          
c     *        zero for all subsequent steps.                        *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine chk_kill_element_now( ielem, debug ) 
c                                  
      use global_data, only : out, iprops, props, nstrs, nstr,
     &                             current_load_time_step
      use main_data, only : elems_to_blocks, cohesive_ele_types                                                              
      use elem_block_data, only : history_blocks, urcs_n_blocks,                
     &                            urcs_n1_blocks, history_blk_list,
     &                            eps_n_blocks, eps_n1_blocks 
      use damage_data, only : dam_ptr, crack_growth_type,
     &                        use_mesh_regularization,
     &                        use_distortion_metric,
     &                        smcs_removed_list_file_flag,
     &                        smcs_removed_list_file_name
      use elem_extinct_data, only : smcs_start_kill_step, Oddy_metrics
      use constants           
c
      implicit none   
c   
      integer, intent(in) :: ielem
      logical, intent(in) :: debug
c 
c              locals
c
      integer :: elem_type, material_type, ngp, blk, rel_elem, device,
     &           hist_size, hist_offset, gp, hisloc, start, end, i,
     &           se, ss, elem_ptr, npts
      double precision :: porosity
      double precision, dimension(:), pointer :: history, urcs_n1,
     &                                           urcs_n, eps_n, eps_n1
      logical :: gurson, regular_material, cohesive_elem, 
     &           output_fully_killed_msg, fexists
      character(len=80) :: file_name
c                                                                               
      if( debug ) write (out,*) '>>>> in chk_kill_element'                           
c                                                                               
c              for regular elements (hex, tet, et.c),                           
c              set young's modulus and poisson's ratio to zero for              
c              element so that any future request for a linear stiffness        
c              computes zero. for cohesive elements, we set a row of            
c              element props table = 1 so the element and cohesive model        
c              routines know element is killed.                                 
c                                                                               
      elem_type     = iprops(1,ielem)   
      material_type = iprops(25,ielem)      
      gurson        = material_type .eq. 3  ! also could just be mises                                      
      regular_material = .not. gurson
      cohesive_elem = cohesive_ele_types(elem_type)    
c                         
      if( .not. cohesive_elem ) then                                           
         props(7,ielem) = 0.0    !   E and nu                                                
         props(8,ielem) = 0.0                                                  
      else  ! cohesive elements                                                                    
         props(7:9,ielem) = 0.0
         props(20,ielem)  = 0.0
         props(21,ielem)  = 0.0
         iprops(32,ielem) = 1                                                   
      end if                                                                    
c                                                                              
c              for each gauss point of the element:                             
c               (a)  reset history data to a null state characteristic          
c                    of an unstressed element                                   
c               (b)  set stresses, streains at n and n+1 to zero.                        
c                                                                               
      ngp           = iprops(6,ielem)                                           
      blk           = elems_to_blocks(ielem,1)                                  
      rel_elem      = elems_to_blocks(ielem,2)                                  
      hist_size     = history_blk_list(blk)                                     
      hist_offset   = (rel_elem-1)*hist_size*ngp                                
      history       => history_blocks(blk)%ptr                                  
      urcs_n1       => urcs_n1_blocks(blk)%ptr                                  
      urcs_n        => urcs_n_blocks(blk)%ptr                                   
      urcs_n1       => urcs_n1_blocks(blk)%ptr                                  
      eps_n         => eps_n_blocks(blk)%ptr                                   
      eps_n1        => eps_n1_blocks(blk)%ptr                                   
c    
      if( gurson ) then ! don't zero porosity
         do gp = 1, ngp                                                            
           hisloc = hist_offset + (gp-1)*hist_size 
           porosity = history(hisloc+5)                          
           history(hisloc+1:hisloc+hist_size) = zero 
           history(hisloc+5) = porosity  
         end do                                         
      end if
c
      if( cohesive_elem ) then
         do gp = 1, ngp                                                            
            hisloc = hist_offset + (gp-1)*hist_size                                
            history(hisloc+1) = zero                                            
            history(hisloc+2) = zero                                            
            history(hisloc+3) = zero                                            
         end do                                                                    
      end if
c
      if( regular_material ) then
         do gp = 1, ngp                                                            
           hisloc = hist_offset + (gp-1)*hist_size 
           history(hisloc+1:hisloc+hist_size) = zero 
         end do
      end if
c                                                                               
      start = (rel_elem-1)*nstrs*ngp + 1                                    
      end   = start + ngp*nstrs - 1                                     
      do i = start, end                                                
        urcs_n1(i) = zero                                                      
        urcs_n(i)  = zero                                                      
      end do  
c                                                                  
      start = (rel_elem-1)*nstr*ngp + 1                                    
      end   = start + ngp*nstr - 1                                     
      do i = start, end                                                
        eps_n1(i) = zero                                                      
        eps_n(i)  = zero                                                      
      end do                                                                    
c 
      output_fully_killed_msg = .true.
      if( crack_growth_type .ne. 3 ) return
      if( .not. use_mesh_regularization ) return
      if( .not. output_fully_killed_msg ) return    
      elem_ptr = dam_ptr(ielem)                                               
      if( elem_ptr == 0 ) then
         write(out,9000) ielem
         call die_abort
      end if
c      
      if( .not. smcs_removed_list_file_flag ) return ! remove elements file
      inquire( file=smcs_removed_list_file_name, exist=fexists )
      open( newunit=device,file=smcs_removed_list_file_name,
     &      form='formatted', status='unknown', 
     &      access='sequential', position='append' )   
      if( .not. fexists ) write(device,9018)
      se = current_load_time_step
      ss = smcs_start_kill_step(elem_ptr)
      if( use_distortion_metric ) then
         write(device,9020) ielem, ss, se, se-ss, 
     &    Oddy_metrics(elem_ptr,2) / Oddy_metrics(elem_ptr,1)
      else
         write(device,9020) ielem, ss, se, se-ss
      end if   
      close(unit=device)    
c                                                                          
      if( debug ) write(out,*) '>>>> leaving chk_kill_element_now'                      
c                                                                               
      return
c
 9000 format('>> FATAL ERROR: routine ',
     & 'chk_kill_element_now',
     & 10x,'inconsistent condition. for element:',i8,
     & ' job terminated.'//)
 9018 format('!',/,
     & '!  SMCS steps numbers for element removal', /,
     & 6x,'elem  start step   end step # steps release',
     & '  max Oddy ratio')  
 9020 format(2x,i8,4x,i6,5x,i6,8x,i6,4x,f9.2)
 9100 format(10x,8f5.2)
c
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                subroutine chk_update_killed_energy           *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 03/9/21 rhd                *          
c     *                                                              *          
c     *        An element is being killed. We need to save the       *          
c     *        internal and plastic work right now so it can be      *          
c     *        included in future work totals. this is because all   *          
c     *        such data is about to be zeroed.                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine chk_update_killed_energy( ielem )  
c                                
      use global_data, only : out, iprops, nstrs, killed_ele_int_work,
     &                        killed_ele_pls_work
      use main_data, only : elems_to_blocks, cohesive_ele_types
      use elem_block_data, only : urcs_n_blocks, element_vol_blocks
      use constants
c            
      implicit none
c
      integer :: ielem
c
      integer :: elem_type, ngp, blk, rel_elem, sig_start, gp
      double precision :: sum_internal_energy, sum_plastic_energy,                       
     &                    element_volume, rgp                                                  
      double precision, pointer :: urcs_n(:) 
      logical :: cohesive_elem   
      logical, parameter :: debug = .false.                                           
c                                                                               
      if( debug ) write (out,*) '>>>> in update_killed_energy'                   
c                                                                               
c              for regular elements (hex, tet, et.c),                           
c              compute the average internal energy and plastic energy           
c              at the gauss points. these values are multipled by element       
c              volume and added to accumluated totals for killed                
c              elements.                                                        
c                                                                               
      elem_type     = iprops(1,ielem)                                           
      cohesive_elem = cohesive_ele_types(elem_type)                             
      if( cohesive_elem ) return                                               
      ngp           = iprops(6,ielem)                                           
      blk           = elems_to_blocks(ielem,1)                                  
      rel_elem      = elems_to_blocks(ielem,2)                                  
      urcs_n        => urcs_n_blocks(blk)%ptr                                   
      sig_start     = (rel_elem-1) * nstrs * ngp                                
c                                                                               
      sum_internal_energy = zero                                                
      sum_plastic_energy  = zero                                                
c                                                                               
      do gp = 1, ngp                                                            
        sum_internal_energy = sum_internal_energy + urcs_n(sig_start+7)         
        sum_plastic_energy  = sum_plastic_energy  + urcs_n(sig_start+8)         
        sig_start           = sig_start + nstrs                                 
      end do                                                                    
c                                                                               
      rgp                 = dble(ngp)                                           
      sum_internal_energy = sum_internal_energy / rgp                           
      sum_plastic_energy  = sum_plastic_energy  / rgp                           
c                                                                               
      element_volume      = element_vol_blocks(blk)%ptr(rel_elem)               
      killed_ele_int_work = killed_ele_int_work +                               
     &                      element_volume * sum_internal_energy                
      killed_ele_pls_work = killed_ele_pls_work +                               
     &                      element_volume * sum_plastic_energy                 
c                                                                               
      if ( debug ) write (out,*) '>>>> leaving update_killed_energy'              
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                   subroutine chk_update_node_elecnt          *          
c     *                                                              *          
c     *                       written by : ag                        *          
c     *                                                              *          
c     *                   last modified : 03/9/21 rhd                *          
c     *                                                              *          
c     *        Updates the vector that holds the number of           *          
c     *        elements connected to each node after a node is       *          
c     *        killed.                                               *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine chk_update_node_elecnt( ielem, debug ) 
c                            
      use global_data, only : out, iprops
      use elem_extinct_data, only : dam_node_elecnt                             
      use main_data, only : incid, incmap, inverse_incidences                   
c                                                                               
      implicit none
c
      integer, intent(in) :: ielem
      logical, intent(in) :: debug                                                             
c 
      integer :: nnode, i, node, elem, ndof                                                                              
c
      if( debug ) write (out,*) '>>>> in update_node_elecnt '                    
c                                                                               
c                loop over nodes connected to killed element.  Get node         
c                number.  Then decrease the count of the 
c                elements connected to that node by one.                                           
c                
      nnode = iprops(2,ielem)                                                   
      do i = 1, nnode                                                           
         node = incid(incmap(ielem)+(i-1))                                      
         elem = inverse_incidences(node)%element_list(1)                        
         ndof = iprops(4,elem)                                                  
         dam_node_elecnt(node)= dam_node_elecnt(node) - 1                       
      end do                                                                    
      if( debug ) write (out,*) '>>>> leaving update_node_elecnt'                  
c                                                                               
      return                                                                    
 9000 format ('>> node ',i2,' of elem ',i7,' is ',i6)                           
 9010 format ('>>   elems connected to node is ',i2)                            
 9020 format('old constraints for node ',i7,' :',3(1x,e14.6)) 
c                  
      end                                                                       
c
c     ****************************************************************
c     *                                                              *
c     *                    subroutine chk_get_d                      *
c     *                                                              *
c     *                   last modified : 4/30/21 rhd                *
c     *                                                              *
c     *     get current value of damage parameter "d" for mesh       *
c     *     regularization                                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine chk_get_d( elem, deps_plas, stress_at_death, d_now,
     &                      out )
      use damage_data, only : tol_regular =>
     &                        tolerance_mesh_regularization,
     &                        regular_npoints, regular_length,
     &                        regular_up_max, regular_points, ! 10x2
     &                        regular_type, regular_alpha,
     &                        regular_Gf, regular_m_power
      use constants
c
      implicit none
c
      integer, intent(in) :: out, elem
      double precision, intent(in) :: deps_plas, stress_at_death
      double precision, intent(out) :: d_now
      double precision, external :: chklint
c
      double precision :: u_plastic, u_normed, Gf, m, up, uf, T0,
     &                    T1, T2, Gtilde, Gratio, alpha
      logical :: local_debug
c
      local_debug = elem == -1
      u_plastic = regular_length * deps_plas
c
      select case( regular_type )
      case( 1 )
        u_normed =  u_plastic / regular_up_max
        d_now = chklint( u_normed, regular_npoints, 
     &                   regular_points(1,1), regular_points(1,2) )
        if( u_plastic >= regular_up_max ) d_now = one 
      case( 2 )
        if( stress_at_death < zero ) then
          write(out,9300) elem, stress_at_death
          call die_abort
        end if 
        Gf = regular_Gf
        m  = regular_m_power
        T0 = stress_at_death
        up = u_plastic
        uf = Gf * (m+one) / m / T0
        T1 = one - (up/uf)**m / (m + one)
        Gtilde = T0 * up * T1
        d_now = Gtilde / Gf
        if( up >= uf ) d_now = one ! happens w/ power law + large jumps
        if( local_debug ) write(out,9100) elem, u_plastic, regular_Gf,
     &                    m, uf, Gtilde, stress_at_death, d_now
      case( 3 )
        if( stress_at_death < zero ) then
          write(out,9300) elem, stress_at_death
          call die_abort
        end if 
        Gf     = regular_Gf
        m      = regular_m_power
        alpha  = regular_alpha
        T0     = stress_at_death
        up     = u_plastic
        uf     = Gf * (m+one) / m / T0
        T1     = one - (up/uf)**m / (m + one)
        Gtilde = T0 * up * T1
        Gratio = Gtilde / Gf
        T2     = one - exp( -alpha * Gratio )
        d_now  = T2 / ( one - exp(-alpha))
        if( up >= uf ) d_now = one ! happens w/ power law + large jumps
        if( local_debug ) write(out,9200) elem, u_plastic, regular_Gf,
     &                    m, alpha, uf, Gtilde, Gratio, stress_at_death,
     &                    d_now
      case default
         write(out,9000)
         call die_abort
      end select
c
      if( d_now < zero ) d_now = zero
      if( d_now > tol_regular ) d_now = one 
      if( local_debug ) write(out,9020) d_now
c
      return
 9000 format('>> FATAL ERROR: routine ',
     & 'chk_get_d.',
     & 10x,'inconsistent condition.',
     & ' job terminated.'//)
 9020 format(10x,'leaving chk_get_d w/ d_now: ',f6.3)
 9100 format(/,10x,'... compute d for element: ',i8,
     & /,15x,'u plastic:         ',f10.6,
     & /,15x,'Gf:                ',f10.6,
     & /,15x,'m-power:           ',f10.6,
     & /,15x,'computed uf:       ',f10.6,
     & /,15x,'computed G-tilde:  ',f10.5,
     & /,15x,'stress @ death:    ',f10.3,
     & /,15x,'d now:             ',f10.3 ) 
 9200 format(/,10x,'... compute d for element: ',i8,
     & /,15x,'u plastic:         ',f10.6,
     & /,15x,'Gf:                ',f10.6,
     & /,15x,'m-power:           ',f10.6,
     & /,15x,'alpha:             ',f10.6,
     & /,15x,'computed uf:       ',f10.6,
     & /,15x,'computed G-tilde:  ',f10.5,
     & /,15x,'Gratio:            ',f10.5,
     & /,15x,'stress @ death:    ',f10.3,
     & /,15x,'d now:             ',f10.3 ) 
 9300 format('>> FATAL ERROR: routine ',
     & 'chk_get_d. the stress @ element death < 0. this is an',
     & 10x,'inconsistent condition.',
     & ' job terminated.'//)
c
      end
c     ****************************************************************
c     *                                                              *
c     *         double precision function chklint                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 4/1/21 rhd                  *
c     *                                                              *
c     *    execute linear interpolation on a tabular function        *
c     *    of a single variable where the x values are sorted in     *
c     *    increasing value but are not req'd to be uniformly spaced *
c     *                                                              *
c     ****************************************************************
c
      function chklint( xvalue, n, x, y ) result( value )
      implicit none
c
      integer :: n
      double precision :: xvalue, x(n), y(n), value
c
      integer :: i, point
      double precision :: x1, x2, y1, y2
c
      if( xvalue .le. x(1) ) then
        value = y(1)
        return
      end if
c
      if( xvalue .ge. x(n) ) then
        value = y(n)
        return
      end if
c
      do point = 2, n
        if( xvalue .gt. x(point) ) cycle
        x1 = x(point-1)
        x2 = x(point)
        y1 = y(point-1)
        y2 = y(point)
        value = y1 + (xvalue-x1)*(y2-y1)/(x2-x1)
        return
      end do
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *          logical function chkcrack_is_elem_deleted           *
c     *                                                              *
c     *                   last modified : 3/28/21 rhd                *
c     *                                                              *
c     *     has element been (fully) deleted from the solution?      *
c     *                                                              *
c     ****************************************************************
c
c
      logical function chkcrack_is_elem_deleted( elem ) result( test )
c
      use global_data, only : out
      use dam_param_code, only : dam_param
      use elem_extinct_data, only : dam_state, smcs_d_values
      use damage_data, only : dam_ptr,
     &                        use_mesh_regularization,
     &                        tol_regular =>
     &                        tolerance_mesh_regularization     
c
      implicit none
c
      integer, intent(in) :: elem
c
      integer :: elem_ptr
      logical :: std_kill  
c
      test = .false.
      elem_ptr = dam_ptr(elem)   ! assumes much checking before here
c                                            
      if( elem_ptr == 0 ) return ! element not killable        
      if( dam_state(elem_ptr) == 0 ) return ! not yet started killing
c
      std_kill = .not. use_mesh_regularization
      if( std_kill ) then
        test = .true.
        return ! dam_state /= 0 element actually deleted
      end if
c
      if( use_mesh_regularization ) then
        test = smcs_d_values(elem_ptr) > tol_regular ! fully deleted 
        return
      end if
c
      write(out,9000) elem
      call die_abort
c
 9000 format('>> FATAL ERROR: routine ',
     & 'chkcrack_is_elem_deleted.',
     & 10x,'inconsistent condition. for element:',i8,
     & ' job terminated.'//)
c
      end function chkcrack_is_elem_deleted
c
c     ****************************************************************
c     *                                                              *
c     *          subroutine chkcrack_Oddy_needed                     *
c     *                                                              *
c     *            last modified : 8/29/21 rhd                       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine chkcrack_Oddy_needed( felem, lflag )
c
      use global_data, only : out
      use elem_extinct_data, only : dam_state
      use damage_data, only : dam_ptr, crack_growth_type, 
     &                        use_distortion_metric
      use constants
c
      implicit none
c
      integer, intent(in) :: felem
      logical, intent(out) :: lflag
c
      integer :: elem, elem_ptr
      logical :: process
c
      lflag = .false.
c
      if( .not. use_distortion_metric ) return
      process = crack_growth_type .eq. 1  .or. crack_growth_type .eq. 3
      if( .not. process ) return  ! no element deletion growth
      elem = felem
      elem_ptr = dam_ptr(elem)                                               
      if( elem_ptr == 0 ) return ! no elements in block killable
      lflag = .true. ! killable element in block
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *          subroutine chkcrack_Oddy_initial_print              *
c     *                                                              *
c     *            last modified : 9/10/21 rhd                       *
c     *                                                              *
c     *     after step = 1, iter = 1, print a file of Oddy metrics   *
c     *     if requested by the user. file name is fixed             *
c     *                                                              *
c     ****************************************************************
c
      subroutine chkcrack_Oddy_initial_print
c
      use global_data, only : out, noelem
      use damage_data, only :  dam_ptr, Oddy_print_initial
      use elem_extinct_data, only : oddy_metrics_initial
      use constants
c
      integer :: elem, elem_ptr, fileno, element_min, element_max,
     &           element_max_ratio
      logical, save :: already_printed = .false.
      real :: Oddy_min, Oddy_max, global_min, global_max, 
     &        global_max_ratio
      character(len=80) :: fname
c
      if( .not. Oddy_print_initial ) return
      if( already_printed ) return ! only print file once
      already_printed = .true.
c
      if( .not. allocated( dam_ptr ) ) then
         write(out,9900) 1
         call die_abort
      end if
      if( .not. allocated( Oddy_metrics_initial ) ) then
         write(out,9900) 2
         call die_abort
      end if
c          
      fname = "Oddy_initial_metrics_killable_elements"
      open(newunit=fileno,file=fname,status="unknown")
      write(fileno,9010) 
c
      global_min = 1000000.0
      global_max = -100000.0
      global_max_ratio = -100000.0
c
      do elem = 1, noelem
       elem_ptr = dam_ptr(elem)                                               
       if( elem_ptr == 0 ) cycle 
       Oddy_min = Oddy_metrics_initial(elem_ptr,1)
       Oddy_max = Oddy_metrics_initial(elem_ptr,2)
       if( Oddy_min < global_min ) then
         global_min = Oddy_min
         element_min = elem
       end if
       if( Oddy_max > global_max ) then
         global_max = Oddy_max
         element_max = elem
       end if
       ele_ratio = Oddy_max / oddy_min
       if( ele_ratio > global_max_ratio ) then
         global_max_ratio = ele_ratio
         element_max_ratio = elem
       end if
       write(fileno,9000) elem, Oddy_min, Oddy_max, ele_ratio
      end do
c
      write(fileno,9020) global_min, element_min, global_max,
     &                   element_max, global_max_ratio, 
     &                   element_max_ratio
      write(fileno,*) " "
c
      close(unit=fileno,status='keep')
      deallocate( Oddy_metrics_initial ) ! not saved in restart
c
      return
c
 9000 format(2x,i8,2x,f12.3,2x,f12.3,2x,f12.3)
 9010 format(/," Oddy metrics for killable elements @ t=0",
     & /,/,'    element      Oddy min      Oddy max      max/min')
 9020 format(/,5x,'minimum:       ',f12.3,' @ element: ',i8,/,
     &         5x,'maximum:       ',f12.3,' @ element: ',i8,/,
     &         5x,'maximum ratio: ',f12.3,' @ element: ',i8,/)
 9900 format('>> FATAL ERROR: @',i2,' in chkcrack_Oddy_initial_print',
     & /,    '                Job terminated',//)           
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *          subroutine chkcrack_Oddy                            *
c     *                                                              *
c     *            last modified : 9/10/21 rhd                       *
c     *                                                              *
c     *     compute/store the Oddy distortion paramter for block     *
c     *     of killable elements. only called if calcs needed        *
c     *                                                              *
c     ****************************************************************
c
      subroutine chkcrack_Oddy( step, itype, span, felem, gpn, detJ,
     &                          jac )
c
      use global_data, only : out, mxvl, iprops, props
      use dam_param_code, only : dam_param
      use elem_extinct_data, only : dam_state, Oddy_metrics, 
     &    oddy_metrics_initial
      use damage_data, only : dam_ptr, Oddy_print_initial
      use constants
c
c              compute/store the current Oddy paramter for a block 
c              of killable elements at gpn.
c              allocate space during step 1 as needed. 
c              called by gtmat while processing blocks in parallel
c              to compute various required deformation gradients F.
c              gtmat1 comoutes coordinate jacobian (and det) for
c              current element (deformed) shape.
c              this routines runs on each thread.
c
      implicit none
c
      integer, intent(in) :: step, itype, span, felem, gpn
      double precision, intent(in) :: detJ(mxvl), jac(mxvl,3,3)
c
      integer :: i, elem, elem_ptr
      logical :: process
      logical, parameter :: ldebug1 = .false., ldebug2 = .false.,
     &                      ldebug3 = .false., print = .false.
      logical, save :: header = .true.
      double precision :: Oddy_value, Oddy_prior_value, Oddy_new
      real :: ymod, old_value, new_value, old_min, new_min, old_max,
     &        new_max
c
      if( ldebug1 ) write(out,9000) step, itype, felem, gpn, span
      if( .not. allocated( dam_ptr ) ) then
         write(out,8900)
         call die_abort
      end if
c
c              use special loop for step 1. undeformed Jac passed
c
      if( itype == 1 ) then
        do i = 1, span
          elem = felem + i - 1
          elem_ptr = dam_ptr(elem)                                               
          if( elem_ptr == 0 ) cycle ! really fatal error
          call chkcrack_Oddy_compute( i, Oddy_value )
          old_value = Oddy_metrics(elem_ptr,1)
          new_value = min( sngl(Oddy_value), old_value )
          if( new_value <= 0.0) new_value = 1.0 ! cubical elems
          Oddy_metrics(elem_ptr,1) = new_value
          Oddy_metrics(elem_ptr,2) = 0.0
          if( ldebug1 .and. elem == 1 ) then
              write(out,9010) elem, detJ(i), Oddy_value, new_value
              write(out,9012) jac(i,1:3,1:3)
          end if
          if( Oddy_print_initial ) then
            old_min = Oddy_metrics_initial(elem_ptr,1)
            new_min = min( old_min, sngl(Oddy_value) )
            old_max = Oddy_metrics_initial(elem_ptr,2)
            new_max = max( old_max, sngl(Oddy_value) )
            Oddy_metrics_initial(elem_ptr,1) = new_min
            Oddy_metrics_initial(elem_ptr,2) = new_max
          end if
          if( .not. print ) cycle
          if( header ) write(out,9310)
          header = .false.
          write(out,9320) elem, gpn, new_min, new_max
        end do
        return
      end if
c
c              regular loop after step 1. all data allocated. put 
c              current Oddy value in row 2 of table. only need to 
c              process elements being killed and not yet fully 
c              deleted from model
c
      if( itype /= 2 ) then
        write(out,9300) 
        call die_abort
      end if
c        
      do i = 1, span
        elem = felem + i - 1
        ymod = props(7,elem)  ! has been zeroed for fully killed element
        if( ymod <= 0.0 ) cycle  ! props is single precision
        elem_ptr = dam_ptr(elem) 
        call chkcrack_Oddy_compute( i, Oddy_value )
        Oddy_prior_value = dble( Oddy_metrics(elem_ptr,2) )
        if( Oddy_value <= zero ) Oddy_value = one
        Oddy_new = max( Oddy_prior_value, Oddy_value )
        Oddy_metrics(elem_ptr,2) = sngl(Oddy_new) 
        if( step <4 .and. ldebug3 .and. elem == 1 ) then
            write(out,9200) elem, gpn, Oddy_prior_value, Oddy_value, 
     &                      Oddy_new
            write(out,9012) jac(i,1:3,1:3)
        end if    
      end do
c
      return 
c
 8900 format(/1x,'>>>>> FATAL ERROR: routine chkcrack_Oddy.'
     &  /,'                    should not be called',
     &  /,'        Execution terminated.',/)
 9000 format(/,'.... Enter chkcrack_Oddy:',/,
     &  10x,'step, itype, felem, gpn, span: ',i7,i3,i7,i3,i4) 
 9010 format(10x,'elem, detJ, Oddy value, Oddy new:',i8,3f10.4)
 9012 format(15x,3f10.5)
 9100 format(10x, i8, 3f10.4 )
 9200 format('.... Oddy trace. elem, gp, Oddy prior, Oddy now, ',
     &   ' Oddy new: ', i8,i4,3f10.4) 
 9300 format(/,">>>>>  FATAL ERROR: routine chkcrack_Oddy",/,
     &         "         Job terminated.")
 9310 format(5x,'--------   Oddy metric values undeformed element  ',
     &   '--------',/,
     &       15x,'Element     Gauss Pt     min Oddy value   ',
     &       'max Oddy value') 
 9320 format(15x, i7, 8x, i2, 9x, f8.3, 9x, f8.3)          
c
      contains
c     --------
c
      subroutine chkcrack_Oddy_compute( relem, Oddy_value )
      implicit none
      integer, intent(in) :: relem
      double precision, intent(out) :: Oddy_value
c
      double precision :: tens(3,3), trans_tens(3,3), C(3,3), 
     &                    Cvec(9), t1, t2
      equivalence ( C, Cvec )
c
      tens(1:3,1:3) = jac(relem,1:3,1:3) / detJ(relem)**third 
      trans_tens = transpose( tens )
      C = matmul( trans_tens, tens )
      t2 = ( C(1,1) + C(2,2) + C(3,3) )**2 / three
      t1 = dot_product( Cvec, Cvec ) ! tensor contraction C:C
      Oddy_value = one +  (t1 - t2)
      if( ldebug2 .and. elem == 1 ) then
         write(out,9000)
         write(out,9010) C
         write(out,9015) t1, t2
      end if
c
      return
 9000 format(10x,"... [C] ....:")
 9010 format(15x,3f10.5)
 9015 format(10x,"... t1,t2: ", 2f10.5)
c
      end subroutine chkcrack_Oddy_compute
      end subroutine chkcrack_Oddy
        



