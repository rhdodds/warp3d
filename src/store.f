c     ****************************************************************
c     *                                                              *
c     *                      subroutine store                        *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 6/28/2018 rhd              *
c     *                                                              *
c     *                  writes analysis restart file                *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine store( savnam, savfil, sbflg1, sbflg2 )
      use global_data ! old common.main
c
      use elem_block_data, only: nonlocal_flags, nonlocal_data_n1,
     &                           initial_state_data
      use elem_extinct_data
      use node_release_data
      use elem_load_data
      use main_data
      use segmental_curves
      use mod_mpc, only : mpcs_exist, num_user_mpc, user_mpc_table,
     &                    tied_con_mpcs_constructed, num_tied_con_mpc,
     &                    tied_con_mpc_table
      use stiffness_data, only : total_lagrange_forces
      use contact
      use damage_data
      use hypre_parameters
      use performance_data
      use distributed_stiffness_data, only: parallel_assembly_allowed,
     &      parallel_assembly_used, distributed_stiffness_used,
     &      initial_map_type, final_map_type
      use mm10_defs, only :
     & one_crystal_hist_size, common_hist_size
      use erflgs
c
      implicit none

c
c               parameters
c
      logical :: sbflg1, sbflg2
      character :: savnam*8
      character(len=*) :: savfil
c
c               locals
c
      integer :: dum, fileno, last, prec_fact, node, local_length,
     &           i,  how_defined, node_count, num_patterns, mpc,
     &           ntrm, itab, num_rows, num_cols, ilist, ulen, nsize
      integer, parameter :: check_data_key=2147483647
      character :: dbname*100
      logical :: nameok, scanms, delfil, wrt_nod_lod, write_table,
     &           initial_stresses_exist
      real :: dumr, rsum
      real, parameter :: rzero = 0.0
      double precision :: dumd
      double precision, parameter :: dzero=0.0d0
      character(len=1) :: dums
c
c                       at least one load step must have been computed
c                       before a restart file is created.  if no load
c                       steps have yet been computed, write and error
c                       message and return.
c
      if ( ltmstp .eq. 0 ) then
         call errmsg (293, dum, dums, dumr, dumd )
         goto 9999
      endif
c
c                                 figure out filename:
c                                       if savfil unset (no filename),
c                                         make a filename from the structure nam
c                                       if savfil set, resolve name.
c                                       if neither are set, error.
c
      fileno = 11
      dbname = ' '
      if ( scanms(savfil,'namenone',8) ) then
        if ( scanms(savnam,'itsblank',8) ) then
          call errmsg( 151, dum, dums, dumr, dumd )
          dbname = 'default_db'
        else
          last = 8
          call name_strip( savnam, last )
          dbname = savnam(1:last) // '_db'
        end if
      else
        call tilde( savfil, dbname, nameok )
        if ( .not.nameok ) dbname = 'default_db'
      end if
c
c                                 check file.  if file there, delete it
c
      inquire( file = dbname, exist = delfil )
      if ( delfil ) then
         open( fileno, file=dbname, status='old', access='sequential',
     &     form='unformatted')
         close( fileno, status= 'delete' )
       end if
c
c                       open the file.
c
      open( fileno, file=dbname, status='new', access='sequential',
     &     form='unformatted', recordtype='segmented' )
c
c                       MPI:
c                         we need to gather all the stresses, strains,
c                         element volumes, initial state arrays  back
c                         to the root processor to store in the restart
c                         database.
c
      call wmpi_get_str ( 2 )
      call wmpi_get_str ( 4 )
      call wmpi_get_str ( 5 )
      if( initial_state_option ) call wmpi_get_initial_state
c
c                       after each block of data written on the file,
c                       we write a record containing a 'check' variable
c                       to catch any inconsistencies between the save
c                       here and reopen
c
c                       write error flags
c
      write(fileno) numnod,numel,fatal,coor,elprop,elinc,constr,block
c
c
c                       write out integer (scalar) variables
c
c
      write(fileno) mathed, noelem, nonode, nodof, csthed, inctop,
     &              crdtop, nummat, numcol, histep,
     &              lowstp, ltmstp, nlibel, numlod, mxiter,
     &              mniter, lodhed, lgnmcn, mxlitr,
     &              nogp, nprs, nplrs, nelblk, numgrp, lgoump,
     &              max_current_pts, max_current_curves,
     &              num_seg_curve_sets,
     &              solver_flag, old_solver_flag, solver_memory,
     &              solver_threads, eq_node_force_len,
     &              precond_type,
     &              hsolver_type, precond_printlevel,
     &              solver_printlevel, max_step_limit,
     &              hypre_max, levels, symme, error_count,
     &              precond_fail_count, ntimes_assembly,
     &              initial_map_type, final_map_type,
     &              coarsening, agg_levels, interpolation, relaxation,
     &              sweeps, cf, cycle_type, max_levels,
     &              one_crystal_hist_size, common_hist_size,
     &              initial_state_step
      write (fileno) check_data_key
c
c
c                       write out logical and character (scalar)
c                       variables.
c
c
      initial_stresses_exist = allocated( initial_stresses )
      write(fileno) prnres,halt,
     &              linmas,ifvcmp,zrocon,growth_k_flag,
     &              newtrn,lsldnm,incflg,
     &              prlres,adaptive_flag,batch_messages,
     &              signal_flag, scalar_blocking, stname,
     &              qbar_flag, solver_out_of_core, solver_scr_dir,
     &              show_details, temperatures, sparse_stiff_output,
     &              sparse_stiff_binary, sparse_research,
     &              solver_mkl_iterative, output_packets,
     &              temperatures_ref, fgm_node_values_defined,
     &              fgm_node_values_used,
     &              hyp_trigger_step, hyp_first_solve,
     &              time_assembly, parallel_assembly_allowed,
     &              parallel_assembly_used,
     &              distributed_stiffness_used, nonlocal_analysis,
     &              umat_serial, umat_used, asymmetric_assembly,
     &              extrapolate, extrap_off_next_step,
     &              divergence_check, diverge_check_strict,
     &              line_search, ls_details, initial_stresses_exist,
     &              initial_stresses_user_routine,
     &              initial_state_option, initial_stresses_input
      write(fileno) sparse_stiff_file_name, packet_file_name,
     &              initial_stresses_file
      write (fileno) check_data_key
c
c
c                       write out double precision (scalar) variables.
c
c
      write(fileno) dt,nbeta,emax,fmax,prdmlt,total_mass,ext_work,
     &              beta_fact,total_model_time,scaling_adapt,eps_bbar,
     &              overshoot_limit, control_load_fact,
     &              old_load_fact, min_load_fact, killed_ele_int_work,
     &              killed_ele_pls_work,hypre_tol,threshold,filter,
     &              loadbal, start_assembly_step,
     &              assembly_total, truncation, relax_wt,
     &              relax_outer_wt, mg_threshold, ls_min_step_length,
     &              ls_max_step_length, ls_rho, ls_slack_tol
      write (fileno) check_data_key
c
c
c                       write out real  (scalar) variables.
c
      write (fileno) time_limit
      write(out,9000)
c
c                       write out integer arrays.
c
      call wrt2d( fileno, outmap,mxlbel,nlibel,mxelmp )
      call wrtbk( fileno, matlst, mxmat )
      call wrtbk( fileno, invdst, nodof )
      call wrtbk( fileno, plrlst, mxlsz )
      call wrtbk( fileno, stprng, mxlc*2 )
      call wrtbk( fileno, gpmap,  nogp )
      call wrtbk( fileno, incmap, noelem )
      call wrtbk( fileno, crdmap, nonode )
      call wrtbk( fileno, dstmap, nonode )
      call wrtbk( fileno, cstmap, nodof )
      call wrtbk( fileno, state,  nogp )
      call wrtbk( fileno, incid,  inctop )
      call wrtbk( fileno, prslst, mxlsz )
      call wrtbk( fileno, lodlst, mxlc )
      write (fileno) check_data_key
      call store_cmplx_int( fileno, 1 )
      call store_cmplx_int( fileno, 2 )
      call store_cmplx_int( fileno, 3 )
      call wrt2d( fileno, elblks(0,1) , 4, 4, nelblk  )
      call wrtbk( fileno, cp, mxedof )
      call wrtbk( fileno, dcp, mxedof )
      call wrt2d( fileno, icp, mxutsz, mxutsz, 2  )
      call wrtbk( fileno, num_seg_points, max_seg_curves )
      call wrtbk( fileno, seg_curves_type, max_seg_curves )
      call wrtbk( fileno, seg_curve_table, (max_seg_curves+1)*
     &            max_seg_curve_sets )
      call wrt2d( fileno, imatprp, mxmtpr, mxmtpr, mxmat )
      write (fileno) check_data_key
      write(out,9010)
c
c
c                       write out logical and character arrays.
c                       use blocked store for large logical vectors.
c                       write character variables/arrays directly.
c                       see some sizes set in main, param_def
c
c
      write(fileno) convrg, trace                  ! short logical vecs
      call wrtbk( fileno, trn, nonode)             ! logical vec
      call wrtbk( fileno, stpchk, max_step_limit ) ! logical vec
      call wrtbk( fileno, repeat_incid, noelem )   ! logical vec
c
      write(fileno) lodnam, lodtyp, matnam, elelib ! short char vecs
      write(fileno) smatprp                        ! char array
      write(fileno) check_data_key
      write(out,9020)
c
c                       write out real arrays.
c
      call wrt2d( fileno, props, mxelpr, mxelpr, noelem )
      call wrt2d( fileno, matprp, mxmtpr, mxmtpr, mxmat )
      call wrtbk( fileno, user_cnstrn_stp_factors, max_step_limit )
      call wrtbk( fileno, actual_cnstrn_stp_factors, max_step_limit )
      write (fileno) check_data_key
      if ( fgm_node_values_defined )
     &  call wrt2d( fileno, fgm_node_values, nonode, nonode,
     &              fgm_node_values_cols )
      write (fileno) check_data_key
      write(out,9030)
c
c
c                       write out double precision data
c                         1) global vectors
c
      prec_fact = 2
      call wrtbk( fileno, u, prec_fact*nodof )
      call wrtbk( fileno, c, prec_fact*nodof )
      call wrtbk( fileno, ifv, prec_fact*nodof )
      call wrtbk( fileno, load, prec_fact*nodof )
      call wrtbk( fileno, dload, prec_fact*nodof )
      call wrtbk( fileno, rload, prec_fact*nodof )
      call wrtbk( fileno, rload_nm1, prec_fact*nodof )
      if( allocated( total_user_nodal_forces ) ) then
         write(fileno) .true.
         call wrtbk( fileno, total_user_nodal_forces, prec_fact*nodof )
      else
         write(fileno) .false.
      end if
      call wrt2d( fileno, load_pattern_factors, mxlc*prec_fact,
     &            numlod*prec_fact, 2 )
      call wrtbk( fileno, cnstrn_in, prec_fact*nodof )
      write (fileno) check_data_key
      call wrtbk( fileno, v, prec_fact*nodof )
      call wrtbk( fileno, a, prec_fact*nodof )
      call wrtbk( fileno, du, prec_fact*nodof )
      call wrtbk( fileno, temper_nodes, prec_fact*nonode )
      call wrtbk( fileno, temper_nodes_ref, prec_fact*nonode )
      call wrtbk( fileno, temper_elems, prec_fact*noelem )
      call wrt2d( fileno, dmatprp, 2*mxmtpr, 2*mxmtpr, mxmat )
      write (fileno) check_data_key
      write(out,9040)
c
c                         2) blocked data structures for
c                             material point histories and [D]
c                             stresses
c                             material point rotations
c                             trial elastic stresses
c                             strains
c                             element volumes
c
      call store_blocks( fileno, 1 )
      write(out,9050)
      call store_blocks( fileno, 5 )
      write (fileno) check_data_key
      write(out,9055)
      call store_blocks( fileno, 2 )
      write (fileno) check_data_key
      write(out,9060)
      call store_blocks( fileno, 4 )
      write (fileno) check_data_key
      write(out,9080)
      call store_blocks( fileno, 6 )
      write (fileno) check_data_key
      write(out,9090)
c
c                         3) coordinate transformation matrices for
c                            node coordinate systems for constraints
c
      do node = 1, nonode
         if (trn(node)) then
            write (fileno) trnmat(node)%mat
         end if
      end do
      write (fileno) check_data_key
c
      call wrtbk( fileno, tol, prec_fact*mxcvtests)
      call wrtbk( fileno, mdiag, prec_fact*nodof )
      call wrt3d( fileno, seg_curves, prec_fact*max_seg_points,2,
     &            prec_fact*max_current_pts, 2, max_current_curves )
      call wrtbk( fileno, seg_curves_min_stress,
     &            prec_fact*max_current_curves )
      call wrtbk( fileno, seg_curves_value,
     &            prec_fact*max_current_curves )
      call wrtbk( fileno, seg_curves_ym,
     &            prec_fact*max_current_curves )
      call wrtbk( fileno, seg_curves_nu,
     &            prec_fact*max_current_curves )
      call wrtbk( fileno, seg_curves_alpha,
     &            prec_fact*max_current_curves )
      call wrtbk( fileno, seg_curves_gp_sigma_0,
     &            prec_fact*max_current_curves )
      call wrtbk( fileno, seg_curves_gp_h_u,
     &            prec_fact*max_current_curves )
      call wrtbk( fileno, seg_curves_gp_beta_u,
     &            prec_fact*max_current_curves )
      call wrtbk( fileno, seg_curves_gp_delta_u,
     &            prec_fact*max_current_curves )
      write (fileno) check_data_key
      write(out,9100)
c
c                       save the nodal loading data. there are
c                       two data structures and some scalars for the
c                       applied nodal loads and nodal temperaturres.
c
      wrt_nod_lod = .false.
      do i = 1, mxlc
       if ( node_load_defs(i)%node_count .gt. 0 ) wrt_nod_lod = .true.
       if ( node_load_defs(i)%how_defined .gt. 0 ) wrt_nod_lod = .true.
      end do
      write(fileno) wrt_nod_lod
      if ( wrt_nod_lod ) then
        do i = 1, mxlc
         node_count = node_load_defs(i)%node_count
         how_defined = node_load_defs(i)%how_defined
         write(fileno) node_count, how_defined
         write(fileno) node_load_defs(i)%user_file_name
         if ( node_count .gt. 0 ) then
          call wrtbk( fileno, node_load_defs(i)%nodal_loads(1,1),
     &                2*node_count )
         end if
        end do
      end if
      write (fileno) check_data_key
c
      write(fileno) num_loddat_blks, next_loddat_col
      do i = 1, num_loddat_blks
        call wrtbk( fileno, loddat_blocks(i)%block(1,1),
     &              mxndldcm*sizeof_loddat_blks )
      end do
      write (fileno) check_data_key
c
c                       save the element loading data, if given.
c
      do i = 1, mxlc
         write(fileno) elem_loads(i)%size
         if ( elem_loads(i)%size .gt. 0 ) then
            call wrt2d( fileno, elem_loads(i)%data(1,1),
     &           elem_loads(i)%size, elem_loads(i)%size, 3 )
            call wrtbk( fileno, elem_loads(i)%vals(1),
     &           elem_loads(i)%size )
            call wrtbk( fileno, elem_loads(i)%piston_tabnum(1),
     &           elem_loads(i)%size )
         endif
      end do
      write (fileno) check_data_key
      write(out,9110)
c
c                       save the element equivalent nodal force vectors
c                       for applied body forces, surface tractions,
c                       pressures.
c
      if ( eq_node_force_len .gt. 0 ) then
         call wrtbk( fileno, eq_node_force_indexes, noelem )
         call wrtbk( fileno, eq_node_forces,
     &               eq_node_force_len*prec_fact )
      end if
      write (fileno) check_data_key
c
c                       save the nonlinear load step data
c
      if ( .not. allocated( step_load_data ) ) call mem_allocate( 19 )
      do i = 1, max_step_limit
          write(fileno) step_load_data(i)%num_load_patterns
      end do
      do i = 1, max_step_limit
       num_patterns = step_load_data(i)%num_load_patterns
       if ( num_patterns .gt. 0 ) then
        write(fileno) step_load_data(i)%load_patt_num(1:num_patterns)
        write(fileno) step_load_data(i)%load_patt_factor(1:num_patterns)
       end if
      end do
      write (fileno) check_data_key
c
c
c                       now save the crack growth parameters, if
c                       needed. 1) integers, 2) doubles. then
c                       data arrays based on user selected
c                       growth options.
c
c
      write(fileno) crack_growth_type, max_dam_state, num_kill_elem,
     &              csttail, print_status, kill_order, no_killed_elems,
     &              num_print_list, num_kill_order_list,
     &              release_type, crk_pln_normal_idx,
     &              num_crack_plane_nodes, crack_front_start,
     &              crack_front_end, crkfrnt_garbage_start,
     &              crkfrnt_garbage_end, growth_by_kill,
     &              growth_by_release, min_steps_for_release,
     &              load_size_control_crk_grth,
     &              overshoot_control_crk_grth, overshoot_allocated,
     &              no_released_nodes, g_stp_cntrl_allocated,
     &              const_front, num_crack_fronts, num_nodes_thick,
     &              num_nodes_back, master_lines_set,
     &              num_nodes_grwinc, num_steps_min, load_reduced,
     &              all_elems_killed, num_elements_killed,
     &              enforce_node_release, num_ctoa_released_nodes
      write (fileno) check_data_key
c
      write(fileno) porosity_limit, gurson_cell_size,
     &              crack_plane_coord, release_fraction,
     &              critical_angle, char_length, release_height,
     &              crack_plane_sign, init_crit_ang, smcs_alpha,
     &              smcs_beta, CTOA_range, perm_load_fact,
     &              max_porosity_change, max_plast_strain_change,
     &              init_ctoa_dist, ctoa_dist, crkpln_srch_tol,
     &              max_deff_change, critical_cohes_deff_fract,
     &              ppr_kill_displ_fraction
      write (fileno) check_data_key
      call wrtbk( fileno, dam_ptr, noelem )
      write (fileno) check_data_key
      write(out,9120)
c
c                             save element extinction variables
c
      if ( growth_by_kill ) then
c
         call wrt2d( fileno, dam_ifv, prec_fact*mxedof,
     &               prec_fact*mxedof, num_kill_elem )
         call wrtbk( fileno, dam_state, num_kill_elem )
         call wrtbk( fileno, dam_blk_killed, nelblk )
c
         if ( print_status ) then
            call wrtbk( fileno, dam_print_list, num_print_list )
            call wrtbk( fileno, old_mises, num_print_list*prec_fact )
            call wrtbk( fileno, old_mean, num_print_list*prec_fact )
         end if
         write (fileno) check_data_key
c
         if ( kill_order ) then
            call wrtbk( fileno, kill_order_list, num_kill_order_list )
         end if
         write (fileno) check_data_key
c
         if ( .not. no_killed_elems ) then
            call wrtbk( fileno, dam_node_elecnt, nonode )
         end if
         write (fileno) check_data_key
c
         if ( release_type .eq. 2  ) then
            if ( .not. no_killed_elems ) then
               call wrtbk( fileno, dam_face_nodes, 4*num_kill_elem )
               call wrtbk( fileno, dam_dbar_elems,
     &              (2*prec_fact)*num_kill_elem )
            end if
         end if
         write (fileno) check_data_key
         write(out,9130)
c
c                                   if load step size reduction is on,
c                                   save porosity or plastic strain from
c                                   last step.
c
         if ( load_size_control_crk_grth ) then
            if ( crack_growth_type .eq. 1 ) then
               call wrtbk( fileno, old_porosity,
     &              prec_fact * num_kill_elem )
               call wrtbk( fileno, del_poros,
     &              prec_fact * mxstp_store)
            else if ( crack_growth_type .eq. 3) then
               call wrtbk( fileno, old_plast_strain,
     &              prec_fact * num_kill_elem )
            else if ( crack_growth_type .eq. 4 ) then
               call wrtbk( fileno, old_deff,
     &              prec_fact * num_kill_elem )
               call wrtbk( fileno, del_deff,
     &              prec_fact * mxstp_store )
            end if
         end if
         write (fileno) check_data_key
c
c                             save node release variables
c
      else if ( growth_by_release ) then
c
         call wrtbk( fileno, crack_plane_nodes, num_crack_plane_nodes )
         call wrtbk( fileno, inv_crkpln_nodes, nonode )
         call wrtbk( fileno, num_neighbors, num_crack_plane_nodes )
         call wrt2d( fileno, neighbor_nodes, mxconn, mxconn,
     &        num_crack_plane_nodes )
         call wrtbk( fileno, crkpln_nodes_state,
     &        num_crack_plane_nodes )
         call wrtbk( fileno, crkpln_nodes_react,
     &        num_crack_plane_nodes * prec_fact )
         call wrt2d( fileno, crack_front_nodes, num_crack_plane_nodes,
     &     num_crack_plane_nodes, 2 )
         write (fileno) check_data_key

         if ( release_type .eq. 2 ) call wrtbk ( fileno,
     &        node_release_frac,num_crack_plane_nodes * prec_fact )
         write (fileno) check_data_key

      if ( overshoot_control_crk_grth ) call wrt2d( fileno,
     &        old_angles_at_front, num_crack_plane_nodes * prec_fact,
     &        num_crack_plane_nodes * prec_fact, mxconn )
         write (fileno) check_data_key

         if ( const_front ) then
           call wrtbk( fileno, master_nodes, num_crack_fronts )
           call wrt2d( fileno, crack_front_list, num_crack_fronts*
     &                 num_nodes_grwinc,num_crack_fronts*
     &                 num_nodes_grwinc,num_nodes_thick )
           call wrt2d( fileno, master_lines, num_crack_fronts,
     &                 num_crack_fronts, num_nodes_back + 1 )
         end if
         write (fileno) check_data_key
         write(out,9140)
c
      end if
      write (fileno) check_data_key
c
c                       save data structures for contact
c
      write (fileno) use_contact
c
      if (use_contact) then
         write (fileno) num_contact
         call wrtbk( fileno, contact_shape, num_contact)
         call wrtbk( fileno, contact_outside, num_contact)
         call wrt2d( fileno, contact_cause, maxcontact, maxcontact,
     &        nonode)
         call wrtbk( fileno, contact_stiff, num_contact* prec_fact)
         call wrtbk( fileno, contact_depth, num_contact*prec_fact)
         call wrtbk( fileno, contact_fric, num_contact*prec_fact)
         write (fileno) check_data_key
         call wrt2d( fileno, cshape_param, maxcntprm* prec_fact,
     &        maxcntprm*prec_fact, num_contact)
         call wrt2d( fileno, cshape_rate, 3*prec_fact, 3*prec_fact,
     &        num_contact)
         call wrt2d( fileno, cshape_pnt, 3*prec_fact, 3*prec_fact,
     &        num_contact)
         call wrt2d( fileno, cshape_norm, 3*prec_fact, 3*prec_fact,
     &        num_contact)
         call wrt3d( fileno, cplane_vec, 3*prec_fact, 2*prec_fact,
     &        3*prec_fact, 2*prec_fact,  num_contact)
         write(out,9150)
      endif
      write (fileno) check_data_key
c
c                       save data structures for mpcs. current values
c                       of lagrange nodal forces
c
      write (fileno) mpcs_exist
      if (mpcs_exist) then
         write (fileno) num_user_mpc
         do mpc = 1, num_user_mpc
            ntrm = user_mpc_table(mpc)%num_terms
            write (fileno) ntrm
            write (fileno) user_mpc_table(mpc)%constant
            call wrtbk( fileno, user_mpc_table(mpc)%node_list, ntrm )
            call wrtbk( fileno, user_mpc_table(mpc)%dof_list, ntrm )
            call wrtbk( fileno, user_mpc_table(mpc)%multiplier_list,
     &                  ntrm )
            write (fileno) check_data_key
         end do
         write(out,9160)
      end if
c
      write (fileno) tied_con_mpcs_constructed
      if (tied_con_mpcs_constructed) then
         write (fileno) num_tied_con_mpc
         do mpc = 1, num_tied_con_mpc
            ntrm = tied_con_mpc_table(mpc)%num_terms
            write (fileno) ntrm
            write (fileno) tied_con_mpc_table(mpc)%constant
            call wrtbk( fileno, tied_con_mpc_table(mpc)%node_list, ntrm)
            call wrtbk( fileno, tied_con_mpc_table(mpc)%dof_list, ntrm )
            call wrtbk( fileno, tied_con_mpc_table(mpc)%multiplier_list,
     &                  ntrm )
            write (fileno) check_data_key
         end do
         write(out,9170)
      end if
c
      if( mpcs_exist .or. tied_con_mpcs_constructed ) then
        call wrtbk( fileno, total_lagrange_forces, nodof* prec_fact )
      end if
      write (fileno) check_data_key
c
c                       save the table definitions, save both single
c                       and double precision values.
c
      do itab = 1, max_tables
         write(fileno) tables(itab)%table_name
         write(fileno) tables(itab)%table_type
         num_rows = tables(itab)%num_rows
         num_cols = tables(itab)%num_cols
         write(fileno) num_rows
         write(fileno) num_cols
         if ( num_rows .ne. 0 ) then
            call wrt2d( fileno, tables(itab)%table_values_sgl(1,1),
     &           num_rows, num_rows, num_cols )
            call wrt2d( fileno, tables(itab)%table_values_dbl(1,1),
     &           num_rows*prec_fact, num_rows*prec_fact, num_cols )
         endif
      end do
      write(out,9180)
      write(fileno) check_data_key
c
c                       save the user list definitions, if given.
c
      do ilist = 1, max_user_lists
         write(fileno) user_lists(ilist)%name
         ulen =  user_lists(ilist)%length_list
         write(fileno) ulen
         if ( ulen .ne. 0 ) then
           call wrtbk(fileno,user_lists(ilist)%list(1),ulen)
         endif
      end do
      write(out,9190)
      write(fileno) check_data_key
c
c                       write the crystal data structures
c
      call write_cry_data( fileno )
      write(out,9200)
      write(fileno) check_data_key
c
c
c                       save nonlocal material state data if it exists.
c                       shared between solids & interface-cohesive
c
      if( nonlocal_analysis ) then
        write(fileno) check_data_key
        nsize = nonlocal_shared_state_size
        do i = 1, noelem
           if( .not. nonlocal_flags(i) ) cycle
           write(fileno) nonlocal_data_n1(i)%state_values(1:nsize)
        end do
        write(out,9195)
      end if
c
c                       save convergence history data for user
c                       routine to adjust next step loading,
c                       solution parameters, etc.
c
      write(fileno) check_data_key
      do i = 1, 5
        write(fileno) convergence_history(i)%step_converged,
     &         convergence_history(i)%adaptive_used,
     &         convergence_history(i)%iterations_for_convergence,
     &         convergence_history(i)%adapt_substeps
      end do
      write(fileno) run_user_solution_routine
      write(fileno) check_data_key
c
c                       save the release_cons_table if we have nodes
c                       undergoing release operations
c
      write_table = .false.
      if( allocated( release_cons_table ) ) write_table = .true.
      write(fileno) write_table
      if( write_table ) then
        do i = 1, nonode
          write(fileno) release_cons_table(1:3,i)%num_release_steps
          write(fileno) release_cons_table(1:3,i)%
     &                  remaining_steps_for_release
          write(fileno) release_cons_table(1:3,i)%reaction_force
        end do
      end if
      write(out,9210)
      write(fileno) check_data_key
c
c                       save for output commands file ... information
c                       save file name and bitmap list if it exists
c
      write(fileno) output_command_file
      local_length = 0
      if( allocated( output_step_bitmap_list ) )
     &   local_length = size( output_step_bitmap_list )
      write(fileno) local_length
      if( local_length > 0 )
     &   write(fileno) output_step_bitmap_list
      write(fileno) check_data_key
c
c                       save initial stress data
c
      if( initial_stresses_exist ) then
        call wrtbk( fileno, initial_stresses, prec_fact*noelem*6 )
        write(out,9230)
      end if
      write(fileno) check_data_key
c
c                       save initial state arrays if they exist yet
c
      if( initial_state_option ) then
          call store_blocks( fileno, 3 )
          write(out,9240)
      end if
      write(fileno) check_data_key
c
c                       USER routines data
c
      call uexternaldb_store( fileno, out )
      write(out,9220)
c
c                       to help indicate if a restart file has been corrupted.
c                       The number is 2**31 - 1

c                       save a final check variable for checking in restart
c                       to help indicate if a restart file has been corrupted.
c                       The number is 2**31 - 1
c
      write (fileno) check_data_key
c
c                       close the file
c
      close( fileno, status='keep' )
c
c          uexternaldb for Abaqus compatible support
c
      douextdb = 5  ! common.main. tells uexter.. what to do
      call wmpi_do_uexternaldb
c
      call errmsg( 193, dum, dbname, dumr, dumd )
c
c
 9999 sbflg1= .false.
      sbflg2= .false.
c
c
 1000 format ( 3x, e16.6 )
 9000 format(15x,'> scalars written...')
 9010 format(15x,'> integer arrays written...')
 9020 format(15x,'> logical & character arrays written...')
 9030 format(15x,'> real arrays written...')
 9040 format(15x,'> double precision arrays written...')
 9050 format(15x,'> histories and [Dt] blocks written...')
 9055 format(15x,'> stress blocks written...')
 9060 format(15x,'> material rotation blocks written...')
 9080 format(15x,'> strain blocks written...')
 9090 format(15x,'> element volume blocks written...')
 9100 format(15x,'> material curves written...')
 9110 format(15x,'> node and element loading data written...')
 9120 format(15x,'> crack growth parameters/status written...')
 9130 format(15x,'> element extinction data written...')
 9140 format(15x,'> crack growth by node release data written...')
 9150 format(15x,'> contact data written...')
 9160 format(15x,'> user-defined mpcs written...')
 9170 format(15x,'> tied mesh mpcs written...')
 9180 format(15x,'> table definitions written...')
 9190 format(15x,'> user list definitions written...')
 9195 format(15x,'> nonlocal material data...')
 9200 format(15x,'> crystal data written...')
 9210 format(15x,'> convergence history written...')
 9220 format(15x,'> user routine data written...')
 9230 format(15x,'> initial stress data written...')
 9240 format(15x,'> initial state data arrays written...')
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine store_blocks                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                   last modified : 6/28/2018 rhd              *
c     *                                                              *
c     *     write data into restart file for the requested data      *
c     *     structures which following the element blocking          *
c     *                                                              *
c     ****************************************************************
c
      subroutine store_blocks( fileno, proc_type )
      use global_data ! old common.main
c
      use elem_block_data, only : history_blocks, rot_n1_blocks,
     &                            eps_n_blocks, initial_state_data,
     &                            urcs_n_blocks, history_blk_list,
     &                            rot_blk_list,
     &                            eps_blk_list, urcs_blk_list,
     &                            element_vol_blocks,
     &                            cep_blocks, cep_blk_list
c
      implicit none
c
      integer :: fileno, proc_type
c
      integer :: blk

      integer, parameter :: check_data_key = 2147483647
c
c
      select case( proc_type )
c
c            element history data and [Dt]s
c
      case( 1 )
      do blk = 1, nelblk
       if ( history_blk_list(blk) .gt. 0 ) then
          write(fileno) history_blocks(blk)%ptr
          write(fileno) check_data_key
       end if
       if ( cep_blk_list(blk) .gt. 0 ) then
          write(fileno) cep_blocks(blk)%vector
          write(fileno) check_data_key
       end if
      end do
c
c            element rotation matrices for gauss (material) points
c            for geometric nonlinear element blocks.
c
      case( 2 )
      do blk = 1, nelblk
       if ( rot_blk_list(blk) .eq. 1 ) then
           write(fileno) rot_n1_blocks(blk)%ptr
           write(fileno) check_data_key
       end if
      end do
c
c            initial state arrays
c
      case( 3 )
      if( .not. allocated( initial_state_data ) ) then
        write(out,9020)
        call die_abort
      end if
c
      associate( x => initial_state_data )
      do blk = 1, nelblk
        if( allocated( x(blk)%W_plastic_nis_block ) ) then
          write(fileno) x(blk)%W_plastic_nis_block
          write(fileno) check_data_key
        else
          write(out,9010) blk
          call die_abort
        end if
      end do
      end associate
c
c            element strains at gauss points - all elements
c
      case( 4 )
      do blk = 1, nelblk
        if ( eps_blk_list(blk) .eq. 1 ) then
          write(fileno) eps_n_blocks(blk)%ptr
          write(fileno) check_data_key
        end if
      end do
c
c            element stresses at gauss points - all elements
c
      case( 5 )
      do blk = 1, nelblk
        if ( urcs_blk_list(blk) .eq. 1 ) then
           write(fileno) urcs_n_blocks(blk)%ptr
           write(fileno) check_data_key
        end if
      end do
c
c            element volumes
c
      case( 6 )
      do blk = 1, nelblk
         write(fileno) element_vol_blocks(blk)%ptr
         write(fileno) check_data_key
      end do
c
      case default
         write(out,9000)
        call die_abort
      end select
c
      return
c
 9000 format(/1x,'>>>>> error: routine store_blocks inconsistency.',
     &   /1x,    '             job terminated...')
 9010 format(/1x,'>>>>> error: routine store_blocks. action 3, blk:',
     &  i10, /1x,    '             job terminated...')
 9020 format(/1x,'>>>>> error: routine store_blocks. action 3',
     &  i10, /1x,    '             job terminated...')
c
      end

c
c     ****************************************************************
c     *                                                              *
c     *                 subroutine store_cmplx_int                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 06/15/01                   *
c     *                                                              *
c     *     write data into restart file for the requested data      *
c     *     structures which additional complexity                   *
c     *                                                              *
c     ****************************************************************
c
      subroutine store_cmplx_int( fileno, proc_type )
      use global_data ! old common.main
c
      use main_data, only : inverse_incidences, inverse_dof_map
c
      implicit integer (a-z)
      data check_data_key / 2147483647 /
c
      go to ( 100, 200, 300 ) proc_type
c
 100  continue
      write(fileno) ( inverse_incidences(stnd)%element_count, stnd = 1,
     &                nonode )
      return
c
 200  continue
      do stnd = 1, nonode
        ecount = inverse_incidences(stnd)%element_count
        write(fileno) ( inverse_incidences(stnd)%element_list(i),
     &                  i = 1, ecount )
      end do
      return
c
 300  continue
      do stnd = 1, nonode
        ecount = inverse_incidences(stnd)%element_count
        write(fileno) (( inverse_dof_map(stnd)%edof_table(i,j),
     &                  i = 1, ecount ), j = 1, 3 )
      end do
      return

      end
c
c     ****************************************************************
c     *                                                              *
c     *               subroutine write_cry_data                      *
c     *                                                              *
c     *                    written by : mcm                          *
c     *                                                              *
c     *                last modified : 12/23/2015 rhd                *
c     *                                                              *
c     *           Write the CP crystal definitions to file           *
c     *                                                              *
c     ****************************************************************
c
      subroutine write_cry_data(fileno)
      use crystal_data
      implicit none
      integer, intent(in) :: fileno
c
      integer :: nelem, mxcry
      intrinsic size
c
      write(fileno) defined_crystal
c
      if( defined_crystal ) then
        write(fileno) cry_multiplier
        write(fileno) c_array
        nelem = size(angle_input,1)
        mxcry = size(angle_input,2)
        write(fileno) nelem
        write(fileno) mxcry
        write(fileno) data_offset
        write(fileno) angle_input
        write(fileno) crystal_input
      end if
c
      write(fileno) srequired
c
      if( srequired)  then
        write(fileno) nangles
        write(fileno) simple_angles
        write(fileno) mc_array
      end if
c
      return
      end subroutine








