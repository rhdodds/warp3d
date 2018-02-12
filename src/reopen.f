c     ****************************************************************
c     *                                                              *
c     *                      subroutine reopen                       *
c     *                                                              *
c     *                      written by : bh                         *
c     *                                                              *
c     *                   last modified : 2/10/2018 rhd              *
c     *                                                              *
c     *          read restart file. get solution start up            *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine reopen( resnam, resfil, sbflg1, sbflg2 )
      use global_data ! old common.main
      use main_data
      use elem_block_data, only : nonlocal_flags, nonlocal_data_n,
     &                            nonlocal_data_n1
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
      character :: resnam*8
      character(len=*) :: resfil
c
c               locals
c
      integer :: dum, last, prec_fact, node, i, node_count, eload_size,
     &           num_patterns, err, dumi, mpc, ntrm, fileno,
     &           how_defined, itab, num_rows, num_cols, ilist,
     &           ulen, nsize, local_length
      character ::  dbname*100, string*80, dums*1
      logical :: flexst, scanms, nameok, initial_stresses_exist,
     &          read_nod_load, msg_flag, read_table, exist_flag
      real :: dumr, restart_time
      real, external :: wcputime
      double precision :: dumd
      double precision, parameter :: dzero=0.0d0
c
      msg_flag = .true.
c
c
c
c                       figure out filename:
c                               if resfil unset (no filename),
c                               make a filename from the structure name.
c                               if resfil set, resolve name.
c                               if neither are set, error.
c
      fileno = 11
      if ( scanms( resfil, 'namenone', 8 ) ) then
         if ( scanms( resnam, 'itsblank', 8 ) ) then
            call errmsg( 148, dum, 'no name', dumr, dumd )
            go to 9999
         else
            last = 8
            call name_strip( resnam, last )
            dbname = resnam(1:last) // '_db'
         end if
      else
         call tilde( resfil, dbname, nameok )
         if (.not.nameok) then
            call errmsg( 192, dum, dums, dumr, dumd )
            go to 9999
         end if
      end if
c
c                      find out if file exists.  if not, error.
c                               if so, open it.
c
      inquire( file=dbname, exist=flexst )
      if( .not. flexst ) then
         call errmsg( 148, dum, dbname, dumr, dumd )
         go to 9999
      end if
c
      open( fileno, file=dbname, status='old', access='sequential',
     &     form='unformatted', recordtype='segmented' )
c
c                       rewind the data base to insure positioning
c                       before the first record.
c
      rewind fileno
c
c                       after each block of data read from the file,
c                       we read a record that should contain a 'check'
c                       variable to catch any inconsistencies
c                       between the previous save and the reopen here
c
c
c                       read error flags
c
c
      read(fileno) numnod,numel,fatal,coor,elprop,elinc,constr,block
c
c
c                       read in integer (scalar) variables
c
c
      read(fileno) mathed, noelem, nonode, nodof, csthed, inctop,
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
     &              one_crystal_hist_size, common_hist_size
      call chk_data_key( fileno, 1, 0 )
      call mem_allocate( 4 ) ! vectors based on # nodes
c
c
c                       read in logical and character (scalar)
c                       variables.
c
c
      read(fileno) prnres,halt,
     &             linmas,ifvcmp,zrocon,growth_k_flag,
     &             newtrn,lsldnm,incflg,
     &             prlres,adaptive_flag,batch_messages,
     &             signal_flag, scalar_blocking, stname,
     &             qbar_flag, solver_out_of_core, solver_scr_dir,
     &             show_details, temperatures, sparse_stiff_output,
     &             sparse_stiff_binary, sparse_research,
     &             solver_mkl_iterative, output_packets,
     &             temperatures_ref, fgm_node_values_defined,
     &             hyp_trigger_step, hyp_first_solve,
     &             time_assembly, parallel_assembly_allowed,
     &             parallel_assembly_used,
     &             distributed_stiffness_used, nonlocal_analysis,
     &             umat_serial, umat_used, asymmetric_assembly,
     &             extrapolate, extrap_off_next_step,
     &             divergence_check, diverge_check_strict,
     &             line_search, ls_details, initial_stresses_exist,
     &             initial_stresses_user_routine
      read(fileno) sparse_stiff_file_name, packet_file_name,
     &             initial_stresses_file
      call chk_data_key( fileno, 1, 1 )
c
c
c                       read in double precision variables.
c
c
      read(fileno) dt,nbeta,emax,fmax,prdmlt,total_mass,ext_work,
     &             beta_fact,total_model_time,scaling_adapt,eps_bbar,
     &             overshoot_limit, control_load_fact,
     &             old_load_fact, min_load_fact, killed_ele_int_work,
     &             killed_ele_pls_work,hypre_tol,threshold,filter,
     &             loadbal, start_assembly_step,
     &             assembly_total, truncation, relax_wt,
     &             relax_outer_wt, mg_threshold, ls_min_step_length,
     &             ls_max_step_length, ls_rho, ls_slack_tol
      call chk_data_key( fileno, 1, 2 )
c
c
c                       read in real variables -- also initialize the
c                          time limit procedures
c
      read (fileno) time_limit
      call steptime(dum,1)
      call steptime(dum,4)
      write(out,9000)
c
c                       read in integer arrays.
c
      call mem_allocate( 3 )
      call rd2d( fileno, outmap,mxlbel,nlibel,mxelmp )
      call rdbk( fileno, matlst, mxmat )
      call init_maps( fileno, 2 )
      call rdbk( fileno, plrlst, mxlsz )
      call rdbk( fileno, stprng, mxlc*2 )
      call rdbk( fileno, gpmap,  nogp )
      call mem_allocate( 9 )
      call rdbk( fileno, incmap, noelem )
      call mem_allocate( 14 )
      call rdbk( fileno, crdmap, nonode )
      call rdbk( fileno, dstmap, nonode )
      call rdbk( fileno, cstmap, nodof )
      call rdbk( fileno, state, nogp )
      call rdbk( fileno, incid, inctop )
      call rdbk( fileno, prslst, mxlsz )
      call rdbk( fileno, lodlst, mxlc )
      call chk_data_key( fileno, 2, 0 )
      call init_maps( fileno, 4 )
      call init_maps( fileno, 5 )
      call init_maps( fileno, 6 )
      call init_maps( fileno, 7 )
      call rd2d( fileno, elblks(0,1) , 4, 4, nelblk )
      call rdbk( fileno, cp, mxedof )
      call rdbk( fileno, dcp, mxedof )
      call rd2d( fileno, icp, mxutsz, mxutsz, 2 )
      call rdbk( fileno, num_seg_points, max_seg_curves )
      call rdbk( fileno, seg_curves_type, max_seg_curves )
      call rdbk( fileno, seg_curve_table, (max_seg_curves+1)*
     &           max_seg_curve_sets )
      call rd2d( fileno, imatprp, mxmtpr, mxmtpr, mxmat )
      call chk_data_key( fileno, 2, 1 )
      write(out,9010)
c
c                       read logical and character arrays.
c
      read(fileno) convrg, trace                  ! short logical vecs
      call rdbk( fileno, trn, nonode )            ! logical vec
      call rdbk( fileno, stpchk, max_step_limit ) ! logical vec
      call init_maps( fileno, 3 )
c
      read(fileno) lodnam, lodtyp, matnam, elelib ! short char vecs
      read(fileno) smatprp                        ! char array
      call chk_data_key( fileno, 2, 2 )
      write(out,9020)
c
c                       read real arrays.
c
      call rd2d( fileno, props, mxelpr, mxelpr, noelem )
      call rd2d( fileno, matprp, mxmtpr, mxmtpr, mxmat )
      call rdbk( fileno, user_cnstrn_stp_factors, max_step_limit )
      call rdbk( fileno, actual_cnstrn_stp_factors, max_step_limit )
      call chk_data_key( fileno, 3, 0 )
      if ( fgm_node_values_defined ) then
        call mem_allocate( 20 )
        call rd2d( fileno, fgm_node_values, nonode, nonode,
     &             fgm_node_values_cols )
      end if
      call chk_data_key( fileno, 3, 1 )
      write(out,9030)
c
c                       read double precision vectors
c
c
      prec_fact = 2
      call rdbk( fileno, u, prec_fact*nodof )
      call rdbk( fileno, c, prec_fact*nodof )
      call rdbk( fileno, ifv, prec_fact*nodof )
      call rdbk( fileno, load, prec_fact*nodof )
      call mem_allocate( 16 )
      call rdbk( fileno, dload, prec_fact*nodof )
      call rdbk( fileno, rload, prec_fact*nodof )
      call rdbk( fileno, rload_nm1, prec_fact*nodof )
      read(fileno) exist_flag ! for total_user_nodal_forces
      if( exist_flag ) then
        allocate( total_user_nodal_forces(nodof) )
        call rdbk( fileno, total_user_nodal_forces, prec_fact*nodof )
      end if
      call rd2d( fileno, load_pattern_factors, mxlc*prec_fact,
     &           numlod*prec_fact, 2 )
      call mem_allocate( 15 )
      call rdbk( fileno, cnstrn_in, prec_fact*nodof )
      call chk_data_key( fileno, 4, 0 )
      call rdbk( fileno, v, prec_fact*nodof )
      call rdbk( fileno, a, prec_fact*nodof )
      call rdbk( fileno, du, prec_fact*nodof )
      call mem_allocate( 1 )
      call rdbk( fileno, temper_nodes, prec_fact*nonode )
      call rdbk( fileno, temper_nodes_ref, prec_fact*nonode )
      call mem_allocate( 2 )
      call rdbk( fileno, temper_elems, prec_fact*noelem )
      call rd2d( fileno, dmatprp, 2*mxmtpr, 2*mxmtpr, mxmat )
      call chk_data_key( fileno, 4, 1 )
      write(out,9040)
c
c                         blocked data structures for
c                             material point histories and [D]
c                             stresses
c                             material point rotations
c                             strains
c                             element volumes
c
      call history_cep_init( fileno, 2 )
      write(out,9050)
      call stresses_init( fileno, 2 )
      write(out,9055)
      call chk_data_key( fileno, 5, 0 )
      call rotation_init( fileno, 2 )
      write(out,9060)
      call chk_data_key( fileno, 5, 1 )
      call chk_data_key( fileno, 5, 2 )
      call strains_init( fileno, 2 )
      write(out,9080)
      call chk_data_key( fileno, 5, 3 )
      call element_volumes_init( fileno, 2 )
      write(out,9090)
      call chk_data_key( fileno, 5, 4 )
c
c                          coordinate transformation matrices for
c                          node coordinate systems for constraints
c
      do node = 1, nonode
         if (trn(node)) call allo_trnmat(node, 3, fileno)
      end do
      call chk_data_key( fileno, 5, 5 )
c
      call rdbk( fileno, tol, prec_fact*mxcvtests )
      call mem_allocate( 12 )
      call rdbk( fileno, mdiag, prec_fact*nodof )
      call rd3d( fileno, seg_curves, prec_fact*max_seg_points,2,
     &             prec_fact*max_current_pts, 2, max_current_curves )
      call rdbk( fileno, seg_curves_min_stress,
     &           prec_fact*max_current_curves )
      call rdbk( fileno, seg_curves_value,
     &           prec_fact*max_current_curves )
      call rdbk( fileno, seg_curves_ym,
     &           prec_fact*max_current_curves )
      call rdbk( fileno, seg_curves_nu,
     &           prec_fact*max_current_curves )
      call rdbk( fileno, seg_curves_alpha,
     &           prec_fact*max_current_curves )
      call rdbk( fileno, seg_curves_gp_sigma_0,
     &            prec_fact*max_current_curves )
      call rdbk( fileno, seg_curves_gp_h_u,
     &            prec_fact*max_current_curves )
      call rdbk( fileno, seg_curves_gp_beta_u,
     &            prec_fact*max_current_curves )
      call rdbk( fileno, seg_curves_gp_delta_u,
     &            prec_fact*max_current_curves )
      call chk_data_key( fileno, 5, 6 )
      write(out,9100)
c
c                       retrieve the nodal loading data. there are
c                       two data structures and some scalars.
c
      call mem_allocate( 5 )
      read(fileno) read_nod_load
      if ( read_nod_load ) then
        do i = 1, mxlc
         read(fileno) node_count, how_defined
         read(fileno) string
         node_load_defs(i)%node_count = node_count
         node_load_defs(i)%how_defined = how_defined
         node_load_defs(i)%user_file_name = string
         if ( node_count .gt. 0 ) then
          allocate(node_load_defs(i)%nodal_loads(1:node_count,1:2) )
          call rdbk( fileno, node_load_defs(i)%nodal_loads(1,1),
     &               2*node_count )
         end if
        end do
      end if
      call chk_data_key( fileno, 6, 0 )
c
      read(fileno) num_loddat_blks, next_loddat_col
      do i = 1, num_loddat_blks
        if ( i .gt. 1 ) then
         allocate(
     &       loddat_blocks(i)%block(1:mxndldcm,1:sizeof_loddat_blks) )
        end if
        call rdbk( fileno, loddat_blocks(i)%block(1,1),
     &             mxndldcm*sizeof_loddat_blks )
      end do
      call chk_data_key( fileno, 6, 1 )
c
c                       retrieve the element loading data
c
      do i = 1, mxlc
         read(fileno) eload_size
         call do_perm_allo( 1, i, eload_size, .false. )
         if ( eload_size .gt.0 ) then
            call rd2d( fileno, elem_loads(i)%data(1,1),
     &           eload_size, eload_size, 3 )
            call rdbk( fileno, elem_loads(i)%vals(1), eload_size )
            call rdbk( fileno, elem_loads(i)%piston_tabnum(1),
     &           eload_size )
            call eloads_rebuild_thread_list( elem_loads(i)%data(1,1),
     &           eload_size, elem_loads(i)%thread_number(1) )
         end if
      end do
      call chk_data_key( fileno, 6, 2 )
c
c                       retrieve the element equivalent nodal force
c                       vectors for step just completed.
c
      if ( eq_node_force_len .gt. 0 ) then
       call mem_allocate( 23 )
       call rdbk( fileno, eq_node_force_indexes, noelem )
       call rdbk( fileno, eq_node_forces,
     &            eq_node_force_len*prec_fact )
      end if
      call chk_data_key( fileno, 6, 3 )
      write(out,9110)
c
c                       retrieve the nonlinear load step data
c
      if ( .not. allocated( step_load_data ) ) call mem_allocate( 19 )
      do i = 1, max_step_limit
         read(fileno) step_load_data(i)%num_load_patterns
      end do
      do i = 1, max_step_limit
       num_patterns = step_load_data(i)%num_load_patterns
       if ( num_patterns .gt. 0 ) then
        allocate( step_load_data(i)%load_patt_num(num_patterns),
     &            step_load_data(i)%load_patt_factor(num_patterns) )
        read(fileno) step_load_data(i)%load_patt_num(1:num_patterns)
        read(fileno) step_load_data(i)%load_patt_factor(1:num_patterns)
       end if
      end do
      call chk_data_key( fileno, 7, 0 )
c
c                       retrieve the crack growth parameters, if
c                       needed.  Also allocate the variables if they have
c                       been set.
c
c
      read(fileno) crack_growth_type, max_dam_state, num_kill_elem,
     &              csttail, print_status, kill_order, no_killed_elems,
     &              num_print_list, num_kill_order_list,
     &              release_type, crk_pln_normal_idx,
     &              num_crack_plane_nodes, crack_front_start,
     &              crack_front_end, crkfrnt_garbage_start,
     &              crkfrnt_garbage_end,  growth_by_kill,
     &              growth_by_release, min_steps_for_release,
     &              load_size_control_crk_grth,
     &              overshoot_control_crk_grth, overshoot_allocated,
     &              no_released_nodes, g_stp_cntrl_allocated,
     &              const_front, num_crack_fronts, num_nodes_thick,
     &              num_nodes_back, master_lines_set,
     &              num_nodes_grwinc, num_steps_min, load_reduced,
     &              all_elems_killed, num_elements_killed,
     &              enforce_node_release, num_ctoa_released_nodes
      call chk_data_key( fileno, 8, 0 )
c
      read(fileno) porosity_limit, gurson_cell_size,
     &              crack_plane_coord, release_fraction,
     &              critical_angle, char_length, release_height,
     &              crack_plane_sign, init_crit_ang, smcs_alpha,
     &              smcs_beta, CTOA_range, perm_load_fact,
     &              max_porosity_change, max_plast_strain_change,
     &              init_ctoa_dist, ctoa_dist, crkpln_srch_tol,
     &              max_deff_change, critical_cohes_deff_fract,
     &              ppr_kill_displ_fraction
      call chk_data_key( fileno, 8, 1 )
      call allocate_damage( 12 )
      call rdbk( fileno, dam_ptr, noelem )
      call chk_data_key( fileno, 8, 2 )
      write(out,9120)
c
c                       allocate and read arrays for element
c                       extinction
c
      if ( growth_by_kill ) then
         call allocate_damage( 1 )
         call read_damage( 1, fileno, prec_fact )
         if ( print_status ) then
            call allocate_damage( 2 )
            call read_damage( 2, fileno, prec_fact )
         end if
         call chk_data_key( fileno, 9, 0 )
c
         if ( kill_order ) then
            call allocate_damage( 3 )
            call read_damage( 3, fileno, prec_fact )
         end if
         call chk_data_key( fileno, 9, 1 )
c
         if ( .not. no_killed_elems ) then
            call allocate_damage( 4 )
            call read_damage( 4, fileno, prec_fact )
         end if
         call chk_data_key( fileno, 9, 2 )
c
         if ( release_type .eq. 2  ) then
          if ( .not. no_killed_elems ) then
            call allocate_damage( 5 )
            call read_damage( 5, fileno, prec_fact )
          end if
         end if
         call chk_data_key( fileno, 9, 3 )

c
         if ( load_size_control_crk_grth ) then
            call allocate_damage( 9 )
            call read_damage( 9, fileno, prec_fact )
         endif
         call chk_data_key( fileno, 9, 4 )
         write(out,9130)
c
      else if ( growth_by_release ) then
c
         call allocate_damage( 6 )
         call read_damage( 6, fileno, prec_fact )
         call chk_data_key( fileno, 9, 5 )

         if ( release_type .eq.2 ) then
            call allocate_damage( 7 )
            call read_damage( 7, fileno, prec_fact )
         end if
         call chk_data_key( fileno, 9, 6 )
c
         if ( overshoot_control_crk_grth ) then
            call allocate_damage( 8 )
            call read_damage( 8, fileno, prec_fact )
         end if
         call chk_data_key( fileno, 9, 7 )
c
         if ( const_front ) then
            call allocate_damage( 10 )
            call allocate_damage( 11 )
            call read_damage( 10, fileno, prec_fact )
         end if
         call chk_data_key( fileno, 9, 7 )
         write(out,9140)
c
      end if
      call chk_data_key( fileno, 9, 8 )
c
c                       read data structures for contact
c
      read (fileno) use_contact
      call mem_allocate( 24 )
c
      if( use_contact ) then
         call mem_allocate( 25 )
         read (fileno) num_contact
         call rdbk( fileno, contact_shape, num_contact)
         call rdbk( fileno, contact_outside, num_contact)
         call rd2d( fileno, contact_cause, maxcontact, maxcontact,
     &        nonode)
         call rdbk( fileno, contact_stiff, num_contact*prec_fact)
         call rdbk( fileno, contact_depth, num_contact*prec_fact)
         call rdbk( fileno, contact_fric, num_contact*prec_fact)
         call chk_data_key( fileno, 10, 0 )
         call rd2d( fileno, cshape_param, maxcntprm* prec_fact,
     &        maxcntprm*prec_fact, num_contact)
         call rd2d( fileno, cshape_rate, 3*prec_fact, 3*prec_fact,
     &        num_contact)
         call rd2d( fileno, cshape_pnt, 3*prec_fact, 3*prec_fact,
     &        num_contact)
         call rd2d( fileno, cshape_norm, 3*prec_fact, 3*prec_fact,
     &        num_contact)
         call rd3d( fileno, cplane_vec, 3*prec_fact, 2*prec_fact,
     &        3*prec_fact, 2*prec_fact,  num_contact)
         write(out,9150)
      endif
      call chk_data_key( fileno, 10, 1 )
c
c                       read data structures for mpcs
c
      read (fileno) mpcs_exist
      if (mpcs_exist) then
         read (fileno) num_user_mpc
         max_mpc = num_user_mpc
         allocate( user_mpc_table(num_user_mpc), stat=err)
         if (err .ne. 0) then
            call errmsg2(64,dumi,dums,dumr,dumd)
            call die_abort
         end if
         do mpc = 1, num_user_mpc
            read (fileno) ntrm
            user_mpc_table(mpc)%num_terms = ntrm
            allocate ( user_mpc_table(mpc)%node_list(ntrm),
     &                 user_mpc_table(mpc)%dof_list(ntrm),
     &                 user_mpc_table(mpc)%multiplier_list(ntrm),
     &                 stat=err )
            if (err .ne. 0) then
               call errmsg2(64,dumi,dums,dumr,dumd)
               call die_abort
            end if
            read (fileno) user_mpc_table(mpc)%constant
            call rdbk( fileno, user_mpc_table(mpc)%node_list, ntrm )
            call rdbk( fileno, user_mpc_table(mpc)%dof_list, ntrm )
            call rdbk( fileno, user_mpc_table(mpc)%multiplier_list,
     &                    ntrm )
            call chk_data_key( fileno, 11, 1)
         end do
         write(out,9160)
      end if
c
      read (fileno) tied_con_mpcs_constructed
      if (tied_con_mpcs_constructed) then
         read (fileno) num_tied_con_mpc
         max_mpc_tied = num_tied_con_mpc
         allocate( tied_con_mpc_table(num_tied_con_mpc), stat=err)
         if (err .ne. 0) then
            call errmsg2(64,dumi,dums,dumr,dumd)
            call die_abort
         end if
         do mpc = 1, num_tied_con_mpc
            read (fileno) ntrm
            tied_con_mpc_table(mpc)%num_terms = ntrm
            allocate ( tied_con_mpc_table(mpc)%node_list(ntrm),
     &                 tied_con_mpc_table(mpc)%dof_list(ntrm),
     &                 tied_con_mpc_table(mpc)%multiplier_list(ntrm),
     &                 stat=err )
            if (err .ne. 0) then
               call errmsg2(64,dumi,dums,dumr,dumd)
               call die_abort
            end if
            read (fileno) tied_con_mpc_table(mpc)%constant
            call rdbk( fileno, tied_con_mpc_table(mpc)%node_list, ntrm )
            call rdbk( fileno, tied_con_mpc_table(mpc)%dof_list, ntrm )
            call rdbk( fileno, tied_con_mpc_table(mpc)%multiplier_list,
     &                    ntrm )
            call chk_data_key( fileno, 11, 2)
         end do
         write(out,9170)
      end if

      if( mpcs_exist .or. tied_con_mpcs_constructed ) then
        allocate( total_lagrange_forces(nodof) )
        call rdbk( fileno, total_lagrange_forces, nodof* prec_fact )
      end if
      call chk_data_key( fileno, 11, 3)
c
c                       retrieve the table definitions
c
      do itab = 1, max_tables
         read(fileno) tables(itab)%table_name
         read(fileno) tables(itab)%table_type
         read(fileno) num_rows
         read(fileno) num_cols
         call do_table_allo( itab, num_rows, num_cols, .false. )
         if ( num_rows .ne. 0 ) then
            call rd2d( fileno, tables(itab)%table_values_sgl(1,1),
     &           num_rows, num_rows, num_cols )
            call rd2d( fileno, tables(itab)%table_values_dbl(1,1),
     &           num_rows*prec_fact, num_rows*prec_fact, num_cols )
         end if
      end do
      write(out,9180)
      call chk_data_key( fileno, 12, 0 )
c
c                       retrieve the user list definitions, if given.
c
      do ilist = 1, max_user_lists
         read(fileno) user_lists(ilist)%name
         read(fileno) ulen
         user_lists(ilist)%length_list = ulen
         if ( ulen .gt. 0 ) then
           allocate( user_lists(ilist)%list(ulen) )
           call rdbk(fileno,user_lists(ilist)%list(1),ulen)
         endif
      end do
      write(out,9190)
      call chk_data_key( fileno, 12, 1 )
c
c                       read the crystal plasticity crystal data structures
c                       which are handled separately from the material params
c
      call read_alloc_cry( fileno )
      write(out,9200)
      call chk_data_key( fileno, 12, 1)
c
c
c
c
c                       set structure name and flags.
c
      newmas = .false.
      newstf = .false.
c
c                       initialize necessary arrays. cdest, edest,
c                       pbar.
c
      call init_eblock_map
      call cdest_init
      call edest_init
      call solid_cohes_init
      call mem_allocate( 13 )
c
c                       retreive nonlocal material data if it exists.
c                       data structure created by solid_cohes_init.
c                       make sets at n and n+1 the same.
c
      if( nonlocal_analysis ) then
        call chk_data_key( fileno, 13, 0 )
        nsize = nonlocal_shared_state_size
        do i = 1, noelem
           if( .not. nonlocal_flags(i) ) cycle
           read(fileno) nonlocal_data_n(i)%state_values(1:nsize)
           nonlocal_data_n1(i)%state_values(1:nsize) =
     &                  nonlocal_data_n(i)%state_values(1:nsize)
        end do
        write(out,9195)
      end if
c
c                       get convergence history data for user
c                       routine to adjust next step loading,
c                       solution parameters, etc.
c
      call chk_data_key( fileno, 14, 0 )
      do i = 1, 5
        read(fileno) convergence_history(i)%step_converged,
     &         convergence_history(i)%adaptive_used,
     &         convergence_history(i)%iterations_for_convergence,
     &         convergence_history(i)%adapt_substeps
      end do
      read(fileno) run_user_solution_routine
c
c                       get the release_cons_table if we have nodes
c                       undergoing release operations
c
      read(fileno) read_table
      if( read_table ) then
        allocate( release_cons_table(3,nonode) )
        do i = 1, nonode
          read(fileno) release_cons_table(1:3,i)%num_release_steps
          read(fileno) release_cons_table(1:3,i)%
     &                  remaining_steps_for_release
          read(fileno) release_cons_table(1:3,i)%reaction_force
        end do
      end if
      write(out,9210)
c
c                       get output commands file ... information
c                       save file name and bitmap list if it exists
c
      read(fileno) output_command_file
      read(fileno) local_length
      if( local_length > 0 ) then
        allocate( output_step_bitmap_list(local_length ) )
        read(fileno) output_step_bitmap_list
      end if
c
c                       get initial stress data
c
      if( initial_stresses_exist ) then
        allocate( initial_stresses(6,noelem) )
        call rdbk( fileno, initial_stresses, prec_fact*noelem*6 )
        write(out,9230)
      end if
c
c                       USER routines data
c
      call uexternaldb_reopen( fileno, out )
      write(out,9220)
c
c                       read final check variable -- check if correct.
c                       if not, the restored data is corrupted.
c
c
      call chk_data_key( fileno, 15, 0 )
c
c
c                       close the restart file
c
      close( fileno, status='keep' )
      restart_time = wcputime( 1 )
      call errmsg( 194, ltmstp, dbname, restart_time,
     &             total_model_time )
c
c                       if we are running MPI:
c                         send basic data, constraint data, and
c                         analysis parameters data to the slave processors.
c                         also send the element information to the
c                         processor who owns the element.
c                        Allocate the diagonal stiffness for non-MPI
c
      call wmpi_send_basic
      call wmpi_send_const
      call wmpi_send_analysis
      call wmpi_send_reopen
      call wmpi_init_owner
      call wmpi_send_contact (.true.)
      call wmpi_send_crystals
      if ( .not. use_mpi ) call mem_allocate ( 11 )

c
c                       initialize the adaptive algorithm so that
c                       it can operate immediately if needed.
c
      if ( adaptive_flag ) call adaptive_save
c
c                        recalculate mass
c
      call cmpmas
c
c                        open the existing binary file of packet
c                        results if name found in restart file.
c
      call open_packets_file(msg_flag)
c
c                        uexternaldb for Abaqus compatible support
c
      douextdb = 6  ! common.main. tells uexgtern... what to do
      call wmpi_do_uexternaldb
c
9999  sbflg1 = .false.
      sbflg2 = .false.
c
c
 1000 format( 3x, e16.6 )
 9000 format(15x,'> scalars read...')
 9010 format(15x,'> integer arrays read...')
 9020 format(15x,'> logical arrays read...')
 9030 format(15x,'> real arrays read...')
 9040 format(15x,'> double precision arrays read...')
 9050 format(15x,'> material histories and [D] blocks read...')
 9055 format(15x,'> stress blocks read...')
 9060 format(15x,'> material rotation blocks read...')
 9080 format(15x,'> strain blocks read...')
 9090 format(15x,'> element volume blocks read...')
 9100 format(15x,'> material curves read...')
 9110 format(15x,'> node and element loading data read...')
 9120 format(15x,'> crack growth parameters/status read...')
 9130 format(15x,'> element extinction data read...')
 9140 format(15x,'> crack growth by node release data read...')
 9150 format(15x,'> contact data read...')
 9160 format(15x,'> user-defined mpcs read...')
 9170 format(15x,'> tied mesh mpcs read...')
 9180 format(15x,'> table definitions read...')
 9190 format(15x,'> user list definitions read...')
 9195 format(15x,'> nonlocal material data read...')
 9200 format(15x,'> crystal data read...')
 9210 format(15x,'> convergence history read...')
 9220 format(15x,'> user routine data read...')
 9230 format(15x,'> initial stress data read...')
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine rdmass                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/20/2015                 *
c     *                                                              *
c     *     create the blocked data structure for storage of element *
c     *     nodal masses. if requested, read the blocks from the     *
c     *     save file.                                               *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine rdmass
      use global_data ! old common.main
      use elem_block_data, only: mass_blocks
      implicit integer (a-z)
      logical myblk
c
      allocate( mass_blocks(nelblk),stat=iok )
      if( iok .ne. 0 ) then
          call iodevn( idummy, iout, dummy, 1 )
          write(iout,9100) iok
          call die_abort
      end if
c
      do blk = 1, nelblk
c
c                       if we are using MPI:
c                         elblks(2,blk) holds which processor owns the
c                         block.  For slave processors, if we don't own
c                         the block, don't allocate it.  For root,
c                         allocate everything.
c                       if we are using the serial version:
c                         the serial process is root, so allocate all blocks.
c
         myblk = myid .eq. elblks(2,blk)
         if( myid .eq. 0 ) myblk = .true.
         if( .not. myblk ) cycle
c
         felem          = elblks(1,blk)
         num_enodes     = iprops(2,felem)
         num_enode_dof  = iprops(4,felem)
         totdof         = num_enodes * num_enode_dof
         span           = elblks(0,blk)
         allocate( mass_blocks(blk)%ptr(num_enodes,span),stat=iok )
         if( iok .ne. 0 ) then
           call iodevn( idummy, iout, dummy, 1 )
           write(iout,9100) iok
           call die_abort
         end if
      end do
c
      return
c
 9100 format('>> FATAL ERROR: rdmass, memory allocate failure',
     &  /,   '                status= ',i5,
     &  /,   '                job terminated' )
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine history_cep_init             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 04/28/12 rhd               *
c     *                                                              *
c     *     create the blocked data structure for storage of element *
c     *     history data at step n and n+1. also current [D] mats    *
c     *     for use in iter = 0.                                     *
c     *     if requested, read the blocks from the save file.        *
c     *                                                              *
c     ****************************************************************
c
      subroutine history_cep_init( fileno, proc_type )
      use global_data ! old common.main
c
      use elem_block_data, only : history_blocks, history1_blocks,
     &                            history_blk_list, gausspts_blk_list,
     &                            cep_blocks, cep_blk_list
c
      implicit integer (a-z)
      double precision
     &    dummy(1)
      integer matl_info(10)
      logical myblk
c
c
c            proc_type :  = 1, create data arrays only and zero
c                         = 2, create and read from save file
c            fileno:      = save file number. needed only for
c                           proc_type = 2
c
c            the data structure is a vector for each element
c            block (dynamically allocated). the vectors are
c            hung from a dynamically allocated pointer vector
c            of length nelblk.
c
      allocate( history_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 1
           call die_abort
           stop
      end if
c
      allocate( history1_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 2
           call die_abort
           stop
      end if
c
      allocate( history_blk_list(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 5
           call die_abort
           stop
      end if
      allocate( gausspts_blk_list(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 6
           call die_abort
           stop
      end if
c
      allocate( cep_blocks(nelblk), stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 7
           call die_abort
           stop
      end if
      allocate( cep_blk_list(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 8
           call die_abort
           stop
      end if
c
      history_blk_list(1:nelblk)  = 0
      gausspts_blk_list(1:nelblk) = 0
      cep_blk_list(1:nelblk)      = 0
c
c            loop over all element blocks. allocate blocks and
c            zero or read from restart file.
c
      do blk = 1, nelblk
c
c                       if we are using MPI:
c                         elblks(2,blk) holds which processor owns the
c                         block.  For slave processors, if we don't own
c                         the block, don't allocate it.  For root,
c                         allocate everything.
c                       if we are using the serial version:
c                         the serial process is root, so allocate all blocks.
c
        myblk = myid .eq. elblks(2,blk)
        if ( root_processor ) myblk = .true.
        if ( .not. myblk ) cycle
c
        felem      = elblks(1,blk)
        span       = elblks(0,blk)
        ngp        = iprops(6,felem)
        call material_model_info( felem, 0, 1, hist_size )
        call material_model_info( felem, 0, 2, cep_size )
c
c                      process history blocks
c                      ----------------------
c
        block_size = span * ngp * hist_size
c
        allocate( history_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 3
           call die_abort
           stop
        end if
c
        allocate( history1_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 4
           call die_abort
           stop
        end if
c
        if ( proc_type .eq. 1 ) then
          call ro_zero_vec( history_blocks(blk)%ptr(1), block_size )
          call ro_zero_vec( history1_blocks(blk)%ptr(1), block_size )
        end if
c
        if ( proc_type .eq. 2 ) then
               read(fileno) history_blocks(blk)%ptr(1:block_size)
               call chk_data_key( fileno, 200, blk )
               call vec_ops( history1_blocks(blk)%ptr(1),
     &                       history_blocks(blk)%ptr(1),
     &                       dummy, block_size, 5 )
        end if
c
        history_blk_list(blk)  = hist_size
        gausspts_blk_list(blk) = ngp
c
c
c                      process [D] (i.e. cep) blocks
c                      -----------------------------
c
        block_size = span * ngp * cep_size
        allocate( cep_blocks(blk)%vector(block_size),
     &             stat = alloc_stat )
          if ( alloc_stat .ne. 0 ) then
             write(out,9900)
             write(out,9910) 5
             call die_abort
             stop
          end if
c
        if ( proc_type .eq. 1 )
     &         call ro_zero_vec( cep_blocks(blk)%vector(1), block_size )
c
        if ( proc_type .eq. 2 ) then
          read(fileno) cep_blocks(blk)%vector(1:block_size)
          call chk_data_key( fileno, 200, blk )
        end if
c
        cep_blk_list(blk)  = cep_size
c
      end do
      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 history_cep_init: @',i2,/,
     &       '>>> Job terminated....')
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine cdest_init                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/1/98                    *
c     *                                                              *
c     *     create the blocked data structure for storage of element *
c     *     coordinate indexes (cdest). fill integer indexes         *
c     *                                                              *
c     ****************************************************************
c
      subroutine cdest_init
      use global_data ! old common.main
c
      use elem_block_data, only : cdest_blocks, cdest_blk_list
      use main_data,       only : incid, incmap, elems_to_blocks,
     &                            crdmap
c
      implicit integer (a-z)
      integer, dimension (:,:), pointer :: cdest
      logical myblk
c
c            the data structure is a vector for each element
c            block (dynamically allocated). the vectors are
c            hung from a dynamically allocated pointer vector
c            of length nelblk.
c
      allocate( cdest_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 1
           call die_abort
           stop
      end if
c
      allocate( cdest_blk_list(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 5
           call die_abort
           stop
      end if
      cdest_blk_list(1:nelblk) = 0
c
c            loop over all element blocks. allocate blocks.
c
      do blk = 1, nelblk
c
c                       if we are using MPI:
c                         elblks(2,blk) holds which processor owns the
c                         block.  For slave processors, if we don't own
c                         the block, don't allocate it.  For root,
c                         allocate everything.
c                       if we are using the serial version:
c                         the serial process is root, so allocate all blocks.
c
         myblk = myid .eq. elblks(2,blk)
         if (root_processor) myblk = .true.
         if (.not. myblk) cycle
c
         felem         = elblks(1,blk)
         span          = elblks(0,blk)
         nnode         = iprops(2,felem)
         num_enode_dof = iprops(4,felem)
         totdof        = nnode * num_enode_dof
c
         allocate( cdest_blocks(blk)%ptr(totdof,span),
     &             stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
            write(out,9900)
            write(out,9910) 3
            call die_abort
            stop
         end if
         cdest_blk_list(blk) = 1
      end do
c
c            loop over all elements. create the coordinate
c            index pointers.
c
      do elem = 1, noelem
c
c                       if we are using MPI:
c                         elblks(2,blk) holds which processor owns the
c                         block.  For slave processors, don't compute
c                         the terms unless we own the block. For root,
c                         compute all the terms.
c                       if we are using the serial version:
c                         the serial process is root, so allocate all blocks.
c
         blk            = elems_to_blocks(elem,1)
         myblk = myid .eq. elblks(2,blk)
         if (root_processor) myblk = .true.
         if (.not. myblk) cycle
c
         rel_elem       = elems_to_blocks(elem,2)
         nnode          = iprops(2,elem)
         incptr         = incmap(elem)
         cdest => cdest_blocks(blk)%ptr
         do j = 1, nnode
          cdest(j,rel_elem)         = crdmap(incid(incptr+j-1))
          cdest(nnode+j,rel_elem)   = crdmap(incid(incptr+j-1)) + 1
          cdest(2*nnode+j,rel_elem) = crdmap(incid(incptr+j-1)) + 2
         end do
      end do

      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 cdest_init: @',i2,/,
     &       '>>> Job terminated....')
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine edest_init                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/1/98                    *
c     *                                                              *
c     *     create the blocked data structure for storage of element *
c     *     displacement indexes (cdest). fill integer indexes       *
c     *                                                              *
c     ****************************************************************
c
      subroutine edest_init
      use global_data ! old common.main
c
      use elem_block_data, only : edest_blocks, edest_blk_list
      use main_data,       only : incid, incmap, elems_to_blocks
c
      implicit integer (a-z)
c
      integer, dimension (:,:), pointer :: edest
      logical myblk
c
c            the data structure is a vector for each element
c            block (dynamically allocated). the vectors are
c            hung from a dynamically allocated pointer vector
c            of length nelblk.
c
      allocate( edest_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 1
           call die_abort
           stop
      end if
c
      allocate( edest_blk_list(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 5
           call die_abort
           stop
      end if
      edest_blk_list(1:nelblk) = 0
c
c            loop over all element blocks. allocate blocks.
c
      do blk = 1, nelblk
c
c                       if we are using MPI:
c                         elblks(2,blk) holds which processor owns the
c                         block.  For slave processors, if we don't own
c                         the block, don't allocate it.  For root,
c                         allocate everything.
c                       if we are using the serial version:
c                         the serial process is root, so allocate all blocks.
c
         myblk = myid .eq. elblks(2,blk)
         if (root_processor) myblk = .true.
         if (.not. myblk) cycle
c
         felem         = elblks(1,blk)
         span          = elblks(0,blk)
         nnode         = iprops(2,felem)
         num_enode_dof = iprops(4,felem)
         totdof        = nnode * num_enode_dof
c
         allocate( edest_blocks(blk)%ptr(totdof,span),
     &             stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
            write(out,9900)
            write(out,9910) 3
            call die_abort
            stop
         end if
         edest_blk_list(blk) = 1
      end do
c
c            loop over all elements. create the displacement
c            index pointers.
c
      do elem = 1, noelem
c
c                       if we are using MPI:
c                         elblks(2,blk) holds which processor owns the
c                         block.  For slave processors, don't compute
c                         the terms unless we own the block. For root,
c                         compute all the terms.
c                       if we are using the serial version:
c                         the serial process is root, so allocate all blocks.
c
         blk            = elems_to_blocks(elem,1)
         myblk = myid .eq. elblks(2,blk)
         if (root_processor) myblk = .true.
         if (.not. myblk) cycle
c
         rel_elem       = elems_to_blocks(elem,2)
         nnode          = iprops(2,elem)
         incptr         = incmap(elem)
         edest => edest_blocks(blk)%ptr
         do j = 1, nnode
          edest(j,rel_elem)         = dstmap(incid(incptr+j-1))
          edest(nnode+j,rel_elem)   = dstmap(incid(incptr+j-1)) + 1
          edest(2*nnode+j,rel_elem) = dstmap(incid(incptr+j-1)) + 2
         end do
      end do
c
      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 edest_init: @',i2,/,
     &       '>>> Job terminated....')
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine strains_init                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/8/98                    *
c     *                                                              *
c     *     create the blocked data structure for storage of element *
c     *     strain data at step n and n+1. if requested, read the    *
c     *     blocks from the save file.                               *
c     *                                                              *
c     ****************************************************************
c
      subroutine strains_init( fileno, proc_type )
      use global_data ! old common.main
c
      use elem_block_data, only : eps_n_blocks, eps_n1_blocks,
     &                            eps_blk_list
c
      implicit integer (a-z)
      double precision
     &    dummy(1)
      logical myblk
c
c
c            proc_type :  = 1, create data arrays only and zero
c                         = 2, create and read from save file
c            fileno:      = save file number. needed only for
c                           proc_type = 2
c
c            the data structure is a vector for each element
c            block (dynamically allocated). the vectors are
c            hung from a dynamically allocated pointer vector
c            of length nelblk.
c
      allocate( eps_n_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 1
           call die_abort
           stop
      end if
c
      allocate( eps_n1_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 2
           call die_abort
           stop
      end if
c
      allocate( eps_blk_list(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 5
           call die_abort
           stop
      end if
      eps_blk_list(1:nelblk) = 0
c
c            loop over all element blocks. allocate blocks and
c            zero or read from restart file.
c
      do blk = 1, nelblk
c
c                       if we are using MPI:
c                         elblks(2,blk) holds which processor owns the
c                         block.  For slave processors, if we don't own
c                         the block, don't allocate it.  For root,
c                         allocate everything.
c                       if we are using the serial version:
c                         the serial process is root, so allocate all blocks.
c
        myblk = myid .eq. elblks(2,blk)
        if (root_processor) myblk = .true.
        if (.not. myblk) cycle
c
        felem      = elblks(1,blk)
        span       = elblks(0,blk)
        ngp        = iprops(6,felem)
        block_size = span * ngp * nstr
c
        allocate( eps_n_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 3
           call die_abort
           stop
        endif
c
        allocate( eps_n1_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 4
           call die_abort
           stop
        endif
c
        if ( proc_type .eq. 1 ) then
          call zero_vector( eps_n_blocks(blk)%ptr(1), block_size )
          call zero_vector( eps_n1_blocks(blk)%ptr(1), block_size )
        end if
        if ( proc_type .eq. 2 ) then
               read(fileno) eps_n_blocks(blk)%ptr(1:block_size)
               call chk_data_key( fileno, 360, blk )
               call vec_ops( eps_n1_blocks(blk)%ptr(1),
     &                       eps_n_blocks(blk)%ptr(1), dummy,
     &                       block_size, 5 )
        end if
c
        eps_blk_list(blk) = 1
c
      end do
      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 strains_init: @',i2,/,
     &       '>>> Job terminated....')
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine stresses_init                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/8/98                    *
c     *                                                              *
c     *     create the blocked data structure for storage of element *
c     *     stress data at step n and n+1. if requested, read the    *
c     *     blocks from the save file.                               *
c     *                                                              *
c     ****************************************************************
c
      subroutine stresses_init( fileno, proc_type )
      use global_data ! old common.main
c
      use elem_block_data, only : urcs_n_blocks, urcs_n1_blocks,
     &                            urcs_blk_list
c
      implicit integer (a-z)
      double precision
     &    dummy(1)
      logical myblk
c
c
c            proc_type :  = 1, create data arrays only and zero
c                         = 2, create and read from save file
c            fileno:      = save file number. needed only for
c                           proc_type = 2
c
c            the data structure is a vector for each element
c            block (dynamically allocated). the vectors are
c            hung from a dynamically allocated pointer vector
c            of length nelblk.
c
      allocate( urcs_n_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 1
           call die_abort
           stop
      end if
c
      allocate( urcs_n1_blocks(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 2
           call die_abort
           stop
      end if
c
      allocate( urcs_blk_list(nelblk),stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 5
           call die_abort
           stop
      end if
      urcs_blk_list(1:nelblk) = 0
c
c            loop over all element blocks. allocate blocks and
c            zero or read from restart file.
c
      do blk = 1, nelblk
c
c                       if we are using MPI:
c                         elblks(2,blk) holds which processor owns the
c                         block.  For slave processors, if we don't own
c                         the block, don't allocate it.  For root,
c                         allocate everything.
c                       if we are using the serial version:
c                         the serial process is root, so allocate all blocks.
c
        myblk = myid .eq. elblks(2,blk)
        if (root_processor) myblk = .true.
        if (.not. myblk) cycle
c
        felem      = elblks(1,blk)
        span       = elblks(0,blk)
        ngp        = iprops(6,felem)
        block_size = span * ngp * nstrs
c
        allocate( urcs_n_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 3
           call die_abort
           stop
        endif
c
        allocate( urcs_n1_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 4
           call die_abort
           stop
        endif
c
        if ( proc_type .eq. 1 ) then
         call zero_vector( urcs_n_blocks(blk)%ptr(1), block_size )
          call zero_vector( urcs_n1_blocks(blk)%ptr(1), block_size )
        end if
        if ( proc_type .eq. 2 ) then
               read(fileno) urcs_n_blocks(blk)%ptr(1:block_size)
               call chk_data_key( fileno, 370, blk )
               call vec_ops( urcs_n1_blocks(blk)%ptr(1),
     &                       urcs_n_blocks(blk)%ptr(1), dummy,
     &                       block_size, 5 )
        end if
c
        urcs_blk_list(blk) = 1
c
      end do
      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 stresses_init: @',i2,/,
     &       '>>> Job terminated....')
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine element_volumes_init         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 12/24/99                   *
c     *                                                              *
c     *     create the blocked data structure for storage of element *
c     *     volumes. if requested, read the blocks from the          *
c     *     save file.                                               *
c     *                                                              *
c     ****************************************************************
c
      subroutine element_volumes_init( fileno, proc_type )
      use global_data ! old common.main
c
      use elem_block_data, only : element_vol_blocks
c
      implicit integer (a-z)
      double precision
     &    model_volume, zero
      logical myblk
      data zero / 0.0 /
c
c
c            proc_type :  = 1, create data arrays only and zero
c                         = 2, create and read from save file
c            fileno:      = save file number. needed only for
c                           proc_type = 2
c
c            the data structure is a vector for each element
c            block (dynamically allocated). the vectors are
c            hung from a dynamically allocated pointer vector
c            of length nelblk.
c
      allocate( element_vol_blocks(nelblk), stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 1
           call die_abort
           stop
      end if
c
c            loop over all element blocks. allocate blocks and
c            zero or read from restart file.
c
      model_volume = zero
      do blk = 1, nelblk
c
c                       if we are using MPI:
c                         elblks(2,blk) holds which processor owns the
c                         block.  For slave processors, if we don't own
c                         the block, don't allocate it.  For root,
c                         allocate everything.
c                       if we are using the serial version:
c                         the serial process is root, so allocate all blocks.
c
        myblk = myid .eq. elblks(2,blk)
        if( root_processor ) myblk = .true.
        if( .not. myblk ) cycle
c
        felem      = elblks(1,blk)
        span       = elblks(0,blk)
        block_size = span
c
        allocate( element_vol_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 3
           call die_abort
           stop
        endif
c
        if ( proc_type .eq. 1 ) then
          call zero_vector( element_vol_blocks(blk)%ptr(1), block_size )
        end if
        if ( proc_type .eq. 2 ) then
           read(fileno) element_vol_blocks(blk)%ptr(1:block_size)
           call chk_data_key( fileno, 380, blk )
        end if
c
      end do

      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 element_volumes_init: @',i2,/,
     &       '>>> Job terminated....')
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine rotation_init                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/1/98                    *
c     *                                                              *
c     *     create the blocked data structure for storage of element *
c     *     rotation matrices at gauss points to support finite      *
c     *     strains/large displacements. allocate space for those    *
c     *     blocks which require space. initialize or read from      *
c     *     restart file                                             *
c     *                                                              *
c     ****************************************************************
c
      subroutine rotation_init( fileno, proc_type )
      use global_data ! old common.main
c
      use elem_block_data, only : rot_n_blocks, rot_n1_blocks,
     &                            rot_blk_list
      implicit integer (a-z)
      double precision
     &    zero, rot_init(9), dummy(1)
      logical geo_non_flg, local_debug
      data local_debug, zero, one / .false., 0.0, 1.0 /
      logical myblk
c
c
c            proc_type :  = 1, create data arrays only and zero
c                         = 2, create and read from save file
c            fileno:      = save file number. needed only for
c                           proc_type = 2
c
c            the data structure is a vector for each element
c            block (dynamically allocated). the vectors are
c            hung from a dynamically allocated pointer vector
c            of length nelblk.
c
c            determine if any blocks of elements have finite strains/
c            large displacements. if so, we have to create the
c            blocked data structures for material point rotation
c            matrices. for a completely small strain analysis,
c            nothing gets done here.
c
      if ( .not. allocated(rot_blk_list) ) then
         allocate( rot_blk_list(nelblk),stat=alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 7
           call die_abort
           stop
          end if
          rot_blk_list(1:nelblk) = 0
      end if
c
      do blk = 1, nelblk
        felem       = elblks(1,blk)
        geo_non_flg = lprops(18,felem)
        if ( geo_non_flg ) go to 100
      end do
c
c            no blocks with geometric nonlinear elements found.
c            this is a small displacement solution. just leave.
c
      return
c
 100  continue
      allocate( rot_n_blocks(nelblk), stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 1
           call die_abort
           stop
      end if
      allocate( rot_n1_blocks(nelblk), stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 2
           call die_abort
           stop
      end if
c
c            loop over all element blocks. allocate blocks which
c            contain geometric nonlinear elements. initialize with
c            unit matrices or read from restart file.
c
      rot_init(1:9) = zero
      rot_init(1)   = one
      rot_init(5)   = one
      rot_init(9)   = one
c
c
      do blk = 1, nelblk
c
c                       if we are using MPI:
c                         elblks(2,blk) holds which processor owns the
c                         block.  For slave processors, if we don't own
c                         the block, don't allocate it.  For root,
c                         allocate everything.
c                       if we are using the serial version:
c                         the serial process is root, so allocate all blocks.
c
        myblk = myid .eq. elblks(2,blk)
        if (root_processor) myblk = .true.
        if (.not. myblk) cycle
c
        felem       = elblks(1,blk)
        geo_non_flg = lprops(18,felem)
        if ( .not. geo_non_flg ) cycle
        span       = elblks(0,blk)
        ngp        = iprops(6,felem)
        block_size = span * ngp * 9
c
        allocate( rot_n_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 4
           call die_abort
           stop
        endif
c
        allocate( rot_n1_blocks(blk)%ptr(block_size),
     &            stat = alloc_stat )
        if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 5
           call die_abort
           stop
        endif
c
        rot_blk_list(blk) = 1
c
        if ( proc_type .eq. 1 ) then
          low = 1
          up  = 9
          do i = 1, span
             do j = 1, ngp
                rot_n_blocks(blk)%ptr(low:up) = rot_init(1:9)
                low = low + 9
                up  = up + 9
             end do
          end do
          call vec_ops( rot_n1_blocks(blk)%ptr(1),
     &                  rot_n_blocks(blk)%ptr(1), dummy, block_size, 5 )
        end if
c
        if ( proc_type .eq. 2 ) then
               read(fileno) rot_n1_blocks(blk)%ptr
               call chk_data_key( fileno, 300, blk )
               rot_n_blocks(blk)%ptr = rot_n1_blocks(blk)%ptr
        end if
c
      end do
      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 rotation_init: @',i2,/,
     &       '>>> Job terminated....')
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine init_maps                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 2/16/2016 rhd              *
c     *                                                              *
c     *     allocate space for various mapping data structures.      *
c     *     if reastarting, read data from the restart file.         *
c     *                                                              *
c     ****************************************************************
c
      subroutine init_maps( fileno, proc_type )
      use global_data ! old common.main
      use main_data, only : invdst, repeat_incid,
     &                      inverse_incidences, inverse_dof_map
      implicit integer (a-z)
c
c
c            fileno:      = file number to read to recover data
c                           proc_type = 0, do not read
c
      go to ( 100, 200, 300, 400, 500, 600, 700 ), proc_type
c
c            available
c
 100  continue
      return
c
c            inverse dest mapping
c
 200  continue
      allocate( invdst(nodof), stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 3
           call die_abort
           stop
      end if
      if ( fileno .eq. 0 ) return
      call rdbk( fileno, invdst, nodof )
      return
c
c            logical vector indicating if element has repeated nodes
c            in the incidences, e.g., collapsed elements used in
c            crack modeling.
c
 300  continue
      allocate( repeat_incid(noelem), stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 4
           call die_abort
           stop
      end if
      if ( fileno .eq. 0 ) return
      call rdbk( fileno, repeat_incid, noelem)
      return
c
c            inverse incidences: list of elements connected to
c            each model node. allocate the data structure for each node.
c            read values from restart file if requested
c
 400  continue
      allocate( inverse_incidences(nonode), stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 5
           call die_abort
           stop
      end if
      if ( fileno .eq. 0 ) return
      read(fileno) (inverse_incidences(i)%element_count,i=1,nonode)
      return
c
c            inverse incidences: list of elements connected to
c            each model node. allocate the list vector for each node
c            read values from researt file if requested.
c
 500  continue
      do i = 1, nonode
       ecount = inverse_incidences(i)%element_count
       allocate( inverse_incidences(i)%element_list(ecount),
     &           stat=alloc_stat )
       if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 6
           call die_abort
           stop
       end if
       if ( fileno .ne. 0 ) then
         read(fileno) (inverse_incidences(i)%element_list(j),
     &                 j = 1, ecount )
       end if
      end do
      return
c
c            inverse dof maps for elements connected to nodes.
c            allocate the data structure for each node.
c            no values to read from file
c
 600  continue
      allocate( inverse_dof_map(nonode), stat=alloc_stat )
      if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 7
           call die_abort
           stop
      end if
      return
c
c            inverse dof maps for elements connected to nodes.
c            allocated the 2-d array for each model node.
c            read each array from disk if required
c
 700  continue
      do stnd = 1, nonode
       ecount = inverse_incidences(stnd)%element_count
       allocate( inverse_dof_map(stnd)%edof_table(ecount,3),
     &           stat=alloc_stat )
       if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 7
           call die_abort
           stop
       end if
       if ( fileno .ne. 0 ) then
         read(fileno) ((inverse_dof_map(stnd)%edof_table(j,k),
     &                 j = 1, ecount ), k = 1, 3 )
       end if
      end do
      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 init_maps: @',i2,/,
     &       '>>> Job terminated....' )
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine init_eblock_map              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/1/98                    *
c     *                                                              *
c     *     create a data structure that provides the block number   *
c     *     for an element and the relative element within the block *
c     *     column 1 provides the block number, column 2 the         *
c     *     relative element within the block                        *
c     *                                                              *
c     ****************************************************************
c
      subroutine init_eblock_map
      use global_data ! old common.main
      use main_data
      implicit integer (a-z)

      if ( allocated( elems_to_blocks ) ) then
           write(out,9900)
           write(out,9910) 1
           call die_abort
           stop
      end if
c
      allocate( elems_to_blocks(noelem,2), stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           write(out,9910) 2
           call die_abort
           stop
         end if
c
      do blk = 1, nelblk
       span  = elblks(0,blk)
       felem = elblks(1,blk)
       do i = 1, span
         elems_to_blocks(felem+i-1,1) = blk
         elems_to_blocks(felem+i-1,2) = i
       end do
      end do
c
      return
c
 9900 format('>>> FATAL ERROR: memory allocate failure...')
 9910 format('                 init_elblock_map: @',i2,/,
     &       '>>> Job terminated....' )
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *               subroutine allo_trnmat                         *
c     *                                                              *
c     *                    written by : asg                          *
c     *                                                              *
c     *                last modified : 7/1/98                        *
c     *                                                              *
c     *     this subroutine provides the allocation/deallocation     *
c     *     of node transformation matricies. If called during a     *
c     *     restart, this routine also drives the reading in of      *
c     *     the saved data.                                          *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine allo_trnmat( node, status, fileno )
      use global_data ! old common.main
c
      use main_data
c
      implicit integer (a-z)
      double precision
     &  zero
      data zero / 0.0 /
c
c              status values:
c                1 -- allocate trnmat(node)%mat
c                2 -- deallocate trnmat(node)%mat
c                3 -- allocate trnmat(node)%mat and read from restart file
c
c
c              allocate 3x3 rotation matrix at node to
c              support non-global constraints.
c
      if ( status .eq. 1 .or. status .eq. 3 ) then
c
         allocate( trnmat(node)%mat(3,3), stat = alloc_stat )
         if ( alloc_stat .ne. 0 ) then
           write(out,9900)
           call die_abort
           stop
         end if
c
         if ( status .eq. 3 ) then
            read (fileno) trnmat(node)%mat
         endif
c
         return
c
c              deallocate 3x3 rotation matrix at node if it exits.
c
      else if ( status .eq. 2 ) then
c
         deallocate( trnmat(node)%mat )
c
      endif
c
 9900 format('>>> FATAL ERROR: memory allocate failure -- nodal trnmat')
c
      end

c
c     ****************************************************************
c     *                                                              *
c     *               subroutine chk_data_key                        *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *                last modified : 2/20/01                       *
c     *                                                              *
c     *     check that next record of restart file is a data key     *
c     *     if not, terminate. caller passes level and sub-level     *
c     *     of where we are in the restart process                   *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine chk_data_key( fileno, level1, level2 )
      implicit integer (a-z)
c
      data check_data_key / 2147483647 /
c
      read (fileno,end=100) checkvar
c
      if ( checkvar .ne. check_data_key ) then
         write(*,9000) level1, level2
         call die_abort
      end if
      return
c
 100  write(*,9100) level1, level2
      call die_abort

c
      return
 9000 format(
     & //,">  FATAL ERROR: unexpected event while reading",
     &  /,"                restart file. the check for embedded data",
     &  /,"                keys failed. point of failure is",
     &  /,"                level: ",i4," sub-level: ",i4,
     &  /,"                job terminated...." ///)
 9100 format(
     & //,">  FATAL ERROR: unexpected end-of-file while reading",
     &  /,"                restart file. the check for embedded data",
     &  /,"                keys failed. point of failure is",
     &  /,"                level: ",i4," sublevel: ",i4,
     &  /,"                job terminated...." ///)
      end


c     ****************************************************************
c     *                                                              *
c     *               subroutine ro_zero_vec                         *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *                last modified : 04/28/12                      *
c     *                                                              *
c     *            zero a vector. written to be inlined  everywhere  *
c     *            in this .f                                        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ro_zero_vec( vec, isize )
      implicit none
c
      double precision
     &  vec(*), zero
c
      integer isize
      data zero / 0.0d00 /
c
      vec(1:isize) = zero
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *               subroutine ro_copy_vec                         *
c     *                                                              *
c     *                    written by : rhd                          *
c     *                                                              *
c     *                last modified : 02/27/13                      *
c     *                                                              *
c     *        copy one vector into another. use integer type        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ro_copy_vec( vec_left, vec_right, isize )
      implicit none
c
      integer vec_left(*), vec_right(*), isize
c
      vec_left(1:isize) = vec_right(1:isize)
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                   subroutine solid_cohes_init                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/27/2013 rhd             *
c     *                                                              *
c     *     build data structure that gives the 2 solid elements     *
c     *     connected to each cohesive-interface element             *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine solid_cohes_init
      use global_data ! old common.main
c
      use main_data, only : cohesive_ele_types, nonlocal_analysis
c
      use elem_block_data, only : solid_interface_lists,
     &                            nonlocal_flags, nonlocal_data_n,
     &                            nonlocal_data_n1
c
      implicit integer (a-z)
c
      double precision
     &  zero
      data zero / 0.0 /
c
      logical local_debug
c
c          global nonlocal flag must have be set in user input to
c          include all data structures for nonlocal analyses of
c          cohesive elements
c
      local_debug = .false.
      if( local_debug ) write(out,9000)
      if( .not. nonlocal_analysis ) return
c
c          find first block with cohesive elements. may not have
c          any work to do
c
      do blk = 1, nelblk
        felem = elblks(1,blk)
        eletype = iprops(1,felem)
        if( cohesive_ele_types(eletype) ) go to 100
      end do
c
      if( nonlocal_analysis ) then   ! inconsistency in data values
        write(iout,9300)
        call die_abort
      end if
      return
c
c          for cohesive-interface elements, build the
c          table of the 2 solid elements connected to
c          each cohesive-interface. Data structure:
c          vector 1:num blocks. each entry
c          is a dynamically allocated array of span rows x
c          2. cols 1 & 2 contain number solid element connected to
c          top & bottom surface.
c
c          then build the nonlocal material data used to support
c          cohesive elements. material models for solid elements are
c          allowed to share a set of state variables with cohesive
c          elements to support nonlocal cohesive models
c
 100  continue
      if( local_debug ) write(out,9010) blk
      allocate( solid_interface_lists(nelblk),
     &          nonlocal_flags(noelem), nonlocal_data_n(noelem),
     &          nonlocal_data_n1(noelem) )
      nonlocal_flags(1:noelem) = .false.
c
c          algorithm:
c            o - loop over all blocks
c            o - skip solid element blocks
c            o - allocate array for block to
c                store 2 solid element numbers for each
c                cohesive element in block
c            o - loop over each cohesive element in block
c            o - use lower level routine for details of incidence
c                matching (actually v. fast)
c            o - save solid element number connected
c                to top and bottom surface.
c
c            Block loop probably not worth running parallel
c            since the number of blocks with cohesive-interface
c            elements likely not large. And only integer
c            computations for setting up solid-interface
c            topology data structure.
c
      do blk = 1, nelblk
        felem = elblks(1,blk)
        span = elblks(0,blk)
        eletype = iprops(1,felem)
        if( .not. cohesive_ele_types(eletype) ) cycle
        if( local_debug ) write(out,9020) blk, felem, span,
     &         cohesive_ele_types(eletype)
        allocate( solid_interface_lists(blk)%list(span,2) )
        nnodes = iprops(2,felem) ! for cohes type in block
        do i = 1, span
          elem_cohes = felem + i - 1
          if( local_debug ) write(out,9030) i, elem_cohes, nnodes
          call solid_cohes_init_a( iprops, mxelpr,
     &            elem_cohes, bottom_solid, top_solid, out )
          solid_interface_lists(blk)%list(i,1) = bottom_solid
          solid_interface_lists(blk)%list(i,2) = top_solid
          if( bottom_solid .gt. 0 )
     &           nonlocal_flags(bottom_solid) = .true.
          if( top_solid .gt. 0 ) nonlocal_flags(top_solid) = .true.
        end do  ! cohes elems in block
      end do ! next block of elements
c
c          build the nonlocal material data used to support
c          cohesive elements. material models for solid elements are
c          allowed to share a set of state variables with cohesive
c          elements to support nonlocal cohesive models.
c
c          any solid element & material model combination is
c          allowed to share state type data with cohesive elements.
c          the cohesive material models must be aware of the contents
c          created by the material model for the solid elements -
c          WARP3D just manages the space.
c
      if( local_debug ) write(out,9300)
      nsize = nonlocal_shared_state_size
      count = 0
      do i = 1, noelem
         if( .not. nonlocal_flags(i) ) cycle
         allocate( nonlocal_data_n(i)%state_values(nsize),
     &          nonlocal_data_n1(i)%state_values(nsize) )
         nonlocal_data_n(i)%state_values(1:nsize) = zero
         nonlocal_data_n1(i)%state_values(1:nsize) = zero
         count = count + 1
      end do
      if( local_debug ) write(out,9350) count
c
c          dump found solid element connected to each cohesive
c          element
c
      if( local_debug ) then
        write(out,*) ' '
        write(out,*) '... summary of solids connected to ',
     &    'cohesive-interface elements'
        do blk = 1, nelblk
         if( .not. allocated( solid_interface_lists(blk)%list ) ) cycle
         felem   = elblks(1,blk)
         span    = elblks(0,blk)
         eletype = iprops(1,felem)
         nnodes  = iprops(2,felem) ! for cohes type in block
         write(out,9100) blk, span, eletype, nnodes
         do i = 1, span
           write(out,9110) i, felem+i-1,
     &      solid_interface_lists(blk)%list(i,1) ,
     &      solid_interface_lists(blk)%list(i,2)
         end do  ! cohes elems in block
         write(out,*) ' '
        end do ! next block of elements
      end if
c
      return
c
 9000 format(/,' ... inside solid_cohes_init ...')
 9010 format(/,'... first block w/ cohesive-interface elements:',i5)
 9020 format(/,'... processing blk, felem, span, is_cohes:',
     &    i5, i8, i4, l3 )
 9030 format(5x,'i, elem_cohes, nnodes:',i4,i8,i3)
 9100 format(5x,'processing blk, span, eletype, nnodes:',i4,i4,i3,i3)
 9110 format(8x,'i, elem, bottom & top elements:',i4,i8,2i9)
 9200 format(5x,'blk: ',i5,' no cohesive-interface elements')
 9300 format(/,'... build nonlocal material data structure ...')
 9310 format(5x,'processing blk, span, eletype:',
     &       i4,i4,i3,i3)
 9320 format(5x,'number of umat blocks created for nonlocal: ',i5)
 9350 format(5x,'number of model elements with nonlocal state vector:',
     &   i8 )
c
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine solid_cohes_init_a                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/8/13 rhd                *
c     *                                                              *
c     *     process at single cohesive element to get solid element  *
c     *     connected to top and bottom surface.                     *
c     *                                                              *
c     *     Note: there may not be a top or bottom solid element.    *
c     *       - cohesive elements connected to a symmetry plane      *
c     *       - a model with only cohesive elements                  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine solid_cohes_init_a( iprops, mxelpr,
     &              elem, bottom_solid, top_solid, iout )
c
      use main_data, only : incid, incmap, inverse_incidences,
     &                      cohesive_ele_types
c
      implicit none
      integer :: mxelpr, elem, bottom_solid, top_solid, iout
      integer iprops(mxelpr,*)
c
c          elem = absolute number for a cohesive-interface element
c
c          return bottom_solid, top_solid element numbers
c            (one or both could be zero)
c
c          locals -- sizes are plenty big
c
      integer cohes_nodes(30), top_nodes(15), bottom_nodes(15)
      integer nnodes, i, nface, coh_node, num_connected,
     &        connected_elem, eletype, num_node_solid
      logical solid_cohes_init_b ! supporting search function
      logical local_debug
c
c
      nnodes = iprops(2,elem)  ! on cohesive element
      do i = 1, nnodes
        cohes_nodes(i) = incid(incmap(elem)+i-1)
      end do
c
      nface = nnodes / 2  ! nnodes always multiple of 2
      do i = 1, nface
        bottom_nodes(i) = cohes_nodes(i)
        top_nodes(i)    = cohes_nodes(nface+i)
      end do
      local_debug = .false.
      if( local_debug ) then
        write(iout,*) ' ';write(iout,*) ' '
        write(iout,*) '.. solid_cohes_init_a, nface:, ',nface
        write(iout,*) '  bottom_nodes, top_nodes: '
          write(iout,*) '    ', bottom_nodes(1:nface)
          write(iout,*) '    ', top_nodes(1:nface)
       end if
c
      bottom_solid = 0
      top_solid    = 0
c
c          find solid element connected to bottom surface
c          algorithm:
c            o - use the list of elements connected to node 1
c                of the cohesive element. the solid element we
c                want must be in that list
c            o - loop over list of elements connected to node 1
c                of cohesive element
c            o - get element number, type. skip if is another
c                cohesive element (could be the cohesive element
c                we are processing). get number of nodes on the
c                solid element
c            o - use support routine that compares incidences of
c                this solid element and nodes on surface of
c                this cohesive element
c
c
c          Consistency check:
c            there can be only 1 solid element connected to
c            surface of cohesive element. otherwise there is an
c            issue in the mesh definition
c
c          A cohesive element with no solid connected is ok.
c
      coh_node      = bottom_nodes(1)
      num_connected = inverse_incidences(coh_node)%element_count
      if( local_debug ) then
         write(iout,*) '    coh_node, num_connected:',coh_node,
     &          num_connected
         write(iout,*) '    connected element loop....'
      end if
c
      do i = 1, num_connected
        connected_elem = inverse_incidences(coh_node)%element_list(i)
        eletype = iprops(1,connected_elem)
        if( local_debug) then
           write(iout,9900) i, connected_elem, eletype,
     &             cohesive_ele_types(eletype)
        end if
        if( cohesive_ele_types(eletype) ) cycle
        num_node_solid = iprops(2,connected_elem)
        if( solid_cohes_init_b( nface, bottom_nodes, num_node_solid,
     &       incid(incmap(connected_elem)) ) ) then
          if( bottom_solid .eq. 0 ) then
              bottom_solid = connected_elem
          else
             write(iout,9000) elem
             call die_abort
          end if
        end if
      end do  ! elements connect to cohesive node
c
c          find solid on top surface. Same algorithm as above.
c
      coh_node      = top_nodes(1)
      num_connected = inverse_incidences(coh_node)%element_count
      do i = 1, num_connected
        connected_elem = inverse_incidences(coh_node)%element_list(i)
        eletype = iprops(1,connected_elem)
        if( cohesive_ele_types(eletype) ) cycle
        num_node_solid = iprops(2,connected_elem)
        if( solid_cohes_init_b( nface, top_nodes, num_node_solid,
     &       incid(incmap(connected_elem)) ) ) then
          if( top_solid .eq. 0 ) then
              top_solid = connected_elem
          else
             write(iout,9010) elem
             call die_abort
          end if
        end if
      end do  ! elements connect to cohesive node
c
      return
c
 9000 format('>>>>>> FATAL ERROR: ',
     &/,8x,'cohesive element: ',i8,' appears to have more',
     &/,8x,'than one solid element attached to the bottom surface.',
     &/,8x,'WARP3D data structures are not designed for this case.',
     &/,8x,'Terminating analysis at this point...',//)
c
 9010 format('>>>>>> FATAL ERROR: ',
     &/,8x,'cohesive element: ',i8,' appears to have more',
     &/,8x,'than one solid element attached to the top surface.',
     &/,8x,'WARP3D data structures are not designed for this case.',
     &/,8x,'Terminating analysis at this point...',//)
c
 9900  format(8x,'i, contd ele, etype, iscoh: ', 3i8,l2)
c
      end
c     ****************************************************************
c     *                                                              *
c     *               function solid_cohes_init_b                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 01/8/13 rhd                *
c     *                                                              *
c     *     low level support function to find solid element         *
c     *     connected to top and bottom surface of cohesive element. *
c     *                                                              *
c     ****************************************************************
c

      logical function solid_cohes_init_b( nface, cohes_list, nsolid,
     &                                     solid_list )
      implicit none
      integer nface, cohes_list(*), nsolid, solid_list(*)
c
c          locals
c
      integer i, cohes_node, j
      logical found(20) ! plenty big enough
c
c          nface:   num nodes on surface of cohes element
c          cohes_list: (structure) nodes on surface
c          nsolid: number of nodes on solid element to be tested
c          solid_list: (structure) nodes on solid
c
c          return .true. if all surface nodes are present in list
c          of solid element nodes.
c
c          algorithm:
c            o init vector. one entry for each node on surface
c              of cohesive element
c            o loop over all nodes on surface of cohesive element
c            o is that node present in the node list for the
c              element passed in
c            o if each surface node of the cohesive element is
c              present in node list for solid - we have a match.
c
      found(1:nface) = .false.
      do i = 1, nface
       cohes_node = cohes_list(i)
       do j = 1, nsolid
         if( cohes_node .eq. solid_list(j) ) then
           found(i) = .true.
           exit   ! from inside loop only
          end if
       end do
      end do
c
      solid_cohes_init_b = .false.
      do i = 1, nface
       if( .not. found(i) ) return
      end do
      solid_cohes_init_b = .true.
c
      return
      end

c
c     ****************************************************************
c     *                                                              *
c     *               subroutine read_alloc_cry                      *
c     *                                                              *
c     *                    written by : mcm                          *
c     *                                                              *
c     *                last modified : 12/23/2105 rhd                *
c     *                                                              *
c     *      Read the crystal data from file while allocating the    *
c     *      appropriate structures.                                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine read_alloc_cry( fileno )
      use global_data ! old common.main
      use crystal_data, only: c_array,data_offset,angle_input,
     &      crystal_input, srequired, nangles, simple_angles,
     &      mc_array, defined_crystal
      implicit integer (a-z)
      integer, intent(in) :: fileno
c
      integer :: nelem
c
      read(fileno) defined_crystal
c
      if( defined_crystal ) then
        read(fileno) cry_multiplier
        read(fileno) c_array
        read(fileno) nelem
        read(fileno) mxcry
        if( .not. allocated(data_offset) ) then
            allocate( data_offset(noelem) )
        end if
        if( .not. allocated(angle_input) ) then
            allocate( angle_input(nelem,mxcry,3) )
        end if
        if( .not. allocated(crystal_input) ) then
            allocate( crystal_input(nelem,mxcry) )
        end if
        read(fileno) data_offset
        read(fileno) angle_input
        read(fileno) crystal_input
      end if
c
      read(fileno) srequired
c
      if( srequired ) then
        read(fileno) nangles
        if( .not. allocated(simple_angles) ) then
          allocate( simple_angles(nangles,3) )
        end if
        read(fileno) simple_angles
        if( .not. allocated(mc_array) ) then
          allocate( mc_array(noelem,max_crystals) )
        end if
        read(fileno) mc_array
      end if
c
      return
      end subroutine


