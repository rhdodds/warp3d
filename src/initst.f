c     ****************************************************************
c     *                                                              *
c     *                      subroutine initst                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 11/26/2018 rhd             *
c     *                                                              *
c     *     at program startup, initializes various variables and    *
c     *     arrays needed to set up the program correctly.           *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine initst( sbflg1, sbflg2 )
      use global_data ! old common.main
c
      use ifport
      use elem_load_data, only : elem_loads
      use segmental_curves
      use main_data, only : output_packets, packet_file_name,
     &                      packet_file_no, ascii_packet_file_name,
     &                      ascii_packet_file_no, temperatures_ref,
     &                      cp_matls_present,
     &                      cohesive_ele_types, linear_displ_ele_types,
     &                      adjust_constants_ele_types,
     &                      axisymm_ele_types, umat_used,
     &                      implemented_ele_types, bar_types,
     &                      tables, user_lists, nonlocal_analysis,
     &                      modified_mpcs, umat_serial,
     &                      convergence_history, link_types,
     &                      run_user_solution_routine, cp_unloading,
     &                      divergence_check, diverge_check_strict,
     &                      asymmetric_assembly, output_command_file,
     &                      material_model_names, batch_mess_fname,
     &                      creep_model_used, extrapolate,
     &                      extrap_off_next_step, line_search,
     &                      ls_details, ls_min_step_length,
     &                      ls_max_step_length, ls_rho, ls_slack_tol,
     &                      initial_stresses_user_routine,
     &                      initial_stresses_file, max_step_limit,
     &                      user_cnstrn_stp_factors, stpchk,
     &                      actual_cnstrn_stp_factors,
     &                      fgm_node_values_defined,
     &                      fgm_node_values_used,
     &                      fgm_node_values_cols,
     &                      initial_state_option, initial_state_step,
     &                      initial_stresses_input
c
      use stiffness_data
      use file_info
      use contact
      use damage_data
      use hypre_parameters
      use performance_data
      use distributed_stiffness_data, only: parallel_assembly_allowed,
     &                                      initial_map_type,
     &                                      final_map_type
c
      use crystal_data, only: c_array, cry_multiplier, defined_crystal,
     &                        nangles, srequired
c
      use erflgs
c
      implicit none
c
      integer :: nblank, reclen, endchr, lword, rword, i, strtop,
     &           strstp, rottop, env_length, k
      integer, external :: omp_get_max_threads
      logical :: promsw, echosw, comsw, atrdsw, eolsw, eofsw, menusw,
     &           ptsw, signsw, sbflg1, sbflg2
      double precision, parameter :: zero = 0.0d0, one=1.0d0,
     &            fourth=0.25d0, ten_billion=1.0e10,
     &            hundredth=0.01d0, twentyth=0.05d0,
     &            thousandth=0.001d0
      character(len=80) :: env_val
c
c                       initialize the file input and output parameters
c
      inlun(1) = 5
      inlun(2) = 80
      inlun(3) = 81
      inlun(4) = 82
      inlun(5) = 83
      inlun(6) = 84
      inlun(7) = 85
      inlun(8) = 86
      inlun(9) = 87
      inlun(10) = 88
c
      outlun(1) = 6
      outlun(2) = 89
c
      filcnt   = 1
      in       = inlun(filcnt)
      out      = outlun(1)
      outing   = .false.
c
      output_packets = .false.
      packet_file_no = 97
      ascii_packet_file_no = 96
      packet_file_name(1:) = ' '
      ascii_packet_file_name(1:) = ' '
      batch_mess_fname(1:) = ' '
c

c                       summary of fortran file numbers used in warp3d
c                       add new ones here....
c
c                       code updates are gradually using the function
c                       warp3d_get_device_number() to obtain an available
c                       file number for use during well-defined block
c                       of execution
c
c         terminal input:               5  (Unix standard input)
c         terminal output:              6  (Unix standard output)
c         ascii energy file:            11
c         other ascii input files
c           for *input from file:       80-88
c         other ascii output files
c           for *output to file:        89
c         ascii data packet file:       96
c         binary data packet file:      97
c         binary results file:          98
c         formatted results file:       99

c
c
c                       initialize scan
c
      call setin(in)
      call setout(out)
      nblank= 80
      reclen= 80
      endchr= 1h$
      promsw= .false.
      echosw= .true.
      comsw= .false.
      atrdsw= .false.
      eolsw= .true.
      eofsw= .true.
      menusw= .false.
      ptsw= .false.
      signsw= .false.
      call check_to_prompt( promsw )
      call scinit(nblank,reclen,endchr,promsw,echosw,comsw,atrdsw,
     &            eolsw,eofsw,menusw,ptsw,signsw)
c
c
c                       set the element library.
c
c
      call setelb( nlibel, elelib, outmap, mxlbel, out)
c
c                       initialize material library vector and its
c                       occupied and open linked lists.
c
c                       initialize strings with names of WARP3D
c                       material models to use, for example in output
c                       messages.
c
c                       cp_matls_present flag. -1 to init. = 0 if
c                       checked and no cp matls in model. = 1 if
c                       checked and cp matls present in model
c
      lword= 32460
      rword= 1
      mathed= lword*two16+rword
c
      do i = 1, mxmat
         matnam(i) = ' '
         lword     = 0
         rword     = i+1
         if( i .eq. mxmat ) rword = 32460
         matlst(i) = lword*two16+rword
      end do
      allocate( material_model_names(mxmat) )
      material_model_names(1:mxmat)(1:) = "not_used"
      material_model_names(1)(1:) = "bilinear"
      material_model_names(2)(1:) = "deformation"
      material_model_names(3)(1:) = "mises_gurson"
      material_model_names(4)(1:) = "cohesive"
      material_model_names(5)(1:) = "cyclic"
      material_model_names(6)(1:) = "creep"
      material_model_names(7)(1:) = "mises_hydrogen"
      material_model_names(8)(1:) = "umat"
      material_model_names(9)(1:) = "not_used"
      material_model_names(10)(1:) = "crystal_plasticity"
      material_model_names(11)(1:) = "interface_damage"
c
      cp_matls_present = -1
c
c
c                       initialize table library object and its
c                       variables.
c
       if ( .not. allocated(tables ) ) then
          allocate( tables(max_tables) )
          do i = 1, max_tables
             tables(i)%table_name = ' '
             tables(i)%table_type = ' '
             tables(i)%num_rows = 0
             tables(i)%num_cols = 0
          end do
       end if
c
c                       initialize library of user-defined integer lists.
c
      do i = 1, max_user_lists ! if changed, fix size in mod_main
        user_lists(i)%name = ' '
        user_lists(i)%length_list = 0
      end do
c
c                       initialize loading library vector  and its
c                       occupied and open linked lists.  Also initialize
c                       data structure for element loads

c
      lword= 32460
      rword= 1
      lodhed= lword*two16+rword
c
      do i = 1, mxlc
         lodnam(i) = ' '
         lodtyp(i) = ' '
         lword     = 0
         rword     = i+1
         if(i .eq. mxlc ) rword = 32460
         lodlst(i) = lword*two16+rword
      end do
c
      if ( .not. allocated(elem_loads) ) then
         allocate (elem_loads(mxlc))
         do i = 1, mxlc
            elem_loads(i)%size = 0
         end do
      endif
c
      new_constraints = .false.
      new_analysis_param = .false.
c
c                       initialize the existence of non-zero, reference
c                       (nodal) tempertures.
c
      temperatures_ref = .false.
      temperatures     = .false.
c
c                       initialize the constraint multipliers for
c                       each load step
c
      max_step_limit = 50000
      k = max_step_limit
      allocate( actual_cnstrn_stp_factors(k),
     &          user_cnstrn_stp_factors(k), stpchk(k) )
      do i = 1, k
         user_cnstrn_stp_factors(i) = 1.0
         stpchk(i) = .false.
      end do
c
c                       initialize pointer and counting variables.
c
      nummat = 0
      numlod = 0
      numcol = 0
      crdtop = 0
      inctop = 0
      strtop = 0
      strstp = 0
      rottop = 0
c
c                       initialize logical flags for structural
c                       size input
c
      numnod  = .false.
      numel   = .false.
      coor    = .false.
      elprop  = .false.
      elinc   = .false.
      constr  = .false.
      block   = .false.
c
c                       set flag indicating proper initialization
c                       to false.
c
      sbflg1 = .false.
      sbflg2 = .false.
      fatal  = .false.
      input_ok = .true.
c
c                       initialize structure name
c
      stname = ' '
c
c
c                       initialize time step counter and the load
c                       name for that step.
c
c
      ltmstp = 0
      lsldnm = ' '
c
c                       initialize bit vector
c
      bits(1) = 1
      do i = 2, 31
         bits(i) = 2*bits(i-1)
      end do
c
c                       initialaize dynamic/nonlinear analysis
c                       parameters to their default values.
c
      trace(1)       =  .true.
      trace(2)       =  .false.
      trace(3)       =  .false.
      trace(4)       =  .false.
      trace(5)       =  .false.
      prnres         =  .false.
      prlres         =  .false.
      growth_k_flag  = .false.
      adaptive_flag  = .false.
      qbar_flag      = .true.
      batch_messages = .false.
      mxiter         =  10
      mniter         =  2
      mxlitr         =  10
      nbeta          =  fourth
      dt             =  one
      total_model_time  =  zero
      linmas         =  .false.
      newstf         =  .false.
      newmas         =  .false.
      newtrn         =  .false.
      halt           =  .true.
      extrapolate    =  .true.
      extrap_off_next_step  = .false.
      line_search    =  .false.
      ls_details     =  .false.
      ls_min_step_length = 0.01d00
      ls_max_step_length = one
      ls_rho         = 0.7d00
      ls_slack_tol   = 0.5d00
      incflg         =  .false.
      ifvcmp         =  .false.
      convrg(1:mxcvtests)   =  .false.
      tol(1:mxcvtests)      =  zero
      emax           =  ten_billion
      fmax           =  ten_billion
      signal_flag    =  .false.
      time_limit     =  0.0
      show_details   = .true.
      sparse_stiff_output = .false.
      sparse_stiff_binary = .true.
      sparse_stiff_file_name(1:) = 'sparse_stiffness_output'
      divergence_check = .true.
      diverge_check_strict = .false.
c
c                       initialize solver parameters
c
c    **********************************************************************
c    *                                                                    *
c    *                     Selection of Solver Type                       *
c    *                (some are now obsolete and removed)                 *
c    *                                                                    *
c    *    The solver flag is an integer:                                  *
c    *      -1 = Lin. Preconditioned Conj. Gradient solver (lnpcg)        *
c    *           no longer allowed via input. code gradually removed      *
c    *      0, 1, 2, 3, 4, 5,  available                                  *
c    *      6 = (not implemented) IBM WSMP                                *
c    *      7 = Intel MKL Pardiso symmetric (Win, Mac, Linux)             *
c    *      8 = Intel MKL Pardiso asymmetric (Win, Mac, Linux)            *
c    *      9 = hypre iterative solver w/ various preconditioners (LLNL)  *
c    *          (Linux only)                                              *
c    *     10 = Cluster Pardiso - symmetric. Linux only                   *
c    *          (Linux only)                                              *
c    *     11 = Cluster Pardiso - asymmetric. Linux only                  *
c    *                                                                    *
c    **********************************************************************
c
c
      old_solver_flag = -2
      solver_flag = 7
      if( use_mpi ) solver_flag = 10
      solver_out_of_core = .false.
      solver_scr_dir(1:) = './warp3d_ooc_solver'
      solver_memory      = 500
      solver_mkl_iterative = .false.
c
c     Never used - MM
c      solver_threads = slv_threads
c
c                       initialize input error flags
c
      num_warn  = 0
      num_error = 0
      num_fatal = 0
c
c                       initialize timing parameters.
c
      times(1:mxtim,1) = zero
      times(1:mxtim,2) = zero
c
c                       initialize domain integral
c
      call initdm
c
c                       initialize the element damage/crack growth
c                       parameters
c
      no_killed_elems       = .true.
      crack_growth_type     = 0
      growth_by_kill        = .false.
      growth_by_release     = .false.
      porosity_limit    = .20d0
      smcs_alpha            = 1.0
      smcs_beta             = 1.5
      max_dam_state         = 5
      num_kill_elem         = 0
      print_status          = .false.
      kill_order            = .false.
      release_type          = 1
      crk_pln_normal_idx    = 0
      gurson_cell_size      = zero
      release_fraction      = 0.1
      critical_angle        = zero
      init_crit_ang         = zero
      num_crack_plane_nodes = 0
      no_released_nodes     = .true.
      char_length           = zero
      list_crkpln_nodes     = .false.
      list_crkfrnt_nodes    = .false.
      overshoot_control_crk_grth = .false.
      overshoot_allocated   = .false.
      control_load_fact     = one
      min_load_fact         = thousandth
      ctoa_range            = hundredth
      overshoot_limit       = twentyth
      perm_load_fact        = one
      min_steps_for_release = 6
      max_porosity_change   = zero
      max_deff_change       = 0.20
      critical_cohes_deff_fract = 5.0
      max_plast_strain_change = hundredth
      g_stp_cntrl_allocated = .false.
      load_size_control_crk_grth = .false.
      const_front           = .false.
      num_nodes_thick       = 0
      num_crack_fronts      = 0
      ctoa_dist             = zero
      init_ctoa_dist        = zero
      num_nodes_back        = 100
      master_lines_set      = .false.
      num_ctoa_released_nodes = 0
      num_steps_min         = 0
      del_poros(1:mxstp_store) = one
      del_deff(1:mxstp_store) = one
      all_elems_killed      = .false.
      load_reduced          = .false.
      killed_ele_int_work   = zero
      killed_ele_pls_work   = zero
      enforce_node_release  = .false.
      num_elements_in_force_release = 0
      ppr_kill_displ_fraction = 0.90
c
c                       initialize the segmental stress strain curve definitions
c
      seg_curves(1:max_seg_points,1:2,1:max_seg_curves) = zero
      seg_curve_table(1:max_seg_curves+1,1:max_seg_curve_sets) = 0
      num_seg_points(1:max_seg_curves) = 0
      max_current_pts    = 0
      max_current_curves = 0
      num_points         = 0
      num_seg_curve_sets = 0
c
c                       initialize the adaptive scaling factor
c
c
      scaling_adapt = one
c
c                       initialize the bbar element stabilization
c
      eps_bbar = zero
c
c                       initialize the contact surface data
c
      use_contact = .false.
      do i = 1, maxcontact
         cplane_vec(1:3,1:2,i)       = zero
         cshape_norm(1:3,i)          = zero
         cshape_pnt(1:3,i)           = zero
         cshape_rate(1:3,i)          = zero
         cshape_param(1:maxcntprm,i) = zero
         contact_stiff(i)            = zero
         contact_fric(i)             = zero
         contact_shape(i)            = 0
         contact_depth(i)            = ten_billion
         contact_outside(i)          = .true.
      end do
c
c                       initialize variables used for fgm material
c                       properties specified at the model nodes.
c
      fgm_node_values_defined = .false.
      fgm_node_values_cols    = 8
      fgm_node_values_used    = .false.
c
c                       initialize variables for support of a user
c                       defined initial state (J-integrals)
c
      initial_state_option = .false.
      initial_state_step = int(1.0e09)
c
c                       has user input initial (residual) stresses
c
      initial_stresses_input = .false.
c
c                       initialize logical flags used throughout code
c                       to indicate various characteristics of elements
c
      do i = 1, 50
        cohesive_ele_types(i)         = .false.
        linear_displ_ele_types(i)     = .false.
        adjust_constants_ele_types(i) = .false.
        axisymm_ele_types(i)          = .false.
        implemented_ele_types(i)      = .false.
        bar_types(i)                  = .false.
        link_types(i)                 = .false.
      end do
c
c                       cohesive element types
c
      cohesive_ele_types(12) = .true.
      cohesive_ele_types(14) = .true.
      cohesive_ele_types(15) = .true.
c
c                       elements considered to have linear displacement
c                       fields for temperature, fgm operations
c
      linear_displ_ele_types(2)  = .true.
      linear_displ_ele_types(5)  = .true.
      linear_displ_ele_types(12) = .true.
      linear_displ_ele_types(13) = .true.
      linear_displ_ele_types(14) = .true.
c
c                       non-hex and non-interface elements require
c                       adjustments in geometry constants (triangles,
c                       tets, wedges, axisymmetric)
c
      adjust_constants_ele_types(6)  = .true.
      adjust_constants_ele_types(7)  = .true.
      adjust_constants_ele_types(8)  = .true.
      adjust_constants_ele_types(10) = .true.
      adjust_constants_ele_types(11) = .true.
      adjust_constants_ele_types(13) = .true.
c
c                       axisymmetric elements
c
      axisymm_ele_types(10) = .true.
      axisymm_ele_types(11) = .true.
c
c                       bar_types, link_types
c
      bar_types(18)  = .true.
      link_types(19) = .true.
c
c                       element types currently implemented
c                       see setelb routine
c
      implemented_ele_types(1:6)   = .true.
      implemented_ele_types(11)    = .true.
      implemented_ele_types(12:15) = .true.
      implemented_ele_types(18:19) = .true.

c
c                       initialize for possible nonlocal_analysis using
c                       cohesive elements.
c
      nonlocal_analysis = .false.
c
c                       initialize mpc and tied contact data structures
c
      call mpc_init()  ! sets lots of flags
      modified_mpcs = .false.
c
c                       initialize stiffness matrix data structures
c                       most of these are only used if mpcs present
c
      ncoeff     = 0
      big_ncoeff = 0
      newcount   = 0
      dep_len    = 0
      ind_len    = 0
      dia_len    = 0
      temp_len   = 0
c
      if( allocated(new_locations) )   deallocate( new_locations )
      if( allocated(k_indexes) )       deallocate( k_indexes )
      if( allocated(ind_temp) )        deallocate( ind_temp )
      if( allocated(new_indexes) )     deallocate( new_indexes )
      if( allocated(new_ptrs) )        deallocate( new_ptrs )
      if( allocated(dep_locations) )   deallocate( dep_locations )
      if( allocated(ind_locations) )   deallocate( ind_locations )
      if( allocated(diag_locations) )  deallocate( diag_locations )
      if( allocated(dep_loc) )         deallocate( dep_loc )
      if( allocated(ind_loc) )         deallocate( ind_loc )
      if( allocated(dia_loc) )         deallocate( dia_loc )
      if( allocated(k_coeffs) )        deallocate( k_coeffs )
      if( allocated(cof_temp) )        deallocate( cof_temp )
      if( allocated(total_lagrange_forces) )
     &     deallocate( total_lagrange_forces )
      if( allocated(d_lagrange_forces) ) deallocate(d_lagrange_forces)
      if( allocated(i_lagrange_forces) ) deallocate(i_lagrange_forces)
c
c
c                        initialize for OMP threaded execution
c                        omp_get_max_threads = omp_num_threads
c
      num_threads = max( 1, omp_get_max_threads() )
      if( num_threads .gt. max_threads ) then
            write(out,9220) max_threads
            call die_abort
      end if
c
c                       global flags for UMAT used in model,
c                       run UMAT element blocks in serial
c
      umat_serial = .false.
      umat_used   = .false.
c
c                       global flags for modeling containing
c                       a material that creep. will cause
c                       iter =0 computations to be run
c
      creep_model_used = .false.
c
c                       hypre solver defaults.  See the hypre reference manual
c                       for rational (if provided) as I just copy
c                       these over from there.  The exception is the
c                       maximum iterations value (hypre_max)
c                       which I increased to 10,000 to allow for our typically
c                       ill-conditioned plasticity problems.
c
      precond_type = 1
      hsolver_type = 1
      precond_printlevel = 0
      solver_printlevel = 0
c
      hypre_tol = 1.0D-8
      hypre_max = 10000
c
      levels = 0
      threshold = 0.25
      filter = 0.1
      symme = 1
      loadbal = 0.0
c
      hyp_first_solve = .true.
c
      precond_fail_count = 0
c
      hyp_trigger_step = .false.
c
      max_levels = 10
      mg_threshold = 0.8
      coarsening = 10
      agg_levels = 1
      interpolation = 5
      truncation = 0.0
      relaxation = 6
      relax_wt = 1.0
      relax_outer_wt = 1.0
      cf = 0
      cycle_type = 1
      sweeps = 1
c
c           Performance data defaults (zero out assembly counter)
c
      ntimes_assembly = 0
      assembly_total = 0.0
      time_assembly = .false.
c
c           Distributed assembly params
c
      parallel_assembly_allowed = .true.
      initial_map_type = 2
      final_map_type = 1
c
c           Crystal properties
c
      do i = 1, max_crystals
            c_array(i)%valid = .false.
      end do
      cry_multiplier = 1.0
      defined_crystal = .false.
      cp_unloading = .false.
      nangles = 0
      srequired = .false.
c
c                       convergence history data for user
c                       routine to adjust next step loading,
c                       solution parameters, etc.
c
      do i = 1, 5
       convergence_history(i)%step_converged = .false.
       convergence_history(i)%adaptive_used = .false.
       convergence_history(i)%iterations_for_convergence = 0
       convergence_history(i)%adapt_substeps = 0
      end do
      run_user_solution_routine = .false.
c
c                       Asymmetric assembly
c
      asymmetric_assembly = .false.
c
c                       initialize the file name for the
c                       "output commands file ... steps ..."
c
      output_command_file(1:) = " "
c
c                       message counters to shut of gazillions of
c                       the same messages for large models
c
      msg_count_1 = 0
      msg_count_2 = 0
c
c                       model initial stresses
c
      initial_stresses_file = " "
      initial_stresses_user_routine = .false.
c
c                       blocking table initial size. now dynamically
c                       allocated
c
      mxnmbl = 0
      ptr_iprops = 0
c
      return
c
 9220 format(/1x,'>>>>> fatal error: the number of threads',
     &          ' exceeds the',/,
     &          '              WARP3D limit of:',i5,/,
     &          '              Job aborted at this point...',/)
c
      end
c
c     ****************************************************************
c     *                                                              *
c     *               subroutine set_element_type                    *
c     *               (requires element type no.)                    *
c     *                                                              *
c     *                       written by: rhd                        *
c     *                   last modified : 08/11/2017 rhd             *
c     *                                                              *
c     ****************************************************************
c
      subroutine set_element_type( element_type, threed_solid_elem,
     &                             hex_elem, wedge_elem, tet_elem,
     &                             twod_elem, quad_elem,
     &                             triangle_elem, axisymm_elem,
     &                             cohesive_elem, bar_elem, link_elem )
      implicit none
c
      integer :: element_type
      logical :: threed_solid_elem, hex_elem, wedge_elem, tet_elem,
     &           twod_elem, quad_elem, triangle_elem,
     &           axisymm_elem, cohesive_elem, bar_elem, link_elem
c
      hex_elem      = element_type .ge. 1
     &                .and. element_type .le. 5
c
      wedge_elem    = element_type .eq. 7
c
      tet_elem      = element_type .eq. 6
     &                .or. element_type .eq. 13
c
      quad_elem     = element_type .eq. 9
     &                .or. element_type .eq. 10
c
      triangle_elem = element_type .eq. 8
     &                .or. element_type .eq. 11
c
      axisymm_elem  = element_type .eq. 10
     &                .or. element_type .eq. 11
c
      cohesive_elem = element_type .eq. 12
     &                .or. element_type .eq. 14
     &                .or. element_type .eq. 15
c
      threed_solid_elem = hex_elem .or. wedge_elem .or. tet_elem
c
      twod_elem = quad_elem .or. triangle_elem .or. axisymm_elem
c
      bar_elem  = element_type .eq. 18
c
      link_elem = element_type .eq. 19
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *            subroutine adjust_scalar_weights                  *
c     *              (requires element type no.)                     *
c     *                       written by: rhd                        *
c     *                   last modified : 02/1/02                    *
c     *                                                              *
c     ****************************************************************
c
      subroutine adjust_scalar_weights( elem_type, scalar )
      implicit none
c
c                  argument variables
c
c
      integer :: elem_type
      double precision :: scalar
c
c                   local variables
c
      logical :: tet_elem, tri_elem, axisymm_elem
c
      double precision :: two, half, pi
      data two, half, pi / 2.0, 0.5, 3.14159265358979323846 /
c
      tet_elem = .false.
      tri_elem = .false.
      axisymm_elem = .false.
      scalar = 1.0d00
c
c
c                    the element type determines the scalar
c                    multiple to get the correct adjustment to |J|.
c
c                    planar triangle, multiply by 0.5.
c                    tetrahedral element, multiply by 1/6.
c
c
      tet_elem = elem_type .eq. 6 .or. elem_type .eq. 13
      tri_elem = elem_type .eq. 8
      axisymm_elem = elem_type .eq. 10 .or. elem_type .eq. 11
c
c                    tetrahedral elements
c
      if( tet_elem ) then
      scalar = 1.0D0/6.0D0
        return
      end if
c
c                     triangular elements
c
      if( tri_elem ) then
        scalar = half
        return
      end if
c
c                      axisymmetric elements
c
c                      axisymmetric elements, multiply
c                      by 2*pi*radius
c
c                      note that for the axisymmetric
c                      triangle 2*pi*0.5 = pi.
c
      if( axisymm_elem ) then
         if( elem_type .eq. 10 ) then
           scalar = two*pi
         else if( elem_type .eq. 11 ) then
           scalar = pi
         end if
         return
      end if
c
      end

c     ****************************************************************
c     *                                                              *
c     *              function  warp3d_matl_num                       *
c     *                                                              *
c     *       return the internal  number for a material model       *
c     *               written by: rhd                                *
c     *                   last modified : 5/28/2016 rhd              *
c     *                                                              *
c     ****************************************************************
c
      integer function warp3d_matl_num( material_model_id, nc )
c
      use main_data, only : material_model_names
      implicit none
      include 'param_def'
c
      integer :: nc, m
      character(len=nc) :: material_model_id
c
      select case( material_model_id )

      case( "bilinear" )
       m = 1
      case( "deformation", "deform" )
       m = 2
      case( "mises", "gurson", "mises-gurson", "mises_gurson" )
       m = 3
      case( "cohesive" )
       m = 4
      case( "cyclic", "adv cyclic", "adv. cyclic" )
       m = 5
      case( "creep", "Norton", "norton" )
       m = 6
      case( "mises_hydrogen", "hydrogen" )
       m = 7
      case( "umat", "UMAT", "um" )
       m = 8
      case( "crystal", "CP", "cp", "crystal_plasticity" )
       m = 10
      case default
         write(*,9000) material_model_id
        call die_abort
      end select
c
      warp3d_matl_num = m
      return
 9000 format('>> FATAL ERROR: routine warp3d_matl_num ',
     & /,    '                invalide material model name: ',a,
     & /,    '                job aborted',//)
      end
