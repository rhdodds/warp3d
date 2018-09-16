c     ****************************************************************
c     *                                                              *
c     *                    f-90 module erflgs                        *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *              last modified : 6/5/2017 rhd                    *
c     *                                                              *
c     *     this small module replaces the old common /erflgs/       *
c     *                                                              *
c     ****************************************************************
c
      module erflgs

      logical :: numnod, numel, fatal, coor, elprop, elinc,
     &           constr, block
c
      end module erflgs
c     ****************************************************************
c     *                                                              *
c     *                    f-90 module erflgs                        *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *              last modified : 2/8/2018 rhd                    *
c     *                                                              *
c     *     this small module replaces the old common.main           *
c     *                                                              *
c     ****************************************************************


      module global_data
      implicit none
c
      include 'param_def'
c
c              --- double precision arrays/vectors ---
c
      double precision :: tol(mxcvtests)
!dir$ attributes align: 64 :: u, v, a, du,idu, load,res, ifv, c, tol
      double precision, allocatable, dimension(:) :: u, v, a, du, idu,
     &                      load, res, ifv, c
c
c              --- double precision scalars ---
c
      double precision ::
     &    dt, nbeta, stplen, aparm, bparm, cparm, eparm,
     &    fparm, emax, fmax, prdmlt, total_mass, ext_work, beta_fact,
     &    total_model_time, eps_bbar, sum_ifv, sum_loads,
     &    internal_energy, plastic_work, killed_ele_pls_work,
     &    killed_ele_int_work, scaling_adapt, scaling_factor
c
c              --- real arrays/vectors/scalars ---
c
!dir$ attributes align: 64:: props, times
      real :: props(mxelpr,mxel)
      real :: times(mxtim,2)
      real :: strtm, time_limit
c
c              --- logical arrays/vectors/scalars ---
c
!dir$ attributes align: 64::trace, convrg
      logical ::
     &    prnres, halt, linmas, newstf, zrocon, newtrn,
     &    newmas, incflg, ifvcmp, prlres, input_ok, adaptive_flag,
     &    new_constraints, batch_messages, signal_flag,
     &    scalar_blocking, growth_k_flag, qbar_flag,
     &    solver_out_of_core, show_details, new_analysis_param,
     &    sparse_stiff_output, sparse_stiff_binary, solver_threads,
     &    sparse_research, solver_mkl_iterative, temperatures,
     &    root_processor, slave_processor, worker_processor,
     &    use_mpi, ! end of scalars
     &    trace(ntrc), convrg(10)
c
!dir$ attributes align: 64 :: lprops
      logical :: lprops(mxelpr,mxel) ! only for equivalencing to props
      equivalence ( lprops, props )
c
c              --- integer arrays/vectors/scalars ---
c
!dir$ attributes align: 64 :: dstmap,cstmap
      integer, allocatable, dimension (:) :: dstmap, cstmap

!dir$ attributes align: 64 :: gpmap, cp, dcp, icp,
     &    matlst, lodlst, prslst, plrlst, stprng, state, bits, outmap,
     &    elblks, blk_ptr_head, MPI_DOF_LOCAL, num_dof_local,
     &    proc_pids
      integer :: gpmap(mxtgp),
     &    cp(mxedof), dcp(mxedof), icp(mxutsz,2), matlst(mxmat),
     &    lodlst(mxlc), prslst(mxlsz), plrlst(mxlsz), stprng(mxlc,2),
     &    state(mxtgp), bits(31),  outmap(mxlbel,mxelmp),
     &    elblks(0:3,mxnmbl), blk_ptr_head(0:max_procs - 1),
     &    MPI_DOF_LOCAL(0:max_procs-1), num_dof_local(0:max_procs-1),
     &    proc_pids(1:max_procs-1)     ! end of arrays
c
      integer :: noelem, nonode, nummat, nogp, numcol,
     &    nodof, nlibel, numlod, nprs, nplrs, numstc, nelblk, numgrp,
     &    lgnmcn, mxiter, mniter, lgoump, mxlitr, num_term_ifv,
     &    num_term_loads, mathed, csthed,  lodhed, inctop, crdtop,
     &    in, out, histep, lowstp, ltmstp, num_warn, num_error,
     &    num_fatal, solver_flag, old_solver_flag, solver_memory,
     &    num_threads,  max_mpc, max_mpc_tied,
     &    myid, numprocs, MPI_VAL, douextdb
c
c                  counters to shut off error messages at some point
c
      integer, save :: msg_count_1, msg_count_2   ! needs atomic update
c
!dir$ attributes align:64 :: iprops
      integer :: iprops(mxelpr,mxel) ! only for equiv with props
      equivalence ( iprops, props )
c
c              --- character arrays/vectors/scalars ---
c
      character(len=8) :: lodnam(mxlc), lodtyp(mxlc), elelib(mxlbel),
     &                    snames(mxstc), stname, lsldnm
      character(len=24) :: matnam(mxmat)
      character(len=80) :: solver_scr_dir, sparse_stiff_file_name
c
c              --- end of common_main file ---
c
      end module global_data
c     ****************************************************************
c     *                                                              *
c     *                     module main_data                         *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *              last modified : 6/28/2018 rhd                   *
c     *                                                              *
c     *     define the data structures for main, large arrays        *
c     *     used in warp3d solutions. also other variables as we     *
c     *     gradually reduce dependence on common.main               *
c     *                                                              *
c     ****************************************************************
c
      module main_data

      logical :: windows_os, linux_os, osx_os
c
c
c                element incidencs and incmap.
c
c
      integer, dimension(:), save, allocatable ::  incmap, incid
c
c
c                inverse of element incidences, i.e., the list
c                of elements connected to each model node. also
c                the inverse dof maps.
c
c
      type :: elements_on_node
        integer element_count
        integer, allocatable, dimension(:) :: element_list
      end type
      type(elements_on_node), save, allocatable,
     &                            dimension(:) :: inverse_incidences
c
      type :: inv_dof_table
        integer, allocatable, dimension(:,:) :: edof_table
      end type
      type(inv_dof_table), save, allocatable,
     &                            dimension(:) :: inverse_dof_map
c
c                temporary storage of element definitions until
c                start of computations.
c
c
      integer, dimension(:,:), save, allocatable ::  elstor
c
c
c                nodal load definitions (flags, dofs, etc.)
c
c
      type :: node_loading_cond
        integer node_count
        integer how_defined
        character(len=80) :: user_file_name
        integer, pointer, dimension(:,:) :: nodal_loads
      end type
      type(node_loading_cond), save, allocatable,
     &                         dimension(:) :: node_load_defs
      integer, save, allocatable, dimension(:) :: temp_nodmap
      integer, save, allocatable, dimension(:,:) :: temp_nodlod
      integer num_loddat_blks, sizeof_loddat_blks, next_loddat_col,
     &        max_loddat_blks
      type :: node_load_block
        real, pointer, dimension(:,:) :: block
      end type
      type(node_load_block), save, allocatable,
     &                       dimension(:) :: loddat_blocks
c
c
c                nodal and element temperatures. current
c                totals and step increments.
c
c
      double precision,
     &       dimension(:), save, allocatable ::
     &        temper_nodes, dtemp_nodes, temper_elems,
     &        dtemp_elems, temper_nodes_ref
      logical temperatures_ref
c
c
c                nodal constraints definitions. release constraints
c
c
      double precision,
     &       dimension(:), save, allocatable :: cnstrn, cnstrn_in
c
      type :: release_con_info
        integer num_release_steps, remaining_steps_for_release
      double precision
     &   reaction_force
      end type
      type(release_con_info), save, allocatable,
     &                 dimension(:,:) :: release_cons_table
      integer :: release_cons_steps
c
c
c                nodal displacement vector for adaptive save/
c                restart
c
c
      double precision,
     &       dimension(:), save, allocatable :: du_nm1
c
c
c                global element-to-block mapping vectors
c                vectors.
c
c
      integer, save, allocatable, dimension(:,:) :: elems_to_blocks
c
c
c                global coordinate mapping
c
c
      integer, save, allocatable, dimension(:) :: crdmap
c
c
c                global mapping vectors vectors.
c
c
      integer, save, allocatable, dimension(:) ::     invdst
      logical, save, allocatable, dimension(:) ::     repeat_incid
c
c
c                global data to store non-global constraint
c                transformations
c
c
      type :: trn_ptr_type
      double precision, dimension(:,:), allocatable :: mat
      end type
c
      type (trn_ptr_type), save, allocatable, dimension(:) :: trnmat
      logical, save, allocatable, dimension(:) ::     trn
c
c
c                global data for diagonal mass, pbar
c
c
      double precision,
     &      save, allocatable, dimension(:) :: mdiag, pbar
c
c
c                global data for rloads, dloads
c
c
      double precision,
     &      save, allocatable, dimension(:) :: rload, dload, rload_nm1,
     &                                         total_user_nodal_forces
      double precision,
     &      save, allocatable, dimension(:,:) :: load_pattern_factors
c
c                global data to store the load pattern numbers and
c                multipliers for each nonlinear load step
c
      integer :: max_step_limit
      real, allocatable, dimension(:) :: actual_cnstrn_stp_factors,
     &                                   user_cnstrn_stp_factors
      logical, allocatable, dimension(:) :: stpchk
      type :: load_data_for_a_step
        integer :: num_load_patterns
        integer, allocatable, dimension(:) :: load_patt_num
        double precision,allocatable,dimension(:)::load_patt_factor
      end type
c
      type(load_data_for_a_step), save, allocatable,
     &                            dimension(:) :: step_load_data
c
c
c                storage for element equivalent force vectors for
c                current load step (applied pressures, tractions,
c                body forces).
c
c
      logical :: elem_equiv_loads_now
      type :: elem_forces
        integer :: ncols
        double precision, pointer, dimension(:,:) :: forces
      end type
c
      type(elem_forces), save, allocatable,
     &                         dimension(:) :: elem_eq_loads
      double precision,
     &      save, allocatable, dimension(:) :: eq_node_forces
      integer, save, allocatable, dimension(:) :: eq_node_force_indexes
      integer :: eq_node_force_len
c
c
c                 global variables being gradually moved from
c                 common.main
c
c
      logical :: nonlocal_analysis, modified_mpcs,
     &           divergence_check, diverge_check_strict
c
c
c                 information for output packets
c
c
      integer :: packet_file_no, ascii_packet_file_no
      logical :: output_packets
      character(len=50) :: packet_file_name
      character(len=80) :: ascii_packet_file_name
      character(len=80) :: batch_mess_fname
c
c
c                 material properties specified at nodes to support
c                 functionally graded materials
c
c
      real, save, allocatable, dimension(:,:) :: fgm_node_values
      logical :: fgm_node_values_defined, fgm_node_values_used
      integer :: fgm_node_values_cols
c
c
c                 logical vectors indicating element types with
c                 specific characteristics. initialized in initst.f
c
c
      logical :: cohesive_ele_types(50),
     &           linear_displ_ele_types(50),
     &           adjust_constants_ele_types(50),
     &           axisymm_ele_types(50),
     &           implemented_ele_types(50), bar_types(50),
     &           link_types(50)
c
c
c                 material properties array. these sizes correspond to
c                 mxmtpr x  mxmat in param_def and must always
c                 be consistent !
c
      integer :: imatprp(300,500)
      real  :: matprp(300,500)
      logical :: lmtprp(300,500)
      double precision :: dmatprp(300,500)
      equivalence (matprp,lmtprp)
      character(len=24), dimension(300,500) :: smatprp
c
c
c
c                 general tables of input for use throughout code.
c                 only "piston" input loading supported at present
c
c
      type :: table_entry
       character(len=24) :: table_name
       character(len=8) ::  table_type
       integer num_rows
       integer num_cols
       real, dimension(:,:), allocatable :: table_values_sgl
       double precision, dimension(:,:),
     &                   allocatable :: table_values_dbl
      end type
c
      type (table_entry), save, allocatable,
     &                          dimension(:) :: tables
c
c
c                 used defined, named lists of integers.
c                 usually node numbers or element numbers
c
c
      type :: ulist
        character(len=24) :: name
        integer :: length_list
        integer, allocatable, dimension(:) :: list
      end type
c
      type (ulist), dimension(100) :: user_lists ! 100 is set in param_def
c
c               UMAT model used ? Force serialization of umats?
c
      logical :: umat_serial, umat_used
c
c
c               creep material appears in solution. will force
c               iter=0 computations
c
      logical :: creep_model_used
c
c                 convergence information for last few load steps.
c                 used for user_solution_paramters
c
      type :: step_convergence_data
        logical :: step_converged
        logical :: adaptive_used
        integer :: iterations_for_convergence
        integer :: adapt_substeps
      end type
      type( step_convergence_data ), dimension(5) ::
     &   convergence_history
      logical :: run_user_solution_routine
c
c                 A CP flag, stick here b/c it's a solution parameter
c
      logical :: cp_unloading
c
c                 Another solution parameter telling us whether
c                 or not to use asymmetric assembly
c
      logical :: asymmetric_assembly
c
c          file name for "output commands file ... after steps <list>'
c          bit map to store expanded list of steps
c
      character(len=80) :: output_command_file
      integer, save, allocatable, dimension(:) ::
     &         output_step_bitmap_list
c
c          string names for WARP3D material models
c
      character(len=20), save,
     &   allocatable, dimension(:) :: material_model_names
c
c                 do we have extrapolated displacement increments
c                 for the step and/or non-zero imposed
c                 (user) displacement increments for the step
c                 global extrapolate flag and no extrapolate next
c                 step only

      logical :: extrapolated_du, non_zero_imposed_du,
     &           extrapolate, extrap_off_next_step
c
c                 line search parameters
c
      logical :: line_search, ls_details
      double precision ::
     & ls_min_step_length, ls_max_step_length, ls_rho,
     & ls_slack_tol
c
c                 does model have crystal plasticity materials.
c                 no need to savein restart file
c                 see init.f
c
      integer :: cp_matls_present
c
c                 options for states results on usual output cmd
c
      integer :: output_states_type_opt1, output_states_type_opt2
c

      double precision, allocatable :: initial_stresses(:,:)
      logical ::  initial_stresses_user_routine, 
     &            initial_stresses_input
      character(len=100) :: initial_stresses_file
c
c                 support for the initial state concept
c
      logical :: initial_state_option
      integer :: initial_state_step
c
      end module





