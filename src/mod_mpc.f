c
c     ****************************************************************
c     *                                                              *
c     *                        module mod_mpc                        *
c     *                                                              *
c     *                       written by : gvt                       *
c     *                                                              *
c     *                    last modified : 09/3/2017 rhd             *
c     *                                                              *
c     *  Module for the multipoint constraint (MPC) and tied contact *
c     *  data structures; use a module to permit use of global       *
c     *  values in the MPC calculations.                             *
c     *                                                              *
c     ****************************************************************
c
c
      module mod_mpc
c
c           Declare Modules
c
      include 'param_def'
c
c           Declare Variables
c
c
c           Variables for the multipoint constraint equations
c
c           num_user_mpc = number of multipoint constraint (MPC)
c               equations specified by the user; stored in the
c               user_mpc_table() array
c           num_tied_con_mpc = number of MPC equations generated
c               by the tied contact processor; stored in the
c               tied_con_mpc_table() array
c
c           Components of the "mpc_eqn" data type:
c
c           num_terms = number of terms in the MPC equation
c           constant = right hand side constant value for the MPC
c               equation; in general the constant can be zero or
c               non-zero, but expecting a constant equal to zero for
c               tied contact
c           node_list = integer pointer component array to store the
c               node ID numbers in the MPC equation; set the dimension
c               of node_list to num_terms
c           dof_list = integer pointer component array to store the
c               node DOF numbers in the MPC equation (1=u, 2=v, 3=w);
c               set the dimension of node_list to num_terms
c           multiplier_list = real pointer component array to store
c               the coefficient multiplier for each node DOF in the
c               MPC equation; set the dimension of node_list to
c               num_terms
c
c           Arrays declared with the "mpc_eqn" data type:
c
c           user_mpc_table = store the MPC equations specified by the
c               user; the MPC equation values are stored in the data
c               type components
c           tied_con_mpc_table = the MPC equations generated and stored
c               by the tied contact processing routine; expect the right
c               hand side constant to be zero for the tied contact MPC
c               equations; the tied contact MPC equation values are
c               stored in the data type components
c
c           Note that the max_mpc array size for user_mpc_table
c           and max_mpc_tied for tied_con_mpc_table are set in
c           mem_allocaste.f
c
      integer  num_user_mpc, num_tied_con_mpc, nmpc
c
      logical  mpcs_exist
c
      type :: mpc_eqn
        integer num_terms
        real constant
        integer, pointer,  dimension (:) :: node_list
        integer, pointer, dimension (:) :: dof_list
        real, pointer, dimension (:) :: multiplier_list
      end type
c
      type (mpc_eqn), allocatable,
     &            dimension (:) :: user_mpc_table
      type (mpc_eqn), allocatable,
     &            dimension (:) :: tied_con_mpc_table
c
c
c           Variables for the tied contact data
c
c           num_surfaces = number of mesh surfaces defined in the input
c               data file; these are all of the master and slave mesh
c               surfaces used in the tied contact sets
c           num_tied_sets = number of tied contact pair sets defined in
c               the input data file; each tied contact set can use one
c               or more pairs of surfaces
c
c           Components of the "surface" data type:
c
c           id = character label for the mesh surface name
c           num_elems = number of elements in the list to define the
c               surface; use to set the size of elem_list and face_list
c           elem_list = integer pointer component array to store the
c               element ID numbers that define the surface
c           face_list = integer pointer component array to store the
c               element face numbers that define the surface
c
c           Components of the "tied_set" data type:
c
c           id = character label for the tied contact set name
c           tolerance = gap tolerance distance between the master and
c               slave surface for slave nodes to be in contact with the
c               master surface
c           adjust_gap = logical flag, if true then adjust the slave
c               node initial position to be located on the master
c               surface
c           num_pairs = number of surface pairs in the tied contact set;
c               a tied contact set is defined with one or more pairs
c               of surfaces (a master and slave surface)
c           master_list = integer pointer component array to store the
c               the master surface ID numbers in the tied contact set;
c               use this integer surface ID number to access the
c               data in the surface_table() array
c           slave_list = integer pointer component array to store the
c               the slave surface ID numbers in the tied contact set;
c               use this integer surface ID number to access the
c               data in the surface_table() array
c
c           Array declared with the "surface" data type:
c
c           surface_table = store the data for all the mesh surfaces
c               defined by the user; the element list for each surface
c               is stored in the data type components
c
c           Array declared with the "tied_set" data type:
c
c           tied_contact_table = store the data for all the tied contact
c               sets defined by the user; the surface lists and data are
c               stored in the data type components
c
c           Note that the max_surfaces and max_tied_sets array sizes
c           are set in param_def
c
c
      integer  num_surfaces, num_tied_sets
c
      logical  tied_sets_exist, tied_con_mpcs_constructed,
     &         display_mpcs
c
      type :: surface
        character(len=16) :: id
        integer num_elems
        integer, pointer, dimension (:) :: elem_list
        integer, pointer, dimension (:) :: face_list
      end type
c
      type :: tied_set
        character(len=16) :: id
        real tolerance
        logical adjust_gap
        integer num_pairs
        integer, pointer, dimension (:) :: master_list
        integer, pointer, dimension (:) :: slave_list
      end type
c
c
      type (surface), allocatable, dimension (:) :: surface_table
      type (tied_set), allocatable, dimension (:) :: tied_contact_table
c
c
c           Components of the "nonzero_locs" data type:
c
c           length = the current length of the integer vector component
c           loc_list = integer vector containing the list of non-zero
c               indices in the equation (similar to k_indexes)
c
c           Arrays declared with the "nonzero_locs" data type:
c
c           eqn_row = used by the MPC solver routines to temporarily store
c               the row portions of all equations in the stiffness matrix
c               while adding terms (later converted k_indexes)
c               It is also used to store the entire dep eqns (row & col).
c               This structure is deallocated after the first MPC solve.
c
      type :: nonzero_locs
        integer length
        integer, pointer, dimension (:) :: loc_list
      end type
c
      type (nonzero_locs), allocatable, dimension(:) :: eqn_row
c
c
c           Variables for MPC solver routines
c
c           dep_check = integer vector used to tell if a given equation
c               number is a dependent term in an MPC equation and which
c               MPC equation it is in, e.g. if dof 13 is a dep term in
c               in MPC #5, dep_check(13) = 5, else = 0
c           dep_dof = lists the dep term of each MPC eqn, e.g. if MPC eqn
c               #5 has dep term 13, dep_dof(5) = 13
c           ind_dof = lists all the independent terms of each MPC eqn
c           num_terms = contains the number of ind terms in each MPC eqn,
c               e.g. if MPC eqn #5 has 6 terms, num_terms(5) = 5 (because
c               the first term is the dep).  This vector is used to reference
c               the correct positions in 'ind_dof' and 'multi_list' for
c               each MPC eqn
c           multi_list = containst the multiplier corresponding to each ind
c               term in each MPC, therefore has same length as ind_dof
c
c           The 'dep_dof', 'ind_dof', 'num_terms', and 'multi_list', data
c               structures are used to prevent excessive access to the f90
c               data structures used to contain the MPC equations.  They are
c               quick to create, but can obviously take a lot of memory, so
c               they are created and then deallocated for each MPC solve.
c
      integer  dep_trms_len
c
      integer, allocatable,
     &        dimension (:) :: dep_check, dep_dof,
     &                         ind_dof, num_terms,
     &                         dep_ptr, num_dep_trms,
     &                         dep_trms, abs_dep_ptr
c
      real, allocatable, dimension (:) :: multi_list
c
      double precision,
     &            allocatable,
     &            dimension (:) :: dep_coef, dep_rhs
c
c
      end module mod_mpc
