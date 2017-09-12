c
c     ****************************************************************
c     *                                                              *
c     *                       module mm10_defs                       *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified: 12/20/2016 tjt              *
c     *                                                              *
c     *       small module to hold crystal update data structs       *
c     *       also hold integer indexes into history vector for      *
c     *       an integration point                                   *
c     *                                                              *
c     ****************************************************************
c
      module mm10_defs
c
      implicit integer (a-z)
      include 'param_def'
c              includes all key size limits for WARP3D (e.g. mxvl)
c
c              Massive list of properties
c
      type :: crystal_props
        double precision :: rate_n, tau_hat_y, G_0_y, burgers,
     &        p_v, q_v, boltzman, theta_0, eps_dot_0_v, eps_dot_0_y,
     &        p_y, q_y, tau_a, tau_hat_v, G_0_v,
     &        k_0, mu_0, D_0, T_0, tau_y, tau_v, voche_m,
     &        u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, iD_v
        double precision :: cp_001, cp_002, cp_003, cp_004, cp_005,
     &                      cp_006, cp_007, cp_008, cp_009, cp_010,
     &                      cp_011, cp_012, cp_013, cp_014, cp_015,
     &                      cp_016, cp_017, cp_018, cp_019, cp_020,
     &                      cp_021, cp_022, cp_023, cp_024, cp_025,
     &                      cp_026, cp_027, cp_028, cp_029, cp_030,
     &                      cp_031, cp_032, cp_033, cp_034, cp_035,
     &                      cp_036, cp_037, cp_038, cp_039, cp_040,
     &                      cp_041, cp_042, cp_043, cp_044, cp_045,
     &                      cp_046, cp_047, cp_048, cp_049, cp_050,
     &                      cp_051, cp_052, cp_053, cp_054, cp_055,
     &                      cp_056, cp_057, cp_058, cp_059, cp_060,
     &                      cp_061, cp_062, cp_063, cp_064, cp_065,
     &                      cp_066, cp_067, cp_068, cp_069, cp_070,
     &                      cp_071, cp_072, cp_073, cp_074, cp_075,
     &                      cp_076, cp_077, cp_078, cp_079, cp_080,
     &                      cp_081, cp_082, cp_083, cp_084, cp_085,
     &                      cp_086, cp_087, cp_088, cp_089, cp_090,
     &                      cp_091, cp_092, cp_093, cp_094, cp_095,
     &                      cp_096, cp_097, cp_098, cp_099, cp_100
        double precision :: atol, atol1, rtol, rtol1, xtol, xtol1
        double precision, dimension(3,3) :: g
        double precision :: ms(6,max_slip_sys), qs(3,max_slip_sys),
     &                      ns(3,max_slip_sys)
        double precision, dimension(6,6) :: stiffness
        integer :: angle_type, angle_convention, nslip,
     &        h_type, miter, gpp, s_type, cnum, method,
     &        st_it(3), num_hard, tang_calc, out
        logical :: solver, strategy, debug, gpall, alter_mode
        ! constants for use in material models
        double precision, dimension(:,:), allocatable :: Gmat,Hmat
      end type
c
      type :: crystal_state
        double precision, dimension(3,3) :: R, Rp
        double precision, dimension(6) :: stress, D, eps
        double precision, dimension(3) :: euler_angles
        double precision, dimension(max_slip_sys) :: tau_l, slip_incs
        double precision, dimension(3,3,3) :: gradFeinv
        double precision, dimension(6,6) :: tangent
        double precision, dimension(max_uhard) :: tau_tilde
        double precision, dimension(max_uhard) :: tt_rate
        double precision :: temp, tinc, dg, tau_v, tau_y,
     &                      mu_harden, work_inc, p_work_inc,
     &                      p_strain_inc
        double precision :: ms(6,max_slip_sys), qs(3,max_slip_sys),
     &                      qc(3,max_slip_sys)
        double precision, dimension(max_uhard) :: u
        double precision, dimension(6) :: ed, ep
        integer :: step, elem, gp, iter
      end type
c
c              store integer indexes into history vector for a an
c              integration point.
c              WARP3D makes sure an mm10 routine is called to
c              create the arraays & set values. No need to save
c              across restarts. These arrays are red-only during
c              threaded processing of element blocks.
c
      integer, allocatable :: indexes_common(:,:),
     &                        index_crys_hist(:,:,:),
     &                        length_crys_hist(:),
     &                        length_comm_hist(:)
      integer :: num_common_indexes, num_crystal_terms,
     &           one_crystal_hist_size, common_hist_size
c
      end module
      
c
c     ****************************************************************
c     *                                                              *
c     *                       module mm10_constants                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 10/13/2016 rhd              *
c     *                                                              *
c     *       small module to hold numerical constants for           *
c     *       connivence                                             *
c     *                                                              *
c     ****************************************************************
c
      module mm10_constants
      implicit none
      
      double precision, parameter ::
     &   zero = 0.0d0, one = 1.d0, two = 2.0, three = 3.0d0, 
     &   four = 4.0d0, ten = 10.0d0, one_eighty = 180.0d0
c
      double precision, parameter ::
     &   ptone = 0.1d0, half = 0.5d0, third = 0.33333333333333333333d0,
     &   pi = 3.141592653589793d0, root3 = 1.7320508075688772d0,
     &   onept5 = 1.5d0, pt1 = 0.1d0
    
      end module      
      
c
c ****************************************************************************
c *                                                                          *
c *    mod_crystals - contains all crystal data, struct for each crystal as  *
c *                   well as the array of crystals and some helper methods  *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 04/01/15 tjt                                     *
c *                                                                          *
c ****************************************************************************
c
      module crystal_data
            implicit integer (a-z)
      include 'param_def'
c                 Crystal array data structures
c
            type :: crystal
                  double precision, dimension(6,6) :: elast_stiff,
     &                                                elast_flex
                  double precision, dimension(max_slip_sys,3) :: ni,bi

c                 Material parameters
                  double precision :: e, nu, mu, harden_n, tau_a,
     &                                tau_hat_y, g_o_y, b, p_v, q_v,
     &                                p_y, q_y, boltz, 
     &                                eps_dot_o_y, t_o,
     &                                theta_o, tau_bar_o,
     &                                tau_hat_v, g_o_v,
     &                                eps_dot_o_v, k_o,
     &                                mu_o, D_o, tau_y, tau_v, iD_v,
     &                                voche_m !  yes it is spelled wrong
                  double precision :: u1, u2, u3, u4, u5, u6,
     &                                u7, u8, u9, u10
                  double precision :: C11, C12, C13, C33, C44, C55
c
                  double precision :: atol, atol1, rtol, rtol1,
     &                                xtol, xtol1
                  integer :: st_it(3)
                  integer ::  slip_type
c                             1) fcc
c                             2) bcc
c                             3) single
c                             6) fcc roters order
c                             7) bcc 12
c                             8) bcc 48
c                             9) hcp 6
c                             10) hcp 18
                  integer :: elastic_type
c                             1) isotropic
c                             2) cubic
c                             3) ti6242
                  integer :: nslip
                  integer :: num_hard
                  integer :: h_type
c                             1) voche or voce
c                             2) mts
c                             3) user
c                             4) ornl
c                             7) roters
c                             8) Armstrong-Frederick
c                             9) DJGM - Ti-6242
                  integer :: tang_calc, miter, gpp, method
c        General parameters
        double precision :: cp_001, cp_002, cp_003, cp_004, cp_005,
     &                      cp_006, cp_007, cp_008, cp_009, cp_010,
     &                      cp_011, cp_012, cp_013, cp_014, cp_015,
     &                      cp_016, cp_017, cp_018, cp_019, cp_020,
     &                      cp_021, cp_022, cp_023, cp_024, cp_025,
     &                      cp_026, cp_027, cp_028, cp_029, cp_030,
     &                      cp_031, cp_032, cp_033, cp_034, cp_035,
     &                      cp_036, cp_037, cp_038, cp_039, cp_040,
     &                      cp_041, cp_042, cp_043, cp_044, cp_045,
     &                      cp_046, cp_047, cp_048, cp_049, cp_050,
     &                      cp_051, cp_052, cp_053, cp_054, cp_055,
     &                      cp_056, cp_057, cp_058, cp_059, cp_060,
     &                      cp_061, cp_062, cp_063, cp_064, cp_065,
     &                      cp_066, cp_067, cp_068, cp_069, cp_070,
     &                      cp_071, cp_072, cp_073, cp_074, cp_075,
     &                      cp_076, cp_077, cp_078, cp_079, cp_080,
     &                      cp_081, cp_082, cp_083, cp_084, cp_085,
     &                      cp_086, cp_087, cp_088, cp_089, cp_090,
     &                      cp_091, cp_092, cp_093, cp_094, cp_095,
     &                      cp_096, cp_097, cp_098, cp_099, cp_100
c                 Solver flags
                  logical :: solver, strategy, gpall, alter_mode
c
                  logical :: valid

            end type crystal
c
            type(crystal), dimension(max_crystals) :: c_array
c
c           Data structures for mm11 Monte Carlo angles
            integer, allocatable :: mc_array(:,:)

            double precision, allocatable :: angle_input(:,:,:)
            double precision, allocatable :: simple_angles(:,:)
            integer, allocatable :: crystal_input(:,:)
            integer, allocatable :: data_offset(:)
            double precision :: cry_multiplier
            logical :: defined_crystal, srequired
            integer :: nangles

      contains
c                 Initializes crystal num with defaults
            subroutine initialize_new_crystal(num, out)
                  integer, intent(in) :: num, out
c                  
                  if (num .gt. max_crystals) then
                        write (out,*)
     &                   "Error: crystal number exceeds max crystals."
                        call die_gracefully
                  end if                  
                  if (num .lt. 1) then
                        write (out,*)
     &                   "Error: invalid crystal number."
                        call die_gracefully
                  end if
c
                  c_array(num)%slip_type = 1
                  c_array(num)%elastic_type = 1
                  c_array(num)%e = 69000.
                  c_array(num)%nu = 0.33
                  c_array(num)%mu = 69000./2./1.33
                  c_array(num)%C11 = 1.3281250d5
                  c_array(num)%C12 = 7.6562500d4
                  c_array(num)%C13 = 6.7187500d4
                  c_array(num)%C33 = 1.5937500d5
                  c_array(num)%C44 = 2.8125000d4
                  c_array(num)%C55 = 3.9843750d4
                  c_array(num)%harden_n = 20
                  c_array(num)%tau_a = 0
                  c_array(num)%tau_hat_y = -1.0
                  c_array(num)%g_o_y = -1.0
                  c_array(num)%tau_hat_v = -1.0
                  c_array(num)%g_o_v = -1.0
                  c_array(num)%b = 2.87e-7
                  c_array(num)%p_v = 0.5
                  c_array(num)%q_v = 2
                  c_array(num)%p_y = 0.5
                  c_array(num)%q_y = 2
                  c_array(num)%boltz = 1.3806E-20
                  c_array(num)%eps_dot_o_y = 1.0E10
                  c_array(num)%eps_dot_o_v = 1.0E10
                  c_array(num)%t_o = 294.
                  c_array(num)%theta_o = 100.0
                  c_array(num)%k_o = 0.0
                  c_array(num)%mu_o =  69000./2./1.33
                  c_array(num)%D_o = 0.0
                  c_array(num)%tau_y = 0.0
                  c_array(num)%tau_v = 0.0
                  c_array(num)%voche_m = 1.0
                  c_array(num)%iD_v = 0.0

                  c_array(num)%u1 = 0.0
                  c_array(num)%u2 = 0.0
                  c_array(num)%u3 = 0.0
                  c_array(num)%u4 = 0.0
                  c_array(num)%u5 = 0.0
                  c_array(num)%u6 = 0.0
c
                  c_array(num)%cp_001 = 0.d0
                  c_array(num)%cp_002 = 0.d0
                  c_array(num)%cp_003 = 0.d0
                  c_array(num)%cp_004 = 0.d0
                  c_array(num)%cp_005 = 0.d0
                  c_array(num)%cp_006 = 0.d0
                  c_array(num)%cp_007 = 0.d0
                  c_array(num)%cp_008 = 0.d0
                  c_array(num)%cp_009 = 0.d0
                  c_array(num)%cp_010 = 0.d0
                  c_array(num)%cp_011 = 0.d0
                  c_array(num)%cp_012 = 0.d0
                  c_array(num)%cp_013 = 0.d0
                  c_array(num)%cp_014 = 0.d0
                  c_array(num)%cp_015 = 0.d0
                  c_array(num)%cp_016 = 0.d0
                  c_array(num)%cp_017 = 0.d0
                  c_array(num)%cp_018 = 0.d0
                  c_array(num)%cp_019 = 0.d0
                  c_array(num)%cp_020 = 0.d0
                  c_array(num)%cp_021 = 0.d0
                  c_array(num)%cp_022 = 0.d0
                  c_array(num)%cp_023 = 0.d0
                  c_array(num)%cp_024 = 0.d0
                  c_array(num)%cp_025 = 0.d0
                  c_array(num)%cp_026 = 0.d0
                  c_array(num)%cp_027 = 0.d0
                  c_array(num)%cp_028 = 0.d0
                  c_array(num)%cp_029 = 0.d0
                  c_array(num)%cp_030 = 0.d0
                  c_array(num)%cp_031 = 0.d0
                  c_array(num)%cp_032 = 0.d0
                  c_array(num)%cp_033 = 0.d0
                  c_array(num)%cp_034 = 0.d0
                  c_array(num)%cp_035 = 0.d0
                  c_array(num)%cp_036 = 0.d0
                  c_array(num)%cp_037 = 0.d0
                  c_array(num)%cp_038 = 0.d0
                  c_array(num)%cp_039 = 0.d0
                  c_array(num)%cp_040 = 0.d0
                  c_array(num)%cp_041 = 0.d0
                  c_array(num)%cp_042 = 0.d0
                  c_array(num)%cp_043 = 0.d0
                  c_array(num)%cp_044 = 0.d0
                  c_array(num)%cp_045 = 0.d0
                  c_array(num)%cp_046 = 0.d0
                  c_array(num)%cp_047 = 0.d0
                  c_array(num)%cp_048 = 0.d0
                  c_array(num)%cp_049 = 0.d0
                  c_array(num)%cp_050 = 0.d0
                  c_array(num)%cp_051 = 0.d0
                  c_array(num)%cp_052 = 0.d0
                  c_array(num)%cp_053 = 0.d0
                  c_array(num)%cp_054 = 0.d0
                  c_array(num)%cp_055 = 0.d0
                  c_array(num)%cp_056 = 0.d0
                  c_array(num)%cp_057 = 0.d0
                  c_array(num)%cp_058 = 0.d0
                  c_array(num)%cp_059 = 0.d0
                  c_array(num)%cp_060 = 0.d0
                  c_array(num)%cp_061 = 0.d0
                  c_array(num)%cp_062 = 0.d0
                  c_array(num)%cp_063 = 0.d0
                  c_array(num)%cp_064 = 0.d0
                  c_array(num)%cp_065 = 0.d0
                  c_array(num)%cp_066 = 0.d0
                  c_array(num)%cp_067 = 0.d0
                  c_array(num)%cp_068 = 0.d0
                  c_array(num)%cp_069 = 0.d0
                  c_array(num)%cp_070 = 0.d0
                  c_array(num)%cp_071 = 0.d0
                  c_array(num)%cp_072 = 0.d0
                  c_array(num)%cp_073 = 0.d0
                  c_array(num)%cp_074 = 0.d0
                  c_array(num)%cp_075 = 0.d0
                  c_array(num)%cp_076 = 0.d0
                  c_array(num)%cp_077 = 0.d0
                  c_array(num)%cp_078 = 0.d0
                  c_array(num)%cp_079 = 0.d0
                  c_array(num)%cp_080 = 0.d0
                  c_array(num)%cp_081 = 0.d0
                  c_array(num)%cp_082 = 0.d0
                  c_array(num)%cp_083 = 0.d0
                  c_array(num)%cp_084 = 0.d0
                  c_array(num)%cp_085 = 0.d0
                  c_array(num)%cp_086 = 0.d0
                  c_array(num)%cp_087 = 0.d0
                  c_array(num)%cp_088 = 0.d0
                  c_array(num)%cp_089 = 0.d0
                  c_array(num)%cp_090 = 0.d0
                  c_array(num)%cp_091 = 0.d0
                  c_array(num)%cp_092 = 0.d0
                  c_array(num)%cp_093 = 0.d0
                  c_array(num)%cp_094 = 0.d0
                  c_array(num)%cp_095 = 0.d0
                  c_array(num)%cp_096 = 0.d0
                  c_array(num)%cp_097 = 0.d0
                  c_array(num)%cp_098 = 0.d0
                  c_array(num)%cp_099 = 0.d0
                  c_array(num)%cp_100 = 0.d0
c
                  c_array(num)%h_type = 1
                  c_array(num)%num_hard = 1
                  c_array(num)%tang_calc = 0
c            Solver parameters
                  c_array(num)%atol = 1.0d-5
                  c_array(num)%atol1 = 1.0d-5
                  c_array(num)%rtol = 5.0d-5
                  c_array(num)%rtol1 = 1.0d-5
                  c_array(num)%xtol = 1.0d-4
                  c_array(num)%xtol1 = 1.0d-4
                  c_array(num)%miter = 30
                  c_array(num)%solver = .true.
                  c_array(num)%strategy = .true.
                  c_array(num)%gpall = .false.
                  c_array(num)%gpp = 0
                  c_array(num)%method = 3
                  c_array(num)%st_it(1) = -2
                  c_array(num)%st_it(2) = -2
                  c_array(num)%st_it(3) = -2
c            Alternative model features, currently Voce only
                  c_array(num)%alter_mode = .false.

            end subroutine
c
c                 Inserts "derived" properties like slip system geometry
c                 and elasticity tensors
            subroutine finalize_new_crystal(num, out)
                  integer, intent(in) :: num, out
                  double precision :: f,e,v,u,a,c,ac1,ac2
                  double precision :: z0,f2,f3,f112,f211,
     &                  f123,f213,f312
                  double precision :: C11,C22,C44,C55,C66
                  double precision :: C12,C33,C13,C23
                  integer :: i,j,info
                  double precision, dimension(600) :: work
                  integer, dimension(600) :: lwork
                  integer, dimension(6) :: ipiv
c
                  if (num .gt. max_crystals) then
                        write (out,*)
     &                   "Error: crystal number exceeds max crystals."
                        call die_gracefully
                  end if                  
                  if (num .lt. 1) then
                        write (out,*)
     &                   "Error: invalid crystal number."
                        call die_gracefully
                  end if
c
c
                  if (c_array(num)%slip_type .eq. 1) then
                        c_array(num)%nslip = 12
                        f = 1/sqrt(2.0D0)
                        
                        c_array(num)%bi(1,1)=0
                        c_array(num)%bi(1,2)=-f
                        c_array(num)%bi(1,3)=f
                        c_array(num)%bi(2,1)=f
                        c_array(num)%bi(2,2)=0
                        c_array(num)%bi(2,3)=-f
                        c_array(num)%bi(3,1)=-f
                        c_array(num)%bi(3,2)=f
                        c_array(num)%bi(3,3)=0
                        c_array(num)%bi(4,1)=0
                        c_array(num)%bi(4,2)=f
                        c_array(num)%bi(4,3)=f
                        c_array(num)%bi(5,1)=f
                        c_array(num)%bi(5,2)=0
                        c_array(num)%bi(5,3)=f
                        c_array(num)%bi(6,1)=f
                        c_array(num)%bi(6,2)=-f
                        c_array(num)%bi(6,3)=0
                        c_array(num)%bi(7,1)=0
                        c_array(num)%bi(7,2)=-f
                        c_array(num)%bi(7,3)=f
                        c_array(num)%bi(8,1)=-f
                        c_array(num)%bi(8,2)=0
                        c_array(num)%bi(8,3)=-f
                        c_array(num)%bi(9,1)=f
                        c_array(num)%bi(9,2)=f
                        c_array(num)%bi(9,3)=0
                        c_array(num)%bi(10,1)=0
                        c_array(num)%bi(10,2)=f
                        c_array(num)%bi(10,3)=f
                        c_array(num)%bi(11,1)=f
                        c_array(num)%bi(11,2)=0
                        c_array(num)%bi(11,3)=-f
                        c_array(num)%bi(12,1)=-f
                        c_array(num)%bi(12,2)=-f
                        c_array(num)%bi(12,3)=0

                        f = 1/sqrt(3.0D0)
                        c_array(num)%ni(1,1)=f
                        c_array(num)%ni(1,2)=f
                        c_array(num)%ni(1,3)=f
                        c_array(num)%ni(2,1)=f
                        c_array(num)%ni(2,2)=f
                        c_array(num)%ni(2,3)=f
                        c_array(num)%ni(3,1)=f
                        c_array(num)%ni(3,2)=f
                        c_array(num)%ni(3,3)=f
                        c_array(num)%ni(4,1)=-f
                        c_array(num)%ni(4,2)=-f
                        c_array(num)%ni(4,3)=f
                        c_array(num)%ni(5,1)=-f
                        c_array(num)%ni(5,2)=-f
                        c_array(num)%ni(5,3)=f
                        c_array(num)%ni(6,1)=-f
                        c_array(num)%ni(6,2)=-f
                        c_array(num)%ni(6,3)=f
                        c_array(num)%ni(7,1)=-f
                        c_array(num)%ni(7,2)=f
                        c_array(num)%ni(7,3)=f
                        c_array(num)%ni(8,1)=-f
                        c_array(num)%ni(8,2)=f
                        c_array(num)%ni(8,3)=f
                        c_array(num)%ni(9,1)=-f
                        c_array(num)%ni(9,2)=f
                        c_array(num)%ni(9,3)=f
                        c_array(num)%ni(10,1)=f
                        c_array(num)%ni(10,2)=-f
                        c_array(num)%ni(10,3)=f
                        c_array(num)%ni(11,1)=f
                        c_array(num)%ni(11,2)=-f
                        c_array(num)%ni(11,3)=f
                        c_array(num)%ni(12,1)=f
                        c_array(num)%ni(12,2)=-f
                        c_array(num)%ni(12,3)=f
                  elseif (c_array(num)%slip_type .eq. 2) then
                        c_array(num)%nslip = 12
                        f = 1/sqrt(2.0D0)
                        
                        c_array(num)%ni(1,1)=0
                        c_array(num)%ni(1,2)=-f
                        c_array(num)%ni(1,3)=f
                        c_array(num)%ni(2,1)=f
                        c_array(num)%ni(2,2)=0
                        c_array(num)%ni(2,3)=-f
                        c_array(num)%ni(3,1)=-f
                        c_array(num)%ni(3,2)=f
                        c_array(num)%ni(3,3)=0
                        c_array(num)%ni(4,1)=0
                        c_array(num)%ni(4,2)=f
                        c_array(num)%ni(4,3)=f
                        c_array(num)%ni(5,1)=f
                        c_array(num)%ni(5,2)=0
                        c_array(num)%ni(5,3)=f
                        c_array(num)%ni(6,1)=f
                        c_array(num)%ni(6,2)=-f
                        c_array(num)%ni(6,3)=0
                        c_array(num)%ni(7,1)=0
                        c_array(num)%ni(7,2)=-f
                        c_array(num)%ni(7,3)=f
                        c_array(num)%ni(8,1)=-f
                        c_array(num)%ni(8,2)=0
                        c_array(num)%ni(8,3)=-f
                        c_array(num)%ni(9,1)=f
                        c_array(num)%ni(9,2)=f
                        c_array(num)%ni(9,3)=0
                        c_array(num)%ni(10,1)=0
                        c_array(num)%ni(10,2)=f
                        c_array(num)%ni(10,3)=f
                        c_array(num)%ni(11,1)=f
                        c_array(num)%ni(11,2)=0
                        c_array(num)%ni(11,3)=-f
                        c_array(num)%ni(12,1)=-f
                        c_array(num)%ni(12,2)=-f
                        c_array(num)%ni(12,3)=0

                        f = 1/sqrt(3.0D0)
                        c_array(num)%bi(1,1)=f
                        c_array(num)%bi(1,2)=f
                        c_array(num)%bi(1,3)=f
                        c_array(num)%bi(2,1)=f
                        c_array(num)%bi(2,2)=f
                        c_array(num)%bi(2,3)=f
                        c_array(num)%bi(3,1)=f
                        c_array(num)%bi(3,2)=f
                        c_array(num)%bi(3,3)=f
                        c_array(num)%bi(4,1)=-f
                        c_array(num)%bi(4,2)=-f
                        c_array(num)%bi(4,3)=f
                        c_array(num)%bi(5,1)=-f
                        c_array(num)%bi(5,2)=-f
                        c_array(num)%bi(5,3)=f
                        c_array(num)%bi(6,1)=-f
                        c_array(num)%bi(6,2)=-f
                        c_array(num)%bi(6,3)=f
                        c_array(num)%bi(7,1)=-f
                        c_array(num)%bi(7,2)=f
                        c_array(num)%bi(7,3)=f
                        c_array(num)%bi(8,1)=-f
                        c_array(num)%bi(8,2)=f
                        c_array(num)%bi(8,3)=f
                        c_array(num)%bi(9,1)=-f
                        c_array(num)%bi(9,2)=f
                        c_array(num)%bi(9,3)=f
                        c_array(num)%bi(10,1)=f
                        c_array(num)%bi(10,2)=-f
                        c_array(num)%bi(10,3)=f
                        c_array(num)%bi(11,1)=f
                        c_array(num)%bi(11,2)=-f
                        c_array(num)%bi(11,3)=f
                        c_array(num)%bi(12,1)=f
                        c_array(num)%bi(12,2)=-f
                        c_array(num)%bi(12,3)=f
                  elseif (c_array(num)%slip_type .eq. 3) then
                        c_array(num)%nslip = 1

                        c_array(num)%bi(1,1) = 1.0
                        c_array(num)%bi(1,2) = 0.0
                        c_array(num)%bi(1,3) = 0.0

                        c_array(num)%ni(1,1) = 0.0
                        c_array(num)%ni(1,2) = 1.0
                        c_array(num)%ni(1,3) = 0.0
                  elseif (c_array(num)%slip_type .eq. 6) then
                        c_array(num)%nslip = 12
                        f = 1/dsqrt(2.0D0)
c                     edge
                        c_array(num)%bi(1,1)=f
                        c_array(num)%bi(1,2)=-f
                        c_array(num)%bi(1,3)=0.d0
                        
                        c_array(num)%bi(2,1)=f
                        c_array(num)%bi(2,2)=0.d0
                        c_array(num)%bi(2,3)=-f
                        
                        c_array(num)%bi(3,1)=0.d0
                        c_array(num)%bi(3,2)=f
                        c_array(num)%bi(3,3)=-f
                        
                        c_array(num)%bi(4,1)=f
                        c_array(num)%bi(4,2)=f
                        c_array(num)%bi(4,3)=0.d0
                        
                        c_array(num)%bi(5,1)=f
                        c_array(num)%bi(5,2)=0.d0
                        c_array(num)%bi(5,3)=f
                        
                        c_array(num)%bi(6,1)=0.d0
                        c_array(num)%bi(6,2)=f
                        c_array(num)%bi(6,3)=-f
                        
                        c_array(num)%bi(7,1)=f
                        c_array(num)%bi(7,2)=f
                        c_array(num)%bi(7,3)=0.d0
                        
                        c_array(num)%bi(8,1)=f
                        c_array(num)%bi(8,2)=0.d0
                        c_array(num)%bi(8,3)=-f
                        
                        c_array(num)%bi(9,1)=0.d0
                        c_array(num)%bi(9,2)=f
                        c_array(num)%bi(9,3)=f
                        
                        c_array(num)%bi(10,1)=f
                        c_array(num)%bi(10,2)=-f
                        c_array(num)%bi(10,3)=0.d0
                        
                        c_array(num)%bi(11,1)=f
                        c_array(num)%bi(11,2)=0.d0
                        c_array(num)%bi(11,3)=f
                        
                        c_array(num)%bi(12,1)=0.d0
                        c_array(num)%bi(12,2)=f
                        c_array(num)%bi(12,3)=f
c                        
                        f = 1/dsqrt(3.0D0)
c                edge
                        c_array(num)%ni(1,1)=f
                        c_array(num)%ni(1,2)=f
                        c_array(num)%ni(1,3)=f
                        
                        c_array(num)%ni(2,1)=-f
                        c_array(num)%ni(2,2)=-f
                        c_array(num)%ni(2,3)=-f
                        
                        c_array(num)%ni(3,1)=f
                        c_array(num)%ni(3,2)=f
                        c_array(num)%ni(3,3)=f
                        
                        c_array(num)%ni(4,1)=f
                        c_array(num)%ni(4,2)=-f
                        c_array(num)%ni(4,3)=-f
                        
                        c_array(num)%ni(5,1)=-f
                        c_array(num)%ni(5,2)=f
                        c_array(num)%ni(5,3)=f
                        
                        c_array(num)%ni(6,1)=f
                        c_array(num)%ni(6,2)=-f
                        c_array(num)%ni(6,3)=-f
                        
                        c_array(num)%ni(7,1)=f
                        c_array(num)%ni(7,2)=-f
                        c_array(num)%ni(7,3)=f
                        
                        c_array(num)%ni(8,1)=f
                        c_array(num)%ni(8,2)=-f
                        c_array(num)%ni(8,3)=f
                        
                        c_array(num)%ni(9,1)=-f
                        c_array(num)%ni(9,2)=f
                        c_array(num)%ni(9,3)=-f
                        
                        c_array(num)%ni(10,1)=-f
                        c_array(num)%ni(10,2)=-f
                        c_array(num)%ni(10,3)=f
                        
                        c_array(num)%ni(11,1)=-f
                        c_array(num)%ni(11,2)=-f
                        c_array(num)%ni(11,3)=f
                        
                        c_array(num)%ni(12,1)=f
                        c_array(num)%ni(12,2)=f
                        c_array(num)%ni(12,3)=-f
                  elseif (c_array(num)%slip_type .eq. 7) then
                        c_array(num)%nslip = 12
                        z0 = 0.d0
                        f2 = 1.d0/sqrt(2.0D0)
                        f3 = 1.d0/sqrt(3.0D0)
                        f112 = 1.d0/sqrt(6.0D0)
                        f211 = 2.d0/sqrt(6.0D0)
                        f123 = 1.d0/sqrt(14.0D0)
                        f213 = 2.d0/sqrt(14.0D0)
                        f312 = 3.d0/sqrt(14.0D0)
                        !{110}<111>
                        c_array(num)%bi( 1,1)=-f3
                        c_array(num)%bi( 1,2)= f3
                        c_array(num)%bi( 1,3)= f3
                        
                        c_array(num)%bi( 2,1)= f3
                        c_array(num)%bi( 2,2)=-f3
                        c_array(num)%bi( 2,3)= f3
                        
                        c_array(num)%bi( 3,1)= f3
                        c_array(num)%bi( 3,2)= f3
                        c_array(num)%bi( 3,3)= f3
                        
                        c_array(num)%bi( 4,1)= f3
                        c_array(num)%bi( 4,2)= f3
                        c_array(num)%bi( 4,3)=-f3
                        
                        c_array(num)%bi( 5,1)= f3
                        c_array(num)%bi( 5,2)= f3
                        c_array(num)%bi( 5,3)=-f3
                        
                        c_array(num)%bi( 6,1)=-f3
                        c_array(num)%bi( 6,2)= f3
                        c_array(num)%bi( 6,3)= f3
                        
                        c_array(num)%bi( 7,1)= f3
                        c_array(num)%bi( 7,2)= f3
                        c_array(num)%bi( 7,3)= f3
                        
                        c_array(num)%bi( 8,1)= f3
                        c_array(num)%bi( 8,2)=-f3
                        c_array(num)%bi( 8,3)= f3
                        
                        c_array(num)%bi( 9,1)= f3
                        c_array(num)%bi( 9,2)= f3
                        c_array(num)%bi( 9,3)=-f3
                        
                        c_array(num)%bi(10,1)= f3
                        c_array(num)%bi(10,2)=-f3
                        c_array(num)%bi(10,3)= f3
                        
                        c_array(num)%bi(11,1)= f3
                        c_array(num)%bi(11,2)= f3
                        c_array(num)%bi(11,3)= f3
                        
                        c_array(num)%bi(12,1)=-f3
                        c_array(num)%bi(12,2)= f3
                        c_array(num)%bi(12,3)= f3
c                        
                        !{110}<111>
                        c_array(num)%ni( 1,1)= f2
                        c_array(num)%ni( 1,2)= f2
                        c_array(num)%ni( 1,3)= z0
                        
                        c_array(num)%ni( 2,1)= f2
                        c_array(num)%ni( 2,2)= f2
                        c_array(num)%ni( 2,3)= z0
                        
                        c_array(num)%ni( 3,1)= f2
                        c_array(num)%ni( 3,2)=-f2
                        c_array(num)%ni( 3,3)= z0
                        
                        c_array(num)%ni( 4,1)= f2
                        c_array(num)%ni( 4,2)=-f2
                        c_array(num)%ni( 4,3)= z0
                        
                        c_array(num)%ni( 5,1)= f2
                        c_array(num)%ni( 5,2)= z0
                        c_array(num)%ni( 5,3)= f2
                        
                        c_array(num)%ni( 6,1)= f2
                        c_array(num)%ni( 6,2)= z0
                        c_array(num)%ni( 6,3)= f2
                        
                        c_array(num)%ni( 7,1)= f2
                        c_array(num)%ni( 7,2)= z0
                        c_array(num)%ni( 7,3)=-f2
                        
                        c_array(num)%ni( 8,1)= f2
                        c_array(num)%ni( 8,2)= z0
                        c_array(num)%ni( 8,3)=-f2
                        
                        c_array(num)%ni( 9,1)= z0
                        c_array(num)%ni( 9,2)= f2
                        c_array(num)%ni( 9,3)= f2
                        
                        c_array(num)%ni(10,1)= z0
                        c_array(num)%ni(10,2)= f2
                        c_array(num)%ni(10,3)= f2
                        
                        c_array(num)%ni(11,1)= z0
                        c_array(num)%ni(11,2)= f2
                        c_array(num)%ni(11,3)=-f2
                        
                        c_array(num)%ni(12,1)= z0
                        c_array(num)%ni(12,2)= f2
                        c_array(num)%ni(12,3)=-f2
                  elseif (c_array(num)%slip_type .eq. 8) then ! BCC 48 sys
                        c_array(num)%nslip = 48
                        z0 = 0.d0
                        f2 = 1.d0/sqrt(2.0D0)
                        f3 = 1.d0/sqrt(3.0D0)
                        f112 = 1.d0/sqrt(6.0D0)
                        f211 = 2.d0/sqrt(6.0D0)
                        f123 = 1.d0/sqrt(14.0D0)
                        f213 = 2.d0/sqrt(14.0D0)
                        f312 = 3.d0/sqrt(14.0D0)
                        !{110}<111>
                        c_array(num)%bi( 1,1)=-f3
                        c_array(num)%bi( 1,2)= f3
                        c_array(num)%bi( 1,3)= f3
                        
                        c_array(num)%bi( 2,1)= f3
                        c_array(num)%bi( 2,2)=-f3
                        c_array(num)%bi( 2,3)= f3
                        
                        c_array(num)%bi( 3,1)= f3
                        c_array(num)%bi( 3,2)= f3
                        c_array(num)%bi( 3,3)= f3
                        
                        c_array(num)%bi( 4,1)= f3
                        c_array(num)%bi( 4,2)= f3
                        c_array(num)%bi( 4,3)=-f3
                        
                        c_array(num)%bi( 5,1)= f3
                        c_array(num)%bi( 5,2)= f3
                        c_array(num)%bi( 5,3)=-f3
                        
                        c_array(num)%bi( 6,1)=-f3
                        c_array(num)%bi( 6,2)= f3
                        c_array(num)%bi( 6,3)= f3
                        
                        c_array(num)%bi( 7,1)= f3
                        c_array(num)%bi( 7,2)= f3
                        c_array(num)%bi( 7,3)= f3
                        
                        c_array(num)%bi( 8,1)= f3
                        c_array(num)%bi( 8,2)=-f3
                        c_array(num)%bi( 8,3)= f3
                        
                        c_array(num)%bi( 9,1)= f3
                        c_array(num)%bi( 9,2)= f3
                        c_array(num)%bi( 9,3)=-f3
                        
                        c_array(num)%bi(10,1)= f3
                        c_array(num)%bi(10,2)=-f3
                        c_array(num)%bi(10,3)= f3
                        
                        c_array(num)%bi(11,1)= f3
                        c_array(num)%bi(11,2)= f3
                        c_array(num)%bi(11,3)= f3
                        
                        c_array(num)%bi(12,1)=-f3
                        c_array(num)%bi(12,2)= f3
                        c_array(num)%bi(12,3)= f3
                        !{112}<111>
                        c_array(num)%bi(13,1)= f3
                        c_array(num)%bi(13,2)= f3
                        c_array(num)%bi(13,3)=-f3
                        
                        c_array(num)%bi(14,1)= f3
                        c_array(num)%bi(14,2)=-f3
                        c_array(num)%bi(14,3)= f3
                        
                        c_array(num)%bi(15,1)=-f3
                        c_array(num)%bi(15,2)= f3
                        c_array(num)%bi(15,3)= f3
                        
                        c_array(num)%bi(16,1)= f3
                        c_array(num)%bi(16,2)= f3
                        c_array(num)%bi(16,3)= f3
                        
                        c_array(num)%bi(17,1)= f3
                        c_array(num)%bi(17,2)=-f3
                        c_array(num)%bi(17,3)= f3
                        
                        c_array(num)%bi(18,1)= f3
                        c_array(num)%bi(18,2)= f3
                        c_array(num)%bi(18,3)=-f3
                        
                        c_array(num)%bi(19,1)= f3
                        c_array(num)%bi(19,2)= f3
                        c_array(num)%bi(19,3)= f3
                        
                        c_array(num)%bi(20,1)=-f3
                        c_array(num)%bi(20,2)= f3
                        c_array(num)%bi(20,3)= f3
                        
                        c_array(num)%bi(21,1)=-f3
                        c_array(num)%bi(21,2)= f3
                        c_array(num)%bi(21,3)= f3
                        
                        c_array(num)%bi(22,1)= f3
                        c_array(num)%bi(22,2)= f3
                        c_array(num)%bi(22,3)= f3
                        
                        c_array(num)%bi(23,1)= f3
                        c_array(num)%bi(23,2)= f3
                        c_array(num)%bi(23,3)=-f3
                        
                        c_array(num)%bi(24,1)= f3
                        c_array(num)%bi(24,2)=-f3
                        c_array(num)%bi(24,3)= f3
                        !{123}<111>
                        c_array(num)%bi(25,1)= f3
                        c_array(num)%bi(25,2)= f3
                        c_array(num)%bi(25,3)=-f3
                        
                        c_array(num)%bi(26,1)= f3
                        c_array(num)%bi(26,2)=-f3
                        c_array(num)%bi(26,3)= f3
                        
                        c_array(num)%bi(27,1)=-f3
                        c_array(num)%bi(27,2)= f3
                        c_array(num)%bi(27,3)= f3
                        
                        c_array(num)%bi(28,1)= f3
                        c_array(num)%bi(28,2)= f3
                        c_array(num)%bi(28,3)= f3
                        
                        c_array(num)%bi(29,1)= f3
                        c_array(num)%bi(29,2)=-f3
                        c_array(num)%bi(29,3)= f3
                        
                        c_array(num)%bi(30,1)= f3
                        c_array(num)%bi(30,2)= f3
                        c_array(num)%bi(30,3)=-f3
                        
                        c_array(num)%bi(31,1)= f3
                        c_array(num)%bi(31,2)= f3
                        c_array(num)%bi(31,3)= f3
                        
                        c_array(num)%bi(32,1)=-f3
                        c_array(num)%bi(32,2)= f3
                        c_array(num)%bi(32,3)= f3
                        
                        c_array(num)%bi(33,1)= f3
                        c_array(num)%bi(33,2)= f3
                        c_array(num)%bi(33,3)=-f3
                        
                        c_array(num)%bi(34,1)= f3
                        c_array(num)%bi(34,2)=-f3
                        c_array(num)%bi(34,3)= f3
                        
                        c_array(num)%bi(35,1)=-f3
                        c_array(num)%bi(35,2)= f3
                        c_array(num)%bi(35,3)= f3
                        
                        c_array(num)%bi(36,1)= f3
                        c_array(num)%bi(36,2)= f3
                        c_array(num)%bi(36,3)= f3
                        
                        c_array(num)%bi(37,1)= f3
                        c_array(num)%bi(37,2)=-f3
                        c_array(num)%bi(37,3)= f3
                        
                        c_array(num)%bi(38,1)= f3
                        c_array(num)%bi(38,2)= f3
                        c_array(num)%bi(38,3)=-f3
                        
                        c_array(num)%bi(39,1)= f3
                        c_array(num)%bi(39,2)= f3
                        c_array(num)%bi(39,3)= f3
                        
                        c_array(num)%bi(40,1)=-f3
                        c_array(num)%bi(40,2)= f3
                        c_array(num)%bi(40,3)= f3
                        
                        c_array(num)%bi(41,1)=-f3
                        c_array(num)%bi(41,2)= f3
                        c_array(num)%bi(41,3)= f3
                        
                        c_array(num)%bi(42,1)= f3
                        c_array(num)%bi(42,2)= f3
                        c_array(num)%bi(42,3)= f3
                        
                        c_array(num)%bi(43,1)= f3
                        c_array(num)%bi(43,2)= f3
                        c_array(num)%bi(43,3)=-f3
                        
                        c_array(num)%bi(44,1)= f3
                        c_array(num)%bi(44,2)=-f3
                        c_array(num)%bi(44,3)= f3
                        
                        c_array(num)%bi(45,1)=-f3
                        c_array(num)%bi(45,2)= f3
                        c_array(num)%bi(45,3)= f3
                        
                        c_array(num)%bi(46,1)= f3
                        c_array(num)%bi(46,2)= f3
                        c_array(num)%bi(46,3)= f3
                        
                        c_array(num)%bi(47,1)= f3
                        c_array(num)%bi(47,2)= f3
                        c_array(num)%bi(47,3)=-f3
                        
                        c_array(num)%bi(48,1)= f3
                        c_array(num)%bi(48,2)=-f3
                        c_array(num)%bi(48,3)= f3
                        !{110}<111>
                        c_array(num)%ni( 1,1)= f2
                        c_array(num)%ni( 1,2)= f2
                        c_array(num)%ni( 1,3)= z0
                        
                        c_array(num)%ni( 2,1)= f2
                        c_array(num)%ni( 2,2)= f2
                        c_array(num)%ni( 2,3)= z0
                        
                        c_array(num)%ni( 3,1)= f2
                        c_array(num)%ni( 3,2)=-f2
                        c_array(num)%ni( 3,3)= z0
                        
                        c_array(num)%ni( 4,1)= f2
                        c_array(num)%ni( 4,2)=-f2
                        c_array(num)%ni( 4,3)= z0
                        
                        c_array(num)%ni( 5,1)= f2
                        c_array(num)%ni( 5,2)= z0
                        c_array(num)%ni( 5,3)= f2
                        
                        c_array(num)%ni( 6,1)= f2
                        c_array(num)%ni( 6,2)= z0
                        c_array(num)%ni( 6,3)= f2
                        
                        c_array(num)%ni( 7,1)= f2
                        c_array(num)%ni( 7,2)= z0
                        c_array(num)%ni( 7,3)=-f2
                        
                        c_array(num)%ni( 8,1)= f2
                        c_array(num)%ni( 8,2)= z0
                        c_array(num)%ni( 8,3)=-f2
                        
                        c_array(num)%ni( 9,1)= z0
                        c_array(num)%ni( 9,2)= f2
                        c_array(num)%ni( 9,3)= f2
                        
                        c_array(num)%ni(10,1)= z0
                        c_array(num)%ni(10,2)= f2
                        c_array(num)%ni(10,3)= f2
                        
                        c_array(num)%ni(11,1)= z0
                        c_array(num)%ni(11,2)= f2
                        c_array(num)%ni(11,3)=-f2
                        
                        c_array(num)%ni(12,1)= z0
                        c_array(num)%ni(12,2)= f2
                        c_array(num)%ni(12,3)=-f2
                        !{112}<111>
                        c_array(num)%ni(13,1)= f112
                        c_array(num)%ni(13,2)= f112
                        c_array(num)%ni(13,3)= f211
                        
                        c_array(num)%ni(14,1)=-f112
                        c_array(num)%ni(14,2)= f112
                        c_array(num)%ni(14,3)= f211
                        
                        c_array(num)%ni(15,1)= f112
                        c_array(num)%ni(15,2)=-f112
                        c_array(num)%ni(15,3)= f211
                        
                        c_array(num)%ni(16,1)= f112
                        c_array(num)%ni(16,2)= f112
                        c_array(num)%ni(16,3)=-f211
                        
                        c_array(num)%ni(17,1)= f112
                        c_array(num)%ni(17,2)= f211
                        c_array(num)%ni(17,3)= f112
                        
                        c_array(num)%ni(18,1)=-f112
                        c_array(num)%ni(18,2)= f211
                        c_array(num)%ni(18,3)= f112
                        
                        c_array(num)%ni(19,1)= f112
                        c_array(num)%ni(19,2)=-f211
                        c_array(num)%ni(19,3)= f112
                        
                        c_array(num)%ni(20,1)= f112
                        c_array(num)%ni(20,2)= f211
                        c_array(num)%ni(20,3)=-f112
                        
                        c_array(num)%ni(21,1)= f211
                        c_array(num)%ni(21,2)= f112
                        c_array(num)%ni(21,3)= f112
                        
                        c_array(num)%ni(22,1)=-f211
                        c_array(num)%ni(22,2)= f112
                        c_array(num)%ni(22,3)= f112
                        
                        c_array(num)%ni(23,1)= f211
                        c_array(num)%ni(23,2)=-f112
                        c_array(num)%ni(23,3)= f112
                        
                        c_array(num)%ni(24,1)= f211
                        c_array(num)%ni(24,2)= f112
                        c_array(num)%ni(24,3)=-f112
                        !{123}<111>
                        c_array(num)%ni(25,1)= f123
                        c_array(num)%ni(25,2)= f213
                        c_array(num)%ni(25,3)= f312
                        
                        c_array(num)%ni(26,1)=-f123
                        c_array(num)%ni(26,2)= f213
                        c_array(num)%ni(26,3)= f312
                        
                        c_array(num)%ni(27,1)= f123
                        c_array(num)%ni(27,2)=-f213
                        c_array(num)%ni(27,3)= f312
                        
                        c_array(num)%ni(28,1)= f123
                        c_array(num)%ni(28,2)= f213
                        c_array(num)%ni(28,3)=-f312
                        
                        c_array(num)%ni(29,1)= f123
                        c_array(num)%ni(29,2)= f312
                        c_array(num)%ni(29,3)= f213
                        
                        c_array(num)%ni(30,1)=-f123
                        c_array(num)%ni(30,2)= f312
                        c_array(num)%ni(30,3)= f213
                        
                        c_array(num)%ni(31,1)= f123
                        c_array(num)%ni(31,2)=-f312
                        c_array(num)%ni(31,3)= f213
                        
                        c_array(num)%ni(32,1)= f123
                        c_array(num)%ni(32,2)= f312
                        c_array(num)%ni(32,3)=-f213
                        
                        c_array(num)%ni(33,1)= f213
                        c_array(num)%ni(33,2)= f123
                        c_array(num)%ni(33,3)= f312
                        
                        c_array(num)%ni(34,1)=-f213
                        c_array(num)%ni(34,2)= f123
                        c_array(num)%ni(34,3)= f312
                        
                        c_array(num)%ni(35,1)= f213
                        c_array(num)%ni(35,2)=-f123
                        c_array(num)%ni(35,3)= f312
                        
                        c_array(num)%ni(36,1)= f213
                        c_array(num)%ni(36,2)= f123
                        c_array(num)%ni(36,3)=-f312

                        c_array(num)%ni(37,1)= f213
                        c_array(num)%ni(37,2)= f312
                        c_array(num)%ni(37,3)= f123
                        
                        c_array(num)%ni(38,1)=-f213
                        c_array(num)%ni(38,2)= f312
                        c_array(num)%ni(38,3)= f123
                        
                        c_array(num)%ni(39,1)= f213
                        c_array(num)%ni(39,2)=-f312
                        c_array(num)%ni(39,3)= f123
                        
                        c_array(num)%ni(40,1)= f213
                        c_array(num)%ni(40,2)= f312
                        c_array(num)%ni(40,3)=-f123
                        
                        c_array(num)%ni(41,1)= f312
                        c_array(num)%ni(41,2)= f123
                        c_array(num)%ni(41,3)= f213
                        
                        c_array(num)%ni(42,1)=-f312
                        c_array(num)%ni(42,2)= f123
                        c_array(num)%ni(42,3)= f213
                        
                        c_array(num)%ni(43,1)= f312
                        c_array(num)%ni(43,2)=-f123
                        c_array(num)%ni(43,3)= f213
                        
                        c_array(num)%ni(44,1)= f312
                        c_array(num)%ni(44,2)= f123
                        c_array(num)%ni(44,3)=-f213
                        
                        c_array(num)%ni(45,1)= f312
                        c_array(num)%ni(45,2)= f213
                        c_array(num)%ni(45,3)= f123
                        
                        c_array(num)%ni(46,1)=-f312
                        c_array(num)%ni(46,2)= f213
                        c_array(num)%ni(46,3)= f123
                        
                        c_array(num)%ni(47,1)= f312
                        c_array(num)%ni(47,2)=-f213
                        c_array(num)%ni(47,3)= f123
                        
                        c_array(num)%ni(48,1)= f312
                        c_array(num)%ni(48,2)= f213
                        c_array(num)%ni(48,3)=-f123

                        
                  elseif (c_array(num)%slip_type .eq. 9) then !Ran HCP
                        c_array(num)%nslip = 6
                        ! material constant of hcp
                        a = 0.295d0
                        c = 0.468d0
                        ac1 = sqrt(c**2+a**2)
                        ac2 = sqrt(4.d0*c**2+3.d0*a**2)
                        
                        ! Basal Slip Systems {0 0 0 1}<1 1 -2 0>
                        ! B1
                        c_array(num)%bi(1,1)=0.5d0
                        c_array(num)%bi(1,2)=-sqrt(3.d0)/2.d0
                        c_array(num)%bi(1,3)=0.d0
                        ! B2
                        c_array(num)%bi(2,1)=0.5d0
                        c_array(num)%bi(2,2)=sqrt(3.d0)/2.d0
                        c_array(num)%bi(2,3)=0.d0
                        ! B3
                        c_array(num)%bi(3,1)=-1.d0
                        c_array(num)%bi(3,2)=0.d0
                        c_array(num)%bi(3,3)=0.d0
                        
                        ! Prismatic Slip Systems {1 0 -1 0}<1 1 -2 0>
                        ! P1
                        c_array(num)%bi(4,1)=1.d0
                        c_array(num)%bi(4,2)=0.d0
                        c_array(num)%bi(4,3)=0.d0
                        ! P2
                        c_array(num)%bi(5,1)=0.5d0
                        c_array(num)%bi(5,2)=sqrt(3.d0)/2.d0
                        c_array(num)%bi(5,3)=0.d0
                        ! P3
                        c_array(num)%bi(6,1)=-0.5d0
                        c_array(num)%bi(6,2)=sqrt(3.d0)/2.d0
                        c_array(num)%bi(6,3)=0.d0
                        
c                         ! Pyramidal <a> Slip Systems {1 0 -1 1}<1 1 -2 0>
c                         ! R1
c                         c_array(num)%bi(7,1)=1
c                         c_array(num)%bi(7,2)=0
c                         c_array(num)%bi(7,3)=0
c                         ! R2
c                         c_array(num)%bi(8,1)=0.5
c                         c_array(num)%bi(8,2)=sqrt(3.d0)/2
c                         c_array(num)%bi(8,3)=0
c                         c R3
c                         c_array(num)%bi(9,1)=-0.5
c                         c_array(num)%bi(9,2)=sqrt(3.d0)/2
c                         c_array(num)%bi(9,3)=0
c                         c R4
c                         c_array(num)%bi(10,1)=-1
c                         c_array(num)%bi(10,2)=0
c                         c_array(num)%bi(10,3)=0
c                         c R5
c                         c_array(num)%bi(11,1)=-0.5
c                         c_array(num)%bi(11,2)=-sqrt(3.d0)/2
c                         c_array(num)%bi(11,3)=0
c                         c R6
c                         c_array(num)%bi(12,1)=0.5
c                         c_array(num)%bi(12,2)=-sqrt(3.d0)/2
c                         c_array(num)%bi(12,3)=0
c                         
c                         c Pyramidal <c+a> Slip Systems {1 0 -1 1}<1 1 -2 3>
c                         c R7
c                         c_array(num)%bi(13,1)=a/(2*ac1)
c                         c_array(num)%bi(13,2)=sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(13,3)=c/ac1
c                         c R8
c                         c_array(num)%bi(14,1)=-a/(2*ac1)
c                         c_array(num)%bi(14,2)=sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(14,3)=c/ac1
c                         c R9
c                         c_array(num)%bi(15,1)=-a/ac1
c                         c_array(num)%bi(15,2)=0
c                         c_array(num)%bi(15,3)=c/ac1
c                         c R10
c                         c_array(num)%bi(16,1)=-a/(2*ac1)
c                         c_array(num)%bi(16,2)=-sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(16,3)=c/ac1
c                         c R11
c                         c_array(num)%bi(17,1)=a/(2*ac1)
c                         c_array(num)%bi(17,2)=-sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(17,3)=c/ac1
c                         c R12
c                         c_array(num)%bi(18,1)=a/ac1
c                         c_array(num)%bi(18,2)=0
c                         c_array(num)%bi(18,3)=c/ac1
c                         c R13
c                         c_array(num)%bi(19,1)=-a/(2*ac1)
c                         c_array(num)%bi(19,2)=sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(19,3)=c/ac1
c                         c R14
c                         c_array(num)%bi(20,1)=-a/ac1
c                         c_array(num)%bi(20,2)=0
c                         c_array(num)%bi(20,3)=c/ac1
c                         c R15
c                         c_array(num)%bi(21,1)=-a/(2*ac1)
c                         c_array(num)%bi(21,2)=-sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(21,3)=c/ac1
c                         c R16
c                         c_array(num)%bi(22,1)=a/(2*ac1)
c                         c_array(num)%bi(22,2)=-sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(22,3)=c/ac1
c                         c R17
c                         c_array(num)%bi(23,1)=a/ac1
c                         c_array(num)%bi(23,2)=0
c                         c_array(num)%bi(23,3)=c/ac1
c                         c R18
c                         c_array(num)%bi(24,1)=a/(2*ac1)
c                         c_array(num)%bi(24,2)=sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(24,3)=c/ac1
c 
c                         c Pyramidal <c+a> Slip Systems {1 1 -2 2}<1 1 -2 3>
c                         c R19
c                         c_array(num)%bi(25,1)=-a/(2*ac1)
c                         c_array(num)%bi(25,2)=sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(25,3)=c/ac1
c                         c R20
c                         c_array(num)%bi(26,1)=-a/ac1
c                         c_array(num)%bi(26,2)=0
c                         c_array(num)%bi(26,3)=c/ac1
c                         c R21
c                         c_array(num)%bi(27,1)=-a/(2*ac1)
c                         c_array(num)%bi(27,2)=-sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(27,3)=c/ac1
c                         c R22
c                         c_array(num)%bi(28,1)=a/(2*ac1)
c                         c_array(num)%bi(28,2)=-sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(28,3)=c/ac1
c                         c R23
c                         c_array(num)%bi(29,1)=a/ac1
c                         c_array(num)%bi(29,2)=0
c                         c_array(num)%bi(29,3)=c/ac1
c                         c R24
c                         c_array(num)%bi(30,1)=a/(2*ac1)
c                         c_array(num)%bi(30,2)=sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(30,3)=c/ac1
c -------------------------------------------------------------------------
                        ! Basal Slip Systems {0 0 0 1}<1 1 -2 0>
                        ! B1
                        c_array(num)%ni(1,1)=0.d0
                        c_array(num)%ni(1,2)=0.d0
                        c_array(num)%ni(1,3)=1.d0
                        ! B2
                        c_array(num)%ni(2,1)=0.d0
                        c_array(num)%ni(2,2)=0.d0
                        c_array(num)%ni(2,3)=1.d0
                        ! B3
                        c_array(num)%ni(3,1)=0.d0
                        c_array(num)%ni(3,2)=0.d0
                        c_array(num)%ni(3,3)=1.d0
                        
                        ! Prismatic Slip Systems {1 0 -1 0}<1 1 -2 0>
                        ! P1
                        c_array(num)%ni(4,1)=0.d0
                        c_array(num)%ni(4,2)=1.d0
                        c_array(num)%ni(4,3)=0.d0
                        ! P2
                        c_array(num)%ni(5,1)=-sqrt(3.d0)/2.d0
                        c_array(num)%ni(5,2)=0.5d0
                        c_array(num)%ni(5,3)=0.d0
                        ! P3
                        c_array(num)%ni(6,1)=-sqrt(3.d0)/2.d0
                        c_array(num)%ni(6,2)=-0.5d0
                        c_array(num)%ni(6,3)=0.d0
                        
c                         c Pyramidal <a> Slip Systems {1 0 -1 1}<1 1 -2 0>
c                         c R1
c                         c_array(num)%ni(7,1)=0
c                         c_array(num)%ni(7,2)=-2*c/ac2
c                         c_array(num)%ni(7,3)=sqrt(3.d0)*a/ac2
c                         c R2
c                         c_array(num)%ni(8,1)=sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(8,2)=-c/ac2
c                         c_array(num)%ni(8,3)=sqrt(3.d0)*a/ac2
c                         c R3
c                         c_array(num)%ni(9,1)=sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(9,2)=c/ac2
c                         c_array(num)%ni(9,3)=sqrt(3.d0)*a/ac2
c                         c R4
c                         c_array(num)%ni(10,1)=0
c                         c_array(num)%ni(10,2)=2*c/ac2
c                         c_array(num)%ni(10,3)=sqrt(3.d0)*a/ac2
c                         c R5
c                         c_array(num)%ni(11,1)=-sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(11,2)=c/ac2
c                         c_array(num)%ni(11,3)=sqrt(3.d0)*a/ac2
c                         c R6
c                         c_array(num)%ni(12,1)=-sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(12,2)=-c/ac2
c                         c_array(num)%ni(12,3)=sqrt(3.d0)*a/ac2
c                         
c                         c Pyramidal <c+a> Slip Systems {1 0 -1 1}<1 1 -2 3>
c                         c R7
c                         c_array(num)%ni(13,1)=0
c                         c_array(num)%ni(13,2)=-2*c/ac2
c                         c_array(num)%ni(13,3)=sqrt(3.d0)*a/ac2
c                         c R8
c                         c_array(num)%ni(14,1)=sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(14,2)=-c/ac2
c                         c_array(num)%ni(14,3)=sqrt(3.d0)*a/ac2
c                         c R9
c                         c_array(num)%ni(15,1)=sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(15,2)=c/ac2
c                         c_array(num)%ni(15,3)=sqrt(3.d0)*a/ac2
c                         c R10
c                         c_array(num)%ni(16,1)=0
c                         c_array(num)%ni(16,2)=2*c/ac2
c                         c_array(num)%ni(16,3)=sqrt(3.d0)*a/ac2
c                         c R11
c                         c_array(num)%ni(17,1)=-sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(17,2)=c/ac2
c                         c_array(num)%ni(17,3)=sqrt(3.d0)*a/ac2
c                         c R12
c                         c_array(num)%ni(18,1)=-sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(18,2)=-c/ac2
c                         c_array(num)%ni(18,3)=sqrt(3.d0)*a/ac2
c                         c R13
c                         c_array(num)%ni(19,1)=0
c                         c_array(num)%ni(19,2)=-2*c/ac2
c                         c_array(num)%ni(19,3)=sqrt(3.d0)*a/ac2
c                         c R14
c                         c_array(num)%ni(20,1)=sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(20,2)=-c/ac2
c                         c_array(num)%ni(20,3)=sqrt(3.d0)*a/ac2
c                         c R15
c                         c_array(num)%ni(21,1)=sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(21,2)=c/ac2
c                         c_array(num)%ni(21,3)=sqrt(3.d0)*a/ac2
c                         c R16
c                         c_array(num)%ni(22,1)=0
c                         c_array(num)%ni(22,2)=2*c/ac2
c                         c_array(num)%ni(22,3)=sqrt(3.d0)*a/ac2
c                         c R17
c                         c_array(num)%ni(23,1)=-sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(23,2)=c/ac2
c                         c_array(num)%ni(23,3)=sqrt(3.d0)*a/ac2
c                         c R18
c                         c_array(num)%ni(24,1)=-sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(24,2)=-c/ac2
c                         c_array(num)%ni(24,3)=sqrt(3.d0)*a/ac2
c 
c                         c Pyramidal <c+a> Slip Systems {1 1 -2 2}<1 1 -2 3>
c                         c R19
c                         c_array(num)%ni(25,1)=c/(2*ac1)
c                         c_array(num)%ni(25,2)=-sqrt(3.d0)*c/(2*ac1)
c                         c_array(num)%ni(25,3)=a/ac1
c                         c R20
c                         c_array(num)%ni(26,1)=c/ac1
c                         c_array(num)%ni(26,2)=0
c                         c_array(num)%ni(26,3)=a/ac1
c                         c R21
c                         c_array(num)%ni(27,1)=c/(2*ac1)
c                         c_array(num)%ni(27,2)=sqrt(3.d0)*c/(2*ac1)
c                         c_array(num)%ni(27,3)=a/ac1
c                         c R22
c                         c_array(num)%ni(28,1)=-c/(2*ac1)
c                         c_array(num)%ni(28,2)=sqrt(3.d0)*c/(2*ac1)
c                         c_array(num)%ni(28,3)=a/ac1
c                         c R23
c                         c_array(num)%ni(29,1)=-c/ac1
c                         c_array(num)%ni(29,2)=0
c                         c_array(num)%ni(29,3)=a/ac1
c                         c R24
c                         c_array(num)%ni(30,1)=-c/(2*ac1)
c                         c_array(num)%ni(30,2)=-sqrt(3.d0)*c/(2*ac1)
c                         c_array(num)%ni(30,3)=a/ac1

c HCP: basal, prismatic, 1st-order pyramidal
                  elseif (c_array(num)%slip_type .eq. 10) then
                        c_array(num)%nslip = 6+12
                        ! material constant of hcp
                        a = 0.295d0
                        c = 0.468d0
                        ac1 = sqrt(c**2+a**2)
                        ac2 = sqrt(4.d0*c**2+3.d0*a**2)
                        
                        ! Basal Slip Systems {0 0 0 1}<1 1 -2 0>
                        ! B1
                        c_array(num)%bi(1,1)=0.5d0
                        c_array(num)%bi(1,2)=-sqrt(3.d0)/2.d0
                        c_array(num)%bi(1,3)=0.d0
                        ! B2
                        c_array(num)%bi(2,1)=0.5d0
                        c_array(num)%bi(2,2)=sqrt(3.d0)/2.d0
                        c_array(num)%bi(2,3)=0.d0
                        ! B3
                        c_array(num)%bi(3,1)=-1.d0
                        c_array(num)%bi(3,2)=0.d0
                        c_array(num)%bi(3,3)=0.d0
                        
                        ! Prismatic Slip Systems {1 0 -1 0}<1 1 -2 0>
                        ! P1
                        c_array(num)%bi(4,1)=1.d0
                        c_array(num)%bi(4,2)=0.d0
                        c_array(num)%bi(4,3)=0.d0
                        ! P2
                        c_array(num)%bi(5,1)=0.5d0
                        c_array(num)%bi(5,2)=sqrt(3.d0)/2.d0
                        c_array(num)%bi(5,3)=0.d0
                        ! P3
                        c_array(num)%bi(6,1)=-0.5d0
                        c_array(num)%bi(6,2)=sqrt(3.d0)/2.d0
                        c_array(num)%bi(6,3)=0.d0
                        
c                         c Pyramidal <a> Slip Systems {1 0 -1 1}<1 1 -2 0>
c                         c R1
c                         c_array(num)%bi(7,1)=1
c                         c_array(num)%bi(7,2)=0
c                         c_array(num)%bi(7,3)=0
c                         c R2
c                         c_array(num)%bi(8,1)=0.5
c                         c_array(num)%bi(8,2)=sqrt(3.d0)/2
c                         c_array(num)%bi(8,3)=0
c                         c R3
c                         c_array(num)%bi(9,1)=-0.5
c                         c_array(num)%bi(9,2)=sqrt(3.d0)/2
c                         c_array(num)%bi(9,3)=0
c                         c R4
c                         c_array(num)%bi(10,1)=-1
c                         c_array(num)%bi(10,2)=0
c                         c_array(num)%bi(10,3)=0
c                         c R5
c                         c_array(num)%bi(11,1)=-0.5
c                         c_array(num)%bi(11,2)=-sqrt(3.d0)/2
c                         c_array(num)%bi(11,3)=0
c                         c R6
c                         c_array(num)%bi(12,1)=0.5
c                         c_array(num)%bi(12,2)=-sqrt(3.d0)/2
c                         c_array(num)%bi(12,3)=0
c                         
                        ! Pyramidal <c+a> Slip Systems {1 0 -1 1}<1 1 -2 3>
                        ! R7
                        c_array(num)%bi(13-6,1)=a/(2.d0*ac1)
                        c_array(num)%bi(13-6,2)=sqrt(3.d0)*a/(2.d0*ac1)
                        c_array(num)%bi(13-6,3)=c/ac1
                        ! R8
                        c_array(num)%bi(14-6,1)=-a/(2.d0*ac1)
                        c_array(num)%bi(14-6,2)=sqrt(3.d0)*a/(2.d0*ac1)
                        c_array(num)%bi(14-6,3)=c/ac1
                        ! R9
                        c_array(num)%bi(15-6,1)=-a/ac1
                        c_array(num)%bi(15-6,2)=0.d0
                        c_array(num)%bi(15-6,3)=c/ac1
                        ! R10
                        c_array(num)%bi(16-6,1)=-a/(2.d0*ac1)
                        c_array(num)%bi(16-6,2)=-sqrt(3.d0)*a/(2.d0*ac1)
                        c_array(num)%bi(16-6,3)=c/ac1
                        ! R11
                        c_array(num)%bi(17-6,1)=a/(2.d0*ac1)
                        c_array(num)%bi(17-6,2)=-sqrt(3.d0)*a/(2.d0*ac1)
                        c_array(num)%bi(17-6,3)=c/ac1
                        ! R12
                        c_array(num)%bi(18-6,1)=a/ac1
                        c_array(num)%bi(18-6,2)=0.d0
                        c_array(num)%bi(18-6,3)=c/ac1
                        ! R13
                        c_array(num)%bi(19-6,1)=-a/(2.d0*ac1)
                        c_array(num)%bi(19-6,2)=sqrt(3.d0)*a/(2.d0*ac1)
                        c_array(num)%bi(19-6,3)=c/ac1
                        ! R14
                        c_array(num)%bi(20-6,1)=-a/ac1
                        c_array(num)%bi(20-6,2)=0.d0
                        c_array(num)%bi(20-6,3)=c/ac1
                        ! R15
                        c_array(num)%bi(21-6,1)=-a/(2.d0*ac1)
                        c_array(num)%bi(21-6,2)=-sqrt(3.d0)*a/(2.d0*ac1)
                        c_array(num)%bi(21-6,3)=c/ac1
                        ! R16
                        c_array(num)%bi(22-6,1)=a/(2.d0*ac1)
                        c_array(num)%bi(22-6,2)=-sqrt(3.d0)*a/(2.d0*ac1)
                        c_array(num)%bi(22-6,3)=c/ac1
                        ! R17
                        c_array(num)%bi(23-6,1)=a/ac1
                        c_array(num)%bi(23-6,2)=0.d0
                        c_array(num)%bi(23-6,3)=c/ac1
                        ! R18
                        c_array(num)%bi(24-6,1)=a/(2.d0*ac1)
                        c_array(num)%bi(24-6,2)=sqrt(3.d0)*a/(2.d0*ac1)
                        c_array(num)%bi(24-6,3)=c/ac1
c 
c                         c Pyramidal <c+a> Slip Systems {1 1 -2 2}<1 1 -2 3>
c                         c R19
c                         c_array(num)%bi(25,1)=-a/(2*ac1)
c                         c_array(num)%bi(25,2)=sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(25,3)=c/ac1
c                         c R20
c                         c_array(num)%bi(26,1)=-a/ac1
c                         c_array(num)%bi(26,2)=0
c                         c_array(num)%bi(26,3)=c/ac1
c                         c R21
c                         c_array(num)%bi(27,1)=-a/(2*ac1)
c                         c_array(num)%bi(27,2)=-sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(27,3)=c/ac1
c                         c R22
c                         c_array(num)%bi(28,1)=a/(2*ac1)
c                         c_array(num)%bi(28,2)=-sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(28,3)=c/ac1
c                         c R23
c                         c_array(num)%bi(29,1)=a/ac1
c                         c_array(num)%bi(29,2)=0
c                         c_array(num)%bi(29,3)=c/ac1
c                         c R24
c                         c_array(num)%bi(30,1)=a/(2*ac1)
c                         c_array(num)%bi(30,2)=sqrt(3.d0)*a/(2*ac1)
c                         c_array(num)%bi(30,3)=c/ac1
c -------------------------------------------------------------------------
                        ! Basal Slip Systems {0 0 0 1}<1 1 -2 0>
                        ! B1
                        c_array(num)%ni(1,1)=0.d0
                        c_array(num)%ni(1,2)=0.d0
                        c_array(num)%ni(1,3)=1.d0
                        ! B2
                        c_array(num)%ni(2,1)=0.d0
                        c_array(num)%ni(2,2)=0.d0
                        c_array(num)%ni(2,3)=1.d0
                        ! B3
                        c_array(num)%ni(3,1)=0.d0
                        c_array(num)%ni(3,2)=0.d0
                        c_array(num)%ni(3,3)=1.d0
                        
                        ! Prismatic Slip Systems {1 0 -1 0}<1 1 -2 0>
                        ! P1
                        c_array(num)%ni(4,1)=0.d0
                        c_array(num)%ni(4,2)=1.d0
                        c_array(num)%ni(4,3)=0.d0
                        ! P2
                        c_array(num)%ni(5,1)=-sqrt(3.d0)/2.d0
                        c_array(num)%ni(5,2)=0.5d0
                        c_array(num)%ni(5,3)=0.d0
                        ! P3
                        c_array(num)%ni(6,1)=-sqrt(3.d0)/2.d0
                        c_array(num)%ni(6,2)=-0.5d0
                        c_array(num)%ni(6,3)=0.d0
                        
c                         c Pyramidal <a> Slip Systems {1 0 -1 1}<1 1 -2 0>
c                         c R1
c                         c_array(num)%ni(7,1)=0.d0
c                         c_array(num)%ni(7,2)=-2.d0*c/ac2
c                         c_array(num)%ni(7,3)=sqrt(3.d0)*a/ac2
c                         c R2
c                         c_array(num)%ni(8,1)=sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(8,2)=-c/ac2
c                         c_array(num)%ni(8,3)=sqrt(3.d0)*a/ac2
c                         c R3
c                         c_array(num)%ni(9,1)=sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(9,2)=c/ac2
c                         c_array(num)%ni(9,3)=sqrt(3.d0)*a/ac2
c                         c R4
c                         c_array(num)%ni(10,1)=0.d0
c                         c_array(num)%ni(10,2)=2.d0*c/ac2
c                         c_array(num)%ni(10,3)=sqrt(3.d0)*a/ac2
c                         c R5
c                         c_array(num)%ni(11,1)=-sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(11,2)=c/ac2
c                         c_array(num)%ni(11,3)=sqrt(3.d0)*a/ac2
c                         c R6
c                         c_array(num)%ni(12,1)=-sqrt(3.d0)*c/ac2
c                         c_array(num)%ni(12,2)=-c/ac2
c                         c_array(num)%ni(12,3)=sqrt(3.d0)*a/ac2
c                         
                        ! Pyramidal <c+a> Slip Systems {1 0 -1 1}<1 1 -2 3>
                        ! R7
                        c_array(num)%ni(13-6,1)=0.d0
                        c_array(num)%ni(13-6,2)=-2.d0*c/ac2
                        c_array(num)%ni(13-6,3)=sqrt(3.d0)*a/ac2
                        ! R8
                        c_array(num)%ni(14-6,1)=sqrt(3.d0)*c/ac2
                        c_array(num)%ni(14-6,2)=-c/ac2
                        c_array(num)%ni(14-6,3)=sqrt(3.d0)*a/ac2
                        ! R9
                        c_array(num)%ni(15-6,1)=sqrt(3.d0)*c/ac2
                        c_array(num)%ni(15-6,2)=c/ac2
                        c_array(num)%ni(15-6,3)=sqrt(3.d0)*a/ac2
                        ! R10
                        c_array(num)%ni(16-6,1)=0.d0
                        c_array(num)%ni(16-6,2)=2.d0*c/ac2
                        c_array(num)%ni(16-6,3)=sqrt(3.d0)*a/ac2
                        ! R11
                        c_array(num)%ni(17-6,1)=-sqrt(3.d0)*c/ac2
                        c_array(num)%ni(17-6,2)=c/ac2
                        c_array(num)%ni(17-6,3)=sqrt(3.d0)*a/ac2
                        ! R12
                        c_array(num)%ni(18-6,1)=-sqrt(3.d0)*c/ac2
                        c_array(num)%ni(18-6,2)=-c/ac2
                        c_array(num)%ni(18-6,3)=sqrt(3.d0)*a/ac2
                        ! R13
                        c_array(num)%ni(19-6,1)=0.d0
                        c_array(num)%ni(19-6,2)=-2.d0*c/ac2
                        c_array(num)%ni(19-6,3)=sqrt(3.d0)*a/ac2
                        ! R14
                        c_array(num)%ni(20-6,1)=sqrt(3.d0)*c/ac2
                        c_array(num)%ni(20-6,2)=-c/ac2
                        c_array(num)%ni(20-6,3)=sqrt(3.d0)*a/ac2
                        ! R15
                        c_array(num)%ni(21-6,1)=sqrt(3.d0)*c/ac2
                        c_array(num)%ni(21-6,2)=c/ac2
                        c_array(num)%ni(21-6,3)=sqrt(3.d0)*a/ac2
                        ! R16
                        c_array(num)%ni(22-6,1)=0.d0
                        c_array(num)%ni(22-6,2)=2.d0*c/ac2
                        c_array(num)%ni(22-6,3)=sqrt(3.d0)*a/ac2
                        ! R17
                        c_array(num)%ni(23-6,1)=-sqrt(3.d0)*c/ac2
                        c_array(num)%ni(23-6,2)=c/ac2
                        c_array(num)%ni(23-6,3)=sqrt(3.d0)*a/ac2
                        ! R18
                        c_array(num)%ni(24-6,1)=-sqrt(3.d0)*c/ac2
                        c_array(num)%ni(24-6,2)=-c/ac2
                        c_array(num)%ni(24-6,3)=sqrt(3.d0)*a/ac2
c 
c                         c Pyramidal <c+a> Slip Systems {1 1 -2 2}<1 1 -2 3>
c                         c R19
c                         c_array(num)%ni(25,1)=c/(2*ac1)
c                         c_array(num)%ni(25,2)=-sqrt(3.d0)*c/(2*ac1)
c                         c_array(num)%ni(25,3)=a/ac1
c                         c R20
c                         c_array(num)%ni(26,1)=c/ac1
c                         c_array(num)%ni(26,2)=0
c                         c_array(num)%ni(26,3)=a/ac1
c                         c R21
c                         c_array(num)%ni(27,1)=c/(2*ac1)
c                         c_array(num)%ni(27,2)=sqrt(3.d0)*c/(2*ac1)
c                         c_array(num)%ni(27,3)=a/ac1
c                         c R22
c                         c_array(num)%ni(28,1)=-c/(2*ac1)
c                         c_array(num)%ni(28,2)=sqrt(3.d0)*c/(2*ac1)
c                         c_array(num)%ni(28,3)=a/ac1
c                         c R23
c                         c_array(num)%ni(29,1)=-c/ac1
c                         c_array(num)%ni(29,2)=0
c                         c_array(num)%ni(29,3)=a/ac1
c                         c R24
c                         c_array(num)%ni(30,1)=-c/(2*ac1)
c                         c_array(num)%ni(30,2)=-sqrt(3.d0)*c/(2*ac1)
c                         c_array(num)%ni(30,3)=a/ac1
                  else
                        write (out,*) "Error: invalid slip type."
                        call die_gracefully
                  end if
c
c ! Allocate number of hardening variables
c
c   Tang_calc: variable to denote Jacobian calculation
c   0 for user supplied, no checking
c   1 for real finite difference (mm10_solveB only)
c   2 for complex finite difference
c   3 for checking with real finite difference
c   4 for checking with complex finite difference
c
c   Note: TJT 10/4/16: tang_calc is now specified through user-input
c                      default value if not specified is 0
c
c ***** START: Add new Constitutive Models into this block *****
                  if (c_array(num)%h_type .eq. 1) then !Voche
                          c_array(num)%num_hard = 1
c                          c_array(num)%tang_calc = 0
                  elseif (c_array(num)%h_type .eq. 2) then !MTS
                          c_array(num)%num_hard = 1
c                          c_array(num)%tang_calc = 0
                  elseif (c_array(num)%h_type .eq. 3) then !User
                          c_array(num)%num_hard = 0
                          c_array(num)%tang_calc = 1
                  elseif (c_array(num)%h_type .eq. 4) then !ORNL
                          c_array(num)%num_hard = c_array(num)%nslip
c                          c_array(num)%tang_calc = 0
                  elseif (c_array(num)%h_type .eq. 7) then !MRR
                          c_array(num)%num_hard = c_array(num)%nslip
c                          c_array(num)%tang_calc = 0
                  elseif (c_array(num)%h_type .eq. 8) then !Armstrong-Frederick
                          c_array(num)%num_hard = c_array(num)%nslip
c                          c_array(num)%tang_calc = 0
                  elseif (c_array(num)%h_type .eq. 9) then !DJGM
                          c_array(num)%num_hard = c_array(num)%nslip
c                          c_array(num)%tang_calc = 0
                  else
                     write(*,101) c_array(num)%h_type
 101  format(
     &      10x,'>> Error: unknown hardening type ', 'i6', '.',
     &    /,10x, 'Aborting...')
      call die_gracefully
                  endif
c ****** END: Add new Constitutive Models into this block ******
c
c                
                  e = c_array(num)%e
                  v = c_array(num)%nu
                  u = c_array(num)%mu
                  if (c_array(num)%elastic_type .eq. 1) then
                        c_array(num)%elast_flex(1,1)=1/e
                        c_array(num)%elast_flex(1,2)=-v/e
                        c_array(num)%elast_flex(1,3)=-v/e
                        c_array(num)%elast_flex(1,4)=0
                        c_array(num)%elast_flex(1,5)=0
                        c_array(num)%elast_flex(1,6)=0
                        c_array(num)%elast_flex(2,1)=-v/e
                        c_array(num)%elast_flex(2,2)=1/e
                        c_array(num)%elast_flex(2,3)=-v/e
                        c_array(num)%elast_flex(2,4)=0
                        c_array(num)%elast_flex(2,5)=0
                        c_array(num)%elast_flex(2,6)=0
                        c_array(num)%elast_flex(3,1)=-v/e
                        c_array(num)%elast_flex(3,2)=-v/e
                        c_array(num)%elast_flex(3,3)=1/e
                        c_array(num)%elast_flex(3,4)=0
                        c_array(num)%elast_flex(3,5)=0
                        c_array(num)%elast_flex(3,6)=0
                        c_array(num)%elast_flex(4,1)=0
                        c_array(num)%elast_flex(4,2)=0
                        c_array(num)%elast_flex(4,3)=0
                        c_array(num)%elast_flex(4,4)=2*(1+v)/e
                        c_array(num)%elast_flex(4,5)=0
                        c_array(num)%elast_flex(4,6)=0
                        c_array(num)%elast_flex(5,1)=0
                        c_array(num)%elast_flex(5,2)=0
                        c_array(num)%elast_flex(5,3)=0
                        c_array(num)%elast_flex(5,4)=0
                        c_array(num)%elast_flex(5,5)=2*(1+v)/e
                        c_array(num)%elast_flex(5,6)=0
                        c_array(num)%elast_flex(6,1)=0
                        c_array(num)%elast_flex(6,2)=0
                        c_array(num)%elast_flex(6,3)=0
                        c_array(num)%elast_flex(6,4)=0
                        c_array(num)%elast_flex(6,5)=0
                        c_array(num)%elast_flex(6,6)=2*(1+v)/e
                  elseif (c_array(num)%elastic_type .eq. 2) then
                        c_array(num)%elast_flex(1,1)=1/e
                        c_array(num)%elast_flex(1,2)=-v/e
                        c_array(num)%elast_flex(1,3)=-v/e
                        c_array(num)%elast_flex(1,4)=0
                        c_array(num)%elast_flex(1,5)=0
                        c_array(num)%elast_flex(1,6)=0
                        c_array(num)%elast_flex(2,1)=-v/e
                        c_array(num)%elast_flex(2,2)=1/e
                        c_array(num)%elast_flex(2,3)=-v/e
                        c_array(num)%elast_flex(2,4)=0
                        c_array(num)%elast_flex(2,5)=0
                        c_array(num)%elast_flex(2,6)=0
                        c_array(num)%elast_flex(3,1)=-v/e
                        c_array(num)%elast_flex(3,2)=-v/e
                        c_array(num)%elast_flex(3,3)=1/e
                        c_array(num)%elast_flex(3,4)=0
                        c_array(num)%elast_flex(3,5)=0
                        c_array(num)%elast_flex(3,6)=0
                        c_array(num)%elast_flex(4,1)=0
                        c_array(num)%elast_flex(4,2)=0
                        c_array(num)%elast_flex(4,3)=0
                        c_array(num)%elast_flex(4,4)=1/u
                        c_array(num)%elast_flex(4,5)=0
                        c_array(num)%elast_flex(4,6)=0
                        c_array(num)%elast_flex(5,1)=0
                        c_array(num)%elast_flex(5,2)=0
                        c_array(num)%elast_flex(5,3)=0
                        c_array(num)%elast_flex(5,4)=0
                        c_array(num)%elast_flex(5,5)=1/u
                        c_array(num)%elast_flex(5,6)=0
                        c_array(num)%elast_flex(6,1)=0
                        c_array(num)%elast_flex(6,2)=0
                        c_array(num)%elast_flex(6,3)=0
                        c_array(num)%elast_flex(6,4)=0
                        c_array(num)%elast_flex(6,5)=0
                        c_array(num)%elast_flex(6,6)=1/u
                  elseif (c_array(num)%elastic_type .eq. 3) then
                        ! do nothing here
                  else
                        write (out,*) "Error: invalid elasticity type."
                        call die_gracefully
                  end if
c                       Generate (one time) the inverse of the flexibility
                  if (c_array(num)%elastic_type .ne. 3) then
                  c_array(num)%elast_stiff = c_array(num)%elast_flex
                  call mm10_invsym(c_array(num)%elast_stiff,6)
                  else
                        C11 = c_array(num)%C11
                        C13 = c_array(num)%C13
                        C44 = c_array(num)%C44
                        C55 = c_array(num)%C55

                        C33 = c_array(num)%C33
                        C12 = c_array(num)%C12
                        C22 = C11
                        C23 = C13
                        C66 = C55               
                        c_array(num)%elast_stiff(1,1)=C11
                        c_array(num)%elast_stiff(1,2)=C12
                        c_array(num)%elast_stiff(1,3)=C13
                        c_array(num)%elast_stiff(1,4)=0.d0
                        c_array(num)%elast_stiff(1,5)=0.d0
                        c_array(num)%elast_stiff(1,6)=0.d0
                        c_array(num)%elast_stiff(2,1)=C12
                        c_array(num)%elast_stiff(2,2)=C22
                        c_array(num)%elast_stiff(2,3)=C23
                        c_array(num)%elast_stiff(2,4)=0.d0
                        c_array(num)%elast_stiff(2,5)=0.d0
                        c_array(num)%elast_stiff(2,6)=0.d0
                        c_array(num)%elast_stiff(3,1)=C13
                        c_array(num)%elast_stiff(3,2)=C23
                        c_array(num)%elast_stiff(3,3)=C33
                        c_array(num)%elast_stiff(3,4)=0.d0
                        c_array(num)%elast_stiff(3,5)=0.d0
                        c_array(num)%elast_stiff(3,6)=0.d0
                        c_array(num)%elast_stiff(4,1)=0.d0
                        c_array(num)%elast_stiff(4,2)=0.d0
                        c_array(num)%elast_stiff(4,3)=0.d0
                        c_array(num)%elast_stiff(4,4)=C44
                        c_array(num)%elast_stiff(4,5)=0.d0
                        c_array(num)%elast_stiff(4,6)=0.d0
                        c_array(num)%elast_stiff(5,1)=0.d0
                        c_array(num)%elast_stiff(5,2)=0.d0
                        c_array(num)%elast_stiff(5,3)=0.d0
                        c_array(num)%elast_stiff(5,4)=0.d0
                        c_array(num)%elast_stiff(5,5)=C55
                        c_array(num)%elast_stiff(5,6)=0.d0
                        c_array(num)%elast_stiff(6,1)=0.d0
                        c_array(num)%elast_stiff(6,2)=0.d0
                        c_array(num)%elast_stiff(6,3)=0.d0
                        c_array(num)%elast_stiff(6,4)=0.d0
                        c_array(num)%elast_stiff(6,5)=0.d0
                        c_array(num)%elast_stiff(6,6)=C66 
c                       write(*,*) c_array(num)%elast_stiff(1:6,1:6)
                    endif
c
                  c_array(num)%valid = .true.
            end subroutine
c
c                 Debug routine, dump the definition to STDOUT
            subroutine print_crystal(num)
                  integer, intent(in) :: num
                  if (num .gt. max_crystals) then
                        write (*,*)
     &                   "Error: crystal number exceeds max crystals."
                        call die_gracefully
                  end if                  
                  if (num .lt. 1) then
                        write (*,*)
     &                   "Error: invalid crystal number."
                        call die_gracefully
                  end if
                  write (*,*) "elip, elastic"
                  write (*,*) c_array(num)%slip_type,
     &                        c_array(num)%elastic_type
                  write (*,*) "e nu"
                  write (*,*) c_array(num)%e,c_array(num)%nu,
     &                        c_array(num)%mu
                  write (*,*) "harden, tau_a"
                  write (*,*) c_array(num)%harden_n,c_array(num)%tau_a,
     &                        c_array(num)%tau_hat_v,c_array(num)%g_o_v,
     &                        c_array(num)%tau_hat_y,c_array(num)%g_o_y,
     &                        c_array(num)%b,c_array(num)%p_v,
     &                        c_array(num)%q_v,c_array(num)%p_y,
     &                        c_array(num)%q_y,c_array(num)%boltz,
     &                        c_array(num)%eps_dot_o_y,
     &                        c_array(num)%eps_dot_o_v,
     &                        c_array(num)%mu_o,c_array(num)%t_o,
     &                        c_array(num)%D_o,
     &                        c_array(num)%theta_o,
     &                        c_array(num)%k_o,
     &                        c_array(num)%tau_y,
     &                        c_array(num)%tau_v,
     &                        c_array(num)%iD_V,
     &                        c_array(num)%voche_m
                  write (*,*) "user"
                  write (*,*) c_array(num)%u1,c_array(num)%u2,
     &                        c_array(num)%u3,c_array(num)%u4,
     &                        c_array(num)%u5,c_array(num)%u6,
     &                        c_array(num)%u7,c_array(num)%u8,
     &                        c_array(num)%u9,c_array(num)%u10
                  write (*,*) "Cxx"
                  write (*,*) c_array(num)%C11,c_array(num)%C12,
     &                        c_array(num)%C13,c_array(num)%C33,
     &                        c_array(num)%C44,c_array(num)%C55
                  write (*,*) "elast_stiff"
                  write (*,*) c_array(num)%elast_stiff
                  write (*,*) "flex"
                  write (*,*) c_array(num)%elast_flex
                  write (*,*) "slip sys"
                  write (*,*) c_array(num)%nslip
                  write (*,*) c_array(num)%ni
                  write (*,*) c_array(num)%bi
                  write (*,*) "htype"
                  write (*,*) c_array(num)%h_type
                  write (*,*) c_array(num)%num_hard
                  write (*,*) "solvertol"
                  write (*,*) c_array(num)%atol,c_array(num)%atol1,
     &                        c_array(num)%rtol,c_array(num)%rtol1,
     &                        c_array(num)%xtol,c_array(num)%xtol1
                  write (*,*) "st_it"
                  write (*,*) c_array(num)%st_it
                  write (*,*) "solver"
                  write (*,*) c_array(num)%solver
                  write (*,*) c_array(num)%strategy
                  write (*,*) c_array(num)%gpall
                  write (*,*) c_array(num)%alter_mode
                  write (*,*) c_array(num)%valid
                  write (*,*) "solver2"
                  write (*,*) c_array(num)%tang_calc,c_array(num)%gpp,
     &                        c_array(num)%method,c_array(num)%miter
                  write (*,*) "user2"
                  write (*,*) c_array(num)%cp_001,c_array(num)%cp_002,
     &                        c_array(num)%cp_003,c_array(num)%cp_004



c
            end subroutine
      end module crystal_data
c
c     ****************************************************************
c     *                                                              *
c     *                       subroutine read_simple_angles          *
c     *                                                              *
c     *                       written by: mcm                        *
c     *                       last modified: 3/10/14                 *
c     *                                                              *
c     *     Read angles to file for the damaged interface material   *
c     *                                                              *
c     ****************************************************************
c
      subroutine read_simple_angles
      use main_data
      use crystal_data, only: simple_angles, nangles, srequired,
     &      mc_array
      use global_data ! old common.main
      implicit integer (a-z)
c
      integer :: inter_mat, nlines, avail_device
      integer, external :: warp3d_get_device_number
      character(len=24) :: filen
      double precision :: t1, t2, t3
      real, dimension(max_crystals) :: rand_reals
c
c           Check to see if the simple angles are required
c
      srequired = .false.
      do el=1,noelem
            inter_mat = elstor(11,el)
            if (inter_mat .ne. -1) then
                  srequired = .true.
                  filen = smatprp(112,inter_mat) 
                  exit
            end if
      end do
c
      if (.not. srequired) return
c
      avail_device = warp3d_get_device_number()
      if (avail_device .eq. -1) then
         write(*,*) "No device available."
         call die_gracefully
      end if
c
c     Count the number of lines
      nlines = 0
      open(avail_device, file=filen)
      do
        read(avail_device,*,END=100)
        nlines = nlines + 1
      end do
 100  close(avail_device)
c
c     Read in all the data
      nangles = nlines
      allocate(simple_angles(nangles,3))
c
      i = 1
      open(avail_device, file=filen)
      do
        read(avail_device, *, end=200) t1, t2, t3
        simple_angles(i,1) = t1
        simple_angles(i,2) = t2
        simple_angles(i,3) = t3
        i = i + 1
      end do
 200  close(avail_device)
c
c
c           This is an overly-expensive way of doing this, but oh well
c           Setup mc_array
      allocate(mc_array(noelem, max_crystals))
      call init_random_seed()
      do i=1,noelem
        call random_number(rand_reals)
        mc_array(i,1:max_crystals) = int(rand_reals*nangles)+1
      end do
c
      call wmpi_send_simple_angles
c
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine read_crystal_data            *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 4/17/12 mcm                *
c     *                                                              *
c     *      read in crystal data from file to in memory structures  *
c     *                                                              *
c     ****************************************************************
c
      subroutine read_crystal_data
      use crystal_data, only : angle_input, crystal_input,
     &            data_offset, defined_crystal
      use main_data
      use global_data ! old common.main
      implicit integer (a-z)
c
      integer :: mxcry, nelem, ncry matnum, el, i, iodev, p_matnum
      integer :: avail_device
      integer, external :: warp3d_get_device_number
      logical :: countme, crystalsl, anglesl, open_file, send_mess
      character(len=24) :: filen
c     
c     Just skip if we don't have a crystal
c
      if (.not. defined_crystal) then
            return
      end if
c
c
c     Forward pass: get numbers to allocate and offsets
c
      allocate(data_offset(noelem))
      data_offset = 0
      mxcry = 0
      nelem = 0
      do el=1,noelem
            matnum = elstor(2,el)
            mattype = matprp(9,matnum)
            if (mattype .eq. 10) then
                  countme = (imatprp(104,matnum) .eq. 2) .or.
     &                  (imatprp(107,matnum) .eq. 2)
                  if (countme) then
                        nelem = nelem +1
                        data_offset(el) = nelem
                        ncry = imatprp(101,matnum)
                        if (ncry .gt. mxcry) mxcry = ncry
                  end if
            end if
      end do
c      
c     Backward pass: actual data
c
      send_mess = .false.
      if( any(data_offset .gt. 0) ) then
            write(out,*) '' 
            write(out,'(A)') '>> Reading crystal definitions ...'
            send_mess = .true.
            allocate(angle_input(nelem,mxcry,3))
            allocate(crystal_input(nelem,mxcry))
            p_matnum = -1
            iodev = warp3d_get_device_number()
            if (iodev .eq. -1) then
              write(*,*) "No available device."
              call die_gracefully
            end if
            open_file = .false.
            do el=1,noelem
                  if (data_offset(el) .ne. 0) then
                        matnum = elstor(2,el)
                        if (matnum .ne. p_matnum) then
                          if (el .ne. 1) then
                              close(iodev)
                          end if
                          filen = smatprp(112,matnum)
                          open(iodev,FILE=filen,READONLY)
                          open_file = .true.
                        end if
                        ncry = imatprp(101,matnum)
c                        filen = smatprp(112,matnum)
c                        open(iodev,FILE=filen,READONLY)
                        if (imatprp(104,matnum) .eq. 2) then
                              crystalsl = .true.
                        else
                              crystalsl = .false.
                        end if
                        if (imatprp(107,matnum) .eq. 2) then
                              anglesl = .true.
                        else
                              anglesl = .false.
                        end if
                        call read_defs(iodev,el,ncry,anglesl,
     &                        crystalsl, angle_input(el,1:ncry,1:3),
     &                        crystal_input(el,1:ncry),out)
c                        close(iodev)
                        p_matnum = matnum
                  end if
            end do
            if (open_file) then
                  close(iodev)
            end if
      end if
c
c           If we're using MPI, send out the arrays (including the
c           general crystal struct).  If not it's a dummy routine
c
      call wmpi_send_crystals

      if( send_mess) then
         write(out,*) '    ... Done'
         write(out,*) 
      end if
c     
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine read_defs                       *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 03/23/2012                 *
c     *                                                              *
c     *     Helper which reads crystal numbers and orientations from *
c     *     a flat file                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine read_defs(ionum,elnum,ncry,angles,crystals,results_ang,
     &            results_cry,out)
            integer, intent(in) :: ionum, elnum, ncry, out
            logical, intent(in) :: angles, crystals
            double precision, dimension(ncry,3), intent(out) :: 
     &            results_ang
            double precision :: a,b,c
            integer :: n, d, i
            integer, dimension(ncry), intent(out) :: results_cry
c
            character(len=24) :: fmat

c           Scan the file until we find the first entry which matches el
            do
c                 Read into dummy variables just in case
                  if (angles .and. crystals) then
                        read(ionum, *,end=10) n,a,b,c,d
                  elseif (angles) then
                        read(ionum,*,end=10) n,a,b,c
                  elseif (crystals) then
                        read(ionum,*,end=10) n,d
                  else
                        exit
                  endif
c                 Test and input
                  i = 1
                  if (n .eq. elnum) then
                        if (angles .and. crystals) then
                              results_ang(i,1) = a
                              results_ang(i,2) = b
                              results_ang(i,3) = c
                              results_cry(i) = d
                        elseif (angles) then
                              results_ang(i,1) = a
                              results_ang(i,2) = b
                              results_ang(i,3) = c
                        elseif (crystals) then
                              results_cry(i) = d
                        else
                              exit
                        endif
c                             Actual loop to read data
                        do i=2,ncry
                              if (angles .and. crystals) then
                                    read(ionum, *,end=12) n,
     &                                    results_ang(i,1),
     &                                    results_ang(i,2),
     &                                    results_ang(i,3),
     &                                    results_cry(i)
                              elseif (angles) then
                                    read(ionum, *,end=12) n,
     &                                    results_ang(i,1),
     &                                    results_ang(i,2),
     &                                    results_ang(i,3)
                              elseif (crystals) then
                                    read(ionum,*,end=12) n,
     &                                    results_cry(i)
                              else
                                    exit
                              endif
                              if (n .ne. elnum) then
c                                   Terrible error
                                    goto 12
                              end if
                        end do
                        exit
                  end if


            end do


            return
10    continue
            write (out,14) elnum
            call die_gracefully

12    continue
            write (out,13) elnum
            call die_gracefully

13    format(/1x,'>>>>> Parse Error: insufficient data for element ',
     &            i4/)
14    format(/1x,'>>>>> Parse Error: element ', i4, ' not found.' /)

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine cleanup_crystal                 *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 04/18/2012                 *
c     *                                                              *
c     *     Call at the end of everything to clean up allocatable    *
c     *     arrays.                                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine cleanup_crystal
      use crystal_data, only : angle_input, crystal_input,
     &            data_offset, simple_angles, mc_array
      use global_data ! old common.main
      implicit integer (a-z)
c
      if (allocated(angle_input)) deallocate(angle_input)
      if (allocated(crystal_input)) deallocate(crystal_input)
      if (allocated(data_offset)) deallocate(data_offset)
      if (allocated(simple_angles)) deallocate(simple_angles)
      if (allocated(mc_array)) deallocate(mc_array)

      call wmpi_dealloc_crystals

      return

      end subroutine


c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine avg_cry_elast_props          *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 4/5/12 mcm                 *
c     *                                                              *
c     *      Set average elastic properties.  For elastic materials  *
c     *           this will be exact, for anisotropic materials we   *
c     *           need to think about it more                        *
c     *                                                              *
c     ****************************************************************
c
      subroutine avg_cry_elast_props
      use crystal_data, only : angle_input, crystal_input,
     &            data_offset, c_array, defined_crystal
      use main_data
      use global_data !  old common.main
      implicit integer (a-z)
c
      integer :: i, j, k, ncrystals, cnum, osn, num, ecount
      real :: e_avg, nu_avg
c
c
c     Just skip if we don't have a crystal
c
      if (.not. defined_crystal) then
            return
      end if
c
c     Run through our materials, find the CP materials, average their elastic
c     constants
      do i=1,nummat
      if (matprp(9,i) .eq. 10) then
        matprp(1,i) = 0.0
        matprp(2,i) = 0.0
        ecount = 0
        do j = 1, noelem
          if (iprops(38,j)  .eq. i) then
            ecount = ecount + 1
            e_avg = 0.0
            nu_avg = 0.0
            ncrystals = imatprp(101, i)
            do k = 1, ncrystals
c             Get the local crystal number
              if (imatprp(104,i) .eq. 1) then
                cnum = imatprp(105,i)
              elseif (imatprp(104,i) .eq. 2) then
                osn = data_offset(j)
                cnum = crystal_input(osn,k)
c               Couldn't do this earlier, so check here
                if ((cnum .gt. max_crystals) .or. 
     &                (cnum .lt. 0)) then
                  write (out,'("Crystal ", i3, " not valid")')
     &                 cnum
                  call die_gracefully
                 end if
              else
                write(out,9502) 
                call die_gracefully
              end if
c              
c             INSERT AVERAGING
c   
              e_avg = e_avg + SNGL(c_array(cnum)%e)
              nu_avg = nu_avg + SNGL(c_array(cnum)%nu)
            end do
            e_avg = e_avg / SNGL(ncrystals)
            nu_avg = nu_avg / SNGL(ncrystals)
            props(7,j) = e_avg
            props(8,j) = nu_avg
            matprp(1,i) = matprp(1,i) + e_avg
            matprp(2,i) = matprp(2,i) + nu_avg
          end if
        end do
        matprp(1,i) = matprp(1,i) / SNGL(ecount)
        matprp(2,i) = matprp(2,i) / SNGL(ecount)
      end if
      end do
c
c     I maintain this is not my fault, we now need to go fix the element props
c
      do i=1,noelem
            if (iprops(25, i) .eq. 10) then
                  num = iprops(38,i)
                  props(7,i) = matprp(1,num)
                  props(8,i) = matprp(2,num)
            end if
      end do

      return

 9502 format(/,1x,
     & '>>>> System error: unexpected input type in avg_elast_props!',
     &       ' Aborting.'/)

      end subroutine
c
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine init_random_see              *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 3/24/14                    *
c     *                                                              *
c     *      Setup the random seed for number generation             *
c     *      This function copied from the GNU documentation         *
c     *                                                              *
c     ****************************************************************
c
      subroutine init_random_seed()
      use ifport, only: getpid
      implicit none
      integer, allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8), pid, t(2), s
      integer(8) :: count, tms
          
      call random_seed(size = n)
      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream", 
     &   form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
         read(un) seed
         close(un)
      else
         ! Fallback to XOR:ing the current time and pid. The PID is
         ! useful in case one launches multiple instances of the same
         ! program in parallel.
         call system_clock(count)
         if (count /= 0) then
            t = transfer(count, t)
         else
            call date_and_time(values=dt)
            tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 
     &           + dt(2) * 31_8 * 24 * 60 * 60 * 1000
     &           + dt(3) * 24 * 60 * 60 * 60 * 1000
     &           + dt(5) * 60 * 60 * 1000 
     &           + dt(6) * 60 * 1000 + dt(7) * 1000
     &           + dt(8)
            t = transfer(tms, t)
         end if
         s = ieor(t(1), t(2))
         pid = getpid() + 1099279 ! Add a prime
         s = ieor(s, pid)
         if (n >= 3) then
            seed(1) = t(1) + 36269
            seed(2) = t(2) + 72551
            seed(3) = pid
            if (n > 3) then
               seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
            end if
         else
            seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
         end if
      end if
      call random_seed(put=seed)
      end subroutine init_random_seed
