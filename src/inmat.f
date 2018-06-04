c     ****************************************************************
c     *                                                              *
c     *                      subroutine inmat                        *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 4/9/2018 rhd               *
C     *                                                              *
C     *     input of properties of the materials in the material     *
c     *     library for the current problem.                         *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine inmat( sbflg1, sbflg2, matnum )
      use global_data ! old common.main
      use main_data, only : matprp, lmtprp, imatprp, dmatprp,
     &  smatprp, nonlocal_analysis, umat_used, creep_model_used
      use erflgs
c
      implicit none
      logical :: sbflg1, sbflg2
      integer :: matnum
c
c
c                       local declarations
c
      integer :: dum, nc, matn, open, param, openxt, premat,
     &           curve_set_no, yl_prop, yl_start, um_prop,
     &           um_start, status
      character :: name*80, mname*24, dums*10
      logical :: matchs, numr, endcrd, true, label, scanms, numd,
     &           matchs_exact
      logical, parameter :: local_debug = .false.
      real ::  dumr, kx_stiff, ky_stiff, kz_stiff, link_mass,
     &        link_stiff
      real, parameter :: fgm_mark = -99.0
      double precision :: dumd
      double precision, parameter :: one=1.0d0
c
      if( sbflg1 ) go to 205
c
c                       if the element properties have already been
c                       defined, then print a warning that this mat-
c                       erial will not be used in the current problem.
c
      if( elprop ) call errmsg(116,dum,dums,dumr,dumd)
c
c
c **********************************************************************
c *                                                                    *
c *                     get the name of the material being input.      *
c *                     each material must have a name for id          *
c *                     purposes. if there is no name, print an        *
c *                     error message and read another high level      *
c *                     command.                                       *
c *                                                                    *
c **********************************************************************
c
c
      if( label(dum) ) then
         name  = ' '
         mname = ' '
         call entits( name, nc )
         if( nc .gt. 24 ) nc = 24
         mname(1:nc) = name(1:nc)
c
c                       search for the specified material name in the
c                       existing material library. if it is found,
c                       print error message and return to the driver
c                       subroutine because an existing material can-
c                       not be overridden, only deleted and then re-
c                       defined.
c
         matn = mathed/two16
c
 201     if( matn.eq.32460 ) go to 202
            if( scanms(matnam(matn),mname,24) ) then
               call errmsg(115,matn,dums,dumr,dumd)
               go to 9998
            end if
            matn = matlst(matn)/two16
         go to 201
c
c                       not in library. check to make sure library is
c                       not full.
c
 202     nummat = nummat + 1
         if( nummat.gt.mxmat ) then
            nummat = mxmat
            call errmsg(25,dum,dums,dumr,dumd)
            go to 9998
         end if
c
c                       find the first open slot in the material lib-
c                       rary vector matnam and assign it to the spec-
c                       ified material.
c
         matn = mathed/two16
         open = mathed-matn*two16
c
c                       make sure there isn't an error in the slot
c                       available.
c
         if( open.le.0.or.open.gt.mxmat ) then
            param = 1
            call errmsg(114,param,dums,dumr,dumd)
            go to 9998
         end if
c
         matnum         = open
         matnam(matnum) = mname
c
c                       update the list of open and occupied slots
c                       in matnam. note that if a slot is occupied,
c                       it cannot also be open. therefore, if the
c                       slot is occupied, set the open linked list
c                       to zero for that slot, and vice-versa.
c
c                       find the next open slot
c
         openxt = matlst(open)-(matlst(open)/two16)*two16
c
c                       check for the condition that the head of the
c                       occupied list is also the end. if so, modify
c                       both lists. in either case, modify the heads
c                       of both lists accordingly.
c
         if( matn.eq.32460 ) then
            mathed       = open*two16+openxt
            matlst(open) = matn*two16
            go to 204
         else
            mathed       = matn*two16+openxt
         end if
c
c                       the head of the occupied linked list is not
c                       the end. find the end and modify both lists.
c
 203     premat = matn
         matn   = matlst(premat)/two16
         if( matn.eq.32460 ) then
            matlst(premat) = open*two16
            matlst(open)   = matn*two16
         else
            go to 203
         end if
c
c                       assign default material values. the current
c                       ordering of material values is:
c
c                  (*)  1  -- young's modulus
c                  (*)  2  -- poisson's ratio
c                       3  -- kinematic hardening ratio (beta)
c                       4  -- tangent modulus for bilinear strain
c                             hardening (et)
c                       5  -- inviscid yield stress (sigma-o)
c                  (*)  6  -- thermal expansion coefficient (isotropic)
c                  (*)  7  -- mass density
c                       8  -- linear elastic material flag
c                       9  -- material model type
c                              = 1 vectorized linear elastic
c                                  and invisicid plasticity model.
c                                  isotropic/kinematic hardening
c                              = 2 nonlinear elastic with linear +
c                                  power law. small strains only.
c                                  rate independent
c                              = 3 general gurson/mises model including
c                                  nucleation, linear hardening,
c                                  power-law hardening, matrix
c                                  viscoplasticity
c                              = 4 Cohesive zone models: linear elastic,
c                                  bilinear, ramp, exponential_1 and
c                                  exponential_2, ppr, cavitation
c                              = 5 advanced cyclic plasticity model
c                                  developed by kristine cochran
c                              = 6 creep model
c                              = 7 advanced mises model with hydrogen
c                                  effects developed by yuiming liang
c                              = 8 general UMAT for warp3d.
c                              = 10 (matching back up with file #s)
c                                   CP model by mark messner
c                              = 11 interface damage model
c                      10  -- viscoplastic m power
c                      11  -- power-law hardening n power
c                      12  -- viscoplastic reference strain rate
c                      13  -- debug material model computations
c                      14  -- initial porosity (f sub 0)
c                      15  -- gurson model parameter q1
c                      16  -- gurson model parameter q2
c                      17  -- gurson model parameter q3
c                      18  -- nucleation flag (true or false)
c                      19  -- nucleation parameter sn
c                      20  -- nucleation parameter en
c                      21  -- nucleation parameter fn
c                      22  -- allow material model to cut step due to
c                      23  -- flag to allow crack growth element killing
c                             excessive reversed plasticity
c                      24  -- segmental stress-strain curve logical flag
c                      25  -- flag to indicate anisotropic thermal
c                             coefficients are defined
c                      26-31  anisotropic thermal expansion coefficients
c                      32  -- interface stiffness in the longitudinal direction
c                      33  -- interface stiffness in the transverse direction
c                      34  -- interface stiffness in the normal direction
c                      35  -- critical normal stress of the interface
c                      36  -- critical shear stress of the interface
c                      37  -- shape parameter for bilinear and ramp
c                      38  -- second ( additional) shape parameter for ramp
c                      39  -- critical separation distance in sliding
c                      40  -- critical separation distance in opening
c                      41  -- equivalent critical separation distance
c                      42  -- a ratio to determine the equivalent separation
c                             under mixed mode loading ( =0 => mode I )
c                      43  -- a flag for identifying the interface element
c                      44  -- type of interface models
c                             1 - linear elastic; 2- bilinear
c                             3 - ramp; 4 - exponential_1; 5 - exponential_2
c                      45  -- for segmental curve model, the segmental curve
c                             set number
c                  (*) 46  -- ductile material volume fracture for fgm
c                      47  -- = 0 homogeneous cohesive material
c                             = 1 functionally graded cohesive material
c                      48  -- critical separation distance ductile (fgm)
c                      49  -- critical separation distance brittle (fgm)
c                      50  -- critical stress ductile (fgm)
c                      51  -- critical stress brittle (fgm)
c                      52  -- beta_ductile (fgm)
c                      53  -- beta_brittle (fgm)
c                      54  -- compression stiffness multiplier for
c                             cohesive materials
c                      55  -- start of props for cyclic plasticity model
c                              see inmat_cyclic
c                      70 -- start of props for yuemin liang's adv.
c                             mises model that includes effects of
c                             staturated hydrogen
c                              70 - yl_1
c                              71 - yl_2
c                              72 - yl_3
c                              73 - yl_4
c                              74 - yl_5
c                              75 - yl_6
c                              76 - yl_7
c                              77 - yl_8
c                              78 - yl_9
c                              79 - yl_10
c
c                      80-89  creep model
c
c                      90-- PPR cohesive model
c                              90-98 currently used. see comments
c                              in inmat_cohesive routine
c                      100-- local tolerance for CP NR loop
c                      101-- number of crystals at e. material point (or max n)
c                      102-- angle convention (1=Kocks)
c                      103-- angle type (1=degree, 2=radians)
c                      104-- crystal input (1=single, 2=file)
c                      105-- crystal number (for single)
c                      106-- crystal offset (for list)
c                      107-- orientation input (1=single, 2=file)
c                      108-110-- psi, theta, phi (for single)
c                      111-- orientation offset (for list)
c                      112-- STRING crystal list (for offset/list)
c                      113-- STRING orientation list (for offset/list)
c
c                      115-- macroscale material model number
c                      116-118-- s vector
c                      119-121-- l vector
c                      122-125-- t vector
c                      126-- l_s
c                      127-- l_l
c                      128-- l_t
c                      129-- alpha_dmg
c                      130-- nstacks (temp, should calculate from element sz)
c                      131-- nfail ("")
c                       132-- macro_sz
c                       133-- cp_sz
c
c                      148-- link2 x-stiffness
c                      149-- link2 y-stiffness
c                      150-- link2 z-stiffness
c
c                      151-200 -- Abaqus compatible UMAT
c                              151 - um_1
c                              151 - um_2
c                              ...
c                              ...
c                              200 - um_50
c
c                      201-230 cavity option for cohesive material
c
c                      the beta_fact is used to assist in construction
c                      of planar models containing one layer of elements.
c                      the stiffness and internal forces are mutiplied
c                      by the beta_fact to simulate the effect of a
c                      non-unit thickness.
c
c                  (*) some material values maybe specified as having the
c                      value 'fgm' (a string). In such cases the user
c                      supplied nodal values for the model are interpolated
c                      at the element gauss points during execution.
c
c
 204     continue
         beta_fact = one
         matprp(1,matnum)  = 30000.0
         matprp(2,matnum)  = 0.3
         matprp(3,matnum)  = 1.0
         matprp(4,matnum)  = 0.0
         matprp(5,matnum)  = 0.0
         matprp(6,matnum)  = 0.0
         matprp(7,matnum)  = 0.0
         lmtprp(8,matnum)  = .false.
         matprp(9,matnum)  = 1
         matprp(10,matnum) = 0.0
         matprp(11,matnum) = 0.0
         matprp(12,matnum) = 0.0
         lmtprp(13,matnum) = .false.
         matprp(14,matnum) = 0.0
         matprp(15,matnum) = 1.5
         matprp(16,matnum) = 1.0
         matprp(17,matnum) = 2.25
         lmtprp(18,matnum) = .false.
         matprp(19,matnum) = 0.1
         matprp(20,matnum) = 0.3
         matprp(21,matnum) = 0.04
         lmtprp(22,matnum) = .true.
         lmtprp(23,matnum) = .false.
         lmtprp(24,matnum) = .false.
         lmtprp(25,matnum) = .false.
         matprp(26,matnum) = 0.0
         matprp(27,matnum) = 0.0
         matprp(28,matnum) = 0.0
         matprp(29,matnum) = 0.0
         matprp(30,matnum) = 0.0
         matprp(31,matnum) = 0.0
         matprp(32,matnum) = 300000.0
         matprp(33,matnum) = 300000.0
         matprp(34,matnum) = 300000.0
         matprp(35,matnum) = 30000.0
         matprp(36,matnum) = 11538.0
         matprp(37,matnum) = .15
         matprp(38,matnum) = .5
         matprp(39,matnum) = 0.01
         matprp(40,matnum) = 0.01
         matprp(41,matnum) = 0.01
         matprp(42,matnum) = 0.0
         lmtprp(43,matnum) = .false.
         matprp(46,matnum) = 0.0
         matprp(47,matnum) = 0.0
         matprp(48,matnum) = 0.0
         matprp(49,matnum) = 0.0
         matprp(50,matnum) = 0.0
         matprp(51,matnum) = 0.0
         matprp(52,matnum) = 0.0
         matprp(53,matnum) = 0.0
         matprp(54,matnum) = 10.0
         matprp(55:99,matnum) = 0.0
         dmatprp(100,matnum) = 1D-06
         imatprp(101,matnum) = 1
         imatprp(102,matnum) = 1
         imatprp(103,matnum) = 1
         imatprp(104,matnum) = 3
         imatprp(105,matnum) = 0
         imatprp(106,matnum) = 0
         imatprp(107,matnum) = 3
         dmatprp(108:110,matnum) = 0.0d00
         imatprp(111,matnum) = 0
         smatprp(112,matnum) = 'ni'
         smatprp(113,matnum) = 'ni'
         matprp(114:mxmtpr,matnum) = 0.0
         dmatprp(151:200,matnum) = 0.0d00   ! UMAT properties
      else
         call errmsg(3,dum,dums,dumr,dumd)
         go to 9998
      end if
c
c
c **********************************************************************
c *                                                                    *
c *                     input properties of the material               *
c *                                                                    *
c **********************************************************************
c
c
      if( endcrd(dum) ) go to 9997
c
c                       properties on next card(s)
c
 205  continue
      if( matchs('properties',10) ) then
         if( endcrd(dum) ) go to 206
         go to 210
      else
         go to 210
      end if
c
c                       read the name of the property and its value.
c                       if there is no match between what is input
c                       and the property list, then print an error
c                       message and continue scanning for a property
c                       name. if an end of card is encountered, print
c                       another error message and return to read a
c                       high level command.
c
 206  call readsc
c
c                       note the matching scheme requires that
c                       some keywords be checked before others,
c                       e.g., nucleation is checked before nu, alphaxy
c                       before alphax and al. use of matchs_exact
c                       helps reduce this problem.
c
 210  continue
      if ( matchs('eps_n',5)            ) go to 420
      if ( matchs_exact('en')           ) go to 420
      if ( matchs_exact('e_n')          ) go to 420
      if ( matchs_exact('e')            ) go to 220
      if ( matchs('nucleation',5)       ) go to 400
      if ( matchs_exact('nu')           ) go to 230
      if ( matchs_exact('beta')         ) go to 240
      if ( matchs('tan_e',5)            ) go to 250
      if ( matchs('yld_pt',6)           ) go to 260
      if ( matchs('alphaxy',7)          ) go to 270
      if ( matchs('alphayz',7)          ) go to 271
      if ( matchs('alphaxz',7)          ) go to 272
      if ( matchs('alphax',6)           ) go to 273
      if ( matchs('alphay',6)           ) go to 274
      if ( matchs('alphaz',6)           ) go to 275
      if ( matchs('alpha',5)            ) go to 276
      if ( matchs_exact('rho')          ) go to 280
      if ( matchs('thickness_ratio',8)  ) go to 285
      if ( matchs_exact('bilinear')     ) go to 300
      if ( matchs_exact('link')         ) go to 300
      if ( matchs_exact('gurson')       ) go to 310
      if ( matchs_exact('mises')        ) go to 312
      if ( matchs('n_power',5)          ) go to 330
      if ( matchs('m_power',5)          ) go to 320
      if ( matchs('ref_eps',7)          ) go to 340
      if ( matchs('deformation',8)      ) go to 315
      if ( matchs('debug',4)            ) go to 350
      if ( matchs_exact('f_0')          ) go to 360
      if ( matchs_exact('f0')           ) go to 360
      if ( matchs_exact('q1')           ) go to 370
      if ( matchs_exact('q2')           ) go to 380
      if ( matchs_exact('q3')           ) go to 390
      if ( matchs_exact('s_n')          ) go to 410
      if ( matchs_exact('f_n')          ) go to 430
      if ( matchs('no_cutback',6)       ) go to 440
      if ( matchs('killable',4)         ) go to 450
      if ( matchs('curve',5)            ) go to 460
      if ( matchs('cohesive',6)         ) go to 465
      if ( matchs_exact('nonlocal')     ) go to 470
      if ( matchs('sig_o',5)            ) go to 260
      if ( matchs('linear_elastic',6)   ) go to 290
      if ( matchs_exact('creep')        ) go to 600
      if ( matchs_exact('B')            ) go to 605
      if ( matchs_exact('cyclic') )         go to 800
      if ( matchs_exact('cylic') )          go to 800
      if ( matchs_exact('mises_hydrogen') ) go to 1000
        if ( matchs_exact('trap_number') )              go to 1010
        if ( matchs_exact('init_disloc_density') )      go to 1010
        if ( matchs_exact('disloc_density_parameter') ) go to 1010
        if ( matchs_exact('lattice_parameter' ) )       go to 1010
        if ( matchs_exact('binding_energy') )           go to 1010
        if ( matchs_exact('nils_number') )              go to 1010
        if ( matchs_exact('molar_volume') )             go to 1010
        if ( matchs_exact('h_partial_volume') )         go to 1010
        if ( matchs_exact('initial_hydrogen_conc') )    go to 1010
        if ( matchs_exact('softening_parameter') )      go to 1010
      if ( matchs('umat',4) )                           go to 1100
         if( matchs_exact( 'um_1') ) go to 1110
         if( matchs_exact( 'um_2') ) go to 1110
         if( matchs_exact( 'um_3') ) go to 1110
         if( matchs_exact( 'um_4') ) go to 1110
         if( matchs_exact( 'um_5') ) go to 1110
         if( matchs_exact( 'um_6') ) go to 1110
         if( matchs_exact( 'um_7') ) go to 1110
         if( matchs_exact( 'um_8') ) go to 1110
         if( matchs_exact( 'um_9') ) go to 1110
         if( matchs_exact( 'um_10') ) go to 1110
         if( matchs_exact( 'um_11') ) go to 1110
         if( matchs_exact( 'um_12') ) go to 1110
         if( matchs_exact( 'um_13') ) go to 1110
         if( matchs_exact( 'um_14') ) go to 1110
         if( matchs_exact( 'um_15') ) go to 1110
         if( matchs_exact( 'um_16') ) go to 1110
         if( matchs_exact( 'um_17') ) go to 1110
         if( matchs_exact( 'um_18') ) go to 1110
         if( matchs_exact( 'um_19') ) go to 1110
         if( matchs_exact( 'um_20') ) go to 1110
         if( matchs_exact( 'um_21') ) go to 1110
         if( matchs_exact( 'um_22') ) go to 1110
         if( matchs_exact( 'um_23') ) go to 1110
         if( matchs_exact( 'um_24') ) go to 1110
         if( matchs_exact( 'um_25') ) go to 1110
         if( matchs_exact( 'um_26') ) go to 1110
         if( matchs_exact( 'um_27') ) go to 1110
         if( matchs_exact( 'um_28') ) go to 1110
         if( matchs_exact( 'um_29') ) go to 1110
         if( matchs_exact( 'um_30') ) go to 1110
         if( matchs_exact( 'um_31') ) go to 1110
         if( matchs_exact( 'um_32') ) go to 1110
         if( matchs_exact( 'um_33') ) go to 1110
         if( matchs_exact( 'um_34') ) go to 1110
         if( matchs_exact( 'um_35') ) go to 1110
         if( matchs_exact( 'um_36') ) go to 1110
         if( matchs_exact( 'um_37') ) go to 1110
         if( matchs_exact( 'um_38') ) go to 1110
         if( matchs_exact( 'um_39') ) go to 1110
         if( matchs_exact( 'um_40') ) go to 1110
         if( matchs_exact( 'um_41') ) go to 1110
         if( matchs_exact( 'um_42') ) go to 1110
         if( matchs_exact( 'um_43') ) go to 1110
         if( matchs_exact( 'um_44') ) go to 1110
         if( matchs_exact( 'um_45') ) go to 1110
         if( matchs_exact( 'um_46') ) go to 1110
         if( matchs_exact( 'um_47') ) go to 1110
         if( matchs_exact( 'um_48') ) go to 1110
         if( matchs_exact( 'um_49') ) go to 1110
         if( matchs_exact( 'um_50') ) go to 1110
         if( matchs_exact('CP')     ) go to 2000
         if( matchs_exact( 'interface_damage')) go to 2100
      if( matchs_exact('kx_stiff') ) go to 2200
      if( matchs_exact('ky_stiff') ) go to 2205
      if( matchs_exact('kz_stiff') ) go to 2210
      if( matchs_exact('mass_link') ) go to 2215
      if( matchs_exact('stiff_link') ) go to 2220

c
      if ( matchs(',',1) ) go to 215

c
c                       there is no match with existing properties.
c                       check for end of card. if not, print error
c                       message.
c
      if( endcrd(dum) ) then
         go to 9998
      else
         call errmsg(4,dum,dums,dumr,dumd)
         if( true(dum) ) go to 210
      end if
c
c                       comma separator. if at the end of a line,
c                       denotes continuation. otherwise, ignore.
c
 215  continue
      if( endcrd(dum) ) then
         go to 206
      else
         go to 210
      end if
c
c * ********************************************************************
c *                                                                    *
c *                     input modulus of elasticity                    *
c *                                                                    *
c **********************************************************************
c
c
 220  continue
c
      call inmat_fgm( status, 'e'  )
      if ( status  .eq. 0 ) then
          matprp(1,matnum) = fgm_mark
          go to 210
      elseif ( status .eq. 1 ) then
          go to 210
      end if
c
      if( .not. numr(matprp(1,matnum)) ) then
         call errmsg(5,dum,'e',dumr,dumd)
      end if
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                     input poisson's ratio                          *
c *                                                                    *
c **********************************************************************
c
c
 230  continue
      call inmat_fgm( status, 'nu'  )
      if ( status  .eq. 0 ) then
          matprp(2,matnum) = fgm_mark
          go to 210
      elseif ( status .eq. 1 ) then
          go to 210
      end if
c
      if( .not. numr(matprp(2,matnum)) ) then
         call errmsg(5,dum,'nu',dumr,dumd)
      end if
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                     input hardening parameter                      *
c *                                                                    *
c **********************************************************************
c
c
 240  continue
      if( .not. numr(matprp(3,matnum)) ) then
         call errmsg(5,dum,'beta', dumr,dumd)
      end if
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                     input tangent modulus of elasticity            *
c *                                                                    *
c **********************************************************************
c
c
 250  continue
      call inmat_fgm( status, 'tan_e'  )
      if ( status  .eq. 0 ) then
          matprp(4,matnum) = fgm_mark
          go to 210
      elseif ( status .eq. 1 ) then
          go to 210
      end if
c
      if( .not. numr(matprp(4,matnum)) ) then
         call errmsg(5,dum,'etan',dumr,dumd)
      end if
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                     input initial yield point                      *
c *                                                                    *
c **********************************************************************
c
c
 260  continue
      call inmat_fgm( status, 'yld_pt'  )
      if ( status  .eq. 0 ) then
          matprp(5,matnum) = fgm_mark
          go to 210
      elseif ( status .eq. 1 ) then
          go to 210
      end if

      if( .not. numr(matprp(5,matnum)) ) then
         call errmsg(5,dum,'ylpt',dumr,dumd)
      end if
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *      input the coefficients of thermal expansion.                  *
c *                                                                    *
c *  user can specify the single alpha value to request isotropic      *
c *  properties. if any other alpha values (alphax, alphay, ...)       *
c *  are input, the material is marked anistropic.                     *
c *                                                                    *
c **********************************************************************
c
c
 270  continue
      if( .not. numr(matprp(29,matnum)) ) then
         call errmsg(5,dum,'alphaxy',dumr,dumd)
         go to 210
      end if
      lmtprp(25,matnum) = .true.
      go to 210
 271  continue
      if( .not. numr(matprp(30,matnum)) ) then
         call errmsg(5,dum,'alphayz',dumr,dumd)
         go to 210
      end if
      lmtprp(25,matnum) = .true.
      go to 210
 272  continue
      if( .not. numr(matprp(31,matnum)) ) then
         call errmsg(5,dum,'alphaxz',dumr,dumd)
         go to 210
      end if
      lmtprp(25,matnum) = .true.
      go to 210
 273  continue
      if( .not. numr(matprp(26,matnum)) ) then
         call errmsg(5,dum,'alphax',dumr,dumd)
         go to 210
      end if
      lmtprp(25,matnum) = .true.
      go to 210
 274  continue
      if( .not. numr(matprp(27,matnum)) ) then
         call errmsg(5,dum,'alphay',dumr,dumd)
         go to 210
      end if
      lmtprp(25,matnum) = .true.
      go to 210
 275  continue
      if( .not. numr(matprp(28,matnum)) ) then
         call errmsg(5,dum,'alphaz',dumr,dumd)
         go to 210
      end if
      lmtprp(25,matnum) = .true.
      go to 210
 276  continue
      call inmat_fgm( status, 'alpha'  )
      if ( status  .eq. 0 ) then
          matprp(6,matnum) = fgm_mark
          go to 210
      elseif ( status .eq. 1 ) then
          go to 210
      end if
      if( .not. numr(matprp(6,matnum)) ) then
         call errmsg(5,dum,'alpha',dumr,dumd)
         go to 210
      end if
      lmtprp(25,matnum) = .false.
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                     input the mass density                         *
c *                                                                    *
c **********************************************************************
c
c
 280  continue
      call inmat_fgm( status, 'rho'  )
      if ( status  .eq. 0 ) then
          matprp(7,matnum) = fgm_mark
          go to 210
      elseif ( status .eq. 1 ) then
          go to 210
      end if
c
      if(.not.numr(matprp(7,matnum))) then
         call errmsg(5,dum,'rho',dumr,dumd)
      end if
      go to 210

c
c
c **********************************************************************
c *                                                                    *
c *                     input the thickness ratio for 2d work          *
c *                                                                    *
c **********************************************************************
c
c
 285  continue
      if(.not.numd(beta_fact)) then
         call errmsg(5,dum,'thick', dumr, dumd)
      end if
      go to 210

c
c **********************************************************************
c *                                                                    *
c *                     linear_elastic option                          *
c *                                                                    *
c **********************************************************************
c
c
 290  continue
      lmtprp(8,matnum) = .true.
      matprp(9,matnum) = 1
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                     vectorized linear and von mises                *
c *                     plasticity model. also used for bar2 and       *
c *                     link2 elements.
c *                                                                    *
c **********************************************************************
c
c
 300  continue
      matprp(9,matnum)= 1
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                 gurson material.                                   *
c *                                                                    *
c **********************************************************************
c
c
 310  continue
      matprp(9,matnum)= 3
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                 mises material.                                    *
c *                                                                    *
c **********************************************************************
c
c
 312  continue
      matprp(9,matnum)  = 3
      matprp(15,matnum) = 0.0
      matprp(16,matnum) = 0.0
      matprp(17,matnum) = 0.0
      lmtprp(18,matnum) = .false.
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 nonlinear elastic model                            *
c *                                                                    *
c **********************************************************************
c
c
 315  continue
      matprp(9,matnum)= 2
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 viscoplasticity power                              *
c *                                                                    *
c **********************************************************************
c
c
 320  continue
      if(.not.numr(matprp(10,matnum))) then
         call errmsg(5,dum,'m_power',dumr,dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 exponent for power-law hardening                   *
c *                                                                    *
c **********************************************************************
c
c
 330  continue
      call inmat_fgm( status, 'n_power'  )
      if ( status  .eq. 0 ) then
          matprp(11,matnum) = fgm_mark
          go to 210
      elseif ( status .eq. 1 ) then
          go to 210
      end if

      if(.not.numr(matprp(11,matnum))) then
         call errmsg(5,dum,'n_power',dumr,dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 reference strain rate                              *
c *                                                                    *
c **********************************************************************
c
c
 340  continue
      if(.not.numr(matprp(12,matnum))) then
         call errmsg(5,dum,'eps_reference',dumr,dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 debug the material model                           *
c *                                                                    *
c **********************************************************************
c
c
 350  continue
      lmtprp(13,matnum) = .true.
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 initial porosity for gurson model                  *
c *                                                                    *
c **********************************************************************
c
c
 360  continue
      if(.not.numr(matprp(14,matnum))) then
         call errmsg(5,dum,'f_o',dumr, dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 q1 for gurson model                                *
c *                                                                    *
c **********************************************************************
c
c
 370  continue
      if(.not.numr(matprp(15,matnum))) then
         call errmsg(5,dum,'q_1',dumr,dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 q2 for gurson model                                *
c *                                                                    *
c **********************************************************************
c
c
 380  continue
      if(.not.numr(matprp(16,matnum))) then
         call errmsg(5,dum,'q_2',dumr,dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 q3 for gurson model                                *
c *                                                                    *
c **********************************************************************
c
c
 390  continue
      if(.not.numr(matprp(17,matnum))) then
         call errmsg(5,dum,'q_3',dumr,dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 nucleation flag for gurson model                   *
c *                                                                    *
c **********************************************************************
c
c
 400  continue
      lmtprp(18,matnum) = .true.
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                 s_n for nucleation in gurson model                 *
c *                                                                    *
c **********************************************************************
c
c
 410  continue
      if(.not.numr(matprp(19,matnum))) then
         call errmsg(5,dum,'s_n',dumr,dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 eps_n for nucleation in gurson model               *
c *                                                                    *
c **********************************************************************
c
c
 420  continue
      if(.not.numr(matprp(20,matnum))) then
         call errmsg(5,dum,'eps_n',dumr,dumd)
      end if
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                 f_n for nucleation in gurson model                 *
c *                                                                    *
c **********************************************************************
c
c
 430  continue
      if(.not.numr(matprp(21,matnum))) then
         call errmsg(5,dum,'f_n',dumr,dumd)
      end if
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                 whether to allow material model to                 *
c *                 cutback the adaptive step size solution            *
c *                 due to excessive reversed plasticity               *
c *                                                                    *
c **********************************************************************
c
c
 440  continue
      lmtprp(22,matnum) = .false.
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                 whether to allow the element to be killed          *
c *                 in a crack growth problem                          *
c *                                                                    *
c **********************************************************************
c
c
 450  continue
      lmtprp(23,matnum) = .true.
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 segmental stress strain curve                      *
c *                    (1)  get the list of curves to use. expand list *
c *                         and save                                   *
c *                    (2)  set segmental flag true                    *
c *                                                                    *
c *                                                                    *
c **********************************************************************
c
c
 460  continue
      call inseg_curve_list( curve_set_no, lmtprp(24,matnum) )
      if ( lmtprp(24,matnum) ) matprp(45,matnum) = curve_set_no
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                 cohesive material for use with interface elements  *
c *                                                                    *
c **********************************************************************
c
c
 465  continue
      matprp(9,matnum) = 4
      call inmat_cohesive( matprp(1,matnum), nonlocal_analysis )
      go to 9998
c
c
c **********************************************************************
c *                                                                    *
c *                 global nonlocal flag - cohesive material           *
c *                                                                    *
c **********************************************************************
c
c
 470  continue
      nonlocal_analysis = .true.
      go to 210
c
c **********************************************************************
c *                                                                    *
c *         properties for creep model                                 *
c *                                                                    *
c **********************************************************************
c
c
 600  continue
      matprp(9,matnum)= 6
      creep_model_used = .true.
      go to 210
 605  continue
      if( numr(matprp(80,matnum)) ) then
      else
         call errmsg(5,dum,'B',dumr,dumd)
         matprp(80,matnum) = 0.0
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *         properties for cyclic plasticity model                     *
c *         model                                                      *
c *                                                                    *
c **********************************************************************
c
c
 800  continue
      matprp(9,matnum)= 5
      call inmat_cyclic( matprp(1,matnum), lmtprp(1,matnum)  )
      go to 9998
c
c **********************************************************************
c *                                                                    *
c *         properties for yueming liang's mises model with            *
c *         hydrogen effects                                           *
c *                                                                    *
c **********************************************************************
c
 1000 continue
      matprp(9,matnum)= 7
      go to 210
c
 1010 continue
      call backsp( 1 )
      yl_prop  = 0
      yl_start = 70
      if( matchs_exact( 'trap_number') ) then
        yl_prop = 1
        go to 1020
      end if
      if( matchs_exact( 'init_disloc_density') ) then
        yl_prop = 2
        go to 1020
      end if
      if( matchs_exact( 'disloc_density_parameter') ) then
        yl_prop = 3
        go to 1020
      end if
      if( matchs_exact( 'lattice_parameter') ) then
        yl_prop = 4
        go to 1020
      end if
      if( matchs_exact( 'binding_energy') ) then
        yl_prop = 5
        go to 1020
      end if
      if( matchs_exact( 'nils_number') ) then
        yl_prop = 6
        go to 1020
      end if
      if( matchs_exact( 'molar_volume') ) then
        yl_prop = 7
        go to 1020
      end if
      if( matchs_exact( 'h_partial_volume') ) then
        yl_prop = 8
        go to 1020
      end if
      if( matchs_exact( 'initial_hydrogen_conc') ) then
        yl_prop = 9
        go to 1020
      end if
      if( matchs_exact( 'softening_parameter') ) then
        yl_prop = 10
        go to 1020
      end if
      write(*,9800) 3  ! really bad if get here
      call die_abort
c
 1020 continue
      if( numr(matprp(yl_start-1+yl_prop,matnum)) ) then
      else
         call errmsg(5,dum,'yl_##',dumr,dumd)
         matprp(yl_start-1+yl_prop,matnum) = 0.0
      end if
      if( local_debug ) then
        write(*,*) '.. yl_start, yl_prop, value: ',
     &     yl_start, yl_prop, matprp(yl_start-1+yl_prop,matnum)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *         properties for general UMAT                                *
c *                                                                    *
c **********************************************************************
c
 1100 continue
      matprp(9,matnum) = 8
      umat_used = .true.
      go to 210
c
 1110 continue
      call backsp( 1 )
      um_prop  = 0
      um_start = 151
      if( matchs_exact( 'um_1') ) then
        um_prop = 1
        go to 1120
      end if
      if( matchs_exact( 'um_2') ) then
        um_prop = 2
        go to 1120
      end if
      if( matchs_exact( 'um_3') ) then
        um_prop = 3
        go to 1120
      end if
      if( matchs_exact( 'um_4') ) then
        um_prop = 4
        go to 1120
      end if
      if( matchs_exact( 'um_5') ) then
        um_prop = 5
        go to 1120
      end if
      if( matchs_exact( 'um_6') ) then
        um_prop = 6
        go to 1120
      end if
      if( matchs_exact( 'um_7') ) then
        um_prop = 7
        go to 1120
      end if
      if( matchs_exact( 'um_8') ) then
        um_prop = 8
        go to 1120
      end if
      if( matchs_exact( 'um_9') ) then
        um_prop = 9
        go to 1120
      end if
      if( matchs_exact( 'um_10') ) then
        um_prop = 10
        go to 1120
      end if
      if( matchs_exact( 'um_11') ) then
        um_prop = 11
        go to 1120
      end if
      if( matchs_exact( 'um_12') ) then
        um_prop = 12
        go to 1120
      end if
      if( matchs_exact( 'um_13') ) then
        um_prop = 13
        go to 1120
      end if
      if( matchs_exact( 'um_14') ) then
        um_prop = 14
        go to 1120
      end if
      if( matchs_exact( 'um_15') ) then
        um_prop = 15
        go to 1120
      end if
      if( matchs_exact( 'um_16') ) then
        um_prop = 16
        go to 1120
      end if
      if( matchs_exact( 'um_17') ) then
        um_prop = 17
        go to 1120
      end if
      if( matchs_exact( 'um_18') ) then
        um_prop = 18
        go to 1120
      end if
      if( matchs_exact( 'um_19') ) then
        um_prop = 19
        go to 1120
      end if
      if( matchs_exact( 'um_20') ) then
        um_prop = 20
        go to 1120
      end if
      if( matchs_exact( 'um_21') ) then
        um_prop = 21
        go to 1120
      end if
      if( matchs_exact( 'um_22') ) then
        um_prop = 22
        go to 1120
      end if
      if( matchs_exact( 'um_23') ) then
        um_prop = 23
        go to 1120
      end if
      if( matchs_exact( 'um_24') ) then
        um_prop = 24
        go to 1120
      end if
      if( matchs_exact( 'um_25') ) then
        um_prop = 25
        go to 1120
      end if
      if( matchs_exact( 'um_26') ) then
        um_prop = 26
        go to 1120
      end if
      if( matchs_exact( 'um_27') ) then
        um_prop = 27
        go to 1120
      end if
      if( matchs_exact( 'um_28') ) then
        um_prop = 28
        go to 1120
      end if
      if( matchs_exact( 'um_29') ) then
        um_prop = 29
        go to 1120
      end if
      if( matchs_exact( 'um_30') ) then
        um_prop = 30
        go to 1120
      end if
      if( matchs_exact( 'um_31') ) then
        um_prop = 31
        go to 1120
      end if
      if( matchs_exact( 'um_32') ) then
        um_prop = 32
        go to 1120
      end if
      if( matchs_exact( 'um_33') ) then
        um_prop = 33
        go to 1120
      end if
      if( matchs_exact( 'um_34') ) then
        um_prop = 34
        go to 1120
      end if
      if( matchs_exact( 'um_35') ) then
        um_prop = 35
        go to 1120
      end if
      if( matchs_exact( 'um_36') ) then
        um_prop = 36
        go to 1120
      end if
      if( matchs_exact( 'um_37') ) then
        um_prop = 37
        go to 1120
      end if
      if( matchs_exact( 'um_38') ) then
        um_prop = 38
        go to 1120
      end if
      if( matchs_exact( 'um_39') ) then
        um_prop = 39
        go to 1120
      end if
      if( matchs_exact( 'um_40') ) then
        um_prop = 40
        go to 1120
      end if
      if( matchs_exact( 'um_41') ) then
        um_prop = 41
        go to 1120
      end if
      if( matchs_exact( 'um_42') ) then
        um_prop = 42
        go to 1120
      end if
      if( matchs_exact( 'um_43') ) then
        um_prop = 43
        go to 1120
      end if
      if( matchs_exact( 'um_44') ) then
        um_prop = 44
        go to 1120
      end if
      if( matchs_exact( 'um_45') ) then
        um_prop = 45
        go to 1120
      end if
      if( matchs_exact( 'um_46') ) then
        um_prop = 46
        go to 1120
      end if
      if( matchs_exact( 'um_47') ) then
        um_prop = 47
        go to 1120
      end if
      if( matchs_exact( 'um_48') ) then
        um_prop = 48
        go to 1120
      end if
      if( matchs_exact( 'um_49') ) then
        um_prop = 49
        go to 1120
      end if
      if( matchs_exact( 'um_50') ) then
        um_prop = 50
        go to 1120
      end if
      write(*,9800) 4 ! really bad if get here
      call die_abort
c
 1120 continue
      if( numd(dmatprp(um_start-1+um_prop,matnum)) ) then
      else
         call errmsg(5,dum,'um_##',dumr,dumd)
         dmatprp(um_start-1+um_prop,matnum) = 0.0d00
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *         properties for CP model                                    *
c *                                                                    *
c **********************************************************************
c
 2000 continue
      matprp(9,matnum)= 10
      call inmat_cp( matnum  )
      go to 9998
c
c ***********************************************************************
c *                                                                     *
c *         properties for the interface damage model                   *
c *                                                                     *
c ***********************************************************************
c
 2100 continue
      write(*,*) '... model not yet included in newest WARP3D'
      call die_abort
      matprp(9, matnum) = 11
      call inmat_inter( matnum)
      go to 9998

c
c **********************************************************************
c *                                                                    *
c *         link material stiffness, mass values                       *
c *                                                                    *
c **********************************************************************
c
 2200 continue
      if( numr(kx_stiff) ) then
         if( kx_stiff < 0.0e0 ) then
           call errmsg( 87, dum, dums, dumr, dumd )
           go to 210
         end if
         matprp(148,matnum) = kx_stiff
         go to 210
      end if
      call errmsg( 124, dum, 'kx_stiff', dumr, dumd )
      go to 210
 2205 continue
      if( numr(ky_stiff) ) then
         if( ky_stiff < 0.0e0 ) then
           call errmsg( 87, dum, dums, dumr, dumd )
           go to 210
         matprp(149,matnum) = ky_stiff
         end if
         go to 210
      end if
      call errmsg( 124, dum, 'ky_stiff', dumr, dumd )
      go to 210
 2210 continue
      if( numr(kz_stiff) ) then
         if( kz_stiff < 0.0e0 ) then
           call errmsg( 87, dum, dums, dumr, dumd )
           go to 210
         end if
         matprp(150,matnum) = kz_stiff
         go to 210
      end if
      call errmsg( 124, dum, 'kz_stiff', dumr, dumd )
      go to 210
 2215 continue
      if( numr(link_mass) ) then
         if( link_mass < 0.0e0 ) then
           call errmsg( 140, dum, dums, dumr, dumd )
           go to 210
         end if
         matprp(7,matnum) = link_mass
         go to 210
      end if
      call errmsg( 352, dum, dums, dumr, dumd )
      go to 210
 2220 continue
      if( numr(link_stiff) ) then
         if( link_stiff < 0.0e0 ) then
           call errmsg( 87, dum, dums, dumr, dumd )
           matprp(148,matnum) = 0.0e0
           matprp(149,matnum) = 0.0e0
           matprp(150,matnum) = 0.0e0
           go to 210
         end if
         matprp(148,matnum) = link_stiff
         matprp(149,matnum) = link_stiff
         matprp(150,matnum) = link_stiff
         go to 210
      end if
      call errmsg( 124, dum, 'link_stiff', dumr, dumd )
      go to 210
c
 9997 sbflg1= .true.
      go to 9999
c
 9998 sbflg1= .false.
c
 9999 sbflg2= .true.
c
      return
 9800 format('>>>> FATAL ERROR: in inmat @',i3,
     &    /, '                  job aborted.'//)
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine inseg_curve_list             *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 2/28/00                    *
c     *                                                              *
c     *     this subroutine supervises and conducts the input of the *
c     *     list of segmental stress-strain curves for use with      *
c     *     this material                                            *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine inseg_curve_list( cur_set_no, ok  )
      use global_data ! old common.main
c
      use segmental_curves
c
      implicit integer (a-z)
c
c                       parameters
c
      logical ok
c                       local declarations
c
      dimension intlst(100), tcurve_list(max_seg_curves)
      character dums*10
      logical local_debug, true
      real dumr
      double precision
     &   dumd, tcurve_value(max_seg_curves), strain_1, strain_2,
     &   tol_strain
      data local_debug, tol_strain / .false., 1.0e-05 /
c
c                 get the integerlist of segmental stress-strain curves
c                 to be associated with this material. we enter with scan
c                 test routines having a last test true. the list
c                 routine does not look at the test flag so we must force
c                 scanner to first look at the list.
c
      ok         = .false.
      cur_set_no = 0
      if ( local_debug ) write(*,*) '>> inside inseg_curve_list'
c
      call scan
      call trlist( intlst, 100, 0, lenlst, errnum )
c
c                       branch on the return code from trlist. a
c                       value of 1 indicates no error. a value of
c                       2 indicates that the parse rules failed in
c                       the list. a value of 3 indicates that the
c                       list overflowed its maximum length of mxlsz.
c                       a value of 4 indicates that no list was found.
c                       in these last 3 cases, the rest of the card
c                       will be ignored and a new card will be sought.
c
      if( errnum .eq. 2 ) then
         param = 1
         call errmsg( 24, param, dums, dumr, dumd )
         call errmsg2( 8, param, dums, dumr, dumd )
         return
      else if( errnum .eq.3 ) then
         param = 2
         call errmsg(24,param,dums,dumr,dumd)
         call errmsg2( 8, param, dums, dumr, dumd )
         return
      else if( errnum .eq. 4 ) then
         param =4
         call errmsg(24,param,dums,dumr,dumd)
         call errmsg2( 8, param, dums, dumr, dumd )
         return
      else
         if( errnum .eq. 1 ) then
            call backsp(1)
            go to 200
         end if
         param= 3
         call errmsg(24,param,dums,dumr,dumd)
         call errmsg2( 8, param, dums, dumr, dumd )
         return
      end if
c
c                 an integer list of curve numbers was found. create the
c                 next "set" of curves in the curve set table. expand the
c                 list of curves and store as the set list
c
 200  continue
      if ( local_debug ) write(*,*) '>> got integer list..'
      num_seg_curve_sets = num_seg_curve_sets + 1
      if ( num_seg_curve_sets .gt. max_seg_curve_sets ) then
        call errmsg2( 9, max_seg_curve_sets, dums, dumr, dumd )
        num_seg_curve_sets = num_seg_curve_sets -1
        return
      end if
c
      count  = 0
      icn    = 0
      iplist = 1
c
 300  continue
      call trxlst( intlst, lenlst, iplist, icn, curve_no )
      if ( local_debug ) write(*,*)
     &           '>> next curve no. in list: ',curve_no
      count = count + 1
c
      if ( curve_no .le. 0 .or. curve_no .gt. max_seg_curves ) then
        call errmsg2( 10, curve_no, dums, dumr, dumd )
        seg_curve_table(1:max_seg_curves,num_seg_curve_sets) = 0
        num_seg_curve_sets = num_seg_curve_sets -1
        return
      end if
c
      if( .not. seg_curve_def(curve_no) ) then
        call errmsg2( 10, curve_no, dums, dumr, dumd )
        seg_curve_table(1:max_seg_curves,num_seg_curve_sets) = 0
        num_seg_curve_sets = num_seg_curve_sets -1
        return
      end if
c
      seg_curve_table(count+1,num_seg_curve_sets) = curve_no
      seg_curve_table(1,num_seg_curve_sets)       = count
      if( iplist .ne. 0 ) go to 300
c
c                 list of curves expanded and checked. define curve
c                 set number and ok flag for calling routine. set
c                 scanner test to true so we will move on to next
c                 data item for material
c
      ok         = .true.
      cur_set_no = num_seg_curve_sets
      if ( true(dumr) )call splunj
      if ( local_debug ) then
          write(*,*) '>> curve set table....'
          do ii = 1, num_seg_curve_sets
            write(*,*) '      > curve set: ', ii
            count = seg_curve_table(1,ii)
            write(*,9100)  seg_curve_table(1:count+1,ii)
          end do
      end if
c
c                 examine the just defined set of curves. if the set
c                 has more than 1 curve, the curves must be either
c                 all temperature dependent or rate dependent. In these
c                 cases, the listed curves should be sorted in order
c                 of increasing value of the dependency variable.
c
      if ( seg_curve_table(1,cur_set_no) .eq. 1 ) return
c
c                 all curves in the set must have the same dependency
c                 type.
c
      curve_no         =  seg_curve_table(2,cur_set_no)
      curve_type       =  seg_curves_type(curve_no)
      number_of_curves =  seg_curve_table(1,cur_set_no)
      do i = 1, number_of_curves
         now_curve = seg_curve_table(i+1,cur_set_no)
         if ( seg_curves_type(now_curve) .eq. curve_type ) cycle
         call errmsg2( 11, curve_no, dums, dumr, dumd )
         ok = .false.
         seg_curve_table(1:max_seg_curves,cur_set_no) = 0
         num_seg_curve_sets = num_seg_curve_sets -1
         return
      end do
c
c                 the curves in the set must all have the same number
c                 of points and the strain values must be the same.
c
c
      curve_no         =  seg_curve_table(2,cur_set_no)
      curve_type       =  seg_curves_type(curve_no)
      number_of_curves =  seg_curve_table(1,cur_set_no)
      do j = 1, number_of_curves
         now_curve = seg_curve_table(j+1,cur_set_no)
         if ( num_seg_points(curve_no) .ne.
     &        num_seg_points(now_curve) ) then
           call errmsg2( 13, now_curve, dums, dumr, dumd )
           ok = .false.
           seg_curve_table(1:max_seg_curves,cur_set_no) = 0
           num_seg_curve_sets = num_seg_curve_sets -1
           return
         end if
         do i = 1, num_seg_points(curve_no)
            strain_1 = seg_curves(i,1,curve_no)
            strain_2 = seg_curves(i,1,now_curve)
            if ( abs(strain_1-strain_2) .gt. tol_strain*strain_1 ) then
              call errmsg2( 14, now_curve, dums, dumr, dumd )
              ok = .false.
              seg_curve_table(1:max_seg_curves,cur_set_no) = 0
              num_seg_curve_sets = num_seg_curve_sets -1
              return
            end if
         end do
      end do
c
c                 all curves of the set are the same dependency type,
c                 they have the same number of points and the
c                 strain values are all the same.
c
c                 sort the list of curves by ascending value of
c                 the dependency variable. if they are curve type 0,
c                 there is no dependency value. just return.
c
      if ( curve_type .eq. 0 ) then
         call errmsg2( 12, curve_no, dums, dumr, dumd )
         return
      end if
      do i = 2, number_of_curves + 1
         curve_no          = seg_curve_table(i,cur_set_no)
         tcurve_list(i-1)  = curve_no
         tcurve_value(i-1) = seg_curves_value(curve_no)
      end do
c
      if ( local_debug ) then
        write(*,*) '>> list of curves and values prior to sort:'
        do i = 1, number_of_curves
           write(*,*) '> curve no, depend. value: ',tcurve_list(i),
     &                tcurve_value(i)
        end do
      end if
c
      call warp3d_sort_float( number_of_curves, tcurve_value,
     &                        tcurve_list )
c
c                 put the sorted list of curve numbers back into
c                 the seg_curve_table for the set.
c
      seg_curve_table(2:number_of_curves+1,cur_set_no) =
     &             tcurve_list(1:number_of_curves)
c
c
      if ( local_debug ) then
        write(*,*) '>> sorted list of curves for the set:'
        do i = 1, number_of_curves
           write(*,*) '> curve no: ',seg_curve_table(i+1,cur_set_no)
        end do
      end if
c
 9100 format(10x,10i3)
      end
c     ****************************************************************
c     *                                                              *
c     *                    subroutine inmat_fgm                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 4/9/2018 rhd               *
c     *                                                              *
c     *     check input stream for 'fgm' and return status           *
c     *      = -1 (no string), = 1 (string but no match), = 0 ok     *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine inmat_fgm( status, prop_id )
      use main_data, only : fgm_node_values_used
c
      implicit none
c
      integer :: status
      character(len=*) :: prop_id
c
      integer :: idummy, ncstring
      logical, external :: string
      logical :: ok
      character(len=80) :: text
      real :: dumr
      double precision :: dumd
c
      status = -1
      if( .not. string(idummy) ) return
      call entits( text, ncstring )
c
      if( ncstring .lt. 3 ) then
         call errmsg2(30,idummy,prop_id,dumr,dumd)
         status = 1
         return
      end if
      ok = text(1:3) .eq. 'fgm' .or. text(1:3) .eq. 'FGM'
      if( ok ) then
         status = 0
         fgm_node_values_used = .true.
      else
         status = 1
         call errmsg2(30,idummy,prop_id,dumr,dumd)
      end if
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine inmat_cp                     *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 3/21/12                    *
c     *                                                              *
c     *     input properties for the crystal plasticity model (#10)  *
c     *                                                              *
c     ****************************************************************
c
      subroutine inmat_cp(matnum)
      use main_data, only : matprp, lmtprp, imatprp, dmatprp, smatprp
      implicit integer (a-z)
      integer, intent(in) :: matnum
      integer :: dumi, nc
      real :: dumr
      double precision :: dumd
      character :: dums, lab*24, filen*24
      logical :: reading
      logical, external :: matchs_exact, isstring, label, numr, numd,
     &                     matchs, numi, endcrd

c           See above for a detailed summary of each material option, in
c           general, handle input for 6,7,9,13,22,25,100-113
c
c           Set defaults
      lmtprp(13,matnum)=.false.
      lmtprp(22,matnum)=.true.
      lmtprp(25,matnum)=.false.
c           Read in properties
      reading = .true.
      do while (reading)
            if ( matchs_exact('angle_convention')) then
                  if (.not. label(dumi)) then
                        call errmsg(356,dumi,dums,dumr,dumd)
                  else
                        lab = ' '
                        call entits(lab,nc)
                  end if
c
                  if (lab(1:nc) .eq. 'kocks') then
                        imatprp(102,matnum) = 1
                  else
                        call errmsg(357,dumi,lab(1:nc),dumr,dumd)
                  end if
            elseif ( matchs_exact('alpha')) then
                  if (.not. numr(matprp(6,matnum))) then
                        call errmsg(5,dumi,'alpha',dumr,dumd)
                  end if
            elseif ( matchs_exact('rho')) then
                  if (.not. numr(matprp(7,matnum))) then
                        call errmsg(5,dumi,'rho',dumr,dumd)
                  end if
            elseif ( matchs_exact('tolerance')) then
                  if (.not. numd(dmatprp(100,matnum))) then
                        call errmsg(5,dumi,'tolerance',dumr,dumd)
                  end if
            elseif ( matchs_exact('n_crystals')) then
                  if (.not. numi(imatprp(101,matnum))) then
                        call errmsg(5,dumi,'n_crystals',dumr,dumd)
                  end if
            elseif ( matchs_exact('angle_type')) then
                  if (.not. label(dumi)) then
                        call errmsg(358,dumi,dums,dumr,dumd)
                  else
                        lab = ' '
                        call entits(lab,nc)
                  end if
c
                  if (lab(1:nc) .eq. 'degrees') then
                        imatprp(103,matnum) = 1
                  elseif (lab(1:nc) .eq. 'radians') then
                        imatprp(103,matnum) = 2
                  else
                        call errmsg(359,dumi,lab(1:nc),dumr,dumd)
                  end if
            elseif ( matchs_exact('crystal_input')) then
                  if (.not. label(dumi)) then
                        call errmsg(360,dumi,dums,dumr,dumd)
                  else
                        lab = ' '
                        call entits(lab,nc)
                  end if
c
                  if (lab(1:nc) .eq. 'single') then
                        imatprp(104,matnum) = 1
                  elseif (lab(1:nc) .eq. 'file') then
                        imatprp(104,matnum) = 2
                  else
                        call errmsg(361,dumi,lab(1:nc),dumr,dumd)
                  end if
            elseif ( matchs_exact('crystal_type')) then
                  if (.not. numi(imatprp(105,matnum))) then
                        call errmsg(5,dumi,'crystal_type',dumr,dumd)
                  end if
            elseif ( matchs_exact('orientation_input')) then
                  if (.not. label(dumi)) then
                        call errmsg(360,dumi,dums,dumr,dumd)
                  else
                        lab = ' '
                        call entits(lab,nc)
                  end if
c
                  if (lab(1:nc) .eq. 'single') then
                        imatprp(107,matnum) = 1
                  elseif (lab(1:nc) .eq. 'file') then
                        imatprp(107,matnum) = 2
                  else
                        call errmsg(361,dumi,lab(1:nc),dumr,dumd)
                  end if
            elseif ( matchs_exact('angles') ) then
                  if (.not. (numd(dmatprp(108,matnum)) .and.
     &                       numd(dmatprp(109,matnum)) .and.
     &                       numd(dmatprp(110,matnum)))) then
                        call errmsg(362,dumi,dums,dumr,dumd)
                  end if
            elseif ( matchs_exact('filename')) then
                  call doscan()
                  if (.not. isstring()) then
                        call errmsg(363,dumi,dums,dumr,dumd)
                  else
                        filen = ' '
                        call entits(filen,nc)
                        if (nc .gt. 24) then
                              call errmsg(365,dumi,dums,dumr,dumd)
                        end if
                        call scan()
                        smatprp(112,matnum) = filen
                  end if
            elseif ( matchs_exact('debug')) then
                  if (.not. label(dumi)) then
                        call errmsg(364,dumi,dums,dumr,dumd)
                  else
                        lab = ' '
                        call entits(lab,nc)
                        if (lab(1:nc) .eq. 'on') then
                              lmtprp(13,matnum) = .true.
                        elseif (lab(1:nc) .eq. 'off') then
                              lmtprp(13,matnum) = .false.
                        else
                              call errmsg(364,dumi,dums,dumr,dumd)
                        end if
                  end if
            elseif ( endcrd(dum) ) then
                  reading = .false.
                  cycle
            elseif ( matchs(',',1) ) then
                  call readsc()
            else
                  call entits(lab,nc)
                  call errmsg(355,dumi,lab(1:nc),dumr,dumd)
                  call scan()
                  cycle
                  end if
            end do

      return
      end subroutine

c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine inmat_inter                  *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 3/10/14                    *
c     *                                                              *
c     *     input properties for the crystal interface damage model  *
c     *                                                              *
c     ****************************************************************
c
      subroutine inmat_inter(matnum)
      use global_data ! old common.main
      use main_data, only : matprp, lmtprp, imatprp, dmatprp, smatprp
      implicit integer (a-z)
      integer, intent(in) :: matnum
      integer :: dumi, nc
      real :: dumr
      double precision :: dumd
      character :: dums, lab*24, filen*24, name*80, mname*24
      logical :: reading, found
      logical, external :: matchs_exact, isstring, label, numr, numd,
     &                     matchs, numi, endcrd, scanms

c           See above for a detailed summary of each material option, in
c           general, handle input for 6,7,13,22,25,100,102,112,115-129
c
c           Set defaults
      lmtprp(13,matnum)=.false.
      lmtprp(22,matnum)=.true.
      lmtprp(25,matnum)=.false.
      imatprp(107,matnum) = 2
      imatprp(101,matnum) = 50
      mname(1:8) = 'nodefalt'
c           Read in properties
      reading = .true.
      do while (reading)
            if ( matchs_exact('angle_convention')) then
                  if (.not. label(dumi)) then
                        call errmsg(356,dumi,dums,dumr,dumd)
                  else
                        lab = ' '
                        call entits(lab,nc)
                  end if
c
                  if (lab(1:nc) .eq. 'kocks') then
                        imatprp(102,matnum) = 1
                  else
                        call errmsg(357,dumi,lab(1:nc),dumr,dumd)
                  end if
            elseif ( matchs_exact('max_crystals')) then
                  if (.not. numi(imatprp(101,matnum))) then
                        call errmsg(5,dumi,'max_crystals',dumr,dumd)
                  end if
            elseif ( matchs_exact('tolerance')) then
                  if (.not. numd(dmatprp(100,matnum))) then
                        call errmsg(5,dumi,'tolerance',dumr,dumd)
                  end if
            elseif ( matchs_exact('angle_type')) then
                  if (.not. label(dumi)) then
                        call errmsg(358,dumi,dums,dumr,dumd)
                  else
                        lab = ' '
                        call entits(lab,nc)
                  end if
c
                  if (lab(1:nc) .eq. 'degrees') then
                        imatprp(103,matnum) = 1
                  elseif (lab(1:nc) .eq. 'radians') then
                        imatprp(103,matnum) = 2
                  else
                        call errmsg(359,dumi,lab(1:nc),dumr,dumd)
                  end if
            elseif ( matchs_exact('crystal_type')) then
                  if (.not. numi(imatprp(105,matnum))) then
                        call errmsg(5,dumi,'crystal_type',dumr,dumd)
                  end if
            elseif ( matchs_exact('o_file')) then
                  call doscan()
                  if (.not. isstring()) then
                        call errmsg(363,dumi,dums,dumr,dumd)
                  else
                        filen = ' '
                        call entits(filen,nc)
                        if (nc .gt. 24) then
                              call errmsg(365,dumi,dums,dumr,dumd)
                        end if
                        call scan()
                        smatprp(112,matnum) = filen
                  end if
            elseif ( matchs_exact('debug')) then
                  if (.not. label(dumi)) then
                        call errmsg(364,dumi,dums,dumr,dumd)
                  else
                        lab = ' '
                        call entits(lab,nc)
                        if (lab(1:nc) .eq. 'on') then
                              lmtprp(13,matnum) = .true.
                        elseif (lab(1:nc) .eq. 'off') then
                              lmtprp(13,matnum) = .false.
                        else
                              call errmsg(364,dumi,dums,dumr,dumd)
                        end if
                  end if
            elseif ( matchs_exact('macroscale')) then
                  if (.not. label(dumi)) then
                        call errmsg(364,dumi,dums,dumr,dumd)
                  else
                        name  = ' '
                        mname = ' '
                        call entits(name,nc)
                        if( nc .gt. 24 ) nc = 24
                        mname(1:nc) = name(1:nc)
                        found = .false.
                        matn = mathed/two16
 5810                   if( matn .eq. 32460 ) go to 5820
                         if( scanms(matnam(matn),mname,24) ) then
                         imatprp(115, matnum) = matn
                         found  = .true.
                         go to 5820
                        end if
                        matn = matlst(matn)/two16
                        go to 5810
 5820                   if( .not. found ) then
                          call errmsg(351,dum,dums,dumr,dumd)
                          mname(1:8) = 'nodefalt'
                        end if
                  end if
            elseif ( matchs_exact('s_v')) then
                  if (.not. numd(dmatprp(116,matnum))) then
                        call errmsg(5,dumi,'s_v(1)',dumr,dumd)
                  end if
                  if (.not. numd(dmatprp(117,matnum))) then
                        call errmsg(5,dumi,'s_v(2)',dumr,dumd)
                  end if
                  if (.not. numd(dmatprp(118,matnum))) then
                        call errmsg(5,dumi,'s_v(3)',dumr,dumd)
                  end if
            elseif ( matchs_exact('l_v')) then
                  if (.not. numd(dmatprp(119,matnum))) then
                        call errmsg(5,dumi,'l_v(1)',dumr,dumd)
                  end if
                  if (.not. numd(dmatprp(120,matnum))) then
                        call errmsg(5,dumi,'l_v(2)',dumr,dumd)
                  end if
                  if (.not. numd(dmatprp(121,matnum))) then
                        call errmsg(5,dumi,'l_v(3)',dumr,dumd)
                  end if
            elseif ( matchs_exact('t_v')) then
                  if (.not. numd(dmatprp(122,matnum))) then
                        call errmsg(5,dumi,'t_v(1)',dumr,dumd)
                  end if
                  if (.not. numd(dmatprp(123,matnum))) then
                        call errmsg(5,dumi,'t_v(2)',dumr,dumd)
                  end if
                  if (.not. numd(dmatprp(124,matnum))) then
                        call errmsg(5,dumi,'t_v(3)',dumr,dumd)
                  end if
            elseif ( matchs_exact('l_s')) then
                  if (.not. numd(dmatprp(126,matnum))) then
                        call errmsg(5,dumi,'l_s',dumr,dumd)
                  end if
            elseif ( matchs_exact('l_l')) then
                  if (.not. numd(dmatprp(127,matnum))) then
                        call errmsg(5,dumi,'l_l',dumr,dumd)
                  end if
            elseif ( matchs_exact('l_t')) then
                  if (.not. numd(dmatprp(128,matnum))) then
                        call errmsg(5,dumi,'l_t',dumr,dumd)
                  end if
            elseif ( matchs_exact('alpha_dmg')) then
                  if (.not. numd(dmatprp(129,matnum))) then
                        call errmsg(5,dumi,'alpha_dmg',dumr,dumd)
                  end if
            elseif ( endcrd(dum) ) then
                  reading = .false.
                  cycle
            elseif ( matchs(',',1) ) then
                  call readsc()
            else
                  call entits(lab,nc)
                  call errmsg(355,dumi,lab(1:nc),dumr,dumd)
                  call scan()
                  cycle
                  end if
            end do
      return
      end subroutine

c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine inmat_cyclic                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 7/30/2011                  *
c     *                                                              *
c     *     input properties for the cyclic plasticty model (#5)     *
c     *                                                              *
c     ****************************************************************
c
      subroutine inmat_cyclic( matprp, lmtprp )
c
      implicit integer (a-z)
c
c                       parameter declarations (the column of the
c                       table is passed for this material)
c
      real :: matprp(*)
      logical :: lmtprp(*)
c
c                       local declarations
c
      logical :: local_debug, true, matchs, matchs_exact, endcrd,
     &           numr
      real ::  dumr
      double precision :: dumd
      real root23, onept5, value, twothrds
      character(len=1) :: dums
      data root23, onept5, twothrds / 0.81649, 1.5, 0.6666667 /
c
c                       enter here having seen the properties and
c                       cyclic keywords at the beginning of a line.
c                       we leave at the end of a logical line having
c                       processed all property values for the
c                       cyclic plasticty model.
c
c                       the model supports two "options":
c                        nonlinear_hardening (formerly known as the
c                                 Frederick-Armstrong option)
c                        generalized_plasticity
c
c
c                  Changes to support major overhaul the the model
c                  in summer 2009, summer 2011 (for temperature
c                  dependent generalized_plasticity).
c
c                  Properties for the nonlinear_hardening option:
c                  ---------------------------------------------
c
c                  Input names: E, nu, rho, yld_pt, sig_tol,
c                               q_u, b_u, h_u, gamma_u
c
c                  mm05 expects as input: q_u, b(3D), h(3D), gamma(3D).
c                  the user's *_u values are changed here to the
c                  expected values before storing into properties.
c                       b(3D) = sqrt(2/3) * b_u
c                       h(3D) = (2/3) * h_u
c                       gamma(3D) = sqrt(2/3) *  gamma_u
c
c
c                  Properties for the generalized_plasticity option:
c                  -------------------------------------------------
c
c                  Input names: E, nu, rho, yld_pt, sig_tol,
c                               gp_h_u, gp_tau, gp_beta_u, gp_delta_u
c                               (changes in 7/2011 add h_u, tau and
c                                delete h_iu, h_ku)
c
c                  mm05 code expects the uniaxial (*_u) values to be
c                  provided from the user input. no conversions are
c                  needed here.
c
c                  storage layout:
c
c                              FA option      GP option
c                     55        q_u             gp_h_u
c                     56        b_u             gp_tau
c                     57        h_u             gp_beta_u
c                     58        1.0              -1.0
c                     59        gamma_u         gp_delta_u
c                     60 sig_tol
c                     61-64 <available>
c
      local_debug = .false.
      if( local_debug ) write(*,*) '.... entered inmat_cyclic ....'
      cyclic_start = 55
c
      matprp(cyclic_start:cyclic_start-1+6) = 0.0
      matprp(cyclic_start-1+4) = 0.0
c
      call reset
      if( true( dummy ) ) call splunj
      if( matchs('properties',10) ) call splunj
      if( matchs('cyclic',6) ) call splunj
      if( local_debug ) write(*,*) '  @ 1'
      option_set = 0
      go to 210
c
 206  continue
      call readsc

 210  continue
      if( local_debug ) write(*,*) '  @ 2'
      if( matchs_exact( 'type' ) ) go to 210
      if( matchs_exact( 'e'  ) ) go to 300
      if( matchs_exact( 'nu' ) ) go to 310
      if( matchs_exact( 'rho' ) ) go to 320
      if( matchs_exact( 'alpha' ) ) go to 325
      if( matchs_exact( 'yld_pt' ) ) go to 330
      if( matchs_exact( 'q_u' ) ) go to 340
      if( matchs_exact( 'b_u' ) ) go to 350
      if( matchs_exact( 'h_u' ) ) go to 360
      if( matchs_exact( 'nonlinear_hardening' ) ) go to 370
      if( matchs_exact( 'generalized_plasticity' ) ) go to 380
      if( matchs_exact( 'gamma_u' ) ) go to 390
      if( matchs_exact( 'sig_tol' ) ) go to 400
      if( matchs_exact( 'gp_h_u' ) ) go to 410
      if( matchs_exact( 'gp_tau' ) ) go to 420
      if( matchs_exact( 'gp_delta_u' ) ) go to 430
      if( matchs_exact( 'gp_beta_u' ) ) go to 440
      if( matchs( 'curves', 5) ) go to 460
c
      if( matchs(',',1)                ) go to 215
c
c                       there is no match with existing properties.
c                       check for end of card. if not, print error
c                       message, skip over bad item and keep going.
c
      if( endcrd(dum) ) then
         if( local_debug ) then
           write(*,9000) matprp(1), matprp(2), matprp(7), matprp(5)
           do i = 1, 6
             write(*,9100) i, matprp(cyclic_start-1+i)
           end do
         end if
         return
      else
         call errmsg(4,dum,dums,dumr,dumd)
         if( true(dum) ) go to 210
      end if
c
c                       comma separator. if at the end of a line,
c                       denotes continuation. otherwise, ignore.
c
 215  continue
      if( endcrd(dum) ) then
         go to 206
      else
         go to 210
      end if
c
c
c                      we have a valid id for a cyclic plasticty
c                      property. get the value, check validity. store.
c
c
c * ********************************************************************
c *                                                                    *
c *                     input modulus of elasticity                    *
c *                                                                    *
c **********************************************************************
c
 300  continue
      if( .not. numr(matprp(1)) ) then
         call errmsg(5,dum,'e',dumr,dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                     input poisson's ratio                          *
c *                                                                    *
c **********************************************************************
c
c
 310  continue
      if( .not. numr(matprp(2)) ) then
         call errmsg(5,dum,'nu',dumr,dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                     input the mass density                         *
c *                                                                    *
c **********************************************************************
c
 320  continue
      if( .not. numr(matprp(7)) ) then
         call errmsg(5,dum,'rho',dumr,dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                     input the thermal expansion                    *
c *                                                                    *
c **********************************************************************
c
 325  continue
      if( .not. numr(matprp(6)) ) then
         call errmsg(5,dum,'alpha',dumr,dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                     input initial yield point                      *
c *                                                                    *
c **********************************************************************
c
 330  continue
      if( .not. numr(matprp(5)) ) then
         call errmsg(5,dum,'yld_pt',dumr,dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                     input q_u for nonlinear hardening option       *
c *                                                                    *
c **********************************************************************
c
 340  continue
      if( .not. numr(value) ) then
         call errmsg(5,dum,'q_u',dumr,dumd)
         go to 210
      end if
      if( option_set .eq. 1 ) then
        matprp(cyclic_start-1+1) = value
      else
        if( option_set .eq. 0 ) then
           call inmat_err_message( 2, 'q_u', 3 )
        else if( option_set .eq. 2 ) then
           call inmat_err_message( 1, 'q_u', 3 )
        end if
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                     input b_u for nonlinear hardening option       *
c *                                                                    *
c **********************************************************************
c
 350  continue
      if( .not. numr(value) ) then
         call errmsg(5,dum,'b_u',dumr,dumd)
         go to 210
      end if
      if( option_set .eq. 1 ) then
          matprp(cyclic_start-1+2) = value * root23
      else
       if( option_set .eq. 0 ) then
           call inmat_err_message( 2, 'b_u', 3 )
        else if( option_set .eq. 2 ) then
           call inmat_err_message( 1, 'b_u', 3 )
        end if
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                     input h_u for nonlinear hardening option       *
c *                                                                    *
c **********************************************************************
c
 360  continue
      if( .not. numr(value) ) then
         call errmsg(5,dum,'h_u',dumr,dumd)
         go to 210
      end if
      if( option_set .eq. 1 ) then
         matprp(cyclic_start-1+3) = value * twothrds
      else
       if( option_set .eq. 0 ) then
           call inmat_err_message( 2, 'h_u', 3 )
        else if( option_set .eq. 2 ) then
           call inmat_err_message( 1, 'h_u', 3 )
        end if
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                     set nonlinear_hardening option flag            *
c *                     (formerly called Frederick-Armstrong)          *
c *                                                                    *
c **********************************************************************
c
 370  continue
      matprp(cyclic_start-1+4) = 1.0
      option_set = 1
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                     set generalized_plasticity option flag         *
c *                                                                    *
c **********************************************************************
c
 380  continue
      matprp(cyclic_start-1+4) = -1.0
      option_set = 2
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                  input gamma_u for nonlininear hardening option    *
c *                                                                    *
c **********************************************************************
c
 390  continue
      if( .not. numr(value) ) then
         call errmsg(5,dum,'gamma_u',dumr,dumd)
         go to 210
      end if
      if( option_set .eq. 1 ) then
         matprp(cyclic_start-1+5) = value * root23
      else
       if( option_set .eq. 0 ) then
           call inmat_err_message( 2, 'gamma_u', 7 )
        else if( option_set .eq. 2 ) then
           call inmat_err_message( 1, 'gamma_u', 7 )
        end if
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                     input sig_tol for both options                 *
c *                                                                    *
c **********************************************************************
c
 400  continue
      if( .not. numr(matprp(cyclic_start-1+6)) ) then
         call errmsg(5,dum,'sig_tol',dumr,dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                     input gp_h_u generalized hardening option      *
c *                                                                    *
c **********************************************************************
c
 410  continue
      if( .not. numr(value) ) then
         call errmsg(5,dum,'gp_h_u',dumr,dumd)
         go to 210
      end if
      if( option_set .eq. 2 ) then
          matprp(cyclic_start-1+1) = value
      else
       if( option_set .eq. 0 ) then
           call inmat_err_message( 2, 'gp_h_u', 6 )
        else if( option_set .eq. 1 ) then
           call inmat_err_message( 3, 'gp_h_u', 6 )
        end if
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                     input gp_tau for generalized hardening option  *
c *                                                                    *
c **********************************************************************
c
 420  continue
      if( .not. numr(value) ) then
         call errmsg(5,dum,'gp_tau',dumr,dumd)
         go to 210
      end if
      if( option_set .eq. 2 ) then
          matprp(cyclic_start-1+2) = value
      else
       if( option_set .eq. 0 ) then
           call inmat_err_message( 2, 'gp_tau', 6 )
        else if( option_set .eq. 1 ) then
           call inmat_err_message( 3, 'gp_tau', 6 )
        end if
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 input gp_delta_u for generalized hardening option  *
c *                                                                    *
c **********************************************************************
c
 430  continue
      if( .not. numr(value) ) then
         call errmsg(5,dum,'gp_delta_u',dumr,dumd)
         go to 210
      end if
      if( option_set .eq. 2 ) then
          matprp(cyclic_start-1+5) = value
      else
       if( option_set .eq. 0 ) then
           call inmat_err_message( 2, 'gp_delta_u', 9 )
        else if( option_set .eq. 1 ) then
           call inmat_err_message( 3, 'gp_delta_u', 9 )
        end if
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *               input gp_beta_u for generalized hardening option     *
c *                                                                    *
c **********************************************************************
c
 440  continue
      if( .not. numr(value) ) then
         call errmsg(5,dum,'gp_beta_u',dumr,dumd)
         go to 210
      end if
      if( option_set .eq. 2 ) then
          matprp(cyclic_start-1+3) = value
      else
       if( option_set .eq. 0 ) then
           call inmat_err_message( 2, 'gp_beta_u', 9 )
        else if( option_set .eq. 1 ) then
           call inmat_err_message( 3, 'gp_beta_u', 9 )
        end if
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 curves used to set temperature dependent props     *
c *                 for generalized_plasticity                         *
c *                    (1)  get the list of curves to use. expand list *
c *                         and save                                   *
c *                    (2)  set segmental flag true                    *
c *                                                                    *
c **********************************************************************
c
 460  continue
      call inseg_curve_list( curve_set_no, lmtprp(24) )
      if ( lmtprp(24) ) matprp(45) = curve_set_no
      go to 210

c
c
 9000 format(/,2x,'** values for cyclic plasticity model:',
     &  /,5x,'e:        ',f10.0,
     &  /,5x,'nu:       ',f10.3,
     &  /,5x,'rho:      ',f10.5,
     &  /,5x,'yld_pt:   ',f10.3,
     & /,5x,'cyclic properties...' )
 9100 format(5x,i2,f15.3)
      end
c     ****************************************************************
c     *                                                              *
c     *                  subroutine inmat_err_mess                   *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 5/16/2017                  *
c     *                                                              *
c     *     output various error messages during scan and checks     *
c     *     of material property values                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine inmat_err_message( ierrno, errstrng, nchars )
      use global_data ! old common.main
      implicit none
c
      integer :: ierrno, nchars
      character(len=*) :: errstrng
c
      character(len=50) :: erprms
c
c
      select case( ierrno )
c
      case( 1 )
         num_error = num_error + 1
         write(out,9001) errstrng(1:nchars)
         input_ok = .false.
c
      case( 2 )
         num_error = num_error + 1
         write(out,9002) errstrng(1:nchars)
         input_ok = .false.
c
      case( 3 )
         num_error = num_error + 1
         write(out,9003) errstrng(1:nchars)
         input_ok = .false.
c
      case( 4 )
         num_error = num_error + 1
         write(out,9004) errstrng(1:nchars)
         input_ok = .false.
c
      case( 5 )
         num_error = num_error + 1
         write(out,9005) errstrng(1:nchars)
         input_ok = .false.
c
      case( 6 )
         num_error = num_error + 1
         write(out,9006) errstrng(1:nchars)
         input_ok = .false.
c
      case( 7 )
         num_error = num_error + 1
         write(out,9007) errstrng(1:nchars)
         input_ok = .false.
c
      case default
        write(out,9999)
        stop
      end select
c
      return
c
 9001 format(/1x,'>>>>> error: invalid property for ',
     & 'generalized_plasticity option: ', a,
     & /14x,'This value ignored...')
 9002 format(/1x,'>>>>> error: no option set prior to defining',
     & /14x,'cyclic model properties. Must specifiy',
     & /14x,'nonlinear_hardening or generalized_plasticity before',
     & /14x,'giving property values. Property: ',a,' ignored...')
 9003 format(/1x,'>>>>> error: invalid property for ',
     & 'nonlinear_hardening option: ', a,
     & /14x,'This value ignored...')
 9004 format(/1x,'>>>>> error: the property: ',a,
     & /14x,'is not definable for the exp1_intf option',
     & /14x,'This value ignored...'/)
 9005 format(/1x,'>>>>> error: the property: ',a,
     & /14x,'is definable only for the ppr option',
     & /14x,'This value ignored...'/)
 9006 format(/1x,'>>>>> error: cohesive option: ',a,
     & /14x,'is not available...'/)
 9007 format(/1x,'>>>>> error: the property: ',a,
     & /14x,'is not definable for this cohesive option...'/)
c
 9999 format(/1x,'>>>>> Fatal Error: routine inmat_err_mess.',
     &   /16x,   'should have not reach this point.')
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine inmat_cohesive               *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/1/20154 rhd              *
c     *                                                              *
c     *     input properties for all the various options of the      *
c     *     cohesive material model (#4)                             *
c     *                                                              *
c     ****************************************************************
c
      subroutine inmat_cohesive( matprp, nonlocal_analysis )
c
      implicit none
c
c                       parameter declarations (the column of the
c                       table is passed for this material)
c
      real :: matprp(*)
      logical nonlocal_analysis
c
c                       local declarations
c
      integer dummy, dum, status, iout, idummy1, idummy2, cavit_loc
      logical  local_debug, true, matchs, matchs_exact, endcrd, numr,
     &         ltrue, lfalse, local_exp1_type, local_ppr_type,
     &         local_linear_type, local_fgm_type, local_cavit_type
      real  dumr, rword, value, fgm_mark, rtrue, rfalse
      double precision
     &   dumd
      equivalence (lfalse, rfalse), (rtrue, ltrue)
      character(len=1) :: dums
      data fgm_mark / -99.0 /
c
c           enter here having seen the properties and
c           cohesive keywords at the beginning of a line.
c           we leave at the end of a logical line having
c           processed all property values for the
c           cohesive-type plasticty model.
c
c           the model supports various options for the cohesive
c           formulation
c                 1 - linear elastic;
c                 2 - bilinear
c                 3 - ramp;
c                 4 - exponential_1;
c                 5 - exponential_2;
c                 6 - PPR (Park, Paulion, Roesler)
c                 7 - nonlocal creep cavtiation
c           At this time, only options 1, 4, 6, 7 are defined
c           FGM option has been removed
c           We still scan and store properties for the yet to be implemented
c
c
c      Changes to support major overhaul the the model
c      in summer 2010. *Note* some properties are defined/translated from input
c      for options not yet implemented.
c
c      The current options are implemented:
c        1: linear, 3: exponential, 6: PPR, 7: creep cavitation
c        The ramp, exponential_2 and bilinear options are not implemented
c        even though most of the input translation and property storage is
c        done here.
c
c           assign default material values. the current
c           ordering of material values is:
c
c           9  -- material model type 9 (set in inmat)
c                  = 4 Cohesive zone models: linear elastic,
c                      bilinear, ramp, exponential_1,
c                      exponential_2, PPR, cavitation
c          13  -- debug material model computations
c          22  -- allow material model to cut step due to
c          23  -- flag to allow crack growth element killing
c                 excessive reversed plasticity
c          32  -- linear option: interface stiffness in
c                 the longitudinal direction
c          33  -- linear option: interface stiffness in
c                 the transverse direction
c          34  -- linear option: interface stiffness in
c                 the normal direction
c          35  -- maximum normal stress of the traction-separation curve
c          36  -- maximum shear stress of the traction-separation curve
c          37  -- shape parameter for bilinear and ramp
c          38  -- second (additional) shape parameter for ramp
c          39  -- separation distance in sliding at peak stress
c          40  -- separation distance in opening at peak stress
c          41  -- "equivalent" separation distance at peak stress --
c                 makes normal + shear into single term via beta_coh(esive)
c                 value (exp1 only)
c          42  -- For exp1: "beta" - a ratio to determine the equivalent
c                           separation. Set =0 for Mode I only behavior.
c          43  -- a flag for identifying the interface element (logical)
c          44  -- type of interface models
c                 1 - linear elastic; 2- bilinear
c                 3 - ramp; 4 - exponential_1; 5 - exponential_2;
c                 6 - PPR, 7 - creep cavitation
c      (*) 46  -- ductile material volume fracture for fgm
c          47  -- = not used now. was used for fgm
c          48  -- separation distance at peak stress for ductile (fgm)
c          49  -- separation distance brittle (fgm)
c          50  -- maximum normal stress on curve ductile (fgm)
c          51  -- maximum shear stress on curve brittle (fgm)
c          52  -- beta_ductile (fgm) [beta cohesive]
c          53  -- beta_brittle (fgm) [beta cohesive]
c          54  -- compression stiffness multiplier for
c                 cohesive materials
c          90-- PPR cohesive model
c                  90 - sig_max
c                  91 - tau_max
c                  92 - G_normal
c                  93 - G_shear
c                  94 - ratio_normal
c                  95 - ratio_shear
c                  96 - shape_normal
c                  97 - shape_shear
c                  98 - debug model computations
c
c
c         see inmat for more details on variables in matprp vector
c         See code below that sets cavity property location
c
      local_debug = .false.
      ltrue       = .true.
      lfalse      = .false.
      call iodevn( idummy1, iout, idummy2, 1 )
c
      cavit_loc = 200   ! cavit props start @ loc +1 row
c
      if( local_debug ) write(iout,*) '.... entered inmat_cohesive ....'
      matprp(13) = rfalse
      matprp(22) = rfalse
      matprp(23) = rfalse
      matprp(44) = 1.0
      matprp(32) = 0.0
      matprp(33) = 0.0
      matprp(34) = 0.0
      matprp(35) = 0.0
      matprp(36) = 0.0
      matprp(37) = 0.0
      matprp(38) = 0.0
      matprp(39) = 0.0
      matprp(40) = 0.0
      matprp(41) = 0.0
      matprp(42) = 0.0
      matprp(43) = rtrue
      matprp(46) = 0.0
      matprp(47) = 0.0
      matprp(48) = 0.0
      matprp(49) = 0.0
      matprp(50) = 0.0
      matprp(51) = 0.0
      matprp(52) = 0.0
      matprp(53) = 0.0
      matprp(54) = 10.0
      matprp(90) = 0.0
      matprp(91) = 0.0
      matprp(92) = 0.0
      matprp(93) = 0.0
      matprp(94) = 0.0
      matprp(95) = 0.0
      matprp(96) = 0.0
      matprp(97) = 0.0
c
      local_linear_type = .false.
      local_exp1_type   = .false.
      local_ppr_type    = .false.
      local_fgm_type    = .false.
      local_cavit_type  = .false.
c
c
      call reset
      if( true( dummy ) ) call splunj
      if( matchs('properties',10) ) call splunj
      if( matchs_exact('nonlocal') ) nonlocal_analysis = .true.
      if( matchs('cohesive',6) ) call splunj
c      if( local_debug ) write(iout,*) '  @ 1'
      go to 210
c
 206  continue
      call readsc

 210  continue
c      if( local_debug ) write(iout,*) '  @ 2'
      if ( matchs_exact( 'type' )           ) go to 210
      if ( matchs_exact('exp1_intf')        ) go to 500
      if ( matchs_exact('exp2_intf')        ) go to 510
      if ( matchs_exact('beta_ductile')     ) go to 680
      if ( matchs_exact('beta_brittle')     ) go to 690
      if ( matchs_exact('beta_coh')         ) go to 620
      if ( matchs_exact('beta')             ) go to 620
      if ( matchs_exact('ramp_intf')        ) go to 490
      if ( matchs_exact('bilinear_intf')    ) go to 480
      if ( matchs_exact('delta_crit_ductile') ) go to 640
      if ( matchs_exact('delta_peak_ductile') ) go to 640
      if ( matchs_exact('delta_crit_brittle') ) go to 650
      if ( matchs_exact('delta_peak_brittle') ) go to 650
      if ( matchs_exact('deltat')           ) go to 590
      if ( matchs_exact('deltan')           ) go to 600
      if ( matchs_exact('delta_crit')       ) go to 610
      if ( matchs_exact('delta_peak')       ) go to 610
      if ( matchs_exact('killable')         ) go to 450
      if ( matchs_exact('stifft')           ) go to 520
      if ( matchs_exact('stiffs')           ) go to 530
      if ( matchs_exact('stiffn')           ) go to 540
      if ( matchs_exact('sig_max_ductile') ) go to 660
      if ( matchs_exact('sig_max_brittle') ) go to 670
      if ( matchs_exact('sig_max')          ) go to 550
      if ( matchs_exact('sig_peak')          ) go to 550
      if ( matchs_exact('tau_max')           ) go to 560
      if ( matchs_exact('tau_peak')          ) go to 560
      if ( matchs_exact('linear_intf')      ) go to 470
      if ( matchs_exact('lambda1')          ) go to 570
      if ( matchs_exact('lambda2')          ) go to 580
      if ( matchs_exact('vol_fract_ductile') ) go to 630
      if ( matchs_exact('compression_multiplier') ) go to 700
      if ( matchs_exact('PPR')         ) go to 800
      if ( matchs_exact('G_normal')         ) go to 805
      if ( matchs_exact('G_shear')          ) go to 810
      if ( matchs_exact('ratio_normal')    ) go to 820
      if ( matchs_exact('ratio_shear')     ) go to 830
      if ( matchs_exact('shape_normal')    ) go to 840
      if ( matchs_exact('shape_shear')     ) go to 850
      if ( matchs_exact('debug')            ) go to 860
      if ( matchs('cavitation',5)           ) go to 870
      if ( matchs_exact('lambda_1')         ) go to 880
      if ( matchs_exact('lambda_2')         ) go to 890
      if ( matchs_exact('lambda_3')         ) go to 900
      if ( matchs_exact('lambda_4')         ) go to 910
      if ( matchs_exact('lambda_5')         ) go to 920
      if ( matchs_exact('lambda_6')         ) go to 930
      if ( matchs_exact('delta_star')       ) go to 940
      if ( matchs_exact('sigma_star')       ) go to 950
      if ( matchs_exact('lambda_7')         ) go to 960
      if ( matchs_exact('eta_b')            ) go to 965
      if ( matchs_exact(' r_i')              ) go to 970 ! see note
      if ( matchs_exact('sigma_0')          ) go to 975
      if ( matchs_exact('t_c')              ) go to 980
      if ( matchs_exact('d')                ) go to 985
      if ( matchs_exact('a_0')              ) go to 990
      if ( matchs_exact('b_0')              ) go to 995
      if ( matchs_exact('psi_angle')        ) go to 1000
      if ( matchs_exact(' n_i')             ) go to 1005 ! see note
      if ( matchs_exact('f_n')              ) go to 1010
      if ( matchs_exact('n')                ) go to 1015
      if ( matchs_exact('n_max')            ) go to 1020
      if ( matchs_exact('local_debug')      ) go to 1900
      if ( matchs(',',1)                    ) go to 215
c
c                       change made on 3/1/2015. R_I, N_I no longer
c                       allowed on input.
c
c
c                       there is no match with existing properties.
c                       check for end of card. if not, print error
c                       message, skip over bad item and keep going.
c
      if( endcrd(dum) ) then
         if( local_ppr_type ) call inmat_ppr_print( matprp, iout )
         if( local_debug ) then
           write(iout,9000) matprp(44)
           write(iout,9002) local_linear_type, local_exp1_type,
     &                   local_ppr_type, local_cavit_type
           if( local_linear_type ) then
             write(iout,9100) matprp(32)
             write(iout,9102) matprp(33)
             write(iout,9104) matprp(34)
           end if
           if( .not. local_cavit_type ) then
             write(iout,9106) matprp(35)
             write(iout,9108) matprp(36)
             write(iout,9110) matprp(37)
             write(iout,9112) matprp(38)
             write(iout,9114) matprp(39)
             write(iout,9116) matprp(40)
             write(iout,9118) matprp(41)
             write(iout,9120) matprp(42)
             write(iout,9121) matprp(43)
             write(iout,9138) matprp(54)
             write(iout,9140) matprp(90)
             if( local_ppr_type ) then
               write(iout,9142) matprp(91)
               write(iout,9144) matprp(92)
               write(iout,9146) matprp(93)
               write(iout,9148) matprp(94)
               write(iout,9150) matprp(95)
               write(iout,9152) matprp(96)
               write(iout,9154) matprp(97)
             end if
           end if
           if( local_cavit_type ) then
             write(iout,9202) matprp(cavit_loc+1)
             write(iout,9204) matprp(cavit_loc+2)
             write(iout,9206) matprp(cavit_loc+3)
             write(iout,9208) matprp(cavit_loc+4)
             write(iout,9210) matprp(cavit_loc+5)
             write(iout,9212) matprp(cavit_loc+6)
             write(iout,9214) matprp(cavit_loc+7)
             write(iout,9216) matprp(cavit_loc+8)
             write(iout,9218) matprp(cavit_loc+9)
             write(iout,9220) matprp(cavit_loc+10)
c             write(iout,9222) matprp(cavit_loc+11)
             write(iout,9224) matprp(cavit_loc+12)
             write(iout,9226) matprp(cavit_loc+13)
             write(iout,9228) matprp(cavit_loc+14)
             write(iout,9230) matprp(cavit_loc+15)
             write(iout,9232) matprp(cavit_loc+16)
             write(iout,9234) matprp(cavit_loc+17)
c             write(iout,9236) matprp(cavit_loc+18)
             write(iout,9238) matprp(cavit_loc+19)
             write(iout,9240) matprp(cavit_loc+20)
             write(iout,9242) matprp(cavit_loc+21)
           end if
           write(iout,9156) matprp(98)
         end if
         return
      else
         call errmsg(4,dum,dums,dumr,dumd)
         if( true(dum) ) go to 210
      end if
c
c                       comma separator. if at the end of a line,
c                       denotes continuation. otherwise, ignore.
c
 215  continue
      if( endcrd(dum) ) then
         go to 206
      else
         go to 210
      end if
c
c
c **********************************************************************
c *                                                                    *
c *                 whether to allow the element to be killed          *
c *                 in a crack growth problem                          *
c *                                                                    *
c **********************************************************************
c
 450  continue
      matprp(23) = rtrue
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                 cohesive zone : linear elastic interface           *
c *                                                                    *
c **********************************************************************
c
c
 470  continue
      matprp(43) = rtrue
      matprp(44) = 1.0
      matprp(54) = 10.0
      local_linear_type = .true.
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                 cohesive zone : bilinear interface                 *
c *                      ** not implemented **                         *
c *                                                                    *
c **********************************************************************
c
c
 480  continue
      matprp(43) = rtrue
      matprp(44) = 2.0
      call inmat_err_message(6,'bilinear_intf',13)
      go to 210

c
c **********************************************************************
c *                                                                    *
c *                 cohesive zone : ramp curve interface               *
c *                      ** not implemented **                         *
c *                                                                    *
c **********************************************************************
c
c
 490  continue
      matprp(43) = rtrue
      matprp(44) = 3.0
      call inmat_err_message(6,'ramp_intf',9)
      go to 210

c
c **********************************************************************
c *                                                                    *
c *                 cohesive zone : exponential type 1                 *
c *                                                                    *
c **********************************************************************
c
c
 500  continue
      matprp(43) = rtrue
      matprp(44) = 4.0
      matprp(54) = 10.0
      local_exp1_type = .true.
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                 cohesivce zone : exponential type 2                *
c *                      ** not implemented **                         *
c *                                                                    *
c **********************************************************************
c
c
 510  continue
      matprp(43) = rtrue
      matprp(44) = 5.0
      call inmat_err_message(6,'exp2_intf',9)
      go to 210
c
c **********************************************************************
c *                                                                    *
c *                 interface stiffness in longitudinal direction       *
c *                      linear elastic interface behaviour            *
c *                                                                    *
c **********************************************************************
c
c
 520  continue
      if( .not. numr(matprp(32)) ) then
         call errmsg(5,dum,'stifft',dumr,dumd)
      end if
      if( .not. local_linear_type )
     &  call inmat_err_message(7,'stifft',6)
      go to 210               !   1234567890123
c
c
 530  continue
      if( .not. numr(matprp(33)) ) then
         call errmsg(5,dum,'stiffs',dumr,dumd)
      end if
      if( .not. local_linear_type )
     &  call inmat_err_message(7,'stiffs',6)
      go to 210               !   1234567890123
c
c
c
 540  continue
      if( .not. numr(matprp(34)) ) then
         call errmsg(5,dum,'stiffn',dumr,dumd)
      end if
      if( .not. local_linear_type )
     &  call inmat_err_message(7,'stiffn',6)
      go to 210               !   1234567890123
c
c
c **********************************************************************
c *                                                                    *
c *                 peak normal traction on curve                      *
c *  save in both locations for both regular CZM and PPR               *
c *                                                                    *
c **********************************************************************
c
c
 550  continue
      if( .not. numr(matprp(35)) ) then
         call errmsg(5,dum,'sig_peak',dumr,dumd)
      end if
      matprp(90) = matprp(35)
      go to 210
c
c
c **********************************************************************
c *                                                                    *
c *                 peak shear traction on curve                       *
c *  save in both locations for both regular CZM and PPR               *
c *                                                                    *
c **********************************************************************
c
c
 560  continue
      if( .not. numr(matprp(36)) ) then
         call errmsg(5,dum,'tau_peak',dumr,dumd)
      end if
      if( .not. local_ppr_type )
     &  call inmat_err_message(7,'tau_peak',8)
      matprp(91) = matprp(36)
      go to 210               !   1234567890123
      go to 210
c
 570  continue
      if( .not. numr(matprp(37)) ) then
         call errmsg(5,dum,'lambda1',dumr,dumd)
      end if
      go to 210
c
c
 580  continue
      if( .not. numr(matprp(38)) ) then
         call errmsg(5,dum,'lambda2',dumr,dumd)
      end if
      go to 210
c
c
 590  continue
      if( .not. numr(matprp(39)) ) then
         call errmsg(5,dum,'deltat',dumr,dumd)
      end if
      go to 210
c
c
 600  continue
      if(.not.numr(matprp(40))) then
         call errmsg(5,dum,'deltan',dumr,dumd)
      end if
      go to 210
c
c
 610  continue
      if( .not. numr(matprp(41)) ) then
         call errmsg(5,dum,'delta_peak',dumr,dumd)
      end if
      go to 210
c
c
 620  continue
      if( .not. numr(matprp(42)) ) then
         call errmsg(5,dum,'beta',dumr,dumd)
      end if
      if( .not. local_exp1_type )
     &  call inmat_err_message(7,'beta',4)
      go to 210               !   1234567890123
c
c
 630  continue
c
      call inmat_fgm( status, 'volume_fraction'  )
      if ( status  .eq. 0 ) then
          matprp(46) = fgm_mark
          matprp(47) = 1.0
          go to 210
      elseif ( status .eq. 1 ) then
          go to 210
      end if
c
      if( .not. numr(matprp(46)) ) then
         call errmsg(5,dum,'vol_fraction',dumr,dumd)
      end if
      matprp(47) = 1.0
      go to 210
c
c
 640  continue
      if( .not. numr(matprp(48)) ) then
         call errmsg(5,dum,'delta_crit_ductile',dumr,dumd)
      end if
      matprp(47) = 1.0
      go to 210
c
c
 650  continue
      if( .not. numr(matprp(49)) ) then
         call errmsg(5,dum,'delta_crit_brittle',dumr,dumd)
      end if
      matprp(47) = 1.0
      go to 210
c
c
 660  continue
      if( .not. numr(matprp(50)) ) then
         call errmsg(5,dum,'critial_stress_ductile',dumr,dumd)
      end if
      matprp(47) = 1.0
      go to 210
c
c
 670  continue
      if( .not. numr(matprp(51)) ) then
         call errmsg(5,dum,'critial_stress_brittle',dumr,dumd)
      end if
      matprp(47) = 1.0
      go to 210
c
c
 680  continue
      if( .not. numr(matprp(52)) ) then
         call errmsg(5,dum,'beta_ductile',dumr,dumd)
      end if
      matprp(47) = 1.0
      go to 210
c
c
 690  continue
      if( .not. numr(matprp(53)) ) then
         call errmsg(5,dum,'beta_brittle',dumr,dumd)
      end if
      matprp(47) = 1.0
      go to 210
c
c
 700  continue
      if( .not. numr(matprp(54)) ) then
         call errmsg(5,dum,'compression_multiplier',dumr,dumd)
      end if
      go to 210
c
c **********************************************************************
c *                                                                    *
c *         properties for PPR type cohesive material                  *
c *                                                                    *
c **********************************************************************
c
 800  continue
      matprp(42) = 0.0
      matprp(43) = rtrue
      matprp(44) = 6.0
      matprp(47) = 0.0
      matprp(54) = 2.0
      matprp(98) = rfalse
      local_ppr_type = .true.
      go to 210
c
 805  continue
      if( .not. numr(matprp(92)) ) then
         call errmsg(5,dum,'G_normal',dumr,dumd)
      end if
      if( .not. local_ppr_type )
     &  call inmat_err_message(5,'G_normal',8)
      go to 210
c
 810  continue
      if( .not. numr(matprp(93)) ) then
         call errmsg(5,dum,'G_shear',dumr,dumd)
      end if
      if( .not. local_ppr_type )
     &  call inmat_err_message(5,'G_shear',7)
      go to 210
c
 820  continue
      if( .not. numr(matprp(94)) ) then
         call errmsg(5,dum,'ratio_normal',dumr,dumd)
      end if
      if( .not. local_ppr_type )
     &  call inmat_err_message(5,'ratio_normal',12)
      go to 210
c
 830  continue
      if( .not. numr(matprp(95)) ) then
         call errmsg(5,dum,'ratio_shear',dumr,dumd)
      end if
      if( .not. local_ppr_type )
     &  call inmat_err_message(5,'ratio_shear',11)
      go to 210
c
 840  continue
      if( .not. numr(matprp(96)) ) then
         call errmsg(5,dum,'shape_normal',dumr,dumd)
      end if
      if( .not. local_ppr_type )
     &  call inmat_err_message(5,'shape_normal',12)
      go to 210
c
 850  continue
      if( .not. numr(matprp(97)) ) then
         call errmsg(5,dum,'shape_shear',dumr,dumd)
      end if
      if( .not. local_ppr_type )
     &  call inmat_err_message(5,'shape_shear',11)
      go to 210
c
 860  continue
      matprp(13) = rtrue
      matprp(98) = rtrue
      go to 210
c
c **********************************************************************
c *                                                                    *
c *           properties for the creep cavitation model                *
c *                                                                    *
c **********************************************************************
c
c
 870  continue
      matprp(43) = rtrue
      matprp(44) = 7.0
      matprp(47) = 0.0
      matprp(54) = 10.0
      matprp(98) = rfalse
      local_cavit_type = .true.
      matprp(cavit_loc+1:cavit_loc+21) = -1.0
      go to 210
c
 880  continue
      if( .not. numr(matprp(cavit_loc+1)) ) then
         call errmsg(5,dum,'lambda_1',dumr,dumd)
      end if
      go to 210
 890  continue
      if( .not. numr(matprp(cavit_loc+2)) ) then
         call errmsg(5,dum,'lambda_2',dumr,dumd)
      end if
      go to 210
 900  continue
      if( .not. numr(matprp(cavit_loc+3)) ) then
         call errmsg(5,dum,'lambda_3',dumr,dumd)
      end if
      go to 210
 910  continue
      if( .not. numr(matprp(cavit_loc+4)) ) then
         call errmsg(5,dum,'lambda_4',dumr,dumd)
      end if
      go to 210
 920  continue
      if( .not. numr(matprp(cavit_loc+5)) ) then
         call errmsg(5,dum,'lambda_5',dumr,dumd)
      end if
      go to 210
 930  continue
      if( .not. numr(matprp(cavit_loc+6)) ) then
         call errmsg(5,dum,'lambda_6',dumr,dumd)
      end if
      go to 210
 940  continue
      if( .not. numr(matprp(cavit_loc+7)) ) then
         call errmsg(5,dum,'delta_star',dumr,dumd)
      end if
      go to 210
 950  continue
      if( .not. numr(matprp(cavit_loc+8)) ) then
         call errmsg(5,dum,'sigma_star',dumr,dumd)
      end if
 960  continue
      if( .not. numr(matprp(cavit_loc+9)) ) then
         call errmsg(5,dum,'lambda_7',dumr,dumd)
      end if
      go to 210
 965  continue
      if( .not. numr(matprp(cavit_loc+10)) ) then
         call errmsg(5,dum,'eta_b',dumr,dumd)
      end if
      go to 210
 970  continue ! available for re-use. R_I not allowed on
c                 input
      if( .not. numr(matprp(cavit_loc+11)) ) then
         call errmsg(5,dum,'r_i',dumr,dumd)
      end if
      go to 210
 975  continue
      if( .not. numr(matprp(cavit_loc+12)) ) then
         call errmsg(5,dum,'sigma_0',dumr,dumd)
      end if
      go to 210
 980  continue
      if( .not. numr(matprp(cavit_loc+13)) ) then
         call errmsg(5,dum,'t_c',dumr,dumd)
      end if
      go to 210
 985  continue
      if( .not. numr(matprp(cavit_loc+14)) ) then
         call errmsg(5,dum,'d',dumr,dumd)
      end if
      go to 210
 990  continue
      if( .not. numr(matprp(cavit_loc+15)) ) then
         call errmsg(5,dum,'a_0',dumr,dumd)
      end if
      go to 210
 995  continue
      if( .not. numr(matprp(cavit_loc+16)) ) then
         call errmsg(5,dum,'b_0',dumr,dumd)
      end if
      go to 210
1000  continue
      if( .not. numr(matprp(cavit_loc+17)) ) then
         call errmsg(5,dum,'psi_angle',dumr,dumd)
      end if
      go to 210
1005  continue  ! available for re-use. N_I not allowed on
c                 input
      if( .not. numr(matprp(cavit_loc+18)) ) then
         call errmsg(5,dum,'n_i',dumr,dumd)
      end if
      go to 210
1010  continue
      if( .not. numr(matprp(cavit_loc+19)) ) then
         call errmsg(5,dum,'f_n',dumr,dumd)
      end if
      go to 210
1015  continue
      if( .not. numr(matprp(cavit_loc+20)) ) then
         call errmsg(5,dum,'n',dumr,dumd)
      end if
      go to 210
1020  continue
      if( .not. numr(matprp(cavit_loc+21)) ) then
         call errmsg(5,dum,'n_max',dumr,dumd)
      end if
      go to 210
c
c
 1900 continue
      local_debug = .true.
      go to 210

c
 9000 format(/,2x,'** values for cohesive plasticity model type:',
     &    f4.1,' **'/)
 9002 format(2x,'local_linear_type: ',l1,
     &      /2x,'local_exp1_type:   ',l1,
     &      /2x,'local_ppr_type:    ',l1,
     &      /2x,'local_cavit_type:  ',l1)
 9100 format(2x,'interface stiffness in the longitudinal direction: ',
     &    f10.2 )
 9102 format(2x,'interface stiffness in the transverse direction: ',
     &    f10.2 )
 9104 format(2x,'interface stiffness in the normal direction: ',
     &    f10.2 )
 9106 format(2x,'peak normal stress of the interface: ',f10.2 )
 9108 format(2x,'peak shear stress of the interface : ',f10.2 )
 9110 format(2x,'shape parameter for bilinear and ramp: ',f10.2 )
 9112 format(2x,'second (additional) shape parameter for ramp: ',
     &    f10.2 )
 9114 format(2x,'separation distance @ peak in sliding: ',f10.6 )
 9116 format(2x,'separation distance @ peak in opening: ',f10.6 )
 9118 format(2x,'equivalent separation distance @ peak: ',f10.6 )
 9120 format(2x,'for exp1 option: a ratio to determine the ',
     &          'equivalent separation: ',
     &     /,18x,' under mixed mode loading ( =0 => mode I ) ',
     &     /2x,'value= ', f10.3)
 9121 format(2x,'a flag for identifying the interface element: ',l1 )
 9138 format(2x,'compression stiffness multiplier: ',f10.2)
 9140 format(2x,'sig_max(sig_peak): ',f10.2)
 9142 format(2x,'tau_max(tau_peak): ',f10.2)
 9144 format(2x,'G_normal: ',f10.2)
 9146 format(2x,'G_shear:  ',f10.2)
 9148 format(2x,'ratio_normal: ',f10.2)
 9150 format(2x,'ratio_shear:  ',f10.2)
 9152 format(2x,'shape_normal: ',f10.2)
 9154 format(2x,'shape_shear:  ',f10.2)
 9156 format(2x,'debug model computations: ',l1)
 9202 format(2x,'lambda_1:     ',f15.6)
 9204 format(2x,'lambda_2:     ',f15.6)
 9206 format(2x,'lambda_3:     ',f15.6)
 9208 format(2x,'lambda_4:     ',f15.6)
 9210 format(2x,'lambda_5:     ',f15.6)
 9212 format(2x,'lambda_6:     ',f15.6)
 9214 format(2x,'delta_star:   ',f15.6)
 9216 format(2x,'sigma_star:   ',f15.6)
 9218 format(2x,'lambda_7:     ',f15.6)
 9220 format(2x,'eta_b:        ',f15.6)
 9222 format(2x,'r_i:          ',f15.6)
 9224 format(2x,'sigma_0:      ',f15.6)
 9226 format(2x,'t_c:          ',f15.6)
 9228 format(2x,'d:            ',f15.6)
 9230 format(2x,'a_0:          ',f15.6)
 9232 format(2x,'b_0:          ',f15.6)
 9234 format(2x,'psi_angle:    ',f15.6)
 9236 format(2x,'n_i:          ',f15.6)
 9238 format(2x,'f_n:          ',f15.6)
 9240 format(2x,'n:            ',f15.6)
 9242 format(2x,'n_max:        ',f15.6)
c
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine inmat_ppr_print              *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 11/28/2010                 *
c     *                                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine inmat_ppr_print( matprp, iout )
      implicit none
c
      integer :: iout
      real :: matprp(*)
      double precision
     &   alph, beta, ln, lt, m, n, Gn, Gt, Tn_m, Tt_m, dn, dt, one
      data one
     &    / 1.0d00 /
c
      alph   = matprp(96)
      beta   = matprp(97)
      ln     = matprp(94)
      lt     = matprp(95)
      Gn     = matprp(92)
      Gt     = matprp(93)
      Tn_m   = matprp(90)
      Tt_m   = matprp(91)
c
      m = (alph-one)*alph*ln**2/(one-alph*ln**2)
      n = (beta-one)*beta*lt**2/(one-beta*lt**2)
      dn = alph*Gn/(m*Tn_m)*(one-ln)**(alph-one)
     &     * (alph/m*ln+one)**(m-one)*(alph+m)*ln
      dt = beta*Gt/(n*Tt_m)*(one-lt)**(beta-one)
     &     * (beta/n*lt+one)**(n-one)*(beta+n)*lt

      write(iout,9000) dn, dt, dn*ln, dt*lt
      return
 9000 format(/,
     &   '>>>>      information from ppr option. computed values: ',/,
     &   '                delta_n: ',f12.9,'        delta_t: ',f12.9,/,
     &   '          delta_n(peak): ',f12.9,'  delta_t(peak): ',f12.9 / )
      end
















