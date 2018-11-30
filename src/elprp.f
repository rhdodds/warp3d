c
c     Catalogue of Element Numbers (types) in WARP3D
c     ----------------------------------------------
c
c     Type    Name        Description
c     ----    ----        -----------
c    1        q3disop      20 node hex
c    2        l3disop       8 nodes hex
c    3        ts12isop     12 node hex transition
c    4        ts15isop     15 node hex transition
c    5        ts9isop       9 node hex transition
c    6        tet10        10 node tetrahedron
c    7        wedge15      15 node wedge element
c    8        tri6          6 node 2-D triangle
c    9        quad8         8 node 2-D quadralateral
c    10       axiquad8      8 node axisymmetric quadralateral
c    11       axitri6       6 node triangular, axisymmetric
c    12       inter_8       8 node quadralateral interface
c    13       tet_4         4 node tetrahedron
c    14       trint6        6 node triangular interface
c    15       trint12      12 node triangular interface
c    16       ***           reserved to use in obtaining
c                           integration points, derivatives,
c                           shape functions
c    17       ***           reserved to use in obtaining
c                           integration points, derivatives,
c                           shape functions
c    18       bar2          2 node bar
c    19       link2         2 node link element
c
c
c         initst routine must be updated when new elements are added
c         to modify the logical vectors that flag certain
c         characteristics, e.g. is it a cohesive-interface element
c         the vectors are stored in mod_main.f
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp18                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 8/11/2017 rhd              *
c     *                                                              *
c     * put the properties of the bar2 element into global storage   *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine elprp18( elem, type )
c
      use global_data ! old common.main
      use main_data
      use segmental_curves
c
      implicit none
c
      integer :: elem, type
c
      integer :: intord, outloc, iloc, geonl, matnum, iword
      character(len=1) :: dums
      real :: e, et, h,  dumr, rword
      real, parameter :: rzero = 0.0, one = 1.0
      double precision :: dumd
      logical :: fgm
      equivalence( iword, rword )
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       see routine inamt for definition of rows of
c                       matprp table.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 2             ! number nodes
      iprops(3,elem)  = 24*two16+6    ! same as solid nodes. not used in code
      iprops(4,elem)  = 3             ! number dof each element node
      iprops(11,elem) = 6             ! not used in code
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
c
      iprops(42,elem) = elstor(11,elem)
c
c
c                       store data for the order of integration.
c                       the default value - only 1 point for bar
c                       for spring.
c
      intord = elstor(3,elem)
      iprops(5,elem) = 1   ! overrides any user value
      iprops(6,elem) = 1   !    ""
c
c                       store data for output location
c                       the default value is gausspts. other options
c                       are element or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       store data concerning the output format for
c                       output. the default value is short.
c
      lprops(16,elem) = .true.
c
c                       store data concerning the geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c                       store data concerning the bbar flag.
c                       the default value is off.
c
      lprops(19,elem) = .false.
C
C                       STORE DATA CONCERNING MATERIAL PROPERTIES.
C                       yes, this wastes space, but we maintain
c                       consistency in that all data for an element
c                       is here.
C                       Row            Quantity
c                       ---            --------
c                        7           young's modulus
c                        8           poisson's ratio. set = 0 for bar
c                        9           thermal expansion coeff.
c                                    (alpha for isotropic, alphax for
c                                     anisotropic)
c                        10          mass density
c                        13          thermal alphay for anisotropic
c                        14          kinematic hardening parameter
c                                      (beta)
c                        15          h-prime for linear hardening
c                        17          request for linear material
c                                    formulation
c                        20          viscopls m-power
c                        21          power-law hardening n-power or
c                                    set number for segmental
c                                    stress-strain curves
c                        22          viscoplastic ref strain rate
c                                    ( the user inputs D which the inverse
c                                      of eta, here we change D to eta
c                                      before passing the value to
c                                      props(22,elem) )
c                        23          uniaxial yield stress (sig-0)
c                        24          bit 1:  request for debug of material
c                                            computations
c                                    bit 2:  allow material model to
c                                            cutback adaptive solution
c                                            due to reversed plasticity
c                                    bit 3:  set of segmental stress-strain
c                                            curves specified for element's
c                                            material
c                        25          material model type
c                        26          initial porosity for gurson model
c                        27          q1 for gurson model
c                        28          q2 for gurson model
c                        29          q3 for gurson model
c                        30          bit 1:  logical nucleation flag for
c                                            gurson model
c                                    bit 2:  flag if element is
c                                            killable
c                        31          s_n parameter for gurson model
c                        32          eps_n parameter for gurson model
c                        33          f_n parameter for gurson model
c                        34          thermal alphaz for anisotropic
c                        35          thermal alphaxy for anisotropic
c                        36          thermal alphayz for anisotropic
c                        37          thermal alphaxz for anisotropic
c                        38          material number during input
c                        43          area for bar elements
c                        44-46       stiffness values for link elements
c
c
c
      matnum          = elstor(2,elem)
      iprops(38,elem) = matnum
      e               = matprp(1,matnum)
      et              = matprp(4,matnum)
      h               = (e*et)/(e-et)
      props(7,elem)   = e
      props(8,elem)   = rzero
c
c                       set thermal expansion coeeficients
c                       isotropic only  -- default
c
      props(9,elem)  = matprp(6,matnum)
      props(13,elem) = matprp(6,matnum)
      props(34,elem) = matprp(6,matnum)
      props(35,elem) = rzero
      props(36,elem) = rzero
      props(37,elem) = rzero
c
      props(10,elem)  = matprp(7,matnum)
c
      props(14,elem)  = matprp(3,matnum)
      props(15,elem)  = h
      lprops(17,elem) = lmtprp(8,matnum)
      props(20,elem)  = matprp(10,matnum)
      props(21,elem)  = matprp(11,matnum)
      props(22,elem)  = rzero
      if ( matprp(12,matnum) .ne. rzero )
     &  props(22,elem) = one / matprp(12,matnum)
      props(23,elem)  = 1.0e10 !  matprp(5,matnum) for linear elastic bar
      iprops(24,elem)  = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 ) ! bit 0
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 ) ! bit 1
      if ( lmtprp(24,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 4 ) ! bit 2
      iprops(25,elem) = matprp(9,matnum)
      props(26,elem)  = rzero
      props(27,elem)  = rzero
      props(28,elem)  = rzero
      props(29,elem)  = rzero
      iprops(30,elem) = 0
      props(31,elem)  = rzero
      props(32,elem)  = rzero
      props(33,elem)  = rzero
      props(34,elem)  = rzero
      props(35,elem)  = rzero
      props(36,elem)  = rzero
      props(37,elem)  = rzero
      iword           = elstor(12,elem)
      props(43,elem)  = rword ! element area
c
c                       if the material associated with the element has
c                       a set of stress-strain curves, store the set number
c                       for the curves and set the
c                       yield stress for the element material using the
c                       first curve in the set.
c
      call elem_set_segmental( elem, matnum )
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp19                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 8/20/2017 rhd              *
c     *                                                              *
c     * put the properties of the link2 element into global storage  *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine elprp19( elem, type )
c
      use global_data ! old common.main
      use main_data
      use segmental_curves
c
      implicit none
c
      integer :: elem, type
c
      integer :: intord, outloc, iloc, geonl, matnum, iword
      character(len=1) :: dums
      real :: e, et, h,  dumr, rword
      real, parameter :: rzero = 0.0, one = 1.0, rbig = -1.0e20
      double precision :: dumd
      logical :: fgm
      equivalence( iword, rword )
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       see routine inamt for definition of rows of
c                       matprp table.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 2             ! number nodes
      iprops(3,elem)  = 24*two16+6    ! same as solid nodes. not used in code
      iprops(4,elem)  = 3             ! number dof each element node
      iprops(11,elem) = 6             ! not used in code
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
c
      iprops(42,elem) = elstor(11,elem)
c
c
c                       store data for the order of integration.
c                       the default value - only 1 point for bar
c                       for spring.
c
      intord = elstor(3,elem)
      iprops(5,elem) = 1   ! overrides any user value
      iprops(6,elem) = 1   !    ""
c
c                       store data for output location. gausspts,
c                       nodes, center are all same
c
      outloc = elstor(4,elem)
      iloc   = 1
      iprops(12,elem) = iloc
c
c                       store data concerning the output format for
c                       output. the default value is short.
c
      lprops(16,elem) = .true.
c
c                       store data concerning the geometric nonlinearity flag.
c                       the default value is on.
c                       link2 elements have small displacement only
c                       formulation.
c
      geonl= elstor(7,elem)
      lprops(18,elem) = .false.

c
c                       store data concerning the bbar flag.
c                       the default value is off.
c
      lprops(19,elem) = .false.
C
C                       STORE DATA CONCERNING MATERIAL PROPERTIES.
C                       yes, this wastes space, but we maintain
c                       consistency in that all data for an element
c                       is here.
C                       Row            Quantity
c                       ---            --------
c                        7           young's modulus
c                                    -> kx_stiff for link2 element
c                        8           poisson's ratio. set = 0 for bar
c                                    -> ky_stiff for link2 element
c                        9           thermal expansion coeff.
c                                    -> kz_stiff for link2 element
c                        10          mass density
c                        13          thermal alphay for anisotropic
c                                    ** not used for link element
c                        14          kinematic hardening parameter
c                                      (beta)
c                                    ** not used for link element
c                        15          h-prime for linear hardening
c                                    ** not used for link element
c                        17          request for linear material
c                                    formulation
c                                    ** not used for link element
c                        20          viscopls m-power
c                                    ** not used for link element
c                        21          power-law hardening n-power or
c                                    ** not used for link element
c                                    set number for segmental
c                                    stress-strain curves
c                                    ** not used for link element
c                        22          viscoplastic ref strain rate
c                                    ( the user inputs D which the inverse
c                                      of eta, here we change D to eta
c                                      before passing the value to
c                                      props(22,elem) )
c                                    ** not used for link element
c                        23          uniaxial yield stress (sig-0)
c                                    ** not used for link element
c                        24          bit 1:  request for debug of material
c                                            computations
c                                    bit 2:  allow material model to
c                                            cutback adaptive solution
c                                            due to reversed plasticity
c                                    bit 3:  set of segmental stress-strain
c                                            curves specified for element's
c                                            material
c                                    ** not used for link element
c                        25          material model type
c                                    ** not used for link element
c                        26          initial porosity for gurson model
c                                    ** not used for link element
c                        27          q1 for gurson model
c                                    ** not used for link element
c                        28          q2 for gurson model
c                                    ** not used for link element
c                        29          q3 for gurson model
c                                    ** not used for link element
c                        30          bit 1:  logical nucleation flag for
c                                            gurson model
c                                    bit 2:  flag if element is
c                                            killable
c                                    ** not used for link element
c                        31          s_n parameter for gurson model
c                                    ** not used for link element
c                        32          eps_n parameter for gurson model
c                                    ** not used for link element
c                        33          f_n parameter for gurson model
c                                    ** not used for link element
c                        34          thermal alphaz for anisotropic
c                                    ** not used for link element
c                        35          thermal alphaxy for anisotropic
c                                    ** not used for link element
c                        36          thermal alphayz for anisotropic
c                                    ** not used for link element
c                        37          thermal alphaxz for anisotropic
c                                    ** not used for link element
c                        38          material number during input
c                                    ** not used for link element
c                        43          area for bar elements
c                                    ** not used for link element
c                                    ** not used for link element
c
c
c
      matnum          = elstor(2,elem)
      iprops(38,elem) = matnum
      props(7,elem)   = matprp(148,matnum)
      props(8,elem)   = matprp(149,matnum)
      props(9,elem)   = matprp(150,matnum)
c
c                       set unused real values to big negative to
c                       catch places where link2 element not handled
c
      props(13,elem) = rbig
      props(34,elem) = rbig
      props(35,elem) = rbig
      props(36,elem) = rbig
      props(37,elem) = rbig
c
      props(10,elem)  = matprp(7,matnum)
c
      props(14,elem)  = rbig
      props(15,elem)  = rbig
      lprops(17,elem) = .true.
      props(20,elem)  = rbig
      props(21,elem)  = rbig
      props(22,elem)  = rbig
      props(23,elem)  = rbig
      iprops(24,elem)  = 0
      iprops(25,elem) = 1 ! must be = 1 for link -> bilinear
      props(26,elem)  = rbig
      props(27,elem)  = rbig
      props(28,elem)  = rbig
      props(29,elem)  = rbig
      iprops(30,elem) = 0
      props(31,elem)  = rbig
      props(32,elem)  = rbig
      props(33,elem)  = rbig
      props(34,elem)  = rbig
      props(35,elem)  = rbig
      props(36,elem)  = rbig
      props(37,elem)  = rbig
      props(43,elem)  = rbig
c
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp1                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 7/10/01                    *
c     *                                                              *
c     *     this subroutine places the properties of the q3disop el- *
c     *     ement given into global storage.                         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine elprp1( elem, type )
      use global_data ! old common.main
      use main_data
      use segmental_curves
      implicit integer (a-z)
c
      character :: dums
      real e, et, h,  dumr
      double precision
     &    zero, one, dumd
      logical fgm
      data zero, one
     &    / 0.0, 1.0 /
c
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       see routine inamt for definition of rows of
c                       matprp table.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 20
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
c
c                       store data concerning the order of integration.
c                       the default value is 2x2x2. the other
c                       ordering permitted is the 14 pt rule.
c
      intord = elstor(3,elem)
c
      if( intord.eq.4HDEFA.or.intord.eq.4HO222 ) then
         iprops(5,elem) = 8
         iprops(6,elem) = 8
      else if( intord.eq.4HO14P ) then
         iprops(5,elem) = 9
         iprops(6,elem) = 14
      else if( intord.eq.4HO09P ) then
         iprops(5,elem) = 2
         iprops(6,elem) = 9
      else
         param = elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem) = 8
         iprops(6,elem) = 8
      end if
c
c                       store data concerning the location of output.
c                       the default value is gausspts. other options
c                       are element or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       store data concerning the output format for
c                       output. the default value is short.
c
      outfmt= elstor(6,elem)
      if( outfmt .eq. 4HDEFA .or. outfmt .eq. 4HSHRT ) then
         lprops(16,elem) = .false.
      else
         lprops(16,elem) = .true.
      end if
c
c                       store data concerning the geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c                       store data concerning the bbar flag.
c                       the default value is off. bbar not
c                       defined for 20-node brick
c
      lprops(19,elem) = .false.
C
C                       STORE DATA CONCERNING MATERIAL PROPERTIES.
C                       yes, this wastes space, but we maintain
c                       consistency in that all data for an element
c                       is here.
C                       Row            Quantity
c                       ---            --------
c                        7           young's modulus
c                        8           poisson's ratio
c                        9           thermal expansion coeff.
c                                    (alpha for isotropic, alphax for
c                                     anisotropic)
c                        10          mass density
c                        13          thermal alphay for anisotropic
c                        14          kinematic hardening parameter
c                                      (beta)
c                        15          h-prime for linear hardening
c                        17          request for linear material
c                                    formulation
c                        20          viscopls m-power
c                        21          power-law hardening n-power or
c                                    set number for segmental
c                                    stress-strain curves
c                        22          viscoplastic ref strain rate
c                                    ( the user inputs D which the inverse
c                                      of eta, here we change D to eta
c                                      before passing the value to
c                                      props(22,elem) )
c                        23          uniaxial yield stress (sig-0)
c                        24          bit 1:  request for debug of material
c                                            computations
c                                    bit 2:  allow material model to
c                                            cutback adaptive solution
c                                            due to reversed plasticity
c                                    bit 3:  set of segmental stress-strain
c                                            curves specified for element's
c                                            material
c                        25          material model type
c                        26          initial porosity for gurson model
c                        27          q1 for gurson model
c                        28          q2 for gurson model
c                        29          q3 for gurson model
c                        30          bit 1:  logical nucleation flag for
c                                            gurson model
c                                    bit 2:  flag if element is
c                                            killable
c                        31          s_n parameter for gurson model
c                        32          eps_n parameter for gurson model
c                        33          f_n parameter for gurson model
c                        34          thermal alphaz for anisotropic
c                        35          thermal alphaxy for anisotropic
c                        36          thermal alphayz for anisotropic
c                        37          thermal alphaxz for anisotropic
c                        38          material number during input
c
      matnum          = elstor(2,elem)
      iprops(38,elem) = matnum
      e               = matprp(1,matnum)
      et              = matprp(4,matnum)
      h               = (e*et)/(e-et)
      props(7,elem)   = e
      props(8,elem)   = matprp(2,matnum)
c
c                       set thermal expansion coeeficients. storage is
c                       always anisotopic. for istropic the three
c                       normal values are alpha, and shear values zero.
c
c                      alpha_{ij} is symmetric
c
      if ( lmtprp(25,matnum) ) then
        props(9,elem)  = matprp(26,matnum)  !  xx
        props(13,elem) = matprp(27,matnum)  !  yy
        props(34,elem) = matprp(28,matnum)  !  zz
        props(35,elem) = matprp(29,matnum)  !  xy
        props(36,elem) = matprp(30,matnum)  !  yz
        props(37,elem) = matprp(31,matnum)  !  xz
      else
        props(9,elem)  = matprp(6,matnum)
        props(13,elem) = matprp(6,matnum)
        props(34,elem) = matprp(6,matnum)
        props(35,elem) = 0.0
        props(36,elem) = 0.0
        props(37,elem) = 0.0
      end if
c
      props(10,elem)  = matprp(7,matnum)
      props(14,elem)  = matprp(3,matnum)
      props(15,elem)  = h
      lprops(17,elem) = lmtprp(8,matnum)
      props(20,elem)  = matprp(10,matnum)
      props(21,elem)  = matprp(11,matnum)
      props(22,elem)  = zero
      if ( matprp(12,matnum) .ne. zero )
     &  props(22,elem) = one / matprp(12,matnum)
      props(23,elem)  = matprp(5,matnum)
      iprops(24,elem)  = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      if ( lmtprp(24,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 4 )
      iprops(25,elem) = matprp(9,matnum)
      props(26,elem)  = matprp(14,matnum)
      props(27,elem)  = matprp(15,matnum)
      props(28,elem)  = matprp(16,matnum)
      props(29,elem)  = matprp(17,matnum)
      iprops(30,elem) = 0
      if ( lmtprp(18,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 1 )
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
      props(31,elem)  = matprp(19,matnum)
      props(32,elem)  = matprp(20,matnum)
      props(33,elem)  = matprp(21,matnum)
c
c                       if the material associated with the element has
c                       a set of stress-strain curves, store the set number
c                       for the curves and set the
c                       yield stress for the element material using the
c                       first curve in the set.
c
      call elem_set_segmental( elem, matnum )
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp2                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   last modified : 02/28/00                   *
c     *                                                              *
c     *     this subroutine places the properties of the l3disop el- *
c     *     ement given into global storage.                         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine elprp2( elem, type )
      use global_data ! old common.main
      use main_data
      use segmental_curves
      implicit integer (a-z)
c
      character :: dums
      real e, et, h,  dumr
      double precision
     &    zero, one, dumd
      logical local_debug
      data zero, one, local_debug
     &    / 0.0, 1.0, .true. /
c
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       see routine inmat for definition of rows of
c                       matprp table.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 8
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
c                       store data concerning the order of integration.
c                       the default value is 2x2x2.
c
      intord= elstor(3,elem)
c
      if(intord.eq.4HDEFA.OR.INTORD.EQ.4HO222) then
         iprops(5,elem)= 1
         iprops(6,elem)= 8
      else if(intord.eq.4HO06P) then
         iprops(5,elem)= 2
         iprops(6,elem)= 6
      else
         param= elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem)= 1
         iprops(6,elem)= 8
      end if
c
c                       store data concerning the location of output.
c                       the default value is gausspts. other options
c                       are element or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       store data concerning the output format for
c                       output. the default value is short.
c
      outfmt= elstor(6,elem)
c
      if(outfmt.eq.4HDEFA.OR.OUTFMT.EQ.4HSHRT) then
         lprops(16,elem)= .false.
      else
         lprops(16,elem)= .true.
      end if
c
c                       store data concerning the geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
c
      if(geonl.eq.4HDEFA.OR.GEONL.EQ.4HTRUE) then
         lprops(18,elem)= .true.
      else
         lprops(18,elem)= .false.
      end if
c
c                       store data concerning the bbar flag.
c                       the default value is on.
c
      bbar = elstor(8,elem)
c
      if(bbar.eq.4HTRUE) THEN
         LPROPS(19,ELEM)= .TRUE.
      ELSE
         LPROPS(19,ELEM)= .FALSE.
      END IF
C
C                       STORE DATA CONCERNING MATERIAL PROPERTIES.
C                       yes, this wastes space, but we maintain
c                       consistency in that all data for an element
c                       is here.
C                       Row            Quantity
c                       ---            --------
c                        7           young's modulus
c                        8           poisson's ratio
c                        9           thermal expansion coeff.
c                                    (alpha for isotropic, alphax for
c                                     anisotropic)
c                        10          mass density
c                        13          thermal alphay for anisotropic
c                        14          kinematic hardening parameter
c                                      (beta)
c                        15          h-prime for linear hardening
c                        17          request for linear material
c                                    formulation
c                        20          viscopls m-power
c                        21          power-law hardening n-power or
c                                    set number for segmental
c                                    stress-strain curves
c                        22          viscoplastic ref strain rate
c                                    ( the user inputs D which the inverse
c                                      of eta, here we change D to eta
c                                      before passing the value to
c                                      props(22,elem) )
c                        23          uniaxial yield stress (sig-0)
c                        24          bit 1:  request for debug of material
c                                            computations
c                                    bit 2:  allow material model to
c                                            cutback adaptive solution
c                                            due to reversed plasticity
c                                    bit 3:  set of segmental stress-strain
c                                            curves specified for element's
c                                            material
c                        25          material model type
c                        26          initial porosity for gurson model
c                        27          q1 for gurson model
c                        28          q2 for gurson model
c                        29          q3 for gurson model
c                        30          bit 1:  logical nucleation flag for
c                                            gurson model
c                                    bit 2:  flag if element is
c                                            killable
c                        31          s_n parameter for gurson model
c                        32          eps_n parameter for gurson model
c                        33          f_n parameter for gurson model
c                        34          thermal alphaz for anisotropic
c                        35          thermal alphaxy for anisotropic
c                        36          thermal alphayz for anisotropic
c                        37          thermal alphaxz for anisotropic
c                        38          material number during input
c
c
c
c
      matnum          = elstor(2,elem)
      iprops(38,elem) = matnum
      e               = matprp(1,matnum)
      et              = matprp(4,matnum)
      h               = (e*et)/(e-et)
      props(7,elem)   = e
      props(8,elem)   = matprp(2,matnum)
c
c                       set thermal expansion coeeficients. storage is
c                       always anisotopic. for istropic the three
c                       normal values are alpha, and shear values zero.
c
      if ( lmtprp(25,matnum) ) then
        props(9,elem)  = matprp(26,matnum)
        props(13,elem) = matprp(27,matnum)
        props(34,elem) = matprp(28,matnum)
        props(35,elem) = matprp(29,matnum)
        props(36,elem) = matprp(30,matnum)
        props(37,elem) = matprp(31,matnum)
      else
        props(9,elem)  = matprp(6,matnum)
        props(13,elem) = matprp(6,matnum)
        props(34,elem) = matprp(6,matnum)
        props(35,elem) = 0.0
        props(36,elem) = 0.0
        props(37,elem) = 0.0
      end if
c
      props(10,elem)  = matprp(7,matnum)
      props(14,elem)  = matprp(3,matnum)
      props(15,elem)  = h
      lprops(17,elem) = lmtprp(8,matnum)
      props(20,elem)  = matprp(10,matnum)
      props(21,elem)  = matprp(11,matnum)
      props(22,elem)  = zero
      if ( matprp(12,matnum) .ne. zero )
     &  props(22,elem) = one / matprp(12,matnum)
      props(23,elem)  = matprp(5,matnum)
      iprops(24,elem)  = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      if ( lmtprp(24,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 4 )
      iprops(25,elem) = matprp(9,matnum)
      props(26,elem)  = matprp(14,matnum)
      props(27,elem)  = matprp(15,matnum)
      props(28,elem)  = matprp(16,matnum)
      props(29,elem)  = matprp(17,matnum)
      iprops(30,elem) = 0
      if ( lmtprp(18,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 1 )
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
      props(31,elem)  = matprp(19,matnum)
      props(32,elem)  = matprp(20,matnum)
      props(33,elem)  = matprp(21,matnum)
c
c                       if the material associated with the element has
c                       a set of stress-strain curves, store the set number
c                       for the curves and set the
c                       yield stress for the element material using the
c                       first curve in the set.
c
      call elem_set_segmental( elem, matnum )
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp3                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/28/00                   *
c     *                                                              *
c     *     this subroutine places the properties of the ts12isop    *
c     *     element given into global storage.                       *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine elprp3( elem, type )
      use global_data ! old common.main
      use main_data
      use segmental_curves
      implicit integer (a-z)
c
      character :: dums
      real e, et, h,  dumr
      double precision
     &    zero, one, dumd
      data zero, one
     &    / 0.0, 1.0 /
c
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       see routine inamt for definition of rows of
c                       matprp table.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 12
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
c                       store data concerning the order of integration.
c                       the default value is 2x2x2.
c
      intord = elstor(3,elem)
c
      if( intord.eq.4HDEFA.or.intord.eq.4HO222 ) then
         iprops(5,elem) = 8
         iprops(6,elem) = 8
      else if( intord.eq.4HO14P ) then
         iprops(5,elem) = 9
         iprops(6,elem) = 14
      else if( intord.eq.4HO09P ) then
         iprops(5,elem) = 2
         iprops(6,elem) = 9
      else
         param = elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem) = 8
         iprops(6,elem) = 8
      end if
c
c                       store data concerning the location of output.
c                       the default value is gausspts. other options
c                       are element or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       store data concerning the output format for
c                       output. the default value is short.
c
      outfmt= elstor(6,elem)
      if( outfmt .eq. 4HDEFA .or. outfmt .eq. 4HSHRT ) then
         lprops(16,elem) = .false.
      else
         lprops(16,elem) = .true.
      end if
c
c                       store data concerning the geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c                       store data concerning the bbar flag.
c                       the default value is off. bbar not
c                       defined for 12-node brick
c
      lprops(19,elem) = .false.
C
C                       STORE DATA CONCERNING MATERIAL PROPERTIES.
C                       yes, this wastes space, but we maintain
c                       consistency in that all data for an element
c                       is here.
C                       Row            Quantity
c                       ---            --------
c                        7           young's modulus
c                        8           poisson's ratio
c                        9           thermal expansion coeff.
c                                    (alpha for isotropic, alphax for
c                                     anisotropic)
c                        10          mass density
c                        13          thermal alphay for anisotropic
c                        14          kinematic hardening parameter
c                                      (beta)
c                        15          h-prime for linear hardening
c                        17          request for linear material
c                                    formulation
c                        20          viscopls m-power
c                        21          power-law hardening n-power or
c                                    set number for segmental
c                                    stress-strain curves
c                        22          viscoplastic ref strain rate
c                                    ( the user inputs D which the inverse
c                                      of eta, here we change D to eta
c                                      before passing the value to
c                                      props(22,elem) )
c                        23          uniaxial yield stress (sig-0)
c                        24          bit 1:  request for debug of material
c                                            computations
c                                    bit 2:  allow material model to
c                                            cutback adaptive solution
c                                            due to reversed plasticity
c                                    bit 3:  set of segmental stress-strain
c                                            curves specified for element's
c                                            material
c                        25          material model type
c                        26          initial porosity for gurson model
c                        27          q1 for gurson model
c                        28          q2 for gurson model
c                        29          q3 for gurson model
c                        30          bit 1:  logical nucleation flag for
c                                            gurson model
c                                    bit 2:  flag if element is
c                                            killable
c                        31          s_n parameter for gurson model
c                        32          eps_n parameter for gurson model
c                        33          f_n parameter for gurson model
c                        34          thermal alphaz for anisotropic
c                        35          thermal alphaxy for anisotropic
c                        36          thermal alphayz for anisotropic
c                        37          thermal alphaxz for anisotropic
c                        38          material number during input
c
c
c
      matnum          = elstor(2,elem)
      iprops(38,elem) = matnum
      e               = matprp(1,matnum)
      et              = matprp(4,matnum)
      h               = (e*et)/(e-et)
      props(7,elem)   = e
      props(8,elem)   = matprp(2,matnum)
c
c                       set thermal expansion coeeficients. storage is
c                       always anisotopic. for istropic the three
c                       normal values are alpha, and shear values zero.
c
      if ( lmtprp(25,matnum) ) then
        props(9,elem)  = matprp(26,matnum)
        props(13,elem) = matprp(27,matnum)
        props(34,elem) = matprp(28,matnum)
        props(35,elem) = matprp(29,matnum)
        props(36,elem) = matprp(30,matnum)
        props(37,elem) = matprp(31,matnum)
      else
        props(9,elem)  = matprp(6,matnum)
        props(13,elem) = matprp(6,matnum)
        props(34,elem) = matprp(6,matnum)
        props(35,elem) = 0.0
        props(36,elem) = 0.0
        props(37,elem) = 0.0
      end if
c
      props(10,elem)  = matprp(7,matnum)
      props(14,elem)  = matprp(3,matnum)
      props(15,elem)  = h
      lprops(17,elem) = lmtprp(8,matnum)
      props(20,elem)  = matprp(10,matnum)
      props(21,elem)  = matprp(11,matnum)
      props(22,elem)  = zero
      if ( matprp(12,matnum) .ne. zero )
     &  props(22,elem) = one / matprp(12,matnum)
      props(23,elem)  = matprp(5,matnum)
      iprops(24,elem)  = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      if ( lmtprp(24,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 4 )
      iprops(25,elem) = matprp(9,matnum)
      props(26,elem)  = matprp(14,matnum)
      props(27,elem)  = matprp(15,matnum)
      props(28,elem)  = matprp(16,matnum)
      props(29,elem)  = matprp(17,matnum)
      iprops(30,elem) = 0
      if ( lmtprp(18,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 1 )
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
      props(31,elem)  = matprp(19,matnum)
      props(32,elem)  = matprp(20,matnum)
      props(33,elem)  = matprp(21,matnum)
c
c                       if the material associated with the element has
c                       a set of stress-strain curves, store the set number
c                       for the curves and set the
c                       yield stress for the element material using the
c                       first curve in the set.
c
      call elem_set_segmental( elem, matnum )
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp4                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/28/00                   *
c     *                                                              *
c     *     this subroutine places the properties of the ts15isop    *
c     *     element into global storage.                             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine elprp4( elem, type )
      use global_data ! old common.main
      use main_data
      use segmental_curves
      implicit integer (a-z)
c
      character :: dums
      real e, et, h,  dumr
      double precision
     &    zero, one, dumd
      data zero, one
     &    / 0.0, 1.0 /
c
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       see routine inamt for definition of rows of
c                       matprp table.
c
c                       store data directly tied only to element type.
c
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 15
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c                       store data concerning the order of integration.
c                       the default value is 2x2x2.
c
      intord = elstor(3,elem)
c
      if( intord.eq.4HDEFA.or.intord.eq.4HO222 ) then
         iprops(5,elem) = 8
         iprops(6,elem) = 8
      else if( intord.eq.4HO14P ) then
         iprops(5,elem) = 9
         iprops(6,elem) = 14
      else if( intord.eq.4HO09P ) then
         iprops(5,elem) = 2
         iprops(6,elem) = 9
      else
         param = elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem) = 8
         iprops(6,elem) = 8
      end if
c
c                       store data concerning the location of output.
c                       the default value is gausspts. other options
c                       are element or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       store data concerning the output format for
c                       output. the default value is short.
c
      outfmt= elstor(6,elem)
      if( outfmt .eq. 4HDEFA .or. outfmt .eq. 4HSHRT ) then
         lprops(16,elem) = .false.
      else
         lprops(16,elem) = .true.
      end if
c
c                       store data concerning the geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c                       store data concerning the bbar flag.
c                       the default value is off. bbar not
c                       defined for 15-node brick
c
      lprops(19,elem) = .false.
C
C                       STORE DATA CONCERNING MATERIAL PROPERTIES.
C                       yes, this wastes space, but we maintain
c                       consistency in that all data for an element
c                       is here.
C                       Row            Quantity
c                       ---            --------
c                        7           young's modulus
c                        8           poisson's ratio
c                        9           thermal expansion coeff.
c                                    (alpha for isotropic, alphax for
c                                     anisotropic)
c                        10          mass density
c                        13          thermal alphay for anisotropic
c                        14          kinematic hardening parameter
c                                      (beta)
c                        15          h-prime for linear hardening
c                        17          request for linear material
c                                    formulation
c                        20          viscopls m-power
c                        21          power-law hardening n-power or
c                                    set number for segmental
c                                    stress-strain curves
c                        22          viscoplastic ref strain rate
c                                    ( the user inputs D which the inverse
c                                      of eta, here we change D to eta
c                                      before passing the value to
c                                      props(22,elem) )
c                        23          uniaxial yield stress (sig-0)
c                        24          bit 1:  request for debug of material
c                                            computations
c                                    bit 2:  allow material model to
c                                            cutback adaptive solution
c                                            due to reversed plasticity
c                                    bit 3:  set of segmental stress-strain
c                                            curves specified for element's
c                                            material
c                        25          material model type
c                        26          initial porosity for gurson model
c                        27          q1 for gurson model
c                        28          q2 for gurson model
c                        29          q3 for gurson model
c                        30          bit 1:  logical nucleation flag for
c                                            gurson model
c                                    bit 2:  flag if element is
c                                            killable
c                        31          s_n parameter for gurson model
c                        32          eps_n parameter for gurson model
c                        33          f_n parameter for gurson model
c                        34          thermal alphaz for anisotropic
c                        35          thermal alphaxy for anisotropic
c                        36          thermal alphayz for anisotropic
c                        37          thermal alphaxz for anisotropic
c                        38          material number during input
c
c
      matnum          = elstor(2,elem)
      iprops(38,elem) = matnum
      e               = matprp(1,matnum)
      et              = matprp(4,matnum)
      h               = (e*et)/(e-et)
      props(7,elem)   = e
      props(8,elem)   = matprp(2,matnum)
c
c                       set thermal expansion coeeficients. storage is
c                       always anisotopic. for istropic the three
c                       normal values are alpha, and shear values zero.
c
      if ( lmtprp(25,matnum) ) then
        props(9,elem)  = matprp(26,matnum)
        props(13,elem) = matprp(27,matnum)
        props(34,elem) = matprp(28,matnum)
        props(35,elem) = matprp(29,matnum)
        props(36,elem) = matprp(30,matnum)
        props(37,elem) = matprp(31,matnum)
      else
        props(9,elem)  = matprp(6,matnum)
        props(13,elem) = matprp(6,matnum)
        props(34,elem) = matprp(6,matnum)
        props(35,elem) = 0.0
        props(36,elem) = 0.0
        props(37,elem) = 0.0
      end if
c
      props(10,elem)  = matprp(7,matnum)
      props(14,elem)  = matprp(3,matnum)
      props(15,elem)  = h
      lprops(17,elem) = lmtprp(8,matnum)
      props(20,elem)  = matprp(10,matnum)
      props(21,elem)  = matprp(11,matnum)
      props(22,elem)  = zero
      if ( matprp(12,matnum) .ne. zero )
     &  props(22,elem) = one / matprp(12,matnum)
      props(23,elem)  = matprp(5,matnum)
      iprops(24,elem)  = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      if ( lmtprp(24,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 4 )
      iprops(25,elem) = matprp(9,matnum)
      props(26,elem)  = matprp(14,matnum)
      props(27,elem)  = matprp(15,matnum)
      props(28,elem)  = matprp(16,matnum)
      props(29,elem)  = matprp(17,matnum)
      iprops(30,elem) = 0
      if ( lmtprp(18,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 1 )
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
      props(31,elem)  = matprp(19,matnum)
      props(32,elem)  = matprp(20,matnum)
      props(33,elem)  = matprp(21,matnum)
c
c                       if the material associated with the element has
c                       a set of stress-strain curves, store the set number
c                       for the curves and set the
c                       yield stress for the element material using the
c                       first curve in the set.
c
      call elem_set_segmental( elem, matnum )

      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp5                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/28/00                   *
c     *                                                              *
c     *     this subroutine places the properties of the ts9isop     *
c     *     element into global storage.                             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine elprp5( elem, type )
      use global_data ! old common.main
      use main_data
      use segmental_curves
      implicit integer (a-z)
c
      character :: dums
      real e, et, h,  dumr
      double precision
     &    zero, one, dumd
      data zero, one
     &    / 0.0, 1.0 /
c
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       see routine inamt for definition of rows of
c                       matprp table.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 9
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
c                       store data concerning the order of integration.
c                       the default value is 2x2x2.
c
      intord = elstor(3,elem)
c
      if( intord.eq.4HDEFA.or.intord.eq.4HO222 ) then
         iprops(5,elem) = 8
         iprops(6,elem) = 8
      else if( intord.eq.4HO14P ) then
         iprops(5,elem) = 9
         iprops(6,elem) = 14
      else if( intord.eq.4HO09P ) then
         iprops(5,elem) = 2
         iprops(6,elem) = 9
      else
         param = elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem) = 8
         iprops(6,elem) = 8
      end if
c
c                       store data concerning the location of output.
c                       the default value is gausspts. other options
c                       are element or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       store data concerning the output format for
c                       output. the default value is short.
c
      outfmt= elstor(6,elem)
      if( outfmt .eq. 4HDEFA .or. outfmt .eq. 4HSHRT ) then
         lprops(16,elem) = .false.
      else
         lprops(16,elem) = .true.
      end if
c
c                       store data concerning the geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c                       store data concerning the bbar flag.
c                       the default value is off. bbar not
c                       defined for 9-node brick
c
      lprops(19,elem) = .false.
C
C                       STORE DATA CONCERNING MATERIAL PROPERTIES.
C                       yes, this wastes space, but we maintain
c                       consistency in that all data for an element
c                       is here.
C                       Row            Quantity
c                       ---            --------
c                        7           young's modulus
c                        8           poisson's ratio
c                        9           thermal expansion coeff.
c                                    (alpha for isotropic, alphax for
c                                     anisotropic)
c                        10          mass density
c                        13          thermal alphay for anisotropic
c                        14          kinematic hardening parameter
c                                      (beta)
c                        15          h-prime for linear hardening
c                        17          request for linear material
c                                    formulation
c                        20          viscopls m-power
c                        21          power-law hardening n-power or
c                                    set number for segmental
c                                    stress-strain curves
c                        22          viscoplastic ref strain rate
c                                    ( the user inputs D which the inverse
c                                      of eta, here we change D to eta
c                                      before passing the value to
c                                      props(22,elem) )
c                        23          uniaxial yield stress (sig-0)
c                        24          bit 1:  request for debug of material
c                                            computations
c                                    bit 2:  allow material model to
c                                            cutback adaptive solution
c                                            due to reversed plasticity
c                                    bit 3:  set of segmental stress-strain
c                                            curves specified for element's
c                                            material
c                        25          material model type
c                        26          initial porosity for gurson model
c                        27          q1 for gurson model
c                        28          q2 for gurson model
c                        29          q3 for gurson model
c                        30          bit 1:  logical nucleation flag for
c                                            gurson model
c                                    bit 2:  flag if element is
c                                            killable
c                        31          s_n parameter for gurson model
c                        32          eps_n parameter for gurson model
c                        33          f_n parameter for gurson model
c                        34          thermal alphaz for anisotropic
c                        35          thermal alphaxy for anisotropic
c                        36          thermal alphayz for anisotropic
c                        37          thermal alphaxz for anisotropic
c                        38          material number during input
c
c
c
      matnum          = elstor(2,elem)
      iprops(38,elem) = matnum
      e               = matprp(1,matnum)
      et              = matprp(4,matnum)
      h               = (e*et)/(e-et)
      props(7,elem)   = e
      props(8,elem)   = matprp(2,matnum)
c
c                       set thermal expansion coeeficients. storage is
c                       always anisotopic. for istropic the three
c                       normal values are alpha, and shear values zero.
c
      if ( lmtprp(25,matnum) ) then
        props(9,elem)  = matprp(26,matnum)
        props(13,elem) = matprp(27,matnum)
        props(34,elem) = matprp(28,matnum)
        props(35,elem) = matprp(29,matnum)
        props(36,elem) = matprp(30,matnum)
        props(37,elem) = matprp(31,matnum)
      else
        props(9,elem)  = matprp(6,matnum)
        props(13,elem) = matprp(6,matnum)
        props(34,elem) = matprp(6,matnum)
        props(35,elem) = 0.0
        props(36,elem) = 0.0
        props(37,elem) = 0.0
      end if
c
      props(10,elem)  = matprp(7,matnum)
      props(14,elem)  = matprp(3,matnum)
      props(15,elem)  = h
      lprops(17,elem) = lmtprp(8,matnum)
      props(20,elem)  = matprp(10,matnum)
      props(21,elem)  = matprp(11,matnum)
      props(22,elem)  = zero
      if ( matprp(12,matnum) .ne. zero )
     &  props(22,elem) = one / matprp(12,matnum)
      props(23,elem)  = matprp(5,matnum)
      iprops(24,elem)  = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      if ( lmtprp(24,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 4 )
      iprops(25,elem) = matprp(9,matnum)
      props(26,elem)  = matprp(14,matnum)
      props(27,elem)  = matprp(15,matnum)
      props(28,elem)  = matprp(16,matnum)
      props(29,elem)  = matprp(17,matnum)
      iprops(30,elem) = 0
      if ( lmtprp(18,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 1 )
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
      props(31,elem)  = matprp(19,matnum)
      props(32,elem)  = matprp(20,matnum)
      props(33,elem)  = matprp(21,matnum)
c
c                       if the material associated with the element has
c                       a set of stress-strain curves, store the set number
c                       for the curves and set the
c                       yield stress for the element material using the
c                       first curve in the set.
c
      call elem_set_segmental( elem, matnum )
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp6                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/28/00                   *
c     *                                                              *
c     *     this subroutine places the properties of the tet10       *
c     *     element into global storage.                             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine elprp6( elem, type )
      use global_data ! old common.main
      use main_data
      use segmental_curves
      implicit integer (a-z)
c
      character :: dums
      real e, et, h,  dumr
      double precision
     &    zero, one, dumd
      data zero, one
     &    / 0.0, 1.0 /
c
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 10
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
c                       store data concerning the order of integration.
c                       the default value is a 4-point rule.
c
c                       see routine inmat for definition of rows of
c                       matprp table.
c
c                       for tetrahedral elements the possible integration
c                       orders have 1, 4, or 5 points
c
      intord = elstor(3,elem)
c
      if( intord.eq.4HDEFA.or.intord.eq.4HO04P ) then
         iprops(5,elem) = 4
         iprops(6,elem) = 4
      else if( intord.eq.4HO01P ) then
         iprops(5,elem) = 1
         iprops(6,elem) = 1
      else if( intord.eq.4HO05P ) then
         iprops(5,elem) = 5
         iprops(6,elem) = 5
      else
         param = elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem) = 4
         iprops(6,elem) = 4
      end if
c
c                       store data concerning the location of output.
c                       the default value is gausspts. other options
c                       are element or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       store data concerning the output format for
c                       output. the default value is short.
c
      outfmt= elstor(6,elem)
      if( outfmt .eq. 4HDEFA .or. outfmt .eq. 4HSHRT ) then
         lprops(16,elem) = .false.
      else
         lprops(16,elem) = .true.
      end if
c
c                       store data concerning the geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c                       store data concerning the bbar flag.
c                       the default value is off. bbar not
c                       defined for 10-node tet
c
      lprops(19,elem) = .false.
C
C                       STORE DATA CONCERNING MATERIAL PROPERTIES.
C                       yes, this wastes space, but we maintain
c                       consistency in that all data for an element
c                       is here.
C                       Row            Quantity
c                       ---            --------
c                        7           young's modulus
c                        8           poisson's ratio
c                        9           thermal expansion coeff.
c                                    (alpha for isotropic, alphax for
c                                     anisotropic)
c                        10          mass density
c                        13          thermal alphay for anisotropic
c                        14          kinematic hardening parameter
c                                      (beta)
c                        15          h-prime for linear hardening
c                        17          request for linear material
c                                    formulation
c                        20          viscopls m-power
c                        21          power-law hardening n-power or
c                                    set number for segmental
c                                    stress-strain curves
c                        22          viscoplastic ref strain rate
c                                    ( the user inputs D which the inverse
c                                      of eta, here we change D to eta
c                                      before passing the value to
c                                      props(22,elem) )
c                        23          uniaxial yield stress (sig-0)
c                        24          bit 1:  request for debug of material
c                                            computations
c                                    bit 2:  allow material model to
c                                            cutback adaptive solution
c                                            due to reversed plasticity
c                                    bit 3:  set of segmental stress-strain
c                                            curves specified for element's
c                                            material
c                        25          material model type
c                        26          initial porosity for gurson model
c                        27          q1 for gurson model
c                        28          q2 for gurson model
c                        29          q3 for gurson model
c                        30          bit 1:  logical nucleation flag for
c                                            gurson model
c                                    bit 2:  flag if element is
c                                            killable
c                        31          s_n parameter for gurson model
c                        32          eps_n parameter for gurson model
c                        33          f_n parameter for gurson model
c                        34          thermal alphaz for anisotropic
c                        35          thermal alphaxy for anisotropic
c                        36          thermal alphayz for anisotropic
c                        37          thermal alphaxz for anisotropic
c                        38          material number during input
c
c
      matnum          = elstor(2,elem)
      iprops(38,elem) = matnum
      e               = matprp(1,matnum)
      et              = matprp(4,matnum)
      h               = (e*et)/(e-et)
      props(7,elem)   = e
      props(8,elem)   = matprp(2,matnum)
c
c                       set thermal expansion coeeficients. storage is
c                       always anisotopic. for istropic the three
c                       normal values are alpha, and shear values zero.
c
      if ( lmtprp(25,matnum) ) then
        props(9,elem)  = matprp(26,matnum)
        props(13,elem) = matprp(27,matnum)
        props(34,elem) = matprp(28,matnum)
        props(35,elem) = matprp(29,matnum)
        props(36,elem) = matprp(30,matnum)
        props(37,elem) = matprp(31,matnum)
      else
        props(9,elem)  = matprp(6,matnum)
        props(13,elem) = matprp(6,matnum)
        props(34,elem) = matprp(6,matnum)
        props(35,elem) = 0.0
        props(36,elem) = 0.0
        props(37,elem) = 0.0
      end if
c
      props(10,elem)  = matprp(7,matnum)
      props(14,elem)  = matprp(3,matnum)
      props(15,elem)  = h
      lprops(17,elem) = lmtprp(8,matnum)
      props(20,elem)  = matprp(10,matnum)
      props(21,elem)  = matprp(11,matnum)
      props(22,elem)  = zero
      if ( matprp(12,matnum) .ne. zero )
     &  props(22,elem) = one / matprp(12,matnum)
      props(23,elem)  = matprp(5,matnum)
      iprops(24,elem)  = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      if ( lmtprp(24,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 4 )
      iprops(25,elem) = matprp(9,matnum)
      props(26,elem)  = matprp(14,matnum)
      props(27,elem)  = matprp(15,matnum)
      props(28,elem)  = matprp(16,matnum)
      props(29,elem)  = matprp(17,matnum)
      iprops(30,elem) = 0
      if ( lmtprp(18,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 1 )
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
      props(31,elem)  = matprp(19,matnum)
      props(32,elem)  = matprp(20,matnum)
      props(33,elem)  = matprp(21,matnum)
c
c                       if the material associated with the element has
c                       a set of stress-strain curves, store the set number
c                       for the curves and set the
c                       yield stress for the element material using the
c                       first curve in the set.
c
      call elem_set_segmental( elem, matnum )
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp7                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/28/00                   *
c     *                                                              *
c     *     this subroutine places the properties of the wedge15     *
c     *     element into global storage.                             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine elprp7( elem, type )
      use global_data ! old common.main
      use main_data
      use segmental_curves
      implicit integer (a-z)
c
      character :: dums
      real e, et, h,  dumr
      double precision
     &    zero, one, dumd
      data zero, one
     &    / 0.0, 1.0 /
c
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       see routine inmat for definition of rows of
c                       matprp table.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 15
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
c                       store data concerning the order of integration.
c                       the default value is 2x2x2.
c
      intord = elstor(3,elem)
c
      if( intord.eq.4HDEFA.or.intord.eq.4HO222 ) then
         iprops(5,elem) = 8
         iprops(6,elem) = 8
      else if( intord.eq.4HO14P ) then
         iprops(5,elem) = 9
         iprops(6,elem) = 14
      else if( intord.eq.4HO09P ) then
         iprops(5,elem) = 2
         iprops(6,elem) = 9
      else
         param = elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem) = 8
         iprops(6,elem) = 8
      end if
c
c                       store data concerning the location of output.
c                       the default value is gausspts. other options
c                       are element or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       store data concerning the output format for
c                       output. the default value is short.
c
      outfmt= elstor(6,elem)
      if( outfmt .eq. 4HDEFA .or. outfmt .eq. 4HSHRT ) then
         lprops(16,elem) = .false.
      else
         lprops(16,elem) = .true.
      end if
c
c                       store data concerning the geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c                       store data concerning the bbar flag.
c                       the default value is off. bbar not
c                       defined for 15-node wedge
c
      lprops(19,elem) = .false.
C
C                       STORE DATA CONCERNING MATERIAL PROPERTIES.
C                       yes, this wastes space, but we maintain
c                       consistency in that all data for an element
c                       is here.
C                       Row            Quantity
c                       ---            --------
c                        7           young's modulus
c                        8           poisson's ratio
c                        9           thermal expansion coeff.
c                                    (alpha for isotropic, alphax for
c                                     anisotropic)
c                        10          mass density
c                        13          thermal alphay for anisotropic
c                        14          kinematic hardening parameter
c                                      (beta)
c                        15          h-prime for linear hardening
c                        17          request for linear material
c                                    formulation
c                        20          viscopls m-power
c                        21          power-law hardening n-power or
c                                    set number for segmental
c                                    stress-strain curves
c                        22          viscoplastic ref strain rate
c                                    ( the user inputs D which the inverse
c                                      of eta, here we change D to eta
c                                      before passing the value to
c                                      props(22,elem) )
c                        23          uniaxial yield stress (sig-0)
c                        24          bit 1:  request for debug of material
c                                            computations
c                                    bit 2:  allow material model to
c                                            cutback adaptive solution
c                                            due to reversed plasticity
c                                    bit 3:  set of segmental stress-strain
c                                            curves specified for element's
c                                            material
c                        25          material model type
c                        26          initial porosity for gurson model
c                        27          q1 for gurson model
c                        28          q2 for gurson model
c                        29          q3 for gurson model
c                        30          bit 1:  logical nucleation flag for
c                                            gurson model
c                                    bit 2:  flag if element is
c                                            killable
c                        31          s_n parameter for gurson model
c                        32          eps_n parameter for gurson model
c                        33          f_n parameter for gurson model
c                        34          thermal alphaz for anisotropic
c                        35          thermal alphaxy for anisotropic
c                        36          thermal alphayz for anisotropic
c                        37          thermal alphaxz for anisotropic
c                        38          material number during input
c
c
      matnum          = elstor(2,elem)
      iprops(38,elem) = matnum
      e               = matprp(1,matnum)
      et              = matprp(4,matnum)
      h               = (e*et)/(e-et)
      props(7,elem)   = e
      props(8,elem)   = matprp(2,matnum)
c
c                       set thermal expansion coeeficients. storage is
c                       always anisotopic. for istropic the three
c                       normal values are alpha, and shear values zero.
c
      if ( lmtprp(25,matnum) ) then
        props(9,elem)  = matprp(26,matnum)
        props(13,elem) = matprp(27,matnum)
        props(34,elem) = matprp(28,matnum)
        props(35,elem) = matprp(29,matnum)
        props(36,elem) = matprp(30,matnum)
        props(37,elem) = matprp(31,matnum)
      else
        props(9,elem)  = matprp(6,matnum)
        props(13,elem) = matprp(6,matnum)
        props(34,elem) = matprp(6,matnum)
        props(35,elem) = 0.0
        props(36,elem) = 0.0
        props(37,elem) = 0.0
      end if
c
      props(10,elem)  = matprp(7,matnum)
      props(14,elem)  = matprp(3,matnum)
      props(15,elem)  = h
      lprops(17,elem) = lmtprp(8,matnum)
      props(20,elem)  = matprp(10,matnum)
      props(21,elem)  = matprp(11,matnum)
      props(22,elem)  = zero
      if ( matprp(12,matnum) .ne. zero )
     &  props(22,elem) = one / matprp(12,matnum)
      props(23,elem)  = matprp(5,matnum)
      iprops(24,elem)  = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      if ( lmtprp(24,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 4 )
      iprops(25,elem) = matprp(9,matnum)
      props(26,elem)  = matprp(14,matnum)
      props(27,elem)  = matprp(15,matnum)
      props(28,elem)  = matprp(16,matnum)
      props(29,elem)  = matprp(17,matnum)
      iprops(30,elem) = 0
      if ( lmtprp(18,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 1 )
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
      props(31,elem)  = matprp(19,matnum)
      props(32,elem)  = matprp(20,matnum)
      props(33,elem)  = matprp(21,matnum)
c
c                       if the material associated with the element has
c                       a set of stress-strain curves, store the set number
c                       for the curves and set the
c                       yield stress for the element material using the
c                       first curve in the set.
c
      call elem_set_segmental( elem, matnum )
c
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp8                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/28/00                   *
c     *                                                              *
c     *     this subroutine places the properties of the tri6        *
c     *     element into global storage.                             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine elprp8( elem, type )
      use global_data ! old common.main
      use main_data
      use segmental_curves
      implicit integer (a-z)
c
      character :: dums
      real e, et, h,  dumr
      double precision
     &    zero, one, dumd
      data zero, one
     &    / 0.0, 1.0 /
c
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       see routine inmat for definition of rows of
c                       matprp table.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 6
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
c                       store data concerning the order of integration.
c                       the default value is 2x2x2.
c
      intord = elstor(3,elem)
c
      if( intord.eq.4HDEFA.or.intord.eq.4HO222 ) then
         iprops(5,elem) = 8
         iprops(6,elem) = 8
      else if( intord.eq.4HO14P ) then
         iprops(5,elem) = 9
         iprops(6,elem) = 14
      else if( intord.eq.4HO09P ) then
         iprops(5,elem) = 2
         iprops(6,elem) = 9
      else
         param = elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem) = 8
         iprops(6,elem) = 8
      end if
c
c                       store data concerning the location of output.
c                       the default value is gausspts. other options
c                       are element or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       store data concerning the output format for
c                       output. the default value is short.
c
      outfmt= elstor(6,elem)
      if( outfmt .eq. 4HDEFA .or. outfmt .eq. 4HSHRT ) then
         lprops(16,elem) = .false.
      else
         lprops(16,elem) = .true.
      end if
c
c                       store data concerning the geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c                       store data concerning the bbar flag.
c                       the default value is off. bbar not
c                       defined for 6-node triangle
c
      lprops(19,elem) = .false.
C
C                       STORE DATA CONCERNING MATERIAL PROPERTIES.
C                       yes, this wastes space, but we maintain
c                       consistency in that all data for an element
c                       is here.
C                       Row            Quantity
c                       ---            --------
c                        7           young's modulus
c                        8           poisson's ratio
c                        9           thermal expansion coeff.
c                                    (alpha for isotropic, alphax for
c                                     anisotropic)
c                        10          mass density
c                        13          thermal alphay for anisotropic
c                        14          kinematic hardening parameter
c                                      (beta)
c                        15          h-prime for linear hardening
c                        17          request for linear material
c                                    formulation
c                        20          viscopls m-power
c                        21          power-law hardening n-power or
c                                    set number for segmental
c                                    stress-strain curves
c                        22          viscoplastic ref strain rate
c                                    ( the user inputs D which the inverse
c                                      of eta, here we change D to eta
c                                      before passing the value to
c                                      props(22,elem) )
c                        23          uniaxial yield stress (sig-0)
c                        24          bit 1:  request for debug of material
c                                            computations
c                                    bit 2:  allow material model to
c                                            cutback adaptive solution
c                                            due to reversed plasticity
c                                    bit 3:  set of segmental stress-strain
c                                            curves specified for element's
c                                            material
c                        25          material model type
c                        26          initial porosity for gurson model
c                        27          q1 for gurson model
c                        28          q2 for gurson model
c                        29          q3 for gurson model
c                        30          bit 1:  logical nucleation flag for
c                                            gurson model
c                                    bit 2:  flag if element is
c                                            killable
c                        31          s_n parameter for gurson model
c                        32          eps_n parameter for gurson model
c                        33          f_n parameter for gurson model
c                        34          thermal alphaz for anisotropic
c                        35          thermal alphaxy for anisotropic
c                        36          thermal alphayz for anisotropic
c                        37          thermal alphaxz for anisotropic
c                        38          material number during input
c
      matnum          = elstor(2,elem)
      iprops(38,elem) = matnum
      e               = matprp(1,matnum)
      et              = matprp(4,matnum)
      h               = (e*et)/(e-et)
      props(7,elem)   = e
      props(8,elem)   = matprp(2,matnum)
c
c                       set thermal expansion coeeficients. storage is
c                       always anisotopic. for istropic the three
c                       normal values are alpha, and shear values zero.
c
      if ( lmtprp(25,matnum) ) then
        props(9,elem)  = matprp(26,matnum)
        props(13,elem) = matprp(27,matnum)
        props(34,elem) = matprp(28,matnum)
        props(35,elem) = matprp(29,matnum)
        props(36,elem) = matprp(30,matnum)
        props(37,elem) = matprp(31,matnum)
      else
        props(9,elem)  = matprp(6,matnum)
        props(13,elem) = matprp(6,matnum)
        props(34,elem) = matprp(6,matnum)
        props(35,elem) = 0.0
        props(36,elem) = 0.0
        props(37,elem) = 0.0
      end if
c
      props(10,elem)  = matprp(7,matnum)
      props(14,elem)  = matprp(3,matnum)
      props(15,elem)  = h
      lprops(17,elem) = lmtprp(8,matnum)
      props(20,elem)  = matprp(10,matnum)
      props(21,elem)  = matprp(11,matnum)
      props(22,elem)  = zero
      if ( matprp(12,matnum) .ne. zero )
     &  props(22,elem) = one / matprp(12,matnum)
      props(23,elem)  = matprp(5,matnum)
      iprops(24,elem)  = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      if ( lmtprp(24,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 4 )
      iprops(25,elem) = matprp(9,matnum)
      props(26,elem)  = matprp(14,matnum)
      props(27,elem)  = matprp(15,matnum)
      props(28,elem)  = matprp(16,matnum)
      props(29,elem)  = matprp(17,matnum)
      iprops(30,elem) = 0
      if ( lmtprp(18,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 1 )
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
      props(31,elem)  = matprp(19,matnum)
      props(32,elem)  = matprp(20,matnum)
      props(33,elem)  = matprp(21,matnum)
c
c                       if the material associated with the element has
c                       a set of stress-strain curves, store the set number
c                       for the curves and set the
c                       yield stress for the element material using the
c                       first curve in the set.
c
      call elem_set_segmental( elem, matnum )
c
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp9                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/28/00                   *
c     *                                                              *
c     *     this subroutine places the properties of the quad8       *
c     *     element into global storage.                             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine elprp9( elem, type )
      use global_data ! old common.main
      use main_data
      use segmental_curves
      implicit integer (a-z)
c
      character :: dums
      real e, et, h,  dumr
      double precision
     &    zero, one, dumd
      data zero, one
     &    / 0.0, 1.0 /
c
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       see routine inmat for definition of rows of
c                       matprp table.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 8
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
c                       store data concerning the order of integration.
c                       the default value is 2x2x2.
c
      intord = elstor(3,elem)
c
      if( intord.eq.4HDEFA.or.intord.eq.4HO222 ) then
         iprops(5,elem) = 8
         iprops(6,elem) = 8
      else if( intord.eq.4HO14P ) then
         iprops(5,elem) = 9
         iprops(6,elem) = 14
      else if( intord.eq.4HO09P ) then
         iprops(5,elem) = 2
         iprops(6,elem) = 9
      else
         param = elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem) = 8
         iprops(6,elem) = 8
      end if
c
c                       store data concerning the location of output.
c                       the default value is gausspts. other options
c                       are element or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       store data concerning the output format for
c                       output. the default value is short.
c
      outfmt= elstor(6,elem)
      if( outfmt .eq. 4HDEFA .or. outfmt .eq. 4HSHRT ) then
         lprops(16,elem) = .false.
      else
         lprops(16,elem) = .true.
      end if
c
c                       store data concerning the geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c                       store data concerning the bbar flag.
c                       the default value is off. bbar not
c                       defined for 8-node quad
c
      lprops(19,elem) = .false.
C
C                       STORE DATA CONCERNING MATERIAL PROPERTIES.
C                       yes, this wastes space, but we maintain
c                       consistency in that all data for an element
c                       is here.
C                       Row            Quantity
c                       ---            --------
c                        7           young's modulus
c                        8           poisson's ratio
c                        9           thermal expansion coeff.
c                                    (alpha for isotropic, alphax for
c                                     anisotropic)
c                        10          mass density
c                        13          thermal alphay for anisotropic
c                        14          kinematic hardening parameter
c                                      (beta)
c                        15          h-prime for linear hardening
c                        17          request for linear material
c                                    formulation
c                        20          viscopls m-power
c                        21          power-law hardening n-power or
c                                    set number for segmental
c                                    stress-strain curves
c                        22          viscoplastic ref strain rate
c                                    ( the user inputs D which the inverse
c                                      of eta, here we change D to eta
c                                      before passing the value to
c                                      props(22,elem) )
c                        23          uniaxial yield stress (sig-0)
c                        24          bit 1:  request for debug of material
c                                            computations
c                                    bit 2:  allow material model to
c                                            cutback adaptive solution
c                                            due to reversed plasticity
c                                    bit 3:  set of segmental stress-strain
c                                            curves specified for element's
c                                            material
c                        25          material model type
c                        26          initial porosity for gurson model
c                        27          q1 for gurson model
c                        28          q2 for gurson model
c                        29          q3 for gurson model
c                        30          bit 1:  logical nucleation flag for
c                                            gurson model
c                                    bit 2:  flag if element is
c                                            killable
c                        31          s_n parameter for gurson model
c                        32          eps_n parameter for gurson model
c                        33          f_n parameter for gurson model
c                        34          thermal alphaz for anisotropic
c                        35          thermal alphaxy for anisotropic
c                        36          thermal alphayz for anisotropic
c                        37          thermal alphaxz for anisotropic
c                        38          material number during input
c
c
c
      matnum          = elstor(2,elem)
      iprops(38,elem) = matnum
      e               = matprp(1,matnum)
      et              = matprp(4,matnum)
      h               = (e*et)/(e-et)
      props(7,elem)   = e
      props(8,elem)   = matprp(2,matnum)
c
c                       set thermal expansion coeeficients. storage is
c                       always anisotopic. for istropic the three
c                       normal values are alpha, and shear values zero.
c
      if ( lmtprp(25,matnum) ) then
        props(9,elem)  = matprp(26,matnum)
        props(13,elem) = matprp(27,matnum)
        props(34,elem) = matprp(28,matnum)
        props(35,elem) = matprp(29,matnum)
        props(36,elem) = matprp(30,matnum)
        props(37,elem) = matprp(31,matnum)
      else
        props(9,elem)  = matprp(6,matnum)
        props(13,elem) = matprp(6,matnum)
        props(34,elem) = matprp(6,matnum)
        props(35,elem) = 0.0
        props(36,elem) = 0.0
        props(37,elem) = 0.0
      end if
c
      props(10,elem)  = matprp(7,matnum)
      props(14,elem)  = matprp(3,matnum)
      props(15,elem)  = h
      lprops(17,elem) = lmtprp(8,matnum)
      props(20,elem)  = matprp(10,matnum)
      props(21,elem)  = matprp(11,matnum)
      props(22,elem)  = zero
      if ( matprp(12,matnum) .ne. zero )
     &  props(22,elem) = one / matprp(12,matnum)
      props(23,elem)  = matprp(5,matnum)
      iprops(24,elem)  = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      if ( lmtprp(24,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 4 )
      iprops(25,elem) = matprp(9,matnum)
      props(26,elem)  = matprp(14,matnum)
      props(27,elem)  = matprp(15,matnum)
      props(28,elem)  = matprp(16,matnum)
      props(29,elem)  = matprp(17,matnum)
      iprops(30,elem) = 0
      if ( lmtprp(18,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 1 )
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
      props(31,elem)  = matprp(19,matnum)
      props(32,elem)  = matprp(20,matnum)
      props(33,elem)  = matprp(21,matnum)
c
c                       if the material associated with the element has
c                       a set of stress-strain curves, store the set number
c                       for the curves and set the
c                       yield stress for the element material using the
c                       first curve in the set.
c
      call elem_set_segmental( elem, matnum )
c
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp10                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/28/00                   *
c     *                                                              *
c     *     this subroutine places the properties of the axiquad8    *
c     *     element into global storage.                             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine elprp10( elem, type )
      use global_data ! old common.main
      use main_data
      use segmental_curves
      implicit integer (a-z)
c
      character :: dums
      real e, et, h,  dumr
      double precision
     &    zero, one, dumd
      data zero, one
     &    / 0.0, 1.0 /
c
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       see routine inmat for definition of rows of
c                       matprp table.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 8
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
c                       store data concerning the order of integration.
c                       the default value is 3x3 (9 int points).
c
c                       for quadrilateral elements the
c                       possible integration order = 1,2,3
c                       for 1x1, 2x2, 3x3 quadrillateral
c                       gauss integration.
c
      intord = elstor(3,elem)
c
      if( intord.eq.4HDEFA.or.intord.eq.4HO09P ) then
         iprops(5,elem) = 3
         iprops(6,elem) = 9
      else if( intord.eq.4HO01P ) then
         iprops(5,elem) = 1
         iprops(6,elem) = 1
      else if( intord.eq.4HO04P ) then
         iprops(5,elem) = 2
         iprops(6,elem) = 4
      else
         param = elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem) = 3
         iprops(6,elem) = 9
      end if
c
c                       store data concerning the location of output.
c                       the default value is gausspts. other options
c                       are element or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       store data concerning the output format for
c                       output. the default value is short.
c
      outfmt= elstor(6,elem)
      if( outfmt .eq. 4HDEFA .or. outfmt .eq. 4HSHRT ) then
         lprops(16,elem) = .false.
      else
         lprops(16,elem) = .true.
      end if
c
c                       store data concerning the geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c                       store data concerning the bbar flag.
c                       the default value is off. bbar not
c                       defined for 8-node axi quad
c
      lprops(19,elem) = .false.
C
C                       STORE DATA CONCERNING MATERIAL PROPERTIES.
C                       yes, this wastes space, but we maintain
c                       consistency in that all data for an element
c                       is here.
C                       Row            Quantity
c                       ---            --------
c                        7           young's modulus
c                        8           poisson's ratio
c                        9           thermal expansion coeff.
c                                    (alpha for isotropic, alphax for
c                                     anisotropic)
c                        10          mass density
c                        13          thermal alphay for anisotropic
c                        14          kinematic hardening parameter
c                                      (beta)
c                        15          h-prime for linear hardening
c                        17          request for linear material
c                                    formulation
c                        20          viscopls m-power
c                        21          power-law hardening n-power or
c                                    set number for segmental
c                                    stress-strain curves
c                        22          viscoplastic ref strain rate
c                                    ( the user inputs D which the inverse
c                                      of eta, here we change D to eta
c                                      before passing the value to
c                                      props(22,elem) )
c                        23          uniaxial yield stress (sig-0)
c                        24          bit 1:  request for debug of material
c                                            computations
c                                    bit 2:  allow material model to
c                                            cutback adaptive solution
c                                            due to reversed plasticity
c                                    bit 3:  set of segmental stress-strain
c                                            curves specified for element's
c                                            material
c                        25          material model type
c                        26          initial porosity for gurson model
c                        27          q1 for gurson model
c                        28          q2 for gurson model
c                        29          q3 for gurson model
c                        30          bit 1:  logical nucleation flag for
c                                            gurson model
c                                    bit 2:  flag if element is
c                                            killable
c                        31          s_n parameter for gurson model
c                        32          eps_n parameter for gurson model
c                        33          f_n parameter for gurson model
c                        34          thermal alphaz for anisotropic
c                        35          thermal alphaxy for anisotropic
c                        36          thermal alphayz for anisotropic
c                        37          thermal alphaxz for anisotropic
c                        38          material number during input
c
c
      matnum          = elstor(2,elem)
      iprops(38,elem) = matnum
      e               = matprp(1,matnum)
      et              = matprp(4,matnum)
      h               = (e*et)/(e-et)
      props(7,elem)   = e
      props(8,elem)   = matprp(2,matnum)
c
c                       set thermal expansion coeeficients. storage is
c                       always anisotopic. for istropic the three
c                       normal values are alpha, and shear values zero.
c
      if ( lmtprp(25,matnum) ) then
        props(9,elem)  = matprp(26,matnum)
        props(13,elem) = matprp(27,matnum)
        props(34,elem) = matprp(28,matnum)
        props(35,elem) = matprp(29,matnum)
        props(36,elem) = matprp(30,matnum)
        props(37,elem) = matprp(31,matnum)
      else
        props(9,elem)  = matprp(6,matnum)
        props(13,elem) = matprp(6,matnum)
        props(34,elem) = matprp(6,matnum)
        props(35,elem) = 0.0
        props(36,elem) = 0.0
        props(37,elem) = 0.0
      end if
c
      props(10,elem)  = matprp(7,matnum)
      props(14,elem)  = matprp(3,matnum)
      props(15,elem)  = h
      lprops(17,elem) = lmtprp(8,matnum)
      props(20,elem)  = matprp(10,matnum)
      props(21,elem)  = matprp(11,matnum)
      props(22,elem)  = zero
      if ( matprp(12,matnum) .ne. zero )
     &  props(22,elem) = one / matprp(12,matnum)
      props(23,elem)  = matprp(5,matnum)
      iprops(24,elem)  = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      if ( lmtprp(24,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 4 )
      iprops(25,elem) = matprp(9,matnum)
      props(26,elem)  = matprp(14,matnum)
      props(27,elem)  = matprp(15,matnum)
      props(28,elem)  = matprp(16,matnum)
      props(29,elem)  = matprp(17,matnum)
      iprops(30,elem) = 0
      if ( lmtprp(18,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 1 )
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
      props(31,elem)  = matprp(19,matnum)
      props(32,elem)  = matprp(20,matnum)
      props(33,elem)  = matprp(21,matnum)
c
c                       if the material associated with the element has
c                       a set of stress-strain curves, store the set number
c                       for the curves and set the
c                       yield stress for the element material using the
c                       first curve in the set.
c
      call elem_set_segmental( elem, matnum )
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp11                      *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/28/00                   *
c     *                                                              *
c     *     this subroutine places the properties of the axitri6     *
c     *     element into global storage.                             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine elprp11( elem, type )
      use global_data ! old common.main
      use main_data
      use segmental_curves
      implicit integer (a-z)
c
      character :: dums
      real e, et, h,  dumr
      double precision
     &    zero, one, dumd
      logical local_debug
      data zero, one, local_debug
     &    / 0.0, 1.0, .false. /
c
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       see routine inmat for definition of rows of
c                       matprp table.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 6
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
c                       store data concerning the order of integration.
c                       the default value is the 3 point rule.
c
c                       for triangle elements the
c                       possible integration order = 1,3,4,6,7
c
      intord = elstor(3,elem)
c
      if( intord.eq.4HDEFA.or.intord.eq.4HO03P ) then
         iprops(5,elem) = 3
         iprops(6,elem) = 3
      else if( intord.eq.4HO01P ) then
         iprops(5,elem) = 1
         iprops(6,elem) = 1
      else if( intord.eq.4HO04P ) then
         iprops(5,elem) = 4
         iprops(6,elem) = 4
      else if( intord.eq.4HO06P ) then
         iprops(5,elem) = 6
         iprops(6,elem) = 6
       else if( intord.eq.4HO07P ) then
         iprops(5,elem) = 7
         iprops(6,elem) = 7
      else
         param = elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem) = 3
         iprops(6,elem) = 3
      end if

      if( local_debug ) then
        write(*,500)elem,type,intord,iprops(5,elem),iprops(6,elem)
      end if
c
c                       store data concerning the location of output.
c                       the default value is gausspts. other options
c                       are element or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       store data concerning the output format for
c                       output. the default value is short.
c
      outfmt= elstor(6,elem)
      if( outfmt .eq. 4HDEFA .or. outfmt .eq. 4HSHRT ) then
         lprops(16,elem) = .false.
      else
         lprops(16,elem) = .true.
      end if
c
c                       store data concerning the geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c                       store data concerning the bbar flag.
c                       the default value is off. bbar not
c                       defined for 6-node axi triangle
c
      lprops(19,elem) = .false.
C
C                       STORE DATA CONCERNING MATERIAL PROPERTIES.
C                       yes, this wastes space, but we maintain
c                       consistency in that all data for an element
c                       is here.
C                       Row            Quantity
c                       ---            --------
c                        7           young's modulus
c                        8           poisson's ratio
c                        9           thermal expansion coeff.
c                                    (alpha for isotropic, alphax for
c                                     anisotropic)
c                        10          mass density
c                        13          thermal alphay for anisotropic
c                        14          kinematic hardening parameter
c                                      (beta)
c                        15          h-prime for linear hardening
c                        17          request for linear material
c                                    formulation
c                        20          viscopls m-power
c                        21          power-law hardening n-power or
c                                    set number for segmental
c                                    stress-strain curves
c                        22          viscoplastic ref strain rate
c                                    ( the user inputs D which the inverse
c                                      of eta, here we change D to eta
c                                      before passing the value to
c                                      props(22,elem) )
c                        23          uniaxial yield stress (sig-0)
c                        24          bit 1:  request for debug of material
c                                            computations
c                                    bit 2:  allow material model to
c                                            cutback adaptive solution
c                                            due to reversed plasticity
c                                    bit 3:  set of segmental stress-strain
c                                            curves specified for element's
c                                            material
c                        25          material model type
c                        26          initial porosity for gurson model
c                        27          q1 for gurson model
c                        28          q2 for gurson model
c                        29          q3 for gurson model
c                        30          bit 1:  logical nucleation flag for
c                                            gurson model
c                                    bit 2:  flag if element is
c                                            killable
c                        31          s_n parameter for gurson model
c                        32          eps_n parameter for gurson model
c                        33          f_n parameter for gurson model
c                        34          thermal alphaz for anisotropic
c                        35          thermal alphaxy for anisotropic
c                        36          thermal alphayz for anisotropic
c                        37          thermal alphaxz for anisotropic
c                        38          material number during input
c
      matnum          = elstor(2,elem)
      iprops(38,elem) = matnum
      e               = matprp(1,matnum)
      et              = matprp(4,matnum)
      h               = (e*et)/(e-et)
      props(7,elem)   = e
      props(8,elem)   = matprp(2,matnum)
c
c                       set thermal expansion coeeficients. storage is
c                       always anisotopic. for istropic the three
c                       normal values are alpha, and shear values zero.
c
      if ( lmtprp(25,matnum) ) then
        props(9,elem)  = matprp(26,matnum)
        props(13,elem) = matprp(27,matnum)
        props(34,elem) = matprp(28,matnum)
        props(35,elem) = matprp(29,matnum)
        props(36,elem) = matprp(30,matnum)
        props(37,elem) = matprp(31,matnum)
      else
        props(9,elem)  = matprp(6,matnum)
        props(13,elem) = matprp(6,matnum)
        props(34,elem) = matprp(6,matnum)
        props(35,elem) = 0.0
        props(36,elem) = 0.0
        props(37,elem) = 0.0
      end if
c
      props(10,elem)  = matprp(7,matnum)
      props(14,elem)  = matprp(3,matnum)
      props(15,elem)  = h
      lprops(17,elem) = lmtprp(8,matnum)
      props(20,elem)  = matprp(10,matnum)
      props(21,elem)  = matprp(11,matnum)
      props(22,elem)  = zero
      if ( matprp(12,matnum) .ne. zero )
     &  props(22,elem) = one / matprp(12,matnum)
      props(23,elem)  = matprp(5,matnum)
      iprops(24,elem)  = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      if ( lmtprp(24,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 4 )
      iprops(25,elem) = matprp(9,matnum)
      props(26,elem)  = matprp(14,matnum)
      props(27,elem)  = matprp(15,matnum)
      props(28,elem)  = matprp(16,matnum)
      props(29,elem)  = matprp(17,matnum)
      iprops(30,elem) = 0
      if ( lmtprp(18,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 1 )
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
      props(31,elem)  = matprp(19,matnum)
      props(32,elem)  = matprp(20,matnum)
      props(33,elem)  = matprp(21,matnum)
c
c                       if the material associated with the element has
c                       a set of stress-strain curves, store the set number
c                       for the curves and set the
c                       yield stress for the element material using the
c                       first curve in the set.
c
      call elem_set_segmental( elem, matnum )
c
c
      return
c
500     format(1x,//,'    In elprp11 (axitri6 element):',/,
     &               '          element number =',i7,/,
     &               '            element type =',i7,/,
     &               '    intord=elstor(3,elem)=',i14,/,
     &               '           iprops(5,elem)=',i7,/,
     &               '           iprops(6,elem)=',i7,/)
c
       end
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp12                      *
c     *                                                              *
c     *                       written by : aroy                      *
c     *                                                              *
c     *                   last modified : 04/16/2013 rhd             *
c     *                                 : add cavit option           *
c     *                                                              *
c     *     this subroutine places the properties of the inter_8     *
c     *     element given into global storage.                       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine elprp12( elem, type )
      use global_data ! old common.main
      use main_data
      implicit integer (a-z)
c
      character :: dums
      real  dumr
      double precision
     &    zero, one, dumd
      logical exponential_type, ppr_type, cavit_type
      data zero, one
     &    / 0.0, 1.0 /
c
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 8
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
c                       order of integration.
c                       the default value is 2x2 nodal lumping. the other
c                       rules permitted are 2x2 and  1 point Gauss
c                       Legendre rule.
c
      intord = elstor(3,elem)
c
c
      if( intord.eq.4HDEFA .or. intord.eq.4HO22G ) then
         iprops(5,elem) = 5
         iprops(6,elem) = 4
      else if( intord.eq.4HO22N ) then
         iprops(5,elem) = 4
         iprops(6,elem) = 4
      else if( intord.eq.4HO01P ) then
         iprops(5,elem) = 6
         iprops(6,elem) = 1
      else
         param = elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem) = 5
         iprops(6,elem) = 4
      end if
c
c                       location of output.
c                       the default value is gausspts. other options
c                       are element or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       output format. the default value is short.
c                       -- cohesive element can output only in short fmt.
c
      outfmt = elstor(6,elem)
      if( outfmt .eq. 4HDEFA .or. outfmt .eq. 4HSHRT ) then
         lprops(16,elem) = .false.
      else
         lprops(16,elem) = .false.
      end if
c
c                       geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c                       interface-cohesive elements do not currently
c                       use a bbar type formulation (as in solids)
c
      lprops(19,elem) = .false.
c
c                       reference surface
c                       for geometrical update for the cohesive element
c
      surf = elstor(10, elem)
      if( surf.eq.4HDEFA .or. surf.eq.4HO222 ) then
         iprops(26,elem) = 2
      else if( surf.eq.4HO111 ) then
         iprops(26,elem) = 1
      else if( surf.eq.4HO333 ) then
         iprops(26,elem) = 3
      else
         iprops(26,elem) = 2
      end if
c
c             Move material properties into (rprops,iprops,lprops) table.
c             yes, this wastes space, but we maintain
c             consistency in that all data for an element
c             is here.
c
c             props row      Quantity
c             ---------      --------
c              7    interface stiffness in the longitudinal direction
c              8    interface stiffness in the transverse direction
c              9    interface stiffness in the normal direction
c             10    a flag for identifying the interface element
c             13    peak normal stress of the interface
c             14    peak shear stress of the interface
c             15    shape parameter for bilinear and ramp
c             28    second ( additional) shape parameter for ramp]
c             17    request for linear material formulation
c             20    separation distance @ peak stress in sliding
c             21    separation distance @ peak stress in opening
c             22    equivalent separation distance @ peak stress
c             23    For exponential form: a ratio [beta_cohesive]to
c                       determine the equivalent separation under
c                       mixed mode loading ( =0 => mode I )
c                   For PPR: only_normal_mode = 2.0
c                            only_shear_mode = -2.0
c                            real mixed mode =  0.0
c                          use tests of > 1.0 or < -1.0; not "exact"
c                          equalities since we have floats
c             24    bit 1:  request for debug of material
c                           computations
c                   bit 2:  allow material model to
c                           cutback adaptive solution
c                           due to reversed plasticity
c                   bit 3:  segmental stress-strain def
c             25    material model type
c             27    type of interface models
c                   1 - linear elastic; 2- bilinear
c                   3 - ramp; 4 - exponential_1; 5 - exponential_2;
c                   6 - PPR, 7 - cavit
c             28    second ( additional) shape parameter for ramp
c             29    compression (normal) stiffness multiplier
c             30    = 0 element is not marked killable;
c                     1 element is marked as killable
c             32    = 0 element is not yet killed
c                   = 1 element is currently killed
c             33    volume fraction of ductile matl in fgm cohesive
c             34    = 0 homogeneous cohesive
c                   = 1 fgm cohesive
c       ***   35  -- critical separation distance ductile (fgm)
c                 -- G_normal
c       ***   36  -- critical separation distance brittle (fgm)
c                 -- G_shear
c       ***   37  -- critical stress ductile (fgm)
c                 -- ratio_nomal
c             38  -- material number
c       ***   39  -- critical stress brittle (fgm)
c                 -- ratio_shear
c       ***   40  -- beta_ductile (fgm)
c                 -- shape_normal
c       ***   41  -- beta_brittle (fgm)
c                 -- shape_shear
c       (***)     space used by PPR option [cannot be used with FGM]
c
c
c                       data in the material properties table to be moved into p
c                       --------------------------------------------------------
c                       the model supports various options for the cohesive
c                       formulation
c                             1 - linear elastic;
c                             2 - bilinear
c                             3 - ramp;
c                             4 - exponential_1;
c                             5 - exponential_2;
c                             6 - PPR (Park, Paulino, Roesler)
c                             7 - creep cavitation
c
c                       At this time, only options 1, 4, 6, 7 are defined
c                       For options 1 and 4, the model also supports a FGM
c                       formulation through additional material properties.
c                       We still scan and store properties for the yet
c                       to be implemented
c
c      Changes to support major overhaul the the model
c      in summer 2010.
c
c      *Note* some properties are defined/translated from input
c             for options not yet implemented.
c
c      The current options are implemented:
c        1: linear, 3: exponential including fgm, 6: PPR, 7: cavitation
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
c                 makes normal + shear into single term via
c                 beta_coh(esive)  value
c          42  -- For exponential form: a ratio [beta_cohesive] to
c                       determine the equivalent separation under
c                       mixed mode loading ( =0 => mode I )
c                 For PPR: only_normal_mode = 2.0
c                          only_shear_mode = -2.0
c                          real mixed mode =  0.0
c                          use tests of > 1.0 or < -1.0; not "exact"
c                          equalities since we have floats
c          43  -- a flag for identifying the interface element (logical)
c          44  -- type of interface models
c                 1 - linear elastic; 2- bilinear
c                 3 - ramp; 4 - exponential_1; 5 - exponential_2;
c                 6 - PPR
c      (*) 46  -- ductile material volume fracture for fgm
c          47  -- = 0 homogeneous cohesive material (real)
c                 = 1 functionally graded cohesive material
c          48  -- separation distance at peak stress for ductile (fgm)
c          49  -- separation distance brittle (fgm)
c          50  -- maximum normal stress on curve ductile (fgm)
c          51  -- maximum shear stress on curve brittle (fgm)
c          52  -- beta_ductile (fgm) [beta cohesive]
c          53  -- beta_brittle (fgm) [beta cohesive]
c          54  -- compression stiffness multiplier for
c                 cohesive materials
c          90  -- PPR cohesive model
c                  90 - sig_max
c                  91 - tau_max
c                  92 - G_normal
c                  93 - G_shear
c                  94 - ratio_normal
c                  95 - ratio_shear
c                  96 - shape_normal
c                  97 - shape_shear
c          98 --  debug flag for material model computations
c
c          130-139 -- cavitation option for cohesive
c
      matnum          = elstor(2,elem)
      props(7,elem)   = matprp(32,matnum)
      props(8,elem)   = matprp(33,matnum)
      props(9,elem)   = matprp(34,matnum)
      props(13,elem)  = matprp(35,matnum)
      props(14,elem)  = matprp(36,matnum)
      props(15,elem)  = matprp(37,matnum)
      props(28,elem)  = matprp(38,matnum)
      props(20,elem)  = matprp(39,matnum)
      props(21,elem)  = matprp(40,matnum)
      props(22,elem)  = matprp(41,matnum)
c
c                 for symmetric problems ie, when reference surface is top
c                 or bottom, beta_coh (props(23)) should be zero for the
c                 "exponential" model. For the PPR option, it should
c                 be -2.0 or 2.0 (for only shear or only normal mode) but
c                 not zero.
c                 if user has input a non-consistent value, serve a warning
c                 and set a correct default value.
c
      props(23,elem)   = matprp(42,matnum)
      iprops(27,elem)  = matprp(44,matnum)
      lprops(10,elem)  = .true.
      exponential_type = iprops(27,elem) .eq. 4
      ppr_type         = iprops(27,elem) .eq. 6
      cavit_type       = iprops(27,elem) .eq. 7
c
c      if( exponential_type ) then
c        if ( surf.eq.4HO111 .and. matprp(42,matnum) .ne. zero ) then
c       call errmsg2(33,elem,'top',matprp(42,matnum),dumd)
c         props(23,elem) = zero
c        end if
c        if ( surf.eq.4HO333 .and. matprp(42,matnum) .ne. zero ) then
c       call errmsg2(33,elem,'bottom',matprp(42,matnum),dumd)
c         props(23,elem) = zero
c        end if
c      end if
c
c      if( ppr_type ) then
c        if ( surf.eq.4HO111 .and. matprp(42,matnum) .eq. zero ) then
c       call errmsg2(78,elem,'top',matprp(42,matnum),dumd)
c        end if
c        if ( surf.eq.4HO333 .and. matprp(42,matnum) .eq. zero ) then
c       call errmsg2(78,elem,'bottom',matprp(42,matnum),dumd)
c        end if
c      end if
c
c                  debug flag, allow cutback request from material,
c                  no segmental-type stress-strain curve.
c
      iprops(24,elem) = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      iprops(25,elem) = matprp(9,matnum)
c
c                      mark element killable if user has that
c                      property set for material.
c
      iprops(30,elem) = 0
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
c
c                      compression multiplier, element is not yet killed,
c                      fgm properties
c
      props(29,elem)  = matprp(54,matnum)
      iprops(32,elem) = 0
      props(33,elem)  = matprp(46,matnum)
      props(34,elem)  = matprp(47,matnum)
      props(35,elem)  = matprp(48,matnum)
      props(36,elem)  = matprp(49,matnum)
      props(37,elem)  = matprp(50,matnum)
      iprops(38,elem) = matnum
      props(39,elem)  = matprp(51,matnum)
      props(40,elem)  = matprp(52,matnum)
      props(41,elem)  = matprp(53,matnum)
c
c                     for PPR formulation, use spaces for FGM properties
c                     to store the PPR properties. Just keeps from
c                     making props table longer.
c
      if ( ppr_type ) then
          props(35,elem)  = matprp(92,matnum)
          props(36,elem)  = matprp(93,matnum)
          props(37,elem)  = matprp(94,matnum)
          props(39,elem)  = matprp(95,matnum)
          props(40,elem)  = matprp(96,matnum)
          props(41,elem)  = matprp(97,matnum)
      end if
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp13                      *
c     *                                                              *
c     *                       written by : rau                       *
c     *                                                              *
c     *                   last modified : 11/29/00                   *
c     *                                                              *
c     *     this subroutine places the properties of the tet4        *
c     *     element into global storage.                             *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine elprp13( elem, type )
      use global_data ! old common.main
      use main_data
      use segmental_curves
      implicit integer (a-z)
c
      character :: dums
c
      real e, et, h,  dumr
      double precision
     &    zero, one, dumd
      data zero, one
     &    / 0.0, 1.0 /
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
c
      iprops(2,elem)  = 4
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
c                       store data concerning the order of integration.
c                       the default value is a 1-point rule.
c
c                       for tetrahedral elements the possible integration
c                       orders have 1, 4, or 5 points
c
      intord = elstor(3,elem)
c
      if( intord.eq.4HDEFA.or.intord.eq.4HO01P ) then
         iprops(5,elem) = 1
         iprops(6,elem) = 1
      else if( intord.eq.4HO04P ) then
         iprops(5,elem) = 4
         iprops(6,elem) = 4
      else if( intord.eq.4HO05P ) then
         iprops(5,elem) = 5
         iprops(6,elem) = 5
      else
         param = elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem) = 1
         iprops(6,elem) = 1
      end if
c
c                       store data concerning the location of output.
c                       the default value is gausspts. other options
c                       are element or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       store data concerning the output format for
c                       output. the default value is short.
c
      outfmt= elstor(6,elem)
      if( outfmt .eq. 4HDEFA .or. outfmt .eq. 4HSHRT ) then
         lprops(16,elem) = .false.
      else
         lprops(16,elem) = .true.
      end if
c
c                       store data concerning the geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c                       store data concerning the bbar flag.
c                       the default value is off. bbar not
c                       defined for 10-node tet
c
      lprops(19,elem) = .false.
C
C                       STORE DATA CONCERNING MATERIAL PROPERTIES.
C                       yes, this wastes space, but we maintain
c                       consistency in that all data for an element
c                       is here.
C                       Row            Quantity
c                       ---            --------
c                        7           young's modulus
c                        8           poisson's ratio
c                        9           thermal expansion coeff.
c                                    (alpha for isotropic, alphax for
c                                     anisotropic)
c                        10          mass density
c                        13          thermal alphay for anisotropic
c                        14          kinematic hardening parameter
c                                      (beta)
c                        15          h-prime for linear hardening
c                        17          request for linear material
c                                    formulation
c                        20          viscopls m-power
c                        21          power-law hardening n-power or
c                                    set number for segmental
c                                    stress-strain curves
c                        22          viscoplastic ref strain rate
c                                    ( the user inputs D which the inverse
c                                      of eta, here we change D to eta
c                                      before passing the value to
c                                      props(22,elem) )
c                        23          uniaxial yield stress (sig-0)
c                        24          bit 1:  request for debug of material
c                                            computations
c                                    bit 2:  allow material model to
c                                            cutback adaptive solution
c                                            due to reversed plasticity
c                                    bit 3:  set of segmental stress-strain
c                                            curves specified for element's
c                                            material
c                        25          material model type
c                        26          initial porosity for gurson model
c                        27          q1 for gurson model
c                        28          q2 for gurson model
c                        29          q3 for gurson model
c                        30          bit 1:  logical nucleation flag for
c                                            gurson model
c                                    bit 2:  flag if element is
c                                            killable
c                        31          s_n parameter for gurson model
c                        32          eps_n parameter for gurson model
c                        33          f_n parameter for gurson model
c                        34          thermal alphaz for anisotropic
c                        35          thermal alphaxy for anisotropic
c                        36          thermal alphayz for anisotropic
c                        37          thermal alphaxz for anisotropic
c                        38          material number during input
c
c
      matnum          = elstor(2,elem)
      iprops(38,elem) = matnum
      e               = matprp(1,matnum)
      et              = matprp(4,matnum)
      h               = (e*et)/(e-et)
      props(7,elem)   = e
      props(8,elem)   = matprp(2,matnum)
c
c                       set thermal expansion coeeficients. storage is
c                       always anisotopic. for istropic the three
c                       normal values are alpha, and shear values zero.
c
      if ( lmtprp(25,matnum) ) then
        props(9,elem)  = matprp(26,matnum)
        props(13,elem) = matprp(27,matnum)
        props(34,elem) = matprp(28,matnum)
        props(35,elem) = matprp(29,matnum)
        props(36,elem) = matprp(30,matnum)
        props(37,elem) = matprp(31,matnum)
      else
        props(9,elem)  = matprp(6,matnum)
        props(13,elem) = matprp(6,matnum)
        props(34,elem) = matprp(6,matnum)
        props(35,elem) = 0.0
        props(36,elem) = 0.0
        props(37,elem) = 0.0
      end if
c
      props(10,elem)  = matprp(7,matnum)
      props(14,elem)  = matprp(3,matnum)
      props(15,elem)  = h
      lprops(17,elem) = lmtprp(8,matnum)
      props(20,elem)  = matprp(10,matnum)
      props(21,elem)  = matprp(11,matnum)
      props(22,elem)  = zero
      if ( matprp(12,matnum) .ne. zero )
     &  props(22,elem) = one / matprp(12,matnum)
      props(23,elem)  = matprp(5,matnum)
      iprops(24,elem)  = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      if ( lmtprp(24,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 4 )
      iprops(25,elem) = matprp(9,matnum)
      props(26,elem)  = matprp(14,matnum)
      props(27,elem)  = matprp(15,matnum)
      props(28,elem)  = matprp(16,matnum)
      props(29,elem)  = matprp(17,matnum)
      iprops(30,elem) = 0
      if ( lmtprp(18,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 1 )
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
      props(31,elem)  = matprp(19,matnum)
      props(32,elem)  = matprp(20,matnum)
      props(33,elem)  = matprp(21,matnum)
c
c                       if the material associated with the element has
c                       a set of stress-strain curves, store the set number
c                       for the curves and set the
c                       yield stress for the element material using the
c                       first curve in the set.
c
      call elem_set_segmental( elem, matnum )
c
      return
      end
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp14                      *
c     *                                                              *
c     *                       written by : sushovan                  *
c     *                                                              *
c     *                   last modified : 04/16/2013 rhd             *
c     *                                 : add cavit option           *
c     *                                                              *
c     *     this subroutine places the properties of the trint6      *
c     *     element given into global storage.                       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine elprp14( elem, type )
      use global_data ! old common.main
      use main_data
      implicit integer (a-z)
c
      character :: dums
      real dumr
      double precision
     &    dumd,zero
      logical exponential_type, ppr_type, cavit_type
c
      data zero /0.0/
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 6
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
c                       the order of integration.
c                       the default value is 3-midpoint rule (mid-point
c                       of elem sides). other rules permitted are a
c                       second 3-point rule and a 1 point rule.
c
      intord = elstor(3,elem)
c
c
      if( intord.eq.4HDEFA .or. intord.eq.4HO3MP ) then
         iprops(5,elem) = 3
         iprops(6,elem) = 3
      else if( intord.eq.4HO03P ) then
         iprops(5,elem) = 2
         iprops(6,elem) = 3
      else if( intord.eq.4HO01P ) then
         iprops(5,elem) = 1
         iprops(6,elem) = 1
      else
         param = elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem) = 3
         iprops(6,elem) = 3
      end if
c
c                       location of output.
c                       the default value is gausspts. other options
c                       are nodes or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       output format. the default value is short.
c                       -- cohesive element can output only in short fmt.
c
      outfmt = elstor(6,elem)
      if( outfmt .eq. 4HDEFA .or. outfmt .eq. 4HSHRT ) then
         lprops(16,elem) = .false.
      else
         lprops(16,elem) = .false.
      end if
c
c                       geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c
c                       interface-cohesive elements do not currently
c                       use a bbar type formulation (as in solids)
c
      lprops(19,elem) = .false.
c
c
c                       reference surface
c                       for geometrical update for the cohesive element
c
      surf = elstor(10, elem)
      if( surf.eq.4HDEFA .or. surf.eq.4HO222 ) then
         iprops(26,elem) = 2
      else if( surf.eq.4HO111 ) then
         iprops(26,elem) = 1
      else if( surf.eq.4HO333 ) then
         iprops(26,elem) = 3
      else
         iprops(26,elem) = 2
      end if
c
c                       See all comments in elprp12. They apply here.
c
      matnum          = elstor(2,elem)
      props(7,elem)   = matprp(32,matnum)
      props(8,elem)   = matprp(33,matnum)
      props(9,elem)   = matprp(34,matnum)
      props(13,elem)  = matprp(35,matnum)
      props(14,elem)  = matprp(36,matnum)
      props(15,elem)  = matprp(37,matnum)
      props(28,elem)  = matprp(38,matnum)
      props(20,elem)  = matprp(39,matnum)
      props(21,elem)  = matprp(40,matnum)
      props(22,elem)  = matprp(41,matnum)
c
c                 for symmetric problems ie, when reference surface is top
c                 or bottom, beta_coh (props(23)) should be zero for the
c                 "exponential" model. For the PPR option, it should
c                 be -2.0 or 2.0 (for only shear or only normal mode) but
c                 not zero.
c                 if user has input a non-consistent value, serve a warning
c                 and set a correct default value.
c
      props(23,elem)   = matprp(42,matnum)
      iprops(27,elem)  = matprp(44,matnum)
      lprops(10,elem)  = .true.
      exponential_type = iprops(27,elem) .eq. 4
      ppr_type         = iprops(27,elem) .eq. 6
      cavit_type       = iprops(27,elem) .eq. 7
c
      if( exponential_type ) then
        if ( surf.eq.4HO111 .and. matprp(42,matnum) .ne. zero ) then
       call errmsg2(33,elem,'top',matprp(42,matnum),dumd)
         props(23,elem) = zero
        end if
        if ( surf.eq.4HO333 .and. matprp(42,matnum) .ne. zero ) then
       call errmsg2(33,elem,'bottom',matprp(42,matnum),dumd)
         props(23,elem) = zero
        end if
      end if
c
      if( ppr_type ) then
        if ( surf.eq.4HO111 .and. matprp(42,matnum) .eq. zero ) then
       call errmsg2(78,elem,'top',matprp(42,matnum),dumd)
        end if
        if ( surf.eq.4HO333 .and. matprp(42,matnum) .eq. zero ) then
       call errmsg2(78,elem,'bottom',matprp(42,matnum),dumd)
        end if
      end if
c
c                  debug flag, allow cutback request from material,
c                  no segmental-type stress-strain curve.
c
      iprops(24,elem) = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      iprops(25,elem) = matprp(9,matnum)
c
c                      mark element killable if user has that
c                      property set for material.
c
      iprops(30,elem) = 0
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
      iprops(25,elem) = matprp(9,matnum)
c
c                      mark element killable if user has that
c                      property set for material.
c
      iprops(30,elem) = 0
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
c
c                      compression multiplier, element is not yet killed,
c                      fgm properties
c
      props(29,elem)  = matprp(54,matnum)
      iprops(32,elem) = 0
      props(33,elem)  = matprp(46,matnum)
      props(34,elem)  = matprp(47,matnum)
      props(35,elem)  = matprp(48,matnum)
      props(36,elem)  = matprp(49,matnum)
      props(37,elem)  = matprp(50,matnum)
      iprops(38,elem) = matnum
      props(39,elem)  = matprp(51,matnum)
      props(40,elem)  = matprp(52,matnum)
      props(41,elem)  = matprp(53,matnum)
c
c                     for PPR formulation, use spaces for FGM properties
c                     to store the PPR properties. Just keeps from
c                     making props table longer.
c
      if ( ppr_type ) then
          props(35,elem)  = matprp(92,matnum)
          props(36,elem)  = matprp(93,matnum)
          props(37,elem)  = matprp(94,matnum)
          props(39,elem)  = matprp(95,matnum)
          props(40,elem)  = matprp(96,matnum)
          props(41,elem)  = matprp(97,matnum)
      end if
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine elprp15                      *
c     *                                                              *
c     *                       written by : sushovan                  *
c     *                                                              *
c     *                   last modified : 08/11/2010 RDH             *
c     *                                 : add PPR option             *
c     *                                                              *
c     *     this subroutine places the properties of the trint12     *
c     *     element given into global storage.                       *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine elprp15( elem, type )
      use global_data ! old common.main
      use main_data
      implicit integer (a-z)
c
      character :: dums
      real dumr
      double precision
     &    dumd,zero
      logical exponential_type, ppr_type, cavit_type
c
      data zero /0.0/
c
c                       note: data stored in the element properties
c                       table is all single precision.
c
c                       store data directly tied only to element type.
c
      iprops(1,elem)  = type
      iprops(2,elem)  = 12
      iprops(3,elem)  = 24*two16+6
      iprops(4,elem)  = 3
      iprops(11,elem) = 6
c
c                       Store the flag/number indicating if this is
c                       an interface damage material
      iprops(42,elem) = elstor(11,elem)
c
c                       order of integration.
c                       the default value is 7 point rule. Other
c                       rules permitted are two 3 point rules and
c                       a 1 point rule.
c
      intord = elstor(3,elem)
c
      if( intord.eq.4HDEFA .or. intord.eq.4HO07P ) then
         iprops(5,elem) = 7
         iprops(6,elem) = 7
      else if( intord.eq.4HO06P ) then
         iprops(5,elem) = 6
         iprops(6,elem) = 6
      else if( intord.eq.4HO04P ) then
         iprops(5,elem) = 4
         iprops(6,elem) = 4
      else if( intord.eq.4HO3MP ) then
         iprops(5,elem) = 3
         iprops(6,elem) = 3
      else if( intord.eq.4HO03P ) then
         iprops(5,elem) = 2
         iprops(6,elem) = 3
      else if( intord.eq.4HO01P ) then
         iprops(5,elem) = 1
         iprops(6,elem) = 1
      else
         param = elem
         call errmsg(33,elem,dums,dumr,dumd)
         iprops(5,elem) = 7
         iprops(6,elem) = 7
      end if
c
c                       location of output.
c                       the default value is gausspts. other options
c                       are nodes or element center.
c
      outloc = elstor(4,elem)
      iloc   = 1
      if ( outloc .eq. 4hNODE ) then
         iloc = 2
      else if ( outloc .eq. 4hCENT ) then
         iloc = 3
      end if
      iprops(12,elem) = iloc
c
c                       output format. the default value is short.
c                       -- cohesive element can output only in short fmt.
c
      outfmt = elstor(6,elem)
      if( outfmt .eq. 4HDEFA .or. outfmt .eq. 4HSHRT ) then
         lprops(16,elem) = .false.
      else
         lprops(16,elem) = .false.
      end if
c
c                       geometric nonlinearity flag.
c                       the default value is on.
c
      geonl= elstor(7,elem)
      if( geonl .eq. 4HDEFA .or. geonl .eq. 4HTRUE ) then
         lprops(18,elem) = .true.
      else
         lprops(18,elem) = .false.
      end if
c
c                       interface-cohesive elements do not currently
c                       use a bbar type formulation (as in solids)
c
      lprops(19,elem) = .false.
c
c
c                       reference surface
c                       for geometrical update for the cohesive element
c
      surf = elstor(10, elem)
      if( surf.eq.4HDEFA .or. surf.eq.4HO222 ) then
         iprops(26,elem) = 2
      else if( surf.eq.4HO111 ) then
         iprops(26,elem) = 1
      else if( surf.eq.4HO333 ) then
         iprops(26,elem) = 3
      else
         iprops(26,elem) = 2
      end if
c
c                       See all comments in elprp12. They apply here.
c
      matnum          = elstor(2,elem)
      props(7,elem)   = matprp(32,matnum)
      props(8,elem)   = matprp(33,matnum)
      props(9,elem)   = matprp(34,matnum)
      props(13,elem)  = matprp(35,matnum)
      props(14,elem)  = matprp(36,matnum)
      props(15,elem)  = matprp(37,matnum)
      props(28,elem)  = matprp(38,matnum)
      props(20,elem)  = matprp(39,matnum)
      props(21,elem)  = matprp(40,matnum)
      props(22,elem)  = matprp(41,matnum)
c
c                 for symmetric problems ie, when reference surface is top
c                 or bottom, beta_coh (props(23)) should be zero for the
c                 "exponential" model. For the PPR option, it should
c                 be -2.0 or 2.0 (for only shear or only normal mode) but
c                 not zero.
c                 if user has input a non-consistent value, serve a warning
c                 and set a correct default value.
c
      props(23,elem)   = matprp(42,matnum)
      iprops(27,elem)  = matprp(44,matnum)
      lprops(10,elem)  = .true.
      exponential_type = iprops(27,elem) .eq. 4
      ppr_type         = iprops(27,elem) .eq. 6
      cavit_type       = iprops(27,elem) .eq. 7
c
      if( exponential_type ) then
        if ( surf.eq.4HO111 .and. matprp(42,matnum) .ne. zero ) then
       call errmsg2(33,elem,'top',matprp(42,matnum),dumd)
         props(23,elem) = zero
        end if
        if ( surf.eq.4HO333 .and. matprp(42,matnum) .ne. zero ) then
       call errmsg2(33,elem,'bottom',matprp(42,matnum),dumd)
         props(23,elem) = zero
        end if
      end if
c
      if( ppr_type ) then
        if ( surf.eq.4HO111 .and. matprp(42,matnum) .eq. zero ) then
       call errmsg2(78,elem,'top',matprp(42,matnum),dumd)
        end if
        if ( surf.eq.4HO333 .and. matprp(42,matnum) .eq. zero ) then
       call errmsg2(78,elem,'bottom',matprp(42,matnum),dumd)
        end if
      end if
c
c                  debug flag, allow cutback request from material,
c                  no segmental-type stress-strain curve.
c
      iprops(24,elem) = 0
      if ( lmtprp(13,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 1 )
      if ( lmtprp(22,matnum) ) iprops(24,elem) =
     &                         ior( iprops(24,elem), 2 )
      iprops(25,elem) = matprp(9,matnum)
c
c                      mark element killable if user has that
c                      property set for material.
c
      iprops(30,elem) = 0
      if ( lmtprp(23,matnum) ) iprops(30,elem) =
     &                         ior( iprops(30,elem), 2 )
c
c                      compression multiplier, element is not yet killed,
c                      fgm properties
c
      props(29,elem)  = matprp(54,matnum)
      iprops(32,elem) = 0
      props(33,elem)  = matprp(46,matnum)
      props(34,elem)  = matprp(47,matnum)
      props(35,elem)  = matprp(48,matnum)
      props(36,elem)  = matprp(49,matnum)
      props(37,elem)  = matprp(50,matnum)
      iprops(38,elem) = matnum
      props(39,elem)  = matprp(51,matnum)
      props(40,elem)  = matprp(52,matnum)
      props(41,elem)  = matprp(53,matnum)
c
c                     for PPR formulation, use spaces for FGM properties
c                     to store the PPR properties. Just keeps from
c                     making props table longer.
c
      if ( ppr_type ) then
          props(35,elem)  = matprp(92,matnum)
          props(36,elem)  = matprp(93,matnum)
          props(37,elem)  = matprp(94,matnum)
          props(39,elem)  = matprp(95,matnum)
          props(40,elem)  = matprp(96,matnum)
          props(41,elem)  = matprp(97,matnum)
      end if
c
      return
      end
c
c     ****************************************************************
c     *                                                              *
c     *                subroutine elem_set_segmental                 *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 02/28/00                   *
c     *                                                              *
c     *     this subroutine sets up the element definition when the  *
c     *     associated material has properties given by a set of     *
c     *     segmental stress strain curves                           *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine elem_set_segmental( elem, matnum )
      use global_data ! old common.main
      use main_data
      use segmental_curves
      implicit integer (a-z)
c
c                 local variables
c
      character :: dums
      real dumr
      double precision
     &    dumd
      logical local_debug
      data local_debug / .false. /

c
c
c                 if the material associated with the element has
c                 a set of stress-strain curves, store the set number
c                 for the curves and set the
c                 yield stress for the element material using the
c                 first curve in the set.
c
      if ( .not. lmtprp(24,matnum) )  return
c
      curve_set = matprp(45,matnum)
      if ( curve_set .le. 0 .or. curve_set .gt.
     &                  num_seg_curve_sets ) then
         call errmsg( 223, curve_set, dums, dumr, dumd )
         call die_gracefully
         stop
      end if
c
      iprops(21,elem)  = curve_set
      curve_one        = seg_curve_table(2,curve_set)
      props(23,elem)   = seg_curves(1,2,curve_one)
      if ( local_debug ) then
           write(*,*) '>> set up segmental curve set. element: ',elem
           write(*,*) '   > curve_set, curve_one, props(23,elem): ',
     &                  curve_set, curve_one, props(23,elem)
      end if
c
      return
      end

